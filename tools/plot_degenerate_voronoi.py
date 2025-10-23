#!/usr/bin/env python3
"""
Utility to visualize degenerate Voronoi facet incidents detected in log_fuel.txt.

The script parses the diagnostic blocks emitted by vdc when Voronoi facets collapse,
extracts the Delaunay edge endpoints, incident tetrahedra, and circumcenters, and
produces 3D plots that highlight the geometric configuration responsible for
each degeneracy.
"""

import argparse
import itertools
import math
import re
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  Imported for side-effects


Point3 = Tuple[float, float, float]


def parse_log(path: Path) -> List[dict]:
    """Parse the log and return a list of degenerate facet cases."""
    cases: List[dict] = []
    current = None
    cell = None

    edge_pattern = re.compile(
        r"v([12]): ([\d\.\-]+) ([\d\.\-]+) ([\d\.\-]+) \(index: (\d+)\)"
    )
    cell_header_pattern = re.compile(r"Delaunay Cell (\d+)")
    vertex_pattern = re.compile(
        r"\[(\d)\]: \(([\d\.\-]+), ([\d\.\-]+), ([\d\.\-]+)\) index: (\d+)"
    )
    circumcenter_pattern = re.compile(
        r"Circumcenter: \(([\d\.\-]+), ([\d\.\-]+), ([\d\.\-]+)\)"
    )
    voronoi_center_pattern = re.compile(
        r"Delaunay vertex \(Voronoi cell center\): ([\d\.\-]+) ([\d\.\-]+) ([\d\.\-]+) \(index: (\d+)\)"
    )
    voronoi_cell_index_pattern = re.compile(r"- Voronoi cell index: (\d+)")

    with path.open("r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")

            if "[ERROR] Degenerate Voronoi facet detected:" in line:
                if current:
                    cases.append(current)
                current = {
                    "edge": {},
                    "cells": [],
                }
                cell = None
                continue

            if current is None:
                continue

            if "Delaunay edge endpoints" in line:
                continue

            match = edge_pattern.search(line)
            if match:
                key = f"v{match.group(1)}"
                current["edge"][key] = {
                    "point": tuple(map(float, match.group(2, 3, 4))),  # type: ignore[arg-type]
                    "index": int(match.group(5)),
                }
                continue

            match = voronoi_center_pattern.search(line)
            if match:
                current["voronoi_center"] = {
                    "point": tuple(map(float, match.group(1, 2, 3))),  # type: ignore[arg-type]
                    "index": int(match.group(4)),
                }
                continue

            match = cell_header_pattern.search(line)
            if match:
                cell = {"id": int(match.group(1)), "vertices": [], "circumcenter": None}
                current["cells"].append(cell)
                continue

            match = vertex_pattern.search(line)
            if match and cell is not None:
                cell["vertices"].append(
                    {
                        "slot": int(match.group(1)),
                        "point": tuple(map(float, match.group(2, 3, 4))),  # type: ignore[arg-type]
                        "index": int(match.group(5)),
                    }
                )
                continue

            match = circumcenter_pattern.search(line)
            if match and cell is not None:
                cell["circumcenter"] = tuple(
                    map(float, match.group(1, 2, 3))
                )  # type: ignore[arg-type]
                continue

            match = voronoi_cell_index_pattern.search(line)
            if match:
                current["voronoi_cell_index"] = int(match.group(1))

        if current:
            cases.append(current)

    return cases


def set_equal_aspect(ax: Axes3D, points: List[Point3]) -> None:
    """Set equal aspect ratio for 3D axes."""
    xs, ys, zs = zip(*points)
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)
    max_range = max(max_x - min_x, max_y - min_y, max_z - min_z)
    if max_range == 0:
        max_range = 1.0
    mid_x = (max_x + min_x) / 2.0
    mid_y = (max_y + min_y) / 2.0
    mid_z = (max_z + min_z) / 2.0
    half = max_range / 2.0
    ax.set_xlim(mid_x - half, mid_x + half)
    ax.set_ylim(mid_y - half, mid_y + half)
    ax.set_zlim(mid_z - half, mid_z + half)


def plot_case(case: dict, case_id: int, out_dir: Path) -> Path:
    """Create a 3D visualization for a single degenerate case."""
    fig = plt.figure(figsize=(8, 6))
    ax: Axes3D = fig.add_subplot(111, projection="3d")

    colors = ["tab:blue", "tab:orange", "tab:green", "tab:red", "tab:purple"]

    # Collect unique vertices over all cells
    unique_vertices: Dict[Point3, dict] = {}
    for cell in case["cells"]:
        for v in cell["vertices"]:
            unique_vertices.setdefault(v["point"], v)

    # Plot tetrahedra edges
    for cell_idx, cell in enumerate(case["cells"]):
        verts = [v["point"] for v in cell["vertices"]]
        color = colors[cell_idx % len(colors)]
        for p, q in itertools.combinations(verts, 2):
            ax.plot(
                [p[0], q[0]],
                [p[1], q[1]],
                [p[2], q[2]],
                color=color,
                linewidth=1.5,
                alpha=0.8,
            )
        # Circumcenter
        if cell.get("circumcenter"):
            c = cell["circumcenter"]
            ax.scatter(
                [c[0]], [c[1]], [c[2]], color=color, marker="X", s=80, label=f"Cell {cell_idx} circumcenter"
            )

    # Plot edge endpoints prominently
    for key, style in zip(["v1", "v2"], ["o", "^"]):
        vinfo = case["edge"].get(key)
        if vinfo:
            p = vinfo["point"]
            ax.scatter([p[0]], [p[1]], [p[2]], color="red", marker=style, s=80)
            ax.text(p[0], p[1], p[2], f"{key} (#{vinfo['index']})", color="red")

    # Draw the problematic Delaunay edge
    if "v1" in case["edge"] and "v2" in case["edge"]:
        p = case["edge"]["v1"]["point"]
        q = case["edge"]["v2"]["point"]
        ax.plot(
            [p[0], q[0]],
            [p[1], q[1]],
            [p[2], q[2]],
            color="red",
            linewidth=3,
            label="Problematic edge",
        )

    # Plot Voronoi center vertex
    voronoi_center = case.get("voronoi_center")
    if voronoi_center:
        p = voronoi_center["point"]
        ax.scatter([p[0]], [p[1]], [p[2]], color="black", marker="*", s=120, label="Voronoi cell center")
        ax.text(p[0], p[1], p[2], f"Voronoi center #{voronoi_center['index']}", color="black")

    # Annotate remaining vertices if not already
    for point, vinfo in unique_vertices.items():
        if point in [case["edge"].get("v1", {}).get("point"), case["edge"].get("v2", {}).get("point")]:
            continue
        ax.scatter([point[0]], [point[1]], [point[2]], color="gray", marker="s", s=40)
        ax.text(point[0], point[1], point[2], f"#{vinfo['index']}", color="gray")

    title = f"Degenerate Voronoi Facet Case {case_id}"
    if "voronoi_cell_index" in case:
        title += f" (Voronoi cell #{case['voronoi_cell_index']})"
    ax.set_title(title)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    set_equal_aspect(ax, list(unique_vertices.keys()) + [case["edge"]["v1"]["point"], case["edge"]["v2"]["point"]])
    ax.legend(loc="best")
    ax.view_init(elev=18, azim=35)
    fig.tight_layout()

    out_dir.mkdir(parents=True, exist_ok=True)
    output_path = out_dir / f"degenerate_case_{case_id:02d}.png"
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description="Visualize degenerate Voronoi facet cases.")
    parser.add_argument("--log", type=Path, default=Path("buildLocal/log_fuel.txt"), help="Path to the log file.")
    parser.add_argument(
        "--cases",
        type=int,
        nargs="+",
        default=[0, 1, 2],
        help="Indices of cases to plot (0-based).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path("buildLocal/degenerate_plots"),
        help="Where to store the generated images.",
    )
    args = parser.parse_args()

    cases = parse_log(args.log)
    if not cases:
        raise SystemExit(f"No degenerate facet cases found in {args.log}")

    outputs = []
    for case_id in args.cases:
        if case_id < 0 or case_id >= len(cases):
            print(f"[WARN] Case index {case_id} is out of range (0..{len(cases)-1}); skipping.")
            continue
        output_path = plot_case(cases[case_id], case_id, args.out_dir)
        outputs.append(output_path)
        print(f"[INFO] Wrote plot for case {case_id} to {output_path}")

    if not outputs:
        raise SystemExit("No plots were generated. Please supply valid case indices.")


if __name__ == "__main__":
    main()
