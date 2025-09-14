#!/usr/bin/env python3
"""
Batch runner for the VDC executable across multiple isovalues, storing all
outputs and logs in a fresh output directory.

Example:
  python tools/run_vdc_batch.py \
    --input data/sphere-32.nrrd \
    --isovalues -1.0 0.0 1.0 \
    --vdc-bin buildLocal/vdc \
    --format off \
    --output-dir outputs/sphere_batch

You can pass additional vdc flags by appending them after known args, e.g.:
  python tools/run_vdc_batch.py --input data/sphere-32.nrrd --isovalues 0 \
    -- -supersample 2 -ply
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import os
import shlex
import subprocess
import sys
from typing import List


def _iso_tag(val: float) -> str:
    """Return a filesystem-safe tag for an isovalue.

    Examples:
      -1.25 -> "m1p25"; 0.0 -> "0"; 2.5 -> "2p5"
    """
    s = f"{val:g}"
    return s.replace("-", "m").replace(".", "p")


def _ensure_abs(path: str) -> str:
    return os.path.abspath(os.path.expanduser(path))


def _default_outdir(input_file: str) -> str:
    """Default output dir: next to the input file.

    Example: for /path/to/data/vol.nrrd -> /path/to/data/vol_vdc_<timestamp>
    """
    ts = _dt.datetime.now().strftime("%Y%m%d-%H%M%S")
    base = os.path.splitext(os.path.basename(input_file))[0]
    parent = os.path.dirname(_ensure_abs(input_file))
    return os.path.join(parent, f"{base}_vdc_{ts}")


def run_once(
    vdc_bin: str,
    input_file: str,
    isovalue: float,
    out_dir: str,
    fmt: str,
    vdc_args: List[str],
) -> int:
    """Run vdc for a single isovalue, writing outputs/logs into out_dir.

    Returns the vdc process return code.
    """
    tag = _iso_tag(isovalue)
    input_base = os.path.splitext(os.path.basename(input_file))[0]
    output_basename = f"{input_base}_iso_{tag}"

    # Prefer running inside out_dir so any auxiliary outputs land there.
    cwd = out_dir

    fmt_flag = ["-off"] if fmt.lower() == "off" else ["-ply"]
    # vdc expects full filename when using -o; it will NOT append extension then.
    output_filename = output_basename + (".off" if fmt.lower() == "off" else ".ply")

    cmd = [
        vdc_bin,
        *fmt_flag,
        *vdc_args,
        "-o",
        output_filename,
        f"{isovalue:g}",
        input_file,
    ]

    # Log file per isovalue
    log_path = os.path.join(out_dir, f"log_iso_{tag}.txt")

    with open(log_path, "w", encoding="utf-8") as lf:
        lf.write("# Command\n")
        lf.write(" ".join(shlex.quote(c) for c in cmd) + "\n\n")
        lf.write("# Output\n")
        lf.flush()

        proc = subprocess.run(
            cmd,
            cwd=cwd,
            text=True,
            stdout=lf,
            stderr=subprocess.STDOUT,
        )

    return proc.returncode


def main(argv: List[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description=(
            "Run ./vdc across multiple isovalues and store outputs/logs "
            "in a new directory. Unknown flags are forwarded to vdc."
        )
    )
    p.add_argument(
        "--vdc-bin",
        default="buildLocal/vdc",
        help="Path to vdc executable (default: buildLocal/vdc)",
    )
    p.add_argument(
        "--input",
        required=True,
        help="Input NRRD/NHDR file path",
    )
    p.add_argument(
        "--isovalues",
        nargs="+",
        type=float,
        required=True,
        help="List of isovalues (e.g., --isovalues -1 0 1)",
    )
    p.add_argument(
        "--format",
        choices=["off", "ply"],
        default="off",
        help="Output mesh format (default: off)",
    )
    p.add_argument(
        "--output-dir",
        help="Directory to create for this run (default: outputs/<input>_<timestamp>)",
    )
    p.add_argument(
        "--allow-existing",
        action="store_true",
        help="Allow using an existing output directory (append new results)",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned commands without executing",
    )

    # Accept and forward any additional flags to vdc without validation.
    args, vdc_args = p.parse_known_args(argv)

    vdc_bin = _ensure_abs(args.vdc_bin)
    input_file = _ensure_abs(args.input)

    if not os.path.isfile(vdc_bin):
        print(f"Error: vdc binary not found at: {vdc_bin}", file=sys.stderr)
        return 2
    if not os.path.isfile(input_file):
        print(f"Error: input file not found: {input_file}", file=sys.stderr)
        return 2

    out_dir = args.output_dir or _default_outdir(input_file)
    out_dir = _ensure_abs(out_dir)

    if os.path.exists(out_dir) and not args.allow_existing:
        print(
            f"Error: output directory already exists: {out_dir}\n"
            f"Use --allow-existing to reuse it, or choose a new path.",
            file=sys.stderr,
        )
        return 2

    os.makedirs(out_dir, exist_ok=True)

    # Save run metadata
    meta = {
        "timestamp": _dt.datetime.now().isoformat(timespec="seconds"),
        "vdc_bin": vdc_bin,
        "input": input_file,
        "isovalues": args.isovalues,
        "format": args.format,
        "vdc_args": vdc_args,
        "cwd": out_dir,
    }
    with open(os.path.join(out_dir, "run_config.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)

    # Execute per isovalue
    failures = []
    planned_cmds = []
    for iso in args.isovalues:
        tag = _iso_tag(iso)
        input_base = os.path.splitext(os.path.basename(input_file))[0]
        output_basename = f"{input_base}_iso_{tag}"
        fmt_flag = ["-off"] if args.format == "off" else ["-ply"]
        output_filename = output_basename + (".off" if args.format == "off" else ".ply")
        cmd = [
            vdc_bin,
            *fmt_flag,
            *vdc_args,
            "-o",
            output_filename,
            f"{iso:g}",
            input_file,
        ]
        planned_cmds.append(" ".join(shlex.quote(c) for c in cmd))

    if args.dry_run:
        print("Planned commands (one per isovalue):")
        for c in planned_cmds:
            print("  " + c)
        print(f"Output directory: {out_dir}")
        return 0

    for iso in args.isovalues:
        rc = run_once(vdc_bin, input_file, iso, out_dir, args.format, vdc_args)
        if rc != 0:
            failures.append((iso, rc))

    # Write a brief summary
    summary_path = os.path.join(out_dir, "run_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as sf:
        sf.write("VDC batch run summary\n")
        sf.write(f"Output dir: {out_dir}\n")
        sf.write(f"Format: {args.format}\n")
        sf.write("Commands:\n")
        for c in planned_cmds:
            sf.write("  " + c + "\n")
        if failures:
            sf.write("\nFailures:\n")
            for iso, rc in failures:
                sf.write(f"  iso={iso:g} -> rc={rc}\n")

    if failures:
        print(
            f"Completed with {len(failures)} failures. See {summary_path}.",
            file=sys.stderr,
        )
        return 1

    print(f"All runs completed. Outputs in: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
