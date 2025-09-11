#!/usr/bin/env python3

import json
import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection

def load_case(path):
    with open(path, 'r') as f:
        return json.load(f)

def facet_normal(verts):
    # verts: Nx3 array
    if len(verts) < 3:
        return np.array([0.0,0.0,1.0])
    v0 = verts[0]
    v1 = verts[1]
    v2 = verts[2]
    n = np.cross(v1 - v0, v2 - v0)
    n_norm = np.linalg.norm(n)
    return n / n_norm if n_norm > 1e-12 else np.array([0.0,0.0,1.0])

def main():
    if len(sys.argv) < 2:
        print("Usage: tools/plot_modcyc.py modcyc_case.json")
        sys.exit(1)

    data = load_case(sys.argv[1])

    # Build vertices array
    vid_to_pt = {}
    for v in data["vertices"]:
        vid_to_pt[v["id"]] = np.array([v["x"], v["y"], v["z"]], dtype=float)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot facets and segments
    for f in data["facets"]:
        verts = [vid_to_pt[i] for i in f["verts"]]
        V = np.array(verts)
        # facet polygon
        poly = Poly3DCollection([V], alpha=0.2)
        poly.set_facecolor(np.random.rand(3,))
        ax.add_collection3d(poly)
        # boundary edges
        edges = []
        n = len(V)
        for i in range(n):
            edges.append([V[i], V[(i+1)%n]])
        ax.add_collection3d(Line3DCollection(edges, colors='k', linewidths=1.0, linestyles=':'))

        # before iso-segments
        sb = [np.array(seg, dtype=float) for seg in f.get("segments_before", [])]
        if sb:
            ax.add_collection3d(Line3DCollection(sb, colors='r', linewidths=2.0, label='before'))
        # after iso-segments
        sa = [np.array(seg, dtype=float) for seg in f.get("segments_after", [])]
        if sa:
            ax.add_collection3d(Line3DCollection(sa, colors='g', linewidths=2.0, label='after'))

        # show two incident cells as offset facet copies
        nrm = facet_normal(V)
        off = 0.05
        Vp = V + off*nrm
        Vm = V - off*nrm
        ax.add_collection3d(Line3DCollection([[Vp[i], Vp[(i+1)%n]] for i in range(n)], colors='b', linewidths=1.5))
        ax.add_collection3d(Line3DCollection([[Vm[i], Vm[(i+1)%n]] for i in range(n)], colors='b', linewidths=1.5))

    # Plot cycles per cell (optional), draw as polylines
    def draw_cycles(key, color, linestyle):
        cycles = data.get(key, [])
        for entry in cycles:
            for cyc in entry.get("cycles", []):
                P = np.array(cyc, dtype=float)
                segs = [[P[i], P[i+1]] for i in range(len(P)-1)]
                ax.add_collection3d(Line3DCollection(segs, colors=color, linewidths=2.0, linestyles=linestyle))

    draw_cycles("cell_cycles_before", '#cc0000', ':')
    draw_cycles("cell_cycles_after",  '#00aa00', ':')

    # Draw 3D prisms for the two incident cells per facet
    prisms = data.get("facet_prisms", [])
    for pr in prisms:
        top = np.array(pr.get("top", []), dtype=float)
        bot = np.array(pr.get("bottom", []), dtype=float)
        if len(top) == 0 or len(bot) == 0 or len(top) != len(bot):
            continue
        n = len(top)
        faces = []
        # top and bottom
        faces.append(top)
        faces.append(bot[::-1])
        # sides
        for i in range(n):
            a = top[i]
            b = top[(i+1) % n]
            c = bot[(i+1) % n]
            d = bot[i]
            faces.append(np.array([a, b, c, d]))
        poly = Poly3DCollection(faces, alpha=0.2)
        poly.set_facecolor(np.random.rand(3,))
        ax.add_collection3d(poly)

    # Overall aesthetics
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('Modify Cycles Test: iso-segments before (red) vs after (green), cells as prisms')
    # autoscale
    all_pts = np.array([p for p in vid_to_pt.values()])
    mins = all_pts.min(axis=0)
    maxs = all_pts.max(axis=0)
    center = 0.5*(mins+maxs)
    span = (maxs - mins).max()
    lims = np.array([center - 0.6*span, center + 0.6*span])
    ax.set_xlim(lims[:,0])
    ax.set_ylim(lims[:,1])
    ax.set_zlim(lims[:,2])

    plt.show()

if __name__ == '__main__':
    main()
