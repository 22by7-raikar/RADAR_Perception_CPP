#!/usr/bin/env python3
"""Render a bird's-eye view of tracking_demo output.

    ./build/tracking_demo --scenario intersection --out out_tracking
    python3 tools/plot_tracking.py out_tracking docs/images/tracking_intersection.png
"""

import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
CLASS_COLOR = {"vehicle": "#2a78d6", "pedestrian": "#eb6834", "unknown": "#8a8984"}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK_2, "axes.labelcolor": INK, "axes.titlecolor": INK,
    "xtick.color": INK_2, "ytick.color": INK_2, "grid.color": GRID,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 11, "axes.titlesize": 13, "axes.titleweight": "bold", "legend.frameon": False,
})


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "out_tracking")
    img = Path(sys.argv[2] if len(sys.argv) > 2 else "docs/images/tracking.png")
    title = sys.argv[3] if len(sys.argv) > 3 else "Multi-object tracking"
    img.parent.mkdir(parents=True, exist_ok=True)

    dets = read(out / "detections_world.csv")
    truth = read(out / "truth_objects.csv")
    tracks = read(out / "tracks.csv")
    ego = read(out / "ego.csv")

    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.scatter([float(d["x"]) for d in dets], [float(d["y"]) for d in dets], s=2,
               color=GRID, linewidths=0, zorder=1, label="Detections (all frames)")

    gt_paths = defaultdict(list)
    for g in truth:
        if g["class"] != "static":
            gt_paths[g["id"]].append((float(g["x"]), float(g["y"])))
    for i, path in enumerate(gt_paths.values()):
        xs, ys = zip(*path)
        ax.plot(xs, ys, color=INK_2, lw=6, alpha=0.15, solid_capstyle="round", zorder=2,
                label="Ground truth" if i == 0 else None)

    trk_paths = defaultdict(list)
    trk_class = {}
    for t in tracks:
        if t["confirmed"] == "1" and t["class"] != "static":
            trk_paths[t["id"]].append((float(t["x"]), float(t["y"])))
            trk_class[t["id"]] = t["class"]
    seen = set()
    for tid, path in trk_paths.items():
        if len(path) < 5:
            continue
        cls = trk_class[tid]
        xs, ys = zip(*path)
        ax.plot(xs, ys, color=CLASS_COLOR[cls], lw=1.8, zorder=3,
                label=f"Track: {cls}" if cls not in seen else None)
        ax.scatter([xs[-1]], [ys[-1]], s=22, color=CLASS_COLOR[cls], zorder=4)
        seen.add(cls)

    ex, ey = [float(e["x"]) for e in ego], [float(e["y"]) for e in ego]
    ax.plot(ex, ey, color=INK, lw=2, zorder=3)
    ax.scatter([ex[-1]], [ey[-1]], marker="^", s=110, color=INK, zorder=5, label="Sensor")

    ax.set_aspect("equal")
    ax.set(xlabel="x [m]", ylabel="y [m]", title=title)
    ax.grid(True, lw=0.5)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0))
    fig.savefig(img, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {img}")


if __name__ == "__main__":
    main()
