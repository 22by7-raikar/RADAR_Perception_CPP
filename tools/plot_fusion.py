#!/usr/bin/env python3
"""Camera view with camera detections, projected radar tracks and radar points.

    ./build/fusion_demo --out out_fusion
    python3 tools/plot_fusion.py out_fusion docs/images/fusion_camera_view.png
"""

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt

INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
CLASS_COLOR = {"vehicle": "#2a78d6", "pedestrian": "#eb6834", "unknown": "#8a8984", "static": "#8a8984"}
W, H = 1920, 1080

plt.rcParams.update({"font.size": 11, "axes.titlesize": 13, "axes.titleweight": "bold",
                     "figure.facecolor": SURFACE, "savefig.facecolor": SURFACE,
                     "axes.titlecolor": INK, "legend.frameon": False})


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "out_fusion")
    img = Path(sys.argv[2] if len(sys.argv) > 2 else "docs/images/fusion_camera_view.png")
    with open(out / "camera_frame.csv") as f:
        rows = list(csv.DictReader(f))

    fig, ax = plt.subplots(figsize=(13, 5.5))
    ax.add_patch(patches.Rectangle((0, 0), W, 520, color="#eef2f6", zorder=0))     # sky
    ax.add_patch(patches.Rectangle((0, 520), W, H - 520, color="#ecebe7", zorder=0))  # ground

    pts = [r for r in rows if r["kind"] == "point"]
    sc = ax.scatter([float(r["u0"]) for r in pts], [float(r["v0"]) for r in pts],
                    c=[float(r["label"]) for r in pts], cmap="Blues_r", vmin=0, vmax=160, s=14,
                    zorder=2, label="Radar detections")
    for r in rows:
        u0, v0, u1, v1 = (float(r[k]) for k in ("u0", "v0", "u1", "v1"))
        if r["kind"] == "camera":
            ax.add_patch(patches.Rectangle((u0, v0), u1 - u0, v1 - v0, fill=False, lw=1.6,
                                           ls=(0, (4, 3)), ec=INK_2, zorder=3))
        elif r["kind"] == "track" and (r["matched"] == "1" or r["class"] not in ("static", "unknown")):
            color = CLASS_COLOR.get(r["class"], INK_2)
            ax.add_patch(patches.Rectangle((u0, v0), u1 - u0, v1 - v0, fill=False, lw=2.2, ec=color,
                                           zorder=4))
            ax.text(u0, v0 - 8, f"{r['class']}  {r['label']}", color=INK, fontsize=9, zorder=5,
                    bbox=dict(boxstyle="round,pad=0.2", fc=SURFACE, ec=color, lw=1))

    handles = [patches.Patch(fill=False, ec=INK_2, ls="--", label="Camera detection"),
               patches.Patch(fill=False, ec=CLASS_COLOR["vehicle"], lw=2, label="Radar track: vehicle"),
               patches.Patch(fill=False, ec=CLASS_COLOR["pedestrian"], lw=2, label="Radar track: pedestrian")]
    ax.legend(handles=handles + [sc], loc="lower left", ncol=4, fontsize=10)
    ax.set_xlim(0, W)
    tracks = [r for r in rows if r["kind"] in ("camera", "track")]
    top = min(float(r["v0"]) for r in tracks) - 60
    bottom = max(float(r["v1"]) for r in tracks) + 80
    ax.set_ylim(min(H, bottom), max(0, top))
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.set_title("Camera view: radar tracks projected through the calibrated extrinsic")
    cb = fig.colorbar(sc, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("Radar range [m]")
    fig.savefig(img, dpi=140, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {img}")


if __name__ == "__main__":
    main()
