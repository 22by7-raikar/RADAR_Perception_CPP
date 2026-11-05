#!/usr/bin/env python3
"""Trajectory and drift figure from odometry_demo output.

    ./build/odometry_demo --out out_odometry
    python3 tools/plot_odometry.py out_odometry docs/images/odometry.png
"""

import csv
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
SERIES = {"imu": ("IMU only", "#8a8984"), "radar": ("Radar only", "#eb6834"),
          "fused": ("Radar + IMU (ESKF)", "#2a78d6")}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK_2, "axes.labelcolor": INK, "axes.titlecolor": INK,
    "xtick.color": INK_2, "ytick.color": INK_2, "grid.color": GRID,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 11, "axes.titlesize": 13, "axes.titleweight": "bold", "legend.frameon": False,
})


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "out_odometry")
    img = Path(sys.argv[2] if len(sys.argv) > 2 else "docs/images/odometry.png")
    img.parent.mkdir(parents=True, exist_ok=True)
    with open(out / "trajectory.csv") as f:
        rows = list(csv.DictReader(f))
    col = lambda k: [float(r[k]) for r in rows]
    t = col("time")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={"width_ratios": [1.15, 1]})
    ax1.plot(col("true_x"), col("true_y"), color=INK, lw=5, alpha=0.18, label="Ground truth",
             solid_capstyle="round")
    for key, (label, color) in SERIES.items():
        ax1.plot(col(key + "_x"), col(key + "_y"), color=color, lw=1.8, label=label)
        err = [math.hypot(float(r[key + "_x"]) - float(r["true_x"]),
                          float(r[key + "_y"]) - float(r["true_y"])) for r in rows]
        ax2.plot(t, err, color=color, lw=1.8, label=label)
        ax2.text(t[-1] + 0.4, err[-1], f"{err[-1]:.1f} m", va="center", color=INK, fontsize=10)
    ax1.set_aspect("equal")
    ax1.set(xlabel="x [m]", ylabel="y [m]", title="Ego trajectory, 30 s drive")
    ax1.grid(True, lw=0.5)
    ax1.legend(loc="upper left")
    ax2.set(xlabel="Time [s]", ylabel="Position error [m]", title="Drift vs. ground truth",
            xlim=(0, t[-1] * 1.12))
    ax2.grid(True, lw=0.5)
    fig.tight_layout()
    fig.savefig(img, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {img}")


if __name__ == "__main__":
    main()
