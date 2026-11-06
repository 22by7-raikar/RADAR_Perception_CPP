#!/usr/bin/env python3
"""Figures for mimo_demo output.

    ./build/mimo_demo --out out_mimo
    python3 tools/plot_mimo.py out_mimo docs/images
"""

import csv
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, GREEN, GREY = "#2a78d6", "#eb6834", "#2f9e6e", "#8a8984"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK_2, "axes.labelcolor": INK, "axes.titlecolor": INK,
    "xtick.color": INK_2, "ytick.color": INK_2, "grid.color": GRID,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 11, "axes.titlesize": 13, "axes.titleweight": "bold", "legend.frameon": False,
})

TDM_SPAN = 23.2  # m/s, +/- Doppler span with 3 Tx


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def num(x):
    v = float(x)
    return v if math.isfinite(v) else None


def velocity_figure(out, img):
    rows = read(out / "velocity_sweep.csv")
    v = [float(r["velocity"]) for r in rows]
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.6))

    for ax in (a, b):
        ax.axvspan(-TDM_SPAN, TDM_SPAN, color=GRID, alpha=0.6, lw=0)
        ax.grid(True, axis="y")
        ax.set_xlabel("True range rate [m/s]")

    floor = 1e-3
    for key, color, label in (("az_err_none", GREY, "No motion compensation"),
                              ("az_err_comp", ORANGE, "Compensation, no unfolding"),
                              ("az_err_full", BLUE, "Compensation + unfolding")):
        ys = [max(floor, num(r[key])) if num(r[key]) is not None else None for r in rows]
        a.plot(v, ys, color=color, lw=2, label=label)
    a.set_yscale("log")
    a.set_ylabel("Azimuth error [deg]")
    a.set_title("Azimuth error, target at 25°")
    a.text(0, 1.6e-3, "TDM Doppler span", ha="center", color=INK_2, fontsize=9)
    a.legend(loc="lower right", fontsize=9)

    b.plot(v, v, color=GRID, lw=6, solid_capstyle="round", label="Truth")
    b.plot(v, [num(r["v_comp"]) for r in rows], color=ORANGE, lw=2, label="No unfolding (aliased)")
    b.plot(v, [num(r["v_full"]) for r in rows], color=BLUE, lw=2, ls="--",
           label="Unfolded")
    b.set_ylabel("Reported range rate [m/s]")
    b.set_title("Velocity beyond the TDM span")
    b.legend(loc="upper left", fontsize=9)

    fig.tight_layout()
    fig.savefig(img / "mimo_velocity.png", dpi=150)
    plt.close(fig)


def overhead_figure(out, img):
    rows = read(out / "overhead_detections.csv")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6), sharey=True)
    titles = {"3d": "3D radar (no elevation): height unknown",
              "4d": "4D radar (3 Tx × 4 Rx): height ± 2σ"}
    for ax, radar in zip(axes, ("3d", "4d")):
        ax.axhspan(0.25, 2.5, color=GRID, alpha=0.6, lw=0)
        ax.text(20, 1.35, "vehicle height band", color=INK_2, fontsize=9, va="center")
        pts = [r for r in rows if r["radar"] == radar and abs(float(r["y"])) < 1.8]
        for flag, color, label in (("1", ORANGE, "In-lane obstacle"),
                                   ("0", BLUE, "Pass under / over")):
            sel = [r for r in pts if r["obstacle"] == flag]
            xs = [float(r["x"]) for r in sel]
            if radar == "3d":
                ys = [0.7 for _ in sel]  # assumed at sensor height
                ax.scatter(xs, ys, s=18, color=color, label=label, zorder=3)
            else:
                ys = [float(r["height"]) for r in sel]
                err = [2 * float(r["height_std"]) for r in sel]
                ax.errorbar(xs, ys, yerr=err, fmt="o", ms=4, color=color, ecolor=color,
                            elinewidth=1, alpha=0.9, label=label, zorder=3)
        ax.set_title(titles[radar])
        ax.set_xlabel("Distance ahead [m]")
        ax.set_xlim(15, 115)
        ax.grid(True, axis="y")
        ax.legend(loc="upper right", fontsize=9)
    axes[0].set_ylabel("Height above road [m]")
    axes[0].set_ylim(-1.0, 7.5)
    # 3D: every return is placed at sensor height. 4D: at its measured height.
    labels = {
        axes[0]: (((63, 0.7), (48, 3.2)), ((99, 0.7), (86, 3.4)), ((36, 0.7), (22, 2.9))),
        axes[1]: (((66, 5.8), (40, 6.6)), ((100, 0.75), (83, 3.4)), ((40, 0.0), (22, -0.75))),
    }
    for ax, spots in labels.items():
        for text, (xy, xytext) in zip(("sign gantry", "stopped car", "metal plate"), spots):
            ax.annotate(text, xy, xytext=xytext, color=INK_2, fontsize=9,
                        arrowprops=dict(arrowstyle="-", color=GREY))
    fig.tight_layout()
    fig.savefig(img / "mimo_overhead.png", dpi=150)
    plt.close(fig)


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "out_mimo")
    img = Path(sys.argv[2] if len(sys.argv) > 2 else "docs/images")
    img.mkdir(parents=True, exist_ok=True)
    velocity_figure(out, img)
    overhead_figure(out, img)


if __name__ == "__main__":
    main()
