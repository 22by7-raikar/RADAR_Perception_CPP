#!/usr/bin/env python3
"""Figures for impairments_demo output.

    ./build/impairments_demo --out out_impairments
    python3 tools/plot_impairments.py out_impairments docs/images
"""

import csv
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE, RED, GREY = "#2a78d6", "#eb6834", "#c23b3b", "#8a8984"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK_2, "axes.labelcolor": INK, "axes.titlecolor": INK,
    "xtick.color": INK_2, "ytick.color": INK_2, "grid.color": GRID,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 11, "axes.titlesize": 13, "axes.titleweight": "bold", "legend.frameon": False,
})

VELOCITY_RES = 1.0864  # m/s per Doppler bin, default radar
RANGE_RES = 0.9993     # m per range bin


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def interference(out, img):
    chirp = read(out / "interference_chirp.csv")
    prof = read(out / "interference_profile.csv")
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.4))

    env = np.array([float(r["envelope"]) for r in chirp])
    res = np.array([float(r["residual"]) for r in chirp])
    blanked = np.array([r["blanked"] == "1" for r in chirp])
    n = np.arange(len(env))
    for i in np.flatnonzero(blanked):
        a.axvspan(i - 0.5, i + 0.5, color=ORANGE, alpha=0.25, lw=0)
    a.plot(n, 10 * np.log10(np.maximum(env, 1e-9) / np.median(env)), color=GREY, lw=1,
           label="Raw samples (strong echo dominates)")
    a.plot(n, 10 * np.log10(np.maximum(res, 1e-9) / np.median(res)), color=INK, lw=1,
           label="Residual, echoes notched out")
    a.plot([], [], color=ORANGE, alpha=0.4, lw=8, label="Blanked")
    a.set_xlabel("ADC sample in chirp")
    a.set_ylabel("Power vs. median [dB]")
    a.set_title("One chirp: interference burst")
    a.set_ylim(-25, 45)
    a.legend(loc="upper left", fontsize=9)
    a.grid(True, axis="y")

    rng = np.array([float(r["range_bin"]) for r in prof]) * RANGE_RES
    for key, color, label, lw in (("interfered_db", RED, "Interfered", 1.5),
                                  ("mitigated_db", BLUE, "Blanked + IMAT", 2.0),
                                  ("clean_db", GREY, "Clean", 1.2)):
        a_ = np.array([float(r[key]) for r in prof])
        b.plot(rng, a_ - np.median([float(r["clean_db"]) for r in prof]), color=color, lw=lw,
               label=label, ls="--" if key == "clean_db" else "-")
    b.set_xlabel("Range [m]")
    b.set_ylabel("Mean power over Doppler [dB]")
    b.set_title("Noise floor after the range-Doppler FFT")
    b.legend(loc="upper right", fontsize=9)
    b.grid(True, axis="y")
    fig.tight_layout()
    fig.savefig(img / "impairments_interference.png", dpi=150)
    plt.close(fig)


def clutter(out, img):
    with open(out / "clutter_power.csv") as f:
        header = f.readline().strip().split(",")
        rows, cols = int(header[1]), int(header[3])
        power = np.loadtxt(f, delimiter=",").reshape(rows, cols)
    dets = read(out / "clutter_detections.csv")
    truth = {t["name"]: (float(t["velocity"]), float(t["range"])) for t in read(out / "clutter_truth.csv")}
    r0, r1 = 55, 115
    vel = (np.arange(cols) - cols // 2) * VELOCITY_RES
    vsel = (vel > -14) & (vel < 11)
    crop = power[r0:r1, vsel]
    crop = crop - np.median(power)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, name in zip(axes, ("CA", "OS")):
        ax.imshow(crop, origin="lower", aspect="auto", cmap="Greys", vmin=0, vmax=45,
                  extent=[vel[vsel][0] - VELOCITY_RES / 2, vel[vsel][-1] + VELOCITY_RES / 2,
                          r0 * RANGE_RES, r1 * RANGE_RES])
        sel = [d for d in dets if d["cfar"] == name]
        ax.scatter([float(d["velocity"]) for d in sel], [float(d["range"]) for d in sel], s=28,
                   facecolors="none", edgecolors=ORANGE, linewidths=1.5, label="Detections")
        ax.set_title(f"{name}-CFAR")
        ax.set_xlabel("Range rate [m/s]")
        labels = {"vegetation": ((0.0, 106), (3.0, 110)),
                  "pedestrian": (truth["pedestrian"], (-6.5, 60)),
                  "car": (truth["car"], (-13.5, 74)),
                  "motorbike (weak)": (truth["motorbike"], (-13.5, 95)),
                  "same motorbike, alone": (truth["reference"], (2.0, 64))}
        for text, (xy, xytext) in labels.items():
            ax.annotate(text, xy, xytext=xytext, color=BLUE, fontsize=9,
                        arrowprops=dict(arrowstyle="-", color=BLUE))
    axes[0].set_ylabel("Range [m]")
    axes[1].legend(loc="lower left", fontsize=9)
    fig.tight_layout()
    fig.savefig(img / "impairments_cfar.png", dpi=150)
    plt.close(fig)


def multipath(out, img):
    dets = read(out / "multipath_detections.csv")
    rail = read(out / "multipath_rail.csv")[0]
    fig, ax = plt.subplots(figsize=(12, 4.8))
    styles = {
        "rail": dict(color=GREY, s=10, label="Guardrail (static)"),
        "real": dict(color=BLUE, s=60, label="Real object"),
    }
    for src, st in styles.items():
        sel = [d for d in dets if d["source"] == src]
        ax.scatter([float(d["x"]) for d in sel], [float(d["y"]) for d in sel], zorder=3, **st)
    ghosts = [d for d in dets if d["source"].startswith("ghost")]
    for flagged, color, label in (("1", ORANGE, "Ghost, flagged"), ("0", RED, "Ghost, missed")):
        sel = [d for d in ghosts if d["flagged"] == flagged]
        ax.scatter([float(d["x"]) for d in sel], [float(d["y"]) for d in sel], marker="x", s=70,
                   color=color, linewidths=2, label=label, zorder=4)
    if rail["valid"] == "1":
        xs = np.array([0, 150])
        ax.plot(xs, float(rail["slope"]) * xs + float(rail["offset"]), color=INK_2, lw=1.5,
                ls="--", label="Fitted guardrail")
    ax.axhspan(-3.5, 5, color=GRID, alpha=0.4, lw=0)
    ax.scatter([0], [0], marker=">", s=120, color=INK, zorder=5, label="Radar")
    ax.set_xlim(-3, 150)
    ax.set_ylim(-22, 8)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("Multipath off a guardrail (one frame)")
    ax.legend(loc="lower left", fontsize=9, ncol=3)
    fig.tight_layout()
    fig.savefig(img / "impairments_multipath.png", dpi=150)
    plt.close(fig)


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "out_impairments")
    img = Path(sys.argv[2] if len(sys.argv) > 2 else "docs/images")
    img.mkdir(parents=True, exist_ok=True)
    interference(out, img)
    clutter(out, img)
    multipath(out, img)


if __name__ == "__main__":
    main()
