#!/usr/bin/env python3
"""Render figures from radar_demo output.

    ./build/radar_demo --out out
    python3 tools/plot_results.py out docs/images
"""

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Radar parameters of the default RadarConfig (keep in sync with config.hpp).
C = 299_792_458.0
FC, B, T_CHIRP, T_IDLE = 77e9, 150e6, 12e-6, 2e-6
NR, ND = 256, 128
R_RES = C / (2 * B)
V_RES = (C / FC) / (2 * ND * (T_CHIRP + T_IDLE))

INK, INK_2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
BLUE, ORANGE = "#2a78d6", "#eb6834"

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK_2, "axes.labelcolor": INK, "axes.titlecolor": INK,
    "xtick.color": INK_2, "ytick.color": INK_2, "grid.color": GRID,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 11, "axes.titlesize": 13, "axes.titleweight": "bold",
    "legend.frameon": False,
})


def read_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


def column(rows, key, cast=float):
    return np.array([cast(r[key]) for r in rows])


def plot_rdm(out, img, dets):
    rdm = np.loadtxt(out / "power_map_db.csv", delimiter=",")
    rdm -= np.median(rdm)  # dB relative to the noise floor
    vel = (np.arange(ND) - ND // 2) * V_RES
    rng = np.arange(NR) * R_RES

    fig, ax = plt.subplots(figsize=(8, 6))
    im = ax.pcolormesh(vel, rng, rdm, cmap="Blues", shading="nearest", vmin=-5, vmax=45)
    f0 = [d for d in dets if int(d["frame"]) == 0]
    ax.scatter([float(d["velocity_mps"]) for d in f0], [float(d["range_m"]) for d in f0],
               s=160, facecolors="none", edgecolors=ORANGE, linewidths=2, label="Detection")
    ax.set(xlabel="Radial velocity [m/s]  (negative = approaching)", ylabel="Range [m]",
           title="Range-Doppler map with detections (frame 0)")
    ax.legend(loc="upper right")
    fig.colorbar(im, ax=ax, label="Power above noise floor [dB]")
    fig.savefig(img / "range_doppler_map.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_point_cloud(img, dets, truth):
    fig, ax = plt.subplots(figsize=(7, 7))
    tr_frames = column(truth, "frame", int)
    tr_r, tr_az = column(truth, "range_m"), np.radians(column(truth, "azimuth_deg"))
    tx, ty = tr_r * np.cos(tr_az), tr_r * np.sin(tr_az)
    # Plot y (left) on the horizontal axis flipped, so the view matches a top-down
    # picture with the sensor at the bottom looking up.
    ax.scatter(-ty, tx, s=90, color=INK_2, alpha=0.18, linewidths=0, label="Ground truth")
    ax.scatter(-column(dets, "y_m"), column(dets, "x_m"), s=12, color=BLUE, linewidths=0,
               label="Detections")
    last = tr_frames.max()
    ax.scatter([0], [0], marker="^", s=120, color=INK, label="Sensor")
    for r in (50, 100, 150, 200, 250):
        a = np.linspace(-np.pi / 2, np.pi / 2, 200)
        ax.plot(r * np.sin(a), r * np.cos(a), color=GRID, lw=0.8, zorder=0)
        ax.text(r * np.sin(np.radians(48)), r * np.cos(np.radians(48)), f"{r} m", color=INK_2,
                fontsize=9, ha="left", va="bottom")
    ax.set_aspect("equal")
    ax.set(xlim=(-210, 210), ylim=(-5, 260), xlabel="Lateral [m]  (right +)",
           ylabel="Forward [m]", title=f"Point cloud over {last + 1} frames")
    ax.legend(loc="upper left")
    fig.savefig(img / "point_cloud.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def plot_latency(out, img):
    rows = read_csv(out / "timings.csv")
    stages = [("range_doppler_ms", "Range-Doppler FFT"), ("integration_ms", "Integration"),
              ("cfar_ms", "CA-CFAR"), ("angle_ms", "Angle FFT"), ("peak_grouping_ms", "Peak grouping")]
    means = [column(rows, k).mean() for k, _ in stages]
    total = column(rows, "total_ms")

    fig, ax = plt.subplots(figsize=(8, 3.2))
    y = np.arange(len(stages))[::-1]
    ax.barh(y, means, color=BLUE, height=0.6)
    for yi, m in zip(y, means):
        ax.text(m + max(means) * 0.01, yi, f"{m:.3f} ms", va="center", color=INK, fontsize=10)
    ax.set_yticks(y, [s for _, s in stages])
    ax.set(xlabel="Mean time per frame [ms]", xlim=(0, max(means) * 1.25),
           title=f"Per-stage latency (total mean {total.mean():.2f} ms, "
                 f"p99 {np.percentile(total, 99):.2f} ms)")
    ax.grid(True, axis="x", lw=0.6)
    fig.savefig(img / "latency.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "out")
    img = Path(sys.argv[2] if len(sys.argv) > 2 else "docs/images")
    img.mkdir(parents=True, exist_ok=True)
    dets = read_csv(out / "detections.csv")
    truth = read_csv(out / "truth.csv")
    plot_rdm(out, img, dets)
    plot_point_cloud(img, dets, truth)
    plot_latency(out, img)
    print(f"Figures written to {img}/")


if __name__ == "__main__":
    main()
