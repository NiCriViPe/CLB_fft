"""Reproduce the README figures from the sample measurement.

    python examples/method_comparison.py [archive] [output_dir]

1. before_after_f1.png - the legacy spectrum (one 280 s rectangular FFT,
   amplitude in 10*log10) vs. the current one (Welch/Hann, 20*log10) around f1.
2. supply_frequency.png - supply frequency per event record, showing the
   drift that limits useful resolution.
"""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from clb_fft import averaged_amplitude_spectrum, decimate, load_hioki_archive
from clb_fft.quality import frequency_track

BLUE, ORANGE, INK, MUTED = "#257de0", "#f47230", "#1f1f1f", "#6b6b6b"


def style(ax):
    ax.grid(True, color="#e3e3e3", lw=0.6)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#bdbdbd")
    ax.tick_params(colors=MUTED, labelsize=9)


def main(archive="Measurements/NH-I/MotorConVibracion1.tar.xz", out="docs/img"):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rec = load_hioki_archive(archive, channels=["I1"])
    x, fs = rec["I1"], rec.fs

    # Legacy: rectangular window over the full record, 10*log10 of an amplitude
    n = len(x)
    legacy = 2.0 / n * np.abs(np.fft.rfft(x))
    f_legacy = np.fft.rfftfreq(n, 1 / fs)
    legacy_db = 10 * np.log10(legacy / legacy.max())

    xd, fsd = decimate(x, fs, 20)
    new = averaged_amplitude_spectrum(xd, fsd, resolution=0.2)
    new_db = new.db()

    fig, axes = plt.subplots(2, 1, figsize=(9, 5.6), sharex=True, sharey=True, layout="constrained")
    band = (40, 60)
    m = (f_legacy >= band[0]) & (f_legacy <= band[1])
    axes[0].plot(f_legacy[m], legacy_db[m], color=ORANGE, lw=0.7)
    axes[0].set_title(
        "Before: one 280 s FFT, rectangular window, 10·log10 of amplitude",
        loc="left",
        fontsize=10,
        color=INK,
    )
    b = new.band(*band)
    axes[1].plot(b.f, new_db[(new.f >= band[0]) & (new.f <= band[1])], color=BLUE, lw=1.2)
    axes[1].set_title(
        "After: Welch average (0.2 Hz, Hann), anti-aliased decimation, 20·log10",
        loc="left",
        fontsize=10,
        color=INK,
    )
    for ax in axes:
        style(ax)
        ax.set_ylabel("dB rel. f1", color=MUTED, fontsize=9)
        ax.set_ylim(-110, 5)
    axes[1].set_xlabel("Frequency [Hz]", color=MUTED, fontsize=9)
    axes[1].set_xlim(*band)
    fig.suptitle(
        "Stator current I1 around the 50 Hz fundamental (MotorConVibracion1)",
        x=0.01,
        ha="left",
        fontsize=11,
        color=INK,
    )
    fig.savefig(out / "before_after_f1.png", dpi=150)

    bnd = rec.record_boundaries()
    track = frequency_track(x, fs, window_s=0.2, boundaries=bnd)
    t = np.r_[0, bnd][: len(track)] / fs
    fig, ax = plt.subplots(figsize=(9, 2.8), layout="constrained")
    ax.plot(t, track, color=BLUE, lw=0.8)
    style(ax)
    ax.set_xlabel("Time [s]", color=MUTED, fontsize=9)
    ax.set_ylabel("f1 [Hz]", color=MUTED, fontsize=9)
    ax.set_title(
        f"Supply frequency per 0.28 s event record: drift {np.ptp(track):.2f} Hz over {rec.duration:.0f} s",
        loc="left",
        fontsize=10,
        color=INK,
    )
    fig.savefig(out / "supply_frequency.png", dpi=150)
    print(f"figures written to {out}/")


if __name__ == "__main__":
    main(*sys.argv[1:])
