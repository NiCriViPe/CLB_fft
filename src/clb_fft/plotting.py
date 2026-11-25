"""Matplotlib figures. Functions return the Figure and never call ``show()``,
so they work in scripts, notebooks, tests (Agg backend) and GUIs alike."""

from __future__ import annotations

from collections.abc import Iterable

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.figure import Figure

from .faults import FaultFrequency
from .spectrum import Spectrum

__all__ = ["plot_spectrum", "plot_frequency_track", "FAMILY_COLORS"]

FAMILY_COLORS = {
    "grid": "#7f7f7f",
    "broken_bar": "#d62728",
    "eccentricity": "#1f77b4",
    "eccentricity_dynamic": "#17becf",
    "rotor_slot": "#9467bd",
    "stator_turns": "#ff7f0e",
    "bearing": "#2ca02c",
}


def plot_spectrum(
    spectra: Spectrum | Iterable[tuple[str, Spectrum]],
    fmin: float = 0.0,
    fmax: float | None = None,
    markers: Iterable[FaultFrequency] = (),
    title: str = "",
    ylim: tuple[float, float] = (-100, 5),
) -> Figure:
    """Plot one or several spectra in dB (each relative to its own maximum)."""
    spectra = [("", spectra)] if isinstance(spectra, Spectrum) else list(spectra)
    fig, ax = plt.subplots(figsize=(10, 4.5), layout="constrained")
    for name, s in spectra:
        b = s.band(fmin, fmax if fmax is not None else s.f[-1])
        ref = float(np.max(s.values))
        ax.plot(b.f, b.db(ref=ref), lw=0.8, label=name or None)
    seen = set()
    for m in markers:
        if fmin <= m.frequency <= (fmax or np.inf):
            c = FAMILY_COLORS.get(m.family, "k")
            ax.axvline(
                m.frequency,
                color=c,
                ls="--",
                lw=0.8,
                alpha=0.7,
                label=m.family if m.family not in seen else None,
            )
            seen.add(m.family)
    unit = "PSD [dB rel. max]" if spectra[0][1].kind == "psd" else "Amplitude [dB rel. max]"
    ax.set(xlabel="Frequency [Hz]", ylabel=unit, title=title, ylim=ylim)
    ax.grid(True, which="major", alpha=0.6)
    ax.grid(True, which="minor", ls=":", alpha=0.3)
    ax.minorticks_on()
    if any(name for name, _ in spectra) or seen:
        ax.legend(loc="upper right", fontsize=8)
    return fig


def plot_frequency_track(track: np.ndarray, window_s: float, title: str = "") -> Figure:
    fig, ax = plt.subplots(figsize=(10, 3), layout="constrained")
    ax.plot(np.arange(len(track)) * window_s, track, lw=0.8)
    ax.set(xlabel="Time [s]", ylabel="Supply frequency [Hz]", title=title)
    ax.grid(True, alpha=0.5)
    return fig
