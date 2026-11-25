"""Spectral estimation for motor current signature analysis (MCSA).

Conventions
-----------
* Amplitude spectra are returned in **peak units** (a sine of amplitude
  ``A`` shows a peak of height ``A``), corrected for the window coherent
  gain.
* Decibels of an *amplitude* use ``20*log10``; decibels of a *power*
  (PSD) use ``10*log10``.
* Frequency axes come from :func:`numpy.fft.rfftfreq` (bin spacing exactly
  ``fs/N``).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import signal

__all__ = [
    "Spectrum",
    "amplitude_spectrum",
    "averaged_amplitude_spectrum",
    "psd",
    "decimate",
    "to_db",
    "estimate_supply_frequency",
]


@dataclass
class Spectrum:
    """A one-sided spectrum ``f`` (Hz) vs ``values``."""

    f: np.ndarray
    values: np.ndarray
    kind: str  # "amplitude" (peak units) or "psd" (units^2/Hz)

    @property
    def resolution(self) -> float:
        return float(self.f[1] - self.f[0])

    def db(self, ref: float | None = None) -> np.ndarray:
        """Spectrum in dB relative to ``ref`` (default: its maximum)."""
        power = self.kind == "psd"
        return to_db(self.values, ref=ref, power=power)

    def band(self, fmin: float, fmax: float) -> Spectrum:
        m = (self.f >= fmin) & (self.f <= fmax)
        return Spectrum(self.f[m], self.values[m], self.kind)

    def value_at(self, f0: float, tol: float) -> tuple[float, float]:
        """Largest value within ``f0 +/- tol``; returns ``(frequency, value)``."""
        m = np.abs(self.f - f0) <= tol
        if not m.any():
            i = int(np.argmin(np.abs(self.f - f0)))
            return float(self.f[i]), float(self.values[i])
        idx = np.flatnonzero(m)
        i = idx[np.argmax(self.values[idx])]
        return float(self.f[i]), float(self.values[i])


def to_db(x: np.ndarray, ref: float | None = None, power: bool = False) -> np.ndarray:
    """Convert to decibels.

    ``power=False`` (amplitudes) uses ``20*log10``; ``power=True`` uses
    ``10*log10``. ``ref`` defaults to ``max(x)`` so the peak is 0 dB.
    """
    x = np.asarray(x, dtype=float)
    ref = float(np.max(x)) if ref is None else float(ref)
    k = 10.0 if power else 20.0
    with np.errstate(divide="ignore"):
        return k * np.log10(np.maximum(x, np.finfo(float).tiny) / ref)


def amplitude_spectrum(
    x: np.ndarray, fs: float, window: str | tuple = "hann", detrend: bool = True
) -> Spectrum:
    """Single FFT of the whole signal, scaled to peak amplitude.

    Resolution is ``fs/len(x)``. Use a window: with the rectangular window
    the leakage skirt of the supply fundamental can easily hide fault
    sidebands 40-60 dB below it.
    """
    x = np.asarray(x, dtype=float)
    if detrend:
        x = x - x.mean()
    n = len(x)
    w = signal.get_window(window, n, fftbins=True)
    spec = np.abs(np.fft.rfft(x * w)) * 2.0 / w.sum()
    spec[0] /= 2.0
    if n % 2 == 0:
        spec[-1] /= 2.0
    return Spectrum(np.fft.rfftfreq(n, 1.0 / fs), spec, "amplitude")


def averaged_amplitude_spectrum(
    x: np.ndarray,
    fs: float,
    resolution: float = 0.05,
    window: str | tuple = "hann",
    overlap: float = 0.5,
) -> Spectrum:
    """Welch-averaged amplitude spectrum with a chosen resolution (Hz).

    Averaging over segments of ``1/resolution`` seconds trades resolution
    for robustness: in field recordings the supply frequency drifts
    (typically +/-0.1 Hz over minutes), so one huge FFT does not give the
    resolution it promises and smears the fundamental.
    """
    nperseg = int(round(fs / resolution))
    nperseg = min(nperseg, len(x))
    f, p = signal.welch(
        x,
        fs=fs,
        window=window,
        nperseg=nperseg,
        noverlap=int(nperseg * overlap),
        scaling="spectrum",
        detrend="constant",
    )
    # scaling="spectrum" gives squared RMS per bin -> peak amplitude
    return Spectrum(f, np.sqrt(2.0 * p), "amplitude")


def psd(
    x: np.ndarray,
    fs: float,
    resolution: float = 0.5,
    window: str | tuple = "hann",
    average: str = "mean",
) -> Spectrum:
    """Welch power spectral density (units^2/Hz)."""
    nperseg = min(int(round(fs / resolution)), len(x))
    f, p = signal.welch(x, fs=fs, window=window, nperseg=nperseg, detrend="constant", average=average)
    return Spectrum(f, p, "psd")


def decimate(x: np.ndarray, fs: float, q: int) -> tuple[np.ndarray, float]:
    """Downsample by ``q`` **with** an anti-aliasing filter.

    Taking every q-th sample (``x[::q]``) folds everything above the new
    Nyquist frequency (inverter switching, slot harmonics...) back into the
    low-frequency band, creating false peaks.
    """
    if q < 1:
        raise ValueError("q must be >= 1")
    # Decimate in stages of at most 10 (keeps each anti-alias filter short).
    stages, r = [], q
    for s in range(10, 1, -1):
        while r % s == 0:
            stages.append(s)
            r //= s
    if r != 1:
        raise ValueError(f"q={q} has a prime factor > 10; choose another factor")
    y = np.asarray(x, dtype=float)
    for s in stages:
        y = signal.decimate(y, s, ftype="fir", zero_phase=True)
    return y, fs / q


def estimate_supply_frequency(x: np.ndarray, fs: float, nominal: float = 50.0, search: float = 5.0) -> float:
    """Estimate the supply (fundamental) frequency with sub-bin accuracy.

    Uses a Hann-windowed FFT and quadratic interpolation of the
    log-magnitude around the largest peak in ``nominal +/- search`` Hz.
    """
    s = amplitude_spectrum(x, fs)
    m = np.flatnonzero(np.abs(s.f - nominal) <= search)
    k = m[np.argmax(s.values[m])]
    a, b, c = np.log(s.values[k - 1 : k + 2] + 1e-300)
    denom = a - 2 * b + c
    delta = 0.5 * (a - c) / denom if denom != 0 else 0.0
    return float(s.f[k] + delta * s.resolution)
