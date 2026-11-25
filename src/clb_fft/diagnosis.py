"""Peak detection, labelling against fault families, and severity indicators."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd
from scipy.signal import find_peaks

from .faults import FaultFrequency, Motor, all_fault_frequencies, broken_bar
from .spectrum import Spectrum, to_db

__all__ = [
    "Peak",
    "detect_peaks",
    "label_peaks",
    "BrokenBarAssessment",
    "assess_broken_bars",
    "diagnose",
    "peaks_table",
]


@dataclass
class Peak:
    frequency: float
    db: float  # relative to the fundamental
    family: str = "unclassified"
    label: str = ""
    expected: float | None = None


def detect_peaks(
    spec: Spectrum,
    fmin: float = 1.0,
    fmax: float | None = None,
    prominence_db: float = 6.0,
    floor_db: float = -90.0,
    min_distance_hz: float | None = None,
) -> list[Peak]:
    """Find spectral peaks, in dB relative to the largest component.

    ``prominence_db`` is measured against the local noise floor, which is
    more robust than an absolute height threshold.
    """
    s = spec.band(fmin, fmax if fmax is not None else spec.f[-1])
    db = to_db(s.values, ref=float(np.max(spec.values)), power=spec.kind == "psd")
    distance = None
    if min_distance_hz:
        distance = max(1, int(min_distance_hz / spec.resolution))
    idx, _ = find_peaks(db, prominence=prominence_db, height=floor_db, distance=distance)
    return [Peak(float(s.f[i]), float(db[i])) for i in idx]


def label_peaks(peaks: list[Peak], expected: list[FaultFrequency], tol_hz: float) -> list[Peak]:
    """Assign each peak to the nearest expected frequency within ``tol_hz``.

    A peak is labelled at most once (with its closest match); peaks with no
    match stay ``"unclassified"``.
    """
    if not expected:
        return peaks
    ef = np.array([e.frequency for e in expected])
    for pk in peaks:
        d = np.abs(ef - pk.frequency)
        j = int(np.argmin(d))
        if d[j] <= tol_hz:
            pk.family, pk.label, pk.expected = expected[j].family, expected[j].label, float(ef[j])
    return peaks


def peaks_table(peaks: list[Peak]) -> pd.DataFrame:
    return pd.DataFrame([asdict(p) for p in peaks])


# Indicative thresholds for the dB difference between f1 and the larger
# (1-2s)f1 / (1+2s)f1 sideband, after W. T. Thomson & M. Fenger (2001).
# Valid only for a motor at (near) full load; trending against a baseline
# of the *same* machine is always more reliable than absolute thresholds.
_BRB_TABLE = [
    (54.0, "excellent"),
    (48.0, "good"),
    (42.0, "moderate"),
    (36.0, "rotor bar crack or high-resistance joint likely"),
    (30.0, "two bars cracked/broken likely"),
    (0.0, "multiple broken bars / end-ring fault likely"),
]


@dataclass
class BrokenBarAssessment:
    f1: float
    slip: float
    lower: tuple[float, float, bool]  # (Hz, dB rel f1, detected as a real peak)
    upper: tuple[float, float, bool]
    delta_db: float  # dB between f1 and the larger *detected* sideband (or noise bound)
    condition: str
    resolvable: bool

    @property
    def graded(self) -> bool:
        """True if a sideband was detected and graded against the severity table."""
        return self.resolvable and (self.lower[2] or self.upper[2])

    def summary(self) -> str:
        def fmt(name, sb):
            what = "peak" if sb[2] else "no peak, level <="
            return f"{name} = {sb[0]:.3f} Hz: {what} {sb[1]:.1f} dB"

        return "\n".join(
            [
                f"f1 = {self.f1:.3f} Hz, s = {self.slip:.4f}",
                fmt("(1-2s)f1", self.lower),
                fmt("(1+2s)f1", self.upper),
                f"difference = {self.delta_db:.1f} dB -> {self.condition}"
                + (" (indicative, full-load thresholds)" if self.graded else ""),
            ]
        )


def _sideband(spec: Spectrum, f0: float, tol: float, ref: float, prominence_db: float):
    """Largest *local peak* within f0 +/- tol, else the max level there (noise bound)."""
    m = np.flatnonzero(np.abs(spec.f - f0) <= tol)
    lo, hi = max(m[0] - 5, 0), min(m[-1] + 6, len(spec.f))
    db = to_db(spec.values[lo:hi], ref=ref)
    idx, _ = find_peaks(db, prominence=prominence_db)
    idx = [i for i in idx if lo + i in set(m)]
    if idx:
        i = max(idx, key=lambda i: db[i])
        return float(spec.f[lo + i]), float(db[i]), True
    j = m[np.argmax(spec.values[m])]
    return float(spec.f[j]), float(to_db(spec.values[j : j + 1], ref=ref)[0]), False


def assess_broken_bars(
    spec: Spectrum,
    f1: float,
    slip: float,
    tol_hz: float | None = None,
    prominence_db: float = 6.0,
) -> BrokenBarAssessment:
    """Broken-rotor-bar indicator from the first pair of sidebands f1(1 +/- 2s).

    Guards against the classic false alarm of reading the fundamental's own
    leakage skirt as a sideband:

    * the sidebands must be at least ~4 bins away from f1 (``resolvable``);
    * a sideband only counts if it is a local peak with ``prominence_db``;
      otherwise its level is reported as an upper bound (noise floor).
    """
    if spec.kind != "amplitude":
        raise ValueError("Use an amplitude spectrum")
    sep = 2 * slip * f1
    tol = tol_hz if tol_hz is not None else max(1.5 * spec.resolution, 0.02)
    tol = min(tol, 0.3 * sep)
    f_fund, a_fund = spec.value_at(f1, max(2 * spec.resolution, 0.05))
    resolvable = sep >= 4 * spec.resolution and tol >= spec.resolution / 2
    lo_f, hi_f = (b.frequency for b in broken_bar(f1, slip, n=1))
    lower = _sideband(spec, lo_f, max(tol, spec.resolution / 2), a_fund, prominence_db)
    upper = _sideband(spec, hi_f, max(tol, spec.resolution / 2), a_fund, prominence_db)
    detected = [sb[1] for sb in (lower, upper) if sb[2]]
    if not resolvable:
        delta = float("nan")
        cond = (
            f"not resolvable: sidebands {sep:.3f} Hz from f1 need resolution "
            f"<= {sep / 4:.3f} Hz (have {spec.resolution:.3f} Hz)"
        )
    elif detected:
        delta = -max(detected)
        cond = next(c for thr, c in _BRB_TABLE if delta >= thr)
    else:
        delta = -max(lower[1], upper[1])
        cond = f"no sideband detected (> {delta:.0f} dB below f1)"
    return BrokenBarAssessment(f_fund, slip, lower, upper, delta, cond, resolvable)


def diagnose(spec: Spectrum, motor: Motor, tol_hz: float | None = None, **peak_kw) -> pd.DataFrame:
    """Detect peaks and label them with every applicable fault family."""
    tol = tol_hz if tol_hz is not None else max(3 * spec.resolution, 0.1)
    peaks = detect_peaks(spec, **peak_kw)
    return peaks_table(label_peaks(peaks, all_fault_frequencies(motor), tol))
