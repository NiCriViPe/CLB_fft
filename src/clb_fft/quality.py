"""Signal-quality checks to run *before* trusting a spectrum.

Two issues seen in real PQ ONE exports:

1. The waveform is a chain of short event records (e.g. 1000 x 0.28 s).
   Some joints are discontinuous; each discontinuity is a broadband
   impulse in the spectrum and, being roughly periodic, can create spurious
   lines spaced by ``1/record_length`` Hz.
2. The supply frequency drifts during long captures. A spectrum whose
   resolution is finer than that drift is meaningless around f1.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .hioki import Recording

__all__ = ["QualityReport", "check_quality", "boundary_errors", "frequency_track"]


@dataclass
class QualityReport:
    n_records: int
    record_length_s: float | None
    discontinuous_joints: int
    f1_min: float
    f1_max: float
    f1_drift_hz: float

    @property
    def recommended_resolution(self) -> float:
        """Finest resolution (Hz) worth asking for, given the supply drift."""
        return max(self.f1_drift_hz, 0.01)

    def summary(self) -> str:
        lines = [
            f"event records          : {self.n_records}"
            + (f" x {self.record_length_s:.3f} s" if self.record_length_s else ""),
            f"discontinuous joints   : {self.discontinuous_joints}",
            f"supply frequency       : {self.f1_min:.3f} .. {self.f1_max:.3f} Hz "
            f"(drift {self.f1_drift_hz:.3f} Hz)",
            f"suggested resolution   : >= {self.recommended_resolution:.3f} Hz",
        ]
        return "\n".join(lines)


def boundary_errors(x: np.ndarray, boundaries: np.ndarray) -> tuple[np.ndarray, float]:
    """Linear-prediction error at record joints vs. the 99.9th percentile inside.

    Returns ``(errors_at_joints, interior_threshold)``.
    """
    x = np.asarray(x, dtype=float)
    pred_err = np.abs(x[2:] - (2 * x[1:-1] - x[:-2]))
    threshold = float(np.percentile(pred_err, 99.9))
    b = boundaries[(boundaries >= 2)]
    joint_err = np.abs(x[b] - (2 * x[b - 1] - x[b - 2]))
    return joint_err, threshold


def _zc_frequency(seg: np.ndarray, fs: float) -> float | None:
    seg = seg - np.mean(seg)
    s = np.signbit(seg)
    zc = np.flatnonzero(s[:-1] & ~s[1:])  # rising crossings
    if len(zc) < 3:
        return None
    t = zc + seg[zc] / (seg[zc] - seg[zc + 1])  # sub-sample crossing instants
    return fs * (len(t) - 1) / (t[-1] - t[0])


def frequency_track(
    x: np.ndarray, fs: float, window_s: float = 1.0, boundaries: np.ndarray | None = None
) -> np.ndarray:
    """Fundamental frequency per window, from interpolated zero crossings.

    If ``boundaries`` (record starts) are given and records hold at least
    ``window_s`` seconds, the estimate is made per record so that joints
    between records never fall inside an analysis window.
    """
    if boundaries is not None and len(boundaries):
        segs = np.split(np.asarray(x, dtype=float), boundaries)
        if np.median([len(g) for g in segs]) >= window_s * fs:
            return np.array([f for g in segs if (f := _zc_frequency(g, fs)) is not None])
    n = int(window_s * fs)
    out = [_zc_frequency(x[i : i + n], fs) for i in range(0, len(x) - n + 1, n)]
    return np.array([f for f in out if f is not None])


def check_quality(rec: Recording, channel: str = "I1") -> QualityReport:
    """Run all checks on one channel of a recording."""
    x = rec[channel]
    bnd = rec.record_boundaries()
    n_records = len(bnd) + 1
    rec_len = None
    if len(bnd):
        rec_len = float(np.median(np.diff(np.r_[0, bnd, rec.n_samples]))) / rec.fs
    joint_err, thr = boundary_errors(x, bnd)
    # >= 0.2 s (~10 cycles) per window for a stable zero-crossing estimate
    track = frequency_track(x, rec.fs, window_s=0.2, boundaries=bnd)
    if track.size == 0:
        raise ValueError("Signal too short or without zero crossings to track frequency.")
    return QualityReport(
        n_records=n_records,
        record_length_s=rec_len,
        discontinuous_joints=int(np.sum(joint_err > 2 * thr)),
        f1_min=float(track.min()),
        f1_max=float(track.max()),
        f1_drift_hz=float(track.max() - track.min()),
    )
