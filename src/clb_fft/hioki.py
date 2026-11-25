"""Loading of HIOKI PQ ONE waveform exports.

PQ ONE splits a long waveform export into several files::

    MotorConVibracion1.csv
    MotorConVibracion1_001.csv
    MotorConVibracion1_002.csv
    ...

Every file starts with a header row (``EventNo,Date,Time,U1,I1,...``) and
the first file carries a UTF-8 byte-order mark. The waveform is stored as a
sequence of *event records* (``EventNo``) of a fixed number of samples; the
records are close to, but not guaranteed to be, contiguous in time (see
:mod:`clb_fft.quality`).
"""

from __future__ import annotations

import re
import tarfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

__all__ = ["Recording", "find_parts", "load_hioki_csv", "load_hioki_archive"]

_PART_RE = re.compile(r"^(?P<stem>.+)_(?P<idx>\d{3})$")
_META_COLUMNS = ("EventNo", "Date", "Time")


@dataclass
class Recording:
    """A multi-channel waveform recording.

    Attributes
    ----------
    signals:
        Mapping channel name (``"I1"``, ``"U1"``...) to a 1-D float array.
    fs:
        Sampling frequency in Hz.
    event_no:
        Event record number of every sample (``None`` if not available).
    start_time:
        Timestamp of the first sample (``None`` if not available).
    source:
        Files the recording was read from.
    """

    signals: dict[str, np.ndarray]
    fs: float
    event_no: np.ndarray | None = None
    start_time: pd.Timestamp | None = None
    source: list[str] = field(default_factory=list)

    def __getitem__(self, channel: str) -> np.ndarray:
        return self.signals[channel]

    @property
    def channels(self) -> list[str]:
        return list(self.signals)

    @property
    def n_samples(self) -> int:
        return len(next(iter(self.signals.values())))

    @property
    def duration(self) -> float:
        """Duration in seconds."""
        return self.n_samples / self.fs

    @property
    def time(self) -> np.ndarray:
        """Time axis in seconds, starting at zero."""
        return np.arange(self.n_samples) / self.fs

    def record_boundaries(self) -> np.ndarray:
        """Indices where a new event record starts (excluding index 0)."""
        if self.event_no is None:
            return np.array([], dtype=int)
        return np.flatnonzero(np.diff(self.event_no) != 0) + 1


def find_parts(first: str | Path) -> list[Path]:
    """Return ``first`` followed by its ``_001``, ``_002``... continuation files.

    ``first`` may also be any of the continuation files; the search always
    starts from the base file. Numbering may go beyond ``_009``.
    """
    first = Path(first)
    m = _PART_RE.match(first.stem)
    stem = m.group("stem") if m else first.stem
    base = first.with_name(stem + first.suffix)
    if not base.is_file():
        raise FileNotFoundError(base)
    parts = [base]
    i = 1
    while (p := first.with_name(f"{stem}_{i:03d}{first.suffix}")).is_file():
        parts.append(p)
        i += 1
    return parts


def _infer_fs(time_col: pd.Series) -> float:
    """Infer the sampling rate from the millisecond-resolution ``Time`` column.

    PQ ONE repeats the same timestamp for every sample inside one
    millisecond, so the run length of identical timestamps is ``fs / 1000``.
    """
    values = time_col.to_numpy()
    change = np.flatnonzero(values[1:] != values[:-1]) + 1
    if len(change) < 3:
        raise ValueError("Cannot infer sampling rate: pass fs explicitly.")
    runs = np.diff(change)
    return float(np.median(runs)) * 1000.0


def _frames_to_recording(
    frames: list[pd.DataFrame],
    sources: list[str],
    channels: list[str] | None,
    fs: float | None,
) -> Recording:
    df = pd.concat(frames, ignore_index=True)
    df.columns = [c.strip().lstrip("﻿") for c in df.columns]
    data_cols = [c for c in df.columns if c not in _META_COLUMNS]
    if channels is None:
        channels = data_cols
    missing = set(channels) - set(df.columns)
    if missing:
        raise KeyError(f"Channels not found: {sorted(missing)}; available: {data_cols}")

    if fs is None:
        if "Time" not in df.columns:
            raise ValueError("No 'Time' column: pass fs explicitly.")
        fs = _infer_fs(df["Time"])

    start = None
    if {"Date", "Time"} <= set(df.columns):
        try:
            start = pd.to_datetime(
                df["Date"].iloc[0] + " " + df["Time"].iloc[0], format="%d-%b-%y %H:%M:%S.%f"
            )
        except (ValueError, TypeError):
            start = None

    return Recording(
        signals={c: df[c].to_numpy(dtype=float) for c in channels},
        fs=fs,
        event_no=df["EventNo"].to_numpy() if "EventNo" in df.columns else None,
        start_time=start,
        source=sources,
    )


def load_hioki_csv(
    path: str | Path,
    channels: list[str] | None = None,
    fs: float | None = None,
) -> Recording:
    """Load a (possibly multi-part) PQ ONE CSV export.

    Parameters
    ----------
    path:
        Path to the first file of the export (or any of its parts).
    channels:
        Columns to load, e.g. ``["I1"]``. Defaults to all data columns.
    fs:
        Sampling frequency in Hz. Inferred from the ``Time`` column if omitted.
    """
    parts = find_parts(path)
    frames = [pd.read_csv(p, encoding="utf-8-sig") for p in parts]
    return _frames_to_recording(frames, [str(p) for p in parts], channels, fs)


def load_hioki_archive(
    path: str | Path,
    channels: list[str] | None = None,
    fs: float | None = None,
) -> Recording:
    """Load a PQ ONE export packed in a ``.tar``/``.tar.xz`` archive.

    All ``*.csv`` members are read in name order (``x.csv``, ``x_001.csv``...)
    without extracting them to disk.
    """
    frames, names = [], []
    with tarfile.open(path) as tar:
        members = sorted(
            (m for m in tar.getmembers() if m.isfile() and m.name.lower().endswith(".csv")),
            key=lambda m: m.name,
        )
        if not members:
            raise ValueError(f"No CSV files in {path}")
        for m in members:
            fh = tar.extractfile(m)
            assert fh is not None
            frames.append(pd.read_csv(fh, encoding="utf-8-sig"))
            names.append(f"{path}::{m.name}")
    return _frames_to_recording(frames, names, channels, fs)
