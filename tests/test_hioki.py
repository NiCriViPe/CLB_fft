import tarfile

import numpy as np
import pytest

from clb_fft.hioki import find_parts, load_hioki_archive, load_hioki_csv
from clb_fft.quality import check_quality

FS = 20_000
REC = 400  # samples per event record (20 ms)


def write_export(folder, stem="Meas", n_parts=3, records_per_part=5, f1=50.0, gap_after=None):
    """Write a fake PQ ONE export (BOM, 1 ms timestamps, event records)."""
    t0 = 13 * 3600 + 44 * 60
    sample = 0
    paths = []
    for part in range(n_parts):
        name = f"{stem}.csv" if part == 0 else f"{stem}_{part:03d}.csv"
        lines = ["﻿EventNo,Date,Time,U1,I1"]
        for r in range(records_per_part):
            ev = part * records_per_part + r + 1
            if gap_after is not None and ev == gap_after + 1:
                sample += 37  # drop samples -> discontinuity
            for _ in range(REC):
                tsec = t0 + sample / FS
                ms = int(round((tsec % 1) * 1000 - 0.5)) % 1000
                hh, mm, ss = int(tsec // 3600), int(tsec % 3600 // 60), int(tsec % 60)
                i1 = 100 * np.sin(2 * np.pi * f1 * sample / FS)
                lines.append(f"{ev},17-Feb-23,{hh:02d}:{mm:02d}:{ss:02d}.{ms:03d},+0.05E+00,{i1:+.2f}E+00")
                sample += 1
        p = folder / name
        p.write_text("\n".join(lines) + "\n", encoding="utf-8")
        paths.append(p)
    return paths


def test_find_parts_any_entry_point(tmp_path):
    paths = write_export(tmp_path, n_parts=3)
    assert find_parts(paths[0]) == paths
    assert find_parts(paths[2]) == paths


def test_find_parts_beyond_nine(tmp_path):
    (tmp_path / "a.csv").write_text("x\n")
    for i in range(1, 13):
        (tmp_path / f"a_{i:03d}.csv").write_text("x\n")
    assert len(find_parts(tmp_path / "a.csv")) == 13


def test_load_infers_fs_and_strips_bom(tmp_path):
    paths = write_export(tmp_path)
    rec = load_hioki_csv(paths[0])
    assert rec.fs == FS
    assert rec.channels == ["U1", "I1"]
    assert rec.n_samples == 3 * 5 * REC
    assert len(rec.record_boundaries()) == 14
    assert rec.start_time is not None


def test_load_from_archive(tmp_path):
    paths = write_export(tmp_path)
    arc = tmp_path / "m.tar.xz"
    with tarfile.open(arc, "w:xz") as tar:
        for p in paths:
            tar.add(p, arcname=p.name)
    a = load_hioki_archive(arc, channels=["I1"])
    b = load_hioki_csv(paths[0], channels=["I1"])
    np.testing.assert_array_equal(a["I1"], b["I1"])


def test_unknown_channel(tmp_path):
    paths = write_export(tmp_path, n_parts=1)
    with pytest.raises(KeyError):
        load_hioki_csv(paths[0], channels=["I4"])


def test_quality_detects_gap(tmp_path):
    (tmp_path / "a").mkdir()
    (tmp_path / "b").mkdir()
    clean = check_quality(load_hioki_csv(write_export(tmp_path / "a")[0]))
    assert clean.discontinuous_joints == 0
    assert clean.f1_min == pytest.approx(50.0, abs=0.01)
    assert clean.f1_drift_hz < 0.01
    gappy = check_quality(load_hioki_csv(write_export(tmp_path / "b", gap_after=7)[0]))
    assert gappy.discontinuous_joints == 1
