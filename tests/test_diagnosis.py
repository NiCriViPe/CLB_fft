import numpy as np
import pytest

from clb_fft.diagnosis import assess_broken_bars, detect_peaks, diagnose, label_peaks
from clb_fft.faults import Motor, broken_bar
from clb_fft.spectrum import averaged_amplitude_spectrum

FS = 5_000.0


def motor_current(f1=50.0, slip=0.03, sideband_db=-40.0, dur=40.0, noise=1e-4, seed=0):
    """Synthetic stator current with broken-bar sidebands at f1(1 +/- 2s)."""
    rng = np.random.default_rng(seed)
    t = np.arange(int(dur * FS)) / FS
    a = 10 ** (sideband_db / 20)
    x = np.sin(2 * np.pi * f1 * t)
    for f in (f1 * (1 - 2 * slip), f1 * (1 + 2 * slip)):
        x += a * np.sin(2 * np.pi * f * t + rng.uniform(0, 6))
    return x + noise * rng.standard_normal(len(t))


@pytest.mark.parametrize(
    "level, condition",
    [(-60, "excellent"), (-45, "moderate"), (-33, "two bars cracked/broken likely")],
)
def test_broken_bar_severity(level, condition):
    spec = averaged_amplitude_spectrum(motor_current(sideband_db=level), FS, resolution=0.05)
    r = assess_broken_bars(spec, 50.0, 0.03)
    assert r.delta_db == pytest.approx(-level, abs=1.0)
    assert r.condition == condition


def test_peaks_are_labelled():
    spec = averaged_amplitude_spectrum(motor_current(), FS, resolution=0.05)
    peaks = label_peaks(detect_peaks(spec, 40, 60), broken_bar(50, 0.03, 1), tol_hz=0.1)
    labels = {p.label for p in peaks if p.family == "broken_bar"}
    assert labels == {"f1(1-2s)", "f1(1+2s)"}


def test_each_peak_labelled_once_and_unknown_kept():
    spec = averaged_amplitude_spectrum(motor_current(), FS, resolution=0.05)
    table = diagnose(spec, Motor(50, 2, 0.03), fmin=1, fmax=200, prominence_db=10)
    assert table.frequency.is_unique
    assert (table.family == "grid").sum() == 1
    assert (table.family == "broken_bar").sum() == 2


def test_no_false_alarm_from_fundamental_skirt():
    # Healthy machine (no sidebands) with a supply frequency drifting 0.2 Hz
    t = np.arange(int(60 * FS)) / FS
    f_inst = 50 + 0.1 * np.sin(2 * np.pi * t / 60)
    x = np.sin(2 * np.pi * np.cumsum(f_inst) / FS)
    spec = averaged_amplitude_spectrum(x, FS, resolution=0.05)
    r = assess_broken_bars(spec, 50.0, 0.01)
    assert r.resolvable
    assert not r.lower[2] and not r.upper[2]
    assert r.condition.startswith("no sideband detected")


def test_unresolvable_sidebands_are_reported():
    spec = averaged_amplitude_spectrum(motor_current(slip=0.004), FS, resolution=0.2)
    r = assess_broken_bars(spec, 50.0, 0.004)
    assert not r.resolvable
    assert r.condition.startswith("not resolvable")
