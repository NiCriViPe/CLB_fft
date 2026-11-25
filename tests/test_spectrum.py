import numpy as np
import pytest

from clb_fft.spectrum import (
    amplitude_spectrum,
    averaged_amplitude_spectrum,
    decimate,
    estimate_supply_frequency,
    psd,
    to_db,
)

FS = 20_000.0


def sine(f, a=1.0, dur=2.0, fs=FS, phase=0.3):
    t = np.arange(int(dur * fs)) / fs
    return a * np.sin(2 * np.pi * f * t + phase)


def test_frequency_axis_is_exact():
    s = amplitude_spectrum(np.zeros(1000), FS)
    assert s.resolution == pytest.approx(FS / 1000, rel=1e-12)
    assert s.f[-1] == pytest.approx(FS / 2)


@pytest.mark.parametrize("window", ["hann", "flattop", "boxcar"])
def test_amplitude_of_bin_centred_sine(window):
    s = amplitude_spectrum(sine(50.0, a=100.0), FS, window=window)
    f, a = s.value_at(50.0, 0.01)
    assert f == pytest.approx(50.0)
    assert a == pytest.approx(100.0, rel=1e-3)


def test_flattop_amplitude_off_bin():
    # Off-bin tone: flattop keeps amplitude error < 0.1 %
    s = amplitude_spectrum(sine(50.23, a=10.0), FS, window="flattop")
    assert s.value_at(50.23, 2.0)[1] == pytest.approx(10.0, rel=1e-3)


def test_amplitude_db_uses_20log10():
    assert to_db(np.array([1.0, 0.1]))[1] == pytest.approx(-20.0)
    assert to_db(np.array([1.0, 0.1]), power=True)[1] == pytest.approx(-10.0)


def test_window_reveals_sideband_that_rectangular_hides():
    # Off-bin fundamental + a sideband 60 dB below, 2 Hz away
    x = sine(50.27, 100.0, dur=4) + sine(48.0, 0.1, dur=4, phase=1.0)
    for window, visible in [("boxcar", False), ("hann", True)]:
        s = amplitude_spectrum(x, FS, window=window)
        _, a = s.value_at(48.0, 0.05)
        level = 20 * np.log10(a / 100.0)
        assert bool(abs(level - -60) < 3) is visible, window


def test_averaged_spectrum_amplitude_and_resolution():
    s = averaged_amplitude_spectrum(sine(50.0, 5.0, dur=20), FS, resolution=0.1)
    assert s.resolution == pytest.approx(0.1)
    assert s.value_at(50.0, 0.2)[1] == pytest.approx(5.0, rel=0.02)


def test_psd_integrates_to_power():
    x = sine(123.0, 2.0, dur=10)
    s = psd(x, FS, resolution=1.0)
    power = np.sum(s.values) * s.resolution
    assert power == pytest.approx(2.0**2 / 2, rel=0.01)


def test_decimate_suppresses_aliasing():
    # 990 Hz tone: with q=20 (new fs = 1 kHz) naive x[::q] aliases it to 10 Hz
    x = sine(990.0, 1.0, dur=4) + sine(50.0, 1.0, dur=4)
    naive = amplitude_spectrum(x[::20], FS / 20)
    good_x, fs2 = decimate(x, FS, 20)
    good = amplitude_spectrum(good_x, fs2)
    assert fs2 == 1000.0
    assert naive.value_at(10.0, 0.5)[1] > 0.5
    assert good.value_at(10.0, 0.5)[1] < 1e-3
    assert good.value_at(50.0, 0.5)[1] == pytest.approx(1.0, rel=0.01)


def test_decimate_rejects_large_prime():
    with pytest.raises(ValueError):
        decimate(np.zeros(1000), FS, 13)


@pytest.mark.parametrize("f1", [49.87, 50.0, 50.123, 59.9])
def test_estimate_supply_frequency(f1):
    nominal = 60.0 if f1 > 55 else 50.0
    assert estimate_supply_frequency(sine(f1, dur=10), FS, nominal) == pytest.approx(f1, abs=2e-3)
