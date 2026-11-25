import pytest

from clb_fft.faults import (
    BearingGeometry,
    Motor,
    all_fault_frequencies,
    bearing_characteristic,
    broken_bar,
    eccentricity_mixed,
    rotor_slot_harmonics,
    slip_from_speed,
)


def test_slip_from_speed():
    # 4-pole 50 Hz machine: ns = 1500 rpm
    assert slip_from_speed(1470, 50, 2) == pytest.approx(0.02)


def test_motor_rotation_frequency():
    m = Motor.from_speed(50, 2, 1470)
    assert m.fr == pytest.approx(24.5)
    assert m.speed_rpm == pytest.approx(1470)


def test_broken_bar_uses_given_supply_frequency():
    f = sorted(x.frequency for x in broken_bar(60.0, 0.02, n=2))
    assert f == pytest.approx([55.2, 57.6, 62.4, 64.8])


def test_broken_bar_count():
    assert len(broken_bar(50, 0.01, n=5)) == 10


def test_mixed_eccentricity():
    m = Motor(50, 2, 0.02)
    f = sorted(x.frequency for x in eccentricity_mixed(m, n=1))
    assert f == pytest.approx([25.5, 74.5])


def test_rotor_slot_harmonic_principal():
    # R = 28, p = 2, s = 0.02 -> 28*0.98/2*50 = 686 Hz; +/- 50
    m = Motor(50, 2, 0.02, rotor_bars=28)
    f = sorted(x.frequency for x in rotor_slot_harmonics(m, nd=(0,), nw=(1,)))
    assert f == pytest.approx([636.0, 736.0])


def test_rotor_slot_requires_bar_count():
    assert rotor_slot_harmonics(Motor(50, 2, 0.02)) == []


def test_bearing_frequencies_6205():
    # SKF 6205: 9 balls, d = 7.94 mm, D = 39.04 mm. Known ratios to shaft speed:
    # BPFO 3.585, BPFI 5.415, BSF 2.357 (approx), FTF 0.398
    b = bearing_characteristic(BearingGeometry(9, 7.94, 39.04), fr=1.0)
    assert b["BPFO"] == pytest.approx(3.585, abs=0.01)
    assert b["BPFI"] == pytest.approx(5.415, abs=0.01)
    assert b["FTF"] == pytest.approx(0.398, abs=0.005)
    assert b["BSF"] == pytest.approx(2.357, abs=0.02)


def test_all_frequencies_sorted():
    m = Motor(50, 2, 0.02, rotor_bars=28, bearing=BearingGeometry(9, 7.94, 39.04))
    f = [x.frequency for x in all_fault_frequencies(m)]
    assert f == sorted(f)
    assert {x.family for x in all_fault_frequencies(m)} >= {"grid", "broken_bar", "bearing", "rotor_slot"}
