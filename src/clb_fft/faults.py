"""Characteristic fault frequencies of induction machines in the stator current.

All functions return a list of :class:`FaultFrequency`, each tagged with the
fault family and harmonic indices, so detected peaks can be labelled.

Notation
--------
f1 : supply frequency [Hz]
s  : slip [p.u.]
p  : pole **pairs**
fr : mechanical rotation frequency = f1 (1 - s) / p  [Hz]

References
----------
* W. T. Thomson, M. Fenger, "Current signature analysis to detect induction
  motor faults", IEEE Industry Applications Magazine, 7(4), 2001.
* M. E. H. Benbouzid, "A review of induction motors signature analysis as a
  medium for faults detection", IEEE Trans. Ind. Electron., 47(5), 2000.
* R. R. Schoen et al., "Motor bearing damage detection using stator current
  monitoring", IEEE Trans. Ind. Appl., 31(6), 1995.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, radians

__all__ = [
    "FaultFrequency",
    "Motor",
    "BearingGeometry",
    "slip_from_speed",
    "grid_harmonics",
    "broken_bar",
    "eccentricity_mixed",
    "rotor_slot_harmonics",
    "stator_shorted_turns",
    "bearing_characteristic",
    "bearing_current",
    "all_fault_frequencies",
]


@dataclass(frozen=True)
class FaultFrequency:
    family: str
    frequency: float
    label: str


@dataclass(frozen=True)
class BearingGeometry:
    """Rolling element bearing geometry (from the bearing catalogue)."""

    n_balls: int
    ball_diameter: float
    pitch_diameter: float
    contact_angle_deg: float = 0.0


@dataclass(frozen=True)
class Motor:
    """Operating point and construction data of the machine under test."""

    f1: float
    pole_pairs: int
    slip: float
    rotor_bars: int | None = None
    bearing: BearingGeometry | None = None

    @property
    def fr(self) -> float:
        """Mechanical rotation frequency [Hz]."""
        return self.f1 * (1.0 - self.slip) / self.pole_pairs

    @property
    def speed_rpm(self) -> float:
        return 60.0 * self.fr

    @classmethod
    def from_speed(cls, f1: float, pole_pairs: int, speed_rpm: float, **kw) -> Motor:
        return cls(f1=f1, pole_pairs=pole_pairs, slip=slip_from_speed(speed_rpm, f1, pole_pairs), **kw)


def slip_from_speed(speed_rpm: float, f1: float, pole_pairs: int) -> float:
    """Slip from measured shaft speed: s = (ns - n) / ns."""
    ns = 60.0 * f1 / pole_pairs
    return (ns - speed_rpm) / ns


def grid_harmonics(f1: float, n: int = 10) -> list[FaultFrequency]:
    """Supply harmonics k*f1, k = 1..n."""
    return [FaultFrequency("grid", k * f1, f"{k}f1") for k in range(1, n + 1)]


def broken_bar(f1: float, slip: float, n: int = 3) -> list[FaultFrequency]:
    """Broken rotor bar sidebands f1 (1 +/- 2 k s), k = 1..n."""
    out = []
    for k in range(1, n + 1):
        out.append(FaultFrequency("broken_bar", f1 * (1 - 2 * k * slip), f"f1(1-{2 * k}s)"))
        out.append(FaultFrequency("broken_bar", f1 * (1 + 2 * k * slip), f"f1(1+{2 * k}s)"))
    return out


def eccentricity_mixed(motor: Motor, n: int = 3) -> list[FaultFrequency]:
    """Mixed eccentricity / misalignment / load oscillation: f1 +/- k fr."""
    out = []
    for k in range(1, n + 1):
        out.append(FaultFrequency("eccentricity", abs(motor.f1 - k * motor.fr), f"f1-{k}fr"))
        out.append(FaultFrequency("eccentricity", motor.f1 + k * motor.fr, f"f1+{k}fr"))
    return out


def rotor_slot_harmonics(
    motor: Motor, nd: tuple[int, ...] = (0, 1), nw: tuple[int, ...] = (1, 3), k: int = 1
) -> list[FaultFrequency]:
    """Rotor slot / eccentricity harmonics (Thomson & Fenger 2001).

    f = f1 [ (k R +/- nd)(1 - s)/p +/- nw ]

    ``nd = 0`` static eccentricity (principal slot harmonics),
    ``nd = 1, 2...`` dynamic eccentricity; ``nw`` = stator MMF harmonic order.
    Needs the number of rotor bars. Also useful to *estimate the slip* when
    no tachometer is available.
    """
    if motor.rotor_bars is None:
        return []
    R, s, p, f1 = motor.rotor_bars, motor.slip, motor.pole_pairs, motor.f1
    out = []
    for d in nd:
        for sign_d in (+1, -1) if d else (+1,):
            for w in nw:
                for sign_w in (+1, -1):
                    f = f1 * ((k * R + sign_d * d) * (1 - s) / p + sign_w * w)
                    if f > 0:
                        lbl = f"RSH nd={sign_d * d:+d} nw={sign_w * w:+d}"
                        fam = "rotor_slot" if d == 0 else "eccentricity_dynamic"
                        out.append(FaultFrequency(fam, f, lbl))
    return out


def stator_shorted_turns(
    motor: Motor, n: tuple[int, ...] = (1, 2, 3), k: tuple[int, ...] = (1, 3)
) -> list[FaultFrequency]:
    """Shorted turns (Thomson): f = f1 [ n (1 - s)/p +/- k ]."""
    out = []
    for ni in n:
        for ki in k:
            for sign in (+1, -1):
                f = motor.f1 * (ni * (1 - motor.slip) / motor.pole_pairs + sign * ki)
                if f > 0:
                    out.append(FaultFrequency("stator_turns", f, f"n={ni} k={sign * ki:+d}"))
    return out


def bearing_characteristic(geom: BearingGeometry, fr: float) -> dict[str, float]:
    """Mechanical bearing defect frequencies (BPFO, BPFI, BSF, FTF) in Hz."""
    ratio = geom.ball_diameter / geom.pitch_diameter * cos(radians(geom.contact_angle_deg))
    nb = geom.n_balls
    return {
        "BPFO": nb / 2 * fr * (1 - ratio),
        "BPFI": nb / 2 * fr * (1 + ratio),
        "BSF": geom.pitch_diameter / (2 * geom.ball_diameter) * fr * (1 - ratio**2),
        "FTF": fr / 2 * (1 - ratio),
    }


def bearing_current(motor: Motor, m: int = 2) -> list[FaultFrequency]:
    """Bearing defects seen in the current: |f1 +/- m f_v| (Schoen 1995).

    Needs ``motor.bearing`` (returns ``[]`` otherwise). Without the full
    geometry, ``BearingGeometry`` can approximate BPFO ~ 0.4 Nb fr and
    BPFI ~ 0.6 Nb fr by using ``ball_diameter/pitch_diameter = 0.2``.
    Note that ``fr`` is the shaft speed, **not** the supply frequency.
    """
    if motor.bearing is None:
        return []
    out = []
    for name, fv in bearing_characteristic(motor.bearing, motor.fr).items():
        for mi in range(1, m + 1):
            out.append(FaultFrequency("bearing", abs(motor.f1 - mi * fv), f"f1-{mi}{name}"))
            out.append(FaultFrequency("bearing", motor.f1 + mi * fv, f"f1+{mi}{name}"))
    return out


def all_fault_frequencies(motor: Motor, n: int = 3) -> list[FaultFrequency]:
    """Every family that the motor data allows, sorted by frequency."""
    out = (
        grid_harmonics(motor.f1, 7)
        + broken_bar(motor.f1, motor.slip, n)
        + eccentricity_mixed(motor, n)
        + rotor_slot_harmonics(motor)
        + stator_shorted_turns(motor)
        + bearing_current(motor)
    )
    return sorted(out, key=lambda x: x.frequency)
