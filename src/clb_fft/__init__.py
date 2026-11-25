"""clb_fft - motor current signature analysis (MCSA) of HIOKI PQ ONE exports."""

__version__ = "0.2.0"

from .faults import BearingGeometry, Motor, all_fault_frequencies, slip_from_speed  # noqa: E402
from .hioki import Recording, load_hioki_archive, load_hioki_csv  # noqa: E402
from .spectrum import (  # noqa: E402
    amplitude_spectrum,
    averaged_amplitude_spectrum,
    decimate,
    estimate_supply_frequency,
    psd,
)

__all__ = [
    "BearingGeometry",
    "Motor",
    "Recording",
    "all_fault_frequencies",
    "amplitude_spectrum",
    "averaged_amplitude_spectrum",
    "decimate",
    "estimate_supply_frequency",
    "load_hioki_archive",
    "load_hioki_csv",
    "psd",
    "slip_from_speed",
]
