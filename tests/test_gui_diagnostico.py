"""The GUI's 'Diagnóstico' button logic, tested without opening a window."""

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pytest  # noqa: E402
from test_hioki import write_export  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "Predictive_maintenance_module"))
import diagnostico  # noqa: E402


@pytest.fixture(scope="module")
def medicion(tmp_path_factory):
    folder = tmp_path_factory.mktemp("med")
    return write_export(folder, n_parts=3, records_per_part=100)[0]  # 6 s at 20 kHz


def test_analizar_end_to_end(medicion):
    texto, figuras, tabla = diagnostico.analizar(medicion, polos=4, velocidad_rpm=1470)
    assert "CALIDAD DE LA MEDICIÓN" in texto and "BARRAS ROTÓRICAS" in texto
    assert "s = 0.0200" in texto
    assert len(figuras) == 2
    assert (tabla.family == "grid").any()
    assert (medicion.parent / "Meas_DIAGNOSTICO.txt").is_file()
    assert (medicion.parent / "Meas_PICOS.csv").is_file()


@pytest.mark.parametrize("polos, rpm", [(3, 1470), (4, 1600)])
def test_analizar_rejects_bad_motor_data(medicion, polos, rpm):
    with pytest.raises(ValueError):
        diagnostico.analizar(medicion, polos=polos, velocidad_rpm=rpm)
