"""Diagnóstico automático para la GUI, usando el paquete ``clb_fft``.

No requiere instalar el paquete: se agrega ``<repo>/src`` al ``sys.path``,
así la GUI sigue abriéndose con un simple ``python gui.py``.
"""

import sys
from pathlib import Path

_SRC = Path(__file__).resolve().parents[1] / "src"
if _SRC.is_dir() and str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from clb_fft import (  # noqa: E402
    Motor,
    all_fault_frequencies,
    averaged_amplitude_spectrum,
    decimate,
    estimate_supply_frequency,
    load_hioki_archive,
    load_hioki_csv,
)
from clb_fft.diagnosis import assess_broken_bars, diagnose  # noqa: E402
from clb_fft.plotting import plot_spectrum  # noqa: E402
from clb_fft.quality import check_quality  # noqa: E402


def _cargar(ruta, canal):
    if ruta.lower().endswith((".tar", ".tar.xz", ".tar.gz", ".tgz")):
        return load_hioki_archive(ruta, channels=[canal])
    return load_hioki_csv(ruta, channels=[canal])


def _salida(ruta, sufijo):
    """Archivo de salida junto a la medición: Medicion.csv -> Medicion_<sufijo>."""
    p = Path(ruta)
    nombre = p.name
    for ext in (".tar.xz", ".tar.gz", ".tgz", ".tar", ".csv"):
        if nombre.lower().endswith(ext):
            nombre = nombre[: -len(ext)]
            break
    return p.with_name(f"{nombre}_{sufijo}")


def analizar(ruta, polos, velocidad_rpm, canal="I1", f_nominal=50.0):
    """Ejecuta el diagnóstico completo de una medición.

    Parámetros
    ----------
    ruta: primer archivo .csv de la exportación (o un .tar.xz).
    polos: número de polos del motor (ej. 4).
    velocidad_rpm: velocidad medida del eje [rpm].

    Retorna ``(texto_resumen, [figuras], tabla_de_picos)`` y guarda
    ``<medición>_DIAGNOSTICO.txt`` y ``<medición>_PICOS.csv`` junto a los datos.
    """
    if polos < 2 or polos % 2:
        raise ValueError("El número de polos debe ser par (2, 4, 6, ...).")
    rec = _cargar(str(ruta), canal)
    x, fs = rec[canal], rec.fs

    calidad = check_quality(rec, canal)
    f1 = estimate_supply_frequency(x, fs, nominal=f_nominal)
    n_sinc = 120 * f1 / polos
    if not 0 < velocidad_rpm < n_sinc:
        raise ValueError(
            f"Velocidad {velocidad_rpm:.0f} rpm fuera de rango: con {polos} polos y "
            f"f1 = {f1:.2f} Hz la velocidad síncrona es {n_sinc:.0f} rpm."
        )
    motor = Motor.from_speed(f1, polos // 2, velocidad_rpm)

    resolucion = max(0.05, calidad.recommended_resolution)
    xd, fsd = decimate(x, fs, 20)
    espectro = averaged_amplitude_spectrum(xd, fsd, resolution=resolucion)
    barras = assess_broken_bars(espectro, f1, motor.slip)
    tabla = diagnose(espectro, motor, fmin=1.0, fmax=min(10 * f1, fsd / 2), prominence_db=10)

    marcas = all_fault_frequencies(motor)
    nombre = Path(ruta).name
    figuras = [
        plot_spectrum(
            espectro,
            f1 - 15,
            f1 + 15,
            marcas,
            f"{nombre} - entorno de f1 (s = {motor.slip:.4f})",
        ),
        plot_spectrum(
            espectro,
            0,
            min(4 * f1, fsd / 2),
            marcas,
            f"{nombre} - espectro de amplitud (resolución {resolucion:.3f} Hz)",
        ),
    ]

    texto = "\n".join(
        [
            f"Archivo: {nombre}",
            f"Duración: {rec.duration:.1f} s   fs = {fs:g} Hz",
            "",
            "CALIDAD DE LA MEDICIÓN",
            calidad.summary(),
            "",
            f"MOTOR: {polos} polos, {velocidad_rpm:.0f} rpm, deslizamiento s = {motor.slip:.4f}",
            "",
            "BARRAS ROTÓRICAS",
            barras.summary(),
            "",
            f"Picos etiquetados: {(tabla.family != 'unclassified').sum()} de {len(tabla)}",
            "Interpretación final: ingeniero/a eléctrico/a.",
        ]
    )
    _salida(ruta, "DIAGNOSTICO.txt").write_text(texto, encoding="utf-8")
    tabla.to_csv(_salida(ruta, "PICOS.csv"), index=False)
    return texto, figuras, tabla
