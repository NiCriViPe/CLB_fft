"""Command line interface: ``clb-fft analyze <file> [options]``."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

from . import __version__


def _load(path: str, channel: str, fs: float | None):
    from .hioki import load_hioki_archive, load_hioki_csv

    if path.endswith((".tar", ".tar.xz", ".tar.gz", ".tgz")):
        return load_hioki_archive(path, channels=[channel], fs=fs)
    return load_hioki_csv(path, channels=[channel], fs=fs)


def cmd_quality(args) -> int:
    from .quality import check_quality

    rec = _load(args.file, args.channel, args.fs)
    print(f"{rec.n_samples} samples, fs = {rec.fs:g} Hz, duration = {rec.duration:.1f} s")
    print(check_quality(rec, args.channel).summary())
    return 0


def cmd_analyze(args) -> int:
    matplotlib.use("Agg")
    from .diagnosis import assess_broken_bars, diagnose
    from .faults import Motor, all_fault_frequencies
    from .plotting import plot_frequency_track, plot_spectrum
    from .quality import check_quality, frequency_track
    from .spectrum import averaged_amplitude_spectrum, decimate, estimate_supply_frequency, psd

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rec = _load(args.file, args.channel, args.fs)
    x, fs = rec[args.channel], rec.fs
    q = check_quality(rec, args.channel)
    f1 = estimate_supply_frequency(x, fs, nominal=args.f_nominal)
    print(q.summary())
    print(f"estimated f1           : {f1:.3f} Hz")

    if args.speed is not None:
        motor = Motor.from_speed(f1, args.pole_pairs, args.speed, rotor_bars=args.rotor_bars)
    else:
        motor = Motor(f1, args.pole_pairs, args.slip, rotor_bars=args.rotor_bars)

    res = max(args.resolution, q.recommended_resolution)
    # Low-frequency band: decimate (with anti-alias filter) for speed.
    xd, fsd = decimate(x, fs, args.decimate) if args.decimate > 1 else (x, fs)
    lf = averaged_amplitude_spectrum(xd, fsd, resolution=res)
    hf = psd(x, fs, resolution=1.0)

    marks = all_fault_frequencies(motor)
    name = Path(args.file).name
    figs = {
        "spectrum_low.png": plot_spectrum(
            lf, 0, min(4 * f1, fsd / 2), marks, f"{name} - amplitude spectrum (res {res:.3f} Hz)"
        ),
        "spectrum_fundamental.png": plot_spectrum(
            lf, f1 - 15, f1 + 15, marks, f"{name} - around f1 (s = {motor.slip:.4f})"
        ),
        "psd_full.png": plot_spectrum(hf, 0, fs / 2, (), f"{name} - Welch PSD", ylim=(-140, 5)),
        "frequency_track.png": plot_frequency_track(
            frequency_track(x, fs, 1.0), 1.0, f"{name} - supply frequency"
        ),
    }
    for fname, fig in figs.items():
        fig.savefig(out / fname, dpi=120)

    table = diagnose(lf, motor, fmin=1.0, fmax=min(10 * f1, fsd / 2), prominence_db=args.prominence)
    table.to_csv(out / "peaks.csv", index=False)
    brb = assess_broken_bars(lf, f1, motor.slip)
    print(brb.summary())
    summary = {
        "file": args.file,
        "fs": fs,
        "duration_s": rec.duration,
        "f1": f1,
        "slip": motor.slip,
        "resolution_hz": res,
        "quality": q.__dict__,
        "broken_bar": brb.__dict__,
        "labelled_peaks": table[table.family != "unclassified"].to_dict("records"),
    }
    (out / "summary.json").write_text(json.dumps(summary, indent=2, default=float))
    print(f"results written to {out}/")
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="clb-fft", description="Motor current signature analysis of HIOKI PQ ONE exports"
    )
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("file", help="first CSV of the export, or a .tar.xz archive")
        sp.add_argument("--channel", default="I1")
        sp.add_argument("--fs", type=float, default=None, help="sampling rate [Hz] (default: inferred)")

    q = sub.add_parser("quality", help="check record continuity and supply-frequency drift")
    common(q)
    q.set_defaults(func=cmd_quality)

    a = sub.add_parser("analyze", help="spectra, fault-frequency labelling and broken-bar indicator")
    common(a)
    a.add_argument("--pole-pairs", type=int, required=True)
    g = a.add_mutually_exclusive_group(required=True)
    g.add_argument("--speed", type=float, help="shaft speed [rpm]")
    g.add_argument("--slip", type=float, help="slip [p.u.]")
    a.add_argument("--rotor-bars", type=int, default=None)
    a.add_argument("--f-nominal", type=float, default=50.0)
    a.add_argument("--resolution", type=float, default=0.05, help="target resolution [Hz]")
    a.add_argument("--decimate", type=int, default=20, help="decimation factor for the low band")
    a.add_argument("--prominence", type=float, default=10.0, help="peak prominence [dB]")
    a.add_argument("--out", default="results")
    a.set_defaults(func=cmd_analyze)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
