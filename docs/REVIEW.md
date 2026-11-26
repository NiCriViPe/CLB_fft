# Repository review and roadmap

*Review date: November 2026 · scope: the whole repository as of commit `b990113`.*

This document has three parts:

1. **Findings**: what was wrong or risky, by severity, and what was done about it.
2. **Best-practice changes** made to the repository.
3. **Roadmap**: how to extend the work so it is interesting to electrical
   engineers, to academia, and to software/industry reviewers.

Line references (`file:line`) point to the code **before** this review
(`git show b990113:<file>`).

---

## 1. Findings

### 1.1 Signal-processing correctness (high impact: these change diagnoses)

| # | Finding | Where | Status |
|---|---|---|---|
| S1 | **Amplitude spectra plotted with `10·log10`.** An FFT magnitude is an amplitude, so dB is `20·log10`. With `10·log10` every level looks half as far below f1 as it really is: a −60 dB sideband reads −30 dB. This directly breaks severity thresholds such as Thomson's broken-bar table, which is defined in amplitude dB. | `POST_proceso.py:236,344` | Fixed (legacy + new core) |
| S2 | **Decimation without an anti-aliasing filter.** `a[range(0,N,20)]` keeps every 20th sample, so every component above 500 Hz (harmonics, slot harmonics, inverter switching) folds into 0–500 Hz as false peaks. A unit test reproduces this: a 990 Hz tone appears at 10 Hz with the naive method and is suppressed by 60+ dB with `scipy.signal.decimate`. | `POST_proceso.py:176` | Fixed |
| S3 | **No window (rectangular).** Leakage from the fundamental (≈ 100 A) masks sidebands 40–60 dB below it, which is exactly where early broken-bar and eccentricity signatures sit. A unit test shows a −60 dB sideband 2 Hz from an off-bin fundamental reading as −40 dB with a rectangular window and correctly with Hann. | `POST_proceso.py:192` | Fixed (Hann; flattop available for amplitude accuracy) |
| S4 | **Frequency axis built with `linspace(0, fs/2, N//2)`.** That spacing is `fs/2/(N/2−1)`, not `fs/N`, so the axis is slightly stretched. It is tiny for 4 M points but wrong by construction. `numpy.fft.rfftfreq` is exact. | `POST_proceso.py:193,302` | Fixed |
| S5 | **Low-frequency FFT normalised by `len(xf)`** instead of the maximum, so the dB axis has no physical reference. | `POST_proceso.py:196-198` | Fixed (0 dB = fundamental) |
| S6 | **The "continuous" waveform is 1000 event records of 0.28 s.** In the sample measurement, 151 of the 999 joints are discontinuous. Each discontinuity is a small impulse, and because they recur about every 0.28 s they can produce spurious lines spaced about 3.6 Hz apart. Nothing in the original code knows the records exist. | data | New `clb_fft.quality` detects and counts them |
| S7 | **Supply-frequency drift sets the real resolution.** f1 wanders 49.89 → 50.08 Hz during the 280 s capture. A single 280 s FFT promises 0.0036 Hz resolution but smears the fundamental over ≈ 50 bins, producing a broad "skirt" where sidebands should be read (see `docs/img/before_after_f1.png`). | `POST_proceso.py:157` (4 M-sample FFT) | New Welch-averaged spectrum; quality report recommends a resolution ≥ drift. Order tracking proposed (roadmap A1). |
| S8 | **False-alarm risk in automatic broken-bar grading.** When sidebands at f1(1±2s) are not resolvable from f1, reading the maximum near them returns the fundamental's own skirt. A first version of the new indicator fell into this trap on the real data: it reported "two bars broken" for a machine that shows no sidebands. | new code | Fixed: the indicator now requires resolvable sidebands **and** a real local peak, and otherwise reports "not resolvable" or "no sideband detected (> X dB)". Covered by regression tests. |
| S9 | `FreqGen.broken_bar_frequencies` ignored its `fgrid` argument (hard-coded 50 Hz) and returned 2(N−1) components instead of 2N. | `FreqGen.py:30-33` | Fixed |
| S10 | The notebook's bearing frequencies use the **supply** frequency (60 Hz) as shaft speed and an undefined variable `n` instead of `nb`. The diagnostic loop appends a peak to `unclassified` once for **each** family it does not match. | `Inline_testingV2.ipynb` | Re-implemented correctly in `clb_fft.faults` / `clb_fft.diagnosis` |
| S11 | Sampling rate hard-coded (`Fs = 20e3`). | `POST_proceso.py:75` | New loader infers it from the 1 ms timestamps (20 rows/ms → 20 kHz) |

### 1.2 Robustness and usability

| # | Finding | Status |
|---|---|---|
| R1 | Cancelling any file dialog passes `''` and crashes; cancelling the "N° de archivos" dialog passes `None` to `range()`. | Fixed (legacy) |
| R2 | "FFT baja/alta frecuencia" crash with `FileNotFoundError` unless "Ejecutar" was pressed first (they read `_y-t_DATA.csv`); only `PSD` checked. | Fixed (legacy) |
| R3 | Multi-part search stops at `_009`; a long PQ ONE export can have more parts. | Fixed (legacy + new loader) |
| R4 | `welch(nperseg=sy/Ntrunc)` passes a float; recent SciPy requires an int. | Fixed |
| R5 | `s.replace('.csv', ...)` rewrites any `.csv` in directory names; paths split on `"/"`; `aux[-4]` raises `IndexError` for shallow paths. | New code uses `pathlib`; legacy unchanged (low risk on Windows) |
| R6 | Heavy work runs on the Tk main thread, so the window freezes for tens of seconds. `show(block=True)` in `PSD` blocks the GUI. | Open (roadmap B4) |
| R7 | An intermediate `_y-t_DATA.csv` (time + current, about 150 MB) is written next to the data on every run. | Open in legacy; the new core works in memory |
| R8 | The GUI writes "PSD [dB/Hz]" on a spectrum normalised to its maximum (dB relative). | Fixed (label) |

### 1.3 Repository hygiene

| # | Finding | Status |
|---|---|---|
| H1 | `__pycache__/*.pyc` committed even though `.gitignore` lists one of them. | Untracked; `.gitignore` now covers `__pycache__/` |
| H2 | `APP.lnk` is a Windows shortcut to `C:\Users\<user>\...\gui.py`: not portable, and it leaks a local path. | Untracked; `*.lnk` ignored |
| H3 | No dependency declaration, tests, CI, or license. | `pyproject.toml`, 37 tests, GitHub Actions added. **License: your decision** (MIT or BSD-3 are usual for portfolio code). |
| H4 | About 40 % of `POST_proceso.py` is commented-out code; `low_freq_fft` and `high_freq_fft` are ~95 % duplicated. | Superseded by the `clb_fft` package; delete the dead code once the GUI calls the package (B1). |
| H5 | Mixed Spanish/English identifiers and docs. | README is now in English with a Spanish summary. Pick one language for code identifiers. |
| H6 | **Data governance.** `Measurements/NH-I/` and the README name a real plant and operator. | **Action for you:** confirm you are allowed to publish this measurement. If not, anonymise it (asset tag, plant name) or replace it with a synthetic/public dataset. |

---

## 2. What was added

```
src/clb_fft/
  hioki.py      PQ ONE loader: any number of parts, BOM, fs inference, .tar.xz without extracting
  spectrum.py   calibrated amplitude spectrum, Welch amplitude/PSD, anti-aliased decimation,
                sub-bin supply-frequency estimate
  faults.py     fault frequencies with references: broken bars, mixed / static / dynamic
                eccentricity (rotor slot harmonics), shorted turns, bearings (BPFO/BPFI/BSF/FTF)
  diagnosis.py  peak detection on prominence, labelling, broken-bar indicator with false-alarm guards
  quality.py    record-joint discontinuities, per-record f1 tracking, resolution recommendation
  plotting.py   figure builders that never call show()
  cli.py        `clb-fft quality|analyze`  -> PNGs + peaks.csv + summary.json
tests/          37 tests on synthetic signals with known answers
examples/method_comparison.py   regenerates docs/img/*
.github/workflows/ci.yml        ruff + pytest on 3.10/3.12/3.13 + CLI smoke test on the real data
```

Design choices worth defending in an interview:

* **Tests use physics, not snapshots.** A sine of amplitude A must read A; a
  sideband injected at −45 dB must be graded "moderate"; a drifting healthy
  supply must **not** raise an alarm.
* **The library never shows figures or opens dialogs.** The GUI, CLI and
  notebooks are thin shells, so everything can be tested headlessly.
* **The indicator can say "I can't tell".** "Not resolvable" is a valid
  output, and that is what separates a diagnostic tool from a plotting script.

---

## 3. Roadmap: adding value for three audiences

Items are ordered by value per effort within each track.

### A. Electrical engineering / condition monitoring

| ID | Proposal | Why it matters |
|---|---|---|
| A1 | **Supply-synchronous resampling (order tracking).** Use the per-record f1 track (already computed) to resample the current so f1 is constant, then FFT. | Removes the drift skirt in `before_after_f1.png` and recovers the real resolution. Needed for large, low-slip machines (s ≈ 0.3–0.8 %, sidebands 0.3–0.8 Hz from f1). |
| A2 | **Slip estimation without a tachometer** from rotor slot harmonics (`faults.rotor_slot_harmonics` already generates them) or from f1 ± fr. | Field measurements rarely include speed, and every fault frequency depends on slip. |
| A3 | **Use all three phases.** Park's vector / extended Park's vector approach (EPVA), and negative-sequence current. | Distinguishes stator asymmetry and supply unbalance from rotor faults, which a single-phase analysis cannot do. |
| A4 | **Use the voltage channel (U1).** Subtract supply-borne harmonics (present in U) from current harmonics. Compute instantaneous power spectra. | Separates faults from grid pollution. Power-spectrum MCSA shows the fr component directly. |
| A5 | **Load normalisation and trending.** Store (date, load %, slip, sideband dB) per asset, and alarm on change from that asset's baseline rather than on absolute thresholds. | Matches IEEE 1415 / ISO 20958 practice. Absolute thresholds are load-dependent. |
| A6 | **Envelope / Hilbert and spectral kurtosis** for bearings, and a vibration channel where available (the file is called *MotorConVibracion*). | Bearing signatures in current are weak; envelope methods improve SNR. |
| A7 | **Inverter-fed motors**: variable f1, PWM switching sidebands. Already partly supported via `--f-nominal`. | Most new installations are VFD-fed. |

### B. Software engineering / industry portfolio

| ID | Proposal | Why it matters |
|---|---|---|
| B1 | **Rewire the Tk GUI onto `clb_fft`** (started: the "Diagnóstico" button uses it) and delete the duplicated legacy functions. | Removes two parallel implementations. |
| B2 | **Fast storage**: convert CSV to Parquet/NPZ once (about 10× smaller and faster), with metadata (asset, load, speed) in a sidecar JSON. | The 250 MB of CSV per measurement is the main bottleneck. |
| B3 | **Report generation**: one HTML/PDF per measurement with the spectra, the labelled peak table, the quality report and the indicators. | This is what a maintenance planner actually consumes. |
| B4 | **Responsive GUI**: run jobs in a worker thread or process, add a progress bar, keep the plots interactive. Or move to a small web UI (Streamlit/Panel) that is easier to show in a portfolio. | Fixes R6, and works without a Windows install. |
| B5 | **Packaging and releases**: semantic versions, a CHANGELOG, a PyPI release, `pre-commit` with ruff, `mypy --strict` on `src/`. | Shows engineering discipline at a glance. |
| B6 | **Asset database and trend dashboard**: SQLite plus a dashboard of sideband dB vs. time per motor. | Turns one-off analyses into condition monitoring (A5). |
| B7 | **Documentation site** (MkDocs + mkdocstrings) with a theory page per fault family. | Reaches both audiences at once. |

### C. Academic / research

| ID | Proposal | Why it matters |
|---|---|---|
| C1 | **Benchmark on public datasets** (for example the IEEE DataPort broken-rotor-bar test-bench database from the University of São Paulo, or the Paderborn University bearing dataset, which includes motor currents). Report detection rates against load and severity. | Turns a tool into a result; makes the work citable and reproducible. |
| C2 | **Synthetic signal generator** with controllable faults, drift, noise and record gaps (the test helpers are a seed). | Ground truth for every method; supports ablation studies. |
| C3 | **Classical vs. learned features**: compare the rule-based indicator with a small classifier (gradient boosting on sideband features, or a 1-D CNN on order-tracked spectra) under train/test splits by machine. | A clear research question with honest baselines. |
| C4 | **High-resolution estimators**: zoom-FFT, MUSIC/ESPRIT, or matrix pencil for sidebands closer than the FFT resolution. | Directly addresses the low-slip / short-record limitation shown here. |
| C5 | **Effect of record stitching**: quantify how PQ-analyser event-record joints bias MCSA, and propose a joint-aware estimator. | A small, novel, publishable contribution specific to this data source. Conference paper scale (e.g. IEEE SDEMPED, ICEM). |
| C6 | **Uncertainty**: confidence intervals for sideband levels from Welch segment variance. | "−48 ± 2 dB" is more defensible than "−48 dB". |

### Suggested next three steps

1. **A1 + A2** (order tracking and slip estimation): highest engineering value, and they reuse code that now exists.
2. **B1 + B3** (GUI on the new core, plus an HTML report): makes the improvements visible to users.
3. **C1** (benchmark on a public dataset): makes the project credible outside the plant.
