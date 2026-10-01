# HOWTO: installing, building, and running FASER

This is a single, consolidated reference for the current (V10) repository
structure: a from-scratch install, the build, where simulation/reconstruction
data actually lives, and every option of every Python wrapper script. It
complements rather than replaces the existing docs — see "Further reading"
at the end for the deeper, topic-specific write-ups this points to.

The overall pipeline is: a GENIE or official FASERMC neutrino-interaction
sample is converted into the common `TPOEvent` format by **ConvertGENIE** or
**ConvertFASERMC**; **FASERG4** (`faserps`) runs that through Geant4 and
writes per-event `TcalEvent` truth files; **Batch** (`batchreco.exe`)
reconstructs those into `TPORecoEvent` files; and **EvDisplay** or the
analysis/TauSearch code consumes the result. `run_regression_tests.py` wraps
`faserps` and `batchreco.exe` end to end for an automated
simulate-reconstruct-compare check.

## 1. Install from scratch

```bash
git clone https://github.com/rubbiaa/FASER.git
cd FASER
source setup.sh
```

`setup.sh` auto-detects which known site you're on — André's Mac, the Ubuntu
(ryzen01) box, or lxplus — and sources the right `thisroot.sh`/`geant4.sh`
for that site, then hands off to `common_setup.sh` for everything that's
identical everywhere (`CLHEPINSTALL`/`RAVEINSTALL`/`GENFITINSTALL`,
`LD_LIBRARY_PATH`, `CMAKE_PREFIX_PATH`, `PATH`, and a sanity check that
reports any required variable that's missing or stale). It also defines a
shortcut, `fb` ("faser build"), that reruns `cmake --build $HOMEFASER/build
-j` from anywhere. On a machine it doesn't recognize, it prints exactly what
to add (an `elif` branch near the top of `setup.sh` that sources the right
ROOT/Geant4 setup scripts for that site).

For a completely new machine — installing ROOT and Geant4 themselves, OS
package prerequisites, and troubleshooting a clean build — see
[`docs/INSTALL.md`](INSTALL.md) and
[`docs/GEANT4_INSTALL.md`](GEANT4_INSTALL.md). This HOWTO assumes
ROOT and Geant4 are already installed somewhere `setup.sh` can find them.

## 2. Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j
```

The first build also fetches and compiles CLHEP, Rave, GenFit, googletest,
and Pythia8 as part of the same CMake superbuild, so it can take a while;
later builds are incremental. After `source setup.sh`, the `fb` shortcut
reruns the second line from anywhere — handy after a quick source edit.
Every built executable lands in `build/bin/` (`faserps`, `batchreco.exe`,
`batchreco_detresp.exe`, `dumphits.exe`, `evDisplay.exe`, `ConvertGENIE.exe`,
`Convert.exe`, `CombineFluxes.exe`, the gtest binaries, ...).

The most commonly useful `-D<OPTION>=<value>` flags at the `cmake -S . -B
build` step (full list and reuse-an-existing-dependency examples in
`docs/INSTALL.md`):

| Option | Default | Meaning |
|---|---|---|
| `FASER_BUILD_EXTERNALS` | `ON` | Build CLHEP/Rave/GenFit/googletest/Pythia8 from source |
| `FASER_BUILD_GEANT4_PACKAGES` | `ON` | Build FASERG4 and FASERCalProtoG4 (requires Geant4) |
| `FASER_BUILD_DISPLAY` | `ON` | Build the GenFit-based event display / muon spectrometer package |
| `FASER_BUILD_TAUSEARCH` | `ON` | Build TauSearch |
| `FASER_BUILD_TESTS` | `ON` | Build FASER's own gtest suite (`Tests/`), registered with CTest |

## 3. Where the data lives (`$FASERDATA`)

Every executable that reads or writes FASERG4/Batch data resolves the
location through `$FASERDATA` (`FASER::GetDataDir()`, in
`CoreUtils/FaserDataDir.hh`) instead of a hardcoded relative path or an
inter-directory symlink — `GetDataDir()` creates a requested subdirectory on
demand if it doesn't already exist. `common_setup.sh` defaults `FASERDATA`
to `$HOMEFASER/data` (gitignored) unless a per-site branch in `setup.sh` has
already exported it to somewhere else (scratch, EOS, a checkout dedicated to
one run/target, ...) before `common_setup.sh` runs.

| Subdirectory | Contents | Written by | Read by |
|---|---|---|---|
| `$FASERDATA/faserG4/` | Per-event `TcalEvent` truth files, `FASERG4-Tcalevent_<run>_<event>[_<mask>].root` | `faserps` (via `run_faserps.py`) | `batchreco.exe`, `dumphits.exe`, `evDisplay.exe` |
| `$FASERDATA/batch/` | Reconstructed `TPORecoEvent` files, `Batch-TPORecevent_<run>_<min>_<max>[_<mask>].root` | `batchreco.exe` (via `run_batchreco.py`) | Analysis code, `summarize_output.py` |
| `$FASERDATA/GENIE/` | GENIE-generated input samples (too large for git) | Fetched automatically — see §4, `fetch_data.py` | `faserps` (via `run_faserps.py`'s `--input-file`) |
| `$FASERDATA/GDML/` | The detector geometry FASERG4 exports on every run | `faserps` | `run_batchreco.py`'s default `--geometry-file` |

## 4. The Python wrapper scripts

### `run_faserps.py` — runs the Geant4 simulation

Builds the V10-geometry macro in memory and pipes it straight to `faserps`'
stdin — there's no `.mac` file to keep in sync by hand. Full prose
documentation (custom geometry parameters, interactive UI, MuonDIS physics)
is in `run_faserps.md`; every flag:

| Flag | Default | Meaning |
|---|---|---|
| `--build-dir` | `build/` | CMake build directory containing `bin/faserps` |
| `--vis` | off | Launch `faserps`' interactive Geant4 UI instead of piping a macro (no macro used) |
| `--print-macro` | off | Print the resolved macro text and exit, without running `faserps` |
| `--dry-run` | off | Print the command (and macro) but don't execute it |
| `--input-file` | `$FASERDATA/GENIE/FASERMC-PO-Run10000-0_53954_3DCAL.root` | Value for `/generator/rootinputfilename`. A relative custom value resolves against `FASERG4/` (faserps' cwd). Ignored with `--muons`/`--muondis` |
| `--start-event` | `0` | Value for `/generator/startevent`. Ignored with `--muons`/`--muondis` |
| `--n-events` | `100` | Value for `/run/beamOn` |
| `--muons` | off | Generate single fixed-momentum muons instead of reading `--input-file` |
| `--muon-momentum-gev` | `100.0` | Muon momentum in GeV (`--muons`/`--muondis` only) |
| `--muondis` | off | Enable MuonDIS (on-the-fly Pythia8 DIS per muon interaction); implies `--muons`. See `docs/README_MuonDIS.md` |
| `--muondis-cross-section-bias` | `150.0` | `/physics/muondis/crossSectionBias` (`--muondis` only; `1` = unbiased) |
| `--muondis-q2min` | `1.0` | `/physics/muondis/q2min` (GeV², `--muondis` only) |
| `--muondis-interaction-log` | `""` (disabled) | `/physics/muondis/interactionLog` CSV path (`--muondis` only) |
| `--muondis-pdf-set` | `""` (Pythia8's built-in proton PDF) | `/physics/muondis/pdfSet`, a path under `FASERG4/input/`, e.g. `input/NNPDF40_nnlo_as_01180_charmasy_0000.dat` (`--muondis` only) |
| `--muondis-xbjmin` | `0.0` (no cut) | `/physics/muondis/xbjmin` (`--muondis` only) |
| `--muondis-debug` | off | Add `/physics/muondis/debug true` (`--muondis` only) |
| `--tilt-deg` | `-4.5` | `/FASER/tiltY` (degrees) |
| `--shift-x-cm` | `45.0` | `/FASER/LOS/shiftX` (cm) |
| `--shift-y-cm` | `24.0` | `/FASER/LOS/shiftY` (cm) |

The default `--input-file` is fetched automatically from a public CERNBox
link the first time it's needed (sha256-verified), so a fresh `git clone`
works with no separate download step.

### `run_batchreco.py` — runs the reconstruction

Wraps `Batch/batchreco.exe` with named, validated flags instead of
positional arguments. Full prose documentation in `run_batchreco.md`; every
flag:

| Flag | Default | Meaning |
|---|---|---|
| `--run` | *(required)* | Run number — deliberately no default, since silently reconstructing the wrong run is worse than being forced to say which one |
| `--min-event` | `0` | First event index |
| `--max-event` | `99999999` | Event index to stop before (exclusive), matching `BatchReco.cc`'s own loop |
| `--mask` | none (all events) | Process only this interaction type: `nueCC`, `numuCC`, `nutauCC`, `nuNC`, or `nuES` |
| `--multi-thread` | off | Pass `-mt` to enable multi-threading |
| `--geometry-file` | `$FASERDATA/GDML/FASERCAL_V10.gdml` | GDML geometry to load |
| `--build-dir` | `build/` | CMake build directory containing `bin/batchreco.exe` |
| `--split` | `1` | Split `[min-event, max-event)` into this many contiguous chunks, each run as a separate background job with its own log under `data/batch/logs/` |
| `--dry-run` | off | Print the command(s) that would run, but don't execute them |

### `run_regression_tests.py` — simulate → reconstruct → compare

Runs the real pipeline end to end for a handful of cases (the default
neutrino sample reconstructed four ways — `nueCC`/`numuCC`/`nutauCC`/`nuNC`
— plus `--muons` and `--muondis`) and compares aggregate output statistics
against a committed golden JSON baseline. Full design in
`docs/REGRESSION_TESTS.md`; every flag:

| Flag | Default | Meaning |
|---|---|---|
| `--case` | all cases | Restrict to this golden case name (repeatable); see `--list` for names |
| `--record` | off | Write the golden file(s) for the selected case(s) instead of comparing |
| `--rel-tol` | `1e-9` | Relative tolerance for floating-point aggregate comparisons |
| `--build-dir` | `build/` | CMake build directory |
| `--python` | the interpreter running this script | Python interpreter used for the `run_faserps.py`/`run_batchreco.py`/`summarize_output.py` subprocesses — point this at a PyROOT-enabled interpreter if it differs from the default one |
| `--list` | — | List available case names and exit |

```bash
python3 run_regression_tests.py --record   # first time, or after a deliberate behavior change
python3 run_regression_tests.py            # every other time: compare against golden/
```

### `fetch_data.py` — downloads large input files from CERNBox

A small, growable manifest of public-CERNBox-link downloads into
`$FASERDATA`, sha256-verified, used automatically by `run_faserps.py` for
the default GENIE sample.

| Flag | Default | Meaning |
|---|---|---|
| `--force` | off | Re-download even if already present |
| `--dry-run` | off | Show what would be fetched, don't fetch |
| `--list` | — | List manifest entries and exit |

### `Tests/regression/summarize_output.py` — PyROOT output summarizer

Invoked internally by `run_regression_tests.py`'s comparison step (and
callable standalone for debugging a single case by hand).

| Flag | Default | Meaning |
|---|---|---|
| `--dict-path` | CoreUtils' built ROOT dictionary | Path to `CoreUtils`' built ROOT dictionary |
| `--truth-dir` | none (skip truth summary) | Directory of `FASERG4-Tcalevent_<run>_<event>.root` files, typically `$FASERDATA/faserG4` — requires `--run` and `--n-events` |
| `--run` | none | Run number, required with `--truth-dir` |
| `--n-events` | none | Number of event indices (0..N-1) to look for, required with `--truth-dir` |
| `--reco-file` | none (skip reco summary) | Path to one `Batch-TPORecevent_*.root` file |

At least one of `--truth-dir` or `--reco-file` must be given.

## 5. A typical end-to-end run

```bash
source setup.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo && cmake --build build -j

python3 run_faserps.py --n-events 500                          # simulate (auto-fetches the default sample)
python3 run_batchreco.py --run 10000 --mask numuCC              # reconstruct just the numuCC events
build/bin/evDisplay.exe 10000 numuCC                             # look at the result
```

or, to exercise the whole pipeline automatically and check it against the
committed baseline:

```bash
python3 run_regression_tests.py
```

## 6. Tests

FASER's own C++ gtest suite (`MuonSpectrometerFieldTest`,
`MagnetGeometryProbeTest`, ...) builds automatically whenever
`FASER_BUILD_TESTS` and `FASER_BUILD_GEANT4_PACKAGES` are both `ON` (the
default), and runs with:

```bash
ctest --test-dir build --output-on-failure
```

The heavier simulate-reconstruct-compare regression suite is
`run_regression_tests.py`, described above and in `docs/REGRESSION_TESTS.md`.

## Further reading

- [`docs/INSTALL.md`](INSTALL.md) — the full CMake build: prerequisites, every option, troubleshooting a clean build.
- [`docs/GEANT4_INSTALL.md`](GEANT4_INSTALL.md) — building/installing Geant4 itself.
- [`run_faserps.md`](../run_faserps.md) / [`run_batchreco.md`](../run_batchreco.md) — prose documentation for those two wrappers.
- [`docs/README_MuonDIS.md`](README_MuonDIS.md) — MuonDIS physics and every `/physics/muondis/...` option.
- [`docs/REGRESSION_TESTS.md`](REGRESSION_TESTS.md) — the simulate-reconstruct-compare suite's full design, what's verified vs. not, and open items.
