# FASERCAL - electronic calorimeter for FASER Run4

FASERCAL code to simulate and analyse events in the FASERCAL detector.

# Event data flow

This project follows a structured workflow for simulating and processing neutrino interactions within the FASERCAL experiment using GEANT4. The process begins with GENIE-generated neutrino-interaction samples, converted into the TPOEvent class format by **ConvertGENIE** (official FASERMC Monte Carlo files can be converted the same way with **ConvertFASERMC** instead). Next, the **FASERG4** module uses GEANT4 to simulate events from the TPOEvent class, generating TcalEvent objects. These simulated events can be visualized using the **EvDisplay** module or processed through the **Batch** module for batch reconstruction, providing a comprehensive analysis pipeline for neutrino interaction events.

![Diagram of the project](images/eventchainflow.png)

# Quick start

The full build is CMake-based; see `docs/INSTALL.md` for prerequisites, build
options and troubleshooting. This section just covers the everyday
build/run/test loop once your machine is already set up.

## Build

```bash
git clone https://github.com/rubbiaa/FASER.git
cd FASER
source setup.sh                                      # sets up ROOT/Geant4/Pythia8 for known sites
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j
```

The first build also fetches and compiles CLHEP, Rave, GenFit, googletest
and Pythia8 (see `docs/INSTALL.md`), so it can take a while; later builds are
incremental. `source setup.sh` also defines a shortcut, `fb`, that reruns
`cmake --build $HOMEFASER/build -j` from anywhere - useful after a quick
source edit.

All executables land in `build/bin/` (`faserps`, `batchreco.exe`,
`evDisplay.exe`, `ConvertGENIE.exe`, `Convert.exe`, ...).

## Where simulation/reconstruction data lives

Every executable that reads or writes FASERG4/batchreco data resolves the
location through `$FASERDATA` (`FASER::GetDataDir()`,
`CoreUtils/FaserDataDir.hh`) rather than a hardcoded path or symlink.
`setup.sh`/`common_setup.sh` default `FASERDATA` to `$HOMEFASER/data`
(gitignored) and create it automatically; a specific site or checkout can
point it elsewhere (scratch, EOS, ...) by exporting it before
`common_setup.sh` runs - see the comments in `setup.sh`. `FASERG4` writes
under `$FASERDATA/faserG4/`, `Batch` reads that and writes under
`$FASERDATA/batch/`. Two more subdirectories hold input/derived data too
large or too generated-on-the-fly to check into git: `$FASERDATA/GENIE/`
(the default GENIE-generated input sample) and `$FASERDATA/GDML/` (the
detector geometry FASERG4 exports every run, via `FASER::GetDataDir("GDML")`
- read back by `run_batchreco.py`). The GENIE sample is fetched
automatically from a public CERNBox link the first time it's needed - see
`fetch_data.py` and "Run faserps" below.

## Run faserps (the GEANT4 simulation)

```bash
python3 run_faserps.py                       # 100 events, V10 geometry, the current default sample
python3 run_faserps.py --n-events 500 --input-file some_other_sample.root
python3 run_faserps.py --muons --n-events 1000 --muon-momentum-gev 250
```

`run_faserps.py` builds the V10-geometry macro in memory and pipes it
straight to `faserps`' stdin - there's no `.mac` file to keep in sync by
hand. See `run_faserps.md` for every option (custom geometry parameters,
`--vis` for the interactive Geant4 UI, `--print-macro`/`--dry-run`, muon
mode, MuonDIS, ...).

The default input sample is too large to commit to git; if it's missing
under `$FASERDATA/GENIE/`, `run_faserps.py` fetches it automatically from
a public CERNBox link before running (sha256-verified - see
`fetch_data.py`), so the commands above work right after a fresh
`git clone` with no separate download step. Run `python3 fetch_data.py`
directly to fetch it (or anything else in its manifest) ahead of time
instead, e.g. before going offline.

## MuonDIS

MuonDIS replaces the primary background muon's Geant4 nuclear interaction
with an on-the-fly, Pythia8-driven muon-nucleon deep-inelastic-scattering
event, generated per interaction from the muon's actual energy and
direction. It's disabled by default; enable it with:

```bash
python3 run_faserps.py --muondis --n-events 1000
```

See `run_faserps.md` (MuonDIS mode) for the `--muondis-*` options, and
[`docs/README_MuonDIS.md`](docs/README_MuonDIS.md) for the physics
details (target treatment, PDF choice, truth-level output added to
`TPOEvent`, every `/physics/muondis/...` command).

## Run batchreco (reconstruction)

```bash
python3 run_batchreco.py --run 10000                          # reconstruct all of run 10000
python3 run_batchreco.py --run 10000 --mask numuCC             # only numuCC events
python3 run_batchreco.py --run 10000 --max-event 12000 --split 6   # split into 6 background jobs
```

`run_batchreco.py` wraps `Batch/batchreco.exe` with named, validated
flags instead of positional arguments. See `run_batchreco.md` for the full
option list (event ranges, `--multi-thread`, a non-default `--geometry-file`,
splitting a run into parallel background jobs, ...).

## Run the tests (gtests)

FASER has its own small C++ regression test suite under `Tests/`, built on
Google Test (fetched automatically as part of the superbuild, see
`docs/INSTALL.md`) and registered with CTest:

- **MuonSpectrometerFieldTest** - cross-checks that the GEANT4
  simulation's muon-spectrometer magnetic field
  (`FASERG4/src/MuonDetMagneticField.cc`) and GenFit reconstruction's
  field (`CoreUtils/GenMagneticField.hh`) actually agree. Both are thin
  adapters around the one shared field model in
  `CoreUtils/MuonSpectrometerField.hh`; this test is what proves the two
  adapters wire it up identically, rather than just trusting that by
  inspection.
- **MagnetGeometryProbeTest** - unit tests for `FASER::ProbeMagnetSlit`
  (`CoreUtils/MagnetGeometryProbe.hh`), the geometry/field consistency
  probe: builds small ROOT `TGeo` shapes by hand to prove the probe both
  agrees with the real MDT magnet's numbers and actually detects a
  genuine mismatch when block and slit are deliberately built
  inconsistent.

They're built automatically whenever `FASER_BUILD_TESTS` and
`FASER_BUILD_GEANT4_PACKAGES` are both `ON` (the default for both). Run
them with:

```bash
ctest --test-dir build --output-on-failure
```

or build/run one directly, e.g. `build/bin/MuonSpectrometerFieldTest`.
Run these after touching either magnetic-field adapter, the shared field
model, or the magnet geometry probe.

## Run the regression tests (simulate -> reconstruct -> compare)

A second, heavier kind of test alongside the gtests above: `run_regression_tests.py`
runs the *real* pipeline end to end (`faserps` then `batchreco.exe`) for a
handful of event types - the default neutrino sample reconstructed four
ways (`nueCC`/`numuCC`/`nutauCC`/`nuNC`), plus `--muons` and `--muondis` -
and compares aggregate output statistics against a committed golden JSON
baseline (`Tests/regression/golden/*.json`), relying on `faserps.cc`'s
hardcoded random seed for exact reproducibility rather than a statistical
tolerance.

```bash
python3 run_regression_tests.py --record   # first time, or after a deliberate behavior change
python3 run_regression_tests.py            # every other time: compare against golden/
```

See [`docs/REGRESSION_TESTS.md`](docs/REGRESSION_TESTS.md) for the full
design, what's verified from source vs. not yet exercised against a real
build, and open items (CI isn't wired up to run this yet).

# EvDisplay

Interactive event display of FASERG4 output.

Usage: `evDisplay.exe [-g <geometryfile>] [-r] <run> [mask]`
   <run>                     Run number
   mask                      To process only specific events (def=none):   nueCC, numuCC, nutauCC, or nuNC

It reads from `$FASERDATA/faserG4/` like every other tool above - no
`input` symlink to set up first.

- run the event display

   to display nueCC events
   ```bash
   $ build/bin/evDisplay.exe 200026 nueCC
   ````

   to display numuCC events
   ```bash
   $ build/bin/evDisplay.exe 200025 numuCC
   ````

   to display nutauCC events
   ```bash
   $ build/bin/evDisplay.exe 200035 nutauCC
   ````

![Diagram of the project](images/numuCC_ev1.jpg)

# DumpHits (in Batch directory)

Very simple app to read and dump all hits from events.

Usage: `dumphits.exe <run> [maxevent] [mask]`
   <run>                     Run number
   maxevent                  Maximum number of events to process (def=-1)
   mask                      To process only specific events (def=none):   nueCC, numuCC, nutauCC, nuNC or nuES

for example to get all the hits of nueCC events from the kaon decay flux:

   ```bash
   $ build/bin/dumphits.exe 200026 10 nueCC > dump.log
   ```

# ConvertFASERMC

Converts official FASER MC files into FASERCAL PO files (generator level),
producing the TPOEvent-format input `faserps`/FASERG4 consumes. Built as
`Convert.exe`.

# ConvertGENIE

Converts GENIE-generated event files into the same TPOEvent format, as an
alternative to ConvertFASERMC for GENIE-based samples. Built as
`ConvertGENIE.exe`; `CombineFluxes.exe` (same directory) combines flux
files upstream of it.

# TauSearch
A generator level tau search analysis code

- t.C : code to convert FASER ntuple into event summary tuples
- s.C : analyse event summary tuples for each tau decay channel and create sig/background tuples
- a.C : read sig/bkg tuples for each decay channel and perform BDT analysis

# Installation

See `docs/INSTALL.md` for the full CMake build: prerequisites, build options,
and troubleshooting a clean build. In short:

```bash
git clone https://github.com/rubbiaa/FASER.git
cd FASER
source setup.sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j
```

`setup.sh` auto-detects which known site you're on (including lxplus) and
sets up ROOT/Geant4/Pythia8 accordingly, so the same command works
everywhere; on a machine it doesn't recognize, it prints what to do (add
an `elif` branch for your site - see the comments at the top of the
script).

## Further documentation

Longer, topic-specific write-ups (build/installation details, physics
reviews, one-off studies, ...) live under [`docs/`](docs/) rather than
cluttering the repo root or the subdirectory they're about:

- [`docs/INSTALL.md`](docs/INSTALL.md) - the full CMake build: prerequisites, options, troubleshooting.
- [`docs/GEANT4_INSTALL.md`](docs/GEANT4_INSTALL.md) - building/installing Geant4 itself.
- [`docs/README_MuonDIS.md`](docs/README_MuonDIS.md) - MuonDIS physics and every `/physics/muondis/...` option (see "MuonDIS" above for the quick start).
- [`docs/HYPERON_DECAY_REVIEW.md`](docs/HYPERON_DECAY_REVIEW.md) - why long-lived hyperons have their Pythia8 decay switched off (see `CoreUtils/TPOEvent.cc`'s `initialize_pythia()`).
- [`docs/NEUTRON_SPECTRUM_RUN.md`](docs/NEUTRON_SPECTRUM_RUN.md) - notes on a neutron-spectrum run with `FASERCalProtoG4`.
- [`docs/MuonSpectrometerReport.md`](docs/MuonSpectrometerReport.md) - the muon-spectrometer magnetic field consistency work (shared Geant4/GenFit field model, geometry probe, gtests).
- [`docs/REGRESSION_TESTS.md`](docs/REGRESSION_TESTS.md) - the simulate->reconstruct->compare regression suite (see "Run the regression tests" above).

`run_faserps.md` and `run_batchreco.md` stay at the repo root, next to the
scripts they document, so script and doc can't drift out of sync with each
other.

 # Event masks

 - event masks are used to select only a type of events when running a job

 - the currently available event masks are:

   nueCC - nue charged currents
 
   numuCC  - numu charged currents
 
   nutauCC - nutau charged currents
 
   nuNC - all neutrinos neutral currents
 
   nuES - elastic scattering off target electrons 


# Run numbers

 - run numbers are taken from the official FASER conventions

    200025 flux from pion decay (i.e. basically numu)
 
    200026 flux from kaon decay (i.e. mainly numu and nue)
 
    200035 flux from charm decay (i.e. numu, nue and some nutau)

# Instructions for Reading ROOT Files using PyROOT

Please check the directory `Python_io`
