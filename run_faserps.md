# Running faserps with `run_faserps.py`

`run_faserps.py` (at the repo root) runs FASERG4's `faserps` simulation
with the V10 geometry macro, built as a Python string and piped straight
to the executable's stdin. No `.mac` file is ever written to, or read
from, disk.

## Why this exists

`FASERG4/` used to carry about 18 `run*.mac` / `RunFASER_*.mac` files
(`runFASER_1.mac`, `runFASER_muon.mac`, `runFASER_V5.mac`, `runFASER_V10.mac`,
...), almost all of them leftovers from earlier detector geometries. Nothing
stopped `faserps` from being run with a stale one, which silently simulates
the wrong detector with no error to catch it.

All of those files have been removed. The one macro that matches the current
V10 geometry now lives in this script as `build_v10_macro()`, a parametrized
Python function, and its text is piped to `faserps -` at run time. A future
run variant (different momentum, different input file, a tweak to the tilt
angle, ...) is a new keyword argument or a new function, not a new file to
keep in sync by hand.

## Prerequisites

- `faserps` must be built, and built from a `faserps.cc` that has the `-`
  (stdin) mode. If you see `ERROR: Can not open a macro file <->` when
  running this script, the built binary predates that change and needs
  rebuilding:

  ```
  cmake --build build --target faserps
  ```

- `FASERDATA` must be set - `source setup.sh` (or `common_setup.sh`) before
  running this script; it defaults `FASERDATA` to `$HOMEFASER/data` (see
  "Where the output goes" below and `common_setup.sh`). The script checks
  this itself and exits with a clear message if it's missing, rather than
  letting `faserps` fail partway through with a confusing ROOT/C++ error.

- Run the script from the repo root (`python3 run_faserps.py`), or with a
  full/relative path to it from anywhere - it resolves its own location via
  `__file__` and finds `FASERG4/` under it either way. The `faserps`
  subprocess itself still always runs with `FASERG4/` as its working
  directory, but this no longer matters for the default `--input-file`
  (an absolute `$FASERDATA/GENIE/...` path - see "Where the input sample
  comes from" below) or the GDML geometry `faserps` writes out (also an
  absolute `$FASERDATA/GDML/...` path now, via `FASER::GetDataDir("GDML")`
  - see `run_batchreco.md`). It only still matters for a custom
  `--input-file` value that isn't already absolute, which is resolved
  relative to `FASERG4/` the same way it always was.

## Where the input sample comes from

The default `--input-file`, `FASERMC-PO-Run10000-0_53954_3DCAL.root`, is
too large to commit to git and lives under `$FASERDATA/GENIE/` (gitignored,
same as the rest of `$FASERDATA`). If it's missing there, `run_faserps.py`
fetches it automatically - before invoking `faserps` - from a public
CERNBox link (`fetch_data.py`'s `REMOTE_FILES` manifest), verifying its
sha256 after download. This only happens when `--input-file` is left at
its default; a custom `--input-file` is your own file, and is never
auto-fetched. You can also run `python3 fetch_data.py` directly (see its
own docstring) to fetch everything in the manifest up front, e.g. before
going offline, or `--force` to re-fetch regardless of what's already
there.

## Where the output goes

`faserps` writes its TcalEvent ROOT output under `$FASERDATA/faserG4/` -
via `FASER::GetDataDir("faserG4")` (`CoreUtils/FaserDataDir.hh`), called
from `CoreUtils/TcalEvent.cc` - instead of a hardcoded relative `output/`
path. `GetDataDir()` creates `$FASERDATA/faserG4` itself if it doesn't
exist yet, so nothing needs to `mkdir` it up front.

`FASERDATA` itself defaults to `$HOMEFASER/data` (set in
`common_setup.sh`, only if not already exported - see `setup.sh` for how a
site or a specific checkout can point it somewhere else, e.g. scratch/EOS
space instead of inside the git checkout). This is the same consolidated
top-level `data/` directory `run_batchreco.py` writes into
(`$FASERDATA/batch/` - see `run_batchreco.md`); the two scripts, and every
C++ executable that reads or writes FASERG4/batchreco data, all resolve
through the same `FASER::GetDataDir()` helper, so there's exactly one
place - `$FASERDATA` - that decides where the data actually lives. There
are no longer any `input`/`output` symlinks between subdirectories for
this (there used to be: `FASERG4/output`, `Batch/input`, `EvDisplay/input`,
`data/batch/input` - all removed).

`data/` (or wherever `FASERDATA` points, if overridden) is gitignored; no
`.root` file under it is ever tracked.

## Basic usage

```
python3 run_faserps.py
```

This runs `faserps` in batch mode with the same settings as the old
`runFASER_V10.mac`: the tilted V10 geometry, 100 neutrino-interaction events
read starting from event 0 of `$FASERDATA/GENIE/FASERMC-PO-Run10000-0_53954_3DCAL.root`.

## Options

| Flag | Default | Meaning |
|---|---|---|
| `--build-dir PATH` | `../build` | CMake build directory containing `bin/faserps`. |
| `--input-file NAME` | `$FASERDATA/GENIE/FASERMC-PO-Run10000-0_53954_3DCAL.root` | Value for `/generator/rootinputfilename`. The default is an absolute path (the sample lives under `$FASERDATA/GENIE/`, not `FASERG4/`); a custom value that isn't already absolute is still resolved relative to `FASERG4/`, since that's faserps' cwd. Ignored if `--muons` or `--muondis` is given. |
| `--start-event N` | `0` | Value for `/generator/startevent`. Ignored if `--muons` or `--muondis` is given. |
| `--n-events N` | `100` | Value for `/run/beamOn`. |
| `--muons` | off | Generate single fixed-momentum muons instead of reading neutrino-interaction events (see below). |
| `--muon-momentum-gev X` | `100` | Muon momentum in GeV. Only used with `--muons`/`--muondis`. |
| `--muondis` | off | Enable MuonDIS (see below). Implies `--muons`. |
| `--muondis-cross-section-bias X` | `150` | Value for `/physics/muondis/crossSectionBias`. Only used with `--muondis`; use `1` for an unbiased cross section. |
| `--muondis-q2min X` | `1.0` | Value for `/physics/muondis/q2min` (GeV²). Only used with `--muondis`. |
| `--muondis-interaction-log PATH` | off | Value for `/physics/muondis/interactionLog`. Only used with `--muondis`. |
| `--muondis-pdf-set PATH` | off | Value for `/physics/muondis/pdfSet`. Only used with `--muondis`. |
| `--muondis-xbjmin X` | `0` | Value for `/physics/muondis/xbjmin`. Only used with `--muondis`. |
| `--muondis-debug` | off | Add `/physics/muondis/debug true`. Only used with `--muondis`. |
| `--tilt-deg X` | `-4.5` | Value for `/FASER/tiltY`. |
| `--shift-x-cm X` | `45` | Value for `/FASER/LOS/shiftX`. |
| `--shift-y-cm X` | `24` | Value for `/FASER/LOS/shiftY`. |
| `--print-macro` | off | Print the resolved macro text and exit. Nothing is run. |
| `--dry-run` | off | Print the command and the macro, but don't execute `faserps`. |
| `--vis` | off | Launch `faserps`' interactive Geant4 UI (`faserps vis`) instead of a batch run. No macro is used in this mode. |

The other geometry constants baked into the old `runFASER_V10.mac`
(scintillator size, target size, voxel size, number of layers, ...) are also
keyword arguments of `build_v10_macro()`, just not exposed as command-line
flags since they're rarely changed. Call the function yourself from a small
script if you need to override one of them.

## Muon mode

`--muons` adds these two lines to the macro instead of the
`/generator/rootinputfilename` / `/generator/startevent` pair:

```
/generator/wantMuonBackground true
/generator/singleMomentum 100 GeV
```

`PrimaryGeneratorAction::GeneratePrimaries` (see
`FASERG4/src/PrimaryGeneratorAction.cc`) only opens the neutrino-interaction
ROOT input file when both `wantMuonBackground` and `wantSingleParticle` are
false, so `--input-file`/`--start-event` are simply unused once `--muons` is
given - the script leaves those two lines out of the generated macro
entirely rather than printing them alongside the muon lines misleadingly.

## MuonDIS mode

`--muondis` additionally enables **MuonDIS**: the primary muon's nuclear
interaction is replaced by an on-the-fly Pythia8 deep-inelastic-scattering
event, generated per interaction from the muon's actual Geant4 energy and
direction, instead of Geant4's standard muon-nuclear final state. It implies
`--muons` (MuonDIS only applies to muon-background primaries) and adds these
lines to the macro, before `/run/initialize` (MuonDIS's own UI commands only
take effect if set before then):

```
/physics/muondis/enable true
/physics/muondis/crossSectionBias 150
/physics/muondis/q2min 1
```

`crossSectionBias`/`q2min` default to the same values
`FASERG4/include/MuonDISPhysics.hh` itself defaults to (150, biasing the
interaction cross section up so DIS events are frequent enough to study
without huge statistics; and 1 GeV², respectively), so the printed macro
stays self-documenting even when you don't override them -- pass
`--muondis-cross-section-bias 1` for an unbiased cross section.
`--muondis-interaction-log`/`--muondis-pdf-set`/`--muondis-xbjmin`/`--muondis-debug`
are all off unless given explicitly. See `docs/README_MuonDIS.md` for the
full physics (target treatment, PDF choice, truth-level output added to
`TPOEvent`, ...) and the complete list of `/physics/muondis/...` options.

## Examples

```
# Default: 100 neutrino-interaction events, V10 geometry
python3 run_faserps.py

# 500 events starting from event 100 of the same input file
python3 run_faserps.py --n-events 500 --start-event 100

# A different input file
python3 run_faserps.py --input-file some_other_sample.root

# 1000 single 100 GeV muons instead of neutrino events
python3 run_faserps.py --muons --n-events 1000

# Single muons at 250 GeV
python3 run_faserps.py --muons --muon-momentum-gev 250

# MuonDIS: 1000 muon-background events with the DIS interaction enabled
python3 run_faserps.py --muondis --n-events 1000

# MuonDIS with an unbiased cross section and a non-default nucleon PDF
python3 run_faserps.py --muondis --muondis-cross-section-bias 1 \
    --muondis-pdf-set input/NNPDF40_nnlo_as_01180_charmasy_0000.dat

# See the macro without running anything
python3 run_faserps.py --muons --print-macro

# See exactly what would be run, without running it
python3 run_faserps.py --muons --dry-run

# Interactive Geant4 UI
python3 run_faserps.py --vis

# Point at a non-default build directory
python3 run_faserps.py --build-dir /path/to/other/build
```

## Under the hood

`faserps.cc` gained a third invocation mode alongside its existing
`faserps <macro-file>` and `faserps vis`:

```
faserps -
```

which reads G4 UI commands from stdin, one per line, skipping blank lines
and `#`-prefixed comments exactly like a real macro file would, instead of
`/control/execute`-ing a file on disk. `run_faserps.py` builds the full
macro text in memory with `build_v10_macro(...)` and pipes it to this mode
via `subprocess.run(..., input=macro_text.encode())`. The old
`faserps <macro-file>` mode still works unchanged, in case anything else
still needs to run a real macro file.

## Location

Both this file and `run_faserps.py` live at the repo root, not inside
`FASERG4/` - they used to, but the script is the standard way to run a
FASERG4 simulation project-wide (the same role `setup.sh` plays for
configuring the build), and keeping a script and its doc apart invites them
to drift out of sync.
