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

- Run the script from the repo root (`python3 run_faserps.py`), or with a
  full/relative path to it from anywhere - it resolves its own location via
  `__file__` and finds `FASERG4/` under it either way. Regardless of where
  you invoke it *from*, the `faserps` subprocess itself always runs with
  `FASERG4/` as its working directory, because the macro's
  `/generator/rootinputfilename` and the `output/` directory it writes into
  are relative paths that only resolve correctly from there - same as the
  GDML file `faserps` writes out.

- The script creates `FASERG4/output/` if it doesn't already exist. ROOT's
  `TFile` does not create missing directories on its own, so without this a
  run fails partway through with a confusing ROOT error instead of a clear
  one up front. On a normal checkout `FASERG4/output` is a symlink to
  `../data/faserG4` (see "Where the output goes" below), so this is a no-op
  in practice - it only creates a real directory on a fresh checkout that
  hasn't set that symlink up yet.

## Where the output goes

`faserps` itself always writes to the relative path `output/FASERG4-Tcalevent_<run>_<event>.root` (hardcoded in `CoreUtils/TcalEvent.cc`), which resolves to `FASERG4/output/` because that's the subprocess's working directory (see above).

`FASERG4/output` is a symlink to `../data/faserG4` - the repo keeps all simulation/reconstruction output consolidated under a top-level `data/` directory (`data/faserG4/` for this script, `data/batch/` for `run_batchreco.py` - see `run_batchreco.md`) instead of scattered inside `FASERG4/` and `Batch/` themselves. Existing consumers of the old path (`Batch/input`, `EvDisplay/input`, both symlinked to `../FASERG4/output`) keep working unchanged, since they resolve through this symlink transparently - nothing in the C++ needed to change.

Both `output/` and `data/` are gitignored; no `.root` file or these symlinks are ever tracked.

## Basic usage

```
python3 run_faserps.py
```

This runs `faserps` in batch mode with the same settings as the old
`runFASER_V10.mac`: the tilted V10 geometry, 100 neutrino-interaction events
read starting from event 0 of `FASERMC-PO-Run10000-0_53954_3DCAL.root`.

## Options

| Flag | Default | Meaning |
|---|---|---|
| `--build-dir PATH` | `../build` | CMake build directory containing `bin/faserps`. |
| `--input-file NAME` | `FASERMC-PO-Run10000-0_53954_3DCAL.root` | Value for `/generator/rootinputfilename` (relative to `FASERG4/`). Ignored if `--muons` is given. |
| `--start-event N` | `0` | Value for `/generator/startevent`. Ignored if `--muons` is given. |
| `--n-events N` | `100` | Value for `/run/beamOn`. |
| `--muons` | off | Generate single fixed-momentum muons instead of reading neutrino-interaction events (see below). |
| `--muon-momentum-gev X` | `100` | Muon momentum in GeV. Only used with `--muons`. |
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
