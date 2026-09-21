# Running batchreco with `run_batchreco.py`

`run_batchreco.py` (at the repo root) runs `Batch/batchreco.exe` with named,
validated command-line arguments, instead of remembering `batchreco.exe`'s
positional argument order or relying on a hand-copied shell script.

## Why this exists

`Batch/go` used to be the way to launch a large reconstruction job in
parallel:

```
./batchreco.exe 1 0 1999 &
./batchreco.exe 1 2000 3999 &
./batchreco.exe 1 4000 5999 &
./batchreco.exe 1 6000 7999 &
./batchreco.exe 1 8000 9999 &
./batchreco.exe 1 10000 10999 &
```

Three problems with this, in increasing order of severity:

1. No explanation of what "1", "0", "1999", etc. mean - you had to go read
   `BatchReco.cc`'s usage text to know these are `<run> [minevent] [maxevent]`.
2. All six processes wrote their output straight to the same terminal, with
   nothing to tell one job's log lines from another's.
3. A real off-by-one gap: `BatchReco.cc`'s event loop is `ievent < max_event`
   (max_event is **exclusive**), so "0 1999" processes events 0-1998 and
   "2000 3999" processes 2000-3998. Event **1999 was never processed by any
   chunk** - and likewise 3999, 5999, 7999, 9999. Five events silently
   dropped, every time this script was used with these boundaries.

`run_batchreco.py` replaces it. `--split N` divides `[min-event, max-event)`
into `N` chunks whose boundaries are computed so each chunk's end is exactly
the next chunk's start - no gaps, by construction - and gives each chunk its
own log file instead of a shared terminal. `Batch/go` has been removed.

It also does not reuse `BatchReco.cc`'s own hardcoded default geometry file.

## The geometry file default is a separate bug this script works around

`BatchReco.cc` (and `BatchReco_DetResp.cc`, and `DumpHits.cc`) default to
loading `../GeomGDML/geometry.gdml` if you don't pass `-g` yourself. That
file **does not exist** in this checkout - `GeomGDML/` has `geometry_v5.gdml`,
`geometry_15_noshiftLOS.gdml`, `geometry_tilted_5degree.gdml` and
`FaserNu3.gdml`, but no plain `geometry.gdml`. Relying on that default is
itself a silent-wrong-geometry footgun, the same class of problem this whole
`run_faserps.py`/`run_batchreco.py` effort exists to close.

This script does not touch that C++ default (fixing it is out of scope
here - it's used by three different executables) but it never relies on it
either: it always passes an explicit `-g`, defaulting to
`FASERG4/FASERCAL_V10.gdml` - the geometry FASERG4 currently actually
exports. Override with `--geometry-file` if you need a different one.

## Prerequisites

- `batchreco.exe` must be built:

  ```
  cmake --build build --target batchreco.exe
  ```

- Run the script from anywhere - it resolves its own location via
  `__file__`. The `batchreco.exe` subprocess itself always runs with
  `Batch/` as its working directory, because it reads from a relative
  `input/` path (normally a symlink to `../FASERG4/output`) and writes its
  output ROOT file into the current directory. The script warns if
  `Batch/input` is missing.

## Basic usage

```
python3 run_batchreco.py --run 10000
```

Reconstructs every event of run 10000 using the default geometry
(`FASERG4/FASERCAL_V10.gdml`), writing
`Batch/Batch-TPORecevent_10000_0_99999999.root`.

## Options

| Flag | Default | Meaning |
|---|---|---|
| `--run N` | *(required)* | Run number. No default - unlike `Batch/go`'s hardcoded `1`, silently reconstructing the wrong run is worse than being forced to say which one you mean. |
| `--min-event N` | `0` | First event index. |
| `--max-event N` | `99999999` | Event index to stop before (**exclusive**), matching `BatchReco.cc`'s own `ievent < max_event` loop. |
| `--mask {nueCC,numuCC,nutauCC,nuNC,nuES}` | none | Process only events of this interaction type. |
| `--multi-thread` | off | Passes `-mt` to `batchreco.exe`. |
| `--geometry-file PATH` | `FASERG4/FASERCAL_V10.gdml` | GDML geometry to load (`-g`). |
| `--build-dir PATH` | `./build` | CMake build directory containing `bin/batchreco.exe`. |
| `--split N` | `1` | Split `[min-event, max-event)` into `N` contiguous, gap-free chunks and launch them as background jobs, each with its own log file under `Batch/logs/`. `1` means a single synchronous run whose output you see directly. |
| `--dry-run` | off | Print the command(s) that would run, but don't execute them. |

`batchreco.exe`'s own `-r` flag ("open reconstructed files"), listed in its
`--help`-style usage text, is not wired to anything in `BatchReco.cc`'s
actual argument parsing - it's dead text, not a real option - so it has no
equivalent flag here.

## Examples

```
# Reconstruct all of run 10000 with the default geometry
python3 run_batchreco.py --run 10000

# Just the first 2000 events
python3 run_batchreco.py --run 10000 --min-event 0 --max-event 2000

# Only numuCC events
python3 run_batchreco.py --run 10000 --mask numuCC

# Split events 0..12000 of run 10000 into 6 parallel background jobs
# (replaces Batch/go, without the output-file gaps)
python3 run_batchreco.py --run 10000 --max-event 12000 --split 6

# Multi-threaded, non-default geometry
python3 run_batchreco.py --run 10000 --multi-thread --geometry-file /path/to/other.gdml

# See the command(s) without running anything
python3 run_batchreco.py --run 10000 --max-event 12000 --split 6 --dry-run
```

When `--split` is used, the script launches every chunk in the background
and returns immediately (it does not wait for them), printing each job's
PID and log path:

```
[run_batchreco] launched 6 background job(s):
  pid=12345  log=Batch/logs/batchreco_run10000_0_2000.log
  ...
[run_batchreco] not waiting for them -- monitor with `tail -f <log>` or `ps`.
```

## Location

Both this file and `run_batchreco.py` live at the repo root, alongside
`run_faserps.py`/`run_faserps.md`, for the same reason: a script and its
doc are easiest to keep in sync when they're next to each other, and both
scripts play the same project-wide "how do I run this" role rather than
being an implementation detail of one package.
