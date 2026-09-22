# Running batchreco with `run_batchreco.py`

`run_batchreco.py` (at the repo root) runs `Batch/batchreco.exe` with named,
validated command-line arguments, instead of remembering `batchreco.exe`'s
positional argument order or relying on a hand-copied shell script.

## Prerequisites

- `batchreco.exe` must be built:

  ```
  cmake --build build --target batchreco.exe
  ```

- `FASERDATA` must be set - `source setup.sh` (or `common_setup.sh`)
  before running this script; it defaults `FASERDATA` to `$HOMEFASER/data`
  (see "Where the output goes" below and `common_setup.sh`). The script
  checks this itself and exits with a clear message if it's missing,
  rather than letting `batchreco.exe` fail with a confusing C++ exception.

- Run the script from anywhere - it resolves its own location via
  `__file__`. Unlike the old `Batch/go` (or running `batchreco.exe`
  directly from `Batch/`), this script does not depend on any particular
  working directory for input/output: `batchreco.exe` reads and writes
  under `$FASERDATA` directly (see "Where the output goes" below), so the
  script no longer needs to set `cwd=` for the subprocess at all.

## Basic usage

```
python3 run_batchreco.py --run 10000
```

Reconstructs every event of run 10000 using the default geometry
(`$FASERDATA/GDML/FASERCAL_V10.gdml`), writing
`data/batch/Batch-TPORecevent_10000_0_99999999.root`.

## Options

| Flag | Default | Meaning |
|---|---|---|
| `--run N` | *(required)* | Run number. No default - unlike `Batch/go`'s hardcoded `1`, silently reconstructing the wrong run is worse than being forced to say which one you mean. |
| `--min-event N` | `0` | First event index. |
| `--max-event N` | `99999999` | Event index to stop before (**exclusive**), matching `BatchReco.cc`'s own `ievent < max_event` loop. |
| `--mask {nueCC,numuCC,nutauCC,nuNC,nuES}` | none | Process only events of this interaction type. |
| `--multi-thread` | off | Passes `-mt` to `batchreco.exe`. |
| `--geometry-file PATH` | `$FASERDATA/GDML/FASERCAL_V10.gdml` | GDML geometry to load (`-g`). |
| `--build-dir PATH` | `./build` | CMake build directory containing `bin/batchreco.exe`. |
| `--split N` | `1` | Split `[min-event, max-event)` into `N` contiguous, gap-free chunks and launch them as background jobs, each with its own log file under `data/batch/logs/`. `1` means a single synchronous run whose output you see directly. |
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
  pid=12345  log=/path/to/data/batch/logs/batchreco_run10000_0_2000.log
  ...
[run_batchreco] not waiting for them -- monitor with `tail -f <log>` or `ps`.
```

## Where the output goes

`batchreco.exe` reads its input and writes its output ROOT file under
`$FASERDATA/faserG4/` and `$FASERDATA/batch/` respectively - via
`FASER::GetDataDir("faserG4")` / `FASER::GetDataDir("batch")`
(`CoreUtils/FaserDataDir.hh`), called from `Batch/BatchReco.cc` - instead
of a relative `input/` path and cwd-dependent output. `GetDataDir()`
creates these directories itself if they don't exist yet.

`FASERDATA` defaults to `$HOMEFASER/data` (set in `common_setup.sh`, only
if not already exported - see `setup.sh` for how a site or a specific
checkout can point it elsewhere, e.g. scratch/EOS space instead of inside
the git checkout). This is the same consolidated top-level `data/`
directory `run_faserps.py` writes into (`$FASERDATA/faserG4/` - see
`run_faserps.md`); every C++ executable that reads or writes
FASERG4/batchreco data - `BatchReco.cc`, `BatchReco_DetResp.cc`,
`DumpHits.cc`, `AnalyReco.cc`, `FileMask.cc`, `MyMainFrame.cc` - resolves
through the same `FASER::GetDataDir()` helper, so there's exactly one
place, `$FASERDATA`, that decides where the data actually lives. There
are no longer any `input`/`output` symlinks between subdirectories for
this (there used to be: `FASERG4/output`, `Batch/input`, `EvDisplay/input`,
`data/batch/input` - all removed), and no more cwd-dependent behaviour
either.

`data/` (or wherever `FASERDATA` points, if overridden) is gitignored; no
`.root` file or log under it is ever tracked.

## Location

Both this file and `run_batchreco.py` live at the repo root, alongside
`run_faserps.py`/`run_faserps.md`, for the same reason: a script and its
doc are easiest to keep in sync when they're next to each other, and both
scripts play the same project-wide "how do I run this" role rather than
being an implementation detail of one package.
