# Running batchreco with `run_batchreco.py`

`run_batchreco.py` (at the repo root) runs `Batch/batchreco.exe` with named,
validated command-line arguments, instead of remembering `batchreco.exe`'s
positional argument order or relying on a hand-copied shell script.

## Prerequisites

- `batchreco.exe` must be built:

  ```
  cmake --build build --target batchreco.exe
  ```

- Run the script from anywhere - it resolves its own location via
  `__file__`. The `batchreco.exe` subprocess itself always runs with
  `data/batch/` as its working directory (not `Batch/`, where the source
  and build files live), because it reads from a relative `input/` path
  and writes its output ROOT file into the current directory - see
  "Where the output goes" below. The script warns if `data/batch/input`
  is missing.

## Basic usage

```
python3 run_batchreco.py --run 10000
```

Reconstructs every event of run 10000 using the default geometry
(`FASERG4/FASERCAL_V10.gdml`), writing
`data/batch/Batch-TPORecevent_10000_0_99999999.root`.

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
  pid=12345  log=data/batch/logs/batchreco_run10000_0_2000.log
  ...
[run_batchreco] not waiting for them -- monitor with `tail -f <log>` or `ps`.
```

## Where the output goes

`batchreco.exe` itself always writes its output ROOT file straight into its current working directory (no `output/` subdirectory of its own, unlike `faserps`/`FASERG4`). This script runs it with `data/batch/` as that working directory instead of `Batch/`, so the output lands in the same consolidated top-level `data/` tree as `run_faserps.py`'s output (`data/faserG4/` - see `run_faserps.md`).

`data/batch/input` is a symlink to `../faserG4` - the same underlying files `Batch/input` already points at via `../FASERG4/output` - so `batchreco.exe`'s own relative `input/` lookup keeps resolving to the right files with no change to `BatchReco.cc`.

`Batch/input` itself is untouched and still works if you run `batchreco.exe` directly from `Batch/` (bypassing this script) - but its output would then land back in `Batch/` rather than `data/batch/`, so prefer this script for anything you want consolidated. One known rough edge: `EvDisplay/src/MyMainFrame.cc` has a hardcoded `../Batch/Batch-TPORecevent_*` lookup that assumes output still lands directly in `Batch/` - it won't see new output produced via this script until that one hardcoded path is updated too (not done as part of this change, to keep it scoped to the two wrapper scripts).

`data/` is gitignored in its entirety; no `.root` file, log, or symlink under it is ever tracked.

## Location

Both this file and `run_batchreco.py` live at the repo root, alongside
`run_faserps.py`/`run_faserps.md`, for the same reason: a script and its
doc are easiest to keep in sync when they're next to each other, and both
scripts play the same project-wide "how do I run this" role rather than
being an implementation detail of one package.
