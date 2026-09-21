#!/usr/bin/env python3
"""
run_batchreco.py -- run Batch's `batchreco.exe` reconstruction with named,
validated arguments instead of remembering positional argv order and
hand-copied shell scripts.

Why: Batch/go hardcoded six parallel invocations of batchreco.exe (run 1,
event chunks "0 1999", "2000 3999", ..., "10000 10999") with no
explanation of what the numbers meant, no output separation (all six
processes wrote straight to the same terminal), and - because
BatchReco.cc's event loop is `ievent < max_event` (max_event is
EXCLUSIVE) - a real off-by-one gap: event 1999 fell between the "0 1999"
and "2000 3999" chunks and was silently never processed, and likewise
3999, 5999, 7999, 9999. This script replaces it: --split N divides
[min-event, max-event) into N contiguous, gap-free chunks and launches
them as background jobs with their own log files; a single run (the
default, no --split) runs synchronously so you see its output directly.

It also does not reuse BatchReco.cc's own hardcoded default geometry file
("../GeomGDML/geometry.gdml", relative to Batch/) - that file does not
actually exist in this checkout (GeomGDML/ has geometry_v5.gdml,
geometry_15_noshiftLOS.gdml, geometry_tilted_5degree.gdml and
FaserNu3.gdml, but no plain geometry.gdml), so relying on it is itself a
silent-wrong-geometry footgun (the same class of bug this whole
run_faserps.py/run_batchreco.py effort is about). This script's own
default instead points at FASERG4/FASERCAL_V10.gdml, the geometry FASERG4
currently actually exports - override with --geometry-file if you need a
different one.

Usage:
    python3 run_batchreco.py --run 10000
    python3 run_batchreco.py --run 10000 --min-event 0 --max-event 2000
    python3 run_batchreco.py --run 10000 --mask numuCC
    python3 run_batchreco.py --run 10000 --max-event 12000 --split 6
    python3 run_batchreco.py --run 10000 --multi-thread
    python3 run_batchreco.py --run 10000 --geometry-file /path/to/other.gdml
    python3 run_batchreco.py --run 10000 --dry-run
    python3 run_batchreco.py --run 10000 --max-event 12000 --split 6 --dry-run
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_BUILD_DIR = REPO_ROOT / "build"

# BatchReco.cc's own hardcoded default ("../GeomGDML/geometry.gdml",
# relative to Batch/) does not exist in this checkout -- see the module
# docstring. Always pass an explicit -g instead of relying on it.
DEFAULT_GEOMETRY_FILE = REPO_ROOT / "FASERG4" / "FASERCAL_V10.gdml"

VALID_MASKS = ("nueCC", "numuCC", "nutauCC", "nuNC", "nuES")


def build_command(binary_path, *, run, min_event, max_event, mask, multi_thread, geometry_file):
    """Builds the batchreco.exe argv list. Order matches BatchReco.cc's own
    usage: [-mt] [-g <geometryFile>] <run> [minevent] [maxevent] [mask]."""
    command = [str(binary_path)]
    if multi_thread:
        command.append("-mt")
    command += ["-g", str(geometry_file)]
    command += [str(run), str(min_event), str(max_event)]
    if mask:
        command.append(mask)
    return command


def split_ranges(min_event, max_event, n_splits):
    """Splits [min_event, max_event) into n_splits contiguous, gap-free
    chunks -- each chunk's max is exactly the next chunk's min, unlike the
    old Batch/go script (see module docstring)."""
    total = max_event - min_event
    if total <= 0:
        raise ValueError(f"max-event ({max_event}) must be greater than min-event ({min_event})")
    base, extra = divmod(total, n_splits)
    ranges = []
    start = min_event
    for i in range(n_splits):
        size = base + (1 if i < extra else 0)
        end = start + size
        ranges.append((start, end))
        start = end
    return ranges


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run Batch's batchreco.exe with named, validated arguments "
                     "instead of Batch/go's hardcoded parallel invocations.",
    )
    parser.add_argument("--run", type=int, required=True,
                         help="Run number (required -- deliberately no default, unlike Batch/go's "
                              "hardcoded 1, since silently reconstructing the wrong run is worse "
                              "than being forced to say which one you mean).")
    parser.add_argument("--min-event", type=int, default=0, help="First event index (default: 0).")
    parser.add_argument("--max-event", type=int, default=99999999,
                         help="Event index to stop before (EXCLUSIVE), matching BatchReco.cc's own "
                              "'ievent < max_event' loop (default: 99999999, i.e. no limit).")
    parser.add_argument("--mask", choices=VALID_MASKS, default=None,
                         help="Process only events of this interaction type (default: none/all).")
    parser.add_argument("--multi-thread", action="store_true", help="Pass -mt (enable multi-threading).")
    parser.add_argument("--geometry-file", type=Path, default=DEFAULT_GEOMETRY_FILE,
                         help=f"GDML geometry to load (default: {DEFAULT_GEOMETRY_FILE}).")
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD_DIR,
                         help=f"CMake build directory containing bin/batchreco.exe (default: {DEFAULT_BUILD_DIR}).")
    parser.add_argument("--split", type=int, default=1,
                         help="Split [min-event, max-event) into this many contiguous chunks and run "
                              "them as separate background jobs, each with its own log file under "
                              "data/batch/logs/ (default: 1, i.e. a single synchronous run).")
    parser.add_argument("--dry-run", action="store_true",
                         help="Print the command(s) that would run, but don't execute them.")
    return parser.parse_args()


def main():
    args = parse_args()

    binary_path = args.build_dir / "bin" / "batchreco.exe"
    if not binary_path.is_file():
        sys.exit(
            f"error: batchreco.exe not found: {binary_path}\n"
            f"       (build it first, e.g. `cmake --build {args.build_dir} --target batchreco.exe`)"
        )
    if not (binary_path.stat().st_mode & 0o111):
        sys.exit(f"error: {binary_path} is not executable")

    geometry_file = args.geometry_file.resolve()
    if not geometry_file.is_file():
        sys.exit(f"error: geometry file not found: {geometry_file}")

    # batchreco.exe reads and writes under $FASERDATA (see
    # CoreUtils/FaserDataDir.hh / Batch/BatchReco.cc), not a relative
    # "input/" path or its own cwd, so FASERDATA must be set in the
    # environment this subprocess inherits (source setup.sh first).
    faserdata = os.environ.get("FASERDATA")
    if not faserdata:
        sys.exit(
            "error: FASERDATA is not set.\n"
            "       Source setup.sh first (`source setup.sh`), or export FASERDATA\n"
            "       yourself to point at FASER's consolidated data directory."
        )
    log_dir = Path(faserdata) / "batch" / "logs"

    if args.split < 1:
        sys.exit("error: --split must be >= 1")

    ranges = split_ranges(args.min_event, args.max_event, args.split)

    print(f"[run_batchreco] FASERDATA: {faserdata}")

    if args.split == 1:
        command = build_command(
            binary_path, run=args.run, min_event=ranges[0][0], max_event=ranges[0][1],
            mask=args.mask, multi_thread=args.multi_thread, geometry_file=geometry_file,
        )
        print(f"[run_batchreco] command:           {' '.join(command)}")
        if args.dry_run:
            print("[run_batchreco] --dry-run: not executing.")
            return 0
        result = subprocess.run(command)
        return result.returncode

    # --split > 1: launch each chunk as its own background job, logging to
    # its own file -- unlike Batch/go, which sent every job's output
    # straight to the same terminal with nothing to tell them apart.
    if not args.dry_run:
        log_dir.mkdir(parents=True, exist_ok=True)

    jobs = []
    for i, (chunk_min, chunk_max) in enumerate(ranges):
        command = build_command(
            binary_path, run=args.run, min_event=chunk_min, max_event=chunk_max,
            mask=args.mask, multi_thread=args.multi_thread, geometry_file=geometry_file,
        )
        log_path = log_dir / f"batchreco_run{args.run}_{chunk_min}_{chunk_max}.log"
        print(f"[run_batchreco] chunk {i}: {' '.join(command)}  (log: {log_path})")
        if not args.dry_run:
            log_file = open(log_path, "w")
            proc = subprocess.Popen(command, stdout=log_file, stderr=subprocess.STDOUT)
            jobs.append((proc, log_path))

    if args.dry_run:
        print("[run_batchreco] --dry-run: not executing.")
        return 0

    print(f"[run_batchreco] launched {len(jobs)} background job(s):")
    for proc, log_path in jobs:
        print(f"  pid={proc.pid}  log={log_path}")
    print("[run_batchreco] not waiting for them -- monitor with `tail -f <log>` or `ps`.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
