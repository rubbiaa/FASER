#!/usr/bin/env python3
"""
run_convertgenie.py -- run ConvertGENIE's `ConvertGENIE.exe` with named,
auto-discovered arguments instead of hand-typing a file list and
remembering positional argv order.

Why: ConvertGENIE.exe used to take exactly one <genieroot> pattern and a
positional <run> number (`ConvertGENIE.exe [options] <genieroot> <run>
[detector]`). That made it awkward to convert a production that is split
across several raw GENIE files (e.g. the Aki2024 v10 "charm" and "light"
channels fetched by fetch_data.py into $FASERDATA/GENIE/ as
fasercal.Aki2024.v10.charm.0.gfaser.root and
fasercal.Aki2024.v10.light.0.gfaser.root) without either running the
executable twice (two separate output files to merge by hand) or
constructing a shell glob yourself. ConvertGENIE.cc was changed alongside
this script: <genieroot> is now one or more positional arguments (each
added to the same TChain, so multiple files become one converted output),
and <run> moved from a positional argument to a required `-r <run>`
option -- see ConvertGENIE/ConvertGENIE.cc's own usage text.

This script does the two tedious parts for you:
  - discovers which *.gfaser.root files actually exist under
    $FASERDATA/GENIE/ and passes all of them (override with
    --input-files to convert only specific files/patterns);
  - locates the GDML geometry file under $FASERDATA/GDML/ automatically
    (override with --geometry-file), the same directory faserps itself
    exports FASERCAL_V10.gdml into (see run_batchreco.py's
    _default_geometry_file() for the same convention).

Usage:
    python3 run_convertgenie.py --run 10000
    python3 run_convertgenie.py --run 10000 --detector 3DCAL
    python3 run_convertgenie.py --run 10000 --charmonly
    python3 run_convertgenie.py --run 10000 --input-files data/GENIE/fasercal.Aki2024.v10.charm.0.gfaser.root
    python3 run_convertgenie.py --run 10000 --geometry-file /path/to/other.gdml
    python3 run_convertgenie.py --run 10000 --dry-run
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_BUILD_DIR = REPO_ROOT / "build"

VALID_DETECTORS = ("3DCAL", "ECAL", "AHCAL", "MuonSpec", "ALL")


# Mirrors run_faserps.py's _default_genie_input_file() / run_batchreco.py's
# _default_geometry_file(): FASERDATA falls back to REPO_ROOT/data so this
# resolves whether or not setup.sh has been sourced yet.
def _faserdata_dir() -> Path:
    return Path(os.environ.get("FASERDATA", str(REPO_ROOT / "data")))


def discover_genie_files(genie_dir: Path):
    """Returns every *.gfaser.root file under genie_dir, sorted for a
    deterministic, reproducible TChain order."""
    if not genie_dir.is_dir():
        return []
    return sorted(genie_dir.glob("*.gfaser.root"))


def discover_geometry_file(gdml_dir: Path) -> Path:
    """Finds the single *.gdml file under gdml_dir. Unlike
    run_batchreco.py (which hardcodes FASERCAL_V10.gdml, the one file that
    directory has always held), this looks at whatever is actually there --
    the user asked for "-g comes from the file in data/GDML/", not a
    hardcoded name -- so it only succeeds when exactly one candidate
    exists; with zero or several it errors out and tells the caller to
    pass --geometry-file explicitly."""
    if not gdml_dir.is_dir():
        sys.exit(f"error: GDML directory not found: {gdml_dir}\n"
                  f"       Pass --geometry-file to point at a GDML file explicitly.")
    candidates = sorted(gdml_dir.glob("*.gdml"))
    if len(candidates) == 0:
        sys.exit(f"error: no .gdml file found under {gdml_dir}\n"
                  f"       Pass --geometry-file to point at one explicitly.")
    if len(candidates) > 1:
        names = ", ".join(c.name for c in candidates)
        sys.exit(f"error: multiple .gdml files found under {gdml_dir} ({names})\n"
                  f"       Pass --geometry-file to pick one explicitly.")
    return candidates[0]


def build_command(binary_path, *, run, genie_files, geometry_file, detector,
                   charm_only, taucc_only):
    """Builds the ConvertGENIE.exe argv list, matching its own usage text:
    [options] <genieroot> [<genieroot> ...] [detector], with -r/-g/
    -charmonly/-tauCConly as options (see ConvertGENIE/ConvertGENIE.cc)."""
    command = [str(binary_path), "-r", str(run), "-g", str(geometry_file)]
    if charm_only:
        command.append("-charmonly")
    if taucc_only:
        command.append("-tauCConly")
    command += [str(f) for f in genie_files]
    if detector != "ALL":
        command.append(detector)
    return command


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run ConvertGENIE's ConvertGENIE.exe with auto-discovered input "
                     "files and geometry instead of hand-typed positional arguments.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--run", type=int, required=True,
                         help="Output run number, passed through as -r <run> (required -- "
                              "deliberately no default, same reasoning as run_batchreco.py's --run).")
    parser.add_argument("--input-files", nargs="+", default=None, metavar="FILE",
                         help="Explicit list of genieroot files/glob patterns to convert, in "
                              "place of auto-discovering every *.gfaser.root file under "
                              "$FASERDATA/GENIE/.")
    parser.add_argument("--genie-dir", type=Path, default=None,
                         help="Directory to auto-discover *.gfaser.root files in when "
                              "--input-files isn't given (default: $FASERDATA/GENIE).")
    parser.add_argument("--geometry-file", type=Path, default=None,
                         help="GDML geometry file, passed through as -g <gdmlfile> (default: "
                              "auto-discovered as the single *.gdml file under $FASERDATA/GDML).")
    parser.add_argument("--detector", choices=VALID_DETECTORS, default="ALL",
                         help="Detector selection (default: ALL).")
    parser.add_argument("--charmonly", action="store_true", help="Pass -charmonly.")
    parser.add_argument("--taucconly", action="store_true", help="Pass -tauCConly.")
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD_DIR,
                         help=f"CMake build directory containing bin/ConvertGENIE.exe "
                              f"(default: {DEFAULT_BUILD_DIR}).")
    parser.add_argument("--dry-run", action="store_true",
                         help="Print the command that would run, but don't execute it.")
    return parser.parse_args()


def print_run_summary(args, *, genie_files, geometry_file):
    print("[run_convertgenie] ==================== run summary ====================")
    print(f"[run_convertgenie] run:                      {args.run}")
    print(f"[run_convertgenie] build dir:                {args.build_dir}")
    print(f"[run_convertgenie] geometry file:             {geometry_file}")
    print(f"[run_convertgenie] detector:                  {args.detector}")
    print(f"[run_convertgenie] charmonly:                 {args.charmonly}")
    print(f"[run_convertgenie] taucconly:                 {args.taucconly}")
    print(f"[run_convertgenie] input file(s) ({len(genie_files)}):")
    for f in genie_files:
        print(f"[run_convertgenie]   {f}")
    print(f"[run_convertgenie] dry-run:                   {args.dry_run}")
    print("[run_convertgenie] ======================================================")


def main():
    args = parse_args()

    binary_path = args.build_dir / "bin" / "ConvertGENIE.exe"
    if not binary_path.is_file():
        sys.exit(
            f"error: ConvertGENIE.exe not found: {binary_path}\n"
            f"       (build it first, e.g. `cmake --build {args.build_dir} --target ConvertGENIE.exe`)"
        )
    if not (binary_path.stat().st_mode & 0o111):
        sys.exit(f"error: {binary_path} is not executable")

    faserdata = _faserdata_dir()

    if args.input_files:
        genie_files = [Path(f) for f in args.input_files]
        for f in genie_files:
            if not f.is_file() and "*" not in str(f) and "?" not in str(f):
                sys.exit(f"error: input file not found: {f}")
    else:
        genie_dir = args.genie_dir if args.genie_dir is not None else faserdata / "GENIE"
        genie_files = discover_genie_files(genie_dir)
        if not genie_files:
            sys.exit(
                f"error: no *.gfaser.root files found under {genie_dir}\n"
                f"       Fetch them first (see fetch_data.py), or pass --input-files/--genie-dir "
                f"explicitly."
            )

    if args.geometry_file is not None:
        geometry_file = args.geometry_file.resolve()
        if not geometry_file.is_file():
            sys.exit(f"error: geometry file not found: {geometry_file}")
    else:
        geometry_file = discover_geometry_file(faserdata / "GDML").resolve()

    print_run_summary(args, genie_files=genie_files, geometry_file=geometry_file)

    command = build_command(
        binary_path, run=args.run, genie_files=genie_files, geometry_file=geometry_file,
        detector=args.detector, charm_only=args.charmonly, taucc_only=args.taucconly,
    )
    print(f"[run_convertgenie] command: {' '.join(command)}")

    if args.dry_run:
        print("[run_convertgenie] --dry-run: not executing.")
        return 0

    # ConvertGENIE.exe writes its FASERMC-PO-Run*.root output to its own
    # cwd (see ConvertGENIE.cc's ROOTOutputFile, a bare relative filename --
    # no FaserDataDir.hh-style FASERDATA redirection like faserps/batchreco
    # use), so run it from the repo root rather than FASERG4/ or Batch/,
    # consistent with where this script itself lives.
    result = subprocess.run(command, cwd=REPO_ROOT)
    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
