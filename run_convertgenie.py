#!/usr/bin/env python3
"""
run_convertgenie.py -- run ConvertGENIE's `ConvertGENIE.exe` once per
detector (3DCAL, ECAL, AHCAL by default), with named, auto-discovered
arguments instead of hand-typing a file list and remembering positional
argv order.

Why: ConvertGENIE.exe used to take exactly one <genieroot> pattern and a
positional <run> number (`ConvertGENIE.exe [options] <genieroot> <run>
[detector]`), and had to be invoked by hand once per detector to get
separate 3DCAL/ECAL/AHCAL output files -- with no record afterwards of
how many events of each interaction type (nueCC, numuCC, nutauCC, NC, ES)
ended up in each one. ConvertGENIE.cc was changed alongside this script:
<genieroot> is now one or more positional file/pattern arguments (each
added to the same TChain, so a production split across several raw GENIE
files converts to a single output), and <run> moved from a positional
argument to a required -r <run> option -- see
ConvertGENIE/ConvertGENIE.cc's own usage text.

This script does the tedious parts for you:
  - discovers which *.gfaser.root files actually exist under
    $FASERDATA/GENIE/ and passes all of them (override with
    --input-files to convert only specific files/patterns);
  - locates the GDML geometry file under $FASERDATA/GDML/ automatically
    (override with --geometry-file), the same directory faserps itself
    exports FASERCAL_V10.gdml into (see run_batchreco.py's
    _default_geometry_file() for the same convention);
  - runs ConvertGENIE.exe once per detector in --detectors (default:
    3DCAL ECAL AHCAL -- deliberately not MuonSpec/ALL, since those are
    niche/debugging selections), each invocation producing its own
    FASERMC-PO-Run<run>-...-_<detector>.root (ConvertGENIE.cc appends
    the _<detector> suffix itself whenever detector != ALL). ConvertGENIE.cc
    writes that file (and its own interaction_plots.root / *_distribution.C
    histogram dump) to its own cwd rather than redirecting via FASERDATA
    itself, so this script runs it with cwd set to --output-dir (default:
    $FASERDATA/CVGENIE/Run<run>) instead of the repo root;
  - captures each run's "nueCC = ... numuCC = ... nutauCC = ... NC = ..."
    / "ES = ..." counters (TPOEvent::dump_stats(), already wired into
    ConvertGENIE.cc's event loop and already scoped to just the events
    written for that detector -- see convert_FASERMC's per-event
    want_<detector> filter, applied before update_stats()) and writes
    them into one consolidated run-summary log, plus keeps each
    detector's full raw console output as its own log file -- both
    under --output-dir alongside the ROOT outputs;
  - asks for confirmation before reusing an --output-dir that already has
    files in it (a re-run with the same --run number, most likely) -- and
    if you confirm (or pass --force), deletes that directory completely
    first rather than writing over/alongside what's there, so no stale
    file from a previous run (an old detector subset, old-format plots,
    ...) can survive into this one;
  - after all detectors finish, overlays each detector's reference
    histograms (h_z_all, h_x_all, h_y_all, h_xy_all, h_xz_all, h_yz_all --
    each detector's own full set is kept in its own
    interaction_plots_<detector>.root, not just the single-detector PNGs)
    onto one canvas per quantity via ConvertGENIE/make_overlay_plots.C,
    producing z_distribution_overlay.png etc. in --output-dir so the
    detectors' distributions can be compared directly. Skipped (with a
    warning, not an error) if the `root` executable isn't on PATH, or
    disabled outright with --no-overlay;
  - looks up the integrated luminosity (in fb^-1) each input GENIE file was
    generated at, parsed out of the gevgen_faser -l/-o flags in whichever
    run*.sh script under the GENIE directory produced it (e.g.
    runFASERCAL_v10_Aki2024.sh's 'gevgen_faser -l 1000.0 ... -o
    fasercal.Aki2024.v10.light'), and records it in the run-summary log --
    luminosity sets the scale (the event rate normalization) of everything
    converted, so it needs to travel with the output, not just live in a
    script someone has to go find again.

Usage:
    python3 run_convertgenie.py --run 10000
    python3 run_convertgenie.py --run 10000 --detectors 3DCAL
    python3 run_convertgenie.py --run 10000 --detectors 3DCAL ECAL AHCAL MuonSpec
    python3 run_convertgenie.py --run 10000 --charmonly
    python3 run_convertgenie.py --run 10000 --input-files data/GENIE/fasercal.Aki2024.v10.charm.0.gfaser.root
    python3 run_convertgenie.py --run 10000 --geometry-file /path/to/other.gdml
    python3 run_convertgenie.py --run 10000 --dry-run
"""
import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
DEFAULT_BUILD_DIR = REPO_ROOT / "build"

VALID_DETECTORS = ("3DCAL", "ECAL", "AHCAL", "MuonSpec", "ALL")

# The three detectors ConvertGENIE.exe is run against by default -- see the
# module docstring. MuonSpec/ALL are still available via --detectors but
# are not run unless asked for.
DEFAULT_DETECTORS = ("3DCAL", "ECAL", "AHCAL")

# Matches TPOEvent::dump_stats()'s own two printed lines:
#   " nueCC = <n> numuCC = <n> nutauCC = <n> NC = <n>\n ES = <n>\n"
# re.DOTALL lets ".*?" cross that line break between NC and ES.
STATS_PATTERN = re.compile(
    r"nueCC\s*=\s*(\d+).*?numuCC\s*=\s*(\d+).*?nutauCC\s*=\s*(\d+).*?NC\s*=\s*(\d+).*?ES\s*=\s*(\d+)",
    re.DOTALL,
)
STATS_FIELDS = ("nueCC", "numuCC", "nutauCC", "NC", "ES")

# Matches ConvertGENIE.cc main()'s own "The output file is <name>" line --
# parsed out rather than recomputed here so the two can't silently drift
# apart if ConvertGENIE.cc's naming logic ever changes.
OUTPUT_FILE_PATTERN = re.compile(r"The output file is\s+(\S+\.root)")


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
    _no_geometry_hint = (
        "       GDML geometry isn't fetched or checked in -- FASERG4 exports it as a\n"
        "       side effect of running faserps. Generate one cheaply, with no input\n"
        "       sample needed, by running:\n"
        "           python3 run_faserps.py --muons --n-events 1\n"
        "       or pass --geometry-file to point at an existing GDML file explicitly.")
    if not gdml_dir.is_dir():
        sys.exit(f"error: GDML directory not found: {gdml_dir}\n{_no_geometry_hint}")
    candidates = sorted(gdml_dir.glob("*.gdml"))
    if len(candidates) == 0:
        sys.exit(f"error: no .gdml file found under {gdml_dir}\n{_no_geometry_hint}")
    if len(candidates) > 1:
        names = ", ".join(c.name for c in candidates)
        sys.exit(f"error: multiple .gdml files found under {gdml_dir} ({names})\n"
                  f"       Pass --geometry-file to pick one explicitly.")
    return candidates[0]


# Matches a gevgen_faser invocation's -l (luminosity, fb^-1) and -o (output
# file prefix) flags, in either order, on the same line -- e.g.
# "gevgen_faser -l 1000.0 -r 0 -g ... -o fasercal.Aki2024.v10.light".
GEVGEN_LUMI_PATTERN = re.compile(r"-l\s+([0-9.eE+-]+)")
GEVGEN_OUTPUT_PATTERN = re.compile(r"-o\s+(\S+)")


def discover_luminosities(search_dirs):
    """Scans every *.sh script in search_dirs for gevgen_faser invocations
    and returns {output_prefix: luminosity_fb}, e.g. scanning
    runFASERCAL_v10_Aki2024.sh's two gevgen_faser lines gives
    {"fasercal.Aki2024.v10.light": 1000.0, "fasercal.Aki2024.v10.charm": 1000.0}.
    Luminosity (fb^-1) is the integrated luminosity that sample was
    generated at -- it sets the scale of every event count downstream, so
    this is what lets the run-summary log say what statistics it
    actually represents instead of just a raw event count."""
    lumi_by_prefix = {}
    seen_dirs = set()
    for d in search_dirs:
        d = Path(d)
        if d in seen_dirs or not d.is_dir():
            continue
        seen_dirs.add(d)
        for script in sorted(d.glob("*.sh")):
            try:
                text = script.read_text()
            except OSError:
                continue
            for line in text.splitlines():
                if "gevgen_faser" not in line:
                    continue
                lum_match = GEVGEN_LUMI_PATTERN.search(line)
                out_match = GEVGEN_OUTPUT_PATTERN.search(line)
                if lum_match and out_match:
                    lumi_by_prefix[out_match.group(1)] = float(lum_match.group(1))
    return lumi_by_prefix


def lookup_luminosity(genie_file, lumi_by_prefix):
    """Matches a genie file's name against the discovered -o prefixes (the
    actual file is e.g. "fasercal.Aki2024.v10.light.0.gfaser.root", with
    ".0.gfaser.root" appended by GENIE's own run-numbering convention onto
    the "-o fasercal.Aki2024.v10.light" prefix), preferring the longest/most
    specific prefix match. Returns None if no script accounts for this file."""
    name = Path(genie_file).name
    best = None
    for prefix, lumi in lumi_by_prefix.items():
        if name.startswith(prefix) and (best is None or len(prefix) > len(best[0])):
            best = (prefix, lumi)
    return best[1] if best else None


def format_luminosity_lines(genie_files, luminosities):
    """Shared by print_plan() and format_summary() so the console plan and
    the persisted run-summary log agree on how this is reported."""
    lines = ["Luminosity (fb^-1) of input GENIE sample(s) -- sets the overall event-rate scale:"]
    values = []
    for f in genie_files:
        lumi = luminosities.get(str(f))
        if lumi is None:
            lines.append(f"  {Path(f).name:<45} UNKNOWN")
        else:
            lines.append(f"  {Path(f).name:<45} {lumi:g}")
            values.append(lumi)
    distinct = sorted(set(values))
    if len(values) == len(genie_files) and len(distinct) == 1:
        lines.append(f"  -> consistent: {distinct[0]:g} fb^-1")
    elif len(distinct) > 1:
        lines.append("  -> WARNING: input files were generated at different luminosities (see "
                      "above) -- the combined sample's statistics are NOT a single clean "
                      "luminosity.")
    else:
        lines.append("  -> WARNING: luminosity not found for one or more input files -- check "
                      "the run*.sh script(s) that generated them under the GENIE directory.")
    return lines


def build_command(binary_path, *, run, genie_files, geometry_file, detector,
                   charm_only, taucc_only, luminosity_fb):
    """Builds the ConvertGENIE.exe argv list for a single detector,
    matching its own usage text: [options] <genieroot> [<genieroot> ...]
    [detector], with -r/-g/-l/-charmonly/-tauCConly as options (see
    ConvertGENIE/ConvertGENIE.cc). luminosity_fb is omitted (ConvertGENIE.cc
    then records "L=unknown" in its histogram titles) when the input files'
    luminosities aren't known or aren't all the same -- see
    resolve_consistent_luminosity()."""
    command = [str(binary_path), "-r", str(run), "-g", str(geometry_file)]
    if luminosity_fb is not None:
        command += ["-l", str(luminosity_fb)]
    if charm_only:
        command.append("-charmonly")
    if taucc_only:
        command.append("-tauCConly")
    command += [str(f) for f in genie_files]
    if detector != "ALL":
        command.append(detector)
    return command


def resolve_consistent_luminosity(genie_files, luminosities):
    """Returns the single luminosity value (fb^-1) to pass ConvertGENIE.exe
    via -l, or None if it can't say one cleanly -- either because some
    input file's luminosity is unknown, or because the inputs were
    generated at different luminosities (format_luminosity_lines() already
    warns about both cases to the console/log; this just decides what, if
    anything, is safe to hand to ConvertGENIE.exe for its histogram
    titles)."""
    values = [luminosities.get(str(f)) for f in genie_files]
    if any(v is None for v in values):
        return None
    distinct = set(values)
    if len(distinct) != 1:
        return None
    return distinct.pop()


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run ConvertGENIE's ConvertGENIE.exe once per detector, with "
                     "auto-discovered input files and geometry instead of hand-typed "
                     "positional arguments.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--run", type=int, required=True,
                         help="Output run number, passed through as -r <run> (required -- "
                              "deliberately no default, same reasoning as run_batchreco.py's --run).")
    parser.add_argument("--detectors", nargs="+", choices=VALID_DETECTORS,
                         default=list(DEFAULT_DETECTORS),
                         help="Detector(s) to convert; ConvertGENIE.exe is run once per entry, "
                              f"each producing its own _<detector>-suffixed output file "
                              f"(default: {' '.join(DEFAULT_DETECTORS)}).")
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
    parser.add_argument("--charmonly", action="store_true", help="Pass -charmonly.")
    parser.add_argument("--taucconly", action="store_true", help="Pass -tauCConly.")
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD_DIR,
                         help=f"CMake build directory containing bin/ConvertGENIE.exe "
                              f"(default: {DEFAULT_BUILD_DIR}).")
    parser.add_argument("--output-dir", type=Path, default=None,
                         help="Directory for ConvertGENIE.exe's own output -- the "
                              "FASERMC-PO-Run*.root files, interaction_plots.root, and the "
                              "*_distribution.C histogram dumps (ConvertGENIE.exe writes all of "
                              "these to its own cwd, which this script sets to this directory) "
                              "-- as well as the per-detector raw-output logs and the "
                              "consolidated run-summary log this script writes itself "
                              "(default: $FASERDATA/CVGENIE/Run<run>).")
    parser.add_argument("--dry-run", action="store_true",
                         help="Print the command(s) that would run, but don't execute them.")
    parser.add_argument("--force", "-f", action="store_true",
                         help="Don't ask for confirmation when --output-dir already exists and "
                              "has files in it -- just proceed (e.g. for a non-interactive/"
                              "scripted invocation).")
    parser.add_argument("--no-overlay", action="store_true",
                         help="Skip generating the overlaid-detectors reference PNGs after all "
                              "detector runs finish.")
    return parser.parse_args()


def print_plan(args, *, genie_files, geometry_file, output_dir, luminosities):
    print("[run_convertgenie] ==================== run plan ====================")
    print(f"[run_convertgenie] run:                      {args.run}")
    print(f"[run_convertgenie] build dir:                {args.build_dir}")
    print(f"[run_convertgenie] output dir:                {output_dir}")
    print(f"[run_convertgenie] geometry file:             {geometry_file}")
    print(f"[run_convertgenie] detectors ({len(args.detectors)}):              {' '.join(args.detectors)}")
    print(f"[run_convertgenie] charmonly:                 {args.charmonly}")
    print(f"[run_convertgenie] taucconly:                 {args.taucconly}")
    print(f"[run_convertgenie] input file(s) ({len(genie_files)}):")
    for f in genie_files:
        print(f"[run_convertgenie]   {f}")
    for line in format_luminosity_lines(genie_files, luminosities):
        print(f"[run_convertgenie] {line}")
    print(f"[run_convertgenie] dry-run:                   {args.dry_run}")
    print("[run_convertgenie] ===================================================")


def run_one_detector(command, *, detector, cwd, log_path):
    """Runs ConvertGENIE.exe for a single detector (with cwd set to
    --output-dir, since ConvertGENIE.exe writes its ROOT output relative
    to its own cwd -- see the module docstring), streaming its stdout to
    the console live (a conversion can take a while) while also capturing
    it -- both to parse out the stats/output-file lines afterwards, and
    to save verbatim as this detector's own log file."""
    print(f"[run_convertgenie] ---- {detector}: {' '.join(command)}")
    lines = []
    with subprocess.Popen(command, cwd=cwd, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT, text=True, bufsize=1) as proc:
        for line in proc.stdout:
            print(f"[{detector}] {line}", end="")
            lines.append(line)
        returncode = proc.wait()

    output_text = "".join(lines)
    log_path.write_text(output_text)

    stats = None
    m = STATS_PATTERN.search(output_text)
    if m:
        stats = {field: int(value) for field, value in zip(STATS_FIELDS, m.groups())}

    output_file = None
    m = OUTPUT_FILE_PATTERN.search(output_text)
    if m:
        output_file = m.group(1)

    return {
        "detector": detector,
        "returncode": returncode,
        "stats": stats,
        "output_file": output_file,
        "log_path": log_path,
    }


OVERLAY_MACRO_PATH = REPO_ROOT / "ConvertGENIE" / "make_overlay_plots.C"


def interaction_plots_path(output_dir, detector):
    """The per-detector histogram file ConvertGENIE.cc itself writes --
    see its own plot_suffix logic (mirrored here): interaction_plots.root
    with a _<detector> suffix whenever detector != ALL."""
    suffix = "" if detector == "ALL" else f"_{detector}"
    return output_dir / f"interaction_plots{suffix}.root"


def make_overlays(results, *, output_dir, run):
    """Runs ConvertGENIE/make_overlay_plots.C over every detector that
    actually produced an interaction_plots_<detector>.root, overlaying
    their reference histograms onto one canvas per quantity. Returns
    without error (just a printed warning) if `root` isn't on PATH or no
    detector's output file is found -- this is a nice-to-have on top of
    the per-detector PNGs, not something that should fail the whole run."""
    root_bin = shutil.which("root")
    if root_bin is None:
        print("[run_convertgenie] skipping detector overlay: `root` executable not found on "
              "PATH (source setup.sh?)")
        return

    manifest_lines = []
    for r in results:
        if r["returncode"] != 0:
            continue
        plots_path = interaction_plots_path(output_dir, r["detector"])
        if not plots_path.is_file():
            print(f"[run_convertgenie] skipping {r['detector']} in overlay: "
                  f"{plots_path.name} not found")
            continue
        manifest_lines.append(f"{r['detector']} {plots_path}")

    if not manifest_lines:
        print("[run_convertgenie] skipping detector overlay: no detector produced an "
              "interaction_plots_<detector>.root to overlay")
        return

    manifest_path = output_dir / f"_overlay_manifest_Run{run}.txt"
    manifest_path.write_text("\n".join(manifest_lines) + "\n")

    macro_call = f'{OVERLAY_MACRO_PATH}("{manifest_path}", "{output_dir}")'
    print(f"[run_convertgenie] ---- overlay: {root_bin} -l -b -q '{macro_call}'")
    result = subprocess.run([root_bin, "-l", "-b", "-q", macro_call],
                             cwd=REPO_ROOT, capture_output=True, text=True)
    print(result.stdout, end="")
    if result.stderr:
        print(result.stderr, end="")
    if result.returncode != 0:
        print(f"[run_convertgenie] warning: overlay step exited with code {result.returncode}")


def format_summary(args, results, *, genie_files, luminosities):
    col_w = 10
    fields = STATS_FIELDS + ("Total",)
    header = "Detector".ljust(12) + "".join(f.rjust(col_w) for f in fields)
    rule = "-" * len(header)
    lines = []
    lines.append("==================== ConvertGENIE run summary ====================")
    lines.append(f"Run:       {args.run}")
    lines.append(f"Detectors: {' '.join(args.detectors)}")
    lines.extend(format_luminosity_lines(genie_files, luminosities))
    lines.append(rule)
    lines.append(header)
    lines.append(rule)
    totals = {field: 0 for field in STATS_FIELDS}
    for r in results:
        if r["stats"] is None:
            lines.append(f"{r['detector'].ljust(12)}  (no stats parsed -- returncode {r['returncode']}, "
                          f"see {r['log_path'].name})")
            continue
        row_total = sum(r["stats"].values())
        row = r["detector"].ljust(12) + "".join(str(r["stats"][f]).rjust(col_w) for f in STATS_FIELDS)
        row += str(row_total).rjust(col_w)
        lines.append(row)
        for field in STATS_FIELDS:
            totals[field] += r["stats"][field]
    lines.append(rule)
    grand_total = sum(totals.values())
    total_row = "TOTAL".ljust(12) + "".join(str(totals[f]).rjust(col_w) for f in STATS_FIELDS)
    total_row += str(grand_total).rjust(col_w)
    lines.append(total_row)
    lines.append(rule)
    lines.append("Output files:")
    for r in results:
        name = r["output_file"] or "(unknown -- not found in console output)"
        lines.append(f"  {r['detector']:<10} {name}")
    lines.append("====================================================================")
    return "\n".join(lines) + "\n"


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

    output_dir = args.output_dir if args.output_dir is not None else faserdata / "CVGENIE" / f"Run{args.run}"

    # Luminosity search dirs: each genie file's own parent directory (glob
    # patterns' ".parent" still resolves fine, e.g. Path("x/*.root").parent
    # == Path("x")), plus the default GENIE dir in case --input-files
    # pointed elsewhere but the generating run*.sh script is still there.
    lumi_search_dirs = {Path(f).resolve().parent if "*" not in str(f) and "?" not in str(f)
                         else Path(f).parent for f in genie_files}
    lumi_search_dirs.add(faserdata / "GENIE")
    lumi_by_prefix = discover_luminosities(lumi_search_dirs)
    luminosities = {str(f): lookup_luminosity(f, lumi_by_prefix) for f in genie_files}

    print_plan(args, genie_files=genie_files, geometry_file=geometry_file, output_dir=output_dir,
               luminosities=luminosities)

    resolved_luminosity_fb = resolve_consistent_luminosity(genie_files, luminosities)

    commands = {
        detector: build_command(
            binary_path, run=args.run, genie_files=genie_files, geometry_file=geometry_file,
            detector=detector, charm_only=args.charmonly, taucc_only=args.taucconly,
            luminosity_fb=resolved_luminosity_fb,
        )
        for detector in args.detectors
    }

    if args.dry_run:
        for detector in args.detectors:
            print(f"[run_convertgenie] ---- {detector}: {' '.join(commands[detector])}")
        print("[run_convertgenie] --dry-run: not executing.")
        return 0

    if output_dir.is_dir():
        existing = sorted(p.name for p in output_dir.iterdir())
        if existing:
            preview = ", ".join(existing[:5]) + (", ..." if len(existing) > 5 else "")
            print(f"[run_convertgenie] warning: --output-dir already exists and has "
                  f"{len(existing)} file(s) in it: {preview}")
            print(f"[run_convertgenie]   {output_dir}")
            if not args.force:
                reply = input("[run_convertgenie] Delete it and regenerate from scratch? [y/N] ")
                if reply.strip().lower() not in ("y", "yes"):
                    print("[run_convertgenie] aborted -- nothing was run.")
                    return 1
            # Confirmed (interactively, or via --force): remove the whole
            # directory rather than writing over/alongside what's there, so
            # no stale file from a previous run can survive into this one.
            print(f"[run_convertgenie] removing existing output dir before regenerating: {output_dir}")
            shutil.rmtree(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    results = []
    for detector in args.detectors:
        log_path = output_dir / f"ConvertGENIE_Run{args.run}_{detector}.log"
        result = run_one_detector(commands[detector], detector=detector, cwd=output_dir, log_path=log_path)
        results.append(result)
        if result["returncode"] != 0:
            print(f"[run_convertgenie] warning: {detector} run exited with code "
                  f"{result['returncode']} (see {log_path})")

    summary_text = format_summary(args, results, genie_files=genie_files, luminosities=luminosities)
    print()
    print(summary_text)

    summary_path = output_dir / f"ConvertGENIE_Run{args.run}_summary.log"
    summary_path.write_text(summary_text)
    print(f"[run_convertgenie] wrote run summary: {summary_path}")

    if not args.no_overlay:
        make_overlays(results, output_dir=output_dir, run=args.run)

    return 0 if all(r["returncode"] == 0 for r in results) else 1


if __name__ == "__main__":
    sys.exit(main())
