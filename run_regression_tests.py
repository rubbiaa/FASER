#!/usr/bin/env python3
"""
run_regression_tests.py -- simulate (faserps) + reconstruct (batchreco.exe)
a handful of event types and compare the output against a stored golden
baseline (Tests/regression/golden/*.json). See docs/REGRESSION_TESTS.md
for the full design, what's verified-from-source vs. not yet run, and the
open items (a missing CI input sample, in particular).

Deterministic by construction: faserps.cc hardcodes its random seed, so a
run should be exactly reproducible given the same code, same input, same
--n-events, and single-threaded Geant4 simulation (run_faserps.py's own
default -- this script never overrides its --n-threads). A mismatch here
means something actually changed, not statistical noise -- see
docs/REGRESSION_TESTS.md's "Determinism this relies on" for why Geant4 MT
(a different, unrelated knob) is deliberately out of scope.

batchreco.exe's own -mt flag (TPORecoEvent::multiThread, parallelizing
Reconstruct3DPS_2's per-module voxel reconstruction across one
std::thread per detector module) is a separate axis, on by default here
-- see --multi-thread/--no-multi-thread. The existing golden/*.json were
recorded single-threaded, and the one comparison run with -mt on so far
passed against them, which is evidence (not proof) that the per-module
work is independent and race-free. If a future run ever FAILs under the
default but passes again with --no-multi-thread, that's a real lead (a
race in Reconstruct3DPS_2/reconstruct3DPS_module) worth reporting, not
noise to retry away.

STATUS: written but not yet run anywhere with ROOT/Geant4 available (see
docs/REGRESSION_TESTS.md). The first `--record` run on your Mac is the
real validation of this script, not this docstring.

Usage:
    fb                                    # build first, always
    python3 run_regression_tests.py --record          # record all cases' baselines
    python3 run_regression_tests.py                    # compare all cases against golden/
    python3 run_regression_tests.py --case muondis            # just one case
    python3 run_regression_tests.py --case muondis --record   # (re-)record just one case
    python3 run_regression_tests.py --case muondis --skip-faserps  # reco-only, reuse truth sample
    python3 run_regression_tests.py --case muondis --skip-reco     # summarize-only, reuse reco file
    python3 run_regression_tests.py --no-multi-thread  # old sequential batchreco.exe (-mt is on by default)
    python3 run_regression_tests.py --list             # list available case names and exit

Every run also writes each case's freshly-computed summary under
Tests/regression/results/<case>.json -- the actual numbers (mean_*/rms_*
etc.), regardless of PASS/FAIL/--record, for inspection without needing
--record (which would overwrite golden/'s baseline). Not a committed
baseline itself -- see .gitignore.

A plain comparison run (no --record) also writes
Tests/regression/results/comparison.csv: one row per field per case --
golden, current, abs_diff, rel_diff, status -- for every field, not just
the FAILs the terminal prints. rel_tol defaults to 1e-9 (near machine
precision, not a statistical tolerance): this suite expects *exact*
reproducibility (faserps' seed is pinned, see "Deterministic by
construction" below), so PASS means "the same to ~15 significant
figures", not "close enough given expected run-to-run noise" -- there
isn't supposed to be any.

faserps/batchreco.exe's own (often long) stdout/stderr no longer streams
to the terminal either -- it's written to
Tests/regression/work/<key>/logs/{faserps,batchreco*}.log instead, with
only a one-line "output -> <path>" pointer printed; a failure still
prints the last few dozen lines immediately so you're not stuck opening
the file to see what broke.
"""
import argparse
import csv
import json
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent

# Sibling script, same repo root as this one -- reused below so the
# neutrino case's input-file resolution (CVGENIE/Run<run>/ discovery)
# can't drift out of sync with run_faserps.py's own logic.
sys.path.insert(0, str(REPO_ROOT))
import run_faserps  # noqa: E402
REGRESSION_DIR = REPO_ROOT / "Tests" / "regression"
GOLDEN_DIR = REGRESSION_DIR / "golden"
WORK_DIR = REGRESSION_DIR / "work"
# Every case's freshly-computed summary, written on every run (--record or
# plain compare) -- unlike golden/, this isn't a committed baseline, it's
# just "what did the last run actually compute", so you can look at the
# real numbers (mean_*/rms_* etc.) without needing --record to see them,
# and without that --record clobbering golden/'s baseline. See .gitignore.
RESULTS_DIR = REGRESSION_DIR / "results"
DEFAULT_BUILD_DIR = REPO_ROOT / "build"

# Default relative tolerance for floating-point aggregates. Tight on
# purpose -- see docs/REGRESSION_TESTS.md's "Determinism" section for why
# this suite expects exact reproducibility, not statistical agreement.
DEFAULT_REL_TOL = 1e-9

# One entry per faserps simulation run. --muons and --muondis both
# hardcode run_number=999 on the C++ side (documented in
# docs/REGRESSION_TESTS.md), which is exactly why each group below gets
# its own isolated $FASERDATA under work/<key>/data -- without that,
# muons and muondis would silently overwrite each other's per-event
# output files.
#
# The "neutrino" group's nueCC/numuCC/nutauCC/nuNC golden cases are NOT
# separate batchreco.exe invocations filtered by --mask: faserps/TcalEvent
# never tag their own truth output with a mask in this (CVGENIE-based,
# unbiased-sample) workflow -- TPOEvent::SetEventMask() is only ever
# called from ConvertFASERMC.cc's conversion path, not from
# FASERG4/FASERCalProtoG4's ParticleManager.cc, which is what actually
# writes FASERG4-Tcalevent_*.root here. So every truth file's event_mask
# is 0, every masked Load_event() lookup in Batch/BatchReco.cc would miss
# for every single event (not just "most", as an earlier version of this
# comment assumed), and reconstructing with --mask nueCC/numuCC/nutauCC/
# nuNC would each produce an empty reco file, not "this flavor's events".
# Instead, "split_by_reaction" below means: run faserps once and
# batchreco.exe once, UNMASKED, over the whole unbiased sample, then split
# the single truth+reco summary into these four buckets in Python by each
# event's own actual interaction type (TPOEvent::reaction_desc()) -- see
# Tests/regression/summarize_output.py's --split-by-reaction and
# docs/REGRESSION_TESTS.md.
SIMULATION_GROUPS = [
    {
        "key": "neutrino",
        "faserps_args": [],  # --input-file is appended at run time -- see
                             # run_group()'s cvgenie_detector handling below
        "cvgenie_detector": "3DCAL",  # which $FASERDATA/CVGENIE/Run10000/ variant to use
        "run_number": 10000,  # from the input file's own stored TPOEvent.run_number -- see docs/REGRESSION_TESTS.md
        "n_events": 100,
        "split_by_reaction": {
            "neutrino_nueCC": "nueCC",
            "neutrino_numuCC": "numuCC",
            "neutrino_nutauCC": "nutauCC",
            "neutrino_nuNC": "nuNC",
        },
    },
    {
        "key": "muons",
        "faserps_args": ["--muons"],
        "run_number": 999,
        "n_events": 20,
        "reco_cases": [{"golden_name": "muons", "batchreco_args": []}],
    },
    {
        "key": "muondis",
        "faserps_args": ["--muondis"],
        "run_number": 999,
        "n_events": 20,
        "reco_cases": [{"golden_name": "muondis", "batchreco_args": []}],
    },
]


def group_golden_names(group):
    """Every golden case name a single SIMULATION_GROUPS entry produces,
    regardless of which shape it uses (split_by_reaction vs. reco_cases)."""
    if "split_by_reaction" in group:
        return list(group["split_by_reaction"].keys())
    return [rc["golden_name"] for rc in group["reco_cases"]]


def all_golden_names():
    return [name for group in SIMULATION_GROUPS for name in group_golden_names(group)]


_LOG_TAIL_LINES = 40


def run(cmd, *, env, label, capture=True, log_path=None):
    """capture=True (the default) buffers stdout/stderr so the caller can
    read result.stdout -- needed for the summarize_output.py step, whose
    JSON output we parse.

    log_path, when given, redirects the child's stdout+stderr (merged)
    straight to that file instead -- used for faserps/batchreco, whose
    own Geant4/BatchReco diagnostic output is voluminous (MDT geometry
    scans, per-event progress, ...) and was otherwise flooding the
    terminal on every run with nothing kept afterward. On failure, the
    last _LOG_TAIL_LINES lines are printed immediately (so you don't have
    to go open the file just to see what broke), with the full path
    alongside for the rest.

    capture=False (and log_path=None) lets the child inherit this
    process's own stdout/stderr instead, so its output streams live as it
    happens: capture_output=True doesn't just delay *our* printing, the
    child's own Python stdout also silently switches from line-buffered
    to block-buffered the moment it isn't a real terminal, so nothing
    appears until the whole subprocess exits or its buffer fills.
    Inheriting a real terminal fd avoids both at once -- no caller
    actually uses this mode any more (log_path replaced it for
    faserps/batchreco), but it's kept as the capture=False fallback."""
    print(f"[run_regression_tests] {label}: {' '.join(str(c) for c in cmd)}")
    if log_path is not None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"[run_regression_tests] {label}: output -> {log_path}")
        with open(log_path, "w") as f:
            result = subprocess.run(cmd, cwd=REPO_ROOT, env=env, stdout=f,
                                     stderr=subprocess.STDOUT, text=True)
    elif capture:
        result = subprocess.run(cmd, cwd=REPO_ROOT, env=env, capture_output=True, text=True)
    else:
        result = subprocess.run(cmd, cwd=REPO_ROOT, env=env, text=True)
    if result.returncode != 0:
        if log_path is not None:
            tail = log_path.read_text().splitlines()[-_LOG_TAIL_LINES:]
            print(f"error: {label} failed (exit {result.returncode}) -- "
                  f"last {len(tail)} line(s) of {log_path}:", file=sys.stderr)
            for line in tail:
                print(f"    {line}", file=sys.stderr)
        elif capture:
            print(result.stdout)
            print(result.stderr, file=sys.stderr)
        suffix = f" -- see {log_path}" if log_path is not None else ""
        sys.exit(f"error: {label} failed (exit {result.returncode}){suffix}")
    return result


def resolve_neutrino_input_file(group):
    """Finds the real, already-converted CVGENIE PO file for the neutrino
    case (group["run_number"]/group["cvgenie_detector"]) under the REAL
    $FASERDATA -- i.e. before run_group() below isolates FASERDATA to a
    scratch directory for this case's own simulation/reco output. Reuses
    run_faserps.py's own discovery helper rather than reimplementing the
    CVGENIE/Run<run>/ layout a second time, so the two can't drift apart.
    Exits with a clear, actionable error (same style run_faserps.py itself
    uses) if the sample isn't there -- this needs a real converted sample
    on disk, not something run_regression_tests.py can generate itself."""
    faserdata_dir = run_faserps._faserdata_dir()
    run_number = group["run_number"]
    detector = group["cvgenie_detector"]
    po_file = run_faserps.resolve_cvgenie_po_file(faserdata_dir, run_number, detector)
    if po_file is None:
        sys.exit(
            f"error: run_regression_tests.py's neutrino case needs a converted "
            f"GENIE sample at {faserdata_dir}/CVGENIE/Run{run_number}/ for detector "
            f"{detector}, but none was found.\n"
            f"       Run `python3 run_convertgenie.py --run {run_number}` first (see "
            f"docs/HOWTO.md, \"run_convertgenie.py\") to produce it."
        )
    return po_file


def run_group(group, *, build_dir, python_exe, skip_faserps=False, skip_reco=False, multi_thread=False):
    """Runs one faserps simulation and every batchreco/summarize pass that
    reads it, isolated in its own $FASERDATA. Returns {golden_name: summary_dict}.

    skip_faserps=True reuses whatever truth sample is already sitting in
    this group's work/<key>/data/faserG4 from an earlier (non-skipped) run
    instead of re-simulating -- useful while iterating on reconstruction-
    only code (BatchReco.cc, TPORecoEvent, ...), where re-running Geant4
    every time is pure overhead and the truth sample hasn't changed. It's
    on the caller to know that's actually true: this doesn't hash/compare
    faserps_args or detect a stale sample, it just checks the truth files
    are there at all (see the error below if they aren't).

    skip_reco=True is the same idea one step further down the pipeline:
    reuses whatever reco file is already sitting in this group's
    work/<key>/data/batch from an earlier (non-skipped) run instead of
    re-running batchreco.exe -- useful while iterating on summarize_output.py
    or the golden-comparison logic itself, where the reco file hasn't
    changed. Independent of skip_faserps: combine both to jump straight to
    the summarize step, reusing both the truth sample and the reco file.

    multi_thread=True passes batchreco.exe's own -mt flag through
    run_batchreco.py's --multi-thread (see Batch/BatchReco.cc,
    TPORecoEvent::multiThread -- it only affects Reconstruct3DPS_2's
    per-module voxel reconstruction, parallelized with one std::thread per
    detector module instead of a sequential loop). Ignored when
    skip_reco=True, since then batchreco.exe isn't run at all. Has no
    effect on faserps/the truth sample. Useful for checking -mt doesn't
    change results (it shouldn't, if Reconstruct3DPS_2's per-module work
    is actually independent) against the existing (single-threaded)
    golden files -- no reason to --record a second set of goldens just to
    check that."""
    key = group["key"]
    faserps_args = list(group["faserps_args"])
    if "cvgenie_detector" in group:
        # File-input mode (not --muons/--muondis): resolve the real sample
        # now, from the real $FASERDATA, before the isolation below hides
        # it -- see resolve_neutrino_input_file()'s docstring.
        faserps_args += ["--input-file", str(resolve_neutrino_input_file(group))]

    faserdata = WORK_DIR / key / "data"
    faserdata.mkdir(parents=True, exist_ok=True)
    # Raw stdout+stderr of each faserps/batchreco.exe subprocess call below,
    # one file per step -- see run()'s log_path. Kept (not printed to the
    # terminal) so a normal run's output stays to the PASS/FAIL lines,
    # while the full Geant4/BatchReco diagnostic firehose is still there
    # to open if something needs triaging.
    logs_dir = WORK_DIR / key / "logs"
    env = dict(os.environ)
    env["FASERDATA"] = str(faserdata)

    if skip_faserps:
        truth_dir = faserdata / "faserG4"
        if not truth_dir.is_dir() or not any(truth_dir.glob(f"FASERG4-Tcalevent_{group['run_number']}_*.root")):
            sys.exit(
                f"error: --skip-faserps given for '{key}', but no truth files found at "
                f"{truth_dir} (FASERG4-Tcalevent_{group['run_number']}_*.root). Run once "
                f"without --skip-faserps first to produce them."
            )
        print(f"[run_regression_tests] {key}: faserps SKIPPED (--skip-faserps) -- "
              f"reusing truth sample already in {truth_dir}")
    else:
        run([
            python_exe, "run_faserps.py",
            *faserps_args,
            "--n-events", str(group["n_events"]),
            "--build-dir", str(build_dir),
        ], env=env, label=f"{key}: faserps", log_path=logs_dir / "faserps.log")

    run_number = group["run_number"]
    n_events = group["n_events"]

    if "split_by_reaction" in group:
        # Single unbiased, unmasked reconstruction pass -- faserps/TcalEvent
        # never tag their own truth output with a mask in this workflow (see
        # SIMULATION_GROUPS' comment above), so there is exactly one
        # batchreco.exe run here, and the nueCC/numuCC/nutauCC/nuNC split
        # happens afterward in Python, from each event's own actual
        # interaction type (TPOEvent::reaction_desc()), not from a
        # mask-tagged filename.
        reco_file = faserdata / "batch" / f"Batch-TPORecevent_{run_number}_0_{n_events}.root"
        if skip_reco:
            if not reco_file.is_file():
                sys.exit(
                    f"error: --skip-reco given for '{key}', but no reco file found at "
                    f"{reco_file}. Run once without --skip-reco first to produce it."
                )
            print(f"[run_regression_tests] {key}: batchreco SKIPPED (--skip-reco) -- "
                  f"reusing reco file already at {reco_file}")
        else:
            batchreco_cmd = [
                python_exe, "run_batchreco.py",
                "--run", str(run_number),
                "--max-event", str(n_events),
                "--build-dir", str(build_dir),
            ]
            if multi_thread:
                batchreco_cmd.append("--multi-thread")
            run(batchreco_cmd, env=env, label=f"{key}: batchreco", log_path=logs_dir / "batchreco.log")

        summarize_cmd = [
            python_exe, str(REGRESSION_DIR / "summarize_output.py"),
            "--truth-dir", str(faserdata / "faserG4"),
            "--run", str(run_number),
            "--n-events", str(n_events),
            "--reco-file", str(reco_file),
            "--split-by-reaction",
        ]
        result = run(summarize_cmd, env=env, label=f"{key}: summarize")
        try:
            by_reaction = json.loads(result.stdout)
        except json.JSONDecodeError as e:
            sys.exit(f"error: summarize_output.py for {key} did not print valid JSON: {e}\n"
                      f"stdout was:\n{result.stdout}")

        results = {}
        for golden_name, reaction in group["split_by_reaction"].items():
            summary = by_reaction.get(reaction, {})
            summary.setdefault("truth", {"n_events": 0})
            summary.setdefault("reco", {"n_events_reconstructed": 0})
            summary["meta"] = {
                "case": golden_name,
                "reaction": reaction,
                "faserps_args": group["faserps_args"],
                "run_number": run_number,
                "n_events_simulated": n_events,
            }
            results[golden_name] = summary
        return results

    # Legacy path (--muons/--muondis): a single, already-homogeneous sample
    # -- nothing to split by reaction, one reco_case per golden name.
    results = {}
    for reco_case in group["reco_cases"]:
        mask = None
        if "--mask" in reco_case["batchreco_args"]:
            mask = reco_case["batchreco_args"][reco_case["batchreco_args"].index("--mask") + 1]
        reco_filename = f"Batch-TPORecevent_{run_number}_0_{n_events}"
        if mask:
            reco_filename += f"_{mask}"
        reco_filename += ".root"
        reco_file = faserdata / "batch" / reco_filename

        if skip_reco:
            if not reco_file.is_file():
                sys.exit(
                    f"error: --skip-reco given for '{reco_case['golden_name']}', but no reco "
                    f"file found at {reco_file}. Run once without --skip-reco first to produce it."
                )
            print(f"[run_regression_tests] {reco_case['golden_name']}: batchreco SKIPPED "
                  f"(--skip-reco) -- reusing reco file already at {reco_file}")
        else:
            batchreco_cmd = [
                python_exe, "run_batchreco.py",
                "--run", str(run_number),
                "--max-event", str(n_events),
                *reco_case["batchreco_args"],
                "--build-dir", str(build_dir),
            ]
            if multi_thread:
                batchreco_cmd.append("--multi-thread")
            run(batchreco_cmd, env=env, label=f"{reco_case['golden_name']}: batchreco",
                log_path=logs_dir / f"batchreco_{reco_case['golden_name']}.log")

        summarize_cmd = [
            python_exe, str(REGRESSION_DIR / "summarize_output.py"),
            "--truth-dir", str(faserdata / "faserG4"),
            "--run", str(run_number),
            "--n-events", str(n_events),
            "--reco-file", str(reco_file),
        ]
        result = run(summarize_cmd, env=env, label=f"{reco_case['golden_name']}: summarize")
        try:
            summary = json.loads(result.stdout)
        except json.JSONDecodeError as e:
            sys.exit(f"error: summarize_output.py for {reco_case['golden_name']} did not print valid JSON: {e}\n"
                      f"stdout was:\n{result.stdout}")

        summary["meta"] = {
            "case": reco_case["golden_name"],
            "faserps_args": group["faserps_args"],
            "batchreco_args": reco_case["batchreco_args"],
            "run_number": run_number,
            "n_events_simulated": n_events,
        }
        results[reco_case["golden_name"]] = summary
    return results


def _compare_value(cur_val, gold_val, rel_tol):
    """Single-field tolerance check, shared by compare() (reports only the
    FAILs, as short strings) and diff_rows() (reports every field, FAIL or
    not) -- one place for the actual numeric comparison so the two can't
    quietly drift apart and disagree about what counts as a match.

    Returns (abs_diff, rel_diff, ok). abs_diff/rel_diff are None when the
    two aren't both present numeric values (missing/appeared-disappeared,
    or a non-numeric field) -- "difference" isn't meaningful there, only
    equality is. rel_diff is also None when gold_val == 0 (relative
    difference is undefined; ok instead falls back to an absolute check
    against rel_tol, same as compare() always did)."""
    if cur_val is None or gold_val is None:
        return None, None, cur_val == gold_val
    if isinstance(cur_val, (int, float)) and isinstance(gold_val, (int, float)):
        abs_diff = cur_val - gold_val
        if gold_val == 0:
            return abs_diff, None, abs(cur_val) < rel_tol
        rel_diff = abs_diff / abs(gold_val)
        return abs_diff, rel_diff, abs(rel_diff) <= rel_tol
    return None, None, cur_val == gold_val


def _diff_fields(current, golden):
    """Every (section, key, cur_val, gold_val) triple across both dicts'
    truth/reco sub-dicts -- the field union compare() and diff_rows() both
    walk, pulled out once so they can't disagree about which fields to
    look at either."""
    for section in ("truth", "reco"):
        cur_section = current.get(section, {})
        gold_section = golden.get(section, {})
        for key in sorted(set(cur_section) | set(gold_section)):
            yield section, key, cur_section.get(key), gold_section.get(key)


def compare(golden_name, current, golden, rel_tol):
    """Compares the "truth"/"reco" sub-dicts field by field. Returns a list
    of human-readable mismatch strings; empty means a pass."""
    mismatches = []
    for section, key, cur_val, gold_val in _diff_fields(current, golden):
        _, _, ok = _compare_value(cur_val, gold_val, rel_tol)
        if ok:
            continue
        if cur_val is None or gold_val is None:
            mismatches.append(f"{section}.{key}: golden={gold_val!r} current={cur_val!r} (field appeared/disappeared)")
        else:
            mismatches.append(f"{section}.{key}: golden={gold_val} current={cur_val}")
    return mismatches


def diff_rows(golden_name, current, golden, rel_tol):
    """Like compare(), but yields one dict per field -- PASS and FAIL alike
    -- instead of only the mismatch strings: {case, section, field,
    golden, current, abs_diff, rel_diff, status}. Meant to be collected
    across every case into one table (see main()'s comparison.csv) so you
    can see each value against its golden and by how much, not just
    whether the case as a whole passed."""
    for section, key, cur_val, gold_val in _diff_fields(current, golden):
        abs_diff, rel_diff, ok = _compare_value(cur_val, gold_val, rel_tol)
        yield {
            "case": golden_name,
            "section": section,
            "field": key,
            "golden": gold_val,
            "current": cur_val,
            "abs_diff": abs_diff,
            "rel_diff": rel_diff,
            "status": "PASS" if ok else "FAIL",
        }


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--case", action="append", dest="cases", default=None,
                         help="Restrict to this golden case name (repeatable). Default: all cases. "
                              "Use --list to see the names.")
    parser.add_argument("--record", action="store_true",
                         help="Write the golden file(s) for the selected case(s) instead of comparing.")
    parser.add_argument("--rel-tol", type=float, default=DEFAULT_REL_TOL,
                         help=f"Relative tolerance for floating-point aggregates (default: {DEFAULT_REL_TOL}).")
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD_DIR,
                         help=f"CMake build directory (default: {DEFAULT_BUILD_DIR}).")
    parser.add_argument("--python", default=os.environ.get("FASER_PYTHON", sys.executable),
                         help="Python interpreter to use for run_faserps.py/run_batchreco.py/"
                              "summarize_output.py subprocesses (default: $FASER_PYTHON if set "
                              "-- see setup.sh -- otherwise this interpreter). Needs a "
                              "PyROOT-enabled interpreter matching the Python your ROOT build was "
                              "linked against (a conda/venv-shadowed \"python3\" is a common way "
                              "to get this wrong -- PyROOT's import fails loudly with a "
                              "major.minor mismatch if so, see docs/REGRESSION_TESTS.md).")
    parser.add_argument("--skip-faserps", action="store_true",
                         help="Skip the faserps simulation step and reuse whatever truth sample "
                              "is already sitting in work/<key>/data/faserG4 from an earlier run "
                              "of the same case(s). Fails with a clear error if that sample isn't "
                              "there yet. Useful while iterating on reconstruction-only code -- "
                              "no reason to re-run Geant4 if the truth sample hasn't changed.")
    parser.add_argument("--skip-reco", action="store_true",
                         help="Skip the batchreco.exe step and reuse whatever reco file is "
                              "already sitting in work/<key>/data/batch from an earlier run of "
                              "the same case(s). Fails with a clear error if that file isn't "
                              "there yet. Useful while iterating on summarize_output.py or the "
                              "golden-comparison logic -- no reason to re-run reconstruction if "
                              "the reco file hasn't changed. Independent of --skip-faserps; pass "
                              "both to jump straight to the summarize step.")
    parser.add_argument("--multi-thread", action="store_true", dest="multi_thread", default=True,
                         help="Pass batchreco.exe's -mt flag (via run_batchreco.py's "
                              "--multi-thread) -- parallelizes Reconstruct3DPS_2's per-module "
                              "voxel reconstruction across one std::thread per detector module "
                              "instead of a sequential loop. Ignored together with --skip-reco. "
                              "On by default; pass --no-multi-thread for the old sequential "
                              "behavior (e.g. to isolate whether a mismatch is -mt-related).")
    parser.add_argument("--no-multi-thread", action="store_false", dest="multi_thread",
                         help="Opposite of --multi-thread -- run batchreco.exe single-threaded.")
    parser.add_argument("--list", action="store_true", help="List available case names and exit.")
    args = parser.parse_args()
    if args.list:
        for name in all_golden_names():
            print(name)
        sys.exit(0)
    if args.cases:
        unknown = set(args.cases) - set(all_golden_names())
        if unknown:
            parser.error(f"unknown case(s): {', '.join(sorted(unknown))} -- see --list")
    return args


def main():
    args = parse_args()
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    selected_names = set(args.cases) if args.cases else set(all_golden_names())
    groups_to_run = [g for g in SIMULATION_GROUPS if selected_names & set(group_golden_names(g))]

    any_failure = False
    all_diff_rows = []
    for group in groups_to_run:
        results = run_group(group, build_dir=args.build_dir, python_exe=args.python,
                            skip_faserps=args.skip_faserps, skip_reco=args.skip_reco,
                            multi_thread=args.multi_thread)
        for golden_name, summary in results.items():
            if golden_name not in selected_names:
                continue

            results_path = RESULTS_DIR / f"{golden_name}.json"
            with open(results_path, "w") as f:
                json.dump(summary, f, indent=2, sort_keys=True)
                f.write("\n")

            golden_path = GOLDEN_DIR / f"{golden_name}.json"
            if args.record:
                with open(golden_path, "w") as f:
                    json.dump(summary, f, indent=2, sort_keys=True)
                    f.write("\n")
                print(f"[run_regression_tests] {golden_name}: recorded -> {golden_path}")
                continue

            if not golden_path.is_file():
                print(f"[run_regression_tests] {golden_name}: FAIL -- no golden file at {golden_path} "
                      f"(run with --record first)")
                any_failure = True
                continue
            with open(golden_path) as f:
                golden = json.load(f)
            all_diff_rows.extend(diff_rows(golden_name, summary, golden, args.rel_tol))
            mismatches = compare(golden_name, summary, golden, args.rel_tol)
            if mismatches:
                print(f"[run_regression_tests] {golden_name}: FAIL")
                for m in mismatches:
                    print(f"    {m}")
                any_failure = True
            else:
                print(f"[run_regression_tests] {golden_name}: PASS")

    print(f"[run_regression_tests] full per-case results (the actual computed numbers, "
          f"not just PASS/FAIL) written under {RESULTS_DIR}/")

    if all_diff_rows:
        comparison_path = RESULTS_DIR / "comparison.csv"
        with open(comparison_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["case", "section", "field", "golden",
                                                     "current", "abs_diff", "rel_diff", "status"])
            writer.writeheader()
            writer.writerows(all_diff_rows)
        n_fail = sum(1 for row in all_diff_rows if row["status"] == "FAIL")
        print(f"[run_regression_tests] field-by-field comparison against golden/ for every "
              f"case ({len(all_diff_rows)} fields, {n_fail} FAIL) written to {comparison_path}")

    if args.record:
        return 0
    return 1 if any_failure else 0


if __name__ == "__main__":
    sys.exit(main())
