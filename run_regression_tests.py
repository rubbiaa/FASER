#!/usr/bin/env python3
"""
run_regression_tests.py -- simulate (faserps) + reconstruct (batchreco.exe)
a handful of event types and compare the output against a stored golden
baseline (Tests/regression/golden/*.json). See docs/REGRESSION_TESTS.md
for the full design, what's verified-from-source vs. not yet run, and the
open items (a missing CI input sample, in particular).

Deterministic by construction: faserps.cc hardcodes its random seed, so a
run should be exactly reproducible given the same code, same input, same
--n-events, and single-threaded mode (this script's default). A mismatch
here means something actually changed, not statistical noise -- see
docs/REGRESSION_TESTS.md's "Determinism this relies on".

STATUS: written but not yet run anywhere with ROOT/Geant4 available (see
docs/REGRESSION_TESTS.md). The first `--record` run on your Mac is the
real validation of this script, not this docstring.

Usage:
    fb                                    # build first, always
    python3 run_regression_tests.py --record          # record all cases' baselines
    python3 run_regression_tests.py                    # compare all cases against golden/
    python3 run_regression_tests.py --case muondis            # just one case
    python3 run_regression_tests.py --case muondis --record   # (re-)record just one case
    python3 run_regression_tests.py --list             # list available case names and exit
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
REGRESSION_DIR = REPO_ROOT / "Tests" / "regression"
GOLDEN_DIR = REGRESSION_DIR / "golden"
WORK_DIR = REGRESSION_DIR / "work"
DEFAULT_BUILD_DIR = REPO_ROOT / "build"

# Default relative tolerance for floating-point aggregates. Tight on
# purpose -- see docs/REGRESSION_TESTS.md's "Determinism" section for why
# this suite expects exact reproducibility, not statistical agreement.
DEFAULT_REL_TOL = 1e-9

# One entry per faserps simulation run. nueCC/numuCC/nutauCC/nuNC share a
# single "neutrino" simulation and differ only in which --mask batchreco.py
# reconstructs with -- see docs/REGRESSION_TESTS.md's "run_number per case"
# section for why (mask filtering happens at reconstruction time, not
# simulation time). --muons and --muondis both hardcode run_number=999 on
# the C++ side (also documented there), which is exactly why each group
# below gets its own isolated $FASERDATA under work/<key>/data -- without
# that, muons and muondis would silently overwrite each other's per-event
# output files.
SIMULATION_GROUPS = [
    {
        "key": "neutrino",
        "faserps_args": [],  # default mode: reads the default input sample
        "run_number": 10000,  # from the input file's own stored TPOEvent.run_number -- see docs/REGRESSION_TESTS.md
        "n_events": 100,
        "reco_cases": [
            {"golden_name": "neutrino_nueCC", "batchreco_args": ["--mask", "nueCC"]},
            {"golden_name": "neutrino_numuCC", "batchreco_args": ["--mask", "numuCC"]},
            {"golden_name": "neutrino_nutauCC", "batchreco_args": ["--mask", "nutauCC"]},
            {"golden_name": "neutrino_nuNC", "batchreco_args": ["--mask", "nuNC"]},
        ],
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


def all_golden_names():
    return [rc["golden_name"] for group in SIMULATION_GROUPS for rc in group["reco_cases"]]


def run(cmd, *, env, label):
    print(f"[run_regression_tests] {label}: {' '.join(str(c) for c in cmd)}")
    result = subprocess.run(cmd, cwd=REPO_ROOT, env=env, capture_output=True, text=True)
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        sys.exit(f"error: {label} failed (exit {result.returncode})")
    return result


def run_group(group, *, build_dir, python_exe):
    """Runs one faserps simulation and every batchreco/summarize pass that
    reads it, isolated in its own $FASERDATA. Returns {golden_name: summary_dict}."""
    key = group["key"]
    faserdata = WORK_DIR / key / "data"
    faserdata.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["FASERDATA"] = str(faserdata)

    run([
        python_exe, "run_faserps.py",
        *group["faserps_args"],
        "--n-events", str(group["n_events"]),
        "--build-dir", str(build_dir),
    ], env=env, label=f"{key}: faserps")

    results = {}
    for reco_case in group["reco_cases"]:
        run_number = group["run_number"]
        n_events = group["n_events"]
        run([
            python_exe, "run_batchreco.py",
            "--run", str(run_number),
            "--max-event", str(n_events),
            *reco_case["batchreco_args"],
            "--build-dir", str(build_dir),
        ], env=env, label=f"{reco_case['golden_name']}: batchreco")

        mask = None
        if "--mask" in reco_case["batchreco_args"]:
            mask = reco_case["batchreco_args"][reco_case["batchreco_args"].index("--mask") + 1]
        reco_filename = f"Batch-TPORecevent_{run_number}_0_{n_events}"
        if mask:
            reco_filename += f"_{mask}"
        reco_filename += ".root"
        reco_file = faserdata / "batch" / reco_filename

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


def compare(golden_name, current, golden, rel_tol):
    """Compares the "truth"/"reco" sub-dicts field by field. Returns a list
    of human-readable mismatch strings; empty means a pass."""
    mismatches = []
    for section in ("truth", "reco"):
        cur_section = current.get(section, {})
        gold_section = golden.get(section, {})
        keys = set(cur_section) | set(gold_section)
        for key in sorted(keys):
            cur_val = cur_section.get(key)
            gold_val = gold_section.get(key)
            if cur_val is None or gold_val is None:
                if cur_val != gold_val:
                    mismatches.append(f"{section}.{key}: golden={gold_val!r} current={cur_val!r} (field appeared/disappeared)")
                continue
            if isinstance(cur_val, (int, float)) and isinstance(gold_val, (int, float)):
                if gold_val == 0:
                    ok = abs(cur_val) < rel_tol
                else:
                    ok = abs(cur_val - gold_val) / abs(gold_val) <= rel_tol
                if not ok:
                    mismatches.append(f"{section}.{key}: golden={gold_val} current={cur_val}")
            elif cur_val != gold_val:
                mismatches.append(f"{section}.{key}: golden={gold_val!r} current={cur_val!r}")
    return mismatches


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
    parser.add_argument("--python", default=sys.executable,
                         help="Python interpreter to use for run_faserps.py/run_batchreco.py/"
                              "summarize_output.py subprocesses (default: this interpreter). "
                              "Use a PyROOT-enabled interpreter if it differs from the default one.")
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

    selected_names = set(args.cases) if args.cases else set(all_golden_names())
    groups_to_run = [g for g in SIMULATION_GROUPS if any(rc["golden_name"] in selected_names for rc in g["reco_cases"])]

    any_failure = False
    for group in groups_to_run:
        results = run_group(group, build_dir=args.build_dir, python_exe=args.python)
        for golden_name, summary in results.items():
            if golden_name not in selected_names:
                continue
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
            mismatches = compare(golden_name, summary, golden, args.rel_tol)
            if mismatches:
                print(f"[run_regression_tests] {golden_name}: FAIL")
                for m in mismatches:
                    print(f"    {m}")
                any_failure = True
            else:
                print(f"[run_regression_tests] {golden_name}: PASS")

    if args.record:
        return 0
    return 1 if any_failure else 0


if __name__ == "__main__":
    sys.exit(main())
