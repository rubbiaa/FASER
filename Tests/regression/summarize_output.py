#!/usr/bin/env python3
"""
summarize_output.py -- read a FASERG4 truth output directory (one
FASERG4-Tcalevent_<run>_<event>.root file per event) and/or one BatchReco
output file, and print a small JSON summary of aggregate quantities to
stdout. This is the PyROOT half of the simulate->reconstruct->compare
regression suite (see docs/REGRESSION_TESTS.md) -- run_regression_tests.py
calls it as a subprocess per case and loads its stdout as the case's
"truth"/"reco" summary.

Deliberately kept to aggregate numbers (counts, means), not a full
per-event dump: the golden JSON files this feeds are meant to be readable
in a PR diff.

Loads CoreUtils/libCoreUtilsDict.so via PyROOT, the same dictionary
CoreUtils/dumpReco.C already uses from a plain ROOT macro (gSystem->Load +
TChain("RecoEvent") + SetBranchAddress("TPORecoEvent", ...)) -- this
script is a PyROOT translation of that same, already-working pattern, not
a new one. See docs/REGRESSION_TESTS.md's "Status" note: this has NOT
been run yet (no ROOT/PyROOT available where this was written), so a
first --record run on your Mac is the real validation step.

Usage:
    # Truth only (summarize FASERG4's own output for a run):
    python3 summarize_output.py --truth-dir /path/to/FASERDATA/faserG4 \\
        --run 999 --n-events 20

    # Truth + reco:
    python3 summarize_output.py --truth-dir /path/to/FASERDATA/faserG4 \\
        --run 10000 --n-events 100 \\
        --reco-file /path/to/FASERDATA/batch/Batch-TPORecevent_10000_0_100_numuCC.root

    # Reco only:
    python3 summarize_output.py \\
        --reco-file /path/to/FASERDATA/batch/Batch-TPORecevent_999_0_20.root
"""
import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_DICT_PATH = REPO_ROOT / "CoreUtils" / "libCoreUtilsDict.so"


def load_dictionary(dict_path: Path):
    """Loads CoreUtils' ROOT dictionary so PyROOT can read TPOEvent/
    TPORecoEvent objects out of a tree, exactly like CoreUtils/dumpReco.C
    does via gSystem->Load(). Built automatically by the normal CMake
    build (CoreUtils/CMakeLists.txt's CoreUtilsDict target) -- if this is
    missing, build first."""
    if not dict_path.is_file():
        sys.exit(
            f"error: dictionary not found: {dict_path}\n"
            f"       (build first, e.g. `cmake --build {REPO_ROOT / 'build'} -j` "
            f"or `fb` -- CoreUtils/libCoreUtilsDict.so is produced as part of "
            f"the normal build, see CoreUtils/CMakeLists.txt)"
        )
    import ROOT  # noqa: E402 -- imported lazily so --help works without PyROOT installed
    ok = ROOT.gSystem.Load(str(dict_path))
    if ok < 0:
        sys.exit(f"error: ROOT.gSystem.Load({dict_path}) failed (return code {ok})")
    return ROOT


def mean(values):
    return sum(values) / len(values) if values else None


def summarize_truth(ROOT, truth_dir: Path, run: int, n_events: int):
    """Reads FASERG4's own per-event output files (FASERG4-Tcalevent_<run>_
    <event>.root, one per event -- see CoreUtils/TcalEvent.cc), tree
    "calEvent", branch "event" (a TPOEvent). Missing event indices are
    skipped, not treated as an error -- a small/sparse run is expected for
    some cases (see docs/REGRESSION_TESTS.md's mask-coverage note)."""
    evis, n_particles, nuE, q2, xbj, yinel = [], [], [], [], [], []
    found = 0
    for ievent in range(n_events):
        f = truth_dir / f"FASERG4-Tcalevent_{run}_{ievent}.root"
        if not f.is_file():
            continue
        tfile = ROOT.TFile.Open(str(f))
        if not tfile or tfile.IsZombie():
            print(f"warning: could not open {f}", file=sys.stderr)
            continue
        tree = tfile.Get("calEvent")
        if not tree:
            print(f"warning: {f} has no 'calEvent' tree", file=sys.stderr)
            tfile.Close()
            continue
        po_event = ROOT.TPOEvent()
        tree.SetBranchAddress("event", ROOT.AddressOf(po_event) if hasattr(ROOT, "AddressOf") else po_event)
        # PyROOT can usually bind a branch of object type directly without
        # AddressOf(); the above covers both bindings depending on PyROOT
        # version. If neither works on your setup, see the fallback note
        # in docs/REGRESSION_TESTS.md and adjust this one line.
        tree.GetEntry(0)
        found += 1
        evis.append(po_event.Evis)
        n_particles.append(po_event.n_particles())
        # DIS kinematics are only meaningful for interaction-type events
        # (Q2 > 0 is the simplest available signal that kinematics_event()
        # actually ran for this event); skip them otherwise rather than
        # padding the mean with structural zeros.
        if po_event.Q2 > 0:
            nuE.append(po_event.nuE)
            q2.append(po_event.Q2)
            xbj.append(po_event.xBj)
            yinel.append(po_event.yInel)
        tfile.Close()

    summary = {
        "n_events": found,
        "mean_Evis": mean(evis),
        "mean_n_particles": mean(n_particles),
    }
    if q2:
        summary["mean_nuE"] = mean(nuE)
        summary["mean_Q2"] = mean(q2)
        summary["mean_xBj"] = mean(xbj)
        summary["mean_yInel"] = mean(yinel)
    return summary


def summarize_reco(ROOT, reco_file: Path):
    """Reads one BatchReco output file (Batch-TPORecevent_<run>_<min>_
    <max>[_<mask>].root -- see Batch/BatchReco.cc), tree "RecoEvent",
    branch "TPORecoEvent" -- the exact pattern CoreUtils/dumpReco.C already
    uses (TChain + SetBranchAddress), just via PyROOT instead of a ROOT
    macro."""
    if not reco_file.is_file():
        sys.exit(f"error: reco file not found: {reco_file}")

    chain = ROOT.TChain("RecoEvent")
    chain.Add(str(reco_file))
    n_entries = chain.GetEntries()
    if n_entries == 0:
        return {"n_events_reconstructed": 0}

    reco = ROOT.TPORecoEvent()
    chain.SetBranchAddress("TPORecoEvent", ROOT.AddressOf(reco) if hasattr(ROOT, "AddressOf") else reco)

    n_porecs, n_tktracks, n_tkvertices, n_mutracks = [], [], [], []
    total_evis, total_ecompensated = [], []

    for ientry in range(n_entries):
        chain.GetEntry(ientry)
        porecs = reco.fPORecs
        n_porecs.append(len(porecs))
        n_tktracks.append(len(reco.fTKTracks))
        n_tkvertices.append(len(reco.fTKVertices))
        n_mutracks.append(len(reco.fMuTracks))
        evis_sum = 0.0
        ecomp_sum = 0.0
        for porec in porecs:
            evis_sum += porec.TotalEvis()
            ecomp_sum += porec.fTotal.Ecompensated
        total_evis.append(evis_sum)
        total_ecompensated.append(ecomp_sum)

    return {
        "n_events_reconstructed": n_entries,
        "mean_n_PORecs": mean(n_porecs),
        "mean_n_TKTracks": mean(n_tktracks),
        "mean_n_TKVertices": mean(n_tkvertices),
        "mean_n_MuTracks": mean(n_mutracks),
        "mean_total_Evis_reco": mean(total_evis),
        "mean_total_Ecompensated": mean(total_ecompensated),
    }


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dict-path", type=Path, default=DEFAULT_DICT_PATH,
                         help=f"Path to CoreUtils' built ROOT dictionary (default: {DEFAULT_DICT_PATH})")
    parser.add_argument("--truth-dir", type=Path, default=None,
                         help="Directory containing FASERG4-Tcalevent_<run>_<event>.root files "
                              "(typically $FASERDATA/faserG4). Omit to skip the truth summary.")
    parser.add_argument("--run", type=int, default=None, help="Run number, required with --truth-dir.")
    parser.add_argument("--n-events", type=int, default=None,
                         help="Number of event indices (0..N-1) to look for, required with --truth-dir.")
    parser.add_argument("--reco-file", type=Path, default=None,
                         help="Path to one Batch-TPORecevent_*.root file. Omit to skip the reco summary.")
    args = parser.parse_args()
    if args.truth_dir and (args.run is None or args.n_events is None):
        parser.error("--truth-dir requires --run and --n-events")
    if not args.truth_dir and not args.reco_file:
        parser.error("give at least one of --truth-dir or --reco-file")
    return args


def main():
    args = parse_args()
    ROOT = load_dictionary(args.dict_path)
    ROOT.gROOT.SetBatch(True)

    result = {}
    if args.truth_dir:
        result["truth"] = summarize_truth(ROOT, args.truth_dir, args.run, args.n_events)
    if args.reco_file:
        result["reco"] = summarize_reco(ROOT, args.reco_file)

    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
