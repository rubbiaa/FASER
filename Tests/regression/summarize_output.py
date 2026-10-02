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
a new one.

Two output modes:

  * Flat (default): one {"truth": {...}, "reco": {...}} dict, aggregated
    over every event found. This is what --muons/--muondis use -- there's
    nothing to split, every event in those samples is the same "type".

  * --split-by-reaction: faserps/TcalEvent never tag truth files with a
    mask -- SetEventMask() is only ever called from ConvertFASERMC.cc's
    conversion path, not from FASERG4/FASERCalProtoG4's ParticleManager.cc
    (which is what actually writes FASERG4-Tcalevent_*.root for the
    CVGENIE-based "neutrino" regression case). So a truth/reco sample
    produced by the normal faserps->batchreco.exe pipeline is always
    unbiased/unmasked, and splitting it into nueCC/numuCC/nutauCC/nuNC
    has to happen here, in Python, by reading each event's *actual*
    interaction type -- TPOEvent::reaction_desc() -- not by filename or
    by re-running batchreco.exe once per flavor with --mask (which can
    never match anything in this workflow: see docs/REGRESSION_TESTS.md).
    Prints {"<reaction>": {"truth": {...}, "reco": {...}}, ...}, one entry
    per distinct reaction_desc() string actually present in the sample
    (typically "nueCC", "numuCC", "nutauCC", "nuNC", but also
    "antinueCC"/"antinumuCC"/"antinutauCC"/"ES" if the sample has any).

This script needs to run under a PyROOT-enabled Python matching the major.
minor version your ROOT build was linked against — `import ROOT` fails
loudly, naming both versions, if invoked under the wrong one (a conda/venv
`python3` easily shadows the right one without you noticing). Run it with
run_regression_tests.py (which resolves this via `$FASER_PYTHON`/`--python`
— see docs/REGRESSION_TESTS.md's "Running it"), or invoke the matching
interpreter directly if calling this script standalone.

Usage:
    # Truth only (summarize FASERG4's own output for a run):
    python3 summarize_output.py --truth-dir /path/to/FASERDATA/faserG4 \\
        --run 999 --n-events 20

    # Truth + reco, flat:
    python3 summarize_output.py --truth-dir /path/to/FASERDATA/faserG4 \\
        --run 999 --n-events 20 \\
        --reco-file /path/to/FASERDATA/batch/Batch-TPORecevent_999_0_20.root

    # Truth + reco, split by each event's actual interaction type:
    python3 summarize_output.py --truth-dir /path/to/FASERDATA/faserG4 \\
        --run 10000 --n-events 100 \\
        --reco-file /path/to/FASERDATA/batch/Batch-TPORecevent_10000_0_100.root \\
        --split-by-reaction

    # Reco only:
    python3 summarize_output.py \\
        --reco-file /path/to/FASERDATA/batch/Batch-TPORecevent_999_0_20.root
"""
import argparse
import contextlib
import json
import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
# CMake names the shared library by platform (.dylib on macOS, .so
# elsewhere -- see cmake/Externals.cmake's _faser_shlib_suffix), so a
# hardcoded ".so" here silently never matches on a Mac. Same convention
# ConvertNPZ/root_io.py's dictionary_path() already uses for the same
# file, kept in sync with it rather than reinventing it independently.
_DICT_SUFFIX = ".dylib" if sys.platform == "darwin" else ".so"
DEFAULT_DICT_PATH = REPO_ROOT / "CoreUtils" / f"libCoreUtilsDict{_DICT_SUFFIX}"


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


@contextlib.contextmanager
def _silence_native_stdout():
    """Temporarily redirects OS-level fd 1 (stdout) to /dev/null.

    Needed because this script's contract is "stdout is exactly one JSON
    document" (run_regression_tests.py parses the whole of it as one) --
    but PyROOT/GenFit's own C++-side std::cout calls write straight to
    the real file descriptor, bypassing Python's sys.stdout entirely, so
    patching sys.stdout (e.g. contextlib.redirect_stdout) can't catch
    them. In particular, TPORecoEvent's constructor (CoreUtils/
    TPORecoEvent.cc) unconditionally prints "Initializing Genfit" /
    "[GenFit] Material effects..." / "[GenFit] FieldManager initialized..."
    the first time one is built -- here, the first TPORecoEvent() in
    _iter_reco_entries() -- landing ahead of the JSON and breaking it as a
    single parseable document. Warnings this script prints to stderr
    elsewhere are unaffected (stderr is a different fd)."""
    sys.stdout.flush()
    saved_fd = os.dup(1)
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull_fd, 1)
        yield
    finally:
        sys.stdout.flush()
        os.dup2(saved_fd, 1)
        os.close(devnull_fd)
        os.close(saved_fd)


def mean(values):
    return sum(values) / len(values) if values else None


def rms(values):
    """Spread around the mean -- population standard deviation, ROOT's
    TH1::GetRMS() convention (sqrt(mean((x - mean(x))**2)), not
    sqrt(mean(x**2))). None when mean() would also be None."""
    if not values:
        return None
    m = mean(values)
    return (sum((v - m) ** 2 for v in values) / len(values)) ** 0.5


def _bind_po_event(ROOT, tree):
    po_event = ROOT.TPOEvent()
    tree.SetBranchAddress("event", ROOT.AddressOf(po_event) if hasattr(ROOT, "AddressOf") else po_event)
    return po_event


def _iter_truth_events(ROOT, truth_dir: Path, run: int, n_events: int):
    """Opens each FASERG4-Tcalevent_<run>_<event>.root under truth_dir (one
    per event, see CoreUtils/TcalEvent.cc), yielding each event's TPOEvent.
    Missing indices are skipped, not treated as an error."""
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
        po_event = _bind_po_event(ROOT, tree)
        # PyROOT can usually bind a branch of object type directly without
        # AddressOf(); the above covers both bindings depending on PyROOT
        # version. If neither works on your setup, see the fallback note
        # in docs/REGRESSION_TESTS.md and adjust this one line.
        tree.GetEntry(0)
        yield po_event
        tfile.Close()


def _truth_bucket_update(bucket, po_event):
    bucket["evis"].append(po_event.Evis)
    bucket["n_particles"].append(po_event.n_particles())
    # DIS kinematics are only meaningful for interaction-type events (Q2 >
    # 0 is the simplest available signal that kinematics_event() actually
    # ran for this event); skip them otherwise rather than padding the
    # mean with structural zeros.
    if po_event.Q2 > 0:
        bucket["nuE"].append(po_event.nuE)
        bucket["q2"].append(po_event.Q2)
        bucket["xbj"].append(po_event.xBj)
        bucket["yinel"].append(po_event.yInel)


def _truth_bucket_to_summary(bucket, n_found):
    summary = {
        "n_events": n_found,
        "mean_Evis": mean(bucket["evis"]),
        "rms_Evis": rms(bucket["evis"]),
        "mean_n_particles": mean(bucket["n_particles"]),
        "rms_n_particles": rms(bucket["n_particles"]),
    }
    if bucket["q2"]:
        summary["mean_nuE"] = mean(bucket["nuE"])
        summary["rms_nuE"] = rms(bucket["nuE"])
        summary["mean_Q2"] = mean(bucket["q2"])
        summary["rms_Q2"] = rms(bucket["q2"])
        summary["mean_xBj"] = mean(bucket["xbj"])
        summary["rms_xBj"] = rms(bucket["xbj"])
        summary["mean_yInel"] = mean(bucket["yinel"])
        summary["rms_yInel"] = rms(bucket["yinel"])
    return summary


def _new_truth_bucket():
    return {"evis": [], "n_particles": [], "nuE": [], "q2": [], "xbj": [], "yinel": []}


def summarize_truth(ROOT, truth_dir: Path, run: int, n_events: int):
    """Flat aggregate over every truth event found, regardless of its
    actual interaction type -- used by --muons/--muondis, and by the
    plain (non-split) mode."""
    bucket = _new_truth_bucket()
    found = 0
    for po_event in _iter_truth_events(ROOT, truth_dir, run, n_events):
        found += 1
        _truth_bucket_update(bucket, po_event)
    return _truth_bucket_to_summary(bucket, found)


def classify_truth(ROOT, truth_dir: Path, run: int, n_events: int):
    """Like summarize_truth(), but buckets each event by its own
    TPOEvent::reaction_desc() (e.g. "nueCC", "numuCC", "nutauCC", "nuNC")
    instead of lumping them all together -- see this script's module
    docstring. Returns {reaction: truth_summary_dict}."""
    buckets = {}
    counts = {}
    for po_event in _iter_truth_events(ROOT, truth_dir, run, n_events):
        reaction = str(po_event.reaction_desc())
        bucket = buckets.setdefault(reaction, _new_truth_bucket())
        _truth_bucket_update(bucket, po_event)
        counts[reaction] = counts.get(reaction, 0) + 1
    return {reaction: _truth_bucket_to_summary(bucket, counts[reaction])
            for reaction, bucket in buckets.items()}


def _new_reco_bucket():
    return {"n_porecs": [], "n_tktracks": [], "n_tkvertices": [], "n_mutracks": [],
            "total_evis": [], "total_ecompensated": []}


def _reco_bucket_update(bucket, reco):
    # fPORecs is private in TPORecoEvent (CoreUtils/TPORecoEvent.hh) --
    # unlike fTKTracks/fTKVertices/fMuTracks below, which are public and so
    # readable directly through PyROOT -- so it needs its getter instead.
    porecs = reco.GetPORecs()
    bucket["n_porecs"].append(len(porecs))
    bucket["n_tktracks"].append(len(reco.fTKTracks))
    bucket["n_tkvertices"].append(len(reco.fTKVertices))
    bucket["n_mutracks"].append(len(reco.fMuTracks))
    evis_sum = 0.0
    ecomp_sum = 0.0
    for porec in porecs:
        evis_sum += porec.TotalEvis()
        ecomp_sum += porec.fTotal.Ecompensated
    bucket["total_evis"].append(evis_sum)
    bucket["total_ecompensated"].append(ecomp_sum)


def _reco_bucket_to_summary(bucket, n_entries):
    return {
        "n_events_reconstructed": n_entries,
        "mean_n_PORecs": mean(bucket["n_porecs"]),
        "rms_n_PORecs": rms(bucket["n_porecs"]),
        "mean_n_TKTracks": mean(bucket["n_tktracks"]),
        "rms_n_TKTracks": rms(bucket["n_tktracks"]),
        "mean_n_TKVertices": mean(bucket["n_tkvertices"]),
        "rms_n_TKVertices": rms(bucket["n_tkvertices"]),
        "mean_n_MuTracks": mean(bucket["n_mutracks"]),
        "rms_n_MuTracks": rms(bucket["n_mutracks"]),
        "mean_total_Evis_reco": mean(bucket["total_evis"]),
        "rms_total_Evis_reco": rms(bucket["total_evis"]),
        "mean_total_Ecompensated": mean(bucket["total_ecompensated"]),
        "rms_total_Ecompensated": rms(bucket["total_ecompensated"]),
    }


def _iter_reco_entries(ROOT, reco_file: Path):
    """Reads one BatchReco output file (Batch-TPORecevent_<run>_<min>_
    <max>.root -- see Batch/BatchReco.cc), tree "RecoEvent", branch
    "TPORecoEvent" -- the exact pattern CoreUtils/dumpReco.C already uses
    (TChain + SetBranchAddress), just via PyROOT instead of a ROOT macro.
    Yields each entry's bound TPORecoEvent."""
    if not reco_file.is_file():
        sys.exit(f"error: reco file not found: {reco_file}")

    chain = ROOT.TChain("RecoEvent")
    chain.Add(str(reco_file))
    n_entries = chain.GetEntries()
    if n_entries == 0:
        return

    reco = ROOT.TPORecoEvent()
    chain.SetBranchAddress("TPORecoEvent", ROOT.AddressOf(reco) if hasattr(ROOT, "AddressOf") else reco)
    for ientry in range(n_entries):
        chain.GetEntry(ientry)
        yield reco


def summarize_reco(ROOT, reco_file: Path):
    """Flat aggregate over every reconstructed event in reco_file,
    regardless of its actual interaction type."""
    bucket = _new_reco_bucket()
    n_entries = 0
    for reco in _iter_reco_entries(ROOT, reco_file):
        n_entries += 1
        _reco_bucket_update(bucket, reco)
    if n_entries == 0:
        return {"n_events_reconstructed": 0}
    return _reco_bucket_to_summary(bucket, n_entries)


def classify_reco(ROOT, reco_file: Path):
    """Like summarize_reco(), but buckets each reconstructed event by its
    own truth TPOEvent's reaction_desc() -- BatchReco.cc constructs every
    TPORecoEvent as `new TPORecoEvent(fTcalEvent, fTcalEvent->fTPOEvent)`
    (see Batch/BatchReco.cc), and that per-event TPOEvent* is a persisted
    (non-"//!") member of TPORecoEvent (CoreUtils/TPORecoEvent.hh), so
    each reco entry carries its own truth event right along with it --
    accessible via TPORecoEvent::GetPOEvent(). No need to also cross
    reference the truth-dir files to classify a reco entry. Returns
    {reaction: reco_summary_dict}."""
    buckets = {}
    counts = {}
    for reco in _iter_reco_entries(ROOT, reco_file):
        po_event = reco.GetPOEvent()
        reaction = str(po_event.reaction_desc()) if po_event else "unknown"
        bucket = buckets.setdefault(reaction, _new_reco_bucket())
        _reco_bucket_update(bucket, reco)
        counts[reaction] = counts.get(reaction, 0) + 1
    return {reaction: _reco_bucket_to_summary(bucket, counts[reaction])
            for reaction, bucket in buckets.items()}


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
    parser.add_argument("--split-by-reaction", action="store_true",
                         help="Instead of one flat {truth, reco} summary, bucket events by each "
                              "event's own TPOEvent::reaction_desc() (nueCC/numuCC/nutauCC/nuNC/...) "
                              "and print {reaction: {truth, reco}, ...}. Use this for a single "
                              "unbiased/unmasked sample that mixes interaction types -- see this "
                              "script's module docstring and docs/REGRESSION_TESTS.md.")
    args = parser.parse_args()
    if args.truth_dir and (args.run is None or args.n_events is None):
        parser.error("--truth-dir requires --run and --n-events")
    if not args.truth_dir and not args.reco_file:
        parser.error("give at least one of --truth-dir or --reco-file")
    return args


def main():
    args = parse_args()

    with _silence_native_stdout():
        ROOT = load_dictionary(args.dict_path)
        ROOT.gROOT.SetBatch(True)

        if args.split_by_reaction:
            truth_by_reaction = (classify_truth(ROOT, args.truth_dir, args.run, args.n_events)
                                  if args.truth_dir else {})
            reco_by_reaction = classify_reco(ROOT, args.reco_file) if args.reco_file else {}
            reactions = set(truth_by_reaction) | set(reco_by_reaction)
            result = {}
            for reaction in reactions:
                entry = {}
                if reaction in truth_by_reaction:
                    entry["truth"] = truth_by_reaction[reaction]
                if reaction in reco_by_reaction:
                    entry["reco"] = reco_by_reaction[reaction]
                result[reaction] = entry
        else:
            result = {}
            if args.truth_dir:
                result["truth"] = summarize_truth(ROOT, args.truth_dir, args.run, args.n_events)
            if args.reco_file:
                result["reco"] = summarize_reco(ROOT, args.reco_file)

    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
