#!/usr/bin/env python3
"""Preflight check to run on the NEW muon-spectrometer-fixed ROOT files
BEFORE launching the full npz conversion.

Why this exists
---------------
TMuTrack is linked WITHOUT '+' in both the producer and the reader, which
uses FASER's CoreUtilsDict library built from CoreUtils/LinkDef.h. It uses the
legacy rootcint Streamer with NO schema evolution: members are streamed in
declaration order and the on-file bytes are interpreted using the COMPILED
class layout. Upstream has repeatedly changed TMuTrack's members WITHOUT
bumping ClassDef, so several mutually-incompatible layouts all call themselves
"version 2" and are distinguishable only by TClass checksum. Observed:

    1631166315  May-2026 layout, ends at fipErr
    3748019090  rubbiaa/FASER 5b39075, ClassDef 2, ends at ffit_ok + fpAnalytic
     663777236  rubbiaa/FASER a1bb833, ClassDef 3, adds `int fchargeMode`
                after fcharge (commit 2ee2913 "Fix charge identification in
                MDT")                                              <-- expected

If the compiled checksum does not match the on-file one, fMuTracks reads
produce silent garbage (wild momenta, vector::length_error), not an exception.
That is exactly the failure this script exists to catch, because catching it
after converting ~2M events is expensive.

Usage
-----
    source setup.sh
    python3 preflight_check_reco.py <one_new_RECO_file.root> [--caldata <FASERG4-Tcalevent_*.root>]

Exits 0 if everything matches, 1 otherwise.
"""

import argparse
import sys

from root_io import dictionary_path, load_dictionary as load_root_dictionary

ROOT = None

# Layout of rubbiaa/FASER a1bb833 (ClassDef 3 = ClassDef 2 + fchargeMode).
EXPECTED_TMUTRACK_CHECKSUM = 663777236

# Members the converter reads off each TMuTrack.
REQUIRED_TMUTRACK_MEMBERS = (
    "ftrackID", "fPDG", "fcharge", "fchargeMode", "fpos",
    "fpx", "fpy", "fpz", "fp", "fpErr",
    "fchi2", "fnDoF", "fpval", "fipErr",
    "fQOverP", "fQOverPErr", "ffit_ok", "fpAnalytic",
)


def load_dictionary():
    global ROOT
    ROOT = load_root_dictionary()
    ROOT.gErrorIgnoreLevel = ROOT.kError
    return str(dictionary_path())


def on_file_streamerinfo(path, clsname):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        sys.exit(f"FAIL: cannot open {path}")
    found = None
    lst = f.GetStreamerInfoList()
    if lst:
        for s in lst:
            if s.GetName() == clsname:
                found = (s.GetClassVersion(), s.GetCheckSum(),
                         [(e.GetTypeName(), e.GetName()) for e in s.GetElements()])
                break
    f.Close()
    return found


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("reco", help="one NEW Batch-TPORecevent_*.root file")
    ap.add_argument("--caldata", help="matching FASERG4-Tcalevent_*.root (optional)")
    ap.add_argument("--events", type=int, default=25, help="events to sample")
    args = ap.parse_args()

    lib = load_dictionary()
    print(f"dictionary: {lib}\n")

    failures = []

    # ---- 1. compiled vs on-file TMuTrack layout -----------------------------
    cl = ROOT.TClass.GetClass("TMuTrack")
    compiled_sum = cl.GetCheckSum()
    compiled_members = [m.GetName() for m in cl.GetListOfDataMembers()]
    print(f"compiled TMuTrack: version={cl.GetClassVersion()} checksum={compiled_sum}")

    if compiled_sum != EXPECTED_TMUTRACK_CHECKSUM:
        failures.append(
            f"compiled TMuTrack checksum {compiled_sum} != expected "
            f"{EXPECTED_TMUTRACK_CHECKSUM}. CoreUtils/TMuTrack.hh does not match "
            f"the expected production layout. Build the dictionary from the same "
            f"CoreUtils revision as the ROOT file producer."
        )

    missing = [m for m in REQUIRED_TMUTRACK_MEMBERS if m not in compiled_members]
    if missing:
        failures.append(f"compiled TMuTrack missing members used by the converter: {missing}")

    info = on_file_streamerinfo(args.reco, "TMuTrack")
    if info is None:
        print("NOTE: no TMuTrack StreamerInfo on file (no muon tracks written?)")
    else:
        version, checksum, elements = info
        print(f"on-file TMuTrack:  version={version} checksum={checksum}")
        if checksum != compiled_sum:
            failures.append(
                f"ON-FILE TMuTrack checksum {checksum} != compiled {compiled_sum}. "
                f"The legacy streamer has no schema evolution, so fMuTracks would be "
                f"read as SILENT GARBAGE. On-file members:\n      "
                + "\n      ".join(f"{t:<24} {n}" for t, n in elements)
            )

    # ---- 2. actually read some events ---------------------------------------
    print()
    f = ROOT.TFile.Open(args.reco)
    tree = f.Get("RecoEvent") or f.Get("Reco") or None
    if tree is None:
        keys = [k.GetName() for k in f.GetListOfKeys()]
        obj = None
        for k in keys:
            o = f.Get(k)
            if o and o.InheritsFrom("TTree"):
                tree, obj = o, k
                break
        if tree is not None:
            print(f"using tree '{obj}'")
    if tree is None:
        print("WARN: no TTree found; skipping value checks")
    else:
        ev = ROOT.TPORecoEvent()
        try:
            tree.SetBranchAddress("TPORecoEvent", ev)
        except Exception as exc:
            print(f"WARN: could not bind TPORecoEvent branch ({exc}); skipping value checks")
            tree = None

    if tree is not None:
        n = min(args.events, tree.GetEntries())
        ntracks = nfit_ok = nbad = ndup = 0
        pmax = 0.0
        for i in range(n):
            tree.GetEntry(i)
            seen = set()
            for t in ev.fMuTracks:
                ntracks += 1
                if t.ftrackID in seen:
                    ndup += 1
                seen.add(t.ftrackID)
                if t.ffit_ok:
                    nfit_ok += 1
                    p = t.fp
                    if not (0.0 < p < 1e5):
                        nbad += 1
                    pmax = max(pmax, p)
        print(f"sampled {n} events: {ntracks} fMuTracks, {nfit_ok} with ffit_ok, "
              f"{ndup} duplicate-trackID rows, max p={pmax:.1f} GeV")
        if nbad:
            failures.append(
                f"{nbad}/{nfit_ok} accepted fits have non-physical p -- "
                f"a classic symptom of a streamer layout mismatch."
            )
        if ndup:
            print("  NOTE: duplicate ftrackID rows are EXPECTED (upstream two-track "
                  "rescue). Deduplicate before counting muons.")
    f.Close()

    # ---- 3. CALDATA voxel response ------------------------------------------
    if args.caldata:
        print()
        vinfo = on_file_streamerinfo(args.caldata, "TcalEvent::FASERCALVOXELRESPONSE")
        if vinfo:
            names = [n for _, n in vinfo[2]]
            print(f"on-file FASERCALVOXELRESPONSE members: {names}")
            if "totalEnergyDeposit" in names:
                print("  totalEnergyDeposit present (new upstream field, '+' streamer "
                      "handles it via schema evolution)")
        mdt = on_file_streamerinfo(args.caldata, "MDTTrack")
        if mdt:
            mv, mc, _ = mdt
            print(f"on-file MDTTrack: version={mv} checksum={mc}")
            mcl = ROOT.TClass.GetClass("MDTTrack")
            if mcl and mcl.GetCheckSum() != mc:
                print(f"  NOTE: compiled MDTTrack checksum {mcl.GetCheckSum()} differs. "
                      f"MDTTrack IS linked with '+' so schema evolution applies, but "
                      f"check the members if you read mdttracks.")

    # ---- verdict -------------------------------------------------------------
    print()
    if failures:
        print("=" * 72)
        print("PREFLIGHT FAILED -- do NOT start the conversion")
        print("=" * 72)
        for i, msg in enumerate(failures, 1):
            print(f"  {i}. {msg}")
        return 1

    print("=" * 72)
    print("PREFLIGHT PASSED -- layouts match, safe to convert")
    print("=" * 72)
    return 0


if __name__ == "__main__":
    sys.exit(main())
