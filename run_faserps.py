#!/usr/bin/env python3
"""
run_faserps.py -- run FASERG4's `faserps` simulation with a macro built
as an in-memory Python string, piped straight to the executable's stdin.
No .mac file is ever written to disk.

Why: FASERG4/ used to carry ~17 run*.mac files (RunFASER_newgeo.mac,
runFASER_1.mac, runFASER_muon.mac, runFASER_V5.mac, ...), almost all
leftovers from earlier detector geometries -- and nothing stopped
faserps from being run with a stale one, silently simulating the wrong
detector. This script replaces all of them: the one macro that matches
the current V10 geometry lives here as build_v10_macro(), a parametrized
function, and its text is piped to `faserps -` (a new mode added to
faserps.cc that reads G4 UI commands from stdin instead of
/control/execute'ing a file). A future "run variant" is a new function
(or new keyword arguments to this one), not a new file to keep in sync
by hand. This includes MuonDIS (--muondis, below), ported the same way
from the standalone runFASER_muondis.mac example that used to live under
FASERG4/.

Usage:
    python3 run_faserps.py                        # defaults matching the old runFASER_V10.mac
    python3 run_faserps.py --n-events 500 --start-event 100
    python3 run_faserps.py --input-file some_other_sample.root
    python3 run_faserps.py --muons --n-events 1000 # single 100 GeV muons instead of neutrino events
    python3 run_faserps.py --muons --muon-momentum-gev 250
    python3 run_faserps.py --muondis --n-events 1000  # muons with MuonDIS (Pythia8 DIS) enabled
    python3 run_faserps.py --print-macro           # print the macro, don't run it
    python3 run_faserps.py --vis                   # interactive G4 UI (unrelated to the macro)
    python3 run_faserps.py --build-dir PATH        # override the build/ location
"""
import argparse
import os
import subprocess
import sys
from pathlib import Path

# This script lives at the repo root (moved out of FASERG4/ so the
# runner and its doc, run_faserps.md, stay together) - FASERG4_DIR is
# therefore an explicit subdirectory, not Path(__file__).parent.
REPO_ROOT = Path(__file__).resolve().parent
FASERG4_DIR = REPO_ROOT / "FASERG4"
DEFAULT_BUILD_DIR = REPO_ROOT / "build"

# FASERMC-PO-Run10000-0_53954_3DCAL.root (the default --input-file sample)
# lives under $FASERDATA/GENIE/, not next to faserps in FASERG4/ anymore --
# it's an input data file, not a source file, so it belongs alongside
# faserps' own faserG4/batch output rather than checked out with the code.
# This mirrors common_setup.sh's own FASERDATA fallback
# (`${FASERDATA:=$HOMEFASER/data}`) so the default resolves correctly
# whether or not setup.sh has been sourced yet -- unlike the faserG4/batch
# output dirs, we can't rely on FASER::GetDataDir() (C++-only) to create
# this one, since it's an input we ship, not output faserps generates.
def _default_genie_input_file():
    faserdata = Path(os.environ.get("FASERDATA", str(REPO_ROOT / "data")))
    return str(faserdata / "GENIE" / "FASERMC-PO-Run10000-0_53954_3DCAL.root")


def build_v10_macro(
    *,
    scint_size_x_cm: float = 48.0,
    scint_size_y_cm: float = 48.0,
    scint_size_z_cm: float = 20.0,
    scint_voxel_mm: float = 10.0,
    targetw_size_x_cm: float = 48.0,
    targetw_size_y_cm: float = 48.0,
    targetw_size_z_mm: float = 5.0,
    n_layers: int = 10,
    los_shift_x_cm: float = 45.0,
    los_shift_y_cm: float = 24.0,
    tilt_deg: float = -4.5,
    tracking_verbose: int = 0,
    run_verbose: int = 0,
    control_verbose: int = 0,
    n_threads: int = 1,
    input_root_file: str = None,  # resolved below -- see _default_genie_input_file()
    start_event: int = 0,
    n_events: int = 100,
    muon_mode: bool = False,
    muon_momentum_gev: float = 100.0,
    muondis_mode: bool = False,
    muondis_cross_section_bias: float = 150.0,
    muondis_q2min: float = 1.0,
    muondis_interaction_log: str = "",
    muondis_pdf_set: str = "",
    muondis_xbjmin: float = 0.0,
    muondis_debug: bool = False,
) -> str:
    """Builds the V10-geometry macro as a string. Defaults match the old
    runFASER_V10.mac exactly; pass keyword overrides for a variant instead
    of hand-editing a copy of a .mac file.

    muon_mode=True switches the generator from reading neutrino-interaction
    events out of `input_root_file` to firing single fixed-momentum muons
    instead (PrimaryGeneratorAction::GeneratePrimaries only opens the ROOT
    input file when both wantMuonBackground and wantSingleParticle are
    false - see FASERG4/src/PrimaryGeneratorAction.cc - so
    input_root_file/start_event are simply unused in this mode, and left
    out of the generated macro rather than printed misleadingly).

    muondis_mode=True additionally enables MuonDIS (see
    docs/README_MuonDIS.md): the primary muon's nuclear interaction is
    replaced by an on-the-fly Pythia8 deep-inelastic-scattering event
    instead of Geant4's standard muon-nuclear final state. This only makes
    sense for muon-background primaries, so callers should also set
    muon_mode=True (run_faserps.py's --muondis flag does this for you via
    main()). The /physics/muondis/... commands must be issued before
    /run/initialize, so they're emitted just above it. crossSectionBias and
    q2min default to the same values MuonDISPhysics.hh itself defaults to
    (150.0, 1.0) so the macro stays self-documenting even when nothing is
    overridden; interactionLog/pdfSet/xbjmin/debug default to "off" (empty
    string / 0.0 / false) and are only emitted when explicitly given, since
    unlike crossSectionBias/q2min they have no meaningful non-empty
    default on the C++ side either."""
    if input_root_file is None:
        input_root_file = _default_genie_input_file()
    if muon_mode:
        generator_lines = "\n".join([
            "/generator/wantMuonBackground true",
            f"/generator/singleMomentum {muon_momentum_gev:g} GeV",
        ])
    else:
        generator_lines = "\n".join([
            f"/generator/rootinputfilename {input_root_file}",
            f"/generator/startevent {start_event}",
        ])

    if muondis_mode:
        muondis_cmds = [
            "/physics/muondis/enable true",
            f"/physics/muondis/crossSectionBias {muondis_cross_section_bias:g}",
            f"/physics/muondis/q2min {muondis_q2min:g}",
        ]
        if muondis_interaction_log:
            muondis_cmds.append(f"/physics/muondis/interactionLog {muondis_interaction_log}")
        if muondis_pdf_set:
            muondis_cmds.append(f"/physics/muondis/pdfSet {muondis_pdf_set}")
        if muondis_xbjmin:
            muondis_cmds.append(f"/physics/muondis/xbjmin {muondis_xbjmin:g}")
        if muondis_debug:
            muondis_cmds.append("/physics/muondis/debug true")
        muondis_block = (
            "#\n"
            "# MuonDIS: replace the primary muon's nuclear interaction with an\n"
            "# on-the-fly Pythia8 DIS event (see docs/README_MuonDIS.md).\n"
            "# Must be set before /run/initialize.\n"
            + "\n".join(muondis_cmds) + "\n"
        )
    else:
        muondis_block = ""

    return f"""# V10 geometry, tilted -- generated by run_faserps.build_v10_macro(), not a .mac file
/FASER/scint/sizeX {scint_size_x_cm:g} cm
/FASER/scint/sizeY {scint_size_y_cm:g} cm
/FASER/scint/sizeZ {scint_size_z_cm:g} cm
/FASER/scint/voxel {scint_voxel_mm:g} mm
/FASER/targetW/sizeX {targetw_size_x_cm:g} cm
/FASER/targetW/sizeY {targetw_size_y_cm:g} cm
/FASER/targetW/sizeZ {targetw_size_z_mm:g} mm
/FASER/layers {n_layers}
/FASER/LOS/shiftX {los_shift_x_cm:g} cm
/FASER/LOS/shiftY {los_shift_y_cm:g} cm
/FASER/tiltY {tilt_deg:g} deg
#
/tracking/verbose {tracking_verbose}
/run/verbose {run_verbose}
/control/verbose {control_verbose}
/run/numberOfThreads {n_threads}
{muondis_block}/run/initialize
{generator_lines}
/run/beamOn {n_events}
"""


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run FASERG4's faserps with the V10 macro, built as a "
                     "Python string and piped to stdin -- no .mac file involved.",
    )
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD_DIR,
                         help=f"CMake build directory containing bin/faserps (default: {DEFAULT_BUILD_DIR})")
    parser.add_argument("--vis", action="store_true",
                         help="Launch faserps' interactive Geant4 UI (`faserps vis`) instead of "
                              "piping a macro. No macro is used in this mode.")
    parser.add_argument("--print-macro", action="store_true",
                         help="Print the resolved macro text and exit without running faserps.")
    parser.add_argument("--dry-run", action="store_true",
                         help="Print the command (and the macro) but don't execute it.")

    # The handful of parameters actually worth varying run-to-run; the
    # rest of build_v10_macro()'s keyword arguments cover the geometry
    # constants and are left at their V10 defaults unless you call the
    # function yourself from a custom script.
    parser.add_argument("--input-file", default=_default_genie_input_file(),
                         help="Value for /generator/rootinputfilename. Defaults to "
                              "$FASERDATA/GENIE/FASERMC-PO-Run10000-0_53954_3DCAL.root (an absolute "
                              "path, so it resolves regardless of cwd). A custom value that isn't "
                              "already absolute is still resolved relative to FASERG4/, since that's "
                              "faserps' cwd when run_faserps.py invokes it. "
                              "Ignored if --muons or --muondis is given.")
    parser.add_argument("--start-event", type=int, default=0,
                         help="Value for /generator/startevent. Ignored if --muons or --muondis is given.")
    parser.add_argument("--n-events", type=int, default=100,
                         help="Value for /run/beamOn.")
    parser.add_argument("--muons", action="store_true",
                         help="Generate single fixed-momentum muons instead of reading "
                              "neutrino-interaction events from --input-file (adds "
                              "/generator/wantMuonBackground true and /generator/singleMomentum).")
    parser.add_argument("--muon-momentum-gev", type=float, default=100.0,
                         help="Muon momentum in GeV, only used with --muons/--muondis (default: 100).")
    parser.add_argument("--muondis", action="store_true",
                         help="Enable MuonDIS: replace the primary muon's nuclear interaction with "
                              "an on-the-fly Pythia8 deep-inelastic-scattering event, generated per "
                              "interaction from the muon's actual Geant4 energy/direction. Implies "
                              "--muons (MuonDIS only applies to muon-background primaries). See "
                              "docs/README_MuonDIS.md for the physics and every /physics/muondis/... "
                              "option.")
    parser.add_argument("--muondis-cross-section-bias", type=float, default=150.0,
                         help="Value for /physics/muondis/crossSectionBias, only used with --muondis "
                              "(default: 150, matching MuonDISPhysics.hh's own default -- biases the "
                              "interaction cross section up so DIS events are frequent enough to study "
                              "without huge statistics; use 1 for an unbiased cross section).")
    parser.add_argument("--muondis-q2min", type=float, default=1.0,
                         help="Value for /physics/muondis/q2min (GeV^2), only used with --muondis "
                              "(default: 1.0, matching MuonDISPhysics.hh's own default).")
    parser.add_argument("--muondis-interaction-log", default="",
                         help="Value for /physics/muondis/interactionLog, only used with --muondis. "
                              "Leave unset to disable the CSV interaction log (default).")
    parser.add_argument("--muondis-pdf-set", default="",
                         help="Value for /physics/muondis/pdfSet (path under FASERG4/input/, e.g. "
                              "input/NNPDF40_nnlo_as_01180_charmasy_0000.dat), only used with "
                              "--muondis. Leave unset to use Pythia8's built-in proton PDF (default).")
    parser.add_argument("--muondis-xbjmin", type=float, default=0.0,
                         help="Value for /physics/muondis/xbjmin, only used with --muondis. "
                              "Leave at 0 to disable the cut (default).")
    parser.add_argument("--muondis-debug", action="store_true",
                         help="Add /physics/muondis/debug true, only used with --muondis.")
    parser.add_argument("--tilt-deg", type=float, default=-4.5,
                         help="Value for /FASER/tiltY (degrees).")
    parser.add_argument("--shift-x-cm", type=float, default=45.0,
                         help="Value for /FASER/LOS/shiftX (cm).")
    parser.add_argument("--shift-y-cm", type=float, default=24.0,
                         help="Value for /FASER/LOS/shiftY (cm).")
    return parser.parse_args()


def main():
    args = parse_args()

    macro_text = build_v10_macro(
        input_root_file=args.input_file,
        start_event=args.start_event,
        n_events=args.n_events,
        tilt_deg=args.tilt_deg,
        los_shift_x_cm=args.shift_x_cm,
        los_shift_y_cm=args.shift_y_cm,
        muon_mode=args.muons or args.muondis,
        muon_momentum_gev=args.muon_momentum_gev,
        muondis_mode=args.muondis,
        muondis_cross_section_bias=args.muondis_cross_section_bias,
        muondis_q2min=args.muondis_q2min,
        muondis_interaction_log=args.muondis_interaction_log,
        muondis_pdf_set=args.muondis_pdf_set,
        muondis_xbjmin=args.muondis_xbjmin,
        muondis_debug=args.muondis_debug,
    )

    if args.print_macro:
        print(macro_text, end="")
        return 0

    binary_path = args.build_dir / "bin" / "faserps"
    if not binary_path.is_file():
        sys.exit(
            f"error: faserps binary not found: {binary_path}\n"
            f"       (build it first, e.g. `cmake --build {args.build_dir} --target faserps`)"
        )
    if not (binary_path.stat().st_mode & 0o111):
        sys.exit(f"error: {binary_path} is not executable")

    # faserps writes its TcalEvent ROOT output under $FASERDATA/faserG4
    # (see CoreUtils/FaserDataDir.hh / CoreUtils/TcalEvent.cc), not a
    # relative "output/" path -- so FASERDATA must be set in the
    # environment this subprocess inherits (source setup.sh first).
    # The C++ side creates $FASERDATA/faserG4 itself if it's missing, but
    # checking FASERDATA here too gives a clear error up front instead of
    # a ROOT/C++ exception partway through the run.
    if not os.environ.get("FASERDATA"):
        sys.exit(
            "error: FASERDATA is not set.\n"
            "       Source setup.sh first (`source setup.sh`), or export FASERDATA\n"
            "       yourself to point at FASER's consolidated data directory."
        )

    command = [str(binary_path), "vis" if args.vis else "-"]

    print(f"[run_faserps] working directory: {FASERG4_DIR}")
    print(f"[run_faserps] command:           {' '.join(command)}")
    if not args.vis:
        print("[run_faserps] macro (piped via stdin, not written to disk):")
        print("  " + macro_text.strip().replace("\n", "\n  "))

    if args.dry_run:
        print("[run_faserps] --dry-run: not executing.")
        return 0

    # cwd=FASERG4_DIR because the macro's /generator/rootinputfilename
    # is a relative path that only resolves correctly from here (same as
    # the GDML faserps writes out). Output no longer depends on cwd at
    # all -- see the FASERDATA check above.
    if args.vis:
        result = subprocess.run(command, cwd=FASERG4_DIR)
    else:
        result = subprocess.run(command, cwd=FASERG4_DIR, input=macro_text.encode())
    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
