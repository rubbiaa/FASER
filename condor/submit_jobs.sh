#!/bin/bash
###############################################################################
# Build a chunk list and submit faserps / batchreco jobs to HTCondor.
# Runs ON lxplus, from inside the FASER checkout's condor/ directory (or
# from the Mac via submit_from_mac.sh, which just calls this over ssh).
#
#   submit_jobs.sh faserps   (--input-file FILE | --run N --detector D) \
#                            --total-events T [--chunk-size C] [--seed S] \
#                            [--faserdata DIR] [--dry-run] [-- run_faserps.py args]
#
#   submit_jobs.sh faserps   --total-events T [...] -- --muons --muon-momentum-gev 250
#
#   submit_jobs.sh batchreco --run N [--total-events T] [--chunk-size C] \
#                            [--faserdata DIR] [--dry-run] [-- run_batchreco.py args]
#
# --faserdata DIR   shared FASERDATA the jobs read inputs from and publish
#                   results to. Default: $FASER_SHARED_DATA, else $FASERDATA,
#                   else <checkout>/data. Put it on EOS, not in your AFS
#                   checkout: e.g. /eos/user/r/rubbiaa/FASERDATA. You can set
#                   FASER_SHARED_DATA once in ~/.faser_condor.conf on lxplus.
# --total-events    events to process, split into chunks of --chunk-size
#                   (default 5000). For batchreco it defaults to "everything
#                   faserps produced for --run" (highest event id + 1).
# --seed S          faserps base seed (default 123456789); chunk k uses S+k.
# --dry-run         create the job list and show the condor_submit command,
#                   but don't submit.
#
# Each submission gets a spool directory $FASER_CONDOR_SPOOL/<mode>_<timestamp>/
# (default ~/faser_condor, must be on AFS - CERN's standard schedds reject /eos
# paths in the submit file) holding the .sub copy, the job list and the logs.
###############################################################################
set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }
usage() { sed -n '2,/^####/p' "${BASH_SOURCE[0]}" | sed -e '1d' -e '$d' -e 's/^# \{0,1\}//'; }

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
HOMEFASER="$(cd "$HERE/.." && pwd)"

# Optional per-user settings on this machine (FASER_SHARED_DATA=...)
# shellcheck disable=SC1090
[[ -f "$HOME/.faser_condor.conf" ]] && source "$HOME/.faser_condor.conf"

[[ $# -ge 1 ]] || { usage; exit 1; }
MODE="$1"; shift
case "$MODE" in
  faserps|batchreco) ;;
  -h|--help|help) usage; exit 0 ;;
  *) usage; die "first argument must be 'faserps' or 'batchreco' (got '$MODE')" ;;
esac

SHARED="${FASER_SHARED_DATA:-${FASERDATA:-$HOMEFASER/data}}"
RUN=""; DETECTOR=""; INPUT_FILE=""; TOTAL=""; CHUNK_SIZE=5000; SEED0=123456789; DRY=0
EXTRA=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --faserdata)    SHARED="$2"; shift 2 ;;
    --run)          RUN="$2"; shift 2 ;;
    --detector)     DETECTOR="$2"; shift 2 ;;
    --input-file)   INPUT_FILE="$2"; shift 2 ;;
    --total-events) TOTAL="$2"; shift 2 ;;
    --chunk-size)   CHUNK_SIZE="$2"; shift 2 ;;
    --seed)         SEED0="$2"; shift 2 ;;
    --dry-run)      DRY=1; shift ;;
    -h|--help)      usage; exit 0 ;;
    --)             shift; EXTRA=("$@"); break ;;
    *)              die "unknown option '$1' (arguments for run_*.py go after a literal --)" ;;
  esac
done

num() { [[ "$2" =~ ^[0-9]+$ ]] || die "$1 must be a non-negative integer (got '$2')"; }
[[ -d "$SHARED" ]] || die "shared FASERDATA does not exist: $SHARED (create it, or pass --faserdata)"
num --chunk-size "$CHUNK_SIZE"; num --seed "$SEED0"
[[ "$CHUNK_SIZE" -ge 1 ]] || die "--chunk-size must be >= 1"

[[ -x "$HOMEFASER/build/bin/$([[ $MODE == faserps ]] && echo faserps || echo batchreco.exe)" ]] \
  || die "$HOMEFASER/build/bin/ has no $MODE binary - build on lxplus first (cd $HOMEFASER; source setup.sh; fb)"

INPUT="none"
if [[ "$MODE" == "faserps" ]]; then
  [[ -n "$TOTAL" ]] || die "faserps needs --total-events"
  num --total-events "$TOTAL"

  muon_mode=0
  for a in "${EXTRA[@]+"${EXTRA[@]}"}"; do
    [[ "$a" == "--muons" || "$a" == "--muondis" ]] && muon_mode=1
  done

  if [[ "$muon_mode" -eq 1 ]]; then
    INPUT="none"
  elif [[ -n "$INPUT_FILE" ]]; then
    [[ -f "$INPUT_FILE" ]] || die "--input-file not found: $INPUT_FILE"
    INPUT="$(readlink -f "$INPUT_FILE")"
  else
    [[ -n "$RUN" && -n "$DETECTOR" ]] || die "faserps needs --input-file FILE, or --run N --detector {3DCAL,AHCAL,ECAL} (looked up under \$FASERDATA/CVGENIE), or -- --muons"
    num --run "$RUN"
    shopt -s nullglob
    cands=( "$SHARED/CVGENIE/Run$RUN"/FASERMC-PO-Run"$RUN"-*_"$DETECTOR".root
            "$SHARED/CVGENIE/Run$RUN"/FASERMC-PO-Run"$RUN"-*_"$DETECTOR"_*.root )
    shopt -u nullglob
    [[ ${#cands[@]} -ge 1 ]] || die "no FASERMC-PO-Run$RUN-*_$DETECTOR*.root under $SHARED/CVGENIE/Run$RUN (copy the converted CVGENIE sample there, or run run_convertgenie.py on lxplus)"
    INPUT="$(readlink -f "${cands[0]}")"
  fi
else
  [[ -n "$RUN" ]] || die "batchreco needs --run N"
  num --run "$RUN"
  if [[ -z "$TOTAL" ]]; then
    TOTAL_ID="$(find "$SHARED/faserG4" -maxdepth 1 -name "FASERG4-Tcalevent_${RUN}_*.root" -printf '%f\n' 2>/dev/null \
                | sed -E "s/^FASERG4-Tcalevent_${RUN}_([0-9]+)\.root$/\1/" | sort -n | tail -1)"
    [[ -n "$TOTAL_ID" ]] || die "no FASERG4-Tcalevent_${RUN}_*.root under $SHARED/faserG4 - nothing to reconstruct (or pass --total-events)"
    TOTAL=$((TOTAL_ID + 1))
    echo "batchreco: highest event id for run $RUN is $TOTAL_ID -> --total-events $TOTAL"
  fi
  num --total-events "$TOTAL"
  [[ -s "$SHARED/GDML/FASERCAL_V10.gdml" ]] || die "no $SHARED/GDML/FASERCAL_V10.gdml - run faserps first (it publishes it) or copy it there"
fi
[[ "$TOTAL" -ge 1 ]] || die "--total-events must be >= 1"

# CERN's standard batch schedds refuse a submit file whose paths are on /eos
# (the executable, logs and job list must be on AFS), and the checkout may well
# live on EOS. So each submission gets its own small spool directory on AFS with
# a copy of the .sub file and the wrapper, the job list and the logs. The jobs
# themselves still read the checkout and the shared data on /eos from the worker.
SPOOL_BASE="${FASER_CONDOR_SPOOL:-$HOME/faser_condor}"
mkdir -p "$SPOOL_BASE" || die "cannot create spool directory $SPOOL_BASE (set FASER_CONDOR_SPOOL to a directory on AFS)"
case "$(readlink -f "$SPOOL_BASE")" in
  /eos/*) die "spool directory $SPOOL_BASE is on /eos; condor_submit needs it on AFS (set FASER_CONDOR_SPOOL to a directory on AFS, e.g. your ~/faser_condor)" ;;
esac
RUNDIR="$SPOOL_BASE/${MODE}_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RUNDIR/logs"
cp "$HERE/submit_${MODE}.sub" "$HERE/${MODE}_chunk.sh" "$RUNDIR/"
chmod +x "$RUNDIR/${MODE}_chunk.sh"
JOBS="$RUNDIR/jobs.list"
python3 "$HERE/make_jobs_list.py" --total-events "$TOTAL" --chunk-size "$CHUNK_SIZE" --out "$JOBS"

args=( -append "homefaser = $HOMEFASER"
       -append "shared = $SHARED"
       -append "jobslist = jobs.list"
       -append "extra = ${EXTRA[*]:-}" )
if [[ "$MODE" == "faserps" ]]; then
  args+=( -append "inputroot = $INPUT" -append "seed0 = $SEED0" )
else
  args+=( -append "run = $RUN" )
fi

echo "mode:        $MODE"
echo "checkout:    $HOMEFASER"
echo "shared data: $SHARED"
[[ "$MODE" == "faserps" ]] && echo "input:       $INPUT"
echo "spool dir:   $RUNDIR"
echo "job list:    $JOBS ($(wc -l < "$JOBS" | tr -d ' ') chunks)"
echo "command:     (cd $RUNDIR && condor_submit ${args[*]} submit_${MODE}.sub)"

if [[ "$DRY" -eq 1 ]]; then
  echo "--dry-run: not submitting."
  exit 0
fi

command -v condor_submit >/dev/null || die "condor_submit not found - run this on lxplus (or use submit_from_mac.sh)"
cd "$RUNDIR"
condor_submit "${args[@]}" "submit_${MODE}.sub"
echo "Submitted. Monitor with: condor_q   (logs in $RUNDIR/logs)"
