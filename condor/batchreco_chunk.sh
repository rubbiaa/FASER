#!/bin/bash
###############################################################################
# HTCondor worker wrapper: reconstruct one event range with the FASER repo's
# own setup.sh and run_batchreco.py (as with faserps_chunk.sh, only this
# script is transferred - the build and the data are read from the checkout
# and the shared FASERDATA area).
#
# usage:
#   batchreco_chunk.sh <homefaser> <shared_faserdata> <run> <start_evt> <nevt> <chunk> \
#                      [run_batchreco.py args...]
#
# Reconstructs events [start_evt, start_evt + nevt) of <run> (batchreco's
# --max-event is exclusive, matching this), reading
#   <shared>/faserG4/FASERG4-Tcalevent_<run>_<evt>.root   (written by faserps_chunk.sh)
#   <shared>/GDML/FASERCAL_V10.gdml                       (published by faserps_chunk.sh)
# and publishing
#   <shared>/batch/Batch-TPORecevent_<run>_<start>_<start+nevt>.root
#
# Inputs are staged to job-local scratch first (as the old xrdcp-based
# script did), batchreco runs against a private job-local FASERDATA, and the
# result is only copied to the shared area once it has been verified.
###############################################################################

set -eo pipefail

usage() {
  echo "Usage: $0 <homefaser> <shared_faserdata> <run> <start_evt> <nevt> <chunk> [run_batchreco.py args...]" >&2
}

if [[ $# -lt 6 ]]; then usage; exit 1; fi

HOMEFASER="$1"; SHARED="$2"; RUN="$3"; START="$4"; NEVT="$5"; CHUNK="$6"
shift 6
EXTRA=("$@")

for name in RUN START NEVT CHUNK; do
  [[ "${!name}" =~ ^[0-9]+$ ]] || { echo "ERROR: $name must be a non-negative integer (got '${!name}')" >&2; exit 2; }
done
if [[ "$NEVT" -lt 1 ]]; then echo "ERROR: nevt must be >= 1" >&2; exit 2; fi
MIN=$START
MAX=$((START + NEVT))   # exclusive

echo "=== batchreco chunk start ==="
echo "Date:      $(date)"
echo "Host:      $(hostname)"
echo "HOMEFASER: $HOMEFASER"
echo "SHARED:    $SHARED"
echo "run=$RUN chunk=$CHUNK events=[$MIN,$MAX) (n=$NEVT)"
echo "extra:     ${EXTRA[*]:-(none)}"

[[ -f "$HOMEFASER/setup.sh" ]]               || { echo "ERROR: no setup.sh in HOMEFASER=$HOMEFASER" >&2; exit 3; }
[[ -x "$HOMEFASER/build/bin/batchreco.exe" ]] || { echo "ERROR: $HOMEFASER/build/bin/batchreco.exe missing - build on lxplus first (source setup.sh; fb)" >&2; exit 4; }
[[ -f "$HOMEFASER/run_batchreco.py" ]]       || { echo "ERROR: no run_batchreco.py in HOMEFASER=$HOMEFASER" >&2; exit 3; }
[[ -d "$SHARED/faserG4" ]]                   || { echo "ERROR: no $SHARED/faserG4 - nothing to reconstruct" >&2; exit 5; }
[[ -s "$SHARED/GDML/FASERCAL_V10.gdml" ]]    || { echo "ERROR: no $SHARED/GDML/FASERCAL_V10.gdml (faserps publishes it; copy it there by hand otherwise)" >&2; exit 22; }

###############################################################################
# Job-local scratch + FASERDATA
###############################################################################
SCRATCH="${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}"
WORKDIR="$SCRATCH/batchreco_run${RUN}_chunk${CHUNK}_$$"
mkdir -p "$WORKDIR"
cd "$WORKDIR"
echo "WORKDIR=$WORKDIR"

cleanup() {
  local rc=$?
  cd / 2>/dev/null || true
  if [[ "${KEEP_WORKDIR:-0}" != "1" ]]; then rm -rf "$WORKDIR"; else echo "KEEP_WORKDIR=1: leaving $WORKDIR"; fi
  exit "$rc"
}
trap cleanup EXIT

# See faserps_chunk.sh: reuse the login node's Pythia8 shim if visible,
# otherwise make a private one in scratch.
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${HOME:-/nonexistent}/.cache}"
if ! compgen -G "$XDG_CACHE_HOME/faser/pythia8-shim/libpythia8-*.so" > /dev/null; then
  export XDG_CACHE_HOME="$WORKDIR/cache"
fi

export FASERDATA="$WORKDIR/data"            # before setup.sh, see faserps_chunk.sh
export HEARTBEAT_FILE="$WORKDIR/heartbeat.txt"  # batchreco's liveness probe, default is ./heartbeat.txt

###############################################################################
# FASER environment
###############################################################################
set +e +o pipefail
source "$HOMEFASER/setup.sh"
setup_rc=$?
set -eo pipefail
[[ $setup_rc -eq 0 ]]              || { echo "ERROR: setup.sh failed (rc=$setup_rc)" >&2; exit 7; }
command -v root-config >/dev/null  || { echo "ERROR: ROOT not on PATH after setup.sh" >&2; exit 7; }
[[ "$FASERDATA" == "$WORKDIR/data" ]] || { echo "ERROR: setup.sh changed FASERDATA to '$FASERDATA'" >&2; exit 7; }

###############################################################################
# Stage inputs (geometry + the event files of this chunk)
###############################################################################
mkdir -p "$FASERDATA/faserG4" "$FASERDATA/GDML" "$FASERDATA/batch"
cp "$SHARED/GDML/FASERCAL_V10.gdml" "$FASERDATA/GDML/FASERCAL_V10.gdml"

echo "staging $NEVT event file(s) from $SHARED/faserG4 ..."
missing=0
for ((evt = MIN; evt < MAX; evt++)); do
  f="FASERG4-Tcalevent_${RUN}_${evt}.root"
  if [[ -s "$SHARED/faserG4/$f" ]]; then
    cp "$SHARED/faserG4/$f" "$FASERDATA/faserG4/$f"
  else
    missing=$((missing + 1))
    [[ "$missing" -le 10 ]] && echo "  missing input: $SHARED/faserG4/$f" >&2
  fi
done
if [[ "$missing" -ne 0 ]]; then
  echo "ERROR: $missing of $NEVT input event file(s) missing - run faserps for this range first" >&2
  exit 21
fi
n_in=$(find "$FASERDATA/faserG4" -maxdepth 1 -type f -name "FASERG4-Tcalevent_${RUN}_*.root" | wc -l | tr -d ' ')
[[ "$n_in" -eq "$NEVT" ]] || { echo "ERROR: staged $n_in input file(s), expected $NEVT" >&2; exit 21; }

###############################################################################
# Run
###############################################################################
cmd=(python3 "$HOMEFASER/run_batchreco.py" --run "$RUN" --min-event "$MIN" --max-event "$MAX")
cmd+=("${EXTRA[@]}")

echo "=== running: ${cmd[*]}"
"${cmd[@]}"
echo "=== run_batchreco.py finished ==="

###############################################################################
# Verify, then publish
###############################################################################
out="Batch-TPORecevent_${RUN}_${MIN}_${MAX}.root"
if [[ ! -s "$FASERDATA/batch/$out" ]]; then
  echo "ERROR: expected output $FASERDATA/batch/$out is missing or empty" >&2
  exit 30
fi

mkdir -p "$SHARED/batch"
cp -f "$FASERDATA/batch/$out" "$SHARED/batch/$out"
if [[ "$(stat -c %s "$FASERDATA/batch/$out")" != "$(stat -c %s "$SHARED/batch/$out" 2>/dev/null || echo -1)" ]]; then
  echo "ERROR: copy of $out to $SHARED/batch is incomplete" >&2
  exit 31
fi

echo "=== published $SHARED/batch/$out ==="
echo "=== batchreco chunk done: $(date) ==="
