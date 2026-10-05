#!/bin/bash
###############################################################################
# HTCondor worker wrapper: one chunk of a faserps simulation, using the
# FASER repo's own setup.sh and run_faserps.py (nothing is transferred to
# the worker except this script - the build, setup.sh and the input sample
# are read straight from the checkout / shared data area, which lxplus
# batch nodes see through /afs, /eos and /cvmfs).
#
# usage:
#   faserps_chunk.sh <homefaser> <shared_faserdata> <inputroot|none> \
#                    <start_evt> <nevt> <chunk> <seed0> [run_faserps.py args...]
#
#   homefaser         the FASER checkout on lxplus (contains setup.sh, build/)
#   shared_faserdata  FASERDATA the results are published into
#                     (<shared>/faserG4/FASERG4-Tcalevent_<run>_<evt>.root)
#   inputroot         absolute path of the converted PO input file
#                     ($FASERDATA/CVGENIE/Run<N>/FASERMC-PO-...root), or
#                     "none" for --muons/--muondis runs
#   seed0             base random seed; this chunk uses seed0 + chunk
#
# The job itself runs with a private, job-local FASERDATA in the scratch
# directory, so concurrent jobs never write the same GDML/file, and only
# verified results are copied to the shared area at the end.
###############################################################################

# No -u (setup.sh is not nounset-clean); -e/pipefail are switched off while
# sourcing it (see below) and back on afterwards.
set -eo pipefail

usage() {
  echo "Usage: $0 <homefaser> <shared_faserdata> <inputroot|none> <start_evt> <nevt> <chunk> <seed0> [run_faserps.py args...]" >&2
}

if [[ $# -lt 7 ]]; then usage; exit 1; fi

HOMEFASER="$1"; SHARED="$2"; INPUTROOT="$3"; START="$4"; NEVT="$5"; CHUNK="$6"; SEED0="$7"
shift 7
EXTRA=("$@")

for name in START NEVT CHUNK SEED0; do
  [[ "${!name}" =~ ^[0-9]+$ ]] || { echo "ERROR: $name must be a non-negative integer (got '${!name}')" >&2; exit 2; }
done
SEED=$((SEED0 + CHUNK))

echo "=== faserps chunk start ==="
echo "Date:      $(date)"
echo "Host:      $(hostname)"
echo "HOMEFASER: $HOMEFASER"
echo "SHARED:    $SHARED"
echo "chunk=$CHUNK start=$START nevt=$NEVT seed=$SEED (seed0=$SEED0)"
echo "input:     $INPUTROOT"
echo "extra:     ${EXTRA[*]:-(none)}"

[[ -f "$HOMEFASER/setup.sh" ]]            || { echo "ERROR: no setup.sh in HOMEFASER=$HOMEFASER" >&2; exit 3; }
[[ -x "$HOMEFASER/build/bin/faserps" ]]   || { echo "ERROR: $HOMEFASER/build/bin/faserps missing - build on lxplus first (source setup.sh; fb)" >&2; exit 4; }
[[ -f "$HOMEFASER/run_faserps.py" ]]      || { echo "ERROR: no run_faserps.py in HOMEFASER=$HOMEFASER" >&2; exit 3; }
[[ -d "$SHARED" ]]                        || { echo "ERROR: shared FASERDATA does not exist: $SHARED" >&2; exit 5; }
if [[ "$INPUTROOT" != "none" ]]; then
  [[ -f "$INPUTROOT" ]]                   || { echo "ERROR: input file not found: $INPUTROOT" >&2; exit 6; }
fi

###############################################################################
# Job-local scratch + FASERDATA
###############################################################################
SCRATCH="${_CONDOR_SCRATCH_DIR:-${TMPDIR:-/tmp}}"
WORKDIR="$SCRATCH/faserps_chunk${CHUNK}_$$"
mkdir -p "$WORKDIR"
cd "$WORKDIR"
echo "WORKDIR=$WORKDIR"

# Clean up on exit (success or failure). KEEP_WORKDIR=1 to debug.
cleanup() {
  local rc=$?
  cd / 2>/dev/null || true
  if [[ "${KEEP_WORKDIR:-0}" != "1" ]]; then rm -rf "$WORKDIR"; else echo "KEEP_WORKDIR=1: leaving $WORKDIR"; fi
  exit "$rc"
}
trap cleanup EXIT

# setup.sh (lxplus branch) keeps a small Pythia8 library shim under
# $XDG_CACHE_HOME/faser/pythia8-shim. Reuse the one made on the login node
# if the worker can see it; otherwise let each job make its own in scratch
# (a one-off CVMFS lookup) instead of racing on a shared directory.
export XDG_CACHE_HOME="${XDG_CACHE_HOME:-${HOME:-/nonexistent}/.cache}"
if ! compgen -G "$XDG_CACHE_HOME/faser/pythia8-shim/libpythia8-*.so" > /dev/null; then
  export XDG_CACHE_HOME="$WORKDIR/cache"
fi

# Must be exported BEFORE setup.sh: common_setup.sh only defaults
# FASERDATA if it is not already set.
export FASERDATA="$WORKDIR/data"

###############################################################################
# FASER environment (the repo's own setup.sh)
###############################################################################
set +e +o pipefail
source "$HOMEFASER/setup.sh"
setup_rc=$?
set -eo pipefail
[[ $setup_rc -eq 0 ]]              || { echo "ERROR: setup.sh failed (rc=$setup_rc)" >&2; exit 7; }
command -v root-config >/dev/null  || { echo "ERROR: ROOT not on PATH after setup.sh" >&2; exit 7; }
[[ "$FASERDATA" == "$WORKDIR/data" ]] || { echo "ERROR: setup.sh changed FASERDATA to '$FASERDATA'" >&2; exit 7; }

###############################################################################
# Run
###############################################################################
cmd=(python3 "$HOMEFASER/run_faserps.py" --start-event "$START" --n-events "$NEVT" --seed "$SEED")
if [[ "$INPUTROOT" != "none" ]]; then cmd+=(--input-file "$INPUTROOT"); fi
cmd+=("${EXTRA[@]}")

echo "=== running: ${cmd[*]}"
"${cmd[@]}"
echo "=== run_faserps.py finished ==="

###############################################################################
# Verify, then publish to the shared FASERDATA
###############################################################################
OUTDIR="$FASERDATA/faserG4"
# (find fails if faserps never even created $OUTDIR; don't let pipefail+set -e
# turn that into a silent exit 1 before the explicit check below)
n_out=$({ find "$OUTDIR" -maxdepth 1 -type f -name 'FASERG4-Tcalevent_*.root' 2>/dev/null || true; } | wc -l | tr -d ' ')
if [[ "$n_out" -eq 0 ]]; then
  echo "ERROR: faserps produced no FASERG4-Tcalevent_*.root files in $OUTDIR" >&2
  exit 40
fi
echo "faserps produced $n_out event file(s) (requested $NEVT)"
if [[ "$n_out" -ne "$NEVT" ]]; then
  echo "WARNING: expected $NEVT event files, found $n_out" >&2
fi

mkdir -p "$SHARED/faserG4"
find "$OUTDIR" -maxdepth 1 -type f -name 'FASERG4-Tcalevent_*.root' -exec cp -f -t "$SHARED/faserG4/" {} +

bad=0
while IFS= read -r f; do
  b=$(basename "$f")
  if [[ "$(stat -c %s "$f")" != "$(stat -c %s "$SHARED/faserG4/$b" 2>/dev/null || echo -1)" ]]; then
    echo "ERROR: copy mismatch for $b" >&2
    bad=$((bad + 1))
  fi
done < <(find "$OUTDIR" -maxdepth 1 -type f -name 'FASERG4-Tcalevent_*.root')
if [[ "$bad" -ne 0 ]]; then
  echo "ERROR: $bad file(s) did not copy correctly to $SHARED/faserG4" >&2
  exit 41
fi

# faserps exports the detector geometry on every run; batchreco needs it.
# Publish it once (first finished job wins) so reconstruction jobs find
# $SHARED/GDML/FASERCAL_V10.gdml. Delete it there if you change geometry.
GDML_SRC="$FASERDATA/GDML/FASERCAL_V10.gdml"
if [[ -s "$GDML_SRC" && ! -e "$SHARED/GDML/FASERCAL_V10.gdml" ]]; then
  mkdir -p "$SHARED/GDML"
  tmp="$SHARED/GDML/.FASERCAL_V10.gdml.$$"
  cp "$GDML_SRC" "$tmp" && mv -n "$tmp" "$SHARED/GDML/FASERCAL_V10.gdml"
  rm -f "$tmp"
  echo "published geometry to $SHARED/GDML/FASERCAL_V10.gdml"
fi

echo "=== $n_out file(s) published to $SHARED/faserG4 ==="
echo "=== faserps chunk done: $(date) ==="
