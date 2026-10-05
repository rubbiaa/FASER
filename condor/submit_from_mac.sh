#!/bin/bash
###############################################################################
# Drive the lxplus HTCondor submission from your Mac.
#
# HTCondor jobs can only be submitted from lxplus, so this just (1) rsyncs
# this condor/ directory into your FASER checkout on lxplus and (2) runs
# submit_jobs.sh there over ssh. Everything is built and run on lxplus from
# that checkout (sources, build/ and setup.sh come from git, not from here).
#
#   ./submit_from_mac.sh faserps   --run 10000 --detector 3DCAL --total-events 350000
#   ./submit_from_mac.sh faserps   --total-events 100000 -- --muons --muon-momentum-gev 250
#   ./submit_from_mac.sh batchreco --run 10000
#   ./submit_from_mac.sh status                # condor_q
#   ./submit_from_mac.sh sync                  # only copy condor/ over
#   ./submit_from_mac.sh pull                  # git pull --ff-only in the lxplus checkout
#   ./submit_from_mac.sh ssh <command...>      # run anything in the lxplus checkout
#
# (faserps / batchreco take the same options as submit_jobs.sh - add
#  --dry-run to see the job list and the condor_submit command first.)
#
# Settings, from the environment or ~/.faser_condor.conf on the Mac:
#   LXPLUS_HOST     ssh host alias        (default: lxplus)
#   LXPLUS_FASER    FASER checkout there  (default: /afs/cern.ch/work/r/rubbiaa/FASER)
# The shared data area on lxplus (EOS) is FASER_SHARED_DATA in
# ~/.faser_condor.conf *on lxplus*, or pass --faserdata DIR here.
#
# Tip: with `ControlMaster auto` / `ControlPersist` for the lxplus host in
# ~/.ssh/config you only do the password + 2FA once for the several ssh and
# rsync calls this makes.
###############################################################################
set -euo pipefail

# shellcheck disable=SC1090
[[ -f "$HOME/.faser_condor.conf" ]] && source "$HOME/.faser_condor.conf"
LXPLUS_HOST="${LXPLUS_HOST:-lxplus}"
LXPLUS_FASER="${LXPLUS_FASER:-/afs/cern.ch/work/r/rubbiaa/FASER}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

die() { echo "ERROR: $*" >&2; exit 1; }
usage() { sed -n '2,/^####/p' "${BASH_SOURCE[0]}" | sed -e '1d' -e '$d' -e 's/^# \{0,1\}//'; }

[[ $# -ge 1 ]] || { usage; exit 1; }
CMD="$1"; shift

sync_condor() {
  echo "[mac] syncing $HERE -> $LXPLUS_HOST:$LXPLUS_FASER/condor/"
  ssh "$LXPLUS_HOST" "test -d $(printf '%q' "$LXPLUS_FASER") && mkdir -p $(printf '%q' "$LXPLUS_FASER/condor")" \
    || die "no FASER checkout at $LXPLUS_HOST:$LXPLUS_FASER (set LXPLUS_FASER; clone https://github.com/rubbiaa/FASER there and build it first)"
  # logs/ and jobs/ live on lxplus only; keep them out of the sync both ways
  rsync -az --delete --exclude 'logs/' --exclude 'jobs/' --exclude '.DS_Store' \
    "$HERE/" "$LXPLUS_HOST:$LXPLUS_FASER/condor/"
}

remote() {  # remote <cwd-relative-to-checkout> <command...>
  local dir="$1"; shift
  ssh "$LXPLUS_HOST" "cd $(printf '%q' "$LXPLUS_FASER/$dir") && $*"
}

quote_args() { local out=""; for a in "$@"; do out+="$(printf '%q ' "$a")"; done; printf '%s' "$out"; }

case "$CMD" in
  faserps|batchreco)
    sync_condor
    remote condor "bash ./submit_jobs.sh $CMD $(quote_args "$@")"
    ;;
  status)  ssh "$LXPLUS_HOST" "condor_q $(quote_args "$@")" ;;
  sync)    sync_condor ;;
  pull)    remote . "git pull --ff-only" ;;
  ssh)     remote . "$(quote_args "$@")" ;;
  -h|--help|help) usage ;;
  *)       usage; die "unknown command '$CMD'" ;;
esac
