# FASER environment setup - auto-detects which known site you're on and
# sets up ROOT/Geant4/Pythia8 accordingly, then hands off to
# common_setup.sh for the logic that's identical everywhere (CLHEP/Rave/
# GenFit paths, LD_LIBRARY_PATH, CMAKE_PREFIX_PATH, PATH, and a sanity
# check). One script, sourced the same way regardless of machine:
#
#     source setup.sh
#
# HOMEFASER is derived from this script's own location (works in both
# bash and zsh), not $PWD - so it works whether you source it from the
# repo root or from somewhere else, including a checkout whose absolute
# path is fixed (e.g. an AFS work area) rather than wherever you happen
# to `cd` from.
#
# Adding a new site: add another `elif` branch below that recognizes your
# machine (a distinctive path, hostname, or username works) and sets up
# ROOT (source thisroot.sh) / GEANT4_INSTALL (export + source geant4.sh) /
# PYTHIA8 for it, the same way the existing branches do.

HOMEFASER="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
export HOMEFASER

# Finds the newest ROOT 6.x release under CVMFS's LCG software area, for
# lxplus accounts that don't keep their own ROOT build inside the
# checkout. Echoes the resolved install directory (the one containing
# bin/thisroot.sh) on success, or nothing if none can be found.
_faser_find_latest_cvmfs_root6() {
  releases_dir=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT
  [ -d "$releases_dir" ] || return 0

  version=$(ls "$releases_dir" 2>/dev/null | grep -E '^6\.[0-9]+\.[0-9]+$' | sort -V | tail -1)
  [ -n "$version" ] || return 0

  # Several platform builds usually exist per version - one per OS/
  # compiler combination CERN builds for - and which one matches this
  # machine can vary node to node. Any optimized (non-debug) x86_64 one
  # works for our purposes; pick the last alphabetically for a stable,
  # deterministic choice rather than whatever `ls` happens to return first.
  platform=$(ls "$releases_dir/$version" 2>/dev/null | grep -E '^x86_64-.*-opt$' | sort | tail -1)
  [ -n "$platform" ] || return 0

  candidate="$releases_dir/$version/$platform"
  [ -f "$candidate/bin/thisroot.sh" ] && echo "$candidate"
}

if [ -d /Users/rubbiaa/Documents/GitHub/GEANT4/geant4-v11.4.2-install ]; then
  # André's Mac (Apple Silicon)
  echo "FASER setup: detected site = André's Mac"
  source /Users/rubbiaa/Documents/GitHub/ROOT/root_install/bin/thisroot.sh

  export PYTHIA8=$HOMEFASER/pythia8312

  export GEANT4_INSTALL=/Users/rubbiaa/Documents/GitHub/GEANT4/geant4-v11.4.2-install/
  source $GEANT4_INSTALL/bin/geant4.sh

elif [ -d /home/rubbiaa/geant4-install ]; then
  # Ubuntu box (rubbiaa, Ryzen)
  echo "FASER setup: detected site = Ubuntu (Ryzen)"
  echo "Current working directory: $HOMEFASER"

  source /home/rubbiaa/ROOT/root_install_v6.32.02/bin/thisroot.sh
  echo "Root installed in $HOMEFASER/ROOT/root_install"

  export GEANT4_INSTALL=/home/rubbiaa/geant4-install/
  source $GEANT4_INSTALL/bin/geant4.sh
  echo "GEANT4 installed in $GEANT4_INSTALL"

  export PYTHIA8=/home/rubbiaa/ROOT/pythia8312
  echo "Pythia8 installed in $PYTHIA8"

elif [ -d /cvmfs/geant4.cern.ch ]; then
  # lxplus: prefer a ROOT build local to the checkout (root-install/) if
  # one is there, since that's what building your own ROOT from source
  # into the repo implies you want used. Otherwise, fall back to
  # auto-discovering the latest ROOT 6 release CVMFS itself publishes,
  # rather than requiring one to be built locally at all.
  echo "FASER setup: detected site = lxplus"
  echo "Current working directory: $HOMEFASER"

  if [ -f "$HOMEFASER/root-install/bin/thisroot.sh" ]; then
    echo "ROOT: using the local build at $HOMEFASER/root-install"
    source $HOMEFASER/root-install/bin/thisroot.sh
  else
    _faser_cvmfs_root=$(_faser_find_latest_cvmfs_root6)
    if [ -n "$_faser_cvmfs_root" ]; then
      echo "ROOT: no local build at $HOMEFASER/root-install - using the latest CVMFS release, $_faser_cvmfs_root"
      source "$_faser_cvmfs_root/bin/thisroot.sh"
    else
      echo "FASER setup: no ROOT found - neither a local build at"
      echo "  $HOMEFASER/root-install nor a ROOT 6 release under"
      echo "  /cvmfs/sft.cern.ch/lcg/app/releases/ROOT. Source ROOT yourself"
      echo "  before this script, or build one into root-install/."
      return 1 2>/dev/null || exit 1
    fi
    unset _faser_cvmfs_root
  fi

  export GEANT4_INSTALL=/cvmfs/geant4.cern.ch/geant4/11.2.p01/x86_64-el9-gcc11-optdeb
  pushd . > /dev/null
  cd $GEANT4_INSTALL/bin
  source geant4.sh
  popd > /dev/null
  echo "GEANT4 installed in $GEANT4_INSTALL"

  export PYTHIA8=$HOMEFASER/pythia8312
  echo "Pythia8 installed in $PYTHIA8"

else
  echo "FASER setup: could not auto-detect a known site."
  echo "  Checked for: André's Mac, the Ubuntu (Ryzen) box, and lxplus (/cvmfs/geant4.cern.ch)."
  echo "  If this is a new machine, add an elif branch for it near the top of"
  echo "  setup.sh (source thisroot.sh, export + source GEANT4_INSTALL,"
  echo "  export PYTHIA8) - or set those by hand right now and just source"
  echo "  common_setup.sh yourself:"
  echo "    source $HOMEFASER/common_setup.sh"
  return 1 2>/dev/null || exit 1
fi

source $HOMEFASER/common_setup.sh
