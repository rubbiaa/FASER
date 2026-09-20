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
# ROOT (source thisroot.sh) / GEANT4_INSTALL (export + source geant4.sh) the
# same way the existing branches do. PYTHIA8 doesn't need to be set here at
# all unless your site keeps its own standalone Pythia8 build outside the
# checkout (see common_setup.sh for the default) - the Ubuntu branch below
# is the one example of that.

HOMEFASER="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
export HOMEFASER

# Best-effort detection of a regex fragment matching this machine's OS in
# CVMFS release directory names. Echoes nothing if it can't be determined
# - callers treat that as "no OS preference", not as an error.
#
# Different CVMFS/LCG release areas spell the same OS differently: Geant4's
# own releases use the short "el9" form for RHEL/Alma/Rocky 9, but ROOT's
# release area names it after the actual rebuild distro instead, e.g.
# "almalinux9.8" (observed on lxplus9, an el9/RHEL9.8 machine) rather than
# "el9" - so a single fixed spelling isn't enough. Match any of the common
# spellings for this major version instead of guessing one.
_faser_os_platform_tag() {
  if [ -f /etc/os-release ]; then
    ( . /etc/os-release
      major="${VERSION_ID%%.*}"
      case "$ID" in
        rhel|almalinux|rocky|centos)
          echo "(el${major}|almalinux${major}|rocky${major}|centos${major})"
          ;;
        ubuntu)
          # Observed spelling keeps the dot, e.g. "ubuntu22.04" - not "ubuntu2204".
          echo "ubuntu${VERSION_ID}"
          ;;
      esac
    )
  fi
}

# Finds a ROOT 6.x release under CVMFS's LCG software area, for lxplus
# accounts that don't keep their own ROOT build inside the checkout.
# Echoes the resolved install directory (the one containing
# bin/thisroot.sh) on success, or nothing if none can be found.
_faser_find_latest_cvmfs_root6() {
  releases_dir=/cvmfs/sft.cern.ch/lcg/app/releases/ROOT
  [ -d "$releases_dir" ] || return 0

  # Newest-to-oldest, so the first match found below is also the newest.
  versions=$(ls "$releases_dir" 2>/dev/null | grep -E '^6\.[0-9]+\.[0-9]+$' | sort -Vr)
  [ -n "$versions" ] || return 0

  # Several platform builds usually exist per version - one per OS/
  # compiler combination CERN builds for. Picking whichever sorts last
  # alphabetically (the old approach) is a trap: CVMFS also publishes
  # builds for other OSes (e.g. Ubuntu), and their platform strings can
  # sort AFTER "el9" (u > e) even on an actual el9 lxplus node - silently
  # picking an ABI-incompatible ROOT that then fails to configure against
  # system/CVMFS dependencies (VDT and friends) further down the line.
  # Prefer a same-OS build instead, and fall back across older versions
  # (not just the newest one) until one is found, since the newest
  # version doesn't always have a build for every OS yet.
  os_tag=$(_faser_os_platform_tag)
  if [ -n "$os_tag" ]; then
    for version in $versions; do
      platform=$(ls "$releases_dir/$version" 2>/dev/null | grep -E "^x86_64-${os_tag}.*-opt$" | sort | tail -1)
      if [ -n "$platform" ]; then
        candidate="$releases_dir/$version/$platform"
        if [ -f "$candidate/bin/thisroot.sh" ]; then
          echo "$candidate"
          return 0
        fi
      fi
    done
  fi

  # No same-OS build found for any version (or the OS couldn't be
  # detected) - fall back to the newest version's alphabetically-last
  # optimized x86_64 build, same as the old behaviour. Better than
  # failing outright, though it may still be for a different OS's ABI.
  version=$(echo "$versions" | head -1)
  platform=$(ls "$releases_dir/$version" 2>/dev/null | grep -E '^x86_64-.*-opt$' | sort | tail -1)
  [ -n "$platform" ] || return 0

  candidate="$releases_dir/$version/$platform"
  [ -f "$candidate/bin/thisroot.sh" ] && echo "$candidate"
}

if [ -d /Users/rubbiaa/Documents/GitHub/GEANT4/geant4-v11.4.2-install ]; then
  # André's Mac (Apple Silicon)
  echo "FASER setup: detected site = André's Mac"
  source /Users/rubbiaa/Documents/GitHub/ROOT/root_install/bin/thisroot.sh

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

else
  echo "FASER setup: could not auto-detect a known site."
  echo "  Checked for: André's Mac, the Ubuntu (Ryzen) box, and lxplus (/cvmfs/geant4.cern.ch)."
  echo "  If this is a new machine, add an elif branch for it near the top of"
  echo "  setup.sh (source thisroot.sh, export + source GEANT4_INSTALL) -"
  echo "  or set those by hand right now and just source common_setup.sh"
  echo "  yourself:"
  echo "    source $HOMEFASER/common_setup.sh"
  return 1 2>/dev/null || exit 1
fi

source $HOMEFASER/common_setup.sh
