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
#
# FASERDATA (where FASERG4/batchreco output lives - see
# CoreUtils/FaserDataDir.hh) works the same way: common_setup.sh defaults
# it to $HOMEFASER/data, but only if it isn't already set, so a site
# branch here can `export FASERDATA=/some/other/path` before
# common_setup.sh runs to point a specific machine - or a checkout
# dedicated to a specific run/target - at its own data area instead.

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

if [ -d /Users/rubbiaa/Documents/GitHub/GEANT4/geant4-v11.4.3-install ]; then
  # André's Mac (Apple Silicon)
  echo "FASER setup: detected site = André's Mac"
  source /Users/rubbiaa/Documents/GitHub/ROOT/root_install/bin/thisroot.sh

  export GEANT4_INSTALL=/Users/rubbiaa/Documents/GitHub/GEANT4/geant4-v11.4.3-install/
  source $GEANT4_INSTALL/bin/geant4.sh

  export PYTHIA8=/Users/rubbiaa/Documents/GitHub/ROOT/pythia8312
  echo "Pythia8 installed in $PYTHIA8"

  # The ROOT build sourced above (root_install) is linked against Homebrew's
  # Python 3.14, not whatever "python3" happens to resolve to interactively
  # (e.g. an active conda "(base)" env shadows it with a different Python
  # minor version -- PyROOT refuses to import across a minor-version
  # mismatch). run_regression_tests.py's --python defaults to this when set,
  # so PyROOT scripts (summarize_output.py and friends) pick up the right
  # interpreter without having to pass --python by hand every time.
  export FASER_PYTHON=/opt/homebrew/bin/python3.14
  echo "PyROOT scripts will default to FASER_PYTHON=$FASER_PYTHON"

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

  # This particular Geant4 CVMFS release is itself dynamically linked
  # against a standalone CLHEP from this LCG view (confirmed via
  # `ldd $GEANT4_INSTALL/lib64/libG4global.so`), not a copy bundled inside
  # GEANT4_INSTALL - so cmake/Externals.cmake's bundled-CLHEP auto-detect
  # doesn't find it on its own. Exporting $CLHEP_ROOT here makes
  # Externals.cmake default to reusing this exact CLHEP instead of
  # building FASER's own from source, which would otherwise be a second,
  # ABI-incompatible CLHEP build that fails to link in any test/executable
  # pulling in both GenFit and Geant4/ROOT at once (observed: a
  # CLHEP_SINGLE_THREAD-vs-not TLS/.bss mismatch on
  # CLHEP::RandGaussZiggurat's internal state).
  export CLHEP_ROOT=/cvmfs/sft.cern.ch/lcg/views/LCG_104b_geant4ext20231106/x86_64-el9-gcc11-opt
  echo "CLHEP: reusing the standalone CLHEP at $CLHEP_ROOT (the one this Geant4 release is itself linked against)"

  # Work around a system-Python quirk observed on lxplus (AlmaLinux9):
  # once thisroot.sh has added ROOT's own lib/ to $PYTHONPATH, the system
  # /usr/bin/python3 (3.9, a non-default/secondary interpreter on EL9)
  # starts resolving its own compiled stdlib extensions (lib-dynload) from
  # /usr/lib/python3.9/lib-dynload/ - the 32-bit (i686) multilib copy -
  # instead of /usr/lib64/python3.9/lib-dynload/, the real 64-bit one this
  # binary needs. Pure-Python stdlib modules (e.g. csv.py) still resolve
  # fine either way, but any module backed by a compiled extension (e.g.
  # csv's own `_csv`) then fails with "ModuleNotFoundError: No module
  # named '_csv'" - with no PYTHONHOME/PYTHONPLATLIBDIR override and no
  # shadowing file anywhere involved; a bare `python3 -c "import csv"`
  # reproduces it as soon as that one PYTHONPATH entry is present.
  # Prepending the correct 64-bit lib-dynload directory to $PYTHONPATH
  # fixes it (confirmed on lxplus969); compute it rather than hardcoding
  # "3.9" so this keeps working if lxplus's system Python version changes.
  _faser_py_ver=$(python3 -c 'import sys; print("%d.%d" % sys.version_info[:2])' 2>/dev/null)
  if [ -n "$_faser_py_ver" ] && [ -d "/usr/lib64/python${_faser_py_ver}/lib-dynload" ]; then
    export PYTHONPATH="/usr/lib64/python${_faser_py_ver}/lib-dynload:${PYTHONPATH}"
  fi
  unset _faser_py_ver

  # Flag checked below, after common_setup.sh has had its turn at
  # $LD_LIBRARY_PATH too - see the matching block at the end of this file.
  _faser_site_is_lxplus=1

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

if [ -n "$_faser_site_is_lxplus" ]; then
  # Second half of the lxplus Python fix above: common_setup.sh (just
  # sourced) prepends GenFit/Rave/CLHEP onto $LD_LIBRARY_PATH, but by then
  # it already contains the CLHEP_ROOT LCG view's own lib/ dir (added
  # while resolving the reused CLHEP, earlier in this branch) - and that
  # view bundles its own standalone libpython3.9.so.1.0 (a plain upstream
  # build, CVMFS "Python" release 3.9.12-9a1bc). Found ahead of the
  # system's own /usr/lib64/libpython3.9.so.1.0 (which carries a Red
  # Hat-only backported symbol, _PyModule_AddObjectRef), that view's copy
  # gets dlopen'd instead for any later python3 extension module needing
  # libpython - e.g. hashlib's _hashlib.so - and fails with
  # "undefined symbol: _PyModule_AddObjectRef" (confirmed on lxplus969).
  # Pure stdlib Python doesn't care which libpython it runs under, so
  # putting the system's own lib64 first is safe, and doesn't affect any
  # FASER/CLHEP/Rave/GenFit/Geant4/ROOT library resolution - none of
  # those names collide with anything under /usr/lib64.
  if [ -d /usr/lib64 ]; then
    export LD_LIBRARY_PATH="/usr/lib64:${LD_LIBRARY_PATH}"
  fi
  unset _faser_site_is_lxplus
fi
