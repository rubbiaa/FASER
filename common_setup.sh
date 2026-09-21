# FASER environment setup - shared logic
#
# Sourced by setup.sh once it has auto-detected which site you're on and
# handled the handful of things that genuinely differ per site (where
# ROOT's thisroot.sh lives, GEANT4_INSTALL, HOMEFASER). Everything past
# that point - PYTHIA8's default, where CLHEP/Rave/GenFit's own
# from-source builds actually land, LD_LIBRARY_PATH, CMAKE_PREFIX_PATH,
# PATH, and the sanity check below - is identical across every site, so it
# lives here once instead of inside each of setup.sh's per-site branches.
# PYTHIA8 only needs to be set in a setup.sh site branch if that machine
# keeps its own standalone Pythia8 build (see the default below).
#
# Adding a new site: see setup.sh - add an elif branch there that sets
# ROOT/GEANT4_INSTALL, it already ends with
# `source $HOMEFASER/common_setup.sh`.

: "${HOMEFASER:=$PWD}"
export HOMEFASER

# Where FASERG4/batchreco/etc. read and write simulation/
# reconstruction data - see CoreUtils/FaserDataDir.hh, which every
# consumer (TcalEvent, BatchReco.cc, DumpHits.cc, BatchReco_DetResp.cc,
# AnalyReco.cc, FileMask.cc, MyMainFrame.cc) resolves through instead of
# a hardcoded relative "input/"/"output/" path or a symlink between
# subdirectories. Defaults to a `data/` directory inside the checkout,
# but - same pattern as PYTHIA8 below - only if a per-site branch in
# setup.sh hasn't already exported it to something else first (e.g. a
# site that wants large output on scratch/EOS instead of inside the git
# checkout, or a checkout dedicated to a specific run/target that wants
# its own data area). GetDataDir() creates missing subdirectories
# itself, so nothing here needs to mkdir them.
: "${FASERDATA:=$HOMEFASER/data}"
export FASERDATA
# Created eagerly (not left for the first executable that needs it) so
# the sanity check below can treat FASERDATA as required, the same way
# HOMEFASER is - CoreUtils/FaserDataDir.hh still creates faserG4/batch/
# etc. subdirectories on demand.
mkdir -p "$FASERDATA"

# CLHEP/Rave/GenFit are built by FASER's own CMake superbuild
# (cmake/Externals.cmake) into build/external-install/, not into
# top-level *-install directories - point at the real thing.
export CLHEPINSTALL=$HOMEFASER/build/external-install/CLHEP
export RAVEINSTALL=$HOMEFASER/build/external-install/rave
export GENFITINSTALL=$HOMEFASER/build/external-install/GenFit
export LD_LIBRARY_PATH=$GENFITINSTALL/lib:$GENFITINSTALL/lib64:$RAVEINSTALL/lib:$CLHEPINSTALL/lib:$CLHEPINSTALL/lib64:$LD_LIBRARY_PATH

# Pythia8 is built by the same CMake superbuild, but unlike CLHEP/Rave/
# GenFit it's never `make install`-ed into build/external-install/ - its
# own build places headers/libs directly under the extracted source tree,
# which ExternalProject_Add stages under build/external/pythia8/src/. Point
# PYTHIA8 there by default, but only if a per-site branch in setup.sh
# hasn't already set it to something else (e.g. a standalone Pythia8 build
# that predates the CMake superbuild, kept around on a specific machine).
: "${PYTHIA8:=$HOMEFASER/build/external/pythia8/src/pythia8_external}"
export PYTHIA8

# Let `cmake` (find_package(Geant4) in FASERG4/FASERCalProtoG4) locate the
# Geant4 install without needing -DGeant4_DIR=... by hand every time.
if [ -n "$GEANT4_INSTALL" ]; then
  export CMAKE_PREFIX_PATH=$GEANT4_INSTALL:$CMAKE_PREFIX_PATH
fi

# Executables built by the top-level CMake build (build/bin -
# CMAKE_RUNTIME_OUTPUT_DIRECTORY in CMakeLists.txt) - AnalyReco.exe,
# batchreco.exe, evDisplay.exe, etc.
export PATH=$HOMEFASER/build/bin:$PATH

# -----------------------------------------------------------------------------
# Sanity check
# -----------------------------------------------------------------------------
# Catches the two failure modes that have repeatedly cost real debugging
# time on this project: a variable that's silently unset (a script forgets
# `export`, or a copy-pasted script skips a step - this is exactly how the
# Ubuntu box's missing `export GEANT4_INSTALL` went unnoticed for a while),
# and a path that's set but stale or wrong (pointing at an install that's
# since moved, or that was never built). Report clearly right here, instead
# of it surfacing 10 minutes later as a confusing CMake or link failure.
#
# CLHEPINSTALL/RAVEINSTALL/GENFITINSTALL are deliberately "optional": they
# won't exist yet on a fresh checkout before the first `cmake --build`, and
# that's expected, not an error.
faser_setup_ok=1
faser_check() {
  # $1 = label, $2 = value, $3 = required|optional
  if [ -z "$2" ]; then
    if [ "$3" = "required" ]; then
      echo "  [MISSING]   $1 is not set"
      faser_setup_ok=0
    else
      echo "  [ -- ]      $1 not set (optional)"
    fi
  elif [ ! -e "$2" ]; then
    echo "  [NOT FOUND] $1 = $2  (path does not exist)"
    if [ "$3" = "required" ]; then faser_setup_ok=0; fi
  else
    echo "  [ok]        $1 = $2"
  fi
}

echo "FASER environment check:"
faser_check "HOMEFASER    " "$HOMEFASER"      required
faser_check "FASERDATA    " "$FASERDATA"      required
faser_check "GEANT4_INSTALL" "$GEANT4_INSTALL" optional
faser_check "PYTHIA8      " "$PYTHIA8"        optional
faser_check "CLHEPINSTALL " "$CLHEPINSTALL"   optional
faser_check "RAVEINSTALL  " "$RAVEINSTALL"    optional
faser_check "GENFITINSTALL" "$GENFITINSTALL"  optional

if command -v root-config >/dev/null 2>&1; then
  echo "  [ok]        ROOT          = $(root-config --version)  ($(command -v root-config))"
else
  echo "  [MISSING]   ROOT not found on PATH - source thisroot.sh before this script"
  faser_setup_ok=0
fi

if [ "$faser_setup_ok" = "1" ]; then
  echo "FASER environment OK."
else
  echo "FASER environment INCOMPLETE - see [MISSING] lines above."
fi
unset faser_setup_ok
