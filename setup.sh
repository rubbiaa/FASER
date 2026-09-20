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
# repo root or from somewhere else. That also removes the old need for a
# hardcoded absolute HOMEFASER on Anna's lxplus/AFS checkout: it's
# whatever directory this file actually lives in.
#
# Adding a new site: add another `elif` branch below that recognizes your
# machine (a distinctive path, hostname, or username works) and sets up
# ROOT (source thisroot.sh) / GEANT4_INSTALL (export + source geant4.sh) /
# PYTHIA8 for it, the same way the existing branches do.

HOMEFASER="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
export HOMEFASER

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

elif [ "$(whoami)" = "amascell" ] && [ -d /cvmfs/geant4.cern.ch ]; then
  # Anna's lxplus/AFS checkout - a CVMFS-provided ROOT release, rather
  # than a ROOT built locally into the checkout like generic lxplus below.
  echo "FASER setup: detected site = lxplus (Anna)"
  source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/bin/thisroot.sh

  export GEANT4_INSTALL=/cvmfs/geant4.cern.ch/geant4/11.2.p01/x86_64-el9-gcc11-optdeb
  source $GEANT4_INSTALL/bin/geant4.sh

  export PYTHIA8=$HOMEFASER/pythia8312

elif [ -d /cvmfs/geant4.cern.ch ]; then
  # Generic lxplus: ROOT built locally into the checkout (root-install/),
  # Geant4 from CVMFS.
  echo "FASER setup: detected site = lxplus"
  echo "Current working directory: $HOMEFASER"

  source $HOMEFASER/root-install/bin/thisroot.sh
  echo "Root installed in $HOMEFASER/ROOT/root_install"

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
