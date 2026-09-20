# Site setup: Anna's lxplus/AFS checkout (CVMFS ROOT + Geant4, fixed repo
# path since this doesn't run from the repo root). Sets the handful of
# things that are genuinely specific to this checkout, then hands off to
# common_setup.sh for everything shared across all sites.
export HOMEFASER=/afs/cern.ch/work/a/amascell/FASERCal/debug/FASER

source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/bin/thisroot.sh

# Previously this sourced geant4.sh without ever exporting GEANT4_INSTALL
# itself, so cmake/Externals.cmake's smart CLHEP-reuse-from-Geant4 default
# (which keys off $GEANT4_INSTALL) never triggered here - CLHEP was always
# built from source instead of reusing CVMFS's own bundled copy, unlike on
# the Mac/Ubuntu setups. Exporting it properly fixes that.
export GEANT4_INSTALL=/cvmfs/geant4.cern.ch/geant4/11.2.p01/x86_64-el9-gcc11-optdeb
source $GEANT4_INSTALL/bin/geant4.sh

export PYTHIA8=$HOMEFASER/pythia8312

source $HOMEFASER/common_setup.sh
