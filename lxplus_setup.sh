# Site setup: lxplus (CVMFS ROOT + Geant4). Sets the handful of things that
# are genuinely specific to this site, then hands off to common_setup.sh
# for everything shared across all sites.
export HOMEFASER=$PWD

echo "Setting up environment for FASER simulation"
echo "Current working directory: $HOMEFASER"

source $HOMEFASER/root-install/bin/thisroot.sh
echo "Root installed in $HOMEFASER/ROOT/root_install"

# Previously this cd'd into the CVMFS bin/ dir and sourced geant4.sh
# without ever exporting GEANT4_INSTALL itself, so the echo below printed
# an empty value and cmake/Externals.cmake's smart CLHEP-reuse-from-Geant4
# default (which keys off $GEANT4_INSTALL) never triggered on lxplus,
# unlike on the Mac/Ubuntu setups - CLHEP was always built from source
# here instead of reusing CVMFS's own bundled copy. Exporting it properly
# fixes both.
export GEANT4_INSTALL=/cvmfs/geant4.cern.ch/geant4/11.2.p01/x86_64-el9-gcc11-optdeb
pushd . > /dev/null
cd $GEANT4_INSTALL/bin
source geant4.sh
popd > /dev/null
echo "GEANT4 installed in $GEANT4_INSTALL"

export PYTHIA8=$HOMEFASER/pythia8312
echo "Pythia8 installed in $PYTHIA8"

source $HOMEFASER/common_setup.sh
