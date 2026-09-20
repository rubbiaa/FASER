export HOMEFASER=$PWD

echo "Setting up environment for FASER simulation"
echo "Current working directory: $HOMEFASER"

source /home/rubbiaa/ROOT/root_install_v6.32.02/bin/thisroot.sh

echo "Root installed in $HOMEFASER/ROOT/root_install"

export GEANT4_INSTALL=/home/rubbiaa/geant4-install/

source $GEANT4_INSTALL/bin/geant4.sh
echo "GEANT4 installed in $GEANT4_INSTALL"

# Let `cmake` (find_package(Geant4)) locate this install without needing
# -DGeant4_DIR=... by hand every time.
export CMAKE_PREFIX_PATH=$GEANT4_INSTALL:$CMAKE_PREFIX_PATH

export PYTHIA8=/home/rubbiaa/ROOT/pythia8312
echo "Pythia8 installed in $PYTHIA8"

# CLHEP/Rave/GenFit are built by FASER's own CMake superbuild
# (cmake/Externals.cmake) into build/external-install/, not into
# top-level *-install directories - point at the real thing.
export CLHEPINSTALL=$HOMEFASER/build/external-install/CLHEP
export RAVEINSTALL=$HOMEFASER/build/external-install/rave
export GENFITINSTALL=$HOMEFASER/build/external-install/GenFit
export LD_LIBRARY_PATH=$GENFITINSTALL/lib:$GENFITINSTALL/lib64:$RAVEINSTALL/lib:$CLHEPINSTALL/lib:$CLHEPINSTALL/lib64:$LD_LIBRARY_PATH

# Executables built by the top-level CMake build (build/bin -
# CMAKE_RUNTIME_OUTPUT_DIRECTORY in CMakeLists.txt) - AnalyReco.exe,
# batchreco.exe, evDisplay.exe, etc.
export PATH=$HOMEFASER/build/bin:$PATH
