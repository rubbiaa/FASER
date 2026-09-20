export HOMEFASER=$PWD

source /Users/rubbiaa/Documents/GitHub/ROOT/root_install/bin/thisroot.sh

export PYTHIA8=$HOMEFASER/pythia8312

export CLHEPINSTALL=$HOMEFASER/CLHEP-install
export RAVEINSTALL=$HOMEFASER/rave-install
export GENFITINSTALL=$HOMEFASER/GenFit-install
export LD_LIBRARY_PATH=$GENFITINSTALL/lib:$GENFITINSTALL/lib64:$RAVEINSTALL/lib:$CLHEPINSTALL/lib:$LD_LIBRARY_PATH

export GEANT4_INSTALL=/Users/rubbiaa/Documents/GitHub/GEANT4/geant4-v11.4.2-install/
source $GEANT4_INSTALL/bin/geant4.sh

# Let `cmake` (find_package(Geant4)) locate this install without needing
# -DGeant4_DIR=... by hand every time.
export CMAKE_PREFIX_PATH=$GEANT4_INSTALL:$CMAKE_PREFIX_PATH

# Executables built by the top-level CMake build (build/bin -
# CMAKE_RUNTIME_OUTPUT_DIRECTORY in CMakeLists.txt) - AnalyReco.exe,
# batchreco.exe, evDisplay.exe, etc.
export PATH=$HOMEFASER/build/bin:$PATH
