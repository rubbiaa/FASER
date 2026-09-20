export HOMEFASER=/afs/cern.ch/work/a/amascell/FASERCal/debug/FASER

source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.32.02/x86_64-almalinux9.4-gcc114-opt/bin/thisroot.sh
source /cvmfs/geant4.cern.ch/geant4/11.2.p01/x86_64-el9-gcc11-optdeb/bin/geant4.sh

export PYTHIA8=$HOMEFASER/pythia8312

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
