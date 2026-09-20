
source /Users/ukose/sw/kits/root-install/bin/thisroot.sh

source /usr/local/bin/geant4.sh

export PYTHIA8=/Users/ukose/sw/kits/pythia8312
export PYTHIA8DATA=$PYTHIA8/xmldoc  # Where PYTHIA’s XML data is stored
export LD_LIBRARY_PATH=$PYTHIA8/lib:$LD_LIBRARY_PATH  # Or DYLD_LIBRARY_PATH on macOS

export CLHEPINSTALL=/Users/ukose/sw/kits/CLHEP-install
export RAVEINSTALL=/Users/ukose/sw/kits/rave-install
export GENFITINSTALL=/Users/ukose/sw/kits/GenFit-install

export DYLD_LIBRARY_PATH=/opt/homebrew/opt/boost/lib:$DYLD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$CLHEPINSTALL/lib:$DYLD_LIBRARY_PATH
export LD_LIBRARY_PATH=$GENFITINSTALL/lib:$RAVEINSTALL/lib:$CLHEPINSTALL/lib:$LD_LIBRARY_PATH


# Executables built by the top-level CMake build (build/bin -
# CMAKE_RUNTIME_OUTPUT_DIRECTORY in CMakeLists.txt) - AnalyReco.exe,
# batchreco.exe, evDisplay.exe, etc. This script doesn't define a
# HOMEFASER-style repo-root variable like the other setup scripts, so
# this assumes it's sourced from the repo root (as the others expect too).
export PATH=$PWD/build/bin:$PATH
