#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.13...3.27)
project(FASERPS)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all gdml)
else()
  find_package(Geant4 REQUIRED)
endif()

# Set the optimization level to O3 and C++ standard to C++17
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -g")

# Specify the external directories
set(EXTERNAL_INCLUDE_DIR "../CoreUtils")
set(EXTERNAL_SOURCE_DIR "../CoreUtils")

include_directories(${EXTERNAL_INCLUDE_DIR})

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

# Ensure Geant4 is compiled with C++17 support
set(GEANT4_COMPILE_FEATURES cxx_std_17)

find_package(ROOT REQUIRED COMPONENTS Geom)

# --- MuonDIS: Pythia8 (on-the-fly muon-nucleon DIS event generation) -------
# Built via the top-level Makefile's `pythia8` target into ../pythia8312
# relative to this CMakeLists.txt (i.e. fasermuondis/pythia8312). Adjust
# PYTHIA8_DIR below if your Pythia8 install lives elsewhere.
set(PYTHIA8_DIR "${CMAKE_SOURCE_DIR}/../pythia8312" CACHE PATH "Path to the Pythia8 install used by MuonDIS")
include_directories(${PYTHIA8_DIR}/include)
set(PYTHIA8_LIBRARY "${PYTHIA8_DIR}/lib/libpythia8.a")
if(NOT EXISTS ${PYTHIA8_LIBRARY})
  message(WARNING "MuonDIS: Pythia8 static library not found at ${PYTHIA8_LIBRARY} -- "
                   "build it first (see fasermuondis/Makefile's `pythia8` target), or set "
                   "-DPYTHIA8_DIR=/path/to/pythia8312 when configuring CMake.")
endif()
# ----------------------------------------------------------------------------

include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})

set(SOURCES
    src/Track.cc
    src/ActionInitialization.cc
    src/DetectorConstruction.cc
    src/MuonDetMagneticField.cc
    src/DetectorMessenger.cc
    src/EventAction.cc
    src/ParticleManager.cc
    src/PrimaryGeneratorAction.cc
    src/PrimaryGeneratorMessenger.cc
    src/RunAction.cc
    src/TrackerSD.cc
    src/MDTSD.cc
    src/TauDecayPhysics.cc
    src/CustomTauDecay.cc
    src/MuonDISPhysics.cc
    src/MuonDISMessenger.cc
    src/MuonDISMuonNuclearModel.cc
    src/MuonDISNuclearWrapperProcess.cc
    src/MuonDISInteractionRecorder.cc
    src/MuonDISPythiaGenerator.cc
    ${EXTERNAL_SOURCE_DIR}/TcalEvent.cc
    ${EXTERNAL_SOURCE_DIR}/TPOEvent.cc
    faserps.cc
    )

file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Generate ROOT dictionary
ROOT_GENERATE_DICTIONARY(G__PS ${EXTERNAL_SOURCE_DIR}/TPOEvent.hh ${EXTERNAL_SOURCE_DIR}/TcalEvent.hh ${PROJECT_SOURCE_DIR}/include/Geant4Process.hh LINKDEF ${PROJECT_SOURCE_DIR}/include/LinkDef.h)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(faserps ${SOURCES} ${os_headers} G__PS.cxx)
target_link_libraries(faserps ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} EG ROOTDataFrame ROOTNTuple ${PYTHIA8_LIBRARY} ${CMAKE_DL_LIBS})

#----------------------------------------------------------------------------
# Build the shared library for the classes mentioned in the LinkDef.h file
#
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# add_library(FaserPS SHARED G__PS.cxx)

#----------------------------------------------------------------------------
# Copy the stuff from the ./input directory to the build directory
#
file(GLOB input_files ${PROJECT_SOURCE_DIR}/input/*)
file(COPY ${input_files} DESTINATION ${PROJECT_BINARY_DIR})
