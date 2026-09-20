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
target_link_libraries(faserps ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} EG ROOTDataFrame ROOTNTuple)

# Geant4's own libraries (specifically the analysis-manager component
# bundled into ${Geant4_LIBRARIES}) can be built against a completely
# different ROOT install than the one FASER's own dictionaries are
# compiled against. On lxplus, Geant4 11.2.p01 is built alongside (and
# links its optional ROOT/HDF5/Xerces-C analysis output drivers against)
# the LCG "geant4ext" view's own bundled ROOT, while FASER's own
# find_package(ROOT) above resolves a separate, newer standalone ROOT
# release for faserps's own TPOEvent/TcalEvent dictionaries. faserps never
# calls into Geant4's built-in ROOT analysis-manager output driver, so
# the two never need to agree with each other at runtime - but the GNU
# linker still refuses to produce the executable unless it can fully
# resolve every symbol referenced from within that unused ROOT library,
# which fails here because the two ROOT versions' TGenericClassInfo
# constructor signature differs (undefined reference to
# ROOT::TGenericClassInfo::TGenericClassInfo(...) from libROOTHist.so).
# --allow-shlib-undefined tells the linker that's fine: an unresolved
# symbol inside a transitively-linked shared library only matters if
# something we actually call needs it. Linux (GNU ld) only - ld64 on
# macOS doesn't understand this flag and doesn't need it, since the Mac
# build uses a single self-built ROOT/Geant4 pair with no such clash.
if(UNIX AND NOT APPLE)
  target_link_options(faserps PRIVATE -Wl,--allow-shlib-undefined)
endif()

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
