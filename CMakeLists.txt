# =============================================================================
# FASERCAL / FASER build system
#
# Modernized from the original top-level Makefile. Build like any CMake
# project:
#
#     cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
#     cmake --build build -j
#
# The first build will also fetch and compile the bundled third-party
# dependencies (CLHEP, Rave, GenFit, googletest, Pythia8) unless you turn
# FASER_BUILD_EXTERNALS off - see cmake/Externals.cmake for the individual
# FASER_BUILD_<NAME> switches and *_ROOT overrides.
#
# ROOT and Geant4 are expected to already be installed on your machine
# (Homebrew, a module system, CVMFS, ...) - source your usual setup.sh /
# mac_setup.sh / lxplus_setup.sh (or just have ROOTSYS / thisroot.sh /
# geant4.sh sourced) before running cmake so find_package() can locate them.
# =============================================================================
cmake_minimum_required(VERSION 3.20)
project(FASER LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

# All package executables land in one place, like the old per-package .exe
# files did in their own directories, but now discoverable in one spot.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

include(ProcessorCount)
ProcessorCount(FASER_BUILD_PARALLEL_JOBS)
if(FASER_BUILD_PARALLEL_JOBS EQUAL 0)
  set(FASER_BUILD_PARALLEL_JOBS 4)
endif()

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# -----------------------------------------------------------------------------
# Top-level switches for which parts of the tree to build
# -----------------------------------------------------------------------------
option(FASER_BUILD_EXTERNALS         "Build bundled third-party deps (CLHEP/Rave/GenFit/googletest/Pythia8) from source" ON)
option(FASER_BUILD_GEANT4_PACKAGES   "Build the Geant4 simulation packages (FASERG4, FASERCalProtoG4)" ON)
option(FASER_BUILD_DISPLAY           "Build the GenFit-based event display / muon spectrometer package (Display)" ON)
option(FASER_BUILD_TAUSEARCH         "Build the generator-level tau search package (TauSearch)" ON)

find_package(ROOT REQUIRED COMPONENTS
  Core RIO Tree Hist Physics EG Geom EGPythia8 RooFit RooFitCore
  Eve GenVector Gpad Graf RGL Rint Gui GuiHtml Thread Gdml Matrix
)
include(${ROOT_USE_FILE})

if(FASER_BUILD_EXTERNALS)
  include(Externals)
else()
  # Fall back to system/module-provided copies. Point these at existing
  # installs, e.g. -DCLHEP_ROOT=/opt/CLHEP -DRAVE_ROOT=... -DGENFIT_ROOT=...
  # -DPYTHIA8_ROOT=$PYTHIA8
  set(FASER_BUILD_CLHEP OFF)
  set(FASER_BUILD_RAVE OFF)
  set(FASER_BUILD_GENFIT OFF)
  set(FASER_BUILD_GOOGLETEST OFF)
  set(FASER_BUILD_PYTHIA8 OFF)
  include(Externals)
endif()

# -----------------------------------------------------------------------------
# FASER packages
# -----------------------------------------------------------------------------
add_subdirectory(CoreUtils)
add_subdirectory(ConvertFASERMC)
add_subdirectory(ConvertGENIE)
add_subdirectory(Analysis)
add_subdirectory(Batch)
add_subdirectory(EvDisplay)
add_subdirectory(FileMask)

if(FASER_BUILD_TAUSEARCH)
  add_subdirectory(TauSearch)
endif()

if(FASER_BUILD_GEANT4_PACKAGES)
  add_subdirectory(FASERG4)
  add_subdirectory(FASERCalProtoG4)
endif()

if(FASER_BUILD_DISPLAY)
  add_subdirectory(Display)
endif()

message(STATUS "")
message(STATUS "FASER configuration summary")
message(STATUS "  Build type              : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Build externals         : ${FASER_BUILD_EXTERNALS}")
message(STATUS "  Build Geant4 packages    : ${FASER_BUILD_GEANT4_PACKAGES}")
message(STATUS "  Build Display            : ${FASER_BUILD_DISPLAY}")
message(STATUS "  Build TauSearch          : ${FASER_BUILD_TAUSEARCH}")
message(STATUS "  Binaries will be in      : ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message(STATUS "")
