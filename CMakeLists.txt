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
# (Homebrew, a module system, CVMFS, ...) - source setup.sh (it
# auto-detects which known site you're on - Mac, Ubuntu, lxplus - and
# sets ROOT/Geant4/Pythia8 up accordingly; or just have ROOTSYS /
# thisroot.sh / geant4.sh sourced yourself) before running cmake so
# find_package() can locate them.
# =============================================================================
cmake_minimum_required(VERSION 3.20)
project(FASER LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Several modules (Analysis, Batch, ...) compile their CoreUtils sources
# once into an OBJECT library that's then linked into BOTH a plain
# executable AND a small SHARED library for interactive ROOT/PyROOT use
# (e.g. AnalysisCore -> AnalyReco.exe + Analysis_TPORec/libTPORec.so).
# OBJECT libraries don't get -fPIC by default, only SHARED ones do - so
# without this, linking the .so half fails on Linux (GNU ld enforces
# this strictly: "relocation ... can not be used when making a shared
# object; recompile with -fPIC"), even though it happened to link fine
# on macOS. Setting this globally means every target compiles as
# position-independent from the start, so this whole class of bug can't
# recur as new OBJECT-library-into-SHARED-library targets get added.
# FASERG4/FASERCalProtoG4 already set this locally for themselves - this
# just makes it universal.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

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
# Environment sanity check
# -----------------------------------------------------------------------------
# A missing or stale setup.sh environment has
# repeatedly caused confusing failures well downstream of here - a cryptic
# Rave/CLHEP autotools error, or a plain find_package(Geant4) failure
# inside FASERG4/CMakeLists.txt with no further context. Report what CMake
# actually sees, up front, before anything else runs, rather than making
# someone reconstruct it from a failure several minutes into the build.
#
# GEANT4_INSTALL itself is optional (see cmake/Externals.cmake's smart
# CLHEP-reuse default) - CLHEP will just be built from source instead of
# reusing Geant4's bundled copy if it's unset, which is slower but not
# wrong - so this only warns, it never fails the configure. A genuinely
# missing Geant4 install still fails loudly and separately, on its own,
# the normal CMake way, via the find_package(Geant4 REQUIRED) calls in
# FASERG4/FASERCalProtoG4.
message(STATUS "FASER environment check:")
if(DEFINED ENV{GEANT4_INSTALL})
  if(EXISTS "$ENV{GEANT4_INSTALL}")
    message(STATUS "  GEANT4_INSTALL = $ENV{GEANT4_INSTALL}")
  else()
    message(WARNING "GEANT4_INSTALL is set to '$ENV{GEANT4_INSTALL}', but that "
      "path does not exist. This is very likely a stale environment - re-source "
      "your setup script (source setup.sh) before running cmake.")
  endif()
else()
  message(STATUS "  GEANT4_INSTALL is not set - source your setup script first "
    "(or pass -DCLHEP_ROOT=... yourself); otherwise FASER will build CLHEP "
    "from source unnecessarily instead of reusing Geant4's bundled copy")
endif()
if(CMAKE_PREFIX_PATH)
  message(STATUS "  CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}")
endif()

# -----------------------------------------------------------------------------
# Top-level switches for which parts of the tree to build
# -----------------------------------------------------------------------------
option(FASER_BUILD_EXTERNALS         "Build bundled third-party deps (CLHEP/Rave/GenFit/googletest/Pythia8) from source" ON)
option(FASER_BUILD_GEANT4_PACKAGES   "Build the Geant4 simulation packages (FASERG4, FASERCalProtoG4)" ON)
option(FASER_BUILD_DISPLAY           "Build the GenFit-based event display / muon spectrometer package (Display)" ON)
option(FASER_BUILD_TAUSEARCH         "Build the generator-level tau search package (TauSearch)" ON)
option(FASER_BUILD_TESTS             "Build FASER's own regression tests (Tests/, e.g. the muon-spectrometer field consistency check)" ON)

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
# Some CVMFS/LCG-distributed ROOT releases (observed with ROOT 6.40.04's
# x86_64-almalinux9.8-gcc115-opt build on lxplus) export their ROOT::EGPythia8
# imported target with an INTERFACE_LINK_LIBRARIES entry that is a hardcoded
# absolute path to whatever Pythia8 install happened to exist on the machine
# that built that ROOT release (e.g. "/usr/lib64/libpythia8.so") - a system
# RPM path that doesn't exist here. Because it's a plain absolute file path
# rather than an imported target, CMake has no choice but to treat it as a
# literal Makefile prerequisite, and the build fails with
# "No rule to make target '/usr/lib64/libpythia8.so'" for anything that
# links ROOT::EGPythia8 (ConvertFASERMC's Convert.exe, ConvertGENIE's
# ConvertGENIE.exe/countinteractions_tpoevent.exe), even though our own
# Pythia8 (Pythia8::pythia8, from this same superbuild) is perfectly capable
# of satisfying whatever symbols ROOT::EGPythia8 actually needs. Scrub any
# such stale absolute path from every ROOT component target we use, and swap
# in Pythia8::pythia8 in its place so the dependency is still satisfied.
foreach(_faser_root_pythia8_dep_target IN ITEMS
    Core RIO Tree Hist Physics EG Geom EGPythia8 RooFit RooFitCore
    Eve GenVector Gpad Graf RGL Rint Gui GuiHtml Thread Gdml Matrix)
  if(TARGET ROOT::${_faser_root_pythia8_dep_target})
    get_target_property(_faser_root_pythia8_dep_libs
      ROOT::${_faser_root_pythia8_dep_target} INTERFACE_LINK_LIBRARIES)
    if(_faser_root_pythia8_dep_libs)
      set(_faser_root_pythia8_dep_changed FALSE)
      set(_faser_root_pythia8_dep_new "")
      foreach(_faser_root_pythia8_dep_lib IN LISTS _faser_root_pythia8_dep_libs)
        if(_faser_root_pythia8_dep_lib MATCHES "/libpythia8\.(a|so)(\.[0-9.]+)?$"
           AND NOT _faser_root_pythia8_dep_lib STREQUAL "Pythia8::pythia8"
           AND NOT EXISTS "${_faser_root_pythia8_dep_lib}")
          message(STATUS "FASER: dropping stale ${_faser_root_pythia8_dep_lib} "
                         "from ROOT::${_faser_root_pythia8_dep_target}'s link "
                         "interface (doesn't exist here) - using "
                         "Pythia8::pythia8 instead")
          list(APPEND _faser_root_pythia8_dep_new "Pythia8::pythia8")
          set(_faser_root_pythia8_dep_changed TRUE)
        else()
          list(APPEND _faser_root_pythia8_dep_new "${_faser_root_pythia8_dep_lib}")
        endif()
      endforeach()
      if(_faser_root_pythia8_dep_changed)
        list(REMOVE_DUPLICATES _faser_root_pythia8_dep_new)
        set_target_properties(ROOT::${_faser_root_pythia8_dep_target} PROPERTIES
          INTERFACE_LINK_LIBRARIES "${_faser_root_pythia8_dep_new}")
      endif()
      unset(_faser_root_pythia8_dep_new)
      unset(_faser_root_pythia8_dep_changed)
    endif()
    unset(_faser_root_pythia8_dep_libs)
  endif()
endforeach()
unset(_faser_root_pythia8_dep_target)

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

if(FASER_BUILD_TESTS AND FASER_BUILD_GEANT4_PACKAGES)
  # The only FASER-authored test right now (MuonSpectrometerFieldTest)
  # cross-checks FASERG4's field against GenFit's, so it needs the
  # Geant4 packages built too - see Tests/CMakeLists.txt.
  enable_testing()
  add_subdirectory(Tests)
endif()

message(STATUS "")
message(STATUS "FASER configuration summary")
message(STATUS "  Build type              : ${CMAKE_BUILD_TYPE}")
message(STATUS "  Build externals         : ${FASER_BUILD_EXTERNALS}")
message(STATUS "  Build Geant4 packages    : ${FASER_BUILD_GEANT4_PACKAGES}")
message(STATUS "  Build Display            : ${FASER_BUILD_DISPLAY}")
message(STATUS "  Build TauSearch          : ${FASER_BUILD_TAUSEARCH}")
message(STATUS "  Build tests              : ${FASER_BUILD_TESTS}")
message(STATUS "  Binaries will be in      : ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message(STATUS "")
