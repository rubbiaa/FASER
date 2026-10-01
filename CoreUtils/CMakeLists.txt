# =============================================================================
# CoreUtils
#
# CoreUtils is not an independent library that other packages link against:
# each downstream package (Analysis, Batch, EvDisplay, FileMask, ...) has
# always compiled these sources directly into its own executable together
# with its own ROOT dictionary, because each package's LinkDef.h selects a
# different subset of CoreUtils classes for I/O. This file preserves that
# behaviour and just exposes the source/include list as CMake variables so
# every package's CMakeLists.txt can reuse them instead of re-deriving a
# `wildcard ../CoreUtils/*.cc` by hand, as the old Makefiles did.
#
# It also builds libCoreUtilsDict, a standalone dictionary library covering
# every CoreUtils class (as `make dict` used to, from the top-level
# Makefile) - CoreUtils/dumpReco.C loads it with gSystem->Load(...), and
# so does Tests/regression/summarize_output.py via PyROOT.
# =============================================================================

# ======== RPATH for macOS runtime discovery ========
# Same issue (and same fix) as Batch/CMakeLists.txt's "RPATH for macOS
# runtime discovery" section, which has the full explanation: GenFit's
# IMPORTED target doesn't give CMake's automatic RPATH computation enough
# to go on, so without this, libCoreUtilsDict.so links fine but fails to
# load at runtime (dyld: "Library not loaded: @rpath/libgenfit2...") --
# whether loaded by a FASER executable or, as here, by `python3` itself
# via PyROOT's gSystem.Load(), which has no FASER-aware RPATH of its own
# to fall back on, so this dylib must carry a correct one itself.
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
list(APPEND CMAKE_INSTALL_RPATH
  "${ROOT_LIBRARY_DIR}"
  "${GENFIT_INSTALL_DIR}/lib"
  "${RAVE_INSTALL_DIR}/lib"
  "${CLHEP_INSTALL_DIR}/lib"
)

file(GLOB FASER_COREUTILS_SOURCES CONFIGURE_DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")
set(FASER_COREUTILS_SOURCES ${FASER_COREUTILS_SOURCES} PARENT_SCOPE)
set(FASER_COREUTILS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} PARENT_SCOPE)

root_generate_dictionary(G__CoreUtilsDict
  TPORecoEvent.hh
  TTauSearch.hh
  TParticleGun.hh
  TPSCluster.hh
  TTKTrack.hh
  TMuTrack.hh
  TMuonSpect.hh
  LINKDEF LinkDef.h
  OPTIONS -I${GENFIT_INSTALL_DIR}/include
  MODULE CoreUtilsDict
)

add_library(CoreUtilsDict SHARED ${FASER_COREUTILS_SOURCES} G__CoreUtilsDict.cxx)
target_include_directories(CoreUtilsDict PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${GENFIT_INSTALL_DIR}/include
)
target_link_libraries(CoreUtilsDict PRIVATE
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Physics ROOT::EG ROOT::Geom
  GenFit::genfit2
)
set_target_properties(CoreUtilsDict PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
# PyROOT finds dictionary metadata next to the loaded library.
faser_copy_dict_pcm(CoreUtilsDict CoreUtilsDict)
