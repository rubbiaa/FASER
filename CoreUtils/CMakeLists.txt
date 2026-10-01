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
# Makefile) - CoreUtils/dumpReco.C loads it with gSystem->Load(...).
# =============================================================================
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
