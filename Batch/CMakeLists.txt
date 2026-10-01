# =============================================================================
# Batch: batchreco.exe, batchreco_detresp.exe, dumphits.exe (+ libTPORec.so)
# =============================================================================

# ======== RPATH for macOS runtime discovery ========
# GenFit::genfit2/Rave::RaveBase/CLHEP (hand-rolled IMPORTED targets from
# the superbuild in cmake/Externals.cmake, not a proper find_package()
# Config) don't give CMake's automatic RPATH computation enough to go on,
# so without this, every Batch executable links fine but fails at
# runtime with dyld errors like:
#   Library not loaded: @rpath/libgenfit2.2.3.dylib
# one per external dependency, each only surfacing once the previous one
# is fixed. Display/CMakeLists.txt hit the exact same GenFit issue first
# and fixed it the same way (see its own "RPATH for macOS Runtime
# Discovery" section) -- this mirrors that, plus Rave/CLHEP, which
# batchreco.exe also links against transitively via GenFit::genfit2 and
# BatchCore.
# CMAKE_BUILD_WITH_INSTALL_RPATH TRUE means CMAKE_INSTALL_RPATH below
# REPLACES CMake's automatic per-target RPATH computation for the
# build-tree binaries, not add to it -- so ROOT_LIBRARY_DIR has to be
# listed explicitly too, even though ROOT::Core's own rpath currently
# gets added automatically without this (found via find_package(ROOT)'s
# proper exported Config, unlike GenFit/Rave/CLHEP below) -- otherwise
# turning this on for GenFit/Rave/CLHEP would silently break that.
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
list(APPEND CMAKE_INSTALL_RPATH
  "${ROOT_LIBRARY_DIR}"
  "${GENFIT_INSTALL_DIR}/lib"
  "${RAVE_INSTALL_DIR}/lib"
  "${CLHEP_INSTALL_DIR}/lib"
)

root_generate_dictionary(G__BatchDict
  ${FASER_COREUTILS_INCLUDE_DIR}/TPORecoEvent.hh ${FASER_COREUTILS_INCLUDE_DIR}/TTauSearch.hh ${FASER_COREUTILS_INCLUDE_DIR}/TParticleGun.hh ${FASER_COREUTILS_INCLUDE_DIR}/TPSCluster.hh ${FASER_COREUTILS_INCLUDE_DIR}/TTKTrack.hh ${FASER_COREUTILS_INCLUDE_DIR}/TMuonSpect.hh
  LINKDEF LinkDef.h
  OPTIONS -I${FASER_COREUTILS_INCLUDE_DIR} -I${GENFIT_INSTALL_DIR}/include
  MODULE BatchDict
)

# All four Batch targets need the same CoreUtils objects + the same
# dictionary. Compiling them once into an object library - instead of
# listing them directly as sources in four separate targets - avoids a real
# CMake/Makefiles-generator pitfall: multiple targets that all directly list
# the same generated dictionary source can end up invoking rootcling
# concurrently under a parallel build, racing on the same output file.
add_library(BatchCore OBJECT ${FASER_COREUTILS_SOURCES} G__BatchDict.cxx)
target_include_directories(BatchCore PUBLIC ${FASER_COREUTILS_INCLUDE_DIR} ${GENFIT_INSTALL_DIR}/include)
target_link_libraries(BatchCore PUBLIC
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Physics ROOT::EG ROOT::Geom
  GenFit::genfit2
)
target_compile_options(BatchCore PUBLIC -O3)

function(faser_batch_executable name source)
  add_executable(${name} ${source})
  target_link_libraries(${name} PRIVATE BatchCore)
  faser_copy_dict_pcm(${name} BatchDict)
endfunction()

faser_batch_executable(batchreco.exe          BatchReco.cc)
faser_batch_executable(dumphits.exe           DumpHits.cc)
faser_batch_executable(batchreco_detresp.exe  BatchReco_DetResp.cc)

# Target name is scoped (CMake target names are global) but the emitted
# file is still Batch/libTPORec.so, as before.
add_library(Batch_TPORec SHARED)
target_link_libraries(Batch_TPORec PRIVATE BatchCore)
set_target_properties(Batch_TPORec PROPERTIES
  OUTPUT_NAME TPORec
  PREFIX "lib"
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
