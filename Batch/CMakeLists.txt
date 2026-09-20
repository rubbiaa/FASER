# =============================================================================
# Batch: batchreco.exe, batchreco_detresp.exe, dumphits.exe (+ libTPORec.so)
# =============================================================================
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
