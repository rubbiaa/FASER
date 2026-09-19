# =============================================================================
# FileMask: FileMask.exe
# =============================================================================
root_generate_dictionary(G__FileMaskDict
  ${FASER_COREUTILS_INCLUDE_DIR}/TPORecoEvent.hh ${FASER_COREUTILS_INCLUDE_DIR}/TTauSearch.hh ${FASER_COREUTILS_INCLUDE_DIR}/TParticleGun.hh ${FASER_COREUTILS_INCLUDE_DIR}/TPSCluster.hh ${FASER_COREUTILS_INCLUDE_DIR}/TTKTrack.hh ${FASER_COREUTILS_INCLUDE_DIR}/TPSTrack.hh
  LINKDEF LinkDef.h
  OPTIONS -I${FASER_COREUTILS_INCLUDE_DIR} -I${GENFIT_INSTALL_DIR}/include
  MODULE FileMaskDict
)

add_executable(FileMask.exe FileMask.cc ${FASER_COREUTILS_SOURCES} G__FileMaskDict.cxx)
target_include_directories(FileMask.exe PRIVATE ${FASER_COREUTILS_INCLUDE_DIR} ${GENFIT_INSTALL_DIR}/include)
target_link_libraries(FileMask.exe PRIVATE
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Physics ROOT::EG ROOT::Geom
  GenFit::genfit2
)
