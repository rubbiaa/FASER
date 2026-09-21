# =============================================================================
# EvDisplay: evDisplay.exe
# =============================================================================
root_generate_dictionary(G__EvDisplayDict
  MyMainFrame.h ${FASER_COREUTILS_INCLUDE_DIR}/TPORecoEvent.hh ${FASER_COREUTILS_INCLUDE_DIR}/TTauSearch.hh ${FASER_COREUTILS_INCLUDE_DIR}/TParticleGun.hh ${FASER_COREUTILS_INCLUDE_DIR}/TMuonSpect.hh ${FASER_COREUTILS_INCLUDE_DIR}/TPSCluster.hh ${FASER_COREUTILS_INCLUDE_DIR}/TTKTrack.hh ${FASER_COREUTILS_INCLUDE_DIR}/TPSTrack.hh
  LINKDEF LinkDef.h
  OPTIONS -I${FASER_COREUTILS_INCLUDE_DIR} -I${GENFIT_INSTALL_DIR}/include
  MODULE EvDisplayDict
)

add_executable(evDisplay.exe evDisplay.C MyMainFrame.cc ${FASER_COREUTILS_SOURCES} G__EvDisplayDict.cxx)
target_include_directories(evDisplay.exe PRIVATE ${CMAKE_CURRENT_SOURCE_DIR} ${FASER_COREUTILS_INCLUDE_DIR} ${GENFIT_INSTALL_DIR}/include)
target_link_libraries(evDisplay.exe PRIVATE
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Physics ROOT::EG ROOT::Geom
  ROOT::Gui ROOT::GuiHtml
  GenFit::genfit2
)
faser_copy_dict_pcm(evDisplay.exe EvDisplayDict)
