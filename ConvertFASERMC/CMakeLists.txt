# =============================================================================
# ConvertFASERMC: Convert.exe
# Converts official FASER MC files into FASERCAL PO files. Needs Pythia8,
# but not GenFit/CLHEP/Rave.
# =============================================================================
root_generate_dictionary(G__ConvertFASERMCDict
  ${FASER_COREUTILS_INCLUDE_DIR}/TPOEvent.hh
  LINKDEF LinkDef.h
  OPTIONS -I${FASER_COREUTILS_INCLUDE_DIR}
  MODULE ConvertFASERMCDict
)

add_executable(Convert.exe Convert.cc ${FASER_COREUTILS_INCLUDE_DIR}/TPOEvent.cc G__ConvertFASERMCDict.cxx)
target_include_directories(Convert.exe PRIVATE ${FASER_COREUTILS_INCLUDE_DIR})
target_compile_definitions(Convert.exe PRIVATE _INCLUDE_PYTHIA_)
target_link_libraries(Convert.exe PRIVATE
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Physics ROOT::EG ROOT::EGPythia8 ROOT::Geom
  Pythia8::pythia8
)
