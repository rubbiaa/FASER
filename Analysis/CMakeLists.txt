# =============================================================================
# Analysis: AnalyReco.exe (+ libTPORec.so) and plots.exe
# =============================================================================
root_generate_dictionary(G__AnalysisDict
  ${FASER_COREUTILS_INCLUDE_DIR}/TPORecoEvent.hh ${FASER_COREUTILS_INCLUDE_DIR}/TTauSearch.hh ${FASER_COREUTILS_INCLUDE_DIR}/TParticleGun.hh ${FASER_COREUTILS_INCLUDE_DIR}/TPSCluster.hh ${FASER_COREUTILS_INCLUDE_DIR}/TTKTrack.hh ${FASER_COREUTILS_INCLUDE_DIR}/TPSTrack.hh
  LINKDEF LinkDef.h
  OPTIONS -I${FASER_COREUTILS_INCLUDE_DIR} -I${GENFIT_INSTALL_DIR}/include
  MODULE AnalysisDict
)

# AnalyReco.exe and libTPORec.so both need the same CoreUtils objects + the
# same dictionary. Compiling them once into an object library - instead of
# listing them directly as sources in two separate targets - avoids a real
# CMake/Makefiles-generator pitfall: two targets that both directly list the
# same generated dictionary source can end up invoking rootcling twice
# concurrently under a parallel build, racing on the same output file.
add_library(AnalysisCore OBJECT ${FASER_COREUTILS_SOURCES} G__AnalysisDict.cxx)
target_include_directories(AnalysisCore PUBLIC ${FASER_COREUTILS_INCLUDE_DIR} ${GENFIT_INSTALL_DIR}/include)
target_link_libraries(AnalysisCore PUBLIC
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Physics ROOT::EG ROOT::Geom
  GenFit::genfit2
)

add_executable(AnalyReco.exe AnalyReco.cc)
target_link_libraries(AnalyReco.exe PRIVATE AnalysisCore)

# The old Makefile also produced a small shared library, from the same
# dictionary + CoreUtils objects, alongside AnalyReco.exe (for interactive
# ROOT/PyROOT use). Target name is scoped (CMake target names are global)
# but the emitted file is still Analysis/libTPORec.so as before.
add_library(Analysis_TPORec SHARED)
target_link_libraries(Analysis_TPORec PRIVATE AnalysisCore)
set_target_properties(Analysis_TPORec PROPERTIES
  OUTPUT_NAME TPORec
  PREFIX "lib"
  LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)

# plots.exe draws canvases/legends (TCanvas, TLegend, ...), which live in
# ROOT's Gpad/Graf libraries - the old Makefile linked with the broader
# `root-config --glibs`, which pulls those in automatically.
add_executable(plots.exe plots.C)
target_link_libraries(plots.exe PRIVATE
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::Graf ROOT::RooFit ROOT::RooFitCore
)
