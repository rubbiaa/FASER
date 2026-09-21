# =============================================================================
# ConvertGENIE: ConvertGENIE.exe, CombineFluxes.exe
# (+ a handful of standalone analysis executables that used to be built by
#  hand with the mycompile.sh script)
# =============================================================================
root_generate_dictionary(G__ConvertGENIEDict
  ${FASER_COREUTILS_INCLUDE_DIR}/TPOEvent.hh
  LINKDEF LinkDef.h
  OPTIONS -I${FASER_COREUTILS_INCLUDE_DIR}
  MODULE ConvertGENIEDict
)

# ConvertGENIE.exe and countinteractions_tpoevent.exe both need this same
# dictionary. Compiling the generated G__ConvertGENIEDict.cxx once into a
# small object library - instead of listing it directly as a source in two
# separate add_executable() calls - avoids a real CMake/Makefiles-generator
# pitfall: two targets that both directly list the same generated source
# can end up with rootcling invoked twice concurrently under a parallel
# build, racing on the same output file.
add_library(ConvertGENIEDictObj OBJECT G__ConvertGENIEDict.cxx ${FASER_COREUTILS_INCLUDE_DIR}/TPOEvent.cc)
target_include_directories(ConvertGENIEDictObj PUBLIC ${FASER_COREUTILS_INCLUDE_DIR})
target_compile_definitions(ConvertGENIEDictObj PUBLIC _INCLUDE_PYTHIA_)
target_link_libraries(ConvertGENIEDictObj PUBLIC
  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Physics ROOT::EG ROOT::EGPythia8 ROOT::Geom
  Pythia8::pythia8
)

add_executable(ConvertGENIE.exe ConvertGENIE.cc)
target_link_libraries(ConvertGENIE.exe PRIVATE ConvertGENIEDictObj)
faser_copy_dict_pcm(ConvertGENIE.exe ConvertGENIEDict)

# CombineFluxes only needs ROOT, no Pythia8 (it is a plain .C file, as the
# old Makefile compiled it: `g++ ... CombineFluxes.C -o CombineFluxes.exe`).
add_executable(CombineFluxes.exe CombineFluxes.C)
target_link_libraries(CombineFluxes.exe PRIVATE ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::EG)

# Small standalone utilities previously built ad hoc via mycompile.sh. They
# draw canvases/legends/lines, which live in ROOT's Gpad/Graf libraries (the
# original script used the broader `root-config --libs`, which pulls those
# in automatically).
add_executable(countInteractions countInteractions.cc)
target_link_libraries(countInteractions PRIVATE ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::Graf)

add_executable(plotInteractionVertices plotInteractionVertices.cc)
target_link_libraries(plotInteractionVertices PRIVATE ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::Graf)

add_executable(etaCoverage etaCoverage.cc)
target_link_libraries(etaCoverage PRIVATE ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::Graf)

add_executable(etaCoverage_geometry etaCoverage_geometry.cc)
target_link_libraries(etaCoverage_geometry PRIVATE ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::Gpad ROOT::Graf)

add_executable(countinteractions_tpoevent.exe countinteractions_tpoevent.cc)
target_link_libraries(countinteractions_tpoevent.exe PRIVATE ConvertGENIEDictObj)
faser_copy_dict_pcm(countinteractions_tpoevent.exe ConvertGENIEDict)
