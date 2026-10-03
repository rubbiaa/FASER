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

# Some CVMFS/LCG-distributed ROOT releases (observed with ROOT 6.40.04's
# x86_64-almalinux9.8-gcc115-opt build on lxplus) ship a libEGPythia8.so
# that was itself linked against a *shared* Pythia8 build (e.g.
# libpythia8-8.3.17.so) that doesn't exist on this machine - we vendor our
# own, separately-versioned, *static* Pythia8 (Pythia8::pythia8, built from
# source above) instead of that LCG release. The linker reports this as
# "libpythia8-8.3.17.so, needed by .../libEGPythia8.so, not found" followed
# by undefined references to Pythia8 symbols (Settings::mode, HadronLevel::
# moreDecays, ...) *as seen from inside libEGPythia8.so* - i.e. ROOT's own
# internal use of Pythia8, which ConvertGENIE.exe/countinteractions_tpoevent.exe
# never call into. The actual Pythia8 API calls this code makes are
# satisfied by our own Pythia8::pythia8 archive right above. PUBLIC (not
# PRIVATE) because ConvertGENIEDictObj is an OBJECT library - it is never
# itself linked, so a PRIVATE link option here would simply never take
# effect; PUBLIC makes it a usage requirement of whatever executable links
# this object library (ConvertGENIE.exe, countinteractions_tpoevent.exe),
# same as the Pythia8::pythia8 link above. Same root cause and same fix as
# faserps's ROOT/Geant4 clash in FASERG4/CMakeLists.txt - see the comment
# there for the general pattern. Linux (GNU ld) only; ld64 on macOS doesn't
# understand this flag and the Mac build doesn't hit this (single
# self-built ROOT with no separate Pythia8 release underneath it).
if(UNIX AND NOT APPLE)
  target_link_options(ConvertGENIEDictObj PUBLIC -Wl,--allow-shlib-undefined)
endif()

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
