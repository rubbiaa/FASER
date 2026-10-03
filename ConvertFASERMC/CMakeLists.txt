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

# Some CVMFS/LCG-distributed ROOT releases (observed with ROOT 6.40.04's
# x86_64-almalinux9.8-gcc115-opt build on lxplus) ship a libEGPythia8.so
# that was itself linked against a *shared* Pythia8 build (e.g.
# libpythia8-8.3.17.so) that doesn't exist on this machine - we vendor our
# own, separately-versioned, *static* Pythia8 (Pythia8::pythia8, built from
# source above) instead of that LCG release. The linker reports this as
# "libpythia8-8.3.17.so, needed by .../libEGPythia8.so, not found" followed
# by undefined references to Pythia8 symbols (Settings::mode, HadronLevel::
# moreDecays, ...) *as seen from inside libEGPythia8.so* - i.e. ROOT's own
# internal use of Pythia8, which Convert.exe never calls into. The actual
# Pythia8 API calls Convert.exe itself makes are satisfied by our own
# Pythia8::pythia8 archive right above. Same root cause and same fix as
# faserps's ROOT/Geant4 clash in FASERG4/CMakeLists.txt - see the comment
# there for the general pattern. Linux (GNU ld) only; ld64 on macOS doesn't
# understand this flag and the Mac build doesn't hit this (single
# self-built ROOT with no separate Pythia8 release underneath it).
if(UNIX AND NOT APPLE)
  target_link_options(Convert.exe PRIVATE -Wl,--allow-shlib-undefined)
endif()

faser_copy_dict_pcm(Convert.exe ConvertFASERMCDict)
