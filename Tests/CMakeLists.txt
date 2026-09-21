# =============================================================================
# Tests
#
# FASER's own regression tests, as opposed to the bundled externals'
# (GenFit's, CLHEP's, ...) internal test suites, which FASER never builds
# or runs on its own account. Add more tests as their own
# add_executable()/add_test() pair below.
# =============================================================================

# Cross-checks that the Geant4 simulation's magnetic field
# (MuonMagneticField, FASERG4/src/MuonDetMagneticField.cc) and the
# GenFit reconstruction's magnetic field (GenMagneticField, CoreUtils/
# GenMagneticField.hh) agree - both are thin adapters around the one
# shared field model in CoreUtils/MuonSpectrometerField.hh, but this
# test is what actually proves the two adapters wire it up identically,
# rather than just trusting that by inspection. Run it after touching
# either adapter or the shared field model itself.
find_package(Geant4 REQUIRED ui_all vis_all gdml)
include(${Geant4_USE_FILE})

add_executable(MuonSpectrometerFieldTest
  MuonSpectrometerFieldTest.cc
  ${FASER_COREUTILS_INCLUDE_DIR}/MuonSpectrometerField.cc
  ${FASER_COREUTILS_INCLUDE_DIR}/GenMagneticField.cc
  ${CMAKE_SOURCE_DIR}/FASERG4/src/MuonDetMagneticField.cc
)
target_include_directories(MuonSpectrometerFieldTest PRIVATE
  ${FASER_COREUTILS_INCLUDE_DIR}
  ${CMAKE_SOURCE_DIR}/FASERG4/include
  ${GENFIT_INSTALL_DIR}/include
  ${ROOT_INCLUDE_DIRS}
)
target_link_libraries(MuonSpectrometerFieldTest PRIVATE
  ${Geant4_LIBRARIES}
  ROOT::Core ROOT::RIO ROOT::Physics ROOT::Geom
  GenFit::genfit2
  GTest::gtest GTest::gtest_main
)

if(UNIX AND NOT APPLE)
  # Same Geant4/ROOT ABI clash as FASERG4/FASERCalProtoG4 - this test
  # links Geant4 too, so it needs the same workaround on lxplus. See
  # FASERG4/CMakeLists.txt for the full explanation.
  target_link_options(MuonSpectrometerFieldTest PRIVATE -Wl,--allow-shlib-undefined)
endif()

add_test(NAME MuonSpectrometerFieldTest COMMAND MuonSpectrometerFieldTest)

# Unit tests for FASER::ProbeMagnetSlit (CoreUtils/MagnetGeometryProbe.hh)
# - the "full paranoia" geometry<->field cross-check: builds small ROOT
# TGeo shapes by hand (no Geant4/GDML needed for this one) to prove the
# probe both agrees with the real MDT magnet's numbers AND actually
# detects a genuine mismatch when the block and slit are deliberately
# built inconsistent, rather than just trusting matching literals like
# GenMagneticField.cc's SetupMDTMagneticField used to.
add_executable(MagnetGeometryProbeTest
  MagnetGeometryProbeTest.cc
  ${FASER_COREUTILS_INCLUDE_DIR}/MagnetGeometryProbe.cc
)
target_include_directories(MagnetGeometryProbeTest PRIVATE
  ${FASER_COREUTILS_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIRS}
)
target_link_libraries(MagnetGeometryProbeTest PRIVATE
  ROOT::Core ROOT::RIO ROOT::Geom
  GTest::gtest GTest::gtest_main
)

if(UNIX AND NOT APPLE)
  target_link_options(MagnetGeometryProbeTest PRIVATE -Wl,--allow-shlib-undefined)
endif()

add_test(NAME MagnetGeometryProbeTest COMMAND MagnetGeometryProbeTest)
