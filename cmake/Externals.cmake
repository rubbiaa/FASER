# =============================================================================
# FASER third-party dependency superbuild
#
# This mirrors, step for step, what the old top-level Makefile did by hand
# with wget/git-clone/configure recipes, but as a proper CMake ExternalProject
# graph so `cmake --build .` builds everything in dependency order:
#
#     CLHEP  ->  Rave  ->  GenFit  (+ googletest, Linux only)
#     Pythia8 (independent)
#
# ROOT and Geant4 are deliberately NOT built here. They are expected to
# already be installed (Homebrew, a module system, CVMFS, ...) and are
# located with find_package() in the top-level CMakeLists.txt, exactly as
# FASERG4/CMakeLists.txt and Display/CMakeLists.txt already did before this
# modernization. Building ROOT itself from source is a multi-hour build
# that nobody on the team actually relied on (setup.sh / mac_setup.sh /
# lxplus_setup.sh all point at a pre-existing ROOT install) so it is out of
# scope for this superbuild.
#
# Each dependency can be switched off individually (FASER_BUILD_<NAME>) to
# point at an already-installed copy instead - set the matching *_ROOT /
# *_DIR variable when doing so.
# =============================================================================
include(ExternalProject)
include(GNUInstallDirs)

set(FASER_EXTERNAL_INSTALL_DIR "${CMAKE_BINARY_DIR}/external-install" CACHE PATH
    "Install prefix for FASER's bundled third-party dependencies")
set(FASER_EXTERNAL_STAGE_DIR "${CMAKE_BINARY_DIR}/external" CACHE PATH
    "Scratch/build directory for FASER's bundled third-party dependencies")

if(APPLE)
  set(_faser_shlib_suffix ".dylib")
else()
  set(_faser_shlib_suffix ".so")
endif()

# -----------------------------------------------------------------------------
# CLHEP
# -----------------------------------------------------------------------------
option(FASER_BUILD_CLHEP "Build CLHEP from source (gitlab.cern.ch/CLHEP/CLHEP)" ON)
set(CLHEP_ROOT "" CACHE PATH "Pre-installed CLHEP prefix (used when FASER_BUILD_CLHEP=OFF)")

if(FASER_BUILD_CLHEP)
  set(CLHEP_INSTALL_DIR "${FASER_EXTERNAL_INSTALL_DIR}/CLHEP")

  ExternalProject_Add(clhep_external
    GIT_REPOSITORY    https://gitlab.cern.ch/CLHEP/CLHEP.git
    GIT_TAG           develop  # CLHEP's default branch (verified via the GitLab UI)
    GIT_SHALLOW       TRUE
    PREFIX            "${FASER_EXTERNAL_STAGE_DIR}/CLHEP"
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CLHEP_INSTALL_DIR}
      -DCLHEP_SINGLE_THREAD=ON
      -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    BUILD_COMMAND     ${CMAKE_COMMAND} --build <BINARY_DIR> --parallel
    INSTALL_COMMAND   ${CMAKE_COMMAND} --build <BINARY_DIR> --target install
    BUILD_BYPRODUCTS  "${CLHEP_INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/libCLHEP${_faser_shlib_suffix}"
  )
else()
  if(NOT CLHEP_ROOT)
    message(FATAL_ERROR "FASER_BUILD_CLHEP=OFF but CLHEP_ROOT was not set to a pre-installed CLHEP prefix")
  endif()
  set(CLHEP_INSTALL_DIR "${CLHEP_ROOT}")
endif()

# The include dir must exist at CMake generate time (CMake validates
# INTERFACE_INCLUDE_DIRECTORIES on imported targets eagerly) even though
# ExternalProject only populates it later, during the build step.
#
# Named FASER::CLHEP rather than the more natural CLHEP::CLHEP: on Linux,
# find_package(Geant4 ...) below (in FASERG4/FASERCalProtoG4) transitively
# calls find_dependency(CLHEP) because conda-forge's Geant4 is built
# against system/conda CLHEP, not a bundled copy. That real find_package()
# call defines its own genuine CLHEP::CLHEP imported target - if we had
# already claimed that exact name for our own from-source CLHEP build,
# CMake errors at configure time with "add_library cannot create imported
# target CLHEP::CLHEP because another target with that name already
# exists". Using our own namespaced target name sidesteps the collision
# entirely (this never showed up on macOS because the Mac Geant4 install
# here is built with its own internal CLHEP, so it never calls
# find_dependency(CLHEP) in the first place).
file(MAKE_DIRECTORY "${CLHEP_INSTALL_DIR}/include")
add_library(FASER::CLHEP SHARED IMPORTED GLOBAL)
set_target_properties(FASER::CLHEP PROPERTIES
  IMPORTED_LOCATION             "${CLHEP_INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/libCLHEP${_faser_shlib_suffix}"
  INTERFACE_INCLUDE_DIRECTORIES "${CLHEP_INSTALL_DIR}/include"
)
if(TARGET clhep_external)
  add_dependencies(FASER::CLHEP clhep_external)
endif()

# -----------------------------------------------------------------------------
# Rave (autotools; source already vendored in-tree under ./rave)
# -----------------------------------------------------------------------------
option(FASER_BUILD_RAVE "Build Rave from the vendored ./rave source tree" ON)
set(RAVE_ROOT "" CACHE PATH "Pre-installed Rave prefix (used when FASER_BUILD_RAVE=OFF)")

if(FASER_BUILD_RAVE)
  set(RAVE_INSTALL_DIR "${FASER_EXTERNAL_INSTALL_DIR}/rave")
  set(RAVE_SOURCE_DIR  "${CMAKE_SOURCE_DIR}/rave")

  if(APPLE)
    # The vendored ./rave tree ships a checked-in `configure` script but no
    # Makefile.in files anywhere (autotools was apparently never run on this
    # copy - every Makefile.am is missing its generated Makefile.in), so
    # config.status has nothing to instantiate and fails with "cannot find
    # input file: `Makefile.in'". autoreconf regenerates them first, exactly
    # like the non-Darwin branch below already does. Requires autoconf/
    # automake/libtool (e.g. `brew install autoconf automake libtool`).
    set(_rave_configure_cmd
      sh -c "autoreconf -fi ${RAVE_SOURCE_DIR} && ${RAVE_SOURCE_DIR}/configure --prefix=${RAVE_INSTALL_DIR} --disable-java --with-boost=/opt/homebrew --with-boost-libdir=/opt/homebrew/lib --with-clhep=${CLHEP_INSTALL_DIR}")
  else()
    # Matches the non-Darwin branch: autoreconf first, then configure.
    set(_rave_configure_cmd
      sh -c "autoreconf -fi ${RAVE_SOURCE_DIR} && ${RAVE_SOURCE_DIR}/configure --prefix=${RAVE_INSTALL_DIR} --disable-java --with-clhep=${CLHEP_INSTALL_DIR}")
  endif()

  # Rave's autotools build is run in-source (as the old Makefile did: `cd
  # rave && ./configure ...`), so it leaves build artifacts inside the
  # tracked ./rave directory. Run `git clean -fdx rave/` if you ever need a
  # pristine rebuild.
  ExternalProject_Add(rave_external
    DEPENDS           clhep_external
    SOURCE_DIR        ${RAVE_SOURCE_DIR}
    BUILD_IN_SOURCE   1
    CONFIGURE_COMMAND ${_rave_configure_cmd}
    BUILD_COMMAND     make CXXFLAGS=-g\ -std=c++11 LHEPINCPATH=. -j${FASER_BUILD_PARALLEL_JOBS}
    INSTALL_COMMAND   make install
    BUILD_BYPRODUCTS  "${RAVE_INSTALL_DIR}/lib/libRaveBase${_faser_shlib_suffix}"
  )
else()
  if(NOT RAVE_ROOT)
    message(FATAL_ERROR "FASER_BUILD_RAVE=OFF but RAVE_ROOT was not set to a pre-installed Rave prefix")
  endif()
  set(RAVE_INSTALL_DIR "${RAVE_ROOT}")
endif()

file(MAKE_DIRECTORY "${RAVE_INSTALL_DIR}/include")
add_library(Rave::RaveBase SHARED IMPORTED GLOBAL)
set_target_properties(Rave::RaveBase PROPERTIES
  IMPORTED_LOCATION             "${RAVE_INSTALL_DIR}/lib/libRaveBase${_faser_shlib_suffix}"
  INTERFACE_INCLUDE_DIRECTORIES "${RAVE_INSTALL_DIR}/include"
)
if(TARGET rave_external)
  add_dependencies(Rave::RaveBase rave_external)
endif()

# -----------------------------------------------------------------------------
# googletest (Linux-only build-time dependency of GenFit's CMake build)
# -----------------------------------------------------------------------------
if(NOT APPLE)
  option(FASER_BUILD_GOOGLETEST "Build googletest from source (needed by GenFit on Linux)" ON)
  set(GOOGLETEST_ROOT "" CACHE PATH "Pre-installed googletest prefix (used when FASER_BUILD_GOOGLETEST=OFF)")

  if(FASER_BUILD_GOOGLETEST)
    set(GOOGLETEST_INSTALL_DIR "${FASER_EXTERNAL_INSTALL_DIR}/googletest")

    ExternalProject_Add(googletest_external
      GIT_REPOSITORY   https://github.com/google/googletest.git
      GIT_TAG          main  # googletest's default branch
      GIT_SHALLOW      TRUE
      PREFIX           "${FASER_EXTERNAL_STAGE_DIR}/googletest"
      CMAKE_ARGS
        -DCMAKE_INSTALL_PREFIX=${GOOGLETEST_INSTALL_DIR}
      BUILD_BYPRODUCTS "${GOOGLETEST_INSTALL_DIR}/lib/libgtest.a"
                        "${GOOGLETEST_INSTALL_DIR}/lib/libgtest_main.a"
    )
  else()
    if(NOT GOOGLETEST_ROOT)
      message(FATAL_ERROR "FASER_BUILD_GOOGLETEST=OFF but GOOGLETEST_ROOT was not set")
    endif()
    set(GOOGLETEST_INSTALL_DIR "${GOOGLETEST_ROOT}")
  endif()
endif()

# -----------------------------------------------------------------------------
# GenFit
# -----------------------------------------------------------------------------
option(FASER_BUILD_GENFIT "Build GenFit from source (github.com/GenFit/GenFit)" ON)
set(GENFIT_ROOT "" CACHE PATH "Pre-installed GenFit prefix (used when FASER_BUILD_GENFIT=OFF)")

if(FASER_BUILD_GENFIT)
  set(GENFIT_INSTALL_DIR "${FASER_EXTERNAL_INSTALL_DIR}/GenFit")

  set(_genfit_cmake_args
    -DCMAKE_BUILD_TYPE=Debug
    -DCMAKE_INSTALL_PREFIX=${GENFIT_INSTALL_DIR}
    "-DRave_CFLAGS=-DRaveDllExport= -DWITH_FLAVORTAGGING -DWITH_KINEMATICS"
    -DRave_INCLUDE_DIRS=${RAVE_INSTALL_DIR}/include/
  )

  set(_genfit_deps clhep_external rave_external)

  if(APPLE)
    list(APPEND _genfit_cmake_args
      "-DRave_LDFLAGS=-L${RAVE_INSTALL_DIR}/lib/ -lRaveBase -L${CLHEP_INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/ -lCLHEP")
    # Plain parallel build; no gtest workaround needed on Darwin (matches
    # the old Makefile's Darwin branch).
    set(_genfit_build_cmd ${CMAKE_COMMAND} --build <BINARY_DIR> --parallel)
  else()
    list(APPEND _genfit_cmake_args
      -DGTEST_LIBRARY=${GOOGLETEST_INSTALL_DIR}/lib/libgtest.a
      -DGTEST_INCLUDE_DIR=${GOOGLETEST_INSTALL_DIR}/include
      -DGTEST_MAIN_LIBRARY=${GOOGLETEST_INSTALL_DIR}/lib/libgtest_main.a
      "-DRave_LDFLAGS=-Wl,-rpath-link,${RAVE_INSTALL_DIR}/lib/ -L${RAVE_INSTALL_DIR}/lib/ -lRaveBase -L${CLHEP_INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/ -lCLHEP")
    list(APPEND _genfit_deps googletest_external)
    # GenFit's own CMake build has a known issue linking its gtest binaries
    # on Linux; the old Makefile worked around it with:
    #   make CXXFLAGS="-g" -j; sh CMakeFiles/gtests.dir/link.txt; make -j
    # Reproduce that exact three-step dance as a single shell command
    # (ExternalProject_Add's BUILD_COMMAND takes one command line).
    set(_genfit_build_cmd sh -c
      "cd <BINARY_DIR> && (${CMAKE_COMMAND} --build . --parallel || true) && (sh CMakeFiles/gtests.dir/link.txt || true) && ${CMAKE_COMMAND} --build . --parallel")
  endif()

  ExternalProject_Add(genfit_external
    DEPENDS           ${_genfit_deps}
    GIT_REPOSITORY    https://github.com/GenFit/GenFit.git
    GIT_TAG           main  # GenFit's default branch
    GIT_SHALLOW       TRUE
    PREFIX            "${FASER_EXTERNAL_STAGE_DIR}/GenFit"
    # -N: skip hunks that are already applied instead of interactively
    # asking "Assume -R?" - without it, this step fails on every build
    # after the first, because GenFit tracks the `main` branch and is
    # re-fetched (see the update step above) without reverting the
    # previous run's already-applied local patch. `|| true` covers the
    # remaining case where -N still exits non-zero purely because every
    # hunk was already applied (nothing left to do is not a failure).
    #
    # GenFit also vendors its own pre-Eigen-3.3 FindEigen3.cmake, which
    # parses EIGEN_WORLD_VERSION/EIGEN_MAJOR_VERSION/EIGEN_MINOR_VERSION
    # out of Eigen's Macros.h with a regex that no longer matches modern
    # Eigen headers (Homebrew's included) - it silently resolves the
    # version to "..", which always fails the minimum-version check.
    # Deleting it makes find_package(Eigen3) fall through to Eigen's own
    # proper exported CMake config instead (Eigen3Config.cmake, shipped
    # since Eigen 3.3), which handles versioning correctly.
    PATCH_COMMAND     sh -c "patch -p0 -N -u -i ${CMAKE_SOURCE_DIR}/genfit.patch || true; rm -f cmake/FindEigen3.cmake"
    CMAKE_ARGS        ${_genfit_cmake_args}
    BUILD_COMMAND     ${_genfit_build_cmd}
    INSTALL_COMMAND   ${CMAKE_COMMAND} --build <BINARY_DIR> --target install
    BUILD_BYPRODUCTS  "${GENFIT_INSTALL_DIR}/lib/libgenfit2${_faser_shlib_suffix}"
  )
else()
  if(NOT GENFIT_ROOT)
    message(FATAL_ERROR "FASER_BUILD_GENFIT=OFF but GENFIT_ROOT was not set to a pre-installed GenFit prefix")
  endif()
  set(GENFIT_INSTALL_DIR "${GENFIT_ROOT}")
endif()

file(MAKE_DIRECTORY "${GENFIT_INSTALL_DIR}/include")
add_library(GenFit::genfit2 SHARED IMPORTED GLOBAL)
set_target_properties(GenFit::genfit2 PROPERTIES
  IMPORTED_LOCATION             "${GENFIT_INSTALL_DIR}/lib/libgenfit2${_faser_shlib_suffix}"
  INTERFACE_INCLUDE_DIRECTORIES "${GENFIT_INSTALL_DIR}/include"
  INTERFACE_LINK_LIBRARIES      "Rave::RaveBase;FASER::CLHEP"
)
if(TARGET genfit_external)
  add_dependencies(GenFit::genfit2 genfit_external)
endif()

# -----------------------------------------------------------------------------
# Pythia8
# -----------------------------------------------------------------------------
option(FASER_BUILD_PYTHIA8 "Build Pythia8 from source" ON)
set(PYTHIA8_ROOT "" CACHE PATH "Pre-installed/pre-built Pythia8 tree (used when FASER_BUILD_PYTHIA8=OFF)")

if(FASER_BUILD_PYTHIA8)
  # Same source mirror the old Makefile's `pythia8_tar` target used. Pythia8's
  # own build places headers/libs directly under the source tree - it is
  # never `make install`-ed, matching the previous behaviour.
  ExternalProject_Add(pythia8_external
    URL                        https://cernbox.cern.ch/s/lOhu8P3H0beVnb0/download
    DOWNLOAD_NAME              pythia8312.tgz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    PREFIX                     "${FASER_EXTERNAL_STAGE_DIR}/pythia8"
    BUILD_IN_SOURCE   1
    CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<SOURCE_DIR>
    BUILD_COMMAND     make -j${FASER_BUILD_PARALLEL_JOBS}
    INSTALL_COMMAND   ""
    BUILD_BYPRODUCTS  "<SOURCE_DIR>/lib/libpythia8.a"
  )
  ExternalProject_Get_Property(pythia8_external SOURCE_DIR)
  set(PYTHIA8_INSTALL_DIR "${SOURCE_DIR}")
else()
  if(NOT PYTHIA8_ROOT)
    message(FATAL_ERROR "FASER_BUILD_PYTHIA8=OFF but PYTHIA8_ROOT was not set to a Pythia8 tree with include/ and lib/")
  endif()
  set(PYTHIA8_INSTALL_DIR "${PYTHIA8_ROOT}")
endif()

file(MAKE_DIRECTORY "${PYTHIA8_INSTALL_DIR}/include")
add_library(Pythia8::pythia8 STATIC IMPORTED GLOBAL)
set_target_properties(Pythia8::pythia8 PROPERTIES
  IMPORTED_LOCATION             "${PYTHIA8_INSTALL_DIR}/lib/libpythia8.a"
  INTERFACE_INCLUDE_DIRECTORIES "${PYTHIA8_INSTALL_DIR}/include"
)
if(TARGET pythia8_external)
  add_dependencies(Pythia8::pythia8 pythia8_external)
endif()
