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
# GNU Make jobserver inheritance for nested external builds
# -----------------------------------------------------------------------------
# Rave and Pythia8 are both built with a raw, hand-invoked `make` inside
# ExternalProject_Add's BUILD_COMMAND, rather than via `cmake --build`. When
# the top-level generator is itself "Unix Makefiles", the outer
# `cmake --build . --parallel N` / `make -jN` already runs as a GNU Make
# jobserver (a shared pipe of N tokens that coordinates concurrency across
# the whole build tree). Passing an explicit `-j${FASER_BUILD_PARALLEL_JOBS}`
# to Rave/Pythia8's nested `make` - as this file used to do - ignores that
# shared pool and starts a second, uncoordinated one instead, which is why
# builds print:
#     warning: -jN forced in submake: resetting jobserver mode.
#     warning: jobserver unavailable: using -j1.  Add '+' to parent make rule.
#
# GNU Make only preserves the jobserver's file descriptors across fork/exec
# for a recipe line whose *unexpanded* text contains the literal substring
# "$(MAKE)" (its heuristic for "this line recursively invokes make"). CMake's
# own variable syntax is ${...}, not $(...), so writing the literal token
# $(MAKE) in a BUILD_COMMAND passes through CMake substitution untouched and
# survives into the generated Makefile recipe for GNU Make to recognize -
# letting Rave/Pythia8's nested make share the outer jobserver's token pool
# instead of racing it with a fixed -jN of its own.
#
# This only makes sense for the "Unix Makefiles" generator (Ninja has no
# jobserver concept), and only for BUILD_COMMANDs that invoke `make`
# directly - CLHEP and GenFit both build via
# `${CMAKE_COMMAND} --build <dir> --parallel N`, a different invocation
# shape where the $(MAKE)-token trick doesn't apply, and neither has ever
# been observed emitting this warning, so both are deliberately left alone
# here rather than reworked speculatively.
#
# CI explicitly opts back out of this (-DFASER_EXTERNALS_USE_JOBSERVER=OFF
# in .github/workflows/build.yml) to keep its already-tested, deliberately
# serial -DFASER_BUILD_PARALLEL_JOBS=1 OOM-avoidance behavior byte-for-byte
# unchanged rather than silently switching a constrained CI runner over to
# jobserver-shared concurrency without being able to re-verify it there.
option(FASER_EXTERNALS_USE_JOBSERVER
       "Let Rave/Pythia8's nested `make` inherit the outer build's GNU Make jobserver instead of using a fixed -j" ON)

if(FASER_EXTERNALS_USE_JOBSERVER AND CMAKE_GENERATOR STREQUAL "Unix Makefiles")
  # Literal, unexpanded "$(MAKE)" - see explanation above. No explicit -j:
  # the recursive make inherits concurrency from the outer jobserver.
  set(_faser_external_make "$(MAKE)")
else()
  set(_faser_external_make make -j${FASER_BUILD_PARALLEL_JOBS})
endif()

# -----------------------------------------------------------------------------
# CLHEP
# -----------------------------------------------------------------------------
# Smart default: if $GEANT4_INSTALL is set (mac_setup.sh / lxplus_setup.sh
# both export it, the latter via CVMFS's own geant4.sh) and that install
# has a bundled CLHEP *complete enough for FASER's own use* (see the
# detection logic and explanation in the FASER_BUILD_CLHEP=OFF branch
# below), default to reusing it instead of building CLHEP from source -
# so a plain `cmake -S . -B build` with no extra -D flags at all already
# skips CLHEP's own from-source build on a normal dev machine. This only
# ever *shortens* the default build: if GEANT4_INSTALL isn't set, or its
# bundled CLHEP isn't complete enough (see below - this is actually the
# common case), this is a no-op and FASER_BUILD_CLHEP still defaults to ON
# as before. An explicit -DFASER_BUILD_CLHEP=ON/OFF or -DCLHEP_ROOT=... on
# the command line always overrides this (ordinary CMake cache-variable
# precedence: option()/set(... CACHE ...) never overwrites a cache entry
# that's already defined, whether from the command line or a previous
# configure of the same build directory).
#
# "Complete enough" specifically means it has CLHEP's Matrix module, not
# just Vector: Geant4 only bundles the CLHEP modules *it* uses internally
# (Vector, Random, Geometry, Evaluator, Units, Utility - confirmed against
# an actual Geant4 11.4.2 install/source tree), which does NOT include
# Matrix - and Rave's DataFormats/CLHEP/interface/AlgebraicObjects.h needs
# CLHEP::HepMatrix/HepSymMatrix for track/vertex covariance algebra. A
# first version of this auto-detection checked only for libG4clhep itself
# and broke Rave's build with "'CLHEP/Matrix/Matrix.h' file not found";
# checking for the Matrix module too means a plain Geant4-bundled CLHEP
# (which never has it) safely falls through to building from source
# instead of silently producing an incomplete CLHEP.
set(_faser_clhep_default_build ON)
set(_faser_clhep_default_root  "")
if(DEFINED ENV{GEANT4_INSTALL})
  foreach(_faser_libdir lib lib64)
    if(EXISTS "$ENV{GEANT4_INSTALL}/${_faser_libdir}/libG4clhep${_faser_shlib_suffix}" AND
       EXISTS "$ENV{GEANT4_INSTALL}/include/Geant4/CLHEP/Matrix")
      set(_faser_clhep_default_build OFF)
      set(_faser_clhep_default_root  "$ENV{GEANT4_INSTALL}")
      break()
    endif()
  endforeach()
endif()

option(FASER_BUILD_CLHEP "Build CLHEP from source (gitlab.cern.ch/CLHEP/CLHEP)" ${_faser_clhep_default_build})
set(CLHEP_ROOT "${_faser_clhep_default_root}" CACHE PATH "Pre-installed CLHEP prefix (used when FASER_BUILD_CLHEP=OFF)")

if(FASER_BUILD_CLHEP)
  message(STATUS "CLHEP: building from source (gitlab.cern.ch/CLHEP/CLHEP) - pass -DCLHEP_ROOT=... or set $GEANT4_INSTALL to a Geant4 install with a bundled CLHEP to skip this")
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
    BUILD_COMMAND     ${CMAKE_COMMAND} --build <BINARY_DIR> --parallel ${FASER_BUILD_PARALLEL_JOBS}
    INSTALL_COMMAND   ${CMAKE_COMMAND} --build <BINARY_DIR> --target install
    BUILD_BYPRODUCTS  "${CLHEP_INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/libCLHEP${_faser_shlib_suffix}"
  )
else()
  if(NOT CLHEP_ROOT)
    message(FATAL_ERROR "FASER_BUILD_CLHEP=OFF but CLHEP_ROOT was not set to a pre-installed CLHEP prefix")
  endif()

  # CLHEP_ROOT can point at two different kinds of pre-installed CLHEP:
  #
  #  - a normal, standalone CLHEP install (lib/libCLHEP.*,
  #    include/CLHEP/...) - a Homebrew `clhep` formula, a conda-forge
  #    `clhep` package, or the CLHEP_INSTALL_DIR of another FASER build
  #    that built CLHEP from source itself; or
  #
  #  - a Geant4 install. By default (GEANT4_USE_SYSTEM_CLHEP=OFF, which is
  #    Geant4's own default and how this project's Geant4 installs are
  #    built) Geant4 compiles its own bundled copy of CLHEP as an ordinary
  #    part of its own build - source/externals/clhep in the Geant4 source
  #    tree - and installs it renamed (lib/libG4clhep.*) under a nested
  #    include path (include/Geant4/CLHEP/...) so it never collides with a
  #    real system CLHEP install. Since FASER requires a Geant4 install
  #    anyway (see the top of this file), that bundled CLHEP is already
  #    sitting right there for anyone building FASER - reusing it skips
  #    CLHEP's own from-source build (by far the heaviest single external:
  #    a fresh git clone plus a full separate CMake configure+build)
  #    entirely.
  #
  #  Rave's `--with-clhep=` autotools flag and GenFit's `-lCLHEP` link
  #  flag (both further down in this file) assume the normal, standalone
  #  layout, so a bundled Geant4 CLHEP is exposed through a small
  #  directory of symlinks (under the build tree) presenting it with that
  #  same shape - everything past this point then just uses
  #  CLHEP_INSTALL_DIR uniformly, however CLHEP was actually obtained.
  #
  #  Geant4 only bundles the CLHEP modules it uses internally (Vector,
  #  Random, Geometry, Evaluator, Units, Utility) - never Matrix, which
  #  Rave needs for CLHEP::HepMatrix/HepSymMatrix. Require Matrix to be
  #  present too, not just libG4clhep itself, so an incomplete bundled
  #  CLHEP is correctly rejected here rather than failing later with
  #  "'CLHEP/Matrix/Matrix.h' file not found" partway through Rave's build.
  if(EXISTS "${CLHEP_ROOT}/${CMAKE_INSTALL_LIBDIR}/libG4clhep${_faser_shlib_suffix}" AND
     EXISTS "${CLHEP_ROOT}/include/Geant4/CLHEP/Matrix")
    set(CLHEP_INSTALL_DIR "${FASER_EXTERNAL_STAGE_DIR}/clhep-from-geant4")
    file(MAKE_DIRECTORY "${CLHEP_INSTALL_DIR}/include" "${CLHEP_INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}")
    if(NOT EXISTS "${CLHEP_INSTALL_DIR}/include/CLHEP")
      file(CREATE_LINK "${CLHEP_ROOT}/include/Geant4/CLHEP" "${CLHEP_INSTALL_DIR}/include/CLHEP" SYMBOLIC)
    endif()
    if(NOT EXISTS "${CLHEP_INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/libCLHEP${_faser_shlib_suffix}")
      file(CREATE_LINK "${CLHEP_ROOT}/${CMAKE_INSTALL_LIBDIR}/libG4clhep${_faser_shlib_suffix}" "${CLHEP_INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/libCLHEP${_faser_shlib_suffix}" SYMBOLIC)
    endif()
    message(STATUS "CLHEP: reusing the CLHEP bundled inside the Geant4 install at ${CLHEP_ROOT} (via symlink shim ${CLHEP_INSTALL_DIR})")
  elseif(EXISTS "${CLHEP_ROOT}/${CMAKE_INSTALL_LIBDIR}/libCLHEP${_faser_shlib_suffix}")
    # A standalone CLHEP install - but don't hand CLHEP_ROOT to Rave/GenFit
    # as-is even here: CLHEP_ROOT can be a *shared* environment prefix
    # rather than a CLHEP-only one - e.g. CI passes -DCLHEP_ROOT=$CONDA_PREFIX,
    # a conda env with both `clhep` and `root` installed side by side under
    # the same include/ and lib/. Rave's own `./configure --with-clhep=...`
    # then adds that whole prefix's include/ to its build, and Rave's
    # vendored ROOT/smatrix headers fall through (via a generic #include
    # that isn't satisfied within Rave's own vendored tree) to the *real*
    # ROOT headers sitting in that same prefix - which are built for a much
    # newer C++ standard than the -std=c++11 Rave itself is compiled with,
    # producing errors like "'constexpr' constructor does not have empty
    # body" deep inside ROOT's GenVector headers when Rave's own tests/
    # get built. Isolating CLHEP behind a symlink shim, exactly like the
    # Geant4-bundled branch above already does, means Rave/GenFit only
    # ever see CLHEP's own headers, never whatever else happens to live in
    # the same prefix.
    set(CLHEP_INSTALL_DIR "${FASER_EXTERNAL_STAGE_DIR}/clhep-shim")
    file(MAKE_DIRECTORY "${CLHEP_INSTALL_DIR}/include" "${CLHEP_INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}")
    if(NOT EXISTS "${CLHEP_INSTALL_DIR}/include/CLHEP")
      file(CREATE_LINK "${CLHEP_ROOT}/include/CLHEP" "${CLHEP_INSTALL_DIR}/include/CLHEP" SYMBOLIC)
    endif()
    if(NOT EXISTS "${CLHEP_INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/libCLHEP${_faser_shlib_suffix}")
      file(CREATE_LINK "${CLHEP_ROOT}/${CMAKE_INSTALL_LIBDIR}/libCLHEP${_faser_shlib_suffix}" "${CLHEP_INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/libCLHEP${_faser_shlib_suffix}" SYMBOLIC)
    endif()
    message(STATUS "CLHEP: using the pre-installed standalone CLHEP at ${CLHEP_ROOT} (via isolating symlink shim ${CLHEP_INSTALL_DIR}, so anything else sharing that prefix - e.g. a conda env's ROOT - can't leak into Rave/GenFit's own builds)")
  else()
    # CLHEP_ROOT matched neither known layout (no libG4clhep+Matrix, and
    # no standalone libCLHEP either) - most commonly a *stale cache*: this
    # variable is a CACHE PATH, so once a build directory has been
    # configured once, CMake will keep reusing whatever value is already
    # in CMakeCache.txt on every subsequent reconfigure, even after the
    # detection logic above (or GEANT4_INSTALL, or the Geant4 install
    # itself) changes - it is never silently recomputed. Failing loudly
    # here, instead of silently accepting a bad path (as a prior version
    # of this file did), turns that into a clear, actionable CMake error
    # instead of a confusing failure deep inside Rave's own `./configure`
    # ("configure: error: required clhep not found").
    message(FATAL_ERROR
      "CLHEP_ROOT (${CLHEP_ROOT}) is neither a Geant4 install with a "
      "Matrix-complete bundled CLHEP (missing "
      "${CLHEP_ROOT}/${CMAKE_INSTALL_LIBDIR}/libG4clhep${_faser_shlib_suffix} "
      "and/or ${CLHEP_ROOT}/include/Geant4/CLHEP/Matrix) nor a standalone "
      "CLHEP install (missing "
      "${CLHEP_ROOT}/${CMAKE_INSTALL_LIBDIR}/libCLHEP${_faser_shlib_suffix}). "
      "If you recently rebuilt/moved your Geant4 install, or changed "
      "GEANT4_INSTALL, this is very likely a STALE CACHED VALUE from an "
      "earlier configure of this same build directory: CLHEP_ROOT/"
      "FASER_BUILD_CLHEP are CACHE variables and are never recomputed on "
      "their own. Fix: delete the build directory (or at least "
      "CMakeCache.txt) and reconfigure from scratch, e.g. "
      "`rm -rf build && cmake -S . -B build ...`.")
  endif()
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
  set(RAVE_INSTALL_DIR  "${FASER_EXTERNAL_INSTALL_DIR}/rave")
  # The vendored, git-tracked copy - never built in place (see below).
  set(RAVE_VENDORED_DIR "${CMAKE_SOURCE_DIR}/rave")
  # A private copy of it under the build tree, staged fresh by the
  # DOWNLOAD_COMMAND below. autoreconf/configure/make all run against
  # *this* directory, not RAVE_VENDORED_DIR.
  set(RAVE_SOURCE_DIR   "${FASER_EXTERNAL_STAGE_DIR}/rave-src")

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
  # rave && ./configure ...`) - but *in-source relative to the staged copy
  # above*, not the tracked ./rave directory. `autoreconf -fi` regenerates
  # configure/config.h.in/aclocal.m4/libtool's scaffolding using whatever
  # autoconf/automake/libtool happen to be installed on the machine doing
  # the build; running that against the tracked directory directly used to
  # leave `git status` permanently dirty with toolchain-version churn in
  # generated files (different autoconf/automake versions on different
  # machines regenerating slightly different boilerplate). Staging a fresh
  # copy under the build tree first means all of that noise lands there
  # instead, and the tracked ./rave tree is never written to by the build.
  #
  # This copy is only made once (DOWNLOAD_COMMAND runs on the first
  # configure and is not re-triggered automatically). If you edit files
  # under the vendored ./rave source itself, do a clean rebuild
  # (`rm -rf build`) to pick the changes up in a fresh staged copy.
  #
  # Only add a build-order dependency on clhep_external when it actually
  # exists - it doesn't when FASER_BUILD_CLHEP=OFF (e.g. reusing the CLHEP
  # bundled in a Geant4 install, or a conda/Homebrew CLHEP), and
  # ExternalProject_Add errors at configure time if DEPENDS names a
  # nonexistent target.
  set(_rave_deps "")
  if(TARGET clhep_external)
    list(APPEND _rave_deps clhep_external)
  endif()

  ExternalProject_Add(rave_external
    DEPENDS           ${_rave_deps}
    DOWNLOAD_COMMAND  sh -c "rm -rf ${RAVE_SOURCE_DIR} && ${CMAKE_COMMAND} -E copy_directory ${RAVE_VENDORED_DIR} ${RAVE_SOURCE_DIR}"
    SOURCE_DIR        ${RAVE_SOURCE_DIR}
    BUILD_IN_SOURCE   1
    CONFIGURE_COMMAND ${_rave_configure_cmd}
    # -Wno-deprecated-declarations: Rave's own vendored code (std::auto_ptr,
    # std::unary_function throughout) is responsible for the overwhelming
    # majority of this build's warning output otherwise - both are still
    # fully functional in this C++11 build, just deprecated in later
    # standards, so silencing them costs nothing and isn't worth patching
    # dozens of files in vendored code for.
    BUILD_COMMAND     ${_faser_external_make} CXXFLAGS=-g\ -std=c++11\ -Wno-deprecated-declarations LHEPINCPATH=.
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

  # Same reasoning as Rave's DEPENDS above: only reference clhep_external
  # when FASER actually built it.
  set(_genfit_deps rave_external)
  if(TARGET clhep_external)
    list(APPEND _genfit_deps clhep_external)
  endif()

  if(APPLE)
    list(APPEND _genfit_cmake_args
      "-DRave_LDFLAGS=-L${RAVE_INSTALL_DIR}/lib/ -lRaveBase -L${CLHEP_INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/ -lCLHEP")
    # Plain parallel build; no gtest workaround needed on Darwin (matches
    # the old Makefile's Darwin branch).
    set(_genfit_build_cmd ${CMAKE_COMMAND} --build <BINARY_DIR> --parallel ${FASER_BUILD_PARALLEL_JOBS})
  else()
    list(APPEND _genfit_cmake_args
      -DGTEST_LIBRARY=${GOOGLETEST_INSTALL_DIR}/lib/libgtest.a
      -DGTEST_INCLUDE_DIR=${GOOGLETEST_INSTALL_DIR}/include
      -DGTEST_MAIN_LIBRARY=${GOOGLETEST_INSTALL_DIR}/lib/libgtest_main.a
      "-DRave_LDFLAGS=-Wl,-rpath-link,${RAVE_INSTALL_DIR}/lib/ -L${RAVE_INSTALL_DIR}/lib/ -lRaveBase -L${CLHEP_INSTALL_DIR}/${CMAKE_INSTALL_LIBDIR}/ -lCLHEP")
    list(APPEND _genfit_deps googletest_external)

    # GenFit's own CMakeLists.txt defaults BUILD_TESTING to ON on every
    # non-Apple platform (it's already OFF on Darwin, in the APPLE branch
    # above - matched here rather than left as an upstream default FASER
    # never asked for). FASER only ever uses GenFit through the genfit2
    # shared library (GenFit::genfit2 further down), never its bundled
    # test/example binaries (fitterTests, minimalFittingExample(2),
    # measurementFactoryExample, ...) - and on Linux those extra binaries
    # additionally link ROOT::Geom (genfit2 itself only needs
    # ROOT::Core/Physics/Eve), which on a conda-forge ROOT install pulls
    # in libGeom.so/libGraf.so built against conda's own newer libstdc++
    # (providing symbols like GLIBCXX_3.4.31/CXXABI_1.3.15). The system
    # g++ doing the actual link still resolves its own implicit
    # `-lstdc++` against Ubuntu's older system libstdc++.so.6 first, which
    # doesn't have those symbols, so those extra binaries fail to link
    # with e.g. "libCore.so.6.40.04: undefined reference to
    # `__cxa_call_terminate@CXXABI_1.3.15'" - a build failure in code
    # FASER was never going to use in the first place. GenFit's own
    # ADD_GENFIT_TEST macro adds these targets EXCLUDE_FROM_ALL when
    # BUILD_TESTING is OFF, so they simply aren't built at all - sidestepping
    # the ABI mismatch entirely instead of trying to out-guess the linker.
    list(APPEND _genfit_cmake_args -DBUILD_TESTING=OFF)

    # GenFit's own CMake build has a known issue linking its gtest binaries
    # on Linux; the old Makefile worked around it with:
    #   make CXXFLAGS="-g" -j; sh CMakeFiles/gtests.dir/link.txt; make -j
    # Reproduce that exact three-step dance as a single shell command
    # (ExternalProject_Add's BUILD_COMMAND takes one command line). With
    # BUILD_TESTING=OFF above, the `gtests` target this targets doesn't
    # even get defined any more (it lives inside GenFit's own
    # IF(BUILD_TESTING) block, googletest deps and all), so the
    # `sh CMakeFiles/gtests.dir/link.txt` step is now a guaranteed no-op -
    # left in place (harmlessly, via the existing `|| true`) rather than
    # ripped out along with the now-unused GOOGLETEST_INSTALL_DIR/GTEST_*
    # args above, to keep this change scoped to the actual CI failure.
    set(_genfit_build_cmd sh -c
      "cd <BINARY_DIR> && (${CMAKE_COMMAND} --build . --parallel ${FASER_BUILD_PARALLEL_JOBS} || true) && (sh CMakeFiles/gtests.dir/link.txt || true) && ${CMAKE_COMMAND} --build . --parallel ${FASER_BUILD_PARALLEL_JOBS}")
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
    BUILD_COMMAND     ${_faser_external_make}
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
