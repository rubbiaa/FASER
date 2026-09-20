# FASER Software Installation Guide

This guide covers building FASER with the CMake build system. A legacy
Makefile-based flow is still present (see "Legacy Makefile build" below)
but CMake is now the supported path, and the steps below have been
validated with a real, end-to-end clean build (clone -> configure ->
build -> all binaries) on macOS (Apple Silicon).

---

## Prerequisites

- `git`, `cmake` (>= 3.20), `make`, a C++17 compiler (`g++`/`clang++`)
- Autotools (`automake`, `autoconf`, `libtool`, `m4`, `perl`) - needed
  because Rave is still autotools-based. On macOS specifically you need
  Homebrew's `libtool` (it provides `glibtoolize`, which `autoreconf`
  requires - macOS's own built-in `libtool` is a different, unrelated
  tool and does not satisfy this).
- Boost
- **Eigen3** - required by GenFit's own CMake build. Easy to miss because
  nothing prints a friendly error for it until GenFit gets around to
  configuring itself.
- **ROOT** (>= 6.20) and **Geant4**, already built/installed - these are
  *not* built by this project; they're expected to come from Homebrew, a
  module system, CVMFS, or a from-source build of your own, exactly as
  before.

### macOS (Homebrew)

```bash
brew install cmake boost automake autoconf libtool eigen
# ROOT and Geant4: brew install root geant4 (or use your existing install)
```

If you build ROOT and/or Geant4 from source yourself rather than via
Homebrew, see "Toolchain / stale-compiler-path issues" below - it's the
single most confusing class of error you can hit here, and it isn't
really a FASER problem.

### Ubuntu/Debian

```bash
sudo apt update
sudo apt install build-essential git cmake automake autoconf libtool m4 perl \
                  libboost-all-dev libeigen3-dev
```

Either way, source your usual environment setup first so `root-config`,
`$ROOTSYS`, and Geant4's `geant4.sh` are on your PATH/environment - e.g.
`source setup.sh` (edit it to point at your ROOT/Geant4 install first), or
`source mac_setup.sh` / `source lxplus_setup.sh`. `mac_setup.sh` also
exports `GEANT4_INSTALL` and adds it to `CMAKE_PREFIX_PATH`, so
`find_package(Geant4)` can locate it without you having to pass
`-DGeant4_DIR=...` by hand every time - edit the paths at the top of that
script to match where your own ROOT/Geant4 actually live before sourcing
it.

---

## Build

```bash
git clone https://github.com/rubbiaa/FASER.git
cd FASER
source mac_setup.sh   # or setup.sh / lxplus_setup.sh - see Prerequisites
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j
```

The first configure+build will also fetch and compile the bundled
third-party dependencies (CLHEP, Rave, GenFit, googletest on Linux,
Pythia8) - this replaces the old `make clhep`, `make rave`, `make genfit`,
`make pythia8` steps and can take a while the first time (expect Rave's
autotools build in particular to take several minutes).

All FASER executables land in `build/bin/`, for example:

```
build/bin/AnalyReco.exe            build/bin/FileMask.exe
build/bin/plots.exe                build/bin/t.exe / s.exe        (TauSearch)
build/bin/batchreco.exe            build/bin/Convert.exe          (ConvertFASERMC)
build/bin/batchreco_detresp.exe    build/bin/ConvertGENIE.exe
build/bin/dumphits.exe             build/bin/CombineFluxes.exe
build/bin/evDisplay.exe            build/bin/faserps               (FASERG4)
build/bin/FaserCalDisplayApp       build/bin/faserps_proto          (FASERCalProtoG4)
build/bin/FaserCalAnalyzer(2nd)
```

(`faserps_proto` is FASERCalProtoG4's simulation binary, renamed from its
original `faserps` only to avoid clashing with FASERG4's own `faserps` now
that both are built from one unified CMake project - they used to be
configured completely independently.)

Package-local shared libraries some workflows load interactively
(PyROOT/`gSystem->Load(...)`) are still written next to their sources, not
into `build/bin/`: `Analysis/libTPORec.so`, `Batch/libTPORec.so`,
`CoreUtils/libCoreUtilsDict.so`.

### Useful options

Pass these with `-D<OPTION>=<value>` at the `cmake -S . -B build` step:

| Option                          | Default | Meaning |
|----------------------------------|---------|---------|
| `FASER_BUILD_EXTERNALS`          | `ON`    | Build CLHEP/Rave/GenFit/googletest/Pythia8 from source |
| `FASER_BUILD_GEANT4_PACKAGES`    | `ON`    | Build FASERG4 and FASERCalProtoG4 (requires Geant4) |
| `FASER_BUILD_DISPLAY`            | `ON`    | Build the Display package (event display / muon spectrometer analysis) |
| `FASER_BUILD_TAUSEARCH`          | `ON`    | Build TauSearch |
| `FASER_BUILD_CLHEP` / `_RAVE` / `_GENFIT` / `_PYTHIA8` | `ON` | Build each dependency individually; set to `OFF` and pass `-D<NAME>_ROOT=/path/to/existing/install` to reuse an already-built copy instead |
| `FASER_BUILD_GOOGLETEST`         | `ON` on Linux only | GenFit's CMake build needs googletest on Linux (not on macOS); set `OFF` + `-DGOOGLETEST_ROOT=...` to reuse an existing one |

Already have CLHEP/Rave/GenFit/Pythia8 built (e.g. from a previous
Makefile-based build)? Point CMake at them instead of rebuilding:

```bash
cmake -S . -B build \
  -DFASER_BUILD_EXTERNALS=OFF \
  -DCLHEP_ROOT=$PWD/CLHEP-install \
  -DRAVE_ROOT=$PWD/rave-install \
  -DGENFIT_ROOT=$PWD/GenFit-install \
  -DPYTHIA8_ROOT=$PWD/pythia8312
```

### Skip building CLHEP - reuse the copy bundled inside your Geant4 install

CLHEP's own from-source build (a fresh git clone plus a full separate
`cmake`+build of its own) is the heaviest single piece of
`FASER_BUILD_EXTERNALS`. You very likely don't need it: unless your
Geant4 was built with `-DGEANT4_USE_SYSTEM_CLHEP=ON`, Geant4 already
compiled and installed its own bundled copy of CLHEP as an ordinary part
of its own build (`source/externals/clhep` in the Geant4 source tree) -
installed as `lib/libG4clhep.*` plus `include/Geant4/CLHEP/...` under
your Geant4 install prefix, purely renamed so it never collides with a
real system CLHEP.

Point `CLHEP_ROOT` at that same Geant4 install and FASER will detect the
bundled layout automatically and reuse it - no separate CLHEP build at
all:

```bash
cmake -S . -B build \
  -DFASER_BUILD_CLHEP=OFF \
  -DCLHEP_ROOT=$GEANT4_INSTALL   # mac_setup.sh already exports this
```

(Everything else - `FASER_BUILD_RAVE`, `FASER_BUILD_GENFIT`,
`FASER_BUILD_PYTHIA8` - is unaffected and still builds from source unless
you turn those off too.) Configure prints which layout it found:

```
-- CLHEP: reusing the CLHEP bundled inside the Geant4 install at ... (via symlink shim ...)
```

If your Geant4 *was* built with `GEANT4_USE_SYSTEM_CLHEP=ON` (true of the
conda-forge Geant4 package the Linux CI workflow uses, for example),
`CLHEP_ROOT` can instead point straight at wherever that system/package
CLHEP lives (e.g. `$CONDA_PREFIX`) - FASER detects that it's a normal,
standalone CLHEP install (not a Geant4-bundled one) and uses it as-is.

### Clean rebuild

```bash
rm -rf build
```

(Rave is still built in-source under `./rave/`, since its autotools build
was always run that way; `git clean -fdx rave/` clears that out if
needed. Rave's `CONFIGURE_COMMAND` runs `autoreconf -fi` itself before
`configure` on every platform now, so a clean rebuild doesn't require you
to do anything extra there - see "What the CMake modernization changed".)

---

## Troubleshooting a clean build

These are real errors hit while validating this exact process end-to-end;
if your first build fails, check here before opening an issue.

**`Could not find a package configuration file ... ROOT` / `find_package(ROOT)` fails**
ROOT's environment isn't sourced. Run `root-config --prefix` to find your
install, then either source `mac_setup.sh`/`setup.sh` (recommended) or
pass `-DROOT_DIR=<that prefix>/cmake` directly.

**`find_package(Geant4)` fails**
Same idea: locate it with `geant4-config --prefix` or `brew --prefix
geant4`, then pass `-DGeant4_DIR=<prefix>/lib/Geant4-<version>`, or set
`-DFASER_BUILD_GEANT4_PACKAGES=OFF` if you don't need the Geant4-based
packages at all.

**`Cannot find header ....hh to generate dictionary` / `fatal error: '....hh' file not found` during dictionary generation**
This means a `root_generate_dictionary()` call is missing an `OPTIONS
-I<path>` entry for wherever that header actually lives (CoreUtils,
GenFit, etc). Every package's `CMakeLists.txt` in this repo already has
the right `OPTIONS -I...` flags for its own dependencies as of this
writing; if you add a new header include to a class that's part of a
dictionary, you may need to add its directory the same way.

**`Could NOT find Eigen3 (missing: EIGEN3_INCLUDE_DIR EIGEN3_VERSION_OK)`**
Eigen3 isn't installed (`brew install eigen` / `apt install
libeigen3-dev`). If it's installed and this still fails with a *blank*
found version (e.g. "Eigen3 version .. found ... but at least version
2.91.0 is required"), that's GenFit's own bundled, pre-Eigen-3.3
`cmake/FindEigen3.cmake` failing to parse a modern Eigen header - our
`Externals.cmake` already deletes that file from GenFit's checkout as
part of its patch step specifically to avoid this, so a fresh clone
shouldn't hit it; if you do, your `build/external/GenFit` checkout is
probably stale from before that fix - `rm -rf build/external/GenFit
build/external-install/GenFit` and rebuild.

**`Can't exec "glibtoolize"` during Rave's configure step**
Install Homebrew's `libtool` (`brew install libtool`). This is GNU
libtool, not the same thing as the `libtool` macOS ships by default.

**`Reversed (or previously applied) patch detected!` during GenFit's patch step**
This used to happen on every rebuild after the first, because GenFit
tracks the `main` branch (re-fetched on every build) without reverting
the previous run's already-applied local patch. The `PATCH_COMMAND` in
`cmake/Externals.cmake` now uses `patch -N` plus a tolerant `|| true` so
this is a non-issue on a fresh clone; you should never see it.

**`config.status: error: cannot find input file: 'Makefile.in'` during Rave's configure step**
The vendored `./rave` source ships `configure` but not the `Makefile.in`
files it expects (autotools was never run on this particular checkout).
`cmake/Externals.cmake` runs `autoreconf -fi` before `configure` on every
platform to regenerate them, so this shouldn't happen on a fresh clone
either - if it does, you're likely missing `automake`/`autoconf`/`libtool`
(see Prerequisites).

### Toolchain / stale-compiler-path issues

If you build ROOT or Geant4 from source yourself (rather than via
Homebrew) and later switch your active Xcode/Command Line Tools version,
you can hit errors like:

```
No rule to make target '/Applications/Xcode-beta.app/.../libexpat.tbd'
```

or dictionary-generation failures mentioning a sysroot under a specific
Xcode version that no longer matches your current `xcode-select -p`.
This happens because `rootcling`/`cling` and some CMake `find_library()`
results bake in *absolute* paths recorded at the time ROOT/Geant4 were
themselves built or configured - they don't automatically re-resolve
against whatever Xcode is active now. If you hit this:

- For a dependency's own exported CMake cache (e.g. Geant4's
  `lib/cmake/Geant4/Geant4PackageCache.cmake`), you can often just patch
  the stale absolute path directly to the correct one for your current
  SDK (`xcrun --show-sdk-path`) - much faster than a full rebuild, but
  you'll need to fix it again if the cached copy in an existing `build/`
  directory also has the stale value baked into its own `CMakeCache.txt`.
- For `rootcling` itself, the compiler path is baked in at the time ROOT
  was built, and there's no equivalent cache file to patch - you need to
  rebuild ROOT under your currently-active Xcode.

This is a pre-existing toolchain hazard on macOS, unrelated to the CMake
modernization itself - it would affect the old Makefile build exactly the
same way, and it only ever comes up if you build ROOT/Geant4 from source
and then change Xcode versions afterward.

---

## What the CMake modernization changed

- One `cmake -S . -B build && cmake --build build -j` replaces the old
  `make clhep && make rave && make genfit && make pythia8`, followed by a
  separate `make` inside every package directory.
- `cmake/Externals.cmake` builds CLHEP, Rave, GenFit, googletest, and
  Pythia8 as a proper dependency graph (`ExternalProject_Add`), instead of
  hand-written `wget`/`git clone`/`configure` recipes. The git-based ones
  are pinned to explicit branches (CLHEP: `develop`, GenFit/googletest:
  `main`) rather than relying on an implicit default, since GenFit and
  googletest no longer default to `master`.
- Rave's `CONFIGURE_COMMAND` runs `autoreconf -fi` on every platform
  (not just Linux) before `./configure`, since the vendored `./rave`
  source ships `configure` but not the `Makefile.in` files it needs.
- GenFit's `PATCH_COMMAND` is idempotent (`patch -N ... || true`) so it
  survives being re-run on every build without erroring on an
  already-applied patch, and it also strips GenFit's own bundled
  pre-Eigen-3.3 `cmake/FindEigen3.cmake`, which can't parse the version
  out of modern Eigen headers - this lets `find_package(Eigen3)` fall
  through to Eigen's own proper exported CMake config instead.
- Every package that had a raw Makefile (`CoreUtils`, `Analysis`, `Batch`,
  `ConvertFASERMC`, `ConvertGENIE`, `EvDisplay`, `FileMask`, `TauSearch`)
  now has a `CMakeLists.txt` that reproduces the same dictionaries,
  libraries and executables. `FASERG4`, `FASERCalProtoG4` and `Display`
  already had their own `CMakeLists.txt`; they're now wired into the same
  top-level build instead of being configured independently.
- `Display/CMakeLists.txt` had a GenFit path hardcoded to one developer's
  `/eos/user/u/ukose/...` install; it now picks up GenFit from the
  superbuild (or `-DGENFIT_ROOT=...`).
- `FASERCalProtoG4`'s executable was renamed `faserps` -> `faserps_proto`
  (see above), to resolve both an executable-name and a dictionary-target
  collision with FASERG4's own `faserps` once both packages moved into
  one unified CMake project.
- Every `root_generate_dictionary()` call across the affected packages
  now passes explicit `OPTIONS -I<path>` flags (for CoreUtils headers and
  for GenFit's headers, as needed) so `rootcling` can resolve `#include`s
  written inside each package's own headers and `LinkDef.h` - previously
  these relied on compiler include paths that only applied to normal
  compilation, not to dictionary generation, which is a separate
  invocation of `rootcling`.
- `Analysis/LinkDef.h` and `FileMask/LinkDef.h` were missing
  `#pragma link C++ class TMuTrack;` (present in `CoreUtils/LinkDef.h` and
  `Batch/LinkDef.h`) even though both packages compile code that uses
  `TMuTrack` - this meant `TMuTrack`'s vtable was never emitted for those
  two packages' link, causing "vtable for TMuTrack" undefined-symbol
  errors. Fixed by adding the missing dictionary entries.
- `EvDisplay/MyMainFrame.h` had no include guard, which caused a
  redefinition error once it was (correctly) included both directly and
  transitively via `LinkDef.h`. Fixed with `#pragma once`.
- `EvDisplay/CMakeLists.txt` compiled with `-fsanitize=address` but never
  linked it in (an old Makefile quirk, kept only for parity) - this left
  every AddressSanitizer runtime symbol unresolved at link time. Removed.

### Validated

This has been build-tested end-to-end on macOS (Apple Silicon): a clean
`git clone` followed by the exact commands in "Build" above successfully
configures and builds every package, including the full CLHEP/Rave/
GenFit/Pythia8 superbuild, GenFit-based tracking, and the event display.
The Linux/lxplus path uses the same `CMakeLists.txt`/`Externals.cmake`
logic (with its own already-established branches, e.g. googletest and the
gtest-link workaround for GenFit) but has not been independently
re-verified in this same end-to-end pass - if you hit something there,
please report it.

---

## Legacy Makefile build

The original hand-rolled `make clhep && make rave && make genfit &&
make pythia8`, followed by `make` inside each package directory, still
works and is unchanged. It's kept around during the transition to CMake
but is no longer the recommended path - see the top-level `Makefile` for
its targets.

---

## Support

If you encounter problems, open an issue on GitHub:
https://github.com/rubbiaa/FASER/issues
