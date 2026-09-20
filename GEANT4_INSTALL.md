# Installing Geant4 for FASER

FASER doesn't build Geant4 itself (see `cmake/Externals.cmake` - it's
treated the same way as ROOT: expected to already be installed, and
located via `find_package(Geant4 ...)`). This is a short guide for
building a Geant4 install that actually has everything FASER needs.

FASER requires the `ui_all`, `vis_all` and `gdml` Geant4 components
(`FASERG4/CMakeLists.txt` / `FASERCalProtoG4/CMakeLists.txt`):

```cmake
find_package(Geant4 REQUIRED ui_all vis_all gdml)
```

A "plain" Geant4 build without extra flags typically won't have GDML
support, which will fail at FASER's `cmake` configure step (or, if a
stale `CMakeCache.txt` masks it, later at compile time with something
like `fatal error: 'G4GDMLParser.hh' file not found`).

## 1. Prerequisites

### macOS (Homebrew)

```bash
brew install cmake xerces-c qt
```

- `xerces-c` - required for GDML support
- `qt` - gives Geant4 a real interactive UI shell and OpenGL visualization
  driver (the standard modern choice on macOS)

### Linux (Ubuntu / Debian)

```bash
sudo apt-get update
sudo apt-get install -y \
  build-essential cmake \
  libxerces-c-dev \
  qtbase5-dev qttools5-dev \
  libx11-dev libxpm-dev libxft-dev libxext-dev libxmu-dev libxi-dev \
  freeglut3-dev libglu1-mesa-dev
```

- `libxerces-c-dev` - required for GDML support (same role as `xerces-c`
  above)
- `qtbase5-dev` / `qttools5-dev` - Qt5 UI + OpenGL visualization driver.
  On Ubuntu 24.04+, where apt defaults to Qt6, use `qt6-base-dev` instead
  and add `-DGEANT4_USE_QT6=ON` to the configure step below.
- the `libx*`/`freeglut3-dev`/`libglu1-mesa-dev` set - headers needed to
  build Geant4's OpenGL/X11 visualization driver
  (`-DGEANT4_USE_OPENGL_X11=ON` below). Already validated working in this
  project's own CI (`.github/workflows/build.yml`, `ubuntu-22.04`).

Or skip building Geant4 on Linux entirely - see the conda-forge note at
the bottom, which is what CI actually uses.

## 2. Get the source

Download/clone a Geant4 release from <https://github.com/Geant4/geant4>
(pick a tagged release, e.g. `v11.4.2`), or use a copy you already have.

## 3. Configure

```bash
cmake -S /path/to/geant4-source -B /path/to/geant4-build \
  -DCMAKE_INSTALL_PREFIX=/path/to/geant4-install \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DGEANT4_INSTALL_DATA=ON \
  -DGEANT4_USE_GDML=ON \
  -DGEANT4_USE_QT=ON \
  -DGEANT4_USE_OPENGL_X11=ON \
  -DGEANT4_USE_SYSTEM_EXPAT=OFF
```

What these do:

| Flag | Why |
|---|---|
| `GEANT4_INSTALL_DATA=ON` | Downloads and installs the physics data libraries (G4NDL, G4EMLOW, PhotonEvaporation, ...) as part of the build, so you don't have to fetch/set them up separately. Needs network access during the build. |
| `GEANT4_USE_GDML=ON` | Satisfies FASER's `gdml` component request. Needs `xerces-c`. |
| `GEANT4_USE_QT=ON` / `GEANT4_USE_OPENGL_X11=ON` | Give the `ui_all`/`vis_all` components an actual interactive UI + visualization driver to find. |
| `GEANT4_USE_SYSTEM_EXPAT=OFF` | Builds Geant4's own bundled expat instead of linking the macOS SDK's copy. Avoids a real bug already hit in this project: a hardcoded `libexpat.tbd` path baked in at build time that breaks later when Xcode is updated/moved. |

Deliberately **not** set: `GEANT4_USE_SYSTEM_CLHEP` (left at its default,
`OFF`). Geant4's own bundled CLHEP only includes the modules Geant4 uses
internally (Vector, Random, Geometry, Evaluator, Units, Utility) - never
Matrix - so it's not a substitute for FASER's own CLHEP either way; no
need to change this default.

## 4. Build and install

```bash
cmake --build /path/to/geant4-build -j
cmake --build /path/to/geant4-build --target install
```

This takes a while - Geant4 is a large codebase.

## 5. Point FASER at it

Update `GEANT4_INSTALL` in whichever setup script you use
(`mac_setup.sh` / `lxplus_setup.sh`):

```bash
export GEANT4_INSTALL=/path/to/geant4-install
source $GEANT4_INSTALL/bin/geant4.sh
export CMAKE_PREFIX_PATH=$GEANT4_INSTALL:$CMAKE_PREFIX_PATH
```

Then `source` that script before configuring FASER.

## 6. Verify

Check what got baked into a Geant4 build's own cache, or what FASER
actually picked up:

```bash
# What was this Geant4 build actually configured with?
grep -E "GEANT4_USE_GDML|GEANT4_USE_SYSTEM_CLHEP|GEANT4_USE_QT|CMAKE_INSTALL_PREFIX" \
  /path/to/geant4-build/CMakeCache.txt

# Which Geant4 install did FASER's own build configure against?
grep -E "^Geant4_DIR" /path/to/FASER/build/CMakeCache.txt
```

FASER's own `cmake -S . -B build` configure step should print a line
confirming Geant4 was found, e.g.:

```
-- Found Geant4: /path/to/geant4-install/lib/cmake/Geant4/Geant4Config.cmake (found version "11.4.2")
```

## Alternative: skip building Geant4 entirely

- **CVMFS** (CERN/lxplus): source a pre-built Geant4 from `/cvmfs/geant4.cern.ch/...`
  (see `lxplus_setup.sh`) - already has GDML/vis/data.
- **conda-forge**: `mamba install -c conda-forge geant4` gets you a
  working install (with GDML) with no manual build at all - this is what
  the project's GitHub Actions CI uses (`.github/workflows/build.yml`).
