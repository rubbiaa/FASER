# MuonDIS (FASERG4 port)

Ports the muon deep-inelastic-scattering enhancement from the Calypso
`Simulation/G4Extensions/MuonDIS` package (Athena/ATLAS-FASER offline
framework) into this standalone Geant4 app.

## What changed vs. the Calypso version

Calypso's version wraps a Gaudi `AthAlgTool`/`ComponentAccumulator` around the
physics, and reads a pre-generated POWHEG event library via HepMC3
(`MuonDISPowhegEventIndex`), picking the library event whose incoming muon
energy is closest to the actual Geant4 muon and rotating its final state into
the true track direction.

This app has no Gaudi/Athena and no HepMC3 dependency, so instead:

- `MuonDISPhysics` is a plain `G4VPhysicsConstructor`, registered directly in
  `faserps.cc` (`physicsList->RegisterPhysics(new MuonDISPhysics())`),
  following the same pattern as `TauDecayPhysics` in this directory.
- Configuration is via Geant4 UI macro commands (`MuonDISMessenger`, under
  `/physics/muondis/...`) instead of Gaudi properties -- see below.
- `MuonDISPythiaGenerator` generates the DIS final state **on the fly**, one
  Pythia8 call per interaction, using the actual Geant4 muon energy and
  direction as the incoming beam (`Beams:frameType = 3`). This is more exact
  than the library-lookup approach (no closest-energy matching, no
  generator-frame rotation) but means Pythia8 is now a real link-time
  dependency of `faserps` (see Build below).
- `MuonDISNuclearWrapperProcess` and `MuonDISInteractionRecorder` are ported
  essentially unchanged -- they only ever depended on Geant4/STL, not on
  Athena or HepMC3.

## Physics settings -- please review before trusting the output

This was written without a local Geant4/ROOT/Pythia8 build to test against
(this environment doesn't have them installed), so treat the Pythia8
configuration in `MuonDISPythiaGenerator.cc` as a first cut, cross-checked
against the online Pythia8 manual but not run:

- Process: `WeakBosonExchange:ff2ff(t:gmZ)` (neutral current, photon/Z
  exchange). Charged current (`t:W`) is off -- negligible at these energies
  for a charged lepton beam.
- `PhaseSpace:Q2Min` defaults to 1.0 GeV^2 (`/physics/muondis/q2min`) to stay
  away from the divergent photoproduction region. Pick whatever your analysis
  needs.
- Nuclear target: for each interaction, the struck nucleon (Pythia8 beam B)
  is randomly chosen as a free proton or neutron with probability Z/A and
  (A-Z)/A (isoscalar mix), at rest -- no Fermi motion, binding energy,
  shadowing or EMC effect. Pythia8 ships no nuclear PDFs itself; see
  `MuonDISPythiaGenerator::chooseNucleon()` if you want to refine this later.
- The recoil nucleus added to the Geant4 secondaries is at rest (no recoil
  kinematics against the struck nucleon) -- see the comment on
  `GetFatalEnergyCheckLevels()` in `MuonDISMuonNuclearModel.cc`.
- Every interaction currently does a full Pythia8 `init()` (beam
  species/kinematics differ essentially every call). If this becomes a
  performance issue, look at `Beams:allowVariableEnergy` and the
  `setKinematics()` overload for `frameType=3` in your installed Pythia8's
  `Pythia.h`.

## Build

MuonDIS needs Pythia8 built at `fasermuondis/pythia8312` (the top-level
`Makefile`'s `pythia8` target does this: `make pythia8` from the
`fasermuondis` root). `FASERG4/CMakeLists.txt` points at
`../pythia8312/{include,lib/libpythia8.a}` relative to the CMake source dir;
override with `-DPYTHIA8_DIR=/path/to/pythia8312` if yours lives elsewhere.
CMake will warn (not fail) at configure time if the library isn't found yet.

## Enabling it in a run macro

All `/physics/muondis/...` commands configure the physics constructor and
must be issued **before** `/run/initialize` (that's when `ConstructProcess()`
actually installs the process on mu-/mu+); commands issued after have no
effect. See `runFASER_muondis.mac` for a full example combined with the
existing `/generator/wantMuonBackground` single-muon mode.

```
/physics/muondis/enable true
/physics/muondis/crossSectionBias 150
/physics/muondis/q2min 1.0
/physics/muondis/interactionLog muondis_interactions.csv
#/physics/muondis/debug true
/run/initialize
```

MuonDIS is disabled by default -- nothing changes for existing macros that
don't set `/physics/muondis/enable true`.
