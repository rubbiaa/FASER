# Review: Turning off Pythia8 decay for Λ, Σ±, Ξ0, Ξ±, Ω⁻ in `CoreUtils/TPOEvent.cc`

Review of a proposed 7-line diff to `initialize_pythia()` that adds:

```cpp
fPythia8->ReadString("3222:onMode = off");  // Sigma+
fPythia8->ReadString("3112:onMode = off");  // Sigma-
fPythia8->ReadString("3122:onMode = off");  // Lambda0
fPythia8->ReadString("3322:onMode = off");  // Xi0
fPythia8->ReadString("3312:onMode = off");  // Xi-
fPythia8->ReadString("3334:onMode = off");  // Omega-
```

alongside the existing `111`/`310`/`130` (π0/K0S/K0L) lines.

## 1. c·τ table

Values cross-checked between the current PDG Review of Particle Physics and the actual Geant4 v11.4.2 source FASER links against (`particles/hadrons/barions/src/G4*.cc`).

| Particle | PDG code | Mass (MeV) | τ (PDG) | cτ (PDG) | Geant4's hardcoded τ | Geant4's cτ |
|---|---|---|---|---|---|---|
| Λ⁰ | 3122 | 1115.683 | 2.617×10⁻¹⁰ s | 7.845 cm | 2.631×10⁻¹⁰ s | 7.89 cm |
| Σ⁺ | 3222 | 1189.37 | 0.8018×10⁻¹⁰ s | 2.404 cm | 0.8018×10⁻¹⁰ s | 2.40 cm |
| Σ⁻ | 3112 | 1197.449 | 1.479×10⁻¹⁰ s | 4.434 cm | 1.479×10⁻¹⁰ s | 4.43 cm |
| Ξ⁰ | 3322 | 1314.82 | 2.90×10⁻¹⁰ s | 8.71 cm | 2.90×10⁻¹⁰ s | 8.69 cm |
| Ξ⁻ | 3312 | 1321.70 | 1.639×10⁻¹⁰ s | 4.91 cm | 1.639×10⁻¹⁰ s | 4.91 cm |
| Ω⁻ | 3334 | 1672.43 | 0.821×10⁻¹⁰ s | 2.461 cm | 0.821×10⁻¹⁰ s | 2.46 cm |

Geant4 matches PDG essentially exactly for five of six; only Λ⁰'s lifetime is ~0.5% high (Geant4 hasn't picked up the latest PDG average). Not a defect in the reviewed PR — just a footnote. All six have macroscopic (cm-scale) decay lengths, i.e. exactly the class of particle that should be tracked and decayed *in the detector*, not force-decayed at the primary vertex.

## 2. Does Geant4 process these particles correctly?

Verified from source, not from general recollection:

- `FASERG4/faserps.cc` builds its physics with `FTFP_BERT` (`new FTFP_BERT`), Geant4's standard reference list.
- `FTFP_BERT.cc` registers `G4DecayPhysics`.
- `G4DecayPhysics::ConstructParticle()` → `G4BaryonConstructor`, which explicitly instantiates `G4Lambda`, `G4SigmaPlus`, `G4SigmaMinus`, `G4XiZero`, `G4XiMinus`, `G4OmegaMinus` (+ antiparticles) with the mass/lifetime above and `G4PhaseSpaceDecayChannel` decay tables whose branching ratios match PDG's dominant modes (e.g. Ω⁻ → ΛK⁻ 67.8%, Ξ⁰π⁻ 23.6%, Ξ⁻π⁰ 8.6%).
- `G4DecayPhysics::ConstructProcess()` loops over *every* particle with no exclusion list and attaches `G4Decay` to any particle the process reports applicable — all six get in-flight decay automatically.
- `FASERG4/src/PrimaryGeneratorAction.cc` turns `TPOEvent::PO` records into primaries via `particleTable->FindParticle(aPO.m_pdg_id)` — a generic PDG lookup, not a hardcoded switch. The only PDG-restrictive filter in that file is scoped to a debug-only single-pion mode, not the normal event loop.
- `ParticleManager.cc` records tracks/hits by PDG code generically.
- `TrackerSD.cc`'s only track-kill logic is scoped to `opticalphoton`.

**Conclusion: no FASER-side code changes are needed.** The existing PDG-driven plumbing already covers these six particles as both primaries and secondaries.

## 3. What the diff actually does

`initialize_pythia()`'s `onMode = off` list only runs inside `TPOEvent::perform_taulepton_decay()` / `TPOEvent::perform_charmhadron_decay()`, both compiled only under `_INCLUDE_PYTHIA_` — defined solely for `ConvertGENIE.exe` and `Convert.exe` (not `FASERG4`/`Batch`). Those two functions inject a single tau lepton or charm hadron into a bare Pythia8 instance and let it decay per Pythia's tables. The existing π0/K0S/K0L lines already stop Pythia at those particles inside a tau/charm-hadron decay chain so Geant4 does the in-flight decay at the correct boosted vertex. The new lines extend exactly that pattern to Λ/Σ±/Ξ0/Ξ±/Ω⁻ — physically well-motivated, since real charm-baryon decays (Λc⁺→Λπ⁺, Ξc→Ξπ, Ωc⁰→Ω⁻π⁺, ...) do produce these hyperons, and without the fix Pythia would force-decay them at the primary vertex instead of letting Geant4 displace them correctly.

The "+-" in the accompanying comment is charge-conjugate shorthand — Pythia's `onMode` is shared between a particle and its antiparticle, so no separate `-3122`/`-3222`/etc. lines are needed. The six lines are the complete, correct set.

**Verdict: correct, complete, consistent with the existing π0/K0S/K0L precedent. Recommend merging as is.**

## 4. Do Pythia8 and Geant4 agree on the decay kinematics themselves?

This is the deeper question: even with the *lifetime* correct on both sides, do the two codes' decay *matrix elements* (the physics that shapes daughter angles/energies) agree? The answer is different for hyperons than for tau/charm.

### Hyperons (Λ, Σ±, Ξ0, Ξ±, Ω⁻): yes, they agree

- Geant4's `G4PhaseSpaceDecayChannel` (used for all six) is flat, isotropic phase space — no polarization handling anywhere in its source.
- Pythia8's manual (`ParticleDecays.html`, `meMode` table) shows no hyperon-specific or decay-asymmetry (α parameter, e.g. α_Λ ≈ 0.73 for Λ→pπ⁻) physics; hyperons fall through to the same generic `meMode 0` flat phase space as most tabulated decays.
- All six primary decay modes are strictly two-body (Λ→pπ⁻/nπ⁰, Σ⁺→pπ⁰/nπ⁺, Σ⁻→nπ⁻, Ξ⁰→Λπ⁰, Ξ⁻→Λπ⁻, Ω⁻→ΛK⁻/Ξ⁰π⁻/Ξ⁻π⁰). Two-body kinematics are fully fixed by mass — no Dalitz-plot shape freedom exists to get "wrong." Only the *angular* distribution could differ, and only if the parent is polarized.
- `struct PO` (`CoreUtils/TPOEvent.hh`) carries no spin/polarization field at all, so neither decayer receives polarization info to act on regardless of what it supports.

**Same model (flat 2-body phase space), same masses, no polarization info flowing through either path → statistically equivalent results.** No matrix-element discrepancy for this PR.

### Tau and charm: no, they genuinely disagree — confirmed from source

| | Geant4 | Pythia8 |
|---|---|---|
| **Charm hadrons** (D±, D⁰, Λc⁺, ...) | **No decay table at all.** Checked `G4LambdacPlus.cc` and `G4DMesonPlus.cc` directly — both construct the particle with `nullptr` for the decay table and never register one. Geant4 cannot decay these on its own. | Generic flat phase space for hadronic channels (`meMode 0`); a generic V−A spectator treatment (not full form factors) for semileptonic channels (`meMode 22/23`). Not fully realistic, but present. |
| **Tau** | Has a decay table (`G4TauMinus.cc`) with roughly correct branching ratios, but its own source comment on the leptonic channel (`G4TauLeptonicDecayChannel`) says: *"this version neglects muon polarization... assumes the pure V-A coupling... gives incorrect energy spectrum for neutrinos."* The angular distribution is drawn isotropically (`costheta = 2*G4UniformRand()-1`). Hadronic channels (π⁻ν, ρ→ππ⁰ν, a1→3πν) use plain `G4PhaseSpaceDecayChannel` — no ρ(770)/a1(1260) resonance lineshape at all. | Since v8.150, a dedicated `TauDecays` package (Ilten, arXiv:1211.6730; on by default, `TauDecays:mode = 1`) implements proper V−A helicity amplitudes with full spin correlations, plus Kühn–Santamaria hadronic currents giving realistic ρ/a1 resonance shapes. |

For charm it isn't "different matrix element" so much as "Geant4 has zero decay physics for it" — exactly why `perform_charmhadron_decay` exists at all. For tau, Geant4's model is a real but crude approximation (roughly correct branching ratios, wrong angular/spin structure, wrong neutrino spectrum, no resonance shapes), while Pythia8's is a proper matrix-element calculation — a genuine, well-documented discrepancy, distinct from the hyperon case above.

## 5. Separate, pre-existing gap worth flagging

`perform_taulepton_decay` passes `0.0` for the tau's polarization into Pythia8's `event.append()`, with an explicit `// TODO: specify polarization!!` comment in the code. So even Pythia8's superior spin-correlation machinery isn't fully exploited yet for taus that are physically polarized (as they are in real ν_τ CC interactions). Unrelated to the hyperon PR, but relevant if FASER's tau-polarimetry program cares about the production/decay kinematic correlation — would need its own follow-up.

## Sources

- PDG live particle listings (pdgLive) for Λ0/Σ+/Σ-/Ξ0/Ξ-/Ω-, and the PDG Baryon Summary Table.
- Geant4 v11.4.2 source: `particles/hadrons/barions/src/G4{Lambda,SigmaPlus,SigmaMinus,XiZero,XiMinus,OmegaMinus,LambdacPlus}.cc`, `particles/management/src/G4TauLeptonicDecayChannel.cc`, `particles/leptons/src/G4TauMinus.cc`, `particles/hadrons/mesons/src/G4DMesonPlus.cc`, `physics_lists/constructors/decay/src/G4DecayPhysics.cc`, `physics_lists/lists/src/FTFP_BERT.cc`.
- FASER source: `CoreUtils/TPOEvent.cc` / `.hh`, `FASERG4/src/PrimaryGeneratorAction.cc`, `FASERG4/src/ParticleManager.cc`, `FASERG4/src/TrackerSD.cc`, `FASERG4/faserps.cc`.
- Ilten, *"Tau Decays in Pythia 8,"* arXiv:1211.6730.
- Sjöstrand et al., *"An Introduction to PYTHIA 8.2,"* arXiv:1410.3012.
- Pythia 8 online manual, `ParticleDecays.html` (`meMode` table).
