# =============================================================================
# TauSearch: t.exe, s.exe
# Standalone generator-level tau search code; only needs ROOT (no CoreUtils,
# no GenFit).
# =============================================================================
add_executable(t.exe t.C faserntuplib.cc)
target_link_libraries(t.exe PRIVATE ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::EG)

add_executable(s.exe
  s.C
  TKinFitter.cc
  TFitParticleESpher.cc
  TAbsFitParticle.cc
  TFitConstraintEp.cc
  TAbsFitConstraint.cc
)
target_link_libraries(s.exe PRIVATE ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::EG)
