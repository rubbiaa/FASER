// GeomGDML/ (obsolete, hand-copied .gdml snapshots) is gone -- this now
// imports the single current geometry FASERG4/src/DetectorConstruction.cc
// itself writes out, under $FASERDATA/GDML/ (see CoreUtils/FaserDataDir.hh).
// Run with setup.sh/common_setup.sh sourced first (so $FASERDATA is set).
gSystem->Load("libGeom");
{
  const char* faserdata = gSystem->Getenv("FASERDATA");
  if (!faserdata || !*faserdata) {
    Error("browseGeom", "FASERDATA is not set - source setup.sh (or "
          "common_setup.sh) before running this macro.");
  } else {
    TGeoManager::Import(TString::Format("%s/GDML/FASERCAL_V10.gdml", faserdata));
    gGeoManager->GetTopVolume()->Draw("ogl");
  }
}
