// Moved here from GeomGDML/ (which held obsolete, hand-copied .gdml
// snapshots) - this now imports the single current geometry that
// FASERG4/src/DetectorConstruction.cc itself writes out, under
// $FASERDATA/GDML/ (see CoreUtils/FaserDataDir.hh). Run this after a
// faserps run has generated that file, e.g. `root -l FASERG4/PlotGeometry.C`
// from a shell with setup.sh/common_setup.sh sourced (so $FASERDATA is set).
{
    const char* faserdata = gSystem->Getenv("FASERDATA");
    if (!faserdata || !*faserdata) {
        Error("PlotGeometry", "FASERDATA is not set - source setup.sh (or "
              "common_setup.sh) before running this macro.");
    } else {
        // Load geometry
        TGeoManager::Import(TString::Format("%s/GDML/FASERCAL_V10.gdml", faserdata));
        gGeoManager->GetTopVolume()->Draw("ogl");
    }
}
