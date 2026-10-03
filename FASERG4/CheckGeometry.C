// Moved here from GeomGDML/ (which held obsolete, hand-copied .gdml
// snapshots) - this now imports the single current geometry that
// FASERG4/src/DetectorConstruction.cc itself writes out, under
// $FASERDATA/GDML/ (see CoreUtils/FaserDataDir.hh). Run this after a
// faserps run has generated that file, e.g. `root -l FASERG4/CheckGeometry.C`
// from a shell with setup.sh/common_setup.sh sourced (so $FASERDATA is set).
{
  const char* faserdata = gSystem->Getenv("FASERDATA");
  if (!faserdata || !*faserdata) {
    Error("CheckGeometry", "FASERDATA is not set - source setup.sh (or "
          "common_setup.sh) before running this macro.");
    return;
  }
  TString gdmlPath = TString::Format("%s/GDML/FASERCAL_V10.gdml", faserdata);
  TGeoManager::Import(gdmlPath);
  gGeoManager->GetTopVolume()->Print();        // quick overview


  TGeoNode* world = gGeoManager->GetTopNode();
  std::cout << "Top node: " << world->GetName()
          << "  vol=" << world->GetVolume()->GetName() << std::endl;

  TGeoNode* detAsm = world->GetDaughter(0);  // this is the one you printed
  std::cout << "DetectorAssembly node: " << detAsm->GetName()
          << " vol=" << detAsm->GetVolume()->GetName() << std::endl;

  int nd = detAsm->GetNdaughters();
  for (int i = 0; i < nd; ++i) {
    TGeoNode* d = detAsm->GetDaughter(i);
    std::cout << "  child " << i
              << " node=" << d->GetName()
              << "  vol="  << d->GetVolume()->GetName()
              << std::endl;
  }

















  gGeoManager->GetListOfVolumes()->ls();       // list volumes










  TObjArray* vols = gGeoManager->GetListOfVolumes();
  for (int i = 0; i < vols->GetEntriesFast(); ++i) {
    TGeoVolume* v = (TGeoVolume*)vols->At(i);
    if (!v) continue;
    TString n = v->GetName();
    if (n.Contains("HCal", TString::kIgnoreCase) ||
	n.Contains("HCAL", TString::kIgnoreCase) ||
	n.Contains("rear", TString::kIgnoreCase)) {
      std::cout << i << "  " << n << std::endl;
    }
  }

  
}
