{
  TGeoManager::Import("geometry_tilted_5degree.gdml");
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
