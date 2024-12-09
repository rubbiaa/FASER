void plotcharm() {
    gSystem->Load("libTPORec.so");
    TFile f("Batch-TPORecevent_1.root");

    TH1D* h_charm_Enuall;
    TH1D* h_charm_Enucharmed;
    f.GetObject("h_charm_Enuall",h_charm_Enuall);
    f.GetObject("h_charm_Enucharmed",h_charm_Enucharmed);


    TH1F *hRatio = (TH1F*)h_charm_Enucharmed->Clone("hRatio");
    hRatio->SetTitle("Charmed events fraction");
    hRatio->Divide(h_charm_Enuall);

TCanvas* c1 = new TCanvas("c1", "Charm Energy", 800, 600);
 
hRatio->SetMarkerStyle(21);
hRatio->SetMarkerColor(kBlue);
hRatio->SetLineColor(kBlue);

TLine *line = new TLine(hRatio->GetXaxis()->GetXmin(), 1, hRatio->GetXaxis()->GetXmax(), 1);
line->SetLineColor(kRed);
line->SetLineStyle(2);
line->Draw("same");
    hRatio->Draw();
 c1->Modified();
 c1->Update();
c1->SaveAs("charm.png");
}
