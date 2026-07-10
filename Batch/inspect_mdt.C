// inspect_mdt.C
// Inspect MDT reconstruction output from Batch-TPORecevent_*.root
//
// Usage (from FASER/Batch/):
//   root -l inspect_mdt.C
// or with a specific file:
//   root -l 'inspect_mdt.C("Batch-TPORecevent_6000_0_500.root")'

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLine.h"
#include "TLatex.h"
#include <iostream>
#include <cmath>

void inspect_mdt(const char* fname = "Batch-TPORecevent_6000_0_500.root")
{
    gStyle->SetOptStat("emr");
    gStyle->SetOptTitle(1);

    // -----------------------------------------------------------------------
    // Open file
    // -----------------------------------------------------------------------
    TFile* f = TFile::Open(fname);
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << fname << "\n";
        return;
    }
    std::cout << "\n=== File: " << fname << " ===\n";
    f->ls();

    // -----------------------------------------------------------------------
    // MuonSpectrometer tree (flat branches, no dictionary needed)
    // -----------------------------------------------------------------------
    TTree* ms = (TTree*)f->Get("MuonSpectrometer");
    if (!ms) { std::cerr << "No MuonSpectrometer tree!\n"; return; }

    const int MAXTR = 10;
    int   ntracks;
    float charge[MAXTR], p[MAXTR], px[MAXTR], py[MAXTR], pz[MAXTR];
    float chi2[MAXTR], pval[MAXTR], fpErr[MAXTR], fipErr[MAXTR];
    int   nDoF[MAXTR], npoints[MAXTR], fit_ok[MAXTR];
    float p_analytic[MAXTR];

    ms->SetBranchAddress("ntracks",    &ntracks);
    ms->SetBranchAddress("charge",     charge);
    ms->SetBranchAddress("p",          p);
    ms->SetBranchAddress("px",         px);
    ms->SetBranchAddress("py",         py);
    ms->SetBranchAddress("pz",         pz);
    ms->SetBranchAddress("chi2",       chi2);
    ms->SetBranchAddress("nDoF",       nDoF);
    ms->SetBranchAddress("pval",       pval);
    ms->SetBranchAddress("fpErr",      fpErr);
    ms->SetBranchAddress("fipErr",     fipErr);
    ms->SetBranchAddress("npoints",    npoints);

    // fit_ok and p_analytic were added in ClassDef v3; older files won't have them.
    bool has_fit_ok    = (ms->GetBranch("fit_ok")     != nullptr);
    bool has_panalytic = (ms->GetBranch("p_analytic") != nullptr);
    if (!has_fit_ok || !has_panalytic)
        std::cout << "\n[WARNING] Old file format (pre-ClassDef v3): "
                  << "fit_ok / p_analytic branches not present.\n"
                  << "  Re-run batchreco.exe on the input events to get new branches.\n"
                  << "  Convergence rate and p_analytic plots will be skipped.\n\n";
    if (has_fit_ok)    ms->SetBranchAddress("fit_ok",     fit_ok);
    if (has_panalytic) ms->SetBranchAddress("p_analytic", p_analytic);

    // -----------------------------------------------------------------------
    // Histograms
    // -----------------------------------------------------------------------
    // pAnalytic — use log-uniform binning so sub-GeV and TeV-scale tracks both visible
    const int nLogBins = 60;
    double logEdges[nLogBins+1];
    for (int i = 0; i <= nLogBins; ++i)
        logEdges[i] = std::pow(10.0, -1.0 + 5.0 * i / nLogBins); // 0.1 GeV – 10 TeV
    TH1F* h_pAna_all  = new TH1F("h_pAna_all",  "p_{analytic} (all tracks);p_{analytic} [GeV/c];tracks", nLogBins, logEdges);
    TH1F* h_pAna_ok   = new TH1F("h_pAna_ok",   "p_{analytic} (fit_ok==1);p_{analytic} [GeV/c];tracks",  nLogBins, logEdges);
    TH1F* h_pAna_fail = new TH1F("h_pAna_fail",  "p_{analytic} (fit_ok==0);p_{analytic} [GeV/c];tracks", nLogBins, logEdges);
    TH1F* h_pAna_log  = new TH1F("h_pAna_log",  "p_{analytic} (log scale);log_{10}(p_{analytic}/GeV);tracks", 50, -1, 4);

    // GenFit fitted momentum (fit_ok==1 only)
    TH1F* h_pFit      = new TH1F("h_pFit",      "p_{fit} (fit_ok==1);p_{fit} [GeV/c];tracks", 50, 0, 2000);
    TH1F* h_pFit_log  = new TH1F("h_pFit_log",  "p_{fit} (log scale);log_{10}(p_{fit}/GeV);tracks", 50, -1, 4);

    // chi2 / nDoF  (fit_ok==1)
    TH1F* h_chi2ndf   = new TH1F("h_chi2ndf",   "#chi^{2}/nDoF (fit_ok==1);#chi^{2}/nDoF;tracks", 50, 0, 10);
    TH1F* h_nDoF      = new TH1F("h_nDoF",      "nDoF (fit_ok==1);nDoF;tracks",                    15, 0, 15);
    TH1F* h_pval      = new TH1F("h_pval",      "p-value (fit_ok==1);p-value;tracks",              25, 0, 1);

    // Charge
    TH1F* h_charge    = new TH1F("h_charge",    "Fitted charge (fit_ok==1);charge;tracks",          3,-1.5,1.5);

    // pFit vs pAnalytic 2D (fit_ok==1)
    TH2F* h2_paVpf    = new TH2F("h2_paVpf",
                                  "p_{analytic} vs p_{fit};p_{analytic} [GeV/c];p_{fit} [GeV/c]",
                                  40, 0, 600, 40, 0, 2000);

    // relative momentum error 1/p * sigma_1/p  (fit_ok==1)
    TH1F* h_relErr    = new TH1F("h_relErr",    "#sigma(1/p)/(1/p) (fit_ok==1);#sigma_{1/p}/(1/p);tracks", 50, 0, 2);

    // npoints per track (all)
    TH1F* h_npts      = new TH1F("h_npts",      "npoints per track (all);npoints;tracks", 20, 0, 20);

    // fit_ok per event
    TH1F* h_ntracks   = new TH1F("h_ntracks",   "ntracks per event;ntracks;events",        8, 0, 8);
    TH1F* h_nok       = new TH1F("h_nok",       "fit_ok==1 per event;n fit_ok tracks;events", 6, 0, 6);

    // -----------------------------------------------------------------------
    // Event loop
    // -----------------------------------------------------------------------
    Long64_t nev = ms->GetEntries();
    std::cout << "\nMuonSpectrometer entries: " << nev << "\n";

    int total_tracks = 0, total_ok = 0, total_fail = 0;
    int cnt_lt1 = 0, cnt_1_10 = 0, cnt_10_100 = 0, cnt_gt100 = 0;

    for (Long64_t ev = 0; ev < nev; ++ev) {
        ms->GetEntry(ev);
        h_ntracks->Fill(ntracks);
        int nok_ev = 0;
        for (int t = 0; t < ntracks && t < MAXTR; ++t) {
            ++total_tracks;
            float pa = p_analytic[t];
            if (pa > 0 && std::isfinite(pa))
                h_pAna_log->Fill(std::log10(pa));

            if (fit_ok[t] == 1) {
                ++total_ok;
                ++nok_ev;
                h_pAna_ok->Fill(pa);
                h_pFit->Fill(p[t]);
                if (p[t] > 0 && std::isfinite(p[t]))
                    h_pFit_log->Fill(std::log10(p[t]));
                if (nDoF[t] > 0)
                    h_chi2ndf->Fill(chi2[t] / nDoF[t]);
                h_nDoF->Fill(nDoF[t]);
                h_pval->Fill(pval[t]);
                h_charge->Fill(charge[t]);
                h2_paVpf->Fill(pa, p[t]);
                // relative inverse-p error:  sigma(1/p) / (1/p)  =  sigma(1/p) * p
                if (p[t] > 0 && fipErr[t] > 0)
                    h_relErr->Fill(fipErr[t] * p[t]);
            } else {
                ++total_fail;
                h_pAna_fail->Fill(pa);
            }
            h_pAna_all->Fill(pa);
            h_npts->Fill(npoints[t]);
            if (pa > 0 && std::isfinite(pa)) {
                if      (pa <   1.0) ++cnt_lt1;
                else if (pa <  10.0) ++cnt_1_10;
                else if (pa < 100.0) ++cnt_10_100;
                else                 ++cnt_gt100;
            }
        }
        h_nok->Fill(nok_ev);
    }

    // -----------------------------------------------------------------------
    // Summary printout
    // -----------------------------------------------------------------------
    std::cout << "\n=== MDT Reconstruction Summary (" << nev << " events) ===\n";
    std::cout << "  Total tracks reconstructed : " << total_tracks << "\n";
    std::cout << "  GenFit SUCCESS  (fit_ok==1): " << total_ok << "\n";
    std::cout << "  GenFit FAILED   (fit_ok==0): " << total_fail << "\n";
    if (total_tracks > 0)
        std::cout << "  Convergence rate           : "
                  << 100.0 * total_ok / total_tracks << " %\n";

    std::cout << "\n  p_analytic distribution (all " << total_tracks << " tracks):\n";
    std::cout << "    p < 1  GeV : " << cnt_lt1    << "\n";
    std::cout << "    1-10  GeV  : " << cnt_1_10   << "\n";
    std::cout << "    10-100 GeV : " << cnt_10_100 << "\n";
    std::cout << "    100+ GeV   : " << cnt_gt100  << "\n";

    if (total_ok > 0) {
        std::cout << "\n  Fitted momentum (fit_ok==1):\n";
        std::cout << "    mean p_fit  = " << h_pFit->GetMean()   << " GeV/c\n";
        std::cout << "    mean chi2   = " << h_chi2ndf->GetMean() << " /nDoF\n";
        std::cout << "    mean pval   = " << h_pval->GetMean()   << "\n";
        std::cout << "    mean relErr = " << h_relErr->GetMean() << "\n";
    }

    // -----------------------------------------------------------------------
    // Plotting
    // -----------------------------------------------------------------------
    gStyle->SetPalette(kBird);

    // Canvas 1: pAnalytic
    TCanvas* c1 = new TCanvas("c1", "MDT – Analytic Momentum", 1200, 500);
    c1->Divide(2, 1);

    c1->cd(1);
    h_pAna_all->SetLineColor(kBlue+1);
    gPad->SetLogx();
    h_pAna_ok->SetLineColor(kGreen+2);
    h_pAna_ok->SetFillColor(kGreen-9);
    h_pAna_fail->SetLineColor(kRed+1);
    h_pAna_all->Draw("hist");
    h_pAna_ok->Draw("hist same");
    h_pAna_fail->Draw("hist same");
    TLegend* leg1 = new TLegend(0.55,0.65,0.88,0.88);
    leg1->SetBorderSize(0);
    leg1->AddEntry(h_pAna_all,  "All tracks",     "l");
    leg1->AddEntry(h_pAna_ok,   "fit_ok == 1",    "lf");
    leg1->AddEntry(h_pAna_fail, "fit_ok == 0",    "l");
    leg1->Draw();

    c1->cd(2);
    gPad->SetLogy();
    h_pAna_log->SetLineColor(kBlue+1);
    h_pAna_log->Draw("hist");
    TLatex lat;
    lat.SetNDC(); lat.SetTextSize(0.04);
    lat.DrawLatex(0.15, 0.85, Form("All %d tracks", total_tracks));

    c1->SaveAs("mdt_panalytic.png");
    std::cout << "\nSaved: mdt_panalytic.png\n";

    // Canvas 2: Fitted momentum and fit quality (fit_ok==1 only)
    TCanvas* c2 = new TCanvas("c2", "MDT – GenFit Results (fit_ok==1)", 1400, 600);
    c2->Divide(3, 2);

    c2->cd(1); gPad->SetLogy();
    h_pFit->SetLineColor(kBlue+1); h_pFit->Draw("hist");
    lat.DrawLatex(0.15, 0.85, Form("%d tracks", total_ok));

    c2->cd(2); gPad->SetLogy();
    h_pFit_log->SetLineColor(kBlue+1); h_pFit_log->Draw("hist");

    c2->cd(3);
    h_chi2ndf->SetLineColor(kRed+1); h_chi2ndf->Draw("hist");
    TLine* lchi = new TLine(1, 0, 1, h_chi2ndf->GetMaximum());
    lchi->SetLineStyle(2); lchi->SetLineColor(kGray+1); lchi->Draw();

    c2->cd(4);
    h_nDoF->SetLineColor(kMagenta+1); h_nDoF->Draw("hist");

    c2->cd(5);
    h_pval->SetLineColor(kGreen+2); h_pval->Draw("hist");

    c2->cd(6);
    h_charge->SetLineColor(kOrange+1); h_charge->Draw("hist");

    c2->SaveAs("mdt_fitquality.png");
    std::cout << "Saved: mdt_fitquality.png\n";

    // Canvas 3: pAnalytic vs pFit 2D + relative error
    TCanvas* c3 = new TCanvas("c3", "MDT – pAnalytic vs pFit", 1200, 500);
    c3->Divide(2, 1);

    c3->cd(1);
    gPad->SetLogx(); gPad->SetLogy();
    h2_paVpf->Draw("colz");

    c3->cd(2);
    h_relErr->SetLineColor(kBlue+1); h_relErr->Draw("hist");
    lat.DrawLatex(0.15, 0.85, Form("mean = %.2f", h_relErr->GetMean()));

    c3->SaveAs("mdt_pana_vs_pfit.png");
    std::cout << "Saved: mdt_pana_vs_pfit.png\n";

    // Canvas 4: per-event counters
    TCanvas* c4 = new TCanvas("c4", "MDT – Per-event", 1000, 450);
    c4->Divide(2, 1);

    c4->cd(1);
    h_ntracks->SetLineColor(kBlue+1); h_ntracks->Draw("hist");

    c4->cd(2);
    h_nok->SetLineColor(kGreen+2); h_nok->Draw("hist");
    TLatex lat2;
    lat2.SetNDC(); lat2.SetTextSize(0.04);
    if (total_tracks > 0)
        lat2.DrawLatex(0.15, 0.85,
            Form("Conv. rate: %.1f%%", 100.0*total_ok/total_tracks));

    c4->SaveAs("mdt_per_event.png");
    std::cout << "Saved: mdt_per_event.png\n";

    std::cout << "\nDone.\n";
}
