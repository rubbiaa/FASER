// plot_mdt_validation.C
// Produces three validation plots from a MuonSpectrometer TTree (ReconstructMDT_fin output):
//   1. qp_scatter.png   -- q*p_truth vs q*p_fit  (signed momentum, tests charge+momentum together)
//   2. p_truth_reco.png -- p_truth vs p_fit      (unsigned momentum resolution)
//   3. efficiency.png   -- fit-success rate and charge-accuracy rate vs p_truth, binned
//
// Usage: root -l -b -q 'plot_mdt_validation.C("Batch-TPORecevent_6000_0_5210.root")'

void plot_mdt_validation(const char* filename, const char* outPrefix = "mdtval") {
    TFile* f = TFile::Open(filename);
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", filename); return; }
    TTree* t = (TTree*)f->Get("MuonSpectrometer");
    if (!t) { printf("No MuonSpectrometer tree\n"); return; }

    const int MAXTRK = 10;
    int   ntracks;
    int   fit_ok[MAXTRK] = {};
    int   nDoF[MAXTRK] = {};
    float p[MAXTRK] = {};
    float p_truth[MAXTRK] = {};
    float chi2[MAXTRK] = {};
    float charge[MAXTRK] = {};
    float charge_truth[MAXTRK] = {};

    t->SetBranchAddress("ntracks", &ntracks);
    t->SetBranchAddress("fit_ok", fit_ok);
    t->SetBranchAddress("nDoF", nDoF);
    t->SetBranchAddress("p", p);
    t->SetBranchAddress("p_truth", p_truth);
    t->SetBranchAddress("chi2", chi2);
    t->SetBranchAddress("charge", charge);
    t->SetBranchAddress("charge_truth", charge_truth);

    gStyle->SetOptStat(0);

    // ---- Plot 1: signed q*p scatter ----
    TGraph* gQP = new TGraph();
    int nQP = 0;

    // ---- Plot 2: unsigned p truth vs reco ----
    TGraph* gP = new TGraph();
    int nP = 0;

    // ---- Plot 4: signed qp distributions, truth vs reco, overlaid ----
    const double qpHistMax = 2000.0;
    TH1D* hQPTruth = new TH1D("hQPTruth", "", 100, -qpHistMax, qpHistMax);
    TH1D* hQPFit   = new TH1D("hQPFit",   "", 100, -qpHistMax, qpHistMax);

    // ---- Plot 3: efficiency vs p_truth, log-spaced bins ----
    const int NBINS_P = 8;
    double pEdges[NBINS_P+1];
    for (int i = 0; i <= NBINS_P; ++i) pEdges[i] = 5.0 * std::pow(1000.0, i/(double)NBINS_P);
    int nAttempt_p[NBINS_P] = {}; // tracks with truth momentum in bin (denominator for fit efficiency, from all fit_ok+truth attempts)
    int nFitOk_p[NBINS_P] = {};
    int nAssigned_p[NBINS_P] = {};
    int nCorrect_p[NBINS_P] = {};

    Long64_t nentries = t->GetEntries();
    for (Long64_t ev = 0; ev < nentries; ++ev) {
        t->GetEntry(ev);
        for (int tr = 0; tr < std::min(ntracks, MAXTRK); ++tr) {
            if (p_truth[tr] <= 0 || !std::isfinite(p_truth[tr])) continue;
            if (std::fabs(charge_truth[tr]) < 0.5) continue;

            int pbin = -1;
            for (int b = 0; b < NBINS_P; ++b) if (p_truth[tr] >= pEdges[b] && p_truth[tr] < pEdges[b+1]) { pbin = b; break; }
            if (pbin >= 0) nAttempt_p[pbin]++;

            if (!fit_ok[tr] || nDoF[tr] <= 0) continue;
            if (pbin >= 0) nFitOk_p[pbin]++;

            const double cndf = chi2[tr] / nDoF[tr];

            gP->SetPoint(nP++, p_truth[tr], p[tr]);

            // Only plot tracks with an actual assigned charge -- charge=0
            // (ambiguous, Mode 0) tracks would otherwise show up as a
            // spurious horizontal qpFit=0 band with no physical meaning.
            if (std::fabs(charge[tr]) > 0.5) {
                const double qpTruth = charge_truth[tr] * p_truth[tr];
                const double qpFit   = charge[tr] * p[tr];
                gQP->SetPoint(nQP++, qpTruth, qpFit);
                if (std::fabs(qpFit) > 2000.0) {
                    printf("  [QP OUTLIER] event=%lld track=%d qpTruth=%.2f qpFit=%.2f p_fit=%.2f p_truth=%.2f chi2/ndf=%.3f nDoF=%d\n",
                           ev, tr, qpTruth, qpFit, p[tr], p_truth[tr], cndf, nDoF[tr]);
                }
                // Same track population as the qp scatter above (assigned charge only),
                // so the two plots are directly comparable.
                hQPTruth->Fill(qpTruth);
                hQPFit->Fill(qpFit);
            }

            if (std::fabs(charge[tr]) > 0.5) {
                if (pbin >= 0) nAssigned_p[pbin]++;
                if (charge[tr] * charge_truth[tr] > 0 && pbin >= 0) nCorrect_p[pbin]++;
            }
        }
    }

    // ---- Draw Plot 1: qp scatter ----
    {
        TCanvas* c1 = new TCanvas("c1", "qp", 800, 800);
        c1->SetGrid();
        gQP->SetMarkerStyle(20);
        gQP->SetMarkerSize(0.6);
        gQP->SetMarkerColor(kAzure+2);
        gQP->SetTitle("Signed momentum: truth vs reco;q_{truth} #times p_{truth}  [GeV/c];q_{fit} #times p_{fit}  [GeV/c]");
        gQP->Draw("AP");
        // Fixed, physics-motivated range rather than autoscaling to the max point:
        // a handful of very-high-momentum (low-curvature, hardest to charge-ID)
        // outliers would otherwise compress the whole well-behaved bulk into an
        // unreadable cluster at the origin. Those outliers are listed separately
        // below (see [QP OUTLIER] lines) rather than hidden by rescaling around them.
        const double qpMax = 2000.0;
        gQP->GetXaxis()->SetLimits(-qpMax, qpMax);
        gQP->SetMinimum(-qpMax);
        gQP->SetMaximum(qpMax);
        TLine* diag = new TLine(-qpMax, -qpMax, qpMax, qpMax);
        diag->SetLineColor(kRed);
        diag->SetLineStyle(2);
        diag->Draw();
        c1->SaveAs(Form("%s_qp_scatter.png", outPrefix));
    }

    // ---- Draw Plot 2: p truth vs reco ----
    {
        TCanvas* c2 = new TCanvas("c2", "p", 800, 800);
        c2->SetGrid();
        c2->SetLogx(); c2->SetLogy();
        gP->SetMarkerStyle(20);
        gP->SetMarkerSize(0.6);
        gP->SetMarkerColor(kAzure+2);
        gP->SetTitle("Unsigned momentum: truth vs reco;p_{truth}  [GeV/c];p_{fit}  [GeV/c]");
        gP->Draw("AP");
        gP->GetXaxis()->SetLimits(3, 6000);
        gP->SetMinimum(3);
        gP->SetMaximum(6000);
        TLine* diag2 = new TLine(3, 3, 6000, 6000);
        diag2->SetLineColor(kRed);
        diag2->SetLineStyle(2);
        diag2->Draw();
        c2->SaveAs(Form("%s_p_truth_reco.png", outPrefix));
    }

    // ---- Draw Plot 4: qp distributions, truth vs reco overlaid ----
    {
        TCanvas* c4 = new TCanvas("c4", "qp_hist", 900, 650);
        c4->SetGrid();
        hQPTruth->SetLineColor(kRed+1);
        hQPTruth->SetLineWidth(2);
        hQPTruth->SetTitle("Signed momentum distribution: truth vs reco;q #times p  [GeV/c];tracks");
        hQPFit->SetLineColor(kAzure+2);
        hQPFit->SetLineWidth(2);
        double ymax = std::max(hQPTruth->GetMaximum(), hQPFit->GetMaximum());
        hQPTruth->SetMaximum(1.15 * ymax);
        hQPTruth->Draw("HIST");
        hQPFit->Draw("HIST SAME");
        TLegend* legQP = new TLegend(0.62, 0.72, 0.88, 0.88);
        legQP->AddEntry(hQPTruth, "truth (q_{truth} #times p_{truth})", "l");
        legQP->AddEntry(hQPFit, "reco (q_{fit} #times p_{fit})", "l");
        legQP->SetBorderSize(0);
        legQP->Draw();
        c4->SaveAs(Form("%s_qp_hist.png", outPrefix));
    }

    // ---- Draw Plot 3: efficiency vs p_truth ----
    {
        TCanvas* c3 = new TCanvas("c3", "eff", 900, 650);
        c3->SetGrid();
        c3->SetLogx();
        TH1D* hFitEff = new TH1D("hFitEff", "MDT reconstruction performance vs truth momentum;p_{truth} [GeV/c];efficiency", NBINS_P, pEdges);
        TH1D* hChargeAcc = new TH1D("hChargeAcc", "", NBINS_P, pEdges);
        for (int b = 0; b < NBINS_P; ++b) {
            double effFit = (nAttempt_p[b] > 0) ? nFitOk_p[b] / (double)nAttempt_p[b] : 0;
            double effFitErr = (nAttempt_p[b] > 0) ? std::sqrt(effFit*(1-effFit)/nAttempt_p[b]) : 0;
            hFitEff->SetBinContent(b+1, effFit);
            hFitEff->SetBinError(b+1, effFitErr);
            double accQ = (nAssigned_p[b] > 0) ? nCorrect_p[b] / (double)nAssigned_p[b] : 0;
            double accQErr = (nAssigned_p[b] > 0) ? std::sqrt(accQ*(1-accQ)/nAssigned_p[b]) : 0;
            hChargeAcc->SetBinContent(b+1, accQ);
            hChargeAcc->SetBinError(b+1, accQErr);
        }
        hFitEff->SetLineColor(kAzure+2);
        hFitEff->SetMarkerColor(kAzure+2);
        hFitEff->SetMarkerStyle(20);
        hFitEff->SetLineWidth(2);
        hFitEff->SetMinimum(0);
        hFitEff->SetMaximum(1.1);
        hFitEff->Draw("E1");
        hChargeAcc->SetLineColor(kGreen+2);
        hChargeAcc->SetMarkerColor(kGreen+2);
        hChargeAcc->SetMarkerStyle(21);
        hChargeAcc->SetLineWidth(2);
        hChargeAcc->Draw("E1 SAME");
        TLegend* leg = new TLegend(0.55, 0.15, 0.88, 0.32);
        leg->AddEntry(hFitEff, "fit success (of MDT tracks)", "lep");
        leg->AddEntry(hChargeAcc, "charge accuracy (of assigned)", "lep");
        leg->SetBorderSize(0);
        leg->Draw();
        c3->SaveAs(Form("%s_efficiency.png", outPrefix));
    }

    printf("\nSummary: nQP=%d nP=%d\n", nQP, nP);
    for (int b = 0; b < NBINS_P; ++b) {
        printf("  p_truth [%7.1f,%7.1f)  attempt=%4d fitOk=%4d(%5.1f%%)  assigned=%4d correct=%4d(%5.1f%%)\n",
               pEdges[b], pEdges[b+1], nAttempt_p[b], nFitOk_p[b],
               nAttempt_p[b] ? 100.0*nFitOk_p[b]/nAttempt_p[b] : 0.0,
               nAssigned_p[b], nCorrect_p[b],
               nAssigned_p[b] ? 100.0*nCorrect_p[b]/nAssigned_p[b] : 0.0);
    }
}
