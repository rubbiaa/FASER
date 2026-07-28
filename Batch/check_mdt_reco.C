// check_mdt_reco.C
// Usage:  root -l -b -q 'check_mdt_reco.C("Batch-TPORecevent_6000_0_5210.root")'
//
// Prints MDT reconstruction statistics and plots truth vs reco momentum resolution.
// Requires the p_truth branch (added after Jul 2026 rebuild).

void check_mdt_reco(const char* fname = "Batch-TPORecevent_6000_0_5210.root")
{
    TFile* f = TFile::Open(fname);
    if (!f || f->IsZombie()) { std::cerr << "Cannot open " << fname << "\n"; return; }

    TTree* t = (TTree*)f->Get("MuonSpectrometer");
    if (!t) { std::cerr << "TTree 'MuonSpectrometer' not found\n"; return; }

    const bool hasTruth = (t->GetBranch("p_truth") != nullptr);
    if (!hasTruth) printf("[WARNING] p_truth branch not found — truth resolution plots disabled.\n");

    const Long64_t nEntries = t->GetEntries();

    // --- branch variables ---
    const int MAXTRK = 10;
    int   ntracks = 0;
    int   fit_ok [MAXTRK] = {};
    int   nDoF   [MAXTRK] = {};
    float p      [MAXTRK] = {};
    float p_analytic[MAXTRK] = {};
    float p_truth   [MAXTRK] = {};
    float chi2   [MAXTRK] = {};
    float charge [MAXTRK] = {};
    float fpErr  [MAXTRK] = {};

    t->SetBranchAddress("ntracks",     &ntracks);
    t->SetBranchAddress("fit_ok",       fit_ok);
    t->SetBranchAddress("nDoF",         nDoF);
    t->SetBranchAddress("p",            p);
    t->SetBranchAddress("p_analytic",   p_analytic);
    t->SetBranchAddress("chi2",         chi2);
    t->SetBranchAddress("charge",       charge);
    t->SetBranchAddress("fpErr",        fpErr);
    if (hasTruth) t->SetBranchAddress("p_truth", p_truth);

    // --- counters ---
    int nEventsTotal     = 0;
    int nEventsNoTrack   = 0;
    int nEventsAttempted = 0;
    int nFitFailed       = 0;
    int nSecondTrack     = 0;
    int nNDF0            = 0;

    std::map<int,int> ndfHist;

    // momentum data vectors (fit_ok, physics range)
    // "all"  = all fit_ok, NDF>0 tracks in physics range
    // "good" = additionally require chi2/NDF < chi2Cut (well-reconstructed)
    const double chi2Cut = 10.0;
    std::vector<double> pFit, pAna, pTru;           // all
    std::vector<double> pFitG, pAnaG, pTruG;        // good chi2
    std::vector<double> dp_over_p_genfit,   dp_over_p_genfit_g;
    std::vector<double> dp_over_p_analytic, dp_over_p_analytic_g;
    std::vector<double> dinvp_genfit,       dinvp_genfit_g;
    std::vector<double> dinvp_analytic,     dinvp_analytic_g;
    std::vector<double> chi2ndfVec;

    const float PMAX = 1e4f; // upper sanity cut

    for (Long64_t ev = 0; ev < nEntries; ++ev) {
        t->GetEntry(ev);
        ++nEventsTotal;

        if (ntracks == 0) { ++nEventsNoTrack; continue; }
        ++nEventsAttempted;
        if (ntracks >= 2) ++nSecondTrack;

        bool anyOk = false;
        for (int tr = 0; tr < std::min(ntracks, MAXTRK); ++tr) {
            if (!fit_ok[tr]) continue;
            anyOk = true;

            int ndf = nDoF[tr];
            ndfHist[ndf]++;
            if (ndf == 0) { ++nNDF0; continue; }

            const double cndf = (ndf > 0) ? chi2[tr] / ndf : 1e9;
            if (cndf < 1e6) chi2ndfVec.push_back(cndf);
            const bool goodChi2 = (cndf < chi2Cut);

            const float pf = p[tr];
            const float pa = p_analytic[tr];
            const float pt = hasTruth ? p_truth[tr] : -999.f;

            // GenFit vs truth
            if (hasTruth
             && std::isfinite(pf) && pf > 0 && pf < PMAX
             && std::isfinite(pt) && pt > 0 && pt < PMAX) {
                pFit.push_back(pf);
                pTru.push_back(pt);
                dp_over_p_genfit.push_back((pf - pt) / pt);
                dinvp_genfit.push_back(1.0/pf - 1.0/pt);
                if (goodChi2) {
                    pFitG.push_back(pf); pTruG.push_back(pt);
                    dp_over_p_genfit_g.push_back((pf - pt) / pt);
                    dinvp_genfit_g.push_back(1.0/pf - 1.0/pt);
                }
            }
            // Analytic vs truth
            if (hasTruth
             && std::isfinite(pa) && pa > 0 && pa < PMAX
             && std::isfinite(pt) && pt > 0 && pt < PMAX) {
                pAna.push_back(pa);
                dp_over_p_analytic.push_back((pa - pt) / pt);
                dinvp_analytic.push_back(1.0/pa - 1.0/pt);
                if (goodChi2) {
                    pAnaG.push_back(pa);
                    dp_over_p_analytic_g.push_back((pa - pt) / pt);
                    dinvp_analytic_g.push_back(1.0/pa - 1.0/pt);
                }
            }
        }
        if (!anyOk) ++nFitFailed;
    }

    // compute mean and RMS of a vector
    auto vecMean = [](const std::vector<double>& v) -> double {
        if (v.empty()) return 0;
        double s = 0; for (double x : v) s += x; return s / v.size();
    };
    auto vecRMS = [&vecMean](const std::vector<double>& v) -> double {
        if (v.empty()) return 0;
        double mu = vecMean(v), s = 0;
        for (double x : v) s += (x-mu)*(x-mu);
        return std::sqrt(s / v.size());
    };

    double mu_dpp_gf      = vecMean(dp_over_p_genfit),    sig_dpp_gf      = vecRMS(dp_over_p_genfit);
    double mu_dpp_ana     = vecMean(dp_over_p_analytic),  sig_dpp_ana     = vecRMS(dp_over_p_analytic);
    double mu_dinvp_gf    = vecMean(dinvp_genfit),        sig_dinvp_gf    = vecRMS(dinvp_genfit);
    double mu_dinvp_ana   = vecMean(dinvp_analytic),      sig_dinvp_ana   = vecRMS(dinvp_analytic);
    double mu_dpp_gf_g    = vecMean(dp_over_p_genfit_g),  sig_dpp_gf_g    = vecRMS(dp_over_p_genfit_g);
    double mu_dpp_ana_g   = vecMean(dp_over_p_analytic_g),sig_dpp_ana_g   = vecRMS(dp_over_p_analytic_g);
    double mu_dinvp_gf_g  = vecMean(dinvp_genfit_g),      sig_dinvp_gf_g  = vecRMS(dinvp_genfit_g);
    double mu_dinvp_ana_g = vecMean(dinvp_analytic_g),    sig_dinvp_ana_g = vecRMS(dinvp_analytic_g);
    double mu_chi2        = vecMean(chi2ndfVec);

    // --- Print summary ---
    printf("\n========================================\n");
    printf("  MDT Reco Statistics: %s\n", fname);
    printf("========================================\n");
    printf("  Total events                  : %d\n", nEventsTotal);
    printf("  Events with no MDT track      : %d\n", nEventsNoTrack);
    printf("  Events with MDT track(s)      : %d  (%.1f%%)\n",
           nEventsAttempted, 100.0*nEventsAttempted/nEventsTotal);
    printf("  Events with 2nd track (B)     : %d\n", nSecondTrack);
    printf("\n");
    printf("  Fit SUCCESS events            : %d  (%.1f%% of attempted)\n",
           nEventsAttempted - nFitFailed,
           nEventsAttempted ? 100.0*(nEventsAttempted-nFitFailed)/nEventsAttempted : 0.0);
    printf("  Fit FAILED  events            : %d\n", nFitFailed);
    printf("\n");
    printf("  NDF distribution (fit_ok tracks):\n");
    for (auto& kv : ndfHist)
        printf("    NDF = %3d : %d tracks\n", kv.first, kv.second);
    printf("  fit_ok tracks with NDF=0      : %d  (DAF garbage)\n", nNDF0);
    printf("  Mean chi2/NDF (NDF>0)         : %.2f\n", mu_chi2);
    printf("\n");
    if (hasTruth) {
        printf("  --- GenFit vs truth  ALL (N=%zu) ---\n", pFit.size());
        printf("  (p_fit-p_truth)/p_truth    : mean=%+.3f  RMS=%.3f\n", mu_dpp_gf, sig_dpp_gf);
        printf("  1/p_fit-1/p_truth [1/GeV]  : mean=%+.3e  RMS=%.3e\n", mu_dinvp_gf, sig_dinvp_gf);
        printf("  --- GenFit vs truth  chi2/NDF<%.0f (N=%zu) ---\n", chi2Cut, pFitG.size());
        printf("  (p_fit-p_truth)/p_truth    : mean=%+.3f  RMS=%.3f\n", mu_dpp_gf_g, sig_dpp_gf_g);
        printf("  1/p_fit-1/p_truth [1/GeV]  : mean=%+.3e  RMS=%.3e\n", mu_dinvp_gf_g, sig_dinvp_gf_g);
        printf("\n");
        printf("  --- Analytic vs truth  ALL (N=%zu) ---\n", pAna.size());
        printf("  (p_ana-p_truth)/p_truth    : mean=%+.3f  RMS=%.3f\n", mu_dpp_ana, sig_dpp_ana);
        printf("  1/p_ana-1/p_truth [1/GeV]  : mean=%+.3e  RMS=%.3e\n", mu_dinvp_ana, sig_dinvp_ana);
        printf("  --- Analytic vs truth  chi2/NDF<%.0f (N=%zu) ---\n", chi2Cut, pAnaG.size());
        printf("  (p_ana-p_truth)/p_truth    : mean=%+.3f  RMS=%.3f\n", mu_dpp_ana_g, sig_dpp_ana_g);
        printf("  1/p_ana-1/p_truth [1/GeV]  : mean=%+.3e  RMS=%.3e\n", mu_dinvp_ana_g, sig_dinvp_ana_g);
    }
    printf("========================================\n\n");

    // --- Histograms ---
    TCanvas* c = new TCanvas("c_mdt","MDT Reco Check",1800,1000);
    int ncols = hasTruth ? 4 : 3;
    c->Divide(ncols, 2);
    int pad = 1;

    // 1) NDF
    c->cd(pad++);
    TH1D* hNDF = new TH1D("hNDF","NDF of fit_ok tracks;NDF;Tracks",25,-0.5,24.5);
    for (auto& kv : ndfHist)
        for (int i = 0; i < kv.second; ++i) hNDF->Fill(kv.first);
    hNDF->Draw();

    // 2) chi2/NDF
    c->cd(pad++);
    TH1D* hChi2 = new TH1D("hChi2","#chi^{2}/NDF (fit_ok, NDF>0);#chi^{2}/NDF;Tracks",60,0,60);
    for (double v : chi2ndfVec) hChi2->Fill(v);
    hChi2->Draw();

    // 3) Fitted p spectrum — log x
    c->cd(pad++);
    gPad->SetLogx();
    const int NBINS = 50;
    double logBins[NBINS+1];
    for (int i = 0; i <= NBINS; ++i) logBins[i] = std::pow(10.0, 0.0 + i * 4.0/NBINS); // 1..10000 GeV
    TH1D* hP = new TH1D("hP","Fitted |p| (fit_ok, NDF>0);p [GeV/c];Tracks",NBINS,logBins);
    for (double v : pFit) hP->Fill(v);
    hP->Draw();

    if (hasTruth) {
        // 4) p_fit vs p_truth 2D — log-log
        c->cd(pad++);
        gPad->SetLogx(); gPad->SetLogy();
        double logBins2[NBINS+1];
        for (int i = 0; i <= NBINS; ++i) logBins2[i] = std::pow(10.0, 0.0 + i * 4.0/NBINS);
        TH2D* h2P = new TH2D("h2P","p_{fit} vs p_{truth};p_{truth} [GeV/c];p_{fit} [GeV/c]",
                              NBINS,logBins2, NBINS,logBins2);
        for (size_t i = 0; i < pFit.size(); ++i) h2P->Fill(pTru[i], pFit[i]);
        h2P->Draw("COLZ");
        TLine* diag = new TLine(1,1,1e4,1e4);
        diag->SetLineColor(kRed); diag->SetLineWidth(2); diag->Draw();

        // 5) GenFit (p_fit-p_truth)/p_truth — all vs good chi2
        c->cd(pad++);
        TH1D* hDp  = new TH1D("hDp_gf_all","GenFit #Deltap/p (all);(p_{fit}-p_{truth})/p_{truth};Tracks",60,-3,3);
        TH1D* hDpG = new TH1D("hDp_gf_good",Form("#chi^{2}/NDF<%.0f",chi2Cut),60,-3,3);
        for (double v : dp_over_p_genfit)   hDp->Fill(v);
        for (double v : dp_over_p_genfit_g) hDpG->Fill(v);
        hDp->SetLineColor(kBlue-7);  hDp->Draw();
        hDpG->SetLineColor(kBlue+2); hDpG->SetLineWidth(2); hDpG->Draw("SAME");
        gPad->BuildLegend(0.55,0.7,0.98,0.92);

        // 6) Analytic (p_ana-p_truth)/p_truth — all vs good chi2
        c->cd(pad++);
        TH1D* hDpA  = new TH1D("hDp_ana_all","Analytic #Deltap/p (all);(p_{ana}-p_{truth})/p_{truth};Tracks",60,-3,3);
        TH1D* hDpAG = new TH1D("hDp_ana_good",Form("#chi^{2}/NDF<%.0f",chi2Cut),60,-3,3);
        for (double v : dp_over_p_analytic)   hDpA->Fill(v);
        for (double v : dp_over_p_analytic_g) hDpAG->Fill(v);
        hDpA->SetLineColor(kGreen-7);  hDpA->Draw();
        hDpAG->SetLineColor(kGreen+3); hDpAG->SetLineWidth(2); hDpAG->Draw("SAME");
        gPad->BuildLegend(0.55,0.7,0.98,0.92);

        // 7) GenFit 1/p curvature resolution — all vs good chi2
        c->cd(pad++);
        TH1D* hIp  = new TH1D("hInvP_gf_all","GenFit #Delta(1/p) (all);1/p_{fit}-1/p_{truth} [GeV^{-1}];Tracks",60,-0.1,0.1);
        TH1D* hIpG = new TH1D("hInvP_gf_good",Form("#chi^{2}/NDF<%.0f",chi2Cut),60,-0.1,0.1);
        for (double v : dinvp_genfit)   hIp->Fill(v);
        for (double v : dinvp_genfit_g) hIpG->Fill(v);
        hIp->SetLineColor(kBlue-7);  hIp->Draw();
        hIpG->SetLineColor(kBlue+2); hIpG->SetLineWidth(2); hIpG->Draw("SAME");
        gPad->BuildLegend(0.55,0.7,0.98,0.92);

        // 8) Analytic 1/p curvature resolution — all vs good chi2
        c->cd(pad++);
        TH1D* hIpA  = new TH1D("hInvP_ana_all","Analytic #Delta(1/p) (all);1/p_{ana}-1/p_{truth} [GeV^{-1}];Tracks",60,-0.1,0.1);
        TH1D* hIpAG = new TH1D("hInvP_ana_good",Form("#chi^{2}/NDF<%.0f",chi2Cut),60,-0.1,0.1);
        for (double v : dinvp_analytic)   hIpA->Fill(v);
        for (double v : dinvp_analytic_g) hIpAG->Fill(v);
        hIpA->SetLineColor(kGreen-7);  hIpA->Draw();
        hIpAG->SetLineColor(kGreen+3); hIpAG->SetLineWidth(2); hIpAG->Draw("SAME");
        gPad->BuildLegend(0.55,0.7,0.98,0.92);
    } else {
        // fallback without truth: p_fit vs p_analytic scatter
        c->cd(pad++);
        std::vector<double> pAnaFallback;
        for (Long64_t ev = 0; ev < nEntries; ++ev) {
            t->GetEntry(ev);
            for (int tr = 0; tr < std::min(ntracks,MAXTRK); ++tr) {
                if (!fit_ok[tr] || nDoF[tr]<=0) continue;
                if (std::isfinite(p[tr]) && p[tr]>0 && p[tr]<PMAX
                 && std::isfinite(p_analytic[tr]) && p_analytic[tr]>0 && p_analytic[tr]<PMAX)
                    pAnaFallback.push_back(p_analytic[tr]);
            }
        }
        TH1D* hPA = new TH1D("hPA","Analytic |p|;p_{analytic} [GeV/c];Tracks",50,0,500);
        for (double v : pAnaFallback) hPA->Fill(v);
        hPA->Draw();
    }

    c->SaveAs("mdt_reco_check.pdf");
    printf("  Plots saved to mdt_reco_check.pdf\n\n");
}
