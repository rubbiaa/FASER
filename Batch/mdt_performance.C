// mdt_performance.C
// Momentum + charge resolution vs truth momentum and vs number of hits,
// plus truth-vs-reco comparisons (momentum and charge).
// Requires the p_truth and charge_truth branches (charge_truth added Jul 2026).
// Usage: root -l -b -q 'mdt_performance.C("Batch-TPORecevent_6000_0_5210.root")'

void mdt_performance(const char* fname = "Batch-TPORecevent_6000_0_5210.root",
                     const char* outpdf = "mdt_performance.pdf")
{
    TFile* f = TFile::Open(fname);
    if (!f || f->IsZombie()) { std::cerr << "Cannot open " << fname << "\n"; return; }
    TTree* t = (TTree*)f->Get("MuonSpectrometer");
    if (!t) { std::cerr << "TTree 'MuonSpectrometer' not found\n"; return; }

    const bool hasTruth  = (t->GetBranch("p_truth") != nullptr);
    const bool hasQTruth = (t->GetBranch("charge_truth") != nullptr);
    if (!hasTruth)  { std::cerr << "[FATAL] p_truth branch not found\n"; return; }
    if (!hasQTruth) { std::cerr << "[WARNING] charge_truth branch not found — charge plots disabled\n"; }

    const int M = 10;
    int   ntracks = 0, fit_ok[M]={}, nDoF[M]={}, npoints[M]={};
    float p[M]={}, p_analytic[M]={}, p_truth[M]={}, chi2[M]={}, charge[M]={}, charge_truth[M]={};

    t->SetBranchAddress("ntracks",      &ntracks);
    t->SetBranchAddress("fit_ok",        fit_ok);
    t->SetBranchAddress("nDoF",          nDoF);
    t->SetBranchAddress("npoints",       npoints);
    t->SetBranchAddress("p",             p);
    t->SetBranchAddress("p_analytic",    p_analytic);
    t->SetBranchAddress("chi2",          chi2);
    t->SetBranchAddress("charge",        charge);
    t->SetBranchAddress("p_truth",       p_truth);
    if (hasQTruth) t->SetBranchAddress("charge_truth", charge_truth);

    const float PMAX     = 1e4f;
    const double chi2Cut = 10.0;

    // ── log-p bins: 5 → 5000 GeV, 8 bins (coarser than mdt_resolution.C so
    // charge mis-ID counts per bin stay statistically meaningful) ───────────
    const int NBINS_P = 8;
    double pEdges[NBINS_P+1];
    for (int i = 0; i <= NBINS_P; ++i)
        pEdges[i] = std::pow(10.0, std::log10(5.0) + i*(std::log10(5000.0)-std::log10(5.0))/NBINS_P);

    // ── nhits bins: 6..12 (2..4 stations x 3 planes) ────────────────────────
    const int NHITS_MIN = 6, NHITS_MAX = 12;
    const int NBINS_H = NHITS_MAX - NHITS_MIN + 1;

    // per-p-bin accumulators
    std::vector<double> dinvp_bin[NBINS_P], dpp_bin[NBINS_P];
    int   nTot_p[NBINS_P]={}, nMis_p[NBINS_P]={};
    // per-nhits-bin accumulators
    std::vector<double> dinvp_hbin[NBINS_H], dpp_hbin[NBINS_H];
    int   nTot_h[NBINS_H]={}, nMis_h[NBINS_H]={};

    // global vectors for truth-vs-reco plots
    std::vector<double> pFit, pTru, dpp_all, dinvp_all;
    int nQTot=0, nQMis=0;

    TH2D* h2Charge = new TH2D("h2Charge","Charge: reco vs truth;q_{truth};q_{reco}",
                               3,-1.5,1.5, 3,-1.5,1.5);
    TH2D* h2P = new TH2D("h2P","p_{fit} vs p_{truth} (fit_ok, #chi^{2}/NDF<%.0f cut);p_{truth} [GeV/c];p_{fit} [GeV/c]",
                          50,std::log10(5.0),std::log10(5000.0), 50,std::log10(5.0),std::log10(5000.0));

    for (Long64_t ev = 0; ev < t->GetEntries(); ++ev) {
        t->GetEntry(ev);
        for (int tr = 0; tr < std::min(ntracks, M); ++tr) {
            if (!fit_ok[tr] || nDoF[tr] <= 0) continue;

            const double cndf = chi2[tr]/nDoF[tr];
            if (cndf >= chi2Cut) continue;   // "good" sample only, for cleaner resolution/mis-ID curves

            const float pt = p_truth[tr];
            if (!(pt > 0 && pt < PMAX)) continue;
            const int nh = npoints[tr];

            // -- momentum resolution --
            const float pf = p[tr];
            if (std::isfinite(pf) && pf > 0 && pf < PMAX) {
                const double dp  = (pf - pt)/pt;
                const double dip = 1.0/pf - 1.0/pt;
                pFit.push_back(pf); pTru.push_back(pt);
                dpp_all.push_back(dp); dinvp_all.push_back(dip);
                h2P->Fill(std::log10(pt), std::log10(pf));

                for (int b = 0; b < NBINS_P; ++b)
                    if (pt >= pEdges[b] && pt < pEdges[b+1])
                        { dpp_bin[b].push_back(dp); dinvp_bin[b].push_back(dip); }

                if (nh >= NHITS_MIN && nh <= NHITS_MAX) {
                    int hb = nh - NHITS_MIN;
                    dpp_hbin[hb].push_back(dp); dinvp_hbin[hb].push_back(dip);
                }
            }

            // -- charge resolution / mis-ID --
            if (hasQTruth && std::abs(charge_truth[tr]) > 0.5) {
                const double qt = charge_truth[tr];
                const double qf = charge[tr];
                h2Charge->Fill(qt, qf);
                const bool misID = (qf*qt < 0);
                ++nQTot; if (misID) ++nQMis;

                for (int b = 0; b < NBINS_P; ++b)
                    if (pt >= pEdges[b] && pt < pEdges[b+1])
                        { ++nTot_p[b]; if (misID) ++nMis_p[b]; }
                if (nh >= NHITS_MIN && nh <= NHITS_MAX) {
                    int hb = nh - NHITS_MIN;
                    ++nTot_h[hb]; if (misID) ++nMis_h[hb];
                }
            }
        }
    }

    auto vmean = [](const std::vector<double>& v)->double{
        if (v.empty()) return 0; double s=0; for(double x:v) s+=x; return s/v.size(); };
    auto vrms = [&](const std::vector<double>& v)->double{
        if (v.size()<2) return 0; double mu=vmean(v),s=0;
        for(double x:v) s+=(x-mu)*(x-mu); return std::sqrt(s/v.size()); };
    auto binomErr = [](int k,int n)->double{
        if (n<=0) return 0; double p=double(k)/n; return std::sqrt(p*(1-p)/n); };

    printf("\n========================================\n");
    printf("  MDT Performance Summary: %s\n", fname);
    printf("========================================\n");
    printf("  Momentum: N=%zu  mean(dp/p)=%+.3f RMS=%.3f  mean(d(1/p))=%+.3e RMS=%.3e\n",
           dpp_all.size(), vmean(dpp_all), vrms(dpp_all), vmean(dinvp_all), vrms(dinvp_all));
    if (hasQTruth)
        printf("  Charge mis-ID rate: %d/%d = %.2f%% +/- %.2f%%\n",
               nQMis, nQTot, nQTot? 100.0*nQMis/nQTot:0.0, nQTot? 100.0*binomErr(nQMis,nQTot):0.0);
    printf("========================================\n\n");

    // ── build resolution-vs-p graphs ────────────────────────────────────────
    std::vector<double> xc,ex,yrInvP,eyrInvP,yrPP;
    std::vector<double> xcMis,eyMis_lo,eyMis_hi,yMis;
    for (int b=0;b<NBINS_P;++b) {
        if (dinvp_bin[b].size() < 2) continue;
        double xmid = std::sqrt(pEdges[b]*pEdges[b+1]);
        double hw   = (pEdges[b+1]-pEdges[b])/2.0;
        xc.push_back(xmid); ex.push_back(hw);
        double r = vrms(dinvp_bin[b]);
        yrInvP.push_back(r); eyrInvP.push_back(r/std::sqrt(dinvp_bin[b].size()));
        yrPP.push_back(vrms(dpp_bin[b]));
    }
    if (hasQTruth) {
        for (int b=0;b<NBINS_P;++b) {
            if (nTot_p[b] < 5) continue;
            double xmid = std::sqrt(pEdges[b]*pEdges[b+1]);
            xcMis.push_back(xmid);
            double rate = double(nMis_p[b])/nTot_p[b];
            yMis.push_back(rate);
            double e = binomErr(nMis_p[b], nTot_p[b]);
            eyMis_lo.push_back(e); eyMis_hi.push_back(e);
        }
    }
    int nGp = (int)xc.size();
    TGraphErrors* grInvP = new TGraphErrors(nGp, xc.data(), yrInvP.data(), ex.data(), eyrInvP.data());
    TGraph*       grPP   = new TGraph(nGp, xc.data(), yrPP.data());
    int nGm = (int)xcMis.size();
    TGraphErrors* grMisP = hasQTruth ? new TGraphErrors(nGm, xcMis.data(), yMis.data(), nullptr, eyMis_lo.data()) : nullptr;

    // ── resolution / mis-ID vs nhits ────────────────────────────────────────
    std::vector<double> xh,yrInvP_h,eyrInvP_h,yrPP_h,xhMis,yMis_h,eyMis_h;
    for (int b=0;b<NBINS_H;++b) {
        if (dinvp_hbin[b].size() < 2) continue;
        xh.push_back(NHITS_MIN+b);
        double r = vrms(dinvp_hbin[b]);
        yrInvP_h.push_back(r); eyrInvP_h.push_back(r/std::sqrt(dinvp_hbin[b].size()));
        yrPP_h.push_back(vrms(dpp_hbin[b]));
    }
    if (hasQTruth) {
        for (int b=0;b<NBINS_H;++b) {
            if (nTot_h[b] < 5) continue;
            xhMis.push_back(NHITS_MIN+b);
            yMis_h.push_back(double(nMis_h[b])/nTot_h[b]);
            eyMis_h.push_back(binomErr(nMis_h[b], nTot_h[b]));
        }
    }
    int nHp = (int)xh.size();
    TGraphErrors* grInvP_h = new TGraphErrors(nHp, xh.data(), yrInvP_h.data(), nullptr, eyrInvP_h.data());
    TGraph*       grPP_h   = new TGraph(nHp, xh.data(), yrPP_h.data());
    int nHm = (int)xhMis.size();
    TGraphErrors* grMisP_h = hasQTruth ? new TGraphErrors(nHm, xhMis.data(), yMis_h.data(), nullptr, eyMis_h.data()) : nullptr;

    auto styleG = [](TGraph* g, int col, int mk){
        g->SetLineColor(col); g->SetLineWidth(2);
        g->SetMarkerColor(col); g->SetMarkerStyle(mk); g->SetMarkerSize(1.1); };
    styleG(grInvP,  kBlue+1, 20);
    styleG(grPP,    kBlue+1, 20);
    styleG(grInvP_h,kBlue+1, 20);
    styleG(grPP_h,  kBlue+1, 20);
    if (grMisP)   styleG(grMisP,   kRed+1, 21);
    if (grMisP_h) styleG(grMisP_h, kRed+1, 21);

    // =========================================================================
    TCanvas* c = new TCanvas("c","MDT Performance",1400,800);
    c->Print(Form("%s[", outpdf));

    // ── Page 1: truth vs reco momentum (2D) + charge confusion matrix ───────
    c->Clear(); c->Divide(2,1);
    c->cd(1);
    h2P->GetXaxis()->SetTitle("log_{10}(p_{truth} / GeV)");
    h2P->GetYaxis()->SetTitle("log_{10}(p_{fit} / GeV)");
    h2P->SetTitle(Form("p_{fit} vs p_{truth} (#chi^{2}/NDF<%.0f)",chi2Cut));
    h2P->Draw("COLZ");
    TLine* diag = new TLine(std::log10(5.0),std::log10(5.0),std::log10(5000.0),std::log10(5000.0));
    diag->SetLineColor(kRed); diag->SetLineWidth(2); diag->Draw();

    c->cd(2);
    if (hasQTruth) {
        h2Charge->SetStats(0);
        h2Charge->Draw("COLZ TEXT");
    }
    c->Print(Form("%s", outpdf));


    // ── Page 2: RMS(Delta(1/p)) and RMS(Deltap/p) vs p_truth ────────────────
    c->Clear(); c->Divide(1,2);
    c->cd(1); gPad->SetLogx();
    grInvP->SetTitle("RMS(#Delta(1/p)) vs p_{truth};p_{truth} [GeV/c];RMS(1/p_{fit}-1/p_{truth}) [GeV^{-1}]");
    grInvP->Draw("APL");
    c->cd(2); gPad->SetLogx();
    grPP->SetTitle("RMS(#Deltap/p) vs p_{truth};p_{truth} [GeV/c];RMS((p_{fit}-p_{truth})/p_{truth})");
    grPP->Draw("APL");
    c->Print(Form("%s", outpdf));


    // ── Page 3: RMS resolution vs N_hits ────────────────────────────────────
    c->Clear(); c->Divide(1,2);
    c->cd(1);
    grInvP_h->SetTitle("RMS(#Delta(1/p)) vs N_{hits};N_{hits} in fit;RMS(1/p_{fit}-1/p_{truth}) [GeV^{-1}]");
    grInvP_h->Draw("APL");
    c->cd(2);
    grPP_h->SetTitle("RMS(#Deltap/p) vs N_{hits};N_{hits} in fit;RMS((p_{fit}-p_{truth})/p_{truth})");
    grPP_h->Draw("APL");
    c->Print(Form("%s", outpdf));


    // ── Page 4: charge mis-ID rate vs p_truth and vs N_hits ─────────────────
    if (hasQTruth) {
        c->Clear(); c->Divide(2,1);
        c->cd(1); gPad->SetLogx();
        grMisP->SetTitle("Charge mis-ID rate vs p_{truth};p_{truth} [GeV/c];mis-ID rate");
        grMisP->GetYaxis()->SetRangeUser(0.0, std::max(0.05, 1.2*TMath::MaxElement(grMisP->GetN(),grMisP->GetY())));
        grMisP->Draw("APL");
        c->cd(2);
        grMisP_h->SetTitle("Charge mis-ID rate vs N_{hits};N_{hits} in fit;mis-ID rate");
        grMisP_h->GetYaxis()->SetRangeUser(0.0, std::max(0.05, 1.2*TMath::MaxElement(grMisP_h->GetN(),grMisP_h->GetY())));
        grMisP_h->Draw("APL");
        c->Print(Form("%s", outpdf));

    }

    // ── Page 5: Deltap/p and Delta(1/p) 1D distributions ────────────────────
    c->Clear(); c->Divide(2,1);
    TH1D* hDpp  = new TH1D("hDpp","#Deltap/p;(p_{fit}-p_{truth})/p_{truth};Tracks",80,-2,2);
    TH1D* hDinv = new TH1D("hDinv","#Delta(1/p);1/p_{fit}-1/p_{truth} [GeV^{-1}];Tracks",80,-0.1,0.1);
    for (double v: dpp_all)   hDpp->Fill(v);
    for (double v: dinvp_all) hDinv->Fill(v);
    hDpp->SetLineColor(kBlue+1);  hDpp->SetLineWidth(2);
    hDinv->SetLineColor(kBlue+1); hDinv->SetLineWidth(2);
    c->cd(1); hDpp->Draw("HIST");
    c->cd(2); hDinv->Draw("HIST");
    c->Print(Form("%s", outpdf));


    c->Print(Form("%s]", outpdf));
    printf("Saved 5-page PDF: %s\n\n", outpdf);
}
