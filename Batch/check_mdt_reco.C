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

    const bool hasTruth  = (t->GetBranch("p_truth") != nullptr);
    const bool hasQTruth = (t->GetBranch("charge_truth") != nullptr);
    const bool hasNPoints = (t->GetBranch("npoints") != nullptr);
    if (!hasTruth) printf("[WARNING] p_truth branch not found — truth resolution plots disabled.\n");
    if (!hasQTruth) printf("[WARNING] charge_truth branch not found — charge correctness counters disabled.\n");
    if (!hasNPoints) printf("[WARNING] npoints branch not found — fit-failure hit-count breakdown disabled.\n");

    const Long64_t nEntries = t->GetEntries();

    // --- branch variables ---
    const int MAXTRK = 10;
    int   ntracks = 0;
    int   fit_ok [MAXTRK] = {};
    int   nDoF   [MAXTRK] = {};
    int   npoints[MAXTRK] = {};
    float p      [MAXTRK] = {};
    float p_analytic[MAXTRK] = {};
    float p_truth   [MAXTRK] = {};
    float chi2   [MAXTRK] = {};
    float charge [MAXTRK] = {};
    int   charge_mode[MAXTRK] = {};
    float charge_truth[MAXTRK] = {};
    float fpErr  [MAXTRK] = {};

    t->SetBranchAddress("ntracks",     &ntracks);
    t->SetBranchAddress("fit_ok",       fit_ok);
    t->SetBranchAddress("nDoF",         nDoF);
    if (hasNPoints) t->SetBranchAddress("npoints", npoints);
    t->SetBranchAddress("p",            p);
    t->SetBranchAddress("p_analytic",   p_analytic);
    t->SetBranchAddress("chi2",         chi2);
    t->SetBranchAddress("charge",       charge);
    if (t->GetBranch("charge_mode")) t->SetBranchAddress("charge_mode", charge_mode);
    t->SetBranchAddress("fpErr",        fpErr);
    if (hasTruth) t->SetBranchAddress("p_truth", p_truth);
    if (hasQTruth) t->SetBranchAddress("charge_truth", charge_truth);

    // --- counters ---
    int nEventsTotal     = 0;
    int nEventsNoTrack   = 0;
    int nEventsAttempted = 0;
    int nFitFailed       = 0;
    int nSecondTrack     = 0;
    int nNDF0            = 0;

    // fit-failure diagnostic (events with tracks but no fit_ok)
    int nFailMaxHitsLE5  = 0;
    int nFailMaxHits6to8 = 0;
    int nFailMaxHitsGE9  = 0;

    const int NBINS_P = 8;
    double pEdges[NBINS_P+1];
    for (int i = 0; i <= NBINS_P; ++i)
        pEdges[i] = std::pow(10.0, std::log10(5.0) + i*(std::log10(5000.0)-std::log10(5.0))/NBINS_P);

    // charge diagnostics (track-level, fit_ok + NDF>0)
    int nChargeTracks          = 0;
    int nChargeAssignedTracks  = 0;
    int nChargeAmbiguousTracks = 0;
    int nChargeTruthTracks     = 0;
    int nChargeCorrectTracks   = 0;
    int nChargeWrongTracks     = 0;
    int nChargeTot_ndf[8]    = {};
    int nChargeAss_ndf[8]    = {};
    int nChargeWrong_ndf[8]  = {};
    int nChargeTot_mode[6]   = {};
    int nChargeAss_mode[6]   = {};
    int nChargeWrong_mode[6] = {};
    int nChargeTot_pndf[NBINS_P][8] = {};
    int nChargeAss_pndf[NBINS_P][8] = {};
    int nChargeWrong_pndf[NBINS_P][8] = {};

    // charge diagnostics (event-level, leading track tr=0 if fit_ok + NDF>0)
    int nLeadChargeEvents          = 0;
    int nLeadChargeAssignedEvents  = 0;
    int nLeadChargeAmbiguousEvents = 0;
    int nLeadChargeTruthEvents     = 0;
    int nLeadChargeCorrectEvents   = 0;
    int nLeadChargeWrongEvents     = 0;
    int nLeadChargeTot_ndf[8]      = {};
    int nLeadChargeAss_ndf[8]      = {};
    int nLeadChargeWrong_ndf[8]    = {};
    int nLeadChargeTot_mode[6]     = {};
    int nLeadChargeAss_mode[6]     = {};
    int nLeadChargeWrong_mode[6]    = {};
    int nLeadChargeTot_pndf[NBINS_P][8] = {};
    int nLeadChargeAss_pndf[NBINS_P][8] = {};
    int nLeadChargeWrong_pndf[NBINS_P][8] = {};
    int nChargeTot_p[NBINS_P]      = {};
    int nChargeAss_p[NBINS_P]      = {};
    int nChargeWrong_p[NBINS_P]    = {};
    int nLeadChargeTot_p[NBINS_P]   = {};
    int nLeadChargeAss_p[NBINS_P]   = {};
    int nLeadChargeWrong_p[NBINS_P] = {};

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
        int maxHitsThisEvent = -1;
        for (int tr = 0; tr < std::min(ntracks, MAXTRK); ++tr) {
            if (hasNPoints) maxHitsThisEvent = std::max(maxHitsThisEvent, npoints[tr]);
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

            // Charge diagnostics on all successful tracks with valid NDF.
            const bool assigned = (std::abs(charge[tr]) > 0.5f);
            ++nChargeTracks;
            if (assigned) ++nChargeAssignedTracks;
            else          ++nChargeAmbiguousTracks;

            if (hasQTruth && std::abs(charge_truth[tr]) > 0.5f) {
                ++nChargeTruthTracks;
                int mode = (charge_mode[tr] >= 0 && charge_mode[tr] <= 5) ? charge_mode[tr] : 0;
                ++nChargeTot_mode[mode];
                if (assigned) ++nChargeAss_mode[mode];
                if (ndf >= 1 && ndf <= 7) {
                    ++nChargeTot_ndf[ndf];
                    if (assigned) ++nChargeAss_ndf[ndf];
                }
                if (std::isfinite(pt) && pt > 0 && pt < PMAX && ndf >= 1 && ndf <= 7) {
                    for (int b = 0; b < NBINS_P; ++b) {
                        if (pt >= pEdges[b] && pt < pEdges[b+1]) {
                            ++nChargeTot_pndf[b][ndf];
                            if (assigned) ++nChargeAss_pndf[b][ndf];
                            break;
                        }
                    }
                }
                if (std::isfinite(pt) && pt > 0 && pt < PMAX) {
                    for (int b = 0; b < NBINS_P; ++b) {
                        if (pt >= pEdges[b] && pt < pEdges[b+1]) {
                            ++nChargeTot_p[b];
                            if (assigned) ++nChargeAss_p[b];
                            break;
                        }
                    }
                }
                if (!assigned) {
                    // Count as ambiguous, not as wrong-sign.
                } else if (charge[tr] * charge_truth[tr] > 0.f) {
                    ++nChargeCorrectTracks;
                } else {
                    ++nChargeWrongTracks;
                    ++nChargeWrong_mode[mode];
                    if (ndf >= 1 && ndf <= 7) ++nChargeWrong_ndf[ndf];
                        if (std::isfinite(pt) && pt > 0 && pt < PMAX && ndf >= 1 && ndf <= 7) {
                            for (int b = 0; b < NBINS_P; ++b) {
                                if (pt >= pEdges[b] && pt < pEdges[b+1]) {
                                    ++nChargeWrong_pndf[b][ndf];
                                    break;
                                }
                            }
                        }
                    if (std::isfinite(pt) && pt > 0 && pt < PMAX) {
                        for (int b = 0; b < NBINS_P; ++b) {
                            if (pt >= pEdges[b] && pt < pEdges[b+1]) {
                                ++nChargeWrong_p[b];
                                break;
                            }
                        }
                    }
                }
            }

            // Event-level leading-track charge diagnostics.
            if (tr == 0) {
                ++nLeadChargeEvents;
                if (assigned) ++nLeadChargeAssignedEvents;
                else          ++nLeadChargeAmbiguousEvents;

                if (hasQTruth && std::abs(charge_truth[tr]) > 0.5f) {
                    ++nLeadChargeTruthEvents;
                    int mode = (charge_mode[tr] >= 0 && charge_mode[tr] <= 5) ? charge_mode[tr] : 0;
                    ++nLeadChargeTot_mode[mode];
                    if (assigned) ++nLeadChargeAss_mode[mode];
                    if (ndf >= 1 && ndf <= 7) {
                        ++nLeadChargeTot_ndf[ndf];
                        if (assigned) ++nLeadChargeAss_ndf[ndf];
                    }
                    if (std::isfinite(pt) && pt > 0 && pt < PMAX && ndf >= 1 && ndf <= 7) {
                        for (int b = 0; b < NBINS_P; ++b) {
                            if (pt >= pEdges[b] && pt < pEdges[b+1]) {
                                ++nLeadChargeTot_pndf[b][ndf];
                                if (assigned) ++nLeadChargeAss_pndf[b][ndf];
                                break;
                            }
                        }
                    }
                    if (std::isfinite(pt) && pt > 0 && pt < PMAX) {
                        for (int b = 0; b < NBINS_P; ++b) {
                            if (pt >= pEdges[b] && pt < pEdges[b+1]) {
                                ++nLeadChargeTot_p[b];
                                if (assigned) ++nLeadChargeAss_p[b];
                                break;
                            }
                        }
                    }
                    if (!assigned) {
                        // ambiguous
                    } else if (charge[tr] * charge_truth[tr] > 0.f) {
                        ++nLeadChargeCorrectEvents;
                    } else {
                        ++nLeadChargeWrongEvents;
                        ++nLeadChargeWrong_mode[mode];
                        if (ndf >= 1 && ndf <= 7) ++nLeadChargeWrong_ndf[ndf];
                        if (std::isfinite(pt) && pt > 0 && pt < PMAX && ndf >= 1 && ndf <= 7) {
                            for (int b = 0; b < NBINS_P; ++b) {
                                if (pt >= pEdges[b] && pt < pEdges[b+1]) {
                                    ++nLeadChargeWrong_pndf[b][ndf];
                                    break;
                                }
                            }
                        }
                        if (std::isfinite(pt) && pt > 0 && pt < PMAX) {
                            for (int b = 0; b < NBINS_P; ++b) {
                                if (pt >= pEdges[b] && pt < pEdges[b+1]) {
                                    ++nLeadChargeWrong_p[b];
                                    break;
                                }
                            }
                        }
                    }
                }
            }

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
        if (!anyOk) {
            ++nFitFailed;
            if (hasNPoints) {
                if (maxHitsThisEvent <= 5)      ++nFailMaxHitsLE5;
                else if (maxHitsThisEvent <= 8) ++nFailMaxHits6to8;
                else                            ++nFailMaxHitsGE9;
            }
        }
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
    if (hasNPoints && nFitFailed > 0) {
        printf("    - failure with max npoints <=5 : %d\n", nFailMaxHitsLE5);
        printf("    - failure with max npoints 6-8 : %d\n", nFailMaxHits6to8);
        printf("    - failure with max npoints >=9 : %d\n", nFailMaxHitsGE9);
    }
    printf("\n");
    printf("  NDF distribution (fit_ok tracks):\n");
    for (auto& kv : ndfHist)
        printf("    NDF = %3d : %d tracks\n", kv.first, kv.second);
    printf("  fit_ok tracks with NDF=0      : %d  (DAF garbage)\n", nNDF0);
    printf("  Mean chi2/NDF (NDF>0)         : %.2f\n", mu_chi2);
        printf("\n");
        printf("  Charge (track-level, fit_ok and NDF>0):\n");
        printf("    assigned (|q|>0.5)          : %d/%d (%.1f%%)\n",
            nChargeAssignedTracks, nChargeTracks,
            nChargeTracks ? 100.0*nChargeAssignedTracks/nChargeTracks : 0.0);
        printf("    ambiguous (q=0)             : %d/%d (%.1f%%)\n",
            nChargeAmbiguousTracks, nChargeTracks,
            nChargeTracks ? 100.0*nChargeAmbiguousTracks/nChargeTracks : 0.0);
        if (hasQTruth) {
         printf("    correct sign (assigned only): %d/%d (%.1f%%)\n",
             nChargeCorrectTracks, nChargeAssignedTracks,
             nChargeAssignedTracks ? 100.0*nChargeCorrectTracks/nChargeAssignedTracks : 0.0);
         printf("    wrong sign   (assigned only): %d/%d (%.1f%%)\n",
             nChargeWrongTracks, nChargeAssignedTracks,
             nChargeAssignedTracks ? 100.0*nChargeWrongTracks/nChargeAssignedTracks : 0.0);
         printf("    truth-tagged tracks          : %d\n", nChargeTruthTracks);
        }
        printf("\n");
        printf("  Charge (event-level, leading track tr=0):\n");
        printf("    assigned events             : %d/%d (%.1f%%)\n",
            nLeadChargeAssignedEvents, nLeadChargeEvents,
            nLeadChargeEvents ? 100.0*nLeadChargeAssignedEvents/nLeadChargeEvents : 0.0);
        printf("    ambiguous events            : %d/%d (%.1f%%)\n",
            nLeadChargeAmbiguousEvents, nLeadChargeEvents,
            nLeadChargeEvents ? 100.0*nLeadChargeAmbiguousEvents/nLeadChargeEvents : 0.0);
        if (hasQTruth) {
         printf("    correct sign (assigned only): %d/%d (%.1f%%)\n",
             nLeadChargeCorrectEvents, nLeadChargeAssignedEvents,
             nLeadChargeAssignedEvents ? 100.0*nLeadChargeCorrectEvents/nLeadChargeAssignedEvents : 0.0);
         printf("    wrong sign   (assigned only): %d/%d (%.1f%%)\n",
             nLeadChargeWrongEvents, nLeadChargeAssignedEvents,
             nLeadChargeAssignedEvents ? 100.0*nLeadChargeWrongEvents/nLeadChargeAssignedEvents : 0.0);
         printf("    truth-tagged leading tracks : %d\n", nLeadChargeTruthEvents);
        }
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
    if (hasQTruth) {
        auto rate = [](int k, int n) -> double { return (n > 0) ? 100.0 * double(k) / double(n) : 0.0; };
        printf("\n  Charge vs truth momentum (track-level, fit_ok and NDF>0):\n");
        printf("    p_truth bin [GeV/c]         assigned  wrong-sign  wrong/assigned\n");
        for (int b = 0; b < NBINS_P; ++b) {
            printf("    [%7.1f, %7.1f)   %8d  %9d  %6.1f%%\n",
                   pEdges[b], pEdges[b+1], nChargeAss_p[b], nChargeWrong_p[b], rate(nChargeWrong_p[b], nChargeAss_p[b]));
        }
        printf("\n  Charge vs truth momentum (event-level, leading track tr=0):\n");
        printf("    p_truth bin [GeV/c]         assigned  wrong-sign  wrong/assigned\n");
        for (int b = 0; b < NBINS_P; ++b) {
            printf("    [%7.1f, %7.1f)   %8d  %9d  %6.1f%%\n",
                   pEdges[b], pEdges[b+1], nLeadChargeAss_p[b], nLeadChargeWrong_p[b], rate(nLeadChargeWrong_p[b], nLeadChargeAss_p[b]));
        }

        printf("\n  Charge vs NDF (track-level, fit_ok and truth-tagged):\n");
        printf("    NDF   assigned  wrong-sign  wrong/assigned\n");
        for (int ndf = 1; ndf <= 7; ++ndf) {
            printf("    %3d   %8d  %9d  %6.1f%%\n",
                   ndf, nChargeAss_ndf[ndf], nChargeWrong_ndf[ndf], rate(nChargeWrong_ndf[ndf], nChargeAss_ndf[ndf]));
        }

        printf("\n  Charge vs NDF (event-level, leading track tr=0):\n");
        printf("    NDF   assigned  wrong-sign  wrong/assigned\n");
        for (int ndf = 1; ndf <= 7; ++ndf) {
            printf("    %3d   %8d  %9d  %6.1f%%\n",
                   ndf, nLeadChargeAss_ndf[ndf], nLeadChargeWrong_ndf[ndf], rate(nLeadChargeWrong_ndf[ndf], nLeadChargeAss_ndf[ndf]));
        }

        printf("\n  Charge vs charge_mode (track-level, truth-tagged):\n");
        printf("    [mode key: 0=both-hyp ambiguous, 1=clear chi2 winner, 2=only mu- ok,\n");
        printf("               3=only mu+ ok (charge flipped), 4=slope tiebreaker, 5=rescue]\n");
        printf("    mode  assigned  wrong-sign  wrong/assigned\n");
        for (int mode = 0; mode <= 5; ++mode) {
            printf("    %4d  %8d  %9d  %6.1f%%\n",
                   mode, nChargeAss_mode[mode], nChargeWrong_mode[mode], rate(nChargeWrong_mode[mode], nChargeAss_mode[mode]));
        }
        printf("    -- what-if flip sign within mode --\n");
        printf("    mode  correct_now/assigned  correct_if_flip/assigned\n");
        for (int mode = 0; mode <= 5; ++mode) {
            const int ass = nChargeAss_mode[mode];
            const int wrong = nChargeWrong_mode[mode];
            const int correctNow = ass - wrong;
            const int correctIfFlip = wrong;
            printf("    %4d  %8d/%d (%.1f%%)     %8d/%d (%.1f%%)\n",
                   mode,
                   correctNow, ass, rate(correctNow, ass),
                   correctIfFlip, ass, rate(correctIfFlip, ass));
        }

        printf("\n  Charge vs charge_mode (event-level, leading track tr=0):\n");
        printf("    mode  assigned  wrong-sign  wrong/assigned\n");
        for (int mode = 0; mode <= 5; ++mode) {
            printf("    %4d  %8d  %9d  %6.1f%%\n",
                   mode, nLeadChargeAss_mode[mode], nLeadChargeWrong_mode[mode], rate(nLeadChargeWrong_mode[mode], nLeadChargeAss_mode[mode]));
        }
        printf("    -- what-if flip sign within mode --\n");
        printf("    mode  correct_now/assigned  correct_if_flip/assigned\n");
        for (int mode = 0; mode <= 5; ++mode) {
            const int ass = nLeadChargeAss_mode[mode];
            const int wrong = nLeadChargeWrong_mode[mode];
            const int correctNow = ass - wrong;
            const int correctIfFlip = wrong;
            printf("    %4d  %8d/%d (%.1f%%)     %8d/%d (%.1f%%)\n",
                   mode,
                   correctNow, ass, rate(correctNow, ass),
                   correctIfFlip, ass, rate(correctIfFlip, ass));
        }

        printf("\n  Charge vs p_truth x NDF (track-level):\n");
        printf("    p_truth bin [GeV/c]      NDF  assigned  wrong-sign  wrong/assigned\n");
        for (int b = 0; b < NBINS_P; ++b) {
            for (int ndf = 1; ndf <= 7; ++ndf) {
                if (nChargeAss_pndf[b][ndf] == 0 && nChargeWrong_pndf[b][ndf] == 0) continue;
                printf("    [%7.1f, %7.1f)     %3d   %8d  %9d  %6.1f%%\n",
                       pEdges[b], pEdges[b+1], ndf,
                       nChargeAss_pndf[b][ndf], nChargeWrong_pndf[b][ndf],
                       rate(nChargeWrong_pndf[b][ndf], nChargeAss_pndf[b][ndf]));
            }
        }

        printf("\n  Charge vs p_truth x NDF (event-level, leading track tr=0):\n");
        printf("    p_truth bin [GeV/c]      NDF  assigned  wrong-sign  wrong/assigned\n");
        for (int b = 0; b < NBINS_P; ++b) {
            for (int ndf = 1; ndf <= 7; ++ndf) {
                if (nLeadChargeAss_pndf[b][ndf] == 0 && nLeadChargeWrong_pndf[b][ndf] == 0) continue;
                printf("    [%7.1f, %7.1f)     %3d   %8d  %9d  %6.1f%%\n",
                       pEdges[b], pEdges[b+1], ndf,
                       nLeadChargeAss_pndf[b][ndf], nLeadChargeWrong_pndf[b][ndf],
                       rate(nLeadChargeWrong_pndf[b][ndf], nLeadChargeAss_pndf[b][ndf]));
            }
        }
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
