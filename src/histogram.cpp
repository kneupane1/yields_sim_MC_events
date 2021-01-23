/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "histogram.hpp"

Histogram::Histogram(const std::string &output_file) {
        RootOutputFile = std::make_shared<TFile>(output_file.c_str(), "RECREATE");
        def = std::make_shared<TCanvas>("def");
        if (getenv("BEAM_E") != NULL) {
                if (atof(getenv("BEAM_E")) < 3) {
                        q2_max = 1.0;
                        w_max = 3.5;
                        p_max = 3.0;
                } else if (atof(getenv("BEAM_E")) < 8) {
                        q2_max = 4.0;
                        w_max = 4.0;
                        p_max = 4.0;
                } else if (atof(getenv("BEAM_E")) < 9) {
                        q2_max = 7.0;
                        w_max = 7.0;
                        p_max = 7.0;
                }
        }
// background substraction
        eff_check_mPim_all_W = std::make_shared<TH1D>("MMSQ_missingpim_all_W", "MMSQ_missingpim_all_W", bins,
                                                      -0.4, 0.4);
        eff_check_mPim_all_W_with_mmsq_cuts = std::make_shared<TH1D>("MMSQ_missingpim_all_W_with_MMSQ_Cuts", "MMSQ_missingpim_all_W_with_MMSQ_Cuts", bins,
                                                                     -0.4, 0.4);
        eff_check_exclusive_all_W = std::make_shared<TH1D>("MMSQ_exclusive_all_W", "MMSQ_exclusive_all_W", bins,
                                                           -0.4, 0.4);
        eff_check_exclusive_MMSQ_cuts_all_W = std::make_shared<TH1D>("after_MMSQ_cuts_exclusive_all_W", "after_MMSQ_cuts_exclusive_all_W", bins, -0.4, 0.4);

        // const Int_t ndims_MMSQ = 1;
        // Int_t bins_MMSQ[ndims_MMSQ] = {80};
        // Double_t xmin_MMSQ[ndims_MMSQ] = {-0.4};
        // Double_t xmax_MMSQ[ndims_MMSQ] = {0.4};
        for (size_t e = 0; e < EFF_CONDITIONS_NUM_MMSQ; e++) {
                for (short w = 0; w < 64; w++) {
                        for (size_t t = 0; t < THETA_BINS_NUM; t++) {
                                ostringstream qqq;
                                qqq.str("");

                                // Defining the histogram using all the above ingredients
                                qqq << "MMSQ_"
                                        //  <<(1.0 + 1.0 * q2) << "<=Q2<=" << (1.0 + 1.0 * q2+ 1.0) << " GeV2_"
                                    <<EFF_CONDITIONS_NAME_MMSQ[e].c_str()
                                    <<(1.0 + 0.05 * w)<<"<=W<="
                                    << (1.0 + 0.05 * w + 0.05)<<" GeV"
                                    << THETA_BINS_NAME[t].c_str();

                                // eff_check_mPim[EFF_CONDITIONS_NUM][w_bin] = std::make_shared<TH1D>(qqq.str().c_str(), qqq.str().c_str(), 80, -0.4, 0.4);
                                eff_check_mPim[e][w][t] = std::make_shared<TH1D>(qqq.str().c_str(), qqq.str().c_str(), 80, -0.4, 0.4);
                        }
                }

        }
        const Int_t ndims_4D = 4;

        Int_t bins_4D[ndims_4D] = {100, 12, 10, 10};
        Double_t xmin_4D[ndims_4D] = {-0.4, 0, 0., 0.};
        Double_t xmax_4D[ndims_4D] = { 0.4, p_max, 180, 360};

        for (size_t e = 0; e < EFF_CONDITIONS_NUM; e++) {
                // Defining the histogram using all the above ingredients
                background_4DHist_mpim[e] = std::make_shared<THnSparseD>(
                        Form("eff_checks_4D_%6.20s ", EFF_CONDITIONS_NAME[e].c_str()),
                        Form("eff_checks_4D_%6.20s ", EFF_CONDITIONS_NAME[e].c_str()), ndims_4D, bins_4D, xmin_4D, xmax_4D);
                background_4DHist_mpim[e]->Sumw2();
        }
        // const Int_t ndims = 7;
        //
        // Int_t bins_7D[ndims] = {6, 6, 30, 30, 30, 30, 30};
        // Double_t xmin[ndims] = {1.0, 1.0, 1.1, 0.2, 0, 0, 0};
        // Double_t xmax[ndims] = {4.2, 9.0, 2.1, 1.2, 180, 360, 360};

        const Int_t ndims_5D = 5;

        Int_t bins_5D[ndims_5D] = {12, 12, 10, 6, 8};
        Double_t xmin_5D[ndims_5D] = {(0.938272 + 0.13957), (0.13957 + 0.13957), 0., 0., 0.};
        Double_t xmax_5D[ndims_5D] = { 2.1, 1.2, 180, 360, 360};
        //  Int_t bins_3D[4] = {76, 22, 200, 200};
        //  Double_t xmin_3D[4] = {0.5, 1.0, -0.02, -0.03};
        //  Double_t xmax_3D[4] = {4.3, 12.0, 0.02, 0.03};

        //  threeDHist = new THnSparseD("threeD_hist", "threeD_hist", 4, bins_3D,
        //  xmin_3D,
        //                            xmax_3D);
        // sevenDHist_pim = new THnSparseD("sevenD_hist_pim", "sevenD_hist_pim", ndims,
        //                                 bins_7D, xmin, xmax);
        // sevenD_Hist_thrown_pim =
        //         new THnSparseD("sevenD_hist_thrown_pim", "sevenD_hist_thrown_pim", ndims,
        //                        bins_7D, xmin, xmax);
        // sevenDHist_pip = new THnSparseD("sevenD_hist_pip", "sevenD_hist_pip", ndims,
        //                                 bins_7D, xmin, xmax);
        // sevenD_Hist_thrown_pip =
        //         new THnSparseD("sevenD_hist_thrown_pip", "sevenD_hist_thrown_pip", ndims,
        //                        bins_7D, xmin, xmax);
        for (short q2 = 0; q2 < 11; q2++) {
                for (short w = 0; w < 64; w++) {


                        ostringstream qqq;
                        qqq.str("");

                        // Defining the histogram using all the above ingredients
                        qqq << "h_5dim_"
                            <<(1.0 + 1.0 * q2) << "<=Q2<=" << (1.0 + 1.0 * q2+ 1.0) << " GeV2_"<<(1.0 + 0.05 * w)<<"<=W<="
                            << (1.0 + 0.05 * w + 0.05)<<" GeV";


                        sevenDHist_prot[q2][w] = new THnSparseD(qqq.str().c_str(), qqq.str().c_str(),
                                                                ndims_5D, bins_5D, xmin_5D, xmax_5D);
                        sevenDHist_prot[q2][w] -> Sumw2();

                        sevenD_Hist_thrown_prot[q2][w] = new THnSparseD(qqq.str().c_str(), qqq.str().c_str(),
                                                                        ndims_5D, bins_5D, xmin_5D, xmax_5D);

                        sevenD_Hist_thrown_prot[q2][w] -> Sumw2();
                        sevenDHist_pip[q2][w] = new THnSparseD(qqq.str().c_str(), qqq.str().c_str(),
                                                               ndims_5D, bins_5D, xmin_5D, xmax_5D);

                        sevenDHist_pip[q2][w] -> Sumw2();

                        sevenD_Hist_thrown_pip[q2][w] = new THnSparseD(qqq.str().c_str(), qqq.str().c_str(),
                                                                       ndims_5D, bins_5D, xmin_5D, xmax_5D);

                        sevenD_Hist_thrown_pip[q2][w] -> Sumw2();
                        sevenDHist_pim[q2][w] =new THnSparseD(qqq.str().c_str(), qqq.str().c_str(),
                                                              ndims_5D, bins_5D, xmin_5D, xmax_5D);

                        sevenDHist_pim[q2][w] -> Sumw2();
                        sevenD_Hist_thrown_pim[q2][w] =new THnSparseD(qqq.str().c_str(), qqq.str().c_str(),
                                                                      ndims_5D, bins_5D, xmin_5D, xmax_5D);
                        sevenD_Hist_thrown_pim[q2][w] -> Sumw2();
                }
        }
        // sevenD_Hist_thrown_prot =
        //         new THnSparseD("sevenD_hist_thrown_prot", "sevenD_hist_thrown_prot",
        //                        ndims, bins_7D, xmin, xmax);
        weight_hist = std::make_shared<TH1D>("weight", "weight", 1500, -0.001, 0.005);
        momentum = std::make_shared<TH1D>("mom", "mom", bins, p_min, p_max);
        W_hist = std::make_shared<TH1D>("W", "W", bins, zero, w_max);

        // Theta_vs_mom_x_mu = std::make_shared<TH2D>("theta_vs_mom_x_mu_all_sec",
        // "theta_vs_mom_x_mu_all_sec", bins, zero,
        //                                            p_max, bins, zero, 120);

        W_P2pi_hist = std::make_shared<TH1D>("W_P2pi", "W_P2pi", bins, zero, w_max);

        Q2_hist = std::make_shared<TH1D>("Q2", "Q2", bins, zero, q2_max);
        W_vs_q2 = std::make_shared<TH2D>("W_vs_q2", "W_vs_q2", bins, zero, w_max,
                                         bins, zero, q2_max);

        W_thrown = std::make_shared<TH1D>("W_thrown", "W_thrown", bins, zero, w_max);
        W_vs_Q2_thrown =
                std::make_shared<TH2D>("W_vs_q2_thrown", "W_vs_q2_thrown", bins, zero,
                                       w_max, bins, zero, q2_max);

        theta_prot =
                std::make_shared<TH1D>("theta_prot", "theta_prot", bins, 0.5, 180);
        theta_pip = std::make_shared<TH1D>("theta_pip", "theta_pip", bins, 0.5, 180);
        theta_pim = std::make_shared<TH1D>("theta_pim", "theta_pim", bins, 0.5, 180);

        Phi_gamma = std::make_shared<TH1D>("Phi_gamma", "Phi_gamma", bins, 0, 360);
        Phi_prot = std::make_shared<TH1D>("Phi_prot", "Phi_prot", bins, 0, 360);
        Phi_pip = std::make_shared<TH1D>("Phi_pip", "Phi_pip", bins, 0, 360);
        Phi_pim = std::make_shared<TH1D>("Phi_pim", "Phi_pim", bins, 0, 360);

        alpha_prot =
                std::make_shared<TH1D>("alpha_prot", "alpha_prot", bins, 0.5, 360);
        alpha_pip = std::make_shared<TH1D>("alpha_pip", "alpha_pip", bins, 0.5, 360);
        alpha_pim = std::make_shared<TH1D>("alpha_pim", "alpha_pim", bins, 0.5, 360);
        //
        // theta_prot_mc = std::make_shared<TH1D>("theta_prot_mc", "theta_prot_mc",
        // bins, 0, 180); theta_pip_mc = std::make_shared<TH1D>("theta_pip_mc",
        // "theta_pip_mc", bins, 0, 180); theta_pim_mc =
        // std::make_shared<TH1D>("theta_pim_mc", "theta_pim_mc", bins, 0, 180);
        //
        // Phi_gamma_mc = std::make_shared<TH1D>("Phi_gamma_mc", "Phi_gamma_mc", bins,
        // 0, 360); Phi_prot_mc = std::make_shared<TH1D>("Phi_prot_mc", "Phi_prot_mc",
        // bins, 0, 360); Phi_pip_mc = std::make_shared<TH1D>("Phi_pip_mc",
        // "Phi_pip_mc", bins, 0, 360); Phi_pim_mc =
        // std::make_shared<TH1D>("Phi_pim_mc", "Phi_pim_mc", bins, 0, 360);
        //
        // alpha_prot_mc = std::make_shared<TH1D>("alpha_prot_mc", "alpha_prot_mc",
        // bins, 0, 360); alpha_pip_mc = std::make_shared<TH1D>("alpha_pip_mc",
        // "alpha_pip_mc", bins, 0, 360); alpha_pim_mc =
        // std::make_shared<TH1D>("alpha_pim_mc", "alpha_pim_mc", bins, 0, 360);

        theta_prot_thrown = std::make_shared<TH1D>("theta_prot_thrown",
                                                   "theta_prot_thrown", bins, 0, 180);
        theta_pip_thrown = std::make_shared<TH1D>("theta_pip_thrown",
                                                  "theta_pip_thrown", bins, 0, 180);
        theta_pim_thrown = std::make_shared<TH1D>("theta_pim_thrown",
                                                  "theta_pim_thrown", bins, 0, 180);

        Phi_gamma_thrown = std::make_shared<TH1D>("Phi_gamma_thrown",
                                                  "Phi_gamma_thrown", bins, 0, 360);
        Phi_prot_thrown = std::make_shared<TH1D>("Phi_prot_thrown", "Phi_prot_thrown",
                                                 bins, 0, 360);
        Phi_pip_thrown =
                std::make_shared<TH1D>("Phi_pip_thrown", "Phi_pip_thrown", bins, 0, 360);
        Phi_pim_thrown =
                std::make_shared<TH1D>("Phi_pim_thrown", "Phi_pim_thrown", bins, 0, 360);

        alpha_prot_thrown = std::make_shared<TH1D>("alpha_prot_thrown",
                                                   "alpha_prot_thrown", bins, 0, 360);
        alpha_pip_thrown = std::make_shared<TH1D>("alpha_pip_thrown",
                                                  "alpha_pip_thrown", bins, 0, 360);
        alpha_pim_thrown = std::make_shared<TH1D>("alpha_pim_thrown",
                                                  "alpha_pim_thrown", bins, 0, 360);

        // MM_neutron = std::make_shared<TH1D>("missMass", "missMass", bins,
        // zero, 4.0);
        MM_twoPi = std::make_shared<TH1D>("MMsq_missingpim", "MMsq_missingpim", bins,
                                          -0.2, 0.2);
        MM2_twoPi = std::make_shared<TH1D>("MMSQ_missingpim", "MMSQ_missingPim", bins,
                                           -0.2, 0.2);

        MM2_twoPi_missingPip = std::make_shared<TH1D>(
                "MMSQ_missingpip", "MMSQ_missingPip", bins, -0.2, 0.2);

        MM2_twoPi_missingProt = std::make_shared<TH1D>(
                "MMSQ_missingProt", "MMSQ_missingProt", bins, 0.6, 1.2);
        W_hist_twoPi =
                std::make_shared<TH1D>("W_twoPi", "W_twoPi", bins, zero, w_max);
        Q2_hist_twoPi =
                std::make_shared<TH1D>("Q2_twoPi", "Q2_twoPi", bins, zero, q2_max);
        W_vs_q2_twoPi = std::make_shared<TH2D>("W_vs_q2_twoPi", "W_vs_q2_twoPi", bins,
                                               zero, w_max, bins, zero, q2_max);
        W_vs_q2_twoPi_thrown = std::make_shared<TH2D>("W_vs_q2_twoPi_thrown", "W_vs_q2_twoPi_thrown", bins,
                                                      zero, w_max, bins, zero, q2_max);
        Theta_prot_cm_vs_mom_prot = std::make_shared<TH2D>(
                "Theta_prot_cm_vs_mom_prot", "Theta_prot_cm_vs_mom_prot", bins, zero, 8.5,
                bins, zero, 180);
        Theta_pip_cm_vs_mom_pip = std::make_shared<TH2D>(
                "Theta_pip_cm_vs_mom_pip", "Theta_pip_cm_vs_mom_pip", bins, zero, 8.5,
                bins, zero, 180);
        Theta_pim_cm_vs_mom_pim = std::make_shared<TH2D>(
                "Theta_pim_cm_vs_mom_pim", "Theta_pim_cm_vs_mom_pim", bins, zero, 8.5,
                bins, zero, 180);
        Theta_prot_lab_vs_mom_prot = std::make_shared<TH2D>(
                "Theta_prot_lab_vs_mom_prot", "Theta_prot_lab_vs_mom_prot", bins, zero,
                8.5, bins, zero, 80);
        Theta_pip_lab_vs_mom_pip = std::make_shared<TH2D>(
                "Theta_pip_lab_vs_mom_pip", "Theta_pip_lab_vs_mom_pip", bins, zero, 8.5,
                bins, zero, 130);
        Theta_pim_lab_vs_mom_pim = std::make_shared<TH2D>(
                "Theta_pim_lab_vs_mom_pim", "Theta_pim_lab_vs_mom_pim", bins, zero, 8.5,
                bins, zero, 130);

        Theta_prot_thrown_cm_vs_mom_prot = std::make_shared<TH2D>(
                "Theta_prot_thrown_cm_vs_mom_prot", "Theta_prot_thrown_cm_vs_mom_prot",
                bins, zero, 8.5, bins, zero, 180);
        Theta_pip_thrown_cm_vs_mom_pip = std::make_shared<TH2D>(
                "Theta_pip_thrown_cm_vs_mom_pip", "Theta_pip_thrown_cm_vs_mom_pip", bins,
                zero, 8.5, bins, zero, 180);
        Theta_pim_thrown_cm_vs_mom_pim = std::make_shared<TH2D>(
                "Theta_pim_thrown_cm_vs_mom_pim", "Theta_pim_thrown_cm_vs_mom_pim", bins,
                zero, 8.5, bins, zero, 180);

        Theta_prot_thrown_lab_vs_mom_prot = std::make_shared<TH2D>(
                "Theta_prot_thrown_lab_vs_mom_prot", "Theta_prot_thrown_cm_vs_mom_prot",
                bins, zero, 8.5, bins, zero, 80);
        Theta_pip_thrown_lab_vs_mom_pip = std::make_shared<TH2D>(
                "Theta_pip_thrown_lab_vs_mom_pip", "Theta_pip_thrown_cm_vs_mom_pip", bins,
                zero, 8.5, bins, zero, 130);
        Theta_pim_thrown_lab_vs_mom_pim = std::make_shared<TH2D>(
                "Theta_pim_thrown_lab_vs_mom_pim", "Theta_pim_thrown_cm_vs_mom_pim", bins,
                zero, 8.5, bins, zero, 130);
        // W_hist_singlePi = std::make_shared<TH1D>("W_singlePi", "W_singlePi", bins,
        // zero, w_max); Q2_hist_singlePi = std::make_shared<TH1D>("Q2_singlePi",
        // "Q2_singlePi", bins, zero, q2_max); W_vs_q2_singlePi =
        //     std::make_shared<TH2D>("W_vs_q2_singlePi", "W_vs_q2_singlePi", bins,
        //     zero, w_max, bins, zero, q2_max);

        makeHists_sector();
        makeHists_deltat();
        makeHists_MomVsBeta();
        //  makeHists_electron_cuts();
        makeHists_x_mu();
        makeHistTheta_pim_measured();
        makeHist_efficiency_check();
        makeHist_efficiency_check_1d();


}

Histogram::~Histogram() {
        this->Write();
}

void Histogram::Write() {
        std::cout << GREEN << "Writting" << DEF << std::endl;
        // THnSparse
        // std::cerr << BOLDBLUE << " Hists_3D()" << DEF << std::endl;
        // TDirectory *THnSparse_3D_folder = RootOutputFile->mkdir("THnSparse_3D");
        // THnSparse_3D_folder->cd();
        // writeHists3D();
        //std::cerr << BOLDBLUE << " Hists_7D()" << DEF << std::endl;
        TDirectory *eff_check_mPim =
                RootOutputFile->mkdir("eff_check_mPim");
        eff_check_mPim->cd();
        writeHist_eff_check_mPim();

        // TDirectory *eff_check_MMSQ_cuts_mPim =
        //         RootOutputFile->mkdir("eff_check_MMSQ_cuts_mPim");
        // eff_check_MMSQ_cuts_mPim->cd();
        // writeHist_eff_check_mPim_after_MMSQ_cuts();
        //
        // TDirectory *eff_check_excusive =
        //         RootOutputFile->mkdir("eff_check_excusive");
        // eff_check_excusive->cd();
        // writeHist_eff_check_exclusive();
        //
        // TDirectory *eff_check_after_exclusive_MMSQ_cuts =
        //         RootOutputFile->mkdir("eff_check_after_exclusive_MMSQ_cuts");
        // eff_check_after_exclusive_MMSQ_cuts->cd();
        // writeHist_eff_check_exclusive_MMSQ_cuts();
        // TDirectory *THnSparse_4D_background =
        //         RootOutputFile->mkdir("THnSparse_4D_background");
        // THnSparse_4D_background->cd();
        // writehists4D_background();

//         TDirectory *THnSparse_7D_prot_folder =
//                 RootOutputFile->mkdir("THnSparse_7D_prot");
//         THnSparse_7D_prot_folder->cd();
//         writeHists7D_prot();
//         TDirectory *THnSparse_7D_thrown_prot_folder =
//                 RootOutputFile->mkdir("THnSparse_7D_thrown_prot");
//         THnSparse_7D_thrown_prot_folder->cd();
//         writeHists7D_thrown_prot();
//         TDirectory *THnSparse_7D_pim_folder =
//                 RootOutputFile->mkdir("THnSparse_7D_pim");
//         THnSparse_7D_pim_folder->cd();
//         writeHists7D_pim();
//         TDirectory *THnSparse_7D_thrown_pim_folder =
//                 RootOutputFile->mkdir("THnSparse_7D_thrown_pim");
//         THnSparse_7D_thrown_pim_folder->cd();
//         writeHists7D_thrown_pim();
//         TDirectory *THnSparse_7D_pip_folder =
//                 RootOutputFile->mkdir("THnSparse_7D_pip");
//         THnSparse_7D_pip_folder->cd();
//         writeHists7D_pip();
//         TDirectory *THnSparse_7D_thrown_pip_folder =
//                 RootOutputFile->mkdir("THnSparse_7D_thrown_pip");
//         THnSparse_7D_thrown_pip_folder->cd();
//         writeHists7D_thrown_pip();
// //  Write_EC();
        // std::cerr << BOLDBLUE << "WvsQ2()" << DEF << std::endl;
        // TDirectory *WvsQ2_folder = RootOutputFile->mkdir("W vs Q2");
        // WvsQ2_folder->cd();
        // Write_WvsQ2();
        // std::cerr << BOLDBLUE << "write_hist_x_mu()" << DEF << std::endl;
        // TDirectory *hists_x_mu = RootOutputFile->mkdir("hists_x_mu");
        // hists_x_mu->cd();
        // write_hist_x_mu();

        // std::cerr << BOLDBLUE << "Write_MomVsBeta()" << DEF << std::endl;
        // TDirectory *Write_MomVsBeta_folder = RootOutputFile->mkdir("Mom Vs Beta");
        // Write_MomVsBeta_folder->cd();
        // Write_MomVsBeta();
        //
        // std::cerr << BOLDBLUE << "Write_Electron_cuts()" << DEF << std::endl;
        // TDirectory* Electron_Cuts = RootOutputFile->mkdir("Electron_Cuts");
        // Electron_Cuts->cd();
        // Write_Electron_cuts();

        // std::cerr << BOLDBLUE << "Write_deltat()" << DEF << std::endl;
        // TDirectory* Write_deltat_folder = RootOutputFile->mkdir("Delta_t");
        // Write_deltat_folder->cd();
        // Write_deltat();
        // //
        // // std::cerr << BOLDBLUE << "SF_CUTS()" << DEF << std::endl;
        // // TDirectory* SF_CUTS = RootOutputFile->mkdir("SF_CUTS");
        // // SF_CUTS->cd();
        // // Write_SF();
        // //
        //
        // std::cerr << BOLDBLUE << "write_hist_theta_pim_measured()" << DEF << std::endl;
        // TDirectory* theta_pim_measured = RootOutputFile->mkdir("theta_pim_measured");
        // theta_pim_measured->cd();
        // write_hist_theta_pim_measured();
        //
        // std::cerr << BOLDBLUE << "write_hist_efficiency_check_CD()" << DEF << std::endl;
        // TDirectory* efficiency_check_CD = RootOutputFile->mkdir("efficiency_check_CD");
        // efficiency_check_CD->cd();
        // write_hist_efficiency_check_CD();
        //
        // std::cerr << BOLDBLUE << "write_hist_efficiency_check_FD()" << DEF << std::endl;
        // TDirectory* efficiency_check_FD = RootOutputFile->mkdir("efficiency_check_FD");
        // efficiency_check_FD->cd();
        // write_hist_efficiency_check_FD();
        //
        // std::cerr << BOLDBLUE << "write_hist_efficiency_check_1d()" << DEF << std::endl;
        // TDirectory* efficiency_check_1d_hist = RootOutputFile->mkdir("efficiency_check_1d_hist");
        // efficiency_check_1d_hist->cd();
        // write_hist_efficiency_check_1d();

        std::cerr << BOLDBLUE << "Done Writing!!!" << DEF << std::endl;
}
//


// void Histogram::Fill_histthreeD(const std::shared_ptr<Reaction> &_e) {
//         // fill it
//         const Int_t ndims = 5;
//         Double_t x[ndims];
//         x[0] = _e->W();
//         x[1] = _e->Q2();
//         x[2] = _e->MM2();
//         x[3] = _e->MM2_mpip();
//         x[4] = _e->MM2_mprot();
//         //if (_e->MM_cut()) { // abs(mmsq<0.03)
//         TThread::Lock();
//         threeDHist->Fill(x, _e->weight());
//         TThread::UnLock();
//         threeDHist->GetNbins();
//         // }
// }
// void Histogram::writeHists3D() {
//         threeDHist->Write();
// }





void Histogram::populate_eff_check_mPim(const std::shared_ptr<Reaction> &_e, double min_w, double max_w, double min_theta, double max_theta, double min_phi, double max_phi, short index_w, short index_theta, short index_phi) {
        if (_e->W() > min_w && _e->W() < max_w) {
                if (_e->pim_theta_lab() > min_theta && _e->pim_theta_lab() < max_theta) {
                        //        if (_e->pim_Phi_lab() > min_phi && _e->pim_Phi_lab() < max_phi) {

                        if(_e->TwoPion_missingPim()) {
                                eff_check_mPim[0][index_w][index_theta] /*[index_phi]*/ -> Fill(_e->MM2(), _e->weight());
                                // if (_e->MM_cut()) {
                                //         eff_check_mPim[1][index_w][index_theta] -> Fill(_e->MM2(), _e->weight());
                                // }
                        }
                        // if(_e->TwoPion_exclusive()) {
                        //         eff_check_mPim[2][index_w] -> Fill(_e->MM2(), _e->weight());
                        //         if(abs(_e->MM2_exclusive()) < 0.03 )
                        //                 eff_check_mPim[3][index_w] -> Fill(_e->MM2(), _e->weight());
                        // }
                }
        }
        //}

}

void Histogram::Fill_eff_ckeck_mPim(const std::shared_ptr<Reaction> &_e) {
        short index_w = 0;
        short index_theta = 0;
        short index_phi = 0;
        if(_e->TwoPion_missingPim()) {

                eff_check_mPim_all_W->Fill(_e->MM2(), _e->weight());
                if (_e->MM_cut()) {
                        eff_check_mPim_all_W_with_mmsq_cuts->Fill(_e->MM2(), _e->weight());
                }
        }
        for (float x = 1.; x <= 4.2; x = x + 0.05) {
                float y = x + 0.05;
                //populate_eff_check_mPim(_e, x, y, index_w);
                for (float x_theta = 0.; x_theta <= 150; x_theta = x_theta + 15) {
                        float y_theta = x_theta + 15;
                        for (float x_phi = 0.; x_phi <= 360; x_phi = x_phi + 45) {
                                float y_phi = x_phi + 45;
                                populate_eff_check_mPim(_e, x, y, x_theta, y_theta, x_phi, y_phi, index_w, index_theta, index_phi);
                                index_phi ++;
                        }
                        index_theta ++;
                }

                index_w ++;
        }
}

// void Histogram::populate_eff_check_exclusive(const std::shared_ptr<Reaction> &_e, double min, double max, short index_w) {
//         if(_e->TwoPion_exclusive()) {
//                 eff_check_exclusive_all_W->Fill(_e->MM2(), _e->weight());
//                 if(abs(_e->MM2_exclusive()) < 0.03 ) {
//                         eff_check_exclusive_MMSQ_cuts_all_W->Fill(_e->MM2(), _e->weight());
//                 }
//                 if (_e->W() > min && _e->W() < max) {
//                         eff_check_mPim[2][index_w] -> Fill(_e->MM2(), _e->weight());
//                         // if(abs(_e->MM2_exclusive()) < 0.03 )
//                         //         eff_check_mPim[3][index_w] -> Fill(_e->MM2(), _e->weight());
//                 }
//         }
//
// }
//
// void Histogram::Fill_eff_ckeck_exclusive(const std::shared_ptr<Reaction> &_e) {
//         short index_w = 0;
//
//         for (float x = 1.; x <= 4.2; x = x + 0.05) {
//                 float y = x + 0.05;
//                 populate_eff_check_exclusive(_e, x, y, index_w);
//                 index_w ++;
//         }
// }
void Histogram::writeHist_eff_check_mPim() {
        eff_check_mPim_all_W->Write();
        //for (size_t e = 0; e < EFF_CONDITIONS_NUM; e ++) {

        for (size_t w = 0; w < 64; w ++) {
                for (size_t t = 0; t < THETA_BINS_NUM; t ++) {

                        eff_check_mPim[0][w][t] -> Write();
                        //}
                }
        }
}
// void Histogram::writeHist_eff_check_mPim_after_MMSQ_cuts() {
//         eff_check_mPim_all_W_with_mmsq_cuts->Write();
//
//         //for (size_t e = 0; e < EFF_CONDITIONS_NUM; e ++) {
//
//         for (size_t w = 0; w < 64; w ++) {
//                 eff_check_mPim[1][w] -> Write();
//                 //}
//         }
// }
// void Histogram::writeHist_eff_check_exclusive() {
//         eff_check_exclusive_all_W->Write();
//         //for (size_t e = 0; e < EFF_CONDITIONS_NUM; e ++) {
//
//         for (size_t w = 0; w < 64; w ++) {
//                 eff_check_mPim[2][w] -> Write();
//                 //}
//         }
// }
// void Histogram::writeHist_eff_check_exclusive_MMSQ_cuts() {
//         eff_check_exclusive_MMSQ_cuts_all_W->Write();
//         //for (size_t e = 0; e < EFF_CONDITIONS_NUM; e ++) {
//
//         for (size_t w = 0; w < 64; w ++) {
//                 eff_check_mPim[3][w] -> Write();
//                 //}
//         }
// }

void Histogram::Fill_hists4D_background(const std::shared_ptr<Reaction> &_e) {
        // fill it
        const Int_t ndims = 4;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->MM2();
        x[1] = _e->pim_momentum();
        x[2] = _e->pim_theta_lab();
        x[3] = _e->pim_Phi_lab();

        if (_e->W() <= 4.2 && _e->W() >= 1.0 && _e->Q2() >= 1.0 && _e->Q2() <= 12.0) {
                TThread::Lock();
                background_4DHist_mpim[0] -> Fill(x, _e->weight());
                TThread::UnLock();
                background_4DHist_mpim[0] -> GetNbins();
                if (_e->MM_cut()) {  // abs(mmsq<0.03)
                        TThread::Lock();
                        background_4DHist_mpim[1] -> Fill(x, _e->weight());
                        TThread::UnLock();
                        background_4DHist_mpim[1] -> GetNbins();

                        if(abs(_e->MM2_exclusive()) < 0.03 ) {

                                TThread::Lock();
                                background_4DHist_mpim[2] -> Fill(x, _e->weight());
                                TThread::UnLock();
                                background_4DHist_mpim[2] -> GetNbins();
                        }
                }
        }
}
void Histogram::writehists4D_background() {
        for (size_t e = 0; e < EFF_CONDITIONS_NUM; e ++) {

                background_4DHist_mpim[e] -> Write();
        }
}

void Histogram::Fill_histSevenD_prot(const std::shared_ptr<Reaction> &_e) {
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppip();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->prot_theta();
        x[3] = _e->prot_Phi();
        x[4] = _e->alpha_pippim_pipf();
        if (_e->MM_cut()) {  // abs(mmsq<0.03)
                if (_e->W() <= 4.2 && _e->W() >= 1.0 && _e->Q2() >= 1.0 && _e->Q2() <= 12.0) {
                        TThread::Lock();
                        sevenDHist_prot[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Fill(x, _e->weight());
                        //sevenDHist_prot[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                        TThread::UnLock();
                        sevenDHist_prot[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> GetNbins();
                }
        }
}
void Histogram::writeHists7D_prot() {
        for (size_t q2 = 0; q2 < 11; q2 ++) {
                for (size_t w = 0; w < 64; w ++) {
                        sevenDHist_prot[q2][w] -> Write();
                }
        }
}
void Histogram::Fill_histSevenD_thrown_prot(const std::shared_ptr<MCReaction> &_e) {
        // fill it
        const Int_t ndims = 5;
        Double_t x_thrown[ndims];
        //         x_thrown[0] = _e->W_mc();
        // x_thrown[1] = _e->Q2_mc();
        x_thrown[0] = _e->MCinv_Ppip();
        x_thrown[1] = _e->MCinv_pip_pim();
        x_thrown[2] = _e->MCprot_theta_thrown();
        x_thrown[3] = _e->MCprot_Phi_thrown();
        x_thrown[4] = _e->MCalpha_pippim_pipf_thrown();
        if (_e->W_mc() <= 4.2 && _e->W_mc() >= 1.0 ) {
                if ( _e->Q2_mc() >= 1.0 && _e->Q2_mc() <= 12.0) {
                        TThread::Lock();
                        sevenD_Hist_thrown_prot[int(_e->Q2_mc() - 1.0)][int((_e->W_mc()-1.0)/0.05)] -> Fill(x_thrown, _e->weight());
                        //sevenD_Hist_thrown_prot[int(_e->Q2_mc() - 1.0)][int((_e->W_mc()-1.0)/0.05)] -> Sumw2();
                        TThread::UnLock();
                        sevenD_Hist_thrown_prot[int(_e->Q2_mc() - 1.0)][int((_e->W_mc()-1.0)/0.05)] -> GetNbins();
                }
        }

}
// hN1->Sumw2();
// gDirectory->Append(hN1);
// ret->Add(hN1);

void Histogram::writeHists7D_thrown_prot() {
        for (size_t q2 = 0; q2 < 11; q2 ++) {
                for (size_t w = 0; w < 64; w ++) {
                        sevenD_Hist_thrown_prot[q2][w] -> Write();
                }
        }
}
void Histogram::Fill_histSevenD_pip(const std::shared_ptr<Reaction> &_e) {
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppim();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->pip_theta();
        x[3] = _e->pip_Phi();
        x[4] = _e->alpha_ppim_pipip();
        if (_e->MM_cut()) {  // abs(mmsq<0.03)
                if (_e->W() <= 4.2 && _e->W() >= 1.0 && _e->Q2() >= 1.0 && _e->Q2() <= 12.0) {
                        TThread::Lock();
                        sevenDHist_pip[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Fill(x, _e->weight());
                        //sevenDHist_pip[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                        TThread::UnLock();
                        sevenDHist_pip[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> GetNbins();
                }
        }
}
void Histogram::writeHists7D_pip() {
        for (size_t q2 = 0; q2 < 11; q2 ++) {
                for (size_t w = 0; w < 64; w ++) {
                        sevenDHist_pip[q2][w] -> Write();
                }
        }
}
void Histogram::Fill_histSevenD_thrown_pip(const std::shared_ptr<MCReaction> &_e) {
        // fill it
        const Int_t ndims = 5;
        Double_t x_thrown[ndims];
        //         x_thrown[0] = _e->W_mc();
        // x_thrown[1] = _e->Q2_mc();
        x_thrown[0] = _e->MCinv_Ppim();
        x_thrown[1] = _e->MCinv_pip_pim();
        x_thrown[2] = _e->MCpip_theta_thrown();
        x_thrown[3] = _e->MCpip_Phi_thrown();
        x_thrown[4] = _e->MCalpha_ppim_pipip_thrown();
        if (_e->W_mc() <= 4.2 && _e->W_mc() >= 1.0 ) {
                if ( _e->Q2_mc() >= 1.0 && _e->Q2_mc() <= 12.0) {
                        TThread::Lock();
                        sevenD_Hist_thrown_pip[int(_e->Q2_mc() - 1.0)][int((_e->W_mc()-1.0)/0.05)] -> Fill(x_thrown, _e->weight());
                        //sevenD_Hist_thrown_pip[int(_e->Q2_mc() - 1.0)][int((_e->W_mc()-1.0)/0.05)] -> Sumw2();

                        TThread::UnLock();
                        sevenD_Hist_thrown_pip[int(_e->Q2_mc() - 1.0)][int((_e->W_mc()-1.0)/0.05)] -> GetNbins();
                }
        }

}

void Histogram::writeHists7D_thrown_pip() {
        for (size_t q2 = 0; q2 < 11; q2 ++) {
                for (size_t w = 0; w < 64; w ++) {
                        sevenD_Hist_thrown_pip[q2][w] -> Write();
                }
        }
}
void Histogram::Fill_histSevenD_pim(const std::shared_ptr<Reaction> &_e) {
        // fill it
        const Int_t ndims = 5;
        Double_t x[ndims];
        // x[0] = _e->W();
        //  x[1] = _e->Q2();
        x[0] = _e->inv_Ppip();
        x[1] = _e->inv_pip_pim();
        x[2] = _e->pim_theta();
        x[3] = _e->pim_Phi();
        x[4] =  _e->alpha_ppip_pipim();
        if (_e->MM_cut()) {  // abs(mmsq<0.03)
                if (_e->W() <= 4.2 && _e->W() >= 1.0 && _e->Q2() >= 1.0 && _e->Q2() <= 12.0) {
                        TThread::Lock();
                        sevenDHist_pim[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Fill(x, _e->weight());
                        //sevenDHist_pim[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> Sumw2();
                        TThread::UnLock();
                        sevenDHist_pim[int((_e->Q2() - 1.0)/1.0)][int((_e->W()-1.0)/0.05)] -> GetNbins();
                }
        }
}
void Histogram::writeHists7D_pim() {
        for (size_t q2 = 0; q2 < 11; q2 ++) {
                for (size_t w = 0; w < 64; w ++) {
                        sevenDHist_pim[q2][w] -> Write();
                }
        }
}
void Histogram::Fill_histSevenD_thrown_pim(const std::shared_ptr<MCReaction> &_e) {
        // fill it
        const Int_t ndims = 5;
        Double_t x_thrown[ndims];
        //         x_thrown[0] = _e->W_mc();
        // x_thrown[1] = _e->Q2_mc();
        x_thrown[0] = _e->MCinv_Ppip();
        x_thrown[1] = _e->MCinv_pip_pim();
        x_thrown[2] = _e->MCpim_theta_thrown();
        x_thrown[3] = _e->MCpim_Phi_thrown();
        x_thrown[4] = _e->MCalpha_ppip_pipim_thrown();
        if (_e->W_mc() <= 4.2 && _e->W_mc() >= 1.0 ) {
                if ( _e->Q2_mc() >= 1.0 && _e->Q2_mc() <= 12.0) {
                        TThread::Lock();
                        sevenD_Hist_thrown_pim[int(_e->Q2_mc() - 1.0)][int((_e->W_mc()-1.0)/0.05)] -> Fill(x_thrown, _e->weight());
                        //sevenD_Hist_thrown_pim[int(_e->Q2_mc() - 1.0)][int((_e->W_mc()-1.0)/0.05)] -> Sumw2();
                        TThread::UnLock();
                        sevenD_Hist_thrown_pim[int(_e->Q2_mc() - 1.0)][int((_e->W_mc()-1.0)/0.05)] -> GetNbins();
                }
        }

}

void Histogram::writeHists7D_thrown_pim() {
        for (size_t q2 = 0; q2 < 11; q2 ++) {
                for (size_t w = 0; w < 64; w ++) {
                        sevenD_Hist_thrown_pim[q2][w] -> Write();
                }
        }
}
//
// void Histogram::Fill_histSevenD_prot(const std::shared_ptr<Reaction> &_e) {
//         // fill it
//         const Int_t ndims = 7;
//         Double_t x[ndims];
//         x[0] = _e->W();
//         x[1] = _e->Q2();
//         x[2] = _e->inv_Ppip();
//         x[3] = _e->inv_pip_pim();
//         x[4] = _e->prot_theta();
//         x[5] = _e->prot_Phi();
//         x[6] = _e->alpha_pippim_pipf();
//         if (_e->MM_cut()) { // abs(mmsq<0.03)
//                 TThread::Lock();
//                 sevenDHist_prot->Fill(x, _e->weight());
//                 TThread::UnLock();
//                 sevenDHist_prot->GetNbins();
//         }
// }
// void Histogram::writeHists7D_prot() {
//         sevenDHist_prot->Write();
// }
// void Histogram::Fill_histSevenD_thrown_prot(
//         const std::shared_ptr<MCReaction> &_e) {
//         // fill it
//         const Int_t ndims = 7;
//         Double_t x_thrown[ndims];
//         x_thrown[0] = _e->W_mc();
//         x_thrown[1] = _e->Q2_mc();
//         x_thrown[2] = _e->MCinv_Ppip();
//         x_thrown[3] = _e->MCinv_pip_pim();
//         x_thrown[4] = _e->MCprot_theta_thrown();
//         x_thrown[5] = _e->MCprot_Phi_thrown();
//         x_thrown[6] = _e->MCalpha_pippim_pipf_thrown();
//         //  if (_e->MM_cut()) {  // abs(mmsq<0.03)
//
//         TThread::Lock();
//         sevenD_Hist_thrown_prot->Fill(x_thrown, _e->weight());
//         TThread::UnLock();
//         sevenD_Hist_thrown_prot->GetNbins();
//         //}
// }
// void Histogram::writeHists7D_thrown_prot() {
//         sevenD_Hist_thrown_prot->Write();
// }
//
// void Histogram::Fill_histSevenD_pim(const std::shared_ptr<Reaction> &_e) {
//         // fill it
//         const Int_t ndims = 7;
//         Double_t x[ndims];
//         x[0] = _e->W();
//         x[1] = _e->Q2();
//         x[2] = _e->inv_Ppip();
//         x[3] = _e->inv_pip_pim();
//         x[4] = _e->pim_theta();
//         x[5] = _e->pim_Phi();
//         x[6] = _e->alpha_ppip_pipim();
//         if (_e->MM_cut()) { // abs(mmsq<0.03)
//
//                 TThread::Lock();
//                 sevenDHist_pim->Fill(x, _e->weight());
//                 TThread::UnLock();
//                 sevenDHist_pim->GetNbins();
//         }
// }
// void Histogram::writeHists7D_pim() {
//         sevenDHist_pim->Write();
// }
// void Histogram::Fill_histSevenD_thrown_pim(
//         const std::shared_ptr<MCReaction> &_e) {
//         // fill it
//         const Int_t ndims = 7;
//         Double_t x_thrown[ndims];
//         x_thrown[0] = _e->W_mc();
//         x_thrown[1] = _e->Q2_mc();
//         x_thrown[2] = _e->MCinv_Ppip();
//         x_thrown[3] = _e->MCinv_pip_pim();
//         x_thrown[4] = _e->MCpim_theta_thrown();
//         x_thrown[5] = _e->MCpim_Phi_thrown();
//         x_thrown[6] = _e->MCalpha_ppip_pipim_thrown();
//         //  if (_e->MM_cut()) {  // abs(mmsq<0.03)
//
//         TThread::Lock();
//         sevenD_Hist_thrown_pim->Fill(x_thrown, _e->weight());
//         TThread::UnLock();
//         sevenD_Hist_thrown_pim->GetNbins();
//         //}
// }
// void Histogram::writeHists7D_thrown_pim() {
//         sevenD_Hist_thrown_pim->Write();
// }
//
// void Histogram::Fill_histSevenD_pip(const std::shared_ptr<Reaction> &_e) {
//         // fill it
//         const Int_t ndims = 7;
//         Double_t x[ndims];
//         x[0] = _e->W();
//         x[1] = _e->Q2();
//         x[2] = _e->inv_Ppim();
//         x[3] = _e->inv_pip_pim();
//         x[4] = _e->pip_theta();
//         x[5] = _e->pip_Phi();
//         x[6] = _e->alpha_ppim_pipip();
//         if (_e->MM_cut()) { // abs(mmsq<0.03)
//                 TThread::Lock();
//                 sevenDHist_pip->Fill(x, _e->weight());
//                 TThread::UnLock();
//                 sevenDHist_pip->GetNbins();
//         }
// }
// void Histogram::writeHists7D_pip() {
//         sevenDHist_pip->Write();
// }
// void Histogram::Fill_histSevenD_thrown_pip(
//         const std::shared_ptr<MCReaction> &_e) {
//         // fill it
//         const Int_t ndims = 7;
//         Double_t x_thrown[ndims];
//         x_thrown[0] = _e->W_mc();
//         x_thrown[1] = _e->Q2_mc();
//         x_thrown[2] = _e->MCinv_Ppim();
//         x_thrown[3] = _e->MCinv_pip_pim();
//         x_thrown[4] = _e->MCpip_theta_thrown();
//         x_thrown[5] = _e->MCpip_Phi_thrown();
//         x_thrown[6] = _e->MCalpha_ppim_pipip_thrown();
//         //  if (_e->MM_cut()) {  // abs(mmsq<0.03)
//
//         TThread::Lock();
//         sevenD_Hist_thrown_pip->Fill(x_thrown, _e->weight());
//         TThread::UnLock();
//         sevenD_Hist_thrown_pip->GetNbins();
//         //}
// }
// void Histogram::writeHists7D_thrown_pip() {
//         sevenD_Hist_thrown_pip->Write();
// }

void Histogram::makeHists_x_mu() {
        for (short i = 0; i < NUM_CONDITIONS; i ++) {
                E_x_mu_hist[i] = std::make_shared<TH1D>(
                        Form("E_x_mu_ %4.12s ", NUM_CONDITIONS_NAME[i].c_str()),
                        Form("E_x_mu %4.12s ", NUM_CONDITIONS_NAME[i].c_str()), bins, -2.0,
                        10.0);
                E_x_mu_hist[i] -> Sumw2();
                diff_E2_P2_x_mu_hist[i] = std::make_shared<TH1D>(
                        Form("diff_E2_P2_x_mu_ %4.12s ", NUM_CONDITIONS_NAME[i].c_str()),
                        Form("diff_E2_P2_x_mu %4.12s ", NUM_CONDITIONS_NAME[i].c_str()), bins,
                        -0.2, 0.2);
                diff_E_P_x_mu_hist[i] = std::make_shared<TH1D>(
                        Form("diff_E_P_x_mu_ %4.12s ", NUM_CONDITIONS_NAME[i].c_str()),
                        Form("diff_E_P_x_mu %4.12s)", NUM_CONDITIONS_NAME[i].c_str()), bins,
                        -2.0, 10.0);
                mom_vs_E_x_mu_hist[i] = std::make_shared<TH2D>(
                        Form("mom_vs_E_ %4.12s ", NUM_CONDITIONS_NAME[i].c_str()),
                        Form("mom_vs_E %4.12s ", NUM_CONDITIONS_NAME[i].c_str()), bins, -1.0,
                        10.0, bins, 0.0, 10.0);

                theta_elec_hist[i] = std::make_shared<TH1D>(
                        Form("theta_elec_scattering_ %4.12s ", NUM_CONDITIONS_NAME[i].c_str()),
                        Form("theta_elec_scattering %4.12s ", NUM_CONDITIONS_NAME[i].c_str()),
                        bins, 0.0, 50.0);
                theta_x_mu_hist[i] = std::make_shared<TH1D>(
                        Form("theta_x_mu_ %4.12s ", NUM_CONDITIONS_NAME[i].c_str()),
                        Form("theta_x_mu %4.12s ", NUM_CONDITIONS_NAME[i].c_str()), bins, 0.0,
                        180.0);

                diff_theta_elec_x_mu_hist[i] = std::make_shared<TH1D>(
                        Form("diff_theta_elec_x_mu_ %4.12s ", NUM_CONDITIONS_NAME[i].c_str()),
                        Form("diff_theta_elec_x_mu %4.12s ", NUM_CONDITIONS_NAME[i].c_str()),
                        bins, -50.0, 180.0);

                MM2_VS_W_x_mu_hist[i] = std::make_shared<TH2D>(
                        Form("MM2_vs_W_ %4.12s ", NUM_CONDITIONS_NAME[i].c_str()),
                        Form("MM2_vs_W %4.12s ", NUM_CONDITIONS_NAME[i].c_str()), bins, 0.0,
                        w_max, bins, -0.2, 0.2);
        }
}
void Histogram::Fill_x_mu(const std::shared_ptr<Reaction> &_e) {

        if (_e->MM_cut())
                if (_e->TwoPion_missingPim()) {
                        P_x_mu->Fill(_e->P_x_mu(), _e->weight());
                        //P_x_mu->Sumw2();
                        // diff_theta_in_x_mu->Fill(-_e->theta_beam() + _e->theta_x_mu());
                        E_x_mu_hist[0] -> Fill(_e->E_x_mu(), _e->weight());
                        diff_E2_P2_x_mu_hist[0] -> Fill(_e->M2_x_mu(),
                                                        _e->weight()); //, _e->weight());
                        diff_E_P_x_mu_hist[0] -> Fill(_e->M_x_mu(), _e->weight()); //(_e->M_x_mu());
                        diff_E_P_x_mu_hist_->Fill(_e->M_x_mu(), _e->weight());
                        theta_elec_hist[0] -> Fill(_e->theta_elec(), _e->weight());
                        theta_x_mu_hist[0] -> Fill(_e->theta_x_mu(), _e->weight());
                        mom_vs_E_x_mu_hist[0] -> Fill(_e->P_x_mu(), _e->E_x_mu(), _e->weight());
                        diff_theta_elec_x_mu_hist[0] -> Fill(_e->theta_x_mu() - _e->theta_elec(),
                                                             _e->weight());
                        MM2_VS_W_x_mu_hist[0] -> Fill(_e->W(), _e->M2_x_mu(), _e->weight());
                }

        // else if (_e->missingPim()) {
        //   E_x_mu_hist[1]->Fill(_e->E_x_mu());
        //   diff_E2_P2_x_mu_hist[1]->Fill(_e->M2_x_mu());
        //   diff_E_P_x_mu_hist[1]->Fill(_e->M_x_mu());
        //   mom_vs_E_x_mu_hist[1]->Fill(_e->P_x_mu(), _e->E_x_mu());
        //   theta_x_mu_hist[1]->Fill(_e->theta_x_mu());
        //   theta_elec_hist[1]->Fill(_e->theta_elec());
        //   diff_theta_elec_x_mu_hist[1]->Fill(_e->theta_x_mu() - _e->theta_elec());
        //   MM2_VS_W_x_mu_hist[1]->Fill(_e->W(), _e->M2_x_mu());
        // }

        // if (_e->onePositive_at180()) {
        //   E_x_mu_hist[2]->Fill(_e->E_x_mu());
        //   diff_E2_P2_x_mu_hist[2]->Fill(_e->M2_x_mu());
        //   diff_E_P_x_mu_hist[2]->Fill(_e->M_x_mu());
        //   mom_vs_E_x_mu_hist[2]->Fill(_e->P_x_mu(), _e->E_x_mu());
        //   theta_x_mu_hist[2]->Fill(_e->theta_x_mu());
        //   theta_elec_hist[2]->Fill(_e->theta_elec());
        //   diff_theta_elec_x_mu_hist[2]->Fill(_e->theta_x_mu() - _e->theta_elec());
        //   MM2_VS_W_x_mu_hist[2]->Fill(_e->W(), _e->M2_x_mu());
        // }
        //  if (abs(_e->M2_x_mu()) < 0.05) {
}
void Histogram::write_hist_x_mu() {
        diff_E_P_x_mu_hist_->SetXTitle("E - P (GeV)");
        diff_E_P_x_mu_hist_->Write();
        P_x_mu->SetXTitle("3 Mom comp (GeV)");
        P_x_mu->Write();
        for (short i; i < NUM_CONDITIONS; i ++) {
                E_x_mu_hist[i] -> SetXTitle("Energy comp (GeV)");
                E_x_mu_hist[i] -> Write();
                diff_E2_P2_x_mu_hist[i] -> SetXTitle("E2-P2 (GeV)");
                diff_E2_P2_x_mu_hist[i] -> Write();
                diff_E_P_x_mu_hist[i] -> SetXTitle("E-P (GeV)");
                diff_E_P_x_mu_hist[i] -> Write();
                mom_vs_E_x_mu_hist[i] -> SetXTitle("Energy (GeV)");
                mom_vs_E_x_mu_hist[i] -> SetYTitle("Mom (GeV)");
                mom_vs_E_x_mu_hist[i] -> SetOption("COLZ");
                mom_vs_E_x_mu_hist[i] -> Write();
                theta_elec_hist[i] -> SetXTitle("theta (deg)");
                theta_elec_hist[i] -> Write();
                theta_x_mu_hist[i] -> SetXTitle("theta (deg)");
                theta_x_mu_hist[i] -> Write();
                diff_theta_elec_x_mu_hist[i] -> SetXTitle("theta (deg)");
                diff_theta_elec_x_mu_hist[i] -> Write();
                MM2_VS_W_x_mu_hist[i] -> SetXTitle("W (GeV)");
                MM2_VS_W_x_mu_hist[i] -> SetYTitle("MM2 (GeV)");
                MM2_VS_W_x_mu_hist[i] -> SetOption("COLZ");
                MM2_VS_W_x_mu_hist[i] -> Write();
        }
}

void Histogram::Fill_pi0(const std::shared_ptr<Reaction> &_e) {
        if (_e->pi0_mass() < 0.0001)
                return;
        short sec = _e->sec();
        if (sec > 0 && sec <= 6) {
                mass_pi0_any_event[sec - 1] -> Fill(_e->pi0_mass(), _e->weight());
                //  std::cout << "pi0 mass : " << _e->pi0_mass() << '\n';
                if (_e->TwoPion_missingPim()) {
                        mass_pi0_hist_before_mmsq_cut[sec - 1] -> Fill(_e->pi0_mass(),
                                                                       _e->weight());
                        if (_e->MM_cut()) { // abs(mmsq<0.03
                                mass_pi0_hist_after_mmsq_cut[sec - 1] -> Fill(_e->pi0_mass(),
                                                                              _e->weight());
                        }
                }
        }
}
// void Histogram::Fill_x_mu(const std::shared_ptr<Reaction>& _e) {
//   E_x_mu->Fill(_e->E_x_mu());
//   P_x_mu->Fill(_e->P_x_mu());
//   diff_E_P_x_mu->Fill(-_e->E_x_mu() + _e->P_x_mu());
//   mom_vs_E_x_mu->Fill(_e->E_x_mu(), _e->P_x_mu());
//   diff_theta_in_x_mu->Fill(-_e->theta_beam() + _e->theta_x_mu());
//   Theta_vs_mom_x_mu->Fill(_e->P_x_mu(), _e->theta_x_mu());
// }
// void Histogram::Fill_SF(const std::shared_ptr<Branches12>& _d, int part) {
//   sf_hist->Fill(_d->p(0), _d->ec_tot_energy(0) / _d->p(0));
// }
// void Histogram::Write_SF() {
//   // set stat
//
//   sf_hist->SetOption("COLZ");
//   gStyle->SetOptFit(1111);
//   sf_hist->Write();
//
//   E_x_mu->SetXTitle("Energy comp (GeV)");
//   E_x_mu->Write();
//
//   P_x_mu->SetXTitle("3 momentum comp (GeV)");
//   P_x_mu->Write();
//
//   diff_E_P_x_mu->SetXTitle("diif (E-P_mag) (GeV)");
//   diff_E_P_x_mu->Write();
//
//   diff_theta_in_x_mu->SetXTitle("theta_x_mu - theta_beam");
//   diff_theta_in_x_mu->Write();
//
//   mom_vs_E_x_mu->SetXTitle("E comp (GeV)");
//   mom_vs_E_x_mu->SetYTitle("P comp (GeV)");
//   mom_vs_E_x_mu->SetOption("COLZ");
//   mom_vs_E_x_mu->Write();
//
//   Dthtea_vs_Dphi->SetXTitle("#Delata#Theta (deg)");
//   Dthtea_vs_Dphi->SetYTitle("#Delata#Phi (deg)");
//   Dthtea_vs_Dphi->SetOption("COLZ");
//   Dthtea_vs_Dphi->Write();
// }

void Histogram::makeHistTheta_pim_measured() {
        for (size_t i = 0; i < theta_bin_NUM; i ++) {
                theta_pim_measured_3_sigma[i] = std::make_shared<TH1D>(Form("theta_pim_measured_%1.12s (deg)", theta_bin_NAME[i].c_str()),
                                                                       Form("theta_pim_measured_%1.12s (deg)", theta_bin_NAME[i].c_str()), bins, 0.0, 180);
        }
}
void Histogram::populate_theta_pim_measured(const std::shared_ptr<Reaction> &_e, double min, double max, short index_theta_pim) {
        if (_e->TwoPion_missingPim()) {
                if (_e->MM_cut()) {
                        if(_e->pim_momentum_measured()>0) {
                                if (_e->pim_theta_lab() > min && _e->pim_theta_lab() < max) {
                                        theta_pim_measured_3_sigma[index_theta_pim]->Fill(_e->pim_theta_lab_measured(),_e->weight());

                                }
                        }
                }
        }
}
void Histogram::Fill_theta_pim_measured(const std::shared_ptr<Reaction> &_e) {
        short index_theta_pim = 0;
        // populate_SF(_d, 0., 0.5, 0);
        // if (_d->p(0) > 1. && _d->p(0) < 1.5) SF_1D[2]->Fill(_d->ec_tot_energy(0) / _d->p(0));
        for (float x = 0.; x <= 170; x = x + 10) {
                float y = x + 10;
                populate_theta_pim_measured(_e, x, y, index_theta_pim);
                index_theta_pim++;
        }
}

void Histogram::write_hist_theta_pim_measured() {
        for (size_t i = 0; i < theta_bin_NUM; i++) {
                // theta_pim_measured_3_sigma[i]->SetXTitle("sf (E/P)");
                // theta_pim_measured_3_sigma[i]->Fit("gaus", "QMR+", "QMR+", 0.224, 0.28);
                //     //  gROOT->SetStyle("Plain");
                gStyle->SetOptFit(1111);
                //     // if (SF_1D[i]->GetEntries())
                theta_pim_measured_3_sigma[i]->Write();
        }
}

//
// void Histogram::Fill_deltat_pi(const std::shared_ptr<Branches12> &data,
//                                const std::shared_ptr<Delta_T> &dt, int part) {
//         auto _cuts = std::make_unique<Cuts>(data, dt);
//         int charge = data->charge(part);
//         int status = abs(data->status(part));
//         bool fc = dt->ctof();
//         int pid = data->pid(part);
//         float mom = data->p(part);
//         float time = NAN;
//         if (fc)
//                 time = dt->dt_ctof_Pi();
//         else
//                 time = dt->dt_Pi();
//
//         if (charge == 1) {
//                 delta_t_hist[1][0][0][fc] -> Fill(mom, time);
//                 if (_cuts->IsPip(part))
//                         delta_t_hist[1][0][1][fc] -> Fill(mom, time);
//                 else
//                         delta_t_hist[1][0][2][fc] -> Fill(mom, time);
//         } else if (charge == -1) {
//                 delta_t_hist[1][1][0][fc] -> Fill(mom, time);
//                 //if (_cuts->IsPim(part))
//                 if (pid !=ELECTRON)
//                         delta_t_hist[1][1][1][fc] -> Fill(mom, time);
//                 else
//                         delta_t_hist[1][1][2][fc] -> Fill(mom, time);
//         }
// }
void Histogram::makeHist_efficiency_check() {
        for (size_t i = 0; i < EFF_CONDITIONS_NUM; i++) {
                pim_theta_rec_vs_mom_CD[i] = std::make_shared<TH2D>(
                        Form("pim_theta_rec_vs_mom_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_theta_rec_vs_mom_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        p_max, bins, zero, 140.0);

                pim_phi_rec_vs_mom_CD[i]  = std::make_shared<TH2D>(
                        Form("pim_Phi_rec_vs_mom_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_Phi_rec_vs_mom_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        p_max, bins, zero, 360.0);

                pim_theta_measured_vs_mom_CD[i] = std::make_shared<TH2D>(
                        Form("pim_theta_measured_vs_mom_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_theta_measured_vs_mom_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        p_max, bins, zero, 140.0);

                pim_phi_measured_vs_mom_CD[i]  = std::make_shared<TH2D>(
                        Form("pim_Phi_measured_vs_mom_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_Phi_measured_vs_mom_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        p_max, bins, zero, 360.0);

                pim_phi_vs_theta_rec_CD[i]  = std::make_shared<TH2D>(
                        Form("pim_Phi_vs_theta_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_Phi_vs_theta_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        180.0, bins, zero, 360.0);

                pim_phi_vs_theta_measured_CD[i]  = std::make_shared<TH2D>(
                        Form("pim_Phi_vs_theta_measured_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_Phi_vs_theta_measured_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        180.0, bins, zero, 360.0);

                pim_E_vs_mom_rec_CD[i]  = std::make_shared<TH2D>(
                        Form("pim_E_vs_mom_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_E_vs_mom_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        p_max, bins, zero, p_max);

                pim_E_vs_mom_measured_CD[i]  = std::make_shared<TH2D>(
                        Form("pim_E_vs_mom_measured_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_E_vs_mom_measured_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        p_max, bins, zero, p_max);

                pim_mom_measured_vs_rec_CD[i]  = std::make_shared<TH2D>(
                        Form("pim_mom_measured_vs_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_mom_measured_vs_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        p_max, bins, zero, p_max);

                pim_E_measured_vs_rec_CD[i]  = std::make_shared<TH2D>(
                        Form("pim_E_measured_vs_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_E_measured_vs_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        p_max, bins, zero, p_max);

                pim_theta_measured_vs_rec_CD[i]  = std::make_shared<TH2D>(
                        Form("pim_theta_measured_vs_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_theta_measured_vs_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        180.0, bins, zero, 180.0);

                pim_phi_measured_vs_rec_CD[i]  = std::make_shared<TH2D>(
                        Form("pim_phi_measured_vs_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_phi_measured_vs_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        360.0, bins, zero, 360.0);


                pim_theta_measured_minus_rec_vs_measured_CD[i]  = std::make_shared<TH2D>(
                        Form("pim_theta_measured_minus_rec_vs_measured_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_theta_measured_minus_rec_vs_measured_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, 30.0,
                        140.0, bins, -140.0, 140.0);

                // pim_theta_measured_minus_rec_vs_rec_CD[i]  = std::make_shared<TH2D>(
                //         Form("pim_theta_measured_minus_rec_vs_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                //         Form("pim_theta_measured_minus_rec_vs_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, 30.0,
                //         140.0, bins, -140.0, 140.0);

                pim_theta_measured_vs_thrown_CD[i]  = std::make_shared<TH2D>(
                        Form("pim_theta_measured_vs_thrown_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_theta_measured_vs_thrown_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, 30.0,
                        140.0, bins, -140.0, 140.0);
                pim_theta_rec_vs_thrown_CD[i]  = std::make_shared<TH2D>(
                        Form("pim_theta_rec_vs_thrown_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_theta_rec_vs_thrown_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),  bins, 30.0,
                        140.0, bins, -140.0, 140.0);

                MM2_exclusive_hist_CD[i]  = std::make_shared<TH1D>(
                        Form("MM2_exclusive_hist_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("MM2_exclusive_hist_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, -0.7, 0.7);

                // pim_theta_measured_minus_rec_CD[i] = std::make_shared<TH1D>(
                //         Form("pim_theta_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                //         Form("pim_theta_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins,-20, 20);
                //
                // pim_phi_measured_minus_rec_CD[i] = std::make_shared<TH1D>(
                //         Form("pim_phi_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                //         Form("pim_phi_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins,-10, 10);
                //
                // pim_mom_measured_minus_rec_CD[i] = std::make_shared<TH1D>(
                //         Form("pim_mom_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                //         Form("pim_mom_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins,-0.4, 0.4);




                pim_theta_rec_vs_mom_FD[i] = std::make_shared<TH2D>(
                        Form("pim_theta_rec_vs_mom_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_theta_rec_vs_mom_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        p_max, bins, zero, 60.0);

                pim_phi_rec_vs_mom_FD[i]  = std::make_shared<TH2D>(
                        Form("pim_Phi_rec_vs_mom_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_Phi_rec_vs_mom_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        p_max, bins, zero, 360.0);

                pim_theta_measured_vs_mom_FD[i] = std::make_shared<TH2D>(
                        Form("pim_theta_measured_vs_mom_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_theta_measured_vs_mom_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        p_max, bins, zero, 60.0);

                pim_phi_measured_vs_mom_FD[i]  = std::make_shared<TH2D>(
                        Form("pim_Phi_measured_vs_mom_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_Phi_measured_vs_mom_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        p_max, bins, zero, 360.0);

                pim_phi_vs_theta_rec_FD[i]  = std::make_shared<TH2D>(
                        Form("pim_Phi_vs_theta_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_Phi_vs_theta_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        60.0, bins, zero, 360.0);

                pim_phi_vs_theta_measured_FD[i]  = std::make_shared<TH2D>(
                        Form("pim_Phi_vs_theta_measured_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_Phi_vs_theta_measured_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        60.0, bins, zero, 360.0);

                pim_E_vs_mom_rec_FD[i]  = std::make_shared<TH2D>(
                        Form("pim_E_vs_mom_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_E_vs_mom_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        p_max, bins, zero, p_max);

                pim_E_vs_mom_measured_FD[i]  = std::make_shared<TH2D>(
                        Form("pim_E_vs_mom_measured_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_E_vs_mom_measured_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        p_max, bins, zero, p_max);

                pim_mom_measured_vs_rec_FD[i]  = std::make_shared<TH2D>(
                        Form("pim_mom_measured_vs_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_mom_measured_vs_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        p_max, bins, zero, p_max);

                pim_E_measured_vs_rec_FD[i]  = std::make_shared<TH2D>(
                        Form("pim_E_measured_vs_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_E_measured_vs_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        p_max, bins, zero, p_max);

                pim_theta_measured_vs_rec_FD[i]  = std::make_shared<TH2D>(
                        Form("pim_theta_measured_vs_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_theta_measured_vs_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        60.0, bins, zero, 60.0);

                pim_phi_measured_vs_rec_FD[i]  = std::make_shared<TH2D>(
                        Form("pim_phi_measured_vs_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_phi_measured_vs_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        360.0, bins, zero, 360.0);

                pim_theta_measured_minus_rec_vs_measured_FD[i]  = std::make_shared<TH2D>(
                        Form("pim_theta_measured_minus_rec_vs_measured_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_theta_measured_minus_rec_vs_measured_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        60.0, bins, -60.0, 60.0);

                // pim_theta_measured_minus_rec_vs_rec_FD[i]  = std::make_shared<TH2D>(
                //         Form("pim_theta_measured_minus_rec_vs_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                //         Form("pim_theta_measured_minus_rec_vs_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                //         60.0, bins, -60, 60.0);

                pim_theta_measured_vs_thrown_FD[i]  = std::make_shared<TH2D>(
                        Form("pim_theta_measured_vs_thrown_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_theta_measured_vs_thrown_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        60.0, bins, -60, 60.0);

                pim_theta_rec_vs_thrown_FD[i]  = std::make_shared<TH2D>(
                        Form("pim_theta_rec_vs_thrown_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("pim_theta_rec_vs_thrown_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, zero,
                        60.0, bins, -60, 60.0);

                MM2_exclusive_hist_FD[i]  = std::make_shared<TH1D>(
                        Form("MM2_exclusive_hist_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                        Form("MM2_exclusive_hist_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins, -0.7, 0.7);

                // pim_theta_measured_minus_rec_FD[i] = std::make_shared<TH1D>(
                //         Form("pim_theta_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                //         Form("pim_theta_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins,-7, 7);
                //
                // pim_phi_measured_minus_rec_FD[i] = std::make_shared<TH1D>(
                //         Form("pim_phi_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                //         Form("pim_phi_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins,-10, 10);
                //
                // pim_mom_measured_minus_rec_FD[i] = std::make_shared<TH1D>(
                //         Form("pim_mom_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()),
                //         Form("pim_mom_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME[i].c_str()), bins,-0.4, 0.4);

                //
                // // missingPiP

                pip_theta_rec_vs_mom[i] = std::make_shared<TH2D>(
                        Form("pip_theta_rec_vs_mom_%6.20s ", EFF_CONDITIONS_NAME_ALL[i].c_str()),
                        Form("pip_theta_rec_vs_mom_%6.20s ", EFF_CONDITIONS_NAME_ALL[i].c_str()), bins, zero,
                        p_max, bins, zero, 140.0);

                pip_theta_measured_vs_mom[i] = std::make_shared<TH2D>(
                        Form("pip_theta_measured_vs_mom_%6.20s ", EFF_CONDITIONS_NAME_ALL[i].c_str()),
                        Form("pip_theta_measured_vs_mom_%6.20s ", EFF_CONDITIONS_NAME_ALL[i].c_str()), bins, zero,
                        p_max, bins, zero, 140.0);
                //
                // pip_theta_measured_minus_rec_CD[i] = std::make_shared<TH1D>(
                //         Form("pip_theta_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME_PIP[i].c_str()),
                //         Form("pip_theta_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME_PIP[i].c_str()), bins,-30, 30);
                //
                // pip_phi_measured_minus_rec_CD[i] = std::make_shared<TH1D>(
                //         Form("pip_phi_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME_PIP[i].c_str()),
                //         Form("pip_phi_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME_PIP[i].c_str()), bins,-10, 10);
                //
                // pip_mom_measured_minus_rec_CD[i] = std::make_shared<TH1D>(
                //         Form("pip_mom_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME_PIP[i].c_str()),
                //         Form("pip_mom_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME_PIP[i].c_str()), bins,-0.4, 0.4);
                //
                //
                // pip_theta_measured_minus_rec_FD[i] = std::make_shared<TH1D>(
                //         Form("pip_theta_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME_PIP[i].c_str()),
                //         Form("pip_theta_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME_PIP[i].c_str()), bins,-6, 6);
                //
                // pip_phi_measured_minus_rec_FD[i] = std::make_shared<TH1D>(
                //         Form("pip_phi_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME_PIP[i].c_str()),
                //         Form("pip_phi_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME_PIP[i].c_str()), bins,-10, 10);
                //
                // pip_mom_measured_minus_rec_FD[i] = std::make_shared<TH1D>(
                //         Form("pip_mom_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME_PIP[i].c_str()),
                //         Form("pip_mom_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME_PIP[i].c_str()), bins,-0.4, 0.4);
                //
                //
                //
                // // missingProt
                prot_theta_rec_vs_mom[i] = std::make_shared<TH2D>(
                        Form("prot_theta_rec_vs_mom_%6.20s ", EFF_CONDITIONS_NAME_ALL[i].c_str()),
                        Form("prot_theta_rec_vs_mom_%6.20s ", EFF_CONDITIONS_NAME_ALL[i].c_str()), bins, zero,
                        p_max, bins, zero, 140.0);

                prot_theta_measured_vs_mom[i] = std::make_shared<TH2D>(
                        Form("prot_theta_measured_vs_mom_%6.20s ", EFF_CONDITIONS_NAME_ALL[i].c_str()),
                        Form("prot_theta_measured_vs_mom_%6.20s ", EFF_CONDITIONS_NAME_ALL[i].c_str()), bins, zero,
                        p_max, bins, zero, 140.0);
                // prot_theta_measured_minus_rec_CD[i] = std::make_shared<TH1D>(
                //         Form("prot_theta_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME_PROT[i].c_str()),
                //         Form("prot_theta_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME_PROT[i].c_str()), bins,-20, 20);
                //
                // prot_phi_measured_minus_rec_CD[i] = std::make_shared<TH1D>(
                //         Form("prot_phi_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME_PROT[i].c_str()),
                //         Form("prot_phi_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME_PROT[i].c_str()), bins,-10, 10);
                //
                // prot_mom_measured_minus_rec_CD[i] = std::make_shared<TH1D>(
                //         Form("prot_mom_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME_PROT[i].c_str()),
                //         Form("prot_mom_measured_minus_rec_CD_%6.20s ", EFF_CONDITIONS_NAME_PROT[i].c_str()), bins,-0.4, 0.4);
                //
                //
                // prot_theta_measured_minus_rec_FD[i] = std::make_shared<TH1D>(
                //         Form("prot_theta_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME_PROT[i].c_str()),
                //         Form("prot_theta_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME_PROT[i].c_str()), bins,-10, 10);
                //
                // prot_phi_measured_minus_rec_FD[i] = std::make_shared<TH1D>(
                //         Form("prot_phi_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME_PROT[i].c_str()),
                //         Form("prot_phi_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME_PROT[i].c_str()), bins,-10, 10);
                //
                // prot_mom_measured_minus_rec_FD[i] = std::make_shared<TH1D>(
                //         Form("prot_mom_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME_PROT[i].c_str()),
                //         Form("prot_mom_measured_minus_rec_FD_%6.20s ", EFF_CONDITIONS_NAME_PROT[i].c_str()), bins,-0.4, 0.4);

        }
}

void Histogram::makeHist_efficiency_check_1d() {
        for (size_t h = 0; h < HADRON_NUM; h++) {
                for (size_t e = 0; e < EFF_CONDITIONS_NUM_ALL; e++) {
                        for (size_t pt = 0; pt < DETECTOR_NUM_PROT; pt++) {
                                for (size_t pp = 0; pp < DETECTOR_NUM_PIP; pp++) {
                                        for (size_t pm = 0; pm < DETECTOR_NUM_PIM; pm++) {
                                                mom_measured_minus_rec[h][e][pt][pp][pm] = std::make_shared<TH1D>(
                                                        Form("mom measured - rec %s_%s_%s_%s_%s ", HADRON_NAME[h].c_str(), EFF_CONDITIONS_NAME_ALL[e].c_str(),DETECTOR_NAME_PROT[pt].c_str(),DETECTOR_NAME_PIP[pp].c_str(),DETECTOR_NAME_PIM[pm].c_str()),
                                                        Form("mom measured - rec %s_%s_%s_%s_%s ", HADRON_NAME[h].c_str(), EFF_CONDITIONS_NAME_ALL[e].c_str(),DETECTOR_NAME_PROT[pt].c_str(),DETECTOR_NAME_PIP[pp].c_str(),DETECTOR_NAME_PIM[pm].c_str()),
                                                        bins,-0.4, 0.4);

                                                theta_measured_minus_rec[h][e][pt][pp][pm] = std::make_shared<TH1D>(
                                                        Form("#theta measured - rec %s_%s_%s_%s_%s ", HADRON_NAME[h].c_str(), EFF_CONDITIONS_NAME_ALL[e].c_str(),DETECTOR_NAME_PROT[pt].c_str(),DETECTOR_NAME_PIP[pp].c_str(),DETECTOR_NAME_PIM[pm].c_str()),
                                                        Form("#theta measured - rec %s_%s_%s_%s_%s ", HADRON_NAME[h].c_str(), EFF_CONDITIONS_NAME_ALL[e].c_str(),DETECTOR_NAME_PROT[pt].c_str(),DETECTOR_NAME_PIP[pp].c_str(),DETECTOR_NAME_PIM[pm].c_str()),
                                                        bins, -20, 20);
                                                phi_measured_minus_rec[h][e][pt][pp][pm] = std::make_shared<TH1D>(
                                                        Form("#phi measured - rec %s_%s_%s_%s_%s ", HADRON_NAME[h].c_str(), EFF_CONDITIONS_NAME_ALL[e].c_str(),DETECTOR_NAME_PROT[pt].c_str(),DETECTOR_NAME_PIP[pp].c_str(),DETECTOR_NAME_PIM[pm].c_str()),
                                                        Form("#phi measured - rec %s_%s_%s_%s_%s ", HADRON_NAME[h].c_str(), EFF_CONDITIONS_NAME_ALL[e].c_str(),DETECTOR_NAME_PROT[pt].c_str(),DETECTOR_NAME_PIP[pp].c_str(),DETECTOR_NAME_PIM[pm].c_str()),
                                                        bins, -10, 10);

                                                // theta_rec_vs_mom[h][e][pt][pp][pm] = std::make_shared<TH2D>(
                                                //         Form("theta_rec_vs_mom %s_%s_%s_%s_%s ", HADRON_NAME[h].c_str(), EFF_CONDITIONS_NAME_ALL[e].c_str(),DETECTOR_NAME_PROT[pt].c_str(),DETECTOR_NAME_PIP[pp].c_str(),DETECTOR_NAME_PIM[pm].c_str()),
                                                //         Form("theta_rec_vs_mom %s_%s_%s_%s_%s ", HADRON_NAME[h].c_str(), EFF_CONDITIONS_NAME_ALL[e].c_str(),DETECTOR_NAME_PROT[pt].c_str(),DETECTOR_NAME_PIP[pp].c_str(),DETECTOR_NAME_PIM[pm].c_str()),
                                                //         bins, zero, p_max, bins, zero, 140.0);

                                        }
                                }
                        }
                }
        }
}
void Histogram::Fill_efficiency_check_CD_rec(const std::shared_ptr<Reaction> &_e){
        if(_e->pim_momentum()!= NAN) {

                if (_e->TwoPion_missingPim()) {
                        MM2_exclusive_hist_CD[0]->Fill(_e->MM2(),_e->weight());
                        if (_e->MM_cut()) {
                                pim_theta_rec_vs_mom_CD[0]->Fill(_e->pim_momentum(),_e->pim_theta_lab(),_e->weight());
                                pim_phi_rec_vs_mom_CD[0]->Fill(_e->pim_momentum(),_e->pim_Phi_lab(),_e->weight());
                                pim_phi_vs_theta_rec_CD[0]->Fill(_e->pim_theta_lab(),_e->pim_Phi_lab(),_e->weight());
                                pim_E_vs_mom_rec_CD[0]->Fill(_e->pim_momentum(),_e->pim_E(),_e->weight());
                                if (_e->TwoPion_exclusive())  {
                                        MM2_exclusive_hist_CD[1]->Fill(_e->MM2(),_e->weight());

                                        pim_theta_rec_vs_mom_CD[1]->Fill(_e->pim_momentum(),_e->pim_theta_lab(),_e->weight());
                                        pim_phi_rec_vs_mom_CD[1]->Fill(_e->pim_momentum(),_e->pim_Phi_lab(),_e->weight());
                                        pim_phi_vs_theta_rec_CD[1]->Fill(_e->pim_theta_lab(),_e->pim_Phi_lab(),_e->weight());
                                        pim_E_vs_mom_rec_CD[1]->Fill(_e->pim_momentum(),_e->pim_E(),_e->weight());
                                }
                        }
                }
        }
}
// void Histogram::Fill_efficiency_check_CD_thrown(const std::shared_ptr<MCReaction> &_mc_e,const std::shared_ptr<Reaction> &_e) {
// //        short sec = _e->sec();
//         if (_e->TwoPion_missingPim()) {
//                 if(_e->pim_momentum_measured()!= NAN) {
//                         pim_theta_measured_vs_thrown_CD[0]->Fill(_mc_e->MCpim_theta_lab(), (_e->pim_theta_lab_measured()-_mc_e->MCpim_theta_lab()),_e->weight());
//                         if (_e->MM_cut()) {
//                                 pim_theta_measured_vs_thrown_CD[1]->Fill(_mc_e->MCpim_theta_lab(), (_e->pim_theta_lab_measured()-_mc_e->MCpim_theta_lab()),_e->weight());
//                         }
//                         if(abs(_mc_e->MCpim_theta_lab() -_e->pim_theta_lab_measured()) < 10 )
//                                 pim_theta_measured_vs_thrown_CD[2]->Fill(_mc_e->MCpim_theta_lab(), (_e->pim_theta_lab_measured()-_mc_e->MCpim_theta_lab()),_e->weight());
//                 }
//         }
// }

void Histogram::Fill_efficiency_check_1d_hist_prot(const std::shared_ptr<Reaction> &_e){


        if (_e->TwoPion_missingProt()) {
                if (_e->MM2_mprot()>0.7 &&_e->MM2_mprot()<1.2) {
                        if(_e->reconstructed_prot_mom_lab() > 0.0)
                                prot_theta_rec_vs_mom[0]->Fill(_e->reconstructed_prot_mom_lab(),_e->reconstructed_prot_theta_lab(), _e->weight());
                        if(_e->measured_prot_mom_lab_no_missing() > 0.0)
                                prot_theta_measured_vs_mom[0]->Fill(_e->measured_prot_mom_lab_no_missing(),_e->measured_prot_theta_lab_no_missing(), _e->weight());

                        if(_e->ctof_prot()) {
                                if(_e->ctof_pip()) {
                                        if(_e->ctof_pim()) {
                                                mom_measured_minus_rec[0][0][0][0][0] -> Fill(_e->diff_prot_momentum(),_e->weight());
                                                theta_measured_minus_rec[0][0][0][0][0] -> Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                phi_measured_minus_rec[0][0][0][0][0] -> Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                //theta_rec_vs_mom[0][0][0][0][0]->Fill(_e->diff_prot_momentum(),_e->weight());

                                        }
                                        else if(_e->ftof_pim()) {
                                                mom_measured_minus_rec[0][0][0][0][1] -> Fill(_e->diff_prot_momentum(),_e->weight());
                                                theta_measured_minus_rec[0][0][0][0][1] -> Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                phi_measured_minus_rec[0][0][0][0][1] -> Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                //theta_rec_vs_mom[0][0][0][0][1]->Fill(_e->diff_prot_momentum(),_e->weight());
                                        }
                                }
                                else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                  mom_measured_minus_rec[0][0][0][1][0] -> Fill(_e->diff_prot_momentum(),_e->weight());
                                                                  theta_measured_minus_rec[0][0][0][1][0] -> Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                                  phi_measured_minus_rec[0][0][0][1][0] -> Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                                  //theta_rec_vs_mom[0][0][0][1][0]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                          }
                                                          else if(_e->ftof_pim()) {
                                                                  mom_measured_minus_rec[0][0][0][1][1] -> Fill(_e->diff_prot_momentum(),_e->weight());
                                                                  theta_measured_minus_rec[0][0][0][1][1] -> Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                                  phi_measured_minus_rec[0][0][0][1][1] -> Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                                  //theta_rec_vs_mom[0][0][0][1][1]->Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                          }}
                        }
                        else if(_e->ftof_prot()) {
                                if(_e->ctof_pip()) {
                                        if(_e->ctof_pim()) {
                                                mom_measured_minus_rec[0][0][1][0][0] -> Fill(_e->diff_prot_momentum(),_e->weight());
                                                phi_measured_minus_rec[0][0][1][0][0] -> Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                //theta_rec_vs_mom[0][0][1][0][0]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                theta_measured_minus_rec[0][0][1][0][0] -> Fill(_e->diff_prot_theta_lab(),_e->weight());

                                        }
                                        else if(_e->ftof_pim()) {
                                                mom_measured_minus_rec[0][0][1][0][1] -> Fill(_e->diff_prot_momentum(),_e->weight());
                                                phi_measured_minus_rec[0][0][1][0][1] -> Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                //theta_rec_vs_mom[0][0][1][0][1]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                theta_measured_minus_rec[0][0][1][0][1] -> Fill(_e->diff_prot_theta_lab(),_e->weight());
                                        }
                                }
                                else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                  mom_measured_minus_rec[0][0][1][1][0] -> Fill(_e->diff_prot_momentum(),_e->weight());
                                                                  phi_measured_minus_rec[0][0][1][1][0] -> Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                                  //theta_rec_vs_mom[0][0][1][1][0]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                  theta_measured_minus_rec[0][0][1][1][0] -> Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                          }
                                                          else if(_e->ftof_pim()) {
                                                                  mom_measured_minus_rec[0][0][1][1][1] -> Fill(_e->diff_prot_momentum(),_e->weight());
                                                                  phi_measured_minus_rec[0][0][1][1][1] -> Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                                  //theta_rec_vs_mom[0][0][1][1][0]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                  theta_measured_minus_rec[0][0][1][1][1] -> Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                          }}
                        }

                        if (_e->TwoPion_exclusive())  {
                                if(_e->reconstructed_prot_mom_lab() > 0.0)
                                        prot_theta_rec_vs_mom[1] -> Fill(_e->reconstructed_prot_mom_lab(),_e->reconstructed_prot_theta_lab(), _e->weight());
                                if(_e->measured_prot_mom_lab_no_missing() > 0.0)
                                        prot_theta_measured_vs_mom[1] -> Fill(_e->measured_prot_mom_lab_no_missing(),_e->measured_prot_theta_lab_no_missing(), _e->weight());

                                if(_e->ctof_prot()) {
                                        if(_e->ctof_pip()) {
                                                if(_e->ctof_pim()) {
                                                        mom_measured_minus_rec[0][1][0][0][0] -> Fill(_e->diff_prot_momentum(),_e->weight());
                                                        theta_measured_minus_rec[0][1][0][0][0] -> Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                        phi_measured_minus_rec[0][1][0][0][0] -> Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                        //theta_rec_vs_mom[0][1][0][0][0]->Fill(_e->diff_prot_momentum(),_e->weight());

                                                }
                                                else if(_e->ftof_pim()) {
                                                        mom_measured_minus_rec[0][1][0][0][1] -> Fill(_e->diff_prot_momentum(),_e->weight());
                                                        theta_measured_minus_rec[0][1][0][0][1] -> Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                        phi_measured_minus_rec[0][1][0][0][1] -> Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                        //theta_rec_vs_mom[0][1][0][0][1]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                }
                                        }
                                        else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                          mom_measured_minus_rec[0][1][0][1][0] -> Fill(_e->diff_prot_momentum(),_e->weight());
                                                                          theta_measured_minus_rec[0][1][0][1][0] -> Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                                          phi_measured_minus_rec[0][1][0][1][0] -> Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                                          //theta_rec_vs_mom[0][1][0][1][0]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                  }
                                                                  else if(_e->ftof_pim()) {
                                                                          mom_measured_minus_rec[0][1][0][1][1] -> Fill(_e->diff_prot_momentum(),_e->weight());
                                                                          theta_measured_minus_rec[0][1][0][1][1] -> Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                                          phi_measured_minus_rec[0][1][0][1][1] -> Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                                          //theta_rec_vs_mom[0][1][0][1][1]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                  }}
                                }
                                else if(_e->ftof_prot()) {
                                        if(_e->ctof_pip()) {
                                                if(_e->ctof_pim()) {
                                                        mom_measured_minus_rec[0][1][1][0][0] -> Fill(_e->diff_prot_momentum(),_e->weight());
                                                        phi_measured_minus_rec[0][1][1][0][0] -> Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                        //theta_rec_vs_mom[0][1][1][0][0]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                        theta_measured_minus_rec[0][1][1][0][0] -> Fill(_e->diff_prot_theta_lab(),_e->weight());

                                                }
                                                else if(_e->ftof_pim()) {
                                                        mom_measured_minus_rec[0][1][1][0][1] -> Fill(_e->diff_prot_momentum(),_e->weight());
                                                        phi_measured_minus_rec[0][1][1][0][1] -> Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                        //theta_rec_vs_mom[0][1][1][0][1]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                        theta_measured_minus_rec[0][1][1][0][1] -> Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                }
                                        }
                                        else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                          mom_measured_minus_rec[0][1][1][1][0] -> Fill(_e->diff_prot_momentum(),_e->weight());
                                                                          phi_measured_minus_rec[0][1][1][1][0] -> Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                                          //theta_rec_vs_mom[0][1][1][1][0]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                          theta_measured_minus_rec[0][1][1][1][0] -> Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                                  }
                                                                  else if(_e->ftof_pim()) {
                                                                          mom_measured_minus_rec[0][1][1][1][1] -> Fill(_e->diff_prot_momentum(),_e->weight());
                                                                          phi_measured_minus_rec[0][1][1][1][1] -> Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                                          //theta_rec_vs_mom[0][1][1][1][0]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                          theta_measured_minus_rec[0][1][1][1][1] -> Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                                  }}

                                }
                                if(abs(_e->MM2_exclusive()) < 0.03 ) {
                                        if(_e->reconstructed_prot_mom_lab()>0)
                                                prot_theta_rec_vs_mom[2] -> Fill(_e->reconstructed_prot_mom_lab(),_e->reconstructed_prot_theta_lab(), _e->weight());
                                        if(_e->measured_prot_mom_lab_no_missing() > 0.0)
                                                prot_theta_measured_vs_mom[2] -> Fill(_e->measured_prot_mom_lab_no_missing(),_e->measured_prot_theta_lab_no_missing(), _e->weight());

                                        if(_e->ctof_prot()) {
                                                if(_e->ctof_pip()) {
                                                        if(_e->ctof_pim()) {
                                                                mom_measured_minus_rec[0][2][0][0][0]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                theta_measured_minus_rec[0][2][0][0][0]->Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                                phi_measured_minus_rec[0][2][0][0][0]->Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                                //theta_rec_vs_mom[0][2][0][0][0]->Fill(_e->diff_prot_momentum(),_e->weight());

                                                        }
                                                        else if(_e->ftof_pim()) {
                                                                mom_measured_minus_rec[0][2][0][0][1]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                theta_measured_minus_rec[0][2][0][0][1]->Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                                phi_measured_minus_rec[0][2][0][0][1]->Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                                //theta_rec_vs_mom[0][2][0][0][1]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                        }
                                                }
                                                else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                                  mom_measured_minus_rec[0][2][0][1][0]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                                  theta_measured_minus_rec[0][2][0][1][0]->Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                                                  phi_measured_minus_rec[0][2][0][1][0]->Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                                                  //theta_rec_vs_mom[0][2][0][1][0]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                          }
                                                                          else if(_e->ftof_pim()) {
                                                                                  mom_measured_minus_rec[0][2][0][1][1]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                                  theta_measured_minus_rec[0][2][0][1][1]->Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                                                  phi_measured_minus_rec[0][2][0][1][1]->Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                                                  //theta_rec_vs_mom[0][2][0][1][1]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                          }}
                                        }
                                        else if(_e->ftof_prot()) {
                                                if(_e->ctof_pip()) {
                                                        if(_e->ctof_pim()) {
                                                                mom_measured_minus_rec[0][2][1][0][0]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                phi_measured_minus_rec[0][2][1][0][0]->Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                                //theta_rec_vs_mom[0][2][1][0][0]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                theta_measured_minus_rec[0][2][1][0][0]->Fill(_e->diff_prot_theta_lab(),_e->weight());

                                                        }
                                                        else if(_e->ftof_pim()) {
                                                                mom_measured_minus_rec[0][2][1][0][1]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                phi_measured_minus_rec[0][2][1][0][1]->Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                                //theta_rec_vs_mom[0][2][1][0][1]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                theta_measured_minus_rec[0][2][1][0][1]->Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                        }
                                                }
                                                else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                                  mom_measured_minus_rec[0][2][1][1][0]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                                  phi_measured_minus_rec[0][2][1][1][0]->Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                                                  //theta_rec_vs_mom[0][2][1][1][0]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                                  theta_measured_minus_rec[0][2][1][1][0]->Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                                          }
                                                                          else if(_e->ftof_pim()) {
                                                                                  mom_measured_minus_rec[0][2][1][1][1]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                                  phi_measured_minus_rec[0][2][1][1][1]->Fill(_e->diff_prot_Phi_lab(),_e->weight());
                                                                                  //theta_rec_vs_mom[0][2][1][1][0]->Fill(_e->diff_prot_momentum(),_e->weight());
                                                                                  theta_measured_minus_rec[0][2][1][1][1]->Fill(_e->diff_prot_theta_lab(),_e->weight());
                                                                          }}
                                        }
                                }


                        }
                }
        }
}
void Histogram::Fill_efficiency_check_1d_hist_pip(const std::shared_ptr<Reaction> &_e){

        if(_e->TwoPion_missingPip()) {
                if (_e->MM2_mpip()>-0.06 &&_e->MM2_mpip()<0.08) {
                        if(_e->reconstructed_pip_mom_lab() > 0.0)
                                pip_theta_rec_vs_mom[0]->Fill(_e->reconstructed_pip_mom_lab(),_e->reconstructed_pip_theta_lab(), _e->weight());
                        if(_e->measured_pip_mom_lab_no_missing() > 0.0)
                                pip_theta_measured_vs_mom[0]->Fill(_e->measured_pip_mom_lab_no_missing(),_e->measured_pip_theta_lab_no_missing(), _e->weight());

                        if(_e->ctof_prot()) {
                                if(_e->ctof_pip()) {
                                        if(_e->ctof_pim()) {
                                                mom_measured_minus_rec[1][0][0][0][0] -> Fill(_e->diff_pip_momentum(),_e->weight());
                                                theta_measured_minus_rec[1][0][0][0][0] -> Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                phi_measured_minus_rec[1][0][0][0][0] -> Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                //theta_rec_vs_mom[0][0][0][0][0]->Fill(_e->diff_pip_momentum(),_e->weight());

                                        }
                                        else if(_e->ftof_pim()) {
                                                mom_measured_minus_rec[1][0][0][0][1] -> Fill(_e->diff_pip_momentum(),_e->weight());
                                                theta_measured_minus_rec[1][0][0][0][1] -> Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                phi_measured_minus_rec[1][0][0][0][1] -> Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                //theta_rec_vs_mom[0][0][0][0][1]->Fill(_e->diff_pip_momentum(),_e->weight());
                                        }
                                }
                                else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                  mom_measured_minus_rec[1][0][0][1][0] -> Fill(_e->diff_pip_momentum(),_e->weight());
                                                                  theta_measured_minus_rec[1][0][0][1][0] -> Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                                  phi_measured_minus_rec[1][0][0][1][0] -> Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                                  //theta_rec_vs_mom[0][0][0][1][0]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                          }
                                                          else if(_e->ftof_pim()) {
                                                                  mom_measured_minus_rec[1][0][0][1][1] -> Fill(_e->diff_pip_momentum(),_e->weight());
                                                                  theta_measured_minus_rec[1][0][0][1][1] -> Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                                  phi_measured_minus_rec[1][0][0][1][1] -> Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                                  //theta_rec_vs_mom[0][0][0][1][1]->Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                          }}
                        }
                        else if(_e->ftof_prot()) {
                                if(_e->ctof_pip()) {
                                        if(_e->ctof_pim()) {
                                                mom_measured_minus_rec[1][0][1][0][0] -> Fill(_e->diff_pip_momentum(),_e->weight());
                                                phi_measured_minus_rec[1][0][1][0][0] -> Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                //theta_rec_vs_mom[0][0][1][0][0]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                theta_measured_minus_rec[1][0][1][0][0] -> Fill(_e->diff_pip_theta_lab(),_e->weight());

                                        }
                                        else if(_e->ftof_pim()) {
                                                mom_measured_minus_rec[1][0][1][0][1] -> Fill(_e->diff_pip_momentum(),_e->weight());
                                                phi_measured_minus_rec[1][0][1][0][1] -> Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                //theta_rec_vs_mom[0][0][1][0][1]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                theta_measured_minus_rec[1][0][1][0][1] -> Fill(_e->diff_pip_theta_lab(),_e->weight());
                                        }
                                }
                                else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                  mom_measured_minus_rec[1][0][1][1][0] -> Fill(_e->diff_pip_momentum(),_e->weight());
                                                                  phi_measured_minus_rec[1][0][1][1][0] -> Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                                  //theta_rec_vs_mom[0][0][1][1][0]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                  theta_measured_minus_rec[1][0][1][1][0] -> Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                          }
                                                          else if(_e->ftof_pim()) {
                                                                  mom_measured_minus_rec[1][0][1][1][1] -> Fill(_e->diff_pip_momentum(),_e->weight());
                                                                  phi_measured_minus_rec[1][0][1][1][1] -> Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                                  //theta_rec_vs_mom[0][0][1][1][0]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                  theta_measured_minus_rec[1][0][1][1][1] -> Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                          }}
                        }

                        if (_e->TwoPion_exclusive())  {
                                pip_theta_rec_vs_mom[1] -> Fill(_e->reconstructed_pip_mom_lab(),_e->reconstructed_pip_theta_lab(), _e->weight());
                                if(_e->measured_pip_mom_lab_no_missing() > 0.0)
                                        pip_theta_measured_vs_mom[1] -> Fill(_e->measured_pip_mom_lab_no_missing(),_e->measured_pip_theta_lab_no_missing(), _e->weight());

                                if(_e->ctof_prot()) {
                                        if(_e->ctof_pip()) {
                                                if(_e->ctof_pim()) {
                                                        mom_measured_minus_rec[1][1][0][0][0] -> Fill(_e->diff_pip_momentum(),_e->weight());
                                                        theta_measured_minus_rec[1][1][0][0][0] -> Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                        phi_measured_minus_rec[1][1][0][0][0] -> Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                        //theta_rec_vs_mom[0][1][0][0][0]->Fill(_e->diff_pip_momentum(),_e->weight());

                                                }
                                                else if(_e->ftof_pim()) {
                                                        mom_measured_minus_rec[1][1][0][0][1] -> Fill(_e->diff_pip_momentum(),_e->weight());
                                                        theta_measured_minus_rec[1][1][0][0][1] -> Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                        phi_measured_minus_rec[1][1][0][0][1] -> Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                        //theta_rec_vs_mom[0][1][0][0][1]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                }
                                        }
                                        else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                          mom_measured_minus_rec[1][1][0][1][0] -> Fill(_e->diff_pip_momentum(),_e->weight());
                                                                          theta_measured_minus_rec[1][1][0][1][0] -> Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                                          phi_measured_minus_rec[1][1][0][1][0] -> Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                                          //theta_rec_vs_mom[0][1][0][1][0]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                  }
                                                                  else if(_e->ftof_pim()) {
                                                                          mom_measured_minus_rec[1][1][0][1][1] -> Fill(_e->diff_pip_momentum(),_e->weight());
                                                                          theta_measured_minus_rec[1][1][0][1][1] -> Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                                          phi_measured_minus_rec[1][1][0][1][1] -> Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                                          //theta_rec_vs_mom[0][1][0][1][1]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                  }}
                                }
                                else if(_e->ftof_prot()) {
                                        if(_e->ctof_pip()) {
                                                if(_e->ctof_pim()) {
                                                        mom_measured_minus_rec[1][1][1][0][0] -> Fill(_e->diff_pip_momentum(),_e->weight());
                                                        phi_measured_minus_rec[1][1][1][0][0] -> Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                        //theta_rec_vs_mom[0][1][1][0][0]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                        theta_measured_minus_rec[1][1][1][0][0] -> Fill(_e->diff_pip_theta_lab(),_e->weight());

                                                }
                                                else if(_e->ftof_pim()) {
                                                        mom_measured_minus_rec[1][1][1][0][1] -> Fill(_e->diff_pip_momentum(),_e->weight());
                                                        phi_measured_minus_rec[1][1][1][0][1] -> Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                        //theta_rec_vs_mom[0][1][1][0][1]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                        theta_measured_minus_rec[1][1][1][0][1] -> Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                }
                                        }
                                        else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                          mom_measured_minus_rec[1][1][1][1][0] -> Fill(_e->diff_pip_momentum(),_e->weight());
                                                                          phi_measured_minus_rec[1][1][1][1][0] -> Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                                          //theta_rec_vs_mom[0][1][1][1][0]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                          theta_measured_minus_rec[1][1][1][1][0] -> Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                                  }
                                                                  else if(_e->ftof_pim()) {
                                                                          mom_measured_minus_rec[1][1][1][1][1] -> Fill(_e->diff_pip_momentum(),_e->weight());
                                                                          phi_measured_minus_rec[1][1][1][1][1] -> Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                                          //theta_rec_vs_mom[0][1][1][1][0]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                          theta_measured_minus_rec[1][1][1][1][1] -> Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                                  }}
                                }

                                if(abs(_e->MM2_exclusive()) < 0.03 ) {
                                        pip_theta_rec_vs_mom[2] -> Fill(_e->reconstructed_pip_mom_lab(),_e->reconstructed_pip_theta_lab(), _e->weight());
                                        if(_e->measured_pip_mom_lab_no_missing() > 0.0)
                                                pip_theta_measured_vs_mom[2] -> Fill(_e->measured_pip_mom_lab_no_missing(),_e->measured_pip_theta_lab_no_missing(), _e->weight());

                                        if(_e->ctof_prot()) {
                                                if(_e->ctof_pip()) {
                                                        if(_e->ctof_pim()) {
                                                                mom_measured_minus_rec[1][2][0][0][0]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                theta_measured_minus_rec[1][2][0][0][0]->Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                                phi_measured_minus_rec[1][2][0][0][0]->Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                                //theta_rec_vs_mom[0][2][0][0][0]->Fill(_e->diff_pip_momentum(),_e->weight());

                                                        }
                                                        else if(_e->ftof_pim()) {
                                                                mom_measured_minus_rec[1][2][0][0][1]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                theta_measured_minus_rec[1][2][0][0][1]->Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                                phi_measured_minus_rec[1][2][0][0][1]->Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                                //theta_rec_vs_mom[0][2][0][0][1]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                        }
                                                }
                                                else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                                  mom_measured_minus_rec[1][2][0][1][0]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                                  theta_measured_minus_rec[1][2][0][1][0]->Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                                                  phi_measured_minus_rec[1][2][0][1][0]->Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                                                  //theta_rec_vs_mom[0][2][0][1][0]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                          }
                                                                          else if(_e->ftof_pim()) {
                                                                                  mom_measured_minus_rec[1][2][0][1][1]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                                  theta_measured_minus_rec[1][2][0][1][1]->Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                                                  phi_measured_minus_rec[1][2][0][1][1]->Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                                                  //theta_rec_vs_mom[0][2][0][1][1]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                          }}
                                        }
                                        else if(_e->ftof_prot()) {
                                                if(_e->ctof_pip()) {
                                                        if(_e->ctof_pim()) {
                                                                mom_measured_minus_rec[1][2][1][0][0]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                phi_measured_minus_rec[1][2][1][0][0]->Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                                //theta_rec_vs_mom[0][2][1][0][0]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                theta_measured_minus_rec[1][2][1][0][0]->Fill(_e->diff_pip_theta_lab(),_e->weight());

                                                        }
                                                        else if(_e->ftof_pim()) {
                                                                mom_measured_minus_rec[1][2][1][0][1]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                phi_measured_minus_rec[1][2][1][0][1]->Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                                //theta_rec_vs_mom[0][2][1][0][1]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                theta_measured_minus_rec[1][2][1][0][1]->Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                        }
                                                }
                                                else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                                  mom_measured_minus_rec[1][2][1][1][0]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                                  phi_measured_minus_rec[1][2][1][1][0]->Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                                                  //theta_rec_vs_mom[0][2][1][1][0]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                                  theta_measured_minus_rec[1][2][1][1][0]->Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                                          }
                                                                          else if(_e->ftof_pim()) {
                                                                                  mom_measured_minus_rec[1][2][1][1][1]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                                  phi_measured_minus_rec[1][2][1][1][1]->Fill(_e->diff_pip_Phi_lab(),_e->weight());
                                                                                  //theta_rec_vs_mom[0][2][1][1][0]->Fill(_e->diff_pip_momentum(),_e->weight());
                                                                                  theta_measured_minus_rec[1][2][1][1][1]->Fill(_e->diff_pip_theta_lab(),_e->weight());
                                                                          }}
                                        }
                                }
                        }
                }
        }
}
void Histogram::Fill_efficiency_check_1d_hist_pim(const std::shared_ptr<Reaction> &_e){
        if(_e->TwoPion_missingPim()) {
                if (_e->MM_cut()) {
                        if(_e->ctof_prot()) {
                                if(_e->ctof_pip()) {
                                        if(_e->ctof_pim()) {
                                                mom_measured_minus_rec[2][0][0][0][0] -> Fill(_e->diff_pim_momentum(),_e->weight());
                                                theta_measured_minus_rec[2][0][0][0][0] -> Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                phi_measured_minus_rec[2][0][0][0][0] -> Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                //theta_rec_vs_mom[0][0][0][0][0]->Fill(_e->diff_pim_momentum(),_e->weight());

                                        }
                                        else if(_e->ftof_pim()) {
                                                mom_measured_minus_rec[2][0][0][0][1] -> Fill(_e->diff_pim_momentum(),_e->weight());
                                                theta_measured_minus_rec[2][0][0][0][1] -> Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                phi_measured_minus_rec[2][0][0][0][1] -> Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                //theta_rec_vs_mom[0][0][0][0][1]->Fill(_e->diff_pim_momentum(),_e->weight());
                                        }
                                }
                                else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                  mom_measured_minus_rec[2][0][0][1][0] -> Fill(_e->diff_pim_momentum(),_e->weight());
                                                                  theta_measured_minus_rec[2][0][0][1][0] -> Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                                  phi_measured_minus_rec[2][0][0][1][0] -> Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                                  //theta_rec_vs_mom[0][0][0][1][0]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                          }
                                                          else if(_e->ftof_pim()) {
                                                                  mom_measured_minus_rec[2][0][0][1][1] -> Fill(_e->diff_pim_momentum(),_e->weight());
                                                                  theta_measured_minus_rec[2][0][0][1][1] -> Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                                  phi_measured_minus_rec[2][0][0][1][1] -> Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                                  //theta_rec_vs_mom[0][0][0][1][1]->Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                          }}
                        }
                        else if(_e->ftof_prot()) {
                                if(_e->ctof_pip()) {
                                        if(_e->ctof_pim()) {
                                                mom_measured_minus_rec[2][0][1][0][0] -> Fill(_e->diff_pim_momentum(),_e->weight());
                                                phi_measured_minus_rec[2][0][1][0][0] -> Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                //theta_rec_vs_mom[0][0][1][0][0]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                theta_measured_minus_rec[2][0][1][0][0] -> Fill(_e->diff_pim_theta_lab(),_e->weight());

                                        }
                                        else if(_e->ftof_pim()) {
                                                mom_measured_minus_rec[2][0][1][0][1] -> Fill(_e->diff_pim_momentum(),_e->weight());
                                                phi_measured_minus_rec[2][0][1][0][1] -> Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                //theta_rec_vs_mom[0][0][1][0][1]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                theta_measured_minus_rec[2][0][1][0][1] -> Fill(_e->diff_pim_theta_lab(),_e->weight());
                                        }
                                }
                                else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                  mom_measured_minus_rec[2][0][1][1][0] -> Fill(_e->diff_pim_momentum(),_e->weight());
                                                                  phi_measured_minus_rec[2][0][1][1][0] -> Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                                  //theta_rec_vs_mom[0][0][1][1][0]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                  theta_measured_minus_rec[2][0][1][1][0] -> Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                          }
                                                          else if(_e->ftof_pim()) {
                                                                  mom_measured_minus_rec[2][0][1][1][1] -> Fill(_e->diff_pim_momentum(),_e->weight());
                                                                  phi_measured_minus_rec[2][0][1][1][1] -> Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                                  //theta_rec_vs_mom[0][0][1][1][0]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                  theta_measured_minus_rec[2][0][1][1][1] -> Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                          }}
                        }

                        if (_e->TwoPion_exclusive())  {
                                if(_e->ctof_prot()) {
                                        if(_e->ctof_pip()) {
                                                if(_e->ctof_pim()) {
                                                        mom_measured_minus_rec[2][1][0][0][0] -> Fill(_e->diff_pim_momentum(),_e->weight());
                                                        theta_measured_minus_rec[2][1][0][0][0] -> Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                        phi_measured_minus_rec[2][1][0][0][0] -> Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                        //theta_rec_vs_mom[0][1][0][0][0]->Fill(_e->diff_pim_momentum(),_e->weight());

                                                }
                                                else if(_e->ftof_pim()) {
                                                        mom_measured_minus_rec[2][1][0][0][1] -> Fill(_e->diff_pim_momentum(),_e->weight());
                                                        theta_measured_minus_rec[2][1][0][0][1] -> Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                        phi_measured_minus_rec[2][1][0][0][1] -> Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                        //theta_rec_vs_mom[0][1][0][0][1]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                }
                                        }
                                        else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                          mom_measured_minus_rec[2][1][0][1][0] -> Fill(_e->diff_pim_momentum(),_e->weight());
                                                                          theta_measured_minus_rec[2][1][0][1][0] -> Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                                          phi_measured_minus_rec[2][1][0][1][0] -> Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                                          //theta_rec_vs_mom[0][1][0][1][0]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                  }
                                                                  else if(_e->ftof_pim()) {
                                                                          mom_measured_minus_rec[2][1][0][1][1] -> Fill(_e->diff_pim_momentum(),_e->weight());
                                                                          theta_measured_minus_rec[2][1][0][1][1] -> Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                                          phi_measured_minus_rec[2][1][0][1][1] -> Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                                          //theta_rec_vs_mom[0][1][0][1][1]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                  }}
                                }
                                else if(_e->ftof_prot()) {
                                        if(_e->ctof_pip()) {
                                                if(_e->ctof_pim()) {
                                                        mom_measured_minus_rec[2][1][1][0][0] -> Fill(_e->diff_pim_momentum(),_e->weight());
                                                        phi_measured_minus_rec[2][1][1][0][0] -> Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                        //theta_rec_vs_mom[0][1][1][0][0]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                        theta_measured_minus_rec[2][1][1][0][0] -> Fill(_e->diff_pim_theta_lab(),_e->weight());

                                                }
                                                else if(_e->ftof_pim()) {
                                                        mom_measured_minus_rec[2][1][1][0][1] -> Fill(_e->diff_pim_momentum(),_e->weight());
                                                        phi_measured_minus_rec[2][1][1][0][1] -> Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                        //theta_rec_vs_mom[0][1][1][0][1]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                        theta_measured_minus_rec[2][1][1][0][1] -> Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                }
                                        }
                                        else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                          mom_measured_minus_rec[2][1][1][1][0] -> Fill(_e->diff_pim_momentum(),_e->weight());
                                                                          phi_measured_minus_rec[2][1][1][1][0] -> Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                                          //theta_rec_vs_mom[0][1][1][1][0]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                          theta_measured_minus_rec[2][1][1][1][0] -> Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                                  }
                                                                  else if(_e->ftof_pim()) {
                                                                          mom_measured_minus_rec[2][1][1][1][1] -> Fill(_e->diff_pim_momentum(),_e->weight());
                                                                          phi_measured_minus_rec[2][1][1][1][1] -> Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                                          //theta_rec_vs_mom[0][1][1][1][0]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                          theta_measured_minus_rec[2][1][1][1][1] -> Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                                  }}
                                }

                                if(abs(_e->MM2_exclusive()) < 0.03 ) {
                                        if(_e->ctof_prot()) {
                                                if(_e->ctof_pip()) {
                                                        if(_e->ctof_pim()) {
                                                                mom_measured_minus_rec[2][2][0][0][0]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                theta_measured_minus_rec[2][2][0][0][0]->Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                                phi_measured_minus_rec[2][2][0][0][0]->Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                                //theta_rec_vs_mom[0][2][0][0][0]->Fill(_e->diff_pim_momentum(),_e->weight());

                                                        }
                                                        else if(_e->ftof_pim()) {
                                                                mom_measured_minus_rec[2][2][0][0][1]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                theta_measured_minus_rec[2][2][0][0][1]->Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                                phi_measured_minus_rec[2][2][0][0][1]->Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                                //theta_rec_vs_mom[0][2][0][0][1]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                        }
                                                }
                                                else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                                  mom_measured_minus_rec[2][2][0][1][0]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                                  theta_measured_minus_rec[2][2][0][1][0]->Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                                                  phi_measured_minus_rec[2][2][0][1][0]->Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                                                  //theta_rec_vs_mom[0][2][0][1][0]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                          }
                                                                          else if(_e->ftof_pim()) {
                                                                                  mom_measured_minus_rec[2][2][0][1][1]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                                  theta_measured_minus_rec[2][2][0][1][1]->Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                                                  phi_measured_minus_rec[2][2][0][1][1]->Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                                                  //theta_rec_vs_mom[0][2][0][1][1]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                          }}
                                        }
                                        else if(_e->ftof_prot()) {
                                                if(_e->ctof_pip()) {
                                                        if(_e->ctof_pim()) {
                                                                mom_measured_minus_rec[2][2][1][0][0]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                phi_measured_minus_rec[2][2][1][0][0]->Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                                //theta_rec_vs_mom[0][2][1][0][0]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                theta_measured_minus_rec[2][2][1][0][0]->Fill(_e->diff_pim_theta_lab(),_e->weight());

                                                        }
                                                        else if(_e->ftof_pim()) {
                                                                mom_measured_minus_rec[2][2][1][0][1]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                phi_measured_minus_rec[2][2][1][0][1]->Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                                //theta_rec_vs_mom[0][2][1][0][1]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                theta_measured_minus_rec[2][2][1][0][1]->Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                        }
                                                }
                                                else if(_e->ftof_pip()) { if(_e->ctof_pim()) {
                                                                                  mom_measured_minus_rec[2][2][1][1][0]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                                  phi_measured_minus_rec[2][2][1][1][0]->Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                                                  //theta_rec_vs_mom[0][2][1][1][0]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                                  theta_measured_minus_rec[2][2][1][1][0]->Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                                          }
                                                                          else if(_e->ftof_pim()) {
                                                                                  mom_measured_minus_rec[2][2][1][1][1]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                                  phi_measured_minus_rec[2][2][1][1][1]->Fill(_e->diff_pim_Phi_lab(),_e->weight());
                                                                                  //theta_rec_vs_mom[0][2][1][1][0]->Fill(_e->diff_pim_momentum(),_e->weight());
                                                                                  theta_measured_minus_rec[2][2][1][1][1]->Fill(_e->diff_pim_theta_lab(),_e->weight());
                                                                          }}
                                        }
                                }
                        }
                }
        }
}
// if (_e->TwoPion_missingPip()) {
//
//         if (_e->MM2_mpip()>-0.06 &&_e->MM2_mpip()<0.08) {
//                 //if(_e->pip_momentum_measured()!= NAN) {
//                 pip_theta_measured_minus_rec_CD[0]->Fill(_e->diff_pip_theta_lab(),_e->weight());
//                 pip_phi_measured_minus_rec_CD[0]->Fill(_e->diff_pip_Phi_lab(),_e->weight());
//                 pip_mom_measured_minus_rec_CD[0]->Fill(_e->diff_pip_momentum(),_e->weight());
//                 if (_e->TwoPion_exclusive())  {
//                         pip_theta_measured_minus_rec_CD[1]->Fill(_e->diff_pip_theta_lab(),_e->weight());
//                         pip_phi_measured_minus_rec_CD[1]->Fill(_e->diff_pip_Phi_lab(),_e->weight());
//                         pip_mom_measured_minus_rec_CD[1]->Fill(_e->diff_pip_momentum(),_e->weight());
//
//                         if(abs(_e->MM2_exclusive()) < 0.03 ) {
//                                 pip_theta_measured_minus_rec_CD[2]->Fill(_e->diff_pip_theta_lab(),_e->weight());
//                                 pip_phi_measured_minus_rec_CD[2]->Fill(_e->diff_pip_Phi_lab(),_e->weight());
//                                 pip_mom_measured_minus_rec_CD[2]->Fill(_e->diff_pip_momentum(),_e->weight());
//                         }
//                 }
//         }
//         //}
// }
//}
//void Histogram::Fill_efficiency_check_missingPip_FD(const std::shared_ptr<Reaction> &_e){
// if (_e->TwoPion_missingPip()) {
//
//         if (_e->MM2_mpip()>-0.06 &&_e->MM2_mpip()<0.08) {
//                 //if(_e->pip_momentum_measured()!= NAN) {
//                 pip_theta_measured_minus_rec_FD[0]->Fill(_e->diff_pip_theta_lab(),_e->weight());
//                 pip_phi_measured_minus_rec_FD[0]->Fill(_e->diff_pip_Phi_lab(),_e->weight());
//                 pip_mom_measured_minus_rec_FD[0]->Fill(_e->diff_pip_momentum(),_e->weight());
//                 if (_e->TwoPion_exclusive())  {
//                         pip_theta_measured_minus_rec_FD[1]->Fill(_e->diff_pip_theta_lab(),_e->weight());
//                         pip_phi_measured_minus_rec_FD[1]->Fill(_e->diff_pip_Phi_lab(),_e->weight());
//                         pip_mom_measured_minus_rec_FD[1]->Fill(_e->diff_pip_momentum(),_e->weight());
//
//                         if(abs(_e->MM2_exclusive()) < 0.03 ) {
//                                 pip_theta_measured_minus_rec_FD[2]->Fill(_e->diff_pip_theta_lab(),_e->weight());
//                                 pip_phi_measured_minus_rec_FD[2]->Fill(_e->diff_pip_Phi_lab(),_e->weight());
//                                 pip_mom_measured_minus_rec_FD[2]->Fill(_e->diff_pip_momentum(),_e->weight());
//                         }
//
//                 }
//         }
//         //}
// }
//}


//void Histogram::Fill_efficiency_check_missingProt_CD(const std::shared_ptr<Reaction> &_e){
// if (_e->TwoPion_missingProt()) {
//         if (_e->TwoPion_missingProt())
//
//                 if (_e->MM2_mprot()>0.7 &&_e->MM2_mprot()<1.2) {
//                         //if(_e->prot_momentum_measured()!= NAN) {
//                         prot_theta_measured_minus_rec_CD[0]->Fill(_e->diff_prot_theta_lab(),_e->weight());
//                         prot_phi_measured_minus_rec_CD[0]->Fill(_e->diff_prot_Phi_lab(),_e->weight());
//                         prot_mom_measured_minus_rec_CD[0]->Fill(_e->diff_prot_momentum(),_e->weight());
//                         if (_e->TwoPion_exclusive())  {
//                                 prot_theta_measured_minus_rec_CD[1]->Fill(_e->diff_prot_theta_lab(),_e->weight());
//                                 prot_phi_measured_minus_rec_CD[1]->Fill(_e->diff_prot_Phi_lab(),_e->weight());
//                                 prot_mom_measured_minus_rec_CD[1]->Fill(_e->diff_prot_momentum(),_e->weight());
//
//                                 if(abs(_e->MM2_exclusive()) < 0.03 ) {
//                                         prot_theta_measured_minus_rec_CD[2]->Fill(_e->diff_prot_theta_lab(),_e->weight());
//                                         prot_phi_measured_minus_rec_CD[2]->Fill(_e->diff_prot_Phi_lab(),_e->weight());
//                                         prot_mom_measured_minus_rec_CD[2]->Fill(_e->diff_prot_momentum(),_e->weight());
//                                 }
//                         }
//                 }
//         //}
// }
//}
//void Histogram::Fill_efficiency_check_missingProt_FD(const std::shared_ptr<Reaction> &_e){
// if (_e->TwoPion_missingProt()) {
//         if (_e->MM2_mprot()>0.7 &&_e->MM2_mprot()<1.2) {
//                 //if(_e->prot_momentum_measured()!= NAN) {
//                 prot_theta_measured_minus_rec_FD[0]->Fill(_e->diff_prot_theta_lab(),_e->weight());
//                 prot_phi_measured_minus_rec_FD[0]->Fill(_e->diff_prot_Phi_lab(),_e->weight());
//                 prot_mom_measured_minus_rec_FD[0]->Fill(_e->diff_prot_momentum(),_e->weight());
//                 if (_e->TwoPion_exclusive())  {
//                         prot_theta_measured_minus_rec_FD[1]->Fill(_e->diff_prot_theta_lab(),_e->weight());
//                         prot_phi_measured_minus_rec_FD[1]->Fill(_e->diff_prot_Phi_lab(),_e->weight());
//                         prot_mom_measured_minus_rec_FD[1]->Fill(_e->diff_prot_momentum(),_e->weight());
//
//                         if(abs(_e->MM2_exclusive()) < 0.03 ) {
//                                 prot_theta_measured_minus_rec_FD[2]->Fill(_e->diff_prot_theta_lab(),_e->weight());
//                                 prot_phi_measured_minus_rec_FD[2]->Fill(_e->diff_prot_Phi_lab(),_e->weight());
//                                 prot_mom_measured_minus_rec_FD[2]->Fill(_e->diff_prot_momentum(),_e->weight());
//                         }
//
//                 }
//         }
//         //}
// }
//}

void Histogram::Fill_efficiency_check_CD(const std::shared_ptr<Reaction> &_e){
        if (_e->TwoPion_missingPim()) {
                if (_e->MM_cut()) {
                        if(_e->pim_momentum_measured()!= NAN) {

                                pim_theta_measured_vs_mom_CD[0]->Fill(_e->pim_momentum_measured(),_e->pim_theta_lab_measured(),_e->weight());
                                pim_phi_measured_vs_mom_CD[0]->Fill(_e->pim_momentum_measured(),_e->pim_Phi_lab_measured(),_e->weight());
                                pim_E_vs_mom_measured_CD[0]->Fill(_e->pim_momentum_measured(),_e->pim_E_measured(),_e->weight());
                                pim_mom_measured_vs_rec_CD[0]->Fill(_e->pim_momentum(),_e->pim_momentum_measured(),_e->weight());
                                pim_E_measured_vs_rec_CD[0]->Fill(_e->pim_E(),_e->pim_E_measured(),_e->weight());
                                pim_theta_measured_vs_rec_CD[0]->Fill(_e->pim_theta_lab(),_e->pim_theta_lab_measured(),_e->weight());
                                pim_phi_measured_vs_rec_CD[0]->Fill(_e->pim_Phi_lab(),_e->pim_Phi_lab_measured(),_e->weight());
                                pim_phi_vs_theta_measured_CD[0]->Fill(_e->pim_theta_lab_measured(),_e->pim_Phi_lab_measured(),_e->weight());
                                //if(_e->pim_momentum_measured()>0)
                                pim_theta_measured_minus_rec_vs_measured_CD[0]->Fill(_e->pim_theta_lab_measured(), (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                //        pim_theta_measured_minus_rec_vs_rec_CD[0]->Fill(_e->pim_theta_lab(), (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                // pim_theta_measured_minus_rec_CD[0]->Fill( (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                // pim_phi_measured_minus_rec_CD[0]->Fill( (_e->pim_Phi_lab_measured()-_e->pim_Phi_lab()),_e->weight());
                                // pim_mom_measured_minus_rec_CD[0]->Fill( (_e->pim_momentum_measured()-_e->pim_momentum()),_e->weight());
                                // pim_theta_measured_minus_rec_vs_rec_CD[0]->Fill(_e->pim_theta_lab(), (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                // //pim_theta_measured_vs_thrown_CD[0]->Fill(_e->pim_theta_thrown(), _e->pim_theta_lab_measured());
                                // //pim_theta_rec_vs_thrown_CD[0]->Fill(_e->pim_theta_thrown(),_e->pim_theta_lab());
                                //
                                // //if(_e->pim_momentum_measured()>0) {
                                if (_e->TwoPion_exclusive())  {
                                        pim_theta_measured_vs_mom_CD[1]->Fill(_e->pim_momentum_measured(),_e->pim_theta_lab_measured(),_e->weight());
                                        pim_phi_measured_vs_mom_CD[1]->Fill(_e->pim_momentum_measured(),_e->pim_Phi_lab_measured(),_e->weight());
                                        pim_E_vs_mom_measured_CD[1]->Fill(_e->pim_momentum_measured(),_e->pim_E_measured(),_e->weight());
                                        pim_mom_measured_vs_rec_CD[1]->Fill(_e->pim_momentum(),_e->pim_momentum_measured(),_e->weight());
                                        pim_E_measured_vs_rec_CD[1]->Fill(_e->pim_E(),_e->pim_E_measured(),_e->weight());
                                        pim_theta_measured_vs_rec_CD[1]->Fill(_e->pim_theta_lab(),_e->pim_theta_lab_measured(),_e->weight());
                                        pim_phi_measured_vs_rec_CD[1]->Fill(_e->pim_Phi_lab(),_e->pim_Phi_lab_measured(),_e->weight());
                                        pim_phi_vs_theta_measured_CD[1]->Fill(_e->pim_theta_lab_measured(),_e->pim_Phi_lab_measured(),_e->weight());

                                        pim_theta_measured_minus_rec_vs_measured_CD[1]->Fill(_e->pim_theta_lab_measured(),( _e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                        //        pim_theta_measured_minus_rec_vs_rec_CD[1]->Fill(_e->pim_theta_lab(), (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                        // pim_theta_measured_minus_rec_CD[1]->Fill( (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                        // pim_phi_measured_minus_rec_CD[1]->Fill( (_e->pim_Phi_lab_measured()-_e->pim_Phi_lab()),_e->weight());
                                        // pim_mom_measured_minus_rec_CD[1]->Fill( (_e->pim_momentum_measured()-_e->pim_momentum()),_e->weight());

                                        // pim_theta_measured_minus_rec_vs_rec_CD[1]->Fill(_e->pim_theta_lab(),( _e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                        // //pim_theta_measured_vs_thrown_CD[1]->Fill(_e->pim_theta_thrown(), _e->pim_theta_lab_measured());
                                        // //pim_theta_rec_vs_thrown_CD[1]->Fill(_e->pim_theta_thrown(),_e->pim_theta_lab());
                                        if(abs(_e->MM2_exclusive()) < 0.03 ) {
                                                MM2_exclusive_hist_CD[2]->Fill(_e->MM2(),_e->weight());
                                                pim_theta_measured_vs_mom_CD[2]->Fill(_e->pim_momentum_measured(),_e->pim_theta_lab_measured(),_e->weight());
                                                pim_phi_measured_vs_mom_CD[2]->Fill(_e->pim_momentum_measured(),_e->pim_Phi_lab_measured(),_e->weight());
                                                pim_E_vs_mom_measured_CD[2]->Fill(_e->pim_momentum_measured(),_e->pim_E_measured(),_e->weight());
                                                pim_theta_measured_vs_rec_CD[2]->Fill(_e->pim_theta_lab(),_e->pim_theta_lab_measured(),_e->weight());
                                                pim_theta_measured_minus_rec_vs_measured_CD[2]->Fill(_e->pim_theta_lab_measured(),( _e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                                //pim_theta_measured_minus_rec_vs_rec_CD[2]->Fill(_e->pim_theta_lab(), (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                                pim_phi_measured_vs_rec_CD[2]->Fill(_e->pim_Phi_lab(),_e->pim_Phi_lab_measured(),_e->weight());
                                                pim_E_measured_vs_rec_CD[2]->Fill(_e->pim_E(),_e->pim_E_measured(),_e->weight());
                                                pim_mom_measured_vs_rec_CD[2]->Fill(_e->pim_momentum(),_e->pim_momentum_measured(),_e->weight());
                                                // pim_theta_measured_minus_rec_CD[2]->Fill( (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                                // pim_phi_measured_minus_rec_CD[2]->Fill( (_e->pim_Phi_lab_measured()-_e->pim_Phi_lab()),_e->weight());
                                                // pim_mom_measured_minus_rec_CD[2]->Fill( (_e->pim_momentum_measured()-_e->pim_momentum()),_e->weight());


                                        }
                                }

                        }
                }
        }

}

void Histogram::Fill_efficiency_check_FD_rec(const std::shared_ptr<Reaction> &_e){
        short sec = _e->pim_sec();
        if (_e->TwoPion_missingPim()) {
                if (_e->MM_cut()) {
                        MM2_exclusive_hist_FD[0]->Fill(_e->MM2(),_e->weight());

                        if(_e->pim_momentum() != NAN) {

                                if (sec > 0 && sec <= 6) {
                                        theta_vs_mom_elec[sec - 1]->Fill(_e->elec_momentum(),_e->theta_elec(), _e->weight());
                                        theta_vs_mom_prot[sec - 1]->Fill(_e->prot_momentum(),_e->prot_theta_lab(), _e->weight());
                                        theta_vs_mom_pip[sec - 1]->Fill(_e->pip_momentum(),_e->pip_theta_lab(), _e->weight());
                                        theta_vs_mom_pim[sec - 1]->Fill(_e->pim_momentum(),_e->pim_theta_lab(), _e->weight());


                                        pim_phi_vs_theta_rec_FD_sec[sec - 1]->Fill(_e->pim_theta_lab(),_e->pim_Phi_lab(),_e->weight());
                                        //if((sec -1) == 4) {
                                        //        float phi_value_for_sec_4 = _e->pim_Phi_lab();
                                        //        if (phi_value_for_sec_4 > 180) {phi_value_for_sec_4 = phi_value_for_sec_4 -360;}
                                        //pim_phi_vs_theta_rec_FD_sec[sec - 1]->Fill(_e->pim_theta_lab(),/* phi_value_for_sec_4*/ _e->pim_Phi_lab(),_e->weight());
                                        //}

                                        if (_e->TwoPion_exclusive())  {
                                                MM2_exclusive_hist_FD[1]->Fill(_e->MM2(),_e->weight());

                                                pim_phi_vs_theta_rec_FD_after_exclusive_sec[sec - 1]->Fill(_e->pim_theta_lab(),_e->pim_Phi_lab(),_e->weight());
                                                //if((sec -1) == 4) {
                                                //        float phi_value_for_sec_4_ = _e->pim_Phi_lab();
                                                //        if (phi_value_for_sec_4_ > 180) {phi_value_for_sec_4_ = phi_value_for_sec_4_ -360;}
                                                //pim_phi_vs_theta_rec_FD_after_exclusive_sec[sec - 1]->Fill(_e->pim_theta_lab(), /*phi_value_for_sec_4_ */ _e->pim_Phi_lab(),_e->weight());
                                                //}
                                        }
                                }

                                pim_theta_rec_vs_mom_FD[0]->Fill(_e->pim_momentum(),_e->pim_theta_lab(),_e->weight());
                                pim_phi_rec_vs_mom_FD[0]->Fill(_e->pim_momentum(),_e->pim_Phi_lab(),_e->weight());
                                pim_phi_vs_theta_rec_FD[0]->Fill(_e->pim_theta_lab(),_e->pim_Phi_lab(),_e->weight());
                                pim_E_vs_mom_rec_FD[0]->Fill(_e->pim_momentum(),_e->pim_E(),_e->weight());

                                if (_e->TwoPion_exclusive())  {
                                        pim_theta_rec_vs_mom_FD[1]->Fill(_e->pim_momentum(),_e->pim_theta_lab(),_e->weight());
                                        pim_phi_rec_vs_mom_FD[1]->Fill(_e->pim_momentum(),_e->pim_Phi_lab(),_e->weight());
                                        pim_phi_vs_theta_rec_FD[1]->Fill(_e->pim_theta_lab(),_e->pim_Phi_lab(),_e->weight());
                                        pim_E_vs_mom_rec_FD[1]->Fill(_e->pim_momentum(),_e->pim_E(),_e->weight());
                                }

                        }
                }
        }
}
// void Histogram::Fill_efficiency_check_FD_thrown(const std::shared_ptr<MCReaction> &_mc_e,const std::shared_ptr<Reaction> &_e) {
// //        short sec = _e->sec();
//         if (_e->TwoPion_missingPim()) {
//                 if(_e->pim_momentum_measured()!= NAN) {
//                         pim_theta_measured_vs_thrown_FD[0]->Fill(_mc_e->MCpim_theta_lab(), (_e->pim_theta_lab_measured()-_mc_e->MCpim_theta_lab()),_e->weight());
//                         if (_e->MM_cut()) {
//                                 pim_theta_measured_vs_thrown_FD[1]->Fill(_mc_e->MCpim_theta_lab(), (_e->pim_theta_lab_measured()-_mc_e->MCpim_theta_lab()),_e->weight());
//                         }
//                         if(abs(_mc_e->MCpim_theta_lab() -_e->pim_theta_lab_measured()) < 10 )
//                                 pim_theta_measured_vs_thrown_FD[2]->Fill(_mc_e->MCpim_theta_lab(), (_e->pim_theta_lab_measured()-_mc_e->MCpim_theta_lab()),_e->weight());
//                 }
//         }
// }
void Histogram::Fill_efficiency_check_FD(const std::shared_ptr<Reaction> &_e) {
        if (_e->TwoPion_missingPim()) {
                if (_e->MM_cut()) {
                        short sec = _e->pim_sec();
                        if (sec > 0 && sec <= 6) {
                                if(_e->pim_momentum_measured()> 0.0) {
                                        pim_phi_vs_theta_measured_FD_sec[sec - 1]->Fill(_e->pim_theta_lab_measured(),_e->pim_Phi_lab_measured(),_e->weight());
                                        //if((sec -1) == 4) {
                                        //        float phi_value_for_sec_4_measured = _e->pim_Phi_lab();
                                        //        if (phi_value_for_sec_4_measured > 180) {phi_value_for_sec_4_measured = phi_value_for_sec_4_measured -360;}
                                        //pim_phi_vs_theta_measured_FD_sec[sec - 1]->Fill(_e->pim_theta_lab_measured(),/* phi_value_for_sec_4_measured */ _e->pim_Phi_lab_measured(),_e->weight());
                                        //}
                                        // mmsq > -0.06 && mmsq < 0.08)

                                        pim_phi_vs_theta_measured_FD_after_exclusive_sec[sec - 1]->Fill(_e->pim_theta_lab_measured(),_e->pim_Phi_lab_measured(),_e->weight());
                                        //if((sec -1) == 4) {
                                        //        float phi_value_for_sec_4_measured_  = _e->pim_Phi_lab();
                                        //        if (phi_value_for_sec_4_measured_  > 180) {phi_value_for_sec_4_measured_ = phi_value_for_sec_4_measured_  -360;}
                                        //pim_phi_vs_theta_measured_FD_after_exclusive_sec[sec - 1]->Fill(_e->pim_theta_lab_measured(),/* phi_value_for_sec_4_measured_ */ _e->pim_Phi_lab_measured(),_e->weight());
                                        //}
                                }
                        }
                        if(_e->pim_momentum_measured()!= NAN) {
                                pim_theta_measured_vs_mom_FD[0]->Fill(_e->pim_momentum_measured(),_e->pim_theta_lab_measured(),_e->weight());
                                pim_phi_measured_vs_mom_FD[0]->Fill(_e->pim_momentum_measured(),_e->pim_Phi_lab_measured(),_e->weight());

                                pim_E_vs_mom_measured_FD[0]->Fill(_e->pim_momentum_measured(),_e->pim_E_measured(),_e->weight());

                                pim_mom_measured_vs_rec_FD[0]->Fill(_e->pim_momentum(),_e->pim_momentum_measured(),_e->weight());
                                pim_E_measured_vs_rec_FD[0]->Fill(_e->pim_E(),_e->pim_E_measured(),_e->weight());


                                pim_theta_measured_vs_rec_FD[0]->Fill(_e->pim_theta_lab(),_e->pim_theta_lab_measured(),_e->weight());
                                pim_phi_measured_vs_rec_FD[0]->Fill(_e->pim_Phi_lab(),_e->pim_Phi_lab_measured(),_e->weight());
                                pim_phi_vs_theta_measured_FD[0]->Fill(_e->pim_theta_lab_measured(),_e->pim_Phi_lab_measured(),_e->weight());
                                pim_theta_measured_minus_rec_vs_measured_FD[0]->Fill(_e->pim_theta_lab_measured(), (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                //pim_theta_measured_minus_rec_vs_rec_FD[0]->Fill(_e->pim_theta_lab(), (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                // pim_theta_measured_minus_rec_FD[0]->Fill( (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                // pim_phi_measured_minus_rec_FD[0]->Fill( (_e->pim_Phi_lab_measured()-_e->pim_Phi_lab()),_e->weight());
                                // pim_mom_measured_minus_rec_FD[0]->Fill( (_e->pim_momentum_measured()-_e->pim_momentum()),_e->weight());

                                // //if(_e->pim_momentum_measured()>0) {
                                if (_e->TwoPion_exclusive())  {

                                        pim_theta_measured_vs_mom_FD[1]->Fill(_e->pim_momentum_measured(),_e->pim_theta_lab_measured(),_e->weight());
                                        pim_phi_measured_vs_mom_FD[1]->Fill(_e->pim_momentum_measured(),_e->pim_Phi_lab_measured(),_e->weight());
                                        pim_E_vs_mom_measured_FD[1]->Fill(_e->pim_momentum_measured(),_e->pim_E_measured(),_e->weight());

                                        pim_mom_measured_vs_rec_FD[1]->Fill(_e->pim_momentum(),_e->pim_momentum_measured(),_e->weight());
                                        pim_E_measured_vs_rec_FD[1]->Fill(_e->pim_E(),_e->pim_E_measured(),_e->weight());

                                        pim_theta_measured_vs_rec_FD[1]->Fill(_e->pim_theta_lab(),_e->pim_theta_lab_measured(),_e->weight());
                                        pim_phi_measured_vs_rec_FD[1]->Fill(_e->pim_Phi_lab(),_e->pim_Phi_lab_measured(),_e->weight());
                                        pim_phi_vs_theta_measured_FD[1]->Fill(_e->pim_theta_lab_measured(),_e->pim_Phi_lab_measured(),_e->weight());
                                        pim_theta_measured_minus_rec_vs_measured_FD[1]->Fill(_e->pim_theta_lab_measured(), (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                        //pim_theta_measured_minus_rec_vs_rec_FD[1]->Fill(_e->pim_theta_lab(), (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                        // pim_theta_measured_minus_rec_FD[1]->Fill( (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                        // pim_phi_measured_minus_rec_FD[1]->Fill( (_e->pim_Phi_lab_measured()-_e->pim_Phi_lab()),_e->weight());
                                        // pim_mom_measured_minus_rec_FD[1]->Fill( (_e->pim_momentum_measured()-_e->pim_momentum()),_e->weight());

                                        if(abs(_e->MM2_exclusive()) < 0.03 ) {
                                                MM2_exclusive_hist_FD[2]->Fill(_e->MM2(),_e->weight());
                                                pim_theta_measured_vs_mom_FD[2]->Fill(_e->pim_momentum_measured(),_e->pim_theta_lab_measured(),_e->weight());
                                                pim_phi_measured_vs_mom_FD[2]->Fill(_e->pim_momentum_measured(),_e->pim_Phi_lab_measured(),_e->weight());
                                                pim_E_vs_mom_measured_FD[2]->Fill(_e->pim_momentum_measured(),_e->pim_E_measured(),_e->weight());
                                                pim_theta_measured_vs_rec_FD[2]->Fill(_e->pim_theta_lab(),_e->pim_theta_lab_measured(),_e->weight());
                                                pim_theta_measured_minus_rec_vs_measured_FD[2]->Fill(_e->pim_theta_lab_measured(), (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                                //pim_theta_measured_minus_rec_vs_rec_FD[2]->Fill(_e->pim_theta_lab(), (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                                // pim_theta_measured_minus_rec_FD[2]->Fill( (_e->pim_theta_lab_measured()-_e->pim_theta_lab()),_e->weight());
                                                // pim_phi_measured_minus_rec_FD[2]->Fill( (_e->pim_Phi_lab_measured()-_e->pim_Phi_lab()),_e->weight());
                                                // pim_mom_measured_minus_rec_FD[2]->Fill( (_e->pim_momentum_measured()-_e->pim_momentum()),_e->weight());
                                                pim_phi_measured_vs_rec_FD[2]->Fill(_e->pim_Phi_lab(),_e->pim_Phi_lab_measured(),_e->weight());
                                                pim_E_measured_vs_rec_FD[2]->Fill(_e->pim_E(),_e->pim_E_measured(),_e->weight());
                                                pim_mom_measured_vs_rec_FD[2]->Fill(_e->pim_momentum(),_e->pim_momentum_measured(),_e->weight());

                                        }
                                }

                        }

                }

        }
}
void Histogram::write_hist_efficiency_check_1d() {
        for (size_t h = 0; h < HADRON_NUM; h++) {
                for (size_t e = 0; e < EFF_CONDITIONS_NUM_ALL; e++) {
                        for (size_t pt = 0; pt < DETECTOR_NUM_PROT; pt++) {
                                for (size_t pp = 0; pp < DETECTOR_NUM_PIP; pp++) {
                                        for (size_t pm = 0; pm < DETECTOR_NUM_PIM; pm++) {
                                                mom_measured_minus_rec[h][e][pt][pp][pm]->SetXTitle("P (GeV)");
                                                mom_measured_minus_rec[h][e][pt][pp][pm]->Write();

                                                theta_measured_minus_rec[h][e][pt][pp][pm]->SetXTitle("#theta (deg)");
                                                theta_measured_minus_rec[h][e][pt][pp][pm]->Write();

                                                phi_measured_minus_rec[h][e][pt][pp][pm]->SetXTitle("#phi (deg)");
                                                phi_measured_minus_rec[h][e][pt][pp][pm]->Write();

                                                //theta_rec_vs_mom[h][e][pt][pp][pm]->SetTitle("#theta vs momentum");
                                                //theta_rec_vs_mom[h][e][pt][pp][pm]->SetXTitle("P (GeV)");
                                                //theta_rec_vs_mom[h][e][pt][pp][pm]->SetYTitle("#theta (deg)");
                                                //theta_rec_vs_mom[h][e][pt][pp][pm]->Write();
                                        }
                                }
                        }
                }
        }
}



// for (size_t i = 0; i < EFF_CONDITIONS_NUM; i++) {
//
//         // pip_theta_measured_minus_rec_CD[i]->Fit("gaus", "QMR+", "QMR+", -6, 5);
//         // gStyle->SetOptFit(1111);
//         pip_theta_measured_minus_rec_CD[i]->SetXTitle("#theta (deg)");
//         pip_theta_measured_minus_rec_CD[i]->Write();
//
//         // pip_phi_measured_minus_rec_FD[i]->Fit("gaus", "QMR+", "QMR+", -2, 2);
//         pip_phi_measured_minus_rec_CD[i]->SetXTitle("#phi (deg)");
//         pip_phi_measured_minus_rec_CD[i]->Write();
//
//         // pip_mom_measured_minus_rec_FD[i]->Fit("gaus", "QMR+", "QMR+", -2, 2);
//         pip_mom_measured_minus_rec_CD[i]->SetXTitle("P (GeV)");
//         pip_mom_measured_minus_rec_CD[i]->Write();
//
//
//         pip_theta_measured_minus_rec_FD[i]->SetXTitle("#theta (deg)");
//         pip_theta_measured_minus_rec_FD[i]->Write();
//
//         // pip_phi_measured_minus_rec_FD[i]->Fit("gaus", "QMR+", "QMR+", -2, 2);
//         pip_phi_measured_minus_rec_FD[i]->SetXTitle("#phi (deg)");
//         pip_phi_measured_minus_rec_FD[i]->Write();
//
//         // pip_mom_measured_minus_rec_FD[i]->Fit("gaus", "QMR+", "QMR+", -2, 2);
//         pip_mom_measured_minus_rec_FD[i]->SetXTitle("P (GeV)");
//         pip_mom_measured_minus_rec_FD[i]->Write();
//
//
//
//         // prot_theta_measured_minus_rec_CD[i]->Fit("gaus", "QMR+", "QMR+", -6, 5);
//         prot_theta_measured_minus_rec_CD[i]->SetXTitle("#theta (deg)");
//         prot_theta_measured_minus_rec_CD[i]->Write();
//
//         // prot_phi_measured_minus_rec_FD[i]->Fit("gaus", "QMR+", "QMR+", -2, 2);
//         prot_phi_measured_minus_rec_CD[i]->SetXTitle("#phi (deg)");
//         prot_phi_measured_minus_rec_CD[i]->Write();
//
//         // prot_mom_measured_minus_rec_FD[i]->Fit("gaus", "QMR+", "QMR+", -2, 2);
//         prot_mom_measured_minus_rec_CD[i]->SetXTitle("P (GeV)");
//         prot_mom_measured_minus_rec_CD[i]->Write();
//
//
//         prot_theta_measured_minus_rec_FD[i]->SetXTitle("#theta (deg)");
//         prot_theta_measured_minus_rec_FD[i]->Write();
//
//         // prot_phi_measured_minus_rec_FD[i]->Fit("gaus", "QMR+", "QMR+", -2, 2);
//         prot_phi_measured_minus_rec_FD[i]->SetXTitle("#phi (deg)");
//         prot_phi_measured_minus_rec_FD[i]->Write();
//
//         // prot_mom_measured_minus_rec_FD[i]->Fit("gaus", "QMR+", "QMR+", -2, 2);
//         prot_mom_measured_minus_rec_FD[i]->SetXTitle("P (GeV)");
//         prot_mom_measured_minus_rec_FD[i]->Write();

//         }
// }
void Histogram::write_hist_efficiency_check_CD() {

        for (size_t i = 0; i < EFF_CONDITIONS_NUM; i++) {

                pip_theta_rec_vs_mom[i]->SetOption("COLZ");
                pip_theta_rec_vs_mom[i]->Write();
                pip_theta_measured_vs_mom[i]->SetOption("COLZ");
                pip_theta_measured_vs_mom[i]->Write();

                prot_theta_rec_vs_mom[i]->SetOption("COLZ");
                prot_theta_rec_vs_mom[i]->Write();
                prot_theta_measured_vs_mom[i]->SetOption("COLZ");
                prot_theta_measured_vs_mom[i]->Write();

                MM2_exclusive_hist_CD[i]->SetXTitle("MMSQ (GeV)");
                MM2_exclusive_hist_CD[i]->Write();
                // pim_theta_measured_minus_rec_CD[i]->Fit("gaus", "QMR+", "QMR+", -6, 5);
                // gStyle->SetOptFit(1111);
                // pim_theta_measured_minus_rec_CD[i]->SetXTitle("#theta (deg)");
                // pim_theta_measured_minus_rec_CD[i]->Write();
                //
                // // pim_phi_measured_minus_rec_FD[i]->Fit("gaus", "QMR+", "QMR+", -2, 2);
                // pim_phi_measured_minus_rec_CD[i]->SetXTitle("#phi (deg)");
                // pim_phi_measured_minus_rec_CD[i]->Write();
                //
                // // pim_mom_measured_minus_rec_FD[i]->Fit("gaus", "QMR+", "QMR+", -2, 2);
                // pim_mom_measured_minus_rec_CD[i]->SetXTitle("P (GeV)");
                // pim_mom_measured_minus_rec_CD[i]->Write();


                pim_theta_rec_vs_mom_CD[i]->SetOption("COLZ");
                pim_theta_rec_vs_mom_CD[i]->Write();

                pim_phi_rec_vs_mom_CD[i]->SetOption("COLZ");
                pim_phi_rec_vs_mom_CD[i]->Write();


                pim_theta_measured_vs_mom_CD[i]->SetOption("COLZ");
                pim_theta_measured_vs_mom_CD[i]->Write();
                pim_phi_measured_vs_mom_CD[i]->SetOption("COLZ");
                pim_phi_measured_vs_mom_CD[i]->Write();

                pim_E_vs_mom_rec_CD[i]->SetOption("COLZ");
                pim_E_vs_mom_rec_CD[i]->Write();

                pim_E_vs_mom_measured_CD[i]->SetOption("COLZ");
                pim_E_vs_mom_measured_CD[i]->Write();

                pim_phi_vs_theta_rec_CD[i]->SetOption("COLZ");
                pim_phi_vs_theta_rec_CD[i]->SetTitle("#phi_{#pi^{-}} vs #theta_{#pi^{-}} measured ");
                pim_phi_vs_theta_rec_CD[i]->SetXTitle("#theta_{#pi^{-}} (deg)");
                pim_phi_vs_theta_rec_CD[i]->SetYTitle("#phi_{#pi^{-}} (deg)");
                pim_phi_vs_theta_rec_CD[i]->Write();


                pim_E_measured_vs_rec_CD[i]->SetOption("COLZ");
                pim_E_measured_vs_rec_CD[i]->SetTitle("E_{#pi^{-}} measured vs E_{#pi^{-}} rec ");
                pim_E_measured_vs_rec_CD[i]->SetXTitle("E_{#pi^{-}} rec (GeV)");
                pim_E_measured_vs_rec_CD[i]->SetYTitle("E_{#pi^{-}} measured (GeV)");
                if (pim_E_measured_vs_rec_CD[i]->GetEntries())
                        pim_E_measured_vs_rec_CD[i]->Write();


                pim_mom_measured_vs_rec_CD[i]->SetOption("COLZ");
                pim_mom_measured_vs_rec_CD[i]->SetTitle("P_{#pi^{-}} measured vs P_{#pi^{-}} rec ");
                pim_mom_measured_vs_rec_CD[i]->SetXTitle("P_{#pi^{-}} rec (GeV)");
                pim_mom_measured_vs_rec_CD[i]->SetYTitle("P_{#pi^{-}} measured (GeV)");
                if (pim_mom_measured_vs_rec_CD[i]->GetEntries())
                        pim_mom_measured_vs_rec_CD[i]->Write();

                pim_theta_measured_vs_rec_CD[i]->SetOption("COLZ");
                pim_theta_measured_vs_rec_CD[i]->SetTitle("#theta_{#pi^{-}} measured vs #theta_{#pi^{-}} rec ");
                pim_theta_measured_vs_rec_CD[i]->SetXTitle("#theta_{#pi^{-}} rec (deg)");
                pim_theta_measured_vs_rec_CD[i]->SetYTitle("#theta_{#pi^{-}} measured (deg)");
                if (pim_theta_measured_vs_rec_CD[i]->GetEntries())
                        pim_theta_measured_vs_rec_CD[i]->Write();

                pim_phi_measured_vs_rec_CD[i]->SetOption("COLZ");
                pim_phi_measured_vs_rec_CD[i]->SetTitle("#phi_{#pi^{-}} measured vs #phi_{#pi^{-}} rec ");
                pim_phi_measured_vs_rec_CD[i]->SetXTitle("#phi_{#pi^{-}} rec (deg)");
                pim_phi_measured_vs_rec_CD[i]->SetYTitle("#phi_{#pi^{-}} measured (deg)");
                if (pim_phi_measured_vs_rec_CD[i]->GetEntries())
                        pim_phi_measured_vs_rec_CD[i]->Write();

                pim_phi_vs_theta_measured_CD[i]->SetOption("COLZ");
                pim_phi_vs_theta_measured_CD[i]->SetOption("COLZ");
                pim_phi_vs_theta_measured_CD[i]->SetTitle("#phi_{#pi^{-}} vs #theta_{#pi^{-}} measured ");
                pim_phi_vs_theta_measured_CD[i]->SetXTitle("#theta_{#pi^{-}} (deg)");
                pim_phi_vs_theta_measured_CD[i]->SetYTitle("#phi_{#pi^{-}} (deg)");
                pim_phi_vs_theta_measured_CD[i]->Write();


                pim_theta_measured_minus_rec_vs_measured_CD[i]->SetOption("COLZ");
                pim_theta_measured_minus_rec_vs_measured_CD[i]->SetTitle("#theta_{#pi^{-}} measured - #theta_{#pi^{-}} rec vs #theta_{#pi^{-}} rec ");
                pim_theta_measured_minus_rec_vs_measured_CD[i]->SetXTitle("#theta_{#pi^{-}} rec  (deg)");
                pim_theta_measured_minus_rec_vs_measured_CD[i]->SetYTitle("#theta_{#pi^{-}} measured - #theta_{#pi^{-}} rec (deg)");
                //if (pim_theta_measured_minus_rec_vs_measured_CD[i]->GetEntries())
                pim_theta_measured_minus_rec_vs_measured_CD[i]->Write();

                // pim_theta_measured_minus_rec_vs_rec_CD[i]->SetOption("COLZ");
                // pim_theta_measured_minus_rec_vs_rec_CD[i]->SetTitle("#theta_{#pi^{-}} measured - #theta_{#pi^{-}} rec vs #theta_{#pi^{-}} rec ");
                // pim_theta_measured_minus_rec_vs_rec_CD[i]->SetXTitle("#theta_{#pi^{-}} measured  (deg)");
                // pim_theta_measured_minus_rec_vs_rec_CD[i]->SetYTitle("#theta_{#pi^{-}} measured - #theta_{#pi^{-}} rec (deg)");
                // //        if (pim_theta_measured_minus_rec_vs_rec_CD[i]->GetEntries())
                // pim_theta_measured_minus_rec_vs_rec_CD[i]->Write();


                pim_theta_measured_vs_thrown_CD[i]->SetOption("COLZ");
                pim_theta_measured_vs_thrown_CD[i]->SetTitle("#theta_{#pi^{-}} measured - #theta_{#pi^{-}} thrown vs #theta_{#pi^{-}} thrown ");
                pim_theta_measured_vs_thrown_CD[i]->SetXTitle("#theta_{#pi^{-}} thrown  (deg)");
                pim_theta_measured_vs_thrown_CD[i]->SetYTitle("#theta_{#pi^{-}} measured (deg)");
                //if (pim_theta_measured_vs_thrown_CD[i]->GetEntries())
                pim_theta_measured_vs_thrown_CD[i]->Write();
                //
                // pim_theta_rec_vs_thrown_CD[i]->SetOption("COLZ");
                // pim_theta_rec_vs_thrown_CD[i]->SetTitle("#theta_{#pi^{-}} rec vs #theta_{#pi^{-}} thrown ");
                // pim_theta_rec_vs_thrown_CD[i]->SetXTitle("#theta_{#pi^{-}} thrown  (deg)");
                // pim_theta_rec_vs_thrown_CD[i]->SetYTitle("#theta_{#pi^{-}} rec (deg)");
                // //if (pim_theta_rec_vs_thrown_CD[i]->GetEntries())
                // pim_theta_rec_vs_thrown_CD[i]->Write();
        }
}

void Histogram::write_hist_efficiency_check_FD() {


        auto pim_phi_vs_theta_rec_FD_can =
                std::make_unique<TCanvas>("pim_phi_vs_theta_rec_FD", "pim_phi_vs_theta_rec_FD sectors", 1920, 1080);
        pim_phi_vs_theta_rec_FD_can->Divide(3, 2);
        for (short i = 0; i < num_sectors; i++) {
                pim_phi_vs_theta_rec_FD_sec[i]->SetTitle("#phi_{#pi^{-}} vs #theta_{#pi^{-}} rec ");
                pim_phi_vs_theta_rec_FD_sec[i]->SetYTitle("#phi (deg)");
                pim_phi_vs_theta_rec_FD_sec[i]->SetXTitle("#theta (deg)");
                pim_phi_vs_theta_rec_FD_sec[i]->SetOption("COLZ");
                pim_phi_vs_theta_rec_FD_can->cd(i + 1);
                pim_phi_vs_theta_rec_FD_sec[i]->Draw("same");
        }
        pim_phi_vs_theta_rec_FD_can->Write();
        auto pim_phi_vs_theta_measured_FD_can =
                std::make_unique<TCanvas>("pim_phi_vs_theta_measured_FD", "pim_phi_vs_theta_measured_FD sectors", 1920, 1080);
        pim_phi_vs_theta_measured_FD_can->Divide(3, 2);
        for (short i = 0; i < num_sectors; i++) {
                pim_phi_vs_theta_measured_FD_sec[i]->SetTitle("#phi_{#pi^{-}} vs #theta_{#pi^{-}} measured ");
                pim_phi_vs_theta_measured_FD_sec[i]->SetYTitle("#phi (deg)");
                pim_phi_vs_theta_measured_FD_sec[i]->SetXTitle("#theta (deg)");
                pim_phi_vs_theta_measured_FD_sec[i]->SetOption("COLZ");
                pim_phi_vs_theta_measured_FD_can->cd(i + 1);
                pim_phi_vs_theta_measured_FD_sec[i]->Draw("same");
        }
        pim_phi_vs_theta_measured_FD_can->Write();

        auto pim_phi_vs_theta_rec_FD_after_MMSQ_can =
                std::make_unique<TCanvas>("pim_phi_vs_theta_rec_FD_after_MMSQ", "pim_phi_vs_theta_rec_FD_after_MMSQ sectors", 1920, 1080);
        pim_phi_vs_theta_rec_FD_after_MMSQ_can->Divide(3, 2);
        for (short i = 0; i < num_sectors; i++) {
                pim_phi_vs_theta_rec_FD_after_exclusive_sec[i]->SetTitle("#phi_{#pi^{-}} vs #theta_{#pi^{-}} rec ");
                pim_phi_vs_theta_rec_FD_after_exclusive_sec[i]->SetYTitle("#phi (deg)");
                pim_phi_vs_theta_rec_FD_after_exclusive_sec[i]->SetXTitle("#theta (deg)");
                pim_phi_vs_theta_rec_FD_after_exclusive_sec[i]->SetOption("COLZ");
                pim_phi_vs_theta_rec_FD_after_MMSQ_can->cd(i + 1);
                pim_phi_vs_theta_rec_FD_after_exclusive_sec[i]->Draw("same");
        }
        pim_phi_vs_theta_rec_FD_after_MMSQ_can->Write();
        auto pim_phi_vs_theta_measured_FD_after_MMSQ_can =
                std::make_unique<TCanvas>("pim_phi_vs_theta_measured_FD_after_MMSQ", "pim_phi_vs_theta_measured_FD_after_MMSQ sectors", 1920, 1080);
        pim_phi_vs_theta_measured_FD_after_MMSQ_can->Divide(3, 2);
        for (short i = 0; i < num_sectors; i++) {
                pim_phi_vs_theta_measured_FD_after_exclusive_sec[i]->SetTitle("#phi_{#pi^{-}} vs #theta_{#pi^{-}} measured ");
                pim_phi_vs_theta_measured_FD_after_exclusive_sec[i]->SetYTitle("#phi (deg)");
                pim_phi_vs_theta_measured_FD_after_exclusive_sec[i]->SetXTitle("#theta (deg)");
                pim_phi_vs_theta_measured_FD_after_exclusive_sec[i]->SetOption("COLZ");
                pim_phi_vs_theta_measured_FD_after_MMSQ_can->cd(i + 1);
                pim_phi_vs_theta_measured_FD_after_exclusive_sec[i]->Draw("same");
        }
        pim_phi_vs_theta_measured_FD_after_MMSQ_can->Write();

        for (size_t i = 0; i < EFF_CONDITIONS_NUM; i++) {

                MM2_exclusive_hist_FD[i]->SetXTitle("MMSQ (GeV)");
                MM2_exclusive_hist_FD[i]->Write();
                // pim_theta_measured_minus_rec_FD[i]->Fit("gaus", "QMR+", "QMR+", -2, 2);
                // gStyle->SetOptFit(1111);
                // pim_theta_measured_minus_rec_FD[i]->SetXTitle("#theta (deg)");
                // pim_theta_measured_minus_rec_FD[i]->Write();
                //
                // // pim_phi_measured_minus_rec_FD[i]->Fit("gaus", "QMR+", "QMR+", -2, 2);
                // pim_phi_measured_minus_rec_FD[i]->SetXTitle("#phi (deg)");
                // pim_phi_measured_minus_rec_FD[i]->Write();
                //
                // // pim_mom_measured_minus_rec_FD[i]->Fit("gaus", "QMR+", "QMR+", -2, 2);
                // pim_mom_measured_minus_rec_FD[i]->SetXTitle("P (GeV)");
                // pim_mom_measured_minus_rec_FD[i]->Write();


                pim_theta_rec_vs_mom_FD[i]->SetOption("COLZ");
                pim_theta_rec_vs_mom_FD[i]->Write();

                pim_phi_rec_vs_mom_FD[i]->SetOption("COLZ");
                pim_phi_rec_vs_mom_FD[i]->Write();


                pim_theta_measured_vs_mom_FD[i]->SetOption("COLZ");
                pim_theta_measured_vs_mom_FD[i]->Write();
                pim_phi_measured_vs_mom_FD[i]->SetOption("COLZ");
                pim_phi_measured_vs_mom_FD[i]->Write();

                pim_E_vs_mom_rec_FD[i]->SetOption("COLZ");
                pim_E_vs_mom_rec_FD[i]->Write();

                pim_E_vs_mom_measured_FD[i]->SetOption("COLZ");
                pim_E_vs_mom_measured_FD[i]->Write();

                pim_phi_vs_theta_rec_FD[i]->SetOption("COLZ");
                pim_phi_vs_theta_rec_FD[i]->SetTitle("#phi_{#pi^{-}} vs #theta_{#pi^{-}} measured ");
                pim_phi_vs_theta_rec_FD[i]->SetXTitle("#theta_{#pi^{-}} (deg)");
                pim_phi_vs_theta_rec_FD[i]->SetYTitle("#phi_{#pi^{-}} (deg)");
                pim_phi_vs_theta_rec_FD[i]->Write();


                pim_E_measured_vs_rec_FD[i]->SetOption("COLZ");
                pim_E_measured_vs_rec_FD[i]->SetTitle("E_{#pi^{-}} measured vs E_{#pi^{-}} rec ");
                pim_E_measured_vs_rec_FD[i]->SetXTitle("E_{#pi^{-}} rec (GeV)");
                pim_E_measured_vs_rec_FD[i]->SetYTitle("E_{#pi^{-}} measured (GeV)");
                if (pim_E_measured_vs_rec_FD[i]->GetEntries())
                        pim_E_measured_vs_rec_FD[i]->Write();


                pim_mom_measured_vs_rec_FD[i]->SetOption("COLZ");
                pim_mom_measured_vs_rec_FD[i]->SetTitle("P_{#pi^{-}} measured vs P_{#pi^{-}} rec ");
                pim_mom_measured_vs_rec_FD[i]->SetXTitle("P_{#pi^{-}} rec (GeV)");
                pim_mom_measured_vs_rec_FD[i]->SetYTitle("P_{#pi^{-}} measured (GeV)");
                if (pim_mom_measured_vs_rec_FD[i]->GetEntries())
                        pim_mom_measured_vs_rec_FD[i]->Write();

                pim_theta_measured_vs_rec_FD[i]->SetOption("COLZ");
                pim_theta_measured_vs_rec_FD[i]->SetTitle("#theta_{#pi^{-}} measured vs #theta_{#pi^{-}} rec ");
                pim_theta_measured_vs_rec_FD[i]->SetXTitle("#theta_{#pi^{-}} rec (deg)");
                pim_theta_measured_vs_rec_FD[i]->SetYTitle("#theta_{#pi^{-}} measured (deg)");
                if (pim_theta_measured_vs_rec_FD[i]->GetEntries())
                        pim_theta_measured_vs_rec_FD[i]->Write();

                pim_phi_measured_vs_rec_FD[i]->SetOption("COLZ");
                pim_phi_measured_vs_rec_FD[i]->SetTitle("#phi_{#pi^{-}} measured vs #phi_{#pi^{-}} rec ");
                pim_phi_measured_vs_rec_FD[i]->SetXTitle("#phi_{#pi^{-}} rec (deg)");
                pim_phi_measured_vs_rec_FD[i]->SetYTitle("#phi_{#pi^{-}} measured (deg)");
                if (pim_phi_measured_vs_rec_FD[i]->GetEntries())
                        pim_phi_measured_vs_rec_FD[i]->Write();

                pim_phi_vs_theta_measured_FD[i]->SetOption("COLZ");
                pim_phi_vs_theta_measured_FD[i]->SetOption("COLZ");
                pim_phi_vs_theta_measured_FD[i]->SetTitle("#phi_{#pi^{-}} vs #theta_{#pi^{-}} measured ");
                pim_phi_vs_theta_measured_FD[i]->SetXTitle("#theta_{#pi^{-}} (deg)");
                pim_phi_vs_theta_measured_FD[i]->SetYTitle("#phi_{#pi^{-}} (deg)");
                pim_phi_vs_theta_measured_FD[i]->Write();


                pim_theta_measured_minus_rec_vs_measured_FD[i]->SetOption("COLZ");
                pim_theta_measured_minus_rec_vs_measured_FD[i]->SetTitle("#theta_{#pi^{-}} measured - #theta_{#pi^{-}} rec vs #theta_{#pi^{-}} rec ");
                pim_theta_measured_minus_rec_vs_measured_FD[i]->SetXTitle("#theta_{#pi^{-}} rec  (deg)");
                pim_theta_measured_minus_rec_vs_measured_FD[i]->SetYTitle("#theta_{#pi^{-}} measured - #theta_{#pi^{-}} rec (deg)");
                //if (pim_theta_measured_minus_rec_vs_measured_FD[i]->GetEntries())
                pim_theta_measured_minus_rec_vs_measured_FD[i]->Write();

                // pim_theta_measured_minus_rec_vs_rec_FD[i]->SetOption("COLZ");
                // pim_theta_measured_minus_rec_vs_rec_FD[i]->SetTitle("#theta_{#pi^{-}} measured - #theta_{#pi^{-}} rec vs #theta_{#pi^{-}} rec ");
                // pim_theta_measured_minus_rec_vs_rec_FD[i]->SetXTitle("#theta_{#pi^{-}} measured  (deg)");
                // pim_theta_measured_minus_rec_vs_rec_FD[i]->SetYTitle("#theta_{#pi^{-}} measured - #theta_{#pi^{-}} rec (deg)");
                // //if (pim_theta_measured_minus_rec_vs_rec_FD[i]->GetEntries())
                // pim_theta_measured_minus_rec_vs_rec_FD[i]->Write();


                pim_theta_measured_vs_thrown_FD[i]->SetOption("COLZ");
                pim_theta_measured_vs_thrown_FD[i]->SetTitle("#theta_{#pi^{-}} measured - #theta_{#pi^{-}} thrown vs #theta_{#pi^{-}} thrown ");
                pim_theta_measured_vs_thrown_FD[i]->SetXTitle("#theta_{#pi^{-}} thrown  (deg)");
                pim_theta_measured_vs_thrown_FD[i]->SetYTitle("#theta_{#pi^{-}} measured (deg)");
                //if (pim_theta_measured_vs_thrown_FD[i]->GetEntries())
                pim_theta_measured_vs_thrown_FD[i]->Write();

                // pim_theta_rec_vs_thrown_FD[i]->SetOption("COLZ");
                // pim_theta_rec_vs_thrown_FD[i]->SetTitle("#theta_{#pi^{-}} rec vs #theta_{#pi^{-}} thrown ");
                // pim_theta_rec_vs_thrown_FD[i]->SetXTitle("#theta_{#pi^{-}} thrown  (deg)");
                // pim_theta_rec_vs_thrown_FD[i]->SetYTitle("#theta_{#pi^{-}} rec (deg)");
                // //if (pim_theta_rec_vs_thrown_FD[i]->GetEntries())
                // pim_theta_rec_vs_thrown_FD[i]->Write();

        }
}
void Histogram::Fill_WvsQ2(const std::shared_ptr<Reaction> &_e) {
        short sec = _e->sec();

        //if(_e->W() > 1.75 && _e->W() < 1.8) {
        //        if( _e->Q2()> 4.0 && _e->Q2()< 5.0) {
        // if (_e->TwoPion_missingPim()) {
        //
        //         //std::cout << " measured =  "<<_e->pim_theta_measured() << " Reconstructed =  "<< _e->pim_theta()<<endl;
        //         theta_pim_rec->Fill(_e->pim_theta_lab(),_e->weight());
        //         pim_theta_rec_vs_mom->Fill(_e->pim_momentum(),_e->pim_theta_lab(),_e->weight());
        //         phi_pim_rec->Fill(_e->pim_Phi_lab(),_e->weight());
        //         pim_phi_rec_vs_mom->Fill(_e->pim_momentum(),_e->pim_Phi_lab(),_e->weight());
        //         pim_theta_vs_phi_rec->Fill(_e->pim_Phi_lab(),_e->pim_theta_lab(),_e->weight());
        //
        //         //}
        //         if(_e->pim_momentum_measured()>0) {
        //                 theta_pim_measured->Fill(_e->pim_theta_lab_measured(),_e->weight());
        //                 pim_theta_measured_vs_rec->Fill(_e->pim_theta_lab(),_e->pim_theta_lab_measured(),_e->weight());
        //                 pim_theta_measured_vs_mom->Fill(_e->pim_momentum_measured(),_e->pim_theta_lab_measured(),_e->weight());
        //                 phi_pim_measured->Fill(_e->pim_Phi_lab_measured(),_e->weight());
        //                 pim_phi_measured_vs_rec->Fill(_e->pim_Phi_lab(),_e->pim_Phi_lab_measured(),_e->weight());
        //                 pim_phi_measured_vs_mom->Fill(_e->pim_momentum_measured(),_e->pim_Phi_lab_measured(),_e->weight());
        //                 pim_theta_vs_phi_measured->Fill(_e->pim_Phi_lab_measured(),_e->pim_theta_lab_measured(),_e->weight());
        //
        //         }
        //         if (_e->MM_cut()) {                 // mmsq > -0.06 && mmsq < 0.08)
        //                 theta_pim_rec_after_mmsq_applied->Fill(_e->pim_theta_lab(),_e->weight());
        //                 pim_theta_rec_vs_mom_after_mmsq_applied->Fill(_e->pim_momentum(),_e->pim_theta_lab(),_e->weight());
        //
        //                 phi_pim_rec_after_mmsq_applied->Fill(_e->pim_Phi_lab(),_e->weight());
        //                 pim_phi_rec_vs_mom_after_mmsq_applied->Fill(_e->pim_momentum(),_e->pim_Phi_lab(),_e->weight());
        //
        //                 pim_theta_vs_phi_rec_after_mmsq_applied->Fill(_e->pim_Phi_lab(),_e->pim_theta_lab(),_e->weight());
        //
        //                 //}
        //                 if(_e->pim_momentum_measured()>0) {
        //                         theta_pim_measured_after_mmsq_applied->Fill(_e->pim_theta_lab_measured(),_e->weight());
        //                         pim_theta_measured_vs_rec_after_mmsq_applied->Fill(_e->pim_theta_lab(),_e->pim_theta_lab_measured(),_e->weight());
        //                         if(abs(_e->pim_theta_lab() -_e->pim_theta_lab_measured()) < 10 ) {
        //                                 pim_theta_measured_vs_rec_after_y_equal_mc_cuts->Fill(_e->pim_theta_lab(),_e->pim_theta_lab_measured(),_e->weight());
        //                         }
        //                         pim_theta_measured_vs_mom_after_mmsq_applied->Fill(_e->pim_momentum_measured(),_e->pim_theta_lab_measured(),_e->weight());
        //                         phi_pim_measured_after_mmsq_applied->Fill(_e->pim_Phi_lab_measured(),_e->weight());
        //                         pim_phi_measured_vs_rec_after_mmsq_applied->Fill(_e->pim_Phi_lab(),_e->pim_Phi_lab_measured(),_e->weight());
        //                         if(abs(_e->pim_Phi_lab() -_e->pim_Phi_lab_measured()) < 20 ) {
        //
        //                                 pim_phi_measured_vs_rec_after_y_equal_mc_cuts->Fill(_e->pim_Phi_lab(),_e->pim_Phi_lab_measured(),_e->weight());
        //                         }
        //                         pim_phi_measured_vs_mom_after_mmsq_applied->Fill(_e->pim_momentum_measured(),_e->pim_Phi_lab_measured(),_e->weight());
        //                         pim_theta_vs_phi_measured_after_mmsq_applied->Fill(_e->pim_Phi_lab_measured(),_e->pim_theta_lab_measured(),_e->weight());
        //                 }
        //
        //         }
        // }
        //        }
        //}
        if (_e->TwoPion_missingPip())
                MM2_twoPi_missingPip->Fill(_e->MM2_mpip(), _e->weight());
        if (_e->TwoPion_missingProt())
                MM2_twoPi_missingProt->Fill(_e->MM2_mprot(), _e->weight());
        if (sec > 0 && sec <= 6) {
                W_vs_q2_sec[sec - 1]->Fill(_e->W(), _e->Q2(), _e->weight());
                W_sec[sec - 1]->Fill(_e->W()), _e->weight();
                if (_e->TwoPion_missingPip())
                        MM2_twoPi_missingPip_sec[sec - 1]->Fill(_e->MM2_mpip(), _e->weight());
                if (_e->TwoPion_missingProt())
                        MM2_twoPi_missingProt_sec[sec - 1]->Fill(_e->MM2_mprot(), _e->weight());
        }

        short det = _e->det();
        if (det == 1 && _e->W() <= 3.5) {
                W_det[0]->Fill(_e->W(), _e->weight());
                WQ2_det[0]->Fill(_e->W(), _e->Q2(), _e->weight());
        } else if (det == 2) {
                W_det[1]->Fill(_e->W(), _e->weight());
                WQ2_det[1]->Fill(_e->W(), _e->Q2(), _e->weight());
        } else {
                W_det[2]->Fill(_e->W(), _e->weight());
                WQ2_det[2]->Fill(_e->W(), _e->Q2(), _e->weight());
        }
}
// void Histogram::Fill_WvsQ2(const std::shared_ptr<MCReaction>& _e) {
//   W_vs_q2->Fill(_e->W(), _e->Q2(), _e->weight());
//   W_hist->Fill(_e->W(), _e->weight());
//   Q2_hist->Fill(_e->Q2(), _e->weight());
//
//   short sec = _e->sec();
//   if (sec > 0 && sec <= 6) {
//     W_vs_q2_sec[sec - 1]->Fill(_e->W(), _e->Q2(), _e->weight());
//     W_sec[sec - 1]->Fill(_e->W(), _e->weight());
//   }
//
//   short det = _e->det();
//   if (det == 1 && _e->W() <= 3.5) {
//     W_det[0]->Fill(_e->W(), _e->weight());
//     WQ2_det[0]->Fill(_e->W(), _e->Q2(), _e->weight());
//   } else if (det == 2) {
//     W_det[1]->Fill(_e->W(), _e->weight());
//     WQ2_det[1]->Fill(_e->W(), _e->Q2(), _e->weight());
//   } else {
//     W_det[2]->Fill(_e->W(), _e->weight());
//     WQ2_det[2]->Fill(_e->W(), _e->Q2(), _e->weight());
//   }
// }

// void Histogram::Fill_WvsQ2_singlePi(const std::shared_ptr<Reaction>& _e) {
//   short sec = _e->sec();
//
//   W_vs_q2_singlePi->Fill(_e->W(), _e->Q2(), _e->weight());
//   W_hist_singlePi->Fill(_e->W(), _e->weight());
//   Q2_hist_singlePi->Fill(_e->Q2(), _e->weight());
//   MM_neutron->Fill(_e->MM(), _e->weight());
//   if (sec > 0 && sec <= 6) {
//     W_vs_MM_singlePi[sec - 1]->Fill(_e->W(), _e->MM(), _e->weight());
//     W_vs_q2_singlePi_sec[sec - 1]->Fill(_e->W(), _e->Q2(), _e->weight());
//     W_singlePi_sec[sec - 1]->Fill(_e->W(), _e->weight());
//     MM_neutron_sec[sec - 1]->Fill(_e->MM(), _e->weight());
//   }
// }

void Histogram::Fill_WvsQ2_twoPi(const std::shared_ptr<Reaction> &_e) {
        short sec = _e->sec();
        weight_hist->Fill(_e->weight());
        if (_e->MM_cut()) { // abs(mmsq<0.03)

                // theta_prot_mc->Fill(_e->prot_theta(), _e->weight());
                // theta_pip_mc->Fill(_e->pip_theta(), _e->weight());
                // theta_pim_mc->Fill(_e->pim_theta(), _e->weight());
                //
                // Phi_gamma_mc->Fill(_e->gamma_Phi(), _e->weight());
                // Phi_prot_mc->Fill(_e->prot_Phi(), _e->weight());
                // Phi_pip_mc->Fill(_e->pip_Phi(), _e->weight());
                // Phi_pim_mc->Fill(_e->pim_Phi(), _e->weight());
                //
                // alpha_prot_mc->Fill(_e->alpha_pippim_pipf(), _e->weight());
                // alpha_pip_mc->Fill(_e->alpha_ppim_pipip(), _e->weight());
                // alpha_pim_mc->Fill(_e->alpha_ppip_pipim(), _e->weight());
                if (sec > 0 && sec <= 6) {
                        W_vs_MM_twoPi[sec - 1]->Fill(_e->W(), _e->MM(), _e->weight());
                        W_vs_q2_twoPi_sec[sec - 1]->Fill(_e->W(), _e->Q2(), _e->weight());
                        W_twoPi_sec[sec - 1]->Fill(_e->W(), _e->weight());
                        MM_twoPi_sec[sec - 1]->Fill(_e->MM(), _e->weight());
                        MM2_twoPi_sec[sec - 1]->Fill(_e->MM2(), _e->weight());
                }
                W_vs_q2_twoPi->Fill(_e->W(), _e->Q2(), _e->weight());
                W_hist_twoPi->Fill(_e->W(), _e->weight());
                Q2_hist_twoPi->Fill(_e->Q2(), _e->weight());

                MM_twoPi->Fill(_e->MM2(), _e->weight());

                MM2_twoPi->Fill(_e->MM2(), _e->weight());
                W_P2pi_hist->Fill(_e->w_P2pi_rec(), _e->weight());
                Q2_hist->Fill(_e->Q2(), _e->weight());
                W_vs_q2->Fill(_e->W(), _e->Q2(), _e->weight());

                theta_prot->Fill(_e->prot_theta(), _e->weight());
                theta_pip->Fill(_e->pip_theta(), _e->weight());
                theta_pim->Fill(_e->pim_theta(), _e->weight());
                if (_e->inv_pip_pim() > 0.7 && _e->inv_pip_pim() < 0.9) {
                        Theta_prot_cm_vs_mom_prot->Fill(_e->prot_momentum(), _e->prot_theta(),
                                                        _e->weight());
                        Theta_prot_lab_vs_mom_prot->Fill(_e->prot_momentum(),
                                                         _e->prot_theta_lab(), _e->weight());
                }
                if (_e->inv_Ppim() > 1.1 && _e->inv_Ppim() < 1.35) {
                        Theta_pip_cm_vs_mom_pip->Fill(_e->pip_momentum(), _e->pip_theta(),
                                                      _e->weight());
                        Theta_pip_lab_vs_mom_pip->Fill(_e->pip_momentum(), _e->pip_theta_lab(),
                                                       _e->weight());
                }
                if (_e->inv_Ppip() > 1.1 && _e->inv_Ppip() < 1.35) {
                        Theta_pim_cm_vs_mom_pim->Fill(_e->pim_momentum(), _e->pim_theta(),
                                                      _e->weight());
                        Theta_pim_lab_vs_mom_pim->Fill(_e->pim_momentum(), _e->pim_theta_lab(),
                                                       _e->weight());
                }

                if (_e->gamma_Phi() != 0)
                        Phi_gamma->Fill(_e->gamma_Phi(), _e->weight());
                if (_e->prot_Phi() != 0)
                        Phi_prot->Fill(_e->prot_Phi(), _e->weight());
                if (_e->pip_Phi() != 0)
                        Phi_pip->Fill(_e->pip_Phi(), _e->weight());
                if (_e->pim_Phi() != 0)
                        Phi_pim->Fill(_e->pim_Phi(), _e->weight());

                alpha_prot->Fill(_e->alpha_pippim_pipf(), _e->weight());
                alpha_pip->Fill(_e->alpha_ppim_pipip(), _e->weight());
                alpha_pim->Fill(_e->alpha_ppip_pipim(), _e->weight());

                inv_mass_P_pip[0][0] -> Fill(_e->inv_Ppip(), _e->weight());
                inv_mass_P_pim[0][0] -> Fill(_e->inv_Ppim(), _e->weight());
                inv_mass_pip_pim[0][0] -> Fill(_e->inv_pip_pim(), _e->weight());
                // std::cout << "theta_P_cm_ " << _e->p_mu_prime_cm().Theta() * 180 / PI <<
                // '\n';
                theta_P_vs_mass_pip_pim[0][0] -> Fill(_e->inv_pip_pim(), _e->prot_theta(),
                                                      _e->weight());
                theta_pip_vs_mass_Ppim[0][0] -> Fill(_e->inv_Ppim(), _e->pip_theta(),
                                                     _e->weight());
                theta_pim_vs_mass_Ppip[0][0] -> Fill(_e->inv_Ppip(), _e->pim_theta(),
                                                     _e->weight());

                theta_P_lab_vs_mass_pip_pim[0][0] -> Fill(_e->inv_pip_pim(),
                                                          _e->prot_theta_lab(), _e->weight());
                theta_pip_lab_vs_mass_Ppim[0][0] -> Fill(_e->inv_Ppim(), _e->pip_theta_lab(),
                                                         _e->weight());
                theta_pim_lab_vs_mass_Ppip[0][0] -> Fill(_e->inv_Ppip(), _e->pim_theta_lab(),
                                                         _e->weight());

                if (_e->Q2() < 4.5) {
                        inv_mass_P_pip[0][1] -> Fill(_e->inv_Ppip(), _e->weight());
                        inv_mass_P_pim[0][1] -> Fill(_e->inv_Ppim(), _e->weight());
                        inv_mass_pip_pim[0][1] -> Fill(_e->inv_pip_pim(), _e->weight());
                        // std::cout << "theta_P_cm_ " << _e->p_mu_prime_cm().Theta() * 180 / PI
                        // <<  '\n';
                        theta_P_vs_mass_pip_pim[0][1] -> Fill(_e->inv_pip_pim(), _e->prot_theta(),
                                                              _e->weight());
                        theta_pip_vs_mass_Ppim[0][1] -> Fill(_e->inv_Ppim(), _e->pip_theta(),
                                                             _e->weight());
                        theta_pim_vs_mass_Ppip[0][1] -> Fill(_e->inv_Ppip(), _e->pim_theta(),
                                                             _e->weight());

                        theta_P_lab_vs_mass_pip_pim[0][1] -> Fill(
                                _e->inv_pip_pim(), _e->prot_theta_lab(), _e->weight());
                        theta_pip_lab_vs_mass_Ppim[0][1] -> Fill(_e->inv_Ppim(),
                                                                 _e->pip_theta_lab(), _e->weight());
                        theta_pim_lab_vs_mass_Ppip[0][1] -> Fill(_e->inv_Ppip(),
                                                                 _e->pim_theta_lab(), _e->weight());
                } else if (_e->Q2() > 4.5) {
                        inv_mass_P_pip[0][2] -> Fill(_e->inv_Ppip(), _e->weight());
                        inv_mass_P_pim[0][2] -> Fill(_e->inv_Ppim(), _e->weight());
                        inv_mass_pip_pim[0][2] -> Fill(_e->inv_pip_pim(), _e->weight());
                        // std::cout << "theta_P_cm_ " << _e->p_mu_prime_cm().Theta() * 180 / PI
                        // <<  '\n';
                        theta_P_vs_mass_pip_pim[0][2] -> Fill(_e->inv_pip_pim(), _e->prot_theta(),
                                                              _e->weight());
                        theta_pip_vs_mass_Ppim[0][2] -> Fill(_e->inv_Ppim(), _e->pip_theta(),
                                                             _e->weight());
                        theta_pim_vs_mass_Ppip[0][2] -> Fill(_e->inv_Ppip(), _e->pim_theta(),
                                                             _e->weight());

                        theta_P_lab_vs_mass_pip_pim[0][2] -> Fill(
                                _e->inv_pip_pim(), _e->prot_theta_lab(), _e->weight());
                        theta_pip_lab_vs_mass_Ppim[0][2] -> Fill(_e->inv_Ppim(),
                                                                 _e->pip_theta_lab(), _e->weight());
                        theta_pim_lab_vs_mass_Ppip[0][2] -> Fill(_e->inv_Ppip(),
                                                                 _e->pim_theta_lab(), _e->weight());
                }
                if (_e->W() < 2.5) {
                        if (_e->Q2() < 4.5) {
                                inv_mass_P_pip[1][1] -> Fill(_e->inv_Ppip(), _e->weight());
                                inv_mass_P_pim[1][1] -> Fill(_e->inv_Ppim(), _e->weight());
                                inv_mass_pip_pim[1][1] -> Fill(_e->inv_pip_pim(), _e->weight());
                                // std::cout << "theta_P_cm_ " << _e->p_mu_prime_cm().Theta() * 180 / PI
                                // <<  '\n';
                                theta_P_vs_mass_pip_pim[1][1] -> Fill(_e->inv_pip_pim(), _e->prot_theta(),
                                                                      _e->weight());
                                theta_pip_vs_mass_Ppim[1][1] -> Fill(_e->inv_Ppim(), _e->pip_theta(),
                                                                     _e->weight());
                                theta_pim_vs_mass_Ppip[1][1] -> Fill(_e->inv_Ppip(), _e->pim_theta(),
                                                                     _e->weight());

                                theta_P_lab_vs_mass_pip_pim[1][1] -> Fill(
                                        _e->inv_pip_pim(), _e->prot_theta_lab(), _e->weight());
                                theta_pip_lab_vs_mass_Ppim[1][1] -> Fill(
                                        _e->inv_Ppim(), _e->pip_theta_lab(), _e->weight());
                                theta_pim_lab_vs_mass_Ppip[1][1] -> Fill(
                                        _e->inv_Ppip(), _e->pim_theta_lab(), _e->weight());
                        } else if (_e->Q2() > 4.5) {
                                inv_mass_P_pip[1][2] -> Fill(_e->inv_Ppip(), _e->weight());
                                inv_mass_P_pim[1][2] -> Fill(_e->inv_Ppim(), _e->weight());
                                inv_mass_pip_pim[1][2] -> Fill(_e->inv_pip_pim(), _e->weight());
                                // std::cout << "theta_P_cm_ " << _e->p_mu_prime_cm().Theta() * 180 / PI
                                // <<  '\n';
                                theta_P_vs_mass_pip_pim[1][2] -> Fill(_e->inv_pip_pim(), _e->prot_theta(),
                                                                      _e->weight());
                                theta_pip_vs_mass_Ppim[1][2] -> Fill(_e->inv_Ppim(), _e->pip_theta(),
                                                                     _e->weight());
                                theta_pim_vs_mass_Ppip[1][2] -> Fill(_e->inv_Ppip(), _e->pim_theta(),
                                                                     _e->weight());

                                theta_P_lab_vs_mass_pip_pim[1][2] -> Fill(
                                        _e->inv_pip_pim(), _e->prot_theta_lab(), _e->weight());
                                theta_pip_lab_vs_mass_Ppim[1][2] -> Fill(
                                        _e->inv_Ppim(), _e->pip_theta_lab(), _e->weight());
                                theta_pim_lab_vs_mass_Ppip[1][2] -> Fill(
                                        _e->inv_Ppip(), _e->pim_theta_lab(), _e->weight());
                        }
                } else if (_e->W() > 2.5) {
                        if (_e->Q2() < 4.5) {
                                inv_mass_P_pip[2][1] -> Fill(_e->inv_Ppip(), _e->weight());
                                inv_mass_P_pim[2][1] -> Fill(_e->inv_Ppim(), _e->weight());
                                inv_mass_pip_pim[2][1] -> Fill(_e->inv_pip_pim(), _e->weight());
                                theta_P_vs_mass_pip_pim[2][1] -> Fill(_e->inv_pip_pim(), _e->prot_theta(),
                                                                      _e->weight());
                                theta_pip_vs_mass_Ppim[2][1] -> Fill(_e->inv_Ppim(), _e->pip_theta(),
                                                                     _e->weight());
                                theta_pim_vs_mass_Ppip[2][1] -> Fill(_e->inv_Ppip(), _e->pim_theta(),
                                                                     _e->weight());

                                theta_P_lab_vs_mass_pip_pim[2][1] -> Fill(
                                        _e->inv_pip_pim(), _e->prot_theta_lab(), _e->weight());
                                theta_pip_lab_vs_mass_Ppim[2][1] -> Fill(
                                        _e->inv_Ppim(), _e->pip_theta_lab(), _e->weight());
                                theta_pim_lab_vs_mass_Ppip[2][1] -> Fill(
                                        _e->inv_Ppip(), _e->pim_theta_lab(), _e->weight());
                        } else if (_e->Q2() > 4.5) {
                                inv_mass_P_pip[2][2] -> Fill(_e->inv_Ppip(), _e->weight());
                                inv_mass_P_pim[2][2] -> Fill(_e->inv_Ppim(), _e->weight());
                                inv_mass_pip_pim[2][2] -> Fill(_e->inv_pip_pim(), _e->weight());
                                // std::cout << "theta_P_cm_ " << _e->p_mu_prime_cm().Theta() * 180 / PI
                                // <<  '\n';
                                theta_P_vs_mass_pip_pim[2][2] -> Fill(_e->inv_pip_pim(), _e->prot_theta(),
                                                                      _e->weight());
                                theta_pip_vs_mass_Ppim[2][2] -> Fill(_e->inv_Ppim(), _e->pip_theta(),
                                                                     _e->weight());
                                theta_pim_vs_mass_Ppip[2][2] -> Fill(_e->inv_Ppip(), _e->pim_theta(),
                                                                     _e->weight());

                                theta_P_lab_vs_mass_pip_pim[2][2] -> Fill(
                                        _e->inv_pip_pim(), _e->prot_theta_lab(), _e->weight());
                                theta_pip_lab_vs_mass_Ppim[2][2] -> Fill(
                                        _e->inv_Ppim(), _e->pip_theta_lab(), _e->weight());
                                theta_pim_lab_vs_mass_Ppip[2][2] -> Fill(
                                        _e->inv_Ppip(), _e->pim_theta_lab(), _e->weight());
                        }
                }
        }
}
void Histogram::Fill_WvsQ2_twoPi_thrown(const std::shared_ptr<MCReaction> &_e) {
        short sec = _e->sec();

        W_vs_q2_twoPi_thrown->Fill(_e->W_mc(), _e->Q2_mc(), _e->weight());
        // W_hist_twoPi_thrown->Fill(_e->W_mc(), _e->weight());
        // Q2_hist_twoPi_thrown->Fill(_e->Q2_mc(), _e->weight());
        // MM_twoPi_thrown->Fill(_e->MM_mc(), _e->weight());
        // MM2_twoPi_thrown->Fill(_e->MM2_mc(), _e->weight());
        // W_vs_Q2_thrown->Fill(_e->W_mc(), _e->Q2_mc(), _e->weight());
        // W_thrown->Fill(_e->W_mc(), _e->weight());

        theta_prot_thrown->Fill(_e->MCprot_theta_thrown(), _e->weight());
        theta_pip_thrown->Fill(_e->MCpip_theta_thrown(), _e->weight());
        theta_pim_thrown->Fill(_e->MCpim_theta_thrown(), _e->weight());

        if (_e->MCinv_pip_pim() > 0.7 && _e->MCinv_pip_pim() < 0.9) {
                Theta_prot_thrown_cm_vs_mom_prot->Fill(
                        _e->prot_momentum_thrown(), _e->MCprot_theta_thrown(), _e->weight());
                Theta_prot_thrown_lab_vs_mom_prot->Fill(
                        _e->prot_momentum_thrown(), _e->MCprot_theta_lab(), _e->weight());
        }
        if (_e->MCinv_Ppim() > 1.1 && _e->MCinv_Ppim() < 1.35) {
                Theta_pip_thrown_cm_vs_mom_pip->Fill(
                        _e->pip_momentum_thrown(), _e->MCpip_theta_thrown(), _e->weight());
                Theta_pip_thrown_lab_vs_mom_pip->Fill(_e->pip_momentum_thrown(),
                                                      _e->MCpip_theta_lab(), _e->weight());
        }
        if (_e->MCinv_Ppip() > 1.1 && _e->MCinv_Ppip() < 1.35) {
                Theta_pim_thrown_cm_vs_mom_pim->Fill(
                        _e->pim_momentum_thrown(), _e->MCpim_theta_thrown(), _e->weight());
                Theta_pim_thrown_lab_vs_mom_pim->Fill(_e->pim_momentum_thrown(),
                                                      _e->MCpim_theta_lab(), _e->weight());
        }
        Phi_gamma_thrown->Fill(_e->MCgamma_Phi_thrown(), _e->weight());
        Phi_prot_thrown->Fill(_e->MCprot_Phi_thrown(), _e->weight());
        Phi_pip_thrown->Fill(_e->MCpip_Phi_thrown(), _e->weight());
        Phi_pim_thrown->Fill(_e->MCpim_Phi_thrown(), _e->weight());

        alpha_prot_thrown->Fill(_e->MCalpha_pippim_pipf_thrown(), _e->weight());
        alpha_pip_thrown->Fill(_e->MCalpha_ppim_pipip_thrown(), _e->weight());
        alpha_pim_thrown->Fill(_e->MCalpha_ppip_pipim_thrown(), _e->weight());

        if (sec > 0 && sec <= 6) {
                W_vs_MM_twoPi_thrown[sec - 1] -> Fill(_e->W_mc(), _e->MM_mc(),
                                                      _e->weight());
                W_vs_q2_twoPi_sec_thrown[sec - 1] -> Fill(_e->W_mc(),
                                                          _e->Q2_mc(), _e->weight());
                W_twoPi_sec_thrown[sec - 1] -> Fill(_e->W_mc(),
                                                    _e->weight());
                MM_twoPi_sec_thrown[sec - 1] -> Fill(_e->MM_mc(),
                                                     _e->weight());
                MM2_twoPi_sec_thrown[sec - 1] -> Fill(_e->MM2_mc(),
                                                      _e->weight());
        }
        // inv_mass_P_pip_thrown[0][0]->Fill(_e->MCinv_Ppip(), _e->weight());
        // inv_mass_P_pim_thrown[0][0]->Fill(_e->MCinv_Ppim(), _e->weight());
        // inv_mass_pip_pim_thrown[0][0]->Fill(_e->MCinv_pip_pim(), _e->weight());
        // theta_P_vs_mass_pip_pim_thrown[0][0]->Fill(_e->MCinv_pip_pim(),
        // _e->MCprot_theta_thrown(), _e->weight());
        // theta_pip_vs_mass_Ppim_thrown[0][0]->Fill(_e->MCinv_Ppim(),
        // _e->MCpip_theta_thrown(), _e->weight());
        // theta_pim_vs_mass_Ppip_thrown[0][0]->Fill(_e->MCinv_Ppip(),
        // _e->MCpim_theta_thrown(), _e->weight());
        // theta_P_lab_vs_mass_pip_pim_thrown[0][0]->Fill(_e->MCinv_pip_pim(),
        // _e->MCprot_theta_lab(), _e->weight());
        // theta_pip_lab_vs_mass_Ppim_thrown[0][0]->Fill(_e->MCinv_Ppim(),
        // _e->MCpip_theta_lab(), _e->weight());
        // theta_pim_lab_vs_mass_Ppip_thrown[0][0]->Fill(_e->MCinv_Ppip(),
        // _e->MCpim_theta_lab(), _e->weight());
        //
        // if (_e->Q2_mc() < 4.5) {
        //   inv_mass_P_pip_thrown[0][1]->Fill(_e->MCinv_Ppip(), _e->weight());
        //   inv_mass_P_pim_thrown[0][1]->Fill(_e->MCinv_Ppim(), _e->weight());
        //   inv_mass_pip_pim_thrown[0][1]->Fill(_e->MCinv_pip_pim(), _e->weight());
        //   theta_P_vs_mass_pip_pim_thrown[0][1]->Fill(_e->MCinv_pip_pim(),
        //   _e->MCprot_theta_thrown(), _e->weight());
        //   theta_pip_vs_mass_Ppim_thrown[0][1]->Fill(_e->MCinv_Ppim(),
        //   _e->MCpip_theta_thrown(), _e->weight());
        //   theta_pim_vs_mass_Ppip_thrown[0][1]->Fill(_e->MCinv_Ppip(),
        //   _e->MCpim_theta_thrown(), _e->weight());
        //   theta_P_lab_vs_mass_pip_pim_thrown[0][1]->Fill(_e->MCinv_pip_pim(),
        //   _e->MCprot_theta_lab(), _e->weight());
        //   theta_pip_lab_vs_mass_Ppim_thrown[0][1]->Fill(_e->MCinv_Ppim(),
        //   _e->MCpip_theta_lab(), _e->weight());
        //   theta_pim_lab_vs_mass_Ppip_thrown[0][1]->Fill(_e->MCinv_Ppip(),
        //   _e->MCpim_theta_lab(), _e->weight());
        // } else if (_e->Q2_mc() > 4.5) {
        //   inv_mass_P_pip_thrown[0][2]->Fill(_e->MCinv_Ppip(), _e->weight());
        //   inv_mass_P_pim_thrown[0][2]->Fill(_e->MCinv_Ppim(), _e->weight());
        //   inv_mass_pip_pim_thrown[0][2]->Fill(_e->MCinv_pip_pim(), _e->weight());
        //   theta_P_vs_mass_pip_pim_thrown[0][2]->Fill(_e->MCinv_pip_pim(),
        //   _e->MCprot_theta_thrown(), _e->weight());
        //   theta_pip_vs_mass_Ppim_thrown[0][2]->Fill(_e->MCinv_Ppim(),
        //   _e->MCpip_theta_thrown(), _e->weight());
        //   theta_pim_vs_mass_Ppip_thrown[0][2]->Fill(_e->MCinv_Ppip(),
        //   _e->MCpim_theta_thrown(), _e->weight());
        //   theta_P_lab_vs_mass_pip_pim_thrown[0][2]->Fill(_e->MCinv_pip_pim(),
        //   _e->MCprot_theta_lab(), _e->weight());
        //   theta_pip_lab_vs_mass_Ppim_thrown[0][2]->Fill(_e->MCinv_Ppim(),
        //   _e->MCpip_theta_lab(), _e->weight());
        //   theta_pim_lab_vs_mass_Ppip_thrown[0][2]->Fill(_e->MCinv_Ppip(),
        //   _e->MCpim_theta_lab(), _e->weight());
        // }
        //
        // if (_e->W_mc() < 2.5) {
        //   if (_e->Q2_mc() < 4.5) {
        //     inv_mass_P_pip_thrown[1][1]->Fill(_e->MCinv_Ppip(), _e->weight());
        //     inv_mass_P_pim_thrown[1][1]->Fill(_e->MCinv_Ppim(), _e->weight());
        //     inv_mass_pip_pim_thrown[1][1]->Fill(_e->MCinv_pip_pim(), _e->weight());
        //     theta_P_vs_mass_pip_pim_thrown[1][1]->Fill(_e->MCinv_pip_pim(),
        //     _e->MCprot_theta_thrown(), _e->weight());
        //     theta_pip_vs_mass_Ppim_thrown[1][1]->Fill(_e->MCinv_Ppim(),
        //     _e->MCpip_theta_thrown(), _e->weight());
        //     theta_pim_vs_mass_Ppip_thrown[1][1]->Fill(_e->MCinv_Ppip(),
        //     _e->MCpim_theta_thrown(), _e->weight());
        //     theta_P_lab_vs_mass_pip_pim_thrown[1][1]->Fill(_e->MCinv_pip_pim(),
        //     _e->MCprot_theta_lab(), _e->weight());
        //     theta_pip_lab_vs_mass_Ppim_thrown[1][1]->Fill(_e->MCinv_Ppim(),
        //     _e->MCpip_theta_lab(), _e->weight());
        //     theta_pim_lab_vs_mass_Ppip_thrown[1][1]->Fill(_e->MCinv_Ppip(),
        //     _e->MCpim_theta_lab(), _e->weight());
        //   } else if (_e->Q2_mc() > 4.5) {
        //     inv_mass_P_pip_thrown[1][2]->Fill(_e->MCinv_Ppip(), _e->weight());
        //     inv_mass_P_pim_thrown[1][2]->Fill(_e->MCinv_Ppim(), _e->weight());
        //     inv_mass_pip_pim_thrown[1][2]->Fill(_e->MCinv_pip_pim(), _e->weight());
        //     theta_P_vs_mass_pip_pim_thrown[1][2]->Fill(_e->MCinv_pip_pim(),
        //     _e->MCprot_theta_thrown(), _e->weight());
        //     theta_pip_vs_mass_Ppim_thrown[1][2]->Fill(_e->MCinv_Ppim(),
        //     _e->MCpip_theta_thrown(), _e->weight());
        //     theta_pim_vs_mass_Ppip_thrown[1][2]->Fill(_e->MCinv_Ppip(),
        //     _e->MCpim_theta_thrown(), _e->weight());
        //     theta_P_lab_vs_mass_pip_pim_thrown[1][2]->Fill(_e->MCinv_pip_pim(),
        //     _e->MCprot_theta_lab(), _e->weight());
        //     theta_pip_lab_vs_mass_Ppim_thrown[1][2]->Fill(_e->MCinv_Ppim(),
        //     _e->MCpip_theta_lab(), _e->weight());
        //     theta_pim_lab_vs_mass_Ppip_thrown[1][2]->Fill(_e->MCinv_Ppip(),
        //     _e->MCpim_theta_lab(), _e->weight());
        //   }
        //
        // } else if (_e->W_mc() > 2.5) {
        //   if (_e->Q2_mc() < 4.5) {
        //     inv_mass_P_pip_thrown[2][1]->Fill(_e->MCinv_Ppip(), _e->weight());
        //     inv_mass_P_pim_thrown[2][1]->Fill(_e->MCinv_Ppim(), _e->weight());
        //     inv_mass_pip_pim_thrown[2][1]->Fill(_e->MCinv_pip_pim(), _e->weight());
        //     // std::cout << "inv_macc_pip_pim_mc " << _e->MCinv_pip_pim() << '\n';
        //     theta_P_vs_mass_pip_pim_thrown[2][1]->Fill(_e->MCinv_pip_pim(),
        //     _e->MCprot_theta_thrown(), _e->weight());
        //     theta_pip_vs_mass_Ppim_thrown[2][1]->Fill(_e->MCinv_Ppim(),
        //     _e->MCpip_theta_thrown(), _e->weight());
        //     theta_pim_vs_mass_Ppip_thrown[2][1]->Fill(_e->MCinv_Ppip(),
        //     _e->MCpim_theta_thrown(), _e->weight());
        //     theta_P_lab_vs_mass_pip_pim_thrown[2][1]->Fill(_e->MCinv_pip_pim(),
        //     _e->MCprot_theta_lab(), _e->weight());
        //     theta_pip_lab_vs_mass_Ppim_thrown[2][1]->Fill(_e->MCinv_Ppim(),
        //     _e->MCpip_theta_lab(), _e->weight());
        //     theta_pim_lab_vs_mass_Ppip_thrown[2][1]->Fill(_e->MCinv_Ppip(),
        //     _e->MCpim_theta_lab(), _e->weight());
        //   } else if (_e->Q2_mc() > 4.5) {
        //     inv_mass_P_pip_thrown[2][2]->Fill(_e->MCinv_Ppip(), _e->weight());
        //     inv_mass_P_pim_thrown[2][2]->Fill(_e->MCinv_Ppim(), _e->weight());
        //     inv_mass_pip_pim_thrown[2][2]->Fill(_e->MCinv_pip_pim(), _e->weight());
        //     theta_P_vs_mass_pip_pim_thrown[2][2]->Fill(_e->MCinv_pip_pim(),
        //     _e->MCprot_theta_thrown(), _e->weight());
        //     theta_pip_vs_mass_Ppim_thrown[2][2]->Fill(_e->MCinv_Ppim(),
        //     _e->MCpip_theta_thrown(), _e->weight());
        //     theta_pim_vs_mass_Ppip_thrown[2][2]->Fill(_e->MCinv_Ppip(),
        //     _e->MCpim_theta_thrown(), _e->weight());
        //     theta_P_lab_vs_mass_pip_pim_thrown[2][2]->Fill(_e->MCinv_pip_pim(),
        //     _e->MCprot_theta_lab(), _e->weight());
        //     theta_pip_lab_vs_mass_Ppim_thrown[2][2]->Fill(_e->MCinv_Ppim(),
        //     _e->MCpip_theta_lab(), _e->weight());
        //     theta_pim_lab_vs_mass_Ppip_thrown[2][2]->Fill(_e->MCinv_Ppip(),
        //     _e->MCpim_theta_lab(), _e->weight());
        //   }
        // }
}

// W and Q^2
// void Histogram::Fill_WvsQ2_Npip(const std::shared_ptr<Reaction>& _e) {
//   short sec = _e->sec();
//   if (sec > 0 && sec <= 6) {
//     W_vs_q2_Npip_sec[sec - 1]->Fill(_e->W(), _e->Q2(), _e->weight());
//     W_Npip_sec[sec - 1]->Fill(_e->W(), _e->weight());
//     MM_Npip_sec[sec - 1]->Fill(_e->MM(), _e->weight());
//   }
// }

void Histogram::Write_WvsQ2() {

        weight_hist->SetXTitle("weight");
        weight_hist->Write();
        W_vs_q2->SetYTitle("Q^{2} (GeV^{2})");
        W_vs_q2->SetXTitle("W (GeV)");
        W_vs_q2->SetOption("COLZ");
        W_vs_q2->Write();
        for (short i = 0; i < 3; i ++) {
                WQ2_det[i] -> SetXTitle("W (GeV)");
                WQ2_det[i] -> SetYTitle("Q^{2} (GeV^2)");
                WQ2_det[i] -> SetOption("COLZ");
                if (WQ2_det[i] -> GetEntries())
                        WQ2_det[i] -> Write();
                W_det[i] -> SetXTitle("W (GeV)");
                if (W_det[i] -> GetEntries())
                        W_det[i] -> Write();
        }



        auto WvsQ2_can =
                std::make_unique<TCanvas>("WvsQ2_can", "W vs Q2 sectors", 1920, 1080);
        WvsQ2_can->Divide(3, 2);
        for (short i = 0; i < num_sectors; i ++) {
                W_vs_q2_sec[i] -> SetYTitle("Q^{2} (GeV^{2})");
                W_vs_q2_sec[i] -> SetXTitle("W (GeV)");
                W_vs_q2_sec[i] -> SetOption("COLZ1");
                WvsQ2_can->cd(i + 1);
                W_vs_q2_sec[i] -> Draw("same");
        }
        WvsQ2_can->Write();

        auto W_can = std::make_unique<TCanvas>("W_can", "W sectors", 1920, 1080);
        W_can->Divide(3, 2);
        for (short i = 0; i < num_sectors; i ++) {
                W_sec[i] -> SetXTitle("W (GeV)");
                W_can->cd(i + 1);

                //  W_sec[i]->Fit("gaus", "QMR+", "QMR+", 0.85, 1.05);
                // gStyle->SetOptFit(01);
                W_sec[i] -> Draw("same");
        }
        W_can->Write();

        W_vs_q2->SetXTitle("W (GeV)");
        W_vs_q2->SetYTitle("Q^{2} (GeV^{2})");
        W_vs_q2->SetOption("COLZ1");
        if (W_vs_q2->GetEntries())
                W_vs_q2->Write();

        W_hist->SetXTitle("W (GeV)");
        if (W_hist->GetEntries())
                W_hist->Write();

        W_P2pi_hist->SetXTitle("W_P2pi (GeV)");
        if (W_P2pi_hist->GetEntries())
                W_P2pi_hist->Write();

        Q2_hist->SetXTitle("Q^{2} (GeV^{2})");
        if (Q2_hist->GetEntries())
                Q2_hist->Write();

        W_vs_Q2_thrown->SetXTitle("W (GeV)");
        W_vs_Q2_thrown->SetYTitle("Q^{2} (GeV^{2})");
        W_vs_Q2_thrown->SetOption("COLZ1");
        if (W_vs_Q2_thrown->GetEntries())
                W_vs_Q2_thrown->Write();

        W_thrown->SetXTitle("W_thrown (GeV)");
        if (W_thrown->GetEntries())
                W_thrown->Write();

        W_vs_q2_twoPi->SetXTitle("W (GeV)");
        W_vs_q2_twoPi->SetYTitle("Q^{2} (GeV^{2})");
        W_vs_q2_twoPi->SetOption("COLZ1");
        if (W_vs_q2_twoPi->GetEntries())
                W_vs_q2_twoPi->Write();

        W_vs_q2_twoPi_thrown->SetXTitle("W_thrown (GeV)");
        W_vs_q2_twoPi_thrown->SetYTitle("Q^{2} (GeV^{2})");
        W_vs_q2_twoPi_thrown->SetOption("COLZ1");
        if (W_vs_q2_twoPi_thrown->GetEntries())
                W_vs_q2_twoPi_thrown->Write();

        W_hist_twoPi->SetXTitle("W (GeV)");
        if (W_hist_twoPi->GetEntries())
                W_hist_twoPi->Write();

        Q2_hist_twoPi->SetXTitle("Q^{2} (GeV^{2})");
        if (Q2_hist_twoPi->GetEntries())
                Q2_hist_twoPi->Write();

        MM_twoPi->SetXTitle("MM2 (GeV2)");
        if (MM_twoPi->GetEntries())
                MM_twoPi->Write();
        MM2_twoPi->SetXTitle("MM2 (GeV2)");
        if (MM2_twoPi->GetEntries())
                MM2_twoPi->Write();

        MM2_twoPi_missingPip->SetXTitle("MM2 (GeV2)");
        if (MM2_twoPi_missingPip->GetEntries())
                MM2_twoPi_missingPip->Write();

        MM2_twoPi_missingProt->SetXTitle("MM2 (GeV2)");
        if (MM2_twoPi_missingProt->GetEntries())
                MM2_twoPi_missingProt->Write();

        auto wvsq2_sec = RootOutputFile->mkdir("wvsq2_sec");
        wvsq2_sec->cd();
        for (short i = 0; i < num_sectors; i ++) {
                W_vs_q2_sec[i] -> SetYTitle("Q^{2} (GeV^{2})");
                W_vs_q2_sec[i] -> SetXTitle("W (GeV)");
                W_vs_q2_sec[i] -> SetOption("COLZ1");
                W_vs_q2_sec[i] -> Write();
        }
        auto w_sec = RootOutputFile->mkdir("w_sec");
        w_sec->cd();
        for (short i = 0; i < num_sectors; i ++) {
                W_sec[i] -> SetXTitle("W (GeV)");

                //  W_sec[i]->Fit("gaus", "QMR+", "QMR+", 0.85, 1.05);
                // gStyle->SetOptFit(01);
                W_sec[i] -> Write();
        }
        // W_vs_q2_singlePi->SetXTitle("W (GeV)");
        // W_vs_q2_singlePi->SetYTitle("Q^{2} (GeV^{2})");
        // W_vs_q2_singlePi->SetOption("COLZ1");
        // if (W_vs_q2_singlePi->GetEntries()) W_vs_q2_singlePi->Write();
        //
        // W_hist_singlePi->SetXTitle("W (GeV)");
        // if (W_hist_singlePi->GetEntries()) W_hist_singlePi->Write();
        //
        // Q2_hist_singlePi->SetXTitle("Q^{2} (GeV^{2})");
        // if (Q2_hist_singlePi->GetEntries()) Q2_hist_singlePi->Write();
        //
        // if (MM_neutron->GetEntries()) MM_neutron->Write();
        // auto efficiency = RootOutputFile->mkdir("efficiency");
        // efficiency->cd();
        // theta_pim_rec->SetTitle("#theta_{#pi^{-}} rec ");
        // theta_pim_rec->SetXTitle("#theta_{#pi^{-}} rec (deg)");
        // theta_pim_rec->Write();
        // theta_pim_measured->SetTitle("#theta_{#pi^{-}} measured ");
        // theta_pim_measured->SetXTitle("#theta_{#pi^{-}} measured (deg)");
        // theta_pim_measured->Write();
        //
        // if (pim_theta_measured_vs_rec->GetEntries()) {
        //         pim_theta_measured_vs_rec->SetTitle("#theta_{#pi^{-}} measured vs #theta_{#pi^{-}} rec");
        //         pim_theta_measured_vs_rec->SetXTitle("#theta_{#pi^{-}} rec (deg)");
        //         pim_theta_measured_vs_rec->SetYTitle("#theta_{#pi^{-}} measured (deg)");
        //         pim_theta_measured_vs_rec->SetOption("COLZ1");
        //         pim_theta_measured_vs_rec->Write();
        // }
        // if (pim_theta_rec_vs_mom->GetEntries()) {
        //         pim_theta_rec_vs_mom->SetTitle("#theta_{#pi^{-}} rec vs momentum");
        //         pim_theta_rec_vs_mom->SetXTitle("P_{#pi^{-}} rec (GeV)");
        //         pim_theta_rec_vs_mom->SetYTitle("#theta_{#pi^{-}} rec (deg)");
        //         pim_theta_rec_vs_mom->SetOption("COLZ1");
        //         pim_theta_rec_vs_mom->Write();
        // }
        //
        // if (pim_theta_measured_vs_mom->GetEntries()) {
        //         pim_theta_measured_vs_mom->SetTitle("#theta_{#pi^{-}} measured vs momentum");
        //         pim_theta_measured_vs_mom->SetXTitle("P_{#pi^{-}} measured (GeV)");
        //         pim_theta_measured_vs_mom->SetYTitle("#theta_{#pi^{-}} measured (deg)");
        //         pim_theta_measured_vs_mom->SetOption("COLZ1");
        //         pim_theta_measured_vs_mom->Write();
        // }
        //
        // phi_pim_rec->SetTitle("#phi_{#pi^{-}} rec ");
        // phi_pim_rec->SetXTitle("#phi_{#pi^{-}} rec (deg)");
        // phi_pim_rec->Write();
        // phi_pim_measured->SetTitle("#phi_{#pi^{-}} measured (deg)");
        // phi_pim_measured->SetXTitle("#phi_{#pi^{-}} measured (deg)");
        // phi_pim_measured->Write();
        //
        // if (pim_phi_measured_vs_rec->GetEntries()) {
        //         pim_phi_measured_vs_rec->SetTitle("#phi_{#pi^{-}} measured vs #phi_{#pi^{-}} rec ");
        //         pim_phi_measured_vs_rec->SetXTitle("#phi_{#pi^{-}} rec (deg)");
        //         pim_phi_measured_vs_rec->SetYTitle("#phi_{#pi^{-}} measured (deg)");
        //         pim_phi_measured_vs_rec->SetOption("COLZ1");
        //         pim_phi_measured_vs_rec->Write();
        // }
        //
        // if (pim_phi_rec_vs_mom->GetEntries()) {
        //
        //         pim_phi_rec_vs_mom->SetTitle("#phi_{#pi^{-}} rec vs P_{#pi^{-}}");
        //         pim_phi_rec_vs_mom->SetXTitle("P_{#pi^{-}} rec (GeV)");
        //         pim_phi_rec_vs_mom->SetYTitle("#phi_{#pi^{-}} rec (deg)");
        //         pim_phi_rec_vs_mom->SetOption("COLZ1");
        //         pim_phi_rec_vs_mom->Write();
        // }
        //
        // if (pim_phi_measured_vs_mom->GetEntries()) {
        //
        //         pim_phi_measured_vs_mom->SetTitle("#phi_{#pi^{-}} measured vs P_{#pi^{-}}");
        //         pim_phi_measured_vs_mom->SetXTitle("P_{#pi^{-}} measured (GeV)");
        //         pim_phi_measured_vs_mom->SetYTitle("#phi_{#pi^{-}} measured (deg)");
        //         pim_phi_measured_vs_mom->SetOption("COLZ1");
        //         pim_phi_measured_vs_mom->Write();
        // }
        //
        // pim_theta_vs_phi_rec->SetTitle("#theta_{#pi^{-}} vs #phi_{#pi^{-}} rec");
        // pim_theta_vs_phi_rec->SetXTitle("#phi_{#pi^{-}} rec (GeV)");
        // pim_theta_vs_phi_rec->SetYTitle("#theta_{#pi^{-}} rec (deg)");
        // pim_theta_vs_phi_rec->SetOption("COLZ1");
        // pim_theta_vs_phi_rec->Write();
        //
        // pim_theta_vs_phi_measured->SetTitle("#theta_{#pi^{-}} vs #phi_{#pi^{-}} measured");
        // pim_theta_vs_phi_measured->SetXTitle("#phi_{#pi^{-}} measured (GeV)");
        // pim_theta_vs_phi_measured->SetYTitle("#theta_{#pi^{-}} measured (deg)");
        // pim_theta_vs_phi_measured->SetOption("COLZ1");
        // pim_theta_vs_phi_measured->Write();
        //
        //
        // theta_pim_rec_after_mmsq_applied->SetTitle("#theta_{#pi^{-}} rec ");
        // theta_pim_rec_after_mmsq_applied->SetXTitle("#theta_{#pi^{-}} rec (deg)");
        // theta_pim_rec_after_mmsq_applied->Write();
        // theta_pim_measured_after_mmsq_applied->SetTitle("#theta_{#pi^{-}} measured ");
        // theta_pim_measured_after_mmsq_applied->SetXTitle("#theta_{#pi^{-}} measured (deg)");
        // theta_pim_measured_after_mmsq_applied->Write();
        //
        // pim_theta_measured_vs_rec_after_mmsq_applied->SetTitle("#theta_{#pi^{-}} measured vs #theta_{#pi^{-}} rec");
        // pim_theta_measured_vs_rec_after_mmsq_applied->SetXTitle("#theta_{#pi^{-}} rec (deg)");
        // pim_theta_measured_vs_rec_after_mmsq_applied->SetYTitle("#theta_{#pi^{-}} measured (deg)");
        // pim_theta_measured_vs_rec_after_mmsq_applied->SetOption("COLZ1");
        // pim_theta_measured_vs_rec_after_mmsq_applied->Write();
        //
        // if (pim_theta_measured_vs_rec_after_y_equal_mc_cuts->GetEntries()) {
        //         pim_theta_measured_vs_rec_after_y_equal_mc_cuts->SetTitle("#theta_{#pi^{-}} measured vs #theta_{#pi^{-}} rec");
        //         pim_theta_measured_vs_rec_after_y_equal_mc_cuts->SetXTitle("#theta_{#pi^{-}} rec (deg)");
        //         pim_theta_measured_vs_rec_after_y_equal_mc_cuts->SetYTitle("#theta_{#pi^{-}} measured (deg)");
        //         pim_theta_measured_vs_rec_after_y_equal_mc_cuts->SetOption("COLZ1");
        //         pim_theta_measured_vs_rec_after_y_equal_mc_cuts->Write();
        // }
        //
        // pim_theta_rec_vs_mom_after_mmsq_applied->SetTitle("#theta_{#pi^{-}} rec vs momentum");
        // pim_theta_rec_vs_mom_after_mmsq_applied->SetXTitle("P_{#pi^{-}} rec (GeV)");
        // pim_theta_rec_vs_mom_after_mmsq_applied->SetYTitle("#theta_{#pi^{-}} rec (deg)");
        // pim_theta_rec_vs_mom_after_mmsq_applied->SetOption("COLZ1");
        // pim_theta_rec_vs_mom_after_mmsq_applied->Write();
        //
        // pim_theta_measured_vs_mom_after_mmsq_applied->SetTitle("#theta_{#pi^{-}} measured vs momentum");
        // pim_theta_measured_vs_mom_after_mmsq_applied->SetXTitle("P_{#pi^{-}} measured (GeV)");
        // pim_theta_measured_vs_mom_after_mmsq_applied->SetYTitle("#theta_{#pi^{-}} measured (deg)");
        // pim_theta_measured_vs_mom_after_mmsq_applied->SetOption("COLZ1");
        // pim_theta_measured_vs_mom_after_mmsq_applied->Write();
        //
        // phi_pim_rec_after_mmsq_applied->SetTitle("#phi_{#pi^{-}} rec ");
        // phi_pim_rec_after_mmsq_applied->SetXTitle("#phi_{#pi^{-}} rec (deg)");
        // phi_pim_rec_after_mmsq_applied->Write();
        // phi_pim_measured_after_mmsq_applied->SetTitle("#phi_{#pi^{-}} measured (deg)");
        // phi_pim_measured_after_mmsq_applied->SetXTitle("#phi_{#pi^{-}} measured (deg)");
        // phi_pim_measured_after_mmsq_applied->Write();
        //
        // pim_phi_measured_vs_rec_after_mmsq_applied->SetTitle("#phi_{#pi^{-}} measured vs #phi_{#pi^{-}} rec ");
        // pim_phi_measured_vs_rec_after_mmsq_applied->SetXTitle("#phi_{#pi^{-}} rec (deg)");
        // pim_phi_measured_vs_rec_after_mmsq_applied->SetYTitle("#phi_{#pi^{-}} measured (deg)");
        // pim_phi_measured_vs_rec_after_mmsq_applied->SetOption("COLZ1");
        // pim_phi_measured_vs_rec_after_mmsq_applied->Write();
        //
        // pim_phi_measured_vs_rec_after_y_equal_mc_cuts->SetTitle("#phi_{#pi^{-}} measured vs #phi_{#pi^{-}} rec ");
        // pim_phi_measured_vs_rec_after_y_equal_mc_cuts->SetXTitle("#phi_{#pi^{-}} rec (deg)");
        // pim_phi_measured_vs_rec_after_y_equal_mc_cuts->SetYTitle("#phi_{#pi^{-}} measured (deg)");
        // pim_phi_measured_vs_rec_after_y_equal_mc_cuts->SetOption("COLZ1");
        // pim_phi_measured_vs_rec_after_y_equal_mc_cuts->Write();
        //
        // pim_phi_rec_vs_mom_after_mmsq_applied->SetTitle("#phi_{#pi^{-}} rec vs P_{#pi^{-}}");
        // pim_phi_rec_vs_mom_after_mmsq_applied->SetXTitle("P_{#pi^{-}} rec (GeV)");
        // pim_phi_rec_vs_mom_after_mmsq_applied->SetYTitle("#phi_{#pi^{-}} rec (deg)");
        // pim_phi_rec_vs_mom_after_mmsq_applied->SetOption("COLZ1");
        // pim_phi_rec_vs_mom_after_mmsq_applied->Write();
        //
        // pim_phi_measured_vs_mom_after_mmsq_applied->SetTitle("#phi_{#pi^{-}} measured vs P_{#pi^{-}}");
        // pim_phi_measured_vs_mom_after_mmsq_applied->SetXTitle("P_{#pi^{-}} measured (GeV)");
        // pim_phi_measured_vs_mom_after_mmsq_applied->SetYTitle("#phi_{#pi^{-}} measured (deg)");
        // pim_phi_measured_vs_mom_after_mmsq_applied->SetOption("COLZ1");
        // pim_phi_measured_vs_mom_after_mmsq_applied->Write();
        //
        // pim_theta_vs_phi_rec_after_mmsq_applied->SetTitle("#theta_{#pi^{-}} vs #phi_{#pi^{-}} rec");
        // pim_theta_vs_phi_rec_after_mmsq_applied->SetXTitle("#phi_{#pi^{-}} rec (GeV)");
        // pim_theta_vs_phi_rec_after_mmsq_applied->SetYTitle("#theta_{#pi^{-}} rec (deg)");
        // pim_theta_vs_phi_rec_after_mmsq_applied->SetOption("COLZ1");
        // pim_theta_vs_phi_rec->Write();
        //
        // pim_theta_vs_phi_measured_after_mmsq_applied->SetTitle("#theta_{#pi^{-}} vs #phi_{#pi^{-}} measured");
        // pim_theta_vs_phi_measured_after_mmsq_applied->SetXTitle("#phi_{#pi^{-}} measured (GeV)");
        // pim_theta_vs_phi_measured_after_mmsq_applied->SetYTitle("#theta_{#pi^{-}} measured (deg)");
        // pim_theta_vs_phi_measured_after_mmsq_applied->SetOption("COLZ1");
        // pim_theta_vs_phi_measured_after_mmsq_applied->Write();


        auto theta_vs_mom = RootOutputFile->mkdir("theta_vs_mom");
        theta_vs_mom->cd();
        Theta_prot_cm_vs_mom_prot->SetXTitle("mom_prot (GeV)");
        Theta_prot_cm_vs_mom_prot->SetYTitle("theta_prot (Deg)");
        Theta_prot_cm_vs_mom_prot->SetOption("COLZ1");
        Theta_prot_cm_vs_mom_prot->Write("");
        Theta_pip_cm_vs_mom_pip->SetXTitle("mom_pip (GeV)");
        Theta_pip_cm_vs_mom_pip->SetYTitle("theta_pip (Deg)");
        Theta_pip_cm_vs_mom_pip->SetOption("COLZ1");
        Theta_pip_cm_vs_mom_pip->Write("");
        Theta_pim_cm_vs_mom_pim->SetXTitle("mom_pim (GeV)");
        Theta_pim_cm_vs_mom_pim->SetYTitle("theta_pim (Deg)");
        Theta_pim_cm_vs_mom_pim->SetOption("COLZ1");
        Theta_pim_cm_vs_mom_pim->Write("");

        Theta_prot_lab_vs_mom_prot->SetXTitle("mom_prot (GeV)");
        Theta_prot_lab_vs_mom_prot->SetYTitle("theta_prot (Deg)");
        Theta_prot_lab_vs_mom_prot->SetOption("COLZ1");
        Theta_prot_lab_vs_mom_prot->Write("");
        Theta_pip_lab_vs_mom_pip->SetXTitle("mom_pip (GeV)");
        Theta_pip_lab_vs_mom_pip->SetYTitle("theta_pip (Deg)");
        Theta_pip_lab_vs_mom_pip->SetOption("COLZ1");
        Theta_pip_lab_vs_mom_pip->Write("");
        Theta_pim_lab_vs_mom_pim->SetXTitle("mom_pim (GeV)");
        Theta_pim_lab_vs_mom_pim->SetYTitle("theta_pim (Deg)");
        Theta_pim_lab_vs_mom_pim->SetOption("COLZ1");
        Theta_pim_lab_vs_mom_pim->Write("");

        auto theta_vs_mom_thrown = RootOutputFile->mkdir("theta_vs_mom_thrown");
        theta_vs_mom_thrown->cd();
        Theta_prot_thrown_cm_vs_mom_prot->SetXTitle("mom_prot (GeV)");
        Theta_prot_thrown_cm_vs_mom_prot->SetYTitle("theta_prot (Deg)");
        Theta_prot_thrown_cm_vs_mom_prot->SetOption("COLZ1");
        Theta_prot_thrown_cm_vs_mom_prot->Write("");
        Theta_pip_thrown_cm_vs_mom_pip->SetXTitle("mom_pip (GeV)");
        Theta_pip_thrown_cm_vs_mom_pip->SetYTitle("theta_pip (Deg)");
        Theta_pip_thrown_cm_vs_mom_pip->SetOption("COLZ1");
        Theta_pip_thrown_cm_vs_mom_pip->Write("");
        Theta_pim_thrown_cm_vs_mom_pim->SetXTitle("mom_pim (GeV)");
        Theta_pim_thrown_cm_vs_mom_pim->SetYTitle("theta_pim (Deg)");
        Theta_pim_thrown_cm_vs_mom_pim->SetOption("COLZ1");
        Theta_pim_thrown_cm_vs_mom_pim->Write("");

        Theta_prot_thrown_lab_vs_mom_prot->SetXTitle("mom_prot (GeV)");
        Theta_prot_thrown_lab_vs_mom_prot->SetYTitle("theta_prot (Deg)");
        Theta_prot_thrown_lab_vs_mom_prot->SetOption("COLZ1");
        Theta_prot_thrown_lab_vs_mom_prot->Write("");
        Theta_pip_thrown_lab_vs_mom_pip->SetXTitle("mom_pip (GeV)");
        Theta_pip_thrown_lab_vs_mom_pip->SetYTitle("theta_pip (Deg)");
        Theta_pip_thrown_lab_vs_mom_pip->SetOption("COLZ1");
        Theta_pip_thrown_lab_vs_mom_pip->Write("");
        Theta_pim_thrown_lab_vs_mom_pim->SetXTitle("mom_pim (GeV)");
        Theta_pim_thrown_lab_vs_mom_pim->SetYTitle("theta_pim (Deg)");
        Theta_pim_thrown_lab_vs_mom_pim->SetOption("COLZ1");
        Theta_pim_thrown_lab_vs_mom_pim->Write("");

        auto Angles = RootOutputFile->mkdir("Angles");
        Angles->cd();

        theta_prot->SetXTitle("theta_prot (Deg)");
        if (theta_prot->GetEntries())
                theta_prot->Write();
        theta_pip->SetXTitle("theta_pip (Deg)");
        if (theta_pip->GetEntries())
                theta_pip->Write();
        theta_pim->SetXTitle("theta_pim (Deg)");
        if (theta_pim->GetEntries())
                theta_pim->Write();

        Phi_gamma->SetXTitle("Phi_gamma (Deg)");
        if (Phi_gamma->GetEntries())
                Phi_gamma->Write();
        Phi_prot->SetXTitle("Phi_prot (Deg)");
        if (Phi_prot->GetEntries())
                Phi_prot->Write();
        Phi_pip->SetXTitle("Phi_pip (Deg)");
        if (Phi_pip->GetEntries())
                Phi_pip->Write();
        Phi_pim->SetXTitle("Phi_pim (Deg)");
        if (Phi_pim->GetEntries())
                Phi_pim->Write();

        alpha_prot->SetXTitle("alpha_prot (Deg)");
        if (alpha_prot->GetEntries())
                alpha_prot->Write();
        alpha_pip->SetXTitle("alpha_pip (Deg)");
        if (alpha_pip->GetEntries())
                alpha_pip->Write();
        alpha_pim->SetXTitle("alpha_pim (Deg)");
        if (alpha_pim->GetEntries())
                alpha_pim->Write();
        //
        // theta_prot_mc->SetXTitle("theta_prot_mc (Deg)");
        // if (theta_prot_mc->GetEntries()) theta_prot_mc->Write();
        // theta_pip_mc->SetXTitle("theta_pip_mc (Deg)");
        // if (theta_pip_mc->GetEntries()) theta_pip_mc->Write();
        // theta_pim_mc->SetXTitle("theta_pim_mc (Deg)");
        // if (theta_pim_mc->GetEntries()) theta_pim_mc->Write();
        //
        // Phi_gamma_mc->SetXTitle("Phi_gamma_mc (Deg)");
        // if (Phi_gamma_mc->GetEntries()) Phi_gamma_mc->Write();
        // Phi_prot_mc->SetXTitle("Phi_prot_mc (Deg)");
        // if (Phi_prot_mc->GetEntries()) Phi_prot_mc->Write();
        // Phi_pip_mc->SetXTitle("Phi_pip_mc (Deg)");
        // if (Phi_pip_mc->GetEntries()) Phi_pip_mc->Write();
        // Phi_pim_mc->SetXTitle("Phi_pim_mc (Deg)");
        // if (Phi_pim_mc->GetEntries()) Phi_pim_mc->Write();
        //
        // alpha_prot_mc->SetXTitle("alpha_prot_mc (Deg)");
        // if (alpha_prot_mc->GetEntries()) alpha_prot_mc->Write();
        // alpha_pip_mc->SetXTitle("alpha_pip_mc (Deg)");
        // if (alpha_pip_mc->GetEntries()) alpha_pip_mc->Write();
        // alpha_pim_mc->SetXTitle("alpha_pim_mc (Deg)");
        // if (alpha_pim_mc->GetEntries()) alpha_pim_mc->Write();
        auto Angles_thrown = RootOutputFile->mkdir("Angles_thrown");
        Angles_thrown->cd();
        theta_prot_thrown->SetXTitle("theta_prot_thrown (Deg)");
        if (theta_prot_thrown->GetEntries())
                theta_prot_thrown->Write();
        theta_pip_thrown->SetXTitle("theta_pip_thrown (Deg)");
        if (theta_pip_thrown->GetEntries())
                theta_pip_thrown->Write();
        theta_pim_thrown->SetXTitle("theta_pim_thrown (Deg)");
        if (theta_pim_thrown->GetEntries())
                theta_pim_thrown->Write();

        Phi_gamma_thrown->SetXTitle("Phi_gamma_thrown (Deg)");
        if (Phi_gamma_thrown->GetEntries())
                Phi_gamma_thrown->Write();
        Phi_prot_thrown->SetXTitle("Phi_prot_thrown (Deg)");
        if (Phi_prot_thrown->GetEntries())
                Phi_prot_thrown->Write();
        Phi_pip_thrown->SetXTitle("Phi_pip_thrown (Deg)");
        if (Phi_pip_thrown->GetEntries())
                Phi_pip_thrown->Write();
        Phi_pim_thrown->SetXTitle("Phi_pim_thrown (Deg)");
        if (Phi_pim_thrown->GetEntries())
                Phi_pim_thrown->Write();

        alpha_prot_thrown->SetXTitle("alpha_prot_thrown (Deg)");
        if (alpha_prot_thrown->GetEntries())
                alpha_prot_thrown->Write();
        alpha_pip_thrown->SetXTitle("alpha_pip_thrown (Deg)");
        if (alpha_pip_thrown->GetEntries())
                alpha_pip_thrown->Write();
        alpha_pim_thrown->SetXTitle("alpha_pim_thrown (Deg)");
        if (alpha_pim_thrown->GetEntries())
                alpha_pim_thrown->Write();
        //
        // Theta_vs_mom_x_mu->SetYTitle("Theta_x_mu");
        // Theta_vs_mom_x_mu->SetXTitle("mom_x_mu");
        // Theta_vs_mom_x_mu->SetOption("COLZ1");
        // Theta_vs_mom_x_mu->Write();

        // auto singlePi_sec = RootOutputFile->mkdir("singlePi_sec");
        // singlePi_sec->cd();
        // for (short i = 0; i < num_sectors; i++) {
        //   W_vs_q2_singlePi_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
        //   W_vs_q2_singlePi_sec[i]->SetXTitle("W (GeV)");
        //   W_vs_q2_singlePi_sec[i]->SetOption("COLZ1");
        //   if (W_vs_q2_singlePi_sec[i]->GetEntries())
        //   W_vs_q2_singlePi_sec[i]->Write();
        // }
        //
        // for (short i = 0; i < num_sectors; i++) {
        //   W_vs_MM_singlePi[i]->SetOption("COLZ1");
        //   W_vs_MM_singlePi[i]->SetYTitle("MM (GeV)");
        //   W_vs_MM_singlePi[i]->SetXTitle("W (GeV)");
        //   if (W_vs_MM_singlePi[i]->GetEntries()) W_vs_MM_singlePi[i]->Write();
        // }
        //
        // for (short i = 0; i < num_sectors; i++) {
        //   W_singlePi_sec[i]->SetXTitle("W (GeV)");
        //   if (W_singlePi_sec[i]->GetEntries()) W_singlePi_sec[i]->Write();
        // }
        //
        // for (short i = 0; i < num_sectors; i++) {
        //   //    if (MM_neutron_sec[i]->GetEntries()) MM_neutron_sec[i]->Fit("gaus",
        //   "QMR+", "QMR+", 0.7, 1.1); MM_neutron_sec[i]->SetXTitle("Mass (GeV)"); if
        //   (MM_neutron_sec[i]->GetEntries()) MM_neutron_sec[i]->Write();
        // }

        auto twoPi_sec = RootOutputFile->mkdir("twoPi_sec");
        twoPi_sec->cd();
        for (short i = 0; i < num_sectors; i ++) {
                W_vs_q2_twoPi_sec[i] -> SetYTitle("Q^{2} (GeV^{2})");
                W_vs_q2_twoPi_sec[i] -> SetXTitle("W (GeV)");
                W_vs_q2_twoPi_sec[i] -> SetOption("COLZ1");
                if (W_vs_q2_twoPi_sec[i] -> GetEntries())
                        W_vs_q2_twoPi_sec[i] -> Write();
        }

        for (short i = 0; i < num_sectors; i ++) {
                W_vs_MM_twoPi[i] -> SetOption("COLZ1");
                W_vs_MM_twoPi[i] -> SetYTitle("MM (GeV)");
                W_vs_MM_twoPi[i] -> SetXTitle("W (GeV)");
                if (W_vs_MM_twoPi[i] -> GetEntries())
                        W_vs_MM_twoPi[i] -> Write();
        }

        for (short i = 0; i < num_sectors; i ++) {
                W_twoPi_sec[i] -> SetXTitle("W (GeV)");
                if (W_twoPi_sec[i] -> GetEntries())
                        W_twoPi_sec[i] -> Write();
        }

        for (short i = 0; i < num_sectors; i ++) {
                //  if (MM_twoPi_sec[i]->GetEntries()) MM_twoPi_sec[i]->Fit("gaus", "QMR+",
                //  "QMR+", -0.1, 0.1);
                MM_twoPi_sec[i] -> SetXTitle("MisingMass (GeV)");
                MM_twoPi_sec[i] -> Write();
        }

        for (short i = 0; i < num_sectors; i ++) {
                // if (MM2_twoPi_sec[i]->GetEntries()) MM2_twoPi_sec[i]->Fit("gaus", "QMR+",
                // "QMR+", -0.1, 0.1);
                MM2_twoPi_sec[i] -> SetXTitle("MM2 (GeV2)");
                MM2_twoPi_sec[i] -> Write();
        }

        for (short i = 0; i < num_sectors; i ++) {
                MM2_twoPi_missingPip_sec[i] -> SetXTitle("MM2 (GeV2)");
                MM2_twoPi_missingPip_sec[i] -> Write();
        }

        for (short i = 0; i < num_sectors; i ++) {
                MM2_twoPi_missingProt_sec[i] -> SetXTitle("MM2 (GeV2)");
                MM2_twoPi_missingProt_sec[i] -> Write();
        }
        auto pi0_mass_hist = RootOutputFile->mkdir("pi0_mass_hist");
        pi0_mass_hist->cd();
        for (short i = 0; i < num_sectors; i ++) {
                mass_pi0_any_event[i] -> SetXTitle("pi0 mass (GeV)");
                mass_pi0_any_event[i] -> Write();
        }
        for (short i = 0; i < num_sectors; i ++) {
                mass_pi0_hist_before_mmsq_cut[i] -> SetXTitle("pi0 mass (GeV)");
                mass_pi0_hist_before_mmsq_cut[i] -> Write();
        }
        for (short i = 0; i < num_sectors; i ++) {
                mass_pi0_hist_after_mmsq_cut[i] -> SetXTitle("pi0 mass (GeV)");
                mass_pi0_hist_after_mmsq_cut[i] -> Write();
        }

        auto theta_vs_mom_hist = RootOutputFile->mkdir("theta_vs_mom_hist");
        theta_vs_mom_hist->cd();
        for (short i = 0; i < num_sectors; i ++) {
                theta_vs_mom_elec[i] -> SetTitle("#theta_{e} vs mom");
                theta_vs_mom_elec[i] -> SetXTitle("P (GeV)");
                theta_vs_mom_elec[i] -> SetYTitle("#theta_{e} (deg)");
                theta_vs_mom_elec[i] -> SetOption("COLZ1");
                theta_vs_mom_elec[i] -> Write();
        }
        for (short i = 0; i < num_sectors; i ++) {
                theta_vs_mom_prot[i] -> SetTitle("#theta_{p} vs mom");
                theta_vs_mom_prot[i] -> SetXTitle("P (GeV)");
                theta_vs_mom_prot[i] -> SetYTitle("#theta_{p} (deg)");
                theta_vs_mom_prot[i] -> SetOption("COLZ1");
                theta_vs_mom_prot[i] -> Write();
        }
        for (short i = 0; i < num_sectors; i ++) {
                theta_vs_mom_pip[i] -> SetTitle("#theta_{#pi^{+}} vs mom");
                theta_vs_mom_pip[i] -> SetXTitle("P (GeV)");
                theta_vs_mom_pip[i] -> SetYTitle("#theta_{#pi^{+}} (deg)");
                theta_vs_mom_pip[i] -> SetOption("COLZ1");
                theta_vs_mom_pip[i] -> Write();
        }
        for (short i = 0; i < num_sectors; i ++) {
                theta_vs_mom_pim[i] -> SetTitle("#theta_{#pi^{-}} vs mom");
                theta_vs_mom_pim[i] -> SetXTitle("P (GeV)");
                theta_vs_mom_pim[i] -> SetYTitle("#theta_{#pi^{-}} (deg)");
                theta_vs_mom_pim[i] -> SetOption("COLZ1");
                theta_vs_mom_pim[i] -> Write();
        }

        auto Inv_Mass_hist = RootOutputFile->mkdir("Inv_Mass_hist");
        Inv_Mass_hist->cd();
        for (short i = 0; i < w_range_num; i ++) {
                for (short j = 0; j < q2_range_num; j ++) {
                        inv_mass_P_pip[i][j] -> SetXTitle("inv_mass (GeV)");
                        inv_mass_P_pip[i][j] -> Write();
                        inv_mass_P_pim[i][j] -> SetXTitle("inv_mass (GeV)");
                        inv_mass_P_pim[i][j] -> Write();
                        inv_mass_pip_pim[i][j] -> SetXTitle("inv_mass (GeV)");
                        inv_mass_pip_pim[i][j] -> Write();

                        theta_P_vs_mass_pip_pim[i][j] -> SetYTitle("theta_P_cm");
                        theta_P_vs_mass_pip_pim[i][j] -> SetXTitle("inv_mass_pip_pim (GeV)");
                        theta_P_vs_mass_pip_pim[i][j] -> SetOption("COLZ1");
                        theta_P_vs_mass_pip_pim[i][j] -> Write();

                        theta_pip_vs_mass_Ppim[i][j] -> SetYTitle("theta_pip_cm");
                        theta_pip_vs_mass_Ppim[i][j] -> SetXTitle("inv_mass_Ppim (GeV)");
                        theta_pip_vs_mass_Ppim[i][j] -> SetOption("COLZ1");
                        theta_pip_vs_mass_Ppim[i][j] -> Write();

                        theta_pim_vs_mass_Ppip[i][j] -> SetYTitle("theta_pim_cm");
                        theta_pim_vs_mass_Ppip[i][j] -> SetXTitle("inv_mass_Ppip (GeV)");
                        theta_pim_vs_mass_Ppip[i][j] -> SetOption("COLZ1");
                        theta_pim_vs_mass_Ppip[i][j] -> Write();

                        theta_P_lab_vs_mass_pip_pim[i][j] -> SetYTitle("theta_P_cm");
                        theta_P_lab_vs_mass_pip_pim[i][j] -> SetXTitle("inv_mass_pip_pim (GeV)");
                        theta_P_lab_vs_mass_pip_pim[i][j] -> SetOption("COLZ1");
                        theta_P_lab_vs_mass_pip_pim[i][j] -> Write();

                        theta_pip_lab_vs_mass_Ppim[i][j] -> SetYTitle("theta_pip_cm");
                        theta_pip_lab_vs_mass_Ppim[i][j] -> SetXTitle("inv_mass_Ppim (GeV)");
                        theta_pip_lab_vs_mass_Ppim[i][j] -> SetOption("COLZ1");
                        theta_pip_lab_vs_mass_Ppim[i][j] -> Write();

                        theta_pim_lab_vs_mass_Ppip[i][j] -> SetYTitle("theta_pim_cm");
                        theta_pim_lab_vs_mass_Ppip[i][j] -> SetXTitle("inv_mass_Ppip (GeV)");
                        theta_pim_lab_vs_mass_Ppip[i][j] -> SetOption("COLZ1");
                        theta_pim_lab_vs_mass_Ppip[i][j] -> Write();
                }
        }
        auto twoPi_sec_thrown = RootOutputFile->mkdir("twoPi_sec_thrown");
        twoPi_sec_thrown->cd();

        for (short i = 0; i < num_sectors; i ++) {
                W_vs_q2_twoPi_sec_thrown[i] -> SetYTitle("Q^{2} (GeV^{2})");
                W_vs_q2_twoPi_sec_thrown[i] -> SetXTitle("W (GeV)");
                W_vs_q2_twoPi_sec_thrown[i] -> SetOption("COLZ1");
                if (W_vs_q2_twoPi_sec_thrown[i] -> GetEntries())
                        W_vs_q2_twoPi_sec_thrown[i] -> Write();
        }

        for (short i = 0; i < num_sectors; i ++) {
                W_vs_MM_twoPi_thrown[i] -> SetOption("COLZ1");
                W_vs_MM_twoPi_thrown[i] -> SetYTitle("MM (GeV)");
                W_vs_MM_twoPi_thrown[i] -> SetXTitle("W (GeV)");
                if (W_vs_MM_twoPi[i] -> GetEntries()) W_vs_MM_twoPi_thrown[i] -> Write();
        }
        //
        // for (short i = 0; i < num_sectors; i++) {
        //   W_twoPi_sec_thrown[i]->SetXTitle("W (GeV)");
        //   if (W_twoPi_sec_thrown[i]->GetEntries()) W_twoPi_sec_thrown[i]->Write();
        // }
        //
        // for (short i = 0; i < num_sectors; i++) {
        //   //  if (MM_twoPi_sec_thrown[i]->GetEntries())
        //   MM_twoPi_sec_thrown[i]->Fit("gaus", "QMR+", "QMR+", -0.1, 0.1);
        //   MM_twoPi_sec_thrown[i]->SetXTitle("MisingMass (GeV)");
        //   MM_twoPi_sec_thrown[i]->Write();
        // }
        //
        // for (short i = 0; i < num_sectors; i++) {
        //   // if (MM2_twoPi_sec_thrown[i]->GetEntries())
        //   MM2_twoPi_sec_thrown[i]->Fit("gaus", "QMR+", "QMR+", -0.1, 0.1);
        //   MM2_twoPi_sec_thrown[i]->SetXTitle("MisingMass (GeV)");
        //   MM2_twoPi_sec_thrown[i]->Write();
        // }
        //
        // auto Inv_Mass_hist_thrown = RootOutputFile->mkdir("Inv_Mass_hist_thrown");
        // Inv_Mass_hist_thrown->cd();
        // for (short i = 0; i < w_range_num; i++) {
        //   for (short j = 0; j < q2_range_num; j++) {
        //     inv_mass_P_pip_thrown[i][j]->SetXTitle("inv_mass (GeV)");
        //     inv_mass_P_pip_thrown[i][j]->Write();
        //     inv_mass_P_pim_thrown[i][j]->SetXTitle("inv_mass (GeV)");
        //     inv_mass_P_pim_thrown[i][j]->Write();
        //     inv_mass_pip_pim_thrown[i][j]->SetXTitle("inv_mass (GeV)");
        //     inv_mass_pip_pim_thrown[i][j]->Write();
        //
        //     theta_P_vs_mass_pip_pim_thrown[i][j]->SetYTitle("theta_P_cm");
        //     theta_P_vs_mass_pip_pim_thrown[i][j]->SetXTitle("inv_mass_pip_pim
        //     (GeV)"); theta_P_vs_mass_pip_pim_thrown[i][j]->SetOption("COLZ1");
        //     theta_P_vs_mass_pip_pim_thrown[i][j]->Write();
        //
        //     theta_pip_vs_mass_Ppim_thrown[i][j]->SetYTitle("theta_pip_cm");
        //     theta_pip_vs_mass_Ppim_thrown[i][j]->SetXTitle("inv_mass_Ppim (GeV)");
        //     theta_pip_vs_mass_Ppim_thrown[i][j]->SetOption("COLZ1");
        //     theta_pip_vs_mass_Ppim_thrown[i][j]->Write();
        //
        //     theta_pim_vs_mass_Ppip_thrown[i][j]->SetYTitle("theta_pim_cm");
        //     theta_pim_vs_mass_Ppip_thrown[i][j]->SetXTitle("inv_mass_Ppip (GeV)");
        //     theta_pim_vs_mass_Ppip_thrown[i][j]->SetOption("COLZ1");
        //     theta_pim_vs_mass_Ppip_thrown[i][j]->Write();
        //
        //     theta_P_lab_vs_mass_pip_pim_thrown[i][j]->SetYTitle("theta_P_cm");
        //     theta_P_lab_vs_mass_pip_pim_thrown[i][j]->SetXTitle("inv_mass_pip_pim
        //     (GeV)"); theta_P_lab_vs_mass_pip_pim_thrown[i][j]->SetOption("COLZ1");
        //     theta_P_lab_vs_mass_pip_pim_thrown[i][j]->Write();
        //
        //     theta_pip_lab_vs_mass_Ppim_thrown[i][j]->SetYTitle("theta_pip_cm");
        //     theta_pip_lab_vs_mass_Ppim_thrown[i][j]->SetXTitle("inv_mass_Ppim
        //     (GeV)"); theta_pip_lab_vs_mass_Ppim_thrown[i][j]->SetOption("COLZ1");
        //     theta_pip_lab_vs_mass_Ppim_thrown[i][j]->Write();
        //
        //     theta_pim_lab_vs_mass_Ppip_thrown[i][j]->SetYTitle("theta_pim_cm");
        //     theta_pim_lab_vs_mass_Ppip_thrown[i][j]->SetXTitle("inv_mass_Ppip
        //     (GeV)"); theta_pim_lab_vs_mass_Ppip_thrown[i][j]->SetOption("COLZ1");
        //     theta_pim_lab_vs_mass_Ppip_thrown[i][j]->Write();
        //   }
        // }
        // auto Npip_sec = RootOutputFile->mkdir("Npip_sec");
        // Npip_sec->cd();
        //
        // for (short i = 0; i < num_sectors; i++) {
        //   W_vs_q2_Npip_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
        //   W_vs_q2_Npip_sec[i]->SetXTitle("W (GeV)");
        //   W_vs_q2_Npip_sec[i]->SetOption("COLZ1");
        //   W_vs_q2_Npip_sec[i]->Write();
        // }
        // for (short i = 0; i < num_sectors; i++) {
        //   W_Npip_sec[i]->SetXTitle("W (GeV)");
        //   W_Npip_sec[i]->Write();
        // }
        // for (short i = 0; i < num_sectors; i++) {
        //   MM_Npip_sec[i]->SetXTitle("Mass (GeV)");
        //   MM_Npip_sec[i]->Write();
        // }
}
// void Histogram::makeHists_electron_cuts() {
//   for (auto&& cut : before_after_cut) {
//     int c = cut.first;
//     auto type = cut.second.c_str();
//     EC_sampling_fraction[c] =
//         std::make_shared<TH2D>(Form("EC_sampling_fraction%s", type),
//         Form("EC_sampling_fraction%s", type), bins, p_min,
//                                p_max, bins, zero, 0.5);
//     vz_position[c] = std::make_shared<TH1D>(Form("vz_position%s", type),
//     Form("vz_position%s", type), bins, -40, 40); pcal_sec[c] =
//         std::make_shared<TH2D>(Form("pcal_sec%s", type), Form("pcal_sec%s",
//         type), bins, -420, 420, bins, -420, 420);
//     dcr1_sec[c] =
//         std::make_shared<TH2D>(Form("dcr1_sec%s", type), Form("dcr1_sec%s",
//         type), bins, -180, 180, bins, -180, 180);
//     dcr2_sec[c] =
//         std::make_shared<TH2D>(Form("dcr2_sec%s", type), Form("dcr2_sec%s",
//         type), bins, -270, 270, bins, -270, 270);
//     dcr3_sec[c] =
//         std::make_shared<TH2D>(Form("dcr3_sec%s", type), Form("dcr3_sec%s",
//         type), bins, -320, 320, bins, -320, 320);
//   }
// }
// void Histogram::FillHists_electron_cuts(const std::shared_ptr<Branches12>&
// _d) {
//   vz_position[before_cut]->Fill(_d->vz(0));
//   pcal_sec[before_cut]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
//   dcr1_sec[before_cut]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0));
//   dcr2_sec[before_cut]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0));
//   dcr3_sec[before_cut]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0));
//   EC_sampling_fraction[before_cut]->Fill(_d->p(0), _d->ec_tot_energy(0) /
//   _d->p(0));
// }
//
// void Histogram::FillHists_electron_with_cuts(const
// std::shared_ptr<Branches12>& _d) {
//   vz_position[after_cut]->Fill(_d->vz(0));
//   pcal_sec[after_cut]->Fill(_d->ec_pcal_x(0), _d->ec_pcal_y(0));
//   dcr1_sec[after_cut]->Fill(_d->dc_r1_x(0), _d->dc_r1_y(0));
//   dcr2_sec[after_cut]->Fill(_d->dc_r2_x(0), _d->dc_r2_y(0));
//   dcr3_sec[after_cut]->Fill(_d->dc_r3_x(0), _d->dc_r3_y(0));
//   EC_sampling_fraction[after_cut]->Fill(_d->p(0), _d->ec_tot_energy(0) /
//   _d->p(0));
// }
//
// void Histogram::Write_Electron_cuts() {
//   for (auto&& cut : before_after_cut) {
//     int c = cut.first;
//     vz_position[c]->SetXTitle("vz (GeV)");
//     vz_position[c]->Fit("gaus", "QMR+", "QMR+", -7.089, 2.0);
//     gStyle->SetOptFit(1111);
//     if (vz_position[c]->GetEntries()) vz_position[c]->Write();
//     pcal_sec[c]->SetXTitle("x (cm)");
//     pcal_sec[c]->SetYTitle("y (cm)");
//     pcal_sec[c]->SetOption("COLZ1");
//     if (pcal_sec[c]->GetEntries()) pcal_sec[c]->Write();
//
//     dcr1_sec[c]->SetXTitle("x (cm)");
//     dcr1_sec[c]->SetYTitle("y (cm)");
//     dcr1_sec[c]->SetOption("COLZ1");
//     if (dcr1_sec[c]->GetEntries()) dcr1_sec[c]->Write();
//
//     dcr2_sec[c]->SetXTitle("x (cm)");
//     dcr2_sec[c]->SetYTitle("y (cm)");
//     dcr2_sec[c]->SetOption("COLZ1");
//     if (dcr2_sec[c]->GetEntries()) dcr2_sec[c]->Write();
//
//     dcr3_sec[c]->SetXTitle("x (cm)");
//     dcr3_sec[c]->SetYTitle("y (cm)");
//     dcr3_sec[c]->SetOption("COLZ1");
//     if (dcr3_sec[c]->GetEntries()) dcr3_sec[c]->Write();
//
//     EC_sampling_fraction[c]->SetXTitle("Momentum (GeV)");
//     EC_sampling_fraction[c]->SetYTitle("Sampling Fraction");
//     EC_sampling_fraction[c]->SetOption("COLZ1");
//     EC_sampling_fraction[c]->Write();
//   }
// }
void Histogram::makeHists_sector() {
        for (short i = 0; i < w_range_num; i ++) {
                for (short j = 0; j < q2_range_num; j ++) {
                        inv_mass_P_pip[i][j] = std::make_shared<TH1D>(
                                Form("inv_mass_P_pip_%6.12s_%6.12s", w_range_name[i].c_str(),
                                     q2_range_name[j].c_str()),
                                Form("inv_mass_P_pip: %6.12s_%6.12s", w_range_name[i].c_str(),
                                     q2_range_name[j].c_str()),
                                bins, 0., 3.5);
                        inv_mass_P_pim[i][j] = std::make_shared<TH1D>(
                                Form("inv_mass_P_pim_%6.12s_%6.12s", w_range_name[i].c_str(),
                                     q2_range_name[j].c_str()),
                                Form("inv_mass_P_pim: %6.12s_%6.12s", w_range_name[i].c_str(),
                                     q2_range_name[j].c_str()),
                                bins, 0., 3.5);
                        inv_mass_pip_pim[i][j] = std::make_shared<TH1D>(
                                Form("inv_mass_pip_pim_%6.12s_%6.12s", w_range_name[i].c_str(),
                                     q2_range_name[j].c_str()),
                                Form("inv_mass_pip_pim: %6.12s_%6.12s", w_range_name[i].c_str(),
                                     q2_range_name[j].c_str()),
                                bins, 0.0, 2.5);
                        theta_P_vs_mass_pip_pim[i][j] = std::make_shared<TH2D>(
                                Form("theta_P_vs_inv_mass_pip_pim_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                Form("theta_P_vs_inv_mass_pip_pim_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                bins, 0.2, 3.5, bins, zero, 180);
                        theta_pim_vs_mass_Ppip[i][j] = std::make_shared<TH2D>(
                                Form("theta_pim_vs_inv_mass_Ppip_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                Form("theta_pim_vs_inv_mass_Ppip_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                bins, 1., 3.5, bins, zero, 180);
                        theta_pip_vs_mass_Ppim[i][j] = std::make_shared<TH2D>(
                                Form("theta_pip_vs_inv_mass_Ppim_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                Form("theta_pip_vs_inv_mass_Ppim_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                bins, 1., 3.5, bins, zero, 180);
                        theta_P_lab_vs_mass_pip_pim[i][j] = std::make_shared<TH2D>(
                                Form("theta_P_lab_vs_inv_mass_pip_pim_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                Form("theta_P_lab_vs_inv_mass_pip_pim_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                bins, 0.2, 2.5, bins, zero, 60);
                        theta_pim_lab_vs_mass_Ppip[i][j] = std::make_shared<TH2D>(
                                Form("theta_pim_lab_vs_inv_mass_Ppip_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                Form("theta_pim_lab_vs_inv_mass_Ppip_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                bins, 1., 3.5, bins, zero, 120);
                        theta_pip_lab_vs_mass_Ppim[i][j] = std::make_shared<TH2D>(
                                Form("theta_pip_lab_vs_inv_mass_Ppim_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                Form("theta_pip_lab_vs_inv_mass_Ppim_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                bins, 1., 3.5, bins, zero, 120);

                        inv_mass_P_pip_thrown[i][j] = std::make_shared<TH1D>(
                                Form("inv_mass_P_pip_thrown_%6.12s_%6.12s", w_range_name[i].c_str(),
                                     q2_range_name[j].c_str()),
                                Form("inv_mass_P_pip_thrown_: %6.12s_%6.12s", w_range_name[i].c_str(),
                                     q2_range_name[j].c_str()),
                                bins, 1., 3.5);
                        inv_mass_P_pim_thrown[i][j] = std::make_shared<TH1D>(
                                Form("inv_mass_P_pim_thrown_%6.12s_%6.12s", w_range_name[i].c_str(),
                                     q2_range_name[j].c_str()),
                                Form("inv_mass_P_pim_thrown_: %6.12s_%6.12s", w_range_name[i].c_str(),
                                     q2_range_name[j].c_str()),
                                bins, 1., 3.5);
                        inv_mass_pip_pim_thrown[i][j] = std::make_shared<TH1D>(
                                Form("inv_mass_pip_pim_thrown%6.12s_%6.12s", w_range_name[i].c_str(),
                                     q2_range_name[j].c_str()),
                                Form("inv_mass_pip_pim_thrown: %6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                bins, 0.2, 2.5);
                        theta_P_vs_mass_pip_pim_thrown[i][j] = std::make_shared<TH2D>(
                                Form("theta_P_vs_inv_mass_pip_pim_thrown_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                Form("theta_P_vs_inv_mass_pip_pim_thrown_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                bins, 0.2, 3.5, bins, zero, 180);
                        theta_pim_vs_mass_Ppip_thrown[i][j] = std::make_shared<TH2D>(
                                Form("theta_pim_vs_inv_mass_Ppip_thrown_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                Form("theta_pim_vs_inv_mass_Ppip_thrown_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                bins, 1., 3.5, bins, zero, 180);
                        theta_pip_vs_mass_Ppim_thrown[i][j] = std::make_shared<TH2D>(
                                Form("theta_pip_vs_inv_mass_Ppim_thrown_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                Form("theta_pip_vs_inv_mass_Ppim_thrown_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                bins, 1., 3.5, bins, zero, 180);
                        theta_P_lab_vs_mass_pip_pim_thrown[i][j] = std::make_shared<TH2D>(
                                Form("theta_P_lab_vs_inv_mass_pip_pim_thrown_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                Form("theta_P_lab_vs_inv_mass_pip_pim_thrown_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                bins, 0.2, 2.5, bins, zero, 70);
                        theta_pim_lab_vs_mass_Ppip_thrown[i][j] = std::make_shared<TH2D>(
                                Form("theta_pim_lab_vs_inv_mass_Ppip_thrown_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                Form("theta_pim_lab_vs_inv_mass_Ppip_thrown_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                bins, 1., 3.5, bins, zero, 120);
                        theta_pip_lab_vs_mass_Ppim_thrown[i][j] = std::make_shared<TH2D>(
                                Form("theta_pip_lab_vs_inv_mass_Ppim_thrown_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                Form("theta_pip_lab_vs_inv_mass_Ppim_thrown_%6.12s_%6.12s",
                                     w_range_name[i].c_str(), q2_range_name[j].c_str()),
                                bins, 1., 3.5, bins, zero, 120);
                }
        }
        for (short i = 0; i < 3; i ++) {
                W_det[i] = std::make_shared<TH1D>(Form("W_det_%d", i + 1),
                                                  Form("W detector: %d", i + 1), bins, zero,
                                                  w_max);
                if (i == 0)
                        WQ2_det[i] = std::make_shared<TH2D>(
                                Form("WQ2_det_%d", i + 1), Form("W vs Q^{2} detector: %d", i + 1),
                                bins, zero, w_max, bins, zero, 0.5);
                else
                        WQ2_det[i] = std::make_shared<TH2D>(
                                Form("WQ2_det_%d", i + 1), Form("W vs Q^{2} detector: %d", i + 1),
                                bins, zero, w_max, bins, zero, q2_max);
        }

        for (short i = 0; i < num_sectors; i ++) {
                W_vs_q2_sec[i] = std::make_shared<TH2D>(
                        Form("wvsq2_sec_%d", i + 1), Form("W vs Q^{2} Sector: %d", i + 1), bins,
                        zero, w_max, bins, zero, q2_max);

                W_sec[i] =
                        std::make_shared<TH1D>(Form("w_sec_%d", i + 1),
                                               Form("W Sector: %d", i + 1), bins, zero, w_max);

                W_vs_q2_twoPi_sec[i] =
                        std::make_shared<TH2D>(Form("wvsq2_sec_twoPi_%d", i + 1),
                                               Form("W vs Q^{2} W_twoPi Sector: %d", i + 1),
                                               bins, zero, w_max, bins, zero, q2_max);
                W_twoPi_sec[i] = std::make_shared<TH1D>(Form("w_sec_twoPi_%d", i + 1),
                                                        Form("W twoPi Sector: %d", i + 1),
                                                        bins, zero, w_max);

                W_vs_MM_twoPi[i] = std::make_shared<TH2D>(
                        Form("W_vs_MM_twoPi_%d", i + 1), Form("W_vs_MM_twoPi_%d", i + 1), bins,
                        zero, w_max, bins, -q2_max, q2_max);

                MM_twoPi_sec[i] = std::make_shared<TH1D>(
                        Form("MM_sec_%d", i + 1), Form("MM missingPim Sector: %d", i + 1), bins,
                        -0.4, 0.4);
                MM2_twoPi_sec[i] = std::make_shared<TH1D>(
                        Form("MM_SQ_sec_%d", i + 1), Form("MM_SQ missingPim Sector: %d", i + 1),
                        bins, -0.2, 0.2);
                MM2_twoPi_missingPip_sec[i] = std::make_shared<TH1D>(
                        Form("MM_SQ_missingPip_sec_%d", i + 1),
                        Form("MM_SQ missingPiP Sector: %d", i + 1), bins, -0.2, 0.2);
                MM2_twoPi_missingProt_sec[i] = std::make_shared<TH1D>(
                        Form("MM_SQt_missingProt_sec_%d", i + 1),
                        Form("MM_SQ missingProt Sector: %d", i + 1), bins, 0.6, 1.2);

                W_vs_MM_twoPi_thrown[i] =
                        std::make_shared<TH2D>(Form("W_vs_MM_missingPim_thrown_%d", i + 1),
                                               Form("W_vs_MM_missingPim_thrown_%d", i + 1),
                                               bins, zero, w_max, bins, -q2_max, q2_max);

                W_vs_q2_twoPi_sec_thrown[i] = std::make_shared<TH2D>(
                        Form("wvsq2_sec_missingPim_thrown_%d", i + 1),
                        Form("W vs Q^{2} W_missingPim Sector_thrown: %d", i + 1), bins, zero,
                        w_max, bins, zero, q2_max);
                W_twoPi_sec_thrown[i] = std::make_shared<TH1D>(
                        Form("w_sec_missingPim_thrown_%d", i + 1),
                        Form("W missingPim Sector_thrown: %d", i + 1), bins, zero, w_max);

                MM_twoPi_sec_thrown[i] = std::make_shared<TH1D>(
                        Form("MM_sec_thrown_%d", i + 1),
                        Form("MM twoPi Sector_thrown: %d", i + 1), bins, -0.4, 0.4);
                MM2_twoPi_sec_thrown[i] = std::make_shared<TH1D>(
                        Form("MM_SQ_sec_thrown_%d", i + 1),
                        Form("MM_SQ twoPi Sector_thrown: %d", i + 1), bins, -0.4, 0.4);

                W_vs_q2_singlePi_sec[i] =
                        std::make_shared<TH2D>(Form("wvsq2_sec_singlePi_%d", i + 1),
                                               Form("W vs Q^{2} W_singlePi Sector: %d", i + 1),
                                               bins, zero, w_max, bins, zero, q2_max);

                W_singlePi_sec[i] = std::make_shared<TH1D>(
                        Form("w_sec_singlePi_%d", i + 1), Form("W singlePi Sector: %d", i + 1),
                        bins, zero, w_max);

                W_vs_q2_Npip_sec[i] =
                        std::make_shared<TH2D>(Form("wvsq2_sec_Npip_%d", i + 1),
                                               Form("W vs Q^{2} W_Npip Sector: %d", i + 1),
                                               bins, zero, w_max, bins, zero, q2_max);

                W_Npip_sec[i] = std::make_shared<TH1D>(Form("w_sec_Npip_%d", i + 1),
                                                       Form("W Npip Sector: %d", i + 1),
                                                       bins, zero, w_max);

                MM_neutron_sec[i] = std::make_shared<TH1D>(
                        Form("MM_Sec_%d", i + 1), Form("MM neutron Sector: %d", i + 1), bins,
                        zero, 4.0);

                MM_Npip_sec[i] = std::make_shared<TH1D>(
                        Form("MM_Npip_Sec_%d", i + 1),
                        Form("MM^{2} neutron pip Sector: %d", i + 1), bins, -2.0, 2.0);

                W_vs_MM_singlePi[i] = std::make_shared<TH2D>(
                        Form("W_vs_MM_singlePi_%d", i + 1), Form("W_vs_MM_singlePi_%d", i + 1),
                        bins, zero, w_max, bins, -q2_max, q2_max);

                mass_pi0_any_event[i] = std::make_shared<TH1D>(
                        Form("mass_pi0_hist_any_event%d", i + 1),
                        Form("mass_pi0_hist_any_event_%d", i + 1), bins, 0, 1.1);
                mass_pi0_hist_before_mmsq_cut[i] = std::make_shared<TH1D>(
                        Form("mass_pi0_hist_2pi_before_MMSQ_cut%d", i + 1),
                        Form("mass_pi0_hist_TwoPi_before_MMSQ_cut_%d", i + 1), bins, 0, 1.1);
                mass_pi0_hist_after_mmsq_cut[i] = std::make_shared<TH1D>(
                        Form("mass_pi0_hist_2pi_after_MMSQ_cut_%d", i + 1),
                        Form("mass_pi0_TwoPi_hist_after_MMSQ_cuts_%d", i + 1), bins, 0, 1.1);

                // theta vs mom histograms added Nov Nov_2020
                theta_vs_mom_elec[i] = std::make_shared<TH2D>(
                        Form("theta_vs_mom_elec_%d", i + 1), Form("theta_vs_mom_elec_%d", i + 1), bins,
                        zero, p_max, bins, zero, 60.0);
                theta_vs_mom_prot[i] = std::make_shared<TH2D>(
                        Form("theta_vs_mom_prot_%d", i + 1), Form("theta_vs_mom_prot_%d", i + 1), bins,
                        zero, p_max, bins, zero, 140.0);
                theta_vs_mom_pip[i] = std::make_shared<TH2D>(
                        Form("theta_vs_mom_pip_%d", i + 1), Form("theta_vs_mom_pip_%d", i + 1), bins,
                        zero, p_max, bins, zero, 140.0);
                theta_vs_mom_pim[i] = std::make_shared<TH2D>(
                        Form("theta_vs_mom_pim_%d", i + 1), Form("theta_vs_mom_pim_%d", i + 1), bins,
                        zero, p_max, bins, zero, 180.0);

                //int Phi_min_max_val = 360.0;
                //Phi_min_max_val = 180 + i* 60;
                // if (i>=3)
                //         Phi_min_max_val = (i-3)* 60;
                pim_phi_vs_theta_rec_FD_sec[i] = std::make_shared<TH2D>(
                        Form("pim_phi_vs_theta_rec_FD_%d", i + 1), Form("pim_phi_vs_theta_rec_FD_%d", i + 1), bins,
                        zero, 60, bins, 0, 360);

                pim_phi_vs_theta_measured_FD_sec[i] = std::make_shared<TH2D>(
                        Form("pim_phi_vs_theta_measured_FD_%d", i + 1), Form("pim_phi_vs_theta_measured_FD_%d", i + 1), bins,
                        zero, 60, bins, 0, 360);

                pim_phi_vs_theta_rec_FD_after_exclusive_sec[i] = std::make_shared<TH2D>(
                        Form("pim_phi_vs_theta_rec_FD_after_MMSQ_%d", i + 1), Form("pim_phi_vs_theta_rec_FD_after_MMSQ_%d", i + 1), bins,
                        zero, 60, bins, 0, 360);

                pim_phi_vs_theta_measured_FD_after_exclusive_sec[i] = std::make_shared<TH2D>(
                        Form("pim_phi_vs_theta_measured_FD_after_MMSQ_%d", i + 1), Form("pim_phi_vs_theta_measured_FD_after_MMSQ_%d", i + 1), bins,
                        zero, 60, bins, 0,  360);


        }
}

void Histogram::makeHists_deltat() {
        std::string tof = "";
        for (short p = 0; p < particle_num; p ++) {
                for (short c = 0; c < charge_num; c ++) {
                        for (short i = 0; i < with_id_num; i ++) {
                                tof = "ftof";
                                delta_t_hist[p][c][i][0] = std::make_shared<TH2D>(
                                        Form("delta_t_%s_%s_%s_%s", tof.c_str(), particle_name[p].c_str(),
                                             charge_name[c].c_str(), id_name[i].c_str()),
                                        Form("#Deltat %s %s %s %s", tof.c_str(), particle_name[p].c_str(),
                                             charge_name[c].c_str(), id_name[i].c_str()),
                                        bins, p_min, p_max, bins, Dt_min, Dt_max);

                                tof = "ctof";
                                delta_t_hist[p][c][i][1] = std::make_shared<TH2D>(
                                        Form("delta_t_%s_%s_%s_%s", tof.c_str(), particle_name[p].c_str(),
                                             charge_name[c].c_str(), id_name[i].c_str()),
                                        Form("#Deltat %s %s %s %s", tof.c_str(), particle_name[p].c_str(),
                                             charge_name[c].c_str(), id_name[i].c_str()),
                                        bins, 0, 3.0, bins, -6.0, 6.0);
                        }
                }
        }
}

void Histogram::Fill_deltat_pi(const std::shared_ptr<Branches12> &data,
                               const std::shared_ptr<Delta_T> &dt, int part) {
        auto _cuts = std::make_unique<Cuts>(data, dt);
        int charge = data->charge(part);
        bool fc = dt->ctof();
        int pid = data->pid(part);
        float mom = data->p(part);
        float time = NAN;
        if (fc)
                time = dt->dt_ctof_Pi();
        else
                time = dt->dt_Pi();

        if (charge == 1) {
                delta_t_hist[1][0][0][fc] -> Fill(mom, time);
                if (_cuts->IsPip(part))
                        delta_t_hist[1][0][1][fc] -> Fill(mom, time);
                else
                        delta_t_hist[1][0][2][fc] -> Fill(mom, time);
        } else if (charge == -1) {
                delta_t_hist[1][1][0][fc] -> Fill(mom, time);
                //if (_cuts->IsPim(part))
                if (pid !=ELECTRON)
                        delta_t_hist[1][1][1][fc] -> Fill(mom, time);
                else
                        delta_t_hist[1][1][2][fc] -> Fill(mom, time);
        }
}

void Histogram::Fill_deltat_prot(const std::shared_ptr<Branches12> &data,
                                 const std::shared_ptr<Delta_T> &dt, int part) {
        auto _cuts = std::make_unique<Cuts>(data, dt);
        int status = abs(data->status(part));
        int charge = data->charge(part);
        bool fc = dt->ctof();
        int pid = data->pid(part);
        float mom = data->p(part);
        float time = NAN;

        if (fc)
                //if(status >=4000 && status < 6000){
                time = dt->dt_ctof_P();
        //}
        else
                time = dt->dt_P();

        if (charge == 1) {
                delta_t_hist[2][0][0][fc] -> Fill(mom, time);
                if (_cuts->IsProton(part))
                        delta_t_hist[2][0][1][fc] -> Fill(mom, time);
                else
                        delta_t_hist[2][0][2][fc] -> Fill(mom, time);

                delta_t_hist[2][1][0][fc] -> Fill(mom, time);
                if (pid == PROTON)
                        delta_t_hist[2][1][1][fc] -> Fill(mom, time);
                else
                        delta_t_hist[2][1][2][fc] -> Fill(mom, time);
        }
}

// do
// pid == 11 at first;
// skim vs pid at 0 is 11. compare.
void Histogram::Write_deltat() {
        TDirectory *ftof_folder = RootOutputFile->mkdir("ftof");
        ftof_folder->cd();
        for (short p = 0; p < particle_num; p ++) {
                for (short c = 0; c < charge_num; c ++) {
                        for (short i = 0; i < with_id_num; i ++) {
                                delta_t_hist[p][c][i][0] -> SetXTitle("Momentum (GeV)");
                                delta_t_hist[p][c][i][0] -> SetYTitle("#Deltat");
                                delta_t_hist[p][c][i][0] -> SetOption("COLZ1");
                                if (delta_t_hist[p][c][i][0]->GetEntries() > 1)
                                        delta_t_hist[p][c][i][0]->Write();
                        }
                }
        }
        TDirectory *ctof_folder = RootOutputFile->mkdir("ctof");
        ctof_folder->cd();
        for (short p = 0; p < particle_num; p ++) {
                for (short c = 0; c < charge_num; c ++) {
                        for (short i = 0; i < with_id_num; i ++) {
                                delta_t_hist[p][c][i][1] -> SetXTitle("Momentum (GeV)");
                                delta_t_hist[p][c][i][1] -> SetYTitle("#Deltat");
                                delta_t_hist[p][c][i][1] -> SetOption("COLZ1");
                                if (delta_t_hist[p][c][i][1]->GetEntries() > 1)
                                        delta_t_hist[p][c][i][1]->Write();
                        }
                }
        }
}

void Histogram::makeHists_MomVsBeta() {
        for (short p = 0; p < particle_num; p ++) {
                for (short c = 0; c < charge_num; c ++) {
                        for (short i = 0; i < with_id_num; i ++) {
                                momvsbeta_hist[p][c][i] = std::make_shared<TH2D>(
                                        Form("mom_vs_beta_%s_%s_%s", particle_name[p].c_str(),
                                             charge_name[c].c_str(), id_name[i].c_str()),
                                        Form("Momentum vs #beta %s %s %s", particle_name[p].c_str(),
                                             charge_name[c].c_str(), id_name[i].c_str()),
                                        bins, p_min, p_max, bins, zero, 1.2);
                        }
                }
        }
}

void Histogram::Fill_MomVsBeta(const std::shared_ptr<Branches12> &data,
                               int part) {
        int good_ID = 0;
        float beta = data->beta(part);
        float mom = data->p(part);
        int charge = data->charge(part);
        int pid = data->pid(part);
        if (beta != 0) {
                momentum->Fill(mom);
                for (short p = 0; p < particle_num; p ++) {
                        switch (p) {
                        case 0:
                                good_ID = ELECTRON;
                                break;
                        case 1:
                                good_ID = PIP;
                                break;
                        case 2:
                                good_ID = PROTON;
                                break;
                        case 3:
                                good_ID = KP;
                                break;
                        }

                        momvsbeta_hist[p][0][0] -> Fill(mom, beta);
                        if (good_ID == abs(pid)) {
                                momvsbeta_hist[p][0][1] -> Fill(mom, beta);
                        } else {
                                momvsbeta_hist[p][0][2] -> Fill(mom, beta);
                        }

                        if (charge == -1) {
                                momvsbeta_hist[p][2][0] -> Fill(mom, beta);
                                if (-good_ID == pid) {
                                        momvsbeta_hist[p][2][1]->Fill(mom, beta);
                                } else {
                                        momvsbeta_hist[p][2][2]->Fill(mom, beta);
                                }
                        } else if (charge == 1) {
                                momvsbeta_hist[p][1][0] -> Fill(mom, beta);
                                if (good_ID == pid) {
                                        momvsbeta_hist[p][1][1]->Fill(mom, beta);
                                } else {
                                        momvsbeta_hist[p][1][2]->Fill(mom, beta);
                                }
                        }
                }
        }
}

void Histogram::Write_MomVsBeta() {
        momentum->SetXTitle("Momentum (GeV)");
        momentum->Write();
        for (short p = 0; p < particle_num; p++) {
                for (short c = 0; c < charge_num; c++) {
                        for (short i = 0; i < with_id_num; i++) {
                                momvsbeta_hist[p][c][i]->SetXTitle("Momentum (GeV)");
                                momvsbeta_hist[p][c][i]->SetYTitle("#beta");
                                momvsbeta_hist[p][c][i]->SetOption("COLZ1");
                                momvsbeta_hist[p][c][i]->Write();
                        }
                }
        }
}
