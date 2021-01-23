/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef HIST_H_GUARD
#define HIST_H_GUARD
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TPaveStats.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TThread.h"
#include "colors.hpp"
#include "constants.hpp"
#include "cuts.hpp"
#include "deltat.hpp"
#include "reaction.hpp"
#include <mutex>

using namespace std;

using TH2D_ptr = std::shared_ptr<TH2D>;
using TH1D_ptr = std::shared_ptr<TH1D>;
using THnSparse_ptr = std::shared_ptr<THnSparse>;
using TGraph_ptr = std::shared_ptr<TGraph>;

class Histogram {
protected:
std::shared_ptr<TFile> RootOutputFile;
std::shared_ptr<TCanvas> def;

int bins = 500;
double p_min = 0.0;
double p_max = 10.0;
double Dt_max = 10.0;
double Dt_min = -Dt_max;
double q2_max = 12.0;
double w_max = 5.50;

double zero = 0.0;

static const short particle_num = 4;   // 0-e 1-Pi 2-P 3-K
std::string particle_name[particle_num] = {"e", "pi", "P", "K"};
static const short charge_num = 2;   // 0-pos 1-neg
std::string charge_name[charge_num] = {"positive", "negative"};
static const short with_id_num = 3;   // 0-without 1-with 2-anti
std::string id_name[with_id_num] = {"withoutID", "withID", "antiID"};
static const short num_sectors = 6;
std::string sec_name[num_sectors] = {"1", "2", "3", "4", "5", "6"};

static const short CUTS = 2;
enum cuts { before_cut, after_cut };
std::mutex mutex;

static const short w_range_num = 3;
std::string w_range_name[w_range_num] = {
        "all_W_range", " W < 2.5 ",
        " W > 2.5 "
};                  //{" W<2.0 ", " 2.0<W<2.5 ", " 2.5<W<3.0 ", " 3.0<W<3.5 "};
static const short q2_range_num = 3;
std::string q2_range_name[q2_range_num] =   /*{" Q2<1.0 ",     " 1.0<Q2<2.0 ",
                                               " 2.0<Q2<3.0 ", " 3.0<Q2<4.0 ",
                                               " 4.0<Q2<5.0 ", " 5.0<Q2<6.0 ",
                                               " 6.0<Q2<7.0 ", " 7.0<Q2<8.0 ",
                                               " 8.0<Q2<9.0 ", " Q2>9.0 "};*/
{"all_Q2_range ", " Q2 < 4.5 ", " Q2 > 4.5 "};

static const short inv_Ppip_range_num = 2;
std::string inv_Ppip_range_name[inv_Ppip_range_num] = {" M[Ppip]<1.5 GeV ",
                                                       " M[Ppip]>1.5 GeV "};

static const short inv_pip_pim_range_num = 2;
std::string inv_pip_pim_range_name[inv_pip_pim_range_num] = {
        " M[pi+pi-]<0.9 GeV ", " M[pi+pi-]>0.9 GeV "
};

static const short theta_pim_range_num = 2;
std::string theta_pim_range_name[theta_pim_range_num] = {
        " theta[pi-]<90 deg ", " theta[pi-]>90 deg "
};

static const short phi_pim_range_num = 2;
std::string phi_pim_range_name[phi_pim_range_num] = {" phi[pi-]<90 deg ",
                                                     " phi[pi-]>90 deg "};
static const short alpha_pim_range_num = 2;
std::string alpha_pim_range_name[alpha_pim_range_num] = {
        " alpha[pi-]<180 deg ", " alpha[pi-]>180 deg "
};

static const short NUM_CONDITIONS = 1;
std::string NUM_CONDITIONS_NAME[NUM_CONDITIONS] = {
        "twoPi_event" /*, "missingPim"*/
};
// Kinematics

static const short q2_bin = 11;
static const short w_bin = 64;
THnSparse *threeDHist[q2_bin][w_bin];;
THnSparse *sevenDHist_pim[q2_bin][w_bin];;
THnSparse *sevenD_Hist_thrown_pim[q2_bin][w_bin];;
THnSparse *sevenDHist_pip[q2_bin][w_bin];;
THnSparse *sevenD_Hist_thrown_pip[q2_bin][w_bin];;
THnSparse *sevenDHist_prot[q2_bin][w_bin];
THnSparse *sevenD_Hist_thrown_prot[q2_bin][w_bin];
static const short NUM_CUT = 2;

static const short theta_bin_NUM = 18;
TH1D_ptr theta_pim_measured_3_sigma[theta_bin_NUM];
std::string theta_bin_NAME[theta_bin_NUM] = {
        "0-10","10-20","20-30","30-40","40-50","50-60","60-70","70-80","80-90","90-100","100-110","110-120","120-130","130-140","140-150","150-160","160-170","170-180"
};

//  TH2D_ptr sf_hist = std::make_shared<TH2D>("SF", "SF", 500, 0, 10.5, 500,
//  0, 0.5);
//
TH1D_ptr diff_E_P_x_mu_hist_ = std::make_shared<TH1D>(
        "diff_E_P_x_mu ", "diff_E_P_x_mu ", 500, -1.0, 1.0);
TH1D_ptr P_x_mu = std::make_shared<TH1D>("Mom P ", "Mom P ", 500, -2.0, 10.0);
// // TH1D_ptr diff_E_P_x_mu =
// //     std::make_shared<TH1D>("Energy - mom (x_mu = e(p,p')e')", "Energy - mom
// //     (x_mu = e(p,p')e')", 500, -0.50, 1.0);
// // TH2D_ptr mom_vs_E_x_mu =
// //     std::make_shared<TH2D>("Mom_vs_Energy", " Mom_vs_Energy_component",
// //     500, -1.0, 2.0, 500, 0.0, 2.0);
// // TH1D_ptr diff_theta_in_x_mu = std::make_shared<TH1D>("#Delata#Theta x_mu
// // and initial electron",
// //                                                      "#Delata#Theta x_mu
// //                                                      and initial electron",
// //                                                      500, 0, 180);
// // TH2D_ptr Dthtea_vs_Dphi = std::make_shared<TH2D>("#Delata#Theta vs
// // #Delata#Phi", "#Delata#Theta vs #Delata#Phi", 500,
// //                                                  -180, 180, 500, -180,
// //
// TH1D_ptr theta_pim_rec = std::make_shared<TH1D>("theta pi- rec ", "theta pi- rec ", 500, 0.0, 180.0);
// TH1D_ptr theta_pim_measured = std::make_shared<TH1D>("theta pi- measured ", "theta pi- measured ", 500, 0.0, 180.0);
// TH2D_ptr pim_theta_measured_vs_rec =
//         std::make_shared<TH2D>("theta pi- rec vs measured", " theta pi- rec vs measured",
//                                500, 0.0, 180.0, 500, 0.0, 180.0);
// TH2D_ptr pim_theta_rec_vs_mom =
//         std::make_shared<TH2D>("theta_pi- rec vs momentum", " theta_pi- rec vs momentum",
//                                500, 0.0, 9.0, 500, 0.0, 180.0);
// TH2D_ptr pim_theta_measured_vs_mom =
//         std::make_shared<TH2D>("theta_pi- measured vs momentum", " theta_pi- measured vs momentum",
//                                500, 0.0, 9.0, 500, 0.0, 180.0);
//
// TH1D_ptr phi_pim_rec = std::make_shared<TH1D>("phi pi- rec ", "phi pi- rec ", 500, 0.0, 360.0);
// TH1D_ptr phi_pim_measured = std::make_shared<TH1D>("phi pi- measured ", "phi pi- measured ", 500, 0.0, 360.0);
// TH2D_ptr pim_phi_measured_vs_rec =
//         std::make_shared<TH2D>("phi pi- rec vs measured", " phi pi- rec vs measured",
//                                500, 0.0, 360.0, 500, 0.0, 360.0);
// TH2D_ptr pim_phi_rec_vs_mom =
//         std::make_shared<TH2D>("phi_pi- rec vs momentum", " phi_pi- rec vs momentum",
//                                500, 0.0, 9.0, 500, 0.0, 360.0);
// TH2D_ptr pim_phi_measured_vs_mom =
//         std::make_shared<TH2D>("phi_pi- measured vs momentum", " phi_pi- measured vs momentum",
//                                500, 0.0, 9.0, 500, 0.0, 360.0);
//
// TH2D_ptr pim_theta_vs_phi_rec =
//         std::make_shared<TH2D>("theta_pi- vs phi pi- rec", " theta_pi- ve phi pi- rec",
//                                500, 0.0, 360.0, 500, 0.0, 180.0);
//
// TH2D_ptr pim_theta_vs_phi_measured =
//         std::make_shared<TH2D>("theta_pi- vs phi pi- measured", " theta_pi- ve phi pi- measured",
//                                500, 0.0, 360.0, 500, 0.0, 180.0);
//
//
// TH1D_ptr theta_pim_rec_after_mmsq_applied = std::make_shared<TH1D>("theta pi- rec with MMSQ cut ", "theta pi- rec with MMSQ cut ", 500, 0.0, 180.0);
// TH1D_ptr theta_pim_measured_after_mmsq_applied = std::make_shared<TH1D>("theta pi- measured with MMSQ cut ", "theta pi- measured with MMSQ cut ", 500, 0.0, 180.0);
// TH2D_ptr pim_theta_measured_vs_rec_after_mmsq_applied =
//         std::make_shared<TH2D>("theta pi- rec vs measured with MMSQ cut", " theta pi- rec vs measured with MMSQ cut",
//                                500, 0.0, 180.0, 500, 0.0, 180.0);
//
// TH2D_ptr pim_theta_measured_vs_rec_after_y_equal_mc_cuts =
//         std::make_shared<TH2D>("theta pi- rec vs measured with cuts", " theta pi- rec vs measured with cuts",
//                                500, 0.0, 180.0, 500, 0.0, 180.0);
//
// TH2D_ptr pim_theta_rec_vs_mom_after_mmsq_applied =
//         std::make_shared<TH2D>("theta_pi- rec vs momentum with MMSQ cut", " theta_pi- rec vs momentum with MMSQ cut",
//                                500, 0.0, 9.0, 500, 0.0, 180.0);
// TH2D_ptr pim_theta_measured_vs_mom_after_mmsq_applied =
//         std::make_shared<TH2D>("theta_pi- measured vs momentum with MMSQ cut", " theta_pi- measured vs momentum with MMSQ cut",
//                                500, 0.0, 9.0, 500, 0.0, 180.0);
//
// TH1D_ptr phi_pim_rec_after_mmsq_applied = std::make_shared<TH1D>("phi pi- rec with MMSQ cut ", "phi pi- rec with MMSQ cut ", 500, 0.0, 360.0);
// TH1D_ptr phi_pim_measured_after_mmsq_applied = std::make_shared<TH1D>("phi pi- measured with MMSQ cut ", "phi pi- measured with MMSQ cut ", 500, 0.0, 360.0);
// TH2D_ptr pim_phi_measured_vs_rec_after_mmsq_applied =
//         std::make_shared<TH2D>("phi pi- rec vs measured with MMSQ cut", " phi pi- rec vs measured with MMSQ cut",
//                                500, 0.0, 360.0, 500, 0.0, 360.0);
//
// TH2D_ptr pim_phi_measured_vs_rec_after_y_equal_mc_cuts =
//         std::make_shared<TH2D>("phi pi- rec vs measured with cuts", " phi pi- rec vs measured with cuts",
//                                500, 0.0, 360.0, 500, 0.0, 360.0);
//
// TH2D_ptr pim_phi_rec_vs_mom_after_mmsq_applied =
//         std::make_shared<TH2D>("phi_pi- rec vs momentum with MMSQ cut", " phi_pi- rec vs momentum with MMSQ cut",
//                                500, 0.0, 9.0, 500, 0.0, 360.0);
// TH2D_ptr pim_phi_measured_vs_mom_after_mmsq_applied =
//         std::make_shared<TH2D>("phi_pi- measured vs momentum with MMSQ cut", " phi_pi- measured vs momentum with MMSQ cut",
//                                500, 0.0, 9.0, 500, 0.0, 360.0);
//
// TH2D_ptr pim_theta_vs_phi_rec_after_mmsq_applied =
//         std::make_shared<TH2D>("theta_pi- vs phi pi- rec with MMSQ cut", " theta_pi- ve phi pi- rec with MMSQ cut",
//                                500, 0.0, 360.0, 500, 0.0, 180.0);
//
// TH2D_ptr pim_theta_vs_phi_measured_after_mmsq_applied =
//         std::make_shared<TH2D>("theta_pi- vs phi pi- measured with MMSQ cut", " theta_pi- ve phi pi- measured with MMSQ cut",
//                                500, 0.0, 360.0, 500, 0.0, 180.0);

TH2D_ptr theta_vs_mom_elec[num_sectors];
TH2D_ptr theta_vs_mom_prot[num_sectors];
TH2D_ptr theta_vs_mom_pip[num_sectors];
TH2D_ptr theta_vs_mom_pim[num_sectors];
TH2D_ptr pim_phi_vs_theta_rec_FD_sec[num_sectors];
TH2D_ptr pim_phi_vs_theta_rec_FD_after_exclusive_sec[num_sectors];
TH2D_ptr pim_phi_vs_theta_measured_FD_sec[num_sectors];
TH2D_ptr pim_phi_vs_theta_measured_FD_after_exclusive_sec[num_sectors];

//missingPiP
static const short HADRON_NUM = 3;
std::string HADRON_NAME[HADRON_NUM] = {
        "mProt","mPip","mPim"
};
static const short EFF_CONDITIONS_NUM_ALL = 3;
std::string EFF_CONDITIONS_NAME_ALL[EFF_CONDITIONS_NUM_ALL] = {
        "missing","exclusive","after_MMSQ_exclusive_cuts"
};
TH2D_ptr pip_theta_rec_vs_mom[EFF_CONDITIONS_NUM_ALL];
TH2D_ptr prot_theta_rec_vs_mom[EFF_CONDITIONS_NUM_ALL];
TH2D_ptr pip_theta_measured_vs_mom[EFF_CONDITIONS_NUM_ALL];
TH2D_ptr prot_theta_measured_vs_mom[EFF_CONDITIONS_NUM_ALL];

static const short DETECTOR_NUM_PROT = 2;
std::string DETECTOR_NAME_PROT[DETECTOR_NUM_PROT] = {"prot_CD","prot_FD"};
static const short DETECTOR_NUM_PIP = 2;
std::string DETECTOR_NAME_PIP[DETECTOR_NUM_PIP] = {"pip_CD","pip_FD"};
static const short DETECTOR_NUM_PIM = 2;
std::string DETECTOR_NAME_PIM[DETECTOR_NUM_PIM] = {"pim_CD","pim_FD"};
TH1D_ptr theta_measured_minus_rec[HADRON_NUM][EFF_CONDITIONS_NUM_ALL][DETECTOR_NUM_PROT][DETECTOR_NUM_PIP][DETECTOR_NUM_PIM];
TH1D_ptr phi_measured_minus_rec[HADRON_NUM][EFF_CONDITIONS_NUM_ALL][DETECTOR_NUM_PROT][DETECTOR_NUM_PIP][DETECTOR_NUM_PIM];
TH1D_ptr mom_measured_minus_rec[HADRON_NUM][EFF_CONDITIONS_NUM_ALL] [DETECTOR_NUM_PROT][DETECTOR_NUM_PIP][DETECTOR_NUM_PIM];
TH2D_ptr theta_rec_vs_mom[HADRON_NUM][EFF_CONDITIONS_NUM_ALL][DETECTOR_NUM_PIM][DETECTOR_NUM_PIP][DETECTOR_NUM_PIM];

// TH1D_ptr theta_measured_minus_rec_FD[HADRON_NUM][EFF_CONDITIONS_NUM][DETECTOR_NUM_PROT][DETECTOR_NUM_PIM];
// TH1D_ptr phi_measured_minus_rec_FD[HADRON_NUM][EFF_CONDITIONS_NUM][DETECTOR_NUM_PROT][DETECTOR_NUM_PIM];
// TH1D_ptr mom_measured_minus_rec_FD[HADRON_NUM][EFF_CONDITIONS_NUM][DETECTOR_NUM_PROT][DETECTOR_NUM_PIM];
// TH2D_ptr theta_rec_vs_mom_FD[HADRON_NUM][EFF_CONDITIONS_NUM][DETECTOR_NUM_PIM][DETECTOR_NUM_PIM];
//
// //missingProt
// std::string EFF_CONDITIONS_NAME_PROT[EFF_CONDITIONS_NUM] = {
//         "missingProt","exclusive","after_MMSQ_exclusive_cuts"
// };
// TH1D_ptr prot_theta_measured_minus_rec_CD[EFF_CONDITIONS_NUM][DETECTOR_NUM_PIM][DETECTOR_NUM_PIP];
// TH1D_ptr prot_phi_measured_minus_rec_CD[EFF_CONDITIONS_NUM][DETECTOR_NUM_PIM][DETECTOR_NUM_PIP];
// TH1D_ptr prot_mom_measured_minus_rec_CD[EFF_CONDITIONS_NUM][DETECTOR_NUM_PIM][DETECTOR_NUM_PIP];
// TH2D_ptr prot_theta_rec_vs_mom_CD[EFF_CONDITIONS_NUM][DETECTOR_NUM_PIM][DETECTOR_NUM_PIP];
//
// TH1D_ptr prot_theta_measured_minus_rec_FD[EFF_CONDITIONS_NUM][DETECTOR_NUM_PIM][DETECTOR_NUM_PIP];
// TH1D_ptr prot_phi_measured_minus_rec_FD[EFF_CONDITIONS_NUM][DETECTOR_NUM_PIM][DETECTOR_NUM_PIP];
// TH1D_ptr prot_mom_measured_minus_rec_FD[EFF_CONDITIONS_NUM][DETECTOR_NUM_PIM][DETECTOR_NUM_PIP];
// TH2D_ptr prot_theta_rec_vs_mom_FD[EFF_CONDITIONS_NUM][DETECTOR_NUM_PIM][DETECTOR_NUM_PIP];
static const short EFF_CONDITIONS_NUM = 3;
std::string EFF_CONDITIONS_NAME[EFF_CONDITIONS_NUM] = {
        //"missingPim","exclusive","after_MMSQ_exclusive_cuts"
        "missingPim","exclusive","after_MMSQ_exclusive_cuts"
};
THnSparse_ptr background_4DHist_mpim[EFF_CONDITIONS_NUM];

static const short EFF_CONDITIONS_NUM_MMSQ = 4;
std::string EFF_CONDITIONS_NAME_MMSQ[EFF_CONDITIONS_NUM_MMSQ] = {
        //"missingPim","exclusive","after_MMSQ_exclusive_cuts"
        "missingPim"," missingPim_after_MMSQ_cuts","exclusive","after_MMSQ_exclusive_cuts"
};
static const short THETA_BINS_NUM = 10;
std::string THETA_BINS_NAME[THETA_BINS_NUM] = {
        "0-15","15-30","30-45","45-60","60-75","75-90","90-105","105-120","120-135","135-150"
};
TH1D_ptr eff_check_mPim_all_W;
TH1D_ptr eff_check_mPim_all_W_with_mmsq_cuts;
TH1D_ptr eff_check_exclusive_all_W;
TH1D_ptr eff_check_exclusive_MMSQ_cuts_all_W;
TH1D_ptr eff_check_mPim[EFF_CONDITIONS_NUM_MMSQ][w_bin][THETA_BINS_NUM];
//TH1D_ptr eff_check_mPim[EFF_CONDITIONS_NUM][w_bin];

//
// TH1D_ptr pim_theta_measured_minus_rec_CD[EFF_CONDITIONS_NUM][DETECTOR_NUM_PROT][DETECTOR_NUM_PIP];
// TH1D_ptr pim_phi_measured_minus_rec_CD[EFF_CONDITIONS_NUM][DETECTOR_NUM_PROT][DETECTOR_NUM_PIP];
// TH1D_ptr pim_mom_measured_minus_rec_CD[EFF_CONDITIONS_NUM][DETECTOR_NUM_PROT][DETECTOR_NUM_PIP];
// TH1D_ptr pim_theta_measured_minus_rec_FD[EFF_CONDITIONS_NUM][DETECTOR_NUM_PROT][DETECTOR_NUM_PIP];
// TH1D_ptr pim_phi_measured_minus_rec_FD[EFF_CONDITIONS_NUM][DETECTOR_NUM_PROT][DETECTOR_NUM_PIP];
// TH1D_ptr pim_mom_measured_minus_rec_FD[EFF_CONDITIONS_NUM][DETECTOR_NUM_PROT][DETECTOR_NUM_PIP];
//missingPim
//CD

TH1D_ptr MM2_exclusive_hist_CD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_theta_rec_vs_mom_CD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_phi_rec_vs_mom_CD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_theta_measured_vs_mom_CD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_phi_measured_vs_mom_CD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_phi_vs_theta_rec_CD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_E_vs_mom_rec_CD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_E_vs_mom_measured_CD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_E_measured_vs_rec_CD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_mom_measured_vs_rec_CD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_theta_measured_vs_rec_CD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_phi_measured_vs_rec_CD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_phi_vs_theta_measured_CD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_theta_measured_minus_rec_vs_measured_CD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_theta_measured_minus_rec_vs_rec_CD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_theta_measured_vs_thrown_CD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_theta_rec_vs_thrown_CD[EFF_CONDITIONS_NUM];
//FD

TH1D_ptr MM2_exclusive_hist_FD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_theta_rec_vs_mom_FD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_phi_rec_vs_mom_FD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_theta_measured_vs_mom_FD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_phi_measured_vs_mom_FD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_phi_vs_theta_rec_FD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_E_vs_mom_rec_FD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_E_vs_mom_measured_FD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_E_measured_vs_rec_FD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_mom_measured_vs_rec_FD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_theta_measured_vs_rec_FD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_phi_measured_vs_rec_FD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_phi_vs_theta_measured_FD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_theta_measured_minus_rec_vs_measured_FD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_theta_measured_minus_rec_vs_rec_FD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_theta_measured_vs_thrown_FD[EFF_CONDITIONS_NUM];
TH2D_ptr pim_theta_rec_vs_thrown_FD[EFF_CONDITIONS_NUM];


TH1D_ptr E_x_mu_hist[NUM_CONDITIONS];
TH1D_ptr diff_E2_P2_x_mu_hist[NUM_CONDITIONS];
TH1D_ptr diff_E_P_x_mu_hist[NUM_CONDITIONS];
TH2D_ptr mom_vs_E_x_mu_hist[NUM_CONDITIONS];
TH1D_ptr theta_elec_hist[NUM_CONDITIONS];
TH1D_ptr theta_x_mu_hist[NUM_CONDITIONS];
TH1D_ptr diff_theta_elec_x_mu_hist[NUM_CONDITIONS];
TH2D_ptr MM2_VS_W_x_mu_hist[NUM_CONDITIONS];

TH1D_ptr mass_pi0_any_event[num_sectors];
TH1D_ptr mass_pi0_hist_before_mmsq_cut[num_sectors];
TH1D_ptr mass_pi0_hist_after_mmsq_cut[num_sectors];

TH1D_ptr weight_hist;
TH1D_ptr momentum;
TH1D_ptr W_hist;
TH1D_ptr Q2_hist;
TH2D_ptr W_vs_q2;
TH1D_ptr W_P2pi_hist;

TH1D_ptr W_hist_ftof;
TH2D_ptr W_vs_q2_ftof;
TH1D_ptr W_hist_ctof;
TH2D_ptr W_vs_q2_ctof;

TH1D_ptr W_thrown;
TH2D_ptr W_vs_Q2_thrown;

TH1D_ptr vz_position[CUTS];
TH2D_ptr pcal_sec[CUTS];
TH2D_ptr dcr1_sec[CUTS];
TH2D_ptr dcr2_sec[CUTS];
TH2D_ptr dcr3_sec[CUTS];

TH1D_ptr theta_prot;
TH1D_ptr theta_pip;
TH1D_ptr theta_pim;

TH2D_ptr Theta_prot_cm_vs_mom_prot;
TH2D_ptr Theta_pip_cm_vs_mom_pip;
TH2D_ptr Theta_pim_cm_vs_mom_pim;

TH2D_ptr Theta_prot_lab_vs_mom_prot;
TH2D_ptr Theta_pip_lab_vs_mom_pip;
TH2D_ptr Theta_pim_lab_vs_mom_pim;

TH2D_ptr Theta_prot_thrown_cm_vs_mom_prot;
TH2D_ptr Theta_pip_thrown_cm_vs_mom_pip;
TH2D_ptr Theta_pim_thrown_cm_vs_mom_pim;

TH2D_ptr Theta_prot_thrown_lab_vs_mom_prot;
TH2D_ptr Theta_pip_thrown_lab_vs_mom_pip;
TH2D_ptr Theta_pim_thrown_lab_vs_mom_pim;

TH1D_ptr Phi_gamma;
TH1D_ptr Phi_prot;
TH1D_ptr Phi_pip;
TH1D_ptr Phi_pim;

TH1D_ptr alpha_pim;
TH1D_ptr alpha_pip;
TH1D_ptr alpha_prot;

TH1D_ptr theta_prot_mc;
TH1D_ptr theta_pip_mc;
TH1D_ptr theta_pim_mc;

TH1D_ptr Phi_gamma_mc;
TH1D_ptr Phi_prot_mc;
TH1D_ptr Phi_pip_mc;
TH1D_ptr Phi_pim_mc;

TH1D_ptr alpha_pim_mc;
TH1D_ptr alpha_pip_mc;
TH1D_ptr alpha_prot_mc;

TH1D_ptr theta_prot_thrown;
TH1D_ptr theta_pip_thrown;
TH1D_ptr theta_pim_thrown;

TH1D_ptr Phi_gamma_thrown;
TH1D_ptr Phi_prot_thrown;
TH1D_ptr Phi_pip_thrown;
TH1D_ptr Phi_pim_thrown;

TH1D_ptr alpha_pim_thrown;
TH1D_ptr alpha_pip_thrown;
TH1D_ptr alpha_prot_thrown;

TH2D_ptr W_vs_q2_sec[num_sectors];
TH1D_ptr W_sec[num_sectors];

TH1D_ptr W_det[3];
TH2D_ptr WQ2_det[3];

TH1D_ptr W_hist_singlePi;
TH1D_ptr Q2_hist_singlePi;
TH2D_ptr W_vs_q2_singlePi;
TH2D_ptr W_vs_q2_singlePi_sec[num_sectors];
TH1D_ptr W_singlePi_sec[num_sectors];
TH2D_ptr W_vs_MM_singlePi[num_sectors];

TH2D_ptr W_vs_q2_Npip_sec[num_sectors];
TH1D_ptr W_Npip_sec[num_sectors];
TH1D_ptr MM_Npip_sec[num_sectors];

TH1D_ptr MM_neutron;
TH1D_ptr MM_neutron_sec[num_sectors];

TH2D_ptr W_vs_q2_twoPi_sec[num_sectors];
TH1D_ptr W_twoPi_sec[num_sectors];
TH2D_ptr W_vs_MM_twoPi[num_sectors];
TH1D_ptr MM_twoPi;
TH1D_ptr MM_twoPi_sec[num_sectors];
TH1D_ptr MM2_twoPi;
TH1D_ptr MM2_twoPi_sec[num_sectors];
TH1D_ptr MM2_twoPi_missingPip;
TH1D_ptr MM2_twoPi_missingPip_sec[num_sectors];
TH1D_ptr MM2_twoPi_missingProt;
TH1D_ptr MM2_twoPi_missingProt_sec[num_sectors];
TH1D_ptr W_hist_twoPi;
TH1D_ptr Q2_hist_twoPi;
TH2D_ptr W_vs_q2_twoPi;
TH2D_ptr W_vs_q2_twoPi_thrown;

TH1D_ptr inv_mass_P_pip[w_range_num][q2_range_num];
TH1D_ptr inv_mass_P_pim[w_range_num][q2_range_num];
TH1D_ptr inv_mass_pip_pim[w_range_num][q2_range_num];

TH2D_ptr theta_P_vs_mass_pip_pim[w_range_num][q2_range_num];
TH2D_ptr theta_pim_vs_mass_Ppip[w_range_num][q2_range_num];
TH2D_ptr theta_pip_vs_mass_Ppim[w_range_num][q2_range_num];
TH2D_ptr theta_P_lab_vs_mass_pip_pim[w_range_num][q2_range_num];
TH2D_ptr theta_pim_lab_vs_mass_Ppip[w_range_num][q2_range_num];
TH2D_ptr theta_pip_lab_vs_mass_Ppim[w_range_num][q2_range_num];

TH2D_ptr W_vs_q2_twoPi_sec_thrown[num_sectors];
TH1D_ptr W_twoPi_sec_thrown[num_sectors];
TH2D_ptr W_vs_MM_twoPi_thrown[num_sectors];
TH1D_ptr MM_twoPi_sec_thrown[num_sectors];
TH1D_ptr MM2_twoPi_sec_thrown[num_sectors];

TH1D_ptr inv_mass_P_pip_thrown[w_range_num][q2_range_num];
TH1D_ptr inv_mass_P_pim_thrown[w_range_num][q2_range_num];
TH1D_ptr inv_mass_pip_pim_thrown[w_range_num][q2_range_num];

TH2D_ptr theta_P_vs_mass_pip_pim_thrown[w_range_num][q2_range_num];
TH2D_ptr theta_pim_vs_mass_Ppip_thrown[w_range_num][q2_range_num];
TH2D_ptr theta_pip_vs_mass_Ppim_thrown[w_range_num][q2_range_num];
TH2D_ptr theta_P_lab_vs_mass_pip_pim_thrown[w_range_num][q2_range_num];
TH2D_ptr theta_pim_lab_vs_mass_Ppip_thrown[w_range_num][q2_range_num];
TH2D_ptr theta_pip_lab_vs_mass_Ppim_thrown[w_range_num][q2_range_num];

// EC Sampling Fraction
TH2D_ptr EC_sampling_fraction[CUTS];
// EC Sampling Fraction

// Mom vs Beta
TH2D_ptr momvsbeta_hist[particle_num][charge_num][with_id_num];
// Mom vs Beta

// Delta T
TH2D_ptr delta_t_hist[particle_num][charge_num][with_id_num][2];
// Delta T

public:
Histogram(const std::string &output_file);
~Histogram();

void populate_eff_check_mPim(const std::shared_ptr<Reaction> &_e, double min_w, double max_w, double min_theta, double max_theta, double min_phi, double max_phi, short index_w, short index_theta, short index_phi);
// void populate_eff_check_exclusive(const std::shared_ptr<Reaction> &_e, double min, double max, short index_w);
void Fill_eff_ckeck_mPim(const std::shared_ptr<Reaction> &_e);
//void Fill_eff_ckeck_exclusive(const std::shared_ptr<Reaction> &_e);
void writeHist_eff_check_mPim();
// void writeHist_eff_check_mPim_after_MMSQ_cuts();
// void writeHist_eff_check_exclusive();
// void writeHist_eff_check_exclusive_MMSQ_cuts();



void Fill_histthreeD(const std::shared_ptr<Reaction> &_e);
void writeHists3D();
void Fill_hists4D_background(const std::shared_ptr<Reaction> &_e);
void writehists4D_background();


void Fill_histSevenD_prot(const std::shared_ptr<Reaction> &_e);
void Fill_histSevenD_thrown_prot(const std::shared_ptr<MCReaction> &_e);
void writeHists7D_prot();
void writeHists7D_thrown_prot();
void writeHists7D_pim();
void writeHists7D_thrown_pim();
void Fill_histSevenD_pim(const std::shared_ptr<Reaction> &_e);
void Fill_histSevenD_thrown_pim(const std::shared_ptr<MCReaction> &_e);
void Fill_histSevenD_pip(const std::shared_ptr<Reaction> &_e);
void Fill_histSevenD_thrown_pip(const std::shared_ptr<MCReaction> &_e);
void writeHists7D_pip();
void writeHists7D_thrown_pip();

// W and Q^2
void makeHists_sector();
void Fill_WvsQ2(const std::shared_ptr<Reaction> &_e);
void Fill_WvsQ2(const std::shared_ptr<MCReaction> &_e);
void Fill_WvsQ2_singlePi(const std::shared_ptr<Reaction> &_e);
void Fill_WvsQ2_Npip(const std::shared_ptr<Reaction> &_e);
void Fill_WvsQ2_twoPi(const std::shared_ptr<Reaction> &_e);
// void Fill_WvsQ2_twoPi(const std::shared_ptr<MCReaction>& _e);
void Fill_WvsQ2_twoPi_thrown(const std::shared_ptr<MCReaction> &_e);
void Write_WvsQ2();

void Fill_W_vs_Q2_all_sec();
void Fill_W_vs_Q2_thrown();
void Fill_inv_mass_hist();
// P and E
// ecectron cuts
void makeHists_electron_cuts();
void FillHists_electron_cuts(const std::shared_ptr<Branches12>& _d);
// void FillHists_electron_with_cuts(const std::shared_ptr<Branches12>& _d);

void Write_Electron_cuts();

// sampling Fraction
// void Fill_SF(const std::shared_ptr<Branches12>& _d, int part);
// void Write_SF();
void makeHists_x_mu();
void Fill_x_mu(const std::shared_ptr<Reaction> &_e);
void write_hist_x_mu();
void Fill_pi0(const std::shared_ptr<Reaction> &_e);

void makeHists_MomVsBeta();
void Fill_MomVsBeta(const std::shared_ptr<Branches12> &data, int part);
void Write_MomVsBeta();

// Delta T
void makeHists_deltat();
void Fill_deltat_pi(const std::shared_ptr<Branches12> &data,
                    const std::shared_ptr<Delta_T> &dt, int part);
void Fill_deltat_prot(const std::shared_ptr<Branches12> &data,
                      const std::shared_ptr<Delta_T> &dt, int part);

void Write_deltat();

// EC Sampling Fraction
//  void Fill_EC(double etot, double momentum);
// void Write_EC();

void makeHistTheta_pim_measured();
void populate_theta_pim_measured(const std::shared_ptr<Reaction> &_e, double min, double max, short index_theta_pim);
void Fill_theta_pim_measured(const std::shared_ptr<Reaction> &_e);
void write_hist_theta_pim_measured();
void makeHist_efficiency_check();
void Fill_efficiency_check_CD_rec(const std::shared_ptr<Reaction> &_e);
void Fill_efficiency_check_FD_rec(const std::shared_ptr<Reaction> &_e);
void Fill_efficiency_check_CD(const std::shared_ptr<Reaction> &_e);
void Fill_efficiency_check_FD(const std::shared_ptr<Reaction> &_e);
void write_hist_efficiency_check_CD();
void write_hist_efficiency_check_FD();
void Fill_efficiency_check_CD_thrown(const std::shared_ptr<MCReaction> &_mc_e,const std::shared_ptr<Reaction> &_e);
void Fill_efficiency_check_FD_thrown(const std::shared_ptr<MCReaction> &_mc_e,const std::shared_ptr<Reaction> &_e);
// void Fill_efficiency_check_missingPip_CD(const std::shared_ptr<Reaction> &_e);
// void Fill_efficiency_check_missingPip_FD(const std::shared_ptr<Reaction> &_e);
// void Fill_efficiency_check_missingProt_CD(const std::shared_ptr<Reaction> &_e);
// void Fill_efficiency_check_missingProt_FD(const std::shared_ptr<Reaction> &_e);
void makeHist_efficiency_check_1d();
void write_hist_efficiency_check_1d();
void Fill_efficiency_check_1d_hist_prot(const std::shared_ptr<Reaction> &_e);
void Fill_efficiency_check_1d_hist_pip(const std::shared_ptr<Reaction> &_e);
void Fill_efficiency_check_1d_hist_pim(const std::shared_ptr<Reaction> &_e);
//
void Write();
};

#endif
