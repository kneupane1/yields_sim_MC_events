
#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include <iostream>
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "branches.hpp"
#include "constants.hpp"
#include "physics.hpp"
#include"mom_corr.hpp"
class Reaction {
 protected:
  std::shared_ptr<Branches12> _data;

  double _beam_energy = 10.6;
  std::unique_ptr<TLorentzVector> _beam;
  std::unique_ptr<TLorentzVector> _elec;
  std::unique_ptr<TLorentzVector> _gamma;
  std::unique_ptr<TLorentzVector> _target;
  std::unique_ptr<TLorentzVector> _prot;
  std::unique_ptr<TLorentzVector> _pip;
  std::unique_ptr<TLorentzVector> _pim;
  std::unique_ptr<TLorentzVector> _other;
  std::unique_ptr<TLorentzVector> _neutron;
  std::unique_ptr<TLorentzVector> _boosted_gamma;
  std::unique_ptr<TLorentzVector> _boosted_prot;
  std::unique_ptr<TLorentzVector> _boosted_pip;
  std::unique_ptr<TLorentzVector> _boosted_pim;

  std::unique_ptr<TLorentzVector> _mom_corr_elec;
  std::unique_ptr<TLorentzVector> _mom_corr_pim_th;
  std::unique_ptr<TLorentzVector> _mom_corr_pim_ph;
  std::unique_ptr<TLorentzVector> _mom_corr_pim;
  std::unique_ptr<TLorentzVector> _mom_corr_pip_th;
  std::unique_ptr<TLorentzVector> _mom_corr_pip_ph;
  std::unique_ptr<TLorentzVector> _mom_corr_pip;
  std::unique_ptr<TLorentzVector> _mom_corr_prot_th;
  std::unique_ptr<TLorentzVector> _mom_corr_prot_ph;
  std::unique_ptr<TLorentzVector> _mom_corr_prot;
  std::unique_ptr<TLorentzVector> _Energy_loss_uncorr_pim;
  std::unique_ptr<TLorentzVector> _Energy_loss_uncorr_pip;
  std::unique_ptr<TLorentzVector> _Energy_loss_uncorr_prot;
  std::unique_ptr<TLorentzVector> _pim_tmt;
  std::unique_ptr<TLorentzVector> _pip_tmt;

  std::unique_ptr<TLorentzVector> _rotated_prot;
  std::unique_ptr<TLorentzVector> _rotated_pip;
  std::unique_ptr<TLorentzVector> _rotated_pim;
  std::unique_ptr<TLorentzVector> _rotated_pim_measured;

  TVector3 _prot_Vect3;
  TVector3 _pip_Vect3;
  TVector3 _pim_Vect3;

  // std::unique_ptr<TLorentzVector> _missingPim;
  std::unique_ptr<TLorentzVector> _boosted_pim_measured;

  bool _is_boosted = false;

  bool _mc = false;

  bool _hasE = false;
  bool _hasP = false;
  bool _hasPip = false;
  bool _hasPim = false;
  bool _hasOther = false;
  bool _hasNeutron = false;

  bool _is_FD = false;
  bool _is_CD = false;
  bool _is_lower_band = false;

  short _numProt = 0;
  short _numPip = 0;
  short _numPim = 0;
  short _numPos = 0;
  short _numNeg = 0;
  short _numNeutral = 0;
  short _numOther = 0;

  short _sector = -1;

  float _MM = NAN;
  float _MM2 = NAN;
  float _MM2_exclusive = NAN;
  float _excl_Energy = NAN;
  float _MM2_mPip = NAN;
  float _MM2_mProt = NAN;

  float _W = NAN;
  float _Q2 = NAN;
  float _elec_mom = NAN;
  float _elec_E = NAN;
  float _theta_e = NAN;

  float _W_after = NAN;

  float _P_elec = NAN;

  float _theta_star = NAN;
  float _phi_star = NAN;

  float _rec_pim_mom = NAN;
  float _rec_pim_theta = NAN;
  float _rec_pim_phi = NAN;

  float _elec_phi = NAN;
  float _x_mu_phi = NAN;
  float _diff_elec_x_mu_theta = NAN;
  float _diff_elec_x_mu_phi = NAN;

  float _beam_phi = NAN;
  float _diff_beam_x_mu_theta = NAN;
  float _diff_beam_x_mu_phi = NAN;

  float _prot_status = NAN;
  float _pip_status = NAN;
  float _pim_status = NAN;

  int _sectorPim = -1;
  int _sectorPip = -1;
  int _sectorProt = -1;

  float _inv_Ppip = NAN;
  float _inv_Ppim = NAN;
  float _inv_pip_pim = NAN;

  void SetElec();



  /// finished momentum corrections earlier

  double _elec_mom_corrected = NAN;

  double _cx = NAN;
  double _cy = NAN;
  double _cz = NAN;

  double _px_prime_elec = NAN;
  double _py_prime_elec = NAN;
  double _pz_prime_elec = NAN;

  double fe = NAN;
  double fpro = NAN;
  double fpip = NAN;
  double fpim = NAN;
  float _thetaDC_r1_Prot = NAN;;
  float _thetaDC_r1_Pip = NAN;;
  float _thetaDC_r1_Pim = NAN;;

  //
  static const int CD_SEC = 3;
  static const int FD_SEC = 6;

  // Now for prot mom corrections

  double _px_prime_prot_th = NAN;
  double _py_prime_prot_th = NAN;
  double _pz_prime_prot_th = NAN;
  double _E_prime_prot_th = NAN;

  double _px_prime_prot_ph = NAN;
  double _py_prime_prot_ph = NAN;
  double _pz_prime_prot_ph = NAN;

  double _px_prime_prot_mom = NAN;
  double _py_prime_prot_mom = NAN;
  double _pz_prime_prot_mom = NAN;

  double _px_prime_prot_E = NAN;
  double _py_prime_prot_E = NAN;
  double _pz_prime_prot_E = NAN;

  // float alpha_prot_mom_corr_FD[2] = {0.6, 0.9};
  // float alpha_prot_mom_corr_CD[5] = {1.0, 0.5, 0.95};

  double _prot_mom_prime = NAN;
  double _prot_mom = NAN;
  double _prot_mom_tmt = NAN;
  double _prot_theta_tmt = NAN;
  double _prot_phi_tmt = NAN;

  double _prot_mom_prime_2nd = NAN;
  double _prot_mom_2nd = NAN;
  double _prot_mom_uncorr = NAN;
  float _E_corr_val_prot = NAN;
  double _prot_theta_uncorr = NAN;
  float _prot_phi_uncorr = NAN;

  // static const int Prot_mom_bins_CD = 11;
  // static const int Prot_mom_bins_FD = 22;

  // float min_prot_mom_values_CD[Prot_mom_bins_CD] = {0, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.4, 1.55, 1.75};
  // float max_prot_mom_values_CD[Prot_mom_bins_CD] = {0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.4, 1.55, 1.75, 3.0};
  // float prot_mom_corr_CD[Prot_mom_bins_CD] = {
  //     0.027, 0.0048, -0.0048, -0.0144, -0.024, -0.0336, -0.049, -0.077, -0.088, -0.104, -0.13
  // };

  // float prot_mom_corr_CD[CD_SEC][Prot_mom_bins_CD] = {
  //     {-0.027, -0.024, -0.0336, -0.0432, -0.0528, -0.0624, -0.0918, -0.0918, -0.078, -0.078, -0.035},
  //     {0.015, 0.0144, 0.024, 0.0144, 0.024, 0.0144, 0.0054, 0.0162, 0.006, 0.018, 0.049},
  //     {0.009, 0.024, 0.0432, 0.0336, 0.0336, 0.0432, 0.0378, 0.0486, 0.054, 0.078, 0.119}};
  // float min_prot_mom_values_FD[Prot_mom_bins_FD] = {0,    0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65,
  //                                                   1.75, 1.85, 1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.6,  2.8,  3.1};
  // float max_prot_mom_values_FD[Prot_mom_bins_FD] = {0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45,
  //                                                   1.55, 1.65, 1.75, 1.9,  2.1,  2.5,  3.0,  5.0};

  // double prot_mom_corr_FD[2][Prot_mom_bins_FD] = {
  //     {0.027, 0.024, 0.0144, 0.0048, 0.0144, 0.0048, -0.007, -0.007, -0.007, -0.008, -0.008, -0.008, -0.01, -0.01,
  //      -0.03, -0.03, -0.03},
  //     {0.039, 0.0336, 0.024, 0.024, 0.024, 0.024, 0.021, 0.021, 0.021, 0.024, 0.024, 0.024, 0.03, 0.03, 0.03, 0.01,
  //     0.03}};  // first is for theta < 27 and second is for theta > 27 degrees,

  // FD corrections seperated by theta < 27 or > 27, sectors and mom bins
  // double prot_mom_corr_FD[2][FD_SEC][Prot_mom_bins_FD] = {
  //     {{0.0144, 0.0048, -0.0048, 0.0048, -0.0048, -0.007, -0.007, -0.007, -0.008, -0.008, -0.008,
  //       -0.008, -0.01,  -0.01,   -0.01,  -0.01,   -0.01,  -0.01,  -0.01,  -0.01,  -0.01,  -0.03},
  //      {0.0144, 0.0144, 0.0144, 0.0048, 0.0144, 0.007, 0.007, 0.007, 0.008, 0.008, 0.008,
  //       0.008,  0.01,   -0.01,  0.01,   -0.01,  -0.01, -0.01, -0.01, -0.01, -0.01, -0.01},
  //      {0.0144, 0.0144, 0.0144, 0.0144, 0.0144, 0.007, 0.007, 0.007, 0.008, -0.008, 0.008,
  //       -0.008, 0.01,   -0.01,  -0.01,  -0.01,  -0.01, -0.01, -0.01, -0.01, -0.01,  -0.01},
  //      {0.0144, 0.0144, 0.0144, 0.0048, 0.0048, 0.007, 0.007, -0.007, 0.008, -0.008, -0.008,
  //       0.008,  -0.01,  -0.01,  -0.01,  -0.01,  -0.01, -0.01, -0.01,  -0.01, -0.01,  -0.01},
  //      {0.0048, 0.0144, 0.0144, 0.0048, 0.0048, -0.007, -0.007, -0.007, -0.008, -0.008, -0.008,
  //       -0.008, -0.01,  -0.01,  -0.01,  -0.01,  -0.01,  -0.01,  -0.01,  -0.01,  -0.01,  -0.03},
  //      {0.0048, 0.0144, 0.0048, 0.0048, -0.0048, -0.007, -0.021, -0.021, -0.024, -0.008, -0.024,
  //       -0.024, -0.03,  -0.03,  -0.03,  -0.03,   -0.03,  -0.03,  -0.03,  -0.03,  -0.03,  -0.03}},

  //     {{-0.0048, 0.0048, 0.0144, -0.0048, -0.0048, -0.007, 0.007, -0.007, 0.008, 0.008, 0.008,
  //       0.008,   0.01,   0.01,   0.01,    0.01,    0.01,   0.03,  0.01,   0.03,  0.03,  0.05},
  //      {0.0048, 0.0048, 0.0048, 0.0048, 0.0048, 0.007, 0.007, 0.007, 0.008, 0.008, 0.008,
  //       0.008,  0.01,   0.01,   0.01,   0.03,   0.01,  0.01,  0.01,  0.03,  0.03,  0.05},
  //      {0.0144, 0.024, 0.0144, 0.024, 0.024, 0.021, 0.021, 0.021, 0.024, 0.024, 0.024,
  //       0.024,  0.03,  0.03,   0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  0.03,  0.05},
  //      {0.0048, 0.0048, 0.0144, 0.0144, 0.0144, 0.007, 0.007, 0.021, 0.008, 0.008, 0.024,
  //       0.024,  0.01,   0.01,   0.03,   0.01,   0.01,  0.03,  0.01,  0.03,  0.03,  0.03},
  //      {0.0144, 0.0144, 0.0144, 0.0144, 0.0048, 0.007, 0.007, 0.007, 0.008, 0.024, 0.008,
  //       0.024,  0.03,   0.03,   0.03,   0.03,   0.03,  0.03,  0.03,  0.05,  0.05,  0.07},
  //      {0.0048, 0.0048, -0.0048, 0.0048, -0.0048, 0.007, -0.007, 0.007, -0.008, 0.008, 0.008,
  //       0.008,  0.01,   0.01,    0.01,   0.01,    0.01,  0.01,   0.03,  0.01,   0.03,  0.05}}};

  // float prot_mom_corr_CD_2nd[Prot_mom_bins_CD] = {
  //     0.003, 0.0048, -0.0144, -0.0048, -0.0048, -0.0144, -0.021, -0.035, -0.056, -0.056, -0.09
  // };
  // double prot_mom_corr_FD_2nd[2][Prot_mom_bins_FD] = {
  //     {0.003, 0.0144, -0.0048, 0.0048, -0.0048, -0.0048, 0.007, 0.007, 0.007, 0.008, 0.008, -0.008, 0.03, 0.01,
  //     -0.01, -0.01, 0.01}, {0.003, 0.0048, 0.0048, 0.0048, 0.0048, 0.0048, 0.007, 0.007, 0.007, 0.008, 0.008, 0.008,
  //     0.01, 0.01, 0.01, 0.01, 0.01}};  // first is for theta < 27 and second is for theta > 27 degrees,

  // // double prot_mom_corr_sim_FD[Prot_mom_bins_FD][2] = {
  // // };

  static const int Prot_theta_bins = 10;
  float alpha_prot_theta_corr = 0.5;
  double _prot_theta = NAN;
  double _prot_theta_prime = NAN;

  // float min_prot_theta_values[Prot_theta_bins] = {0, 7, 10, 13, 16, 19, 22, 25, 27, 30, 33, 37, 42, 50, 60};
  // float max_prot_theta_values[Prot_theta_bins] = {7, 10, 13, 16, 19, 22, 25, 27, 30, 33, 37, 42, 50, 60, 180};
  // double prot_theta_corr[Prot_theta_bins] = {0.675, 0.225, 0.225, 0.225, 0.225, 0.225, 0.225, 0.075,
  //                                            0.15,  -0.15, -0.45, 0.15,  -0.45, -3.15, -3.15};
  // ,
  // double prot_theta_corr_sim[Prot_theta_bins] = {
  //     0.1875 ,0.1125, 0.1125, 0.1125, 0.1125, 0.1125, 0.1125, 0.1125, 0.075, 0.075, 0.075, 0.075, -0.225, -1.125,
  //     -1.125
  // // };

  static const int Prot_phi_bins = 11;
  float alpha_prot_phi_corr = 0.5;
  double _prot_phi = NAN;
  double _prot_phi_prime = NAN;
  // float min_prot_phi_values[Prot_phi_bins] = {0,   15,  30,  45,  60,  75,  90,  105, 120, 135, 150, 165,
  //                                             180, 195, 210, 225, 240, 255, 270, 285, 300, 315, 330, 345};
  // float max_prot_phi_values[Prot_phi_bins] = {15,  30,  45,  60,  75,  90,  105, 120, 135, 150, 165, 180,
  //                                             195, 210, 225, 240, 255, 270, 285, 300, 315, 330, 345, 360};
  // double prot_phi_corr[Prot_phi_bins] = {0.7,  1.1, 0.9, 0.7, 0.5, 0.5, 0.5, 0.5, 0.1, 0.3, -0.1, -0.1,
  //                                        -0.3, 0.1, 0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.5, 0.9, 0.9,  0.7};
  // // double prot_phi_corr_sim[Prot_phi_bins] = {0.075,- 0.075,-0.075, 0.075, 0.075, 0.075, -0.075, 0.075, 0.075,
  // -0.075,
  // // -0.075, 0.075, 0.075, -0.075, -0.075, 0.075, 0.075,0.075, -0.075, 0.075, 0.075, -0.075, -0.075, -0.075,};

  // Now for pip mom corrections

  double _px_prime_pip_th = NAN;
  double _py_prime_pip_th = NAN;
  double _pz_prime_pip_th = NAN;
  double _E_prime_pip_th = NAN;

  double _px_prime_pip_ph = NAN;
  double _py_prime_pip_ph = NAN;
  double _pz_prime_pip_ph = NAN;

  double _px_prime_pip_mom = NAN;
  double _py_prime_pip_mom = NAN;
  double _pz_prime_pip_mom = NAN;

  double _px_prime_pip_E = NAN;
  double _py_prime_pip_E = NAN;
  double _pz_prime_pip_E = NAN;

  double _px_prime_pip_E_tmt = NAN;
  double _py_prime_pip_E_tmt = NAN;
  double _pz_prime_pip_E_tmt = NAN;

  // float alpha_pip_mom_corr_FD[2] = {0.5, 0.7};
  // float alpha_pip_mom_corr_CD[3] = {0.9, 0.45, 0.9};

  double _pip_mom = NAN;
  double _pip_mom_prime = NAN;

  double _pip_mom_2nd = NAN;
  double _pip_mom_prime_2nd = NAN;

  double _pip_mom_tmt = NAN;
  double _pip_theta_tmt = NAN;
  double _pip_phi_tmt = NAN;

  double _pip_mom_tmt2 = NAN;
  double _pip_mom_uncorr = NAN;
  double _pip_theta_uncorr = NAN;
  double _pip_phi_uncorr = NAN;

  float _E_corr_val_pip = NAN;
  float _E_corr_val_pip_th = NAN;
  float _E_corr_val_pip2 = NAN;

  // static const int Pip_mom_bins_FD = 18;
  // static const int Pip_mom_bins_CD = 15;

  // float min_pip_mom_values_CD[Pip_mom_bins_CD] = {0,   0.25, 0.3, 0.35, 0.4,  0.45, 0.5, 0.55,
  //                                                 0.6, 0.65, 0.7, 0.75, 0.85, 1.0,  1.2};
  // float max_pip_mom_values_CD[Pip_mom_bins_CD] = {0.25, 0.3, 0.35, 0.4,  0.45, 0.5, 0.55, 0.6,
  //                                                 0.65, 0.7, 0.75, 0.85, 1.0,  1.2, 3.0};
  // float pip_mom_corr_CD[Pip_mom_bins_CD] = {0.045,  0.015,   0.003,   -0.003,  -0.003,  -0.009, -0.009, -0.015,
  // -0.021,
  //                                           -0.027, -0.0336, -0.0432, -0.0528, -0.0624, -0.072, -0.105, -0.13};

  // float pip_mom_corr_CD[CD_SEC][Pip_mom_bins_CD] = {
  //     {-0.0084, -0.014, -0.0084, -0.014, -0.014, -0.018, -0.0252, -0.0252, -0.0324, -0.0324, -0.0396, -0.044, -0.052,
  //      -0.065, -0.077},
  //     {0.0028, 0.0084, 0.0084, 0.014, 0.014, 0.018, 0.018, 0.018, 0.018, 0.018, 0.018, 0.02, 0.02, 0.025, 0.021},
  //     {-0.0028, -0.0028, 0.0028, 0.0028, 0.0084, 0.0108, 0.0108, 0.018, 0.0108, 0.018, 0.018, 0.028, 0.036, 0.045,
  //      0.077}};

  // float min_pip_mom_values_FD[2][Pip_mom_bins_FD] = {
  //     {0, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.9, 2.05, 2.25, 2.45, 2.7, 3.1},
  //     {0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.95, 2.1, 2.25}};
  // float max_pip_mom_values_FD[2][Pip_mom_bins_FD] = {
  //     {0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.9, 2.05, 2.25, 2.45, 2.7, 3.1, 5.0},
  //     {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.95, 2.1, 2.25, 5.0}};
  // double pip_mom_corr_FD[2][Pip_mom_bins_FD] = {
  //     {0.021, 0.024, 0.0144, -0.0048, -0.0048, 0.0048, -0.007, -0.007, 0.007, 0.008, 0.008, 0.008, 0.01, -0.01,
  //     -0.01,
  //      0.01, 0.01},
  //     {0.025, 0.021, 0.027, 0.027, 0.027, 0.027, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.036,
  //     0.028,
  //      0.036}};  // first is for theta < 27 and second is for theta > 27 degrees,
  //                // last one for >27 is adjusted to second last one
  // float pip_mom_corr_CD_2nd[Pip_mom_bins_CD] = {0.033, 0.003,   0.009,   0.009,   0.009,   0.009,   0.021, 0.009,
  // 0.003,
  //                                               0.009, -0.0048, -0.0048, -0.0048, -0.0048, -0.0336, 0.007, -0.09};

  // double pip_mom_corr_FD[2][FD_SEC][Pip_mom_bins_FD] = {
  //     {{-0.0144, -0.0144, -0.0048, -0.0048, -0.0048, -0.005, -0.005, -0.005, -0.006, -0.006, -0.006, -0.007, -0.007,
  //       -0.007, 0.007, 0.007, 0.007, 0.021},
  //      {0.0048, 0.0048, 0.0048, 0.0048, 0.0048, -0.005, 0.005, 0.005, 0.006, 0.006, 0.006, 0.007, 0.007, 0.007, 0.007,
  //       0.007, 0.007, 0.021},
  //      {-0.0048, 0.0048, 0.0048, 0.0048, 0.0048, 0.015, 0.015, 0.005, 0.006, 0.018, 0.006, 0.007, 0.007, 0.007, 0.021,
  //       0.021, 0.021, 0.021},
  //      {0.0048, 0.0048, 0.0048, 0.0048, 0.0048, 0.005, 0.005, 0.005, 0.018, 0.006, 0.006, 0.007, 0.007, 0.007, 0.007,
  //       0.007, 0.007, 0.007},
  //      {0.0048, 0.0048, 0.0048, 0.0048, 0.0048, -0.005, -0.005, -0.005, -0.006, 0.006, 0.006, -0.007, 0.007, -0.007,
  //       0.007, 0.007, 0.007, 0.007},
  //      {-0.0048, -0.0048, -0.0048, 0.0048, -0.0048, -0.005, -0.005, -0.005, -0.006, -0.006, -0.006, -0.007, -0.007,
  //       -0.007, -0.007, 0.007, 0.007, 0.007}},

  //     {{-0.015, -0.009, -0.0048, -0.0048, -0.0048, 0.0048, 0.0048, -0.005, 0.005, 0.005, 0.006, 0.006, 0.018, 0.018,
  //       0.021, 0.021, 0.021, 0.021},
  //      {-0.021, -0.015, -0.0144, -0.0144, -0.0048, -0.0048, -0.0048, 0.005, 0.015, 0.015, 0.006, 0.018, 0.006, 0.006,
  //       0.021, 0.021, 0.021, 0.021},
  //      {-0.009, -0.003, -0.0048, 0.0048, 0.0048, 0.0048, 0.0144, 0.005, 0.015, 0.015, 0.018, 0.018, 0.018, 0.03, 0.021,
  //       0.021, 0.035, 0.035},
  //      {-0.009, 0.009, -0.0144, -0.0048, 0.0048, 0.0048, 0.0048, 0.015, 0.015, 0.015, 0.006, 0.006, 0.018, 0.018, 0.021,
  //       0.021, 0.021, 0.021},
  //      {0.003, 0.009, 0.0048, 0.0048, 0.0048, 0.0048, 0.0144, 0.015, 0.005, 0.025, 0.018, 0.018, 0.018, 0.03, 0.021,
  //       0.021, 0.021, 0.035},
  //      {-0.015, -0.003, 0.0048, 0.0048, 0.0048, 0.0048, 0.0144, 0.005, 0.005, 0.005, 0.006, 0.006, 0.006, 0.018, 0.007,
  //       0.021, 0.021, 0.035}}};

  static const int Pip_theta_bins = 10;
  float alpha_pip_theta_corr = 0.5;
  double _pip_theta = NAN;
  double _pip_theta_prime = NAN;

  // float min_pip_theta_values[Pip_theta_bins] = {0, 19.5, 24, 28, 33, 38, 45, 52, 60, 75};
  // float max_pip_theta_values[Pip_theta_bins] = {19.5, 24, 28, 33, 38, 45, 52, 60, 75, 180};
  // double pip_theta_corr[Pip_theta_bins] = {0.075, 0.075, 0.075, -0.15, -0.45, 1.0, -0.2, -1.2, -1.8, -5.4};
  // double pip_theta_corr_sim[Pip_theta_bins] = {0.1875, 0.1875, 0.1875, 0.375, 0.225, 0.3, 0.3, 0.375, 0.6, 0.2};

  static const int Pip_phi_bins = 11;
  float alpha_pip_phi_corr = 1.0;
  double _pip_phi = NAN;
  double _pip_phi_prime = NAN;
  // float min_pip_phi_values[Pip_phi_bins] = {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300};
  // float max_pip_phi_values[Pip_phi_bins] = {30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 360};
  // double pip_phi_corr[Pip_phi_bins] = {0.9, 0.7, 0.5, 0.3, 0.1, -0.1, -0.1, 0.1, 0.3, 0.3, 0.7};
  // double pip_phi_corr_sim[Pip_phi_bins] = {0.1, -0.1, 0.1, -0.1, 0.1, -0.1, 0.1, 0.1, 0.1, -0.1, -0.1};

  // for pim mom corrections ...............

  double _px_prime_pim_th = NAN;
  double _py_prime_pim_th = NAN;
  double _pz_prime_pim_th = NAN;
  double _E_prime_pim_th = NAN;

  double _px_prime_pim_ph = NAN;
  double _py_prime_pim_ph = NAN;
  double _pz_prime_pim_ph = NAN;

  double _px_prime_pim_mom = NAN;
  double _py_prime_pim_mom = NAN;
  double _pz_prime_pim_mom = NAN;

  double _px_prime_pim_E = NAN;
  double _py_prime_pim_E = NAN;
  double _pz_prime_pim_E = NAN;

  double _px_prime_pim_E_tmt = NAN;
  double _py_prime_pim_E_tmt = NAN;
  double _pz_prime_pim_E_tmt = NAN;

  // float alpha_pim_mom_corr_FD[2] = {0.3, 0.6};
  // float alpha_pim_mom_corr_CD[4] = {0.85, 1.0, 0.4};

  double _pim_mom = NAN;
  double _pim_mom_2nd = NAN;
  double _pim_mom_tmt = NAN;
  double _pim_theta_tmt = NAN;
  double _pim_phi_tmt = NAN;
  double _pim_mom_tmt2 = NAN;
  double _pim_mom_prime = NAN;
  double _pim_mom_prime_2nd = NAN;
  double _pim_mom_uncorr = NAN;
  double _pim_theta_uncorr = NAN;
  double _pim_phi_uncorr = NAN;

  float _E_corr_val_pim = NAN;
  float _E_corr_val_pim_th = NAN;
  float _E_corr_val_pim2 = NAN;

  // static const int Pim_mom_bins_FD = 21;
  // static const int Pim_mom_bins_CD = 17;

  // float min_pim_mom_values_CD[Pim_mom_bins_CD] = {0,   0.2,  0.25, 0.3,  0.35, 0.4,  0.45, 0.5, 0.55,
  //                                                 0.6, 0.65, 0.7,  0.75, 0.85, 0.95, 1.1,  1.3};
  // float max_pim_mom_values_CD[Pim_mom_bins_CD] = {0.2,  0.25, 0.3,  0.35, 0.4,  0.45, 0.5, 0.55, 0.6,
  //                                                 0.65, 0.7,  0.75, 0.85, 0.95, 1.1,  1.3, 3.0};
  // float pim_mom_corr_CD[Pim_mom_bins_CD] = {0.039,   0.009,   0.003,   -0.003,  -0.0048, -0.0048,
  //                                           -0.0048, -0.0048, -0.0048, -0.0144, -0.0144, -0.0144,
  //                                           -0.0144, -0.024,  -0.024,  -0.035,  -0.05};

  // float pim_mom_corr_CD[CD_SEC][Pim_mom_bins_CD] = {
  //     {0.033, 0.009, -0.003, -0.009, -0.0048, -0.0048, -0.0048, -0.0048, -0.0048, -0.0144, -0.0144, -0.0144, -0.0144,
  //      -0.0144, -0.024, -0.021, -0.03},
  //     {0.039, 0.015, 0.003, 0.003, 0.0048, 0.0048, 0.0144, 0.0144, 0.0144, 0.0048, 0.0048, 0.0144, 0.0144, 0.0144,
  //      0.0144, 0.021, 0.03},
  //     {0.045, 0.009, 0.003, -0.003, -0.0048, -0.0144, -0.0144, -0.0144, -0.024, -0.0336, -0.0336, -0.0336, -0.0432,
  //      -0.0528, -0.0624, -0.077, -0.09}};

  // // binning in measured pim
  // float pim_mom_corr_CD[CD_SEC][Pim_mom_bins_CD] = {
  //     {0.0032, -0.0096, -0.0096, -0.0032, -0.0032, -0.0032, -0.0036, 0.0036, 0.0036, 0.0108, 0.0036, 0.0108, 0.012,
  //      0.02, 0.028, 0.055, 0.105},
  //     {-0.0288, -0.0032, 0.0032, 0.0032, 0.0032, 0.0096, 0.0108, 0.018, 0.018, 0.0252, 0.0252, 0.0324, 0.036, 0.044,
  //      0.06, 0.075, 0.119},
  //     {-0.0096, -0.0032, -0.0032, -0.0032, -0.0096, -0.0096, -0.0036, -0.0036, -0.0036, -0.0036, -0.0108, -0.0108,
  //      -0.012, -0.012, -0.012, -0.005, 0.021}};

  // float min_pim_mom_values_FD[2][Pim_mom_bins_FD] = {{0,    0.85, 0.95, 1.1,  1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85,
  //                                                     1.95, 2.05, 2.15, 2.25, 2.35, 2.45, 2.6,  2.75, 3.0,  3.4},
  //                                                    {0.5,  0.7,  0.8,  0.9, 1.0, 1.08, 1.164, 1.25, 1.32, 1.4, 1.48,
  //                                                     1.56, 1.64, 1.72, 1.8, 1.9, 2.05, 2.2,   2.4,  2.6,  2.8}};
  // float max_pim_mom_values_FD[2][Pim_mom_bins_FD] = {{0.85, 0.95, 1.1,  1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95,
  //                                                     2.05, 2.15, 2.25, 2.35, 2.45, 2.6,  2.75, 3.0,  3.4,  5.0},
  //                                                    {0.7,  0.8,  0.9, 1.0, 1.08, 1.164, 1.25, 1.32, 1.4, 1.48, 1.56,
  //                                                     1.64, 1.72, 1.8, 1.9, 2.05, 2.2,   2.4,  2.6,  2.8, 5.0}};
  // double pim_mom_corr_FD[2][Pim_mom_bins_FD] = {
  //     {0.087, 0.024, 0.0144, 0.0048, 0.0144, 0.0048, 0.007, 0.007, 0.007, -0.008, 0.008, 0.008, -0.01, -0.01, -0.01,
  //      -0.01, 0.01},
  //     {0.069, 0.0432, 0.024, 0.024, 0.024, 0.024, 0.021, 0.021, 0.021, 0.024, 0.008, 0.008, 0.01, 0.01, 0.01, 0.01,
  //      0.01}};  // first is for theta < 27 and second is for theta > 27 degrees, last bin for > 27 is
  //               // adjusted to second last bin

  // double pim_mom_corr_FD[2][FD_SEC][Pim_mom_bins_FD] = {
  //     {{0.024,  0.024,  0.0144, 0.0162, 0.0162, 0.0162, 0.006, 0.018, 0.006, 0.018, 0.007,
  //       -0.007, -0.007, -0.007, -0.007, -0.007, -0.007, 0.007, 0.007, 0.007, 0.007},
  //      {0.0144, 0.024, 0.0144, 0.0162, 0.0162, 0.0162, 0.006, 0.006, 0.018, 0.018, 0.021,
  //       0.021,  0.007, 0.007,  0.007,  0.007,  0.021,  0.021, 0.007, 0.007, 0.007},
  //      {0.0144, 0.0144, 0.024, 0.0162, 0.027, 0.0162, 0.018, 0.018, 0.018, 0.018, 0.021,
  //       0.007,  0.007,  0.007, 0.007,  0.007, 0.007,  0.007, 0.021, 0.021, 0.021},
  //      {0.0144, 0.0144, 0.0144, 0.0162, 0.0054, 0.0054, 0.006, 0.006, 0.006, 0.006, 0.007,
  //       0.007,  0.007,  0.007,  0.007,  0.007,  0.007,  0.007, 0.007, 0.007, 0.007},
  //      {0.0048, 0.024, 0.0144, 0.0054, 0.0162, 0.0054, 0.006, -0.006, 0.006, -0.006, -0.007,
  //       0.007,  0.007, 0.021,  0.007,  0.007,  0.007,  0.007, 0.007,  0.007, 0.007},
  //      {0.0048, 0.024,  0.0048, 0.0054, 0.0054, 0.0054, -0.006, -0.006, 0.006, 0.006, -0.007,
  //       -0.007, -0.007, -0.007, -0.007, -0.007, 0.007,  -0.007, -0.007, 0.007, 0.007}},

  //     {{-0.0144, 0.024, 0.0144, 0.024,  0.024,  0.0144, 0.015,  0.005,  0.015, 0.015, -0.006,
  //       0.006,   0.006, 0.006,  -0.006, -0.006, -0.006, -0.006, -0.006, 0.006, -0.018},
  //      {-0.0144, 0.0048, 0.0048, 0.0144, 0.0144, 0.0144, 0.015, 0.015, 0.005, 0.015, 0.018,
  //       0.018,   0.018,  0.018,  0.006,  0.006,  0.018,  0.018, 0.018, 0.006, 0.006},
  //      {-0.024, 0.0144, 0.024, 0.0144, 0.0144, 0.0048, 0.015, 0.015, 0.015, 0.015, 0.018,
  //       0.018,  0.006,  0.018, 0.018,  0.006,  0.006,  0.006, 0.018, 0.018, 0.006},
  //      {-0.0144, 0.0144, 0.0144, 0.0144, -0.0048, 0.0048, 0.015,  0.005, 0.015,  0.005, 0.006,
  //       0.006,   0.006,  0.006,  0.006,  -0.006,  -0.006, -0.006, 0.006, -0.006, 0.018},
  //      {-0.02, -0.0048, -0.0048, 0.0048, -0.0048, 0.0048, 0.015,  -0.005, -0.005, -0.005, 0.018,
  //       0.006, -0.006,  0.006,   -0.006, -0.006,  -0.006, -0.006, -0.006, -0.006, 0.006},
  //      {-0.0048, 0.0048, 0.0144, 0.0048, -0.0048, -0.0048, -0.005, -0.005, 0.005, -0.005, 0.006,
  //       -0.006,  -0.006, -0.006, -0.006, -0.006,  -0.006,  -0.006, -0.006, 0.006, -0.006}}};

  // float pim_mom_corr_CD_2nd[Pim_mom_bins_CD] = {
  //     0.039,  0.009,  0.003,  0.003,   0.0048,  0.0048,  0.0048, 0.0048, 0.0048,
  //     0.0048, 0.0048, 0.0048, -0.0048, -0.0048, -0.0048, -0.007, -0.03,
  // };

  // double pim_mom_corr_FD_2nd[2][Pim_mom_bins_FD] = {{0.069, 0.0048, -0.0048, -0.0048, -0.0144, -0.0144, -0.007,
  // -0.021,
  //                                                    -0.021, -0.008, -0.008, -0.008, -0.01, -0.01, -0.01, -0.01,
  //                                                    -0.01},
  //                                                   {}};

  // double pim_mom_corr_sim_FD[Pim_mom_bins_FD][2] = {
  // };

  static const int Pim_theta_bins = 10;
  float alpha_pim_theta_corr = 0.5;
  double _pim_theta = NAN;
  double _pim_theta_prime = NAN;

  float min_pim_theta_values[Pim_theta_bins] = {0, 19.5, 24, 28, 33, 38, 45, 52, 60, 75};
  float max_pim_theta_values[Pim_theta_bins] = {19.5, 24, 28, 33, 38, 45, 52, 60, 75, 180};
  double pim_theta_corr[Pim_theta_bins] = {
      0.15, 0.05, 0.05, 0.1, -0.3, 0.75, 0.15, -0.2, -0.3, -2.1,
  };
  double pim_theta_corr_sim[Pim_theta_bins] = {
      0.075, 0.075, 0.075, 0.15, 0.15, 0.225, -0.075, 0.1, 0.15, -0.15,
  };

  static const int Pim_phi_bins = 11;
  float alpha_pim_phi_corr = 0.5;
  double _pim_phi = NAN;
  double _pim_phi_prime = NAN;
  float min_pim_phi_values[Pim_phi_bins] = {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300};
  float max_pim_phi_values[Pim_phi_bins] = {30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 360};
  double pim_phi_corr[Pim_phi_bins] = {
      0.5, 0.1, 0.1, 0.3, -0.1, -0.3, -0.3, -0.1, -0.1, -0.1, 0.3,
  };
  double pim_phi_corr_sim[Pim_phi_bins] = {
      -0.1, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1,
  };

 public:
  Reaction(){};
  Reaction(const std::shared_ptr<Branches12> &data, float beam_energy);
  ~Reaction();
  inline float weight() {
    // return _data->mc_weight();
    return 1.0;
  }
  // Check lists when you swich from mc to exp or vice-versa
  // 1. inline weight function above
  // 2. gamma, _w, _q2 and dpp function in electron four vector set up at reaction.cpp because of momentum corrections
  // for elec included only for exp data
  // 3. turn on the SetMomCorrElec() function on clas12_yields.hpp
  // 4. clas12_yields: auto data = std::make_shared<Branches12>(_chain, true);  turn off true for data
  // 5. from if (data->mc_npart() < 1) to all particle set up im mc events.
  // 6. all mc bank related (generated) output parameters will not work in exp data

  // momentum correction
  void SetMomCorrElec();
  double dpp(float px, float py, float pz, int sec_mom_corr, int ivec);
  double dppC(float px, float py, float pz, int sec, int ivec);
  double Corr_elec_mom();

  inline bool mc() { return _mc; }
  void SetProton(int i);
  void SetPip(int i);
  void SetPim(int i);
  void SetOther(int i);
  void SetNeutron(int i);

  float rec_pim_px();
  float rec_pim_py();
  float rec_pim_pz();
  float rec_pim_E();
  float rec_pim_P();
  float rec_pim_mm2();

  float Diff_elec_x_mu_theta();
  float Diff_elec_x_mu_phi();

  float Diff_beam_x_mu_theta();
  float Diff_beam_x_mu_phi();

  float pim_px();
  float pim_py();
  float pim_pz();
  float pim_E();
  float pim_P();

  float pip_px();
  float pip_py();
  float pip_pz();
  float pip_E();

  float prot_px();
  float prot_py();
  float prot_pz();
  float prot_E();

  float elec_px();
  float elec_py();
  float elec_pz();
  inline float elec_mom() { return _elec_mom; }
  inline float elec_En() { return _elec_E; }
  inline float Theta_Elec() { return _theta_e; }

  float beam_px();
  float beam_py();
  float beam_pz();
  float beam_E();

  float target_px();
  float target_py();
  float target_pz();
  float target_E();

  float pim_momentum_corrected();
  float pim_theta_corrected();
  float pim_Phi_corrected();

  float pip_momentum_corrected();
  float pip_theta_corrected();
  float pip_Phi_corrected();

  float prot_momentum_corrected();
  float prot_theta_corrected();
  float prot_Phi_corrected();

  // missingPim
  float pim_momentum();
  float pim_theta_lab();
  float pim_Phi_lab();
  float pim_momentum_measured();
  float pim_theta_lab_measured();
  float pim_Phi_lab_measured();

  float pim_theta_cm();
  float pim_Phi_cm();
  float pim_momentum_cm();
  float pim_theta_cm_measured();
  float pim_Phi_cm_measured();
  float pim_momentum_cm_measured();

  // missingPip
  float pip_momentum();
  float pip_theta_lab();
  float pip_Phi_lab();
  float pip_momentum_measured();
  float pip_theta_lab_measured();
  float pip_Phi_lab_measured();

  // missingProt
  float prot_momentum();
  float prot_theta_lab();
  float prot_Phi_lab();
  float prot_momentum_measured();
  float prot_theta_lab_measured();
  float prot_Phi_lab_measured();

  void boost();

  inline float Theta_star() { return _theta_star; }
  inline float Phi_star() { return _phi_star; }

  void CalcMissMass();
  float MM();
  float MM2();
  float MM2_exclusive();
  float Energy_excl();
  float MM2_mPip();
  float MM2_mProt();
  float MM2_mPim_corr();
  float MM2_mPip_corr();
  float MM2_mProt_corr();

  float w_hadron();
  float w_difference();
  float w_hadron_corr();
  float w_difference_corr();

  float inv_Ppip();
  float inv_Ppim();
  float inv_Pippim();

  void invMassPpip();
  void invMassPpim();
  void invMasspippim();

  virtual std::string CsvHeader();
  virtual std::string ReacToCsv();

  inline float thetaDCr1Prot(){ return _thetaDC_r1_Prot;}
  inline float thetaDCr1Pip() { return _thetaDC_r1_Pip; }
  inline float thetaDCr1Pim() { return _thetaDC_r1_Pim; }

  inline float W() { return _W; }
  inline float Q2() { return _Q2; }

  inline float W_after() { return _W_after; }

  float_t scalar_triple_product();

  inline short sec() { return _data->dc_sec(0); }
  inline short pimSec() { return _sectorPim; }
  inline short pipSec() { return _sectorPip; }
  inline short protSec() { return _sectorProt; }

  inline int det() { return abs(_data->status(0) / 1000); }

  inline bool Inclusive() {
    bool _channelIncl = true;
    _channelIncl &= (_hasE);
    return _channelIncl;
  }

  inline bool TwoPion_missingPim() {
    bool _channelTwoPi = true;
    _channelTwoPi &= ((_numProt == 1 && _numPip == 1) && (_hasE && _hasP && _hasPip));
    return _channelTwoPi;
  }

  inline bool TwoPion_exclusive() {
    bool _channelTwoPi_excl = true;

    _channelTwoPi_excl &= ((_numProt == 1 && _numPip == 1 && _numPim == 1) &&
                           (_hasE && _hasP && _hasPip && _hasPim /*&& !_hasNeutron && !_hasOther*/));
    return _channelTwoPi_excl;
  }
  inline bool TwoPion_missingPip() {
    bool _channelTwoPi_mpip = true;

    _channelTwoPi_mpip &=
        ((_numProt == 1 && _numPim == 1) && (_hasE && _hasP && _hasPim /*&&!_hasPip && !_hasNeutron && !_hasOther*/));
    return _channelTwoPi_mpip;
  }
  inline bool TwoPion_missingProt() {
    bool _channelTwoPi_mprot = true;
    _channelTwoPi_mprot &=
        ((_numPip == 1 && _numPim == 1) && (_hasE && _hasPip && _hasPim /*&&!_hasP  && !_hasOther*/));
    return _channelTwoPi_mprot;
  }

  const TLorentzVector &e_mu() { return *_beam; }
  const TLorentzVector &e_mu_prime() { return *_elec; }
  const TLorentzVector &gamma() { return *_gamma; }
};

class MCReaction : public Reaction {
 private:
  float _weight_mc = NAN;
  float _W_mc = NAN;
  float _Q2_mc = NAN;

  std::unique_ptr<TLorentzVector> _elec_mc;
  std::unique_ptr<TLorentzVector> _gamma_mc;
  std::unique_ptr<TLorentzVector> _prot_mc;
  std::unique_ptr<TLorentzVector> _pip_mc;
  std::unique_ptr<TLorentzVector> _pim_mc;
  std::unique_ptr<TLorentzVector> _other_mc;

  float _MM2_exclusive_mc = NAN;
  float _excl_Energy_mc = NAN;

  float _rec_x_mu_mom_mc = NAN;
  float _rec_x_mu_theta_mc = NAN;
  float _x_mu_phi_mc = NAN;

  float _beam_phi_mc = NAN;
  float _elec_phi_mc = NAN;
  float _diff_elec_x_mu_theta_mc = NAN;
  float _diff_elec_x_mu_phi_mc = NAN;

  float _diff_beam_x_mu_theta_mc = NAN;
  float _diff_beam_x_mu_phi_mc = NAN;

  float _elec_mom_mc = NAN;
  float _elec_E_mc = NAN;
  float _theta_e_mc = NAN;

 public:
  void SetMCProton(int i);
  void SetMCPip(int i);
  void SetMCPim(int i);
  void SetMCOther(int i);

  MCReaction(const std::shared_ptr<Branches12> &data, float beam_energy);
  void SetMCElec();
  inline float weight() { return _data->mc_weight(); }

  inline float elec_mom_mc() { return _elec_mom_mc; }
  inline float elec_En_mc() { return _elec_E_mc; }
  inline float Theta_Elec_mc() { return _theta_e_mc; }
  inline float W_mc() { return _W_mc; }
  inline float Q2_mc() { return _Q2_mc; }

  float pim_mom_mc_gen();
  float pip_mom_mc_gen();
  float prot_mom_mc_gen();

  float pim_theta_mc_gen();
  float pip_theta_mc_gen();
  float prot_theta_mc_gen();

  float pim_phi_mc_gen();
  float pip_phi_mc_gen();
  float prot_phi_mc_gen();

  void CalcMissMass_mc();

  float Diff_elec_x_mu_theta_mc();
  float Diff_elec_x_mu_phi_mc();
  float Diff_beam_x_mu_theta_mc();
  float Diff_beam_x_mu_phi_mc();
  float MM2_exclusive_mc();
  float Energy_excl_mc();
  float x_mu_momentum_mc();
  float x_mu_theta_lab_mc();
  float x_mu_Phi_lab_mc();

  std::string CsvHeader();
  std::string ReacToCsv();
};

#endif
