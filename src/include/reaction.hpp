
#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include <TRandom3.h>
#include <cmath>
#include <iostream>
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "branches.hpp"
#include "constants.hpp"
#include "mom_corr.hpp"
#include "physics.hpp"

class Reaction {
 protected:
  std::shared_ptr<Branches12> _data;

  double _beam_energy = 10.6041;
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

  std::unique_ptr<TLorentzVector> _elecUnSmear;
  std::unique_ptr<TLorentzVector> _protUnSmear;
  std::unique_ptr<TLorentzVector> _pipUnSmear;
  std::unique_ptr<TLorentzVector> _pimUnSmear;

  std::unique_ptr<TLorentzVector> _elecSmear;
  std::unique_ptr<TLorentzVector> _protSmear;
  std::unique_ptr<TLorentzVector> _pipSmear;
  std::unique_ptr<TLorentzVector> _pimSmear;

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

  float _MM_mPim = NAN;
  float _MM2_mPim = NAN;
  float _MM2_exclusive = NAN;
  float _excl_Energy = NAN;
  float _MM2_mPip = NAN;
  float _MM2_mProt = NAN;

  float _W = NAN;
  float _Q2 = NAN;
  float _elec_mom = NAN;
  float _elec_E = NAN;
  float _theta_e = NAN;
  float _phi_elec = NAN;

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
  int _sectorElec = -1;

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
  float _thetaDC_r1_Prot = NAN;
  ;
  float _thetaDC_r1_Pip = NAN;
  ;
  float _thetaDC_r1_Pim = NAN;
  ;

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

  static const int Prot_theta_bins = 10;
  float alpha_prot_theta_corr = 0.5;
  double _prot_theta = NAN;
  double _prot_theta_prime = NAN;

  static const int Prot_phi_bins = 11;
  float alpha_prot_phi_corr = 0.5;
  double _prot_phi = NAN;
  double _prot_phi_prime = NAN;

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

  static const int Pip_theta_bins = 10;
  float alpha_pip_theta_corr = 0.5;
  double _pip_theta = NAN;
  double _pip_theta_prime = NAN;

  static const int Pip_phi_bins = 11;
  float alpha_pip_phi_corr = 1.0;
  double _pip_phi = NAN;
  double _pip_phi_prime = NAN;

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
    return _data->mc_weight();
    // return 1.0;
  }
  // Check lists when you swich from mc to exp or vice-versa
  // 1. inline weight function above
  // 2. gamma, _w, _q2 and dpp function in electron four vector set up at reaction.cpp because of momentum corrections
  // for elec included only for exp data
  // 3. turn on the SetMomCorrElec() function on clas12_yields.hpp
  // 4. clas12_yields: auto data = std::make_shared<Branches12>(_chain, true);  turn off true for data
  // 5. from if (data->mc_npart() < 1) to all particle set up im mc events.
  // 6. all mc bank related (generated) output parameters will not work in exp data

  /// smearing fx's function
  void SmearingFunc(int part_id, int status_part, double p, double theta, double phi, double &pNew, double &thetaNew,
                    double &phiNew) {
    // Constants
    const double pS1 = 0.0184291 - 0.0110083 * theta + 0.00227667 * pow(theta, 2) - 0.000140152 * pow(theta, 3) +
                       3.07424e-6 * pow(theta, 4);
    const double pR = 0.02 * sqrt(pow(pS1 * p, 2) + pow(0.02 * theta, 2));
    const double thetaR = 2.5 * sqrt(pow((0.004 * theta + 0.1) * (pow(p, 2) + 0.13957 * 0.13957) / pow(p, 2), 2));
    const double phiS1 = 0.85 - 0.015 * theta;
    const double phiS2 = 0.17 - 0.003 * theta;
    const double phiR = 3.5 * sqrt(pow(phiS1 * sqrt(pow(p, 2) + 0.13957 * 0.13957) / pow(p, 2), 2) + pow(phiS2, 2));

    // Generate new values
    if (part_id == ELECTRON) {
      phiNew = phi + 0.38 * phiR * gRandom->Gaus(0, 1);
      thetaNew = theta + 0.38 * thetaR * gRandom->Gaus(0, 1);
      pNew = p + 0.38 * pR * gRandom->Gaus(0, 1) * p;

    } else if (part_id == PROTON) {
      double fact_cd = 0;
      double fact_fd = 0;
      if (status_part > 4000) {
        fact_cd = (0.008622 * pow(p, 2) + (0.034953) * p + 0.160969);

        phiNew = phi + 1 / (fact_cd)*phiR * gRandom->Gaus(0, 1);
        thetaNew = theta + 1 / (fact_cd)*thetaR * gRandom->Gaus(0, 1);
        pNew = p + 1 / (fact_cd)*pR * gRandom->Gaus(0, 1) * p;
        std::cout << "mom " << p << "prot fact_cd : " << 1 / fact_cd << std::endl;

      } else if (status_part <= 4000) {
        fact_fd = (0.000055 * pow(p, 3) + (-0.004178) * pow(p, 2) + (0.08473) * p + 0.170);

        phiNew = phi + 1 / (fact_fd)*phiR * gRandom->Gaus(0, 1);
        thetaNew = theta + 1 / (fact_fd)*thetaR * gRandom->Gaus(0, 1);
        pNew = p + 1 / (fact_fd)*pR * gRandom->Gaus(0, 1) * p;

        std::cout << "mom " << p << "prot fact_fd : " << 1 / fact_fd << std::endl;
      }
    }

    else if (part_id == PIP) {
      double fact_cd = 0;
      double fact_fd = 0;
      if (status_part > 4000) {
        fact_cd = (-0.000414) * pow(p, 2) + (0.038539) * p + 0.355856;
        std::cout << "mom " << p << "pip fact_cd : " << 1 / fact_cd << std::endl;

        phiNew = phi + 1 / (fact_cd)*phiR * gRandom->Gaus(0, 1);
        thetaNew = theta + 1 / (fact_cd)*thetaR * gRandom->Gaus(0, 1);
        pNew = p + 1 / (fact_cd)*pR * gRandom->Gaus(0, 1) * p;
      } else if (status_part <= 4000) {
        fact_fd = 0.000075 * pow(p, 3) + (-0.000582) * pow(p, 2) + (-0.015386) * p + 0.772866;
        std::cout << "mom " << p << "pip fact_fd : " << 1 / fact_fd << std::endl;

        phiNew = phi + 1 / (fact_fd)*phiR * gRandom->Gaus(0, 1);
        thetaNew = theta + 1 / (fact_fd)*thetaR * gRandom->Gaus(0, 1);
        pNew = p + 1 / (fact_fd)*pR * gRandom->Gaus(0, 1) * p;
      }
    }

    else if (part_id == PIM) {
      double fact_cd = 0;
      double fact_fd = 0;
      if (status_part > 4000) {
        fact_cd = 0.002389 * pow(p, 2) + (0.007391) * p + 0.322767;
        std::cout << "mom " << p << "pim fact_cd : " << 1 / fact_cd << std::endl;

        phiNew = phi + 1 / (fact_cd)*phiR * gRandom->Gaus(0, 1);
        thetaNew = theta + 1 / (fact_cd)*thetaR * gRandom->Gaus(0, 1);
        pNew = p + 1 / (fact_cd)*pR * gRandom->Gaus(0, 1) * p;
      } else if (status_part <= 4000) {
        fact_fd = 0.000266 * pow(p, 3) + (-0.007222) * pow(p, 2) + (0.058862) * p + 0.544308;
        std::cout << "mom " << p << "pim fact_cd : " << 1 / fact_fd << std::endl;

        phiNew = phi + 1 / (fact_fd)*phiR * gRandom->Gaus(0, 1);
        thetaNew = theta + 1 / (fact_fd)*thetaR * gRandom->Gaus(0, 1);
        pNew = p + 1 / (fact_fd)*pR * gRandom->Gaus(0, 1) * p;
      }
    }
  }
  // /// smearing valerii's function
  // void SmearingFunc1(double p_gen, double p_rec, double &pNew, double frac) {
  //   // Generate new rec momentum
  //   pNew = p_rec + frac * (p_rec - p_gen);
  // }

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
  inline float Phi_Elec() { return _phi_elec; }

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
  float MM_mPim();
  float MM2_mPim();
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

  inline float thetaDCr1Prot() { return _thetaDC_r1_Prot; }
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

    /////////_channelTwoPi_excl &= (_hasE && (_hasP || _hasPip || _hasPim )); // just a trick to see everything
    /// reconstructed dont use it in normal condition

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
  std::unique_ptr<TLorentzVector> _beam_mc;
  std::unique_ptr<TLorentzVector> _elec_mc;
  std::unique_ptr<TLorentzVector> _gamma_mc;
  std::unique_ptr<TLorentzVector> _prot_mc;
  std::unique_ptr<TLorentzVector> _pip_mc;
  std::unique_ptr<TLorentzVector> _pim_mc;
  std::unique_ptr<TLorentzVector> _other_mc;

  std::unique_ptr<TLorentzVector> _boosted_gamma_mc;
  std::unique_ptr<TLorentzVector> _boosted_prot_mc;
  std::unique_ptr<TLorentzVector> _boosted_pip_mc;
  std::unique_ptr<TLorentzVector> _boosted_pim_mc;
  std::unique_ptr<TLorentzVector> _boosted_pim_measured_mc;

  bool _is_boosted_mc = false;
  float _MCinv_Ppip = NAN;
  float _MCinv_Ppim = NAN;
  float _MCinv_pip_pim = NAN;

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

  TVector3 _prot_Vect3_mc;
  TVector3 _pip_Vect3_mc;
  TVector3 _pim_Vect3_mc_measured;

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

  void boost_mc();
  void MCinvMassPpip();
  void MCinvMassPpim();
  void MCinvMasspippim();
  float MCinv_Ppip();
  float MCinv_Ppim();
  float MCinv_pip_pim();

  float scalar_triple_product_thrown();

  std::string CsvHeader();
  std::string ReacToCsv();
};

#endif
