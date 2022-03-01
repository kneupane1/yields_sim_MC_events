/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include <iostream>
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "branches.hpp"
#include "constants.hpp"
#include "physics.hpp"

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

  float _W_after = NAN;

  float _theta_e = NAN;
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

  void SetElec();

  // momentum corrections

  double xx[54] = {
      0.0263375, 0.0158871,  0.0130852,  -0.00366006, 0.00694866,  0.0197195, 0.00767067, 0.00480921,  -0.0175756,
      0.0252757, 0.0156601,  0.00984872, 0.00244435,  0.00681414,  0.0294068, 0.0059881,  0.00286992,  0.0179319,
      0.0171495, 0.00359637, -0.0046115, 0.00314739,  0.0136338,   0.0768753, 0.00675454, -0.0118234,  -0.0288654,
      0.0189465, 0.0131816,  0.0262004,  0.00375165,  0.00907457,  0.0486894, 0.00806305, 0.0006999,   0.00527513,
      0.0116485, 0.0105681,  0.0149848,  0.000318094, -0.00480124, 0.0395545, 0.00824216, -0.00070659, -0.0057075,
      0.0213057, 0.0112999,  0.0100216,  0.000653685, 0.0093174,   0.0822385, 0.00808384, 0.000898799, -0.0172692,
  };
  double pars[6][3][3];
  int ipar = 0;

  double _elec_mom_corrected = NAN;
  double _elec_mom = NAN;

  double _cx = NAN;
  double _cy = NAN;
  double _cz = NAN;

  double _px_prime_elec = NAN;
  double _py_prime_elec = NAN;
  double _pz_prime_elec = NAN;

// for pim mom corrections 

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

  static const int Pim_mom_bins = 10;
  float alpha_pim_mom_corr = 0.5;
  double _pim_mom = NAN;
  double _pim_mom_prime = NAN;
  double _pim_mom_uncorr = NAN;
  float _E_corr_val_pim = NAN;

  float min_pim_mom_values[Pim_mom_bins] = {0, 0.36, 0.5, 0.65, 0.8, 0.95, 1.15, 1.45, 1.9, 2.5};
  float max_pim_mom_values[Pim_mom_bins] = {0.36, 0.5, 0.65, 0.8, 0.95, 1.15, 1.45, 1.9, 2.5, 10};
  double pim_mom_corr[Pim_mom_bins] = {-0.004, -0.004, -0.004, -0.012, -0.025, -0.035, -0.035, -0.03, -0.03, -0.03};
  double pim_mom_corr_sim[Pim_mom_bins] = { -0.002, -0.006, -0.002, -0.002, -0.0025, -0.0075, -0.0125, -0.015, -0.015, -0.015};

  static const int Pim_theta_bins = 10;
  float alpha_pim_theta_corr = 0.5;
  double _pim_theta = NAN;
  double _pim_theta_prime = NAN;

  float min_pim_theta_values[Pim_theta_bins] = {0, 19.5, 24, 28, 33, 38, 45, 52, 60, 75};
  float max_pim_theta_values[Pim_theta_bins] = {19.5, 24, 28, 33, 38, 45, 52, 60, 75, 180};
  double pim_theta_corr[Pim_theta_bins] = {0.25, 0.25, 0.25, 0.3, -0.1, 1.05, 0.75, 0.2, 0.3, -1.5};
  double pim_theta_corr_sim[Pim_theta_bins] = {0.075, 0.125, 0.175, 0.25, 0.25, 0.375, 0.225, 0.3, 0.45, 0.75,};

  static const int Pim_phi_bins = 11;
  float alpha_pim_phi_corr = 0.5;
  double _pim_phi = NAN;
  double _pim_phi_prime = NAN;
  float min_pim_phi_values[Pim_phi_bins] = {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300};
  float max_pim_phi_values[Pim_phi_bins] = {30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 360};
  double pim_phi_corr[Pim_phi_bins] = {0.5, 0.1, 0.1, 0.3, -0.1, -0.3, -0.3, -0.1, -0.1, -0.1, 0.3};
  double pim_phi_corr_sim[Pim_phi_bins] = { -0.1, -0.1, -0.1, -0.1, -0.1, 0.1, -0.1, 0.1, -0.1, 0.1, -0.1};

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

  static const int Pip_mom_bins = 10;
  float alpha_pip_mom_corr = 0.5;
  double _pip_mom = NAN;
  double _pip_mom_prime = NAN;
  double _pip_mom_uncorr = NAN;
  float _E_corr_val_pip = NAN;

  float min_pip_mom_values[Pip_mom_bins] = {0, 0.36, 0.5, 0.65, 0.8, 0.95, 1.15, 1.45, 1.9, 2.5};
  float max_pip_mom_values[Pip_mom_bins] = {0.36, 0.5, 0.65, 0.8, 0.95, 1.15, 1.45, 1.9, 2.5, 10};
  double pip_mom_corr[Pip_mom_bins] = { -0.003, -0.009, -0.021, -0.021, -0.036, -0.028, -0.028, -0.025, -0.025, -0.025,};
  double pip_mom_corr_sim[Pip_mom_bins] = {
      -0.003, -0.003, -0.003, -0.009, -0.012, -0.012, -0.012, -0.015, -0.015, -0.015

  };

  static const int Pip_theta_bins = 10;
  float alpha_pip_theta_corr = 0.5;
  double _pip_theta = NAN;
  double _pip_theta_prime = NAN;

  float min_pip_theta_values[Pip_theta_bins] = {0, 19.5, 24, 28, 33, 38, 45, 52, 60, 75};
  float max_pip_theta_values[Pip_theta_bins] = {19.5, 24, 28, 33, 38, 45, 52, 60, 75, 180};
  double pip_theta_corr[Pip_theta_bins] = { 0.225, 0.225, 0.075, 0.15, -0.15, 1.4, 0.2, -0.75, -2, -6};
  double pip_theta_corr_sim[Pip_theta_bins] = {
      0.1875, 0.1875, 0.1875, 0.375, 0.225, 0.3, 0.3, 0.375, 0.6, 0.2
  };

  static const int Pip_phi_bins = 11;
  float alpha_pip_phi_corr = 0.5;
  double _pip_phi = NAN;
  double _pip_phi_prime = NAN;
  float min_pip_phi_values[Pip_phi_bins] = {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300};
  float max_pip_phi_values[Pip_phi_bins] = {30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 360};
  double pip_phi_corr[Pip_phi_bins] = {0.9, 0.7, 0.5, 0.3, 0.1, -0.1, -0.1, 0.1, 0.3, 0.3, 0.7};
  double pip_phi_corr_sim[Pip_phi_bins] = {
      0.1, -0.1, 0.1, -0.1, 0.1, -0.1, 0.1, 0.1, 0.1, -0.1, -0.1
  };

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


  static const int Prot_mom_bins = 8;
  float alpha_prot_mom_corr = 0.5;
  double _prot_mom_prime = NAN;
  double _prot_mom = NAN;
  double _prot_mom_uncorr = NAN;
  float _E_corr_val_prot = NAN;

      float min_prot_mom_values[Prot_mom_bins] = {0, 0.65, 0.8, 0.95, 1.15, 1.45, 1.9, 2.5};
  float max_prot_mom_values[Prot_mom_bins] = { 0.65, 0.8, 0.95, 1.15, 1.45, 1.9, 2.5,10};
  double prot_mom_corr[Prot_mom_bins] = {-0.021, -0.021, -0.025, -0.025, -0.025, -0.03, -0.03, -0.042};
  double prot_mom_corr_sim[Prot_mom_bins] = {
      -0.0135, -0.0135, -0.0125, -0.0125, -0.0125, -0.009, -0.009, -0.009,
  };

  static const int Prot_theta_bins = 10;
  float alpha_prot_theta_corr = 0.5;
  double _prot_theta = NAN;
  double _prot_theta_prime = NAN;

  float min_prot_theta_values[Prot_theta_bins] = {0, 10, 15, 19.5, 24, 28, 33, 38};
  float max_prot_theta_values[Prot_theta_bins] = {10, 15, 19.5, 24, 28, 33, 38, 180};
  double prot_theta_corr[Prot_theta_bins] = { 0.375, 0.225, 0.225, 0.225, 0.225, 0.15, -0.45, -0.15};
  double prot_theta_corr_sim[Prot_theta_bins] = {
      0.1125, 0.1125, 0.1125, 0.1125, 0.1125, 0.225, 0.075, 0.075,
  };

  static const int Prot_phi_bins = 11;
  float alpha_prot_phi_corr = 0.5;
  double _prot_phi = NAN;
  double _prot_phi_prime = NAN;
  float min_prot_phi_values[Prot_phi_bins] = {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300};
  float max_prot_phi_values[Prot_phi_bins] = {30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 360};
  double prot_phi_corr[Prot_phi_bins] = {0.9, 0.7, 0.5, 0.5, 0.1, -0.1, -0.3, 0.1, 0.1, 0.3, 0.7};
  double prot_phi_corr_sim[Prot_phi_bins] = {
      -0.1, -0.1, -0.1, 0.1, 0.1, -0.1, 0.1, 0.1, 0.1, 0.1, -0.1,
  };

 public:
  Reaction(){};
  Reaction(const std::shared_ptr<Branches12> &data, float beam_energy);
  ~Reaction();
  inline float weight() {
    return _data->mc_weight();
    // return 1.0;
  }

  // momentum correction
  double dpp(float px, float py, float pz, int sec, int ivec);
  double Corr_elec_mom();
  double elec_mom();

  inline bool mc() { return _mc; }
  void SetProton(int i);
  void SetPip(int i);
  void SetPim(int i);
  void SetOther(int i);
  void SetNeutron(int i);

  bool ctof_prot();
  bool ftof_prot();
  bool ctof_pip();
  bool ftof_pip();
  bool ctof_pim();
  bool ftof_pim();

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
  float elec_E();

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
  inline float Theta_E() { return _theta_e; }

  void CalcMissMass();
  float MM();
  float MM2();
  float MM2_exclusive();
  float Energy_excl();
  float MM2_mPip();
  float MM2_mProt();

  virtual std::string CsvHeader();
  virtual std::string ReacToCsv();

  inline float W() { return _W; }
  inline float Q2() { return _Q2; }

  // inline float W_after() { return _W_after; }

  float_t scalar_triple_product();

  inline short sec() { return _data->dc_sec(0); }
  inline int det() { return abs(_data->status(0) / 1000); }

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

 public:
  void SetMCProton(int i);
  void SetMCPip(int i);
  void SetMCPim(int i);
  void SetMCOther(int i);

  MCReaction(const std::shared_ptr<Branches12> &data, float beam_energy);
  void SetMCElec();
  inline float weight() { return _data->mc_weight(); }
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
