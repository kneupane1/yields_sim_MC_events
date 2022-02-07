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

 public:
  Reaction(){};
  Reaction(const std::shared_ptr<Branches12> &data, float beam_energy);
  ~Reaction();
  inline float weight() {
    // return _data->mc_weight();
    return 1.0;
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

  inline float W_after() { return _W_after; }

  float_t scalar_triple_product();

  inline short sec() { return _data->dc_sec(0); }
  inline int det() { return abs(_data->status(0) / 1000); }

  inline bool TwoPion_missingPim() {
    bool _channelTwoPi = true;
    _channelTwoPi &= ((_numProt == 1 ) && (_hasE && _hasP));
    return _channelTwoPi;
  }

  // inline bool TwoPion_missingPim() {
  //   bool _channelTwoPi = true;
  //   _channelTwoPi &= ((_numProt == 1 && _numPip == 1 ) && (_hasE && _hasP && _hasPip));
  //   return _channelTwoPi;
  // }

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
