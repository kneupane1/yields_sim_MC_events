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

  float _W = NAN;
  float _Q2 = NAN;

  float _theta_e = NAN;
  float _theta_star = NAN;
  float _phi_star = NAN;

  void SetElec();

 public:
  Reaction(){};
  Reaction(const std::shared_ptr<Branches12> &data, float beam_energy);
  ~Reaction();
  inline float weight() {
    return _data->mc_weight();
    // return 1.0;
  }

  inline bool mc() { return _mc; }
  void SetProton(int i);
  void SetPip(int i);
  void SetPim(int i);
  void SetOther(int i);
  void SetNeutron(int i);
  // missingPim
  float pim_theta_lab();
  float pim_Phi_lab();
  float pim_momentum();
  float pim_theta_lab_measured();
  float pim_Phi_lab_measured();
  float pim_momentum_measured();

  float pim_theta_cm();
  float pim_Phi_cm();
  float pim_momentum_cm();
  float pim_theta_cm_measured();
  float pim_Phi_cm_measured();
  float pim_momentum_cm_measured();

  void boost();

  inline float Theta_star() { return _theta_star; }
  inline float Phi_star() { return _phi_star; }
  inline float Theta_E() { return _theta_e; }

  void CalcMissMass();
  float MM();
  float MM2();
  float MM2_exclusive();
  virtual std::string CsvHeader();
  virtual std::string ReacToCsv();

  inline float W() { return _W; }
  inline float Q2() { return _Q2; }
  float_t scalar_triple_product();

  inline short sec() { return _data->dc_sec(0); }
  inline int det() { return abs(_data->status(0) / 1000); }

  inline bool TwoPion_missingPim() {
    bool _channelTwoPi = true;
    _channelTwoPi &= ((_numProt == 1 && _numPip == 1 /*&&  _numPim == 1*/) &&
                          (_hasE && _hasP && _hasPip
                           /*&&!_hasPim && !_hasNeutron
                              &&!_hasOther*/));
    return _channelTwoPi;
  }
  inline bool TwoPion_exclusive() {
    bool _channelTwoPi_excl = true;

    return ((_numProt == 1 && _numPip == 1 && _numPim == 1) &&
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

  // float _alpha_ppim_pipip_thrown_mc = NAN;

 public:
  void SetMCProton(int i);
  void SetMCPip(int i);
  void SetMCPim(int i);
  void SetMCOther(int i);

  MCReaction(const std::shared_ptr<Branches12> &data, float beam_energy);
  void SetMCElec();
  inline float weight() { return _data->mc_weight(); }
  inline float W() { return _W_mc; }
  inline float Q2() { return _Q2_mc; }
  std::string CsvHeader();
  std::string ReacToCsv();
};

#endif
