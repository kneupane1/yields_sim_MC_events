/************************************************************************/
/*  Created by Nick Tyler*/
/*	University Of South Carolina*/
/************************************************************************/

#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include "TLorentzVector.h"
#include "branches.hpp"
#include "constants.hpp"
#include "physics.hpp"
#include <iostream>

class Reaction {
protected:
std::shared_ptr<Branches12> _data;

double _beam_energy = 10.6;
std::unique_ptr<TLorentzVector> _beam;
std::unique_ptr<TLorentzVector> _elec;
std::unique_ptr<TLorentzVector> _gamma;
std::unique_ptr<TLorentzVector> _target;

// std::unique_ptr<TLorentzVector> _x_mu;
// std::unique_ptr<TLorentzVector> _photons;
std::vector<std::shared_ptr<TLorentzVector> > _photons;

std::unique_ptr<TLorentzVector> _prot;
std::unique_ptr<TLorentzVector> _pip;
std::unique_ptr<TLorentzVector> _pim;
std::unique_ptr<TLorentzVector> _boosted_gamma;
std::unique_ptr<TLorentzVector> _boosted_prot;
std::unique_ptr<TLorentzVector> _boosted_pip;
std::unique_ptr<TLorentzVector> _boosted_pim;
std::unique_ptr<TLorentzVector> _other;
std::unique_ptr<TLorentzVector> _neutron;

std::unique_ptr<TLorentzVector> _missingPim;
std::unique_ptr<TLorentzVector> _boosted_pim_measured;


float _weight = NAN;

bool _mc = false;

bool _is_boosted = false;

bool _hasE = false;
bool _hasP = false;

bool _hasPip = false;
bool _hasPim = false;
bool _hasOther = false;
bool _hasNeutron = false;

short _numPart = 0;
short _numProt = 0;
short _numPip = 0;
short _numPim = 0;
short _numPos = 0;
short _numNeg = 0;
short _numNeutral = 0;
short _numOther = 0;
short _numPhoton = 0;
short _sector = -1;

float _MM = NAN;
float _MM2 = NAN;

float _MM_mpip = NAN;
float _MM2_mpip = NAN;
float _MM_mprot = NAN;
float _MM2_mprot = NAN;
float _pi0_mass = NAN;
float _MM2_exclusive = NAN;
float _W = NAN;
float _Q2 = NAN;

float _W_e_Prot = NAN;
float _Q2_e_Prot = NAN;

float _inv_Ppip = NAN;
float _inv_Ppim = NAN;
float _inv_pip_pim = NAN;
float _W_P2pi = NAN;

float _phi_gamma = NAN;
float _phi_prot = NAN;
float _phi_pip = NAN;
float _phi_pim = NAN;

float _alpha_ppip_pipim = NAN;
float _alpha_pippim_pipf = NAN;
float _alpha_ppim_pipip = NAN;

float _x_mu_Px = NAN;
float _x_mu_Py = NAN;
float _x_mu_Pz = NAN;

float _x_mu_m = NAN;
float _x_mu_m2 = NAN;
float _x_mu_P = NAN;
float _x_mu_E = NAN;
float _x_mu_theta = NAN;
float _beam_theta = NAN;
float _elec_theta = NAN;
float _E_elec = NAN;
float _pim_theta_measured = NAN;
short _pim_sec = -9999;
float _prot_status = NAN;
float _pip_status = NAN;
float _pim_status = NAN;


void SetElec();

public:
Reaction(){
};
Reaction(const std::shared_ptr<Branches12> &data, float beam_energy);
~Reaction();

inline bool mc() {
        return _mc;
}
void SetProton(int i);
void SetPip(int i);
void SetPim(int i);
void SetOther(int i);
void SetNeutron(int i);
short pim_sec();
bool ctof_prot();
bool ftof_prot();
bool ctof_pip();
bool ftof_pip();
bool ctof_pim();
bool ftof_pim();


// float q_3_();
// TLorentzVector p_mu_prime_cm();
// TLorentzVector pip_mu_prime_cm();
// TLorentzVector pim_mu_prime_cm();
// float theta_();
virtual std::string CsvHeader();
virtual std::string ReacToCsv();

void boost();
void CalcMissMass();
float MM();
float MM2();
float MM2_mpip();
float MM2_mprot();
void CalcMassPi0();
float pi0_mass();
float MM2_exclusive();
inline float weight() {

        return _data->mc_weight();  //
        //return 1.0;
}

float inv_Ppip();
float inv_Ppim();
float inv_pip_pim();
float w_P2pi_rec();

void W_2pi_P();
void invMassPpip();
void invMassPpim();
void invMasspippim();

float prot_theta_lab();
float pip_theta_lab();
float pim_theta_lab();
float pim_theta_lab_measured();

float prot_theta();
float pip_theta();
float pim_theta();

float elec_momentum();
float prot_momentum();
float pip_momentum();
float pim_momentum();
float pim_momentum_measured();
float pim_E();
float pim_E_measured();

float theta_beam();
float theta_elec();
float E_elec();

float theta_x_mu();
float P_x_mu();
float E_x_mu();
float Px_x_mu();
float Py_x_mu();
float Pz_x_mu();
float M_x_mu();
float M2_x_mu();

void AlphaCalc();

float gamma_Phi();
float prot_Phi();
float pip_Phi();
float pim_Phi();
float prot_Phi_lab();
float pip_Phi_lab();
float pim_Phi_lab();
float pim_Phi_lab_measured();

float alpha_ppip_pipim();
float alpha_pippim_pipf();
float alpha_ppim_pipip();


//missingPim
float diff_pim_theta_lab();
float diff_pim_Phi_lab();
float diff_pim_momentum();
//missingPip
float reconstructed_pip_theta_lab();
float reconstructed_pip_mom_lab();
float measured_pip_theta_lab_no_missing();
float measured_pip_mom_lab_no_missing();
float diff_pip_theta_lab();
float diff_pip_Phi_lab();
float diff_pip_momentum();

//missingProt
float reconstructed_prot_theta_lab();
float reconstructed_prot_mom_lab();
float measured_prot_theta_lab_no_missing();
float measured_prot_mom_lab_no_missing();
float diff_prot_theta_lab();
float diff_prot_Phi_lab();
float diff_prot_momentum();


inline float W() {

        return _W;
}
inline float Q2() {
        return _Q2;
}
//
// float theta_x_mu();
// float E_x_mu();
// float P_x_mu();
// float theta_beam();

inline short sec() {
        return _data->dc_sec(0);
}
inline int det() {
        return abs(_data->status(0) / 1000);
}

inline bool MM_cut() {
        return // abs(Reaction::MM2()) < 0.03 // MM2_cut;
               (Reaction::MM2() < 0.08 && Reaction::MM2() > -0.06);
}   // missingPim ko
    // lagi
    // two pi exclusive ko lagi ho

inline bool TwoPion_missingPim() {
        bool _channelTwoPi = true;
        _channelTwoPi &= ((_numProt == 1 && _numPip == 1 /*&&  _numPim == 1*/) &&
                          (_hasE && _hasP && _hasPip
                           /*&&!_hasPim /*&& !_hasNeutron
                              &&!_hasOther*/));
        return _channelTwoPi;
}
inline bool TwoPion_exclusive() {
        return ((_numProt == 1 && _numPip == 1  && _numPim == 1) &&
                (_hasE && _hasP &&
                 _hasPip && _hasPim /*&& !_hasNeutron && !_hasOther*/));
}
inline bool TwoPion_missingPip() {
        bool _channelTwoPi_mpip = true;

        _channelTwoPi_mpip &=((_numProt == 1 && _numPim == 1) &&
                              (_hasE && _hasP &&
                               _hasPim /*&&!_hasPip /*&& !_hasNeutron && !_hasOther*/));
        return _channelTwoPi_mpip;
}
inline bool TwoPion_missingProt() {
        bool _channelTwoPi_mprot = true;
        _channelTwoPi_mprot &= ((_numPip == 1 && _numPim == 1) &&
                                (_hasE && _hasPip && _hasPim /*&&!_hasP /*&& !_hasOther*/));
        return _channelTwoPi_mprot;
}
// inline bool SingleP() {
//   return ((_numProt == 1) && (_hasE && _hasP && !_hasPip && !_hasPim &&
//   !_hasNeutron && !_hasOther));
// }
//
// inline bool NeutronPip() {
//   bool _channel = true;
//   _channel &= ((_numPip == 1) && (_hasE && !_hasP && _hasPip && !_hasPim &&
//   _hasNeutron)) ||
//               (Reaction::SinglePip() && Reaction::MM() >= 0.85 &&
//               Reaction::MM() <= 1.0);
//   return _channel;
// }

const TLorentzVector &e_mu() {
        return *_beam;
}
const TLorentzVector &e_mu_prime() {
        return *_elec;
}
const TLorentzVector &gamma() {
        return *_gamma;
}
};

class MCReaction : public Reaction {
private:
float _weight_mc = NAN;
float _W_mc = NAN;
float _Q2_mc = NAN;

float _MCinv_Ppip = NAN;
float _MCinv_Ppim = NAN;
float _MCinv_pip_pim = NAN;

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

bool _is_boosted_mc = false;

float _MM_mc = NAN;
float _MM2_mc = NAN;
float _MM2_exclusive_mc = NAN;


float _alpha_ppip_pipim_mc = NAN;
float _alpha_pippim_pipf_mc = NAN;
float _alpha_ppim_pipip_mc = NAN;

float _alpha_ppip_pipim_thrown_mc = NAN;
float _alpha_pippim_pipf_thrown_mc = NAN;
float _alpha_ppim_pipip_thrown_mc = NAN;

public:
MCReaction(const std::shared_ptr<Branches12> &data, float beam_energy);
void SetMCElec();
// inline float weight() {
//         return _data->mc_weight();
// }
inline float W_mc() {
        return _W_mc;
}
inline float Q2_mc() {
        return _Q2_mc;
}
void CalcMissMass_mc();
float MM_mc();
float MM2_mc();
float MM2_exclusive_mc();

void SetMCProton(int i);
void SetMCPip(int i);
void SetMCPim(int i);
void SetMCOther(int i);

std::string CsvHeader();
std::string ReacToCsv();

void boost_mc();
void MCinvMassPpip();
void MCinvMassPpim();
void MCinvMasspippim();
float MCinv_Ppip();
float MCinv_Ppim();
float MCinv_pip_pim();
//
// float MCprot_theta();
// float MCpip_theta();
// float MCpim_theta();
//
float MCprot_theta_lab();
float MCpip_theta_lab();
float MCpim_theta_lab();

float prot_momentum_thrown();
float pip_momentum_thrown();
float pim_momentum_thrown();
// float MCgamma_Phi();
// float MCprot_Phi();
// float MCpip_Phi();
// float MCpim_Phi();
//
// float MCalpha_ppip_pipim();
// float MCalpha_pippim_pipf();
// float MCalpha_ppim_pipip();

float MCprot_theta_thrown();
float MCpip_theta_thrown();
float MCpim_theta_thrown();

float MCgamma_Phi_thrown();
float MCprot_Phi_thrown();
float MCpip_Phi_thrown();
float MCpim_Phi_thrown();

void MCAlphaCalc();

float MCalpha_ppip_pipim_thrown();
float MCalpha_pippim_pipf_thrown();
float MCalpha_ppim_pipip_thrown();
};

#endif
