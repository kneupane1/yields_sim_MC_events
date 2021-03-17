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

double _beam_energy = 7.5;
std::unique_ptr<TLorentzVector> _beam;
std::unique_ptr<TLorentzVector> _elec;
std::unique_ptr<TLorentzVector> _gamma;
std::unique_ptr<TLorentzVector> _target;
std::unique_ptr<TLorentzVector> _prot;
std::unique_ptr<TLorentzVector> _pip;
std::unique_ptr<TLorentzVector> _pim;
std::unique_ptr<TLorentzVector> _other;
std::unique_ptr<TLorentzVector> _neutron;

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

float _W = NAN;
float _Q2 = NAN;

float _theta_e = NAN;
float _theta_star = NAN;
float _phi_star = NAN;

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
void boost();
inline float Theta_star() {
        return _theta_star;
}
inline float Phi_star() {
        return _phi_star;
}
inline float Theta_E() {
        return _theta_e;
}

void CalcMissMass();
float MM();
float MM2();
virtual std::string CsvHeader();
virtual std::string ReacToCsv();

inline float W() {
        return _W;
}
inline float Q2() {
        return _Q2;
}
inline short sec() {
        return _data->dc_sec(0);
}
inline int det() {
        return abs(_data->status(0) / 1000);
}

inline bool TwoPion() {
        return ((_numPip == 1 && _numPim == 1) && (_hasE && !_hasP && _hasPip && _hasPim && !_hasNeutron && !_hasOther));
}
inline bool ProtonPim() {
        return ((_numProt == 1 && _numPim == 1) && (_hasE && _hasP && !_hasPip && _hasPim && !_hasNeutron && !_hasOther));
}
inline bool SinglePip() {
        return ((_numPip == 1) && (_hasE && !_hasP && _hasPip && !_hasPim));
}
inline bool SingleP() {
        return ((_numProt == 1) && (_hasE && _hasP && !_hasPip && !_hasPim && !_hasNeutron && !_hasOther));
}

inline bool NeutronPip() {
        bool _channel = true;
        _channel &= (Reaction::SinglePip() && Reaction::MM() >= 0.85 && Reaction::MM() <= 1.0);
        return _channel;
}

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

std::unique_ptr<TLorentzVector> _elec_mc;
std::unique_ptr<TLorentzVector> _gamma_mc;
std::unique_ptr<TLorentzVector> _prot_mc;
std::unique_ptr<TLorentzVector> _pip_mc;
std::unique_ptr<TLorentzVector> _pim_mc;
std::unique_ptr<TLorentzVector> _other_mc;

float _alpha_ppim_pipip_thrown_mc = NAN;

public:
void SetMCProton(int i);
void SetMCPip(int i);
void SetMCPim(int i);
void SetMCOther(int i);

MCReaction(const std::shared_ptr<Branches12> &data, float beam_energy);
void SetMCElec();
inline float weight() {
        return _data->mc_weight();
}
inline float W() {
        return _W_mc;
}
inline float Q2() {
        return _Q2_mc;
}
std::string CsvHeader();
std::string ReacToCsv();
};

#endif
