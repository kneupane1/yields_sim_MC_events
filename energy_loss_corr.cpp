

#include <iostream>
#include "TMath.h"
#include "TLorentzRotation.h"
#include "TLorentzVector.h"

std::unique_ptr<TLorentzVector> _Energy_loss_uncorr_pim;
std::unique_ptr<TLorentzVector> _Energy_loss_uncorr_pip;
std::unique_ptr<TLorentzVector> _Energy_loss_uncorr_prot;

std::unique_ptr<TLorentzVector> _prot;
std::unique_ptr<TLorentzVector> _pip;
std::unique_ptr<TLorentzVector> _pim;

_Energy_loss_uncorr_prot = std::make_unique<TLorentzVector>();
_Energy_loss_uncorr_pip = std::make_unique<TLorentzVector>();
_Energy_loss_uncorr_pim = std::make_unique<TLorentzVector>();

_prot = std::make_unique<TLorentzVector>();
_pip = std::make_unique<TLorentzVector>();
_pim = std::make_unique<TLorentzVector>();

// PDG particle masses in GeV/c2
static const float MASS_P = 0.93827203;
static const float MASS_PIP = 0.13957018;
static const float MASS_PIM = 0.13957018;

static const float PI = TMath::Pi();

void SetProton(int i);
void SetPip(int i);
void SetPim(int i);

double _prot_mom_prime = NAN;
double _prot_mom = NAN;
double _prot_theta = NAN;
double _prot_mom_tmt = NAN;
double _prot_mom_uncorr = NAN;
float _E_corr_val_prot = NAN;

double _px_prime_prot_E = NAN;
double _py_prime_prot_E = NAN;
double _pz_prime_prot_E = NAN;

double _pip_mom_prime = NAN;
double _pip_mom = NAN;
double _pip_theta = NAN;
double _pip_mom_tmt = NAN;
double _pip_mom_uncorr = NAN;
float _E_corr_val_pip = NAN;

double _px_prime_pip_E = NAN;
double _py_prime_pip_E = NAN;
double _pz_prime_pip_E = NAN;

double _pim_mom_prime = NAN;
double _pim_mom = NAN;
double _pim_theta = NAN;
double _pim_mom_tmt = NAN;
double _pim_mom_uncorr = NAN;
float _E_corr_val_pim = NAN;

double _px_prime_pim_E = NAN;
double _py_prime_pim_E = NAN;
double _pz_prime_pim_E = NAN;

////////////// For Prot energy loss corrections

void SetProton(int i) {
  _Energy_loss_uncorr_prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
  _prot_mom_uncorr = _Energy_loss_uncorr_prot->P();
  _prot_theta = _Energy_loss_uncorr_prot->Theta() * 180 / PI;

  if (abs(_data->status(i)) < 4000) {  // FD Corrections
    if (_prot_theta <= 27) {
      _E_corr_val_prot = -0.00078846 * pow(_prot_mom_uncorr, 5) + 0.0093734 * pow(_prot_mom_uncorr, 4) -
                         0.04277868 * pow(_prot_mom_uncorr, 3) + 0.09421284 * pow(_prot_mom_uncorr, 2) -
                         0.10095842 * (_prot_mom_uncorr) + 0.04567203;
    } else {
      _E_corr_val_prot = -0.0023389 * pow(_prot_mom_uncorr, 5) + 0.02838603 * pow(_prot_mom_uncorr, 4) -
                         0.13214962 * pow(_prot_mom_uncorr, 3) + 0.29609571 * pow(_prot_mom_uncorr, 2) -
                         0.32307424 * (_prot_mom_uncorr) + 0.14742569;
    }
  } else if (abs(_data->status(i)) >= 4000) {  // CD Corrections
    _E_corr_val_prot = 0.0;
    _E_corr_val_prot = ((-9.30990933e-05) * pow(_prot_theta, 3) + (1.23584235e-02) * pow(_prot_theta, 2) +
                        (-5.42538215e-01) * (_prot_theta) + 7.87921215e+00) *
                           pow(_prot_mom_uncorr, 3) +

                       (4.17955911e-04 * pow(_prot_theta, 3) + (-5.53676478e-02) * pow(_prot_theta, 2) +
                        (2.42642631e+00) * (_prot_theta) + (-3.51829220e+01)) *
                           pow(_prot_mom_uncorr, 2) +

                       ((-5.58084320e-04) * pow(_prot_theta, 3) + (7.38670367e-02) * pow(_prot_theta, 2) +
                        (-3.23723227e+00) * (_prot_theta) + 4.69456718e+01) *
                           (_prot_mom_uncorr) +

                       ((2.40014720e-04) * pow(_prot_theta, 3) + (-3.17071405e-02) * pow(_prot_theta, 2) +
                        (1.38769727e+00 * (_prot_theta)) + (-2.01072704e+01));
  }

  _prot_mom_tmt = _prot_mom_uncorr + _E_corr_val_prot;

  _px_prime_prot_E = _data->px(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
  _py_prime_prot_E = _data->py(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
  _pz_prime_prot_E = _data->pz(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));

  _prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);  // energy loss corrected prot
}

////////////// For Pip energy loss corrections

void SetPip(int i) {
  _Energy_loss_uncorr_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
  _pip_mom_uncorr = _Energy_loss_uncorr_pip->P();
  _pip_theta = _Energy_loss_uncorr_pip->Theta() * 180 / PI;

  if (abs(_data->status(i)) < 4000) {  // FD Corrections

    if (_pip_theta <= 27) {
      _E_corr_val_pip = 9.21970527e-05 * pow(_pip_mom_uncorr, 3) - 3.70500143e-04 * pow(_pip_mom_uncorr, 2) +
                        2.78880101e-04 * (_pip_mom_uncorr) + 2.66040566e-03;

    } else {
      _E_corr_val_pip = -0.00010482 * pow(_pip_mom_uncorr, 3) + 0.00080463 * pow(_pip_mom_uncorr, 2) -
                        0.0022871 * (_pip_mom_uncorr) + 0.00831496;
    }
  } else if (abs(_data->status(i)) >= 4000) {  // CD Corrections
    _E_corr_val_pip = 0.0;

    _E_corr_val_pip = (-6.50509539e-07 * pow(_pip_theta, 3) + 1.31547371e-04 * pow(_pip_theta, 2) +
                       (-7.99024673e-03) * (_pip_theta) + 1.60563630e-01) *
                          pow(_pip_mom_uncorr, 3) +

                      (2.48202211e-06 * pow(_pip_theta, 3) + (-5.15757241e-04) * pow(_pip_theta, 2) +
                       3.19833135e-02 * (_pip_theta) + (-6.53476057e-01)) *
                          pow(_pip_mom_uncorr, 2) +

                      (-2.71923009e-06 * pow(_pip_theta, 3) + 5.80375203e-04 * pow(_pip_theta, 2) +
                       (-3.75941898e-02) * (_pip_theta) + 7.80443724e-01) *
                          (_pip_mom_uncorr) +

                      4.62456800e-07 * pow(_pip_theta, 3) + (-1.08401698e-04) * pow(_pip_theta, 2) +
                      8.09261138e-03 * (_pip_theta)-2.05315604e-01;
  }
  _pip_mom_tmt = _pip_mom_uncorr + _E_corr_val_pip;

  _px_prime_pip_E = _data->px(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
  _py_prime_pip_E = _data->py(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
  _pz_prime_pip_E = _data->pz(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));

  _pip->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);  // energy loss corrected pip
}

////////////// For Pim energy loss corrections
void SetPim(int i) {
  _Energy_loss_uncorr_pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
  _pim_mom_uncorr = _Energy_loss_uncorr_pim->P();
  _pim_theta = _Energy_loss_uncorr_pim->Theta() * 180 / PI;

  if (abs(_data->status(i)) < 4000) {  // FD Corrections
    _sectorPim = _data->dc_sec(i);

    if (_pim_theta <= 27) {
      _E_corr_val_pim = -0.00035275 * pow(_pim_mom_uncorr, 3) + 0.00291237 * pow(_pim_mom_uncorr, 2) -
                        0.00681058 * (_pim_mom_uncorr) + 0.00736721;

    } else {
      _E_corr_val_pim = 0.00019358 * pow(_pim_mom_uncorr, 3) - 0.00103456 * pow(_pim_mom_uncorr, 2) +
                        0.00024772 * (_pim_mom_uncorr) + 0.00735159;
    }
  } else if (abs(_data->status(i)) >= 4000) {  // CD Corrections
    _E_corr_val_pim = 0.0;

    _E_corr_val_pim = (-4.94426765e-07 * pow(_pim_theta, 3) + 9.85729368e-05 * pow(_pim_theta, 2) +
                       (-5.85778699e-03) * (_pim_theta) + 1.17447168e-01) *
                          pow(_pim_mom_uncorr, 3) +

                      (1.75953956e-06 * pow(_pim_theta, 3) + (-3.63382515e-04) * pow(_pim_theta, 2) +
                       2.21447425e-02 * (_pim_theta) + (-4.54844509e-01)) *
                          pow(_pim_mom_uncorr, 2) +

                      (-1.90446515e-06 * pow(_pim_theta, 3) + 4.08768480e-04 * pow(_pim_theta, 2) +
                       (-2.65277055e-02) * (_pim_theta) + 5.57286393e-01) *
                          (_pim_mom_uncorr) +

                      2.05653097e-07 * pow(_pim_theta, 3) + (-5.44018546e-05) * pow(_pim_theta, 2) +
                      4.61561853e-03 * (_pim_theta)-1.35303212e-01;
  }

  // _pim_mom = _pim_mom_uncorr + _E_corr_val_pim; // first iteration

  _pim_mom_tmt = _pim_mom_uncorr + _E_corr_val_pim;  // first iteration

  _px_prime_pim_E = _data->px(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
  _py_prime_pim_E = _data->py(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
  _pz_prime_pim_E = _data->pz(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));

  _pim->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, MASS_PIM);  // energy loss corrected pim
}