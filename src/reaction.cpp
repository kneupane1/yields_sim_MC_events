/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "reaction.hpp"

Reaction::Reaction(const std::shared_ptr<Branches12>& data, float beam_energy) {
  _data = data;
  _beam = std::make_unique<TLorentzVector>();
  _beam_energy = beam_energy;
  _sector = data->dc_sec(0);

  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  _gamma = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  _elec = std::make_unique<TLorentzVector>();

  _mom_corr_elec = std::make_unique<TLorentzVector>();
  _mom_corr_pim = std::make_unique<TLorentzVector>();
  _mom_corr_pim_th = std::make_unique<TLorentzVector>();
  _mom_corr_pim_ph = std::make_unique<TLorentzVector>();
  _mom_corr_pip = std::make_unique<TLorentzVector>();
  _mom_corr_pip_th = std::make_unique<TLorentzVector>();
  _mom_corr_pip_ph = std::make_unique<TLorentzVector>();
  _mom_corr_prot = std::make_unique<TLorentzVector>();
  _mom_corr_prot_th = std::make_unique<TLorentzVector>();
  _mom_corr_prot_ph = std::make_unique<TLorentzVector>();
  _Energy_loss_uncorr_prot = std::make_unique<TLorentzVector>();
  _Energy_loss_uncorr_pip = std::make_unique<TLorentzVector>();
  _Energy_loss_uncorr_pim = std::make_unique<TLorentzVector>();
  _pim_tmt = std::make_unique<TLorentzVector>();
  _pip_tmt = std::make_unique<TLorentzVector>();

  this->SetElec();

  _prot = std::make_unique<TLorentzVector>();
  _pip = std::make_unique<TLorentzVector>();
  _pim = std::make_unique<TLorentzVector>();
  _other = std::make_unique<TLorentzVector>();
  _neutron = std::make_unique<TLorentzVector>();

  for (int isec = 0; isec < 6; isec++) {
    for (int ivec = 0; ivec < 3; ivec++) {
      double dp1 = xx[ipar++], dp5 = xx[ipar++], dp9 = xx[ipar++];

      pars[isec][ivec][0] = (dp1 - 2 * dp5 + dp9) / 32.;
      pars[isec][ivec][1] = (-7 * dp1) / 16. + (5 * dp5) / 8. - (3 * dp9) / 16.;
      pars[isec][ivec][2] = (45 * dp1) / 32. - (9 * dp5) / 16. + (5 * dp9) / 32.;
    }
  }
}

Reaction::~Reaction() {}

double Reaction::dpp(float px, float py, float pz, int sec, int ivec) {
  double pp = sqrt(px * px + py * py + pz * pz);

  double a = pars[sec - 1][ivec][0], b = pars[sec - 1][ivec][1], c = pars[sec - 1][ivec][2];

  double dp = a * pp * pp + b * pp + c;  // pol2 corr func

  // electron pol1 corr func for each sec and each phi bins
  if (ivec == 0) {
    if (sec == 1) {
      dp = 0.45 * b * (pp - 9) + 0.1 * c;

      // ep 3 phi bins
      // dp = -0.01*b*(pp-9)+1.35*c; //phi<-5
      // dp = 0.6*b*(pp-9)-0.3*c; //-5<phi<5
      // dp = 1.7*b*(pp-9)-1.5*c; //phi>5
    }
    if (sec == 2) {
      dp = -0.15 * b * (pp - 8.0) - 0.3 * c;

      // ep 3 phi bins
      // dp = -0.7*b*(pp-8.0)+0.4*c; //phi<-5
      // dp = -0.05*b*(pp-8.0)-0.4*c; //-5<phi<5
      // dp = 0.01*b*(pp-8.0)-1.5*c; //phi>5
    }
    if (sec == 3) {
      dp = 3. * b * (pp - 5.4) - 0.5 * c;

      // ep 3 phi bins
      // dp = 0.04*b*(pp-5.4)-3.5*c; //phi<-5
      // dp = 0.06*b*(pp-5.4)-3.*c; //-5<phi<5
      // dp = 1.1*b*(pp-5.4)-0.7*c; //phi>5
    }
    if (sec == 4) {
      dp = 0.25 * b * (pp - 9.25) - 0.3 * c;

      // ep 3 phi bins
      // dp = 0.25*b*(pp-9.25)-0.7*c; //phi<-5
      // dp = 0.25*b*(pp-9.25)+0.05*c; //-5<phi<5
      // dp = 0.1*b*(pp-9.25)+1.1*c; //phi>5
    }
    if (sec == 5) {
      dp = 2.2 * b * (pp - 7.5) - 0.5 * c;

      // ep 3 phi bins
      // dp = 2.2*b*(pp-7.5)+0.5*c; //phi<-5
      // dp = 2.2*b*(pp-7.5)-0.1*c; //-5<phi<5
      // dp = 2.2*b*(pp-7.5)-0.6*c; //phi>5
    }
    if (sec == 6) {
      dp = 0.5 * b * (pp - 7) - 0.6 * c;

      // ep 3 phi bins
      // dp = 1.263*b*(pp-7)+0.5*c; //phi<-5
      // dp = 1.*b*(pp-7)-0.5*c; //-5<phi<5
      // dp = 0.5*b*(pp-7)-1.45*c; //phi>5
    }
  }
  return dp / pp;
};

// double fe = dpp(ex, ey, ez, esec, 0) + 1;
// double fpip = dpp(pipx,pipy,pipz,pipsec,1) + 1;
// double fpim = dpp(pimx,pimy,pimz,pimsec,2) + 1;

void Reaction::SetElec() {
  _hasE = true;
  _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);
  *_gamma += *_beam - *_elec;  // be careful you are commenting this only to include the momentum correction

  // Can calculate W and Q2 here (useful for simulations as sim do not have elec mom corrections)
  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);

  // // //One way of  calculating mom - corrected four vectors
  // //   // // _cx = _data->px(0)/_elec->P();
  // //   // // _cy = _data->py(0) / _elec->P();
  // //   // // _cz = _data->pz(0) / _elec->P();
  // //   // _elec_mom_corrected = _elec->P() * (dpp(_data->px(0), _data->py(0), _data->pz(0), _data->dc_sec(0), 0) + 1);

  // //   // _px_prime_elec = _cx * _elec_mom_corrected;
  // //   // _py_prime_elec = _cy * _elec_mom_corrected;
  // //   // _pz_prime_elec = _cz * _elec_mom_corrected; // _mom_corr_elec->SetXYZM(_px_prime_elec, _py_prime_elec,
  // //   // _pz_prime_elec, MASS_E);

  // // //   // mom correction another way
  //   _elec_mom = _elec->P();

  //   _elec_mom_corrected = (dpp(_data->px(0), _data->py(0), _data->pz(0), _data->dc_sec(0), 0) + 1);

  //   _mom_corr_elec->SetPxPyPzE(_data->px(0) * _elec_mom_corrected, _data->py(0) * _elec_mom_corrected,
  //                              _data->pz(0) * _elec_mom_corrected, _elec_mom * _elec_mom_corrected);

  //   *_gamma += *_beam - *_mom_corr_elec;

  //   _W = physics::W_calc(*_beam, *_mom_corr_elec);
  //   _Q2 = physics::Q2_calc(*_beam, *_mom_corr_elec);
}

// double Reaction::Corr_elec_mom() {
//   if (_elec_mom_corrected != _elec_mom_corrected) SetElec();
//   // std::cout << " emec mom corrected " << _elec_mom_corrected << std::endl;

//   return _elec_mom_corrected;
// }

// double Reaction::elec_mom() {
//   if (_elec_mom != _elec_mom) SetElec();
//   // std::cout << " emec mom " << _elec_mom << std::endl;

//   return _elec_mom;
// }

void Reaction::SetProton(int i) {
  _numProt++;
  _numPos++;
  _hasP = true;
  // _prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
  _Energy_loss_uncorr_prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
  // _prot_status = abs(_data->status(i));
  _prot_mom_uncorr = _Energy_loss_uncorr_prot->P();
  _prot_theta = _Energy_loss_uncorr_prot->Theta() * 180 / PI;
  // std::cout << "prot ststus " << _data->status(i) << "   prot theta " << _prot_theta << " prot  mom   "
  //           << _prot_mom_uncorr<< std::endl;
  if (abs(_data->status(i)) < 4000) {
    if (_prot_theta <= 27) {
      _E_corr_val_prot = -0.00078846 * pow(_prot_mom_uncorr, 5) + 0.0093734 * pow(_prot_mom_uncorr, 4) -
                         0.04277868 * pow(_prot_mom_uncorr, 3) + 0.09421284 * pow(_prot_mom_uncorr, 2) -
                         0.10095842 * (_prot_mom_uncorr) + 0.04567203;
    } else {
      _E_corr_val_prot = -0.0023389 * pow(_prot_mom_uncorr, 5) + 0.02838603 * pow(_prot_mom_uncorr, 4) -
                         0.13214962 * pow(_prot_mom_uncorr, 3) + 0.29609571 * pow(_prot_mom_uncorr, 2) -
                         0.32307424 * (_prot_mom_uncorr) + 0.14742569;
    }
  } else if (abs(_data->status(i)) >= 4000) {
    _E_corr_val_prot = 0.01066342 * pow(_prot_mom_uncorr, 2) - 0.05379427 * (_prot_mom_uncorr) + 0.02530928;
  }

  _prot_mom = _prot_mom_uncorr + _E_corr_val_prot;

  _px_prime_prot_E = _data->px(i) * ((_prot_mom) / (_prot_mom_uncorr));
  _py_prime_prot_E = _data->py(i) * ((_prot_mom) / (_prot_mom_uncorr));
  _pz_prime_prot_E = _data->pz(i) * ((_prot_mom) / (_prot_mom_uncorr));
  _prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);

  // if (_prot->Phi() > 0)
  //   _prot_phi = _prot->Phi() * 180 / PI;
  // else if (_prot->Phi() < 0)
  //   _prot_phi = (_prot->Phi() + 2 * PI) * 180 / PI;

  // for (size_t t = 0; t < Prot_theta_bins; t++) {
  //   double theta_min = min_prot_theta_values[t];
  //   double theta_max = max_prot_theta_values[t];
  //   if (_prot_theta > theta_min && _prot_theta < theta_max) {
  //     // for experimental dsta
  //     _prot_theta_prime = _prot_theta - prot_theta_corr[t] * alpha_prot_theta_corr;
  //     // //for simulation data
  //     //       _prot_theta_prime = _prot_theta - prot_theta_corr_sim[t] * alpha_prot_theta_corr;

  //     _px_prime_prot_th = _data->px(i) * (sin(DEG2RAD * _prot_theta_prime) / sin(DEG2RAD * _prot_theta));
  //     _py_prime_prot_th = _data->py(i) * (sin(DEG2RAD * _prot_theta_prime) / sin(DEG2RAD * _prot_theta));
  //     _pz_prime_prot_th = _data->pz(i) * (cos(DEG2RAD * _prot_theta_prime) / cos(DEG2RAD * _prot_theta));
  //     _E_prime_prot_th = sqrt(abs(_px_prime_prot_th * _px_prime_prot_th + _py_prime_prot_th * _py_prime_prot_th +
  //                                 _pz_prime_prot_th * _pz_prime_prot_th));
  //     _mom_corr_prot_th->SetPxPyPzE(_px_prime_prot_th, _py_prime_prot_th, _pz_prime_prot_th, _E_prime_prot_th);
  //   }
  // }

  //   for (size_t p = 0; p < Prot_phi_bins; p++) {
  //     double phi_min = min_prot_phi_values[p];
  //     double phi_max = max_prot_phi_values[p];
  //     if (_prot_phi > phi_min && _prot_phi < phi_max) {
  //       // for experimental data
  //       _prot_phi_prime = _prot_phi - prot_phi_corr[p] * alpha_prot_phi_corr;

  //       // // for simulation data
  //       // _prot_phi_prime = _prot_phi - prot_phi_corr_sim[p] * alpha_prot_phi_corr;

  //       _px_prime_prot_ph = _mom_corr_prot_th->Px() * (cos(DEG2RAD * _prot_phi_prime) / cos(DEG2RAD * _prot_phi));
  //       _py_prime_prot_ph = _mom_corr_prot_th->Py() * (sin(DEG2RAD * _prot_phi_prime) / sin(DEG2RAD * _prot_phi));
  //       _pz_prime_prot_ph = _mom_corr_prot_th->Pz();

  //       _mom_corr_prot_ph->SetXYZM(_px_prime_prot_ph, _py_prime_prot_ph, _pz_prime_prot_ph, MASS_P);
  //     }
  // }

  // for (size_t m = 0; m < Prot_mom_bins; m++) {
  //   double mom_min = min_prot_mom_values[m];
  //   double mom_max = max_prot_mom_values[m];
  //   if (_prot_mom > mom_min && _prot_mom < mom_max) {
  //     //For experimental data
  //     _prot_mom_prime = _prot_mom - prot_mom_corr[m] * alpha_prot_mom_corr;
  //     // // For simulation data
  //     // _prot_mom_prime = _prot_mom - prot_mom_corr_sim[m] * alpha_prot_mom_corr;

  //     _px_prime_prot_mom = _mom_corr_prot_ph->Px() * ((_prot_mom_prime) / (_prot_mom));
  //     _py_prime_prot_mom = _mom_corr_prot_ph->Py() * ((_prot_mom_prime) / (_prot_mom));
  //     _pz_prime_prot_mom = _mom_corr_prot_ph->Pz() * ((_prot_mom_prime) / (_prot_mom));
  //     _mom_corr_prot->SetXYZM(_px_prime_prot_mom, _py_prime_prot_mom, _pz_prime_prot_mom, MASS_P);
  //   }
  // }
}
// bool Reaction::ctof_prot() {
//   bool _prot_ctof = true;
//   _prot_ctof &= (4000 <= _prot_status && _prot_status < 6000);
//   return _prot_ctof;
// }
// bool Reaction::ftof_prot() {
//   bool _prot_ftof = true;
//   _prot_ftof &= (2000 <= _prot_status && _prot_status < 4000);
//   return _prot_ftof;
// }
void Reaction::SetPip(int i) {
  _numPip++;
  _numPos++;
  _hasPip = true;
  // _pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
  //   // _pip_status = abs(_data->status(i));
  _Energy_loss_uncorr_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
  _pip_mom_uncorr = _Energy_loss_uncorr_pip->P();
  _pip_theta = _Energy_loss_uncorr_pip->Theta() * 180 / PI;
  // std::cout << "pip ststus " << _data->status(i) << "   pip theta " << _pip_theta << " pip  mom   "
  //           << _pip_mom_uncorr<< std::endl;
  if (abs(_data->status(i)) < 4000) {
    if (_pip_theta <= 27) {
      _E_corr_val_pip = 9.21970527e-05 * pow(_pip_mom_uncorr, 3) - 3.70500143e-04 * pow(_pip_mom_uncorr, 2) +
                        2.78880101e-04 * (_pip_mom_uncorr) + 2.66040566e-03;

    } else {
      _E_corr_val_pip = -0.00010482 * pow(_pip_mom_uncorr, 3) + 0.00080463 * pow(_pip_mom_uncorr, 2) -
                        0.0022871 * (_pip_mom_uncorr) + 0.00831496;

    }
  } else if (abs(_data->status(i)) >= 4000) {
    _E_corr_val_pip = -0.00631413  * pow(_pip_mom_uncorr, 5) + 0.04713584  * pow(_pip_mom_uncorr, 4) -
                      0.12554256 * pow(_pip_mom_uncorr, 3) + 0.15622077 * pow(_pip_mom_uncorr, 2) -
                      0.11467851 * (_pip_mom_uncorr) + 0.01917004;


    // _E_corr_val_pip =  -0.00279293 * pow(_pip_mom_uncorr, 3) + 0.0206818 * pow(_pip_mom_uncorr, 2) -
    //                   0.05257802 * pow(_pip_mom_uncorr, 2) + 0.00996933;
  }

  _pip_mom_tmt = _pip_mom_uncorr + _E_corr_val_pip;  // first iteration
  _px_prime_pip_E = _data->px(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
  _py_prime_pip_E = _data->py(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
  _pz_prime_pip_E = _data->pz(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));

  _pip->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);

  // if (abs(_data->status(i)) < 4000) {
    
  //     _E_corr_val_pip_th = 0.00000000;
    
  // } else if (abs(_data->status(i)) >= 4000) {

  //   _E_corr_val_pip_th = -7.08389160e-11 * pow(_pip_theta, 5) + 3.75704402e-08 * pow(_pip_theta, 4) -
  //                        7.26740433e-06 * pow(_pip_theta, 3) + 6.45415606e-04 * pow(_pip_theta, 2) -
  //                        2.60057363e-02 * (_pip_theta) + 3.78387868e-01;
  // }
  // _pip_mom_tmt2 = _pip_mom_tmt + _E_corr_val_pip_th;  // theta iteration

  // _px_prime_pip_E_tmt = _data->px(i) * ((_pip_mom_tmt2) / (_pip_mom_uncorr));
  // _py_prime_pip_E_tmt = _data->py(i) * ((_pip_mom_tmt2) / (_pip_mom_uncorr));
  // _pz_prime_pip_E_tmt = _data->pz(i) * ((_pip_mom_tmt2) / (_pip_mom_uncorr));

  // _pip->SetXYZM(_px_prime_pip_E_tmt, _py_prime_pip_E_tmt, _pz_prime_pip_E_tmt, MASS_PIP);

  // second iterations

  // _pip_mom_tmt2 = _pip_tmt->P();  // for second iteration

  // // let's do second iteration for cd pip
  // if (abs(_data->status(i)) < 4000) {
  //   _E_corr_val_pip2 = 0.0;
  // } else if (abs(_data->status(i)) >= 4000) {
  //   _E_corr_val_pip2 = -0.00125164 * pow(_pip_mom_tmt2, 5) + 0.01272027 * pow(_pip_mom_tmt2, 4) -
  //                      0.04457356 * pow(_pip_mom_tmt2, 3) + 0.06272048 * pow(_pip_mom_tmt2, 2) -
  //                      0.03798534 * (_pip_mom_tmt2)-0.00716495;
  // }
  // _pip_mom = _pip_mom_tmt2 + _E_corr_val_pip2;
  // _px_prime_pip_E = _data->px(i) * ((_pip_mom) / (_pip_mom_tmt2));
  // _py_prime_pip_E = _data->py(i) * ((_pip_mom) / (_pip_mom_tmt2));
  // _pz_prime_pip_E = _data->pz(i) * ((_pip_mom) / (_pip_mom_tmt2));

  // _pip->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);

  // //   _pip_mom = _pip->P();
  // //   _pip_theta = _pip->Theta() * 180 / PI;

  //   if (_pip->Phi() > 0)
  //     _pip_phi = _pip->Phi() * 180 / PI;
  //   else if (_pip->Phi() < 0)
  //     _pip_phi = (_pip->Phi() + 2 * PI) * 180 / PI;

  //   for (size_t t = 0; t < Pip_theta_bins; t++) {
  //     double theta_min = min_pip_theta_values[t];
  //     double theta_max = max_pip_theta_values[t];
  //     if (_pip_theta > theta_min && _pip_theta < theta_max) {
  //       // For experimental data
  //       _pip_theta_prime = _pip_theta - pip_theta_corr[t] * alpha_pip_theta_corr;
  //       // // For simulation data
  //       // _pip_theta_prime = _pip_theta - pip_theta_corr_sim[t] * alpha_pip_theta_corr;

  //       _px_prime_pip_th = _data->px(i) * (sin(DEG2RAD * _pip_theta_prime) / sin(DEG2RAD * _pip_theta));
  //       _py_prime_pip_th = _data->py(i) * (sin(DEG2RAD * _pip_theta_prime) / sin(DEG2RAD * _pip_theta));
  //       _pz_prime_pip_th = _data->pz(i) * (cos(DEG2RAD * _pip_theta_prime) / cos(DEG2RAD * _pip_theta));
  //       _E_prime_pip_th = sqrt(abs(_px_prime_pip_th * _px_prime_pip_th + _py_prime_pip_th * _py_prime_pip_th +
  //                                  _pz_prime_pip_th * _pz_prime_pip_th));
  //       _mom_corr_pip_th->SetPxPyPzE(_px_prime_pip_th, _py_prime_pip_th, _pz_prime_pip_th, _E_prime_pip_th);
  //     }
  //   }

  //   for (size_t p = 0; p < Pip_phi_bins; p++) {
  //     double phi_min = min_pip_phi_values[p];
  //     double phi_max = max_pip_phi_values[p];
  //     if (_pip_phi > phi_min && _pip_phi < phi_max) {
  //       //For experimantal data
  //       _pip_phi_prime = _pip_phi - pip_phi_corr[p] * alpha_pip_phi_corr;
  // // For simulations data
  //       // _pip_phi_prime = _pip_phi - pip_phi_corr_sim[p] * alpha_pip_phi_corr;

  //       _px_prime_pip_ph = _mom_corr_pip_th->Px() * (cos(DEG2RAD * _pip_phi_prime) / cos(DEG2RAD * _pip_phi));
  //       _py_prime_pip_ph = _mom_corr_pip_th->Py() * (sin(DEG2RAD * _pip_phi_prime) / sin(DEG2RAD * _pip_phi));
  //       _pz_prime_pip_ph = _mom_corr_pip_th->Pz();

  //       _mom_corr_pip_ph->SetXYZM(_px_prime_pip_ph, _py_prime_pip_ph, _pz_prime_pip_ph, MASS_PIP);
  //     }
  //   }

  //   for (size_t m = 0; m < Pip_mom_bins; m++) {
  //     double mom_min = min_pip_mom_values[m];
  //     double mom_max = max_pip_mom_values[m];
  //     if (_pip_mom > mom_min && _pip_mom < mom_max) {
  //       // For experimantal data
  //       _pip_mom_prime = _pip_mom - pip_mom_corr[m] * alpha_pip_mom_corr;
  //       // // For simulation data
  //       // _pip_mom_prime = _pip_mom - pip_mom_corr_sim[m] * alpha_pip_mom_corr;

  //       _px_prime_pip_mom = _mom_corr_pip_ph->Px() * ((_pip_mom_prime) / (_pip_mom));
  //       _py_prime_pip_mom = _mom_corr_pip_ph->Py() * ((_pip_mom_prime) / (_pip_mom));
  //       _pz_prime_pip_mom = _mom_corr_pip_ph->Pz() * ((_pip_mom_prime) / (_pip_mom));
  //       _mom_corr_pip->SetXYZM(_px_prime_pip_mom, _py_prime_pip_mom, _pz_prime_pip_mom, MASS_PIP);
  //     }
  //   }
}

// bool Reaction::ctof_pip() {
//   bool _pip_ctof = true;
//   _pip_ctof &= (4000 <= _pip_status && _pip_status < 6000);
//   return _pip_ctof;
// }
// bool Reaction::ftof_pip() {
//   bool _pip_ftof = true;
//   _pip_ftof &= (2000 <= _pip_status && _pip_status < 4000);
//   return _pip_ftof;
// }

void Reaction::SetPim(int i) {
  _numPim++;
  _numNeg++;
  _hasPim = true;
  // _pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
  // // // _pim_status = abs(_data->status(i));
  _Energy_loss_uncorr_pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
  _pim_mom_uncorr = _Energy_loss_uncorr_pim->P();
  _pim_theta = _Energy_loss_uncorr_pim->Theta() * 180 / PI;
  // std::cout << "pim ststus " << _data->status(i) << "   pim theta " << _pim_theta << " pim  mom   "
  //           << _pim_mom_uncorr<< std::endl;
  if (abs(_data->status(i)) < 4000) {
    if (_pim_theta <= 27) {
      _E_corr_val_pim = -0.00035275 * pow(_pim_mom_uncorr, 3) + 0.00291237 * pow(_pim_mom_uncorr, 2) -
                        0.00681058 * (_pim_mom_uncorr) + 0.00736721;

    } else {
      _E_corr_val_pim = 0.00019358 * pow(_pim_mom_uncorr, 3) - 0.00103456 * pow(_pim_mom_uncorr, 2) +
                        0.00024772 * (_pim_mom_uncorr) + 0.00735159;

    }
  } else if (abs(_data->status(i)) >= 4000) {
    _E_corr_val_pim = (0.02153442) * pow(_pim_mom_uncorr, 5) -
                      (0.13271424) * pow(_pim_mom_uncorr, 4) +
                      (0.27140262) * pow(_pim_mom_uncorr, 3) -
                      (0.23266059) * pow(_pim_mom_uncorr, 2) +
                      (0.04031421) * (_pim_mom_uncorr) + 0.0036634;

  }
  // _pim_mom = _pim_mom_uncorr + _E_corr_val_pim; // first iteration

  _pim_mom_tmt = _pim_mom_uncorr + _E_corr_val_pim;  // first iteration

  _px_prime_pim_E = _data->px(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
  _py_prime_pim_E = _data->py(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
  _pz_prime_pim_E = _data->pz(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));

  _pim->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, MASS_PIM);

  // if (abs(_data->status(i)) < 4000) {

  //     _E_corr_val_pim_th = 0.00000000;
    
  // } else if (abs(_data->status(i)) >= 4000) {
  //   // -2.07141609e-10 8.81758359e-08 - 1.46534798e-05 1.17681655e-03 - 4.50634123e-02 6.54748237e-01;
  //   _E_corr_val_pim_th = (-2.07141609e-10) * pow(_pim_theta, 5) + (8.81758359e-08) * pow(_pim_theta, 4) +
  //                        (-1.46534798e-05) * pow(_pim_theta, 3) + (1.17681655e-03) * pow(_pim_theta, 2) +
  //                        (-4.50634123e-02) * (_pim_theta) + 6.54748237e-01;
  // }

  //     _pim_mom_tmt2 = _pim_mom_tmt + _E_corr_val_pim_th;

  // _px_prime_pim_E_tmt = _data->px(i) * ((_pim_mom_tmt2) / (_pim_mom_uncorr));
  // _py_prime_pim_E_tmt = _data->py(i) * ((_pim_mom_tmt2) / (_pim_mom_uncorr));
  // _pz_prime_pim_E_tmt = _data->pz(i) * ((_pim_mom_tmt2) / (_pim_mom_uncorr));

  // _pim->SetXYZM(_px_prime_pim_E_tmt, _py_prime_pim_E_tmt, _pz_prime_pim_E_tmt, MASS_PIM);

  // std::cout << "_E_corr_val_pim " << _E_corr_val_pim << "  _E_corr_val_pim_th " << _E_corr_val_pim_th
  //           << "   pim mom tmt  " << _pim_mom_tmt << "   pim mom tmt2  " << _pim_mom_tmt2 << " diff "
  //           << _pim_mom_tmt - _pim_mom_tmt2 << std::endl;

  // _pim_tmt->SetXYZM(_px_prime_pim_E_tmt, _py_prime_pim_E_tmt, _pz_prime_pim_E_tmt, MASS_PIM);

  // _pim_mom_tmt2 = _pim_tmt->P();  // for second iteration

  // // std::cout << " diff  " << _pim_tmt->P() - _pim_mom_tmt << std::endl;
  // // let's do second iteration for cd pim
  // if (abs(_data->status(i)) < 4000) {
  //   _E_corr_val_pim2 = 0.0;
  // } else if (abs(_data->status(i)) >= 4000) {
  //   _E_corr_val_pim2 = 0.07604229 * pow(_pim_mom_tmt2, 7) - 0.69056865 * pow(_pim_mom_tmt2, 6) +
  //                      2.42244641 * pow(_pim_mom_tmt2, 5) - 4.26630462 * pow(_pim_mom_tmt2, 4) +
  //                      4.07033382 * pow(_pim_mom_tmt2, 3) - 2.09075715 * pow(_pim_mom_tmt2, 2) +
  //                      0.52748137 * (_pim_mom_tmt2)-0.04274812;
  // }
  // _pim_mom = _pim_mom_tmt2 + _E_corr_val_pim2;
  // _px_prime_pim_E = _data->px(i) * ((_pim_mom) / (_pim_mom_tmt2));
  // _py_prime_pim_E = _data->py(i) * ((_pim_mom) / (_pim_mom_tmt2));
  // _pz_prime_pim_E = _data->pz(i) * ((_pim_mom) / (_pim_mom_tmt2));

  // _pim->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, MASS_PIM);

  // // _pim_mom = _pim->P();
  // // _pim_theta = _pim->Theta() * 180 / PI;

  // if (_pim->Phi() > 0)
  //   _pim_phi = _pim->Phi()* 180 / PI;
  // else if (_pim->Phi() < 0)
  //   _pim_phi = (_pim->Phi() + 2 * PI)* 180 / PI;

  // for (size_t t = 0; t < Pim_theta_bins; t++) {
  //   double theta_min = min_pim_theta_values[t];
  //   double theta_max = max_pim_theta_values[t];
  //   if (_pim_theta > theta_min && _pim_theta < theta_max) {
  //     //For experimental data
  //     _pim_theta_prime = _pim_theta - pim_theta_corr[t] * alpha_pim_theta_corr;

  //     // // For simulations data
  //     // _pim_theta_prime = _pim_theta - pim_theta_corr_sim[t] * alpha_pim_theta_corr;

  //     _px_prime_pim_th = _data->px(i) * (sin(DEG2RAD * _pim_theta_prime) / sin(DEG2RAD * _pim_theta));
  //     _py_prime_pim_th = _data->py(i) * (sin(DEG2RAD * _pim_theta_prime) / sin(DEG2RAD * _pim_theta));
  //     _pz_prime_pim_th = _data->pz(i) * (cos(DEG2RAD * _pim_theta_prime) / cos(DEG2RAD * _pim_theta));
  //     _E_prime_pim_th =
  //         sqrt(abs(_px_prime_pim_th * _px_prime_pim_th + _py_prime_pim_th* _py_prime_pim_th + _pz_prime_pim_th *
  //         _pz_prime_pim_th));
  //     _mom_corr_pim_th->SetPxPyPzE(_px_prime_pim_th, _py_prime_pim_th, _pz_prime_pim_th, _E_prime_pim_th);
  //   }
  // }

  // for (size_t p = 0; p < Pim_phi_bins; p++) {
  //   double phi_min = min_pim_phi_values[p];
  //   double phi_max = max_pim_phi_values[p];
  //   if (_pim_phi > phi_min && _pim_phi < phi_max) {
  //     // For experimantal data
  //     _pim_phi_prime = _pim_phi - pim_phi_corr[p] * alpha_pim_phi_corr;

  //     // // For simulations data
  //     // _pim_phi_prime = _pim_phi - pim_phi_corr_sim[p] * alpha_pim_phi_corr;

  //     _px_prime_pim_ph = _mom_corr_pim_th->Px() * (cos(DEG2RAD * _pim_phi_prime) / cos(DEG2RAD * _pim_phi));
  //     _py_prime_pim_ph = _mom_corr_pim_th->Py() * (sin(DEG2RAD * _pim_phi_prime) / sin(DEG2RAD * _pim_phi));
  //     _pz_prime_pim_ph = _mom_corr_pim_th->Pz();

  //     _mom_corr_pim_ph->SetXYZM(_px_prime_pim_ph, _py_prime_pim_ph, _pz_prime_pim_ph, MASS_PIM);
  //   }
  // }

  // for (size_t m = 0; m < Pim_mom_bins; m++) {
  //   double mom_min = min_pim_mom_values[m];
  //   double mom_max = max_pim_mom_values[m];
  //   if (_pim_mom > mom_min && _pim_mom < mom_max) {
  //     // For experimantal data
  //     _pim_mom_prime = _pim_mom - pim_mom_corr[m] * alpha_pim_mom_corr;

  //     // // For simulations data
  //     // _pim_mom_prime = _pim_mom - pim_mom_corr_sim[m] * alpha_pim_mom_corr;

  //     _px_prime_pim_mom = _mom_corr_pim_ph->Px() * ((_pim_mom_prime) / (_pim_mom));
  //     _py_prime_pim_mom = _mom_corr_pim_ph->Py() * ((_pim_mom_prime) / (_pim_mom));
  //     _pz_prime_pim_mom = _mom_corr_pim_ph->Pz() * ((_pim_mom_prime) / (_pim_mom));
  //     _mom_corr_pim->SetXYZM(_px_prime_pim_mom, _py_prime_pim_mom, _pz_prime_pim_mom, MASS_PIM);
  //   }
  // }
}
// bool Reaction::ctof_pim() {
//   bool _pim_ctof = true;
//   _pim_ctof &= (4000 <= _pim_status && _pim_status < 6000);
//   return _pim_ctof;
// }
// bool Reaction::ftof_pim() {
//   bool _pim_ftof = true;
//   _pim_ftof &= (2000 <= _pim_status && _pim_status < 4000);
//   return _pim_ftof;
// }

// float Reaction::rec_pim_px() {
//   return _beam->Px() - _elec->Px() + _target->Px() - _pip->Px() - _prot->Px() - _pim->Px();
// }
// float Reaction::rec_pim_py() {
//   return _beam->Py() - _elec->Py() + _target->Py() - _pip->Py() - _prot->Py() - _pim->Py();
// }
// float Reaction::rec_pim_pz() {
//   return _beam->Pz() - _elec->Pz() + _target->Pz() - _pip->Pz() - _prot->Pz() - _pim->Pz();
// }
// float Reaction::rec_pim_E() { return _beam->E() - _elec->E() + _target->E() - _pip->E() - _prot->E() -
// _pim->E(); } float Reaction::rec_pim_P() {
//   return sqrt(abs(pow((_beam->Px() - _elec->Px() + _target->Px() - _pip->Px() - _prot->Px() - _pim->Px()),
//   2)
//   +
//                   pow((_beam->Py() - _elec->Py() + _target->Py() - _pip->Py() - _prot->Py() - _pim->Py()),
//                   2)
//                   + pow((_beam->Pz() - _elec->Pz() + _target->Pz() - _pip->Pz() - _prot->Pz() -
//                   _pim->Pz()), 2)));
// }

// float Reaction::rec_pim_mm2() {
//   return abs(pow(_beam->Px() - _elec->Px() + _target->Px() - _pip->Px() - _prot->Px() - _pim->Px(), 2) +
//              pow(_beam->Py() - _elec->Py() + _target->Py() - _pip->Py() - _prot->Py() - _pim->Py(), 2) +
//              pow(_beam->Pz() - _elec->Pz() + _target->Pz() - _pip->Pz() - _prot->Pz() - _pim->Pz(), 2) -
//              pow(_beam->E() - _elec->E() + _target->E() - _pip->E() - _prot->E() - _pim->E(), 2));
// }

// float Reaction::beam_px() { return _beam->Px(); }
// float Reaction::beam_py() { return _beam->Py(); }
// float Reaction::beam_pz() { return _beam->Pz(); }
// float Reaction::beam_E() { return _beam->E(); }

// float Reaction::elec_px() { return _elec->Px(); }
// float Reaction::elec_py() { return _elec->Py(); }
// float Reaction::elec_pz() { return _elec->Pz(); }
// float Reaction::elec_E() { return _elec->E(); }

// float Reaction::target_px() { return _target->Px(); }
// float Reaction::target_py() { return _target->Py(); }
// float Reaction::target_pz() { return _target->Pz(); }
// float Reaction::target_E() { return _target->E(); }

// float Reaction::pim_px() { return _pim->Px(); }
// float Reaction::pim_py() { return _pim->Py(); }
// float Reaction::pim_pz() { return _pim->Pz(); }
// float Reaction::pim_E() { return _pim->E(); }
// float Reaction::pim_P() { return _pim->P(); }

// float Reaction::pip_px() { return _pip->Px(); }
// float Reaction::pip_py() { return _pip->Py(); }
// float Reaction::pip_pz() { return _pip->Pz(); }
// float Reaction::pip_E() { return _pip->E(); }

// float Reaction::prot_px() { return _prot->Px(); }
// float Reaction::prot_py() { return _prot->Py(); }
// float Reaction::prot_pz() { return _prot->Pz(); }
// float Reaction::prot_E() { return _prot->E(); }

void Reaction::SetNeutron(int i) {
  _numNeutral++;
  _hasNeutron = true;
  _neutron->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_N);
}

void Reaction::SetOther(int i) {
  if (_data->pid(i) == NEUTRON) {
    SetNeutron(i);
  } else {
    _numOther++;
    _hasOther = true;
    _other->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), mass[_data->pid(i)]);
  }
}

void Reaction::CalcMissMass() {
  auto mm = std::make_unique<TLorentzVector>();
  auto mm_mpip = std::make_unique<TLorentzVector>();
  auto mm_mprot = std::make_unique<TLorentzVector>();
  auto mm_excl = std::make_unique<TLorentzVector>();

  *mm += (*_gamma + *_target);

  if (TwoPion_missingPim()) {
    *mm -= *_prot;
    *mm -= *_pip;
    // *mm -= *_pim;
    _MM = mm->M();
    _MM2 = mm->M2();

    // _rec_pim_mom = mm->P();
    // _rec_pim_theta = mm->Theta() * 180 / PI;

    // if (mm->Phi() >= 0)
    //   _rec_pim_phi = (mm->Phi() * 180 / PI);
    // else if (mm->Phi() < 0)
    //   _rec_pim_phi = ((mm->Phi() + 2 * PI) * 180 / PI);

  //   // // // _x_mu_E = mm->E();
  //   // // // _x_mu_P = mm->P();
  //   // // // _x_mu_Px = mm->Px();
  //   // // // _x_mu_Py = mm->Py();
  //   // // // _x_mu_Pz = mm->Pz();
  //   // // // _x_mu_theta = mm->Theta() * RAD2DEG;
  //   // // // _x_mu_m2 = mm->E() * mm->E() - mm->P() * mm->P();
  //   // // // _x_mu_m = mm->E() - mm->P();
  //   // // //   //
  }
  if (TwoPion_exclusive()) {
    // *mm -= *_prot;
    // *mm -= *_pip;
    // // *mm -= *_pim;
    // _MM = mm->M();
    // _MM2 = mm->M2();

    *mm_excl += (*_gamma + *_target);
    *mm_excl -= *_prot;
    *mm_excl -= *_pip;
    *mm_excl -= *_pim;

    _MM2_exclusive = mm_excl->M2();
    _excl_Energy = mm_excl->E();

    // _rec_pim_mom = mm->P();
    // _rec_pim_theta = mm->Theta() * 180 / PI;

    // if (mm->Phi() >= 0)
    //   _rec_pim_phi = (mm->Phi() * 180 / PI);
    // else if (mm->Phi() < 0)
    //   _rec_pim_phi = ((mm->Phi() + 2 * PI) * 180 / PI);

    // // //   // //////// for x_mu - elec/beam theta phi
    // // //   // if (mm_excl->Phi() >= 0)
    // // //   //   _x_mu_phi = (mm_excl->Phi() * 180 / PI);
    // // //   // else if (mm_excl->Phi() < 0)
    // // //   //   _x_mu_phi = ((mm_excl->Phi() + 2 * PI) * 180 / PI);

    // // //   // if (_elec->Phi() >= 0)
    // // //   //   _elec_phi = (_elec->Phi() * 180 / PI);
    // // //   // else if (_elec->Phi() < 0)
    // // //   //   _elec_phi = ((_elec->Phi() + 2 * PI) * 180 / PI);

    // // //   // if (_beam->Phi() >= 0)
    // // //   //   _beam_phi = (_beam->Phi() * 180 / PI);
    // // //   // else if (_beam->Phi() < 0)
    // // //   //   _beam_phi = ((_beam->Phi() + 2 * PI) * 180 / PI);

    // // //   // _diff_elec_x_mu_theta = (_elec->Theta() * 180 / PI);  // - (mm_excl->Theta() * 180 / PI);
    // // //   // _diff_elec_x_mu_phi = (_elec_phi - _x_mu_phi);

    // // //   // _diff_beam_x_mu_theta = (_beam->Theta() * 180 / PI);  //-(mm_excl->Theta() * 180 / PI);
    // // //   // _diff_beam_x_mu_phi = (_beam_phi - _x_mu_phi);

    // // //   // // std::cout << " beam_theta " << _diff_beam_x_mu_theta << std::endl;
    // // //   // // std::cout << " rec_pim_energy " << mm->E() << std::endl;

    // //   // for mPip peak with exclusive events
    *mm_mpip += (*_gamma + *_target);
    *mm_mpip -= *_prot;
    *mm_mpip -= *_pim;
    _MM2_mPip = mm_mpip->M2();

    // //   // for mProt peak with exclusive events
    *mm_mprot += (*_gamma + *_target);
    *mm_mprot -= *_pip;
    *mm_mprot -= *_pim;
    _MM2_mProt = mm_mprot->M2();
  }
  // if (TwoPion_missingPip()) {
  //   *mm_mpip += (*_gamma + *_target);
  //   *mm_mpip -= *_prot;
  //   *mm_mpip -= *_pim;
  //   _MM2_mPip = mm_mpip->M2();
  // }
  // if (TwoPion_missingProt()) {
  //   *mm_mprot += (*_gamma + *_target);
  //   *mm_mprot -= *_pip;
  //   *mm_mprot -= *_pim;
  //   _MM2_mProt = mm_mprot->M2();
  // }
}
// float Reaction::Diff_elec_x_mu_theta() {
//   if (_diff_elec_x_mu_theta != _diff_elec_x_mu_theta) CalcMissMass();
//   return _diff_elec_x_mu_theta;
// }

// float Reaction::Diff_elec_x_mu_phi() {
//   if (_diff_elec_x_mu_phi != _diff_elec_x_mu_phi) CalcMissMass();
//   return _diff_elec_x_mu_phi;
// }

// float Reaction::Diff_beam_x_mu_theta() {
//   if (_diff_beam_x_mu_theta != _diff_beam_x_mu_theta) CalcMissMass();
//   return _diff_beam_x_mu_theta;
// }

// float Reaction::Diff_beam_x_mu_phi() {
//   if (_diff_beam_x_mu_phi != _diff_beam_x_mu_phi) CalcMissMass();
//   return _diff_beam_x_mu_phi;
// }

float Reaction::MM() {
  if (_MM != _MM) CalcMissMass();
  return _MM;
}
float Reaction::MM2() {
  if (_MM2 != _MM2) CalcMissMass();
  return _MM2;
}
float Reaction::MM2_exclusive() {
  if (_MM2_exclusive != _MM2_exclusive) CalcMissMass();
  return _MM2_exclusive;
}
float Reaction::MM2_mPip() {
  if (_MM2_mPip != _MM2_mPip) CalcMissMass();
  return _MM2_mPip;
}
float Reaction::MM2_mProt() {
  if (_MM2_mProt != _MM2_mProt) CalcMissMass();
  return _MM2_mProt;
}
float Reaction::Energy_excl() {
  if (_excl_Energy != _excl_Energy) CalcMissMass();
  //  std::cout << "_x_mu_p  " << _x_mu->E() << '\n';
  //  if (_x_mu_E > 0)
  return _excl_Energy;
  // else
  // return NAN;
}
float Reaction::pim_momentum() {
  // if (_rec_pim_mom != _rec_pim_mom) CalcMissMass();

  if (TwoPion_missingPim()) {
  // if (TwoPion_exclusive()) {
    auto missingpim_ = std::make_unique<TLorentzVector>();
    *missingpim_ += *_gamma + *_target - *_prot - *_pip;
    // *missingpim_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pip;

    return missingpim_->P();
    // return _rec_pim_mom;

  } else
    return NAN;
}
float Reaction::pim_theta_lab() {
  // if (_rec_pim_theta != _rec_pim_theta) CalcMissMass();

  if (TwoPion_missingPim()) {
  // if (TwoPion_exclusive()) {
    auto missingpim_ = std::make_unique<TLorentzVector>();
    *missingpim_ += *_gamma + *_target - *_prot - *_pip;
    // *missingpim_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pip;

    return missingpim_->Theta() * 180.0 / PI;
    // return _rec_pim_theta;
  } else
    return NAN;
}
float Reaction::pim_Phi_lab() {
  // if (_rec_pim_phi != _rec_pim_phi) CalcMissMass();

  if (TwoPion_missingPim()) {
  // if (TwoPion_exclusive()) {
    auto missingpim_ = std::make_unique<TLorentzVector>();
    *missingpim_ += *_gamma + *_target - *_prot - *_pip;
    // *missingpim_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pip;

    if (missingpim_->Phi() > 0)
      return missingpim_->Phi() * 180 / PI;
    else if (missingpim_->Phi() < 0)
      return (missingpim_->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
    // return _rec_pim_phi;
  } else
    return NAN;
}
float Reaction::pim_momentum_measured() {
  if (TwoPion_exclusive())
    return _pim->P();
  else
    return NAN;
}

float Reaction::pim_theta_lab_measured() {
  if (TwoPion_exclusive())
    return _pim->Theta() * 180.0 / PI;
  else
    return NAN;
}

float Reaction::pim_Phi_lab_measured() {
  if (TwoPion_exclusive()) {
    if (_pim->Phi() > 0)
      return _pim->Phi() * 180 / PI;
    else if (_pim->Phi() < 0)
      return (_pim->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
// float Reaction::pim_momentum_corrected() {
//   if (TwoPion_exclusive())
//     return _mom_corr_pim->P();
//   else
//     return NAN;
// }

// float Reaction::pim_theta_corrected() {
//   if (TwoPion_exclusive())
//     return _mom_corr_pim->Theta() * 180.0 / PI;
//   else
//     return NAN;
// }

// float Reaction::pim_Phi_corrected() {
//   if (TwoPion_exclusive()) {
//     if (_mom_corr_pim->Phi() > 0)
//       return _mom_corr_pim->Phi() * 180 / PI;
//     else if (_mom_corr_pim->Phi() < 0)
//       return (_mom_corr_pim->Phi() + 2 * PI) * 180 / PI;
//     else
//       return NAN;
//   } else
//     return NAN;
// }
////////////////mPip
float Reaction::pip_momentum() {
  // if (TwoPion_missingPip()) {
  if (TwoPion_exclusive()) {
    auto missingpip_ = std::make_unique<TLorentzVector>();
    *missingpip_ += *_gamma + *_target - *_prot - *_pim;
    // *missingpip_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pim;

    return missingpip_->P();
  } else
    return NAN;
}
float Reaction::pip_theta_lab() {
  // if (TwoPion_missingPip()) {
  if (TwoPion_exclusive()) {
    auto missingpip_ = std::make_unique<TLorentzVector>();
    *missingpip_ += *_gamma + *_target - *_prot - *_pim;
    // *missingpip_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pim;
    return missingpip_->Theta() * 180.0 / PI;
  } else
    return NAN;
}
float Reaction::pip_Phi_lab() {
  // if (TwoPion_missingPip()) {
  if (TwoPion_exclusive()) {
    auto missingpip_ = std::make_unique<TLorentzVector>();
    *missingpip_ += *_gamma + *_target - *_prot - *_pim;
    // *missingpip_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pim;

    if (missingpip_->Phi() > 0)
      return missingpip_->Phi() * 180 / PI;
    else if (missingpip_->Phi() < 0)
      return (missingpip_->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
float Reaction::pip_momentum_measured() {
  if (TwoPion_exclusive())
    return _pip->P();
  else
    return NAN;
}

float Reaction::pip_theta_lab_measured() {
  if (TwoPion_exclusive())
    return _pip->Theta() * 180.0 / PI;
  else
    return NAN;
}

float Reaction::pip_Phi_lab_measured() {
  if (TwoPion_exclusive()) {
    if (_pip->Phi() > 0)
      return _pip->Phi() * 180 / PI;
    else if (_pip->Phi() < 0)
      return (_pip->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}

// float Reaction::pip_momentum_corrected() {
//   if (TwoPion_exclusive())
//     return _mom_corr_pip->P();
//   else
//     return NAN;
// }
// float Reaction::pip_theta_corrected() {
//   if (TwoPion_exclusive())
//     return _mom_corr_pip->Theta() * 180.0 / PI;
//   else
//     return NAN;
// }

// float Reaction::pip_Phi_corrected() {
//   if (TwoPion_exclusive()) {
//     if (_mom_corr_pip->Phi() > 0)
//       return _mom_corr_pip->Phi() * 180 / PI;
//     else if (_mom_corr_pip->Phi() < 0)
//       return (_mom_corr_pip->Phi() + 2 * PI) * 180 / PI;
//     else
//       return NAN;
//   } else
//     return NAN;
// }

////////////////mProt
float Reaction::prot_momentum() {
  // if (TwoPion_missingProt()) {
  if (TwoPion_exclusive()) {
    auto missingprot_ = std::make_unique<TLorentzVector>();
    *missingprot_ += *_gamma + *_target - *_pip - *_pim;
    // *missingprot_ += *_gamma + *_target - *_mom_corr_pip - *_mom_corr_pim;

    return missingprot_->P();
  } else
    return NAN;
}
float Reaction::prot_theta_lab() {
  // if (TwoPion_missingProt()) {
  if (TwoPion_exclusive()) {
    auto missingprot_ = std::make_unique<TLorentzVector>();
    *missingprot_ += *_gamma + *_target - *_pip - *_pim;
    // *missingprot_ += *_gamma + *_target - *_mom_corr_pip - *_mom_corr_pim;

    return missingprot_->Theta() * 180.0 / PI;
  } else
    return NAN;
}
float Reaction::prot_Phi_lab() {
  // if (TwoPion_missingProt()) {
  if (TwoPion_exclusive()) {
    auto missingprot_ = std::make_unique<TLorentzVector>();
    *missingprot_ += *_gamma + *_target - *_pip - *_pim;
    // *missingprot_ += *_gamma + *_target - *_mom_corr_pip - *_mom_corr_pim;

    if (missingprot_->Phi() > 0)
      return missingprot_->Phi() * 180 / PI;
    else if (missingprot_->Phi() < 0)
      return (missingprot_->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
float Reaction::prot_momentum_measured() {
  if (TwoPion_exclusive())
    return _prot->P();
  else
    return NAN;
}

float Reaction::prot_theta_lab_measured() {
  if (TwoPion_exclusive())
    return _prot->Theta() * 180.0 / PI;
  else
    return NAN;
}

float Reaction::prot_Phi_lab_measured() {
  if (TwoPion_exclusive()) {
    if (_prot->Phi() > 0)
      return _prot->Phi() * 180 / PI;
    else if (_prot->Phi() < 0)
      return (_prot->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}

// float Reaction::prot_momentum_corrected() {
//   if (TwoPion_exclusive())
//     return _mom_corr_prot->P();
//   else
//     return NAN;
// }
// float Reaction::prot_theta_corrected() {
//   if (TwoPion_exclusive())
//     return _mom_corr_prot->Theta() * 180.0 / PI;
//   else
//     return NAN;
// }

// float Reaction::prot_Phi_corrected() {
//   if (TwoPion_exclusive()) {
//     if (_mom_corr_prot->Phi() > 0)
//       return _mom_corr_prot->Phi() * 180 / PI;
//     else if (_mom_corr_prot->Phi() < 0)
//       return (_mom_corr_prot->Phi() + 2 * PI) * 180 / PI;
//     else
//       return NAN;
//   } else
//     return NAN;
// }

/////////////////////
std::string Reaction::CsvHeader() { return "e_rec_p,e_rec_theta,e_rec_phi,e_sec\n"; }
std::string Reaction::ReacToCsv() {
  // e_rec_p,e_rec_theta,e_rec_phi,e_sec
  std::string out = "";
  out += std::to_string(_elec->P()) + ",";
  out += std::to_string(_elec->Theta()) + ",";
  out += std::to_string(_elec->Phi()) + ",";
  out += std::to_string(_sector) + "\n";

  return out;
}

void Reaction::boost() {
  _is_boosted = true;
  _boosted_prot = std::make_unique<TLorentzVector>(*_prot);
  _boosted_pip = std::make_unique<TLorentzVector>(*_pip);
  _boosted_pim = std::make_unique<TLorentzVector>(*_gamma + *_target - *_prot - *_pip);  //(*_pim);
  _boosted_gamma = std::make_unique<TLorentzVector>(*_gamma);
  _boosted_pim_measured = std::make_unique<TLorentzVector>(*_pim);

  _rotated_prot = std::make_unique<TLorentzVector>(*_prot);
  _rotated_pip = std::make_unique<TLorentzVector>(*_pip);
  _rotated_pim = std::make_unique<TLorentzVector>(*_gamma + *_target - *_prot - *_pip);  //(*_pim);
  _rotated_pim_measured = std::make_unique<TLorentzVector>(*_pim);

  TRotation rot;
  _boosted_gamma->Transform(rot);
  float_t beta_1 = ((sqrt(_boosted_gamma->E() * _boosted_gamma->E() + _Q2)) / (_boosted_gamma->E() + MASS_P));
  TVector3 uz = _boosted_gamma->Vect().Unit();                  // uit vector along virtual photon
  TVector3 ux = ((_beam->Vect()).Cross(_elec->Vect())).Unit();  // unit vector along e cross e'
  ux.Rotate(3. * PI / 2, uz);                                   // rotating ux by 3pi/2 with uz as axis of roration
  rot.SetZAxis(uz, ux).Invert();                                // setting TRotation rot

  _boosted_prot->Transform(rot);
  _rotated_prot->Transform(rot);
  _boosted_prot->Boost(0, 0, -beta_1);

  _boosted_pip->Transform(rot);
  _rotated_pip->Transform(rot);
  _boosted_pip->Boost(0, 0, -beta_1);

  _boosted_pim->Transform(rot);
  _rotated_pim->Transform(rot);
  _boosted_pim->Boost(0, 0, -beta_1);

  _boosted_gamma->Boost(0, 0, -beta_1);

  _boosted_pim_measured->Transform(rot);
  _rotated_pim_measured->Transform(rot);
  _boosted_pim_measured->Boost(0, 0, -beta_1);
  // -beta ko value (0.5 to -0.5 huda
  // samma value aauchha nattra aaudyna)

  _prot_Vect3 = _boosted_prot->Vect();
  _pip_Vect3 = _boosted_pip->Vect();
  _pim_Vect3 = _boosted_pim_measured->Vect();
}

float_t Reaction::scalar_triple_product() {
  if (!_is_boosted) boost();
  if (TwoPion_exclusive()) {
    return (_prot_Vect3.Dot(_pip_Vect3.Cross(_pim_Vect3)));

  } else
    return NAN;
}

// // float Reaction::pim_momentum_cm() {
// //         if (!_is_boosted)
// //                 boost();
// //         if (TwoPion_missingPim())
// //                 return _boosted_pim->P();
// //         else
// //                 return NAN;
// // }

// float Reaction::pim_theta_cm() {
//   if (!_is_boosted) boost();
//   if (TwoPion_missingPim())
//     return _rotated_pim->Theta() * 180.0 / PI;
//   else
//     return NAN;
// }

// float Reaction::pim_Phi_cm() {
//   if (!_is_boosted) boost();
//   if (TwoPion_missingPim()) {
//     if (_rotated_pim->Phi() > 0)
//       return _rotated_pim->Phi() * 180 / PI;
//     else if (_rotated_pim->Phi() < 0)
//       return (_rotated_pim->Phi() + 2 * PI) * 180 / PI;
//     else
//       return NAN;
//   } else
//     return NAN;
// }

// // float Reaction::pim_momentum_cm_measured() {
// //         if (!_is_boosted)
// //                 boost();
// //         if (TwoPion_exclusive())
// //                 return _boosted_pim_measured->P();
// //         else
// //                 return NAN;
// // }

// float Reaction::pim_theta_cm_measured() {
//   if (!_is_boosted) boost();
//   if (TwoPion_exclusive())
//     return _rotated_pim_measured->Theta() * 180.0 / PI;
//   else
//     return NAN;
// }

// float Reaction::pim_Phi_cm_measured() {
//   if (!_is_boosted) boost();
//   if (TwoPion_exclusive()) {
//     if (_rotated_pim_measured->Phi() > 0)
//       return _rotated_pim_measured->Phi() * 180 / PI;
//     else if (_rotated_pim_measured->Phi() < 0)
//       return (_rotated_pim_measured->Phi() + 2 * PI) * 180 / PI;
//     else
//       return NAN;
//   } else
//     return NAN;
// }

MCReaction::MCReaction(const std::shared_ptr<Branches12>& data, float beam_enrgy) {
  _data = data;
  if (!_data->mc()) _data->mc_branches();
  _beam = std::make_unique<TLorentzVector>();
  _beam_energy = beam_enrgy;
  _weight_mc = _data->mc_weight();
  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  //_gamma = std::make_unique<TLorentzVector>();  // do i need this?
  _gamma_mc = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  //_elec = std::make_unique<TLorentzVector>();  // do i need this?
  _elec_mc = std::make_unique<TLorentzVector>();
  // this->SetElec();  // do i need this?
  this->SetMCElec();
  _prot_mc = std::make_unique<TLorentzVector>();
  _pip_mc = std::make_unique<TLorentzVector>();
  _pim_mc = std::make_unique<TLorentzVector>();
  //_other = std::make_unique<TLorentzVector>();  // do i need this?
  _other_mc = std::make_unique<TLorentzVector>();
  //_neutron = std::make_unique<TLorentzVector>();
}
// Reaction::~Reaction() {} // why this is not here
void MCReaction::SetMCElec() {
  //  _hasE = true;  //??
  _elec_mc->SetXYZM(_data->mc_px(0), _data->mc_py(0), _data->mc_pz(0), MASS_E);

  *_gamma_mc += *_beam - *_elec_mc;

  // Can calculate W and Q2 here
  _W_mc = physics::W_calc(*_beam, *_elec_mc);
  _Q2_mc = physics::Q2_calc(*_beam, *_elec_mc);
}

void MCReaction::SetMCProton(int i) { _prot_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_P); }

void MCReaction::SetMCPip(int i) { _pip_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIP); }

void MCReaction::SetMCPim(int i) { _pim_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIM); }
// void MCReaction::SetMCOther(int i) {
//   _other_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i),
//   mass[_data->pid(i)]);
// }

float MCReaction::pim_mom_mc_gen() {
  // if (Reaction::TwoPion_exclusive())
  return _pim_mc->P();
  // else
  //   return NAN;
}
float MCReaction::pip_mom_mc_gen() {
  // if (Reaction::TwoPion_exclusive())
  return _pip_mc->P();
  // else
  //   return NAN;
}
float MCReaction::prot_mom_mc_gen() {
  // if (Reaction::TwoPion_exclusive())
  return _prot_mc->P();
  // else
  //   return NAN;
}

float MCReaction::pim_theta_mc_gen() {
  if (Reaction::TwoPion_exclusive())
    return _pim_mc->Theta() * 180 / PI;
  else
    return NAN;
}
float MCReaction::pip_theta_mc_gen() {
  if (Reaction::TwoPion_exclusive())
    return _pip_mc->Theta() * 180 / PI;
  else
    return NAN;
}
float MCReaction::prot_theta_mc_gen() {
  if (TwoPion_exclusive())
    return _prot_mc->Theta() * 180 / PI;
  else
    return NAN;
}

float MCReaction::pim_phi_mc_gen() {
  if (_pim_mc->Phi() >= 0)
    return (_pim_mc->Phi() * 180 / PI);
  else if (_pim_mc->Phi() < 0)
    return ((_pim_mc->Phi() + 2 * PI) * 180 / PI);
  else
    return NAN;
}
float MCReaction::pip_phi_mc_gen() {
  if (_pip_mc->Phi() >= 0)
    return (_pip_mc->Phi() * 180 / PI);
  else if (_pip_mc->Phi() < 0)
    return ((_pip_mc->Phi() + 2 * PI) * 180 / PI);
  else
    return NAN;
}
float MCReaction::prot_phi_mc_gen() {
  if (_prot_mc->Phi() >= 0)
    return (_prot_mc->Phi() * 180 / PI);
  else if (_prot_mc->Phi() < 0)
    return ((_prot_mc->Phi() + 2 * PI) * 180 / PI);
  else
    return NAN;
}

// void MCReaction::CalcMissMass_mc() {
//   auto mm_excl_mc = std::make_unique<TLorentzVector>();

//   *mm_excl_mc += (*_gamma_mc + *_target);
//   *mm_excl_mc -= *_prot_mc;
//   *mm_excl_mc -= *_pip_mc;
//   *mm_excl_mc -= *_pim_mc;
//   _MM2_exclusive_mc = mm_excl_mc->M2();
//   _excl_Energy_mc = mm_excl_mc->E();

// _rec_x_mu_mom_mc = mm_excl_mc->P();
// _rec_x_mu_theta_mc = mm_excl_mc->Theta() * 180 / PI;

// if (mm_excl_mc->Phi() >= 0)
//   _x_mu_phi_mc = (mm_excl_mc->Phi() * 180 / PI);
// else if (mm_excl_mc->Phi() < 0)
//   _x_mu_phi_mc = ((mm_excl_mc->Phi() + 2 * PI) * 180 / PI);

// if (_elec_mc->Phi() >= 0)
//   _elec_phi_mc = (_elec_mc->Phi() * 180 / PI);
// else if (_elec_mc->Phi() < 0)
//   _elec_phi_mc = ((_elec_mc->Phi() + 2 * PI) * 180 / PI);

// if (_beam->Phi() >= 0)
//   _beam_phi_mc = (_beam->Phi() * 180 / PI);
// else if (_beam->Phi() < 0)
//   _beam_phi_mc = ((_beam->Phi() + 2 * PI) * 180 / PI);

// _diff_elec_x_mu_theta_mc = (_elec_mc->Theta() * 180 / PI) - (mm_excl_mc->Theta() * 180 / PI);
// _diff_elec_x_mu_phi_mc = (_elec_phi_mc - _x_mu_phi_mc);

// _diff_beam_x_mu_theta_mc = (mm_excl_mc->Theta() * 180 / PI);
// _diff_beam_x_mu_phi_mc = (_beam_phi_mc - _x_mu_phi_mc);
// }

// float MCReaction::Diff_elec_x_mu_theta_mc() {
//   if (_diff_elec_x_mu_theta_mc != _diff_elec_x_mu_theta_mc) CalcMissMass_mc();
//   return _diff_elec_x_mu_theta_mc;
// }

// float MCReaction::Diff_elec_x_mu_phi_mc() {
//   if (_diff_elec_x_mu_phi_mc != _diff_elec_x_mu_phi_mc) CalcMissMass_mc();
//   return _diff_elec_x_mu_phi_mc;
// }

// float MCReaction::Diff_beam_x_mu_theta_mc() {
//   if (_diff_beam_x_mu_theta_mc != _diff_beam_x_mu_theta_mc) CalcMissMass_mc();
//   return _diff_beam_x_mu_theta_mc;
// }

// float MCReaction::Diff_beam_x_mu_phi_mc() {
//   if (_diff_beam_x_mu_phi_mc != _diff_beam_x_mu_phi_mc) CalcMissMass_mc();
//   return _diff_beam_x_mu_phi_mc;
// }

// float MCReaction::MM2_exclusive_mc() {
//   if (_MM2_exclusive_mc != _MM2_exclusive_mc) CalcMissMass_mc();
//   return _MM2_exclusive_mc;
// }
// float MCReaction::Energy_excl_mc() {
//   if (_excl_Energy_mc != _excl_Energy_mc) CalcMissMass_mc();
//   return _excl_Energy_mc;
// }
// float MCReaction::x_mu_momentum_mc() {
//   if (_rec_x_mu_mom_mc != _rec_x_mu_mom_mc) CalcMissMass_mc();
//   return _rec_x_mu_mom_mc;
// }
// float MCReaction::x_mu_theta_lab_mc() {
//   if (_rec_x_mu_theta_mc != _rec_x_mu_theta_mc) CalcMissMass_mc();
//   return _rec_x_mu_theta_mc;
// }
// float MCReaction::x_mu_Phi_lab_mc() {
//   if (_x_mu_phi_mc != _x_mu_phi_mc) CalcMissMass_mc();
//   return _x_mu_phi_mc;
// }

std::string MCReaction::CsvHeader() {
  return "e_rec_p,e_rec_theta,e_rec_phi,e_sec,e_thrown_p,e_thrown_theta,e_thrown_phi\n";
}
std::string MCReaction::ReacToCsv() {
  // e_rec_p,e_rec_theta,e_rec_phi,e_sec,e_thrown_p,e_thrown_theta,e_thrown_phi
  std::string out = "";
  out += std::to_string(_elec->P()) + ",";
  out += std::to_string(_elec->Theta()) + ",";
  out += std::to_string(_elec->Phi()) + ",";
  out += std::to_string(_sector) + ",";
  out += std::to_string(_elec_mc->P()) + ",";
  out += std::to_string(_elec_mc->Theta()) + ",";
  out += std::to_string(_elec_mc->Phi()) + "\n";

  return out;
}
