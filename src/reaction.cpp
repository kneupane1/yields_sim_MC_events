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
  this->SetElec();

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

  _prot = std::make_unique<TLorentzVector>();
  _pip = std::make_unique<TLorentzVector>();
  _pim = std::make_unique<TLorentzVector>();
  _other = std::make_unique<TLorentzVector>();
  _neutron = std::make_unique<TLorentzVector>();
}

Reaction::~Reaction() {}
auto objMomCorr = std::make_shared<mom_corr>();

void Reaction::SetElec() {
  _hasE = true;
  _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);
  *_gamma += *_beam - *_elec;  // be careful you are commenting this only to include the momentum correction

  // // // // // Can calculate W and Q2 here (useful for simulations as sim do not have elec mom corrections)
  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);
  _elec_mom = _elec->P();
}
// void Reaction::SetMomCorrElec() {
//   // Below shows how the corrections are to be applied using the ROOT momentum 4-vector using the above code:

//   // New electron momentum corrections
//   fe = objMomCorr->dppC(_data->px(0), _data->py(0), _data->pz(0), _data->dc_sec(0), 0) + 1;
//   _mom_corr_elec->SetXYZM(_data->px(0) * fe, _data->py(0) * fe, _data->pz(0) * fe,
//                           MASS_E);  // this is new electron mom corrections aug 2022

//   // _elec_mom_corrected = (dpp(_data->px(0), _data->py(0), _data->pz(0), _data->dc_sec(0), 0) + 1);
//   // _mom_corr_elec->SetXYZM(_data->px(0) * _elec_mom_corrected, _data->py(0) * _elec_mom_corrected,
//   //                         _data->pz(0) * _elec_mom_corrected, MASS_E);

//   // _mom_corr_elec->SetPxPyPzE(_data->px(0) * _elec_mom_corrected, _data->py(0) * _elec_mom_corrected,
//   //                            _data->pz(0) * _elec_mom_corrected, _elec_mom * _elec_mom_corrected);

//   *_gamma += *_beam - *_mom_corr_elec;
//   // _W_after = physics::W_calc(*_beam, *_mom_corr_elec);
//   _W = physics::W_calc(*_beam, *_mom_corr_elec);
//   _Q2 = physics::Q2_calc(*_beam, *_mom_corr_elec);

//   _P_elec = _mom_corr_elec->P();
//   // _E_elec = _mom_corr_elec->E();
//   _theta_e = _mom_corr_elec->Theta() * 180 / PI;
// }
// double Reaction::Corr_elec_mom() {
//   if (_P_elec != _P_elec) SetMomCorrElec();
//   // std::cout << " elec mom corrected " << _elec_mom_corrected << std::endl;

//   return _P_elec;
// }

void Reaction::SetProton(int i) {
  _numProt++;
  _numPos++;
  _hasP = true;
  _prot_status = abs(_data->status(i));
  _Energy_loss_uncorr_prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);

  // _thetaDC_r1_Prot = RAD2DEG * (atan2(sqrt(pow(_data->dc_r1_x(i), 2) + pow(_data->dc_r1_y(i), 2)), _data->dc_r1_z(i)));

  _prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
  _mom_corr_prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);

  // _sectorProt = _data->dc_sec(i);

  // _prot_mom_uncorr = _Energy_loss_uncorr_prot->P();

  // _prot_theta_uncorr = _Energy_loss_uncorr_prot->Theta() * 180 / PI;

  // if (_Energy_loss_uncorr_prot->Phi() > 0)
  //   _prot_phi_uncorr = _Energy_loss_uncorr_prot->Phi() * 180 / PI;
  // else if (_Energy_loss_uncorr_prot->Phi() < 0)
  //   _prot_phi_uncorr = (_Energy_loss_uncorr_prot->Phi() + 2 * PI) * 180 / PI;

  // _is_FD_Prot = objMomCorr->is_FD(_prot_status);
  // _is_CD_Prot = objMomCorr->is_CD(_prot_status);
  // // _is_lower_band = objMomCorr->is_lower_band(_prot_mom_uncorr, _thetaDC_r1_Prot, _prot_status);

  // if (_is_CD_Prot) {
  //   // _prot_mom_tmt = _prot_mom_uncorr;
  //   // _prot_theta_tmt = _prot_theta_uncorr;
  //   // _prot_phi_tmt = _prot_phi_uncorr;

  //   _prot_mom_tmt = objMomCorr->CD_prot_Emom_corr(_prot_mom_uncorr, _prot_theta_uncorr);
  //   // _prot_theta_tmt = objMomCorr->CD_prot_Eth_corr(_prot_mom_uncorr, _prot_theta_uncorr);
  //   // _prot_phi_tmt = objMomCorr->CD_prot_Eph_corr(_prot_mom_uncorr, _prot_theta_uncorr, _prot_phi_uncorr);
  // }
  // if (_is_FD_Prot) {
  //   // // these are Andrey's corrections
  //   if (_prot_theta_uncorr < 27) {
  //     // _prot_theta_tmt = _prot_theta_uncorr;
  //     // _prot_phi_tmt = _prot_phi_uncorr;
  //     _prot_mom_tmt = _prot_mom_uncorr + exp(-2.739 - 3.932 * _prot_theta_uncorr) + 0.002907;
  //   } else {
  //     // _prot_theta_tmt = _prot_theta_uncorr;
  //     // _prot_phi_tmt = _prot_phi_uncorr;
  //     _prot_mom_tmt = _prot_mom_uncorr + exp(-1.2 - 4.228 * _prot_mom_uncorr) + 0.007502;
  //   }

  //   // if (_is_lower_band) {
  //   //   _prot_theta_tmt = objMomCorr->FD_prot_Eth_corr_lower(_prot_mom_uncorr, _prot_theta_uncorr);
  //   //   _prot_phi_tmt = objMomCorr->FD_prot_Eph_corr_lower(_prot_mom_uncorr, _prot_theta_uncorr, _prot_phi_uncorr);
  //   //   if (_prot_mom_uncorr >= 1.0)
  //   //     _prot_mom_tmt = objMomCorr->FD_prot_Emom_corr_lower(_prot_mom_uncorr, _prot_theta_uncorr);
  //   //   else
  //   //     _prot_mom_tmt =
  //   //         _prot_mom_uncorr + objMomCorr->A_p(_prot_mom_uncorr, _prot_theta_uncorr, _thetaDC_r1_Prot, _sectorProt) +
  //   //         objMomCorr->B_p(_prot_mom_uncorr, _prot_theta_uncorr, _thetaDC_r1_Prot, _sectorProt) / _prot_mom_uncorr;

  //   // } else {
  //   //   _prot_theta_tmt = objMomCorr->FD_prot_Eth_corr_upper(_prot_mom_uncorr, _prot_theta_uncorr);
  //   //   _prot_phi_tmt = objMomCorr->FD_prot_Eph_corr_upper(_prot_mom_uncorr, _prot_theta_uncorr, _prot_phi_uncorr);

  //   //   if (_prot_mom_uncorr >= 1.0)
  //   //     _prot_mom_tmt = objMomCorr->FD_prot_Emom_corr_upper(_prot_mom_uncorr, _prot_theta_uncorr);
  //   //   else
  //   //     _prot_mom_tmt =
  //   //         _prot_mom_uncorr + objMomCorr->A_p(_prot_mom_uncorr, _prot_theta_uncorr, _thetaDC_r1_Prot, _sectorProt) +
  //   //         objMomCorr->B_p(_prot_mom_uncorr, _prot_theta_uncorr, _thetaDC_r1_Prot, _sectorProt) / _prot_mom_uncorr;
  //   // }
  // }
  // // _px_prime_prot_E = _data->px(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr)) *
  // //                    (sin(DEG2RAD * _prot_theta_tmt) / sin(DEG2RAD * _prot_theta_uncorr)) *
  // //                    (cos(DEG2RAD * _prot_phi_tmt) / cos(DEG2RAD * _prot_phi_uncorr));
  // // _py_prime_prot_E = _data->py(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr)) *
  // //                    (sin(DEG2RAD * _prot_theta_tmt) / sin(DEG2RAD * _prot_theta_uncorr)) *
  // //                    (sin(DEG2RAD * _prot_phi_tmt) / sin(DEG2RAD * _prot_phi_uncorr));
  // // _pz_prime_prot_E = _data->pz(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr)) *
  // //                    (cos(DEG2RAD * _prot_theta_tmt) / cos(DEG2RAD * _prot_theta_uncorr));

  // _px_prime_prot_E = _data->px(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
  // _py_prime_prot_E = _data->py(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
  // _pz_prime_prot_E = _data->pz(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));

  // /////// _prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P); // energy loss corrected
  // /////// _mom_corr_prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);  // energy loss
  // /// corrected

  // // // Below shows how the corrections are to be applied using the ROOT momentum 4-vector using the above code:
  // if (_is_FD) {
  //   fpro = objMomCorr->dppC(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, _data->dc_sec(i), 3) + 1;
  // } else {
  //   fpro = 1.0;
  // }
  // // // one question here are these corrections good for all FD protons or just for FD prot with FD pip, FD pim???
  // // // // // _px_prime_prot_E = _data->px(i) * fpro * ((_prot_mom_tmt) / (_prot_mom_uncorr));
  // // // // // _py_prime_prot_E = _data->py(i) * fpro * ((_prot_mom_tmt) / (_prot_mom_uncorr));
  // // // // // _pz_prime_prot_E = _data->pz(i) * fpro * ((_prot_mom_tmt) / (_prot_mom_uncorr));
  // // // // // _prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);

  // _prot->SetXYZM(_px_prime_prot_E * fpro, _py_prime_prot_E * fpro, _pz_prime_prot_E * fpro,
  //                MASS_P);  // energy loss + FD had corr

  // // _mom_corr_prot->SetXYZM(_px_prime_prot_E * fpro, _py_prime_prot_E * fpro, _pz_prime_prot_E * fpro, MASS_P);
}
void Reaction::SetPip(int i) {
  _numPip++;
  _numPos++;
  _hasPip = true;
  _pip_status = abs(_data->status(i));
  _sectorPip = _data->dc_sec(i);
  // _thetaDC_r1_Pip = RAD2DEG * (atan2(sqrt(pow(_data->dc_r1_x(i), 2) + pow(_data->dc_r1_y(i), 2)), _data->dc_r1_z(i)));

  _Energy_loss_uncorr_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
  _pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
  _mom_corr_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);

  // _pip_mom_uncorr = _Energy_loss_uncorr_pip->P();
  // _pip_theta_uncorr = _Energy_loss_uncorr_pip->Theta() * 180 / PI;
  // if (_Energy_loss_uncorr_pip->Phi() > 0)
  //   _pip_phi_uncorr = _Energy_loss_uncorr_pip->Phi() * 180 / PI;
  // else if (_Energy_loss_uncorr_pip->Phi() < 0)
  //   _pip_phi_uncorr = (_Energy_loss_uncorr_pip->Phi() + 2 * PI) * 180 / PI;

  // _is_FD = objMomCorr->is_FD(_pip_status);
  // _is_CD = objMomCorr->is_CD(_pip_status);
  // // _is_lower_band = objMomCorr->is_lower_band(_pip_mom_uncorr, _thetaDC_r1_Pip, _pip_status);

  // if (_is_CD) {
  //   // _pip_mom_tmt = _pip_mom_uncorr;
  //   // _pip_theta_tmt = _pip_theta_uncorr;
  //   // _pip_phi_tmt = _pip_phi_uncorr;

  //   _pip_mom_tmt = objMomCorr->CD_pip_Emom_corr(_pip_mom_uncorr, _pip_theta_uncorr);
  //   // _pip_theta_tmt = objMomCorr->CD_pip_Eth_corr(_pip_mom_uncorr, _pip_theta_uncorr);
  //   // _pip_phi_tmt = objMomCorr->CD_pip_Eph_corr(_pip_mom_uncorr, _pip_theta_uncorr, _pip_phi_uncorr);
  // }
  // if (_is_FD) {
  //   _pip_mom_tmt = _pip_mom_uncorr;
  //   //   // _pip_theta_tmt = _pip_theta_uncorr;
  //   //   // _pip_phi_tmt = _pip_phi_uncorr;

  //   //   if (_is_lower_band) {
  //   //     _pip_theta_tmt = objMomCorr->FD_pip_Eth_corr_lower(_pip_mom_uncorr, _pip_theta_uncorr);
  //   //     _pip_phi_tmt = objMomCorr->FD_pip_Eph_corr_lower(_pip_mom_uncorr, _pip_theta_uncorr, _pip_phi_uncorr);
  //   //     _pip_mom_tmt = objMomCorr->FD_pip_Emom_corr_lower(_pip_mom_uncorr, _pip_theta_uncorr);

  //   //   } else {
  //   //     _pip_theta_tmt = objMomCorr->FD_pip_Eth_corr_upper(_pip_mom_uncorr, _pip_theta_uncorr);
  //   //     _pip_phi_tmt = objMomCorr->FD_pip_Eph_corr_upper(_pip_mom_uncorr, _pip_theta_uncorr, _pip_phi_uncorr);
  //   //     _pip_mom_tmt = objMomCorr->FD_pip_Emom_corr_upper(_pip_mom_uncorr, _pip_theta_uncorr);
  //   //   }
  // }
  // // _px_prime_pip_E = _data->px(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr)) *
  // //                   (sin(DEG2RAD * _pip_theta_tmt) / sin(DEG2RAD * _pip_theta_uncorr)) *
  // //                   (cos(DEG2RAD * _pip_phi_tmt) / cos(DEG2RAD * _pip_phi_uncorr));
  // // _py_prime_pip_E = _data->py(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr)) *
  // //                   (sin(DEG2RAD * _pip_theta_tmt) / sin(DEG2RAD * _pip_theta_uncorr)) *
  // //                   (sin(DEG2RAD * _pip_phi_tmt) / sin(DEG2RAD * _pip_phi_uncorr));
  // // _pz_prime_pip_E = _data->pz(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr)) *
  // //                   (cos(DEG2RAD * _pip_theta_tmt) / cos(DEG2RAD * _pip_theta_uncorr));

  // _px_prime_pip_E = _data->px(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
  // _py_prime_pip_E = _data->py(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
  // _pz_prime_pip_E = _data->pz(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));

  // // _pip->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);
  // // _mom_corr_pip->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);

  // if (_is_FD) {
  //   // _sectorPip = _data->dc_sec(i);
  //   fpip = objMomCorr->dppC(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, _data->dc_sec(i), 1) + 1;
  //   // fpip = objMomCorr->dppC(_data->px(i), _data->py(i), _data->pz(i), _data->dc_sec(i), 1) + 1;

  // } else {
  //   fpip = 1.0;
  // }
  // _pip->SetXYZM(_px_prime_pip_E * fpip, _py_prime_pip_E * fpip, _pz_prime_pip_E * fpip, MASS_PIP);
  // // _mom_corr_pip->SetXYZM(_px_prime_pip_E * fpip, _py_prime_pip_E * fpip, _pz_prime_pip_E * fpip, MASS_PIP);

  // // // _pip->SetXYZM(_data->px(i) * fpip, _data->py(i) * fpip, _data->pz(i) * fpip, MASS_PIP);
  // // // // _mom_corr_pip->SetXYZM(_data->px(i) * fpip, _data->py(i) * fpip, _data->pz(i) * fpip, MASS_PIP);
}

void Reaction::SetPim(int i) {
  _numPim++;
  _numNeg++;
  _hasPim = true;
  _pim_status = abs(_data->status(i));
  _sectorPim = _data->dc_sec(i);
  _thetaDC_r1_Pim = RAD2DEG * (atan2(sqrt(pow(_data->dc_r1_x(i), 2) + pow(_data->dc_r1_y(i), 2)), _data->dc_r1_z(i)));

  _Energy_loss_uncorr_pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
  _pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
  _mom_corr_pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);

  // _pim_mom_uncorr = _Energy_loss_uncorr_pim->P();
  // _pim_theta_uncorr = _Energy_loss_uncorr_pim->Theta() * 180 / PI;
  // if (_Energy_loss_uncorr_pim->Phi() > 0)
  //   _pim_phi_uncorr = _Energy_loss_uncorr_pim->Phi() * 180 / PI;
  // else if (_Energy_loss_uncorr_pim->Phi() < 0)
  //   _pim_phi_uncorr = (_Energy_loss_uncorr_pim->Phi() + 2 * PI) * 180 / PI;

  // _is_FD = objMomCorr->is_FD(_pim_status);
  // _is_CD = objMomCorr->is_CD(_pim_status);
  // // _is_lower_band = objMomCorr->is_lower_band(_pim_mom_uncorr, _thetaDC_r1_Pim, _pim_status);

  // if (_is_CD) {
  //   // _pim_mom_tmt = _pim_mom_uncorr;
  //   // _pim_theta_tmt = _pim_theta_uncorr;
  //   // _pim_phi_tmt = _pim_phi_uncorr;

  //   _pim_mom_tmt = objMomCorr->CD_pim_Emom_corr(_pim_mom_uncorr, _pim_theta_uncorr);
  //   // _pim_theta_tmt = objMomCorr->CD_pim_Eth_corr(_pim_mom_uncorr, _pim_theta_uncorr);
  //   // _pim_phi_tmt = objMomCorr->CD_pim_Eph_corr(_pim_mom_uncorr, _pim_theta_uncorr, _pim_phi_uncorr);
  // }
  // if (_is_FD) {
  //   _pim_mom_tmt = _pim_mom_uncorr;
  //   // // _pim_theta_tmt = _pim_theta_uncorr;
  //   // // _pim_phi_tmt = _pim_phi_uncorr;

  //   //   if (_is_lower_band) {
  //   //     _pim_theta_tmt = objMomCorr->FD_pim_Eth_corr_lower(_pim_mom_uncorr, _pim_theta_uncorr);
  //   //     _pim_phi_tmt = objMomCorr->FD_pim_Eph_corr_lower(_pim_mom_uncorr, _pim_theta_uncorr, _pim_phi_uncorr);
  //   //     _pim_mom_tmt = objMomCorr->FD_pim_Emom_corr_lower(_pim_mom_uncorr, _pim_theta_uncorr);

  //   //   } else {
  //   //     _pim_theta_tmt = objMomCorr->FD_pim_Eth_corr_upper(_pim_mom_uncorr, _pim_theta_uncorr);
  //   //     _pim_phi_tmt = objMomCorr->FD_pim_Eph_corr_upper(_pim_mom_uncorr, _pim_theta_uncorr, _pim_phi_uncorr);
  //   //     _pim_mom_tmt = objMomCorr->FD_pim_Emom_corr_upper(_pim_mom_uncorr, _pim_theta_uncorr);
  //   //   }
  // }
  // // _px_prime_pim_E = _data->px(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr)) *
  // //                   (sin(DEG2RAD * _pim_theta_tmt) / sin(DEG2RAD * _pim_theta_uncorr)) *
  // //                   (cos(DEG2RAD * _pim_phi_tmt) / cos(DEG2RAD * _pim_phi_uncorr));
  // // _py_prime_pim_E = _data->py(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr)) *
  // //                   (sin(DEG2RAD * _pim_theta_tmt) / sin(DEG2RAD * _pim_theta_uncorr)) *
  // //                   (sin(DEG2RAD * _pim_phi_tmt) / sin(DEG2RAD * _pim_phi_uncorr));
  // // _pz_prime_pim_E = _data->pz(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr)) *
  // //                   (cos(DEG2RAD * _pim_theta_tmt) / cos(DEG2RAD * _pim_theta_uncorr));

  // // _pim->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, MASS_PIM);
  // // _mom_corr_pim->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, MASS_PIM);

  // _px_prime_pim_E = _data->px(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
  // _py_prime_pim_E = _data->py(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
  // _pz_prime_pim_E = _data->pz(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));

  // if (_is_FD) {
  //   fpim = objMomCorr->dppC(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, _data->dc_sec(i), 2) + 1;
  //   // fpim = objMomCorr->dppC(_data->px(i), _data->py(i), _data->pz(i), _data->dc_sec(i), 2) + 1;
  // } else {
  //   fpim = 1.0;
  // }
  // _pim->SetXYZM(_px_prime_pim_E * fpim, _py_prime_pim_E * fpim, _pz_prime_pim_E * fpim, MASS_PIM);
  // // _mom_corr_pim->SetXYZM(_px_prime_pim_E * fpim, _py_prime_pim_E * fpim, _pz_prime_pim_E * fpim, MASS_PIM);

  // // // _pim->SetXYZM(_data->px(i) * fpim, _data->py(i) * fpim, _data->pz(i) * fpim, MASS_PIM);
  // // // _mom_corr_pim->SetXYZM(_data->px(i) * fpim, _data->py(i) * fpim, _data->pz(i) * fpim, MASS_PIM);
}

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
// // // //// Now Our version of Momentum corrections based on aug task force but now fd in not separated in more cases

// void Reaction::Prot_HMom_corr(int status_prot, int sector_Prot, float alPFD, float alPCD[3]) {
//   auto uncorr_prot = std::make_unique<TLorentzVector>();

//   *uncorr_prot += (*_prot);
//   _is_FD_Prot = objMomCorr->is_FD(status_prot);
//   _is_CD_Prot = objMomCorr->is_CD(status_prot);

//   _prot_mom = uncorr_prot->P();

//   if (uncorr_prot->Phi() > 0)
//     _prot_phi = uncorr_prot->Phi() * 180 / PI;
//   else if (_prot->Phi() < 0)
//     _prot_phi = (uncorr_prot->Phi() + 2 * PI) * 180 / PI;

//   if (_is_CD_Prot) {
//     _prot_mom_prime = objMomCorr->CD_prot_Hmom_corr(_prot_mom, _prot_phi, alPCD);
//   }
//   if (_is_FD_Prot) {
//     _prot_mom_prime = objMomCorr->FD_prot_Hmom_corr(_prot_mom, sector_Prot, alPFD);
//   }

//   _px_prime_prot_mom = uncorr_prot->Px() * ((_prot_mom_prime) / (_prot_mom));
//   _py_prime_prot_mom = uncorr_prot->Py() * ((_prot_mom_prime) / (_prot_mom));
//   _pz_prime_prot_mom = uncorr_prot->Pz() * ((_prot_mom_prime) / (_prot_mom));
//   _mom_corr_prot->SetXYZM(_px_prime_prot_mom, _py_prime_prot_mom, _pz_prime_prot_mom, MASS_P);
// }

// void Reaction::Pip_HMom_corr(int status_pip, int sector_Pip, float alPipFD, float alPipCD[3]) {
//   auto uncorr_pip = std::make_unique<TLorentzVector>();
//   *uncorr_pip += (*_pip);
//   _is_FD_Pip = objMomCorr->is_FD(status_pip);
//   _is_CD_Pip = objMomCorr->is_CD(status_pip);

//   _pip_mom = uncorr_pip->P();

//   if (uncorr_pip->Phi() > 0)
//     _pip_phi = uncorr_pip->Phi() * 180 / PI;
//   else if (_pip->Phi() < 0)
//     _pip_phi = (uncorr_pip->Phi() + 2 * PI) * 180 / PI;

//   if (_is_CD_Pip) {
//     _pip_mom_prime = objMomCorr->CD_pip_Hmom_corr(_pip_mom, _pip_phi, alPipCD);
//   }
//   if (_is_FD_Pip) {
//     _pip_mom_prime = objMomCorr->FD_pip_Hmom_corr(_pip_mom, sector_Pip, alPipFD);
//   }

//   _px_prime_pip_mom = uncorr_pip->Px() * ((_pip_mom_prime) / (_pip_mom));
//   _py_prime_pip_mom = uncorr_pip->Py() * ((_pip_mom_prime) / (_pip_mom));
//   _pz_prime_pip_mom = uncorr_pip->Pz() * ((_pip_mom_prime) / (_pip_mom));
//   _mom_corr_pip->SetXYZM(_px_prime_pip_mom, _py_prime_pip_mom, _pz_prime_pip_mom, MASS_PIP);
// }

// void Reaction::Pim_HMom_corr(int status_pim, int sector_Pim, float alPimFD, float alPimCD[3]) {
//   auto uncorr_pim = std::make_unique<TLorentzVector>();
//   *uncorr_pim += (*_pim);
//   _is_FD_Pim = objMomCorr->is_FD(status_pim);
//   _is_CD_Pim = objMomCorr->is_CD(status_pim);

//   _pim_mom = uncorr_pim->P();

//   if (uncorr_pim->Phi() > 0)
//     _pim_phi = uncorr_pim->Phi() * 180 / PI;
//   else if (_pim->Phi() < 0)
//     _pim_phi = (uncorr_pim->Phi() + 2 * PI) * 180 / PI;

//   if (_is_CD_Pim) {
//     _pim_mom_prime = objMomCorr->CD_pim_Hmom_corr(_pim_mom, _pim_phi, alPimCD);
//   }
//   if (_is_FD_Pim) {
//     _pim_mom_prime = objMomCorr->FD_pim_Hmom_corr(_pim_mom, sector_Pim, alPimFD);
//   }

//   _px_prime_pim_mom = uncorr_pim->Px() * ((_pim_mom_prime) / (_pim_mom));
//   _py_prime_pim_mom = uncorr_pim->Py() * ((_pim_mom_prime) / (_pim_mom));
//   _pz_prime_pim_mom = uncorr_pim->Pz() * ((_pim_mom_prime) / (_pim_mom));
//   _mom_corr_pim->SetXYZM(_px_prime_pim_mom, _py_prime_pim_mom, _pz_prime_pim_mom, MASS_PIM);
// }

// // // // //// Now Our version of Momentum corrections based on Aug task fc mom corr

void Reaction::Prot_HMom_corr(int status_prot, int status_pip, int status_pim, int sector_Prot, float alPFD[4], float
alPCD[3]) {
  auto uncorr_prot = std::make_unique<TLorentzVector>();

  *uncorr_prot += (*_prot);
  _is_FD_Prot = objMomCorr->is_FD(status_prot);
  _is_CD_Prot = objMomCorr->is_CD(status_prot);
  _is_FD_Pip = objMomCorr->is_FD(status_pip);
  _is_FD_Pim = objMomCorr->is_FD(status_pim);

  _prot_mom = uncorr_prot->P();
  _prot_theta = uncorr_prot->Theta() * 180 / PI;

  if (uncorr_prot->Phi() > 0)
    _prot_phi = uncorr_prot->Phi() * 180 / PI;
  else if (_prot->Phi() < 0)
    _prot_phi = (uncorr_prot->Phi() + 2 * PI) * 180 / PI;

  if (_is_CD_Prot) {
    _prot_mom_prime = objMomCorr->CD_prot_Hmom_corr(_prot_mom, _prot_phi, alPCD);
  }
  if (_is_FD_Prot) {
    if (_prot_theta < 27) {
      if ((_is_FD_Pip) && (_is_FD_Pim)) {
        _prot_mom_prime = objMomCorr->FD_prot_Hmom_corr_lower_All_FD(_prot_mom, sector_Prot, alPFD[0]);
      } else {
        _prot_mom_prime = objMomCorr->FD_prot_Hmom_corr_lower_Except_All_FD(_prot_mom, sector_Prot, alPFD[1]);
      }
    } else {
      if ((_is_FD_Pip) && (_is_FD_Pim)) {
        _prot_mom_prime = objMomCorr->FD_prot_Hmom_corr_upper_All_FD(_prot_mom, sector_Prot, alPFD[2]);
      } else {
        _prot_mom_prime = objMomCorr->FD_prot_Hmom_corr_upper_Except_All_FD(_prot_mom, sector_Prot, alPFD[3]);
      }
    }
  }

  _px_prime_prot_mom = uncorr_prot->Px() * ((_prot_mom_prime) / (_prot_mom));
  _py_prime_prot_mom = uncorr_prot->Py() * ((_prot_mom_prime) / (_prot_mom));
  _pz_prime_prot_mom = uncorr_prot->Pz() * ((_prot_mom_prime) / (_prot_mom));
  _mom_corr_prot->SetXYZM(_px_prime_prot_mom, _py_prime_prot_mom, _pz_prime_prot_mom, MASS_P);
}

void Reaction::Pip_HMom_corr(int status_prot, int status_pip, int status_pim, int sector_Pip, float alPipFD[4],
                             float alPipCD[3]) {
  auto uncorr_pip = std::make_unique<TLorentzVector>();
  *uncorr_pip += (*_pip);
  _is_FD_Prot = objMomCorr->is_FD(status_prot);
  _is_FD_Pip = objMomCorr->is_FD(status_pip);
  _is_CD_Pip = objMomCorr->is_CD(status_pip);
  _is_FD_Pim = objMomCorr->is_FD(status_pim);

  _pip_mom = uncorr_pip->P();
  _pip_theta = uncorr_pip->Theta() * 180 / PI;

  if (uncorr_pip->Phi() > 0)
    _pip_phi = uncorr_pip->Phi() * 180 / PI;
  else if (_pip->Phi() < 0)
    _pip_phi = (uncorr_pip->Phi() + 2 * PI) * 180 / PI;

  if (_is_CD_Pip) {
    _pip_mom_prime = objMomCorr->CD_pip_Hmom_corr(_pip_mom, _pip_phi, alPipCD);
  }
  if (_is_FD_Pip) {
    if (_pip_theta < 27) {
      if ((_is_FD_Prot) && (_is_FD_Pim)) {
        _pip_mom_prime = objMomCorr->FD_pip_Hmom_corr_lower_All_FD(_pip_mom, sector_Pip, alPipFD[0]);
      } else {
        _pip_mom_prime = objMomCorr->FD_pip_Hmom_corr_lower_Except_All_FD(_pip_mom, sector_Pip, alPipFD[1]);
      }
    } else {
      if ((_is_FD_Prot) && (_is_FD_Pim)) {
        _pip_mom_prime = objMomCorr->FD_pip_Hmom_corr_upper_All_FD(_pip_mom, sector_Pip, alPipFD[2]);
      } else {
        _pip_mom_prime = objMomCorr->FD_pip_Hmom_corr_upper_Except_All_FD(_pip_mom, sector_Pip, alPipFD[3]);
      }
    }
  }

  _px_prime_pip_mom = uncorr_pip->Px() * ((_pip_mom_prime) / (_pip_mom));
  _py_prime_pip_mom = uncorr_pip->Py() * ((_pip_mom_prime) / (_pip_mom));
  _pz_prime_pip_mom = uncorr_pip->Pz() * ((_pip_mom_prime) / (_pip_mom));
  _mom_corr_pip->SetXYZM(_px_prime_pip_mom, _py_prime_pip_mom, _pz_prime_pip_mom, MASS_PIP);
}

void Reaction::Pim_HMom_corr(int status_prot, int status_pip, int status_pim, int sector_Pim, float alPimFD[4],
                             float alPimCD[3]) {
  auto uncorr_pim = std::make_unique<TLorentzVector>();
  *uncorr_pim += (*_pim);
  _is_FD_Prot = objMomCorr->is_FD(status_prot);
  _is_FD_Pip = objMomCorr->is_FD(status_pip);
  _is_FD_Pim = objMomCorr->is_FD(status_pim);
  _is_CD_Pim = objMomCorr->is_CD(status_pim);

  _pim_mom = uncorr_pim->P();
  _pim_theta = uncorr_pim->Theta() * 180 / PI;

  if (uncorr_pim->Phi() > 0)
    _pim_phi = uncorr_pim->Phi() * 180 / PI;
  else if (_pim->Phi() < 0)
    _pim_phi = (uncorr_pim->Phi() + 2 * PI) * 180 / PI;

  if (_is_CD_Pim) {
    _pim_mom_prime = objMomCorr->CD_pim_Hmom_corr(_pim_mom, _pim_phi, alPimCD);
  }
  if (_is_FD_Pim) {
    if (_pim_theta < 27) {
      if ((_is_FD_Pip) && (_is_FD_Prot)) {
        _pim_mom_prime = objMomCorr->FD_pim_Hmom_corr_lower_All_FD(_pim_mom, sector_Pim, alPimFD[0]);
      } else {
        _pim_mom_prime = objMomCorr->FD_pim_Hmom_corr_lower_Except_All_FD(_pim_mom, sector_Pim, alPimFD[1]);
      }
    } else {
      if ((_is_FD_Pip) && (_is_FD_Prot)) {
        _pim_mom_prime = objMomCorr->FD_pim_Hmom_corr_upper_All_FD(_pim_mom, sector_Pim, alPimFD[2]);
      } else {
        _pim_mom_prime = objMomCorr->FD_pim_Hmom_corr_upper_Except_All_FD(_pim_mom, sector_Pim, alPimFD[3]);
      }
    }
  }

  _px_prime_pim_mom = uncorr_pim->Px() * ((_pim_mom_prime) / (_pim_mom));
  _py_prime_pim_mom = uncorr_pim->Py() * ((_pim_mom_prime) / (_pim_mom));
  _pz_prime_pim_mom = uncorr_pim->Pz() * ((_pim_mom_prime) / (_pim_mom));
  _mom_corr_pim->SetXYZM(_px_prime_pim_mom, _py_prime_pim_mom, _pz_prime_pim_mom, MASS_PIM);
}

void Reaction::CalcMissMass() {
  auto mm_mpim = std::make_unique<TLorentzVector>();
  auto mm_mpip = std::make_unique<TLorentzVector>();
  auto mm_mprot = std::make_unique<TLorentzVector>();
  auto mm_excl = std::make_unique<TLorentzVector>();
  auto mm_excl_corr = std::make_unique<TLorentzVector>();

  *mm_mpim += (*_gamma + *_target);

  // if (TwoPion_missingPim()) {
  //   *mm -= *_prot;
  //   *mm -= *_pip;
  //   // *mm -= *_pim;
  //   _MM = mm->M();
  //   _MM2 = mm->M2();

  // //   // _rec_pim_mom = mm->P();
  // //   // _rec_pim_theta = mm->Theta() * 180 / PI;

  // //   // if (mm->Phi() >= 0)
  // //   //   _rec_pim_phi = (mm->Phi() * 180 / PI);
  // //   // else if (mm->Phi() < 0)
  // //   //   _rec_pim_phi = ((mm->Phi() + 2 * PI) * 180 / PI);

  // // //   // // // _x_mu_E = mm->E();
  // // //   // // // _x_mu_P = mm->P();
  // // //   // // // _x_mu_Px = mm->Px();
  // // //   // // // _x_mu_Py = mm->Py();
  // // //   // // // _x_mu_Pz = mm->Pz();
  // // //   // // // _x_mu_theta = mm->Theta() * RAD2DEG;
  // // //   // // // _x_mu_m2 = mm->E() * mm->E() - mm->P() * mm->P();
  // // //   // // // _x_mu_m = mm->E() - mm->P();
  // // //   // // //   //
  // }
  if (TwoPion_exclusive()) {
    // *mm -= *_mom_corr_prot;
    // *mm -= *_mom_corr_pip;
    // // *mm -= *_pim;
    // _MM = mm->M();
    // _MM2 = mm->M2();

    // *mm_excl += (*_gamma + *_target);
    // *mm_excl -= *_mom_corr_prot;
    // *mm_excl -= *_mom_corr_pip;
    // *mm_excl -= *_mom_corr_pim;

    *mm_mpim -= *_prot;
    *mm_mpim -= *_pip;
    _MM_mPim = mm_mpim->M();
    _MM2_mPim = mm_mpim->M2();

    *mm_excl += (*_gamma + *_target);
    *mm_excl -= *_prot;
    *mm_excl -= *_pip;
    *mm_excl -= *_pim;

    _MM2_exclusive = mm_excl->M2();
    _excl_Energy = mm_excl->E();
    _mom_exclusive = mm_excl->P();

    *mm_excl_corr += (*_gamma + *_target);
    *mm_excl_corr -= *_mom_corr_prot;
    *mm_excl_corr -= *_mom_corr_pip;
    *mm_excl_corr -= *_mom_corr_pim;

    _MM2_exclusive_corr = mm_excl_corr->M2();
    _excl_Energy_corr = mm_excl_corr->E();
    // _mom_exclusive_corr = mm_excl_corr->P();

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

    //   // //   // for mPip peak with exclusive events
    //   *mm_mpip += (*_gamma + *_target);
    //   *mm_mpip -= *_mom_corr_prot;
    //   *mm_mpip -= *_mom_corr_pim;
    //   _MM2_mPip = mm_mpip->M2();

    //   // //   // for mProt peak with exclusive events
    //   *mm_mprot += (*_gamma + *_target);
    //   *mm_mprot -= *_mom_corr_pip;
    //   *mm_mprot -= *_mom_corr_pim;
    //   _MM2_mProt = mm_mprot->M2();
    // }
    // // if (TwoPion_missingPip()) {
    *mm_mpip += (*_gamma + *_target);
    *mm_mpip -= *_prot;
    *mm_mpip -= *_pim;
    _MM2_mPip = mm_mpip->M2();
    // // }
    // // if (TwoPion_missingProt()) {
    *mm_mprot += (*_gamma + *_target);
    *mm_mprot -= *_pip;
    *mm_mprot -= *_pim;
    _MM2_mProt = mm_mprot->M2();
  }
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

float Reaction::MM_mPim() {
  if (_MM_mPim != _MM_mPim) CalcMissMass();
  return _MM_mPim;
}
float Reaction::MM2_mPim() {
  if (_MM2_mPim != _MM2_mPim) CalcMissMass();
  return _MM2_mPim;
}
float Reaction::MM2_exclusive() {
  if (_MM2_exclusive != _MM2_exclusive) CalcMissMass();
  return _MM2_exclusive;
}
float Reaction::MM2_exclusive_corr() {
  if (_MM2_exclusive_corr != _MM2_exclusive_corr) CalcMissMass();
  return _MM2_exclusive_corr;
}
float Reaction::MM2_mPip() {
  if (_MM2_mPip != _MM2_mPip) CalcMissMass();
  return _MM2_mPip;
}
float Reaction::MM2_mProt() {
  if (_MM2_mProt != _MM2_mProt) CalcMissMass();
  return _MM2_mProt;
}

float Reaction::MM2_mPim_corr() {
  // if (_rec_pim_mom != _rec_pim_mom) CalcMissMass();

  // if (TwoPion_missingPim()) {
  if (TwoPion_exclusive()) {
    auto missingpim_ = std::make_unique<TLorentzVector>();
    // *missingpim_ += *_gamma + *_target - *_prot - *_pip;
    *missingpim_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pip;

    return missingpim_->M2();
    // return _rec_pim_mom;

  } else
    return NAN;
}

float Reaction::MM2_mPip_corr() {
  // if (TwoPion_missingPip()) {
  if (TwoPion_exclusive()) {
    auto missingpip_ = std::make_unique<TLorentzVector>();
    // *missingpip_ += *_gamma + *_target - *_prot - *_pim;
    *missingpip_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pim;

    return missingpip_->M2();
  } else
    return NAN;
}

float Reaction::MM2_mProt_corr() {
  // if (TwoPion_missingProt()) {
  if (TwoPion_exclusive()) {
    auto missingprot_ = std::make_unique<TLorentzVector>();
    // *missingprot_ += *_gamma + *_target - *_pip - *_pim;
    *missingprot_ += *_gamma + *_target - *_mom_corr_pip - *_mom_corr_pim;

    return missingprot_->M2();
  } else
    return NAN;
}

float Reaction::Energy_excl() {
  if (_excl_Energy != _excl_Energy) CalcMissMass();
  return _excl_Energy;
  // else
  // return NAN;
}

float Reaction::Mom_excl() {
  if (_mom_exclusive != _mom_exclusive) CalcMissMass();
  return _mom_exclusive;
  // else
  // return NAN;
}

float Reaction::Energy_excl_corr() {
  if (_excl_Energy_corr != _excl_Energy_corr) CalcMissMass();
  return _excl_Energy_corr;
  // else
  // return NAN;
}

// float Reaction::Mom_excl_corr() {
//   if (_mom_exclusive_corr != _mom_exclusive_corr) CalcMissMass();
//   return _mom_exclusive_corr;
//   // else
//   // return NAN;
// }
float Reaction::pim_momentum() {
  // if (_rec_pim_mom != _rec_pim_mom) CalcMissMass();

  // if (TwoPion_missingPim()) {
  if (TwoPion_exclusive()) {
    auto missingpim_ = std::make_unique<TLorentzVector>();
    // *missingpim_ += *_gamma + *_target - *_prot - *_pip;
    *missingpim_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pip;
    /// 4-vect approach
    // *missingpim_ += *_pim - (*_gamma + *_target - *_prot - *_pip);
    return missingpim_->P();
    // return _rec_pim_mom;

  } else
    return NAN;
}
float Reaction::pim_theta_lab() {
  // if (_rec_pim_theta != _rec_pim_theta) CalcMissMass();

  // if (TwoPion_missingPim()) {
  if (TwoPion_exclusive()) {
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

  // if (TwoPion_missingPim()) {
  if (TwoPion_exclusive()) {
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
    if (_pim->Phi() > 0) {
      // std::cout << "phi root >0 is " << _pim->Phi() * 180 / PI << std::endl;
      return _pim->Phi() * 180 / PI;
    } else if (_pim->Phi() < 0) {
      // std::cout << "phi root < 0 is " << (_pim->Phi() + 2 * PI) * 180 / PI << std::endl;
      return (_pim->Phi() + 2 * PI) * 180 / PI;
    } else
      return NAN;
  } else
    return NAN;
}
float Reaction::pim_momentum_corrected() {
  if (TwoPion_exclusive())
    return _mom_corr_pim->P();
  else
    return NAN;
}
float Reaction::w_hadron() {
  if (TwoPion_exclusive())
    return ((*_prot) + (*_pip) + (*_pim)).Mag();
  else
    return NAN;
}
// float Reaction::w_difference() {
//   if (TwoPion_exclusive())
//     return (physics::W_calc(*_beam, *_mom_corr_elec) - ((*_prot) + (*_pip) + (*_pim)).Mag());
//   else
//     return NAN;
// }

// float Reaction::w_hadron_corr() {
//   if (TwoPion_exclusive())
//     return ((*_mom_corr_prot) + (*_mom_corr_pip) + (*_mom_corr_pim)).Mag();
//   else
//     return NAN;
// }
// float Reaction::w_difference_corr() {
//   if (TwoPion_exclusive())
//     return (physics::W_calc(*_beam, *_mom_corr_elec) -
//             ((*_mom_corr_prot) + (*_mom_corr_pip) + (*_mom_corr_pim)).Mag());
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
    // *missingpip_ += *_gamma + *_target - *_prot - *_pim;
    *missingpip_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pim;
    // 4-vect approach
    // *missingpip_ += *_pip - (*_gamma + *_target - *_prot - *_pim);

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

float Reaction::pip_momentum_corrected() {
  if (TwoPion_exclusive())
    return _mom_corr_pip->P();
  else
    return NAN;
}
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
    // *missingprot_ += *_gamma + *_target - *_pip - *_pim;
    *missingprot_ += *_gamma + *_target - *_mom_corr_pip - *_mom_corr_pim;
    // 4-vect approach
    // *missingprot_ += *_prot -(*_gamma + *_target - *_pip - *_pim);

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

float Reaction::prot_momentum_corrected() {
  if (TwoPion_exclusive())
    return _mom_corr_prot->P();
  else
    return NAN;
}
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

void Reaction::invMassPpim() {
  auto inv_Ppim = std::make_unique<TLorentzVector>();
  *inv_Ppim += *_mom_corr_prot;
  *inv_Ppim += *_mom_corr_pim;
  if (TwoPion_exclusive()) _inv_Ppim = inv_Ppim->M();
}
void Reaction::invMasspippim() {
  auto inv_pip_pim = std::make_unique<TLorentzVector>();
  *inv_pip_pim += *_mom_corr_pip;
  *inv_pip_pim += *_mom_corr_pim;
  if (TwoPion_exclusive()) _inv_pip_pim = inv_pip_pim->M();
}
void Reaction::invMassPpip() {
  auto inv_Ppip = std::make_unique<TLorentzVector>();
  *inv_Ppip += *_mom_corr_prot;
  *inv_Ppip += *_mom_corr_pip;
  if (TwoPion_exclusive()) _inv_Ppip = inv_Ppip->M();
}
float Reaction::inv_Ppip() {
  if (_inv_Ppip != _inv_Ppip) invMassPpip();
  return _inv_Ppip;
}
float Reaction::inv_Ppim() {
  if (_inv_Ppim != _inv_Ppim) invMassPpim();
  return _inv_Ppim;
}
float Reaction::inv_Pippim() {
  if (_inv_pip_pim != _inv_pip_pim) invMasspippim();
  return _inv_pip_pim;
}

// Function to rotate a point around Z and Y
TVector3 Reaction::getRotTiltPoint(TVector3& point, int sec) {
  TRotation rotationZ;
  TRotation rotationY;
  rotationZ.RotateZ(-TMath::Pi() / 3.0 * (sec - 1));
  rotationY.RotateY(-TMath::Pi() / 180.0 * 25);

  TVector3 rotatedPoint = rotationY * (rotationZ * point);
  return rotatedPoint;
}
// Function to get the residual as a TVector3
void Reaction::calculateResidualpcal() {
  int sec = _data->dc_sec(0);
  Double_t hx = _data->ec_pcal_hx(0);
  Double_t hy = _data->ec_pcal_hy(0);
  Double_t hz = _data->ec_pcal_hz(0);

  TVector3 xyz(hx - _data->ec_pcal_x(0), hy - _data->ec_pcal_y(0), hz - _data->ec_pcal_z(0));
  TVector3 residual = getRotTiltPoint(xyz, sec);

  _residualXpcal = residual.X();
  _residualYpcal = residual.Y();
  _residualZpcal = residual.Z();

  TVector3 xyz_(_data->ec_pcal_x(0), _data->ec_pcal_y(0), _data->ec_pcal_z(0));
  TVector3 xyz_rot = getRotTiltPoint(xyz_, sec);
  _Xpcal_rot = xyz_rot.X();
  _Ypcal_rot = xyz_rot.Y();
  // return getRotTiltPoint(xyz, sec);
}

// Function to get the x value of the residual
Double_t Reaction::getResidualXpcal() {
  if (std::isnan(_residualXpcal)) calculateResidualpcal();
  return _residualXpcal;
}
Double_t Reaction::getResidualYpcal() {
  if (std::isnan(_residualYpcal)) calculateResidualpcal();
  return _residualYpcal;
}
Double_t Reaction::getResidualZpcal() {
  if (std::isnan(_residualZpcal)) calculateResidualpcal();
  return _residualZpcal;
}
// Function to get the x value of the residual
Double_t Reaction::Xpcal() { return _data->ec_pcal_x(0); }
Double_t Reaction::Ypcal() { return _data->ec_pcal_y(0); }

Double_t Reaction::Xpcal_rot() {
  if (std::isnan(_Xpcal_rot)) calculateResidualpcal();
  return _Xpcal_rot;
}
Double_t Reaction::Ypcal_rot() {
  if (std::isnan(_Ypcal_rot)) calculateResidualpcal();
  return _Ypcal_rot;
}

// Function to get the residual as a TVector3
void Reaction::calculateResidualecin() {
  int sec = _data->dc_sec(0);
  Double_t hx = _data->ec_ecin_hx(0);
  Double_t hy = _data->ec_ecin_hy(0);
  Double_t hz = _data->ec_ecin_hz(0);

  TVector3 xyz_ecin(hx - _data->ec_ecin_x(0), hy - _data->ec_ecin_y(0), hz - _data->ec_ecin_z(0));
  TVector3 residual_ecin = getRotTiltPoint(xyz_ecin, sec);

  _residualXecin = residual_ecin.X();
  _residualYecin = residual_ecin.Y();
  _residualZecin = residual_ecin.Z();

  TVector3 xyz_ecin_(_data->ec_ecin_x(0), _data->ec_ecin_y(0), _data->ec_ecin_z(0));
  TVector3 xyz_rot_ecin = getRotTiltPoint(xyz_ecin_, sec);

  _Xecin_rot = xyz_rot_ecin.X();
  _Yecin_rot = xyz_rot_ecin.Y();
  // return getRotTiltPoint(xyz, sec);
}

// Function to get the x value of the residual
Double_t Reaction::getResidualXecin() {
  if (std::isnan(_residualXecin)) calculateResidualecin();
  return _residualXecin;
}
Double_t Reaction::getResidualYecin() {
  if (std::isnan(_residualYecin)) calculateResidualecin();
  return _residualYecin;
}
Double_t Reaction::getResidualZecin() {
  if (std::isnan(_residualZecin)) calculateResidualecin();
  return _residualZecin;
}

// Function to get the x value of the residual
Double_t Reaction::Xecin() { return _data->ec_ecin_x(0); }
Double_t Reaction::Yecin() { return _data->ec_ecin_y(0); }

Double_t Reaction::Xecin_rot() {
  if (std::isnan(_Xecin_rot)) calculateResidualecin();
  return _Xecin_rot;
}
Double_t Reaction::Yecin_rot() {
  if (std::isnan(_Yecin_rot)) calculateResidualecin();
  return _Yecin_rot;
}

float Reaction::sampling_fraction() { return _data->ec_tot_energy(0) / _data->p(0); }

//////////////////////////
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
  // // original one withiout going through momentum corrections
  // _boosted_prot = std::make_unique<TLorentzVector>(*_prot);
  // _boosted_pip = std::make_unique<TLorentzVector>(*_pip);
  // _boosted_pim = std::make_unique<TLorentzVector>(*_gamma + *_target - *_prot - *_pip);  //(*_pim);
  // _boosted_gamma = std::make_unique<TLorentzVector>(*_gamma);
  // _boosted_pim_measured = std::make_unique<TLorentzVector>(*_pim);

  // _rotated_prot = std::make_unique<TLorentzVector>(*_prot);
  // _rotated_pip = std::make_unique<TLorentzVector>(*_pip);
  // _rotated_pim = std::make_unique<TLorentzVector>(*_gamma + *_target - *_prot - *_pip);  //(*_pim);
  // _rotated_pim_measured = std::make_unique<TLorentzVector>(*_pim);

  // // new one which are already gone through momentum corrections

  _boosted_prot = std::make_unique<TLorentzVector>(*_mom_corr_prot);
  _boosted_pip = std::make_unique<TLorentzVector>(*_mom_corr_pip);
  _boosted_pim = std::make_unique<TLorentzVector>(*_gamma + *_target - *_mom_corr_prot - *_mom_corr_pip);  //(*_pim); // careful here because it is not in exclusive set up sa says by git branch
  _boosted_gamma = std::make_unique<TLorentzVector>(*_gamma);
  _boosted_pim_measured = std::make_unique<TLorentzVector>(*_mom_corr_pim);

  _rotated_prot = std::make_unique<TLorentzVector>(*_mom_corr_prot);
  _rotated_pip = std::make_unique<TLorentzVector>(*_mom_corr_pip);
  _rotated_pim = std::make_unique<TLorentzVector>(*_gamma + *_target - *_mom_corr_prot - *_mom_corr_pip);  //(*_pim);// careful here because it is not in exclusive set up sa says by git branch
  _rotated_pim_measured = std::make_unique<TLorentzVector>(*_mom_corr_pim);

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

float Reaction::pim_momentum_cm() {
  if (!_is_boosted) boost();
  // if (TwoPion_missingPim())
  if (TwoPion_exclusive())
    return _boosted_pim->P();
  else
    return NAN;
}

float Reaction::pim_theta_cm() {
  if (!_is_boosted) boost();
  // if (TwoPion_missingPim())
  if (TwoPion_exclusive())

    return _rotated_pim->Theta() * 180.0 / PI;
  else
    return NAN;
}

float Reaction::pim_Phi_cm() {
  if (!_is_boosted) boost();
  // if (TwoPion_missingPim()) {
  if (TwoPion_exclusive()) {
    if (_rotated_pim->Phi() > 0)
      return _rotated_pim->Phi() * 180 / PI;
    else if (_rotated_pim->Phi() < 0)
      return (_rotated_pim->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}

float Reaction::pim_momentum_cm_measured() {
  if (!_is_boosted) boost();
  if (TwoPion_exclusive())
    return _boosted_pim_measured->P();
  else
    return NAN;
}

float Reaction::pim_theta_cm_measured() {
  if (!_is_boosted) boost();
  if (TwoPion_exclusive())
    return _rotated_pim_measured->Theta() * 180.0 / PI;
  else
    return NAN;
}

float Reaction::pim_Phi_cm_measured() {
  if (!_is_boosted) boost();
  if (TwoPion_exclusive()) {
    if (_rotated_pim_measured->Phi() > 0)
      return _rotated_pim_measured->Phi() * 180 / PI;
    else if (_rotated_pim_measured->Phi() < 0)
      return (_rotated_pim_measured->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}

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

  _elec_mom_mc = _elec_mc->P();
  _elec_E_mc = _elec_mc->E();
  _theta_e_mc = _elec_mc->Theta() * 180 / PI;
}

void MCReaction::SetMCProton(int i) { _prot_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_P); }

void MCReaction::SetMCPip(int i) { _pip_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIP); }

void MCReaction::SetMCPim(int i) { _pim_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIM); }
// void MCReaction::SetMCOther(int i) {
//   _other_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i),
//   mass[_data->pid(i)]);
// }

float MCReaction::pim_mom_mc_gen() { return _pim_mc->P(); }
float MCReaction::pip_mom_mc_gen() { return _pip_mc->P(); }
float MCReaction::prot_mom_mc_gen() { return _prot_mc->P(); }

float MCReaction::pim_theta_mc_gen() { return _pim_mc->Theta() * 180 / PI; }
float MCReaction::pip_theta_mc_gen() { return _pip_mc->Theta() * 180 / PI; }
float MCReaction::prot_theta_mc_gen() { return _prot_mc->Theta() * 180 / PI; }

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
