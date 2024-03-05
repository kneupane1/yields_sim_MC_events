
#include "reaction.hpp"

Reaction::Reaction(const std::shared_ptr<Branches12>& data, float beam_energy) {
  _data = data;
  _beam = std::make_unique<TLorentzVector>();
  _beam_energy = beam_energy;
  _sector = data->dc_sec(0);

  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  _gamma = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  _elecUnSmear = std::make_unique<TLorentzVector>();
  // _elecSmear = std::make_unique<TLorentzVector>();
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

  _protUnSmear = std::make_unique<TLorentzVector>();
  _pipUnSmear = std::make_unique<TLorentzVector>();
  _pimUnSmear = std::make_unique<TLorentzVector>();

  _protSmear = std::make_unique<TLorentzVector>();
  _pipSmear = std::make_unique<TLorentzVector>();
  _pimSmear = std::make_unique<TLorentzVector>();

  _prot = std::make_unique<TLorentzVector>();
  _pip = std::make_unique<TLorentzVector>();
  _pim = std::make_unique<TLorentzVector>();
  _other = std::make_unique<TLorentzVector>();
  _neutron = std::make_unique<TLorentzVector>();
}

Reaction::~Reaction() {}

void Reaction::SetElec() {
  _hasE = true;
  _sectorElec = _data->dc_sec(0);
  _elec_status = abs(_data->status(0));

  _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);

  *_gamma += *_beam - *_elec;  // be careful you are commenting this only to include the momentum correction

  // // // // // Can calculate W and Q2 here (useful for simulations as sim do not have elec mom corrections)
  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);

  _elec_mom = _elec->P();
  _elec_E = _elec->E();
  _theta_e = _elec->Theta() * 180 / PI;

  // //////////////////////////////////////////////////////////////

  // _elecUnSmear->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);

  // double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear,
  // phiSmear;

  // pUnSmear = _elecUnSmear->P();

  // thetaUnSmear = _elecUnSmear->Theta() * 180 / PI;

  // if (_elecUnSmear->Phi() > 0)
  //   phiUnSmear = _elecUnSmear->Phi() * 180 / PI;
  // else if (_elecUnSmear->Phi() < 0)
  //   phiUnSmear = (_elecUnSmear->Phi() + 2 * PI) * 180 / PI;

  // ////////////////////////////////////////////////////////////////

  // // Generate new values
  // Reaction::SmearingFunc(ELECTRON, _elec_status, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear);

  // _pxPrimeSmear = _elecUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
  //                 sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
  // _pyPrimeSmear = _elecUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
  //                 sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
  // _pzPrimeSmear =
  //     _elecUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);

  // // _elecSmear->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_E);  // smeared
  // _elec->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_E);  // smeared
  // // _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);  // smeared

  // //////////////////////////////////////////////////////////////

  // // double pGen, pUnSmear, pSmear;
  // // Reaction::SmearingFunc1(pGen, pUnSmear, pSmear, 0.5);
  // // _pxPrimeSmear = _data->px(0) * ((pSmear) / (pUnSmear));
  // // _pyPrimeSmear = _data->py(0) * ((pSmear) / (pUnSmear));
  // // _pzPrimeSmear = _data->pz(0) * ((pSmear) / (pUnSmear));

  // // _elec->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_E);  // smeared

  // // //////////////////////////////////////////////////////////////
  // *_gamma += *_beam - *_elec;  // be careful you are commenting this only to include the momentum correction

  // // // // // // Can calculate W and Q2 here (useful for simulations as sim do not have elec mom corrections)
  // _W = physics::W_calc(*_beam, *_elec);
  // _Q2 = physics::Q2_calc(*_beam, *_elec);

  // _elec_mom = _elec->P();
  // _elec_E = _elec->E();
  // _theta_e = _elec->Theta() * 180 / PI;

  // if (_elec->Phi() > 0)
  //   _phi_elec = _elec->Phi() * 180 / PI;
  // else if (_elec->Phi() < 0)
  //   _phi_elec = (_elec->Phi() + 2 * PI) * 180 / PI;
}
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

void Reaction::SetProton(int i) {
  _numProt++;
  _numPos++;
  _hasP = true;
  _sectorProt = _data->dc_sec(i);
  _prot_status = abs(_data->status(i));

  _Energy_loss_uncorr_prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
  // _prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);

  _prot_mom_uncorr = _Energy_loss_uncorr_prot->P();

  _prot_theta_uncorr = _Energy_loss_uncorr_prot->Theta() * 180 / PI;

  // if (_Energy_loss_uncorr_prot->Phi() > 0)
  //   _prot_phi_uncorr = _Energy_loss_uncorr_prot->Phi() * 180 / PI;
  // else if (_Energy_loss_uncorr_prot->Phi() < 0)
  //   _prot_phi_uncorr = (_Energy_loss_uncorr_prot->Phi() + 2 * PI) * 180 / PI;

  // _thetaDC_r1_Prot = RAD2DEG * (acos(_data->dc_r1_z(i) / sqrt(pow(_data->dc_r1_x(i), 2) + pow(_data->dc_r1_y(i),
  // 2) +
  //                                                            pow(_data->dc_r1_z(i), 2))));

  // _thetaDC_r1_Prot = RAD2DEG * (atan2(sqrt(pow(_data->dc_r1_x(i), 2) + pow(_data->dc_r1_y(i), 2)),
  // _data->dc_r1_z(i)));

  _is_FD = mom_corr::is_FD(_prot_status);
  _is_CD = mom_corr::is_CD(_prot_status);
  // _is_lower_band = mom_corr::is_lower_band(_prot_mom_uncorr, _thetaDC_r1_Prot, _prot_status);

  if (_is_CD) {
    // std::cout << " this is cd " << _sectorProt << std::endl;

    // _prot_mom_tmt = _prot_mom_uncorr;
    // _prot_theta_tmt = _prot_theta_uncorr;
    // _prot_phi_tmt = _prot_phi_uncorr;

    _prot_mom_tmt = mom_corr::CD_prot_Emom_corr(_prot_mom_uncorr, _prot_theta_uncorr);
    // _prot_theta_tmt = mom_corr::CD_prot_Eth_corr(_prot_mom_uncorr, _prot_theta_uncorr);
    // _prot_phi_tmt = mom_corr::CD_prot_Eph_corr(_prot_mom_uncorr, _prot_theta_uncorr, _prot_phi_uncorr);
  }
  if (_is_FD) {
    // std::cout << " this is fd " << _sectorProt << std::endl;
    // _prot_mom_tmt = _prot_mom_uncorr;
    // // // these are Andrey's corrections
    if (_prot_theta_uncorr < 27) {
      // _prot_theta_tmt = _prot_theta_uncorr;
      // _prot_phi_tmt = _prot_phi_uncorr;
      _prot_mom_tmt = _prot_mom_uncorr + exp(-2.739 - 3.932 * _prot_theta_uncorr) + 0.002907;
    } else {
      // _prot_theta_tmt = _prot_theta_uncorr;
      // _prot_phi_tmt = _prot_phi_uncorr;
      _prot_mom_tmt = _prot_mom_uncorr + exp(-1.2 - 4.228 * _prot_mom_uncorr) + 0.007502;
    }
  }

  _px_prime_prot_E = _data->px(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
  _py_prime_prot_E = _data->py(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
  _pz_prime_prot_E = _data->pz(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));

  _prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);  // energy loss corrected
  // _protUnSmear->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);  // energy loss corrected

  // //////////////////////////////////////////////////////////////
  // double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear,
  // phiSmear;

  // pUnSmear = _protUnSmear->P();

  // thetaUnSmear = _protUnSmear->Theta() * 180 / PI;

  // if (_protUnSmear->Phi() > 0)
  //   phiUnSmear = _protUnSmear->Phi() * 180 / PI;
  // else if (_protUnSmear->Phi() < 0)
  //   phiUnSmear = (_protUnSmear->Phi() + 2 * PI) * 180 / PI;

  // // Generate new values

  // Reaction::SmearingFunc(PROTON, _prot_status, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear);

  // _pxPrimeSmear = _protUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
  //                 sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
  // _pyPrimeSmear = _protUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
  //                 sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
  // _pzPrimeSmear =
  //     _protUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);

  // // _protSmear->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_P);  // smeared
  // _prot->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_P);  // smeared
}
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

void Reaction::SetPip(int i) {
  _numPip++;
  _numPos++;
  _hasPip = true;
  _pip_status = abs(_data->status(i));
  _sectorPip = _data->dc_sec(i);
  // _thetaDC_r1_Pip = RAD2DEG * (atan2(sqrt(pow(_data->dc_r1_x(i), 2) + pow(_data->dc_r1_y(i), 2)),
  // _data->dc_r1_z(i)));

  _Energy_loss_uncorr_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
  // // // _mom_corr_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
  // _pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);

  _pip_mom_uncorr = _Energy_loss_uncorr_pip->P();
  _pip_theta_uncorr = _Energy_loss_uncorr_pip->Theta() * 180 / PI;
  // // if (_Energy_loss_uncorr_pip->Phi() > 0)
  // //   _pip_phi_uncorr = _Energy_loss_uncorr_pip->Phi() * 180 / PI;
  // // else if (_Energy_loss_uncorr_pip->Phi() < 0)
  // //   _pip_phi_uncorr = (_Energy_loss_uncorr_pip->Phi() + 2 * PI) * 180 / PI;

  _is_FD = mom_corr::is_FD(_pip_status);
  _is_CD = mom_corr::is_CD(_pip_status);
  // // _is_lower_band = mom_corr::is_lower_band(_pip_mom_uncorr, _thetaDC_r1_Pip, _pip_status);

  if (_is_CD) {
    // //   _pip_mom_tmt = _pip_mom_uncorr;
    // //   // _pip_theta_tmt = _pip_theta_uncorr;
    // //   // _pip_phi_tmt = _pip_phi_uncorr;

    _pip_mom_tmt = mom_corr::CD_pip_Emom_corr(_pip_mom_uncorr, _pip_theta_uncorr);
    // //   // _pip_theta_tmt = mom_corr::CD_pip_Eth_corr(_pip_mom_uncorr, _pip_theta_uncorr);
    // //   // _pip_phi_tmt = mom_corr::CD_pip_Eph_corr(_pip_mom_uncorr, _pip_theta_uncorr, _pip_phi_uncorr);
  }
  if (_is_FD) {
    _pip_mom_tmt = _pip_mom_uncorr;
    // //   _pip_theta_tmt = _pip_theta_uncorr;
    // //   _pip_phi_tmt = _pip_phi_uncorr;
  }

  // // // _pip->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_P);

  _px_prime_pip_E = _data->px(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
  _py_prime_pip_E = _data->py(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
  _pz_prime_pip_E = _data->pz(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));

  _pip->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);
  // _pipUnSmear->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);

  // // //////////////////////////////////////////////////////////////
  // double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear,
  // phiSmear;

  // pUnSmear = _pipUnSmear->P();

  // thetaUnSmear = _pipUnSmear->Theta() * 180 / PI;

  // if (_pipUnSmear->Phi() > 0)
  //   phiUnSmear = _pipUnSmear->Phi() * 180 / PI;
  // else if (_pipUnSmear->Phi() < 0)
  //   phiUnSmear = (_pipUnSmear->Phi() + 2 * PI) * 180 / PI;

  // // Generate new values
  // Reaction::SmearingFunc(PIP, _pip_status, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear);

  // _pxPrimeSmear = _pipUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
  //                 sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
  // _pyPrimeSmear = _pipUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
  //                 sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
  // _pzPrimeSmear = _pipUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD *
  // thetaUnSmear);

  // // _pipSmear->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIP);  // smeared
  // _pip->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIP);  // smeared
}
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

void Reaction::SetPim(int i) {
  _numPim++;
  _numNeg++;
  _hasPim = true;
  _pim_status = abs(_data->status(i));

  _sectorPim = _data->dc_sec(i);
  // std::cout << " pim status " << _pim_status << std::endl;
  // _thetaDC_r1_Pim = RAD2DEG * (atan2(sqrt(pow(_data->dc_r1_x(i), 2) + pow(_data->dc_r1_y(i), 2)),
  // _data->dc_r1_z(i)));

  _Energy_loss_uncorr_pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
  // _pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);

  _pim_mom_uncorr = _Energy_loss_uncorr_pim->P();
  _pim_theta_uncorr = _Energy_loss_uncorr_pim->Theta() * 180 / PI;
  // // // // if (_Energy_loss_uncorr_pim->Phi() > 0)
  // // // //   _pim_phi_uncorr = _Energy_loss_uncorr_pim->Phi() * 180 / PI;
  // // // // else if (_Energy_loss_uncorr_pim->Phi() < 0)
  // // // //   _pim_phi_uncorr = (_Energy_loss_uncorr_pim->Phi() + 2 * PI) * 180 / PI;

  _is_FD = mom_corr::is_FD(_pim_status);
  _is_CD = mom_corr::is_CD(_pim_status);
  // // // _is_lower_band = mom_corr::is_lower_band(_pim_mom_uncorr, _thetaDC_r1_Pim, _pim_status);

  if (_is_CD) {
    // _pim_mom_tmt = _pim_mom_uncorr;
    // _pim_theta_tmt = _pim_theta_uncorr;
    // _pim_phi_tmt = _pim_phi_uncorr;
    _pim_mom_tmt = mom_corr::CD_pim_Emom_corr(_pim_mom_uncorr, _pim_theta_uncorr);
    // _pim_theta_tmt = mom_corr::CD_pim_Eth_corr(_pim_mom_uncorr, _pim_theta_uncorr);
    // _pim_phi_tmt = mom_corr::CD_pim_Eph_corr(_pim_mom_uncorr, _pim_theta_uncorr, _pim_phi_uncorr);
  }
  if (_is_FD) {
    _pim_mom_tmt = _pim_mom_uncorr;
    // _pim_theta_tmt = _pim_theta_uncorr;
    // _pim_phi_tmt = _pim_phi_uncorr;
  }

  _px_prime_pim_E = _data->px(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
  _py_prime_pim_E = _data->py(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
  _pz_prime_pim_E = _data->pz(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));

  _pim->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, MASS_PIM);
  // _pimUnSmear->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, MASS_PIM);

  // // //////////////////////////////////////////////////////////////
  // double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear,
  // phiSmear;

  // pUnSmear = _pimUnSmear->P();

  // thetaUnSmear = _pimUnSmear->Theta() * 180 / PI;

  // if (_pimUnSmear->Phi() > 0)
  //   phiUnSmear = _pimUnSmear->Phi() * 180 / PI;
  // else if (_pimUnSmear->Phi() < 0)
  //   phiUnSmear = (_pimUnSmear->Phi() + 2 * PI) * 180 / PI;

  // // Generate new values
  // Reaction::SmearingFunc(PIM, _pim_status, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear);

  // _pxPrimeSmear = _pimUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
  //                 sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
  // _pyPrimeSmear = _pimUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
  //                 sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
  // _pzPrimeSmear = _pimUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD *
  // thetaUnSmear);

  // // _pimSmear->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIM);  // smeared
  // _pim->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIM);  // smeared
}
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

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
  auto mm_mPim = std::make_unique<TLorentzVector>();
  auto mm_mpip = std::make_unique<TLorentzVector>();
  auto mm_mprot = std::make_unique<TLorentzVector>();
  auto mm_excl = std::make_unique<TLorentzVector>();

  *mm_mPim += (*_gamma + *_target);

  if (TwoPion_exclusive()) {
    // *mm_mPim -= *_mom_corr_prot;
    // *mm_mPim -= *_mom_corr_pip;
    // // *mm_mPim -= *_pim;
    // _MM_mPim = mm->M();
    // _MM2_mPim = mm->M2();

    // *mm_excl += (*_gamma + *_target);
    // *mm_excl -= *_mom_corr_prot;
    // *mm_excl -= *_mom_corr_pip;
    // *mm_excl -= *_mom_corr_pim;

    *mm_mPim -= *_prot;
    *mm_mPim -= *_pip;
    // *mm_mPim -= *_pim;
    _MM_mPim = mm_mPim->M();
    _MM2_mPim = mm_mPim->M2();

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

  // if (TwoPion_missingPim()) {
  if (TwoPion_exclusive()) {
    auto missingpim_ = std::make_unique<TLorentzVector>();
    *missingpim_ += *_gamma + *_target - *_prot - *_pip;

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
// float Reaction::pim_momentum_corrected() {
//   if (TwoPion_exclusive())
//     return _mom_corr_pim->P();
//   else
//     return NAN;
// }
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

    return missingpip_->P();
  } else
    return NAN;
}
float Reaction::pip_theta_lab() {
  // if (TwoPion_missingPip()) {
  if (TwoPion_exclusive()) {
    auto missingpip_ = std::make_unique<TLorentzVector>();
    *missingpip_ += *_gamma + *_target - *_prot - *_pim;
    return missingpip_->Theta() * 180.0 / PI;
  } else
    return NAN;
}
float Reaction::pip_Phi_lab() {
  // if (TwoPion_missingPip()) {
  if (TwoPion_exclusive()) {
    auto missingpip_ = std::make_unique<TLorentzVector>();
    *missingpip_ += *_gamma + *_target - *_prot - *_pim;

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

////////////////mProt
float Reaction::prot_momentum() {
  // if (TwoPion_missingProt()) {
  if (TwoPion_exclusive()) {
    auto missingprot_ = std::make_unique<TLorentzVector>();
    *missingprot_ += *_gamma + *_target - *_pip - *_pim;

    return missingprot_->P();
  } else
    return NAN;
}
float Reaction::prot_theta_lab() {
  // if (TwoPion_missingProt()) {
  if (TwoPion_exclusive()) {
    auto missingprot_ = std::make_unique<TLorentzVector>();
    *missingprot_ += *_gamma + *_target - *_pip - *_pim;

    return missingprot_->Theta() * 180.0 / PI;
  } else
    return NAN;
}
float Reaction::prot_Phi_lab() {
  // if (TwoPion_missingProt()) {
  if (TwoPion_exclusive()) {
    auto missingprot_ = std::make_unique<TLorentzVector>();
    *missingprot_ += *_gamma + *_target - *_pip - *_pim;

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

/////////////////////

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
  TVector3 uz = _boosted_gamma->Vect().Unit();                  // uit vector along virtual photon
  TVector3 ux = ((_beam->Vect()).Cross(_elec->Vect())).Unit();  // unit vector along e cross e'
  ux.Rotate(3. * PI / 2, uz);                                   // rotating ux by 3pi/2 with uz as axis of roration
  rot.SetZAxis(uz, ux).Invert();                                // setting TRotation rot

  _boosted_gamma->Transform(rot);
  float_t beta_1 = ((sqrt(_boosted_gamma->E() * _boosted_gamma->E() + _Q2)) / (_boosted_gamma->E() + MASS_P));

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
    // return ((_prot->Vect()).Dot((_pip->Vect()).Cross(_pim->Vect())));

  } else
    return NAN;
}
//////////////////
void Reaction::invMassPpim() {
  if (!_is_boosted) boost();

  auto inv_Ppim = std::make_unique<TLorentzVector>();
  *inv_Ppim += *_boosted_prot;
  *inv_Ppim += *_boosted_pim_measured;
  // *inv_Ppim += (*_boosted_gamma + *_target - *_boosted_prot - *_boosted_pip);
  if (TwoPion_exclusive()) _inv_Ppim = inv_Ppim->M();
}
void Reaction::invMasspippim() {
  if (!_is_boosted) boost();

  auto inv_pip_pim = std::make_unique<TLorentzVector>();
  *inv_pip_pim += *_boosted_pip;
  *inv_pip_pim += *_boosted_pim_measured;

  // *inv_pip_pim += (*_boosted_gamma + *_target - *_boosted_prot - *_boosted_pip);
  if (TwoPion_exclusive()) _inv_pip_pim = inv_pip_pim->M();
}
void Reaction::invMassPpip() {
  if (!_is_boosted) boost();

  auto inv_Ppip = std::make_unique<TLorentzVector>();
  *inv_Ppip += *_boosted_prot;
  *inv_Ppip += *_boosted_pip;
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
////////////
float Reaction::pim_momentum_cm() {
  if (!_is_boosted) boost();
  if (TwoPion_exclusive())
    return _boosted_pim->P();
  else
    return NAN;
}

float Reaction::pim_theta_cm() {
  if (!_is_boosted) boost();
  if (TwoPion_exclusive())
    return _rotated_pim->Theta() * 180.0 / PI;
  else
    return NAN;
}

float Reaction::pim_Phi_cm() {
  if (!_is_boosted) boost();
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
  _beam_mc = std::make_unique<TLorentzVector>();
  _beam_energy = beam_enrgy;
  _weight_mc = _data->mc_weight();
  _beam_mc->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

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

  *_gamma_mc += *_beam_mc - *_elec_mc;

  // Can calculate W and Q2 here
  _W_mc = physics::W_calc(*_beam_mc, *_elec_mc);
  _Q2_mc = physics::Q2_calc(*_beam_mc, *_elec_mc);

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

// if (_beam_mc->Phi() >= 0)
//   _beam_phi_mc = (_beam_mc->Phi() * 180 / PI);
// else if (_beam_mc->Phi() < 0)
//   _beam_phi_mc = ((_beam_mc->Phi() + 2 * PI) * 180 / PI);

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

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

void MCReaction::boost_mc() {
  _is_boosted_mc = true;
  _boosted_prot_mc = std::make_unique<TLorentzVector>(*_prot_mc);
  _boosted_pip_mc = std::make_unique<TLorentzVector>(*_pip_mc);
  _boosted_pim_measured_mc = std::make_unique<TLorentzVector>(*_pim_mc);
  _boosted_pim_mc = std::make_unique<TLorentzVector>(*_gamma_mc + *_target - *_prot_mc - *_pip_mc);  //(*_pim);
  // _boosted_pim_mc = std::make_unique<TLorentzVector>(*_pim_mc);
  // std::cout << " P_pz " << (_pim_mc->Pz()) << "    ";
  // std::cout << (_boosted_pim_mc->Pz()) << std::endl;
  // std::cout << " Pz_diff " << (_pim_mc->Pz()) - (_boosted_pim_mc->Pz()) << "    ";
  // std::cout << " P " << (_pim_mc->P()) << "    ";
  // std::cout << " P_cal " <<(_boosted_pim_mc->P()) << "    ";
  // std::cout << " P_diff " << (_pim_mc->P()) - (_boosted_pim_mc->P()) << "    " << std::endl;

  // // if ((_boosted_pim_mc->Pz()) > 0.000001)
  // {
  //         std::cout
  //         //  << _beam_mc->Pz() << " " << _elec_mc->Pz() << " "
  //         //           << _gamma_mc->Pz() << " " << _prot_mc->Pz() << " "
  //         //           << _pip_mc->Pz() << " " << _pim_mc->Pz() << " "
  //                   << _boosted_pim_mc->Px() - _pim_mc->Px() << " "
  //                 //   << _boosted_pim_mc->P() - _pim_mc->P() << " "
  //                 //   << _boosted_pim_mc->E() - _pim_mc->E() << " "
  //                   << std::endl;
  // }

  _boosted_gamma_mc = std::make_unique<TLorentzVector>(*_gamma_mc);

  TRotation rot_mc;

  TVector3 uz_mc = _boosted_gamma_mc->Vect().Unit();                     // uit vector along virtual photon
  TVector3 ux_mc = ((_beam_mc->Vect()).Cross(_elec_mc->Vect())).Unit();  // unit vector along e cross e'
  ux_mc.Rotate(3. * PI / 2, uz_mc);        // rotating ux by 3pi/2 with uz as axis of roration
  rot_mc.SetZAxis(uz_mc, ux_mc).Invert();  // setting TRotation rot

  _boosted_gamma_mc->Transform(rot_mc);
  _boosted_prot_mc->Transform(rot_mc);
  _boosted_pip_mc->Transform(rot_mc);
  _boosted_pim_mc->Transform(rot_mc);
  _boosted_pim_measured_mc->Transform(rot_mc);

  float_t beta_1_mc =
      ((sqrt(_boosted_gamma_mc->E() * _boosted_gamma_mc->E() + _Q2_mc)) / (_boosted_gamma_mc->E() + MASS_P));

  _boosted_prot_mc->Boost(0, 0, -beta_1_mc);
  _boosted_pip_mc->Boost(0, 0, -beta_1_mc);
  _boosted_pim_mc->Boost(0, 0, -beta_1_mc);
  _boosted_gamma_mc->Boost(0, 0, -beta_1_mc);  // -beta ko value (0.5 to -0.5 huda
                                               // samma value aauchha nattra aaudyna)
  _boosted_pim_measured_mc->Boost(0, 0, -beta_1_mc);

  _prot_Vect3_mc = _boosted_prot_mc->Vect();
  _pip_Vect3_mc = _boosted_pip_mc->Vect();
  _pim_Vect3_mc_measured = _boosted_pim_measured_mc->Vect();
}

float MCReaction::scalar_triple_product_thrown() {
  if (!_is_boosted_mc) boost_mc();
  // if (TwoPion_exclusive()) {
  return (_prot_Vect3_mc.Dot(_pip_Vect3_mc.Cross(_pim_Vect3_mc_measured)));

  // lab system
  // return ((_prot_mc->Vect()).Dot((_pip_mc->Vect()).Cross(_pim_mc->Vect())));

  // }
  // return NAN;
}
//
// float MCReaction::MCgamma_Phi() {
//   if (!_is_boosted_mc) boost_mc();
//   // if (TwoPion_missingPim()) { //???? is this correct?
//   if (_gamma_mc->Phi() > 0)
//     return _gamma_mc->Phi() * 180 / PI;
//   else if (_gamma_mc->Phi() < 0)
//     return (_gamma_mc->Phi() + 2 * PI) * 180 / PI;
//   else
//     return NAN;
//   //} else
//   // return NAN;
// }
// float MCReaction::MCprot_Phi() {
//   if (!_is_boosted_mc) boost_mc();
//   // if (TwoPion_missingPim()) {
//   if (_boosted_prot_mc->Phi() > 0)
//     return _boosted_prot_mc->Phi() * 180 / PI;
//   else if (_boosted_prot_mc->Phi() < 0)
//     return (_boosted_prot_mc->Phi() + 2 * PI) * 180 / PI;
//   else
//     return NAN;
//   //} else
//   // return NAN;
// }
// float MCReaction::MCpip_Phi() {
//   // if (TwoPion_missingPim()) {
//   if (_boosted_pip_mc->Phi() > 0)
//     return _boosted_pip_mc->Phi() * 180 / PI;
//   else if (_boosted_pip_mc->Phi() < 0)
//     return (_boosted_pip_mc->Phi() + 2 * PI) * 180 / PI;
//   else
//     return NAN;
//   //} else
//   // return NAN;
// }
// float MCReaction::MCpim_Phi() {
//   if (!_is_boosted_mc) boost_mc();
//   // if (TwoPion_missingPim()) {
//   if (_boosted_pim_mc->Phi() > 0)
//     return _boosted_pim_mc->Phi() * 180 / PI;
//   else if (_boosted_pim_mc->Phi() < 0)
//     return (_boosted_pim_mc->Phi() + 2 * PI) * 180 / PI;
//   else
//     return NAN;
//   //} else
//   // return NAN;
// }
// float MCReaction::MCprot_theta() {
//   if (!_is_boosted_mc) boost_mc();
//   // if (TwoPion_missingPim())
//   return _boosted_prot_mc->Theta() * 180.0 / PI;
//   // else
//   // return NAN;
// }
// float MCReaction::MCpip_theta() {
//   if (!_is_boosted_mc) boost_mc();
//   // if (TwoPion_missingPim())
//   return _boosted_pip_mc->Theta() * 180.0 / PI;
//   // else
//   // return NAN;
// }
// float MCReaction::MCpim_theta() {
//   if (!_is_boosted_mc) boost_mc();
//   // if (TwoPion_missingPim())
//   return _boosted_pim_mc->Theta() * 180.0 / PI;
//   // else
//   // return NAN;
// }

void MCReaction::MCinvMassPpip() {
  if (!_is_boosted_mc) boost_mc();
  auto MCinv_Ppip = std::make_unique<TLorentzVector>();
  *MCinv_Ppip += *_boosted_prot_mc;
  *MCinv_Ppip += *_boosted_pip_mc;
  _MCinv_Ppip = MCinv_Ppip->M();
}

void MCReaction::MCinvMassPpim() {
  // std::cout << *_pim_mc->M() -  << std::endl;

  if (!_is_boosted_mc) boost_mc();
  auto MCinv_Ppim = std::make_unique<TLorentzVector>();
  *MCinv_Ppim += *_boosted_prot_mc;
  *MCinv_Ppim += *_boosted_pim_mc;
  _MCinv_Ppim = MCinv_Ppim->M();
}
void MCReaction::MCinvMasspippim() {
  if (!_is_boosted_mc) boost_mc();
  auto MCinv_pip_pim = std::make_unique<TLorentzVector>();
  *MCinv_pip_pim += *_boosted_pip_mc;
  *MCinv_pip_pim += *_boosted_pim_mc;
  _MCinv_pip_pim = MCinv_pip_pim->M();
}
float MCReaction::MCinv_Ppip() {
  if (_MCinv_Ppip != _MCinv_Ppip) MCinvMassPpip();
  // if (TwoPion_missingPim())
  if (_MCinv_Ppip != NAN) return _MCinv_Ppip;
  // else
  return NAN;
}
float MCReaction::MCinv_Ppim() {
  if (_MCinv_Ppim != _MCinv_Ppim) MCinvMassPpim();
  // if (TwoPion_missingPim())
  if (_MCinv_Ppim != NAN) return _MCinv_Ppim;
  //  else
  return NAN;
}
float MCReaction::MCinv_pip_pim() {
  if (_MCinv_pip_pim != _MCinv_pip_pim) MCinvMasspippim();
  // if (TwoPion_missingPim())
  if (_MCinv_pip_pim != NAN) return _MCinv_pip_pim;
  // else
  return NAN;
}

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
