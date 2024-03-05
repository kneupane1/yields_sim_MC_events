
#include "reaction.hpp"
#include <TMath.h>

Reaction::Reaction(const std::shared_ptr<Branches12>& data, float beam_energy) {
  _data = data;
  _beam = std::make_unique<TLorentzVector>();
  _beam_energy = beam_energy;
  _sector = data->dc_sec(0);

  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  _gamma = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  _elecUnSmear = std::make_unique<TLorentzVector>();
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

  _protUnSmear = std::make_unique<TLorentzVector>();
  _pipUnSmear = std::make_unique<TLorentzVector>();
  _pimUnSmear = std::make_unique<TLorentzVector>();

  _residualXpcal = NAN;
  _residualYpcal = NAN;
  _residualZpcal = NAN;

  _residualXecin = NAN;
  _residualYecin = NAN;
  _residualZecin = NAN;
}

Reaction::~Reaction() {}

void Reaction::SetElec() {
  _hasE = true;
  _sectorElec = _data->dc_sec(0);
  _elec_status = abs(_data->status(0));

  _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);
  *_gamma += *_beam - *_elec;  // be careful you are commenting this only to include the momentum correction

  // // Can calculate W and Q2 here (useful for simulations as sim do not have elec mom corrections)
  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);
  _elec_mom = _elec->P();

  // // //////////////////////////////////////////////////////////////
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

  // *_gamma += *_beam - *_elec;  // be careful you are commenting this only to include the momentum correction

  // // // // // // Can calculate W and Q2 here (useful for simulations as sim do not have elec mom corrections)
  // _W = physics::W_calc(*_beam, *_elec);
  // _Q2 = physics::Q2_calc(*_beam, *_elec);

  // _elec_mom = _elec->P();
  // _elec_E = _elec->E();
  // _theta_e = _elec->Theta() * 180 / PI;

  // // if (_elec->Phi() > 0)
  // //   _phi_elec = _elec->Phi() * 180 / PI;
  // // else if (_elec->Phi() < 0)
  // //   _phi_elec = (_elec->Phi() + 2 * PI) * 180 / PI;
}

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

  // // _thetaDC_r1_Prot = RAD2DEG * (acos(_data->dc_r1_z(i) / sqrt(pow(_data->dc_r1_x(i), 2) + pow(_data->dc_r1_y(i),
  // // 2) +
  // //                                                            pow(_data->dc_r1_z(i), 2))));

  // _thetaDC_r1_Prot = RAD2DEG * (atan2(sqrt(pow(_data->dc_r1_x(i), 2) + pow(_data->dc_r1_y(i), 2)),
  // _data->dc_r1_z(i)));

  _is_FD = mom_corr::is_FD(_prot_status);
  _is_CD = mom_corr::is_CD(_prot_status);
  // _is_lower_band = mom_corr::is_lower_band(_prot_mom_uncorr, _thetaDC_r1_Prot, _prot_status);

  if (_is_CD) {
    // _prot_mom_tmt = _prot_mom_uncorr;
    // _prot_theta_tmt = _prot_theta_uncorr;
    // _prot_phi_tmt = _prot_phi_uncorr;

    _prot_mom_tmt = mom_corr::CD_prot_Emom_corr(_prot_mom_uncorr, _prot_theta_uncorr);
    // _prot_theta_tmt = mom_corr::CD_prot_Eth_corr(_prot_mom_uncorr, _prot_theta_uncorr);
    // _prot_phi_tmt = mom_corr::CD_prot_Eph_corr(_prot_mom_uncorr, _prot_theta_uncorr, _prot_phi_uncorr);
  }
  if (_is_FD) {
    // // these are Andrey's corrections
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
  // // _mom_corr_prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);

  // /////////////////// SMEARING PART ////////////////////////////////////////////////////////////////////////////

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
  // _mom_corr_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
  // _pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);

  _pip_mom_uncorr = _Energy_loss_uncorr_pip->P();
  _pip_theta_uncorr = _Energy_loss_uncorr_pip->Theta() * 180 / PI;
  // if (_Energy_loss_uncorr_pip->Phi() > 0)
  //   _pip_phi_uncorr = _Energy_loss_uncorr_pip->Phi() * 180 / PI;
  // else if (_Energy_loss_uncorr_pip->Phi() < 0)
  //   _pip_phi_uncorr = (_Energy_loss_uncorr_pip->Phi() + 2 * PI) * 180 / PI;

  _is_FD = mom_corr::is_FD(_pip_status);
  _is_CD = mom_corr::is_CD(_pip_status);
  // _is_lower_band = mom_corr::is_lower_band(_pip_mom_uncorr, _thetaDC_r1_Pip, _pip_status);

  if (_is_CD) {
    // _pip_mom_tmt = _pip_mom_uncorr;
    // _pip_theta_tmt = _pip_theta_uncorr;
    // _pip_phi_tmt = _pip_phi_uncorr;

    _pip_mom_tmt = mom_corr::CD_pip_Emom_corr(_pip_mom_uncorr, _pip_theta_uncorr);
    // _pip_theta_tmt = mom_corr::CD_pip_Eth_corr(_pip_mom_uncorr, _pip_theta_uncorr);
    // _pip_phi_tmt = mom_corr::CD_pip_Eph_corr(_pip_mom_uncorr, _pip_theta_uncorr, _pip_phi_uncorr);
  }
  if (_is_FD) {
    _pip_mom_tmt = _pip_mom_uncorr;
    // _pip_theta_tmt = _pip_theta_uncorr;
    // _pip_phi_tmt = _pip_phi_uncorr;
  }

  // // _pip->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_P);

  _px_prime_pip_E = _data->px(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
  _py_prime_pip_E = _data->py(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
  _pz_prime_pip_E = _data->pz(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));

  _pip->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);  // cd eloss corrected

  // // /////////////////////////////////     SMEARING PART  /////////////////////////////

  // _pipUnSmear->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);

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
  // // // // _is_lower_band = mom_corr::is_lower_band(_pim_mom_uncorr, _thetaDC_r1_Pim, _pim_status);

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

  // // /////////////////////////////////     SMEARING PART  /////////////////////////////

  // _pimUnSmear->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, MASS_PIM);

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
  auto mm_mpim = std::make_unique<TLorentzVector>();
  *mm_mpim += (*_gamma + *_target);
  auto mm_mpip = std::make_unique<TLorentzVector>();
  *mm_mpip += (*_gamma + *_target);
  auto mm_mprot = std::make_unique<TLorentzVector>();
  *mm_mprot += (*_gamma + *_target);

  if (TwoPion_missingPim()) {
    *mm_mpim -= *_prot;
    *mm_mpim -= *_pip;
    // *mm -= *_pim;
    _MM_mPim = mm_mpim->M();
    _MM2_mPim = mm_mpim->M2();
    //   // _rec_pim_mom = mm->P();
    //   // _rec_pim_theta = mm->Theta() * 180 / PI;

    //   // if (mm->Phi() >= 0)
    //   //   _rec_pim_phi = (mm->Phi() * 180 / PI);
    //   // else if (mm->Phi() < 0)
    //   //   _rec_pim_phi = ((mm->Phi() + 2 * PI) * 180 / PI);
  }

  if (TwoPion_missingPip()) {
    *mm_mpip -= *_prot;
    *mm_mpip -= *_pim;
    _MM2_mPip = mm_mpip->M2();
  }

  if (TwoPion_missingProt()) {
    *mm_mprot -= *_pip;
    *mm_mprot -= *_pim;
    _MM2_mProt = mm_mprot->M2();
  }
}

float Reaction::MM_mPim() {
  if (_MM_mPim != _MM_mPim) CalcMissMass();
  return _MM_mPim;
}
float Reaction::MM2_mPim() {
  if (_MM2_mPim != _MM2_mPim) CalcMissMass();
  return _MM2_mPim;
}
float Reaction::MM2_mPip() {
  if (_MM2_mPip != _MM2_mPip) CalcMissMass();
  return _MM2_mPip;
}
float Reaction::MM2_mProt() {
  if (_MM2_mProt != _MM2_mProt) CalcMissMass();
  return _MM2_mProt;
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
// float Reaction::w_hadron() {
//   if (TwoPion_exclusive())
//     return ((*_prot) + (*_pip) + (*_pim)).Mag();
//   else
//     return NAN;
// }

/// pip topology

float Reaction::pip_momentum() {
  // if (_rec_pip_mom != _rec_pip_mom) CalcMissMass();

  if (TwoPion_missingPip()) {
    // if (TwoPion_exclusive()) {
    auto missingpip_ = std::make_unique<TLorentzVector>();
    *missingpip_ += *_gamma + *_target - *_prot - *_pim;
    // *missingpip_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pim;

    return missingpip_->P();
    // return _rec_pip_mom;

  } else
    return NAN;
}
float Reaction::pip_theta_lab() {
  // if (_rec_pip_theta != _rec_pip_theta) CalcMissMass();

  if (TwoPion_missingPip()) {
    // if (TwoPion_exclusive()) {
    auto missingpip_ = std::make_unique<TLorentzVector>();
    *missingpip_ += *_gamma + *_target - *_prot - *_pim;
    // *missingpip_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pim;

    return missingpip_->Theta() * 180.0 / PI;
    // return _rec_pip_theta;
  } else
    return NAN;
}
float Reaction::pip_Phi_lab() {
  // if (_rec_pip_phi != _rec_pip_phi) CalcMissMass();

  if (TwoPion_missingPip()) {
    // if (TwoPion_exclusive()) {
    auto missingpip_ = std::make_unique<TLorentzVector>();
    *missingpip_ += *_gamma + *_target - *_prot - *_pim;
    // *missingpip_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pim;

    if (missingpip_->Phi() > 0)
      return missingpip_->Phi() * 180 / PI;
    else if (missingpip_->Phi() < 0)
      return (missingpip_->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
    // return _rec_pip_phi;
  } else
    return NAN;
}

////////////////mProt
float Reaction::prot_momentum() {
  if (TwoPion_missingProt()) {
    // if (TwoPion_exclusive()) {
    auto missingprot_ = std::make_unique<TLorentzVector>();
    *missingprot_ += *_gamma + *_target - *_pip - *_pim;
    // 4-vect approach
    // *missingprot_ += *_prot -(*_gamma + *_target - *_pip - *_pim);

    return missingprot_->P();
  } else
    return NAN;
}
float Reaction::prot_theta_lab() {
  if (TwoPion_missingProt()) {
    // if (TwoPion_exclusive()) {
    auto missingprot_ = std::make_unique<TLorentzVector>();
    *missingprot_ += *_gamma + *_target - *_pip - *_pim;
    // *missingprot_ += *_gamma + *_target - *_mom_corr_pip - *_mom_corr_pim;

    return missingprot_->Theta() * 180.0 / PI;
  } else
    return NAN;
}
float Reaction::prot_Phi_lab() {
  if (TwoPion_missingProt()) {
    // if (TwoPion_exclusive()) {
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
  float_t beta_1 = ((sqrt(_boosted_gamma->E() * _boosted_gamma->E() + _Q2)) / (_boosted_gamma->E() + MASS_P));
  TVector3 uz = _boosted_gamma->Vect().Unit();                  // uit vector along virtual photon
  TVector3 ux = ((_beam->Vect()).Cross(_elec->Vect())).Unit();  // unit vector along e cross e'
  ux.Rotate(3. * PI / 2, uz);                                   // rotating ux by 3pi/2 with uz as axis of roration
  rot.SetZAxis(uz, ux).Invert();                                // setting TRotation rot

  _boosted_gamma->Transform(rot);

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
}

void Reaction::invMassPpim() {
  if (!_is_boosted) boost();
  auto inv_Ppim = std::make_unique<TLorentzVector>();
  *inv_Ppim += *_boosted_prot;
  *inv_Ppim += *_boosted_pim;
  if (TwoPion_missingPip()) _inv_Ppim = inv_Ppim->M();
}
void Reaction::invMasspippim() {
  if (!_is_boosted) boost();
  auto inv_pip_pim = std::make_unique<TLorentzVector>();
  *inv_pip_pim += *_boosted_pip;
  *inv_pip_pim += *_boosted_pim;
  if (TwoPion_missingProt()) _inv_pip_pim = inv_pip_pim->M();
}
void Reaction::invMassPpip() {
  if (!_is_boosted) boost();
  auto inv_Ppip = std::make_unique<TLorentzVector>();
  *inv_Ppip += *_boosted_prot;
  *inv_Ppip += *_boosted_pip;
  if (TwoPion_missingPim()) _inv_Ppip = inv_Ppip->M();
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

float Reaction::pim_momentum_cm() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim())
    return _boosted_pim->P();
  else
    return NAN;
}

float Reaction::pim_theta_cm() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim())
    return _rotated_pim->Theta() * 180.0 / PI;
  else
    return NAN;
}

float Reaction::pim_Phi_cm() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim()) {
    if (_rotated_pim->Phi() > 0)
      return _rotated_pim->Phi() * 180 / PI;
    else if (_rotated_pim->Phi() < 0)
      return (_rotated_pim->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
/////////////
float Reaction::prot_theta() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim()) return _boosted_prot->Theta() * 180.0 / PI;
  // else
  return NAN;
}
float Reaction::pip_theta() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim()) return _boosted_pip->Theta() * 180.0 / PI;
  // else
  return NAN;
}
float Reaction::pim_theta() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim()) return _boosted_pim->Theta() * 180.0 / PI;
  // else
  return NAN;
}

float Reaction::gamma_Phi() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim()) {
    if (_boosted_gamma->Phi() > 0)
      return _boosted_gamma->Phi() * 180 / PI;
    else if (_boosted_gamma->Phi() < 0)
      return (_boosted_gamma->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
float Reaction::prot_Phi() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim()) {
    if (_boosted_prot->Phi() > 0)
      return _boosted_prot->Phi() * 180 / PI;
    else if (_boosted_prot->Phi() < 0)
      return (_boosted_prot->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}

float Reaction::pip_Phi() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim()) {
    if (_boosted_pip->Phi() > 0)
      return _boosted_pip->Phi() * 180 / PI;
    else if (_boosted_pip->Phi() < 0)
      return (_boosted_pip->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}
float Reaction::pim_Phi() {
  if (!_is_boosted) boost();
  if (TwoPion_missingPim()) {
    if (_boosted_pim->Phi() > 0)
      return _boosted_pim->Phi() * 180 / PI;
    else if (_boosted_pim->Phi() < 0)
      return (_boosted_pim->Phi() + 2 * PI) * 180 / PI;
    else
      return NAN;
  } else
    return NAN;
}

void Reaction::AlphaCalc() {
  //  Float_t m_proton, m_pip, beta;
  Float_t a_gamma, b_gamma, a_beta, b_beta;
  TVector3 Vect3_gamma, Vect3_beta, V3_anti_z(0, 0, -1);
  float alpha_PPIp_piPIm;  // proton initial pim
  float alpha_PIpPIm_pipf;
  float alpha_PPIm_piPIp;

  if (!_is_boosted) boost();

  // 1 this one is used for α[π−]
  a_gamma = sqrt(1. / (1 - pow((_boosted_pim->Vect().Unit() * V3_anti_z),
                               2)));  // V3_anti_z(0,0,-1);
  b_gamma = -(_boosted_pim->Vect().Unit() * V3_anti_z) * a_gamma;
  Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _boosted_pim->Vect().Unit();

  a_beta = sqrt(1. / (1 - pow((_boosted_pim->Vect().Unit() * _boosted_pip->Vect().Unit()), 2)));
  b_beta = -(_boosted_pim->Vect().Unit() * _boosted_pip->Vect().Unit()) * a_beta;
  Vect3_beta = a_beta * _boosted_pip->Vect().Unit() + b_beta * _boosted_pim->Vect().Unit();

  alpha_PPIp_piPIm = (180. / PI) * acos(Vect3_gamma * Vect3_beta);
  if (Vect3_gamma.Cross(Vect3_beta) * _boosted_pim->Vect() < 0) alpha_PPIp_piPIm = 360. - alpha_PPIp_piPIm;

  //α[pπ+][p'π−]
  /// 2
  a_gamma = sqrt(1. / (1 - pow((_boosted_prot->Vect().Unit() * V3_anti_z), 2)));
  b_gamma = -(_boosted_prot->Vect().Unit() * V3_anti_z) * a_gamma;
  Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _boosted_prot->Vect().Unit();

  a_beta = sqrt(1. / (1 - pow((_boosted_prot->Vect().Unit() * _boosted_pip->Vect().Unit()), 2)));
  b_beta = -(_boosted_prot->Vect().Unit() * _boosted_pip->Vect().Unit()) * a_beta;
  Vect3_beta = a_beta * _boosted_pip->Vect().Unit() + b_beta * _boosted_prot->Vect().Unit();

  alpha_PIpPIm_pipf = (180. / PI) * acos(Vect3_gamma * Vect3_beta);

  if (Vect3_gamma.Cross(Vect3_beta) * _boosted_prot->Vect() < 0) alpha_PIpPIm_pipf = 360. - alpha_PIpPIm_pipf;
  //α[pp'][π+π−]

  /// 3
  a_gamma = sqrt(1. / (1 - pow((_boosted_pip->Vect().Unit() * V3_anti_z), 2)));
  b_gamma = -(_boosted_pip->Vect().Unit() * V3_anti_z) * a_gamma;
  Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _boosted_pip->Vect().Unit();

  a_beta = sqrt(1. / (1 - pow((_boosted_pip->Vect().Unit() * _boosted_pim->Vect().Unit()), 2)));
  b_beta = -(_boosted_pip->Vect().Unit() * _boosted_pim->Vect().Unit()) * a_beta;
  Vect3_beta = a_beta * _boosted_pim->Vect().Unit() + b_beta * _boosted_pip->Vect().Unit();

  alpha_PPIm_piPIp = (180. / PI) * acos(Vect3_gamma * Vect3_beta);

  if (Vect3_gamma.Cross(Vect3_beta) * _boosted_pip->Vect() < 0) alpha_PPIm_piPIp = 360. - alpha_PPIm_piPIp;

  _alpha_ppip_pipim = alpha_PPIp_piPIm;
  _alpha_pippim_pipf = alpha_PIpPIm_pipf;
  _alpha_ppim_pipip = alpha_PPIm_piPIp;
}

float Reaction::alpha_ppip_pipim() {  // pipim bhaneko proton initial  pim ho?
  if (_alpha_ppip_pipim != _alpha_ppip_pipim) AlphaCalc();
  if (TwoPion_missingPim())
    return _alpha_ppip_pipim;
  else
    return NAN;
}
float Reaction::alpha_pippim_pipf() {  // alpha P (proton initial proton final)
  if (_alpha_pippim_pipf != _alpha_pippim_pipf) AlphaCalc();
  if (TwoPion_missingPim())
    return _alpha_pippim_pipf;
  else
    return NAN;
}
float Reaction::alpha_ppim_pipip() {  // alpha pip (proton initial pip)
  if (_alpha_ppim_pipip != _alpha_ppim_pipip) AlphaCalc();
  if (TwoPion_missingPim())
    return _alpha_ppim_pipip;
  else
    return NAN;
}

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
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////

MCReaction::MCReaction(const std::shared_ptr<Branches12>& data, float beam_energy) {
  _data = data;
  if (!_data->mc()) _data->mc_branches();
  _beam = std::make_unique<TLorentzVector>();
  _beam_energy = beam_energy;
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

//// boost here ...

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
