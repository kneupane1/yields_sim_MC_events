
#include "reaction.hpp"

Reaction::Reaction(const std::shared_ptr<Branches12> &data, float beam_energy) {
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

  _prot = std::vector<std::unique_ptr<TLorentzVector>>();  // Initialize _protons as an empty vector
  _pip = std::vector<std::unique_ptr<TLorentzVector>>();
  _pim = std::vector<std::unique_ptr<TLorentzVector>>();

  _mom_corr_elec = std::make_unique<TLorentzVector>();
  _mom_corr_prot = std::vector<std::unique_ptr<TLorentzVector>>();
  _mom_corr_pip = std::vector<std::unique_ptr<TLorentzVector>>();
  _mom_corr_pim = std::vector<std::unique_ptr<TLorentzVector>>();

  _prot_indices = std::vector<int>();
  _pip_indices = std::vector<int>();
  _pim_indices = std::vector<int>();

  _Energy_loss_uncorr_prot = std::make_unique<TLorentzVector>();
  _Energy_loss_uncorr_pip = std::make_unique<TLorentzVector>();
  _Energy_loss_uncorr_pim = std::make_unique<TLorentzVector>();
  _pim_tmt = std::make_unique<TLorentzVector>();
  _pip_tmt = std::make_unique<TLorentzVector>();

  _other = std::make_unique<TLorentzVector>();
  _neutron = std::make_unique<TLorentzVector>();

  _protUnSmear = std::make_unique<TLorentzVector>();
  _pipUnSmear = std::make_unique<TLorentzVector>();
  _pimUnSmear = std::make_unique<TLorentzVector>();
}

Reaction::~Reaction() {}
auto objMomCorr = std::make_shared<mom_corr>();

void Reaction::SetElec() {
  _hasE = true;
  // _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);
  // *_gamma += *_beam - *_elec;  // be careful you are commenting this only to include the momentum correction

  // // // // // // // Can calculate W and Q2 here (useful for simulations as sim do not have elec mom corrections)
  // _W = physics::W_calc(*_beam, *_elec);
  // _Q2 = physics::Q2_calc(*_beam, *_elec);

  // _elec_mom = _elec->P();
  // _elec_E = _elec->E();
  // _theta_e = _elec->Theta() * 180 / PI;

  // ///////////////////////////// Smearing /////////////////////////////////
  //////////////////////////////// Smearing //////////////////////////////
  if (_mc) {
    _sectorElec = _data->dc_sec(0);
    _elec_status = abs(_data->status(0));

    _elecUnSmear->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);

    double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear,
        phiSmear;

    pUnSmear = _elecUnSmear->P();

    thetaUnSmear = _elecUnSmear->Theta() * 180 / PI;

    if (_elecUnSmear->Phi() > 0)
      phiUnSmear = _elecUnSmear->Phi() * 180 / PI;
    else if (_elecUnSmear->Phi() < 0)
      phiUnSmear = (_elecUnSmear->Phi() + 2 * PI) * 180 / PI;

    ////////////////////////////////////////////////////////////////

    // Generate new values
    Reaction::SmearingFunc(ELECTRON, _elec_status, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear);

    _pxPrimeSmear = _elecUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                    sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
    _pyPrimeSmear = _elecUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                    sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
    _pzPrimeSmear =
        _elecUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);

    // _elecSmear->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_E);  // smeared
    _elec->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_E);  // smeared
    // _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);  // unsmeared

    *_gamma += *_beam - *_elec;  // be careful you are commenting this only to include the momentum correction

    // // // // // Can calculate W and Q2 here (useful for simulations as sim do not have elec mom corrections)
    _W = physics::W_calc(*_beam, *_elec);
    _Q2 = physics::Q2_calc(*_beam, *_elec);

    // _elec_mom = _elec->P();
    // _elec_E = _elec->E();
    // _theta_e = _elec->Theta() * 180 / PI;

    // if (_elec->Phi() > 0)
    //   _phi_elec = _elec->Phi() * 180 / PI;
    // else if (_elec->Phi() < 0)
    //   _phi_elec = (_elec->Phi() + 2 * PI) * 180 / PI;
  }
}

// // ///////////////////////////// MOM CORR /////////////////////////////////
// //////////////////////////////// MOM CORR //////////////////////////////

void Reaction::SetMomCorrElec() {  // New electron momentum corrections
  if (!_mc) {
    fe = objMomCorr->dppC(_data->px(0), _data->py(0), _data->pz(0), _data->dc_sec(0), 0) + 1;
    _mom_corr_elec->SetXYZM(_data->px(0) * fe, _data->py(0) * fe, _data->pz(0) * fe,
                            MASS_E);  // this is new electron mom corrections aug 2022
    _elec->SetXYZM(_data->px(0) * fe, _data->py(0) * fe, _data->pz(0) * fe,
                   MASS_E);  // elec and mom corr elec are SAME !!!!!

    *_gamma += *_beam - *_mom_corr_elec;
    // _W_after = physics::W_calc(*_beam, *_mom_corr_elec);
    _W = physics::W_calc(*_beam, *_mom_corr_elec);
    _Q2 = physics::Q2_calc(*_beam, *_mom_corr_elec);

    _P_elec = _mom_corr_elec->P();
    _elec_E = _mom_corr_elec->E();
    _theta_e = _mom_corr_elec->Theta() * 180 / PI;
  }
}

void Reaction::SetProton(int i) {
  _numProt++;
  _numPos++;
  _hasP = true;
  _sectorProt = _data->dc_sec(i);
  _prot_status = abs(_data->status(i));
  _beta_prot = _data->beta(i);

  _Energy_loss_uncorr_prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
  // _prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);

  auto proton = std::make_unique<TLorentzVector>();
  auto mom_corr_proton = std::make_unique<TLorentzVector>();

  _prot_mom_uncorr = _Energy_loss_uncorr_prot->P();

  _prot_theta_uncorr = _Energy_loss_uncorr_prot->Theta() * 180 / PI;

  if (_Energy_loss_uncorr_prot->Phi() > 0)
    _prot_phi_uncorr = _Energy_loss_uncorr_prot->Phi() * 180 / PI;
  else if (_Energy_loss_uncorr_prot->Phi() < 0)
    _prot_phi_uncorr = (_Energy_loss_uncorr_prot->Phi() + 2 * PI) * 180 / PI;

  _is_FD_Prot = objMomCorr->is_FD(_prot_status);
  _is_CD_Prot = objMomCorr->is_CD(_prot_status);

  // _is_lower_band = mom_corr::is_lower_band(_prot_mom_uncorr, _thetaDC_r1_Prot, _prot_status);

  if (_is_CD_Prot) {
    _prot_mom_tmt = _prot_mom_uncorr;
    // _prot_theta_tmt = _prot_theta_uncorr;
    // _prot_phi_tmt = _prot_phi_uncorr;

    // _prot_mom_tmt = mom_corr::CD_prot_Emom_corr(_prot_mom_uncorr, _prot_theta_uncorr);
    // _prot_theta_tmt = mom_corr::CD_prot_Eth_corr(_prot_mom_uncorr, _prot_theta_uncorr);
    // _prot_phi_tmt = mom_corr::CD_prot_Eph_corr(_prot_mom_uncorr, _prot_theta_uncorr, _prot_phi_uncorr);
  }

  if (_is_FD_Prot) {
    if (_prot_theta_uncorr < 27) {
      _prot_mom_tmt = _prot_mom_uncorr + exp(-2.739 - 3.932 * _prot_theta_uncorr) + 0.002907;
      // if (_prot_mom_tmt < 2.4)
      //   _prot_mom_tmt = _prot_mom_uncorr + (0.000971) * pow(_prot_mom_uncorr, 4) +
      //                   (-0.009736) * pow(_prot_mom_uncorr, 3) + (0.034675) * pow(_prot_mom_uncorr, 2) +
      //                   (-0.052357) * _prot_mom_uncorr + 0.032840;
      // else
      //   _prot_mom_tmt = _prot_mom_uncorr + 0.004530;
    } else {
      _prot_mom_tmt = _prot_mom_uncorr + exp(-1.2 - 4.228 * _prot_mom_uncorr) + 0.007502;
      // if (_prot_mom_tmt < 2.4)
      //   _prot_mom_tmt = _prot_mom_uncorr + (0.002369) * pow(_prot_mom_uncorr, 4) +
      //                   (-0.023523) * pow(_prot_mom_uncorr, 3) + (0.082790) * pow(_prot_mom_uncorr, 2) +
      //                   (-0.124996) * _prot_mom_uncorr + 0.077659;
      // else
      //   _prot_mom_tmt = _prot_mom_uncorr + 0.007968;
    }
  }

  _px_prime_prot_E = _data->px(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
  _py_prime_prot_E = _data->py(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
  _pz_prime_prot_E = _data->pz(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));

  /// // _prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);  // energy loss corrected

  ////////////////////////// Mom CORR /////////////////////////////////
  if (!_mc) {
    // // Below shows how the corrections are to be applied using the ROOT momentum 4-vector using the above code:
    if (_is_FD_Prot) {
      // fpro = 1.0;
      fpro = objMomCorr->dppC(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, _data->dc_sec(i), 3) + 1;
    } else {
      fpro = 1.0;
    }
    proton->SetXYZM(_px_prime_prot_E * fpro, _py_prime_prot_E * fpro, _pz_prime_prot_E * fpro, MASS_P);
    mom_corr_proton->SetXYZM(_px_prime_prot_E * fpro, _py_prime_prot_E * fpro, _pz_prime_prot_E * fpro, MASS_P);

    _prot.push_back(std::move(proton));                    // Add proton to the vector
    _mom_corr_prot.push_back(std::move(mom_corr_proton));  // Add proton to the vector
    _prot_indices.push_back(i);                            // Store the index
  }
  /////////////////// SMEARING PART ////////////////////////////////////////////////////////////////////////////
  if (_mc) {
    _protUnSmear->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);  // energy loss corrected

    //////////////////////////////////////////////////////////////
    double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear,
        phiSmear;

    pUnSmear = _protUnSmear->P();

    thetaUnSmear = _protUnSmear->Theta() * 180 / PI;

    if (_protUnSmear->Phi() > 0)
      phiUnSmear = _protUnSmear->Phi() * 180 / PI;
    else if (_protUnSmear->Phi() < 0)
      phiUnSmear = (_protUnSmear->Phi() + 2 * PI) * 180 / PI;

    // Generate new values

    Reaction::SmearingFunc(PROTON, _prot_status, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear);

    _pxPrimeSmear = _protUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                    sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
    _pyPrimeSmear = _protUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                    sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
    _pzPrimeSmear =
        _protUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);

    // _protSmear->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_P);  // smeared

    proton->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_P);  // smeared

    _prot.push_back(std::move(proton));                    // Add proton to the vector
    _mom_corr_prot.push_back(std::move(mom_corr_proton));  // Add proton to the vector
    _prot_indices.push_back(i);                            // Store the index
  }
}
void Reaction::SetPip(int i) {
  _numPip++;
  _numPos++;
  _hasPip = true;
  _pip_status = abs(_data->status(i));
  _sectorPip = _data->dc_sec(i);
  _beta_pip = _data->beta(i);

  auto pip = std::make_unique<TLorentzVector>();
  auto mom_corr_pip = std::make_unique<TLorentzVector>();

  _Energy_loss_uncorr_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
  // // // _mom_corr_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
  // _pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);

  _pip_mom_uncorr = _Energy_loss_uncorr_pip->P();
  _pip_theta_uncorr = _Energy_loss_uncorr_pip->Theta() * 180 / PI;
  // // // if (_Energy_loss_uncorr_pip->Phi() > 0)
  // // //   _pip_phi_uncorr = _Energy_loss_uncorr_pip->Phi() * 180 / PI;
  // // // else if (_Energy_loss_uncorr_pip->Phi() < 0)
  // // //   _pip_phi_uncorr = (_Energy_loss_uncorr_pip->Phi() + 2 * PI) * 180 / PI;

  _is_FD_Pip = objMomCorr->is_FD(_pip_status);
  _is_CD_Pip = objMomCorr->is_CD(_pip_status);

  // // // _is_lower_band = mom_corr::is_lower_band(_pip_mom_uncorr, _thetaDC_r1_Pip, _pip_status);

  // if (_is_CD_Pip) {
  //   _pip_mom_tmt = _pip_mom_uncorr;
  //   _pip_mom_tmt = _pip_mom_uncorr + objMomCorr->elossPipCD(_pip_mom_uncorr, _pip_theta_uncorr);
  // }
  // if (_is_FD_Pip) {
  //   // _pip_mom_tmt = _pip_mom_uncorr;
  //   _pip_mom_tmt = _pip_mom_uncorr + objMomCorr->elossPipFD(_pip_mom_uncorr, _pip_theta_uncorr);
  // }
  ///// our version of elsoo pip
  // if (_is_CD_Pip) {
  //   _pip_mom_tmt = _pip_mom_uncorr;
  // }
  // if (_is_FD_Pip) {
  //   // _pim_mom_tmt = _pim_mom_uncorr;
  //   if (_pip_theta_uncorr < 27) {
  //     _pip_mom_tmt = _pip_mom_uncorr + 0.0002468543 * _pip_mom_uncorr + 0.00324120;
  //   } else {
  //     _pip_mom_tmt = _pip_mom_uncorr + -0.0004140691 * _pip_mom_uncorr + 0.007524105;
  //   }
  // }

  // eloss used by stefan
  if (_is_CD_Pip) {
    // _pip_mom_tmt = _pip_mom_uncorr;
    _pip_mom_tmt = _pip_mom_uncorr + objMomCorr->elossPipCD(_pip_mom_uncorr, _pip_theta_uncorr);
  }
  if (_is_FD_Pip) {
    // _pip_mom_tmt = _pip_mom_uncorr;
    _pip_mom_tmt = _pip_mom_uncorr + objMomCorr->elossPipFD(_pip_mom_uncorr, _pip_theta_uncorr);
  }
  _px_prime_pip_E = _data->px(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
  _py_prime_pip_E = _data->py(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
  _pz_prime_pip_E = _data->pz(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
  //// _pip->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);  // energy loss corrected

  ////////////////////////// Mom CORR /////////////////////////////////
  if (!_mc) {
    if (_is_FD_Pip) {
      // fpip = 1.0;
      fpip = objMomCorr->dppC(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, _data->dc_sec(i), 1) + 1;
    } else {
      fpip = 1.0;
    }
    //// unique ptr method for one pip
    // _pip->SetXYZM(_px_prime_pip_E * fpip, _py_prime_pip_E * fpip, _pz_prime_pip_E * fpip, MASS_PIP);
    // _mom_corr_pip->SetXYZM(_px_prime_pip_E * fpip, _py_prime_pip_E * fpip, _pz_prime_pip_E * fpip, MASS_PIP);

    /////// vector method for many pip
    pip->SetXYZM(_px_prime_pip_E * fpip, _py_prime_pip_E * fpip, _pz_prime_pip_E * fpip, MASS_PIP);
    mom_corr_pip->SetXYZM(_px_prime_pip_E * fpip, _py_prime_pip_E * fpip, _pz_prime_pip_E * fpip, MASS_PIP);

    _pip.push_back(std::move(pip));                    // Add pip to the vector
    _mom_corr_pip.push_back(std::move(mom_corr_pip));  // Add pip to the vector
    _pip_indices.push_back(i);
  }  // Store the index

  // /////////////////////////////////     SMEARING PART  /////////////////////////////
  if (_mc) {
    _pipUnSmear->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);  // energy loss corrected

    double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear,
        phiSmear;

    pUnSmear = _pipUnSmear->P();

    thetaUnSmear = _pipUnSmear->Theta() * 180 / PI;

    if (_pipUnSmear->Phi() > 0)
      phiUnSmear = _pipUnSmear->Phi() * 180 / PI;
    else if (_pipUnSmear->Phi() < 0)
      phiUnSmear = (_pipUnSmear->Phi() + 2 * PI) * 180 / PI;

    // Generate new values
    Reaction::SmearingFunc(PIP, _pip_status, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear);

    _pxPrimeSmear = _pipUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                    sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
    _pyPrimeSmear = _pipUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                    sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
    _pzPrimeSmear =
        _pipUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);

    // _pipSmear->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIP);  // smeared
    pip->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIP);  // smeared
    _pip.push_back(std::move(pip));
    _pip_indices.push_back(i);  // Store the index
  }
}

void Reaction::SetPim(int i) {
  _numPim++;
  _numNeg++;
  _hasPim = true;
  _pim_status = abs(_data->status(i));
  _sectorPim = _data->dc_sec(i);

  auto pim = std::make_unique<TLorentzVector>();
  auto mom_corr_pim = std::make_unique<TLorentzVector>();

  _Energy_loss_uncorr_pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
  // _pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);

  _pim_mom_uncorr = _Energy_loss_uncorr_pim->P();
  _pim_theta_uncorr = _Energy_loss_uncorr_pim->Theta() * 180 / PI;
  // // // // if (_Energy_loss_uncorr_pim->Phi() > 0)
  // // // //   _pim_phi_uncorr = _Energy_loss_uncorr_pim->Phi() * 180 / PI;
  // // // // else if (_Energy_loss_uncorr_pim->Phi() < 0)
  // // // //   _pim_phi_uncorr = (_Energy_loss_uncorr_pim->Phi() + 2 * PI) * 180 / PI;

  _is_FD_Pim = objMomCorr->is_FD(_pim_status);
  _is_CD_Pim = objMomCorr->is_CD(_pim_status);

  // // // _is_lower_band = mom_corr::is_lower_band(_pim_mom_uncorr, _thetaDC_r1_Pim, _pim_status);
  if (_is_CD_Pim) {
    _pim_mom_tmt = _pim_mom_uncorr;
  }
  if (_is_FD_Pim) {
    // _pim_mom_tmt = _pim_mom_uncorr;
    if (_pim_theta_uncorr < 27) {
      _pim_mom_tmt = _pim_mom_uncorr + 0.00046571 * _pim_mom_uncorr + 0.00322164;
    } else {
      if (_pim_mom_uncorr < 1.7)
        _pim_mom_tmt = _pim_mom_uncorr + (-0.0024313) * pow(_pim_mom_uncorr, 3) +
                       (0.0094416) * pow(_pim_mom_uncorr, 2) + (-0.01257967) * pow(_pim_mom_uncorr, 1) + 0.0122432;
      else
        _pim_mom_tmt = _pim_mom_uncorr + 0.006199071;
    }
  }
  _px_prime_pim_E = _data->px(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
  _py_prime_pim_E = _data->py(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
  _pz_prime_pim_E = _data->pz(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));

  ////////////////////////// Mom CORR /////////////////////////////////
  if (!_mc) {
    if (_is_FD_Pim) {
      // fpim = 1.0;
      fpim = objMomCorr->dppC(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, _data->dc_sec(i), 2) + 1;
    } else {
      fpim = 1.0;
    }
    // _pim->SetXYZM(_px_prime_pim_E * fpim, _py_prime_pim_E * fpim, _pz_prime_pim_E * fpim, MASS_PIM);
    // _mom_corr_pim->SetXYZM(_px_prime_pim_E * fpim, _py_prime_pim_E * fpim, _pz_prime_pim_E * fpim, MASS_PIM);

    /////// vector method for many pim
    pim->SetXYZM(_px_prime_pim_E * fpim, _py_prime_pim_E * fpim, _pz_prime_pim_E * fpim, MASS_PIM);
    mom_corr_pim->SetXYZM(_px_prime_pim_E * fpim, _py_prime_pim_E * fpim, _pz_prime_pim_E * fpim, MASS_PIM);

    _pim.push_back(std::move(pim));                    // Add pim to the vector
    _mom_corr_pim.push_back(std::move(mom_corr_pim));  // Add pim to the vector
    _pim_indices.push_back(i);
  }
  ////// _pim->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, MASS_PIM);  // energy loss corrected
  if (_mc) {
    // /////////////////////////////////     SMEARING PART  /////////////////////////////

    _pimUnSmear->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, MASS_PIM);

    double _pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear,
        phiSmear;

    pUnSmear = _pimUnSmear->P();

    thetaUnSmear = _pimUnSmear->Theta() * 180 / PI;

    if (_pimUnSmear->Phi() > 0)
      phiUnSmear = _pimUnSmear->Phi() * 180 / PI;
    else if (_pimUnSmear->Phi() < 0)
      phiUnSmear = (_pimUnSmear->Phi() + 2 * PI) * 180 / PI;

    // Generate new values
    Reaction::SmearingFunc(PIM, _pim_status, pUnSmear, thetaUnSmear, phiUnSmear, pSmear, thetaSmear, phiSmear);

    _pxPrimeSmear = _pimUnSmear->Px() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                    sin(DEG2RAD * thetaUnSmear) * cos(DEG2RAD * phiSmear) / cos(DEG2RAD * phiUnSmear);
    _pyPrimeSmear = _pimUnSmear->Py() * ((pSmear) / (pUnSmear)) * sin(DEG2RAD * thetaSmear) /
                    sin(DEG2RAD * thetaUnSmear) * sin(DEG2RAD * phiSmear) / sin(DEG2RAD * phiUnSmear);
    _pzPrimeSmear =
        _pimUnSmear->Pz() * ((pSmear) / (pUnSmear)) * cos(DEG2RAD * thetaSmear) / cos(DEG2RAD * thetaUnSmear);

    // _pimSmear->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIM);  // smeared
    pim->SetXYZM(_pxPrimeSmear, _pyPrimeSmear, _pzPrimeSmear, MASS_PIM);  // smeared

    _pim.push_back(std::move(pim));  // Add pim to the vector
    _pim_indices.push_back(i);       // Store the index
  }
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

// /////////////////// new added ////////////////
// void Reaction::CalcMissMassPim(const TLorentzVector &prot, const TLorentzVector &pip) {
//   auto mm_mpim = std::make_unique<TLorentzVector>();

//   *mm_mpim += (*_gamma + *_target);
//   *mm_mpim -= prot;
//   *mm_mpim -= pip;

//   _MM_mPim = mm_mpim->M();
//   _MM2_mPim = mm_mpim->M2();
// }

void Reaction::CalcMissMassExcl(const TLorentzVector &prot, const TLorentzVector &pip, const TLorentzVector &pim)
//     void Reaction::CalcMissMass()
{
  auto mm_mpim = std::make_unique<TLorentzVector>();
  auto mm_mpip = std::make_unique<TLorentzVector>();
  auto mm_mprot = std::make_unique<TLorentzVector>();
  auto mm_excl = std::make_unique<TLorentzVector>();

  if (TwoPion_exclusive()) {
    *mm_excl += (*_gamma + *_target);
    *mm_excl -= prot;
    *mm_excl -= pip;
    *mm_excl -= pim;
    _MM_exclusive = mm_excl->M();
    _MM2_exclusive = mm_excl->M2();
    _excl_Energy = mm_excl->E();
    // }

    // //         if (TwoPion_missingPim())
    // // {
    *mm_mpim += (*_gamma + *_target);
    *mm_mpim -= prot;
    *mm_mpim -= pip;

    _MM_mPim = mm_mpim->M();
    _MM2_mPim = mm_mpim->M2();
    // // }
    // // if (TwoPion_missingPip())
    // // {
    *mm_mpip += (*_gamma + *_target);
    *mm_mpip -= prot;
    *mm_mpip -= pim;
    _MM2_mPip = mm_mpip->M2();
    // }
    // if (TwoPion_missingProt())
    // {
    *mm_mprot += (*_gamma + *_target);
    *mm_mprot -= pip;
    *mm_mprot -= pim;
    _MM2_mProt = mm_mprot->M2();
    // }
  }
}

float Reaction::MM_mPim() {
  // if (_MM_mPim != _MM_mPim)
  //         CalcMissMass(*_prot[0], *_pip[0]); // This is just a default case for first proton/pion pair
  return _MM_mPim;
}

float Reaction::MM2_mPim() {
  // if (_MM2_mPim != _MM2_mPim)
  //         CalcMissMass(*_prot[0], *_pip[0]); // This is just a default case for first proton/pion pair
  return _MM2_mPim;
}
float Reaction::MM2_exclusive() {
  // if (_MM2_exclusive != _MM2_exclusive)
  // CalcMissMass();
  return _MM2_exclusive;
}
float Reaction::MM_exclusive() {
  // if (_MM_exclusive != _MM_exclusive)
  //         CalcMissMass();
  return _MM_exclusive;
}
float Reaction::MM2_mPip() {
  // if (_MM2_mPip != _MM2_mPip)
  //         CalcMissMass();
  return _MM2_mPip;
}
float Reaction::MM2_mProt() {
  // if (_MM2_mProt != _MM2_mProt)
  //         CalcMissMass();
  return _MM2_mProt;
}
float Reaction::Energy_excl() {
  // if (_excl_Energy != _excl_Energy)
  //         CalcMissMass();
  //  std::cout << "_x_mu_p  " << _x_mu->E() << '\n';
  //  if (_x_mu_E > 0)
  return _excl_Energy;
  // else
  // return NAN;
}

float Reaction::prot_momentum(const TLorentzVector &prot) {
  // if (TwoPion_missingPim())
  if (_hasP) return prot.P();
  // return _prot->P();
  else
    return NAN;
}

float Reaction::pip_momentum(const TLorentzVector &pip) {
  // if (TwoPion_missingPim())
  if (_hasPip) return pip.P();
  // return _pip->P();
  else
    return NAN;
}
float Reaction::pim_momentum(const TLorentzVector &prot, const TLorentzVector &pip) {
  if (TwoPion_missingPim()) {
    auto missingpim_ = std::make_unique<TLorentzVector>();
    // *missingpim_ += *_gamma + *_target - *_prot - *_pip;
    *missingpim_ += *_gamma + *_target - prot - pip;
    return missingpim_->P();
  } else
    return NAN;
}

float Reaction::pim_E(const TLorentzVector &prot, const TLorentzVector &pip) {
  if (TwoPion_missingPim()) {
    auto missingpim_ = std::make_unique<TLorentzVector>();
    // *missingpim_ += *_gamma + *_target - *_prot - *_pip;
    *missingpim_ += *_gamma + *_target - prot - pip;
    return missingpim_->E();
  } else
    return NAN;
}
float Reaction::pim_momentum_measured(const TLorentzVector &pim) {
  if (TwoPion_exclusive())
    return pim.P();
  else
    return NAN;
}
float Reaction::pim_E_measured(const TLorentzVector &pim) {
  if (TwoPion_exclusive())
    return pim.E();
  else
    return NAN;
}

// float Reaction::elec_momentum() {
//   if (TwoPion_missingPim())
//     return _elec->P();
//   else
//     return NAN;
// }
// float Reaction::theta_elec() {  /// lab theta mattrai hunchha electron ko case ma
//   // if (TwoPion_missingPim())
//   if (_elec->Theta() > -500)
//     return _elec->Theta() * 180.0 / PI;
//   else
//     return NAN;
// }
// float Reaction::Phi_elec() {  /// lab theta mattrai hunchha electron ko case ma
//   // if (TwoPion_missingPim()) {
//   //         if ((_elec->Phi() * 180.0 / PI) > -150)
//   //                 return _elec->Phi() * 180.0 / PI;
//   //         else
//   //                 return ((_elec->Phi() * 180.0 / PI) + 360);
//   // }
//   // else
//   //         return NAN;

//   if (_elec->Phi() > -500) {
//     if (_elec->Phi() > 0)
//       return _elec->Phi() * 180 / PI;
//     else if (_elec->Phi() < 0)
//       return (_elec->Phi() + 2 * PI) * 180 / PI;
//     else
//       return NAN;
//   } else
//     return NAN;
// }
float Reaction::prot_theta_lab(const TLorentzVector &prot) {
  // if (TwoPion_missingPim())
  if (_hasP)
    return prot.Theta() * 180.0 / PI;
  else
    return NAN;
}
float Reaction::pip_theta_lab(const TLorentzVector &pip) {
  // if (TwoPion_missingPim())
  if (_hasPip)
    return pip.Theta() * 180.0 / PI;
  else
    return NAN;
}
float Reaction::pim_theta_lab(const TLorentzVector &prot, const TLorentzVector &pip) {
  if (TwoPion_missingPim()) {
    auto missingpim_ = std::make_unique<TLorentzVector>();
    *missingpim_ += *_gamma + *_target - prot - pip;
    return missingpim_->Theta() * 180.0 / PI;
  } else
    return NAN;
}
float Reaction::pim_theta_lab_measured(const TLorentzVector &pim) {  /////////////////////////////////////work here

  if (TwoPion_exclusive())
    return pim.Theta() * 180.0 / PI;
  else
    return NAN;
}
float Reaction::prot_Phi_lab(const TLorentzVector &prot) {
  // if (TwoPion_missingPim())
  // {
  if (_hasP) {
    if (prot.Phi() > 0)
      return prot.Phi() * 180 / PI;
    else if (prot.Phi() < 0)
      return (prot.Phi() + 2 * PI) * 180 / PI;
    // else return NAN;
  } else
    return NAN;
}
float Reaction::prot_momT(const TLorentzVector &prot) {
  // if (TwoPion_missingPim())
  if (_hasP)
    return prot.Perp();
  else
    return NAN;
}
float Reaction::pip_momT(const TLorentzVector &pip) {
  // if (TwoPion_missingPim())
  if (_hasPip)
    return pip.Perp();
  else
    return NAN;
}
float Reaction::pim_momT(const TLorentzVector &pim) {
  // if (TwoPion_missingPim())
  return pim.Perp();
  // else
  //         return NAN;
}
float Reaction::pip_Phi_lab(const TLorentzVector &pip) {
  if (_hasPip) {
    if (pip.Phi() > 0)
      return pip.Phi() * 180 / PI;
    else if (pip.Phi() <= 0)
      return (pip.Phi() + 2 * PI) * 180 / PI;
  } else
    return NAN;
}
float Reaction::pim_Phi_lab(const TLorentzVector &prot, const TLorentzVector &pip) {
  auto missingpim_ = std::make_unique<TLorentzVector>();
  *missingpim_ += *_gamma + *_target - prot - pip;
  if (TwoPion_missingPim()) {
    if (missingpim_->Phi() > 0)
      return missingpim_->Phi() * 180 / PI;
    else if (missingpim_->Phi() <= 0)
      return (missingpim_->Phi() + 2 * PI) * 180 / PI;
  } else
    return NAN;
}

float Reaction::pim_Phi_lab_measured(const TLorentzVector &pim) {
  if (TwoPion_exclusive()) {
    if (pim.Phi() > 0)
      return pim.Phi() * 180 / PI;
    else if (pim.Phi() <= 0)
      return (pim.Phi() + 2 * PI) * 180 / PI;
  } else
    return NAN;
}

// // Boost the particles to the center-of-mass system
// // void Reaction::boost(const TLorentzVector &prot, const TLorentzVector &pip, const TLorentzVector &pim)
void Reaction::boost(const TLorentzVector &prot, const TLorentzVector &pip) {
  _is_boosted = true;
  // std::unique_ptr<TLorentzVector> pim = std::make_unique<TLorentzVector>(*_gamma + *_target - prot - pip);
  const TLorentzVector pim = *_gamma + *_target - prot - pip;

  // Boost all particles to the center of mass frame
  _boosted_gamma = std::make_unique<TLorentzVector>(boost_cms::boostToCMS(*_gamma, *_gamma, *_elec, _Q2));
  _boosted_prot = std::make_unique<TLorentzVector>(boost_cms::boostToCMS(prot, *_gamma, *_elec, _Q2));
  _boosted_pip = std::make_unique<TLorentzVector>(boost_cms::boostToCMS(pip, *_gamma, *_elec, _Q2));
  _boosted_pim = std::make_unique<TLorentzVector>(boost_cms::boostToCMS(pim, *_gamma, *_elec, _Q2));
}
// // // Calculate invariant masse
float Reaction::inv_Ppip() { return boost_cms::calculateInvariantMass(*_boosted_prot, *_boosted_pip); }
float Reaction::inv_Ppim() { return boost_cms::calculateInvariantMass(*_boosted_prot, *_boosted_pim); }
float Reaction::inv_pip_pim() { return boost_cms::calculateInvariantMass(*_boosted_pip, *_boosted_pim); }

// float Reaction::w_P2pi_rec()
// {
//         if (_W_P2pi != _W_P2pi)
//                 W_2pi_P();
//         return _W_P2pi;
// }

// //////////////
float Reaction::prot_theta() { return boost_cms::calculateTheta(*_boosted_prot); }
float Reaction::pip_theta() { return boost_cms::calculateTheta(*_boosted_pip); }
float Reaction::pim_theta() { return boost_cms::calculateTheta(*_boosted_pim); }

// //////////////
float Reaction::gamma_Phi() { return boost_cms::calculatePhi(*_boosted_gamma); }

float Reaction::prot_Phi() { return boost_cms::calculatePhi(*_boosted_prot); }

float Reaction::pip_Phi() { return boost_cms::calculatePhi(*_boosted_pip); }
float Reaction::pim_Phi() { return boost_cms::calculatePhi(*_boosted_pim); }

// //////////////
float Reaction::alpha_ppip_pipim() {
  return boost_cms::calculateAlpha(_boosted_pim->Vect().Unit(), _boosted_pip->Vect().Unit(), _boosted_pim->Vect());
}

float Reaction::alpha_pippim_pipf() {
  return boost_cms::calculateAlpha(_boosted_prot->Vect().Unit(), _boosted_pip->Vect().Unit(), _boosted_prot->Vect());
}

float Reaction::alpha_ppim_pipip() {
  return boost_cms::calculateAlpha(_boosted_pip->Vect().Unit(), _boosted_pim->Vect().Unit(), _boosted_pip->Vect());
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

// float_t Reaction::scalar_triple_product() {
//   if (!_is_boosted) boost();
//   if (TwoPion_exclusive()) {
//     return (_prot_Vect3.Dot(_pip_Vect3.Cross(_pim_Vect3)));

//   } else
//     return NAN;
// }

// ////////////
// float Reaction::pim_momentum_cm() {
//   if (!_is_boosted) boost();
//   if (TwoPion_exclusive())
//     return _boosted_pim->P();
//   else
//     return NAN;
// }

// float Reaction::pim_theta_cm() {
//   if (!_is_boosted) boost();
//   if (TwoPion_exclusive())
//     return _rotated_pim->Theta() * 180.0 / PI;
//   else
//     return NAN;
// }

// float Reaction::pim_Phi_cm() {
//   if (!_is_boosted) boost();
//   if (TwoPion_exclusive()) {
//     if (_rotated_pim->Phi() > 0)
//       return _rotated_pim->Phi() * 180 / PI;
//     else if (_rotated_pim->Phi() < 0)
//       return (_rotated_pim->Phi() + 2 * PI) * 180 / PI;
//     else
//       return NAN;
//   } else
//     return NAN;
// }

// float Reaction::pim_momentum_cm_measured() {
//   if (!_is_boosted) boost();
//   if (TwoPion_exclusive())
//     return _boosted_pim_measured->P();
//   else
//     return NAN;
// }

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

MCReaction::MCReaction(const std::shared_ptr<Branches12> &data, float beam_enrgy) {
  _data = data;
  if (!_data->mc()) _data->mc_branches();
  _beam = std::make_unique<TLorentzVector>();
  _beam_energy = beam_enrgy;
  _weight_mc = _data->mc_weight();
  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  //_gamma = std::make_unique<TLorentzVector>();  // do i need this?
  _gamma_mc = std::make_unique<TLorentzVector>();
  _target_mc = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  //_elec = std::make_unique<TLorentzVector>();  // do i need this?
  _elec_mc = std::make_unique<TLorentzVector>();
  // this->SetElec();  // do i need this?
  this->SetMCElec();
  // _prot_mc = std::make_unique<TLorentzVector>();
  // _pip_mc = std::make_unique<TLorentzVector>();
  // _pim_mc = std::make_unique<TLorentzVector>();
  //_other = std::make_unique<TLorentzVector>();  // do i need this?
  _other_mc = std::make_unique<TLorentzVector>();
  //_neutron = std::make_unique<TLorentzVector>();

  _prot_mc = std::vector<std::unique_ptr<TLorentzVector>>();  // Initialize _protons as an empty vector
  _pip_mc = std::vector<std::unique_ptr<TLorentzVector>>();
  _pim_mc = std::vector<std::unique_ptr<TLorentzVector>>();

  _prot_mc_indices = std::vector<int>();
  _pip_mc_indices = std::vector<int>();
  _pim_mc_indices = std::vector<int>();
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

void MCReaction::SetMCProton(int i) {
  // _prot_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_P);
  auto proton_mc = std::make_unique<TLorentzVector>();
  proton_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_P);
  _prot_mc.push_back(std::move(proton_mc));
  _prot_mc_indices.push_back(i);
}

void MCReaction::SetMCPip(int i) {
  // _pip_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIP);
  auto pionp_mc = std::make_unique<TLorentzVector>();
  pionp_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIP);
  _pip_mc.push_back(std::move(pionp_mc));
  _pip_mc_indices.push_back(i);
}

void MCReaction::SetMCPim(int i) {
  auto pionm_mc = std::make_unique<TLorentzVector>();
  pionm_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIM);
  _pim_mc.push_back(std::move(pionm_mc));
  _pim_mc_indices.push_back(i);
}

/////////////////// new added ////////////////
void MCReaction::CalcMissMassPimMC(const TLorentzVector &prot_mc, const TLorentzVector &pip_mc) {
  // void MCReaction::CalcMissMassPimMC() {
  auto mm_mpim_mc = std::make_unique<TLorentzVector>();

  *mm_mpim_mc += (*_gamma_mc + *_target_mc);
  *mm_mpim_mc -= prot_mc;
  *mm_mpim_mc -= pip_mc;
  // *mm_mpim_mc -= *_prot_mc;
  // *mm_mpim_mc -= *_pip_mc;

  _MM_mPim_mc = mm_mpim_mc->M();
  _MM2_mPim_mc = mm_mpim_mc->M2();
}

float MCReaction::MM_mPim_MC() { return _MM_mPim_mc; }

float MCReaction::MM2_mPim_MC() { return _MM2_mPim_mc; }

// float MCReaction::pim_mom_mc_gen() { return _pim_mc->P(); }
// float MCReaction::pip_mom_mc_gen() { return _pip_mc->P(); }
// float MCReaction::prot_mom_mc_gen() { return _prot_mc->P(); }

// float MCReaction::pim_theta_mc_gen() { return _pim_mc->Theta() * 180 / PI; }
// float MCReaction::pip_theta_mc_gen() { return _pip_mc->Theta() * 180 / PI; }
// float MCReaction::prot_theta_mc_gen() { return _prot_mc->Theta() * 180 / PI; }

// float MCReaction::pim_phi_mc_gen() {
//   if (_pim_mc->Phi() >= 0)
//     return (_pim_mc->Phi() * 180 / PI);
//   else if (_pim_mc->Phi() < 0)
//     return ((_pim_mc->Phi() + 2 * PI) * 180 / PI);
//   else
//     return NAN;
// }
// float MCReaction::pip_phi_mc_gen() {
//   if (_pip_mc->Phi() >= 0)
//     return (_pip_mc->Phi() * 180 / PI);
//   else if (_pip_mc->Phi() < 0)
//     return ((_pip_mc->Phi() + 2 * PI) * 180 / PI);
//   else
//     return NAN;
// }
// float MCReaction::prot_phi_mc_gen() {
//   if (_prot_mc->Phi() >= 0)
//     return (_prot_mc->Phi() * 180 / PI);
//   else if (_prot_mc->Phi() < 0)
//     return ((_prot_mc->Phi() + 2 * PI) * 180 / PI);
//   else
//     return NAN;
// }
// ////////////////////////////////  BOOST TO CM SYSTEM ///////////////////////////////

// void MCReaction::boost_mc(const TLorentzVector &prot_mc, const TLorentzVector &pip_mc, const TLorentzVector &pim_mc)
// {
//   _is_boosted_mc = true;

//   // Boost all particles to the center of mass frame
//   _boosted_gamma_mc =
//       std::make_unique<TLorentzVector>(boost_cms::boostToCMS(*_gamma_mc, *_gamma_mc, *_elec_mc, _Q2_mc));
//   _boosted_prot_mc = std::make_unique<TLorentzVector>(boost_cms::boostToCMS(prot_mc, *_gamma_mc, *_elec_mc, _Q2_mc));
//   _boosted_pip_mc = std::make_unique<TLorentzVector>(boost_cms::boostToCMS(pip_mc, *_gamma_mc, *_elec_mc, _Q2_mc));
//   _boosted_pim_mc = std::make_unique<TLorentzVector>(boost_cms::boostToCMS(pim_mc, *_gamma_mc, *_elec_mc, _Q2_mc));
// }
// // // // Calculate invariant masse
// float MCReaction::MCinv_Ppip() { return boost_cms::calculateInvariantMass(*_boosted_prot_mc, *_boosted_pip_mc); }
// float MCReaction::MCinv_Ppim() { return boost_cms::calculateInvariantMass(*_boosted_prot_mc, *_boosted_pim_mc); }
// float MCReaction::MCinv_pip_pim() { return boost_cms::calculateInvariantMass(*_boosted_pip_mc, *_boosted_pim_mc); }

// // //////////////
// float MCReaction::MCprot_theta_thrown() { return boost_cms::calculateTheta(*_boosted_prot_mc); }
// float MCReaction::MCpip_theta_thrown() { return boost_cms::calculateTheta(*_boosted_pip_mc); }
// float MCReaction::MCpim_theta_thrown() { return boost_cms::calculateTheta(*_boosted_pim_mc); }

// // //////////////
// float MCReaction::MCgamma_Phi_thrown() { return boost_cms::calculatePhi(*_boosted_gamma_mc); }

// float MCReaction::MCprot_Phi_thrown() { return boost_cms::calculatePhi(*_boosted_prot_mc); }

// float MCReaction::MCpip_Phi_thrown() { return boost_cms::calculatePhi(*_boosted_pip_mc); }
// float MCReaction::MCpim_Phi_thrown() { return boost_cms::calculatePhi(*_boosted_pim_mc); }

// // //////////////
// float MCReaction::MCalpha_ppip_pipim_thrown() {
//   return boost_cms::calculateAlpha(_boosted_pim_mc->Vect().Unit(), _boosted_pip_mc->Vect().Unit(),
//                                    _boosted_pim_mc->Vect());
// }

// float MCReaction::MCalpha_pippim_pipf_thrown() {
//   return boost_cms::calculateAlpha(_boosted_prot_mc->Vect().Unit(), _boosted_pip_mc->Vect().Unit(),
//                                    _boosted_prot_mc->Vect());
// }

// float MCReaction::MCalpha_ppim_pipip_thrown() {
//   return boost_cms::calculateAlpha(_boosted_pip_mc->Vect().Unit(), _boosted_pim_mc->Vect().Unit(),
//                                    _boosted_pip_mc->Vect());
// }

///////////////////////////////////////////////////////////////////////////
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
