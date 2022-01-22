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
  this->SetElec();
  _prot = std::make_unique<TLorentzVector>();
  _pip = std::make_unique<TLorentzVector>();
  _pim = std::make_unique<TLorentzVector>();
  _other = std::make_unique<TLorentzVector>();
  _neutron = std::make_unique<TLorentzVector>();
}

Reaction::~Reaction() {}

void Reaction::SetElec() {
  _hasE = true;
  _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);

  *_gamma += *_beam - *_elec;

  // Can calculate W and Q2 here
  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);
}

void Reaction::SetProton(int i) {
  _numProt++;
  _numPos++;
  _hasP = true;
  _prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
}
void Reaction::SetPip(int i) {
  _numPip++;
  _numPos++;
  _hasPip = true;
  _pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
}
void Reaction::SetPim(int i) {
  _numPim++;
  _numNeg++;
  _hasPim = true;
  _pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
}

// float Reaction::rec_pim_px() {
//   return _beam->Px() - _elec->Px() + _target->Px() - _pip->Px() - _prot->Px() - _pim->Px();
// }
// float Reaction::rec_pim_py() {
//   return _beam->Py() - _elec->Py() + _target->Py() - _pip->Py() - _prot->Py() - _pim->Py();
// }
// float Reaction::rec_pim_pz() {
//   return _beam->Pz() - _elec->Pz() + _target->Pz() - _pip->Pz() - _prot->Pz() - _pim->Pz();
// }
// float Reaction::rec_pim_E() { return _beam->E() - _elec->E() + _target->E() - _pip->E() - _prot->E() - _pim->E(); }
// float Reaction::rec_pim_P() {
//   return sqrt(abs(pow((_beam->Px() - _elec->Px() + _target->Px() - _pip->Px() - _prot->Px() - _pim->Px()), 2) +
//                   pow((_beam->Py() - _elec->Py() + _target->Py() - _pip->Py() - _prot->Py() - _pim->Py()), 2) +
//                   pow((_beam->Pz() - _elec->Pz() + _target->Pz() - _pip->Pz() - _prot->Pz() - _pim->Pz()), 2)));
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

  // if (TwoPion_missingPim()) {
  //   *mm -= *_prot;
  //   *mm -= *_pip;
  //   *mm -= *_pim;
  //   _MM = mm->E();  /// just for test printing energy
  //   _MM2 = mm->M2();

  // _rec_pim_mom = mm->P();
  // _rec_pim_theta = mm->Theta() * 180 / PI;

  // if (mm->Phi() >= 0)
  //   _rec_pim_phi = (mm->Phi() * 180 / PI);
  // else if (mm->Phi() < 0)
  //   _rec_pim_phi = ((mm->Phi() + 2 * PI) * 180 / PI);

  // _x_mu_E = mm->E();
  // _x_mu_P = mm->P();
  // _x_mu_Px = mm->Px();
  // _x_mu_Py = mm->Py();
  // _x_mu_Pz = mm->Pz();
  // _x_mu_theta = mm->Theta() * RAD2DEG;
  // _x_mu_m2 = mm->E() * mm->E() - mm->P() * mm->P();
  // _x_mu_m = mm->E() - mm->P();
  //   //
  // }
  if (TwoPion_exclusive()) {
    *mm -= *_prot;
    *mm -= *_pip;
    // *mm -= *_pim;
    _MM = mm->M();
    _MM2 = mm->M2();

    *mm_excl += (*_gamma + *_target);
    *mm_excl -= *_prot;
    *mm_excl -= *_pip;
    *mm_excl -= *_pim;


    _MM2_exclusive = mm_excl->M2();
    _excl_Energy = mm_excl->E();

    _rec_pim_mom = mm->P();
    _rec_pim_theta = mm->Theta() * 180 / PI;

    if (mm->Phi() >= 0)
      _rec_pim_phi = (mm->Phi() * 180 / PI);
    else if (mm->Phi() < 0)
      _rec_pim_phi = ((mm->Phi() + 2 * PI) * 180 / PI);

    // if (mm_excl->Phi() >= 0)
    //   _x_mu_phi = (mm_excl->Phi() * 180 / PI);
    // else if (mm_excl->Phi() < 0)
    //   _x_mu_phi = ((mm_excl->Phi() + 2 * PI) * 180 / PI);

    // if (_elec->Phi() >= 0)
    //   _elec_phi = (_elec->Phi() * 180 / PI);
    // else if (_elec->Phi() < 0)
    //   _elec_phi = ((_elec->Phi() + 2 * PI) * 180 / PI);

    // if (_beam->Phi() >= 0)
    //   _beam_phi = (_beam->Phi() * 180 / PI);
    // else if (_beam->Phi() < 0)
    //   _beam_phi = ((_beam->Phi() + 2 * PI) * 180 / PI);

    // _diff_elec_x_mu_theta = (_elec->Theta() * 180 / PI) - (mm_excl->Theta() * 180 / PI);
    // _diff_elec_x_mu_phi = (_elec_phi - _x_mu_phi);

    // _diff_beam_x_mu_theta = (mm_excl->Theta() * 180 / PI);
    // _diff_beam_x_mu_phi = (_beam_phi - _x_mu_phi);

    // std::cout << " beam_theta " << _diff_beam_x_mu_theta << std::endl;
    // std::cout << " rec_pim_energy " << mm->E() << std::endl;

    // for mPip peak with exclusive events
    *mm_mpip += (*_gamma + *_target);
    *mm_mpip -= *_prot;
    *mm_mpip -= *_pim;
    _MM2_mPip = mm_mpip->M2();

    // for mProt peak with exclusive events
    *mm_mprot += (*_gamma + *_target);
    *mm_mprot -= *_pip;
    *mm_mprot -= *_pim;
    _MM2_mProt = mm_mprot->M2();
  }
  //   if (TwoPion_missingPip()) {
  //     *mm_mpip += (*_gamma + *_target);
  //     *mm_mpip -= *_prot;
  //     *mm_mpip -= *_pim;
  //     _MM2_mPip = mm_mpip->M2();
  //   }
  //   if (TwoPion_missingProt()) {
  //     *mm_mprot += (*_gamma + *_target);
  //     *mm_mprot -= *_pip;
  //     *mm_mprot -= *_pim;
  //     _MM2_mProt = mm_mprot->M2();
  //   }
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
  if (_rec_pim_mom != _rec_pim_mom) CalcMissMass();

  // if (TwoPion_missingPim()) {
  // auto missingpim_ = std::make_unique<TLorentzVector>();
  // *missingpim_ += *_gamma + *_target - *_prot - *_pip;
  // return missingpim_->P();
  return _rec_pim_mom;

  // } else
  //   return NAN;
}
float Reaction::pim_theta_lab() {
  if (_rec_pim_theta != _rec_pim_theta) CalcMissMass();

  // if (TwoPion_missingPim()) {
  // auto missingpim_ = std::make_unique<TLorentzVector>();
  // *missingpim_ += *_gamma + *_target - *_prot - *_pip;
  // return missingpim_->Theta() * 180.0 / PI;
  return _rec_pim_theta;
  // } else
  //   return NAN;
}
float Reaction::pim_Phi_lab() {
  if (_rec_pim_phi != _rec_pim_phi) CalcMissMass();

  // if (TwoPion_missingPim()) {
  // auto missingpim_ = std::make_unique<TLorentzVector>();
  // *missingpim_ += *_gamma + *_target - *_prot - *_pip;
  // if (missingpim_->Phi() > 0)
  //   return missingpim_->Phi() * 180 / PI;
  // else if (missingpim_->Phi() < 0)
  //   return (missingpim_->Phi() + 2 * PI) * 180 / PI;
  // else
  //   return NAN;
  return _rec_pim_phi;
  // } else
  //   return NAN;
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
////////////////mPip
float Reaction::pip_momentum() {
  if (TwoPion_missingPip()) {
    auto missingpip_ = std::make_unique<TLorentzVector>();
    *missingpip_ += *_gamma + *_target - *_prot - *_pim;
    return missingpip_->P();
  } else
    return NAN;
}
float Reaction::pip_theta_lab() {
  if (TwoPion_missingPip()) {
    auto missingpip_ = std::make_unique<TLorentzVector>();
    *missingpip_ += *_gamma + *_target - *_prot - *_pim;
    return missingpip_->Theta() * 180.0 / PI;
  } else
    return NAN;
}
float Reaction::pip_Phi_lab() {
  if (TwoPion_missingPip()) {
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
  if (TwoPion_missingProt()) {
    auto missingprot_ = std::make_unique<TLorentzVector>();
    *missingprot_ += *_gamma + *_target - *_pip - *_pim;
    return missingprot_->P();
  } else
    return NAN;
}
float Reaction::prot_theta_lab() {
  if (TwoPion_missingProt()) {
    auto missingprot_ = std::make_unique<TLorentzVector>();
    *missingprot_ += *_gamma + *_target - *_pip - *_pim;
    return missingprot_->Theta() * 180.0 / PI;
  } else
    return NAN;
}
float Reaction::prot_Phi_lab() {
  if (TwoPion_missingProt()) {
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

// float Reaction::pim_momentum_cm() {
//         if (!_is_boosted)
//                 boost();
//         if (TwoPion_missingPim())
//                 return _boosted_pim->P();
//         else
//                 return NAN;
// }

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

// float Reaction::pim_momentum_cm_measured() {
//         if (!_is_boosted)
//                 boost();
//         if (TwoPion_exclusive())
//                 return _boosted_pim_measured->P();
//         else
//                 return NAN;
// }

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
}

void MCReaction::SetMCProton(int i) { _prot_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_P); }

void MCReaction::SetMCPip(int i) { _pip_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIP); }

void MCReaction::SetMCPim(int i) { _pim_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIM); }
// // void MCReaction::SetMCOther(int i) {
// //   _other_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i),
// //   mass[_data->pid(i)]);
// // }
// void MCReaction::CalcMissMass_mc() {
//   auto mm_excl_mc = std::make_unique<TLorentzVector>();

//   *mm_excl_mc += (*_gamma_mc + *_target);
//   *mm_excl_mc -= *_prot_mc;
//   *mm_excl_mc -= *_pip_mc;
//   *mm_excl_mc -= *_pim_mc;
//   _MM2_exclusive_mc = mm_excl_mc->M2();
//   _excl_Energy_mc = mm_excl_mc->E();

//   _rec_x_mu_mom_mc = mm_excl_mc->P();
//   _rec_x_mu_theta_mc = mm_excl_mc->Theta() * 180 / PI;

//   if (mm_excl_mc->Phi() >= 0)
//     _x_mu_phi_mc = (mm_excl_mc->Phi() * 180 / PI);
//   else if (mm_excl_mc->Phi() < 0)
//     _x_mu_phi_mc = ((mm_excl_mc->Phi() + 2 * PI) * 180 / PI);

//   if (_elec_mc->Phi() >= 0)
//     _elec_phi_mc = (_elec_mc->Phi() * 180 / PI);
//   else if (_elec_mc->Phi() < 0)
//     _elec_phi_mc = ((_elec_mc->Phi() + 2 * PI) * 180 / PI);

//   if (_beam->Phi() >= 0)
//     _beam_phi_mc = (_beam->Phi() * 180 / PI);
//   else if (_beam->Phi() < 0)
//     _beam_phi_mc = ((_beam->Phi() + 2 * PI) * 180 / PI);

//   _diff_elec_x_mu_theta_mc = (_elec_mc->Theta() * 180 / PI) - (mm_excl_mc->Theta() * 180 / PI);
//   _diff_elec_x_mu_phi_mc = (_elec_phi_mc - _x_mu_phi_mc);

//   _diff_beam_x_mu_theta_mc = (mm_excl_mc->Theta() * 180 / PI);
//   _diff_beam_x_mu_phi_mc = (_beam_phi_mc - _x_mu_phi_mc);
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
