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

  // // // // Can calculate W and Q2 here (useful for simulations as sim do not have elec mom corrections)
  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);
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

  _thetaDC_r1_Prot = RAD2DEG * (atan2(sqrt(pow(_data->dc_r1_x(i), 2) + pow(_data->dc_r1_y(i), 2)), _data->dc_r1_z(i)));

  // _prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
  // _mom_corr_prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);

  _sectorProt = _data->dc_sec(i);

  _prot_mom_uncorr = _Energy_loss_uncorr_prot->P();

  _prot_theta_uncorr = _Energy_loss_uncorr_prot->Theta() * 180 / PI;

  if (_Energy_loss_uncorr_prot->Phi() > 0)
    _prot_phi_uncorr = _Energy_loss_uncorr_prot->Phi() * 180 / PI;
  else if (_Energy_loss_uncorr_prot->Phi() < 0)
    _prot_phi_uncorr = (_Energy_loss_uncorr_prot->Phi() + 2 * PI) * 180 / PI;

  _is_FD_Prot = objMomCorr->is_FD(_prot_status);
  _is_CD_Prot = objMomCorr->is_CD(_prot_status);
  // _is_lower_band = objMomCorr->is_lower_band(_prot_mom_uncorr, _thetaDC_r1_Prot, _prot_status);

  if (_is_CD_Prot) {
    _prot_mom_tmt = _prot_mom_uncorr;
  }
  if (_is_FD_Prot) {
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

  _prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);  // energy loss + FD had corr
  _mom_corr_prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);
}
void Reaction::SetPip(int i) {
  _numPip++;
  _numPos++;
  _hasPip = true;
  _pip_status = abs(_data->status(i));
  _sectorPip = _data->dc_sec(i);
  _thetaDC_r1_Pip = RAD2DEG * (atan2(sqrt(pow(_data->dc_r1_x(i), 2) + pow(_data->dc_r1_y(i), 2)), _data->dc_r1_z(i)));

  _Energy_loss_uncorr_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
  _pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
  _mom_corr_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
}

void Reaction::SetPim(int i) {
  _numPim++;
  _numNeg++;
  _hasPim = true;
  _pim_status = abs(_data->status(i));
  _sectorPim = _data->dc_sec(i);
  _thetaDC_r1_Pim = RAD2DEG * (atan2(sqrt(pow(_data->dc_r1_x(i), 2) + pow(_data->dc_r1_y(i), 2)), _data->dc_r1_z(i)));

  _pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
  _mom_corr_pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
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

void Reaction::CalcMissMass() {
  auto mm = std::make_unique<TLorentzVector>();
  *mm += (*_gamma + *_target);

  if (TwoPion_missingPim()) {
    *mm -= *_prot;
    *mm -= *_pip;
    // *mm -= *_pim;
    _MM = mm->M();
    _MM2 = mm->M2();

    //   // _rec_pim_mom = mm->P();
    //   // _rec_pim_theta = mm->Theta() * 180 / PI;

    //   // if (mm->Phi() >= 0)
    //   //   _rec_pim_phi = (mm->Phi() * 180 / PI);
    //   // else if (mm->Phi() < 0)
    //   //   _rec_pim_phi = ((mm->Phi() + 2 * PI) * 180 / PI);
  }
}

float Reaction::MM() {
  if (_MM != _MM) CalcMissMass();
  return _MM;
}
float Reaction::MM2() {
  if (_MM2 != _MM2) CalcMissMass();
  return _MM2;
}
float Reaction::pim_momentum() {
  // if (_rec_pim_mom != _rec_pim_mom) CalcMissMass();

  if (TwoPion_missingPim()) {
    // if (TwoPion_exclusive()) {
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

  // // // momentum uncorrected version
  _boosted_prot = std::make_unique<TLorentzVector>(*_prot);
  _boosted_pip = std::make_unique<TLorentzVector>(*_pip);
  _boosted_pim = std::make_unique<TLorentzVector>(*_gamma + *_target - *_prot - *_pip);  //(*_pim);
  _boosted_gamma = std::make_unique<TLorentzVector>(*_gamma);
  _boosted_pim_measured = std::make_unique<TLorentzVector>(*_pim);

  // // momentum corrected version
  // _boosted_prot = std::make_unique<TLorentzVector>(*_mom_corr_prot);
  // _boosted_pip = std::make_unique<TLorentzVector>(*_mom_corr_pip);
  // _boosted_pim = std::make_unique<TLorentzVector>(*_gamma + *_target - *_mom_corr_prot - *_mom_corr_pip); //(*_pim);
  // _boosted_gamma = std::make_unique<TLorentzVector>(*_gamma);
  // _boosted_pim_measured = std::make_unique<TLorentzVector>(*_mom_corr_pim);

  TRotation rot;

  TVector3 uz = _boosted_gamma->Vect().Unit();                  // uit vector along virtual photon
  TVector3 ux = ((_beam->Vect()).Cross(_elec->Vect())).Unit();  // unit vector along e cross e'

  ux.Rotate(3. * PI / 2, uz);     // rotating ux by 3pi/2 with uz as axis of roration
  rot.SetZAxis(uz, ux).Invert();  // setting TRotation rot

  _boosted_gamma->Transform(rot);
  _boosted_prot->Transform(rot);
  _boosted_pip->Transform(rot);
  _boosted_pim->Transform(rot);

  // note beta is calculated only after transforming gamma

  float_t beta_1 = ((sqrt(_boosted_gamma->E() * _boosted_gamma->E() + _Q2)) / (_boosted_gamma->E() + MASS_P));

  _boosted_prot->Boost(0, 0, -beta_1);
  _boosted_pip->Boost(0, 0, -beta_1);
  _boosted_pim->Boost(0, 0, -beta_1);
  _boosted_gamma->Boost(0, 0, -beta_1);

  _boosted_pim_measured->Transform(rot);
  _boosted_pim_measured->Boost(0, 0, -beta_1);
  // -beta ko value (0.5 to -0.5 huda
  // samma value aauchha nattra aaudyna)
}
// float Reaction::weight() { return _weight; }

//////////////////

void Reaction::invMassPpim() {
  if (!_is_boosted) boost();
  auto inv_Ppim = std::make_unique<TLorentzVector>();
  *inv_Ppim += *_boosted_prot;
  *inv_Ppim += *_boosted_pim;

  // *inv_Ppim += (*_boosted_gamma + *_target - *_boosted_prot - *_boosted_pip);
  if (TwoPion_missingPim()) _inv_Ppim = inv_Ppim->M();
}
void Reaction::invMasspippim() {
  if (!_is_boosted) boost();
  auto inv_pip_pim = std::make_unique<TLorentzVector>();
  *inv_pip_pim += *_boosted_pip;
  *inv_pip_pim += *_boosted_pim;

  // *inv_pip_pim += (*_boosted_gamma + *_target - *_boosted_prot - *_boosted_pip);
  if (TwoPion_missingPim()) _inv_pip_pim = inv_pip_pim->M();
}

void Reaction::invMassPpip() {
  if (!_is_boosted) boost();
  auto inv_Ppip = std::make_unique<TLorentzVector>();
  *inv_Ppip += *_boosted_prot;
  *inv_Ppip += *_boosted_pip;

  if (TwoPion_missingPim()) _inv_Ppip = inv_Ppip->M();
}

void Reaction::W_2pi_P() {
  auto W_P2pi = std::make_unique<TLorentzVector>();
  *W_P2pi += *_prot;
  *W_P2pi += *_pip;
  *W_P2pi += (*_gamma + *_target - *_prot - *_pip);

  if (TwoPion_missingPim()) _W_P2pi = W_P2pi->M();
}

float Reaction::inv_Ppip() {
  if (_inv_Ppip != _inv_Ppip) invMassPpip();
  return _inv_Ppip;
}
float Reaction::inv_Ppim() {
  if (_inv_Ppim != _inv_Ppim) invMassPpim();
  return _inv_Ppim;
}
float Reaction::inv_pip_pim() {
  if (_inv_pip_pim != _inv_pip_pim) invMasspippim();
  return _inv_pip_pim;
}
float Reaction::w_P2pi_rec() {
  if (_W_P2pi != _W_P2pi) W_2pi_P();
  return _W_P2pi;
}

//////////////
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

MCReaction::MCReaction(const std::shared_ptr<Branches12>& data, float beam_energy) {
  _data = data;
  if (!_data->mc()) _data->mc_branches();
  _beam = std::make_unique<TLorentzVector>();
  // _beam_energy = 10.6041;
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

/////////////////////////////////////////////////////////////////////////////
void MCReaction::boost_mc() {
  _is_boosted_mc = true;
  _boosted_prot_mc = std::make_unique<TLorentzVector>(*_prot_mc);
  _boosted_pip_mc = std::make_unique<TLorentzVector>(*_pip_mc);
  _boosted_pim_mc = std::make_unique<TLorentzVector>(*_gamma_mc + *_target - *_prot_mc - *_pip_mc);  //(*_pim);
  // _boosted_pim_mc = std::make_unique<TLorentzVector>(*_pim_mc);

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

  // beta is calculated only after transfer in gamma
  float_t beta_1_mc =
      ((sqrt(_boosted_gamma_mc->E() * _boosted_gamma_mc->E() + _Q2_mc)) / (_boosted_gamma_mc->E() + MASS_P));

  _boosted_prot_mc->Boost(0, 0, -beta_1_mc);
  _boosted_pip_mc->Boost(0, 0, -beta_1_mc);
  _boosted_pim_mc->Boost(0, 0, -beta_1_mc);
  _boosted_gamma_mc->Boost(0, 0, -beta_1_mc);  // -beta ko value (0.5 to -0.5 huda
                                               // samma value aauchha nattra aaudyna)
}

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
float MCReaction::MCgamma_Phi_thrown() {
  if (!_is_boosted_mc) boost_mc();
  if (_boosted_gamma_mc->Phi() > 0)
    return _boosted_gamma_mc->Phi() * 180 / PI;
  else if (_boosted_gamma_mc->Phi() < 0)
    return (_boosted_gamma_mc->Phi() + 2 * PI) * 180 / PI;
  // else
  return NAN;
}
float MCReaction::MCprot_Phi_thrown() {
  if (!_is_boosted_mc) boost_mc();

  if (_boosted_prot_mc->Phi() > 0)
    return _boosted_prot_mc->Phi() * 180 / PI;
  else if (_boosted_prot_mc->Phi() < 0)
    return (_boosted_prot_mc->Phi() + 2 * PI) * 180 / PI;
  // else
  return NAN;
}
float MCReaction::MCpip_Phi_thrown() {
  if (!_is_boosted_mc) boost_mc();
  if (_boosted_pip_mc->Phi() > 0)
    return _boosted_pip_mc->Phi() * 180 / PI;
  else if (_boosted_pip_mc->Phi() < 0)
    return (_boosted_pip_mc->Phi() + 2 * PI) * 180 / PI;
  // else
  return NAN;
}
float MCReaction::MCpim_Phi_thrown() {
  if (!_is_boosted_mc) boost_mc();

  if (_boosted_pim_mc->Phi() > 0)
    return _boosted_pim_mc->Phi() * 180 / PI;
  else if (_boosted_pim_mc->Phi() < 0)
    return (_boosted_pim_mc->Phi() + 2 * PI) * 180 / PI;
  // else
  return NAN;
}
float MCReaction::MCprot_theta_thrown() {
  if (!_is_boosted_mc) boost_mc();
  return _boosted_prot_mc->Theta() * 180.0 / PI;

  return NAN;
}
float MCReaction::MCpip_theta_thrown() {
  if (!_is_boosted_mc) boost_mc();
  return _boosted_pip_mc->Theta() * 180.0 / PI;
  return NAN;
}
float MCReaction::MCpim_theta_thrown() {
  if (!_is_boosted_mc) boost_mc();
  return _boosted_pim_mc->Theta() * 180.0 / PI;

  return NAN;
}

void MCReaction::MCAlphaCalc() {
  //  Float_t m_proton, m_pip, beta;
  Float_t a_gamma, b_gamma, a_beta, b_beta;
  TVector3 Vect3_gamma, Vect3_beta, V3_anti_z(0, 0, -1);
  float alpha_PPIp_piPIm_mc;
  float alpha_pippim_pipf_mc;
  float alpha_PPIm_piPIp_mc;

  if (!_is_boosted_mc) boost_mc();
  // 1 this one is used for
  a_gamma = sqrt(1. / (1 - pow((_boosted_pim_mc->Vect().Unit() * V3_anti_z),
                               2)));  // V3_anti_z(0,0,-1);
  b_gamma = -(_boosted_pim_mc->Vect().Unit() * V3_anti_z) * a_gamma;
  Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _boosted_pim_mc->Vect().Unit();

  a_beta = sqrt(1. / (1 - pow((_boosted_pim_mc->Vect().Unit() * _boosted_pip_mc->Vect().Unit()), 2)));
  b_beta = -(_boosted_pim_mc->Vect().Unit() * _boosted_pip_mc->Vect().Unit()) * a_beta;
  Vect3_beta = a_beta * _boosted_pip_mc->Vect().Unit() + b_beta * _boosted_pim_mc->Vect().Unit();

  alpha_PPIp_piPIm_mc = (180. / PI) * acos(Vect3_gamma * Vect3_beta);
  if (Vect3_gamma.Cross(Vect3_beta) * _boosted_pim_mc->Vect() < 0) alpha_PPIp_piPIm_mc = 360. - alpha_PPIp_piPIm_mc;

  //α[pπ+][p'π−]
  /// 2
  a_gamma = sqrt(1. / (1 - pow((_boosted_prot_mc->Vect().Unit() * V3_anti_z), 2)));
  b_gamma = -(_boosted_prot_mc->Vect().Unit() * V3_anti_z) * a_gamma;
  Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _boosted_prot_mc->Vect().Unit();

  a_beta = sqrt(1. / (1 - pow((_boosted_prot_mc->Vect().Unit() * _boosted_pip_mc->Vect().Unit()), 2)));
  b_beta = -(_boosted_prot_mc->Vect().Unit() * _boosted_pip_mc->Vect().Unit()) * a_beta;
  Vect3_beta = a_beta * _boosted_pip_mc->Vect().Unit() + b_beta * _boosted_prot_mc->Vect().Unit();

  alpha_pippim_pipf_mc = (180. / PI) * acos(Vect3_gamma * Vect3_beta);

  if (Vect3_gamma.Cross(Vect3_beta) * _boosted_prot_mc->Vect() < 0) alpha_pippim_pipf_mc = 360. - alpha_pippim_pipf_mc;
  //α[pp'][π+π−]

  /// 3
  a_gamma = sqrt(1. / (1 - pow((_boosted_pip_mc->Vect().Unit() * V3_anti_z), 2)));
  b_gamma = -(_boosted_pip_mc->Vect().Unit() * V3_anti_z) * a_gamma;
  Vect3_gamma = a_gamma * V3_anti_z + b_gamma * _boosted_pip_mc->Vect().Unit();

  a_beta = sqrt(1. / (1 - pow((_boosted_pip_mc->Vect().Unit() * _boosted_pim_mc->Vect().Unit()), 2)));
  b_beta = -(_boosted_pip_mc->Vect().Unit() * _boosted_pim_mc->Vect().Unit()) * a_beta;
  Vect3_beta = a_beta * _boosted_pim_mc->Vect().Unit() + b_beta * _boosted_pip_mc->Vect().Unit();

  alpha_PPIm_piPIp_mc = (180. / PI) * acos(Vect3_gamma * Vect3_beta);

  if (Vect3_gamma.Cross(Vect3_beta) * _boosted_pip_mc->Vect() < 0) alpha_PPIm_piPIp_mc = 360. - alpha_PPIm_piPIp_mc;
  //α[pπ−][p'π+]

  _alpha_ppip_pipim_mc = alpha_PPIp_piPIm_mc;
  _alpha_pippim_pipf_mc = alpha_pippim_pipf_mc;
  _alpha_ppim_pipip_mc = alpha_PPIm_piPIp_mc;
}

float MCReaction::MCalpha_ppip_pipim_thrown() {  // pipim bhaneko proton
  // initial pim ho?
  if (_alpha_ppip_pipim_mc != _alpha_ppip_pipim_mc) MCAlphaCalc();

  return _alpha_ppip_pipim_mc;
  return NAN;
}
float MCReaction::MCalpha_pippim_pipf_thrown() {  // alpha P (proton initial
  // proton final)
  if (_alpha_pippim_pipf_mc != _alpha_pippim_pipf_mc) MCAlphaCalc();

  return _alpha_pippim_pipf_mc;
  return NAN;
}
float MCReaction::MCalpha_ppim_pipip_thrown() {  // alpha pip (proton initial
  // pip)
  if (_alpha_ppim_pipip_mc != _alpha_ppim_pipip_mc) MCAlphaCalc();

  return _alpha_ppim_pipip_mc;
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
