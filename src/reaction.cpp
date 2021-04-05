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

Reaction::~Reaction() {
}

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
        // auto mm_mpip = std::make_unique<TLorentzVector>();
        // auto mm_mprot = std::make_unique<TLorentzVector>();
        auto mm_excl = std::make_unique<TLorentzVector>();

        *mm += (*_gamma + *_target);

        //if (TwoPion_missingPim()) {
        *mm -= *_prot;
        *mm -= *_pip;
        //*mm -= *_pim;
        _MM = mm->M();
        _MM2 = mm->M2();
        //         // _x_mu_E = mm->E();
        //         // _x_mu_P = mm->P();
        //         // _x_mu_Px = mm->Px();
        //         // _x_mu_Py = mm->Py();
        //         // _x_mu_Pz = mm->Pz();
        //         // _x_mu_theta = mm->Theta() * RAD2DEG;
        //         // _x_mu_m2 = mm->E() * mm->E() - mm->P() * mm->P();
        //         // _x_mu_m = mm->E() - mm->P();
        //
//        }
        /*if (TwoPion_exclusive()) {
         * mm -= *_prot;
         * mm -= *_pip;
                _MM = mm->M();
                _MM2 = mm->M2();
                // _x_mu_E = mm->E();
                // _x_mu_P = mm->P();
                // _x_mu_Px = mm->Px();
                // _x_mu_Py = mm->Py();
                // _x_mu_Pz = mm->Pz();
                // _x_mu_theta = mm->Theta() * RAD2DEG;
                // _x_mu_m2 = mm->E() * mm->E() - mm->P() * mm->P();
                // _x_mu_m = mm->E() - mm->P();
                //
         * mm_excl += (*_gamma + *_target);
         * mm_excl -= *_prot;
         * mm_excl -= *_pip;
         * mm_excl -= *_pim;
                _MM2_exclusive =mm_excl->M2();
                // _excl_Energy = mm_excl->E();
                //
                // *mm_mpip += (*_gamma + *_target);
                // *mm_mpip -= *_prot;
                // *mm_mpip -= *_pim;
                // _MM2_mpip = mm_mpip->M2();
                //
                // *mm_mprot += (*_gamma + *_target);
                // *mm_mprot -= *_pip;
                // *mm_mprot -= *_pim;
                // _MM2_mprot = mm_mprot->M2();

           }*/
        // if (TwoPion_missingPip()) {
        //         *mm_mpip += (*_gamma + *_target);
        //         *mm_mpip -= *_prot;
        //         //*mm_mpip -= *_pip;
        //         *mm_mpip -= *_pim;
        //         // _MM2_mpip = mm_mpip->M2();
        // }
        // if (TwoPion_missingProt()) {
        //         *mm_mprot += (*_gamma + *_target);
        //         *mm_mprot -= *_pip;
        //         *mm_mprot -= *_pim;
        //         // _MM2_mprot = mm_mprot->M2();
        //
        //         // } else if (SingleP()) {
        //         //   *mm -= *_prot;
        //         //   _MM = mm->M();
        //         //   _MM2 = mm->M2();
        // }
}

float Reaction::MM() {
        if (_MM != _MM)
                CalcMissMass();
        return _MM;
}
float Reaction::MM2() {
        if (_MM2 != _MM2)
                CalcMissMass();
        return _MM2;
}
// float Reaction::MM2_exclusive() {
//         if (_MM2_exclusive != _MM2_exclusive)
//                 CalcMissMass();
//         return _MM2_exclusive;
// }
// float Reaction::MM2_mpip() {
//         if (_MM2_mpip != _MM2_mpip)
//                 CalcMissMass();
//         return _MM2_mpip;
// }
// float Reaction::MM2_mprot() {
//         if (_MM2_mprot != _MM2_mprot)
//                 CalcMissMass();
//         return _MM2_mprot;
// }
// float Reaction::Energy_excl() {
//         if (_excl_Energy != _excl_Energy)
//                 CalcMissMass();
//         //  std::cout << "_x_mu_p  " << _x_mu->E() << '\n';
//         //  if (_x_mu_E > 0)
//         return _excl_Energy;
//         // else
//         // return NAN;
// }
// float Reaction::pim_momentum() {
//         if (TwoPion_missingPim()) {
//                 auto missingpim_ = std::make_unique<TLorentzVector>();
//                 *missingpim_ += *_gamma + *_target - *_prot - *_pip;
//                 return missingpim_->P();
//         } else
//                 return NAN;
// }
// float Reaction::pim_theta_lab() {
//         if (TwoPion_missingPim()) {
//                 auto missingpim_ = std::make_unique<TLorentzVector>();
//                 *missingpim_ += *_gamma + *_target - *_prot - *_pip;
//                 return missingpim_->Theta() * 180.0 / PI;
//         } else
//                 return NAN;
// }
// float Reaction::pim_Phi_lab() { /////////////////////////////////////work here
//         auto missingpim_ = std::make_unique<TLorentzVector>();
//         *missingpim_ += *_gamma + *_target - *_prot - *_pip;
//         if (TwoPion_missingPim()) {
//                 if (missingpim_->Phi() > 0)
//                         return missingpim_->Phi() * 180 / PI;
//                 else if (missingpim_->Phi() < 0)
//                         return (missingpim_->Phi() + 2 * PI) * 180 / PI;
//                 else
//                         return NAN;
//         } else
//                 return NAN;
// }
// float Reaction::pim_momentum_measured() {
//         if (TwoPion_exclusive())
//                 return _pim->P();
//         else
//                 return NAN;
// }
//
// float Reaction::pim_theta_lab_measured() { /////////////////////////////////////work here
//
//         if (TwoPion_exclusive())
//                 //return _boosted_pim->Theta() * 180.0 / PI;
//                 return _pim->Theta() * 180.0 / PI;
//         else
//                 return NAN;
// }
//
// float Reaction::pim_Phi_lab_measured() {/////////////////////////////////////work here
//         if (TwoPion_exclusive()) {
//                 if (_pim->Phi() > 0)
//                         return _pim->Phi() * 180 / PI;
//                 else if (_pim->Phi() < 0)
//                         return (_pim->Phi() + 2 * PI) * 180 / PI;
//                 else
//                         return NAN;
//         } else
//                 return NAN;
// }


std::string Reaction::CsvHeader() {
        return "e_rec_p,e_rec_theta,e_rec_phi,e_sec\n";
}
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
        _boosted_pim = std::make_unique<TLorentzVector>(*_gamma + *_target - *_prot -
                                                        *_pip); //(*_pim);
        _boosted_gamma = std::make_unique<TLorentzVector>(*_gamma);
        _boosted_pim_measured = std::make_unique<TLorentzVector>(*_pim);

        TRotation rot;
        _boosted_gamma->Transform(rot);
        float_t beta_1 =
                ((sqrt(_boosted_gamma->E() * _boosted_gamma->E() + _Q2)) / (_boosted_gamma->E() + MASS_P));
        TVector3 uz = _boosted_gamma->Vect().Unit(); // uit vector along virtual photon
        TVector3 ux = ((_beam->Vect()).Cross(_elec->Vect())).Unit(); // unit vector along e cross e'
        ux.Rotate(3. * PI / 2, uz); // rotating ux by 3pi/2 with uz as axis of roration
        rot.SetZAxis(uz, ux).Invert(); // setting TRotation rot
        _boosted_prot->Transform(rot);
        _boosted_prot->Boost(0, 0, -beta_1);
        _boosted_pip->Transform(rot);
        _boosted_pip->Boost(0, 0, -beta_1);
        _boosted_pim->Transform(rot);
        _boosted_pim->Boost(0, 0, -beta_1);
        _boosted_gamma->Boost(0, 0, -beta_1);

        _boosted_pim_measured->Transform(rot);
        _boosted_pim_measured->Boost(0, 0, -beta_1);
        // -beta ko value (0.5 to -0.5 huda
        // samma value aauchha nattra aaudyna)

        _prot_Vect3 = _boosted_prot->Vect();
        _pip_Vect3 = _boosted_pip->Vect();
        _pim_Vect3 = _boosted_pim_measured->Vect();
}

float_t Reaction::scalar_triple_product(){
        if (!_is_boosted) boost();
        if(TwoPion_exclusive()) {
                return (_prot_Vect3.Dot(_pip_Vect3.Cross(_pim_Vect3)));

        } else
                return NAN;
}


float Reaction::pim_momentm_cm() {
        if (!_is_boosted)
                boost();
        if (TwoPion_missingPim())
                return _boosted_pim->P();
        else
                return NAN;
}

float Reaction::pim_theta_cm() {
        if (!_is_boosted)
                boost();
        if (TwoPion_missingPim())
                return _boosted_pim->Theta() * 180.0 / PI;
        else
                return NAN;
}

float Reaction::pim_Phi_cm() {
        if (!_is_boosted)
                boost();
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

float Reaction::pim_momentm_cm_measured() {
        if (!_is_boosted)
                boost();
        if (TwoPion_exclusive())
                return _boosted_pim_measured->P();
        else
                return NAN;
}

float Reaction::pim_theta_cm_measured() {
        if (!_is_boosted)
                boost();
        if (TwoPion_exclusive())
                return _boosted_pim_measured->Theta() * 180.0 / PI;
        else
                return NAN;
}

float Reaction::pim_Phi_cm_measured() {
        if (!_is_boosted)
                boost();
        if (TwoPion_exclusive()) {
                if (_boosted_pim_measured->Phi() > 0)
                        return _boosted_pim_measured->Phi() * 180 / PI;
                else if (_boosted_pim_measured->Phi() < 0)
                        return (_boosted_pim_measured->Phi() + 2 * PI) * 180 / PI;
                else
                        return NAN;
        } else
                return NAN;
}

MCReaction::MCReaction(const std::shared_ptr<Branches12> &data,
                       float beam_enrgy) {
        _data = data;
        if (!_data->mc())
                _data->mc_branches();
        _beam = std::make_unique<TLorentzVector>();
        _beam_energy = beam_enrgy;
        _weight_mc = _data->mc_weight();
        _beam->SetPxPyPzE(0.0, 0.0,
                          sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E),
                          _beam_energy);

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

void MCReaction::SetMCProton(int i) {
        _prot_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_P);
}

void MCReaction::SetMCPip(int i) {
        _pip_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIP);
}

void MCReaction::SetMCPim(int i) {
        _pim_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIM);
}
// void MCReaction::SetMCOther(int i) {
//   _other_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i),
//   mass[_data->pid(i)]);
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
