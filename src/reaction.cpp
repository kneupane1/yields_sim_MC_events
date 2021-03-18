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
        *mm += (*_gamma + *_target);
        if (SinglePip() || NeutronPip()) {
                *mm -= *_pip;
                _MM = mm->M();
                _MM2 = mm->M2();
        } else if (TwoPion()) {
                *mm -= *_prot;
                *mm -= *_pip;
                *mm -= *_pim;
                _MM = mm->M();
                _MM2 = mm->M2();
        } else if (ProtonPim()) {
                *mm -= *_prot;
                *mm -= *_pim;
                _MM = mm->M();
                _MM2 = mm->M2();
        } else if (SingleP()) {
                *mm -= *_prot;
                _MM = mm->M();
                _MM2 = mm->M2();
        }
}

void Reaction::boost() {
        // Angles gotten from picture in KPark thesis page 11
        // May be wrong still, need to check

        auto _com_ = *_target + (*_beam - *_elec);

        auto com = std::make_shared<TLorentzVector>(_com_.X(), _com_.Y(), _com_.Z(), _com_.E());
        auto elec_boosted = std::make_shared<TLorentzVector>(_elec->X(), _elec->Y(), _elec->Z(), _elec->E());
        auto gamma_boosted = std::make_shared<TLorentzVector>(_gamma->X(), _gamma->Y(), _gamma->Z(), _gamma->E());

        auto beam_boosted = std::make_shared<TLorentzVector>(_beam->X(), _beam->Y(), _beam->Z(), _beam->E());
        auto pip_boosted = std::make_shared<TLorentzVector>(_pip->X(), _pip->Y(), _pip->Z(), _pip->E());

        //! Calculate rotation: taken from Evan's phys-ana-omega on 08-05-13
        // Copied and modified from Arjun's code

        //  auto uz = _gamma->Vect().Unit();
        // auto ux = _beam->Vect().Cross(_elec->Vect()).Unit();
        // ROOT::Math::VectorUtil::Rotate(ux, uz);

        TVector3 uz = gamma_boosted->Vect().Unit();
        TVector3 ux = (beam_boosted->Vect().Cross(elec_boosted->Vect())).Unit();
        ux.Rotate(-PI / 2, uz);
        TRotation r3; // = new TRotation();
        r3.SetZAxis(uz, ux).Invert();
        //! _w and _q are in z-direction
        TVector3 boost(-1 * com->BoostVector());
        TLorentzRotation r4(r3); //*_boost);
        r4 *= boost;        //*_3rot;

        gamma_boosted->Transform(r4);
        elec_boosted->Transform(r4);
        beam_boosted->Transform(r4);
        pip_boosted->Transform(r4);

        auto _temp = physics::fourVec(pip_boosted->X(), pip_boosted->Y(), pip_boosted->Z(), pip_boosted->M());
        _theta_e = elec_boosted->Theta();
        _theta_star = _temp->Theta();
        _phi_star = physics::phi_boosted(_temp);
}

float Reaction::MM() {
        if (_MM != _MM) CalcMissMass();
        return _MM;
}
float Reaction::MM2() {
        if (_MM2 != _MM2) CalcMissMass();
        return _MM2;
}

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


MCReaction::MCReaction(const std::shared_ptr<Branches12> &data,
                       float beam_enrgy) {
        _data = data;
        if (!_data->mc())
                _data->mc_branches();
        _beam = std::make_unique<TLorentzVector>();
        _beam_energy = beam_enrgy;
        // _weight_mc = _data->mc_weight();
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
