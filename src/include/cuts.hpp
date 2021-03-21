// /**************************************/
// /*																		*/
// /*  Created by Nick Tyler             */
// /*	University Of South Carolina      */
// /**************************************/
//
// #ifndef CUTS_H_GUARD
// #define CUTS_H_GUARD
// #include "branches.hpp"
// #include "constants.hpp"
// #include "deltat.hpp"
// #include "physics.hpp"
//
// class Cuts {
// private:
// std::shared_ptr<Branches12> _data = nullptr;
// std::shared_ptr<Delta_T> _dt = nullptr;
//
// public:
// Cuts(const std::shared_ptr<Branches12>& data);
// Cuts(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt);
// ~Cuts();
//
// bool ElectronCuts();
// bool FiducialCuts();
// bool CC_nphe_cut();
// bool EC_outer_vs_EC_inner_cut();
// bool EC_sampling_fraction_cut();
// bool EC_hit_position_fiducial_cut_homogeneous();
// bool DC_fiducial_cut_XY();
// bool DC_z_vertex_cut();
//
//
//
// bool IsPip(int i);
// bool IsProton(int i);
// bool IsPim(int i);
// bool Hadron_Chi2pid_cut(int i);
// bool Hadron_Delta_vz_cut(int i);
// bool DC_fiducial_cut_theta_phi(int i);
// // bool IsmissingPim(int i);
//
// // bool IsPositive(int i);
// // bool IsFDPositive(int i);
// // bool IsCDPositive(int i);
// };
//
// #endif

#ifndef CUTS_H_GUARD
#define CUTS_H_GUARD
#include "branches.hpp"
#include "constants.hpp"
#include "deltat.hpp"
#include "physics.hpp"

class Cuts {
protected:
std::shared_ptr<Branches12> _data = nullptr;
std::shared_ptr<Delta_T> _dt = nullptr;

public:
Cuts(const std::shared_ptr<Branches12>& data);
Cuts(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt);
~Cuts();

bool ElectronCuts();
bool FiducialCuts();
bool IsPip(int i);
bool IsProton(int i);
bool IsPim(int i);
};

class rga_Cuts : public Cuts {
public:
rga_Cuts(const std::shared_ptr<Branches12>& data) : Cuts(data) {
}
rga_Cuts(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt) : Cuts(data, dt){
};
};

class rgf_Cuts : public Cuts {
public:
rgf_Cuts(const std::shared_ptr<Branches12>& data) : Cuts(data) {
}
rgf_Cuts(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt) : Cuts(data, dt){
};
};

class uconn_Cuts : public Cuts {
public:
uconn_Cuts(const std::shared_ptr<Branches12>& data) : Cuts(data) {
}
uconn_Cuts(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt) : Cuts(data, dt){
};
bool ElectronCuts();

// bool CC_nphe_cut(double nphe);
bool CC_nphe_cut();
bool EC_outer_vs_EC_inner_cut();
bool EC_sampling_fraction_cut();
bool EC_hit_position_fiducial_cut_homogeneous();
bool DC_fiducial_cut_XY();
bool DC_z_vertex_cut();

bool HadronsCuts(int i);
bool DC_fiducial_cut_theta_phi(int i);
bool Hadron_Delta_vz_cut(int i);
bool Hadron_Chi2pid_cut(int i);
};
#endif
