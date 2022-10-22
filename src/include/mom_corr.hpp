#ifndef MOM_CORR_H_GUARD
#define MOM_CORR_H_GUARD
#include "TROOT.h"
#include "constants.hpp"

namespace mom_corr {

bool is_FD(int prot_status);
bool is_CD(int prot_status);
bool is_lower_band(float mom_, float theta_DCr1_, int status_);

float CD_prot_Emom_corr(float mom_, float theta_);
float FD_prot_Emom_corr_lower(float mom_, float theta_);
float FD_prot_Emom_corr_upper(float mom_, float theta_);

float CD_prot_Eth_corr(float mom_, float theta_);
float FD_prot_Eth_corr_lower(float mom_, float theta_);
float FD_prot_Eth_corr_upper(float mom_, float theta_);

float CD_prot_Eph_corr(float mom_, float theta_, float phi_);
float FD_prot_Eph_corr_lower(float mom_, float theta_, float phi_);
float FD_prot_Eph_corr_upper(float mom_, float theta_, float phi_);

// // Calcuating Energy loss corr parameters
float A_p(float mom_, float theta_, float theta_DCr1_p, int dc_sec);
float B_p(float mom_, float theta_, float theta_DCr1_p, int dc_sec);

// float A_th(float mom_, float theta_, float theta_DCr1_p, int dc_sec);
// float B_th(float mom_, float theta_, float theta_DCr1_p, int dc_sec);
// float C_th(float mom_, float theta_, int dc_sec);

// float A_ph(float mom_, float theta_, float theta_DCr1_p, int dc_sec);
// float B_ph(float mom_, float theta_, float theta_DCr1_p, int dc_sec);
// float C_ph(float mom_, float theta_, float theta_DCr1_p, int dc_sec);

// pip energy loss correction functions
float CD_pip_Emom_corr(float mom_, float theta_);
float FD_pip_Emom_corr_lower(float mom_, float theta_);
float FD_pip_Emom_corr_upper(float mom_, float theta_);

float CD_pip_Eth_corr(float mom_, float theta_);
float FD_pip_Eth_corr_lower(float mom_, float theta_);
float FD_pip_Eth_corr_upper(float mom_, float theta_);

float CD_pip_Eph_corr(float mom_, float theta_, float phi_);
float FD_pip_Eph_corr_lower(float mom_, float theta_, float phi_);
float FD_pip_Eph_corr_upper(float mom_, float theta_, float phi_);

// pim energy loss correction functions
float CD_pim_Emom_corr(float mom_, float theta_);
float FD_pim_Emom_corr_lower(float mom_, float theta_);
float FD_pim_Emom_corr_upper(float mom_, float theta_);

float CD_pim_Eth_corr(float mom_, float theta_);
float FD_pim_Eth_corr_lower(float mom_, float theta_);
float FD_pim_Eth_corr_upper(float mom_, float theta_);

float CD_pim_Eph_corr(float mom_, float theta_, float phi_);
float FD_pim_Eph_corr_lower(float mom_, float theta_, float phi_);
float FD_pim_Eph_corr_upper(float mom_, float theta_, float phi_);

// hadron mom corrections
double dppC(float Px, float Py, float Pz, int sec, int ivec);

float CD_prot_Hmom_corr(float mom_, float phi_);
float FD_prot_Hmom_corr_lower(float mom_, float dc_sec);
float FD_prot_Hmom_corr_upper(float mom_, float dc_sec);

float CD_pip_Hmom_corr(float mom_, float phi_);
float FD_pip_Hmom_corr_lower(float mom_, float dc_sec);
float FD_pip_Hmom_corr_upper(float mom_, float dc_sec);

float CD_pim_Hmom_corr(float mom_, float phi_);
float FD_pim_Hmom_corr_lower(float mom_, float dc_sec);
float FD_pim_Hmom_corr_upper(float mom_, float dc_sec);

}  // namespace mom_corr

#endif
