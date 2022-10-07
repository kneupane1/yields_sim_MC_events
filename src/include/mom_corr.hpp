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

float CD_prot_Eph_corr(float mom_, float theta_, float phi_P);
float FD_prot_Eph_corr_lower(float mom_, float theta_, float phi_P);
float FD_prot_Eph_corr_upper(float mom_, float theta_, float phi_P);

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

float CD_pip_Eph_corr(float mom_, float theta_, float phi_P);
float FD_pip_Eph_corr_lower(float mom_, float theta_, float phi_P);
float FD_pip_Eph_corr_upper(float mom_, float theta_, float phi_P);

}  // namespace mom_corr

#endif
