#ifndef MOM_CORR_H_GUARD
#define MOM_CORR_H_GUARD
#include "TROOT.h"
#include "constants.hpp"

namespace mom_corr {

bool is_FD(int dc_sec);
bool is_CD(int dc_sec);
bool is_lower_band(float mom_P, float theta_DCr1_p, int dc_sec);

// Calcuating Energy loss corr parameters
float A_p(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec);
float B_p(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec);

float A_th(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec);
float B_th(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec);
float C_th(float mom_P, float theta_P, int dc_sec);

float A_ph(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec);
float B_ph(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec);
float C_ph(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec);

}  // namespace mom_corr

#endif
