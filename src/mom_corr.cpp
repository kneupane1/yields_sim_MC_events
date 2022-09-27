
#include "mom_corr.hpp"

namespace mom_corr {

bool is_FD(int dc_sec) {
  if (dc_sec >= 1 && dc_sec <= 6)
    return true;
  else
    return false;
}

bool is_CD(int dc_sec) {
  if (dc_sec < 1 || dc_sec > 6)
    return true;
  else
    return false;
}
bool is_lower_band(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec) {
  if (dc_sec >= 1 && dc_sec <= 6) {
    if (theta_DCr1_p < 53.14680163254601 + 79.61307254040804 * pow((mom_P - 0.3), 0.05739232362022314)) {
      return true;
    } else
      return false;
  } else
    return false;
}
// energy loss corrections parameters for momentum of proton
float A_p(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec) {
  if (dc_sec >= 1 && dc_sec <= 6) {
    if (theta_DCr1_p < 53.14680163254601 + 79.61307254040804 * pow((mom_P - 0.3), 0.05739232362022314)) {
      return -0.00051894 - 0.00018104 * theta_P;
      //   Ap = − 0.00051894 − 0.00018104 × θ

    } else
      return -3.03346359e-1 + 1.83368163e-2 * theta_P - 2.86486404e-4 * theta_P * theta_P;
    //   Ap = − 3.03346359 × 10−1 + 1.83368163 × 10−2 × θ − 2.86486404 × 10−4 × θ2
  } else
    return 1.93686914 - 0.116288824 * theta_P + 0.00223685833 * theta_P * theta_P -
           1.40771969e-5 * theta_P * theta_P * theta_P;
  // Ap =1.93686914 − 0.116288824 × θ + 0.00223685833 × θ2 − 1.40771969 × 10−5 × θ3
}

float B_p(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec) {
  if (dc_sec >= 1 && dc_sec <= 6) {
    if (theta_DCr1_p < 53.14680163254601 + 79.61307254040804 * pow((mom_P - 0.3), 0.05739232362022314)) {
      return 3.29466917e-3 + 5.73663160e-4 * theta_P - 1.40807209e-5 * theta_P * theta_P;
      //   Bp =3.29466917 × 10−3 + 5.73663160 × 10−4 × θ − 1.40807209 × 10−5 × θ2.
    } else
      return 2.01023276e-1 - 1.13312215e-2 * theta_P + 1.82487916e-4 * theta_P * theta_P;
    // Bp = 2.01023276 × 10−1 − 1.13312215 × 10−2 × θ + 1.82487916 × 10−4 × θ2.
  } else
    return -0.738047800 + 0.0443343685 * theta_P - 8.50985972e-4 * theta_P * theta_P +
           5.36810280e-6 * theta_P * theta_P * theta_P;
  //   Bp = − 0.738047800 + 0.0443343685 × θ − 8.50985972 × 10−4 × θ2 + 5.36810280 × 10−6 × θ3
}

// energy loss corrections parameters for theta angle of proton

float A_th(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec) {
  if (dc_sec >= 1 && dc_sec <= 6) {
    if (theta_DCr1_p < 53.14680163254601 + 79.61307254040804 * pow((mom_P - 0.3), 0.05739232362022314)) {
      return -0.16742969 + 0.00697925 * theta_P;
    //   Dθ = − 0.16742969 + 0.00697925 × θ
    } else
      return 2.04334532 * 10 - 1.81052405 * theta_P + 5.32556360e-2 * theta_P * theta_P -
             5.23157558e-4 * theta_P * theta_P * theta_P;
    //  Dθ = 2.04334532 × 10 − 1.81052405 × θ + 5.32556360 × 10−2 × θ2 − 5.23157558 × 10−4 × θ3
  } else
    return -1.09849291e2 + 8.86664014 * theta_P - 0.26643881 * theta_P * theta_P +
           3.53814210e-3 * theta_P * theta_P * theta_P - 1.75297107e-5 * theta_P * theta_P * theta_P * theta_P;

  //    Aθ = − 1.09849291 × 102 + 8.86664014 × θ − 0.26643881 × θ2 + 3.53814210 × 10−3 ∗ θ3 − 1.75297107 × 10−5 × θ4
}

float B_th(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec) {
  if (dc_sec >= 1 && dc_sec <= 6) {
    if (theta_DCr1_p < 53.14680163254601 + 79.61307254040804 * pow((mom_P - 0.3), 0.05739232362022314)) {
      return 0.23352115 - 0.01338697 * theta_P;
      //  Eθ = 0.23352115 − 0.01338697 × θ
    } else
      return 8.74233279 - 7.63869344e-1 * theta_P + 2.22376362e-2 * theta_P * theta_P -
             2.16457260e-4 * theta_P * theta_P * theta_P;
    // Eθ = 8.74233279 − 7.63869344 × 10−1 × θ + 2.22376362 × 10−2 × θ2 − 2.16457260 × 10−4 × θ3

  } else
    return 9.52034523e2 - 5.74808292e1 * theta_P + 1.15386949 * theta_P * theta_P -
           7.57970373e-3 * theta_P * theta_P * theta_P;
  //      Bθ = 9.52034523 × 102 − 5.74808292 × 10 × θ +1.15386949 × θ2 − 7.57970373 × 10−3 × θ3
}
float C_th(float mom_P, float theta_P, int dc_sec) {
  if (dc_sec < 1 || dc_sec > 6) {
    return -2.00387313e2 + 1.18979079e1 * theta_P - 2.37730217e-1 * theta_P * theta_P +
           1.55153003e-3 * theta_P * theta_P * theta_P;
    // Cθ = − 2.00387313 × 102 + 1.18979079 × 10 × θ − 2.37730217 × 10−1 × θ2 + 1.55153003 × 10−3 × θ3

  } else
    return NAN;
}
// energy loss corrections parameters for phi angle of proton

float A_ph(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec) {
  if (dc_sec >= 1 && dc_sec <= 6) {
    if (theta_DCr1_p < 53.14680163254601 + 79.61307254040804 * pow((mom_P - 0.3), 0.05739232362022314)) {
      return 0.21192125 - 0.0115175 * theta_P;
      //   Dφ = 0.21192125 − 0.0115175 × θ
    } else
      return 0.54697831 - 0.04896981 * theta_P + 0.00111376 * theta_P * theta_P;
    // Aφ = 0.54697831 − 0.04896981 × θ + 0.00111376 × θ2

  } else
    return 4.94546178 - 3.26662886e-1 * theta_P + 7.39069603e-3 * theta_P * theta_P -
           6.83599356e-5 * theta_P * theta_P * theta_P + 2.12303103e-7 * theta_P * theta_P * theta_P * theta_P;
  //    Aφ = 4.94546178 − 3.26662886 × 10−1 × θ + 7.39069603 × 10−3 × θ2 − 6.83599356 × 10−5 × θ3 +
  //   2.12303103 × 10−7 × θ4;
}

float B_ph(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec) {
  if (dc_sec >= 1 && dc_sec <= 6) {
    if (theta_DCr1_p < 53.14680163254601 + 79.61307254040804 * pow((mom_P - 0.3), 0.05739232362022314)) {
      return -8.94307411e-1 + 1.66349766e-1 * theta_P - 8.90617559e-3 * theta_P * theta_P +
             1.64803754e-4 * theta_P * theta_P * theta_P;
      //   Eφ = − 8.94307411 × 10−1 + 1.66349766 × 10−1 × θ − 8.90617559 × 10−3 × θ2 + 1.64803754 × 10−4 × θ3
    } else
      return -4.06733541e2 + 2.43696202e1 * theta_P - 3.36144736e-1 * theta_P * theta_P;
    // Bφ = − 4.06733541 × 102 + 2.43696202 × 10 × θ − 3.36144736 × 10−1 × θ2;
  } else
    return 1.72181613e5 - 1.36827111e4 * theta_P + 4.00923146e2 * theta_P * theta_P -
           5.12792347 * theta_P * theta_P * theta_P + 2.41793167e-2 * theta_P * theta_P * theta_P * theta_P;
  //   Bφ = 1.72181613 × 105 − 1.36827111 × 104 × θ + 4.00923146 × 102 × θ2 − 5.12792347 × θ3 + 2.41793167 × 10−2 ×
  //   θ4;
}
float C_ph(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec) {
  if (dc_sec >= 1 && dc_sec <= 6) {
    return 2.06378660e1 - 1.42866062 * theta_P + 2.01085440e-2 * theta_P * theta_P;
    //    Cφ = 2.06378660 × 10 − 1.42866062 × θ + 2.01085440 × 10−2 × θ2;

  } else
    return 1.20477219e2 - 5.86630228 * theta_P + 7.44007875e-2 * theta_P * theta_P -
           2.42652473e-4 * theta_P * theta_P * theta_P;
  // Cφ = 1.20477219 × 102 − 5.86630228 × θ + 7.44007875 × 10−2 × θ2 − 2.42652473 × 10−4 × θ3;
}

}  // namespace mom_corr
