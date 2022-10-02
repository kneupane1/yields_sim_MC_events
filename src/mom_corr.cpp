
#include "mom_corr.hpp"
#include "iostream"
namespace mom_corr {

bool is_FD(int prot_status) {
  // if (dc_sec >= 1 && dc_sec <= 6)
  if (prot_status > 2000 && prot_status <= 4000)
    return true;
  else
    return false;
}

bool is_CD(int prot_status) {
  // if (dc_sec < 1 || dc_sec > 6)
  if (prot_status > 4000 && prot_status <= 6000)
    return true;
  else
    return false;
}
bool is_lower_band(float mom_P, float theta_DCr1_p, int prot_status) {
  // if (dc_sec >= 1 && dc_sec <= 6) {
  if (prot_status > 2000 && prot_status <= 4000) {
    if (theta_DCr1_p < -53.14680163254601 + 79.61307254040804 * pow((mom_P - 0.3), 0.05739232362022314)) {
      return true;
    } else
      return false;
  } else
    return false;
}
float CD_prot_Emom_corr(float mom_P, float theta_P) {
  return mom_P +
         ((-4.81194246e-05) * pow(mom_P, 3) + 2.14028275e-04 * pow(mom_P, 2) + (-2.57104043e-04) * mom_P +
          1.02579973e-04) *
             pow(theta_P, 3) +
         (0.00595756 * pow(mom_P, 3) + (-0.02653457) * pow(mom_P, 2) + 0.03182286 * mom_P + (-0.0127522)) *
             pow(theta_P, 2) +
         ((-0.24075865) * pow(mom_P, 3) + 1.07424972 * pow(mom_P, 2) + (-1.28641337) * mom_P + 0.51823688) *
             pow(theta_P, 1) +
         3.18175483 * pow(mom_P, 3) + (-14.22566829) * pow(mom_P, 2) + 16.9859584 * mom_P + (-6.88745671);
}

float FD_prot_Emom_corr_lower(float mom_P, float theta_P) {
  return mom_P +
         (2.41366148e-08 * pow(mom_P, 3) + (-8.48694710e-08) * pow(mom_P, 2) + 2.12520490e-08 * mom_P +
          8.19171862e-11) *
             pow(theta_P, 4) +
         ((-1.79468233e-06) * pow(mom_P, 3) + 6.63527873e-06 * pow(mom_P, 2) + (2.41674379e-06) * mom_P +
          1.93217562e-06) *
             pow(theta_P, 3) +
         (4.60815923e-05 * pow(mom_P, 3) + (-1.84383312e-04) * pow(mom_P, 2) + 1.05318538e-04 * mom_P +
          (-1.15779782e-04)) *
             pow(theta_P, 2) +
         ((-0.00049214) * pow(mom_P, 3) + 0.0022003 * pow(mom_P, 2) + (-0.001929) * mom_P + 0.00218473) * theta_P +
         0.00154294 * pow(mom_P, 3) + (-0.00661294) * pow(mom_P, 2) + 0.00329457 * mom_P + (-0.00185376);
}
float FD_prot_Emom_corr_upper(float mom_P, float theta_P) {
  return mom_P +
         ((-6.33926614e-05) * pow(mom_P, 3) + 3.21255513e-04 * pow(mom_P, 2) + (-4.80918164e-04) * mom_P +
          1.94036549e-04) *
             pow(theta_P, 2) +
         (0.00385508 * pow(mom_P, 3) + (-0.0193179) * pow(mom_P, 2) + 0.0279666 * mom_P + (-0.01032478)) *
             pow(theta_P, 1) +
         (-0.06010495) * pow(mom_P, 3) + 0.30123952 * pow(mom_P, 2) + (-0.43371747) * mom_P + 0.16664826;
}

float CD_prot_Eth_corr(float mom_P, float theta_P) {
  return theta_P +
         (0.01794123 * pow(mom_P, 3) + (-0.09198341) * pow(mom_P, 2) + 0.15148531 * mom_P + (-0.0941657)) *
             pow(theta_P, 1) +
         (-0.7392232) * pow(mom_P, 3) + 3.93194154 * pow(mom_P, 2) + (-6.83838677) * mom_P + 4.5505975;
}

float FD_prot_Eth_corr_lower(float mom_P, float theta_P) {
  return theta_P +
         (2.14391671e-05 * pow(mom_P, 3) + (-1.69415274e-04) * pow(mom_P, 2) + 3.62193361e-04 * mom_P +
          (-1.72672065e-04)) *
             pow(theta_P, 2) +
         ((-0.00014124) * pow(mom_P, 3) + 0.00017366 * pow(mom_P, 2) + 0.00466645 * mom_P + (-0.0111939)) *
             pow(theta_P, 1) +
         (-0.00031486) * pow(mom_P, 3) + 0.00897261 * pow(mom_P, 2) + (-0.05371869) * mom_P + 0.08065691;
}

float FD_prot_Eth_corr_upper(float mom_P, float theta_P) {
  return theta_P +
         (0.00165645 * pow(mom_P, 3) + (-0.00983809) * pow(mom_P, 2) + 0.01821203 * mom_P + (-0.01069836)) *
             pow(theta_P, 2) +
         ((-0.10409645) * pow(mom_P, 3) + 0.61354318 * pow(mom_P, 2) + (-1.12258434) * mom_P + 0.64393271) *
             pow(theta_P, 1) +
         1.66090372 * pow(mom_P, 3) + (-9.75714605) * pow(mom_P, 2) + 17.77247321 * mom_P + (-10.0865238);
}

float CD_prot_Eph_corr(float mom_P, float theta_P, float phi_P) {
  return phi_P +
         (0.0152672 * pow(mom_P, 3) + (-0.07306141) * pow(mom_P, 2) + 0.09932124 * mom_P + (-0.04428166)) *
             pow(theta_P, 1) +
         (-0.71565591) * pow(mom_P, 3) + 3.37273717 * pow(mom_P, 2) + (-4.54191832) * mom_P + 1.87540743;
}
float FD_prot_Eph_corr_lower(float mom_P, float theta_P, float phi_P) {
  return phi_P +
         ((-4.86422409e-05) * pow(mom_P, 4) + 1.21216530e-03 * pow(mom_P, 3) + (-8.15266042e-03) * pow(mom_P, 2) +
          1.93258907e-02 * mom_P + (-1.28009681e-02)) *
             pow(theta_P, 1) +
         0.01081378 * pow(mom_P, 4) + (-0.14401558) * pow(mom_P, 3) + 0.69173611 * pow(mom_P, 2) +
         (-1.3964496) * mom_P + 0.95058901;
}

float FD_prot_Eph_corr_upper(float mom_P, float theta_P, float phi_P) {
  return phi_P +
         ((-0.01255713) * pow(mom_P, 3) + 0.07022673 * pow(mom_P, 2) + (-0.12047137) * mom_P + 0.06254443) *
             pow(theta_P, 1) +
         0.27588214 * pow(mom_P, 3) + (-1.37114604) * pow(mom_P, 2) + 1.82000373 * mom_P + (-0.40190107);
}
// // energy loss corrections parameters for momentum of proton
// float A_p(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec) {
//   if (dc_sec >= 1 && dc_sec <= 6) {
//     if (theta_DCr1_p < -53.14680163254601 + 79.61307254040804 * pow((mom_P - 0.3), 0.05739232362022314)) {
//       return -0.00051894 - 0.00018104 * theta_P;
//       //   Ap = − 0.00051894 − 0.00018104 × θ
//       // CorrectedPp_FD_1 = np.select([df_protonRecFD_1.Pp<1, df_protonRecFD_1.Pp>=1], [const_FD +
//       // coeff_FD/df_protonRecFD_1.loc[:, "Pp"] + df_protonRecFD_1.loc[:, "Pp"], np.exp(-2.739
//       // - 3.932*df_protonRecFD_1.Pp) + 0.002907+df_protonRecFD_1.Pp])
//       // np.exp(-2.739 - 3.932*df_protonRecFD_1.Pp) + 0.002907+df_protonRecFD_1.Pp])
//     } else return -3.03346359e-1 + 1.83368163e-2 * theta_P - 2.86486404e-4 * theta_P * theta_P;
//     //   Ap = − 3.03346359 × 10−1 + 1.83368163 × 10−2 × θ − 2.86486404 × 10−4 × θ2
//     // CorrectedPp_FD_2 = np.select([df_protonRecFD_2.Pp<1, df_protonRecFD_2.Pp>=1], [const_FD +
//     // coeff_FD/df_protonRecFD_2.loc[:, "Pp"] + df_protonRecFD_2.loc[:, "Pp"],
//     //  np.exp(-1.2 - 4.228*df_protonRecFD_2.Pp) + 0.007502+df_protonRecFD_2.Pp])

//   } else
// return  1.93686914 - 0.116288824 * theta_P + 0.00223685833 * theta_P * theta_P -
//              1.40771969e-5 * theta_P * theta_P * theta_P;
//   // Ap =1.93686914 − 0.116288824 × θ + 0.00223685833 × θ2 − 1.40771969 × 10−5 × θ3
// }

// float B_p(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec) {
//   if (dc_sec >= 1 && dc_sec <= 6) {
//     if (theta_DCr1_p < -53.14680163254601 + 79.61307254040804 * pow((mom_P - 0.3), 0.05739232362022314)) {
//       return 3.29466917e-3 + 5.73663160e-4 * theta_P - 1.40807209e-5 * theta_P * theta_P;
//       //   Bp =3.29466917 × 10−3 + 5.73663160 × 10−4 × θ − 1.40807209 × 10−5 × θ2.
//     } else
//       return 2.01023276e-1 - 1.13312215e-2 * theta_P + 1.82487916e-4 * theta_P * theta_P;
//     // Bp = 2.01023276 × 10−1 − 1.13312215 × 10−2 × θ + 1.82487916 × 10−4 × θ2.
//   } else
//     return -0.738047800 + 0.0443343685 * theta_P - 8.50985972e-4 * theta_P * theta_P +
//            5.36810280e-6 * theta_P * theta_P * theta_P;
//   //   Bp = − 0.738047800 + 0.0443343685 × θ − 8.50985972 × 10−4 × θ2 + 5.36810280 × 10−6 × θ3
// }

// // energy loss corrections parameters for theta angle of proton

// float A_th(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec) {
//   if (dc_sec >= 1 && dc_sec <= 6) {
//     if (theta_DCr1_p < -53.14680163254601 + 79.61307254040804 * pow((mom_P - 0.3), 0.05739232362022314)) {
//       return -0.16742969 + 0.00697925 * theta_P;
//       //   Dθ = − 0.16742969 + 0.00697925 × θ
//     } else
//       return 2.04334532 * 10 - 1.81052405 * theta_P + 5.32556360e-2 * theta_P * theta_P -
//              5.23157558e-4 * theta_P * theta_P * theta_P;
//     //  Dθ = 2.04334532 × 10 − 1.81052405 × θ + 5.32556360 × 10−2 × θ2 − 5.23157558 × 10−4 × θ3
//   } else
//     return -1.09849291e2 + 8.86664014 * theta_P - 0.26643881 * theta_P * theta_P +
//            3.53814210e-3 * theta_P * theta_P * theta_P - 1.75297107e-5 * theta_P * theta_P * theta_P * theta_P;

//   //    Aθ = − 1.09849291 × 102 + 8.86664014 × θ − 0.26643881 × θ2 + 3.53814210 × 10−3 ∗ θ3 − 1.75297107 × 10−5 ×
//   θ4
// }

// float B_th(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec) {
//   if (dc_sec >= 1 && dc_sec <= 6) {
//     if (theta_DCr1_p < -53.14680163254601 + 79.61307254040804 * pow((mom_P - 0.3), 0.05739232362022314)) {
//       return 0.23352115 - 0.01338697 * theta_P;
//       //  Eθ = 0.23352115 − 0.01338697 × θ
//     } else
//       return 8.74233279 - 7.63869344e-1 * theta_P + 2.22376362e-2 * theta_P * theta_P -
//              2.16457260e-4 * theta_P * theta_P * theta_P;
//     // Eθ = 8.74233279 − 7.63869344 × 10−1 × θ + 2.22376362 × 10−2 × θ2 − 2.16457260 × 10−4 × θ3

//   } else
//     return 9.52034523e2 - 5.74808292e1 * theta_P + 1.15386949 * theta_P * theta_P -
//            7.57970373e-3 * theta_P * theta_P * theta_P;
//   //      Bθ = 9.52034523 × 102 − 5.74808292 × 10 × θ +1.15386949 × θ2 − 7.57970373 × 10−3 × θ3
// }
// float C_th(float mom_P, float theta_P, int dc_sec) {
//   if (dc_sec < 1 || dc_sec > 6) {
//     return -2.00387313e2 + 1.18979079e1 * theta_P - 2.37730217e-1 * theta_P * theta_P +
//            1.55153003e-3 * theta_P * theta_P * theta_P;
//     // Cθ = − 2.00387313 × 102 + 1.18979079 × 10 × θ − 2.37730217 × 10−1 × θ2 + 1.55153003 × 10−3 × θ3

//   } else
//     return NAN;
// }
// // energy loss corrections parameters for phi angle of proton

// float A_ph(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec) {
//   if (dc_sec >= 1 && dc_sec <= 6) {
//     if (theta_DCr1_p < -53.14680163254601 + 79.61307254040804 * pow((mom_P - 0.3), 0.05739232362022314)) {
//       return 0.21192125 - 0.0115175 * theta_P;
//       //   Dφ = 0.21192125 − 0.0115175 × θ
//     } else
//       return 0.54697831 - 0.04896981 * theta_P + 0.00111376 * theta_P * theta_P;
//     // Aφ = 0.54697831 − 0.04896981 × θ + 0.00111376 × θ2

//   } else
//     return 4.94546178 - 3.26662886e-1 * theta_P + 7.39069603e-3 * theta_P * theta_P -
//            6.83599356e-5 * theta_P * theta_P * theta_P + 2.12303103e-7 * theta_P * theta_P * theta_P * theta_P;
//   //    Aφ = 4.94546178 − 3.26662886 × 10−1 × θ + 7.39069603 × 10−3 × θ2 − 6.83599356 × 10−5 × θ3 +
//   //   2.12303103 × 10−7 × θ4;
// }

// float B_ph(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec) {
//   if (dc_sec >= 1 && dc_sec <= 6) {
//     if (theta_DCr1_p < -53.14680163254601 + 79.61307254040804 * pow((mom_P - 0.3), 0.05739232362022314)) {
//       return -8.94307411e-1 + 1.66349766e-1 * theta_P - 8.90617559e-3 * theta_P * theta_P +
//              1.64803754e-4 * theta_P * theta_P * theta_P;
//       //   Eφ = − 8.94307411 × 10−1 + 1.66349766 × 10−1 × θ − 8.90617559 × 10−3 × θ2 + 1.64803754 × 10−4 × θ3
//     } else
//       return -4.06733541e2 + 2.43696202e1 * theta_P - 3.36144736e-1 * theta_P * theta_P;
//     // Bφ = − 4.06733541 × 102 + 2.43696202 × 10 × θ − 3.36144736 × 10−1 × θ2;
//   } else
//     return 1.72181613e5 - 1.36827111e4 * theta_P + 4.00923146e2 * theta_P * theta_P -
//            5.12792347 * theta_P * theta_P * theta_P + 2.41793167e-2 * theta_P * theta_P * theta_P * theta_P;
//   //   Bφ = 1.72181613 × 105 − 1.36827111 × 104 × θ + 4.00923146 × 102 × θ2 − 5.12792347 × θ3 + 2.41793167 × 10−2 ×
//   //   θ4;
// }
// float C_ph(float mom_P, float theta_P, float theta_DCr1_p, int dc_sec) {
//   if (dc_sec >= 1 && dc_sec <= 6) {
//     return 2.06378660e1 - 1.42866062 * theta_P + 2.01085440e-2 * theta_P * theta_P;
//     //    Cφ = 2.06378660 × 10 − 1.42866062 × θ + 2.01085440 × 10−2 × θ2;

//   } else
//     return 1.20477219e2 - 5.86630228 * theta_P + 7.44007875e-2 * theta_P * theta_P -
//            2.42652473e-4 * theta_P * theta_P * theta_P;
//   // Cφ = 1.20477219 × 102 − 5.86630228 × θ + 7.44007875 × 10−2 × θ2 − 2.42652473 × 10−4 × θ3;
// }

}  // namespace mom_corr
