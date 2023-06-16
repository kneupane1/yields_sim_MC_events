
#include "mom_corr.hpp"
#include <cstdlib>
#include <ctime>
#include "iostream"

mom_corr::~mom_corr() {}

bool mom_corr::is_FD(int part_status) {
  // if (dc_sec >= 1 && dc_sec <= 6)
  if (part_status > 2000 && part_status <= 4000)
    return true;
  else
    return false;
}

// bool mom_corr::is_AllFD(int part1_status, int part2_status, int part3_status) {
//   // if (dc_sec >= 1 && dc_sec <= 6)
//   if ((part1_status > 2000 && part1_status <= 4000) && (part2_status > 2000 && part2_status <= 4000) &&
//       (part3_status > 2000 && part3_status <= 4000))
//     return true;
//   else
//     return false;
// }

bool mom_corr::is_CD(int part_status) {
  // if (dc_sec < 1 || dc_sec > 6)
  if (part_status > 4000 && part_status <= 6000)
    return true;
  else
    return false;
}
bool mom_corr::is_lower_band(float mom_, float theta_DCr1_, int status_) {
  // if (dc_sec >= 1 && dc_sec <= 6) {
  if (status_ > 2000 && status_ <= 4000) {
    if (theta_DCr1_ < -53.14680163254601 + 79.61307254040804 * pow((mom_ - 0.3), 0.05739232362022314)) {
      return true;
    } else
      return false;
  } else
    return false;
}

float mom_corr::CD_prot_Emom_corr(float mom_, float theta_) {
  return mom_ +
         ((-4.81194246e-05) * pow(mom_, 3) + 2.14028275e-04 * pow(mom_, 2) + (-2.57104043e-04) * mom_ +
          1.02579973e-04) *
             pow(theta_, 3) +
         (0.00595756 * pow(mom_, 3) + (-0.02653457) * pow(mom_, 2) + 0.03182286 * mom_ + (-0.0127522)) *
             pow(theta_, 2) +
         ((-0.24075865) * pow(mom_, 3) + 1.07424972 * pow(mom_, 2) + (-1.28641337) * mom_ + 0.51823688) *
             pow(theta_, 1) +
         3.18175483 * pow(mom_, 3) + (-14.22566829) * pow(mom_, 2) + 16.9859584 * mom_ + (-6.88745671);
}

float mom_corr::FD_prot_Emom_corr_lower(float mom_, float theta_) {
  return mom_ +
         (2.41366148e-08 * pow(mom_, 3) + (-8.48694710e-08) * pow(mom_, 2) + 2.12520490e-08 * mom_ + 8.19171862e-11) *
             pow(theta_, 4) +
         ((-1.79468233e-06) * pow(mom_, 3) + 6.63527873e-06 * pow(mom_, 2) + (-2.41674379e-06) * mom_ +
          1.93217562e-06) *
             pow(theta_, 3) +
         (4.60815923e-05 * pow(mom_, 3) + (-1.84383312e-04) * pow(mom_, 2) + 1.05318538e-04 * mom_ +
          (-1.15779782e-04)) *
             pow(theta_, 2) +
         ((-0.00049214) * pow(mom_, 3) + 0.0022003 * pow(mom_, 2) + (-0.001929) * mom_ + 0.00218473) * theta_ +
         0.00154294 * pow(mom_, 3) + (-0.00661294) * pow(mom_, 2) + 0.00329457 * mom_ + (-0.00185376);
}
float mom_corr::FD_prot_Emom_corr_upper(float mom_, float theta_) {
  return mom_ +
         ((-6.33926614e-05) * pow(mom_, 3) + 3.21255513e-04 * pow(mom_, 2) + (-4.80918164e-04) * mom_ +
          1.94036549e-04) *
             pow(theta_, 2) +
         (0.00385508 * pow(mom_, 3) + (-0.0193179) * pow(mom_, 2) + 0.0279666 * mom_ + (-0.01032478)) * pow(theta_, 1) +
         (-0.06010495) * pow(mom_, 3) + 0.30123952 * pow(mom_, 2) + (-0.43371747) * mom_ + 0.16664826;
}

float mom_corr::CD_prot_Eth_corr(float mom_, float theta_) {
  return theta_ +
         (0.01794123 * pow(mom_, 3) + (-0.09198341) * pow(mom_, 2) + 0.15148531 * mom_ + (-0.0941657)) *
             pow(theta_, 1) +
         (-0.7392232) * pow(mom_, 3) + 3.93194154 * pow(mom_, 2) + (-6.83838677) * mom_ + 4.5505975;
}

float mom_corr::FD_prot_Eth_corr_lower(float mom_, float theta_) {
  return theta_ +
         (2.14391671e-05 * pow(mom_, 3) + (-1.69415274e-04) * pow(mom_, 2) + 3.62193361e-04 * mom_ +
          (-1.72672065e-04)) *
             pow(theta_, 2) +
         ((-0.00014124) * pow(mom_, 3) + 0.00017366 * pow(mom_, 2) + 0.00466645 * mom_ + (-0.0111939)) *
             pow(theta_, 1) +
         (-0.00031486) * pow(mom_, 3) + 0.00897261 * pow(mom_, 2) + (-0.05371869) * mom_ + 0.08065691;
}

float mom_corr::FD_prot_Eth_corr_upper(float mom_, float theta_) {
  return theta_ +
         (0.00165645 * pow(mom_, 3) + (-0.00983809) * pow(mom_, 2) + 0.01821203 * mom_ + (-0.01069836)) *
             pow(theta_, 2) +
         ((-0.10409645) * pow(mom_, 3) + 0.61354318 * pow(mom_, 2) + (-1.12258434) * mom_ + 0.64393271) *
             pow(theta_, 1) +
         1.66090372 * pow(mom_, 3) + (-9.75714605) * pow(mom_, 2) + 17.77247321 * mom_ + (-10.0865238);
}

float mom_corr::CD_prot_Eph_corr(float mom_, float theta_, float phi_) {
  return phi_ +
         (0.0152672 * pow(mom_, 3) + (-0.07306141) * pow(mom_, 2) + 0.09932124 * mom_ + (-0.04428166)) *
             pow(theta_, 1) +
         (-0.71565591) * pow(mom_, 3) + 3.37273717 * pow(mom_, 2) + (-4.54191832) * mom_ + 1.87540743;
}
float mom_corr::FD_prot_Eph_corr_lower(float mom_, float theta_, float phi_) {
  return phi_ +
         ((-4.86422409e-05) * pow(mom_, 4) + 1.21216530e-03 * pow(mom_, 3) + (-8.15266042e-03) * pow(mom_, 2) +
          1.93258907e-02 * mom_ + (-1.28009681e-02)) *
             pow(theta_, 1) +
         0.01081378 * pow(mom_, 4) + (-0.14401558) * pow(mom_, 3) + 0.69173611 * pow(mom_, 2) + (-1.3964496) * mom_ +
         0.95058901;
}

float mom_corr::FD_prot_Eph_corr_upper(float mom_, float theta_, float phi_) {
  return phi_ +
         ((-0.01255713) * pow(mom_, 3) + 0.07022673 * pow(mom_, 2) + (-0.12047137) * mom_ + 0.06254443) *
             pow(theta_, 1) +
         0.27588214 * pow(mom_, 3) + (-1.37114604) * pow(mom_, 2) + 1.82000373 * mom_ + (-0.40190107);
}
// // energy loss corrections parameters for momentum of proton
float mom_corr::A_p(float mom_, float theta_, float theta_DCr1_p, int dc_sec) {
  // if (dc_sec >= 1 && dc_sec <= 6) {
  if (theta_DCr1_p < -53.14680163254601 + 79.61307254040804 * pow((mom_ - 0.3), 0.05739232362022314)) {
    return -0.00051894 - 0.00018104 * theta_;
    //   Ap = − 0.00051894 − 0.00018104 × θ
    // CorrectedPp_FD_1 = np.select([df_protonRecFD_1.Pp<1, df_protonRecFD_1.Pp>=1], [const_FD +
    // coeff_FD/df_protonRecFD_1.loc[:, "Pp"] + df_protonRecFD_1.loc[:, "Pp"], np.exp(-2.739
    // - 3.932*df_protonRecFD_1.Pp) + 0.002907+df_protonRecFD_1.Pp])
    // np.exp(-2.739 - 3.932*df_protonRecFD_1.Pp) + 0.002907+df_protonRecFD_1.Pp])
  } else
    return -3.03346359e-1 + 1.83368163e-2 * theta_ - 2.86486404e-4 * theta_ * theta_;
  //   Ap = − 3.03346359 × 10−1 + 1.83368163 × 10−2 × θ − 2.86486404 × 10−4 × θ2
  // CorrectedPp_FD_2 = np.select([df_protonRecFD_2.Pp<1, df_protonRecFD_2.Pp>=1], [const_FD +
  // coeff_FD/df_protonRecFD_2.loc[:, "Pp"] + df_protonRecFD_2.loc[:, "Pp"],
  //  np.exp(-1.2 - 4.228*df_protonRecFD_2.Pp) + 0.007502+df_protonRecFD_2.Pp])

  //   } else
  // return  1.93686914 - 0.116288824 * theta_ + 0.00223685833 * theta_ * theta_ -
  //              1.40771969e-5 * theta_ * theta_ * theta_;
  //   // Ap =1.93686914 − 0.116288824 × θ + 0.00223685833 × θ2 − 1.40771969 × 10−5 × θ3
}

float mom_corr::B_p(float mom_, float theta_, float theta_DCr1_p, int dc_sec) {
  // if (dc_sec >= 1 && dc_sec <= 6) {
  if (theta_DCr1_p < -53.14680163254601 + 79.61307254040804 * pow((mom_ - 0.3), 0.05739232362022314)) {
    return 3.29466917e-3 + 5.73663160e-4 * theta_ - 1.40807209e-5 * theta_ * theta_;
    //   Bp =3.29466917 × 10−3 + 5.73663160 × 10−4 × θ − 1.40807209 × 10−5 × θ2.
  } else
    return 2.01023276e-1 - 1.13312215e-2 * theta_ + 1.82487916e-4 * theta_ * theta_;
  // Bp = 2.01023276 × 10−1 − 1.13312215 × 10−2 × θ + 1.82487916 × 10−4 × θ2.
  // } else
  //   return -0.738047800 + 0.0443343685 * theta_ - 8.50985972e-4 * theta_ * theta_ +
  //          5.36810280e-6 * theta_ * theta_ * theta_;
  // //   Bp = − 0.738047800 + 0.0443343685 × θ − 8.50985972 × 10−4 × θ2 + 5.36810280 × 10−6 × θ3
}

// energy loss corrections for pip

float mom_corr::CD_pip_Emom_corr(float mom_, float theta_) {
  return mom_ +
         ((-6.06092449e-07) * pow(theta_, 3) + 1.32660527e-04 * pow(theta_, 2) + (-9.21399702e-03) * theta_ +
          2.30256661e-01) *
             pow(mom_, 3) +
         (1.99184379e-06 * pow(theta_, 3) + (-4.43181568e-04) * pow(theta_, 2) + 3.15039271e-02 * theta_ +
          (-7.97320779e-01)) *
             pow(mom_, 2) +
         ((-2.00127680e-06) * pow(theta_, 3) + 4.61630337e-04 * pow(theta_, 2) + (-3.41672108e-02) * theta_ +
          8.64527869e-01) *
             pow(mom_, 1) +
         4.14468224e-07 * pow(theta_, 3) + (-1.07089463e-04) * pow(theta_, 2) + 9.25833758e-03 * theta_ +
         (-2.74924349e-01);
}
float mom_corr::FD_pip_Emom_corr_lower(float mom_, float theta_) {
  return mom_ + (-4.67842670e-05) * pow(mom_, 3) + 3.37133020e-04 * pow(mom_, 2) + (-4.79135831e-04) * mom_ +
         2.70872474e-03;
}
float mom_corr::FD_pip_Emom_corr_upper(float mom_, float theta_) {
  return mom_ + (-0.00125149) * pow(mom_, 3) + 0.0053441 * pow(mom_, 2) + (-0.00765213) * mom_ + 0.0102172;
}

float mom_corr::CD_pip_Eth_corr(float mom_, float theta_) {
  if (mom_ <= 0.7) {
    return theta_ +
           (1.50263076e-06 * pow(mom_, 3) + (-4.71834964e-06) * pow(mom_, 2) + 4.19603178e-06 * mom_ +
            (-1.22889036e-06)) *
               pow(theta_, 4) +
           ((-0.00042763) * pow(mom_, 3) + 0.00134022 * pow(mom_, 2) + (-0.00118851) * mom_ + 0.00034294) *
               pow(theta_, 3) +

           (0.04191854 * pow(mom_, 3) + (-0.13037561) * pow(mom_, 2) + 0.11407653 * mom_ + (-0.03178403)) *
               pow(theta_, 2) +
           ((-1.57945065) * pow(mom_, 3) + 4.77845697 * pow(mom_, 2) + (-3.96720052) * mom_ + 0.98986696) *
               pow(theta_, 1) +

           14.28409289 * pow(mom_, 3) + (-37.03568066) * pow(mom_, 2) + 20.96721711 * mom_ + (-0.3402565);
  } else {
    return theta_ + (-0.07926959493130192) * mom_ + 0.29484361324796154;
  }
}
float mom_corr::FD_pip_Eth_corr_lower(float mom_, float theta_) {
  return theta_ +
         (5.82345268e-07 * pow(mom_, 4) + (-6.50577207e-06) * pow(mom_, 3) + 2.69047970e-05 * pow(mom_, 2) +
          (-4.63578237e-05) * pow(mom_, 1) + 2.92063857e-05) *
             pow(theta_, 3) +
         ((-1.70152392e-05) * pow(mom_, 4) + 2.08992182e-04 * pow(mom_, 3) + (-9.71300032e-04) * pow(mom_, 2) +
          1.81681161e-03 * pow(mom_, 1) + (-1.22931209e-03)) *
             pow(theta_, 2) +

         ((-0.0004973) * pow(mom_, 4) + 0.00534567 * pow(mom_, 3) + (-0.01984666) * pow(mom_, 2) +
          0.03332743 * pow(mom_, 1) + (-0.02142081)) *
             pow(theta_, 1) +
         0.00841078 * pow(mom_, 4) + (-0.09350417) * pow(mom_, 3) + 0.36576903 * pow(mom_, 2) +
         (-0.6074988) * pow(mom_, 1) + 0.35290183;
}

float mom_corr::FD_pip_Eth_corr_upper(float mom_, float theta_) {
  return theta_ +
         (0.00094724 * pow(mom_, 3) + (-0.00524101) * pow(mom_, 2) + 0.00919525 * mom_ + (-0.00516691)) *
             pow(theta_, 3) +
         ((-0.09887756) * pow(mom_, 3) + 0.54682169 * pow(mom_, 2) + (-0.95634115) * mom_ + 0.53345618) *
             pow(theta_, 2) +
         (3.44104365 * pow(mom_, 3) + (-19.02680178) * pow(mom_, 2) + 33.17864748 * mom_ + (-18.37813421)) *
             pow(theta_, 1) +
         (-39.82151866) * pow(mom_, 3) + 220.12521819 * pow(mom_, 2) + (-382.61957089) * mom_ + 210.34677439;
}

float mom_corr::CD_pip_Eph_corr(float mom_, float theta_, float phi_) {
  if (mom_ <= 0.7) {
    return phi_ +
           ((-5.02775972e-07) * pow(mom_, 3) + 1.77952733e-06 * pow(mom_, 2) + (-1.91537716e-06) * mom_ +
            8.14069464e-07) *
               pow(theta_, 4) +
           (0.00015302 * pow(mom_, 3) + (-0.00054583) * pow(mom_, 2) + 0.00059431 * mom_ + (-0.00025352)) *
               pow(theta_, 3) +
           ((-0.01619882) * pow(mom_, 3) + 0.05826388 * pow(mom_, 2) + (-0.06424007) * mom_ + 0.02763694) *
               pow(theta_, 2) +
           (0.67677027 * pow(mom_, 3) + (-2.45354146) * pow(mom_, 2) + 2.7393312 * mom_ + (-1.20043622)) *
               pow(theta_, 1) +
           (-8.07766719) * pow(mom_, 3) + 29.66313429 * pow(mom_, 2) + (-33.9669606) * mom_ + 15.78966364;
  } else {
    return phi_ + 0.04826653377945466 * mom_ + (-0.21426965774563544);
  }
}
float mom_corr::FD_pip_Eph_corr_lower(float mom_, float theta_, float phi_) {
  return phi_ +
         ((-4.86422409e-05) * pow(mom_, 4) + 1.21216530e-03 * pow(mom_, 3) + (-8.15266042e-03) * pow(mom_, 2) +
          1.93258907e-02 * mom_ + (-1.28009681e-02)) *
             pow(theta_, 1) +
         0.01081378 * pow(mom_, 4) + (-0.14401558) * pow(mom_, 3) + 0.69173611 * pow(mom_, 2) + (-1.3964496) * mom_ +
         0.95058901;
}
float mom_corr::FD_pip_Eph_corr_upper(float mom_, float theta_, float phi_) {
  return phi_ +
         ((-0.02343664) * pow(mom_, 3) + 0.13264734 * pow(mom_, 2) + (-0.2342437) * mom_ + 0.12601401) *
             pow(theta_, 1) +
         0.50037573 * pow(mom_, 3) + (-2.72628993) * pow(mom_, 2) + 4.48508987 * mom_ + (-2.05446324);
}

// energy loss corrections for pim

float mom_corr::CD_pim_Emom_corr(float mom_, float theta_) {
  // return mom_ + ((-1.66077208e-08) * pow(mom_, 2) + 5.87672135e-08 * mom_ + (-1.35413089e-08)) * pow(theta_, 4) +
  //        (5.15167601e-06 * pow(mom_, 2) + (-1.79444621e-05) * mom_ + 4.06971096e-06) * pow(theta_, 3) +
  //        ((-0.00057812) * pow(mom_, 2) + 0.00197867 * mom_ + (-0.00044994)) * pow(theta_, 2) +
  //        (0.02778557 * pow(mom_, 2) + (-0.09352583) * mom_ + 0.02226586) * pow(theta_, 1) +
  //        (-0.47794319) * pow(mom_, 2) + 1.57678098 * mom_ + (-0.41789067);

  if (theta_ <= 90) {
    return mom_ +
           (-4.94426765e-07 * pow(theta_, 3) + 9.85729368e-05 * pow(theta_, 2) + (-5.85778699e-03) * (theta_) +
            1.17447168e-01) *
               pow(mom_, 3) +

           (1.75953956e-06 * pow(theta_, 3) + (-3.63382515e-04) * pow(theta_, 2) + 2.21447425e-02 * (theta_) +
            (-4.54844509e-01)) *
               pow(mom_, 2) +

           (-1.90446515e-06 * pow(theta_, 3) + 4.08768480e-04 * pow(theta_, 2) + (-2.65277055e-02) * (theta_) +
            5.57286393e-01) *
               (mom_) +

           2.05653097e-07 * pow(theta_, 3) + (-5.44018546e-05) * pow(theta_, 2) +
           4.61561853e-03 * (theta_)-1.35303212e-01;
  } else {
    return mom_ + 2.27546950e-07 * pow(theta_, 3) + (-8.12537308e-05) * pow(theta_, 2) +
           9.10902744e-03 * pow(theta_, 1) + (-3.22464750e-01);
  }
}
float mom_corr::FD_pim_Emom_corr_lower(float mom_, float theta_) { return mom_ + 0.00030448 * mom_ + 0.00232071; }
float mom_corr::FD_pim_Emom_corr_upper(float mom_, float theta_) { return mom_ + (-0.00100881) * mom_ + 0.00780439; }

float mom_corr::CD_pim_Eth_corr(float mom_, float theta_) {
  if (mom_ <= 0.7) {
    return theta_ +
           (7.39231883e-06 * pow(mom_, 3) + (-1.50802473e-05) * pow(mom_, 2) + 9.79813939e-06 * mom_ +
            (-2.16012840e-06)) *
               pow(theta_, 4) +
           ((-0.00222313) * pow(mom_, 3) + 0.00452095 * pow(mom_, 2) + (-0.00291633) * mom_ + 0.00063192) *
               pow(theta_, 3) +

           (0.23465449 * pow(mom_, 3) + (-0.47381883) * pow(mom_, 2) + 0.30115767 * mom_ + (-0.06318632)) *
               pow(theta_, 2) +
           ((-9.96392109) * pow(mom_, 3) + 19.7772132 * pow(mom_, 2) + (-12.14375333) * mom_ + 2.36811829) *
               pow(theta_, 1) +
           130.65299881 * pow(mom_, 3) + (-246.22737915) * pow(mom_, 2) + 135.30865002 * mom_ + (-19.89993903);
  } else {
    return theta_ + (-0.10181687) * mom_ + 0.28868377;
  }
}
float mom_corr::FD_pim_Eth_corr_lower(float mom_, float theta_) {
  return theta_ + ((-1.13685553e-04) * pow(mom_, 4) + 4.19458440e-03 * pow(mom_, 3) + (-3.76566663e-02) * pow(mom_, 2) +
                   1.30733557e-01 * pow(mom_, 1) + (-1.76073418e-01));
}

float mom_corr::FD_pim_Eth_corr_upper(float mom_, float theta_) {
  return theta_ +
         (0.01520214 * pow(mom_, 3) + (-0.08264195) * pow(mom_, 2) + 0.14545703 * mom_ + (-0.0888854)) *
             pow(theta_, 1) +
         (-0.46222418) * pow(mom_, 3) + 2.45741975 * pow(mom_, 2) + (-4.17396135) * mom_ + 2.39541974;
}

float mom_corr::CD_pim_Eph_corr(float mom_, float theta_, float phi_) {
  if (mom_ <= 0.7) {
    return phi_ +
           ((-2.40376620e-06) * pow(mom_, 3) + 5.50564834e-06 * pow(mom_, 2) + (-3.61060685e-06) * mom_ +
            3.18869876e-07) *
               pow(theta_, 4) +
           (8.05453348e-04 * pow(mom_, 3) + (-1.79820994e-03) * pow(mom_, 2) + 1.14727340e-03 * mom_ +
            (-9.70646762e-05)) *
               pow(theta_, 3) +
           ((-0.10033497) * pow(mom_, 3) + 0.21835219 * pow(mom_, 2) + (-0.13602144) * mom_ + 0.01176118) *
               pow(theta_, 2) +
           (5.49119001 * pow(mom_, 3) + (-11.70340094) * pow(mom_, 2) + 7.20324368 * mom_ + (-0.70284274)) *
               pow(theta_, 1) +
           (-110.77154056) * pow(mom_, 3) + 232.99674561 * pow(mom_, 2) + (-143.93904349) * mom_ + 17.1013553;
  } else {
    return phi_ + (-0.08507155) * mom_ + 0.28063752;
  }
}

float mom_corr::FD_pim_Eph_corr_lower(float mom_, float theta_, float phi_) {
  return phi_ +
         ((-2.86749702e-05) * pow(mom_, 4) + 3.03813193e-04 * pow(mom_, 3) + (-1.12379180e-03) * pow(mom_, 2) +
          1.70003187e-03 * pow(mom_, 1) + (-8.77541156e-04)) *
             pow(theta_, 3) +

         (0.00196534 * pow(mom_, 4) + (-0.0209559) * pow(mom_, 3) + 0.07804889 * pow(mom_, 2) +
          (-0.11856395) * pow(mom_, 1) + 0.06067883) *
             pow(theta_, 2) +
         ((-0.04397531) * pow(mom_, 4) + 0.47211422 * pow(mom_, 3) + (-1.76953348) * pow(mom_, 2) +
          2.69302517 * pow(mom_, 1) + (-1.35729049)) *
             pow(theta_, 1) +
         0.32282676 * pow(mom_, 4) + (-3.48574851) * pow(mom_, 3) + 13.11695944 * pow(mom_, 2) +
         (-19.9133663) * pow(mom_, 1) + 9.82183739;
}
float mom_corr::FD_pim_Eph_corr_upper(float mom_, float theta_, float phi_) {
  return phi_ +
         ((-0.00049736) * pow(mom_, 4) + 0.0022372 * pow(mom_, 3) + (-0.00317915) * pow(mom_, 2) +
          0.00218449 * pow(mom_, 1) + (-0.00080044)) *
             pow(theta_, 3) +

         (0.03527214 * pow(mom_, 4) + (-0.13061747) * pow(mom_, 3) + 0.10585015 * pow(mom_, 2) +
          (-0.01732555) * pow(mom_, 1) + 0.013915) *
             pow(theta_, 2) +
         ((-0.80579394) * pow(mom_, 4) + 2.02762174 * pow(mom_, 3) + 1.63055562 * pow(mom_, 2) +
          (-4.17036909) * pow(mom_, 1) + 1.03295098) *
             pow(theta_, 1) +
         6.17963055 * pow(mom_, 4) + (-5.81705813) * pow(mom_, 3) + (-53.39466945) * pow(mom_, 2) +
         77.16020833 * pow(mom_, 1) + (-20.58824011);
}
//////////////////// new mom correction start Aug 2022 Version

double mom_corr::dppC(float Px, float Py, float Pz, int sec, int ivec) {
// auto dppC = [&](float Px, float Py, float Pz, int sec, int ivec) {
// ivec = 0 --> Electron Corrections
// ivec = 1 --> Pi+ Corrections
// ivec = 2 --> Pi- Corrections
// ivec = 3 --> Proton Corrections

// Momentum Magnitude
double pp = sqrt(Px * Px + Py * Py + Pz * Pz);

  // Initializing the correction factor
  double dp = 0;

  // Defining Phi Angle
  double Phi = (180 / 3.1415926) * atan2(Py, Px);

  // (Initial) Shift of the Phi Angle (done to realign sectors whose data is separated when plotted from ±180˚)
  if (((sec == 4 || sec == 3) && Phi < 0) || (sec > 4 && Phi < 90)) {
    Phi += 360;
  }

  // Getting Local Phi Angle
  double PhiLocal = Phi - (sec - 1) * 60;

  // Applying Shift Functions to Phi Angles (local shifted phi = phi)
  double phi = PhiLocal;

  // For Electron Shift
  if (ivec == 0) {
    phi = PhiLocal - 30 / pp;
  }

  // For Pi+ Pion/Proton Shift
  if (ivec == 1 || ivec == 3) {
    phi = PhiLocal + (32 / (pp - 0.05));
  }

  // For Pi- Pion Shift
  if (ivec == 2) {
    phi = PhiLocal - (32 / (pp - 0.05));
  }

  //==========//  PARTICLE = ELECTRON  //==========//

  if (ivec == 0) {
    if (sec == 1) {
      dp = ((1.57e-06) * phi * phi + (5.021e-05) * phi + (-1.74089e-03)) * pp * pp +
           ((-2.192e-05) * phi * phi + (-1.12528e-03) * phi + (0.0146476)) * pp +
           ((8.504e-05) * phi * phi + (2.08012e-03) * phi + (-0.0122501));
    }
    if (sec == 2) {
      dp = ((-3.98e-06) * phi * phi + (1.66e-05) * phi + (-1.55918e-03)) * pp * pp +
           ((2.136e-05) * phi * phi + (-5.7373e-04) * phi + (0.0143591)) * pp +
           ((2.4e-06) * phi * phi + (1.6656e-03) * phi + (-0.0218711));
    }

    if (sec == 3) {
      dp = ((5.57e-06) * phi * phi + (2.3e-07) * phi + (-2.26999e-03)) * pp * pp +
           ((-7.761e-05) * phi * phi + (4.1437e-04) * phi + (0.0152985)) * pp +
           ((2.2542e-04) * phi * phi + (-9.442e-04) * phi + (-0.0231432));
    }

    if (sec == 4) {
      dp = ((3.48e-06) * phi * phi + (2.166e-05) * phi + (-2.29e-04)) * pp * pp +
           ((-2.758e-05) * phi * phi + (7.226e-05) * phi + (-3.38e-03)) * pp +
           ((3.166e-05) * phi * phi + (6.93e-05) * phi + (0.04767));
    }

    if (sec == 5) {
      dp = ((1.19e-06) * phi * phi + (-2.286e-05) * phi + (-1.6332e-04)) * pp * pp +
           ((-1.05e-06) * phi * phi + (7.04e-05) * phi + (-5.0754e-03)) * pp +
           ((-7.22e-06) * phi * phi + (4.1748e-04) * phi + (0.04441));
    }

    if (sec == 6) {
      dp = ((-5.97e-06) * phi * phi + (-3.689e-05) * phi + (5.782e-05)) * pp * pp +
           ((6.573e-05) * phi * phi + (2.1376e-04) * phi + (-9.54576e-03)) * pp +
           ((-1.7732e-04) * phi * phi + (-8.62e-04) * phi + (0.0618975));
    }
  }

  //==========//  PARTICLE = ELECTRON (END)  //==========//

  //==========//  PARTICLE = PI+ PION  //==========//

  if (ivec == 1) {
    if (sec == 1) {
      dp = ((-5.2e-07) * phi * phi + (-1.383e-05) * phi + (4.7179e-04)) * pp * pp +
           ((8.33e-06) * phi * phi + (3.8849e-04) * phi + (-6.81319e-03)) * pp +
           ((-1.645e-05) * phi * phi + (-5.0057e-04) * phi + (1.9902e-02));
    }

    if (sec == 2) {
      dp = ((-1.88e-06) * phi * phi + (3.303e-05) * phi + (1.1331e-03)) * pp * pp +
           ((1.569e-05) * phi * phi + (-3.974e-05) * phi + (-1.25869e-02)) * pp +
           ((-2.903e-05) * phi * phi + (-1.0638e-04) * phi + (2.61529e-02));
    }
    if (sec == 3) {
      dp = ((2.4e-07) * phi * phi + (-1.04e-05) * phi + (7.0864e-04)) * pp * pp +
           ((8.0e-06) * phi * phi + (-5.156e-05) * phi + (-8.12169e-03)) * pp +
           ((-2.42e-05) * phi * phi + (8.928e-05) * phi + (2.13223e-02));
    }
    if (sec == 4) {
      dp = ((-4.0e-08) * phi * phi + (-3.59e-05) * phi + (1.32146e-03)) * pp * pp +
           ((1.023e-05) * phi * phi + (2.2199e-04) * phi + (-1.33043e-02)) * pp +
           ((-2.801e-05) * phi * phi + (-1.576e-04) * phi + (3.27995e-02));
    }
    if (sec == 5) {
      dp = ((2.7e-06) * phi * phi + (5.03e-06) * phi + (1.59668e-03)) * pp * pp +
           ((-1.28e-05) * phi * phi + (-1.99e-06) * phi + (-1.71578e-02)) * pp +
           ((2.091e-05) * phi * phi + (-4.14e-05) * phi + (3.25434e-02));
    }
    if (sec == 6) {
      dp = ((2.13e-06) * phi * phi + (-7.49e-05) * phi + (1.75565e-03)) * pp * pp +
           ((-7.37e-06) * phi * phi + (5.8222e-04) * phi + (-1.27969e-02)) * pp +
           ((4.9e-07) * phi * phi + (-7.2253e-04) * phi + (3.11499e-02));
    }
  }

  //==========//  PARTICLE = PI+ PION (END)  //==========//

  //==========//  PARTICLE = PI- PION  //==========//

  if (ivec == 2) {
    if (sec == 1) {
      dp = ((-4.0192658422317425e-06) * phi * phi - (2.660222128967742e-05) * phi + 0.004774434682983547) * pp * pp;
      dp = dp + ((1.9549520962477972e-05) * phi * phi - 0.0002456062756770577 * phi - 0.03787692408323466) * pp;
      dp = dp + (-2.128953094937459e-05) * phi * phi + 0.0002461708852239913 * phi + 0.08060704449822174 - 0.01;
    }

    if (sec == 2) {
      dp = ((1.193010521758372e-05) * phi * phi - (5.996221756031922e-05) * phi + 0.0009093437955814359) * pp * pp;
      dp = dp + ((-4.89113824430594e-05) * phi * phi + 0.00021676479488147118 * phi - 0.01861892053916726) * pp;
      dp = dp + (4.446394152208071e-05) * phi * phi - (3.6592784167335244e-05) * phi + 0.05498710249944096 - 0.01;
    }

    if (sec == 3) {
      dp = ((-1.6596664895992133e-07) * phi * phi + (6.317189710683516e-05) * phi + 0.0016364212312654086) * pp * pp;
      dp = dp + ((-2.898409777520318e-07) * phi * phi - 0.00014531513577533802 * phi - 0.025456145839203827) * pp;
      dp = dp + (2.6432552410603506e-06) * phi * phi + 0.00018447151306275443 * phi + 0.06442602664627255 - 0.01;
    }

    if (sec == 4) {
      dp = ((2.4035259647558634e-07) * phi * phi - (8.649647351491232e-06) * phi + 0.004558993439848128) * pp * pp;
      dp = dp + ((-5.981498144060984e-06) * phi * phi + 0.00010582131454222416 * phi - 0.033572004651981686) * pp;
      dp = dp + (8.70140266889548e-06) * phi * phi - 0.00020137414379966883 * phi + 0.07258774523336173 - 0.01;
    }

    if (sec == 5) {
      dp = ((2.5817024702834863e-06) * phi * phi + 0.00010132810066914441 * phi + 0.003397314538804711) * pp * pp;
      dp = dp + ((-1.5116941263931812e-05) * phi * phi - 0.00040679799541839254 * phi - 0.028144285760769876) * pp;
      dp = dp + (1.4701931057951464e-05) * phi * phi + 0.0002426350390593454 * phi + 0.06781682510174941 - 0.01;
    }

    if (sec == 6) {
      dp = ((-8.196823669099362e-07) * phi * phi - (5.280412421933636e-05) * phi + 0.0018457238328451137) * pp * pp;
      dp = dp + ((5.2675062282094536e-06) * phi * phi + 0.0001515803461044587 * phi - 0.02294371578470564) * pp;
      dp = dp + (-9.459454671739747e-06) * phi * phi - 0.0002389523716779765 * phi + 0.06428970810739926 - 0.01;
    }
  }

  //==========//  PARTICLE = PI- PION (END)  //==========//

  //==========//  PARTICLE = PROTON  //==========//

  if (ivec == 3) {
    // The following lines should be added up in the order given for the full correction
    // Applying this code as given will give the exact corrections of this analysis
    // These parameters will be combined into a single line at a later point

    if (sec == 1) {
      dp = (5.415e-04) * pp * pp + (-1.0262e-02) * pp + (7.78075e-03);
      dp = dp + ((1.2129e-04) * pp * pp + (1.5373e-04) * pp + (-2.7084e-04));
    }
    if (sec == 2) {
      dp = (-9.5439e-04) * pp * pp + (-2.86273e-03) * pp + (3.38149e-03);
      dp = dp + ((-1.6890e-03) * pp * pp + (4.3744e-03) * pp + (-2.1218e-03));
    }
    if (sec == 3) {
      dp = (-5.5541e-04) * pp * pp + (-7.69739e-03) * pp + (5.7692e-03);
      dp = dp + ((7.6422e-04) * pp * pp + (-1.5425e-03) * pp + (5.4255e-04));
    }
    if (sec == 4) {
      dp = (-1.944e-04) * pp * pp + (-5.77104e-03) * pp + (3.42399e-03);
      dp = dp + ((1.1174e-03) * pp * pp + (-3.2747e-03) * pp + (2.3687e-03));
    }
    if (sec == 5) {
      dp = (1.54009e-03) * pp * pp + (-1.69437e-02) * pp + (1.04656e-02);
      dp = dp + ((-2.1067e-04) * pp * pp + (1.2266e-03) * pp + (-1.0553e-03));
    }
    if (sec == 6) {
      dp = (2.38182e-03) * pp * pp + (-2.07301e-02) * pp + (1.72325e-02);
      dp = dp + ((-3.6002e-04) * pp * pp + (8.9582e-04) * pp + (-1.0093e-03));
    }
  }

  //==========//  PARTICLE = PROTON (END)  //==========//

  return dp / pp;
}


///////////////////////////////////////////// new Momentum Corrections Last Updated: December 23, 2022 (Not Finalized)

// ////// These corrections have not been fully updated yet as the π - pion corrections
// ////// are still under development as of 12 -23 -2022
// /////////////////

// // auto dppC = (float Px, float Py, float Pz, int sec, int ivec) {
// double mom_corr::dppC(float Px, float Py, float Pz, int sec, int ivec) {
//   // ivec = 0 --> Electron Corrections
//   // ivec = 1 --> π+ Corrections
//   // ivec = 2 --> π- Corrections
//   // ivec = 3 --> Proton Corrections (NOT UPDATED YET)

//   // Momentum Magnitude
//   double pp = sqrt(Px * Px + Py * Py + Pz * Pz);

//   // Initializing the correction factor
//   double dp = 0;

//   // Defining Phi Angle
//   double Phi = (180 / 3.1415926) * atan2(Py, Px);

//   // (Initial) Shift of the Phi Angle (done to realign sectors whose data is separated when plotted from ±180˚)
//   if (((sec == 4 || sec == 3) && Phi < 0) || (sec > 4 && Phi < 90)) {
//     Phi += 360;
//   }

//   // Getting Local Phi Angle
//   double PhiLocal = Phi - (sec - 1) * 60;

//   // Applying Shift Functions to Phi Angles (local shifted phi = phi)
//   double phi = PhiLocal;

//   // For Electron Shift
//   if (ivec == 0) {
//     phi = PhiLocal - 30 / pp;
//   }

//   // For π+ Pion/Proton Shift
//   if (ivec == 1 || ivec == 3) {
//     phi = PhiLocal + (32 / (pp - 0.05));
//   }

//   // For π- Pion Shift
//   if (ivec == 2) {
//     phi = PhiLocal - (32 / (pp - 0.05));
//   }

//   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   //==================================================================================================================================//
//   //=======================//=======================//     Electron Corrections
//   ////=======================//=======================//
//   //==================================================================================================================================//
//   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//   if (ivec == 0) {
//     if (sec == 1) {
//       // The CONTINUOUS QUADRATIC function predicted for ∆p_{El} for [Cor = Uncorrected][Sector 1] is:
//       dp = ((-4.3303e-06) * phi * phi + (1.1006e-04) * phi + (-5.7235e-04)) * pp * pp +
//            ((3.2555e-05) * phi * phi + (-0.0014559) * phi + (0.0014878)) * pp +
//            ((-1.9577e-05) * phi * phi + (0.0017996) * phi + (0.025963));
//     }
//     if (sec == 2) {
//       // The CONTINUOUS QUADRATIC function predicted for ∆p_{El} for [Cor = Uncorrected][Sector 2] is:
//       dp = ((-9.8045e-07) * phi * phi + (6.7395e-05) * phi + (-4.6757e-05)) * pp * pp +
//            ((-1.4958e-05) * phi * phi + (-0.0011191) * phi + (-0.0025143)) * pp +
//            ((1.2699e-04) * phi * phi + (0.0033121) * phi + (0.020819));
//     }
//     if (sec == 3) {
//       // The CONTINUOUS QUADRATIC function predicted for ∆p_{El} for [Cor = Uncorrected][Sector 3] is:
//       dp = ((-5.9459e-07) * phi * phi + (-2.8289e-05) * phi + (-4.3541e-04)) * pp * pp +
//            ((-1.5025e-05) * phi * phi + (5.7730e-04) * phi + (-0.0077582)) * pp +
//            ((7.3348e-05) * phi * phi + (-0.001102) * phi + (0.057052));
//     }
//     if (sec == 4) {
//       // The CONTINUOUS QUADRATIC function predicted for ∆p_{El} for [Cor = Uncorrected][Sector 4] is:
//       dp = ((-2.2714e-06) * phi * phi + (-3.0360e-05) * phi + (-8.9322e-04)) * pp * pp +
//            ((2.9737e-05) * phi * phi + (5.1142e-04) * phi + (0.0045641)) * pp +
//            ((-1.0582e-04) * phi * phi + (-5.6852e-04) * phi + (0.027506));
//     }
//     if (sec == 5) {
//       // The CONTINUOUS QUADRATIC function predicted for ∆p_{El} for [Cor = Uncorrected][Sector 5] is:
//       dp = ((-1.1490e-06) * phi * phi + (-6.2147e-06) * phi + (-4.7235e-04)) * pp * pp +
//            ((3.7039e-06) * phi * phi + (-1.5943e-04) * phi + (-8.5238e-04)) * pp +
//            ((4.4069e-05) * phi * phi + (0.0014152) * phi + (0.031933));
//     }
//     if (sec == 6) {
//       // The CONTINUOUS QUADRATIC function predicted for ∆p_{El} for [Cor = Uncorrected][Sector 6] is:
//       dp = ((1.1076e-06) * phi * phi + (4.0156e-05) * phi + (-1.6341e-04)) * pp * pp +
//            ((-2.8613e-05) * phi * phi + (-5.1861e-04) * phi + (-0.0056437)) * pp +
//            ((1.2419e-04) * phi * phi + (4.9084e-04) * phi + (0.049976));
//     }
//   }

//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   //====================================================================================================================================//
//   //======================//======================//     Electron Corrections (End)
//   ////======================//======================//
//   //====================================================================================================================================//
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   //====================================================================================================================================//
//   //=========================//=========================//     π+ Corrections
//   ////=========================//=========================//
//   //====================================================================================================================================//
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   if (ivec == 1) {
//     if (sec == 1) {
//       dp = ((-5.4904e-07) * phi * phi + (-1.4436e-05) * phi + (3.1534e-04)) * pp * pp +
//            ((3.8231e-06) * phi * phi + (3.6582e-04) * phi + (-0.0046759)) * pp +
//            ((-5.4913e-06) * phi * phi + (-4.0157e-04) * phi + (0.010767));
//       dp = dp + ((6.1103e-07) * phi * phi + (5.5291e-06) * phi + (-1.9120e-04)) * pp * pp +
//            ((-3.2300e-06) * phi * phi + (1.5377e-05) * phi + (7.5279e-04)) * pp +
//            ((2.1434e-06) * phi * phi + (-6.9572e-06) * phi + (-7.9333e-05));
//       dp = dp + ((-1.3049e-06) * phi * phi + (1.1295e-05) * phi + (4.5797e-04)) * pp * pp +
//            ((9.3122e-06) * phi * phi + (-5.1074e-05) * phi + (-0.0030757)) * pp +
//            ((-1.3102e-05) * phi * phi + (2.2153e-05) * phi + (0.0040938));
//     }
//     if (sec == 2) {
//       dp = ((-1.0087e-06) * phi * phi + (2.1319e-05) * phi + (7.8641e-04)) * pp * pp +
//            ((6.7485e-06) * phi * phi + (7.3716e-05) * phi + (-0.0094591)) * pp +
//            ((-1.1820e-05) * phi * phi + (-3.8103e-04) * phi + (0.018936));
//       dp = dp + ((8.8155e-07) * phi * phi + (-2.8257e-06) * phi + (-2.6729e-04)) * pp * pp +
//            ((-5.4499e-06) * phi * phi + (3.8397e-05) * phi + (0.0015914)) * pp +
//            ((6.8926e-06) * phi * phi + (-5.9386e-05) * phi + (-0.0021749));
//       dp = dp + ((-2.0147e-07) * phi * phi + (1.1061e-05) * phi + (3.8827e-04)) * pp * pp +
//            ((4.9294e-07) * phi * phi + (-6.0257e-05) * phi + (-0.0022087)) * pp +
//            ((9.8548e-07) * phi * phi + (5.9047e-05) * phi + (0.0022905));
//     }
//     if (sec == 3) {
//       dp = ((8.6722e-08) * phi * phi + (-1.7975e-05) * phi + (4.8118e-05)) * pp * pp +
//            ((2.6273e-06) * phi * phi + (3.1453e-05) * phi + (-0.0015943)) * pp +
//            ((-6.4463e-06) * phi * phi + (-5.8990e-05) * phi + (0.0041703));
//       dp = dp + ((9.6317e-07) * phi * phi + (-1.7659e-06) * phi + (-8.8318e-05)) * pp * pp +
//            ((-5.1346e-06) * phi * phi + (8.3318e-06) * phi + (3.7723e-04)) * pp +
//            ((3.9548e-06) * phi * phi + (-6.9614e-05) * phi + (2.1393e-04));
//       dp = dp + ((5.6438e-07) * phi * phi + (8.1678e-06) * phi + (-9.4406e-05)) * pp * pp +
//            ((-3.9074e-06) * phi * phi + (-6.5174e-05) * phi + (5.4218e-04)) * pp +
//            ((6.3198e-06) * phi * phi + (1.0611e-04) * phi + (-4.5749e-04));
//     }
//     if (sec == 4) {
//       dp = ((4.3406e-07) * phi * phi + (-4.9036e-06) * phi + (2.3064e-04)) * pp * pp +
//            ((1.3624e-06) * phi * phi + (3.2907e-05) * phi + (-0.0034872)) * pp +
//            ((-5.1017e-06) * phi * phi + (2.4593e-05) * phi + (0.0092479));
//       dp = dp + ((6.0218e-07) * phi * phi + (-1.4383e-05) * phi + (-3.1999e-05)) * pp * pp +
//            ((-1.1243e-06) * phi * phi + (9.3884e-05) * phi + (-4.1985e-04)) * pp +
//            ((-1.8808e-06) * phi * phi + (-1.2222e-04) * phi + (0.0014037));
//       dp = dp + ((-2.5490e-07) * phi * phi + (-8.5120e-07) * phi + (7.9109e-05)) * pp * pp +
//            ((2.5879e-06) * phi * phi + (8.6108e-06) * phi + (-5.1533e-04)) * pp +
//            ((-4.4521e-06) * phi * phi + (-1.7012e-05) * phi + (7.4848e-04));
//     }
//     if (sec == 5) {
//       dp = ((2.4292e-07) * phi * phi + (8.8741e-06) * phi + (2.9482e-04)) * pp * pp +
//            ((3.7229e-06) * phi * phi + (7.3215e-06) * phi + (-0.0050685)) * pp +
//            ((-1.1974e-05) * phi * phi + (-1.3043e-04) * phi + (0.0078836));
//       dp = dp + ((1.0867e-06) * phi * phi + (-7.7630e-07) * phi + (-4.4930e-05)) * pp * pp +
//            ((-5.6564e-06) * phi * phi + (-1.3417e-05) * phi + (2.5224e-04)) * pp +
//            ((6.8460e-06) * phi * phi + (9.0495e-05) * phi + (-4.6587e-04));
//       dp = dp + ((8.5720e-07) * phi * phi + (-6.7464e-06) * phi + (-4.0944e-05)) * pp * pp +
//            ((-4.7370e-06) * phi * phi + (5.8808e-05) * phi + (1.9047e-04)) * pp +
//            ((5.7404e-06) * phi * phi + (-1.1105e-04) * phi + (-1.9392e-04));
//     }
//     if (sec == 6) {
//       dp = ((2.1191e-06) * phi * phi + (-3.3710e-05) * phi + (2.5741e-04)) * pp * pp +
//            ((-1.2915e-05) * phi * phi + (2.3753e-04) * phi + (-2.6882e-04)) * pp +
//            ((2.2676e-05) * phi * phi + (-2.3115e-04) * phi + (-0.001283));
//       dp = dp + ((6.0270e-07) * phi * phi + (-6.8200e-06) * phi + (1.3103e-04)) * pp * pp +
//            ((-1.8745e-06) * phi * phi + (3.8646e-05) * phi + (-8.8056e-04)) * pp +
//            ((2.0885e-06) * phi * phi + (-3.4932e-05) * phi + (4.5895e-04));
//       dp = dp + ((4.7349e-08) * phi * phi + (-5.7528e-06) * phi + (-3.4097e-06)) * pp * pp +
//            ((1.7731e-06) * phi * phi + (3.5865e-05) * phi + (-5.7881e-04)) * pp +
//            ((-9.7008e-06) * phi * phi + (-4.1836e-05) * phi + (0.0035403));
//     }
//   }

//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   //====================================================================================================================================//
//   //=========================//=========================//  π+ Corrections (End)
//   ////=========================//=========================//
//   //====================================================================================================================================//
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   //====================================================================================================================================//
//   //==================//==================//    π- Corrections (Updated as of 01-13-2023)
//   ////==================//==================//
//   //====================================================================================================================================//
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//   if (ivec == 2) {
//     if (sec == 1) {
//       dp = (-1.6287E-05 * phi * phi + -4.9698E-04 * phi + -1.5461E-03) * pp * pp +
//            (6.1162E-05 * phi * phi + 1.5151E-03 * phi + -6.2379E-03) * pp + -5.2911E-05 * phi * phi +
//            -1.1023E-03 * phi + 3.3557E-02 - 0.008;
//     }
//     if (sec == 2) {
//       dp = (1.4819E-06 * phi * phi + -3.7222E-05 * phi + 1.3426E-03) * pp * pp +
//            (-1.4990E-06 * phi * phi + 9.6467E-05 * phi + -1.9363E-02) * pp + 1.1426E-06 * phi * phi + 4.5750E-05 * phi +
//            3.7193E-02 - 0.005;
//     }

//     if (sec == 3) {
//       dp = (-1.2521E-05 * phi * phi + -4.0605E-05 * phi + 7.1584E-04) * pp * pp +
//            (5.5105E-05 * phi * phi + 3.5087E-04 * phi + -1.9071E-02) * pp + -4.9659E-05 * phi * phi +
//            -2.9078E-04 * phi + 3.6345E-02 - 0.005;
//     }

//     if (sec == 4) {
//       dp = (-4.6480E-07 * phi * phi + -1.7373E-06 * phi + 3.4723E-03) * pp * pp +
//            (6.8267E-07 * phi * phi + 1.3368E-04 * phi + -2.4534E-02) * pp + 9.9275E-06 * phi * phi + -1.5813E-04 * phi +
//            4.0981E-02 - 0.004;
//     }

//     if (sec == 5) {
//       dp = (-2.7381E-06 * phi * phi + 6.2363E-05 * phi + 1.3284E-03 - 0.004) * pp * pp +
//            (8.3662E-06 * phi * phi + -2.0197E-04 * phi + -1.5436E-02 + 0.01) * pp + -1.3453E-05 * phi * phi +
//            9.2336E-07 * phi + 3.8825E-02 - 0.008;
//     }

//     if (sec == 6) {
//       dp = (4.7125E-06 * phi * phi + 9.8559E-05 * phi + -1.8335E-03) * pp * pp +
//            (-1.8723E-05 * phi * phi + -5.5862E-04 * phi + -5.5286E-03) * pp + 1.7244E-05 * phi * phi +
//            4.8344E-04 * phi + 3.0237E-02 - 0.008;
//     }
//   }

//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   //====================================================================================================================================//
//   //=======================//=======================//      π- Corrections (End)
//   ////=======================//=======================//
//   //====================================================================================================================================//
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   //====================================================================================================================================//
//   //=======================//=======================//     All Proton Corrections
//   ////=======================//=======================//
//   //====================================================================================================================================//
//   ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   if (ivec == 3) {
//     if (sec == 1) {
//       dp = ((1 + TMath::Sign(1, (pp - 1.4))) / 2) * ((4.4034e-03) * pp + (-0.01703)) +
//            ((1 + TMath::Sign(1, -(pp - 1.4))) / 2) *
//                ((-0.10898) * (pp - 1.4) * (pp - 1.4) + (-0.09574) * (pp - 1.4) + ((4.4034e-03) * 1.4 + (-0.01703)));
//     }
//     if (sec == 2) {
//       dp = ((1 + TMath::Sign(1, (pp - 1.5))) / 2) * ((0.01318) * pp + (-0.03403)) +
//            ((1 + TMath::Sign(1, -(pp - 1.5))) / 2) *
//                ((-0.09829) * (pp - 1.5) * (pp - 1.5) + (-0.0986) * (pp - 1.5) + ((0.01318) * 1.5 + (-0.03403)));
//     }
//     if (sec == 3) {
//       dp =
//           ((1 + TMath::Sign(1, (pp - 1.05))) / 2) * ((-4.7052e-03) * pp + (1.2410e-03)) +
//           ((1 + TMath::Sign(1, -(pp - 1.05))) / 2) * ((-0.22721) * (pp - 1.05) * (pp - 1.05) +
//                                                       (-0.09702) * (pp - 1.05) + ((-4.7052e-03) * 1.05 + (1.2410e-03)));
//     }
//     if (sec == 4) {
//       dp = ((1 + TMath::Sign(1, (pp - 1.4))) / 2) * ((-1.0900e-03) * pp + (-4.0573e-03)) +
//            ((1 + TMath::Sign(1, -(pp - 1.4))) / 2) *
//                ((-0.09236) * (pp - 1.4) * (pp - 1.4) + (-0.073) * (pp - 1.4) + ((-1.0900e-03) * 1.4 + (-4.0573e-03)));
//     }
//     if (sec == 5) {
//       dp = ((1 + TMath::Sign(1, (pp - 1.5))) / 2) * ((7.3965e-03) * pp + (-0.02428)) +
//            ((1 + TMath::Sign(1, -(pp - 1.5))) / 2) *
//                ((-0.09539) * (pp - 1.5) * (pp - 1.5) + (-0.09263) * (pp - 1.5) + ((7.3965e-03) * 1.5 + (-0.02428)));
//     }
//     if (sec == 6) {
//       dp =
//           ((1 + TMath::Sign(1, (pp - 1.15))) / 2) * ((-7.6214e-03) * pp + (8.1014e-03)) +
//           ((1 + TMath::Sign(1, -(pp - 1.15))) / 2) * ((-0.12718) * (pp - 1.15) * (pp - 1.15) +
//                                                       (-0.06626) * (pp - 1.15) + ((-7.6214e-03) * 1.15 + (8.1014e-03)));
//     }
//   }

//   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   //=====================================================================================================================================//
//   //=======================//=======================//    End of Proton Corrections
//   ////=======================//=======================//
//   //=====================================================================================================================================//
//   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//   return dp / pp;
// };

// // // The following code is for the Energy Loss Corrections for the proton
// // double dE_loss = 0;
// // // Inbending Energy Loss Correction //
// // if (proth < 27) {
// //   dE_loss = exp(-2.739 - 3.932 * pro) + 0.002907;
// // }
// // if (proth > 27) {
// //   dE_loss = exp(-1.2 - 4.228 * pro) + 0.007502;
// // }
// // double feloss = (pro + dE_loss) / pro;

// // // Below shows how the corrections are to be applied using the ROOT momentum 4-vector using the above code:
// // auto fe = dppC(ex, ey, ez, esec, 0) + 1;
// // auto fpip = dppC(pipx, pipy, pipz, pipsec, 1) + 1;
// // auto fpim = dppC(pimx, pimy, pimz, pimsec, 2) + 1;
// // auto fpro = dppC(prox * feloss, proy* feloss, proz* feloss, prosec, 3) + 1;

// // auto eleC = ROOT::Math::PxPyPzMVector(ex * fe, ey* fe, ez* fe, 0);
// // auto pipC = ROOT::Math::PxPyPzMVector(pipx * fpip, pipy* fpip, pipz* fpip, 0.13957);
// // auto pimC = ROOT::Math::PxPyPzMVector(pimx * fpim, pimy* fpim, pimz* fpim, 0.13957);
// // auto proC = ROOT::Math::PxPyPzMVector(prox * feloss * fpro, proy* feloss* fpro, proz* feloss* fpro, 0.938);

// //// NEW method FOR ALL W RANGE USING TWOPION skim data
// double CDProt[3][5] = {{-0.1469, 0.7793, -1.339, 0.8013, -0.1791},
//                        {0.03082, -0.01495, -0.1974, 0.2642, -0.0895},
//                        {-0.02351, 0.2006, -0.4646, 0.43, -0.1278}};

// float mom_corr::CD_prot_Hmom_corr(float mom_, float phi_, float alpha_prot_mom_corr_CD[3]) {
//   if (phi_ > 270 || phi_ <= 30) {
//     return mom_ - alpha_prot_mom_corr_CD[0] * (CDProt[0][0] * pow(mom_, 4) + CDProt[0][1] * pow(mom_, 3) +
//                                                CDProt[0][2] * pow(mom_, 2) + CDProt[0][3] * mom_ + CDProt[0][4]);
//   } else if (phi_ > 30 && phi_ <= 150) {
//     return mom_ - alpha_prot_mom_corr_CD[1] * (CDProt[1][0] * pow(mom_, 4) + CDProt[1][1] * pow(mom_, 3) +
//                                                CDProt[1][2] * pow(mom_, 2) + CDProt[1][3] * mom_ + CDProt[1][4]);
//   } else if (phi_ > 150 && phi_ <= 270) {
//     return mom_ - alpha_prot_mom_corr_CD[2] * (CDProt[2][0] * pow(mom_, 4) + CDProt[2][1] * pow(mom_, 3) +
//                                                CDProt[2][2] * pow(mom_, 2) + CDProt[2][3] * mom_ + CDProt[2][4]);
//   } else
//     return NAN;
// }

// float FDProt[6][2] = {{-0.005096, 0.002893}, {-0.005005, 0.010056}, {-0.007, 0.01459},
//                       {-0.006756, 0.01022},  {-0.004307, 0.003967}, {-0.010124, 0.006184}};

// float mom_corr::FD_prot_Hmom_corr(float mom_, float dc_sec, float alpha_prot) {
//   if (dc_sec == 1) {
//     return mom_ - alpha_prot * (FDProt[0][0] * pow(mom_, 1) + FDProt[0][1]);
//   } else if (dc_sec == 2) {
//     return mom_ - alpha_prot * (FDProt[1][0] * pow(mom_, 1) + FDProt[1][1]);
//   } else if (dc_sec == 3) {
//     return mom_ - alpha_prot * (FDProt[2][0] * pow(mom_, 1) + FDProt[2][1]);
//   } else if (dc_sec == 4) {
//     return mom_ - alpha_prot * (FDProt[3][0] * pow(mom_, 1) + FDProt[3][1]);
//   } else if (dc_sec == 5) {
//     return mom_ - alpha_prot * (FDProt[4][0] * pow(mom_, 1) + FDProt[4][1]);
//   } else if (dc_sec == 6) {
//     return mom_ - alpha_prot * (FDProt[5][0] * pow(mom_, 1) + FDProt[5][1]);
//   } else
//     return NAN;
// }

// ///// pip
// double CDPip[3][5] = {{0.02008, -0.01386, -0.02003, -0.0848, 0.01132},
//                       {0.003, 0.03384, -0.09216, 0.0508, -0.004814},
//                       {-0.0642, 0.2113, -0.1926, 0.0702, -0.01304}};

// float mom_corr::CD_pip_Hmom_corr(float mom_, float phi_, float alpha_pip_mom_corr_CD[3]) {
//   if (phi_ > 270 || phi_ <= 30) {
//     return mom_ - alpha_pip_mom_corr_CD[0] * (CDPip[0][0] * pow(mom_, 4) + CDPip[0][1] * pow(mom_, 3) +
//                                               CDPip[0][2] * pow(mom_, 2) + CDPip[0][3] * mom_ + CDPip[0][4]);
//   } else if (phi_ > 30 && phi_ <= 150) {
//     return mom_ - alpha_pip_mom_corr_CD[1] * (CDPip[1][0] * pow(mom_, 4) + CDPip[1][1] * pow(mom_, 3) +
//                                               CDPip[1][2] * pow(mom_, 2) + CDPip[1][3] * mom_ + CDPip[1][4]);
//   } else if (phi_ > 150 && phi_ <= 270) {
//     return mom_ - alpha_pip_mom_corr_CD[2] * (CDPip[2][0] * pow(mom_, 4) + CDPip[2][1] * pow(mom_, 3) +
//                                               CDPip[2][2] * pow(mom_, 2) + CDPip[2][3] * mom_ + CDPip[2][4]);
//   } else
//     return NAN;
// }

// float FDPip[6][3] = {{-0.000578, 0.004654, -0.01345},  {-0.0001847, 0.002308, -0.007244},
//                      {-0.000411, 0.001353, -0.000333}, {0.000643, -0.006462, 0.007042},
//                      {0.002375, -0.01843, 0.02043},    {0.001001, -0.00788, 0.004395}};

// float mom_corr::FD_pip_Hmom_corr(float mom_, float dc_sec, float alpha_pip) {
//   if (dc_sec == 1) {
//     return mom_ - alpha_pip * (FDPip[0][0] * pow(mom_, 2) + FDPip[0][1] * mom_ + FDPip[0][2]);
//   } else if (dc_sec == 2) {
//     return mom_ - alpha_pip * (FDPip[1][0] * pow(mom_, 2) + FDPip[1][1] * mom_ + FDPip[1][2]);
//   } else if (dc_sec == 3) {
//     return mom_ - alpha_pip * (FDPip[2][0] * pow(mom_, 2) + FDPip[2][1] * mom_ + FDPip[2][2]);
//   } else if (dc_sec == 4) {
//     return mom_ - alpha_pip * (FDPip[3][0] * pow(mom_, 2) + FDPip[3][1] * mom_ + FDPip[3][2]);
//   } else if (dc_sec == 5) {
//     return mom_ - alpha_pip * (FDPip[4][0] * pow(mom_, 2) + FDPip[4][1] * mom_ + FDPip[4][2]);
//   } else if (dc_sec == 6) {
//     return mom_ - alpha_pip * (FDPip[5][0] * pow(mom_, 2) + FDPip[5][1] * mom_ + FDPip[5][2]);
//   } else
//     return NAN;
// }
// ///// pim
// double CDPim[3][5] = {{-0.09753, 0.394, -0.4932, 0.2395, -0.03833},
//                       {-0.03668, 0.1414, -0.1572, 0.09827, -0.011246},
//                       {-0.1105, 0.4443, -0.566, 0.2224, -0.02582}};

// float mom_corr::CD_pim_Hmom_corr(float mom_, float phi_, float alpha_pim_mom_corr_CD[3]) {
//   if (phi_ > 270 || phi_ <= 30) {
//     return mom_ - alpha_pim_mom_corr_CD[0] * (CDPim[0][0] * pow(mom_, 4) + CDPim[0][1] * pow(mom_, 3) +
//                                               CDPim[0][2] * pow(mom_, 2) + CDPim[0][3] * mom_ + CDPim[0][4]);
//   } else if (phi_ > 30 && phi_ <= 150) {
//     return mom_ - alpha_pim_mom_corr_CD[1] * (CDPim[1][0] * pow(mom_, 4) + CDPim[1][1] * pow(mom_, 3) +
//                                               CDPim[1][2] * pow(mom_, 2) + CDPim[1][3] * mom_ + CDPim[1][4]);
//   } else if (phi_ > 150 && phi_ <= 270) {
//     return mom_ - alpha_pim_mom_corr_CD[2] * (CDPim[2][0] * pow(mom_, 4) + CDPim[2][1] * pow(mom_, 3) +
//                                               CDPim[2][2] * pow(mom_, 2) + CDPim[2][3] * mom_ + CDPim[2][4]);
//   } else
//     return NAN;
// }

// float FDPim[6][3] = {{0.002954, -0.01746, 0.01444},     {-0.0002172, -0.0009136, 0.0006347},
//                      {-0.0009513, 0.004047, -0.003483}, {0.00122, -0.006836, 0.002924},
//                      {0.001049, -0.002333, -0.00816},   {0.0010195, -0.00813, 0.00096}};

// float mom_corr::FD_pim_Hmom_corr(float mom_, float dc_sec, float alpha_pim) {
//   if (dc_sec == 1) {
//     return mom_ - alpha_pim * (FDPim[0][0] * pow(mom_, 2) + FDPim[0][1] * mom_ + FDPim[0][2]);
//   } else if (dc_sec == 2) {
//     return mom_ - alpha_pim * (FDPim[1][0] * pow(mom_, 2) + FDPim[1][1] * mom_ + FDPim[1][2]);
//   } else if (dc_sec == 3) {
//     return mom_ - alpha_pim * (FDPim[2][0] * pow(mom_, 2) + FDPim[2][1] * mom_ + FDPim[2][2]);
//   } else if (dc_sec == 4) {
//     return mom_ - alpha_pim * (FDPim[3][0] * pow(mom_, 2) + FDPim[3][1] * mom_ + FDPim[3][2]);
//   } else if (dc_sec == 5) {
//     return mom_ - alpha_pim * (FDPim[4][0] * pow(mom_, 2) + FDPim[4][1] * mom_ + FDPim[4][2]);
//   } else if (dc_sec == 6) {
//     return mom_ - alpha_pim * (FDPim[5][0] * pow(mom_, 2) + FDPim[5][1] * mom_ + FDPim[5][2]);
//   } else
//     return NAN;
// }

// our hadron momentum correction come from here:

// ////////////// ###########  These are the corrections we used for our w< 2.55 , FD seperated in theta angle cases, as
// // June 04 2023

// double CDProt[3][5] = {{-0.2578, 1.334, -2.3, 1.489, -0.3545},
//                        {-0.0736, 0.4873, -1.048, 0.862, -0.2374},
//                        {-0.0928, 0.5454, -1.074, 0.874, -0.2386}};
// // corrections

// float mom_corr::CD_prot_Hmom_corr(float mom_, float phi_, float alpha_prot_mom_corr_CD[3]) {
//   if (phi_ > 270 || phi_ <= 30) {
//     return mom_ - alpha_prot_mom_corr_CD[0] * (CDProt[0][0] * pow(mom_, 4) + CDProt[0][1] * pow(mom_, 3) +
//                                                CDProt[0][2] * pow(mom_, 2) + CDProt[0][3] * mom_ + CDProt[0][4]);
//   } else if (phi_ > 30 && phi_ <= 150) {
//     return mom_ - alpha_prot_mom_corr_CD[1] * (CDProt[1][0] * pow(mom_, 4) + CDProt[1][1] * pow(mom_, 3) +
//                                                CDProt[1][2] * pow(mom_, 2) + CDProt[1][3] * mom_ + CDProt[1][4]);
//   } else if (phi_ > 150 && phi_ <= 270) {
//     return mom_ - alpha_prot_mom_corr_CD[2] * (CDProt[2][0] * pow(mom_, 4) + CDProt[2][1] * pow(mom_, 3) +
//                                                CDProt[2][2] * pow(mom_, 2) + CDProt[2][3] * mom_ + CDProt[2][4]);
//   } else
//     return NAN;
// }
// float FDProtL[2][6][4] = {{{0.0004573, 0.000176, -0.01131, 0.011406},
//                            {0.001105, -0.0095, 0.01335, 0.003487},
//                            {-0.001555, 0.007084, -0.0193, 0.02263},
//                            {0.000816, -0.004272, -0.00223, 0.01088},
//                            {-0.0002866, 0.004208, -0.02225, 0.0223},
//                            {0.00344, -0.02016, 0.01811, 0.004692}},
//                           {{-0.002356, 0.01585, -0.03143, 0.001087},
//                            {0.003145, -0.01888, 0.0341, -0.02017},
//                            {-0.0005207, 0.002861, -0.004078, 0.001318},
//                            {-0.002945, 0.0183, -0.03568, 0.01697},
//                            {-0.003153, 0.02109, -0.04117, 0.01209},
//                            {-0.00556, 0.0387, -0.08203, 0.02975}}};

// float mom_corr::FD_prot_Hmom_corr_lower_All_FD(float mom_, float dc_sec, float alpha_prot_mom_corr_FD) {
//   if (dc_sec == 1) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtL[0][0][0] * pow(mom_, 3) + FDProtL[0][0][1] * pow(mom_, 2) +
//                                             FDProtL[0][0][2] * mom_ + FDProtL[0][0][3]);
//   } else if (dc_sec == 2) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtL[0][1][0] * pow(mom_, 3) + FDProtL[0][1][1] * pow(mom_, 2) +
//                                             FDProtL[0][1][2] * mom_ + FDProtL[0][1][3]);
//   } else if (dc_sec == 3) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtL[0][2][0] * pow(mom_, 3) + FDProtL[0][2][1] * pow(mom_, 2) +
//                                             FDProtL[0][2][2] * mom_ + FDProtL[0][2][3]);
//   } else if (dc_sec == 4) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtL[0][3][0] * pow(mom_, 3) + FDProtL[0][3][1] * pow(mom_, 2) +
//                                             FDProtL[0][3][2] * mom_ + FDProtL[0][3][3]);
//   } else if (dc_sec == 5) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtL[0][4][0] * pow(mom_, 3) + FDProtL[0][4][1] * pow(mom_, 2) +
//                                             FDProtL[0][4][2] * mom_ + FDProtL[0][4][3]);
//   } else if (dc_sec == 6) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtL[0][5][0] * pow(mom_, 3) + FDProtL[0][5][1] * pow(mom_, 2) +
//                                             FDProtL[0][5][2] * mom_ + FDProtL[0][5][3]);
//   } else
//     return NAN;
// }

// float mom_corr::FD_prot_Hmom_corr_lower_Except_All_FD(float mom_, float dc_sec, float alpha_prot_mom_corr_FD) {
//   if (dc_sec == 1) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtL[1][0][0] * pow(mom_, 3) + FDProtL[1][0][1] * pow(mom_, 2) +
//                                             FDProtL[1][0][2] * mom_ + FDProtL[1][0][3]);
//   } else if (dc_sec == 2) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtL[1][1][0] * pow(mom_, 3) + FDProtL[1][1][1] * pow(mom_, 2) +
//                                             FDProtL[1][1][2] * mom_ + FDProtL[1][1][3]);
//   } else if (dc_sec == 3) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtL[1][2][0] * pow(mom_, 3) + FDProtL[1][2][1] * pow(mom_, 2) +
//                                             FDProtL[1][2][2] * mom_ + FDProtL[1][2][3]);
//   } else if (dc_sec == 4) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtL[1][3][0] * pow(mom_, 3) + FDProtL[1][3][1] * pow(mom_, 2) +
//                                             FDProtL[1][3][2] * mom_ + FDProtL[1][3][3]);
//   } else if (dc_sec == 5) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtL[1][4][0] * pow(mom_, 3) + FDProtL[1][4][1] * pow(mom_, 2) +
//                                             FDProtL[1][4][2] * mom_ + FDProtL[1][4][3]);
//   } else if (dc_sec == 6) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtL[1][5][0] * pow(mom_, 3) + FDProtL[1][5][1] * pow(mom_, 2) +
//                                             FDProtL[1][5][2] * mom_ + FDProtL[1][5][3]);
//   } else
//     return NAN;
// }

// double FDProtH[2][6][4] = {{{-0.01047, 0.0603, -0.08997, 0.0436},
//                             {-0.000767, 0.008934, -0.01267, 0.011955},
//                             {-0.013275, 0.0709, -0.1035, 0.0569},
//                             {0.003572, -0.007786, 0.003706, 0.01051},
//                             {-0.00461, 0.03168, -0.05597, 0.0392},
//                             {-0.00398, 0.02565, -0.04474, 0.03076}},
//                            {{-0.006363, 0.0391, -0.06036, 0.0166},
//                             {0.00442, -0.02858, 0.06744, -0.04587},
//                             {0.000977, -0.00703, 0.02438, -0.00925},
//                             {0.003227, -0.01921, 0.04, -0.01855},
//                             {-0.01031, 0.0607, -0.08954, 0.03084},
//                             {-0.006283, 0.03836, -0.05853, 0.00891}}};

// float mom_corr::FD_prot_Hmom_corr_upper_All_FD(float mom_, float dc_sec, float alpha_prot_mom_corr_FD) {
//   if (dc_sec == 1) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtH[0][0][0] * pow(mom_, 3) + FDProtH[0][0][1] * pow(mom_, 2) +
//                                             FDProtH[0][0][2] * mom_ + FDProtH[0][0][3]);
//   } else if (dc_sec == 2) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtH[0][1][0] * pow(mom_, 3) + FDProtH[0][1][1] * pow(mom_, 2) +
//                                             FDProtH[0][1][2] * mom_ + FDProtH[0][1][3]);
//   } else if (dc_sec == 3) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtH[0][2][0] * pow(mom_, 3) + FDProtH[0][2][1] * pow(mom_, 2) +
//                                             FDProtH[0][2][2] * mom_ + FDProtH[0][2][3]);
//   } else if (dc_sec == 4) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtH[0][3][0] * pow(mom_, 3) + FDProtH[0][3][1] * pow(mom_, 2) +
//                                             FDProtH[0][3][2] * mom_ + FDProtH[0][3][3]);
//   } else if (dc_sec == 5) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtH[0][4][0] * pow(mom_, 3) + FDProtH[0][4][1] * pow(mom_, 2) +
//                                             FDProtH[0][4][2] * mom_ + FDProtH[0][4][3]);
//   } else if (dc_sec == 6) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtH[0][5][0] * pow(mom_, 3) + FDProtH[0][5][1] * pow(mom_, 2) +
//                                             FDProtH[0][5][2] * mom_ + FDProtH[0][5][3]);
//   } else
//     return NAN;
// }

// float mom_corr::FD_prot_Hmom_corr_upper_Except_All_FD(float mom_, float dc_sec, float alpha_prot_mom_corr_FD) {
//   if (dc_sec == 1) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtH[1][0][0] * pow(mom_, 3) + FDProtH[1][0][1] * pow(mom_, 2) +
//                                             FDProtH[1][0][2] * mom_ + FDProtH[1][0][3]);
//   } else if (dc_sec == 2) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtH[1][1][0] * pow(mom_, 3) + FDProtH[1][1][1] * pow(mom_, 2) +
//                                             FDProtH[1][1][2] * mom_ + FDProtH[1][1][3]);
//   } else if (dc_sec == 3) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtH[1][2][0] * pow(mom_, 3) + FDProtH[1][2][1] * pow(mom_, 2) +
//                                             FDProtH[1][2][2] * mom_ + FDProtH[1][2][3]);
//   } else if (dc_sec == 4) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtH[1][3][0] * pow(mom_, 3) + FDProtH[1][3][1] * pow(mom_, 2) +
//                                             FDProtH[1][3][2] * mom_ + FDProtH[1][3][3]);
//   } else if (dc_sec == 5) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtH[1][4][0] * pow(mom_, 3) + FDProtH[1][4][1] * pow(mom_, 2) +
//                                             FDProtH[1][4][2] * mom_ + FDProtH[1][4][3]);
//   } else if (dc_sec == 6) {
//     return mom_ - alpha_prot_mom_corr_FD * (FDProtH[1][5][0] * pow(mom_, 3) + FDProtH[1][5][1] * pow(mom_, 2) +
//                                             FDProtH[1][5][2] * mom_ + FDProtH[1][5][3]);
//   } else
//     return NAN;
// }

// // /// pip hadron corrections

// // float alpha_pip_mom_corr_FD[4] = {0.1, 0.15, 0.5, 0.5};
// // float alpha_pip_mom_corr_CD[3] = {0.8, 0.4, 0.8};

// double CDPip[3][4] = {{0.06775, -0.1256, -0.03055, 0.002312},
//                       {0.0484, -0.11993, 0.0746, -0.00975},
//                       {-0.00775, 0.06445, -0.04684, 0.004112}};

// float mom_corr::CD_pip_Hmom_corr(float mom_, float phi_, float alpha_pip_mom_corr_CD[3]) {
//   if (phi_ > 270 || phi_ <= 30) {
//     return mom_ - alpha_pip_mom_corr_CD[0] *
//                       (CDPip[0][0] * pow(mom_, 3) + CDPip[0][1] * pow(mom_, 2) + CDPip[0][2] * mom_ + CDPip[0][3]);
//   } else if (phi_ > 30 && phi_ <= 150) {
//     return mom_ - alpha_pip_mom_corr_CD[1] *
//                       (CDPip[1][0] * pow(mom_, 3) + CDPip[1][1] * pow(mom_, 2) + CDPip[1][2] * mom_ + CDPip[1][3]);
//   } else if (phi_ > 150 && phi_ <= 270) {
//     return mom_ - alpha_pip_mom_corr_CD[2] *
//                       (CDPip[2][0] * pow(mom_, 3) + CDPip[2][1] * pow(mom_, 2) + CDPip[2][2] * mom_ + CDPip[2][3]);
//   } else
//     return NAN;
// }
// double FDPipL[2][6][4] = {{{-0.001051, 0.004627, 0.006058, -0.01855},
//                            {-0.003084, 0.02007, -0.03488, 0.01718},
//                            {-0.001668, 0.007435, -0.0001147, -0.005516},
//                            {0.0003283, -0.005856, 0.02171, -0.012},
//                            {-0.002243, 0.01291, -0.02052, 0.012505},
//                            {-0.003408, 0.0175, -0.01814, 0.001455}},
//                           {{0.002834, -0.0171, 0.03253, -0.02928},
//                            {0.00416, -0.02376, 0.0383, -0.01701},
//                            {0.00258, -0.01698, 0.0333, -0.01591},
//                            {0.002327, -0.01192, 0.00987, 0.006634},
//                            {0.001894, -0.006947, -0.003706, 0.00831},
//                            {0.001051, -0.004684, 0.001716, -0.006424}}};

// float mom_corr::FD_pip_Hmom_corr_lower_All_FD(float mom_, float dc_sec, float alpha_pip_mom_corr_FD) {
//   if (dc_sec == 1) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipL[0][0][0] * pow(mom_, 3) + FDPipL[0][0][1] * pow(mom_, 2) +
//                                            FDPipL[0][0][2] * mom_ + FDPipL[0][0][3]);
//   } else if (dc_sec == 2) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipL[0][1][0] * pow(mom_, 3) + FDPipL[0][1][1] * pow(mom_, 2) +
//                                            FDPipL[0][1][2] * mom_ + FDPipL[0][1][3]);
//   } else if (dc_sec == 3) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipL[0][2][0] * pow(mom_, 3) + FDPipL[0][2][1] * pow(mom_, 2) +
//                                            FDPipL[0][2][2] * mom_ + FDPipL[0][2][3]);
//   } else if (dc_sec == 4) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipL[0][3][0] * pow(mom_, 3) + FDPipL[0][3][1] * pow(mom_, 2) +
//                                            FDPipL[0][3][2] * mom_ + FDPipL[0][3][3]);
//   } else if (dc_sec == 5) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipL[0][4][0] * pow(mom_, 3) + FDPipL[0][4][1] * pow(mom_, 2) +
//                                            FDPipL[0][4][2] * mom_ + FDPipL[0][4][3]);
//   } else if (dc_sec == 6) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipL[0][5][0] * pow(mom_, 3) + FDPipL[0][5][1] * pow(mom_, 2) +
//                                            FDPipL[0][5][2] * mom_ + FDPipL[0][5][3]);
//   } else
//     return NAN;
// }

// float mom_corr::FD_pip_Hmom_corr_lower_Except_All_FD(float mom_, float dc_sec, float alpha_pip_mom_corr_FD) {
//   if (dc_sec == 1) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipL[1][0][0] * pow(mom_, 3) + FDPipL[1][0][1] * pow(mom_, 2) +
//                                            FDPipL[1][0][2] * mom_ + FDPipL[1][0][3]);
//   } else if (dc_sec == 2) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipL[1][1][0] * pow(mom_, 3) + FDPipL[1][1][1] * pow(mom_, 2) +
//                                            FDPipL[1][1][2] * mom_ + FDPipL[1][1][3]);
//   } else if (dc_sec == 3) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipL[1][2][0] * pow(mom_, 3) + FDPipL[1][2][1] * pow(mom_, 2) +
//                                            FDPipL[1][2][2] * mom_ + FDPipL[1][2][3]);
//   } else if (dc_sec == 4) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipL[1][3][0] * pow(mom_, 3) + FDPipL[1][3][1] * pow(mom_, 2) +
//                                            FDPipL[1][3][2] * mom_ + FDPipL[1][3][3]);
//   } else if (dc_sec == 5) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipL[1][4][0] * pow(mom_, 3) + FDPipL[1][4][1] * pow(mom_, 2) +
//                                            FDPipL[1][4][2] * mom_ + FDPipL[1][4][3]);
//   } else if (dc_sec == 6) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipL[1][5][0] * pow(mom_, 3) + FDPipL[1][5][1] * pow(mom_, 2) +
//                                            FDPipL[1][5][2] * mom_ + FDPipL[1][5][3]);
//   } else
//     return NAN;
// }

// double FDPipH[2][6][4] = {{{-0.00287, 0.00692, 0.01993, -0.02162},
//                            {-0.0057, 0.01073, 0.0326, -0.04004},
//                            {-0.001721, 0.00743, 0.01846, -0.02194},
//                            {-0.006126, 0.00979, 0.03296, -0.03087},
//                            {-0.001957, 0.00483, 0.01576, -0.004658},
//                            {-0.001278, 0.004314, 0.01255, -0.00646}},
//                           {{-0.002602, 0.00621, 0.0158, -0.03085},
//                            {-0.00791, 0.00989, 0.0353, -0.04178},
//                            {-0.003164, 0.006046, 0.01929, -0.01709},
//                            {-0.001335, 0.002363, 0.00827, -0.004005},
//                            {-0.003061, 0.003662, 0.01438, -0.01021},
//                            {-0.001557, 0.002182, 0.006744, -0.01247}}};

// float mom_corr::FD_pip_Hmom_corr_upper_All_FD(float mom_, float dc_sec, float alpha_pip_mom_corr_FD) {
//   if (dc_sec == 1) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipH[0][0][0] * pow(mom_, 3) + FDPipH[0][0][1] * pow(mom_, 2) +
//                                            FDPipH[0][0][2] * mom_ + FDPipH[0][0][3]);
//   } else if (dc_sec == 2) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipH[0][1][0] * pow(mom_, 3) + FDPipH[0][1][1] * pow(mom_, 2) +
//                                            FDPipH[0][1][2] * mom_ + FDPipH[0][1][3]);
//   } else if (dc_sec == 3) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipH[0][2][0] * pow(mom_, 3) + FDPipH[0][2][1] * pow(mom_, 2) +
//                                            FDPipH[0][2][2] * mom_ + FDPipH[0][2][3]);
//   } else if (dc_sec == 4) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipH[0][3][0] * pow(mom_, 3) + FDPipH[0][3][1] * pow(mom_, 2) +
//                                            FDPipH[0][3][2] * mom_ + FDPipH[0][3][3]);
//   } else if (dc_sec == 5) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipH[0][4][0] * pow(mom_, 3) + FDPipH[0][4][1] * pow(mom_, 2) +
//                                            FDPipH[0][4][2] * mom_ + FDPipH[0][4][3]);
//   } else if (dc_sec == 6) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipH[0][5][0] * pow(mom_, 3) + FDPipH[0][5][1] * pow(mom_, 2) +
//                                            FDPipH[0][5][2] * mom_ + FDPipH[0][5][3]);
//   } else
//     return NAN;
// }

// float mom_corr::FD_pip_Hmom_corr_upper_Except_All_FD(float mom_, float dc_sec, float alpha_pip_mom_corr_FD) {
//   if (dc_sec == 1) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipH[1][0][0] * pow(mom_, 3) + FDPipH[1][0][1] * pow(mom_, 2) +
//                                            FDPipH[1][0][2] * mom_ + FDPipH[1][0][3]);
//   } else if (dc_sec == 2) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipH[1][1][0] * pow(mom_, 3) + FDPipH[1][1][1] * pow(mom_, 2) +
//                                            FDPipH[1][1][2] * mom_ + FDPipH[1][1][3]);
//   } else if (dc_sec == 3) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipH[1][2][0] * pow(mom_, 3) + FDPipH[1][2][1] * pow(mom_, 2) +
//                                            FDPipH[1][2][2] * mom_ + FDPipH[1][2][3]);
//   } else if (dc_sec == 4) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipH[1][3][0] * pow(mom_, 3) + FDPipH[1][3][1] * pow(mom_, 2) +
//                                            FDPipH[1][3][2] * mom_ + FDPipH[1][3][3]);
//   } else if (dc_sec == 5) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipH[1][4][0] * pow(mom_, 3) + FDPipH[1][4][1] * pow(mom_, 2) +
//                                            FDPipH[1][4][2] * mom_ + FDPipH[1][4][3]);
//   } else if (dc_sec == 6) {
//     return mom_ - alpha_pip_mom_corr_FD * (FDPipH[1][5][0] * pow(mom_, 3) + FDPipH[1][5][1] * pow(mom_, 2) +
//                                            FDPipH[1][5][2] * mom_ + FDPipH[1][5][3]);
//   } else
//     return NAN;
// }

// // /// pim hadron corrections
// // float alpha_pim_mom_corr_FD[4] = {0.5, 0.15, 0.3, 0.3};
// // float alpha_pim_mom_corr_CD[3] = {0.5, 1.0, 0.5};

// // float alpha_pim_mom_corr_FD[4] = {0.0, 0.0, 0.0, 0.0};
// // float alpha_pim_mom_corr_CD[3] = {0., 0.0, 0.0};

// // double CDPim[3][4] = {
// //     {0.0531, -0.0899, 0.05328, -0.01124}, {0.02277, -0.02846, 0.04657, -0.005}, {0.05997, -0.1099, 0.0093,
// //     0.006657}}; // becareful plot has 1st sector in 3rd place
// double CDPim[3][5] = {{-0.06088, 0.2715, -0.355, 0.1799, -0.03076},
//                       {-0.01833, 0.08844, -0.1082, 0.08466, -0.01088},
//                       {-0.1163, 0.4768, -0.616, 0.2512, -0.03062}};
// float mom_corr::CD_pim_Hmom_corr(float mom_, float phi_, float alpha_pim_mom_corr_CD[3]) {
//   if (phi_ > 270 || phi_ <= 30) {
//     return mom_ - alpha_pim_mom_corr_CD[0] * (CDPim[0][0] * pow(mom_, 4) + CDPim[0][1] * pow(mom_, 3) +
//                                               CDPim[0][2] * pow(mom_, 2) + CDPim[0][3] * mom_ + CDPim[0][4]);
//   } else if (phi_ > 30 && phi_ <= 150) {
//     return mom_ - alpha_pim_mom_corr_CD[1] * (CDPim[1][0] * pow(mom_, 4) + CDPim[1][1] * pow(mom_, 3) +
//                                               CDPim[1][2] * pow(mom_, 2) + CDPim[1][3] * mom_ + CDPim[1][4]);
//   } else if (phi_ > 150 && phi_ <= 270) {
//     return mom_ - alpha_pim_mom_corr_CD[2] * (CDPim[2][0] * pow(mom_, 4) + CDPim[2][1] * pow(mom_, 3) +
//                                               CDPim[2][2] * pow(mom_, 2) + CDPim[2][3] * mom_ + CDPim[2][4]);
//   } else
//     return NAN;
// }

// double FDPimL[2][6][4] = {{{-0.003864, 0.0297, -0.0769, 0.06757},
//                            {-0.00433, 0.03613, -0.1026, 0.0982},
//                            {0.00489, -0.02971, 0.04425, 0.002811},
//                            {-2.235e-05, 0.004494, -0.03032, 0.04724},
//                            {0.00835, -0.0581, 0.1112, -0.0421},
//                            {0.003828, -0.01945, 0.01244, 0.01949}},
//                           {{-0.00321, 0.0282, -0.08203, 0.06024},
//                            {0.003355, -0.02539, 0.05545, -0.03723},
//                            {0.00402, -0.02834, 0.06, -0.0354},
//                            {0.000703, -0.003294, 0.00447, -0.009796},
//                            {-0.002449, 0.0176, -0.03662, 0.00562},
//                            {-0.003347, 0.0308, -0.0873, 0.04517}}};

// float mom_corr::FD_pim_Hmom_corr_lower_All_FD(float mom_, float dc_sec, float alpha_pim_mom_corr_FD) {
//   if (dc_sec == 1) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimL[0][0][0] * pow(mom_, 3) + FDPimL[0][0][1] * pow(mom_, 2) +
//                                            FDPimL[0][0][2] * mom_ + FDPimL[0][0][3]);
//   } else if (dc_sec == 2) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimL[0][1][0] * pow(mom_, 3) + FDPimL[0][1][1] * pow(mom_, 2) +
//                                            FDPimL[0][1][2] * mom_ + FDPimL[0][1][3]);
//   } else if (dc_sec == 3) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimL[0][2][0] * pow(mom_, 3) + FDPimL[0][2][1] * pow(mom_, 2) +
//                                            FDPimL[0][2][2] * mom_ + FDPimL[0][2][3]);
//   } else if (dc_sec == 4) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimL[0][3][0] * pow(mom_, 3) + FDPimL[0][3][1] * pow(mom_, 2) +
//                                            FDPimL[0][3][2] * mom_ + FDPimL[0][3][3]);
//   } else if (dc_sec == 5) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimL[0][4][0] * pow(mom_, 3) + FDPimL[0][4][1] * pow(mom_, 2) +
//                                            FDPimL[0][4][2] * mom_ + FDPimL[0][4][3]);
//   } else if (dc_sec == 6) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimL[0][5][0] * pow(mom_, 3) + FDPimL[0][5][1] * pow(mom_, 2) +
//                                            FDPimL[0][5][2] * mom_ + FDPimL[0][5][3]);
//   } else
//     return NAN;
// }

// float mom_corr::FD_pim_Hmom_corr_lower_Except_All_FD(float mom_, float dc_sec, float alpha_pim_mom_corr_FD) {
//   if (dc_sec == 1) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimL[1][0][0] * pow(mom_, 3) + FDPimL[1][0][1] * pow(mom_, 2) +
//                                            FDPimL[1][0][2] * mom_ + FDPimL[1][0][3]);
//   } else if (dc_sec == 2) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimL[1][1][0] * pow(mom_, 3) + FDPimL[1][1][1] * pow(mom_, 2) +
//                                            FDPimL[1][1][2] * mom_ + FDPimL[1][1][3]);
//   } else if (dc_sec == 3) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimL[1][2][0] * pow(mom_, 3) + FDPimL[1][2][1] * pow(mom_, 2) +
//                                            FDPimL[1][2][2] * mom_ + FDPimL[1][2][3]);
//   } else if (dc_sec == 4) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimL[1][3][0] * pow(mom_, 3) + FDPimL[1][3][1] * pow(mom_, 2) +
//                                            FDPimL[1][3][2] * mom_ + FDPimL[1][3][3]);
//   } else if (dc_sec == 5) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimL[1][4][0] * pow(mom_, 3) + FDPimL[1][4][1] * pow(mom_, 2) +
//                                            FDPimL[1][4][2] * mom_ + FDPimL[1][4][3]);
//   } else if (dc_sec == 6) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimL[1][5][0] * pow(mom_, 3) + FDPimL[1][5][1] * pow(mom_, 2) +
//                                            FDPimL[1][5][2] * mom_ + FDPimL[1][5][3]);
//   } else
//     return NAN;
// }

// double FDPimH[2][6][4] = {{{-0.0002575, -0.002295, -0.003235, 0.01292},
//                            {-0.00398, 0.0052, 0.02222, -0.02208},
//                            {-0.003052, 0.003235, 0.01591, -0.01168},
//                            {0.001004, -0.003637, -0.00909, 0.02184},
//                            {-0.000769, -0.0006876, 0.0008802, 0.002695},
//                            {0.002268, -0.003813, -0.013794, 0.0193}},
//                           {{0.001858, -0.003899, -0.01464, 0.009155},
//                            {0.0002279, -0.000515, -0.002798, -0.003864},
//                            {-0.001788, 0.00237, 0.008934, -0.01624},
//                            {-0.002106, 0.002811, 0.0091, -0.0273},
//                            {-0.003132, 0.00478, 0.01424, -0.04587},
//                            {0.001281, -0.001394, -0.01089, -0.0178}}};

// float mom_corr::FD_pim_Hmom_corr_upper_All_FD(float mom_, float dc_sec, float alpha_pim_mom_corr_FD) {
//   if (dc_sec == 1) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimH[0][0][0] * pow(mom_, 3) + FDPimH[0][0][1] * pow(mom_, 2) +
//                                            FDPimH[0][0][2] * mom_ + FDPimH[0][0][3]);
//   } else if (dc_sec == 2) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimH[0][1][0] * pow(mom_, 3) + FDPimH[0][1][1] * pow(mom_, 2) +
//                                            FDPimH[0][1][2] * mom_ + FDPimH[0][1][3]);
//   } else if (dc_sec == 3) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimH[0][2][0] * pow(mom_, 3) + FDPimH[0][2][1] * pow(mom_, 2) +
//                                            FDPimH[0][2][2] * mom_ + FDPimH[0][2][3]);
//   } else if (dc_sec == 4) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimH[0][3][0] * pow(mom_, 3) + FDPimH[0][3][1] * pow(mom_, 2) +
//                                            FDPimH[0][3][2] * mom_ + FDPimH[0][3][3]);
//   } else if (dc_sec == 5) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimH[0][4][0] * pow(mom_, 3) + FDPimH[0][4][1] * pow(mom_, 2) +
//                                            FDPimH[0][4][2] * mom_ + FDPimH[0][4][3]);
//   } else if (dc_sec == 6) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimH[0][5][0] * pow(mom_, 3) + FDPimH[0][5][1] * pow(mom_, 2) +
//                                            FDPimH[0][5][2] * mom_ + FDPimH[0][5][3]);
//   } else
//     return NAN;
// }

// float mom_corr::FD_pim_Hmom_corr_upper_Except_All_FD(float mom_, float dc_sec, float alpha_pim_mom_corr_FD) {
//   if (dc_sec == 1) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimH[1][0][0] * pow(mom_, 3) + FDPimH[1][0][1] * pow(mom_, 2) +
//                                            FDPimH[1][0][2] * mom_ + FDPimH[1][0][3]);
//   } else if (dc_sec == 2) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimH[1][1][0] * pow(mom_, 3) + FDPimH[1][1][1] * pow(mom_, 2) +
//                                            FDPimH[1][1][2] * mom_ + FDPimH[1][1][3]);
//   } else if (dc_sec == 3) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimH[1][2][0] * pow(mom_, 3) + FDPimH[1][2][1] * pow(mom_, 2) +
//                                            FDPimH[1][2][2] * mom_ + FDPimH[1][2][3]);
//   } else if (dc_sec == 4) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimH[1][3][0] * pow(mom_, 3) + FDPimH[1][3][1] * pow(mom_, 2) +
//                                            FDPimH[1][3][2] * mom_ + FDPimH[1][3][3]);
//   } else if (dc_sec == 5) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimH[1][4][0] * pow(mom_, 3) + FDPimH[1][4][1] * pow(mom_, 2) +
//                                            FDPimH[1][4][2] * mom_ + FDPimH[1][4][3]);
//   } else if (dc_sec == 6) {
//     return mom_ - alpha_pim_mom_corr_FD * (FDPimH[1][5][0] * pow(mom_, 3) + FDPimH[1][5][1] * pow(mom_, 2) +
//                                            FDPimH[1][5][2] * mom_ + FDPimH[1][5][3]);
//   } else
//     return NAN;
// }

// ////////////// ###########  These are the corrections we used for our w< 2.55 , FD seperated in theta angle cases, as
// /// June 04 2023   DONE AT THIS POINT

////// oour final hadron momentum  dp corrections finished

// Below shows how the corrections are to be applied using the ROOT momentum 4-vector using the above code:
// auto fe = dppC(ex, ey, ez, esec, 0) + 1;
// auto fpip = dppC(pipx, pipy, pipz, pipsec, 1) + 1;
// auto fpim = dppC(pimx, pimy, pimz, pimsec, 2) + 1;
// auto fpro = dppC(prox, proy, proz, prosec, 3) + 1;

// auto eleC = ROOT::Math::PxPyPzMVector(ex * fe, ey* fe, ez* fe, 0);
// auto pipC = ROOT::Math::PxPyPzMVector(pipx * fpip, pipy* fpip, pipz* fpip, 0.13957);
// auto pimC = ROOT::Math::PxPyPzMVector(pimx * fpim, pimy* fpim, pimz* fpim, 0.13957);
// auto proC = ROOT::Math::PxPyPzMVector(prox * fpro, proy* fpro, proz* fpro, 0.938);

////////////////// new new mom corr done (Aug-15-2022)

// }  // namespace mom_corr

// ////////////////// old mom corrections (probably better one)

// // momentum corrections earlier

// double xx[54] = {
//     0.0263375, 0.0158871,  0.0130852,  -0.00366006, 0.00694866,  0.0197195, 0.00767067, 0.00480921,  -0.0175756,
//     0.0252757, 0.0156601,  0.00984872, 0.00244435,  0.00681414,  0.0294068, 0.0059881,  0.00286992,  0.0179319,
//     0.0171495, 0.00359637, -0.0046115, 0.00314739,  0.0136338,   0.0768753, 0.00675454, -0.0118234,  -0.0288654,
//     0.0189465, 0.0131816,  0.0262004,  0.00375165,  0.00907457,  0.0486894, 0.00806305, 0.0006999,   0.00527513,
//     0.0116485, 0.0105681,  0.0149848,  0.000318094, -0.00480124, 0.0395545, 0.00824216, -0.00070659, -0.0057075,
//     0.0213057, 0.0112999,  0.0100216,  0.000653685, 0.0093174,   0.0822385, 0.00808384, 0.000898799, -0.0172692,
// };
// double pars[6][3][3];
// int ipar = 0;

// /// this was inside the constructor   // for (int isec_mom_corr = 0; isec_mom_corr < 6; isec_mom_corr++) {
//   for (int ivec = 0; ivec < 3; ivec++) {
//     double dp1 = xx[ipar++], dp5 = xx[ipar++], dp9 = xx[ipar++];

//     pars[isec_mom_corr][ivec][0] = (dp1 - 2 * dp5 + dp9) / 32.;
//     pars[isec_mom_corr][ivec][1] = (-7 * dp1) / 16. + (5 * dp5) / 8. - (3 * dp9) / 16.;
//     pars[isec_mom_corr][ivec][2] = (45 * dp1) / 32. - (9 * dp5) / 16. + (5 * dp9) / 32.;
//   }
// }

// /// now outside the constructor

// double Reaction::dpp(float px, float py, float pz, int sec_mom_corr, int ivec) {
//   double pp = sqrt(px * px + py * py + pz * pz);

//   double a = pars[sec_mom_corr - 1][ivec][0], b = pars[sec_mom_corr - 1][ivec][1],
//          c = pars[sec_mom_corr - 1][ivec][2];

//   // double dp = a * pp * pp + b * pp + c;  // pol2 corr func

//   // electron pol1 corr func for each sec_mom_corr and each phi bins
//   if (ivec == 0) {
//     if (sec_mom_corr == 1) {
//       dp = 0.45 * b * (pp - 9) + 0.1 * c;

//       // ep 3 phi bins
//       // dp = -0.01*b*(pp-9)+1.35*c; //phi<-5
//       // dp = 0.6*b*(pp-9)-0.3*c; //-5<phi<5
//       // dp = 1.7*b*(pp-9)-1.5*c; //phi>5
//     }
//     if (sec_mom_corr == 2) {
//       dp = -0.15 * b * (pp - 8.0) - 0.3 * c;

//       // ep 3 phi bins
//       // dp = -0.7*b*(pp-8.0)+0.4*c; //phi<-5
//       // dp = -0.05*b*(pp-8.0)-0.4*c; //-5<phi<5
//       // dp = 0.01*b*(pp-8.0)-1.5*c; //phi>5
//     }
//     if (sec_mom_corr == 3) {
//       dp = 3. * b * (pp - 5.4) - 0.5 * c;

//       // ep 3 phi bins
//       // dp = 0.04*b*(pp-5.4)-3.5*c; //phi<-5
//       // dp = 0.06*b*(pp-5.4)-3.*c; //-5<phi<5
//       // dp = 1.1*b*(pp-5.4)-0.7*c; //phi>5
//     }
//     if (sec_mom_corr == 4) {
//       dp = 0.25 * b * (pp - 9.25) - 0.3 * c;

//       // ep 3 phi bins
//       // dp = 0.25*b*(pp-9.25)-0.7*c; //phi<-5
//       // dp = 0.25*b*(pp-9.25)+0.05*c; //-5<phi<5
//       // dp = 0.1*b*(pp-9.25)+1.1*c; //phi>5
//     }
//     if (sec_mom_corr == 5) {
//       dp = 2.2 * b * (pp - 7.5) - 0.5 * c;

//       // ep 3 phi bins
//       // dp = 2.2*b*(pp-7.5)+0.5*c; //phi<-5
//       // dp = 2.2*b*(pp-7.5)-0.1*c; //-5<phi<5
//       // dp = 2.2*b*(pp-7.5)-0.6*c; //phi>5
//     }
//     if (sec_mom_corr == 6) {
//       dp = 0.5 * b * (pp - 7) - 0.6 * c;

//       // ep 3 phi bins
//       // dp = 1.263*b*(pp-7)+0.5*c; //phi<-5
//       // dp = 1.*b*(pp-7)-0.5*c; //-5<phi<5
//       // dp = 0.5*b*(pp-7)-1.45*c; //phi>5
//     }
//   }
//   return dp / pp;
// };

// double fe = dpp(ex, ey, ez, esec, 0) + 1;
// double fpip = dpp(pipx,pipy,pipz,pipsec,1) + 1;
// double fpim = dpp(pimx,pimy,pimz,pimsec,2) + 1;

///// old momentum corrections done!!!!!!!!

// my previous 2d in cd and 1d in fd, mom part only energy loss corrections: .............
// _Energy_loss_uncorr_prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
// _prot_status = abs(_data->status(i));

// _prot_mom_uncorr = _Energy_loss_uncorr_prot->P();
// _prot_theta = _Energy_loss_uncorr_prot->Theta() * 180 / PI;
// if (abs(_data->status(i)) < 4000) {
//   _sectorProt = _data->dc_sec(i);
//   if (_prot_theta <= 27) {
//     _E_corr_val_prot = -0.00078846 * pow(_prot_mom_uncorr, 5) + 0.0093734 * pow(_prot_mom_uncorr, 4) -
//                        0.04277868 * pow(_prot_mom_uncorr, 3) + 0.09421284 * pow(_prot_mom_uncorr, 2) -
//                        0.10095842 * (_prot_mom_uncorr) + 0.04567203;
//   } else {
//     _E_corr_val_prot = -0.0023389 * pow(_prot_mom_uncorr, 5) + 0.02838603 * pow(_prot_mom_uncorr, 4) -
//                        0.13214962 * pow(_prot_mom_uncorr, 3) + 0.29609571 * pow(_prot_mom_uncorr, 2) -
//                        0.32307424 * (_prot_mom_uncorr) + 0.14742569;
//   }
// } else if (abs(_data->status(i)) >= 4000) {
//   _E_corr_val_prot = 0.0;
//   // _E_corr_val_prot = ((-9.30990933e-05) * pow(_prot_theta, 3) + (1.23584235e-02) * pow(_prot_theta, 2) +
//   //                     (-5.42538215e-01) * (_prot_theta) + 7.87921215e+00) *
//   //                        pow(_prot_mom_uncorr, 3) +

//   //                    (4.17955911e-04 * pow(_prot_theta, 3) + (-5.53676478e-02) * pow(_prot_theta, 2) +
//   //                     (2.42642631e+00) * (_prot_theta) + (-3.51829220e+01)) *
//   //                        pow(_prot_mom_uncorr, 2) +

//   //                    ((-5.58084320e-04) * pow(_prot_theta, 3) + (7.38670367e-02) * pow(_prot_theta, 2) +
//   //                     (-3.23723227e+00) * (_prot_theta) + 4.69456718e+01) *
//   //                        (_prot_mom_uncorr) +

//   //                    ((2.40014720e-04) * pow(_prot_theta, 3) + (-3.17071405e-02) * pow(_prot_theta, 2) +
//   //                     (1.38769727e+00 * (_prot_theta)) + (-2.01072704e+01));
// }

// _prot_mom_tmt = _prot_mom_uncorr + _E_corr_val_prot;

// _px_prime_prot_E = _data->px(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
// _py_prime_prot_E = _data->py(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
// _pz_prime_prot_E = _data->pz(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));

// // _prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);

// /* pip corrections 1d in fd and 2d in cd mom part only
//     // //   // _pip_status = abs(_data->status(i));
//     // _Energy_loss_uncorr_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
//     // _pip_mom_uncorr = _Energy_loss_uncorr_pip->P();
//     // _pip_theta = _Energy_loss_uncorr_pip->Theta() * 180 / PI;
//     // if (abs(_data->status(i)) < 4000) {
//     //   // _sectorPip = _data->dc_sec(i);
//     //   // fpip = dppC(_data->px(i), _data->py(i), _data->pz(i), _data->dc_sec(i), 1) + 1;

//     //   if (_pip_theta <= 27) {
//     //     _E_corr_val_pip = 9.21970527e-05 * pow(_pip_mom_uncorr, 3) - 3.70500143e-04 * pow(_pip_mom_uncorr, 2) +
//     //                       2.78880101e-04 * (_pip_mom_uncorr) + 2.66040566e-03;

//     //   } else {
//     //     _E_corr_val_pip = -0.00010482 * pow(_pip_mom_uncorr, 3) + 0.00080463 * pow(_pip_mom_uncorr, 2) -
//     //                       0.0022871 * (_pip_mom_uncorr) + 0.00831496;
//     //   }
//     // } else if (abs(_data->status(i)) >= 4000) {
//     //   _E_corr_val_pip = 0.0;

//     //   // _E_corr_val_pip = (-6.50509539e-07 * pow(_pip_theta, 3) + 1.31547371e-04 * pow(_pip_theta, 2) +
//     //   //                    (-7.99024673e-03) * (_pip_theta) + 1.60563630e-01) *
//     //   //                       pow(_pip_mom_uncorr, 3) +

//     //   //                   (2.48202211e-06 * pow(_pip_theta, 3) + (-5.15757241e-04) * pow(_pip_theta, 2) +
//     //   //                    3.19833135e-02 * (_pip_theta) + (-6.53476057e-01)) *
//     //   //                       pow(_pip_mom_uncorr, 2) +

//     //   //                   (-2.71923009e-06 * pow(_pip_theta, 3) + 5.80375203e-04 * pow(_pip_theta, 2) +
//     //   //                    (-3.75941898e-02) * (_pip_theta) + 7.80443724e-01) *
//     //   //                       (_pip_mom_uncorr) +

//     //   //                   4.62456800e-07 * pow(_pip_theta, 3) + (-1.08401698e-04) * pow(_pip_theta, 2) +
//     //   //                   8.09261138e-03 * (_pip_theta)-2.05315604e-01;
//     // }
//     // _pip_mom_tmt = _pip_mom_uncorr + _E_corr_val_pip;

//     // _px_prime_pip_E = _data->px(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
//     // _py_prime_pip_E = _data->py(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
//     // _pz_prime_pip_E = _data->pz(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));

//     // // _pip->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);

/*......................
      _pim_status = abs(_data->status(i));
      _Energy_loss_uncorr_pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
      _pim_mom_uncorr = _Energy_loss_uncorr_pim->P();
      _pim_theta = _Energy_loss_uncorr_pim->Theta() * 180 / PI;

      // // this is for energy loss corrections

      if (abs(_data->status(i)) < 4000) {
        _sectorPim = _data->dc_sec(i);
        // fpim = dppC(_data->px(i), _data->py(i), _data->pz(i), _data->dc_sec(i), 2) + 1;

        if (_pim_theta <= 27) {
          _E_corr_val_pim = -0.00035275 * pow(_pim_mom_uncorr, 3) + 0.00291237 * pow(_pim_mom_uncorr, 2) -
                            0.00681058 * (_pim_mom_uncorr) + 0.00736721;

        } else {
          _E_corr_val_pim = 0.00019358 * pow(_pim_mom_uncorr, 3) - 0.00103456 * pow(_pim_mom_uncorr, 2) +
                            0.00024772 * (_pim_mom_uncorr) + 0.00735159;
        }
      } else if (abs(_data->status(i)) >= 4000) {
        _E_corr_val_pim = 0.0;

        // // _E_corr_val_pim = (0.02153442) * pow(_pim_mom_uncorr, 5) -
        // //                   (0.13271424) * pow(_pim_mom_uncorr, 4) +
        // //                   (0.27140262) * pow(_pim_mom_uncorr, 3) -
        // //                   (0.23266059) * pow(_pim_mom_uncorr, 2) +
        // //                   (0.04031421) * (_pim_mom_uncorr) + 0.0036634;

        // _E_corr_val_pim = (-4.94426765e-07 * pow(_pim_theta, 3) + 9.85729368e-05 * pow(_pim_theta, 2) +
        //                    (-5.85778699e-03) * (_pim_theta) + 1.17447168e-01) *
        //                       pow(_pim_mom_uncorr, 3) +

        //                   (1.75953956e-06 * pow(_pim_theta, 3) + (-3.63382515e-04) * pow(_pim_theta, 2) +
        //                    2.21447425e-02 * (_pim_theta) + (-4.54844509e-01)) *
        //                       pow(_pim_mom_uncorr, 2) +

        //                   (-1.90446515e-06 * pow(_pim_theta, 3) + 4.08768480e-04 * pow(_pim_theta, 2) +
        //                    (-2.65277055e-02) * (_pim_theta) + 5.57286393e-01) *
        //                       (_pim_mom_uncorr) +

        //                   2.05653097e-07 * pow(_pim_theta, 3) + (-5.44018546e-05) * pow(_pim_theta, 2) +
        //                   4.61561853e-03 * (_pim_theta)-1.35303212e-01;
      }

      // _pim_mom = _pim_mom_uncorr + _E_corr_val_pim; // first iteration

      _pim_mom_tmt = _pim_mom_uncorr + _E_corr_val_pim;  // first iteration

      _px_prime_pim_E = _data->px(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
      _py_prime_pim_E = _data->py(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
      _pz_prime_pim_E = _data->pz(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));

      // _pim->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, MASS_PIM);  // energy loss corrected
*/

//////
////..........
// // Sangbaek energy loss corrections parameters for theta angle of proton

// float A_th(float mom_, float theta_, float theta_DCr1_p, int dc_sec) {
//   if (dc_sec >= 1 && dc_sec <= 6) {
//     if (theta_DCr1_p < -53.14680163254601 + 79.61307254040804 * pow((mom_ - 0.3), 0.05739232362022314)) {
//       return -0.16742969 + 0.00697925 * theta_;
//       //   Dθ = − 0.16742969 + 0.00697925 × θ
//     } else
//       return 2.04334532 * 10 - 1.81052405 * theta_ + 5.32556360e-2 * theta_ * theta_ -
//              5.23157558e-4 * theta_ * theta_ * theta_;
//     //  Dθ = 2.04334532 × 10 − 1.81052405 × θ + 5.32556360 × 10−2 × θ2 − 5.23157558 × 10−4 × θ3
//   } else
//     return -1.09849291e2 + 8.86664014 * theta_ - 0.26643881 * theta_ * theta_ +
//            3.53814210e-3 * theta_ * theta_ * theta_ - 1.75297107e-5 * theta_ * theta_ * theta_ * theta_;

//   //    Aθ = − 1.09849291 × 102 + 8.86664014 × θ − 0.26643881 × θ2 + 3.53814210 × 10−3 ∗ θ3 − 1.75297107 × 10−5 ×
//   θ4
// }

// float B_th(float mom_, float theta_, float theta_DCr1_p, int dc_sec) {
//   if (dc_sec >= 1 && dc_sec <= 6) {
//     if (theta_DCr1_p < -53.14680163254601 + 79.61307254040804 * pow((mom_ - 0.3), 0.05739232362022314)) {
//       return 0.23352115 - 0.01338697 * theta_;
//       //  Eθ = 0.23352115 − 0.01338697 × θ
//     } else
//       return 8.74233279 - 7.63869344e-1 * theta_ + 2.22376362e-2 * theta_ * theta_ -
//              2.16457260e-4 * theta_ * theta_ * theta_;
//     // Eθ = 8.74233279 − 7.63869344 × 10−1 × θ + 2.22376362 × 10−2 × θ2 − 2.16457260 × 10−4 × θ3

//   } else
//     return 9.52034523e2 - 5.74808292e1 * theta_ + 1.15386949 * theta_ * theta_ -
//            7.57970373e-3 * theta_ * theta_ * theta_;
//   //      Bθ = 9.52034523 × 102 − 5.74808292 × 10 × θ +1.15386949 × θ2 − 7.57970373 × 10−3 × θ3
// }
// float C_th(float mom_, float theta_, int dc_sec) {
//   if (dc_sec < 1 || dc_sec > 6) {
//     return -2.00387313e2 + 1.18979079e1 * theta_ - 2.37730217e-1 * theta_ * theta_ +
//            1.55153003e-3 * theta_ * theta_ * theta_;
//     // Cθ = − 2.00387313 × 102 + 1.18979079 × 10 × θ − 2.37730217 × 10−1 × θ2 + 1.55153003 × 10−3 × θ3

//   } else
//     return NAN;
// }
// // energy loss corrections parameters for phi angle of proton

// float A_ph(float mom_, float theta_, float theta_DCr1_p, int dc_sec) {
//   if (dc_sec >= 1 && dc_sec <= 6) {
//     if (theta_DCr1_p < -53.14680163254601 + 79.61307254040804 * pow((mom_ - 0.3), 0.05739232362022314)) {
//       return 0.21192125 - 0.0115175 * theta_;
//       //   Dφ = 0.21192125 − 0.0115175 × θ
//     } else
//       return 0.54697831 - 0.04896981 * theta_ + 0.00111376 * theta_ * theta_;
//     // Aφ = 0.54697831 − 0.04896981 × θ + 0.00111376 × θ2

//   } else
//     return 4.94546178 - 3.26662886e-1 * theta_ + 7.39069603e-3 * theta_ * theta_ -
//            6.83599356e-5 * theta_ * theta_ * theta_ + 2.12303103e-7 * theta_ * theta_ * theta_ * theta_;
//   //    Aφ = 4.94546178 − 3.26662886 × 10−1 × θ + 7.39069603 × 10−3 × θ2 − 6.83599356 × 10−5 × θ3 +
//   //   2.12303103 × 10−7 × θ4;
// }

// float B_ph(float mom_, float theta_, float theta_DCr1_p, int dc_sec) {
//   if (dc_sec >= 1 && dc_sec <= 6) {
//     if (theta_DCr1_p < -53.14680163254601 + 79.61307254040804 * pow((mom_ - 0.3), 0.05739232362022314)) {
//       return -8.94307411e-1 + 1.66349766e-1 * theta_ - 8.90617559e-3 * theta_ * theta_ +
//              1.64803754e-4 * theta_ * theta_ * theta_;
//       //   Eφ = − 8.94307411 × 10−1 + 1.66349766 × 10−1 × θ − 8.90617559 × 10−3 × θ2 + 1.64803754 × 10−4 × θ3
//     } else
//       return -4.06733541e2 + 2.43696202e1 * theta_ - 3.36144736e-1 * theta_ * theta_;
//     // Bφ = − 4.06733541 × 102 + 2.43696202 × 10 × θ − 3.36144736 × 10−1 × θ2;
//   } else
//     return 1.72181613e5 - 1.36827111e4 * theta_ + 4.00923146e2 * theta_ * theta_ -
//            5.12792347 * theta_ * theta_ * theta_ + 2.41793167e-2 * theta_ * theta_ * theta_ * theta_;
//   //   Bφ = 1.72181613 × 105 − 1.36827111 × 104 × θ + 4.00923146 × 102 × θ2 − 5.12792347 × θ3 + 2.41793167 × 10−2 ×
//   //   θ4;
// }
// float C_ph(float mom_, float theta_, float theta_DCr1_p, int dc_sec) {
//   if (dc_sec >= 1 && dc_sec <= 6) {
//     return 2.06378660e1 - 1.42866062 * theta_ + 2.01085440e-2 * theta_ * theta_;
//     //    Cφ = 2.06378660 × 10 − 1.42866062 × θ + 2.01085440 × 10−2 × θ2;

//   } else
//     return 1.20477219e2 - 5.86630228 * theta_ + 7.44007875e-2 * theta_ * theta_ -
//            2.42652473e-4 * theta_ * theta_ * theta_;
//   // Cφ = 1.20477219 × 102 − 5.86630228 × θ + 7.44007875 × 10−2 × θ2 − 2.42652473 × 10−4 × θ3;
// }
// }

// // Here are the functions used to do the energy loss corrections
// // pnew =p + Ap + Bp/p;
// // θnew = θ + Dθ + Eθ / p2; or = θ+Aθ +Bθ ×exp(Cθp)
// // φnew = φ + Aφ + Bφ ×exp(Cφ ×p).;
// // if (_is_FD && _prot_mom_uncorr >= 1.0) {
// //     // these are Andrey's corrections
// //     if (_is_lower_band) /// it is not lower band it is for less than 27 degree.......... be careful here.
// //       _prot_mom_tmt = _prot_mom_uncorr + exp(-2.739 - 3.932 * _prot_theta_uncorr) + 0.002907;
// //     else
// //       _prot_mom_tmt = _prot_mom_uncorr + exp(-1.2 - 4.228 * _prot_mom_uncorr) + 0.007502;
// // } else {
// // if (_is_CD || (_is_FD && _prot_mom_uncorr < 1.0)) {
// _prot_mom_tmt = _prot_mom_uncorr +
//                 mom_corr::A_p(_prot_mom_uncorr, _prot_theta_uncorr, _thetaDC_r1_Prot, _sectorProt) +
//                 mom_corr::B_p(_prot_mom_uncorr, _prot_theta_uncorr, _thetaDC_r1_Prot, _sectorProt) /
//                 _prot_mom_uncorr;
// // }

// if (_is_FD) {
//   std::cout <<  " ststus fd " << abs(_data->status(i)) << std::endl;

//   _prot_theta_tmt = _prot_theta_uncorr +
//                     mom_corr::A_th(_prot_mom_uncorr, _prot_theta_uncorr, _thetaDC_r1_Prot, _sectorProt) +
//                     (mom_corr::B_th(_prot_mom_uncorr, _prot_theta_uncorr, _thetaDC_r1_Prot, _sectorProt) /
//                         _prot_mom_uncorr * _prot_mom_uncorr);
// } else if (_is_CD) {
//   std::cout << " ststus cd " << abs(_data->status(i)) << std::endl;

//   _prot_theta_tmt = _prot_theta_uncorr +
//                     mom_corr::A_th(_prot_mom_uncorr, _prot_theta_uncorr, _thetaDC_r1_Prot, _sectorProt) +
//                     mom_corr::B_th(_prot_mom_uncorr, _prot_theta_uncorr, _thetaDC_r1_Prot, _sectorProt) *
//                         exp(mom_corr::C_th(_prot_mom_uncorr, _prot_theta_uncorr, _sectorProt) *
//                         _prot_mom_uncorr);
// }

// if (_is_lower_band) {
//   // std::cout << " dc theta lower band " << _thetaDC_r1_Prot << " ststus " << abs(_data->status(i)) <<
//   std::endl; _prot_phi_tmt = _prot_phi_uncorr + mom_corr::A_ph(_prot_mom_uncorr, _prot_theta_uncorr,
//   _thetaDC_r1_Prot, _sectorProt) +
//       mom_corr::B_ph(_prot_mom_uncorr, _prot_theta_uncorr, _thetaDC_r1_Prot, _sectorProt) /
//           _prot_mom_uncorr* _prot_mom_uncorr;
// } else if (!_is_lower_band) {
//   // std::cout << " dc theta upper band cd " << _thetaDC_r1_Prot << " ststus " << abs(_data->status(i)) <<
//   std::endl;

//   _prot_phi_tmt =
//   _prot_phi_uncorr + mom_corr::A_ph(_prot_mom_uncorr, _prot_theta_uncorr, _thetaDC_r1_Prot, _sectorProt) +
//       mom_corr::B_ph(_prot_mom_uncorr, _prot_theta_uncorr, _thetaDC_r1_Prot, _sectorProt) *
//           exp(mom_corr::C_ph(_prot_mom_uncorr, _prot_theta_uncorr, _thetaDC_r1_Prot, _sectorProt) *
//           _prot_mom_uncorr);
// }

// // // (-0.00051894 - 0.00018104 * _prot_theta_uncorr) +
// // //     (3.29466917e-3 + 5.73663160e-4 * _prot_theta_uncorr − 1.40807209e-5 * _prot_theta * _prot_theta) /
// // _prot_mom_uncorr;
// // // else if () {
// //   // Ap = − 3.03346359 × 10−1 + 1.83368163 × 10−2 × θ − 2.86486404 × 10−4 × θ2(14) Bp =
// //   //     2.01023276 × 10−1 − 1.13312215 × 10−2 × θ + 1.82487916 × 10−4 × θ2.;
// // // }
// // std::cout << " sin theta " << sinf(_prot_theta_tmt) << " cos phi " << cosf(_prot_phi_tmt)<<std::endl;

////// sangbaek corrections done!!!!!!!!!
