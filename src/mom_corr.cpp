
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

float alpha_pip_mom_corr_FD[4] = {0.1, 0.15, 0.5, 0.5};
float alpha_pip_mom_corr_CD[3] = {0.8, 0.4, 0.8};

float alpha_pim_mom_corr_FD[4] = {0.5, 0.15, 0.3, 0.3};
float alpha_pim_mom_corr_CD[3] = {0.5, 1.0, 0.5};

// float alpha_prot_mom_corr_FD[4] = {0.0, 0.0, 0.0, 0.0};
// float alpha_prot_mom_corr_CD[3] = {0.0, 0.0, 0.0};

// float alpha_pip_mom_corr_FD[4] = {0.0, 0.0, 0.0, 0.0};
// float alpha_pip_mom_corr_CD[3] = {0.0, 0.0, 0.0};

// float alpha_pim_mom_corr_FD[4] = {0.0, 0.0, 0.0, 0.0};
// float alpha_pim_mom_corr_CD[3] = {0.0, 0.0, 0.0};

// double CDProt[3][4] = {
//     {-0.05237, 0.3066, -0.5225, 0.1763}, {0.0853, -0.275, 0.2484, -0.0692}, {0.0435, -0.1208, 0.128, -0.0371}};
// 3rd order pol

double CDProt[3][5] = {{-0.2578, 1.334, -2.3, 1.489, -0.3545},
                       {-0.0736, 0.4873, -1.048, 0.862, -0.2374},
                       {-0.0928, 0.5454, -1.074, 0.874, -0.2386}};
// corrections

float mom_corr::CD_prot_Hmom_corr(float mom_, float phi_, float alpha_prot_mom_corr_CD[3]) {
  if (phi_ > 270 || phi_ <= 30) {
    return mom_ - alpha_prot_mom_corr_CD[0] * (CDProt[0][0] * pow(mom_, 4) + CDProt[0][1] * pow(mom_, 3) +
                                               CDProt[0][2] * pow(mom_, 2) + CDProt[0][3] * mom_ + CDProt[0][4]);
  } else if (phi_ > 30 && phi_ <= 150) {
    return mom_ - alpha_prot_mom_corr_CD[1] * (CDProt[1][0] * pow(mom_, 4) + CDProt[1][1] * pow(mom_, 3) +
                                               CDProt[1][2] * pow(mom_, 2) + CDProt[1][3] * mom_ + CDProt[1][4]);
  } else if (phi_ > 150 && phi_ <= 270) {
    return mom_ - alpha_prot_mom_corr_CD[2] * (CDProt[2][0] * pow(mom_, 4) + CDProt[2][1] * pow(mom_, 3) +
                                               CDProt[2][2] * pow(mom_, 2) + CDProt[2][3] * mom_ + CDProt[2][4]);
  } else
    return NAN;
}

float FDProtL[2][6][4] = {{{0.0004573, 0.000176, -0.01131, 0.011406},
                           {0.001105, -0.0095, 0.01335, 0.003487},
                           {-0.001555, 0.007084, -0.0193, 0.02263},
                           {0.000816, -0.004272, -0.00223, 0.01088},
                           {-0.0002866, 0.004208, -0.02225, 0.0223},
                           {0.00344, -0.02016, 0.01811, 0.004692}},
                          {{-0.002356, 0.01585, -0.03143, 0.001087},
                           {0.003145, -0.01888, 0.0341, -0.02017},
                           {-0.0005207, 0.002861, -0.004078, 0.001318},
                           {-0.002945, 0.0183, -0.03568, 0.01697},
                           {-0.003153, 0.02109, -0.04117, 0.01209},
                           {-0.00556, 0.0387, -0.08203, 0.02975}}};

float mom_corr::FD_prot_Hmom_corr_lower_All_FD(float mom_, float dc_sec, float alpha_prot_mom_corr_FD) {
  if (dc_sec == 1) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtL[0][0][0] * pow(mom_, 3) + FDProtL[0][0][1] * pow(mom_, 2) +
                                            FDProtL[0][0][2] * mom_ + FDProtL[0][0][3]);
  } else if (dc_sec == 2) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtL[0][1][0] * pow(mom_, 3) + FDProtL[0][1][1] * pow(mom_, 2) +
                                            FDProtL[0][1][2] * mom_ + FDProtL[0][1][3]);
  } else if (dc_sec == 3) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtL[0][2][0] * pow(mom_, 3) + FDProtL[0][2][1] * pow(mom_, 2) +
                                            FDProtL[0][2][2] * mom_ + FDProtL[0][2][3]);
  } else if (dc_sec == 4) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtL[0][3][0] * pow(mom_, 3) + FDProtL[0][3][1] * pow(mom_, 2) +
                                            FDProtL[0][3][2] * mom_ + FDProtL[0][3][3]);
  } else if (dc_sec == 5) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtL[0][4][0] * pow(mom_, 3) + FDProtL[0][4][1] * pow(mom_, 2) +
                                            FDProtL[0][4][2] * mom_ + FDProtL[0][4][3]);
  } else if (dc_sec == 6) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtL[0][5][0] * pow(mom_, 3) + FDProtL[0][5][1] * pow(mom_, 2) +
                                            FDProtL[0][5][2] * mom_ + FDProtL[0][5][3]);
  } else
    return NAN;
}

float mom_corr::FD_prot_Hmom_corr_lower_Except_All_FD(float mom_, float dc_sec, float alpha_prot_mom_corr_FD) {
  if (dc_sec == 1) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtL[1][0][0] * pow(mom_, 3) + FDProtL[1][0][1] * pow(mom_, 2) +
                                            FDProtL[1][0][2] * mom_ + FDProtL[1][0][3]);
  } else if (dc_sec == 2) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtL[1][1][0] * pow(mom_, 3) + FDProtL[1][1][1] * pow(mom_, 2) +
                                            FDProtL[1][1][2] * mom_ + FDProtL[1][1][3]);
  } else if (dc_sec == 3) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtL[1][2][0] * pow(mom_, 3) + FDProtL[1][2][1] * pow(mom_, 2) +
                                            FDProtL[1][2][2] * mom_ + FDProtL[1][2][3]);
  } else if (dc_sec == 4) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtL[1][3][0] * pow(mom_, 3) + FDProtL[1][3][1] * pow(mom_, 2) +
                                            FDProtL[1][3][2] * mom_ + FDProtL[1][3][3]);
  } else if (dc_sec == 5) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtL[1][4][0] * pow(mom_, 3) + FDProtL[1][4][1] * pow(mom_, 2) +
                                            FDProtL[1][4][2] * mom_ + FDProtL[1][4][3]);
  } else if (dc_sec == 6) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtL[1][5][0] * pow(mom_, 3) + FDProtL[1][5][1] * pow(mom_, 2) +
                                            FDProtL[1][5][2] * mom_ + FDProtL[1][5][3]);
  } else
    return NAN;
}

double FDProtH[2][6][4] = {{{-0.01047, 0.0603, -0.08997, 0.0436},
                            {-0.000767, 0.008934, -0.01267, 0.011955},
                            {-0.013275, 0.0709, -0.1035, 0.0569},
                            {0.003572, -0.007786, 0.003706, 0.01051},
                            {-0.00461, 0.03168, -0.05597, 0.0392},
                            {-0.00398, 0.02565, -0.04474, 0.03076}},
                           {{-0.006363, 0.0391, -0.06036, 0.0166},
                            {0.00442, -0.02858, 0.06744, -0.04587},
                            {0.000977, -0.00703, 0.02438, -0.00925},
                            {0.003227, -0.01921, 0.04, -0.01855},
                            {-0.01031, 0.0607, -0.08954, 0.03084},
                            {-0.006283, 0.03836, -0.05853, 0.00891}}};

float mom_corr::FD_prot_Hmom_corr_upper_All_FD(float mom_, float dc_sec, float alpha_prot_mom_corr_FD) {
  if (dc_sec == 1) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtH[0][0][0] * pow(mom_, 3) + FDProtH[0][0][1] * pow(mom_, 2) +
                                            FDProtH[0][0][2] * mom_ + FDProtH[0][0][3]);
  } else if (dc_sec == 2) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtH[0][1][0] * pow(mom_, 3) + FDProtH[0][1][1] * pow(mom_, 2) +
                                            FDProtH[0][1][2] * mom_ + FDProtH[0][1][3]);
  } else if (dc_sec == 3) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtH[0][2][0] * pow(mom_, 3) + FDProtH[0][2][1] * pow(mom_, 2) +
                                            FDProtH[0][2][2] * mom_ + FDProtH[0][2][3]);
  } else if (dc_sec == 4) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtH[0][3][0] * pow(mom_, 3) + FDProtH[0][3][1] * pow(mom_, 2) +
                                            FDProtH[0][3][2] * mom_ + FDProtH[0][3][3]);
  } else if (dc_sec == 5) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtH[0][4][0] * pow(mom_, 3) + FDProtH[0][4][1] * pow(mom_, 2) +
                                            FDProtH[0][4][2] * mom_ + FDProtH[0][4][3]);
  } else if (dc_sec == 6) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtH[0][5][0] * pow(mom_, 3) + FDProtH[0][5][1] * pow(mom_, 2) +
                                            FDProtH[0][5][2] * mom_ + FDProtH[0][5][3]);
  } else
    return NAN;
}

float mom_corr::FD_prot_Hmom_corr_upper_Except_All_FD(float mom_, float dc_sec, float alpha_prot_mom_corr_FD) {
  if (dc_sec == 1) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtH[1][0][0] * pow(mom_, 3) + FDProtH[1][0][1] * pow(mom_, 2) +
                                            FDProtH[1][0][2] * mom_ + FDProtH[1][0][3]);
  } else if (dc_sec == 2) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtH[1][1][0] * pow(mom_, 3) + FDProtH[1][1][1] * pow(mom_, 2) +
                                            FDProtH[1][1][2] * mom_ + FDProtH[1][1][3]);
  } else if (dc_sec == 3) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtH[1][2][0] * pow(mom_, 3) + FDProtH[1][2][1] * pow(mom_, 2) +
                                            FDProtH[1][2][2] * mom_ + FDProtH[1][2][3]);
  } else if (dc_sec == 4) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtH[1][3][0] * pow(mom_, 3) + FDProtH[1][3][1] * pow(mom_, 2) +
                                            FDProtH[1][3][2] * mom_ + FDProtH[1][3][3]);
  } else if (dc_sec == 5) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtH[1][4][0] * pow(mom_, 3) + FDProtH[1][4][1] * pow(mom_, 2) +
                                            FDProtH[1][4][2] * mom_ + FDProtH[1][4][3]);
  } else if (dc_sec == 6) {
    return mom_ - alpha_prot_mom_corr_FD * (FDProtH[1][5][0] * pow(mom_, 3) + FDProtH[1][5][1] * pow(mom_, 2) +
                                            FDProtH[1][5][2] * mom_ + FDProtH[1][5][3]);
  } else
    return NAN;
}

// /// pip hadron corrections

// float alpha_pip_mom_corr_FD[4] = {0.1, 0.15, 0.5, 0.5};
// float alpha_pip_mom_corr_CD[3] = {0.8, 0.4, 0.8};

double CDPip[3][4] = {{0.06775, -0.1256, -0.03055, 0.002312},
                      {0.0484, -0.11993, 0.0746, -0.00975},
                      {-0.00775, 0.06445, -0.04684, 0.004112}};

float mom_corr::CD_pip_Hmom_corr(float mom_, float phi_, float alpha_pip_mom_corr_CD[3]) {
  if (phi_ > 270 || phi_ <= 30) {
    return mom_ - alpha_pip_mom_corr_CD[0] *
                      (CDPip[0][0] * pow(mom_, 3) + CDPip[0][1] * pow(mom_, 2) + CDPip[0][2] * mom_ + CDPip[0][3]);
  } else if (phi_ > 30 && phi_ <= 150) {
    return mom_ - alpha_pip_mom_corr_CD[1] *
                      (CDPip[1][0] * pow(mom_, 3) + CDPip[1][1] * pow(mom_, 2) + CDPip[1][2] * mom_ + CDPip[1][3]);
  } else if (phi_ > 150 && phi_ <= 270) {
    return mom_ - alpha_pip_mom_corr_CD[2] *
                      (CDPip[2][0] * pow(mom_, 3) + CDPip[2][1] * pow(mom_, 2) + CDPip[2][2] * mom_ + CDPip[2][3]);
  } else
    return NAN;
}
double FDPipL[2][6][4] = {{{-0.001051, 0.004627, 0.006058, -0.01855},
                           {-0.003084, 0.02007, -0.03488, 0.01718},
                           {-0.001668, 0.007435, -0.0001147, -0.005516},
                           {0.0003283, -0.005856, 0.02171, -0.012},
                           {-0.002243, 0.01291, -0.02052, 0.012505},
                           {-0.003408, 0.0175, -0.01814, 0.001455}},
                          {{0.002834, -0.0171, 0.03253, -0.02928},
                           {0.00416, -0.02376, 0.0383, -0.01701},
                           {0.00258, -0.01698, 0.0333, -0.01591},
                           {0.002327, -0.01192, 0.00987, 0.006634},
                           {0.001894, -0.006947, -0.003706, 0.00831},
                           {0.001051, -0.004684, 0.001716, -0.006424}}};

float mom_corr::FD_pip_Hmom_corr_lower_All_FD(float mom_, float dc_sec, float alpha_pip_mom_corr_FD) {
  if (dc_sec == 1) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipL[0][0][0] * pow(mom_, 3) + FDPipL[0][0][1] * pow(mom_, 2) +
                                           FDPipL[0][0][2] * mom_ + FDPipL[0][0][3]);
  } else if (dc_sec == 2) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipL[0][1][0] * pow(mom_, 3) + FDPipL[0][1][1] * pow(mom_, 2) +
                                           FDPipL[0][1][2] * mom_ + FDPipL[0][1][3]);
  } else if (dc_sec == 3) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipL[0][2][0] * pow(mom_, 3) + FDPipL[0][2][1] * pow(mom_, 2) +
                                           FDPipL[0][2][2] * mom_ + FDPipL[0][2][3]);
  } else if (dc_sec == 4) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipL[0][3][0] * pow(mom_, 3) + FDPipL[0][3][1] * pow(mom_, 2) +
                                           FDPipL[0][3][2] * mom_ + FDPipL[0][3][3]);
  } else if (dc_sec == 5) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipL[0][4][0] * pow(mom_, 3) + FDPipL[0][4][1] * pow(mom_, 2) +
                                           FDPipL[0][4][2] * mom_ + FDPipL[0][4][3]);
  } else if (dc_sec == 6) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipL[0][5][0] * pow(mom_, 3) + FDPipL[0][5][1] * pow(mom_, 2) +
                                           FDPipL[0][5][2] * mom_ + FDPipL[0][5][3]);
  } else
    return NAN;
}

float mom_corr::FD_pip_Hmom_corr_lower_Except_All_FD(float mom_, float dc_sec, float alpha_pip_mom_corr_FD) {
  if (dc_sec == 1) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipL[1][0][0] * pow(mom_, 3) + FDPipL[1][0][1] * pow(mom_, 2) +
                                           FDPipL[1][0][2] * mom_ + FDPipL[1][0][3]);
  } else if (dc_sec == 2) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipL[1][1][0] * pow(mom_, 3) + FDPipL[1][1][1] * pow(mom_, 2) +
                                           FDPipL[1][1][2] * mom_ + FDPipL[1][1][3]);
  } else if (dc_sec == 3) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipL[1][2][0] * pow(mom_, 3) + FDPipL[1][2][1] * pow(mom_, 2) +
                                           FDPipL[1][2][2] * mom_ + FDPipL[1][2][3]);
  } else if (dc_sec == 4) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipL[1][3][0] * pow(mom_, 3) + FDPipL[1][3][1] * pow(mom_, 2) +
                                           FDPipL[1][3][2] * mom_ + FDPipL[1][3][3]);
  } else if (dc_sec == 5) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipL[1][4][0] * pow(mom_, 3) + FDPipL[1][4][1] * pow(mom_, 2) +
                                           FDPipL[1][4][2] * mom_ + FDPipL[1][4][3]);
  } else if (dc_sec == 6) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipL[1][5][0] * pow(mom_, 3) + FDPipL[1][5][1] * pow(mom_, 2) +
                                           FDPipL[1][5][2] * mom_ + FDPipL[1][5][3]);
  } else
    return NAN;
}

double FDPipH[2][6][4] = {{{-0.00287, 0.00692, 0.01993, -0.02162},
                           {-0.0057, 0.01073, 0.0326, -0.04004},
                           {-0.001721, 0.00743, 0.01846, -0.02194},
                           {-0.006126, 0.00979, 0.03296, -0.03087},
                           {-0.001957, 0.00483, 0.01576, -0.004658},
                           {-0.001278, 0.004314, 0.01255, -0.00646}},
                          {{-0.002602, 0.00621, 0.0158, -0.03085},
                           {-0.00791, 0.00989, 0.0353, -0.04178},
                           {-0.003164, 0.006046, 0.01929, -0.01709},
                           {-0.001335, 0.002363, 0.00827, -0.004005},
                           {-0.003061, 0.003662, 0.01438, -0.01021},
                           {-0.001557, 0.002182, 0.006744, -0.01247}}};

float mom_corr::FD_pip_Hmom_corr_upper_All_FD(float mom_, float dc_sec, float alpha_pip_mom_corr_FD) {
  if (dc_sec == 1) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipH[0][0][0] * pow(mom_, 3) + FDPipH[0][0][1] * pow(mom_, 2) +
                                           FDPipH[0][0][2] * mom_ + FDPipH[0][0][3]);
  } else if (dc_sec == 2) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipH[0][1][0] * pow(mom_, 3) + FDPipH[0][1][1] * pow(mom_, 2) +
                                           FDPipH[0][1][2] * mom_ + FDPipH[0][1][3]);
  } else if (dc_sec == 3) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipH[0][2][0] * pow(mom_, 3) + FDPipH[0][2][1] * pow(mom_, 2) +
                                           FDPipH[0][2][2] * mom_ + FDPipH[0][2][3]);
  } else if (dc_sec == 4) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipH[0][3][0] * pow(mom_, 3) + FDPipH[0][3][1] * pow(mom_, 2) +
                                           FDPipH[0][3][2] * mom_ + FDPipH[0][3][3]);
  } else if (dc_sec == 5) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipH[0][4][0] * pow(mom_, 3) + FDPipH[0][4][1] * pow(mom_, 2) +
                                           FDPipH[0][4][2] * mom_ + FDPipH[0][4][3]);
  } else if (dc_sec == 6) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipH[0][5][0] * pow(mom_, 3) + FDPipH[0][5][1] * pow(mom_, 2) +
                                           FDPipH[0][5][2] * mom_ + FDPipH[0][5][3]);
  } else
    return NAN;
}

float mom_corr::FD_pip_Hmom_corr_upper_Except_All_FD(float mom_, float dc_sec, float alpha_pip_mom_corr_FD) {
  if (dc_sec == 1) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipH[1][0][0] * pow(mom_, 3) + FDPipH[1][0][1] * pow(mom_, 2) +
                                           FDPipH[1][0][2] * mom_ + FDPipH[1][0][3]);
  } else if (dc_sec == 2) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipH[1][1][0] * pow(mom_, 3) + FDPipH[1][1][1] * pow(mom_, 2) +
                                           FDPipH[1][1][2] * mom_ + FDPipH[1][1][3]);
  } else if (dc_sec == 3) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipH[1][2][0] * pow(mom_, 3) + FDPipH[1][2][1] * pow(mom_, 2) +
                                           FDPipH[1][2][2] * mom_ + FDPipH[1][2][3]);
  } else if (dc_sec == 4) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipH[1][3][0] * pow(mom_, 3) + FDPipH[1][3][1] * pow(mom_, 2) +
                                           FDPipH[1][3][2] * mom_ + FDPipH[1][3][3]);
  } else if (dc_sec == 5) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipH[1][4][0] * pow(mom_, 3) + FDPipH[1][4][1] * pow(mom_, 2) +
                                           FDPipH[1][4][2] * mom_ + FDPipH[1][4][3]);
  } else if (dc_sec == 6) {
    return mom_ - alpha_pip_mom_corr_FD * (FDPipH[1][5][0] * pow(mom_, 3) + FDPipH[1][5][1] * pow(mom_, 2) +
                                           FDPipH[1][5][2] * mom_ + FDPipH[1][5][3]);
  } else
    return NAN;
}

// /// pim hadron corrections
// float alpha_pim_mom_corr_FD[4] = {0.5, 0.15, 0.3, 0.3};
// float alpha_pim_mom_corr_CD[3] = {0.5, 1.0, 0.5};

// float alpha_pim_mom_corr_FD[4] = {0.0, 0.0, 0.0, 0.0};
// float alpha_pim_mom_corr_CD[3] = {0., 0.0, 0.0};

// double CDPim[3][4] = {
//     {0.0531, -0.0899, 0.05328, -0.01124}, {0.02277, -0.02846, 0.04657, -0.005}, {0.05997, -0.1099, 0.0093,
//     0.006657}}; // becareful plot has 1st sector in 3rd place
double CDPim[3][5] = {{-0.06088, 0.2715, -0.355, 0.1799, -0.03076},
                      {-0.01833, 0.08844, -0.1082, 0.08466, -0.01088},
                      {-0.1163, 0.4768, -0.616, 0.2512, -0.03062}};
float mom_corr::CD_pim_Hmom_corr(float mom_, float phi_, float alpha_pim_mom_corr_CD[3]) {
  if (phi_ > 270 || phi_ <= 30) {
    return mom_ - alpha_pim_mom_corr_CD[0] * (CDPim[0][0] * pow(mom_, 4) + CDPim[0][1] * pow(mom_, 3) +
                                              CDPim[0][2] * pow(mom_, 2) + CDPim[0][3] * mom_ + CDPim[0][4]);
  } else if (phi_ > 30 && phi_ <= 150) {
    return mom_ - alpha_pim_mom_corr_CD[1] * (CDPim[1][0] * pow(mom_, 4) + CDPim[1][1] * pow(mom_, 3) +
                                              CDPim[1][2] * pow(mom_, 2) + CDPim[1][3] * mom_ + CDPim[1][4]);
  } else if (phi_ > 150 && phi_ <= 270) {
    return mom_ - alpha_pim_mom_corr_CD[2] * (CDPim[2][0] * pow(mom_, 4) + CDPim[2][1] * pow(mom_, 3) +
                                              CDPim[2][2] * pow(mom_, 2) + CDPim[2][3] * mom_ + CDPim[2][4]);
  } else
    return NAN;
}

double FDPimL[2][6][4] = {{{-0.003864, 0.0297, -0.0769, 0.06757},
                           {-0.00433, 0.03613, -0.1026, 0.0982},
                           {0.00489, -0.02971, 0.04425, 0.002811},
                           {-2.235e-05, 0.004494, -0.03032, 0.04724},
                           {0.00835, -0.0581, 0.1112, -0.0421},
                           {0.003828, -0.01945, 0.01244, 0.01949}},
                          {{-0.00321, 0.0282, -0.08203, 0.06024},
                           {0.003355, -0.02539, 0.05545, -0.03723},
                           {0.00402, -0.02834, 0.06, -0.0354},
                           {0.000703, -0.003294, 0.00447, -0.009796},
                           {-0.002449, 0.0176, -0.03662, 0.00562},
                           {-0.003347, 0.0308, -0.0873, 0.04517}}};

float mom_corr::FD_pim_Hmom_corr_lower_All_FD(float mom_, float dc_sec, float alpha_pim_mom_corr_FD) {
  if (dc_sec == 1) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimL[0][0][0] * pow(mom_, 3) + FDPimL[0][0][1] * pow(mom_, 2) +
                                           FDPimL[0][0][2] * mom_ + FDPimL[0][0][3]);
  } else if (dc_sec == 2) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimL[0][1][0] * pow(mom_, 3) + FDPimL[0][1][1] * pow(mom_, 2) +
                                           FDPimL[0][1][2] * mom_ + FDPimL[0][1][3]);
  } else if (dc_sec == 3) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimL[0][2][0] * pow(mom_, 3) + FDPimL[0][2][1] * pow(mom_, 2) +
                                           FDPimL[0][2][2] * mom_ + FDPimL[0][2][3]);
  } else if (dc_sec == 4) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimL[0][3][0] * pow(mom_, 3) + FDPimL[0][3][1] * pow(mom_, 2) +
                                           FDPimL[0][3][2] * mom_ + FDPimL[0][3][3]);
  } else if (dc_sec == 5) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimL[0][4][0] * pow(mom_, 3) + FDPimL[0][4][1] * pow(mom_, 2) +
                                           FDPimL[0][4][2] * mom_ + FDPimL[0][4][3]);
  } else if (dc_sec == 6) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimL[0][5][0] * pow(mom_, 3) + FDPimL[0][5][1] * pow(mom_, 2) +
                                           FDPimL[0][5][2] * mom_ + FDPimL[0][5][3]);
  } else
    return NAN;
}

float mom_corr::FD_pim_Hmom_corr_lower_Except_All_FD(float mom_, float dc_sec, float alpha_pim_mom_corr_FD) {
  if (dc_sec == 1) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimL[1][0][0] * pow(mom_, 3) + FDPimL[1][0][1] * pow(mom_, 2) +
                                           FDPimL[1][0][2] * mom_ + FDPimL[1][0][3]);
  } else if (dc_sec == 2) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimL[1][1][0] * pow(mom_, 3) + FDPimL[1][1][1] * pow(mom_, 2) +
                                           FDPimL[1][1][2] * mom_ + FDPimL[1][1][3]);
  } else if (dc_sec == 3) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimL[1][2][0] * pow(mom_, 3) + FDPimL[1][2][1] * pow(mom_, 2) +
                                           FDPimL[1][2][2] * mom_ + FDPimL[1][2][3]);
  } else if (dc_sec == 4) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimL[1][3][0] * pow(mom_, 3) + FDPimL[1][3][1] * pow(mom_, 2) +
                                           FDPimL[1][3][2] * mom_ + FDPimL[1][3][3]);
  } else if (dc_sec == 5) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimL[1][4][0] * pow(mom_, 3) + FDPimL[1][4][1] * pow(mom_, 2) +
                                           FDPimL[1][4][2] * mom_ + FDPimL[1][4][3]);
  } else if (dc_sec == 6) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimL[1][5][0] * pow(mom_, 3) + FDPimL[1][5][1] * pow(mom_, 2) +
                                           FDPimL[1][5][2] * mom_ + FDPimL[1][5][3]);
  } else
    return NAN;
}

double FDPimH[2][6][4] = {{{-0.0002575, -0.002295, -0.003235, 0.01292},
                           {-0.00398, 0.0052, 0.02222, -0.02208},
                           {-0.003052, 0.003235, 0.01591, -0.01168},
                           {0.001004, -0.003637, -0.00909, 0.02184},
                           {-0.000769, -0.0006876, 0.0008802, 0.002695},
                           {0.002268, -0.003813, -0.013794, 0.0193}},
                          {{0.001858, -0.003899, -0.01464, 0.009155},
                           {0.0002279, -0.000515, -0.002798, -0.003864},
                           {-0.001788, 0.00237, 0.008934, -0.01624},
                           {-0.002106, 0.002811, 0.0091, -0.0273},
                           {-0.003132, 0.00478, 0.01424, -0.04587},
                           {0.001281, -0.001394, -0.01089, -0.0178}}};

float mom_corr::FD_pim_Hmom_corr_upper_All_FD(float mom_, float dc_sec, float alpha_pim_mom_corr_FD) {
  if (dc_sec == 1) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimH[0][0][0] * pow(mom_, 3) + FDPimH[0][0][1] * pow(mom_, 2) +
                                           FDPimH[0][0][2] * mom_ + FDPimH[0][0][3]);
  } else if (dc_sec == 2) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimH[0][1][0] * pow(mom_, 3) + FDPimH[0][1][1] * pow(mom_, 2) +
                                           FDPimH[0][1][2] * mom_ + FDPimH[0][1][3]);
  } else if (dc_sec == 3) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimH[0][2][0] * pow(mom_, 3) + FDPimH[0][2][1] * pow(mom_, 2) +
                                           FDPimH[0][2][2] * mom_ + FDPimH[0][2][3]);
  } else if (dc_sec == 4) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimH[0][3][0] * pow(mom_, 3) + FDPimH[0][3][1] * pow(mom_, 2) +
                                           FDPimH[0][3][2] * mom_ + FDPimH[0][3][3]);
  } else if (dc_sec == 5) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimH[0][4][0] * pow(mom_, 3) + FDPimH[0][4][1] * pow(mom_, 2) +
                                           FDPimH[0][4][2] * mom_ + FDPimH[0][4][3]);
  } else if (dc_sec == 6) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimH[0][5][0] * pow(mom_, 3) + FDPimH[0][5][1] * pow(mom_, 2) +
                                           FDPimH[0][5][2] * mom_ + FDPimH[0][5][3]);
  } else
    return NAN;
}

float mom_corr::FD_pim_Hmom_corr_upper_Except_All_FD(float mom_, float dc_sec, float alpha_pim_mom_corr_FD) {
  if (dc_sec == 1) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimH[1][0][0] * pow(mom_, 3) + FDPimH[1][0][1] * pow(mom_, 2) +
                                           FDPimH[1][0][2] * mom_ + FDPimH[1][0][3]);
  } else if (dc_sec == 2) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimH[1][1][0] * pow(mom_, 3) + FDPimH[1][1][1] * pow(mom_, 2) +
                                           FDPimH[1][1][2] * mom_ + FDPimH[1][1][3]);
  } else if (dc_sec == 3) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimH[1][2][0] * pow(mom_, 3) + FDPimH[1][2][1] * pow(mom_, 2) +
                                           FDPimH[1][2][2] * mom_ + FDPimH[1][2][3]);
  } else if (dc_sec == 4) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimH[1][3][0] * pow(mom_, 3) + FDPimH[1][3][1] * pow(mom_, 2) +
                                           FDPimH[1][3][2] * mom_ + FDPimH[1][3][3]);
  } else if (dc_sec == 5) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimH[1][4][0] * pow(mom_, 3) + FDPimH[1][4][1] * pow(mom_, 2) +
                                           FDPimH[1][4][2] * mom_ + FDPimH[1][4][3]);
  } else if (dc_sec == 6) {
    return mom_ - alpha_pim_mom_corr_FD * (FDPimH[1][5][0] * pow(mom_, 3) + FDPimH[1][5][1] * pow(mom_, 2) +
                                           FDPimH[1][5][2] * mom_ + FDPimH[1][5][3]);
  } else
    return NAN;
}
////////// our final mom corr finished ////////////////

//////////////////// new mom correction start pass2 2024
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

// Below shows how the corrections are to be applied using the ROOT momentum 4-vector using the above code:
// auto fe = dppC(ex, ey, ez, esec, 0) + 1;
// auto fpip = dppC(pipx, pipy, pipz, pipsec, 1) + 1;
// auto fpim = dppC(pimx, pimy, pimz, pimsec, 2) + 1;
// auto fpro = dppC(prox, proy, proz, prosec, 3) + 1;

// auto eleC = ROOT::Math::PxPyPzMVector(ex * fe, ey* fe, ez* fe, 0);
// auto pipC = ROOT::Math::PxPyPzMVector(pipx * fpip, pipy* fpip, pipz* fpip, 0.13957);
// auto pimC = ROOT::Math::PxPyPzMVector(pimx * fpim, pimy* fpim, pimz* fpim, 0.13957);
// auto proC = ROOT::Math::PxPyPzMVector(prox * fpro, proy* fpro, proz* fpro, 0.938);
////////////////// Eloss corr pip ////////////////////////////////////

// double eloss_pip(double pim_p, double pim_theta, double status_pim) {
//   double dp_pim = 0.0;

//   // INBENDING
//   if (is_FD(status_pion)) {  // Forward Detector
//     if (pim_theta < 27) {
//       dp_pim = 0.00044836 * pim_p + 0.00325965;
//     } else if (pim_theta >= 27) {
//       dp_pim = -0.00208368 * pim_p + 0.00908514;
//     }
//   }
// }
////////////////// Eloss corr pip ////////////////////////////////////

double mom_corr::elossPipFD(double pion_p, double pip_theta) {
  // momentum loss correction for low momentum pions:
  // input: p = pion momentum in GeV, pip_theta = pion theta in degree,
  //        pion_det = pion detector (2 = FD, 3 = CD),  outbending = torus polarity
  // output: dp_pion_fd = generated momentum - reconstructed momentum = momentum loss (+) / gain (-)

  double dp_pion_fd = 0.0;

  // INBENDING
  if (pip_theta < 27) {
    dp_pion_fd = 0.00342646 + (-0.00282934) * pion_p + (0.00205983) * pow(pion_p, 2) + (-0.00043158) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta < 27 && pion_p >= 2.5) {
    dp_pion_fd =
        0.00342646 + (-0.00282934) * 2.5 + (0.00205983) * pow(2.5, 2) + (-0.00043158) * pow(2.5, 3) + (0) * pow(2.5, 4);
  }
  if (pip_theta > 27 && pip_theta < 28) {
    dp_pion_fd = 0.00328565 + (-0.00376042) * pion_p + (0.00433886) * pow(pion_p, 2) + (-0.00141614) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 27 && pip_theta < 28 && pion_p >= 1.83) {
    dp_pion_fd = 0.00328565 + (-0.00376042) * 1.83 + (0.00433886) * pow(1.83, 2) + (-0.00141614) * pow(1.83, 3) +
                 (0) * pow(1.83, 4);
  }
  if (pip_theta > 28 && pip_theta < 29) {
    dp_pion_fd = 0.00328579 + (-0.00281121) * pion_p + (0.00342749) * pow(pion_p, 2) + (-0.000932614) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 28 && pip_theta < 29 && pion_p >= 2) {
    dp_pion_fd =
        0.00328579 + (-0.00281121) * 2 + (0.00342749) * pow(2, 2) + (-0.000932614) * pow(2, 3) + (0) * pow(2, 4);
  }
  if (pip_theta > 29 && pip_theta < 30) {
    dp_pion_fd = 0.00167358 + (0.00441871) * pion_p + (-0.000834667) * pow(pion_p, 2) +
                 (-0.000137968) * pow(pion_p, 3) + (0) * pow(pion_p, 4);
  }
  if (pip_theta > 29 && pip_theta < 30 && pion_p >= 1.9) {
    dp_pion_fd = 0.00167358 + (0.00441871) * 1.9 + (-0.000834667) * pow(1.9, 2) + (-0.000137968) * pow(1.9, 3) +
                 (0) * pow(1.9, 4);
  }
  if (pip_theta > 30 && pip_theta < 31) {
    dp_pion_fd = 0.00274159 + (0.00635686) * pion_p + (-0.00380977) * pow(pion_p, 2) + (0.00071627) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 30 && pip_theta < 31 && pion_p >= 1.9) {
    dp_pion_fd =
        0.00274159 + (0.00635686) * 1.9 + (-0.00380977) * pow(1.9, 2) + (0.00071627) * pow(1.9, 3) + (0) * pow(1.9, 4);
  }
  if (pip_theta > 31 && pip_theta < 32) {
    dp_pion_fd = 0.00450241 + (0.00248969) * pion_p + (-0.00336795) * pow(pion_p, 2) + (0.00111193) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 31 && pip_theta < 32 && pion_p >= 1.8) {
    dp_pion_fd =
        0.00450241 + (0.00248969) * 1.8 + (-0.00336795) * pow(1.8, 2) + (0.00111193) * pow(1.8, 3) + (0) * pow(1.8, 4);
  }
  if (pip_theta > 32 && pip_theta < 33) {
    dp_pion_fd = 0.00505593 + (-0.00246203) * pion_p + (0.00172984) * pow(pion_p, 2) + (-0.000406701) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 32 && pip_theta < 33 && pion_p >= 1.8) {
    dp_pion_fd = 0.00505593 + (-0.00246203) * 1.8 + (0.00172984) * pow(1.8, 2) + (-0.000406701) * pow(1.8, 3) +
                 (0) * pow(1.8, 4);
  }
  if (pip_theta > 33 && pip_theta < 34) {
    dp_pion_fd = 0.00273402 + (0.00440449) * pion_p + (-0.00373488) * pow(pion_p, 2) + (0.000996612) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 33 && pip_theta < 34 && pion_p >= 1.8) {
    dp_pion_fd =
        0.00273402 + (0.00440449) * 1.8 + (-0.00373488) * pow(1.8, 2) + (0.000996612) * pow(1.8, 3) + (0) * pow(1.8, 4);
  }
  if (pip_theta > 34 && pip_theta < 35) {
    dp_pion_fd = 0.00333542 + (0.00439874) * pion_p + (-0.00397776) * pow(pion_p, 2) + (0.00105586) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 34 && pip_theta < 35 && pion_p >= 1.8) {
    dp_pion_fd =
        0.00333542 + (0.00439874) * 1.8 + (-0.00397776) * pow(1.8, 2) + (0.00105586) * pow(1.8, 3) + (0) * pow(1.8, 4);
  }
  if (pip_theta > 35 && pip_theta < 36) {
    dp_pion_fd = 0.00354663 + (0.00565397) * pion_p + (-0.00513503) * pow(pion_p, 2) + (0.00153346) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 35 && pip_theta < 36 && pion_p >= 1.8) {
    dp_pion_fd =
        0.00354663 + (0.00565397) * 1.8 + (-0.00513503) * pow(1.8, 2) + (0.00153346) * pow(1.8, 3) + (0) * pow(1.8, 4);
  }
  if (pip_theta > 36 && pip_theta < 37) {
    dp_pion_fd = 0.00333909 + (0.00842367) * pion_p + (-0.0077321) * pow(pion_p, 2) + (0.0022489) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 36 && pip_theta < 37 && pion_p >= 1.8) {
    dp_pion_fd =
        0.00333909 + (0.00842367) * 1.8 + (-0.0077321) * pow(1.8, 2) + (0.0022489) * pow(1.8, 3) + (0) * pow(1.8, 4);
  }
  if (pip_theta > 37 && pip_theta < 38) {
    dp_pion_fd = 0.00358828 + (0.0112108) * pion_p + (-0.0133854) * pow(pion_p, 2) + (0.00486924) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 37 && pip_theta < 38 && pion_p >= 1.4) {
    dp_pion_fd =
        0.00358828 + (0.0112108) * 1.4 + (-0.0133854) * pow(1.4, 2) + (0.00486924) * pow(1.4, 3) + (0) * pow(1.4, 4);
  }
  if (pip_theta > 38 && pip_theta < 39) {
    dp_pion_fd = 0.00354343 + (0.0117121) * pion_p + (-0.0129649) * pow(pion_p, 2) + (0.00455602) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 38 && pip_theta < 39 && pion_p >= 1.3) {
    dp_pion_fd =
        0.00354343 + (0.0117121) * 1.3 + (-0.0129649) * pow(1.3, 2) + (0.00455602) * pow(1.3, 3) + (0) * pow(1.3, 4);
  }
  if (pip_theta > 39 && pip_theta < 40) {
    dp_pion_fd = -0.00194951 + (0.0409713) * pion_p + (-0.0595861) * pow(pion_p, 2) + (0.0281588) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 39 && pip_theta < 40 && pion_p >= 0.9) {
    dp_pion_fd =
        -0.00194951 + (0.0409713) * 0.9 + (-0.0595861) * pow(0.9, 2) + (0.0281588) * pow(0.9, 3) + (0) * pow(0.9, 4);
  }
  if (pip_theta > 40 && pip_theta < 41) {
    dp_pion_fd = -0.0099217 + (0.0808096) * pion_p + (-0.119836) * pow(pion_p, 2) + (0.0559553) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 40 && pip_theta < 41 && pion_p >= 0.75) {
    dp_pion_fd =
        -0.0099217 + (0.0808096) * 0.75 + (-0.119836) * pow(0.75, 2) + (0.0559553) * pow(0.75, 3) + (0) * pow(0.75, 4);
  }
  if (pip_theta > 41 && pip_theta < 42) {
    dp_pion_fd = 0.00854898 + (0.00025037) * pion_p + (-0.0113992) * pow(pion_p, 2) + (0.0145178) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 41 && pip_theta < 42 && pion_p >= 0.65) {
    dp_pion_fd = 0.00854898 + (0.00025037) * 0.65 + (-0.0113992) * pow(0.65, 2) + (0.0145178) * pow(0.65, 3) +
                 (0) * pow(0.65, 4);
  }
  if (pip_theta > 42) {
    dp_pion_fd = 0.00564818 + (0.00706606) * pion_p + (0.0042602) * pow(pion_p, 2) + (-0.01141) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 42 && pion_p >= 0.65) {
    dp_pion_fd =
        0.00564818 + (0.00706606) * 0.65 + (0.0042602) * pow(0.65, 2) + (-0.01141) * pow(0.65, 3) + (0) * pow(0.65, 4);
  }
  return dp_pion_fd;
}
double mom_corr::elossPipCD(double pion_p, double pip_theta) {
  double dp_pion_cd = 0.0;

  if (pip_theta < 39) {
    dp_pion_cd =
        -0.045 + (-0.102652) + (0.455589) * pion_p + (-0.671635) * pow(pion_p, 2) + (0.303814) * pow(pion_p, 3);
  }
  if (pip_theta < 39 && pion_p >= 0.7) {
    dp_pion_cd = -0.045 + (-0.102652) + (0.455589) * 0.7 + (-0.671635) * pow(0.7, 2) + (0.303814) * pow(0.7, 3);
  }
  if (pip_theta > 39 && pip_theta < 40) {
    dp_pion_cd = 0.0684552 + (-0.766492) * pion_p + (1.73092) * pow(pion_p, 2) + (-1.46215) * pow(pion_p, 3) +
                 (0.420127) * pow(pion_p, 4);
  }
  if (pip_theta > 39 && pip_theta < 40 && pion_p >= 1.4) {
    dp_pion_cd =
        0.0684552 + (-0.766492) * 1.4 + (1.73092) * pow(1.4, 2) + (-1.46215) * pow(1.4, 3) + (0.420127) * pow(1.4, 4);
  }
  if (pip_theta > 40 && pip_theta < 41) {
    dp_pion_cd = 0.751549 + (-7.4593) * pion_p + (26.8037) * pow(pion_p, 2) + (-47.1576) * pow(pion_p, 3) +
                 (43.8527) * pow(pion_p, 4) + (-20.7039) * pow(pion_p, 5) + (3.90931) * pow(pion_p, 6);
  }
  if (pip_theta > 40 && pip_theta < 41 && pion_p >= 1.45) {
    dp_pion_cd = 0.751549 + (-7.4593) * 1.45 + (26.8037) * pow(1.45, 2) + (-47.1576) * pow(1.45, 3) +
                 (43.8527) * pow(1.45, 4) + (-20.7039) * pow(1.45, 5) + (3.90931) * pow(1.45, 6);
  }
  if (pip_theta > 41 && pip_theta < 42) {
    dp_pion_cd = -1.35043 + (10.0788) * pion_p + (-30.4829) * pow(pion_p, 2) + (47.7792) * pow(pion_p, 3) +
                 (-40.996) * pow(pion_p, 4) + (18.2662) * pow(pion_p, 5) + (-3.30449) * pow(pion_p, 6);
  }
  if (pip_theta > 41 && pip_theta < 42 && pion_p >= 1.2) {
    dp_pion_cd = -1.35043 + (10.0788) * 1.2 + (-30.4829) * pow(1.2, 2) + (47.7792) * pow(1.2, 3) +
                 (-40.996) * pow(1.2, 4) + (18.2662) * pow(1.2, 5) + (-3.30449) * pow(1.2, 6);
  }
  if (pip_theta > 42 && pip_theta < 43) {
    dp_pion_cd = -0.0231195 + (0.0744589) * pion_p + (-0.0807029) * pow(pion_p, 2) + (0.0264266) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 42 && pip_theta < 43 && pion_p >= 1.3) {
    dp_pion_cd =
        -0.0231195 + (0.0744589) * 1.3 + (-0.0807029) * pow(1.3, 2) + (0.0264266) * pow(1.3, 3) + (0) * pow(1.3, 4);
  }
  if (pip_theta > 43 && pip_theta < 44) {
    dp_pion_cd = -0.00979928 + (0.0351043) * pion_p + (-0.0365865) * pow(pion_p, 2) + (0.00977218) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 43 && pip_theta < 44 && pion_p >= 1.1) {
    dp_pion_cd =
        -0.00979928 + (0.0351043) * 1.1 + (-0.0365865) * pow(1.1, 2) + (0.00977218) * pow(1.1, 3) + (0) * pow(1.1, 4);
  }
  if (pip_theta > 44 && pip_theta < 45) {
    dp_pion_cd = 0.00108491 + (-0.00924885) * pion_p + (0.0216431) * pow(pion_p, 2) + (-0.0137762) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 44 && pip_theta < 45 && pion_p >= 1.1) {
    dp_pion_cd =
        0.00108491 + (-0.00924885) * 1.1 + (0.0216431) * pow(1.1, 2) + (-0.0137762) * pow(1.1, 3) + (0) * pow(1.1, 4);
  }
  if (pip_theta > 45 && pip_theta < 55) {
    dp_pion_cd = 0.0092263 + (-0.0676178) * pion_p + (0.168778) * pow(pion_p, 2) + (-0.167463) * pow(pion_p, 3) +
                 (0.05661) * pow(pion_p, 4);
  }
  if (pip_theta > 45 && pip_theta < 55 && pion_p >= 1.3) {
    dp_pion_cd =
        0.0092263 + (-0.0676178) * 1.3 + (0.168778) * pow(1.3, 2) + (-0.167463) * pow(1.3, 3) + (0.05661) * pow(1.3, 4);
  }
  if (pip_theta > 55 && pip_theta < 65) {
    dp_pion_cd = 0.00805642 + (-0.0670962) * pion_p + (0.188536) * pow(pion_p, 2) + (-0.20571) * pow(pion_p, 3) +
                 (0.0765) * pow(pion_p, 4);
  }
  if (pip_theta > 55 && pip_theta < 65 && pion_p >= 1.05) {
    dp_pion_cd = 0.00805642 + (-0.0670962) * 1.05 + (0.188536) * pow(1.05, 2) + (-0.20571) * pow(1.05, 3) +
                 (0.0765) * pow(1.05, 4);
  }
  if (pip_theta > 65 && pip_theta < 75) {
    dp_pion_cd = 0.00312202 + (-0.0269717) * pion_p + (0.0715236) * pow(pion_p, 2) + (-0.0545622) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 65 && pip_theta < 75 && pion_p >= 0.75) {
    dp_pion_cd = 0.00312202 + (-0.0269717) * 0.75 + (0.0715236) * pow(0.75, 2) + (-0.0545622) * pow(0.75, 3) +
                 (0) * pow(0.75, 4);
  }
  if (pip_theta > 75 && pip_theta < 85) {
    dp_pion_cd = 0.00424971 + (-0.0367683) * pion_p + (0.10417) * pow(pion_p, 2) + (-0.0899651) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 75 && pip_theta < 85 && pion_p >= 0.65) {
    dp_pion_cd =
        0.00424971 + (-0.0367683) * 0.65 + (0.10417) * pow(0.65, 2) + (-0.0899651) * pow(0.65, 3) + (0) * pow(0.65, 4);
  }
  if (pip_theta > 85 && pip_theta < 95) {
    dp_pion_cd = 0.00654123 + (-0.0517915) * pion_p + (0.147888) * pow(pion_p, 2) + (-0.14253) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 85 && pip_theta < 95 && pion_p >= 0.5) {
    dp_pion_cd =
        0.00654123 + (-0.0517915) * 0.5 + (0.147888) * pow(0.5, 2) + (-0.14253) * pow(0.5, 3) + (0) * pow(0.5, 4);
  }
  if (pip_theta > 95 && pip_theta < 105) {
    dp_pion_cd = -0.00111721 + (0.00478119) * pion_p + (0.0158753) * pow(pion_p, 2) + (-0.052902) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 95 && pip_theta < 105 && pion_p >= 0.45) {
    dp_pion_cd = -0.00111721 + (0.00478119) * 0.45 + (0.0158753) * pow(0.45, 2) + (-0.052902) * pow(0.45, 3) +
                 (0) * pow(0.45, 4);
  }
  if (pip_theta > 105 && pip_theta < 115) {
    dp_pion_cd = -0.00239839 + (0.00790738) * pion_p + (0.0311713) * pow(pion_p, 2) + (-0.104157) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 105 && pip_theta < 115 && pion_p >= 0.35) {
    dp_pion_cd = -0.00239839 + (0.00790738) * 0.35 + (0.0311713) * pow(0.35, 2) + (-0.104157) * pow(0.35, 3) +
                 (0) * pow(0.35, 4);
  }
  if (pip_theta > 115 && pip_theta < 125) {
    dp_pion_cd = -0.00778793 + (0.0256774) * pion_p + (0.0932503) * pow(pion_p, 2) + (-0.32771) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 115 && pip_theta < 125 && pion_p >= 0.35) {
    dp_pion_cd =
        -0.00778793 + (0.0256774) * 0.35 + (0.0932503) * pow(0.35, 2) + (-0.32771) * pow(0.35, 3) + (0) * pow(0.35, 4);
  }
  if (pip_theta > 125 && pip_theta < 135) {
    dp_pion_cd = -0.00292778 + (-0.00536697) * pion_p + (-0.00414351) * pow(pion_p, 2) + (0.0196431) * pow(pion_p, 3) +
                 (0) * pow(pion_p, 4);
  }
  if (pip_theta > 125 && pip_theta < 135 && pion_p >= 0.35) {
    dp_pion_cd = -0.00292778 + (-0.00536697) * 0.35 + (-0.00414351) * pow(0.35, 2) + (0.0196431) * pow(0.35, 3) +
                 (0) * pow(0.35, 4);
  }

  return dp_pion_cd;
}
