
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
bool is_lower_band(float mom_, float theta_DCr1_, int status_) {
  // if (dc_sec >= 1 && dc_sec <= 6) {
  if (status_ > 2000 && status_ <= 4000) {
    if (theta_DCr1_ < -53.14680163254601 + 79.61307254040804 * pow((mom_ - 0.3), 0.05739232362022314)) {
      return true;
    } else
      return false;
  } else
    return false;
}

float CD_prot_Emom_corr(float mom_, float theta_) {
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

float FD_prot_Emom_corr_lower(float mom_, float theta_) {
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
float FD_prot_Emom_corr_upper(float mom_, float theta_) {
  return mom_ +
         ((-6.33926614e-05) * pow(mom_, 3) + 3.21255513e-04 * pow(mom_, 2) + (-4.80918164e-04) * mom_ +
          1.94036549e-04) *
             pow(theta_, 2) +
         (0.00385508 * pow(mom_, 3) + (-0.0193179) * pow(mom_, 2) + 0.0279666 * mom_ + (-0.01032478)) * pow(theta_, 1) +
         (-0.06010495) * pow(mom_, 3) + 0.30123952 * pow(mom_, 2) + (-0.43371747) * mom_ + 0.16664826;
}

float CD_prot_Eth_corr(float mom_, float theta_) {
  return theta_ +
         (0.01794123 * pow(mom_, 3) + (-0.09198341) * pow(mom_, 2) + 0.15148531 * mom_ + (-0.0941657)) *
             pow(theta_, 1) +
         (-0.7392232) * pow(mom_, 3) + 3.93194154 * pow(mom_, 2) + (-6.83838677) * mom_ + 4.5505975;
}

float FD_prot_Eth_corr_lower(float mom_, float theta_) {
  return theta_ +
         (2.14391671e-05 * pow(mom_, 3) + (-1.69415274e-04) * pow(mom_, 2) + 3.62193361e-04 * mom_ +
          (-1.72672065e-04)) *
             pow(theta_, 2) +
         ((-0.00014124) * pow(mom_, 3) + 0.00017366 * pow(mom_, 2) + 0.00466645 * mom_ + (-0.0111939)) *
             pow(theta_, 1) +
         (-0.00031486) * pow(mom_, 3) + 0.00897261 * pow(mom_, 2) + (-0.05371869) * mom_ + 0.08065691;
}

float FD_prot_Eth_corr_upper(float mom_, float theta_) {
  return theta_ +
         (0.00165645 * pow(mom_, 3) + (-0.00983809) * pow(mom_, 2) + 0.01821203 * mom_ + (-0.01069836)) *
             pow(theta_, 2) +
         ((-0.10409645) * pow(mom_, 3) + 0.61354318 * pow(mom_, 2) + (-1.12258434) * mom_ + 0.64393271) *
             pow(theta_, 1) +
         1.66090372 * pow(mom_, 3) + (-9.75714605) * pow(mom_, 2) + 17.77247321 * mom_ + (-10.0865238);
}

float CD_prot_Eph_corr(float mom_, float theta_, float phi_) {
  return phi_ +
         (0.0152672 * pow(mom_, 3) + (-0.07306141) * pow(mom_, 2) + 0.09932124 * mom_ + (-0.04428166)) *
             pow(theta_, 1) +
         (-0.71565591) * pow(mom_, 3) + 3.37273717 * pow(mom_, 2) + (-4.54191832) * mom_ + 1.87540743;
}
float FD_prot_Eph_corr_lower(float mom_, float theta_, float phi_) {
  return phi_ +
         ((-4.86422409e-05) * pow(mom_, 4) + 1.21216530e-03 * pow(mom_, 3) + (-8.15266042e-03) * pow(mom_, 2) +
          1.93258907e-02 * mom_ + (-1.28009681e-02)) *
             pow(theta_, 1) +
         0.01081378 * pow(mom_, 4) + (-0.14401558) * pow(mom_, 3) + 0.69173611 * pow(mom_, 2) + (-1.3964496) * mom_ +
         0.95058901;
}

float FD_prot_Eph_corr_upper(float mom_, float theta_, float phi_) {
  return phi_ +
         ((-0.01255713) * pow(mom_, 3) + 0.07022673 * pow(mom_, 2) + (-0.12047137) * mom_ + 0.06254443) *
             pow(theta_, 1) +
         0.27588214 * pow(mom_, 3) + (-1.37114604) * pow(mom_, 2) + 1.82000373 * mom_ + (-0.40190107);
}
// // energy loss corrections parameters for momentum of proton
float A_p(float mom_, float theta_, float theta_DCr1_p, int dc_sec) {
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

float B_p(float mom_, float theta_, float theta_DCr1_p, int dc_sec) {
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

float CD_pip_Emom_corr(float mom_, float theta_) {
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
float FD_pip_Emom_corr_lower(float mom_, float theta_) {
  return mom_ + (-4.67842670e-05) * pow(mom_, 3) + 3.37133020e-04 * pow(mom_, 2) + (-4.79135831e-04) * mom_ +
         2.70872474e-03;
}
float FD_pip_Emom_corr_upper(float mom_, float theta_) {
  return mom_ + (-0.00125149) * pow(mom_, 3) + 0.0053441 * pow(mom_, 2) + (-0.00765213) * mom_ + 0.0102172;
}

float CD_pip_Eth_corr(float mom_, float theta_) {
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
float FD_pip_Eth_corr_lower(float mom_, float theta_) {
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

float FD_pip_Eth_corr_upper(float mom_, float theta_) {
  return theta_ +
         (0.00094724 * pow(mom_, 3) + (-0.00524101) * pow(mom_, 2) + 0.00919525 * mom_ + (-0.00516691)) *
             pow(theta_, 3) +
         ((-0.09887756) * pow(mom_, 3) + 0.54682169 * pow(mom_, 2) + (-0.95634115) * mom_ + 0.53345618) *
             pow(theta_, 2) +
         (3.44104365 * pow(mom_, 3) + (-19.02680178) * pow(mom_, 2) + 33.17864748 * mom_ + (-18.37813421)) *
             pow(theta_, 1) +
         (-39.82151866) * pow(mom_, 3) + 220.12521819 * pow(mom_, 2) + (-382.61957089) * mom_ + 210.34677439;
}

float CD_pip_Eph_corr(float mom_, float theta_, float phi_) {
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
float FD_pip_Eph_corr_lower(float mom_, float theta_, float phi_) {
  return phi_ +
         ((-4.86422409e-05) * pow(mom_, 4) + 1.21216530e-03 * pow(mom_, 3) + (-8.15266042e-03) * pow(mom_, 2) +
          1.93258907e-02 * mom_ + (-1.28009681e-02)) *
             pow(theta_, 1) +
         0.01081378 * pow(mom_, 4) + (-0.14401558) * pow(mom_, 3) + 0.69173611 * pow(mom_, 2) + (-1.3964496) * mom_ +
         0.95058901;
}
float FD_pip_Eph_corr_upper(float mom_, float theta_, float phi_) {
  return phi_ +
         ((-0.02343664) * pow(mom_, 3) + 0.13264734 * pow(mom_, 2) + (-0.2342437) * mom_ + 0.12601401) *
             pow(theta_, 1) +
         0.50037573 * pow(mom_, 3) + (-2.72628993) * pow(mom_, 2) + 4.48508987 * mom_ + (-2.05446324);
}

// energy loss corrections for pim

float CD_pim_Emom_corr(float mom_, float theta_) {
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
    return 2.27546950e-07 * pow(theta_, 3) + (-8.12537308e-05) * pow(theta_, 2) + 9.10902744e-03 * pow(theta_, 1) +
           (-3.22464750e-01);
  }
}
float FD_pim_Emom_corr_lower(float mom_, float theta_) { return mom_ + 0.00030448 * mom_ + 0.00232071; }
float FD_pim_Emom_corr_upper(float mom_, float theta_) { return mom_ + (-0.00100881) * mom_ + 0.00780439; }

float CD_pim_Eth_corr(float mom_, float theta_) {
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
float FD_pim_Eth_corr_lower(float mom_, float theta_) {
  return theta_ + ((-1.13685553e-04) * pow(mom_, 4) + 4.19458440e-03 * pow(mom_, 3) + (-3.76566663e-02) * pow(mom_, 2) +
                   1.30733557e-01 * pow(mom_, 1) + (-1.76073418e-01));
}

float FD_pim_Eth_corr_upper(float mom_, float theta_) {
  return theta_ +
         (0.01520214 * pow(mom_, 3) + (-0.08264195) * pow(mom_, 2) + 0.14545703 * mom_ + (-0.0888854)) *
             pow(theta_, 1) +
         (-0.46222418) * pow(mom_, 3) + 2.45741975 * pow(mom_, 2) + (-4.17396135) * mom_ + 2.39541974;
}

float CD_pim_Eph_corr(float mom_, float theta_, float phi_) {
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

float FD_pim_Eph_corr_lower(float mom_, float theta_, float phi_) {
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
         0.32282676 * pow(mom_, 3) + (-3.48574851) * pow(mom_, 3) + 13.11695944 * pow(mom_, 2) +
         (-19.9133663) * pow(mom_, 1) + 9.82183739;
}
float FD_pim_Eph_corr_upper(float mom_, float theta_, float phi_) {
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
         6.17963055 * pow(mom_, 3) + (-5.81705813) * pow(mom_, 3) + (-53.39466945) * pow(mom_, 2) +
         77.16020833 * pow(mom_, 1) + (-20.58824011);
}
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

}  // namespace mom_corr
