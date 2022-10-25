
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
    return mom_ + 2.27546950e-07 * pow(theta_, 3) + (-8.12537308e-05) * pow(theta_, 2) +
           9.10902744e-03 * pow(theta_, 1) + (-3.22464750e-01);
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
         0.32282676 * pow(mom_, 4) + (-3.48574851) * pow(mom_, 3) + 13.11695944 * pow(mom_, 2) +
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
         6.17963055 * pow(mom_, 4) + (-5.81705813) * pow(mom_, 3) + (-53.39466945) * pow(mom_, 2) +
         77.16020833 * pow(mom_, 1) + (-20.58824011);
}
//////////////////// new mom correction start

double dppC(float Px, float Py, float Pz, int sec, int ivec) {
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

// our hadron momentum correction come from here:
// proton mom corr
float alpha_prot_mom_corr_FD[2] = {0.6, 0.9};
float alpha_prot_mom_corr_CD[5] = {1.0, 0.5, 0.95};

// float alpha_prot_mom_corr_FD[2] = {1., 1.};
// float alpha_prot_mom_corr_CD[5] = {1.0, 1.0, 1.0};

double CDProt[3][4] = {{0.01887542, -0.02475295, -0.1015926, 0.03270923},
                       {0.08789515, -0.28722038, 0.26654857, -0.07150531},
                       {0.05834911, -0.1821812, 0.21029297, -0.06586461}};

float CD_prot_Hmom_corr(float mom_, float phi_) {
  if (phi_ > 270 || phi_ <= 30) {
    return mom_ - alpha_prot_mom_corr_CD[0] *
                      (CDProt[0][0] * pow(mom_, 3) + CDProt[0][1] * pow(mom_, 2) + CDProt[0][2] * mom_ + CDProt[0][3]);
  }
  else if (phi_ > 30 && phi_ <= 150) {
    return mom_ - alpha_prot_mom_corr_CD[1] *
                      (CDProt[1][0] * pow(mom_, 3) + CDProt[1][1] * pow(mom_, 2) + CDProt[1][2] * mom_ + CDProt[1][3]);
  }
  else if (phi_ > 150 && phi_ <= 270) {
    return mom_ - alpha_prot_mom_corr_CD[2] *
                      (CDProt[2][0] * pow(mom_, 3) + CDProt[2][1] * pow(mom_, 2) + CDProt[2][2] * mom_ + CDProt[2][3]);
  } else
    return NAN;
}

double FDProtL[6][4] = {
    {-0.00295384, 0.02153553, -0.05321597, 0.0246334}, {-0.00223025, 0.01588222, -0.04035392, 0.0281309},
    {-0.00217536, 0.015924, -0.04262257, 0.03281738},  {-0.00223025, 0.01588222, -0.04035392, 0.0281309},
    {-0.00211633, 0.0186769, -0.05749485, 0.03704521}, {-0.00421353, 0.03406236, -0.09387399, 0.05380671}};

float FD_prot_Hmom_corr_lower(float mom_, float dc_sec) {
  if (dc_sec == 1) {
    return mom_ - alpha_prot_mom_corr_FD[0] * (FDProtL[0][0] * pow(mom_, 3) + FDProtL[0][1] * pow(mom_, 2) +
                                               FDProtL[0][2] * mom_ + FDProtL[0][3]);
  }
  else if (dc_sec == 2) {
    return mom_ - alpha_prot_mom_corr_FD[0] * (FDProtL[1][0] * pow(mom_, 3) + FDProtL[1][1] * pow(mom_, 2) +
                                               FDProtL[1][2] * mom_ + FDProtL[1][3]);
  }
  else if (dc_sec == 3) {
    return mom_ - alpha_prot_mom_corr_FD[0] * (FDProtL[2][0] * pow(mom_, 3) + FDProtL[2][1] * pow(mom_, 2) +
                                               FDProtL[2][2] * mom_ + FDProtL[2][3]);
  }
  else if (dc_sec == 4) {
    return mom_ - alpha_prot_mom_corr_FD[0] * (FDProtL[3][0] * pow(mom_, 3) + FDProtL[3][1] * pow(mom_, 2) +
                                               FDProtL[3][2] * mom_ + FDProtL[3][3]);
  }
  else if (dc_sec == 5) {
    return mom_ - alpha_prot_mom_corr_FD[0] * (FDProtL[4][0] * pow(mom_, 3) + FDProtL[4][1] * pow(mom_, 2) +
                                               FDProtL[4][2] * mom_ + FDProtL[4][3]);
  }
  else if (dc_sec == 6) {
    return mom_ - alpha_prot_mom_corr_FD[0] * (FDProtL[5][0] * pow(mom_, 3) + FDProtL[5][1] * pow(mom_, 2) +
                                               FDProtL[5][2] * mom_ + FDProtL[0][3]);
  } else
    return NAN;
}

double FDProth[6][4] = {
    {-0.00647504, 0.05130311, -0.10250439, 0.05724377}, {0.01688694, -0.08105666, 0.12060526, -0.04750832},
    {-0.01177275, 0.06558347, -0.10800544, 0.07067348}, {0.0022607, -0.00595311, 0.00188426, 0.01132344},
    {-0.00149675, 0.0334335, -0.08967463, 0.06174838},  {-0.00422545, 0.037816, -0.08261621, 0.04522374}};

float FD_prot_Hmom_corr_upper(float mom_, float dc_sec) {
  if (dc_sec == 1) {
    return mom_ - alpha_prot_mom_corr_FD[1] * (FDProth[0][0] * pow(mom_, 3) + FDProth[0][1] * pow(mom_, 2) +
                                               FDProth[0][2] * mom_ + FDProth[0][3]);
  }
  else if (dc_sec == 2) {
    return mom_ - alpha_prot_mom_corr_FD[1] * (FDProth[1][0] * pow(mom_, 3) + FDProth[1][1] * pow(mom_, 2) +
                                               FDProth[1][2] * mom_ + FDProth[1][3]);
  }
  else if (dc_sec == 3) {
    return mom_ - alpha_prot_mom_corr_FD[1] * (FDProth[2][0] * pow(mom_, 3) + FDProth[2][1] * pow(mom_, 2) +
                                               FDProth[2][2] * mom_ + FDProth[2][3]);
  }
  else if (dc_sec == 4) {
    return mom_ - alpha_prot_mom_corr_FD[1] * (FDProth[3][0] * pow(mom_, 3) + FDProth[3][1] * pow(mom_, 2) +
                                               FDProth[3][2] * mom_ + FDProth[3][3]);
  }
  else if (dc_sec == 5) {
    return mom_ - alpha_prot_mom_corr_FD[1] * (FDProth[4][0] * pow(mom_, 3) + FDProth[4][1] * pow(mom_, 2) +
                                               FDProth[4][2] * mom_ + FDProth[4][3]);
  }
  else if (dc_sec == 6) {
    return mom_ - alpha_prot_mom_corr_FD[1] * (FDProth[5][0] * pow(mom_, 3) + FDProth[5][1] * pow(mom_, 2) +
                                               FDProth[5][2] * mom_ + FDProth[0][3]);
  } else
    return NAN;
}

/// pip hadron corrections

float alpha_pip_mom_corr_FD[2] = {0.5, 0.7};
float alpha_pip_mom_corr_CD[3] = {0.9, 0.45, 0.9};

// float alpha_pip_mom_corr_FD[2] = {0.50, 0.0};
// float alpha_pip_mom_corr_CD[3] = {0.0, 0.0, 0.0};

double CDPip[3][3] = {
    {0.04719538, -0.1493156, 0.02066552}, {0.01280761, -0.02654298, 0.01043671}, {0.01650979, 0.00771899, -0.00896883}};

float CD_pip_Hmom_corr(float mom_, float phi_) {
  if (phi_ > 270 || phi_ <= 30) {
    return mom_ - alpha_pip_mom_corr_CD[0] * (CDPip[0][0] * pow(mom_, 2) + CDPip[0][1] * mom_ + CDPip[0][2]);
  }
  else if (phi_ > 30 && phi_ <= 150) {
    return mom_ - alpha_pip_mom_corr_CD[1] * (CDPip[1][0] * pow(mom_, 2) + CDPip[1][1] * mom_ + CDPip[1][2]);
  }
  else if (phi_ > 150 && phi_ <= 270) {
    return mom_ - alpha_pip_mom_corr_CD[2] * (CDPip[2][0] * pow(mom_, 2) + CDPip[2][1] * mom_ + CDPip[2][2]);
  } else
    return NAN;
}

double FDPipL[6][4] = {{0.00121648, -0.00639932, 0.01722629, -0.02311194},
                       {0.0006273, -0.00401133, 0.01160091, -0.00776435},
                       {0.0012861, -0.011386, 0.03326134, -0.01889528},
                       {-4.35537266e-05, 2.88048986e-05, 1.16868864e-03, 2.25183481e-03},
                       {-0.00246769, 0.01826193, -0.04038504, 0.02252866},
                       {-0.0008343, 0.00713883, -0.01056021, -0.00563172}};

float FD_pip_Hmom_corr_lower(float mom_, float dc_sec) {
  if (dc_sec == 1) {
    return mom_ - alpha_pip_mom_corr_FD[0] *
                      (FDPipL[0][0] * pow(mom_, 3) + FDPipL[0][1] * pow(mom_, 2) + FDPipL[0][2] * mom_ + FDPipL[0][3]);
  }
  else if (dc_sec == 2) {
    return mom_ - alpha_pip_mom_corr_FD[0] *
                      (FDPipL[1][0] * pow(mom_, 3) + FDPipL[1][1] * pow(mom_, 2) + FDPipL[1][2] * mom_ + FDPipL[1][3]);
  }
  else if (dc_sec == 3) {
    return mom_ - alpha_pip_mom_corr_FD[0] *
                      (FDPipL[2][0] * pow(mom_, 3) + FDPipL[2][1] * pow(mom_, 2) + FDPipL[2][2] * mom_ + FDPipL[2][3]);
  }
  else if (dc_sec == 4) {
    return mom_ - alpha_pip_mom_corr_FD[0] *
                      (FDPipL[3][0] * pow(mom_, 3) + FDPipL[3][1] * pow(mom_, 2) + FDPipL[3][2] * mom_ + FDPipL[3][3]);
  }
  else if (dc_sec == 5) {
    return mom_ - alpha_pip_mom_corr_FD[0] *
                      (FDPipL[4][0] * pow(mom_, 3) + FDPipL[4][1] * pow(mom_, 2) + FDPipL[4][2] * mom_ + FDPipL[4][3]);
  }
  else if (dc_sec == 6) {
    return mom_ - alpha_pip_mom_corr_FD[0] *
                      (FDPipL[5][0] * pow(mom_, 3) + FDPipL[5][1] * pow(mom_, 2) + FDPipL[5][2] * mom_ + FDPipL[0][3]);
  }
  else
  return NAN;
}

double FDPipH[6][4] = {
    {-0.00356278, 0.01426312, 0.01091276, -0.02680428}, {0.01664683, -0.07549475, 0.1280444, -0.07204956},
    {0.00441982, -0.02617576, 0.06874411, -0.04127967}, {-0.0121876, 0.05399212, -0.04808956, 0.00589738},
    {-0.00209284, 0.01315519, -0.00717711, 0.00233726}, {0.02996566, -0.12590414, 0.18225925, -0.08344463}};

float FD_pip_Hmom_corr_upper(float mom_, float dc_sec) {
  if (dc_sec == 1) {
    return mom_ - alpha_pip_mom_corr_FD[1] *
                      (FDPipH[0][0] * pow(mom_, 3) + FDPipH[0][1] * pow(mom_, 2) + FDPipH[0][2] * mom_ + FDPipH[0][3]);
  }
  else if (dc_sec == 2) {
    return mom_ - alpha_pip_mom_corr_FD[1] *
                      (FDPipH[1][0] * pow(mom_, 3) + FDPipH[1][1] * pow(mom_, 2) + FDPipH[1][2] * mom_ + FDPipH[1][3]);
  }
  else if (dc_sec == 3) {
    return mom_ - alpha_pip_mom_corr_FD[1] *
                      (FDPipH[2][0] * pow(mom_, 3) + FDPipH[2][1] * pow(mom_, 2) + FDPipH[2][2] * mom_ + FDPipH[2][3]);
  }
  else if (dc_sec == 4) {
    return mom_ - alpha_pip_mom_corr_FD[1] *
                      (FDPipH[3][0] * pow(mom_, 3) + FDPipH[3][1] * pow(mom_, 2) + FDPipH[3][2] * mom_ + FDPipH[3][3]);
  }
  else if (dc_sec == 5) {
    return mom_ - alpha_pip_mom_corr_FD[1] *
                      (FDPipH[4][0] * pow(mom_, 3) + FDPipH[4][1] * pow(mom_, 2) + FDPipH[4][2] * mom_ + FDPipH[4][3]);
  }
  else if (dc_sec == 6) {
    return mom_ - alpha_pip_mom_corr_FD[1] *
                      (FDPipH[5][0] * pow(mom_, 3) + FDPipH[5][1] * pow(mom_, 2) + FDPipH[5][2] * mom_ + FDPipH[0][3]);
  } else
    return NAN;
}

/// pim hadron corrections
float alpha_pim_mom_corr_FD[2] = {0.3, 0.6};
float alpha_pim_mom_corr_CD[4] = {0.85, 1.0, 0.4};

// float alpha_pim_mom_corr_FD[2] = {1., 1.};
// float alpha_pim_mom_corr_CD[4] = {1., 1.0, 1.};

double CDPim[3][3] = {
    {0.03532859, -0.0343919, 0.00118764}, {0.02215662, 0.00657728, -0.00050019}, {0.03928169, -0.09742097, 0.02256432}};

float CD_pim_Hmom_corr(float mom_, float phi_) {
  if (phi_ > 270 || phi_ <= 30) {
    return mom_ - alpha_pim_mom_corr_CD[0] * (CDPim[0][0] * pow(mom_, 2) + CDPim[0][1] * mom_ + CDPim[0][2]);
  }
  else if (phi_ > 30 && phi_ <= 150) {
    return mom_ - alpha_pim_mom_corr_CD[1] * (CDPim[1][0] * pow(mom_, 2) + CDPim[1][1] * mom_ + CDPim[1][2]);
  }
  else if (phi_ > 150 && phi_ <= 270) {
    return mom_ - alpha_pim_mom_corr_CD[2] * (CDPim[2][0] * pow(mom_, 2) + CDPim[2][1] * mom_ + CDPim[2][2]);
  } else
    return NAN;
}

double FDPimL[6][4] = {
    {-0.00093805, 0.01419687, -0.05424683, 0.05540511}, {0.00101109, -0.00627509, 0.00820368, 0.0083081},
    {0.000657, -0.00192145, -0.00327336, 0.01824693},   {-0.0048234, 0.03892757, -0.09560444, 0.0720411},
    {-0.00028774, 0.00718568, -0.02632964, 0.02415357}, {-0.00022294, 0.00749173, -0.03113873, 0.0249824}};

float FD_pim_Hmom_corr_lower(float mom_, float dc_sec) {
  if (dc_sec == 1) {
    return mom_ - alpha_pim_mom_corr_FD[0] *
                      (FDPimL[0][0] * pow(mom_, 3) + FDPimL[0][1] * pow(mom_, 2) + FDPimL[0][2] * mom_ + FDPimL[0][3]);
  }
  else if (dc_sec == 2) {
    return mom_ - alpha_pim_mom_corr_FD[0] *
                      (FDPimL[1][0] * pow(mom_, 3) + FDPimL[1][1] * pow(mom_, 2) + FDPimL[1][2] * mom_ + FDPimL[1][3]);
  }
  else if (dc_sec == 3) {
    return mom_ - alpha_pim_mom_corr_FD[0] *
                      (FDPimL[2][0] * pow(mom_, 3) + FDPimL[2][1] * pow(mom_, 2) + FDPimL[2][2] * mom_ + FDPimL[2][3]);
  }
  else if (dc_sec == 4) {
    return mom_ - alpha_pim_mom_corr_FD[0] *
                      (FDPimL[3][0] * pow(mom_, 3) + FDPimL[3][1] * pow(mom_, 2) + FDPimL[3][2] * mom_ + FDPimL[3][3]);
  }
  else if (dc_sec == 5) {
    return mom_ - alpha_pim_mom_corr_FD[0] *
                      (FDPimL[4][0] * pow(mom_, 3) + FDPimL[4][1] * pow(mom_, 2) + FDPimL[4][2] * mom_ + FDPimL[4][3]);
  }
  else if (dc_sec == 6) {
    return mom_ - alpha_pim_mom_corr_FD[0] *
                      (FDPimL[5][0] * pow(mom_, 3) + FDPimL[5][1] * pow(mom_, 2) + FDPimL[5][2] * mom_ + FDPimL[0][3]);
  } else
    return NAN;
}

double FDPimH[6][5] = {{-0.03708776, 0.28179957, -0.75582601, 0.82592834, -0.30001301},
                       {-0.02561556, 0.18208368, -0.46698648, 0.51041072, -0.1898823},
                       {-0.02949599, 0.2299059, -0.63540584, 0.72691335, -0.27715446},
                       {-0.02387114, 0.18711438, -0.51217348, 0.57264114, -0.21573554},
                       {-0.0132594, 0.10936398, -0.31834157, 0.38727052, -0.16851244},
                       {-0.00052433, 0.00753084, -0.02917275, 0.03994074, -0.02118046}};

float FD_pim_Hmom_corr_upper(float mom_, float dc_sec) {
  if (dc_sec == 1) {
    return mom_ - alpha_pim_mom_corr_FD[1] * (FDPimH[0][0] * pow(mom_, 4) + FDPimH[0][1] * pow(mom_, 3) +
                                              FDPimH[0][2] * pow(mom_, 2) + FDPimH[0][3] * mom_ + FDPimH[0][4]);
  }
  else if (dc_sec == 2) {
    return mom_ - alpha_pim_mom_corr_FD[1] * (FDPimH[1][0] * pow(mom_, 4) + FDPimH[1][1] * pow(mom_, 3) +
                                              FDPimH[1][2] * pow(mom_, 2) + FDPimH[1][3] * mom_ + FDPimH[1][4]);
  }
  else if (dc_sec == 3) {
    return mom_ - alpha_pim_mom_corr_FD[1] * (FDPimH[2][0] * pow(mom_, 4) + FDPimH[2][1] * pow(mom_, 3) +
                                              FDPimH[2][2] * pow(mom_, 2) + FDPimH[2][3] * mom_ + FDPimH[2][4]);
  }
  else if (dc_sec == 4) {
    return mom_ - alpha_pim_mom_corr_FD[1] * (FDPimH[3][0] * pow(mom_, 4) + FDPimH[3][1] * pow(mom_, 3) +
                                              FDPimH[3][2] * pow(mom_, 2) + FDPimH[3][3] * mom_ + FDPimH[3][4]);
  }
  else if (dc_sec == 5) {
    return mom_ - alpha_pim_mom_corr_FD[1] * (FDPimH[4][0] * pow(mom_, 4) + FDPimH[4][1] * pow(mom_, 3) +
                                              FDPimH[4][2] * pow(mom_, 2) + FDPimH[4][3] * mom_ + FDPimH[4][4]);
  }
  else if (dc_sec == 6) {
    return mom_ - alpha_pim_mom_corr_FD[1] * (FDPimH[5][0] * pow(mom_, 4) + FDPimH[5][1] * pow(mom_, 3) +
                                              FDPimH[5][2] * pow(mom_, 2) + FDPimH[5][3] * mom_ + FDPimH[5][4]);
  } else
    return NAN;
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

////////////////// new new mom corr done (Aug-15-2022)

}  // namespace mom_corr

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
// //     if (_is_lower_band)
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
