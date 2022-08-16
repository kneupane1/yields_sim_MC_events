/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "reaction.hpp"

Reaction::Reaction(const std::shared_ptr<Branches12>& data, float beam_energy) {
  _data = data;
  _beam = std::make_unique<TLorentzVector>();
  _beam_energy = beam_energy;
  _sector = data->dc_sec(0);

  _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

  _gamma = std::make_unique<TLorentzVector>();
  _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
  _elec = std::make_unique<TLorentzVector>();
  this->SetElec();

  _mom_corr_elec = std::make_unique<TLorentzVector>();
  _mom_corr_pim = std::make_unique<TLorentzVector>();
  _mom_corr_pim_th = std::make_unique<TLorentzVector>();
  _mom_corr_pim_ph = std::make_unique<TLorentzVector>();
  _mom_corr_pip = std::make_unique<TLorentzVector>();
  _mom_corr_pip_th = std::make_unique<TLorentzVector>();
  _mom_corr_pip_ph = std::make_unique<TLorentzVector>();
  _mom_corr_prot = std::make_unique<TLorentzVector>();
  _mom_corr_prot_th = std::make_unique<TLorentzVector>();
  _mom_corr_prot_ph = std::make_unique<TLorentzVector>();
  _Energy_loss_uncorr_prot = std::make_unique<TLorentzVector>();
  _Energy_loss_uncorr_pip = std::make_unique<TLorentzVector>();
  _Energy_loss_uncorr_pim = std::make_unique<TLorentzVector>();
  _pim_tmt = std::make_unique<TLorentzVector>();
  _pip_tmt = std::make_unique<TLorentzVector>();

  _prot = std::make_unique<TLorentzVector>();
  _pip = std::make_unique<TLorentzVector>();
  _pim = std::make_unique<TLorentzVector>();
  _other = std::make_unique<TLorentzVector>();
  _neutron = std::make_unique<TLorentzVector>();

  for (int isec_mom_corr = 0; isec_mom_corr < 6; isec_mom_corr++) {
    for (int ivec = 0; ivec < 3; ivec++) {
      double dp1 = xx[ipar++], dp5 = xx[ipar++], dp9 = xx[ipar++];

      pars[isec_mom_corr][ivec][0] = (dp1 - 2 * dp5 + dp9) / 32.;
      pars[isec_mom_corr][ivec][1] = (-7 * dp1) / 16. + (5 * dp5) / 8. - (3 * dp9) / 16.;
      pars[isec_mom_corr][ivec][2] = (45 * dp1) / 32. - (9 * dp5) / 16. + (5 * dp9) / 32.;
    }
  }
}

Reaction::~Reaction() {}

//////////////////// new mom correction start

double Reaction::dppC(float Px, float Py, float Pz, int sec, int ivec) {
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
      // The following lines should be added up in the order given for the full correction
      // Applying this code as given will give the exact corrections of this analysis
      // These parameters will be combined into a single line at a later point

      dp = ((-6.2e-07) * phi * phi + (-2.4e-06) * phi + (3.3867e-04)) * pp * pp +
           ((8.66e-06) * phi * phi + (3.8389e-04) * phi + (-5.732e-03)) * pp +
           ((-1.607e-05) * phi * phi + (-4.9186e-04) * phi + (0.01490708));

      dp = dp + ((2.12e-06) * phi * phi + (-6.48e-06) * phi + (-3.0978e-04)) * pp * pp +
           ((-1.717e-05) * phi * phi + (6.22e-06) * phi + (2.892e-03)) * pp +
           ((2.56e-05) * phi * phi + (6.209e-05) * phi + (-3.692e-03));

      dp = dp + ((-1.3e-06) * phi * phi + (-2.01e-06) * phi + (3.6961e-04)) * pp * pp +
           ((1.103e-05) * phi * phi + (4.46e-06) * phi + (-3.5553e-03)) * pp +
           ((-2.272e-05) * phi * phi + (-1.788e-05) * phi + (8.2909e-03));

      dp = dp + ((-4.1e-07) * phi * phi + (-9.8e-07) * phi + (-1.71e-06)) * pp * pp +
           ((3.13e-06) * phi * phi + (-1.651e-05) * phi + (3.055e-04)) * pp +
           ((2.5e-07) * phi * phi + (-2.298e-05) * phi + (-1.09082e-03));

      dp = dp + ((-3.1e-07) * phi * phi + (-1.96e-06) * phi + (7.5e-05)) * pp * pp +
           ((2.68e-06) * phi * phi + (1.043e-05) * phi + (-7.2339e-04)) * pp +
           ((-3.51e-06) * phi * phi + (-2.994e-05) * phi + (1.48688e-03));
    }

    if (sec == 2) {
      // The following lines should be added up in the order given for the full correction
      // Applying this code as given will give the exact corrections of this analysis
      // These parameters will be combined into a single line at a later point

      dp = ((2.38e-06) * phi * phi + (6.164e-05) * phi + (1.7685e-04)) * pp * pp +
           ((-1.336e-05) * phi * phi + (-1.0139e-04) * phi + (-3.65463e-03)) * pp +
           ((8.22e-06) * phi * phi + (-2.3904e-04) * phi + (6.66967e-03));

      dp = dp + ((-2.99e-06) * phi * phi + (-3.46e-06) * phi + (8.3225e-04)) * pp * pp +
           ((2.08e-05) * phi * phi + (-8.059e-05) * phi + (-7.9052e-03)) * pp +
           ((-2.799e-05) * phi * phi + (3.0292e-04) * phi + (0.0169132));

      dp = dp + ((-1.2e-07) * phi * phi + (-2.241e-05) * phi + (-1.3058e-04)) * pp * pp +
           ((-1.21e-06) * phi * phi + (1.3224e-04) * phi + (1.5422e-03)) * pp +
           ((7.43e-06) * phi * phi + (-1.3034e-04) * phi + (-3.34406e-03));

      dp = dp + ((-9.7e-07) * phi * phi + (-1.14e-06) * phi + (2.3808e-04)) * pp * pp +
           ((8.14e-06) * phi * phi + (3.86e-06) * phi + (-2.41185e-03)) * pp +
           ((-1.481e-05) * phi * phi + (-3.1e-05) * phi + (5.53527e-03));

      dp = dp + ((-1.8e-07) * phi * phi + (-1.6e-06) * phi + (1.65e-05)) * pp * pp +
           ((1.32e-06) * phi * phi + (6.14e-06) * phi + (-1.5742e-04)) * pp +
           ((-1.88e-06) * phi * phi + (-8.92e-06) * phi + (3.7879e-04));
    }

    if (sec == 3) {
      // The following lines should be added up in the order given for the full correction
      // Applying this code as given will give the exact corrections of this analysis
      // These parameters will be combined into a single line at a later point

      dp = ((2.49e-06) * phi * phi + (-4.841e-05) * phi + (2.8506e-04)) * pp * pp +
           ((-9.06e-06) * phi * phi + (1.2071e-04) * phi + (-3.07157e-03)) * pp +
           ((4.81e-06) * phi * phi + (2.85e-06) * phi + (8.91138e-03));

      dp = dp + ((-2.35e-06) * phi * phi + (3.44e-05) * phi + (3.4527e-04)) * pp * pp +
           ((2.068e-05) * phi * phi + (-2.005e-04) * phi + (-4.6252e-03)) * pp +
           ((-4.26e-05) * phi * phi + (2.2642e-04) * phi + (0.01231));

      dp = dp + ((-9e-07) * phi * phi + (1.367e-05) * phi + (1.6897e-04)) * pp * pp +
           ((7.27e-06) * phi * phi + (-9.074e-05) * phi + (-1.39503e-03)) * pp +
           ((-1.275e-05) * phi * phi + (1.316e-04) * phi + (2.2438e-03));

      dp = dp + ((8.7e-07) * phi * phi + (-8e-06) * phi + (-8.044e-05)) * pp * pp +
           ((-9.43e-06) * phi * phi + (9.486e-05) * phi + (8.8623e-04)) * pp +
           ((2.283e-05) * phi * phi + (-2.1403e-04) * phi + (-2.10146e-03));

      dp = dp + ((1.3e-07) * phi * phi + (-2.06e-06) * phi + (-1.022e-05)) * pp * pp +
           ((-1.46e-06) * phi * phi + (2.411e-05) * phi + (8.388e-05)) * pp +
           ((3.51e-06) * phi * phi + (-5.756e-05) * phi + (-4.141e-05));
    }

    if (sec == 4) {
      // The following lines should be added up in the order given for the full correction
      // Applying this code as given will give the exact corrections of this analysis
      // These parameters will be combined into a single line at a later point

      dp = ((3.03e-06) * phi * phi + (-3.98e-06) * phi + (1.24811e-03)) * pp * pp +
           ((-8.87e-06) * phi * phi + (3.0e-07) * phi + (-0.01203987)) * pp +
           ((-8.3e-07) * phi * phi + (1.2562e-04) * phi + (0.02730945));

      dp = dp + ((-2.25e-06) * phi * phi + (-9.05e-06) * phi + (3.837e-05)) * pp * pp +
           ((1.656e-05) * phi * phi + (5.745e-05) * phi + (-1.3013e-03)) * pp +
           ((-3.035e-05) * phi * phi + (-8.443e-05) * phi + (5.295e-03));

      dp = dp + ((3.1e-07) * phi * phi + (-1.313e-05) * phi + (-2.7121e-04)) * pp * pp +
           ((-3.49e-06) * phi * phi + (9.059e-05) * phi + (2.3969e-03)) * pp +
           ((6.43e-06) * phi * phi + (-9.583e-05) * phi + (-4.0804e-03));

      dp = dp + ((-4.8e-07) * phi * phi + (-3.41e-06) * phi + (1.9609e-04)) * pp * pp +
           ((1.56e-06) * phi * phi + (2.632e-05) * phi + (-1.44602e-03)) * pp +
           ((3.17e-06) * phi * phi + (-2.63e-05) * phi + (2.6809e-03));

      dp = dp + ((-6.5e-07) * phi * phi + (-6.33e-06) * phi + (1.101e-04)) * pp * pp +
           ((4.47e-06) * phi * phi + (4.733e-05) * phi + (-9.1397e-04)) * pp +
           ((-6.43e-06) * phi * phi + (-7.666e-05) * phi + (1.59452e-03));
    }

    if (sec == 5) {
      // The following lines should be added up in the order given for the full correction
      // Applying this code as given will give the exact corrections of this analysis
      // These parameters will be combined into a single line at a later point

      dp = ((6.33e-06) * phi * phi + (1.607e-05) * phi + (-2.2161e-04)) * pp * pp +
           ((-5.05e-05) * phi * phi + (-9.049e-05) * phi + (9.6753e-04)) * pp +
           ((1.0928e-04) * phi * phi + (1.2416e-04) * phi + (-0.01023377));

      dp = dp + ((-4.5e-07) * phi * phi + (-6.42e-06) * phi + (5.021e-05)) * pp * pp +
           ((4.9e-06) * phi * phi + (5.907e-05) * phi + (-1.0226e-03)) * pp +
           ((-1.21e-05) * phi * phi + (-1.1515e-04) * phi + (3.9894e-03));

      dp = dp + ((6e-07) * phi * phi + (-1.054e-05) * phi + (-4.725e-05)) * pp * pp +
           ((-5.61e-06) * phi * phi + (1.0622e-04) * phi + (6.8259e-04)) * pp +
           ((1.065e-05) * phi * phi + (-2.1664e-04) * phi + (-1.5396e-03));

      dp = dp + ((5.5e-07) * phi * phi + (6.89e-06) * phi + (-7.51e-06)) * pp * pp +
           ((-4.29e-06) * phi * phi + (-7.935e-05) * phi + (-7.698e-05)) * pp +
           ((7.89e-06) * phi * phi + (1.5742e-04) * phi + (1.52222e-03));

      dp = dp + ((-4.33e-06) * phi * phi + (-9.7e-07) * phi + (1.82284e-03)) * pp * pp +
           ((4.27e-05) * phi * phi + (2.56e-06) * phi + (-0.01770834)) * pp +
           ((-9.481e-05) * phi * phi + (8.81e-06) * phi + (0.03880514));
    }

    if (sec == 6) {
      // The following lines should be added up in the order given for the full correction
      // Applying this code as given will give the exact corrections of this analysis
      // These parameters will be combined into a single line at a later point

      dp = ((4.87e-06) * phi * phi + (-8.83e-05) * phi + (1.94282e-03)) * pp * pp +
           ((-2.412e-05) * phi * phi + (6.0476e-04) * phi + (-0.01370525)) * pp +
           ((2.318e-05) * phi * phi + (-6.665e-04) * phi + (0.03056206));

      dp = dp + ((-5.9e-07) * phi * phi + (1.391e-05) * phi + (-7.5974e-04)) * pp * pp +
           ((2.62e-06) * phi * phi + (-6.794e-05) * phi + (5.3699e-03)) * pp +
           ((-3.38e-06) * phi * phi + (7.579e-05) * phi + (-6.5997e-03));

      dp = dp + ((-6.1e-07) * phi * phi + (-1.245e-05) * phi + (3.0939e-04)) * pp * pp +
           ((4.26e-06) * phi * phi + (9.574e-05) * phi + (-2.08024e-03)) * pp +
           ((-3.25e-06) * phi * phi + (-1.0427e-04) * phi + (9.967e-04));

      dp = dp + ((-7.8e-07) * phi * phi + (8.1e-06) * phi + (8.971e-05)) * pp * pp +
           ((4.56e-06) * phi * phi + (-2.082e-05) * phi + (-1.1163e-03)) * pp +
           ((-7.55e-06) * phi * phi + (-7.662e-05) * phi + (4.17743e-03));

      dp = dp + ((-7.6e-07) * phi * phi + (3.84e-06) * phi + (1.7347e-04)) * pp * pp +
           ((5.31e-06) * phi * phi + (-2.952e-05) * phi + (-1.265e-03)) * pp +
           ((-8.51e-06) * phi * phi + (4.907e-05) * phi + (2.01338e-03));
    }
  }

  //==========//  PARTICLE = PI+ PION (END)  //==========//

  //==========//  PARTICLE = PI- PION  //==========//

  if (ivec == 2) {
    if (sec == 1) {
      // std::cout << " " << sec << " " << Phi << std::endl;
      // std::cout << " " << sec << " " << Phi << std::endl;

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
    if (sec == 1) {
      // The following lines should be added up in the order given for the full correction
      // Applying this code as given will give the exact corrections of this analysis
      // These parameters will be combined into a single line at a later point

      dp = ((3.54e-06) * phi * phi + (0.00012741) * phi + (-0.00169485)) * pp * pp +
           ((-6.8e-06) * phi * phi + (-0.00018409) * phi + (0.00756841)) * pp +
           ((3.49e-06) * phi * phi + (6.297e-05) * phi + (-0.00334762));
      dp = dp + ((-1.857e-05) * phi * phi + (-1.14e-04) * phi + (4.346e-03)) * pp * pp +
           ((5.978e-05) * phi * phi + (3.776e-04) * phi + (-0.0192)) * pp +
           ((-3.864e-05) * phi * phi + (-2.651e-04) * phi + (0.0163));
      dp = dp + ((9.05e-06) * phi * phi + (-0.00011072) * phi + (-0.0019497)) * pp * pp +
           ((-2.473e-05) * phi * phi + (0.00025581) * phi + (0.006382)) * pp +
           ((1.288e-05) * phi * phi + (-0.00012678) * phi + (-0.0034358));
      dp = dp + ((-4.88e-06) * phi * phi + (0.00015119) * phi + (0.00157759)) * pp * pp +
           ((1.467e-05) * phi * phi + (-0.00032185) * phi + (-0.00770597)) * pp +
           ((-1.062e-05) * phi * phi + (0.00016935) * phi + (0.00680068));
      dp = dp + ((-6.48e-06) * phi * phi + (9.682e-05) * phi + (0.00122778)) * pp * pp +
           ((2.14e-05) * phi * phi + (-0.00026216) * phi + (-0.00452799)) * pp +
           ((-1.54e-05) * phi * phi + (0.00013753) * phi + (0.00366475));
    }

    if (sec == 2) {
      // The following lines should be added up in the order given for the full correction
      // Applying this code as given will give the exact corrections of this analysis
      // These parameters will be combined into a single line at a later point

      dp = ((-7.51e-06) * phi * phi + (7.449e-05) * phi + (-0.00193523)) * pp * pp +
           ((1.372e-05) * phi * phi + (-0.00033226) * phi + (0.01605589)) * pp +
           ((-2.27e-06) * phi * phi + (0.00017048) * phi + (-0.01159601));
      dp = dp + ((-2e-07) * phi * phi + (1.393e-04) * phi + (-9.325e-04)) * pp * pp +
           ((3.51e-06) * phi * phi + (-4.582e-04) * phi + (-2.414e-03)) * pp +
           ((1.13e-06) * phi * phi + (3.343e-04) * phi + (4.76e-03));
      dp = dp + ((3.84e-06) * phi * phi + (-1.841e-05) * phi + (0.00085236)) * pp * pp +
           ((-1.049e-05) * phi * phi + (3.949e-05) * phi + (-0.0024575)) * pp +
           ((7.08e-06) * phi * phi + (-2.962e-05) * phi + (0.00084752));
      dp = dp + ((1.51e-06) * phi * phi + (-8.957e-05) * phi + (-0.0029601)) * pp * pp +
           ((-1.07e-05) * phi * phi + (0.00022399) * phi + (0.0077783)) * pp +
           ((1.042e-05) * phi * phi + (-0.00015105) * phi + (-0.00291493));
      dp = dp + ((2.57e-06) * phi * phi + (-1.635e-05) * phi + (-0.00105274)) * pp * pp +
           ((-8.65e-06) * phi * phi + (4.65e-05) * phi + (0.00267734)) * pp +
           ((5.91e-06) * phi * phi + (-2.431e-05) * phi + (-0.00144716));
    }

    if (sec == 3) {
      // The following lines should be added up in the order given for the full correction
      // Applying this code as given will give the exact corrections of this analysis
      // These parameters will be combined into a single line at a later point

      dp = ((-1.192e-05) * phi * phi + (5.745e-05) * phi + (0.00184011)) * pp * pp +
           ((3.226e-05) * phi * phi + (-0.00017688) * phi + (0.00281795)) * pp +
           ((-1.339e-05) * phi * phi + (0.00010705) * phi + (-0.00584473));
      dp = dp + ((1.153e-05) * phi * phi + (-4.931e-05) * phi + (-5.09e-03)) * pp * pp +
           ((-3.071e-05) * phi * phi + (7.061e-05) * phi + (6.413e-03)) * pp +
           ((1.525e-05) * phi * phi + (3.982e-05) * phi + (2.289e-03));
      dp = dp + ((2.49e-06) * phi * phi + (-2.58e-05) * phi + (-0.00047769)) * pp * pp +
           ((-6.8e-06) * phi * phi + (7.153e-05) * phi + (0.0012909)) * pp +
           ((3.31e-06) * phi * phi + (-5.431e-05) * phi + (-0.00023919));
      dp = dp + ((4.49e-06) * phi * phi + (-0.00019179) * phi + (0.00160744)) * pp * pp +
           ((-8.79e-06) * phi * phi + (0.00055794) * phi + (-0.0092487)) * pp +
           ((1.43e-06) * phi * phi + (-0.00031925) * phi + (0.00929164));
      dp = dp + ((7.71e-06) * phi * phi + (-2.683e-05) * phi + (-0.00231227)) * pp * pp +
           ((-2.053e-05) * phi * phi + (3.837e-05) * phi + (0.00664668)) * pp +
           ((1.225e-05) * phi * phi + (-7.91e-06) * phi + (-0.00462877));
    }

    if (sec == 4) {
      // The following lines should be added up in the order given for the full correction
      // Applying this code as given will give the exact corrections of this analysis
      // These parameters will be combined into a single line at a later point

      dp = ((1.02e-06) * phi * phi + (-5.753e-05) * phi + (-0.00116896)) * pp * pp +
           ((6.6e-07) * phi * phi + (0.00027025) * phi + (0.01147826)) * pp +
           ((5.7e-07) * phi * phi + (-0.00018222) * phi + (-0.01011935));
      dp = dp + ((-1.496e-05) * phi * phi + (5.225e-05) * phi + (1.183e-03)) * pp * pp +
           ((4.374e-05) * phi * phi + (-1.675e-04) * phi + (-8.0593e-03)) * pp +
           ((-3.157e-05) * phi * phi + (1.325e-04) * phi + (9.979e-03));
      dp = dp + ((5e-07) * phi * phi + (-1.454e-05) * phi + (-0.0007149)) * pp * pp +
           ((5.47e-06) * phi * phi + (7.743e-05) * phi + (0.00068488)) * pp +
           ((-7.88e-06) * phi * phi + (-8.695e-05) * phi + (7.4e-05));
      dp = dp + ((7.8e-07) * phi * phi + (0.00011795) * phi + (0.00165541)) * pp * pp +
           ((-3.13e-06) * phi * phi + (-0.00042059) * phi + (-0.00798432)) * pp +
           ((6.23e-06) * phi * phi + (0.00033604) * phi + (0.00680092));
      dp = dp + ((-3.98e-06) * phi * phi + (-2.054e-05) * phi + (0.00014266)) * pp * pp +
           ((1.058e-05) * phi * phi + (6.702e-05) * phi + (-0.00055213)) * pp +
           ((-8.63e-06) * phi * phi + (-8.036e-05) * phi + (0.00105203));
    }

    if (sec == 5) {
      // The following lines should be added up in the order given for the full correction
      // Applying this code as given will give the exact corrections of this analysis
      // These parameters will be combined into a single line at a later point

      dp = ((8.52e-06) * phi * phi + (-0.00019083) * phi + (-0.00545608)) * pp * pp +
           ((-2.943e-05) * phi * phi + (0.00051946) * phi + (0.01909462)) * pp +
           ((2.152e-05) * phi * phi + (-0.00028127) * phi + (-0.01224992));
      dp = dp + ((-7.62e-06) * phi * phi + (4.846e-05) * phi + (1.713e-03)) * pp * pp +
           ((1.092e-05) * phi * phi + (-4.886e-05) * phi + (-8.412e-03)) * pp +
           ((5.17e-06) * phi * phi + (-2.448e-05) * phi + (5.572e-03));
      dp = dp + ((4.63e-06) * phi * phi + (8.79e-06) * phi + (-0.0016709)) * pp * pp +
           ((-1.308e-05) * phi * phi + (-3.022e-05) * phi + (0.0046265)) * pp +
           ((7.58e-06) * phi * phi + (1.813e-05) * phi + (-0.0026691));
      dp = dp + ((-1.6e-07) * phi * phi + (-3.7e-06) * phi + (0.00249659)) * pp * pp +
           ((9e-06) * phi * phi + (-2.182e-05) * phi + (-0.012996)) * pp +
           ((-1.66e-05) * phi * phi + (-1.649e-05) * phi + (0.01366649));
      dp = dp + ((-3.89e-06) * phi * phi + (-8.14e-06) * phi + (0.00086004)) * pp * pp +
           ((1.439e-05) * phi * phi + (4.012e-05) * phi + (-0.00281493)) * pp +
           ((-1.111e-05) * phi * phi + (-3.062e-05) * phi + (0.00129624));
    }

    if (sec == 6) {
      // The following lines should be added up in the order given for the full correction
      // Applying this code as given will give the exact corrections of this analysis
      // These parameters will be combined into a single line at a later point

      dp = ((-6.8e-06) * phi * phi + (-9.966e-05) * phi + (0.00217971)) * pp * pp +
           ((1.846e-05) * phi * phi + (0.000334) * phi + (-0.00358955)) * pp +
           ((-7.83e-06) * phi * phi + (-0.00022134) * phi + (0.00280734));
      dp = dp + ((-3.43e-06) * phi * phi + (3.171e-04) * phi + (-1.572e-03)) * pp * pp +
           ((1.799e-05) * phi * phi + (-8.219e-04) * phi + (-2.945e-03)) * pp +
           ((-9.34e-06) * phi * phi + (3.74e-04) * phi + (5.1794e-03));
      dp = dp + ((-3.64e-06) * phi * phi + (1.113e-05) * phi + (0.0014469)) * pp * pp +
           ((1.031e-05) * phi * phi + (-2.765e-05) * phi + (-0.0041404)) * pp +
           ((-5.83e-06) * phi * phi + (-3.08e-06) * phi + (0.0024571));
      dp = dp + ((1.969e-05) * phi * phi + (-6.395e-05) * phi + (-0.00457687)) * pp * pp +
           ((-5.808e-05) * phi * phi + (0.00018976) * phi + (0.01091338)) * pp +
           ((3.344e-05) * phi * phi + (-0.00015701) * phi + (-0.00426004));
      dp = dp + ((2.39e-06) * phi * phi + (3.47e-06) * phi + (0.00087654)) * pp * pp +
           ((-7.18e-06) * phi * phi + (4.082e-05) * phi + (-0.00295533)) * pp +
           ((2.94e-06) * phi * phi + (-4.526e-05) * phi + (0.00265678));
    }
  }

  // ==========//  PARTICLE = PROTON (END)  //==========//

  return dp / pp;
};

// Below shows how the corrections are to be applied using the ROOT momentum 4-vector using the above code:
// auto fe = dppC(ex, ey, ez, esec, 0) + 1;
// auto fpip = dppC(pipx, pipy, pipz, pipsec, 1) + 1;
// auto fpim = dppC(pimx, pimy, pimz, pimsec, 2) + 1;
// auto fpro = dppC(prox, proy, proz, prosec, 3) + 1;

// auto eleC = ROOT::Math::PxPyPzMVector(ex * fe, ey * fe, ez * fe, 0);
// auto pipC = ROOT::Math::PxPyPzMVector(pipx * fpip, pipy * fpip, pipz * fpip, 0.13957);
// auto pimC = ROOT::Math::PxPyPzMVector(pimx * fpim, pimy * fpim, pimz * fpim, 0.13957);
// auto proC = ROOT::Math::PxPyPzMVector(prox * fpro, proy * fpro, proz * fpro, 0.938);

////////////////// new mom corr done

//////////////////// old mom corrections (probably better one)
// double Reaction::dpp(float px, float py, float pz, int sec_mom_corr, int ivec) {
//   double pp = sqrt(px * px + py * py + pz * pz);

//   double a = pars[sec_mom_corr - 1][ivec][0], b = pars[sec_mom_corr - 1][ivec][1],
//          c = pars[sec_mom_corr - 1][ivec][2];

//   double dp = a * pp * pp + b * pp + c;  // pol2 corr func

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

void Reaction::SetElec() {
  _hasE = true;
  _elec->SetXYZM(_data->px(0), _data->py(0), _data->pz(0), MASS_E);
  // *_gamma += *_beam - *_elec;  // be careful you are commenting this only to include the momentum correction

  // // // // Can calculate W and Q2 here (useful for simulations as sim do not have elec mom corrections)
  _W = physics::W_calc(*_beam, *_elec);
  _Q2 = physics::Q2_calc(*_beam, *_elec);

  // // //One way of  calculating mom - corrected four vectors
  // //   // // _cx = _data->px(0)/_elec->P();
  // //   // // _cy = _data->py(0) / _elec->P();
  // //   // // _cz = _data->pz(0) / _elec->P();
  // //   // _elec_mom_corrected = _elec->P() * (dpp(_data->px(0), _data->py(0), _data->pz(0), _data->dc_sec(0),
  // 0) + 1);

  // //   // _px_prime_elec = _cx * _elec_mom_corrected;
  // //   // _py_prime_elec = _cy * _elec_mom_corrected;
  // //   // _pz_prime_elec = _cz * _elec_mom_corrected; // _mom_corr_elec->SetXYZM(_px_prime_elec,
  // _py_prime_elec,
  // //   // _pz_prime_elec, MASS_E);

  _elec_mom = _elec->P();
  // _elec_mom_corrected = (dpp(_data->px(0), _data->py(0), _data->pz(0), _data->dc_sec(0), 0) + 1);

  //   _mom_corr_elec->SetPxPyPzE(_data->px(0) * _elec_mom_corrected, _data->py(0) * _elec_mom_corrected,
  //                              _data->pz(0) * _elec_mom_corrected, _elec_mom * _elec_mom_corrected);

  // *_gamma += *_beam - *_mom_corr_elec;

  // _W = physics::W_calc(*_beam, *_mom_corr_elec);
  // _Q2 = physics::Q2_calc(*_beam, *_mom_corr_elec);
  // _P_elec = _elec->P();
}
void Reaction::SetMomCorrElec() {
  // Below shows how the corrections are to be applied using the ROOT momentum 4-vector using the above code:

  // New electron momentum corrections
  fe = dppC(_data->px(0), _data->py(0), _data->pz(0), _data->dc_sec(0), 0) + 1;
  _mom_corr_elec->SetXYZM(_data->px(0) * fe, _data->py(0) * fe, _data->pz(0) * fe,
                          MASS_E);  // this is new electron mom corrections aug 2022

  // _elec_mom_corrected = (dpp(_data->px(0), _data->py(0), _data->pz(0), _data->dc_sec(0), 0) + 1);
  // _mom_corr_elec->SetXYZM(_data->px(0) * _elec_mom_corrected, _data->py(0) * _elec_mom_corrected,
  //                         _data->pz(0) * _elec_mom_corrected, MASS_E);

  // _mom_corr_elec->SetPxPyPzE(_data->px(0) * _elec_mom_corrected, _data->py(0) * _elec_mom_corrected,
  //                            _data->pz(0) * _elec_mom_corrected, _elec_mom * _elec_mom_corrected);

  *_gamma += *_beam - *_mom_corr_elec;
  _W_after = physics::W_calc(*_beam, *_mom_corr_elec);
  // _W = physics::W_calc(*_beam, *_mom_corr_elec);
  // _Q2 = physics::Q2_calc(*_beam, *_mom_corr_elec);

  _P_elec = _mom_corr_elec->P();

  // _E_elec = _mom_corr_elec->E();
}
double Reaction::Corr_elec_mom() {
  if (_P_elec != _P_elec) SetMomCorrElec();
  // std::cout << " elec mom corrected " << _elec_mom_corrected << std::endl;

  return _P_elec;
}

double Reaction::elec_mom() {
  if (_elec_mom != _elec_mom) SetElec();
  // std::cout << " emec mom " << _elec_mom << std::endl;

  return _elec_mom;
}

void Reaction::SetProton(int i) {
  _numProt++;
  _numPos++;
  _hasP = true;

  // _prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
  _Energy_loss_uncorr_prot->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_P);
  // _prot_status = abs(_data->status(i));

  _prot_mom_uncorr = _Energy_loss_uncorr_prot->P();
  _prot_theta = _Energy_loss_uncorr_prot->Theta() * 180 / PI;
  // std::cout << "prot ststus " << _data->status(i) << "   prot theta " << _prot_theta << " prot  mom   "
  //           << _prot_mom_uncorr<< std::endl;
  if (abs(_data->status(i)) < 4000) {
    _sectorProt = _data->dc_sec(i);

    // fpro = dppC(_data->px(i), _data->py(i), _data->pz(i), _data->dc_sec(i), 3) + 1;

    if (_prot_theta <= 27) {
      _E_corr_val_prot = -0.00078846 * pow(_prot_mom_uncorr, 5) + 0.0093734 * pow(_prot_mom_uncorr, 4) -
                         0.04277868 * pow(_prot_mom_uncorr, 3) + 0.09421284 * pow(_prot_mom_uncorr, 2) -
                         0.10095842 * (_prot_mom_uncorr) + 0.04567203;
    } else {
      _E_corr_val_prot = -0.0023389 * pow(_prot_mom_uncorr, 5) + 0.02838603 * pow(_prot_mom_uncorr, 4) -
                         0.13214962 * pow(_prot_mom_uncorr, 3) + 0.29609571 * pow(_prot_mom_uncorr, 2) -
                         0.32307424 * (_prot_mom_uncorr) + 0.14742569;
    }

    // _E_corr_val_prot = (-2.38530518e-07 * pow(_prot_theta, 3) + 1.04005173e-05 * pow(_prot_theta, 2) +
    //                     (-1.49765839e-04) * (_prot_theta) + (-1.41952603e-05)) *
    //                        pow(_prot_mom_uncorr, 5) +

    //                    (1.57607943e-06 * pow(_prot_theta, 3) + (-4.53383617e-05) * pow(_prot_theta, 2) +
    //                     3.07855716e-04 * (_prot_theta) + (8.04420446e-03)) *
    //                        pow(_prot_mom_uncorr, 4) +

    //                    ((-1.67283666e-06) * pow(_prot_theta, 3) + (-1.39877228e-04) * pow(_prot_theta, 2) +
    //                     5.19206096e-03 * (_prot_theta)-7.29631272e-02) *
    //                        pow(_prot_mom_uncorr, 3) +

    //                    (-7.85196782e-06 * pow(_prot_theta, 3) + 1.04094323e-03 * pow(_prot_theta, 2) +
    //                     (-2.55883006e-02) * (_prot_theta) + 2.40630846e-01) *
    //                        pow(_prot_mom_uncorr, 2) +

    //                    (2.03637619e-05 * pow(_prot_theta, 3) + (-1.88926745e-03) * pow(_prot_theta, 2) +
    //                     4.26626313e-02 * (_prot_theta)-3.46543011e-01) *
    //                        (_prot_mom_uncorr) +

    //                    ((-1.37003888e-05) * pow(_prot_theta, 3) + (1.15294285e-03) * pow(_prot_theta, 2) +
    //                     (-2.53253135e-02) * (_prot_theta) + 1.93782983e-01);

  } else if (abs(_data->status(i)) >= 4000) {
    // fpro = 1.0;
    // _E_corr_val_prot = 0.01066342 * pow(_prot_mom_uncorr, 2) - 0.05379427 * (_prot_mom_uncorr) + 0.02530928;
    // a x ^ 2 bx c[-9.30990933e-05 1.23584235e-02 - 5.42538215e-01 7.87921215e+00]........................ a x
    // ^
    //     2 bx c[4.17955911e-04 - 5.53676478e-02 2.42642631e+00 - 3.51829220e+01]........................ a x ^
    //     2 bx c[-5.58084320e-04 7.38670367e-02 - 3.23723227e+00 4.69456718e+01]........................ a x ^
    //     2 bx c[2.40014720e-04 - 3.17071405e-02 1.38769727e+00 -2.01072704e+01]........................;

    _E_corr_val_prot = ((-9.30990933e-05) * pow(_prot_theta, 3) + (1.23584235e-02) * pow(_prot_theta, 2) +
                        (-5.42538215e-01) * (_prot_theta) + 7.87921215e+00) *
                           pow(_prot_mom_uncorr, 3) +

                       (4.17955911e-04 * pow(_prot_theta, 3) + (-5.53676478e-02) * pow(_prot_theta, 2) +
                        (2.42642631e+00) * (_prot_theta) + (-3.51829220e+01)) *
                           pow(_prot_mom_uncorr, 2) +

                       ((-5.58084320e-04) * pow(_prot_theta, 3) + (7.38670367e-02) * pow(_prot_theta, 2) +
                        (-3.23723227e+00) * (_prot_theta) + 4.69456718e+01) *
                           (_prot_mom_uncorr) +

                       ((2.40014720e-04) * pow(_prot_theta, 3) + (-3.17071405e-02) * pow(_prot_theta, 2) +
                        (1.38769727e+00 * (_prot_theta)) + (-2.01072704e+01));
  }

  _prot_mom_tmt = _prot_mom_uncorr + _E_corr_val_prot;

  _px_prime_prot_E = _data->px(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
  _py_prime_prot_E = _data->py(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));
  _pz_prime_prot_E = _data->pz(i) * ((_prot_mom_tmt) / (_prot_mom_uncorr));

  // _prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);
  // _mom_corr_prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);

  // Below shows how the corrections are to be applied using the ROOT momentum 4-vector using the above code:
  if (abs(_data->status(i)) < 4000) {
    fpro = dppC(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, _data->dc_sec(i), 3) + 1;
  } else {
    fpro = 1.0;
  }

  // _px_prime_prot_E = _data->px(i) * fpro * ((_prot_mom_tmt) / (_prot_mom_uncorr));
  // _py_prime_prot_E = _data->py(i) * fpro * ((_prot_mom_tmt) / (_prot_mom_uncorr));
  // _pz_prime_prot_E = _data->pz(i) * fpro * ((_prot_mom_tmt) / (_prot_mom_uncorr));
  // _prot->SetXYZM(_px_prime_prot_E, _py_prime_prot_E, _pz_prime_prot_E, MASS_P);

  _prot->SetXYZM(_px_prime_prot_E * fpro, _py_prime_prot_E * fpro, _pz_prime_prot_E * fpro, MASS_P);

  if (_prot->Phi() > 0)
    _prot_phi = _prot->Phi() * 180 / PI;
  else if (_prot->Phi() < 0)
    _prot_phi = (_prot->Phi() + 2 * PI) * 180 / PI;

  // for (size_t t = 0; t < Prot_theta_bins; t++) {
  //   double theta_min = min_prot_theta_values[t];
  //   double theta_max = max_prot_theta_values[t];
  //   if (_prot_theta > theta_min && _prot_theta < theta_max) {
  //     // for experimental dsta
  //     _prot_theta_prime = _prot_theta - prot_theta_corr[t] * alpha_prot_theta_corr;
  //     // //for simulation data
  //     //       _prot_theta_prime = _prot_theta - prot_theta_corr_sim[t] * alpha_prot_theta_corr;

  //     _px_prime_prot_th = _data->px(i) * (sin(DEG2RAD * _prot_theta_prime) / sin(DEG2RAD * _prot_theta));
  //     _py_prime_prot_th = _data->py(i) * (sin(DEG2RAD * _prot_theta_prime) / sin(DEG2RAD * _prot_theta));
  //     _pz_prime_prot_th = _data->pz(i) * (cos(DEG2RAD * _prot_theta_prime) / cos(DEG2RAD * _prot_theta));
  //     _E_prime_prot_th = sqrt(abs(_px_prime_prot_th * _px_prime_prot_th + _py_prime_prot_th *
  //     _py_prime_prot_th +
  //                                 _pz_prime_prot_th * _pz_prime_prot_th));
  //     _mom_corr_prot_th->SetPxPyPzE(_px_prime_prot_th, _py_prime_prot_th, _pz_prime_prot_th,
  //     _E_prime_prot_th);
  //   }
  // }

  //   for (size_t p = 0; p < Prot_phi_bins; p++) {
  //     double phi_min = min_prot_phi_values[p];
  //     double phi_max = max_prot_phi_values[p];
  //     if (_prot_phi > phi_min && _prot_phi < phi_max) {
  //       // for experimental data
  //       _prot_phi_prime = _prot_phi - prot_phi_corr[p] * alpha_prot_phi_corr;

  //       // // for simulation data
  //       // _prot_phi_prime = _prot_phi - prot_phi_corr_sim[p] * alpha_prot_phi_corr;

  //       _px_prime_prot_ph = _mom_corr_prot_th->Px() * (cos(DEG2RAD * _prot_phi_prime) / cos(DEG2RAD *
  //       _prot_phi)); _py_prime_prot_ph = _mom_corr_prot_th->Py() * (sin(DEG2RAD * _prot_phi_prime) /
  //       sin(DEG2RAD * _prot_phi)); _pz_prime_prot_ph = _mom_corr_prot_th->Pz();

  //       _mom_corr_prot_ph->SetXYZM(_px_prime_prot_ph, _py_prime_prot_ph, _pz_prime_prot_ph, MASS_P);
  //     }
  // }

  // _px_prime_prot_mom = _mom_corr_prot_ph->Px() * ((_prot_mom_prime) / (_prot_mom));
  // _py_prime_prot_mom = _mom_corr_prot_ph->Py() * ((_prot_mom_prime) / (_prot_mom));
  // _pz_prime_prot_mom = _mom_corr_prot_ph->Pz() * ((_prot_mom_prime) / (_prot_mom));
  // _mom_corr_prot->SetXYZM(_px_prime_prot_mom, _py_prime_prot_mom, _pz_prime_prot_mom, MASS_P);

  // ///---- For Hadron mom corr  --------///

  _prot_mom = _prot->P();

  if (_prot->Phi() > 0)
    _prot_phi = _prot->Phi() * 180 / PI;
  else if (_prot->Phi() < 0)
    _prot_phi = (_prot->Phi() + 2 * PI) * 180 / PI;

  if (abs(_data->status(i)) < 4000) {
    for (size_t m = 0; m < Prot_mom_bins_FD; m++) {
      double mom_min = min_prot_mom_values_FD[m];
      double mom_max = max_prot_mom_values_FD[m];
      if (_prot_mom > mom_min && _prot_mom < mom_max) {
        if (_prot_theta <= 27) {
          //   //   // For experimental data
          if (_data->dc_sec(i) == 1) _prot_mom_prime = _prot_mom - prot_mom_corr_FD[0][0][m] * alpha_prot_mom_corr_FD[0];
          if (_data->dc_sec(i) == 2) _prot_mom_prime = _prot_mom - prot_mom_corr_FD[0][1][m] * alpha_prot_mom_corr_FD[0];
          if (_data->dc_sec(i) == 3) _prot_mom_prime = _prot_mom - prot_mom_corr_FD[0][2][m] * alpha_prot_mom_corr_FD[0];
          if (_data->dc_sec(i) == 4) _prot_mom_prime = _prot_mom - prot_mom_corr_FD[0][3][m] * alpha_prot_mom_corr_FD[0];
          if (_data->dc_sec(i) == 5) _prot_mom_prime = _prot_mom - prot_mom_corr_FD[0][4][m] * alpha_prot_mom_corr_FD[0];
          if (_data->dc_sec(i) == 6) _prot_mom_prime = _prot_mom - prot_mom_corr_FD[0][5][m] * alpha_prot_mom_corr_FD[0];

        } else {
          if (_data->dc_sec(i) == 1) _prot_mom_prime = _prot_mom - prot_mom_corr_FD[1][0][m] * alpha_prot_mom_corr_FD[0];
          if (_data->dc_sec(i) == 2) _prot_mom_prime = _prot_mom - prot_mom_corr_FD[1][1][m] * alpha_prot_mom_corr_FD[0];
          if (_data->dc_sec(i) == 3) _prot_mom_prime = _prot_mom - prot_mom_corr_FD[1][2][m] * alpha_prot_mom_corr_FD[0];
          if (_data->dc_sec(i) == 4) _prot_mom_prime = _prot_mom - prot_mom_corr_FD[1][3][m] * alpha_prot_mom_corr_FD[0];
          if (_data->dc_sec(i) == 5) _prot_mom_prime = _prot_mom - prot_mom_corr_FD[1][4][m] * alpha_prot_mom_corr_FD[0];
          if (_data->dc_sec(i) == 6) _prot_mom_prime = _prot_mom - prot_mom_corr_FD[1][5][m] * alpha_prot_mom_corr_FD[0];
        }
      }
    }
  } else if (abs(_data->status(i)) >= 4000) {
    for (size_t m = 0; m < Prot_mom_bins_CD; m++) {
      double mom_min = min_prot_mom_values_CD[m];
      double mom_max = max_prot_mom_values_CD[m];
      if (_prot_mom > mom_min && _prot_mom < mom_max) {
        //   // For experimental data
        if (_prot_phi > 270 || _prot_phi <= 30)
          _prot_mom_prime = _prot_mom - prot_mom_corr_CD[0][m] * alpha_prot_mom_corr_CD[0];
        else if (_prot_phi > 30 && _prot_phi <= 150) {
          if (_pim_mom < 2.0)
            _prot_mom_prime = _prot_mom - prot_mom_corr_CD[1][m] * alpha_prot_mom_corr_CD[1];
            else
              _prot_mom_prime = _prot_mom - prot_mom_corr_CD[1][m] * alpha_prot_mom_corr_CD[2];
        } else if (_prot_phi > 150 && _prot_phi <= 270) {
          if (_pim_mom < 2.0) _prot_mom_prime = _prot_mom - prot_mom_corr_CD[2][m] * alpha_prot_mom_corr_CD[3];
          else
            _prot_mom_prime = _prot_mom - prot_mom_corr_CD[2][m] * alpha_prot_mom_corr_CD[4];}
        }
    }
  }

  _px_prime_prot_mom = _prot->Px() * ((_prot_mom_prime) / (_prot_mom));
  _py_prime_prot_mom = _prot->Py() * ((_prot_mom_prime) / (_prot_mom));
  _pz_prime_prot_mom = _prot->Pz() * ((_prot_mom_prime) / (_prot_mom));
  _mom_corr_prot->SetXYZM(_px_prime_prot_mom, _py_prime_prot_mom, _pz_prime_prot_mom, MASS_P);

  // // 2nd iteration
  //   _prot_mom_2nd = _prot_mom_prime;
  //   if (abs(_data->status(i)) < 4000) {
  //     // _prot_mom_2nd = _prot_mom_prime;

  //     for (size_t m = 0; m < Prot_mom_bins_FD; m++) {
  //       double mom_min = min_prot_mom_values_FD[m];
  //       double mom_max = max_prot_mom_values_FD[m];
  //       if (_prot_mom_2nd > mom_min && _prot_mom_2nd < mom_max) {
  //         if (_prot_theta <= 27) {
  //           //   //   // For experimental data
  //           _prot_mom_prime_2nd = _prot_mom_2nd - prot_mom_corr_FD_2nd[0][m] * alpha_prot_mom_corr_2nd[1];
  //         } else {
  //           _prot_mom_prime_2nd = _prot_mom_2nd - prot_mom_corr_FD_2nd[1][m] * alpha_prot_mom_corr_2nd[2];
  //         }
  //       }
  //     }
  //   }
  //   else if (abs(_data->status(i)) >= 4000) {
  //     for (size_t m = 0; m < Prot_mom_bins_CD; m++) {
  //       double mom_min = min_prot_mom_values_CD[m];
  //       double mom_max = max_prot_mom_values_CD[m];
  //       if (_prot_mom_2nd > mom_min && _prot_mom_2nd < mom_max) {
  //         //   // For experimental data
  //         // _prot_mom_prime = _prot_mom - prot_mom_corr_CD[m] * alpha_prot_mom_corr;
  //         _prot_mom_prime_2nd = _prot_mom_2nd - prot_mom_corr_CD_2nd[m] * alpha_prot_mom_corr_2nd[0];
  //       }
  //     }
  //   }

  //   _px_prime_prot_mom = _prot->Px() * ((_prot_mom_prime_2nd) / (_prot_mom));
  //   _py_prime_prot_mom = _prot->Py() * ((_prot_mom_prime_2nd) / (_prot_mom));
  //   _pz_prime_prot_mom = _prot->Pz() * ((_prot_mom_prime_2nd) / (_prot_mom));
  //   _mom_corr_prot->SetXYZM(_px_prime_prot_mom, _py_prime_prot_mom, _pz_prime_prot_mom, MASS_P);
}

// bool Reaction::ctof_prot() {
//   bool _prot_ctof = true;
//   _prot_ctof &= (4000 <= _prot_status && _prot_status < 6000);
//   return _prot_ctof;
// }
// bool Reaction::ftof_prot() {
//   bool _prot_ftof = true;
//   _prot_ftof &= (2000 <= _prot_status && _prot_status < 4000);
//   return _prot_ftof;
// }
void Reaction::SetPip(int i) {
  _numPip++;
  _numPos++;
  _hasPip = true;

  // _pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
  //   // _pip_status = abs(_data->status(i));
  _Energy_loss_uncorr_pip->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIP);
  _pip_mom_uncorr = _Energy_loss_uncorr_pip->P();
  _pip_theta = _Energy_loss_uncorr_pip->Theta() * 180 / PI;
  // std::cout << "pip ststus " << _data->status(i) << "   pip theta " << _pip_theta << " pip  mom   "
  //           << _pip_mom_uncorr<< std::endl;
  if (abs(_data->status(i)) < 4000) {
    _sectorPip = _data->dc_sec(i);
    // fpip = dppC(_data->px(i), _data->py(i), _data->pz(i), _data->dc_sec(i), 1) + 1;

    if (_pip_theta <= 27) {
      _E_corr_val_pip = 9.21970527e-05 * pow(_pip_mom_uncorr, 3) - 3.70500143e-04 * pow(_pip_mom_uncorr, 2) +
                        2.78880101e-04 * (_pip_mom_uncorr) + 2.66040566e-03;

    } else {
      _E_corr_val_pip = -0.00010482 * pow(_pip_mom_uncorr, 3) + 0.00080463 * pow(_pip_mom_uncorr, 2) -
                        0.0022871 * (_pip_mom_uncorr) + 0.00831496;
    }
  } else if (abs(_data->status(i)) >= 4000) {
    // fpip = 1.0;

    // _E_corr_val_pip = -0.00631413  * pow(_pip_mom_uncorr, 5) + 0.04713584  * pow(_pip_mom_uncorr, 4) -
    //                   0.12554256 * pow(_pip_mom_uncorr, 3) + 0.15622077 * pow(_pip_mom_uncorr, 2) -
    //                   0.11467851 * (_pip_mom_uncorr) + 0.01917004;
    _E_corr_val_pip = (-6.50509539e-07 * pow(_pip_theta, 3) + 1.31547371e-04 * pow(_pip_theta, 2) +
                       (-7.99024673e-03) * (_pip_theta) + 1.60563630e-01) *
                          pow(_pip_mom_uncorr, 3) +

                      (2.48202211e-06 * pow(_pip_theta, 3) + (-5.15757241e-04) * pow(_pip_theta, 2) +
                       3.19833135e-02 * (_pip_theta) + (-6.53476057e-01)) *
                          pow(_pip_mom_uncorr, 2) +

                      (-2.71923009e-06 * pow(_pip_theta, 3) + 5.80375203e-04 * pow(_pip_theta, 2) +
                       (-3.75941898e-02) * (_pip_theta) + 7.80443724e-01) *
                          (_pip_mom_uncorr) +

                      4.62456800e-07 * pow(_pip_theta, 3) + (-1.08401698e-04) * pow(_pip_theta, 2) +
                      8.09261138e-03 * (_pip_theta)-2.05315604e-01;

    // _E_corr_val_pip =  -0.00279293 * pow(_pip_mom_uncorr, 3) + 0.0206818 * pow(_pip_mom_uncorr, 2) -
    //                   0.05257802 * pow(_pip_mom_uncorr, 2) + 0.00996933;
  }
  _pip_mom_tmt = _pip_mom_uncorr + _E_corr_val_pip;

  _px_prime_pip_E = _data->px(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
  _py_prime_pip_E = _data->py(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));
  _pz_prime_pip_E = _data->pz(i) * ((_pip_mom_tmt) / (_pip_mom_uncorr));

  // _pip->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);
  // _mom_corr_pip->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);

  if (abs(_data->status(i)) < 4000) {
    fpip = dppC(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, _data->dc_sec(i), 1) + 1;
  } else {
    fpip = 1.0;
  }
  _pip->SetXYZM(_px_prime_pip_E * fpip, _py_prime_pip_E * fpip, _pz_prime_pip_E * fpip, MASS_PIP);

  // if (abs(_data->status(i)) < 4000) {

  //     _E_corr_val_pip_th = 0.00000000;

  // } else if (abs(_data->status(i)) >= 4000) {

  //   _E_corr_val_pip_th = -7.08389160e-11 * pow(_pip_theta, 5) + 3.75704402e-08 * pow(_pip_theta, 4) -
  //                        7.26740433e-06 * pow(_pip_theta, 3) + 6.45415606e-04 * pow(_pip_theta, 2) -
  //                        2.60057363e-02 * (_pip_theta) + 3.78387868e-01;
  // }
  // _pip_mom_tmt2 = _pip_mom_tmt + _E_corr_val_pip_th;  // theta iteration

  // _px_prime_pip_E_tmt = _data->px(i) * ((_pip_mom_tmt2) / (_pip_mom_uncorr));
  // _py_prime_pip_E_tmt = _data->py(i) * ((_pip_mom_tmt2) / (_pip_mom_uncorr));
  // _pz_prime_pip_E_tmt = _data->pz(i) * ((_pip_mom_tmt2) / (_pip_mom_uncorr));

  // _pip->SetXYZM(_px_prime_pip_E_tmt, _py_prime_pip_E_tmt, _pz_prime_pip_E_tmt, MASS_PIP);

  // second iterations

  // _pip_mom_tmt2 = _pip_tmt->P();  // for second iteration

  // // let's do second iteration for cd pip
  // if (abs(_data->status(i)) < 4000) {
  //   _E_corr_val_pip2 = 0.0;
  // } else if (abs(_data->status(i)) >= 4000) {
  //   _E_corr_val_pip2 = -0.00125164 * pow(_pip_mom_tmt2, 5) + 0.01272027 * pow(_pip_mom_tmt2, 4) -
  //                      0.04457356 * pow(_pip_mom_tmt2, 3) + 0.06272048 * pow(_pip_mom_tmt2, 2) -
  //                      0.03798534 * (_pip_mom_tmt2)-0.00716495;
  // }
  // _pip_mom = _pip_mom_tmt2 + _E_corr_val_pip2;
  // _px_prime_pip_E = _data->px(i) * ((_pip_mom) / (_pip_mom_tmt2));
  // _py_prime_pip_E = _data->py(i) * ((_pip_mom) / (_pip_mom_tmt2));
  // _pz_prime_pip_E = _data->pz(i) * ((_pip_mom) / (_pip_mom_tmt2));

  // _pip->SetXYZM(_px_prime_pip_E, _py_prime_pip_E, _pz_prime_pip_E, MASS_PIP);

  // //   _pip_mom = _pip->P();
  // //   _pip_theta = _pip->Theta() * 180 / PI;

  //   if (_pip->Phi() > 0)
  //     _pip_phi = _pip->Phi() * 180 / PI;
  //   else if (_pip->Phi() < 0)
  //     _pip_phi = (_pip->Phi() + 2 * PI) * 180 / PI;

  //   for (size_t t = 0; t < Pip_theta_bins; t++) {
  //     double theta_min = min_pip_theta_values[t];
  //     double theta_max = max_pip_theta_values[t];
  //     if (_pip_theta > theta_min && _pip_theta < theta_max) {
  //       // For experimental data
  //       _pip_theta_prime = _pip_theta - pip_theta_corr[t] * alpha_pip_theta_corr;
  //       // // For simulation data
  //       // _pip_theta_prime = _pip_theta - pip_theta_corr_sim[t] * alpha_pip_theta_corr;

  //       _px_prime_pip_th = _data->px(i) * (sin(DEG2RAD * _pip_theta_prime) / sin(DEG2RAD * _pip_theta));
  //       _py_prime_pip_th = _data->py(i) * (sin(DEG2RAD * _pip_theta_prime) / sin(DEG2RAD * _pip_theta));
  //       _pz_prime_pip_th = _data->pz(i) * (cos(DEG2RAD * _pip_theta_prime) / cos(DEG2RAD * _pip_theta));
  //       _E_prime_pip_th = sqrt(abs(_px_prime_pip_th * _px_prime_pip_th + _py_prime_pip_th * _py_prime_pip_th
  //       +
  //                                  _pz_prime_pip_th * _pz_prime_pip_th));
  //       _mom_corr_pip_th->SetPxPyPzE(_px_prime_pip_th, _py_prime_pip_th, _pz_prime_pip_th, _E_prime_pip_th);
  //     }
  //   }

  //   for (size_t p = 0; p < Pip_phi_bins; p++) {
  //     double phi_min = min_pip_phi_values[p];
  //     double phi_max = max_pip_phi_values[p];
  //     if (_pip_phi > phi_min && _pip_phi < phi_max) {
  //       //For experimantal data
  //       _pip_phi_prime = _pip_phi - pip_phi_corr[p] * alpha_pip_phi_corr;
  // // For simulations data
  //       // _pip_phi_prime = _pip_phi - pip_phi_corr_sim[p] * alpha_pip_phi_corr;

  //       _px_prime_pip_ph = _mom_corr_pip_th->Px() * (cos(DEG2RAD * _pip_phi_prime) / cos(DEG2RAD *
  //       _pip_phi)); _py_prime_pip_ph = _mom_corr_pip_th->Py() * (sin(DEG2RAD * _pip_phi_prime) / sin(DEG2RAD
  //       * _pip_phi)); _pz_prime_pip_ph = _mom_corr_pip_th->Pz();

  //       _mom_corr_pip_ph->SetXYZM(_px_prime_pip_ph, _py_prime_pip_ph, _pz_prime_pip_ph, MASS_PIP);
  //     }
  //   }

  // for (size_t m = 0; m < Pip_mom_bins; m++) {
  //   double mom_min = min_pip_mom_values[m];
  //   double mom_max = max_pip_mom_values[m];
  //   if (_pip_mom > mom_min && _pip_mom < mom_max) {
  //     // For experimantal data
  //     _pip_mom_prime = _pip_mom - pip_mom_corr[m] * alpha_pip_mom_corr;
  //     // // For simulation data
  //     // _pip_mom_prime = _pip_mom - pip_mom_corr_sim[m] * alpha_pip_mom_corr;

  //     _px_prime_pip_mom = _mom_corr_pip_ph->Px() * ((_pip_mom_prime) / (_pip_mom));
  //     _py_prime_pip_mom = _mom_corr_pip_ph->Py() * ((_pip_mom_prime) / (_pip_mom));
  //     _pz_prime_pip_mom = _mom_corr_pip_ph->Pz() * ((_pip_mom_prime) / (_pip_mom));
  //     _mom_corr_pip->SetXYZM(_px_prime_pip_mom, _py_prime_pip_mom, _pz_prime_pip_mom, MASS_PIP);
  //   }
  // }

  // ///---- For Hadron mom corr  --------///
  _pip_mom = _pip->P();
  if (_pip->Phi() > 0)
    _pip_phi = _pip->Phi() * 180 / PI;
  else if (_pip->Phi() < 0)
    _pip_phi = (_pip->Phi() + 2 * PI) * 180 / PI;

  if (abs(_data->status(i)) < 4000) {
    for (size_t m = 0; m < Pip_mom_bins_FD; m++) {
      if (_pip_theta <= 27) {
        double mom_min = min_pip_mom_values_FD[0][m];
        double mom_max = max_pip_mom_values_FD[0][m];
        if (_pip_mom > mom_min && _pip_mom < mom_max) {
          //   //   // For experimental data
          if (_data->dc_sec(i) == 1) _pip_mom_prime = _pip_mom - pip_mom_corr_FD[0][0][m] * alpha_pip_mom_corr_FD[0];
          if (_data->dc_sec(i) == 2) _pip_mom_prime = _pip_mom - pip_mom_corr_FD[0][1][m] * alpha_pip_mom_corr_FD[0];
          if (_data->dc_sec(i) == 3) _pip_mom_prime = _pip_mom - pip_mom_corr_FD[0][2][m] * alpha_pip_mom_corr_FD[0];
          if (_data->dc_sec(i) == 4) _pip_mom_prime = _pip_mom - pip_mom_corr_FD[0][3][m] * alpha_pip_mom_corr_FD[0];
          if (_data->dc_sec(i) == 5) _pip_mom_prime = _pip_mom - pip_mom_corr_FD[0][4][m] * alpha_pip_mom_corr_FD[0];
          if (_data->dc_sec(i) == 6) _pip_mom_prime = _pip_mom - pip_mom_corr_FD[0][5][m] * alpha_pip_mom_corr_FD[0];
        }
      } else {
        double mom_min = min_pip_mom_values_FD[1][m];
        double mom_max = max_pip_mom_values_FD[1][m];
        if (_pip_mom > mom_min && _pip_mom < mom_max) {
          if (_data->dc_sec(i) == 1) _pip_mom_prime = _pip_mom - pip_mom_corr_FD[1][0][m] * alpha_pip_mom_corr_FD[1];
          if (_data->dc_sec(i) == 2) _pip_mom_prime = _pip_mom - pip_mom_corr_FD[1][1][m] * alpha_pip_mom_corr_FD[1];
          if (_data->dc_sec(i) == 3) _pip_mom_prime = _pip_mom - pip_mom_corr_FD[1][2][m] * alpha_pip_mom_corr_FD[1];
          if (_data->dc_sec(i) == 4) _pip_mom_prime = _pip_mom - pip_mom_corr_FD[1][3][m] * alpha_pip_mom_corr_FD[1];
          if (_data->dc_sec(i) == 5) _pip_mom_prime = _pip_mom - pip_mom_corr_FD[1][4][m] * alpha_pip_mom_corr_FD[1];
          if (_data->dc_sec(i) == 6) _pip_mom_prime = _pip_mom - pip_mom_corr_FD[1][5][m] * alpha_pip_mom_corr_FD[1];
        }
      }
    }
  } else if (abs(_data->status(i)) >= 4000) {
    for (size_t m = 0; m < Pip_mom_bins_CD; m++) {
      double mom_min = min_pip_mom_values_CD[m];
      double mom_max = max_pip_mom_values_CD[m];
      if (_pip_mom > mom_min && _pip_mom < mom_max) {
        //   // For experimental data
        if (_pip_phi > 270 || _pip_phi <= 30)
          _pip_mom_prime = _pip_mom - pip_mom_corr_CD[0][m] * alpha_pip_mom_corr_CD[0];
        else if (_pip_phi > 30 && _pip_phi <= 150)
          _pip_mom_prime = _pip_mom - pip_mom_corr_CD[1][m] * alpha_pip_mom_corr_CD[1];
        else if (_pip_phi > 150 && _pip_phi <= 270)
          _pip_mom_prime = _pip_mom - pip_mom_corr_CD[2][m] * alpha_pip_mom_corr_CD[2];
      }
    }
  }

  _px_prime_pip_mom = _pip->Px() * ((_pip_mom_prime) / (_pip_mom));
  _py_prime_pip_mom = _pip->Py() * ((_pip_mom_prime) / (_pip_mom));
  _pz_prime_pip_mom = _pip->Pz() * ((_pip_mom_prime) / (_pip_mom));
  _mom_corr_pip->SetXYZM(_px_prime_pip_mom, _py_prime_pip_mom, _pz_prime_pip_mom, MASS_PIP);

  // // 2nd iteration
  // _pip_mom_2nd = _pip_mom_prime;
  // if (abs(_data->status(i)) < 4000) {
  //    for (size_t m = 0; m < Pip_mom_bins_FD; m++) {
  //     double mom_min = min_pip_mom_values_FD[m];
  //     double mom_max = max_pip_mom_values_FD[m];
  //     if (_pip_mom_2nd > mom_min && _pip_mom_2nd < mom_max) {
  //       if (_pip_theta <= 27) {
  //           //   // For experimental data
  //         _pip_mom_prime_2nd = _pip_mom_2nd - pip_mom_corr_FD_2nd[0][m] * alpha_pip_mom_corr_2nd[1];
  //       } else  {
  //          _pip_mom_prime_2nd = _pip_mom_2nd - pip_mom_corr_FD_2nd[1][m] * alpha_pip_mom_corr_2nd[2]; }
  //       }
  //   }
  // }
  // else if (abs(_data->status(i)) >= 4000) {

  //   for (size_t m = 0; m < Pip_mom_bins_CD; m++) {
  //     double mom_min = min_pip_mom_values_CD[m];
  //     double mom_max = max_pip_mom_values_CD[m];
  //     if (_pip_mom_2nd > mom_min && _pip_mom_2nd < mom_max) {
  //       //   // For experimental data
  //       // _pip_mom_prime = _pip_mom - pip_mom_corr_CD[m] * alpha_pip_mom_corr;
  //       _pip_mom_prime_2nd = _pip_mom_2nd - pip_mom_corr_CD_2nd[m] * alpha_pip_mom_corr_2nd[0];
  //     }
  //   }
  // }

  // _px_prime_pip_mom = _pip->Px() * ((_pip_mom_prime_2nd) / (_pip_mom));
  // _py_prime_pip_mom = _pip->Py() * ((_pip_mom_prime_2nd) / (_pip_mom));
  // _pz_prime_pip_mom = _pip->Pz() * ((_pip_mom_prime_2nd) / (_pip_mom));
  // _mom_corr_pip->SetXYZM(_px_prime_pip_mom, _py_prime_pip_mom, _pz_prime_pip_mom, MASS_PIP);
}

// bool Reaction::ctof_pip() {
//   bool _pip_ctof = true;
//   _pip_ctof &= (4000 <= _pip_status && _pip_status < 6000);
//   return _pip_ctof;
// }
// bool Reaction::ftof_pip() {
//   bool _pip_ftof = true;
//   _pip_ftof &= (2000 <= _pip_status && _pip_status < 4000);
//   return _pip_ftof;
// }

void Reaction::SetPim(int i) {
  _numPim++;
  _numNeg++;
  _hasPim = true;

  // _pim->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_PIM);
  // // // _pim_status = abs(_data->status(i));
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
    // fpim = 1.0;

    // _E_corr_val_pim = (0.02153442) * pow(_pim_mom_uncorr, 5) -
    //                   (0.13271424) * pow(_pim_mom_uncorr, 4) +
    //                   (0.27140262) * pow(_pim_mom_uncorr, 3) -
    //                   (0.23266059) * pow(_pim_mom_uncorr, 2) +
    //                   (0.04031421) * (_pim_mom_uncorr) + 0.0036634;

    _E_corr_val_pim = (-4.94426765e-07 * pow(_pim_theta, 3) + 9.85729368e-05 * pow(_pim_theta, 2) +
                       (-5.85778699e-03) * (_pim_theta) + 1.17447168e-01) *
                          pow(_pim_mom_uncorr, 3) +

                      (1.75953956e-06 * pow(_pim_theta, 3) + (-3.63382515e-04) * pow(_pim_theta, 2) +
                       2.21447425e-02 * (_pim_theta) + (-4.54844509e-01)) *
                          pow(_pim_mom_uncorr, 2) +

                      (-1.90446515e-06 * pow(_pim_theta, 3) + 4.08768480e-04 * pow(_pim_theta, 2) +
                       (-2.65277055e-02) * (_pim_theta) + 5.57286393e-01) *
                          (_pim_mom_uncorr) +

                      2.05653097e-07 * pow(_pim_theta, 3) + (-5.44018546e-05) * pow(_pim_theta, 2) +
                      4.61561853e-03 * (_pim_theta)-1.35303212e-01;
  }

  // _pim_mom = _pim_mom_uncorr + _E_corr_val_pim; // first iteration

  _pim_mom_tmt = _pim_mom_uncorr + _E_corr_val_pim;  // first iteration

  _px_prime_pim_E = _data->px(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
  _py_prime_pim_E = _data->py(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));
  _pz_prime_pim_E = _data->pz(i) * ((_pim_mom_tmt) / (_pim_mom_uncorr));

  // _pim->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, MASS_PIM);  // energy loss corrected
  // _mom_corr_pim->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, MASS_PIM);

  if (abs(_data->status(i)) < 4000) {
    // std::cout << "pim sec is " << _data->dc_sec(i) <<std::endl;
    fpim = dppC(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, _data->dc_sec(i), 2) + 1;
  } else {
    fpim = 1.0;
  }
  _pim->SetXYZM(_px_prime_pim_E * fpim, _py_prime_pim_E * fpim, _pz_prime_pim_E * fpim, MASS_PIM);

  // if (abs(_data->status(i)) < 4000) {

  //     _E_corr_val_pim_th = 0.00000000;

  // } else if (abs(_data->status(i)) >= 4000) {

  //   // -2.07141609e-10 8.81758359e-08 - 1.46534798e-05 1.17681655e-03 - 4.50634123e-02 6.54748237e-01;
  //   _E_corr_val_pim_th = (-2.07141609e-10) * pow(_pim_theta, 5) + (8.81758359e-08) * pow(_pim_theta, 4) +
  //                        (-1.46534798e-05) * pow(_pim_theta, 3) + (1.17681655e-03) * pow(_pim_theta, 2) +
  //                        (-4.50634123e-02) * (_pim_theta) + 6.54748237e-01;
  // }

  //     _pim_mom_tmt2 = _pim_mom_tmt + _E_corr_val_pim_th;

  // _px_prime_pim_E_tmt = _data->px(i) * ((_pim_mom_tmt2) / (_pim_mom_uncorr));
  // _py_prime_pim_E_tmt = _data->py(i) * ((_pim_mom_tmt2) / (_pim_mom_uncorr));
  // _pz_prime_pim_E_tmt = _data->pz(i) * ((_pim_mom_tmt2) / (_pim_mom_uncorr));

  // _pim->SetXYZM(_px_prime_pim_E_tmt, _py_prime_pim_E_tmt, _pz_prime_pim_E_tmt, MASS_PIM);

  // std::cout << "_E_corr_val_pim " << _E_corr_val_pim << "  _E_corr_val_pim_th " << _E_corr_val_pim_th
  //           << "   pim mom tmt  " << _pim_mom_tmt << "   pim mom tmt2  " << _pim_mom_tmt2 << " diff "
  //           << _pim_mom_tmt - _pim_mom_tmt2 << std::endl;

  // _pim_tmt->SetXYZM(_px_prime_pim_E_tmt, _py_prime_pim_E_tmt, _pz_prime_pim_E_tmt, MASS_PIM);

  // _pim_mom_tmt2 = _pim_tmt->P();  // for second iteration

  // // std::cout << " diff  " << _pim_tmt->P() - _pim_mom_tmt << std::endl;
  // // let's do second iteration for cd pim
  // if (abs(_data->status(i)) < 4000) {
  //   _E_corr_val_pim2 = 0.0;
  // } else if (abs(_data->status(i)) >= 4000) {
  //   _E_corr_val_pim2 = 0.07604229 * pow(_pim_mom_tmt2, 7) - 0.69056865 * pow(_pim_mom_tmt2, 6) +
  //                      2.42244641 * pow(_pim_mom_tmt2, 5) - 4.26630462 * pow(_pim_mom_tmt2, 4) +
  //                      4.07033382 * pow(_pim_mom_tmt2, 3) - 2.09075715 * pow(_pim_mom_tmt2, 2) +
  //                      0.52748137 * (_pim_mom_tmt2)-0.04274812;
  // }
  // _pim_mom = _pim_mom_tmt2 + _E_corr_val_pim2;
  // _px_prime_pim_E = _data->px(i) * ((_pim_mom) / (_pim_mom_tmt2));
  // _py_prime_pim_E = _data->py(i) * ((_pim_mom) / (_pim_mom_tmt2));
  // _pz_prime_pim_E = _data->pz(i) * ((_pim_mom) / (_pim_mom_tmt2));

  // _pim->SetXYZM(_px_prime_pim_E, _py_prime_pim_E, _pz_prime_pim_E, MASS_PIM);

  // from here it is for mom theta and phi corrections 1-D

  // // _pim_mom = _pim->P();
  // // _pim_theta = _pim->Theta() * 180 / PI;

  // if (_pim->Phi() > 0)
  //   _pim_phi = _pim->Phi()* 180 / PI;
  // else if (_pim->Phi() < 0)
  //   _pim_phi = (_pim->Phi() + 2 * PI)* 180 / PI;

  // for (size_t t = 0; t < Pim_theta_bins; t++) {
  //   double theta_min = min_pim_theta_values[t];
  //   double theta_max = max_pim_theta_values[t];
  //   if (_pim_theta > theta_min && _pim_theta < theta_max) {
  //     //For experimental data
  //     _pim_theta_prime = _pim_theta - pim_theta_corr[t] * alpha_pim_theta_corr;

  //     // // For simulations data
  //     // _pim_theta_prime = _pim_theta - pim_theta_corr_sim[t] * alpha_pim_theta_corr;

  //     _px_prime_pim_th = _data->px(i) * (sin(DEG2RAD * _pim_theta_prime) / sin(DEG2RAD * _pim_theta));
  //     _py_prime_pim_th = _data->py(i) * (sin(DEG2RAD * _pim_theta_prime) / sin(DEG2RAD * _pim_theta));
  //     _pz_prime_pim_th = _data->pz(i) * (cos(DEG2RAD * _pim_theta_prime) / cos(DEG2RAD * _pim_theta));
  //     _E_prime_pim_th =
  //         sqrt(abs(_px_prime_pim_th * _px_prime_pim_th + _py_prime_pim_th* _py_prime_pim_th +
  //         _pz_prime_pim_th * _pz_prime_pim_th));
  //     _mom_corr_pim_th->SetPxPyPzE(_px_prime_pim_th, _py_prime_pim_th, _pz_prime_pim_th, _E_prime_pim_th);
  //   }
  // }

  // for (size_t p = 0; p < Pim_phi_bins; p++) {
  //   double phi_min = min_pim_phi_values[p];
  //   double phi_max = max_pim_phi_values[p];
  //   if (_pim_phi > phi_min && _pim_phi < phi_max) {
  //     // For experimantal data
  //     _pim_phi_prime = _pim_phi - pim_phi_corr[p] * alpha_pim_phi_corr;

  //     // // For simulations data
  //     // _pim_phi_prime = _pim_phi - pim_phi_corr_sim[p] * alpha_pim_phi_corr;

  //     _px_prime_pim_ph = _mom_corr_pim_th->Px() * (cos(DEG2RAD * _pim_phi_prime) / cos(DEG2RAD * _pim_phi));
  //     _py_prime_pim_ph = _mom_corr_pim_th->Py() * (sin(DEG2RAD * _pim_phi_prime) / sin(DEG2RAD * _pim_phi));
  //     _pz_prime_pim_ph = _mom_corr_pim_th->Pz();

  //     _mom_corr_pim_ph->SetXYZM(_px_prime_pim_ph, _py_prime_pim_ph, _pz_prime_pim_ph, MASS_PIM);
  //   }
  // }

  // for (size_t m = 0; m < Pim_mom_bins; m++) {
  //   double mom_min = min_pim_mom_values[m];
  //   double mom_max = max_pim_mom_values[m];
  //   if (_pim_mom > mom_min && _pim_mom < mom_max) {
  //     // For experimantal data
  //     _pim_mom_prime = _pim_mom - pim_mom_corr[m] * alpha_pim_mom_corr;

  //     // // For simulations data
  //     // _pim_mom_prime = _pim_mom - pim_mom_corr_sim[m] * alpha_pim_mom_corr;

  //     _px_prime_pim_mom = _mom_corr_pim_ph->Px() * ((_pim_mom_prime) / (_pim_mom));
  //     _py_prime_pim_mom = _mom_corr_pim_ph->Py() * ((_pim_mom_prime) / (_pim_mom));
  //     _pz_prime_pim_mom = _mom_corr_pim_ph->Pz() * ((_pim_mom_prime) / (_pim_mom));
  //     _mom_corr_pim->SetXYZM(_px_prime_pim_mom, _py_prime_pim_mom, _pz_prime_pim_mom, MASS_PIM);
  //   }
  // }

  // // // Now we are just applying momentum corrections only because we dont have much deviation in theta and
  // phi mes and miss

  // ///---- For Hadron mom corr  --------///

  _pim_mom = _pim->P();
  if (_pim->Phi() > 0)
    _pim_phi = _pim->Phi() * 180 / PI;
  else if (_pim->Phi() < 0)
    _pim_phi = (_pim->Phi() + 2 * PI) * 180 / PI;

  if (abs(_data->status(i)) < 4000) {
    for (size_t m = 0; m < Pim_mom_bins_FD; m++) {
      if (_pim_theta <= 27) {
        double mom_min = min_pim_mom_values_FD[0][m];
        double mom_max = max_pim_mom_values_FD[0][m];
        if (_pim_mom > mom_min && _pim_mom < mom_max) {
          //   //   // For experimental data
          if (_data->dc_sec(i) == 1) _pim_mom_prime = _pim_mom - pim_mom_corr_FD[0][0][m] * alpha_pim_mom_corr_FD[0];
          if (_data->dc_sec(i) == 2) _pim_mom_prime = _pim_mom - pim_mom_corr_FD[0][1][m] * alpha_pim_mom_corr_FD[0];
          if (_data->dc_sec(i) == 3) _pim_mom_prime = _pim_mom - pim_mom_corr_FD[0][2][m] * alpha_pim_mom_corr_FD[0];
          if (_data->dc_sec(i) == 4) _pim_mom_prime = _pim_mom - pim_mom_corr_FD[0][3][m] * alpha_pim_mom_corr_FD[0];
          if (_data->dc_sec(i) == 5) _pim_mom_prime = _pim_mom - pim_mom_corr_FD[0][4][m] * alpha_pim_mom_corr_FD[0];
          if (_data->dc_sec(i) == 6) _pim_mom_prime = _pim_mom - pim_mom_corr_FD[0][5][m] * alpha_pim_mom_corr_FD[0];
        }
      } else {
        double mom_min = min_pim_mom_values_FD[1][m];
        double mom_max = max_pim_mom_values_FD[1][m];
        if (_pim_mom > mom_min && _pim_mom < mom_max) {
          if (_data->dc_sec(i) == 1) _pim_mom_prime = _pim_mom - pim_mom_corr_FD[1][0][m] * alpha_pim_mom_corr_FD[1];
          if (_data->dc_sec(i) == 2) _pim_mom_prime = _pim_mom - pim_mom_corr_FD[1][1][m] * alpha_pim_mom_corr_FD[1];
          if (_data->dc_sec(i) == 3) _pim_mom_prime = _pim_mom - pim_mom_corr_FD[1][2][m] * alpha_pim_mom_corr_FD[1];
          if (_data->dc_sec(i) == 4) _pim_mom_prime = _pim_mom - pim_mom_corr_FD[1][3][m] * alpha_pim_mom_corr_FD[1];
          if (_data->dc_sec(i) == 5) _pim_mom_prime = _pim_mom - pim_mom_corr_FD[1][4][m] * alpha_pim_mom_corr_FD[1];
          if (_data->dc_sec(i) == 6) _pim_mom_prime = _pim_mom - pim_mom_corr_FD[1][5][m] * alpha_pim_mom_corr_FD[1];
        }
      }
    }
  }
    else if (abs(_data->status(i)) >= 4000) {
      for (size_t m = 0; m < Pim_mom_bins_CD; m++) {
        double mom_min = min_pim_mom_values_CD[m];
        double mom_max = max_pim_mom_values_CD[m];
        if (_pim_mom > mom_min && _pim_mom < mom_max) {
          //   // For experimental data
          if (_pim_phi > 270 || _pim_phi <= 30){
            if (_pim_mom < 1.1) _pim_mom_prime = _pim_mom - pim_mom_corr_CD[0][m] * alpha_pim_mom_corr_CD[0];
            else if (_pim_mom >= 1.1) _pim_mom_prime = _pim_mom - pim_mom_corr_CD[0][m] * alpha_pim_mom_corr_CD[1];
          } else if (_pim_phi > 30 && _pim_phi <= 150)
            _pim_mom_prime = _pim_mom - pim_mom_corr_CD[1][m] * alpha_pim_mom_corr_CD[2];
          else if (_pim_phi > 150 && _pim_phi <= 270)
            _pim_mom_prime = _pim_mom - pim_mom_corr_CD[2][m] * alpha_pim_mom_corr_CD[3];
        }
      }
    }

    _px_prime_pim_mom = _pim->Px() * ((_pim_mom_prime) / (_pim_mom));
    _py_prime_pim_mom = _pim->Py() * ((_pim_mom_prime) / (_pim_mom));
    _pz_prime_pim_mom = _pim->Pz() * ((_pim_mom_prime) / (_pim_mom));
    _mom_corr_pim->SetXYZM(_px_prime_pim_mom, _py_prime_pim_mom, _pz_prime_pim_mom, MASS_PIM);

    // // 2nd iteration
    // _pim_mom_2nd = _pim_mom_prime;
    // if (abs(_data->status(i)) < 4000) {
    //   for (size_t m = 0; m < Pim_mom_bins_FD; m++) {
    //     double mom_min = min_pim_mom_values_FD[m];
    //     double mom_max = max_pim_mom_values_FD[m];
    //     if (_pim_mom_2nd > mom_min && _pim_mom_2nd < mom_max) {
    //       if (_pim_theta <= 27) {
    //         //   //   // For experimental data
    //         _pim_mom_prime_2nd = _pim_mom_2nd - pim_mom_corr_FD_2nd[0][m] * alpha_pim_mom_corr_2nd[1];
    //       } else {
    //         _pim_mom_prime_2nd = _pim_mom_2nd - pim_mom_corr_FD_2nd[1][m] * alpha_pim_mom_corr_2nd[2];
    //       }
    //     }
    //   }
    // } else if (abs(_data->status(i)) >= 4000) {
    //   for (size_t m = 0; m < Pim_mom_bins_CD; m++) {
    //     double mom_min = min_pim_mom_values_CD[m];
    //     double mom_max = max_pim_mom_values_CD[m];
    //     if (_pim_mom_2nd > mom_min && _pim_mom_2nd < mom_max) {
    //       //   // For experimental data
    //       // _pim_mom_prime = _pim_mom - pim_mom_corr_CD[m] * alpha_pim_mom_corr;
    //       _pim_mom_prime_2nd = _pim_mom_2nd - pim_mom_corr_CD_2nd[m] * alpha_pim_mom_corr_2nd[0];
    //     }
    //   }
    // }

    // _px_prime_pim_mom = _pim->Px() * ((_pim_mom_prime_2nd) / (_pim_mom));
    // _py_prime_pim_mom = _pim->Py() * ((_pim_mom_prime_2nd) / (_pim_mom));
    // _pz_prime_pim_mom = _pim->Pz() * ((_pim_mom_prime_2nd) / (_pim_mom));
    // _mom_corr_pim->SetXYZM(_px_prime_pim_mom, _py_prime_pim_mom, _pz_prime_pim_mom, MASS_PIM);
  }
  // bool Reaction::ctof_pim() {
  //   bool _pim_ctof = true;
  //   _pim_ctof &= (4000 <= _pim_status && _pim_status < 6000);
  //   return _pim_ctof;
  // }
  // bool Reaction::ftof_pim() {
  //   bool _pim_ftof = true;
  //   _pim_ftof &= (2000 <= _pim_status && _pim_status < 4000);
  //   return _pim_ftof;
  // }

  // float Reaction::rec_pim_px() {
  //   return _beam->Px() - _elec->Px() + _target->Px() - _pip->Px() - _prot->Px() - _pim->Px();
  // }
  // float Reaction::rec_pim_py() {
  //   return _beam->Py() - _elec->Py() + _target->Py() - _pip->Py() - _prot->Py() - _pim->Py();
  // }
  // float Reaction::rec_pim_pz() {
  //   return _beam->Pz() - _elec->Pz() + _target->Pz() - _pip->Pz() - _prot->Pz() - _pim->Pz();
  // }
  // float Reaction::rec_pim_E() { return _beam->E() - _elec->E() + _target->E() - _pip->E() - _prot->E() -
  // _pim->E(); } float Reaction::rec_pim_P() {
  //   return sqrt(abs(pow((_beam->Px() - _elec->Px() + _target->Px() - _pip->Px() - _prot->Px() - _pim->Px()),
  //   2)
  //   +
  //                   pow((_beam->Py() - _elec->Py() + _target->Py() - _pip->Py() - _prot->Py() - _pim->Py()),
  //                   2)
  //                   + pow((_beam->Pz() - _elec->Pz() + _target->Pz() - _pip->Pz() - _prot->Pz() -
  //                   _pim->Pz()), 2)));
  // }

  // float Reaction::rec_pim_mm2() {
  //   return abs(pow(_beam->Px() - _elec->Px() + _target->Px() - _pip->Px() - _prot->Px() - _pim->Px(), 2) +
  //              pow(_beam->Py() - _elec->Py() + _target->Py() - _pip->Py() - _prot->Py() - _pim->Py(), 2) +
  //              pow(_beam->Pz() - _elec->Pz() + _target->Pz() - _pip->Pz() - _prot->Pz() - _pim->Pz(), 2) -
  //              pow(_beam->E() - _elec->E() + _target->E() - _pip->E() - _prot->E() - _pim->E(), 2));
  // }

  // float Reaction::beam_px() { return _beam->Px(); }
  // float Reaction::beam_py() { return _beam->Py(); }
  // float Reaction::beam_pz() { return _beam->Pz(); }
  // float Reaction::beam_E() { return _beam->E(); }

  // float Reaction::elec_px() { return _elec->Px(); }
  // float Reaction::elec_py() { return _elec->Py(); }
  // float Reaction::elec_pz() { return _elec->Pz(); }
  // float Reaction::elec_E() { return _elec->E(); }

  // float Reaction::target_px() { return _target->Px(); }
  // float Reaction::target_py() { return _target->Py(); }
  // float Reaction::target_pz() { return _target->Pz(); }
  // float Reaction::target_E() { return _target->E(); }

  // float Reaction::pim_px() { return _pim->Px(); }
  // float Reaction::pim_py() { return _pim->Py(); }
  // float Reaction::pim_pz() { return _pim->Pz(); }
  // float Reaction::pim_E() { return _pim->E(); }
  // float Reaction::pim_P() { return _pim->P(); }

  // float Reaction::pip_px() { return _pip->Px(); }
  // float Reaction::pip_py() { return _pip->Py(); }
  // float Reaction::pip_pz() { return _pip->Pz(); }
  // float Reaction::pip_E() { return _pip->E(); }

  // float Reaction::prot_px() { return _prot->Px(); }
  // float Reaction::prot_py() { return _prot->Py(); }
  // float Reaction::prot_pz() { return _prot->Pz(); }
  // float Reaction::prot_E() { return _prot->E(); }

  void Reaction::SetNeutron(int i) {
    _numNeutral++;
    _hasNeutron = true;
    _neutron->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), MASS_N);
  }

  void Reaction::SetOther(int i) {
    if (_data->pid(i) == NEUTRON) {
      SetNeutron(i);
    } else {
      _numOther++;
      _hasOther = true;
      _other->SetXYZM(_data->px(i), _data->py(i), _data->pz(i), mass[_data->pid(i)]);
    }
  }

  void Reaction::CalcMissMass() {
    auto mm = std::make_unique<TLorentzVector>();
    auto mm_mpip = std::make_unique<TLorentzVector>();
    auto mm_mprot = std::make_unique<TLorentzVector>();
    auto mm_excl = std::make_unique<TLorentzVector>();

    *mm += (*_gamma + *_target);

    // if (TwoPion_missingPim()) {
    //   *mm -= *_prot;
    //   *mm -= *_pip;
    //   // *mm -= *_pim;
    //   _MM = mm->M();
    //   _MM2 = mm->M2();

    // //   // _rec_pim_mom = mm->P();
    // //   // _rec_pim_theta = mm->Theta() * 180 / PI;

    // //   // if (mm->Phi() >= 0)
    // //   //   _rec_pim_phi = (mm->Phi() * 180 / PI);
    // //   // else if (mm->Phi() < 0)
    // //   //   _rec_pim_phi = ((mm->Phi() + 2 * PI) * 180 / PI);

    // // //   // // // _x_mu_E = mm->E();
    // // //   // // // _x_mu_P = mm->P();
    // // //   // // // _x_mu_Px = mm->Px();
    // // //   // // // _x_mu_Py = mm->Py();
    // // //   // // // _x_mu_Pz = mm->Pz();
    // // //   // // // _x_mu_theta = mm->Theta() * RAD2DEG;
    // // //   // // // _x_mu_m2 = mm->E() * mm->E() - mm->P() * mm->P();
    // // //   // // // _x_mu_m = mm->E() - mm->P();
    // // //   // // //   //
    // }
    if (TwoPion_exclusive()) {
      // *mm -= *_mom_corr_prot;
      // *mm -= *_mom_corr_pip;
      // // *mm -= *_pim;
      // _MM = mm->M();
      // _MM2 = mm->M2();

      // *mm_excl += (*_gamma + *_target);
      // *mm_excl -= *_mom_corr_prot;
      // *mm_excl -= *_mom_corr_pip;
      // *mm_excl -= *_mom_corr_pim;

      *mm -= *_prot;
      *mm -= *_pip;
      // *mm -= *_pim;
      _MM = mm->M();
      _MM2 = mm->M2();

      *mm_excl += (*_gamma + *_target);
      *mm_excl -= *_prot;
      *mm_excl -= *_pip;
      *mm_excl -= *_pim;

      _MM2_exclusive = mm_excl->M2();
      _excl_Energy = mm_excl->E();

      // _rec_pim_mom = mm->P();
      // _rec_pim_theta = mm->Theta() * 180 / PI;

      // if (mm->Phi() >= 0)
      //   _rec_pim_phi = (mm->Phi() * 180 / PI);
      // else if (mm->Phi() < 0)
      //   _rec_pim_phi = ((mm->Phi() + 2 * PI) * 180 / PI);

      // // //   // //////// for x_mu - elec/beam theta phi
      // // //   // if (mm_excl->Phi() >= 0)
      // // //   //   _x_mu_phi = (mm_excl->Phi() * 180 / PI);
      // // //   // else if (mm_excl->Phi() < 0)
      // // //   //   _x_mu_phi = ((mm_excl->Phi() + 2 * PI) * 180 / PI);

      // // //   // if (_elec->Phi() >= 0)
      // // //   //   _elec_phi = (_elec->Phi() * 180 / PI);
      // // //   // else if (_elec->Phi() < 0)
      // // //   //   _elec_phi = ((_elec->Phi() + 2 * PI) * 180 / PI);

      // // //   // if (_beam->Phi() >= 0)
      // // //   //   _beam_phi = (_beam->Phi() * 180 / PI);
      // // //   // else if (_beam->Phi() < 0)
      // // //   //   _beam_phi = ((_beam->Phi() + 2 * PI) * 180 / PI);

      // // //   // _diff_elec_x_mu_theta = (_elec->Theta() * 180 / PI);  // - (mm_excl->Theta() * 180 / PI);
      // // //   // _diff_elec_x_mu_phi = (_elec_phi - _x_mu_phi);

      // // //   // _diff_beam_x_mu_theta = (_beam->Theta() * 180 / PI);  //-(mm_excl->Theta() * 180 / PI);
      // // //   // _diff_beam_x_mu_phi = (_beam_phi - _x_mu_phi);

      // // //   // // std::cout << " beam_theta " << _diff_beam_x_mu_theta << std::endl;
      // // //   // // std::cout << " rec_pim_energy " << mm->E() << std::endl;

      //   // //   // for mPip peak with exclusive events
      //   *mm_mpip += (*_gamma + *_target);
      //   *mm_mpip -= *_mom_corr_prot;
      //   *mm_mpip -= *_mom_corr_pim;
      //   _MM2_mPip = mm_mpip->M2();

      //   // //   // for mProt peak with exclusive events
      //   *mm_mprot += (*_gamma + *_target);
      //   *mm_mprot -= *_mom_corr_pip;
      //   *mm_mprot -= *_mom_corr_pim;
      //   _MM2_mProt = mm_mprot->M2();
      // }
      // // if (TwoPion_missingPip()) {
      *mm_mpip += (*_gamma + *_target);
      *mm_mpip -= *_prot;
      *mm_mpip -= *_pim;
      _MM2_mPip = mm_mpip->M2();
      // // }
      // // if (TwoPion_missingProt()) {
      *mm_mprot += (*_gamma + *_target);
      *mm_mprot -= *_pip;
      *mm_mprot -= *_pim;
      _MM2_mProt = mm_mprot->M2();
    }
  }
  // float Reaction::Diff_elec_x_mu_theta() {
  //   if (_diff_elec_x_mu_theta != _diff_elec_x_mu_theta) CalcMissMass();
  //   return _diff_elec_x_mu_theta;
  // }

  // float Reaction::Diff_elec_x_mu_phi() {
  //   if (_diff_elec_x_mu_phi != _diff_elec_x_mu_phi) CalcMissMass();
  //   return _diff_elec_x_mu_phi;
  // }

  // float Reaction::Diff_beam_x_mu_theta() {
  //   if (_diff_beam_x_mu_theta != _diff_beam_x_mu_theta) CalcMissMass();
  //   return _diff_beam_x_mu_theta;
  // }

  // float Reaction::Diff_beam_x_mu_phi() {
  //   if (_diff_beam_x_mu_phi != _diff_beam_x_mu_phi) CalcMissMass();
  //   return _diff_beam_x_mu_phi;
  // }

  float Reaction::MM() {
    if (_MM != _MM) CalcMissMass();
    return _MM;
  }
  float Reaction::MM2() {
    if (_MM2 != _MM2) CalcMissMass();
    return _MM2;
  }
  float Reaction::MM2_exclusive() {
    if (_MM2_exclusive != _MM2_exclusive) CalcMissMass();
    return _MM2_exclusive;
  }
  float Reaction::MM2_mPip() {
    if (_MM2_mPip != _MM2_mPip) CalcMissMass();
    return _MM2_mPip;
  }
  float Reaction::MM2_mProt() {
    if (_MM2_mProt != _MM2_mProt) CalcMissMass();
    return _MM2_mProt;
  }

  float Reaction::MM2_mPim_corr() {
    // if (_rec_pim_mom != _rec_pim_mom) CalcMissMass();

    // if (TwoPion_missingPim()) {
    if (TwoPion_exclusive()) {
      auto missingpim_ = std::make_unique<TLorentzVector>();
      // *missingpim_ += *_gamma + *_target - *_prot - *_pip;
      *missingpim_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pip;

      return missingpim_->M2();
      // return _rec_pim_mom;

    } else
      return NAN;
  }

  float Reaction::MM2_mPip_corr() {
    // if (TwoPion_missingPip()) {
    if (TwoPion_exclusive()) {
      auto missingpip_ = std::make_unique<TLorentzVector>();
      // *missingpip_ += *_gamma + *_target - *_prot - *_pim;
      *missingpip_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pim;

      return missingpip_->M2();
    } else
      return NAN;
  }

  float Reaction::MM2_mProt_corr() {
    // if (TwoPion_missingProt()) {
    if (TwoPion_exclusive()) {
      auto missingprot_ = std::make_unique<TLorentzVector>();
      // *missingprot_ += *_gamma + *_target - *_pip - *_pim;
      *missingprot_ += *_gamma + *_target - *_mom_corr_pip - *_mom_corr_pim;

      return missingprot_->M2();
    } else
      return NAN;
  }

  float Reaction::Energy_excl() {
    if (_excl_Energy != _excl_Energy) CalcMissMass();
    //  std::cout << "_x_mu_p  " << _x_mu->E() << '\n';
    //  if (_x_mu_E > 0)
    return _excl_Energy;
    // else
    // return NAN;
  }
  float Reaction::pim_momentum() {
    // if (_rec_pim_mom != _rec_pim_mom) CalcMissMass();

    // if (TwoPion_missingPim()) {
    if (TwoPion_exclusive()) {
      auto missingpim_ = std::make_unique<TLorentzVector>();
      // *missingpim_ += *_gamma + *_target - *_prot - *_pip;
      *missingpim_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pip;

      return missingpim_->P();
      // return _rec_pim_mom;

    } else
      return NAN;
  }
  float Reaction::pim_theta_lab() {
    // if (_rec_pim_theta != _rec_pim_theta) CalcMissMass();

    // if (TwoPion_missingPim()) {
    if (TwoPion_exclusive()) {
      auto missingpim_ = std::make_unique<TLorentzVector>();
      *missingpim_ += *_gamma + *_target - *_prot - *_pip;
      // *missingpim_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pip;

      return missingpim_->Theta() * 180.0 / PI;
      // return _rec_pim_theta;
    } else
      return NAN;
  }
  float Reaction::pim_Phi_lab() {
    // if (_rec_pim_phi != _rec_pim_phi) CalcMissMass();

    // if (TwoPion_missingPim()) {
    if (TwoPion_exclusive()) {
      auto missingpim_ = std::make_unique<TLorentzVector>();
      *missingpim_ += *_gamma + *_target - *_prot - *_pip;
      // *missingpim_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pip;

      if (missingpim_->Phi() > 0)
        return missingpim_->Phi() * 180 / PI;
      else if (missingpim_->Phi() < 0)
        return (missingpim_->Phi() + 2 * PI) * 180 / PI;
      else
        return NAN;
      // return _rec_pim_phi;
    } else
      return NAN;
  }
  float Reaction::pim_momentum_measured() {
    if (TwoPion_exclusive())
      return _pim->P();
    else
      return NAN;
  }

  float Reaction::pim_theta_lab_measured() {
    if (TwoPion_exclusive())
      return _pim->Theta() * 180.0 / PI;
    else
      return NAN;
  }

  float Reaction::pim_Phi_lab_measured() {
    if (TwoPion_exclusive()) {
      if (_pim->Phi() > 0) {
        // std::cout << "phi root >0 is " << _pim->Phi() * 180 / PI << std::endl;
        return _pim->Phi() * 180 / PI;
      } else if (_pim->Phi() < 0) {
        // std::cout << "phi root < 0 is " << (_pim->Phi() + 2 * PI) * 180 / PI << std::endl;
        return (_pim->Phi() + 2 * PI) * 180 / PI;
      } else
        return NAN;
    } else
      return NAN;
  }

  float Reaction::pim_momentum_corrected() {
    if (TwoPion_exclusive())
      return _mom_corr_pim->P();
    else
      return NAN;
  }
  float Reaction::w_hadron() {
    if (TwoPion_exclusive())
      return ((*_prot) + (*_pip) + (*_pim)).Mag();
    else
      return NAN;
  }
  // float Reaction::w_difference() {
  //   if (TwoPion_exclusive())
  //     return (physics::W_calc(*_beam, *_mom_corr_elec) - ((*_prot) + (*_pip) + (*_pim)).Mag());
  //   else
  //     return NAN;
  // }

  float Reaction::w_hadron_corr() {
    if (TwoPion_exclusive())
      return ((*_mom_corr_prot) + (*_mom_corr_pip) + (*_mom_corr_pim)).Mag();
    else
      return NAN;
  }
  // float Reaction::w_difference_corr() {
  //   if (TwoPion_exclusive())
  //     return (physics::W_calc(*_beam, *_mom_corr_elec) -
  //             ((*_mom_corr_prot) + (*_mom_corr_pip) + (*_mom_corr_pim)).Mag());
  //   else
  //     return NAN;
  // }

  // float Reaction::pim_theta_corrected() {
  //   if (TwoPion_exclusive())
  //     return _mom_corr_pim->Theta() * 180.0 / PI;
  //   else
  //     return NAN;
  // }

  // float Reaction::pim_Phi_corrected() {
  //   if (TwoPion_exclusive()) {
  //     if (_mom_corr_pim->Phi() > 0)
  //       return _mom_corr_pim->Phi() * 180 / PI;
  //     else if (_mom_corr_pim->Phi() < 0)
  //       return (_mom_corr_pim->Phi() + 2 * PI) * 180 / PI;
  //     else
  //       return NAN;
  //   } else
  //     return NAN;
  // }
  ////////////////mPip
  float Reaction::pip_momentum() {
    // if (TwoPion_missingPip()) {
    if (TwoPion_exclusive()) {
      auto missingpip_ = std::make_unique<TLorentzVector>();
      // *missingpip_ += *_gamma + *_target - *_prot - *_pim;
      *missingpip_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pim;

      return missingpip_->P();
    } else
      return NAN;
  }
  float Reaction::pip_theta_lab() {
    // if (TwoPion_missingPip()) {
    if (TwoPion_exclusive()) {
      auto missingpip_ = std::make_unique<TLorentzVector>();
      *missingpip_ += *_gamma + *_target - *_prot - *_pim;
      // *missingpip_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pim;
      return missingpip_->Theta() * 180.0 / PI;
    } else
      return NAN;
  }
  float Reaction::pip_Phi_lab() {
    // if (TwoPion_missingPip()) {
    if (TwoPion_exclusive()) {
      auto missingpip_ = std::make_unique<TLorentzVector>();
      *missingpip_ += *_gamma + *_target - *_prot - *_pim;
      // *missingpip_ += *_gamma + *_target - *_mom_corr_prot - *_mom_corr_pim;

      if (missingpip_->Phi() > 0)
        return missingpip_->Phi() * 180 / PI;
      else if (missingpip_->Phi() < 0)
        return (missingpip_->Phi() + 2 * PI) * 180 / PI;
      else
        return NAN;
    } else
      return NAN;
  }
  float Reaction::pip_momentum_measured() {
    if (TwoPion_exclusive())
      return _pip->P();
    else
      return NAN;
  }

  float Reaction::pip_theta_lab_measured() {
    if (TwoPion_exclusive())
      return _pip->Theta() * 180.0 / PI;
    else
      return NAN;
  }

  float Reaction::pip_Phi_lab_measured() {
    if (TwoPion_exclusive()) {
      if (_pip->Phi() > 0)
        return _pip->Phi() * 180 / PI;
      else if (_pip->Phi() < 0)
        return (_pip->Phi() + 2 * PI) * 180 / PI;
      else
        return NAN;
    } else
      return NAN;
  }

  float Reaction::pip_momentum_corrected() {
    if (TwoPion_exclusive())
      return _mom_corr_pip->P();
    else
      return NAN;
  }
  // float Reaction::pip_theta_corrected() {
  //   if (TwoPion_exclusive())
  //     return _mom_corr_pip->Theta() * 180.0 / PI;
  //   else
  //     return NAN;
  // }

  // float Reaction::pip_Phi_corrected() {
  //   if (TwoPion_exclusive()) {
  //     if (_mom_corr_pip->Phi() > 0)
  //       return _mom_corr_pip->Phi() * 180 / PI;
  //     else if (_mom_corr_pip->Phi() < 0)
  //       return (_mom_corr_pip->Phi() + 2 * PI) * 180 / PI;
  //     else
  //       return NAN;
  //   } else
  //     return NAN;
  // }

  ////////////////mProt
  float Reaction::prot_momentum() {
    // if (TwoPion_missingProt()) {
    if (TwoPion_exclusive()) {
      auto missingprot_ = std::make_unique<TLorentzVector>();
      // *missingprot_ += *_gamma + *_target - *_pip - *_pim;
      *missingprot_ += *_gamma + *_target - *_mom_corr_pip - *_mom_corr_pim;

      return missingprot_->P();
    } else
      return NAN;
  }
  float Reaction::prot_theta_lab() {
    // if (TwoPion_missingProt()) {
    if (TwoPion_exclusive()) {
      auto missingprot_ = std::make_unique<TLorentzVector>();
      *missingprot_ += *_gamma + *_target - *_pip - *_pim;
      // *missingprot_ += *_gamma + *_target - *_mom_corr_pip - *_mom_corr_pim;

      return missingprot_->Theta() * 180.0 / PI;
    } else
      return NAN;
  }
  float Reaction::prot_Phi_lab() {
    // if (TwoPion_missingProt()) {
    if (TwoPion_exclusive()) {
      auto missingprot_ = std::make_unique<TLorentzVector>();
      *missingprot_ += *_gamma + *_target - *_pip - *_pim;
      // *missingprot_ += *_gamma + *_target - *_mom_corr_pip - *_mom_corr_pim;

      if (missingprot_->Phi() > 0)
        return missingprot_->Phi() * 180 / PI;
      else if (missingprot_->Phi() < 0)
        return (missingprot_->Phi() + 2 * PI) * 180 / PI;
      else
        return NAN;
    } else
      return NAN;
  }
  float Reaction::prot_momentum_measured() {
    if (TwoPion_exclusive())
      return _prot->P();
    else
      return NAN;
  }

  float Reaction::prot_theta_lab_measured() {
    if (TwoPion_exclusive())
      return _prot->Theta() * 180.0 / PI;
    else
      return NAN;
  }

  float Reaction::prot_Phi_lab_measured() {
    if (TwoPion_exclusive()) {
      if (_prot->Phi() > 0)
        return _prot->Phi() * 180 / PI;
      else if (_prot->Phi() < 0)
        return (_prot->Phi() + 2 * PI) * 180 / PI;
      else
        return NAN;
    } else
      return NAN;
  }

  float Reaction::prot_momentum_corrected() {
    if (TwoPion_exclusive())
      return _mom_corr_prot->P();
    else
      return NAN;
  }
  // float Reaction::prot_theta_corrected() {
  //   if (TwoPion_exclusive())
  //     return _mom_corr_prot->Theta() * 180.0 / PI;
  //   else
  //     return NAN;
  // }

  // float Reaction::prot_Phi_corrected() {
  //   if (TwoPion_exclusive()) {
  //     if (_mom_corr_prot->Phi() > 0)
  //       return _mom_corr_prot->Phi() * 180 / PI;
  //     else if (_mom_corr_prot->Phi() < 0)
  //       return (_mom_corr_prot->Phi() + 2 * PI) * 180 / PI;
  //     else
  //       return NAN;
  //   } else
  //     return NAN;
  // }

  /////////////////////
  std::string Reaction::CsvHeader() { return "e_rec_p,e_rec_theta,e_rec_phi,e_sec\n"; }
  std::string Reaction::ReacToCsv() {
    // e_rec_p,e_rec_theta,e_rec_phi,e_sec
    std::string out = "";
    out += std::to_string(_elec->P()) + ",";
    out += std::to_string(_elec->Theta()) + ",";
    out += std::to_string(_elec->Phi()) + ",";
    out += std::to_string(_sector) + "\n";

    return out;
  }

  void Reaction::boost() {
    _is_boosted = true;
    _boosted_prot = std::make_unique<TLorentzVector>(*_prot);
    _boosted_pip = std::make_unique<TLorentzVector>(*_pip);
    _boosted_pim = std::make_unique<TLorentzVector>(*_gamma + *_target - *_prot - *_pip);  //(*_pim);
    _boosted_gamma = std::make_unique<TLorentzVector>(*_gamma);
    _boosted_pim_measured = std::make_unique<TLorentzVector>(*_pim);

    _rotated_prot = std::make_unique<TLorentzVector>(*_prot);
    _rotated_pip = std::make_unique<TLorentzVector>(*_pip);
    _rotated_pim = std::make_unique<TLorentzVector>(*_gamma + *_target - *_prot - *_pip);  //(*_pim);
    _rotated_pim_measured = std::make_unique<TLorentzVector>(*_pim);

    TRotation rot;
    _boosted_gamma->Transform(rot);
    float_t beta_1 = ((sqrt(_boosted_gamma->E() * _boosted_gamma->E() + _Q2)) / (_boosted_gamma->E() + MASS_P));
    TVector3 uz = _boosted_gamma->Vect().Unit();                  // uit vector along virtual photon
    TVector3 ux = ((_beam->Vect()).Cross(_elec->Vect())).Unit();  // unit vector along e cross e'
    ux.Rotate(3. * PI / 2, uz);                                   // rotating ux by 3pi/2 with uz as axis of roration
    rot.SetZAxis(uz, ux).Invert();                                // setting TRotation rot

    _boosted_prot->Transform(rot);
    _rotated_prot->Transform(rot);
    _boosted_prot->Boost(0, 0, -beta_1);

    _boosted_pip->Transform(rot);
    _rotated_pip->Transform(rot);
    _boosted_pip->Boost(0, 0, -beta_1);

    _boosted_pim->Transform(rot);
    _rotated_pim->Transform(rot);
    _boosted_pim->Boost(0, 0, -beta_1);

    _boosted_gamma->Boost(0, 0, -beta_1);

    _boosted_pim_measured->Transform(rot);
    _rotated_pim_measured->Transform(rot);
    _boosted_pim_measured->Boost(0, 0, -beta_1);
    // -beta ko value (0.5 to -0.5 huda
    // samma value aauchha nattra aaudyna)

    _prot_Vect3 = _boosted_prot->Vect();
    _pip_Vect3 = _boosted_pip->Vect();
    _pim_Vect3 = _boosted_pim_measured->Vect();
  }

  float_t Reaction::scalar_triple_product() {
    if (!_is_boosted) boost();
    if (TwoPion_exclusive()) {
      return (_prot_Vect3.Dot(_pip_Vect3.Cross(_pim_Vect3)));

    } else
      return NAN;
  }

  // // float Reaction::pim_momentum_cm() {
  // //         if (!_is_boosted)
  // //                 boost();
  // //         if (TwoPion_missingPim())
  // //                 return _boosted_pim->P();
  // //         else
  // //                 return NAN;
  // // }

  // float Reaction::pim_theta_cm() {
  //   if (!_is_boosted) boost();
  //   if (TwoPion_missingPim())
  //     return _rotated_pim->Theta() * 180.0 / PI;
  //   else
  //     return NAN;
  // }

  // float Reaction::pim_Phi_cm() {
  //   if (!_is_boosted) boost();
  //   if (TwoPion_missingPim()) {
  //     if (_rotated_pim->Phi() > 0)
  //       return _rotated_pim->Phi() * 180 / PI;
  //     else if (_rotated_pim->Phi() < 0)
  //       return (_rotated_pim->Phi() + 2 * PI) * 180 / PI;
  //     else
  //       return NAN;
  //   } else
  //     return NAN;
  // }

  // // float Reaction::pim_momentum_cm_measured() {
  // //         if (!_is_boosted)
  // //                 boost();
  // //         if (TwoPion_exclusive())
  // //                 return _boosted_pim_measured->P();
  // //         else
  // //                 return NAN;
  // // }

  // float Reaction::pim_theta_cm_measured() {
  //   if (!_is_boosted) boost();
  //   if (TwoPion_exclusive())
  //     return _rotated_pim_measured->Theta() * 180.0 / PI;
  //   else
  //     return NAN;
  // }

  // float Reaction::pim_Phi_cm_measured() {
  //   if (!_is_boosted) boost();
  //   if (TwoPion_exclusive()) {
  //     if (_rotated_pim_measured->Phi() > 0)
  //       return _rotated_pim_measured->Phi() * 180 / PI;
  //     else if (_rotated_pim_measured->Phi() < 0)
  //       return (_rotated_pim_measured->Phi() + 2 * PI) * 180 / PI;
  //     else
  //       return NAN;
  //   } else
  //     return NAN;
  // }

  MCReaction::MCReaction(const std::shared_ptr<Branches12>& data, float beam_enrgy) {
    _data = data;
    if (!_data->mc()) _data->mc_branches();
    _beam = std::make_unique<TLorentzVector>();
    _beam_energy = beam_enrgy;
    _weight_mc = _data->mc_weight();
    _beam->SetPxPyPzE(0.0, 0.0, sqrt(_beam_energy * _beam_energy - MASS_E * MASS_E), _beam_energy);

    //_gamma = std::make_unique<TLorentzVector>();  // do i need this?
    _gamma_mc = std::make_unique<TLorentzVector>();
    _target = std::make_unique<TLorentzVector>(0.0, 0.0, 0.0, MASS_P);
    //_elec = std::make_unique<TLorentzVector>();  // do i need this?
    _elec_mc = std::make_unique<TLorentzVector>();
    // this->SetElec();  // do i need this?
    this->SetMCElec();
    _prot_mc = std::make_unique<TLorentzVector>();
    _pip_mc = std::make_unique<TLorentzVector>();
    _pim_mc = std::make_unique<TLorentzVector>();
    //_other = std::make_unique<TLorentzVector>();  // do i need this?
    _other_mc = std::make_unique<TLorentzVector>();
    //_neutron = std::make_unique<TLorentzVector>();
  }
  // Reaction::~Reaction() {} // why this is not here
  void MCReaction::SetMCElec() {
    //  _hasE = true;  //??
    _elec_mc->SetXYZM(_data->mc_px(0), _data->mc_py(0), _data->mc_pz(0), MASS_E);

    *_gamma_mc += *_beam - *_elec_mc;

    // Can calculate W and Q2 here
    _W_mc = physics::W_calc(*_beam, *_elec_mc);
    _Q2_mc = physics::Q2_calc(*_beam, *_elec_mc);
  }

  void MCReaction::SetMCProton(int i) { _prot_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_P); }

  void MCReaction::SetMCPip(int i) { _pip_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIP); }

  void MCReaction::SetMCPim(int i) { _pim_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i), MASS_PIM); }
  // void MCReaction::SetMCOther(int i) {
  //   _other_mc->SetXYZM(_data->mc_px(i), _data->mc_py(i), _data->mc_pz(i),
  //   mass[_data->pid(i)]);
  // }

  float MCReaction::pim_mom_mc_gen() {
    // if (Reaction::TwoPion_exclusive())
    return _pim_mc->P();
    // else
    //   return NAN;
  }
  float MCReaction::pip_mom_mc_gen() {
    // if (Reaction::TwoPion_exclusive())
    return _pip_mc->P();
    // else
    //   return NAN;
  }
  float MCReaction::prot_mom_mc_gen() {
    // if (Reaction::TwoPion_exclusive())
    return _prot_mc->P();
    // else
    //   return NAN;
  }

  float MCReaction::pim_theta_mc_gen() {
    if (Reaction::TwoPion_exclusive())
      return _pim_mc->Theta() * 180 / PI;
    else
      return NAN;
  }
  float MCReaction::pip_theta_mc_gen() {
    if (Reaction::TwoPion_exclusive())
      return _pip_mc->Theta() * 180 / PI;
    else
      return NAN;
  }
  float MCReaction::prot_theta_mc_gen() {
    if (TwoPion_exclusive())
      return _prot_mc->Theta() * 180 / PI;
    else
      return NAN;
  }

  float MCReaction::pim_phi_mc_gen() {
    if (_pim_mc->Phi() >= 0)
      return (_pim_mc->Phi() * 180 / PI);
    else if (_pim_mc->Phi() < 0)
      return ((_pim_mc->Phi() + 2 * PI) * 180 / PI);
    else
      return NAN;
  }
  float MCReaction::pip_phi_mc_gen() {
    if (_pip_mc->Phi() >= 0)
      return (_pip_mc->Phi() * 180 / PI);
    else if (_pip_mc->Phi() < 0)
      return ((_pip_mc->Phi() + 2 * PI) * 180 / PI);
    else
      return NAN;
  }
  float MCReaction::prot_phi_mc_gen() {
    if (_prot_mc->Phi() >= 0)
      return (_prot_mc->Phi() * 180 / PI);
    else if (_prot_mc->Phi() < 0)
      return ((_prot_mc->Phi() + 2 * PI) * 180 / PI);
    else
      return NAN;
  }

  // void MCReaction::CalcMissMass_mc() {
  //   auto mm_excl_mc = std::make_unique<TLorentzVector>();

  //   *mm_excl_mc += (*_gamma_mc + *_target);
  //   *mm_excl_mc -= *_prot_mc;
  //   *mm_excl_mc -= *_pip_mc;
  //   *mm_excl_mc -= *_pim_mc;
  //   _MM2_exclusive_mc = mm_excl_mc->M2();
  //   _excl_Energy_mc = mm_excl_mc->E();

  // _rec_x_mu_mom_mc = mm_excl_mc->P();
  // _rec_x_mu_theta_mc = mm_excl_mc->Theta() * 180 / PI;

  // if (mm_excl_mc->Phi() >= 0)
  //   _x_mu_phi_mc = (mm_excl_mc->Phi() * 180 / PI);
  // else if (mm_excl_mc->Phi() < 0)
  //   _x_mu_phi_mc = ((mm_excl_mc->Phi() + 2 * PI) * 180 / PI);

  // if (_elec_mc->Phi() >= 0)
  //   _elec_phi_mc = (_elec_mc->Phi() * 180 / PI);
  // else if (_elec_mc->Phi() < 0)
  //   _elec_phi_mc = ((_elec_mc->Phi() + 2 * PI) * 180 / PI);

  // if (_beam->Phi() >= 0)
  //   _beam_phi_mc = (_beam->Phi() * 180 / PI);
  // else if (_beam->Phi() < 0)
  //   _beam_phi_mc = ((_beam->Phi() + 2 * PI) * 180 / PI);

  // _diff_elec_x_mu_theta_mc = (_elec_mc->Theta() * 180 / PI) - (mm_excl_mc->Theta() * 180 / PI);
  // _diff_elec_x_mu_phi_mc = (_elec_phi_mc - _x_mu_phi_mc);

  // _diff_beam_x_mu_theta_mc = (mm_excl_mc->Theta() * 180 / PI);
  // _diff_beam_x_mu_phi_mc = (_beam_phi_mc - _x_mu_phi_mc);
  // }

  // float MCReaction::Diff_elec_x_mu_theta_mc() {
  //   if (_diff_elec_x_mu_theta_mc != _diff_elec_x_mu_theta_mc) CalcMissMass_mc();
  //   return _diff_elec_x_mu_theta_mc;
  // }

  // float MCReaction::Diff_elec_x_mu_phi_mc() {
  //   if (_diff_elec_x_mu_phi_mc != _diff_elec_x_mu_phi_mc) CalcMissMass_mc();
  //   return _diff_elec_x_mu_phi_mc;
  // }

  // float MCReaction::Diff_beam_x_mu_theta_mc() {
  //   if (_diff_beam_x_mu_theta_mc != _diff_beam_x_mu_theta_mc) CalcMissMass_mc();
  //   return _diff_beam_x_mu_theta_mc;
  // }

  // float MCReaction::Diff_beam_x_mu_phi_mc() {
  //   if (_diff_beam_x_mu_phi_mc != _diff_beam_x_mu_phi_mc) CalcMissMass_mc();
  //   return _diff_beam_x_mu_phi_mc;
  // }

  // float MCReaction::MM2_exclusive_mc() {
  //   if (_MM2_exclusive_mc != _MM2_exclusive_mc) CalcMissMass_mc();
  //   return _MM2_exclusive_mc;
  // }
  // float MCReaction::Energy_excl_mc() {
  //   if (_excl_Energy_mc != _excl_Energy_mc) CalcMissMass_mc();
  //   return _excl_Energy_mc;
  // }
  // float MCReaction::x_mu_momentum_mc() {
  //   if (_rec_x_mu_mom_mc != _rec_x_mu_mom_mc) CalcMissMass_mc();
  //   return _rec_x_mu_mom_mc;
  // }
  // float MCReaction::x_mu_theta_lab_mc() {
  //   if (_rec_x_mu_theta_mc != _rec_x_mu_theta_mc) CalcMissMass_mc();
  //   return _rec_x_mu_theta_mc;
  // }
  // float MCReaction::x_mu_Phi_lab_mc() {
  //   if (_x_mu_phi_mc != _x_mu_phi_mc) CalcMissMass_mc();
  //   return _x_mu_phi_mc;
  // }

  std::string MCReaction::CsvHeader() {
    return "e_rec_p,e_rec_theta,e_rec_phi,e_sec,e_thrown_p,e_thrown_theta,e_thrown_phi\n";
  }
  std::string MCReaction::ReacToCsv() {
    // e_rec_p,e_rec_theta,e_rec_phi,e_sec,e_thrown_p,e_thrown_theta,e_thrown_phi
    std::string out = "";
    out += std::to_string(_elec->P()) + ",";
    out += std::to_string(_elec->Theta()) + ",";
    out += std::to_string(_elec->Phi()) + ",";
    out += std::to_string(_sector) + ",";
    out += std::to_string(_elec_mc->P()) + ",";
    out += std::to_string(_elec_mc->Theta()) + ",";
    out += std::to_string(_elec_mc->Phi()) + "\n";

    return out;
  }
