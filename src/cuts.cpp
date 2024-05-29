
#include "cuts.hpp"
#include <iostream>
#include "TFile.h"
#include "histogram.hpp"
#include "reaction.hpp"

Cuts::Cuts(const std::shared_ptr<Branches12> &data) : _data(data) { _dt = std::make_shared<Delta_T>(data); }
Cuts::Cuts(const std::shared_ptr<Branches12> &data, const std::shared_ptr<Delta_T> &dt) : _data(data), _dt(dt) {}
// size_t run(std::shared_ptr<TChain> _chain, std::shared_ptr<Histogram> _hists,
//            int thread_id);

Cuts::~Cuts() {}

bool Cuts::ElectronCuts() {
  bool _elec = true;
  // Number of good particles is greater than 0
  // So that we can check at(0) without errors
  _elec &= (_data->gpart() > 0);
  if (!_elec) return false;

  _elec &= (_data->gpart() < 20);

  _elec &= (_data->charge(0) == NEGATIVE);
  _elec &= (_data->pid(0) == ELECTRON);
  // //  _elec &= !std::isnan(_data->cc_nphe_tot(0));
  //
  //_elec &= (_data->beta(0) > 0.05);
  _elec &= (_data->p(0) > 1.50);
  _elec &= (2000 <= abs(_data->status(0)) && abs(_data->status(0)) < 4000);
  _elec &= (_data->vz(0) > -(2.78 + 2 * 2.16) && _data->vz(0) < (-2.78 + 2 * 2.16));  // 3 sigma cut
  // Use the chi2pid instead of straight line cuts on SF
  //_elec &= (abs(_data->chi2pid(0)) < 3);
  _elec &=
      (_data->ec_tot_energy(0) / _data->p(0) < (0.30676 - 0.00111 * _data->p(0) - 0.00031 * _data->p(0) * _data->p(0)));
  _elec &=
      (_data->ec_tot_energy(0) / _data->p(0) > (0.15546 + 0.01714 * _data->p(0) - 0.00151 * _data->p(0) * _data->p(0)));

  //
  // FiducialCuts is the slowest of the cuts because of all the calcuations
  // If it already fails a different cut we will quit before
  // calulating for the FiducialCuts to save time
  if (!_elec) return _elec;
  _elec &= FiducialCuts();

  return _elec;
}

bool Cuts::FiducialCuts() {
  bool _fid_cut = true;
  // DC sector never changes so get it once and store it to use all the time
  short dc_sec = (_data->dc_sec(0) - 1);
  // Same with these values
  float sin_dc_sec = sinf(dc_sec * ROTATE);
  float cos_dc_sec = cosf(dc_sec * ROTATE);

  float x_PCAL_rot = _data->ec_pcal_y(0) * sin_dc_sec + _data->ec_pcal_x(0) * cos_dc_sec;
  float y_PCAL_rot = _data->ec_pcal_y(0) * cos_dc_sec - _data->ec_pcal_x(0) * sin_dc_sec;

  float left_PCAL = (HEIGHT_PCAL - SLOPE * y_PCAL_rot);
  float right_PCAL = (HEIGHT_PCAL + SLOPE * y_PCAL_rot);
  float radius2_PCAL = X_SQUARE_PCAL - (y_PCAL_rot * y_PCAL_rot);  // circle radius r^2 = x^2 + y^2

  // I do this to clean up what is happening and makse sure that the cuts are
  // not ambiguous
  _fid_cut &= (x_PCAL_rot > left_PCAL);
  _fid_cut &= (x_PCAL_rot > right_PCAL);
  _fid_cut &= (x_PCAL_rot * x_PCAL_rot > radius2_PCAL);
  _fid_cut &= (x_PCAL_rot < 372);

  // If it fails pcal cut return before calculating DC cut to save time
  if (!_fid_cut) return _fid_cut;

  float x1_rot = _data->dc_r1_y(0) * sin_dc_sec + _data->dc_r1_x(0) * cos_dc_sec;
  float y1_rot = _data->dc_r1_y(0) * cos_dc_sec - _data->dc_r1_x(0) * sin_dc_sec;
  float left_r1 = (DCR1_HEIGHT - SLOPE * y1_rot);
  float right_r1 = (DCR1_HEIGHT + SLOPE * y1_rot);
  float radius2_DCr1 = DCR1_SQUARE - (y1_rot * y1_rot);

  _fid_cut &= (x1_rot > left_r1);
  _fid_cut &= (x1_rot > right_r1);
  _fid_cut &= (x1_rot * x1_rot > radius2_DCr1);

  // If it fails cut return before calculating cut to save time
  if (!_fid_cut) return _fid_cut;

  float x2_rot = _data->dc_r2_y(0) * sin_dc_sec + _data->dc_r2_x(0) * cos_dc_sec;
  float y2_rot = _data->dc_r2_y(0) * cos_dc_sec - _data->dc_r2_x(0) * sin_dc_sec;
  float left_r2 = (DCR2_HEIGHT - SLOPE * y2_rot);
  float right_r2 = (DCR2_HEIGHT + SLOPE * y2_rot);
  float radius2_DCr2 = DCR2_SQUARE - (y2_rot * y2_rot);

  _fid_cut &= (x2_rot > left_r2);
  _fid_cut &= (x2_rot > right_r2);
  _fid_cut &= ((x2_rot * x2_rot) > radius2_DCr2);

  // If it fails cut return before calculating cut to save time
  if (!_fid_cut) return _fid_cut;

  float x3_rot = _data->dc_r3_y(0) * sin_dc_sec + _data->dc_r3_x(0) * cos_dc_sec;
  float y3_rot = _data->dc_r3_y(0) * cos_dc_sec - _data->dc_r3_x(0) * sin_dc_sec;
  float left_r3 = (DCR3_HEIGHT - SLOPE * y3_rot);
  float right_r3 = (DCR3_HEIGHT + SLOPE * y3_rot);
  float radius2_DCr3 = DCR3_SQUARE - pow(y3_rot, 2);

  _fid_cut &= (x3_rot > left_r3);
  _fid_cut &= (x3_rot > right_r3);
  _fid_cut &= ((x3_rot * x3_rot) > radius2_DCr3);

  return _fid_cut;
}
/////////////////////////// EXP data dt cuts ////////////// prot, pip, pim /////////////////

double dt_cut_fd[3][6] = {
    {-2.46098382e-03, 5.50052526e-02, -4.62941277e-01, 1.81269058e+00, -3.26981471e+00, 2.53628001e+00},
    {-1.37230195e-03, 3.05561520e-02, -2.56302213e-01, 1.00260344e+00, -1.81912637e+00, 1.54216958e+00},
    {-1.28173274e-03, 2.78960131e-02, -2.26379356e-01, 8.41779657e-01, -1.40649467e+00, 1.13639326e+00}};

double dt_cut_cd[3][3] = {{0.06046542, -0.37050291, 0.83190914},
                          {0.01848491, -0.08887598, 0.44818322},
                          {0.02421193, -0.10836307, 0.44539544}};

// ///////////////////////// SIM data dt cuts ////////////// prot, pip, pim /////////////////
// double dt_cut_fd[3][6] = {
//     {-1.75409746e-03, 3.93868741e-02, -3.33545338e-01, 1.32093252e+00, -2.44830422e+00, 2.15632989e+00},
//     {-4.34008708e-04, 9.56377615e-03, -8.31534626e-02, 3.60812599e-01, -7.89520521e-01, 1.08811065e+00},
//     {0.0, 0.0, 0.0, 0.00309318, -0.05048557, 0.50591481}};

// double dt_cut_cd[3][3] = {
//     {0.08992734, -0.50475099, 0.96629123}, {0.01654657, -0.0938437, 0.45075823}, {0.01747117, -0.07756261,
//     0.4276971}};

bool Cuts::IsPip(int i) {
  if (_data->gpart() <= i) return false;
  bool _pip = true;
  //   _pip &= (_data->charge(i) == POSITIVE);
  _pip &= (_data->pid(i) == PIP);
  // _pip &= (abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) < 0.4);
  _pip &=
      (abs(_dt->dt_Pi(i)) < (dt_cut_cd[1][0] * pow(_data->p(i), 2) + dt_cut_cd[1][1] * _data->p(i) + dt_cut_cd[1][2]) ||
       abs(_dt->dt_ctof_Pi(i)) < (dt_cut_fd[1][0] * pow(_data->p(i), 5) + dt_cut_fd[1][1] * pow(_data->p(i), 4) +
                                  dt_cut_fd[1][2] * pow(_data->p(i), 3) + dt_cut_fd[1][3] * pow(_data->p(i), 2) +
                                  dt_cut_fd[1][4] * pow(_data->p(i), 1) + dt_cut_fd[1][4]));
  // _pip &= !(abs(_dt->dt_P(i)) < 0.5 || abs(_dt->dt_ctof_P(i)) < 0.2);
  _pip &= (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 6000);

  // // min/max mom cuts
  if (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000) {
    _pip &= (_data->p(i) > 0.5);
    // _pip &= (_data->p(i) < 4.6);
  } else if (abs(_data->status(i)) >= 4000) {
    _pip &= (_data->p(i) > 0.2);
    // _pip &= (_data->p(i) < 1.7);
  }
  // _pip &= (_data->p(i) > 0.2);

  return _pip;
}
bool Cuts::IsProton(int i) {
  if (_data->gpart() <= i) return false;
  bool _proton = true;
  //   _proton &= (_data->charge(i) == POSITIVE);
  _proton &= (_data->pid(i) == PROTON);
  // _proton &= (abs(_dt->dt_P(i)) < 0.5 || abs(_dt->dt_ctof_P(i)) < 0.4);
  // // _proton &= !(abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) < 0.2);
  _proton &=
      (abs(_dt->dt_P(i)) < (dt_cut_cd[0][0] * pow(_data->p(i), 2) + dt_cut_cd[0][1] * _data->p(i) + dt_cut_cd[0][2]) ||
       abs(_dt->dt_ctof_P(i)) < (dt_cut_fd[0][0] * pow(_data->p(i), 5) + dt_cut_fd[0][1] * pow(_data->p(i), 4) +
                                 dt_cut_fd[0][2] * pow(_data->p(i), 3) + dt_cut_fd[0][3] * pow(_data->p(i), 2) +
                                 dt_cut_fd[0][4] * pow(_data->p(i), 1) + dt_cut_fd[0][4]));
  _proton &= (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 6000);
  // min/max mom cuts
  if (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000) {
    _proton &= (_data->p(i) > 0.4);
    // _proton &= (_data->p(i) < 4.5);
  } else if (abs(_data->status(i)) >= 4000) {
    _proton &= (_data->p(i) > 0.4);  /// this 0.4 look harse when we do missing Pim channel
                                     // _proton &= (_data->p(i) < 2.0);
  }
  // _proton &= (_data->p(i) > 0.2);
  //_proton &= (abs(_data->chi2pid(i)) < 0.5);
  return _proton;
}

bool Cuts::IsPim(int i) {
  if (_data->gpart() <= i) return false;
  bool _pim = true;
  //   _pim &= (_data->charge(i) == NEGATIVE);
  _pim &= (_data->pid(i) == PIM);
  // _pim &= (abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) < 0.5);
  _pim &=
      (abs(_dt->dt_Pi(i)) < (dt_cut_cd[2][0] * pow(_data->p(i), 2) + dt_cut_cd[2][1] * _data->p(i) + dt_cut_cd[2][2]) ||
       abs(_dt->dt_ctof_Pi(i)) < (dt_cut_fd[2][0] * pow(_data->p(i), 5) + dt_cut_fd[2][1] * pow(_data->p(i), 4) +
                                  dt_cut_fd[2][2] * pow(_data->p(i), 3) + dt_cut_fd[2][3] * pow(_data->p(i), 2) +
                                  dt_cut_fd[2][4] * pow(_data->p(i), 1) + dt_cut_fd[2][4]));

  _pim &= (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 6000);
  // min / max mom cuts
  if (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000) {
    _pim &= (_data->p(i) > 0.5);
    // _pim &= (_data->p(i) < 4.5);
  } else if (abs(_data->status(i)) >= 4000) {
    _pim &= (_data->p(i) > 0.2);
    // _pim &= (_data->p(i) < 1.9);
  }
  // _pim &= (_data->p(i) > 0.2);

  return _pim;
}

// /////////////////////// Pass2_Cuts ///////////////////////
bool Pass2_Cuts::ElectronCuts() {
  bool cut = true;
  cut &= (_data->gpart() > 0);
  if (!cut) return false;

  cut &= (_data->gpart() < 20);
  // //
  cut &= (_data->charge(0) == NEGATIVE);
  cut &= (_data->pid(0) == ELECTRON);
  // cut &= (_data->p(0) > 1.50);
  // cut &= (2000 <= abs(_data->status(0)) && abs(_data->status(0)) < 4000);
  // cut &= DC_z_vertex_cut();
  // // cut &= (abs(_data->chi2pid(0)) < 3);  ////////////// check it.......
  // cut &= CC_nphe_cut();
  // cut &= PCAL_Minimum_Energy_cut();
  // cut &= PCAL_fiducial_cut_HX_HY();
  // // cut &= EC_outer_vs_EC_inner_cut();
  // cut &= EC_sampling_fraction_cut();
  // cut &= EC_hit_position_fiducial_cut_homogeneous();
  // cut &= DC_fiducial_cut_XY();
  return cut;
}
bool Pass2_Cuts::HadronsCuts(int i) {
  bool cut = true;
  // if (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000) cut &= DC_fiducial_cut_theta_phi(i);
  // cut &= Hadron_Delta_vz_cut(i);
  // cut &= Hadron_Chi2pid_cut(i);
  return cut;
}

bool Pass2_Cuts::CC_nphe_cut() {
  float nphe_min = 2;
  return (_data->cc_nphe_tot(0) > nphe_min);
}

bool Pass2_Cuts::PCAL_Minimum_Energy_cut() {
  double edep_tight = 0.06, edep_medium = 0.07, edep_loose = 0.09;
  return (_data->ec_pcal_energy(0) > edep_medium);
}
/////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

bool Pass2_Cuts::EC_outer_vs_EC_inner_cut() {
  short isector = (_data->ec_pcal_sec(0) - 1);

  double param_a_exp[9][6] = {
      {0.201232, 0.197423, 0.197381, 0.18777, 0.197711, 0.198306},   // <2 GeV
      {0.207553, 0.20437, 0.209949, 0.198804, 0.208762, 0.209544},   // 2−3 GeV
      {0.213829, 0.215049, 0.215857, 0.209816, 0.218481, 0.221003},  // 3−4 GeV
      {0.217145, 0.218486, 0.21988, 0.215048, 0.221863, 0.22733},    // 4−5 GeV
      {0.220458, 0.218645, 0.220504, 0.218401, 0.222802, 0.230317},  // 5−6 GeV
      {0.22359, 0.219856, 0.226815, 0.221764, 0.222474, 0.23842},    // 6−7 GeV
      {0.226479, 0.219286, 0.22881, 0.225656, 0.22175, 0.242428},    // 7−8 GeV
      {0.226668, 0.219087, 0.226924, 0.228833, 0.21954, 0.24472},    // 8−9 GeV
      {0.225488, 0.221795, 0.21997, 0.228161, 0.216218, 0.246185}    // >9 GeV
  };

  double param_b_exp[9][6] = {
      {-1.00947, -0.975077, -0.988142, -0.892892, -0.97406, -0.981174},  // <2 GeV
      {-1.07716, -1.01246, -1.10285, -1.00481, -1.07982, -1.08578},      // 2−3 GeV
      {-1.10165, -1.08631, -1.1188, -1.0906, -1.15156, -1.17655},        // 3−4 GeV
      {-1.12186, -1.1006, -1.14603, -1.11891, -1.16233, -1.22691},       // 4−5 GeV
      {-1.16741, -1.10554, -1.15965, -1.13694, -1.1562, -1.26216},       // 5−6 GeV
      {-1.22825, -1.13945, -1.26374, -1.16584, -1.14528, -1.38061},      // 6−7 GeV
      {-1.30402, -1.16559, -1.33568, -1.2161, -1.13725, -1.47588},       // 7−8 GeV
      {-1.37011, -1.21059, -1.38956, -1.28444, -1.11712, -1.57429},      // 8−9 GeV
      {-1.45688, -1.30088, -1.4082, -1.34519, -1.09426, -1.67625}        // >9 GeV
  };

  int momRangeIndex = getMomRangeIndex(_data->p(0));
  double a = param_a_exp[momRangeIndex][isector];
  double b = param_b_exp[momRangeIndex][isector];
  double sf_ecin = _data->ec_ecin_energy(0) / _data->p(0);
  double sf_pcal = _data->ec_pcal_energy(0) / _data->p(0);
  // std::cout << "ec sec " << isector << " a  " << a << " b  " << b << '\n';

  return sf_pcal > a + b * sf_ecin;

  // double param_a_sim[9][6] = {
  //     {0.20233, 0.202647, 0.200986, 0.201774, 0.201084, 0.201648},  // <2 GeV
  //     {0.212437, 0.21288, 0.211826, 0.213158, 0.212129, 0.212346},  // 2−3 GeV
  //     {0.219554, 0.220162, 0.219478, 0.220034, 0.219012, 0.219218}, // 3−4 GeV
  //     {0.224078, 0.224914, 0.223814, 0.224357, 0.224371, 0.224192}, // 4−5 GeV
  //     {0.22785, 0.228319, 0.227713, 0.228105, 0.22709, 0.227346},   // 5−6 GeV
  //     {0.230326, 0.230837, 0.230451, 0.230969, 0.229402, 0.229665}, // 6−7 GeV
  //     {0.23258, 0.233087, 0.2319, 0.233145, 0.232148, 0.232072},    // 7−8 GeV
  //     {0.232341, 0.233287, 0.231747, 0.232848, 0.231751, 0.232098}, // 8−9 GeV
  //     {0.22484, 0.233669, 0.236727, 0.233346, 0.233825, 0.23524}    // >9 GeV
  // };

  // double param_b_sim[9][6] = {
  //     {-0.949695, -0.956382, -0.934262, -0.940217, -0.933617, -0.936363}, // <2 GeV
  //     {-1.04219, -1.04446, -1.035, -1.04832, -1.03597, -1.03963},         // 2−3 GeV
  //     {-1.08307, -1.09327, -1.08529, -1.08885, -1.07962, -1.08145},       // 3−4 GeV
  //     {-1.11223, -1.11772, -1.10889, -1.11142, -1.11528, -1.11228},       // 4−5 GeV
  //     {-1.13177, -1.13684, -1.13048, -1.13521, -1.12421, -1.12983},       // 5−6 GeV
  //     {-1.14596, -1.14755, -1.14677, -1.15218, -1.13488, -1.14008},       // 6−7 GeV
  //     {-1.16058, -1.16113, -1.1519, -1.1639, -1.15457, -1.1564},          // 7−8 GeV
  //     {-1.15733, -1.16306, -1.14965, -1.15888, -1.15131, -1.15527},       // 8−9 GeV
  //     {-1.08707, -1.16674, -1.21746, -1.18771, -1.17542, -1.19943}        // >9 GeV
  // };
}
/////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

bool Pass2_Cuts::EC_sampling_fraction_cut() {
  int isec = (_data->ec_pcal_sec(0) - 1);
  double upper_lim_total = 0;
  double lower_lim_total = 0;

  double mean_minus_3_5_sigma[6][3] = {{-0.0001186, 0.0001892, 0.1942}, {-0.000856, 0.01084, 0.1637},
                                       {-0.001184, 0.014046, 0.1593},   {-0.001268, 0.01918, 0.1287},
                                       {-0.0002744, 0.003532, 0.1844},  {-0.001039, 0.012505, 0.1593}};

  double mean_plus_3_5_sigma[6][3] = {{-0.0004027, 0.001746, 0.2903}, {-9.36e-05, -0.000999, 0.2979},
                                      {-0.0003238, 0.00101, 0.2957},  {-4.303e-05, -0.0004702, 0.2954},
                                      {-0.0001818, 0.003223, 0.2742}, {-0.0002906, 0.0015335, 0.2883}};

  // mean -3 *sigma : {{-0.0001389 ,0.0003004 ,0.201 ,}, {-0.0008016 ,0.01 ,0.1733 ,}, {-0.001123 ,0.013115 ,0.169 ,},
  // {-0.001181 ,0.01778 ,0.1405 ,}, {-0.0002677 ,0.00351 ,0.1908 ,}, {-0.000985 ,0.01172 ,0.1685 ,}, }
  //  maen + 3*sigma : {{-0.0003824 ,0.001636 ,0.2834 ,}, {-0.000148 ,-0.0001535 ,0.2883 ,}, {-0.0003853 ,0.001941
  //  ,0.286 ,}, {-0.0001305 ,0.0009336 ,0.2834 ,}, {-0.0001885 ,0.003246 ,0.2678 ,}, {-0.0003443 ,0.002317 ,0.279 ,}, }

  // mean -4 *sigma : {{-9.835e-05 ,7.79e-05 ,0.1874 ,}, {-0.0009108 ,0.01169 ,0.1542 ,}, {-0.001246 ,0.01498 ,0.1495
  // ,}, {-0.001355 ,0.02058 ,0.1167 ,}, {-0.000281 ,0.003553 ,0.178 ,}, {-0.001093 ,0.01329 ,0.15 ,}, }
  //  maen + 4*sigma : {{-0.000423 ,0.001858 ,0.297 ,}, {-3.91e-05 ,-0.001845 ,0.3074 ,}, {-0.0002623 ,7.83e-05 ,0.3054
  //  ,}, {4.447e-05 ,-0.001874 ,0.3074 ,}, {-0.0001752 ,0.003202 ,0.2808 ,}, {-0.0002373 ,0.0007496 ,0.2976 ,}, }

  // ////////////////////////////simulations 3.5 sigma cuts ////////////////////////
  // double mean_minus_3_5_sigma[6][3] = {{-0.00058, 0.00687, 0.19312}, {-0.00088, 0.01022, 0.18360},
  //                                      {-0.00089, 0.00941, 0.18832}, {-0.00066, 0.00888, 0.18466},
  //                                      {-0.00066, 0.00798, 0.18884}, {-0.00055, 0.00685, 0.19319}};

  // double mean_plus_3_5_sigma[6][3] = {{-0.00002, -0.00078, 0.29991}, {0.00023, -0.00396, 0.31026},
  //                                     {0.00010, -0.00156, 0.30077},  {0.00017, -0.00400, 0.31052},
  //                                     {0.00018, -0.00342, 0.30823},  {0.00012, -0.00297, 0.30706}};

  for (Int_t k = 0; k < 6; k++) {
    if (isec == k) {
      upper_lim_total = mean_plus_3_5_sigma[k][0] * pow(_data->p(0), 2) + (mean_plus_3_5_sigma[k][1]) * _data->p(0) +
                        mean_plus_3_5_sigma[k][2];

      lower_lim_total = mean_minus_3_5_sigma[k][0] * pow(_data->p(0), 2) + (mean_minus_3_5_sigma[k][1]) * _data->p(0) +
                        mean_minus_3_5_sigma[k][2];

      // double p0mean[] = {0.112005, 0.113961, 0.111551, 0.114676, 0.112113, 0.112245};
      // double p1mean[] = {-0.103884, -0.0441433, -0.228211, 0.078763, -0.247011, -0.219538};
      // double p2mean[] = {0.00818948, 0.00792426, 0.0111649, 0.0072319, 0.00864546, 0.00970343};
      // double p3mean[] = {-0.000937046, -0.000921012, -0.00131273, -0.000757202, -0.00101063, -0.00122006};

      // double p0sigma[] = {0.0266184, 0.0410105, 0.0246484, 0.0210765, 0.0397108, 0.0200216};
      // double p1sigma[] = {-0.00287631, -0.0163519, -0.00037226, 0.0047182, -0.0149788, 0.00625846};
      // double p2sigma[] = {-0.00354911, -0.00645957, -0.00361653, -0.00303481, -0.00628577, -0.00319593};
      // double p3sigma[] = {0.00026367, 0.000507061, 0.000271306, 0.000234627, 0.000498977, 0.000290974};

      // double sigma_range = 3.5;
      // double mean = 0;
      // double sigma = 0;
      // double upper_lim_total = 0;
      // double lower_lim_total = 0;

      // for (Int_t k = 0; k < 6; k++) {
      //   if (isec == k) {
      //     mean = p0mean[k] * (1 + _data->p(0) / sqrt(_data->p(0) * _data->p(0) + p1mean[k])) + p2mean[k] *
      //     _data->p(0) +
      //            p3mean[k] * _data->p(0) * _data->p(0);
      //     sigma = p0sigma[k] + p1sigma[k] / sqrt(_data->p(0)) + p2sigma[k] * _data->p(0) +
      //             p3sigma[k] * _data->p(0) * _data->p(0);
      //     upper_lim_total = mean + sigma_range * sigma;
      //     lower_lim_total = mean - sigma_range * sigma;
    }
  }

  bool pass_band = _data->ec_tot_energy(0) / _data->p(0) <= upper_lim_total &&
                   _data->ec_tot_energy(0) / _data->p(0) >= lower_lim_total;
  bool pass_triangle = true;

  // if (_data->p(0) < 4.5) {
  //   pass_triangle = true;
  // } else {
  //   pass_triangle = _data->ec_ecin_energy(0) / _data->p(0) > (0.2 - _data->ec_pcal_energy(0) / _data->p(0));
  // }

  if (pass_band && pass_triangle)
    return true;
  else
    return false;
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

bool Pass2_Cuts::EC_hit_position_fiducial_cut_homogeneous() {  //// these are not updated because there is some question
                                                               /// on houw to update

  // Cut using the natural directions of the scintillator bars/ fibers:
  ///////////////////////////////////////////////////////////////////
  /// inbending:
  //
  double min_v_tight_inb[6] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
  double min_v_med_inb[6] = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
  double min_v_loose_inb[6] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
  //
  double max_v_tight_inb[6] = {400, 400, 400, 400, 400, 400};
  double max_v_med_inb[6] = {400, 400, 400, 400, 400, 400};
  double max_v_loose_inb[6] = {400, 400, 400, 400, 400, 400};
  //
  double min_w_tight_inb[6] = {19.0, 19.0, 19.0, 19.0, 19.0, 19.0};
  double min_w_med_inb[6] = {14.0, 14.0, 14.0, 14.0, 14.0, 14.0};
  double min_w_loose_inb[6] = {9.0, 9.0, 9.0, 9.0, 9.0, 9.0};
  //
  double max_w_tight_inb[6] = {400, 400, 400, 400, 400, 400};
  double max_w_med_inb[6] = {400, 400, 400, 400, 400, 400};
  double max_w_loose_inb[6] = {400, 400, 400, 400, 400, 400};

  int isec = (_data->ec_pcal_sec(0) - 1);
  double min_v = min_v_med_inb[isec];
  double max_v = max_v_med_inb[isec];
  double min_w = min_w_med_inb[isec];
  double max_w = max_w_med_inb[isec];
  return (_data->ec_pcal_lv(0) > min_v && _data->ec_pcal_lv(0) < max_v && _data->ec_pcal_lw(0) > min_w &&
          _data->ec_pcal_lw(0) < max_w);
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

bool Pass2_Cuts::PCAL_fiducial_cut_HX_HY() {
  double minparams_pcal_in[6][2] = {{-0.52452, 20.33242}, {-0.51548, 18.38758}, {-0.49609, 19.04455},
                                    {-0.51318, 22.13909}, {-0.50361, 20.48697}, {-0.51821, 19.48394}};

  double maxparams_pcal_in[6][2] = {{0.52494, -20.38030}, {0.50706, -22.01970}, {0.50900, -21.77000},
                                    {0.51967, -19.31667}, {0.52082, -23.08091}, {0.52288, -20.65061}};

  short pcal_sector = (_data->ec_pcal_sec(0) - 1);

  double HX = _data->ec_pcal_hx(0);
  double HY = _data->ec_pcal_hy(0);

  float HX_new = HX * cos(DEG2RAD * (-60 * (pcal_sector))) - HY * sin(DEG2RAD * (-60 * (pcal_sector)));
  HY = HX * sin(DEG2RAD * (-60 * (pcal_sector))) + HY * cos(DEG2RAD * (-60 * (pcal_sector)));

  HX = HX_new;

  double calc_min = minparams_pcal_in[pcal_sector][0] * HX + minparams_pcal_in[pcal_sector][1];
  double calc_max = maxparams_pcal_in[pcal_sector][0] * HX + maxparams_pcal_in[pcal_sector][1];

  return ((HY > calc_min) && (HY < calc_max));
}
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

bool Pass2_Cuts::DC_fiducial_cut_XY() {
  // bool _dc_fid_cut = true;
  // bool isinbending = true;
  // new cut parameters for the linear cut based on x and y coordinates (inbending field):
  // replace it in the function: bool DC_fiducial_cut_XY(int j, int region)
  // (optimized for electrons, do not use it for hadrons)

  // maxparams_in[1][6][3][2] -> [pid][sec][regions][a,b]->a*x+b

  double minparams_in[1][6][3][2] = {{{{-0.56964, 11.63393}, {-0.58683, 18.99917}, {-0.56401, 21.69753}},
                                      {{-0.56179, 10.35179}, {-0.58000, 17.43889}, {-0.55588, 20.13489}},
                                      {{-0.52107, 9.64821}, {-0.57333, 19.19444}, {-0.55379, 22.50467}},
                                      {{-0.56571, 12.33571}, {-0.58350, 20.34917}, {-0.56077, 23.63846}},
                                      {{-0.53714, 10.55000}, {-0.57667, 19.56667}, {-0.54082, 21.52995}},
                                      {{-0.55929, 10.78929}, {-0.57817, 17.77583}, {-0.54918, 19.78929}}}};

  double maxparams_in[1][6][3][2] = {{{{0.56679, -10.93393}, {0.58683, -18.22139}, {0.56231, -21.30000}},
                                      {{0.54964, -11.43393}, {0.58650, -20.73417}, {0.53275, -21.52033}},
                                      {{0.55536, -11.49107}, {0.57667, -19.56667}, {0.54588, -22.23874}},
                                      {{0.57000, -11.27143}, {0.58983, -18.60639}, {0.56038, -21.00577}},
                                      {{0.56893, -12.67321}, {0.58533, -20.64556}, {0.56434, -24.17720}},
                                      {{0.57536, -11.69107}, {0.58683, -18.99917}, {0.56857, -22.52143}}}};

  // // BE CAREFUL HERE

  int pid = 0;
  short dc_sector = (_data->dc_sec(0) - 1);
  // float sin_dc_sec = sinf(-dc_sector * ROTATE);
  // float cos_dc_sec = cosf(-dc_sector * ROTATE);

  double X1 = _data->dc_r1_x(0);
  double Y1 = _data->dc_r1_y(0);
  // double X1_new =
  //         X1 * cos_dc_sec - Y1 * sin_dc_sec;
  // Y1 = X1 * sin_dc_sec + Y1 * cos_dc_sec;

  float X1_new = X1 * cos(DEG2RAD * (-60 * (dc_sector))) - Y1 * sin(DEG2RAD * (-60 * (dc_sector)));
  Y1 = X1 * sin(DEG2RAD * (-60 * (dc_sector))) + Y1 * cos(DEG2RAD * (-60 * (dc_sector)));

  X1 = X1_new;
  int region_1 = 1;

  double calc_min1 = minparams_in[pid][dc_sector][region_1 - 1][0] * X1 + minparams_in[pid][dc_sector][region_1 - 1][1];
  double calc_max1 = maxparams_in[pid][dc_sector][region_1 - 1][0] * X1 + maxparams_in[pid][dc_sector][region_1 - 1][1];
  // _dc_fid_cut &= (Y1 > calc_min1);
  // _dc_fid_cut &= (Y1 < calc_max1);

  double X2 = _data->dc_r2_x(0);
  double Y2 = _data->dc_r2_y(0);
  // double X2_new =
  //         X2 * cos_dc_sec - Y2 * sin_dc_sec;
  // Y2 = X2 * sin_dc_sec + Y2 * cos_dc_sec;

  float X2_new = X2 * cos(DEG2RAD * (-60 * (dc_sector))) - Y2 * sin(DEG2RAD * (-60 * (dc_sector)));
  Y2 = X2 * sin(DEG2RAD * (-60 * (dc_sector))) + Y2 * cos(DEG2RAD * (-60 * (dc_sector)));

  X2 = X2_new;
  int region_2 = 2;

  double calc_min2 = minparams_in[pid][dc_sector][region_2 - 1][0] * X2 + minparams_in[pid][dc_sector][region_2 - 1][1];
  double calc_max2 = maxparams_in[pid][dc_sector][region_2 - 1][0] * X2 + maxparams_in[pid][dc_sector][region_2 - 1][1];
  // _dc_fid_cut &= (Y2 > calc_min2);
  // _dc_fid_cut &= (Y2 < calc_max2);

  double X3 = _data->dc_r3_x(0);
  double Y3 = _data->dc_r3_y(0);
  // double X3_new =
  //         X3 * cos_dc_sec - Y3 * sin_dc_sec;
  // Y3 = X3 * sin_dc_sec + Y3 * cos_dc_sec;

  float X3_new = X3 * cos(DEG2RAD * (-60 * (dc_sector))) - Y3 * sin(DEG2RAD * (-60 * (dc_sector)));
  Y3 = X3 * sin(DEG2RAD * (-60 * (dc_sector))) + Y3 * cos(DEG2RAD * (-60 * (dc_sector)));

  X3 = X3_new;
  int region_3 = 3;

  double calc_min3 = minparams_in[pid][dc_sector][region_3 - 1][0] * X3 + minparams_in[pid][dc_sector][region_3 - 1][1];
  double calc_max3 = maxparams_in[pid][dc_sector][region_3 - 1][0] * X3 + maxparams_in[pid][dc_sector][region_3 - 1][1];
  // _dc_fid_cut &= (Y3 > calc_min3);
  // _dc_fid_cut &= (Y3 < calc_max3);
  // std::cout << "y2 " << Y2  << " calc_max2  " << calc_max2 <<'\n';
  // std::cout << "y3 " << Y3  << " calc_max2  " << calc_max3 <<'\n';

  return ((Y1 > calc_min1) && (Y1 < calc_max1) && (Y2 > calc_min2) && (Y2 < calc_max2) && (Y3 > calc_min3) &&
          (Y3 < calc_max3));
  // _dc_fid_cut &= (Y1 < calc_max1);
  // return _dc_fid_cut;

  // if (!_dc_fid_cut)
  //         return _dc_fid_cut;

  // if (!_dc_fid_cut)
  //         return _dc_fid_cut;
  //

  // if (!_dc_fid_cut)
  //         return _dc_fid_cut;
  // int pid = 0;
  // switch (_data->pid(i)) {
  // case 11:
  //         pid = 0;
  //         break;
  // case 2212:
  //         pid = 1;
  //         break;
  // case 211:
  //         pid = 2;
  //         break;
  // case -211:
  //         pid = 3;
  //         break;
  // case 321:
  //         pid = 4;
  //         break;
  // case -321:
  //         pid = 5;
  //         break;
  // default:
  //         return false;
  // }
  // if(inbending == true) pid = 0; // use only for electrons in inbending case
  // double calc_min = minparams[pid][dc_sector - 1][region - 1][0] + minparams[pid][dc_sector - 1][region - 1][1] *
  // X; double calc_max = maxparams[pid][dc_sector - 1][region - 1][0] + maxparams[pid][dc_sector - 1][region -
  // 1][1] * X; return (Y > calc_min) && (Y < calc_max);
}

bool Pass2_Cuts::DC_z_vertex_cut() {
  int pcal_sector = _data->ec_pcal_sec(0);
  float partvz = _data->vz(0);

  double vz_min_sect_inb[] = {-10, -10, -10, -10, -10, -10};
  double vz_max_sect_inb[] = {5, 5, 5, 5, 5, 5};

  float vz_min_sect[6];
  float vz_max_sect[6];

  for (int i = 0; i < 6; i++) {
    vz_min_sect[i] = vz_min_sect_inb[i];
    vz_max_sect[i] = vz_max_sect_inb[i];
  }

  int isec = pcal_sector - 1;
  float vz_min = vz_min_sect[isec];
  float vz_max = vz_max_sect[isec];

  return partvz > vz_min && partvz < vz_max;
}

// public class HadronCuts {
//
// /**
//  * DC fiducial cut for hadrons
//  * @param dc_sector sector of hits in DC
//  * @param region specify fiducial Pass2_Cuts for which region to use
//  * @param trajx x for region 1 or 2 or 3 from REC::Traj
//  * @param trajy y for region 1 or 2 or 3 from REC::Traj
//  * @param trajz z for region 1 or 2 or 3 from REC::Traj
//  * @param partpid pid assigned to particle candidate
//  * @param isinbending True if magnetic field is inbending
//  */
bool Pass2_Cuts::DC_fiducial_cut_theta_phi(int i) {  ///////////////// use this and vs above xy cuts and compare.
  // new cut parameters for the polynomial cut based on the local theta and phi coordinates (inbending field):
  // replace it in the function: bool DC_fiducial_cut_theta_phi(int j, int region)
  // (optimized for pi+ and pi-, not optimized for Kaons yet)
  //
  short dc_sector = (_data->dc_sec(i));  //_data->dc_sec(i) ??

  float trajx1 = _data->dc_r1_x(i);
  float trajy1 = _data->dc_r1_y(i);
  float trajz1 = _data->dc_r1_z(i);

  float trajx2 = _data->dc_r2_x(i);
  float trajy2 = _data->dc_r2_y(i);
  float trajz2 = _data->dc_r2_z(i);

  float trajx3 = _data->dc_r3_x(i);
  float trajy3 = _data->dc_r3_y(i);
  float trajz3 = _data->dc_r3_z(i);

  int partpid = _data->pid(i);
  bool isinbending = true;

  double maxparams_in[6][6][3][4] = {{{{-37.5489, 27.4543, -1.11484, 0.00522935},
                                       {-29.7228, 26.7512, -1.52592, 0.0122397},
                                       {-20.3559, 23.1586, -1.47441, 0.0133898}},
                                      {{-36.2719, 25.1427, -0.817973, 0.00233912},
                                       {-28.2118, 25.0664, -1.29748, 0.00947493},
                                       {-20.6015, 22.9639, -1.39759, 0.012069}},
                                      {{-34.1013, 25.9343, -1.23555, 0.00959955},
                                       {-24.0285, 22.9346, -1.165, 0.00846331},
                                       {-8.04969, 12.5436, -0.268326, 9.03561e-11}},
                                      {{-48.5546, 36.1076, -2.07362, 0.0161268},
                                       {-24.7284, 22.9355, -1.12754, 0.00796403},
                                       {-22.5292, 24.1624, -1.52361, 0.0137042}},
                                      {{-40.4295, 30.8386, -1.77195, 0.0156563},
                                       {-26.7149, 23.5322, -1.1011, 0.00715825},
                                       {-10.9822, 13.8127, -0.312534, 1.32292e-05}},
                                      {{-38.1396, 28.0524, -1.19166, 0.00613986},
                                       {-26.1238, 24.3235, -1.28254, 0.00950751},
                                       {-19.0376, 22.042, -1.32482, 0.0113948}}},
                                     {{{-1.67037e-08, 12.8334, -0.820443, 0.00818882},
                                       {-6.23823, 14.8659, -0.776403, 0.00624484},
                                       {-5.75713, 11.4787, -0.227124, 6.61281e-10}},
                                      {{-6.09637e-07, 12.7972, -0.813133, 0.00808401},
                                       {-5.51055, 13.9682, -0.639287, 0.00441366},
                                       {-7.90046, 12.5383, -0.271117, 1.86929e-10}},
                                      {{-2.84217e-14, 13.0836, -0.864047, 0.00869759},
                                       {-6.78639, 15.3367, -0.827197, 0.00677168},
                                       {-4.8928, 11.1884, -0.221965, 1.51263e-10}},
                                      {{-3.8595e-09, 12.9673, -0.841224, 0.0083938},
                                       {-4.01784, 12.9989, -0.557548, 0.00367493},
                                       {-1.95023, 9.69687, -0.157901, 5.33239e-09}},
                                      {{-6.43496e-10, 12.9804, -0.850651, 0.00863353},
                                       {-5.10299, 13.9958, -0.671087, 0.00489619},
                                       {-6.03313, 11.7973, -0.249435, 1.2754e-11}},
                                      {{-2.94932e-10, 13.1054, -0.859032, 0.00848181},
                                       {-6.05945, 14.7331, -0.742818, 0.00558374},
                                       {-5.63811, 11.6686, -0.247509, 2.33147e-13}}},
                                     {{{-2.68279e-07, 12.99, -0.846226, 0.00845788},
                                       {-14.6317, 19.3874, -1.09244, 0.00899541},
                                       {-38.1915, 29.8688, -1.59229, 0.0120089}},
                                      {{-0.996514, 13.9379, -0.964686, 0.00982941},
                                       {-15.9613, 20.2461, -1.16106, 0.00955431},
                                       {-35.9455, 29.0996, -1.586, 0.0122175}},
                                      {{-1.14284e-07, 13.6015, -0.966952, 0.0101523},
                                       {-15.5288, 20.3045, -1.20523, 0.0102808},
                                       {-34.2682, 26.4216, -1.20609, 0.0078434}},
                                      {{-1.70075e-08, 13.0005, -0.832325, 0.00817159},
                                       {-7.66776, 15.4526, -0.779727, 0.00585967},
                                       {-26.8035, 23.9995, -1.2322, 0.00942061}},
                                      {{-9.53804e-10, 13.2563, -0.898206, 0.00917629},
                                       {-6.85083, 14.8485, -0.722803, 0.0053221},
                                       {-39.3606, 31.5412, -1.83015, 0.0148302}},
                                      {{-7.66835e-07, 13.937, -1.05153, 0.0118223},
                                       {-9.7913, 16.925, -0.913158, 0.00712552},
                                       {-27.722, 23.9412, -1.1314, 0.00761088}}},
                                     {{{-22.1832, 20.4134, -0.764848, 0.00310923},
                                       {-31.0844, 28.2369, -1.715, 0.0145145},
                                       {-9.52175, 18.7932, -1.38896, 0.0150233}},
                                      {{-21.5849, 20.2457, -0.762109, 0.00305359},
                                       {-19.5601, 21.5945, -1.18955, 0.00939109},
                                       {-1.57084, 13.3989, -0.823161, 0.00795227}},
                                      {{-16.052, 16.6264, -0.444308, 2.82701e-06},
                                       {-13.8291, 18.6541, -1.01549, 0.00825776},
                                       {-1.92223e-05, 13.0305, -0.881089, 0.00925281}},
                                      {{-19.821, 18.4301, -0.516168, 2.17199e-10},
                                       {-30.6295, 28.0989, -1.71897, 0.0146585},
                                       {-9.23709, 17.1589, -1.03955, 0.00943673}},
                                      {{-16.1795, 16.7121, -0.448883, 1.53774e-11},
                                       {-23.6418, 24.5748, -1.48652, 0.01254},
                                       {-4.2626e-09, 12.899, -0.845374, 0.00872171}},
                                      {{-9.74791, 15.0287, -0.531727, 0.00192371},
                                       {-41.0848, 33.1802, -1.97671, 0.0158148},
                                       {-4.12428, 14.3361, -0.820483, 0.00725632}}},
                                     {{{-1.05499e-08, 12.7347, -0.800158, 0.00789171},
                                       {-3.78358, 13.3272, -0.620589, 0.0043452},
                                       {-31.0947, 26.2276, -1.33783, 0.00961276}},
                                      {{-3.20108e-05, 13.2084, -0.89232, 0.00907651},
                                       {-11.5913, 18.4403, -1.08132, 0.00895511},
                                       {-26.4998, 23.4434, -1.09015, 0.00695521}},
                                      {{-1.54979e-07, 13.3849, -0.912541, 0.00919697},
                                       {-4.77271, 14.366, -0.750675, 0.00582608},
                                       {-31.7881, 27.2978, -1.49603, 0.0115217}},
                                      {{-8.46957e-07, 13.135, -0.863007, 0.00850261},
                                       {-5.91254, 14.7345, -0.748863, 0.00564354},
                                       {-27.2818, 24.4544, -1.24541, 0.009006}},
                                      {{-8.97242e-09, 12.8923, -0.825914, 0.00815967},
                                       {-6.91507, 16.0014, -0.917916, 0.00756705},
                                       {-18.1359, 18.5543, -0.695074, 0.00311518}},
                                      {{-2.50141e-08, 13.1356, -0.864227, 0.00854005},
                                       {-6.62648, 15.5703, -0.861224, 0.00697927},
                                       {-19.9356, 18.969, -0.647219, 0.00209364}}},
                                     {{{-31.056, 26.1595, -1.20596, 0.00643836},
                                       {-44.4944, 36.2986, -2.35276, 0.020162},
                                       {-12.2855, 21.0109, -1.61628, 0.0172125}},
                                      {{-27.3898, 25.1282, -1.2366, 0.00728902},
                                       {-24.9794, 23.2357, -1.09342, 0.00656412},
                                       {-16.9519, 23.8236, -1.78734, 0.017541}},
                                      {{-28.7906, 26.9219, -1.49542, 0.0104976},
                                       {-22.0922, 23.6046, -1.37835, 0.0110503},
                                       {-5.24383, 16.5267, -1.15701, 0.0113067}},
                                      {{-3.92728, 12.0692, -0.372323, 0.0011559},
                                       {-23.5702, 22.3459, -1.04378, 0.00649998},
                                       {-17.3561, 24.4119, -1.93535, 0.0204532}},
                                      {{-30.442, 26.0012, -1.2191, 0.00674908},
                                       {-54.5014, 42.354, -2.8256, 0.0242569},
                                       {-0.751452, 13.9234, -0.958253, 0.00952713}},
                                      {{-31.216, 26.1169, -1.20087, 0.00650951},
                                       {-31.0314, 28.4075, -1.70479, 0.0137299},
                                       {-13.8981, 22.326, -1.72999, 0.0176742}}}};

  double minparams_in[6][6][3][4] = {{{{45.6964, -33.9555, 1.83632, -0.0133721},
                                       {16.3132, -19.1709, 0.95922, -0.00719164},
                                       {17.4745, -21.3091, 1.29658, -0.0114378}},
                                      {{34.063, -25.5129, 0.992129, -0.00445872},
                                       {22.4188, -23.1898, 1.33328, -0.011079},
                                       {15.558, -20.779, 1.32969, -0.0122892}},
                                      {{28.8399, -21.4732, 0.662977, -0.00227941},
                                       {15.2776, -18.4944, 0.917128, -0.00703012},
                                       {25.9277, -26.2555, 1.70407, -0.0154587}},
                                      {{43.4091, -32.329, 1.78095, -0.0143066},
                                       {34.8052, -27.7186, 1.43403, -0.0108989},
                                       {26.384, -24.813, 1.4364, -0.0123938}},
                                      {{42.094, -32.8674, 2.12321, -0.0208007},
                                       {39.6248, -33.4591, 2.1938, -0.0196953},
                                       {17.5854, -17.6921, 0.617536, -0.00282672}},
                                      {{24.4957, -19.3118, 0.481099, -6.0729e-07},
                                       {22.7714, -23.2117, 1.31478, -0.0107808},
                                       {16.2955, -21.0448, 1.33876, -0.0123879}}},
                                     {{{2.01913e-05, -13.2206, 0.868885, -0.00845047},
                                       {6.86331, -15.0105, 0.765473, -0.00602765},
                                       {5.15884, -11.18, 0.21433, -1.79763e-09}},
                                      {{3.24593, -15.5188, 1.12128, -0.011555},
                                       {8.61633, -16.3281, 0.913374, -0.00783236},
                                       {4.51456, -11.0507, 0.243113, -0.000607925}},
                                      {{0.905676, -13.3623, 0.85485, -0.00835569},
                                       {6.87062, -14.5731, 0.694399, -0.00526577},
                                       {3.8283, -10.4277, 0.178245, -4.2334e-10}},
                                      {{5.54817e-07, -12.6609, 0.744683, -0.00664861},
                                       {6.25817, -14.6969, 0.728253, -0.00543273},
                                       {6.01169, -11.8105, 0.251251, -3.71394e-10}},
                                      {{9.30801e-09, -13.3207, 0.888792, -0.00873133},
                                       {8.41797, -16.4985, 0.956897, -0.00841779},
                                       {4.36256, -10.8341, 0.202655, -3.44186e-09}},
                                      {{0.27863, -13.1208, 0.833431, -0.0079631},
                                       {7.38412, -15.4188, 0.82054, -0.00681735},
                                       {4.48567, -10.7376, 0.190611, -9.77392e-10}}},
                                     {{{1.59369e-06, -13.8294, 0.990918, -0.0103128},
                                       {20.1273, -23.853, 1.58449, -0.0145959},
                                       {40.8152, -32.8944, 2.00731, -0.0171007}},
                                      {{1.4334, -14.5452, 1.04379, -0.0106791},
                                       {19.9242, -23.3894, 1.5036, -0.0134429},
                                       {45.1348, -34.9897, 2.11238, -0.0175613}},
                                      {{4.48276e-06, -12.6688, 0.757818, -0.006981},
                                       {10.2525, -16.9056, 0.909637, -0.00739798},
                                       {33.2958, -27.7763, 1.53467, -0.0123488}},
                                      {{3.817e-06, -13.2285, 0.856439, -0.0081744},
                                       {12.5356, -19.0801, 1.1686, -0.0102758},
                                       {37.3388, -29.7344, 1.64296, -0.0130658}},
                                      {{3.64842e-07, -14.1631, 1.0771, -0.0118569},
                                       {9.85442, -17.8198, 1.12641, -0.010627},
                                       {34.7, -28.5335, 1.57226, -0.0124004}},
                                      {{0.828721, -13.6429, 0.895665, -0.00866683},
                                       {10.8176, -18.0919, 1.11147, -0.010183},
                                       {29.9288, -24.3389, 1.08973, -0.00703934}}},
                                     {{{15.8302, -16.9632, 0.53561, -0.00136216},
                                       {32.8002, -29.2569, 1.79783, -0.015324},
                                       {1.98393, -13.0099, 0.70788, -0.00615153}},
                                      {{16.0367, -16.5901, 0.470678, -0.000728065},
                                       {32.4005, -29.7403, 1.92286, -0.0171968},
                                       {2.39707, -13.6612, 0.816883, -0.00770837}},
                                      {{22.0623, -21.6319, 1.02811, -0.00680893},
                                       {32.7467, -29.6099, 1.87839, -0.0164223},
                                       {1.19902e-08, -12.972, 0.863127, -0.00884759}},
                                      {{21.5883, -21.198, 0.957819, -0.00575361},
                                       {25.7387, -25.4963, 1.5428, -0.0131855},
                                       {6.06479, -16.6311, 1.16092, -0.0117194}},
                                      {{19.6915, -19.1751, 0.704086, -0.00288768},
                                       {28.6596, -27.3351, 1.70309, -0.0148193},
                                       {5.30096e-08, -11.8562, 0.621373, -0.00541869}},
                                      {{20.6594, -19.8704, 0.786033, -0.00394155},
                                       {20.7612, -22.3774, 1.27116, -0.0104109},
                                       {2.56196, -14.4159, 0.98009, -0.0100214}}},
                                     {{{6.84429e-08, -11.7778, 0.558372, -0.00403519},
                                       {5.88119, -14.1561, 0.630592, -0.00400605},
                                       {22.9399, -21.6066, 0.97379, -0.00604844}},
                                      {{5.49686, -16.3382, 1.10037, -0.0105049},
                                       {9.25791, -16.8955, 0.947447, -0.00774283},
                                       {19.4826, -18.4694, 0.587601, -0.00147216}},
                                      {{0.148482, -12.4191, 0.691879, -0.00595948},
                                       {6.95863, -15.5624, 0.862069, -0.00725014},
                                       {16.6631, -16.746, 0.461105, -0.000520762}},
                                      {{2.64705e-10, -11.8828, 0.574528, -0.00419463},
                                       {5.45746, -13.9134, 0.602948, -0.00360009},
                                       {31.3252, -27.342, 1.51348, -0.0115756}},
                                      {{3.46769, -15.3338, 1.02031, -0.00951104},
                                       {0.368693, -11.8657, 0.574108, -0.0044343},
                                       {39.7957, -32.8529, 2.02652, -0.016978}},
                                      {{0.00207118, -12.0447, 0.602167, -0.00447581},
                                       {3.03476, -12.9176, 0.603586, -0.00440659},
                                       {32.0315, -26.8451, 1.37417, -0.00966969}}},
                                     {{{56.9355, -42.3826, 2.61014, -0.0202986},
                                       {28.8989, -27.1772, 1.63996, -0.0136625},
                                       {4.30155, -15.1455, 0.995784, -0.0100192}},
                                      {{13.4916, -17.1287, 0.681434, -0.0031646},
                                       {32.246, -29.0499, 1.77696, -0.0148718},
                                       {2.22052, -9.65178, 0.133616, -9.0964e-05}},
                                      {{41.8686, -33.5132, 1.92542, -0.0142307},
                                       {0.0645903, -9.74163, 0.217245, -2.22987e-05},
                                       {9.58895e-09, -13.2013, 0.926579, -0.00993616}},
                                      {{34.8087, -28.1804, 1.3547, -0.00784213},
                                       {31.3059, -28.7057, 1.76134, -0.0146575},
                                       {8.66833, -17.8896, 1.20937, -0.0116248}},
                                      {{42.0802, -33.525, 1.91492, -0.0140721},
                                       {36.8805, -31.3893, 1.91131, -0.0157056},
                                       {6.11008, -17.0626, 1.24276, -0.0127673}},
                                      {{39.6762, -31.6354, 1.73354, -0.0123964},
                                       {30.2451, -27.8243, 1.67413, -0.0138583},
                                       {4.78902, -14.9558, 0.912758, -0.00855026}}}};

  double trajr1 = sqrt(pow(trajx1, 2) + pow(trajy1, 2) + pow(trajz1, 2));
  double theta_DCr1 = RAD2DEG * (acos(trajz1 / trajr1));
  double phi_DCr_raw1 = RAD2DEG * (atan2(trajy1 / trajr1, trajx1 / trajr1));

  // std::cout << " trajx1 " << trajx1 << "  trajy1 " << trajy1  <<"  trajz1 "<<trajz1 << "  r1  "<<trajr1<<'\n';
  // std::cout << " acos  " << acos(trajz1/trajr1) <<'\n';
  // std::cout << "  atan2  " <<atan2(trajy1/trajr1, trajx1/trajr1)<< '\n';
  // std::cout << "theta " << theta_DCr1 << '\n';

  double phi_DCr1 = 5000;

  if (dc_sector == 1) phi_DCr1 = phi_DCr_raw1;
  if (dc_sector == 2) phi_DCr1 = phi_DCr_raw1 - 60;
  if (dc_sector == 3) phi_DCr1 = phi_DCr_raw1 - 120;
  if (dc_sector == 4 && phi_DCr_raw1 > 0) phi_DCr1 = phi_DCr_raw1 - 180;
  if (dc_sector == 4 && phi_DCr_raw1 < 0) phi_DCr1 = phi_DCr_raw1 + 180;
  if (dc_sector == 5) phi_DCr1 = phi_DCr_raw1 + 120;
  if (dc_sector == 6) phi_DCr1 = phi_DCr_raw1 + 60;

  // std::cout << "phi " << phi_DCr_raw1<< '\n';

  double trajr2 = sqrt(pow(trajx2, 2) + pow(trajy2, 2) + pow(trajz2, 2));
  double theta_DCr2 = RAD2DEG * (acos(trajz2 / trajr2));
  double phi_DCr_raw2 = RAD2DEG * (atan2(trajy2 / trajr2, trajx2 / trajr2));

  double phi_DCr2 = 5000;

  if (dc_sector == 1) phi_DCr2 = phi_DCr_raw2;
  if (dc_sector == 2) phi_DCr2 = phi_DCr_raw2 - 60;
  if (dc_sector == 3) phi_DCr2 = phi_DCr_raw2 - 120;
  if (dc_sector == 4 && phi_DCr_raw2 > 0) phi_DCr2 = phi_DCr_raw2 - 180;
  if (dc_sector == 4 && phi_DCr_raw2 < 0) phi_DCr2 = phi_DCr_raw2 + 180;
  if (dc_sector == 5) phi_DCr2 = phi_DCr_raw2 + 120;
  if (dc_sector == 6) phi_DCr2 = phi_DCr_raw2 + 60;

  double trajr3 = sqrt(pow(trajx3, 2) + pow(trajy3, 2) + pow(trajz3, 2));
  double theta_DCr3 = RAD2DEG * (acos(trajz3 / trajr3));
  double phi_DCr_raw3 = RAD2DEG * (atan2(trajy3 / trajr3, trajx3 / trajr3));

  double phi_DCr3 = 5000;

  if (dc_sector == 1) phi_DCr3 = phi_DCr_raw3;
  if (dc_sector == 2) phi_DCr3 = phi_DCr_raw3 - 60;
  if (dc_sector == 3) phi_DCr3 = phi_DCr_raw3 - 120;
  if (dc_sector == 4 && phi_DCr_raw3 > 0) phi_DCr3 = phi_DCr_raw3 - 180;
  if (dc_sector == 4 && phi_DCr_raw3 < 0) phi_DCr3 = phi_DCr_raw3 + 180;
  if (dc_sector == 5) phi_DCr3 = phi_DCr_raw3 + 120;
  if (dc_sector == 6) phi_DCr3 = phi_DCr_raw3 + 60;

  int pid = 0;

  switch (partpid) {
    case 11:
      pid = 0;
      break;
    case 2212:
      pid = 1;
      break;
    case 211:
      pid = 2;
      break;
    case -211:
      pid = 3;
      break;
    case 321:
      pid = 4;
      break;
    case -321:
      pid = 5;
      break;
    default:
      return false;
  }
  int region1 = 1;
  int region2 = 2;
  int region3 = 3;
  // std::cout << "pid = " <<pid<< '\n';
  double calc_phi_min1 = minparams_in[pid][dc_sector - 1][region1 - 1][0] +
                         minparams_in[pid][dc_sector - 1][region1 - 1][1] * log(theta_DCr1) +
                         minparams_in[pid][dc_sector - 1][region1 - 1][2] * theta_DCr1 +
                         minparams_in[pid][dc_sector - 1][region1 - 1][3] * theta_DCr1 * theta_DCr1;

  double calc_phi_max1 = maxparams_in[pid][dc_sector - 1][region1 - 1][0] +
                         maxparams_in[pid][dc_sector - 1][region1 - 1][1] * log(theta_DCr1) +
                         maxparams_in[pid][dc_sector - 1][region1 - 1][2] * theta_DCr1 +
                         maxparams_in[pid][dc_sector - 1][region1 - 1][3] * theta_DCr1 * theta_DCr1;
  // std::cout << "  phi dcr1 " <<phi_DCr1 << '\n';
  //
  // std::cout << "phi_min " <<calc_phi_min1  << "  phi_max " <<calc_phi_max1<< '\n';
  // std::cout << "log 10 " <<log(10)<< '\n';
  // std::cout << "calc_phi_min " << calc_phi_min1 <<'\n';
  // std::cout << "calc_phi_max " << calc_phi_max1 <<'\n';

  // return (phi_DCr1 > calc_phi_min1) && (phi_DCr1 < calc_phi_max1);
  //
  double calc_phi_min2 = minparams_in[pid][dc_sector - 1][region2 - 1][0] +
                         minparams_in[pid][dc_sector - 1][region2 - 1][1] * log(theta_DCr2) +
                         minparams_in[pid][dc_sector - 1][region2 - 1][2] * theta_DCr2 +
                         minparams_in[pid][dc_sector - 1][region2 - 1][3] * theta_DCr2 * theta_DCr2;

  double calc_phi_max2 = maxparams_in[pid][dc_sector - 1][region2 - 1][0] +
                         maxparams_in[pid][dc_sector - 1][region2 - 1][1] * log(theta_DCr2) +
                         maxparams_in[pid][dc_sector - 1][region2 - 1][2] * theta_DCr2 +
                         maxparams_in[pid][dc_sector - 1][region2 - 1][3] * theta_DCr2 * theta_DCr2;

  // return (phi_DCr2 > calc_phi_min2) && (phi_DCr2 < calc_phi_max2);

  double calc_phi_min3 = minparams_in[pid][dc_sector - 1][region3 - 1][0] +
                         minparams_in[pid][dc_sector - 1][region3 - 1][1] * log(theta_DCr3) +
                         minparams_in[pid][dc_sector - 1][region3 - 1][2] * theta_DCr3 +
                         minparams_in[pid][dc_sector - 1][region3 - 1][3] * theta_DCr3 * theta_DCr3;

  double calc_phi_max3 = maxparams_in[pid][dc_sector - 1][region3 - 1][0] +
                         maxparams_in[pid][dc_sector - 1][region3 - 1][1] * log(theta_DCr3) +
                         maxparams_in[pid][dc_sector - 1][region3 - 1][2] * theta_DCr3 +
                         maxparams_in[pid][dc_sector - 1][region3 - 1][3] * theta_DCr3 * theta_DCr3;

  return ((phi_DCr1 > calc_phi_min1) && (phi_DCr1 < calc_phi_max1) && (phi_DCr2 > calc_phi_min2) &&
          (phi_DCr2 < calc_phi_max2) && (phi_DCr3 > calc_phi_min3) && (phi_DCr3 < calc_phi_max3));
}

/** Delta VZ cut for hadrons
 * @param pid hadron PID code
 * @param dvz difference between Vz of hadron candidate and electron
 */
bool Pass2_Cuts::Hadron_Delta_vz_cut(int i) {
  int pid = _data->pid(i);
  // if(pid==PROTON){
  float dvz = (_data->vz(i) - _data->vz(0));

  // std::cout<<"dvz  "<<dvz<<std::endl;

  // return dvz > -20 && dvz < 20;}
  switch (pid) {
    case 2212:
      return dvz > -20 && dvz < 20;
    case 22:
      return dvz > -20 && dvz < 20;
    case 2112:
      return dvz > -20 && dvz < 20;
    case 211:
      return dvz > -20 && dvz < 20;
    case -211:
      return dvz > -20 && dvz < 20;
    case 321:
      return dvz > -20 && dvz < 20;
    case -321:
      return dvz > -20 && dvz < 20;
  }
  return false;
}

/** chi2pid cut for hadrons
 * @param chi2pid chi2pid value
 * @param pid hadron PID code
 */
bool Pass2_Cuts::Hadron_Chi2pid_cut(int i) {
  bool isstrict = false;
  float chi2pid = _data->chi2pid(i);
  float p = _data->p(i);
  int pid = _data->pid(i);
  int status = abs(_data->status(i));

  double coef;
  if (pid == 211)
    coef = 0.88;
  else if (pid == -211)
    coef = 0.93;

  else if (pid == 2212) {
    if (status < 4000)
      return abs(chi2pid) < 3.0;  /// please confirm this first 2.64 is given for rga fall 2018
    else {
      return abs(chi2pid) < 6.0;
    }
  }

  else
    return false;

  bool chi2cut = false;
  if (status < 4000) {
    if (isstrict) {
      if (p < 2.44)
        chi2cut = chi2pid < 3 * coef;
      else if (p < 4.6)
        chi2cut = chi2pid < coef * (0.00869 + 14.98587 * exp(-p / 1.18236) + 1.81751 * exp(-p / 4.86394));
      else
        chi2cut = chi2pid < coef * (-1.14099 + 24.14992 * exp(-p / 1.36554) + 2.66876 * exp(-p / 6.80522));
    } else {
      if (p < 2.44)
        chi2cut = chi2pid < 3 * coef;
      else
        chi2cut = chi2pid < coef * (0.00869 + 14.98587 * exp(-p / 1.18236) + 1.81751 * exp(-p / 4.86394));
    }

    return chi2cut && chi2pid > coef * -3;
  } else {
    return abs(chi2pid) < 6.0;
  }
}

//}

///////////////////// Pass2_Cuts ///////////////////////
