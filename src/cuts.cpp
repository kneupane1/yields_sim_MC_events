/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
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
bool Cuts::IsPip(int i) {
  if (_data->gpart() <= i) return false;
  bool _pip = true;
  //   _pip &= (_data->charge(i) == POSITIVE);
  _pip &= (_data->pid(i) == PIP);
  _pip &= (abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) < 0.4);
  // _pip &= !(abs(_dt->dt_P(i)) < 0.5 || abs(_dt->dt_ctof_P(i)) < 0.2);
  //   _pip &= (_data->p(i) > 0.2);
  _pip &= (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000);

  // // min/max mom cuts
  // if (abs(_data->status(i)) < 4000) {
  //   _pip &= (_data->p(i) > 0.3);
  //   _pip &= (_data->p(i) < 5.0);
  // } else if (abs(_data->status(i)) >= 4000) {
  //   _pip &= (_data->p(i) > 0.0);
  //   _pip &= (_data->p(i) < 3.0);
  // }
  ////_pip &= (abs(_data->chi2pid(i)) < 0.5);

  //_pip &= DC_fiducial_cut_theta_phi(i);
  //_pip &= Hadron_Delta_vz_cut(i);
  //_pip&= Hadron_Chi2pid_cut(i);

  return _pip;
}
bool Cuts::IsProton(int i) {
  if (_data->gpart() <= i) return false;
  bool _proton = true;
  //   _proton &= (_data->charge(i) == POSITIVE);
  _proton &= (_data->pid(i) == PROTON);
  _proton &= (abs(_dt->dt_P(i)) < 0.5 || abs(_dt->dt_ctof_P(i)) < 0.4);
  // _proton &= !(abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) < 0.2);
  _proton &= (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000);
  // min/max mom cuts
  // if (abs(_data->status(i)) < 4000) {
  //   _proton &= (_data->p(i) > 0.4);
  //   _proton &= (_data->p(i) < 5.0);
  // } else if (abs(_data->status(i)) >= 4000) {
  //   _proton &= (_data->p(i) > 0.0);
  //   _proton &= (_data->p(i) < 3.0);
  // }
  //   _proton &= (_data->p(i) > 0.2);
  // //_proton &= (abs(_data->chi2pid(i)) < 0.5);
  return _proton;
}
bool Cuts::IsPim(int i) {
  if (_data->gpart() <= i) return false;
  bool _pim = true;
  //   _pim &= (_data->charge(i) == NEGATIVE);
  _pim &= (_data->pid(i) == PIM);
  _pim &= (abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) < 0.5);
  _pim &= (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000);
  // min/max mom cuts
  // if (abs(_data->status(i)) < 4000){
  //   _pim &= (_data->p(i) > 0.4);
  //   _pim &= (_data->p(i) < 5.0);
  // } else if (abs(_data->status(i)) >= 4000) {
  //   _pim &= (_data->p(i) > 0.0);
  //   _pim &= (_data->p(i) < 3.0);
  // }

  // //_pim &= (abs(_data->chi2pid(i)) < 0.5);

  //_pim &= DC_fiducial_cut_theta_phi(i);
  //_pim &= Hadron_Delta_vz_cut(i);
  //_pim &= Hadron_Chi2pid_cut(i);

  return _pim;
}

// // bool Cuts::IsmissingPim(int i) {
// //   if (_data->gpart() <= i) return false;
// //   bool _missingPim = true;
// //   _missingPim &= (_data->charge(i) == NEGATIVE);
// //   //_missingPim &= (abs(_dt->dt_Pi(i)) < 0.5 || abs(_dt->dt_ctof_Pi(i)) <
// //   0.5);
// //   //  _missingPim &= (_data->pid(i) != PIM);
// //   _missingPim &= (_data->p(i) > 0.2);
// //
// //   //_pim &= (abs(_data->chi2pid(i)) < 0.5);
// //   return _missingPim;
// // }
// /////////////////////// uconn_Cuts ///////////////////////
bool uconn_Cuts::ElectronCuts() {
  bool cut = true;
  cut &= (_data->gpart() > 0);
  if (!cut) return false;

  cut &= (_data->gpart() < 20);
  //
  cut &= (_data->charge(0) == NEGATIVE);
  cut &= (_data->pid(0) == ELECTRON);
  cut &= (_data->p(0) > 1.50);
  cut &= (2000 <= abs(_data->status(0)) && abs(_data->status(0)) < 4000);
  cut &= (abs(_data->chi2pid(0)) < 3);  ////////////// check it.......
  cut &= CC_nphe_cut();
  cut &= EC_outer_vs_EC_inner_cut();
  cut &= EC_sampling_fraction_cut();
  cut &= EC_hit_position_fiducial_cut_homogeneous();
  cut &= DC_fiducial_cut_XY();
  cut &= DC_z_vertex_cut();
  return cut;
}
bool uconn_Cuts::HadronsCuts(int i) {
  bool cut = true;
  if (2000 <= abs(_data->status(i)) && abs(_data->status(i)) < 4000) cut &= DC_fiducial_cut_theta_phi(i);
  cut &= Hadron_Delta_vz_cut(i);
  cut &= Hadron_Chi2pid_cut(i);
  return cut;
}
//
//
// // bool uconn_Cuts::CC_nphe_cut(double nphe) {
// //         double nphe_min = 2;
// //         return nphe > nphe_min;
// // }
//
//
bool uconn_Cuts::CC_nphe_cut() {
  float nphe_min = 2;
  return (_data->cc_nphe_tot(0) > nphe_min);
}
//
// bool uconn_Cuts::EC_outer_vs_EC_inner_cut(double pcal_energy) {
//         double edep_tight = 0.06, edep_medium = 0.07, edep_loose = 0.09;
//         return pcal_energy > edep_medium;
// }
bool uconn_Cuts::EC_outer_vs_EC_inner_cut() {
  double edep_tight = 0.06, edep_medium = 0.07, edep_loose = 0.09;
  return (_data->ec_pcal_energy(0) > edep_medium);
}

bool uconn_Cuts::EC_sampling_fraction_cut() {
  double ecal_e_sampl_mu[3][6] = {{0.2531, 0.2550, 0.2514, 0.2494, 0.2528, 0.2521},
                                  {-0.6502, -0.7472, -0.7674, -0.4913, -0.3988, -0.703},
                                  {4.939, 5.350, 5.102, 6.440, 6.149, 4.957}};
  double ecal_e_sampl_sigm[3][6] = {{2.726e-3, 4.157e-3, 5.222e-3, 5.398e-3, 8.453e-3, 6.533e-3},
                                    {1.062, 0.859, 0.5564, 0.6576, 0.3242, 0.4423},
                                    {-4.089, -3.318, -2.078, -2.565, -0.8223, -1.274}};
  double sigma_range = 3.5;
  // double ectotal_energy = pcal_energy + ecin_energy + ecout_energy;
  int isec = (_data->ec_pcal_sec(0) - 1);
  double mean =
      ecal_e_sampl_mu[0][isec] + ecal_e_sampl_mu[1][isec] / 1000 * pow(_data->p(0) - ecal_e_sampl_mu[2][isec], 2);
  double sigma =
      ecal_e_sampl_sigm[0][isec] + ecal_e_sampl_sigm[1][isec] / (10 * (_data->p(0) - ecal_e_sampl_sigm[2][isec]));
  double upper_lim_total = mean + sigma_range * sigma;
  double lower_lim_total = mean - sigma_range * sigma;
  bool pass_band = _data->ec_tot_energy(0) / _data->p(0) <= upper_lim_total &&
                   _data->ec_tot_energy(0) / _data->p(0) >= lower_lim_total;
  bool pass_triangle = false;
  if (_data->p(0) < 4.5) {
    pass_triangle = true;
  } else {
    pass_triangle = _data->ec_ecin_energy(0) / _data->p(0) > (0.2 - _data->ec_pcal_energy(0)) / _data->p(0);
  }
  return (pass_band && pass_triangle);
}

bool uconn_Cuts::EC_hit_position_fiducial_cut_homogeneous() {
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
  //////////////////////////////////////////////////////////////
  int isec = (_data->ec_pcal_sec(0) - 1);
  double min_v = min_v_loose_inb[isec];
  double max_v = max_v_loose_inb[isec];
  double min_w = min_w_loose_inb[isec];
  double max_w = max_w_loose_inb[isec];
  return (_data->ec_pcal_lv(0) > min_v && _data->ec_pcal_lv(0) < max_v && _data->ec_pcal_lw(0) > min_w &&
          _data->ec_pcal_lw(0) < max_w);
}

bool uconn_Cuts::DC_fiducial_cut_XY() {
  // bool _dc_fid_cut = true;
  // bool isinbending = true;
  // new cut parameters for the linear cut based on x and y coordinates (inbending field):
  // replace it in the function: bool DC_fiducial_cut_XY(int j, int region)
  // (optimized for electrons, do not use it for hadrons)
  //

  double maxparams_in[6][6][3][2] = {{{{-14.563, 0.60032}, {-19.6768, 0.58729}, {-22.2531, 0.544896}},
                                      {{-12.7486, 0.587631}, {-18.8093, 0.571584}, {-19.077, 0.519895}},
                                      {{-11.3481, 0.536385}, {-18.8912, 0.58099}, {-18.8584, 0.515956}},
                                      {{-10.7248, 0.52678}, {-18.2058, 0.559429}, {-22.0058, 0.53808}},
                                      {{-16.9644, 0.688637}, {-17.1012, 0.543961}, {-21.3974, 0.495489}},
                                      {{-13.4454, 0.594051}, {-19.4173, 0.58875}, {-22.8771, 0.558029}}},
                                     {{{-6.2928, 0.541828}, {-16.7759, 0.57962}, {-32.5232, 0.599023}},
                                      {{-6.3996, 0.543619}, {-16.7429, 0.578472}, {-32.5408, 0.600826}},
                                      {{-5.49712, 0.53463}, {-16.1294, 0.576928}, {-32.5171, 0.597735}},
                                      {{-6.4374, 0.54839}, {-16.9511, 0.582143}, {-33.0501, 0.59995}},
                                      {{-5.30128, 0.529377}, {-16.1229, 0.579019}, {-30.7768, 0.593861}},
                                      {{-5.89201, 0.541124}, {-16.1245, 0.575001}, {-32.2617, 0.601506}}},
                                     {{{-6.3618, 0.546384}, {-17.0277, 0.582344}, {-34.9276, 0.612875}},
                                      {{-6.36432, 0.546268}, {-15.8404, 0.574102}, {-33.0627, 0.599142}},
                                      {{-6.34357, 0.548411}, {-16.0496, 0.575913}, {-34.8535, 0.610211}},
                                      {{-5.8568, 0.541784}, {-16.1124, 0.576473}, {-32.8547, 0.599033}},
                                      {{-5.91941, 0.536801}, {-15.726, 0.575211}, {-34.0964, 0.606777}},
                                      {{-5.55498, 0.536609}, {-15.9853, 0.579705}, {-33.4886, 0.606439}}},
                                     {{{-12.594, 0.613062}, {-18.4504, 0.588136}, {-16.3157, 0.529461}},
                                      {{-12.3417, 0.61231}, {-18.1498, 0.590748}, {-13.8106, 0.52335}},
                                      {{-12.1761, 0.609307}, {-15.919, 0.572156}, {-13.0598, 0.5194}},
                                      {{-12.5467, 0.612645}, {-16.2129, 0.572974}, {-12.8611, 0.51252}},
                                      {{-13.0976, 0.615928}, {-16.9233, 0.580972}, {-13.0906, 0.519738}},
                                      {{-12.884, 0.622133}, {-17.2566, 0.585572}, {-12.1874, 0.510124}}},
                                     {{{-6.51157, 0.545763}, {-16.4246, 0.583603}, {-32.2001, 0.60425}},
                                      {{-6.21169, 0.541872}, {-16.8484, 0.591172}, {-31.7785, 0.606234}},
                                      {{-5.89452, 0.54464}, {-16.612, 0.591506}, {-29.9143, 0.589656}},
                                      {{-6.68908, 0.553374}, {-16.2993, 0.585165}, {-30.252, 0.59519}},
                                      {{-6.17185, 0.540496}, {-16.7197, 0.591664}, {-31.619, 0.608306}},
                                      {{-5.7526, 0.541761}, {-16.2054, 0.587326}, {-31.3653, 0.604081}}},
                                     {{{-11.8798, 0.62389}, {-20.2212, 0.610786}, {-16.4137, 0.51337}},
                                      {{-12.0817, 0.631621}, {-20.7511, 0.610844}, {-16.9407, 0.522958}},
                                      {{-9.72746, 0.605471}, {-20.4903, 0.622337}, {-15.3363, 0.520589}},
                                      {{-12.4566, 0.627481}, {-20.238, 0.606098}, {-20.7651, 0.56974}},
                                      {{-11.6712, 0.622265}, {-18.2649, 0.591062}, {-19.2569, 0.580894}},
                                      {{-12.0943, 0.630674}, {-22.4432, 0.633366}, {-17.2197, 0.537965}}}};
  double minparams_in[6][6][3][2] = {{{{12.2692, -0.583057}, {17.6233, -0.605722}, {19.7018, -0.518429}},
                                      {{12.1191, -0.582662}, {16.8692, -0.56719}, {20.9153, -0.534871}},
                                      {{11.4562, -0.53549}, {19.3201, -0.590815}, {20.1025, -0.511234}},
                                      {{13.202, -0.563346}, {20.3542, -0.575843}, {23.6495, -0.54525}},
                                      {{12.0907, -0.547413}, {17.1319, -0.537551}, {17.861, -0.493782}},
                                      {{13.2856, -0.594915}, {18.5707, -0.597428}, {21.6804, -0.552287}}},
                                     {{{5.35616, -0.531295}, {16.9702, -0.583819}, {36.3388, -0.612192}},
                                      {{6.41665, -0.543249}, {17.3455, -0.584322}, {37.1294, -0.61791}},
                                      {{6.86336, -0.550492}, {17.2747, -0.575263}, {39.6389, -0.625934}},
                                      {{6.82938, -0.558897}, {17.8618, -0.599931}, {39.3376, -0.631517}},
                                      {{6.05547, -0.54347}, {15.7765, -0.569165}, {35.6589, -0.611349}},
                                      {{6.3468, -0.544882}, {16.7144, -0.578363}, {38.2501, -0.617055}}},
                                     {{{6.70668, -0.558853}, {17.0627, -0.587751}, {36.1194, -0.617417}},
                                      {{6.3848, -0.542992}, {16.6355, -0.581708}, {34.6781, -0.609794}},
                                      {{6.36802, -0.539521}, {15.9829, -0.569165}, {32.5691, -0.59588}},
                                      {{5.94912, -0.546191}, {18.0321, -0.601764}, {36.5238, -0.619185}},
                                      {{5.65108, -0.541684}, {15.5009, -0.567131}, {34.0489, -0.602048}},
                                      {{6.71064, -0.547956}, {16.4449, -0.577051}, {34.4375, -0.602515}}},
                                     {{{12.4734, -0.608063}, {16.1064, -0.575034}, {16.0751, -0.536452}},
                                      {{12.1936, -0.6034}, {15.9302, -0.571271}, {14.2791, -0.520157}},
                                      {{12.216, -0.600017}, {14.8741, -0.56304}, {11.1766, -0.498955}},
                                      {{12.7941, -0.616044}, {17.1516, -0.583616}, {11.6077, -0.500028}},
                                      {{12.7448, -0.611315}, {16.2814, -0.572461}, {13.1033, -0.506663}},
                                      {{12.7949, -0.612051}, {16.1565, -0.569143}, {12.9295, -0.504203}}},
                                     {{{7.19022, -0.562083}, {16.5946, -0.591266}, {31.9033, -0.589167}},
                                      {{7.80002, -0.571429}, {17.8587, -0.595543}, {36.5772, -0.630136}},
                                      {{7.96121, -0.569485}, {17.8085, -0.592936}, {37.553, -0.632848}},
                                      {{7.52041, -0.566112}, {17.3385, -0.603462}, {33.7712, -0.606047}},
                                      {{7.35796, -0.562782}, {15.2865, -0.57433}, {29.8283, -0.574685}},
                                      {{7.80003, -0.571429}, {16.1751, -0.583286}, {39.1972, -0.642803}}},
                                     {{{13.4466, -0.633911}, {22.0097, -0.62205}, {18.8862, -0.519652}},
                                      {{13.0534, -0.626648}, {20.2994, -0.60581}, {19.3973, -0.573994}},
                                      {{12.547, -0.62145}, {18.9322, -0.596491}, {16.2331, -0.546036}},
                                      {{14.5339, -0.64585}, {20.0211, -0.608462}, {19.0405, -0.563914}},
                                      {{12.7388, -0.617954}, {21.1677, -0.621012}, {15.4502, -0.525165}},
                                      {{13.4019, -0.63075}, {16.6584, -0.554797}, {19.0302, -0.55004}}}};
  // // double maxparams_out[6][6][3][2] =
  // double maxparams_out [6][6][3][2] = {{{{-9.86221, 0.565985}, {-16.4397, 0.569087}, {-29.7787, 0.586842}},
  //                                       {{-10.2065, 0.565541}, {-16.5554, 0.571394}, {-28.933, 0.582078}},
  //                                       {{-8.48034, 0.550706}, {-16.4397, 0.569087}, {-27.1037, 0.563767}},
  //                                       {{-6.77188, 0.53062}, {-16.4397, 0.569087}, {-30.485, 0.587534}},
  //                                       {{-8.00705, 0.543502}, {-16.4038, 0.571178}, {-27.7934, 0.573472}},
  //                                       {{-10.3328, 0.571942}, {-16.69, 0.575252}, {-30.8177, 0.592418}}},
  //                                      {{{-5.43811, 0.550931}, {-17.1906, 0.57936}, {-18.552, 0.546789}},
  //                                       {{-5.46281, 0.549659}, {-18.0351, 0.588876}, {-17.6981, 0.549803}},
  //                                       {{-3.26087, 0.531677}, {-16.3762, 0.578005}, {-17.6831, 0.55049}},
  //                                       {{-4.5985, 0.542017}, {-17.2735, 0.581566}, {-16.7013, 0.538853}},
  //                                       {{-6.83053, 0.561019}, {-16.5082, 0.579816}, {-18.0846, 0.553592}},
  //                                       {{-5.67358, 0.5558}, {-18.8196, 0.594965}, {-19.4333, 0.560965}}},
  //                                      {{{-12.6317, 0.611023}, {-16.5644, 0.578978}, {-11.5882, 0.496324}},
  //                                       {{-12.8886, 0.614807}, {-17.0847, 0.584072}, {-14.9561, 0.532125}},
  //                                       {{-11.4504, 0.600574}, {-16.3862, 0.57885}, {-12.3309, 0.515431}},
  //                                       {{-12.2256, 0.609801}, {-16.2134, 0.574306}, {-12.7661, 0.515787}},
  //                                       {{-12.6311, 0.611069}, {-16.2486, 0.577577}, {-12.6783, 0.519597}},
  //                                       {{-12.6937, 0.615423}, {-16.1427, 0.57847}, {-11.5156, 0.509458}}},
  //                                      {{{-5.95834, 0.538479}, {-15.8909, 0.570164}, {-30.2922, 0.586335}},
  //                                       {{-6.15277, 0.542134}, {-16.1129, 0.573794}, {-31.6024, 0.592681}},
  //                                       {{-6.12341, 0.542023}, {-16.1611, 0.575971}, {-29.8604, 0.581528}},
  //                                       {{-6.37691, 0.546536}, {-16.8501, 0.580239}, {-30.0623, 0.580497}},
  //                                       {{-5.96605, 0.537402}, {-15.7154, 0.5704}, {-31.2955, 0.594146}},
  //                                       {{-5.86704, 0.539556}, {-16.2268, 0.580945}, {-31.2345, 0.590849}}},
  //                                      {{{-11.7796, 0.614043}, {-19.0763, 0.595015}, {-18.804, 0.559538}},
  //                                       {{-12.4399, 0.623126}, {-19.1733, 0.600646}, {-17.675, 0.557016}},
  //                                       {{-10.4158, 0.605483}, {-18.0044, 0.595497}, {-17.5441, 0.556504}},
  //                                       {{-12.1552, 0.617782}, {-19.7134, 0.603519}, {-17.3756, 0.549676}},
  //                                       {{-11.3901, 0.612121}, {-18.2429, 0.596796}, {-10.0097, 0.482578}},
  //                                       {{-12.5004, 0.626384}, {-19.9266, 0.60993}, {-16.4668, 0.543148}}},
  //                                      {{{-5.60572, 0.537153}, {-16.3196, 0.582537}, {-32.4336, 0.601487}},
  //                                       {{-5.52369, 0.532985}, {-15.2055, 0.568935}, {-31.9046, 0.600079}},
  //                                       {{-5.78558, 0.546316}, {-16.3328, 0.583765}, {-36.0074, 0.617008}},
  //                                       {{-5.82321, 0.542839}, {-15.9551, 0.580441}, {-31.4304, 0.597132}},
  //                                       {{-5.36526, 0.535923}, {-15.9219, 0.586886}, {-30.4245, 0.599613}},
  //                                       {{-5.14766, 0.53037}, {-14.1986, 0.561504}, {-31.7548, 0.60233}}}};
  // double minparams_out[6][6][3][2] ={{{{8.07831, -0.548881}, {16.4382, -0.569075}, {33.7768, -0.607402}},
  //                                     {{8.51057, -0.551773}, {16.7782, -0.571381}, {32.2613, -0.600686}},
  //                                     {{8.5232, -0.552628}, {16.4274, -0.56775}, {31.1516, -0.584708}},
  //                                     {{7.98845, -0.544571}, {16.4381, -0.569077}, {31.8093, -0.595237}},
  //                                     {{7.46705, -0.538557}, {16.7414, -0.573345}, {31.1888, -0.586751}},
  //                                     {{7.82627, -0.538957}, {16.2409, -0.565872}, {32.1089, -0.596846}}},
  //                                    {{{7.1519, -0.563678}, {16.1038, -0.571795}, {20.0449, -0.559802}},
  //                                     {{6.38228, -0.553174}, {16.4526, -0.576382}, {19.3523, -0.556484}},
  //                                     {{7.11359, -0.561586}, {17.2815, -0.578095}, {14.9667, -0.53314}},
  //                                     {{5.89053, -0.556406}, {17.4946, -0.585038}, {17.3607, -0.545739}},
  //                                     {{7.08253, -0.562099}, {15.1516, -0.569192}, {16.9665, -0.545949}},
  //                                     {{5.53089, -0.546315}, {16.4962, -0.574014}, {17.9593, -0.545788}}},
  //                                    {{{12.4879, -0.610527}, {16.7782, -0.575065}, {11.7704, -0.511182}},
  //                                     {{12.1931, -0.604779}, {15.6443, -0.560967}, {12.7304, -0.515606}},
  //                                     {{12.206, -0.602999}, {16.5979, -0.573274}, {12.3971, -0.513795}},
  //                                     {{11.5538, -0.604186}, {16.6974, -0.576753}, {12.7385, -0.517811}},
  //                                     {{12.9718, -0.611968}, {17.7233, -0.583943}, {10.6601, -0.49233}},
  //                                     {{12.2966, -0.607592}, {15.923, -0.564133}, {13.9314, -0.525363}}},
  //                                    {{{5.92493, -0.539308}, {17.4444, -0.586183}, {31.6974, -0.591988}},
  //                                     {{5.467, -0.525876}, {16.0649, -0.570869}, {30.5937, -0.590071}},
  //                                     {{5.67798, -0.531096}, {16.5072, -0.57205}, {30.7922, -0.586727}},
  //                                     {{6.85795, -0.558336}, {14.9425, -0.545596}, {31.3159, -0.592865}},
  //                                     {{6.0155, -0.545283}, {16.0649, -0.570869}, {30.6644, -0.587002}},
  //                                     {{6.18343, -0.539055}, {17.4516, -0.583221}, {32.6264, -0.594317}}},
  //                                    {{{12.9118, -0.618907}, {19.7061, -0.60171}, {18.9352, -0.559461}},
  //                                     {{13.0612, -0.618743}, {19.0954, -0.595406}, {19.7019, -0.568119}},
  //                                     {{12.4007, -0.613459}, {17.544, -0.581147}, {12.8175, -0.511017}},
  //                                     {{13.3144, -0.625596}, {18.9225, -0.594001}, {15.1524, -0.530046}},
  //                                     {{13.101, -0.620887}, {18.5616, -0.595279}, {14.8807, -0.533111}},
  //                                     {{12.2964, -0.613529}, {19.0686, -0.595276}, {19.2596, -0.562706}}},
  //                                    {{{5.34118, -0.530584}, {16.3015, -0.585185}, {38.7808, -0.641362}},
  //                                     {{6.68051, -0.548747}, {16.4236, -0.583598}, {38.4718, -0.630423}},
  //                                     {{6.87, -0.552602}, {16.4285, -0.57977}, {36.8889, -0.624053}},
  //                                     {{7.15338, -0.565067}, {16.9387, -0.595922}, {37.2398, -0.624177}},
  //                                     {{6.06995, -0.550001}, {15.7376, -0.577755}, {32.6004, -0.601595}},
  //                                     {{6.20459, -0.543148}, {14.6326, -0.561623}, {39.2154, -0.631762}}}};
  // double minparams[6][6][3][2]=  {minparams_in[6][6][3][2]};
  // double maxparams[6][6][3][2]=  {maxparams_in[6][6][3][2]};
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

  double calc_min1 = minparams_in[pid][dc_sector][region_1 - 1][0] + minparams_in[pid][dc_sector][region_1 - 1][1] * X1;
  double calc_max1 = maxparams_in[pid][dc_sector][region_1 - 1][0] + maxparams_in[pid][dc_sector][region_1 - 1][1] * X1;
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

  double calc_min2 = minparams_in[pid][dc_sector][region_2 - 1][0] + minparams_in[pid][dc_sector][region_2 - 1][1] * X2;
  double calc_max2 = maxparams_in[pid][dc_sector][region_2 - 1][0] + maxparams_in[pid][dc_sector][region_2 - 1][1] * X2;
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

  double calc_min3 = minparams_in[pid][dc_sector][region_3 - 1][0] + minparams_in[pid][dc_sector][region_3 - 1][1] * X3;
  double calc_max3 = maxparams_in[pid][dc_sector][region_3 - 1][0] + maxparams_in[pid][dc_sector][region_3 - 1][1] * X3;
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
  // double calc_min = minparams[pid][dc_sector - 1][region - 1][0] + minparams[pid][dc_sector - 1][region - 1][1] * X;
  // double calc_max = maxparams[pid][dc_sector - 1][region - 1][0] + maxparams[pid][dc_sector - 1][region - 1][1] * X;
  // return (Y > calc_min) && (Y < calc_max);
}

bool uconn_Cuts::DC_z_vertex_cut() {
  int pcal_sector = _data->ec_pcal_sec(0);
  float partvz = _data->vz(0);
  bool isinbending = false;

  float vz_min_sect_inb[6] = {-13, -13, -13, -13, -13, -13};
  float vz_max_sect_inb[6] = {12, 12, 12, 12, 12, 12};

  float vz_min_sect_outb[6] = {-18, -18, -18, -18, -18, -18};
  float vz_max_sect_outb[6] = {10, 10, 10, 10, 10, 10};

  float vz_min_sect[6];
  float vz_max_sect[6];

  for (int i = 0; i < 6; i++) {
    if (isinbending) {
      vz_min_sect[i] = vz_min_sect_inb[i];
      vz_max_sect[i] = vz_max_sect_inb[i];
    } else {
      vz_min_sect[i] = vz_min_sect_outb[i];
      vz_max_sect[i] = vz_max_sect_outb[i];
    }
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
//  * @param region specify fiducial uconn_Cuts for which region to use
//  * @param trajx x for region 1 or 2 or 3 from REC::Traj
//  * @param trajy y for region 1 or 2 or 3 from REC::Traj
//  * @param trajz z for region 1 or 2 or 3 from REC::Traj
//  * @param partpid pid assigned to particle candidate
//  * @param isinbending True if magnetic field is inbending
//  */
bool uconn_Cuts::DC_fiducial_cut_theta_phi(int i) {
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

  // //fitted values outbending
  // double maxparams_out[6][6][3][4] = {
  //         {   {{-3.69457, 12.3755, -0.41328, 0.00129631},{-54.3237, 40.3308, -2.39952, 0.0181339},{-39.8661, 27.1428,
  //         -0.907303, 0.00220974}},
  //             {{-37.6199, 26.2865, -0.826366, 0.000862203},{-72.4212, 54.7953, -4.04856,
  //             0.0373308},{-21.1791, 17.0759, -0.391795, 0.00151085}},
  //             {{-0.421685, 10.482, -0.272111, 8.69408e-05},{-43.3635, 32.746, -1.6541, 0.0101454},{-62.6387, 41.1869,
  //             -1.97298, 0.0107022}},
  //             {{-42.0766, 29.6387, -0.993426, 1.97101e-09},{-44.7036, 33.0587, -1.64131,
  //             0.0099416},{-47.2703, 32.6109, -1.46533, 0.00817871}},
  //             {{-22.2035, 20.6894, -0.689051, 0.000592423},{-74.6572, 54.7065, -3.83999,
  //             0.0351952},{-38.9183, 25.7212, -0.711499, 2.5796e-12}},
  //             {{-52.078, 45.571, -3.71942, 0.0376577},{-65.4047, 49.1723, -3.36623, 0.0288435},{-53.9611, 35.9294,
  //             -1.58589, 0.00772417}}},
  //         {   {{-2.20312e-07, 13.0916, -0.864184, 0.0086342},{-6.44026e-08, 12.056, -0.675801,
  //         0.00643464},{-20.2596, 23.5977, -1.545, 0.0141047}},
  //             {{-4.42537e-05, 10.2799, -0.322454, 0.00154825},{-1.63659e-07, 11.0228, -0.451412,
  //             0.00308633},{-8.5382, 15.6903, -0.785315, 0.00602734}},
  //             {{-2.32088, 11.6343, -0.363509, 0.000902217},{-0.301128, 12.0319, -0.643794,
  //             0.00581994},{-22.4378, 25.2772, -1.73656, 0.0164181}},
  //             {{-7.40627, 13.601, -0.382439, 2.45262e-05},{-5.50415e-08, 11.9792, -0.652368,
  //             0.00597647},{-15.1608, 20.6455, -1.33827, 0.0127123}},
  //             {{-0.203913, 10.7032, -0.322123, 0.000691162},{-1.73184e-07, 10.735, -0.379993,
  //             0.00196037},{-0.155443, 10.1794, -0.249841, 6.24278e-05}},
  //             {{-1.87352e-07, 12.4226, -0.730141, 0.0068049},{-1.40236e-07, 12.5356, -0.750615,
  //             0.00719921},{-16.8681, 21.8555, -1.43078, 0.0131935}}},
  //         {   {{-8.89326e-08, 10.0681, -0.240869, 9.9612e-12},{-15.2705, 21.635, -1.55291,
  //         0.0166645},{-10.5976, 17.9928, -1.08432, 0.00950807}},
  //             {{-0.00389562, 10.2092, -0.254082, 4.15737e-06},{-9.16032e-11, 10.527, -0.334641,
  //             0.00129061},{-9.63013e-07, 11.0668, -0.42453, 0.0022955}},
  //             {{-2.40163e-06, 13.4151, -0.949883, 0.0107662},{-1.60937e-07, 10.5128, -0.35046,
  //             0.00173787},{-29.2647, 30.1252, -2.20552, 0.0213809}},
  //             {{-2.69733e-08, 11.7703, -0.589854, 0.00482124},{-3.77564e-08, 11.3764, -0.527037,
  //             0.00416671},{-4.85047, 13.7737, -0.650441, 0.0047428}},
  //             {{-3.90816e-07, 12.2683, -0.692591, 0.00625884},{-9.70203e-10, 11.0335, -0.438323,
  //             0.00275342},{-2.54193, 13.5404, -0.76861, 0.00684486}},
  //             {{-3.23439e-10, 10.7412, -0.348557, 0.00113794},{-1.79623, 11.7499, -0.449432,
  //             0.00247294},{-13.1393, 19.4689, -1.17148, 0.00984086}}},
  //         {   {{-5.07611e-08, 11.7796, -0.516966, 0.00295389},{-4.87018, 12.2727,
  //         -0.322719, 9.12315e-06},{-35.9369, 31.015, -1.95133, 0.0169834}},
  //             {{-1.32385e-07, 11.6454, -0.495467, 0.00272602},{-2.70664, 12.0151, -0.434014,
  //             0.00203292},{-8.97137, 15.0453, -0.646138, 0.00429196}},
  //             {{-7.92247e-09, 12.5189, -0.682231, 0.00539531},{-0.0942499, 10.3465, -0.280521,
  //             0.000405358},{-19.7485, 21.7919, -1.24334, 0.0105088}},
  //             {{-8.50093e-11, 10.739, -0.302295, 5.6862e-11},{-0.184771, 10.4358, -0.285869,
  //             0.000389546},{-21.9469, 24.9675, -1.77893, 0.0183075}},
  //             {{-4.34589, 12.5902, -0.362849, 4.996e-15},{-0.000684493, 10.6055, -0.332363,
  //             0.00104632},{-21.328, 22.0864, -1.20993, 0.00989151}},
  //             {{-0.0202168, 12.0097, -0.539165, 0.00299034},{-0.5239, 10.7167, -0.309141,
  //             0.000535617},{-10.0299, 16.3179, -0.812315, 0.00617078}}},
  //         {   {{-0.169908, 10.902, -0.353938, 0.00100715},{-3.2818, 13.2193, -0.65495,
  //         0.00515117},{-0.013532, 8.51331, -0.070239, 1.755e-05}},
  //             {{-8.51985e-08, 11.6512, -0.56808, 0.00453582},{-1.2381e-07, 10.6653, -0.368149,
  //             0.00181989},{-9.30287e-08, 10.0352, -0.254321, 0.000417053}},
  //             {{-0.150407, 10.6338, -0.308676, 0.000481694},{-0.00186321, 10.4259, -0.303092,
  //             0.00073092},{-21.3328, 28.0803, -2.37912, 0.025101}},
  //             {{-14.4411, 19.817, -1.13705, 0.00894685},{-6.25263e-09, 11.7414, -0.586098,
  //             0.00478932},{-5.49193, 16.1248, -1.11306, 0.0115644}},
  //             {{-1.54761, 12.0015, -0.462506, 0.00204729},{-5.72883, 14.9638, -0.795325,
  //             0.00616222},{-50.229, 45.8456, -3.88803, 0.0414729}},
  //             {{-40.7531, 33.6269, -2.03771, 0.01609},{-1.33363e-09, 11.9894, -0.614358,
  //             0.004924},{-27.2506, 29.2602, -2.1426, 0.0203235}}},
  //         {   {{-1.62999e-10, 14.0422, -1.03609, 0.0107179},{-6.71565, 15.6964, -0.887791,
  //         0.00740777},{-38.9148, 32.9935, -2.09023, 0.0177295}},
  //             {{-1.09078e-05, 13.4131, -0.878092, 0.00825152},{-15.0102, 21.6968, -1.4935,
  //             0.0138851},{-19.5261, 20.3932, -0.969464, 0.00661531}},
  //             {{-1.39619e-08, 12.3593, -0.618488, 0.00415536},{-5.38271e-07, 11.5631, -0.512607,
  //             0.00334452},{-23.0902, 24.7093, -1.57315, 0.0140132}},
  //             {{-1.73908e-08, 12.0348, -0.591608, 0.00423834},{-8.35134, 17.3066, -1.11555,
  //             0.010407},{-2.74909e-07, 9.59202, -0.216455, 0.000527479}},
  //             {{-0.0449157, 10.5243, -0.334389, 0.00134555},{-0.0143489, 10.0993,
  //             -0.2434, 1.57595e-10},{-22.3661, 23.2499, -1.32946, 0.0108047}},
  //             {{-5.83731e-07, 14.5234, -1.14022, 0.0122177},{-1.4586e-08, 11.6946, -0.520935,
  //             0.00324975},{-12.4252, 16.3216, -0.652566, 0.00365791}}}
  // };
  //
  // double minparams_out[6][6][3][4] = {
  //         {   {{3.73672, -12.3584,0.390616, -0.000795415},{51.644, -37.8546,1.99228, -0.0119973},{32.3551,
  //         -22.9742,0.624096, -4.30811e-05}},
  //             {{6.11614, -13.6358,0.491668, -0.0018637},{47.5098, -35.902,1.97535, -0.0134876},{82.9536,
  //             -58.2741,4.12662, -0.0378612}},
  //             {{0.000950108, -7.99619,0.000506416, -0.0020788},{64.0688, -47.8642,3.16007, -0.025878},{70.0064,
  //             -50.3249,3.38975, -0.029639}},
  //             {{37.0145, -35.0316,2.61892, -0.0250306},{14.5954, -15.6554,0.426733, -0.000879865},{28.9035,
  //             -21.5279,0.610475, -0.00087271}},
  //             {{5.65685, -13.3347,0.400781, -1.46612e-11},{67.3504, -50.152,3.33677, -0.0270726},{47.0772,
  //             -32.1506,1.38851, -0.00719898}},
  //             {{8.95987, -15.1646,0.585477, -0.00246174},{41.6154, -29.7967,1.1817, -0.00403765},{61.1631,
  //             -41.6465,2.32522, -0.0175271}}},
  //         {   {{8.80954e-10, -11.0364,0.413853, -0.00210254},{6.50072e-08, -11.2505,0.501571, -0.00380973},{10.9643,
  //         -17.4701,0.989297, -0.00860789}},
  //             {{2.33292e-08, -11.2353,0.470728, -0.00309666},{2.29373e-07, -11.2458,0.50218, -0.00383969},{29.5429,
  //             -29.9965,2.19166, -0.021366}},
  //             {{1.61826e-08, -11.861,0.577321, -0.00433276},{2.9436e-07, -11.5738,0.581015, -0.00503307},{19.5142,
  //             -23.451,1.58724, -0.0151339}},
  //             {{2.07231e-09, -12.7453,0.751184, -0.00664181},{1.77802e-07, -11.4574,0.537367, -0.00422656},{12.5683,
  //             -18.4632,1.05475, -0.00892182}},
  //             {{7.6216e-08, -13.9769,1.01051, -0.0107372},{1.33092e-08, -11.9128,0.628521, -0.00550105},{13.5537,
  //             -20.1708,1.32578, -0.0123213}},
  //             {{9.25941, -19.658,1.51566, -0.0157124},{6.25983e-10, -11.6806,0.599263, -0.00532588},{17.0479,
  //             -22.0046,1.47474, -0.0140475}}},
  //         {   {{4.65436e-08, -11.1925,0.466196, -0.00308992},{18.4968, -22.5122,1.4594, -0.0135962},{18.9488,
  //         -23.3348,1.57414, -0.0146183}},
  //             {{3.67722e-08, -10.9985,0.428395, -0.00257574},{16.3745, -21.0105,1.3093, -0.0119156},{11.4404,
  //             -18.6679,1.15919, -0.010306}},
  //             {{1.46846e-08, -10.865,0.398638, -0.00212392},{20.7337, -23.3738,1.46852, -0.0130115},{28.2098,
  //             -28.9406,2.05908, -0.0197782}},
  //             {{0.237058, -10.4694,0.271985, -1.08731e-07},{2.32759, -11.9354,0.469887, -0.00291497},{13.287,
  //             -20.8621,1.49656, -0.0148999}},
  //             {{0.000149907, -10.4632,0.294713, -0.000431947},{6.96663, -15.3946,0.845078, -0.00724722},{11.0939,
  //             -17.4733,0.944239, -0.00747728}},
  //             {{3.10006e-08, -10.1416,0.247764, -1.36913e-11},{5.41915, -14.6085,0.795369, -0.00684375},{5.89127,
  //             -13.0881,0.453024, -0.0020325}}},
  //         {   {{4.16588e-09, -12.9305,0.749425, -0.00611725},{5.65263, -14.1661,0.637395, -0.00400239},{4.66325,
  //         -12.9519,0.565753, -0.00442033}},
  //             {{8.0428e-08, -13.1625,0.836744, -0.00778246},{12.3243, -18.8718,1.11103, -0.00917354},{7.20312,
  //             -16.0935,0.987223, -0.00930883}},
  //             {{0.00147165, -10.4992,0.280542, -1.79846e-06},{3.20232, -11.6892,0.350774, -0.00101099},{8.14117e-08,
  //             -10.9813,0.524839, -0.00507885}},
  //             {{0.470888, -13.5446,0.820782, -0.00768941},{3.9697, -13.0821,0.540847, -0.00303209},{3.44817,
  //             -12.3932,0.533804, -0.00414144}},
  //             {{1.05038e-08, -10.6539,0.297078, -6.04694e-05},{15.0983, -21.1791,1.38383, -0.0124058},{17.3666,
  //             -20.3986,1.16663, -0.0102393}},
  //             {{8.49365e-07, -13.765,0.964056, -0.00956575},{9.38084, -16.7385,0.904339, -0.00707907},{12.1048,
  //             -17.3704,0.91318, -0.00757461}}},
  //         {   {{10.6378, -19.5017,1.45275, -0.017057},{1.24368e-08, -10.5134,0.338985, -0.00143696},{37.3291,
  //         -35.1606,2.60092, -0.0242728}},
  //             {{19.1614, -24.0851,1.73932, -0.0185466},{14.1293, -19.8382,1.21613, -0.0107037},{20.9629,
  //             -24.0839,1.60283, -0.015173}},
  //             {{0.000450804, -8.15062,0.0103867, -2.00709e-05},{5.72496, -14.338,0.717819, -0.00567964},{16.9428,
  //             -21.8075,1.4216, -0.0131736}},
  //             {{6.15991e-10, -11.5278,0.536105, -0.00402223},{2.17842e-07, -10.5338,0.327427, -0.0010898},{20.7387,
  //             -24.3028,1.65004, -0.0155857}},
  //             {{0.650351, -10.6177,0.275393, -6.4664e-08},{8.05811, -16.1558,0.913735, -0.00788487},{0.308897,
  //             -10.2816,0.275186, -0.000561299}},
  //             {{0.427836, -10.168,0.240458, -5.90042e-06},{2.30661, -12.8686,0.664796, -0.00562626},{0.00499667,
  //             -11.6585,0.62597, -0.00619261}}},
  //         {   {{9.01249e-07, -11.8437,0.494125, -0.00223452},{14.3941, -21.2365,1.46048, -0.0137349},{13.7095,
  //         -15.4704,0.408961, -0.000312145}},
  //             {{0.000251044, -11.3084,0.438545, -0.0020791},{0.00847078, -12.6769,0.804431, -0.00836705},{1.09388,
  //             -9.66797,0.175278, -1.8721e-11}},
  //             {{4.04693e-10, -11.9001,0.585913, -0.00440376},{5.05178, -12.1514,0.31134, -0.000112735},{30.8105,
  //             -28.0795,1.73625, -0.0151639}},
  //             {{3.86607e-11, -13.471,0.889111, -0.0083617},{8.86591e-09, -9.25745,0.163052, -6.08491e-12},{27.1358,
  //             -24.3255,1.23326, -0.00891886}},
  //             {{0.196086, -11.7392,0.480055, -0.00224614},{0.18667, -10.5859,0.287231, -6.53153e-06},{14.8865,
  //             -17.1338,0.653576, -0.00333176}},
  //             {{2.7955e-07, -13.1311,0.848222, -0.00812719},{29.5508, -32.9514,2.77917, -0.0291596},{59.7514,
  //             -47.3033,3.54495, -0.0341802}}}
  // };
  //
  // // double[][][][] minparams = isinbending ? minparams_in : minparams_out;
  // // double[][][][] maxparams = isinbending ? maxparams_in : maxparams_out;

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
bool uconn_Cuts::Hadron_Delta_vz_cut(int i) {
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
bool uconn_Cuts::Hadron_Chi2pid_cut(int i) {
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

///////////////////// uconn_Cuts ///////////////////////
