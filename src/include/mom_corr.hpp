#ifndef REACTION_H_GUARD
#define REACTION_H_GUARD

#include <iostream>
#include "TLorentzRotation.h"
#include "TLorentzVector.h"
#include "branches.hpp"
#include "constants.hpp"
#include "physics.hpp"

class MomCorr {
 protected:
  std::shared_ptr<Branches12> _data;

  double _beam_energy = 10.6;

  std::unique_ptr<TLorentzVector> _mom_corr_pim;

  static const int Pim_mom_bins = 10;
  static const alpha_mom_corr = 0.5;

//   double pim_theta_corr[10] = {0.25, 0.25, 0.25, 0.3, -0.1, 1.05, 0.75, 0.2, 0.3, -1.5};
//   double pim_phi_corr[10] = {0.5, 0.1, 0.1, 0.3, -0.1, -0.3, -0.3, -0.1, -0.1, -0.1, 0.3};
//   double pim_mom_pars[10];  // mom bin 10
  float pim_mom_values[11] = {0, 0.36, 0.5, 0.65, 0.8, 0.95, 1.15, 1.45, 1.9, 2.5, 10};
  float pim_mom_ranges[10] = {0.36, 0.14, 0.15, 0.15, 0.15, 0.20, 0.30, 0.45, 0.6, 7.5};
  double pim_mom_corr[Pim_mom_bins] = {-0.004, -0.004, -0.004, -0.012, -0.025, -0.035, -0.035, -0.03, -0.03, -0.03};

  double _pim_mom_corr_fac = NAN;
  double _pim_mom = NAN;

  void MomCorr::corrfunc(int i) {

    for (size_t j = 0; j < Pim_mom_bins; j++) {
        double min = pim_mom_val[j];
        double max = pim_mom_val[j] + pim_mom_ranges[j];
        pim_mom = _pim->P();
        if (_pim_mom > min && pim_mom < max) {
          _mom_corr_pim->SetPxPyPzE(_data->px(i) * pim_mom_corr[j] * alpha_mom_corr, _data->py(i) * pim_mom_corr[j] * alpha_mom_corr,
                                     _data->pz(i) * pim_mom_corr[j] * alpha_mom_corr, _elec_mom * pim_mom_corr[j] * alpha_mom_corr);
      }
    }
}
    // void Histogram::populate_SF(const std::shared_ptr<Branches12>& _d, double min, double max, short index_sf) {
    //   if (_d->p(0) > min && _d->p(0) < max) {
    //     SF_1D[index_sf]->Fill(_d->ec_tot_energy(0) / _d->p(0));
    //   }
    // }
    // void Histogram::Fill_SF(const std::shared_ptr<Branches12>& _d) {
    //   short index_sf = 0;
    //   // populate_SF(_d, 0., 0.5, 0);
    //   // if (_d->p(0) > 1. && _d->p(0) < 1.5) SF_1D[2]->Fill(_d->ec_tot_energy(0) / _d->p(0));
    //   for (float x = 1; x <= Pim_mom_bins; x = x + pim) {
    //     float y = x + 0.3;
    //     populate_SF(_d, x, y, index_sf);
    //     index_sf++;
    //     j++;
    //   }
    // }

    int q2_bining(float q2) {
      if (q2 < 2.0)
        return 0;
      else if (q2 < 2.4)
        return 1;
      else if (q2 < 3.0)
        return 2;
      else if (q2 < 3.5)
        return 3;
      else if (q2 < 4.2)
        return 4;
      else if (q2 < 5.0)
        return 5;
      else if (q2 < 6.2)
        return 6;
      else if (q2 < 7.4)
        return 7;
      else if (q2 < 8.6)
        return 8;
      else if (q2 < 9.8)
        return 9;
      else
        return 10;
    }
  };
#endif