// #include "mom_corr.hpp"

// MomCorr::MomCorr(const std::shared_ptr<Branches12>& data, float beam_energy) {
//   _data = data;
//   _beam = std::make_unique<TLorentzVector>();
//   _beam_energy = beam_energy;
//   _sector = data->dc_sec(0);

//   _mom_corr_elec = std::make_unique<TLorentzVector>();


//  for (size_t i = 0; i < Pim_mom_bins; i++)
//  {

//  }

//     }
//   }
// }
// MomCorr::~MomCorr() {}

//   _mom_corr_elec->SetPxPyPzE(_data->px(0) * _elec_mom_corrected, _data->py(0) * _elec_mom_corrected,
//                              _data->pz(0) * _elec_mom_corrected, _elec_mom * _elec_mom_corrected);