#ifndef CSV_DATA_H_GUARD
#define CSV_DATA_H_GUARD

#include <fstream>
#include <string>

struct csv_data {
  float pim_mom_mPim;
  float pim_theta_mPim;
  float pim_phi_mPim;
  float mm2_mPim;
  float mm2_mPim_corr;
  float weight_mPim;

  float pim_mom_corr;
  float pim_theta_corr;
  float pim_phi_corr;

  float gen_pim_mom;
  float gen_pim_theta;
  float gen_pim_phi;

  int status_Pim;
  int status_Pip;
  int status_Prot;
  // Static functions can be called without making a new struct
  static std::string header() {
    // Make a string for the header of the csv file mPim case
    // return "pim_mom_mPim,pim_theta_mPim,pim_phi_mPim,mm2_mPim,mm2_mPim_corr,weight";
    return "pim_mom_mPim,pim_theta_mPim,pim_phi_mPim,mm2_mPim,weight";
  }

  friend std ::ostream &operator<<(std::ostream &os, const csv_data &data) {
    // mPim
    os << std::setprecision(7);
    os << data.pim_mom_mPim << ",";
    os << data.pim_theta_mPim << ",";
    os << data.pim_phi_mPim << ",";
    os << std::setprecision(7);
    os << data.mm2_mPim << ",";
    // //  os << data.mm2_mPim_corr << ",";
    // os << data.status_Pim << ",";
    // os << data.status_Pip << ",";
    // os << data.status_Prot << ",";
    os << std::setprecision(1);
    os << data.weight_mPim << ",";

    return os;
  }
};

#endif
