#ifndef CSV_DATA_H_GUARD
#define CSV_DATA_H_GUARD

#include <fstream>
#include <string>

struct csv_data {
  short electron_sector;
  float w;
  float q2;
  int status_pim;
  int status_pip;
  int status_prot;
  float scalar_product;
  float pim_mom_mPim;
  float pim_theta_mPim;
  float pim_phi_mPim;
  // float pim_mom_mPim_cm;
  float pim_theta_mPim_cm;
  float pim_phi_mPim_cm;
  float mm2_mPim;
  float weight_mPim;

  float pim_mom_exclusive;
  float pim_theta_exclusive;
  float pim_phi_exclusive;
  // float pim_mom_exclusive_cm;
  float pim_theta_exclusive_cm;
  float pim_phi_exclusive_cm;
  float mm2_exclusive;
  float energy_excl;
  float mm2_mPip;
  float mm2_mProt;

  float mm2_exclusive_at_zero;
  float weight_exclusive;

  // Static functions can be called without making a new struct
  static std::string header() {
    // Make a string for the header of the csv file
    // return "electron_sector,w,q2,pim_mom_mPim,pim_theta_mPim,pim_phi_mPim,mm2_mPim,weight,pim_theta_mPim_cm,pim_phi_"
    //        "mPim_cm";
    return "w,q2,stp,pim_mom_exclusive,pim_theta_exclusive,pim_phi_exclusive,mm2_exclusive,energy_excl,mm2_mPip,mm2_"
           "mProt,weight";
  }

  friend std ::ostream &operator<<(std::ostream &os, const csv_data &data) {
    os << std::setprecision(7);
    // os << data.electron_sector << ",";
    os << data.w << ",";
    os << data.q2 << ",";
    //     os << data.status_prot << ",";
    //     os << data.status_pip << ",";
    //     os << data.status_pim << ",";

    // os << data.pim_mom_mPim << ",";
    // os << data.pim_theta_mPim << ",";
    // os << data.pim_phi_mPim << ",";
    // os << data.mm2_mPim << ",";
    // os << data.weight_mPim << ",";
    // os << data.pim_theta_mPim_cm << ",";
    // os << data.pim_phi_mPim_cm << ",";

    os << data.scalar_product << ",";
    os << data.pim_mom_exclusive << ",";
    os << data.pim_theta_exclusive << ",";
    os << data.pim_phi_exclusive << ",";
    os << data.mm2_exclusive << ",";
    os << data.energy_excl << ",";
    os << data.mm2_mPip << ",";
    os << data.mm2_mProt << ",";
    // os << data.mm2_exclusive_at_zero<<",";
    os << data.weight_exclusive << ",";
    // os << data.pim_theta_exclusive_cm << ",";
    // os << data.pim_phi_exclusive_cm << ",";

    return os;
  }
};

#endif
