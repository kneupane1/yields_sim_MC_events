#ifndef CSV_DATA_H_GUARD
#define CSV_DATA_H_GUARD

#include <fstream>
#include <string>

struct csv_data {
  float w;
  float q2;
  float w_mc;
  float q2_mc;

  float pim_mom_mPim;
  float pim_theta_mPim;
  float pim_phi_mPim;

  float pim_mom_mPim_cm;
  float pim_theta_mPim_cm;
  float pim_phi_mPim_cm;

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

  float pip_mom_mPip;
  float pip_theta_mPip;
  float pip_phi_mPip;
  float mm2_mPip;
  float mm2_mPip_corr;
  float weight_mPip;

  float prot_mom_mProt;
  float prot_theta_mProt;
  float prot_phi_mProt;
  float mm2_mProt;
  float mm2_mProt_corr;
  float weight_mProt;

  // Static functions can be called without making a new struct
  static std::string header() {
    // Make a string for the header of the csv file mPim case
    return "w_gen,q2_gen,w_rec,q2_rec,weight";

    // return "pim_mom_mPim_cm,pim_theta_mPim_cm,pim_phi_mPim_cm,mm2_mPim,weight";
    // return "pim_mom_mPim,pim_theta_mPim,pim_phi_mPim,mm2_mPim,weight";
    // return "pip_mom_mPip,pip_theta_mPip,pip_phi_mPip,mm2_mPip,weight";
    // return "prot_mom_mProt,prot_theta_mProt,prot_phi_mProt,mm2_mProt,weight";
  }

  friend std ::ostream &operator<<(std::ostream &os, const csv_data &data) {

    os << data.w_mc << ",";
    os << data.q2_mc << ",";
    os << data.w << ",";
    os << data.q2 << ",";
    // mPim
    // os << std::setprecision(7);
    // os << data.pim_mom_mPim << ",";
    // os << data.pim_theta_mPim << ",";
    // os << data.pim_phi_mPim << ",";

    // os << data.pim_mom_mPim_cm << ",";
    // os << data.pim_theta_mPim_cm << ",";
    // os << data.pim_phi_mPim_cm << ",";
    // os << std::setprecision(7);
    // os << data.mm2_mPim << ",";
    // //  os << data.mm2_mPim_corr << ",";
    // // os << std::setprecision(1);

    // // os << data.status_Pim << ",";
    // // os << data.status_Pip << ",";
    // // os << data.status_Prot << ",";
    os << std::setprecision(7);
    os << data.weight_mPim << ",";

    // // mPip
    // os << std::setprecision(7);
    // os << data.pip_mom_mPip << ",";
    // os << data.pip_theta_mPip << ",";
    // os << data.pip_phi_mPip << ",";
    // os << std::setprecision(7);
    // os << data.mm2_mPip << ",";
    // //  os << data.mm2_mPip_corr << ",";
    // // os << std::setprecision(1);
    // os << std::setprecision(7);
    // os << data.weight_mPip << ",";

    // // mProt
    // os << std::setprecision(7);
    // os << data.prot_mom_mProt << ",";
    // os << data.prot_theta_mProt << ",";
    // os << data.prot_phi_mProt << ",";
    // os << std::setprecision(7);
    // os << data.mm2_mProt << ",";
    // //  os << data.mm2_mProt_corr << ",";
    // // os << std::setprecision(1);
    // os << std::setprecision(7);
    // os << data.weight_mProt << ",";

    return os;
  }
};

#endif
