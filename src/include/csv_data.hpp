#ifndef CSV_DATA_H_GUARD
#define CSV_DATA_H_GUARD

#include <fstream>
#include <string>

struct csv_data {
  short electron_sector;
  short pim_sec;
  short pip_sec;
  short prot_sec;
  float w;
  float q2;

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
  float pip_mom_mes;
  float pip_theta_mes;
  float pip_phi_mes;

  float prot_mom_mProt;
  float prot_theta_mProt;
  float prot_phi_mProt;
  float mm2_mProt;
  float mm2_mProt_corr;
  float weight_mProt;
  float prot_mom_mes;
  float prot_theta_mes;
  float prot_phi_mes;
  float prot_mom_corr;

  float inv_ppip;
  float inv_ppim;
  float inv_pip_pim;

  float prot_eff;
  float pip_eff;

  // Static functions can be called without making a new struct
  static std::string header() {
    return "w_rec,q2_rec,mm2_mPim_corr,weight";
    // return "elec_sec,w,q2,pim_mom_mPim,pim_theta_mPim,pim_phi_mPim,mm2_mPim,mm2_mPim_corr,inv_pPip,weight";
    // Make a string for the header of the csv file mPim case
    // return
    // "sec_pim,sec_pip,sec_prot,pim_mom_mPim,pim_theta_mPim,pim_phi_mPim,mm2_mPim,mm2_mPim_corr,status_Pim,status_"
    //        "Pip,status_Prot,pip_mom_mes,pip_theta_mes,pip_phi_"
    //        "mes,prot_mom_mes,prot_theta_mes,prot_phi_mes,prot_eff,pip_eff,weight";
    // return "pip_mom_mPip,pip_theta_mPip,pip_phi_mPip,mm2_mPip,mm2_mPip_corr,weight";
    // return
    // "prot_mom_mProt,prot_theta_mProt,prot_phi_mProt,prot_mom_mes,prot_theta_mes,prot_phi_mes,prot_mom_corr,mm2_"
    //        "mProt,mm2_mProt_corr,weight";
  }

  friend std ::ostream &operator<<(std::ostream &os, const csv_data &data) {
    // // // // mPim
    // os << std::setprecision(1);
    // // os << data.electron_sector << ",";
    // os << data.pim_sec << ",";
    // os << data.pip_sec << ",";
    // os << data.prot_sec << ",";

    os << std::setprecision(7);

    os << data.w << ",";
    os << data.q2 << ",";

    // os << data.pim_mom_mPim << ",";
    // os << data.pim_theta_mPim << ",";
    // os << data.pim_phi_mPim << ",";
    // // os << data.pim_mom_mPim_cm << ",";
    // // os << data.pim_theta_mPim_cm << ",";
    // // os << data.pim_phi_mPim_cm << ",";
    // // os << std::setprecision(7);
    // os << data.mm2_mPim << ",";
    os << data.mm2_mPim_corr << ",";

    // os << std::setprecision(1);
    // os << data.status_Pim << ",";
    // os << data.status_Pip << ",";
    // os << data.status_Prot << ",";

    // // os << std::setprecision(7);
    // // os << data.inv_ppip << ",";
    // // os << data.inv_ppim << ",";
    // // os << data.inv_pip_pim << ",";

    // // // mPip
    // os << std::setprecision(7);
    // // os << data.pip_mom_mPip << ",";
    // // os << data.pip_theta_mPip << ",";
    // // os << data.pip_phi_mPip << ",";
    // os << data.pip_mom_mes << ",";
    // os << data.pip_theta_mes << ",";
    // os << data.pip_phi_mes << ",";
    // // os << data.mm2_mPip << ",";
    // // os << data.mm2_mPip_corr << ",";
    // // os << std::setprecision(1);
    // // os << data.weight_mPip << ",";

    // // // mProt
    // // os << std::setprecision(7);
    // //  os << data.prot_mom_mProt << ",";
    // //  os << data.prot_theta_mProt << ",";
    // //  os << data.prot_phi_mProt << ",";
    // os << data.prot_mom_mes << ",";
    // os << data.prot_theta_mes << ",";
    // os << data.prot_phi_mes << ",";
    // //  os << data.prot_mom_corr << ",";
    // //  os << data.mm2_mProt << ",";
    // //  os << data.mm2_mProt_corr << ",";
    // //  os << std::setprecision(1);
    // //  os << data.weight_mProt << ",";

    // os << std::setprecision(7);
    // os << data.prot_eff << ",";
    // os << data.pip_eff << ",";
    os << std::setprecision(1);
    os << data.weight_mPim << ",";

    return os;
  }
};

#endif
