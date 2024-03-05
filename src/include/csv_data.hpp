#ifndef CSV_DATA_H_GUARD
#define CSV_DATA_H_GUARD

#include <fstream>
#include <string>

struct csv_data {
  short electron_sector;

  float w;
  float q2;
  float w_mc;
  float q2_mc;
  float elec_mom;
  float elec_energy;
  float elec_theta;
  float sf;

  float pim_mom_mPim;
  float pim_theta_mPim;
  float pim_phi_mPim;

  float pim_mom_mPim_cm;
  float pim_theta_mPim_cm;
  float pim_phi_mPim_cm;

  float mm2_mPim;
  float mm2_mPim_corr;
  float weight_mPim_rec;
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

  float prot_mom_mes;
  float prot_theta_mes;
  float prot_phi_mes;

  float inv_ppip;
  float inv_ppim;
  float inv_pip_pim;

  double residualXpcal;
  double residualYpcal;
  double residualZpcal;

  double residualXecin;
  double residualYecin;
  double residualZecin;

  double Xpcal;
  double Ypcal;

  double Xecin;
  double Yecin;

  double Xpcal_rot;
  double Ypcal_rot;

  double Xecin_rot;
  double Yecin_rot;
  // Static functions can be called without making a new struct
  static std::string header() {
    // Make a string for the header of the csv file mPim case
    // return "elec_sec,sf,elec_mom,ResiXpcal,ResiYpcal,pcalX,pcalY,pcalX_rot,pcalY_rot,ResiXecin,"
    //        "ResiYecin,ecinX,ecinY,ecinX_rot,ecinY_rot,weight";
    // return "w_mc,q2_mc,w_ec,q2_rec,mm2_mPim,weight";

    // return "pim_mom_mPim_cm,pim_theta_mPim_cm,pim_phi_mPim_cm,mm2_mPim,weight";
    return "w_rec,q2_rec,mm2_mPim,weight";
    // return "pip_mom_mPip,pip_theta_mPip,pip_phi_mPip,mm2_mPip,weight";
    // return "prot_mom_mProt,prot_theta_mProt,prot_phi_mProt,mm2_mProt,"
    //        "weight";
  }

  friend std ::ostream &operator<<(std::ostream &os, const csv_data &data) {
    // os << std::setprecision(1);

    // os << data.electron_sector << ",";
    os << std::setprecision(7);
    // os << data.sf << ",";
    // os << data.elec_mom << ",";
    // // os << data.elec_energy << ",";
    // // os << data.elec_theta << ",";

    // os << data.w_mc << ",";
    // os << data.q2_mc << ",";
    os << data.w << ",";
    os << data.q2 << ",";

    // // os << data.residualXpcal << ",";
    // // os << data.residualYpcal << ",";
    // // // os << data.residualZpcal << ",";

    // // os << data.Xpcal << ",";
    // // os << data.Ypcal << ",";

    // // os << data.Xpcal_rot << ",";
    // // os << data.Ypcal_rot << ",";

    // // os << data.residualXecin << ",";
    // // os << data.residualYecin << ",";
    // // // os << data.residualZecin << ",";
    // // os << data.Xecin << ",";
    // // os << data.Yecin << ",";

    // // os << data.Xecin_rot << ",";
    // // os << data.Yecin_rot << ",";
    // // // // mPim
    // os << data.pim_mom_mPim << ",";
    // os << data.pim_theta_mPim << ",";
    // os << data.pim_phi_mPim << ",";

    // // // // os << data.pim_mom_mPim_cm << ",";
    // // // // os << data.pim_theta_mPim_cm << ",";
    // // // // os << data.pim_phi_mPim_cm << ",";
    // // // // os << std::setprecision(7);
    os << data.mm2_mPim << ",";
    // // // //  os << data.mm2_mPim_corr << ",";
    // // // // os << std::setprecision(1);

    // // // os << data.status_Pim << ",";
    // // // os << data.status_Pip << ",";
    // // // os << data.status_Prot << ",";

    // // // mPip
    // // os << std::setprecision(7);
    // // os << data.pip_mom_mPip << ",";
    // // os << data.pip_theta_mPip << ",";
    // // os << data.pip_phi_mPip << ",";
    // // os << std::setprecision(7);
    // // os << data.mm2_mPip << ",";
    // // //  os << data.mm2_mPip_corr << ",";
    // // // os << std::setprecision(1);
    // // os << std::setprecision(7);
    // // os << data.weight_mPip << ",";

    // // // // mProt
    // // os << std::setprecision(7);
    // // os << data.prot_mom_mProt << ",";
    // // os << data.prot_theta_mProt << ",";
    // // os << data.prot_phi_mProt << ",";

    // // // os << data.prot_mom_mes << ",";
    // // // os << data.prot_theta_mes << ",";
    // // // os << data.prot_phi_mes << ",";

    // // os << std::setprecision(7);
    // // os << data.mm2_mProt << ",";
    // // //  os << data.mm2_mProt_corr << ",";
    // // // os << std::setprecision(1);
    // // os << std::setprecision(7);
    // // os << data.weight_mProt << ",";

    // // // os << data.inv_ppip << ",";
    // // // os << data.inv_ppim << ",";
    // // // os << data.inv_pip_pim << ",";
    // // // os << std::setprecision(7);
    // os << data.weight_mPim_rec << ",";
    os << data.weight_mPim << ",";

    return os;
  }
};

#endif
