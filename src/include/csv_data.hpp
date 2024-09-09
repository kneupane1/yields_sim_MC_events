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
  float w_mc;
  float q2_mc;
  float weight_mc;
  float mm2_mPim_mc;
  float w_had;
  float w_diff;
  float w_had_corr;
  float w_diff_corr;
  float w_after;

  float elec_mom;
  float elec_energy;
  float elec_theta;
  float elec_mom_mc;
  float elec_energy_mc;
  float elec_theta_mc;
  float elec_mom_rec;
  float elec_energy_rec;
  float elec_theta_rec;

  float prot_mom_mc;
  float prot_theta_mc;
  float prot_mom_rec;
  float prot_theta_rec;

  float pip_mom_mc;
  float pip_theta_mc;
  float pip_mom_rec;
  float pip_theta_rec;

  float pim_mom_mc;
  float pim_theta_mc;
  float pim_mom_rec;
  float pim_theta_rec;

  float corr_elec_mom;

  float scalar_product;
  float pim_mom_mPim;
  float pim_theta_mPim;
  float pim_phi_mPim;
  float pim_mom_mPim_cm;
  float pim_theta_mPim_cm;
  float pim_phi_mPim_cm;

  float mm2_mPim;
  float mm2_mPim_corr;
  float weight_mPim;
  float pim_mom_exclusive;

  float pim_mom_corr;
  float pim_theta_corr;
  float pim_phi_corr;

  float pip_mom_corr;
  float pip_theta_corr;
  float pip_phi_corr;

  float prot_mom_corr;
  float prot_theta_corr;
  float prot_phi_corr;

  float pim_theta_exclusive;
  float pim_phi_exclusive;
  float mm2_exclusive;
  float mm2_exclusive_at_zero;
  float weight_exclusive;

  float pip_mom_mPip;
  float pip_theta_mPip;
  float pip_phi_mPip;
  float mm2_mPip;
  float mm2_mPip_corr;
  float weight_mPip;
  float pip_mom_exclusive;
  float pip_theta_exclusive;
  float pip_phi_exclusive;
  float energy_x_mu;

  float prot_mom_mProt;
  float prot_theta_mProt;
  float prot_phi_mProt;
  float mm2_mProt;
  float mm2_mProt_corr;
  float weight_mProt;

  float prot_mom_exclusive;
  float prot_theta_exclusive;
  float prot_phi_exclusive;
  float prot_dcr1theta_exclusive;
  float pip_dcr1theta_exclusive;
  float pim_dcr1theta_exclusive;

  float diff_ex_theta;
  float diff_ex_phi;
  float diff_bx_theta;
  float diff_bx_phi;

  float x_mu_mom_exclusive;
  float x_mu_theta_exclusive;
  float x_mu_phi_exclusive;

  // float diff_rec_mes_pim_mom;
  // float diff_rec_mes_pim_theta;
  // float diff_rec_mes_pim_phi;

  // float diff_gen_pim_mom;
  // float diff_gen_pim_theta;
  // float diff_gen_pim_phi;

  float gen_pim_mom;
  float gen_pim_theta;
  float gen_pim_phi;

  float gen_pip_mom;
  float gen_pip_theta;
  float gen_pip_phi;

  float gen_prot_mom;
  float gen_prot_theta;
  float gen_prot_phi;

  // float diff_rec_mes_pip_mom;
  // float diff_rec_mes_pip_theta;
  // float diff_rec_mes_pip_phi;

  // float diff_rec_mes_prot_mom;
  // float diff_rec_mes_prot_theta;
  // float diff_rec_mes_prot_phi;

  int status_Pim;
  int status_Pip;
  int status_Prot;

  float beta_Prot;
  float beta_Pip;

  float inv_ppip;
  float inv_ppim;
  float inv_pip_pim;
  int prot_pid_mc, prot_pid_rec, pip_pid_mc, pip_pid_rec, pim_pid_rec;

  // Static functions can be called without making a new struct
  static std::string header() {
    // Make a string for the header of the csv file mPim case
    // return "w_mc,q2_mc,mm2_mPim_mc,weight";

    return "prot_pid_rec,pip_pid_rec,pim_pid_rec,w_rec,q2_rec,mm2_mProt,mm2_mPip,mm2_mPim,status_Pip,status_Prot,"
           "weight";
  }

  friend std ::ostream &operator<<(std::ostream &os, const csv_data &data) {
    ////.......................................
    os << std::setprecision(1);

    // os << data.prot_pid_mc << ",";
    os << data.prot_pid_rec << ",";
    // os << data.pip_pid_mc << ",";
    os << data.pip_pid_rec << ",";
    os << data.pim_pid_rec << ",";

    // // For excl case

    // // //  // os << data.electron_sector << ",";
    // os << data.pim_sec << ",";
    // os << data.pip_sec << ",";
    // os << data.prot_sec << ",";

    os << std::setprecision(7);

    // os << data.w_mc << ",";
    // os << data.q2_mc << ",";
    // os << data.mm2_mPim_mc << ",";

    // os << data.weight_mc << ",";

    /////////////////////////////////////////

    os << data.w << ",";
    os << data.q2 << ",";
    // // //  // // // os << data.w_after << ",";

    // // //  // os << data.w_had << ",";
    // // //  // // // // os << data.w_diff << ",";
    // // //  // // // os << data.w_had_corr << ",";
    // // //  // // // // os << data.w_diff_corr << ",";

    // // //  // // // // os << data.w_after << ",";
    // // //  // // os << data.elec_mom << ",";
    // // //  // os << data.elec_energy << ",";
    // // //  // os << data.elec_theta << ",";

    // os << data.elec_mom_mc << ",";
    // // os << data.elec_energy_mc << ",";
    // os << data.elec_theta_mc << ",";

    // os << data.elec_mom_rec << ",";
    // // os << data.elec_energy_rec << ",";
    // os << data.elec_theta_rec << ",";

    // // //  //  // // os << data.corr_elec_mom << ",";
    // os << data.scalar_product << ",";
    // //  //  // // // // Generated
    // //  //  // // // os << std::setprecision(5);

    // os << data.gen_prot_mom << ",";
    // os << data.gen_prot_theta << ",";
    // os << data.gen_prot_phi << ",";

    // os << data.gen_pip_mom << ",";
    // os << data.gen_pip_theta << ",";
    // os << data.gen_pip_phi << ",";

    // os << data.gen_pim_mom << ",";
    // os << data.gen_pim_theta << ",";
    // os << data.gen_pim_phi << ",";

    // // //  // // // Missing
    // os << data.prot_mom_mProt << ",";
    // os << data.prot_theta_mProt << ",";
    // os << data.prot_phi_mProt << ",";

    // os << data.pip_mom_mPip << ",";
    // os << data.pip_theta_mPip << ",";
    // os << data.pip_phi_mPip << ",";

    // os << data.pim_mom_mPim << ",";
    // os << data.pim_theta_mPim << ",";
    // os << data.pim_phi_mPim << ",";

    // // // //  // // // measured

    // os << data.prot_mom_exclusive << ",";
    // os << data.prot_theta_exclusive << ",";
    // os << data.prot_phi_exclusive << ",";
    // // os << data.prot_dcr1theta_exclusive << ",";

    // os << data.pip_mom_exclusive << ",";
    // os << data.pip_theta_exclusive << ",";
    // os << data.pip_phi_exclusive << ",";
    // // os << data.pip_dcr1theta_exclusive << ",";

    // os << data.pim_mom_exclusive << ",";
    // os << data.pim_theta_exclusive << ",";
    // os << data.pim_phi_exclusive << ",";
    // // os << data.pim_dcr1theta_exclusive << ",";

    os << data.mm2_mProt << ",";
    os << data.mm2_mPip << ",";
    os << data.mm2_mPim << ",";
    // os << data.mm2_mPim_mc << ",";

    os << data.mm2_exclusive_at_zero << ",";
    os << data.energy_x_mu << ",";

    // os << std::setprecision(1);

    os << data.status_Pim << ",";
    os << data.status_Pip << ",";
    os << data.status_Prot << ",";

    // os << data.beta_Pip << ",";
    // os << data.beta_Prot << ",";
    // os << std::setprecision(7);
    // os << data.inv_ppip << ",";
    // os << data.inv_ppim << ",";
    // os << data.inv_pip_pim << ",";

    // os << std::setprecision(1);

    os << data.weight_exclusive << ",";

    //  ///.......................................

    // // for pim eff check

    //      os << std::setprecision(7);

    //     //  os << data.w << ",";
    //     //  os << data.q2 << ",";
    //      os << data.scalar_product << ",";

    // //      //missing
    // //      os << data.pim_mom_mPim << ",";
    // //      os << std::setprecision(5);
    // //      os << data.pim_theta_mPim << ",";
    // //      os << data.pim_phi_mPim << ",";
    // os << data.pim_mom_mPim_cm << ",";
    // os << data.pim_theta_mPim_cm << ",";
    // os << data.pim_phi_mPim_cm << ",";

    //     //  // // // measured
    //     //  os << std::setprecision(7);
    //     //  os << data.pim_mom_exclusive << ",";
    //     //  os << std::setprecision(5);
    //     //  os << data.pim_theta_exclusive << ",";
    //     //  os << data.pim_phi_exclusive << ",";
    //      os << data.mm2_mProt << ",";
    //      os << data.mm2_mPip << ",";
    //      os << data.mm2_mPim << ",";
    //      // os << data.mm2_mPim_corr << ",";
    //      os << std::setprecision(7);
    //      os << data.mm2_exclusive_at_zero << ",";
    //      os << data.energy_x_mu << ",";
    //      os << data.weight_exclusive << ",";

    //  ///.......................................
    return os;
  }
};

#endif
