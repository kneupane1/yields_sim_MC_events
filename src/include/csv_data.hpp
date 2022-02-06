#ifndef CSV_DATA_H_GUARD
#define CSV_DATA_H_GUARD

#include <fstream>
#include <string>

struct csv_data {
  short electron_sector;
  float w;
  float w_after;
  float q2;
  float scalar_product;
  float pim_mom_mPim;
  float pim_theta_mPim;
  float pim_phi_mPim;
  float mm2_mPim;
  float weight_mPim;
  float pim_mom_exclusive;
  float pim_theta_exclusive;
  float pim_phi_exclusive;
  float mm2_exclusive;
  float mm2_exclusive_at_zero;
  float weight_exclusive;

  float pip_mom_mPip;
  float pip_theta_mPip;
  float pip_phi_mPip;
  float mm2_mPip;
  float weight_mPip;
  float pip_mom_exclusive;
  float pip_theta_exclusive;
  float pip_phi_exclusive;
  float energy_x_mu;

  float prot_mom_mProt;
  float prot_theta_mProt;
  float prot_phi_mProt;
  float mm2_mProt;
  float weight_mProt;
  float prot_mom_exclusive;
  float prot_theta_exclusive;
  float prot_phi_exclusive;

  float diff_ex_theta;
  float diff_ex_phi;
  float diff_bx_theta;
  float diff_bx_phi;

  float x_mu_mom_exclusive;
  float x_mu_theta_exclusive;
  float x_mu_phi_exclusive;

  float diff_rec_mes_pim_mom;
  float diff_rec_mes_pim_theta;
  float diff_rec_mes_pim_phi;

  int status_Pim;
  int status_Pip;
  int status_Prot;

  float elec_mom;
  float corr_elec_mom;

  // Static functions can be called without making a new struct
  static std::string header() {
    // Make a string for the header of the csv file mPim case
    // return "w,pim_mom_mPim,pim_theta_mPim,pim_phi_mPim,mm2_mPim,weight";
    // return "w,,stp,pim_mom_exclusive,pim_theta_exclusive,pim_phi_exclusive,mm2_exclusive,weight";
    // return
    // "w,stp,pim_mom_exclusive,pim_theta_exclusive,pim_phi_exclusive,mm2_exclusive,mm2_exclusive_at_zero,energy_x_mu,diff_ex_theta,diff_ex_phi,diff_bx_theta,diff_bx_phi,weight";
    // return
    // "w,x_mu_mom_exclusive,x_mu_theta_exclusive,x_mu_phi_exclusive,mm2_exclusive_at_zero,energy_x_mu,diff_ex_theta,diff_ex_phi,diff_bx_theta,diff_bx_phi,weight";

    // for mom thete phi rec- mes check
    return "sec_ele,w_after,pim_mom_mPim,pim_theta_mPim,pim_phi_mPim,mm2_mPim,status_Pim,status_Pip,status_Prot,weight";

    // mPip case
    // return "w,pip_mom_mPip,pip_theta_mPip,pip_phi_mPip,mm2_mPip,weight";
    // return "w,stp,pip_mom_exclusive,pip_theta_exclusive,pip_phi_exclusive,mm2_exclusive,weight";

    // // mProt case
    // return "w,prot_mom_mProt,prot_theta_mProt,prot_phi_mProt,mm2_mProt,weight";
    // return "w,stp,prot_mom_exclusive,prot_theta_exclusive,prot_phi_exclusive,mm2_exclusive,weight";
  }

  friend std ::ostream &operator<<(std::ostream &os, const csv_data &data) {
    os << std::setprecision(3);
    os << data.electron_sector<< ",";
    // os << data.w << ",";
    os << data.w_after << ",";

    // mPim
    /*     os << data.pim_mom_mPim << ",";
         os << data.pim_theta_mPim << ",";
         os << data.pim_phi_mPim << ",";
         os << data.mm2_mPim << ",";
         os << std::setprecision(5);
         os << data.weight_mPim << ",";


     */
    // os << data.elec_mom << ",";
    // os << data.corr_elec_mom << ",";

    os << data.pim_mom_mPim << ",";
    os << data.pim_theta_mPim << ",";
    os << data.pim_phi_mPim << ",";
    os << data.mm2_mPim << ",";

    // os << data.scalar_product << ",";
    // os << data.pim_mom_exclusive << ",";
    // os << data.pim_theta_exclusive << ",";
    // os << data.pim_phi_exclusive << ",";
    // // os << data.mm2_exclusive << ",";

    // os << data.mm2_exclusive_at_zero << ",";
    // os << data.energy_x_mu << ",";
    // os << data.mm2_mPip << ",";
    // os << data.mm2_mProt << ",";

    // os << data.diff_rec_mes_pim_mom << ",";
    // os << data.diff_rec_mes_pim_theta << ",";
    // os << data.diff_rec_mes_pim_phi << ",";

    // // os << data.diff_ex_theta << ",";
    // // os << data.diff_ex_phi << ",";
    // // os << data.diff_bx_theta << ",";
    // // os << data.diff_bx_phi << ",";

    os << data.status_Pim << ",";
    os << data.status_Pip << ",";
    os << data.status_Prot << ",";

    os << std::setprecision(5);
    os << data.weight_exclusive << ",";

    // os << data.x_mu_mom_exclusive << ",";
    // os << data.x_mu_theta_exclusive << ",";
    // os << data.x_mu_phi_exclusive << ",";
    // // os << data.mm2_exclusive << ",";

    // os << data.mm2_exclusive_at_zero << ",";
    // os << data.energy_x_mu << ",";

    // os << data.diff_ex_theta << ",";
    // os << data.diff_ex_phi << ",";
    // os << data.diff_bx_theta << ",";
    // os << data.diff_bx_phi << ",";

    // os << std::setprecision(5);
    // os << data.weight_exclusive<<",";

    // mPip
    /*  os << data.pip_mom_mPip << ",";
      os << data.pip_theta_mPip << ",";
      os << data.pip_phi_mPip << ",";
      os << data.mm2_mPip << ",";
      os << std::setprecision(5);
      os << data.weight_mPip << ",";
      */

    /*  os << data.scalar_product << ",";
      os << data.pip_mom_exclusive << ",";
      os << data.pip_theta_exclusive << ",";
      os << data.pip_phi_exclusive << ",";
      os << data.mm2_exclusive << ",";
      os << std::setprecision(5);
      os << data.weight_exclusive << ",";
  */
    // mProt
    /* os << data.prot_mom_mProt << ",";
     os << data.prot_theta_mProt << ",";
     os << data.prot_phi_mProt << ",";
     os << data.mm2_mProt << ",";
         os << std::setprecision(5);
     os << data.weight_mProt << ",";
*/
    /*      os << data.scalar_product << ",";
          os << data.prot_mom_exclusive << ",";
          os << data.prot_theta_exclusive << ",";
          os << data.prot_phi_exclusive << ",";
          os << data.mm2_exclusive << ",";
              os << std::setprecision(5);
          os << data.weight_exclusive << ",";
  */
    return os;
  }
};

#endif
