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

  float corr_elec_mom;


  float scalar_product;
  float pim_mom_mPim;
  float pim_theta_mPim;
  float pim_phi_mPim;
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

  float inv_ppip;
  float inv_ppim;
  float inv_pip_pim;

  // Static functions can be called without making a new struct
  static std::string header() {
    // Make a string for the header of the csv file mPim case
    // return "w,q2,elec_en,elec_theta,weight";
    // return "w,q2,elec_mom,elec_en,elec_theta,w_mc,q2_mc,elec_mom_mc,elec_en_mc,elec_theta_mc,weight";

    // return "pim_mom_mPim,pim_theta_mPim,pim_phi_mPim,mm2_mPim,mm2_mPim_corr,weight";
    // return "pim_mom_mPim,pim_theta_mPim,pim_phi_mPim,mm2_mPim,weight";


    // for all 3 particles // // Data/ Simulations excl....    prot_mom_gen,pip_mom_gen,pim_mom_gen,
    // return "prot_mom_gen,prot_theta_gen,prot_phi_gen,pip_mom_gen,pip_theta_gen,pip_phi_gen,pim_mom_gen,pim_theta_gen,"
    //        "pim_phi_gen,prot_mom_mes,prot_theta_mes,prot_phi_mes,"
    //        "pip_mom_mes,pip_theta_mes,pip_phi_mes,pim_mom_mes,pim_theta_mes,pim_phi_mes,"
    //        "mm2_exclusive_at_zero,energy_x_mu,status_Pim,"
    //        "status_Pip,status_Prot,weight";
    // return "prot_sec,prot_mom_gen,prot_theta_gen,prot_phi_gen,"
    //        "prot_mom_mes,prot_theta_mes,prot_phi_mes,dcr1_theta_prot,"
    //        "mm2_exclusive_at_zero,energy_x_mu,status_Prot,weight";

    // return "pip_sec,pip_mom_gen,pip_theta_gen,pip_phi_gen,"
    //        "pip_mom_mes,pip_theta_mes,pip_phi_mes,dcr1_theta_pip,"
    //        "status_Pip,weight";
    // return "pim_sec,pim_mom_gen,pim_theta_gen,pim_phi_gen,"
    //        "pim_mom_mes,pim_theta_mes,pim_phi_mes,dcr1_theta_pim,"
    //        "status_Pim,weight";
    return "sec_pim,sec_pip,sec_prot,prot_mom_miss,prot_theta_miss,prot_phi_"
           "miss,pip_mom_miss,pip_"
           "theta_miss,pip_phi_miss,pim_mom_miss,pim_theta_miss,pim_phi_miss,prot_mom_mes,prot_theta_mes,prot_phi_mes,"
           "dcr1_theta_prot,prot_mom_corr,"
           "pip_mom_mes,pip_theta_mes,pip_phi_"
           "mes,dcr1_theta_pip,pip_mom_corr,pim_mom_mes,pim_theta_mes,pim_phi_mes,dcr1_theta_pim,pim_mom_corr,mm2_"
           "mProt,mm2_mProt_corr,mm2_mPip,mm2_mPip_corr,mm2_mPim,mm2_mPim_corr,mm2_"
           "exclusive_at_zero,energy_x_mu,"
           "status_Pim,"
           "status_Pip,status_Prot,weight";

    // mPip case
    // return "pip_mom_mPip,pip_theta_mPip,pip_phi_mPip,mm2_mPip,mm2_mPip_corr,weight";
    // return "w,stp,pip_mom_exclusive,pip_theta_exclusive,pip_phi_exclusive,mm2_exclusive,weight";

    // // mProt case
    // return "prot_mom_mProt,prot_theta_mProt,prot_phi_mProt,mm2_mProt,mm2_mProt_corr,weight";
    // return "w,stp,prot_mom_exclusive,prot_theta_exclusive,prot_phi_exclusive,mm2_exclusive,weight";

    //for cross-section checks
    // return "w,q2,w_had,w_had_corr_1st_iter,mm2_exclusive_at_zero,energy_x_mu,weight";
    // return "w,q2,w_had,w_mc,q2_mc,mm2_exclusive_at_zero,energy_x_mu,weight";
    // return "w_mc,q2_mc,weight";
  }

  friend std ::ostream &operator<<(std::ostream &os, const csv_data &data) {
    ////.......................................
     os << std::setprecision(1);

    // os << data.electron_sector << ",";
    os << data.pim_sec << ",";
    os << data.pip_sec << ",";
    os << data.prot_sec << ",";

    os << std::setprecision(7);

    // os << data.w << ",";
    // os << data.q2 << ",";
    // // // os << data.w_after << ",";

    // // os << data.w_had << ",";
    // // // // os << data.w_diff << ",";
    // // // os << data.w_had_corr << ",";
    // // // // os << data.w_diff_corr << ",";


    // // // // os << data.w_after << ",";
    // // os << data.elec_mom << ",";
    // os << data.elec_energy << ",";
    // os << data.elec_theta << ",";

    // os << data.w_mc << ",";
    // os << data.q2_mc << ",";
    // os << data.elec_mom_mc << ",";
    // os << data.elec_energy_mc << ",";
    // os << data.elec_theta_mc << ",";

    // // os << data.corr_elec_mom << ",";
    // os << data.scalar_product << ",";
    // // // // Generated
    // // // os << std::setprecision(5);

    // os << data.gen_prot_mom << ",";
    // os << data.gen_prot_theta << ",";
    // os << data.gen_prot_phi << ",";

    // os << data.gen_pip_mom << ",";
    // os << data.gen_pip_theta << ",";
    // os << data.gen_pip_phi << ",";

    // os << data.gen_pim_mom << ",";
    // os << data.gen_pim_theta << ",";
    // os << data.gen_pim_phi << ",";

    // // // Missing
    os << std::setprecision(10);
    os << data.prot_mom_mProt << ",";
    os << std::setprecision(7);

    os << data.prot_theta_mProt << ",";
    os << data.prot_phi_mProt << ",";
    os << std::setprecision(10);

    os << data.pip_mom_mPip << ",";
    os << std::setprecision(7);

    os << data.pip_theta_mPip << ",";
    os << data.pip_phi_mPip << ",";
    os << std::setprecision(10);

    os << data.pim_mom_mPim << ",";
    os << std::setprecision(7);

    os << data.pim_theta_mPim << ",";
    os << data.pim_phi_mPim << ",";

    // // // measured
    os << std::setprecision(10);

    os << data.prot_mom_exclusive << ",";
    os << std::setprecision(7);

    os << data.prot_theta_exclusive << ",";
    os << data.prot_phi_exclusive << ",";
    os << data.prot_dcr1theta_exclusive << ",";
    os << std::setprecision(10);

    os << data.prot_mom_corr << ",";
    // os << data.prot_theta_corr << ",";
    // os << data.prot_phi_corr << ",";
    os << std::setprecision(10);

    os << data.pip_mom_exclusive << ",";
    os << std::setprecision(7);

    os << data.pip_theta_exclusive << ",";
    os << data.pip_phi_exclusive << ",";
    os << data.pip_dcr1theta_exclusive << ",";

    os << std::setprecision(10);
    os << data.pip_mom_corr << ",";
    // os << data.pip_theta_corr << ",";
    // os << data.pip_phi_corr << ",";
    // os << std::setprecision(10);

    os << data.pim_mom_exclusive << ",";
    // os << std::setprecision(7);

    os << data.pim_theta_exclusive << ",";
    os << data.pim_phi_exclusive << ",";
    os << data.pim_dcr1theta_exclusive << ",";

    os << std::setprecision(10);
    os << data.pim_mom_corr << ",";
    // // // os << data.pim_theta_corr << ",";
    // // // os << data.pim_phi_corr << ",";
    // // // os << std::setprecision(10);


    os << data.mm2_mProt << ",";
    os << data.mm2_mProt_corr << ",";

    os << data.mm2_mPip << ",";
    os << data.mm2_mPip_corr << ",";

    os << data.mm2_mPim << ",";
    os << data.mm2_mPim_corr << ",";
    os << std::setprecision(7);
    os << data.mm2_exclusive_at_zero << ",";
    os << data.energy_x_mu << ",";
    os << std::setprecision(1);

    os << data.status_Pim << ",";
    os << data.status_Pip << ",";
    os << data.status_Prot << ",";

    // os << std::setprecision(7);
    // os << data.inv_ppip << ",";
    // os << data.inv_ppim << ",";
    // os << data.inv_pip_pim << ",";

    // os << std::setprecision(1);

    os << data.weight_exclusive << ",";

    ///.......................................

    // os << data.pim_mom_mPim << ",";
    // os << data.pim_mom_exclusive << ",";
    // os << data.pim_mom_corr << ",";

    // os << data.pim_theta_mPim << ",";
    // os << data.pim_theta_exclusive << ",";
    // os << data.pim_theta_corr << ",";

    // os << data.pim_phi_mPim << ",";
    // os << data.pim_phi_exclusive << ",";
    // os << data.pim_phi_corr << ",";

    // os << data.pip_mom_mPip << ",";
    // os << data.pip_mom_exclusive << ",";
    // os << data.pip_mom_corr << ",";

    // os << data.pip_theta_mPip << ",";
    // os << data.pip_theta_exclusive << ",";
    // os << data.pip_theta_corr << ",";

    // os << data.pip_phi_mPip << ",";
    // os << data.pip_phi_exclusive << ",";
    // os << data.pip_phi_corr << ",";

    // os << data.prot_mom_mProt << ",";
    // os << data.prot_mom_exclusive << ",";
    // os << data.prot_mom_corr << ",";

    // os << data.prot_theta_mProt << ",";
    // os << data.prot_theta_exclusive << ",";
    // os << data.prot_theta_corr << ",";

    // os << data.prot_phi_mProt << ",";
    // os << data.prot_phi_exclusive << ",";
    // os << data.prot_phi_corr << ",";

    // os << data.mm2_exclusive_at_zero << ",";
    // os << data.energy_x_mu << ",";
    // os << data.weight_exclusive << ",";

    // // // // // // // // // // mPim
    //  os << data.pim_mom_mPim << ",";
    //  os << data.pim_theta_mPim << ",";
    //  os << data.pim_phi_mPim << ",";
    //  os << std::setprecision(10);
    //  os << data.mm2_mPim << ",";
    // //  os << data.mm2_mPim_corr << ",";
    //  os << std::setprecision(7);
    //  os << data.weight_mPim << ",";

    // // // // // // //
    // // // // // // // os << data.elec_mom << ",";
    // // // // // // // os << data.corr_elec_mom << ",";

    // // os << data.prot_mom_mProt << ",";
    // // os << data.prot_theta_mProt << ",";
    // // os << data.prot_phi_mProt << ",";
    // // os << data.mm2_mProt << ",";

    // os << data.scalar_product << ",";
    // os << data.pim_mom_exclusive << ",";
    // os << data.pim_theta_exclusive << ",";
    // os << data.pim_phi_exclusive << ",";
    // os << data.mm2_exclusive << ",";

    // os << data.mm2_exclusive_at_zero << ",";
    // os << data.energy_x_mu << ",";
    // // // os << data.mm2_mPip << ",";
    // // // os << data.mm2_mProt << ",";

    // // os << data.gen_prot_mom << ",";
    // // os << data.gen_prot_theta << ",";
    // // os << data.gen_prot_phi << ",";

    // // // os << data.diff_ex_theta << ",";
    // // // os << data.diff_ex_phi << ",";
    // // // os << data.diff_bx_theta << ",";
    // // // os << data.diff_bx_phi << ",";

    // // os << data.status_Pim << ",";
    // // os << data.status_Pip << ",";
    // // os << data.status_Prot << ",";

    // // // // os << std::setprecision(10);
    // // os << data.weight_exclusive << ",";

    // // os << data.x_mu_mom_exclusive << ",";
    // // os << data.x_mu_theta_exclusive << ",";
    // // os << data.x_mu_phi_exclusive << ",";
    // // // os << data.mm2_exclusive << ",";

    // // os << data.mm2_exclusive_at_zero << ",";
    // // os << data.energy_x_mu << ",";

    // // os << data.diff_ex_theta << ",";
    // // os << data.diff_ex_phi << ",";
    // // os << data.diff_bx_theta << ",";
    // // os << data.diff_bx_phi << ",";

    // // os << std::setprecision(10);
    // os << data.weight_exclusive<<",";

     // mPip .......................................
   /*    os << data.pip_mom_mPip << ",";
       os << data.pip_theta_mPip << ",";
       os << data.pip_phi_mPip << ",";
       os << data.mm2_mPip << ",";
       os << data.mm2_mPip_corr << ",";
       os << std::setprecision(1);
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
     // mProt .......................................
 /*    os << data.prot_mom_mProt << ",";
     os << data.prot_theta_mProt << ",";
     os << data.prot_phi_mProt << ",";
     os << data.mm2_mProt << ",";
     os << data.mm2_mProt_corr << ",";
     os << std::setprecision(1);
     os << data.weight_mProt << ",";
*/
    /*      os << data.scalar_product << ",";
          os << data.prot_mom_exclusive << ",";
          os << data.prot_theta_exclusive << ",";
          os << data.prot_phi_exclusive << ",";
          os << data.mm2_exclusive << ",";
              os << std::setprecision(10);
          os << data.weight_exclusive << ",";
  */
    return os;
  }
};

#endif
