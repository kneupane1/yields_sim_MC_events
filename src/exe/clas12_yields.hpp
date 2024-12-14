
#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
// #include "QADB.h"
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"
#include "cuts.hpp"
#include "histogram.hpp"
#include "reaction.hpp"
#include "syncfile.hpp"

/////////////////////////////////////////
template <class CutType>
size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<SyncFile>& _sync, int thread_id) {
  // size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<SyncFile>& _sync, const std::shared_ptr<QA::QADB>&
  // _qa,
  //            int thread_id) {
  // Get the number of events in this thread
  size_t num_of_events = (int)_chain->GetEntries();

  float beam_energy = 10.6041;
  if (std::is_same<CutType, Pass2_Cuts>::value) {
    beam_energy = 10.6041;
  }

  if (getenv("BEAM_E") != NULL) beam_energy = atof(getenv("BEAM_E"));

  // Print some information for each thread
  std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
            << num_of_events << " Events " << DEF << "===============\n";

  // Make a data object which all the branches can be accessed from
  // for sim data use it
  // auto data = std::make_shared<Branches12>(_chain, true);
  // for exp data use it
  auto data = std::make_shared<Branches12>(_chain, _mc);

  // Total number of events "Processed"
  size_t total = 0;
  int twoPion_excl = 0;
  int twoPion_mPim = 0;
  int twoPion_mPip = 0;
  int twoPion_mProt = 0;

  int numElec = 0;
  int numPip = 0;
  int numProt = 0;
  int numPim = 0;

  int numElec_mc = 0;
  int numPip_mc = 0;
  int numProt_mc = 0;
  int numPim_mc = 0;
  int Pip_pid_mc = -9999;
  int Prot_pid_mc = -9999;
  int Pim_pid_mc = -9999;
  int Pip_pid_rec = -9999;
  int Prot_pid_rec = -9999;
  int Pim_pid_rec = -9999;
  double dp_prot1 = NAN;
  double dp_pip1 = NAN;
  double dp_prot2 = NAN;
  double dp_pip2 = NAN;
  double dp_prot3 = NAN;
  double dp_pip3 = NAN;
  double dp_pim1 = NAN;
  double dp_pim2 = NAN;
  double dp_pim3 = NAN;

  // For each event
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    // for (size_t current_event = 220857; current_event < 220858; current_event++) {
    // Get current event
    _chain->GetEntry(current_event);

    // If we are the 0th thread print the progress of the thread every 1000 events
    if (thread_id == 0 && current_event % 1000 == 0)
      std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;

    // /////////////////////////////// Generated sim only //////////////////////////////////////
    // if (_mc) {
    if (data->mc_npart() < 1 || data->mc_weight() <= 0) continue;
    numElec_mc++;

    // If we pass electron cuts the event is processed
    total++;
    // std::cout << " event " << current_event << std::endl;

    // Make a reaction class from the data given
    auto mc_event = std::make_shared<MCReaction>(data, beam_energy);
    // if (mc_event->weight() <= 0.0) continue;
    // std::cout << " event " << current_event << std::endl;

    for (int part = 1; part < data->mc_npart(); part++) {
      // Check particle ID's and fill the reaction class

      if (data->mc_pid(part) == PIP) {
        numPip_mc++;
        // mc_event->SetMCProton(part);

        mc_event->SetMCPip(part);
      }
      if (data->mc_pid(part) == PROTON) {
        numProt_mc++;
        // mc_event->SetMCPip(part);

        mc_event->SetMCProton(part);
        // std::cout << mc_event->GetMcProtons().size() << std::endl;

      } else if (data->mc_pid(part) == PIM) {
        numPim_mc++;

        mc_event->SetMCPim(part);
        // } else {
        //   mc_event->SetMCOther(part);
      }
    }
    // // // Retrieve the number of protons and pions in the event
    // size_t num_protons_mc = mc_event->GetMcProtons().size();
    // size_t num_pips_mc = mc_event->GetMcPips().size();
    // for (size_t i = 0; i < num_protons_mc; ++i) {
    //   for (size_t j = 0; j < num_pips_mc; ++j) {
    //     mc_event->CalcMissMassPimMC(*mc_event->GetMcProtons()[i], *mc_event->GetMcPips()[j]);
    //   }
    // }
    // {
    //   csv_data output;

    // // // // for generated case
    // output.w_mc = mc_event->W_mc();
    // output.q2_mc = mc_event->Q2_mc();
    // output.mm2_mPim_mc = mc_event->MM2_mPim_MC();
    // output.weight_mc = mc_event->mc_weight();
    // _sync->write(output);
    // }

    /////////////////////////////// Reconstruction only //////////////////////////////////////
    auto event = std::make_shared<Reaction>(data, beam_energy);
    auto dt = std::make_shared<Delta_T>(data);
    auto cuts = std::make_shared<Pass2_Cuts>(data);
    // // auto cuts = std::make_shared<rga_Cuts>(data);
    // if (!_qa->Golden(data->getRun(), data->getEvent())) continue;

    if (!cuts->ElectronCuts()) continue;
    // std::cout << " chi2pid at 0 " << data->chi2pid(0) << std::endl;
    // event->SetMomCorrElec();

    numElec++;
    /*
        // // For each particle in the event
        for (int part = 0; part < data->gpart(); part++) {
          dt->dt_calc(part);
          // Check particle ID's and fill the reaction class
          if (cuts->IsPip(part)) {
            numPip++;
            event->SetPip(part);
            // event->SetProton(part);
          }
          if (cuts->IsProton(part)) {
            numProt++;
            no_of_protons_in_ths_event++;
            event->SetProton(part);
            // // event->SetPip(part);
          }
          if (cuts->IsPim(part)) {
            numPim++;
            event->SetPim(part);
          } else {
            event->SetOther(part);
          }
        }
        */
    // no_of_protons_in_ths_event = 0;

    // ///////////////////////////////////////  For  dp cut method ////////////////////
    // ///////////////////////////////////////  For  dp cut method ////////////////////
    // // ///////////////////////////////////////  For  dp cut method ////////////////////

    // // Define vectors to store dp values and indices for protons and pions
    std::vector<std::pair<int, double>> proton_dps;  // Pair of (index, dp_prot)
    std::vector<std::pair<int, double>> pip_dps;     // Pair of (index, dp_pip)

    for (int part = 1; part < data->gpart(); part++) {
      if (data->charge(part) != 0) {
        dt->dt_calc(part);

        if (data->charge(part) > 0) {
          // Check if the particle satisfies proton and/or pion conditions
          if (cuts->IsProton(part)) {
            // prot++;
            double dp_prot = pow(mc_event->prot_momX_mc_gen() - data->px(part), 2) +
                             pow(mc_event->prot_momY_mc_gen() - data->py(part), 2) +
                             pow(mc_event->prot_momZ_mc_gen() - data->pz(part), 2);
            proton_dps.push_back(std::make_pair(part, dp_prot));  // Store index and dp value for proton

            // event->SetProton(part);                 // for overlapped proton index
          }

          if (cuts->IsPip(part)) {
            // pip++;
            double dp_pip = pow(mc_event->pip_momX_mc_gen() - data->px(part), 2) +
                            pow(mc_event->pip_momY_mc_gen() - data->py(part), 2) +
                            pow(mc_event->pip_momZ_mc_gen() - data->pz(part), 2);
            pip_dps.push_back(std::make_pair(part, dp_pip));  // Store index and dp value for proton

            // event->SetPip(part);                // for overlapped pip index
          }
        }

        else {
          if (cuts->IsPim(part))

          {
            event->SetPim(part);

            // pim++;
          }
        }
      }
    }

    // // // Now, find the pair of proton and pip with the minimum dp_prot + dp_pip
    double min_dp_sum = std::numeric_limits<double>::max();
    int best_proton_index = -1;
    int best_pip_index = -1;
    std::vector<std::pair<int, int>> non_minimum_pairs;  // Stores all non-minimum proton-pip pairs

    // Loop over all combinations of protons and pions to find the minimum dp_prot + dp_pip
    for (size_t i = 0; i < proton_dps.size(); i++) {
      int prot_index = proton_dps[i].first;
      double dp_prot = proton_dps[i].second;

      for (size_t j = 0; j < pip_dps.size(); j++) {
        int pip_index = pip_dps[j].first;
        double dp_pip = pip_dps[j].second;

        double dp_sum = dp_prot + dp_pip;
        if (dp_sum < min_dp_sum) {
          min_dp_sum = dp_sum;
          best_proton_index = prot_index;
          best_pip_index = pip_index;
        }
      }
    }

    // // Set the proton and pip with the minimum dp_sum for further processing
    if (best_proton_index != -1 && best_pip_index != -1) {
      event->SetProton(best_proton_index);
      event->SetPip(best_pip_index);
    }

    // // Overlapped loop over all combinations of protons and pions
    // for (size_t i = 0; i < proton_dps.size(); i++) {
    //   int prot_index = proton_dps[i].first;
    //   event->SetProton(prot_index);  // for overlapped proton index
    // }

    // for (size_t j = 0; j < pip_dps.size(); j++) {
    //   int pip_index = pip_dps[j].first;
    //   event->SetPip(pip_index);  // for overlapped pip index
    // }

    // if (event->W() > 1.35 && event->W() <= 2.15 && event->Q2() <= 9.0 && event->Q2() >= 1.95 && event->weight() >
    // 0.0)
    {
      // if (event->TwoPion_missingPim() || event->TwoPion_missingPip() || event->TwoPion_missingProt()||
      // event->TwoPion_exclusive()) {

      // if (event->TwoPion_exclusive()) {
      //   for (size_t i = 0; i < event->GetProtons().size(); ++i) {
      //     for (size_t j = 0; j < event->GetPips().size(); ++j) {
      //       for (size_t k = 0; k < event->GetPims().size(); ++k) {
      //         // event->CalcMissMassExcl(*event->GetProtons()[i], *event->GetPips()[j], *event->GetPims()[k]);

      //         //                                 two_pion_Excl_events++;
      //         //                                 _hists->Fill_WvsQ2(event);

      //         // // You should have a similar method for π⁻ if applicable
      //         // dt->dt_calc(event->GetPimIndices()[k]);
      //         // _hists->Fill_deltat_pim_after_cut(data, dt, event->GetPimIndices()[k], event);
      //         // _hists->FillHists_pim_pid_with_cuts(data, event, event->GetPimIndices()[k]);
      //         //       }
      //         //     }
      //         //   }
      //         // }

      if (event->TwoPion_missingPim()) {
        // if (event->TwoPion_missingPip()) {
        // if (event->TwoPion_missingProt()) {
        // if (event->TwoPion_exclusive()) {
        // // twoPion_excl++;
        // // if (event->Inclusive()) {
        // {
        // {

        {  // // Retrieve the number of protons and pions in the event
          size_t num_protons = event->GetProtons().size();
          size_t num_pips = event->GetPips().size();
          int num_combinations = 0;
          // std::cout << "    prot size   :  " << num_protons << std::endl;
          // std::cout << "    pip size   :  " << num_pips << std::endl;

          // First loop: count valid combinations
          for (size_t i = 0; i < num_protons; ++i) {
            for (size_t j = 0; j < num_pips; ++j) {
              if (event->GetProtonIndices()[i] != event->GetPipIndices()[j]) {
                num_combinations++;
              }
            }
          }

          for (size_t i = 0; i < num_protons; ++i) {
            for (size_t j = 0; j < num_pips; ++j) {
              // if (event->GetProtonIndices()[i] == event->GetPipIndices()[j]) no_prot_pip++;

              // Exclude the case where the same particle is assigned as both proton and pip
              if (event->GetProtonIndices()[i] != event->GetPipIndices()[j]) {
                // // // std::cout << "    num_combinations   :  " << num_combinations << std::endl;
                // // // if (num_combinations >= 2)
                // if ((best_proton_index != event->GetProtonIndices()[i]) ||
                //     (best_pip_index != event->GetPipIndices()[j]))
                {
                  ////////////////////////////////////
                  int proton_part_idx = event->GetProtonIndices()[i];
                  int pip_part_idx = event->GetPipIndices()[j];

                  double dp_Sum = proton_dps[i].second + pip_dps[j].second;

                  {
                    event->SetSwappedProton(pip_part_idx);
                    event->SetSwappedPip(proton_part_idx);
                    // Extract velocity components with energy normalization
                    double v_original_x_Prot = event->GetProtons()[i]->Px() / event->GetProtons()[i]->E();
                    double v_original_y_Prot = event->GetProtons()[i]->Py() / event->GetProtons()[i]->E();
                    double v_original_z_Prot = event->GetProtons()[i]->Pz() / event->GetProtons()[i]->E();
                    double v_swapped_x_Prot = event->GetProtonsSwapped()->Px() / event->GetProtonsSwapped()->E();
                    double v_swapped_y_Prot = event->GetProtonsSwapped()->Py() / event->GetProtonsSwapped()->E();
                    double v_swapped_z_Prot = event->GetProtonsSwapped()->Pz() / event->GetProtonsSwapped()->E();

                    // Calculate delta_V^2
                    double dv2_Prot = pow(v_swapped_x_Prot - v_original_x_Prot, 2) +
                                      pow(v_swapped_y_Prot - v_original_y_Prot, 2) +
                                      pow(v_swapped_z_Prot - v_original_z_Prot, 2);
                    // std::cout << "Delta V^2: " << dv2 << std::endl;

                    // // Extract velocity components with energy normalization
                    // double v_original_x_Pip = event->GetPips()[j]->Px() / event->GetPips()[j]->E();
                    // double v_original_y_Pip = event->GetPips()[j]->Py() / event->GetPips()[j]->E();
                    // double v_original_z_Pip = event->GetPips()[j]->Pz() / event->GetPips()[j]->E();
                    // double v_swapped_x_Pip = event->GetPipsSwapped()->Px() / event->GetPipsSwapped()->E();
                    // double v_swapped_y_Pip = event->GetPipsSwapped()->Py() / event->GetPipsSwapped()->E();
                    // double v_swapped_z_Pip = event->GetPipsSwapped()->Pz() / event->GetPipsSwapped()->E();

                    // // Calculate delta_V^2
                    // double dv2_Pip = pow(v_swapped_x_Pip - v_original_x_Pip, 2) +
                    //                  pow(v_swapped_y_Pip - v_original_y_Pip, 2) +
                    //                  pow(v_swapped_z_Pip - v_original_z_Pip, 2);

                    event->CalcMissMassPim(*event->GetProtons()[i], *event->GetPips()[j]);
                    event->boost(*event->GetProtons()[i], *event->GetPips()[j]);
                    // event->CalcMissMassExcl(*event->GetProtons()[i], *event->GetPips()[j], *event->GetPims()[k]);

                    event->CalcMissMassPimSwapped();
                    ////////////  CONTROL OVER HAOW MANY FILLING PER EVENT /////////
                    ////////////  CONTROL OVER HAOW MANY FILLING PER EVENT /////////
                    ////////////  CONTROL OVER HAOW MANY FILLING PER EVENT /////////
                    // if (num_combinations == 2)
                    // if (event->MM2_mPim() < -0.1)
                    // if (dv2_Prot >= 0.0005 && dv2_Prot < 0.001)
                    // if (event->MM2_mPim() > -0.1 && event->MM2_mPim() < 0.1)
                    // if (num_combinations > 1) {

                    //   csv_data output;

                    //   // // // for generated case
                    //   output.w_mc = mc_event->W_mc();
                    //   output.q2_mc = mc_event->Q2_mc();
                    //   output.mm2_mPim_mc = mc_event->MM2_mPim_MC();
                    //   output.weight_mc = mc_event->weight();
                    //   _sync->write(output);
                    // }
                    csv_data output;

                    // // // // // //// using exclusive topology ...................................

                    // // // // // // output.electron_sector = event->sec();
                    // output.pim_sec = event->pimSec();
                    // output.pip_sec = event->pipSec();
                    // output.prot_sec = event->protSec();

                    // // output.prot_pid_mc = Prot_pid_mc;
                    // output.prot_pid_rec = Prot_pid_rec;
                    // // output.pip_pid_mc = Pip_pid_mc;
                    // output.pip_pid_rec = Pip_pid_rec;
                    // output.pim_pid_rec = Pim_pid_rec;

                    // output.w_before = event->W_before();
                    // output.q2_before = event->Q2_before();
                    output.w = event->W();
                    output.q2 = event->Q2();
                    output.dv2_prot = dv2_Prot;
                    output.dp_sum = dp_Sum;

                    // // // output.w_had = event->w_hadron();
                    // // // // output.w_diff = event->w_difference();
                    // // // output.w_had_corr = event->w_hadron_corr();
                    // // // // output.w_diff_corr = event->w_difference_corr();

                    // // // output.elec_mom = event->elec_mom();
                    // // // output.elec_energy = event->elec_En();
                    // // // output.elec_theta = event->Theta_Elec();
                    // // // output.corr_elec_mom = event->Corr_elec_mom();

                    // // // // //   // // for generated case
                    // // // output.w_mc = mc_event->W_mc();
                    // // // output.q2_mc = mc_event->Q2_mc();

                    // output.elec_mom_mc = mc_event->elec_mom_mc();
                    // // output.elec_energy_mc = mc_event->elec_En_mc();
                    // output.elec_theta_mc = mc_event->Theta_Elec_mc();

                    // output.elec_mom_rec = event->elec_mom();
                    // // output.elec_energy_rec = event->elec_En();
                    // output.elec_theta_rec = event->Theta_Elec();

                    // output.scalar_product = event->scalar_triple_product();

                    // //   // output.weight_exclusive = mc_event->weight();

                    // // // // for energy loss corrections : gen
                    // output.gen_prot_mom = (mc_event->prot_mom_mc_gen());
                    // output.gen_prot_theta = (mc_event->prot_theta_mc_gen());
                    // output.gen_prot_phi = (mc_event->prot_phi_mc_gen());

                    // output.gen_pip_mom = (mc_event->pip_mom_mc_gen());
                    // output.gen_pip_theta = (mc_event->pip_theta_mc_gen());
                    // output.gen_pip_phi = (mc_event->pip_phi_mc_gen());

                    // output.gen_pim_mom = (mc_event->pim_mom_mc_gen());
                    // output.gen_pim_theta = (mc_event->pim_theta_mc_gen());
                    // output.gen_pim_phi = (mc_event->pim_phi_mc_gen());

                    // // // // // // // // // // // missing
                    // // auto proton_vector = event->GetProtonIndices()[i];
                    // output.prot_mom_mProt = event->prot_momentum(*event->GetProtons()[i]);
                    // output.prot_theta_mProt = event->prot_theta_lab(*event->GetProtons()[i]);
                    // output.prot_phi_mProt = event->prot_Phi_lab(*event->GetProtons()[i]);

                    // output.pip_mom_mPip = event->pip_momentum(*event->GetPips()[j]);
                    // output.pip_theta_mPip = event->pip_theta_lab(*event->GetPips()[j]);
                    // output.pip_phi_mPip = event->pip_Phi_lab(*event->GetPips()[j]);

                    // output.pim_mom_mPim = event->pim_momentum();
                    // output.pim_theta_mPim = event->pim_theta_lab();
                    // output.pim_phi_mPim = event->pim_Phi_lab();

                    // // output.pim_mom_mPim_cm = event->pim_momentum_cm();
                    // // output.pim_theta_mPim_cm = event->pim_theta_cm();
                    // // output.pim_phi_mPim_cm = event->pim_Phi_cm();

                    // // // // // recon mes

                    output.prot_mom_exclusive = event->prot_momentum(*event->GetProtons()[i]);
                    output.prot_theta_exclusive = event->prot_theta_lab(*event->GetProtons()[i]);
                    output.prot_phi_exclusive = event->prot_Phi_lab(*event->GetProtons()[i]);
                    // // output.prot_dcr1theta_exclusive = event->thetaDCr1Prot();

                    output.pip_mom_exclusive = event->pip_momentum(*event->GetPips()[j]);
                    output.pip_theta_exclusive = event->pip_theta_lab(*event->GetPips()[j]);
                    output.pip_phi_exclusive = event->pip_Phi_lab(*event->GetPips()[j]);
                    // // output.pip_dcr1theta_exclusive = event->thetaDCr1Pip();

                    // output.pim_mom_exclusive = event->pim_momentum_measured();
                    // output.pim_theta_exclusive = event->pim_theta_lab_measured();
                    // output.pim_phi_exclusive = event->pim_Phi_lab_measured();
                    // // output.pim_dcr1theta_exclusive = event->thetaDCr1Pim();

                    // // // output.pim_mom_corr = event->pim_momentum_corrected();
                    // // // // output.pim_theta_corr = event->pim_theta_corrected();
                    // // // // output.pim_phi_corr = event->pim_Phi_corrected();

                    // output.mm2_mProt = event->MM2_mProt();
                    // // output.mm2_mProt_corr = event->MM2_mProt_corr();
                    // output.mm2_mPip = event->MM2_mPip();
                    // // // // output.mm2_mPip_corr = event->MM2_mPip_corr();
                    output.mm2_mPim = event->MM2_mPim();
                    // output.mm2_mPim_mc = mc_event->MM2_mPim_MC();
                    output.mm2_mPim_swapped = event->MM2_mPim_swapped();

                    // // // // output.mm2_mPim_corr = event->MM2_mPim_corr();

                    // output.mm2_exclusive_at_zero = event->MM2_exclusive();
                    // output.energy_x_mu = event->Energy_excl();

                    // // output.status_Pim = event->pimStatus();
                    output.status_Pip = event->pipStatus();
                    output.status_Prot = event->protStatus();

                    output.beta_Pip = event->betaPip();
                    output.beta_Prot = event->betaProt();

                    output.inv_ppip = event->inv_Ppip();
                    output.inv_ppim = event->inv_Ppim();
                    output.inv_pip_pim = event->inv_pip_pim();

                    output.alpha_Prot = event->alpha_pippim_pipf();
                    output.alpha_Pip = event->alpha_ppim_pipip();
                    output.alpha_Pim = event->alpha_ppip_pipim();
                    output.weight_exclusive = event->weight();
                    _sync->write(output);

                    /// ..........................................

                    // _sync->write(output);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
  // // Return the total number of events
  std::cout << " number of events = " << total << " numElec " << numElec << "  ratio  " << numElec / float(total) * 100
            << std::endl;
  std::cout << " number of mc elec = " << numElec_mc << std::endl;

  //  << "  mc  prot = " << numProt_mc
  //       << "  mc pip =
  // // " << numPip_mc
  // //           << "  mc pim  = " << numPim_mc << std::endl;

  // std::cout << " number of elec = " << numElec << "   prot = " << numProt << "  pip = " << numPip
  //           << "  pim  = " << numPim << std::endl;

  return num_of_events;
}
#endif
