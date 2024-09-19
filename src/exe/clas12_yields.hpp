
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
  // For each event
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    // for (size_t current_event = 0; current_event < 100; current_event++) {
    // Get current event
    _chain->GetEntry(current_event);

    // If we are the 0th thread print the progress of the thread every 1000 events
    if (thread_id == 0 && current_event % 1000 == 0)
      std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;

    // /////////////////////////////// Generated sim only //////////////////////////////////////
    // if (_mc) {
    if (data->mc_npart() < 1) continue;
    // numElec_mc++;

    // If we pass electron cuts the event is processed
    total++;
    // std::cout << " event " << current_event << std::endl;

    // Make a reaction class from the data given
    auto mc_event = std::make_shared<MCReaction>(data, beam_energy);
    // if (mc_event->weight() <= 0.0) continue;
    // std::cout << " event " << current_event;

    for (int part = 1; part < data->mc_npart(); part++) {
      // Check particle ID's and fill the reaction class

      if (data->mc_pid(part) == PIP) {
        numPip_mc++;

        mc_event->SetMCPip(part);
      }
      if (data->mc_pid(part) == PROTON) {
        numProt_mc++;

        mc_event->SetMCProton(part);
        // std::cout << mc_event->GetMcProtons().size() << std::endl;

      } else if (data->mc_pid(part) == PIM) {
        numPim_mc++;

        mc_event->SetMCPim(part);
        // } else {
        //   mc_event->SetMCOther(part);
      }
    }
    // // std::cout << mc_event->GetMcProtons().size() << std::endl;
    // // // Retrieve the number of protons and pions in the event
    // size_t num_protons_mc = mc_event->GetMcProtons().size();
    // size_t num_pips_mc = mc_event->GetMcPips().size();
    // for (size_t i = 0; i < num_protons_mc; ++i) {
    //   for (size_t j = 0; j < num_pips_mc; ++j) {
    //     // std::cout << mc_event->GetMcProtons().size() << std::endl;
    //     mc_event->CalcMissMassPimMC(*mc_event->GetMcProtons()[i], *mc_event->GetMcPips()[j]);
    //     // std::cout << "  mc mm2 mPim inside loop " << mc_event->MM2_mPim_MC();
    //   }
    // }
    // mc_event->CalcMissMassPimMC();

    // std::cout << "  mc mass " << mc_event->MM_mPim_MC() << std::endl;
    // std::cout << "  mc mm2 mPim " << mc_event->MM2_mPim_MC() << std::endl;

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

    // // For each particle in the event
    for (int part = 1; part < data->gpart(); part++) {
      dt->dt_calc(part);

      // Check particle ID's and fill the reaction class

      if (cuts->IsPip(part)) {
        // Get the generated pip (π⁺) indices
        // const std::vector<int>& mc_pip_indices = mc_event->GetPipMcIndices();
        // if (mc_pip_indices.size() != 1) std::cout << "Number of gen pip : " << mc_pip_indices.size() << std::endl;
        // int mc_pip = mc_pip_indices[0];  // Access the first (and only) pip index
        // Pip_pid_mc = data->mc_pid(mc_pip);
        Pip_pid_rec = data->pid(part);
        // for (int mc_pip : mc_pip_indices) {
        // if (data->mc_pid(mc_pip) != data->pid(part)) {
        numPip++;
        event->SetPip(part);
        // break;
        //   }
        // }
      }

      if (cuts->IsProton(part)) {
        // Get the generated proton indices
        const std::vector<int>& mc_proton_indices = mc_event->GetProtonMcIndices();
        int mc_proton = mc_proton_indices[0];  // Access the first (and only) proton index
        Prot_pid_mc = data->mc_pid(mc_proton);
        Prot_pid_rec = data->pid(part);
        // for (int mc_proton : mc_proton_indices) {
        //   if (data->mc_pid(mc_proton) != data->pid(part)) {
        numProt++;
        event->SetProton(part);
        // break;  // Stop once a match is found
        // }
        // }
      }
      if (cuts->IsPim(part)) {
        // event->SetPim(part);
        // // Get the generated pim (π⁻) indices
        // const std::vector<int>& mc_pim_indices = mc_event->GetPimMcIndices();
        // int mc_pim = mc_pim_indices[0];  // Access the first (and only) pim index
        // Pim_pid_mc = data->mc_pid(mc_pim);
        Pim_pid_rec = data->pid(part);
        // for (int mc_pim : mc_pim_indices) {
        //   if (data->mc_pid(mc_pim) != data->pid(part)) {
        numPim++;
        event->SetPim(part);
        //     break;
        //   }
        // }
      } else {
        event->SetOther(part);
      }
    }
    // const std::vector<int>& mc_proton_indices = mc_event->GetProtonMcIndices();
    // for (int mc_proton : mc_proton_indices) {
    //   if (mc_proton != 1) std::cout << "    _prot_mc_indices  " << mc_proton << std::endl;
    // }
    // const std::vector<int>& mc_pip_indices = mc_event->GetPipMcIndices();
    // for (int mc_pip : mc_pip_indices) {
    //   // if (mc_pip != 1)
    //   std::cout << "    _pip_mc_indices  " << mc_pip << std::endl;
    // }

    // const std::vector<int>& mc_pim_indices = mc_event->GetPimMcIndices();
    // for (int mc_pim : mc_pim_indices) {
    //   // if (mc_pim != 1)
    //   std::cout << "    _pim_mc_indices  " << mc_pim << std::endl;
    //   std::cout << "    mc pid  " << data->mc_pid(mc_pim) << std::endl;
    //   std::cout << "    rec pid  " << data->pid(mc_pim) << std::endl;
    // }

    if (event->W() > 1.35 && event->W() <= 2.15 && event->Q2() <= 9.0 && event->Q2() >= 1.95 && event->weight() > 0.0) {
      // if (event->TwoPion_missingPim() || event->TwoPion_missingPip() || event->TwoPion_missingProt()||
      // event->TwoPion_exclusive()) {

      if (event->TwoPion_exclusive()) {
        for (size_t i = 0; i < event->GetProtons().size(); ++i) {
          for (size_t j = 0; j < event->GetPips().size(); ++j) {
            for (size_t k = 0; k < event->GetPims().size(); ++k) {
              // event->CalcMissMassExcl(*event->GetProtons()[i], *event->GetPips()[j], *event->GetPims()[k]);

              //                                 two_pion_Excl_events++;
              //                                 _hists->Fill_WvsQ2(event);

              // // You should have a similar method for π⁻ if applicable
              // dt->dt_calc(event->GetPimIndices()[k]);
              // _hists->Fill_deltat_pim_after_cut(data, dt, event->GetPimIndices()[k], event);
              // _hists->FillHists_pim_pid_with_cuts(data, event, event->GetPimIndices()[k]);
              //       }
              //     }
              //   }
              // }

              // if (event->TwoPion_missingPim()) {
              //   // if (event->TwoPion_missingPip()) {
              //   // if (event->TwoPion_missingProt()) {
              //   // if (event->TwoPion_exclusive()) {
              //   // // twoPion_excl++;
              //   // // if (event->Inclusive()) {
              //   // {
              //   // {

              //   {  // // Retrieve the number of protons and pions in the event
              //     size_t num_protons = event->GetProtons().size();
              //     size_t num_pips = event->GetPips().size();

              //     // // Calculate the total number of combinations
              //     // size_t num_combinations = num_protons * num_pips;

              //     // // Get the vector of protons
              //     // const auto &protons = event->GetProtons();
              //     // const auto &pip = event->GetPips();

              //     for (size_t i = 0; i < num_protons; ++i) {
              //       for (size_t j = 0; j < num_pips; ++j) {
              //         // if (event->GetProtonIndices()[i] == event->GetPipIndices()[j]) no_prot_pip++;

              //         // // // Print//////////////////////////////////////
              //         // std::cout << "Event " << current_event << ": "
              //         //           << num_protons << " proton(s), "
              //         //           << num_pips << " pip(s), "
              //         //           << num_combinations << " combination(s)." << std::endl;

              //         // if (both_prot_pip >= 1)
              //         // Access the i-th proton and j-th pip

              //         // Exclude the case where the same particle is assigned as both proton and pip
              if (event->GetProtonIndices()[i] != event->GetPipIndices()[j]) {
                // event->CalcMissMassPim(*event->GetProtons()[i], *event->GetPips()[j]);
                // event->boost(*event->GetProtons()[i], *event->GetPips()[j]);
                event->CalcMissMassExcl(*event->GetProtons()[i], *event->GetPips()[j], *event->GetPims()[k]);

                // // std::cout << "  rec mass mPim " << event->MM_mPim() << std::endl;
                // {
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

                // output.prot_pid_mc = Prot_pid_mc;
                output.prot_pid_rec = Prot_pid_rec;
                // output.pip_pid_mc = Pip_pid_mc;
                output.pip_pid_rec = Pip_pid_rec;
                output.pim_pid_rec = Pim_pid_rec;

                output.w_before = event->W_before();
                output.q2_before = event->Q2_before();
                output.w = event->W();
                output.q2 = event->Q2();
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

                // // // // // // // // // // missing
                // auto proton_vector = event->GetProtonIndices()[i];
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

                // output.prot_mom_exclusive = event->prot_momentum_measured();
                // output.prot_theta_exclusive = event->prot_theta_lab_measured();
                // output.prot_phi_exclusive = event->prot_Phi_lab_measured();
                // // output.prot_dcr1theta_exclusive = event->thetaDCr1Prot();

                // output.pip_mom_exclusive = event->pip_momentum_measured();
                // output.pip_theta_exclusive = event->pip_theta_lab_measured();
                // output.pip_phi_exclusive = event->pip_Phi_lab_measured();
                // // output.pip_dcr1theta_exclusive = event->thetaDCr1Pip();

                // output.pim_mom_exclusive = event->pim_momentum_measured();
                // output.pim_theta_exclusive = event->pim_theta_lab_measured();
                // output.pim_phi_exclusive = event->pim_Phi_lab_measured();
                // // output.pim_dcr1theta_exclusive = event->thetaDCr1Pim();

                // // // output.pim_mom_corr = event->pim_momentum_corrected();
                // // // // output.pim_theta_corr = event->pim_theta_corrected();
                // // // // output.pim_phi_corr = event->pim_Phi_corrected();

                output.mm2_mProt = event->MM2_mProt();
                // output.mm2_mProt_corr = event->MM2_mProt_corr();
                output.mm2_mPip = event->MM2_mPip();
                // // // output.mm2_mPip_corr = event->MM2_mPip_corr();
                output.mm2_mPim = event->MM2_mPim();
                // output.mm2_mPim_mc = mc_event->MM2_mPim_MC();

                // // // output.mm2_mPim_corr = event->MM2_mPim_corr();

                output.mm2_exclusive_at_zero = event->MM2_exclusive();
                output.energy_x_mu = event->Energy_excl();

                output.status_Pim = event->pimStatus();
                output.status_Pip = event->pipStatus();
                output.status_Prot = event->protStatus();

                // output.beta_Pip = event->betaPip();
                // output.beta_Prot = event->betaProt();

                // // output.inv_ppip = event->inv_Ppip();
                // // output.inv_ppim = event->inv_Ppim();
                // // output.inv_pip_pim = event->inv_Pippim();

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
  std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
  // // Return the total number of events
  // std::cout << " number of events = " << total << "   exclusive twoPion = " << twoPion_excl <<
  // std::endl;
  // // std::cout << " number of mc elec = " << numElec_mc << "  mc  prot = " << numProt_mc << "  mc pip =
  // " << numPip_mc
  // //           << "  mc pim  = " << numPim_mc << std::endl;

  // std::cout << " number of elec = " << numElec << "   prot = " << numProt << "  pip = " << numPip
  //           << "  pim  = " << numPim << std::endl;

  return num_of_events;
}
#endif
