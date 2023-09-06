
#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"
#include "cuts.hpp"
#include "histogram.hpp"
#include "reaction.hpp"
#include "syncfile.hpp"

template <class CutType>
size_t run(std::shared_ptr<TChain> _chain, const std::shared_ptr<SyncFile>& _sync, int thread_id) {
  // Get the number of events in this thread
  size_t num_of_events = (int)_chain->GetEntries();

  float beam_energy = 10.6;
  if (std::is_same<CutType, rga_Cuts>::value) {
    beam_energy = 10.6;
  } else if (std::is_same<CutType, uconn_Cuts>::value) {
    beam_energy = 10.6;
    // } else if (std::is_same<CutType, rgf_Cuts>::value) {
    //         beam_energy = rgf_E0;
    // }
    // else if (std::is_same<CutType, rgk_Cuts>::value) {
    //         beam_energy = rgk_E0;
  }

  if (getenv("BEAM_E") != NULL) beam_energy = atof(getenv("BEAM_E"));

  // Print some information for each thread
  std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
            << num_of_events << " Events " << DEF << "===============\n";

  // Make a data object which all the branches can be accessed from
  // for sim data use it
  auto data = std::make_shared<Branches12>(_chain, true);
  // for exp data use it
  // auto data = std::make_shared<Branches12>(_chain);

  // Total number of events "Processed"
  size_t total = 0;
  // For each event
  for (size_t current_event = 0; current_event < num_of_events; current_event++) {
    // for (size_t current_event = 0; current_event < 350; current_event++) {
    // Get current event
    _chain->GetEntry(current_event);

    // If we are the 0th thread print the progress of the thread every 1000 events
    if (thread_id == 0 && current_event % 1000 == 0)
      std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;

    int statusPim = -9999;
    int statusPip = -9999;
    int statusProt = -9999;

    if (data->mc_npart() < 1) continue;

    // If we pass electron cuts the event is processed
    total++;

    // Make a reaction class from the data given
    auto mc_event = std::make_shared<MCReaction>(data, beam_energy);

    for (int part = 1; part < data->mc_npart(); part++) {
      // Check particle ID's and fill the reaction class

      if (data->mc_pid(part) == PIP) {
        mc_event->SetMCPip(part);
      } else if (data->mc_pid(part) == PROTON) {
        mc_event->SetMCProton(part);
      } else if (data->mc_pid(part) == PIM) {
        mc_event->SetMCPim(part);
        // } else {
        //   mc_event->SetMCOther(part);
      }
    }

    auto dt = std::make_shared<Delta_T>(data);
    auto cuts = std::make_shared<uconn_Cuts>(data);
    // auto cuts = std::make_shared<rga_Cuts>(data);
    if (!cuts->ElectronCuts()) continue;

    // Make a reaction class from the data given
    auto event = std::make_shared<Reaction>(data, beam_energy);
    // event->SetMomCorrElec();

    // // For each particle in the event
    for (int part = 1; part < data->gpart(); part++) {
      dt->dt_calc(part);
      if (!isnan(data->sc_ctof_time(part))) std::cout << "sc ctof time part :  " << data->sc_ctof_time(part) << std::endl;

      // Check particle ID's and fill the reaction class
      if (cuts->IsProton(part)) {
        if (cuts->HadronsCuts(part)) {
          event->SetProton(part);
          statusProt = abs(data->status(part));
          // if (statusProt < 4000 && statusProt > 2000) sectorProt = data->dc_sec(part);
        }

      } else if (cuts->IsPip(part)) {
        if (cuts->HadronsCuts(part)) {
          event->SetPip(part);
          statusPip = abs(data->status(part));
          // if (statusPip<4000 && statusPip> 2000) sectorPip = data->dc_sec(part);
        }
      } else if (cuts->IsPim(part)) {
        if (cuts->HadronsCuts(part)) {
          event->SetPim(part);
          statusPim = abs(data->status(part));
          // if (statusPim < 4000 && statusPim > 2000) sectorPim = data->dc_sec(part);
        }
      } else {
        event->SetOther(part);
      }
    }
    // std::cout << "sc ctof time elec :  " << data->sc_ctof_time(0) << std::endl;

    // if (event->TwoPion_missingPim() || event->TwoPion_missingPip() || event->TwoPion_missingProt() ||
    //     event->TwoPion_exclusive()) {
    // if (event->TwoPion_missingPip()) {
      if (event->TwoPion_missingPim()) {
    // if (event->TwoPion_missingProt()) {
      if (event->W() > 1.25 && event->W() < 2.55 && event->Q2() > 1.5 && event->Q2() < 10.5) {


        //   // total++;
        csv_data output;
        output.electron_sector = event->sec();

        // //   // // for generated case
        // output.w_mc = mc_event->W_mc();
        // output.q2_mc = mc_event->Q2_mc();

        output.w = event->W();
        output.q2 = event->Q2();
        // // mPim .......................................

          output.pim_mom_mPim = event->pim_momentum();
          output.pim_theta_mPim = event->pim_theta_lab();
          output.pim_phi_mPim = event->pim_Phi_lab();

        // output.pim_mom_mPim_cm = event->pim_momentum_cm();
        // output.pim_theta_mPim_cm = event->pim_theta_cm();
        // output.pim_phi_mPim_cm = event->pim_Phi_cm();

        output.mm2_mPim = event->MM2_mPim();
        // // output.mm2_mPim_corr = event->MM2_mPim_corr();

        // // output.status_Pim = statusPim;
        // // output.status_Pip = statusPip;
        // // output.status_Prot = statusProt;
        // output.weight_mPim = event->weight();

        // // mPip .......................................

        // output.pip_mom_mPip = event->pip_momentum();
        // output.pip_theta_mPip = event->pip_theta_lab();
        // output.pip_phi_mPip = event->pip_Phi_lab();
        // output.mm2_mPip = event->MM2_mPip();
        // // output.mm2_mPip_corr = event->MM2_mPip_corr();
        // output.weight_mPip = event->weight();

        // mProt .......................................

        // output.prot_mom_mProt = event->prot_momentum();
        // output.prot_theta_mProt = event->prot_theta_lab();
        // output.prot_phi_mProt = event->prot_Phi_lab();
        // output.mm2_mProt = event->MM2_mProt();
        // // output.mm2_mProt_corr = event->MM2_mProt_corr();
        // output.weight_mProt = event->weight();

        output.inv_ppip = event->inv_Ppip();
        // output.inv_ppim = event->inv_Ppim();
        // output.inv_pip_pim = event->inv_Pippim();
        output.weight_mPim = mc_event->weight();

        _sync->write(output);
      }
    }
  }
  std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
  // Return the total number of events
  return num_of_events;
}
#endif
