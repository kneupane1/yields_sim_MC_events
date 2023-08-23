
#ifndef MAIN_H_GUARD
#define MAIN_H_GUARD

#include <iostream>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "branches.hpp"
#include "colors.hpp"
#include "cuts.hpp"
#include "histogram.hpp"
#include "reaction.hpp"
#include "syncfile.hpp"
using namespace std;

float alpha_CD[3][3] = {{0.9, 0.9, 0.95}, {0.8, 0.4, 0.8}, {0.5, 1.0, 0.5}};
float alpha_FD[3][4] = {{0.5, 0.6, 0.5, 0.5}, {0.1, 0.15, 0.5, 0.5}, {0.5, 0.15, 0.3, 0.3}};

// // float alpha_CD[3][3];
// // float alpha_FD[3][4];

// void initialize_alphas() {
//   std::srand(std::chrono::duration_cast<std::chrono::milliseconds>(
//                  std::chrono::high_resolution_clock::now().time_since_epoch())
//                  .count());

//   for (int i = 0; i < 3; i++) {
//     for (int j = 0; j < 4; j++) {
//       float rand_no = static_cast<float>(std::rand()) / RAND_MAX;
//       alpha_FD[i][j] += alpha_FD[i][j] * (rand_no - 0.5) * 0.1;
//       // alpha_FD[i][j] = rand_no;

//       std::cout << "fd rand no: " << rand_no << "  alpha : " << alpha_FD[i][j] << std::endl;
//     }
//   }

//   for (int i = 0; i < 3; i++) {
//     for (int j = 0; j < 3; j++) {
//       float rand_no = static_cast<float>(std::rand()) / RAND_MAX;
//       alpha_CD[i][j] += alpha_CD[i][j] * (rand_no - 0.5) * 0.1;
//       // alpha_CD[i][j] = rand_no;

//       std::cout << " cd rand no: " << rand_no << "  alpha : " << alpha_CD[i][j] << std::endl;
//     }
//   }
// }

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
  // auto data = std::make_shared<Branches12>(_chain, true);
  // for exp data use it
  auto data = std::make_shared<Branches12>(_chain);

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
    int sectorPim = -1;
    int sectorPip = -1;
    int sectorProt = -1;

    // if (data->mc_npart() < 1) continue;

    // // If we pass electron cuts the event is processed
    total++;

    // // Make a reaction class from the data given
    // auto mc_event = std::make_shared<MCReaction>(data, beam_energy);

    // for (int part = 1; part < data->mc_npart(); part++) {
    //   // Check particle ID's and fill the reaction class

    //   if (data->mc_pid(part) == PIP) {
    //     mc_event->SetMCPip(part);
    //   } else if (data->mc_pid(part) == PROTON) {
    //     mc_event->SetMCProton(part);
    //   } else if (data->mc_pid(part) == PIM) {
    //     mc_event->SetMCPim(part);
    //     // } else {
    //     //   mc_event->SetMCOther(part);
    //   }
    // }

    auto dt = std::make_shared<Delta_T>(data);
    auto cuts = std::make_shared<uconn_Cuts>(data);
    // auto cuts = std::make_shared<rga_Cuts>(data);
    if (!cuts->ElectronCuts()) continue;

    // Make a reaction class from the data given
    auto event = std::make_shared<Reaction>(data, beam_energy);
    event->SetMomCorrElec();

    // // For each particle in the event
    for (int part = 1; part < data->gpart(); part++) {
      dt->dt_calc(part);

      // Check particle ID's and fill the reaction class
      if (cuts->IsProton(part)) {
        if (cuts->HadronsCuts(part)) {
          event->SetProton(part);
          statusProt = abs(data->status(part));
          sectorProt = data->dc_sec(part);
          // if (statusProt < 4000 && statusProt > 2000) sectorProt = data->dc_sec(part);
        }

      } else if (cuts->IsPip(part)) {
        if (cuts->HadronsCuts(part)) {
          event->SetPip(part);
          statusPip = abs(data->status(part));
          sectorPip = data->dc_sec(part);
          // if (statusPip<4000 && statusPip> 2000) sectorPip = data->dc_sec(part);
        }
      } else if (cuts->IsPim(part)) {
        if (cuts->HadronsCuts(part)) {
          event->SetPim(part);
          statusPim = abs(data->status(part));
          sectorPim = data->dc_sec(part);
          // if (statusPim < 4000 && statusPim > 2000) sectorPim = data->dc_sec(part);
        }
      } else {
        event->SetOther(part);
      }
    }

    // if (event->TwoPion_missingPim() || event->TwoPion_missingPip() || event->TwoPion_missingProt() ||
    //     event->TwoPion_exclusive()) {
    // if (event->TwoPion_missingPim()) {
    // if (event->TwoPion_missingPip()) {
    // if (event->TwoPion_missingProt()) {
    // if (event->TwoPion_exclusive()) {
    // if (event->W() > 1.25 && event->W() < 2.55 && event->Q2() > 1.5 && event->Q2() < 10.5) {  // &&
    // // abs(event->Energy_excl()) < 0.3) {
    // float deltapCom = NAN;
    // float min_deltapCom = 9999.9;
    // float minimum_alphap = NAN;
    // float minimum_alphapip = NAN;
    // float minimum_alphapim = NAN;
    // // float alpha_universe[30] = {-5.0, -4.0, -3.0, -2.0, -1.0, -0.75, -0.5, -0.25, -0.1, 0, 0.05, 0.10, 0.15,
    // // 0.25, 0.4, 0.5, 0.65, 0.75, 0.85, 0.95, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.5, 3.0, 4.0, 5.0};

    // // float alpha_universe[11] = {-10.0, -7.0, -4.0, -1.0, -0.5,  0,  0.5, 1.0, 4.0, 7.0, 10};
    // float alpha_universe[20] = {-15.0, -10.0, -5.0, -2.0, -1.0, -0.75, -0.5, -0.25, -0.15, 0.0,
    //                             0.15,  0.25,  0.5,   0.75,  1.0,  2.0,  5.0, 7.0,  10.0,  15.0};

    // // float alpha_universe[20] = {-1.0, -0.9, -0.75, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0,
    // //                             0.1,  0.2,  0.3,   0.4,  0.5,  0.6,  0.7,  0.8,  0.9,  1.0};

    // float alpha_proton = NAN;
    // float deltapP = NAN;
    // float deltapInitialP = NAN;
    // float alpha_pip = NAN;
    // float deltapPip = NAN;
    // float deltapInitialPip = NAN;
    // float alpha_pim = NAN;
    // float deltapPim = NAN;
    // float deltapInitialPim = NAN;

    // for (int alpha_countP = 0; alpha_countP < 20; alpha_countP++) {
    //   for (int alpha_countPip = 0; alpha_countPip < 20; alpha_countPip++) {
    //     for (int alpha_countPim = 0; alpha_countPim < 20; alpha_countPim++) {
    //       alpha_proton = alpha_universe[alpha_countP];
    //       alpha_pip = alpha_universe[alpha_countPip];
    //       alpha_pim = alpha_universe[alpha_countPim];

    //       event->Prot_HMom_corr(statusProt, statusPip, statusPim, sectorProt, alpha_proton);
    //       event->Pip_HMom_corr(statusProt, statusPip, statusPim, sectorPip, alpha_pip);
    //       event->Pim_HMom_corr(statusProt, statusPip, statusPim, sectorPim, alpha_pim);

    //       deltapCom = pow((event->prot_momentum_corrected() - event->prot_momentum()), 2) +
    //                   pow((event->pip_momentum_corrected() - event->pip_momentum()), 2) +
    //                   pow((event->pim_momentum_corrected() - event->pim_momentum()), 2);

    //       if (deltapCom < min_deltapCom) {
    //         minimum_alphap = alpha_proton;
    //         minimum_alphapip = alpha_pip;
    //         minimum_alphapim = alpha_pim;
    //         min_deltapCom = deltapCom;
    //       } else
    //         continue;

    //       // if (deltapInitialP < deltapP)
    //       //   deltapInitialP = deltapP;
    //       // else
    //       //   continue;
    //       // if (deltapInitialPip < deltapPip)
    //       //   deltapInitialPip = deltapPip;
    //       // else
    //       //   continue;
    //       // if (deltapInitialPim < deltapPim)
    //       //   deltapInitialPim = deltapPim;
    //       // else
    //       //   continue;
    //     }
    //   }
    // }
    // event->Prot_HMom_corr(statusProt, statusPip, statusPim, sectorProt, minimum_alphap);
    // event->Pip_HMom_corr(statusProt, statusPip, statusPim, sectorPip, minimum_alphapip);
    // event->Pim_HMom_corr(statusProt, statusPip, statusPim, sectorPim, minimum_alphapim);

    // if (event->TwoPion_missingPim() || event->TwoPion_missingPip() || event->TwoPion_missingProt() ||
    // if (event->TwoPion_exclusive()) {
    if (event->TwoPion_missingPim()) {
      // if (event->TwoPion_missingPip()) {
    // if (event->TwoPion_missingProt()) {
      if (event->W() > 1.25 && event->W() < 2.55 && event->Q2() > 1.5 && event->Q2() < 10.5) {
        event->Prot_HMom_corr(statusProt, statusPip, statusPim, sectorProt, alpha_FD[0], alpha_CD[0]);
        event->Pip_HMom_corr(statusProt, statusPip, statusPim, sectorPip, alpha_FD[1], alpha_CD[1]);
        event->Pim_HMom_corr(statusProt, statusPip, statusPim, sectorPim, alpha_FD[2], alpha_CD[2]);

        //   // total++;
        csv_data output;

        // // // mPim .......................................
        output.electron_sector = event->sec();
        output.w = event->W();
        output.q2 = event->Q2();

        output.pim_mom_mPim = event->pim_momentum();
        output.pim_theta_mPim = event->pim_theta_lab();
        output.pim_phi_mPim = event->pim_Phi_lab();
        // output.pim_mom_mPim_cm = event->pim_momentum_cm();
        // output.pim_theta_mPim_cm = event->pim_theta_cm();
        // output.pim_phi_mPim_cm = event->pim_Phi_cm();
        output.mm2_mPim = event->MM2_mPim();
        output.mm2_mPim_corr = event->MM2_mPim_corr();
        // output.status_Pim = statusPim;
        // output.status_Pip = statusPip;
        // output.status_Prot = statusProt;
        output.weight_mPim = event->weight();

        // // mpip .......................................

        // output.pip_mom_mPip = event->pip_momentum();
        // output.pip_theta_mPip = event->pip_theta_lab();
        // output.pip_phi_mPip = event->pip_Phi_lab();
        // output.mm2_mPip = event->MM2_mPip();
        // output.mm2_mPip_corr = event->MM2_mPip_corr();
        // output.weight_mPip = event->weight();

        // // mProt .......................................

        // output.prot_mom_mProt = event->prot_momentum();
        // output.prot_theta_mProt = event->prot_theta_lab();
        // output.prot_phi_mProt = event->prot_Phi_lab();
        // output.mm2_mProt = event->MM2_mProt();
        // output.mm2_mProt_corr = event->MM2_mProt_corr();
        // output.weight_mProt = event->weight();

        _sync->write(output);
      }
    }
    // }

    // std::cout << "event no = " <<total << std::endl;
  }

  std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
  // Return the total number of events
  return num_of_events;
}
#endif
