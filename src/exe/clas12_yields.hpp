/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

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

size_t run(const std::shared_ptr<TChain>& _chain, const std::shared_ptr<SyncFile>& _sync, int thread_id) {
        // Get the number of events in this thread
        size_t num_of_events = (int)_chain->GetEntries();
        float beam_energy = NAN;
        if (getenv("BEAM_E") != NULL) beam_energy = atof(getenv("BEAM_E"));

        // Print some information for each thread
        std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
                  << num_of_events << " Events " << DEF << "===============\n";

        // Make a data object which all the branches can be accessed from
        auto data = std::make_shared<Branches12>(_chain);

        // Total number of events "Processed"
        size_t total = 0;
        // For each event
        for (size_t current_event = 0; current_event < num_of_events; current_event++) {
                // Get current event
                _chain->GetEntry(current_event);

                // If we are the 0th thread print the progress of the thread every 1000 events
                if (thread_id == 0 && current_event % 10000 == 0)
                        std::cerr << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;


                if (data->mc_npart() < 1) continue;
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
// changed new
                if (data->gpart() == 0) continue;
                bool elec = true;
                elec &= (data->charge(0) == NEGATIVE);
                elec &= (data->pid(0) == 11);
                if (!elec) continue;
                auto cuts = std::make_shared<Cuts>(data);
                if (!cuts->ElectronCuts()) continue;
// changed new
                // auto dt = std::make_shared<Delta_T>(data);
                // auto cuts = std::make_shared<Cuts>(data, dt);
                // if (!cuts->ElectronCuts()) continue;

                // Make a reaction class from the data given
                auto event = std::make_shared<Reaction>(data, beam_energy);
                auto dt = std::make_shared<Delta_T>(data);

                // For each particle in the event
                for (int part = 1; part < data->gpart(); part++) {
                        dt->dt_calc(part);

                        // Check particle ID's and fill the reaction class
                        if (cuts->IsPip(part)) {
                                event->SetPip(part);
                        } else if (cuts->IsProton(part)) {
                                event->SetProton(part);
                        } else if (cuts->IsPim(part)) {
                                event->SetPim(part);
                        } else {
                                event->SetOther(part);
                        }
                }

                //event->boost();

                if (event->TwoPion_missingPim()) {
                        // if(event->weight() > 0.5)
                        //         std::cout << "mc_weight from clas12_mc " << event->weight()<< '\n';
                        total++;
                        csv_data output;
                        output.w = event->W();
                        output.q2 = event->Q2();
                        output.pim_mom_mPim = event->pim_momentum();
                        output.pim_theta_mPim = event->pim_theta_lab();
                        output.pim_phi_mPim = event->pim_Phi_lab();
                        output.mm2_mPim = event->MM2();
                        output.weight_mPim = event->weight();
                        //
                        //
                        // output.electron_sector = event->sec();

                        // output.pip_theta = event->();
                        // output.pip_phi = event->Phi_star();
                        // output.mm2 = event->MM2();
// rsync
                        _sync->write(output);
                }
                // if (event->TwoPion_exclusive()) {
                //         total++;
                //         //std::cout << "WEIGHT " << event->weight()<<'\n';
                //
                //         csv_data output_exclusive;
                //         output_exclusive.w = event->W();
                //         output_exclusive.q2 = event->Q2();
                //         output_exclusive.pim_mom_exclusive = event->pim_momentum();
                //         output_exclusive.pim_theta_exclusive = event->pim_theta_lab();
                //         output_exclusive.pim_phi_exclusive = event->pim_Phi_lab();
                //         output_exclusive.mm2_exclusive = event->MM2();
                //         output_exclusive.mm2_exclusive_at_zero = event->MM2_exclusive();
                //         output_exclusive.weight_exclusive = event->weight();
                //
                //         _sync->write(output_exclusive);
                // }
        }
        std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
        // Return the total number of events
        return num_of_events;
}
#endif

/*
   ep -> e x+

   W for all events
   W for 2 particles
   W for 2 Part 2nd positive
   hist Phi_e - Phi_pos ~ 90
   W for cut around 90
   pos_mom vs pos_theta
   theta_p_pos_calc_from_electron - theta_pos_measured

   calc theta from magnitude of pos momentum

 */
