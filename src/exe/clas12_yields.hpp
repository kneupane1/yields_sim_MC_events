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
        float beam_energy = 10.6041;
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
                if (thread_id == 0 && current_event % 1000 == 0)
                        std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;
                std::cout << "mc_npart " << data->mc_npart()<<'\n';

                // if (data->mc_npart() > 1) { //continue;

                // If we pass electron cuts the event is processed
                total++;

                // Make a reaction class from the data given
                auto mc_event = std::make_shared<MCReaction>(data, beam_energy);

                std::cout << "mc_npart after " << data->mc_npart()<<'\n';
                // if (data->mc_npart() > 1) {         //continue;

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
                auto cuts = std::make_shared<Cuts>(data, dt);
                if (!cuts->ElectronCuts()) continue;

                // Make a reaction class from the data given
                auto event = std::make_shared<Reaction>(data, beam_energy);
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

                event->boost();

                if (event->SinglePip()) {
                        total++;
                        csv_data output;
                        output.electron_sector = event->sec();
                        output.w = event->W();
                        output.q2 = event->Q2();
                        output.pip_theta = event->Theta_star();
                        output.pip_phi = event->Phi_star();
                        // output.mm2 = event->MM2();
                        output.mm2 = data->mc_weight();
                        if(data->mc_weight() > 0.5)
                                std::cout << "weight: " << data->mc_weight() <<'\n';
                        _sync->write(output);
                }

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
