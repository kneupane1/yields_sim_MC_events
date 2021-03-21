/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/

#ifndef CSV_H_GUARD
#define CSV_H_GUARD

#include <iostream>
#include "TFile.h"
#include "branches.hpp"
#include "colors.hpp"
#include "cuts.hpp"
#include "deltat.hpp"
#include "reaction.hpp"

template <class CutType>
std::string run(std::shared_ptr<TChain> _chain, int thread_id, bool mc) {
        std::string csv_out;
        // Get the number of events in this thread
        size_t num_of_events = (int)_chain->GetEntries();

        float beam_energy = rga_E0;
        if (std::is_same<CutType, rga_Cuts>::value) {
                beam_energy = rga_E0;
        } else if (std::is_same<CutType, rgf_Cuts>::value) {
                beam_energy = rgf_E0;
        }

        // Print some information for each thread
        std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
                  << num_of_events << " Events " << DEF << "===============\n";

        // Make a data object which all the branches can be accessed from
        auto data = std::make_shared<Branches12>(_chain, mc);

        // Total number of events "Processed"
        size_t total = 0;
        // For each event
        for (size_t current_event = 0; current_event < num_of_events; current_event++) {
                // Get current event
                data->GetEntry(current_event);

                // If we are the 0th thread print the progress of the thread every 1000 events
                if (thread_id == 0 && current_event % 1000 == 0)
                        std::cout << BLUE << "\t" << (100 * current_event / num_of_events) << " %\r" << DEF << std::flush;

                total++;
                auto cuts = std::make_shared<CutType>(data);
                if (!cuts->ElectronCuts()) continue;
                auto event = mc ? std::make_shared<MCReaction>(data, beam_energy) : std::make_shared<Reaction>(data, beam_energy);
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
                // Check the reaction class what kind of even it is and fill the appropriate histograms

                // if (event->NeutronPip()) csv_out += event->ReacToCsv();
        }

        // Return the total number of events
        return csv_out;
}
#endif
