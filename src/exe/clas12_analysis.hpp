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

size_t run(std::shared_ptr<TChain> _chain, std::shared_ptr<Histogram> _hists, int thread_id);
size_t run_files(std::vector<std::string> inputs, std::shared_ptr<Histogram> hists, int thread_id);

size_t run_files(std::vector<std::string> inputs, std::shared_ptr<Histogram> hists, int thread_id) {
        // Called once for each thread
        // Make a new chain to process for this thread
        auto chain = std::make_shared<TChain>("clas12");
        // Add every file to the chain
        for (auto in : inputs) chain->Add(in.c_str());

        // Run the function over each thread
        return run(chain, hists, thread_id);
}

size_t run(std::shared_ptr<TChain> _chain, std::shared_ptr<Histogram> _hists, int thread_id) {
        // Get the number of events in this thread
        size_t num_of_events = (int)_chain->GetEntries();
        float beam_energy = NAN;
        if (getenv("CLAS12_E") != NULL) beam_energy = atof(getenv("CLAS12_E"));

        // Print some information for each thread
        std::cout << "=============== " << RED << "Thread " << thread_id << DEF << " =============== " << BLUE
                  << num_of_events << " Events " << DEF << "===============\n";

        // Make a data object which all the branches can be accessed from
        auto data = std::make_shared<Branches12>(_chain);

        // Total number of events "Processed"
        size_t total = 0;
        int no_of_events = 0, prot = 0, pim = 0, pip = 0, other = 0, twopi = 0;
        int total_part = 0;
        // For each event
        for (size_t current_event = 0; current_event < num_of_events; current_event++) {
                // Get current event
                _chain->GetEntry(current_event);
                // If we are the 0th thread print the progress of the thread every 1000 events
                if (thread_id == 0 && current_event % 1000 == 0)
                        std::cout << "\t" << (100 * current_event / num_of_events) << " %\r" << std::flush;
                //  _hists->FillHists_electron_cuts(event, data->ec_tot_energy(0), data->p(0));
                auto event = std::make_shared<Reaction>(data, beam_energy);

                //  if (data->charge(0) == -1) _hists->FillHists_electron_cuts(data);

                auto cuts = std::make_shared<Cuts>(data);
                if (!cuts->ElectronCuts()) continue;
                no_of_events++;
                //  _hists->FillHists_electron_with_cuts(data);

                //_hists->Fill_EC(data->ec_tot_energy(0), data->p(0));

                // If we pass electron cuts the event is processed
                total++;
                //_hists->Fill_SF(data, part);
                int status_pim = -9999;
                int status_pip = -9999;
                int status_prot = -9999;



                // Make a reaction class from the data given
                auto dt = std::make_shared<Delta_T>(data);

                // For each particle in the event
                //  if (data->gpart() > 1) continue;
                for (int part = 1; part < data->gpart(); part++) {
                        total_part++;
                        dt->dt_calc(part);
                        //dt->dt_calc_ctof(part);

                        _hists->Fill_MomVsBeta(data, part);
                        //if (event->TwoPion_missingPim()) {
                        if (event->TwoPion_missingPim()) {
                                //if(data->status(part) >=2000 && data->status(part) < 6000) {
                                _hists->Fill_deltat_pi(data, dt, part);
                                //}
                                _hists->Fill_deltat_prot(data, dt, part);
                                //}
                        }
//_hists->Fill_deltat_positive(data, dt, event, part);

                        // Check particle ID's and fill the reaction class

                        if (cuts->IsProton(part)) {
                                prot++;
                                event->SetProton(part);
                                status_prot = abs(data->status(part));

                        } else if (cuts->IsPip(part)) {
                                pip++;
                                event->SetPip(part);
                                status_pip = abs(data->status(part));

                        } else if (cuts->IsPim(part)) {
                                //   if (event->MM_cut()) {
                                pim++;
                                event->SetPim(part);
                                status_pim = abs(data->status(part));
                        }

                        // } else if (cuts->IsmissingPim(part)) {
                        //   if (event->MM_cut()) event->SetmissingPim(part);

                        else {
                                other++;
                                event->SetOther(part);
                        }
                }
                // // //  if (event->TwoPion_missingPim()) {
                // // for (int part = 1; part < data->gpart(); part++) {
                // //   if (event->MM_cut()) event->SetmissingPim(part);
                // //   //  }
                // // }
                //
                // // Check the reaction class what kind of even it is and fill the appropriate histograms
                //
                // _hists->Fill_pi0(event);
                // //_hists->Fill_histthreeD(event);
                // _hists->Fill_WvsQ2(event);
                // _hists->Fill_theta_pim_measured(event);
                //
                // _hists->Fill_efficiency_check_CD_rec(event);
                // _hists->Fill_efficiency_check_FD_rec(event);
                //
                //
                // if(status_pim >=4000 && status_pim < 6000) {
                //         _hists->Fill_efficiency_check_CD(event);
                //
                // }
                // else if(status_pim >=2000 && status_pim < 4000) {
                //         _hists->Fill_efficiency_check_FD(event);
                // }
                //
                // if (event->TwoPion_missingProt())
                //         _hists->Fill_efficiency_check_1d_hist_prot(event);
                // if (event->TwoPion_missingPip())
                //         _hists->Fill_efficiency_check_1d_hist_pip(event);
                // if (event->TwoPion_missingPim())
                //         _hists->Fill_efficiency_check_1d_hist_pim(event);
                // // if(status_pip >=4000 && status_pip < 6000) {
                // //         _hists->Fill_efficiency_check_missingPip_CD(event);
                // //
                // // }
                // // else if(status_pip >=2000 && status_pip < 4000) {
                // //         _hists->Fill_efficiency_check_missingPip_FD(event);
                // // }
                // //
                // // if(status_prot >=4000 && status_prot < 6000) {
                // //         _hists->Fill_efficiency_check_missingProt_CD(event);
                // // }
                // // else if(status_prot >=2000 && status_prot < 4000) {
                // //         _hists->Fill_efficiency_check_missingProt_FD(event);
                // }
                // if (event->TwoPion_exclusive()) {
                //         _hists->Fill_eff_ckeck_exclusive(event);
                // }

                // if (event->TwoPion_missingPim()) {
                //         twopi++;
                //         _hists->Fill_eff_ckeck_mPim(event);
                //         //_hists->Fill_hists4D_background(event);
                //         _hists->Fill_WvsQ2_twoPi(event);
                //         // _hists->Fill_pi0(event);
                //         // _hists->Fill_x_mu(event);
                //         // _hists->Fill_histSevenD_pim(event);
                //         // _hists->Fill_histSevenD_pip(event);
                //         // _hists->Fill_histSevenD_prot(event);
                //
                // }
                //
                // // if (event->SinglePip()) {
                // //   _hists->Fill_WvsQ2_singlePi(event);
                // // }
                // // if (event->NeutronPip()) _hists->Fill_WvsQ2_Npip(event);
        }
        std::cout << "Percent = " << 100.0 * total / num_of_events << std::endl;
        std::cout << "no of events with good electron  " << no_of_events << '\n'
                  << "pim  " << pim << '\n'
                  << "prot  " << prot << '\n'
                  << "pip  " << pip << '\n'
                  << "other  " << other << " total_part " << total_part << '\n'
                  << "twopion = " << twopi << '\n';
        // Return the total number of events
        return num_of_events;
}
#endif
