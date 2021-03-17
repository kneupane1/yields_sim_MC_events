/**************************************/
/*																		*/
/*  Created by Nick Tyler             */
/*	University Of South Carolina      */
/**************************************/
#include "histogram.hpp"

Histogram::Histogram(const std::string& output_file) {
  RootOutputFile = std::make_shared<TFile>(output_file.c_str(), "RECREATE");
  def = std::make_shared<TCanvas>("def");
  if (getenv("BEAM_E") != NULL) {
    if (atof(getenv("BEAM_E")) < 3) {
      q2_max = 1.0;
      w_max = 3.5;
      p_max = 3.0;
    } else if (atof(getenv("BEAM_E")) < 8) {
      q2_max = 4.0;
      w_max = 4.0;
      p_max = 4.0;
    } else if (atof(getenv("BEAM_E")) < 9) {
      q2_max = 7.0;
      w_max = 7.0;
      p_max = 7.0;
    }
  }
  // Kinematics
  momentum = std::make_shared<TH1D>("mom", "mom", bins, p_min, p_max);
  W_hist = std::make_shared<TH1D>("W", "W", bins, zero, w_max);
  Q2_hist = std::make_shared<TH1D>("Q2", "Q2", bins, zero, q2_max);
  W_vs_q2 = std::make_shared<TH2D>("W_vs_q2", "W_vs_q2", bins, zero, w_max, bins, zero, q2_max);

  W_MC_hist = std::make_shared<TH1D>("W_MC", "W", bins, zero, w_max);
  Q2_MC_hist = std::make_shared<TH1D>("Q2_MC", "Q2", bins, zero, q2_max);
  W_vs_q2_MC = std::make_shared<TH2D>("W_vs_q2_MC", "W_vs_q2", bins, zero, w_max, bins, zero, q2_max);

  MM_neutron = std::make_shared<TH1D>("missMass", "missMass", bins, zero, 4.0);

  W_hist_singlePip = std::make_shared<TH1D>("W_singlePip", "W_singlePip", bins, zero, w_max);
  Q2_hist_singlePip = std::make_shared<TH1D>("Q2_singlePip", "Q2_singlePip", bins, zero, q2_max);
  W_vs_q2_singlePip =
      std::make_shared<TH2D>("W_vs_q2_singlePip", "W_vs_q2_singlePip", bins, zero, w_max, bins, zero, q2_max);

  EC_sampling_fraction =
      std::make_shared<TH2D>("EC_sampling_fraction", "EC_sampling_fraction", bins, p_min, p_max, bins, zero, 1.0);
  makeHists_sector();
  makeHists_deltat();
  makeHists_MomVsBeta();
}

Histogram::~Histogram() { this->Write(); }

void Histogram::Write() {
  std::cout << GREEN << "Writting" << DEF << std::endl;
  Write_EC();
  std::cerr << BOLDBLUE << "WvsQ2()" << DEF << std::endl;
  TDirectory* WvsQ2_folder = RootOutputFile->mkdir("W vs Q2");
  WvsQ2_folder->cd();
  Write_WvsQ2();

  std::cerr << BOLDBLUE << "Write_MomVsBeta()" << DEF << std::endl;
  TDirectory* Write_MomVsBeta_folder = RootOutputFile->mkdir("Mom Vs Beta");
  Write_MomVsBeta_folder->cd();
  Write_MomVsBeta();

  std::cerr << BOLDBLUE << "Write_deltat()" << DEF << std::endl;
  // TDirectory* Write_deltat_folder = RootOutputFile->mkdir("Delta_t");
  // Write_deltat_folder->cd();
  Write_deltat();

  Write_WvsQ2MC();

  std::cerr << BOLDBLUE << "Done Writing!!!" << DEF << std::endl;
}

void Histogram::Fill_WvsQ2(const std::shared_ptr<Reaction>& _e) {
  W_vs_q2->Fill(_e->W(), _e->Q2());
  W_hist->Fill(_e->W());
  Q2_hist->Fill(_e->Q2());

  short sec = _e->sec();
  if (sec > 0 && sec <= 6) {
    W_vs_q2_sec[sec - 1]->Fill(_e->W(), _e->Q2());
    W_sec[sec - 1]->Fill(_e->W());
  }
  short det = _e->det();
  if (det == 1 && _e->W() <= 3.5) {
    W_det[0]->Fill(_e->W());
    WQ2_det[0]->Fill(_e->W(), _e->Q2());
  } else if (det == 2) {
    W_det[1]->Fill(_e->W());
    WQ2_det[1]->Fill(_e->W(), _e->Q2());
  } else {
    W_det[2]->Fill(_e->W());
    WQ2_det[2]->Fill(_e->W(), _e->Q2());
  }
}
void Histogram::Fill_WvsQ2(const std::shared_ptr<MCReaction>& _e) {
  W_vs_q2_MC->Fill(_e->W(), _e->Q2(), _e->weight());
  W_MC_hist->Fill(_e->W(), _e->weight());
  Q2_MC_hist->Fill(_e->Q2(), _e->weight());

  short sec = _e->sec();
  if (sec > 0 && sec <= 6) {
    W_vs_q2_sec_MC[sec - 1]->Fill(_e->W(), _e->Q2(), _e->weight());
    W_sec_MC[sec - 1]->Fill(_e->W(), _e->weight());
  }

  short det = _e->det();
  if (det == 1 && _e->W() <= 3.5) {
    W_det_MC[0]->Fill(_e->W(), _e->weight());
    WQ2_det_MC[0]->Fill(_e->W(), _e->Q2(), _e->weight());
  } else if (det == 2) {
    W_det_MC[1]->Fill(_e->W(), _e->weight());
    WQ2_det_MC[1]->Fill(_e->W(), _e->Q2(), _e->weight());
  } else {
    W_det_MC[2]->Fill(_e->W(), _e->weight());
    WQ2_det_MC[2]->Fill(_e->W(), _e->Q2(), _e->weight());
  }
}

// W and Q^2
void Histogram::Fill_WvsQ2_singlePip(const std::shared_ptr<Reaction>& _e) {
  short sec = _e->sec();
  W_vs_q2_singlePip->Fill(_e->W(), _e->Q2());
  W_hist_singlePip->Fill(_e->W());
  Q2_hist_singlePip->Fill(_e->Q2());
  MM_neutron->Fill(_e->MM());
  if (sec > 0 && sec <= 6) {
    W_vs_MM_singlePip[sec - 1]->Fill(_e->W(), _e->MM());
    W_vs_q2_singlePip_sec[sec - 1]->Fill(_e->W(), _e->Q2());
    W_singlePip_sec[sec - 1]->Fill(_e->W());
    MM_neutron_sec[sec - 1]->Fill(_e->MM());
  }
}

// W and Q^2
void Histogram::Fill_WvsQ2_Npip(const std::shared_ptr<Reaction>& _e) {
  short sec = _e->sec();
  if (sec > 0 && sec <= 6) {
    W_vs_q2_Npip_sec[sec - 1]->Fill(_e->W(), _e->Q2());
    W_Npip_sec[sec - 1]->Fill(_e->W());
    MM_Npip_sec[sec - 1]->Fill(_e->MM());
  }
}

void Histogram::Write_WvsQ2() {
  for (short i = 0; i < 3; i++) {
    WQ2_det[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    WQ2_det[i]->SetYTitle("Q^{2} [Momentum Transfer] (GeV^{2})");
    WQ2_det[i]->SetOption("COLZ1");
    if (WQ2_det[i]->GetEntries()) WQ2_det[i]->Write();
    W_det[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    if (W_det[i]->GetEntries()) W_det[i]->Write();
  }
  auto WvsQ2_can = std::make_unique<TCanvas>("WvsQ2_can", "W vs Q2 sectors", 1920, 1080);
  WvsQ2_can->Divide(3, 2);
  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_sec[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    W_vs_q2_sec[i]->SetOption("COLZ1");
    WvsQ2_can->cd(i + 1);
    W_vs_q2_sec[i]->Draw("same");
  }
  WvsQ2_can->Write();

  auto W_can = std::make_unique<TCanvas>("W_can", "W sectors", 1920, 1080);
  W_can->Divide(3, 2);
  for (short i = 0; i < num_sectors; i++) {
    W_sec[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    W_can->cd(i + 1);
    // W_sec[i]->Fit("gaus", "QMR+", "QMR+", 0.5, 1.2);
    W_sec[i]->Draw("same");
  }
  W_can->Write();

  W_vs_q2->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
  W_vs_q2->SetYTitle("Q^{2} (GeV^{2})");
  W_vs_q2->SetOption("COLZ1");
  if (W_vs_q2->GetEntries()) W_vs_q2->Write();

  W_hist->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
  if (W_hist->GetEntries()) W_hist->Write();

  Q2_hist->SetXTitle("Q^{2} (GeV^{2})");
  if (Q2_hist->GetEntries()) Q2_hist->Write();

  W_vs_q2_singlePip->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
  W_vs_q2_singlePip->SetYTitle("Q^{2} (GeV^{2})");
  W_vs_q2_singlePip->SetOption("COLZ1");
  if (W_vs_q2_singlePip->GetEntries()) W_vs_q2_singlePip->Write();

  W_hist_singlePip->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
  if (W_hist_singlePip->GetEntries()) W_hist_singlePip->Write();

  Q2_hist_singlePip->SetXTitle("Q^{2} (GeV^{2})");
  if (Q2_hist_singlePip->GetEntries()) Q2_hist_singlePip->Write();

  if (MM_neutron->GetEntries()) MM_neutron->Write();

  auto wvsq2_sec = RootOutputFile->mkdir("wvsq2_sec");
  wvsq2_sec->cd();
  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_sec[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    W_vs_q2_sec[i]->SetOption("COLZ1");
    W_vs_q2_sec[i]->Write();
  }
  auto w_sec = RootOutputFile->mkdir("w_sec");
  w_sec->cd();
  for (short i = 0; i < num_sectors; i++) {
    W_sec[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    W_sec[i]->Write();
  }
  auto singlePip_sec = RootOutputFile->mkdir("singlePip_sec");
  singlePip_sec->cd();
  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_singlePip_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_singlePip_sec[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    W_vs_q2_singlePip_sec[i]->SetOption("COLZ1");
    if (W_vs_q2_singlePip_sec[i]->GetEntries()) W_vs_q2_singlePip_sec[i]->Write();
  }

  for (short i = 0; i < num_sectors; i++) {
    W_vs_MM_singlePip[i]->SetOption("COLZ1");
    W_vs_MM_singlePip[i]->SetYTitle("MM (GeV)");
    W_vs_MM_singlePip[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    if (W_vs_MM_singlePip[i]->GetEntries()) W_vs_MM_singlePip[i]->Write();
  }

  for (short i = 0; i < num_sectors; i++) {
    W_singlePip_sec[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    if (W_singlePip_sec[i]->GetEntries()) W_singlePip_sec[i]->Write();
  }

  for (short i = 0; i < num_sectors; i++) {
    if (MM_neutron_sec[i]->GetEntries()) MM_neutron_sec[i]->Fit("gaus", "QMR+", "QMR+", 0.7, 1.1);
    MM_neutron_sec[i]->SetXTitle("Mass (GeV)");
    if (MM_neutron_sec[i]->GetEntries()) MM_neutron_sec[i]->Write();
  }

  auto Npip_sec = RootOutputFile->mkdir("Npip_sec");
  Npip_sec->cd();

  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_Npip_sec[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_Npip_sec[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    W_vs_q2_Npip_sec[i]->SetOption("COLZ1");
    W_vs_q2_Npip_sec[i]->Write();
  }
  for (short i = 0; i < num_sectors; i++) {
    W_Npip_sec[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    W_Npip_sec[i]->Write();
  }
  for (short i = 0; i < num_sectors; i++) {
    MM_Npip_sec[i]->SetXTitle("Mass (GeV)");
    MM_Npip_sec[i]->Write();
  }
}

void Histogram::Write_WvsQ2MC() {
  auto MC_dir = RootOutputFile->mkdir("MC");
  MC_dir->cd();
  for (short i = 0; i < 3; i++) {
    WQ2_det_MC[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    WQ2_det_MC[i]->SetYTitle("Q^{2} [Momentum Transfer] (GeV^{2})");
    WQ2_det_MC[i]->SetOption("COLZ1");
    if (WQ2_det_MC[i]->GetEntries()) WQ2_det_MC[i]->Write();
    W_det_MC[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    if (W_det_MC[i]->GetEntries()) W_det_MC[i]->Write();
  }
  auto WvsQ2_can = std::make_unique<TCanvas>("WvsQ2_can_MC", "W vs Q2 sectors", 1920, 1080);
  WvsQ2_can->Divide(3, 2);
  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_sec_MC[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_sec_MC[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    W_vs_q2_sec_MC[i]->SetOption("COLZ1");
    WvsQ2_can->cd(i + 1);
    W_vs_q2_sec_MC[i]->Draw("same");
  }
  WvsQ2_can->Write();

  auto W_can = std::make_unique<TCanvas>("W_can_MC", "W sectors", 1920, 1080);
  W_can->Divide(3, 2);
  for (short i = 0; i < num_sectors; i++) {
    W_sec_MC[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    W_can->cd(i + 1);
    // W_sec[i]->Fit("gaus", "QMR+", "QMR+", 0.5, 1.2);
    W_sec_MC[i]->Draw("same");
  }
  W_can->Write();

  W_vs_q2_MC->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
  W_vs_q2_MC->SetYTitle("Q^{2} (GeV^{2})");
  W_vs_q2_MC->SetOption("COLZ1");
  if (W_vs_q2_MC->GetEntries()) W_vs_q2_MC->Write();

  W_MC_hist->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
  if (W_MC_hist->GetEntries()) W_MC_hist->Write();

  Q2_MC_hist->SetXTitle("Q^{2} (GeV^{2})");
  if (Q2_MC_hist->GetEntries()) Q2_MC_hist->Write();

  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_sec_MC[i]->SetYTitle("Q^{2} (GeV^{2})");
    W_vs_q2_sec_MC[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    W_vs_q2_sec_MC[i]->SetOption("COLZ1");
    W_vs_q2_sec_MC[i]->Write();
  }

  for (short i = 0; i < num_sectors; i++) {
    W_sec_MC[i]->SetXTitle("W [Invariant Mass] (GeV/c^{2})");
    W_sec_MC[i]->Write();
  }
}

void Histogram::makeHists_sector() {
  for (short i = 0; i < 3; i++) {
    W_det[i] = std::make_shared<TH1D>(Form("W_det_%d", i + 1), Form("W detector: %d", i + 1), bins, zero, w_max);
    W_det_MC[i] =
        std::make_shared<TH1D>(Form("W_det_MC_%d", i + 1), Form("W detector (MC): %d", i + 1), bins, zero, w_max);
    if (i == 0) {
      WQ2_det[i] = std::make_shared<TH2D>(Form("WQ2_det_%d", i + 1), Form("W vs Q^{2} detector: %d", i + 1), bins, zero,
                                          w_max, bins, zero, 0.5);
      WQ2_det_MC[i] = std::make_shared<TH2D>(Form("WQ2_det_MC_%d", i + 1), Form("W vs Q^{2} detector: %d", i + 1), bins,
                                             zero, w_max, bins, zero, 0.5);
    } else {
      WQ2_det[i] = std::make_shared<TH2D>(Form("WQ2_det_%d", i + 1), Form("W vs Q^{2} detector: %d", i + 1), bins, zero,
                                          w_max, bins, zero, q2_max);
      WQ2_det_MC[i] = std::make_shared<TH2D>(Form("WQ2_det_MC_%d", i + 1), Form("W vs Q^{2} detector: %d", i + 1), bins,
                                             zero, w_max, bins, zero, q2_max);
    }
  }

  for (short i = 0; i < num_sectors; i++) {
    W_vs_q2_sec[i] = std::make_shared<TH2D>(Form("wvsq2_sec_%d", i + 1), Form("W vs Q^{2} Sector: %d", i + 1), bins,
                                            zero, w_max, bins, zero, q2_max);

    W_sec[i] = std::make_shared<TH1D>(Form("w_sec_%d", i + 1), Form("W Sector: %d", i + 1), bins, zero, w_max);

    W_vs_q2_sec_MC[i] = std::make_shared<TH2D>(Form("wvsq2_sec_MC_%d", i + 1), Form("W vs Q^{2} Sector: %d", i + 1),
                                               bins, zero, w_max, bins, zero, q2_max);

    W_sec_MC[i] = std::make_shared<TH1D>(Form("w_sec_MC_%d", i + 1), Form("W Sector: %d", i + 1), bins, zero, w_max);

    W_vs_q2_singlePip_sec[i] =
        std::make_shared<TH2D>(Form("wvsq2_sec_singlePip_%d", i + 1), Form("W vs Q^{2} W_singlePip Sector: %d", i + 1),
                               bins, zero, w_max, bins, zero, q2_max);

    W_singlePip_sec[i] = std::make_shared<TH1D>(Form("w_sec_singlePip_%d", i + 1),
                                                Form("W singlePip Sector: %d", i + 1), bins, zero, w_max);

    W_vs_q2_Npip_sec[i] =
        std::make_shared<TH2D>(Form("wvsq2_sec_Npip_%d", i + 1), Form("W vs Q^{2} W_Npip Sector: %d", i + 1), bins,
                               zero, w_max, bins, zero, q2_max);

    W_Npip_sec[i] =
        std::make_shared<TH1D>(Form("w_sec_Npip_%d", i + 1), Form("W Npip Sector: %d", i + 1), bins, zero, w_max);

    MM_neutron_sec[i] =
        std::make_shared<TH1D>(Form("MM_Sec_%d", i + 1), Form("MM neutron Sector: %d", i + 1), bins, zero, 4.0);

    MM_Npip_sec[i] = std::make_shared<TH1D>(Form("MM_Npip_Sec_%d", i + 1), Form("MM^{2} neutron pip Sector: %d", i + 1),
                                            bins, zero, 4.0);
    W_vs_MM_singlePip[i] =
        std::make_shared<TH2D>(Form("W_vs_MM_singlePip_%d", i + 1), Form("W_vs_MM_singlePip_%d", i + 1), bins, zero,
                               w_max, bins, -q2_max, q2_max);
  }
}

void Histogram::makeHists_deltat() {
  for (short sec = 0; sec < num_sectors; sec++) {
    delta_t_pip[sec] = std::make_shared<TH2D>(Form("delta_t_pip_%d", sec), Form("#Deltat #pi^{+} Sector %d", sec), bins,
                                              0, 6.0, bins, -0.6, 0.6);
  }

  std::string tof = "";
  for (short p = 0; p < particle_num; p++) {
    for (short c = 0; c < charge_num; c++) {
      for (short i = 0; i < with_id_num; i++) {
        tof = "ftof";
        delta_t_hist[p][c][i][0] =
            std::make_shared<TH2D>(Form("delta_t_%s_%s_%s_%s", tof.c_str(), particle_name[p].c_str(),
                                        charge_name[c].c_str(), id_name[i].c_str()),
                                   Form("#Deltat %s %s %s %s", tof.c_str(), particle_name[p].c_str(),
                                        charge_name[c].c_str(), id_name[i].c_str()),
                                   bins, p_min, p_max, bins, Dt_min, Dt_max);

        tof = "ctof";
        delta_t_hist[p][c][i][1] =
            std::make_shared<TH2D>(Form("delta_t_%s_%s_%s_%s", tof.c_str(), particle_name[p].c_str(),
                                        charge_name[c].c_str(), id_name[i].c_str()),
                                   Form("#Deltat %s %s %s %s", tof.c_str(), particle_name[p].c_str(),
                                        charge_name[c].c_str(), id_name[i].c_str()),
                                   bins, 0, 3.0, bins, -6.0, 6.0);
      }
    }
  }
}

void Histogram::Fill_deltat_pi(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt, int part) {
  auto _cuts = std::make_unique<Cuts>(data, dt);
  int charge = data->charge(part);
  bool fc = dt->ctof();
  int pid = data->pid(part);
  float mom = data->p(part);
  float time = NAN;
  if (fc)
    time = dt->dt_ctof_Pi();
  else
    time = dt->dt_Pi();

  if (charge == 1) {
    delta_t_hist[1][0][0][fc]->Fill(mom, time);
    if (_cuts->IsPip(part)) {
      delta_t_hist[1][0][1][fc]->Fill(mom, time);
      if (data->dc_sec(part) >= 1 && data->dc_sec(part) <= 6) delta_t_pip[data->dc_sec(part) - 1]->Fill(mom, time);
    } else
      delta_t_hist[1][0][2][fc]->Fill(mom, time);
  } else if (charge == -1) {
    delta_t_hist[1][1][0][fc]->Fill(mom, time);
    if (_cuts->IsPim(part))
      delta_t_hist[1][1][1][fc]->Fill(mom, time);
    else
      delta_t_hist[1][1][2][fc]->Fill(mom, time);
  }
}

void Histogram::Fill_deltat_prot(const std::shared_ptr<Branches12>& data, const std::shared_ptr<Delta_T>& dt,
                                 int part) {
  auto _cuts = std::make_unique<Cuts>(data, dt);
  int charge = data->charge(part);
  bool fc = dt->ctof();
  int pid = data->pid(part);
  float mom = data->p(part);
  float time = NAN;

  if (fc)
    time = dt->dt_ctof_P();
  else
    time = dt->dt_P();

  if (charge == 1) {
    delta_t_hist[2][0][0][fc]->Fill(mom, time);
    if (_cuts->IsProton(part))
      delta_t_hist[2][0][1][fc]->Fill(mom, time);
    else
      delta_t_hist[2][0][2][fc]->Fill(mom, time);

    delta_t_hist[2][1][0][fc]->Fill(mom, time);
    if (pid == PROTON)
      delta_t_hist[2][1][1][fc]->Fill(mom, time);
    else
      delta_t_hist[2][1][2][fc]->Fill(mom, time);
  }
}

void Histogram::Write_deltat() {
  TDirectory* ftof_folder = RootOutputFile->mkdir("ftof");
  ftof_folder->cd();
  for (short sec = 0; sec < num_sectors; sec++) {
    delta_t_pip[sec]->SetXTitle("Momentum (GeV)");
    delta_t_pip[sec]->SetYTitle("#Deltat");
    delta_t_pip[sec]->SetOption("COLZ1");
    delta_t_pip[sec]->Write();
  }

  for (short p = 0; p < particle_num; p++) {
    for (short c = 0; c < charge_num; c++) {
      for (short i = 0; i < with_id_num; i++) {
        delta_t_hist[p][c][i][0]->SetXTitle("Momentum (GeV)");
        delta_t_hist[p][c][i][0]->SetYTitle("#Deltat");
        delta_t_hist[p][c][i][0]->SetOption("COLZ1");
        if (delta_t_hist[p][c][i][0]->GetEntries() > 1) delta_t_hist[p][c][i][0]->Write();
      }
    }
  }
  TDirectory* ctof_folder = RootOutputFile->mkdir("ctof");
  ctof_folder->cd();
  for (short p = 0; p < particle_num; p++) {
    for (short c = 0; c < charge_num; c++) {
      for (short i = 0; i < with_id_num; i++) {
        delta_t_hist[p][c][i][1]->SetXTitle("Momentum (GeV)");
        delta_t_hist[p][c][i][1]->SetYTitle("#Deltat");
        delta_t_hist[p][c][i][1]->SetOption("COLZ1");
        if (delta_t_hist[p][c][i][1]->GetEntries() > 1) delta_t_hist[p][c][i][1]->Write();
      }
    }
  }
}

void Histogram::makeHists_MomVsBeta() {
  for (short p = 0; p < particle_num; p++) {
    for (short c = 0; c < charge_num; c++) {
      for (short i = 0; i < with_id_num; i++) {
        momvsbeta_hist[p][c][i] = std::make_shared<TH2D>(
            Form("mom_vs_beta_%s_%s_%s", particle_name[p].c_str(), charge_name[c].c_str(), id_name[i].c_str()),
            Form("Momentum vs #beta %s %s %s", particle_name[p].c_str(), charge_name[c].c_str(), id_name[i].c_str()),
            bins, p_min, p_max, bins, zero, 1.2);
      }
    }
  }
}

void Histogram::Fill_MomVsBeta(const std::shared_ptr<Branches12>& data, int part) {
  int good_ID = 0;
  float beta = data->beta(part);
  float mom = data->p(part);
  int charge = data->charge(part);
  int pid = data->pid(part);
  if (beta != 0) {
    momentum->Fill(mom);
    for (short p = 0; p < particle_num; p++) {
      switch (p) {
        case 0:
          good_ID = ELECTRON;
          break;
        case 1:
          good_ID = PIP;
          break;
        case 2:
          good_ID = PROTON;
          break;
        case 3:
          good_ID = KP;
          break;
      }

      momvsbeta_hist[p][0][0]->Fill(mom, beta);
      if (good_ID == abs(pid)) {
        momvsbeta_hist[p][0][1]->Fill(mom, beta);
      } else {
        momvsbeta_hist[p][0][2]->Fill(mom, beta);
      }

      if (charge == -1) {
        momvsbeta_hist[p][2][0]->Fill(mom, beta);
        if (-good_ID == pid) {
          momvsbeta_hist[p][2][1]->Fill(mom, beta);
        } else {
          momvsbeta_hist[p][2][2]->Fill(mom, beta);
        }
      } else if (charge == 1) {
        momvsbeta_hist[p][1][0]->Fill(mom, beta);
        if (good_ID == pid) {
          momvsbeta_hist[p][1][1]->Fill(mom, beta);
        } else {
          momvsbeta_hist[p][1][2]->Fill(mom, beta);
        }
      }
    }
  }
}

void Histogram::Write_MomVsBeta() {
  momentum->SetXTitle("Momentum (GeV)");
  momentum->Write();
  for (short p = 0; p < particle_num; p++) {
    for (short c = 0; c < charge_num; c++) {
      for (short i = 0; i < with_id_num; i++) {
        momvsbeta_hist[p][c][i]->SetXTitle("Momentum (GeV)");
        momvsbeta_hist[p][c][i]->SetYTitle("#beta");
        momvsbeta_hist[p][c][i]->SetOption("COLZ1");
        momvsbeta_hist[p][c][i]->Write();
      }
    }
  }
}

void Histogram::Fill_EC(const std::shared_ptr<Branches12>& data) {
  EC_sampling_fraction->Fill(data->p(0), data->ec_tot_energy(0) / data->p(0));
}

void Histogram::Write_EC() {
  EC_sampling_fraction->SetXTitle("Momentum (GeV)");
  EC_sampling_fraction->SetYTitle("Sampling Fraction");
  EC_sampling_fraction->SetOption("COLZ1");
  EC_sampling_fraction->Write();
}
