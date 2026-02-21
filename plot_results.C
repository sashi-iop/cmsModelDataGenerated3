/*
 * plot_results.C — ROOT plotting macro for CMS Model simulation & digitisation
 *
 * Usage:
 *   root -l plot_results.C                          (uses defaults)
 *   root -l 'plot_results.C("myfile.root")'         (simulation only)
 *   root -l 'plot_results.C("sim.root","digi.root")'(both)
 *
 * Produces:
 *   plots_simulation.pdf  — simulation histograms
 *   plots_digitisation.pdf — digitisation results (if digi file provided)
 */

void plot_results(const char *simFile = "test_run10.root",
                  const char *digiFile = "") {
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(1);
  gStyle->SetPalette(kViridis);
  gStyle->SetTitleSize(0.05, "XY");
  gStyle->SetLabelSize(0.04, "XY");

  // ================================================================
  //  PART 1: SIMULATION OUTPUT (T1 tree + histograms)
  // ================================================================
  TFile *fSim = TFile::Open(simFile, "READ");
  if (!fSim || fSim->IsZombie()) {
    cout << "ERROR: Cannot open simulation file: " << simFile << endl;
    return;
  }
  cout << "Opened simulation file: " << simFile << endl;

  // --- List all keys in the file for reference ---
  cout << "\n--- Contents of " << simFile << " ---" << endl;
  fSim->ls();
  cout << "-----------------------------------\n" << endl;

  TTree *T1 = (TTree *)fSim->Get("T1");

  // ==============================================
  // Page 1: ECAL Histograms
  // ==============================================
  TCanvas *c1 = new TCanvas("c1", "ECAL Histograms", 1200, 900);
  c1->Divide(2, 2);

  TH1F *h_ecal = (TH1F *)fSim->Get("ecalenergy");
  TH2F *h_ecal2d = (TH2F *)fSim->Get("h2decalenergy");
  TH1F *h_posRE = (TH1F *)fSim->Get("pPosRE");
  TH1F *h_posEE = (TH1F *)fSim->Get("pPosEE");

  c1->cd(1);
  if (h_ecal) {
    h_ecal->SetLineColor(kBlue + 1);
    h_ecal->SetLineWidth(2);
    h_ecal->Draw();
  } else
    cout << "  [SKIP] ecalenergy not found" << endl;

  c1->cd(2);
  if (h_ecal2d) {
    h_ecal2d->Draw("COLZ");
  } else
    cout << "  [SKIP] h2decalenergy not found" << endl;

  c1->cd(3);
  if (h_posRE) {
    h_posRE->SetLineColor(kGreen + 2);
    h_posRE->SetLineWidth(2);
    h_posRE->Draw();
  } else
    cout << "  [SKIP] pPosRE not found" << endl;

  c1->cd(4);
  if (h_posEE) {
    h_posEE->SetLineColor(kRed + 1);
    h_posEE->SetLineWidth(2);
    h_posEE->Draw();
  } else
    cout << "  [SKIP] pPosEE not found" << endl;

  c1->SaveAs("plots_simulation.pdf("); // Open multi-page PDF

  // ==============================================
  // Page 2: HCAL per-layer energy histograms
  // ==============================================
  TCanvas *c2 = new TCanvas("c2", "HCAL Layer Energy", 1400, 1000);

  // Try to find how many HCAL layer histograms exist
  int nhcalFound = 0;
  for (int i = 0; i < 20; i++) {
    TH1F *h = (TH1F *)fSim->Get(Form("hcalenergy_L%d", i));
    if (h)
      nhcalFound++;
    else
      break;
  }
  if (nhcalFound == 0)
    nhcalFound = 1;

  int nx = (int)ceil(sqrt((double)nhcalFound));
  int ny = (int)ceil((double)nhcalFound / nx);
  c2->Divide(nx, ny);

  int colors[] = {kBlue + 1,   kRed + 1,    kGreen + 2,  kMagenta + 1,
                  kOrange + 1, kCyan + 1,   kViolet + 1, kTeal + 1,
                  kPink + 1,   kSpring + 2, kAzure + 1,  kYellow + 2,
                  kGray + 1,   kBlue - 4,   kRed - 4,    kGreen - 4,
                  kMagenta - 4};

  for (int i = 0; i < nhcalFound; i++) {
    c2->cd(i + 1);
    TH1F *h = (TH1F *)fSim->Get(Form("hcalenergy_L%d", i));
    if (h) {
      h->SetLineColor(colors[i % 17]);
      h->SetLineWidth(2);
      h->Draw();
    }
  }
  c2->SaveAs("plots_simulation.pdf");

  // ==============================================
  // Page 3: HCAL 2D energy maps (eta-phi)
  // ==============================================
  TCanvas *c2b = new TCanvas("c2b", "HCAL 2D Energy Maps", 1400, 1000);
  c2b->Divide(nx, ny);

  for (int i = 0; i < nhcalFound; i++) {
    c2b->cd(i + 1);
    TH2F *h = (TH2F *)fSim->Get(Form("hcal2denergy_L%d", i));
    if (h) {
      h->Draw("COLZ");
    }
  }
  c2b->SaveAs("plots_simulation.pdf");

  // ==============================================
  // Page 4: HCAL position histograms
  // ==============================================
  TCanvas *c3 = new TCanvas("c3", "HCAL Position", 1200, 400);
  c3->Divide(3, 1);

  TH1F *h_posRH = (TH1F *)fSim->Get("pPosRH");
  TH1F *h_posEH = (TH1F *)fSim->Get("pPosEH");
  TH1F *h_posPH = (TH1F *)fSim->Get("pPosPH");

  c3->cd(1);
  if (h_posRH) {
    h_posRH->SetLineColor(kBlue + 1);
    h_posRH->SetLineWidth(2);
    h_posRH->Draw();
  }

  c3->cd(2);
  if (h_posEH) {
    h_posEH->SetLineColor(kRed + 1);
    h_posEH->SetLineWidth(2);
    h_posEH->Draw();
  }

  c3->cd(3);
  if (h_posPH) {
    h_posPH->SetLineColor(kGreen + 2);
    h_posPH->SetLineWidth(2);
    h_posPH->Draw();
  }

  c3->SaveAs("plots_simulation.pdf");

  // ==============================================
  // Page 5: Tracker layer energy
  // ==============================================
  TCanvas *c4 = new TCanvas("c4", "Tracker Layer Energy", 1400, 800);
  int ntrkFound = 0;
  for (int i = 0; i < 15; i++) {
    TH1F *h = (TH1F *)fSim->Get(Form("trkenergy_L%d", i));
    if (h)
      ntrkFound++;
    else
      break;
  }
  if (ntrkFound > 0) {
    int tnx = (int)ceil(sqrt((double)ntrkFound));
    int tny = (int)ceil((double)ntrkFound / tnx);
    c4->Divide(tnx, tny);
    for (int i = 0; i < ntrkFound; i++) {
      c4->cd(i + 1);
      TH1F *h = (TH1F *)fSim->Get(Form("trkenergy_L%d", i));
      if (h) {
        h->SetLineColor(colors[i % 17]);
        h->SetLineWidth(2);
        h->Draw();
      }
    }
  }
  c4->SaveAs("plots_simulation.pdf");

  // ==============================================
  // Page 6: Calibration Histograms (HCAL-only and ECAL+HCAL)
  // ==============================================
  TCanvas *c5 = new TCanvas("c5", "Calibration", 1200, 500);
  c5->Divide(2, 1);

  TH1F *h_hcal_total = (TH1F *)fSim->Get("hcal_total_energy");
  TH1F *h_ecal_hcal = (TH1F *)fSim->Get("ecal_hcal_total_energy");

  c5->cd(1);
  if (h_hcal_total) {
    h_hcal_total->SetLineColor(kOrange + 1);
    h_hcal_total->SetLineWidth(2);
    h_hcal_total->SetFillColor(kOrange - 9);
    h_hcal_total->SetFillStyle(3004);
    h_hcal_total->Draw();
    // Fit with Gaussian to extract calibration peak
    if (h_hcal_total->GetEntries() > 10) {
      h_hcal_total->Fit("gaus", "Q");
    }
  } else {
    cout << "  [SKIP] hcal_total_energy not found" << endl;
  }

  c5->cd(2);
  if (h_ecal_hcal) {
    h_ecal_hcal->SetLineColor(kBlue + 1);
    h_ecal_hcal->SetLineWidth(2);
    h_ecal_hcal->SetFillColor(kBlue - 9);
    h_ecal_hcal->SetFillStyle(3004);
    h_ecal_hcal->Draw();
    // Fit with Gaussian to extract calibration peak
    if (h_ecal_hcal->GetEntries() > 10) {
      h_ecal_hcal->Fit("gaus", "Q");
    }
  } else {
    cout << "  [SKIP] ecal_hcal_total_energy not found" << endl;
  }

  c5->SaveAs("plots_simulation.pdf");

  // ==============================================
  // Page 7: T1 tree distributions
  // ==============================================
  if (T1) {
    TCanvas *c6 = new TCanvas("c6", "Generated Particle Info", 1200, 800);
    c6->Divide(2, 2);

    c6->cd(1);
    T1->Draw("momin[0]>>h_momin(100, 0, 300)", "", "");
    TH1F *hm = (TH1F *)gDirectory->Get("h_momin");
    if (hm) {
      hm->SetTitle("Generated Momentum;Momentum (GeV);Events");
      hm->SetLineWidth(2);
    }

    c6->cd(2);
    T1->Draw("nsimhtHL>>h_nhits(100, 0, 500)", "", "");
    TH1F *hn = (TH1F *)gDirectory->Get("h_nhits");
    if (hn) {
      hn->SetTitle("HCAL Hits per Event;N hits;Events");
      hn->SetLineWidth(2);
    }

    c6->cd(3);
    T1->Draw("nsimhtEC>>h_necal(100, 0, 500)", "", "");
    TH1F *hne = (TH1F *)gDirectory->Get("h_necal");
    if (hne) {
      hne->SetTitle("ECAL Hits per Event;N hits;Events");
      hne->SetLineWidth(2);
    }

    c6->cd(4);
    T1->Draw("pidin[0]>>h_pid(500, -250, 250)", "", "");
    TH1F *hp = (TH1F *)gDirectory->Get("h_pid");
    if (hp) {
      hp->SetTitle("Particle ID;PID;Events");
      hp->SetLineWidth(2);
    }

    c6->SaveAs("plots_simulation.pdf");

    // Page 8: Energy from packed words
    TCanvas *c7 = new TCanvas("c7", "Packed Word Analysis", 1200, 500);
    c7->Divide(2, 1);

    c7->cd(1);
    // Extract energy from lower 17 bits of detidHL
    T1->Draw("(detidHL & 0x1FFFF)/1000.0>>h_packed_energy(200, 0, 50)", "", "");
    TH1F *hpe = (TH1F *)gDirectory->Get("h_packed_energy");
    if (hpe) {
      hpe->SetTitle("HCAL Hit Energy (from packed word);Energy (GeV);Hits");
      hpe->SetLineColor(kRed + 1);
      hpe->SetLineWidth(2);
    }

    c7->cd(2);
    // Extract depth from bits 0-2 of cell ID (which is bits 17-19 of packed
    // word)
    T1->Draw("(detidHL >> 17) & 0x7>>h_depth(10, 0, 10)", "", "");
    TH1F *hd = (TH1F *)gDirectory->Get("h_depth");
    if (hd) {
      hd->SetTitle("HCAL Depth Distribution (3-bit);Depth;Hits");
      hd->SetLineColor(kGreen + 2);
      hd->SetLineWidth(2);
    }

    c7->SaveAs("plots_simulation.pdf)"); // Close multi-page PDF
  }

  cout << "\n=== Saved: plots_simulation.pdf ===" << endl;

  // ================================================================
  //  PART 2: DIGITISATION OUTPUT (T2 tree)
  // ================================================================
  TString digiStr(digiFile);
  if (digiStr.Length() == 0) {
    // Try to auto-detect a digitisation file
    TString tryFile = "test_pion_klong_xwoecl_40_200_25_100.root";
    if (!gSystem->AccessPathName(tryFile.Data())) {
      digiStr = tryFile;
    }
  }

  if (digiStr.Length() == 0) {
    cout << "\nNo digitisation file found. Skipping digitisation plots."
         << endl;
    cout << "Run: ./ecal_hcal_digitisation 40 200 25 100" << endl;
    return;
  }

  TFile *fDigi = TFile::Open(digiStr.Data(), "READ");
  if (!fDigi || fDigi->IsZombie()) {
    cout << "Cannot open digitisation file: " << digiStr << endl;
    return;
  }
  cout << "\nOpened digitisation file: " << digiStr << endl;

  TTree *T2 = (TTree *)fDigi->Get("T2");
  if (!T2) {
    cout << "ERROR: T2 tree not found in " << digiStr << endl;
    return;
  }

  TCanvas *c8 = new TCanvas("c8", "Digitisation Results", 1200, 900);
  c8->Divide(2, 2);

  // Total energy: simulated vs digitised (ECAL + HCAL)
  c8->cd(1);
  T2->Draw("totsimenr>>h_totsim(200, 0, 300)", "", "");
  TH1F *hts = (TH1F *)gDirectory->Get("h_totsim");
  if (hts) {
    hts->SetTitle("ECAL+HCAL Simulated Energy;Energy (GeV);Events");
    hts->SetLineColor(kBlue + 1);
    hts->SetLineWidth(2);
    if (hts->GetEntries() > 10)
      hts->Fit("gaus", "Q");
  }

  c8->cd(2);
  T2->Draw("totdigienr>>h_totdigi(200, 0, 300)", "", "");
  TH1F *htd = (TH1F *)gDirectory->Get("h_totdigi");
  if (htd) {
    htd->SetTitle("ECAL+HCAL Digitised Energy;Energy (GeV);Events");
    htd->SetLineColor(kRed + 1);
    htd->SetLineWidth(2);
    if (htd->GetEntries() > 10)
      htd->Fit("gaus", "Q");
  }

  // HCAL only: simulated vs digitised
  c8->cd(3);
  T2->Draw("hclsimenr>>h_hclsim(200, 0, 300)", "", "");
  TH1F *hhs = (TH1F *)gDirectory->Get("h_hclsim");
  if (hhs) {
    hhs->SetTitle("HCAL-only Simulated Energy;Energy (GeV);Events");
    hhs->SetLineColor(kGreen + 2);
    hhs->SetLineWidth(2);
    if (hhs->GetEntries() > 10)
      hhs->Fit("gaus", "Q");
  }

  c8->cd(4);
  T2->Draw("hcldigienr>>h_hcldigi(200, 0, 300)", "", "");
  TH1F *hhd = (TH1F *)gDirectory->Get("h_hcldigi");
  if (hhd) {
    hhd->SetTitle("HCAL-only Digitised Energy;Energy (GeV);Events");
    hhd->SetLineColor(kOrange + 1);
    hhd->SetLineWidth(2);
    if (hhd->GetEntries() > 10)
      hhd->Fit("gaus", "Q");
  }

  c8->SaveAs("plots_digitisation.pdf(");

  // Page 2: Overlay comparisons
  TCanvas *c9 = new TCanvas("c9", "Comparison Overlays", 1200, 500);
  c9->Divide(2, 1);

  // ECAL+HCAL: sim vs digi overlay
  c9->cd(1);
  T2->Draw("totsimenr>>h_ov_sim(200, 0, 300)", "", "");
  T2->Draw("totdigienr>>h_ov_digi(200, 0, 300)", "", "same");
  TH1F *hovs = (TH1F *)gDirectory->Get("h_ov_sim");
  TH1F *hovd = (TH1F *)gDirectory->Get("h_ov_digi");
  if (hovs) {
    hovs->SetLineColor(kBlue + 1);
    hovs->SetLineWidth(2);
    hovs->SetTitle("ECAL+HCAL: Simulated vs Digitised;Energy (GeV);Events");
  }
  if (hovd) {
    hovd->SetLineColor(kRed + 1);
    hovd->SetLineWidth(2);
    hovd->SetLineStyle(2);
  }

  TLegend *leg1 = new TLegend(0.6, 0.7, 0.88, 0.85);
  if (hovs)
    leg1->AddEntry(hovs, "Simulated", "l");
  if (hovd)
    leg1->AddEntry(hovd, "Digitised", "l");
  leg1->Draw();

  // HCAL only: sim vs digi overlay
  c9->cd(2);
  T2->Draw("hclsimenr>>h_ov_hsim(200, 0, 300)", "", "");
  T2->Draw("hcldigienr>>h_ov_hdigi(200, 0, 300)", "", "same");
  TH1F *hovhs = (TH1F *)gDirectory->Get("h_ov_hsim");
  TH1F *hovhd = (TH1F *)gDirectory->Get("h_ov_hdigi");
  if (hovhs) {
    hovhs->SetLineColor(kGreen + 2);
    hovhs->SetLineWidth(2);
    hovhs->SetTitle("HCAL-only: Simulated vs Digitised;Energy (GeV);Events");
  }
  if (hovhd) {
    hovhd->SetLineColor(kOrange + 1);
    hovhd->SetLineWidth(2);
    hovhd->SetLineStyle(2);
  }

  TLegend *leg2 = new TLegend(0.6, 0.7, 0.88, 0.85);
  if (hovhs)
    leg2->AddEntry(hovhs, "Simulated", "l");
  if (hovhd)
    leg2->AddEntry(hovhd, "Digitised", "l");
  leg2->Draw();

  c9->SaveAs("plots_digitisation.pdf");

  // Page 3: Resolution plot (digi/sim ratio)
  TCanvas *c10 = new TCanvas("c10", "Energy Resolution", 1200, 500);
  c10->Divide(2, 1);

  c10->cd(1);
  T2->Draw("totdigienr/totsimenr>>h_ratio_tot(200, 0, 2)", "totsimenr>1", "");
  TH1F *hrt = (TH1F *)gDirectory->Get("h_ratio_tot");
  if (hrt) {
    hrt->SetTitle("ECAL+HCAL: Digitised/Simulated;E_{digi}/E_{sim};Events");
    hrt->SetLineColor(kBlue + 1);
    hrt->SetLineWidth(2);
    if (hrt->GetEntries() > 10)
      hrt->Fit("gaus", "Q");
  }

  c10->cd(2);
  T2->Draw("hcldigienr/hclsimenr>>h_ratio_hcl(200, 0, 2)", "hclsimenr>1", "");
  TH1F *hrh = (TH1F *)gDirectory->Get("h_ratio_hcl");
  if (hrh) {
    hrh->SetTitle("HCAL-only: Digitised/Simulated;E_{digi}/E_{sim};Events");
    hrh->SetLineColor(kOrange + 1);
    hrh->SetLineWidth(2);
    if (hrh->GetEntries() > 10)
      hrh->Fit("gaus", "Q");
  }

  c10->SaveAs("plots_digitisation.pdf)"); // Close PDF

  cout << "=== Saved: plots_digitisation.pdf ===" << endl;

  // Print summary
  cout << "\n========== SUMMARY ==========" << endl;
  if (hrt && hrt->GetFunction("gaus")) {
    TF1 *f = hrt->GetFunction("gaus");
    cout << "ECAL+HCAL resolution: mean = " << f->GetParameter(1)
         << ", sigma = " << f->GetParameter(2)
         << ", sigma/mean = " << f->GetParameter(2) / f->GetParameter(1)
         << endl;
  }
  if (hrh && hrh->GetFunction("gaus")) {
    TF1 *f = hrh->GetFunction("gaus");
    cout << "HCAL-only resolution: mean = " << f->GetParameter(1)
         << ", sigma = " << f->GetParameter(2)
         << ", sigma/mean = " << f->GetParameter(2) / f->GetParameter(1)
         << endl;
  }
  cout << "=============================" << endl;
}
