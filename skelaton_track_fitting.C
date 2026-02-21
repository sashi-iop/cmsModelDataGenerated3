/*
 * Track Finder using Hough Transformation
 *
 * Hands-on Session — Complete Implementation
 *
 * Features:
 *   1. Per-event display: signal (black) + noise (red) in (ρ,z) and (x,y)
 *   2. Hough transform in R-Z: slope vs intercept, color-coded
 *   3. Hough transform in R-φ: curvature vs φ₀
 *   4. Resolution plots: ΔpT(GeV), Δφ(mrad), Δθ(mrad) with log scale
 *   5. Noise comparison mode (nNoiseHit)
 *   6. Multiple trajectory finding
 *
 * Usage:
 *   cd build
 *   echo "test_run10.root 10000" > skeleton_track.log
 *   ./track_finder
 *
 * Output:
 *   output_skeleton.root  — T2 tree + histograms
 *   track_event_display.pdf — per-event hit and Hough space plots
 *   track_resolution.pdf  — resolution and summary plots
 */

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TObject.h"
#include "TPad.h"
#include "TPostScript.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TTree.h"
#include "TVector3.h"
#include <TRandom3.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

#define DELTARAY // Removal of hits from delta ray

using namespace std;

struct fourparam {
  float hit1;
  float hit2;
  float slope;
  float inter;
  int type; // 0=signal, 1=mixed, 2=noise
};

struct TrackResult {
  double theta;
  double phi;
  double pt;
  vector<int> hit_indices;
};

const double pival = acos(-1.0);
const double twopi = 2 * pival;
const double pibytwo = pival / 2.;
const double pibyfour = pival / 4.;
const double radtodeg = 180. / pival;

//=========================================================
// CONFIGURABLE PARAMETERS
//=========================================================
int nNoiseHit = 0;               // Runtime: pass as CLI arg. e.g. ./track_finder 13
const int nsilayer = 13;
const int nmnhit = 3;            // Minimum hits for a track
const int nEventDisplay = 4;     // Number of events to display individually
const double hough_nsigma = 2.0; // σ cut for Hough cleanup
const int hough_maxiter = 10;    // Max iterations for cleanup
//=========================================================

double sirad[nsilayer] = {0.029, 0.068, 0.109, 0.160, 0.220, 0.290, 0.370,
                          0.450, 0.520, 0.590, 0.660, 0.730, 0.800};
double momconv = 0.299792 * 3.8; // c * B (3.8 T)
double extrad = 0.0289;

double phireso2fact = 2.5e-9 / 12.0;
double zreso2fact = 4.0e-8 / 12.0;

TRandom3 *gRandom3 = new TRandom3();

//------------------------------------------------------------
// Helper: TVector3 from cylindrical (rho, phi, z)
//------------------------------------------------------------
TVector3 MakeRhoPhiZ(double rho, double phi, double z) {
  return TVector3(rho * cos(phi), rho * sin(phi), z);
}

TVector3 MakeRhoPhiTheta(double rho, double phi, double theta) {
  double z = rho / tan(theta);
  return TVector3(rho * cos(phi), rho * sin(phi), z);
}

double sfdph(double phi1, double phi2) {
  double dphi = phi2 - phi1;
  if (dphi > M_PI)
    dphi -= 2 * M_PI;
  if (dphi < -M_PI)
    dphi += 2 * M_PI;
  return dphi;
}

// Forward declarations
void fill_digipos_array(int nhit, unsigned int *detid, int *pdgid,
                        float *simenr, float *simtime, vector<TVector3> &hitpos,
                        vector<bool> &isSignal);
void straight_line_fit(vector<TVector3> finder, double *par, double *parerr);
bool GetPrediction(double radin, double *StateVector, double *Prediction,
                   int il, double *Q_k_minus, int ix);
void track_move_pt_align(double *b, double *tin, double q, double *x,
                         double *tout, double &dl);

//------------------------------------------------------------
// Hough Transform in R-Z plane
// Returns cleaned (slope, intercept) pairs and reconstructed theta
//------------------------------------------------------------
double hough_rz(vector<TVector3> &hits, vector<bool> &isSignal,
                vector<fourparam> &cleaned, vector<int> &selected_hits) {

  vector<fourparam> allpairs;
  int nhits = hits.size();

  // Generate all pairs
  for (int i = 0; i < nhits; i++) {
    for (int j = i + 1; j < nhits; j++) {
      double r1 = hits[i].Perp(), z1 = hits[i].Z();
      double r2 = hits[j].Perp(), z2 = hits[j].Z();
      if (fabs(r2 - r1) < 1e-6)
        continue;

      fourparam fp;
      fp.hit1 = i;
      fp.hit2 = j;
      fp.slope = (z2 - z1) / (r2 - r1); // cotθ
      fp.inter = z1 - fp.slope * r1;    // z₀

      // Classify: 0=signal(both signal), 1=mixed, 2=noise(both noise)
      if (isSignal[i] && isSignal[j])
        fp.type = 0;
      else if (!isSignal[i] && !isSignal[j])
        fp.type = 2;
      else
        fp.type = 1;

      allpairs.push_back(fp);
    }
  }

  // Iterative 2σ cleanup
  cleaned = allpairs;
  for (int iter = 0; iter < hough_maxiter; iter++) {
    if (cleaned.size() < 2)
      break;

    double ms = 0, mi = 0;
    for (auto &p : cleaned) {
      ms += p.slope;
      mi += p.inter;
    }
    ms /= cleaned.size();
    mi /= cleaned.size();

    double rs = 0, ri = 0;
    for (auto &p : cleaned) {
      rs += (p.slope - ms) * (p.slope - ms);
      ri += (p.inter - mi) * (p.inter - mi);
    }
    rs = sqrt(rs / cleaned.size());
    ri = sqrt(ri / cleaned.size());

    if (rs < 1e-10 && ri < 1e-10)
      break;

    vector<fourparam> tmp;
    for (auto &p : cleaned) {
      bool ok = true;
      if (rs > 1e-10 && fabs(p.slope - ms) > hough_nsigma * rs)
        ok = false;
      if (ri > 1e-10 && fabs(p.inter - mi) > hough_nsigma * ri)
        ok = false;
      if (ok)
        tmp.push_back(p);
    }
    if (tmp.size() == cleaned.size())
      break;
    cleaned = tmp;
  }

  // Collect selected hit indices
  set<int> sel;
  for (auto &p : cleaned) {
    sel.insert((int)p.hit1);
    sel.insert((int)p.hit2);
  }
  selected_hits.assign(sel.begin(), sel.end());

  // Mean slope → theta
  double mean_slope = 0;
  if (cleaned.size() > 0) {
    for (auto &p : cleaned)
      mean_slope += p.slope;
    mean_slope /= cleaned.size();
  }

  // Refine with straight line fit
  vector<TVector3> fitpts;
  for (int idx : selected_hits)
    fitpts.push_back(hits[idx]);
  if (fitpts.size() >= 2) {
    double par[2], parerr[3];
    straight_line_fit(fitpts, par, parerr);
    if (par[0] > -0.99)
      mean_slope = par[0];
  }

  double theta = atan2(1.0, mean_slope);
  if (theta < 0)
    theta += pival;
  return theta;
}

//------------------------------------------------------------
// Hough Transform in R-φ plane
// Returns reconstructed phi and pT
//------------------------------------------------------------
void hough_rphi(vector<TVector3> &hits, vector<bool> &isSignal,
                vector<int> &selected_hits, double &reco_phi, double &reco_pt) {

  vector<fourparam> allpairs;
  int nhits = hits.size();

  for (int i = 0; i < nhits; i++) {
    for (int j = i + 1; j < nhits; j++) {
      double r1 = hits[i].Perp(), phi1 = hits[i].Phi();
      double r2 = hits[j].Perp(), phi2 = hits[j].Phi();
      if (fabs(r2 - r1) < 1e-6)
        continue;

      double dphi = sfdph(phi1, phi2);
      double slope = dphi / (r2 - r1);  // ≈ κ/2
      double inter = phi1 - slope * r1; // ≈ φ₀

      fourparam fp;
      fp.hit1 = i;
      fp.hit2 = j;
      fp.slope = slope;
      fp.inter = inter;
      if (isSignal[i] && isSignal[j])
        fp.type = 0;
      else if (!isSignal[i] && !isSignal[j])
        fp.type = 2;
      else
        fp.type = 1;
      allpairs.push_back(fp);
    }
  }

  // Use only hits from R-Z selection
  set<int> rz_set(selected_hits.begin(), selected_hits.end());
  vector<fourparam> filtered;
  for (auto &p : allpairs) {
    if (rz_set.count((int)p.hit1) && rz_set.count((int)p.hit2))
      filtered.push_back(p);
  }
  if (filtered.empty())
    filtered = allpairs;

  // Iterative cleanup
  for (int iter = 0; iter < hough_maxiter; iter++) {
    if (filtered.size() < 2)
      break;
    double ms = 0, mi = 0;
    for (auto &p : filtered) {
      ms += p.slope;
      mi += p.inter;
    }
    ms /= filtered.size();
    mi /= filtered.size();
    double rs = 0, ri = 0;
    for (auto &p : filtered) {
      rs += (p.slope - ms) * (p.slope - ms);
      ri += (p.inter - mi) * (p.inter - mi);
    }
    rs = sqrt(rs / filtered.size());
    ri = sqrt(ri / filtered.size());
    if (rs < 1e-10 && ri < 1e-10)
      break;
    vector<fourparam> tmp;
    for (auto &p : filtered) {
      bool ok = true;
      if (rs > 1e-10 && fabs(p.slope - ms) > hough_nsigma * rs)
        ok = false;
      if (ri > 1e-10 && fabs(p.inter - mi) > hough_nsigma * ri)
        ok = false;
      if (ok)
        tmp.push_back(p);
    }
    if (tmp.size() == filtered.size())
      break;
    filtered = tmp;
  }

  double mean_slope = 0, mean_inter = 0;
  if (filtered.size() > 0) {
    for (auto &p : filtered) {
      mean_slope += p.slope;
      mean_inter += p.inter;
    }
    mean_slope /= filtered.size();
    mean_inter /= filtered.size();
  }

  reco_phi = mean_inter;
  if (fabs(mean_slope) > 1e-8)
    reco_pt = fabs(momconv / (2.0 * mean_slope));
  else
    reco_pt = 99999;
}

//============================================================
// MAIN — accepts optional CLI argument for noise hits
//   ./track_finder          (no noise)
//   ./track_finder 13       (13 noise hits per event)
//============================================================
int main(int argc, char *argv[]) {

  // Parse noise hits from command line
  if (argc > 1) {
    nNoiseHit = atoi(argv[1]);
    if (nNoiseHit < 0) nNoiseHit = 0;
  }
  cout << "Noise hits per event: " << nNoiseHit << endl;

  for (int ij = 0; ij < nsilayer; ij++)
    sirad[ij] += 0.00015;

  char datafile[100];
  unsigned irun, ievt, ngent, ntrkt;
  float ievt_wt;
  const unsigned int ngenmx = 50;
  int pidin[ngenmx];
  float momin[ngenmx], thein[ngenmx], phiin[ngenmx];
  float posxin[ngenmx], posyin[ngenmx], poszin[ngenmx];

  unsigned int nsimhtTk;
  const unsigned int nsimhtmxTk = 2000;
  unsigned int detidTk[nsimhtmxTk];
  float simtimeTk[nsimhtmxTk], simenrTk[nsimhtmxTk];
  int simpdgidTk[nsimhtmxTk];
  float simvxTk[nsimhtmxTk], simvyTk[nsimhtmxTk], simvzTk[nsimhtmxTk];
  float simpxTk[nsimhtmxTk], simpyTk[nsimhtmxTk], simpzTk[nsimhtmxTk];

  vector<TVector3> digihitpos;
  vector<bool> hitIsSignal; // true=signal, false=noise

  TFile *fileOut = new TFile("output_skeleton.root", "recreate");

  TTree *Tout = new TTree("T2", "Track Reconstruction Output");

  // Output variables
  float reco_theta, reco_phi, reco_pt;
  float gen_theta, gen_phi, gen_pt, gen_mom;
  float delta_theta_mrad, delta_phi_mrad, delta_pt_GeV;
  int nhits_used, nhits_total, ntrack_found;

  Tout->Branch("irun", &irun, "irun/i");
  Tout->Branch("ievt", &ievt, "ievt/i");
  Tout->Branch("gen_theta", &gen_theta, "gen_theta/F");
  Tout->Branch("gen_phi", &gen_phi, "gen_phi/F");
  Tout->Branch("gen_pt", &gen_pt, "gen_pt/F");
  Tout->Branch("gen_mom", &gen_mom, "gen_mom/F");
  Tout->Branch("reco_theta", &reco_theta, "reco_theta/F");
  Tout->Branch("reco_phi", &reco_phi, "reco_phi/F");
  Tout->Branch("reco_pt", &reco_pt, "reco_pt/F");
  Tout->Branch("delta_theta_mrad", &delta_theta_mrad, "delta_theta_mrad/F");
  Tout->Branch("delta_phi_mrad", &delta_phi_mrad, "delta_phi_mrad/F");
  Tout->Branch("delta_pt_GeV", &delta_pt_GeV, "delta_pt_GeV/F");
  Tout->Branch("nhits_used", &nhits_used, "nhits_used/I");
  Tout->Branch("nhits_total", &nhits_total, "nhits_total/I");
  Tout->Branch("ntrack_found", &ntrack_found, "ntrack_found/I");

  //============================================================
  // Resolution histograms (matching slide format exactly)
  //============================================================
  TH1F *h_delta_pt =
      new TH1F("h_delta_pt", ";#Delta p_{T} (GeV);Entries", 200, -50, 50);
  TH1F *h_delta_phi =
      new TH1F("h_delta_phi", ";#Delta#phi (mrad);Entries", 200, -5, 5);
  TH1F *h_delta_theta =
      new TH1F("h_delta_theta", ";#Delta#theta (mrad);Entries", 200, -5, 5);

  h_delta_pt->SetLineColor(kRed + 1);
  h_delta_pt->SetLineWidth(2);
  h_delta_phi->SetLineColor(kRed + 1);
  h_delta_phi->SetLineWidth(2);
  h_delta_theta->SetLineColor(kRed + 1);
  h_delta_theta->SetLineWidth(2);

  // Additional analysis histograms
  TH1F *h_nhits =
      new TH1F("h_nhits", "Hits per event;N_{hits};Events", 30, 0, 30);
  TH1F *h_reco_pt_dist = new TH1F(
      "h_reco_pt_dist", "Reconstructed p_{T};p_{T} (GeV);Events", 100, 0, 100);
  TH1F *h_ntracks = new TH1F(
      "h_ntracks", "Tracks found per event;N_{tracks};Events", 10, 0, 10);

  // Per-event display canvases (only for first N events)
  TCanvas *c_evtdisp = nullptr;

  gStyle->SetOptStat(1111);

  ifstream file_db;
  file_db.open("skeleton_track.log");
  if (!file_db.is_open()) {
    cout << "ERROR: Cannot open skeleton_track.log" << endl;
    cout << "Create: echo 'test_run10.root 10000' > skeleton_track.log" << endl;
    return 1;
  }

  int global_evt = 0;
  bool pdfOpened = false;

  while (!(file_db.eof())) {
    int nentrymx = -1;
    file_db >> datafile >> nentrymx;
    if (strstr(datafile, "#"))
      continue;
    if (file_db.eof())
      break;

    cout << "Opening: " << datafile << endl;

    TFile *fileIn = new TFile(datafile, "read");
    if (!fileIn || fileIn->IsZombie()) {
      cout << "ERROR: Cannot open " << datafile << endl;
      continue;
    }
    TTree *Tin = (TTree *)fileIn->Get("T1");
    if (!Tin) {
      cout << "ERROR: T1 tree not found in " << datafile << endl;
      continue;
    }

    Tin->SetBranchAddress("irun", &irun);
    Tin->SetBranchAddress("ievt", &ievt);
    Tin->SetBranchAddress("ngent", &ngent);
    Tin->SetBranchAddress("pidin", pidin);
    Tin->SetBranchAddress("ievt_wt", &ievt_wt);
    Tin->SetBranchAddress("momin", momin);
    Tin->SetBranchAddress("thein", thein);
    Tin->SetBranchAddress("phiin", phiin);
    Tin->SetBranchAddress("posxin", posxin);
    Tin->SetBranchAddress("posyin", posyin);
    Tin->SetBranchAddress("poszin", poszin);
    Tin->SetBranchAddress("nsimhtTk", &nsimhtTk);
    Tin->SetBranchAddress("detidTk", detidTk);
    Tin->SetBranchAddress("simpdgidTk", simpdgidTk);
    Tin->SetBranchAddress("simtimeTk", simtimeTk);
    Tin->SetBranchAddress("simenrTk", simenrTk);
    Tin->SetBranchAddress("simvxTk", simvxTk);
    Tin->SetBranchAddress("simvyTk", simvyTk);
    Tin->SetBranchAddress("simvzTk", simvzTk);
    Tin->SetBranchAddress("simpxTk", simpxTk);
    Tin->SetBranchAddress("simpyTk", simpyTk);
    Tin->SetBranchAddress("simpzTk", simpzTk);

    int nentries = Tin->GetEntries();
    cout << "Entries: " << nentries << endl;

    for (int iev = 0; iev < min(nentries, nentrymx); iev++) {
      fileIn->cd();
      Tin->GetEntry(iev);

      if (iev % 500 == 0)
        cout << "Event " << iev << "/" << nentries << endl;

      fileOut->cd();

      digihitpos.clear();
      hitIsSignal.clear();

      if (nsimhtTk == 0)
        continue;

      // Fill hit positions with signal/noise tagging
      fill_digipos_array(nsimhtTk, detidTk, simpdgidTk, simenrTk, simtimeTk,
                         digihitpos, hitIsSignal);

      if ((int)digihitpos.size() < nmnhit)
        continue;

      nhits_total = digihitpos.size();
      h_nhits->Fill(nhits_total);

      // Generator truth
      gen_theta = thein[0];
      gen_phi = phiin[0];
      gen_mom = momin[0] / 1000.0; // MeV → GeV
      gen_pt = gen_mom * sin(gen_theta);

      //============================================================
      // Per-event display plots (first N events)
      //============================================================
      if (global_evt < nEventDisplay) {

        c_evtdisp =
            new TCanvas(Form("c_evt%d", global_evt),
                        Form("Event %d Display", global_evt), 1200, 1000);
        c_evtdisp->Divide(2, 2);

        // ---- Pad 1: ρ vs Z (signal=black, noise=red) ----
        c_evtdisp->cd(1);
        gPad->SetGrid();

        vector<double> sig_r, sig_z, noi_r, noi_z;
        vector<double> sig_x, sig_y, noi_x, noi_y;

        for (int ih = 0; ih < (int)digihitpos.size(); ih++) {
          double r = digihitpos[ih].Perp();
          double z = digihitpos[ih].Z();
          double x = digihitpos[ih].X();
          double y = digihitpos[ih].Y();
          if (hitIsSignal[ih]) {
            sig_r.push_back(r);
            sig_z.push_back(z);
            sig_x.push_back(x);
            sig_y.push_back(y);
          } else {
            noi_r.push_back(r);
            noi_z.push_back(z);
            noi_x.push_back(x);
            noi_y.push_back(y);
          }
        }

        // Draw frame
        TH2F *frame_rz = new TH2F(Form("frz%d", global_evt), ";Rho(m);Z(m)", 10,
                                  0, 0.9, 10, -0.5, 0.5);
        frame_rz->Draw();

        if (!noi_r.empty()) {
          TGraph *gnoi = new TGraph(noi_r.size(), &noi_r[0], &noi_z[0]);
          gnoi->SetMarkerStyle(7);
          gnoi->SetMarkerColor(kRed);
          gnoi->Draw("P SAME");
        }
        if (!sig_r.empty()) {
          TGraph *gsig = new TGraph(sig_r.size(), &sig_r[0], &sig_z[0]);
          gsig->SetMarkerStyle(20);
          gsig->SetMarkerSize(0.8);
          gsig->SetMarkerColor(kBlack);
          gsig->Draw("P SAME");
        }

        TLegend *leg1 = new TLegend(0.12, 0.12, 0.35, 0.3);
        TGraph *gdummy_n = new TGraph();
        gdummy_n->SetMarkerStyle(7);
        gdummy_n->SetMarkerColor(kRed);
        TGraph *gdummy_s = new TGraph();
        gdummy_s->SetMarkerStyle(20);
        gdummy_s->SetMarkerColor(kBlack);
        leg1->AddEntry(gdummy_n, "Noise", "p");
        leg1->AddEntry(gdummy_s, "Signal", "p");
        leg1->Draw();

        // ---- Pad 2: X vs Y ----
        c_evtdisp->cd(2);
        gPad->SetGrid();
        TH2F *frame_xy = new TH2F(Form("fxy%d", global_evt), ";X(m);Y(m)", 10,
                                  0, 0.9, 10, -0.6, 0.6);
        frame_xy->Draw();

        if (!noi_x.empty()) {
          TGraph *gnoi2 = new TGraph(noi_x.size(), &noi_x[0], &noi_y[0]);
          gnoi2->SetMarkerStyle(7);
          gnoi2->SetMarkerColor(kRed);
          gnoi2->Draw("P SAME");
        }
        if (!sig_x.empty()) {
          TGraph *gsig2 = new TGraph(sig_x.size(), &sig_x[0], &sig_y[0]);
          gsig2->SetMarkerStyle(20);
          gsig2->SetMarkerSize(0.8);
          gsig2->SetMarkerColor(kBlack);
          gsig2->Draw("P SAME");
        }

        // ---- Pad 3: Hough Space R-Z (slope vs intercept) ----
        c_evtdisp->cd(3);
        gPad->SetGrid();

        // Generate all Hough pairs with classification
        vector<fourparam> all_hough;
        for (int i = 0; i < (int)digihitpos.size(); i++) {
          for (int j = i + 1; j < (int)digihitpos.size(); j++) {
            double r1 = digihitpos[i].Perp(), z1 = digihitpos[i].Z();
            double r2 = digihitpos[j].Perp(), z2 = digihitpos[j].Z();
            if (fabs(r2 - r1) < 1e-6)
              continue;
            fourparam fp;
            fp.hit1 = i;
            fp.hit2 = j;
            fp.slope = (z2 - z1) / (r2 - r1);
            fp.inter = z1 - fp.slope * r1;
            if (hitIsSignal[i] && hitIsSignal[j])
              fp.type = 0;
            else if (!hitIsSignal[i] && !hitIsSignal[j])
              fp.type = 2;
            else
              fp.type = 1;
            all_hough.push_back(fp);
          }
        }

        TH2F *frame_hough = new TH2F(Form("fh%d", global_evt),
                                     ";Slope (cot#theta);Intercept (m)", 10, -2,
                                     2, 10, -0.5, 0.5);
        frame_hough->Draw();

        // Draw noise (red), then mixed (black), then signal (green) on top
        vector<double> ns, ni, ms, mi, ss, si;
        for (auto &p : all_hough) {
          if (p.type == 2) {
            ns.push_back(p.slope);
            ni.push_back(p.inter);
          } else if (p.type == 1) {
            ms.push_back(p.slope);
            mi.push_back(p.inter);
          } else {
            ss.push_back(p.slope);
            si.push_back(p.inter);
          }
        }

        if (!ns.empty()) {
          TGraph *gn = new TGraph(ns.size(), &ns[0], &ni[0]);
          gn->SetMarkerStyle(7);
          gn->SetMarkerColor(kRed);
          gn->Draw("P SAME");
        }
        if (!ms.empty()) {
          TGraph *gm = new TGraph(ms.size(), &ms[0], &mi[0]);
          gm->SetMarkerStyle(7);
          gm->SetMarkerColor(kBlack);
          gm->Draw("P SAME");
        }
        if (!ss.empty()) {
          TGraph *gs = new TGraph(ss.size(), &ss[0], &si[0]);
          gs->SetMarkerStyle(20);
          gs->SetMarkerSize(0.5);
          gs->SetMarkerColor(kGreen + 2);
          gs->Draw("P SAME");
        }

        TLegend *leg3 = new TLegend(0.12, 0.7, 0.4, 0.88);
        TGraph *gd1 = new TGraph();
        gd1->SetMarkerStyle(7);
        gd1->SetMarkerColor(kRed);
        TGraph *gd2 = new TGraph();
        gd2->SetMarkerStyle(7);
        gd2->SetMarkerColor(kBlack);
        TGraph *gd3 = new TGraph();
        gd3->SetMarkerStyle(20);
        gd3->SetMarkerColor(kGreen + 2);
        leg3->AddEntry(gd1, "Noise", "p");
        leg3->AddEntry(gd2, "Mixed", "p");
        leg3->AddEntry(gd3, "Signal", "p");
        leg3->Draw();

        // ---- Pad 4: Zoomed Hough inset (signal cluster) ----
        c_evtdisp->cd(4);
        gPad->SetGrid();

        // Find signal cluster center for zoom
        double zoom_s = 0, zoom_i = 0;
        if (!ss.empty()) {
          for (int k = 0; k < (int)ss.size(); k++) {
            zoom_s += ss[k];
            zoom_i += si[k];
          }
          zoom_s /= ss.size();
          zoom_i /= si.size();
        }
        double zw_s = 0.3, zw_i = 0.1; // zoom window half-widths

        TH2F *frame_zoom = new TH2F(
            Form("fhz%d", global_evt), "Zoomed Signal Region;Slope;Intercept",
            10, zoom_s - zw_s, zoom_s + zw_s, 10, zoom_i - zw_i, zoom_i + zw_i);
        frame_zoom->Draw();

        // Replot only points in zoom window
        for (auto &p : all_hough) {
          if (fabs(p.slope - zoom_s) < zw_s && fabs(p.inter - zoom_i) < zw_i) {
            TGraph *gpt = new TGraph(1);
            gpt->SetPoint(0, p.slope, p.inter);
            if (p.type == 0) {
              gpt->SetMarkerStyle(20);
              gpt->SetMarkerSize(0.6);
              gpt->SetMarkerColor(kGreen + 2);
            } else if (p.type == 1) {
              gpt->SetMarkerStyle(7);
              gpt->SetMarkerColor(kBlack);
            } else {
              gpt->SetMarkerStyle(7);
              gpt->SetMarkerColor(kRed);
            }
            gpt->Draw("P SAME");
          }
        }

        // Save
        if (!pdfOpened) {
          c_evtdisp->SaveAs("track_event_display.pdf(");
          pdfOpened = true;
        } else {
          c_evtdisp->SaveAs("track_event_display.pdf");
        }
      }

      //============================================================
      // TRACK RECONSTRUCTION (single track)
      //============================================================
      vector<fourparam> cleaned_rz;
      vector<int> selected_hits_rz;
      reco_theta =
          hough_rz(digihitpos, hitIsSignal, cleaned_rz, selected_hits_rz);

      double rphi_phi, rphi_pt;
      hough_rphi(digihitpos, hitIsSignal, selected_hits_rz, rphi_phi, rphi_pt);
      reco_phi = rphi_phi;
      reco_pt = rphi_pt;

      nhits_used = selected_hits_rz.size();

      //============================================================
      // MULTIPLE TRAJECTORY FINDING
      //============================================================
      ntrack_found = 1;
      vector<TVector3> remaining_hits = digihitpos;
      vector<bool> remaining_signal = hitIsSignal;

      // Remove track 1 hits
      set<int> used_set(selected_hits_rz.begin(), selected_hits_rz.end());
      vector<TVector3> leftover_hits;
      vector<bool> leftover_signal;
      for (int k = 0; k < (int)remaining_hits.size(); k++) {
        if (!used_set.count(k)) {
          leftover_hits.push_back(remaining_hits[k]);
          leftover_signal.push_back(remaining_signal[k]);
        }
      }

      // Try to find more tracks in remaining hits
      while ((int)leftover_hits.size() >= nmnhit) {
        vector<fourparam> cl;
        vector<int> sel;
        double tht = hough_rz(leftover_hits, leftover_signal, cl, sel);
        if ((int)sel.size() < nmnhit)
          break;
        ntrack_found++;

        // Remove these hits
        set<int> used2(sel.begin(), sel.end());
        vector<TVector3> tmp_hits;
        vector<bool> tmp_signal;
        for (int k = 0; k < (int)leftover_hits.size(); k++) {
          if (!used2.count(k)) {
            tmp_hits.push_back(leftover_hits[k]);
            tmp_signal.push_back(leftover_signal[k]);
          }
        }
        leftover_hits = tmp_hits;
        leftover_signal = tmp_signal;
      }
      h_ntracks->Fill(ntrack_found);

      //============================================================
      // Resolution (first track only, matching slide units)
      //============================================================
      delta_theta_mrad = (reco_theta - gen_theta) * 1000.0; // rad → mrad
      delta_phi_mrad = sfdph(gen_phi, reco_phi) * 1000.0;   // rad → mrad
      delta_pt_GeV = reco_pt - gen_pt;                      // GeV

      h_delta_theta->Fill(delta_theta_mrad);
      h_delta_phi->Fill(delta_phi_mrad);
      h_delta_pt->Fill(delta_pt_GeV);
      h_reco_pt_dist->Fill(reco_pt);

      Tout->Fill();
      global_evt++;

    } // events
    fileIn->cd();
    delete Tin;
    delete fileIn;
  } // files

  //============================================================
  // Close event display PDF
  //============================================================
  if (pdfOpened) {
    // Create a blank closing page
    TCanvas *cclose = new TCanvas("cclose", "", 100, 100);
    cclose->SaveAs("track_event_display.pdf)");
  }

  //============================================================
  // Resolution plots (matching slide 3 format)
  //============================================================
  gStyle->SetOptStat(1111);
  gStyle->SetStatX(0.88);
  gStyle->SetStatY(0.88);

  TCanvas *c_res = new TCanvas("c_res", "Resolution", 1500, 500);
  c_res->Divide(3, 1);

  c_res->cd(1);
  gPad->SetLogy();
  h_delta_pt->Draw();
  TLatex *t1 = new TLatex();
  t1->SetNDC();
  t1->SetTextSize(0.04);
  t1->SetTextColor(kGreen + 2);

  c_res->cd(2);
  gPad->SetLogy();
  h_delta_phi->Draw();

  c_res->cd(3);
  gPad->SetLogy();
  h_delta_theta->Draw();

  c_res->SaveAs("track_resolution.pdf(");

  // Additional summary plots
  TCanvas *c_sum = new TCanvas("c_sum", "Summary", 1200, 400);
  c_sum->Divide(3, 1);

  c_sum->cd(1);
  h_nhits->SetLineWidth(2);
  h_nhits->Draw();
  c_sum->cd(2);
  h_reco_pt_dist->SetLineWidth(2);
  h_reco_pt_dist->Draw();
  c_sum->cd(3);
  h_ntracks->SetLineWidth(2);
  h_ntracks->Draw();

  c_sum->SaveAs("track_resolution.pdf)");

  //============================================================
  // Write output
  //============================================================
  fileOut->cd();
  fileOut->Write();

  // Print summary
  cout << "\n========== TRACK RECONSTRUCTION SUMMARY ==========" << endl;
  cout << "Events processed: " << global_evt << endl;
  cout << "Noise hits per event: " << nNoiseHit << endl;
  cout << "\nResolution:" << endl;
  cout << "  Delta pT:    Mean=" << h_delta_pt->GetMean()
       << " GeV, StdDev=" << h_delta_pt->GetStdDev() << " GeV" << endl;
  cout << "  Delta phi:   Mean=" << h_delta_phi->GetMean()
       << " mrad, StdDev=" << h_delta_phi->GetStdDev() << " mrad" << endl;
  cout << "  Delta theta: Mean=" << h_delta_theta->GetMean()
       << " mrad, StdDev=" << h_delta_theta->GetStdDev() << " mrad" << endl;
  cout << "\nOutput files:" << endl;
  cout << "  output_skeleton.root" << endl;
  cout << "  track_event_display.pdf" << endl;
  cout << "  track_resolution.pdf" << endl;
  cout << "===================================================" << endl;

  return 0;
}

////////////////////////////////////////////////////////////////////////////////////
// Digitisation: decode detid → position, add noise hits
// Now also tags each hit as signal or noise
////////////////////////////////////////////////////////////////////////////////////
void fill_digipos_array(int nhit, unsigned int *detid, int *pdgid,
                        float *simenr, float *simtime, vector<TVector3> &hitpos,
                        vector<bool> &isSignal) {
  TVector3 tmp3vect;

  for (int ij = 0; ij < nhit; ij++) {
    if (simenr[ij] <= 0)
      continue;

#ifdef DELTARAY
    if (abs(pdgid[ij]) != 13)
      continue;
    if (simtime[ij] > 10)
      continue;
#endif

    if (gRandom3->Uniform() > 0.90)
      continue;
    int irad = ((detid[ij] >> 28) & 0xF);
    if (irad < 0 || irad >= nsilayer)
      continue;
    if (irad <= 6 && gRandom3->Uniform() > 0.90)
      continue;

    double rho = sirad[irad];
    int InThe = ((detid[ij] >> 15) & 0x1FFF);
    double zz = (0.2 * (InThe + 0.5) - 750.0);
#define SMEARING
#ifdef SMEARING
    zz += gRandom3->Gaus(0.0, 0.2);
#endif
    zz = (zz / 1000.0);

    int InPhi = (detid[ij] & 0x7FFF);
    double phi = ((InPhi + 0.5) / 20000.0 - pibyfour);
#ifdef SMEARING
    phi += gRandom3->Gaus(0.0, 0.05);
    double xx = gRandom3->Gaus(0.0, 1.e-4);
    phi += xx / rho;
#endif
    tmp3vect = MakeRhoPhiZ(rho, phi, zz);
    hitpos.push_back(tmp3vect);
    isSignal.push_back(true); // This is a signal hit
  }

  // Noise hits (configurable via nNoiseHit constant)
  for (int ij = 0; ij < nNoiseHit; ij++) {
    double phi =
        (int((pibytwo * gRandom3->Uniform() - pibyfour) * 20000 + 0.5)) /
        20000.;
    double the = acos(sqrt(2) * (gRandom3->Uniform() - 0.5));
    double rad = sirad[int(nsilayer * gRandom3->Uniform())];
    tmp3vect = MakeRhoPhiTheta(rad, phi, the);
    hitpos.push_back(tmp3vect);
    isSignal.push_back(false); // This is a noise hit
  }
}

////////////////////////////////////////////////////////
// Straight line fit: z = slope*r + intercept (Luis Lyons)
///////////////////////////////////////////////////////
void straight_line_fit(vector<TVector3> finder, double *par, double *parerr) {

  double zerr2 = 4.e-8;
  const int size = finder.size();
  double srhoz = 0, srho = 0, sz = 0, sn = 0, srho2 = 0;

  for (int ij = 0; ij < size; ij++) {
    srhoz += finder[ij].Perp() * finder[ij].Z() / zerr2;
    srho += finder[ij].Perp() / zerr2;
    srho2 += finder[ij].Perp() * finder[ij].Perp() / zerr2;
    sz += finder[ij].Z() / zerr2;
    sn += 1 / zerr2;
  }

  if (sn > 0 && srho2 * sn - srho * srho != 0) {
    double slope = par[0] =
        (srhoz * sn - srho * sz) / (srho2 * sn - srho * srho);
    par[1] = sz / sn - slope * srho / sn;

    double determ = (sn * srho2 - srho * srho);
    parerr[0] = srho2 / determ;
    parerr[1] = -srho / determ;
    parerr[2] = sn / determ;
  } else {
    par[0] = par[1] = parerr[0] = parerr[1] = parerr[2] = -1;
  }
}

////////////////////////////////////////////////////////////
// Track extrapolation with energy loss
////////////////////////////////////////////////////////////
bool GetPrediction(double radin, double *Statev, double *Prediction, int il,
                   double *Q_k_minus, int iiter) {

  double charge = (Statev[4] > 0) ? 1 : -1;
  double aftin = -charge / momconv;
  double pt0 = abs(1. / Statev[4]);
  double thetx = (abs(Statev[2]) > 1.e-13) ? atan(1. / Statev[2]) : pibytwo;
  if (thetx < 0)
    thetx += pival;

  double ptot = pt0 / sin(thetx);
  double px0 = pt0 * cos(Statev[3]);
  double py0 = pt0 * sin(Statev[3]);
  double pz0 = ptot * cos(thetx);
  double x0 = radin * cos(Statev[0]);
  double y0 = radin * sin(Statev[0]);
  double z0 = Statev[1];

  double aa = sirad[il] * sirad[il] - radin * radin;
  double bb = 2 * aftin * (aftin * pt0 * pt0 - x0 * py0 + y0 * px0);
  double cc = 2 * aftin * (x0 * px0 + y0 * py0);

  double axa = bb * bb + cc * cc;
  double bxb = 2 * bb * (aa - bb);
  double cxc = ((aa - bb) * (aa - bb) - cc * cc);
  double det = bxb * bxb - 4 * axa * cxc;

  if (det < 0) {
    for (int ix = 0; ix < 5; ix++)
      Prediction[ix] = Statev[4];
    return false;
  }

  double cosrs =
      max(-0.999999999999999,
          min(0.99999999999999, (-bxb + sqrt(max(1.e-50, det))) / (2 * axa)));
  double sinrs = -charge * sqrt(1.0 - cosrs * cosrs);

  double x = x0 + aftin * (px0 * sinrs - py0 * (1 - cosrs));
  double y = y0 + aftin * (py0 * sinrs + px0 * (1 - cosrs));
  double ss = acos(cosrs) * ptot * abs(aftin);
  double z = z0 + ss * cos(thetx);
  double px = px0 * cosrs - py0 * sinrs;
  double py = py0 * cosrs + px0 * sinrs;
  double pz = pz0;

  double elos = 0.001826 * (ss / (sirad[il] - radin)) *
                (14.9 + 0.96 * fabs(log(ptot * 2.8)) +
                 0.033 * ptot * (1.0 - pow(ptot, -0.33))) *
                1e-2 * 1.10;

  double pscale = (ptot - elos) / ptot;
  ptot *= pscale;
  px *= pscale;
  py *= pscale;
  pz *= pscale;

  Prediction[0] = atan2(y, x);
  Prediction[1] = z;
  double thet = acos(pz / ptot);
  if (thet < 0)
    thet += pival;
  Prediction[2] = 1 / tan(thet);
  Prediction[3] = atan2(py, px);
  Prediction[4] = charge / (ptot * sin(thet));
  return true;
}

// Track propagation
void track_move_pt_align(double *b, double *tin, double q, double *x,
                         double *tout, double &dl) {
  double c_b = 1.0;
  double brich = pow(b[0] * b[0] + b[1] * b[1] + b[2] * b[2], 0.5);
  if (brich < 0.00000001)
    brich = 0.00001;
  double pt = pow(pow(tin[3], 2.) + pow(tin[4], 2.), 0.5);
  if (pt != 0.0 && q != 0.0) {
    double a = -c_b * brich * q;
    double ptinv = 1. / pt;
    double rho = a * ptinv;
    double delx = tin[0] - x[0], dely = tin[1] - x[1];
    double ainv = 1. / a;
    double px = tin[3], py = tin[4];
    double cosps = 1. - rho * (delx * py - dely * px) * ptinv;
    double sinps = -rho * (delx * px + dely * py) * ptinv;
    double sovp = atan2(sinps, cosps) * ainv;
    double alpha = 1. / sqrt(cosps * cosps + sinps * sinps);
    cosps *= alpha;
    sinps *= alpha;
    double dcosps = 1. - cosps;
    tout[3] = px * cosps - py * sinps;
    tout[4] = py * cosps + px * sinps;
    tout[5] = tin[5];
    tout[0] = tin[0] + (px * sinps - py * dcosps) * ainv;
    tout[1] = tin[1] + (py * sinps + px * dcosps) * ainv;
    tout[2] = tin[2] + abs(sovp) * tin[5];
    dl = fabs(sovp) *
         pow(pow(tin[3], 2.0) + pow(tin[4], 2.) + pow(tin[5], 2.), 0.5);
  } else {
    for (int i = 0; i < 6; i++)
      tout[i] = tin[i];
    dl = 0.;
  }
}
