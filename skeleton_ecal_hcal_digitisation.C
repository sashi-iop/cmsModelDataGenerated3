/*
  ECAL + HCAL Digitisation Macro

  This macro reads simulation output (T1 tree from Geant4),
  applies realistic detector digitisation:
    - ECAL: crystal-level noise (pedwidth) + energy thresholds
    - HCAL: unpacks 32-bit words, applies tower noise + thresholds
  and writes digitised energies to a new ROOT file (T2 tree).

  Usage:
    ./ecal_hcal_digitisation <ecal_noise_MeV> <ecal_threshold_MeV> <hcal_noise_MeV> <hcal_threshold_MeV>

  Example:
    ./ecal_hcal_digitisation 40 200 25 100

  Input: reads file list from "test_pion_klong.log"
         Format per line: <rootfile.root> <max_events>
  Output: creates a new ROOT file with digitised T2 tree

  PbWO4 : 0.7%
  Plastic : ~40% Factor 50, but loss in fibre, only 7% ->3.5%  (overall 5
  times)->10 times
 */
// Note: CLHEP Matrix/LorentzVector not needed for digitisation
// #include "CLHEP/Matrix/Matrix.h"
// #include "CLHEP/Vector/LorentzVector.h"
// #include "CLHEP/Vector/ThreeVector.h"
#include <cmath>
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TObject.h"
#include "TPostScript.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TTree.h"
#include <TRandom3.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#define SMEARING
using namespace std;

const int nhclwedge = 36;
const int nhcalLayer = 7; // Merged depth: 17 layers -> 3 bits -> max 8 depth bins (0-7), but effective ~6 used
const int nhcalEtaDiv = 34;

double thetval[34] = {
    0.652481, 0.684299, 0.718216, 0.755979, 0.791376, 0.834511, 0.879688,
    0.930175, 0.983599, 1.04162,  1.10367,  1.17053,  1.24089,  1.31505,
    1.3919,   1.47074,  1.54993,  1.6316,   1.71117,  1.78886,  1.86516,
    1.93752,  2.00573,  2.06932,  2.12961,  2.18582,  2.23783,  2.28572,
    2.32959,  2.366,    2.40572,  2.44111,  2.47378,  2.50406};

const double pival = acos(-1.0);
const double twopi = 2 * pival;
const double pibytwo = pival / 2.;
const double pibyfour = pival / 4.;
const double degtorad = pival / 180.;
const double radtodeg = 180. / pival;
static unsigned int mypow_2[32];

double sfdph(double phi1, double phi2) {
  double dphi = phi2 - phi1;
  if (dphi > M_PI) {
    dphi = 2 * M_PI - dphi;
  }
  if (dphi < -M_PI) {
    dphi = 2 * M_PI + dphi;
  }
  return dphi;
}

// Helper function (requires CLHEP Hep3Vector, uncomment if needed)
// double dgap(Hep3Vector v1, Hep3Vector v2) {
//   double dthe = v1.theta() - v2.theta();
//   double dphi = sfdph(v1.phi(), v2.phi());
//   return sqrt(dthe * dthe + dphi * dphi);
// }

TRandom3 *gRandom3 = new TRandom3();

int main(int argc, char **argv) {

  if (argc < 5) {
    cout << "Usage: " << argv[0] << " <ecal_noise_MeV> <ecal_threshold_MeV> <hcal_noise_MeV> <hcal_threshold_MeV>" << endl;
    cout << "Give seedthresh, pedwidth and otherthrs in MeV unit" << endl;
    return 0;
  }

  const double pedwidth =
      atof(argv[1]) / 1000.0; // ECAL Noise per crystal in GeV (input in MeV)
  const double otherthrs =
      atof(argv[2]) / 1000.0; // Energy threshold for ECAL crystal in GeV
  const double hclwidth =
      atof(argv[3]) / 1000.0; // HCAL noise per tower in GeV (input in MeV)
  const double hclthrs =
      atof(argv[4]) / 1000.0; // Threshold for HCAL tower energy in GeV

  cout << "ECAL noise=" << pedwidth << " GeV, ECAL threshold=" << otherthrs
       << " GeV, HCAL noise=" << hclwidth << " GeV, HCAL threshold=" << hclthrs
       << " GeV" << endl;

  for (int ij = 0; ij < 32; ij++) {
    mypow_2[ij] = (int)pow(float(2), ij);
  }

  gStyle->SetOptStat(111111);
  char datafile[100];

  unsigned irun; // Run number of these events
  unsigned ievt; // Event number
  unsigned ngent;

  const unsigned int ngenmx = 50;
  int pidin[ngenmx];   // PID of incident particle
  float momin[ngenmx]; // Energy of incident particle
  float thein[ngenmx]; // Initial polar angle of incident particle
  float phiin[ngenmx]; // Initial azimuthal angle of incident particle

  // For Simulation output of ECAL
  unsigned int nsimhtEC;
  static const unsigned int nsimhtmxEC = 2000;
  unsigned int detidEC[nsimhtmxEC];
  float energyEC[nsimhtmxEC];
  float thetaEC[nsimhtmxEC];
  float phiEC[nsimhtmxEC];

  unsigned int nsimhtTk; // Looked only for early showering

  // For Simulation output of HCAL
  // NOTE: detidHL is unsigned long int (64-bit) in the simulation TTree branch ("/l")
  unsigned int nsimhtHL;
  static const unsigned int nsimhtmxHL = 2000;
  unsigned long int detidHL[nsimhtmxHL];  // Must match the "/l" branch format
  unsigned int timeHL[nsimhtmxHL];

  char rootfiles[100];
  char outfile[100];
  char outfilx[100];
  ifstream file_db;

  sprintf(rootfiles, "test_pion_klong.log");

  int len = strlen(rootfiles);
  strncpy(outfilx, rootfiles, len - 4);
  outfilx[len - 4] = '\0';
  sprintf(outfile, "%s_xwoecl_%s_%s_%s_%s.root", outfilx, argv[1], argv[2],
          argv[3], argv[4]);

  float totsimenr = 0;   // Total simulated energy (ECAL + HCAL) in GeV
  float totdigienr = 0;  // Total digitised energy (ECAL + HCAL) in GeV
  float hclsimenr = 0;   // HCAL simulated energy in GeV
  float hcldigienr = 0;  // HCAL digitised energy in GeV

  TFile *fileOut = new TFile(outfile, "recreate");
  TTree *Tout = new TTree("T2", "Digitised ECAL+HCAL");

  Tout->Branch("totsimenr", &totsimenr, "totsimenr/F");
  Tout->Branch("totdigienr", &totdigienr, "totdigienr/F");
  Tout->Branch("hclsimenr", &hclsimenr, "hclsimenr/F");
  Tout->Branch("hcldigienr", &hcldigienr, "hcldigienr/F");

  file_db.open(rootfiles);
  while (!(file_db.eof())) {
    int nentrymx = -1;

    file_db >> datafile >> nentrymx;
    if (strstr(datafile, "#"))
      continue;

    cout << "datafile = " << datafile << endl;

    if (file_db.eof())
      break;

    TFile *fileIn = new TFile(datafile, "read");
    TTree *Tin = (TTree *)fileIn->Get("T1");

    Tin->SetBranchAddress("irun", &irun);
    Tin->SetBranchAddress("ievt", &ievt);

    Tin->SetBranchAddress("ngent", &ngent);
    Tin->SetBranchAddress("pidin", pidin);
    Tin->SetBranchAddress("momin", momin);
    Tin->SetBranchAddress("thein", thein);
    Tin->SetBranchAddress("phiin", phiin);

    // Simulation output of ECAL
    Tin->SetBranchAddress("nsimhtEC", &nsimhtEC);
    Tin->SetBranchAddress("detidEC", detidEC);
    Tin->SetBranchAddress("energyEC", energyEC);
    Tin->SetBranchAddress("thetaEC", thetaEC);
    Tin->SetBranchAddress("phiEC", phiEC);

    // Simulation output of HCAL
    Tin->SetBranchAddress("nsimhtHL", &nsimhtHL);
    Tin->SetBranchAddress("detidHL", detidHL);
    Tin->SetBranchAddress("timeHL", timeHL);

    int nentries = Tin->GetEntries();
    cout << "nentr " << datafile << " " << nentries << endl;

    for (int iev = 0; iev < min(nentries, nentrymx); iev++) {
      fileIn->cd();
      Tin->GetEntry(iev);
      if (nsimhtEC == 0 && nsimhtHL == 0)
        continue;

      totsimenr = 0;
      totdigienr = 0;
      hclsimenr = 0;
      hcldigienr = 0;

      //============================================================
      // ECAL Digitisation
      // - energyEC[] is in native Geant4 units (MeV)
      // - Add Gaussian noise (pedwidth) per crystal
      // - Apply energy threshold (otherthrs) per crystal
      // - Sum surviving crystals for total ECAL digitised energy
      //============================================================
      float ecal_sim_total = 0;
      float ecal_digi_total = 0;

      for (unsigned int ih = 0; ih < nsimhtEC; ih++) {
        float sim_energy_GeV = energyEC[ih] / 1000.0;  // MeV -> GeV
        ecal_sim_total += sim_energy_GeV;

#ifdef SMEARING
        // Add Gaussian noise to each crystal
        float smeared_energy = sim_energy_GeV + gRandom3->Gaus(0, pedwidth);
#else
        float smeared_energy = sim_energy_GeV;
#endif
        // Apply threshold: only include crystals above threshold
        if (smeared_energy > otherthrs) {
          ecal_digi_total += smeared_energy;
        }
      }

      //============================================================
      // HCAL Digitisation
      // - detidHL[] is a packed 32-bit word:
      //     bits 17-31: cell ID (ieta[6] + iphi[6] + depth[3])
      //     bits 0-16:  energy in MeV (17 bits, least count = 1 MeV)
      // - Extract energy from lower 17 bits
      // - Add Gaussian noise (hclwidth) per tower
      // - Apply energy threshold (hclthrs) per tower
      // - Sum surviving towers for total HCAL digitised energy
      //============================================================
      float hcal_sim_total = 0;
      float hcal_digi_total = 0;

      for (unsigned int ih = 0; ih < nsimhtHL; ih++) {
        // Extract energy from lower 17 bits of the packed word
        unsigned int energy_MeV = detidHL[ih] & 0x1FFFF;  // 2^17 - 1 = 131071
        float sim_energy_GeV = energy_MeV / 1000.0;       // MeV -> GeV

        // Extract cell ID from upper bits (for debugging/analysis)
        // unsigned int cellid = (detidHL[ih] >> 17) & 0x7FFF;  // 15-bit cell ID
        // unsigned int ieta   = (cellid >> 9) & 0x3F;
        // unsigned int iphi   = (cellid >> 3) & 0x3F;
        // unsigned int depth  = cellid & 0x7;

        hcal_sim_total += sim_energy_GeV;

#ifdef SMEARING
        // Add Gaussian noise to each tower
        float smeared_energy = sim_energy_GeV + gRandom3->Gaus(0, hclwidth);
#else
        float smeared_energy = sim_energy_GeV;
#endif
        // Apply threshold: only include towers above threshold
        if (smeared_energy > hclthrs) {
          hcal_digi_total += smeared_energy;
        }
      }

      //============================================================
      // Store results
      //============================================================
      // ECAL + HCAL combined
      totsimenr  = ecal_sim_total + hcal_sim_total;
      totdigienr = ecal_digi_total + hcal_digi_total;

      // HCAL only (for calibration without ECAL)
      hclsimenr  = hcal_sim_total;
      hcldigienr = hcal_digi_total;

      fileOut->cd();
      Tout->Fill();
    } // end of events loop

    fileIn->cd();
    delete Tin;
    delete fileIn;

  } // end of file loop
  fileOut->cd();
  fileOut->Write();

  cout << "Digitisation complete. Output file: " << outfile << endl;
}
