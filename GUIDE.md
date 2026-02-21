# CMS Model — Parameter Tuning & ROOT Analysis Guide

> A self-contained reference for modifying simulation parameters and analyzing ROOT output.

---

## Table of Contents

- [Part A: Simulation Parameters](#part-a-simulation-parameters)
  - [1. Changing Particle Type](#1-changing-particle-type)
  - [2. Changing Beam Energy](#2-changing-beam-energy)
  - [3. Changing Number of Events](#3-changing-number-of-events)
  - [4. Changing Beam Direction & Position](#4-changing-beam-direction--position)
  - [5. Changing Energy Smearing](#5-changing-energy-smearing)
  - [6. Changing Angle Smearing](#6-changing-angle-smearing)
  - [7. Changing Vertex Smearing](#7-changing-vertex-smearing)
  - [8. Multiple Particles per Event](#8-multiple-particles-per-event)
  - [9. Random Mode On/Off](#9-random-mode-onoff)
  - [10. Output File Name](#10-output-file-name)
  - [11. Digitisation Parameters](#11-digitisation-parameters)
  - [12. Detector Physics Parameters (Source Code)](#12-detector-physics-parameters-source-code)
  - [13. Quick Reference Table](#13-quick-reference-table)
- [Part B: ROOT File Analysis](#part-b-root-file-analysis)
  - [Level 1–13: Interactive ROOT analysis](#level-1-opening--browsing-root-files)
  - [Available Branches & Histograms](#available-branches--histograms)
- [Part C: Track Finder](#part-c-track-finder)
  - [1. Overview](#c1-overview)
  - [2. How to Run](#c2-how-to-run)
  - [3. Hough Transform Theory](#c3-hough-transform-theory)
  - [4. Configurable Parameters](#c4-configurable-parameters)
  - [5. Noise Comparison](#c5-noise-comparison)
  - [6. Output Files](#c6-output-files)
  - [7. Understanding Resolution Plots](#c7-understanding-resolution-plots)
  - [8. Multi-Track Finding](#c8-multi-track-finding)
  - [9. T2 Tree Branches (Track Finder)](#c9-t2-tree-branches-track-finder)
  - [10. Troubleshooting](#c10-troubleshooting)

---

# Part A: Simulation Parameters

There are **two ways** to change parameters:
- **Method 1 — `run.mac` file** (recommended, no recompilation)
- **Method 2 — Source code** (requires `make` rebuild)

---

## 1. Changing Particle Type

### Via `run.mac` (no rebuild):
```
# Particle ID (PDG code)
/Serc19/gun/pid  211
```

### Via source code (`src/Serc19PrimaryGeneratorAction.cc`, line ~49):
```cpp
SetPartId(211);    // Set default PID
```

### Common PDG codes:
| PID | Particle | Use Case |
|-----|----------|----------|
| `11` | e⁻ (electron) | EM shower calibration |
| `-11` | e⁺ (positron) | EM shower calibration |
| `13` | μ⁻ (muon) | MIP calibration |
| `22` | γ (photon) | EM shower |
| `111` | π⁰ | Decays to 2γ |
| `211` | π⁺ | Hadronic shower (default) |
| `-211` | π⁻ | Hadronic shower |
| `321` | K⁺ | Hadronic comparison |
| `-321` | K⁻ | Hadronic comparison |
| `2212` | proton | Nuclear interactions |
| `2112` | neutron | Missing energy studies |

### To enable random PID mixing:
```
/Serc19/gun/rndmPID  on
```

---

## 2. Changing Beam Energy

### Via `run.mac` (no rebuild):
```
# Energy in GeV
/Serc19/gun/energy 200.0
```

### Via source code (`src/Serc19PrimaryGeneratorAction.cc`, line ~51):
```cpp
SetIncEnergy(200.0*GeV);    // Energy value
```

### Typical energy scan values:
```
# For calibration studies, try these one at a time:
/Serc19/gun/energy 10.0
/Serc19/gun/energy 20.0
/Serc19/gun/energy 50.0
/Serc19/gun/energy 100.0
/Serc19/gun/energy 200.0
```

---

## 3. Changing Number of Events

### Via `run.mac` (last line):
```
# Number of events to simulate
/run/beamOn 1000
```

### Typical values:
| Events | Purpose | Approx Time |
|--------|---------|-------------|
| 10 | Quick test | seconds |
| 100 | Development | ~1 min |
| 1000 | Decent statistics | ~10 min |
| 10000 | Publication quality | ~1-2 hrs |

---

## 4. Changing Beam Direction & Position

### Direction (unit vector components):
```
# Format: /Serc19/gun/incdir  x  y  z
/Serc19/gun/incdir 1.0 0 0      # Along +X axis (default)
/Serc19/gun/incdir 1.0 0 0.6    # Tilted in XZ plane
/Serc19/gun/incdir 0 0 1        # Along +Z axis
```

### Starting position:
```
# Format: /Serc19/gun/incpos  x  y  z (in cm)
/Serc19/gun/incpos  0.  0.  0.    # Origin (default)
/Serc19/gun/incpos  0.  5.  0.    # Shifted 5cm in Y
```

---

## 5. Changing Energy Smearing

```
# Sigma of Gaussian energy smearing in GeV
# Positive = fixed sigma, Negative = uniform random between 0 and |value|
/Serc19/gun/ensmear 0.5       # Gaussian sigma = 0.5 GeV
/Serc19/gun/ensmear -30.0     # Uniform E in [0, 30] GeV (negative = flat)
```

### In source code (`src/Serc19PrimaryGeneratorAction.cc`, line ~52):
```cpp
SetIncEnergySmr(100*MeV);    // Smearing sigma
```

---

## 6. Changing Angle Smearing

```
# Theta smearing in cos(theta) terms
# Negative value = enable uniform random direction
/Serc19/gun/thsmear -700.     # Random in cos(theta)

# Phi smearing in milliradians
/Serc19/gun/phsmear -500.     # Random phi spread
```

---

## 7. Changing Vertex Smearing

```
# Gaussian spread of starting position (in cm)
/Serc19/gun/vxsmear 0.01    # X spread = 0.01 cm
/Serc19/gun/vysmear 0.01    # Y spread = 0.01 cm
/Serc19/gun/vzsmear 5.0     # Z spread = 5.0 cm
```

---

## 8. Multiple Particles per Event

```
# Fire N particles per event
/Serc19/gun/mult 1      # Single particle (default)
/Serc19/gun/mult 5      # 5 particles per event
```

---

## 9. Random Mode On/Off

```
# Enable/disable random generation
/Serc19/gun/rndm  on     # Random energy/direction each event
/Serc19/gun/rndm  off    # Fixed energy/direction every event
```

---

## 10. Output File Name

```
# Base name for output (ROOT file will be <name>_run10.root)
/Serc19/run/output_file   test
/Serc19/run/output_file2  test_geant.root
```

---

## 11. Digitisation Parameters

These are passed as command-line arguments (in **MeV**):
```bash
./ecal_hcal_digitisation <ecal_noise> <ecal_threshold> <hcal_noise> <hcal_threshold>
```

| Parameter | Description | Typical Value | Effect |
|-----------|-------------|---------------|--------|
| Arg 1 | ECAL noise/crystal | 40 MeV | Electronic noise per PbWO₄ crystal |
| Arg 2 | ECAL threshold | 200 MeV | Only crystals above this contribute |
| Arg 3 | HCAL noise/tower | 25 MeV | Electronic noise per HCAL tower |
| Arg 4 | HCAL threshold | 100 MeV | Only towers above this contribute |

### Example parameter scans:
```bash
# Low noise
./ecal_hcal_digitisation 20 100 10 50

# High noise (stress test)
./ecal_hcal_digitisation 100 300 80 200

# No noise (perfect detector)
# Comment out #define SMEARING in the source code, then rebuild
```

---

## 12. Detector Physics Parameters (Source Code)

These require editing source files and rebuilding (`make`):

### Energy Resolution — `src/Serc19HclSD.cc`

```cpp
// Photo-electron yield (p.e. per MeV) — line ~120
G4double npe_mean = 40.0 * edep_MeV;     // Change 40.0 to adjust
//   Higher → better stochastic resolution
//   Lower  → worse stochastic resolution

// Collection efficiency — line ~130
G4double efficiency = 1.0 + 0.1 * G4double(approx_depth) / G4double(nhcalLayer);
//   Change 0.1 to adjust depth dependence
//   0.0 → uniform, 0.2 → strong depth dependence

// Electronic noise sigma (MeV) — line ~137
edep_smeared_MeV += G4RandGauss::shoot(0.0, 200.0);
//   Change 200.0 to adjust noise level
//   Smaller → cleaner signal
//   Larger  → more noise
```

### 32-bit Packing — `src/Serc19HclSD.cc`

```cpp
// Energy least count — currently 1 MeV (17 bits, max 131071 MeV)
// To change to 0.1 MeV: multiply by 10 before packing, divide when reading
G4int energy_units = G4int(edep_smeared_MeV * 10);  // 0.1 MeV least count
// Then unpack: energy_MeV = (detidHL & 0x1FFFF) / 10.0;
```

---

## 13. Quick Reference Table

| Parameter | Where | Value |
|-----------|-------|-------|
| Particle type | `run.mac`: `/Serc19/gun/pid` | 211, 321, 13, etc. |
| Energy | `run.mac`: `/Serc19/gun/energy` | in GeV |
| N events | `run.mac`: `/run/beamOn` | integer |
| Direction | `run.mac`: `/Serc19/gun/incdir` | x y z |
| Energy smear | `run.mac`: `/Serc19/gun/ensmear` | GeV |
| Photo-electrons | `Serc19HclSD.cc` line ~120 | p.e./MeV |
| Noise sigma | `Serc19HclSD.cc` line ~137 | MeV |
| ECAL digi noise | CLI arg 1 | MeV |
| HCAL digi noise | CLI arg 3 | MeV |

---

# Part B: ROOT File Analysis

All examples should be run from the `build/` directory.

---

## Level 1: Opening & Browsing ROOT Files

### Interactive ROOT session:
```cpp
// Start ROOT
root -l

// Open a file
TFile *f = TFile::Open("test_run10.root");

// See what's inside
f->ls();                    // List all objects
f->Print();                 // Detailed info

// Use the GUI browser
TBrowser b;                 // Opens graphical browser (double-click histograms!)
```

### One-liner from shell:
```bash
# List contents without opening ROOT interactively
root -l -b -q -e 'TFile f("test_run10.root"); f.ls();'
```

---

## Level 2: Reading TTree Branches

```cpp
root -l

// Open file and get tree
TFile *f = TFile::Open("test_run10.root");
TTree *T1 = (TTree*)f->Get("T1");

// See all branches
T1->Print();                 // Shows branch names, types, sizes

// Quick scan of first 10 events
T1->Scan("irun:ievt:ngent:nsimhtEC:nsimhtHL", "", "", 10);

// Get number of entries
cout << "Total events: " << T1->GetEntries() << endl;
```

---

## Level 3: Drawing Histograms from TTree

### Basic 1D histogram:
```cpp
// Draw a branch directly — ROOT creates a histogram automatically
T1->Draw("nsimhtHL");                        // HCAL hit count
T1->Draw("nsimhtEC");                        // ECAL hit count
T1->Draw("momin[0]");                        // Incident momentum
T1->Draw("pidin[0]");                        // Particle ID
```

### With binning control:
```cpp
// Syntax: "branch>>histname(nbins, xmin, xmax)"
T1->Draw("nsimhtHL>>h1(100, 0, 500)");
T1->Draw("momin[0]>>h2(200, 0, 300)");
```

### With selection cuts:
```cpp
// Only events with >10 HCAL hits
T1->Draw("nsimhtHL>>h(100,0,500)", "nsimhtHL>10");

// Only pion events
T1->Draw("momin[0]>>h(100,0,300)", "pidin[0]==211");

// Combine cuts with && or ||
T1->Draw("momin[0]", "nsimhtHL>5 && nsimhtEC>3");
```

### 2D scatter plot:
```cpp
// Syntax: "y:x" — NOTE: y comes first!
T1->Draw("nsimhtHL:nsimhtEC", "", "COLZ");   // HCAL vs ECAL hits
T1->Draw("nsimhtHL:momin[0]", "", "COLZ");   // HCAL hits vs momentum
```

---

## Level 4: Reading Pre-made Histograms

The simulation writes standalone histograms (not in the TTree):

```cpp
TFile *f = TFile::Open("test_run10.root");

// Get a 1D histogram
TH1F *h = (TH1F*)f->Get("ecalenergy");
h->Draw();

// Get a 2D histogram
TH2F *h2d = (TH2F*)f->Get("h2decalenergy");
h2d->Draw("COLZ");

// Available histograms:
// "ecalenergy"            — ECAL total energy (log10 GeV)
// "h2decalenergy"         — ECAL energy 2D (eta vs phi)
// "pPosRT/RE/RH"          — Radial position (Tracker/ECAL/HCAL)
// "pPosET/EE/EH"          — Eta position
// "pPosPT/PE/PH"          — Phi position
// "trkenergy_L0..L12"     — Tracker energy per layer
// "hcalenergy_L0..L16"    — HCAL energy per layer
// "hcal2denergy_L0..L16"  — HCAL 2D energy per layer
// "hcal_total_energy"     — HCAL total energy (GeV)
// "ecal_hcal_total_energy"— ECAL+HCAL total energy (GeV)
```

---

## Level 5: Customizing Plot Appearance

```cpp
TFile *f = TFile::Open("test_run10.root");
TH1F *h = (TH1F*)f->Get("hcal_total_energy");

// Colors
h->SetLineColor(kBlue+1);          // Line color
h->SetLineWidth(2);                 // Line thickness
h->SetFillColor(kBlue-9);           // Fill color
h->SetFillStyle(3004);              // Fill pattern (3004 = hatched)
h->SetMarkerColor(kRed);            // Marker color
h->SetMarkerStyle(20);              // Marker shape (20 = filled circle)

// Titles
h->SetTitle("My Title;X axis label;Y axis label");
h->GetXaxis()->SetTitle("Energy (GeV)");
h->GetYaxis()->SetTitle("Events");
h->GetXaxis()->SetTitleSize(0.05);
h->GetYaxis()->SetTitleSize(0.05);

// Range
h->GetXaxis()->SetRangeUser(0, 200);  // Zoom X axis
h->SetMaximum(50);                     // Set Y max

// Draw options
h->Draw();           // Default (line)
h->Draw("E");        // With error bars
h->Draw("HIST");     // Histogram style (no error bars)
h->Draw("PE");       // Points with errors

// Canvas
TCanvas *c = new TCanvas("c", "My Canvas", 800, 600);
gPad->SetLogy();     // Log Y axis
gPad->SetLogx();     // Log X axis
gPad->SetGrid();     // Grid lines
h->Draw();
c->SaveAs("myplot.pdf");
c->SaveAs("myplot.png");

// Colors reference:
// kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kViolet, kYellow
// Add +1/+2 for darker, -1/-9 for lighter
```

---

## Level 6: Multiple Histograms & Overlays

### Overlay on same pad:
```cpp
TFile *f = TFile::Open("test_run10.root");
TH1F *h1 = (TH1F*)f->Get("hcal_total_energy");
TH1F *h2 = (TH1F*)f->Get("ecal_hcal_total_energy");

h1->SetLineColor(kRed+1);   h1->SetLineWidth(2);
h2->SetLineColor(kBlue+1);  h2->SetLineWidth(2);

h1->Draw();           // Draw first
h2->Draw("SAME");     // Overlay second — "SAME" is the key!

// Add legend
TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);  // x1,y1,x2,y2 (NDC)
leg->AddEntry(h1, "HCAL only", "l");
leg->AddEntry(h2, "ECAL + HCAL", "l");
leg->SetBorderSize(0);
leg->Draw();
```

### Split canvas into pads:
```cpp
TCanvas *c = new TCanvas("c", "Multi-pad", 1200, 600);
c->Divide(2, 1);     // 2 columns, 1 row

c->cd(1);            // Go to pad 1
h1->Draw();

c->cd(2);            // Go to pad 2
h2->Draw();

c->SaveAs("comparison.pdf");
```

### More grid options:
```cpp
c->Divide(3, 2);     // 3 columns, 2 rows = 6 pads
c->Divide(2, 2);     // 2x2 = 4 pads
```

---

## Level 7: 2D Histograms & Color Maps

```cpp
TFile *f = TFile::Open("test_run10.root");
TH2F *h = (TH2F*)f->Get("h2decalenergy");

// Draw options for 2D:
h->Draw("COLZ");     // Color map + Z-axis palette (most common)
h->Draw("LEGO");     // 3D lego plot
h->Draw("LEGO2");    // 3D lego with colors
h->Draw("SURF");     // Surface plot
h->Draw("CONT");     // Contour lines
h->Draw("BOX");      // Box plot (area = value)

// Color palettes:
gStyle->SetPalette(kViridis);      // Scientific default
gStyle->SetPalette(kRainBow);      // Rainbow
gStyle->SetPalette(kBird);         // Blue to yellow
gStyle->SetPalette(kTemperatureMap);// Blue-white-red

// Project 2D to 1D:
TH1D *px = h->ProjectionX("px");   // Project onto X axis
TH1D *py = h->ProjectionY("py");   // Project onto Y axis
px->Draw();

// Profile (mean Y vs X):
TProfile *prof = h->ProfileX("prof");
prof->Draw();
```

---

## Level 8: Fitting Distributions

```cpp
TFile *f = TFile::Open("test_run10.root");
TH1F *h = (TH1F*)f->Get("ecal_hcal_total_energy");

// Built-in fits:
h->Fit("gaus");              // Gaussian
h->Fit("expo");              // Exponential
h->Fit("pol1");              // Linear (1st order polynomial)
h->Fit("pol2");              // Quadratic
h->Fit("landau");            // Landau distribution

// Fit in a range:
h->Fit("gaus", "", "", 50, 250);   // Fit only between 50 and 250

// Quiet fit (no printout):
h->Fit("gaus", "Q");

// Access fit results:
TF1 *fit = h->GetFunction("gaus");
double mean  = fit->GetParameter(1);     // Gaussian mean
double sigma = fit->GetParameter(2);     // Gaussian sigma
double chi2  = fit->GetChisquare();
int    ndf   = fit->GetNDF();
cout << "Mean = " << mean << ", Sigma = " << sigma << endl;
cout << "Resolution = " << sigma/mean * 100 << "%" << endl;

// Custom function fit:
TF1 *myFunc = new TF1("myFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3]", 0, 400);
myFunc->SetParameters(100, 200, 30, 5);  // Initial guesses
myFunc->SetParNames("Amplitude", "Mean", "Sigma", "Background");
h->Fit(myFunc);
```

---

## Level 9: Multi-page PDF Output

```cpp
TCanvas *c = new TCanvas("c", "", 800, 600);

// Open PDF (note the parenthesis in filename!)
c->SaveAs("all_plots.pdf(");    // "(" = open

// ... draw something ...
h1->Draw();
c->SaveAs("all_plots.pdf");     // Middle pages (no parenthesis)

// ... draw something else ...
h2->Draw();
c->SaveAs("all_plots.pdf");     // Another page

// Close PDF
h3->Draw();
c->SaveAs("all_plots.pdf)");    // ")" = close and finalize
```

---

## Level 10: Looping Over Events in C++

Save this as `loop_events.C` and run with `root -l loop_events.C`:

```cpp
void loop_events() {
    TFile *f = TFile::Open("test_run10.root");
    TTree *T1 = (TTree*)f->Get("T1");

    // Declare variables matching branch types
    unsigned int nsimhtHL;
    unsigned long int detidHL[2000];
    unsigned int nsimhtEC;
    float energyEC[2000];
    float momin[50];
    int pidin[50];

    // Connect branches
    T1->SetBranchAddress("nsimhtHL", &nsimhtHL);
    T1->SetBranchAddress("detidHL", detidHL);
    T1->SetBranchAddress("nsimhtEC", &nsimhtEC);
    T1->SetBranchAddress("energyEC", energyEC);
    T1->SetBranchAddress("momin", momin);
    T1->SetBranchAddress("pidin", pidin);

    // Create output histograms
    TH1F *h_total = new TH1F("h_total", "Total HCAL Energy;Energy (GeV);Events",
                              200, 0, 300);

    int nentries = T1->GetEntries();
    for (int i = 0; i < nentries; i++) {
        T1->GetEntry(i);

        // Loop over HCAL hits and sum energy
        double total_energy = 0;
        for (unsigned int j = 0; j < nsimhtHL; j++) {
            // Extract energy from lower 17 bits
            unsigned int energy_MeV = detidHL[j] & 0x1FFFF;
            total_energy += energy_MeV / 1000.0;  // MeV -> GeV
        }

        h_total->Fill(total_energy);

        // Print progress every 100 events
        if (i % 100 == 0) {
            cout << "Event " << i << "/" << nentries
                 << ": HCAL hits=" << nsimhtHL
                 << ", ECAL hits=" << nsimhtEC
                 << ", HCAL energy=" << total_energy << " GeV" << endl;
        }
    }

    h_total->Fit("gaus", "Q");
    h_total->Draw();

    TF1 *fit = h_total->GetFunction("gaus");
    if (fit) {
        cout << "\n=== Calibration Result ===" << endl;
        cout << "Mean  = " << fit->GetParameter(1) << " GeV" << endl;
        cout << "Sigma = " << fit->GetParameter(2) << " GeV" << endl;
        cout << "Resolution = " << fit->GetParameter(2)/fit->GetParameter(1)*100
             << "%" << endl;
    }
}
```

---

## Level 11: Comparing Two Files

```cpp
void compare_particles() {
    // Open two simulation files (e.g., pion vs kaon)
    TFile *f1 = TFile::Open("test_pion_run10.root");
    TFile *f2 = TFile::Open("test_kaon_run10.root");

    TH1F *h1 = (TH1F*)f1->Get("ecal_hcal_total_energy");
    TH1F *h2 = (TH1F*)f2->Get("ecal_hcal_total_energy");

    // Normalize to same area for shape comparison
    if (h1->Integral() > 0) h1->Scale(1.0 / h1->Integral());
    if (h2->Integral() > 0) h2->Scale(1.0 / h2->Integral());

    h1->SetLineColor(kBlue+1); h1->SetLineWidth(2);
    h2->SetLineColor(kRed+1);  h2->SetLineWidth(2);

    TCanvas *c = new TCanvas("c", "Comparison", 800, 600);
    h1->SetTitle("Energy Response: #pi^{+} vs K^{+};Energy (GeV);Normalized");
    h1->Draw("HIST");
    h2->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->AddEntry(h1, "#pi^{+} (200 GeV)", "l");
    leg->AddEntry(h2, "K^{+} (200 GeV)", "l");
    leg->Draw();

    c->SaveAs("pion_vs_kaon.pdf");
}
```

---

## Level 12: Unpacking Packed Detector Words

```cpp
void unpack_detid() {
    TFile *f = TFile::Open("test_run10.root");
    TTree *T1 = (TTree*)f->Get("T1");

    unsigned int nsimhtHL;
    unsigned long int detidHL[2000];
    T1->SetBranchAddress("nsimhtHL", &nsimhtHL);
    T1->SetBranchAddress("detidHL", detidHL);

    TH1F *h_eta   = new TH1F("h_eta",   "ieta;ieta;Hits",   64, 0, 64);
    TH1F *h_phi   = new TH1F("h_phi",   "iphi;iphi;Hits",   64, 0, 64);
    TH1F *h_depth = new TH1F("h_depth", "depth;depth;Hits",  8,  0, 8);
    TH1F *h_energy= new TH1F("h_energy","Energy;E (MeV);Hits", 200, 0, 5000);

    for (int i = 0; i < T1->GetEntries(); i++) {
        T1->GetEntry(i);
        for (unsigned int j = 0; j < nsimhtHL; j++) {
            //  Packed word layout:
            //  |  cellid (15 bits)  | energy (17 bits) |
            //  | ieta[6] iphi[6] depth[3] |  MeV      |

            unsigned int energy = detidHL[j] & 0x1FFFF;         // bits 0-16
            unsigned int cellid = (detidHL[j] >> 17) & 0x7FFF;  // bits 17-31
            unsigned int depth  = cellid & 0x7;                  // bits 0-2
            unsigned int iphi   = (cellid >> 3) & 0x3F;          // bits 3-8
            unsigned int ieta   = (cellid >> 9) & 0x3F;          // bits 9-14

            h_eta->Fill(ieta);
            h_phi->Fill(iphi);
            h_depth->Fill(depth);
            h_energy->Fill(energy);
        }
    }

    TCanvas *c = new TCanvas("c", "Unpacked", 1200, 800);
    c->Divide(2, 2);
    c->cd(1); h_eta->Draw();
    c->cd(2); h_phi->Draw();
    c->cd(3); h_depth->Draw();
    c->cd(4); h_energy->Draw();
    c->SaveAs("unpacked_detid.pdf");
}
```

---

## Level 13: Using PyROOT (Python)

If you prefer Python (with conda `hep` environment):

```python
import ROOT

# Open file
f = ROOT.TFile.Open("test_run10.root")
T1 = f.Get("T1")

# Quick draw
c = ROOT.TCanvas("c", "", 800, 600)
T1.Draw("nsimhtHL")
c.SaveAs("hits.pdf")

# Get pre-made histogram
h = f.Get("hcal_total_energy")
h.SetLineColor(ROOT.kBlue)
h.Fit("gaus")
h.Draw()
c.SaveAs("hcal_energy.pdf")

# Loop over events
for i in range(T1.GetEntries()):
    T1.GetEntry(i)
    print(f"Event {i}: HCAL hits = {T1.nsimhtHL}, ECAL hits = {T1.nsimhtEC}")
```

Run with:
```bash
conda activate hep
python my_analysis.py
```

---

## Available Branches & Histograms

### T1 Tree Branches (Simulation Output)

| Branch | Type | Description |
|--------|------|-------------|
| `irun` | `unsigned int` | Run number |
| `ievt` | `unsigned int` | Event number |
| `ngent` | `unsigned int` | Number of generated particles |
| `pidin[ngent]` | `int[]` | PID of generated particles |
| `momin[ngent]` | `float[]` | Momentum of generated particles |
| `thein[ngent]` | `float[]` | Theta of generated particles |
| `phiin[ngent]` | `float[]` | Phi of generated particles |
| `nsimhtTk` | `unsigned int` | Number of tracker hits |
| `detidTk[nsimhtTk]` | `unsigned int[]` | Tracker hit detector ID |
| `nsimhtEC` | `unsigned int` | Number of ECAL hits |
| `detidEC[nsimhtEC]` | `unsigned int[]` | ECAL packed detector ID |
| `energyEC[nsimhtEC]` | `float[]` | ECAL hit energy (MeV) |
| `thetaEC[nsimhtEC]` | `float[]` | ECAL hit theta |
| `phiEC[nsimhtEC]` | `float[]` | ECAL hit phi |
| `nsimhtHL` | `unsigned int` | Number of HCAL hits |
| `detidHL[nsimhtHL]` | `unsigned long[]` | HCAL packed word (cellid<<17 + energy_MeV) |
| `timeHL[nsimhtHL]` | `unsigned int[]` | HCAL hit time |

### T2 Tree Branches (Digitisation Output)

| Branch | Type | Description |
|--------|------|-------------|
| `totsimenr` | `float` | Total ECAL+HCAL simulated energy (GeV) |
| `totdigienr` | `float` | Total ECAL+HCAL digitised energy (GeV) |
| `hclsimenr` | `float` | HCAL-only simulated energy (GeV) |
| `hcldigienr` | `float` | HCAL-only digitised energy (GeV) |

### Standalone Histograms in Simulation File

| Name | Type | Description |
|------|------|-------------|
| `ecalenergy` | TH1F | ECAL energy (log10 GeV) |
| `h2decalenergy` | TH2F | ECAL energy 2D (η, φ) |
| `hcalenergy_L{i}` | TH1F | HCAL energy per layer (log10 keV) |
| `hcal2denergy_L{i}` | TH2F | HCAL 2D per layer (η, φ) |
| `trkenergy_L{i}` | TH1F | Tracker energy per layer |
| `pPosRH/RE/RT` | TH1F | Radial position (HCAL/ECAL/Trk) |
| `pPosEH/EE/ET` | TH1F | Eta position |
| `hcal_total_energy` | TH1F | **HCAL total energy (GeV)** |
| `ecal_hcal_total_energy` | TH1F | **ECAL+HCAL total energy (GeV)** |

---

# Part C: Track Finder

The track finder reconstructs muon trajectories from silicon tracker data using the **Hough Transformation** method.

---

## C1. Overview

The track finder (`skelaton_track_fitting.C`) reconstructs three key parameters:

| Parameter | Plane | Method |
|-----------|-------|--------|
| **θ** (polar angle) | R-Z | Hough transform: `z = cot(θ) × r + z₀` |
| **φ** (azimuthal angle) | R-φ | Hough transform: `φ = φ₀ + (κ/2) × r` |
| **pT** (transverse momentum) | R-φ | From curvature: `pT = 0.3 × B / κ` |

The CMS detector has **13 silicon tracker layers** at radii from 29 mm to 800 mm, with a **3.8 T** magnetic field along the Z direction.

---

## C2. How to Run

### Easiest: Run entire pipeline from scratch
```bash
cd ~/Desktop/cmsModel
./runit.sh scratch
```
This runs **all 6 steps** automatically: Build → Simulate → Digitise → Plot → Track(no noise) → Track(1% noise).

### Individual commands:
```bash
./runit.sh build          # Just build (cmake + make)
./runit.sh run            # Build + simulate (run.mac)
./runit.sh digi           # Digitise simulation output
./runit.sh plot           # Generate plots (sim + digi)
./runit.sh track          # Track finder, no noise
./runit.sh track 13       # Track finder, 13 noise hits/event
./runit.sh track 50       # Track finder, 50 noise hits (stress test)
./runit.sh all            # Build + Sim + Digi + Plot (no track)
./runit.sh scratch        # EVERYTHING
./runit.sh                # Interactive mode (with visualization)
```

### Manual method (if not using `runit.sh`):
```bash
# 1. Build
cd ~/Desktop/cmsModel/build
cmake -DGeant4_DIR=$HOME/software/geant4/geant4-v11.4.0/lib/cmake/Geant4 \
      -DCMAKE_PREFIX_PATH=$HOME/miniforge3/envs/hep ..
make -j$(sysctl -n hw.ncpu)

# 2. Simulate (ensure PID 13 for muons in run.mac)
./cmsmodel run.mac

# 3. Digitise
echo "test_run10.root 10000" > test_pion_klong.log
./ecal_hcal_digitisation 40 200 25 100

# 4. Plot
root -l -b -q '../plot_results.C("test_run10.root","test_pion_klong_xwoecl_40_200_25_100.root")'

# 5. Track finder
echo "test_run10.root 10000" > skeleton_track.log
./track_finder         # no noise
./track_finder 13      # with noise
```

> **⚠️ Important**: Use `./track_finder` (compiled), NOT `root -l skelaton_track_fitting.C`.

---

## C3. Hough Transform Theory

### R-Z Plane (Non-bending, for θ reconstruction)

In a uniform B field along Z, tracks are straight lines in the R-Z plane:
```
z = cot(θ) × r + z₀
```

**Hough Transform idea**: For every pair of hits `(r₁, z₁)` and `(r₂, z₂)`, compute:
```
slope     = (z₂ - z₁) / (r₂ - r₁)   → cot(θ)
intercept = z₁ - slope × r₁          → z₀
```

Signal hits from the same track all converge to `one point (cotθ, z₀)` in Hough space. Noise hits scatter randomly.

**Iterative 2σ cleanup**: Remove (slope, intercept) pairs outside 2σ of the mean, repeat until stable.

**Final**: `θ = atan2(1, mean_slope)`, refined with straight-line fit.

### R-φ Plane (Bending, for φ and pT reconstruction)

In the magnetic field, the track curves:
```
φ_hit ≈ φ₀ + (κ/2) × r
```

From the **slope** (= κ/2):
- `φ₀ = intercept` → reconstructed azimuthal angle
- `pT = (0.3 × 3.8) / (2 × slope)` → transverse momentum in GeV

---

## C4. Configurable Parameters

All parameters are at the top of `skelaton_track_fitting.C` (lines 72–78):

| Constant | Default | Description |
|----------|---------|-------------|
| `nNoiseHit` | `0` | **Runtime CLI arg**: `./track_finder 13` (no recompile!) |
| `nsilayer` | `13` | Number of silicon layers |
| `nmnhit` | `3` | Minimum hits required for a track |
| `nEventDisplay` | `4` | Number of events with per-event plots |
| `hough_nsigma` | `2.0` | σ-cut for iterative Hough cleanup |
| `hough_maxiter` | `10` | Maximum cleanup iterations |

### Position smearing (in `fill_digipos_array`):
| Effect | Code | Description |
|--------|------|-------------|
| Z smearing | `gRandom3->Gaus(0.0, 0.2)` | σ = 200 μm in Z |
| φ smearing | `gRandom3->Gaus(0.0, 0.05)` | σ = 50 mrad in φ |
| R-φ position | `gRandom3->Gaus(0.0, 1.e-4)` | σ = 100 μm in r×φ |

### Hit efficiency:
- Outer layers (irad > 6): 90% efficiency
- Inner layers (irad ≤ 6): 90% × 90% = 81% efficiency

### Delta ray removal:
Enabled via `#define DELTARAY` — keeps only muon hits (PDG=±13) with time < 10 ns.

---

## C5. Noise Comparison

Noise is controlled via **command-line argument** (no recompile needed):

```bash
# Step 1: Run with no noise (default)
./track_finder
cp track_resolution.pdf track_resolution_no_noise.pdf
cp track_event_display.pdf track_event_display_no_noise.pdf

# Step 2: Run with 13 noise hits (~1%)
./track_finder 13
cp track_resolution.pdf track_resolution_1pct_noise.pdf
cp track_event_display.pdf track_event_display_1pct_noise.pdf

# Step 3: Run with more noise (stress test)
./track_finder 50
```

> **No recompile needed!** The noise count is passed directly as a CLI argument.

**Expected effect**: Noise adds random red dots in event displays and widens resolution tails, especially ΔpT.

---

## C6. Output Files

| File | Contents |
|------|----------|
| `output_skeleton.root` | ROOT file with T2 tree + resolution histograms |
| `track_event_display.pdf` | Per-event scatter plots (first 4 events) |
| `track_resolution.pdf` | Resolution (ΔpT, Δφ, Δθ) + summary plots |

### `track_event_display.pdf` pages (per event):
| Pad | X axis | Y axis | Colors |
|-----|--------|--------|--------|
| 1 | ρ (m) | Z (m) | Black=signal, Red=noise |
| 2 | X (m) | Y (m) | Black=signal, Red=noise |
| 3 | Slope (cotθ) | Intercept (m) | Green=signal, Black=mixed, Red=noise |
| 4 | Slope (zoom) | Intercept (zoom) | Same, zoomed on signal cluster |

### `track_resolution.pdf` pages:
| Page | Plots |
|------|-------|
| 1 | ΔpT (GeV), Δφ (mrad), Δθ (mrad) — all with **log Y** |
| 2 | Hits/event, Reconstructed pT, Tracks found |

---

## C7. Understanding Resolution Plots

The resolution histograms show:
```
Δθ = (θ_reco - θ_gen) × 1000     [mrad]
Δφ = (φ_reco - φ_gen) × 1000     [mrad]
ΔpT = pT_reco - pT_gen             [GeV]
```

**Interpreting results**:
- **Mean ≈ 0**: Unbiased reconstruction
- **StdDev**: Resolution width (smaller = better)
- Non-bending plane (θ) has much better resolution than bending plane (φ, pT)
- **Log Y scale** reveals non-Gaussian tails from noise contamination

### Reference values (from assignment slides, ~10k events):
| Quantity | Mean | StdDev |
|----------|------|--------|
| ΔpT | 0.014 GeV | 5.88 GeV |
| Δφ | 0.0009 mrad | 0.80 mrad |
| Δθ | -0.0003 mrad | 0.38 mrad |

---

## C8. Multi-Track Finding

The algorithm finds multiple tracks iteratively:

1. **Run Hough R-Z** on all hits → find first track
2. **Refine with Hough R-φ** → get φ and pT
3. **Remove** found track's hits from the hit list
4. **Repeat** steps 1–3 on remaining hits
5. **Stop** when fewer than `nmnhit` (3) hits remain

The number of found tracks is stored in the `ntrack_found` branch.

---

## C9. T2 Tree Branches (Track Finder)

| Branch | Type | Description |
|--------|------|-------------|
| `irun` | `unsigned int` | Run number |
| `ievt` | `unsigned int` | Event number |
| `gen_theta` | `float` | Generated θ (rad) |
| `gen_phi` | `float` | Generated φ (rad) |
| `gen_pt` | `float` | Generated pT (GeV) |
| `gen_mom` | `float` | Generated momentum (GeV) |
| `reco_theta` | `float` | Reconstructed θ (rad) |
| `reco_phi` | `float` | Reconstructed φ (rad) |
| `reco_pt` | `float` | Reconstructed pT (GeV) |
| `delta_theta_mrad` | `float` | Δθ = θ_reco - θ_gen (mrad) |
| `delta_phi_mrad` | `float` | Δφ = φ_reco - φ_gen (mrad) |
| `delta_pt_GeV` | `float` | ΔpT = pT_reco - pT_gen (GeV) |
| `nhits_used` | `int` | Hits used for first track |
| `nhits_total` | `int` | Total hits in event |
| `ntrack_found` | `int` | Number of tracks found |

### Reading the track tree:
```cpp
root -l output_skeleton.root
T2->Draw("delta_theta_mrad");     // θ resolution
T2->Draw("delta_phi_mrad");       // φ resolution
T2->Draw("delta_pt_GeV");         // pT resolution
T2->Draw("reco_pt:gen_pt", "", "COLZ");  // 2D correlation
T2->Draw("nhits_used");           // Track quality
T2->Draw("ntrack_found");         // Multi-track stat
```

---

## C10. Troubleshooting

| Problem | Solution |
|---------|----------|
| `root -l skelaton_track_fitting.C` fails | Use compiled `./track_finder` instead |
| `Cannot open skeleton_track.log` | Create: `echo "test_run10.root 10000" > skeleton_track.log` |
| All ΔpT values very large | Check `run.mac` uses PID 13 (muons), not 211 (pions) |
| No events processed | Ensure `test_run10.root` exists with T1 tree |
| Very wide Δθ | Check smearing parameters, try smaller σ |
| Build error | Run `cmake ... && make track_finder` from `build/` |
| "T1 tree not found" | Input file may be a digitisation file; use simulation output |
