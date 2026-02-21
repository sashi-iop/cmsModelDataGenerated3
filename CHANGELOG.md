# CMS Model — Complete Changelog & Developer Notes

> **Date**: 21 February 2026  
> **Project**: `~/Desktop/cmsModel` — Geant4 + ROOT CMS Detector Simulation  
> **Geant4 Version**: 11.4.0 (`~/software/geant4/geant4-v11.4.0/`)  
> **ROOT**: via conda `hep` environment (`~/miniforge3/envs/hep/`)

---

## Table of Contents

1. [Build System Fixes](#1-build-system-fixes)
2. [Compilation Error Fixes](#2-compilation-error-fixes)
3. [Runtime Crash Fix](#3-runtime-crash-fix)
4. [Hands-On Session: Primary Particle Config](#4-hands-on-session-primary-particle-configuration)
5. [Hands-On Session: Energy Resolution](#5-hands-on-session-energy-resolution-in-hcal)
6. [Hands-On Session: 32-bit Word Packing](#6-hands-on-session-32-bit-word-packing)
7. [Hands-On Session: Calibration Histograms](#7-hands-on-session-calibration-histograms)
8. [Digitisation Macro](#8-digitisation-macro)
9. [Build Script](#9-build-script)
10. [How to Build & Run](#10-how-to-build--run)
11. [Track Finder: CLHEP → TVector3](#11-track-finder-clhep--tvector3)
12. [Track Finder: Hough Transform R-Z](#12-track-finder-hough-transform-r-z)
13. [Track Finder: Hough Transform R-φ](#13-track-finder-hough-transform-r-φ)
14. [Track Finder: Resolution Plots](#14-track-finder-resolution-plots)
15. [Track Finder: Multi-Track Finding](#15-track-finder-multi-track-finding)
16. [Track Finder: Noise Comparison](#16-track-finder-noise-comparison)
17. [Track Finder: Per-Event Display](#17-track-finder-per-event-display)
18. [Track Finder: Build & Run](#18-track-finder-build--run)
19. [Track Finder: Noise as CLI Argument](#19-track-finder-noise-as-command-line-argument)
20. [Build Safety: Preserve Data on Clean](#20-build-safety-preserve-data-on-clean-rebuild)
21. [Full Pipeline Script (`runit.sh` rewrite)](#21-full-pipeline-script-runitsh-rewrite)
22. [Plot Macro Default Filename Fix](#22-plot-macro-default-filename-fix)

---

## 1. Build System Fixes

### Problem
`cmake` could not find Geant4 or ROOT packages.

### Solution
Pass explicit paths to cmake:
```bash
cmake -DGeant4_DIR=$HOME/software/geant4/geant4-v11.4.0/lib/cmake/Geant4 \
      -DCMAKE_PREFIX_PATH=$HOME/miniforge3/envs/hep \
      ..
```
- `-DGeant4_DIR` points to the Geant4 cmake config
- `-DCMAKE_PREFIX_PATH` points to the conda env so ROOT and its dependency `Vdt` are found

### File Changed
**None** — this is a cmake command-line fix, not a code change.

---

## 2. Compilation Error Fixes

### Problem
Three "variable-sized object may not be initialized" errors.  
In C++, you cannot initialize a variable-length array (VLA) with `= {values}`.  
Accessing `pAnalysis->nsilayer` through a pointer makes the compiler treat it as a runtime value, even though `nsilayer` is `static const`.

### Fix: Use the class constant directly

#### File: `src/Serc19DetectorConstruction.cc` (lines 250, 252)
```diff
-  const int nsilayer=pAnalysis->nsilayer;
-  double sirad[pAnalysis->nsilayer]={2.9, 6.8, ...};
+  const int nsilayer=Serc19SimAnalysis::nsilayer;
+  double sirad[Serc19SimAnalysis::nsilayer]={2.9, 6.8, ...};
```

#### File: `src/Serc19EventAction.cc` (lines 153, 176)
```diff
-  G4double totET[pAnalysis->nsilayer] = {0};
+  G4double totET[Serc19SimAnalysis::nsilayer] = {0};

-  G4double totEH[pAnalysis->nhcalLayer] = {0};
+  G4double totEH[Serc19SimAnalysis::nhcalLayer] = {0};
```

### Why it works
`Serc19SimAnalysis::nsilayer` is a compile-time constant (static const), so the compiler knows the array size at compile time.

---

## 3. Runtime Crash Fix

### Problem
Running `./cmsmodel` without arguments caused a **segmentation fault**.

### Root Cause
Line 41 of `serc19_cmsmodel.cc` unconditionally accessed `argv[1]` and `argv[2]`, which are NULL when `argc == 1`.

### Fix
#### File: `serc19_cmsmodel.cc` (line 41)
```diff
-  G4cout <<"argc "<<argc<<" "<<argv[0]<<" "<<argv[1]<<" "<<argv[2]<<G4endl;
+  G4cout <<"argc "<<argc<<" "<<argv[0];
+  if (argc > 1) G4cout <<" "<<argv[1];
+  if (argc > 2) G4cout <<" "<<argv[2];
+  G4cout <<G4endl;
```

---

## 4. Hands-On Session: Primary Particle Configuration

### Change
Updated beam energy from 15 GeV to 200 GeV for π⁺ study.

#### File: `src/Serc19PrimaryGeneratorAction.cc` (line 51)
```diff
-  SetIncEnergy(15.0*GeV);
+  SetIncEnergy(200.0*GeV);
```

`SetPartId(211)` (pion) was already set — no change needed.

### For kaon comparison (future)
```cpp
SetPartId(321);          // K+ instead of π+
SetIncEnergy(40.0*GeV);  // or 100.0*GeV
```

---

## 5. Hands-On Session: Energy Resolution in HCAL

### What was implemented
Three effects applied in `EndOfEvent()` after cell-level energy accumulation:

#### File: `src/Serc19HclSD.cc` — `EndOfEvent()` method

**a) Photon Statistics (Stochastic Term)**
```cpp
#include "G4Poisson.hh"
// ...
G4double npe_mean = 40.0 * edep_MeV;            // 40 p.e. per MeV
G4long   npe_sampled = G4Poisson(npe_mean);      // Poisson sampling
G4double edep_smeared_MeV = edep_MeV * (G4double(npe_sampled) / npe_mean);
```

**b) Collection Efficiency (depth-dependent)**
```cpp
G4double efficiency = 1.0 + 0.1 * G4double(approx_depth) / G4double(nhcalLayer);
edep_smeared_MeV *= efficiency;
```

**c) Electronic Noise**
```cpp
#include "CLHEP/Random/RandGauss.h"
// ...
edep_smeared_MeV += G4RandGauss::shoot(0.0, 200.0);  // σ = 200 MeV
```

### New includes added to `src/Serc19HclSD.cc`
```cpp
#include "Randomize.hh"
#include "CLHEP/Random/RandGauss.h"
#include "G4Poisson.hh"
```

---

## 6. Hands-On Session: 32-bit Word Packing

### Old layout (17 bits for detid)
```
ieta(6) << 11 | iphi(6) << 5 | idepth(5)
```

### New layout (32 bits total)
```
| ieta (6 bits) | iphi (6 bits) | depth (3 bits) | energy (17 bits) |
|   bits 31-26  |  bits 25-20   |   bits 19-17   |    bits 16-0     |
```

#### File: `src/Serc19HclSD.cc` — `ProcessHits()` method
```cpp
// 17 layers mapped to 3 bits
unsigned int depth3 = (idepth / 3);
if (depth3 > 7) depth3 = 7;  // cap at 3-bit max

unsigned int detid = ieta;
detid <<= 6;
detid += iphi;
detid <<= 3;
detid += depth3;
```

#### File: `src/Serc19HclSD.cc` — `EndOfEvent()` method
```cpp
// Pack into 32-bit word: cellid(15 bits) << 17 | energy(17 bits)
G4int energy_MeV = (edep_smeared_MeV > 0) ? G4int(edep_smeared_MeV) : 0;
if (energy_MeV > 131071) energy_MeV = 131071;  // 2^17 - 1

unsigned long int packed_word = (static_cast<unsigned long int>(cellid) << 17) + energy_MeV;
pAnalysis->detidHL[pAnalysis->nsimhtHL] = packed_word;
```

#### File: `src/Serc19EventAction.cc` — `EndOfEventAction()` (depth extraction updated)
```diff
-  unsigned il = ((*EHC2)[ij]->GetHitId()) & 0x1F;  // old: 5-bit depth
+  unsigned il = detid & 0x7;                         // new: 3-bit depth
+  unsigned layer = il * 3;                            // map back to ~original layer
```

---

## 7. Hands-On Session: Calibration Histograms

### New histogram declarations
#### File: `include/Serc19SimAnalysis.hh`
```cpp
// Calibration histograms (with and without ECAL)
TH1F* h_hcal_total_energy;        // Total HCAL energy only
TH1F* h_ecal_hcal_total_energy;   // Combined ECAL + HCAL energy
```

### Histogram creation
#### File: `src/Serc19SimAnalysis.cc` — in `OpenRootfiles()`
```cpp
h_hcal_total_energy = new TH1F("hcal_total_energy",
    "HCAL Total Energy (GeV)", 200, 0, 400);
h_ecal_hcal_total_energy = new TH1F("ecal_hcal_total_energy",
    "ECAL+HCAL Total Energy (GeV)", 200, 0, 400);
```

### Histogram filling
#### File: `src/Serc19EventAction.cc` — in `EndOfEventAction()`
```cpp
// ECAL energy moved to wider scope
G4double totEE = 0;
// ... ECAL loop fills totEE ...

// HCAL total accumulated
G4double totEH_all = 0;
// ... HCAL loop fills totEH_all ...

// Fill calibration histograms
G4double hcal_GeV = totEH_all / 1.0e6;  // keV -> GeV
G4double ecal_GeV = totEE / GeV;
pAnalysis->h_hcal_total_energy->Fill(hcal_GeV);
pAnalysis->h_ecal_hcal_total_energy->Fill(ecal_GeV + hcal_GeV);
```

---

## 8. Digitisation Macro

### File: `skeleton_ecal_hcal_digitisation.C`

The skeleton was provided with empty digitisation logic. The following was implemented:

#### a) Bug fix: `detidHL` type mismatch
```diff
-  unsigned int detidHL[nsimhtmxHL];     // WRONG: doesn't match /l branch
+  unsigned long int detidHL[nsimhtmxHL]; // CORRECT: matches 64-bit branch
```

#### b) ECAL Digitisation (crystal-level)
```cpp
for (unsigned int ih = 0; ih < nsimhtEC; ih++) {
    float sim_energy_GeV = energyEC[ih] / 1000.0;  // MeV -> GeV
    float smeared_energy = sim_energy_GeV + gRandom3->Gaus(0, pedwidth);
    if (smeared_energy > otherthrs) {
        ecal_digi_total += smeared_energy;
    }
}
```

#### c) HCAL Digitisation (32-bit unpacking + noise)
```cpp
for (unsigned int ih = 0; ih < nsimhtHL; ih++) {
    unsigned int energy_MeV = detidHL[ih] & 0x1FFFF;  // lower 17 bits
    float sim_energy_GeV = energy_MeV / 1000.0;
    float smeared_energy = sim_energy_GeV + gRandom3->Gaus(0, hclwidth);
    if (smeared_energy > hclthrs) {
        hcal_digi_total += smeared_energy;
    }
}
```

#### d) Removed CLHEP dependency
Commented out `CLHEP/Matrix/Matrix.h`, `CLHEP/Vector/LorentzVector.h`, `CLHEP/Vector/ThreeVector.h` and the `dgap()` function — not needed for digitisation.

#### e) Added to `CMakeLists.txt`
```cmake
add_executable(ecal_hcal_digitisation skeleton_ecal_hcal_digitisation.C)
target_include_directories(ecal_hcal_digitisation PRIVATE include)
target_link_libraries(ecal_hcal_digitisation ${ROOT_LIBRARIES})
target_link_libraries(ecal_hcal_digitisation ${Geant4_LIBRARIES})
```

### Output
Produces a ROOT file with tree `T2` containing branches:
- `totsimenr` — ECAL + HCAL simulated energy (GeV)
- `totdigienr` — ECAL + HCAL digitised energy (GeV)
- `hclsimenr` — HCAL-only simulated energy (GeV)
- `hcldigienr` — HCAL-only digitised energy (GeV)

---

## 9. Build Script

### File: `runit.sh`

Complete build/run/digitise script with commands:
```
./runit.sh build   — Build only (incremental)
./runit.sh clean   — Clean + full rebuild
./runit.sh run     — Build + batch simulation
./runit.sh digi    — Run digitisation
./runit.sh all     — Full pipeline: build → simulate → digitise
./runit.sh         — Build + interactive mode
```

---

## 10. How to Build & Run

### Prerequisites
```bash
conda activate hep
```

### Full pipeline
```bash
cd ~/Desktop/cmsModel
./runit.sh all
```

### Manual step-by-step
```bash
# Build
cd ~/Desktop/cmsModel
source ~/software/geant4/geant4-v11.4.0/bin/geant4.sh
cd build
cmake -DGeant4_DIR=$HOME/software/geant4/geant4-v11.4.0/lib/cmake/Geant4 \
      -DCMAKE_PREFIX_PATH=$HOME/miniforge3/envs/hep ..
make -j$(sysctl -n hw.ncpu)

# Simulate (batch mode)
./cmsmodel run.mac

# Digitise
echo "test_run10.root 10000" > test_pion_klong.log
./ecal_hcal_digitisation 40 200 25 100
```

### Digitisation arguments (all in MeV)
| Arg | Meaning | Typical Value |
|-----|---------|---------------|
| 1 | ECAL noise per crystal | 40 |
| 2 | ECAL energy threshold | 200 |
| 3 | HCAL noise per tower | 25 |
| 4 | HCAL energy threshold | 100 |

---

---

## 11. Track Finder: CLHEP → TVector3

### Problem
The skeleton code `skelaton_track_fitting.C` used CLHEP headers (`Hep3Vector`, `CLHEP/Matrix/Matrix.h`, etc.) which are unavailable in ROOT's C++ interpreter and not linked for standalone compilation.

### Fix
Replaced all `Hep3Vector` with ROOT's `TVector3`. Key API differences:
```diff
-#include "CLHEP/Matrix/Matrix.h"
-#include "CLHEP/Vector/LorentzVector.h"
-#include "CLHEP/Vector/ThreeVector.h"
-using namespace CLHEP;
+#include "TVector3.h"

// Hep3Vector methods → TVector3 equivalents:
-  Hep3Vector v;        →  TVector3 v;
-  v.setRhoPhiZ(r,p,z); →  v = MakeRhoPhiZ(r, p, z);  // custom helper
-  v.rho()              →  v.Perp()
-  v.z()                →  v.Z()
-  v.phi()              →  v.Phi()
```

#### Helper functions added (since TVector3 lacks cylindrical setters):
```cpp
TVector3 MakeRhoPhiZ(double rho, double phi, double z) {
  return TVector3(rho * cos(phi), rho * sin(phi), z);
}

TVector3 MakeRhoPhiTheta(double rho, double phi, double theta) {
  double z = rho / tan(theta);
  return TVector3(rho * cos(phi), rho * sin(phi), z);
}
```

### Bug fix: `simenrTk == 0`
```diff
-  if (simenrTk == 0)       // BUG: comparing array pointer to NULL
+  if (nsimhtTk == 0)       // FIX: check hit count instead
```

---

## 12. Track Finder: Hough Transform R-Z

### Algorithm (reconstructs polar angle θ)

In the R-Z plane, a straight track follows: `z = slope × r + intercept`  
where `slope = cot(θ)` and `intercept = z₀` (vertex).

```cpp
// For every pair of hits (i, j):
double slope = (z_j - z_i) / (r_j - r_i);   // = cot(θ)
double intercept = z_i - slope * r_i;         // = z₀
```

#### Iterative 2σ cleanup:
```cpp
for (int iter = 0; iter < 10; iter++) {
  // 1. Calculate mean and RMS of (slope, intercept)
  // 2. Remove points outside 2σ of mean
  // 3. Repeat until stable
}
```

#### Final θ reconstruction:
```cpp
reco_theta = atan2(1.0, mean_slope);  // slope = cot(θ) → θ = atan(1/slope)
if (reco_theta < 0) reco_theta += π;
```

Refined with `straight_line_fit()` using selected hits.

---

## 13. Track Finder: Hough Transform R-φ

### Algorithm (reconstructs φ and pT)

In a 3.8T magnetic field, charged particles curve in the R-φ plane:
```
φ_hit ≈ φ₀ + (κ/2) × r
```
where `κ = charge / (pT × R_factor)` is the curvature.

```cpp
// For every pair of hits:
double dphi = sfdph(phi_i, phi_j);        // Δφ with wrapping
double slope = dphi / (r_j - r_i);         // ≈ κ/2
double intercept = phi_i - slope * r_i;    // ≈ φ₀

// After iterative cleanup:
reco_phi = mean_intercept;                 // = φ₀
reco_pt = momconv / (2.0 * mean_slope);    // momconv = 0.299792 × 3.8T
```

---

## 14. Track Finder: Resolution Plots

Resolution computed in the units matching the assignment slides:

```cpp
delta_theta_mrad = (reco_theta - gen_theta) * 1000.0;  // rad → mrad
delta_phi_mrad = sfdph(gen_phi, reco_phi) * 1000.0;    // rad → mrad
delta_pt_GeV = reco_pt - gen_pt;                        // absolute GeV
```

### Typical results (100 muon events, 45 GeV):
| Quantity | Mean | StdDev |
|----------|------|--------|
| Δθ | -0.019 mrad | 0.637 mrad |
| Δφ | 0.664 mrad | 2.623 mrad |
| ΔpT | 11.74 GeV | 9.478 GeV |

All plotted with **log Y scale** to show tails, matching slide 3.

---

## 15. Track Finder: Multi-Track Finding

After finding the first track:
1. Remove its hits from the list
2. Repeat Hough transform on remaining hits
3. Continue until < 3 hits remain

```cpp
while (leftover_hits.size() >= nmnhit) {
  double theta = hough_rz(leftover_hits, ...);
  if (selected_hits.size() < nmnhit) break;
  ntrack_found++;
  // Remove found hits from leftover list
}
```

---

## 16. Track Finder: Noise Comparison

To test robustness with 1% noise hits, edit line 73:
```cpp
const int nNoiseHit = 0;   // Default: no noise
const int nNoiseHit = 13;  // ~1% noise (for comparison)
```

Noise hits are randomly generated in (ρ, φ, θ) space across all 13 silicon layers.
Signal hits are tagged `isSignal = true`, noise as `isSignal = false`.

---

## 17. Track Finder: Per-Event Display

The first 4 events produce detailed scatter plots:
- **Pad 1**: ρ vs Z — signal (black ●) + noise (red ·)
- **Pad 2**: X vs Y — same color coding
- **Pad 3**: Hough space (slope vs intercept) with signal (green), mixed (black), noise (red)
- **Pad 4**: Zoomed Hough inset around signal cluster

Output: `track_event_display.pdf`

---

## 18. Track Finder: Build & Run

### Added to `CMakeLists.txt`:
```cmake
add_executable(track_finder skelaton_track_fitting.C)
target_include_directories(track_finder PRIVATE include)
target_link_libraries(track_finder ${ROOT_LIBRARIES})
```

### How to run:
```bash
cd ~/Desktop/cmsModel/build
echo "test_run10.root 10000" > skeleton_track.log
./track_finder
```

> **Important**: Use `./track_finder` (compiled executable), NOT `root -l skelaton_track_fitting.C` (the file uses `int main()`, not a ROOT macro function).

### Output files:
- `output_skeleton.root` — T2 tree with gen/reco parameters
- `track_event_display.pdf` — per-event hit + Hough space plots
- `track_resolution.pdf` — ΔpT, Δφ, Δθ resolution + summary

---

## 19. Track Finder: Noise as Command-Line Argument

### Problem
`nNoiseHit` was a compile-time `const`, requiring recompilation to change.
Users editing the value and running without `make track_finder` saw no noise.

### Fix
#### File: `skelaton_track_fitting.C`
```diff
-const int nNoiseHit = 13;  // compile-time constant
+int nNoiseHit = 0;          // runtime, set via CLI
```

```cpp
int main(int argc, char *argv[]) {
  if (argc > 1) {
    nNoiseHit = atoi(argv[1]);
    if (nNoiseHit < 0) nNoiseHit = 0;
  }
```

### New usage:
```bash
./track_finder          # No noise (default)
./track_finder 13       # 13 noise hits per event (~1%)
./track_finder 50       # 50 noise hits (stress test)
```

---

## 20. Build Safety: Preserve Data on Clean Rebuild

### Problem
`./runit.sh clean` ran `rm -rf "$BUILD_DIR"` which **destroyed** all ROOT output (simulation data, digitisation, track finder output, PDFs) in the build directory.

### Fix
#### File: `runit.sh` — `do_clean()` function

Now backs up `.root`, `.pdf`, `.log` files to `/tmp` before deleting, rebuilds, then restores:
```bash
# Before rm -rf:
BACKUP_DIR="/tmp/cmsmodel_backup_$$"
cp "$BUILD_DIR"/*.root "$BACKUP_DIR/"
cp "$BUILD_DIR"/*.pdf  "$BACKUP_DIR/"
cp "$BUILD_DIR"/*.log  "$BACKUP_DIR/"

rm -rf "$BUILD_DIR"
do_build

# After rebuild:
cp "$BACKUP_DIR"/* "$BUILD_DIR/"
rm -rf "$BACKUP_DIR"
```

### What is preserved:
| Extension | Examples |
|-----------|----------|
| `.root` | `test_run10.root`, `output_skeleton.root`, `ecal_hcal_digi.root` |
| `.pdf` | `track_resolution.pdf`, `track_event_display.pdf`, `plots_simulation.pdf` |
| `.log` | `skeleton_track.log`, `test_pion_klong.log` |

### What is NOT preserved (intentionally):
- Compiled binaries (`.o`, executables) — these get rebuilt
- CMake cache files — regenerated by cmake
- Makefiles — regenerated by cmake

> **Safe commands**: `./runit.sh build` and `make` **never** delete any files. Only `./runit.sh clean` triggered the deletion, and it's now safe.

---

## 21. Full Pipeline Script (`runit.sh` rewrite)

### What changed
Complete rewrite of `runit.sh` adding 3 new commands and a full `scratch` pipeline.

### New commands:
```bash
./runit.sh track [N]   # Track finder (N = noise hits, 0 if omitted)
./runit.sh plot        # Generate simulation + digitisation plots
./runit.sh scratch     # FULL pipeline: build → sim → digi → plot → track(0) → track(13)
```

### `scratch` pipeline (runs all 6 steps):
| Step | Command | Output |
|------|---------|--------|
| 1 | Build | `cmsmodel`, `ecal_hcal_digitisation`, `track_finder` |
| 2 | Simulate | `test_run10.root` |
| 3 | Digitise | `test_pion_klong_xwoecl_40_200_25_100.root` |
| 4 | Plot | `plots_simulation.pdf`, `plots_digitisation.pdf` |
| 5 | Track (no noise) | `track_resolution_no_noise.pdf`, `track_event_display_no_noise.pdf` |
| 6 | Track (1% noise) | `track_resolution_noise13.pdf`, `track_event_display_noise13.pdf` |

### Other improvements:
- Visual step banners showing progress
- `do_plot()` auto-detects simulation and digitisation files
- `do_track()` auto-creates `skeleton_track.log` and saves copies with noise label
- Better help message with `./runit.sh help`

---

## 22. Plot Macro Default Filename Fix

### Problem
`plot_results.C` defaulted to `test_run10_smr001.root` which doesn't exist in `build/`.

### Fix
#### File: `plot_results.C` (line 14)
```diff
-void plot_results(const char *simFile = "test_run10_smr001.root",
+void plot_results(const char *simFile = "test_run10.root",
```

Now works from `build/` without arguments:
```bash
cd build
root -l ../plot_results.C
```

---

## Files Modified Summary

| File | Changes |
|------|---------|
| `serc19_cmsmodel.cc` | Guard `argv[1]`/`argv[2]` access |
| `src/Serc19DetectorConstruction.cc` | VLA fix: `Serc19SimAnalysis::nsilayer` |
| `src/Serc19EventAction.cc` | VLA fix + calibration histograms + 3-bit depth mask |
| `src/Serc19PrimaryGeneratorAction.cc` | Beam energy 15→200 GeV |
| `src/Serc19HclSD.cc` | Energy resolution + 32-bit packing (major rewrite) |
| `include/Serc19SimAnalysis.hh` | Added calibration histogram pointers |
| `src/Serc19SimAnalysis.cc` | Created calibration histograms |
| `skeleton_ecal_hcal_digitisation.C` | Completed digitisation logic |
| `skelaton_track_fitting.C` | CLHEP→TVector3 + Hough + noise CLI arg (major rewrite) |
| `CMakeLists.txt` | Added `ecal_hcal_digitisation` + `track_finder` targets |
| `runit.sh` | Full pipeline script: build/run/digi/plot/track/scratch |
| `plot_results.C` | Plotting macro + default filename fix |
| `GUIDE.md` | Complete guide: params + ROOT analysis + track finder |
| `CHANGELOG.md` | This file |

