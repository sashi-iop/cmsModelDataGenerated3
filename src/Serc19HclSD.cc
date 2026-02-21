
#include "Serc19HclSD.hh"
#include "Serc19HclHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "Randomize.hh"
#include "CLHEP/Random/RandGauss.h"
#include "G4Poisson.hh"


Serc19HclSD::Serc19HclSD(G4String name)
  :G4VSensitiveDetector(name),numberInCell(20000), InCell(0) 
{
  G4String HCname;
  collectionName.insert(HCname="HclCollect");
  pAnalysis = Serc19SimAnalysis::AnPointer;

}

Serc19HclSD::~Serc19HclSD(){;}

void Serc19HclSD::Initialize(G4HCofThisEvent* HCE)
{
  static int HCID = -1;
  HclCollection = new Serc19HclHitsCollection
                      (SensitiveDetectorName,collectionName[0]); 
  if(HCID<0)
  { HCID = GetCollectionID(0); }
  HCE->AddHitsCollection(HCID,HclCollection);
}

G4bool Serc19HclSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
  
  G4double edep = aStep->GetTotalEnergyDeposit();
  
  //  if(edep==0.) return false;
  
  G4TouchableHistory* theTouchable = 
    (G4TouchableHistory*)( aStep->GetPreStepPoint()->GetTouchable() );
  
  //  int level = theTouchable->GetHistoryDepth();
  G4ThreeVector glbpos = 0.5*(aStep->GetPreStepPoint()->GetPosition() + aStep->GetPostStepPoint()->GetPosition());

  pAnalysis->pPosRH->Fill(glbpos.rho(), edep);
  pAnalysis->pPosEH->Fill(glbpos.eta(), edep);
  pAnalysis->pPosPH->Fill(glbpos.phi(), edep);  
  
  int iphi = theTouchable->GetCopyNumber( 3 );
  int idepth = theTouchable->GetCopyNumber( 2 );
  int ieta = theTouchable->GetCopyNumber( 0 );

  pAnalysis->h2d_hcalenergy[idepth]->Fill(glbpos.eta(), glbpos.phi(), edep/MeV);

  if (edep <10*eV) return false;
  
  edep /=keV;

  G4int nInT = G4int(10*aStep->GetPreStepPoint()->GetGlobalTime()/ns); //in 100 ps

  //------------------------------------------------------------
  // New 32-bit word packing (assignment spec):
  //   ieta   : 6 bits (bits 9-14 of cellid)
  //   iphi   : 6 bits (bits 3-8 of cellid)
  //   depth  : 3 bits (bits 0-2 of cellid)
  //   17 layers mapped to 3 bits: depth3 = min(idepth/3, 7)
  //------------------------------------------------------------
  unsigned int depth3 = (idepth / 3);
  if (depth3 > 7) depth3 = 7;  // cap at 3-bit max

  unsigned int detid = ieta;
  detid <<= 6;
  detid += iphi;
  detid <<= 3;
  detid += depth3;

  int oldCellId = -1;
  for (int ij=0; ij<InCell; ij++) {
    if (detid ==CellDetID[ij]) {oldCellId = ij; break;}
  }
  
  if (oldCellId ==-1 && InCell <numberInCell -1 ) {
    Serc19HclHit* newHit = new Serc19HclHit();
    newHit->SetHitId(detid);
    newHit->SetEdep( edep );
    newHit->SetPos(glbpos);
    newHit->SetTime(nInT);

    InCell = HclCollection->insert( newHit );
    CellDetID[InCell-1] = detid;

  } else {
    (*HclCollection)[oldCellId]->AddEdep(edep);
    if (nInT <(*HclCollection)[oldCellId]->GetTime()) {
      (*HclCollection)[oldCellId]->SetTime(nInT);
    }
  }

  return true;
}

void Serc19HclSD::EndOfEvent(G4HCofThisEvent*) {

  InCell = 0;
  pAnalysis->nsimhtHL = 0; // HclCollection->entries();
  
  for (int ij=0; ij< int(HclCollection->entries()); ij++) {

    G4double raw_edep_keV = (*HclCollection)[ij]->GetEdep();  // in keV
    G4double edep_MeV = raw_edep_keV / 1000.0;                // convert to MeV

    //------------------------------------------------------------
    // Energy Resolution Implementation (assignment spec):
    //
    // 1) Photon statistics (stochastic term):
    //    Average photo-electrons = 40 per MeV of deposited energy
    //    Sample actual p.e. from Poisson distribution
    //------------------------------------------------------------
    G4double npe_mean = 40.0 * edep_MeV;
    G4long   npe_sampled = G4Poisson(npe_mean);
    G4double edep_smeared_MeV = (npe_mean > 0) ? edep_MeV * (G4double(npe_sampled) / npe_mean) : 0.0;

    //------------------------------------------------------------
    // 2) Collection efficiency: depth-dependent
    //    efficiency = 1.0 + 0.1 * depth / nhcalLayer
    //    Extract depth3 from the stored detid (lower 3 bits)
    //------------------------------------------------------------
    unsigned int cellid = (*HclCollection)[ij]->GetHitId();
    unsigned int depth3 = cellid & 0x7;
    // Map depth3 back to approximate original layer for efficiency calc
    unsigned int approx_depth = depth3 * 3;
    G4double efficiency = 1.0 + 0.1 * G4double(approx_depth) / G4double(Serc19SimAnalysis::nhcalLayer);
    edep_smeared_MeV *= efficiency;

    //------------------------------------------------------------
    // 3) Electronic noise: Gaussian with sigma = 200 MeV per channel
    //------------------------------------------------------------
    edep_smeared_MeV += G4RandGauss::shoot(0.0, 200.0);

    // Energy in MeV, clamp to non-negative for digitization
    G4int energy_MeV = (edep_smeared_MeV > 0) ? G4int(edep_smeared_MeV) : 0;

    //------------------------------------------------------------
    // 32-bit word packing:
    //   cellid (15 bits: ieta[6]+iphi[6]+depth[3]) << 17 | energy (17 bits)
    //   Energy least count = 1 MeV, max = 2^17 - 1 = 131071 MeV
    //------------------------------------------------------------
    if (energy_MeV > 131071) energy_MeV = 131071;  // cap at 17-bit max

    if (energy_MeV > 0 && pAnalysis->nsimhtHL < pAnalysis->nsimhtmxHL) {
      unsigned long int packed_word = (static_cast<unsigned long int>(cellid) << 17) + energy_MeV;
      pAnalysis->detidHL[pAnalysis->nsimhtHL] = packed_word;
      pAnalysis->timeHL[pAnalysis->nsimhtHL] = edep_smeared_MeV;  // store smeared energy for diagnostics
      pAnalysis->nsimhtHL++;
    }
  }
}

void Serc19HclSD::clear()
{
} 

void Serc19HclSD::DrawAll()
{;
} 

void Serc19HclSD::PrintAll()
{;
} 
