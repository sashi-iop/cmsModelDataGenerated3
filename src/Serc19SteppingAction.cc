// $Id: Serc19SteppingAction.cc,v 1.9 2005/02/02 17:11:11 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Serc19SteppingAction.hh"

#include "Serc19EventAction.hh"
#include "G4Track.hh"

////#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19SteppingAction::Serc19SteppingAction(Serc19EventAction* evt)
  : eventaction(evt)
{
  // ical0Field = new Serc19ElectroMagneticField();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19SteppingAction::~Serc19SteppingAction() { 
  cout<<"Closing Serc19SteppingAction "<<endl;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Serc19SteppingAction::UserSteppingAction(const G4Step* aStep) {

  const G4Track* track = aStep->GetTrack();

  G4VPhysicalVolume* volume = track->GetVolume();
  
  G4double edep = aStep->GetTotalEnergyDeposit();

  G4String namex = volume->GetLogicalVolume()->GetName(); // material->GetName();

  G4double stepl = 0.;
  if (track->GetDefinition()->GetPDGCharge() != 0.) {
    stepl = aStep->GetStepLength();
  }

  if (volume->GetName() == "solidTrackSpt") {
    eventaction->Addcal0(edep,stepl);
  } else if (volume->GetName() == "ParfPhysi") {
    eventaction->AddGap(edep,stepl);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......



