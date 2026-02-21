#include "Serc19PhysicsList.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include <iomanip>   

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4HadronElasticPhysics.hh"

//#include "G4DataQuestionaire.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"


#include "G4ProcessTable.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "Serc19DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19PhysicsList::Serc19PhysicsList() 
  :G4VModularPhysicsList() {
  G4LossTableManager::Instance();
  
  // default cut value  (1.0mm) 
  // defaultCutValue = 1.0*mm;
  //GMA older version 29th Nov 2024 :  G4DataQuestionaire it(photon, neutron, no, no, no, neutronxs);
  cout << "<< Geant4 Physics List: Serc19PhysicsList " <<endl;

  defaultCutValue = 0.7*mm;
  G4int ver = 1;
  SetVerboseLevel(ver);
  
  // EM Physics
  RegisterPhysics(new G4EmStandardPhysics(ver));
  
  // Synchroton Radiation & GN Physics
  RegisterPhysics(new G4EmExtraPhysics(ver));
  // Decays
  RegisterPhysics(new G4DecayPhysics(ver));
  
  // Hadron physics
  RegisterPhysics(new G4HadronElasticPhysics(ver) ); 
  RegisterPhysics(new G4HadronPhysicsQGSP_BERT_HP(ver));
  
  // Ion Physics
  RegisterPhysics(new G4IonPhysics(ver));
  
}

Serc19PhysicsList::~Serc19PhysicsList() {
  cout<<"Closing Serc19PhysicsList"<<endl;  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void Serc19PhysicsList::SetCuts() {
  if (verboseLevel >1){
    G4cout << "Serc19PhysicsList::SetCuts: default cut length : "
	   << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }
  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  
  SetCutsWithDefault(); 
  
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(1500*eV, 100*TeV); 

}


