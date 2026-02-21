//
// $Id: Serc19RunAction.cc,v 1.15 2003/11/25 16:50:13 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//#include "Serc19Analysis.hh"

#include "Serc19RunAction.hh"
//#include "g4root.hh"
#include "g4root_defs.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19RunAction::Serc19RunAction() {
  Serc19SimAnalysis *pAnalysis = Serc19SimAnalysis::AnPointer;
  theRunActMessenger = new Serc19RunActionMessenger(this);

  for (int ij=0; ij<10; ij++) {pAnalysis->fNtColId[ij]=-1; }

  SetOutputFile("test");
  SetOutputFile2("g4_test.root");

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19RunAction::~Serc19RunAction() {
  
delete theRunActMessenger; theRunActMessenger=0;
 cout<<"Closing Serc19RunAction"<<endl; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Serc19RunAction::BeginOfRunAction(const G4Run* aRun) {

  Serc19SimAnalysis *pAnalysis = Serc19SimAnalysis::AnPointer;

  G4cout << "### RunID "<< aRun->GetRunID() <<" Nevt " <<aRun->GetNumberOfEvent()<< G4endl;
  
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  //  gRandom->SetSeed(1327512);
  char output_title[100];
  char name[100];
  sprintf(name, r_file_title);
  int ij=aRun->GetRunID();

  sprintf(output_title, "%s_run%d.root", name, ij*10); //body of the file name
  
  pAnalysis->OpenRootfiles(output_title);
  cout <<"output_title = "<<output_title<<endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Serc19RunAction::EndOfRunAction(const G4Run* ) {

  Serc19SimAnalysis *pAnalysis = Serc19SimAnalysis::AnPointer;
  cout<<"Serc19SimAnalysis *pAnalysis = Serc19SimAnalysis::AnPointer;"<<endl;
  pAnalysis->CloseRootfiles();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
