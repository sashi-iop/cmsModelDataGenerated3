// $Id: Serc19PrimaryGeneratorAction.cc,v 1.7 2003/09/15 15:38:18 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Serc19PrimaryGeneratorAction.hh"

#include "Serc19DetectorConstruction.hh"
#include "Serc19PrimaryGeneratorMessenger.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"
#include "G4Box.hh"

#include "math.h"
#include "CLHEP/Random/RandGauss.h"
#include "G4UnitsTable.hh"
#include "G4RunManager.hh"

//using namespace std;
//#include "Serc19DetectorParameterDef.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Serc19PrimaryGeneratorAction *Serc19PrimaryGeneratorAction::AnPointer;
Serc19PrimaryGeneratorAction::Serc19PrimaryGeneratorAction(
							 Serc19SimAnalysis *panalysis)
  :pAnalysis(panalysis)
{ 
  AnPointer =this;
  // G4cout<<" Initialized Serc19PrimaryGeneratorAction Constructor"<<endl;
  G4int n_particle = 1;
  particleGun  = new G4ParticleGun(n_particle);
  
  //Default settings :
  
  SetRunNumber(0);
  //  SetOutputFile("simulation.root");
  //  SetFirstEvt(1);

  SetRndmFlag("off");
  SetRndmPIDFlag("off");

  SetPartId(211);

  SetIncEnergy(200.0*GeV);
  SetIncEnergySmr(100*MeV);

  SetIncDirection(G4ThreeVector(1.0,0.0,0.0));
  SetIncThetaSmr(2*mrad);
  SetIncPhiSmr(2*mrad);

  SetIncPosition(G4ThreeVector(0.0*cm,0.0*cm, 0.0*cm));
  SetIncVxSmr(0.03*mm);
  SetIncVySmr(0.03*mm);
  SetIncVzSmr(4.0*cm);
  nMultiplicity = 1;
  initialise = 0;
  
  //create a messenger for this class
  gunMessenger = new Serc19PrimaryGeneratorMessenger(this);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19PrimaryGeneratorAction::~Serc19PrimaryGeneratorAction() {
  cout<<"Closing Serc19PrimaryGeneratorAction"<<endl;  
  if (particleGun)	{delete particleGun;}
  if (gunMessenger) {delete gunMessenger;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Serc19PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  G4double   vx=0, vy=0, vz=0;
  
  if (initialise==0) {
    gunMessenger = new Serc19PrimaryGeneratorMessenger(this);
    g_nevt=-1;
    initialise = 1;
  }
  
  pAnalysis->irun = RunNumber;  //Keep an option that in a file, one may have more than two run number
  g_nevt++;

  pAnalysis->ievt=g_nevt;
  pAnalysis->ngent = (nMultiplicity <=int(pAnalysis->ngenmx)) ? nMultiplicity : pAnalysis->ngenmx;
  cout<<"multi "<< pAnalysis->ngent<<" "<<nMultiplicity<<endl;

  vx = incPosition.x()*mm;
  vy = incPosition.y()*mm;
  vz = incPosition.z()*mm;
    
  if(incVxSmr>=0) {
    vx += G4RandGauss::shoot(0,incVxSmr*mm);
  } else {
    vx += incVxSmr*(2*G4UniformRand()-1.)*mm;
  }

  if(incVySmr>=0) {
    vy += G4RandGauss::shoot(0,incVySmr*mm);
  } else {
    vy += incVySmr*(2*G4UniformRand()-1.)*mm;
  }

  if(incVzSmr>=0) {
    vz += G4RandGauss::shoot(0,incVzSmr*mm);
  } else {
    vz += incVzSmr*(2*G4UniformRand()-1.)*mm;
  }

  double xyrad=pow(vx*vx + vy*vy, 0.5);
  if (xyrad>25) { vx =vx*25.0/xyrad; vy =vy*25.0/xyrad;}

  double tmpthe = (30.0 + 10.0*G4UniformRand());

  double in_Energy = incEnergy*MeV;
  if(incEnergySmr>=0) {
    in_Energy += G4RandGauss::shoot(0,incEnergySmr);
    if (in_Energy <1*MeV) in_Energy=1*MeV;
  } else if(incEnergySmr<0) {
    in_Energy += incEnergySmr*(2*G4UniformRand()-1);
    if (in_Energy <1*MeV) in_Energy=1*MeV;
  }


  int sign=1;
  if (G4UniformRand()>0.5) sign=-1;

  for (unsigned int ij=0; ij<pAnalysis->ngent; ij++)	{
    //Option to have multiple particle
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

    cout<<"partId "<<partId<<endl; 
    G4ParticleDefinition* particle = particleTable->FindParticle(partId*((ij==0) ? sign : -sign));
    particleGun->SetParticleDefinition(particle);
    
    G4ThreeVector ini_Dir(incDirection);
    
    if (rndmFlag=="on") {
      cout<<"ini_Dir "<<tmpthe<<" "<<ini_Dir.theta()<<" "<<(180.0/acos(-1))*ini_Dir.theta()<<endl;
     
     ini_Dir.setPhi(0.0); // (acos(-1.0)/180.)*30*(2*G4UniformRand()-1));
     ini_Dir.setTheta((acos(-1.0)/2.)); // ..*(90.0+30*(2*G4UniformRand()-1)));
         
     in_Energy = 1000.0*(0.25+99.50*G4UniformRand());
     

     cout <<"init "<< in_Energy<<" "<<ini_Dir<<endl;
    }
    
    cout <<"theta "<< ini_Dir<<endl;
    cout <<"theta2 "<< ini_Dir<<endl;
    particleGun->SetParticleMomentumDirection(ini_Dir);
    particleGun->SetParticleMomentum(in_Energy);
    
    particleGun->SetParticlePosition(G4ThreeVector(vx, vy, vz));
    particleGun->GeneratePrimaryVertex(anEvent);
      cout <<"init2 "<< in_Energy<<" "<<ini_Dir<<endl;
    cout<<"Gen Pid "<<particle->GetPDGEncoding()<<" "<<G4BestUnit(in_Energy, "Energy")<<" Dir "<<ini_Dir<<" Pos "<<G4BestUnit(vx, "Length")<<" "<<G4BestUnit(vy, "Length")<<" "<<G4BestUnit(vz, "Length")<<" "<< particleGun->GetParticleMomentumDirection()<<" "<<endl;
    
    if (ij<(int)pAnalysis->ngenmx) { 
      pAnalysis->pidin[ij] = particle->GetPDGEncoding();
      pAnalysis->posxin[ij] = vx/1000.0; //conver from mm to metre
      pAnalysis->posyin[ij] = vy/1000.0;
      pAnalysis->poszin[ij] = vz/1000.0;
      //can not use charge for neutral particle *(particle->GetPDGCharge());
      if (particle->GetPDGCharge()<0) { 
        pAnalysis->momin[ij] = -in_Energy/1000.0; //convert in GeV 
      } else {
        pAnalysis->momin[ij] = in_Energy/1000.0; ////convert in GeV
      }
      pAnalysis->thein[ij] = particleGun->GetParticleMomentumDirection().theta();
      pAnalysis->phiin[ij] = particleGun->GetParticleMomentumDirection().phi();

      cout<<" pAnalysis->thein[ij] "<< pAnalysis->momin[ij]<<" "<<pAnalysis->thein[ij]<<" "<<pAnalysis->phiin[ij]<<endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void Serc19PrimaryGeneratorAction::SetMultiplicity(G4int p) {
  nMultiplicity = p;
}


void Serc19PrimaryGeneratorAction::SetIncEnergy(G4double p) {
  incEnergy = p;
}

void Serc19PrimaryGeneratorAction::SetPartId(G4int p) {
  partId = p;
}


void Serc19PrimaryGeneratorAction::SetIncEnergySmr(G4double p) {
  incEnergySmr = p;
}

void Serc19PrimaryGeneratorAction::SetIncDirection(G4ThreeVector p){
  incDirection = p;
} 

void Serc19PrimaryGeneratorAction::SetIncThetaSmr(G4double p) {
  incThetaSmr =p;
}

void Serc19PrimaryGeneratorAction::SetIncPhiSmr(G4double p) {
  incPhiSmr =p;
}

void Serc19PrimaryGeneratorAction::SetIncPosition(G4ThreeVector p){
  incPosition = p;
} 

void Serc19PrimaryGeneratorAction::SetIncVxSmr(G4double p) {
  incVxSmr =p;
}

void Serc19PrimaryGeneratorAction::SetIncVySmr(G4double p) {
  incVySmr =p;
}

void Serc19PrimaryGeneratorAction::SetIncVzSmr(G4double p) {
  incVzSmr =p;
}

//void Serc19PrimaryGeneratorAction::SetOutputFile(G4String p) {
//  OutputFile = p; pAnalysis->text_outputFile=p;
//}
