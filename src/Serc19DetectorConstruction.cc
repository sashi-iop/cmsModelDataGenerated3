//////////////////////////////////
//
//  Define "uniformField" to have uniform field and use field value
//  in SetUniformMagField(1*tesla);<-
//  and choose field fdirection in 
//  magField = new G4UniformMagField(G4ThreeVector(fieldValue,0., 0.)); //Field along x-axis
//
///////////////////////
//
// Define "arbField" for arbtraty magnetic field map
// give input field value as a function of x in file_magField.dat, 
// which is read in constructor
// Serc19MagneticField::Serc19MagneticField(const G4String &file_magField)
// Define proper filed value in memberfunction
//void Serc19MagneticField::MagneticField(const double x[3], double B[3]) const
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//#define uniformField
#define arbField 

#include "Serc19DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Trd.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4PVDivision.hh"

#ifdef uniformField
#include "G4UniformMagField.hh"
#endif

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4NistManager.hh"
#include "G4PVParameterised.hh"
#include "G4RotationMatrix.hh"

#include "G4UniformMagField.hh"

#include "G4FieldManager.hh"
#include "Serc19Field.hh"
#define debug

#ifdef debug
#include "G4Timer.hh"
#endif

#include "G4ProductionCuts.hh"
#include "G4MaterialCutsCouple.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19DetectorConstruction::Serc19DetectorConstruction(Serc19SimAnalysis *p)
  : pAnalysis(p),
    Air(0),
    Brass(0),
    SiliconStr(0),
    Scintillator(0),
    G10(0),
    CarbonFRP(0),
    Paraffin(0),
    PbWO4(0),    

    solidWorld(0),
    solidTrack(0),
    solidTrackSpt(0),
    solidEcal(0),
    solidElectronics(0),
    solidHcal(0),
    solidHcalAbs(0),
    solidHcalBox(0),
    solidHcalSci(0),

    logicWorld(0), physiWorld(0),
    logicTrack(0), physiTrack(0),
    logicTrackSpt(0), physiTrackSpt(0),
    logicEcal(0), physiEcal(0),
    logicElectronics(0), physiElectronics(0),
    logicHcal(0), physiHcal(0),
    logicHcalAbs(0), physiHcalAbs(0),
    logicHcalBox(0), physiHcalBox(0),
    logicHcalSci(0), physiHcalSci(0),
    physiHcalSci_div(0), physiParf(0),
    magField(0), 
    visWorld(0),
    visTrack(0),
    visTrackSpt(0),
    visEcal(0),
    visElectronics(0),
    visHcal(0),
    visHcalAbs(0),
    visHcalBox(0),
    visHcalSci(0),
    visParf(0),
    visNull(0),
    fieldMgr(0),
    localfldMgr(0)

{

  // materials
  DefineMaterials();

  SDman = G4SDManager::GetSDMpointer();
  ecalName = "serc19Trk";
  TrkSD = new Serc19TrkSD(ecalName);
  SDman->AddNewDetector(TrkSD);

  aRegion0 = new G4Region("Tracker_block");

  ecalName = "serc19Ecal";
  EcalSD = new Serc19EcalSD(ecalName);
  SDman->AddNewDetector(EcalSD);

  aRegion1 = new G4Region("ECAL_Block");

  ecalName = "serc19Hcl";
  HclSD = new Serc19HclSD(ecalName);
  SDman->AddNewDetector(HclSD);

  aRegion2 = new G4Region("HCAL_Block");


  // create commands for interactive definition of the calorimeter

  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Serc19DetectorConstruction::~Serc19DetectorConstruction() {

  cout<<"xxxxx"<<endl;
  //delete   cal0SD ; cal0SD=0;
}
  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Serc19DetectorConstruction::Construct() {
  
  // SetUniformMagField(0.0*tesla);
  SetUniformMagField(3.8*tesla);
  
  physiWorld = ConstructCalorimeter();

  return physiWorld;

  //return ConstructCalorimeter();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Serc19DetectorConstruction::DefineMaterials()
{ 

	//This function illustrates the possible ways to define materials
  //  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4String name, symbol;             //a=mass of a mole;
  G4double density;      //z=mean number of protons;  
  // n=number of nucleons in an isotope;
  
  G4int ncomponents,natoms;
  G4double fractionmass;
  
  G4NistManager* mat = G4NistManager::Instance();  
  mat->SetVerbose(1);

  // define Elements
  G4Element* C = mat->FindOrBuildElement("C");
  G4Element* H = mat->FindOrBuildElement("H");
  G4Element* O = mat->FindOrBuildElement("O");
  G4Element* Si = mat->FindOrBuildElement("Si");
  G4Element* Cu = mat->FindOrBuildElement("Cu");
  G4Element* Zn = mat->FindOrBuildElement("Zn");
  
  Iron =  mat->FindOrBuildMaterial("G4_Fe");
  Air =  mat->FindOrBuildMaterial("G4_AIR");

  Brass =  new G4Material("Brass", density=8.4*g/cm3, ncomponents=2);
  Brass->AddElement(Cu, fractionmass=70.0*perCent);
  Brass->AddElement(Zn, fractionmass=30.0*perCent);  
  
  SiliconStr = new G4Material("Silicon", density=2.32*g/cm3, ncomponents=2);
  SiliconStr->AddElement(Si, fractionmass=99.99*perCent);
  SiliconStr->AddElement(O, fractionmass=0.01*perCent);

  Scintillator = new G4Material(name="Scintillator", density=1.032*g/cm3, ncomponents=2);
  Scintillator->AddElement(C, 9);
  Scintillator->AddElement(H, 10);

  // //  Scintillator = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");

  G10 = new G4Material("G10", density= 1.09*g/cm3, ncomponents=4);
  
  G10->AddElement(Si, natoms=1);
  G10->AddElement(O , natoms=2);
  G10->AddElement(C , natoms=3);
  G10->AddElement(H , natoms=3);
  
  CarbonFRP = new G4Material("FRPCarbon",density = 2.26*g/cm3,ncomponents=2);
  CarbonFRP->AddElement(C, fractionmass=99.*perCent);
  CarbonFRP->AddElement(H, fractionmass= 1.*perCent);

  Paraffin =  new G4Material("Paraffin",density = 0.9*g/cm3, ncomponents=2);
  Paraffin->AddElement(C, natoms=31);
  Paraffin->AddElement(H, natoms=64);
 
  //  G4Element* elO = new G4Element(name="Oxygen",symbol="O",z=8.,16.00*g/mole);
  G4Element* elW = new G4Element(name="Tungsten",symbol="W",z=74.,183.84*g/mole);
  G4Element* elPb = new G4Element(name="Lead",symbol="Pb",z=82.,207.20*g/mole);
  
  PbWO4 = new G4Material(name="PbWO4",density=8.280*g/cm3, ncomponents=3);
  PbWO4->AddElement(elPb, natoms=1); 
  PbWO4->AddElement(elW, natoms=1);
  PbWO4->AddElement(O, natoms=4);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Serc19DetectorConstruction::ConstructCalorimeter()
{
  // cout <<"G4VPhysicalVolume* Serc19DetectorConstruction::ConstructCalorimeter()"<<endl;
  // Clean old geometry, if any//
  
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();


  //13 layer silicon detector with silicon support structure
  const int nsilayer=Serc19SimAnalysis::nsilayer;
  //  double sirad[Serc19SimAnalysis::nsilayer]={2.9, 6.8, 10.9, 16.0, 23, 32, 41, 50, 60, 70, 80, 90, 100};
  double sirad[Serc19SimAnalysis::nsilayer]={2.9, 6.8, 10.9, 16.0, 22, 29, 37, 45, 52, 59, 66, 73, 80}; 

  double sithickness = 0.3*mm;
  double supportwidth=0.5*mm;
  double epsilon=1.e-6*mm;
  for (int ij=0; ij<nsilayer; ij++) { sirad[ij] *=cm;}

  double pival = acos(-1)*rad;
  //     
  // World
  //
  //  solidWorld = 0; logicWorld=0; physiWorld=0;
  solidWorld = new G4Tubs("WorldSolid", 0.0, 3.2*m, 4.5*m, -pival/2.-5*epsilon, pival+10*epsilon);
  
  logicWorld = new G4LogicalVolume(solidWorld,
                                   Air,
                                   "WorldLogic");
  physiWorld = new G4PVPlacement(0,
                                 G4ThreeVector(0,0,0),
                                 logicWorld,
                                 "WorldPhysi",
                                 0,
                                 false,
                                 0);
  
  for (int ij=0; ij<nsilayer; ij++) { 
    solidTrack = new G4Tubs("Silicontrk", sirad[ij], sirad[ij]+sithickness, 75*cm, -pival/4., pival/2.);
    logicTrack  = new G4LogicalVolume(solidTrack, SiliconStr, "SitrkLog");
    physiTrack  = new G4PVPlacement(0, 
                                    G4ThreeVector(0,0,0), 
                                    logicTrack,
                                    "solidTrack",
                                    logicWorld,
                                    false,
                                    ij);
    logicTrack->SetSensitiveDetector(TrkSD);
    logicTrack->SetRegion(aRegion0);
    aRegion0->AddRootLogicalVolume(logicTrack);
    
  }
  
  
  for (int ij=0; ij<nsilayer; ij++) { 
    solidTrackSpt = new G4Tubs("TrackSpt", sirad[ij]-supportwidth, sirad[ij]-epsilon, 80*cm, -pival/4., pival/2);
    logicTrackSpt  = new G4LogicalVolume(solidTrackSpt, CarbonFRP, "TrackSptLog");
    
    
    physiTrackSpt  = new G4PVPlacement(0, 
                                       G4ThreeVector(0,0,0), 
                                       logicTrackSpt,
                                       "solidTrackSpt",
                                       logicWorld,
                                       false,
                                       nsilayer+ij);
  }
  
  
  //Segmentation fault  
  //Paraffin in between tracker and ECAL
  G4Tubs* solidParf = new G4Tubs("SolParaffin", (sirad[nsilayer-1]+2.0*cm), (sirad[nsilayer-1]+6.0*cm), 1.2*m, -pival/4., pival/2.);
  
  G4LogicalVolume* logicParf = new G4LogicalVolume(solidParf,
                                                   Paraffin,
                                                   "ParfLogic");
  physiParf = new G4PVPlacement(0,
                                G4ThreeVector(0,0,0),
                                logicParf,
                                "ParfPhysi",
                                logicWorld,
                                false,
                                0);  
  
  //  //ECAL : for simplicity used simple sphreical shell, 
  //       but each crystal in different theta/eta should have different trapizium 
  
  solidEcal = new G4Sphere("SolECAL", 123.0*cm, 152.6*cm, -pival/4., pival/2., pival/4., pival/2.);
  logicEcal = new G4LogicalVolume(solidEcal, PbWO4, "EcalLogic");
  physiEcal = new G4PVPlacement(0,
                                G4ThreeVector(0,0,0),
                                logicEcal,
                                "EcalPhysi",
                                logicWorld,
                                false,
                                0); 
  
  logicEcal->SetSensitiveDetector(EcalSD);
  logicEcal->SetRegion(aRegion1);
  aRegion1->AddRootLogicalVolume(logicEcal);  
  
  //Electronics behind ECAL using G10
  //Segmentation fault
  solidElectronics = new G4Sphere("SolElectronics", 155.0*cm, 160*cm, -pival/4., pival/2., pival/4., pival/2.);
  logicElectronics = new G4LogicalVolume( solidElectronics, G10, "ElectronicsLogic");
  physiElectronics = new G4PVPlacement(0,
                                       G4ThreeVector(0,0,0),
                                       logicElectronics,
                                       "ElectronicsPhysi",
                                       logicWorld,
                                       false,
                                       0);
  
  //HCAL : Make volume of Brass and then placed 4mm scintillator inside
  double r1 = 181.0*cm;  //Inner radius
  double r2 = 295.0*cm;
  double dz = 260.0*cm;
  
  int nhclwedge=pAnalysis->nhclwedge; //36
  int nhcalLayer=pAnalysis->nhcalLayer; //17;
  int nhcalEtaDiv=pAnalysis->nhcalEtaDiv; //34;
	
  double active_dr=0.4*cm; //Width of air box is twice the schintillator width
  double passive_dr=5.0*cm; //Thickness of Brass absorber  
  double wedgeangle=5.0*degree; //Width of scintillator supporting wedge
  double scintangle=4.9*degree; //Width of scintillator in phi
	
  localfldMgr = new G4FieldManager(new Serc19Field());
	
  solidHcalBox = new G4Tubs("SolHcalBox", 0.95*r1, 1.05*r2, 4.3*m, -pival/2.-2*epsilon, pival+4*epsilon);
  logicHcalBox = new G4LogicalVolume( solidHcalBox, Brass, "HcalBoxLogic", localfldMgr);
  //  logicHcalBox = new G4LogicalVolume( solidHcalBox, Brass, "HcalBoxLogic");
  //  logicHcalBox->setFieldManager(localFldMgr, true);
  
  
  physiHcalBox = new G4PVPlacement(0,
                                   G4ThreeVector(0,0,0),
                                   logicHcalBox,
                                   "HcalBoxPhysi",
                                   logicWorld,
                                   false,
                                   0);    
  double angle = atan((dz-5*cm)/r1); //Side is also made with 5cm thick iron 
  double dx1 = r1*tan(angle);
  double dx2 = r2*tan(angle);
  double dy1 = r1*tan(scintangle)/2.;
  double dy2 = r2*tan(scintangle)/2.;
  double dr = (r1+r2)/2.; //Middle of the wedge
	
  solidHcalAbs = new G4Trd("SolHcalAbs", dx1, dx2, dy1, dy2, (r2-r1)/2.);
  logicHcalAbs = new G4LogicalVolume(solidHcalAbs, Brass, "HcalAbsLogic");
	
  for (int ij=0; ij<nhclwedge; ij++) {
    double rotang = -pival/2. + (ij+0.5)*wedgeangle;
    G4ThreeVector xAxis(0.0, 0.0, -1.0);    
    G4ThreeVector yAxis(-sin(rotang), cos(rotang), 0.0);
    G4ThreeVector zAxis(cos(rotang), sin(rotang), 0.0);    
    G4RotationMatrix* rot = new G4RotationMatrix();
    rot->rotateAxes(xAxis, yAxis, zAxis);
    rot->invert();
    
    physiHcalAbs = new G4PVPlacement(rot,
                                     G4ThreeVector(dr*cos(rotang), dr*sin(rotang),0),
                                     logicHcalAbs,
                                     "HcalAbsPhysi",
                                     logicHcalBox,
                                     false,
                                     ij);   
  }
  
  for (int ij=0; ij<nhcalLayer; ij++) { 
    double drsc = r1+8.0*cm+passive_dr/2.0+ij*(2*active_dr+passive_dr);
    double dxsc = 0.995*drsc*tan(angle);
    double dysc = 0.95*drsc*tan(scintangle)/2.;
    double tilesize = dxsc/nhcalEtaDiv;
		
    solidHcal = new G4Box("SolHcal", dxsc+epsilon, dysc, active_dr);
    logicHcal = new G4LogicalVolume(solidHcal, Air, "HcalLogic");
		
    solidHcalSci = new G4Box("SolHcalSci", 0.995*dxsc, 0.95*dysc-epsilon, active_dr/2.0);
    logicHcalSci = new G4LogicalVolume(solidHcalSci, Air, "HcalLogicSci");
    
    G4Box* solidHcalSci_div = new G4Box("SolHcalSci_div", tilesize, dysc-epsilon, active_dr/2.0);
    G4LogicalVolume* logicHcalSci_div =  new G4LogicalVolume(solidHcalSci_div, Scintillator, "HcalLogicSci_div");
    physiHcalSci_div = new G4PVDivision("PhyHcalSci_div", logicHcalSci_div, logicHcalSci, kXAxis, nhcalEtaDiv, tilesize);     
    
    physiHcal = new G4PVPlacement(0,
                                  G4ThreeVector(0,0,drsc-dr),
                                  logicHcal,
                                  "HcalPhysi",
                                  logicHcalAbs,
                                  false,
                                  ij);   
    physiHcalSci = new G4PVPlacement(0,
                                     G4ThreeVector(0,0,0),
                                     logicHcalSci,
                                     "HcalPhysi",
                                     logicHcal,
                                     false,
                                     0);     
    logicHcalSci_div->SetVisAttributes(visNull);
    
    logicHcalSci_div->SetSensitiveDetector(HclSD);
    
  }
  
  //   if(!visWorld){visWorld= new G4VisAttributes(true, G4Colour(0.0,0.0,1.0));}//blue
  if(!visTrack){ visTrack= new G4VisAttributes(true, G4Colour(1.0,0.0,0.0));}//red
  visTrack->SetForceAuxEdgeVisible(true);
  
  //   if(!visTrackSpt){ visTrackSpt= new G4VisAttributes(true, G4Colour(1.0,1.0,0.0));}//yellow
  //   //visTrackSpt->SetForceSolid(true);
  if(!visEcal){
    visEcal= new G4VisAttributes(true, G4Colour(0.5,0.5,0.5));
    visEcal->SetForceSolid(true);
  }//gray
  //   if(!visElectronics){visElectronics= new G4VisAttributes(true, G4Colour(0.0,1.0,1.0));}//cyan
  //   if(!visHcal){visHcal= new G4VisAttributes(true, G4Colour(1.0,0.0,1.0));}//cyan
  if(!visHcalAbs){
    visHcalAbs= new G4VisAttributes(true, G4Colour(0.0,1.0,0.0));
    visHcalAbs->SetForceWireframe(1);
    //      visHcalAbs->SetForceSolid(true);
  }//green
  if(!visHcalBox){visHcalBox= new G4VisAttributes(true, G4Colour(0.0,0.0,1.0));}//blue
  //    visHcalBox->SetForceSolid(true);
  
  //   if(!visParf){visParf= new G4VisAttributes(true, G4Colour(1.0,1.0,1.0));}//white
  
  //   //  logicWorld->SetVisAttributes(visWorld);
  //   logicWorld->SetVisAttributes(visNull);
  
  logicTrack->SetVisAttributes(visTrack);
  //   //   logicTrackSpt->SetVisAttributes(visTrackSpt);
  //  logicTrackSpt->SetVisAttributes(visNull);
  logicEcal->SetVisAttributes(visEcal);
  //   //  logicEcal->SetVisAttributes(visNull);
  //   logicElectronics->SetVisAttributes(visNull);
  logicElectronics->SetVisAttributes(visElectronics);
  
  //   logicHcal->SetVisAttributes(visNull);
  logicHcal->SetVisAttributes(visHcal);
  logicHcalAbs->SetVisAttributes(visHcalAbs);
  //    logicHcalAbs->SetVisAttributes(visNull);
  logicHcalBox->SetVisAttributes(visHcalBox);
  //       logicHcalBox->SetVisAttributes(visNull);
  logicHcalSci->SetVisAttributes(visHcalSci);
  //    logicHcalSci->SetVisAttributes(visNull);
  
  
  logicParf->SetVisAttributes(visParf);
  //   logicParf->SetVisAttributes(visNull);
  
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

void Serc19DetectorConstruction::SetUniformMagField(G4double fieldValue) {
  //apply a global uniform magnetic field along Z axis

  fieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  cout<<"fieldvalue "<< fieldValue/tesla<<" Tesla"<<endl;
  if(magField) { delete magField;}		//delete the existing magn field

  magField = new G4UniformMagField(G4ThreeVector(0.0, 0.0, fieldValue)); //Field along z-axis
  fieldMgr->SetDetectorField(magField);
  fieldMgr->CreateChordFinder(magField);

}

