//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file dosesimDetectorConstruction.cc
/// \brief Implementation of the dosesimDetectorConstruction class

#include "dosesimDetectorConstruction.hh"
#include "DoseConstants.hh"
#include "G4UserLimits.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVReplica.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "DoseCellParameterisation.hh"
#include "G4PVParameterised.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

dosesimDetectorConstruction::dosesimDetectorConstruction()
: G4VUserDetectorConstruction(),
  fStepLimit(0),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

dosesimDetectorConstruction::~dosesimDetectorConstruction()
{
  delete fStepLimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* dosesimDetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();


/////////........ Build Materials ..........//////////////// /////////........ Build Materials ..........////////////////

G4double a, z, pressure, temperature, fractionmass, density;
G4String name, symbol;
G4int ncomponents;

//Elements
a = 14.01*g/mole;
G4Element* elN = new G4Element(name="Nitrogen",symbol="N" , z= 7., a);
a = 16.00*g/mole;
G4Element* elO = new G4Element(name="Oxygen" ,symbol="O" , z= 8., a);
a=1.01*g/mole;
G4Element* eH = new G4Element(name="Hydrogen", symbol="H", z = 1., a);
a = 39.95*g/mole;
G4Element* eAr = new G4Element(name="Argon", symbol="H", z = 1., a);
a=12.011*g/mole;
G4Element* eC = new G4Element(name="Carbon",symbol="C", z = 12., a);

//Titanium
G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_Ti");

// Aluminum
G4Material* Alu = nist->FindOrBuildMaterial("G4_Al");
G4Material* nyl = nist->FindOrBuildMaterial("G4_NYLON-6-6");

//ABS plastic
G4Material* Abs = new G4Material("ABS", 1.07*g/cm3,3);
Abs->AddElement(eC,15);
Abs->AddElement(eH,17);
Abs->AddElement(elN,1);

// Water
density = 0.006*mg/cm3;// density calculation given below
G4Material* water = new G4Material(name="Water",density,ncomponents=2);
water->AddElement(eH,1);
water->AddElement(elO,2);


temperature= 273.15*kelvin;  // 20 c

pressure =0.853*atmosphere; //0.853*atmosphere;// Pocatello Altitude ~4400 ft -> ~86,383 kPa ~0.853 atm

density = 1.021*mg/cm3;// density at 4400 ft at 86 kPa calculator used from https://www.omnicalculator.com/physics/air-density

                        // Relative humidity is the measure of the percentage of the total allowed amount of
                        //water vapor to the actual amount of water vapor. There is a limit on the amount of
                        //water vapor dependent on the tempuratre. for example in summer maybe the limit is
                        //30g/cm^3, but there is only 10g/cm^3, that is a relative humidity of 33.33%.
                        //At 20 C the limit is 17.3 g/cm^3 Idaho is around 20-40% humidity in summer so im using~ 6 g/m^3.
                        //Water is ~ 17 g/mole -> 0.353 moles/(1/100)^3*cm^3 = 3.53*10^-7 moles/cm^3;
                        //Air is 0.78*14g/mole (N) + 0.21*16g/mole (O)+ 0.01*39.95g/mole (Ar) = 14.68 g/mole
                        // density = 0.00102 g/cm^3 -> 0.000069863 moles/cm^3; Calculating (see airwatercalcs.nb)
                        //Ar=1%,O=20.9%,N=77.6%,water=0.5%
G4Material* water2 = nist->FindOrBuildMaterial("G4_WATER");

//Synthetic human tissue
G4Material* htissue = nist->FindOrBuildMaterial("G4_A-150_TISSUE");

// "Pocatello" Air
G4Material* Air = new G4Material(name="Air ",density,ncomponents=4,kStateGas,temperature,pressure);
Air->AddElement(elN, fractionmass=77.6*perCent);
Air->AddElement(elO, fractionmass=20.9*perCent);
Air->AddElement(eAr,1.*perCent);
Air->AddMaterial(water,fractionmass=0.5*perCent);

//Dry Air near sea level
G4Material* air2 = nist->ConstructNewGasMaterial("air2","G4_AIR",temperature,pressure);


//Vacuum
pressure = 1.e-19*pascal;
temperature = 0.1*kelvin;
G4Material* world_mat = new G4Material("Galactic", z=1., a=1.01*g/mole, 1.e-25*g/cm3,
                   kStateGas,temperature,pressure);

// //OSL mat
G4Material* AlO = nist->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
G4Material* C12 =   new G4Material("C12", 6., 12.011*g/mole , 2.267*g/cm3);
G4Material* osl_mat = new G4Material("Carbon Doped Aluminum Oxide", 3.9698*g/cm3 , 2);
osl_mat->AddMaterial(AlO, 99.988*perCent);
osl_mat->AddMaterial(C12, 0.012*perCent);

//Sensitive volume offset
G4double offset = 0*cm; //+ for closer to the source - for further away



/////////-------------------------------- Air Envelope parameters ----------//////////////// /////////-------------------------------- Air Envelope parameters ----------////////////////

G4double env_sizeXY = 100*cm;
G4double env_sizeZ = env_sizeXY;
G4bool checkOverlaps = true;


/////////-------------------------------- World parameters ----------//////////////// /////////-------------------------------- World material ----------////////////////

G4double world_sizeXY = 1.2*env_sizeXY;
G4double world_sizeZ  = 1.2*env_sizeZ;


/////////-------------------------------- Build World Volume ----------//////////////// /////////-------------------------------- Build World Volume ----------////////////////
  G4Box* solidWorld =
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking



/////////-------------------------------- Build Air Envelope Volume ----------//////////////// /////////-------------------------------- Build Air Envelope Volume ----------////////////////


  G4Box* solidEnv =
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size

  G4LogicalVolume* logicEnv =
    new G4LogicalVolume(solidEnv,            //its solid
                        air2,             //its material
                        "Envelope");         //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking



/////////-------------------------------- BuildTitanium Plate ----------//////////////// /////////-------------------------------- Titanium Plate ----------////////////////


  G4double shape1_xy = 10*cm;
  G4double shape1_dz = .0254*mm;


  G4Box* solidShape1 =
    new G4Box("Shape1",                    //its name
        0.5*shape1_xy, 0.5*shape1_xy, 0.5*shape1_dz);

  G4ThreeVector pos1 = G4ThreeVector(0, 0, -50*cm-0.5*shape1_dz);//-0.5*shape1_dz


  G4LogicalVolume* logicShape1 =
    new G4LogicalVolume(solidShape1,         //its solid
                        shape1_mat,          //its material
                        "Shape1");           //its name

  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    logicShape1,             //its logical volume
                    "Shape1",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking*/




/////////-------------------------------- Aluminum plate and Abs ----------//////////////// /////////-------------------------------- Aluminum plate and Abs ----------////////////////

// G4double alum_xy = 10*cm;
// G4double alum_z = 0.1875*2.54*cm;
// G4Box* aluminumPlate =  new G4Box("Alum",0.5*alum_xy, 0.5*alum_xy, 0.5*alum_z);
// G4ThreeVector pos3 = G4ThreeVector(0, 0,-0.5*alum_z-0.85*mm+offset);//).5*shape2_dr);//-(3*2.54+50*cm);
//
// G4LogicalVolume* logicalum =
//   new G4LogicalVolume(aluminumPlate,Alu,"Aluminum");
//   new G4PVPlacement(0,    //yRot                   //no rotation
//                     pos3,                    //at position
//                     logicalum,             //its logical volume
//                     "Aluminum",                //its name
//                     logicEnv,                //its mother  volume
//                     false,                   //no boolean operation
//                     0,                       //copy number
//                     checkOverlaps);          //overlaps checking

// G4Box* absbox = new G4Box("ABS_Material",0.5*alum_xy, 0.5*alum_xy, 0.5*0.85*mm);
// G4ThreeVector posabs = G4ThreeVector(0,0,-0.5*0.85*mm+offset);
//
// G4LogicalVolume* logicabs =
//   new G4LogicalVolume(absbox,Abs,"ABS_Material");
//   new G4PVPlacement(0,    //yRot                   //no rotation
//                     posabs,                    //at position
//                     logicabs,             //its logical volume
//                     "ABS_Material",                //its name
//                     logicEnv,                //its mother  volume
//                     false,                   //no boolean operation
//                     0,                       //copy number
//                     checkOverlaps);          //overlaps checking
//
// auto absvis = new G4VisAttributes(G4Colour(0.,0.,1.));
// absvis->SetVisibility(true);
// logicabs->SetVisAttributes(absvis);





/////////-------------------------------- Quartz ----------//////////////// /////////-------------------------------- Quartz ----------////////////////

// a=28.085*g/mole;
// G4Element* eSi = new G4Element(name="Silicon",symbol="Si" , z= 14., a);

// G4Material* shape2_mat = new G4Material(name="Quartz",density=2.202*g/cm3,ncomponents=2);
// shape2_mat->AddElement(eSi,1);
// shape2_mat->AddElement(elO,2);

// G4double shape2_dr = 1*cm;
// G4double shape2_xy = 5.5*cm;
// G4double shape2_dz = 1*cm;

G4double shape2_xy = 13*mm;
G4double shape2_dz = 3*mm;

G4Box* solidShape2 =  new G4Box("Shape2",0.5*shape2_xy, 0.5*shape2_xy, 0.5*shape2_dz);

// G4Tubs* solidShape2 = new G4Tubs("Shape2",  0.*cm,shape2_dr, 0.5*shape2_dz, 0.0*deg,360.*deg); //its size



// G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_SILICON_DIOXIDE"); // 2.648 g/cm3


//Dimensions for quartz cylinder peak 5mm need to rotate with yRot


//Dimensions for quartz cylinder full need to rotate with yRot

//Filter full
//G4double shape2_dr = 25*mm;
//G4double shape2_dz = 3*mm;

//Filter peak 5mm
//G4double shape2_dr = 5*mm;
//G4double shape2_dz = 3*mm;


//film peak 5mm
//G4double shape2_dr = 5*mm;
//G4double shape2_dz = .5*mm;


/////////-------------------------------- Glass Plate ----------//////////////// /////////-------------------------------- Glass Plate ----------////////////////


/*

G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_GLASS_PLATE"); // 2.648 g/cm3
G4double shape2_x = 8.5*2.54*cm;
G4double shape2_y = 11*2.54*cm;
G4double shape2_dz = 3*mm;
G4Box* solidShape2 =  new G4Box("Shape2",0.5*shape2_x, 0.5*shape2_y, 0.5*shape2_dz);

*/

/////////-------------------------------- OSL ----------//////////////// /////////-------------------------------- OSL ----------////////////////

// G4double shape2_dr = 0.25*cm;
// G4double shape2_dz = .02505*cm;
// G4Tubs* solidShape2 = new G4Tubs("Shape2",  0.*cm,shape2_dr, 0.5*shape2_dz, 0.0*deg,360.*deg); //its size

// G4ThreeVector pos2 = G4ThreeVector(0, 0, 1.25*shape2_dz+offset);

// /////////-------------------------------- OSL Array ----------//////////////// /////////-------------------------------- OSL Array ----------////////////////
//
// G4Box* oslbox = new G4Box("oslbox", kNofEmRows*10*cm*0.5, kNofEmColumns*10*cm*0.5, shape2_dz*1.25);
//
//
//
// G4LogicalVolume* oslboxlog =
//  new G4LogicalVolume(oslbox,         //its solid
//                      air2,          //its material
//                      "oslbox");           //its name
//  new G4PVPlacement(0,    //yRot                   //no rotation
//      pos2,                    //at position
//      oslboxlog,             //its logical volume
//      "oslbox",                //its name
//      logicEnv,                //its mother  volume
//      false,                   //no boolean operation
//      0,                       //copy number
//      true);          //overlaps checking
//
//      auto visAttributes = new G4VisAttributes(G4Colour(0.,0.,0.));
//      visAttributes->SetVisibility(false);
//      oslboxlog->SetVisAttributes(visAttributes);
     //fVisAttributes.push_back(visAttributes);


//G4double shape2_xy = 20*cm;
//G4double shape2_dz = 0.02*cm;
//G4Box* solidShape2 =  new G4Box("Shape2",0.5*shape2_xy, 0.5*shape2_xy, 0.5*shape2_dz);


/////////-------------------------------- Build Scoring Volume ----------//////////////// /////////-------------------------------- Build Scoring Volume ----------////////////////



G4RotationMatrix* yRot = new G4RotationMatrix;  // Rotates X and Z axes only
yRot->rotateX(90*deg);                     // Rotates 45 degrees


//G4ThreeVector pos2 = G4ThreeVector(0, 0, 0.5*shape2_dz);//).5*shape2_dr);//-(3*2.54+50*cm);


G4LogicalVolume* logicShape2 =
 new G4LogicalVolume(solidShape2,         //its solid
                     //shape2_mat,          //its material
                      nyl,
                     "Shape2");           //its name



// G4VPVParameterisation* shape2parm = new DoseCellParameterisation();
// new G4PVParameterised("OSL_Grid",
//                       logicShape2,        //logical volume
//                       oslboxlog,          //mother volume
//                       kXAxis,             //Placed along this axis
//                       kNofEmCells,        //number of copies
//                       shape2parm,         //the parameterization
//                       checkOverlaps);     //check for overlaps


G4ThreeVector pos2 = G4ThreeVector(0, 0, 0.5*shape2_dz);


 new G4PVPlacement(0,    //yRot                   //no rotation
      pos2,                    //at position
      logicShape2,             //its logical volume
      "Shape2",                //its name
      logicEnv,                //its mother  volume
      false,                   //no boolean operation
      0,                       //copy number
      checkOverlaps);          //overlaps checking


      // G4double maxStep = shape2_dz/5; // 1/5 of Z
      // fStepLimit = new G4UserLimits(maxStep);
      // logicShape2->SetUserLimits(fStepLimit);

  // Set Shape2 as scoring volume
  //
  fScoringVolume = logicShape2;


  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
