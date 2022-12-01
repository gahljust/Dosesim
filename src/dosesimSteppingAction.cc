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
/// \file dosesimSteppingAction.cc
/// \brief Implementation of the dosesimSteppingAction class

#include "dosesimSteppingAction.hh"
#include "dosesimEventAction.hh"
#include "dosesimDetectorConstruction.hh"
#include "dosesimAnalysis.hh"

#include "G4ThreeVector.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4VProcess.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

dosesimSteppingAction::dosesimSteppingAction(dosesimEventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

dosesimSteppingAction::~dosesimSteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void dosesimSteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) {
    const dosesimDetectorConstruction* detectorConstruction
      = static_cast<const dosesimDetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();
  }

  // get volume of the current step
  G4LogicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();

  G4Track* fTrack = step->GetTrack();

  // check if we are in scoring volume
  if (volume != fScoringVolume) return;

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);
  G4int TrackID = fTrack->GetTrackID();
  G4int ParticleID = fTrack->GetDefinition()->GetPDGEncoding();
  G4ThreeVector prePoint = step->GetPreStepPoint() ->GetPosition();
  G4ThreeVector postPoint = step->GetPostStepPoint()->GetPosition();
  G4ThreeVector point = prePoint + G4UniformRand()*(postPoint - prePoint);
  G4double ke = step->GetPostStepPoint()->GetKineticEnergy();
  G4double px = step->GetPostStepPoint()->GetMomentum().x();
  G4double py = step->GetPostStepPoint()->GetMomentum().y();
  G4double pz = step->GetPostStepPoint()->GetMomentum().z();
  //auto physicsProcess = fTrack->GetCreatorProcess()->GetProcessName();
  const G4VProcess* CurrentProcess = step->GetPreStepPoint()->GetProcessDefinedStep();
  const G4VProcess* LastProcess = step->GetPostStepPoint()->GetProcessDefinedStep();
  const G4String & PostStepProcessName = LastProcess->GetProcessName();
  G4double parentId = fTrack->GetParentID();
  G4double charge = fTrack->GetDefinition()->GetPDGCharge();


  G4double xPos = point[0];
  G4double yPos = point[1];
  G4double zPos = point[2];


  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->FillNtupleDColumn(0, edepStep);
  analysisManager->FillNtupleDColumn(1, xPos);
  analysisManager->FillNtupleDColumn(2, yPos);
  analysisManager->FillNtupleDColumn(3, zPos);
  analysisManager->FillNtupleDColumn(4, ParticleID);
  analysisManager->FillNtupleDColumn(5, TrackID);
  analysisManager->FillNtupleDColumn(6, px);
  analysisManager->FillNtupleDColumn(7, py);
  analysisManager->FillNtupleDColumn(8, pz);
  analysisManager->FillNtupleDColumn(9,parentId);
  analysisManager->FillNtupleDColumn(10,charge);
  analysisManager->FillNtupleDColumn(12,ke);
  if (CurrentProcess != 0)
  {const G4String & StepProcessName = CurrentProcess->GetProcessName();
  analysisManager->FillNtupleSColumn(11,StepProcessName);}
  else
  analysisManager->FillNtupleSColumn(11,PostStepProcessName);
  analysisManager->AddNtupleRow();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
