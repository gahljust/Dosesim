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
/// \file dosesimPrimaryGeneratorAction.hh
/// \brief Definition of the dosesimPrimaryGeneratorAction class

#ifndef dosesimPrimaryGeneratorAction_h
#define dosesimPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"
#include "globals.hh"

class G4GeneralParticleSource;
class G4Event;
class G4Box;

/// The primary generator action class with particle gun.
///
/// The default kinematic is a 6 MeV gamma, randomly distribued
/// in front of the phantom across 80% of the (X,Y) phantom size.

class dosesimPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    dosesimPrimaryGeneratorAction();
    virtual ~dosesimPrimaryGeneratorAction();

    // method from the base class
    virtual void GeneratePrimaries(G4Event*);

    // method to access particle gun
    const G4GeneralParticleSource* GetParticleGun() const { return fParticleGun; }

  private:
    G4GeneralParticleSource*  fParticleGun; // pointer a to G4 gun class
    G4Box* fEnvelopeBox;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
