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
// TETActionInitialization.cc
// 
// Author: Haegin Han
// Reference: ICRP Publication 145. Ann. ICRP 49(3), 2020.
// Geant4 Contributors: J. Allison and S. Guatelli
//

#include "TETActionInitialization.hh"

TETActionInitialization::TETActionInitialization(TETModelImport* _tetData, G4String _output, G4Timer* _initT)
 : G4VUserActionInitialization(), fTetData(_tetData), fOutput(_output), initT(_initT)
{}

void TETActionInitialization::BuildForMaster() const
{
 SetUserAction(new TETRunAction(fTetData, fOutput, initT));
}

void TETActionInitialization::Build() const
{
 // initialise UserAction classes
 SetUserAction(new PrimaryGeneratorAction(fTetData));
 SetUserAction(new TETRunAction(fTetData, fOutput, initT));
 SetUserAction(new TETSteppingAction());
}  

