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
//

#ifndef TETPSSkinDosimeter_h
#define TETPSSkinDosimeter_h 1

#include "G4VPrimitivePlotter.hh"
#include "G4THitsMap.hh"
#include "TETModelImport.hh"

////////////////////////////////////////////////////////////////////////////////
// Description:
//   This is a primitive scorer class for scoring energy deposit.
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura
// 2010-07-22   Introduce Unit specification.
//
// 2020-09-03   Use G4VPrimitivePlotter and fill 1-D histo of energy deposit (x)
//              vs. track weight (y)                   (Makoto Asai)
//
///////////////////////////////////////////////////////////////////////////////

class TETPSSkinDosimeter : public G4VPrimitivePlotter
{
 public:
  TETPSSkinDosimeter(G4String name, TETModelImport* _tetData);  // default unit
  ~TETPSSkinDosimeter() override = default;

 public:
  void Initialize(G4HCofThisEvent*) override;
  void clear() override;
  void PrintAll() override;

  virtual void SetUnit(const G4String& unit);

 protected:
  G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
  virtual G4int GetIndex(G4Step*);

 private:
  G4int HCID, skinDetNum;
  G4THitsMap<G4double>* EvtMap;
  TETModelImport* tetData;
};
#endif