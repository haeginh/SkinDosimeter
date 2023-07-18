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
// TETPSSkinDosimeter
#include "TETPSSkinDosimeter.hh"
#include "G4VScoreHistFiller.hh"
#include "G4UnitsTable.hh"
////////////////////////////////////////////////////////////////////////////////
// Description:
//   This is a primitive scorer class for scoring energy deposit.
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 2010-07-22   Introduce Unit specification.
// 2020-09-03   Use G4VPrimitivePlotter and fill 1-D histo of energy deposit (x)
//              vs. track weight (y)                   (Makoto Asai)
//
///////////////////////////////////////////////////////////////////////////////

TETPSSkinDosimeter::TETPSSkinDosimeter(G4String name, TETModelImport* _tetData)
  : G4VPrimitivePlotter(name, 0)
  , HCID(-1)
  , EvtMap(nullptr), tetData(_tetData)
{
  SetUnit("MeV");
  skinDetNum = tetData->GetNumSkinDet();
}

G4bool TETPSSkinDosimeter::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4int index  = GetIndex(aStep);
  if(index<0) return false;
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep == 0.)
    return false;
  G4double wei = aStep->GetPreStepPoint()->GetWeight();  // (Particle Weight)
  G4double edepwei = edep * wei;
  EvtMap->add(index, edepwei); //edep
  
  edep = -aStep->GetDeltaEnergy();
  edepwei = edep * wei;
  EvtMap->add(skinDetNum+index, edepwei); //kerma approx.

  if(!hitIDMap.empty() && hitIDMap.find(index) != hitIDMap.cend())
  {
    auto filler = G4VScoreHistFiller::Instance();
    if(filler == nullptr)
    {
      G4Exception(
        "TETPSSkinDosimeter::ProcessHits", "SCORER0123", JustWarning,
        "G4TScoreHistFiller is not instantiated!! Histogram is not filled.");
    }
    else
    {
      filler->FillH1(hitIDMap[index], edep, wei);
    }
  }

  return true;
}

void TETPSSkinDosimeter::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<G4double>(GetMultiFunctionalDetector()->GetName(),
                                    GetName());
  if(HCID < 0)
  {
    HCID = GetCollectionID(0);
  }
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*) EvtMap);
}

void TETPSSkinDosimeter::clear() { EvtMap->clear(); }

void TETPSSkinDosimeter::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  for(const auto& [copy, edep] : *(EvtMap->GetMap()))
  {
    G4cout << "  copy no.: " << copy
           << "  energy deposit: " << *(edep) / GetUnitValue() << " ["
           << GetUnit() << "]" << G4endl;
  }
}

void TETPSSkinDosimeter::SetUnit(const G4String& unit)
{
  CheckAndSetUnit(unit, "Energy");
}


G4int TETPSSkinDosimeter::GetIndex(G4Step* aStep)
{
  // return the organ ID (= material index)
  G4int copyNo = aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
  G4int idx = tetData->GetMaterialIndex(copyNo);
  if(idx>0) idx = -1;
  else idx = floor(-idx*0.1);

  return idx;
}
