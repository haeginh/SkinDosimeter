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
// TETRun.cc
//
// Author: Haegin Han
// Reference: ICRP Publication 145. Ann. ICRP 49(3), 2020.
// Geant4 Contributors: J. Allison and S. Guatelli
//

#include "TETRun.hh"

TETRun::TETRun(TETModelImport* tetData)
:G4Run(), eff(0), eff2(0), eff_drf(0), eff_drf2(0), fCollID(-1), fCollID_skin(-1), fCollID_drf(-1)
{
	organ2dose = tetData->GetOrgan2Dose();
	effCoeff = tetData->GetEffCoeff();
	doseMass = tetData->GetDoseMass();
	auto massMap  = tetData->GetMassMap();
	auto rbmRatio = tetData->GetRBMratio();
	auto bsRatio  = tetData->GetBSratio();

	for(auto rbm:rbmRatio)
		rbmFactor[rbm.first] = rbm.second / massMap[rbm.first];
	for(auto bs:bsRatio)
		bsFactor[bs.first] = bs.second / massMap[bs.first];
}

TETRun::~TETRun()
{
 	fEdepMap.clear();
	fSkinDoseMap.clear();
}

void TETRun::RecordEvent(const G4Event* event)
{
	if(fCollID<0)
	{
		fCollID	= G4SDManager::GetSDMpointer()->GetCollectionID("PhantomSD/eDep");
		fCollID_skin	= G4SDManager::GetSDMpointer()->GetCollectionID("PhantomSD/skin");
		fCollID_drf	= G4SDManager::GetSDMpointer()->GetCollectionID("PhantomSD/drf");
	}

	// Hits collections
	//
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
	if(!HCE) return;

	std::map<G4int, G4double> edepSum; // -1: RBM_DRF, -2: BS_DRF, -3: RBM, -4: BS
	auto* evtMap = static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID));
	for (auto itr : *evtMap->GetMap())
	{
		for(auto doseID:organ2dose[itr.first])
			edepSum[doseID]  += *itr.second;
		auto findRBM = rbmFactor.find(itr.first);
		auto findBS = bsFactor.find(itr.first);
		if(findRBM!=rbmFactor.end()) edepSum[-3] += *itr.second * findRBM->second;
		if(findBS!=bsFactor.end()) edepSum[-4] += *itr.second * findBS->second;
	}
	evtMap = static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_drf));
	for (auto itr : *evtMap->GetMap())
	{
		G4double dose = *itr.second * (joule/kg);
		edepSum[-1-itr.first] += dose;                   //sum
		edepSum[-1-itr.first] += dose * dose; //sum square
	}

	for(auto &iter:edepSum)
	{
		if(iter.first<0) continue;
		iter.second = iter.second / doseMass[iter.first];
	}

	for(auto edep:edepSum){
		fEdepMap[edep.first].first += edep.second;                 //sum
		fEdepMap[edep.first].second += edep.second * edep.second;  //square sum
	}
	// effective dose
	double eff0(0), eff_drf0(0);
	for(size_t i=0;i<effCoeff.size();i++)
		eff0+=edepSum[i]*effCoeff[i];
	eff_drf0 = eff0 + 0.12*edepSum[-1] + 0.01*edepSum[-2];
	eff0 += 0.12*edepSum[-3] + 0.01*edepSum[-4];
	eff += eff0;
	eff2 += eff0*eff0;
	eff_drf += eff_drf0;
	eff_drf2 += eff_drf0*eff_drf0;

	evtMap = static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_skin));
	for (auto itr : *evtMap->GetMap())
	{
		fSkinDoseMap[itr.first].first  += *itr.second;                   //sum
		fSkinDoseMap[itr.first].second += (*itr.second) * (*itr.second); //sum square
	}
}

void TETRun::Merge(const G4Run* run)
{
	auto localRun = static_cast<const TETRun*>(run);
 	EDEPMAP localMap = localRun->fEdepMap;
 	for(auto itr : localMap){
	 	fEdepMap[itr.first].first  += itr.second.first;
        fEdepMap[itr.first].second += itr.second.second;
	}

 	localMap = localRun->fSkinDoseMap;
 	for(auto itr : localMap){
	 	fSkinDoseMap[itr.first].first  += itr.second.first;
        fSkinDoseMap[itr.first].second += itr.second.second;
	}

	eff+=localRun->eff;
	eff2+=localRun->eff2;
	eff_drf+=localRun->eff_drf;
	eff_drf2+=localRun->eff_drf2;

	primary = localRun->primary;
	dir = localRun->dir;
	primaryE = localRun->primaryE;
	beamArea = localRun->beamArea;
	isExternal = localRun->isExternal;

	G4Run::Merge(run);
}





