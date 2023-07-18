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
// TETTETRunAction.cc
// \file   MRCP_GEANT4/External/src/TETTETRunAction.cc
// \author Haegin Han
//

#include "G4Timer.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4Neutron.hh"
#include <iostream>
#include "TETRunAction.hh"

TETRunAction::TETRunAction(TETModelImport* _tetData, G4String _output, G4Timer* _init)
:tetData(_tetData), fRun(0), numOfEvent(0), runID(0), outputFile(_output), initTimer(_init), runTimer(0),
 primaryEnergy(-1.), beamArea(-1.), isExternal(true)
{
	if(!isMaster) return;

	runTimer = new G4Timer;
	std::ofstream ofs(outputFile);

    //massMap will be initialized for negative IDs (RBM and BS) in the for loop
	ofs<<"[External: pGycm2 / Internal: SAF (kg-1)]"<<G4endl;
	ofs<<"run#\tnps\tinitT\trunT\tparticle\tsource\tenergy[MeV]\t";
	
	ofs<<"RBM(DRF)\t\tBS(DRF)\t\tRBM\t\tBS\t\t";
	for(size_t i=0;i<tetData->GetDoseName().size();i++)
		ofs<<tetData->GetDoseName()[i]<<"\t"<<tetData->GetDoseMass()[i]/g<<"\t";
	ofs<<"eff. dose (DRF)"<<"\t\t"<< "eff. dose";
	ofs<<G4endl;
	ofs.close();

	std::ofstream ofs2("skin_"+outputFile);
	
    //massMap will be initialized for negative IDs (RBM and BS) in the for loop
	ofs2<<"run#\t";
	for(int i=0;i<tetData->GetNumSkinDet();i++)
		ofs2<<i<<"\t"<<tetData->GetVolume(-i)/cm3<<"\t";
	ofs2<<G4endl;
	ofs2.close();
}

G4Run* TETRunAction::GenerateRun()
{
	// generate run
	fRun = new TETRun(tetData);
	return fRun;
}


void TETRunAction::BeginOfRunAction(const G4Run* aRun)
{
	// print the progress at the interval of 10%
	numOfEvent=aRun->GetNumberOfEventToBeProcessed();
	G4RunManager::GetRunManager()->SetPrintProgress(int(numOfEvent*0.1));
//		    FILE* file = fopen("/proc/self/status", "r");
//		    G4String result;
//		    char line[128];
//
//		    while (fgets(line, 128, file) != NULL){
//		        if (strncmp(line, "VmRSS:", 6) == 0){
//		            result = G4String(line);
//		            break;
//		        }
//		    }
//		    G4cout<<result;
//		    fclose(file);
	if(isMaster){
		initTimer->Stop();
		runTimer->Start();
	}

	const PrimaryGeneratorAction* primary =
			dynamic_cast<const PrimaryGeneratorAction*>(G4RunManager::GetRunManager()
			->GetUserPrimaryGeneratorAction());
	if(!primary) return;
	primaryParticle = primary->GetParticleGun()->GetParticleDefinition()->GetParticleName();
	primarySourceName = primary->GetSourceName();
	primaryEnergy = primary->GetParticleGun()->GetParticleEnergy();
	isExternal = primary-> GetSourceGenerator()->IsExternal();
	if(isExternal) beamArea = primary->GetExternalBeamGenerator()->GetBeamArea();
	fRun->SetPrimary(primaryParticle, primarySourceName, primaryEnergy, beamArea, isExternal);
}

void TETRunAction::EndOfRunAction(const G4Run* aRun)
{
	// print the result only in the Master
	if(!isMaster) return;
	runTimer->Stop();
	// get the run ID
	runID = aRun->GetRunID();

	//get primary info
	primaryParticle = fRun->GetParticleName();
	primarySourceName  = fRun->GetBeamDirName();
	primaryEnergy   = fRun->GetBeamEnergy();
	beamArea        = fRun->GetBeamArea();
	isExternal      = fRun->GetIsExternal();

	G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle(primaryParticle);
	weight = GetRadiationWeighting(particle, primaryEnergy);
	// Print the run result by G4cout and std::ofstream
	//

	// print by G4cout
	if(isExternal) PrintResultExternal(G4cout);
	else           PrintResultInternal(G4cout);

	// print by std::ofstream
	std::ofstream ofs(outputFile.c_str(), std::ios::app);
	if(isExternal)PrintLineExternal(ofs);
	else          PrintLineInternal(ofs);
	ofs.close();

	if(isExternal && tetData->GetNumSkinDet()>0) 
	{
		std::ofstream ofs2("skin_"+outputFile, std::ios::app);
		PrintSkinDoses(ofs2);
		ofs2.close();
	}

	initTimer->Start();
}

void TETRunAction::PrintResultExternal(std::ostream &out)
{
	// Print run result
	//
	using namespace std;

	out << G4endl
	    << "=======================================================================" << G4endl
	    << " Run #" << runID << " / Number of event processed : "<< numOfEvent     << G4endl
	    << "=======================================================================" << G4endl
		<< " Init time: " << initTimer->GetRealElapsed() << " s / Run time: "<< runTimer->GetRealElapsed()<<" s"<< G4endl
	    << "=======================================================================" << G4endl
	    << setw(27) << "organ ID| "
		<< setw(15) << "Organ Mass (g)"
		<< setw(15) << "Dose (pGy*cm2)"
		<< setw(15) << "Relative Error" << G4endl;

	out.precision(3);
	EDEPMAP edepMap = *fRun->GetEdepMap();
	std::vector<G4String> boneNames = {"RBM(DRF)","BS(DRF)","RBM","BS"};
	for(G4int i=0;i<4;i++)
	{
		G4double meanDose = edepMap[-1-i].first / (G4double) numOfEvent;
		G4double squareDoese = edepMap[-1-i].second;
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;
		out << setw(25) << boneNames[i] <<"| ";
		out	<< setw(30) << scientific << meanDose/(joule/kg)*1e12*beamArea/cm2;
		out	<< setw(15) << fixed      << relativeE << G4endl;
	}

	for(size_t i=0;i<tetData->GetDoseName().size();i++){
		G4double meanDose = edepMap[i].first / (G4double) numOfEvent;
		G4double squareDoese = edepMap[i].second;
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;
		out << setw(25) << tetData->GetDoseName()[i] << "| ";
		out	<< setw(15) << fixed      << tetData->GetDoseMass()[i]/g;
        out	<< setw(15) << scientific << meanDose/(joule/kg)*1e12*beamArea/cm2;
		out	<< setw(15) << fixed      << relativeE << G4endl;
	}

	//effective dose
	G4double meanDose = fRun->GetEff_DRF() / (G4double) numOfEvent;
	G4double squareDoese = fRun->GetEff_DRF2();
	G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
	G4double relativeE   = sqrt(variance)/meanDose;
	out << setw(25) << "eff. dose (DRF)" << "| ";
	out	<< setw(15) << " "                ;
    out	<< setw(15) << scientific << meanDose/(joule/kg)*1e12*beamArea/cm2 * weight;
	out	<< setw(15) << fixed      << relativeE << G4endl;

	meanDose = fRun->GetEff() / (G4double) numOfEvent;
	squareDoese = fRun->GetEff2();
	variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
	relativeE   = sqrt(variance)/meanDose;
	out << setw(25) << "eff. dose" << "| ";
	out	<< setw(15) << " "                ;
    out	<< setw(15) << scientific << meanDose/(joule/kg)*1e12*beamArea/cm2 * weight;
	out	<< setw(15) << fixed      << relativeE << G4endl;

	out << "=======================================================================" << G4endl << G4endl;
}

void TETRunAction::PrintResultInternal(std::ostream &out)
{
	// Print run result
	//
	// using namespace std;

	// out << G4endl
	//     << "=======================================================================" << G4endl
	//     << " Run #" << runID << " / Number of event processed : "<< numOfEvent     << G4endl
	//     << "=======================================================================" << G4endl
	// 	<< " Init time: " << initTimer->GetRealElapsed() << " s / Run time: "<< runTimer->GetRealElapsed()<<" s"<< G4endl
	//     << "=======================================================================" << G4endl
	//     << setw(27) << "organ ID| "
	// 	<< setw(15) << "Organ Mass (g)"
	// 	<< setw(15) << "SAF (kg-1)"
	// 	<< setw(15) << "Relative Error" << G4endl;

	// out.precision(3);

	// for(auto itr : massMap){
	// 	if(tetData->DoseWasOrganized()||itr.first<0) out << setw(25) << nameMap[itr.first]<< "| ";
	// 	else                            out << setw(25) << tetData->GetMaterial(itr.first)->GetName()<< "| ";
	// 	out	<< setw(15) << fixed      << itr.second/g;
	// 	out	<< setw(15) << scientific << doses[itr.first].first/primaryEnergy/(1./kg);
	// 	out	<< setw(15) << fixed      << doses[itr.first].second << G4endl;
	// }
	// out << "=======================================================================" << G4endl << G4endl;
}

void TETRunAction::PrintLineExternal(std::ostream &out)
{
	// Print run result
	//
	using namespace std;
	EDEPMAP edepMap = *fRun->GetEdepMap();

	out << runID << "\t" <<numOfEvent<<"\t"<< initTimer->GetRealElapsed() << "\t"<< runTimer->GetRealElapsed()<<"\t"
		<< primaryParticle << "\t" <<primarySourceName<< "\t" << primaryEnergy/MeV << "\t";

	for(G4int i=0;i<4;i++)
	{
		G4double meanDose = edepMap[-1-i].first / (G4double) numOfEvent;
		G4double squareDoese = edepMap[-1-i].second;
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;
		out << meanDose/(joule/kg) * beamArea/cm2 <<"\t" << relativeE << "\t";
	}
	for(size_t i=0;i<tetData->GetDoseName().size();i++)
	{
		G4double meanDose = edepMap[i].first / (G4double) numOfEvent;
		G4double squareDoese = edepMap[i].second;
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;
		out << meanDose/(joule/kg) * beamArea/cm2 <<"\t" << relativeE << "\t";
	}

	G4double meanDose = fRun->GetEff_DRF() / (G4double) numOfEvent;
	G4double squareDoese = fRun->GetEff_DRF2();
	G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
	G4double relativeE   = sqrt(variance)/meanDose;
    out << meanDose/(joule/kg) * beamArea/cm2 *weight<<"\t" << relativeE << "\t";

	meanDose = fRun->GetEff() / (G4double) numOfEvent;
	squareDoese = fRun->GetEff2();
	variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
	relativeE   = sqrt(variance)/meanDose;
    out << meanDose/(joule/kg) * beamArea/cm2 *weight<<"\t" << relativeE << G4endl;
}

void TETRunAction::PrintLineInternal(std::ostream &out)
{
	// // Print run result
	// //
	// using namespace std;
	// EDEPMAP edepMap = *fRun->GetEdepMap();

	// out << runID << "\t" <<numOfEvent<<"\t"<< initTimer->GetRealElapsed() << "\t"<< runTimer->GetRealElapsed()<<"\t"
	// 	<< primaryParticle << "\t" <<primarySourceName<< "\t" << primaryEnergy/MeV << "\t";

	// for(auto itr:doses){
	// 	out << itr.second.first/primaryEnergy/(1./kg) <<"\t" << itr.second.second << "\t";
	// }
	// out<<G4endl;
}

void TETRunAction::PrintSkinDoses(std::ostream &out)
{
	EDEPMAP edepMap = *fRun->GetSkinDoseMap();
	G4double skinDen = tetData->GetMaterial(12200)->GetDensity();
	out<<std::to_string(runID)+"_D"<<"\t";
	for(G4int i=0;i<tetData->GetNumSkinDet();i++)
	{
		G4double meanDose = edepMap[i].first / (G4double) numOfEvent ;
		G4double squareDoese = edepMap[i].second;
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;
		out<<meanDose/(tetData->GetVolume(-i) * skinDen)/(joule/kg)*beamArea/cm2<<"\t"<<relativeE<<"\t";
	}
	out<<G4endl;
	out<<std::to_string(runID)+"_K"<<"\t";
	for(G4int i=tetData->GetNumSkinDet();i<tetData->GetNumSkinDet()*2;i++)
	{
		G4double meanDose = edepMap[i].first / (G4double) numOfEvent ;
		G4double squareDoese = edepMap[i].second;
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;
		out<<meanDose/(tetData->GetVolume(-(i-tetData->GetNumSkinDet())) * skinDen)/(joule/kg)*beamArea/cm2<<"\t"<<relativeE<<"\t";
	}
	out<<G4endl;
}

G4double TETRunAction::GetRadiationWeighting(G4ParticleDefinition* _particle, G4double _energy)
{
	G4double weightingFactor = 1.0; // for Gamma and Electron, Note that Muons need to be considered later.

	if(_particle == G4Proton::Proton()) { //charged pions need to be considered later.
		weightingFactor = 2.0;
	}
	else if(_particle == G4Alpha::Alpha()) { // fission fragments and heavy ions need to be considered later.
		weightingFactor = 20.0;
	}
	else if(_particle == G4Neutron::Neutron()) { //neutron
		if( _energy < 1.0) {// under 1 MeV
			weightingFactor = 2.5 + 18.2*exp(-(pow((log(_energy)),2))/6);
		}
		else if(_energy <= 50) { // From 1 MeV to 50 MeV
			weightingFactor = 5.0 + 17.0*exp(-(pow((log(2.0*_energy)),2))/6);
		}
		else {// more than 50 MeV
			weightingFactor = 2.5 + 3.25*exp(-(pow((log(0.04*_energy)),2))/6);
		}
	}

	return weightingFactor;
}


