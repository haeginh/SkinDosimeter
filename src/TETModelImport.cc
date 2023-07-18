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
// Author: Haegin Han
// Reference: ICRP Publication 145. Ann. ICRP 49(3), 2020.
// Geant4 Contributors: J. Allison and S. Guatelli
//

#include "TETModelImport.hh"

TETModelImport::TETModelImport(G4String _phantomName, G4UIExecutive* ui)
:fPhantomName(_phantomName)
{
  G4cout << "================================================================================"<<G4endl;
  G4cout << "\t" << fPhantomName << " was implemented in this CODE!!   "<< G4endl;
  G4cout << "================================================================================"<<G4endl;

  G4String eleFile      =  fPhantomName + ".ele";
  G4String nodeFile     =  fPhantomName + ".node";
  G4String materialFile =  fPhantomName + ".material";
  G4String boneFile     =  fPhantomName + ".RBMnBS";
	G4String drfFile      =  fPhantomName + ".DRF";
	
  // read dose file (*.dose) -if there is any
	DoseRead();
  // read phantom data files (*. ele, *.node)
  DataRead(eleFile, nodeFile);
  // read material file (*.material)
  MaterialRead(materialFile);
  // read bone file (*.RBMnBS)
	RBMBSRead(boneFile);
	// read bone file (*.DRF)
	DRFRead(drfFile);
  // read colour data file (colour.dat) if this is interactive mode
  if(ui) ColourRead();
  // print the summary of phantom information
  PrintMaterialInfomation();
}

void TETModelImport::DoseRead(){
	//read dose file : PLEASE be careful not to include dose ID 0
	std::ifstream ifs("dose.txt");
	if(!ifs.is_open())
  {
    // exception for the case when there is no *.node file
    G4Exception("TETModelImport::DoseRead","",FatalErrorInArgument,
                G4String("      There is no dose.txt").c_str());
  }
	
	G4String aLine;
  G4int doseID(0);
  G4double effSum(0);
	while(!ifs.eof()){
		getline(ifs, aLine);
		if(aLine.empty()) break;

		std::stringstream ss(aLine);
    double eff;	ss>>eff;
    effSum += eff;
    G4String name; ss>>name; 
    effCoeff.push_back(eff);
		doseName.push_back(name);
		
    G4int organID;
		while(ss>>organID){
			if(organ2dose.find(organID)==organ2dose.end()) organ2dose[organID] = {doseID};
			else	                                         organ2dose[organID].push_back(doseID);
		}
    doseID++;
	}
  doseMass.resize(doseID);
	ifs.close();
  if(fabs(effSum-(1.-0.13)>1e-5))
  {
        // exception for the case when there is no *.node file
    G4Exception("TETModelImport::DoseRead","",FatalErrorInArgument,
                G4String("      effSum is " + std::to_string(effSum+0.13)).c_str());
  }
}

void TETModelImport::DataRead(G4String eleFile, G4String nodeFile)
{
  G4String tempStr;
  G4int tempInt;

  // Read *.node file
  //
  std::ifstream ifpNode(nodeFile);

  if(!ifpNode.is_open())
  {
    // exception for the case when there is no *.node file
    G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
                G4String("      There is no " + nodeFile ).c_str());
  }
  G4cout << "  Opening TETGEN node (vertex points: x y z) file '"
         << nodeFile << "'" <<G4endl;

  G4int numVertex;
  G4double xPos, yPos, zPos;
  G4double xMin(DBL_MAX), yMin(DBL_MAX), zMin(DBL_MAX);
  G4double xMax(DBL_MIN), yMax(DBL_MIN), zMax(DBL_MIN);

  ifpNode >> numVertex >> tempInt >> tempInt >> tempInt;

  for(G4int i=0; i<numVertex; ++i)
  {
    ifpNode >> tempInt >> xPos >> yPos >> zPos;
 
    // set the unit
    xPos*=cm;
    yPos*=cm;
    zPos*=cm;

    // save the node data as the form of std::vector<G4ThreeVector>
    fVertexVector.push_back(G4ThreeVector(xPos, yPos, zPos));

    // to get the information of the bounding box of phantom
    if (xPos < xMin) xMin = xPos;
    if (xPos > xMax) xMax = xPos;
    if (yPos < yMin) yMin = yPos;
    if (yPos > yMax) yMax = yPos;
    if (zPos < zMin) zMin = zPos;
    if (zPos > zMax) zMax = zPos;
  }

  // set the variables for the bounding box and phantom size
  fBoundingBox_Min = G4ThreeVector(xMin,yMin,zMin);
  fBoundingBox_Max = G4ThreeVector(xMax,yMax,zMax);
  fPhantomSize = G4ThreeVector(xMax-xMin,yMax-yMin,zMax-zMin);
  G4ThreeVector center = (fBoundingBox_Min+fBoundingBox_Max)*0.5;

  ifpNode.close();

  // Read *.ele file
  //
  std::ifstream ifpEle(eleFile);

  if(!ifpEle.is_open())
  {
    // exception for the case when there is no *.ele file
    G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
                G4String("      There is no " + eleFile ).c_str());
  }
  G4cout << "  Opening TETGEN elements (tetrahedron with node No.) file '"
         << eleFile << "'" <<G4endl;

  G4int numEle;
  ifpEle >> numEle  >> tempInt >> tempInt;

  for(G4int i=0; i<numEle; ++i)
  {
    ifpEle >> tempInt;
    auto* ele = new G4int[4];
    for(G4int j=0;j<4;j++)
    {
      ifpEle >> tempInt;
      ele[j]=tempInt;
    }
    fEleVector.push_back(ele);
    ifpEle >> tempInt;
    fMaterialVector.push_back(tempInt);

    // save the element (tetrahedron) data as the form of std::vector<G4Tet*>
    fTetVector.push_back(new G4Tet("Tet_Solid",
                                   fVertexVector[ele[0]]-center,
                                   fVertexVector[ele[1]]-center,
                                   fVertexVector[ele[2]]-center,
                                   fVertexVector[ele[3]]-center));

    // calculate the total volume and the number of tetrahedrons for each organ
    int idx = fMaterialVector[i];

    if(idx<=0)
    {
      fVolumeMap[-floor(-idx*0.1)] += fTetVector[i]->GetCubicVolume();
      fNumTetMap[-floor(-idx*0.1)]++;
      if(idx&1) idx = 12201;
      else idx = 12200;
      fVolumeMap[idx] += fTetVector[i]->GetCubicVolume();
      fNumTetMap[idx]++;
    }
    else
    {
      fVolumeMap[idx] += fTetVector[i]->GetCubicVolume();
      fNumTetMap[idx]++;
    }
  }

  ifpEle.close();
}

void TETModelImport::MaterialRead(G4String materialFile)
{
  // Read material file (*.material)
  //
  std::ifstream ifpMat(materialFile);

  if(!ifpMat.is_open())
  {
    // exception for the case when there is no *.material file
    G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
    G4String("      There is no " + materialFile ).c_str());
  }

  G4cout << "  Opening material file '" << materialFile << "'" <<G4endl;

  char read_data[50];
  char* token;
  G4double zaid;
  G4double fraction;
  G4String MaterialName;
  G4double density;
  G4int i=0;

  while(!ifpMat.eof())
  {
    ifpMat >> read_data;                   //ex) 'C' RBM
    ifpMat >> MaterialName;                //ex)  C 'RBM'
    ifpMat >> read_data;
    density = std::atof(read_data);        //ex) 1.30
    ifpMat >> read_data;                   //ex) g/cm3
    ifpMat >> read_data;
    token = std::strtok(read_data,"m");
    G4int matID = std::atoi(token);        //ex) m'10'
    fMaterialIndex.push_back(matID);
    fOrganNameMap[matID]= MaterialName;
    fDensityMap[matID] = density*g/cm3;

    for(i=0 ;  ; ++i)
    {
      ifpMat >> read_data;
      if(std::strcmp(read_data, "C")==0 || ifpMat.eof()) break;

      zaid = (G4int)(std::atoi(read_data)/1000);
      ifpMat >> read_data;
      fraction = -1.0 * std::atof(read_data);
      fMaterialIndexMap[matID].push_back(std::make_pair(G4int(zaid), fraction));
    }
  }
  ifpMat.close();

  // Construct materials for each organ
  //
  G4NistManager* nistManager = G4NistManager::Instance();

  for(i=0;i<(G4int)fMaterialIndex.size();++i)
  {
    G4int idx = fMaterialIndex[i];
    auto* mat = new G4Material(fOrganNameMap[idx], fDensityMap[idx],
                               G4int(fMaterialIndexMap[idx].size()),
                               kStateSolid, NTP_Temperature, STP_Pressure);
    for(G4int j=0;j<G4int(fMaterialIndexMap[idx].size());++j)
    {
      mat->AddElement(nistManager->FindOrBuildElement(fMaterialIndexMap[idx][j].first),
                                                      fMaterialIndexMap[idx][j].second);
    }
    fMaterialMap[idx]=mat;
    fMassMap[idx]=fDensityMap[idx]*fVolumeMap[idx];
  }
}


void TETModelImport::RBMBSRead(G4String bonefile){
	std::ifstream ifs(bonefile);
	if(!ifs.is_open()) {
		// exception for the case when there is no *.material file
		G4Exception("TETModelImport::RBMBSRead","",FatalException,
				G4String("      There is no " + bonefile ).c_str());
		return;
	}
	G4int idx;
	G4double rbm, bs;
	while(ifs>>idx>>rbm>>bs){
        if(rbmRatio.find(idx)!=rbmRatio.end()) {
            G4cerr<<idx<<" is duplicated in RBMBS file.."<<G4endl;
            exit(0);
        }
		rbmRatio[idx]=rbm;
		bsRatio[idx]=bs;
	}
}

void TETModelImport::DRFRead(G4String DRFfile){
	std::ifstream ifp;
	ifp.open(DRFfile.c_str());

	if(!ifp.is_open()) {
		G4Exception("TETModelImport::DRFRead","",FatalException,
				G4String("      There is no " + DRFfile ).c_str());
		return;
	}

	G4int ID;
    G4double DRF;
    while (!ifp.eof()) {
        G4String dump;
        getline(ifp, dump);
        std::stringstream ss(dump); dump.clear();
        ss >> dump;
        if(dump.empty()) continue;
        ID = atoi(dump.c_str());
        if(rbmDRF.find(ID)!=rbmDRF.end()) {
            G4cerr<<ID<<" is duplicated in DRF file.."<<G4endl;
            exit(0);
        }
        rbmDRF[ID]={};
        bsDRF[ID]={};
    	for (int j=0; j<25; j++) {
            ss >> DRF;
    		rbmDRF[ID].push_back(DRF);
    	}
    	for (int j=0; j<25; j++) {
            ss >> DRF;
    		bsDRF[ID].push_back(DRF);
    	}
    }
    ifp.close();
}

void TETModelImport::ColourRead()
{
  // Read colour data file (colour.dat)
  //
  std::ifstream ifpColour("colour.dat");

  if(!ifpColour.is_open())
  {
    // exception for the case when there is no colour.dat file
    G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
                G4String("Colour data file was not found ").c_str());
  }

  G4cout << "  Opening colour data file 'colour.dat'" <<G4endl;

  G4int organID;
  G4double red, green, blue, alpha;
  while( ifpColour >> organID >> red >> green >> blue >> alpha )
  {
    fColourMap[organID] = G4Colour(red, green, blue, alpha);
  }

  ifpColour.close();
}

void TETModelImport::PrintMaterialInfomation()
{
  // Print the overall information for each organ
  //
  G4cout << G4endl
         << std::setw(9)  << "Organ ID"
         << std::setw(11) << "# of Tet"
         << std::setw(11) << "vol [cm3]"
         << std::setw(11) << "d [g/cm3]"
         << std::setw(11) << "mass [g]"
         << "\t" << "organ/tissue"<< G4endl ;

  G4cout << "--------------------------------------------------------------------------------"<<G4endl;

  std::map<G4int, G4Material*>::const_iterator matIter;
  G4cout<<std::setiosflags(std::ios::fixed);
  G4cout.precision(3);
  for(matIter=fMaterialMap.cbegin(); matIter!=fMaterialMap.cend(); ++matIter)
  {
    G4int idx = matIter->first;
    G4cout << std::setw(9)  << idx                          // organ ID
           << std::setw(11) << fNumTetMap[idx]              // # of tetrahedrons
           << std::setw(11) << fVolumeMap[idx]/cm3          // organ volume
           << std::setw(11) << fMaterialMap[idx] ->GetDensity()/(g/cm3)     // organ density
           << std::setw(11) << fMassMap[idx]/g              // organ mass
           << "\t"<<fMaterialMap[idx]->GetName() << G4endl; // organ name
  }
  G4cout << "--------------------------------------------------------------------------------"<<G4endl;
  doseMass.resize(doseName.size());
  for(auto iter:organ2dose)
  {
    for(auto dose:iter.second) doseMass[dose] += fMassMap[iter.first];
  }
  for(size_t i=0;i<doseMass.size();i++)
  {
    G4cout<< std::setw(20) << doseName[i] <<  std::setw(11) << doseMass[i]/g<<" [g]"<<G4endl;
  }
}
