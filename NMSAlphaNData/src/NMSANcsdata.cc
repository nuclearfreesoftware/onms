/*
 * ONMS - Open Neutron Multiplicity Simulation
 *
 *
 * Copyright (C) 2013-2016 Moritz Kütt
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact: moritz@nuclearfreesoftware.org
 */

#include "NMSANcsdata.hh"

#include "G4SystemOfUnits.hh"
#include "G4AutoLock.hh"

namespace
{
  G4Mutex singMutex = G4MUTEX_INITIALIZER; //Protects singleton access
}

const G4int NMSANcsdata::csno = 17;
const G4int NMSANcsdata::csarray[csno][2] = {{3, 6}, {3, 7},
					    {4, 9},
					    {5, 10}, {5, 11},
					    {6, 12}, {6, 13},
					    {7, 14}, {7, 15},
					    {8, 17}, {8, 18},
					    {9, 19},
					    {11, 23},
					    {13, 27},
					    {14, 28}, {14, 29}, {14, 30}
};
const G4String NMSANcsdata::nameString[100] = {"Hydrogen", "Helium",
 "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen", "Fluorine",
 "Neon", "Sodium", "Magnesium", "Aluminum", "Silicon", "Phosphorous", 
 "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium", "Scandium",
 "Titanium", "Vanadium", "Chromium", "Manganese", "Iron", "Cobalt", "Nickel",
 "Copper", "Zinc", "Gallium", "Germanium", "Arsenic", "Selenium", "Bromine",
 "Krypton", "Rubidium", "Strontium", "Yttrium", "Zirconium", "Niobium",
 "Molybdenum", "Technetium", "Ruthenium", "Rhodium", "Palladium", "Silver",
 "Cadmium", "Indium", "Tin", "Antimony", "Tellurium", "Iodine", "Xenon",
 "Cesium", "Barium", "Lanthanum", "Cerium", "Praseodymium", "Neodymium",
 "Promethium", "Samarium", "Europium", "Gadolinium", "Terbium", "Dysprosium",
 "Holmium", "Erbium", "Thulium", "Ytterbium", "Lutetium", "Hafnium",
 "Tantalum", "Tungsten", "Rhenium", "Osmium", "Iridium", "Platinium", "Gold",
 "Mercury", "Thallium", "Lead", "Bismuth", "Polonium", "Astatine", "Radon", 
 "Francium", "Radium", "Actinium", "Thorium", "Protactinium", "Uranium", 
 "Neptunium", "Plutonium", "Americium", "Curium", "Berkelium", "Californium",
 "Einsteinium","Fermium"};

NMSANcsdata *NMSANcsdata::Instance() {
  G4AutoLock lock(&singMutex);
  static NMSANcsdata instance;
  return &instance;
}

void NMSANcsdata::Init() {
    // Load cs data
  if(!getenv("NMSALPHALEDATA")) 
    G4Exception("NMScsdata::init", "No CS found", FatalException, "Please setenv NMSALPHALEDATA to point to the alpha n cross-section files.");
  G4String dirName = getenv("NMSALPHALEDATA");
  G4String filename = "";

    for(  std::vector<G4LPhysicsFreeVector *>::iterator it = csIsoVector.begin() ; it != csIsoVector.end(); ++it) {
    delete(*it);
  }
  csIsoVector.clear();
  
  for(G4int i = 0; i < csno; i++) {
  //for(G4int i = 0; i < 2; i++) {
    G4int Z = csarray[i][0];
    G4int A = csarray[i][1];
    std::stringstream ss;
    ss << dirName << "/Inelastic/CrossSection/" << Z << "_" << A << "_" << nameString[Z - 1];
    ss >> filename;
    ss.clear();
    //	filename = dirName + "/Inelastic/CrossSection/" + std::to_string(Z) + "_" + std::to_string(A) + (*G4Element::GetElementTable())[elI]->GetName();
    if(verboseLevel >= 1) {
      G4cout << "NMSAlphaLE: Try to read from " << filename << G4endl;
    }
    std::ifstream filestream( filename , std::ios::in );
    if(filestream.good()) {
      G4String dummy;
      filestream >> dummy; // Set
      filestream >> dummy; // Source
      G4int nData;
      filestream >> nData;
      if(verboseLevel >= 2) {
	G4cout << "NMSAlphaLE: File has " << nData << " entries" << G4endl;
      }
      G4LPhysicsFreeVector * temp = new G4LPhysicsFreeVector(nData, 0 * MeV, 15 * MeV);
      G4double x,y;

      for (G4int iDat=0;iDat<nData;iDat++) {
	if(filestream.good()) {
	  filestream >> x >> y;
	  temp->PutValues(iDat, x * eV, y * barn);
	}
      }
      csIsoVector.push_back(temp);
      // read into pv + add to data
      filestream.close();
    }
    else {
      if(verboseLevel >= 1) {
	G4cout << "NMSAlphaLE: Could not find file with data, using GG-cs" << G4endl;
      }
    }

    std::stringstream ss2;
    ss2 << dirName << "/Inelastic/F01/" << Z << "_" << A << "_" << nameString[Z - 1];
    ss2 >> filename;
    ss2.clear();
    if(verboseLevel >= 1) {
      G4cout << "NMSAlphaLE: Try to read from " << filename << G4endl;
    }
    std::vector< inelasticdata > idv;
    inelasticdata dat;
    std::ifstream fs( filename , std::ios::in );
    if(fs.good()) {
      G4int count = 0;
      G4int mt = 0;
      while(fs.good() and mt != 91) {
	count++;
	G4String dummy;
	fs >> dummy; // Set
	fs >> dummy; // Source
	G4int Qex;
	G4int nData;
	fs >> mt >> dummy;
	fs >> Qex >> dummy >> nData;
	if(verboseLevel >= 1) {
	  G4cout << "qvalue excited state: "<< Qex << G4endl;
	  G4cout << "entries : "<< nData << G4endl;
	  G4cout << "mt: " << mt << G4endl;
	}
	dat.Qex = Qex * eV;
	dat.mt = mt;
	dat.cs = new G4LPhysicsFreeVector(nData, 0 * MeV, 15 * MeV);
	G4double x,y;
	for (G4int iDat=0;iDat<nData;iDat++) {
	  if(fs.good()) {
	    fs >> x >> y;
	    dat.cs->PutValues(iDat, x * eV, y * barn);
	  }
	}
	idv.push_back(dat);
      }

    }
    else {
      if(verboseLevel >= 1) {
	G4cout << "NMSANcsdata: Could not find file with data" << G4endl;
      }
    }
    inelasticCsVector.push_back(idv);

  }

}

G4bool NMSANcsdata::has(G4int Z, G4int A) {
  for(G4int k = 0; k < csno; k++) {
    if(csarray[k][0] == Z && csarray[k][1] == A) {
      return true;
    }
  }
  return false;
}

G4int NMSANcsdata::index(G4int Z, G4int A) {
  for(G4int k = 0; k < csno; k++) {
    if(csarray[k][0] == Z && csarray[k][1] == A) {
      return k;
    }
  }
  return -1;
}

G4LPhysicsFreeVector * NMSANcsdata::getcsVector(G4int idx) {
  if(idx >= 0 && idx < csno) {
    return csIsoVector[idx];
  }
  return 0;
}

G4int NMSANcsdata::inelasticNo(G4int idx) {
  if(idx >= 0 && idx < csno) {
    return inelasticCsVector[idx].size();
  }
  return 0;
}

G4int NMSANcsdata::getInelasticCsMT(G4int idx, G4int ridx) {
  if(idx >= 0 && idx < csno) {
    if(ridx >= 0 && ridx < inelasticCsVector[idx].size()) {
      return inelasticCsVector[idx][ridx].mt;
    }
  }
  return -1;
}

G4LPhysicsFreeVector * NMSANcsdata::getInelasticCsVector(G4int idx, G4int ridx) {
  if(idx >= 0 && idx < csno) {
    if(ridx >= 0 && ridx < inelasticCsVector[idx].size()) {
      return inelasticCsVector[idx][ridx].cs;
    }
  }
  return 0;
}

G4double NMSANcsdata::getInelasticCsQex(G4int idx, G4int ridx) {
  if(idx >= 0 && idx < csno) {
    if(ridx >= 0 && ridx < inelasticCsVector[idx].size()) {
      return inelasticCsVector[idx][ridx].Qex;
    }
  }
  return 0;
}

NMSANcsdata::NMSANcsdata() {
  verboseLevel = 0;
}

NMSANcsdata::~NMSANcsdata() {
  
}
