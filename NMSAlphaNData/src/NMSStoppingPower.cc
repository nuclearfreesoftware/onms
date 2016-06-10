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

#include "NMSStoppingPower.hh"

NMSStoppingPower::NMSStoppingPower() {
  ds = STOPPINGPOWER_DS_G4EMCALCULATOR;
  alpha = G4Alpha::AlphaDefinition();
  initialized = false;
}

NMSStoppingPower::~NMSStoppingPower() {
  
}

void NMSStoppingPower::initialize() {

  if(!getenv("NMSZIEGLERDATA")) 
    G4Exception("NMSANYield::initialize", "No Ziegler Table Data", FatalException, "Please setenv NMSZIEGLERDATA to point to the file with Ziegler coefficients (table 4 from 1977 publication).");
  G4String filename = getenv("NMSZIEGLERDATA");

  // load file
  G4int maxlines = 92;
  std::ifstream filestream( filename , std::ios::in );
  if(filestream.good()) {
    for(G4int i = 0; i<maxlines; i++) {
      G4String elementname;
      filestream >> elementname
		 >> zieglerdata[i][0]
		 >> zieglerdata[i][1]
		 >> zieglerdata[i][2]
		 >> zieglerdata[i][3]
		 >> zieglerdata[i][4];
    }
  }
  else {
    G4Exception("NMSANYield::initialize", "Could not open Ziegler Table Data", FatalException, "Please make sure that there is a file at the location described by NMSZIEGLERDATA");
  }

  // ASTAR
  astar = new G4ASTARStopping();
  astar->Initialise();

  // EMCalculator
  emprocessname = G4LossTableManager::Instance()->GetEnergyLossProcess(alpha)->GetProcessName();
  initialized = true;
}
  
G4double NMSStoppingPower::getDEDX(G4double E, G4int Z) {
  G4double notyet = E * Z;
  return 0;
}

G4double NMSStoppingPower::getDEDX(G4double E, G4Material * mat) {
  if(!initialized) {
    initialize();
  }

  if(ds == STOPPINGPOWER_DS_ZIEGLER_T3) {
    return getDEDXfromZieglerT3(E, mat);
  }
  else if(ds == STOPPINGPOWER_DS_ASTAR) {
    if(astar->GetIndex(mat) != -1) {
      return astar->GetElectronicDEDX(mat, E) * mat->GetDensity();
    }
    else {
      G4cout << "Material not in ASTAR database" << G4endl;
      return 0;
    }

  }
  else if(ds == STOPPINGPOWER_DS_ASTARZIEGLER) {
    if(astar->GetIndex(mat) != -1) {
      return astar->GetElectronicDEDX(mat, E) * mat->GetDensity();
    }
    else {
      return getDEDXfromZieglerT3(E, mat);
    }
  }

  // Default: G4EMCALCULATOR
  return emcal.ComputeDEDX(E, alpha, emprocessname, mat);
  
}

G4double NMSStoppingPower::getDEDXfromZieglerT3(G4double E, G4Material* mat) {
    // G4cout << "ZIEGLER" << G4endl;

  G4int elements = mat->GetNumberOfElements();
  const G4double* vecab = mat->GetVecNbOfAtomsPerVolume();
  const G4double * vecww = mat->GetFractionVector();
  G4double tots = 0;
  G4double tots2 = 0;
  for(G4int i = 0; i < elements; i++) {
    // G4cout << vecab[i] << " " << vec2[i] << " " << vecww[i] << G4endl;
    G4double slowel = 0;
    G4double shighel = 0;
    G4double sel = 0;
    const G4Element * el = mat->GetElement(i);
    G4int Z = el->GetZ();
    //G4cout << "Z:" << Z << " E:" << E / MeV << G4endl;

    if(E == 0) {
      return 0;
    }
    if(Z <= 92) {
      // for(G4int p = 0; p < 5; p++)
      // 	G4cout << zieglerdata[Z-1][p] << " ";
      // G4cout << G4endl;
      slowel = zieglerdata[Z-1][0] * G4Pow::GetInstance()->powA(E * 1000, zieglerdata[Z-1][1]);
      shighel = (zieglerdata[Z-1][2] / E)
	* G4Log(1 + (zieglerdata[Z-1][3] / E) + (zieglerdata[Z-1][4] * E));
      sel = 1 / ((1 / slowel) + (1 / shighel));
      tots2 += sel * vecww[i];
      tots += sel * vecab[i] ;
      //      G4cout << slowel << "  " << shighel << G4endl;
    }
    if(Z == 94) { // from LA-8869-MS
      G4double pu[5] = {5.1486, -0.171158, -0.272723, 0.100975, -0.0160365};

      G4double newel = pu[0]
	+ pu[1] * G4Log(E)
	+ pu[2] * G4Pow::GetInstance()->powA(G4Log(E), 2.0)
	+ pu[3] * G4Pow::GetInstance()->powA(G4Log(E), 3.0)
	+ pu[4] * G4Pow::GetInstance()->powA(G4Log(E), 4.0);
      newel = G4Exp(newel);
      
      // Uranium for comparison
      slowel = zieglerdata[Z-3][0] * G4Pow::GetInstance()->powA(E * 1000, zieglerdata[Z-3][1]);
      shighel = (zieglerdata[Z-3][2] / E)
	* G4Log(1 + (zieglerdata[Z-3][3] / E) + (zieglerdata[Z-3][4] * E));
      sel = 1 / ((1 / slowel) + (1 / shighel));

      tots += newel * vecab[i];
    }
  }
  //G4cout << "tots" << tots << G4endl;
  //  tots2 = tots2 * eV*cm2*1.0e-15;
  //  G4cout << tots2 * totab << " " << tots2 * totab * (g / (MeV * cm * cm)) / mat->GetDensity() << G4endl;
  return tots * eV*cm2*1.0e-15; // "Ziegler factor"
}

G4double NMSStoppingPower::getDEDXfromEMCal(G4double E, G4Material * mat) {
  return emcal.ComputeDEDX(E, alpha, emprocessname, mat);
}

G4double NMSStoppingPower::getDEDXfromASTAR(G4double E, G4Material * mat) {
  if(astar->GetIndex(mat) != -1) {
    return astar->GetElectronicDEDX(mat, E) * mat->GetDensity();
  }
  else {
    G4cout << "Material not in ASTAR database" << G4endl;
    return 0;
  }
}

void NMSStoppingPower::setDataSource(StoppingDataSource newds) {
  // add check: data is accessible?
  ds = newds;
  initialized = false;
}
