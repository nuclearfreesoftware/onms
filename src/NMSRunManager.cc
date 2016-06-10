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


#include "NMSRunManager.hh"

NMSRunManager::NMSRunManager() : G4RunManager() {
  nmsrmm = new NMSRunManagerMessenger(this);
  
  runtime = 1 * s;
  printModulo = 1;
}

NMSRunManager::~NMSRunManager() {
  delete nmsrmm;
}

void NMSRunManager::SetCLOptions(G4bool newgps) {
  gps = newgps;
}

void NMSRunManager::SetRuntime(G4double rt) {
  if(rt > 0) { // FIXME
    runtime = rt;
    if(!gps) { // If we have normal PGA (NMSMaterialDecaySource), not debugging
      NMSPrimaryGeneratorAction* pga = (NMSPrimaryGeneratorAction*) this->GetUserPrimaryGeneratorAction();
      if(pga != 0) {
	pga->SetRuntime(rt);
      }
      else {
	G4cout << "NMSRunManager ERROR: Could not find PrimaryGenerator Action" << G4endl;
	G4cout << "Aborting." << G4endl;
	exit(1);
      }
    }
  }
}

void NMSRunManager::SetPrintModulo(G4int pm) {
  if(pm > 0) {
    printModulo = pm;
    NMSEventAction* ea = (NMSEventAction*) this->GetUserEventAction();
    if(ea != 0) {
      ea->SetPrintModulo(pm);
    }
  }
  else {
    printModulo = 1;
  }
}

void NMSRunManager::RunFromTime(G4double runs) {
  if(!gps) {
    if(ConfirmBeamOnCondition()) {
      RunInitialization();
      ConstructScoringWorlds();
      G4int runsI = 1;
      G4double events = 0;
      NMSPrimaryGeneratorAction* pga = (NMSPrimaryGeneratorAction*) this->GetUserPrimaryGeneratorAction();
      G4double act = pga->GetActivity();

      if(runs > 0 && runs < 1) {
	events = int(std::floor(runtime * act * runs));
	if(G4double(events) != std::floor(runtime * act * runs)) {
	  G4Exception("NMSRunManager::RunFromTime", "Could not calculate event number", FatalException, "This is a problem.");
	}

	G4cout << "Will start a 'fraction' runs" << G4endl;
	G4cout << "Simulated source activity: " << act * s << G4endl;
	G4cout << "Runtime:                   " << runtime / s << G4endl;
	G4cout << "Fraction:                  " << runs << G4endl;
	G4cout << "Events:                    " << events << G4endl;

      }
      if(runs >= 1) {
	events = int(std::floor(runtime * act));
	runsI = int(std::floor(runs));
	G4cout << "Will start single run / multiple runs" << G4endl;
	G4cout << "Simulated source activity: " << act * s << G4endl;
	G4cout << "Runtime:                   " << runtime / s << G4endl;
	G4cout << "Events:                    " << events << G4endl;
	G4cout << "Number of Runs:            " << runsI << G4endl;
      }
      G4cout << "Starting (1.) Run from Runtime" << G4endl;

      DoEventLoop(events);
      RunTermination();

      if(runs > 1) {
	for(G4int i = 2; i <= runsI; i++) {
	  if(ConfirmBeamOnCondition()) {
	    G4cout << "Starting " << i << ". Run from Runtime" << G4endl;
	    RunInitialization();
	    ConstructScoringWorlds();
	    DoEventLoop(events);
	    RunTermination();
	  }
	}
      }
      
      //      BeamOn(events);
    }

  }
  else {
    G4cout << "WARNING! RunFromTime can not be started with this particle source." << G4endl;
  }
}
