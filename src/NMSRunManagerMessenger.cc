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


#include "NMSRunManagerMessenger.hh"

NMSRunManagerMessenger::NMSRunManagerMessenger(NMSRunManager* newrm) {
  rm = newrm;

  nmsDir = new G4UIdirectory("/NMS/");
  runDir = new G4UIdirectory("/NMS/run/");

  runtimeCmd = new G4UIcmdWithADoubleAndUnit("/NMS/run/runtime", this);
  runtimeCmd->SetGuidance("Set simulated measurement time");
  runtimeCmd->SetParameterName("runtime", false, false);
  runtimeCmd->SetDefaultValue(1);
  runtimeCmd->SetDefaultUnit("s");

  printModuloCmd = new G4UIcmdWithAnInteger("/NMS/run/printmodulo",this);
  printModuloCmd->SetGuidance("Print events modulo n");
  printModuloCmd->SetParameterName("eventmodulo",false);
  printModuloCmd->SetRange("eventmodulo>0");

  runFromTimeCmd = new G4UIcmdWithADouble("/NMS/run/beamOnRuntime", this);
  runFromTimeCmd->SetGuidance("Start calculation with number of events calculated from activity.");
  runFromTimeCmd->SetGuidance("  The activity of NMSMaterialDecaySource is taken.");
  runFromTimeCmd->SetGuidance("  You can also specify a parameter <runs>:");
  runFromTimeCmd->SetGuidance("  0 < runs < 1: Only start a fraction <runs> of the calculated number of events");
  runFromTimeCmd->SetGuidance("  1 <= runs   : Start <runs> different runs with calculated number of events");
  runFromTimeCmd->SetGuidance("                In this case, <runs> will be converted to integer (floor)");

  runFromTimeCmd->SetParameterName("runs", false);
  runFromTimeCmd->SetRange("runs>0");
  runFromTimeCmd->SetDefaultValue(1);

  randomSeedCmd = new G4UIcommand("/NMS/randomSeedList", this);
  randomSeedCmd->SetGuidance("Output some random seeds to file.");

  G4UIparameter * param;
  param = new G4UIparameter("SeedsPerLine",'i',true);
  param->SetParameterRange("SeedsPerLine>0");
  param->SetDefaultValue(2);
  param->SetGuidance("Seeds per line in file");
  randomSeedCmd->SetParameter(param);
  param = new G4UIparameter("SeedLines",'i',true);
  param->SetGuidance("Number of lines");
  param->SetParameterRange("SeedLines>0");
  param->SetDefaultValue(10);
  randomSeedCmd->SetParameter(param);
  param = new G4UIparameter("Filename",'s',true);
  param->SetGuidance("Filename in which seeds will be written");
  param->SetDefaultValue("some_seeds.txt");
  randomSeedCmd->SetParameter(param);
}

NMSRunManagerMessenger::~NMSRunManagerMessenger() {
  delete randomSeedCmd;
  
  delete runFromTimeCmd;
  delete printModuloCmd;
  delete runtimeCmd;

  delete runDir;
  delete nmsDir;
}

void NMSRunManagerMessenger::SetNewValue(G4UIcommand* cmd, G4String newval) {
  if( cmd == runtimeCmd ) {
    rm->SetRuntime( runtimeCmd->GetNewDoubleValue(newval));
  }
  if( cmd == printModuloCmd ) {
    //    rm->GetUserEventAction()->SetPrintModulo(printModuloCmd->GetNewIntValue(newval));
    rm->SetPrintModulo(printModuloCmd->GetNewIntValue(newval));
  }
  if (cmd == runFromTimeCmd) {
    rm->RunFromTime(runFromTimeCmd->GetNewDoubleValue(newval));
  }
  if (cmd == randomSeedCmd) {
    CLHEP::HepRandomEngine* masterRNGEngine = G4Random::getTheEngine();
    std::istringstream is(newval);
    G4int eventsperline;
    G4int lines;
    G4String filename;
    is >> eventsperline >> lines >> filename;
    G4double * randnumbers = new double[eventsperline*lines];
    masterRNGEngine->flatArray(eventsperline*lines,randnumbers); 
    G4RNGHelper* helper = G4RNGHelper::GetInstance();
    helper->Fill(randnumbers,lines,lines,eventsperline);
    std::ofstream outfile;
    outfile.open(filename.c_str());
    for(int i = 0; i < lines; i++) {
      for(int j = 0; j < eventsperline; j++) {
	outfile << helper->GetSeed(i * eventsperline + j);
	if(j + 1 < eventsperline) {
	  outfile << ",";
	}
      }
      outfile << std::endl;
    }
    outfile.close();
    delete randnumbers;
  }
}








