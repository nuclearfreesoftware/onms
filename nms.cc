/*
 * NMS - NeutronMultiplicitySimulation
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


//GEANT4 header files
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4Version.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "QGSP_BIC_HP.hh"
#include "QGSP_BIC.hh"
#include "NMS_QGSP_BIC_HP_Thermal.hh"

//Output tee
#include <boost/iostreams/tee.hpp>
#include <boost/iostreams/stream.hpp>
#include <fstream>
#include <iostream>

//Project header files
#include "NMSDetectorConstruction.hh"
#include "NMSPrimaryGeneratorAction.hh"
#include "NMSPrimaryGeneratorActionGPS.hh"
#include "NMSAnalysisManager.hh"
#include "NMSTrackingAction.hh"
#include "NMSSteppingAction.hh"
#include "NMSEventAction.hh"
#include "NMSRunAction.hh"
#include "NMSRunManager.hh"

typedef boost::iostreams::tee_device<std::ostream, std::ofstream> TeeDevice;
typedef boost::iostreams::stream<TeeDevice> TeeStream;

//Definitions
const G4String LVNAME_DETECTOR_CAVITY = "NMSLVName-DetectorCavity";
const G4String LVNAME_SAMPLE_MOTHER = "NMSLVName-SampleMother";
const G4String GDML_SAMPLE_POSITION = "NMSGDML-SamplePosition";

void usage() {
      G4cout << "Usage:" << G4endl;
      G4cout << "nms [options]" << G4endl;
      G4cout << G4endl;
      G4cout << "Options:" << G4endl;
      G4cout << "   -D <detector>         GDML File with Detector specification (Default: gdml/detector/PSMC.gdml)" << G4endl;
      G4cout << "   -S <sample>           GDML File with Sample specification (Default: gdml/sample/pm1.gdml)" << G4endl;
      G4cout << G4endl;
      G4cout << "   -n                    do NOT use thermal neutron scattering (use standard QGSP_BIC_HP)" << G4endl;
      G4cout << "   -m <macrofilename>    Load <macrofilename> (Default: vis.mac)" << G4endl;
      G4cout << "   -b <macrofilename>    Run <macrofilename> in batch mode" << G4endl;
      G4cout << "   -l <logfile>          Write all output to log file <logfile>" << G4endl;
      G4cout << "   -p <name>             Name for all output files (log, result etc.). Replaces defaults" << G4endl;
      G4cout << "   -r <no1> <no2>        Create file random_seeds.txt with <no2> lines, each line containing <no1> random seeds" << G4endl;
	
      G4cout << G4endl;
      G4cout << "For Debugging/Testing: " << G4endl;
      G4cout << "   -g                    Use GeneralParticleSource instead of NMSMaterialDecaySource" << G4endl;
      G4cout << "   -f                    Use QGSP_BIC Physics List (no High Precision neutron data)" << G4endl;
      G4cout << "                         (no source Messenger!)" << G4endl;
}

G4bool file_exists (const G4String& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}


int main(int argc, char** argv)
{
  // Variables to store command line parameter selections (with default values)
  G4String detectorgdml = "gdml/detector/PSMC.gdml";
  G4String samplegdml = "gdml/samples/pm1.gdml";
  G4String session = "tcsh";
  G4bool batch = false;
  G4String macro = "vis.mac";
  G4bool thermal = true;
  G4bool gps = false;
  G4String filetest;
  G4String outputfilename = "nms-output";
  G4String logfilename = "nms-output.log";
  G4bool nohp = false;

  // Set Environment variable for heavy elements
  setenv("AllowForHeavyElements", "1", 1);

  // Evaluate command line arguments
  for ( int i=1; i<argc; i=i+1 ) {
    if ( std::string(argv[i]) == "-h" || std::string(argv[i]) == "--help") {
      usage();
      exit(0);
    }
    else if ( std::string(argv[i]) == "-n") {
      thermal = false;
    }
    else if ( std::string(argv[i]) == "-p") {
      i++;
      outputfilename = argv[i];
      logfilename = outputfilename;
      logfilename += ".log";
    }
    else if ( std::string(argv[i]) == "-m" ) {
      i++;
      filetest = argv[i];
      if(file_exists(filetest)) {
	macro = filetest;
      }
      else {
	G4cout << "NMS: ERROR - specified macro file does not exist." << G4endl;
	G4cout << "Aborting..." << G4endl;
	exit(1);
      }
    }
    else if ( std::string(argv[i]) == "-u" ) {
      i++;
      filetest = argv[i];
      if(filetest == "qt" || filetest == "tcsh") {
	session = filetest;
      }
      else {
	G4cout << "NMS: ERROR - specified user interface does not exist (only 'qt' or 'tcsh')." << G4endl;
	G4cout << "Aborting..." << G4endl;
	exit(1);
      }
    }
    else if ( std::string(argv[i]) == "-b" ) {
      i++;
      filetest = argv[i];
      if(file_exists(filetest)) {
	batch = true;
	macro = filetest;
      }
      else {
	G4cout << "NMS: ERROR - batch run, specified macro file does not exist." << G4endl;
	G4cout << "Aborting..." << G4endl;
	exit(1);
      }
    }
    else if ( std::string(argv[i]) == "-D" ) {
      i++;
      filetest = argv[i];
      if(file_exists("gdml/detector/" + filetest)) {
	detectorgdml = "gdml/detector/" + filetest;
      }
      else {
	if(file_exists(filetest)) {
	  detectorgdml = filetest;
	}
	else {
	  G4cout << "NMS: ERROR - specified detector file does not exist." << G4endl;
	  G4cout << "Aborting..." << G4endl;
	  exit(1);
	}
      }
    }
    else if ( std::string(argv[i]) == "-S" ) {
      i++;
      filetest = argv[i];
      if(filetest != "none") {
	if(file_exists("gdml/sample/" + filetest)) {
	  samplegdml = "gdml/sample/" + filetest;
	}
	else {
	  if(file_exists(filetest)) {
	    samplegdml = filetest;
	  }
	  else {
	    G4cout << "NMS: ERROR - specified sample file does not exist." << G4endl;
	    G4cout << "Aborting..." << G4endl;
	    exit(1);
	  }
	}
      }
      else {
	samplegdml = filetest;
      }
    }
    else if ( std::string(argv[i]) == "-l" ) {
      i++;
      logfilename = argv[i];
    }

    else if ( std::string(argv[i]) == "-f" ) { 
      nohp = true;
    }

    else if ( std::string(argv[i]) == "-g" ) { 
      gps = true;
    }
    else if ( std::string(argv[i]) == "-r" ) {
      i++;
      G4int eventsperline = std::atoi(argv[i]);
      i++;
      G4int lines = std::atoi(argv[i]);
      CLHEP::HepRandomEngine* masterRNGEngine = G4Random::getTheEngine();
      G4double * randnumbers = new double[eventsperline*lines];
      masterRNGEngine->flatArray(eventsperline*lines,randnumbers); 
      G4RNGHelper* helper = G4RNGHelper::GetInstance();
      helper->Fill(randnumbers,lines,lines,eventsperline);
      std::ofstream outfile;
      outfile.open("random_seeds.txt");
      for(int ii = 0; ii < lines; ii++) {
	for(int j = 0; j < eventsperline; j++) {
	  outfile << helper->GetSeed(ii * eventsperline + j);
	  if(j + 1 < eventsperline) {
	    outfile << ",";
	  }
	}
	outfile << std::endl;
      }
      outfile.close();
      delete randnumbers;
      exit(0);
    }
    else {
      G4cout << "Warning: " << argv[i] << " is not a valid option" << G4endl;
    }
  }

  // Logging
  std::ofstream logFile;
  logFile.open(logfilename, std::ofstream::out | std::ofstream::app);

  std::ostream tmp(std::cout.rdbuf());
  TeeDevice outputDevice(tmp, logFile);
  TeeStream logger(outputDevice);    

  G4cout.rdbuf(logger.rdbuf());

  // Optional feature: Add Date & Time (Version?)
  G4cout << "Start logging to file!" << G4endl;
  G4cout << "#################################################################################" << G4endl;
  G4cout << "Open Neutron Multiplicity Simulation based on GEANT4" << G4endl;
  G4cout << "#################################################################################" << G4endl;

  // Check for NeutronHP / ParticleHP data set
  if(!nohp) {
    if(!getenv("G4NEUTRONHPDATA")) {
      G4Exception("ONMS main()", "No NeutronHP Data", FatalException, "Please setenv G4NEUTRONHPDATA to point to a folder with G4NDL compatible neutron cross section data.");
    }
    G4String nhp = getenv("G4NEUTRONHPDATA");
    if(nhp.find("JEFF31N") == std::string::npos) {
      G4cout << "Warning: Could not find 'JEFF31N' in path for G4NEUTRONHPDATA, however JEFF 3.1 cross sections are recommended." << G4endl;
    }
  }
  
  // Initialise RunManager
  NMSRunManager * runManager = new NMSRunManager();
  runManager->SetCLOptions(gps);

  // Initialise Geometry (DetectorConstruction)
  G4VUserDetectorConstruction* detector = new NMSDetectorConstruction(detectorgdml, samplegdml);
  runManager->SetUserInitialization(detector);

  G4int verbose = 0;
  if(thermal) {
    NMS_QGSP_BIC_HP_Thermal* physlist = new NMS_QGSP_BIC_HP_Thermal(verbose);
    runManager->SetUserInitialization(physlist);
  }
  else {
    if(nohp) {
      QGSP_BIC* physlist = new QGSP_BIC(verbose);
      runManager->SetUserInitialization(physlist);
    }
    else {
      QGSP_BIC_HP* physlist = new QGSP_BIC_HP(verbose);
      runManager->SetUserInitialization(physlist);
    }
  }

  NMSAnalysisManager * analysisManager = NMSAnalysisManager::GetInstance();
  analysisManager->SetOutputFilename(outputfilename);
  NMSTrackingAction * trackingAction = new NMSTrackingAction(analysisManager);
  NMSSteppingAction * steppingAction = new NMSSteppingAction(analysisManager);
  NMSEventAction * eventAction = new NMSEventAction();
  NMSRunAction * runAction = new NMSRunAction();
  if(!gps) {
    NMSPrimaryGeneratorAction* generatorAction = new NMSPrimaryGeneratorAction();
    runManager->SetUserAction(generatorAction);
  }
  else {
    NMSPrimaryGeneratorActionGPS* generatorAction = new NMSPrimaryGeneratorActionGPS();
    runManager->SetUserAction(generatorAction);
  }
  runManager->SetUserAction(trackingAction);
  runManager->SetUserAction(steppingAction);
  runManager->SetUserAction(eventAction);
  runManager->SetUserAction(runAction);

  runManager->Initialize();

  G4VisManager* visManager = new G4VisExecutive();
  visManager->Initialize();

  G4UImanager * UImanager = G4UImanager::GetUIpointer();
  G4UIExecutive * ui = new G4UIExecutive(argc,argv,session);


  UImanager->ApplyCommand("/control/execute " + macro);

  if(!batch) {
    ui->SessionStart();
  }

  delete ui;

  delete runManager;

  return 0;


}
