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

#ifndef NMSAnalysisManager_h
#define NMSAnalysisManager_h 1

enum RunMode { RUNMODE_ALPHAN, RUNMODE_HE_DETECTOR };
enum SDTOutput { SDTperSecond, SDTtotal };
enum NeutronEfficiencyOutput { EFFICIENCY_SOURCE, EFFICIENCY_TOTAL };

#include "NMSDetectedEventSet.hh"
#include "NMSPulsetrainManager.hh"
#include "NMSRunManager.hh"
#include "NMSTrackingAction.hh"
class NMSTrackingAction;

#include "NMSAnalysisManagerMessenger.hh"
class NMSAnalysisManagerMessenger;

class NMSAnalysisManager {
public:
  static NMSAnalysisManager* GetInstance();

  ~NMSAnalysisManager();

  void reset();

  inline void addSourceNeutron(G4double &energy) {
    sourceNeutronCount++;
    totalSourceNeutronEnergy += energy;
    sourceNeutronSpectrum.push_back(energy);
  };
  inline void addSourceAlpha(G4double &energy)  {
    sourceAlphaCount++;
    totalSourceAlphaEnergy += energy;
  };
  inline void addSourceGamma(G4double &energy) {
    sourceGammaCount++;
    totalSourceGammaEnergy += energy;
  };
  inline void addSourceSF() {
    sourceSFcount++;
  };
  inline void addSecondaryNeutron(G4double &energy) {
    secondaryNeutronCount++;
    totalSecondaryNeutronEnergy += energy;
  };

  inline void addLostNeutron(NMSDetectedEvent de) {
    lostneutronlist.push_back(de);
  };

  inline void addDetectorNeutron(NMSDetectedEvent de) {
    nmspm->add(de); 
  }

  inline void addAlphaNReaction(NMSAlphaNReaction ar) {
    alphanset.push_back(ar);
  }

  inline void addSomeNeutron(NMSDetectedEvent de) {
    anyneutronlist.push_back(de);
  }
  void saveAlphaNToFile(G4String filename = "std.alphan");
  //    - (DumpValues)
  //    - dumptofile
  void saveNeutronSpectrumToFile(G4String filename = "std.sns");

  G4int getSourceNeutrons() { return sourceNeutronCount; };
  G4double getSourceNeutronsEnergy();
  G4double getAverageSourceNeutronsPerTime();

  G4int getSourceAlphas() { return sourceAlphaCount; };
  G4double getSourceAlphasEnergy();
  G4int getSourceGammas() { return sourceGammaCount; };
  G4double getSourceGammasEnergy();
  G4int getSourceSF() { return sourceSFcount; };
  G4double getAverageNu();
  G4int getSecondaryNeutrons() { return secondaryNeutronCount; };
  G4double getSecondaryNeutronsEnergy();

  G4double getDetectorNeutronEfficiency( NeutronEfficiencyOutput neo = EFFICIENCY_SOURCE );
  G4double getLostNeutrons() { return lostneutronlist.size(); }
  
  G4double getSingles(SDTOutput output = SDTperSecond);
  G4double getDoubles(SDTOutput output = SDTperSecond);
  G4double getTriples(SDTOutput output = SDTperSecond);
  G4double getEffPuMass();
  G4double getTotalPuMass();

  
  // lifetimehistogram
  void calculateResults();
  void showResults();
  void showDetectorResults();
  void exportToFile(G4String filename);
  void exportToFile();
  void exportCheck(G4int = 1);

  void loadPulsetrainFromFile(G4String filename);
  void loadPulsetrainFromShakesFile(G4String filename);
  void addPulsetrainFromFile(G4String filename, G4double offset = 0);
  void splitPulsetrain(G4String filename, G4int splits = 10);
  
  void SetDetectorVolume(G4String);
  inline G4String GetDetectorVolume() { return detectorVolume; };
  void SetLostNeutronVolume(G4String);
  inline G4String GetLostNeutronVolume() { return lostNeutronVolume; };

  void SetRunMode(RunMode newrm);
  RunMode GetRunMode() { return rm; };

  NMSPulsetrainManager* GetPulsetrainManager() {return nmspm; };

  void SetPulsetrainAnalysisMode(G4int pam);
  G4int GetPulsetrainAnalysisMode() { return pulsetrainAnalysisMode; };

  void SetRuntimeFromPulsetrain(G4bool pfp);
  G4bool GetRuntimeFromPulsetrain() { return runtimeFromPulsetrain; };
  
  void SetOutputFilename(G4String filename);
  G4String GetOuputFilename() {return outputfilename; };

  void SetExportAfterRun(G4bool ear) { exportAfterRun = ear; };
  G4bool GetExportAfterRun() { return exportAfterRun; };

  void SetPulsetrainExport(G4bool ep) { exportPulsetrain = ep; };
  G4bool SetPulseTrainExport() {return exportPulsetrain; };

  void SetEventOffset(G4int eo) { eventIDoffset = eo;  };
  G4int GetEventOffset() { return eventIDoffset; }; 

  void SetSettingsChanged() { calculated = false; }
  
private:
  static NMSAnalysisManager* amInstance;
  
  NMSAnalysisManager();

  G4bool calculated;

  G4int sourceNeutronCount;
  G4double totalSourceNeutronEnergy;
  G4int sourceAlphaCount;
  G4double totalSourceAlphaEnergy;
  G4int sourceGammaCount;
  G4double totalSourceGammaEnergy;
  G4int sourceSFcount;
  std::vector <G4double> sourceNeutronSpectrum;

  G4int secondaryNeutronCount;
  G4double totalSecondaryNeutronEnergy;

  NMSDetectedEventSet lostneutronlist;
  NMSDetectedEventSet anyneutronlist;
  NMSPulsetrainManager* nmspm;
  NMSMultiplicityResult nmsmr;

  NMSAlphaNSet alphanset;

  NMSAnalysisManagerMessenger* nmsamm;

  G4String detectorVolume;
  G4String lostNeutronVolume;

  G4int pulsetrainAnalysisMode;
  G4bool runtimeFromPulsetrain;
  G4int splits;
  
  RunMode rm;

  G4String outputfilename;
  G4bool exportAfterRun;
  G4bool exportPulsetrain;

  G4int eventIDoffset;
};

#endif /* NMSAnalysisManager_h */
