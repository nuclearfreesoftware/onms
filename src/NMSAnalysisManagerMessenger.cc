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

#include "NMSAnalysisManagerMessenger.hh"

NMSAnalysisManagerMessenger::NMSAnalysisManagerMessenger(NMSAnalysisManager* newnmsam) {
  G4String guid;
  
  nmsam = newnmsam;

  analysisDir = new G4UIdirectory("/NMS/analysis/");
  analysisDir->SetGuidance("NMS analysis settings, export/import commands");

  ptdataDir = new G4UIdirectory("/NMS/analysis/pulsetraindata/");
  ptdataDir->SetGuidance("Functions to load, store and manipulate pulsetrain data");
  
  setDetectorVolumeCmd = new G4UIcmdWithAString("/NMS/analysis/detectorvolume", this);
  setDetectorVolumeCmd->SetGuidance("Set volume name of the volume where neutron absorption is added to pulsetrain.");

  setLostNeutronVolumeCmd = new G4UIcmdWithAString("/NMS/analysis/lostvolume", this);
  setLostNeutronVolumeCmd->SetGuidance("Set volume name where neutrons are counted as lost.");

  showResultsCmd = new G4UIcmdWithoutParameter("/NMS/analysis/showresults", this);
  showResultsCmd->SetGuidance("Show results of current run or loaded data.");
  
  saveSourceNeutronSpectrumToFileCmd = new G4UIcmdWithAString("/NMS/analysis/writesourceneutronenergies", this);
  saveSourceNeutronSpectrumToFileCmd->SetGuidance("Write list of source neutron energies to file.");

  eventOffsetCmd = new G4UIcmdWithAnInteger("/NMS/analysis/eventoffset", this);
  eventOffsetCmd->SetGuidance("Set an event-id offset for pulsetrain (e.g. for multiple parallel runs).");
  eventOffsetCmd->SetDefaultValue(0);

  writeResultsToFileCmd = new G4UIcmdWithAString("/NMS/analysis/writeresults", this);
  writeResultsToFileCmd->SetGuidance("Write results to files (detector statistics, multiplicity results and settings)");
  writeResultsToFileCmd->SetGuidance("  If parameter is omitted and the -p option has been");
  writeResultsToFileCmd->SetGuidance("  specified at the start of NMS, that filename is");
  writeResultsToFileCmd->SetGuidance("  used. Otherwise will use a default filename.");
  writeResultsToFileCmd->SetParameterName("outputfilename", true);

  exportResultsAfterRunCmd = new G4UIcmdWithABool("/NMS/analysis/writeresultsafterrun", this);
  exportResultsAfterRunCmd->SetGuidance("Automatically export results after each run."); 
  exportResultsAfterRunCmd->SetGuidance("The filename given with the -p option will be used,");
  exportResultsAfterRunCmd->SetGuidance("the current run no is added");
  exportResultsAfterRunCmd->SetDefaultValue(true);

  exportPulsetrainCmd = new G4UIcmdWithABool("/NMS/analysis/writeincludepulsetrain", this);
  exportPulsetrainCmd->SetGuidance("Turn pulsetrain file export on/off for normal result write.");
  
  loadPulsetrainFromFileCmd = new G4UIcmdWithAString("/NMS/analysis/readpulsetrain", this);
  loadPulsetrainFromFileCmd->SetGuidance("Load a pulsetrain from a pulsetrain file. Do not specify file ending.");

  loadPulsetrainFromShakesFileCmd = new G4UIcmdWithAString("/NMS/analysis/readshakespulsetrain", this);
  loadPulsetrainFromShakesFileCmd->SetGuidance("Load a pulsetrain from a file that contains a list of events in shakes.");

  addPulsetrainFromFileCmd = new G4UIcommand("/NMS/analysis/pulsetraindata/add", this);
  addPulsetrainFromFileCmd->SetGuidance("Add pulsetraindata from file to set in memory.");
  addPulsetrainFromFileCmd->SetGuidance("  (Just loads data if memory is empty)");
  G4UIparameter * param;
  param = new G4UIparameter("filename", 's', false);
  param->SetGuidance("Filename (.pulsetrain will be automatically added)");
  addPulsetrainFromFileCmd->SetParameter(param);
  param = new G4UIparameter("offsettime", 'd', false);
  param->SetGuidance("Time offset to be added to data");
  addPulsetrainFromFileCmd->SetParameter(param);

  splitPulsetrainCmd = new G4UIcommand("/NMS/analysis/pulsetraindata/split", this);
  splitPulsetrainCmd->SetGuidance("Split pulsetrain into smaller time blocks and write to file");
  splitPulsetrainCmd->SetGuidance("Default number of splits is 10");
  splitPulsetrainCmd->SetGuidance("  (Just loads data if memory is empty)");
  param = new G4UIparameter("filename", 's', false);
  param->SetGuidance("Filename (.pulsetrain will be automatically added)");
  splitPulsetrainCmd->SetParameter(param);
  param = new G4UIparameter("splits", 'i', true);
  param->SetGuidance("Number of splits");
  param->SetDefaultValue(10);
  splitPulsetrainCmd->SetParameter(param);

  savePulsetrainCmd = new G4UIcmdWithAString("/NMS/analysis/pulsetraindata/save", this);
  savePulsetrainCmd->SetGuidance("Write pulsetrain to file");
  // sequential split pulsetrain

  resetPulsetrainCmd = new G4UIcmdWithoutParameter("/NMS/analysis/pulsetraindata/reset", this);
  resetPulsetrainCmd->SetGuidance("Delete pulsetraindata from memory.");
  
  recalculateCmd = new G4UIcmdWithoutParameter("/NMS/analysis/recalculate", this);
  recalculateCmd->SetGuidance("Re-do pulsetrain analysis (normally automatic on changed settings)");
  
  setMultiplicityLengthCmd = new G4UIcmdWithAnInteger("/NMS/analysis/multiplicitylength", this);
  setMultiplicityLengthCmd->SetGuidance("Set maximal multiplicity for R / A distributions");

  setRegisterLengthCmd = new G4UIcmdWithAnInteger("/NMS/analysis/registerlength", this);
  setRegisterLengthCmd->SetGuidance("Set length of R / A shift register (in positions!)");
  setRegisterLengthCmd->SetGuidance("  the length in time is calculated by registerlength * registerperiod");
  setRegisterLengthCmd->SetGuidance("Value will also be used for register quantization, if turned on.");
  setRegisterPeriodCmd = new G4UIcmdWithADoubleAndUnit("/NMS/analysis/registerperiod", this);
  setRegisterPeriodCmd->SetGuidance("Set register period of R / A shift register");
  setRegisterPeriodCmd->SetGuidance("  the length in time is calculated by registerlength * registerperiod");
  setRegisterPeriodCmd->SetGuidance("Value will also be used for register quantization, if turned on.");
  setPredelayCmd = new G4UIcmdWithADoubleAndUnit("/NMS/analysis/setpredelay", this);
  setPredelayCmd->SetGuidance("Predelay for pulsetrain analysis");
  setLongdelayCmd = new G4UIcmdWithADoubleAndUnit("/NMS/analysis/setlongdelay", this);
  setLongdelayCmd->SetGuidance("Long delay for pulsetrain analysis");

  setGateCmd = new G4UIcmdWithADoubleAndUnit("/NMS/analysis/setgate", this);
  setGateCmd->SetGuidance("Gate for pulsetrain analysis");
  setGateCmd->SetGuidance("  registerperiod is automatically new calculated as gate / registerlength");

  setDerandomizePeriodCmd = new G4UIcmdWithADoubleAndUnit("/NMS/analysis/setderandomizeperiod", this);
  setDerandomizePeriodCmd->SetGuidance("Set period for derandomize-quantization step.");
  derandomizeDoCmd = new G4UIcmdWithABool("/NMS/analysis/derandomizedo", this);
  derandomizeDoCmd->SetGuidance("Turn on/off derandomize quantization step.");
  registerQuantizationDoCmd = new G4UIcmdWithABool("/NMS/analysis/registerquantizationdo", this);
  registerQuantizationDoCmd->SetGuidance("Turn on/off register quantization step.");

  setEfficiencyCmd = new G4UIcmdWithADouble("/NMS/analysis/efficiency", this);
  setEfficiencyCmd->SetGuidance("Set detector neutron efficiency.");
  dieawayMethodCmd = new G4UIcmdWithAnInteger("/NMS/analysis/dieawaymethod", this);
  guid = "  " + std::to_string(DAM_FIXED) + " - fixed dieaway time";
  dieawayMethodCmd->SetGuidance(guid.c_str());
  guid = "  " + std::to_string(DAM_AVERAGELIFETIME) + " - set dieaway time to average lifetime.";
  dieawayMethodCmd->SetGuidance(guid.c_str());
  dieawayCmd = new G4UIcmdWithADoubleAndUnit("/NMS/analysis/dieaway", this);
  dieawayCmd->SetGuidance("Set detector neutron dieaway time.");
  guid = "Only used, if dieawaymethod is set to " + std::to_string(DAM_FIXED);
  dieawayCmd->SetGuidance(guid);
  infiniteGateCmd = new G4UIcmdWithADoubleAndUnit("/NMS/analysis/infinitegate", this);
  infiniteGateCmd->SetGuidance("Set length for 'infinite' gate to be used for gate fraction calculations.");

  gatefractionMethodCmd = new G4UIcmdWithAnInteger("/NMS/analysis/gatefractionmethod", this);
  gatefractionMethodCmd->SetGuidance("Choose method for calculation of gate fractions");
  guid = "  " + std::to_string(GM_FROM_DIEAWAY) + " - set gate fractions as calculated from dieaway time.";
  gatefractionMethodCmd->SetGuidance(guid.c_str());
  guid = "  " + std::to_string(GM_FROM_LONGGATE) + " - set gate fractions as calculated from 'infinite' gate.";
  gatefractionMethodCmd->SetGuidance(guid.c_str());
  guid = "  " + std::to_string(GM_FROM_EVENTNO) + " - set gate fractions as calculated per eventno method.";
  gatefractionMethodCmd->SetGuidance(guid.c_str());
  guid = "  " + std::to_string(GM_FIXED) + " - set gate fractions to fixed values.";
  gatefractionMethodCmd->SetGuidance(guid.c_str());
  //gatefractionMethodCmd->SetGuidance("  " + std::to_string(GM_FROM_DIEAWAY) + " - set gate fractions as calculated from dieaway time.");
  gateFractionsCmd = new G4UIcommand("/NMS/analysis/gatefractions", this);
  gateFractionsCmd->SetGuidance("Set gate fractions");
  param = new G4UIparameter("doublegatefraction", 'd', false);
  param->SetParameterRange("doublegatefraction>0 && doublegatefraction<=1");
  param->SetGuidance("double gate fraction");
  gateFractionsCmd->SetParameter(param);
  param = new G4UIparameter("triplegatefraction", 'd', false);
  param->SetParameterRange("triplegatefraction>0 && triplegatefraction<=1");
  param->SetGuidance("triple gate fraction");
  gateFractionsCmd->SetParameter(param);
  
  predeadtimeCmd = new G4UIcmdWithADoubleAndUnit("/NMS/analysis/predeadtime", this);
  predeadtimeCmd->SetGuidance("Deadtime for algorithm carried out PRE quantization");
  predeadtimeUpdatingCmd = new G4UIcmdWithABool("/NMS/analysis/predeadtimeupdating", this);
  predeadtimeUpdatingCmd->SetGuidance("Set deadtime to be/not to be updating for algorithm carried out PRE quantization");
  postdeadtimeCmd = new G4UIcmdWithADoubleAndUnit("/NMS/analysis/postdeadtime", this);
  postdeadtimeCmd->SetGuidance("Deadtime for algorithm carried out POST quantization");
  postdeadtimeUpdatingCmd = new G4UIcmdWithABool("/NMS/analysis/postdeadtimeupdating", this);
  postdeadtimeUpdatingCmd->SetGuidance("Set deadtime to be/not to be updating for algorithm carried out POST quantization");
  
  pulsetrainAnalysisModeCmd = new G4UIcmdWithAnInteger("/NMS/analysis/pulsetrainanalysismode", this);
  pulsetrainAnalysisModeCmd->SetGuidance("Set pulsetrain analysis method");
  std::stringstream ss;
  G4String idxstring;
  ss << RESULTS_LL_ALL_EVENTS_NO_SHIFT_BACKWARD << "   backward calculation, all events";
  idxstring = ss.str();
  ss.str("");
  pulsetrainAnalysisModeCmd->SetGuidance(idxstring);
  ss << RESULTS_LL_ALL_EVENTS_NO_SHIFT_FORWARD << "   forward calculation, all events";
  idxstring = ss.str();
  ss.str("");
  pulsetrainAnalysisModeCmd->SetGuidance(idxstring);
  ss << RESULTS_LL_SINGLE_EVENT_NO_SHIFT_BACKWARD << "   forward calculation, single event";
  idxstring = ss.str();
  ss.str("");
  pulsetrainAnalysisModeCmd->SetGuidance(idxstring);
  ss << RESULTS_LL_SINGLE_EVENT_NO_SHIFT_FORWARD << "   forward calculation, single event";
  idxstring = ss.str();
  ss.str("");
  pulsetrainAnalysisModeCmd->SetGuidance(idxstring);
  exportResultsAfterRunCmd->SetDefaultValue(RESULTS_LL_ALL_EVENTS_NO_SHIFT_FORWARD);

  setRuntimeFromPulsetrainCmd = new G4UIcmdWithABool("/NMS/analysis/setruntimefrompulsetrain", this);
  setRuntimeFromPulsetrainCmd->SetGuidance("Set runtime to be last event in currently stored pulsetrain.");

  setRunModeCmd = new G4UIcmdWithAString("/NMS/analysis/setrunmode", this);
  setRunModeCmd->SetGuidance("Set run mode of NMS");
  setRunModeCmd->SetGuidance("Either He-Detector or AlphaN");
  setRunModeCmd->SetDefaultValue("He-Detector");
  setRunModeCmd->SetCandidates("He-Detector AlphaN");
}

NMSAnalysisManagerMessenger::~NMSAnalysisManagerMessenger() {
  delete setRunModeCmd;

  delete setRuntimeFromPulsetrainCmd;
  delete pulsetrainAnalysisModeCmd;
  delete postdeadtimeUpdatingCmd;
  delete postdeadtimeCmd;
  delete predeadtimeUpdatingCmd;
  delete predeadtimeCmd;

  delete gateFractionsCmd;
  delete gatefractionMethodCmd;
  delete infiniteGateCmd;
  delete dieawayCmd;
  delete dieawayMethodCmd;
  delete setEfficiencyCmd;
  
  delete registerQuantizationDoCmd;
  delete derandomizeDoCmd;
  delete setDerandomizePeriodCmd;

  delete setLongdelayCmd;
  delete setPredelayCmd;
  delete setGateCmd;
  delete setRegisterPeriodCmd;
  delete setRegisterLengthCmd;
  delete setMultiplicityLengthCmd;

  delete recalculateCmd;

  delete resetPulsetrainCmd;
  delete savePulsetrainCmd;
  delete splitPulsetrainCmd;
  delete addPulsetrainFromFileCmd;

  delete loadPulsetrainFromShakesFileCmd;
  delete loadPulsetrainFromFileCmd;
  delete exportPulsetrainCmd;
  delete writeResultsToFileCmd;

  delete eventOffsetCmd;

  delete saveSourceNeutronSpectrumToFileCmd;
  
  delete showResultsCmd;

  delete setLostNeutronVolumeCmd;
  delete setDetectorVolumeCmd;

  delete ptdataDir;
  delete analysisDir;
}

void NMSAnalysisManagerMessenger::SetNewValue(G4UIcommand* cmd, G4String newval) {


  if(cmd == setDetectorVolumeCmd) {
    // check if volume exists / some exist
    nmsam->SetDetectorVolume(newval);
  }
  if(cmd == setLostNeutronVolumeCmd) {
    // check if volume exists / some exist
    nmsam->SetLostNeutronVolume(newval);
  }

  if(cmd == showResultsCmd) {
    nmsam->showResults();
  }
  if(cmd == saveSourceNeutronSpectrumToFileCmd) {
    G4String value = newval;
    if(newval == "") {
      nmsam->saveNeutronSpectrumToFile();
    }
    else {
      nmsam->saveNeutronSpectrumToFile(value);
    }
  }

  if(cmd == eventOffsetCmd) {
    nmsam->SetEventOffset(eventOffsetCmd->ConvertToInt(newval));
  }
  if(cmd == writeResultsToFileCmd) {
    G4String value = newval;
    if(newval == "") {
      nmsam->exportToFile();
    }
    else {
      nmsam->exportToFile(value);
    }
  }
  if(cmd == exportPulsetrainCmd) {
    nmsam->SetPulsetrainExport(exportPulsetrainCmd->ConvertToBool(newval));
  }

  if(cmd == exportResultsAfterRunCmd) {
    nmsam->SetExportAfterRun(exportResultsAfterRunCmd->ConvertToBool(newval));
  }

  if(cmd == loadPulsetrainFromFileCmd) {
    nmsam->loadPulsetrainFromFile(newval);
  }

  if(cmd == loadPulsetrainFromShakesFileCmd) {
    nmsam->loadPulsetrainFromShakesFile(newval);
  }

  if(cmd == addPulsetrainFromFileCmd) {
    std::istringstream is(newval);
    G4String filename;
    G4double offset;
    is >> filename >> offset;
    nmsam->addPulsetrainFromFile(filename, offset);
  }

  if(cmd == splitPulsetrainCmd) {
    std::istringstream is(newval);
    G4String filename;
    G4int splits;
    is >> filename >> splits;
    nmsam->splitPulsetrain(filename, splits);
  }

  if(cmd == savePulsetrainCmd) {
    nmsam->GetPulsetrainManager()->savePulsetrainToFile(newval);
  }

  if(cmd == resetPulsetrainCmd) {
    nmsam->GetPulsetrainManager()->reset();
    nmsam->SetSettingsChanged();
  }

  if(cmd == recalculateCmd) {
    nmsam->SetSettingsChanged();
    nmsam->calculateResults();
  }

  if(cmd == setMultiplicityLengthCmd) {
    nmsam->GetPulsetrainManager()->setMultiplicityLength(setMultiplicityLengthCmd->GetNewIntValue(newval));
    nmsam->SetSettingsChanged();
  }
  
  // Shift register settings Commands
  if(cmd == setRegisterLengthCmd) {
    nmsam->GetPulsetrainManager()->setRegisterLength(setRegisterLengthCmd->GetNewIntValue(newval));
    nmsam->SetSettingsChanged();
  }
  if(cmd == setRegisterPeriodCmd) {
    nmsam->GetPulsetrainManager()->setRegisterPeriod(setRegisterPeriodCmd->GetNewDoubleValue(newval) / microsecond);
    nmsam->SetSettingsChanged();
  }
  if(cmd == setPredelayCmd) {
    nmsam->GetPulsetrainManager()->setPredelay(setPredelayCmd->GetNewDoubleValue(newval) / microsecond);
    nmsam->SetSettingsChanged();
  }
  if(cmd == setLongdelayCmd) {
    nmsam->GetPulsetrainManager()->setAdelay(setLongdelayCmd->GetNewDoubleValue(newval) / microsecond);
    nmsam->SetSettingsChanged();
  }
  if(cmd == setGateCmd) {
    nmsam->GetPulsetrainManager()->setGate(setGateCmd->GetNewDoubleValue(newval) / microsecond);
    nmsam->SetSettingsChanged();
  }

  // Quantization Commands
  if(cmd == setDerandomizePeriodCmd) {
    nmsam->GetPulsetrainManager()->setDerandomizePeriod(setDerandomizePeriodCmd->GetNewDoubleValue(newval) / microsecond);
    nmsam->SetSettingsChanged();
  }
  if(cmd == derandomizeDoCmd) {
    nmsam->GetPulsetrainManager()->setDerandomize(derandomizeDoCmd->GetNewBoolValue(newval));
    nmsam->SetSettingsChanged();
  }
  if(cmd == registerQuantizationDoCmd) {
    nmsam->GetPulsetrainManager()->setRegisterQuantization(registerQuantizationDoCmd->GetNewBoolValue(newval));
    nmsam->SetSettingsChanged();
  }

  // Detector parameters
  if(cmd == setEfficiencyCmd) {
    nmsam->GetPulsetrainManager()->setEfficiency(setEfficiencyCmd->GetNewDoubleValue(newval));
    nmsam->SetSettingsChanged();
  }
  if(cmd == dieawayMethodCmd) {
    G4int setting = dieawayMethodCmd->GetNewIntValue(newval);
    if(setting > 1) { 
      G4cout << "This is not a valid dieaway calculation method." << G4endl;
    }
    else {
      nmsam->GetPulsetrainManager()->setDieawayMethod(DieawayMethod(setting));
      nmsam->SetSettingsChanged();
    }
  }
  if(cmd == dieawayCmd) {
    nmsam->GetPulsetrainManager()->setDieaway(dieawayCmd->GetNewDoubleValue(newval) / microsecond);
    nmsam->SetSettingsChanged();
  }
  if(cmd == infiniteGateCmd) {
    nmsam->GetPulsetrainManager()->setInfiniteGate(infiniteGateCmd->GetNewDoubleValue(newval) / microsecond);
    nmsam->SetSettingsChanged();
  }
  if(cmd == gatefractionMethodCmd) {
    G4int setting = gatefractionMethodCmd->GetNewIntValue(newval);
    if(setting > 3) { 
      G4cout << "This is not a valid gate fraction calculation method." << G4endl;
    }
    else {
      nmsam->GetPulsetrainManager()->setGatefractionMethod(GatefractionMethod(setting));
      nmsam->SetSettingsChanged();
    }
  }
  if(cmd == gateFractionsCmd) {
    std::istringstream is(newval);
    G4double fd;
    G4double ft;
    is >> fd >> ft;
    nmsam->GetPulsetrainManager()->setGateFractions(fd, ft);
    nmsam->SetSettingsChanged();
  }


  // Deadtime Commands
  if(cmd == predeadtimeCmd) {
    nmsam->GetPulsetrainManager()->setPreDeadtime(predeadtimeCmd->GetNewDoubleValue(newval) / microsecond);
    nmsam->SetSettingsChanged();
  }
  if(cmd == predeadtimeUpdatingCmd) {
    nmsam->GetPulsetrainManager()->setPreDeadtimeUpdating(predeadtimeUpdatingCmd->GetNewBoolValue(newval));
    nmsam->SetSettingsChanged();
  }
  if(cmd == postdeadtimeCmd) {
    nmsam->GetPulsetrainManager()->setPostDeadtime(postdeadtimeCmd->GetNewDoubleValue(newval) / microsecond);
    nmsam->SetSettingsChanged();
  }
  if(cmd == postdeadtimeUpdatingCmd) {
    nmsam->GetPulsetrainManager()->setPostDeadtimeUpdating(postdeadtimeUpdatingCmd->GetNewBoolValue(newval));
    nmsam->SetSettingsChanged();
  }

  // Pulsetrain Analysis Method
  if(cmd == pulsetrainAnalysisModeCmd) {
    nmsam->SetPulsetrainAnalysisMode(pulsetrainAnalysisModeCmd->GetNewIntValue(newval));
  }

  if(cmd == setRuntimeFromPulsetrainCmd) {
    nmsam->SetRuntimeFromPulsetrain(setRuntimeFromPulsetrainCmd->GetNewBoolValue(newval));
  }
    
  // RunMode
  if(cmd == setRunModeCmd) {
    if(newval == "He-Detector") {
      nmsam->SetRunMode(RUNMODE_HE_DETECTOR);
    }
    else if(newval == "AlphaN") {
      nmsam->SetRunMode(RUNMODE_ALPHAN);
    }
  }


}
















