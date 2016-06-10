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

#ifndef NMSAnalysisManagerMessenger_h
#define NMSAnalysisManagerMessenger_h 1

#include "G4UImessenger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithABool.hh"

class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;

#include "NMSAnalysisManager.hh"
class NMSAnalysisManager;

class NMSAnalysisManagerMessenger : public G4UImessenger {
public:
  NMSAnalysisManagerMessenger(NMSAnalysisManager* );
  ~NMSAnalysisManagerMessenger();

  void SetNewValue(G4UIcommand* cmd, G4String newval);

private:
  NMSAnalysisManager* nmsam;

  G4UIdirectory* analysisDir;
  G4UIdirectory* ptdataDir;

  G4UIcmdWithAString* setDetectorVolumeCmd;
  G4UIcmdWithAString* setLostNeutronVolumeCmd;

  G4UIcmdWithoutParameter* showResultsCmd;

  G4UIcmdWithAString * saveSourceNeutronSpectrumToFileCmd;

  G4UIcmdWithAnInteger * eventOffsetCmd;
  
  G4UIcmdWithAString* writeResultsToFileCmd;
  G4UIcmdWithABool* exportPulsetrainCmd;
  G4UIcmdWithABool* exportResultsAfterRunCmd;
  G4UIcmdWithAString* loadPulsetrainFromFileCmd;
  G4UIcmdWithAString* loadPulsetrainFromShakesFileCmd;
  G4UIcommand* addPulsetrainFromFileCmd;
  G4UIcommand* splitPulsetrainCmd;
  G4UIcmdWithAString* savePulsetrainCmd;
  G4UIcmdWithoutParameter* resetPulsetrainCmd;

  G4UIcmdWithoutParameter* recalculateCmd;
  
  G4UIcmdWithAnInteger* setMultiplicityLengthCmd;

  G4UIcmdWithAnInteger* setRegisterLengthCmd;
  G4UIcmdWithADoubleAndUnit* setRegisterPeriodCmd;
  G4UIcmdWithADoubleAndUnit* setPredelayCmd;
  G4UIcmdWithADoubleAndUnit* setLongdelayCmd;
  G4UIcmdWithADoubleAndUnit* setGateCmd;

  G4UIcmdWithADoubleAndUnit* setDerandomizePeriodCmd;
  G4UIcmdWithABool* derandomizeDoCmd;
  G4UIcmdWithABool* registerQuantizationDoCmd;

  G4UIcmdWithADouble* setEfficiencyCmd;
  G4UIcmdWithAnInteger* dieawayMethodCmd;
  G4UIcmdWithADoubleAndUnit* dieawayCmd;
  G4UIcmdWithADoubleAndUnit* infiniteGateCmd;
  G4UIcmdWithAnInteger* gatefractionMethodCmd;
  G4UIcommand* gateFractionsCmd;
  
  G4UIcmdWithADoubleAndUnit* predeadtimeCmd;
  G4UIcmdWithABool* predeadtimeUpdatingCmd;
  G4UIcmdWithADoubleAndUnit* postdeadtimeCmd;
  G4UIcmdWithABool* postdeadtimeUpdatingCmd;
   
  G4UIcmdWithAnInteger * pulsetrainAnalysisModeCmd;
  G4UIcmdWithABool * setRuntimeFromPulsetrainCmd;
  
  G4UIcmdWithAString* setRunModeCmd;

};

#endif /* NMSAnalysisManagerMessenger_h */
