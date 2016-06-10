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

#ifndef NMSMaterialDecaySourceMessenger_h
#define NMSMaterialDecaySourceMessenger_h 1

#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UImessenger.hh"

#include "NMSMaterialDecaySource.hh"
#include "NMSMaterialDecaySettings.hh"

class G4UIcmdWithADoubleAndUnit;

class NMSMaterialDecaySource;

class NMSMaterialDecaySourceMessenger : public G4UImessenger {
public:
  NMSMaterialDecaySourceMessenger(NMSMaterialDecaySource* );
  ~NMSMaterialDecaySourceMessenger();

  void SetNewValue(G4UIcommand* cmd, G4String newval);
  
private:
  NMSMaterialDecaySource* nmsmds;
  NMSMaterialDecaySettings* nmsmdsSettings;

  G4UIdirectory* sourceDir;
  G4UIdirectory* sourcePosDistDir;
  G4UIdirectory* sourceANDir;

  G4UIcmdWithAnInteger * VerboseCmd;
  
  G4UIcmdWithAString* MaterialCmd;
  G4UIcmdWithAString* SourceFromVolumeCmd;
  G4UIcommand* Cf252Cmd;
  G4UIcmdWithADoubleAndUnit* ActiveVolumeCmd;
  G4UIcmdWithADoubleAndUnit* ActivityCmd;
  G4UIcmdWithABool* ActivityFixedCmd;

  G4UIcmdWithABool* SFnCmd;
  G4UIcmdWithABool* SFgCmd;
  G4UIcmdWithABool* BetaCmd;
  G4UIcmdWithABool* AlphaCmd;
  G4UIcmdWithABool* AlphaNCmd;

  G4UIcmdWithAString* ANWriteYieldFileCmd;
  G4UIcmdWithAnInteger* ANSampleDirectionCmd;
  G4UIcmdWithAnInteger* ANSampleEnergyCmd;
  G4UIcmdWithAnInteger* ANSamplePositionCmd;
  G4UIcmdWithADoubleAndUnit* ANEnergyCmd;
  G4UIcmdWithAString* ANFileCmd;
  G4UIcmdWithAString* ANEnergyFileCmd;
  G4UIcmdWithAString* ANWriteEnergyFileCmd;
  G4UIcmdWithABool* ANSpectrumCalcMT91;
  G4UIcmdWithAnInteger* ANActivityCalcCmd;
  G4UIcmdWithADouble* ANActivityCmd;

  G4UIcmdWithAString* posTypeCmd;
  G4UIcmdWithAString* posShapeCmd;
  G4UIcmdWithAString* posConfineVolumeCmd;
  G4UIcmdWithADoubleAndUnit* posRadiusCmd;
  G4UIcmdWithADoubleAndUnit* posRadius0Cmd;
  G4UIcmdWithADoubleAndUnit* posHalfXCmd;
  G4UIcmdWithADoubleAndUnit* posHalfYCmd;
  G4UIcmdWithADoubleAndUnit* posHalfZCmd;
  G4UIcmdWith3VectorAndUnit* posCentreCoordsCmd;

  G4UIcmdWithoutParameter* dumpSourceStatusCmd;
};

#endif /* NMSMaterialDecaySourceMessenger_h */
