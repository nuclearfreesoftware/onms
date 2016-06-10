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

#ifndef NMSRunManagerMessenger_h
#define NMSRunManagerMessenger_h 1

#include <fstream>

#include "G4UImessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"

#include "G4RNGHelper.hh"

class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;

#include "NMSRunManager.hh"
class NMSRunManager;


class NMSRunManagerMessenger : public G4UImessenger {
public:
  NMSRunManagerMessenger(NMSRunManager*);
  ~NMSRunManagerMessenger();

  void SetNewValue(G4UIcommand* cmd, G4String newval);

private:
  NMSRunManager* rm;

  G4UIdirectory* nmsDir;
  G4UIdirectory* runDir;

  G4UIcmdWithADoubleAndUnit* runtimeCmd;
  G4UIcmdWithAnInteger* printModuloCmd;
  G4UIcmdWithADouble* runFromTimeCmd;
  // alphaN???

  G4UIcommand* randomSeedCmd;
};

#endif /* NMSRunManagerMessenger_h */
