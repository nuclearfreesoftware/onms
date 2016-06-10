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

#ifndef NMSRunManager_h
#define NMSRunManager_h 1

#include "G4RunManager.hh"

#include "NMSPrimaryGeneratorAction.hh"
#include "NMSEventAction.hh"

#include "NMSRunManagerMessenger.hh"

class NMSRunManagerMessenger;

class NMSRunManager : public G4RunManager {
public:
  NMSRunManager();
  virtual ~NMSRunManager();

  void SetRuntime(G4double rt);
  inline G4double GetRuntime() { return runtime; };

  void SetPrintModulo(G4int pm);
  inline G4int GetPrintModulo(G4int pm);

  void SetCLOptions(G4bool newgps);

  void RunFromTime(G4double);
private:
  NMSRunManagerMessenger* nmsrmm;
  G4double runtime;
  G4int printModulo;
  G4bool gps;
};


#endif /* NMSRunManager_h */
