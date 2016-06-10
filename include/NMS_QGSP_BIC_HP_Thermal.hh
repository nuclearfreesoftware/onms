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

// File derived from QGSP_BIC_HP from Geant4

#ifndef NMS_QGSP_BIC_HP_Thermal_h
#define NMS_QGSP_BIC_HP_Thermal_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4VModularPhysicsList.hh"
#include "CompileTimeConstraints.hh"

class NMS_QGSP_BIC_HP_Thermal: public G4VModularPhysicsList
{
public:
  NMS_QGSP_BIC_HP_Thermal(G4int ver = 1);
  virtual ~NMS_QGSP_BIC_HP_Thermal();
  
public:
  // SetCuts() 
  virtual void SetCuts();

private:
  //  enum {ok = CompileTimeConstraints::IsA<G4VModularPhysicsList>::ok };
};

#endif



