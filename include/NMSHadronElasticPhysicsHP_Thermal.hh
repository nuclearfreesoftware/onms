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


// This file is a modification of G4HadronElasticPhysicsHP.hh of Geant4 10.2
//
// DERIVED CLASS

#ifndef NMSHadronElasticPhysicsHP_Thermal_h
#define NMSHadronElasticPhysicsHP_Thermal_h 1

#include "globals.hh"
#include "G4VPhysicsConstructor.hh"

class G4HadronElasticPhysics;

class NMSHadronElasticPhysicsHP_Thermal : public G4VPhysicsConstructor
{
public: 

  NMSHadronElasticPhysicsHP_Thermal(G4int ver = 1); 

  virtual ~NMSHadronElasticPhysicsHP_Thermal();

  // This method will be invoked in the Construct() method. 
  // each particle type will be instantiated
  virtual void ConstructParticle();
 
  // This method will be invoked in the Construct() method.
  // each physics process will be instantiated and
  // registered to the process manager of each particle type 
  virtual void ConstructProcess();

private:

  NMSHadronElasticPhysicsHP_Thermal(NMSHadronElasticPhysicsHP_Thermal &);
  NMSHadronElasticPhysicsHP_Thermal & operator=(const NMSHadronElasticPhysicsHP_Thermal &right);

  G4int    verbose;
  static G4ThreadLocal G4bool   wasActivated;
  static G4ThreadLocal G4HadronElasticPhysics* mainElasticBuilder;
};


#endif








