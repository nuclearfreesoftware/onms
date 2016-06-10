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

#ifndef NMSANYield_h
#define NMSANYield_h 1

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4LPhysicsFreeVector.hh"
#include "G4SystemOfUnits.hh"

#include "G4IonTable.hh"
#include "G4RadioactiveDecay.hh"
#include "G4AlphaDecayChannel.hh"
#include "G4AlphaDecay.hh"
#include "G4DecayProducts.hh"

#include "NMSANcsdata.hh"

#include "NMSStoppingPower.hh"

struct yielddata {
  G4double activity;
  G4double energy;
  G4double yieldfore;
  G4double yield;
};

class NMSANYield
{
public:
  NMSANYield();
  virtual ~NMSANYield();

  G4double fromSource(G4Material * mat);
  G4double alphaActivityFromSource(G4Material * mat);
  G4double fromEnergy(G4double E, G4Material * mat);

  void initialize(); // store cs data in vector, use method from NMSAlphaLE

  void SetVerboseLevel(G4int i) { verboseLevel = i; }
  G4int GetVerboseLevel() { return verboseLevel; }

public:
  static const G4int csno;
  static const G4int csarray[17][2];
  static const G4String nameString[100];
  
  // Yield Table
  // G4PhysicsVector
private:
  G4int verboseLevel;
  std::vector<G4LPhysicsFreeVector *> csIsoVector;
  std::vector<yielddata> yieldvector;

  NMSANcsdata * csData;
  NMSStoppingPower * stoppingData;
};


#endif /* NMSANYield_h */
