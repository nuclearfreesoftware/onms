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

#ifndef NMSANCSDATA_H
#define NMSANCSDATA_H

#include "G4LPhysicsFreeVector.hh"

struct inelasticdata {
  G4int mt;
  G4double Qex;
  G4LPhysicsFreeVector * cs;
};

class NMSANcsdata {
public:
  static NMSANcsdata* Instance();
  
  NMSANcsdata();
  ~NMSANcsdata();

  void Init();

  void SetVerboseLevel(G4int i) { verboseLevel = i; }
  G4int GetVerboseLevel() { return verboseLevel; }

  G4int no() {return csno; }
  G4bool has(G4int Z, G4int A);
  G4int index(G4int Z, G4int A);
  G4LPhysicsFreeVector * getcsVector(G4int idx);

  G4int inelasticNo(G4int idx);
  G4LPhysicsFreeVector * getInelasticCsVector(G4int idx, G4int ridx);
  G4int getInelasticCsMT(G4int idx, G4int ridx);
  G4double getInelasticCsQex(G4int idx, G4int ridx);
  
public:
  static const G4int csno;
  static const G4int csarray[17][2];
  static const G4String nameString[100];

private:
  G4int verboseLevel;
  std::vector<G4LPhysicsFreeVector *> csIsoVector;
  std::vector< std::vector< inelasticdata > > inelasticCsVector;

};

#endif /* NMSANCSDATA_H */
