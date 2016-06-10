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



#include "NMSNewSingleDecaySource.hh"

NMSNewSingleDecaySource::NMSNewSingleDecaySource() : G4SingleParticleSource() {
  sourcetype = "Undefined";
  verboseLevel = 0;
  initialized = false;
}


NMSNewSingleDecaySource::NMSNewSingleDecaySource(G4String st) : G4SingleParticleSource(){
  sourcetype = st;
  verboseLevel = 0;
  initialized = false;
}

NMSNewSingleDecaySource::~NMSNewSingleDecaySource() {

}

void NMSNewSingleDecaySource::SetVerboseLevel(G4int vl) {
  verboseLevel = vl;
}

G4PrimaryVertex * NMSNewSingleDecaySource::getNewGeneralVertex() {
  G4ThreeVector sourcePosition = GetPosDist()->GenerateOne();

  if(verboseLevel >= 2) {
    G4cout << "New Source Event" << G4endl;
    G4cout << "Time:     " << GetParticleTime() / second << G4endl;
    G4cout << "Position: " << sourcePosition / cm << G4endl;
    G4cout << "Position Distribution Type: " << GetPosDist()->GetPosDisType() << G4endl;
  }

  return new G4PrimaryVertex(sourcePosition, GetParticleTime());

}
