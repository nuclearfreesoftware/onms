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

#ifndef NMSSFSingleIsoDecaySource_h
#define NMSSFSingleIsoDecaySource_h 1

#include "G4DataVector.hh"
#include "Randomize.hh"

#include "G4IonTable.hh"
#include "G4RadioactiveDecay.hh"
#include "G4VDecayChannel.hh"
#include "G4DecayProducts.hh"
#include "G4AlphaDecayChannel.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4SPSPosDistribution.hh"
#include "G4SPSAngDistribution.hh"

#include "G4SingleParticleSource.hh"

#include "NMSPrimaryUserInformation.hh"
#include "NMSNewSingleDecaySource.hh"

class NMSSFSingleIsoDecaySource : public NMSNewSingleDecaySource 
{
public:
  NMSSFSingleIsoDecaySource();
  ~NMSSFSingleIsoDecaySource();
  static double UniformRand();

  void SetCf252n(G4int, G4int);
  
  void GeneratePrimaryVertex(G4Event* anEvent);

  void Init();

  void setIsotope(G4int iso);
  void setDecayType(G4int type);

  void setIsoSourceType(G4int iso);

  G4bool neutronSource() {
    return true;
  }
  G4double neutronPerEvent();

private:
  G4int DecayIsotope;
  G4int DecayType;

  G4int cf252ndist;
  G4int cf252neng;
  
  G4double list;
};

#endif /* NMSSFSingleIsoDecaySource_h */
