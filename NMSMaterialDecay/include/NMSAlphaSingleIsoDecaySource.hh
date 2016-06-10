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

#ifndef NMSAlphaSingleIsoDecaySource_h
#define NMSAlphaSingleIsoDecaySource_h 1

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

class NMSAlphaSingleIsoDecaySource : public NMSNewSingleDecaySource 
{
public:
  NMSAlphaSingleIsoDecaySource();
  ~NMSAlphaSingleIsoDecaySource();

  void GeneratePrimaryVertex(G4Event* anEvent);

  void Init();

  void setIsotope(G4int iso);

private:
  G4int DecayIsotope;

  std::vector<G4double> energyList;
  std::vector<G4double> branchingList;

};

#endif /* NMSAlphaSingleIsoDecaySource_h */
