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
#include "NMSAlphaSingleIsoDecaySource.hh"

#include "G4AlphaDecay.hh"

NMSAlphaSingleIsoDecaySource::NMSAlphaSingleIsoDecaySource() : NMSNewSingleDecaySource("alpha-decay") {
  initialized = false;
  setIsotope(982520);
}

NMSAlphaSingleIsoDecaySource::~NMSAlphaSingleIsoDecaySource() {
  
}

void NMSAlphaSingleIsoDecaySource::GeneratePrimaryVertex(G4Event* anEvent) {

  G4PrimaryVertex* vertex = getNewGeneralVertex();

  G4ParticleDefinition* alpha_definition = G4Alpha::AlphaDefinition();
  G4DynamicParticle* alpha;

  if(!initialized) {
    Init();
  }

  G4double ran = G4UniformRand();
  G4double decayenergy = 0;
  G4int i=0;
  while(ran > branchingList[i])
    {
      i++;
    }
  decayenergy = energyList[i];

  //Calculate Momentum
  alpha = new G4DynamicParticle();
  alpha->SetDefinition(alpha_definition);
  alpha->SetKineticEnergy(decayenergy);
  G4double mom = alpha->GetTotalMomentum();

  //Random, Isotropic direction
  G4double cosTheta = 2*G4UniformRand() - 1., phi = twopi*G4UniformRand();
  G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
  G4double ux = sinTheta*std::cos(phi),
    uy = sinTheta*std::sin(phi),
    uz = cosTheta;

  G4double momx = mom*ux;
  G4double momy = mom*uy;
  G4double momz = mom*uz;


  G4PrimaryParticle* particle = new G4PrimaryParticle(alpha_definition,
                                                      momx, momy, momz,
                                                      decayenergy);

  particle->SetMass(alpha_definition->GetPDGMass());
  particle->SetCharge(alpha_definition->GetPDGCharge());
  // Fix (where should polarization come from)
  // particle->SetPolarization(polar.x(), polar.y(), polar.z());

  if(verboseLevel >= 2) {
    G4cout << " Alpha Decay of " << DecayIsotope << G4endl;
    G4cout << "     Energy:   " << decayenergy << G4endl;
    G4cout << "     Momentum: "<<particle->GetMomentum() << G4endl;
    G4cout << "     Position: "<<vertex->GetPosition() << G4endl;
  }

  vertex->SetPrimary(particle);
  anEvent->AddPrimaryVertex(vertex);

}

void NMSAlphaSingleIsoDecaySource::Init() {
  G4int A;
  G4int Z;
  G4double energy;
  G4ParticleDefinition* ion = 0;
  G4RadioactiveDecay* decay = new G4RadioactiveDecay();
  G4DecayTable* dtable;
  G4VDecayChannel* dc;
  G4DecayProducts* products;

  energyList.clear();
  branchingList.clear();

  Z = DecayIsotope / 10000;
  A = (DecayIsotope - Z * 10000) / 10;
  if (verboseLevel >= 1) {
    G4cout << "================================================================" << G4endl;
  }

  ion =  G4IonTable::GetIonTable()->GetIon(Z,A,0);
  if(verboseLevel >= 1) {
    G4cout << "Loading Alpha Decay for " << ion->GetParticleName() << G4endl;
  }
  dtable = decay->GetDecayTable(ion);
  for(G4int i = 0; i < dtable->entries(); i++) {
    dc = dtable->GetDecayChannel(i);
    // Option A: G4AlphaDecay (was G4AlphaDecayChannel until November 2014)
    // if( G4AlphaDecayChannel* testalpha = dynamic_cast< G4AlphaDecay* >( dc ) ) {
    // Option B: Check for kinematics name
    if(dc->GetKinematicsName() == "alpha decay") {
      products = dc->DecayIt(0);
      for(G4int j = 0; j < products->entries(); j++) {
	if((*products)[j]->GetParticleDefinition()->GetParticleName() == "alpha") {
	  energy = (*products)[j]->GetKineticEnergy();
	  energyList.push_back(energy);
	  branchingList.push_back(dc->GetBR());
	  if(verboseLevel >= 1) {
	    G4cout << " Decay energy: " << energy / MeV << " MeV (Branching Ratio: " << dc->GetBR() << ")" << G4endl;
	  }
	}
	else {
	  //FIX AlphaChannel without alpha?
	}
      }

    }
  }
 

  //normalization + convert to probability
  G4double total = 0.;
  std::vector<G4double> tempList(branchingList.size());
  for (size_t i = 0; i < branchingList.size(); i++) {
    total += branchingList[i];
  }

  branchingList[0] = branchingList[0] / total;
  
  if (verboseLevel >= 1) {
    G4cout << " Branching Ratios are normalized internally: " << G4endl;
  }

  for (size_t i = 1 ;  i < branchingList.size(); i++) {
    if (verboseLevel >= 1) {
      G4cout << " " << branchingList[i] / total << ", summed: ";
    }
    branchingList[i] = branchingList[i] / total + branchingList[i-1];
    if (verboseLevel >= 1) {
      G4cout << " " << branchingList[i] << G4endl;
    }
  }

  if (verboseLevel >= 1) {
    G4cout << "================================================================" << G4endl;
  }

  initialized = true;
  
}

void NMSAlphaSingleIsoDecaySource::setIsotope(G4int iso) {
  if(verboseLevel >= 3) {
    G4cout << "NMSAlphaSingleIsoDecaySource::setIsotope" << G4endl;
  }

  //check CheckIsotope(iso, DecayType);


  //check if we can do isotope + type
  DecayIsotope = iso;
  initialized = false;
}

