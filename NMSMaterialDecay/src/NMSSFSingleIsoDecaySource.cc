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
#include "NMSSFSingleIsoDecaySource.hh"
#include "fissionEvent.h"
#include "fissionEventData.hh"

NMSSFSingleIsoDecaySource::NMSSFSingleIsoDecaySource() : NMSNewSingleDecaySource("neutron - from SF") {
  setIsotope(982520); // Cf-255
  DecayType = NMSDECAY_SF_N;
  cf252ndist = 0;
  cf252neng = 0;
  
  initialized = false;

  fissionEvent::setRNGd(UniformRand);
}


NMSSFSingleIsoDecaySource::~NMSSFSingleIsoDecaySource() {
  
}

double NMSSFSingleIsoDecaySource::UniformRand() {
  return CLHEP::HepRandom::getTheEngine()->flat();
}

void NMSSFSingleIsoDecaySource::SetCf252n(G4int ndist = 0, G4int neng = 0) {
  cf252ndist = ndist;
  cf252neng = neng;
}

void NMSSFSingleIsoDecaySource::GeneratePrimaryVertex(G4Event* anEvent) {

  G4PrimaryVertex* vertex = getNewGeneralVertex();

  G4ParticleDefinition* neutron_definition = G4Neutron::Neutron();
  G4ParticleDefinition* photon_definition = G4Gamma::Gamma();
  G4double mom, momx, momy, momz, eng;
  G4int nPrompt, gPrompt;

  G4int isotope = DecayIsotope / 10;
  G4double timesf = 0.;

  G4DynamicParticle* it;

  // zero polarisation for particles (as in example of llnl fission)
  G4ThreeVector polar;

  if(verboseLevel >= 1) {
    G4cout << "Spontaneous Fission!" << G4endl;
  }

  fissionEvent * fissionevent = new fissionEvent(isotope, timesf, -1., 0., 0);
  fissionevent->setCf252Option(cf252ndist, cf252neng);

  // this can be carried out without check for decay type
  nPrompt = fissionevent->getNeutronNu();
  gPrompt = fissionevent->getPhotonNu();
  
  if(nPrompt == -1) {
    G4cout << "================================================================" << G4endl;
    G4cout << "Error: NMSSingleDecaySource" << G4endl;
    G4cout << "Isotope is not available in Spontaneous Fission library: " << isotope << G4endl;
    exit(0);
  }

  if(verboseLevel >= 2)
    G4cout << "Creating primaries and assigning to vertex" << G4endl;

  // neutrons
  if(DecayType == NMSDECAY_SF || DecayType == NMSDECAY_SF_N) {
    for(G4int i=0; i<nPrompt; i++)
      {
	it = new G4DynamicParticle();
	it->SetDefinition(neutron_definition);
	//eng = getneng_(&i)*MeV;
	eng = fissionevent->getNeutronEnergy(i);
	it->SetKineticEnergy(eng);
	mom = it->GetTotalMomentum();

	momx = mom * fissionevent->getNeutronDircosu(i);
	momy = mom * fissionevent->getNeutronDircosv(i);
	momz = mom * fissionevent->getNeutronDircosw(i);
	
	G4PrimaryParticle* particle = new G4PrimaryParticle(
							    neutron_definition,
							    momx, momy, momz,
							    eng);
	particle->SetMass(neutron_definition->GetPDGMass());
	particle->SetCharge(neutron_definition->GetPDGCharge());
	particle->SetPolarization(polar.x(), polar.y(), polar.z());

	if(verboseLevel >= 2){
	  G4cout << "Particle name: "<<particle->GetG4code()->GetParticleName() << G4endl;
	  G4cout << "     Momentum: "<<particle->GetMomentum() << G4endl;
	  G4cout << "     Position: "<<vertex->GetPosition() << G4endl;
	}
	vertex->SetPrimary(particle);
      }

  }

  // photons
  if(DecayType == NMSDECAY_SF || DecayType == NMSDECAY_SF_GAMMA) {
    for(G4int i=0; i<gPrompt; i++) {
      it = new G4DynamicParticle();
      it->SetDefinition(photon_definition);
      eng = fissionevent->getPhotonEnergy(i);
      it->SetKineticEnergy(eng);
      mom = it->GetTotalMomentum();

      momx = mom * fissionevent->getPhotonDircosu(i);
      momy = mom * fissionevent->getPhotonDircosv(i);
      momz = mom * fissionevent->getPhotonDircosw(i);

      G4PrimaryParticle* particle = new G4PrimaryParticle(
							  photon_definition,
							  momx, momy, momz,
							  eng);
      particle->SetMass(photon_definition->GetPDGMass());
      particle->SetCharge(photon_definition->GetPDGCharge());
      particle->SetPolarization(polar.x(), polar.y(), polar.z());

      if(verboseLevel >= 2){
	G4cout << "Particle name: "<<particle->GetG4code()->GetParticleName() << G4endl;
	G4cout << "     Momentum: "<<particle->GetMomentum() << G4endl;
	G4cout << "     Position: "<<vertex->GetPosition() << G4endl;
      }
      vertex->SetPrimary(particle);

    }
  }

  anEvent->AddPrimaryVertex( vertex );

  delete fissionevent;
}

void NMSSFSingleIsoDecaySource::Init() {
  // Currently nothing to do...
  
}

void NMSSFSingleIsoDecaySource::setIsoSourceType(G4int iso) {
  std::stringstream ss;
  ss << "neutron - from SF (of " << iso << ")";
  G4String stype = ss.str();
  setSourceType(stype);
}

void NMSSFSingleIsoDecaySource::setIsotope(G4int iso) {
  DecayIsotope = iso;
  setIsoSourceType(iso);
  initialized = false;
}

void NMSSFSingleIsoDecaySource::setDecayType(G4int type) {
  DecayType = type;
}

double NMSSFSingleIsoDecaySource::neutronPerEvent() {
  double nb = 0;
  nb = fissionEventData::getSfNubar(DecayIsotope / 10);
  if(nb != -1) {
    return nb;
  }
  return 0;
}
