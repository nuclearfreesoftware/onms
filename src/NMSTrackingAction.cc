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

#include "NMSTrackingAction.hh"

NMSTrackingAction::NMSTrackingAction(NMSAnalysisManager* newnmsam) : nmsam(newnmsam) {
  rm = RUNMODE_HE_DETECTOR;
}

NMSTrackingAction::~NMSTrackingAction() {

}

void NMSTrackingAction::SetRunMode(RunMode newrm) {
  rm = newrm;
}


void NMSTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{

  G4int parent = aTrack->GetParentID();
  G4ParticleDefinition * particleType = aTrack->GetDefinition();

  if (parent == 0) { // primary particle
    G4double energy = aTrack->GetKineticEnergy();
    if(particleType==G4Neutron::NeutronDefinition()) {
      nmsam->addSourceNeutron(energy);
    }
    if(particleType==G4Alpha::AlphaDefinition()) {
      nmsam->addSourceAlpha(energy);
    }
    if(particleType==G4Gamma::GammaDefinition()) {
      nmsam->addSourceGamma(energy);
    }
  }
  else {
    if(particleType==G4Neutron::NeutronDefinition()) {
      G4double energy = aTrack->GetKineticEnergy();
      nmsam->addSecondaryNeutron(energy);
      if(rm == RUNMODE_ALPHAN) {
	G4String processname = "Unknown";
	if(aTrack->GetCreatorProcess() != 0) {
	  processname = aTrack->GetCreatorProcess()->GetProcessName();
	}
	if(processname == "alphaInelastic") {
	  //	fRunAction->addProducedNeutron();
	  //	fRunAction->producedEnergyNeutron(energy);
	  G4cout << "Neutron produced by: " << processname << G4endl;
	  G4cout << "Neutron energy: " << energy << G4endl;
	  NMSAlphaNReaction ar;
	  ar.energy = energy / MeV;
	  ar.time = aTrack->GetGlobalTime() / microsecond;
	  ar.position = aTrack->GetPosition() / cm;
	  ar.alphaDirection = aTrack->GetMomentum();   
	  nmsam->addAlphaNReaction(ar);
	  fpTrackingManager->EventAborted();
	}
      }
    }
    else {
      G4Track* tr = (G4Track*) aTrack;
      tr->SetTrackStatus(fStopAndKill);
    }
  }
}

void NMSTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{

}
