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

#include "NMSPrimaryGeneratorAction.hh"

NMSPrimaryGeneratorAction::NMSPrimaryGeneratorAction()
{
  nmsmdSource = new NMSMaterialDecaySource();
  //  nmspgm = new NMSPrimaryGeneratorMessenger(this);
  verboseLevel = 0;
}

NMSPrimaryGeneratorAction::~NMSPrimaryGeneratorAction()
{
  delete nmsmdSource;
}

void NMSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  if(verboseLevel >= 2) {
    G4cout << "** NMS PGA: Generate Primaries" << G4endl;
    G4cout << G4endl;
    G4cout << "Source information" << G4endl;
    G4cout << "Material / Isotope:" << G4endl;
    if(nmsmdSource->GetAlphaDecay()) {
      G4cout << "AlphaDecay: On" << G4endl;
    }
    else {
      G4cout << "AlphaDecay: Off" << G4endl;
    }
    if(nmsmdSource->GetBetaDecay()) {
      G4cout << "Beta decay: On" << G4endl;
    }
    else {
      G4cout << "Beta decay: Off" << G4endl;
    }
    if(nmsmdSource->GetSpontaneousFissionNeutron() || nmsmdSource->GetSpontaneousFissionGamma()) {
      if(nmsmdSource->GetSpontaneousFissionGamma() && nmsmdSource->GetSpontaneousFissionNeutron()) {
	G4cout << "Spontaneous Fission decay: On (n,g)" << G4endl;
      }
      else if(nmsmdSource->GetSpontaneousFissionGamma() && !nmsmdSource->GetSpontaneousFissionNeutron()) {
	G4cout << "Spontaneous Fission decay: On (g)" << G4endl;
      }
      else if(!nmsmdSource->GetSpontaneousFissionGamma() && nmsmdSource->GetSpontaneousFissionNeutron()) {
	G4cout << "Spontaneous Fission decay: On (n)" << G4endl;
      }
    }
    else {
      G4cout << "Spontaneous Fission decay: Off" << G4endl;
    }
  }
  nmsmdSource->GeneratePrimaryVertex(anEvent);
}


void NMSPrimaryGeneratorAction::SetVerboseLevel(G4int level) {
  verboseLevel = level;

  // Set Verbose Level for Sources...

  if(verboseLevel >= 2) {
    G4cout << "** NMS PGA: Set verbosity to " << verboseLevel << "." << G4endl;
  }

  nmsmdSource->SetVerboseLevel(level);
}

void NMSPrimaryGeneratorAction::SetSourceMaterial(G4Material* sourceMat) {
  nmsmdSource->SetSourceMaterial(sourceMat);
}

void NMSPrimaryGeneratorAction::SetActivity(G4double act) {
  nmsmdSource->SetActivity(act);
}

void NMSPrimaryGeneratorAction::SetActiveVolume(G4double vol) {
  nmsmdSource->SetActiveVolume(vol);
}

void NMSPrimaryGeneratorAction::SetSpontaneousFission(G4bool neutron, G4bool gamma) {
  nmsmdSource->SetSpontaneousFission(neutron, gamma);
}

void NMSPrimaryGeneratorAction::SetAlphaDecay(G4bool status) {
  nmsmdSource->SetAlphaDecay(status);
}

void NMSPrimaryGeneratorAction::SetBetaDecay(G4bool status) {
  nmsmdSource->SetBetaDecay(status);
} 

void NMSPrimaryGeneratorAction::SetAlphaNSource(G4bool status) {
  nmsmdSource->SetAlphaNSource(status);
}

void NMSPrimaryGeneratorAction::SetAlphaNFile(G4String filename) {
  nmsmdSource->SetAlphaNFile(filename);
}

void NMSPrimaryGeneratorAction::SetRuntime(G4double rt) {
  nmsmdSource->SetEventTimeLimits(0 * s, rt);
}
