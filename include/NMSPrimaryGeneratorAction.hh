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

#ifndef NMSPrimaryGeneratorAction_h
#define NMSPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

#include "NMSRunManager.hh"

#include "Fission.hh"
#include "NMSMaterialDecaySource.hh"
#include "NMSPosDistribution.hh"

class G4Event;

class NMSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  NMSPrimaryGeneratorAction();
  ~NMSPrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event* anEvent);
  void SetVerboseLevel(G4int level);

  void SetSourceMaterial(G4Material* sourceMat);
  void SetActivity(G4double act);
  G4double GetActivity() { return nmsmdSource->GetActivity(); }
  void SetActiveVolume(G4double vol);

  void SetSpontaneousFission(G4bool neutron, G4bool gamma);
  void SetAlphaDecay(G4bool status);
  void SetBetaDecay(G4bool status);
  void SetAlphaNSource(G4bool status);
  void SetAlphaNFile(G4String filename);

  void SetRuntime(G4double rt);
  G4double GetRuntime() { return nmsmdSource->GetRuntime(); }

  inline NMSPosDistribution* GetPosDist() { return nmsmdSource->GetPosDist(); };
  inline bool GetSpontaneousFissionGamma() { return nmsmdSource->GetSpontaneousFissionGamma(); };
  inline bool GetSpontaneousFissionNeutron() { return nmsmdSource->GetSpontaneousFissionNeutron(); };

private:
  G4int verboseLevel;

  NMSMaterialDecaySource* nmsmdSource;

};

#endif
