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

//******************************************************************************
// NMSMaterialDecaySource.hh
//
// Wrapper Class for G4GeneralParticleSource
//******************************************************************************

#ifndef NMSMaterialDecaySource_h
#define NMSMaterialDecaySource_h 1

#include "G4GeneralParticleSource.hh"
#include "G4SPSPosDistribution.hh"
#include "G4Material.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"

#include "G4AlphaDecayChannel.hh"
#include "G4AlphaDecay.hh"
#include "G4BetaMinusDecayChannel.hh"

#include "G4RunManager.hh"

#include "NMSMultipleDecaySource.hh"
#include "NMSPosDistribution.hh"
#include "NMSAlphaNSet.hh"
#include "NMSMaterialDecaySettings.hh"
#include "NMSMaterialDecaySourceMessenger.hh"

#include "NMSANYield.hh"

class NMSMaterialDecaySourceMessenger;

class NMSMaterialDecaySource
{
public:
  NMSMaterialDecaySource();
  NMSMaterialDecaySource(G4Material* sourceMat);
  ~NMSMaterialDecaySource();

  // for common constructor instructions;
  void Init();

  G4Material* GetSourceMaterial() {return nmsmdsSettings->GetSourceMaterial(); };
  void SetSourceMaterial(G4Material* sourceMat);

  std::vector< G4VPhysicalVolume* > PhysicalVolumeRelationshipVector(G4VPhysicalVolume*, G4LogicalVolume*);
  void SetSourceFromPhysicalVolumeName(G4String volname);

  void SetActivity(G4double act);
  void SetActiveVolume(G4double vol);
  G4double GetActivity();
  
  NMSPosDistribution* GetPosDist() { return posGenerator; };

  // Set / Unset particles produced
  G4bool GetSpontaneousFissionNeutron() {return spontaneousFissionNeutron; };
  G4bool GetSpontaneousFissionGamma() {return spontaneousFissionGamma;}
  G4bool GetAlphaDecay() {return alphaDecay; };
  G4bool GetBetaDecay() {return betaDecay; };
  G4bool GetAlphaN() {return alphaN; };

  void SetSpontaneousFission(G4bool neutron, G4bool gamma);
  void SetAlphaDecay(G4bool status);
  void SetBetaDecay(G4bool status);
  void SetAlphaNSource(G4bool status);

  void SetAlphaNFile(G4String filename = "std.alphan");
  void SetEnergyFile(G4String filename = "std.energy");

  void GeneratePrimaryVertex(G4Event*);

  void SetEventTimeLimits(G4double start, G4double end);
  G4double GetStartTime() {return startSourceTimeDistribution; };
  G4double GetEndTime() {return endSourceTimeDistribution; };
  G4double GetRuntime() {return endSourceTimeDistribution - startSourceTimeDistribution; };

  void LoadSources();

  void SetVerboseLevel(G4int i);
  void DumpIsotopeData();
  void DumpSourceStatus();

  void WriteANSpectrum(G4String filename);
  void WriteANYield(G4String filename);

private:
  G4bool spontaneousFissionNeutron;
  G4bool spontaneousFissionGamma;
  G4bool alphaDecay;
  G4bool betaDecay;
  G4bool alphaN;

  G4bool sourcesloaded;

  //G4Material* currentSourceMaterial;

  G4double activity;
  G4double activeVolume;

  G4double materialIntensity; // Intensity per ???

  NMSMultipleDecaySource* sourceGenerator;
  NMSPosDistribution* posGenerator;

  G4String alphaNFilename;
  G4String energyFilename;

  G4double startSourceTimeDistribution;
  G4double endSourceTimeDistribution;

  G4int verboseLevel;

  NMSMaterialDecaySettings * nmsmdsSettings;
  NMSMaterialDecaySourceMessenger * nmsmdsMessenger;

  NMSANYield * yieldCalculator;
  
private:
  // Timing issues
  G4double GetNextEventTime();
  void SetTime(G4double);

  G4double GetSFBranching(G4int a, G4int z);
  G4double GetAlphaBranching(G4int a, G4int z);
  G4double GetBetaBranching(G4int a, G4int z);
};


#endif
