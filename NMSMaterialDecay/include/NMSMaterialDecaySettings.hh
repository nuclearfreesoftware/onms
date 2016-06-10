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

#ifndef NMSMaterialDecaySettings_h
#define NMSMaterialDecaySettings_h 1

#include "G4Material.hh"
#include "G4SystemOfUnits.hh"

#include "NMSPosDistribution.hh"

enum NMSANEnergySamplingMode {
  AN_E_FIXED,
  AN_E_PRECALCMC,
  AN_E_CALC,
  AN_E_SPECTRUMFILE
};

enum NMSANPositionSamplingMode {
  AN_POS_SOURCEVOLUME,
  AN_POS_PRECALCMC
};

enum NMSANDirectionSamplingMode {
  AN_DIR_ISOTROPIC,
  AN_DIR_PRECALCMC,
  AN_DIR_CS
};

enum NMSANActivityCalcMode {
  AN_ACT_FIXED,
  AN_ACT_PRECALCMC,
  AN_ACT_CALC
};

class NMSMaterialDecaySettings
{
public:

  static NMSMaterialDecaySettings* Instance();

  G4bool GetSourcesloaded() { return sourcesloaded; }
  void SetSourcesloaded(G4bool sl) { sourcesloaded = sl; }

  G4bool GetActivityFixed() { return activityfixed; }
  void SetActivityFixed(G4bool af);
  
  NMSANDirectionSamplingMode GetANDirectionSamplingMode() { return andsm; }
  void SetANDirectionSamplingMode(NMSANDirectionSamplingMode);
  
  NMSANEnergySamplingMode GetANEnergySamplingMode() { return anesm; }
  void SetANEnergySamplingMode(NMSANEnergySamplingMode);

  NMSANPositionSamplingMode GetANPositionSamplingMode() { return anpsm; }
  void SetANPositionSamplingMode(NMSANPositionSamplingMode);

  G4double GetANEnergy() { return anEnergy; }
  void SetANEnergy(G4double);

  G4String GetANFilename() { return anFilename; }
  void SetANFilename(G4String filename);

  G4String GetANEnergyFilename() { return anEnergyFilename; }
  void SetANEnergyFilename(G4String filename);

  G4bool GetANSpectrumMT91() { return anmt91; };
  void SetANSpectrumMT91(G4bool inc);

  NMSANActivityCalcMode GetANActivityCalcMode() { return anacm; }
  void SetANActivityCalcMode(NMSANActivityCalcMode);

  G4double GetANActivity() { return anActivity; }
  void SetANActivity(G4double);

  void SetCf252n(G4int, G4int);
  G4int GetCf252ndist() { return ndist; }
  G4int GetCf252neng() { return neng; }

  void SetSourceMaterial(G4Material * mat) { currentSourceMaterial = mat; }
  G4Material * GetSourceMaterial() { return currentSourceMaterial; }
  
  G4bool ANwithoutMaterialPossible();

  // Lock Functions
private:
  NMSMaterialDecaySettings();
  virtual ~NMSMaterialDecaySettings();

  void checkAndSetANAllFile();

private:
  // Sources on/off
  G4bool spontaneousFissionNeutron;
  G4bool spontaneousFissionGamma;
  G4bool alphaDecay;
  G4bool betaDecay;
  G4bool alphaN;

  // General sources settings
  G4Material* currentSourceMaterial;

  G4bool activityfixed;
  G4double activity;
  G4double activeVolume;

  G4double materialIntensity; // Intensity per ???

  NMSPosDistribution* posGenerator;

  //Cf252 SF settings
  G4int ndist;
  G4int neng;

  // Alpha N Settings
  NMSANDirectionSamplingMode andsm;
  NMSANEnergySamplingMode anesm;
  NMSANPositionSamplingMode anpsm;
  G4double anEnergy;
  G4String anFilename;
  G4String anEnergyFilename;
  NMSANActivityCalcMode anacm;
  G4double anActivity;
  G4bool anAllFile;
  G4bool anmt91;

  // Timing
  G4double startSourceTimeDistribution;
  G4double endSourceTimeDistribution;

  // Runtime
  G4bool sourcesloaded;
  

  // Debug
  G4int verboseLevel;

};




#endif /* NMSMaterialDecaySettings_h */
