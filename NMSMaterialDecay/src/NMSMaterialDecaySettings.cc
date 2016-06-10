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

#include "NMSMaterialDecaySettings.hh"

#include "G4AutoLock.hh"

#include <sys/stat.h>

namespace
{
  G4Mutex singMutex = G4MUTEX_INITIALIZER; //Protects singleton access
}

bool file_exists (const G4String& name); // defined elsewhere

NMSMaterialDecaySettings* NMSMaterialDecaySettings::Instance() {
  G4AutoLock lock(&singMutex);
  static NMSMaterialDecaySettings instance;
  return &instance;
}

void NMSMaterialDecaySettings::SetActivityFixed(G4bool af) {
  activityfixed = af;
  sourcesloaded = false;
}

void NMSMaterialDecaySettings::SetANDirectionSamplingMode(NMSANDirectionSamplingMode ds) {
  andsm = ds;
  checkAndSetANAllFile();
  sourcesloaded = false;
}

void NMSMaterialDecaySettings::SetANEnergySamplingMode(NMSANEnergySamplingMode es) {
  anesm = es;
  checkAndSetANAllFile();
  sourcesloaded = false;
}

void NMSMaterialDecaySettings::SetANPositionSamplingMode(NMSANPositionSamplingMode ps) {
  anpsm = ps;
  checkAndSetANAllFile();
  sourcesloaded = false;
}

void NMSMaterialDecaySettings::checkAndSetANAllFile() {
  if(anesm == AN_E_PRECALCMC && anpsm == AN_POS_PRECALCMC && andsm == AN_DIR_PRECALCMC) {
    anAllFile = true;
  }
  else {
    anAllFile = false;
  }
}

void NMSMaterialDecaySettings::SetANEnergy(G4double energy) {
  if(anesm == AN_E_FIXED) {
    anEnergy = energy;
  }
  else {
    G4cout << "WARNING: Energy could not be set. Energy sampling is not set to fixed energy mode" << G4endl;
  }
  sourcesloaded = false;
}

void NMSMaterialDecaySettings::SetANFilename(G4String filename) {
  // Add Check if file exists!
  if(file_exists(filename)) {
    anFilename = filename;
    sourcesloaded = false;
  }
  else {
    G4cout << "WARNING: Could not find specified file (" << filename << "), hence no file will be used" << G4endl;
  }
}

void NMSMaterialDecaySettings::SetANEnergyFilename(G4String filename) {
  // Add Check if file exists!
  if(file_exists(filename)) {
    anEnergyFilename = filename;
    sourcesloaded = false;
  }
  else {
    G4cout << "WARNING: Could not find specified file (" << filename << "), hence no file will be used" << G4endl;
  }
}

void NMSMaterialDecaySettings::SetANSpectrumMT91(G4bool inc) {
  anmt91 = inc;
  sourcesloaded = false;
}

void NMSMaterialDecaySettings::SetANActivityCalcMode(NMSANActivityCalcMode acm) {
  anacm = acm;
  sourcesloaded = false;
}

void NMSMaterialDecaySettings::SetANActivity(G4double newactivity) {
  if(anacm = AN_ACT_FIXED) {
    anActivity = newactivity;
    sourcesloaded = false;
  }
  else {
    G4cout << "WARNING: (alpha, n) activity could not be set. Calculation method is not set to fixed activity." << G4endl;
  }
}

void NMSMaterialDecaySettings::SetCf252n(G4int nd, G4int ne) {
  ndist = nd;
  neng = ne;
  sourcesloaded = false;
}

G4bool NMSMaterialDecaySettings::ANwithoutMaterialPossible() {
  if(GetANDirectionSamplingMode() == AN_DIR_CS) {
    return false;
  }
  if(GetANEnergySamplingMode() == AN_E_CALC) {
    return false;
  }
  if(GetANActivityCalcMode() == AN_ACT_CALC) {
    return false;
  }
  if(GetANFilename() == "") {
    return false;
  }
  return true;
}

NMSMaterialDecaySettings::

NMSMaterialDecaySettings::NMSMaterialDecaySettings() {
  spontaneousFissionNeutron = true;
  spontaneousFissionGamma = false;
  alphaDecay = false;
  betaDecay = false;
  alphaN = true;

  activityfixed = false;
  
  andsm = AN_DIR_ISOTROPIC;
  anesm = AN_E_FIXED;
  anpsm = AN_POS_SOURCEVOLUME;
  anacm = AN_ACT_FIXED;

  anEnergy = 1 * MeV;
  anActivity = 1 / (s * cm3);

  anmt91 = false;

  sourcesloaded = false;

  ndist = 0;
  neng = 0;
}

NMSMaterialDecaySettings::~NMSMaterialDecaySettings() {
  
}
