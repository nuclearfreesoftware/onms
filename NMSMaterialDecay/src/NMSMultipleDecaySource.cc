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
 *
 *
 * This file incorporates work covered by the following copyright and  
 * permission notice:
 * 
 *     Copyright (c) 2006-2010 Lawrence Livermore National Security, LLC.
 *     Produced at the Lawrence Livermore National Laboratory 
 *     UCRL-CODE-224807.
 *     
 *     All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *     
 *     o   Redistributions of source code must retain the above copyright notice, this list of conditions and the disclaimer below.
 *     
 *     o  Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the disclaimer (as noted below) in the documentation and/or other materials provided with the distribution.
 *     
 *     o  Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *     
 *     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *     Additional BSD Notice
 *     
 *     1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE. 
 *     
 *     2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or process disclosed, or represents that its use would not infringe privately-owned rights. 
 *     
 *     3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
 *
 * 
 */

#include "NMSMultipleDecaySource.hh"

NMSMultipleDecaySource::NMSMultipleDecaySource() {
  verboseLevel = 0;
  posGenerator = new NMSPosDistribution();
  sourceVector.clear();
  sourceIntensity.clear();
  sourceProbability.clear();

  currentSourceIdx = -1;
  currentSource = 0;
  normalised = false;

}

NMSMultipleDecaySource::NMSMultipleDecaySource(NMSNewSingleDecaySource* src, G4double strength) {
  verboseLevel = 0;

  posGenerator = new NMSPosDistribution();
  sourceVector.clear();
  sourceIntensity.clear();
  sourceProbability.clear();

  currentSourceIdx = -1;
  currentSource = 0;

  sourceVector.push_back(src);
  sourceIntensity.push_back(strength);
  currentSource = src;
  currentSourceIdx = G4int(sourceVector.size() - 1);

  normalised = false;
  IntensityNormalization();
}

NMSMultipleDecaySource::~NMSMultipleDecaySource(){
  delete posGenerator;
}

void NMSMultipleDecaySource::IntensityNormalization() {
  G4double total  = 0.;
  size_t i;

  if(verboseLevel >= 2) {
    G4cout << "**** Intensity normalization" << G4endl;
  }


  for (i = 0; i < sourceIntensity.size(); i++) {
    total += sourceIntensity[i] ;
  }
  sourceProbability.clear();
  sourceProbability.push_back(sourceIntensity[0] / total);
  for (i = 1; i < sourceIntensity.size(); i++) {
    sourceProbability.push_back(sourceIntensity[i] / total + sourceProbability[i-1]);
  }

  if(verboseLevel >= 2) {
    DumpSourceTable();
  }


  normalised = true;
}

void NMSMultipleDecaySource::DumpSourceTable(G4double volume) {
  G4cout << "Relative Source Table ==========================================" << G4endl;
  G4cout << " Number of sources: " << sourceIntensity.size() << G4endl;
  G4cout.width(15);
  G4cout << "Cum. Prob.";
  G4cout.width(15);
  G4cout << "Activ./(s*cm3)";
  G4cout.width(15);
  G4cout << "n / (s*cm3)";
  //  G4cout << "Intens./(g*s)";
  G4cout << " Type";
  G4cout << G4endl;
  G4double totneutron = 0;
  for (G4int i = 0; i < sourceProbability.size(); i++) {
    G4cout.width(15);
    G4cout << sourceProbability[i];
    G4cout.width(15);
    G4cout << sourceIntensity[i] * s * cm3;
    G4cout.width(15);
    if(sourceVector[i]->neutronSource()) {
      G4cout << sourceIntensity[i] * s * cm3 * sourceVector[i]->neutronPerEvent();
      totneutron += sourceIntensity[i] * s * cm3 * sourceVector[i]->neutronPerEvent();
    }
    else {
      G4cout << " ";
    }
    //    G4cout << sourceIntensity[i] * s * g;
    G4cout << " " << sourceVector[i]->getSourceType();
    G4cout << G4endl;
  }
  G4cout << "n / (s*cm3) only for neutron emitting decays" << G4endl;
  G4cout << "Total neutrons: " << totneutron << " n / (s * cm3)" << G4endl;
  G4cout << "================================================================" << G4endl;

}

void NMSMultipleDecaySource::GeneratePrimaryVertex(G4Event* anEvent) {
  if(verboseLevel >= 2) {
    G4cout << "================================================================" << G4endl;
    G4cout << "*** NMSMultipleDecaySource GeneratePrimaryVertex" << G4endl;
  }

  size_t i = 0;

  if(!normalised) {
    IntensityNormalization();
  }

  if(verboseLevel >= 2) {
    DumpSourceTable();
  }

  G4double ran = G4UniformRand();

  i=0;
  while(ran > sourceProbability[i]) {
    i++;
  }

  if(verboseLevel >= 2) {
    G4cout << "**** Selected Source: " << sourceVector[i]->getSourceType() << G4endl;
  }

  SetCurrentSourceto(i);

  currentSource->GeneratePrimaryVertex(anEvent);
}

void NMSMultipleDecaySource::AddaSource(NMSNewSingleDecaySource* src, G4double strength){
  currentSource = src;
  sourceVector.push_back(currentSource);
  sourceIntensity.push_back(strength);
  currentSourceIdx = G4int(sourceVector.size() - 1);
  currentSource->SetVerboseLevel(verboseLevel);
  posGenerator->AddaPosDist(currentSource->GetPosDist());

  normalised = false;
  //  IntensityNormalization();
}

void NMSMultipleDecaySource::DeleteaSource(G4int idx) {
  size_t id = size_t (idx);

  if( id <= sourceIntensity.size() ) {
    sourceVector.erase(sourceVector.begin() + idx);
    sourceIntensity.erase(sourceIntensity.begin() + idx);

    if ( currentSourceIdx == idx ) {
      if (sourceVector.size() > 0) {
	currentSource = sourceVector[0];
	currentSourceIdx = 1;
      }
      else {
	currentSource = 0;
	currentSourceIdx = -1;
      }
    }
  }
  else {
    G4cout << " source index is invalid " << G4endl;
    G4cout << "    it shall be <= " << sourceIntensity.size() << G4endl;
  }

  posGenerator->DeleteaPosDist(idx);

  normalised = false;
  //  IntensityNormalization();
}

void NMSMultipleDecaySource::ClearAll() {
  currentSourceIdx = 0;
  currentSource = 0;
  sourceVector.clear();
  sourceIntensity.clear();
  sourceProbability.clear();

  posGenerator->ClearAll();

  normalised = false;
}

void NMSMultipleDecaySource::SetCurrentSourceto(G4int idx) {
  size_t id = size_t (idx);

  if ( id < sourceVector.size() ) {
    currentSourceIdx = idx;
    currentSource = sourceVector[id];
  }
  else {
    G4cout << " source index is invalid " << G4endl;
    G4cout << "    it shall be < " << sourceIntensity.size() << G4endl;
  }
}

void NMSMultipleDecaySource::SetCurrentSourceIntensity(G4double strength) {
  sourceIntensity[currentSourceIdx] = strength;

  normalised = false;
  //  IntensityNormalization();
}

void NMSMultipleDecaySource::SetVerboseLevel(G4int verbosity) {
  verboseLevel = verbosity;

  for(size_t i = 0; i < sourceVector.size(); i++) {
    sourceVector[i]->SetVerboseLevel(verbosity);
  }
}

void NMSMultipleDecaySource::SetParticleTime(G4double time) {
  for(size_t i = 0; i < sourceVector.size(); i++) {
    sourceVector[i]->SetParticleTime(time);
  }
}
