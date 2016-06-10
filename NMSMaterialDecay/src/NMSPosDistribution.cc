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
 * This file incorporates work covered by the following copyright and  
 * permission notice:
 * 
 *
 *     ********************************************************************
 *     * License and Disclaimer                                           *
 *     *                                                                  *
 *     * The  Geant4 software  is  copyright of the Copyright Holders  of *
 *     * the Geant4 Collaboration.  It is provided  under  the terms  and *
 *     * conditions of the Geant4 Software License,  included in the file *
 *     * LICENSE and available at  http://cern.ch/geant4/license .  These *
 *     * include a list of copyright holders.                             *
 *     *                                                                  *
 *     * Neither the authors of this software system, nor their employing *
 *     * institutes,nor the agencies providing financial support for this *
 *     * work  make  any representation or  warranty, express or implied, *
 *     * regarding  this  software system or assume any liability for its *
 *     * use.  Please see the license in the file  LICENSE  and URL above *
 *     * for the full disclaimer and the limitation of liability.         *
 *     *                                                                  *
 *     * This  code  implementation is the result of  the  scientific and *
 *     * technical work of the GEANT4 collaboration.                      *
 *     * By using,  copying,  modifying or  distributing the software (or *
 *     * any work based  on the software)  you  agree  to acknowledge its *
 *     * use  in  resulting  scientific  publications,  and indicate your *
 *     * acceptance of all terms of the Geant4 Software license.          *
 *     ********************************************************************
 *
 */

#include "NMSPosDistribution.hh"

NMSPosDistribution::NMSPosDistribution(){
  Confine = false;
}

NMSPosDistribution::~NMSPosDistribution(){
  posVector.clear();
}

void NMSPosDistribution::SetVerboseLevel(G4int verbose) {
  verboseLevel = verbose;
}


void NMSPosDistribution::ClearAll(){
  posVector.clear();
}

void NMSPosDistribution::AddaPosDist(G4SPSPosDistribution* posDist){
  posVector.push_back(posDist);

  posDist->SetPosDisType(SourcePosType);
  posDist->SetPosDisShape(Shape);
  posDist->SetCentreCoords(CentreCoords);
  posDist->SetRadius(Radius);
  posDist->SetRadius0(Radius0);
  posDist->SetHalfX(halfx);
  posDist->SetHalfY(halfy);
  posDist->SetHalfZ(halfz);

  if(Confine) {
    if(verboseLevel >= 2) {
      G4cout << "NMSPosDistribution: Confine to Volume " << VolName << G4endl;
    }
    posDist->ConfineSourceToVolume(VolName);
  }
    
  // Set all Settings for posDist
}

void NMSPosDistribution::DeleteaPosDist(G4int idx){
  size_t sizeIdx = size_t (idx);
  if(sizeIdx <= posVector.size()){
    posVector.erase(posVector.begin() + idx);
  }
  else {
    G4cout << " source index is invalid " << G4endl;
    G4cout << "    it shall be <= " << posVector.size() << G4endl;
  }
}

void NMSPosDistribution::SetPosDisType(G4String type){
  SourcePosType = type;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->SetPosDisType(type);
  }
}

void NMSPosDistribution::SetPosDisShape(G4String shape){
  Shape = shape;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->SetPosDisShape(shape);
  }
}


void NMSPosDistribution::ConfineSourceToVolume(G4String volume){
  Confine = true;
  VolName = volume;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->ConfineSourceToVolume(volume);
  }
}

void NMSPosDistribution::SetCentreCoords(G4ThreeVector ccoords){
  CentreCoords = ccoords;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->SetCentreCoords(ccoords);
  }

}

// FIX Add PosRot

// FIX implement
void NMSPosDistribution::SetPosRot2(G4ThreeVector posrot){

}

void NMSPosDistribution::SetHalfX(G4double hx){
  halfx = hx;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->SetHalfX(hx);
  }
}

void NMSPosDistribution::SetHalfY(G4double hy){
  halfy = hy;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->SetHalfY(hy);
  }

}

void NMSPosDistribution::SetHalfZ(G4double hz){
  halfz = hz;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->SetHalfZ(hz);
  }
}

void NMSPosDistribution::SetRadius(G4double r){
  Radius = r;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->SetRadius(r);
  }

}

void NMSPosDistribution::SetRadius0(G4double r0){
  Radius0 = r0;
  for(size_t i = 0; i < posVector.size(); i++) {
    posVector[i]->SetRadius0(r0);
  }
}
