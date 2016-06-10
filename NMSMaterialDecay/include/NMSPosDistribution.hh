/* Copyright (C) 2014-2016, Moritz Kütt
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

#ifndef NMSPosDistribution_h
#define NMSPosDistribution_h 1

#include "G4SPSPosDistribution.hh"
#include <vector>

class NMSPosDistribution
{
public:
  NMSPosDistribution();
  ~NMSPosDistribution();

  void SetVerboseLevel(G4int i);

  void ClearAll();
  void AddaPosDist(G4SPSPosDistribution* posDist);
  void DeleteaPosDist(G4int);

  void SetPosDisType(G4String); // Point, Plane, Surface, Volume
  inline G4String GetPosDisType() { return SourcePosType; };
  void SetPosDisShape(G4String);
  inline G4String GetPosDisShape() { return Shape; };

  void ConfineSourceToVolume(G4String);
  void SetCentreCoords(G4ThreeVector);
  inline G4ThreeVector GetCentreCoords() { return CentreCoords; } ;
  void SetPosRot2(G4ThreeVector);
  void SetHalfX(G4double);
  inline G4double GetHalfX() { return halfx; } ;
  void SetHalfY(G4double);
  inline G4double GetHalfY()  { return halfy; } ;
  void SetHalfZ(G4double);
  inline G4double GetHalfZ()  { return halfz; } ;
  void SetRadius(G4double);
  inline G4double GetRadius()  { return Radius; };
  void SetRadius0(G4double);
  inline G4double GetRadius0() { return Radius0; };

private:
  G4int verboseLevel;

  std::vector<G4SPSPosDistribution*> posVector;

  G4String SourcePosType; //Point,Plane,Surface,Volume
  G4String Shape; //Circle,Square,Rectangle etc..

  G4double halfx, halfy, halfz; //half lengths
  G4double Radius; //Radius for circles or spheres
  G4double Radius0; // The inner radius of an annulus

  G4ThreeVector CentreCoords; // Coords of centre of input shape

  G4bool Confine; //If true confines source distribution to VolName
  G4String VolName;

};


#endif
