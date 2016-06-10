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

#ifndef NMSMultipleDecaySource_h
#define NMSMultipleDecaySource_h 1

#include "G4Event.hh"
#include "Randomize.hh"

#include "NMSPosDistribution.hh"
#include "NMSNewSingleDecaySource.hh"
#include "NMSAlphaSingleIsoDecaySource.hh"
#include "NMSSFSingleIsoDecaySource.hh"
#include "NMSANSource.hh"

class NMSMultipleDecaySource : public G4VPrimaryGenerator
{
public:
  NMSMultipleDecaySource();
  NMSMultipleDecaySource(NMSNewSingleDecaySource* src, G4double strength);

  ~NMSMultipleDecaySource();

  G4int GetNumberofSource() { return G4int(sourceVector.size()); };

  void IntensityNormalization();
  void GeneratePrimaryVertex(G4Event* anEvent);

  void AddaSource (NMSNewSingleDecaySource* src, G4double strength);

  void DeleteaSource(G4int);
  void ClearAll();

  void SetCurrentSourceto(G4int) ;
  void SetCurrentSourceIntensity(G4double);

  NMSNewSingleDecaySource* GetCurrentSource() {return currentSource;};
  G4int GetCurrentSourceIndex() { return currentSourceIdx; };
  G4double GetCurrentSourceIntensity() { return sourceIntensity[currentSourceIdx]; };

  void SetParticleTime(G4double time);

  // Wrapper Class for all Position Distributions - allows for easier adjustment
  inline NMSPosDistribution* GetAllPosDist() {return posGenerator; };

  void SetVerboseLevel(G4int i);
  void DumpSourceTable(G4double volume = 0);
  
private:
  NMSNewSingleDecaySource* currentSource;
  G4int currentSourceIdx;

  std::vector <NMSNewSingleDecaySource*> sourceVector;
  std::vector <G4double> sourceIntensity;
  std::vector <G4double> sourceProbability;

  G4bool normalised;

  NMSPosDistribution* posGenerator;

  G4int verboseLevel;
};
#endif
