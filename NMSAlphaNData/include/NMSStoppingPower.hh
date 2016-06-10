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

#ifndef NMSStoppingPower_h
#define NMSStoppingPower_h 1

#include "G4Material.hh"
#include "G4ASTARStopping.hh"
#include "G4EmCalculator.hh"
#include "G4Alpha.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"

enum StoppingDataSource { STOPPINGPOWER_DS_G4EMCALCULATOR,
			  STOPPINGPOWER_DS_ZIEGLER_T3,
			  STOPPINGPOWER_DS_ASTAR,
			  STOPPINGPOWER_DS_ASTARZIEGLER
};


class NMSStoppingPower
{
public:
  NMSStoppingPower();
  virtual ~NMSStoppingPower();
  void initialize();
  
  G4double getDEDX(G4double E, G4int Z);
  G4double getDEDX(G4double E, G4Material * mat);
  G4double getDEDXfromZieglerT3(G4double E, G4Material* mat);
  G4double getDEDXfromEMCal(G4double E, G4Material* mat);
  G4double getDEDXfromASTAR(G4double E, G4Material* mat);

  G4bool ASTARapplicable(G4double E, G4Material* mat);
  
  void setDataSource(StoppingDataSource newds);
  StoppingDataSource getDataSource() { return ds; }

private:
  G4bool initialized;
  StoppingDataSource ds;
  G4ParticleDefinition * alpha;
  G4ASTARStopping * astar;
  G4EmCalculator emcal;
  G4String emprocessname;

  G4double zieglerdata[92][5];
};



#endif /* NMSStoppingPower_h */
