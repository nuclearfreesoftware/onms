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

#ifndef NMSANSPECTRUM_H
#define NMSANSPECTRUM_H 1

#include "G4Material.hh"
#include "G4LPhysicsFreeVector.hh"

#include "NMSStoppingPower.hh"
#include "NMSANcsdata.hh"

struct spectrumdata {
  G4double intensity;
  G4double energy;
  G4LPhysicsFreeVector * spectrum;
};

class NMSANSpectrum {
public:
  NMSANSpectrum();
  ~NMSANSpectrum();

  void SetUseMT91(G4bool mt91b = true) { mt91 = mt91b; };
  G4bool GetUseMT91() { return mt91; };

  G4LPhysicsFreeVector * fromSource(G4Material * mat);
  G4LPhysicsFreeVector * fromEnergy(G4double E, G4Material * mat);
  G4LPhysicsFreeVector * fromEnergy(G4double E, G4Material * mat, const G4Isotope * iso);
  G4LPhysicsFreeVector * fromEnergy(G4double E, G4Material * mat, const G4Isotope * iso, G4int lidx);

  void clear();
  
private:
  G4int verboseLevel;
  
  G4int neutronbins;
  G4double neutronEmin;
  G4double neutronEmax;
  G4double ngridbinwidth;

  G4bool mt91;
  
  std::vector<spectrumdata> spectrumvector;
  NMSStoppingPower * stoppingData;
  NMSANcsdata * csData;
};

#endif /* NMSANSPECTRUM_H */

