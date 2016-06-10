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


#ifndef NMSPrimaryUserInformation_h
#define NMSPrimaryUserInformation_h 1

#include "G4VUserPrimaryParticleInformation.hh"

enum NMSOrigin { ORIGIN_ALPHA_N, ORIGIN_SF, ORIGIN_ALPHA };

class NMSPrimaryUserInformation : public G4VUserPrimaryParticleInformation {
public:
  NMSPrimaryUserInformation();
  ~NMSPrimaryUserInformation();

  void Print() const;

  NMSOrigin GetOrigin();
  void SetOrigin(NMSOrigin newor);

private:
  NMSOrigin ori;

};

#endif /* NMSPrimaryUserInformation_h */
