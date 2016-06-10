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

#ifndef NMSDetectedEventSet_h
#define NMSDetectedEventSet_h 1

#include <vector>
#include "NMSDetectedEvent.hh"

class NMSDetectedEventSet : public std::vector<NMSDetectedEvent> {
public:
  NMSDetectedEventSet();
  virtual ~NMSDetectedEventSet();

  // retrieve and store

};

#endif /* NMSDetectedEventSet_h */
