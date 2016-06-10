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

/** \brief Class for single pulse train events
 *
 */

#ifndef NMSDetectedEvent_h
#define NMSDetectedEvent_h 1

#include <string>

class NMSDetectedEvent {
 public:
  NMSDetectedEvent(double et = 0, double lt = 0, int eid = 0, std::string einfo = "");

  inline double getStarttime() {
    return eventtime - lifetime;
  }

  inline bool operator < (const NMSDetectedEvent& e) const {
    return (eventtime < e.eventtime);
  }

public:
  double eventtime;
  double lifetime;
  int eventid;
  //  std::string info;
};

#endif
