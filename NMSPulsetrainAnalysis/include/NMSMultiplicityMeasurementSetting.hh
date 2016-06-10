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

#ifndef NMSMultiplicityMeasurementSetting_h
#define NMSMultiplicityMeasurementSetting_h 1

enum GatefractionMethod {
  GM_FROM_DIEAWAY,
  GM_FROM_LONGGATE,
  GM_FROM_EVENTNO,
  GM_FIXED
};

enum DieawayMethod {
  DAM_FIXED,
  DAM_AVERAGELIFETIME
};
  
struct NMSMultiplicityMeasurementSetting {
  double predelay;
  double gate;
  double adelay;

  int registerLength;
  double registerPeriod;
  int quantizeRegisterBuffer;
  bool quantizeRegister;

  double derandomizePeriod;
  int derandomizeBuffer;
  bool derandomize;
  
  double efficiency;
  DieawayMethod dm;
  double dieaway;
  double infinitegate;
  GatefractionMethod gm;
  double doublegatefraction;
  double triplegatefraction;

  double predeadtime;
  bool predeadtimeUpdating;
  double postdeadtime;
  bool postdeadtimeUpdating;

  int multiplicityLength;

  double measurementLength;
};

#endif /* NMSMultiplicityMeasurementSetting_h */
