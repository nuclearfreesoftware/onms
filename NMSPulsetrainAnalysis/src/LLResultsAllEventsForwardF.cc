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

#include "NMSMultiplicityResult.hh"
#include "NMSPulsetrainManager.hh"

#include <iomanip>

NMSMultiplicityResult NMSPulsetrainManager::LLResultsAllEventsForwardF() {

  NMSMultiplicityResult mr(mms.multiplicityLength);
  mr.SetSettings(mms);
  
  unsigned long long lastevent = 0;
  unsigned long long pd = llround(mms.predelay * 1000000);
  unsigned long long ad = llround(mms.adelay * 1000000);
  unsigned long long gate = llround(mms.gate * 1000000);

  // Timing of process
  double start;
  double stop;
  start = clock();

  // No Events => return empty result
  if( eventsprocessed.size() == 0) {
    return mr;
  }

  // Verbose output
  int verbositySteps = eventsprocessed.size() / 50;
  if(verbositySteps == 0) {
    verbositySteps = 1;
  }
  if(verboseLevel >= 1) {
    std::cout << "Calculating";
    std::cout.flush();
  }

  // Loop / count variables
  int j = 0;
  int shiftcount = 0;

  // Main loop
  for(int i = 0; i < eventsprocessed.size(); i++) {
  //for(int i = 0; i < 3000; i++) {
    if(i % verbositySteps == 0) {
      if(verboseLevel >= 1) {
	std::cout << ".";
	std::cout.flush();
      }
    }

    // jump over predelay events
    shiftcount = 0;
    j = i + 1;
    while(j < eventsprocessed.size()
	  && eventsprocessed[i].eventtime + pd > eventsprocessed[j].eventtime ) {
      j++;
    }

    // count events in gate
    while(j < eventsprocessed.size()
	  && eventsprocessed[i].eventtime + pd + gate >= eventsprocessed[j].eventtime) {
      j++;
      //      std::cout << eventsprocessed[j].eventtime << " ";
      shiftcount++;
    }
    //    std::cout << std::endl;
    mr.addRA(shiftcount, 1);

    // jump over long delay events
    shiftcount = 0;
    j = i + 1;
    while(j < eventsprocessed.size()
	  && eventsprocessed[i].eventtime + ad > eventsprocessed[j].eventtime) {
      j++;
    }

    // count events in gate
    while(j < eventsprocessed.size()
	  && eventsprocessed[i].eventtime + ad + gate >= eventsprocessed[j].eventtime) {
      j++;
      shiftcount++;
    }
    mr.addA(shiftcount, 1);
  }

  if(verboseLevel >= 1) {
    std::cout << std::endl;
  }

  mr.setLastEvent(eventsprocessed.back().eventtime);

  stop = clock();
  if(verboseLevel >= 1) {
    std::cout << "WALL: Calculation took " << ( stop - start ) / CLOCKS_PER_SEC << " second" << std::endl;
  }

  return mr;
}
