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
#include <map>

NMSMultiplicityResult NMSPulsetrainManager::LLResultsSingleEventForward(bool infinitegate) {

  NMSMultiplicityResult gateResult(mms.multiplicityLength);
  gateResult.SetSettings(mms);
  NMSMultiplicityResult noGateResult(mms.multiplicityLength);
  noGateResult.SetSettings(mms);
  
  unsigned long long lastevent = 0;
  unsigned long long pd = llround(mms.predelay * 1000000);
  unsigned long long gate = llround(mms.gate * 1000000);
  unsigned long long ad = llround(mms.adelay * 1000000);

  // No Events => return empty result
  if( eventsprocessed.size() == 0) {
    return gateResult;
  }

  // Map pulsetrain by event groups
  std::map<int, size_t> eventmap;
  std::vector< std::vector< int > > eventgroups;
  std::vector< int > eventidx;
  for(size_t i = 0; i < eventsprocessed.size(); i++ ) {
    if ( eventmap.find(eventsprocessed[i].eventid) == eventmap.end() ) {
      eventmap[eventsprocessed[i].eventid] = eventgroups.size();
      eventidx.clear();
      eventidx.push_back(i);
      eventgroups.push_back(eventidx);
    } else {
      eventgroups[eventmap[eventsprocessed[i].eventid]].push_back(i);
    }
  }

  // Calculate results
  for(size_t i = 0; i < eventgroups.size(); i++) {
    for(size_t j = 0; j < eventgroups[i].size(); j++ ) {
      size_t k = j + 1;
      size_t jidx = eventgroups[i][j];
      size_t kidx = eventgroups[i][k];

      // Predelay
      int pcount = 0;
      while(k < eventgroups[i].size()
	    && eventsprocessed[kidx].eventtime < eventsprocessed[jidx].eventtime + pd) {
	k++; kidx = eventgroups[i][k];
	pcount++;
      }

      // Gate
      int gcount = 0;
      while(k < eventgroups[i].size()
	    && eventsprocessed[kidx].eventtime < eventsprocessed[jidx].eventtime + pd + gate) {
	k++; kidx = eventgroups[i][k];
	gcount++;
      }
      
      size_t remainingevents = eventgroups[i].size() - k;
      if(gcount < mms.multiplicityLength) {
	gateResult.addRA(gcount);
      }
      else {
	std::cout << "Gate has " << gcount << " events, that is more than the length of the virtual shift register" << std::endl;
      }
      if(remainingevents + gcount < mms.multiplicityLength) {
	noGateResult.addRA(remainingevents + gcount + pcount);
      }
      else {
	std::cout << "'Infinite' Gate has " << remainingevents + gcount + pcount<< " events, that is more than the length of the virtual shift register" << std::endl;
      }
    }
  }

  if(infinitegate) {
    return noGateResult;
  }
  else {
    return gateResult;
  }
}
