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

#include "NMSPulsetrainManager.hh"
#include <algorithm>
#include <vector>
#include <queue>
#include <map>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

NMSPulsetrainManager::NMSPulsetrainManager(int registerLengthInit)  {
  // Default Times (these times are in microseconds)
  mms.predelay = 4.5;
  mms.gate = 64;
  mms.adelay = 4096;

  // Multiplicity distribution arrays are per default as long as imaginary shift register
  mms.multiplicityLength = registerLengthInit;
  
  // Shiftregister parameters - period is calculated from gate length
  mms.registerLength = registerLengthInit;
  mms.registerPeriod = mms.gate / registerLengthInit;
  
  mms.quantizeRegisterBuffer = 100;
  mms.quantizeRegister = false;

  mms.derandomizePeriod = 0.1; // 100 ns, defined in microseconds
  mms.derandomizeBuffer = 100;
  mms.derandomize = false;

  mms.predeadtime = 0; // defined in microseconds
  mms.predeadtimeUpdating = false; // default: non-updating deadtime
  mms.postdeadtime = 0; // defined in microseconds
  mms.postdeadtimeUpdating = false; // default: non-updating deadtime

  // default: infintegate is twice as long as normal gate
  mms.infinitegate = mms.gate * 2;

  // default: fixed dieaway of 50 microseconds
  mms.dm = DAM_FIXED;
  mms.dieaway = 50;

  // Set default efficiency to 1
  mms.efficiency = 1;

  // gatefractions per default from dieaway (single exponential function)
  mms.gm = GM_FROM_DIEAWAY;
  mms.doublegatefraction = simplegatefraction();
  mms.triplegatefraction = simplegatefraction() * simplegatefraction();

  processed = false;

  verboseLevel = 0;
}

NMSPulsetrainManager::~NMSPulsetrainManager() {
}

void NMSPulsetrainManager::setVerboseLevel(int vlevel) {
  verboseLevel = vlevel;
}

void NMSPulsetrainManager::sortOriginal() {
  double start = 0;
  double stop = 0;
  NMSSimpleEvent te;
  //if(required) {
    start = clock();
    if(verboseLevel >= 1) {
      std::cout << "Sorting data" << std::endl;
    }
    //sort(processedevents.begin(), processedevents.end());
    //sort(eventsprocessed.begin(), eventsprocessed.end());
    sort(events.begin(), events.end());
    stop = clock();
    if (verboseLevel >= 1) {
      std::cout << "WALL: Sort took " << ( stop - start ) / CLOCKS_PER_SEC << " second" << std::endl;
    }
  //}

  if(verboseLevel >= 1) {
    std::cout << "Check if sorting is required." << std::endl;
  }
  start = clock();
  bool required = false;
  //processedevents.clear();
  eventsprocessed.clear();
  if(events.size() > 0) {
    //processedevents.push_back(events[0].eventtime);
    te.eventtime = llround(events[0].eventtime * 1000000);
    te.lifetime = llround(events[0].lifetime * 1000000);
    te.eventid = events[0].eventid;
    eventsprocessed.push_back(te);
  }
  for(int i = 1; i < events.size(); i++) {
    //processedevents.push_back(events[i].eventtime);
    te.eventtime = llround(events[i].eventtime * 1000000);
    te.lifetime = llround(events[i].lifetime * 1000000);
    te.eventid = events[i].eventid;
    eventsprocessed.push_back(te);
    if(events[i].eventtime < events[i - 1].eventtime) {
      required = true;
    }
  }
  stop = clock();
  if (verboseLevel >= 1) {
    std::cout << "WALL: Vector copy / check took " << ( stop - start ) / CLOCKS_PER_SEC << " second" << std::endl;
  }

}

void NMSPulsetrainManager::reset() {
  eventsprocessed.clear();
  eventslostpredeadtime.clear();
  eventslostpostdeadtime.clear();
  eventslostderandomize.clear();
  eventslostregisterquantization.clear();
  events.clear();

  processed = false;
}

void NMSPulsetrainManager::setMultiplicityLength(int ml) {
  mms.multiplicityLength = ml;
}

double NMSPulsetrainManager::simplegatefraction() {
  return exp(- mms.predelay / mms.dieaway) * (1 - exp(-mms.gate / mms.dieaway));
}


void NMSPulsetrainManager::setEfficiency(double eff) {
  mms.efficiency = eff;
}

void NMSPulsetrainManager::setMeasurementLength(double newtime) {
  mms.measurementLength = newtime;
}

void NMSPulsetrainManager::setMeasurementLengthFromPulsetrain() {
  double max = 0;
  for(int i = 0; i < events.size(); i++) {
    if(events[i].eventtime > max) {
      max = events[i].eventtime;
    }
  }
  std::cout << "The last event happened at an eventtime of " << max << " microseconds (" << max / 1000000 << " seconds)." <<  std::endl;
  mms.measurementLength = max;
}

void NMSPulsetrainManager::setPredelay(double t) {
  mms.predelay = t;
  processed = false;
}

void NMSPulsetrainManager::setGate(double t) {
  mms.gate = t;
  mms.registerPeriod = mms.gate / mms.registerLength;
  processed = false;
}

void NMSPulsetrainManager::setAdelay(double t) {
  mms.adelay = t;
  processed = false;
}

void NMSPulsetrainManager::setRegisterLength(int length) {
  mms.registerLength = length;
  processed = false;
}

// in micro second
void NMSPulsetrainManager::setRegisterPeriod(double t) {
  mms.registerPeriod = t;
  mms.gate = mms.registerPeriod * mms.registerLength;
  processed = false;
}

// in MHz
void NMSPulsetrainManager::setFrequency(double frequency) {
  if(frequency != 0) {
    mms.registerPeriod = 1 / frequency;
    mms.gate = mms.registerPeriod * mms.registerLength;
    processed = false;
  }
  else {
    std::cout << "Error: Can not set register frequency to 0." << std::endl;
  }
}

void NMSPulsetrainManager::setRegisterQuantization(bool rq) {
  mms.quantizeRegister = rq;
  processed = false;
}

void NMSPulsetrainManager::setRegisterQuantizationBuffer(int buffer) {
  mms.quantizeRegisterBuffer = buffer;
  processed = false;
}

void NMSPulsetrainManager::setDieawayMethod(DieawayMethod dm) {
  mms.dm = dm;
}

void NMSPulsetrainManager::setDieaway(double da) {
  mms.dieaway = da;
}

void NMSPulsetrainManager::setInfiniteGate(double t) {
  mms.infinitegate = t;
}

void NMSPulsetrainManager:: setGatefractionMethod(GatefractionMethod gm) {
  mms.gm = gm;
}

void NMSPulsetrainManager::setGateFractions(double fd, double ft) {
  mms.doublegatefraction = fd;
  mms.triplegatefraction = ft;
}


void NMSPulsetrainManager::setDerandomizePeriod(double t) {
  mms.derandomizePeriod = t;
  processed = false;
}

void NMSPulsetrainManager::setDerandomizeBuffer(int buffer) {
  mms.derandomizeBuffer = buffer;
  processed = false;
}

void NMSPulsetrainManager::setDerandomize(bool dr) {
  mms.derandomize = dr;
  processed = false;
}

void NMSPulsetrainManager::setPreDeadtime(double dt) {
  mms.predeadtime = dt;
  processed = false;
}

void NMSPulsetrainManager::setPreDeadtimeUpdating(bool dtu) {
  mms.predeadtimeUpdating = dtu;
  processed = false;
}

void NMSPulsetrainManager::setPostDeadtime(double dt) {
  mms.postdeadtime = dt;
  processed = false;
}

void NMSPulsetrainManager::setPostDeadtimeUpdating(bool dtu) {
  mms.postdeadtimeUpdating = dtu;
  processed = false;
}

void NMSPulsetrainManager::deadtime(unsigned long long dt, bool dtu, std::vector<NMSSimpleEvent>& lost) {
  std::vector<NMSSimpleEvent> temp;
  temp = eventsprocessed;
  eventsprocessed.clear();

  unsigned long long deadtimewindow = 0;
  for(int i = 0; i < temp.size(); i++) {
    if(temp[i].eventtime > deadtimewindow) {
      eventsprocessed.push_back(temp[i]);
      deadtimewindow = temp[i].eventtime + dt;
    }
    else {
      lost.push_back(temp[i]);
      if(dtu) {
	deadtimewindow = temp[i].eventtime + dt;
      }
    }
  }
}

  
void NMSPulsetrainManager::quantize(unsigned long long dp, int buffer, std::vector<NMSSimpleEvent>& quantizationLost) {
  std::queue<NMSSimpleEvent> quantizationQueue;
  std::vector<NMSSimpleEvent> temp;
  NMSSimpleEvent te;
  quantizationLost.clear();

  unsigned long long nextPossibleEvent = 0; 

  // Store last possible slot
  unsigned long long lastslot;
  unsigned long long leftover = eventsprocessed.back().eventtime % dp;
  if (leftover == 0) {
    lastslot = eventsprocessed.back().eventtime;
  }
  else {
    lastslot = eventsprocessed.back().eventtime + (dp - leftover);
  }

  temp = eventsprocessed;
  eventsprocessed.clear();
  
  while(!quantizationQueue.empty()) { quantizationQueue.pop(); }
  for(int i = 0; i < temp.size(); i++) {
    // put events on queue into vector, if free quantization slots before current event
    while(!quantizationQueue.empty() && nextPossibleEvent < temp[i].eventtime) {
      if(quantizationQueue.front().eventtime <= nextPossibleEvent) {
	te.eventtime = nextPossibleEvent;
	te.lifetime = quantizationQueue.front().lifetime;
	te.eventid = quantizationQueue.front().eventid;
	eventsprocessed.push_back(te);
	quantizationQueue.pop();
	nextPossibleEvent + dp;
      }
      else {
	unsigned long long leftover = quantizationQueue.front().eventtime % dp;
	if (leftover == 0) {
	  nextPossibleEvent = quantizationQueue.front().eventtime;
	}
	else {
	  nextPossibleEvent = quantizationQueue.front().eventtime + (dp - leftover);
	}
      }
    }
    if(quantizationQueue.size() == buffer) {
      quantizationLost.push_back(quantizationQueue.front());
      quantizationQueue.pop();
    }
    quantizationQueue.push(temp[i]);
  }

  if(verboseLevel >= 2) {
    std::cout << "Cleaning up " << quantizationQueue.size() << " remaining elements" << std::endl;
    std::cout << "Last possible slot " << lastslot << std::endl;
  }
  
  // empty queue if fitting
  while(!quantizationQueue.empty() && nextPossibleEvent <= lastslot) {
    if(quantizationQueue.front().eventtime <= nextPossibleEvent) {
      te.eventtime = nextPossibleEvent;
      te.lifetime = quantizationQueue.front().lifetime;
      te.eventid = quantizationQueue.front().eventid;
      eventsprocessed.push_back(te);
      quantizationQueue.pop();
      nextPossibleEvent += dp;
    }
    else {
      unsigned long long leftover = quantizationQueue.front().eventtime % dp;
      if (leftover == 0) {
	nextPossibleEvent = quantizationQueue.front().eventtime;
      }
      else {
	nextPossibleEvent = quantizationQueue.front().eventtime + (dp - leftover);
      }
    }
  }

  // empty remaining in lost
  while(!quantizationQueue.empty()) {
    quantizationLost.push_back(quantizationQueue.front());
    quantizationQueue.pop();
  }
}

void NMSPulsetrainManager::processEvents() {
  sortOriginal();
  if(mms.predeadtime > 0) {
    deadtime(llround(mms.predeadtime * 1000000), mms.predeadtimeUpdating, eventslostpostdeadtime);
  }
  if(mms.derandomize) {
    quantize(llround(mms.derandomizePeriod * 1000000), mms.derandomizeBuffer, eventslostderandomize);
  }
  if(mms.quantizeRegister) {
    quantize(llround(mms.registerPeriod * 1000000), mms.quantizeRegisterBuffer, eventslostregisterquantization);
  }
  if(mms.predeadtime > 0) {
    deadtime(llround(mms.predeadtime * 1000000), mms.predeadtimeUpdating, eventslostpostdeadtime);
  }
  processed = true;
}

NMSMultiplicityResult NMSPulsetrainManager::getLongGateResults(int mode) {
  NMSMultiplicityResult resi;
  double tmpgate = mms.gate;
  mms.gate = mms.infinitegate;
  if (mode == RESULTS_LL_ALL_EVENTS_NO_SHIFT_BACKWARD) {
    resi = LLResultsAllEventsBackward();
  }
  else if (mode == RESULTS_LL_SINGLE_EVENT_NO_SHIFT_BACKWARD) {
    resi = LLResultsSingleEventBackward();
  }
  else if (mode == RESULTS_LL_SINGLE_EVENT_NO_SHIFT_FORWARD) {
    resi = LLResultsSingleEventForward();
  }
  else if (mode == RESULTS_LL_ALL_EVENTS_NO_SHIFT_FORWARD_FULL) {
    resi = LLResultsAllEventsForwardF();
  }
  else { // Default
    resi = LLResultsAllEventsForward();
  }
  mms.gate = tmpgate;
  
  return resi;
}

NMSGateFraction NMSPulsetrainManager::gatefractionEventNo() {
  unsigned long long pd = llround(mms.predelay * 1000000);
  unsigned long long ad = llround(mms.adelay * 1000000);
  unsigned long long gate = llround(mms.gate * 1000000);

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
  // for(size_t i = 0; i < eventgroups.size(); i++) {
  //   std::cout << "Group " << i;
  //   for(size_t j = 0; j < eventgroups[i].size(); j++ ) {
  //     std::cout << " " << eventsprocessed[eventgroups[i][j]].eventtime;
  //   }
  //   std::cout << std::endl;
  // }
  NMSMultiplicityResult gateResult(mms.registerLength);
  NMSMultiplicityResult noGateResult(mms.registerLength);
  for(size_t i = 0; i < eventgroups.size(); i++) {
    for(size_t j = 0; j < eventgroups[i].size(); j++ ) {
      size_t k = j + 1;
      size_t jidx = eventgroups[i][j];
      size_t kidx = eventgroups[i][k];
      // int pcount = 0;
      // while(k < eventgroups[i].size()
      // 	    && eventsprocessed[kidx].eventtime < eventsprocessed[jidx].eventtime + pd) {
      // 	k++; kidx = eventgroups[i][k];
      // 	pcount++;
      // }
      int gcount = 0;
      while(k < eventgroups[i].size()
	    && eventsprocessed[kidx].eventtime < eventsprocessed[jidx].eventtime + gate) {
	k++; kidx = eventgroups[i][k];
	gcount++;
      }
      size_t remainingevents = eventgroups[i].size() - k;
      if(gcount < mms.registerLength) {
	gateResult.addRA(gcount);
      }
      else {
	std::cout << "Gate has " << gcount << " events, that is more than the length of the virtual shift register" << std::endl;
      }
      if(remainingevents + gcount < mms.registerLength) {
	noGateResult.addRA(remainingevents + gcount);
      }
      else {
	std::cout << "'Infinite' Gate has " << remainingevents + gcount << " events, that is more than the length of the virtual shift register" << std::endl;
      }
      //std::cout << " GateCount=" << gcount << " InfiniteGateCount=" << remainingevents + gcount << std::endl;
    }
  }

  gateResult.normalize();
  noGateResult.normalize();
  // std::cout << gateResult << std::endl;
  // std::cout << noGateResult << std::endl;
  // std::cout << "singles=" << gateResult.singles(1) << " singles infinite=" << noGateResult.singles(1) << std::endl;
  // std::cout << "doubles=" << gateResult.doubles(1) << " doubles infinite=" << noGateResult.doubles(1) << std::endl;
  // std::cout << "triples=" << gateResult.triples(1) << " triples infinite=" << noGateResult.triples(1) << std::endl;
  NMSGateFraction ret;
  ret.doubles = gateResult.doubles(1) / noGateResult.doubles(1);
  ret.triples = gateResult.triples(1) / noGateResult.triples(1);
  return ret;
}

NMSMultiplicityResult NMSPulsetrainManager::getResults(int mode) {
  if(verboseLevel >= 1) {
    std::cout << "Calculate Results" << std::endl;
  }
  if(!processed) {
    processEvents();
  }
  
  switch(mms.dm) {
  case DAM_FIXED:
    if(verboseLevel >= 1) {
      std::cout << "Use fixed dieaway time." << std::endl;
    }
    if(mms.dieaway == 0) {
      mms.dieaway = 50;
    }
    break;
  case DAM_AVERAGELIFETIME:
    if(verboseLevel >= 1) {
      std::cout << "Dieaway calculated as average neutron life time in detector." << std::endl;
    }
    mms.dieaway = getAverageLifetime();
    break;
  }

  NMSMultiplicityResult mr;
  if(mode == RESULTS_ALL_TIME_STEPS_SHIFT) {
    mr = ResultsAllTimeSteps();
  }
  else if (mode == RESULTS_FAST_FORWARD_SHIFT) {
    mr = ResultsFastForward();
  }
  else if (mode == RESULTS_ALL_EVENTS_NO_SHIFT) {
    mr = ResultsAllEvents();
  }
  else if (mode == RESULTS_ALL_EVENTS_NO_SHIFT_FORWARD) {
    mr = ResultsAllEventsForward();
  }
  else if (mode == RESULTS_LL_ALL_EVENTS_NO_SHIFT_FORWARD_FULL) {
    mr = LLResultsAllEventsForwardF();
  }
  else if (mode == RESULTS_LL_ALL_EVENTS_NO_SHIFT_BACKWARD) {
    if(verboseLevel >= 1) {
        std::cout << "Backward Mode, event independent" << std::endl;
    }
    mr = LLResultsAllEventsBackward();
  }
  else if (mode == RESULTS_LL_ALL_EVENTS_NO_SHIFT_FORWARD) {
    if(verboseLevel >= 1) {
        std::cout << "Forward Mode, event independent" << std::endl;
    }
    mr = LLResultsAllEventsForward();
  }
  else if (mode == RESULTS_LL_SINGLE_EVENT_NO_SHIFT_BACKWARD) {
    if(verboseLevel >= 1) {
        std::cout << "Backward Mode, event based" << std::endl;
    }
    mr = LLResultsSingleEventBackward();
  }
  else if (mode == RESULTS_LL_SINGLE_EVENT_NO_SHIFT_FORWARD) {
    if(verboseLevel >= 1) {
        std::cout << "Forward Mode, event based" << std::endl;
    }
    mr = LLResultsSingleEventForward();
  }

  else { // Default
    mr = LLResultsAllEventsForward();
  }

  NMSGateFraction gf;
  switch(mms.gm) {
  case GM_FROM_DIEAWAY:
    mms.doublegatefraction = simplegatefraction();
    mms.triplegatefraction = simplegatefraction() * simplegatefraction();
    break;
  case GM_FROM_LONGGATE: {
    NMSMultiplicityResult lgmr = getLongGateResults(mode);
    mms.doublegatefraction = mr.doubles(1) / lgmr.doubles(1);
    mms.triplegatefraction = mr.triples(1) / lgmr.triples(1);
    break;
  }
  case GM_FROM_EVENTNO: 
    if(mode == RESULTS_LL_ALL_EVENTS_NO_SHIFT_BACKWARD) {
      NMSMultiplicityResult eventresult = LLResultsSingleEventBackward();
      NMSMultiplicityResult infiniteresult = LLResultsSingleEventBackward(true);
      mms.doublegatefraction = eventresult.doubles(1) / infiniteresult.doubles(1);
      mms.triplegatefraction = eventresult.triples(1) / infiniteresult.triples(1);
    }
    else if (mode == RESULTS_LL_ALL_EVENTS_NO_SHIFT_FORWARD
	     or mode == RESULTS_LL_ALL_EVENTS_NO_SHIFT_FORWARD_FULL) {
      NMSMultiplicityResult eventresult = LLResultsSingleEventForward();
      NMSMultiplicityResult infiniteresult = LLResultsSingleEventForward(true);
      mms.doublegatefraction = eventresult.doubles(1) / infiniteresult.doubles(1);
      mms.triplegatefraction = eventresult.triples(1) / infiniteresult.triples(1);
    }
    else if (mode == RESULTS_LL_SINGLE_EVENT_NO_SHIFT_BACKWARD) {
      NMSMultiplicityResult infiniteresult = LLResultsSingleEventBackward(true);
      mms.doublegatefraction = mr.doubles(1) / infiniteresult.doubles(1);
      mms.triplegatefraction = mr.triples(1) / infiniteresult.triples(1);
    }
    else if (mode == RESULTS_LL_SINGLE_EVENT_NO_SHIFT_FORWARD) {
      NMSMultiplicityResult infiniteresult = LLResultsSingleEventForward(true);
      mms.doublegatefraction = mr.doubles(1) / infiniteresult.doubles(1);
      mms.triplegatefraction = mr.triples(1) / infiniteresult.triples(1);
    }
    break;
  case GM_FIXED:
    if(mms.doublegatefraction == 0 || mms.triplegatefraction == 0) {
      mms.doublegatefraction = simplegatefraction();
      mms.triplegatefraction = simplegatefraction() * simplegatefraction();
    }
    break;
  }

  mr.setGateFractions(mms.doublegatefraction, mms.triplegatefraction);

  if(verboseLevel >= 1) {
    std::cout << "PulsetrainManager Memory Usage:" << std::endl;
    std::cout << "Event number: " << events.size() << std::endl;
    std::cout << "Size of events vector: " << sizeof(NMSDetectedEvent) * events.capacity() / 1024  << " kB" << std::endl;
    std::cout << "Size of NMSSimpleEvent eventtime vector: " << sizeof(NMSSimpleEvent) * eventsprocessed.capacity() / 1024 << " kB" << std::endl;
  }
  return mr;
}

double NMSPulsetrainManager::getAverageLifetime() {
  if(events.size() > 0) {
    double totallifetimes;
    for(int i = 0; i < events.size(); i++) {
      totallifetimes += events[i].lifetime;
    }
    return 1.0 * totallifetimes / events.size();
  }
  else {
    return 0.;
  }

}

void NMSPulsetrainManager::saveSettingsToFile(std::string filename) {
  std::ofstream outfile;
  filename = filename + ".settings";

  if(verboseLevel >= 1) {
    std::cout << "Saving to file " << filename << std::endl;
  }
  outfile.open(filename.c_str());
  outfile << "predelay," << mms.predelay << std::endl;
  outfile << "gate," << mms.gate << std::endl;
  outfile << "adelay," << mms.adelay << std::endl;

  outfile << "registerLength," << mms.registerLength << std::endl;
  outfile << "registerPeriod," << mms.registerPeriod << std::endl;
  outfile << "quantizeRegisterBuffer," << mms.quantizeRegisterBuffer << std::endl;
  outfile << "quantizeRegister," << mms.quantizeRegister << std::endl;

  outfile << "derandomizePeriod," << mms.derandomizePeriod << std::endl;
  outfile << "derandomizeBuffer," << mms.derandomizeBuffer << std::endl;
  outfile << "derandomize," << mms.derandomize << std::endl;

  outfile << "efficiency," << mms.efficiency << std::endl;
  outfile << "dieawaymethod," << mms.dm << std::endl;
  outfile << "dieaway," << mms.dieaway << std::endl;
  outfile << "infinitegate," << mms.infinitegate << std::endl;
  outfile << "gatefractionmethod," << mms.gm << std::endl;
  outfile << "doublegatefraction," << mms.doublegatefraction << std::endl;
  outfile << "triplegatefraction," << mms.triplegatefraction << std::endl;

  outfile << "predeadtime," << mms.predeadtime << std::endl;
  outfile << "predeadtimeUpdating," << mms.predeadtimeUpdating << std::endl;
  outfile << "postdeadtime," << mms.postdeadtime << std::endl;
  outfile << "postdeadtimeUpdating," << mms.postdeadtimeUpdating << std::endl;

  outfile << "measurementLength," << mms.measurementLength << std::endl;

  outfile.close();
}

void NMSPulsetrainManager::savePulsetrainToFile(std::string filename) {
  std::ofstream outfile;
  double start;
  double stop;

  // Always sort before saving
  start = clock();
  if(verboseLevel >= 1) {
    std::cout << "Sorting data" << std::endl;
  }
  sort(events.begin(), events.end());
  stop = clock();
  if (verboseLevel >= 1) {
    std::cout << "WALL: Sort took " << ( stop - start ) / CLOCKS_PER_SEC << " second" << std::endl;
  }

  filename = filename + ".pulsetrain";
  if(verboseLevel >= 1) {
    std::cout << "Saving to file " << filename << std::endl;
  }
  outfile.open(filename.c_str());
  for(int i = 0; i < events.size(); i++) {
    outfile << std::setprecision(15) << events[i].eventtime << "," 
	    << events[i].lifetime << ","
	    << events[i].eventid << ","
	    // << events[i].info << std::endl;
	    << std::endl;
  }
  outfile.close();
}

void NMSPulsetrainManager::loadPulsetrainFromFile(std::string filename) {
  std::ifstream input;
  std::string test;
  std::string fieldstr;
  NMSDetectedEvent tempevent;

  filename = filename + ".pulsetrain";
  this->reset();

  if(verboseLevel >= 1) {
    std::cout << "Loading from file " << filename << std::endl;
  }
  if(file_exists(filename)) {
    input.open(filename.c_str());
    while ( input.good() ) {
      getline(input, test);
      //      input >> test;
      if(test.length() > 0) {
	std::istringstream iss( test );
	getline(iss, fieldstr, ',');
	std::istringstream( fieldstr ) >> tempevent.eventtime;
	getline(iss, fieldstr, ',');
	std::istringstream( fieldstr ) >> tempevent.lifetime;
	getline(iss, fieldstr, ',');
	std::istringstream( fieldstr ) >> tempevent.eventid;
	// getline(iss, fieldstr, ',');
	// tempevent.info = fieldstr;
	events.push_back(tempevent);
      }
    }
    input.close();
    processed = false;
  }
  else {
    std::cout << "ERROR: could not find file " << filename << std::endl;
  }

}

void NMSPulsetrainManager::loadPulsetrainFromShakesFile(std::string filename) {
  std::ifstream input;
  double test;
  this->reset();
  
  double start = 0;
  double stop = 0;
  start = clock();
  
  if(verboseLevel >= 1) {
    std::cout << "Loading from file (in shakes) " << filename << std::endl;
  }
  if(file_exists(filename)) {
    NMSDetectedEvent tempevent;
    //    tempevent.info = "";
    tempevent.lifetime = 0;
    tempevent.eventid = 0;
    input.open(filename.c_str());
    while ( input.good() ) {
      input >> test;
      tempevent.eventtime = test / 100;
      events.push_back(tempevent);
    }
    input.close();
    processed = false;
  }
  else {
    std::cout << "ERROR: could not find file " << filename << std::endl;
  }
  stop = clock();

  if(verboseLevel >= 1) {
    std::cout << "WALL: Calculation took " << ( stop - start ) / CLOCKS_PER_SEC << " second" << std::endl;
  }

}

void NMSPulsetrainManager::addPulsetrainFromFile(std::string filename, double offset) {
  // offset in seconds
  std::ifstream input;
  std::string test;
  std::string fieldstr;
  NMSDetectedEvent tempevent;

  filename = filename + ".pulsetrain";
  if(verboseLevel >= 1) {
    std::cout << "Loading from file " << filename << std::endl; 
    std::cout << "Offset " << offset << "s" << std::endl;
  }
  if(file_exists(filename)) {
    // Get biggest current eventid
    int eventoffset = -1;
    for(size_t i; i < events.size(); i++) {
      if(events[i].eventid > eventoffset)
	eventoffset = events[i].eventid;
    }
    eventoffset += 1;
    input.open(filename.c_str());
    while ( input.good() ) {
      getline(input, test);
      //      input >> test;
      if(test.length() > 0) {
	std::istringstream iss( test );
	getline(iss, fieldstr, ',');
	std::istringstream( fieldstr ) >> tempevent.eventtime;
	tempevent.eventtime += offset * 1000000;
	getline(iss, fieldstr, ',');
	std::istringstream( fieldstr ) >> tempevent.lifetime;
	getline(iss, fieldstr, ',');
	std::istringstream( fieldstr ) >> tempevent.eventid;
	tempevent.eventid += eventoffset;
	// getline(iss, fieldstr, ',');
	// tempevent.info = fieldstr;
	events.push_back(tempevent);
      }
    }
    input.close();
    processed = false;
  }
  else {
    std::cout << "ERROR: could not find file " << filename << std::endl;
  }

}

void NMSPulsetrainManager::splitPulsetrain(std::string filename, int splits) {
  if(mms.measurementLength <= 0) {
    std::cout << "Measurement time (runtime) not set, can not split" << std::endl;
    return;
  } 

  std::ofstream outfile;
  double start;
  double stop;

  // Always sort before saving
  start = clock();
  if(verboseLevel >= 1) {
    std::cout << "Sorting data" << std::endl;
  }
  sort(events.begin(), events.end());
  stop = clock();
  if (verboseLevel >= 1) {
    std::cout << "WALL: Sort took " << ( stop - start ) / CLOCKS_PER_SEC << " second" << std::endl;
  }

  if(verboseLevel >= 1) {
    std::cout << "Saving to file " << filename << std::endl;
  }

  int i = 0;
  int spl = 0;
  std::stringstream ss;
  ss << filename << "_" << spl  <<  ".pulsetrain";
  std::string fn = ss.str();

  outfile.open(fn.c_str());
  double limit = 1.0 * mms.measurementLength / splits;
  while(i < events.size() && events[i].eventtime <= mms.measurementLength) {
    if(events[i].eventtime >= limit) {
      spl++;
      limit = 1.0 * mms.measurementLength / splits * (spl + 1);
      ss.str("");
      ss.clear();
      ss << filename << "_" << spl  <<  ".pulsetrain";
      fn = ss.str();
      outfile.close();
      outfile.open(fn.c_str());
    }
    i++;
    outfile << std::setprecision(15) << events[i].eventtime << "," 
	    << events[i].lifetime << ","
	    << events[i].eventid << ","
	    // << events[i].info << std::endl;
	    << std::endl;
  }
  outfile.close();
}

void NMSPulsetrainManager::debug() {
  for(int i = 0; i < events.size(); i++) {
    std::cout.width(20);
    std::cout << events[i].eventtime;
  }
  std::cout << std::endl;

  std::cout << "=====================================================" << std::endl;
  std::cout << "** Processed Events" << std::endl;
  std::cout << "=====================================================" << std::endl;

  for(int i = 0; i < processedevents.size(); i++) {
    std::cout.width(20);
    std::cout << processedevents[i];
  }

  std::cout << std::endl;

}

void NMSPulsetrainManager::debugToFile(std::string filename) {
  std::ofstream original;
  std::ofstream processedfh;
  std::ofstream lost;

  std::string ofilename = filename;
  ofilename.append(".original.ptrain");
  std::string pfilename = filename;
  pfilename.append(".processedfh.ptrain");
  std::string lfilename = filename;
  lfilename.append(".lost.ptrain");

  original.open(ofilename.c_str());
  processedfh.open(pfilename.c_str());
  lost.open(lfilename.c_str());

  for(int i = 0; i < events.size(); i++) {
    original << std::setprecision(15) << events[i].eventtime << std::endl;
  }
  for(int i = 0; i < eventsprocessed.size(); i++) {
    processedfh << std::setprecision(15) << eventsprocessed[i].eventtime << std::endl;
  }
  // for(int i = 0; i < lostevents.size(); i++) {
  //   lost << std::setprecision(15) << lostevents[i] << std::endl;
  // }

  original.close();
  processedfh.close();
  lost.close();

}
