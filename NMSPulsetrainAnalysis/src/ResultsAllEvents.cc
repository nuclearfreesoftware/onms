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

#include <queue>
#include <iomanip>

// Comments starting with ### give hints on intended conditionals, replaced in code by floating point equivalents

NMSMultiplicityResult NMSPulsetrainManager::ResultsAllEvents() {

  NMSMultiplicityResult mr(mms.registerLength);
  mr.SetSettings(mms);

  double lastevent = 0;
  std::queue<double> predelayqueue;
  std::queue<double> longdelayqueue;
  std::queue<bool> shiftregister;
  int onshiftregister = 0;

  double start;
  double stop;

  start = clock();

  for (int i = 0; i < mms.registerLength; i++) {
    shiftregister.push(0);
  }
  while(!predelayqueue.empty()) {
    predelayqueue.pop();
  }
  if( processedevents.size() == 0) {
    return mr;
  }

  double epsilon = mms.registerPeriod / 1000;
  int j = 0;
  int pcount = 0;
  int oldj = 0;
  int shiftcount = 0;
  int verbositySteps = processedevents.size() / 50;
  if(verbositySteps == 0) {
    verbositySteps = 1;
  }

  if(verboseLevel >= 1) {
    std::cout << "Calculating";
    std::cout.flush();
  }

  // ofstream tmpdebug;
  // string filename = "4-debug";

  // tmpdebug.open(filename.c_str());

  for(int i = 0; i < processedevents.size(); i++) {
    j = i - 1;
    if(i % verbositySteps == 0) {
      if(verboseLevel >= 1) {
	std::cout << ".";
	std::cout.flush();
      }
    }
    pcount = 0;
    shiftcount = 0;

    // while(j >= 0) {
    //   pcount = j;
    //   if(processedevents[j] <= (processedevents[i] - predelay)
    // 	 && processedevents[j] >= (processedevents[i] - predelay - (registerPeriod * registerLength))) {
    // 	shiftcount++;
    //   }
    //   if(processedevents[j] < (processedevents[i] - predelay - (registerPeriod * registerLength))) {
    // 	oldj = j;
    // 	j = -1;
    //   }
    //   else {
    // 	j--;
    //   }
    // }


    // while(j >= 0
    // // ### processedevents[i] - (registerPeriod * registerLength) < processedevents[j]
    // 	  && (processedevents[j] >= (processedevents[i] - predelay - (registerPeriod * registerLength)) ) ) {
    //   j--;
    //   shiftcount ++;
    // }

    while(j >= 0
    	  // ### processedevents[i] - predelay < processedevents[j]
    	  && (processedevents[j] - (processedevents[i] - mms.predelay) > epsilon) ) {
      j--;
      pcount++;
    }

    while(j >= 0
    // ### processedevents[i] - (registerPeriod * registerLength) < processedevents[j]
    	  && (processedevents[j] - (processedevents[i] - mms.predelay - (mms.registerPeriod * mms.registerLength)) > epsilon) ) {
      oldj = j;
      j--;
      shiftcount++;
    }

    // while(j >= 0
    // 	  // ### processedevents[j] >= (processedevents[i] - predelay - (registerPeriod * registerLength))
    // 	  && ((processedevents[j] - (processedevents[i] - predelay - (registerPeriod * registerLength))) > -epsilon)
    // 	  && (processedevents.back() - (processedevents[i] + adelay) > -epsilon)
    // 	  ){
    //   j--;
    //   // ### processedevents[j] <= (processedevents[i] - predelay)
    //   if((processedevents[i] - predelay) - processedevents[j] > -epsilon) {
    // 	shiftcount++;
    //   }
    // }
    // j++;
    mr.addRA(shiftcount, 1);
    // if(shiftcount == 4) {
    //   tmpdebug << setprecision(14) << processedevents[i] << " " << (processedevents[i] - predelay) << " " << (processedevents[i] - predelay - (registerPeriod * registerLength)) << std::endl;
    //   for(int k = oldj; k <= i; k++) {
    // 	tmpdebug << setw(15) << processedevents[k] << " ";
    //   }
    //   tmpdebug << std::endl;
    // }

    pcount = 0;
    j = i + 1;
    while(j < processedevents.size()
	  // ### processedevents[i] + adelay > processedevents[j]
	  && ( ((processedevents[i] + mms.adelay) - processedevents[j]) > epsilon) ) {
      j++;
      pcount++;
    }

    while(j < processedevents.size()
	  // ### processedevents[i] + delay + gatelength > processedevents[j]
	  && ( ((processedevents[i] + mms.adelay + (mms.registerPeriod * mms.registerLength))- processedevents[j]) > epsilon) ) {
      j++;
    }
    j--;
    mr.addA((j - i - pcount), 1);

  }

  // tmpdebug.close();

  if(verboseLevel >= 1) {
    std::cout << std::endl;
  }

  mr.setLastEvent(processedevents.back());

  stop = clock();
  if(verboseLevel >= 1) {
    std::cout << "WALL: Calculation took " << ( stop - start ) / CLOCKS_PER_SEC << " second" << std::endl;
  }

  return mr;

}
