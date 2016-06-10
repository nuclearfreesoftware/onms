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

// Comments starting with ### give hints on intended conditionals, replaced in code by floating point equivalents

NMSMultiplicityResult NMSPulsetrainManager::ResultsFastForward() {

  NMSMultiplicityResult mr(mms.registerLength);
  mr.SetSettings(mms);

  double lastevent = 0;
  std::queue<double> predelayqueue;
  std::queue<double> longdelayqueue;
  std::queue<bool> shiftregister;
  int onshiftregister = 0;

  for (int i = 0; i < mms.registerLength; i++) {
    shiftregister.push(0);
  }

  while(!predelayqueue.empty()) {
    predelayqueue.pop();
  }
  if( processedevents.size() == 0) {
    return mr;
  }

  double t = 0;
  double t2 = 0;

  double epsilon = mms.registerPeriod / 10000;

  for (int i = 0; i < processedevents.size(); i++) {
    t2 = lastevent;

    // fast forward shift register

    predelayqueue.push(processedevents[i] + mms.predelay);
    longdelayqueue.push(processedevents[i] + mms.adelay);

    while(
	  //	  (onshiftregister != 0 || !predelayqueue.empty() || !longdelayqueue.empty())
	  // ### ( t < processedevents[i] )
	   (( processedevents[i] - t) > epsilon)) {

      if(shiftregister.front()) {
	onshiftregister--;
      }
      shiftregister.pop();

      if((!predelayqueue.empty())
	 // ### predelayqueue.front() == t
	 && (
	     ((t - predelayqueue.front()) > -epsilon) && ((t - predelayqueue.front()) < epsilon ))
	 ) {
	shiftregister.push(1);
	onshiftregister++;
	predelayqueue.pop();
      }
      else {
	shiftregister.push(0);
      }
      if(!longdelayqueue.empty()) {
	// ### t > longdelayqueue.front()
	if((t - longdelayqueue.front()) > epsilon) {
	  //	cout << onshiftregister << endl;
	  mr.addA(onshiftregister, 1);
	  longdelayqueue.pop();
	}
      }


      t += mms.registerPeriod;

    }
    t = processedevents[i];
    if(shiftregister.size() != mms.registerLength) {
      std::cout << shiftregister.size() << std::endl;
    }

    // Safety checks
    if(((t - processedevents[i]) < - epsilon)) {
      std::cout << "Fast forward not complete... (" << t - processedevents[i] << ")" << std::endl;
    }
    if(((t - processedevents[i]) >  epsilon)) {
      std::cout << "Fast forward not complete... (" << t - processedevents[i] << ")" << std::endl;
    }


    mr.addRA(onshiftregister, 1);
    //    lastevent = processedevents[i];
  }

  for (double t = processedevents.back(); t <= processedevents.back() + mms.adelay + mms.registerPeriod; t += mms.registerPeriod) {
    //    cout << t << std::endl;
    if(!longdelayqueue.empty()) {
      if(longdelayqueue.front() < t) {
	//	cout << onshiftregister << std::endl;
	mr.addA(onshiftregister, 1);
	longdelayqueue.pop();
      }
    }
    if(shiftregister.front()) {
      onshiftregister--;
    }
    shiftregister.pop();
    if((!predelayqueue.empty()) && (predelayqueue.front() < t)) {
      shiftregister.push(1);
      onshiftregister++;
      predelayqueue.pop();
    }
    else {
      shiftregister.push(0);
    }
    //    cout << setw(5) << t << setw(5) << onshiftregister << std::endl;

  }

  mr.setLastEvent(processedevents.back());
  return mr;

}
