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

NMSMultiplicityResult NMSPulsetrainManager::ResultsAllTimeSteps() {

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

  double epsilon = mms.registerPeriod / 100;
  double oldt = 0;
  int i = 0;
  double t = 0;
  for(int j = 0; j <= (processedevents.back()/mms.registerPeriod); j+=1) {
    t = j * mms.registerPeriod;
    if(processedevents[i] >= oldt && processedevents[i] <= t) {

      //      mr.addRA(onshiftregister,1);
    }

    // if (((t - processedevents[i]) < epsilon) && ((t - processedevents[i]) > -epsilon)){
    // 	mr.addRA(onshiftregister, 1);
    // 	predelayqueue.push(t + predelay);
    // 	i++;
    //   }

      if(shiftregister.front()) {
	onshiftregister--;
      }
      shiftregister.pop();

      //      if (((t - predelayqueue.front()) < epsilon) && ((t - predelayqueue.front()) > -epsilon)){
      // t >= predelayqueue.front()
      if (!predelayqueue.empty() &&  ((t - predelayqueue.front()) > -epsilon)) {
	predelayqueue.pop();
	shiftregister.push(1);
	onshiftregister++;
      }
      else {
	shiftregister.push(0);
      }

      if(!longdelayqueue.empty() && ((t - longdelayqueue.front()) > -epsilon)) {
	mr.addA(onshiftregister,1 );
	longdelayqueue.pop();
      }

      if (((t - processedevents[i]) < epsilon) && ((t - processedevents[i]) > -epsilon)){
	mr.addRA(onshiftregister, 1);
	predelayqueue.push(t + mms.predelay);
	longdelayqueue.push(t + mms.adelay);
	i++;
      }

      /*      // ### (processedevents[i] + predelay) < t
      if((t - (processedevents[i] + predelay)) > epsilon) {

	shiftregister.push(1);
	onshiftregister++;

	//	i++;
      }
      else {
	shiftregister.push(0);
	}*/

      oldt = t;

  }
  //  lastevent = processedevents[9999];
  //  cout << setprecision(14) << lastevent << " " << processedevents[9999] << std::endl;

  mr.setLastEvent(processedevents.back());
  return mr;

}
