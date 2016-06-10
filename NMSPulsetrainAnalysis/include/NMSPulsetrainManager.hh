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

/** \brief Class to store and evaluate pulse trains
 *
 *  This class should handle pulse trains, e.g. produced by neutron counting devices.
 *  Therefore, different methods are available, and different options
 *  regarding predelay, register length and register frequency can be set.
 */

#ifndef NMSPulsetrainManager_h
#define NMSPulsetrainManager_h 1

#include <vector>
#include <sys/stat.h>
#include <sstream>

#include "NMSMultiplicityResult.hh"
#include "NMSDetectedEvent.hh"
#include "NMSDetectedEventSet.hh"
#include "NMSMultiplicityMeasurementSetting.hh"
#include "NMSSimpleEvent.hh"

#define RESULTS_ALL_EVENTS_NO_SHIFT 1
#define RESULTS_ALL_EVENTS_NO_SHIFT_FORWARD 2
#define RESULTS_ALL_TIME_STEPS_SHIFT 3
#define RESULTS_FAST_FORWARD_SHIFT 4

#define RESULTS_LL_ALL_EVENTS_NO_SHIFT_BACKWARD 11
#define RESULTS_LL_ALL_EVENTS_NO_SHIFT_FORWARD 12
#define RESULTS_LL_SINGLE_EVENT_NO_SHIFT_BACKWARD 13
#define RESULTS_LL_SINGLE_EVENT_NO_SHIFT_FORWARD 14

#define RESULTS_LL_ALL_EVENTS_NO_SHIFT_FORWARD_FULL 22

struct NMSGateFraction {
  double doubles;
  double triples;
};

class NMSPulsetrainManager {

public:
  NMSPulsetrainManager(int registerLengthInit = 128);

  ~NMSPulsetrainManager();

  void setVerboseLevel(int vlevel);

  inline void add(NMSDetectedEvent event) {
    events.push_back(event);
  };
  void sortOriginal();
  void reset();

  void setMultiplicityLength(int ml);
  
  void setEfficiency(double eff);
  void setMeasurementLength(double newtime);
  void setMeasurementLengthFromPulsetrain();
  
  void setPredelay(double t);
  void setGate(double t);
  void setAdelay(double t);
  
  void setRegisterLength(int length);
  void setRegisterPeriod(double t);
  void setFrequency(double frequency);
  void setRegisterQuantizationBuffer(int buffer = 1);
  void setRegisterQuantization(bool rq);

  void setDieawayMethod(DieawayMethod dm);
  void setDieaway(double da);
  void setInfiniteGate(double g);
  void setGatefractionMethod(GatefractionMethod gm);
  void setGateFractions(double fd, double ft);

  void setDerandomizePeriod(double t);
  void setDerandomizeBuffer(int storage = 1);
  void setDerandomize(bool dr);

  void setPreDeadtime(double);
  void setPreDeadtimeUpdating(bool);
  void setPostDeadtime(double);
  void setPostDeadtimeUpdating(bool);

  NMSMultiplicityMeasurementSetting getSettings() { return mms; };

  void deadtime(unsigned long long deadtime, bool updating, std::vector<NMSSimpleEvent>& lost);
  void quantize(unsigned long long period, int buffer, std::vector<NMSSimpleEvent>& lost);
  void processEvents();

  NMSMultiplicityResult getLongGateResults(int mode);
  NMSGateFraction gatefractionEventNo();

  NMSMultiplicityResult getResults(int mode = RESULTS_LL_ALL_EVENTS_NO_SHIFT_FORWARD);
  int getEventCount() { return events.size(); };
  double getAverageLifetime();

  void saveSettingsToFile(std::string filename);
  void savePulsetrainToFile(std::string filename);
  void loadPulsetrainFromFile(std::string filename);
  void loadPulsetrainFromShakesFile(std::string filename);
  void addPulsetrainFromFile(std::string filename, double offset = 0);
  void splitPulsetrain(std::string filename, int splits);
  
  void debug();
  void debugToFile(std::string filename);

private:
  NMSMultiplicityResult ResultsAllTimeSteps();
  NMSMultiplicityResult ResultsAllEvents();
  NMSMultiplicityResult ResultsAllEventsForward();
  NMSMultiplicityResult ResultsFastForward();
  NMSMultiplicityResult LLResultsAllEventsBackward();
  NMSMultiplicityResult LLResultsAllEventsForward();
  NMSMultiplicityResult LLResultsAllEventsForwardF();
  NMSMultiplicityResult LLResultsSingleEventBackward(bool infinitegate = false);
  NMSMultiplicityResult LLResultsSingleEventForward(bool infinitegate = false);
  
  double simplegatefraction();

  inline bool file_exists(std::string filename) {
    struct stat buffer;
    return (stat (filename.c_str(), &buffer) == 0);
  }

private:

  NMSMultiplicityMeasurementSetting mms;

  bool processed;
  std::vector<double> processedevents;
  std::vector<double> lostevents;
  std::vector<NMSSimpleEvent> eventsprocessed;
  std::vector<NMSSimpleEvent> eventslostpredeadtime;
  std::vector<NMSSimpleEvent> eventslostpostdeadtime;
  std::vector<NMSSimpleEvent> eventslostderandomize;
  std::vector<NMSSimpleEvent> eventslostregisterquantization;

  int verboseLevel;
  NMSDetectedEventSet events;

};

#endif
