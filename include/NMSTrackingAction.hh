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

#ifndef NMSTrackingAction_h
#define NMSTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "G4ParticleDefinition.hh"

#include "NMSAnalysisManager.hh"
class NMSAnalysisManager;

class NMSTrackingAction : public G4UserTrackingAction {
public:
  NMSTrackingAction(NMSAnalysisManager*);
  virtual ~NMSTrackingAction();
  virtual void PreUserTrackingAction(const G4Track*);
  virtual void PostUserTrackingAction(const G4Track*);

  void SetRunMode(RunMode newrm = RUNMODE_HE_DETECTOR);
  RunMode GetRunMode() { return rm; };

private:
  NMSAnalysisManager* nmsam;
  RunMode rm;
};



#endif /* NMSTrackingAction_h */
