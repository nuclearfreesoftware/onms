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

#include "NMSSteppingAction.hh"

#include "G4EventManager.hh"

NMSSteppingAction::NMSSteppingAction(NMSAnalysisManager* newnmsam)
{
  nmsam = newnmsam;
}

NMSSteppingAction::~NMSSteppingAction()
{

}

void NMSSteppingAction::UserSteppingAction(const G4Step * theStep)
{
  G4Track * theTrack = theStep->GetTrack();
  // check if it is not alive anymore, otherwise return
  if(theTrack->GetTrackStatus()==fAlive) {
    return;
  }

  G4ParticleDefinition * particleType = theTrack->GetDefinition();
  if(particleType==G4Neutron::NeutronDefinition())
  {
    G4String lostneutronvolumename = nmsam->GetLostNeutronVolume();

    G4StepPoint * thePostPoint = theStep->GetPostStepPoint();
    G4VPhysicalVolume * thePostPV = thePostPoint->GetPhysicalVolume();
    G4String thePostPVname = "Unknown";
    if(thePostPV != 0) {
      thePostPVname = thePostPV->GetName();
    }
    G4String processname = "Unknown";
    if(thePostPoint->GetProcessDefinedStep() != 0) {
      processname = thePostPoint->GetProcessDefinedStep()->GetProcessName();
    }
    G4int eid = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID() + nmsam->GetEventOffset();
    if(thePostPVname(0,nmsam->GetDetectorVolume().length()) == nmsam->GetDetectorVolume()) {
      if((processname != "inelastic") && (processname != "NeutronInelastic") && (processname != "neutronInelastic")) {
        G4cout << "** ERROR: Unknown process in detector volume, event will be counted anyway: " << processname << G4endl;
      }
      // convert to non-unit value in microseconds
      nmsam->addDetectorNeutron(NMSDetectedEvent(theTrack->GetGlobalTime() / microsecond, theTrack->GetLocalTime() / microsecond, eid));
    }
    else if(thePostPVname(0,lostneutronvolumename.length()) == lostneutronvolumename) {
      nmsam->addLostNeutron(NMSDetectedEvent(theTrack->GetGlobalTime() / microsecond, theTrack->GetLocalTime() / microsecond, eid));
    }
    //Only for debug purposes
    /* else {
      G4cout << "Neutron absorbed in " << thePostPVname << " with '" << processname << "'" << G4endl;
      nmsam->addSomeNeutron(NMSDetectedEvent(theTrack->GetGlobalTime() / microsecond, theTrack->GetLocalTime() / microsecond, eid, thePostPVname));
      }*/

  }
}
