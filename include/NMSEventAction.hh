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

#ifndef NMSEventAction_h
#define NMSEventAction_h 1

#include "globals.hh"
#include "G4UserEventAction.hh"
#include "G4Event.hh"

#include "NMSPrimaryUserInformation.hh"

class G4Event;

class NMSEventAction : public G4UserEventAction {
public:
  NMSEventAction();

  inline virtual void BeginOfEventAction(const G4Event*);

  void SetPrintModulo(G4int val);

private:
  G4int printModulo;

};

inline void NMSEventAction::BeginOfEventAction(const G4Event* evt) {
  G4PrimaryVertex* pv = evt->GetPrimaryVertex();
  for(G4int i = 0; i < pv->GetNumberOfParticle(); i++) {
    NMSPrimaryUserInformation* pui = (NMSPrimaryUserInformation*) pv->GetPrimary(i)->GetUserInformation();
    if(pui != 0) {
      if(pui->GetOrigin() == ORIGIN_ALPHA_N) {
	G4cout << " From alpha  n" << G4endl;
      }
    }
  }
  G4int nEvt = evt->GetEventID();
  if(G4int(nEvt/printModulo)*printModulo == nEvt) {
      G4cout << "** Event # " << nEvt << " started" << G4endl;
  }
}


#endif /* NMSEventAction_h */
