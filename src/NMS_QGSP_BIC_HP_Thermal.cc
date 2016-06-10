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

#include <iomanip>   

#include "globals.hh"
#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"

#include "G4DataQuestionaire.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"

#include "NMS_QGSP_BIC_HP_Thermal.hh"
#include "NMSHadronElasticPhysicsHP_Thermal.hh"

NMS_QGSP_BIC_HP_Thermal::NMS_QGSP_BIC_HP_Thermal(G4int ver):  G4VModularPhysicsList()
{
  // default cut value  (1.0mm) 
  // defaultCutValue = 1.0*CLHEP::mm;

  G4DataQuestionaire it(photon, neutron);
  G4cout << "<<< Geant4 Physics List simulation engine: Thermal n-scattering (based on QGSP_BIC_HP 2.0)"<<G4endl;
  G4cout <<G4endl;

  this->defaultCutValue = 0.7*CLHEP::mm;  
  this->SetVerboseLevel(ver);


  // EM Physics
  this->RegisterPhysics( new G4EmStandardPhysics(ver) );

  // Synchroton Radiation & GN Physics
  this->RegisterPhysics( new G4EmExtraPhysics(ver) );

  // Decays
  this->RegisterPhysics( new G4DecayPhysics(ver) );

  // This has been changed!
  // Hadron Elastic scattering
  this->RegisterPhysics( new NMSHadronElasticPhysicsHP_Thermal(ver) );

  // Hadron Physics
  this->RegisterPhysics(  new G4HadronPhysicsQGSP_BIC_HP(ver));

  // Stopping Physics
  this->RegisterPhysics( new G4StoppingPhysics(ver) );

  // Ion Physics
  this->RegisterPhysics( new G4IonPhysics(ver));
  
}

NMS_QGSP_BIC_HP_Thermal::~NMS_QGSP_BIC_HP_Thermal()
{
}

void NMS_QGSP_BIC_HP_Thermal::SetCuts()
{
  if (this->verboseLevel >1){
    G4cout << "NMS_QGSP_BIC_HP_Thermal::SetCuts:";
  }  
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 

  this->SetCutsWithDefault();   
  
//  if (this->verboseLevel >0)
//    G4VUserPhysicsList::DumpCutValuesTable();  

}



