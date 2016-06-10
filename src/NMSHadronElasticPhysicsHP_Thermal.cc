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

// Modified file, based on G4HadronElasticPhysicsHP.cc, Geant4

// HP model for n with E < 20 MeV

#include "NMSHadronElasticPhysicsHP_Thermal.hh"
#include "G4SystemOfUnits.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4Neutron.hh"
#include "G4HadronicProcess.hh"
#include "G4HadronElastic.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPThermalScattering.hh"
#include "G4NeutronHPThermalScatteringData.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(NMSHadronElasticPhysicsHP_Thermal);

G4ThreadLocal G4bool NMSHadronElasticPhysicsHP_Thermal::wasActivated = false;
G4ThreadLocal G4HadronElasticPhysics* NMSHadronElasticPhysicsHP_Thermal::mainElasticBuilder = 0;

NMSHadronElasticPhysicsHP_Thermal::NMSHadronElasticPhysicsHP_Thermal(G4int ver)
  : G4VPhysicsConstructor("hElasticWEL_CHIPS_HP"), verbose(ver)
{
  if(verbose > 1) { 
    G4cout << "### NMSHadronElasticPhysicsHP_Thermal: " << GetPhysicsName() 
	   << G4endl; 
  }
  mainElasticBuilder = new G4HadronElasticPhysics(verbose);
}

NMSHadronElasticPhysicsHP_Thermal::~NMSHadronElasticPhysicsHP_Thermal()
{
  delete mainElasticBuilder;
}

void NMSHadronElasticPhysicsHP_Thermal::ConstructParticle()
{
  // G4cout << "G4HadronElasticPhysics::ConstructParticle" << G4endl;
  mainElasticBuilder->ConstructParticle();
}

void NMSHadronElasticPhysicsHP_Thermal::ConstructProcess()
{
  if(wasActivated) { return; }
  wasActivated = true;
  //Needed because this is a TLS object and this method is called by all threads
  if ( ! mainElasticBuilder ) mainElasticBuilder = new G4HadronElasticPhysics(verbose);
  mainElasticBuilder->ConstructProcess();

  mainElasticBuilder->GetNeutronModel()->SetMinEnergy(19.5*MeV);

  G4HadronicProcess* hel = mainElasticBuilder->GetNeutronProcess();
  G4NeutronHPElastic* hp = new G4NeutronHPElastic();
  hp->SetMinEnergy(4.0*eV);
  hel->RegisterMe(hp);
  hel->AddDataSet(new G4NeutronHPElasticData());
  
  G4NeutronHPThermalScattering* theThermalModel = new G4NeutronHPThermalScattering();
  theThermalModel->SetMaxEnergy(4.0*eV);
  hel->RegisterMe(theThermalModel);
  hel->AddDataSet(new G4NeutronHPThermalScatteringData());

  if (verbose > 1)
    {
      G4cout << "### NMSHadronElasticPhysicsHP_Thermal is constructed " << G4endl;
    }
}


