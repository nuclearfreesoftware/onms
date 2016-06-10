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


#include "G4PhysicalVolumeStore.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Ellipsoid.hh"

#include "NMSMaterialDecaySource.hh"

NMSMaterialDecaySource::NMSMaterialDecaySource() {
  verboseLevel = 0;
  
  Init();
  nmsmdsSettings->SetSourceMaterial(0);
}

NMSMaterialDecaySource::NMSMaterialDecaySource(G4Material* sourceMat) {
  verboseLevel = 0;

  Init();

  if(verboseLevel >= 2) {
    G4cout << "NMSMaterialDecaySource: Set source material to " << sourceMat << G4endl;
  }
  nmsmdsSettings->SetSourceMaterial(sourceMat);

  activity = -1;
  materialIntensity = -1;
  sourcesloaded = false;

}

NMSMaterialDecaySource::~NMSMaterialDecaySource() {
  G4cout << "Activity: " << activity * second << G4endl;
  sourceGenerator->ClearAll();
  delete sourceGenerator;
  delete nmsmdsMessenger;
  delete yieldCalculator;
}

void NMSMaterialDecaySource::Init() {
  nmsmdsSettings = NMSMaterialDecaySettings::Instance();
  nmsmdsMessenger = new NMSMaterialDecaySourceMessenger(this);

  yieldCalculator = new NMSANYield();
  
  startSourceTimeDistribution = 0;
  endSourceTimeDistribution = 0;

  spontaneousFissionNeutron = true;
  spontaneousFissionGamma = true;
  alphaDecay = false;
  betaDecay = false;
  alphaN = false;
  sourceGenerator = new NMSMultipleDecaySource();
  posGenerator = sourceGenerator->GetAllPosDist();
  posGenerator->SetPosDisType("Point");
  posGenerator->SetCentreCoords(G4ThreeVector(0,0,0));

  activity = 1;
  materialIntensity = 1;
  activeVolume = 1 * cm * cm * cm;

  sourcesloaded = false;

}

void NMSMaterialDecaySource::SetSourceMaterial(G4Material* sourceMat) {
  nmsmdsSettings->SetSourceMaterial(sourceMat);
  sourcesloaded = false;
  // define isotopes as in G4Material
}

std::vector< G4VPhysicalVolume* > NMSMaterialDecaySource::PhysicalVolumeRelationshipVector(G4VPhysicalVolume* physvol, G4LogicalVolume* logicvol) {
  for(G4int i = 0; i <logicvol->GetNoDaughters(); i++) {
    if(logicvol->GetDaughter(i) == physvol) {
      std::vector<G4VPhysicalVolume*> physhierarchy;
      physhierarchy.push_back(physvol);
      return physhierarchy;
    }
    std::vector<G4VPhysicalVolume*> below = PhysicalVolumeRelationshipVector(physvol, logicvol->GetDaughter(i)->GetLogicalVolume());
    if(below.size() > 0) {
      below.push_back(logicvol->GetDaughter(i));
      return below;
    }
  }
  std::vector<G4VPhysicalVolume*> a;
  return a;
}

void NMSMaterialDecaySource::SetSourceFromPhysicalVolumeName(G4String volname) {
  G4PhysicalVolumeStore * pvs = G4PhysicalVolumeStore::GetInstance(); 
  G4VPhysicalVolume * sourcevolumeP = pvs->GetVolume(volname);
  if(sourcevolumeP != 0) {
    G4LogicalVolume * sourcevolumeL = sourcevolumeP->GetLogicalVolume();
    if(sourcevolumeL != 0) {
      G4VSolid * sourcevolumeS = sourcevolumeL->GetSolid();
      if(sourcevolumeS != 0) {
	G4ThreeVector sourceposition = G4ThreeVector(0, 0, 0);
	std::vector<G4VPhysicalVolume* > vec = PhysicalVolumeRelationshipVector(sourcevolumeP, G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->GetWorldVolume()->GetLogicalVolume());
	if(vec.size() > 0) {
	  for(size_t i = 0; i < vec.size(); i++) {
	    sourceposition = sourceposition + vec[i]->GetTranslation();
	  }
	}
	else {
	  G4cout << "Error: Could not find the volume in the geometry tree - that might be due to an error." << G4endl;
	  exit(-1);
	}	
	if(verboseLevel >= 1) {
	  G4cout << "Found source volume solid" << G4endl;
	  //	sourcevolumeS->DumpInfo();
	  G4cout << "It is of the following material: " << sourcevolumeL->GetMaterial()->GetName() << G4endl;
	  G4cout << "Center of source volume: " << sourceposition / cm << G4endl;
	}
	if(G4Tubs * tubs = dynamic_cast< G4Tubs * >(sourcevolumeS)) {
	  SetSourceMaterial(sourcevolumeL->GetMaterial());
	  posGenerator->SetPosDisType("Volume");
	  if(verboseLevel >= 1) {
	    G4cout << "It is a cylinder!" << G4endl;
	    G4cout << "Volume: " << tubs->GetCubicVolume() / cm3 << " cm3" << G4endl;
	  }
	  posGenerator->SetPosDisShape("Cylinder");
	  posGenerator->SetRadius(tubs->GetOuterRadius());
	  posGenerator->SetRadius0(tubs->GetInnerRadius());
	  posGenerator->SetHalfZ(tubs->GetZHalfLength());
	  //position
	  posGenerator->SetCentreCoords(sourceposition);
	  //volume
	  SetActiveVolume(tubs->GetCubicVolume());
	}
	else if(G4Sphere * sphere = dynamic_cast< G4Sphere * >(sourcevolumeS)) {
	  SetSourceMaterial(sourcevolumeL->GetMaterial());
	  posGenerator->SetPosDisType("Volume");
	  if(verboseLevel >= 1) {
	    G4cout << "It is a sphere!" << G4endl;
	    G4cout << "Volume: " << sphere->GetCubicVolume() / cm3 << " cm3" << G4endl;
	  }
	  posGenerator->SetPosDisShape("Sphere");
	  posGenerator->SetRadius(sphere->GetOuterRadius());
	  posGenerator->SetRadius0(sphere->GetInnerRadius());
	  //position
	  posGenerator->SetCentreCoords(sourceposition);
	  //volume
	  SetActiveVolume(sphere->GetCubicVolume());
	}
	else if(G4Ellipsoid * ellipsoid = dynamic_cast< G4Ellipsoid * >(sourcevolumeS)) {
	  SetSourceMaterial(sourcevolumeL->GetMaterial());
	  posGenerator->SetPosDisType("Volume");
	  if(verboseLevel >= 1) {
	    G4cout << "It is an ellipsoid!" << G4endl;
	    G4cout << "Volume: " << ellipsoid->GetCubicVolume() / cm3 << " cm3" << G4endl;
	  }
	  posGenerator->SetPosDisShape("Ellipsoid");
	  posGenerator->SetHalfX(ellipsoid->GetSemiAxisMax(0));
	  posGenerator->SetHalfY(ellipsoid->GetSemiAxisMax(1));
	  posGenerator->SetHalfZ(ellipsoid->GetSemiAxisMax(2));
	  //position
	  posGenerator->SetCentreCoords(sourceposition);
	  //volume
	  SetActiveVolume(ellipsoid->GetCubicVolume());
	}
	else if(G4Box * box = dynamic_cast< G4Box * >(sourcevolumeS)) {
	  SetSourceMaterial(sourcevolumeL->GetMaterial());
	  posGenerator->SetPosDisType("Volume");
	  if(verboseLevel >= 1) {
	    G4cout << "It is a box!" << G4endl;
	    G4cout << "Volume: " << box->GetCubicVolume() / cm3 << " cm3" << G4endl;
	  }
	  posGenerator->SetPosDisShape("Para");
	  posGenerator->SetHalfX(box->GetZHalfLength());
	  posGenerator->SetHalfY(box->GetZHalfLength());
	  posGenerator->SetHalfZ(box->GetZHalfLength());
	  //position
	  posGenerator->SetCentreCoords(sourceposition);
	  //volume
	  SetActiveVolume(box->GetCubicVolume());
	  G4cout << "Box" << G4endl;
	}
	else {
	  G4cout << "Unknown source volume solid type. Will not set a source" << G4endl;
	}
	
      }
    }
  }
  else {
    G4cout << "Could not find the specified volume." << G4endl;
  }
}


void NMSMaterialDecaySource::SetActiveVolume(G4double vol){
  activeVolume = vol;
  if(verboseLevel >= 1) {
    G4cout << "MaterialDecaySource: Volume set to " << vol / cm3 << " cm3" << G4endl;
  }
  if(materialIntensity != -1.) {
    activity = materialIntensity * activeVolume;
  }
}

void NMSMaterialDecaySource::SetActivity(G4double act){
  activity = act;
}

G4double NMSMaterialDecaySource::GetActivity() {
  if(sourcesloaded == false) {
    LoadSources();
  }
  return activity;
}

void NMSMaterialDecaySource::SetSpontaneousFission(G4bool neutron = true, G4bool gamma = true) {
  if ( spontaneousFissionNeutron != neutron ) {
    spontaneousFissionNeutron = neutron;
    sourcesloaded = false;
  }
  if ( spontaneousFissionGamma != gamma ) {
    spontaneousFissionGamma = gamma;
    sourcesloaded = false;
  }
}

void NMSMaterialDecaySource::SetAlphaDecay(G4bool status) {
  if ( alphaDecay != status ) {
    alphaDecay = status;
    sourcesloaded = false;
  }
}

void NMSMaterialDecaySource::SetBetaDecay(G4bool status) {
  if ( betaDecay != status ) {
    betaDecay = status;
    sourcesloaded = false;
  }
}

void NMSMaterialDecaySource::SetAlphaNSource(G4bool status) {
  if ( alphaN != status) {
    alphaN = status;
    sourcesloaded = false;
  }
}

void NMSMaterialDecaySource::SetAlphaNFile(G4String filename){
  alphaNFilename = filename;
  sourcesloaded = false;
}

void NMSMaterialDecaySource::SetEnergyFile(G4String filename){
  energyFilename = filename;
  sourcesloaded = false;
}

void NMSMaterialDecaySource::GeneratePrimaryVertex(G4Event* anEvent) {
  //if(currentSourceMaterial == 0 && !alphaN) {
  if(nmsmdsSettings->GetSourceMaterial() == 0 && !alphaN) {
    G4cout << "NMSMaterialDecaySource - ERROR: No source material has been selected. Please use messenger commands to select a source material. If you haven't defined an appropriate source material in your detector construction, please do so as well" << G4endl;
    G4cout << "Aborting run!" << G4endl;
    G4RunManager::GetRunManager()-> AbortRun();
    return;
  }
  if(verboseLevel >= 2) {
    G4cout << "NMSMaterialDecaySource: GeneratePrimaryVertex" << G4endl; // VERBOSEF
  }
  if(!sourcesloaded or !nmsmdsSettings->GetSourcesloaded()) {
    if(verboseLevel >= 2) {
      G4cout << "NMSMaterialDecaySource: No loaded sources on call of GeneratePrimaryVertex, will load sources." << G4endl; // VERBOSEF
    }
    LoadSources();
  }
  G4double time = GetNextEventTime();
  //  G4AnalysisManager* anaman = G4AnalysisManager::Instance();
  //  anaman->FillH1(1, time * second);

  SetTime(time);
  sourceGenerator->GeneratePrimaryVertex(anEvent);
}

void NMSMaterialDecaySource::SetEventTimeLimits(G4double start, G4double end) {
  startSourceTimeDistribution = start;
  endSourceTimeDistribution = end;

  //fix check ranges
}

G4double NMSMaterialDecaySource::GetNextEventTime(){
  G4double time;
  G4double ran = G4UniformRand();

  if((startSourceTimeDistribution == 0) && (endSourceTimeDistribution == 0)) {
    time = - 1 / activity * std::log(ran) * second;
  }
  else {
    time = startSourceTimeDistribution + ran * (endSourceTimeDistribution - startSourceTimeDistribution);    
  }
  if(verboseLevel >= 2) {
    G4cout << " time of next source event: " << time / second << G4endl;
  }

  return time;
}

void NMSMaterialDecaySource::SetTime(G4double time){
  sourceGenerator->SetParticleTime(time);
}

void NMSMaterialDecaySource::SetVerboseLevel(G4int verbose) {
  verboseLevel = verbose;
  sourceGenerator->SetVerboseLevel(verbose);
  posGenerator->SetVerboseLevel(verbose);
}

void NMSMaterialDecaySource::LoadSources() {
  G4ParticleDefinition* ion = 0;
  G4IsotopeVector* ivec;
  G4double* relabvec;
  G4double elementatoms;
  G4int a; G4int z;
  G4double lifetime = -1;
  //G4double halflife;
  G4double sfbranching = -1;  G4double alphabranching = -1;  G4double betabranching = -1;
  G4double tmpactivity; G4double totalactivity = 0;
  G4double alphanactivity = 0;

  // Check Settings
  if(nmsmdsSettings->GetSourceMaterial() == 0) {
    if(!(alphaN == true and nmsmdsSettings->ANwithoutMaterialPossible())) {
      G4cout << "Error: No Source Material defined - Can not load sources." << G4endl;
      return; // FIX
    }
  }

  
  sourceGenerator->ClearAll();
  materialIntensity = 0;

  if(verboseLevel >= 1) {
    G4cout << "Load Source ====================================================" << G4endl;
    if(verboseLevel >= 2) {
      G4cout << "  Source Material" << G4endl;
      if(nmsmdsSettings->GetSourceMaterial() == 0) {
	G4cout << "Not defined!" << G4endl;
      }
      else {
	G4cout << nmsmdsSettings->GetSourceMaterial() << G4endl;
      }
    }
    G4cout << "Current source volume is: " << activeVolume / cm3 << " cm^3" << G4endl;
    G4cout << "Isotope data" << G4endl;
    DumpIsotopeData();
  }

  if(nmsmdsSettings->GetSourceMaterial() != 0) {
    const G4ElementVector* evec = nmsmdsSettings->GetSourceMaterial()->GetElementVector();
    // cycle through elements
    for(size_t i=0; i < evec->size(); i++) {
      elementatoms = nmsmdsSettings->GetSourceMaterial()->GetVecNbOfAtomsPerVolume()[i];
      ivec = (*evec)[i]->GetIsotopeVector();
      // Geant Documentation: a vector of relative abundances referring to such isotopes (where relative abundance means the number of atoms per volume)
      relabvec = (*evec)[i]->GetRelativeAbundanceVector();

      // cycle through isotopes
      for(size_t j=0; j < ivec->size(); j++) {
	a = (*ivec)[j]->GetN();
	z = (*ivec)[j]->GetZ();
	ion = G4IonTable::GetIonTable()->GetIon(z,a,0);
	lifetime = ion->GetPDGLifeTime();
	//	halflife = lifetime*log(2);
	if(lifetime != -1) {
	  sfbranching = GetSFBranching(a, z);
	  alphabranching = GetAlphaBranching(a, z);
	  betabranching = GetBetaBranching(a, z);
	  //adjust alpha and betabranching - geants own libraries do not include sf values up to now
	  if(sfbranching != 0.) {
	    alphabranching *= (1 - sfbranching);
	    betabranching *= (1 - sfbranching);
	  }
	  totalactivity += relabvec[j] * elementatoms / lifetime;
	  if(spontaneousFissionNeutron || spontaneousFissionGamma) {
	    // N * branching * lambda = N * branching * tao
	    tmpactivity = relabvec[j] * elementatoms * sfbranching / lifetime;
	    materialIntensity += tmpactivity;
	    NMSSFSingleIsoDecaySource * nsfds = new NMSSFSingleIsoDecaySource();
	    if(a == 252 and z == 98) {
	      nsfds->SetCf252n(nmsmdsSettings->GetCf252ndist(), nmsmdsSettings->GetCf252neng());
	    }
	    nsfds->setIsotope(z * 10000 + a * 10);
	    sourceGenerator->AddaSource(nsfds, tmpactivity);

	    // FIX: SF source should know about type
	    // if(spontaneousFissionNeutron && spontaneousFissionGamma) {
	    //   sourceGenerator->GetCurrentSource()->setDecayType(NMSDECAY_SF);
	    // }
	    // else if (spontaneousFissionNeutron) {
	    //   sourceGenerator->GetCurrentSource()->setDecayType(NMSDECAY_SF_N);
	    // }
	    // else if (spontaneousFissionGamma) {
	    //   sourceGenerator->GetCurrentSource()->setDecayType(NMSDECAY_SF_GAMMA);
	    // }
	  }
	  if(alphaDecay) {
	    // N * branching * lambda = N * branching * tao
	    tmpactivity = relabvec[j] * elementatoms * alphabranching / lifetime;
	    materialIntensity += tmpactivity;
	    NMSAlphaSingleIsoDecaySource * nads = new NMSAlphaSingleIsoDecaySource();
	    nads->setIsotope(z * 10000 + a * 10);
	    sourceGenerator->AddaSource(nads, tmpactivity);
	  }
	  if(betaDecay) {
	  }
	}
	else {
	  // nothing to do
        } // end lifetime
      } // end isotopes
    } // end elements
    if(alphaN) { // A
      tmpactivity = 1 / s / cm3;
      if(nmsmdsSettings->GetANActivityCalcMode() == AN_ACT_FIXED) {
	tmpactivity =  nmsmdsSettings->GetANActivity();
      }
      if(nmsmdsSettings->GetANActivityCalcMode() == AN_ACT_PRECALCMC) {
	// Get from file
      }
      if(nmsmdsSettings->GetANActivityCalcMode() == AN_ACT_CALC) {
	// Get from Yield
	G4double yield = yieldCalculator->fromSource(nmsmdsSettings->GetSourceMaterial());
	tmpactivity = yield;
      }
      // Add Source
      NMSANSource * nands = new NMSANSource();
      sourceGenerator->AddaSource(nands, tmpactivity);
      materialIntensity += tmpactivity;
      alphanactivity += tmpactivity;
      
    } // end alphaN part

  }
  else {
    if(alphaN) {
      // FIX
      NMSANSource * nands = new NMSANSource();
      sourceGenerator->AddaSource(nands, 1);
      // sourceGenerator->GetCurrentSource()->setDecayType(NMSDECAY_ALPHA_N);
      // sourceGenerator->GetCurrentSource()->setAlphaNFile(alphaNFilename);
    }
    else {
      G4cout << "Error: No Source Material defined - Can not load sources." << G4endl;
      return;
    }
  }
  if(!nmsmdsSettings->GetActivityFixed()) {
    activity = materialIntensity * activeVolume;
  }
  if(verboseLevel >= 1) {
    G4cout << "================================================================" << G4endl;
    G4cout << "Activity Analysis" << G4endl;
    G4cout << "Simulated source activity:                        " << activity * second << "/s" << G4endl;
    G4cout << "  - neutrons from (alpha, n):                     " << alphanactivity * activeVolume * second << "/s" << G4endl;
    G4cout << "  - spontaneouos fissions                         " << (materialIntensity - alphanactivity) * activeVolume * second << "/s" << G4endl;
    //FIX Level
    if(verboseLevel >= 2) {
      G4cout << "Simulated source activity per volume:            " << materialIntensity * second * ( cm * cm * cm ) << "/(s * cm^3)"<< G4endl;
      G4cout << "  - neutrons from (alpha, n):                   " << alphanactivity * second * ( cm * cm * cm ) << "/(s * cm^3)" << G4endl;
      G4cout << "Active volume:                                  " << activeVolume / ( cm * cm * cm ) << " cm^3 " << G4endl;
      G4cout << "Simulated activity per weight:                  " << materialIntensity / nmsmdsSettings->GetSourceMaterial()->GetDensity() * g * s << "/(s * g)" << G4endl;
      G4cout << "Total activity (except alpha, n) of the sample: " << activeVolume * totalactivity *second << "/s" << G4endl;
    }
    G4cout << "================================================================" << G4endl;
  }
  sourcesloaded = true;
  nmsmdsSettings->SetSourcesloaded(true);
}

G4double NMSMaterialDecaySource::GetSFBranching(G4int a, G4int z) {
  G4int zaid = z * 1000 + a;
  G4double br = 0;

  switch(zaid) {
  case 90232:
    br = 1.1e-11;
    // RadioactiveDecayData 4.3.1
    // consistent with
    // Nuclear Data Sheets Article from 2006 (10.1016/j.nds.2006.09.001)
    // and most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
      break;
  case 92232:
    br =     2.7e-14;
    // EXCEPTION: NOT RadioactiveDecayData 4.3.1!
    // taken from most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    // Nuclear Data Sheets Article from 2006 (10.1016/j.nds.2006.09.001) lists 8.5e-22,
    // which is most likely wrong. It lists the same half-life in years for SF as in
    // the ENSDF file, so calculation based on this give 2.7e-14 (actually 2.65e-14)
    // RadioactiveDecayData 4.3.1 lists 3e-14;
    break;
  case 92233:
    br = 6e-13;
    // RadioactiveDecayData 4.3.1
    // consistent with
    // most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    // (lists < 6e-13)
    // Nuclear Data Sheets Article from 2005 (10.1016/j.nds.2005.05.002) lists < 6e-11
    // which is most likely wrong. Using SF half-life listed in article gives ~ 6e-13
    break;
  case 92234:
    br = 1.64e-11;
    // EXCEPTION: NOT RadioactiveDecayData 4.3.1!
    // taken from Nuclear Data Sheets Article from 2007 (10.1016/j.nds.2007.02.003)
    // same value can be found in most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    // RadioactiveDecayData 4.3.1 lists 1.6e-11
    break;
  case 92235:
    br = 7.0e-11;
    // RadioactiveDecayData 4.3.1
    // consistent with
    // Nuclear Data Sheets Article from 2014 (10.1016/j.nds.2014.11.002)
    // and most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    break;
  case 92236:
    br = 9.4e-10;
    // RadioactiveDecayData 4.3.1
    // consistent with
    // Nuclear Data Sheets Article from 2006 (10.1016/j.nds.2006.09.002)
    // and most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    break;
  case 92238:
    br = 5.45e-7;
    // EXCEPTION: NOT RadioactiveDecayData 4.3.1!
    // taken from Nuclear Data Sheets Article from 2015 (10.1016/j.nds.2015.07.003)
    // same value can be found in most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    // RadioactiveDecayData 4.3.1 lists 5.5e-7, same as the NNDC Website/Chart of Nuclides (03.04.2016)
    break;
  case 93237:
    br = 2e-12;
    // RadioactiveDecayData 4.3.1
    // consistent (missing <) with
    // Nuclear Data Sheets Article from 2006 (10.1016/j.nds.2006.07.001)
    // and most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    break;
  case 94236:
    br = 1.9e-9;
    // RadioactiveDecayData 4.3.1
    // consistent with
    // Nuclear Data Sheets Article from 2006 (10.1016/j.nds.2006.09.002)
    // and most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    break;
  case 94238:
    br = 1.9e-9;
    // RadioactiveDecayData 4.3.1
    // consistent with
    // Nuclear Data Sheets Article from 2015 (10.1016/j.nds.2015.07.003)
    // and most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    break;
  case 94239:
    br = 3.1e-12;
    // EXCEPTION: NOT RadioactiveDecayData 4.3.1!
    // taken from Nuclear Data Sheets Article from 2014 (10.1016/j.nds.2014.11.003)
    // same value can be found in most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    // RadioactiveDecayData 4.3.1 lists 3e-12
    break;
  case 94240:
    br = 5.7e-8;
    // RadioactiveDecayData 4.3.1
    // consistent with
    // Nuclear Data Sheets Article from 2008 (10.1016/j.nds.2008.09.002)
    // and most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    break;
  case 94241:
    br = 2.4e-16;
    // EXCEPTION: NOT RadioactiveDecayData 4.3.1!
    // taken from Nuclear Data Sheets Article from 2015 (10.1016/j.nds.2015.11.004)
    // same value can be found in most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    // (ignoring <)
    // RadioactiveDecayData 4.3.1 lists 2e-16;
    break;
  case 94242:
    br = 5.5e-6;
    // RadioactiveDecayData 4.3.1
    // consistent with
    // Nuclear Data Sheets Article from 2002 (10.1006/ndsh.2002.0011)
    // and most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    break;
  case 95241:
    br = 3.6e-12;
    // EXCEPTION: NOT RadioactiveDecayData 4.3.1!
    // taken from Nuclear Data Sheets Article from 2015 (10.1016/j.nds.2015.11.004)
    // same value can be found in most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    // RadioactiveDecayData 4.3.1 lists 4e-12, same as the NNDC Website/Chart of Nuclides (03.04.2016)
    break;
  case 96242: 
    br = 6.2e-8;
    // RadioactiveDecayData 4.3.1
    // consistent with
    // Nuclear Data Sheets Article from 2002 (10.1006/ndsh.2002.0011)
    // and most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    break;
  case 96244:
    br = 1.37e-6;
    // EXCEPTION: NOT RadioactiveDecayData 4.3.1!
    // taken from Nuclear Data Sheets Article from 2003 (10.1006/ndsh.2003.0008)
    // same value can be found in most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    // RadioactiveDecayData 4.3.1 has 1.4e-6, as has the NNDC Website/Chart of Nuclides (03.04.2016)
    break;
  case 96246:
    br = 0.02615;
    // EXCEPTION: NOT RadioactiveDecayData 4.3.1!
    // taken from Nuclear Data Sheets Article from 2011 (10.1016/j.nds.2011.06.002)
    // same value can be found in most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    // RadioactiveDecayData 4.3.1 lists 0.0002999
    // 0.0003 is listed on the NNDC Website/Chart of Nuclides (03.04.2016)
    break;
  case 96248: 
    br = 0.0839;
    // RadioactiveDecayData 4.3.1
    // consistent with
    // Nuclear Data Sheets Article from 2014 (10.1016/j.nds.2014.11.004)
    // and most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    break;
  case 97249:
    br = 4.7e-10;
    // does not exist in RadioactiveDecayData 4.3.1!
    // taken from Nuclear Data Sheets Article from 2011 (10.1016/j.nds.2011.08.003)
    // same value can be found in most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    break;
  case 98250: 
    br = 0.00077;
    // EXCEPTION: NOT RadioactiveDecayData 4.3.1!
    // taken from Nuclear Data Sheets Article from 2001 (10.1006/ndsh.2001.0018)
    // same value can be found in most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    // RadioactiveDecayData 4.3.1 lists 0.00079958
    break;
  case 98246:
    br = 2.4e-6;
    // RadioactiveDecayData 4.3.1
    // consistent (rounding from 2.3999e-6) with
    // Nuclear Data Sheets Article from 2011 (10.1016/j.nds.2011.06.002)
    // and most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)    // ++
    break;
  case 98252: 
    br = 3.092e-2;
    // EXCEPTION: NOT RadioactiveDecayData 4.3.1!
    // taken from Nuclear Data Sheets Article from 2005 (10.1016/j.nds.2005.11.003)
    // same value can be found in most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    // RadioactiveDecayData 4.3.1 lists 0.030901
    break;
  case 98254: 
    br = 0.9969;
    // RadioactiveDecayData 4.3.1
    // consistent with
    // Nuclear Data Sheets Article from 2005 (10.1016/j.nds.2005.10.002)
    // and most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)    // ++
    break;
  case 98256: 
    br = 1.; 
    // RadioactiveDecayData 4.3.1
    // consistent with
    // Nuclear Data Sheets Article from 1999 (10.1006/ndsh.1999.0021)
    // and most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)    // ++
    break;
  case 100257:
    br = 0.0021;
    // RadioactiveDecayData 4.3.1
    // consistent with
    // Nuclear Data Sheets Article from 2013 (10.1016/j.nds.2013.08.002)
    // and most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)    // ++
    break;
  case 102252:
    br = 0.322;
    // does not exist in RadioactiveDecayData 4.3.1!
    // taken from Nuclear Data Sheets Article from 2005 (10.1016/j.nds.2005.11.003)
    // same value can be found in most recent ENSDF files (08.12.2015,
    // http://www.nndc.bnl.gov/ensarchivals/distributions/dist15/ensdf_151208_294.zip)
    break;
  default:
    if(verboseLevel >= 1) {
      G4cout << "No spontaneous fission data available for isotope " << zaid << G4endl;
    }
    br = 0;
  }
  return br;
}

G4double NMSMaterialDecaySource::GetAlphaBranching(G4int a, G4int z) {
  G4ParticleDefinition* ion =  G4IonTable::GetIonTable()->GetIon(z,a,0);
  G4RadioactiveDecay* decay = new G4RadioactiveDecay();
  G4DecayTable* dtable = decay->GetDecayTable(ion);
  G4VDecayChannel* dc;

  G4double br = 0;
  for(G4int i = 0; i < dtable->entries(); i++) {
    dc = dtable->GetDecayChannel(i);
    // for backwards compatibility
    {
      G4AlphaDecayChannel* testalpha = dynamic_cast< G4AlphaDecayChannel* >( dc );
      if( testalpha != 0) {
	br += dc->GetBR();
      }
    }
    // newer geant version have G4AlphaDecay as child class of G4NuclearDecayChannel
    {
      G4AlphaDecay* testalpha = dynamic_cast< G4AlphaDecay* >( dc );
      if( testalpha != 0 ) {
	br += dc->GetBR();
      }
    }

  }

  return br;
}

G4double NMSMaterialDecaySource::GetBetaBranching(G4int a, G4int z) {
  G4ParticleDefinition* ion =  G4IonTable::GetIonTable()->GetIon(z,a,0);
  G4RadioactiveDecay* decay = new G4RadioactiveDecay();
  G4DecayTable* dtable = decay->GetDecayTable(ion);
  G4VDecayChannel* dc;

  G4double br = 0;
  for(G4int i = 0; i < dtable->entries(); i++) {
    dc = dtable->GetDecayChannel(i);
    if( G4BetaMinusDecayChannel* testchannel = dynamic_cast< G4BetaMinusDecayChannel* >( dc ) ) {
      br += dc->GetBR();
    }
  }

  return br;

}

void NMSMaterialDecaySource::DumpIsotopeData() {
  G4ParticleDefinition* ion = 0;
  G4IsotopeVector* ivec;
  G4double* relabvec;
  G4double elementatoms;
  G4int a; G4int z;
  G4double lifetime = -1; G4double halflife;

  const G4ElementVector* evec = nmsmdsSettings->GetSourceMaterial()->GetElementVector();
  // cycle through elments
  for(size_t i=0; i < evec->size(); i++) {
    elementatoms = nmsmdsSettings->GetSourceMaterial()->GetVecNbOfAtomsPerVolume()[i];
    ivec = (*evec)[i]->GetIsotopeVector();
      // Geant Documenation: a vector of relative abundances referring to such isotopes (where relative abundance means the number of atoms per volume)
    relabvec = (*evec)[i]->GetRelativeAbundanceVector();

    // cycle through isotopes
    for(size_t j=0; j < ivec->size(); j++) {
      a = (*ivec)[j]->GetN();
      z = (*ivec)[j]->GetZ();
      ion = G4IonTable::GetIonTable()->GetIon(z,a,0);
      lifetime = ion->GetPDGLifeTime();
      halflife = lifetime*log(2);
      G4cout << G4endl;
      G4cout << "      Isotope data for Z: " << z << ", A: " << a << ", ZAID: " << z * 10000 + a * 10 << G4endl;
      G4cout << "        Fraction: ";
      G4cout.width(10);
      G4cout << relabvec[j];
      G4cout << " Atomic Density (1/cm^3): ";
      G4cout.width(10);
      G4cout << relabvec[j] * elementatoms * (cm3) << G4endl;
      if(lifetime != -1) {
	G4cout << "        Lifetime: ";
	G4cout.width(10);
	G4cout << lifetime / second;
	G4cout << "                Halflife: ";
	G4cout.width(10);
	G4cout << halflife / second;
	G4cout << G4endl;
	G4cout << "        Branching Ratios: " << G4endl;
	G4cout << "        SF: ";
	G4cout.width(8);
	G4cout << GetSFBranching(a, z);
	G4cout << " alpha: ";
	G4cout.width(8);
	G4cout << GetAlphaBranching(a, z);
	G4cout << " beta: ";
	G4cout.width(8);
	G4cout << GetBetaBranching(a, z) << G4endl;
      }
      else {
	G4cout << "        Stable isotope, will not be used as source." << G4endl;
      }
    }
  }
  G4cout << G4endl;
}

void NMSMaterialDecaySource::DumpSourceStatus() {
  G4cout << G4endl;
  G4cout << "================================================================================" << G4endl;
  G4cout << "NMSMaterialDecaySource --- Source Status Dump                              start" << G4endl;
  G4cout << "================================================================================" << G4endl;
  G4cout << "  Settings" << G4endl;

  G4cout << "  Alpha Decay:        " << ((alphaDecay) ? "active" : "inactive") << G4endl;
  G4cout << "  Beta Decay:         " << ((betaDecay) ? "active" : "inactive") << G4endl;
  G4cout << "  SF - neutrons:      " << ((spontaneousFissionNeutron) ? "active" : "inactive") << G4endl;
  G4cout << "  SF - gammas:        " << ((spontaneousFissionGamma) ? "active" : "inactive") << G4endl;
  G4cout << "  (alpha,n) neutrons: " << ((alphaN) ? "active" : "inactive") << G4endl;

  DumpIsotopeData();

  sourceGenerator->DumpSourceTable();

  G4cout << "  Source Volume: " << activeVolume / cm3 << " cm3" << G4endl;

  G4cout << "================================================================================" << G4endl;
  G4cout << "NMSMaterialDecaySource --- Source Status Dump                                end" << G4endl;
  G4cout << "================================================================================" << G4endl;
  G4cout << G4endl;

  
}

void NMSMaterialDecaySource::WriteANYield(G4String filename) {
  G4int bins = 1000;
  G4double emax = 10 * MeV;
  G4double estep = emax / bins;
  G4String fn = filename + ".alphanyield";
  std::ofstream outfile;
  outfile.open(fn.c_str());
  // output 0 
  outfile << 0 << " " << 0 << std::endl;
  for(G4int i = 1; i < bins; i++) {
    outfile << estep * i  << " " << yieldCalculator->fromEnergy(estep * i, nmsmdsSettings->GetSourceMaterial()) << std::endl;
  }
  outfile.close();
}

void NMSMaterialDecaySource::WriteANSpectrum(G4String filename) {
  for(G4int i = 0; i < sourceGenerator->GetNumberofSource(); i++) {
    sourceGenerator->SetCurrentSourceto(i);
    NMSNewSingleDecaySource * src = sourceGenerator->GetCurrentSource();
    NMSANSource * asrc = dynamic_cast< NMSANSource *>(src);
    if(asrc != 0) {
      asrc->WriteANSpectrum(filename);
    }
  }
}
