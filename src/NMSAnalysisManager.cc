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
#include <fstream>


#include "NMSAnalysisManager.hh"

NMSAnalysisManager* NMSAnalysisManager::amInstance = 0;


NMSAnalysisManager::NMSAnalysisManager() {
  nmsamm = new NMSAnalysisManagerMessenger(this);
  nmspm = new NMSPulsetrainManager();

  outputfilename = "nms-output";
  exportAfterRun = false;
  exportPulsetrain = false;
  pulsetrainAnalysisMode = RESULTS_LL_ALL_EVENTS_NO_SHIFT_FORWARD;
  runtimeFromPulsetrain = false;
  splits = 10;
  rm = RUNMODE_HE_DETECTOR;
  sourceNeutronSpectrum.clear();
  eventIDoffset = 0;
}

NMSAnalysisManager::~NMSAnalysisManager() {
  delete nmspm;
  delete nmsamm;
}

void NMSAnalysisManager::reset() {
  sourceNeutronCount = 0;
  totalSourceNeutronEnergy = 0;
  sourceAlphaCount = 0;
  totalSourceAlphaEnergy = 0;
  sourceGammaCount = 0;
  totalSourceGammaEnergy = 0;
  sourceSFcount = 0;
  sourceNeutronSpectrum.clear();

  secondaryNeutronCount = 0;
  totalSecondaryNeutronEnergy = 0;

  lostneutronlist.clear();
  nmspm->reset();
  
  calculated = false;
}


NMSAnalysisManager* NMSAnalysisManager::GetInstance() {
  if(amInstance == 0) {
    amInstance = new NMSAnalysisManager();
  }
  return amInstance;
}

void NMSAnalysisManager::saveAlphaNToFile(G4String filename) {
  alphanset.saveToFile(filename);
}

void NMSAnalysisManager::saveNeutronSpectrumToFile(G4String filename) {
  std::ofstream outfile;
  G4String  nspectrumfilename;
  nspectrumfilename = filename + ".nspectrum";
  outfile.open(nspectrumfilename.c_str());
  for(int i = 0; i < sourceNeutronSpectrum.size(); i++) {
    outfile  << sourceNeutronSpectrum[i] / MeV << std::endl;
  }
  outfile.close();
}

G4double NMSAnalysisManager::getSourceNeutronsEnergy() {
  if(sourceNeutronCount > 0) {
    return totalSourceNeutronEnergy / sourceNeutronCount;
  }
  else {
    return 0;
  }
}

G4double NMSAnalysisManager::getAverageSourceNeutronsPerTime() {
  G4double runtime = static_cast<NMSRunManager*>(NMSRunManager::GetRunManager())->GetRuntime();
  if(runtime > 0) {
    return sourceNeutronCount / runtime * s;
  }
  else {
    return 0;
  }
}

G4double NMSAnalysisManager::getSourceAlphasEnergy() {
  if(sourceAlphaCount > 0) {
    return totalSourceAlphaEnergy / sourceAlphaCount;
  }
  else {
    return 0;
  }
}

G4double NMSAnalysisManager::getSourceGammasEnergy() {
  if(sourceGammaCount > 0) {
    return totalSourceGammaEnergy / sourceGammaCount;
  }
  else {
    return 0;
  }
}

G4double NMSAnalysisManager::getAverageNu() {
  return 0;
}

G4double NMSAnalysisManager::getSecondaryNeutronsEnergy() {
  if(secondaryNeutronCount > 0) {
    return totalSecondaryNeutronEnergy / secondaryNeutronCount;
  }
  else {
    return 0;
  }
}

G4double NMSAnalysisManager::getDetectorNeutronEfficiency( NeutronEfficiencyOutput neo ) {
  if(sourceNeutronCount > 0) {
    if( neo == EFFICIENCY_SOURCE) {
      return 1.0 * nmspm->getEventCount() / sourceNeutronCount;
    }
    if( neo == EFFICIENCY_TOTAL) {
      return 1.0 * nmspm->getEventCount() / (sourceNeutronCount + secondaryNeutronCount);
    }
  }
  return 0;
}

G4double NMSAnalysisManager::getSingles(SDTOutput output) {
  if(!calculated) {
    calculateResults();
  }
  if(output == SDTperSecond) {
    G4double runtime = static_cast<NMSRunManager*>(NMSRunManager::GetRunManager())->GetRuntime();
    if(runtime == 0) {
      runtime = 1 * s;
    }
    return nmsmr.singles(1 / runtime * s);
  }
  else {
    return nmsmr.singles(1);
  }
}

G4double NMSAnalysisManager::getDoubles(SDTOutput output) {
  if(!calculated) {
    calculateResults();
  }
  if(output == SDTperSecond) {
    G4double runtime = static_cast<NMSRunManager*>(NMSRunManager::GetRunManager())->GetRuntime();
    if(runtime == 0) {
      runtime = 1 * s;
    }
    return nmsmr.doubles(1 / runtime * s);
  }
  else {
    return nmsmr.doubles(1);
  }
}

G4double NMSAnalysisManager::getTriples(SDTOutput output) {
  if(!calculated) {
    calculateResults();
  }
  if(output == SDTperSecond) {
    G4double runtime = static_cast<NMSRunManager*>(NMSRunManager::GetRunManager())->GetRuntime();
    if(runtime == 0) {
      runtime = 1 * s;
    }
    return nmsmr.triples(1 / runtime * s);
  }
  else {
    return nmsmr.triples(1);
  }

}

G4double NMSAnalysisManager::getEffPuMass() {
  if(!calculated) {
    calculateResults();
  }
  return nmsmr.pu240eff();
}

G4double NMSAnalysisManager::getTotalPuMass() {  // FIXME add Material
  if(!calculated) {
    calculateResults();
  }
  G4Material * sourcemat = NMSMaterialDecaySettings::Instance()->GetSourceMaterial();
  if(sourcemat == 0) {
    return 0.0;
  }
  else {
    const G4ElementVector* evec = sourcemat->GetElementVector();
    // Find plutonium and isotopes
    G4double mass = 0;
    G4double pu238 = 0;
    G4double pu240 = 0;
    G4double pu242 = 0;
    for(size_t e=0; e < evec->size(); e++) {
      if((*evec)[e]->GetZ() == 94) {
	G4IsotopeVector* ivec;
	ivec = (*evec)[e]->GetIsotopeVector();
	G4double* relabvec;
	relabvec = (*evec)[e]->GetRelativeAbundanceVector();
	G4double totw = 0;
	G4double * wtvec = new G4double[ivec->size()];
	for(size_t i=0; i < ivec->size(); i++) {
	  G4double w;
	  w = relabvec[i] * (*ivec)[i]->GetA();
	  wtvec[i] = w;
	  totw += w;
	}
	for(size_t i=0; i < ivec->size(); i++) {
	  wtvec[i] /= totw;
	  if((*ivec)[i]->GetN() == 238) {
	    pu238 = wtvec[i]; 
	  }
	  else if ((*ivec)[i]->GetN() == 240) {
	    pu240 = wtvec[i]; 
	  }
	  else if ((*ivec)[i]->GetN() == 242) {
	    pu242 = wtvec[i]; 
	  }
	}
	return nmsmr.pumass(pu238, pu240, pu242);
      }

    }
    return 0.0;
  }
}

void NMSAnalysisManager::calculateResults() {
  if(runtimeFromPulsetrain) {
    nmspm->setMeasurementLengthFromPulsetrain();
  }
  else {
    G4double runtime = static_cast<NMSRunManager*>(NMSRunManager::GetRunManager())->GetRuntime();
    if(runtime == 0) {
      runtime = 1 * s;
    }
    nmspm->setMeasurementLength(runtime / microsecond);
  }
  
  // Settings
  nmsmr = nmspm->getResults(pulsetrainAnalysisMode);
  nmsmr.normalize();

  calculated = true;
}

void NMSAnalysisManager::showResults() {
  if(!calculated) {
    calculateResults();
  }
  G4double runtime = 0;
  if(runtimeFromPulsetrain) {
    runtime = nmsmr.getLastEvent() * picosecond;
  }
  else {
    runtime = static_cast<NMSRunManager*>(NMSRunManager::GetRunManager())->GetRuntime();
  }
  if(runtime == 0) {
    runtime = 1 * s;
  }

  G4cout << "********************************************************************************" << G4endl;
  G4cout << "** Runtime: " << runtime / s << " s" << G4endl;
  G4cout << "** Last Event: " << nmsmr.getLastEvent() * picosecond / s << G4endl;
  //  G4cout << "** End Run:    " << id << G4endl;
  //  G4cout << "** No of Events: " << events << G4endl;
  //  G4cout << "** Number of spontaneous fissions: " << sfcount << G4endl;
  G4cout << "** Primary alpha particles started: " << sourceAlphaCount << G4endl;
  G4cout << "** Primary gamma particles started: " << sourceGammaCount << G4endl;
  //  G4cout << "** Primary beta particles started: " << sourceAlphaCount << G4endl;
  G4cout << "**" << G4endl;
  G4cout << "** Primary neutron particles started: " << sourceNeutronCount << G4endl;

  //  G4cout << "** SF branching ratio: " << 1. * sfcount / events << G4endl;
  //  G4cout << "** Alpha decay branching ratio: " << 1. * alphacount / events << G4endl;

  G4cout << "** Average neutrons per SF: " << getAverageNu() << G4endl;
  G4cout << "** Average neutrons per time: " << getAverageSourceNeutronsPerTime() << " 1/s" << G4endl;
  G4cout << "** Secondary neutrons started: " << secondaryNeutronCount << G4endl;

  //  G4cout << "** Average secondary neutrons per SF: " << 1. * secondaryneutroncount / sfcount << G4endl;

  G4cout << "** Total neutrons started: " << sourceNeutronCount + secondaryNeutronCount << G4endl;
  //  G4cout << "** Average total neutrons per SF: " << 1. * totalneutroncount / sfcount << G4endl;
  G4cout << "**" << G4endl;
  G4cout << "** Average source neutron energy: " << getSourceNeutronsEnergy() / MeV << " MeV" << G4endl;
  G4cout << "** Average secondary neutron energy: " << getSecondaryNeutronsEnergy() / MeV << " MeV" << G4endl;
  //  G4cout << "** Average total neutron energy: " << (totalneutronenergy + secondarytotalneutronenergy) / totalneutroncount / MeV << " MeV" << G4endl;
  G4cout << "**" << G4endl;
  G4cout << "** Number of neutrons absorbed in Detector (" << detectorVolume << "): " << nmspm->getEventCount() << G4endl;
  G4cout << "** Number of neutrons absorbed in 'Lost Volume' (" << lostNeutronVolume << "): " << getLostNeutrons() << G4endl;
  G4cout << "** Simulated detector efficiency (Source neutrons): " << getDetectorNeutronEfficiency(EFFICIENCY_SOURCE) << G4endl;
  G4cout << "** Simulated detector efficiency (Total neutrons): " << getDetectorNeutronEfficiency(EFFICIENCY_TOTAL) << G4endl;
  G4cout << "** Average neutron lifetime: " << nmspm->getAverageLifetime() << " microseconds" << G4endl;


  std::cout << nmsmr << std::endl;
  G4cout << "== Singles:        " << std::setw(15) << nmsmr.singles(1) << std::setw(15) << nmsmr.singles(1 / runtime * s) << " 1/s" << G4endl;
  G4cout << "== Doubles:        " << std::setw(15) << nmsmr.doubles(1) << std::setw(15) << nmsmr.doubles(1 / runtime * s) << " 1/s" << G4endl;
  G4cout << "== Triples:        " << std::setw(15) << nmsmr.triples(1) << std::setw(15) << nmsmr.triples(1 / runtime * s) << " 1/s" << G4endl;
  G4cout << "== For the pulsetrain split into " << splits << " equal parts:" << std::endl;
  G4cout << "== Avg Singles:    " << std::setw(15) << nmsmr.avgsingles(1) << std::setw(15) << nmsmr.avgsingles(1 / (runtime / splits) * s) << " 1/s" << G4endl;
  G4cout << "== Avg Doubles:    " << std::setw(15) << nmsmr.avgdoubles(1) << std::setw(15) << nmsmr.avgdoubles(1 / (runtime / splits) * s) << " 1/s" << G4endl;
  G4cout << "== Avg Triples:    " << std::setw(15) << nmsmr.avgtriples(1) << std::setw(15) << nmsmr.avgtriples(1 / (runtime / splits ) * s) << " 1/s" << G4endl;
  G4cout << "== StdDev Singles: " << std::setw(15) << nmsmr.uncsingles(1) << std::setw(15) << nmsmr.uncsingles(1 / (runtime / splits) * s) << " 1/s" << G4endl;
  G4cout << "== StdDev Doubles: " << std::setw(15) << nmsmr.uncdoubles(1) << std::setw(15) << nmsmr.uncdoubles(1 / (runtime / splits) * s) << " 1/s" << G4endl;
  G4cout << "== StdDev Triples: " << std::setw(15) << nmsmr.unctriples(1) << std::setw(15) << nmsmr.unctriples(1 / (runtime / splits ) * s) << " 1/s" << G4endl;

  G4cout << "**" << G4endl;
  G4cout << "** Settings for mass calculation" << G4endl;
  G4cout << "** efficiency=" << nmsmr.GetSettings().efficiency
	 << ", dieaway=" << nmsmr.GetSettings().dieaway
	 << ", doublegatefraction=" << nmsmr.GetSettings().doublegatefraction
	 << ", triplegatefraction=" << nmsmr.GetSettings().triplegatefraction
	 << G4endl;
  switch(nmsmr.GetSettings().dm) {
  case DAM_FIXED:
    G4cout << "** Fixed dieaway time" << G4endl;
    break;
  case DAM_AVERAGELIFETIME:
    G4cout << "** Dieaway time calculated from average neutron lifetime" << G4endl;
    break;
  }
  switch(nmsmr.GetSettings().gm) {
  case GM_FROM_DIEAWAY:
    G4cout << "** Gate fractions calculated from dieaway time (single exponential)" << G4endl; break;
  case GM_FROM_LONGGATE:
    G4cout << "** Gate fractions calculated by comparing gate to 'infinite' gate (length " << nmsmr.GetSettings().infinitegate << " microseconds" << G4endl; break;
  case GM_FROM_EVENTNO:
    G4cout << "** Gate fractions calculated per single event" << G4endl; break;
  case GM_FIXED:
    G4cout << "** Fixed gate fractions" << G4endl; break;
  }
  G4cout << "== Pu240eff-Mass : " << nmsmr.pu240eff() << G4endl;
  G4cout << "== Pu-Mass, calculated from source material: " << getTotalPuMass() << G4endl;

  std::streamsize ss = G4cout.precision();
  G4cout << "** Alpha: " << std::setprecision(12) << nmsmr.alpha() << std::setprecision(ss) << G4endl;
  G4cout << "** Multiplication: " << std::setprecision(12) << nmsmr.mul() << std::setprecision(ss) << G4endl;
  G4cout << "** Fission rate: " << std::setprecision(12) << nmsmr.fissionRate() << std::setprecision(ss) << G4endl;
}

void NMSAnalysisManager::showDetectorResults() {

}


void NMSAnalysisManager::exportToFile(G4String filename) {
  if(!calculated) {
    calculateResults();
  }
  std::ofstream fs;
  G4double runtime = 0;
  if(runtimeFromPulsetrain) {
    runtime = nmsmr.getLastEvent() * picosecond;
  }
  else {
    runtime = static_cast<NMSRunManager*>(NMSRunManager::GetRunManager())->GetRuntime();
  }
  if(runtime == 0) {
    runtime = 1 * s;
  }

  G4String  resultsfilename;
  resultsfilename = filename + ".results";
  fs.open(resultsfilename.c_str());
  fs << "runtime," << runtime / s << ",s" << std::endl;
  fs << "lastevent," << nmsmr.getLastEvent() * picosecond / s << ",s" << std::endl;
  fs << "primaryalpha," << sourceAlphaCount << ",1" << std::endl;
  //  fs << "primarybeta" << std::endl;
  fs << "primarygamma," << sourceGammaCount << ",1" << std::endl;
  fs << "primaryneutron," << sourceNeutronCount << ",1" << std::endl;
  fs << "secondaryneutron," << secondaryNeutronCount << ",1" << std::endl;
  fs << "totalneutron," << sourceNeutronCount + secondaryNeutronCount << ",1" << std::endl;
  fs << "sfnubar," << getAverageNu() << ",n/SF" << std::endl;
  fs << "neutronpertime," << getAverageSourceNeutronsPerTime() << ",1/s" << std::endl;
  fs << "primaryneutronenergy," << getSourceNeutronsEnergy() / MeV << ",MeV" << std::endl;
  fs << "secondaryneutronenergy," << getSecondaryNeutronsEnergy() / MeV << ",MeV" << std::endl;
  fs << "neutronabsorbeddetector," << nmspm->getEventCount() << ",1" << std::endl;
  fs << "detectorvolname," << detectorVolume << ",-" << std::endl;
  fs << "neutronabsorbedlostvolume," << getLostNeutrons() << ",1" << std::endl;
  fs << "lostvolumename," << lostNeutronVolume << ",-" << std::endl;
  fs << "efficiencysource," << getDetectorNeutronEfficiency(EFFICIENCY_SOURCE) << ",1" << std::endl;
  fs << "efficiencytotal," << getDetectorNeutronEfficiency(EFFICIENCY_TOTAL) << ",1" << std::endl;
  fs << "neutronlifetime," << nmspm->getAverageLifetime() << ",ms" << std::endl;
  fs << "singlestotal," << nmsmr.singles(1) << ",1" << std::endl;
  fs << "doublestotal," << nmsmr.doubles(1) << ",1" << std::endl;
  fs << "triplestotal," << nmsmr.triples(1) << ",1" << std::endl;
  fs << "singlespertime," << nmsmr.singles(1 / runtime * s) << ",1/s" << std::endl;
  fs << "doublespertime," << nmsmr.doubles(1 / runtime * s) << ",1/s" << std::endl;
  fs << "triplespertime," << nmsmr.triples(1 / runtime * s) << ",1/s" << std::endl;
  fs << "avgsinglespertime," << nmsmr.avgsingles(1 / (runtime / splits) * s) << ",1/s" << std::endl;
  fs << "avgdoublespertime," << nmsmr.avgdoubles(1 / (runtime / splits) * s) << ",1/s" << std::endl;
  fs << "avgtriplespertime," << nmsmr.avgtriples(1 / (runtime / splits) * s) << ",1/s" << std::endl;
  fs << "uncsinglespertime," << nmsmr.uncsingles(1 / (runtime / splits) * s) << ",1/s" << std::endl;
  fs << "uncdoublespertime," << nmsmr.uncdoubles(1 / (runtime / splits) * s) << ",1/s" << std::endl;
  fs << "unctriplespertime," << nmsmr.unctriples(1 / (runtime / splits) * s) << ",1/s" << std::endl;
  std::streamsize ss = G4cout.precision();
  fs << std::setprecision(12) << "multiplication," << nmsmr.mul() << ",1" << std::endl;
  fs << std::setprecision(12) << "alpha," << nmsmr.alpha() << ",1" << std::endl;
  fs << std::setprecision(12) << "fissionrate," << nmsmr.fissionRate() << ",1/s" << std::endl;
  fs << std::setprecision(ss);
  fs << "pu240effmass," << nmsmr.pu240eff() << ",g" << std::endl;
  fs << "pumass," << getTotalPuMass() << ",g" << std::endl; // FIX
  fs.close();

  G4String  tablefilename;
  tablefilename = filename + ".table";
  fs.open(tablefilename.c_str());
  fs << std::setw(13) << "Multiplicity" << std::setw(16) << "R + A" << std::setw(16) << "A" << std::setw(16) << "R + A (norm)" << std::setw(16) << "A (norm)" << std::endl;;
  for(int i = 0; i < nmsmr.getMultiplicities(); i++) {
    fs << std::setw(13) << i
	  << std::setw(16) << nmsmr.getRA(i)
	  << std::setw(16) << nmsmr.getA(i)
	  << std::setw(16) << nmsmr.getNormRA(i)
	  << std::setw(16) << nmsmr.getNormA(i) << std::endl;
  }
  fs << std::setw(13) << "i" << std::setw(16) << "R + A moment" << std::setw(16) << "A moment" << std::endl;
  for(int i = 0; i < 4; i++) {
    fs << std::setw(13) << i
       << std::setw(16) << nmsmr.getNuRA(i)
       << std::setw(16) << nmsmr.getNuA(i)
       << std::endl;
  }
  fs.close();
  
  if(exportPulsetrain) {
      nmspm->savePulsetrainToFile(filename);
  }
  nmspm->saveSettingsToFile(filename);
}

void NMSAnalysisManager::exportToFile() {
  exportToFile(outputfilename);
}

void NMSAnalysisManager::exportCheck(G4int runno) {
  if(exportAfterRun) {
    std::stringstream ss;
    G4String filename = "";
    ss << outputfilename << "." << runno;
    ss >> filename;
    ss.clear();
    exportToFile(filename);
  }
}

void NMSAnalysisManager::loadPulsetrainFromFile(G4String filename) {
  calculated = false;
  nmspm->loadPulsetrainFromFile(filename);
}


void NMSAnalysisManager::loadPulsetrainFromShakesFile(G4String filename) {
  calculated = false;
  nmspm->loadPulsetrainFromShakesFile(filename);
}

void NMSAnalysisManager::addPulsetrainFromFile(G4String filename, G4double offset) {
  calculated = false;
  nmspm->addPulsetrainFromFile(filename, offset);
}

void NMSAnalysisManager::splitPulsetrain(G4String filename, G4int splits) {
  calculated = false;
  nmspm->splitPulsetrain(filename, splits);
}

void NMSAnalysisManager::SetDetectorVolume(G4String newvolume) {
  detectorVolume = newvolume;
}

void NMSAnalysisManager::SetLostNeutronVolume(G4String newvolume) {
  lostNeutronVolume = newvolume;
}

void NMSAnalysisManager::SetRunMode(RunMode newrm) {
  rm = newrm;
  NMSTrackingAction* nmsta = (NMSTrackingAction*) NMSRunManager::GetRunManager()->GetUserTrackingAction();
  nmsta->SetRunMode(newrm);
}

void NMSAnalysisManager::SetPulsetrainAnalysisMode(G4int pam) {
  pulsetrainAnalysisMode = pam;
  calculated = false;
}

void NMSAnalysisManager::SetRuntimeFromPulsetrain(bool pfp) {
  runtimeFromPulsetrain = pfp;
  calculated = false;
}

void NMSAnalysisManager::SetOutputFilename(G4String filename) {
  if(filename != "") {
    outputfilename = filename;
  }
  else {
    G4cout << "Warning! NMSAnalysisManager::SetOutputFilename - given filename was empty. Set to default value" << G4endl;
    outputfilename = "nms-output";
  }
}
