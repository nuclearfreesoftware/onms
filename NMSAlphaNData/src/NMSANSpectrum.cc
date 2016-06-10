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

#include "NMSANSpectrum.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ParticleDefinition.hh"
#include "G4DecayProducts.hh"
#include "G4IonTable.hh"
#include "G4RadioactiveDecay.hh"
#include "G4AlphaDecay.hh"
#include "G4NucleiProperties.hh"

#include "NMSANcsdata.hh"

NMSANSpectrum::NMSANSpectrum() {
  neutronbins = 1000;
  neutronEmin = 0;
  neutronEmax = 10 * MeV;
  ngridbinwidth = (neutronEmax - neutronEmin) / neutronbins;

  mt91 = true;
  
  verboseLevel = 0;
  stoppingData = new NMSStoppingPower();
  stoppingData->setDataSource(STOPPINGPOWER_DS_ZIEGLER_T3);

  csData = NMSANcsdata::Instance();
  csData->Init();
}

NMSANSpectrum::~NMSANSpectrum() {
  delete stoppingData;
}

G4LPhysicsFreeVector * NMSANSpectrum::fromSource(G4Material * mat) {
  clear(); // clear spectrumvector
  
  G4ParticleDefinition* ion = 0;
  G4double lifetime;
  G4int a;
  G4int z;
  G4int elements = mat->GetNumberOfElements();
  G4DecayProducts* products;
  G4double energy;
  G4double elementatoms;
  G4double * relabvec;

  for(G4int i = 0; i < elements; i++) {
    const G4Element * el = mat->GetElement(i);
    elementatoms = mat->GetVecNbOfAtomsPerVolume()[i];
    relabvec = el->GetRelativeAbundanceVector();

    for(size_t j = 0; j < el->GetNumberOfIsotopes(); j++) {
      const G4Isotope * iso = el->GetIsotope(j);
      a = iso->GetN();
      z = iso->GetZ();
      ion = G4IonTable::GetIonTable()->GetIon(z,a,0);
      lifetime = ion->GetPDGLifeTime();

      if(lifetime != -1) {
	G4double isoact = 0;
	isoact = elementatoms * relabvec[j] / lifetime;
	//	G4cout << "Found r-active, activity: " << isoact * cm3 * s << G4endl;
	//	G4cout << "other way: " << elementatoms * relabvec[j] / lifetime * cm3 * s << G4endl;
	//G4double chact = 0;
	G4RadioactiveDecay* decay = new G4RadioactiveDecay();
	G4DecayTable* dtable = decay->GetDecayTable(ion);
	G4VDecayChannel* dc;

	for(G4int cidx = 0; cidx < dtable->entries(); cidx++) {
	  dc = dtable->GetDecayChannel(cidx);
	  //	  G4cout << "BR: " << dc->GetBR() << G4endl;

	  // newer geant version have G4AlphaDecay as child class of G4NuclearDecayChannel
	  G4AlphaDecay* testalpha = dynamic_cast< G4AlphaDecay* >( dc );
	  if( testalpha != 0) {
	  // for backwards compatibility should there be...
	      //	      or  G4AlphaDecayChannel* testalpha = dynamic_cast< G4AlphaDecayChannel* >( dc ) ) {
	    products = dc->DecayIt(0);
	    for(G4int k = 0; k < products->entries(); k++) {
	      if((*products)[k]->GetParticleDefinition()->GetParticleName() == "alpha") {
		energy = (*products)[k]->GetKineticEnergy();
	      }
	    }
	    spectrumdata sd;
	    sd.intensity = isoact * dc->GetBR();
	    //G4cout << sd.intensity << "    E = " << energy << G4endl;
	    sd.spectrum = fromEnergy(energy, mat);
	    sd.energy = energy;
	    spectrumvector.push_back(sd);
	  }
	}
      }
    }
  }

  // sum all spectrum
  G4LPhysicsFreeVector * totalspectrum;
  totalspectrum = new G4LPhysicsFreeVector(neutronbins, neutronEmin, neutronEmax);
  
  for(G4int eidx = 0; eidx < neutronbins; eidx++) {
    G4double val = 0;
    for(size_t sidx = 0; sidx < spectrumvector.size(); sidx++) {
      val += spectrumvector[sidx].spectrum->operator[](eidx) * spectrumvector[sidx].intensity;
    }
    totalspectrum->PutValues(eidx, neutronEmin + eidx * (neutronEmax - neutronEmin) / neutronbins, val);
  }

  return totalspectrum;
}

G4LPhysicsFreeVector * NMSANSpectrum::fromEnergy(G4double E, G4Material * mat) {
  G4LPhysicsFreeVector * result = new G4LPhysicsFreeVector(neutronbins, neutronEmin, neutronEmax);
  for(G4int i = 0; i < neutronbins; i++) {
    result->PutValues(i, neutronEmin + i * (neutronEmax - neutronEmin) / neutronbins, 0);
  }
  stoppingData->initialize();

  G4int elements = mat->GetNumberOfElements();
  for(G4int i = 0; i < elements; i++) {
    const G4Element * el = mat->GetElement(i);
    G4double elementatoms = mat->GetVecNbOfAtomsPerVolume()[i];
    G4double* relabvec = el->GetRelativeAbundanceVector();
    for(size_t j = 0; j < el->GetNumberOfIsotopes(); j++) {
      G4double atden = relabvec[j] * elementatoms;
      const G4Isotope * iso = el->GetIsotope(j);
      G4LPhysicsFreeVector * tmp = fromEnergy(E, mat, iso);
      for(G4int k = 0; k < neutronbins; k++) {
	result->PutValue(k, result->operator[](k) + atden * tmp->operator[](k));
      }
      delete tmp;
    }
  }

  return result;
}

G4LPhysicsFreeVector * NMSANSpectrum::fromEnergy(G4double E, G4Material * mat, const G4Isotope * iso) {

  G4LPhysicsFreeVector * isoresult = new G4LPhysicsFreeVector(neutronbins, neutronEmin, neutronEmax);
  for(G4int i = 0; i < neutronbins; i++) {
    isoresult->PutValues(i, neutronEmin + i * (neutronEmax - neutronEmin) / neutronbins, 0);
  }

  //G4cout << "mt " << mt / MeV << "  mr " << mr / MeV << G4endl;
  //mt = G4IonTable::GetIonTable()->GetIon(iso->GetZ(), iso->GetN())->GetAtomicMass();
  //mr = G4IonTable::GetIonTable()->GetIon(iso->GetZ() + 3, iso->GetN() + 2)->GetAtomicMass();
  //G4cout << "mt " << mt / MeV << "  mr " << mr / MeV << G4endl;

  G4int csidx = csData->index(iso->GetZ(), iso->GetN());
  if(csidx != - 1 ) {
    G4int levels = csData->inelasticNo(csidx);
    //G4cout << levels << " levels" << G4endl;
    if(levels > 0) {
      G4LPhysicsFreeVector * tmp;
      for(G4int lidx = 0; lidx < levels; lidx++) {
	tmp = fromEnergy(E, mat, iso, lidx);
	for(G4int k = 0; k < neutronbins; k++) {
	  isoresult->PutValue(k, isoresult->operator[](k) + tmp->operator[](k));
	}
      }
      delete tmp;
    }
  }
  return isoresult;  
}

G4LPhysicsFreeVector * NMSANSpectrum::fromEnergy(G4double E, G4Material * mat, const G4Isotope * iso, G4int levelidx) {
  G4LPhysicsFreeVector * levelresult = new G4LPhysicsFreeVector(neutronbins, neutronEmin, neutronEmax);
  for(G4int i = 0; i < neutronbins; i++) {
    levelresult->PutValues(i, neutronEmin + i * (neutronEmax - neutronEmin) / neutronbins, 0);
  }

  G4double mn = neutron_mass_c2;
  G4double ma = G4Alpha::Alpha()->GetPDGMass();
  G4double tm = G4NucleiProperties::GetNuclearMass(iso->GetN(), iso->GetZ());
  G4double rm = G4NucleiProperties::GetNuclearMass(iso->GetN() + 3, iso->GetZ() + 2);

  G4int csidx = csData->index(iso->GetZ(), iso->GetN());

  G4int mt = csData->getInelasticCsMT(csidx, levelidx);
  if((mt >= 50 and mt <= 90) or (mt91 and mt == 91)) {
    if(verboseLevel >= 1) {
      G4cout << "Read MT=" << mt << " data" << G4endl;
    }
    G4LPhysicsFreeVector * cs = csData->getInelasticCsVector(csidx, levelidx);

    G4int idx = 0;
    G4double Qm = csData->getInelasticCsQex(csidx, levelidx); // Q-Value - Excitation energy
    G4double stepcs = 0;
    G4double deltae = 0;
    G4double yieldone = 0;
    G4double Enmin = 0; G4double Enmax = 0;
    G4double emin = cs->Energy(idx); G4double emax = cs->Energy(idx + 1);
    G4double csmin = cs->operator[](idx); G4double csmax = cs->operator[](idx + 1);
    G4double smin = 0;
    G4double smax = stoppingData->getDEDX(emax, mat);

    G4int tc2 = 0;
    while(emax <= E) {
      //		G4cout << emax / MeV << G4endl;
      //		G4cout << smax << G4endl;
      G4double stepstopping = (smin + smax) / 2.0;
      if(stepstopping != 0.0) {
	stepcs = (csmin + csmax) / 2;
	deltae = emax - emin;
	yieldone = stepcs / stepstopping * deltae;
      }
      else {
	yieldone = 0;
      }
      G4double sum1 = sqrt(emax * mn / ma) / ( 1 + tm / ma );
      G4double sum2a = (Qm / ( 1 + mn / rm));
      // G4cout << "tm / ma " << (tm / ma) << " mn / rm" << mn / rm << G4endl;
      G4double sum2b = (emax * (tm / ma) / (1 + tm / ma) / (1 + mn / rm ));
      // G4cout << sum2a << "  " << sum2b << G4endl;
      G4double sum2 = sqrt(sum2a + sum2b);
      Enmin = pow(( sum1 - sum2 ), 2);
      Enmax = pow(( sum1 + sum2 ), 2);
      // G4cout << sum1 << " " << sum2 << "  es: ";
      G4int tc = 0;
      for(G4int eidx = 0; eidx < neutronbins; eidx++) {
	//		  G4cout << levelresult->Energy(eidx) << G4endl;
	if(levelresult->Energy(eidx) >= Enmin && levelresult->Energy(eidx) <= Enmax) {
	  tc += 1;
	  G4double add = levelresult->operator[](eidx) + yieldone * ngridbinwidth / (Enmax - Enmin);
	  levelresult->PutValue(eidx, add);
	}
      }
      tc2 += 1;

      idx += 1;
      emin = cs->Energy(idx); emax = cs->Energy(idx + 1);
      csmin = cs->operator[](idx); csmax = cs->operator[](idx + 1);
      smin = smax;
      smax = stoppingData->getDEDX(emax, mat);
    }
  }

  return levelresult;
  
}


void NMSANSpectrum::clear() {
  for(size_t i; i < spectrumvector.size(); i++) {
    delete spectrumvector[i].spectrum;
  }
  spectrumvector.clear();
}
