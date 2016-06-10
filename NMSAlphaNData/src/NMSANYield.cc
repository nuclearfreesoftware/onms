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

#include "NMSANYield.hh"

const G4int NMSANYield::csno = 17;
const G4int NMSANYield::csarray[csno][2] = {{3, 6}, {3, 7},
					    {4, 9},
					    {5, 10}, {5, 11},
					    {6, 12}, {6, 13},
					    {7, 14}, {7, 15},
					    {8, 17}, {8, 18},
					    {9, 19},
					    {11, 23},
					    {13, 27},
					    {14, 28}, {14, 29}, {14, 30}
};
const G4String NMSANYield::nameString[100] = {"Hydrogen", "Helium",
 "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen", "Fluorine",
 "Neon", "Sodium", "Magnesium", "Aluminum", "Silicon", "Phosphorous", 
 "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium", "Scandium",
 "Titanium", "Vanadium", "Chromium", "Manganese", "Iron", "Cobalt", "Nickel",
 "Copper", "Zinc", "Gallium", "Germanium", "Arsenic", "Selenium", "Bromine",
 "Krypton", "Rubidium", "Strontium", "Yttrium", "Zirconium", "Niobium",
 "Molybdenum", "Technetium", "Ruthenium", "Rhodium", "Palladium", "Silver",
 "Cadmium", "Indium", "Tin", "Antimony", "Tellurium", "Iodine", "Xenon",
 "Cesium", "Barium", "Lanthanum", "Cerium", "Praseodymium", "Neodymium",
 "Promethium", "Samarium", "Europium", "Gadolinium", "Terbium", "Dysprosium",
 "Holmium", "Erbium", "Thulium", "Ytterbium", "Lutetium", "Hafnium",
 "Tantalum", "Tungsten", "Rhenium", "Osmium", "Iridium", "Platinium", "Gold",
 "Mercury", "Thallium", "Lead", "Bismuth", "Polonium", "Astatine", "Radon", 
 "Francium", "Radium", "Actinium", "Thorium", "Protactinium", "Uranium", 
 "Neptunium", "Plutonium", "Americium", "Curium", "Berkelium", "Californium",
 "Einsteinium","Fermium"};


NMSANYield::NMSANYield() {
  verboseLevel = 0;
  stoppingData = new NMSStoppingPower();
  stoppingData->setDataSource(STOPPINGPOWER_DS_ZIEGLER_T3);

  csData = NMSANcsdata::Instance();
  csData->Init();
  //  initialize();
}

NMSANYield::~NMSANYield() {
  delete stoppingData;
}

G4double NMSANYield::fromSource(G4Material * mat) {

  G4ParticleDefinition* ion = 0;
  G4double lifetime;
  //G4double halflife;
  G4int a;
  G4int z;
  G4int elements = mat->GetNumberOfElements();
  G4DecayProducts* products;
  G4double energy;
  G4double elementatoms;
  G4double * relabvec;

  yieldvector.clear();
  
  for(G4int i = 0; i < elements; i++) {
    const G4Element * el = mat->GetElement(i);
    elementatoms = mat->GetVecNbOfAtomsPerVolume()[i];
    relabvec = el->GetRelativeAbundanceVector();

    for(size_t j = 0; j < el->GetNumberOfIsotopes(); j++) {
      //      G4cout << relabvec[j] << "  " << elementatoms << G4endl;
      const G4Isotope * iso = el->GetIsotope(j);
      a = iso->GetN();
      z = iso->GetZ();
      ion = G4IonTable::GetIonTable()->GetIon(z,a,0);
      lifetime = ion->GetPDGLifeTime();
      //      halflife = lifetime*log(2);


      if(lifetime != -1) {
	G4double isoact = 0;
	isoact = elementatoms * relabvec[j] / lifetime;
	//	G4cout << "Found r-active, activity: " << isoact * cm3 * s << G4endl;
	//	G4cout << "other way: " << elementatoms * relabvec[j] / lifetime * cm3 * s << G4endl;
	// G4double chact = 0;
	G4RadioactiveDecay* decay = new G4RadioactiveDecay();
	G4DecayTable* dtable = decay->GetDecayTable(ion);
	G4VDecayChannel* dc;

	// isoact = relabvec[j] * elementatoms / lifetime;

	//	G4cout << "Decay Table entries: " << dtable->entries() << G4endl;
	// G4double br = 0;
	for(G4int cidx = 0; cidx < dtable->entries(); cidx++) {
	  dc = dtable->GetDecayChannel(cidx);
	  //	  G4cout << "BR: " << dc->GetBR() << G4endl;

	  // newer geant version have G4AlphaDecay as child class of G4NuclearDecayChannel
	  if( G4AlphaDecay* testalpha = dynamic_cast< G4AlphaDecay* >( dc )) {
	  // for backwards compatibility could be the following ...
	      //	      or  G4AlphaDecayChannel* testalpha = dynamic_cast< G4AlphaDecayChannel* >( dc ) ) {
	    products = dc->DecayIt(0);
	    for(G4int k = 0; k < products->entries(); k++) {
	      if((*products)[k]->GetParticleDefinition()->GetParticleName() == "alpha") {
		energy = (*products)[k]->GetKineticEnergy();
	      }
	    }
	    yielddata ny;
	    ny.activity = isoact * dc->GetBR();
	    ny.energy = energy;
	    ny.yieldfore = fromEnergy(energy, mat);
	    ny.yield = fromEnergy(energy, mat) * ny.activity;
	    yieldvector.push_back(ny);
	  }
	}
	// somewhere should be: delete decay;

      }

    }
  }
  G4double fullyield = 0;
  for(size_t i = 0; i< yieldvector.size(); i++) {
    // G4cout << yieldvector[i].activity * cm3 * s<< " "
    // 	   << yieldvector[i].energy << " "
    // 	   << yieldvector[i].yieldfore << " "
    // 	   << yieldvector[i].yield * cm3 * s<< G4endl;
    fullyield += yieldvector[i].yield;
  }

  return fullyield;
}

G4double NMSANYield::alphaActivityFromSource(G4Material * mat) {
  G4ParticleDefinition* ion = 0;

  G4int elements = mat->GetNumberOfElements();
  G4double elementatoms = 0;
  G4double * relabvec;
  G4int a;
  G4int z;
  G4double lifetime = 0;
  //G4double halflife = 0;

  G4double totact = 0;
  for(G4int i = 0; i < elements; i++) {
    const G4Element * el = mat->GetElement(i);
    elementatoms = mat->GetVecNbOfAtomsPerVolume()[i];
    relabvec = el->GetRelativeAbundanceVector();
    for(size_t j = 0; j < el->GetNumberOfIsotopes(); j++) {
      //      G4cout << relabvec[j] << "  " << elementatoms << G4endl;
      const G4Isotope * iso = el->GetIsotope(j);
      a = iso->GetN();
      z = iso->GetZ();
      ion = G4IonTable::GetIonTable()->GetIon(z,a,0);
      lifetime = ion->GetPDGLifeTime();
      //halflife = lifetime*log(2);

      if(lifetime != -1) {
	G4double isoact = 0;
	isoact = elementatoms * relabvec[j] / lifetime;
	totact += isoact;
      }
    }
  }
  return totact;
}

G4double NMSANYield::fromEnergy(G4double E, G4Material * mat) {

  stoppingData->initialize();

  G4double totalyield = 0;
  
  G4int elements = mat->GetNumberOfElements();
  // G4double totab = mat->GetTotNbOfAtomsPerVolume();
  //  G4cout << "material: number of atoms " << totab * cm3 << G4endl; 
  for(G4int i = 0; i < elements; i++) {
    const G4Element * el = mat->GetElement(i);
    //    G4cout << el->GetName() << G4endl;

    G4double elementatoms = mat->GetVecNbOfAtomsPerVolume()[i];
    //    G4cout << "Element atoms: " << elementatoms * cm3 << G4endl;
    G4double* relabvec = el->GetRelativeAbundanceVector();
    for(size_t j = 0; j < el->GetNumberOfIsotopes(); j++) {
      G4double atden = relabvec[j] * elementatoms;
      const G4Isotope * iso = el->GetIsotope(j);
      G4int csidx = csData->index(iso->GetZ(), iso->GetN());

      // //      G4cout << iso->GetName() << " Z: " << iso->GetZ() << ", A: " << iso->GetN() << " Den:" << atden * cm3 << G4endl;
      // for(G4int k = 0; k < csData->no(); k++) {
      // 	if(csarray[k][0] == iso->GetZ() && csarray[k][1] == iso->GetN()) {
      // 	  csidx = k;
      // 	}
      // }
      if(csidx != -1) {
	G4double isoyield = 0;
	G4int idx = 0;
	G4LPhysicsFreeVector * pv = csData->getcsVector(csidx);
	G4double lowe = pv->Energy(idx);
	G4double highe = pv->Energy(idx + 1);
	G4double lowecs = (*pv)[idx];
	G4double highecs = (*pv)[idx + 1];
	G4double highestopping = 0;
	
	G4double lowestopping = 0;
	G4double stepfull = 0;
	G4double stoppingmean = 0;
	while(highe < E) {
	  lowestopping = highestopping;
	  highestopping = stoppingData->getDEDX(highe, mat);
	  stoppingmean = (lowestopping + highestopping) / 2;

	  if(stoppingmean != 0.0) {
	    G4double stepcsint = ((lowecs + highecs) / 2 * (highe - lowe));
	    stepfull = stepcsint / (stoppingmean);
	    //	    G4cout << lowe / MeV << " " << stoppingmean / (MeV * cm * cm) * g << " " << stepcsint / barn / MeV << G4endl;
	    //	    stepfull = stepfull / 1e24;
	  }
	  else {
	    stepfull = 0;
	  }

	  isoyield += stepfull;

	  //	  G4cout << idx << G4endl;
	  idx += 1;
	  lowe = pv->Energy(idx);
	  lowecs = (*pv)[idx];
	  highe = pv->Energy(idx + 1);
	  highecs = (*pv)[idx + 1];
	}
	isoyield = isoyield * atden;
	totalyield += isoyield;
      }

    }
  }

  return totalyield;
}



void NMSANYield::initialize() {
  // Initialize Stopping power data
  //stoppingData->initialize();

  // Load cs data
  if(!getenv("NMSALPHALEDATA")) 
    G4Exception("NMSANYield::initialize", "No CS found", FatalException, "Please setenv NMSALPHALEDATA to point to the alpha n cross-section files.");
  G4String dirName = getenv("NMSALPHALEDATA");
  G4String filename = "";

  for(  std::vector<G4LPhysicsFreeVector *>::iterator it = csIsoVector.begin() ; it != csIsoVector.end(); ++it) {
    delete(*it);
  }
  csIsoVector.clear();
  
  for(G4int i = 0; i < csno; i++) {
    G4int Z = csarray[i][0];
    G4int A = csarray[i][1];
    std::stringstream ss;
    ss << dirName << "/Inelastic/CrossSection/" << Z << "_" << A << "_" << nameString[Z - 1];
    ss >> filename;
    ss.clear();
    //	filename = dirName + "/Inelastic/CrossSection/" + std::to_string(Z) + "_" + std::to_string(A) + (*G4Element::GetElementTable())[elI]->GetName();
    if(verboseLevel >= 1) {
      G4cout << "NMSAlphaLE: Try to read from " << filename << G4endl;
    }
    std::ifstream filestream( filename , std::ios::in );
    if(filestream.good()) {
      G4String dummy;
      filestream >> dummy; // Set
      filestream >> dummy; // Source
      G4int nData;
      filestream >> nData;
      if(verboseLevel >= 2) {
	G4cout << "NMSAlphaLE: File has " << nData << " entries" << G4endl;
      }
      G4LPhysicsFreeVector * temp = new G4LPhysicsFreeVector(nData, 0 * MeV, 15 * MeV);
      G4double x,y;
      for (G4int iDat=0;iDat<nData;iDat++) {
	if(filestream.good()) {
	  filestream >> x >> y;
	  temp->PutValues(iDat, x * eV, y * barn);
	}
      }
      csIsoVector.push_back(temp);
      // read into pv + add to data
      filestream.close();
    }
    else {
      if(verboseLevel >= 1) {
	G4cout << "NMSAlphaLE: Could not find file with data, using GG-cs" << G4endl;
      }
    }

  }
  
}


