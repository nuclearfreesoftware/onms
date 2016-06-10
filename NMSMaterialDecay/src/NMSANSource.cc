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

#include "NMSANSource.hh"
#include "NMSANSpectrum.hh"

NMSANSource::NMSANSource() : NMSNewSingleDecaySource("neutron from (alpha,n)") {
  verboseLevel = 0; // FIX

  nmsmdsSettings = NMSMaterialDecaySettings::Instance();
  nspectrum = 0;
}

NMSANSource::~NMSANSource() {
  
}

void NMSANSource::GeneratePrimaryVertex(G4Event* anEvent) {

  if(verboseLevel >= 2) {
    G4cout << "SDS-Derived-GPV (Alpha-N) called!" << G4endl;
  }
  
  if(!initialized) {
    Init();
  }

  G4ParticleDefinition* neutronDefinition = G4Neutron::NeutronDefinition();
  G4DynamicParticle* neutron;

  G4double sourceEnergy = 0;

  G4ThreeVector sourcePosition = SamplePosition();
  if(nmsmdsSettings->GetANPositionSamplingMode() == AN_POS_PRECALCMC &&  nmsmdsSettings->GetANEnergySamplingMode() == AN_E_PRECALCMC) {
    sourceEnergy = GetEnergyFromSampleIndex();    
  } else {
    sourceEnergy = SampleEnergy();
  }

  G4PrimaryVertex * vertex = new G4PrimaryVertex(sourcePosition, GetParticleTime());
  neutron = new G4DynamicParticle();
  neutron->SetDefinition(neutronDefinition);
  neutron->SetKineticEnergy(sourceEnergy);
  G4double mom = neutron->GetTotalMomentum();

  G4ThreeVector sourceDirection = SampleDirection();

  G4PrimaryParticle* particle = new G4PrimaryParticle(neutronDefinition,
						      mom * sourceDirection.getX(),
						      mom * sourceDirection.getY(),
						      mom * sourceDirection.getZ(),
						      sourceEnergy);

  vertex->SetPrimary(particle);
  anEvent->AddPrimaryVertex(vertex);

  if(verboseLevel >= 2) {
    G4cout << " Neutron from (alpha, n) reaction "  << G4endl;
    G4cout << "     Energy:   " << sourceEnergy << G4endl;
    G4cout << "     Momentum: "<<particle->GetMomentum() << G4endl;
    G4cout << "     Position: "<<vertex->GetPosition() << G4endl;
  }

}

void NMSANSource::Init() {
  NMSANEnergySamplingMode em = nmsmdsSettings->GetANEnergySamplingMode();

  if(em == AN_E_FIXED) {
    // nothing to do
  }
  if(em == AN_E_PRECALCMC) {
    
  }
  if(em == AN_E_CALC) {
    G4Material * mat = nmsmdsSettings->GetSourceMaterial();
    if(mat != 0) {
      G4cout << mat << G4endl;
    }
    else {
      G4cout << "No material defined" << G4endl;
    }
    if(nspectrum != 0) {
      if(nspectrum->GetVectorLength() != NMSANCalcNeutronBins) {
	nspectrum = new G4LPhysicsFreeVector(NMSANCalcNeutronBins, 0 * MeV, 10 * MeV);
      }
    }
    else {
      nspectrum = new G4LPhysicsFreeVector(NMSANCalcNeutronBins, 0 * MeV, 10 * MeV);
    }

    NMSANSpectrum * ans = new NMSANSpectrum();
    ans->SetUseMT91(nmsmdsSettings->GetANSpectrumMT91());

    G4double total = 0;
    G4LPhysicsFreeVector * spectrum = ans->fromSource(mat);
    for(G4int eidx = 0; eidx < NMSANCalcNeutronBins; eidx++) {
      nspectrum->PutValues(eidx, spectrum->Energy(eidx), total);
      total += spectrum->operator[](eidx);
    }
    if(total == 0) {
      G4Exception("NMSANSource::Init", "Spectrum was initialized with total=0", FatalException, "Somehow, neutron energy pre-calculation was unsuccessful - possibly wrong CS?");
    }
    nspectrum->ScaleVector(1, 1 / total);

    delete ans;
  }
  if(em == AN_E_SPECTRUMFILE) {
    G4String filename = nmsmdsSettings->GetANEnergyFilename();
    G4cout << "Reading Alpha N File: " << filename << G4endl;
    G4int countpoints;
    G4double energy = 0;
    G4double total = 0;
    G4double value = 0;
    if(nspectrum != 0) {
      delete nspectrum;
      nspectrum = 0;
    }
    std::ifstream filestream( filename , std::ios::in );
    if(filestream.good()) {
      filestream >> countpoints;
      G4cout << "countpoins" << countpoints;
      if(countpoints != 0) {
	nspectrum = new G4LPhysicsFreeVector(countpoints, 0 * MeV, 20 * MeV);
	for (G4int iDat=0;iDat<countpoints;iDat++) {
	  if(filestream.good()) {
	    filestream >> energy >> value;
	    nspectrum->PutValues(iDat, energy * MeV, total);
	    total += value;
	  }
	}
	nspectrum->ScaleVector(1, 1 / total);
	G4cout << "Reading successful" << G4endl;
      }
      else {
	G4cout << "Found file, seems to have wrong format!" << G4endl;
      }
    }
    else {
      G4cout << "Could not read from file!" << G4endl;
    }

  }
  initialized = true;
}

G4ThreeVector NMSANSource::SampleDirection() {
  NMSANDirectionSamplingMode dm = nmsmdsSettings->GetANDirectionSamplingMode();
  if(dm == AN_DIR_ISOTROPIC) {
      //Random, Isotropic direction
    G4double cosTheta = 2*G4UniformRand() - 1., phi = twopi*G4UniformRand();
    G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
    G4double ux = sinTheta*std::cos(phi), uy = sinTheta*std::sin(phi), uz = cosTheta;

    return G4ThreeVector(ux, uy, uz);
  }
  if(dm == AN_DIR_PRECALCMC) {

  }
  // default return value
  return G4ThreeVector(1, 0, 0);
}

G4double NMSANSource::SampleEnergy() {
  NMSANEnergySamplingMode em = nmsmdsSettings->GetANEnergySamplingMode();
  if(verboseLevel >= 2) {
      G4cout << "Sampling energy" << G4endl;
  }
  if(em == AN_E_FIXED) {
    return nmsmdsSettings->GetANEnergy();
  }
  if(em == AN_E_PRECALCMC) {
    
  }
  if(em == AN_E_CALC) {
    if(nspectrum != 0) {
      G4double ran = G4UniformRand();
      G4int i = 0;
      while(ran > nspectrum->operator[](i) && i < NMSANCalcNeutronBins) {
	i++;
      }
      G4double energy = nspectrum->Energy(i);
      //G4cout << "Source E-idx: " << i << ", energy=" <<  energy / MeV << "MeV" << G4endl;
      return energy;
    }
    else {
      G4Exception("NMSANSource::SampleEnergy", "Spectrum not initialized", FatalException, "Somehow, there was no neutron energy distribution pre-calculated.");
    }
  }
  if(em == AN_E_SPECTRUMFILE) {
    if(nspectrum != 0) {
      G4double ran = G4UniformRand();
      G4int i = 0;
      while(ran > nspectrum->operator[](i) && i < nspectrum->GetVectorLength()) {
	i++;
      }
      G4double energy = nspectrum->Energy(i);
      //G4cout << "Source E-idx: " << i << ", energy=" <<  energy / MeV << "MeV" << G4endl;
      return energy;
    }
    else {
      G4Exception("NMSANSource::SampleEnergy", "Spectrum not initialized", FatalException, "Somehow, there was no neutron energy distribution that could be read from a file.");
    }
  }
  
  return 1.0;
}

G4ThreeVector NMSANSource::SamplePosition() {
  NMSANPositionSamplingMode pm = nmsmdsSettings->GetANPositionSamplingMode();
  if(pm == AN_POS_SOURCEVOLUME) {
      return GetPosDist()->GenerateOne();
  }
  if(pm == AN_POS_PRECALCMC) {
    // FIX Sample from File
    return G4ThreeVector(0, 0, 0);
  }
  return G4ThreeVector(0, 0, 0);
}

G4double NMSANSource::GetEnergyFromSampleIndex() {
  return 1.0;
}

void NMSANSource::WriteANSpectrum(G4String filename) {
  if(!initialized) {
    Init();
  }

  if(nspectrum != 0) {
    G4String fnCDF = filename + ".alphancdf";
    std::ofstream outfile;
    outfile.open(fnCDF.c_str());
    outfile << nspectrum->GetVectorLength() << std::endl;
    for(size_t i = 0; i < nspectrum->GetVectorLength(); i++) {
      outfile  << nspectrum->Energy(i) << " " << nspectrum->operator[](i) << std::endl;
    }
    outfile.close();
    G4cout << "Wrote cumulative density function to " << fnCDF << G4endl;

    G4String fnDF = filename + ".alphan";
    outfile.open(fnDF.c_str());
    outfile << nspectrum->GetVectorLength() << std::endl;
    G4double total = 0;
    for(size_t i = 0; i < nspectrum->GetVectorLength(); i++) {
      outfile  << nspectrum->Energy(i) << " " << nspectrum->operator[](i) - total << std::endl;
      total = nspectrum->operator[](i);
    }
    outfile.close();
    G4cout << "Wrote spectrum to " << fnDF << G4endl;
  }
  else {
    G4cout << "Warning: Could not write anspectrum to file - none was loaded previously." << G4endl;
  }

}
