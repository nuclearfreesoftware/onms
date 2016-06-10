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


#include "NMSMaterialDecaySourceMessenger.hh"

NMSMaterialDecaySourceMessenger::NMSMaterialDecaySourceMessenger(NMSMaterialDecaySource * newnmsmds) {
  nmsmds = newnmsmds;
  nmsmdsSettings = NMSMaterialDecaySettings::Instance();

  // Directories
  sourceDir = new G4UIdirectory("/NMS/source/");
  sourceDir->SetGuidance("NMSMaterialDecaySource parameter set for source definition");

  sourcePosDistDir = new G4UIdirectory("/NMS/source/geometry/");
  sourcePosDistDir->SetGuidance("NMSMaterialDecaySource source geometry setup");

  sourceANDir = new G4UIdirectory("/NMS/source/alphan/");
  sourceANDir->SetGuidance("NMSMaterialDecaySource setup component for neutrons from (alpha, n)-reactions");

  // Verbose
  VerboseCmd = new G4UIcmdWithAnInteger("/NMS/source/verbose", this);
  VerboseCmd->SetGuidance("Set verbose level for MaterialDecaySource");
  VerboseCmd->SetDefaultValue(0);

  // Source Material, Volume and Activity
  MaterialCmd = new G4UIcmdWithAString("/NMS/source/material", this);
  MaterialCmd->SetGuidance("Select material for NMSMaterialDecaySource");
  MaterialCmd->SetParameterName("material", false, false);
  MaterialCmd->SetDefaultValue("NULL");

  SourceFromVolumeCmd = new G4UIcmdWithAString("/NMS/source/sourcefromvolume", this);
  SourceFromVolumeCmd->SetGuidance("Specify a physical volume name that should be used as source");
  
  Cf252Cmd = new G4UIcommand("/NMS/source/cf252sfoptions", this);
  Cf252Cmd->SetGuidance("Set options for multiplicity and energy distribution of Cf252 spontaneous fission.");
  Cf252Cmd->SetGuidance("(put through to fission library)");
  G4UIparameter * param;
  param = new G4UIparameter("ndist", 'i', true);
  param->SetParameterRange("ndist>=0 && ndist<=1");
  param->SetGuidance("  0 to sample the spontaneous fission neutron multiplicity using tabulated data from Spencer");
  param->SetGuidance("  1 to sample the spontaneous fission neutron multiplicity using tabulated data from Boldeman");
  param->SetDefaultValue(0);
  Cf252Cmd->SetParameter(param);
  param = new G4UIparameter("neng", 'i', true);
  param->SetParameterRange("neng>=0 && neng<=2");
  param->SetGuidance("  0 to sample the Mannhart corrected Maxwellian spectrum");
  param->SetGuidance("  1 to sample the Madland-Nix theoretical spectrum");
  param->SetGuidance("  2 to sample the Froehner Watt spectrum (a=1.175, b=1.04)");
  param->SetDefaultValue(0);
  Cf252Cmd->SetParameter(param);

  ActiveVolumeCmd = new G4UIcmdWithADoubleAndUnit("/NMS/source/activevolume", this);
  ActiveVolumeCmd->SetGuidance("Set active volume of source (density as in material specification)");
  ActiveVolumeCmd->SetParameterName("activevolume", false, false);
  ActiveVolumeCmd->SetDefaultValue(1);
  ActiveVolumeCmd->SetDefaultUnit("cm3");

  ActivityCmd = new G4UIcmdWithADoubleAndUnit("/NMS/source/activity", this);
  ActivityCmd->SetGuidance("Set activity of source");
  ActivityCmd->SetParameterName("activity", false, false);
  ActivityCmd->SetDefaultValue(1);
  ActivityCmd->SetDefaultUnit("Hz");

  ActivityFixedCmd = new G4UIcmdWithABool("/NMS/source/activityfixed", this);
  ActivityFixedCmd->SetGuidance("Set to true to override automatic activity calculation.");
  ActivityFixedCmd->SetDefaultValue(false);
  //timing issues

  // Turn on/off Sourcetypes
  SFnCmd = new G4UIcmdWithABool("/NMS/source/sfn", this);
  SFnCmd->SetGuidance("Enable/Disable spontaneous fission (neutron emission) as a decay mode for source");
  SFnCmd->SetParameterName("sfn", false);

  SFgCmd = new G4UIcmdWithABool("/NMS/source/sfg", this);
  SFgCmd->SetGuidance("Enable/Disable spontaneous fission (gamma emission) as a decay mode for source");
  SFgCmd->SetParameterName("sfg", false);

  BetaCmd = new G4UIcmdWithABool("/NMS/source/beta", this);
  BetaCmd->SetGuidance("Enable/Disable beta decay for source");
  BetaCmd->SetParameterName("beta", false);

  AlphaCmd = new G4UIcmdWithABool("/NMS/source/alpha", this);
  AlphaCmd->SetGuidance("Enable/Disable alpha decay for source");
  AlphaCmd->SetParameterName("alpha", false);

  AlphaNCmd = new G4UIcmdWithABool("/NMS/source/neutronAlphaN", this);
  AlphaNCmd->SetGuidance("Enable/Disable neutrons from (alpha,n) reactions for source");
  AlphaNCmd->SetGuidance("  If enabled, more settings are needed in /NMS/source/alphan/!");
  AlphaNCmd->SetParameterName("alphaN", false);

  ANWriteYieldFileCmd = new G4UIcmdWithAString("/NMS/source/alphan/writeYieldFile", this);
  ANWriteYieldFileCmd->SetGuidance("Output alpha,n reaction yield for source material to file");

  ANSampleDirectionCmd = new G4UIcmdWithAnInteger("/NMS/source/alphan/directionSampling", this);
  ANSampleDirectionCmd->SetGuidance("Set method for sampling the position");
  ANSampleDirectionCmd->SetGuidance(" of neutrons from (alpha, n) reactions");
  ANSampleDirectionCmd->SetGuidance("0 = Isotropic");
  ANSampleDirectionCmd->SetGuidance("1 = use precalculated (alpha, n) source file");
  ANSampleDirectionCmd->SetGuidance("2 = From CS - EnergyAngular Distribution");
  //FIX Does not make sense...
  ANSampleDirectionCmd->SetDefaultValue(AN_DIR_ISOTROPIC);

  ANSampleEnergyCmd = new G4UIcmdWithAnInteger("/NMS/source/alphan/energySampling", this);
  ANSampleEnergyCmd->SetGuidance("Set method for sampling the energy");
  ANSampleEnergyCmd->SetGuidance(" of neutrons from (alpha, n) reactions");
  ANSampleEnergyCmd->SetGuidance("0 = Fixed Energy");
  ANSampleEnergyCmd->SetGuidance("1 = use precalculated (alpha, n) source file");
  ANSampleEnergyCmd->SetGuidance("2 = calculate spectrum using Sources4c-like method");
  //ANSampleEnergyCmd->SetGuidance("3 = From CS - EnergyAngular Distribution");
  ANSampleEnergyCmd->SetGuidance("3 = spectrum read from simple file");
  ANSampleEnergyCmd->SetDefaultValue(AN_E_FIXED);

  ANSamplePositionCmd = new G4UIcmdWithAnInteger("/NMS/source/alphan/positionSampling", this);
  ANSamplePositionCmd->SetGuidance("Set method for sampling the position");
  ANSamplePositionCmd->SetGuidance(" of neutrons from (alpha, n) reactions");
  ANSamplePositionCmd->SetGuidance("0 = isotropic in source volume");
  ANSamplePositionCmd->SetGuidance("1 = use precalculated (alpha, n) source file");
  ANSamplePositionCmd->SetDefaultValue(AN_POS_SOURCEVOLUME);

  ANEnergyCmd = new G4UIcmdWithADoubleAndUnit("/NMS/source/alphan/energy", this);
  ANEnergyCmd->SetGuidance("Set energy of neutrons emitted by (alpha, n) reactions");
  ANEnergyCmd->SetGuidance("only effective if energySampling is set to 0 (default)");
  ANEnergyCmd->SetDefaultValue(1);
  ANEnergyCmd->SetDefaultUnit("MeV");
  
  ANFileCmd = new G4UIcmdWithAString("/NMS/source/alphan/filename", this);
  ANFileCmd->SetGuidance("Set file for input alpha,n data");

  ANEnergyFileCmd = new G4UIcmdWithAString("/NMS/source/alphan/energyFilename", this);
  ANEnergyFileCmd->SetGuidance("Set file with neutron energy spectrum for alpha,n reaction");

  ANWriteEnergyFileCmd = new G4UIcmdWithAString("/NMS/source/alphan/writeEnergyFile", this);
  ANWriteEnergyFileCmd->SetGuidance("Output alpha,n reaction spectrum to file");

  ANSpectrumCalcMT91 = new G4UIcmdWithABool("/NMS/source/alphan/energySpectrumIncludeMT91", this);
  ANSpectrumCalcMT91->SetGuidance("Turn on/off inclusion of MT91 cross section (and q-value) for spectrum calculation");
  ANSpectrumCalcMT91->SetGuidance("Only relevant when energySampling set to 2");

  ANActivityCalcCmd = new G4UIcmdWithAnInteger("/NMS/source/alphan/activityCalculation", this);
  ANActivityCalcCmd->SetGuidance("Set method for calculation of activity for");
  ANActivityCalcCmd->SetGuidance(" (alpha, n) reactions. Activity is the");
  ANActivityCalcCmd->SetGuidance(" number of neutrons emitted per (s * cm3)");
  ANActivityCalcCmd->SetGuidance("0 = Fixed Activity");
  ANActivityCalcCmd->SetGuidance("1 = use precalculated (alpha, n) source file");
  ANActivityCalcCmd->SetGuidance("2 = use (alpha,n) neutron yield calculation");
  
  ANActivityCmd = new G4UIcmdWithADouble("/NMS/source/alphan/activity", this);
  ANActivityCmd->SetGuidance("Set activity of (alpha, n) neutron source");
  ANActivityCmd->SetGuidance("only effective if activityCalculation is set to 0");

  posTypeCmd = new G4UIcmdWithAString("/NMS/source/geometry/type", this);
  posTypeCmd->SetGuidance("Set NMSMaterialDecaySource distribution type");
  posTypeCmd->SetGuidance("Either Point, Surface or Volume");
  posTypeCmd->SetParameterName("Type",false,false);
  posTypeCmd->SetDefaultValue("Point");
  posTypeCmd->SetCandidates("Point Surface Volume");

  posShapeCmd = new G4UIcmdWithAString("/NMS/source/geometry/shape", this);
  posShapeCmd->SetGuidance("Set NMSMaterialDecaySource shape for Surface or Volume source");
  posShapeCmd->SetGuidance("Either Circle, Annulus, Ellipse, Square, Rectangle, Sphere, Ellipsoid, Cylinder or Para");
  posTypeCmd->SetParameterName("Shape",false,false);
  posShapeCmd->SetDefaultValue("NULL");
  posShapeCmd->SetCandidates("Circle Annulus Ellipse Square Rectangle Sphere Ellipsoid Cylinder Para");

  posConfineVolumeCmd = new G4UIcmdWithAString("/NMS/source/geometry/volume", this);
  posConfineVolumeCmd->SetGuidance("Confine source to volume (NULL to unset)");
  posConfineVolumeCmd->SetParameterName("VolName", false, false);
  posConfineVolumeCmd->SetDefaultValue("NULL");

  posRadiusCmd = new G4UIcmdWithADoubleAndUnit("/NMS/source/geometry/radius", this);
  posRadiusCmd->SetGuidance("Set radius.");
  posRadiusCmd->SetParameterName("Radius", false, false);
  posRadiusCmd->SetDefaultUnit("cm");

  posRadius0Cmd = new G4UIcmdWithADoubleAndUnit("/NMS/source/geometry/inner_radius", this);
  posRadius0Cmd->SetGuidance("Set inner radius when required.");
  posRadius0Cmd->SetParameterName("Radius0", false, false);
  posRadius0Cmd->SetDefaultUnit("cm");

  posHalfXCmd = new G4UIcmdWithADoubleAndUnit("/NMS/source/geometry/halfx", this);
  posHalfXCmd->SetGuidance("Set inner radius when required.");
  posHalfXCmd->SetParameterName("HalfX", false, false);
  posHalfXCmd->SetDefaultUnit("cm");

  posHalfYCmd = new G4UIcmdWithADoubleAndUnit("/NMS/source/geometry/halfy", this);
  posHalfYCmd->SetGuidance("Set inner radius when required.");
  posHalfYCmd->SetParameterName("HalfY", false, false);
  posHalfYCmd->SetDefaultUnit("cm");

  posHalfZCmd = new G4UIcmdWithADoubleAndUnit("/NMS/source/geometry/halfz", this);
  posHalfZCmd->SetGuidance("Set z half length.");
  posHalfZCmd->SetParameterName("HalfZ", false, false);
  posHalfZCmd->SetDefaultUnit("cm");

  posCentreCoordsCmd = new G4UIcmdWith3VectorAndUnit("/NMS/source/geometry/centre", this);
  posCentreCoordsCmd->SetGuidance("Set centre coordinates of NMSMaterialDecaySource");
  posCentreCoordsCmd->SetParameterName("X", "Y", "Z", false, false);
  posCentreCoordsCmd->SetDefaultUnit("cm");

  dumpSourceStatusCmd = new G4UIcmdWithoutParameter("/NMS/source/dumpstatus", this);
  dumpSourceStatusCmd->SetGuidance("Output current source status");

}

NMSMaterialDecaySourceMessenger::~NMSMaterialDecaySourceMessenger() {
  delete dumpSourceStatusCmd;
  
  delete posCentreCoordsCmd;
  delete posHalfZCmd;
  delete posHalfYCmd;
  delete posHalfXCmd;
  delete posRadius0Cmd;
  delete posRadiusCmd;
  delete posConfineVolumeCmd;
  delete posShapeCmd;
  delete posTypeCmd;

  delete ANActivityCmd;
  delete ANActivityCalcCmd;
  delete ANSpectrumCalcMT91;
  delete ANWriteEnergyFileCmd;
  delete ANEnergyFileCmd;
  delete ANFileCmd;
  delete ANEnergyCmd;
  delete ANSamplePositionCmd;
  delete ANSampleEnergyCmd;
  delete ANSampleDirectionCmd;
  delete ANWriteYieldFileCmd;
    
  delete AlphaNCmd;
  delete AlphaCmd;
  delete BetaCmd;
  delete SFgCmd;
  delete SFnCmd;

  delete ActivityFixedCmd;
  delete ActivityCmd;
  delete ActiveVolumeCmd;
  delete Cf252Cmd;
  delete SourceFromVolumeCmd;
  delete MaterialCmd;

  delete sourceANDir;
  delete sourcePosDistDir;
  delete sourceDir;
 
}

void NMSMaterialDecaySourceMessenger::SetNewValue(G4UIcommand* cmd, G4String newval) {

  if(cmd == VerboseCmd) {
    nmsmds->SetVerboseLevel(VerboseCmd->GetNewIntValue(newval));
  }
  if(cmd == MaterialCmd) {
    G4String materialname = newval;
    G4Material* sourcematerial = G4Material::GetMaterial(materialname, true);
    if(sourcematerial != 0) {
      nmsmds->SetSourceMaterial(sourcematerial);
    }
  }
  if(cmd == SourceFromVolumeCmd) {
    G4String volumename = newval;
    nmsmds->SetSourceFromPhysicalVolumeName(volumename);
  }
  if(cmd == Cf252Cmd) {
    std::istringstream is(newval);
    G4int ndist;
    G4int neng;
    is >> ndist >> neng;
    nmsmdsSettings->SetCf252n(ndist, neng);
  }

  if(cmd == ActivityCmd) {
    nmsmds->SetActivity(ActivityCmd->GetNewDoubleValue(newval));
    nmsmdsSettings->SetActivityFixed(true);
  }
  if(cmd == ActivityFixedCmd) {
    nmsmdsSettings->SetActivityFixed(ActivityFixedCmd->GetNewBoolValue(newval));
  }
  if(cmd == ActiveVolumeCmd) {
    nmsmds->SetActiveVolume(ActiveVolumeCmd->GetNewDoubleValue(newval));
  }

  if(cmd == SFnCmd) {
    nmsmds->SetSpontaneousFission(SFnCmd->GetNewBoolValue(newval), nmsmds->GetSpontaneousFissionGamma());
  }
  if(cmd == SFgCmd) {
    nmsmds->SetSpontaneousFission(nmsmds->GetSpontaneousFissionGamma(), SFgCmd->GetNewBoolValue(newval));
  }
  if(cmd == BetaCmd) {
    nmsmds->SetBetaDecay(BetaCmd->GetNewBoolValue(newval));
  }
  if(cmd == AlphaCmd) {
    nmsmds->SetAlphaDecay(AlphaCmd->GetNewBoolValue(newval));
  }
  if(cmd == AlphaNCmd) {
    nmsmds->SetAlphaNSource(AlphaNCmd->GetNewBoolValue(newval));
  }

  if(cmd == ANWriteYieldFileCmd) {
    nmsmds->WriteANYield(newval);
  }
  if(cmd == ANSampleDirectionCmd) {
    G4int setting = ANSampleDirectionCmd->GetNewIntValue(newval);
    if(setting > 2) { //FIX
      G4cout << "This is not a valid direction sampling method." << G4endl;
    }
    else {
      nmsmdsSettings->SetANDirectionSamplingMode(NMSANDirectionSamplingMode(setting));
    }
  }
  if(cmd == ANSampleEnergyCmd) {
    G4int setting = ANSampleEnergyCmd->GetNewIntValue(newval);
    if(setting > 3) { //FIX
      G4cout << "This is not a valid energy sampling method." << G4endl;
    }
    else {
      nmsmdsSettings->SetANEnergySamplingMode(NMSANEnergySamplingMode(setting));
    }
    
  }
  if(cmd == ANSamplePositionCmd) {
    G4int setting = ANSamplePositionCmd->GetNewIntValue(newval);
    if(setting > 1) { //FIX
      G4cout << "This is not a valid position sampling method." << G4endl;
    }
    else {
      nmsmdsSettings->SetANPositionSamplingMode(NMSANPositionSamplingMode(setting));
    }
    
  }
  if(cmd == ANFileCmd) {
    nmsmdsSettings->SetANFilename(newval);
  }
  if(cmd == ANEnergyFileCmd) {
    nmsmdsSettings->SetANEnergyFilename(newval);
  }
  if(cmd == ANWriteEnergyFileCmd) {
    nmsmds->WriteANSpectrum(newval);
  }
  if(cmd == ANSpectrumCalcMT91) {
    nmsmdsSettings->SetANSpectrumMT91(ANSpectrumCalcMT91->GetNewBoolValue(newval));
  }
  if(cmd == ANEnergyCmd) {
    nmsmdsSettings->SetANEnergy(ANEnergyCmd->GetNewDoubleValue(newval));
  }
  if(cmd == ANActivityCalcCmd) {
    G4int setting = ANActivityCalcCmd->GetNewIntValue(newval);
    if(setting > 2) { //FIX
      G4cout << "This is not a valid Activity Calc Value" << G4endl;
    }
    else {
      nmsmdsSettings->SetANActivityCalcMode(NMSANActivityCalcMode(setting));
    }
  }
  if(cmd == ANActivityCmd) {
    nmsmdsSettings->SetANActivity(ANActivityCmd->GetNewDoubleValue(newval));
  }



  if(cmd == posTypeCmd) {
    nmsmds->GetPosDist()->SetPosDisType(newval);
  }
  if(cmd == posShapeCmd) {
    nmsmds->GetPosDist()->SetPosDisShape(newval);
  }
  if(cmd == posRadiusCmd) {
    nmsmds->GetPosDist()->SetRadius(posRadiusCmd->GetNewDoubleValue(newval));
  }
  if(cmd == posRadius0Cmd) {
    nmsmds->GetPosDist()->SetRadius0(posRadius0Cmd->GetNewDoubleValue(newval));
  }
  if(cmd == posHalfXCmd) {
    nmsmds->GetPosDist()->SetHalfX(posHalfXCmd->GetNewDoubleValue(newval));
  }
  if(cmd == posHalfYCmd) {
    nmsmds->GetPosDist()->SetHalfY(posHalfYCmd->GetNewDoubleValue(newval));
  }
  if(cmd == posHalfZCmd) {
    nmsmds->GetPosDist()->SetHalfZ(posHalfZCmd->GetNewDoubleValue(newval));
  }
  if(cmd == posCentreCoordsCmd) {
    nmsmds->GetPosDist()->SetCentreCoords(posCentreCoordsCmd->GetNew3VectorValue(newval));
  }

  if(cmd == dumpSourceStatusCmd) {
    nmsmds->DumpSourceStatus();
  }

}
