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

#include "NMSDetectorConstruction.hh"

G4VPhysicalVolume* NMSDetectorConstruction::Construct() {
  G4VPhysicalVolume* worldPhysical;

  detectorparser.Read(detectorfilename);
  worldPhysical = detectorparser.GetWorldVolume();
  G4LogicalVolume* worldLogical = detectorparser.GetVolume(LVNAME_DETECTOR_CAVITY);

  if(samplefilename != "none") {
    sampleparser.Read(samplefilename);
    sampleposition = sampleparser.GetPosition(GDML_SAMPLE_POSITION);
    G4LogicalVolume* sampleLogical = sampleparser.GetVolume(LVNAME_SAMPLE_MOTHER);
    G4VisAttributes * sampleVisAtt = new G4VisAttributes(G4Colour(1.0,0.,0.));
    sampleLogical->SetVisAttributes(sampleVisAtt);
    G4PVPlacement * samplePhysical = new G4PVPlacement(0, sampleposition, sampleLogical, "Sample", worldLogical, false, 700);
  }


  G4cout<<*(G4Material::GetMaterialTable()) <<G4endl;

  return worldPhysical;
}

NMSDetectorConstruction::NMSDetectorConstruction(G4String detector, G4String sample) : detectorfilename(detector), samplefilename(sample) {

}

NMSDetectorConstruction::~NMSDetectorConstruction() {

}
