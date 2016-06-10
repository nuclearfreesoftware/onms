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

#ifndef NMSAlphaNSet_h
#define NMSAlphaNSet_h 1

#include <vector>
#include <fstream>
#include <iomanip>
#include <sys/stat.h>

#include "G4SystemOfUnits.hh"

#include "NMSAlphaNReaction.hh"

class NMSAlphaNSet : public std::vector<NMSAlphaNReaction> {
public:  
  NMSAlphaNSet();
  virtual ~NMSAlphaNSet();

  void saveToFile(G4String filename);
  void loadFromFile(G4String filename);

private:  
  inline bool file_exists(G4String filename) {
    struct stat buffer;
    return (stat (filename.c_str(), &buffer) == 0);
  }

};

#endif
