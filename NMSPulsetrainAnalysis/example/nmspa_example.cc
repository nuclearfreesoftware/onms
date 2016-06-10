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

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <sys/stat.h>

using namespace std;

#include "NMSMultiplicityResult.hh"
#include "NMSPulsetrainManager.hh"


void usage() {
      cout << "Usage:" << endl;
      cout << "<app> [options] <pulsetrain list file>" << endl;
      cout << "    File entries should be in microseconds, floating point numbers are allowed." << endl;
      cout << "    Special options can be set for MCNP Ptrac Files or other base units." << endl << endl;
      cout << "Options:" << endl;
      cout << "   -p <num>    Set pre-delay to <num> µs" << endl;
      cout << "   -l <num>    Set number of shift register steps to <num>" << endl;
      cout << "   -p <num>    Set length of one shift register step to <num> µs" << endl;
      cout << "   -a <num>    Set length of accidential gate to <num> µs" << endl;
      cout << "   -d <num>    Set length of derandomizer period to <num> µs" << endl;
      cout << "   -D          No derandomization" << endl;
      cout << "   -R          No quantization for registers" << endl;
      cout << "   -P          Read file as MCNP Ptrac file" << endl;
      cout << "   -S          Read data entries in unit of shakes" << endl;
      cout << "   -m <num>    Set length of measurement to <num> µs (Default: 10^8 µs)" << endl;
}

inline bool file_exists (const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

int main(int argc,char** argv)
{

  double predelay = 4.5;
  double registerlength = 128;
  double registerperiod = 0.5;
  double agate = 4096;
  double derandomizeperiod = 0.1;
  double measurementlength = 100000000;

  bool gateS = false;
  bool ptracfile = false;
  bool unitshakes = false;
  bool derandomizeDo = true;
  bool registerquantizationDo = true;

  string filename = "";

  cout << "NMSPulsetrainAnalysis Library (version 0.5) - Example runtime" << endl;

  if(argc <= 1) {
      usage();
      exit(0);
  }
  for ( int i=1; i<argc; i=i+1 ) {
    if      ( string(argv[i]) == "-p" ) {
      i++;
      predelay = atof(argv[i]);
    }
    else if ( string(argv[i]) == "-l" ) {
      i++;
      registerlength = atof(argv[i]);
    }
    else if ( string(argv[i]) == "-r" ) {
      i++;
      registerperiod = atof(argv[i]);
    }
    else if ( string(argv[i]) == "-a" ) {
      i++;
      agate = atof(argv[i]);
    }
    else if ( string(argv[i]) == "-d" ) {
      i++;
      derandomizeperiod = atof(argv[i]);
    }
    else if ( string(argv[i]) == "-P" ) ptracfile = true;
    else if ( string(argv[i]) == "-S" ) unitshakes = true;
    else if ( string(argv[i]) == "-D" ) derandomizeDo = false;
    else if ( string(argv[i]) == "-R" ) registerquantizationDo = false;
    else if ( string(argv[i]) == "-m" ) {
      i++;
      measurementlength = atof(argv[i]);
    }
    else {
      if(file_exists(string(argv[i]))) {
	filename = argv[i];
      }
      else {
	cout << "Given pulsetrain list file does not exist." << endl;
	cout << "Filename: " << argv[i] << endl;
	return 1;
      }
    }
  }

  // Shift register length (in number of positions, not time!!!)
  NMSPulsetrainManager * pm = new NMSPulsetrainManager(registerlength);

  // set predelay in µs
  pm->setPredelay(predelay);

  // set accidentials gate delay in µs
  pm->setAdelay(agate);

  // set register period (R+A Gate length) / Shift Register Length
  pm->setRegisterPeriod(registerperiod);

  // set period for derandomize function in µs
  pm->setDerandomizePeriod(derandomizeperiod);

  ifstream input;
  string line;
  double test;

  if(unitshakes) {
    pm->loadPulsetrainFromShakesFile(filename);
  }
  else {
    pm->loadPulsetrainFromFile(filename);
  }

  pm->sortOriginal();

  pm->processEvents();

  cout << "Results using Fast Forward Method" << endl;
  NMSMultiplicityResult mureff = pm->getResults(RESULTS_FAST_FORWARD_SHIFT);

  mureff.normalize();
  cout << mureff;
  cout << endl;
  cout << mureff.singles(1000000 / measurementlength) << "   " << mureff.doubles(1000000 / measurementlength) << "    " << mureff.triples(1000000 / measurementlength) << endl;
  cout << "Results using All Time Steps Method" << endl;
  NMSMultiplicityResult mureats = pm->getResults(RESULTS_ALL_TIME_STEPS_SHIFT);

  mureats.normalize();
  cout << mureats;
  cout << endl;
  cout << mureats.singles(1000000 / measurementlength) << "   " << mureats.doubles(1000000 / measurementlength) << "    " << mureats.triples(1000000 / measurementlength) << endl;

  cout << "Results using Default Method (All Events, No Shift queue)" << endl;
  NMSMultiplicityResult mureae = pm->getResults();

  mureae.normalize();
  cout << mureae;
  cout << endl;
  cout << mureae.singles(1000000 / measurementlength) << "   " << mureae.doubles(1000000 / measurementlength) << "    " << mureae.triples(1000000 / measurementlength) << endl;


  cout << endl;
  cout << "Used data" << endl;
  cout << "Gate width: " << pm->getSettings().gate << endl;
  cout << "Register period: " << registerperiod << endl;
  return 0;

}
