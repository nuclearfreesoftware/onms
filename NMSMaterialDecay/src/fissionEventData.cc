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

#include "fissionEventData.hh"

// Modified line from fission library SmpSpNuDistData.cc
double fissionEventData::sfnu [SP_FISSION_ISOTOPES][SP_FISSION_N] = { 
// BEGIN from fission library SmpSpNuDistData.cc
     {0.0481677,0.2485215,0.4253044,0.2284094,0.0423438,0.0072533,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000},
     {0.0631852,0.2319644,0.3333230,0.2528207,0.0986461,0.0180199,0.0020407,0.0000000,0.0000000,0.0000000,0.0000000},
     {0.0679423,0.2293159,0.3341228,0.2475507,0.0996922,0.0182398,0.0031364,0.0000000,0.0000000,0.0000000,0.0000000},
     {0.0212550,0.1467407,0.3267531,0.3268277,0.1375090,0.0373815,0.0025912,0.0007551,0.0001867,0.0000000,0.0000000},
     {0.0150050,0.1161725,0.2998427,0.3331614,0.1837748,0.0429780,0.0087914,0.0002744,0.0000000,0.0000000,0.0000000},
     {0.0562929,0.2106764,0.3797428,0.2224395,0.1046818,0.0261665,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000},
     {0.0021100,0.0246700,0.1229000,0.2714400,0.3076300,0.1877000,0.0677000,0.0140600,0.0016700,0.0001000,0.0000000},
     {0.0020900,0.0262100,0.1262000,0.2752000,0.3018000,0.1846000,0.0668000,0.0150000,0.0021000,0.0000000,0.0000000},
     {0.0802878,0.2126177,0.3773740,0.2345049,0.0750387,0.0201770,0.0000000,0.0000000,0.0000000,0.0000000,0.0000000},
     {0.0152182,0.0762769,0.2627039,0.3449236,0.2180653,0.0755895,0.0072227,0.0000000,0.0000000,0.0000000,0.0000000},
     {0.0067352,0.0596495,0.2205536,0.3509030,0.2543767,0.0893555,0.0167386,0.0016888,0.0000000,0.0000000,0.0000000},
     {0.0005084,0.1135987,0.2345989,0.2742853,0.2208697,0.1259660,0.0301731,0.0000000,0.0000000,0.0000000,0.0000000},
     {0.0038191,0.0365432,0.1673371,0.2945302,0.2982732,0.1451396,0.0472215,0.0040174,0.0031188,0.0000000,0.0000000},
     {0.0001979,0.0190236,0.1126406,0.2638883,0.3183439,0.1941768,0.0745282,0.0150039,0.0021968,0.0000000,0.0000000},
     {0.0205736,0.0520335,0.1172580,0.1997003,0.2627898,0.2007776,0.1061661,0.0333033,0.0073979,0.0000000,0.0000000},
     {0.0569148,0.0576845,0.0924873,0.1437439,0.1832482,0.1831510,0.1455905,0.0962973,0.0382048,0.0026776,0.0000000}
      };
//END from fission library SmpSpNuDistData.cc

// Modified line from fission library SmpSpNubarData.cc
int fissionEventData::spzaid [SP_FISSION_NUBAR_ISOTOPES] = {
// BEGIN from fission library SmpSpNubarData.cc
  90232, 92232, 92233, 92234, 92235,
  92236, 92238, 93237, 94238, 94239, 
  94240, 94241, 94242, 95241, 96242, 
  96244, 97249, 98252 };
// END from fission library SmpSpNubarData.cc

// Modified line from fission library SmpSpNubarData.cc
double fissionEventData::spnubar [SP_FISSION_NUBAR_ISOTOPES] = {
// BEGIN from fission library SmpSpNubarData.cc
  2.14,  1.71, 1.76,  1.81, 1.86,
  1.91,  2.01, 2.05,  2.21, 2.16, 
  2.156, 2.25, 2.145, 3.22, 2.54, 
  2.72,  3.40, 3.757
};
// END from fission library SmpSpNubarData.cc

fissionEventData::fissionEventData() {
  
} 

fissionEventData::~fissionEventData() {
  
}


int fissionEventData::sfDataIndex(int isotope, int Cf252option) {
  int index = -1;

  // BEGIN from fission library SmpSpNuDistData.cc
  if (isotope == 92238) index = 0;
  else if (isotope == 94240) index = 1;
  else if (isotope == 94242) index = 2;
  else if (isotope == 96242) index = 3;
  else if (isotope == 96244) index = 4;
  else if (isotope == 94238) index = 5;
  else if (isotope == 98252 && Cf252option == 0) index = 6;
  else if (isotope == 98252 && Cf252option == 1) index = 7;
  else if (isotope == 94236) index = 8;
  else if (isotope == 96246) index = 9;
  else if (isotope == 96248) index = 10;
  else if (isotope == 98246) index = 11;
  else if (isotope == 98250) index = 12;
  else if (isotope == 98254) index = 13;
  else if (isotope == 100257) index = 14;
  else if (isotope == 102252) index = 15;
  // END from fission library SmpSpNuDistData.cc

  return index;
}

double fissionEventData::fissionEventNu(int isotope, int n, int Cf252option) {
  if(n >= 0 || n < SP_FISSION_N) {
    int index = sfDataIndex(isotope, Cf252option);
    if (index != -1) {
      return sfnu[index][n];
    }
  }
  return -1.;
}

double fissionEventData::getSfNubar(int isotope, int Cf252option) {
  double weightedtotal = 0;
  int index = sfDataIndex(isotope, Cf252option);

  // Similar to fission library files SmpSpNuDistData.cc / SmpSpNubarData.cc
  if (index != -1) { 
    for(int i = 0; i < SP_FISSION_N; i++) {
      weightedtotal += i * sfnu[index][i];
    }
    return weightedtotal;
  }
  else {
    for(int i = 0; i < SP_FISSION_NUBAR_ISOTOPES; i++) {
      if( spzaid[i] == isotope) {
	return spnubar[i];
      }
    } 
  }
  return -1.;
}
