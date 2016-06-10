/* 
Copyright (c) 2006-2015 Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory 
UCRL-CODE-224807.

All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

o   Redistributions of source code must retain the above copyright notice, this list of conditions and the disclaimer below.

o  Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the disclaimer (as noted below) in the documentation and/or other materials provided with the distribution.

o  Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
Additional BSD Notice

1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE. 

2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or process disclosed, or represents that its use would not infringe privately-owned rights. 

3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence Livermore National Security, LLC, and shall not be used for advertising or product endorsement purposes.
*/


#include "fissionEvent.h"

#define nSPfissIso 16
#define nSPfissn 11

int fissionEvent::SmpSpNuDistData(int isotope, int Cf252option) {

/*
  Description
    Sample Number of Neutrons from spontaneous fission 
    (a) from the neutron multiplicity data for 
        U-238, Pu-240, Pu-242, Cm-242, Cm-244
           using Holden and Zucker's tabulated data
        Pu-236, Pu-238
           using the tabulated data from P.Santi, D.H. 
           Beddingfield and D.R. Mayo: "Revised prompt 
           neutron emission multiplicity distributions 
           for 236,238Pu"
        Cm-246, Cm-248
           using the tabulated data from BNL-36467
        Cf-246
           using the tabulated data from Dakavoski, 
           Sov.Atom.Erg. 17, 360, BNL-36467, (1973).
        Cf-250, Cf-254, Fm-257
           using the tabulated data from Hoffman, 
           Phys.Rev.C. 21, 637, BNL-36467, (1980).
        Cf-252 using either Spencer's tabulated data or 
           Boldeman's data
        No-252
           using the tabulated data from Lazarev, 
           Phys.Lett. 52B, 321, BNL-36467, (1974).
    (b) from Terrell's approximation using nubar for
        Th-232, 
        U-232, U-233, U-234, U-235, U-236,
        Np-237, 
        Pu-239, Pu-241, 
        Am-241, 
        Bk-249
           using Ensslin's data.
*/

/*
  Input
    iso          - isotope
    Cf252option  - 0 to use Spencer's tabulated data
                   1 to use Boldeman's data
  Output
    SmpSpNuDistData - sampled multiplicity
                      -1 is the isotope has 
                         no multiplicity data,
                         nor any nubar data
*/
 
  int i, index;
  double sum, nubar;
  double r;

  static double sfnu [nSPfissIso][nSPfissn] = { 
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
/*
  sample the spontaneous fission neutron number distribution
*/
  index = -1;

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

  if (index != -1) { 
    r=fisslibrng();

    sum=0.;
    for(i=0; i<nSPfissn; i++) {
       sum=sum+sfnu[index][i];
       if (r <= sum) return i;
       if (sum>0. && sfnu[index][i] == 0.) return i-1;
    }
    return nSPfissn-1;
  } else {
// There is no full multiplicity distribution data available
// for that isotope, let's try to find a nubar for it in
// N. Ensslin, et.al., "Application Guide to Neutron
// Multiplicity Counting," LA-13422-M (November 1998)
// and use Terrell's approximation
    nubar = SmpSpNubarData(isotope);
    if (nubar != -1.) {
      return (int) SmpTerrell(nubar);
    } else {
// There is no nubar information for that isotope, return -1,
// meaning no data available for that isotope
      return -1;
    }
  }
  return 0;
}
