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


#include <math.h>
#include "fissionEvent.h"

int fissionEvent::SmpNuDistDataPu239_241_MC(double nubar) {

/*
  Description
    Sample Number of Neutrons from fission in Pu-239 and Pu-241 using 
    Zucker and Holden's tabulated data for Pu-239
    The 11 P(nu) distributions are given as a function of nubar, 
    the average number of neutrons from induced fission for the 
    11 different energies (0 to 10 MeV), based on the Pu-239 data 
    from Zucker and Holden.
*/

/*
  Input
    nubar    - average number of neutrons per fission
  Output
    SmpNuDistDataPu239_241_MC  - sampled multiplicity
    
*/

  static double Pu239nu [11] [9] = {
     {.0108826, .0994916, .2748898, .3269196, .2046061, .0726834, .0097282, .0006301, .0001685},
     {.0084842, .0790030, .2536175, .3289870, .2328111, .0800161, .0155581, .0011760, .0003469},
     {.0062555, .0611921, .2265608, .3260637, .2588354, .0956070, .0224705, .0025946, .0005205},
     {.0045860, .0477879, .1983002, .3184667, .2792811, .1158950, .0301128, .0048471, .0007233},
     {.0032908, .0374390, .1704196, .3071862, .2948565, .1392594, .0386738, .0078701, .0010046},
     {.0022750, .0291416, .1437645, .2928006, .3063902, .1641647, .0484343, .0116151, .0014149},
     {.0014893, .0222369, .1190439, .2756297, .3144908, .1892897, .0597353, .0160828, .0029917},
     {.0009061, .0163528, .0968110, .2558524, .3194566, .2134888, .0729739, .0213339, .0020017},
     {.0004647, .0113283, .0775201, .2335926, .3213289, .2356614, .0886183, .0274895, .0039531},
     {.0002800, .0071460, .0615577, .2089810, .3200121, .2545846, .1072344, .0347255, .0054786},
     {.0002064, .0038856, .0492548, .1822078, .3154159, .2687282, .1295143, .0432654, .0075217}
    };
  static double Pu239nubar [11] = {
      2.8760000,
      3.0088800,
      3.1628300,
      3.3167800,
      3.4707300,
      3.6246800,
      3.7786300,
      3.9325800,
      4.0865300,
      4.2404900,
      4.3944400
    };
  double fraction, r, cum;
  int engind, nu;

/* 
  Check if nubar is within the range of experimental values
*/
  if(nubar >= Pu239nubar[0] && nubar <= Pu239nubar[10]) {
/*
     Use Zucker and Holden Data
*/
     engind = 1;
     while (nubar > Pu239nubar[engind]){ engind++;}
     fraction = (nubar-Pu239nubar[engind-1])/(Pu239nubar[engind]-Pu239nubar[engind-1]);
     if(fisslibrng() > fraction) engind--;

     r = fisslibrng();
     nu = 0;
     cum = Pu239nu[engind][0];
     while (r > cum && nu < 8){ 
       nu++;
       cum += Pu239nu[engind][nu];
     }
     return nu;
  } else {
/*
     Use Terrell's formula
*/
     return (int) SmpTerrell(nubar);
  }
}
