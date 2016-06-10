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


#define nZAfis 40    /* 38 fissionable isotopes in ENDL + U-232 and Pu-236 from ENDF.B-VII */
#define WATTEMIN 1.0e-6
#define WATTEMAX 20.0

#include <string>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include "fissionEvent.h"

double fissionEvent::SmpWatt(double ePart, int iso) {

/*
  Description
    Sample Watt Spectrum as in TART (Kalos algorithm)
*/

/*
  Input
    ePart     - energy of incoming particle
    iso       - isotope
  Output
              - energy of incoming particle
*/

   static int nZA [nZAfis]= {
                      90231, 90232, 90233,
                      91233,
                      92232, 92233, 92234, 92235, 92236, 92237, 92238, 92239, 92240,
                      93235, 93236, 93237, 93238,
                      94236, 94237, 94238, 94239, 94240, 94241, 94242, 94243,
                      95241, 95242, 95243,
                      96242, 96243, 96244, 96245, 96246, 96247, 96248,
                      97249,
                      98249, 98250, 98251, 98252
                      };

   static double Watta [nZAfis][3] = {
                           {6.00949285e-05, -8.36695381e-03,  9.50939496e-01}, // 90231
                           {6.54348443e-05, -8.86574327e-03,  9.55404490e-01}, // 90232
                           {7.08173682e-05, -9.22676286e-03,  9.50088329e-01}, // 90233
                           {6.35839062e-05, -8.63645973e-03,  9.24583535e-01}, // 91233
                           {2.12324523e-05, -8.27743356e-03,  9.18556312e-01}, // 92232
                           {6.21335718e-05, -8.45651858e-03,  9.14717276e-01}, // 92233
                           {6.81386135e-05, -8.99142394e-03,  9.21954824e-01}, // 92234
                           {7.32627297e-05, -9.36908697e-03,  9.20107976e-01}, // 92235
                           {8.06505279e-05, -9.95416671e-03,  9.27890410e-01}, // 92236
                           {8.33208285e-05, -1.01073057e-02,  9.17691654e-01}, // 92237
                           {8.96944680e-05, -1.06491070e-02,  9.25496030e-01}, // 92238
                           {9.44608097e-05, -1.08940419e-02,  9.17795511e-01}, // 92239
                           {1.01395704e-04, -1.15098159e-02,  9.29395462e-01}, // 92240
                           {6.81110009e-05, -8.91619352e-03,  9.00047566e-01}, // 93235
                           {7.21126359e-05, -9.20179363e-03,  8.95722889e-01}, // 93236
                           {7.82371142e-05, -9.67050621e-03,  8.99574933e-01}, // 93237
                           {8.27256297e-05, -9.99353009e-03,  8.97461897e-01}, // 93238
                           {1.31388919e-04, -8.01060484e-03,  8.91083550e-01}, // 94236
                           {7.29458059e-05, -9.22415170e-03,  8.80996165e-01}, // 94237
                           {8.02383914e-05, -9.78291439e-03,  8.88964070e-01}, // 94238
                           {8.50641730e-05, -1.01099145e-02,  8.87304833e-01}, // 94239
                           {9.10537157e-05, -1.05303084e-02,  8.89438514e-01}, // 94240
                           {9.43014320e-05, -1.07133543e-02,  8.82632055e-01}, // 94241
                           {1.02655616e-04, -1.13154691e-02,  8.91617174e-01}, // 94242
                           {1.06118094e-04, -1.14971777e-02,  8.85181637e-01}, // 94243
                           {9.08474473e-05, -1.04296303e-02,  8.71942958e-01}, // 95241
                           {9.35633054e-05, -1.05612167e-02,  8.63930371e-01}, // 95242
                           {1.01940441e-04, -1.11573929e-02,  8.73153437e-01}, // 95243
                           {9.19501202e-05, -1.04229157e-02,  8.58681822e-01}, // 96242
                           {9.42991674e-05, -1.05098872e-02,  8.49103546e-01}, // 96243
                           {1.02747171e-04, -1.11371417e-02,  8.60434431e-01}, // 96244
                           {1.05024967e-04, -1.12138980e-02,  8.51101942e-01}, // 96245
                           {1.14130011e-04, -1.18692049e-02,  8.62838259e-01}, // 96246
                           {1.15163673e-04, -1.18553822e-02,  8.51306646e-01}, // 96247
                           {1.27169055e-04, -1.27033210e-02,  8.68623539e-01}, // 96248
                           {1.24195213e-04, -1.24047085e-02,  8.48974077e-01}, // 97249
                           {1.12616150e-04, -1.15135023e-02,  8.19708800e-01}, // 98249
                           {1.23637465e-04, -1.22869889e-02,  8.35392018e-01}, // 98250
                           {1.22724317e-04, -1.21677963e-02,  8.22569523e-01}, // 98251
                           {1.33891595e-04, -1.29267762e-02,  8.37122909e-01}  // 98252
                            };

   double a;  /* Watt Parameters */
   double b = 1.;

   double rand1,rand2;
   double x,y,z;
   double eSmp;
   int i;


/*
   Find Watt parameters for isotope
*/
   int isoindex=-1;
   for (i=0; isoindex == -1 && i<nZAfis; i++) {
      if (iso == nZA[i]) isoindex = i;
   }
   if (isoindex == -1) {
      std::ostringstream o;
      o << iso;
      std::string errMsg = "No Watt spectrum available for iso " + o.str();
      fissionerr(6, "SmpWatt", errMsg);
   }
   
   a= Watta[isoindex][2] + ePart*(Watta[isoindex][1] + ePart*Watta[isoindex][0]);

   x= 1. + (b/(8.*a));
   y= (x + sqrt(x*x-1.))/a;
   z= a*y - 1.;

   do {

      rand1= -log(fisslibrng());
      rand2= -log(fisslibrng());
 
      eSmp= y*rand1;

   } while ((rand2-z*(rand1+1.))*(rand2-z*(rand1+1.)) > b*y*rand1 ||
             eSmp < WATTEMIN || eSmp > WATTEMAX);
 
   
   return eSmp;
}
