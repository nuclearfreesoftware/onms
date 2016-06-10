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


#define nZAEngConsAllActN 73 /* 73 isotopes for neutrons in paper 
                               "Energy-Dependent Fission Q Values 
                                Generalized for All Actinides" by R. Vogt */

#include <string>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include "fissionEvent.h"

double fissionEvent::getTotEngNEnergyConsAllActinides(double ePart, int iso) {

/*
  Description
    Determine the average total energy of prompt photons emitted by fission.

    The data in this function comes from the paper "Energy-Dependent 
    Fission Q Values Generalized for All Actinides" by R. Vogt. The paper 
    gives the total energy for all emitted prompt fission neutons. Data 
    in the paper is given for all actinides, major and minor, in the 
    Evaluated Nuclear Data Library 2008 release, ENDL2008.
    The average total prompt fission neutron energy is given as a 
    function of the incident neutron energy by the formula:
      NEng = coeffNeutron[2]+coeffNeutron[1]*ePart+coeffNeutron[0]*ePart^2
*/

/*
  Input
    ePart     - energy of incoming particle
    iso       - isotope
  Output
    getTotEngNEnergyConsAllActinides - average total energy of prompt fission
                                       gamma-rays
*/

   static int ZAEngConsAllActN [nZAEngConsAllActN]= {
                   89225, // Ac-225
                   89226, // Ac-226
                   89227, // Ac-227

                   90227, // Th-227
                   90228, // Th-228
                   90229, // Th-229
                   90230, // Th-230
                   90231, // Th-231
                   90232, // Th-232
                   90233, // Th-233
                   90234, // Th-234

                   91229, // Pa-229
                   91230, // Pa-230
                   91231, // Pa-231
                   91232, // Pa-232
                   91233, // Pa-233

                   92230, // U-230
                   92231, // U-231
                   92232, // U-232
                   92233, // U-233
                   92234, // U-234
                   92235, // U-235
                   92236, // U-236
                   92237, // U-237
                   92238, // U-238
                   92239, // U-239
                   92240, // U-240
                   92241, // U-241

                   93234, // Np-234
                   93235, // Np-235
                   93236, // Np-236
                   93237, // Np-237
                   93238, // Np-238
                   93239, // Np-239

                   94236, // Pu-236
                   94237, // Pu-237
                   94238, // Pu-238
                   94239, // Pu-239
                   94240, // Pu-240
                   94241, // Pu-241
                   94242, // Pu-242
                   94243, // Pu-243
                   94244, // Pu-244
                   94246, // Pu-246

                   95240, // Am-240
                   95241, // Am-241
                   95242, // Am-242
                   95243, // Am-243
                   95244, // Am-244

                   96240, // Cm-240
                   96241, // Cm-241
                   96242, // Cm-242
                   96243, // Cm-243
                   96244, // Cm-244
                   96245, // Cm-245
                   96246, // Cm-246
                   96247, // Cm-247
                   96248, // Cm-248
                   96249, // Cm-249
                   96250, // Cm-250

                   97245, // Bk-245
                   97246, // Bk-246
                   97247, // Bk-247
                   97248, // Bk-248
                   97249, // Bk-249
                   97250, // Bk-250

                   98246, // Cf-246
                   98248, // Cf-248
                   98249, // Cf-249
                   98250, // Cf-250
                   98251, // Cf-251
                   98252, // Cf-252
                   98253  // Cf-253
                      };

   static double coeffNeutron [nZAEngConsAllActN][3] = {
                  {-0.001317, 0.1937, 3.478}, // Ac-225
                  {0.004442, 0.1231, 3.635},  // Ac-226
                  {-0.000144, 0.1888, 3.396}, // Ac-227
                   
                  {0.006569, 0.1225, 4.275},  // Th-227
                  {0.003449, 0.2181, 3.787},  // Th-228
                  {0.006267, 0.1339, 4.216},  // Th-229
                  {0.007380, 0.1422, 3.847},  // Th-230
                  {0.006487, 0.1196, 4.095},  // Th-231
                  {-0.000431, 0.3465, 3.401}, // Th-232
                  {0.000663, 0.2566, 3.736},  // Th-233
                  {0.003476, 0.2290, 3.387},  // Th-234
                   
                  {0.005433, 0.1744, 4.605},  // Pa-229
                  {0.005562, 0.1879, 4.720},  // Pa-230
                  {0.006436, 0.1726, 4.524},  // Pa-231
                  {0.006763, 0.1683, 4.699},  // Pa-232
                  {0.000639, 0.3671, 4.076},  // Pa-233
                   
                  {0.006792, 0.1832, 4.977},  // U-230
                  {0.005808, 0.2127, 5.196},  // U-231
                  {0.003243, 0.2782, 6.082},  // U-232
                  {0.002915, 0.2540, 5.141},  // U-233
                  {0.002704, 0.2339, 4.728},  // U-234
                  {-0.001424, 0.3114, 4.864}, // U-235
                  {0.004555, 0.2969, 4.505},  // U-236
                  {0.001783, 0.2680, 4.999},  // U-237
                  {-0.004351, 0.3574, 4.509}, // U-238
                  {0.004266, 0.3647, 4.580},  // U-239
                  {0.000273, 0.3596, 4.561},  // U-240
                  {0.002821, 0.3998, 4.268},  // U-241
                   
                  {0.007642, 0.2311, 5.880},  // Np-234
                  {0.007751, 0.2484, 5.576},  // Np-235
                  {0.008116, 0.2446, 5.080},  // Np-236
                  {0.005819, 0.2768, 5.330},  // Np-237
                  {0.007559, 0.2650, 5.214},  // Np-238
                  {0.004159, 0.2489, 5.416},  // Np-239
                   
                  {0.009279, 0.2240, 6.112},  // Pu-236
                  {0.006790, 0.2599, 6.177},  // Pu-237
                  {0.008211, 0.2189, 6.087},  // Pu-238
                  {-0.002495, 0.3707, 6.092}, // Pu-239
                  {0.008608, 0.2477, 5.906},  // Pu-240
                  {0.009310, 0.2356, 6.161},  // Pu-241
                  {0.008356, 0.2192, 5.926},  // Pu-242
                  {0.005751, 0.4692, 5.781},  // Pu-243
                  {0.008807, 0.2557, 5.655},  // Pu-244
                  {0.007922, 0.3155, 5.145},  // Pu-246
                   
                  {0.002294, 0.3473, 7.150},  // Am-240
                  {-0.004504, 0.4243, 6.957}, // Am-241
                  {0.002294, 0.3473, 7.150},  // Am-242
                  {-0.002387, 0.3523, 7.422}, // Am-243
                  {0.0000, 0.3837, 6.543},    // Am-244
                   
                  {0.011040, 0.2786, 7.525},  // Cm-240
                  {0.007316, 0.3648, 7.699},  // Cm-241
                  {0.011400, 0.2683, 7.701},  // Cm-242
                  {0.005492, 0.2363, 8.104},  // Cm-243
                  {0.010830, 0.2061, 7.103},  // Cm-244
                  {0.005426, 0.2279, 7.984},  // Cm-245
                  {0.009390, 0.2245, 6.939},  // Cm-246
                  {0.008595, 0.3896, 8.216},  // Cm-247
                  {0.013550, 0.2499, 7.295},  // Cm-248
                  {0.008907, 0.3777, 7.124},  // Cm-249
                  {0.006831, 0.4062, 6.973},  // Cm-250
                   
                  {0.009615, 0.3643, 8.210},  // Bk-245
                  {0.005445, 0.4764, 8.274},  // Bk-246
                  {0.008129, 0.4266, 7.831},  // Bk-247
                  {0.006656, 0.4796, 8.145},  // Bk-248
                  {0.010130, 0.4021, 7.519},  // Bk-249
                  {0.008308, 0.4204, 7.879},  // Bk-250
                   
                  {0.009000, 0.4323, 8.900},  // Cf-246
                  {0.010700, 0.3877, 8.661},  // Cf-248
                  {0.007067, 0.4746, 9.428},  // Cf-249
                  {0.007397, 0.4980, 8.226},  // Cf-250
                  {0.010790, 0.4454, 9.407},  // Cf-251
                  {0.007184, 0.5190, 8.627},  // Cf-252
                  {0.018650, 0.2396, 8.449}   // Cf-253
                            };

   double muEn;

   int i;

// Find neutron parameters for isotope
   int isoindexn=-1;
   for (i=0; isoindexn == -1 && i<nZAEngConsAllActN; i++) {
      if (iso == ZAEngConsAllActN[i]) {
         isoindexn = i;
         break;
      }
   }
   if (isoindexn == -1) {
      std::string errMsg = "No total fission neutron energy available for iso "; errMsg += iso;
      fissionerr(6, "getTotEngNEnergyConsAllActinides", errMsg);
   }
   
// Computing total energies for neutrons
   muEn = coeffNeutron[isoindexn][2] + ePart*(coeffNeutron[isoindexn][1] + ePart*coeffNeutron[isoindexn][0]);

   return muEn;
}
