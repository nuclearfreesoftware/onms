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


#define nZANSepEnergy 74 /* Neutron separation energy for the 73"+1" isotopes 
                            in the paper "Energy-Dependent Fission Q 
                            Values Generalized for All Actinides" by
                            R. Vogt, in which we have average neutron and
                            photon multiplicity distributions as a function
                            of the incident neutron energy. 73"+1" because
                            we added Pa-234 (this isotope was needed to 
                            complete the list of isotopes required for the 
                            Watt spectrum. We have Watt spectrum data for 
                            Pa-233, so we needed the neutron separation 
                            energy Sn for Pa-234).

                            The reference Richard B. Firestone, "Table of 
                            Isotopes," Eigth edition, John Wiley and Sons, 
                            Inc. (1996) did not have the neutron separation
                            energy for U-241, and this is the only isotope
                            for which we had Watt spectrum coefficients for 
                            the corresponding neutron induced fission 
                            isotope U-240 but no separation energy, among 
                            the 40 isotopes for which we have Watt spectrum 
                            coefficients. So we decided to look for another
                            reference where to find the separation energy of
                            U-241 for the sake of having a complete list of 
                            40 isotopes. Fortunately, the two references
                            A.H. Wapstra, G. Audi, C. Thibault, "The 
                            AME2003 atomic mass evaluation (I). Evaluation 
                            of input data, adjustment procedures," Nuclear 
                            physics A729, 129 (2003) and 
                            G. Audi, A.H. Wapstra, C. Thibault, "The 
                            AME2003 atomic mass evaluation (II). Tables, 
                            graphs, and references," Nuclear physics A729, 
                            337 (2003) 
                            had the estimated mass of U-241. This estimated
                            mass was used to produced the neutron 
                            separation energy of U-241. */

#include <string>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include "fissionEvent.h"

double fissionEvent::getNSepEng(int iso) {

/*
  Description
    Returns the neutron separation energy S_n [M(N)-M(N-1)-M(n)] in units of MeV.
    The S_n are taken directly from the reference Richard B. Firestone, "Table of 
    Isotopes," Eigth edition, John Wiley and Sons, Inc. (1996), except for U-241
    which comes from
     A.H. Wapstra, G. Audi, C. Thibault, "The AME2003 atomic mass evaluation (I). 
     Evaluation of input data, adjustment procedures," Nuclear physics A729, 129 
     (2003) and 
     G. Audi, A.H. Wapstra, C. Thibault, "The AME2003 atomic mass evaluation (II). 
     Tables, graphs, and references," Nuclear physics A729, 337 (2003) 
*/

/*
  Input
    iso       - isotope ZA (e.g. 92235 for U-235)
  Output
    neutron separation energy S_n in MeV
*/
   static int ZANSepEng [nZANSepEnergy]= {
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
                   91234, // Pa-234

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

   static double NSepEng [nZANSepEnergy] = {
                   6.663, // Ac-225
                   5.399,  // Ac-226
                   6.527,  // Ac-227
                   
                   5.455,  // Th-227
                   7.1101, // Th-228
                   5.255,  // Th-229
                   6.7936, // Th-230
                   5.11804,// Th-231
                   6.4381, // Th-232
                   4.78635,// Th-233
                   6.189,  // Th-234
                   
                   7.051,  // Pa-229
                   5.800,  // Pa-230
                   6.817,  // Pa-231
                   5.553,  // Pa-232
                   6.527,  // Pa-233
                   5.217,  // Pa-234
                   
                   7.672,  // U-230
                   5.900,  // U-231
                   7.250,  // U-232
                   5.760,  // U-233
                   6.8437, // U-234
                   5.29784,// U-235
                   6.5448, // U-236
                   5.1259, // U-237
                   6.1520, // U-238
                   4.80626,// U-239
                   5.933,  // U-240
                   4.589,  // U-241
                   
                   6.270,  // Np-234
                   6.984,  // Np-235
                   5.730,  // Np-236
                   6.580,  // Np-237
                   5.48809,// Np-238
                   6.2168, // Np-239
                   
                   7.380,  // Pu-236
                   5.8775, // Pu-237
                   7.0005, // Pu-238
                   5.6465, // Pu-239
                   6.5335, // Pu-240
                   5.24160,// Pu-241
                   6.3094, // Pu-242
                   5.034,  // Pu-243
                   6.021,  // Pu-244
                   5.780,  // Pu-246
                   
                   5.957,  // Am-240
                   6.641,  // Am-241
                   5.53757,// Am-242
                   6.3670, // Am-243
                   5.3637, // Am-244
                   
                   7.440,  // Cm-240
                   6.0898, // Cm-241
                   6.9699, // Cm-242
                   5.6933, // Cm-243
                   6.8007, // Cm-244
                   5.5198, // Cm-245
                   6.4580, // Cm-246
                   5.156,  // Cm-247
                   6.213,  // Cm-248
                   4.7135, // Cm-249
                   5.832,  // Cm-250
                   
                   6.970,  // Bk-245
                   5.920,  // Bk-246
                   6.550,  // Bk-247
                   5.451,  // Bk-248
                   6.331,  // Bk-249
                   4.970,  // Bk-250
                   
                   7.364,  // Cf-246
                   6.967,  // Cf-248
                   5.585,  // Cf-249
                   6.6247, // Cf-250
                   5.109,  // Cf-251
                   6.172,  // Cf-252
                   4.806   // Cf-253
                            };

// Find neutron separation energy for nuclide in table
   for (int i=0; i<nZANSepEnergy; i++) {
      if (iso == ZANSepEng[i]) {
         return NSepEng[i];
      }
   }
   std::string errMsg = "No neutron separation energy Sn for nuclide "; errMsg += iso;
   fissionerr(6, "NSepEnergy", errMsg);
   return -1.;;
}
