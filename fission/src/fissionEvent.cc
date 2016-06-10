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
#include <stdio.h>
#include <stdlib.h>

int fissionEvent::delayoption=0;
int fissionEvent::correlationoption=0;
int fissionEvent::nudistoption=3;
int fissionEvent::Cf252ndistoption=0;
int fissionEvent::Cf252nengoption=0;

fissionEvent::fissionEvent(int isotope, double time, double nubar, double eng, int fissiontype) {
   /*
    * Constructs a fission event with neutronNu neutrons and photonNu
    * photons.
    */
   int i;

   bool spontaneous = (fissiontype==0);

   neutronNu = 0;
   photonNu = 0;

   effectivecorrelationoption = correlationoption;

   if (effectivecorrelationoption == 3) {
      bool success = false;
#ifdef FREYA
      success = SmpFreya(eng, isotope, fissiontype);
#endif
      if(success) {
         for (i=0; i<photonNu; i++) photonAges[i] = time;
         for (i=0; i<neutronNu; i++) neutronAges[i] = time;
         return;
      }
      effectivecorrelationoption = 0;
   }

   // For photofission, adjust isotope ZA and excitation energy (except for FREYA)
   if (fissiontype==2 && effectivecorrelationoption!=3) {
     eng-=fissionEvent::getNSepEng(isotope); // subtract neutron separation energy
     isotope--;                              // use Z(A-1)
     if (eng<0) eng=0.;
   }

   if (spontaneous) {
      /* spontaneous fission */
      neutronNu = SmpSpNuDistData(isotope, Cf252ndistoption);
      photonNu = SmpSpNugDistData(isotope);
   } else {
     /* induced fission */
     switch (nudistoption) {
       case 0: // nudistoption=0
       case 1: // nudistoption=1
         switch (isotope) {
           case 92235:
              neutronNu = SmpNuDistDataU235(eng,nudistoption);
              break;
           case 92238:
              neutronNu = SmpNuDistDataU238(eng);
              break;
           case 94239:
              neutronNu = SmpNuDistDataPu239(eng);
              break;
           default:
              neutronNu = (int) SmpTerrell(nubar);
              break;
         } 
         break;
       case 2: // nudistoption=2
         switch (isotope) {
           case 92232:
           case 92234:
           case 92236:
           case 92238:
              neutronNu = SmpNuDistDataU232_234_236_238(nubar);
              break;
           case 92233:
           case 92235:
              neutronNu = (int) SmpNuDistDataU233_235(nubar);
              break;
           case 94239:
           case 94241:
              neutronNu = SmpNuDistDataPu239_241(nubar);
              break;
           default:
              neutronNu = (int) SmpTerrell(nubar);
              break;
         }
         break;
       case 3: // nudistoption=3
         switch (isotope) {
           case 92232:
           case 92234:
           case 92236:
           case 92238:
              neutronNu = SmpNuDistDataU232_234_236_238_MC(nubar);
              break;
           case 92233:
           case 92235:
              neutronNu = (int) SmpNuDistDataU233_235_MC(nubar);
              break;
           case 94239:
           case 94241:
              neutronNu = SmpNuDistDataPu239_241_MC(nubar);
              break;
           default:
              neutronNu = (int) SmpTerrell(nubar);
              break;
         } 
         break;
       default: // nudistoption...
         break;
     }
     photonNu = SmpNugDist(isotope, nubar, eng, spontaneous);
   }
   allocateMem(neutronNu, photonNu);
   bool success=false;
   switch (effectivecorrelationoption) {
     case 1: // effectivecorrelationoption=1
       switch (isotope) {
         case 92235:
         case 92238:
         case 94239:
           SmpNPEnergyCons(eng, nubar, isotope, spontaneous);
           success=true;
           break;
         default:
           break;
       }
       break;
     case 2: // effectivecorrelationoption=2
       switch (isotope) {
         case 89225:  // Ac-225
         case 89226:  // Ac-226
         case 89227:  // Ac-227
  
         case 90227:  // Th-227
         case 90228:  // Th-228
         case 90229:  // Th-229
         case 90230:  // Th-230
         case 90231:  // Th-231
         case 90232:  // Th-232
         case 90233:  // Th-233
         case 90234:  // Th-234
  
         case 91229:  // Pa-229
         case 91230:  // Pa-230
         case 91231:  // Pa-231
         case 91232:  // Pa-232
         case 91233:  // Pa-233
  
         case 92230:  // U-230
         case 92231:  // U-231
         case 92232:  // U-232
         case 92233:  // U-233
         case 92234:  // U-234
         case 92235:  // U-235
         case 92236:  // U-236
         case 92237:  // U-237
         case 92238:  // U-238
         case 92239:  // U-239
         case 92240:  // U-240
         case 92241:  // U-241

         case 93234:  // Np-234
         case 93235:  // Np-235
         case 93236:  // Np-236
         case 93237:  // Np-237
         case 93238:  // Np-238
         case 93239:  // Np-239

         case 94236:  // Pu-236
         case 94237:  // Pu-237
         case 94238:  // Pu-238
         case 94239:  // Pu-239
         case 94240:  // Pu-240
         case 94241:  // Pu-241
         case 94242:  // Pu-242
         case 94243:  // Pu-243
         case 94244:  // Pu-244
         case 94246:  // Pu-246

         case 95240:  // Am-240
         case 95241:  // Am-241
         case 95242:  // Am-242
         case 95243:  // Am-243
         case 95244:  // Am-244

         case 96240:  // Cm-240
         case 96241:  // Cm-241
         case 96242:  // Cm-242
         case 96243:  // Cm-243
         case 96244:  // Cm-244
         case 96245:  // Cm-245
         case 96246:  // Cm-246
         case 96247:  // Cm-247
         case 96248:  // Cm-248
         case 96249:  // Cm-249
         case 96250:  // Cm-250

         case 97245:  // Bk-245
         case 97246:  // Bk-246
         case 97247:  // Bk-247
         case 97248:  // Bk-248
         case 97249:  // Bk-249
         case 97250:  // Bk-250

         case 98246:  // Cf-246
         case 98248:  // Cf-248
         case 98249:  // Cf-249
         case 98250:  // Cf-250
         case 98251:  // Cf-251
         case 98252:  // Cf-252
         case 98253:  // Cf-253
           SmpNPEnergyCons(eng, nubar, isotope, spontaneous);
           success = true;
           break;
         default:     // isotope
           break;
       }
       break;
     default: // effectivecorrelationoption...
       break;
   }
   if (!success) {
     if (neutronNu > 0) {
       for (i=0; i<neutronNu; i++) {
         if (spontaneous) {
            /* spontaneous fission */
            if (isotope == 98252) neutronEnergies[i] = SmpNEngCf252(Cf252nengoption);
            else neutronEnergies[i] = SmpSpWatt(isotope);
         } else {
            /* induced fission */
            neutronEnergies[i] = SmpWatt(eng, isotope);
         }
       }
     }
     if (photonNu > 0) {
        for (i=0; i<photonNu; i++) photonEnergies[i] = SmpGEng();
     }
   }
   for (i=0; i<neutronNu; i++) {
      neutronVelocities[i] = SmpNVel(neutronEnergies[i]);
      SmpIsoDir(&(neutronDircosu[i]),
                &(neutronDircosv[i]),
                &(neutronDircosw[i])
               );
   }
   for (i=0; i<photonNu; i++) {
      photonVelocities[i] = SmpPVel();
      SmpIsoDir(&(photonDircosu[i]), 
                &(photonDircosv[i]), 
                &(photonDircosw[i])
               );
   }
   for (i=0; i<photonNu; i++) photonAges[i] = time;
   for (i=0; i<neutronNu; i++) neutronAges[i] = time;
};

fissionEvent::~fissionEvent() {
   if (neutronNu > 0) {
      delete [] neutronEnergies;
      delete [] neutronVelocities;
      delete [] neutronDircosu;
      delete [] neutronDircosv;
      delete [] neutronDircosw;
      delete [] neutronAges;
   }

   if (photonNu > 0) {
      delete [] photonEnergies;
      delete [] photonVelocities;
      delete [] photonDircosu;
      delete [] photonDircosv;
      delete [] photonDircosw;
      delete [] photonAges;
   }
};
