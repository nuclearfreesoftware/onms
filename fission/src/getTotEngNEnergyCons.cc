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


#define nZAEngConsN 3   /* 3 isotopes in paper "Implementation of Energy-
                           Dependent Q Values for Fission" by R. Vogt, 
                           B. Beck, D.A. Brown, F. Daffin, and J. Hedstrom */
#include <string>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include "fissionEvent.h"

double fissionEvent::getTotEngNEnergyCons(double ePart, int iso) {

/*
  Description
    Determine the average total energy of prompt neutrons emitted by fission.

    The data in this function comes from the paper "Implementation of Energy-
    Dependent Q Values for Fission" by R. Vogt, B. Beck, D.A. Brown, F. Daffin,
    and J. Hedstrom. The paper gives the total energy for all emitted prompt 
    fission neutrons. Data in the paper is given for the following 3 isotopes:
    U-235, U-238, Pu-239.
    The average total prompt fission neutron energy is given as a function of 
    the incident neutron energy with the formula:
       NEng = VogtNeutron[0] + VogtNeutron[1] * ePart
*/

/*
  Input
    ePart     - energy of incoming particle
    iso       - isotope
  Output
    getTotEngNEnergyCons - the average total energy available to the 
                           prompt fission gamma-rays
*/
   static int ZAEngConsN [nZAEngConsN]= {
                      92235, 92238, 
                      94239
                      };

   static double VogtNeutron [nZAEngConsN][2] = {
                           {4.838, 0.3004}, // 92235
                           {4.558, 0.3070}, // 92238
                           {6.128, 0.3428}  // 94239
                            };

   double muEn;

   int i;

// Find parameters for isotope
   int isoindex=-1;
   for (i=0; isoindex == -1 && i<nZAEngConsN; i++) {
      if (iso == ZAEngConsN[i]) {
         isoindex = i;
         break;
      }
   }
   if (isoindex == -1) {
      std::string errMsg = "No total fission neutron energy available for iso "; errMsg += iso;
      fissionerr(6, "getTotEngNEnergyCons", errMsg);
   }
   
// Computing average total energy for prompt fission neutrons
   muEn = VogtNeutron[isoindex][0] + ePart*VogtNeutron[isoindex][1];

   return muEn;
}
