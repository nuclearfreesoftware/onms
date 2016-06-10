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


#include <string>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include "fissionEvent.h"

void fissionEvent::SmpNPEnergyCons(double ePart, double nubar, int iso, bool spontaneous) {

/*
  Description
    Determine the energy of prompt neutrons and photons emitted by fission.
    For neutrons, sample Watt Spectrum as in TART (Kalos algorithm) 
    under the constraint that the total emitted prompt neutron energy is bound.
    For photons, sample photon energy spectrum from Maienschein's measurements
    under the constraint that the total emitted prompt photon energy is bound.
*/

/*
  Input
    ePart       - energy of incoming particle
    nubar       - average number of fission neutrons
    iso         - isotope
    spontaneous -  set to true for spontaneous fissions
  Output
    neutronEnergies (implicitly)
              - energies of the emitted fission neutrons
    photonEnergies (implicitly)
              - energies of the emitted fission photons
*/

   bool scaling=true;

   double Etotn,Etotp;
   double muEn, muEp;
   double sigman, sigmap;

   int i;

// Computing total energies for neutrons and photons separately
   muEn = getTotEngN(ePart, nubar, iso, spontaneous);
   sigman = muEn*.25;
   do {
      Etotn = muEn+sigman*normsinv(fisslibrng());
//   } while (Etotn < muEn/32); // We do not want too low an energy
   } while (Etotn < 0.01); // We do not want too low an energy

   muEp = getTotEngP(ePart, nubar, iso, -1, -1, spontaneous);
   sigmap = muEp*.125;
   do {
      Etotp = muEp+sigmap*normsinv(fisslibrng());
//   } while (Etotp < muEp/2); // We do not want too low an energy
   } while (Etotp < 0.100); // We do not want too low an energy

   if (!scaling) {
//    Sampling the Watt spectrum neutronNu-1 times and then add a neutron with 
//    an energy such that the total of the neutron energies adds up to Etotn
      double sumNEng = 0.;
      int resample = 0;
      for (i=0; i<neutronNu-1; i++) {
         neutronEnergies[i] = SmpWatt(ePart, iso);
         sumNEng += neutronEnergies[i];
         if (sumNEng >= Etotn) { // We already have too much energy, resample a neutron energy
            sumNEng -= neutronEnergies[i];
            resample++;
            if (resample > 10) {
               // try resampling 10 times, otherwise, split the energy among the left over neutrons
               double remainder = (Etotn-sumNEng)/(neutronNu-i);
               sigman = .25*remainder;
               for (int j=i; j<neutronNu; j++) {
                  do {
                     neutronEnergies[j] = remainder+sigman*normsinv(fisslibrng());;
                  } while (neutronEnergies[j] < 0);
                  sumNEng += neutronEnergies[j];
               }
               break;
            }
            i--;
         }
      }
      if (neutronNu != 0 && resample <= 10) {
         neutronEnergies[neutronNu-1] = Etotn-sumNEng;
         sumNEng += neutronEnergies[neutronNu-1];
      }
   } else { // scaling method
// We sample the Watt spectrum neutronNu times and then scale the neutron
// energies so that the total of the neutron energies adds up to Etotn
      double sumNEng = 0.;
      for (i=0; i<neutronNu; i++) {
         neutronEnergies[i] = SmpWatt(ePart, iso);
         sumNEng += neutronEnergies[i];
      }
      double ratio = Etotn/sumNEng;
      for (i=0; i<neutronNu; i++) neutronEnergies[i] *= ratio;
   }

   if (!scaling) {
//    Sampling the fission gamma-ray spectrum photonNu-1 times and then add a gamma-ray with 
//    an energy such that the gamma-ray energies add up to Etotp
      double sumPEng = 0.;
      int resample = 0;
      for (i=0; i<photonNu-1; i++) {
         photonEnergies[i] = SmpGEng();
         sumPEng += photonEnergies[i];
         if (sumPEng >= Etotp) { // We already have too much energy, resample a gamma-ray energy
            sumPEng -= photonEnergies[i];
            resample++;
            if (resample > 10) {
               // try resampling 10 times, otherwise, split the energy among the left over gamma-ray
               double remainder = (Etotp-sumPEng)/(photonNu-i);
               sigmap = .25*remainder;
               for (int j=i; j<photonNu; j++) {
                  do {
                     photonEnergies[j] = remainder+sigmap*normsinv(fisslibrng());;
                  } while (photonEnergies[j] < 0);
                  sumPEng += photonEnergies[j];
               }
               break;
            }
            i--;
         }
      }
      if (photonNu != 0 && resample <= 10) {
         photonEnergies[photonNu-1] = Etotp-sumPEng;
         sumPEng += photonEnergies[photonNu-1];
      }
   } else { // scaling method
// We sample the fission gamma-ray spectrum photonNu times and then 
// scale the gamma-ray energies so that the total of the gamma-ray 
// energies adds up to Etotp
      double sumPEng = 0.;
      for (i=0; i<photonNu; i++) {
         photonEnergies[i] = SmpGEng();
         sumPEng += photonEnergies[i];
      }
      double ratio = Etotp/sumPEng;
      for (i=0; i<photonNu; i++) photonEnergies[i] *= ratio;
   }

   return;
}
