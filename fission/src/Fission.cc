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

fissionEvent* fe;

extern "C" {
   double rngdptr_() {
      return fissionEvent::fisslibrng();
   };

   void genspfissevt_(int *isotope, double *time) {
      if (fe != 0) delete fe;
      fe = new fissionEvent(*isotope, *time, -1., 0., 0);
   };

   void genfissevt_(int *isotope, double *time, double *nubar, double *eng) {
      if (fe != 0) delete fe;
      fe = new fissionEvent(*isotope, *time, *nubar, *eng, 1);
   };

   void genphotofissevt_(int *isotope, double *time, double *nubar, double *eng) {
      if (fe != 0) delete fe;
      // The operations below are now performed in fissionEvent constructor
      // int fissionIsotope = *isotope-1;
      // double neutronEnergy = *eng-fissionEvent::getNSepEng(*isotope);
      // if (neutronEnergy<0) neutronEnergy=0.;
      // fe = new fissionEvent(fissionIsotope, *time, *nubar, neutronEnergy, 2);
      fe = new fissionEvent(*isotope, *time, *nubar, *eng, 2);
   };

   int getnnu_() {
      return (*fe).getNeutronNu();
   };

   int getpnu_() {
      return (*fe).getPhotonNu();
   };

   double getneng_(int *index) {
      return (*fe).getNeutronEnergy(*index);
   };

   double getnvel_(int *index) {
      return (*fe).getNeutronVelocity(*index);
   };

   double getndircosu_(int *index) {
      return (*fe).getNeutronDircosu(*index);
   };

   double getndircosv_(int *index) {
      return (*fe).getNeutronDircosv(*index);
   };

   double getndircosw_(int *index) {
      return (*fe).getNeutronDircosw(*index);
   };

   double getpeng_(int *index) {
      return (*fe).getPhotonEnergy(*index);
   };

   double getpvel_(int *index) {
      return (*fe).getPhotonVelocity(*index);
   };

   double getpdircosu_(int *index) {
      return (*fe).getPhotonDircosu(*index);
   };

   double getpdircosv_(int *index) {
      return (*fe).getPhotonDircosv(*index);
   };

   double getpdircosw_(int *index) {
      return (*fe).getPhotonDircosw(*index);
   };

   double getnage_(int *index) {
      return (*fe).getNeutronAge(*index);
   };

   double getpage_(int *index) {
      return (*fe).getPhotonAge(*index);
   };

   void setdelay_(int *delay) {
      (*fe).setDelayOption(*delay);
   };

   void setcorrel_(int *correlation) {
      (*fe).setCorrelationOption(*correlation);
   };

   void setnudist_(int *nudist) {
/*
      where the argument *nudist affects induced fissions only, it
      is set to
         0 for sampling Zucker and Holden probability distributions 
           for U-235,238 and Pu-239. Terrell for other isotopes.
         1 same as above, but using Gwin, Spencer and Ingle 
           tabulated distributions for thermal energies (0 MeV) 
           for U-235. Terrell for other isotopes.
         2 for sampling fission-induced neutron multiplicity in 
           (a) U-232, U-234, U-236 and U-238 using Zucker and 
               Holden's tabulated data for U-238
           (b) U-233 and U-235 using 
               (1) Gwin, Spencer and Ingle tabulated data at 0 MeV
               (2) Zucker and Holden's tabulated data at 1 MeV and 
                   above
               for U-235
           (c) Pu-239 and Pu-241 using Zucker and Holden's tabulated 
               data for Pu-239
           The P(nu) distributions for *nudist=2 are given as a 
           function of the average number of neutrons from fission, 
           based on interpolation of the data above.
           Terrell for other isotopes.
         3 for sampling fission-induced neutron multiplicity in 
           (a) U-232, U-234, U-236 and U-238 using Zucker and 
               Holden's tabulated data for U-238
           (b) U-233 and U-235 using 
               (1) Gwin, Spencer and Ingle tabulated data at 0 MeV
               (2) Zucker and Holden's tabulated data at 1 MeV and 
                   above
               for U-235
           (c) Pu-239 and Pu-241 using Zucker and Holden's tabulated 
               data for Pu-239
           The tables have P(nu) distributions for 11 energies 
           (0 MeV through 10 MeV), along with their nubars. For 
           *nudist=3, we select the P(nu) distribution that has
           a nubar closest either from above, or from below, to the 
           to the nubar entered for the induced fission, based on a 
           random number and fractional distances to the end of the 
           nubar interval thus formed.
           Terrell for other isotopes.
*/

      (*fe).setNudistOption(*nudist);
   };

   void setcf252_(int *ndist, int *neng) {
/*
      where the argument
      *ndist is set to 
         0 to sample the spontaneous fission neutron multiplicity 
           using tabulated data from Spencer
         1 to sample the spontaneous fission neutron multiplicity
           using tabulated data from Boldeman
      *neng is set to 
         0 to sample the Mannhart corrected Maxwellian spectrum
         1 to sample the Madland-Nix theoretical spectrum
         2 to sample the Froehner Watt spectrum
*/
      (*fe).setCf252Option(*ndist, *neng);
   }; 

   double getnsepeng_(int *nuclide) {
      return fissionEvent::getNSepEng(*nuclide);
   }

   void setrngf_(float (*funcptr) (void)) {
      fissionEvent::setRNGf(funcptr);
   }

   void setrngd_(double (*funcptr) (void)) {
      fissionEvent::setRNGd(funcptr);
   }

   void getfreya_errors_(int* length, char* error) {
#ifdef FREYA
      return (*fe).getFREYAerrors(length, error);
#endif
   }
}
