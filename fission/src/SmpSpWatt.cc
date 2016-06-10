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


#define nZAspfis 19
#define WATTEMIN 1.0e-6
#define WATTEMAX 20.0

#include <string>
#include <math.h>
#include <stdio.h>
#include "fissionEvent.h"

double fissionEvent::SmpSpWatt(int iso) {

/*
  Description
    Sample Watt Spectrum for spontaneous fissions as in TART (Kalos algorithm) for 
        Th-232, 
        U-232, U-233, U-234, U-235, U-236, U-238
        Np-237, 
        Pu-238, Pu-239, Pu-240, Pu-241,  Pu-242
        Am-241, 
        Cm-242, Cm-244, 
        Bk-249,
        Cf-252
    based on Ensslin's data.
    N. Ensslin, et.al., "Application Guide to Neutron Multiplicity Counting," 
    LA-13422-M (November 1998) and
        Pu-236
    based on the assumption that Pu-236 has the same Watt spectrum as Np-237 
    since they have approximately the same nubar (2.05 versus 2.07). We think
    this assumption is valid since Dermott Cullen showed in "Sampling ENDL Watt 
    Fission Spectra," UCRL-TR-203251 that the Watt fission spectra can very well
    be approximated with a single parameter a by setting b equal to 1, instead 
    of 2 parameters a and b. He showed this for neutron-induced fissions. Since
    there is only 1 parameter characterizing a Watt spectrum, Watt spectra with
    identical nubars must have the same value for that parameter a (that is 
    because the integral of the spectrum with respect to the energy gives nubar, 
    normalization excluded). If we assume that Watt spectra can be approximated 
    by a single parameter a for spontaneous fissions as well (which I verified 
    and it seems to be a valid assumption), there can only be a single Watt 
    spectrum for a given nubar_sp. I thus concluded that the Watt spectrum for 
    Pu-236 is close to the Watt spectrum for Np-237.
*/

/*
  Input
    iso       - isotope
  Output
              - energy of incoming particle
*/

   static int nZAsp [nZAspfis]= {
                      90232,
                      92232, 92233, 92234, 92235, 92236, 92238,
                      93237,
                      94236, 94238, 94239, 94240, 94241, 94242,
                      95241,
                      96242, 96244,
                      97249,
                      98252
                      };

   static double WattSpa [nZAspfis] = {
                            1.25000e+00, // 90232
                            1.12082e+00, // 92232
                            1.16986e+00, // 92233
                            1.29661e+00, // 92234
                            1.29080e+00, // 92235
                            1.36024e+00, // 92236
                            1.54245e+00, // 92238
                            1.19985e+00, // 93237
                            1.19985e+00, // 94236
                            1.17948e+00, // 94238
                            1.12963e+00, // 94239
                            1.25797e+00, // 94240
                            1.18698e+00, // 94241
                            1.22078e+00, // 94242
                            1.07179e+00, // 95241
                            1.12695e+00, // 96242
                            1.10801e+00, // 96244
                            1.12198e+00, // 97249
                            8.47458e-01  // 98252
                            };

   static double WattSpb [nZAspfis] = {
                            4.00000e+00, // 90232
                            3.72278e+00, // 92232
                            4.03210e+00, // 92233
                            4.92449e+00, // 92234
                            4.85231e+00, // 92235
                            5.35746e+00, // 92236
                            6.81057e+00, // 92238
                            4.24147e+00, // 93237
                            4.24147e+00, // 94236
                            4.16933e+00, // 94238
                            3.80269e+00, // 94239
                            4.68927e+00, // 94240
                            4.15150e+00, // 94241
                            4.36668e+00, // 94242
                            3.46195e+00, // 95241
                            3.89176e+00, // 96242
                            3.72033e+00, // 96244
                            3.79405e+00, // 97249
                            1.03419e+00  // 98252
                            };

   double a;  /* Watt Parameters */
   double b;

   double rand1,rand2;
   double x,y,z;
   double eSmp;
   int i;


/*
   Find Watt parameters for isotope
*/
   int isoindex=-1;
   for (i=0; isoindex == -1 && i<nZAspfis; i++) {
      if (iso == nZAsp[i]) isoindex = i;
   }
   if (isoindex == -1) {
      std::string errMsg = "No Watt spectrum available for iso "; errMsg += iso;
      fissionerr(6, "SmpSpWatt", errMsg);
   }
   
   a= WattSpa[isoindex];
   b= WattSpb[isoindex];

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
