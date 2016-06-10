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

#define nSPfissg 30

int fissionEvent::SmpSpNugDistData(int isotope) {

/*
  Description
    Sample Number of Photons from spontaneous fission in 
    (a) Cf-252 using the double Poisson model from Brunson;
    (b) Th-232, 
        U-232, U-233, U-234, U-235, U-236, U-238*,
        Np-237, 
        Pu-236**, Pu-238*, Pu-239, Pu-240*, Pu-241, Pu-242*,
        Am-241, 
        Cm-242*, Cm-244*, Cm-246***, Cm-248***,
        Bk-249,
        Cf-246****, Cf-250*****, Cf-254*****,
        Fm-257*****,
        No-252******
        using the negative binomial distribution based on the
        spontaneous fission neutron nubar from Ensslin's 
        tabulated data, Holden and Zucker's tabulated data 
        (for isotopes denoted with asterix *), P.Santi, 
        D.H. Beddingfield and D.R. Mayo tabulated data
        (for isotopes denoted with a double asterix **),
        BNL-36467's tabulated data (for isotopes denoted 
        with ***), Dakavoski's tabulated data (for isotopes 
        denoted with ****), Hoffman's tabulated data (for 
        isotopes denoted with *****), Lazarev's tabulated
        data (for isotopes denoted with ******).
*/

/*
  Input
    isotope          - isotope

  Output
    SmpSpNugDistData - sampled multiplicity
                       -1 if there is no multiplicity data for that isotope
*/
 
  int i;
  double sum, nubar;
  double r;

  static double Cf252spdist [nSPfissg] = { 
         5.162699e-4,3.742057e-3,1.360482e-2,3.312786e-2,6.090540e-2,
         9.043537e-2,1.133984e-1,1.240985e-1,1.216759e-1,1.092255e-1,
         9.137106e-2,7.219960e-2,5.438060e-2,3.923091e-2,2.714690e-2,
         1.800781e-2,1.143520e-2,6.942099e-3,4.025720e-3,2.229510e-3,
         1.179602e-3,5.966936e-4,2.888766e-4,1.340137e-4,5.965291e-5,
         2.551191e-5,1.049692e-5,4.160575e-6,1.590596e-6,0.000000e+0
      };

/*
  sample the spontaneous fission photon number distribution 
*/
  nubar=0.;
//  Cf-252 using the double Poisson model from Brunson;
  if (isotope == 98252) {
    r=fisslibrng();

    sum=0.;
    for(i=0; i<nSPfissg; i++) {
       sum=sum+Cf252spdist[i];
       if (r <= sum || Cf252spdist[i+1] == 0.) return i;
    }
/*
    using the spontaneous fission nubar from
    Holden and Zucker's tabulated data.
*/
  } else if (isotope == 92238) {
    nubar = 1.9900002;
  } else if (isotope == 94240) {
    nubar = 2.1540006;
  } else if (isotope == 94242) {
    nubar = 2.1489998;
  } else if (isotope == 96242) {
    nubar = 2.54;
  } else if (isotope == 96244) {
    nubar = 2.7200005;
/*
    using the spontaneous fission nubar from
    P.Santi, D.H. Beddingfield and D.R. Mayo 
    tabulated data.
*/
  } else if (isotope == 94236) {
    nubar = 2.07;
  } else if (isotope == 94238) {
    nubar = 2.19;
/*
    using the spontaneous fission nubar from
    BNL-36467's tabulated data.
*/
  } else if (isotope == 96246) {
    nubar = 2.93;
  } else if (isotope == 96248) {
    nubar = 3.13;
/*
    using the spontaneous fission nubar from
    Dakavoski's tabulated data.
*/
  } else if (isotope == 98246) {
    nubar = 3.1;
/*
    using the spontaneous fission nubar from
    Hoffman's tabulated data.
*/
  } else if (isotope == 98250) {
    nubar = 3.51;
  } else if (isotope == 98254) {
    nubar = 3.85;
  } else if (isotope == 100257) {
    nubar = 3.87;
/*
    using the spontaneous fission nubar from
    Lazarev's tabulated data.
*/
  } else if (isotope == 102252) {
    nubar = 4.2;
  }

  if (nubar != 0.) {
    return SmpNugDist(isotope, nubar, 0., true);
  } else {
/*
    using the spontaneous fission nubar from
    N. Ensslin, et.al., "Application Guide to Neutron
    Multiplicity Counting," LA-13422-M (November 1998)
*/
    nubar = SmpSpNubarData(isotope);
    if (nubar != -1.) {
      return SmpNugDist(isotope, nubar, 0., true);
    } else {
// There is no nubar information for that isotope, return -1,
// meaning no data available for that isotope
      return -1;
    }
  }
}
