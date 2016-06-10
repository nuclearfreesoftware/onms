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

#define nSPfissNubarIso 18

double fissionEvent::SmpSpNubarData(int isotope) {

/*
  Description
    Determine average number of neutrons from spontaneous fission for
        Th-232, 
        U-232, U-233, U-234, U-235, U-236, U-238
        Np-237, 
        Pu-239, Pu-240, Pu-241,  Pu-242
        Am-241, 
        Cm-242, Cm-244, 
        Bk-249,
        Cf-252
    Based on Ensslin's data.
    N. Ensslin, et.al., "Application Guide to Neutron Multiplicity Counting," 
    LA-13422-M (November 1998)
*/

/*
  Input
    iso          - isotope
  Output
    SmpSpNubarData - average number of neutrons
                     -1. is the isotope has 
                         no nubar data
*/
 
  int i;

  static int spzaid [nSPfissNubarIso] = {
      90232, 92232, 92233, 92234, 92235,
      92236, 92238, 93237, 94238, 94239, 
      94240, 94241, 94242, 95241, 96242, 
      96244, 97249, 98252 };
  static double spnubar [nSPfissNubarIso] = {
      2.14,  1.71, 1.76,  1.81, 1.86,
      1.91,  2.01, 2.05,  2.21, 2.16, 
      2.156, 2.25, 2.145, 3.22, 2.54, 
      2.72,  3.40, 3.757
      };

// Find nubar
  for (i=0; i<nSPfissNubarIso; i++) {
    if (isotope == spzaid[i]) {
      return spnubar[i];
    }
  }
// no nubar available for that isotope
  return -1.;
}
