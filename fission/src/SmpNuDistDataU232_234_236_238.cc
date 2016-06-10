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


#include <math.h>
#include "fissionEvent.h"

int fissionEvent::SmpNuDistDataU232_234_236_238(double nubar) {

/*
  Description
    Sample Number of Neutrons from fission in U-232, U-234, U-236 
    and U-238 using Zucker and Holden's tabulated data for U-238
    The P(nu) distribution is given as a function of the average 
    number of neutrons from fission, based on interpolation of the 
    U-238 data from Zucker and Holden.
*/

/*
  Input
    nubar    - average number of neutrons per fission
  Output
    SmpNuDistDataU232_234_236_238  - sampled multiplicity
    
*/

  double pnu[9], cpnu, sum;
  double r;

/* 
  Check if nubar is within the range of experimental values
*/
  if (nubar >= 2.25 && nubar <= 3.80) {
/*
     Use Zucker and Holden Data
*/
     pnu[0]=-7.705432e-3*pow(nubar,3)+8.904671e-2*pow(nubar,2)-3.488123e-1*nubar+4.627291e-1;
     pnu[1]=-2.879938e-2*pow(nubar,3)+3.629189e-1*pow(nubar,2)-1.545284*nubar+2.229503;
     pnu[2]=6.543684e-2*pow(nubar,3)-6.673117e-1*pow(nubar,2)+2.087358*nubar-1.771396;
     pnu[3]=1.412971e-2*pow(nubar,3)-2.309842e-1*pow(nubar,2)+1.022451*nubar-1.032235;
     pnu[4]=-5.163167e-2*pow(nubar,3)+4.457516e-1*pow(nubar,2)-1.114981*nubar+9.484241e-1;
     pnu[5]=8.758841e-4*pow(nubar,3)+3.707461e-2*pow(nubar,2)-1.565149e-1*nubar+1.851039e-1;
     pnu[6]=-3.871089e-5*pow(nubar,3)+1.936524e-2*pow(nubar,2)-8.091057e-2*nubar+9.019871e-2;
     pnu[7]=3.945995e-3*pow(nubar,3)-2.697509e-2*pow(nubar,2)+6.237296e-2*nubar-4.820745e-2;
     pnu[8]=1.708054e-3*pow(nubar,4)-1.706039e-2*pow(nubar,3)+6.550213e-2*pow(nubar,2)-1.135e-1*nubar+7.443828e-2;

     sum=pnu[0]+pnu[1]+pnu[2]+pnu[3]+pnu[4]+pnu[5]+pnu[6]+pnu[7]+pnu[8];

     pnu[0]=pnu[0]/sum;
     pnu[1]=pnu[1]/sum;
     pnu[2]=pnu[2]/sum;
     pnu[3]=pnu[3]/sum;
     pnu[4]=pnu[4]/sum;
     pnu[5]=pnu[5]/sum;
     pnu[6]=pnu[6]/sum;
     pnu[7]=pnu[7]/sum;
     pnu[8]=pnu[8]/sum;

     r=fisslibrng();

     if(r <= pnu[0]) return 0;

     cpnu=pnu[0]+pnu[1];
     if(r <= cpnu) return 1;

     cpnu=cpnu+pnu[2];
     if(r <= cpnu) return 2;

     cpnu=cpnu+pnu[3];
     if(r <= cpnu) return 3;

     cpnu=cpnu+pnu[4];
     if(r <= cpnu) return 4;

     cpnu=cpnu+pnu[5];
     if(r <= cpnu) return 5;

     cpnu=cpnu+pnu[6];
     if(r <= cpnu) return 6;

     cpnu=cpnu+pnu[7];
     if(r <= cpnu) return 7;
     else return 8;

  } else {
/*
     Use Terrell's formula
*/
     return (int) SmpTerrell(nubar);
  }
}
