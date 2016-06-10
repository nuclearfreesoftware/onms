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

int fissionEvent::SmpNuDistDataU233_235(double nubar) {

/*
  Description
    Sample Number of Neutrons from fission in U-233 and U-235 using 
    (a) Gwin, Spencer and Ingle tabulated data at thermal 
        energies (0 MeV),
    (b) Zucker and Holden's tabulated data for U-235 at 1 MeV and 
        higher.
    The P(nu) distribution is given as a function of the average 
    number of neutrons from fission, based on interpolation of the 
    U-235 data above.
*/

/*
  Input
    nubar    - average number of neutrons per fission
  Output
    SmpNuDistDataU233_235  - sampled multiplicity
    
*/

  double pnu[8], cpnu, sum;
  double r;

/* 
  Check if nubar is within the range of experimental values
*/
  if(nubar >= 2.25 && nubar <= 4.0) {
/*
     Use Zucker and Holden Data
*/
     if(nubar <= 2.8738) pnu[0]=-9.279554e-02*pow(nubar,3)+8.036687e-01*pow(nubar,2)-2.342684*nubar+2.309035;
     else if(nubar > 2.8738 && nubar <= 3.4272) pnu[0]=1.50072e-2*pow(nubar,2)-1.109109e-1*nubar+2.063133e-1;
     else pnu[0]=1.498897e+3*exp(-3.883864*nubar);

     if(nubar <= 3.2316) pnu[1]=3.531126e-2*pow(nubar,3)-2.787213e-1*pow(nubar,2)+5.824072e-1*nubar-1.067136e-1;
     else pnu[1]=6.574492e-2*pow(nubar,2)-5.425741e-1*nubar+1.123199;

     pnu[2]=1.274643e-2*pow(nubar,3)-1.387954e-1*pow(nubar,2)+3.264669e-1*nubar+1.77148e-1;

     pnu[3]=5.473738e-2*pow(nubar,5)-8.835826e-1*pow(nubar,4)+5.657201*pow(nubar,3)-1.802669e+1*pow(nubar,2)+2.867937e+1*nubar-1.794296e+1;

     pnu[4]=-3.591076e-2*pow(nubar,3)+3.092624e-1*pow(nubar,2)-7.184805e-1*nubar+5.649400e-1;

     if(nubar <= 2.8738) pnu[5]=1.699374e-2*pow(nubar,2)-1.069558e-3*nubar-6.981430e-2;
     else pnu[5]=2.100175e-2*pow(nubar,3)-1.705788e-1*pow(nubar,2)+5.575467e-1*nubar-6.245873e-1;

     if(nubar <= 3.0387) pnu[6]=9.431919e-7*pow(nubar,8.958848);
     else pnu[6]=4.322428e-3*pow(nubar,3)-2.094790e-2*pow(nubar,2)+4.449671e-2*nubar-4.435987e-2;

     pnu[7]=5.689084e-3*pow(nubar,4)-6.591895e-2*pow(nubar,3)+2.886861e-1*pow(nubar,2)-5.588146e-1*nubar+4.009166e-1;

     sum=pnu[0]+pnu[1]+pnu[2]+pnu[3]+pnu[4]+pnu[5]+pnu[6]+pnu[7];

     pnu[0]=pnu[0]/sum;
     pnu[1]=pnu[1]/sum;
     pnu[2]=pnu[2]/sum;
     pnu[3]=pnu[3]/sum;
     pnu[4]=pnu[4]/sum;
     pnu[5]=pnu[5]/sum;
     pnu[6]=pnu[6]/sum;
     pnu[7]=pnu[7]/sum;

     r=fisslibrng();

     if(r <= pnu[0]) return (int) 0;

     cpnu=pnu[0]+pnu[1];
     if(r <= cpnu) return (int) 1;

     cpnu=cpnu+pnu[2];
     if(r <= cpnu) return (int) 2;

     cpnu=cpnu+pnu[3];
     if(r <= cpnu) return (int) 3;

     cpnu=cpnu+pnu[4];
     if(r <= cpnu) return (int) 4;

     cpnu=cpnu+pnu[5];
     if(r <= cpnu) return (int) 5;

     cpnu=cpnu+pnu[6];
     if(r <= cpnu) return (int) 6;
     else return (int) 7;
  } else {
/*
     Use Terrell's formula
*/
     return (int) SmpTerrell(nubar);
  }
}
