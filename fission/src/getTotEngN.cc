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

double fissionEvent::getTotEngN(double ePart, double nubar, int isotope, bool spontaneous) {

/*
  Description
    Returns the average total energy available to prompt fission neutrons.
    Depending on the value of fissionEvent::effectivecorrelationoption,
    this energy will be computed 1 of 3 ways:
    switch (fissionEvent::effectivecorrelationoption)
       case(1): the average total energy is determined by sampling a 
                normal distribution of mean <Etot_n> and standard 
                deviation <Etot_n>/4, where <Etot_n> is given in B. Beck, 
                D.A. Brown, F. Daffin, J. Hedstrom, R. Vogt, 
                "Implementation of Energy-Dependent Q Values for 
                 Fission," UCRL-TR-234617, Lawrence Livermore National 
                Laboratory (2007).
       case(2): the average total energy is determined by sampling a 
                normal distribution of mean <Etot_n> and standard 
                deviation <Etot_n>/4, where <Etot_n> is given in R. Vogt, 
                "Energy-Dependent Fission Q Values Generalized for All 
                 Actinides," LLNL-TR-407620, Lawrence Livermore 
                National Laboratory (2008).
       default: the average total energy is irrelevant and set to 0.
    For spontaneous fissions, the average total energy available to 
    prompt fission neutrons is computed via T.E. Valentine's method.
*/

/*
  Input
    ePart       - energy of incoming neutron causing fission
    nubar       - average number of fission neutrons
    isotope     - isotope
    spontaneous -  set to true for spontaneous fissions

  Output
    getTotEngN - average total energy available to prompt fission 
                 gamma-rays
*/
 
  double energy; // average total energy of the prompt fission gamma-rays

  // We first treat the spontaneous fission case
  if (spontaneous) {
     energy = getTotEngNEnergyConsAllActinides(0., isotope);
     return energy;
  }

  switch(effectivecorrelationoption) {
     case(1): 
        energy = getTotEngNEnergyCons(ePart, isotope);
        break;
     case(2): 
        energy = getTotEngNEnergyConsAllActinides(ePart, isotope);
        break;
     default: 
        energy = 0;
        break;
  }
  return energy;
}
