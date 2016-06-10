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

double fissionEvent::getTotEngP(double ePart, double nubar, int isotope, int A, int Z, bool spontaneous) {

/*
  Description
    Returns the average total energy available to prompt fission gamma-rays.
    Depending on the value of fissionEvent::effectivecorrelationoption,
    this energy will be computed 1 of 4 ways:
    switch (fissionEvent::effectivecorrelationoption)
       case(0): the average total energy is the one of Tim Valentine's 
                model described in T.E. Valentine, "Evaluation of Prompt 
                 Fission Gamma Rays for Use in Simulating Nuclear 
                 Safeguard Measurements," Ann. Nucl. Eng., 28, 191 
                (2001).
       case(1): the average total energy is determined by sampling a 
                normal distribution of mean <Etot_p> and standard 
                deviation <Etot_p>/8, where <Etot_p> is given in B. Beck, 
                D.A. Brown, F. Daffin, J. Hedstrom, R. Vogt, 
                "Implementation of Energy-Dependent Q Values for 
                 Fission," UCRL-TR-234617, Lawrence Livermore National 
                Laboratory (2007).
       case(2): the average total energy is determined by sampling a 
                normal distribution of mean <Etot_p> and standard 
                deviation <Etot_p>/8, where <Etot_p> is given in R. Vogt, 
                "Energy-Dependent Fission Q Values Generalized for All 
                 Actinides," LLNL-TR-407620, Lawrence Livermore 
                National Laboratory (2008).
       case(3): the average total energy is irrelevant and set to 0, 
                energies are determined independently by FREYA in this 
                case.
    For spontaneous fissions, the average total energy available to 
    prompt fission gamma-rays is computed via T.E. Valentine's method.
*/

/*
  Input
    ePart       - energy of incoming neutron causing fission. If it is 
                  equal to -1, it means this is a spontaneous fission
    nubar       - average number of fission neutrons
    isotope     - isotope
    A           - number of nuclei in atom, if set to -1, A and Z need
                  to be computed
    Z           - number of protons in atom, if set to -1, Z and Z need
                  to be computed
    spontaneous - set to true for spontaneous fissions

  Output
    getTotEngP - average total energy available to prompt fission 
                 gamma-rays
*/
 
  double energy; // average total energy of the prompt fission gamma-rays

  if (A == -1 || Z == -1) {
     A = (int) (isotope-1000*((int)(isotope/1000)));
     Z = (int) ((isotope-A)/1000);
  }

  // We first treat the spontaneous fission case
  if (spontaneous) {
     energy = (2.51-1.13e-5*pow((double)Z,2.)*sqrt((double)A))*nubar+4.0;
     return energy;
  }

  switch(effectivecorrelationoption) {
     case(0): 
        energy = (2.51-1.13e-5*pow((double)Z,2.)*sqrt((double)A))*nubar+4.0;
        break;
     case(1): 
        energy = getTotEngPEnergyCons(ePart, isotope);
        break;
     case(2): 
        energy = getTotEngPEnergyConsAllActinides(ePart, isotope);
        break;
     case(3): 
        energy = 0;
        break;
     default: 
        energy = 0;
        break;
  }
  return energy;
}
