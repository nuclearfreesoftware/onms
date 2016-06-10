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


#ifndef FISSIONEVENT_H
#define FISSIONEVENT_H

#include <string>
#include <vector>
#include <stdlib.h>

class fissionEvent {
   private:
      int neutronNu; // number of neutrons in this fission event
      double* neutronEnergies; 
      double* neutronVelocities; 
      double* neutronDircosu; 
      double* neutronDircosv; 
      double* neutronDircosw; 
      double* neutronAges; 

      int photonNu; // number of photons in this fission event
      double* photonEnergies; 
      double* photonVelocities; 
      double* photonDircosu; 
      double* photonDircosv; 
      double* photonDircosw; 
      double* photonAges; 

      // options
      static int delayoption;
      static int correlationoption;
      int effectivecorrelationoption;
      static int nudistoption;
      static int Cf252ndistoption;
      static int Cf252nengoption;
      static double (*rngdptr)(void);
      static float (*rngfptr)(void);

   public:
      // These are all the methods of this class accessible to the caller of the object 
      fissionEvent(int isotope, double time, double nubar, double eng, int fissiontype);
      ~fissionEvent();
      int getNeutronNu() {
         return neutronNu;
      }
      int getPhotonNu() {
         return photonNu;
      }
      double getNeutronEnergy(int index) {
         if (index >= 0 && index < neutronNu) return neutronEnergies[index];
         else return -1;
      }
      double getNeutronVelocity(int index) {
         if (index >= 0 && index < neutronNu) return neutronVelocities[index];
         else return -1;
      }
      double getNeutronDircosu(int index) {
         if (index >= 0 && index < neutronNu) return neutronDircosu[index];
         else return -1;
      }
      double getNeutronDircosv(int index) {
         if (index >= 0 && index < neutronNu) return neutronDircosv[index];
         else return -1;
      }
      double getNeutronDircosw(int index) {
         if (index >= 0 && index < neutronNu) return neutronDircosw[index];
         else return -1;
      }
      double getPhotonEnergy(int index) {
         if (index >= 0 && index < photonNu) return photonEnergies[index];
         else return -1;
      }
      double getPhotonVelocity(int index) {
         if (index >= 0 && index < photonNu) return photonVelocities[index];
         else return -1;
      }
      double getPhotonDircosu(int index) {
         if (index >= 0 && index < photonNu) return photonDircosu[index];
         else return -1;
      }
      double getPhotonDircosv(int index) {
         if (index >= 0 && index < photonNu) return photonDircosv[index];
         else return -1;
      }
      double getPhotonDircosw(int index) {
         if (index >= 0 && index < photonNu) return photonDircosw[index];
         else return -1;
      }
      double getNeutronAge(int index) {
         if (index >= 0 && index < neutronNu) return neutronAges[index];
         else return -1;
      }
      double getPhotonAge(int index) {
         if (index >= 0 && index < photonNu) return photonAges[index];
         else return -1;
      }
      void getFREYAerrors(int* length, char* errors);
      static void setDelayOption(int delay) {
         delayoption = delay;
      };
      static void setCorrelationOption(int correlation) {
         correlationoption = correlation;
      };
      static int getCorrelationOption() {
         return correlationoption;
      };
      static void setNudistOption(int nudist) {
         nudistoption = nudist;
      };
      static int getNudistOption() {
         return nudistoption;
      };
      static void setCf252Option(int ndist, int neng) {
         Cf252ndistoption = ndist;
         Cf252nengoption = neng;
      };
      static int getCf252ndist() {
         return Cf252ndistoption;
      };
      static int getCf252neng() {
         return Cf252nengoption;
      };
      static void setRNGf(float (*funcptr) (void)) {
         rngfptr = funcptr;
         rngdptr = rngf2d;
      }
      static void setRNGd(double (*funcptr) (void)) {
         rngdptr = funcptr;
      }
      static double erf(double x);
      static double erfc(double x);
      static double getNSepEng(int iso);
      inline static double fisslibrng(void) {
         if (rngdptr == 0) {
#ifdef LLNLBUILD
            return drand48();
#endif
         } else {
            return rngdptr();
         }
         return 0.;
      }

   private:
      int SmpNuDistDataU232_234_236_238(double nubar);
      int SmpNuDistDataU232_234_236_238_MC(double nubar);
      int SmpNuDistDataU233_235(double nubar);
      int SmpNuDistDataU233_235_MC(double nubar);
      int SmpNuDistDataU235(double erg, int option);
      int SmpNuDistDataPu239(double erg);
      double SmpNVel(double eng);
      double SmpNEngCf252(int option);
      void SmpIsoDir(double* cosdiru, double* cosdirv, double* cosdirw);
      double SmpGEng();
      int SmpNuDistDataPu239_241(double nubar);
      int SmpNuDistDataPu239_241_MC(double nubar);
      int SmpNuDistDataU238(double erg);
      int SmpNugDist(int isotope, double nubar, double ePart, bool spontaneous);
      double getNubarg(int isotope, double nubar, double ePart, bool spontaneous);
      double SmpPVel();
      int SmpSpNuDistData(int isotope, int Cf252option);
      double SmpSpNubarData(int isotope);
      int SmpSpNugDistData(int isotope);
      double SmpTerrell(double nubar);
      double SmpSpWatt(int iso);
      double SmpWatt(double ePart, int iso);
      void SmpNPEnergyCons(double ePart, double nubar, int iso, bool spontaneous);
      double getTotEngN(double ePart, double nubar, int isotope, bool spontaneous);
      double getTotEngP(double ePart, double nubar, int isotope, int A, int Z, bool spontaneous);
      double getTotEngNEnergyCons(double ePart, int isotope);
      double getTotEngPEnergyCons(double ePart, int isotope);
      double getTotEngNEnergyConsAllActinides(double ePart, int isotope);
      double getTotEngPEnergyConsAllActinides(double ePart, int isotope);
#ifdef FREYA
      bool SmpFreya(double ePart, int iso, int fissiontype);
      void saveFREYANeutron(float* pp);
      void saveFREYAPhoton(float* pp);
#endif
      static void fissionerr(int iSever, std::string chSubNam, std::string chMsg);
      void allocateMem(int nu, int nug);
      double normsinv(double p);
      static double rngf2d(void);
#ifdef FREYA
      bool handle_freya_error();
#endif

      static double erf1(double x);
      static double erfc1(double x);
      static double erf2(double x);
      static double erfc2(double x);
      static double erf3(double x);
      static double erfc3(double x);
      static double erf4(double x);
      static double erfc4(double x);
      static double erf5(double x);
      static double erfc5( double x);
};
#endif
