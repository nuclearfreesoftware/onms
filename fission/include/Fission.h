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


#ifdef __cplusplus
//....This is needed for C++, otherwise this header works only for C
extern "C" {
#endif
extern void genspfissevt_(int *isotope, double *time);
/*
 * This function is called to trigger a spontaneous fission.
 * Multiple neutrons and photons are generated and stored
 * in a stack along with their energies, directions and 
 * emission times.
 * The arguments of this function are
 *      isotope:        94239 for Pu-239 for instance
 *      time:           the time of the spontaneous fission 
 *                      in seconds
 */

extern void genfissevt_(int *isotope, double *time, double *nubar, double *eng);
/*
 * This function is called to trigger a neutron-induced fission.
 * Multiple neutrons and photons are generated and stored
 * in a stack along with their energies, directions and 
 * emission times. In addition to the arguments above, this
 * function needs
 *      nubar:          user-specified average number of neutrons emitted 
 *                      per fission (e.g. as tabulated in the cross-section 
 *                      libraries used by the particle transport code)
 *      eng:            energy of the neutron inducing fission
 */

extern void genphotofissevt_(int *isotope, double *time, double *nubar, double *eng);
/*
 * This function is called to trigger a photofission.
 * Multiple neutrons and photons are generated and stored
 * in a stack along with their energies, directions and 
 * emission times. The arguments of this function are
 *      isotope:        94239 for Pu-239 for instance
 *      time:           the time of the spontaneous fission in seconds
 *      nubar:          user-specified average number of neutrons emitted 
 *                      per photofission (e.g. as tabulated in the 
 *                      cross-section libraries used by the particle 
 *                      transport code)
 *      eng:            energy of the photon inducing fission
 */

extern int getnnu_();
/*
 * This function returns the number of neutrons emitted by the
 * fission, -1 if there is no neutron data for that isotope in 
 * the fission library.
 */

extern int getpnu_();
/*
 * This function returns the number of photons emitted by the
 * fission, -1 if there is no photon data for that isotope in 
 * the fission library.
 */

extern double getneng_(int *index);
/*
 * Given the index of the emitted neutron, this function returns
 * its energy, -1 if index isout of range.
 */

extern double getnvel_(int *index);
/*
 * Given the index of the emitted neutron, this function returns
 * the amplitude of its velocity, -1 if index is out of range.
 */

extern double getndircosu_(int *index);
extern double getndircosv_(int *index);
extern double getndircosw_(int *index);
/*
 * Given the index of the emitted neutron, this function returns
 * the direction cosines of its velocity vector on the x, y and z 
 * axes.
 */

extern double getpeng_(int *index);
/*
 * Given the index of the emitted photon, this function returns
 * its energy, -1 if index is out of range.
 */

extern double getpvel_(int *index);
/*
 * Given the index of the emitted photon, this function returns
 * the amplitude of its velocity, -1 if index is out of range.
 */

extern double getpdircosu_(int *index);
extern double getpdircosv_(int *index);
extern double getpdircosw_(int *index);
/*
 * Given the index of the emitted photon, this function returns
 * the direction cosines of its velocity.
 */

extern double getnage_(int *index);
/*
 * Given the index of the emitted neutron, this function returns
 * its age, -1 if index is out of range.
 * This age will be different from the time specified
 * in generateFissionEvent and generateSpontaneousFissionEvent
 * for non-prompt neutrons, i.e. delayed neutrons. 
 */

extern double getpage_(int *index);
/*
 * Given the index of the emitted photon, this function returns
 * its age, -1 of index is out of range.
 *  This age will be different from the time specified
 * in generateFissionEvent and generateSpontaneousFissionEvent
 * for photons that are emitted by beta-decay of the fission
 * fragments.
 */

extern void setdelay_(int *delay);
/*
 * This function is called to enable delayed neutrons and photons
 * Input
 *      delay:
 *              0 (default) for strictly prompt neutrons and 
 *                photons
 *              1 (n/a) for prompt neutrons, prompt and delayed 
 *                photons
 *              2 (n/a) for prompt and delayed neutrons, prompt 
 *                photons
 *              3 (n/a) for prompt and delayed neutrons, prompt 
 *                and delayed photons
 */

extern void setcorrel_(int *correlation);
/*
 * This function is called to set the type of neutron photon correlation
 * Input
 *      correlation:
 *              0 (default) for no correlation between neutrons and
 *                photons
 *              1 for neutron energy conservation and photon energy 
 *                conservation using data in paper "Implementation 
 *                of Energy-Dependent Q Values for Fission" by R. Vogt, 
 *                B. Beck, D.A. Brown, F. Daffin, and J. Hedstrom.
 *              2 for neutron energy conservation and photon energy 
 *                conservation using data in paper "Energy-Dependent 
 *                Fission Q Values Generalized for All Actinides" by 
 *                R. Vogt.
 *              3 for neutron and gamma-ray sampled from the FREYA code.
 *                If isotope included in FREYA, options set by 
 *                setnudist_() and setcf252_() are ignored. Otherwise, 
 *                reverts back to correlation option 0.
 */

extern void setnudist_(int *nudist);
/*
 * This function is called to set the data to be sampled for the neutron
 * number distributions in induced fissions
 * Input
 *      nudist:
 *               0 to use the fit to the Zucker and Holden tabulated 
 *                 P(nu) distributions as a function of energy for 
 *                 U235, U238 and Pu239. Terrell for other isotopes.
 *               1 to use fits to the Zucker and Holden tabulated 
 *                 P(nu) distribution as a function of energy for 
 *                 U238 and Pu239, and a fit to the Zucker and Holden 
 *                 data (for E greater than 1MeV) as well as the 
 *                 Gwin, Spencer and Ingle data (at thermal energies, 
 *                 i.e. 0 MeV) as a function of energy for U235. 
 *                 Terrell for other isotopes.
 *               2 to use the fit to the Zucker and Holden tabulated
 *                 P(nu) distributions as a function of nubar for
 *                 U238 and Pu239, and a fit to the Zucker and Holden
 *                 data (for E greater than 1 MeV) and the Gwin, 
 *                 Spencer and Ingle tabulated P(nu) distributions 
 *                 (at 0 MeV) as a function of nubar for U235. 
 *                 The U238 fit is used for the U232, U234, U236 and 
 *                 U238 isotopes, the U235 fit for U233 and U235, the 
 *                 Pu239 fit for Pu239 and Pu241. Terrell for other 
 *                 isotopes.
 *               3 (default) to use the Zucker and Holden as well as
 *                 the Gwin, Spencer and Ingle tabulated
 *                 P(nu) distributions as a function of nubar. The 
 *                 tables have P(nu) distributions for 11 
 *                 energies (0 MeV through 10 MeV), along with their
 *                 nubars. We select the P(nu) distribution that has
 *                 a nubar closest either from above, or from below, 
 *                 to the nubar entered for the induced fission, 
 *                 based on a random number and fractional distances 
 *                 to the end of the nubar interval thus formed.
 *                 The U238 fit is used for the U232, U234, U236 and 
 *                 U238 isotopes, the U235 fit for U233 and U235, the 
 *                 Pu239 fit for Pu239 and Pu241. Terrell for other 
 *                 isotopes.
 * For correlation option 3, FREYA computes the number of sampled 
 * fission neutrons independently, and ignores this setting.
 */


extern void setcf252_(int *ndist, int *neng);
/*
 * This function is called to set the data to be sampled for the 
 * (a) Cf252 spontaneous fission number distribution, and 
 * (b) Cf252 spontaneous fission neutron energy spectrum
 * Input
 *      ndist:
 *              0 (default) to sample the number of neutrons from the 
 *                tabulated data measured by Spencer
 *              1 to sample the number of neutrons from Boldeman's data
 *      neng:
 *              0 to sample the spontaneous fission neutron energy from 
 *                Mannhart corrected Maxwellian spectrum
 *              1 to sample the spontaneous fission neutron energy from 
 *                Madland-Nix theoretical spectrum
 *              2 to sample the spontaneous fission neutron energy from 
 *                the Froehner Watt spectrum
 * For correlation option 3, FREYA computes the number and energies of 
 * sampled fission neutrons independently, and ignores this setting.
 */

extern double getnsepeng_(int *nuclide);
/*
 * This function returns the neutron separation energy for the
 * nuclide ZA specified, that is the excitation energy that a nuclide 
 * ZA would be in if a neutron with no kinetic energy was added to a 
 * nuclide with with one less neutron than the input nuclide [Z(A-1)].
 * Returns -1 if the neutron separation energy is not available for 
 * that nuclide.
*/

extern void setrngf_(float (*funcptr) (void));
/*
 * This function sets the random number generator to the user-defined
 * one specified in the argument. If either setrngf_ or setrngd_ are
 * not specified, the default system call srand48 will be called.
 * Input
 *      funcptr:
 *               a random number generator function that returns a
 *               variable of type float
 */

extern void setrngd_(double (*funcptr) (void));
/*
 * This function sets the random number generator to the user-defined
 * one specified in the argument. If either setrngf_ or setrngd_ are
 * not specified, the default system call srand48 will be called.
 * Input
 *      funcptr:
 *               a random number generator function that returns a
 *               variable of type double
 */

#ifdef FREYA
extern void getfreya_errors_(int *length, char* errors);
/*
 * This function returns an array of characters with errors that 
 * were generated by FREYA. The length of the array errors passed as
 * an argument is given in the first argument length. The length of
 * the array returned is given in the first argument length.
 * Input
 *      errors:
 *              pointer to an allocated array of characters
 *      length:
 *              size of character array errors passed in
 * Output
 *      errors:
 *              error message (if any)
 *      length:
 *              size of error message (if any)
 *              =1 if no errors
 *              >1 if errors
 */
#endif
#ifdef __cplusplus
}
#endif
