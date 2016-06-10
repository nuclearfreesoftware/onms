! Copyright (c) 2013-2015, Lawrence Livermore National Security, LLC.
! Produced at the Lawrence Livermore National Laboratory
! Written by Ramona Voga <vogt2@llnl.gov>, 
!            Jørgen Randrup <jrandrup@lbl.gov>, 
!            Christian Hagmann <hagmann1@llnl.gov>, 
!            Jérôme Verbeke <verbeke2@llnl.gov>.
! LLNL-CODE-636753.
! OCEC-13-161
! All rights reserved.
! 
! This file is part of FREYA, Version:  1.0. For details, see <http://nuclear.llnl.gov/simulations>. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
! 
! o Redistributions of source code must retain the above copyright notice, this list of conditions and the disclaimer below.
! 
! o Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the disclaimer (as noted below) in the documentation and/or other materials provided with the distribution.
! 
! o Neither the name of the LLNS/LLNL nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
! 
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
! IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! 
! Additional BSD Notice
! 
! 1. This notice is required to be provided under our contract with the U.S. Department of Energy (DOE). This work was produced at Lawrence Livermore National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
! 
! 2. Neither the United States Government nor Lawrence Livermore National Security, LLC nor any of their employees, makes any warranty, express or implied, or assumes any liability or responsibility for the accuracy, completeness, or usefulness of any information, apparatus, product, or process disclosed, or represents that its use would not infringe privately-owned rights.
! 
! 3. Also, reference herein to any specific commercial products, process, or services by trade name, trademark, manufacturer or otherwise does not necessarily constitute or imply its endorsement, recommendation, or favoring by the United States Government or Lawrence Livermore National Security, LLC. The views and opinions of authors expressed herein do not necessarily state or reflect those of the United States Government or Lawrence Livermore National Security, LLC, 
!
       module freyaparameters
!      parameters
!      ----------
       implicit none
       save

       integer mxK            ! Number of cases
       integer Mxth
       parameter (Mxth=3)     ! Up to Mxth prefission neutrons
       integer MxEGn
       parameter (MxEGn=1001) ! Max entries in table of Gn/(Gn+Gf)
       integer Nk
       parameter (Nk=1)       ! Number of ejectile types (photons:0)
       integer mxE
       parameter (mxE=300)    ! Maximum number of energy bins (of width=dE)
       integer Mxg
       parameter (Mxg=2)      ! Maximum number of asymmetric fission modes
       integer mMax
       parameter (mMax=30)    ! Max n+g mult in each chain: pre & post1&2
       end module freyaparameters

       module freyaG
!      comG common block
!      -----------------
       use freyaparameters
       implicit none
       save

       real Gnf
       real dEnf
       real,dimension(:,:,:),allocatable :: Pnf   ! Pnf=Gn/(Gn+Gf) [not needed here]
       end module freyaG

       module freyaEv
!      comEv common block
!      ------------------
       implicit none
       save

       integer KK, kEv                         ! Event number (just for monitoring)
       end module freyaEv

       module freyaout
!      COMout common block
!      -------------------
       implicit none
       save

       integer L5,L6,L7,L8,L9                  ! Output channels
       data L5,L6,L7,L8,L9/5,6,7,8,9/
       end module freyaout

       module freyaB
!      comB common block
!      -----------------
!      Thresholds for evaporation and fission
       use freyaparameters
       implicit none
       save

       real,dimension(:,:),allocatable :: Sn   ! Neutron separation energy
       real,dimension(:,:),allocatable :: Bf   ! Fission barrier
       end module freyaB

       module freyaMth
!      cMth common block
!      -----------------
       implicit none
       save

       integer,dimension(:),allocatable :: Mthk! Max number of pre-fission neutrons
       end module freyaMth

       module freyaK
!      comK common block
!      -----------------
       implicit none
       save

       integer,dimension(:),allocatable :: iZk             ! Z of up to mK species
       integer,dimension(:),allocatable :: iAk             ! A of up to mK species
       integer,dimension(:),allocatable :: ntkek           ! number of TKE distributions
       real :: EnMax                                       ! maximum incident neutron energy
       real,dimension(:),allocatable :: Emax               ! Maximum excitation
       character*100,dimension(:),allocatable :: pAfk      ! name of the p(Af) distribution file
       character*100,dimension(:),allocatable ::  preprob  ! names of the files giving the probability 
                                                           ! for pre-equillibrium emission
       character*100,dimension(:),allocatable :: prespec   ! names of the files giving the pre-
                                                           ! equillibrium emission neutron spectra
       character*100,dimension(:,:),allocatable :: tkek    ! names of the total kinetic energy files
       character*30,dimension(:),allocatable :: label      ! Label for process
       character*10,dimension(:),allocatable :: reactlab   ! Label for reaction ( sf, (n,f), ...)
       integer,dimension(:),allocatable :: react           ! reaction index:
                                                           ! 0:sf
                                                           ! 1:(n,f)
       end module freyaK

       module freyaZ
!      comZ common block
!      -----------------
       implicit none
       save

       integer mdZ
       data mdZ/5/                     ! mdZ*Zdis is max Zf dev from <Zf>
       real,dimension(:),allocatable :: Zdis
       end module freyaZ

       module freyaC
!      comC common block
!      -----------------
       implicit none
       save

       real esq,c13,c23,c53
       end module freyaC

       module freyacmass
!      cmass common block
!      ------------------
       implicit none
       save

       real wA,Dn,DH
       data wA/931.49386/              ! [MeV/c2] Audi & Wapstra: NPA595 (1995) 409
       data Dn/8.0713228/              ! [MeV/c2] Audi & Wapstra, Table A
       data DH/7.2889694/              ! [MeV/c2] Audi & Wapstra, Table A
       end module freyacmass

       module freyaconsts
!      consts common block
!      -------------------
       implicit none
       save

       integer iseed
       data iseed/123457/              ! Seed for RAN0
       real pi,twopi
       data    pi/3.141592653/
       data twopi/6.283185308/
       end module freyaconsts

       module freyakmass
!      kmass common block
!      ------------------
       implicit none
       save

       real,dimension(:),allocatable :: wk
       real,dimension(:),allocatable :: Dk
       end module freyakmass

       module freyaLD
!      LD common block
!      ---------------
       implicit none
       save

       real r0
       end module freyaLD

       module freyaCD
!      comCD common block
!      ------------------
       implicit none
       save

       real,dimension(:,:),allocatable :: cNZ
       real,dimension(:,:),allocatable :: dNZ
       character,dimension(:,:),allocatable :: char
       end module freyaCD

       module freyaWg
!      comWg common block
!      ------------------
       implicit none
       save

       real dEf                            ! Width of energy bins in compound nucleus
       data dEf/1.0/
       real,dimension(:,:,:,:),allocatable :: CnNE
       real,dimension(:,:,:,:),allocatable :: FnNE
       real,dimension(:,:,:,:),allocatable :: AnNE
       real,dimension(:,:,:,:),allocatable :: WnNE
       end module freyaWg

       module freyacSission
!      cSission common block
!      ---------------------
       implicit none
       save

       integer maxN1, maxZ1
       parameter (maxN1=130,maxZ1=60)        ! Max N & Z of light fission fragment
       real,dimension(:,:,:,:),allocatable :: Qsf          ! Q-value for spontaneous fission (gs)
       real,dimension(:,:,:,:),allocatable :: Esciss       ! Interaction energy at scission
       real,dimension(:,:,:,:,:),allocatable :: Escissk    ! Interaction energy at scission
       real,dimension(:,:,:,:),allocatable :: Q1sciss      ! Def erg of light fragm at scission
       real,dimension(:,:,:,:),allocatable :: Q2sciss      ! Def erg of heavy fragm at scission
#if defined(WRITEL6) || defined(WRITEL8)
       real,dimension(:,:,:,:),allocatable :: Qsciss       ! Energy available at scission
       real,dimension(:,:,:,:,:),allocatable :: Qscissk    ! Energy available at scission
#endif
       ! real,dimension(:,:,:,:),allocatable :: s12sciss   ! Center sep at scission: s12=c1+c2+d
       ! real,dimension(:,:,:,:,:),allocatable :: s12scissk! Center sep at scission: s12=c1+c2+d
       end module freyacSission

       module freyaSAMPLE
!      SAMPLE common block
!      -------------------
       implicit none
       save

       LOGICAL,dimension(:),allocatable :: sampleAf    ! Get Af by direct sampling of P(Af)?
       end module freyaSAMPLE

       module freyacYAf
!      cYAf common block
!      -----------------
       implicit none
       save

       real,dimension(:,:),allocatable :: YAfk
       real,dimension(:,:),allocatable :: YYAfk
       integer,dimension(:),allocatable :: minAfk
       integer,dimension(:),allocatable :: maxAfk
       end module freyacYAf

       module freyacLOOK
!      cLOOK common block
!      ------------------
       implicit none
       save

       LOGICAL LOOK/.FALSE./               ! controls test printouts
       end module freyacLOOK

       module freyaQmin
!      comQmin common block
!      --------------------
       implicit none
       save

       real Qmin                           ! Min Q-value & ejectile kin erg (numerics)
       data Qmin/0.010/                    ! Minimum accepted neutron kinetic energy
       end module freyaQmin

       module freyaPath
!      ----------------
       implicit none
       save

       integer maxDATAPATH
       parameter (maxDATAPATH=2000)        ! max length of environment variable
                                           ! FREYADATAPATH
       character(len=maxDATAPATH) freyadir ! path to freya data directory

       end module freyaPath

       module freyaparams
!      params common block
!      -------------------
       implicit none
       save

       real alevel0,xeps,dTKE
       real,dimension(:),allocatable :: alevelk        ! alevel
       real,dimension(:),allocatable :: xepsk          ! xeps
       real,dimension(:),allocatable :: dTKE0k         ! value of dTKE (if no file)
       character*100,dimension(:),allocatable :: dTKEfk! names of the files
                                                       ! containing dTKE vs energy
       logical,dimension(:),allocatable :: dTKEe       ! existence of the dTKE files
       integer,dimension(:),allocatable :: iKtodTKEfk  ! mapping between iK index and 
                                                       ! index of (Ein,dTKE) array
       integer,dimension(:),allocatable :: nEin        ! number of pairs (Ein,dTKE)
                                                       ! in the dTKE file
       real,dimension(:,:),allocatable :: dTKEk        ! dTKE versus 
       real,dimension(:,:),allocatable :: dTKE_Ek      ! energies
       end module freyaparams

       module freyaerrors
!      common block for errors
!      -----------------------
       implicit none
       save

!     static variable storing whether initerrors has already been called
      logical :: freyaerrorsinitialized=.FALSE. 

!
!      Table of error codes
!      --------------------
!      code    meaning
!         1    data file react.dat could not be found
!         2    error code out of bounds
!         3    case-specific error code out of bounds
!         4    data file SEPn.dat could not be found
!         5    data file fisbar.dat could not be found
!         6    data file alevel.dat could not be found
!         7    data file MassMNMS.dat could not be found
!         8    data file MassAudi.dat could not be found
!         9    data file inputparameters.dat could not be found
!        10    data file acorrection.dat could not be found
!        11    data file Zdis.dat could not be found
!        12    data file gaussfit.dat could not be found
!        13    environment variable FREYADATAPATH too long
!
       integer maxErrorCodes
       parameter (maxErrorCodes=13)
       character*256,dimension(:),allocatable :: errors! list of errors
       integer,dimension(:),allocatable :: severities  ! severity of errors
       integer,dimension(:),allocatable :: errorctr    ! error counter
!
!      Table of case-specific error codes
!      ----------------------------------
!      code    meaning
!         1    incoming neutron energy Einc out of bounds
!         2    excitation energy eps0 out of bounds
!         3    data file *-Af*.dat could not be found
!         4    data file *-Ktot*.dat could not be found
!         5    data file *.xs could not be found
!         6    data file *.PreEq could not be found
!
       integer maxCaseErrorCodes
       parameter (maxCaseErrorCodes=6)
       character*256,dimension(:),allocatable :: caseerrors! list of errors
       integer,dimension(:),allocatable :: caseseverities  ! severity of errors
       integer,dimension(:),allocatable :: caseerrorctr    ! error counter
!
       integer nerrors
       integer ncaseerrors

       logical errorflag
       logical casespec
       integer lasterrorcode

       character*256 error                             ! error message
       integer ctr                                     ! error counter
       end module freyaerrors

       module freyarng
!      common block for random number generator
!      ----------------------------------------
       use iso_c_binding, only: C_BOOL

       implicit none
       save

       logical (kind=c_bool) :: user_rng=.FALSE.

       end module freyarng
