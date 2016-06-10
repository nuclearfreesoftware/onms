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
       SUBROUTINE msfreya_setup () &
       bind (C, name="msfreya_setup_")

! called from msFREYA_main to set up for event generation with msFREYA_event:

!**************************************************************************
       use freyaK
       use freyaG
       use freyaB
       use freyaout
       use freyaMth
       use freyaEv
       use freyaparams
       use freyaZ
       use freyaC
       use freyacmass
       use freyaconsts
       use freyakmass
       use freyaLD
       use freyaWg
       use freyacSission
       use freyaSAMPLE
       use freyacYAf
       use freyaPath
       use freyaerrors

       use freyainterfaces, only: msfreya_errorflagset

       implicit none

       logical exitonerror

       real,allocatable,dimension(:,:):: Ethr    ! E* threshold for m pre-fiss neutrons
       integer,allocatable,dimension(:,:):: mxEps! Max E* after emission of Nth neutrons

       character*3 element             ! Element (U, Pu, ...)

       character*70 headline           ! Headline of data file
       character*300 line              ! line in data file

       integer mxEin
       parameter (mxEin=200)           ! Max number of incident neutron 
                                       ! energies in files dTKE vs erg
       character*10 react0             ! name of reaction (sf,(n,f))
       character*100 dTKEf             ! name of temporary file where dTKE 
                                       ! vs energy is saved
       integer ndTKEfiles              ! number of files dTKE vs energy
       integer maxnEin                 ! maximum number of entries in
                                       ! the files dTKE vs energy
       character*10 iAstr
       real hbarc,D,gsM,wN,alevel,aA,U,epsmx &
       ,BARRIER,Qk,E,En,R,T,rho0,Gammanf,Gn,Gf &
       ,eps,Pn,alog10,alev,xeps0,dTKE0,msFREYA_SEPn
       integer iost,i,k,l,iK,iZ,iA,iA0,Nth,m
       character(len=maxDATAPATH+100) dataf 
                                       ! name of data file
       logical exists
       integer length

!************************************************************************
!      write (L6,"(80('='))")
#ifdef WRITEL6
       write (L6,*) 'Setting up ms FREYA ...'
#endif
       EnMax=20.0                      ! Max incoming neutron energy

       call initerrors()

!========================================================================
! find the path to the freya data directory
! use local 'data' directory by default
!       
       freyadir="data"
       dataf=trim(freyadir)//"/react.dat"
       inquire(FILE=dataf, EXIST=exists)
       if (exists.eqv..FALSE.) then
! check if FREYADATAPATH specified
         call get_environment_variable('FREYADATAPATH',freyadir,length)
         if (length.gt.maxDATAPATH) then
           write(error, 11) maxDATAPATH
11         format("FREYA: path in environment variable FREYADATAPATH ",&
                  "exceeds maximum allowable length (maxDATAPATH=", &
                  i4,").")
           call seterror(13,error)
           RETURN
         endif
         dataf=trim(freyadir)//"/react.dat"
         inquire(FILE=dataf, EXIST=exists)
         if ((exists.eqv..FALSE.).or.(trim(freyadir).eq.'')) then
           write(error, 12)
12         format("FREYA: data file react.dat could not be found. ", &
                  "Use either local data directory or specify Freya ", &
                  "data directory using variable FREYADATAPATH.")
           call seterror(1,error)
           RETURN
         endif
       endif

! include all processes available in file react.dat
       mxK=1                           ! # mxK identifies the case
       OPEN  (L8,file=trim(freyadir)//'/react.dat',status='old')
       read  (L8,*) headline
       do 
         read (L8,'(a300)', IOSTAT=iost) line
         if (iost.eq.0) then
           mxK=mxK+1
         else
           ! something went wrong reading file react.dat,
           ! or reached end of file
           exit
         end if
       end do  ! loop over file 'react.dat'
       mxK=mxK-1
       CLOSE (L8)

! allocate initial memory to arrays
       allocate(iZk(1:mxK))
       iZk=0
       allocate(iAk(1:mxK))
       iAk=0
       allocate(ntkek(1:mxK))
       ntkek=0
       allocate(Emax(1:mxK))
       Emax=0
       allocate(pAfk(1:mxK))
       allocate(preprob(1:mxK))
       allocate(prespec(1:mxK))
       allocate(label(1:mxK))
       allocate(reactlab(1:mxK))
       allocate(react(1:mxK))
       allocate(tkek(1:mxK,1:3))

       allocate(Pnf(0:MxEGn,0:Mxth,1:mxK))
       Pnf=0
       allocate(Sn(0:Mxth,1:mxK))
       Sn=0
       allocate(Bf(0:Mxth,1:mxK))
       Bf=0

       allocate(Mthk(1:mxK))
       Mthk=0

       allocate(Zdis(1:mxK))
       Zdis=0

       allocate(CnNE(0:Mxg,0:Mxth,mxE,mxK))
       CnNE=0
       allocate(FnNE(0:Mxg,0:Mxth,mxE,mxK))
       FnNE=0
       allocate(AnNE(0:Mxg,0:Mxth,mxE,mxK))
       AnNE=0
       allocate(WnNE(0:Mxg,0:Mxth,mxE,mxK))
       WnNE=0

       allocate(Qsf(maxN1,maxZ1,0:Mxth,mxK))
       Qsf=0
       allocate(Q1sciss(maxN1,maxZ1,0:Mxth,mxK))
       Q1sciss=0
       allocate(Q2sciss(maxN1,maxZ1,0:Mxth,mxK))
       Q2sciss=0
#if defined(WRITEL6) || defined(WRITEL8)
       allocate(Qsciss(maxN1,maxZ1,0:Mxth,mxK))
       Qsciss=0
#endif
       ! allocate(s12sciss(maxN1,maxZ1,0:Mxth,mxK))
       ! s12sciss=0
       allocate(Esciss(maxN1,maxZ1,0:Mxth,mxK))
       Esciss=0

       allocate(sampleAf(mxK))
       sampleAf=.FALSE.

       allocate(YAfk(300,mxK))
       YAfk=0
       allocate(YYAfk(300,mxK))
       YYAfk=0
       allocate(minAfk(mxK))
       minAfk=0
       allocate(maxAfk(mxK))
       maxAfk=0

       allocate(Ethr(0:Mxth,mxK))
       Ethr=0
       allocate(mxEps(0:Mxth,mxK))
       mxEps=0

       OPEN  (L8,file=trim(freyadir)//'/react.dat',status='old')
       read  (L8,*) headline
       do k=1,mxK
         read (L8,'(a300)', IOSTAT=iost) line
         if (iost.eq.0) then
           read (line,*) element,iZk(k),iAk(k),reactlab(k), &
                         Mthk(k),pAfk(k),preprob(k),prespec(k), &
                         ntkek(k),(tkek(k,i),i=1,ntkek(k))
           write(iAstr,'(i3)') iAk(k)
           label(k)=trim(iAstr)//trim(element)//trim(reactlab(k))
         else
           ! something went wrong reading file react.dat,
           ! or reached end of file
           exit
         end if
       end do  ! loop over file 'react.dat'
!
! Find Emax and react type for each case
!
       do k=1,mxK
         if (reactlab(k).eq."sf") then
           Emax(k)=0
           react(k)=0
         else
           call msFREYA_reseterrorflag()
           Emax(k)=msFREYA_SEPn(k,iZk(k),iAk(k))+EnMax
           if (msFREYA_errorflagset().and.exitonerror()) RETURN
           react(k)=1
         end if
       end do
       CLOSE (L8)

! Choose which processes to include in this run:
#ifdef WRITEL6
       write (L6,"(80('-'))")
       write (L6,*) 'These cases are included in this run:'
       do iK=1,mxK
         write (L6,"('Case',i3,':  #',i2,': ',a10,2i4,f8.3,i3)") &
           iK,iK,label(iK),iZk(iK),iAk(iK),Emax(iK),Mthk(iK)
       enddo
       write (L6,"(80('-'))")
#endif

!************************************************************************
! Constants:
       c13      =1.0/3.0
       c23      =2.0/3.0
       c53      =5.0/3.0
       esq      =1.44                  ! alpha*hbar*c = e**2/4*pi*eps0
       hbarc    =197.33
       r0       =1.2                   ! Nuclear radius constant (Lysekil)

       call SetupRandom (0,0)          ! set up & check random generators
#ifdef WRITEL6
       write (L6,*) 'START value of iseed:',iseed!!
#endif

       call msFREYA_reseterrorflag()
       D = gsM (0,0)                   ! Initialize gsM (iZ,iA)
       if (msFREYA_errorflagset().and.exitonerror()) RETURN

! allocate initial memory to arrays
       allocate(wk(0:Nk))
       wk=0
       allocate(Dk(0:Nk))
       Dk=0

! Ejectile masses:
       wN=wA+Dn                        ! Neutron mass (MeV/c2): 939.565128
       Dk(0)=0.0; wk(0)=0.0            ! gamma
       Dk(1)=gsM(0,1); wk(1)=wN        ! neutron

!========================================================================
! These default parameters could be species dependent:
       alevel0  = 10.0724      ! 7.25  !7.010512
       xeps     = 1.23389      ! 1.0   !0.241473
       dTKE     = 1.53729      !-2.0   !1.45725	Should depend on (E*)!
!
! allocate memory
!
       allocate(alevelk(mxK))
       allocate(xepsk(mxK))
       allocate(dTKE0k(mxK))
       allocate(dTKEfk(mxK))
       allocate(iKtodTKEfk(mxK))
       allocate(dTKEe(mxK))
!
! supply default values for alevel0, xeps and dTKE, in case no entry 
! for Z,A,reaction in file 'inputparameters.dat'
!
       do iK=1,mxK
         dTKEe(iK)=.FALSE.
         alevelk(iK)=alevel0
         xepsk(iK)=xeps
         dTKE0k(iK)=dTKE
       end do
!
! read all alevel0, xeps, dTKE values, and
! determine how many dTKE vs erg files we have
!
       dataf=trim(freyadir)//'/inputparameters.dat'
       inquire(FILE=dataf, EXIST=exists)
       if (exists.eqv..FALSE.) then
#ifdef WRITEL6
         write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
         write(error, 13) trim(dataf)
13       format("FREYA: data file ", a, " could not be found.")
         call seterror(9,error)
         RETURN
       endif
       ndTKEfiles=0
       OPEN  (L8,file=dataf,status='old')
       read  (L8,*) headline
       do
         read (L8,'(a300)', IOSTAT=iost) line
         if (iost.ne.0) exit
         read (line,*) iZ,iA,react0,alev,xeps0,dTKEf,dTKE0
         do iK=1,mxK
           if (iZ.eq.iZk(iK).and.iA.eq.iAk(iK).and. &
               react0.eq.reactlab(iK)) then
             alevelk(iK)=alev
             xepsk(iK)=xeps0
             dTKE0k(iK)=dTKE0
             dTKEfk(iK)=dTKEf
             dataf=trim(freyadir)//"/"//trim(dTKEfk(iK))
             inquire(FILE=dataf, EXIST=dTKEe(iK))
             if (dTKEe(iK).eqv..TRUE.) then
               ndTKEfiles=ndTKEfiles+1
             endif
           endif
         end do
       end do
       close(L8)
!
! allocate memory for the arrays that will contain the dTKE's
!
       allocate(nEin(ndTKEfiles))
       nEin=0
       allocate(dTKEk(ndTKEfiles,mxEin))
       allocate(dTKE_Ek(ndTKEfiles,mxEin))
!
! read all 'ndTKEfiles' files
!
       i=0
       do iK=1,mxK
         if (dTKEe(iK).eqv..FALSE.) cycle
         i=i+1
         iKtodTKEfk(iK)=i
         dataf=trim(freyadir)//"/"//trim(dTKEfk(iK))
         OPEN (L8,file=dataf,status='old')
         read (L8,901) headline
  901    format (a70)
! save energies and dTKE in arrays
         read (L8,*) nEin(i)                  ! Number of values Ein>0
         do l=1,nEin(i)
           read (L8,*,IOSTAT=iost) dTKE_Ek(i,l), dTKEk(i,l)
           if (iost.ne.0) then
             exit
           endif 
         enddo
         CLOSE (L8)
       end do
!
! allocate memory to remaining arrays
!
       maxnEin=maxval(nEin)
#if defined(WRITEL6) || defined(WRITEL8)
       allocate(Qscissk(maxN1,maxZ1,0:Mxth,ndTKEfiles,maxnEin))
       Qscissk=0
#endif
       ! allocate(s12scissk(maxN1,maxZ1,0:Mxth,ndTKEfiles,maxnEin))
       ! s12scissk=0
       allocate(Escissk(maxN1,maxZ1,0:Mxth,ndTKEfiles,maxnEin))
       Escissk=0
!
#ifdef WRITEL6
       write (L6,*) 'Preparing fission-fragment level densities:'
#endif
       call msFREYA_reseterrorflag()
       aA = alevel (0,0,0,0.0,U)                ! Initialize alevel (iZ,iA,eps,e,U)
       if (msFREYA_errorflagset().and.exitonerror()) RETURN

       DO iK=1,mxK               ! LOOP OVER mxK CASES:	=================
         iZ=iZk(iK); iA0=iAk(iK); epsmx=Emax(iK)
         do Nth=0,Mthk(iK)
           call msFREYA_reseterrorflag()
           Sn(Nth,iK) = msFREYA_SEPn (iK,iZ,iA0-Nth)
                                                ! Separation energy
           if (msFREYA_errorflagset().and.exitonerror()) RETURN
           call msFREYA_reseterrorflag()
           Bf(Nth,iK) = BARRIER (iK,iZ,iA0-Nth) ! Fission barrier
           if (msFREYA_errorflagset().and.exitonerror()) RETURN
           Qk=gsM(iZ,iA0-Nth)-gsM(iZ,iA0-Nth-1)-Dn
           E=Bf(Nth,iK)
           do m=Nth-1,0,-1
             E=E+Sn(m,iK)
           enddo
           Ethr(Nth,iK)=E; En=E-Sn(0,iK)
#ifdef WRITEL6
           write (L6,*) iZ,iA0-Nth,Bf(Nth,iK),Sn(Nth,iK),Sn(Nth,iK)+Qk
           if (Bf(Nth,iK).lt.Sn(Nth,iK)) then
             write (L6,811) Nth,iZ,iA0,Bf(Nth,iK),'<', &
                           Sn(Nth,iK),Sn(Nth,iK)+Qk,E,En
           else if (Bf(Nth,iK).gt.Sn(Nth,iK)) then
             write (L6,811) Nth,iZ,iA0,Bf(Nth,iK),'>', &
                           Sn(Nth,iK),Sn(Nth,iK)+Qk,E,En
           else
             write (L6,811) Nth,iZ,iA0,Bf(Nth,iK),'=', &
                           Sn(Nth,iK),Sn(Nth,iK)+Qk,E,En
           endif
#endif
           iA=iA-1                              ! Emit one neutron
         enddo
  811    format (3i5,f10.5,2x,a1,5f10.5)
#ifdef WRITEL6
         m=0; write (L6,"('   Nth    E*thr   <   E*max:      above')")
         do while (Ethr(m,iK).lt.epsmx.AND.m.le.Mthk(iK))
           write (L6,"(i5,3f12.5)") m,Ethr(m,iK),epsmx,epsmx-Ethr(m,iK)
           m=m+1
         enddo
         write (L6, 923) m
  923    format('Max number of pre-fiss neutrons possible', &
                ' (if Pnf(E*) had a SHARP cut-off):',i3)
#endif

! Calculate chances for pre-fission neutron emission ("N'th chance fission"):
!         Pnf is tabulated versus the EXCITATION in the emitting nucleus
         dEnf=0.2
         do Nth=0,Mthk(iK)
           iA=iA0-Nth                   ! A of current compound nucleus
#ifdef WRITEL6
           write (L6,*) 'Nth =',Nth,' =>  E* =',epsmx,'in A=',iA,':'
#endif
           R=1.15*float(iA)**c13        !.NE. Lysekil value 1.2 sigma=pi*R**2
!          Gnf:=2*gs*mN*sigma/pi*hbar*2, 
           Gnf=4*wN*(R/hbarc)**2
#ifdef WRITEL6
           write (L6,"(' Gnf:',f8.4,'/MeV = 1/(',f7.4,' MeV)')") Gnf, &
               1/Gnf
#endif
!          write (L8,*) '  Nth   E*     T       Gn        GnX  ' &
!           ,'      Gf        GfX      Pn',iA
           m=0; Pnf(0,Nth,iK)=0.0; eps=0.0
           do while (eps.le.epsmx)
             m=m+1; eps=m*dEnf
             aA=alevel(iK,iZ,iA,eps,U)
             T=sqrt(U/aA)               ! Temperature
             call msFREYA_reseterrorflag()
             rho0 = Gammanf (iK,iZ,iA,eps,Sn(Nth,iK),Bf(Nth,iK),Gn,Gf)
             if (msFREYA_errorflagset().and.exitonerror()) RETURN
             Pn=Gn+Gf; if (Pn.gt.0.0) Pn=Gn/Pn
             if (Gn.gt.0.0) Gn=alog10(Gn)
             if (Gf.gt.0.0) Gf=alog10(Gf)
             Pnf(m,Nth,iK)=Pn
           enddo
           Pnf(m+1,Nth,iK)=Pnf(m,Nth,iK)! for interpolation
           mxEps(Nth,iK)=m              ! Number of energy bins for this Nth
           epsmx=epsmx-Sn(Nth,iK)       ! Reduce max excitation
         enddo

         call msFREYA_reseterrorflag()
         call SCISSION (iK)             ! Get scission info: TKEL and Q1 & Q2.
         if (msFREYA_errorflagset().and.exitonerror()) RETURN
         call msFREYA_reseterrorflag()
         call MASSES   (iK,Emax(iK))    ! Prepare fragment mass distribution
         if (msFREYA_errorflagset().and.exitonerror()) RETURN
       ENDDO                    ! LOOP OVER mxK CASES.	=================

       aA = alevel (0,0,0,0.0,U)        ! Initialize alevel (iZ,iA,eps,e,U)

#ifdef WRITEL6
       write (L6,"(80('-'))")
       write (L6,*)'SETUP SUMMARY: #   Z   A    E*max  Mth   label'
       do iK=1,mxK
         write (L6,"('Case',i4,':',i8,i4,i5,f8.3,1x,i4,3x,a13)") &
               iK,iK,iZk(iK),iAk(iK),Emax(iK),Mthk(iK),label(iK)
       enddo
       write (L6,"(i2,' cases are considered in this msFREYA run.')") &
         mxK
       write (L6,"(80('='))")
#endif
       kEv=0                            ! Event number (just for monitoring)
       RETURN
       END

!************************************************************************
       SUBROUTINE msfreya_getniso(nisosf,nisoif) &
       bind (C, name="msfreya_getniso_")
! must be called after msfreya_setup, returns numbers of
! spontaneous and induced fission isotopes
!
! Output:  nisosf (number of spontaneous fission isotopes)
!          nisoif (number of induced fission isotopes)
!
       use freyaparameters
       use freyaK
       use iso_c_binding, only: C_INT

       implicit none

       integer (kind=c_int) :: nisosf,nisoif

       integer k

       nisosf=0
       nisoif=0
       do k=1,mxK
         if (react(k).eq.0) then
           nisosf=nisosf+1
         end if
         if (react(k).eq.1) then
           nisoif=nisoif+1
         end if
       end do

       RETURN 
       end

!************************************************************************
       SUBROUTINE msfreya_getzas(zas,fistypes) &
       bind (C, name="msfreya_getzas_")
! must be called after msfreya_setup, returns 
! - one array with the ZAs of the fission isotopes
! - one array with the fission type 
!   o spontaneous (0)
!   o induced (1)
!
! Output:  zas      - array of ZAs of fission isotopes
!          fistypes - array of fission type
!
       use freyaparameters
       use freyaK
       use iso_c_binding, only: C_INT, C_FLOAT

       implicit none

       integer (kind=c_int), dimension(*) :: zas
       integer (kind=c_int), dimension(*) :: fistypes

       integer k

       do k=1,mxK
         if (react(k).eq.0) then
           zas(k)=1000*iZk(k)+iAk(k)
         elseif (react(k).eq.1) then
           zas(k)=1000*iZk(k)+iAk(k)-1
         endif
         fistypes(k)=react(k)
       end do

       RETURN 
       end

!************************************************************************
       function msFREYA_gsMassN (iZ,iA) &  ! Ground-state mass of nucleus
       bind (C, name="msfreya_gsmassn_")
!
!      returns the ground state mass of nucleus iZ,iA
!

       use freyacmass
       use iso_c_binding, only: C_INT, C_FLOAT

       implicit none

       real (kind=c_float) :: msFREYA_gsMassN
       integer (kind=c_int), value :: iZ, iA

       real gsM

       msFREYA_gsMassN = iA*wA+gsM(iZ,iA)

       RETURN
       end

!************************************************************************
       FUNCTION DROP (iZ,iA)
! CALCULATES THE LYSEKIL NUCLEAR MASS DEFECT + WIGNER + PAIRING + SHIFT
! Change of notation: WN -> Dn, WP -> DH, W -> D for clarity [JR:30jan08].
! The ground-state mass of M(Z,A) = A*wA + DROP(Z,A) with wA=931.504,
! so the total binding energy is B(Z,A) = Z*MH + (A-Z)*Mn - M(Z,A) =
! Z*(MH-wA) + (A-Z)*(Mn-wA) - DROP(Z,A) = Z*WH + (A-Z)*Wn + DROP(Z,A).
! The neutron separation energy is
! Note: This mass 'defect' is often called mass 'excess'.

       use freyacmass

       implicit none

       real drop
       integer iZ,iA

       real A1,A2,C3,CAPPA,C4,WIG,D1,D2,SHIFT
       DATA A1,A2,C3,CAPPA/15.4941,17.9439,0.7053,1.7826/
       DATA C4/1.1533/
       DATA WIG,D1,D2,SHIFT/30.,12.,10.,50./

        integer IA2,IZ2
        real A,Z,X,T,A3,EV,ES,EC,D,E,EW,EP

       A=float(iA)
       Z=float(iZ)
       X=0.0
       DROP=0.0
       IF (A.LE.0.) RETURN

       T=ABS(1.-2.*Z/A)
       A3=A**0.333333
       EV=-A1*A*(1.-CAPPA*T**2)         ! Volume energy
       ES=A2*A3**2*(1.-CAPPA*T**2)      ! Surface energy
       EC=(C3/A3-C4/A)*Z**2             ! Coulomb energy
       D=(A-Z)*Dn+Z*DH
       DROP=D
       X=0.5*EC/ES                      ! Fissility (WJS)
       IF (A.LE.4.) RETURN

       E=EV+ES+EC+D
       EW=WIG*T
       IA2=int(A/2.+0.1)
       IF (A-2.*IA2-0.1) 10,10,20
   10  IZ2=int(Z/2.+0.1)
       EP=D1/SQRT(A)-D2/A
       IF (Z-2.*IZ2-0.1) 11,11,15
   11  EW=EW-EP
       GO TO 30
   15  EW=EW+EP
       IF (IA2.EQ.2*IZ2) EW=EW+WIG/A
       GO TO 30
   20  EW=EW+D2/A
   30  E=E+EW+SHIFT/A
       DROP=E

       RETURN
       END
     
!************************************************************************
       function gsM (iZ,iA)    ! Ground-state mass defect D(Z,A):
! Ground-state mass DEFECT yields M(Z,A) [MeV] = A + D(Z,A)/931.49386
! whereas we want M = A*wA + gsM with wA=931.504 MeV.
!NOTE: The function must be initialized before usage by c =gsM (0,0)!
! iA<0:	Liquid-Drop mass is used.
! iA=0:	Setup: reading values from MassMNMS.dat & MassAudi.dat
! iA>0:	Looks up the value in the array dNZ: gsM = dNZ(iN=iA-iZ,iZ)
!------------------------------------------------------------------------

       use freyaout
       use freyaC
       use freyaCD
       use freyaPath
       use freyaPath
       use freyaerrors

       implicit none

       integer iZ,iA
       real gsM
      
       integer mxN,mxZ
       parameter (mxN=200,mxZ=100)

       logical MAP

       real flag
       data flag/1000./

       integer nN,nZ,nA
       real drop,wA,W,D,B
       character*70 headline   ! Headline of data file

       character*100 dataf             ! names of data file containing
                                       ! the mass tables
       logical exists

! allocate initial memory to arrays
       if (.not.allocated(cNZ)) then
         allocate(cNZ(0:mxN,0:mxZ))
         cNZ=0
       endif
       if (.not.allocated(dNZ)) then
         allocate(dNZ(0:mxN,0:mxZ))
         dNZ=0
       endif
       if (.not.allocated(char)) then
         allocate(char(0:mxN,0:mxZ))
       endif

!------------------------------------------------------------------------
       gsM=0.            ! default initialization
       MAP=.FALSE.       ! default initialization

       IF (iA.lt.0) THEN ! use Liquid-Drop value:
         gsM = DROP (iZ,iabs(iA))
         RETURN
       ENDIF
!------------------------------------------------------------------------
       IF (iA.eq.0) THEN ! initialize (read in the mass table):
! First we fill up the table with flags so usage can be checked:
         do nN=0,mxN
           do nZ=0,mxZ
             dNZ (nN,nZ) = flag
             char(nN,nZ)=' ' !1h 	! blank character
             if (mod(nZ,10).eq.0) char(nN,nZ)=' ' !1h|
           enddo
         enddo

#ifdef WRITEL6
         write (L6,*) 'Reading the mass tables ...'
#endif
! Then we read in the theoretical mass table:
         dataf=trim(freyadir)//'/MassMNMS.dat'
         inquire(FILE=dataf, EXIST=exists)
         if (exists.eqv..FALSE.) then
#ifdef WRITEL6
           write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
           write(error, 13) trim(dataf)
13         format("FREYA: data file ", a, " could not be found.")
           call seterror(7,error)
           RETURN
         endif
         OPEN  (L7,file=dataf,status='old')
         read  (L7,901) headline
#ifdef WRITEL6
         write (L6,902) headline
#endif
         read  (L7,901) headline
#ifdef WRITEL6
         write (L6,902) headline
#endif
         read  (L7,*) wA
#ifdef WRITEL6
         write (L6,*) 'gsM:',wA
#endif
         read  (L7,901) headline
#ifdef WRITEL6
         write (L6,902) headline
#endif
         read  (L7,901) headline
#ifdef WRITEL6
         write (L6,902) headline
#endif
         nA=1
         do while (nA.gt.0)
           read  (L7,*) nZ,nA,W,D,B
           if (nA.gt.0) then
             nN=nA-nZ
             if (nZ.le.mxZ.AND.nN.le.mxN) then
!              write (L6,*) nZ,nA,W,D,B
               dNZ (nN,nZ) = D
               char(nN,nZ)='0' !1h0
             endif
           endif
         enddo
         CLOSE (L7)
! and finally overwrite it with the experimental table
! in order to use experimental masses where available:
         dataf=trim(freyadir)//'/MassAudi.dat'
         inquire(FILE=dataf, EXIST=exists)
         if (exists.eqv..FALSE.) then
#ifdef WRITEL6
           write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
           write(error, 13) trim(dataf)
           call seterror(8,error)
           RETURN
         endif
         OPEN  (L7,file=dataf,status='old')
         read  (L7,901) headline
#ifdef WRITEL6
         write (L6,902) headline
#endif
         read  (L7,901) headline
#ifdef WRITEL6
         write (L6,902) headline
#endif
         read  (L7,*) wA
#ifdef WRITEL6
         write (L6,*) 'gsM: wA =',wA
#endif
         read  (L7,901) headline
#ifdef WRITEL6
         write (L6,902) headline
#endif
         read  (L7,901) headline
#ifdef WRITEL6
         write (L6,902) headline
#endif
         nA=1
         do while (nA.gt.0)
           read  (L7,*) nZ,nA,W,D,B
           if (nA.gt.0) then
             nN=nA-nZ
             if (nZ.le.mxZ.AND.nN.le.mxN) then
               dNZ (nN,nZ) = D
               char(nN,nZ) ='#' !1h#
             endif
           endif
         enddo
         CLOSE (L7)
#ifdef WRITEL6
         write (L6,*) 'Done reading the mass tables.'
#endif
         gsM = flag
! Chart of dNZ(N,Z):
         if (MAP) then
         else
#ifdef WRITEL9
           write (L9,*) 'NZ chart of nuclear masses:'
           write (L9,*) '  0: Theoretical masses from MNMS'
           write (L9,*) '  #: Experimental masses from Audi'
           write (L9,*) 'N            P R O T O N  N U M B E R  Z'
           do nN=0,mxN
             write (L9,101) nN,(char(nN,nZ),nZ=0,mxZ)
           enddo
#endif
         endif
  101    format (i3,1x,101a1)

! We thus use exp masses where available and theo values elsewhere.
!------------------------------------------------------------------------
       ELSE
         D = dNZ (iA-iZ,iZ)
         if (D.eq.flag) then ! this bin is unfilled:
#ifdef WRITEL6
           write (L6,*) 'gsM:',iZ,iA-iZ,iA,' has NOT been set!'
! Check:			OK up to KK=1,000,000
           write (L6,*) 'gsM:',iZ,iA-iZ,iA,': Continue with DROP?'
!			PAUSE!!
#endif
           D = DROP (iZ,iabs(iA))
         endif
         gsM = D
       ENDIF
  901  format (a70)
  902  format ('gsM: ',a70)

       RETURN
       END

!************************************************************************
! looks up the neutron separation energy from the table below
! which is taken from www.nndc.bnl.gov/masses/rct2.mas03:
! iA<0:	Separation energy is taken from MNMS mass table;
! iA>0:	Separation energy is taken from the table below:
!	The integer values are "estimated values" given in the data base.
!------------------------------------------------------------------------

       function msFREYA_SEPn (iK,iZ,iA)

       use freyaparameters
       use freyaK
       use freyaout
       use freyacmass
       use freyaPath
       use freyaerrors
 
       implicit none
       real :: msFREYA_SEPn
       integer :: iK,iZ,iA

       character*100 dataf               ! name of data file 'SEPn.dat'
       logical exists
       logical, save :: dataread=.FALSE. ! static variable storing whether file
                                         ! SEPn.dat was read
       character*70 headline             ! Headline of data file
       integer,save,dimension(:),allocatable:: iKtoSEPnk
                                         ! mapping between iK and index 
                                         ! of element in SEPnk array
       integer,save,dimension(:),allocatable::Zk 
                                         ! array of the isotope Z's
       integer,save,dimension(:),allocatable::Ak 
                                         ! array of the isotope A's
       real,save,dimension(:),allocatable::SEPnk 
                                         ! array of the separation energies
                                         ! corresponding to the ZA above
       integer, save :: nisotopes        ! number of separation energies/
                                         ! isotopes in file 'SEPn.dat'
       logical foundZ
       integer nZ,nA,iAf,isoindex,iost,iZA
       real w0,W,D,B,Wf,Wi,drop,sep

!*******
       msFREYA_SEPn = 0

       if (iA.lt.0) then                 ! Separation energy from MNMS mass table:
         dataf=trim(freyadir)//'/MassMNMS.dat'
         inquire(FILE=dataf, EXIST=exists)
         if (exists.eqv..FALSE.) then
#ifdef WRITEL6
           write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
           write(error, 13) trim(dataf)
13         format("FREYA: data file ", a, " could not be found.")
           call seterror(8,error)
           RETURN
         endif
         OPEN  (L7,file=dataf,status='old')
         read  (L7,901) headline
         read  (L7,901) headline
         read  (L7,*) w0                 ! wA
         read  (L7,901) headline
         read  (L7,901) headline
         nZ=-1
         do while (nZ.ne.iZ)
           read  (L7,*) nZ,nA,W,D,B
         enddo
         iAf=iabs(iA)-1
         do while (nA.ne.iAf)
           read  (L7,*) nZ,nA,Wf,D,B
         enddo
         read  (L7,*) nZ,nA,Wi,D,B
!          write (L6,*) 'Sn(Z,A) =',Wf+Wn-Wi,':',Wf,Wn,Wi
         msFREYA_SEPn = Wf + (wA+Dn) - Wi
         CLOSE (L7)
         RETURN
  901    format (a70)
       endif                    ! Separation energy from MNMS mass table.
!*******
       if (iA.ge.0) then        ! Separation energy from file SEPn.dat:
         if (dataread.eqv..FALSE.) then
           dataf=trim(freyadir)//"/SEPn.dat"
           inquire(FILE=dataf, EXIST=exists)
           if (exists.eqv..FALSE.) then
#ifdef WRITEL6
             write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
             write(error, 13) trim(dataf)
             call seterror(4,error)
             RETURN
           endif
           OPEN  (L7,file=dataf,status='old')
           read  (L7,901) headline
           isoindex=1
           do 
             read  (L7,*,IOSTAT=iost) nZ,nA,sep
             if (iost.ne.0) then
               ! something went wrong reading file
               ! SEPn.dat, or reached end of file
               exit
             endif 
             isoindex=isoindex+1
           enddo
           CLOSE (L7)
           nisotopes = isoindex-1
! allocate memory to iKtoSEPnk, Zk, Ak and SEPnk
           if (.not.allocated(iKtoSEPnk)) then
             allocate(iKtoSEPnk(mxK))
             iKtoSEPnk=0
           endif
           if (.not.allocated(Zk)) then
             allocate(Zk(nisotopes))
             Zk=0
           endif
           if (.not.allocated(Ak)) then
             allocate(Ak(nisotopes))
             Ak=0
           endif
           if (.not.allocated(SEPnk)) then
             allocate(SEPnk(nisotopes))
             SEPnk=0
           endif
!
           OPEN  (L7,file=dataf,status='old')
           read  (L7,901) headline
           do isoindex=1,nisotopes
             read  (L7,*,IOSTAT=iost) nZ,nA,sep
             if (iost.ne.0) then
               ! something went wrong reading file
               ! SEPn.dat, or reached end of file
               exit
             endif 
             ! The line was read properly
             Zk(isoindex) = nZ
             Ak(isoindex) = nA
             SEPnk(isoindex) = sep
             ! find mapping from iK to index of SEPn
             do iZA=1,mxK
               if (nZ.eq.iZk(iZA).and.nA.eq.iAk(iZA)) then
                 iKtoSEPnk(iZA)=isoindex
               end if
             end do
           enddo
           CLOSE (L7)
           dataread = .TRUE.
         endif ! dataread == false
!
! Look for the separation energy in array SEPnk
!
! See if we have this iK mapped onto one of the elements in
! the array SEPn
!
         isoindex=iKtoSEPnk(iK)
         if (isoindex.ne.0) then
           if (iZk(iK).eq.iZ.and.iAk(iK).eq.iA) then
             msFREYA_SEPn = SEPnk(isoindex)/1000
             RETURN
           endif
         endif
!
         foundZ=.FALSE.
         do isoindex=1,nisotopes
           if (iZ.eq.Zk(isoindex)) then
             foundZ=.TRUE.
             if (iA.eq.Ak(isoindex)) then
               msFREYA_SEPn = SEPnk(isoindex)/1000
               RETURN
             endif
           endif
         enddo ! loop over all ZA's
         if (foundZ.eqv..TRUE.) then
           ! Could find Z but not ZA in file SEPn.dat. In this 
           ! case use the liquid drop value for Sn. This part 
           ! of the code should no longer be invoked.
           msFREYA_SEPn = Dn + DROP (iZ,iA-1) - DROP (iZ,iA)
#ifdef WRITEL6
           write (L6,*) '** SEPn:',iZ,iA
           write (L6,*) '** This isotope has no tabulated value!'
           write (L6,*) '** Using the Liquid Drop value for Sn:', &
                         msFREYA_SEPn
#endif
           RETURN
         end if
       endif   ! Use separation energy from file SEPn.dat.
!*******
       if (iZ.ne.92.and.iZ.ne.94.and.iZ.ne.96.and.iZ.ne.98) then
         ! Use the liquid drop value for the Z's above, 0 otherwise.
         ! This part of the code should no longer be invoked.
#ifdef WRITEL6
         write (L6,*) '**',iZ,': No separation energy for this element!'
#endif
         msFREYA_SEPn = 0.0
       endif
       RETURN
       END


!************************************************************************     
!      This function serves as an "overloaded" function for msFREYA_SEPn      
!      where this function's parameters are type integer (kind=c_int) 
!      instead of the type integer as function msFREYA_SEPn requires.
!
       function msFREYA_SEPn_CVT (p1, p2, p3) &
         bind (C, name="msfreya_sepn_cvt_")
       use iso_c_binding, only: C_INT, C_FLOAT

       implicit none
       real (kind=c_float) :: msFREYA_SEPn_CVT
       real :: msFREYA_SEPn
       integer (kind=c_int), value :: p1, p2, p3
       integer :: a1, a2, a3
  
!      Convert the parameter data types from integer(kind=c_int) to           
!      integer                                                                
       a1 = int(p1)
       a2 = int(p2)
       a3 = int(p3)

!      Call msFREYA_SEPn with arguments of type integer                       
       msFREYA_SEPn_CVT = real(msFREYA_SEPn(a1, a2, a3),kind=c_float)
       RETURN
       END


!************************************************************************
       function BARRIER (iK,iZ,iA)
! looks up the fission barrier from the table below which is based on
!      http://iaeand.iaea.or.at/ripl/
       use freyaparameters
       use freyaK
       use freyaout
       use freyaPath
       use freyaerrors

       implicit none

       real BARRIER
       integer iK,iZ,iA

       logical, save :: dataread=.FALSE.
       character*100 dataf
       logical exists
       integer,save,dimension(:),allocatable::iKtofissbark
                                         ! mapping between iK and index 
                                         ! of element in fissbark array
       integer,save,dimension(:),allocatable::Zk 
                                         ! array of the isotope Z's
       integer,save,dimension(:),allocatable::Ak 
                                         ! array of the isotope A's
       real,save,dimension(:),allocatable::fissbark
                                         ! array of fission barriers
                                         ! corresponding to the ZA above
       integer, save :: nisotopes        ! number of fission barriers/
                                         ! isotopes in file 'fisbar.dat'
       logical foundZ
       integer isoindex,iost,nZ,nA,iZA
       real fissbar
       character*70 headline             ! Headline of data file
!      -----------------------------------------------------------------
       BARRIER=0.                        ! default initialization

       if (dataread.eqv..FALSE.) then
! read in file fisbar.dat, if it exists
         dataf=trim(freyadir)//"/fisbar.dat"
         inquire(FILE=dataf, EXIST=exists)
         if (exists.eqv..FALSE.) then
#ifdef WRITEL6
           write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
           write(error, 14) trim(dataf)
14         format("FREYA: data file ", a, " could not be found.")
           call seterror(5,error)
           RETURN
         end if
         OPEN  (L7,file=dataf,status='old')
         read  (L7,901) headline
         isoindex=1
         do 
           read  (L7,*,IOSTAT=iost) nZ,nA,fissbar
           if (iost.ne.0) then
             ! something went wrong reading file
             ! fisbar.dat, or reached end of file
             exit
           endif 
           isoindex=isoindex+1
         enddo
         CLOSE (L7)
         nisotopes = isoindex-1
! allocate memory to iKtofissbark, Zk, Ak and fissbark
         if (.not.allocated(iKtofissbark)) then
           allocate(iKtofissbark(mxK))
           iKtofissbark=0
         endif
         if (.not.allocated(Zk)) then
           allocate(Zk(nisotopes))
           Zk=0
         endif
         if (.not.allocated(Ak)) then
           allocate(Ak(nisotopes))
           Ak=0
         endif
         if (.not.allocated(fissbark)) then
           allocate(fissbark(nisotopes))
           fissbark=0
         endif
!
         OPEN  (L7,file=dataf,status='old')
         read  (L7,901) headline
         do isoindex=1,nisotopes
           read  (L7,*,IOSTAT=iost) nZ,nA,fissbar
           if (iost.ne.0) then
             ! something went wrong reading file
             ! fisbar.dat, or reached end of file
             exit
           endif 
           ! The line was read properly
           Zk(isoindex) = nZ
           Ak(isoindex) = nA
           fissbark(isoindex) = fissbar
           ! find mapping from iK to index of fissbark
           do iZA=1,mxK
             if (nZ.eq.iZk(iZA).and.nA.eq.iAk(iZA)) then
               iKtofissbark(iZA)=isoindex
             end if
           end do
         end do ! loop over entire file 'fisbar.dat'
         CLOSE  (L7)
         dataread=.TRUE.
       endif ! dataread == false
!
! Look for the fission barrier in array fissbark
!
! See if we have this iK mapped onto one of the elements in
! the array fissbark
!          
       isoindex=iKtofissbark(iK)
       if (isoindex.ne.0) then
         if (iZk(iK).eq.iZ.and.iAk(iK).eq.iA) then
           barrier = fissbark(isoindex)
           RETURN
         endif
       endif
!
       foundZ=.FALSE.
       do isoindex=1,nisotopes
         if (iZ.eq.Zk(isoindex)) then
           foundZ=.TRUE.
           if (iA.eq.Ak(isoindex)) then
             barrier = fissbark(isoindex)
             RETURN
           endif
         endif
       enddo ! loop over all ZA's
       if(foundZ.eqv..TRUE.) then
#ifdef WRITEL6
         write (L6,*) & 
           '**',iZ,iA,': No fission barrier for this isotope!'
#endif
         barrier = 0.0
       else
#ifdef WRITEL6
         write (L6,*) '**',iZ,': No fission barriers for this element!'
#endif
         barrier = 0.0
       endif
       RETURN
  901  format (a70)
       END

!************************************************************************
       FUNCTION alevel (iK,iZ0,iA0,eps,U)   ! Level density parameter 
! based on calculations by [Koura042] provided by RV 11sep0.
! but using the back-shifted U(eps).
! OBS: PM has reservations about the shell effects by Koura!
! -----------------------------------------------------------------------
! Note that the data file obtained from the RIPL-3 people is inconsistent
! with the parameters in the RIPL-3 paper, on which the fits of parameters
! alpha, beta, gamma, delta are based.
! -----------------------------------------------------------------------
! iK    index identifying the isotope considered (iK=1,...,mxK);
! iZ0   charge number of the nucleus
! iA0   mass number of the nucleus
! eps   excitation energy of the nucleus (MeV)
! e     shape eccentricity (not used here)
! U     effective excitation U=eps-Delta, splined smoothly to zero at low eps
! The data file 'alevel.dat' must contain {iZ,iA,deltaW,Delta} in order of
! uninterrupted increasing iZ, and for each iZ uninterrupted increasing iA.
!--------------------------------------------------------------------------
       use freyaout
       use freyaparams
       use freyaparameters
       use freyaPath
       use freyaK
       use freyaerrors

       implicit none

       real alevel
       integer iK,iZ0,iA0
       real eps,U

       character*100 dataf
       logical exists

       integer MAX
       parameter (MAX=7076)                     !Koura: 7050 + 26 additions
       integer,save:: iZn,iAn,noZ,nZ,nA
       real,save:: deltaW,Delta
       dimension iZn(MAX),iAn(MAX),deltaW(MAX),Delta(MAX)
       dimension noZ(0:120),nZ(121),nA(121)

       real gamma
       data gamma/0.05/
       integer n0,nm,np
       data n0,nm,np/0,0,0/                     ! count backshifts
       character*70 headline                    ! Headline of data file

       integer iZ,iA,n,i,kA,kZ,m
       real dW,D,atilde,aA,damp

  901  format (a70)
  902  format ('Headline: ',a70)

!      -----------------------------------------------------------------
       alevel = 0.0             ! default initialization

       if (iA0.lt.0) then       ! finish up:
#ifdef WRITEL6
         do iZ=1,120
           if (noZ(iZ).gt.0) write (L6,*) 'alevel:', &
             ' missing element Z =',iZ,' called',noZ(iZ),' times!'
         enddo
         if (noZ(0).gt.0) write (L6,*) 'alevel:', &
             noZ(0),' missing isotopes in total!'
         write (L6,"('alevel: Backshift counts:',3i8)") n0,nm,np
#endif
         RETURN
       endif
!      -----------------------------------------------------------------
       IF (iA0.eq.0) THEN       ! INITIALIZE:
!        nZ(iZ) points to location in table of first isotope of iZ
!        nA(iZ) is the number of (consequtive) isotopes of iZ
         do iZ=0,120
           noZ(iZ)=0; nZ(iZ)=0
         enddo
         dataf=trim(freyadir)//"/alevel.dat"
         inquire(FILE=dataf, EXIST=exists)
         if (exists.eqv..FALSE.) then
#ifdef WRITEL6
           write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
           write(error, 13) trim(dataf)
13         format("FREYA: data file ", a, " could not be found.")
           call seterror(6,error)
           RETURN
         end if
         OPEN  (L7,file=dataf,status='old')
         read  (L7,901) headline
#ifdef WRITEL6
         write (L6,902) headline
#endif
         n=0
   10    read (L7,*) iz, ia, dW, D
         if (iz.ne.0) then
           if (n.ge.MAX) then
#ifdef WRITEL6
             write (L6,*) '** File  alevel.dat  is too large!'
#endif
             do while (iz.ne.0)
               n=n+1
               read (L7,*) iz, ia, dW, D
             enddo
#ifdef WRITEL6
             write (L6,*) '** Suggestion: Increase MAX to',n,'!'
             write (L6,*) '** Continuing with MAX =',MAX,'!'
#endif
             n=MAX
           else
             n=n+1              ! iza=1000*iz+ia =>
             iZn(n)=iz          ! iza/1000
             iAn(n)=ia          ! MOD(iza,1000)
             deltaW(n)=dW       ! Shell correction to the mass
             Delta(n)=D         ! Pairing energy
             goto 10
           endif
         endif
         CLOSE (L7)
#ifdef WRITEL6
         write (L6,*) 'alevel:',n,' entries in file alevel.dat'
#endif
! Done reading data file - now take a look at it:
         do i=1,121
           nZ(i)=0              ! will point to start of charge Z
         enddo
         kZ=0
         do i=1,n
           if (iZn(i).lt.kZ) then
#ifdef WRITEL6
             write (L6,*) '** Z out of order in alevel.dat!'
             write (L6,*) '** ', i, iZn(i), kZ
             write (L6,*) '** Aborting!'
#endif
             STOP
           else if (iZn(i).gt.kZ) then
             kZ=iZn(i)          ! Next charge in the table
             nZ(kZ)=i           ! Its location in table
             nA(kZ)=1
             kA=iAn(i)
!             write (L6,*) i,kZ,kA
           else if (iAn(i).le.kA) then
#ifdef WRITEL6
             write (L6,*) '** A out of order in alevel.dat!'
             write (L6,*) '** ', i, iZn(i), iAn(i), kA
             write (L6,*) '** Aborting!'
#endif
             STOP
           else if (iAn(i).gt.kA+1) then
#ifdef WRITEL6
             write (L6,*) '** Missing A in alevel.dat!'
             write (L6,*) '** ', i, iZn(i), iAn(i), kA
             write (L6,*) '** Aborting!'
#endif
             STOP
           else
             kA=iAn(i); nA(kZ)=nA(kZ)+1
           endif
         enddo
         nZ(iZn(n)+1)=n+1; nA(iZn(n)+1)=0
         alevel = 0.0
         RETURN
       ENDIF                              ! INITIALIZE.
!      -----------------------------------------------------------------
!
! Find value of alevel0 (form inputparameters.dat)
!
       alevel0=alevelk(iK)
 
!
! CALCULATE the level density parameter:
       n=nZ(iZ0)                         ! jump to the start of this Z=iZ0
       m=n+iA0-iAn(n)
! Check if the specified iZ0 is in the table:
       if (n.eq.0) then
#ifdef WRITEL6
         if (noZ(iZ0).eq.0) write (L6,*) &
            '** Z = ',iZ0,' is not in alevel.dat!'
!     c      ,'  Using a = A/(',alevel0,' MeV).'
#endif
         alevel = float(iA0)/alevel0
         noZ(iZ0)=noZ(iZ0)+1
         RETURN
       endif
! Check if the specified A is in the table for this Z:
       if (iA0.lt.iAn(n).OR.iA0.ge.iAn(n)+nA(iZ0)) then
#ifdef WRITEL6
         write (L6,*) n,':',iAn(n),iA0,iAn(n)+nA(iZ0)-1,iZ0
         write (L6,*) &
           '** (Z,A) = (',iZ0,',',iA0,') is not in alevel.dat!'
!     c    ,'  Using a = A/(',alevel0,' MeV).'
#endif
         alevel = float(iA0)/alevel0; ! pause!!
         noZ(0)=noZ(0)+1
         RETURN
       else    ! calculate the location of iA0:
         n=n+iA0-iAn(n)
       endif
       dW=deltaW(n)                             ! Shell
       D=Delta(n)                               ! Pairing
! Calculate U(eps) by splining for low eps:
       U=eps                                    ! if Delta=0
       if (D.gt.0.0) then                       ! if Delta>0:
         if (eps.ge.2*D) then
           U=eps-D; n0=n0+1
         else
           U=0.25*eps**2/D; np=np+1
         endif
       else                                     ! if Delta<0:
         if (eps.ge.abs(D)) then
           U=eps-D; n0=n0+1
         else
           U=2*sqrt(eps*abs(D)); nm=nm+1
         endif
       endif
! Cancel:      U=eps!!! for RV 27feb11
! Calculate alevel:
       atilde=float(iA0)/alevel0                ! atilde
       aA=atilde*(1.0+gamma*dW)                 ! Use this for U=0
       if (aA.le.0.0) then                      ! happens if dW < -gamma
#ifdef WRITEL6
         write (L6,*) 'a<0:',iZ0,iA0,eps,aA,dW
         write (L6,*) n,iZn(n),iAn(n),gamma
         write (L6,*) Delta(n),deltaW(n),atilde
#endif
         ! pause
       endif
       damp= 1.0-exp(-gamma*U)                  ! 1-exp(-gamma*U)
       if (U.gt.0.0) aA = atilde*(1+damp*dW/U)  ! Use this for U > 0
!      write (L6,101) iZ0,iA0,eps, aA,atilde,dW,D,g
!  101  format ('##',2i5,f10.5,':',5f10.5)
       alevel = aA
!c     alevel=atilde    ! no shell effects!
       RETURN
       END

!************************************************************************
       FUNCTION Gammanf (iK,iZ,iA,eps,Bn,Bf,Gn,Gf)
! calculates the widths for neutron evaporation and fission:
!  ->    iK:       index identifying the isotope considered (iK=1,...,mxK);
!  ->    iZ,iA:    charge & mass numbers of the decaying nucleus;
!  ->    eps:      its excitation energy;
!  ->    Bn & Bf:  barriers for neutron evaporation and fission.
! <-    Gammanf:   2*pi*rho0 = 2*pi*exp(2*sqrt(a0*eps))
! <-    Gn & Gf:   partial widths for neutron evaporation & fission
!------------------------------------------------------------------------
       use freyaout
       use freyaG
       use freyaconsts
       use freyaPath
       use freyaK
       use freyaerrors

       implicit none

       real Gammanf
       integer iK,iZ,iA
       real eps,Bn,Bf,Gn,Gf

       real c1,c2,c3,r0
       data c1,c2,c3/0.04543,0.1355,0.1426/
       data r0/1.15/                ! differs from the Lysekil value 1.2 used elsewhere!
       logical, save:: dataread=.FALSE.
       character*100 dataf
       logical exists

       integer,save,dimension(:),allocatable::iKtoCorrk
                                         ! mapping between iK and index
                                         ! of element in An_corr/Af_corr
                                         ! arrays
       integer,save,dimension(:),allocatable::Zk 
                                         ! array of the isotope Z's
       integer,save,dimension(:),allocatable::Ak 
                                         ! array of the isotope A's
       real,save,dimension(:),allocatable::An_corr
                                         ! array of correction factors An
       real,save,dimension(:),allocatable::Af_corr
                                         ! array of correction factors Af

       integer,save::nisotopes           ! number of An/Af correction
                                         ! factors in file 'acorrection.dat'
       real anco,afco,ancorr,afcorr,R,a0,alevel,U0,rho0,Xf &
       ,af,an,dx,x,xeff,Uf,Un,Xn
       character*300 headline            ! Headline of data file 
       character*300 line                ! line of data
       integer nZ,nA,n,n2,i,iost,isoindex,iZA

       Gammanf=0.                        ! default initialization

       if (dataread.eqv..FALSE.) then
! read in file acorrection.dat, if it exists
         dataf=trim(freyadir)//"/acorrection.dat"
         inquire(FILE=dataf, EXIST=exists)
         if (exists.eqv..FALSE.) then
#ifdef WRITEL6
           write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
           write(error, 13) trim(dataf)
13         format("FREYA: data file ", a, " could not be found.")
           call seterror(10,error)
           RETURN
         end if
         OPEN (L7,file=dataf,status='old')
         read  (L7,*,IOSTAT=iost) headline
         isoindex=1
         do
           read (L7,'(a300)',IOSTAT=iost) line
           if (iost.ne.0) then
             ! something went wrong reading file acorrection.dat
             ! or reached end of file
             exit
           end if
           isoindex=isoindex+1
         enddo
         CLOSE (L7)
         nisotopes = isoindex-1        
! allocate memory to iKtoCorrk, Zk, Ak, An_corr and Af_corr
         if (.not.allocated(iKtoCorrk)) then  
           allocate(iKtoCorrk(mxK))
            iKtoCorrk=0
         endif
         if (.not.allocated(Zk)) then  
           allocate(Zk(nisotopes))
            Zk=0
         endif
         if (.not.allocated(Ak)) then
           allocate(Ak(nisotopes))
            Ak=0
         endif
         if (.not.allocated(An_corr)) then
           allocate(An_corr(nisotopes))
            An_corr=0
         endif
         if (.not.allocated(Af_corr)) then
           allocate(Af_corr(nisotopes))
            Af_corr=0
         endif
!
         OPEN (L7,file=dataf,status='old')
         read  (L7,*,IOSTAT=iost) headline
         do isoindex=1,nisotopes
           read (L7,*,IOSTAT=iost) nZ,nA,Anco,Afco
           if (iost.ne.0) then
             ! something went wrong reading file acorrection.dat
             ! or reached end of file
             exit
           end if
           ! The line was read properly
           Zk(isoindex) = nZ
           Ak(isoindex) = nA
           An_corr(isoindex)=Anco
           Af_corr(isoindex)=Afco
           ! find mapping from iK to index of An_corr/Af_corr arrays
           do iZA=1,mxK
             if (nZ.eq.iZk(iZA).and.nA.eq.iAk(iZA)) then
               iKtoCorrk(iZA)=isoindex
             end if
           end do
         end do ! loop over entire file 'acorrection.dat'
         dataread=.TRUE.
         CLOSE (L7)
       endif ! dataread == false
!
! Look for the an and af correction factors in array An_corr and Af_corr
!
! Set default values
       ancorr = 1.00
       afcorr = 1.00
!
! See if we have this iK mapped onto one of the elements in
! the array An_corr/Af_corr
       isoindex=iKtoCorrk(iK)
       if (isoindex.ne.0.and.iZk(iK).eq.iZ.and.iAk(iK).eq.iA) then
         ancorr = An_corr(isoindex)
         afcorr = Af_corr(isoindex)
       else
         do isoindex=1,nisotopes
           if (iZ.eq.Zk(isoindex)) then
             if (iA.eq.Ak(isoindex)) then
               ancorr = An_corr(isoindex)
               afcorr = Af_corr(isoindex)
               exit
             endif
           endif
         enddo ! loop over all ZA's
       endif
!

! Since the shape effect on the level-density parameter is very small,
! we shall neglect it.
!      write (L6,"('Gammanf:'2i4,3f10.3)") iZ,iA,eps,Bn,Bf

       R=r0*float(iA)**0.333333       ! Nuclear radius
!      a0=((c1*R+c2)*R+c3)*R          ! We put Bsurf & Bcurv to unity
       a0=alevel(iK,iZ,iA,eps,U0)
!       write (L6,*) iZ,iA,eps,': A/a0=',iA/a0,' MeV'
       rho0=twopi*exp(2*sqrt(a0*eps))
       n=200; n2=n*2
! Calculate the fission width:    -----------------------------------------
       Gf=0.0; Xf=eps-Bf              ! Max of excitation energy x at barrier top
       if (Xf.gt.0.0) then            ! integrate over excitation energy x:
         af=alevel(iK,iZ,iA,Xf,Uf)    ! at fission barrier
! Change af for fission of U-Nth:
         af = af*afcorr
         dx=Xf/n2                     ! Step size in exc energy: x=(1,3,5,..)*dx
         do i=1,n2,2
           x=i*dx; xeff=x
           Gf=Gf+exp(2*sqrt(af*xeff))
!          Generally af should be af(x) in the integrand!
         enddo
         Gf=2*dx*Gf/rho0
       else
         Gf=0.0
       endif
! Calculate the evaporation width:-----------------------------------------
       Gn=0.0; Xn=eps-Bn              ! Max excitation energy in n-evap daughter
       if (Xn.gt.0.0) then            ! integrate over excitation energy x:
         an=alevel(iK,iZ,iA-1,Xn,Un)  ! after evaporation
! Change an for evaporation from Pu-Nth:
         an = an*ancorr
         dx=Xn/n2                     ! Step size in exc energy: x=(1,3,5,..)*dx
         do i=1,n2,2
           x=i*dx; xeff=x
           Gn=Gn+(Xn-x)*exp(2*sqrt(an*xeff))
!          Generally an should be an(x) in the integrand!
         enddo
         Gn=2*dx*Gn*Gnf/rho0
!c       if (iA.eq.235) Gn=10.0*Gn!!
       else
         Gn=0.0
       endif
       Gammanf = rho0

       RETURN
       END

!************************************************************************
       SUBROUTINE SCISSION (iK)
! tabulates the Q-value Qsciss(Nlight,Zlight,Nth) for fission of species iK
! at the excitation E*=eps0; the dependence on eps0 is ignored, except that
! eps0=0 gives spontaneous fission.
! This routine adjusts the scission configurations to Ktot(Af) data.

!------------------------------------------------------------------------
       use freyaout
       use freyaK
       use freyaMth
       use freyaB
       use freyaparams
       use freyaZ
       use freyaC
       use freyaLD
       use freyacSission
       use freyaPath
       use freyaerrors

       implicit none

       integer iK

       integer minN1, minZ1

       integer,dimension(:,:),allocatable:: nc ! Number of channels open for Nth
       real Edata(300,0:3)                     ! Data on TKE(Af) [n>0: set n; 0: ave]

       character*80 headline                   ! Headline of data file	
       character*100 dataf
       logical exists
       integer iZ0,iZ,iA,iAf,iA0,iA1,iA2,iAhalf,l,m, &
               i1,i2,jA1,jA2,i, &
               Nth,iN1,iZ1,minA1,maxA1,iZa,iZb,iZ2, &
               iN2,minTKE,maxTKE

       integer iTKEfidx                 ! maps index iK to index of 
                                        ! file (Ein,dTKE)
       integer nergs                    ! number of energies in file
                                        ! (Ein,dTKE)
       real Af,Etot,Eave,D0,ZoverA,dV1,dV2, &
            gsM,Edat,Ec,Z,Z1,Z2,s,D2,Q12, &
            D1,Q

       data minN1,minZ1/8,8/

       if (.not.allocated(nc)) then
         allocate(nc(0:Mxth,mxK))
         nc=0
       endif

       iZ0=iZk(iK); iA0=iAk(iK)
#ifdef WRITEL6
       write (L6,"(80('-'))")
       write (L6,"('Setting up SCISSION for iK=',i3,':',2i5)") &
              iK,iZ0,iA0
#endif

!========================================================================
! Copy over the Ktot(Af) data for the specified specie (iZ0,iA0):
       do iA=1,300
         do l=0,3
           Edata(iA,l)=0.0
         enddo
       enddo
! TKE(Af) data for ZZAAA:       -----------------------------------------
#ifdef WRITEL6
       write (L6,"('TKE(Af) data used for',i4,i5,':')") iZ0,iA0
#endif
       do l=1,ntkek(iK)
         dataf=trim(freyadir)//'/'//trim(tkek(iK,l))
         inquire(FILE=dataf, EXIST=exists)
         if ((exists.eqv..FALSE.).and.(trim(tkek(iK,l)).ne.'-')) then
#ifdef WRITEL6
           write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
           write(error, 13) trim(dataf)
13         format("FREYA: data file ", a, " could not be found. ", &
                  "Check react.dat file in FREYA data directory.")
           call seterrorcase(4,error)
           RETURN
         endif
         OPEN  (L7,file=dataf,status='old')
         read  (L7,901) headline ! Data from Nishio
#ifdef WRITEL6
         write (L6,902) headline
#endif
         read  (L7,*) Af,Etot
         do while (Af.ne.0.0)
           iAf=int(Af+0.5); Edata(iAf,l)=Etot
!          write (L6,'(i5,f10.5)') iAf,Etot
           read  (L7,*) Af,Etot
         enddo
         CLOSE (L7)
       enddo
!------------------------------------------------------------------------
! Calculate the average of the experimental TKE values for each Af,
! counting each data point equally (so Af & A0-Af count equally):
#ifdef WRITEL8
       write (L8,"(80('-'))"); write (L8,"('SCISSION: TKE(Af)')")
#endif
       minTKE=0         ! Smallest mass for which TKE data exist
       maxTKE=0         ! Largest mass for which TKE data exist
       iAhalf=iA0/2
       do iA1=1,iAhalf  ! Fragment Af
         iA2=iA0-iA1    ! Af's partner fragment
         m=0; Eave=0.0
         do l=1,ntkek(iK)
           if (Edata(iA1,l).NE.0.0) then
             m=m+1
             Eave=((m-1)*Eave+Edata(iA1,l))/m
           endif
           if (iA1.NE.iA2.AND.Edata(iA2,l).NE.0.0) then
             m=m+1
             Eave=((m-1)*Eave+Edata(iA2,l))/m
           endif
         enddo
         Edata(iA1,0)=Eave; Edata(iA2,0)=Eave
         if (Eave.eq.0.0.AND.maxTKE.eq.0) minTKE=iA1+1
         if (Eave.gt.0.0.AND.minTKE.gt.0) maxTKE=iA1
#ifdef WRITEL8
         write (L8,"(i5,4f10.4)") iA1,(Edata(iA1,l),l=0,3)
#endif
       enddo
       maxTKE=iA0-minTKE
#ifdef WRITEL8
       do iA1=iAhalf+1,maxTKE
         write (L8,"(i5,4f10.4)") iA1,(Edata(iA1,l),l=0,3)
       enddo
#endif
!------------------------------------------------------------------------
! Construct TKE data by interpolation for any missing fragment masses:
#ifdef WRITEL8
       write (L8,"('SCISSION: TKE(Af):',2i5)") minTKE,maxTKE
#endif
       iA1=minTKE
       do while (iA1.lt.maxTKE)! -> iA1: lower interpolation point, E(iA1)>0
#ifdef WRITEL8
         write (L8,"(i5,':',f10.3)") iA1,Edata(iA1,0)
#endif
         if (Edata(iA1+1,0).eq.0) then
           iA2=iA1+2           ! -> iA2: upper interpolation point, E(iA2)>0
           do while (Edata(iA2,0).eq.0.0)
             iA2=iA2+1
           enddo
           do iA=iA1+1,iA2-1
             Eave= &
              ((iA2-iA)*Edata(iA1,0)+(iA-iA1)*Edata(iA2,0))/(iA2-iA1)
#ifdef WRITEL8
             write (L8,"(i5,':',2f10.3)") iA,Eave,Edata(iA,0)
#endif
             Edata(iA,0)=Eave
           enddo
           iA1=iA2-1
         endif
         iA1=iA1+1
       enddo
#ifdef WRITEL8
       write (L8,"(i5,':',f10.3)") iA1,Edata(iA1,0)
#endif
       call FLUSH (L8)

!========================================================================
       D0=gsM(iZ0,iA0)          ! Mother mass excess
       ZoverA=float(iZ0)/iA0    ! Overall charge-to-mass ratio
       dV1=0.0; dV2=0.0         ! Distortion energies at scission: none here!
#ifdef WRITEL8
       write (L8,*) 'SCISSION: dTKE =',dTKE
#endif
!------------------------------------------------------------------------
! Adjust TKE(Af) to data (use same center separation for all Zf for each Af)
       minA1=minN1+minZ1
       maxA1=iA0/2
       jA1=0; jA2=iA0-jA1
       do iA1=minA1,maxA1
         iA2=iA0-iA1
! Calculate the average experimental kinetic energy for this iA1, TKE(Af):
         Edat=Edata(iA1,0)       ! Exp TKE for this split
         if (Edat.eq.0.0) goto 29! No TKE data -> skip this split!
         if (dTKEe(iK).eqv..FALSE.) then
           nergs = 1
         else
! Compute the quantities s12sciss, Esciss and Qsciss for all the
! different incoming neutron energies
           nergs = nEin(iKtodTKEfk(iK))
         endif
         do l=1,nergs
           iTKEfidx = 0
           if (dTKEe(iK).eqv..FALSE.) then
             dTKE=dTKE0k(iK)
           else
             iTKEfidx = iKtodTKEfk(iK)
             dTKE = dTKEk(iTKEfidx,l)
           endif
           Eave=Edat+dTKE          ! to match both <n> and <TKE>!!
           do Nth=0,Mthk(iK)
             Ec=Eave
             Z1=iA1*iZ0/float(iA0-Nth)
             Z2=iZ0-Z1; iA2=iA0-Nth-iA1
             s=esq*Z1*Z2/Ec               ! Center separation that matches Ec
             iZ1=int(Z1+0.5); iZa=iZ1-6; iZb=iZ1+6
             do iZ1=iZa,iZb
               iN1=iA1-iZ1
               iZ2=iZ0-iZ1; iN2=iA2-iZ2
! Calculate the Q value for this split:
               D1=gsM(iZ1,iA1)                   ! Fragm #1 mass excess
               D2=gsM(iZ2,iA2)                   ! Fragm #2 mass excess
               Q12=D0-D1-D2                      ! sf Q-value [0 -> 1+2]
               Qsf     (iN1,iZ1,Nth,iK) = Q12
               Q1sciss (iN1,iZ1,Nth,iK) = dV1    ! =0 here
               Q2sciss (iN1,iZ1,Nth,iK) = dV2    ! =0 here
               Q=Q12-Ec-dV1-dV2                  ! Internal heat in 1 & 2
               Q=max(0.0,Q)
! Q>0  : Channel is open
! Q==0 : Channel is closed
               if (dTKEe(iK).eqv..FALSE.) then
                 ! s12sciss(iN1,iZ1,Nth,iK) = s
                 Esciss (iN1,iZ1,Nth,iK) = esq*iZ1*iZ2/s
#if defined(WRITEL6) || defined(WRITEL8)
                 Qsciss(iN1,iZ1,Nth,iK)=Q
#endif
               else
                 ! s12scissk(iN1,iZ1,Nth,iTKEfidx,l) = s
                 Escissk (iN1,iZ1,Nth,iTKEfidx,l) = esq*iZ1*iZ2/s
#if defined(WRITEL6) || defined(WRITEL8)
                 Qscissk(iN1,iZ1,Nth,iTKEfidx,l)=Q
#endif
               endif
               if (Q.gt.0.0) then
                 nc(Nth,iK)=nc(Nth,iK)+1         ! # channels open
               endif
#ifdef WRITEL6
               if (iA1.eq.jA1.AND.iA2.eq.jA2) write (L6,"(2i4,5f10.4)") & 
                 iZ1,iZ2,Q12,esq*iZ1*iZ2/s,Qsciss(iN1,iZ1,Nth,iK)
#endif
             enddo
           enddo
         enddo
   29    if (iA1.eq.0.AND.iA2.eq.iA0-0) then
           ! pause!!
         endif
       enddo

!------------------------------------------------------------------------
       continue
#ifdef WRITEL8
       write (L8,*) 'Tables of Qsciss(Zf) for given Af in [100 keV] for'
       write (L8,*) & 
         'Zf = [Z1ave]-6,..., [Z1ave]+6; Zfave = (Z0/A0)*Af ..'
       do Nth=0,Mthk(iK)
         write (L8,*) 'Qsciss for Nth =',Nth,':'
         write (L8,*) ' A1 Z1ave  Z1min                  [Z1ave]' &
           ,'                  Z1max'
         do iA1=60,iA0/2
           Z=ZoverA*iA1; iZ=int(Z+0.5)
           i1=iZ-6; i2=iZ+6
           write (L8,"(i4,f6.2,1x,17i4)") &
             iA1,Z,(int(10*Qsciss(iA1-i,i,Nth,iK)+.5),i=i1,i2)
         enddo
         write (L8,*) nc(0,iK),' fission channels are open for Nth =', &
                      Nth
       enddo

       write (L8,*) 'SCISSION: Number of open fission channels:'
       do Nth=0,Mthk(iK)
         write (L8,*) 'Nth =',Nth,':',nc(Nth,iK)
       enddo
#endif

#ifdef WRITEL6
       write (L6,*) 'SCISSION: Number of open fission channels:'
       do Nth=0,Mthk(iK)
         write (L6,*) 'Nth =',Nth,':',nc(Nth,iK),Bf(Nth,iK),Sn(Nth,iK)
       enddo
#endif

       RETURN
  901  format (a80)
  902  format (a80)
       END

!************************************************************************
       SUBROUTINE MASSES (iK,eps0)
! prepares the fission-fragment mass distribution for specie iK:
       use freyaparameters
       use freyaout
       use freyaK
       use freyaMth
       use freyaZ
       use freyaWg
       use freyaSAMPLE
       use freyacYAf
       use freyaPath
       use freyaerrors

       implicit none

       integer iK
       real eps0

       logical msFREYA_errorflagset
       logical exitonerror

       character*70 char70

       logical, save :: dataread=.FALSE.
       character*100 dataf
       logical exists
       character*100 headline                  ! Headline of data file	

       real,dimension(:,:),allocatable:: CC0, & 
            CCm, & 
            CCp                                ! Auxiliary arrays
       real,dimension(:,:,:),allocatable:: CCg ! Bin fission modes n=0:Mxg for Nth=0:Mth
       real,dimension(:,:),allocatable:: Cf0   ! Number of fissions with given Nth=0:Mth.
       real,dimension(:,:),allocatable:: Cf1   ! Bin fission energies for Nth=0:Mth;
       real,dimension(:,:),allocatable:: Cf2   ! bin the corresponding squares

       integer iZ0, iA0, iA
       integer m, n1, n2, n, i, Nth, iost, ind
       real s, Y, zdisv
       integer minAf, maxAf
       logical exist

       if (.not.allocated(CC0)) then
         allocate(CC0(0:Mxg,mxK))
         CC0=0
       endif
       if (.not.allocated(CCm)) then
         allocate(CCm(0:Mxg,mxK))
         CCm=0
       endif
       if (.not.allocated(CCp)) then
         allocate(CCp(0:Mxg,mxK))
         CCp=0
       endif
       if (.not.allocated(CCg)) then
         allocate(CCg(0:Mxg,0:Mxth,mxK))
         CCg=0
       endif
       if (.not.allocated(Cf0)) then
         allocate(Cf0(0:Mxg,mxK))
         Cf0=0
       endif
       if (.not.allocated(Cf1)) then
         allocate(Cf1(0:Mxg,mxK))
         Cf1=0
       endif
       if (.not.allocated(Cf2)) then
         allocate(Cf2(0:Mxg,mxK))
         Cf2=0
       endif

       iZ0=iZk(iK); iA0=iAk(iK)
#ifdef WRITEL6
       write (L6,"(80('-'))"); write (L8,"(80('-'))")
       write (L6,"('Preparing P(Af) for #',i2,':',i4,i5)") iK, iZ0,iA0
#endif
       sampleAf(iK)=.FALSE.             ! Use 5-g fit
       exist = .FALSE.
        ! If we find a file *-Af*.dat, sample it
       dataf=trim(freyadir)//'/'//trim(pAfk(iK))
       inquire(FILE=dataf, EXIST=exists)
       if ((exists.eqv..FALSE.).and.(trim(pAfk(iK)).ne.'-')) then
#ifdef WRITEL6
         write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
         write(error, 13) trim(dataf)
13       format("FREYA: data file ", a, " could not be found.", &
                "Check react.dat file in FREYA data directory.")
         call seterrorcase(3,error)
         RETURN 
       end if
       sampleAf(iK)=exists
       IF (sampleAf(iK)) THEN           ! sample P(Af) directly:
#ifdef WRITEL6
         write (L6,*) 'Sampling directly from P(Af):'
#endif
         OPEN (L7,file=dataf,status='old')
         read (L7,'(a70)') char70
#ifdef WRITEL6
         write (L6,'(a70)') char70
#endif
         read (L7,*) iA,Y; minAfk(iK)=iA! Minimum Af in data table
         do while (iA.ne.0)
           YYAfk(iA,iK)=Y
           maxAfk(iK)=iA                ! Maximum Af in data table
           read (L7,*) iA,Y
         enddo
         CLOSE (L7)
! Symmetrize Y(Af):
         minAf=minAfk(iK)
         maxAf=maxAfk(iK)
         m=(iA0+1)/2; n1=m-minAf; n2=maxAf-m
         n=min(n1,n2); minAf=m-n; maxAf=m+n
#ifdef WRITEL6
         write (L6,"('minAf=',f10.5,'; maxAf=',2f10.5,i4)") &
               1.*minAf,1.*maxAf,1.*(minAf+maxAf),iA0
#endif
         s=0.0
         do i=minAfk(iK),maxAfk(iK)
                YAfk(i,iK)=0.5*(YYAfk(i,iK)+YYAfk(maxAf-i+minAf,iK))
!               write (6,"(i3,3f10.5)")
!      c      i,YAfk(i,iK),YYAfk(i,iK),YYAfk(maxAf-i+minAf,iK)
           s=s+YYAfk(i,iK)
         enddo
!NOTE: The sampling of Af can be sped up by a factor of two by sampling only
!      the heavy (or light) fragment!
! Normalize Y(Af) to TWICE times unity (to allow selection of Alight):
         YYAfk(minAf-1,iK)=0.0
         do iA=minAfk(iK),maxAfk(iK)
           YAfk(iA,iK)=2*YAfk(iA,iK)/s
           YYAfk(iA,iK)=YYAfk(iA-1,iK)+YAfk(iA,iK)
!           write (L6,'(i5,2f10.5)') iA,YAfk(iA,iK),YYAfk(iA,iK)
         enddo
         YYAfk(maxAf,iK)=2.0
       ELSE                             ! sample from 5-gaussian fit:
         call msFREYA_reseterrorflag()
         call SplitFits (iK,eps0)
         if (msFREYA_errorflagset().and.exitonerror()) RETURN
       ENDIF
       do Nth=0,Mthk(iK)
         do n=0,Mxg
           CCg(n,Nth,iK)=0.0            ! Count various fission modes n
         enddo
         Cf0(Nth,iK)=0.0
         Cf1(Nth,iK)=0.0
         Cf2(Nth,iK)=0.0
       enddo

! Charge distribution:
!      ZDis=sqrt(Zdis0**2 + 1./12)     ! Why 1/12?
! The values of Zdis0 are from Reisdorf et al, NPA 177 (1971) 337.
       if (dataread.eqv..FALSE.) then
! read in file Zdis.dat, if it exists
         dataf=trim(freyadir)//"/Zdis.dat"
         inquire(FILE=dataf, EXIST=exists)
         if (exists.eqv..FALSE.) then
#ifdef WRITEL6
           write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
           write(error, 13) trim(dataf)
           call seterror(11,error)
           RETURN 
         end if
         OPEN (L7,file=dataf,status='old')
         read (L7,*,IOSTAT=iost) headline
         if (iost.eq.0) then
           do 
             read  (L7,*,IOSTAT=iost) iA,zdisv
             if (iost.ne.0) then
               ! something went wrong reading file
               ! Zdis.dat, or reached end of file
               exit
             endif 
             ! The line was read properly
             do ind=1,mxK
               if (iA.eq.iAk(ind)) then
                 Zdis(ind)=zdisv
               end if
             end do
           end do ! loop over entire file 'Zdis.dat'
           dataread=.TRUE.
         end if
         CLOSE (L7)
       endif ! dataread == false
       RETURN
       END

!************************************************************************
       SUBROUTINE SplitFits (iK,epsmx)
!-----------------------------------------------------------------------*
! This routine tabulates the multi-gaussian fits to the fission-fragment*
! mass distribution P(Af) as a function of the excitation energy for the*
! original compound nucleus (Z00,A00) and up to Mth neutron-evaporation *
! daughters; the tabulation is done for E*(iE)=(iE-0.5)*dEf, iE=1,..,mE,*
! with mE=2+epsmx/dEf to ensure that E*(mx)>epsmx.                      *
!-----------------------------------------------------------------------*
! iZ0   Charge number           !>
! iA0   Mass number             | >    of the initial compound nucleus
! epsmx Max excitation energy   !>
! Mth   Maximum number of prefission neutron emissions: Nth=0,...,Mth

!-----------------------------------------------------------------------*
       use freyaparameters
       use freyaout
       use freyaK
       use freyaMth
       use freyaWg
       use freyaPath
       use freyaerrors

       use freyainterfaces

       implicit none

       integer iK
       real epsmx

       integer mode, m0, i0            ! for plot of P(Af) in AfDist
       real E0

       real,save,dimension(:,:,:),allocatable::Cn! -> Weight of each fission mode n=0:Mxg
       real,save,dimension(:,:),allocatable::Dn  ! Centroid shift of fission mode n=0:Mxg
       real,save,dimension(:,:,:),allocatable::Wn! -> Width of Gaussian for mode n=0:Mxg

       logical, save:: dataread=.FALSE.
       character*100 dataf
       logical exists
       character*70 headline           ! Headline of data file
       character*300 line,line2        ! line in data file
       integer iZ,iA,i,iZAindex,iline,iZA,n,iost,iZ0,iA0,iE &
       ,Mth,mE,m
       real sqrt,Tm,aA,En,eps,U,Sn,epsAve,epsm,msFREYA_SEPn
       logical eof
       logical foundZA
!-----------------------------------------------------------------------*
       if (.not.allocated(Cn)) then
         allocate(Cn(0:Mxg,0:2,mxK))
         Cn=0
       endif
       if (.not.allocated(Dn)) then
         allocate(Dn(0:Mxg,mxK))
         Dn=0
       endif
       if (.not.allocated(Wn)) then
         allocate(Wn(0:Mxg,0:2,mxK))
         Wn=0
       endif

       if (dataread.eqv..FALSE.) then
! read in file gaussfit.dat, if it exists
         dataf=trim(freyadir)//"/gaussfit.dat"
         inquire(FILE=dataf, EXIST=exists)
         if (exists.eqv..FALSE.) then
#ifdef WRITEL6
           write (L6,*) 'File ',trim(dataf),' does not exist'
#endif
           write(error, 13) trim(dataf)
13         format("FREYA: data file ", a, " could not be found.")
           call seterror(12,error)
           RETURN
         end if
         eof = .FALSE.
         OPEN (L7,file=dataf,status='old')
         do i=1,7
           read  (L7,*,IOSTAT=iost) headline
           if (iost.ne.0) then
             ! reached end of file
             eof=.TRUE.
             exit
           end if
         end do
         if (eof.neqv..TRUE.) then
           do
             foundZA = .FALSE.
             do iline=1,6
               read (L7,'(a300)',IOSTAT=iost) line
               if (iost.ne.0) then
                 ! something went wrong reading file gaussfit.dat
                 ! or reached end of file
                 eof=.TRUE.
                 exit
               end if
               if (iline.eq.1) then
                 ! Jerome, I do not quite like this
                 read (line,'(i2,x,i3,a)') iZ,iA,line2
                 ! read (line,'(2i,a)') iZ,iA,line2
                 ! Find Z,A in the array iZk,iAk
                 do iZA=1,mxK
                   if (iZ.eq.iZk(iZA).and.iA.eq.iAk(iZA)) then
                     foundZA = .TRUE.
                     line = line2
                     exit
                   end if
                 end do
               endif
               ! if (iZAindex.eq.-1) then
               if (foundZA.eqv..FALSE.) then
                 ! ZZAAA does not exist in file react.dat, 
                 ! skip lines related to ZZAAA
                 cycle
               endif
               do iZA=1,mxK
                 if (iZ.eq.iZk(iZA).and.iA.eq.iAk(iZA)) then
                   iZAindex = iZA
                   if (iline.eq.1.or.iline.eq.2) then
                     read (line,*) Cn(iline,0,iZAindex), &
                                   Cn(iline,1,iZAindex), &
                                   Cn(iline,2,iZAindex)
                   elseif (iline.eq.3) then
                     read (line,*) Dn(1,iZAindex), &
                                   Dn(2,iZAindex)
                   elseif (iline.ge.4.and.iline.le.6) then
                     read (line,*) Wn(iline-4,0,iZAindex), &
                                   Wn(iline-4,1,iZAindex), &
                                   Wn(iline-4,2,iZAindex)
                   endif
                 end if
               end do
             end do ! read 6 lines of data per isotope
             if (eof.eqv..TRUE.) exit
           end do ! read until eof
         end if ! eof != true
         dataread=.TRUE.
         CLOSE (L7)
       endif ! dataread == false
!
#ifdef WRITEL6
       write (L6,*) iZk(iK),iAk(iK),epsmx,Mthk(iK)
#endif
       iZ0=iZk(iK); iA0=iAk(iK); Mth=Mthk(iK)           ! Z0,A0,Mth
#ifdef WRITEL8
       write (L8,*) 'SplitFits:',iZ0,iA0,epsmx,Mth
!
       write (L8,*) 'Fit parameter values:'
#endif
#ifdef WRITEL6
       do n=0,Mxg
         if (n.gt.0) write (L6,101) n,'Weight:',(Cn(n,i,iK),i=0,2)
         if (n.gt.0) write (L6,101) n,'Shift:',Dn(n,iK)
         write (L6,101) n,'Width:',(Wn(n,i,iK),i=0,2)
       enddo
#endif

       iZ=iZ0                             ! Charge number of fissioning nucleus;
       iA=iA0                             ! its mass number;
       epsAve=epsmx                       ! its average excitation, since:
! OBS:                                    ! epsAve=epsmx=En+Sn=fixed for now
       epsm=epsmx                         ! Maximum excitation needed
       mE=int(2+epsm/dEf)                 ! Number of energy entries in table
       mode=0; m0=0; i0=int(0.5+epsm/dEf) !>  for plot
       E0=epsmx-msFREYA_SEPn(iK,iZ,iA)    !>  in AfDist
! Consider emission of up to Mth pre-fission neutrons:
       DO m=0,Mth
#ifdef WRITEL8
         write (L8,*) 'SplitFits for',m,' pre-fission neutrons: ======='
#endif
         Sn=msFREYA_SEPn(iK,iZ,iA)        ! Neutron separation energy in A=A0-m
         do iE=1,mE                       ! Tabulate for all excitations eps:
           eps=(iE-0.5)*dEf               ! Nuclear excitation energy
           En=eps-Sn                      ! Equiv neutron kin erg (may be <0),
! OBS:                                    ! since RV formula is in terms of En
! Calculate mode weights Cn & centroids An & widths Wn:
           CnNE(0,m,iE,iK)=1.0            ! n=0: weight before subtraction of n>0
           AnNE(0,m,iE,iK)=0.5*iA         ! n=0: centroid = half of total
           WnNE(0,m,iE,iK)=Wn(0,0,iK)     ! n=0: width
           do n=1,Mxg                     ! n>0:
             CnNE(n,m,iE,iK) &
               = Cn(n,0,iK)/(1.0+exp((En-Cn(n,1,iK))/Cn(n,2,iK))) ! Cn
             CnNE(0,m,iE,iK)=CnNE(0,m,iE,iK)-CnNE(n,m,iE,iK)      ! C0
             AnNE(n,m,iE,iK)=AnNE(0,m,iE,iK)+Dn(n,iK)             ! An
             WnNE(n,m,iE,iK) &
               = Wn(n,0,iK)+En*(Wn(n,1,iK)+En*Wn(n,2,iK))         ! Wn
           enddo
           FnNE(0,m,iE,iK)=CnNE(0,m,iE,iK)! Cumulative wgt (Fn:C0->1)! n=0
           do n=1,Mxg
             FnNE(n,m,iE,iK)=FnNE(n-1,m,iE,iK)+CnNE(n,m,iE,iK)!+(n>0)
           enddo                          ! FnNE(Mxg,m,iE) should now be unity but we put
           FnNE(Mxg,m,iE,iK)=1.0          ! to avoid numerical trouble
#ifdef WRITEL8
           write (L8,*) 'Gaussian Fits for Nth=',m,' E*=',eps,En,':'
           write (L8,102) 'Weight:',(CnNE(n,m,iE,iK),n=0,Mxg)
           write (L8,102) 'Sum(C):',(FnNE(n,m,iE,iK),n=0,Mxg)
           write (L8,102) 'Center:',(AnNE(n,m,iE,iK),n=0,Mxg)
           write (L8,102) 'Widths:',(WnNE(n,m,iE,iK),n=0,Mxg)
           write (L8,*) 'modes:      n=0       n=1       n=2'
  102      format (a7,5f10.4)
#endif
         enddo
         if (m.lt.Mth) then
           iA=iA-1                        ! Mass number of the daughter;
           epsm=epsm-Sn                   ! its maximum excitation.
           epsAve=epsAve-Sn               ! Average maximum excitation;
!           aA=alevel(iK,iZ,iA,epsAve,U)  ! average maximum level density;
           aA=iA/8.0; U=epsAve!!
           Tm=sqrt(U/aA)                  ! average maximum temperature;
           epsAve=epsAve-2*Tm             ! average mean excitation.
         endif 
       ENDDO
       call FLUSH (L8)

  101  format (i3,2x,a7,5f10.4)
       RETURN
       END

!************************************************************************
!               BEGIN MATHEMATICAL ROUTINES:
!************************************************************************
       subroutine SetupRandom (mix,N)
! sets up & checks the random-number generators ran1
! mix  number of times ran1 is called to scramble iseed
! N    sample size for testing ran1

       use freyaout
       use freyaconsts

       implicit none

       integer mix,N

       integer i,j,k,Nran,nextN
       real a,ran1,ran,ranC,s1,s2,s3,s4,eta,xnormal

       parameter (Nran=8)
       dimension ran(Nran),nextN(Nran),ranC(Nran)

#ifdef WRITEL6
       write (L6,*) 'Default seed for random number generator:',iseed
#endif
       do j=1,mix
         a=ran1(iseed)
!         write (L6,*) j,iseed
       enddo
#ifdef WRITEL6
       write (L6,"(' iseed is now',i10,' after',i8,' calls')") iseed,mix
       write (L6,*) 'Current seed for random number generator:',iseed
#endif
       IF (N.gt.0) THEN
! Check the random number generator ran1:
         do i=1,Nran
           nextN(i)=i+1
           ranC(i)=0.0
           ran(i)=ran1(iseed)
         enddo
         nextN(Nran)=1
         j=1
         s1=0.0
         s2=0.0
         do i=1,N
           eta=ran1(iseed)
           s1=s1+eta
           s2=s2+eta**2
           do k=1,Nran
             ranC(k)=ranC(k)+eta*ran(j)
             j=nextN(j)
           enddo
           ran(j)=eta
           j=nextN(j)
         enddo
! Check that <eta>=1/2 and <eta**2>-<eta>**2=1/12 for P(eta) uniform on (0,1):
         s1=s1/N
         s2=s2/N-s1**2
#ifdef WRITEL6
         write (L6,"('Check RAN1:',f10.7,'  ?=1=? ',f10.7)") 2*s1, 12*s2
#endif
! Correlations between a random number and an earlier one (k removed):
         do k=1,Nran
           ranC(k)=ranC(k)/N-s1**2
         enddo
#ifdef WRITEL6
         write (L6,"('Corr:',8f9.5)") (ranC(k),k=1,Nran)
#endif

!------------------------------------------------------------------------
! Check the normal distribution generator xnormal:
#ifdef WRITEL6
         s1=0.0
         s2=0.0
         s3=0.0
         s4=0.0
         do i=1,N
           eta=xnormal()
           s1=s1+eta
           s2=s2+eta**2
           s3=s3+eta**3
           s4=s4+eta**4
         enddo
         s1=s1/N
         s2=s2/N
         s3=s3/N
         s4=s4/N
         write (L6,"(10f8.4)") s1,s2,s3,s4,3*s2**2
         write (L6,*) '   0       1       0       3   =?  3'
#endif
       ENDIF

       RETURN
       END

!************************************************************************
       subroutine msFREYA_usehostrng() &
       bind(C, name="msfreya_usehostrng_")
!      sets the code to use the host random number generator rngdptr
       
       use freyarng

       implicit none

       user_rng=.TRUE.

       RETURN
       END

!************************************************************************
       function rng (iseed)
!      returns a random number generated by
!                 user_rng
!        ran1     false
!        rngdptr  true

       use freyarng
       use iso_c_binding, only: C_DOUBLE

       implicit none

       real rng
       integer iseed

       real (kind=c_double) :: rngdptr
       external rngdptr
       real ran1

       if (user_rng.eqv..TRUE.) then
          rng = real(rngdptr())
       else
          rng = ran1(iseed)
       endif

       RETURN
       END

!************************************************************************
       function ran0 (iseed)   ! This is ran2 from Numerical Recipes p197.
!      it returns a uniform random number in the half-open interval [0,1);
!      set iseed to any negative value to (re)initialize the sequence.
!
!      parameter (m=714025, ia=1366, ic=150889, rm=1./m)
       use freyaout
       use freyaMth

       implicit none

       real ran0
       integer iseed

       integer, save :: iran
       dimension iran(97)
       integer, save :: iy

       integer ia,ic,m
       data m,ia,ic /714025,1366,150889/
       integer iff
       real rm
       data rm,iff /0.0,0/

       integer j

       if (iseed.lt.0.OR.iff.eq.0) then
         iff=1
         rm=1./m
         iseed=MOD(ic-iseed,m)
         do j=1,97              ! Initialize the shuffle table
           iseed=MOD(ia*iseed+ic,m)
           iran(j)=iseed
         enddo
         iseed=MOD(ia*iseed+ic,m)
         iy=iseed
       endif
       j=1+(97*iy)/m            ! Starts here except on initialization
       if (j.gt.97.OR.j.lt.1) then
#ifdef WRITEL6
         write (L6,*) 'ran2: j out of range (1-97):',j
#endif
         ! PAUSE
       endif
       iy=iran(j)
       ran0=iy*rm
! This function ran0 (= ran2 from NR) maps onto [0,1) so it never becomes 1;
! to instead map onto (0,1], do for example 0 -> 1 as follows:
! change: if (iy*rm.eq.0.0) ran2=1.0      ! Change to avoid ran2=0.0
! or use the function ran1 = 1-ran0 (is defined below).
       iseed=MOD(ia*iseed+ic,m)
       iran(j)=iseed
       RETURN
       END

!************************************************************************
       function ran1 (iseed)    !       See Numerical Recipes p197.
!      it returns a uniform random number in the half-open interval (0,1];
!      set iseed to any negative value to (re)initialize the sequence.

       implicit none

       integer iseed
       real ran1

       real ran0

       ran1 = 1.0-ran0(iseed)
       RETURN
       END

!***********************  NORMAL:  *************************************
       function xnormal ()      !	normal-distributed number:
!      x has a normal distribution with zero mean and unit variance.
       use freyaconsts

       implicit none

       real xnormal

       integer m
       real x,y
       data x,y,m/0.0,0.0,2/

       real a,rng,phi,cos,sin

       m=3-m
       if (m.eq.1) then
    1    a=rng(iseed)           !>
         if (a.eq.0.) goto 1    ! >     choose modulus
         a=sqrt(-2.*alog(a))    !>
         phi=twopi*rng(iseed)   !       choose angle
         x=a*cos(phi)           !       x component
         y=a*sin(phi)           !       y component
         xnormal=x
       else
         xnormal=y
       endif
       RETURN
       end

!************************************************************************
