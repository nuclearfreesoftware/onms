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
       module freyainterfaces
!
!      interfaces between FORTRAN and C for bindings
!      ---------------------------------------------
!
       interface
         subroutine msfreya_event(iK,Einc,eps0,PP0,iZ1,iA1, &
          PP1,iZ2,iA2,PP2,n,p,id) bind (C, name="msfreya_event_")
           use, intrinsic :: iso_c_binding
           use freyaparameters
           implicit none
           integer (kind=c_int), value :: iK
           real (kind=c_float), value :: Einc,eps0
           ! Four-momentum of the fissioning nucleus
           real (kind=c_float), dimension(0:4) :: PP0
           ! Four-momentum of the 1st evaporating nucleus
           real (kind=c_float), dimension(0:4) :: PP1
           ! Four-momentum of the 2nd evaporating nucleus
           real (kind=c_float), dimension(0:4) :: PP2
           integer (kind=c_int) :: iZ1,iA1,iZ2,iA2
           integer (kind=c_int) :: n
           ! Ejectile momenta of all n=n0+n1+n2 neutrons
           real (kind=c_float), dimension(4,3*mMax) :: p
           ! Ejectile type 
           integer (kind=c_int), dimension(3*mMax) :: id
         end subroutine msfreya_event
       end interface
!
       interface
         real (kind=c_float) function msFREYA_SEPn_CVT (iK,iZ,iA)  &
         bind (C, name="msfreya_sepn_cvt_")
           use, intrinsic :: iso_c_binding
           implicit none
           integer (kind=c_int), value :: iK,iZ,iA
         end function msFREYA_SEPn_CVT
       end interface
!
       interface
         subroutine msfreya_setup () &
         bind (C, name="msfreya_setup_")
           use, intrinsic :: iso_c_binding
           implicit none
         end subroutine msfreya_setup
       end interface
!
       interface
         subroutine msfreya_getniso(nisosf,nisoif) &
         bind (C, name="msfreya_getniso_")
           use, intrinsic :: iso_c_binding
           implicit none
           integer (kind=c_int) :: nisosf,nisoif
         end subroutine msfreya_getniso
       end interface
!
       interface
         subroutine msfreya_getzas(zas,fistypes) &
         bind (C, name="msfreya_getzas_")
           use, intrinsic :: iso_c_binding
           implicit none
           integer (kind=c_int), dimension(*) :: zas,fistypes
         end subroutine msfreya_getzas
       end interface
!
       interface
         real (kind=c_float) function msFREYA_gsMassN (iZ,iA)  &
         bind (C, name="msfreya_gsmassn_")
           use, intrinsic :: iso_c_binding
           implicit none
           integer (kind=c_int), value :: iZ,iA
         end function msFREYA_gsMassN
       end interface
!
       interface
         subroutine msFREYA_usehostrng() &
         bind(C, name="msfreya_usehostrng_")
           use, intrinsic :: iso_c_binding
           implicit none
         end subroutine msFREYA_usehostrng
       end interface
!
       interface
         subroutine msFREYA_getlasterror(message,counter) &
         bind(C, name="msfreya_getlasterror_")
           use, intrinsic :: iso_c_binding
           implicit none
           character (kind=c_char,len=1), dimension(256) :: message
           integer (kind=c_int) :: counter
         end subroutine msFREYA_getlasterror
       end interface
!
       interface
         subroutine msFREYA_geterrors(messages,length) &
         bind(C, name="msfreya_geterrors_")
           use, intrinsic :: iso_c_binding
           implicit none
           character (kind=c_char,len=1), dimension(length) :: messages
           integer (kind=c_int) :: length
         end subroutine msFREYA_geterrors
       end interface
!
       interface
         subroutine msFREYA_reseterrorflag() &
         bind(C, name="msfreya_reseterrorflag_")
           use, intrinsic :: iso_c_binding
           implicit none
         end subroutine msFREYA_reseterrorflag
       end interface
!
       interface
         logical (kind=c_bool) function msFREYA_errorflagset ()  &
         bind (C, name="msfreya_errorflagset_")
           use, intrinsic :: iso_c_binding
           implicit none
         end function msFREYA_errorflagset
       end interface
!
       end module freyainterfaces
