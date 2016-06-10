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
!************************************************************************
!               ERROR ROUTINES:
!************************************************************************
      SUBROUTINE initerrors ()
!
!      severity=1 will cause the code to stop
!
      use freyaerrors

      implicit none

!
! allocate memory for the errors
!
      allocate(errors(maxErrorCodes))
      allocate(severities(maxErrorCodes))
      severities=0
      severities(1)=1
      severities(4)=1
      severities(5)=1
      severities(6)=1
      severities(7)=1
      severities(8)=1
      severities(9)=1
      severities(10)=1
      severities(11)=1
      severities(12)=1
      severities(13)=1

      allocate(errorctr(maxErrorCodes))
      errorctr=0
      nerrors=0

      allocate(caseerrors(maxCaseErrorCodes))
      allocate(caseseverities(maxCaseErrorCodes))
      caseseverities=0
      caseseverities(3)=1
      caseseverities(4)=1
      caseseverities(5)=1
      caseseverities(6)=1

      allocate(caseerrorctr(maxCaseErrorCodes))
      caseerrorctr=0
      ncaseerrors=0

      lasterrorcode=0
      casespec=.FALSE.

      freyaerrorsinitialized=.TRUE.

      RETURN
      END

      SUBROUTINE seterror (errorcode,message)
!     - if this error code is called for first time, add that error to 
!       the list. Otherwise, replace error message.
!     - increment error counter
!
      use freyaerrors

      implicit none

      integer errorcode
      character*256 message

      if ((errorcode.ge.1).and.(errorcode.le.maxErrorCodes)) then
        errors(errorcode)=message
        errorctr(errorcode)=errorctr(errorcode)+1
        nerrors=nerrors+1
        lasterrorcode=errorcode
        casespec=.FALSE.
        errorflag=.TRUE.
      else
        write(error, 20) 
20      format("FREYA: Error code out of bound")
        errorcode=2
        errors(errorcode)=error
        errorctr(errorcode)=errorctr(errorcode)+1
        nerrors=nerrors+1
        lasterrorcode=errorcode
        casespec=.FALSE.
        errorflag=.TRUE.
      endif

      RETURN
      END

      SUBROUTINE seterrorcase (errorcode,message)
!     - if this error code is called for first time, add that error to 
!       the list. Otherwise, replace error message.
!     - increment counter
!
      use freyaerrors

      implicit none

      integer errorcode
      character*256 message

      if ((errorcode.ge.1).and.(errorcode.le.maxCaseErrorCodes)) then
        caseerrors(errorcode)=message
        caseerrorctr(errorcode)=caseerrorctr(errorcode)+1
        ncaseerrors=ncaseerrors+1
        lasterrorcode=errorcode
        casespec=.TRUE.
        errorflag=.TRUE.
      else
        write(error, 20) 
20      format("FREYA: Case-specific error code out of bound")
        call seterror (3,error)
      endif

      RETURN
      END

      FUNCTION exitonerror()
!     if one error has severity of 1, exits
       
      use freyaerrors

      implicit none

      logical exitonerror

      integer i

      exitonerror=.FALSE.
      do i=1,maxErrorCodes
        if ((errorctr(i).gt.0).and.(severities(i).eq.1)) then
          exitonerror=.TRUE.
          RETURN
        endif
      enddo

      do i=1,maxCaseErrorCodes
        if ((caseerrorctr(i).gt.0).and.(caseseverities(i).eq.1)) then
          exitonerror=.TRUE.
          RETURN
        endif
      enddo

      RETURN
      END

      SUBROUTINE msFREYA_getlasterror(message,counter) &
       bind (C, name="msfreya_getlasterror_")
!     returns last error and the counter of that error
      
      use iso_c_binding, only: C_CHAR, C_INT, c_null_char
      use freyaerrors

      implicit none

      character (kind=c_char,len=1), dimension(256) :: message
      integer (kind=c_int) :: counter

      if (.not.freyaerrorsinitialized) then
        counter=0
        RETURN
      endif

      if (lasterrorcode.ne.0) then
        if (casespec.eqv..FALSE.) then
          message(1)=errors(lasterrorcode)
          counter=errorctr(lasterrorcode)
        else
          message(1)=caseerrors(lasterrorcode)
          counter=caseerrorctr(lasterrorcode)
        endif
      else
        counter=0
      endif
      RETURN
      END

      SUBROUTINE msFREYA_geterrors(messages,length) &
       bind (C, name="msfreya_geterrors_")
!     returns all errors and their number
      
      use iso_c_binding, only: C_CHAR, C_INT, C_NEW_LINE, C_NULL_CHAR
      use freyaerrors

      implicit none

      character (kind=c_char,len=1), dimension(length) :: messages
      integer (kind=c_int) :: length
      character*10 tmp
      integer counter,i
      integer messageindex

      if (.not.freyaerrorsinitialized) then
        length=1
        RETURN
      endif

      messageindex=1
      do counter=1,maxErrorCodes
        if (errorctr(counter).gt.0) then
          write(tmp,"(i7)") errorctr(counter)
          do i=1,len_trim(tmp)
            messages(messageindex)=tmp(i:i)
            messageindex=messageindex+1
            if (messageindex.gt.length) RETURN
          enddo
          messages(messageindex)=':'
          messageindex=messageindex+1
          if (messageindex.gt.length) RETURN
          do i=1,len_trim(errors(counter))
            messages(messageindex)=errors(counter)(i:i)
            messageindex=messageindex+1
            if (messageindex.gt.length) RETURN
          enddo
          messages(messageindex)=c_new_line
          messageindex=messageindex+1
          if (messageindex.gt.length) RETURN
        endif
      enddo

      do counter=1,maxCaseErrorCodes
        if (caseerrorctr(counter).gt.0) then
          write(tmp,"(i7)") caseerrorctr(counter)
          do i=1,len_trim(tmp)
            messages(messageindex)=tmp(i:i)
            messageindex=messageindex+1
            if (messageindex.gt.length) RETURN
          enddo
          messages(messageindex)=':'
          messageindex=messageindex+1
          if (messageindex.gt.length) RETURN
          do i=1,len_trim(caseerrors(counter))
            messages(messageindex)=caseerrors(counter)(i:i)
            messageindex=messageindex+1
            if (messageindex.gt.length) RETURN
          enddo
          messages(messageindex)=c_new_line
          messageindex=messageindex+1
          if (messageindex.gt.length) RETURN
        endif
      enddo
      messages(messageindex)=c_null_char
      length=messageindex

      RETURN
      END

      SUBROUTINE msFREYA_reseterrorflag() &
       bind (C, name="msfreya_reseterrorflag_")
!
!     resets the error flag
!
      use freyaerrors

      implicit none

      errorflag=.FALSE.
      RETURN
      END

      FUNCTION msFREYA_errorflagset() &
       bind (C, name="msfreya_errorflagset_")
!     returns the value of the error flag

      use iso_c_binding, only: C_BOOL
      use freyaerrors

      implicit none

      logical (kind=c_bool) :: msFREYA_errorflagset
      
      msFREYA_errorflagset=errorflag

      RETURN
      END

