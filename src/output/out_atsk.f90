MODULE out_atsk
!
!
!**********************************************************************************
!*  OUT_ATSK                                                                      *
!**********************************************************************************
!* This module writes a file in the binary format of atomsk.                      *
!* This file format is specific to atomsk and should *NOT* be used to             *
!* store or archive data!!                                                        *
!**********************************************************************************
!* (C) Oct. 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 30 July 2015                                     *
!**********************************************************************************
!* This program is free software: you can redistribute it and/or modify           *
!* it under the terms of the GNU General Public License as published by           *
!* the Free Software Foundation, either version 3 of the License, or              *
!* (at your option) any later version.                                            *
!*                                                                                *
!* This program is distributed in the hope that it will be useful,                *
!* but WITHOUT ANY WARRANTY; without even the implied warranty of                 *
!* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  *
!* GNU General Public License for more details.                                   *
!*                                                                                *
!* You should have received a copy of the GNU General Public License              *
!* along with this program.  If not, see <http://www.gnu.org/licenses/>.          *
!**********************************************************************************
!
USE comv
USE constants
USE functions
USE messages
USE files
USE subroutines
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE WRITE_ATSK(H,P,comment,AUXNAMES,AUX,outputfile,S_in)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=32):: atsktxt
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
INTEGER:: NS, Naux, Nauxnames, Ncomment
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN),OPTIONAL:: S_in
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
IF(ALLOCATED(S)) DEALLOCATE(S)
IF( PRESENT(S_in) .AND. ALLOCATED(S_in) .AND. SIZE(S_in,1)>0 ) THEN
  ALLOCATE( S(SIZE(S_in,1),SIZE(S_in,2)) )
  S(:,:) = S_in(:,:)
ENDIF
NS=0
Naux=0
Nauxnames=0
!
msg = 'entering WRITE_ATSK'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
!
100 CONTINUE
OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN',FORM="UNFORMATTED",ERR=500)
!
!Write a header to indicate what the file is
WRITE(atsktxt,'(a32)') '0.8 Atomsk binary file'
atsktxt = ADJUSTL(atsktxt)
WRITE(40) atsktxt
!
!Write the length of each array
IF( ALLOCATED(S) .AND. SIZE(S,1)>0 ) THEN
  NS = SIZE(S,1)
ELSE
  NS = 0
ENDIF
IF( ALLOCATED(AUX) .AND. SIZE(AUX,1)>0 ) THEN
  Naux = SIZE(AUX,1)
  Nauxnames = SIZE(AUX,2)
ELSE
  Naux = 0
  Nauxnames = 0
ENDIF
IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
  Ncomment = SIZE(comment)
ELSE
  Ncomment = 0
ENDIF
WRITE(40) SIZE(P,1), NS, Naux, Nauxnames, Ncomment
!
!Write arrays
WRITE(40) H
WRITE(40) P
IF(NS>0) THEN
  WRITE(40) S
ENDIF
IF(Naux>0) THEN
  WRITE(40) AUXNAMES
  WRITE(40) AUX
ENDIF
IF(Ncomment>0) THEN
  WRITE(40) comment
ENDIF
!
!
!
500 CONTINUE
CLOSE(40)
msg = "ATSK"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
1000 CONTINUE
!
END SUBROUTINE WRITE_ATSK
!
END MODULE out_atsk
