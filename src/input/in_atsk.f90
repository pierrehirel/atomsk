MODULE in_atsk
!
!
!**********************************************************************************
!*  IN_ATSK                                                                       *
!**********************************************************************************
!* This module reads a file in the binary format of atomsk.                       *
!**********************************************************************************
!* (C) Oct. 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 22 July 2014                                     *
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
SUBROUTINE READ_ATSK(inputfile,H,P,S,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=32):: atsktxt
CHARACTER(LEN=4096):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: comment
INTEGER:: atskversion  !version of atomsk this binary file was written by
INTEGER:: Ncomment, NP, NS, Naux, Nauxnames
REAL(dp),DIMENSION(3,3),INTENT(OUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: P, S
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: AUX !auxiliary properties
!
!
!Initialize variables
atskversion = 0
Ncomment=0
NP=0
NS=0
Naux=0
Nauxnames=0
!
msg = 'entering READ_ATSK'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
!
100 CONTINUE
OPEN(UNIT=30,FILE=inputfile,STATUS='UNKNOWN',FORM="UNFORMATTED",ERR=800)
!
READ(30,ERR=800) atsktxt
READ(atsktxt,*,ERR=120) atskversion
!
120 CONTINUE
!Read the length of each array
READ(30,ERR=800) NP, NS, Naux, Nauxnames, Ncomment
!
!Allocate arrays and read them
READ(30,ERR=800) H
!
ALLOCATE(P(NP,4))
READ(30,ERR=800) P
!
IF(NS>0) THEN
  ALLOCATE(S(NS,4))
  READ(30,ERR=800) S
ENDIF
!
IF(Naux>0 .AND. Nauxnames>0) THEN
  ALLOCATE( AUXNAMES(Nauxnames) )
  ALLOCATE( AUX(NP,Nauxnames) )
  READ(30,ERR=800) AUXNAMES
  READ(30,ERR=800) AUX
ENDIF
!
IF(Ncomment>0) THEN
  ALLOCATE(comment(Ncomment))
  READ(30,ERR=800) comment
ENDIF
!
GOTO 1000
!
!
!
800 CONTINUE
nerr=nerr+1
!
!
!
1000 CONTINUE
CLOSE(30)
!
!
END SUBROUTINE READ_ATSK
!
END MODULE in_atsk
