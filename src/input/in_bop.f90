MODULE in_bop
!
!
!**********************************************************************************
!*  IN_BOP                                                                        *
!**********************************************************************************
!* This module reads files used by the Bond-Order Potential code (BOP).           *
!* The BOP format is not really documented.                                       *
!**********************************************************************************
!* (C) July 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 19 March 2014                                    *
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
USE messages
USE files
USE subroutines
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE READ_BOP(inputfile,H,P,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: test
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
INTEGER:: i, NP, NPd, NPinert
REAL(dp):: a0
REAL(dp):: len1, len2, len3
REAL(dp), DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P
REAL(dp), DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
i = 0
NP = 0
NPd = 0
NPinert = 0
a0 = 1.d0
len1 = 1.d0
len2 = 1.d0
len3 = 1.d0
H(:,:) = 0.d0
ALLOCATE(AUXNAMES(3))
AUXNAMES(1) = "fixx"
AUXNAMES(2) = "fixy"
AUXNAMES(3) = "fixz"
!
!
100 CONTINUE
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
REWIND(30)
DO
  READ(30,'(a128)',ERR=110,END=110) test
  test = TRIM(ADJUSTL(test))
  IF(test(1:1)=='A') THEN
    READ(30,*,ERR=410,END=410) H(1,1), H(1,2), H(1,3)
    READ(30,*,ERR=410,END=410) H(2,1), H(2,2), H(2,3)
    READ(30,*,ERR=410,END=410) H(3,1), H(3,2), H(3,3)
  ELSEIF(test(1:3)=='LEN') THEN
    READ(30,*,ERR=400,END=400) len1, len2, len3
  ELSEIF(test(1:2)=='ND') THEN
    READ(30,*,ERR=420,END=420) NPd
  ELSEIF(test(1:6)=='NINERT') THEN
    READ(30,*,ERR=420,END=420) NPinert
  ENDIF
ENDDO
!
110 CONTINUE
H(1,:) = H(1,:)*len1
H(2,:) = H(2,:)*len2
H(3,:) = H(3,:)*len3
!
!
!
200 CONTINUE
NP = NPd+NPinert
!If number of particles is zero we are in trouble
IF(NP==0) GOTO 400
ALLOCATE(P(NP,4))
ALLOCATE(AUX(NP,3))
AUX(:,:) = 0.d0
!
!Read atomic positions and save them to P
REWIND(30)
DO
  READ(30,'(a128)',END=250,ERR=400) test
  test = TRIM(ADJUSTL(test))
  IF(test(1:3)=='D  ') THEN
    DO i=1,NPd
      READ(30,*,ERR=400,END=400) P(i,1), P(i,2), P(i,3), species
      CALL ATOMNUMBER(species,P(i,4))
      P(i,1) = len1*P(i,1)
      P(i,2) = len2*P(i,2)
      P(i,3) = len3*P(i,3)
    ENDDO
  ELSEIF(test(1:6)=='DINERT') THEN
    DO i=NPd+1,NP
      READ(30,*,ERR=400,END=400) P(i,1), P(i,2), P(i,3), species
      CALL ATOMNUMBER(species,P(i,4))
      P(i,1) = len1*P(i,1)
      P(i,2) = len2*P(i,2)
      P(i,3) = len3*P(i,3)
      AUX(i,1:3) = 1.d0
    ENDDO
  ENDIF
ENDDO
!
250 CONTINUE
GOTO 500
!
!
!
400 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
GOTO 500
!
410 CONTINUE
CALL ATOMSK_MSG(1809,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 500
!
420 CONTINUE
CALL ATOMSK_MSG(1810,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 500
!
!
!
500 CONTINUE
CLOSE(30)
!
END SUBROUTINE READ_BOP
!
END MODULE in_bop
