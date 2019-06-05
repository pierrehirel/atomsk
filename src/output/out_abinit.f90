MODULE out_abinit
!
!**********************************************************************************
!*  OUT_ABINIT                                                                    *
!**********************************************************************************
!* This module writes a file in the ABINIT format.                                *
!* The ABINIT documentation is available here:                                    *
!*     https://www.abinit.org/                                                    *
!**********************************************************************************
!* (C) May 2019 - Pierre Hirel                                                    *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 04 June 2019                                     *
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
SUBROUTINE WRITE_ABINIT(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=12):: test
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: isreduced
INTEGER:: i, j
INTEGER:: typecol  !column of AUX containing the atom types
REAL(dp):: a, b, c, alpha, beta, gamma
REAL(dp):: smass, snumber
REAL(dp):: P1, P2, P3
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: atypes
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: G
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
isreduced = .FALSE.
typecol = 0
G(:,:) = 0.d0
i=1
!
IF( ALLOCATED(AUXNAMES) .AND. SIZE(AUXNAMES)>0 .AND. ALLOCATED(AUX) ) THEN
  DO i=1,SIZE(AUXNAMES)
    IF( AUXNAMES(i)=="type" ) THEN
      typecol = i
      EXIT
    ENDIF
  ENDDO
ENDIF
!
100 CONTINUE
msg = 'entering WRITE_ABINIT'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Find the number of atoms of each type
CALL FIND_NSP(P(:,4),atypes)
!
OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN',ERR=250)
!
!Check if coordinates are reduced or not
CALL FIND_IF_REDUCED(P,isreduced)
WRITE(msg,*) 'isreduced:', isreduced
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
IF( .NOT.isreduced ) THEN
  !Calculate the inverse of matrix H
  msg = 'inverting matrix H'
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  CALL INVMAT(H,G)
ENDIF
!
! Write header of ABINIT file
CONTINUE
!Write comments
DO i=1,SIZE(comment)
  WRITE(40,*) TRIM(ADJUSTL(comment(i)))
ENDDO
!
!Write lattice constant
WRITE(40,'(a12)') "acell    3*1"
!
!Write cell parameters
WRITE(40,'(a9,3(f12.8,2X))') "rprim    ", H(1,1), H(2,1), H(3,1)
WRITE(40,'(9X,3(f12.8,2X))') H(1,2), H(2,2), H(3,2)
WRITE(40,'(9X,3(f12.8,2X))') H(1,3), H(2,3), H(3,3)
!
!Write number of atoms
WRITE(temp,*) SIZE(P,1)
WRITE(40,'(a)') "natom    "//TRIM(ADJUSTL(temp))
!
!Write number of types of atoms
WRITE(temp,*) SIZE(atypes,1)
WRITE(40,'(a)') "ntypat   "//TRIM(ADJUSTL(temp))
!
!Write the type of each atom
temp = " "
IF( typecol>0 ) THEN
  !The "type" is defined as auxiliary property: use it
  DO j=1,SIZE(AUX,1)
    WRITE(msg,*) NINT(AUX(j,typecol))
    temp = TRIM(ADJUSTL(temp))//" "//TRIM(ADJUSTL(msg))
  ENDDO
  !
ELSE
  !The "type" is not defined: use the index of this atom species in aentries
  DO i=1,SIZE(P,1)
    DO j=1,SIZE(atypes,1)
      IF( NINT(P(i,4))==NINT(atypes(j,1)) ) THEN
        WRITE(msg,*) j
        temp = TRIM(ADJUSTL(temp))//" "//TRIM(ADJUSTL(msg))
        EXIT
      ENDIF
    ENDDO
  ENDDO
ENDIF
WRITE(40,'(a)') "typat    "//TRIM(ADJUSTL(temp))
!
!Write atomic number for each type
temp = " "
DO j=1,SIZE(atypes,1)
  WRITE(msg,*) NINT(atypes(j,1))
  temp = TRIM(ADJUSTL(temp))//" "//TRIM(ADJUSTL(msg))
ENDDO
WRITE(40,'(a)') "znucl    "//TRIM(ADJUSTL(temp))
!
!
200 CONTINUE
! Write atom coordinates
IF(isreduced) THEN
  WRITE(40,'(a9,3(f12.8,2X))') "xred     ", P(1,1), P(1,2), P(1,3)
ELSE
  WRITE(40,'(a9,3(f12.8,2X))') "xcart    ", P(1,1), P(1,2), P(1,3)
ENDIF
IF(SIZE(P,1)>1) THEN
  DO i=2,SIZE(P,1)
    WRITE(40,'(9X,3(f12.8,2X))') P(i,1), P(i,2), P(i,3)
  END DO
ENDIF
!
GOTO 300
!
250 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
300 CONTINUE
msg = "ABINIT"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
1000 CONTINUE
CLOSE(40)
!
END SUBROUTINE WRITE_ABINIT
!
!
END MODULE out_abinit
