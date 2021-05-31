MODULE out_crystal
!
!**********************************************************************************
!*  OUT_CRYSTAL                                                                   *
!**********************************************************************************
!* This module writes a file in the CRYSTAL format.                               *
!* The CRYSTAL format is officially described here:                               *
!*     http://www.crystal.unito.it/                                               *
!**********************************************************************************
!* (C) May 2019 - Pierre Hirel                                                    *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 31 May 2021                                      *
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
SUBROUTINE WRITE_CRYSTAL(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=12):: test
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: isreduced
INTEGER:: i, j
REAL(dp):: a, b, c, alpha, beta, gamma
REAL(dp):: smass, snumber
REAL(dp):: P1, P2, P3
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: G
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
isreduced = .FALSE.
G(:,:) = 0.d0
i=1
!
100 CONTINUE
msg = 'entering WRITE_CRYSTAL'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
IF(ofu.NE.6) THEN
  OPEN(UNIT=ofu,FILE=outputfile,STATUS='UNKNOWN',ERR=250)
ENDIF
!
!Check if coordinates are already reduced or not
CALL FIND_IF_REDUCED(H,P,isreduced)
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
! Write header of CRYSTAL file
CONTINUE
!First line is a comment
WRITE(ofu,*) comment(1)
!Second line can be CYSTAL, SLAB, POLYMER, MOLECULE, etc.
!Here "CRYSTAL" is always assumed
WRITE(ofu,'(a7)') "CRYSTAL"
!Then
WRITE(ofu,'(a5)') "0 0 0"
!Space group number: here P1 is always assumed
WRITE(ofu,'(a1)') "1"
!Convert cell into conventional notation
CALL MATCONV(H,a,b,c,alpha,beta,gamma)
!Convert angles into degrees
alpha = RAD2DEG(alpha)
beta  = RAD2DEG(beta)
gamma = RAD2DEG(gamma)
WRITE(ofu,'(6(f12.6,2X))') a, b, c, alpha, beta, gamma
WRITE(temp,*) SIZE(P,1)
WRITE(ofu,'(a)') TRIM(ADJUSTL(temp))
!
!
200 CONTINUE
! Write atom coordinates
DO i=1,SIZE(P,1)
  !
  IF(isreduced) THEN
    WRITE(temp,'(i3,2X,3(f12.8,2X))') NINT(P(i,4)), P(i,1), P(i,2), P(i,3)
  ELSE
    P1 = P(i,1)
    P2 = P(i,2)
    P3 = P(i,3)
    WRITE(temp,'(i3,2X,3(f12.8,2X))')  NINT(P(i,4)), &
                              &  P1*G(1,1) + P2*G(2,1) + P3*G(3,1),     &
                              &  P1*G(1,2) + P2*G(2,2) + P3*G(3,2),     &
                              &  P1*G(1,3) + P2*G(2,3) + P3*G(3,3)
  ENDIF
  !
  !Write the line to the file
  WRITE(ofu,'(a)') TRIM(ADJUSTL(temp))
END DO
! Write keyword "END" to end section
WRITE(ofu,'(a3)') "END"
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
msg = "CRYSTAL"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
1000 CONTINUE
IF(ofu.NE.6) THEN
  CLOSE(ofu)
ENDIF
!
END SUBROUTINE WRITE_CRYSTAL
!
!
END MODULE out_crystal
