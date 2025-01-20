MODULE out_bop
!
!
!**********************************************************************************
!*  IN_BOP                                                                        *
!**********************************************************************************
!* This module writes files used by the Bond-Order Potential code (BOP).          *
!* The BOP format is not really documented.                                       *
!**********************************************************************************
!* (C) July 2011 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 14 Jan. 2025                                     *
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
SUBROUTINE WRITE_BOP(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=2):: species
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: isreduced
INTEGER:: fixx, fixy, fixz
INTEGER:: i, NPd, NPinert
REAL(dp):: len1, len2, len3
REAL(dp):: P1, P2, P3
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: G   !Inverse of H
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
fixx=0
fixy=0
fixz=0
NPd = 0
NPinert = 0
len1 = 1.d0
len2 = 1.d0
len3 = 1.d0
!
!
100 CONTINUE
IF(ofu.NE.6) THEN
  OPEN(UNIT=ofu,FILE=outputfile,FORM='FORMATTED',STATUS='UNKNOWN')
ENDIF
!
!Check if coordinates are already reduced or not
CALL FIND_IF_REDUCED(H,P,isreduced)
WRITE(msg,*) 'isreduced:', isreduced
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Inverse matrix H
CALL INVMAT(H,G)
!
!Check if some atoms are fixed
IF(ALLOCATED(AUX)) THEN
  DO i=1,SIZE(AUXNAMES)
    IF(AUXNAMES(i)=="freeze_x") fixx=i
    IF(AUXNAMES(i)=="freeze_y") fixy=i
    IF(AUXNAMES(i)=="freeze_z") fixz=i
  ENDDO
ENDIF
!
!Count the number of fixed atoms
IF( fixx.NE.0 .AND. fixy.NE.0 .AND. fixz.NE.0 ) THEN
  DO i=1,SIZE(P(:,1))
    IF( AUX(i,fixx)>0.1d0 .OR. AUX(i,fixy)>0.1d0 .OR. AUX(i,fixz)>0.1d0 ) THEN
      NPinert = NPinert+1
    ENDIF
  ENDDO
ELSE
  NPinert = 0
ENDIF
!
!Set the number of mobile atoms
NPd = SIZE(P(:,1))-NPinert
!
!
!
200 CONTINUE
201 FORMAT(3(f16.8,2X),a2,2X,f6.3)
WRITE(ofu,'(a1)') "A"
WRITE(msg,201) H(1,1), H(1,2), H(1,3)
WRITE(ofu,'(a)') TRIM(ADJUSTL(msg))
WRITE(msg,201) H(2,1), H(2,2), H(2,3)
WRITE(ofu,'(a)') TRIM(ADJUSTL(msg))
WRITE(msg,201) H(3,1), H(3,2), H(3,3)
WRITE(ofu,'(a)') TRIM(ADJUSTL(msg))
!
WRITE(ofu,'(a3)') "LEN"
WRITE(msg,201) len1, len2, len3
WRITE(ofu,'(a)') TRIM(ADJUSTL(msg))
!
WRITE(ofu,'(a6)') "LATPAR"
WRITE(ofu,'(a5)') "1.000"
!
!Write coordinates of non-fixed atoms
WRITE(ofu,'(a2)') "ND"
WRITE(msg,*) NPd
WRITE(ofu,'(a)') TRIM(ADJUSTL(msg))
!
IF(NPd>0) THEN
  WRITE(ofu,'(a1)') "D"
  DO i=1,SIZE(P(:,1))
    !Note: internally if AUX(fix)==1 then atom is fixed
    IF( fixx==0 .OR. fixy==0 .OR. fixz==0 .OR.                               &
      & .NOT.(AUX(i,fixx)>0.5d0 .OR. AUX(i,fixy)>0.5d0 .OR. AUX(i,fixz)>0.5d0) ) THEN
      CALL ATOMSPECIES(P(i,4),species)
      IF(isreduced) THEN
        WRITE(msg,201) P(i,1), P(i,2), P(i,3), species, zero
      ELSE
        P1 = P(i,1)
        P2 = P(i,2)
        P3 = P(i,3)
        WRITE(msg,201) P1*G(1,1) + P2*G(2,1) + P3*G(3,1),     &
                     & P1*G(1,2) + P2*G(2,2) + P3*G(3,2),     &
                     & P1*G(1,3) + P2*G(2,3) + P3*G(3,3),     &
                     & species, zero
      ENDIF
      WRITE(ofu,'(a)') TRIM(ADJUSTL(msg))
    ENDIF
  ENDDO
ENDIF
!
!Write coordinates of fixed atoms
WRITE(ofu,'(a6)') "NINERT"
WRITE(msg,*) NPinert
WRITE(ofu,'(a)') TRIM(ADJUSTL(msg))
!
IF(NPinert>0) THEN
  WRITE(ofu,'(a6)') "DINERT"
  DO i=1,SIZE(P(:,1))
    IF( AUX(i,fixx)>0.5d0 .OR. AUX(i,fixy)>0.5d0 .OR. AUX(i,fixz)>0.5d0 ) THEN
      CALL ATOMSPECIES(P(i,4),species)
      IF(isreduced) THEN
        WRITE(msg,201) P(i,1), P(i,2), P(i,3), species, zero
      ELSE
        P1 = P(i,1)
        P2 = P(i,2)
        P3 = P(i,3)
        WRITE(msg,201) P1*G(1,1) + P2*G(2,1) + P3*G(3,1),     &
                     & P1*G(1,2) + P2*G(2,2) + P3*G(3,2),     &
                     & P1*G(1,3) + P2*G(2,3) + P3*G(3,3),     &
                     & species, zero
      ENDIF
      WRITE(ofu,'(a)') TRIM(ADJUSTL(msg))
    ENDIF
  ENDDO
ENDIF
GOTO 1000
!
!
!
400 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
!
!
!
1000 CONTINUE
IF(ofu.NE.6) THEN
  CLOSE(ofu)
ENDIF
!
END SUBROUTINE WRITE_BOP
!
END MODULE out_bop
