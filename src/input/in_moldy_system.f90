MODULE in_moldy
!
!**********************************************************************************
!*  IN_MOLDY                                                                      *
!**********************************************************************************
!* This module reads configuration files in the MOLDY format                      *
!* (by default named "system.in").                                                *
!* The format is described for instance at:                                       *
!*     https://www.wiki.ed.ac.uk/display/ComputerSim/MOLDY+Information            *
!**********************************************************************************
!* (C) Oct. 2011 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
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
!
USE comv
USE constants
USE messages
USE files
USE subroutines
!
!
CONTAINS
!
SUBROUTINE READ_MOLDY(inputfile,H,P,comment,AUXNAMES,AUX)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
INTEGER:: ex, ey, ez  !number of replicas along X, Y, Z
INTEGER:: i, j, o, n, m, NP
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, Q  !atom positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
ex = 0
ey = 0
ez = 0
i = 0
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
!
!
msg = 'entering READ_MOLDY'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!
100 CONTINUE
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!First line is always the number of atoms
READ(30,*,ERR=820,END=820) NP
!
IF(NP==0) GOTO 800
ALLOCATE(P(NP,4))
P(:,:) = 0.d0
!
!2nd line tells how many replicas of the system there are along X, Y and Z
READ(30,*,ERR=800,END=800) ex, ey, ez
IF(ex<=0) ex = 1
IF(ey<=0) ey = 1
IF(ez<=0) ez = 1
!
!Three next lines are the supercell vectors
DO i=1,3
  READ(30,*,ERR=810,END=810) (H(i,j),j=1,3)
ENDDO
!
!
!
200 CONTINUE
!Read coordinates
DO i=1,NP
  READ(30,*,ERR=800,END=800) (P(i,j), j=1,4)
ENDDO
!
!Coordinates are fractional: convert to cartesian
CALL FRAC2CART(P,H)
!
!
!
300 CONTINUE
!If system must be duplicated, do it
IF(ex.NE.1 .OR. ey.NE.1 .OR. ez.NE.1) THEN
  NP = SIZE(P(:,1))*ex*ey*ez
  ALLOCATE(Q(NP,4))
  NP=0
  DO o=1,ez
    DO n=1,ey
      DO m=1,ex
        DO i=1,SIZE(P(:,1))
          NP = NP+1
          Q(NP,1) = P(i,1)+REAL(m-1)*H(1,1)+REAL(n-1)*H(2,1)+REAL(o-1)*H(3,1)
          Q(NP,2) = P(i,2)+REAL(m-1)*H(1,2)+REAL(n-1)*H(2,2)+REAL(o-1)*H(3,2)
          Q(NP,3) = P(i,3)+REAL(m-1)*H(1,3)+REAL(n-1)*H(2,3)+REAL(o-1)*H(3,3)
          Q(NP,4) = P(i,4)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !Replace old P by new Q
  DEALLOCATE(P)
  ALLOCATE(P(SIZE(Q,1),4))
  P(:,:) = Q(:,:)
  DEALLOCATE(Q)
  !
  !Resize the cell dimensions
  H(1,:) = ex*H(1,:)
  H(2,:) = ey*H(2,:)
  H(3,:) = ez*H(3,:)
ENDIF
!
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
GOTO 1000
!
810 CONTINUE
CALL ATOMSK_MSG(1809,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
820 CONTINUE
CALL ATOMSK_MSG(1810,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
830 CONTINUE
CALL ATOMSK_MSG(1804,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
1000 CONTINUE
CLOSE(30)
!
END SUBROUTINE READ_MOLDY
!
!
END MODULE in_moldy
