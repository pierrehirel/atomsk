MODULE deterH
!
!**********************************************************************************
!*  DETERH                                                                        *
!**********************************************************************************
!* This module tries to determine the supercell vectors                           *
!* when they were not found in the input file.                                    *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 23 Oct. 2013                                     *
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
! WARNING: for now, this works only for rectangular or parallelepipedic boxes
!
!
CONTAINS
!
SUBROUTINE DETERMINE_H(H,P)
!
USE comv
USE constants
USE messages
USE subroutines
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
INTEGER:: a2, a3
INTEGER:: i, j
INTEGER,DIMENSION(3):: Natoms  !number of atoms in each direction
REAL(dp),PARAMETER:: tol=0.1d0
REAL(dp),DIMENSION(3):: baseunit !unit length in each direction
REAL(dp),DIMENSION(4):: Pfirst
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P
!
!
!Initialize variables
Natoms(:) = 0
baseunit(:) = 0.d0
H(:,:) = 0.d0
!
!
100 CONTINUE
msg = 'Entering DETERMINE_H...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!
200 CONTINUE
!The atom with smallest X coordinate is used as a starting point
Pfirst(:) = P(1,:)
!
!The we go along each axis and find atoms identical to Pfirst
!(this works for parallelepipedic boxes only)
DO j=1,3  !loop on x, y, z
  IF(j==1) THEN
    a2 = 2
    a3 = 3
  ELSEIF(j==2) THEN
    a2 = 3
    a3 = 1
  ELSE
    a2 = 1
    a3 = 2
  ENDIF
  !
  DO i=2,SIZE(P,1)
    !Go along current direction j and find all atoms that
    !have the same coordinates in the other 2 directions
    IF( P(i,a2)-Pfirst(a2)<=tol .AND. P(i,a3)-Pfirst(a3)<=tol              &
      & .AND. P(i,4)==Pfirst(4)                               ) THEN
      !Increase the number of atoms found in that direction
      Natoms(j) = Natoms(j)+1
      !Let baseunit = distance between first atom and
      !the farest one along the current direction
      IF( P(i,j)-Pfirst(j)>baseunit(j) ) THEN
        baseunit(j) = P(i,j)-Pfirst(j)
      ENDIF
    ENDIF
  ENDDO
  !Normalize baseunit to the number of atoms found along direction j
  !This should give the average distance between atoms along direction j
  baseunit(j) = baseunit(j) / MAX( 1.d0 , DBLE(Natoms(j)+1) )
ENDDO
!
!
300 CONTINUE
!The base vectors of H are the maximum coordinates + baseunit
DO j=1,3
  H(j,j) = MAXVAL(P(:,j)) - MINVAL(P(:,j)) + baseunit(j)
  !Special case for mono-atomic slabs
  IF( H(j,j)<1.0 ) THEN
    H(j,j) = H(j,j)*2.d0
  ENDIF
  !Avoid zero or negative dimensions
  IF( H(j,j)<=0.d0 ) THEN
    H(j,j) = 2.d0
  ENDIF
ENDDO
!
CALL ATOMSK_MSG(2095,(/''/),(/0.d0/))
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE DETERMINE_H
!
END MODULE deterH
