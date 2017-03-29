MODULE deterH
!
!**********************************************************************************
!*  DETERH                                                                        *
!**********************************************************************************
!* This module tries to determine the supercell vectors                           *
!* when they were not found in the input file.                                    *
!* WARNING: this option tries to provide a rectangular bounding box.              *
!* It cannot be expected to be accuratefor any system, and especially not         *
!* for non-rectangular systems (e.g. triclinic) or for molecules.                 *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 29 March 2017                                    *
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
INTEGER:: i, j, n
INTEGER,DIMENSION(3):: Natoms  !number of atoms in each direction
REAL(dp),PARAMETER:: tol=0.3d0   !tolerance on atom positions (allow for small perturbations)
REAL(dp):: maxd        !maximum distance between to atoms along a Cartesian direction
REAL(dp):: maxP, minP  !maximum and minimum coordinates along an axis
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
!The idea is to start from the atom with the smallest coordinate,
!and then along each Cartesian direction, find the periodicity at which
!equivalent atoms are found. Then, the box vector is an integer number times
!this periodicity.
!NOTE: this works only for crystalline materials, and only for parallelepipedic boxes.
!If no suitable period is found, then a surrounding box is still provided, but
!it cannot be expected to be well suited for the system.
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
  !The atom with smallest coordinate along current direction is used as reference
  n = MINLOC(P(:,j),1)
  Pfirst(:) = P(n,:)
  minP = Pfirst(j)
  !
  DO i=2,SIZE(P,1)
    !Go along current direction j and find all atoms that
    !have the same coordinates in the other 2 directions
    IF( DABS(P(i,a2)-Pfirst(a2))<=tol .AND. DABS(P(i,a3)-Pfirst(a3))<=tol ) THEN
      !Increase the number of atoms found in that direction
      Natoms(j) = Natoms(j)+1
      !Let baseunit = distance between first atom and
      !the farest one along the current direction
      IF( P(i,j)-Pfirst(j)>baseunit(j) ) THEN
        baseunit(j) = P(i,j)-Pfirst(j)
      ENDIF
      IF( DABS(P(i,4)-Pfirst(4)) < 1.d-12 ) THEN
        maxP = P(i,j)
      ENDIF
    ENDIF
  ENDDO
  WRITE(msg,*) j, ' Period (A) = ', baseunit(j)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,*) j, ' Number of equivalent atoms found = ', Natoms(j)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !
  IF( Natoms(j) > 1 ) THEN
    !Normalize baseunit to the number of atoms found along direction j
    !This should give the average distance between atoms along direction j
    baseunit(j) = baseunit(j) / MAX( 1.d0 , DBLE(Natoms(j)) )
    !
    !The base vectors of H are the minimum coordinates + ( n*baseunit )
    n = CEILING( (maxP-minP) / baseunit(j) )
    H(j,j) = minP + n*baseunit(j)
    WRITE(msg,*) ' Renormalized period = ', baseunit(j)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  ENDIF
  !
  !Avoid extremely small or negative dimensions
  !Moreover, box vector should be MAXVAL-MINVAL, plus or minus 5 A
  ! if it is not the case then the previous determination was wrong
  maxd = MAXVAL(P(:,j)) - MINVAL(P(:,j))
  IF( H(j,j)<2.d0 .OR. DABS(H(j,j)-maxd)>=5.d0 ) THEN
    WRITE(msg,*) ' Re-setting H = ', H(j,j)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    H(j,j) = maxd + 2.d0
  ENDIF
  !
  ! IF( 
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
