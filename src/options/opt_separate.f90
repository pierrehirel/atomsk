MODULE separate
!
!**********************************************************************************
!*  SEPARATE                                                                      *
!**********************************************************************************
!* This module separates atoms that are too close.                                *
!**********************************************************************************
!* (C) February 2017 - Pierre Hirel                                               *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 21 Feb. 2017                                     *
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
USE atoms
USE comv
USE constants
USE messages
USE files
USE neighbors
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE SEPARATE_XYZ(H,P,S,sep_radius,sep_dist,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: msg
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT   !mask for atom list
INTEGER:: i, iat, j
INTEGER:: Nsep     !number of pairs of atoms that were separated
INTEGER,DIMENSION(:,:),ALLOCATABLE:: NeighList  !list of neighbors
REAL(dp),INTENT(IN):: sep_radius  !atoms closer than this distance will be pulled apart
REAL(dp),INTENT(IN):: sep_dist    !atoms will be pulled apart by this distance
REAL(dp):: vl  !length of a vector
REAL(dp),DIMENSION(3):: Pvector   !unit vector between the two atoms to separateclose
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: PosList  !final positions of neighbors
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S  !positions of atoms, shells
!
species = ''
i = 0
Nsep = 0
!
WRITE(msg,*) 'Entering SEPARATE_XYZ: ', sep_radius, sep_dist
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
CALL ATOMSK_MSG(2140,(/""/),(/sep_radius,sep_dist/))
!
!
!
100 CONTINUE
!Construct neighbor list
CALL ATOMSK_MSG(11,(/""/),(/0.d0/))
CALL NEIGHBOR_LIST(H,P,sep_radius,NeighList)
!
!Find atoms that must be pulled apart
DO i=1,SIZE(P,1)-1  !Loop on all atoms
  !
  !Find positions of neighbours of P(i,:) that are closer than rmd_radius
  CALL NEIGHBOR_POS(H,P(:,:),P(i,1:3),NeighList(i,:),ALLOCATED(NeighList),sep_radius,PosList)
  !Now PosList(:,:) contains the cartesian positions of all neighbors in the rmd_radius,
  !their distance to the atom #i, and their indices.
  !
  !If a neighboring atom is within the rmd_radius, displace the 2 atoms
  IF( ALLOCATED(PosList) .AND. SIZE(PosList,1)>0 ) THEN
    DO j=1,SIZE(PosList,1)
      iat = NINT(PosList(j,5))  !index of current neighbor
      IF( iat>0 .AND. iat<=SIZE(P,1) ) THEN
        IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(iat) ) THEN  !remove only selected atoms
          Pvector(:) = P(i,1:3) - PosList(j,1:3)
          vl = VECLENGTH(Pvector(:))
          IF( vl > 1d-6 ) THEN
            !Normalize vector length
            Pvector(:) = Pvector(:) / VECLENGTH(Pvector(:))
          ELSE
            !Atoms are extremely close
            !Make sure that Pvector has a length 1
            Pvector(:) = DSQRT(3.d0)
          ENDIF
          P(i,1:3) = P(i,1:3) + sep_dist*Pvector(:)
          P(iat,1:3) = P(iat,1:3) - sep_dist*Pvector(:)
          Nsep = Nsep+1
        ENDIF
      ENDIF
    ENDDO
    DEALLOCATE(PosList)
  ENDIF
  !
ENDDO
!
!
!
300 CONTINUE
!
!
CALL ATOMSK_MSG(2141,(/""/),(/DBLE(Nsep)/))
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
!
!
!
1000 CONTINUE
!
!
!
END SUBROUTINE SEPARATE_XYZ
!
END MODULE separate
