MODULE remdoubles
!
!**********************************************************************************
!*  REMDOUBLES                                                                    *
!**********************************************************************************
!* This module searches for atoms that have the same position or are within       *
!* a given distance, and removes all of them but the first one.                   *
!**********************************************************************************
!* (C) March 2011 - Pierre Hirel                                                  *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 01 March 2017                                    *
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
USE neighbors
USE files
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE REMDOUBLES_XYZ(H,P,S,AUX,rmd_radius,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: msg
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT  !mask for atom list
INTEGER:: i, j, iat
INTEGER:: Nremoved
INTEGER,DIMENSION(:,:),ALLOCATABLE:: NeighList  !list of neighbors
REAL(dp):: distance
REAL(dp),INTENT(IN):: rmd_radius  !radius below which atoms are considered too close
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: PosList  !final positions of neighbors
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Q, T               !positions of atoms, shells (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX !auxiliary properties of atoms
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX            !auxiliary properties of atoms (temporary)
!
!
!Initialize variables
i = 0
Nremoved = 0
distance = 0.d0
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(T)) DEALLOCATE(T)
IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
!
!
CALL ATOMSK_MSG(2079,(/''/),(/rmd_radius/))
IF( rmd_radius<0.d0 ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2733,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
100 CONTINUE
!Construct neighbor list
CALL ATOMSK_MSG(11,(/""/),(/0.d0/))
CALL NEIGHBOR_LIST(H,P,rmd_radius,NeighList)
!
!Find atoms that must be removed
!For now, atoms to be removed are marked by setting their P(i,4) to zero
DO i=1,SIZE(P,1)-1  !Loop on all atoms
  IF( P(i,4)>0.1d0 ) THEN  !Ignore atoms that were already eliminated
    !
    !Find positions of neighbours of P(i,:) that are closer than rmd_radius
    CALL NEIGHBOR_POS(H,P(:,:),P(i,1:3),NeighList(i,:),ALLOCATED(NeighList),rmd_radius,PosList)
    !Now PosList(:,:) contains the cartesian positions of all neighbors in the rmd_radius,
    !their distance to the atom #i, and their indices.
    !
    !If neighboring atoms are within the rmd_radius, mark them for elimination
    !(they will be effectively removed from the array later)
    IF( ALLOCATED(PosList) .AND. SIZE(PosList,1)>0 ) THEN
      DO j=1,SIZE(PosList,1)
        iat = NINT(PosList(j,5))  !index of current neighbor
        IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(iat) ) THEN  !remove only selected atoms
          !Only consider atoms that were not removed before
          IF( P(iat,4) > 0.1d0 ) THEN
            Nremoved = Nremoved+1
            CALL ATOMSPECIES(P(iat,4),species)
            P(iat,4) = 0.d0
          ENDIF
        ENDIF
      ENDDO
    ENDIF
    !
  ENDIF
  IF(ALLOCATED(PosList)) DEALLOCATE(PosList)
ENDDO
!
!
WRITE(msg,*) 'Nremoved: ', Nremoved
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!Atoms were marked for termination, now eliminate them for good
IF(Nremoved>0) THEN
  IF( Nremoved < SIZE(P,1) ) THEN
    !Remaining atom positions will be stored in Q
    ALLOCATE( Q( SIZE(P,1)-Nremoved, SIZE(P,2) ) )
    Q(:,:) = 0.d0
    !If auxiliary properties exist they must also be modified
    IF(ALLOCATED(AUX)) THEN
      ALLOCATE( newAUX( SIZE(P,1)-Nremoved, SIZE(AUX,2) ) )
      newAUX(:,:) = 0.d0
    ENDIF
    !If shells exist they must also match P
    IF( ALLOCATED(S) .AND. SIZE(S,2)>0 ) THEN
      ALLOCATE( T( SIZE(P,1)-Nremoved, SIZE(P,2) ) )
      T(:,:) = 0.d0
    ENDIF
    !
    j = 0
    DO i=1,SIZE(P,1)
      IF(P(i,4)>0.1d0) THEN
        j = j+1
        IF( j<=SIZE(Q,1) ) THEN
          Q(j,:) = P(i,:)
          IF(ALLOCATED(AUX)) newAUX(j,:) = AUX(i,:)
          IF( ALLOCATED(S) .AND. SIZE(S,2)>0 ) T(j,:) = S(i,:)
        ELSE
          GOTO 800
        ENDIF
      ENDIF
    ENDDO
    !
    !Replace old P with new Q
    DEALLOCATE(P)
    ALLOCATE( P( SIZE(Q,1), SIZE(Q,2) ) )
    P = Q
    DEALLOCATE(Q)
    !Replace old AUX by newAUX
    IF(ALLOCATED(AUX)) THEN
      DEALLOCATE(AUX)
      ALLOCATE( AUX( SIZE(newAUX,1), SIZE(newAUX,2) ) )
      AUX = newAUX
      DEALLOCATE(newAUX)
    ENDIF
    !Replace old S with new T
    IF( ALLOCATED(S) .AND. SIZE(S,2).NE.0 ) THEN
      DEALLOCATE(S)
      ALLOCATE( S( SIZE(T,1), SIZE(T,2) ) )
      S = T
      DEALLOCATE(T)
    ENDIF
  !
  ELSE !i.e. if Nremoved>=SIZE(P,1), all atoms must be removed
    Nremoved = SIZE(P,1)
    DEALLOCATE(P)
    ALLOCATE(P(0,4))
  ENDIF
ENDIF
!
!
!
300 CONTINUE
IF( Nremoved==1 ) THEN
  CALL ATOMSK_MSG(2080,(/species/),(/DBLE(Nremoved),DBLE(SIZE(P,1))/))
ELSE
  CALL ATOMSK_MSG(2080,(/''/),(/DBLE(Nremoved),DBLE(SIZE(P,1))/))
ENDIF
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
END SUBROUTINE REMDOUBLES_XYZ
!
!
END MODULE remdoubles
