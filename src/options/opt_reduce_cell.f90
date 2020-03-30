MODULE reduce_cell
!
!**********************************************************************************
!*  REDUCE_CELL                                                                   *
!**********************************************************************************
!* This module takes in a cell and atomic coordinates, and attempts to reduce     *
!* its size to provide a "unit cell".                                             *
!**********************************************************************************
!* (C) March 2020 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 30 March 2020                                    *
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
!
CONTAINS
!
!
SUBROUTINE REDUCECELL(H,P,S,AUX,SELECT)
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: SELECT  !mask for atom list
LOGICAL,DIMENSION(:),ALLOCATABLE:: newSELECT             !mask for atom list (temporary)
INTEGER,DIMENSION(3):: reduce
INTEGER:: i, j, k, NP
REAL(dp):: vp, tempreal
REAL(dp),DIMENSION(3):: V, cp
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H      !box vectors
REAL(dp),DIMENSION(3,3):: Hnew
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S  !positions of cores, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Q, T                !positions of cores, shells (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX   !auxiliary properties of atoms
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX              !auxiliary properties of atoms (temporary)
!
!
!
!Initialize variables
i = 0
reduce(:) = 0
Hnew(:,:) = 1.d12
IF(ALLOCATED(newSELECT)) DEALLOCATE(newSELECT)
IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(T)) DEALLOCATE(T)
!
!
!
100 CONTINUE
CALL ATOMSK_MSG(2149,(/''/),(/0.d0/))
NP=0
!
! Loop on all atoms
DO i=2,SIZE(P,1)
  !
  IF( P(i,4) == P(1,4) ) THEN
    !Initialize variables
    k = 0
    vp = 0.d0
    !
    !Atom #i is the same spieces as atom #1
    !Compute vector between those two atoms
    V(:) = P(i,1:3) - P(1,1:3)
    !
    !Check if this vector is aligned along a cell vector
    DO j=1,3
      !Compute vector product and cross product
      vp = DOT_PRODUCT( H(j,:) , V(:) )
      cp(:) = CROSS_PRODUCT( H(j,:) , V(:) )
      IF( DABS(vp)>1.d0 .AND. VECLENGTH(cp)<0.1d0 ) THEN
        !Vector V is aligned with a cell vector
        k = j
        EXIT
      ENDIF
    ENDDO
    !
    IF( k>0 ) THEN
      !Vector V is aligned with cell vector #aligned
      !Check if it is shorter than cell vector, AND shorter than previous vector Hnew
      IF( VECLENGTH(V(:)) < VECLENGTH(H(k,:)) .AND. VECLENGTH(V(:)) < VECLENGTH(Hnew(k,:)) ) THEN
        !This vector is shorter than initial cell vector
        !Save it as new cell vector
        Hnew(k,:) = V(:)
      ENDIF
    ENDIF
  ENDIF
  !
ENDDO
!
!
!
200 CONTINUE
!Replace old cell vectors with new ones
reduce(:)=0
DO i=1,3
  IF( VECLENGTH(Hnew(i,:)) < VECLENGTH(H(i,:)) ) THEN
    H(i,:) = H(i,:)*VECLENGTH(Hnew(i,:)) / VECLENGTH(H(i,:))
    reduce(i) = 1
    !Eliminate approximations in cell vectors components
    DO j=1,3
      IF( DABS(H(i,j)) < 1.d-12 ) H(i,j) = 0.d0
    ENDDO
  ENDIF
ENDDO
!
IF( ANY( reduce(:)>0 ) ) THEN
  !One or more cell vectors were changed: remove all atoms that are not inside the new cell
  !Convert atom positions into reduced coordinates
  CALL CART2FRAC(P,H)
  IF( ALLOCATED(S) ) THEN
    CALL CART2FRAC(S,H)
  ENDIF
  !
  !Allocate temporary arrays
  ALLOCATE( Q( SIZE(P,1),4 ) )
  Q(:,:) = 0.d0
  IF(ALLOCATED(SELECT)) ALLOCATE( newSELECT( SIZE(P,1) ) )
  IF(ALLOCATED(S)) ALLOCATE( T( SIZE(P,1),4 ) )
  IF(ALLOCATED(AUX)) ALLOCATE( newAUX( SIZE(AUX,1), SIZE(AUX,2) ) )
  !
  !Copy positions of atoms that are inside the new cell into temporary arrays (ignore atoms that are outside)
  !(also copy their properties if any are defined)
  NP = 0
  DO i=1,SIZE(P,1)
    IF( P(i,1)>-1.d-12 .AND. P(i,2)>-1.d-12 .AND. P(i,3)>-1.d-12 .AND.      &
      & P(i,1)<1.d0-1.d-12 .AND. P(i,2)<1.d0-1.d-12 .AND. P(i,3)<1.d0-1.d-12   ) THEN
      NP = NP+1
      Q(NP,:) = P(i,:)
      IF(ALLOCATED(S)) T(NP,:) = S(i,:)
      IF(ALLOCATED(AUX)) newAUX(NP,:) = AUX(i,:)
      IF(ALLOCATED(SELECT)) newSELECT(NP) = SELECT(i)
    ENDIF
  ENDDO
  !
  !Re-allocate arrays into their final sizes
  DEALLOCATE(P)
  ALLOCATE(P(NP,4))
  IF( ALLOCATED(S) ) THEN
    DEALLOCATE(S)
    ALLOCATE(S(NP,4))
  ENDIF
  IF( ALLOCATED(AUX) ) THEN
    j = SIZE(AUX,2)
    DEALLOCATE(AUX)
    ALLOCATE(AUX(NP,j))
  ENDIF
  IF( ALLOCATED(SELECT) ) THEN
    DEALLOCATE(SELECT)
    ALLOCATE(SELECT(NP))
  ENDIF
  !
  !Copy data into final arrays
  DO i=1,NP
    P(i,:) = Q(i,:)
    IF(ALLOCATED(S)) S(i,:) = T(i,:)
    IF(ALLOCATED(AUX)) AUX(i,:) = newAUX(i,:)
    IF(ALLOCATED(SELECT)) SELECT(i) = newSELECT(i)
  ENDDO
  !
  !Free temporary arrays
  DEALLOCATE(Q)
  IF(ALLOCATED(T)) DEALLOCATE(T)
  IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
  IF(ALLOCATED(newSELECT)) DEALLOCATE(newSELECT)
  !
  !Restore coordinates back into Cartesian
  CALL FRAC2CART(P,H)
  IF( ALLOCATED(S) ) THEN
    CALL FRAC2CART(S,H)
  ENDIF
  !
  CALL ATOMSK_MSG(2150,(/''/),(/ DBLE(reduce(1)), DBLE(reduce(2)), DBLE(reduce(3)), DBLE(NP) /))
  !
ELSE
  !No cell vector was modified: display warning and go on
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2764,(/''/),(/0.d0/))
  !
ENDIF
!
!
!
1000 CONTINUE
!
!
!
END SUBROUTINE REDUCECELL
!
END MODULE reduce_cell

