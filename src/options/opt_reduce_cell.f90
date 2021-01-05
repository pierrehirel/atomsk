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
!* Last modification: P. Hirel - 05 Jan. 2021                                     *
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
SUBROUTINE REDUCECELL(H,P,S,AUX,dir,SELECT)
!
IMPLICIT NONE
CHARACTER(LEN=1),INTENT(IN):: dir  !direction to reduce (if empty, reduce all directions)
CHARACTER(LEN=128):: msg
LOGICAL,DIMENSION(3):: reducedir !=TRUE if this direction must be reduced
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: SELECT  !mask for atom list
LOGICAL,DIMENSION(:),ALLOCATABLE:: newSELECT             !mask for atom list (temporary)
INTEGER,DIMENSION(3):: reduced !=1 if cell was reduced along X,Y,Z , =0 otherwise
INTEGER:: i, iref, j, k, NP, m, n, o
REAL(dp):: cp, vp, tempreal
REAL(dp),PARAMETER:: rtol=1.d-3  !relative tolerance on atom positions
REAL(dp),DIMENSION(3):: tol      !tolerance on atom positions along X,Y,Z (angströms)
REAL(dp),DIMENSION(3):: V
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
NP=0
tol(:) = 0.d0
IF(ALLOCATED(newSELECT)) DEALLOCATE(newSELECT)
IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(T)) DEALLOCATE(T)
!
reducedir(:) = .FALSE.
IF( dir=="x" .OR. dir=="X" ) THEN
  reducedir(1) = .TRUE.
ELSEIF( dir=="y" .OR. dir=="Y" ) THEN
  reducedir(2) = .TRUE.
ELSEIF( dir=="z" .OR. dir=="Z" ) THEN
  reducedir(3) = .TRUE.
ELSE
  reducedir(:) = .TRUE.
ENDIF
!
!
!
100 CONTINUE
CALL ATOMSK_MSG(2149,(/''/),(/0.d0/))
!By default, use atom closest to the origin (0,0,0) as "reference" and search for its periodic images
iref = 1
DO i=2,SIZE(P,1)
  IF( VECLENGTH(P(i,1:3)) < VECLENGTH(P(iref,1:3)) ) THEN
    iref = i
  ENDIF
ENDDO
!
!
120 CONTINUE
reduced(:) = 0
Hnew(:,:) = 1.d12
IF( verbosity==4 ) THEN
  WRITE(msg,*) "Reference atom #", iref, ": ", P(iref,1), P(iref,2), P(iref,3)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ENDIF
!
! Loop on all atoms to find an atom equivalent to atom #iref
! NOTE: here "equivalent" just means it has the same species.
!      The environment of atoms is not checked
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,k,cp,vp,V)
DO i=1,SIZE(P,1)
  !
  IF( P(i,4)==P(iref,4) .AND. i.NE.iref ) THEN
    !Atom #i is the same spieces as atom #iref
    !
    !Initialize variables
    j = 0
    k = -1
    cp = 0.d0
    vp = 0.d0
    !
    !Compute vector between those two atoms
    V(:) = P(i,1:3) - P(iref,1:3)
    !
    !Check if this vector is aligned along a cell vector
    DO WHILE(j<3 .AND. k.NE.j)
      j = j+1
      !Compute vector product and cross product
      vp = DOT_PRODUCT( H(j,:) , V(:) )
      cp = VECLENGTH( CROSS_PRODUCT( H(j,:) , V(:) ) )
      IF( DABS(vp)>1.d0 .AND. cp<0.1d0 ) THEN
        !Vector V is aligned with a cell vector
        k = j
      ENDIF
    ENDDO
    !
    IF( k>0 ) THEN
      !Vector V is aligned with cell vector #k
      !Check if it is shorter than cell vector, AND shorter than previous vector Hnew
      IF( VECLENGTH(V(:)) < VECLENGTH(H(k,:)) .AND. VECLENGTH(V(:)) < VECLENGTH(Hnew(k,:)) ) THEN
        !This vector is shorter than initial cell vector
        !Save it as new cell vector
        !$OMP CRITICAL
        Hnew(k,:) = V(:)
        reduced(k) = 1
        !$OMP END CRITICAL
      ENDIF
    ENDIF
  ENDIF
  !
ENDDO
!$OMP END PARALLEL DO
!
!
150 CONTINUE
IF( .NOT.ANY(reduced(:)==1) .AND. SIZE(P,1)>100 .AND. iref==1 ) THEN
  !Unable to find a reduced cell using first atom
  !System has a lot of atoms, itmay be complex (e.g. defect at the origin), try with a different reference atom
  iref = SIZE(P,1) / 5
  WRITE(msg,*) "No new cell found, trying with atom #", iref, " as reference"
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  IF(iref>1) GOTO 120
ENDIF
!
!
!
200 CONTINUE
!Replace old cell vectors with new ones
reduced(:)=0  !so far the cell was not reduced
!Loop on X, Y, Z
DO i=1,3
  !Check if user asked to reduce cell along this direction
  IF( reducedir(i) ) THEN
    IF( VECLENGTH(Hnew(i,:)) < VECLENGTH(H(i,:)) .AND. VECLENGTH(Hnew(i,:)) > 0.5d0 ) THEN
      !Replace cell vector H(i,:) with new one
      H(i,:) = H(i,:)*VECLENGTH(Hnew(i,:)) / VECLENGTH(H(i,:))
      !Mark it as modified
      reduced(i) = 1
      !Eliminate approximations in cell vectors components
      DO j=1,3
        IF( DABS(H(i,j)) < rtol ) H(i,j) = 0.d0
      ENDDO
    ENDIF
  ENDIF
ENDDO
!
IF( verbosity==4 ) THEN
  msg = "New cell vectors:"
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  DO i=1,3
    WRITE(msg,'(6X,3f12.6)') H(i,1), H(i,2), H(i,3)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  ENDDO
ENDIF
!
!Set tolerance on atom position along X, Y, Z
DO i=1,3
  tol(i) = MIN( 1.d-3 , rtol / VECLENGTH(H(i,:)) )
ENDDO
!
IF( ANY( reduced(:)>0 ) ) THEN
  !One or more cell vectors were changed: remove all atoms that are not inside the new cell
  !Convert atom positions into reduced coordinates
  CALL CART2FRAC(P,H)
  IF( ALLOCATED(S) ) THEN
    CALL CART2FRAC(S,H)
  ENDIF
  !
  !Allocate temporary arrays
  !NOTE: cell is reduced so there should be FEWER atoms in the final cell
  !Actual array sizes will be adjusted later
  ALLOCATE( Q( SIZE(P,1),4 ) )
  Q(:,:) = 0.d0
  IF(ALLOCATED(SELECT)) ALLOCATE( newSELECT( SIZE(P,1) ) )
  IF(ALLOCATED(S)) ALLOCATE( T( SIZE(P,1),4 ) )
  IF(ALLOCATED(AUX)) ALLOCATE( newAUX( SIZE(AUX,1), SIZE(AUX,2) ) )
  !
  !Copy positions of atoms that are inside the new cell into temporary arrays (ignore atoms that are outside)
  !(also copy their properties if any are defined)
  !Place reference atom at origin (0,0,0) and copy its properties
  NP = 1
  Q(1,1:3) = 0.d0
  Q(1,4) = P(iref,4)
  IF(ALLOCATED(S)) THEN
    T(1,1:3) = S(iref,1:3) - P(iref,1:3)
    T(1,4) = S(iref,4)
  ENDIF
  IF(ALLOCATED(AUX)) newAUX(1,:) = AUX(iref,:)
  IF(ALLOCATED(SELECT)) newSELECT(1) = SELECT(iref)
  !Loop on all atoms
  DO i=1,SIZE(P,1)
    !Only consider atoms that are different from reference atom
    IF( i.NE.iref ) THEN
      !Compute relative position of this atom
      V(:) = P(i,1:3) - P(iref,1:3)
      !Check if this atom is in the new box
      !Along each direction, atom is counted as being inside the box if
      !the box was not reduced along that direction, or if atom is within the new boundaries
      IF( ( reduced(1)==0 .OR. (V(1)>-tol(1) .AND. V(1)<=1.d0-tol(1)) ) .AND. &
        & ( reduced(2)==0 .OR. (V(2)>-tol(2) .AND. V(2)<=1.d0-tol(2)) ) .AND. &
        & ( reduced(3)==0 .OR. (V(3)>-tol(3) .AND. V(3)<=1.d0-tol(3)) )       ) THEN
        !Yes, it is in the box: add it to the list
        NP = NP+1
        IF(NP<=SIZE(Q,1)) THEN
          Q(NP,1:3) = V(:)
          Q(NP,4) = P(i,4)
          IF(ALLOCATED(S)) THEN
            T(NP,1:3) = S(i,1:3) - P(iref,1:3)
            T(NP,4) = S(i,4)
          ENDIF
          IF(ALLOCATED(AUX)) newAUX(NP,:) = AUX(i,:)
          IF(ALLOCATED(SELECT)) newSELECT(NP) = SELECT(i)
        ENDIF
      ENDIF
      !
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
  CALL ATOMSK_MSG(2150,(/''/),(/ DBLE(reduced(1)), DBLE(reduced(2)), DBLE(reduced(3)), DBLE(NP) /))
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

