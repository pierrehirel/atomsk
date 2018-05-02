MODULE orthocell
!
!**********************************************************************************
!*  ORTHOCELL                                                                     *
!**********************************************************************************
!* Provided a cell with arbitrary vectors H and containing atom positions P,      *
!* this module finds the smallest orthogonal cell that maintains periodic         *
!* boundary conditions.                                                           *
!**********************************************************************************
!* (C) Feb. 2018 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@izbs.uni-karlsruhe.de                                         *
!* Last modification: P. Hirel - 24 April 2018                                    *
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
USE resize
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE ORTHOCELL_XYZ(H,P,S,AUX,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
LOGICAL:: doshells, doaux !are shells/auxiliary properties present?
LOGICAL:: new  !is the duplicated atom new?
LOGICAL,DIMENSION(3):: aligned, reversed !is box vector aligned along X,Y,Z? Did we just reverse it?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT, newSELECT  !mask for atom list
INTEGER:: i, j, k, m, n, o
INTEGER:: mminmax, nminmax, ominmax, lminmax
INTEGER:: mfinal, nfinal, ofinal
INTEGER:: NP
REAL(dp):: vlen  !length of a vector
REAL(dp),DIMENSION(3):: vector !just a vector
REAL(dp),DIMENSION(4):: tempP  !temporary atom position and species
REAL(dp),DIMENSION(3,3):: uv   !new vectors along each Cartesian direction
REAL(dp),DIMENSION(3,3):: H    !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, S, Q, T !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX   !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX              !auxiliary properties (temporary)
!
!Initialize variables
doshells = .FALSE.
doaux = .FALSE.
aligned(:) = .FALSE.
reversed(:) = .FALSE.
mminmax = 10
nminmax = 10
ominmax = 10
uv(:,:) = 1.d12
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(T)) DEALLOCATE(T)
!
WRITE(msg,*) 'Entering ORTHOCELL_XYZ'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
CALL ATOMSK_MSG(2143,(/msg/),(/0.d0/))
!
!
!
100 CONTINUE
IF( verbosity>=4 ) THEN
  !Print some debug messages
  WRITE(msg,*) "Old cell vectors:"
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,*) "  H(1,:) = ", H(1,:)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,*) "  H(2,:) = ", H(2,:)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,*) "  H(3,:) = ", H(3,:)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ENDIF
!
!Check if cell vectors are already aligned with Cartesian axes
IF( DABS(VECLENGTH(H(1,:))-DABS(H(1,1))) < 1.d-3 .AND. H(1,1)>0.d0 .AND. &
  & DABS(VECLENGTH(H(2,:))-DABS(H(2,2))) < 1.d-3 .AND. H(2,2)>0.d0 .AND. &
  & DABS(VECLENGTH(H(3,:))-DABS(H(3,3))) < 1.d-3 .AND. H(3,3)>0.d0        ) THEN
  !Cell vectors already orthogonal => skip to the end
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2760,(/msg/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Verify if shells are present
IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
  doshells = .TRUE.
ENDIF
!Verify if auxiliary properties are present
IF( ALLOCATED(AUX) .AND. SIZE(AUX,1)==SIZE(P,1) ) THEN
  doaux = .TRUE.
ENDIF
!
110 CONTINUE
!Along each Cartesian direction, search for the smallest
!linear combination of vectors H that produce an orthogonal base
!NOTE: the triple loop is first parsed with small values of mminmax, nminmax and ominmax
!  This should ensure a fast search for simple systems/orientations
!  If it fails, these values are increased and the loops are parsed again (see label 200 below)
DO i=1,3
  mfinal = 0
  nfinal = 0
  ofinal = 0
  IF( i==1 ) THEN
    j = 2
    k = 3
  ELSEIF( i==2 ) THEN
    j = 3
    k = 1
  ELSE
    j = 1
    k = 2
  ENDIF
  !
  IF( DABS(VECLENGTH(H(i,:))-DABS(H(i,i))) < 1.d-3 ) THEN
    !Vector H is already aligned with Cartesian axis => save it in uv
    aligned(i) = .TRUE.
    IF( H(i,i)<0.d0 ) THEN
      !The new cell vector will just be uv(i,:) = -H(i,:)
      reversed = .TRUE.
    ENDIF
    uv(i,:) = DABS(H(i,:))
    WRITE(msg,'(a3,i1,a17)') "uv(", i, ") already aligned"
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  ELSE
    !This vector is not aligned along correct Cartesian axis
    !Search for the smallest linear combination that produces a vector aligned with i
    DO m=-mminmax,mminmax
      DO n=-nminmax,nminmax
        DO o=-ominmax,ominmax
          vector(:) = m*H(1,:) + n*H(2,:) + o*H(3,:)
          vlen = VECLENGTH(vector)
          IF( DABS(vector(j))<1.d-3 .AND. DABS(vector(k))<1.d-3 .AND. vlen>1.d0 ) THEN
            !This vector is (almost) aligned with the Cartesian axis
            IF( vlen < VECLENGTH(uv(i,:)) ) THEN
              !Vector is shorter => save it into uv
              uv(i,:) = DABS(vector(:))
              mfinal = m
              nfinal = n
              ofinal = o
            !ELSEIF( VECLENGTH(vector) > 3.d0*VECLENGTH(uv(i,:)) ) THEN
              !Vector is much longer => we won't find a shorter vector, exit loop
             ! GOTO 150
            ENDIF
          ELSEIF( DABS(vector(j))<1.d-1 .AND. DABS(vector(k))<1.d-1 .AND. vlen>1.d0 ) THEN
            IF( DABS(vector(j))<uv(i,j) .AND. DABS(vector(k))<uv(i,k) ) THEN
              !This vector is better aligned with Cartesian axis than previously found vector
              !(although this may mean that it is longer)
              uv(i,:) = DABS(vector(:))
              mfinal = m
              nfinal = n
              ofinal = o
            ENDIF
          ENDIF
          !
        ENDDO !o
      ENDDO !n
    ENDDO !m
    150 CONTINUE
    WRITE(msg,'(a3,i1,a4,3(i4,a6),3e16.6)') "uv(", i, ") = ", mfinal, "*H1 + ", nfinal, "*H2 + ", ofinal, "*H3 = ", uv(i,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
  ENDIF
ENDDO
!
!
!
200 CONTINUE
!At this point we should have new cell vectors
IF( VECLENGTH(uv(1,:))<1.d-3 .OR. VECLENGTH(uv(2,:))<1.d-3 .OR. VECLENGTH(uv(3,:))<1.d-3 .OR. &
  & VECLENGTH(uv(1,:))>1.d10 .OR. VECLENGTH(uv(2,:))>1.d10 .OR. VECLENGTH(uv(3,:))>1.d10     ) THEN
  !We don't have appropriate cell vectors
  IF( mminmax+nminmax+ominmax < 80 ) THEN
    !The first attempt with small loops did not work => try to expand the loop
    mminmax = 200
    nminmax = 200
    ominmax = 200
    GOTO 110
  ELSEIF( mminmax+nminmax+ominmax < 620 ) THEN
    !The second attempt with larger loops did not work => try to expand the loop again
    !(this may take a long time)
    CALL ATOMSK_MSG(3,(/""/),(/0.d0/))
    mminmax = 500
    nminmax = 500
    ominmax = 500
    GOTO 110
  ELSE
    !Even the larger loops did not work => abort
    GOTO 800
  ENDIF
ENDIF
!
!Make sure vectors are oriented positively along each axis
!and that non-diagonal elements are all zero
DO i=1,3
  DO j=1,3
    IF( i==j ) THEN
      uv(i,j) = DABS(uv(i,j))
    ELSE
      uv(i,j) = 0.d0
    ENDIF
  ENDDO
ENDDO
IF( verbosity>=4 ) THEN
  !Print some debug messages
  WRITE(msg,*) "New cell vectors:"
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,*) "  uv(1,:) = ", uv(1,:)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,*) "  uv(2,:) = ", uv(2,:)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,*) "  uv(3,:) = ", uv(3,:)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ENDIF
!
IF( aligned(1) .AND. aligned(2) .AND. aligned(3) ) THEN
  !Very special case: all box vectors were already aligned,
  !but maybe some of them pointed towards negative directions
  !In this case, the number of atoms will remain the same,
  !but atoms will be shifted by the box vector
  DO i=1,3
    IF( reversed(i) ) THEN
      !This box vector was just reversed => shift all atoms by this vector
      P(:,i) = P(:,i) + uv(i,i)
    ENDIF
  ENDDO
  !
ELSE
  !Estimate new number of particles NP by comparing volumes of old and new cells
  !Allow for +50% and +25 atoms. Actual size of arrays will be adjusted later
  NP = 1.5d0*CEILING( SIZE(P,1) * DABS( DABS(uv(1,1)*uv(2,2)*uv(3,3)) / &
      & DABS(VECLENGTH(H(1,:))*VECLENGTH(H(2,:))*VECLENGTH(H(3,:))) ) ) + 25
  WRITE(msg,*) "Estimated new number of atoms : ", NP
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  IF(ALLOCATED(Q)) DEALLOCATE(Q)
  ALLOCATE( Q(NP,4) )
  Q(:,:) = 0.d0
  IF(doshells) THEN
    ALLOCATE( T(NP,4) )
    T(:,:) = 0.d0
  ENDIF
  IF(doaux) THEN
    ALLOCATE( newAUX(NP,SIZE(AUX,2)) )
    newAUX(:,:) = 0.d0
  ENDIF
  !Verify if some atoms are selected
  IF( ALLOCATED(SELECT) ) THEN
    ALLOCATE( newSELECT(NP) )
    newSELECT(:) = .FALSE.
  ENDIF
  !
  !Set min/max replica search
  IF( mminmax.NE.1 ) THEN
    mminmax = 10 * NINT( MAXVAL(uv(:,:)) / MIN(VECLENGTH(H(:,1)),VECLENGTH(H(:,2)),VECLENGTH(H(:,3))) )
  ENDIF
  IF( nminmax.NE.1 ) THEN
    nminmax = 10 * NINT( MAXVAL(uv(:,:)) / MIN(VECLENGTH(H(:,1)),VECLENGTH(H(:,2)),VECLENGTH(H(:,3))) )
  ENDIF
  IF( ominmax.NE.1 ) THEN
    ominmax = 10 * NINT( MAXVAL(uv(:,:)) / MIN(VECLENGTH(H(:,1)),VECLENGTH(H(:,2)),VECLENGTH(H(:,3))) )
  ENDIF
  WRITE(msg,*) "Atom duplication min/max = ", mminmax, nminmax, ominmax
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !
  !Loop over all replica in a wide range
  NP = 0
  DO i=1,SIZE(P,1)
    DO m=-mminmax,mminmax
      DO n=-nminmax,nminmax
        DO o=-ominmax,ominmax
          !Compute cartesian position of this replica
          tempP(1) = P(i,1) + DBLE(m)*H(1,1) + DBLE(n)*H(2,1) + DBLE(o)*H(3,1)
          tempP(2) = P(i,2) + DBLE(m)*H(1,2) + DBLE(n)*H(2,2) + DBLE(o)*H(3,2)
          tempP(3) = P(i,3) + DBLE(m)*H(1,3) + DBLE(n)*H(2,3) + DBLE(o)*H(3,3)
          tempP(4) = P(i,4)
          IF( tempP(1)>-1.d-12 .AND. tempP(1)<=uv(1,1)-1.d-12 .AND.             &
            & tempP(2)>-1.d-12 .AND. tempP(2)<=uv(2,2)-1.d-12 .AND.             &
            & tempP(3)>-1.d-12 .AND. tempP(3)<=uv(3,3)-1.d-12       ) THEN
            !This replica is inside the new cell, mark it as new
            new = .TRUE.
            !Verify that its position is different from all previous atoms
            DO k=1,NP
              IF( VECLENGTH(tempP(1:3)-Q(k,1:3)) < 0.1d0 ) THEN
                new = .FALSE.
              ENDIF
            ENDDO
            !
            IF( new ) THEN
              NP = NP+1
              IF(NP>SIZE(Q,1)) THEN
                !Resize array Q
                CALL RESIZE_DBLEARRAY2(Q,SIZE(Q,1)+10,SIZE(Q,2))
              ENDIF
              Q(NP,:) = tempP(:)
              !
              IF( doshells ) THEN
                !Compute position of the replica of this shell
                tempP(1) = S(i,1) + DBLE(m)*H(1,1) + DBLE(n)*H(2,1) + DBLE(o)*H(3,1)
                tempP(2) = S(i,2) + DBLE(m)*H(1,2) + DBLE(n)*H(2,2) + DBLE(o)*H(3,2)
                tempP(3) = S(i,3) + DBLE(m)*H(1,3) + DBLE(n)*H(2,3) + DBLE(o)*H(3,3)
                tempP(4) = S(i,4)
                T(NP,:) = tempP(:)
              ENDIF
              !
              IF( doaux ) THEN
                newAUX(NP,:) = AUX(i,:)
              ENDIF
              !
              IF( ALLOCATED(SELECT) ) THEN
                newSELECT(NP) = SELECT(i)
              ENDIF
              !
            ENDIF
          ENDIF
          !
        ENDDO !o
      ENDDO !n
    ENDDO !m
  ENDDO !i
  !
  !Replace old P with the new Q
  IF(ALLOCATED(P)) DEALLOCATE(P)
  ALLOCATE(P(NP,4))
  DO i=1,NP
    P(i,:) = Q(i,:)
  ENDDO
  IF(ALLOCATED(Q)) DEALLOCATE(Q)
  !Same with shell if they are present
  IF( doshells ) THEN
    IF(ALLOCATED(S)) DEALLOCATE(S)
    ALLOCATE(S(NP,4))
    DO i=1,NP
      S(i,:) = T(i,:)
    ENDDO
    IF(ALLOCATED(T)) DEALLOCATE(T)
  ENDIF
  !Same with auxiliary properties if they are present
  IF( doaux ) THEN
    IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
    ALLOCATE( AUX(NP,SIZE(newAUX,2)) )
    DO i=1,NP
      AUX(i,:) = newAUX(i,:)
    ENDDO
    IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
  ENDIF
  !Same with SELECT array if it is defined
  IF( ALLOCATED(SELECT) ) THEN
    DEALLOCATE(SELECT)
    ALLOCATE(SELECT(NP))
    SELECT(:) = .FALSE.
    DO i=1,NP
      SELECT(i) = newSELECT(i)
    ENDDO
    IF(ALLOCATED(newSELECT)) DEALLOCATE(newSELECT)
  ENDIF
  !
ENDIF
!
!Save new cell vectors to H (making sure non-diagonal components are all zero)
H(:,:) = 0.d0
DO i=1,3
  H(i,i) = uv(i,i)
ENDDO
!
!
CALL ATOMSK_MSG(2144,(/msg/),(/DBLE(SIZE(P,1))/))
GOTO 1000
!
!
!
800 CONTINUE
nerr = nerr+1
CALL ATOMSK_MSG(2819,(/msg/),(/0.d0/))
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE ORTHOCELL_XYZ
!
!
!
END MODULE orthocell