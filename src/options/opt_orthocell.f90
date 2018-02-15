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
!* Last modification: P. Hirel - 15 Feb. 2018                                     *
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
SUBROUTINE ORTHOCELL_XYZ(H,P,S,AUX)
!
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
LOGICAL:: doshells, doaux !are shells/auxiliary properties present?
LOGICAL:: new  !is the duplicated atom new?
INTEGER:: i, j, k, m, n, o
INTEGER:: lminmax
INTEGER:: NP
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
uv(:,:) = 0.d0
DO i=1,3
  uv(i,i) = 1.d12
ENDDO
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
!Check if cell vectors already form an orthogonal cell
IF( VECLENGTH(H(1,:))-DABS(H(1,1)) < 1.d-3 .AND. VECLENGTH(H(2,:))-DABS(H(2,2)) < 1.d-3 .AND. &
  & VECLENGTH(H(3,:))-DABS(H(3,3)) < 1.d-3 ) THEN
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
!Along each Cartesian direction, search for the smallest
!linear combination of vectors H that produce an orthogonal base
DO i=1,3
  IF( VECLENGTH(H(i,:))-DABS(H(i,i)) < 1.d-3 ) THEN
    !Vector H is already aligned with Cartesian axis => save it in uv
    uv(i,:) = DABS(H(i,:))
  ELSE
    !This vector is not aligned along correct Cartesian axis
    !Search for the smallest linear combination that produces a vector aligned with i
    DO m=-200,200
      DO n=-200,200
        DO o=-200,200
          vector(:) = m*H(1,:) + n*H(2,:) + o*H(3,:)
          IF( DABS(vector(i))>1.d0 .AND. VECLENGTH(vector)-DABS(vector(i)) < 1.d-3 ) THEN
            !This vector is (almost) aligned with the Cartesian axis
            IF( VECLENGTH(vector)<VECLENGTH(uv(i,:)) ) THEN
              !Vector is shorter => save it into uv
              uv(i,:) = DABS(vector(:))
            ELSEIF( VECLENGTH(vector)>3.d0*VECLENGTH(uv(i,:)) ) THEN
              !Vector is much longer => we won't find a shorter vector, exit loop
              GOTO 150
            ENDIF
          ENDIF
          !
        ENDDO !o
      ENDDO !n
    ENDDO !m
    !
  ENDIF
  150 CONTINUE
ENDDO
!
!
!
200 CONTINUE
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
!At this point we should have new cell vectors
!If not then something went wrong => abort
IF( VECLENGTH(uv(1,:))<1.d-3 .OR. VECLENGTH(uv(2,:))<1.d-3 .OR. VECLENGTH(uv(3,:))<1.d-3 .OR. &
  & VECLENGTH(uv(1,:))>1.d10 .OR. VECLENGTH(uv(2,:))>1.d10 .OR. VECLENGTH(uv(3,:))>1.d10     ) THEN
  GOTO 800
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
!
!Estimate new number of particles NP by comparing volumes of old and new cells
!Allow for +20% and +20 atoms. Actual size of arrays will be adjusted later
NP = 1.2d0*CEILING( SIZE(P,1) * DABS( DABS(uv(1,1)*uv(2,2)*uv(3,3)) / &
    & DABS(VECLENGTH(H(1,:))*VECLENGTH(H(2,:))*VECLENGTH(H(3,:))) ) ) + 20
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
!
!Set lminmax = 10 * 
lminmax = 10 * NINT( MAXVAL(uv(:,:)) / MIN(VECLENGTH(H(:,1)),VECLENGTH(H(:,2)),VECLENGTH(H(:,3))) )
!
!Loop over all replica in a wide range
NP = 0
DO i=1,SIZE(P,1)
  DO m=-lminmax,lminmax
    DO n=-lminmax,lminmax
      DO o=-lminmax,lminmax
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
              nerr = nerr+1
              CALL ATOMSK_MSG(4821,(/""/),(/DBLE(NP),DBLE(SIZE(Q,1))/))
              GOTO 1000
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
  ALLOCATE(AUX(NP,4))
  DO i=1,NP
    AUX(i,:) = newAUX(i,:)
  ENDDO
  IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
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
CALL ATOMSK_MSG(2819,(/msg/),(/0.d0/))
CALL DISPLAY_MSG(verbosity,msg,logfile)
nerr = nerr+1
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