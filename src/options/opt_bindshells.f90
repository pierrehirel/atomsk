MODULE bindshells
!
!**********************************************************************************
!*  BINDSHELLS                                                                    *
!**********************************************************************************
!* This module re-associates shells (in the sens of a core-shell model) with      *
!* their core.                                                                    *
!**********************************************************************************
!* (C) February 2012 - Pierre Hirel                                               *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 23 Sept. 2019                                    *
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
SUBROUTINE BSHELLS_XYZ(H,P,S,AUXNAMES,AUX,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: AUXNAMES !names of auxiliary prop.
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: newAUXNAMES !names of auxiliary prop.
INTEGER:: i, j, k, l
INTEGER:: Nbound !number of shells that were re-bound
INTEGER:: NP     !number of cores detected
INTEGER:: mass, q, qs !columns of AUX where mass, charge of cores and shells
INTEGER,DIMENSION(:),ALLOCATABLE:: Nlist !index of nearest shell
LOGICAL:: exceeds100 !does the number of neighbours exceed 100?
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: SELECT
LOGICAL,DIMENSION(:),ALLOCATABLE:: newSELECT
REAL(dp):: distance  !distance between 2 particles
REAL(dp),PARAMETER:: maxCSdistance=1.5d0  !maximum allowed core-shell distance
REAL(dp),DIMENSION(4):: Stemp !temporary position of a shell
REAL(dp),DIMENSION(3,3),INTENT(IN):: H    !supercell parameters
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P     !positions of cores
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: S     !positions of shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX   !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newP, newS, newAUX  !auxiliary properties (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: V_NN  !position of nearest shell
!
!Initialize variables
Nbound=0
mass=0
q=0
qs=0
IF(ALLOCATED(newP)) DEALLOCATE(newP)
IF(ALLOCATED(newAUX)) DEALLOCATE(newS)
IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
!
!
msg = 'Entering BSHELLS_XYZ'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
CALL ATOMSK_MSG(2105,(/''/),(/0.d0/))
!
IF( ALLOCATED(AUXNAMES) ) THEN
  DO i=1,SIZE(AUXNAMES)
    IF( AUXNAMES(i)=="mass" ) THEN
      mass = i
    ELSEIF( AUXNAMES(i)=="q" ) THEN
      q = i
    ENDIF
  ENDDO
ENDIF
!
!
!
100 CONTINUE
IF( ALLOCATED(SELECT) ) THEN
  !Selected atoms must be converted into shells
  IF( .NOT.ALLOCATED(S) ) THEN
    !Array S is not allocated => allocate it with the same size as P
    ALLOCATE(newS(SIZE(P,1),4))
    newS(:,:) = 0.d0
    !
    j=0
    DO i=1,SIZE(P,1)
      IF( SELECT(i) ) THEN
        !This atom is selected => make it a shell
        j=j+1
        newS(j,:) = P(i,:)
        !Remove it from list of atoms
        P(i,:) = 0.d0
        !IF(q>0) THEN
        
        !ENDIF
      ENDIF
    ENDDO
    !
    !Resize old arrays
    IF( ALLOCATED(AUX) .AND.(ALLOCATED(AUXNAMES)) ) THEN
      ALLOCATE(newAUX(NP,SIZE(AUXNAMES)))
      newAUX(:,:) = 0.d0
      j=0
      DO i=1,SIZE(AUX,1)
        IF( P(i,4)>0.1d0 ) THEN
          newAUX(j,:) = AUX(i,:)
        ENDIF
      ENDDO
      DEALLOCATE(AUX)
      ALLOCATE(AUX(SIZE(newAUX,1),SIZE(newAUX,2)))
      AUX(:,:) = newAUX(:,:)
      DEALLOCATE(newAUX)
    ENDIF
    !
    ALLOCATE(newP(NP,4))
    newP(:,:) = 0.d0
    j=0
    DO i=1,SIZE(P,1)
      IF( P(i,4)>0.1d0 ) THEN
        newP(j,:) = P(i,:)
      ENDIF
    ENDDO
    DEALLOCATE(P)
    ALLOCATE(P(SIZE(newP,1),4))
    P(:,:) = newP(:,:)
    DEALLOCATE(newP)
    !
    ALLOCATE(S(NP,4))
    P(:,:) = 0.d0
    j=0
    DO i=1,SIZE(newS,1)
      IF( newS(i,4)>0.1d0 ) THEN
        S(j,:) = P(i,:)
      ENDIF
    ENDDO
    DEALLOCATE(newS)
    !
  ELSE
    !Array S is already allocated, it means that new shells will be added to it
    ALLOCATE(newS(SIZE(P,1),4))
    newS(:,:) = 0.d0
    
  ENDIF
!
!
ELSEIF( .NOT.ALLOCATED(S) ) THEN
  !No selection is defined
  !Array S is not allocated => allocate it with the same size as P
  ALLOCATE(newS(SIZE(P,1),4))
  newS(:,:) = 0.d0
  !
  IF( q>0 ) THEN
    CALL RESIZE_DBLEARRAY2(AUX,SIZE(AUX,1),SIZE(AUX,2)+1,l)
    qs = SIZE(AUX,2)
  ENDIF
  !
  NP = 0
  k = 0
  !Loop on all particles in P
  DO i=1,SIZE(P,1)-1
    IF( P(i,4)>0.1d0 ) THEN
      !Loop on all other particles in P
      DO j=i+1,SIZE(P,1)
        IF( P(j,4)>0.1d0 ) THEN
          !
          !Compute distance between particles i and j
          distance = VECLENGTH( P(j,1:3) - P(i,1:3) )
          !
          IF( distance < maxCSdistance ) THEN
            !Particles i and j are very close to each other
            !Determine which one is the core and which the shell
            !and remove the shell from P (set all =0, it will be actually deleted later)
            IF( mass>0 ) THEN
              IF( AUX(i,mass)>AUX(j,mass) ) THEN
                Stemp(:) = P(i,:)
                P(i,:) = 0.d0
              ELSE
                Stemp(:) = P(j,:)
                P(j,:) = 0.d0
              ENDIF
            ELSEIF( q>0 ) THEN
              IF( AUX(i,q)>0.d0 ) THEN
                Stemp(:) = P(i,:)
                P(i,:) = 0.d0
              ELSE
                !Particle i has negative charge
                IF( AUX(j,q)<AUX(i,q) ) THEN
                  Stemp(:) = P(j,:)
                  P(j,:) = 0.d0
                ELSE
                  Stemp(:) = P(i,:)
                  P(i,:) = 0.d0
                ENDIF
              ENDIF
            ELSE
              !No property to rely on => assume that i is the core and j the shell
              Stemp(:) = P(i,:)
              P(i,:) = 0.d0
            ENDIF
            !Save particle in the list of shells
            k = k+1
            newS(k,:) = Stemp(:)
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  !
  !Resize old arrays
  IF( ALLOCATED(AUX) .AND.(ALLOCATED(AUXNAMES)) ) THEN
    ALLOCATE(newAUX(NP,SIZE(AUXNAMES)))
    newAUX(:,:) = 0.d0
    j=0
    DO i=1,SIZE(AUX,1)
      IF( P(i,4)>0.1d0 ) THEN
        newAUX(j,:) = AUX(i,:)
      ENDIF
    ENDDO
    DEALLOCATE(AUX)
    ALLOCATE(AUX(SIZE(newAUX,1),SIZE(newAUX,2)))
    AUX(:,:) = newAUX(:,:)
    DEALLOCATE(newAUX)
  ENDIF
  !
  ALLOCATE(newP(NP,4))
  newP(:,:) = 0.d0
  j=0
  DO i=1,SIZE(P,1)
    IF( P(i,4)>0.1d0 ) THEN
      newP(j,:) = P(i,:)
    ENDIF
  ENDDO
  DEALLOCATE(P)
  ALLOCATE(P(SIZE(newP,1),4))
  P(:,:) = newP(:,:)
  DEALLOCATE(newP)
  !
  ALLOCATE(S(NP,4))
  P(:,:) = 0.d0
  j=0
  DO i=1,SIZE(newS,1)
    IF( newS(i,4)>0.1d0 ) THEN
      S(j,:) = P(i,:)
    ENDIF
  ENDDO
  DEALLOCATE(newS)
  !
ENDIF
!
!
!
200 CONTINUE
DO i=1,SIZE(P,1)
  !
  IF( VECLENGTH(S(i,1:3)-P(i,1:3)) > maxCSdistance ) THEN
    !core-shell distance is very large
    !
    !First, let us assume that core #i and shell #i are associated
    !but that for some reason the shell was shifted by a combination
    !of box vectors because of PBC, and the core was not
    Stemp(:) = S(i,:)
    DO j=-3,3
      DO k=-3,3
        DO l=-3,3
          Stemp(1:3) = S(i,1:3) + DBLE(j)*H(1,:) + DBLE(k)*H(2,:) + DBLE(l)*H(3,:)
          IF( VECLENGTH(Stemp(1:3)-P(i,1:3)) <= maxCSdistance ) THEN
            !The position of this periodic image is suitable
            !Save it in S
            S(i,:) = Stemp(:)
            !One more shell was associated
            Nbound = Nbound+1
            EXIT
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
    !
    IF( VECLENGTH(S(i,1:3)-P(i,1:3)) > maxCSdistance ) THEN
      !Mere translation by box vectors was unsuccessful
      !=> shell #i does not seem to be associated with core #i
      !Parse S to see if another shell is better suited to core #i
      !
      !Search the shell that is closest to that atom
      CALL FIND_NNN(H,S,P(i,:),1,V_NN,Nlist,exceeds100)
      !V_NN(1,:) = position of closest shell
      !Nlist(1) = index of that shell
      !
      !Check that we found a shell that has a different index than i,
      !and that it is actually a shell
      IF( Nlist(1).NE.i .AND. NINT(S(Nlist(1),4)).NE.0 ) THEN
        !This shell was not associated with this core before
        IF( VECLENGTH( V_NN(1,1:3)-P(i,1:3) ) <= maxCSdistance ) THEN
          !If distance is smaller than 1.5 A,then
          !this shell belongs to that core
          IF( NINT(S(i,4)).NE.0 ) THEN
            !Present core already had a shell: save this shell position
            Stemp(:) = S(i,:)
          ENDIF
          !
          !Save closest shell position at current index i
          S(i,1:3) = V_NN(1,1:3)
          S(i,4) = P(i,4)
          !
          IF( NINT(Stemp(4)).NE.0 ) THEN
            !Save previous shell of atom i to index Nlist(1)
            S(Nlist(1),:) = Stemp(:)
            IF( VECLENGTH( S(Nlist(1),1:3)-P(Nlist(1),1:3) ) <= maxCSdistance ) THEN
              !This shell belongs to the core at Nlist(1)
              Nbound = Nbound+1
            ENDIF
          ELSE
            !Otherwise, simply delete previous shell position
            S(Nlist(1),:) = 0.d0
          ENDIF
          !
          !One more shell was associated
          Nbound = Nbound+1
        ELSE
          !Otherwise (distance>1.5A) this core has no shell
          S(i,:) = 0.d0
        ENDIF
      ENDIF
    ENDIF
    !
    !
  ENDIF
ENDDO
!
!
CALL ATOMSK_MSG(2106,(/''/),(/DBLE(Nbound)/))
!
!
!
1000 CONTINUE
IF(ALLOCATED(Nlist)) DEALLOCATE(Nlist)
IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
!
!
!
END SUBROUTINE BSHELLS_XYZ
!
!
!
END MODULE bindshells
