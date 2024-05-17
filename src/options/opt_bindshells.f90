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
!* Last modification: P. Hirel - 17 May 2024                                      *
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
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: newAUXNAMES !names of auxiliary prop. (temporary)
LOGICAL:: doaux
INTEGER:: icore, ishell
INTEGER:: i, j, k, l, m, n, o
INTEGER:: Nbound !number of shells that were re-bound
INTEGER:: NP, NS !number of cores detected
INTEGER:: mass, q, qs !columns of AUX where mass, charge of cores and shells
INTEGER,DIMENSION(:),ALLOCATABLE:: Nlist !index of nearest shell
LOGICAL:: exceeds100 !does the number of neighbours exceed 100?
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: SELECT
REAL(dp):: distance  !distance between 2 particles
REAL(dp),PARAMETER:: maxCSdistance=1.5d0  !maximum allowed core-shell distance
REAL(dp),DIMENSION(3):: shift
REAL(dp),DIMENSION(4):: Ptemp, Stemp !temporary position of a core, shell
REAL(dp),DIMENSION(3,3),INTENT(IN):: H    !supercell parameters
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P     !positions of cores
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: S     !positions of shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX   !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newP, newS, newAUX  !auxiliary properties (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: V_NN  !position of nearest shell
!
!Initialize variables
doaux=.FALSE.
Nbound=0
NP=0
mass=0
q=0
qs=0
IF(ALLOCATED(newP)) DEALLOCATE(newP)
IF(ALLOCATED(newS)) DEALLOCATE(newS)
IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
IF(ALLOCATED(newAUXNAMES)) DEALLOCATE(newAUXNAMES)
!
!
msg = 'Entering BSHELLS_XYZ'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
CALL ATOMSK_MSG(2105,(/''/),(/0.d0/))
!
IF( ALLOCATED(AUXNAMES) .AND. ALLOCATED(AUX) ) THEN
  doaux=.TRUE.
  DO i=1,SIZE(AUXNAMES)
    IF( AUXNAMES(i)=="mass" ) THEN
      mass = i
    ELSEIF( AUXNAMES(i)=="q" ) THEN
      q = i
    ELSEIF( AUXNAMES(i)=="qs" ) THEN
      qs = i
    ENDIF
  ENDDO
ENDIF
!
!
!
100 CONTINUE
IF( ALLOCATED(SELECT) ) THEN
  msg = 'SELECT = .TRUE.'
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !Selected atoms must be converted into shells
  IF( .NOT.ALLOCATED(S) ) THEN
    msg = 'No shells: array S is not allocated'
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !Array S is not allocated => allocate it with the same size as P
    ALLOCATE(newS(SIZE(P,1),4))
    newS(:,:) = 0.d0
    !
    j=0
    NP=SIZE(P,1)
    DO i=1,SIZE(P,1)
      IF( IS_SELECTED(SELECT,i) ) THEN
        !This atom is selected => make it a shell
        NP = NP-1
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
          j=j+1
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
  msg = 'SELECT = .FALSE. , and array S is not allocated'
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !No selection is defined, and array S is not allocated => allocate it with the same size as P
  ALLOCATE(newP(SIZE(P,1),4))
  newP(:,:) = 0.d0
  ALLOCATE(newS(SIZE(P,1),4))
  newS(:,:) = 0.d0
  !
  IF( doaux ) THEN
    k = SIZE(AUXNAMES) !number of auxiliary properties
    IF( q>0 ) THEN
      !Ion charges are defined in AUX
      IF( qs==0 ) THEN
        !Ion charges are defined, but not shell charges
        !We will add them later: add a column in AUX
        k = k+1
        ALLOCATE(newAUXNAMES(k))
        newAUXNAMES(1:k-1) = AUXNAMES(1:k-1)
        newAUXNAMES(k) = "qs"
        DEALLOCATE(AUXNAMES)
        ALLOCATE(AUXNAMES(k))
        AUXNAMES(:) = newAUXNAMES(:)
        DEALLOCATE(newAUXNAMES)
      ENDIF
    ENDIF
    !
    ALLOCATE(newAUX(SIZE(P,1),k))
    newAUX(:,:) = 0.d0
  ENDIF
  !
  NP = 0 !counter for cores
  NS = 0 !counter for shells
  !Loop on all particles in P
  DO i=1,SIZE(P,1)-1
    l=0 !counter =0 as long as no shell is found for atom #i
    IF( P(i,4)>0.1d0 ) THEN
      !Loop on all other particles in P
      DO j=i+1,SIZE(P,1)
        !
        IF( l==0 ) THEN
          !Loop on periodic replica of particle #j
          DO o=-1,1
            DO n=-1,1
              DO m=-1,1
                shift(:) = DBLE(m)*H(1,:) + DBLE(n)*H(2,:) + DBLE(o)*H(3,:)
                Stemp(1:3) = P(j,1:3) + shift(:)  !assuming j is shell
                !Compute distance between particles i and j
                distance = VECLENGTH( Stemp(1:3) - P(i,1:3) )
                !
                IF( distance < maxCSdistance ) THEN
                  !Particles i and j are very close to each other, they form a core-shell pair
                  !Core will stay at same position, shell may be shifted using PBC
                  !Determine which one is the core and which the shell
                  l=1 !to indicate a pair was found
                  NP=NP+1 !counter for cores
                  NS=NS+1 !counter for shells
                  IF( q>0 .AND. qs>0 ) THEN
                    IF( AUX(i,q)>0.d0 ) THEN
                      !Charge of i is positive
                      IF( AUX(i,q)>AUX(j,q) .OR. AUX(i,q)>AUX(j,qs) ) THEN
                        !qi > qj => i is the core, j is the shell
                        icore = i
                        ishell = j
                      ELSE
                        !Otherwise j is the core, i is the shell
                        icore = j
                        ishell = i
                        Stemp(1:3) = P(i,1:3) - shift(:)
                      ENDIF
                    ELSE
                      !Particle i has negative charge
                      IF( AUX(j,qs)<AUX(i,q) ) THEN
                        !qi > qj => i is the core, j is the shell
                        icore = i
                        ishell = j
                      ELSE
                        !Otherwise j is the core, i is the shell
                        icore = j
                        ishell = i
                        Stemp(1:3) = P(i,1:3) - shift(:)
                      ENDIF
                    ENDIF
                  ELSEIF( mass>0 ) THEN
                    IF( AUX(i,mass)>AUX(j,mass) ) THEN
                      !i is the core, j is the shell
                      icore = i
                      ishell = j
                    ELSE
                      !j is the core, i is the shell
                      icore = j
                      ishell = i
                      Stemp(1:3) = P(i,1:3) - shift(:)
                    ENDIF
                  ELSEIF( q>0 ) THEN
                    IF( AUX(i,q)>0.d0 ) THEN
                      !i is the core, j is the shell
                      icore = i
                      ishell = j
                    ELSE
                      !Particle i has negative charge
                      IF( AUX(j,q)<AUX(i,q) ) THEN
                        icore = i
                        ishell = j
                      ELSE
                        icore = j
                        ishell = i
                        Stemp(1:3) = P(i,1:3) - shift(:)
                      ENDIF
                    ENDIF
                  ELSE
                    !No property to rely on => assume that i is the core and j the shell
                    icore = i
                    ishell = j
                  ENDIF
                  !Save position of core and shell in new arrays
                  newP(NP,:) = P(icore,:)
                  newS(NS,1:3) = Stemp(1:3)
                  newS(NS,4) = P(ishell,4)
                  !Delete data from P (set all =0, arrays will be rewritten later)
                  P(icore,:) = 0.d0
                  P(ishell,:) = 0.d0
                  IF(doaux) THEN
                    DO k=1,SIZE(AUX,2)
                      newAUX(NP,k) = AUX(icore,k)
                    ENDDO
                    IF(q>0) THEN
                      IF( qs>0 ) THEN
                        newAUX(NP,qs) = AUX(ishell,q)
                      ELSE
                        newAUX(NP,SIZE(newAUX,2)) = AUX(ishell,q)
                      ENDIF
                    ENDIF
                  ENDIF
                  WRITE(msg,*) "Bound particles # ", icore, " (core) and # ", ishell, " (shell)"
                  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
                  !We have found a core/shell pair, we can break the loop on atoms j
                  GOTO 130
                ENDIF
              ENDDO  !m
            ENDDO  !n
          ENDDO  !o
          !
        ENDIF  !end if (j>i .AND. ...)
        !
      ENDDO !end loop on j
      !
      130 CONTINUE
      IF( l==0 ) THEN
        !No particle was found close to particle #i
        !=> consider that particle #i is a core
        !Particle #j is not a shell and will be dealt with in another loop
        NP=NP+1
        newP(NP,:) = P(i,:)
        newAUX(NP,:) = AUX(i,:)
      ENDIF
      !
    ENDIF
  ENDDO !end loop on i
  !
  CALL ATOMSK_MSG(2154,(/''/),(/DBLE(NP),DBLE(NS)/))
  !
  IF( NP.NE.SIZE(P,1) ) THEN
    !Among all particles in P, some were detected to be shells
    !Resize old arrays
    IF( ALLOCATED(AUX) .AND. (ALLOCATED(AUXNAMES)) ) THEN
      DEALLOCATE(AUX)
      ALLOCATE( AUX( NP , SIZE(newAUX,2) ) )
      AUX(:,:) = 0.d0
      DO i=1,NP
        AUX(i,:) = newAUX(i,:)
      ENDDO
      DEALLOCATE(newAUX)
    ENDIF
    !
    !Re-write P so it contains only cores
    IF(ALLOCATED(P)) DEALLOCATE(P)
    ALLOCATE(P(NP,4))
    P(:,:) = 0.d0
    j=0
    DO i=1,SIZE(newP,1)
      IF( newP(i,4)>0.1d0 ) THEN
        j=j+1
        P(j,:) = newP(i,:)
      ENDIF
    ENDDO
    DEALLOCATE(newP)
    !
    !Save shells (newS) into array S
    ALLOCATE(S(NP,4))
    S(:,:) = 0.d0
    j=0
    DO i=1,SIZE(newS,1)
      IF( newS(i,4)>0.1d0 ) THEN
        j=j+1
        S(j,:) = newS(i,:)
      ENDIF
    ENDDO
    DEALLOCATE(newS)
  ENDIF
  !
ENDIF
!
!
!
200 CONTINUE
IF( ALLOCATED(S) .AND. SIZE(S,1)>0 ) THEN
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
  CALL ATOMSK_MSG(2106,(/''/),(/DBLE(Nbound)/))
  !
ENDIF
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
