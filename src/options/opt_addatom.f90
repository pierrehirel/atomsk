MODULE addatom
!
!**********************************************************************************
!*  ADDATOM                                                                       *
!**********************************************************************************
!* This module adds a new atom of the given species, either at the given          *
!* position, of near an existing atom, or N new atoms at random places.           *
!**********************************************************************************
!* (C) March 2014 - Pierre Hirel                                                  *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 28 Sep. 2016                                     *
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
USE subroutines
USE neighbors
!
!
CONTAINS
!
!
SUBROUTINE ADDATOM_XYZ(H,P,S,AUX,addatom_species,addatom_type,addatom_prop,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=2):: species
CHARACTER(LEN=2),INTENT(IN):: addatom_species !species of atom(s) to add
CHARACTER(LEN=8),INTENT(IN):: addatom_type    !"at" or "near" or "random"
CHARACTER(LEN=128):: msg
LOGICAL:: exceeds100 !does the number of neighbors exceed 100?
LOGICAL:: hasShells  !does this type of atom have shells?
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT  !mask for atom list
INTEGER:: atomindex
INTEGER:: NP !number of particles
INTEGER:: i, j, k, m, n
INTEGER:: addedatoms !number of atoms added
INTEGER,DIMENSION(:),ALLOCATABLE:: Nlist  !list of indices of neighbors
INTEGER,DIMENSION(:,:),ALLOCATABLE:: NeighList !list of index of neighbors
REAL(dp):: distance, distance2, dmax
REAL(dp):: snumber !atomic number of the new atoms
REAL(dp):: x, y, z
REAL(dp),DIMENSION(4),INTENT(IN):: addatom_prop  !properties of atom(s) to add
                                                 !if addatom_type=="at", position x,y,z of new atom
                                                 !if addatom_type=="relative", index of atom and x,y,z
                                                 !if addatom_type=="near", index of atom
                                                 !if addatom-type=="random", number of atoms to add
REAL(dp),DIMENSION(3):: V !a vector
REAL(dp),DIMENSION(:),ALLOCATABLE:: randarray  !array for storing random numbers
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !vectors of supercell
REAL(dp),DIMENSION(4,3):: atomsite  !positions of the 4 atoms forming the interstitial site
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S  !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newP, newS          !positions of atoms, shells (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX   !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX              !auxiliary properties (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: V_NN                !positions of neighbors
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: PosList             !positions of neighbors
!
species = ''
exceeds100 = .FALSE.
hasShells = .FALSE.
i = 0
NP = 0
addedatoms = 0
snumber = 0.d0
x = 0.d0
y = 0.d0
z = 0.d0
IF(ALLOCATED(newP)) DEALLOCATE(newP)
IF(ALLOCATED(newS)) DEALLOCATE(newS)
IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
!
WRITE(msg,*) 'Entering ADDATOM_XYZ: '//TRIM(ADJUSTL(addatom_species))//","// &
             & TRIM(ADJUSTL(addatom_type))
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!
100 CONTINUE
CALL ATOMSK_MSG(2116,(/addatom_species//"      ",addatom_type/),addatom_prop(:))
!
!Convert atom species into atomic number
CALL ATOMNUMBER(addatom_species,snumber)
IF( snumber < 1.d-12 ) THEN
  !Unrecognized species => abandon ship
  CALL ATOMSK_MSG(801,(/addatom_species/),(/0.d0/))
  nerr=nerr+1
  GOTO 1000
ENDIF
!
!
!
200 CONTINUE
SELECT CASE(addatom_type)
!
CASE("at","AT","@")
  !Simplest case: a new atom must be inserted at position x, y, z
  ALLOCATE( newP( SIZE(P,1)+1 , 4 ) )
  DO i=1,SIZE(P,1)
    newP(i,:) = P(i,:)
    IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
      IF( NINT(S(i,4))==NINT(snumber) ) THEN
        hasShells = .TRUE.
      ENDIF
    ENDIF
  ENDDO
  newP(SIZE(newP,1),1) = addatom_prop(1)
  newP(SIZE(newP,1),2) = addatom_prop(2)
  newP(SIZE(newP,1),3) = addatom_prop(3)
  newP(SIZE(newP,1),4) = snumber
  addedatoms = 1
  !
  !
CASE("relative","rel")
  !A new atom must be added near atom with the index addatom_prop(1)
  !Get index of atom
  atomindex = NINT(addatom_prop(1))
  !
  !Place atom at the given vector relatively to given atom
  x = P(atomindex,1) + addatom_prop(2)
  y = P(atomindex,2) + addatom_prop(3)
  z = P(atomindex,3) + addatom_prop(4)
  !
  !Save position of new atom in newP
  ALLOCATE( newP( SIZE(P,1)+1 , 4 ) )
  DO i=1,SIZE(P,1)
    newP(i,:) = P(i,:)
    IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
      IF( NINT(S(i,4))==NINT(snumber) ) THEN
        hasShells = .TRUE.
      ENDIF
    ENDIF
  ENDDO
  newP(SIZE(newP,1),1) = x
  newP(SIZE(newP,1),2) = y
  newP(SIZE(newP,1),3) = z
  newP(SIZE(newP,1),4) = snumber
  addedatoms = 1
  !
  !
CASE("near","NEAR")
  !A new atom must be added near atom with the index addatom_prop(1)
  !Get index of atom
  atomindex = NINT(addatom_prop(1))
  !
  !Perform a neighbor search around that atom
  CALL FIND_NNN(H,P,P(atomindex,1:3),8,V_NN,Nlist,exceeds100)
  !CALL FIND_1NN(H,P,P(atomindex,1:3),V_NN,Nlist,exceeds100)
  !
  !
  IF( SIZE(V_NN,1) >= 3 ) THEN
    !V_NN(:,:) now contains the positions of the nearest neighbors
    !Consider all possible tetrahedral sites that have P(atomindex,:) in a corner
    !Find distance to the closest neighbor to V_NN(1,:)
    dmax = 6.d0
    DO i=1,SIZE(V_NN,1)
      distance = VECLENGTH( P(atomindex,1:3) - V_NN(i,1:3) )
      IF( distance<dmax .AND. distance >1.d-12 ) THEN
        dmax = VECLENGTH( P(atomindex,1:3) - V_NN(i,1:3) )
      ENDIF
    ENDDO
    !
    !Compute the distance between the site center and P(atomindex,:)
    !Put new atom at the position (x,y,z) that maximizes its distance to P(atomindex,:)
    !(but is still closer than the first nearest neighbor)
    distance = 0.d0
    DO i=1,SIZE(V_NN,1)-2
      DO j=i+1,SIZE(V_NN,1)-1
        DO k=j+1,SIZE(V_NN,1)
          !Compute position of the center of this site
          V(1) = ( P(atomindex,1) + V_NN(i,1) + V_NN(j,1) + V_NN(k,1) ) / 4.d0
          V(2) = ( P(atomindex,2) + V_NN(i,2) + V_NN(j,2) + V_NN(k,2) ) / 4.d0
          V(3) = ( P(atomindex,3) + V_NN(i,3) + V_NN(j,3) + V_NN(k,3) ) / 4.d0
          !Compute distance between site center and P(atomindex,:)
          distance2 = VECLENGTH( P(atomindex,1:3) - V(:) )
          IF( distance2>distance .AND. distance2<dmax ) THEN
            !Save this position as a suitable candidate for the new atom
            x = V(1)
            y = V(2)
            z = V(3)
            distance = distance2
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    !
  ELSE
    !No suitable neighbors were found
    !This can only mean that P(atomindex,:) is far from any atom in the system
    !Just place the new atom near P(atomindex,:)
    !(this is completely arbitrary, but something has to be done)
    x = P(atomindex,1) + 0.5d0
    y = P(atomindex,2) + 0.5d0
    z = P(atomindex,3) + 0.5d0
  ENDIF
  IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
  !
  !Save position of new atom in newP
  ALLOCATE( newP( SIZE(P,1)+1 , 4 ) )
  DO i=1,SIZE(P,1)
    newP(i,:) = P(i,:)
    IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
      IF( NINT(S(i,4))==NINT(snumber) ) THEN
        hasShells = .TRUE.
      ENDIF
    ENDIF
  ENDDO
  newP(SIZE(newP,1),1) = x
  newP(SIZE(newP,1),2) = y
  newP(SIZE(newP,1),3) = z
  newP(SIZE(newP,1),4) = snumber
  addedatoms = 1
  !
  !
CASE("random","RANDOM","rand","RAND")
  !Insert NP atoms at random positions (but not too close to existing atom)
  NP = NINT(addatom_prop(1))
  !
  ALLOCATE(newP(SIZE(P,1)+NP,4))
  newP(:,:) = 0.d0
  DO i=1,SIZE(P,1)
    newP(i,:) = P(i,:)
  ENDDO
  !
  !Generate 3*NP random numbers
  CALL GEN_NRANDNUMBERS(3*NP,randarray)
  !
  !randarray(:) now contains 3*NP random numbers between 0 and 1
  !multiply them by the box dimensions to have cartesian coordinates
  DO n=1,NP
   randarray(n)      = H(1,1) * randarray(n)
   randarray(NP+n)   = H(2,2) * randarray(NP+n)
   randarray(2*NP+n) = H(3,3) * randarray(2*NP+n)
  ENDDO
  !
  !Construct neighbor list of new system with all atoms
  CALL ATOMSK_MSG(11,(/""/),(/0.d0/))
  CALL NEIGHBOR_LIST(H,newP(:,:),6.d0,NeighList)
  !
  !For each random position, search for the 4 nearest neighbors
  !and replace the position by the center of the 4 neighbors positions
  m = SIZE(P,1)
  DO n=1,NP
    !Gather the random coordinates of the atom to insert
    x = randarray(n)
    y = randarray(NP+n)
    z = randarray(2*NP+n)
    !
    !x,y,z are random and may be too close to an existing atom
    !To avoid that, replace x,y,z by the closest suitable tetrahedral site
    !Search for the nearest neighbors of the position (x,y,z)
    !Note: in order to avoid introducing two new atoms in the same site,
    !     perform the neighbor search in the system newP containing atoms previously introduced
    !CALL FIND_NNN(H,newP(1:m,:),(/x,y,z/),4,V_NN,Nlist,exceeds100)
    !
    !Generate list of positions of neighbors of atom #n
    CALL NEIGHBOR_POS(H,newP,(/x,y,z/),NeighList(m+n,:),6.d0,PosList)
    !
    IF( SIZE(PosList,1) >= 4 ) THEN
      !Atom #m+n has more than 4 neighbors => try to adjust its position
      !Sort neighbors by increasing distance
      CALL BUBBLESORT(PosList,4,'up  ')
      !Determine the equidistance of the 4 nearest atoms
      x = SUM( PosList(1:4,1) ) / 4.d0
      y = SUM( PosList(1:4,2) ) / 4.d0
      z = SUM( PosList(1:4,3) ) / 4.d0
    ENDIF
    !Save final position to newP
    m=m+1
    newP(m,1) = x
    newP(m,2) = y
    newP(m,3) = z
    newP(m,4) = snumber
    addedatoms = addedatoms+1
    !
    IF(ALLOCATED(PosList)) DEALLOCATE(PosList)
    IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
    !
    !Determine if this type of atoms has shells
    IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
      IF( NINT(S(i,4))==NINT(snumber) ) THEN
        hasShells = .TRUE.
      ENDIF
    ENDIF
  ENDDO
  !
  !
CASE DEFAULT
  !error
  GOTO 800
END SELECT
!
!
!Replace old P by newP; same with S and AUX if necessary
IF( addedatoms>0 ) THEN
  DEALLOCATE(P)
  ALLOCATE(P(SIZE(newP,1),4))
  P(:,:) = 0.d0
  IF(ALLOCATED(S)) THEN
    IF(ALLOCATED(newS)) DEALLOCATE(newS)
    ALLOCATE(newS(SIZE(newP,1),4))
    newS(:,:) = 0.d0
    DO i=1,SIZE(S,1)
      newS(i,:) = S(i,:)
    ENDDO
    DEALLOCATE(S)
    ALLOCATE(S(SIZE(P,1),4))
    S(:,:) = newS(:,:)
    DEALLOCATE(newS)
  ENDIF
  IF(ALLOCATED(AUX)) THEN
    DEALLOCATE(AUX)
    ALLOCATE( AUX( SIZE(newP,1),SIZE(newAUX,2) ) )
    AUX(:,:) = 0.d0
  ENDIF
  !
  DO i=1,SIZE(P,1)
    P(i,:) = newP(i,:)
    IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
      IF( i>SIZE(S,1)-addedatoms .AND. hasShells ) THEN
        !Other atoms of this type have shells
        !=> also add shells to new atoms
        S(i,1:4) = P(i,1:4)
      ENDIF
    ENDIF
    IF(ALLOCATED(AUX)) AUX(i,:) = newAUX(i,:)
  ENDDO
ENDIF
IF(ALLOCATED(newP)) DEALLOCATE(newP)
IF(ALLOCATED(newS)) DEALLOCATE(newS)
IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
!
!
!
300 CONTINUE
CALL ATOMSK_MSG(2117,(/''/),(/DBLE(addedatoms),DBLE(SIZE(P,1))/))
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
END SUBROUTINE ADDATOM_XYZ
!
END MODULE addatom