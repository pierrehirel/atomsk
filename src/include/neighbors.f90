MODULE neighbors
!
!**********************************************************************************
!*  ATOMS                                                                         *
!**********************************************************************************
!* This module contains subroutines concerning atoms.                             *
!**********************************************************************************
!* (C) Feb. 2014 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 20 June 2014                                     *
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
!* List of subroutines in this module:                                            *
!* NEIGHBOR_LIST       constructs a neighbor list                                 *
!* NEIGHBOR_POS        given a neighbor list, gives the coordinates of neighbors  *
!* FIND_NNN            finds the N nearest neighbors of a coordinate              *
!* FIND_NNRdR          finds the neighbors in a skin of radius R and width dR     *
!* FIND_1NN            finds the first nearest neighbors                          *
!**********************************************************************************
!
!
USE comv
USE functions
USE subroutines
!
!
CONTAINS
!
!
!********************************************************
! NEIGHBOR_LIST
! Considering a list of atom positions A, base
! vectors H defining periodic conditions, and a radius R,
! this routine builds a list of neighbors for all atoms,
! i.e. a list that looks like the following:
!        2 4 6 0
!        1 3 4 5
!        2 4 5 0
!        1 2 3 5
!    etc.
! meaning that neighbors of atom #1 are atoms #2, 4, 6,
! neighbors of atom #2 are #1, 3, 4, 5, etc.
! NOTE: atoms are assumed to be wrapped, i.e. all 
!      in the box (no atom outside of the box).
! NOTE2: the neighbor list NeighList(i,:) of the atom #i
!      may contain trailing zeros, as illustrated above.
!      These zeros are meant to be ignored. The first zero
!      encountered means there will be no more neighbors
! NOTE3: an atom is *never* counted as its own neighbor.
!      A neighbor appears only once in the neighbor list.
!      If the system is such that
!      replicas of atom #i are neighbors of atom #i,
!      or such that several replica of an atom #j can
!      be neighbors to an atom #i, then the replica do
!      not appear in the list, i.e. the atom #j is
!      counted only once. In such cases that must be
!      corrected *after* calling this routine.
!********************************************************
SUBROUTINE NEIGHBOR_LIST(H,A,R,NeighList)
!
IMPLICIT NONE
REAL(dp),INTENT(IN):: R  !radius in which Neighbors are searched
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),INTENT(IN):: A  !array of all atom positions
!
LOGICAL,DIMENSION(6):: IsCloseToBorder !is the atom close to the borders of the box?
INTEGER:: i, j, k, l, m, n, u
INTEGER:: kmin, kmax, lmin, lmax, mmin, mmax !Boundaries for neighbour search
INTEGER:: Nneighbours  !number of neighbors of atom i
INTEGER,PARAMETER:: NNincrement=2  !whenever list is full, increase its size by that much
INTEGER,DIMENSION(:,:),ALLOCATABLE:: tempList  !list of neighbours
REAL(dp):: distance
REAL(dp),DIMENSION(3):: d_border !atoms close to a border will be searched for periodic replica
REAL(dp),DIMENSION(1,3):: Vfrac  !position of an atom in reduced coordinates
!
INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: NeighList  !list of neighbours
!
!
!Initialize variables
IF(ALLOCATED(NeighList)) DEALLOCATE(NeighList)
ALLOCATE(NeighList(SIZE(A,1),4))  !initially, allow for 4 neighbors
NeighList(:,:) = 0
IF(ALLOCATED(tempList)) DEALLOCATE(tempList)
!
!Define "close to border" in reduced units
DO i=1,3
  !default=atoms closer than radius R from a border will be duplicated
  d_border(i) = MAX( 0.1d0 , R/VECLENGTH(H(i,:)) )
ENDDO
!
DO i=1,SIZE(A,1)
  !Count the neighbours that were already detected for atom i
  Nneighbours = 0
  n=1
  DO WHILE ( n<=SIZE(NeighList,2) .AND. NeighList(i,n).NE.0 )
    n = n+1
    Nneighbours = Nneighbours+1
  ENDDO
  !
  !Save fractional coordinate of atom i in Vfrac
  Vfrac(1,:) = A(i,1:3)
  CALL CART2FRAC(Vfrac,H)
  !
  !Check if atom i is near a border of the box,
  !"near" meaning "within the radius R"
  IsCloseToBorder(:) = .FALSE.
  IF( DABS(Vfrac(1,1))<d_border(1) ) THEN
    IsCloseToBorder(1) = .TRUE.
  ENDIF
  IF( DABS(1.d0-Vfrac(1,1))<d_border(1) ) THEN
    IsCloseToBorder(2) = .TRUE.
  ENDIF
  IF( DABS(Vfrac(1,2))<d_border(2) ) THEN
    IsCloseToBorder(3) = .TRUE.
  ENDIF
  IF( DABS(1.d0-Vfrac(1,2))<d_border(2) ) THEN
    IsCloseToBorder(4) = .TRUE.
  ENDIF
  IF( DABS(Vfrac(1,3))<d_border(3) ) THEN
    IsCloseToBorder(5) = .TRUE.
  ENDIF
  IF( DABS(1.d0-Vfrac(1,3))<d_border(3) ) THEN
    IsCloseToBorder(6) = .TRUE.
  ENDIF
  !
  DO j=i+1,SIZE(A,1)
    !Save fractional coordinate of atom j in Vfrac
    Vfrac(1,:) = A(j,1:3)
    CALL CART2FRAC(Vfrac,H)
    !
    !If atom i is close to a border, then
    !look for periodic replica of atom j across that border
    kmin = 0
    kmax = 0
    lmin = 0
    lmax = 0
    mmin = 0
    mmax = 0
    IF( IsCloseToBorder(1) ) THEN
      !Atom #i is close to the left border along x
      IF( DABS(1.d0-Vfrac(1,1))<d_border(1) ) THEN
        !Atom #j is close to the right border along x
        kmin = -1
      ENDIF
    ENDIF
    IF( IsCloseToBorder(2) ) THEN
      !Atom #i is close to the right border along x
      IF( DABS(Vfrac(1,1))<d_border(1) ) THEN
        !Atom #j is close to the left border along x
        kmax = 1
      ENDIF
    ENDIF
    IF( IsCloseToBorder(3) ) THEN
      !Atom #i is close to the bottom border along y
      IF( DABS(1.d0-Vfrac(1,2))<d_border(2) ) THEN
        !Atom #j is close to the top border along y
        lmin = -1
      ENDIF
    ENDIF
    IF( IsCloseToBorder(4) ) THEN
      !Atom #i is close to the top border along y
      IF( DABS(Vfrac(1,2))<d_border(2) ) THEN
        !Atom #j is close to the bottom border along y
        lmax = 1
      ENDIF
    ENDIF
    IF( IsCloseToBorder(5) ) THEN
      !Atom #i is close to the bottom border along z
      IF( DABS(1.d0-Vfrac(1,3))<d_border(3) ) THEN
        !Atom #j is close to the top border along z
        mmin = -1
      ENDIF
    ENDIF
    IF( IsCloseToBorder(6) ) THEN
      !Atom #i is close to the top border along z
      IF( DABS(Vfrac(1,3))<d_border(3) ) THEN
        !Atom #j is close to the bottom border along z
        mmax = 1
      ENDIF
    ENDIF
    !
    !Loop on periodic images
    DO k=kmin,kmax
      DO l=lmin,lmax
        DO m=mmin,mmax
          !Compute distance between atom i and this replica of atom j
          Vfrac(1,:) = A(j,1:3) + DBLE(k)*H(1,:) + DBLE(l)*H(2,:) + DBLE(m)*H(3,:)
          distance = VECLENGTH( A(i,1:3) - Vfrac(1,1:3) )
          !
          IF ( distance < R ) THEN
            !Atom j is neighbour of atom i
            Nneighbours = Nneighbours+1
            !If total number of neighbours exceeds size of NeighList, expand NeighList
            IF( Nneighbours > SIZE(NeighList,2) ) THEN
              !The neighbour list of this atom is full
              !=> Increase the size of NeighList by NNincrement
              IF( ALLOCATED(tempList) ) DEALLOCATE(tempList)
              ALLOCATE( tempList (SIZE(NeighList,1) , SIZE(NeighList,2)+NNincrement ) )
              tempList(:,:) = 0
              DO n=1,SIZE(NeighList,2)
                tempList(:,n) = NeighList(:,n)
              ENDDO
              tempList(i,SIZE(NeighList,2)+1) = j
              DEALLOCATE(NeighList)
              ALLOCATE( NeighList( SIZE(tempList,1) , SIZE(tempList,2) ) )
              NeighList(:,:) = tempList(:,:)
              DEALLOCATE(tempList)
              !
            ELSE
              !Add atom j to the list of neighbours of atom i
              NeighList(i,Nneighbours) = j
            ENDIF
            !
            !i is also a neighbour of atom j => save this
            n=1
            DO WHILE ( n<=SIZE(NeighList,2) .AND. NeighList(j,n).NE.0 )
              n=n+1
            ENDDO
            IF( n>SIZE(NeighList,2) ) THEN
              !The neighbour list of this atom is full
              !=> Increase the size of NeighList by NNincrement
              IF( ALLOCATED(tempList) ) DEALLOCATE(tempList)
              ALLOCATE( tempList (SIZE(NeighList,1) , SIZE(NeighList,2)+NNincrement ) )
              tempList(:,:) = 0
              DO u=1,SIZE(NeighList,2)
                tempList(:,u) = NeighList(:,u)
              ENDDO
              tempList(j,n) = i
              DEALLOCATE(NeighList)
              ALLOCATE( NeighList( SIZE(tempList,1) , SIZE(tempList,2) ) )
              NeighList(:,:) = tempList(:,:)
              DEALLOCATE(tempList)
            ELSE
              NeighList(j,n) = i
            ENDIF
            !
            !We already found that j is neighbour of i
            !=> no need to keep on looking for replica of atom j
            !=> exit the loops on replica
            GOTO 200
          ENDIF
          !
        ENDDO !m
      ENDDO  !l
    ENDDO   !k
    !
    200 CONTINUE
    !
  ENDDO !j
  !
ENDDO !i
!
!
IF( ALLOCATED(NeighList) ) THEN
  !If the neighbor list contains only zeros, then no neighbor was found
  !=> deallocate NeighList
  IF( .NOT.ANY(NeighList(:,:).NE.0) ) THEN
    DEALLOCATE(NeighList)
  ENDIF
ENDIF
!
!
END SUBROUTINE NEIGHBOR_LIST
!
!
!********************************************************
! NEIGHBOR_POS
! Given an atom position, the list of indices of its
! neighbors, and the box vectors H, this subroutine
! generates a list of the positions of the neighbors,
! taking periodic boundary conditions into account.
! The position V(:) of the central atom must be provided
! because it must not be counted as his own neighbor,
! but its periodic images may be neighbors.
! This routine returns an array PosList(:,:) containing,
! for each neighbor, its position x, y, z, its distance
! to the position V(:), and the index of the
! neighboring atom.
!********************************************************
SUBROUTINE NEIGHBOR_POS(H,A,V,NeighList,radius,PosList)
!
IMPLICIT NONE
LOGICAL:: selfneighbor  !is central atom a neighbor of itself? (because of PBC)
INTEGER:: i, m, n, o
INTEGER:: Nneighbors
INTEGER,DIMENSION(:),INTENT(IN):: NeighList !list of index of neighbors
REAL(dp):: distance
REAL(dp):: P1, P2, P3
REAL(dp),INTENT(IN):: radius
REAL(dp),DIMENSION(3):: vector
REAL(dp),DIMENSION(3),INTENT(IN):: V !position of central atom
REAL(dp),DIMENSION(3,3),INTENT(IN):: H  !vectors of the box
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: PosList
REAL(dp),DIMENSION(:,:),INTENT(IN):: A !positions of all atoms
!
selfneighbor=.FALSE.
IF(ALLOCATED(PosList)) DEALLOCATE(PosList)
IF(SIZE(NeighList)<=0) RETURN
!
Nneighbors=0
i=1
DO WHILE( i<=SIZE(NeighList) .AND. NeighList(i)>0 )
  !NeighList(i) is the index of a neighboring atom
  !Find which periodic image(s) are actually neighbors
  DO o=-1,1
    DO n=-1,1
      DO m=-1,1
        !Position of the periodic image of neighbor #i
        P1 = A(NeighList(i),1) + DBLE(m)*H(1,1) + DBLE(n)*H(2,1) + DBLE(o)*H(3,1)
        P2 = A(NeighList(i),2) + DBLE(m)*H(1,2) + DBLE(n)*H(2,2) + DBLE(o)*H(3,2)
        P3 = A(NeighList(i),3) + DBLE(m)*H(1,3) + DBLE(n)*H(2,3) + DBLE(o)*H(3,3)
        vector = (/ P1 , P2 , P3 /)
        !This image is a neighbor if distance is smaller than radius
        distance = VECLENGTH( vector(:) - V(:) )
        IF( distance>1.d-12 .AND. distance <= radius ) THEN
          Nneighbors = Nneighbors+1
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  i=i+1
ENDDO
!
DO o=-1,1
  DO n=-1,1
    DO m=-1,1
      !Position of the periodic image of central atom
      P1 = V(1) + DBLE(m)*H(1,1) + DBLE(n)*H(2,1) + DBLE(o)*H(3,1)
      P2 = V(2) + DBLE(m)*H(1,2) + DBLE(n)*H(2,2) + DBLE(o)*H(3,2)
      P3 = V(3) + DBLE(m)*H(1,3) + DBLE(n)*H(2,3) + DBLE(o)*H(3,3)
      vector = (/ P1 , P2 , P3 /)
      !This image is a neighbor if distance is smaller than radius
      distance = VECLENGTH( vector(:) - V(:) )
      IF( distance>1.d-12 .AND. distance <= radius ) THEN
        Nneighbors = Nneighbors+1
        selfneighbor=.TRUE.
      ENDIF
    ENDDO
  ENDDO
ENDDO
!
IF( Nneighbors>0 ) THEN
  ALLOCATE(PosList(Nneighbors,5))
  PosList(:,:) = 0.d0
  !
  Nneighbors=0
  i=1
  DO WHILE( i<=SIZE(NeighList) .AND. NeighList(i)>0 )
    DO o=-1,1
      DO n=-1,1
        DO m=-1,1
          !Position of the periodic image of neighbor #i
          P1 = A(NeighList(i),1) + DBLE(m)*H(1,1) + DBLE(n)*H(2,1) + DBLE(o)*H(3,1)
          P2 = A(NeighList(i),2) + DBLE(m)*H(1,2) + DBLE(n)*H(2,2) + DBLE(o)*H(3,2)
          P3 = A(NeighList(i),3) + DBLE(m)*H(1,3) + DBLE(n)*H(2,3) + DBLE(o)*H(3,3)
          vector = (/ P1 , P2 , P3 /)
          !This image is a neighbor if distance is smaller than radius
          distance = VECLENGTH( vector(:) - V(:) )
          IF( distance>1.d-12 .AND. distance <= radius ) THEN
            Nneighbors = Nneighbors+1
            PosList(Nneighbors,1:3) = vector(:)
            PosList(Nneighbors,4) = distance
            PosList(Nneighbors,5) = DBLE(NeighList(i))
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    i=i+1
  ENDDO
  !
  IF( selfneighbor ) THEN
    DO o=-1,1
      DO n=-1,1
        DO m=-1,1
          !Position of the periodic image of central atom
          P1 = V(1) + DBLE(m)*H(1,1) + DBLE(n)*H(2,1) + DBLE(o)*H(3,1)
          P2 = V(2) + DBLE(m)*H(1,2) + DBLE(n)*H(2,2) + DBLE(o)*H(3,2)
          P3 = V(3) + DBLE(m)*H(1,3) + DBLE(n)*H(2,3) + DBLE(o)*H(3,3)
          vector = (/ P1 , P2 , P3 /)
          !This image is a neighbor if distance is smaller than radius
          distance = VECLENGTH( vector(:) - V(:) )
          IF( distance>1.d-12 .AND. distance <= radius ) THEN
            Nneighbors = Nneighbors+1
            PosList(Nneighbors,1:3) = vector(:)
            PosList(Nneighbors,4) = distance
            PosList(Nneighbors,5) = 0.d0
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  !
ENDIF
!
!
END SUBROUTINE NEIGHBOR_POS
!
!
!********************************************************
! FIND_NNN
! This subroutine reads a real array A of dimension QxM
! with M>=3, assumed to contain atom positions in the
! first three columns, and finds the NNN (N nearest
! neighbours) to a position V (array of size M>=3).
! It returns their coordinates in an array V_NN,
! from the closest to V up to the farest. If all the
! nearest neighbours are at the same distance from
! point V, then they will appear in V_NN in their
! order of appearence in A.
! It also returns the neighbours list in Nlist, which
! contains the indices of neighbours in the same order
! as they appear in V_NN.
! Note that atom positions in A and V must be given as
! reduced coordinates.
! This subroutine also takes periodic boundary conditions
! into account (using the supercell vectors H).
!********************************************************
SUBROUTINE FIND_NNN(H,A,V,NNN,V_NN,Nlist,exceeds100)
!
IMPLICIT NONE
INTEGER:: i, j, k, l, m, n
INTEGER:: kmin, kmax, lmin, lmax, mmin, mmax !Boundaries for neighbour search
INTEGER,INTENT(IN):: NNN     !Number of Nearest Neighbours to find
INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(OUT):: Nlist !index of neighbours
LOGICAL:: diff_from_prev_NN !is current atom different from neighbours previously found?
LOGICAL:: exceeds100 !does the number of neighbours exceed 100?
LOGICAL,DIMENSION(6):: VIsCloseToBorder  !is the position V close to the borders?
REAL(dp):: distance, distance2 !distance between current atom and V
REAL(dp):: tempreal
REAL(dp),DIMENSION(3):: d_border !atoms close to a border will be searched for periodic replica
REAL(dp),DIMENSION(:),INTENT(IN):: V  !position around which NN must be searched
REAL(dp),DIMENSION(1,3):: Vfrac  !copy of V in reduced coordinates
REAL(dp),DIMENSION(:),ALLOCATABLE:: currA !position of current atom
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),INTENT(IN):: A  !array of all atom positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Afrac  !copy of A in reduced coordinates
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: V_100 !positions of first 100 neighbours
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: V_NN  !final positions of 1st NN
!
!Initialize variables
exceeds100 = .FALSE.
DO i=1,3
  !default=atoms closer than 0.1 cell vector or 10A from a border
  !        will be duplicated
  d_border(i) = MAX(0.1d0,10.d0/H(i,i))
ENDDO
IF(ALLOCATED(currA)) DEALLOCATE(currA)
 ALLOCATE( currA( SIZE(A(1,:)) ) )
 currA(:) = 0.d0
IF(ALLOCATED(V_100)) DEALLOCATE(V_100)
 ALLOCATE( V_100( 100,SIZE(A(1,:)) ) )
 V_100(:,:) = 0.d0
IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
IF(ALLOCATED(Nlist)) DEALLOCATE(Nlist)
IF(ALLOCATED(Afrac)) DEALLOCATE(Afrac)
ALLOCATE( Afrac( SIZE(A,1),SIZE(A,2) ) )
Afrac = A
CALL CART2FRAC(Afrac,H)
Vfrac(1,:) = V(1:3)
CALL CART2FRAC(Vfrac,H)
!
!Check if the position V is close to a border
!Note: we have to use reduced coordinates here in case
!the box is not cubic
VIsCloseToBorder(:) = .FALSE.
IF( DABS(Vfrac(1,1))<d_border(1) ) THEN
  VIsCloseToBorder(1) = .TRUE.
ENDIF
IF( DABS(1.d0-Vfrac(1,1))<d_border(1) ) THEN
  VIsCloseToBorder(2) = .TRUE.
ENDIF
IF( DABS(Vfrac(1,2))<d_border(2) ) THEN
  VIsCloseToBorder(3) = .TRUE.
ENDIF
IF( DABS(1.d0-Vfrac(1,2))<d_border(2) ) THEN
  VIsCloseToBorder(4) = .TRUE.
ENDIF
IF( DABS(Vfrac(1,3))<d_border(3) ) THEN
  VIsCloseToBorder(5) = .TRUE.
ENDIF
IF( DABS(1.d0-Vfrac(1,3))<d_border(3) ) THEN
  VIsCloseToBorder(6) = .TRUE.
ENDIF
!
!Start the neighbours search
!
IF(NNN>100.d0) THEN
  !If NNN is greater than 100 we are in trouble
  exceeds100 = .TRUE.
!
ELSE
  !If NNN is positive then it corresponds to the number of neighbours
  !the program has to find
  ALLOCATE( V_NN( NNN,SIZE(A(1,:)) ) )
  V_NN(:,:) = 0.d0
  ALLOCATE(Nlist(NNN))
  Nlist(:) = 0
  !
  i=0
  !Loop on the N nearest neighbours to find
  DO WHILE( i<MAX(1,NNN) )  !we always find at least 1 neighbour
    !
    i=i+1
    tempreal = 1.d9
    !
    DO j=1,SIZE(A(:,1))  !Loop on all atoms
      !
      !We need to look for periodic images only for atoms close to a boundary
      !The purpose is to minimize the ranges [min,max] in which periodic
      !images will be generated in the big embedded loops below
      kmin = 0
      kmax = 0
      lmin = 0
      lmax = 0
      mmin = 0
      mmax = 0
      IF( VIsCloseToBorder(1) ) THEN
        IF( DABS(1.d0-Afrac(j,1))<d_border(1) ) THEN
          kmin = -1
        ENDIF
      ENDIF
      IF( VIsCloseToBorder(2) ) THEN
        IF( DABS(Afrac(j,1))<d_border(1) ) THEN
          kmax = 1
        ENDIF
      ENDIF
      IF( VIsCloseToBorder(3) ) THEN
        IF( DABS(1.d0-Afrac(j,2))<d_border(2) ) THEN
          lmin = -1
        ENDIF
      ENDIF
      IF( VIsCloseToBorder(4) ) THEN
        IF( DABS(Afrac(j,2))<d_border(2) ) THEN
          lmax = 1
        ENDIF
      ENDIF
      IF( VIsCloseToBorder(5) ) THEN
        IF( DABS(1.d0-Afrac(j,3))<d_border(3) ) THEN
          mmin = -1
        ENDIF
      ENDIF
      IF( VIsCloseToBorder(6) ) THEN
        IF( DABS(Afrac(j,3))<d_border(3) ) THEN
          mmax = 1
        ENDIF
      ENDIF
      !
      !Loop on periodic images
      DO k=kmin,kmax
        DO l=lmin,lmax
          DO m=mmin,mmax
            diff_from_prev_NN = .TRUE.
            !Position of the image of current atom
            currA(1:3) = A(j,1:3) + DBLE(k)*H(1,:) + DBLE(l)*H(2,:) + DBLE(m)*H(3,:)
            !Compute distance between current atom and V
            distance = VECLENGTH( currA(1:3)-V(1:3) )
            !
            !Compute minimum distance from neighbours previously found
            DO n=1,i
              distance2 = VECLENGTH( currA(1:3)-V_100(n,1:3) )
              !If distance 2 is very small then we already found
              !that neighbour, skip it
              IF( DABS(distance2)<=1.d-6 ) diff_from_prev_NN = .FALSE.
            ENDDO
            !
            !If distance is smaller than before we save it to V_100,
            !but only if it is also different from the previous neighbour found
            !If distance is close to zero, ignore it: it means that we found the atom
            !that is exactly at position V(:), we don't want to count itself as its own neighbor
            IF( distance<tempreal .AND. diff_from_prev_NN .AND. distance>1.d-12 ) THEN
              tempreal = distance
              !Copy contents of currA to V_100 (coordinates)
              V_100(i,:) = currA(:)
              !Save index of atom into Nlist
              Nlist(i) = j
            ENDIF
            !
          ENDDO !m
        ENDDO !l
      ENDDO !k
    ENDDO !j
    !
  ENDDO !i
  !
  !Fill V_NN with the coordinates of neighbours
  DO j=1,NNN
    V_NN(j,:) = V_100(j,:)
  ENDDO
  !
ENDIF
!
DEALLOCATE(Afrac)
!
END SUBROUTINE FIND_NNN
!
!
!
!********************************************************
! FIND_NNRdR
! This subroutine reads a real array A of dimension QxM
! with M>=3, assumed to contain atom positions in the
! first three columns, and finds all neighbours within a
! distance R < d <= R+dR around a position V
! (array of size M>=3).
! In order to find atoms in a sphere of a certain radius,
! this subroutine can be called with R<=0 and
! dR=radius of the sphere.
! Coordinates of neighbours are returned in an array V_NN,
! from the closest to V up to the farest. If all the
! nearest neighbours are at the same distance from
! point V, then they will appear in V_NN
! in their order of appearence in A.
! The coordinates in A must be given as cartesian
! coordinates.
! This subroutine also takes periodic boundary
! conditions into account (using the supercell vectors H).
!********************************************************
SUBROUTINE FIND_NNRdR(H,A,V,R,dR,V_NN,Nlist,exceeds100)
!
IMPLICIT NONE
INTEGER:: i, j, k, l, m
INTEGER:: kmin, kmax, lmin, lmax, mmin, mmax !Boundaries for neighbour search
INTEGER:: Vsize  !Size of array for storing positions of nearest neighbours
INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(OUT):: Nlist !index of neighbours
INTEGER,DIMENSION(:),ALLOCATABLE:: Nlist100 !index of first 100 neighbours
LOGICAL:: exceeds100 !does the number of neighbours exceed 100?
LOGICAL,DIMENSION(6):: VIsCloseToBorder  !is the position V close to the borders?
REAL(dp):: distance !distance between current atom and V
REAL(dp),INTENT(IN):: R, dR  !radius in which Neighbour are searched
REAL(dp),DIMENSION(3):: d_border !atoms close to a border will be searched for periodic replica
REAL(dp),DIMENSION(:),INTENT(IN):: V  !position around which NN must be searched
REAL(dp),DIMENSION(1,3):: Vfrac  !copy of V in reduced coordinates
REAL(dp),DIMENSION(:),ALLOCATABLE:: currA !position of current atom
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),INTENT(IN):: A  !array of all atom positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Afrac  !copy of A in reduced coordinates
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: V_100 !positions of first 100 neighbours
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: V_NN  !final positions of 1st NN
!
!Initialize variables
exceeds100 = .FALSE.
DO i=1,3
  !default=atoms closer than R+dR from a border
  !        will be duplicated
  d_border(i) = (R+dR)/H(i,i)
ENDDO
IF(ALLOCATED(currA)) DEALLOCATE(currA)
 ALLOCATE( currA( SIZE(A(1,:)) ) )
 currA(:) = 0.d0
IF(ALLOCATED(V_100)) DEALLOCATE(V_100)
 ALLOCATE( V_100( 1000,SIZE(A(1,:)) ) )
 V_100(:,:) = 0.d0
IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
IF(ALLOCATED(Nlist)) DEALLOCATE(Nlist)
IF(ALLOCATED(Nlist100)) DEALLOCATE(Nlist100)
 ALLOCATE( Nlist100(1000) )
 Nlist100(:) = 0
IF(ALLOCATED(Afrac)) DEALLOCATE(Afrac)
ALLOCATE( Afrac( SIZE(A,1),SIZE(A,2) ) )
Afrac = A
CALL CART2FRAC(Afrac,H)
Vfrac(1,:) = V(1:3)
CALL CART2FRAC(Vfrac,H)
!
!Check if the position V is close to a border
!Note: we have to use reduced coordinates here in case
!the box is not cubic
VIsCloseToBorder(:) = .FALSE.
IF( DABS(Vfrac(1,1))<d_border(1) ) THEN
  VIsCloseToBorder(1) = .TRUE.
ENDIF
IF( DABS(1.d0-Vfrac(1,1))<d_border(1) ) THEN
  VIsCloseToBorder(2) = .TRUE.
ENDIF
IF( DABS(Vfrac(1,2))<d_border(2) ) THEN
  VIsCloseToBorder(3) = .TRUE.
ENDIF
IF( DABS(1.d0-Vfrac(1,2))<d_border(2) ) THEN
  VIsCloseToBorder(4) = .TRUE.
ENDIF
IF( DABS(Vfrac(1,3))<d_border(3) ) THEN
  VIsCloseToBorder(5) = .TRUE.
ENDIF
IF( DABS(1.d0-Vfrac(1,3))<d_border(3) ) THEN
  VIsCloseToBorder(6) = .TRUE.
ENDIF
!
!First we have to count the neighbours
Vsize=0
!
DO j=1,SIZE(A(:,1))  !Loop on all atoms (incl. periodic images)
  !We need to look for periodic images only for atoms close to a boundary
  !The purpose is to minimize the ranges [min,max] in which periodic
  !images will be generated in the big embedded loops below
  kmin = 0
  kmax = 0
  lmin = 0
  lmax = 0
  mmin = 0
  mmax = 0
  IF( VIsCloseToBorder(1) ) THEN
    IF( DABS(1.d0-Afrac(j,1))<d_border(1) ) THEN
      kmin = -1
    ENDIF
  ENDIF
  IF( VIsCloseToBorder(2) ) THEN
    IF( DABS(Afrac(j,1))<d_border(1) ) THEN
      kmax = 1
    ENDIF
  ENDIF
  IF( VIsCloseToBorder(3) ) THEN
    IF( DABS(1.d0-Afrac(j,2))<d_border(2) ) THEN
      lmin = -1
    ENDIF
  ENDIF
  IF( VIsCloseToBorder(4) ) THEN
    IF( DABS(Afrac(j,2))<d_border(2) ) THEN
      lmax = 1
    ENDIF
  ENDIF
  IF( VIsCloseToBorder(5) ) THEN
    IF( DABS(1.d0-Afrac(j,3))<d_border(3) ) THEN
      mmin = -1
    ENDIF
  ENDIF
  IF( VIsCloseToBorder(6) ) THEN
    IF( DABS(Afrac(j,3))<d_border(3) ) THEN
      mmax = 1
    ENDIF
  ENDIF
  !
  !Loop on periodic images
  DO k=kmin,kmax  !Loop on periodic images
    DO l=lmin,lmax
      DO m=mmin,mmax
        !Position of the image of current atom
        currA(1:3) = A(j,1:3) + DBLE(k)*H(1,:) + DBLE(l)*H(2,:) + DBLE(m)*H(3,:)
        !Compute distance between current atom and V
        distance = VECLENGTH( currA(1:3)-V(1:3) )
        !If distance is in the skin we consider it a neighbour
        IF( distance>R .AND. distance<=R+dR ) THEN
          Vsize=Vsize+1
          IF(Vsize>1000) THEN
            exceeds100 = .TRUE.
            EXIT
          ELSE
            V_100(Vsize,:) = currA(:)
            Nlist100(Vsize) = j
          ENDIF
        ENDIF
      ENDDO !m
    ENDDO !l
  ENDDO !k
ENDDO !j
!
!Now we know the number of neighbours, allocate V_NN
ALLOCATE( V_NN( Vsize,SIZE(A(1,:)) ) )
V_NN(:,:) = 0.d0
ALLOCATE( Nlist(Vsize) )
Nlist(:) = 0
!
!Fill V_NN with the coordinates of neighbours
DO j=1,SIZE(V_NN(:,1))
  V_NN(j,:) = V_100(j,:)
  Nlist(j) = Nlist100(j)
ENDDO
!
DEALLOCATE(Afrac)
!
END SUBROUTINE FIND_NNRdR
!
!
!********************************************************
! FIND_1NN
! This subroutine reads a real array A of dimension QxM
! with M>=3, assumed to contain atom positions in the
! first three columns, and finds the first nearest
! neighbours to a position V (array of size M>=3).
!********************************************************
SUBROUTINE FIND_1NN(H,A,V,V_NN,Nlist,exceeds100)
!
IMPLICIT NONE
INTEGER:: i, j, k, l, m
INTEGER:: kmin, kmax, lmin, lmax, mmin, mmax !Boundaries for neighbour search
INTEGER:: Vsize  !Size of array for storing positions of nearest neighbours
INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(OUT):: Nlist !index of neighbours
INTEGER,DIMENSION(:),ALLOCATABLE:: Nlist100 !index of first 100 neighbours
LOGICAL:: exceeds100 !does the number of neighbours exceed 100?
LOGICAL,DIMENSION(6):: VIsCloseToBorder  !is the position V close to the borders?
REAL(dp):: distance, distance2 !distance between current atom and V
REAL(dp):: tempreal
REAL(dp),DIMENSION(3):: d_border !atoms close to a border will be searched for periodic replica
REAL(dp),DIMENSION(:),INTENT(IN):: V  !position around which NN must be searched
REAL(dp),DIMENSION(1,3):: Vfrac  !copy of V in reduced coordinates
REAL(dp),DIMENSION(:),ALLOCATABLE:: currA !position of current atom
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),INTENT(IN):: A  !array of all atom positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Afrac  !copy of A in reduced coordinates
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: V_100 !positions of first 100 neighbours
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: V_NN  !final positions of 1st NN
!
!Initialize variables
exceeds100 = .FALSE.
DO i=1,3
  !default=atoms closer than R+dR from a border
  !        will be duplicated
  d_border(i) = MAX(0.1d0,10.d0/H(i,i))
ENDDO
IF(ALLOCATED(currA)) DEALLOCATE(currA)
 ALLOCATE( currA( SIZE(A(1,:)) ) )
 currA(:) = 0.d0
IF(ALLOCATED(V_100)) DEALLOCATE(V_100)
 ALLOCATE( V_100( 100,SIZE(A(1,:)) ) )
 V_100(:,:) = 0.d0
IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
IF(ALLOCATED(Nlist100)) DEALLOCATE(Nlist100)
 ALLOCATE( Nlist100(100) )
 Nlist100(:) = 0
IF(ALLOCATED(Nlist)) DEALLOCATE(Nlist)
IF(ALLOCATED(Afrac)) DEALLOCATE(Afrac)
ALLOCATE( Afrac( SIZE(A,1),SIZE(A,2) ) )
Afrac = A
CALL CART2FRAC(Afrac,H)
Vfrac(1,:) = V(1:3)
CALL CART2FRAC(Vfrac,H)
!
!
!Check if the position V is close to a border
!Note: we have to use reduced coordinates here in case
!the box is not cubic
VIsCloseToBorder(:) = .FALSE.
IF( DABS(Vfrac(1,1))<d_border(1) ) THEN
  VIsCloseToBorder(1) = .TRUE.
ENDIF
IF( DABS(1.d0-Vfrac(1,1))<d_border(1) ) THEN
  VIsCloseToBorder(2) = .TRUE.
ENDIF
IF( DABS(Vfrac(1,2))<d_border(2) ) THEN
  VIsCloseToBorder(3) = .TRUE.
ENDIF
IF( DABS(1.d0-Vfrac(1,2))<d_border(2) ) THEN
  VIsCloseToBorder(4) = .TRUE.
ENDIF
IF( DABS(Vfrac(1,3))<d_border(3) ) THEN
  VIsCloseToBorder(5) = .TRUE.
ENDIF
IF( DABS(1.d0-Vfrac(1,3))<d_border(3) ) THEN
  VIsCloseToBorder(6) = .TRUE.
ENDIF
!
!First, find the closest neighbour
tempreal = 1.d9
!
!Loop on all atoms
DO j=1,SIZE(A(:,1))
  !We need to look for periodic images only for atoms close to a boundary
  !The purpose is to minimize the ranges [min,max] in which periodic
  !images will be generated in the big embedded loops below
  kmin = 0
  kmax = 0
  lmin = 0
  lmax = 0
  mmin = 0
  mmax = 0
  IF( VIsCloseToBorder(1) ) THEN
    IF( DABS(1.d0-Afrac(j,1))<d_border(1) ) THEN
      kmin = -1
    ENDIF
  ENDIF
  IF( VIsCloseToBorder(2) ) THEN
    IF( DABS(Afrac(j,1))<d_border(1) ) THEN
      kmax = 1
    ENDIF
  ENDIF
  IF( VIsCloseToBorder(3) ) THEN
    IF( DABS(1.d0-Afrac(j,2))<d_border(2) ) THEN
      lmin = -1
    ENDIF
  ENDIF
  IF( VIsCloseToBorder(4) ) THEN
    IF( DABS(Afrac(j,2))<d_border(2) ) THEN
      lmax = 1
    ENDIF
  ENDIF
  IF( VIsCloseToBorder(5) ) THEN
    IF( DABS(1.d0-Afrac(j,3))<d_border(3) ) THEN
      mmin = -1
    ENDIF
  ENDIF
  IF( VIsCloseToBorder(6) ) THEN
    IF( DABS(Afrac(j,3))<d_border(3) ) THEN
      mmax = 1
    ENDIF
  ENDIF
  !
  !Loop on periodic images
  DO k=kmin,kmax
    DO l=lmin,lmax
      DO m=mmin,mmax
        !Position of the image of current atom
        currA(1:3) = A(j,1:3) + DBLE(k)*H(1,:) + DBLE(l)*H(2,:) + DBLE(m)*H(3,:)
        !Compute distance between current atom and V
        distance = VECLENGTH( currA(1:3)-V(1:3) )
        !
        !If distance is smaller than before we save it to V_100(1,:)
        IF( distance<tempreal .AND. distance>1.d-6 ) THEN
          tempreal = distance
          V_100(1,:) = currA(:)
          Nlist100(1) = j
        ENDIF
      ENDDO !m
    ENDDO !l
  ENDDO !k
ENDDO !j
!
distance = VECLENGTH(V_100(1,1:3)-V(1:3))
!
!Now find atoms that are at the same distance +/- epsilon
Vsize=1
DO j=1,SIZE(A(:,1))
  kmin = 0
  kmax = 0
  lmin = 0
  lmax = 0
  mmin = 0
  mmax = 0
  IF( VIsCloseToBorder(1) ) THEN
    IF( DABS(1.d0-Afrac(j,1))<d_border(1) ) THEN
      kmin = -1
    ENDIF
  ENDIF
  IF( VIsCloseToBorder(2) ) THEN
    IF( DABS(Afrac(j,1))<d_border(1) ) THEN
      kmax = 1
    ENDIF
  ENDIF
  IF( VIsCloseToBorder(3) ) THEN
    IF( DABS(1.d0-Afrac(j,2))<d_border(2) ) THEN
      lmin = -1
    ENDIF
  ENDIF
  IF( VIsCloseToBorder(4) ) THEN
    IF( DABS(Afrac(j,2))<d_border(2) ) THEN
      lmax = 1
    ENDIF
  ENDIF
  IF( VIsCloseToBorder(5) ) THEN
    IF( DABS(1.d0-Afrac(j,3))<d_border(3) ) THEN
      mmin = -1
    ENDIF
  ENDIF
  IF( VIsCloseToBorder(6) ) THEN
    IF( DABS(Afrac(j,3))<d_border(3) ) THEN
      mmax = 1
    ENDIF
  ENDIF
  !
  !Loop on periodic images
  DO k=kmin,kmax
    DO l=lmin,lmax
      DO m=mmin,mmax
        !Position of the image of current atom
        currA(1:3) = A(j,1:3) + k*H(1,:) + l*H(2,:) + m*H(3,:)
        !Compute distance between current atom and V
        distance2 = VECLENGTH(currA(1:3)-V(1:3))
        !
        !If "distance2" is approximately the same as "distance", save it to V_100
        IF( distance2<distance*1.1d0 .AND. VECLENGTH(currA(1:3)-V_100(1,1:3))>1.d-6 ) THEN
          Vsize=Vsize+1
          IF(Vsize>100) THEN
            exceeds100 = .TRUE.
            EXIT
          ELSE
            V_100(Vsize,:) = currA(:)
            Nlist100(Vsize) = j
          ENDIF
        ENDIF
      ENDDO !m
    ENDDO !l
  ENDDO !k
ENDDO !j
!
!Now we know the number of neighbours, allocate V_NN
ALLOCATE( V_NN( Vsize,SIZE(A(1,:)) ) )
V_NN(:,:) = 0.d0
ALLOCATE( Nlist(Vsize) )
Nlist(:) = 0
!
!Fill V_NN with the coordinates of neighbours
DO j=1,Vsize
  V_NN(j,:) = V_100(j,:)
  Nlist(j) = Nlist100(j)
ENDDO
!
IF(ALLOCATED(currA)) DEALLOCATE(currA)
DEALLOCATE(Afrac)
!
!
END SUBROUTINE FIND_1NN
!
!
!
END MODULE neighbors