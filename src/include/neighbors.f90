MODULE neighbors
!
!**********************************************************************************
!*  NEIGHBORS                                                                     *
!**********************************************************************************
!* This module contains subroutines performing neighbor search.                   *
!**********************************************************************************
!* (C) Feb. 2014 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 12 Nov. 2019                                     *
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
USE messages
USE resize
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
!        2 4 6 0 0 0 0
!        1 3 4 5 0 0 0
!        2 4 5 0 0 0 0
!        1 2 3 5 6 0 0
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
! NOTE4: if the number of atoms is small, then a
!      simplistic Verlet search is performed, which
!      scales as N². Otherwise a cell-list
!      algorithm is used, which scales as N.
! NOTE5: if no neighbor is found, then the array
!      NeighList is returned as UNALLOCATED.
!********************************************************
SUBROUTINE NEIGHBOR_LIST(H,A,R,NeighList)
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
REAL(dp),INTENT(IN):: R  !radius in which Neighbors are searched
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),INTENT(IN):: A  !array of all atom positions
!
LOGICAL,DIMENSION(6):: IsCloseToBorder !is the atom close to the borders of the box?
INTEGER:: a1, a2, a3
INTEGER:: Ix, Iy, Iz
INTEGER:: i, j, k, l, m, n, u
INTEGER:: iCell  !index of a cell (cell-list algorithm)
INTEGER:: kmin, kmax, lmin, lmax, mmin, mmax !Boundaries for neighbour search
INTEGER:: Maxcells !max number of cells along one direction (cell-list algorithm)
INTEGER:: Ncells   !total number of cells (cell-list algorithm)
INTEGER:: OMP_GET_NUM_THREADS
INTEGER,DIMENSION(3):: shift
INTEGER,PARAMETER:: NNincrement=2  !whenever list is full, increase its size by that much
INTEGER,DIMENSION(3):: NcellsX  !number of cells along X, Y, Z (cell-list algorithm)
INTEGER,DIMENSION(:),ALLOCATABLE:: Cell_NP       !number of atoms in each cell (cell-list algorithm)
INTEGER,DIMENSION(:,:),ALLOCATABLE:: Cell_AtomID !index of atoms for each cell (cell-list algorithm)
INTEGER,DIMENSION(:,:),ALLOCATABLE:: tempList    !list of neighbours (temporary)
INTEGER,DIMENSION(:,:,:),ALLOCATABLE:: Cell_Neigh  !neighbors of each cell (cell-list algorithm)
REAL(dp):: distance
REAL(dp):: tempreal
REAL(dp),DIMENSION(3):: d_border !atoms close to a border will be searched for periodic replica
REAL(dp),DIMENSION(3):: Cell_L   !length of cell along X, Y, Z (cell-list algorithm)
INTEGER,DIMENSION(:),ALLOCATABLE:: NNeigh      !number of neighbors of atom #i
INTEGER,DIMENSION(:),ALLOCATABLE:: Atom_Cell   !for each atom, index of the cell it belongs to
REAL(dp),DIMENSION(1,3):: Vfrac  !position of an atom in reduced coordinates
!
INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: NeighList  !the neighbor list
!
!
!Initialize variables
NcellsX(:) = 3
IF(ALLOCATED(NNeigh)) DEALLOCATE(NNeigh)
ALLOCATE(NNeigh(SIZE(A,1)) , STAT=i)
IF( i>0 ) THEN
  ! Allocation failed (not enough memory)
  nerr = nerr+1
  CALL ATOMSK_MSG(819,(/''/),(/-1.d0/))
  RETURN
ENDIF
NNeigh(:) = 0
IF(ALLOCATED(NeighList)) DEALLOCATE(NeighList)
!Assume a continuous atom density to estimate initial size of NeighList
!NOTE: if an atom has a greater number of neighbors the size of NeighList will be changed later
CALL VOLUME_PARA(H,distance)
n = MAX( 100 , NINT((DBLE(SIZE(A,1))/distance)*(2.d0/3.d0)*pi*(R**3)) )
ALLOCATE(NeighList(SIZE(A,1),n) , STAT=i)
IF( i>0 ) THEN
  ! Allocation failed (not enough memory)
  nerr = nerr+1
  CALL ATOMSK_MSG(819,(/''/),(/-1.d0/))
  RETURN
ENDIF
NeighList(:,:) = 0
IF(ALLOCATED(tempList)) DEALLOCATE(tempList)
IF(ALLOCATED(Cell_NP)) DEALLOCATE(Cell_NP)
IF(ALLOCATED(Cell_AtomID)) DEALLOCATE(Cell_AtomID)
IF(ALLOCATED(Cell_Neigh)) DEALLOCATE(Cell_Neigh)
!
!
msg = 'entering NEIGHBOR_LIST'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,*) 'Radius for neighbor search (angstroms) = ', R
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
IF( (VECLENGTH(H(1,:))<1.2d0*R .OR. VECLENGTH(H(2,:))<1.2d0*R .OR. VECLENGTH(H(3,:))<1.2d0*R) &
  & .AND. SIZE(A,1) < 2000 ) THEN
  !
  !System is pseudo-2D or contains a small number of atoms
  !=> a simplistic Verlet neighbor search will suffice
  msg = 'algorithm: VERLET'
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !
  !Define "close to border" in reduced units
  DO i=1,3
    !default=atoms closer than radius R from a border will be duplicated
    d_border(i) = MAX( 0.1d0 , R/VECLENGTH(H(i,:)) )
  ENDDO
  !
  DO i=1,SIZE(A,1)-1
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
              !Atom j is neighbor of atom i
              NNeigh(i) = NNeigh(i)+1
              !If total number of neighbors exceeds size of NeighList, expand NeighList
              IF( NNeigh(i) > SIZE(NeighList,2) ) THEN
                !The neighbor list of this atom is full
                !=> Increase the size of NeighList by NNincrement
                CALL RESIZE_INTARRAY2(NeighList,SIZE(NeighList,1),SIZE(NeighList,2)+NNincrement,l)
              ENDIF
              !Add atom j to the list of neighbors of atom i
              NeighList(i,NNeigh(i)) = j
              !
              !Atom i is also a neighbor of atom j => save this
              Nneigh(j) = NNeigh(j)+1
              IF( NNeigh(j) > SIZE(NeighList,2) ) THEN
                !The neighbor list of this atom is full
                !=> Increase the size of NeighList by NNincrement
                CALL RESIZE_INTARRAY2(NeighList,SIZE(NeighList,1),SIZE(NeighList,2)+NNincrement,l)
              ENDIF
              !Add atom i to the list of neighbors of atom j
              NeighList(j,NNeigh(j)) = i
              !
              !We already found that j is neighbor of i
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
  IF(ALLOCATED(NNeigh)) DEALLOCATE(NNeigh)
  !
  !
ELSE
  !Large system => use a cell list algorithm
  msg = 'algorithm: CELL DECOMPOSITION'
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !Define max. number of cells along each direction
#if defined(OPENMP)
  !Parallel version: use Nthreads to define max. number of cells along each direction
  !$OMP PARALLEL
  Maxcells = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  Maxcells = MAX( 4 , Maxcells*2 )
#else
  !Serial version: use max. 8 cells along any given direction
  Maxcells = MIN( 8 , NINT( SIZE(A,1)**(1.d0/3.d0) ) )
#endif
  WRITE(msg,*) 'Max. allowed number of cells along any direction: ', Maxcells
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !
  !Determine the number of cells needed along each dimension X, Y, Z
  DO i=1,3
    distance = SUM(DABS(H(:,i)))
    tempreal = 1.9d0*distance/R
    IF( distance < 2.d0*R .OR. NINT(tempreal) < 2 ) THEN
      !Ensure that there is always at least one cell along any given direction
      NcellsX(i) = 1
    ELSE
      !Large cell or small R => do not create zillions of small cells
      !=> No more than Maxcell cells along any direction
      NcellsX(i) = MIN( NINT(tempreal) , Maxcells )
    ENDIF
    !Length of a cell along each direction
    Cell_L(i) = distance / DBLE(NcellsX(i))
  ENDDO
  Ncells = PRODUCT( NcellsX(:) )  !total number of cells
  WRITE(msg,*) 'Actual number of cells along X, Y, Z: ', NcellsX(1), NcellsX(2), NcellsX(3)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !WRITE(*,*) '    Maxcells = ', Maxcells
  !WRITE(*,*) ' Total cells = ', Ncells
  !
  !For each atom, find the index of the cell it belongs to
  ALLOCATE( Atom_Cell(SIZE(A,1)) )  !index of cell each atom belongs to
  Atom_Cell(:) = 0
  ALLOCATE( Cell_NP(Ncells) )       !number of atoms in each cell
  Cell_NP(:) = 0
  !$OMP PARALLEL DO DEFAULT(SHARED) &
  !$OMP& PRIVATE(i,Ix,Iy,Iz) REDUCTION(+:Atom_Cell,Cell_NP)
  DO i=1,SIZE(A,1)
    Ix = MAX( CEILING(A(i,1)/Cell_L(1)) , 1 )
    Iy = MAX( CEILING(A(i,2)/Cell_L(2)) , 1 )
    Iz = MAX( CEILING(A(i,3)/Cell_L(3)) , 1 )
    !Link this atom to the corresponding cell
    Atom_Cell(i) = (Iz-1)*NcellsX(1)*NcellsX(2) + (Iy-1)*NcellsX(1) + Ix
    IF( Atom_Cell(i)<=0 ) Atom_Cell(i) = 1
    IF( Atom_Cell(i)>SIZE(Cell_NP) ) Atom_Cell(i) = SIZE(Cell_NP)
    !Increment number of atoms of this cell
    Cell_NP(Atom_Cell(i)) = Cell_NP(Atom_Cell(i)) + 1
  ENDDO
  !$OMP END PARALLEL DO
  !
  IF( verbosity==4 ) THEN
    OPEN(UNIT=56,FILE="atomsk_atomCell.txt",STATUS="UNKNOWN")
    DO i=1,SIZE(A,1)
      WRITE(56,*) i, Atom_Cell(i)
    ENDDO
    CLOSE(56)
  ENDIF
  !
  !Save the indexes of atoms in each cell
  ALLOCATE( Cell_AtomID( Ncells , MAXVAL(Cell_NP) ) )
  Cell_AtomID(:,:) = 0
  !$OMP PARALLEL DO DEFAULT(SHARED) &
  !$OMP& PRIVATE(i,j,k)
  DO j=1,Ncells
    k=0
    DO i=1,SIZE(A,1)
      IF( Atom_Cell(i)==j ) THEN
        !Atom #i belongs to cell #j
        k=k+1
        Cell_AtomID(j,k) = i
      ENDIF
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  !
  IF( verbosity==4 ) THEN
    OPEN(UNIT=56,FILE="atomsk_cellAtomID.txt",STATUS="UNKNOWN")
    DO j=1,Ncells
      WRITE(56,*) j, Cell_AtomID(j,:)
    ENDDO
    CLOSE(56)
  ENDIF
  !
  !Make a neighbor list for the cells
  !For each cell we consider a cube of 3*3*3 = 27 cells (central cell + neighboring cells)
  !For each of the 27 cells, we store 4 numbers: the cell index, and the coordinates of the
  !appropriate periodic image (in units of box vectors)
  ALLOCATE( Cell_Neigh(Ncells,27,4) )
  Cell_Neigh(:,:,:) = 0
  DO i=1,NcellsX(1)
    DO j=1,NcellsX(2)
      DO k=1,NcellsX(3)
        u=0
        iCell = (k-1)*NcellsX(1)*NcellsX(2) + (j-1)*NcellsX(1) + i
        !
        DO l=i-1,i+1
          a1=l
          shift(1) = 0
          IF( a1>NcellsX(1) ) THEN
            !No cell exists to the right of cell #iCell
            !The neighboring cell is actually the periodic image of cell #1 shifted by +H(1,:)
            a1 = 1
            shift(1) = 1
          ELSEIF( a1<=0 ) THEN
            !No cell exists to the left of cell #iCell
            !The neighboring cell is actually the periodic image of cell #NcellsX(1) shifted by -H(1,:)
            a1 = NcellsX(1)
            shift(1) = -1
          ENDIF
          !
          DO m=j-1,j+1
            a2=m
            shift(2) = 0
            IF( a2>NcellsX(2) ) THEN
              !No cell exists on top of cell #iCell
              !The neighboring cell is actually the periodic image of cell #1 shifted by +H(2,:)
              a2 = 1
              shift(2) = 1
            ELSEIF( a2<=0 ) THEN
              !No cell exists below cell #iCell
              !The neighboring cell is actually the periodic image of cell #NcellsX(2) shifted by -H(2,:)
              a2 = NcellsX(2)
              shift(2) = -1
            ENDIF
            !
            DO n=k-1,k+1
              a3=n
              shift(3) = 0
              IF( a3>NcellsX(3) ) THEN
                !No cell exists in front of cell #iCell
                !The neighboring cell is actually the periodic image of cell #1 shifted by +H(3,:)
                a3 = 1
                shift(3) = 1
              ELSEIF( a3<=0 ) THEN
                !No cell exists behind cell #iCell
                !The neighboring cell is actually the periodic image of cell #NcellsX(3) shifted by -H(3,:)
                a3 = NcellsX(3)
                shift(3) = -1
              ENDIF
              !
              u = u+1
              Cell_Neigh(iCell,u,1) = (a3-1)*NcellsX(1)*NcellsX(2) + (a2-1)*NcellsX(1) + a1
              Cell_Neigh(iCell,u,2) = shift(1)
              Cell_Neigh(iCell,u,3) = shift(2)
              Cell_Neigh(iCell,u,4) = shift(3)
            ENDDO !n
          ENDDO !m
        ENDDO !l
        !
      ENDDO !k
    ENDDO !j
  ENDDO !i
  WRITE(msg,*) 'Neighbor list for cells complete'
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !
  !Construct the neighbor list for atoms
  !$OMP PARALLEL DO DEFAULT(SHARED) &
  !$OMP& PRIVATE(i,iCell,j,k,n,u,distance,Vfrac)
  DO i=1,SIZE(A,1)
    !iCell = index of the cell atom #i belongs to
    iCell = Atom_Cell(i)
    !
    !Parse atoms in cell #iCell and its neighbors
    DO j=1,27  !SIZE(Cell_Neigh,2)
      !Parse atoms in cell #j (there are Cell_NP(j) atoms in it)
      DO k=1,Cell_NP(Cell_Neigh(iCell,j,1))
        !Cell_Neigh(iCell,j,1) is the index of the j-th neighboring cell of cell #iCell
        !n = actual index of the k-th atom in that cell
        n = Cell_AtomID( Cell_Neigh(iCell,j,1) , k )
        !Do not count atom #i as its own neighbor (i.e. n!=i)
        IF( n>0 .AND. n<=SIZE(A,1) .AND. n.NE.i ) THEN
          !Check if this atom was already counted as neighbor
          IF( .NOT. ANY(NeighList(i,:)==n) ) THEN
            !Use the correct periodic image of that cell
            !The position of the atom #n has to be translated by Cell_Neigh(iCell,j,2:4) * (box vectors)
            Vfrac(1,1:3) = A(n,1:3) + DBLE(Cell_Neigh(iCell,j,2))*H(1,:)   &
                         &          + DBLE(Cell_Neigh(iCell,j,3))*H(2,:)   &
                         &          + DBLE(Cell_Neigh(iCell,j,4))*H(3,:)
            !Compute distance between atom #i and the periodic image of atom #n
            distance = VECLENGTH( A(i,1:3) - Vfrac(1,1:3) )
            IF( distance < R ) THEN
              !Add atom #n as neighbor of atom #i
              NNeigh(i) = NNeigh(i)+1
              IF( Nneigh(i) <= SIZE(NeighList,2) ) THEN
                NeighList(i,Nneigh(i)) = n
              ENDIF
              !
              IF( .NOT. ANY(NeighList(n,:)==i) ) THEN
                !Also add atom #i as neighbor of atom #n
                !$OMP CRITICAL
                NNeigh(n) = NNeigh(n)+1
                !$OMP END CRITICAL
                IF( Nneigh(n) <= SIZE(NeighList,2) ) THEN
                  NeighList(n,Nneigh(n)) = i
                ENDIF
              ENDIF
            ENDIF
            !
          ENDIF
        ENDIF
        !
      ENDDO  !k
      !
    ENDDO  !j
    !
  ENDDO  !i
  !$OMP END PARALLEL DO
  !
  !Free memory
  IF(ALLOCATED(Cell_NP)) DEALLOCATE(Cell_NP)
  IF(ALLOCATED(Cell_AtomID)) DEALLOCATE(Cell_AtomID)
  IF(ALLOCATED(Cell_Neigh)) DEALLOCATE(Cell_Neigh)
  !
  !
ENDIF
!
!
IF(ALLOCATED(Nneigh)) DEALLOCATE(Nneigh)
!
IF( ALLOCATED(NeighList) ) THEN
  !If the neighbor list contains only zeros, then no neighbor was found
  !=> deallocate NeighList
  IF( .NOT.ANY(NeighList(:,:).NE.0) ) THEN
    DEALLOCATE(NeighList)
  ENDIF
ENDIF
!
IF( ALLOCATED(NeighList) ) THEN
  WRITE(msg,*) 'Max. n. neighbors = ', SIZE(NeighList,2)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ELSE
  msg = 'NeighList UNALLOCATED'
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ENDIF
!
IF( verbosity==4 ) THEN
  !Verbose mode: write neighbor list in a text file
  IF( ALLOCATED(NeighList) .AND. SIZE(NeighList,1)>0 ) THEN
    OPEN(UNIT=51,FILE="atomsk_neighborlist.txt",STATUS="UNKNOWN",FORM="FORMATTED")
    DO i=1,SIZE(NeighList,1)
      WRITE(51,*) (NeighList(i,j),j=1,SIZE(NeighList,2))
    ENDDO
    CLOSE(51)
  ENDIF
ENDIF
!
msg = 'exiting NEIGHBOR_LIST'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
END SUBROUTINE NEIGHBOR_LIST
!
!
!********************************************************
! NEIGHBOR_POS
! Given an atom position V, the list of indices of its
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
! NOTE: as for the subroutine NEIGHBOR_LIST, atom
!      coordinates are assumed to be wrapped.
!********************************************************
SUBROUTINE NEIGHBOR_POS(H,A,V,NeighList,Use_NeighList,radius,PosList)
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
LOGICAL:: selfneighbor  !is central atom a neighbor of itself? (because of PBC)
LOGICAL,INTENT(IN):: Use_NeighList !use provided neighbor list?
INTEGER:: i, m, n, o
INTEGER:: Nm, Nn, No !number of replicas to search along X, Y, Z
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
!
selfneighbor=.FALSE.
IF(ALLOCATED(PosList)) DEALLOCATE(PosList)
!
Nm = MAX( 1 , CEILING(radius/VECLENGTH(H(1,:))) +1 )
Nn = MAX( 1 , CEILING(radius/VECLENGTH(H(2,:))) +1 )
No = MAX( 1 , CEILING(radius/VECLENGTH(H(3,:))) +1 )
!
Nneighbors=0
IF( Use_NeighList ) THEN
  i=1
  DO WHILE( i<=SIZE(NeighList) .AND. NeighList(i)>0 )
    !NeighList(i) is the index of a neighboring atom
    !Find which periodic image(s) are actually neighbors
    DO o=-No,No
      DO n=-Nn,Nn
        DO m=-Nm,Nm
          !Position of the periodic image of neighbor #i
          P1 = A(NeighList(i),1) + DBLE(m)*H(1,1) + DBLE(n)*H(2,1) + DBLE(o)*H(3,1)
          P2 = A(NeighList(i),2) + DBLE(m)*H(1,2) + DBLE(n)*H(2,2) + DBLE(o)*H(3,2)
          P3 = A(NeighList(i),3) + DBLE(m)*H(1,3) + DBLE(n)*H(2,3) + DBLE(o)*H(3,3)
          vector = (/ P1 , P2 , P3 /)
          !This image is a neighbor if distance is smaller than radius
          distance = VECLENGTH( vector(:) - V(:) )
          IF( distance <= radius ) THEN
            Nneighbors = Nneighbors+1
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    i=i+1
  ENDDO
ENDIF
!
DO o=-No,No
  DO n=-Nn,Nn
    DO m=-Nm,Nm
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
  IF( Use_NeighList ) THEN
    i=1
    DO WHILE( i<=SIZE(NeighList) .AND. NeighList(i)>0 )
      DO o=-No,No
        DO n=-Nn,Nn
          DO m=-Nm,Nm
            !Position of the periodic image of neighbor #i
            P1 = A(NeighList(i),1) + DBLE(m)*H(1,1) + DBLE(n)*H(2,1) + DBLE(o)*H(3,1)
            P2 = A(NeighList(i),2) + DBLE(m)*H(1,2) + DBLE(n)*H(2,2) + DBLE(o)*H(3,2)
            P3 = A(NeighList(i),3) + DBLE(m)*H(1,3) + DBLE(n)*H(2,3) + DBLE(o)*H(3,3)
            vector = (/ P1 , P2 , P3 /)
            !This image is a neighbor if distance is smaller than radius
            distance = VECLENGTH( vector(:) - V(:) )
            IF( distance <= radius ) THEN
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
  ENDIF
  !
  IF( selfneighbor ) THEN
    DO o=-No,No
      DO n=-Nn,Nn
        DO m=-Nm,Nm
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
