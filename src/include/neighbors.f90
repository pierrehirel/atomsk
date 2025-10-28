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
!* Last modification: P. Hirel - 01 Oct. 2025                                     *
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
!* NEIGHBOR_LIST       runs Verlet or cell neighbor list algorithm                *
!* VERLET_LIST         constructs a neighbor list (Verlet algorithm)              *
!* CELL_LIST           constructs a neighbor list (cell decomposition)            *
!* NEIGHBOR_POS        given a neighbor list, gives the coordinates of neighbors  *
!* FIND_NNN            finds the N nearest neighbors of a coordinate              *
!* FIND_NNRdR          finds the neighbors in a skin of radius R and width dR     *
!* FIND_1NN            finds the first nearest neighbors                          *
!**********************************************************************************
!
!
USE comv
USE subroutines
USE messages
USE resize
!
!
CONTAINS
!
!
!
!********************************************************
! NEIGHBOR_LIST
! Considering a list of atom positions A, base
! vectors H defining periodic conditions, and a radius R,
! a list of neighbors is constructed for all atoms,
! i.e. a list that looks like the following:
!        2 4 6 0 0 0 0
!        1 3 4 5 0 0 0
!        2 4 5 0 0 0 0
!        1 2 3 5 6 0 0
!    etc.
! meaning that neighbors of atom #1 are atoms #2, 4, 6,
! neighbors of atom #2 are #1, 3, 4, 5, etc.
!     NOTES:
! (1) This subroutine is just an alias that automatically
!     calls VERLET_LIST or CELL_LIST, depending on criteria.
!     It remains possible to call VERLET_LIST or CELL_LIST
!     directly, if desired/necessary.
! (2) atoms are assumed to be wrapped, i.e. all 
!     in the box (no atom outside of the box).
! (3) the neighbor list NeighList(i,:) of the atom #i
!     may contain trailing zeros, as illustrated above.
!     These zeros are meant to be ignored. The first zero
!     encountered means there will be no more neighbors.
! (4) an atom is *never* counted as its own neighbor.
!     A neighbor appears only once in the neighbor list.
!     If the system is such that
!     replicas of atom #i are neighbors of atom #i,
!     or such that several replica of an atom #j can
!     be neighbors to an atom #i, then the replica do
!     not appear in the list, i.e. the atom #j is
!     counted only once. If several replica of atom #j
!     are neighbors of atom #i, that must be corrected for
!     *after* calling this routine (see NEIGHBOR_POS).
! (5) if the number of atoms is small, then a
!     simplistic Verlet search is performed, which
!     scales as N². Otherwise a cell decomposition
!     algorithm is used, which scales as N.
! (6) if no neighbor is found, then the array
!     NeighList is returned as UNALLOCATED.
!********************************************************
SUBROUTINE NEIGHBOR_LIST(H,A,R,NeighList)
!
IMPLICIT NONE
REAL(dp):: distance
REAL(dp),INTENT(IN):: R  !radius in which Neighbors are searched
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),INTENT(IN):: A  !array of all atom positions
INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: NeighList  !the neighbor list
!
IF(ALLOCATED(NeighList)) DEALLOCATE(NeighList)
!
!Compute minimum cell length
distance = MIN( VECLENGTH(H(1,:)) , VECLENGTH(H(2,:)) , VECLENGTH(H(3,:)) )
!
neighsearch = "verlet"  !TEMPORARY: force use of Verlet algorithm
IF( STRDNCASE(neighsearch)=="verlet" ) THEN
  !User forces use of Verlet algorithm
  CALL VERLET_LIST(H,A,R,NeighList)
  !
ELSEIF( STRDNCASE(neighsearch)=="cell" ) THEN
  !User forces use of CELL algorithm
  CALL CELL_LIST(H,A,R,NeighList)
  !
ELSE
  !User did not specify anything: automatically choose VERLET or CELL depending on system size
  !IF( distance>10.d0 .OR. R>10.d0 .OR. SIZE(A,1)>10000 ) THEN
  IF( R>VECLENGTH(H(1,:)) .AND. R>VECLENGTH(H(2,:)) .AND. R>VECLENGTH(H(3,:)) &
    & .AND. SIZE(A,1)>20000 )  THEN
    !Very large system => use cell decomposition algorithm
    CALL CELL_LIST(H,A,R,NeighList)
  ELSE
    !Default: use Verlet algorithm
    CALL VERLET_LIST(H,A,R,NeighList)
  ENDIF
  !
ENDIF
!
END SUBROUTINE NEIGHBOR_LIST
!
!
!
!********************************************************
! VERLET_LIST
! This subroutine constructs a neighbor list using
! the most simple Verlet algorithm. It is expected
! to be robust, but scales poorly with system size (O(N²)).
!********************************************************
SUBROUTINE VERLET_LIST(H,A,R,NeighList)
!
IMPLICIT NONE
CHARACTER(LEN=4096):: msg, temp
REAL(dp),INTENT(IN):: R  !radius in which Neighbors are searched
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),INTENT(IN):: A  !array of all atom positions
!
LOGICAL:: frac 
INTEGER:: i, j, k, m, n
INTEGER:: progress      !To show calculation progress
INTEGER,PARAMETER:: NNincrement=2  !whenever list is full, increase its size by that much
INTEGER,DIMENSION(:,:),ALLOCATABLE:: tempList    !list of neighbours (temporary)
REAL(dp):: distance   !distance between two atoms
REAL(dp):: dx         !difference of coordinates along a direction
REAL(dp):: rho        !average density of the system
REAL(dp):: Vsystem   !volume of the box defined by H(:,:)
REAL(dp),DIMENSION(3):: shift      !shift due to periodic boundary conditions
REAL(dp),DIMENSION(3):: V, Vi, Vj  !fractional coordinates
REAL(dp),DIMENSION(3,3):: G      
INTEGER,DIMENSION(:),ALLOCATABLE:: NNeigh      !number of neighbors of atom #i
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Afrac    !positions in reduced coordinates
!
INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: NeighList  !the neighbor list
!
!
msg = 'entering VERLET_LIST'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,*) 'Radius for neighbor search (angstroms) = ', R
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Initialize variables
frac=.FALSE.
IF(ALLOCATED(tempList)) DEALLOCATE(tempList)
IF(ALLOCATED(NNeigh)) DEALLOCATE(NNeigh)
IF(ALLOCATED(NeighList)) DEALLOCATE(NeighList)
ALLOCATE(NNeigh(SIZE(A,1)) , STAT=i)
IF( i>0 ) THEN
  ! Allocation failed (not enough memory)
  nerr = nerr+1
  CALL ATOMSK_MSG(819,(/''/),(/-1.d0/))
  RETURN
ENDIF
NNeigh(:) = 0
!Assume a continuous atom density to estimate initial size of NeighList
!NOTE: if an atom has a greater number of neighbors the size of NeighList will be changed later
Vsystem = VOLUME_PARA(H)
rho = DBLE(SIZE(A,1))/Vsystem
n = MAX( 100 , CEILING(rho*1.5d0*pi*(R**3)) )
WRITE(msg,*) "N. neighbors/atom estimated from density = ", n
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
IF(n<=0) THEN
  nerr = nerr+1
  RETURN
ENDIF
ALLOCATE(NeighList(SIZE(A,1),n) , STAT=i)
IF( i>0 ) THEN
  ! Allocation failed (not enough memory)
  nerr = nerr+1
  CALL ATOMSK_MSG(819,(/''/),(/-1.d0/))
  RETURN
ENDIF
NeighList(:,:) = 0
!
!
!NOTE: storing reduced coordinates in Afrac(:,:) is faster, but uses more memory.
!If allocation fails (not enough memory/too many atoms), then reduced coordinates
!will be computed on-the-fly in the loop below, but this is slower
ALLOCATE( Afrac(SIZE(A,1),SIZE(A,2)) , STAT=i )
IF( i==0 ) THEN
  ! Allocation succeeded
  frac = .TRUE.
  !Convert atom positions into reduced coordinates
  Afrac(:,:) = A(:,:)
  CALL CART2FRAC(Afrac,H)
ELSE
  ! Allocation failed (not enough memory)
  frac = .FALSE.
  !We will need the inverted cell vector matrix
  CALL INVMAT(H,G)
ENDIF
!
!Loop on all atoms to find their neighbors
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j,k,V,Vi,Vj,dx,shift,distance)
DO i=1,SIZE(A,1)-1
  !
  progress = progress+1
  !
  IF( SIZE(A,1) > 10000 ) THEN
    !If there are many atoms, display a fancy progress bar
    CALL ATOMSK_MSG(10,(/""/),(/DBLE(progress),DBLE(SIZE(A,1))/))
  ENDIF
  !
  !Convert position of atom #i into fractional coordinates
  IF( frac ) THEN
    Vi(:) = Afrac(i,1:3)
  ELSE
    V(:) = A(i,1:3)
    DO k=1,3
      Vi(k) = V(1)*G(1,k) + V(2)*G(2,k) + V(3)*G(3,k)
    ENDDO
  ENDIF
  !
  DO j=i+1,SIZE(A,1)
    !Convert position of atom #i into fractional coordinates
    IF( frac ) THEN
      Vj(:) = Afrac(j,1:3)
    ELSE
      V(:) = A(j,1:3)
      DO k=1,3
        Vj(k) = V(1)*G(1,k) + V(2)*G(2,k) + V(3)*G(3,k)
      ENDDO
    ENDIF
    !
    shift(:) = 0.d0
    !
    !Compute distance (reduced coordinates) between atoms #i and #j
    !Account for periodic boundary conditions
    DO k=1,3  !loop on all 3 directions
      dx = DABS( Vj(k) - Vi(k) )
      IF( DABS(1.d0-dx) < dx ) THEN
        IF( A(j,k)>A(i,k) ) THEN
          shift(k) = -1.d0
        ELSE
          shift(k) = 1.d0
        ENDIF
      ENDIF
    ENDDO
    !
    !Save shift of atom #j in Cartesian coordinates
    shift(:) = shift(1)*H(1,:) + shift(2)*H(2,:) + shift(3)*H(3,:)
    !
    !Compute distance (Cartesian coordinates) between atom #i and atom #j (or replica of #j)
    distance = VECLENGTH( A(j,1:3) + shift(:) - A(i,1:3) )
    !
    IF ( distance <= R ) THEN
      !$OMP CRITICAL
      !IF( .NOT.ANY(NeighList(i,:)==j) ) THEN
        !Atom #j is neighbor of atom #i
        NNeigh(i) = NNeigh(i)+1
        !Add atom #j to the list of neighbors of atom #i
        IF( NNeigh(i) <= SIZE(NeighList,2) ) THEN
          NeighList(i,NNeigh(i)) = j
        ENDIF
      !ENDIF
      !
      !Atom #i is also neighbor of atom #j
      !IF( .NOT.ANY(NeighList(j,:)==i) ) THEN
        NNeigh(j) = NNeigh(j)+1
        !Add atom #i to the list of neighbors of atom #j
        IF( NNeigh(j) <= SIZE(NeighList,2) ) THEN
          NeighList(j,NNeigh(j)) = i
        ENDIF
      !ENDIF
      !$OMP END CRITICAL
    ENDIF
    !
  ENDDO !j
  !
ENDDO !i
!$OMP END PARALLEL DO
!
IF(ALLOCATED(Afrac)) DEALLOCATE(Afrac)
!
IF( ALLOCATED(NeighList) .AND. SIZE(NeighList,1)>0 ) THEN
  IF( .NOT.ANY(NeighList(:,:).NE.0) ) THEN
    !Neighbor list contains only zeros, i.e. no neighbor was found
    !=> deallocate NeighList
    DEALLOCATE(NeighList)
  ENDIF
ENDIF
!
IF( ALLOCATED(NeighList) .AND. SIZE(NeighList,1)>0 ) THEN
  !Count max. number of neighbors
  n = MAXVAL(NNeigh(:))
  WRITE(msg,*) "Max. N. neighbors found = ", n
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  IF( n>0 .AND. n<SIZE(NeighList,2) ) THEN
    !Reduce NeighList
    WRITE(msg,*) "Reducing size of NeighList ..."
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    m = SIZE(NeighList,1)
    CALL RESIZE_INTARRAY2(NeighList,m,n,i)
  ENDIF
  IF( verbosity==4 ) THEN
    WRITE(msg,*) SIZE(NeighList,1)
    WRITE(temp,*) SIZE(NeighList,2)
    WRITE(msg,*) "Final size of NeighList :  "//TRIM(ADJUSTL(msg))//" x "//TRIM(ADJUSTL(temp))
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  ENDIF
ELSE
  msg = 'NeighList UNALLOCATED'
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ENDIF
!
IF( verbosity==4 ) THEN
  !Verbose mode: write neighbor list in a text file
  IF( ALLOCATED(NeighList) .AND. SIZE(NeighList,1)>0 ) THEN
    OPEN(UNIT=51,FILE="atomsk_verletneighborlist.txt",STATUS="UNKNOWN",FORM="FORMATTED")
    DO i=1,SIZE(NeighList,1)
      WRITE(temp,*) i
      msg = "#"//TRIM(ADJUSTL(temp))
      j=1
      DO WHILE( j<=MIN(SIZE(NeighList,2),10) .AND. NeighList(i,j)>0 )
        WRITE(temp,*) NeighList(i,j)
        msg = TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(temp))
        j=j+1
      ENDDO
      WRITE(51,*) TRIM(msg)
    ENDDO
    CLOSE(51)
  ENDIF
ENDIF
!
IF(ALLOCATED(NNeigh)) DEALLOCATE(NNeigh)
!
msg = 'exiting VERLET_LIST'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
END SUBROUTINE VERLET_LIST
!
!
!
!********************************************************
! CELL_LIST
! This subroutine constructs a neighbor list using
! a cell decomposition algorithm. It is expected to be
! much faster than Verlet algorithm for large systems.
!********************************************************
SUBROUTINE CELL_LIST(H,A,R,NeighList)
!
IMPLICIT NONE
CHARACTER(LEN=4096):: msg, temp
REAL(dp),INTENT(IN):: R  !radius in which Neighbors are searched
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),INTENT(IN):: A  !array of all atom positions
!
INTEGER:: Ix, Iy, Iz
INTEGER:: i, j, k, l, m, n
INTEGER:: iCell  !index of a cell
INTEGER:: Maxcells !max number of cells along one direction
INTEGER:: Ncells   !total number of cells
INTEGER:: OMP_GET_NUM_THREADS  !number of OpenMP threads
INTEGER,PARAMETER:: NNincrement=2  !whenever list is full, increase its size by that much
INTEGER,DIMENSION(3):: NcellsX  !number of cells along X, Y, Z (cell-list algorithm)
INTEGER,DIMENSION(SIZE(A,1)):: LinkedList
INTEGER,DIMENSION(:),ALLOCATABLE:: LinkedHead
INTEGER,DIMENSION(:),ALLOCATABLE:: Cell_NP       !number of atoms in each cell
INTEGER,DIMENSION(:),ALLOCATABLE:: NNeigh        !number of neighbors of atom #i
INTEGER,DIMENSION(:),ALLOCATABLE:: Atom_Cell     !for each atom, index of the cell it belongs to
INTEGER,DIMENSION(:,:),ALLOCATABLE:: Cell_AtomID !index of atoms for each cell
INTEGER,DIMENSION(:,:),ALLOCATABLE:: Cell_ijk    !index of each cell
INTEGER,DIMENSION(:,:),ALLOCATABLE:: tempList    !list of neighbours (temporary)
INTEGER,DIMENSION(:,:),ALLOCATABLE:: Cell_Neigh  !neighbors of each cell
REAL(dp):: distance, dx
REAL(dp):: tempreal
REAL(dp):: Vsystem   !volume of the box defined by H(:,:)
REAL(dp),DIMENSION(3):: Cell_L   !length of cell along X, Y, Z (cell-list algorithm)
REAL(dp),DIMENSION(3):: shift   !replica number along X, Y, Z
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Cell_P  !positions of cells
REAL(dp),DIMENSION(SIZE(A,1),SIZE(A,2)):: Afrac  !atom positions in reduced coordinates
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
LinkedList(:) = 0
IF(ALLOCATED(NeighList)) DEALLOCATE(NeighList)
!Assume a continuous atom density to estimate initial size of NeighList
Vsystem = VOLUME_PARA(H)
n = MAX( 100 , NINT((DBLE(SIZE(A,1))/Vsystem)*1.5d0*pi*(R**3))+100 )
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
IF( verbosity==4 ) THEN
  msg = "entering CELL_LIST"
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,*) "Radius for neighbor search (angstroms) = ", R
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,*) "Max. allowed number of neighbors for each atom = ", SIZE(NeighList,2)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,*) "Current box vectors:"
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  DO i=1,3
    WRITE(msg,'(3(3X,f9.3))') H(i,1), H(i,2), H(i,3)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  ENDDO
ENDIF
!
!Define max. number of cells along each direction
!Use max. 8 cells along any given direction
Maxcells = MIN( 16 , NINT( SIZE(A,1)**(1.d0/3.d0) ) )
!Make sure that Maxcells is even
IF( MOD(Maxcells,2).NE.0 .AND. Maxcells>1 ) Maxcells = Maxcells-1
WRITE(msg,*) "Max. allowed number of cells along any direction: ", Maxcells
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Determine the number of cells along each dimension X, Y, Z
DO i=1,3
  !Compute max. distance between two ends of the box
  distance = SUM(DABS(H(:,i)))  !MAXVAL(DABS(H(:,i)))
  !Compute how many cells of length R we can fit in that distance
  tempreal = 0.9d0*distance/R
  IF( distance < R .OR. NINT(tempreal) < 2 ) THEN
    !Cell is roughly the same size as cutoff radius
    !Ensure that there is always at least one cells along any given direction
    NcellsX(i) = 1
  ELSE
    !Large cell or small R => use R as cell size
    !BUT do not create zillions of small cells
    !=> No more than Maxcell cells along any direction
    NcellsX(i) = MIN( NINT(tempreal) , Maxcells )
  ENDIF
  !Make sure that number of cells is even (or 1)
  IF( MOD(NcellsX(i),2).NE.0 .AND. NcellsX(i)>1 ) NcellsX(i) = NcellsX(i)-1
  !Make sure it does not exceed Maxcells
  IF( NcellsX(i)>Maxcells ) NcellsX(i) = Maxcells
  !
  !TEMPORARY: ENFORCE NCELLS ALONG EACH DIRECTION
  !NcellsX(:) = 2
  !
  !Length of a cell along current direction
  Cell_L(i) = distance / DBLE(NcellsX(i))
ENDDO
!
IF( verbosity==4 ) THEN
  !Debug message
  WRITE(msg,*) "Actual number of cells along X, Y, Z: ", NcellsX(1), NcellsX(2), NcellsX(3)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,*) "Actual cell size along X, Y, Z: ", Cell_L(1), Cell_L(2), Cell_L(3)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ENDIF
!PRINT*, "CELL DECOMPOSITION: ", NcellsX(1), NcellsX(2), NcellsX(3)
!
!Compute total number of cells
Ncells = PRODUCT( NcellsX(:) )  !total number of cells
!
!Compute position of each cell
ALLOCATE( Cell_P(Ncells,3) )
Cell_P(:,:) = 0.d0
ALLOCATE( Cell_ijk(Ncells,3) )
Cell_ijk(:,:) = 0
iCell=0
DO i=1,NcellsX(1)
  DO j=1,NcellsX(2)
    DO k=1,NcellsX(3)
      iCell=iCell+1
      Cell_ijk(iCell,1) = i
      Cell_ijk(iCell,2) = j
      Cell_ijk(iCell,3) = k
      Cell_P(iCell,1) = (i-1)*Cell_L(1)
      Cell_P(iCell,2) = (j-1)*Cell_L(2)
      Cell_P(iCell,3) = (k-1)*Cell_L(3)
    ENDDO
  ENDDO
ENDDO
!
!For each cell, find its neighbouring cells in a radius R+(max.cell length)
!Number of cells is quite small, use Verlet list
distance = 1.1d0*DSQRT( Cell_L(1)**2 + Cell_L(2)**2 + Cell_L(3)**2 )
CALL VERLET_LIST(H,Cell_P,R+distance,Cell_Neigh)
!For each cell, append its own index at beginning of neighbour list
IF( ALLOCATED(Cell_Neigh) .AND. SIZE(Cell_Neigh,2)>0 ) THEN
  !Some neighbours were found
  ALLOCATE( tempList( SIZE(Cell_Neigh,1) , SIZE(Cell_Neigh,2)+1 ) )
  tempList(:,:) = 0
  DO i=1,SIZE(Cell_Neigh,1)
    tempList(i,1) = i
    DO j=1,SIZE(Cell_Neigh,2)
      tempList(i,j+1) = Cell_Neigh(i,j)
    ENDDO
  ENDDO
  DEALLOCATE(Cell_Neigh)
  ALLOCATE( Cell_Neigh( SIZE(tempList,1) , SIZE(tempList,2) ) )
  Cell_Neigh(:,:) = tempList(:,:)
  DEALLOCATE(tempList)
ELSE
  !No neighbour was found (e.g. because there is only 1 cell)
  !Simply add cell to its own cell list
  IF(ALLOCATED(Cell_Neigh)) DEALLOCATE(Cell_Neigh)
  ALLOCATE(Cell_Neigh(Ncells,1))
  DO i=1,Ncells
    Cell_Neigh(i,1) = i
  ENDDO
ENDIF
!
IF( verbosity==4 ) THEN
  !debug messages
  WRITE(msg,*) 'Cell neighbor list complete'
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !debug: for each cell, write its position in an XYZ file
  OPEN(UNIT=56,FILE="atomsk_cellP.xyz",STATUS="UNKNOWN")
  WRITE(56,*) Ncells
  WRITE(56,*) "# Positions of cells in CELL algorithm"
  DO j=1,Ncells
    WRITE(56,*) "1  ", Cell_P(j,1:3)
  ENDDO
  CLOSE(56)
  !debug: for each cell, write its neighbor list in a file
  OPEN(UNIT=56,FILE="atomsk_cellNeighList.txt",STATUS="UNKNOWN")
  DO j=1,Ncells
    WRITE(msg,*) j
    msg = "[cell #"//TRIM(ADJUSTL(msg))//"]"
    DO n=1,SIZE(Cell_Neigh,2)
      WRITE(temp,'(i4,a1,3(i2,a1))') Cell_Neigh(j,n)
      msg = TRIM(msg)//"  "//TRIM(ADJUSTL(temp))
    ENDDO
    WRITE(56,'(a)') TRIM(msg)
  ENDDO
  CLOSE(56)
ENDIF
!
!We don't need cell positions anymore: free memory
IF(ALLOCATED(Cell_P)) DEALLOCATE(Cell_P)
!
WRITE(msg,*) 'Constructing Linked List ...'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!Associate each atom with its cell
!This loop is expected to be fast and to scale as O(N)
ALLOCATE( LinkedHead(Ncells) )
LinkedHead(:) = 0
DO i=1,SIZE(A,1)
  !iCell = index of cell atom #i belongs to
  Ix = MAX( CEILING(A(i,1)/Cell_L(1)) , 1 )
  Iy = MAX( CEILING(A(i,2)/Cell_L(2)) , 1 )
  Iz = MAX( CEILING(A(i,3)/Cell_L(3)) , 1 )
  iCell = (Iz-1)*NcellsX(1)*NcellsX(2) + (Iy-1)*NcellsX(1) + Ix
  !PRINT*, "[ Atom #", i, " , Cell #", iCell, " ]"
  !In the LinkedList of atom #i, save 
  LinkedList(i) = LinkedHead(iCell)
  !Make LinkedHead(iCell) point to current atom #i
  LinkedHead(iCell) = i
ENDDO
WRITE(msg,*) 'Linked List finished.'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
IF( verbosity==4 ) THEN
  !debug: write LinkedHead into a file
  OPEN(UNIT=56,FILE="atomsk_LinkedHead.txt",STATUS="UNKNOWN")
  DO i=1,SIZE(LinkedHead)
    WRITE(56,*) i, LinkedHead(i)
  ENDDO
  CLOSE(56)
  !debug: write LinkedList into a file
  OPEN(UNIT=56,FILE="atomsk_LinkedList.txt",STATUS="UNKNOWN")
  DO i=1,SIZE(LinkedList)
    WRITE(56,*) i, LinkedList(i)
  ENDDO
  CLOSE(56)
ENDIF
!
WRITE(msg,*) "Constructing Atom Neighbour List ..."
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Save atom positions in fractional coordinates
Afrac(:,:) = A(:,:)
CALL CART2FRAC(Afrac,H)
!
!Using the previous linked list, construct actual neighbor list
!!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,iCell,Ix,Iy,Iz,j,k,n,shift,Vfrac,distance)
DO i=1,SIZE(A,1)  !Loop on all atoms
  !Atom #i belongs to cell #iCell
  Ix = MAX( CEILING(A(i,1)/Cell_L(1)) , 1 )
  Iy = MAX( CEILING(A(i,2)/Cell_L(2)) , 1 )
  Iz = MAX( CEILING(A(i,3)/Cell_L(3)) , 1 )
  iCell = (Iz-1)*NcellsX(1)*NcellsX(2) + (Iy-1)*NcellsX(1) + Ix
  !PRINT*, "[ Atom #", i, " , Cell #", iCell, " ]"
  !
  l=0
  !n = index of current neighbouring cell
  !Begin by searching for neighbours of atom #i in its own cell #iCell
  n = iCell  !Cell_Neigh(iCell,l)
  DO WHILE( l<SIZE(Cell_Neigh,2) .AND. n>0 )   !l=1,SIZE(Cell_Neigh,2) !Loop on all neighbouring cells
    !Get index of first atom in cell #n
    j = LinkedHead(n)
    !
    !Loop on all atoms in LinkedList, until we meet an index equal to 0
    DO WHILE(j>0)
      !PRINT*, "                Cell #", n, ":", j
      !IF(j==0) EXIT
      IF( j>i .AND. .NOT.ANY(NeighList(i,:)==j) ) THEN
        !Account for periodic boundary conditions
        shift(:) = 0.d0
        DO k=1,3  !loop on all 3 directions
          dx = DABS( Afrac(j,k) - Afrac(i,k) )
          IF( DABS(1.d0-dx) < dx ) THEN
            IF( A(j,k)>A(i,k) ) THEN
              shift(k) = -1.d0
            ELSE
              shift(k) = 1.d0
            ENDIF
          ENDIF
        ENDDO
        !
        !Save shift of atom #j in Cartesian coordinates
        shift(:) = shift(1)*H(1,:) + shift(2)*H(2,:) + shift(3)*H(3,:)
        !
        !Compute distance (Cartesian coordinates) between atom #i and atom #j (or replica of #j)
        distance = VECLENGTH( A(j,1:3) + shift(:) - A(i,1:3) )
        !
        IF ( distance <= R ) THEN
          !PRINT*, "                     is neighbour"
          !Atom #j is neighbour of atom #i
          NNeigh(i) = Nneigh(i)+1
          IF(NNeigh(i)<=SIZE(NeighList,2)) NeighList(i,NNeigh(i)) = j
          !Atom #i is also neighbour of atom #j
          NNeigh(j) = Nneigh(j)+1
          IF(NNeigh(j)<=SIZE(NeighList,2)) NeighList(j,NNeigh(j)) = i
        ENDIF
      ENDIF
      j = LinkedList(j)
    ENDDO  !j
    !
    !After parsing atoms in cell #iCell, parse atoms in neighbouring cells
    l = l+1
    !n = index of next neighbouring cell
    n = Cell_Neigh(iCell,l)
    !
  ENDDO  !l
ENDDO  !i
!!!$OMP END PARALLEL DO
WRITE(msg,*) "Atom Neighbour List COMPLETE ..."
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Free memory
IF(ALLOCATED(Cell_NP)) DEALLOCATE(Cell_NP)
IF(ALLOCATED(Cell_AtomID)) DEALLOCATE(Cell_AtomID)
IF(ALLOCATED(Cell_Neigh)) DEALLOCATE(Cell_Neigh)
!
!
IF( ALLOCATED(NeighList) ) THEN
  IF( .NOT.ANY(NeighList(:,:).NE.0) ) THEN
    !Neighbor list contains only zeros, i.e. no neighbor was found
    !=> deallocate NeighList
    DEALLOCATE(NeighList)
  ENDIF
ENDIF
!
!
IF( ALLOCATED(NeighList) .AND. SIZE(NeighList,1)>0 ) THEN
  !Count max. number of neighbors
  n = MAXVAL(NNeigh(:))
  WRITE(msg,*) "Max. N. neighbors found = ", n
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  IF( n>0 .AND. n<SIZE(NeighList,2) ) THEN
    !Reduce NeighList
    WRITE(msg,*) "Reducing size of NeighList ..."
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    m = SIZE(NeighList,1)
    CALL RESIZE_INTARRAY2(NeighList,m,n,i)
  ENDIF
  IF( verbosity==4 ) THEN
    WRITE(msg,*) SIZE(NeighList,1)
    WRITE(temp,*) SIZE(NeighList,2)
    WRITE(msg,*) "Final size of NeighList :  "//TRIM(ADJUSTL(msg))//" x "//TRIM(ADJUSTL(temp))
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  ENDIF
ELSE
  msg = 'NeighList UNALLOCATED'
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ENDIF
!
IF( verbosity==4 ) THEN
  !Verbose mode: write neighbor list in a text file
  IF( ALLOCATED(NeighList) .AND. SIZE(NeighList,1)>0 ) THEN
    OPEN(UNIT=51,FILE="atomsk_cellneighborlist.txt",STATUS="UNKNOWN",FORM="FORMATTED")
    DO i=1,SIZE(NeighList,1)
      WRITE(temp,*) i
      msg = "[Atom #"//TRIM(ADJUSTL(temp))//"]"
      j=1
      DO WHILE( j<=SIZE(NeighList,2) .AND. NeighList(i,j)>0 )
        WRITE(temp,*) NeighList(i,j)
        msg = TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(temp))
        j=j+1
      ENDDO
      WRITE(51,*) TRIM(msg)
    ENDDO
    CLOSE(51)
  ENDIF
ENDIF
!
IF(ALLOCATED(NNeigh)) DEALLOCATE(NNeigh)
IF(ALLOCATED(LinkedHead)) DEALLOCATE(LinkedHead)
!
msg = "exiting CELL_LIST"
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
END SUBROUTINE CELL_LIST
!
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
      IF( distance>1.d-3 .AND. distance <= radius ) THEN
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
          IF( distance>1.d-3 .AND. distance <= radius ) THEN
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
  ALLOCATE( V_NN( NNN,SIZE(A,2) ) )
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
    DO j=1,SIZE(A,1)  !Loop on all atoms
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
            IF( distance>1.d-12 .AND. distance<tempreal .AND. diff_from_prev_NN ) THEN
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
