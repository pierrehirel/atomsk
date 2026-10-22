MODULE mode_polycrystal
!
!**********************************************************************************
!*  MODE_POLYCRYSTAL                                                              *
!**********************************************************************************
!* This module constructs a polycrystal using a Voronoi tessellation.             *
!* Input: - one or several "seeds" (=unit cells, or more complex atomic systems)  *
!*          defined in Hs,Ps,Ss,AUXs.                                             *
!*        - number of atoms for each seed, NPs(:)                                 *
!*        - seed index for each node, idseed(:)                                   *
!*        - nodes positions in vnodes(:,:)                                        *
!*        - grain orientations (rotation matrices) in vorient(:,:,:)              *
!* Output:- final polycrystal defined in H, P, S,AUX                              *
!**********************************************************************************
!* (C) May 2013 - Pierre Hirel                                                    *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 21 Oct. 2026                                     *
!**********************************************************************************
!* OUTLINE:                                                                       *
!* 100        Read atom positions of seed (usually a unit cell) from ucfile       *
!* 200        Read parameters from vfile                                          *
!* 300        For each node, search neighboring vertices                          *
!* 400        Construct template grain                                            *
!* 500        Construct grains using Voronoi tesselation                          *
!* 600        Apply options and write final result to output file(s)              *
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
USE math
USE messages
USE files
USE subroutines
USE options
USE center
USE orient
USE writeout
USE voronoi
!
CONTAINS
!
!
SUBROUTINE POLYCRYS(prefix,Hs,Ps,Ss,AUXNAMES,AUXs,NPs,idseed,idmatrix,vnodes,vorient,clearance,H,P,S,AUX)
!
!
IMPLICIT NONE
!Input parameters
INTEGER,INTENT(IN):: idmatrix                  !Index of seed to use as a matrix (=0 if none)
INTEGER,DIMENSION(:),ALLOCATABLE:: NPs         !Number of atoms for each seed
INTEGER,DIMENSION(:),ALLOCATABLE:: idseed      !Seed index for each node
REAL(dp),INTENT(IN):: clearance                !clear atoms that close to the GB
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE:: Hs     !Base vectors of seeds
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE:: Ps, Ss !positions of atoms, shells in seeds
REAL(dp),DIMENSION(:,:,:),INTENT(IN):: AUXs    !auxiliary properties of atoms in seeds
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: vnodes   !cartesian coordinate of each node
REAL(dp),DIMENSION(:,:,:):: vorient            !rotation matrix for each node
!Output parameters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of aux.prop. of atoms
REAL(dp),DIMENSION(3,3):: H                            !Base vectors of final polycrystal
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: P, S !positions of atoms, shells in final polycrystal
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: AUX  !aux.prop. of atoms in final polycrystal
!
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: msg, prefix, temp
CHARACTER(LEN=4096):: distfile, idsizefile !name of file containing grain size distribution, grain sizes
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to output
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: newAUXNAMES !names of auxiliary properties of atoms (temporary)
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: doshells, doaux !are there shells, auxiliary properties in initial seed?
LOGICAL:: domatrix        !is this the matrix?
LOGICAL:: Hset            !are the box vectors H(:,:) defined?
LOGICAL:: exceeds100      !does the number of neighbours exceed 100?
LOGICAL:: keepatom        !keep an atom and save it to final system?
LOGICAL:: protectuc       !protect unit cell integrity?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT
LOGICAL,DIMENSION(:),ALLOCATABLE:: ingrain  !is an atom in a grain?
INTEGER,DIMENSION(:),ALLOCATABLE:: newindex !list of sorted indexes
INTEGER:: twodim        !=0 if system is 3-D, =1,2,3 if system is thin along x, y, z
INTEGER:: grainID       !position of the auxiliary property "grainID" in AUX
INTEGER:: seedID        !position of aux.prop. "seedID" in AUX
INTEGER:: i, j, k
INTEGER:: i2, j2, k2
INTEGER:: istart
INTEGER:: linenumber    !line number when reading a file
INTEGER:: m, n, o
INTEGER:: maxvertex     !max. number of vertices to look for (defined for 3-D or 2-D)
INTEGER:: NPexpect      !expected total number of atoms in the final system
INTEGER:: NPnew
INTEGER:: NP          !total number of atoms in the final system
INTEGER:: qi   !used to count atoms in a grain
INTEGER:: inode, jnode
INTEGER:: Nnodes  !number of nodes
INTEGER:: status
INTEGER,DIMENSION(3):: expandmatrix
INTEGER,DIMENSION(3):: Nsurv  !number of atoms that survived a loop
INTEGER,DIMENSION(:),ALLOCATABLE:: Nvertices !number of vertices for each node
INTEGER,DIMENSION(:),ALLOCATABLE:: NPgrains  !number of atoms in each grain
INTEGER,DIMENSION(:,:,:),ALLOCATABLE:: NeighListMask  !mask for neighboring vertices
REAL(dp):: distance    !distance between two points
REAL(dp):: maxdvertices !maximum distance between 2 vertices
REAL(dp):: H1, H2, H3  !temporary position
REAL(dp):: P1, P2, P3  !temporary position
REAL(dp):: seed_density !density of the seed (N.atoms/Volume)
REAL(dp):: Volume, Vmin, Vmax, Vstep  !min, max. volume occupied by a grain, step for grain size distribution
REAL(dp),DIMENSION(3):: GrainCenter !position of center of grain (differs from node position)
REAL(dp),DIMENSION(3):: shift       !shift vector
REAL(dp),DIMENSION(3):: vector      !vector between an atom and a node
REAL(dp),DIMENSION(3):: vnormal     !vector normal to grain boundary
REAL(dp),DIMENSION(3,3):: Huc       !Base vectors of the oriented unit cell or template
REAL(dp),DIMENSION(3,3):: Ht        !Base vectors of the oriented unit cell or template
REAL(dp),DIMENSION(3,3):: ORIENT    !crystalographic orientation
REAL(dp),DIMENSION(3,3):: rotmat    !rotation matrix
REAL(dp),DIMENSION(4,3):: verp      !positions (x,y,z) of some vertices
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:),ALLOCATABLE:: cell_volumes !volume of each Voronoi cell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Q, T      !positions of atoms, shells in a grain
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX_Q     !auxiliary properties of atoms in a grain
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX    !auxiliary properties of atoms (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: shiftgrain !translation vectors for each grain
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: PosList   !positions of nearest neighbours
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE:: vvertex !for each node, position of neighboring vertices
!
!
!Initialize variables
distfile = TRIM(TRIM_EXT(prefix))//"_size-dist.txt"
idsizefile = TRIM(TRIM_EXT(prefix))//"_id-size.txt"
doaux = .FALSE.
doshells = .FALSE.
protectuc = .FALSE.
m=0
Nnodes=SIZE(vnodes,1)  !number of nodes
twodim = 0    !assume system will be 3-D
expandmatrix(:) = 0
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
IF(ALLOCATED(NPgrains)) DEALLOCATE(NPgrains)
 C_tensor(:,:) = 0.d0
IF(ALLOCATED(P)) DEALLOCATE(P)
IF(ALLOCATED(shiftgrain)) DEALLOCATE(shiftgrain)
IF(ALLOCATED(Nvertices)) DEALLOCATE(Nvertices)
IF(ALLOCATED(ingrain)) DEALLOCATE(ingrain)
!
!
!
!
100 CONTINUE
!Check if polycrystal is 2-D
twodim = 0
DO i=1,3
  IF( H(i,i) < 1.1d0*MAXVAL(Hs(:,i,i)) ) THEN
    twodim = i
    EXIT
  ENDIF
ENDDO
!
!
CALL ATOMSK_MSG(4054,(/''/),(/DBLE(twodim),H(1,1),H(2,2),H(3,3),DBLE(Nnodes)/))
!
!Check that the user provided enough information
IF( .NOT. ANY( NINT(H).NE.0 ) ) THEN
  !User did not provide box, or a box with negative or zero size
  nerr=nerr+1
  CALL ATOMSK_MSG(4820,(/""/),(/0.d0/))
  GOTO 1000
ENDIF
!Check that number of nodes is not zero
IF(Nnodes<1) THEN
  CALL ATOMSK_MSG(4831,(/""/),(/0.d0/))
  nerr=nerr+1
  GOTO 1000
ENDIF
!Check that user provided at least one seed
IF( .NOT.ALLOCATED(Ps) ) THEN
  !No seed provided => can't do anything
  PRINT*, "X!X ERROR: no seed provided"
  nerr = nerr+1
  GOTO 1000
ENDIF
!Check if shells are defined
IF( ALLOCATED(Ss) .AND. SIZE(Ss,1)==SIZE(Ps,1) ) THEN
  doshells = .TRUE.
ENDIF
!Check if auxiliary properties are defined
IF( SIZE(AUXs,1)==SIZE(Ps,1) .AND. SIZE(AUXs,2)==SIZE(Ps,2) .AND. SIZE(AUXs,3)==SIZE(AUXNAMES) ) THEN
  doaux = .TRUE.
ENDIF
!Check idseed(:): size must be Nnodes, and values must be < number of seeds
IF( .NOT.ALLOCATED(idseed) ) THEN
  !No seed id defined: assume all nodes use the first seed
  WRITE(msg,*) "/!\ WARNING: no seed defined, allocating idseed(:)"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ALLOCATE(idseed(Nnodes))
  idseed(:) = 1
ELSEIF( SIZE(idseed).NE.Nnodes ) THEN
  !Seeds are defined but do not match number of nodes
  WRITE(msg,*) "/!\ WARNING: idseed(:) size is not Nnodes, resizing"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ALLOCATE(Nvertices(Nnodes))
  Nvertices(:) = 1
  DO i=1,MIN(Nnodes,SIZE(idseed))
    Nvertices(i) = idseed(i)
  ENDDO
  DEALLOCATE(idseed)
  ALLOCATE(idseed(SIZE(Nvertices)))
  idseed(:) = Nvertices(:)
  DEALLOCATE(Nvertices)
ELSE
  !Seeds are defined, check that seed IDs are within range
  DO i=1,SIZE(idseed)
    IF( idseed(i)>SIZE(Ps,1) ) THEN
      WRITE(msg,*) "/!\ WARNING: idseed(", i, ") out of bounds"
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      idseed(i) = SIZE(Ps,1)
    ENDIF
  ENDDO
ENDIF
!Check NPs(:): size must be = number of seeds, and values must be < max.Natoms for a seed
IF( .NOT.ALLOCATED(NPs) ) THEN
  msg = "/!\ WARNING: NPs not allocated, allocating it now"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !Number of atoms for each seed not defined: try to count them
  !For each seed, parse atomis number, the first =0 marks the end of the list
  ALLOCATE(NPs(SIZE(Ps,1)))
  NPs(:) = SIZE(Ps,2)
  DO i=1,SIZE(Ps,1)
    j = 1
    DO WHILE( j<=SIZE(Ps,2) )
      NPs(i) = j
      IF( Ps(i,j,4) > 1.d-3 ) EXIT
    ENDDO
  ENDDO
ELSEIF( SIZE(NPs).NE.SIZE(Ps,1) ) THEN
  WRITE(msg,*) "/!\ WARNING: NPs size does not match Ps", SIZE(NPs), SIZE(Ps,1)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !Natoms are defined but do not match number of seeds
  ALLOCATE(Nvertices(SIZE(Ps,1)))
  Nvertices(:) = SIZE(Ps,2)
  DO i=1,SIZE(Ps,1)
    j = 1
    DO WHILE( j<=SIZE(Ps,2) )
      Nvertices(i) = j
      IF( Ps(i,j,4) > 1.d-3 ) EXIT
    ENDDO
  ENDDO
  DEALLOCATE(NPs)
  ALLOCATE(NPs(SIZE(Nvertices)))
  NPs(:) = Nvertices(:)
  DEALLOCATE(Nvertices)
ELSE
  !Natoms are defined, check that they are within range
  DO i=1,SIZE(NPs)
    IF( NPs(i)<=0 ) THEN
      NPs(i) = 1
    ELSEIF( NPs(i)>SIZE(Ps,2) ) THEN
      NPs(i) = SIZE(Ps,2)
    ENDIF
  ENDDO
ENDIF
!
!For 2-D polycrystal, check that user rotated grains only around shortest axis
IF( twodim>0 ) THEN
  m=0
  DO i=1,SIZE(vorient,1)
    DO j=1,3
      IF( twodim==j .AND. DABS(DABS(vorient(i,j,j))-1.d0)>1.d-12 ) THEN
        m=1
      ENDIF
    ENDDO
  ENDDO
  IF( m>0 ) THEN
    nwarn = nwarn+1
    CALL ATOMSK_MSG(4715,(/''/),(/DBLE(twodim)/))
  ENDIF
ENDIF
!
IF(verbosity==4) THEN
  !Debug messages
  msg = "(((   Seed index for each node   )))"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(idseed,1)
    WRITE(msg,'(6X,i4,a6,i4)') i, "  ->  ", idseed(i)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
  msg = "(((  Atomic positions for seeds  )))"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(Ps,1)
    WRITE(msg,'(6X,a6,i4,a3,i9,a7)') "SEED #", i, "  (", NPs(i), " atoms)"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO j=1,3
      WRITE(msg,'(10X,3f9.3)') Hs(i,j,1), Hs(i,j,2), Hs(i,j,3)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
    DO j=1,NPs(i)
      WRITE(msg,'(6X,i4,6X,3f9.3)') NINT(Ps(i,j,4)), Ps(i,j,1), Ps(i,j,2), Ps(i,j,3)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
  ENDDO
  WRITE(msg,*) "(((  Number of nodes:", Nnodes
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  msg = "Positions of nodes:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(vnodes,1)
    WRITE(msg,'(i4,3f9.3)') i, vnodes(i,:)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
  msg = "Rotation matrices of nodes:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(vorient,1)
    WRITE(msg,'(a2,i4,a5,3f9.3,a1)') 'R[', i, '] = [', (vorient(i,1,j),j=1,3), ']'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,'(a2,i4,a5,3f9.3,a1)') 'R[', i, '] = [', (vorient(i,2,j),j=1,3), ']'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,'(a2,i4,a5,3f9.3,a1)') 'R[', i, '] = [', (vorient(i,3,j),j=1,3), ']'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
ENDIF
!
!
!
300 CONTINUE
!For each node, construct list of neighboring vertices (include/voronoi.f90)
CALL VORONEIGH(H,twodim,vnodes,Nvertices,vvertex,cell_volumes)
!
IF( nerr>0 ) THEN
  GOTO 1000
ENDIF
!
IF(verbosity==4) THEN
  !Debug messages
  WRITE(msg,*) "List of neighboring vertices for each node:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(vnodes,1)
    WRITE(msg,*) "    NODE #", i, ", vol = ", cell_volumes(i), " A^3, ", Nvertices(i), " vertices:"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO j=1,Nvertices(i)
      WRITE(msg,'(6X,i4,3f9.3,a7,f9.3)') j, vvertex(i,j,1:3), " | d = ", vvertex(i,j,4)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
  ENDDO
ENDIF
!
!For each node/grain, apply inverse rotation matrix of that grain to all vertices
rotmat(:,:) = Id_Matrix(:,:)
ALLOCATE(shiftgrain(Nnodes,1:3))
shiftgrain(:,:) = 0.d0
DO inode=1,Nnodes
  n = Nvertices(inode)
  IF( n>0 ) THEN
    !Sort neighboring vertices by increasing distance
    CALL BUBBLESORT(vvertex(inode,1:n,1:4),4,"up  ",newindex)
    !Compute the inverse of rotation matrix for this node
    CALL INVMAT( vorient(inode,1:3,1:3) , rotmat )
    !Rotate node position vector
    P1 = vnodes(inode,1)
    P2 = vnodes(inode,2)
    P3 = vnodes(inode,3)
    vnodes(inode,1) = P1*rotmat(1,1) + P2*rotmat(1,2) + P3*rotmat(1,3)
    vnodes(inode,2) = P1*rotmat(2,1) + P2*rotmat(2,2) + P3*rotmat(2,3)
    vnodes(inode,3) = P1*rotmat(3,1) + P2*rotmat(3,2) + P3*rotmat(3,3)
    !Rotate all vertices
    DO j=1,n
      P1 = vvertex(inode,j,1)
      P2 = vvertex(inode,j,2)
      P3 = vvertex(inode,j,3)
      vvertex(inode,j,1) = P1*rotmat(1,1) + P2*rotmat(1,2) + P3*rotmat(1,3)
      vvertex(inode,j,2) = P1*rotmat(2,1) + P2*rotmat(2,2) + P3*rotmat(2,3)
      vvertex(inode,j,3) = P1*rotmat(3,1) + P2*rotmat(3,2) + P3*rotmat(3,3)
    ENDDO
    !Find sphere that includes all vertices: vector(:)=center, P1=radius of sphere
    CALL BOUNDINGSPHERE(vvertex(inode,1:n,1:3),vector,P1,k)
    IF(twodim>0) THEN
      vector(twodim) = 0.5d0*H(twodim,twodim)
    ENDIF
    !Shift node and all vertices so they are inside [0,xmax], [0,ymax], [0,zmax]
    !NOTE: this shift vector will be used later to place atoms back into position
    IF( Nnodes<6 ) THEN
      P2 = 1.55d0*P1 + 20.d0
    ELSE
      P2 = 1.25d0*P1 + 12.d0
    ENDIF
    IF( twodim>0 ) THEN
      P2 = 1.5d0*P2 + 20.d0
    ENDIF
    DO i=1,3
      IF(i.NE.twodim) THEN
        shiftgrain(inode,i) = -1.d0*vector(i) + P2
      ENDIF
    ENDDO
    vnodes(inode,1:3) = vnodes(inode,1:3) + shiftgrain(inode,1:3)
    DO j=1,n
      vvertex(inode,j,1:3) = vvertex(inode,j,1:3) + shiftgrain(inode,1:3)
    ENDDO
  ENDIF
ENDDO
!
IF(verbosity>=4) THEN
  !Debug messages
  WRITE(msg,*) "List of neighboring vertices after rotation and shift:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO inode=1,SIZE(vnodes,1)
    WRITE(msg,*) "    NODE #", inode, ", vol = ", cell_volumes(inode), " A^3, ", Nvertices(inode), " vertices:"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO j=1,Nvertices(inode)
      WRITE(msg,'(6X,i4,3f9.3,a7,f9.3)') j, vvertex(inode,j,1:3), " | d = ", vvertex(inode,j,4)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
    !Write rotated positions of node and vertices into file
    k=40+inode
    WRITE(msg,*) inode
    OPEN(UNIT=k,FILE="atomsk_grain"//TRIM(ADJUSTL(msg))//"_vertices-rotated.xyz")
    WRITE(k,*) Nvertices(inode)+1
    WRITE(k,*) "# Position of node # "//TRIM(ADJUSTL(msg))//" and its vertices after rotation"
    WRITE(k,'(i4,6f12.3)') 2, vnodes(inode,:)
    DO i=1,Nvertices(inode)
      IF( vvertex(inode,i,4)>1.d-6 ) THEN
        WRITE(k,'(i4,3f12.3)') 1, vvertex(inode,i,1:3)
      ENDIF
    ENDDO
    CLOSE(k)
  ENDDO
ENDIF
!
!
!
500 CONTINUE
!Estimate final number of particles NP = (density of unit cell) * (volume of final cell)
!NOTE: final polycrystal is expected to count fewer atoms than a perfect crystal of same size,
!     however here we allocate for more atoms to avoid having to re-allocate arrays during the loop.
!     Actual final array sizes will be optimized later (see label 600)
!Compute seed density (=max. density for any given seed)
seed_density = 0.d0
DO i=1,SIZE(Hs,1)
  P1 = VOLUME_PARA(Hs(i,:,:))
  P2 = NPs(i) / P1
  IF( P2>seed_density ) THEN
    seed_density = P2
  ENDIF
ENDDO
WRITE(msg,*) " Seed density (max.) = ", seed_density
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
NP = CEILING( seed_density * VOLUME_PARA(H) )
NPexpect = NP
WRITE(msg,*) " Estimated max. number of atoms in this box: NP = ", NP
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!Allow for +10% or +1000 atoms to allocate arrays
NP = MIN( NINT(1.1d0*NP) , NP+1000 )
!Allocate arrays for final polycrystal
WRITE(msg,*) "                Allocating arrays with size: NP = ", NP
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!Array for atom positions
ALLOCATE(Q(NP,4) , STAT=i)
IF( i>0 ) THEN
  ! Allocation failed (not enough memory)
  nerr = nerr+1
  CALL ATOMSK_MSG(819,(/''/),(/DBLE(Nnodes)/))
  GOTO 1000
ENDIF
Q(:,:) = 0.d0
!Array for shells (if any)
IF(doshells) THEN
  ALLOCATE(T(NP,4) , STAT=i)
  IF( i>0 ) THEN
    ! Allocation failed (not enough memory)
    nerr = nerr+1
    CALL ATOMSK_MSG(819,(/''/),(/DBLE(Nnodes)/))
    GOTO 1000
  ENDIF
  T(:,:) = 0.d0
ENDIF
!Array for auxiliary properties (there is always at least one, grainID)
grainID = 0
seedID = 0
IF( ALLOCATED(AUXNAMES) ) THEN
  DO i=1,SIZE(AUXNAMES)
    IF( AUXNAMES(i)=="grainID" ) THEN
      grainID = i
    ELSEIF( AUXNAMES(i)=="seedID" ) THEN
      seedID = i
    ENDIF
  ENDDO
  n=0
  IF( grainID==0 .OR. seedID==0 ) THEN
    n = SIZE(AUXNAMES)
    IF( grainID==0 ) n=n+1
    IF( seedID==0  ) THEN
      !Define "seedID" only if several seeds or a matrix are defined
      IF( MAXVAL(idseed(:))>1 .OR. idmatrix>0 ) THEN
        n=n+1
      ENDIF
    ENDIF
  ENDIF
  IF( n>0 ) THEN
    m = SIZE(AUXNAMES)
    ALLOCATE(newAUXNAMES(n))
    newAUXNAMES(1:m) = AUXNAMES(1:m)
    DEALLOCATE(AUXNAMES)
    n = SIZE(newAUXNAMES)
    ALLOCATE(AUXNAMES(n))
    AUXNAMES(:) = newAUXNAMES(:)
    DEALLOCATE(newAUXNAMES)
    IF( MAXVAL(idseed(:))>1 .OR. idmatrix>0 ) THEN
      grainID = n-1
      seedID = n
    ELSE
      grainID = n
    ENDIF
  ENDIF
ELSE
  IF( MAXVAL(idseed(:))>1 .OR. idmatrix>0 ) THEN
    ALLOCATE(AUXNAMES(2))
    seedID = 2
  ELSE
    ALLOCATE(AUXNAMES(1))
  ENDIF
  AUXNAMES(:) = ""
  grainID = 1
ENDIF
AUXNAMES(grainID) = "grainID"
IF(seedID>0) AUXNAMES(seedID) = "seedID"
ALLOCATE(newAUX(NP,SIZE(AUXNAMES)) , STAT=i)
IF( i>0 ) THEN
  ! Allocation failed (not enough memory)
  nerr = nerr+1
  CALL ATOMSK_MSG(819,(/''/),(/DBLE(Nnodes)/))
  GOTO 1000
ENDIF
newAUX(:,:) = 0.d0
!
!Allocate array containing number of atoms in each grain
IF( idmatrix>0 ) THEN
  ALLOCATE( NPgrains(Nnodes+1) , STAT=i )
ELSE
  ALLOCATE( NPgrains(Nnodes) , STAT=i )
ENDIF
IF( i>0 ) THEN
  ! Allocation failed (not enough memory)
  nerr = nerr+1
  CALL ATOMSK_MSG(819,(/''/),(/DBLE(Nnodes)/))
  GOTO 1000
ENDIF
!
!
!Construct grains using Voronoi tesselation
!For each node, list of neighboring vertices was established above
NPgrains(:) = 0 !count atoms in each grain
NP=0            !count all atoms that will make it in the final polycrystal
rotmat(:,:) = Id_Matrix(:,:)
DO inode=1,Nnodes
  !This is grain # inode
  CALL ATOMSK_MSG(4055,(/''/),(/DBLE(inode),DBLE(Nvertices(inode))/))
  WRITE(msg,*) "  _  _  _  _  _  _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,'(a,i4)') " |             GRAIN #", inode
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,'(a,i4,a3,i9,a7)') " |        USING SEED #", idseed(inode), "  (", NPs(idseed(inode)), " atoms)"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Don't do the matrix now, it will be dealt with later
  !IF( idseed(inode) == idmatrix ) THEN
  !  GOTO 540
  !ENDIF
  !
  !Get number of neighboring vertices
  n = Nvertices(inode)
  !
  !Huc will be that seed's periodic vectors
  Huc(:,:) = Hs(idseed(inode),:,:)
  !
  !Compute current seed density
  seed_density = NPs(idseed(inode)) / VOLUME_PARA(Huc(:,:))
  !
  !Estimate duplication factors
  Ht(:,:) = Huc(:,:)
  DO i=1,3
    distance = 0.d0
    n = Nvertices(inode)
    P1 = MAX( 0.d0 , MAXVAL(vvertex(inode,1:n,i)) ) - MIN( 0.d0 , MINVAL(vvertex(inode,1:n,i)) )
    IF( P1>distance ) THEN
      distance = MAX(P1,20.d0)
    ENDIF
    IF( Nnodes<6 ) THEN
      Ht(i,i) = 1.55d0*distance + 20.d0
    ELSE
      Ht(i,i) = 1.25d0*distance + 12.d0
    ENDIF
  ENDDO
  IF( twodim>0 ) THEN
    !System is pseudo-2D: template won't be duplicated along the short dimension
    !We can afford to make it larger along other directions
    Ht(:,:) = 2.d0*Ht(:,:) + 20.d0
    Ht(twodim,twodim) = VECLENGTH(H(twodim,:))
  ENDIF
  !Estimate how many particles the template will contain = (seed density)*(template volume)
  P1 = 1.05d0 * seed_density * Ht(1,1)*Ht(2,2)*Ht(3,3)
  !
  !Compute duplication factors along X, Y, Z
  !If the system is 2-D, do not expand along the shortest axis
  !NOTE: these duplication factors will be used to duplicate seed to fill the grain
  expandmatrix(:) = 0
  WRITE(msg,*) "Determining expansion factors:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,3  !Loop on X, Y, Z
    IF( i.NE.twodim ) THEN
      !Round up to nearest greater integer
      H2 = MAXVAL(Huc(:,i)) + MINVAL(Huc(:,i))
      expandmatrix(i) = MAX(1 , CEILING( Ht(i,i)/H2 ) )
    ENDIF
  ENDDO
  IF(twodim>0) THEN
    expandmatrix(twodim) = 0
  ENDIF
  WRITE(msg,*) "Initial expansion factors:", expandmatrix(:)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Seed will be duplicated inside the grain
  istart = NP+1
  !
  !Duplicate seed to fill the grain
  !Loop on all atoms in seed/unit cell and their periodic images
  !$OMP PARALLEL DO DEFAULT(SHARED) &
  !$OMP& PRIVATE(i,j,m,n,o,P1,P2,P3,keepatom,vector,vnormal)
  DO o=0,expandmatrix(3)
    DO n=0,expandmatrix(2)
      DO m=0,expandmatrix(1)
        DO i=1,NPs(idseed(inode))
          !Compute (cartesian) position of the replica of this atom
          P1 = Ps(idseed(inode),i,1) + DBLE(m)*Huc(1,1) + DBLE(n)*Huc(2,1) + DBLE(o)*Huc(3,1)
          P2 = Ps(idseed(inode),i,2) + DBLE(m)*Huc(1,2) + DBLE(n)*Huc(2,2) + DBLE(o)*Huc(3,2)
          P3 = Ps(idseed(inode),i,3) + DBLE(m)*Huc(1,3) + DBLE(n)*Huc(2,3) + DBLE(o)*Huc(3,3)
          vector(:) = (/P1,P2,P3/) - vnodes(inode,1:3)
          !Check if atom is inside the grain #inode
          j = 1
          keepatom = .TRUE.
          !Atoms of the seed that are *outside* polyhedra must be removed
          DO WHILE( keepatom .AND. j <= Nvertices(inode) )
            vnormal(:) = vvertex(inode,j,1:3) - vnodes(inode,1:3)
            IF( VEC_PLANE(vnormal,VECLENGTH(vnormal),vector) > -1.d0*clearance ) THEN
              !Atom is above this plane of cut, hence outside of polyhedron
              keepatom = .FALSE.
              EXIT
            ENDIF
            j=j+1
          ENDDO
          !
          IF(keepatom) THEN
            !$OMP CRITICAL
            !Atom #i survived all cuts => save its position and species into Q
            NPgrains(inode) = NPgrains(inode)+1
            NP = NP+1
            IF( NP > SIZE(Q,1) ) THEN
              !Increase array size by 10%
              PRINT*, "     --->   RESIZE Q   <---"
              NPnew = CEILING(1.1d0*SIZE(Q,1))
              CALL RESIZE_DBLEARRAY2(Q,NPnew,4,status)
              IF(doshells) THEN
                CALL RESIZE_DBLEARRAY2(T,NPnew,4,status)
              ENDIF
              IF(doaux) THEN
                CALL RESIZE_DBLEARRAY2(newAUX,NPnew,SIZE(newAUX,2),status)
              ENDIF
            ENDIF
            !Shift atom position and apply rotation
            P1 = P1 - shiftgrain(inode,1)
            P2 = P2 - shiftgrain(inode,2)
            P3 = P3 - shiftgrain(inode,3)
            Q(NP,1) = P1*vorient(inode,1,1) + P2*vorient(inode,1,2) + P3*vorient(inode,1,3)
            Q(NP,2) = P1*vorient(inode,2,1) + P2*vorient(inode,2,2) + P3*vorient(inode,2,3)
            Q(NP,3) = P1*vorient(inode,3,1) + P2*vorient(inode,3,2) + P3*vorient(inode,3,3)
            Q(NP,4) = Ps(idseed(inode),i,4)
            IF(doshells) THEN
              T(NP,4) = Ss(idseed(inode),i,4)
              IF( T(NP,4)>0.1d0 ) THEN
                P1 = Ss(idseed(inode),i,1) - shiftgrain(inode,1)
                P2 = Ss(idseed(inode),i,2) - shiftgrain(inode,2)
                P3 = Ss(idseed(inode),i,3) - shiftgrain(inode,3)
                T(NP,1) = P1*vorient(inode,1,1) + P2*vorient(inode,1,2) + P3*vorient(inode,1,3)
                T(NP,2) = P1*vorient(inode,2,1) + P2*vorient(inode,2,2) + P3*vorient(inode,2,3)
                T(NP,3) = P1*vorient(inode,3,1) + P2*vorient(inode,3,2) + P3*vorient(inode,3,3)
              ENDIF
            ENDIF
            IF(doaux) THEN
              !Copy auxiliary properties of atom #i in newAUX
              newAUX(NP,:) = AUXs(idseed(inode),i,:)
            ENDIF
            !Save grainID and seedID for this atom
            newAUX(NP,grainID) = DBLE(inode)
            IF(seedID>0) newAUX(NP,seedID) = idseed(inode)
            !$OMP END CRITICAL
          ENDIF
        ENDDO !loop on i
      ENDDO !loop on m
    ENDDO !loop on n
  ENDDO !loop on o
  !$OMP END PARALLEL DO
  !
  IF(nerr>0) GOTO 1000
  !
  CALL ATOMSK_MSG(4056,(/''/),(/DBLE(NPgrains(inode)),cell_volumes(inode)/))
  !
  IF(NPgrains(inode)==0) THEN
    !There are zero atoms in this grain, which is strange
    !Display a warning message
    nwarn=nwarn+1
    CALL ATOMSK_MSG(4708,(/""/),(/0.d0/))
  ENDIF
  !
  IF( verbosity==4 ) THEN
    !Debug: write positions of atoms of current grain into an XYZ file
    WRITE(temp,*) inode
    msg = 'atomsk_grain'//TRIM(ADJUSTL(temp))//'.xyz'
    OPEN(UNIT=36,FILE=msg,STATUS="UNKNOWN",FORM="FORMATTED")
    WRITE(36,*) NPgrains(inode)
    msg = '# Debug file for Atomsk, mode --polycrystal, grain # '//TRIM(ADJUSTL(temp))
    WRITE(36,*) TRIM(msg)
    IF(NPgrains(inode)>0) THEN
      DO i=istart,NP
        CALL ATOMSPECIES(Q(i,4),species)
        WRITE(36,'(a2,2X,3(f16.8,1X))') species, Q(i,1:3)
      ENDDO
    ENDIF
!     WRITE(36,'(a4)') 'alat'
!     WRITE(36,'(a3)') '1.0'
!     WRITE(36,'(a9)') 'supercell'
!     WRITE(36,'(3f16.6)') H(1,:)
!     WRITE(36,'(3f16.6)') H(2,:)
!     WRITE(36,'(3f16.6)') H(3,:)
    CLOSE(36)
    !
    WRITE(msg,*) " |_  _  _  _  _  END OF GRAIN #", inode, "  _ _ _ _ _ _ _ _ _ _|"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
  !
  540 CONTINUE
  !
ENDDO  !end loop on inode
!
!
550 CONTINUE
!Now that all grains are filled, deal with the matrix (if any)
IF( idmatrix>0 .AND. idmatrix<=SIZE(Ps,1) ) THEN
  WRITE(msg,*) "  _  _  _  _  _  _  _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "  |      NOW FILLING SPACE BETWEEN GRAINS WITH INTEGRANULAR PHASE:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  CALL ATOMSK_MSG(4082,(/''/),(/0.d0/))
  !
  !Restore the original node positions and orientations
  DO inode=1,Nnodes
    n = Nvertices(inode)
    IF( n>0 ) THEN
      !Define rotation matrix for node #inode
      rotmat(:,:) = vorient(inode,1:3,1:3)
      !Shift and then rotate node #inode
      P1 = vnodes(inode,1) - shiftgrain(inode,1)
      P2 = vnodes(inode,2) - shiftgrain(inode,2)
      P3 = vnodes(inode,3) - shiftgrain(inode,3)
      vnodes(inode,1) = P1*rotmat(1,1) + P2*rotmat(1,2) + P3*rotmat(1,3)
      vnodes(inode,2) = P1*rotmat(2,1) + P2*rotmat(2,2) + P3*rotmat(2,3)
      vnodes(inode,3) = P1*rotmat(3,1) + P2*rotmat(3,2) + P3*rotmat(3,3)
      !Shift and then rotate vertices of node #inode
      DO j=1,n
        P1 = vvertex(inode,j,1) - shiftgrain(inode,1)
        P2 = vvertex(inode,j,2) - shiftgrain(inode,2)
        P3 = vvertex(inode,j,3) - shiftgrain(inode,3)
        vvertex(inode,j,1) = P1*rotmat(1,1) + P2*rotmat(1,2) + P3*rotmat(1,3)
        vvertex(inode,j,2) = P1*rotmat(2,1) + P2*rotmat(2,2) + P3*rotmat(2,3)
        vvertex(inode,j,3) = P1*rotmat(3,1) + P2*rotmat(3,2) + P3*rotmat(3,3)
      ENDDO
    ENDIF
  ENDDO
  !
  ALLOCATE(ingrain(Nnodes*27))
  ingrain(:) = .TRUE.
  !
  !Compute duplication factors along X, Y, Z
  !If the system is 2-D, do not expand along the shortest axis
  !NOTE: these duplication factors will be used to duplicate seed to fill the grain
  expandmatrix(:) = 0
  Huc(:,:) = Hs(idmatrix,:,:)
  WRITE(msg,*) "Determining expansion factors:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,3  !Loop on X, Y, Z
    IF( i.NE.twodim ) THEN
      !Round up to nearest greater integer
      H2 = MAXVAL(Huc(:,i)) + MINVAL(Huc(:,i))
      expandmatrix(i) = MAX(1 , CEILING( H(i,i)/H2 ) )
    ENDIF
  ENDDO
  IF(twodim>0) THEN
    expandmatrix(twodim) = 0
  ENDIF
  WRITE(msg,*) "Initial expansion factors:", expandmatrix(:)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  IF( NPs(idmatrix)<6 ) THEN
    !Resize matrix cell to match supercell size
    DO i=1,3
      IF( expandmatrix(i)>10 ) THEN
        !New cell size = supercell size / integer
        H1 = H(i,i) / expandmatrix(i)
        !Update all atom positions
        DO j=1,NPs(idmatrix)
          Ps(idmatrix,j,i) = Ps(idmatrix,j,i) * H1 / Huc(i,i) + 0.01d0
        ENDDO
        !Update cell size
        Huc(i,i) = H1
      ENDIF
    ENDDO
  ENDIF
  !
  !
  !Duplicate matrix to fill the whole box
  !Atoms that are inside a grain are removed
  !Loop on all atoms in seed/unit cell
  !$OMP PARALLEL DO DEFAULT(SHARED) &
  !$OMP& PRIVATE(inode,i,j,k,i2,j2,k2,m,n,o,H1,H2,H3,P1,P2,P3,ingrain,keepatom,vector,vnormal)
  DO o=0,expandmatrix(3)
    DO n=0,expandmatrix(2)
      DO m=0,expandmatrix(1)
        DO i=1,NPs(idmatrix)  !Loop on all atoms in the matrix
          !Compute (cartesian) position of the replica of this atom
          P1 = Ps(idmatrix,i,1) + DBLE(m)*Huc(1,1) + DBLE(n)*Huc(2,1) + DBLE(o)*Huc(3,1)
          P2 = Ps(idmatrix,i,2) + DBLE(m)*Huc(1,2) + DBLE(n)*Huc(2,2) + DBLE(o)*Huc(3,2)
          P3 = Ps(idmatrix,i,3) + DBLE(m)*Huc(1,3) + DBLE(n)*Huc(2,3) + DBLE(o)*Huc(3,3)
          !
          !Only consider atoms that are inside the final supercell
          IF( P1>=0.d0 .AND. P1<H(1,1) .AND. &
            & P2>=0.d0 .AND. P2<H(2,2) .AND. &
            & P3>=0.d0 .AND. P3<H(3,3)       ) THEN
            !Check if atom #i is inside a grain, if it's not in any grain then we keep it
            ingrain(:) = .FALSE.
            !Loop on all grains and their periodic replicas
            k = 1
            DO inode=1,Nnodes
              DO k2=-1,1,1
                DO j2=-1,1,1
                  DO i2=-1,1,1
                    ingrain(k) = .TRUE.
                    H1 = DBLE(i2)*H(1,1)
                    H2 = DBLE(j2)*H(2,2)
                    H3 = DBLE(k2)*H(3,3)
                    vector(:) = (/P1,P2,P3/) - vnodes(inode,1:3) - (/H1,H2,H3/)
                    !Loop on all vertices of that grain
                    j = 1
                    DO WHILE( ingrain(k) .AND. j <= Nvertices(inode) )
                      vnormal(:) = vvertex(inode,j,1:3) - vnodes(inode,1:3)
                      IF( VEC_PLANE(vnormal,VECLENGTH(vnormal),vector) > -1.d0*clearance+0.2d0 ) THEN
                        !Atom #i is above this plane of cut, hence outside grain #inode
                        ingrain(k) = .FALSE.
                        EXIT
                      ENDIF
                      j=j+1
                    ENDDO
                    IF(ingrain(k)) EXIT
                    k=k+1
                  ENDDO !loop on i2
                  IF(ingrain(k)) EXIT
                ENDDO !loop on j2
                IF(ingrain(k)) EXIT
              ENDDO !loop on k2
              IF(ingrain(k)) EXIT
            ENDDO !loop on inode
            !
            !If any ingrain(:) is TRUE, it means atom #i is inside one of the grains: discard it
            IF( ANY(ingrain(:)) ) THEN
              keepatom = .FALSE.
            ELSE
              !Otherwise, atom #i is outside of all grains: keep it
              keepatom = .TRUE.
            ENDIF
            !
            IF(keepatom) THEN
              !$OMP CRITICAL
              !Atom #i survived all cuts => save its position and species into Q
              NPgrains(Nnodes+1) = NPgrains(Nnodes+1)+1
              NP = NP+1
              IF( NP > SIZE(Q,1) ) THEN
                !Increase array size by 10%
                PRINT*, "     --->   RESIZE Q   <---"
                NPnew = CEILING(1.1d0*SIZE(Q,1))
                CALL RESIZE_DBLEARRAY2(Q,NPnew,4,status)
                IF(doshells) THEN
                  CALL RESIZE_DBLEARRAY2(T,NPnew,4,status)
                ENDIF
                IF(doaux) THEN
                  CALL RESIZE_DBLEARRAY2(newAUX,NPnew,SIZE(newAUX,2),status)
                ENDIF
              ENDIF
              !Save final atom position
              Q(NP,1) = P1
              Q(NP,2) = P2
              Q(NP,3) = P3
              Q(NP,4) = Ps(idmatrix,i,4)
              IF(doshells) THEN
                T(NP,4) = Ss(idmatrix,i,4)
                IF( T(NP,4)>0.1d0 ) THEN
                  T(NP,1) = Ss(idmatrix,i,1)
                  T(NP,2) = Ss(idmatrix,i,2)
                  T(NP,3) = Ss(idmatrix,i,3)
                ENDIF
              ENDIF
              IF(doaux) THEN
                newAUX(NP,:) = AUXs(idmatrix,i,:)
              ENDIF
              newAUX(NP,grainID) = Nnodes+1
              IF(seedID>0) newAUX(NP,seedID) = MAXVAL(idseed)+1
              !$OMP END CRITICAL
            ENDIF
            !
          ENDIF
          !
        ENDDO !loop on i
      ENDDO !loop on m
    ENDDO !loop on n
  ENDDO !loop on o
  !$OMP END PARALLEL DO
  !
  IF(ALLOCATED(ingrain)) DEALLOCATE(ingrain)
  !
  CALL ATOMSK_MSG(4056,(/''/),(/DBLE(NPgrains(Nnodes+1)),0.d0/))
  !
  WRITE(msg,*) " |_  _  _  _  _  _   END OF MATRIX   _ _ _ _ _ _ _ _ _ _ _|"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!
!
!
600 CONTINUE
WRITE(msg,*) "Loop on grains finished, size of Q = ", SIZE(Q,1)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
IF(ALLOCATED(Hs)) DEALLOCATE(Hs)
IF(ALLOCATED(Ps)) DEALLOCATE(Ps)
IF(ALLOCATED(Ss)) DEALLOCATE(Ss)
IF(ALLOCATED(Nvertices)) DEALLOCATE(Nvertices)
IF(ALLOCATED(vnodes)) DEALLOCATE(vnodes)
IF(ALLOCATED(vvertex)) DEALLOCATE(vvertex)
IF(ALLOCATED(shiftgrain)) DEALLOCATE(shiftgrain)
!
!NP is now the actual number of atoms in the final polycrystal
WRITE(msg,*) "Now allocating arrays with NP = ", NP
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
IF( NP .NE. SUM(NPgrains(:)) ) THEN
  PRINT*, "ERROR while counting atoms: ", NP, " vs ", SUM(NPgrains(:))
  nerr = nerr+1
  GOTO 1000
ENDIF
!Q now contains positions of all atoms in all the grains, but may be oversized
!(and T contains the positions of shells, and newAUX the aux.prop. if relevant)
!Copy atom positions into final array P (and S and AUX if needed) with appropriate size
IF(ALLOCATED(P)) DEALLOCATE(P)
ALLOCATE( P(NP,4) , STAT=n )
IF( n>0 ) THEN
  ! Allocation failed (not enough memory)
  nerr = nerr+1
  CALL ATOMSK_MSG(819,(/''/),(/DBLE(NP)/))
  GOTO 1000
ENDIF
P(1:NP,:) = Q(1:NP,:)
IF(ALLOCATED(Q)) DEALLOCATE(Q)
!Copy shell positions into final array S
IF(ALLOCATED(S)) DEALLOCATE(S)
IF(doshells) THEN
  ALLOCATE( S(NP,4) , STAT=n )
  IF( n>0 ) THEN
    ! Allocation failed (not enough memory)
    nerr = nerr+1
    CALL ATOMSK_MSG(819,(/''/),(/DBLE(NP)/))
    GOTO 1000
  ENDIF
  S(1:NP,:) = T(1:NP,:)
  IF(ALLOCATED(T)) DEALLOCATE(T)
ENDIF
!Copy aux.prop. into final array AUX
!NOTE: AUX contains at least the grainID
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
ALLOCATE( AUX(NP,SIZE(AUXNAMES)) , STAT=n )
IF( n>0 ) THEN
  ! Allocation failed (not enough memory)
  nerr = nerr+1
  CALL ATOMSK_MSG(819,(/''/),(/DBLE(NP*SIZE(AUXNAMES))/))
  GOTO 1000
ENDIF
AUX(1:NP,:) = newAUX(1:NP,:)
IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
!
!
!
700 CONTINUE
!From now on, some arrays are not needed anymore => deallocate them
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(S)) DEALLOCATE(S)
IF(ALLOCATED(T)) DEALLOCATE(T)
!
!
!Find the max. volume occupied by a grain
Vmin = MINVAL(cell_volumes)
Vmax = MAXVAL(cell_volumes)
m=0
ALLOCATE(Q(Nnodes,4))
DO i=1,Nnodes
  !Determine the position of the center of mass of this grain
  P1 = 0.d0
  vector(:) = 0.d0 !position of center of mass
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,m,species) REDUCTION(+:P1,vector)
  DO j=1,NPgrains(i)
    m = m+1
    CALL ATOMSPECIES(P(m,4),species)
    CALL ATOMMASS(species,P2)
    vector(:) = vector(:) + P2*P(m,1:3)
    P1 = P1 + P2
  ENDDO
  !$OMP END PARALLEL DO
  vector(:) = vector(:) / P1
  Q(i,1:3) = vector(:)
  Q(i,4) = 1.d0  !positions of com are given the atomic number of hydrogen
ENDDO
!
IF( ofu.NE.6 ) THEN
  !Write positions of center of mass of each grain into a file
  ALLOCATE(comment(1))
  comment(1) = "# Positions of grains centers of mass"
  ALLOCATE(outfileformats(1))
  outfileformats(1) = "xsf"
  msg = TRIM(TRIM_EXT(prefix))//"_grains-com.xsf"
  CALL WRITE_AFF(msg,outfileformats,H,Q,S,comment,newAUXNAMES,newAUX)
  IF(ALLOCATED(Q)) DEALLOCATE(Q)
  IF(ALLOCATED(comment)) DEALLOCATE(comment)
  IF(ALLOCATED(outfileformats)) DEALLOCATE(outfileformats)
  !
  WRITE(msg,*) "Min, max. volume occupied by a grain: ", Vmin, Vmax
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  IF( Vmax > 0.d0 ) THEN
    !Write grain ID and their sizes into a file
    IF(.NOT.overw) CALL CHECKFILE(idsizefile,'writ')
    OPEN(UNIT=41,FILE=idsizefile,FORM="FORMATTED",STATUS="UNKNOWN")
    IF( twodim>0 ) THEN
      WRITE(41,*) "# grainID   ;   N.atoms   ;   Grain area (A^2)"
    ELSE
      WRITE(41,*) "# grainID   ;   N.atoms   ;   Grain volume (A^3)"
    ENDIF
    DO i=1,Nnodes  !loop on all grains
      WRITE(41,*) i, NPgrains(i), cell_volumes(i)
    ENDDO
    CLOSE(41)
    !
    !Compute the grain size distribution and write it to a file
    !with the format:  ( Grain size ; Number of grains with that size )
    Vstep = (Vmax-Vmin)/20.d0  !step for grain size distribution
    IF(.NOT.overw) CALL CHECKFILE(distfile,'writ')
    OPEN(UNIT=41,FILE=distfile,FORM="FORMATTED",STATUS="UNKNOWN")
    IF( twodim>0 ) THEN
      WRITE(41,*) "# Grain size distribution: area (A^2)   ; No. of grains"
    ELSE
      WRITE(41,*) "# Grain size distribution: volume (A^3) ; No. of grains"
    ENDIF
    DO j=0,20  !loop on grain size
      k = 0
      DO i=1,Nnodes !loop on all grains
        IF( cell_volumes(i) >= Vmin+DBLE(j)*Vstep .AND. &
          & cell_volumes(i) < Vmin+DBLE(j+1)*Vstep      ) THEN
          !This grain has the appropriate size => increment counter
          k = k+1
        ENDIF
      ENDDO
      !Write to file
      WRITE(41,*) Vmin+DBLE(j)*Vstep, k
    ENDDO
    CLOSE(41)
    msg = "DATA"
    CALL ATOMSK_MSG(3002,(/distfile,msg/),(/0.d0/))
  ENDIF
ENDIF
GOTO 1000
!
!
!
800 CONTINUE
nerr=nerr+1
GOTO 1000
!
810 CONTINUE
nerr=nerr+1
CALL ATOMSK_MSG(2813,(/TRIM(temp)/),(/0.d0/))
GOTO 1000
!
820 CONTINUE
nerr=nerr+1
CALL ATOMSK_MSG(4820,(/TRIM(temp)/),(/0.d0/))
GOTO 1000
!
830 CONTINUE
nerr=nerr+1
CALL ATOMSK_MSG(1801,(/TRIM(" ")/),(/DBLE(linenumber)/))
GOTO 1000
!
!
!
1000 CONTINUE
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(S)) DEALLOCATE(S)
IF(ALLOCATED(T)) DEALLOCATE(T)
IF(ALLOCATED(Nvertices)) DEALLOCATE(Nvertices)
IF(ALLOCATED(NPgrains)) DEALLOCATE(NPgrains)
!
!
!
END SUBROUTINE POLYCRYS
!
!
!
END MODULE mode_polycrystal
