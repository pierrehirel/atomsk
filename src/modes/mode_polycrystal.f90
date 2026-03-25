MODULE mode_polycrystal
!
!**********************************************************************************
!*  MODE_POLYCRYSTAL                                                              *
!**********************************************************************************
!* This module constructs a polycrystal using a Voronoi tessellation, given       *
!* a unit cell defined by vectors H and atom positions P, and given the           *
!* position of the node and lattice orientation for each grain.                   *
!**********************************************************************************
!* (C) May 2013 - Pierre Hirel                                                    *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 19 Feb. 2026                                     *
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
USE readin
USE options
USE center
USE orient
USE writeout
USE polyx_readparam
USE voronoi
!
CONTAINS
!
!
SUBROUTINE POLYCRYS(ucfile,vfile,options_array,prefix,outfileformats,wof,H,P)
!
!
IMPLICIT NONE
!Input parameters
CHARACTER(LEN=*),INTENT(IN):: ucfile  !name of file containing seed (usually a unit cell, but can be anything)
CHARACTER(LEN=*),INTENT(IN):: vfile   !name of file containing parameters for Voronoi construction
CHARACTER(LEN=*),INTENT(IN):: prefix  !name or prefix for output file (polycrystal)
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to output
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
LOGICAL,INTENT(IN):: wof !write output file?
!
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=4096):: outparamfile  !file where grain parameters are written (if some parameters equal "random")
CHARACTER(LEN=4096):: distfile, idsizefile !name of file containing grain size distribution, grain sizes
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES    !names of auxiliary properties of atoms
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: newAUXNAMES !names of auxiliary properties of atoms (temporary)
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: doshells, doaux !are there shells, auxiliary properties in initial seed?
LOGICAL:: Hset            !are the box vectors H(:,:) defined?
LOGICAL:: exceeds100      !does the number of neighbours exceed 100?
LOGICAL:: isinpolyhedron  !is atom inside the polyhedron?
LOGICAL:: protectuc       !protect unit cell integrity?
LOGICAL:: use_template    !use a template or no?
LOGICAL,DIMENSION(4):: outparam  !Are keywords (1) "node", (2) "lattice", (3) "random" used?
                                 !(4) Must parameters be saved in a text file?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT
LOGICAL,DIMENSION(:),ALLOCATABLE:: Ptmask  !mask for template
INTEGER,DIMENSION(:),ALLOCATABLE:: newindex !list of sorted indexes
INTEGER:: twodim        !=0 if system is 3-D, =1,2,3 if system is thin along x, y, z
INTEGER:: grainID       !position of the auxiliary property "grainID" in AUX
INTEGER:: i, j, k
INTEGER:: istart
INTEGER:: linenumber    !line number when reading a file
INTEGER:: m, n, o
INTEGER:: maxvertex   !max. number of vertices to look for (defined for 3-D or 2-D)
INTEGER:: NPexpect   !expected total number of atoms in the final system
INTEGER:: NP         !total number of atoms in the final system
INTEGER:: qi   !used to count atoms in a grain
INTEGER:: inode, jnode
INTEGER:: Nnodes      !number of nodes
INTEGER:: status
INTEGER,DIMENSION(3):: expandmatrix
INTEGER,DIMENSION(3):: Nsurv  !number of atoms that survived a loop
INTEGER,DIMENSION(:),ALLOCATABLE:: Nvertices !number of vertices for each node
INTEGER,DIMENSION(:),ALLOCATABLE:: NPgrains  !number of atoms in each grain
INTEGER,DIMENSION(:,:,:),ALLOCATABLE:: NeighListMask  !mask for neighboring vertices
REAL(dp):: distance    !distance between two points
REAL(dp):: clearance   !clear atoms that close to the GB
REAL(dp):: maxdvertices !maximum distance between 2 vertices
REAL(dp):: H1, H2, H3  !temporary position
REAL(dp):: P1, P2, P3  !temporary position
REAL(dp):: seed_density !density of the seed (N.atoms/Volume)
REAL(dp):: Volume, Vmin, Vmax, Vstep  !min, max. volume occupied by a grain, step for grain size distribution
REAL(dp),DIMENSION(3):: GrainCenter !position of center of grain (differs from node position)
REAL(dp),DIMENSION(3):: shift       !shift vector
REAL(dp),DIMENSION(3):: vector      !vector between an atom and a node
REAL(dp),DIMENSION(3):: vnormal     !vector normal to grain boundary
REAL(dp),DIMENSION(3,3):: Huc       !Base vectors of the unit cell (seed)
REAL(dp),DIMENSION(3,3):: Ht        !Base vectors of the oriented unit cell or template
REAL(dp),DIMENSION(3,3):: H         !Base vectors of the final supercell
REAL(dp),DIMENSION(3,3):: ORIENT    !crystalographic orientation
REAL(dp),DIMENSION(3,3):: rotmat    !rotation matrix
REAL(dp),DIMENSION(4,3):: verp      !positions (x,y,z) of some vertices
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:),ALLOCATABLE:: cell_volumes !volume of each Voronoi cell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Puc, Suc  !positions of atoms, shells in unit cell (seed)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Pt, St    !positions of atoms, shells in template supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, S      !positions of atoms, shells in final supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Q, T      !positions of atoms, shells in a grain
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUXuc     !auxiliary properties of atoms in the unit cell (seed)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX_Q     !auxiliary properties of atoms in a grain
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX       !auxiliary properties of atoms in the final supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX    !auxiliary properties of atoms (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: shiftgrain !translation vectors for each grain
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: vnodes    !cartesian coordinate of each node
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: PosList   !positions of nearest neighbours
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE:: vvertex !for each node, position of neighboring vertices
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE:: vorient !rotation matrix for each node
!
!
!Initialize variables
CALL NAME_OUTFILE(prefix,temp,"     ")
i=SCAN(temp,".",BACK=.TRUE.)  !remove trailing dot
temp = temp(1:i-1)
distfile = TRIM(ADJUSTL(temp))//"_size-dist.txt"
idsizefile = TRIM(ADJUSTL(temp))//"_id-size.txt"
outparamfile = TRIM(ADJUSTL(temp))//"_param.txt"
outparam(:) = .FALSE.
protectuc = .FALSE.
use_template = .TRUE. !By default use a template, unless its size would be too large (see below)
m=0
Nnodes = 0
twodim = 0    !assume system will be 3-D
expandmatrix(:) = 0
clearance = 0.1d0  !atoms this close to the GB plane will be deleted
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
IF(ALLOCATED(NPgrains)) DEALLOCATE(NPgrains)
 C_tensor(:,:) = 0.d0
IF(ALLOCATED(P)) DEALLOCATE(P)
IF(ALLOCATED(shiftgrain)) DEALLOCATE(shiftgrain)
IF(ALLOCATED(vnodes)) DEALLOCATE(vnodes)
IF(ALLOCATED(vorient)) DEALLOCATE(vorient)
!
!
CALL ATOMSK_MSG(4054,(/''/),(/0.d0/))
!
!
100 CONTINUE
!Read initial seed from file
!NOTE: usually when constructing a polycrystal the seed is a unit cell,
!     e.g. a unit cell of fcc or bcc crystal. Here, no such assumption is made,
!     and the seed can be anything: a unit cell, a supercell, a large system containing
!     defects, dislocations or whatsoever. Also, the seed may be smaller or larger
!     than the final polycrystal, i.e. Huc(:,:) may be smaller than H(:,:), but
!     it might also be larger. Keep that in mind if you modify this routine
CALL READ_AFF(ucfile,Huc,Puc,Suc,comment,AUXNAMES,AUXuc)
!
IF(nerr>0 ) GOTO 1000
IF( .NOT.ALLOCATED(Puc) .OR. SIZE(Puc,1)<1 ) THEN
  CALL ATOMSK_MSG(804,(/''/),(/0.d0/))
  nerr=nerr+1
  GOTO 1000
ENDIF
!
!Check if seed contains shells (in the sense of core-shell model) and/or auxiliary properties
IF( ALLOCATED(Suc) .AND. SIZE(Suc,1)>0 ) THEN
  doshells = .TRUE.
ELSE
  doshells = .FALSE.
ENDIF
IF( ALLOCATED(AUXuc) .AND. SIZE(AUXuc,1)>0 ) THEN
  doaux = .TRUE.
  grainID = SIZE(AUXNAMES)+1
  ALLOCATE( newAUXNAMES(grainID) )
  DO i=1,SIZE(AUXNAMES)
    newAUXNAMES(i) = AUXNAMES(i)
  ENDDO
  DEALLOCATE(AUXNAMES)
  ALLOCATE( AUXNAMES(grainID) )
  AUXNAMES(:) = newAUXNAMES(:)
  DEALLOCATE(newAUXNAMES)
ELSE
  doaux = .FALSE.
  grainID = 1
  ALLOCATE(AUXNAMES(1))
ENDIF
AUXNAMES(grainID) = "grainID"
!
!Compute seed density
P1 = VOLUME_PARA(Huc)
seed_density = SIZE(Puc,1) / P1
!
IF( ALLOCATED(comment) ) DEALLOCATE(comment)
!
!
!
200 CONTINUE
!Read parameter file (modes/polyx_readparam.f90)
CALL ATOMSK_MSG(4057,(/vfile/),(/0.d0/))
CALL POLYX_READ_PARAM(vfile,Huc,Puc,H,Nnodes,vnodes,vorient,twodim,clearance,outparam,status)
!If return status is not zero there was an error => exit
IF( status>0 ) THEN
  GOTO 1000
ENDIF
!
!Check that number of nodes is not zero
IF(Nnodes<1) THEN
  CALL ATOMSK_MSG(4831,(/vfile/),(/0.d0/))
  nerr=nerr+1
  GOTO 1000
ENDIF
!
CALL ATOMSK_MSG(4058,(/''/),(/DBLE(twodim),H(1,1),H(2,2),H(3,3),DBLE(Nnodes)/))
!
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
IF( twodim > 0 ) THEN
  !System is pseudo 2-D => place all nodes on the same plane
  !in the middle of the cell along the short direction
  vnodes(:,twodim) = 0.5d0*H(twodim,twodim)
ENDIF
!
IF( outparam(4) .AND. ofu.NE.6 ) THEN
  !Write positions and orientations in a parameter file
  OPEN(41,FILE=outparamfile,STATUS="UNKNOWN")
  WRITE(41,'(a62)') "# Random positions and rotations of grains generated by Atomsk"
  WRITE(41,'(a66)') "# This parameter file can be used to generate the same polycrystal"
  WRITE(41,'(a4,3f16.6)') "box ", H(1,1), H(2,2), H(3,3)
  DO i=1,SIZE(vorient,1)
    !Compute corresponding rotation vectors and write them into parameter file
    CALL MAT2EULER_ZYX(vorient(i,:,:),P1,P2,P3)
    WRITE(41,'(a5,6f16.6)') "node ", vnodes(i,:), RAD2DEG(P1), RAD2DEG(P2), RAD2DEG(P3)
  ENDDO
  CLOSE(41)
ENDIF
!
IF(verbosity==4) THEN
  !Debug messages
  WRITE(msg,*) "Number of nodes:", Nnodes
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
!Check that the user provided enough information
IF( .NOT. ANY( NINT(H).NE.0 ) ) THEN
  !User did not provide box, or a box with negative or zero size
  nerr=nerr+1
  CALL ATOMSK_MSG(4820,(/""/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Estimate final number of particles NP = (density of unit cell) * (volume of final cell)
NP = CEILING( seed_density * VOLUME_PARA(H) )
NPexpect = NP
WRITE(msg,*) "Estimated number of atoms in this box: NP = ", NP
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!Allow for +10% or +1000 atoms to allocate arrays.
!Actual size of arrays will be adjusted later
NP = MIN( NINT(1.1d0*NP) , NP+1000 )
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
      P2 = 1.5d0*P1 + 10.d0
    ELSE
      P2 = 1.1d0*P1 + 6.d0
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
400 CONTINUE
!Estimate the size of template
Ht(:,:) = 0.d0
DO i=1,3
  distance = 0.d0
  DO inode=1,Nnodes
    n = Nvertices(inode)
    P1 = MAXVAL( vvertex(inode,1:n,i) )
    IF( P1>distance ) THEN
      distance = P1
    ENDIF
  ENDDO
  IF( Nnodes<6 ) THEN
    Ht(i,i) = 1.3d0*distance + 10.d0
  ELSE
    Ht(i,i) = 1.1d0*distance + 6.d0
  ENDIF
ENDDO
IF( twodim>0 ) THEN
  !System is pseudo-2D: template won't be duplicated along the short dimension
  !We can afford to make it larger along other directions
  Ht(:,:) = 1.5d0*Ht(:,:) + 20.d0
  Ht(twodim,twodim) = VECLENGTH(H(twodim,:))
ENDIF
!Estimate how many particles the template will contain = (seed density)*(template volume)
P1 = 1.05d0 * seed_density * Ht(1,1)*Ht(2,2)*Ht(3,3)
IF(verbosity>=4) THEN
  WRITE(msg,'(a,3f12.3)') "Template box size:", Ht(1,1), Ht(2,2), Ht(3,3)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,'(a25,f18.0)') "Expected NP for template:", P1
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!
!If number of nodes is small, or template size is small (<20,000 atoms),
!or if template size is very very large (i.e. several millions of atoms
!*and* much larger than expected NP in final polycrystal),
!then don't use a template
IF( Nnodes<5 .AND. P1>2.d6 .AND. P1>1.4d0*NP ) THEN
  use_template = .FALSE.
  WRITE(msg,*) "Template will NOT be constructed"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ELSEIF( P1<2.d4 ) THEN
  use_template = .TRUE.
  WRITE(msg,*) "Template WILL be constructed"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!
!Compute duplication factors along X, Y, Z
!If the system is 2-D, do not expand along the shortest axis
!NOTE: these duplication factors will be used to duplicate seed to fill grains
expandmatrix(:) = 0
WRITE(msg,*) "Determining expansion factors:"
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
DO i=1,3  !Loop on X, Y, Z
  IF( i.NE.twodim ) THEN
    !Round up to nearest greater integer
    expandmatrix(i) = MAX(1 , CEILING( Ht(i,i)/Huc(i,i) ) )
  ENDIF
ENDDO
WRITE(msg,*) "Initial expansion factors:", expandmatrix(:)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
420 CONTINUE
!If size (memory) is reasonable, construct template grain
!This template grain will be cut later to construct each grain
IF( use_template ) THEN
  ! Allocate array Pt for template grain (full duplicated crystal, not oriented or truncated)
  m = SIZE(Puc,1)*(expandmatrix(1)+1)*(expandmatrix(2)+1)*(expandmatrix(3)+1)
  !
  !Construct template supercell Pt(:,:)
  CALL ATOMSK_MSG(4081,(/''/),(/DBLE(m)/))
  WRITE(msg,*) "ALLOCATE  Pt, SIZE = ", m
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  ALLOCATE( Pt(m,4) , STAT=i )
  IF( i>0 ) THEN
    ! Allocation failed (not enough memory)
    ! Fall back to no-template method
    use_template = .FALSE.
    GOTO 500
  ENDIF
  Pt(:,:) = 0.d0
  IF( doshells ) THEN
    ALLOCATE( St( SIZE(Pt,1) , 4 ) , STAT=i )
    IF( i>0 ) THEN
      ! Allocation failed (not enough memory)
      ! Fall back to no-template method
      IF(ALLOCATED(Pt)) DEALLOCATE(Pt)
      use_template = .FALSE.
      GOTO 500
    ENDIF
    St(:,:) = 0.d0
  ENDIF
  IF( doaux ) THEN
    ALLOCATE( AUX_Q( SIZE(Pt,1) , SIZE(AUXuc,2)+1 ) , STAT=i )
  ELSE
    ALLOCATE( AUX_Q( SIZE(Pt,1) , 1 ) , STAT=i )
  ENDIF
  IF( i>0 ) THEN
    ! Allocation of AUX_Q failed (not enough memory)
    ! Fall back to no-template method
    IF(ALLOCATED(Pt)) DEALLOCATE(Pt)
    IF(ALLOCATED(St)) DEALLOCATE(St)
    use_template = .FALSE.
    GOTO 500
  ENDIF
  AUX_Q(:,:) = 0.d0
  !
  ! Fill array Pt(:,:)
  qi=0
  DO o = 0 , expandmatrix(3)
    DO n = 0 , expandmatrix(2)
      DO m = 0 , expandmatrix(1)
        DO i=1,SIZE(Puc,1)
          !Compute (cartesian) position of the replica of this atom
          qi = qi+1
          Pt(qi,1) = Puc(i,1) + DBLE(m)*Huc(1,1) + DBLE(n)*Huc(2,1) + DBLE(o)*Huc(3,1)
          Pt(qi,2) = Puc(i,2) + DBLE(m)*Huc(1,2) + DBLE(n)*Huc(2,2) + DBLE(o)*Huc(3,2)
          Pt(qi,3) = Puc(i,3) + DBLE(m)*Huc(1,3) + DBLE(n)*Huc(2,3) + DBLE(o)*Huc(3,3)
          Pt(qi,4) = Puc(i,4)
          IF(doshells) THEN
            !Compute (cartesian) position of the replica of this shell
            St(qi,1) = Suc(i,1) + DBLE(m)*Huc(1,1) + DBLE(n)*Huc(2,1) + DBLE(o)*Huc(3,1)
            St(qi,2) = Suc(i,2) + DBLE(m)*Huc(1,2) + DBLE(n)*Huc(2,2) + DBLE(o)*Huc(3,2)
            St(qi,3) = Suc(i,3) + DBLE(m)*Huc(1,3) + DBLE(n)*Huc(2,3) + DBLE(o)*Huc(3,3)
            St(qi,4) = Suc(i,4)
          ENDIF
          IF(doaux) THEN
            !Copy auxiliary properties of that atom
            AUX_Q(qi,:) = AUXuc(i,:)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  !Save template cell vectors (make sure it's not zero)
  Ht(1,:) = Huc(1,:) * MAX(1 , expandmatrix(1) )
  Ht(2,:) = Huc(2,:) * MAX(1 , expandmatrix(2) )
  Ht(3,:) = Huc(3,:) * MAX(1 , expandmatrix(3) )
  !
  IF( verbosity==4 ) THEN
    OPEN(UNIT=40,FILE="atomsk_template.xyz",STATUS="UNKNOWN")
    qi=0
    DO i=1,SIZE(Pt,1)
      IF(NINT(Pt(i,4)).NE.0) qi=qi+1
    ENDDO
    WRITE(40,'(i9)') qi
    WRITE(40,*) "# Template used by Atomsk to construct polycrystal"
    DO i=1,qi
      WRITE(40,'(i4,2X,3(f12.6))') NINT(Pt(i,4)), Pt(i,1), Pt(i,2), Pt(i,3)
    ENDDO
    CLOSE(40)
  ENDIF
  !
  !Allocate array for mask: will be .TRUE. for atoms inside polyhedron/grain, .FALSE. otherwise
  IF(ALLOCATED(Ptmask)) DEALLOCATE(Ptmask)
  ALLOCATE( Ptmask(SIZE(Pt,1)) )
  !
  !
ENDIF !end if(use_template)
!
!
!
500 CONTINUE
!Allocate arrays for final polycrystal
WRITE(msg,*) "             Allocating arrays with size: NP = ", NP
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ALLOCATE(Q(NP,4))
Q(:,:) = 0.d0
IF(doshells) THEN
  ALLOCATE(T(NP,4))
  T(:,:) = 0.d0
ENDIF
ALLOCATE(newAUX(NP,SIZE(AUXNAMES)))
newAUX(:,:) = 0.d0
!
!Allocate array containing number of atoms in each grain
ALLOCATE( NPgrains(Nnodes) , STAT=i )
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
  WRITE(msg,*) "  _  _  _  _  _  _  _"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,'(a,i4)') " |                     GRAIN #", inode
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Get center of grain = average positions of node + vertices
  n = Nvertices(inode)
  !
  IF( use_template ) THEN
    !For each atom of the template Pt, find out if it is located inside the grain
    !If not, set mask Ptmask(i) to .FALSE.
    Ptmask(:) = .TRUE.
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,vnormal,vector)
    DO i=1,SIZE(Pt,1)
      IF( Ptmask(i) ) THEN
        !Shift oriented supercell so that its center of mass is at the position of the node
        !Pt(i,1:3) = Pt(i,1:3) + shift(1:3)
        !Compute vector between atom #i and current node
        vector(:) = Pt(i,1:3) - vnodes(inode,1:3)
        !Check if atom #i is inside the grain #inode, loop on all vertices
        DO j=1,Nvertices(inode)
          !Compute vector between vertex #j and current node #inode
          !By definition this vector is normal to the grain boundary
          vnormal(:) = vvertex(inode,j,1:3) - vnodes(inode,1:3)
          IF( VEC_PLANE(vnormal,VECLENGTH(vnormal),vector) > -1.d0*clearance ) THEN
            !Atom #i is above plane of cut, hence out of polyhedron
            !Mark it for termination and exit loop on j (effectively close )
            Ptmask(i) = .FALSE.
            EXIT
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    !$OMP END PARALLEL DO
    !
    !Now Ptmask(:) is .TRUE. for template atoms that belong in grain #inode, .FALSE. otherwise
    !Copy all .TRUE. atoms into array Q
    !Apply rotation matrix and shift them into their final position
    DO i=1,SIZE(Pt,1)
      IF( Ptmask(i) ) THEN
        !Increment number of atoms belonging to grain #inode
        NPgrains(inode) = NPgrains(inode)+1
        !Increment total number of particles in the system
        NP = NP+1
        IF( NP>SIZE(Q,1) ) THEN
          !Total number of atoms exceeds size of array Q
          !This should not happen, yet here we are: resize Q and all related arrays
          m = MIN( SIZE(Q,1)+NPgrains(inode) , CEILING(1.1*SIZE(Q,1)) )
          CALL RESIZE_DBLEARRAY2( Q , m , 4 , n )
          IF( n>0 ) THEN
            ! Allocation failed (not enough memory)
            nerr = nerr+1
            CALL ATOMSK_MSG(819,(/''/),(/DBLE(m)/))
            GOTO 1000
          ENDIF
          !Resize array for shells with same size
          IF(doshells) CALL RESIZE_DBLEARRAY2( T , SIZE(Q,1) , 4 , n )
          IF( n>0 ) THEN
            ! Allocation failed (not enough memory)
            nerr = nerr+1
            CALL ATOMSK_MSG(819,(/''/),(/DBLE(SIZE(Q,1))/))
            GOTO 1000
          ENDIF
          !Resize array for aux. prop. with same size
          CALL RESIZE_DBLEARRAY2( newAUX , SIZE(Q,1) , SIZE(AUXNAMES) , n )
          IF( n>0 ) THEN
            ! Allocation failed (not enough memory)
            nerr = nerr+1
            CALL ATOMSK_MSG(819,(/''/),(/DBLE(SIZE(Q,1))/))
            GOTO 1000
          ENDIF
        ENDIF
        !Shift atom position and apply rotation
        P1 = Pt(i,1) - shiftgrain(inode,1)
        P2 = Pt(i,2) - shiftgrain(inode,2)
        P3 = Pt(i,3) - shiftgrain(inode,3)
        Q(NP,1) = P1*vorient(inode,1,1) + P2*vorient(inode,1,2) + P3*vorient(inode,1,3)
        Q(NP,2) = P1*vorient(inode,2,1) + P2*vorient(inode,2,2) + P3*vorient(inode,2,3)
        Q(NP,3) = P1*vorient(inode,3,1) + P2*vorient(inode,3,2) + P3*vorient(inode,3,3)
        Q(NP,4) = Pt(i,4)
        IF(doshells) THEN
          P1 = St(i,1) - shiftgrain(inode,1)
          P2 = St(i,2) - shiftgrain(inode,2)
          P3 = St(i,3) - shiftgrain(inode,3)
          T(NP,1) = P1*vorient(inode,1,1) + P2*vorient(inode,1,2) + P3*vorient(inode,1,3)
          T(NP,2) = P1*vorient(inode,2,1) + P2*vorient(inode,2,2) + P3*vorient(inode,2,3)
          T(NP,3) = P1*vorient(inode,3,1) + P2*vorient(inode,3,2) + P3*vorient(inode,3,3)
          T(NP,4) = St(i,4)
        ENDIF
        IF(doaux) THEN
          DO j=1,SIZE(AUXNAMES)-1
            newAUX(NP,j) = AUX_Q(i,j)
          ENDDO
        ENDIF
        newAUX(NP,grainID) = DBLE(inode)
      ENDIF
    ENDDO
    !
    IF( verbosity>=4 ) THEN
      WRITE(temp,'(f12.1)') 100.d0*DBLE(NPgrains(inode)) / DBLE(SIZE(Pt,1))
      WRITE(msg,*) " | INFO: template contains ", SIZE(Pt,1), " atoms, only ", NPgrains(inode), &
                & " (", TRIM(ADJUSTL(temp)), "%) were used in this grain"
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDIF
    !
    !
  ELSE !i.e. if there is no template
    !Seed will be duplicated inside the grain
    istart = NP+1
    !
    !Duplicate seed to fill the grain
    DO o=0,expandmatrix(3)
      DO n=0,expandmatrix(2)
        DO m=0,expandmatrix(1)
          !Loop on all atoms in seed/unit cell
          DO i=1,SIZE(Puc,1)
            !Compute (cartesian) position of the replica of this atom
            P1 = Puc(i,1) + DBLE(m)*Huc(1,1) + DBLE(n)*Huc(2,1) + DBLE(o)*Huc(3,1)
            P2 = Puc(i,2) + DBLE(m)*Huc(1,2) + DBLE(n)*Huc(2,2) + DBLE(o)*Huc(3,2)
            P3 = Puc(i,3) + DBLE(m)*Huc(1,3) + DBLE(n)*Huc(2,3) + DBLE(o)*Huc(3,3)
            vector(:) = (/P1,P2,P3/) - vnodes(inode,1:3)
            !Check if atom is inside the grain #inode
            j = 1
            isinpolyhedron = .TRUE.
            DO WHILE( isinpolyhedron .AND. j <= Nvertices(inode) )
              vnormal(:) = vvertex(inode,j,1:3) - vnodes(inode,1:3)
              IF( VEC_PLANE(vnormal,VECLENGTH(vnormal),vector) > -1.d0*clearance ) THEN
                !Atom is above this plane of cut, hence outside of polyhedron
                isinpolyhedron = .FALSE.
                EXIT
              ENDIF
              j=j+1
            ENDDO
            IF(isinpolyhedron) THEN
              !Atom #i survived all cuts => save its position and species into Q
              NPgrains(inode) = NPgrains(inode)+1
              NP = NP+1
              IF( NP > SIZE(Q,1) ) THEN
                !Increase array size by 10%
                PRINT*, "RESIZE Q"
                CALL RESIZE_DBLEARRAY2(Q,CEILING(1.1d0*NP),4,status)
                IF(doshells) THEN
                  CALL RESIZE_DBLEARRAY2(T,SIZE(Q,1),4,status)
                ENDIF
                IF(doaux) THEN
                  CALL RESIZE_DBLEARRAY2(newAUX,SIZE(Q,1),SIZE(AUXNAMES),status)
                ENDIF
              ENDIF
              !Shift atom position and apply rotation
              P1 = P1 - shiftgrain(inode,1)
              P2 = P2 - shiftgrain(inode,2)
              P3 = P3 - shiftgrain(inode,3)
              Q(NP,1) = P1*vorient(inode,1,1) + P2*vorient(inode,1,2) + P3*vorient(inode,1,3)
              Q(NP,2) = P1*vorient(inode,2,1) + P2*vorient(inode,2,2) + P3*vorient(inode,2,3)
              Q(NP,3) = P1*vorient(inode,3,1) + P2*vorient(inode,3,2) + P3*vorient(inode,3,3)
              Q(NP,4) = Puc(i,4)
              IF(doshells) THEN
                P1 = Suc(i,1) - shiftgrain(inode,1)
                P2 = Suc(i,2) - shiftgrain(inode,2)
                P3 = Suc(i,3) - shiftgrain(inode,3)
                T(NP,1) = P1*vorient(inode,1,1) + P2*vorient(inode,1,2) + P3*vorient(inode,1,3)
                T(NP,2) = P1*vorient(inode,2,1) + P2*vorient(inode,2,2) + P3*vorient(inode,2,3)
                T(NP,3) = P1*vorient(inode,3,1) + P2*vorient(inode,3,2) + P3*vorient(inode,3,3)
                T(NP,4) = Suc(i,4)
              ENDIF
              IF(doaux) THEN
                newAUX(NP,:) = AUXuc(i,:)
              ENDIF
              newAUX(NP,grainID) = DBLE(inode)
            ENDIF
          ENDDO !loop on i
        ENDDO !loop on m
      ENDDO !loop on n
    ENDDO !loop on o
    !
  ENDIF !end if(use_template)
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
    WRITE(msg,*) " |_  _  _  _  _  END OF GRAIN #", inode
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
  !
ENDDO  !end loop on inode
!
!
!
600 CONTINUE
WRITE(msg,*) "Loop on grains finished, size of Q = ", SIZE(Q,1)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
IF(ALLOCATED(Puc)) DEALLOCATE(Puc)
IF(ALLOCATED(Suc)) DEALLOCATE(Suc)
IF(ALLOCATED(Pt)) DEALLOCATE(Pt)
IF(ALLOCATED(St)) DEALLOCATE(St)
IF(ALLOCATED(Ptmask)) DEALLOCATE(Ptmask)
IF(ALLOCATED(Nvertices)) DEALLOCATE(Nvertices)
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
!Generate comment for output file(s)
ALLOCATE(comment(1))
WRITE(temp,*) Nnodes
 comment(1) = "# Voronoi polycrystal with "//TRIM(ADJUSTL(temp))//" grains"
!
!Apply options to the final system
CALL OPTIONS_AFF(options_array,Huc,H,P,S,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
IF(nerr>0) GOTO 1000
!
!Write final system to file(s)
IF(wof) THEN
  CALL WRITE_AFF(prefix,outfileformats,H,P,S,comment,AUXNAMES,AUX)
ENDIF
IF(nerr>0) GOTO 1000
!
!
!
700 CONTINUE
!From now on, some arrays are not needed anymore => deallocate them
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(S)) DEALLOCATE(S)
IF(ALLOCATED(T)) DEALLOCATE(T)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF(ALLOCATED(comment)) DEALLOCATE(comment)
IF(ALLOCATED(outfileformats)) DEALLOCATE(outfileformats)
!
IF( outparam(4) .AND. ofu.NE.6 ) THEN
  !Random parameters were written into a file => inform user
  temp = "parameter"
  CALL ATOMSK_MSG(3002,(/outparamfile,temp/),(/0.d0/))
ENDIF
!
IF( wof .AND. ofu.NE.6 ) THEN
  !Write positions of nodes into a file
  ALLOCATE(comment(1))
  comment(1) = "# Positions of nodes"
  ALLOCATE(outfileformats(1))
  outfileformats(1) = "xsf"
  ALLOCATE(Q(SIZE(vnodes,1),4))
  DO i=1,SIZE(vnodes,1)
    Q(i,1:3) = vnodes(i,1:3)
    Q(i,4) = 1.d0  !nodes are given the atomic number of hydrogen
  ENDDO
  CALL NAME_OUTFILE(prefix,temp,"     ")
  i=SCAN(temp,".",BACK=.TRUE.)  !remove trailing dot
  temp = temp(1:i-1)
  outparamfile = TRIM(ADJUSTL(temp))//"_nodes.xsf"
  CALL WRITE_AFF(outparamfile,outfileformats,H,Q,S,comment,AUXNAMES,AUX)
  IF(ALLOCATED(Q)) DEALLOCATE(Q)
  IF(ALLOCATED(comment)) DEALLOCATE(comment)
  IF(ALLOCATED(outfileformats)) DEALLOCATE(outfileformats)
ENDIF
!
!Find the max. volume occupied by a grain
Vmin = MINVAL(cell_volumes)
Vmax = MAXVAL(cell_volumes)
m=0
ALLOCATE(Q(SIZE(NPgrains),4))
DO i=1,SIZE(NPgrains)
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
IF( wof .AND. ofu.NE.6 ) THEN
  !Write positions of center of mass of each grain into a file
  ALLOCATE(comment(1))
  comment(1) = "# Positions of grains centers of mass"
  ALLOCATE(outfileformats(1))
  outfileformats(1) = "xsf"
  CALL NAME_OUTFILE(prefix,temp,"     ")
  i=SCAN(temp,".",BACK=.TRUE.)  !remove trailing dot
  temp = temp(1:i-1)
  outparamfile = TRIM(ADJUSTL(temp))//"_grains-com.xsf"
  CALL WRITE_AFF(outparamfile,outfileformats,H,Q,S,comment,AUXNAMES,AUX)
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
    DO i=1,SIZE(NPgrains)  !loop on all grains
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
      DO i=1,SIZE(NPgrains) !loop on all grains
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
CALL ATOMSK_MSG(1801,(/TRIM(vfile)/),(/DBLE(linenumber)/))
GOTO 1000
!
!
!
1000 CONTINUE
IF(ALLOCATED(Puc)) DEALLOCATE(Puc)
IF(ALLOCATED(Suc)) DEALLOCATE(Suc)
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(S)) DEALLOCATE(S)
IF(ALLOCATED(T)) DEALLOCATE(T)
IF(ALLOCATED(Nvertices)) DEALLOCATE(Nvertices)
IF(ALLOCATED(NPgrains)) DEALLOCATE(NPgrains)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
!
!
!
END SUBROUTINE POLYCRYS
!
!
!
END MODULE mode_polycrystal
