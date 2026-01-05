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
!* Last modification: P. Hirel - 05 Jan. 2026                                     *
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
INTEGER:: NP   !total number of atoms in the final system
INTEGER:: qi   !used to count atoms in a grain
INTEGER:: inode, jnode
INTEGER:: Nnodes      !number of nodes
INTEGER:: status
INTEGER,DIMENSION(3):: expandmatrix
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
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newP, newS !positions of atoms, shells (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUXuc     !auxiliary properties of atoms in the unit cell (seed)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX_Q     !auxiliary properties of atoms in a grain
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX       !auxiliary properties of atoms in the final supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX    !auxiliary properties of atoms (temporary)
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
use_template = .TRUE.
m=0
Nnodes = 0
twodim = 0         !assume system will be 3-D
clearance = 0.1d0  !atoms this close to the GB plane will be deleted
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
IF(ALLOCATED(NPgrains)) DEALLOCATE(NPgrains)
 C_tensor(:,:) = 0.d0
IF(ALLOCATED(P)) DEALLOCATE(P)
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
CALL ATOMSK_MSG(4058,(/''/),(/DBLE(Nnodes),DBLE(twodim)/))
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
!
!
400 CONTINUE
!Construct template supercell Pt(:,:)
IF( use_template ) THEN
  CALL ATOMSK_MSG(4081,(/''/),(/0.d0/))
  IF( SIZE(vvertex,1)<=3 ) THEN
    Ht(1,1) = MAXVAL(DABS(H(:,:))) * 1.45d0
    Ht(2,2) = Ht(1,1)
    Ht(3,3) = Ht(1,1)
  ELSE
    !Determine max. distance between two vertices of a given node
    maxdvertices = 0.d0
    DO inode=1,SIZE(vvertex,1)  !Loop on all nodes
      n = Nvertices(inode)
      DO j=1,n-1   !Loop on all vertices of that node
        DO k=j+1,n  !Loop on all other vertices of that node
          distance = VECLENGTH( vvertex(inode,j,1:3) - vvertex(inode,k,1:3) )
          IF( distance > maxdvertices ) THEN
            maxdvertices = distance
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    WRITE(msg,'(a,f12.3)') "Max distance between two vertices:", maxdvertices
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !
    !Save box into Ht(:,:), add a few angströms for good measure
    Ht(:,:) = 0.d0
    DO i=1,3
      Ht(i,i) = MIN( maxdvertices+10.d0 , 1.1d0*maxdvertices )
    ENDDO
  ENDIF
  IF( twodim>0 ) THEN
    !System is pseudo-2D: don't duplicate template along the short dimension
    Ht(:,:) = Ht(:,:)*2.25d0
    Ht(twodim,twodim) = VECLENGTH(H(twodim,:))
  ENDIF
  !Compute how many particles the template will contain = (seed density)*(template volume)
  P1 = 1.05d0 * seed_density * Ht(1,1)*Ht(2,2)*Ht(3,3)
  WRITE(msg,'(a25,f18.0)') "Expected NP for template:", P1
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !If expected number of particles is too large, reduce template size
  IF( P1 > 2.1d9 ) THEN
    Ht(:,:) = 0.8d0*Ht(:,:)
    IF( twodim>0 ) THEN
      Ht(twodim,twodim) = VECLENGTH(H(twodim,:))
    ENDIF
    !Re-compute number of particles
    P1 = 1.05d0 * seed_density * Ht(1,1)*Ht(2,2)*Ht(3,3)
  ENDIF
  !If expected number of particles still is too large, abort completely
  IF( P1 > 2.1d9 ) THEN
    CALL ATOMSK_MSG(821,(/""/),(/P1/))
    nerr = nerr+1
    GOTO 1000
  ENDIF
  !
  !By default the template grain is a bit larger than max.cell size * sqrt(3)
  !This template grain will be cut later to construct each grain
  expandmatrix(:) = 1
  WRITE(msg,*) "Determining expansion factors:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,3
    IF( VECLENGTH(Huc(i,:)) < 0.8d0*MAXVAL(Ht(:,:)) ) THEN
      !
      !Compute sum of seed vectors components along direction i
      P2 = DBLE( FLOOR( DABS( SUM(Huc(:,i)) )))
      P2 = MAX( P2 , 1.d0 )
      !
      !P3 = number of times the seed will be duplicated along each base vector direction
      P3 = 1.05d0 * ( Ht(i,i) / P2 )
      !
      !Make sure duplication factors are not crazy
      IF(P3<=0) THEN
        P3 = 1
      ELSEIF(P3>2000) THEN
        P3=2000
      ENDIF
      expandmatrix(i) = CEILING(P3)
    ENDIF
  ENDDO
  !If the system is 2-D, do not expand along the shortest axis
  IF( twodim>0 ) THEN
    expandmatrix(twodim) = 0
  ENDIF
  WRITE(msg,*) "Initial expansion factors:", expandmatrix(:)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  ! Allocate array newP for template grain (full duplicated crystal, not oriented or truncated)
  m = CEILING(P1)
  WRITE(msg,*) "ALLOCATE  newP, SIZE = ", m
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ALLOCATE( newP(m,4) , STAT=i )
  IF( i>0 ) THEN
    ! Allocation failed (not enough memory)
    nerr = nerr+1
    CALL ATOMSK_MSG(819,(/''/),(/DBLE(m)/))
    GOTO 1000
  ENDIF
  newP(:,:) = 0.d0
  IF( doshells ) THEN
    ALLOCATE( newS( SIZE(newP,1) , 4 ) , STAT=i )
    IF( i>0 ) THEN
      ! Allocation failed (not enough memory)
      nerr = nerr+1
      CALL ATOMSK_MSG(819,(/''/),(/DBLE(m)/))
      GOTO 1000
    ENDIF
    newS(:,:) = 0.d0
  ENDIF
  IF( doaux ) THEN
    ALLOCATE( newAUX( SIZE(newP,1) , SIZE(AUXuc,2)+1 ) , STAT=i )
    IF( i>0 ) THEN
      ! Allocation failed (not enough memory)
      nerr = nerr+1
      CALL ATOMSK_MSG(819,(/''/),(/DBLE(m)/))
      GOTO 1000
    ENDIF
  ELSE
    ALLOCATE( newAUX( SIZE(newP,1) , 1 ) , STAT=i )
    IF( i>0 ) THEN
      ! Allocation failed (not enough memory)
      nerr = nerr+1
      CALL ATOMSK_MSG(819,(/''/),(/DBLE(m)/))
      GOTO 1000
    ENDIF
  ENDIF
  newAUX(:,:) = 0.d0
  !
  ! Fill array newP(:,:)
  NP=0
  DO o = 0 , expandmatrix(3)
    DO n = 0 , expandmatrix(2)
      DO m = 0 , expandmatrix(1)
        DO i=1,SIZE(Puc,1)
          !Compute (cartesian) position of the replica of this atom
          P1 = Puc(i,1) + DBLE(m)*Huc(1,1) + DBLE(n)*Huc(2,1) + DBLE(o)*Huc(3,1)
          P2 = Puc(i,2) + DBLE(m)*Huc(1,2) + DBLE(n)*Huc(2,2) + DBLE(o)*Huc(3,2)
          P3 = Puc(i,3) + DBLE(m)*Huc(1,3) + DBLE(n)*Huc(2,3) + DBLE(o)*Huc(3,3)
          distance = VECLENGTH( (/P1,P2,P3/) - 0.5d0*(/Ht(1,1),Ht(2,2),Ht(3,3)/) )
          !Check if this position is inside template box
          IF( P1>-1.d0 .AND. P1<=Ht(1,1) .AND.   &
            & P2>-1.d0 .AND. P2<=Ht(2,2) .AND.   &
            & P3>-1.d0 .AND. P3<=Ht(3,3) .AND.   &
            & ( twodim>0 .OR. distance<0.72d0*MAXVAL(Ht(:,:)) )  ) THEN
            !Yes it is: save atom position and species into newP
            NP = NP+1
            IF( NP > SIZE(newP,1) ) THEN
              !Increase array size by 10%
              CALL RESIZE_DBLEARRAY2(newP,CEILING(1.1d0*NP),4,status)
              IF(doshells) THEN
                CALL RESIZE_DBLEARRAY2(newS,SIZE(newP,1),4,status)
              ENDIF
            ENDIF
            newP(NP,1) = P1
            newP(NP,2) = P2
            newP(NP,3) = P3
            newP(NP,4) = Puc(i,4)
            IF(doshells) THEN
              !Compute (cartesian) position of the replica of this shell
              newS(NP,1) = Suc(i,1) + DBLE(m)*Huc(1,1) + DBLE(n)*Huc(2,1) + DBLE(o)*Huc(3,1)
              newS(NP,2) = Suc(i,2) + DBLE(m)*Huc(1,2) + DBLE(n)*Huc(2,2) + DBLE(o)*Huc(3,2)
              newS(NP,3) = Suc(i,3) + DBLE(m)*Huc(1,3) + DBLE(n)*Huc(2,3) + DBLE(o)*Huc(3,3)
              newS(NP,4) = Puc(i,4)
            ENDIF
            IF(doaux) THEN
              !Copy auxiliary properties of that atom
              newAUX(NP,:) = AUXuc(i,:)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  !Save template into Pt with correct array size
  IF(ALLOCATED(Pt)) DEALLOCATE(Pt)
  ALLOCATE( Pt(NP,4) )
  Pt(:,:) = 0.d0
  IF(doshells) THEN
    IF(ALLOCATED(St)) DEALLOCATE(St)
    ALLOCATE( St(NP,4) )
    St(:,:) = 0.d0
  ENDIF
  IF(doaux) THEN
    IF(ALLOCATED(AUX_Q)) DEALLOCATE(AUX_Q)
    ALLOCATE( AUX_Q(NP,4) )
    AUX_Q(:,:) = 0.d0
  ENDIF
  DO i=1,NP
    Pt(i,:) = newP(i,:)
    IF(doshells) THEN
      St(i,:) = newS(i,:)
    ENDIF
    IF(doaux) THEN
      AUX_Q(i,:) = newAUX(i,:)
    ENDIF
  ENDDO
  DEALLOCATE(newP)
  IF(ALLOCATED(newS)) DEALLOCATE(newS)
  IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
  !
  !Save template cell vectors
  Ht(1,:) = Huc(1,:) * MAX(1 , expandmatrix(1) )
  Ht(2,:) = Huc(2,:) * MAX(1 , expandmatrix(2) )
  Ht(3,:) = Huc(3,:) * MAX(1 , expandmatrix(3) )
  !
  !
  IF( verbosity==4 ) THEN
    OPEN(UNIT=40,FILE="atomsk_template.xyz",STATUS="UNKNOWN")
    NP=0
    DO i=1,SIZE(Pt,1)
      IF(NINT(Pt(i,4)).NE.0) NP=NP+1
    ENDDO
    WRITE(40,'(i9)') NP
    WRITE(40,*) "# Template used by Atomsk to construct polycrystal"
    DO i=1,NP
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
!Estimate final number of particles NP = (density of unit cell) * (volume of final cell)
NP = CEILING( seed_density * DABS(H(1,1)*H(2,2)*H(3,3)) )
WRITE(msg,*) "Estimated number of atoms in polycrystal: NP = ", NP
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!Allow for +10% or +1000 atoms to allocate arrays.
!Actual size of arrays will be adjusted later
NP = MIN( NINT(1.2d0*NP) , NP+2000 )
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
Ht(:,:) = Huc(:,:)
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
  vector(:) = 0.d0
  DO jnode=1,n
    vector(:) = vector(:) + vvertex(inode,jnode,1:3)
  ENDDO
  GrainCenter(:) = ( vnodes(inode,1:3) + vector(:) ) / DBLE(n+1)
  !
  IF( verbosity==4 ) THEN
    WRITE(msg,'(a23,3f9.3)') "  | Node position:      ", vnodes(inode,1:3)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,'(a23,3f9.3)') "  | Center of vertices: ", vector(1:3) / DBLE(n)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,'(a23,3f9.3)') "  | Grain center:       ", GrainCenter(:)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
  !
  !This node is surrounded by Nvertices(inode) vertices
  !Positions of vertices are in array vvertex(inode,:,:)
  CALL BUBBLESORT(vvertex(inode,1:n,1:4),4,"up  ",newindex)
  !
  IF( use_template ) THEN
    !Rotate template grain to match target grain orientation
    !NOTE: unless inode=1, the template was already rotated before
    !so we have to un-rotate it (apply the transpose of previous rotation matrix)
    !before applying the rotation of current inode
    rotmat(:,:) = MATMUL( TRANSPOSE(rotmat) , vorient(inode,1:3,1:3) )
    DO i=1,SIZE(Pt,1)
      H1 = Pt(i,1)
      H2 = Pt(i,2)
      H3 = Pt(i,3)
      Pt(i,1) = H1*rotmat(1,1) + H2*rotmat(1,2) + H3*rotmat(1,3)
      Pt(i,2) = H1*rotmat(2,1) + H2*rotmat(2,2) + H3*rotmat(2,3)
      Pt(i,3) = H1*rotmat(3,1) + H2*rotmat(3,2) + H3*rotmat(3,3)
    ENDDO
    IF( doshells ) THEN
      DO i=1,SIZE(St,1)
        H1 = St(i,1)
        H2 = St(i,2)
        H3 = St(i,3)
        St(i,1) = H1*rotmat(1,1) + H2*rotmat(1,2) + H3*rotmat(1,3)
        St(i,2) = H1*rotmat(2,1) + H2*rotmat(2,2) + H3*rotmat(2,3)
        St(i,3) = H1*rotmat(3,1) + H2*rotmat(3,2) + H3*rotmat(3,3)
      ENDDO
    ENDIF
    IF( verbosity==4 ) THEN
      WRITE(msg,'(a18,3f9.3)') "  | vorient(1):   ", vorient(inode,1,1:3)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      WRITE(msg,'(a18,3f9.3)') "  | vorient(2):   ", vorient(inode,2,1:3)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      WRITE(msg,'(a18,3f9.3)') "  | vorient(3):   ", vorient(inode,3,1:3)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      WRITE(msg,'(a18,3f9.3)') "  | rotmat(1):    ", rotmat(1,1:3)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      WRITE(msg,'(a18,3f9.3)') "  | rotmat(2):    ", rotmat(2,1:3)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      WRITE(msg,'(a18,3f9.3)') "  | rotmat(3):    ", rotmat(3,1:3)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDIF
    !CALL ATOMSK_MSG(2072,(/''/),(/0.d0/))
    !
    !Determine current position of the center of mass of template
    n = SIZE(Pt,1)
    vector(:) = 0.5d0*( Pt(1,1:3) + Pt(n,1:3) )
    !
    !vector(:) is the position of center of mass of template
    !Compute shift vector so that center of mass is at the center of grain
    shift(1:3) = GrainCenter(1:3) - vector(1:3)
    !
    !For each atom of the template Pt, find out if it is located inside the grain
    !If not, set mask Ptmask(i) to .FALSE.
    Ptmask(:) = .TRUE.
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,vnormal,vector)
    DO i=1,SIZE(Pt,1)
      IF( Ptmask(i) ) THEN
        !Shift oriented supercell so that its center of mass is at the position of the node
        Pt(i,1:3) = Pt(i,1:3) + shift(1:3)
        !Compute vector between atom #i and current node
        vector(:) = Pt(i,1:3) - vnodes(inode,1:3)
        !Check if atom is inside the grain #inode
        j = 1
        DO WHILE( Ptmask(i) .AND. j <= Nvertices(inode) )
          !Compute vector between vertex #j and current node #inode
          !By definition this vector is normal to the grain boundary
          vnormal(:) = vvertex(inode,j,1:3) - vnodes(inode,1:3)
          IF( VEC_PLANE(vnormal,VECLENGTH(vnormal),vector) > -1.d0*clearance ) THEN
            !Atom is above plane of cut, hence out of polyhedron => exit loop on j
            Ptmask(i) = .FALSE.
            EXIT
          ENDIF
          j=j+1
        ENDDO
      ENDIF
      !
    ENDDO
    !$OMP END PARALLEL DO
    !
    !Now Ptmask(:) is .TRUE. for atoms that belong to grain #inode, .FALSE. otherwise
    !Copy all .TRUE. atoms into array Q
    DO i=1,SIZE(Pt,1)
      IF( Ptmask(i) ) THEN
        !Increment number of atoms belonging to grain #inode
        NPgrains(inode) = NPgrains(inode)+1
        !Increment total number of particles in the system
        NP = NP+1
        IF( NP>SIZE(Q,1) ) THEN
          !Array Q was too small: resize it
          CALL RESIZE_DBLEARRAY2( Q , NP+1000 , 4 , n )
          IF( n>0 ) THEN
            ! Allocation failed (not enough memory)
            nerr = nerr+1
            CALL ATOMSK_MSG(819,(/''/),(/DBLE(NP+1000)/))
            GOTO 1000
          ENDIF
          IF(doshells) CALL RESIZE_DBLEARRAY2( T , NP+1000 , 4 , n )
          IF( n>0 ) THEN
            ! Allocation failed (not enough memory)
            nerr = nerr+1
            CALL ATOMSK_MSG(819,(/''/),(/DBLE(NP+1000)/))
            GOTO 1000
          ENDIF
          CALL RESIZE_DBLEARRAY2( newAUX , NP+1000 , SIZE(AUXNAMES) , n )
          IF( n>0 ) THEN
            ! Allocation failed (not enough memory)
            nerr = nerr+1
            CALL ATOMSK_MSG(819,(/''/),(/DBLE(NP+1000)/))
            GOTO 1000
          ENDIF
        ENDIF
        Q(NP,:) = Pt(i,:)
        IF(doshells) T(NP,:) = St(i,:)
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
    istart = NP+1
    !Rotate unit cell to match target grain orientation
    rotmat(:,:) = MATMUL( TRANSPOSE(rotmat) , vorient(inode,1:3,1:3) )
    DO i=1,3
      H1 = Ht(i,1)
      H2 = Ht(i,2)
      H3 = Ht(i,3)
      Ht(i,1) = H1*rotmat(1,1) + H2*rotmat(1,2) + H3*rotmat(1,3)
      Ht(i,2) = H1*rotmat(2,1) + H2*rotmat(2,2) + H3*rotmat(2,3)
      Ht(i,3) = H1*rotmat(3,1) + H2*rotmat(3,2) + H3*rotmat(3,3)
    ENDDO
    DO i=1,SIZE(Puc,1)
      H1 = Puc(i,1)
      H2 = Puc(i,2)
      H3 = Puc(i,3)
      Puc(i,1) = H1*rotmat(1,1) + H2*rotmat(1,2) + H3*rotmat(1,3)
      Puc(i,2) = H1*rotmat(2,1) + H2*rotmat(2,2) + H3*rotmat(2,3)
      Puc(i,3) = H1*rotmat(3,1) + H2*rotmat(3,2) + H3*rotmat(3,3)
    ENDDO
    IF( doshells ) THEN
      DO i=1,SIZE(Suc,1)
        H1 = Suc(i,1)
        H2 = Suc(i,2)
        H3 = Suc(i,3)
        Suc(i,1) = H1*rotmat(1,1) + H2*rotmat(1,2) + H3*rotmat(1,3)
        Suc(i,2) = H1*rotmat(2,1) + H2*rotmat(2,2) + H3*rotmat(2,3)
        Suc(i,3) = H1*rotmat(3,1) + H2*rotmat(3,2) + H3*rotmat(3,3)
      ENDDO
    ENDIF
    !CALL ATOMSK_MSG(2072,(/''/),(/0.d0/))
    !
    !Place origin of unit cell at the node position
    !Save it
    DO i=1,SIZE(Puc,1)
      Puc(i,1:3) = vnodes(inode,1:3) - Puc(i,1:3)
      NP = NP+1
      NPgrains(inode) = NPgrains(inode)+1
      Q(NP,:) = Puc(i,:)
      IF(doshells) T(NP,:) = Suc(i,:)
      IF(doaux) newAUX(NP,:) = AUXuc(i,:)
    ENDDO
    !
    !Duplicate unit cell to fill the grain
    !expandmatrix(:) is used to give the sign (direction) of duplication along X,Y,Z
    m=0
    n=0
    o=0
    k=1
    DO WHILE(k>0)
      o=o+1
      DO WHILE(k>0)
        n=n+1
        DO WHILE(k>0)  !loop on m
          IF(twodim>0) expandmatrix(twodim)=0
          PRINT*, "ITERATION: ", m, n, o
          k=0  !counter for atoms inside the grain
          DO i=1,SIZE(Puc,1)
            !Compute (cartesian) position of the replica of this atom
            P1 = Puc(i,1) + DBLE(m)*Ht(1,1) + DBLE(n)*Ht(2,1) + DBLE(o)*Ht(3,1)
            P2 = Puc(i,2) + DBLE(m)*Ht(1,2) + DBLE(n)*Ht(2,2) + DBLE(o)*Ht(3,2)
            P3 = Puc(i,3) + DBLE(m)*Ht(1,3) + DBLE(n)*Ht(2,3) + DBLE(o)*Ht(3,3)
            vector(:) = (/P1,P2,P3/) - vnodes(inode,1:3)
            !Check if atom is inside the grain #inode
            j = 1
            isinpolyhedron = .TRUE.
            DO WHILE( isinpolyhedron .AND. j <= Nvertices(inode) )
              vnormal(:) = vvertex(inode,j,1:3) - vnodes(inode,1:3)
              IF( VEC_PLANE(vnormal,VECLENGTH(vnormal),vector) > -1.d0*clearance ) THEN
                !Atom is above this plane of cut, hence outside of polyhedron
                isinpolyhedron = .FALSE.
              ENDIF
              j=j+1
            ENDDO
            IF(isinpolyhedron) THEN
              !Atom #i survived all cuts => save its position and species into Q
              k = k+1
              NPgrains(inode) = NPgrains(inode)+1
              NP = NP+1
              IF( NP > SIZE(Q,1) ) THEN
                !Increase array size by 10%
                CALL RESIZE_DBLEARRAY2(Q,CEILING(1.1d0*NP),4,status)
                IF(doshells) THEN
                  CALL RESIZE_DBLEARRAY2(T,SIZE(Q,1),4,status)
                ENDIF
                IF(doaux) THEN
                  CALL RESIZE_DBLEARRAY2(newAUX,SIZE(Q,1),SIZE(AUXNAMES),status)
                ENDIF
              ENDIF
              Q(NP,1) = P1
              Q(NP,2) = P2
              Q(NP,3) = P3
              Q(NP,4) = Puc(i,4)
              IF(doshells) THEN
                T(NP,1) = Suc(i,1) + DBLE(m)*Ht(1,1) + DBLE(n)*Ht(2,1) + DBLE(o)*Ht(3,1)
                T(NP,2) = Suc(i,2) + DBLE(m)*Ht(1,2) + DBLE(n)*Ht(2,2) + DBLE(o)*Ht(3,2)
                T(NP,3) = Suc(i,3) + DBLE(m)*Ht(1,3) + DBLE(n)*Ht(2,3) + DBLE(o)*Ht(3,3)
                T(NP,4) = Suc(i,4)
              ENDIF
              IF(doaux) THEN
                newAUX(NP,:) = AUXuc(i,:)
              ENDIF
            ENDIF
          ENDDO !loop on i
          IF( k==0 ) THEN
            !If k==0 it will terminate the loop on m
            k=1
          ELSE
            m=m+1
          ENDIF
        ENDDO !loop on m
        IF( k>0 ) THEN
          n=n+1
        ENDIF
        !If k==0 it will terminate the loop on m
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
DO i=1,NP
  P(i,:) = Q(i,:)
ENDDO
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
  DO i=1,NP
    S(i,:) = T(i,:)
  ENDDO
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
DO i=1,NP
  AUX(i,:) = newAUX(i,:)
ENDDO
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
