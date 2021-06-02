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
!* Last modification: P. Hirel - 02 June 2021                                     *
!**********************************************************************************
!* OUTLINE:                                                                       *
!* 100        Read atom positions of seed (usually a unit cell) from ucfile       *
!* 200        Read parameters for Voronoi construction from vfile                 *
!* 300        Construct grains using Voronoi method                               *
!* 500        Apply options and write final result to output file(s)              *
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
USE neighbors
USE files
USE subroutines
USE readin
USE options
USE center
USE orient
USE writeout
!
CONTAINS
!
!
SUBROUTINE POLYCRYS(ucfile,vfile,options_array,prefix,outfileformats)
!
!
IMPLICIT NONE
!Input parameters
CHARACTER(LEN=*),INTENT(IN):: ucfile  !name of file containing seed (usually a unit cell, but can be anything)
CHARACTER(LEN=*),INTENT(IN):: vfile   !name of file containing parameters for Voronoi construction
CHARACTER(LEN=*),INTENT(IN):: prefix  !name or prefix for output file (polycrystal)
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to output
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
!
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: line
CHARACTER(LEN=128):: or1, or2, or3
CHARACTER(LEN=128):: lattice  !if grains are organized according to a lattice
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=4096):: outparamfile  !file where grain parameters are written (if some parameters equal "random")
CHARACTER(LEN=4096):: distfile, idsizefile !name of file containing grain size distribution, grain sizes
CHARACTER(LEN=32),DIMENSION(3):: oldvec, newvec !new vectors after orient
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES    !names of auxiliary properties of atoms
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: newAUXNAMES !names of auxiliary properties of atoms (temporary)
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: doshells, doaux !are there shells, auxiliary properties in initial seed?
LOGICAL:: Hset            !are the box vectors H(:,:) defined?
LOGICAL:: isinpolyhedron  !is atom inside the polyhedron?
LOGICAL:: miller          !are Miller indices given? (if no then angles are given)
LOGICAL:: paramnode, paramrand, paramlatt !are the keywords "node", "random", "lattice" used in parameter file?
LOGICAL:: outparam        !were parameters saved in a text file?
LOGICAL:: sameplane       !are all nodes in the same plane?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT
INTEGER:: twodim        !=0 if system is 3-D, =1,2,3 if system is thin along x, y, z
INTEGER:: grainID       !position of the auxiliary property "grainID" in AUX
INTEGER:: i, j, k
INTEGER:: linenumber    !line number when reading a file
INTEGER:: m, n, o
INTEGER:: maxvertex   !max. number of vertices to look for (defined for 3-D or 2-D)
INTEGER:: Nvertexmax  !max number of neighboring vertices found
INTEGER:: NP   !total number of atoms in the final system
INTEGER:: qi   !used to count atoms in a grain
INTEGER:: inode, jnode
INTEGER:: Nnodes, Nvertices !number of nodes, of vertices
INTEGER:: status
INTEGER,DIMENSION(:),ALLOCATABLE:: NPgrains  !number of atoms in each grain
INTEGER,DIMENSION(3):: expandmatrix
INTEGER,DIMENSION(:),ALLOCATABLE:: newindex  !list of index after sorting
INTEGER,DIMENSION(:,:),ALLOCATABLE:: vnodesNeighList  !list of neighbours for nodes
REAL(dp):: boxmax      !max. distance from one end of the box to another
REAL(dp):: distance    !distance between two points
REAL(dp):: maxdnodes   !maximum distance between 2 nodes
REAL(dp):: P1, P2, P3  !temporary position
REAL(dp):: Volume, Vmin, Vmax, Vstep  !min, max. volume occupied by a grain, step for grain size distribution
REAL(dp),DIMENSION(3):: GrainCenter !position of center of grain (different from node position)
REAL(dp),DIMENSION(3):: vector    !vector between an atom and a node
REAL(dp),DIMENSION(3):: vnormal   !vector normal to grain boundary
REAL(dp),DIMENSION(3,3):: Huc       !Base vectors of the unit cell (seed)
REAL(dp),DIMENSION(3,3):: Ht        !Base vectors of the oriented unit cell
REAL(dp),DIMENSION(3,3):: Gt        !Inverse of Ht
REAL(dp),DIMENSION(3,3):: Hn, Hend  !Normalized Ht
REAL(dp),DIMENSION(3,3):: H         !Base vectors of the final supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystalographic orientation
REAL(dp),DIMENSION(3,3):: rotmat  !rotation matrix
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:),ALLOCATABLE:: randarray   !random numbers
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Puc, Suc  !positions of atoms, shells in unit cell (seed)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Pt, St    !positions of atoms, shells in template supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Pt2, St2  !copy of template supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, S      !positions of atoms, shells in final supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Q, T      !positions of atoms, shells in a grain
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newP, newS !positions of atoms, shells (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUXuc     !auxiliary properties of atoms in the unit cell (seed)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX_Q     !auxiliary properties of atoms in a grain
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX       !auxiliary properties of atoms in the final supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX    !auxiliary properties of atoms (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: vnodes    !cartesian coordinate of each node
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: vvertex   !cartesian coordinate of each vertex
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE:: vorient !crystallographic orientation of each node
!
!
!Initialize variables
CALL NAME_OUTFILE(prefix,temp,"     ")
i=SCAN(temp,".",BACK=.TRUE.)  !remove trailing dot
temp = temp(1:i-1)
distfile = TRIM(ADJUSTL(temp))//"_size-dist.txt"
idsizefile = TRIM(ADJUSTL(temp))//"_id-size.txt"
outparamfile = TRIM(ADJUSTL(temp))//"_param.txt"
outparam = .FALSE.
sameplane = .FALSE.
Nnodes = 0
twodim = 0  !assume system will be 3-D
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
IF( ALLOCATED(comment) ) DEALLOCATE(comment)
!
!
!
!
200 CONTINUE
CALL ATOMSK_MSG(4057,(/vfile/),(/0.d0/))
!Read the file containing parameters for Voronoi construction
CALL CHECKFILE(vfile,'read')
OPEN(UNIT=31,FILE=vfile)
!Parse the file a first time to count number of nodes
Nnodes = 0
Hset=.FALSE.
paramnode = .FALSE.
paramrand = .FALSE.
paramlatt = .FALSE.
linenumber = 0
DO
  READ(31,'(a)',END=210,ERR=210) line
  line = TRIM(ADJUSTL(line))
  linenumber = linenumber+1
  !
  !Ignore empty lines and lines starting with #
  IF( line(1:1).NE."#" .AND. LEN_TRIM(line)>0 ) THEN
    IF( line(1:3)=="box" ) THEN
      !Read size of the final box
      READ(line(4:),*,END=830,ERR=830) P1, P2, P3
      !Set final box vectors
      H(:,:) = 0.d0
      H(1,1) = DABS(P1)
      H(2,2) = DABS(P2)
      H(3,3) = DABS(P3)
      !
      IF( H(1,1)<1.1d0*Huc(1,1) .OR. H(2,2)<1.1d0*Huc(2,2) .OR. H(3,3)<1.1d0*Huc(3,3) ) THEN
        !The user asked for a final cell with at least one small dimension (smaller than seed)
        !=> Look along which dimension the final box is "small"
        twodim=0  !=1 if cell is small along X; =2 along Y; =3 along Z
        m=0  !counter for how many dimensions are small (only one is allowed to be small)
        DO i=1,3
          IF( i==1 ) msg="X"
          IF( i==2 ) msg="Y"
          IF( i==3 ) msg="Z"
          IF( VECLENGTH(H(i,:)) < 1.1d0*VECLENGTH(Huc(i,:)) ) THEN
            !The final box is zero along this dimension
            m=m+1
            IF( m>=2 ) THEN
              !Final box is small in many directions => error
              nerr=nerr+1
              CALL ATOMSK_MSG(4828,(/""/),(/0.d0/))
              GOTO 1000
            ELSE
              !The final box is small in only one dimension => pseudo-2D system
              IF( VECLENGTH(H(i,:)) < 1.d-12 ) THEN
                !Cell size was actually zero along that dimension => warn that it will be resized
                nwarn=nwarn+1
                CALL ATOMSK_MSG(4714,(/TRIM(msg)/),(/Huc(i,i)/))
              ENDIF
              twodim = i
              !Make sure that the final box dimension matches the seed dimension in that direction
              H(i,i) = Huc(i,i)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      WRITE(msg,*) "twodim = ", twodim
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      Hset=.TRUE.
      !
      !Check that there will be enough memory for final system
      IF(ALLOCATED(Q)) DEALLOCATE(Q)
      CALL VOLUME_PARA(Huc,Vmin)  ! Volume of seed
      CALL VOLUME_PARA(H,Volume)  ! Volume of final box
      m = CEILING( SIZE(Puc,1) * Volume/Vmin )  !estimate of number of atoms in final box
      IF( m<=0 ) THEN
        ! if m is negative, it's probably because it exceeded the accuracy of INTEGER
        nerr = nerr+1
        CALL ATOMSK_MSG(819,(/''/),(/0.d0/))
        GOTO 1000
      ELSE
        ALLOCATE(Q(m,4),STAT=i)
        IF( i>0 ) THEN
          ! Allocation failed (not enough memory)
          nerr = nerr+1
          CALL ATOMSK_MSG(819,(/''/),(/0.d0/))
          GOTO 1000
        ENDIF
        IF(ALLOCATED(Q)) DEALLOCATE(Q)
      ENDIF
      !
    ELSEIF( line(1:5)=="node " .OR. line(1:5)=="grain" ) THEN
      !Check that the box was defined
      IF( .NOT.Hset ) THEN
        GOTO 820
      ENDIF
      Nnodes = Nnodes+1
      paramnode = .TRUE.
      !
    ELSEIF( line(1:7)=="lattice" ) THEN
      !Check that the box was defined
      IF( .NOT.Hset ) THEN
        GOTO 820
      ENDIF
      paramlatt = .TRUE.
      lattice = TRIM(ADJUSTL(line(8:)))
      !Set the number of nodes according to the lattice type
      !Beware of pseudo-2D systems
      IF( lattice=="sc" ) THEN
        Nnodes = 1
      ELSEIF( lattice=="bcc" ) THEN
        Nnodes = 2
      ELSEIF( lattice=="fcc" ) THEN
        IF(twodim>0) THEN
          Nnodes = 3
        ELSE
          Nnodes = 4
        ENDIF
      ELSEIF( lattice=="diamond" .OR. lattice=="dia"  ) THEN
        IF(twodim>0) THEN
          Nnodes = 6
        ELSE
          Nnodes = 8
        ENDIF
      ELSEIF( lattice=="hcp" ) THEN
        IF(twodim>0) THEN
          Nnodes = 2
        ELSE
          Nnodes = 4
        ENDIF
      ELSE
        !Unrecognized lattice type => abort
        nerr=nerr+1
        GOTO 1000
      ENDIF
      !
    ELSEIF( line(1:6)=="random" ) THEN
      !Check that the box was defined
      IF( .NOT.Hset ) THEN
        GOTO 820
      ENDIF
      !Read total number of grains
      READ(line(7:),*,END=800,ERR=800) Nnodes
      paramrand = .TRUE.
      !
    ELSE
      !Unknown command => display warning
      nwarn=nwarn+1
      temp(1:64) = CHARLONG2SHRT(vfile)
      CALL ATOMSK_MSG(1702,(/line,temp/),(/0.d0/))
    ENDIF
    !
  ENDIF
ENDDO
!
210 CONTINUE
IF( verbosity==4 ) THEN
  WRITE(msg,'(a5,3(f9.3,a3))') "Box: ", H(1,1), " x ", H(2,2), " x ",  H(3,3)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "Nnodes = ", Nnodes
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!Keywords "lattice", "node" and "random" are mutually exclusive
!If two of them appear in the param file, display an error message and exit
IF( paramnode .AND. paramrand ) THEN
  CALL ATOMSK_MSG(4823,(/"node   ","random "/),(/0.d0/))
  nerr=nerr+1
  GOTO 1000
ELSEIF( paramnode .AND. paramlatt ) THEN
  CALL ATOMSK_MSG(4823,(/"node   ","lattice"/),(/0.d0/))
  nerr=nerr+1
  GOTO 1000
ELSEIF( paramrand .AND. paramlatt ) THEN
  CALL ATOMSK_MSG(4823,(/"random ","lattice"/),(/0.d0/))
  nerr=nerr+1
  GOTO 1000
ENDIF
!
ALLOCATE(comment(1))
WRITE(temp,*) Nnodes 
 comment(1) = "# Voronoi polycrystal with "//TRIM(ADJUSTL(temp))//" grains"
!
!Final positions of nodes will be stored in array vnodes(:,:)
!Final crystal orientations of grains will be stored in array vorient(:,:)
ALLOCATE(vnodes(Nnodes,3))
vnodes(:,:) = 0.d0
ALLOCATE(vorient(Nnodes,3,3))
vorient(:,:,:) = 0.d0
!
REWIND(31)
linenumber = 0
Nnodes = 0
DO
  READ(31,'(a)',END=250,ERR=250) line
  line = TRIM(ADJUSTL(line))
  linenumber = linenumber+1
  !Ignore empty lines and lines starting with #
  IF( line(1:1).NE."#" .AND. LEN_TRIM(line)>0 ) THEN
    !
    IF( line(1:3)=="box" ) THEN
      !This part was already dealt with before
      !
      !
    ELSEIF( line(1:7)=="lattice" ) THEN
      !Nodes will be placed according to a pattern
      !Patterns are defined by "lattice", can be fcc, bcc, etc.
      !The crystallographic orientation of the grains will be random
      !
      !Check that the box was defined
      IF( .NOT.Hset ) THEN
        GOTO 820
      ENDIF
      !
      lattice = TRIM(ADJUSTL(line(8:)))
      IF( lattice=="sc" ) THEN
        Nnodes = 1
        vnodes(1,:) = 0.d0
      ELSEIF( lattice=="bcc" ) THEN
        Nnodes = 2
        vnodes(1,:) = 0.d0          !(0,0,0)
        vnodes(2,1) = 0.5d0*H(1,1)  !(1/2,1/2,1/2)
        vnodes(2,2) = 0.5d0*H(2,2)
        vnodes(2,3) = 0.5d0*H(3,3)
        IF(twodim>0) THEN
          vnodes(2,twodim) = 0.d0
        ENDIF
      ELSEIF( lattice=="fcc" ) THEN
        IF(twodim>0) THEN
          !System is 2-D => define only 3 nodes
          Nnodes = 3
          vnodes(1,:) = 0.d0          !(0,0)
          IF(twodim==1) THEN
            vnodes(2,2) = 0.d0*H(2,2)  !(0,0)
            vnodes(2,3) = 0.5d0*H(3,3)
            vnodes(3,2) = 0.5d0*H(2,2) !(1/2,0)
            vnodes(3,3) = 0.d0*H(3,3)
          ELSEIF(twodim==2) THEN
            vnodes(2,1) = 0.d0*H(1,1)  !(0,0)
            vnodes(2,3) = 0.5d0*H(3,3)
            vnodes(3,1) = 0.5d0*H(1,1) !(1/2,0)
            vnodes(3,3) = 0.d0*H(3,3)
          ELSEIF(twodim==3) THEN
            vnodes(2,1) = 0.d0*H(1,1)  !(0,0)
            vnodes(2,2) = 0.5d0*H(2,2)
            vnodes(3,1) = 0.5d0*H(1,1) !(1/2,0)
            vnodes(3,2) = 0.d0*H(2,2)
          ENDIF
          vnodes(:,twodim) = 0.5d0
        ELSE
          !System is 3-D => define fcc lattice
          Nnodes = 4
          vnodes(:,:) = 0.d0          !(0,0,0)
          vnodes(2,1) = 0.5d0*H(1,1)  !(1/2,1/2,0)
          vnodes(2,2) = 0.5d0*H(2,2)
          vnodes(3,1) = 0.5d0*H(1,1)  !(1/2,0,1/2)
          vnodes(3,3) = 0.5d0*H(3,3)
          vnodes(4,2) = 0.5d0*H(2,2)  !(0,1/2,1/2)
          vnodes(4,3) = 0.5d0*H(3,3)
        ENDIF
      ELSEIF( lattice=="diamond" .OR. lattice=="dia" ) THEN
        IF(twodim>0) THEN
          !System is 2-D => define only 6 nodes
          Nnodes = 6
          vnodes(:,:) = 0.d0
          vnodes(2,:) = 0.5d0*H(1,1)  !(1/2,1/2)
          vnodes(2,2) = 0.5d0*H(2,2)
          vnodes(2,3) = 0.5d0*H(3,3)
          vnodes(3,1) = 0.25d0*H(1,1) !(1/4,1/4)
          vnodes(3,2) = 0.25d0*H(2,2)
          vnodes(3,3) = 0.25d0*H(3,3)
          vnodes(4,1) = 0.75d0*H(1,1) !(3/4,3/4)
          vnodes(4,2) = 0.75d0*H(2,2)
          vnodes(4,3) = 0.75d0*H(3,3)
          IF(twodim==1) THEN
            vnodes(5,2) = 0.25d0*H(2,2)  !(1/4,3/4)
            vnodes(5,3) = 0.75d0*H(3,3)
            vnodes(6,2) = 0.75d0*H(2,2)  !(3/4,1/4)
            vnodes(6,3) = 0.25d0*H(3,3)
          ELSEIF(twodim==2) THEN
            vnodes(5,1) = 0.25d0*H(1,1)  !(1/4,3/4)
            vnodes(5,3) = 0.75d0*H(3,3)
            vnodes(6,1) = 0.75d0*H(1,1)  !(3/4,1/4)
            vnodes(6,3) = 0.25d0*H(3,3)
          ELSE   !i.e. twodim==3
            vnodes(5,1) = 0.25d0*H(1,1)  !(1/4,3/4)
            vnodes(5,2) = 0.75d0*H(2,2)
            vnodes(6,1) = 0.75d0*H(1,1)  !(3/4,1/4)
            vnodes(6,2) = 0.25d0*H(2,2)
          ENDIF
          vnodes(:,twodim) = 0.5d0
        ELSE
          !System is 3-D => define diamond lattice
          Nnodes = 8
          vnodes(:,:) = 0.d0
          vnodes(2,1) = 0.5d0*H(1,1)  !(1/2,1/2,0)
          vnodes(2,2) = 0.5d0*H(2,2)
          vnodes(3,1) = 0.5d0*H(1,1)  !(1/2,0,1/2)
          vnodes(3,3) = 0.5d0*H(3,3)
          vnodes(4,2) = 0.5d0*H(2,2)  !(0,1/2,1/2)
          vnodes(4,3) = 0.5d0*H(3,3)
          vnodes(5,1) = 0.25d0*H(1,1) !(1/4,1/4,1/4)
          vnodes(5,2) = 0.25d0*H(2,2)
          vnodes(5,3) = 0.25d0*H(3,3)
          vnodes(6,1) = 0.75d0*H(1,1) !(3/4,3/4,1/4)
          vnodes(6,2) = 0.75d0*H(2,2)
          vnodes(6,3) = 0.25d0*H(3,3)
          vnodes(7,1) = 0.75d0*H(1,1) !(3/4,1/4,3/4)
          vnodes(7,2) = 0.25d0*H(2,2)
          vnodes(7,3) = 0.75d0*H(3,3)
          vnodes(8,1) = 0.25d0*H(1,1) !(1/4,3/4,3/4)
          vnodes(8,2) = 0.75d0*H(2,2)
          vnodes(8,3) = 0.75d0*H(3,3)
        ENDIF
      ELSEIF( lattice=="hcp" ) THEN
        IF(twodim>0) THEN
          !System is 2-D => define only 2 nodes
          Nnodes = 2
          vnodes(:,:) = 0.d0          !(0,0,0)
          vnodes(2,1) = 0.5d0*H(1,1)  !(1/2,1/2,1/2)
          vnodes(2,2) = 0.5d0*H(2,2)
          vnodes(2,3) = 0.5d0*H(3,3)
          vnodes(2,twodim) = 0.d0
        ELSE
          !System is 3-D => define hcp lattice
          Nnodes = 4
          vnodes(:,:) = 0.d0          !(0,0,0)
          vnodes(2,2) = H(2,2)/3.d0   !(0,1/3,1/2)
          vnodes(2,3) = 0.5d0*H(3,3)
          vnodes(3,1) = 0.5d0*H(1,1)  !(1/2,1/2,0)
          vnodes(3,2) = 0.5d0*H(2,2)
          vnodes(4,1) = 0.5d0*H(1,1)  !(1/2,5/6,1/2)
          vnodes(4,2) = 5.d0*H(2,2)/6.d0
          vnodes(4,3) = 0.5d0*H(3,3)
        ENDIF
      ELSE
        !unrecognized lattice type (already dealt with before)
      ENDIF
      !
      !Orientation of each grain will be random
      !Generate a list of 3*Nnodes random numbers
      CALL GEN_NRANDNUMBERS( 3*Nnodes , randarray )
      !randarray now contains 3*Nnodes real numbers between 0 and 1
      !They are used to generate rotation matrices
      IF( twodim>0 ) THEN
        !Only one random number will be used => Multiply all of them by 2*pi and subtract pi
        randarray(:) = randarray(:)*2.d0*pi - pi
      ELSE
        DO i=1,Nnodes
          m = 3*(i-1) + 1
          n = 3*(i-1) + 2
          o = 3*(i-1) + 3
          randarray(m) = randarray(m)*2.d0*pi - pi
          randarray(n) = DACOS(2.d0*randarray(n) - 1.d0)
          randarray(o) = randarray(o)*2.d0*pi - pi
        ENDDO
      ENDIF
      DO i=1,Nnodes
        P1 = randarray(3*(i-1)+1)
        P2 = randarray(3*(i-1)+2)
        P3 = randarray(3*(i-1)+3)
        !
        vorient(i,:,:) = Id_Matrix(:,:)  !unity matrix
        IF( twodim==0 .OR. twodim==1 ) THEN
          !Construct the rotation matrix around X
          rotmat(:,:) = 0.d0
          rotmat(1,1) = 1.d0
          rotmat(2,2) = DCOS(P1)
          rotmat(2,3) = -1.d0*DSIN(P1)
          rotmat(3,2) = DSIN(P1)
          rotmat(3,3) = DCOS(P1)
          vorient(i,:,:) = rotmat(:,:)
        ENDIF
        IF( twodim==0 .OR. twodim==2 ) THEN
          !Construct the rotation matrix around Y
          rotmat(:,:) = 0.d0
          rotmat(2,2) = 1.d0
          rotmat(3,3) = DCOS(P2)
          rotmat(3,1) = -1.d0*DSIN(P2)
          rotmat(1,3) = DSIN(P2)
          rotmat(1,1) = DCOS(P2)
          vorient(i,:,:) = MATMUL( rotmat(:,:) , vorient(i,:,:) )
        ENDIF
        IF( twodim==0 .OR. twodim==3 ) THEN
          !Construct the rotation matrix around Z
          rotmat(:,:) = 0.d0
          rotmat(3,3) = 1.d0
          rotmat(1,1) = DCOS(P3)
          rotmat(1,2) = -1.d0*DSIN(P3)
          rotmat(2,1) = DSIN(P3)
          rotmat(2,2) = DCOS(P3)
          vorient(i,:,:) = MATMUL( rotmat(:,:) , vorient(i,:,:) )
        ENDIF
      ENDDO
      !
      GOTO 250
      !
      !
    ELSEIF( line(1:5)=="node " .OR. line(1:5)=="grain" ) THEN
      !Check that the box was defined
      IF( .NOT.Hset ) THEN
        GOTO 820
      ENDIF
      Nnodes = Nnodes+1
      !
      IF( .NOT.ALLOCATED(randarray) ) THEN
        !Generate N random numbers, just in case one or more nodes have random orientation
        CALL GEN_NRANDNUMBERS( 3*SIZE(vnodes,1) , randarray )
        !randarray now contains N real numbers between 0 and 1
        !Multiply them by 2*pi and subtract pi to generate 3 angles alpha, beta and gamma
        IF( twodim>0 ) THEN
          randarray(:) = randarray(:)*2.d0*pi - pi
        ELSE
          DO i=1,Nnodes
            m = 3*(i-1) + 1
            n = 3*(i-1) + 2
            randarray(m) = randarray(m)*2.d0*pi - pi
            randarray(n) = randarray(n)*2.d0*pi - pi
          ENDDO
        ENDIF
      ENDIF
      !
      !Read position of that grain
      !Note: position may be given with respect to box dimension, e.g. "box/2"
      line = TRIM(ADJUSTL(line(6:)))
      i=SCAN(line,' ')
      temp = line(:i)
      CALL BOX2DBLE( H(1,:),temp,vnodes(Nnodes,1),status )
      IF( status>0 ) GOTO 830
      line = TRIM(ADJUSTL(line(i+1:)))
      i=SCAN(line,' ')
      temp = line(:i)
      CALL BOX2DBLE( H(2,:),temp,vnodes(Nnodes,2),status )
      IF( status>0 ) GOTO 830
      line = TRIM(ADJUSTL(line(i+1:)))
      i=SCAN(line,' ')
      temp = line(:i)
      CALL BOX2DBLE( H(3,:),temp,vnodes(Nnodes,3),status )
      line = TRIM(ADJUSTL(line(i+1:)))
      IF( status>0 ) GOTO 830
      !
      !Read crystallographic orientation of that grain
      !(can be explicitely given as Miller indices, or random)
      IF( line(1:6)=="random" ) THEN
        !Generated random parameters will be written into a file later
        outparam = .TRUE.
        !
        IF( twodim>0 ) THEN
          !Pick a random angle from randarray
          P1 = randarray(3*(Nnodes-1)+1)
          !Get indices
          m = twodim
          n = twodim+1
          IF(n>3) n=n-3
          o = twodim+2
          IF(o>3) o=o-3
          !Construct the rotation matrix around short axis
          rotmat(:,:) = 0.d0
          rotmat(m,m) = 1.d0
          rotmat(n,n) = DCOS(P1)
          rotmat(n,o) = -1.d0*DSIN(P1)
          rotmat(o,n) = DSIN(P1)
          rotmat(o,o) = DCOS(P1)
          !Save rotation matrix in vorient
          vorient(Nnodes,:,:) = rotmat(:,:)
        ELSE
          !Generate a random rotation matrix for this grain, using random numbers generated before
          !Method is from "Fast random rotation matrices", James Arvo, Cornell University
          !NOTE: the distribution of orientations will be completely random only at the
          !     condition that all grain orientations are specified as "random".
          !     If specific angles or Miller indices are given explicitely by the user
          !     for some grains, then the distribution will not be completely random.
          !Compute vector V
          P1 = randarray(3*(Nnodes-1)+3)
          vector(1) = DSQRT(P1)*DCOS(randarray(3*(Nnodes-1)+1))
          vector(2) = DSQRT(P1)*DSIN(randarray(3*(Nnodes-1)+1))
          vector(3) = DSQRT(1.d0-P1)
          !Compute matrix R
          rotmat(:,:) = 0.d0
          rotmat(1,1) = DCOS(randarray(3*(Nnodes-1)+2))
          rotmat(1,2) = DSIN(randarray(3*(Nnodes-1)+2))
          rotmat(2,1) = -1.d0*DSIN(randarray(3*(Nnodes-1)+2))
          rotmat(2,2) = DCOS(randarray(3*(Nnodes-1)+2))
          rotmat(3,3) = 1.d0
          !Compute final rotation matrix:  M = ( 2*V^T*V - I ) * R
          vorient(Nnodes,:,:) = MATMUL( 2.d0*VECMAT(vector,vector) - Id_Matrix , rotmat )
        ENDIF
        !
      ELSE
        !The user provides 3 angles or 3 Miller indices
        READ(line,*,END=830,ERR=830) or1, or2, or3
        miller=.TRUE.
        IF( SCAN(or1,"[")>0 .OR. SCAN(or2,"[")>0 .OR. SCAN(or3,"[")>0 .OR.          &
          & SCAN(or1,"_")>0 .OR. SCAN(or2,"_")>0 .OR. SCAN(or3,"_")>0 ) THEN
          !No ambiguity: it should be Miller indices (given in brackets and/or with underscores)
          miller=.TRUE.
        ELSEIF( SCAN(or1,"°")>0 .OR. SCAN(or2,"°")>0 .OR. SCAN(or3,"°")>0 .OR.      &
              & SCAN(or1,".")>0 .OR. SCAN(or2,".")>0 .OR. SCAN(or3,".")>0 ) THEN
          !No ambiguity: it is angles
          miller=.FALSE.
        ELSE
          !Ambiguous data: the user entered something like "110", is it an angle or Miller indices?
          IF( LEN_TRIM(or1)>=3 .AND. LEN_TRIM(or2)>=3 .AND. LEN_TRIM(or3)>=3 ) THEN
            !Try to interpret it as Miller indices, if it fails then it is angles
            miller=.TRUE.
            CALL INDEX_MILLER(or1,rotmat,j)
            IF(j>0) miller=.FALSE.
            CALL INDEX_MILLER(or2,rotmat,j)
            IF(j>0) miller=.FALSE.
            CALL INDEX_MILLER(or3,rotmat,j)
            IF(j>0) miller=.FALSE.
          ELSE
            !or1, or2 and/or or3 contain only 2 digits => consider they are angles
            miller=.FALSE.
          ENDIF
        ENDIF
        !
        IF( miller ) THEN
          !Read and interpret the Miller indices,
          !save the rotation matrix in vorient(Nnodes,:,:)
          CALL INDEX_MILLER(or1,vorient(Nnodes,1,:),j)
          IF(j>0) GOTO 800
          CALL INDEX_MILLER(or2,vorient(Nnodes,2,:),j)
          IF(j>0) GOTO 800
          CALL INDEX_MILLER(or3,vorient(Nnodes,3,:),j)
          IF(j>0) GOTO 800
        ELSE
          !Read and interpret the angles,
          !save the rotation matrix in vorient(Nnodes,:,:)
          !Rotation will be done in the order X, Y, Z
          !Read the first angle from or1
          j=SCAN(or1,"°")
          IF(j>0) or1(j:j)=" "
          READ(or1,*,END=830,ERR=830) P1
          !Make sure angle lies between -180° and 180°
          DO WHILE( P1<=-180.d0 )
            P1 = P1+360.d0
          ENDDO
          DO WHILE( P1>180.d0 )
            P1 = P1-360.d0
          ENDDO
          !Convert into radians
          P1 = DEG2RAD(P1)
          !Read the second angle from or2
          j=SCAN(or2,"°")
          IF(j>0) or2(j:j)=" "
          READ(or2,*,END=830,ERR=830) P2
          !Make sure angle lies between -180° and 180°
          DO WHILE( P2<=-180.d0 )
            P2 = P2+360.d0
          ENDDO
          DO WHILE( P2>180.d0 )
            P2 = P2-360.d0
          ENDDO
          !Convert into radians
          P2 = DEG2RAD(P2)
          !Read the third angle from or3
          j=SCAN(or3,"°")
          IF(j>0) or3(j:j)=" "
          READ(or3,*,END=830,ERR=830) P3
          !Make sure angle lies between -180° and 180°
          DO WHILE( P3<=-180.d0 )
            P3 = P3+360.d0
          ENDDO
          DO WHILE( P3>180.d0 )
            P3 = P3-360.d0
          ENDDO
          !Convert into radians
          P3 = DEG2RAD(P3)
          !Construct the rotation matrix around X
          rotmat(:,:) = 0.d0
          rotmat(1,1) = 1.d0
          rotmat(2,2) = DCOS(P1)
          rotmat(2,3) = -1.d0*DSIN(P1)
          rotmat(3,2) = DSIN(P1)
          rotmat(3,3) = DCOS(P1)
          vorient(Nnodes,:,:) = rotmat(:,:)
          !Construct the rotation matrix around Y
          rotmat(:,:) = 0.d0
          rotmat(2,2) = 1.d0
          rotmat(3,3) = DCOS(P2)
          rotmat(3,1) = -1.d0*DSIN(P2)
          rotmat(1,3) = DSIN(P2)
          rotmat(1,1) = DCOS(P2)
          vorient(Nnodes,:,:) = MATMUL( rotmat(:,:) , vorient(Nnodes,:,:) )
          !Construct the rotation matrix around Z
          rotmat(:,:) = 0.d0
          rotmat(3,3) = 1.d0
          rotmat(1,1) = DCOS(P3)
          rotmat(1,2) = -1.d0*DSIN(P3)
          rotmat(2,1) = DSIN(P3)
          rotmat(2,2) = DCOS(P3)
          vorient(Nnodes,:,:) = MATMUL( rotmat(:,:) , vorient(Nnodes,:,:) )
        ENDIF
      ENDIF
      !
      !
    ELSEIF( line(1:6)=="random" ) THEN
      !Position and orientations of grains are random
      !Generated random parameters will be written into a file later
      outparam = .TRUE.
      !Check that the box was defined
      IF( .NOT.Hset ) THEN
        GOTO 820
      ENDIF
      !Read total number of grains
      READ(line(7:),*,END=830,ERR=830) Nnodes
      IF(Nnodes<1) THEN
        CALL ATOMSK_MSG(4831,(/vfile/),(/0.d0/))
        nerr=nerr+1
        GOTO 1000
      ENDIF
      !Generate a list of 6*Nnodes random numbers
      CALL GEN_NRANDNUMBERS( 6*Nnodes , randarray )
      !randarray now contains 6*Nnodes real numbers between 0 and 1
      !The 3*Nnodes first random numbers are used to generate positions
      DO i=1,Nnodes
        vnodes(i,1) = randarray(3*(i-1)+1) * H(1,1)
        vnodes(i,2) = randarray(3*(i-1)+2) * H(2,2)
        vnodes(i,3) = randarray(3*(i-1)+3) * H(3,3)
      ENDDO
      !
      !The last 3*Nnodes random numbers are used to generate rotation matrices
      !Modify them to generate angles
      IF( twodim>0 ) THEN
        !Grains are rotated only around one axis
        !Multiply all random numbers by 2*pi to generate random angles
        randarray(3*Nnodes:) = randarray(3*Nnodes:)*2.d0*pi - pi
        DO i=1,Nnodes
          !Pick a random angle
          P1 = randarray(3*Nnodes+3*(i-1)+1)
          !Get indices
          m = twodim
          n = twodim+1
          IF(n>3) n=n-3
          o = twodim+2
          IF(o>3) o=o-3
          !Construct the rotation matrix around short axis
          rotmat(:,:) = 0.d0
          rotmat(m,m) = 1.d0
          rotmat(n,n) = DCOS(P1)
          rotmat(n,o) = -1.d0*DSIN(P1)
          rotmat(o,n) = DSIN(P1)
          rotmat(o,o) = DCOS(P1)
          !Save rotation matrix in vorient
          vorient(i,:,:) = rotmat(:,:)
        ENDDO
      ELSE
        !Generate a random rotation matrix for each grain
        !Method is from "Fast random rotation matrices", James Arvo, Cornell University
        DO i=1,Nnodes
          !Generate two random angles; third random number will be used below
          m = 3*Nnodes + 3*(i-1) + 1
          n = 3*Nnodes + 3*(i-1) + 2
          randarray(m) = randarray(m)*2.d0*pi - pi
          randarray(n) = randarray(n)*2.d0*pi - pi
          !Compute vector V
          P1 = randarray(3*Nnodes+3*(i-1)+3)
          vector(1) = DSQRT(P1)*DCOS(randarray(3*Nnodes+3*(i-1)+1))
          vector(2) = DSQRT(P1)*DSIN(randarray(3*Nnodes+3*(i-1)+1))
          vector(3) = DSQRT(1.d0-P1)
          !Compute matrix R
          rotmat(:,:) = 0.d0
          rotmat(1,1) = DCOS(randarray(3*Nnodes+3*(i-1)+2))
          rotmat(1,2) = DSIN(randarray(3*Nnodes+3*(i-1)+2))
          rotmat(2,1) = -1.d0*DSIN(randarray(3*Nnodes+3*(i-1)+2))
          rotmat(2,2) = DCOS(randarray(3*Nnodes+3*(i-1)+2))
          rotmat(3,3) = 1.d0
          !Compute final rotation matrix:  M = ( 2*V^T*V - I ) * R
          vorient(i,:,:) = MATMUL( 2.d0*VECMAT(vector,vector) - Id_Matrix , rotmat )
        ENDDO
      ENDIF
      !
      !
      IF( verbosity==4 ) THEN
        !Write angles into a file for visualization/debug purposes
        !Rotation matrices are applied to vector [100]
        !File is in XYZ format to be visualized with VESTA, gnuplot or other softwares
        OPEN(UNIT=43,FILE="atomsk_angles.xyz",STATUS="UNKNOWN")
        WRITE(43,*) Nnodes
        WRITE(43,*) "# Distribution of random angles generated by Atomsk"
        DO i=1,Nnodes
          vector(1) = vorient(i,1,1)
          vector(2) = vorient(i,2,1)
          vector(3) = vorient(i,3,1)
          WRITE(43,'(a3,3(f16.3,2X))') "1  ", vector(1), vector(2), vector(3)
        ENDDO
        CLOSE(43)
      ENDIF
      !
      GOTO 250
      !
      !
    ENDIF
  ENDIF
ENDDO
250 CONTINUE
CLOSE(31)
!
!Check that number of nodes is not zero
IF(Nnodes<1) THEN
  CALL ATOMSK_MSG(4831,(/vfile/),(/0.d0/))
  nerr=nerr+1
  GOTO 1000
ENDIF
!
!Make sure that all nodes are inside the final box H(:,:)
CALL CART2FRAC(vnodes,H)
k = 0
DO i=1,SIZE(vnodes,1) !loop on all nodes
  m=0
  DO j=1,3  !loop on xyz
    DO WHILE( vnodes(i,j)>=1.d0 )
      vnodes(i,j) = vnodes(i,j)-1.d0
      m=m+1
    ENDDO
    DO WHILE( vnodes(i,j)<0.d0 )
      vnodes(i,j) = vnodes(i,j)+1.d0
      m=m+1
    ENDDO
  ENDDO
  IF(m>0) k=k+1
ENDDO
CALL FRAC2CART(vnodes,H)
IF(k>0) THEN
  !Some nodes were wrapped: display message
  nwarn=nwarn+1
  CALL ATOMSK_MSG(4716,(/''/),(/DBLE(k)/))
ENDIF
!
CALL ATOMSK_MSG(4058,(/''/),(/DBLE(Nnodes),DBLE(twodim)/))
!
!
!For 2-D polycrystal, check that the user rotated the grains
!only around the rotation axis
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
ELSE
  !Otherwise (3-D case), check if all nodes are on the same plane
  IF( SIZE(vnodes,1)>1 ) THEN
    sameplane = .TRUE.
    DO j=1,3 !loop on xyz
      DO i=2,SIZE(vnodes,1)
        IF( DABS(vnodes(i,j)-vnodes(1,j)) > 0.1d0 ) THEN
          sameplane = .FALSE.
        ENDIF
      ENDDO
    ENDDO
    WRITE(msg,*) "All nodes belong to the same plane: ", sameplane
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
ENDIF
!
IF( outparam .AND. ofu.NE.6 ) THEN
  !Write positions and orientations in a parameter file
  OPEN(41,FILE=outparamfile,STATUS="UNKNOWN")
  WRITE(41,'(a62)') "# Random positions and rotations of grains generated by Atomsk"
  WRITE(41,'(a66)') "# This parameter file can be used to generate the same polycrystal"
  WRITE(41,'(a4,3f16.6)') "box ", H(1,1), H(2,2), H(3,3)
  DO i=1,SIZE(vorient,1)
    !Compute corresponding rotation vectors and write them into parameter file
    P1 = DATAN2( vorient(i,3,2) , vorient(i,3,3) )
    P2 = DATAN2( -1.d0*vorient(i,3,1) , DSQRT(vorient(i,3,2)**2 + vorient(i,3,3)**2) )
    P3 = DATAN2( vorient(i,2,1) , vorient(i,1,1) )
    WRITE(41,'(a5,6f16.6)') "node ", vnodes(i,:), RAD2DEG(P1), RAD2DEG(P2), RAD2DEG(P3)
  ENDDO
  CLOSE(41)
ENDIF
!
!The maximum number of faces of any polyhedron should be 11 in 2-D,
!and 59 in 3-D, for proof see e.g.:
!   http://www.ericharshbarger.org/voronoi.html
!To that we will add 26 self-neighbors in 3-D
!This will be used to limit the number of iterations in the loops on jnode below
IF( twodim > 0 ) THEN
  maxvertex = 21
ELSE
  maxvertex = 59
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
  DO i=1,SIZE(vnodes,1)
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
!Compute boxmax = maximum distance from one end of the box to another
boxmax = 1.2d0*VECLENGTH( (/ H(1,1) , H(2,2) , H(3,3) /)  ) + 2.d0
WRITE(msg,*) "Max. distance for neighbor search:", boxmax
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!
300 CONTINUE
!Construct a neighbor list of nodes
CALL NEIGHBOR_LIST(H,vnodes,boxmax,vnodesNeighList)
IF(nerr>0) GOTO 1000
!WARNING: if user asks for only one grain, then vnodesNeighList is NOT ALLOCATED!
IF(verbosity==4) THEN
  !Debug messages
  IF( ALLOCATED(vnodesNeighList) .AND. SIZE(vnodesNeighList,1)>1 ) THEN
    msg = "Neighbor list for nodes:"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO i=1,SIZE(vnodesNeighList,1)
      WRITE(msg,'(32i4)') i, (vnodesNeighList(i,j),j=1,MIN(SIZE(vnodesNeighList,2),20))
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
  ENDIF
ENDIF
!
!Allocate array containing number of atoms in each grain
ALLOCATE( NPgrains(Nnodes) , STAT=i )
IF( i>0 ) THEN
  ! Allocation failed (not enough memory)
  nerr = nerr+1
  CALL ATOMSK_MSG(819,(/''/),(/DBLE(Nnodes)/))
  GOTO 1000
ENDIF
NPgrains(:) = 0
!
!Construct template supercell Pt(:,:)
!By default the template grain is a bit larger than max.cell size * sqrt(3)
!This template grain will be cut later to construct each grain
expandmatrix(:) = 1
WRITE(msg,*) "Determining expansion factors:"
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
DO i=1,3
  IF( VECLENGTH(Huc(i,:)) < 0.8d0*VECLENGTH(H(1,:)) .OR.  &
    & VECLENGTH(Huc(i,:)) < 0.8d0*VECLENGTH(H(2,:)) .OR.  &
    & VECLENGTH(Huc(i,:)) < 0.8d0*VECLENGTH(H(3,:))       ) THEN
    !
    P1 = CEILING( MAX( VECLENGTH(H(1,:))/VECLENGTH(Huc(:,i)) , &
               & VECLENGTH(H(2,:))/VECLENGTH(Huc(:,i)) , VECLENGTH(H(3,:))/VECLENGTH(Huc(:,i)) ) )
    !
    WRITE(msg,'(a11,i1,a4,i6)') "    expand(", i, ") = ", NINT(P1)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !If the number of grains is small, the template grain may not be large enough
    IF( sameplane ) THEN
      !Special case: all nodes are in the same plane
      P1 = 1.8d0*P1
    ELSEIF( Nnodes<=4 ) THEN
      !Relatively few nodes: increase size to ensure it covers all grains
      P1 = 1.2d0*P1
    ENDIF
    !Make sure duplication factors are not crazy
    IF(P1==0) THEN
      P1 = 1
    ELSEIF(P1==2) THEN
      P1=3
    ELSEIF(P1>2000) THEN
      P1=1999
    ENDIF
    expandmatrix(i) = NINT(P1)
  ENDIF
ENDDO
WRITE(msg,*) "Initial expansion factors:", expandmatrix(:)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!If the system is 2-D, do not expand along the shortest axis
IF( twodim>0 ) THEN
  expandmatrix(twodim) = 0
ENDIF
!Compute how many particles the template will contain
P1 = DBLE(MAX(1,expandmatrix(1))) * DBLE(MAX(1,expandmatrix(2))) * DBLE(MAX(1,expandmatrix(3))) * DBLE(SIZE(Puc,1))
WRITE(msg,'(a25,f18.0)') "Expected NP for template:", P1
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!If expected number of particles is very large, reduce some values in expandmatrix(:)
IF( P1 > 2.147d9 ) THEN
  WRITE(msg,'(a25,f18.0)') "NP is too large, reducing expansion factors...:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  IF( Nnodes>8 ) THEN
    !Many nodes in the box => reduce drastically the size of template grain
    expandmatrix(:) = NINT(0.7d0 * expandmatrix(:))
  ELSE
    !There are not many grains => do not reduce too much
    expandmatrix(:) = NINT( 0.9d0 * expandmatrix(:) )
  ENDIF
ENDIF
WRITE(msg,'(a47,3i6)') "Final corrected expansion factors for template:", expandmatrix(:)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!Re-calculate expected number of atoms
P1 = DBLE(MAX(1,expandmatrix(1))) * DBLE(MAX(1,expandmatrix(2))) * DBLE(MAX(1,expandmatrix(3))) * DBLE(SIZE(Puc,1))
!If expected NP is too large, abort completely
IF( P1 > 2.147d9 ) THEN
  CALL ATOMSK_MSG(821,(/""/),(/P1/))
  nerr = nerr+1
  GOTO 1000
ENDIF
!
! Allocate array Pt for template grain (full duplicated crystal, not oriented or truncated)
m = PRODUCT(expandmatrix(:)+1)*SIZE(Puc,1)
WRITE(msg,*) "ALLOCATE  Pt, SIZE = ", m
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ALLOCATE( Pt(m,4) , STAT=i )
IF( i>0 ) THEN
  ! Allocation failed (not enough memory)
  nerr = nerr+1
  CALL ATOMSK_MSG(819,(/''/),(/DBLE(m)/))
  GOTO 1000
ENDIF
Pt(:,:) = 0.d0
IF( doshells ) THEN
  ALLOCATE( St( SIZE(Pt,1) , 4 ) , STAT=i )
  IF( i>0 ) THEN
    ! Allocation failed (not enough memory)
    nerr = nerr+1
    CALL ATOMSK_MSG(819,(/''/),(/DBLE(m)/))
    GOTO 1000
  ENDIF
  St(:,:) = 0.d0
ENDIF
IF( doaux ) THEN
  ALLOCATE( AUX_Q( SIZE(Pt,1) , SIZE(AUXuc,2)+1 ) , STAT=i )
  IF( i>0 ) THEN
    ! Allocation failed (not enough memory)
    nerr = nerr+1
    CALL ATOMSK_MSG(819,(/''/),(/DBLE(m)/))
    GOTO 1000
  ENDIF
ELSE
  ALLOCATE( AUX_Q( SIZE(Pt,1) , 1 ) , STAT=i )
  IF( i>0 ) THEN
    ! Allocation failed (not enough memory)
    nerr = nerr+1
    CALL ATOMSK_MSG(819,(/''/),(/DBLE(m)/))
    GOTO 1000
  ENDIF
ENDIF
AUX_Q(:,:) = 0.d0
!
! Fill array Pt(:,:)
qi=0
DO o = 0 , expandmatrix(3)
  DO n = 0 , expandmatrix(2)
    DO m = 0 , expandmatrix(1)
      DO i=1,SIZE(Puc,1)
        !qi = i + ( m + n*expandmatrix(1) + o*expandmatrix(2)*expandmatrix(1) ) * SIZE(Puc,1)
        qi = qi+1
        !Compute (cartesian) position of the replica of this atom
        Pt(qi,1) = Puc(i,1) + DBLE(m)*Huc(1,1) + DBLE(n)*Huc(2,1) + DBLE(o)*Huc(3,1)
        Pt(qi,2) = Puc(i,2) + DBLE(m)*Huc(1,2) + DBLE(n)*Huc(2,2) + DBLE(o)*Huc(3,2)
        Pt(qi,3) = Puc(i,3) + DBLE(m)*Huc(1,3) + DBLE(n)*Huc(2,3) + DBLE(o)*Huc(3,3)
        Pt(qi,4) = Puc(i,4)
        IF(doshells) THEN
          !Compute (cartesian) position of the replica of this shell
          St(qi,1) = Suc(i,1) + DBLE(m)*Huc(1,1) + DBLE(n)*Huc(2,1) + DBLE(o)*Huc(3,1)
          St(qi,2) = Suc(i,2) + DBLE(m)*Huc(1,2) + DBLE(n)*Huc(2,2) + DBLE(o)*Huc(3,2)
          St(qi,3) = Suc(i,3) + DBLE(m)*Huc(1,3) + DBLE(n)*Huc(2,3) + DBLE(o)*Huc(3,3)
          St(qi,4) = Puc(i,4)
        ENDIF
        IF(doaux) THEN
          AUX_Q(qi,:) = AUXuc(i,:)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
ENDDO
!
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
!
!Estimate new number of particles NP = (density of unit cell) / (volume of final cell)
NP = CEILING( SIZE(Puc,1) * DABS( DABS(H(1,1)*H(2,2)*H(3,3)) / &
      & DABS(VECLENGTH(Huc(1,:))*VECLENGTH(Huc(2,:))*VECLENGTH(Huc(3,:))) ) )
WRITE(msg,*) "Estimated number of atoms in polycrystal: NP = ", NP
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!Allow for +50% and +1000 atoms to allocate arrays.
!Actual size of arrays will be adjusted later
NP = NINT(1.74d0*NP) + 1000
WRITE(msg,*) "Allocating arrays with size: NP = ", NP
ALLOCATE(Q(NP,4))
Q(:,:) = 0.d0
IF(doshells) THEN
  ALLOCATE(T(NP,4))
  T(:,:) = 0.d0
ENDIF
ALLOCATE(newAUX(NP,SIZE(AUXNAMES)))
newAUX(:,:) = 0.d0
!
!
!Construct grains using Voronoi tesselation
NP=0 !counter for atoms that will make it in the final polycrystal
DO inode=1,Nnodes
  !This is grain # inode
  CALL ATOMSK_MSG(4055,(/''/),(/DBLE(inode)/))
  !
  !Compute number of vertices surrounding current node
  Nvertices = 1000
  maxdnodes = 0.d0
  expandmatrix(:) = 1
  !
  !Search neighbors from neighbor list
  IF( twodim>0 ) THEN
    !System is pseudo 2-D => do not look along the short distance
    expandmatrix(twodim) = 0
  ENDIF
  !
  !Allocate memory for vertices
  !NOTE: at this stage Nvertices=1000 is expected to over-estimate the number
  !      of neighboring vertices. This array will be resized later
  ALLOCATE( vvertex(Nvertices,4) )
  vvertex(:,:) = 0.d0
  !
  Nvertices = 0
  IF( ALLOCATED(vnodesNeighList) .AND. SIZE(vnodesNeighList,1)>1 ) THEN
    !
    !Loop on all neighboring nodes
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(m,n,o,P1,P2,P3,vector,distance,i,j,k,jnode)
    DO i=1,SIZE(vnodesNeighList,2)
      !Save index of the i-th neighbor as jnode
      jnode = vnodesNeighList(inode,i)
      !Check if jnode is a valid node index
      IF( jnode>0 ) THEN
        !Node #jnode is neighbor of node #inode
        !Check which periodic image(s) are actually neighbors
        DO o=-expandmatrix(3),expandmatrix(3)
          DO n=-expandmatrix(2),expandmatrix(2)
            DO m=-expandmatrix(1),expandmatrix(1)
              !Position of the periodic image of neighbor #jnode along this Cartesian axis
              P1 = vnodes(jnode,1) + DBLE(m)*H(1,1) + DBLE(n)*H(2,1) + DBLE(o)*H(3,1)
              P2 = vnodes(jnode,2) + DBLE(m)*H(1,2) + DBLE(n)*H(2,2) + DBLE(o)*H(3,2)
              P3 = vnodes(jnode,3) + DBLE(m)*H(1,3) + DBLE(n)*H(2,3) + DBLE(o)*H(3,3)
              vector = (/ P1 , P2 , P3 /)
              !Compute distance between current node and this vertex
              distance = VECLENGTH( vector(:) - vnodes(inode,:) )
              !This image is a neighbor if distance is smaller than max. box size
              IF( distance>1.d-3 .AND. distance <= boxmax ) THEN
                !$OMP CRITICAL
                Nvertices = Nvertices+1
                k = Nvertices
                !$OMP END CRITICAL
                IF( k>SIZE(vvertex,1) ) THEN
                  !Increase size of array vvertex
                  CALL RESIZE_DBLEARRAY2(vvertex,k+10,4)
                ENDIF
                !Save vertex position = middle point between nodes #inode and #jnode
                vvertex(k,1:3) = vnodes(inode,:) + (vector(:)-vnodes(inode,:))/2.d0
                vvertex(k,4) = distance
                !
              ENDIF  !end if distance<boxmax
              !
            ENDDO  ! end loop on m
          ENDDO    ! end loop on n
        ENDDO      ! end loop on o
        !
      ENDIF
    ENDDO
    !$OMP END PARALLEL DO
    !
  ENDIF
  !
  WRITE(msg,'(a,i6)') "N vertices for this grain:", Nvertices
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Fill the rest of the list with very large distances before sorting.
  !This way, after sorting they will be at the end of the list
  !(those artefacts will be removed later)
  IF( Nvertices<SIZE(vvertex,1) ) THEN
    DO i=Nvertices+1,SIZE(vvertex,1)
      vvertex(i,4) = 1.d12
    ENDDO
  ENDIF
  !
  !Sort vertices by increasing distance
  CALL BUBBLESORT(vvertex,4,"up  ",newindex)
  !
  !All neighboring vertices will not be used, only the maxvertex first ones
  !If total number of neighboring vertices is greater than maxvertex, then correct maxdnodes
  IF( SIZE(vvertex,1)>maxvertex ) THEN
    maxdnodes = vvertex(MIN(maxvertex,Nvertices),4)
    WRITE(msg,*) maxdnodes
    msg = "Keep vertices only up to max.distance: "//TRIM(ADJUSTL(msg))
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
  !
  !Get index of last vertex that is closer than maxdnodes
  Nvertices=0
  DO WHILE( Nvertices<=SIZE(vvertex,1) .AND. Nvertices<=maxvertex .AND. vvertex(Nvertices+1,4)<=maxdnodes )
    Nvertices = Nvertices+1
  ENDDO
  !Resize array vvertex to get rid of unused vertices
  IF(twodim>0) THEN
    CALL RESIZE_DBLEARRAY2(vvertex,Nvertices+8,4,status)
  ELSE
    CALL RESIZE_DBLEARRAY2(vvertex,Nvertices+26,4,status)
  ENDIF
  !
  !Append vertices corresponding to the replicas of current node
  WRITE(msg,'(a,i6)') "Adding self neighbors (=periodic replica of current grain)"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO o=-expandmatrix(3),expandmatrix(3)
    DO n=-expandmatrix(2),expandmatrix(2)
      DO m=-expandmatrix(1),expandmatrix(1)
        IF( o.NE.0 .OR. n.NE.0 .OR. m.NE.0 ) THEN
          !Position of the periodic image of node #inode
          P1 = vnodes(inode,1) + DBLE(m)*H(1,1) + DBLE(n)*H(2,1) + DBLE(o)*H(3,1)
          P2 = vnodes(inode,2) + DBLE(m)*H(1,2) + DBLE(n)*H(2,2) + DBLE(o)*H(3,2)
          P3 = vnodes(inode,3) + DBLE(m)*H(1,3) + DBLE(n)*H(2,3) + DBLE(o)*H(3,3)
          vector = (/ P1 , P2 , P3 /)
          !Compute distance
          distance = VECLENGTH( vector(:) - vnodes(inode,:) )
          IF( distance>1.d-3 ) THEN
            Nvertices = Nvertices+1
            vvertex(Nvertices,1:3) = vnodes(inode,:) + (vector(:)-vnodes(inode,:))/2.d0
            vvertex(Nvertices,4) = distance
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  Nvertexmax = Nvertices
  !
  IF(verbosity==4) THEN
    !Debug messages
    msg = "List of neighboring vertices after cleanup:"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    msg = "   N       x           y          z        dist.to current node"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO i=1,SIZE(vvertex,1)
      WRITE(msg,'(i4,6f12.3)') i, vvertex(i,:)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
    WRITE(msg,*) inode
    OPEN(UNIT=40,FILE="atomsk_grain"//TRIM(ADJUSTL(msg))//"_vertices.xyz")
    WRITE(40,*) SIZE(vvertex,1)+1
    WRITE(40,*) "# Position of node # "//TRIM(ADJUSTL(msg))//" and its vertices"
    WRITE(40,*) 2, vnodes(inode,:)
    DO i=1,SIZE(vvertex,1)
      WRITE(40,'(i4,6f12.3)') 1, vvertex(i,:)
    ENDDO
    CLOSE(40)
  ENDIF
  !
  !Copy template supercell Pt(:,:) into Pt2(:,:)
  !Convert to reduced coordinates
  Ht(:,:) = H(:,:)
  IF(ALLOCATED(Pt2)) DEALLOCATE(Pt2)
  IF(ALLOCATED(St2)) DEALLOCATE(St2)
  ALLOCATE( Pt2( SIZE(Pt,1) , SIZE(Pt,2) ) , STAT=n )
  IF( n>0 ) THEN
    ! Allocation failed (not enough memory)
    nerr = nerr+1
    CALL ATOMSK_MSG(819,(/''/),(/DBLE(2*SIZE(Pt,1))/))
    GOTO 1000
  ENDIF
  Pt2(:,:) = Pt(:,:)
  CALL CART2FRAC(Pt2,Ht)
  IF(doshells) THEN
    !Copy shells from template into St2(:,:)
    ALLOCATE( St2( SIZE(St,1) , SIZE(St,2) ) , STAT=n )
    IF( n>0 ) THEN
      ! Allocation failed (not enough memory)
      nerr = nerr+1
      CALL ATOMSK_MSG(819,(/''/),(/DBLE(3*SIZE(Pt,1))/))
      GOTO 1000
    ENDIF
    St2(:,:) = St(:,:)
    CALL CART2FRAC(St2,Ht)
  ENDIF
  !
  !Rotate this cell to obtain the desired crystallographic orientation
  CALL ATOMSK_MSG(2071,(/''/),(/0.d0/))
  DO i=1,3
    IF( VECLENGTH(Ht(i,:)) > 1.d-12 ) THEN
      Hn(i,:) = Ht(i,:)/VECLENGTH(Ht(i,:))
    ELSE
      !we have a problem
      nerr = nerr+1
      GOTO 1000
    ENDIF
    IF( VECLENGTH(vorient(inode,i,:)) > 1.d-12 ) THEN
      Hend(i,:) = vorient(inode,i,:)/VECLENGTH(vorient(inode,i,:))
    ELSE
      !we have a problem
      nerr = nerr+1
      GOTO 1000
    ENDIF
  ENDDO
  CALL INVMAT(Hn,Gt)
  DO i=1,3
    P1 = Ht(i,1)
    P2 = Ht(i,2)
    P3 = Ht(i,3)
    Ht(i,1) = Gt(1,1)*P1 + Gt(1,2)*P2 + Gt(1,3)*P3
    Ht(i,2) = Gt(2,1)*P1 + Gt(2,2)*P2 + Gt(2,3)*P3
    Ht(i,3) = Gt(3,1)*P1 + Gt(3,2)*P2 + Gt(3,3)*P3
    P1 = Ht(i,1)
    P2 = Ht(i,2)
    P3 = Ht(i,3)
    Ht(i,1) = P1*Hend(1,1) + P2*Hend(1,2) + P3*Hend(1,3)
    Ht(i,2) = P1*Hend(2,1) + P2*Hend(2,2) + P3*Hend(2,3)
    Ht(i,3) = P1*Hend(3,1) + P2*Hend(3,2) + P3*Hend(3,3)
  ENDDO
  !
  !Convert back to Cartesian coordinates
  CALL FRAC2CART(Pt2,Ht)
  IF(doshells) THEN
    CALL FRAC2CART(St2,Ht)
  ENDIF
  CALL ATOMSK_MSG(2072,(/''/),(/0.d0/))
  !
  !Shift oriented supercell so that its center of mass is at the center of the box
  CALL CENTER_XYZ(H,Pt2,St2,0,SELECT)
  !
  !Get center of grain = barycenter of the closest vertices
  !The actual node itself is given a "weight" in this calculation, but
  !the center of the grain will be different from the position of the node itself
  GrainCenter(:) = 0.d0
  n=MIN( SIZE(vvertex,1) , 6 )
  DO jnode=1,n
    GrainCenter(:) = GrainCenter(:) + vvertex(jnode,:)
  ENDDO
  GrainCenter(:) = ( 2.d0*vnodes(inode,1:3) + GrainCenter(:) ) / DBLE(n+2)
  IF( verbosity==4 ) THEN
    WRITE(msg,'(a15,3f9.3)') 'Node position: ', vnodes(inode,1:3)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,'(a15,3f9.3)') 'Grain center:  ', GrainCenter(:)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
  !
  qi=NP !save current number of particles
  !For each atom of the template Pt2, find out if it is located inside the grain
  !If so, then save it to the array Q; if not, discard it (change its atomic number to -1)
  !$OMP PARALLEL DO DEFAULT(SHARED) &
  !$OMP& PRIVATE(i,isinpolyhedron,jnode,vnormal,vector) &
  !$OMP& REDUCTION(+:NP,NPgrains)
  DO i=1,SIZE(Pt2,1)
    !Shift oriented supercell so that its center of mass is at the position of the node
    Pt2(i,1:3) = Pt2(i,1:3) - 0.5d0*(/H(:,1)+H(:,2)+H(:,3)/) + GrainCenter(1:3)
    !Determine if this position is inside of the polyhedron
    isinpolyhedron = .TRUE.
    DO jnode = 1 , SIZE(vvertex,1)
      !Compute vector between the vertex and current node
      !By definition this vector is normal to the grain boundary
      vnormal(:) = vvertex(jnode,:) - vnodes(inode,:)
      !Compute vector between atom and current node
      vector(:) = Pt2(i,1:3) - vnodes(inode,:)
      IF( VEC_PLANE(vnormal,VECLENGTH(vnormal),vector) > -1.d-12 ) THEN
        !Atom is above this plane of cut, hence out of the polyhedron
        !=> exit the loop on jnode
        isinpolyhedron = .FALSE.
        EXIT
      ENDIF
    ENDDO
    !
    IF( isinpolyhedron ) THEN
      !Atom is inside the polyhedron
      !Increment number of atoms belonging to grain #inode
      NPgrains(inode) = NPgrains(inode)+1
      !Increment total number of particles in the system
      NP = NP+1
    ELSE
      !Atom is outside of the polyhedron -> mark it for termination
      Pt2(i,4) = -1.d0
    ENDIF
    !
  ENDDO
  !$OMP END PARALLEL DO
  !
  DO i=1,SIZE(Pt2,1)
    IF( Pt2(i,4) > 0.1d0 ) THEN
      qi = qi+1
      Q(qi,:) = Pt2(i,:)
      IF(doshells) T(qi,:) = St2(i,:)
      IF(doaux) THEN
        DO j=1,SIZE(AUXNAMES)-1
          newAUX(qi,j) = AUX_Q(i,j)
        ENDDO
      ENDIF
      newAUX(qi,grainID) = DBLE(inode)
    ENDIF
  ENDDO
  !
  IF(nerr>0) GOTO 1000
  !
  CALL ATOMSK_MSG(4056,(/''/),(/DBLE(NPgrains(inode))/))
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
      DO i=1,SIZE(Pt2,1)
        IF( Pt2(i,4)>0.d0 ) THEN
          CALL ATOMSPECIES(Pt2(i,4),species)
          WRITE(36,'(a2,2X,3(f16.8,1X))') species, Pt2(i,1:3)
        ENDIF
      ENDDO
    ENDIF
    WRITE(36,'(a4)') 'alat'
    WRITE(36,'(a3)') '1.0'
    WRITE(36,'(a9)') 'supercell'
    WRITE(36,'(3f16.6)') H(1,:)
    WRITE(36,'(3f16.6)') H(2,:)
    WRITE(36,'(3f16.6)') H(3,:)
    CLOSE(36)
    WRITE(msg,*) "Wrote atom positions"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,*) "   o  o  o  o  o  o  o  o  o  o  o  o  o  o  o"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
  !
  IF(ALLOCATED(vvertex)) DEALLOCATE(vvertex)
  IF(ALLOCATED(Pt2)) DEALLOCATE(Pt2)
  IF(ALLOCATED(St2)) DEALLOCATE(St2)
  !
ENDDO  !end loop on inode
!
!
!
500 CONTINUE
!NP is now the actual number of atoms in the final polycrystal
!Q now contains positions of all atoms in all the grains, but may be oversized
!(and T contains the positions of shells, and newAUX the aux.prop. if relevant)
!Copy atom positions into final array P with appropriate size
IF(ALLOCATED(P)) DEALLOCATE(P)
ALLOCATE(P(NP,4))
DO i=1,NP
  P(i,:) = Q(i,:)
ENDDO
DEALLOCATE(Q)
!Copy shell positions into final array S
IF(doshells) THEN
  IF(ALLOCATED(S)) DEALLOCATE(S)
  ALLOCATE(S(NP,4))
  DO i=1,NP
    S(i,:) = T(i,:)
  ENDDO
  DEALLOCATE(T)
ENDIF
!Copy aux.prop. into final array AUX
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
ALLOCATE(AUX(NP,SIZE(newAUX,2)))
DO i=1,NP
  AUX(i,:) = newAUX(i,:)
ENDDO
DEALLOCATE(newAUX)
!
!Apply options to the final system
CALL OPTIONS_AFF(options_array,Huc,H,P,S,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
!
!Write final system to file(s)
CALL WRITE_AFF(prefix,outfileformats,H,P,S,comment,AUXNAMES,AUX)
!
!
!
600 CONTINUE
!From now on, some arrays are not needed anymore => deallocate them
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(S)) DEALLOCATE(S)
IF(ALLOCATED(T)) DEALLOCATE(T)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF(ALLOCATED(comment)) DEALLOCATE(comment)
IF(ALLOCATED(outfileformats)) DEALLOCATE(outfileformats)
!
IF( outparam .AND. ofu.NE.6 ) THEN
  !Random parameters were written into a file => inform user
  temp = "parameter"
  CALL ATOMSK_MSG(3002,(/temp,outparamfile/),(/0.d0/))
ENDIF
!
IF( ofu.NE.6 ) THEN
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
Vmin = HUGE(1.d0)
Vmax = 0.d0
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
  !
  !Estimate volume occupied by the grain
  CALL VOLUME_PARA(Huc,Volume)
  Volume = DBLE(NPgrains(i)) * Volume / DBLE(SIZE(Puc,1))
  IF( Volume > Vmax ) THEN
    Vmax = Volume
  ENDIF
  IF( Volume < Vmin ) THEN
    Vmin = Volume
  ENDIF
ENDDO
!
IF( ofu.NE.6 ) THEN
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
    WRITE(41,*) "# grainID   ;   N.atoms   ;   Grain volume (A^3)"
    DO i=1,SIZE(NPgrains)  !loop on all grains
      !Estimate volume occupied by that grain
      CALL VOLUME_PARA(Huc,Volume)
      Volume = DBLE(NPgrains(i)) * Volume / DBLE(SIZE(Puc,1))
      WRITE(41,*) i, NPgrains(i), Volume
    ENDDO
    CLOSE(41)
    !
    !Compute the grain size distribution and write it to a file
    !with the format:  ( Grain size ; Number of grains with that size )
    Vstep = (Vmax-Vmin)/20.d0  !step for grain size distribution
    IF(.NOT.overw) CALL CHECKFILE(distfile,'writ')
    OPEN(UNIT=41,FILE=distfile,FORM="FORMATTED",STATUS="UNKNOWN")
    WRITE(41,*) "# Grain size distribution: volume (A^3) ; No. of grains"
    DO j=0,20  !loop on grain size
      Nnodes = 0
      DO i=1,SIZE(NPgrains) !loop on all grains
        !Estimate volume occupied by that grain
        CALL VOLUME_PARA(Huc,Volume)
        Volume = DBLE(NPgrains(i)) * Volume / DBLE(SIZE(Puc,1))
        IF( Volume >= Vmin+DBLE(j)*Vstep .AND. Volume < Vmin+DBLE(j+1)*Vstep ) THEN
          !This grain has the appropriate size => increment counter
          Nnodes = Nnodes+1
        ENDIF
      ENDDO
      !Write to file
      WRITE(41,*) Vmin+DBLE(j)*Vstep, Nnodes
    ENDDO
    CLOSE(41)
    msg = "DAT"
    CALL ATOMSK_MSG(3002,(/msg,distfile/),(/0.d0/))
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
IF(ALLOCATED(P)) DEALLOCATE(P)
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(S)) DEALLOCATE(S)
IF(ALLOCATED(T)) DEALLOCATE(T)
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
