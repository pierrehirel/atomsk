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
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 22 Sept. 2015                                    *
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
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE,INTENT(IN):: outfileformats !list of formats to output
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
!
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: distfile   !name of file containing grain size distribution
CHARACTER(LEN=128):: line
CHARACTER(LEN=128):: or1, or2, or3
CHARACTER(LEN=128):: lattice  !if grains are organized according to a lattice
CHARACTER(LEN=128):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES    !names of auxiliary properties of atoms
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: newAUXNAMES !names of auxiliary properties of atoms (temporary)
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: doshells, doaux !are there shells, auxiliary properties in initial seed?
LOGICAL:: Hset            !are the box vectors H(:,:) defined?
LOGICAL:: isinpolyhedron  !is atom inside the polyhedron?
LOGICAL:: miller          !are Miller indices given? (if no then angles are given)
LOGICAL:: paramnode, paramrand, paramlatt !are the keywords "node", "random", "lattice" used in parameter file?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT
INTEGER:: twodim        !=0 if system is 3-D, =1,2,3 if system is thin along x, y, z
INTEGER:: grainID       !position of the auxiliary property "grainID" in AUX
INTEGER:: i, j
INTEGER:: m, n, o
INTEGER:: maxvertex   !max. number of vertices to look for
INTEGER:: NP   !total number of atoms in the final system
INTEGER:: qi   !used to count atoms in a grain
INTEGER:: inode, jnode
INTEGER:: Nnodes, Nvertices !number of nodes, of vertices
INTEGER:: status
INTEGER,DIMENSION(:),ALLOCATABLE:: NPgrains  !number of atoms in each grain
INTEGER,DIMENSION(3):: expandmatrix
INTEGER,DIMENSION(:,:),ALLOCATABLE:: vnodesNeighList  !list of neighbours for nodes
REAL(dp):: boxmax      !max. distance from one end of the box to another
REAL(dp):: distance, distance2    !distance between two points
REAL(dp):: maxdnodes   !maximum distance between 2 nodes
REAL(dp):: P1, P2, P3  !temporary position
REAL(dp):: Volume, Vmin, Vmax, Vstep  !min, max. volume occupied by a grain, step for grain size distribution
REAL(dp),DIMENSION(3):: GrainCenter !position of center of grain (different from node position)
REAL(dp),DIMENSION(3):: vector    !vector between an atom and a node
REAL(dp),DIMENSION(3):: vnormal   !vector normal to grain boundary
REAL(dp),DIMENSION(3,3):: Huc       !Base vectors of the unit cell (seed)
REAL(dp),DIMENSION(3,3):: Ht        !Base vectors of the oriented unit cell
REAL(dp),DIMENSION(3,3):: Hunity    !unit matrix
REAL(dp),DIMENSION(3,3):: H         !Base vectors of the final supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystalographic orientation
REAL(dp),DIMENSION(3,3):: rotmat  !rotation matrix
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:),ALLOCATABLE:: randarray   !random numbers
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Puc, Suc  !positions of atoms, shells in unit cell (seed)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Pt, St    !positions of atoms, shells in template supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, S      !positions of atoms, shells in final supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Q, T      !positions of atoms, shells in a grain
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newP, newQ, newS !positions of atoms, shells (temporary)
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
CALL NAME_OUTFILE(prefix,distfile,"dat  ")
Nnodes = 0
twodim = 0  !assume system will be 3-D
Hunity(:,:) = 0.d0
DO i=1,3
  Hunity(i,i) = 1.d0
ENDDO
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
!     than the final polycrystal, i.e. Huc(:,:) may be larger than H(:,:).
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
DO
  READ(31,'(a)',END=210,ERR=210) line
  line = TRIM(ADJUSTL(line))
  !
  !Ignore empty lines and lines starting with #
  IF( line(1:1).NE."#" .AND. LEN_TRIM(line)>0 ) THEN
    IF( line(1:3)=="box" ) THEN
      !Read size of the final box
      READ(line(4:),*,END=800,ERR=800) P1, P2, P3
      !Set box vectors
      H(:,:) = 0.d0
      H(1,1) = P1
      H(2,2) = P2
      H(3,3) = P3
      IF( SIZE(Puc,1)<100 .OR. ANY(H(:,:)<20.d0) ) THEN
        !The user provided a small cell (<100 atoms), or asked for a final cell with a small dimension (<20A)
        !=> If the final box is smaller than 2 times the unit cell along one dimension,
        !   then consider that it is a 2-D system and use 2-D Voronoi construction
        DO i=1,3
          IF( H(i,i)<2.1d0*VECLENGTH(Huc(i,:)) ) THEN
            !The final box is "small" along this dimension
            twodim = i
            !Make sure that the final box dimension matches the unit cell length
            H(i,i) = Huc(i,i)
          ENDIF
        ENDDO
      ENDIF
      WRITE(msg,*) "twodim = ", twodim
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      Hset=.TRUE.
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
      IF( lattice=="bcc" ) THEN
        Nnodes = 2
      ELSEIF( lattice=="fcc" ) THEN
        IF(twodim>0) THEN
          Nnodes = 2
        ELSE
          Nnodes = 4
        ENDIF
      ELSEIF( lattice=="diamond" ) THEN
        IF(twodim>0) THEN
          Nnodes = 6
        ELSE
          Nnodes = 8
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
Nnodes = 0
DO
  READ(31,'(a)',END=250,ERR=250) line
  line = TRIM(ADJUSTL(line))
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
      IF( lattice=="bcc" ) THEN
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
          !System is 2-D => define only 2 nodes
          Nnodes = 2
          vnodes(1,:) = 0.d0          !(0,0)
          vnodes(2,1) = 0.5d0*H(1,1)  !(1/2,1/2)
          vnodes(2,2) = 0.5d0*H(2,2)
          vnodes(2,3) = 0.5d0*H(3,3)
          vnodes(2,twodim) = 0.d0
        ELSE
          !System is 3-D => define fcc lattice
          Nnodes = 4
          vnodes(1,:) = 0.d0          !(0,0,0)
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
          vnodes(2,1) = 0.5d0*H(1,1)  !(1/2,1/2)
          vnodes(2,2) = 0.5d0*H(2,2)
          vnodes(2,3) = 0.5d0*H(3,3)
          vnodes(2,twodim) = 0.d0
          vnodes(3,1) = 0.25d0*H(1,1) !(1/4,1/4)
          vnodes(3,2) = 0.25d0*H(2,2)
          vnodes(3,3) = 0.25d0*H(3,3)
          vnodes(3,twodim) = 0.d0
          vnodes(4,1) = 0.75d0*H(1,1) !(3/4,3/4)
          vnodes(4,2) = 0.75d0*H(2,2)
          vnodes(4,3) = 0.75d0*H(3,3)
          vnodes(4,twodim) = 0.d0
          IF(twodim==1) THEN
            vnodes(5,2) = 0.25d0*H(2,2)  !(1/4,3/4)
            vnodes(5,3) = 0.75d0*H(3,3)
            vnodes(6,2) = 0.75d0*H(2,2)  !(1/4,3/4)
            vnodes(6,3) = 0.25d0*H(3,3)
          ELSEIF(twodim==2) THEN
            vnodes(5,1) = 0.25d0*H(1,1)  !(1/4,3/4)
            vnodes(5,3) = 0.75d0*H(3,3)
            vnodes(6,1) = 0.75d0*H(1,1)  !(1/4,3/4)
            vnodes(6,3) = 0.25d0*H(3,3)
          ELSE   !i.e. twodim==3
            vnodes(5,1) = 0.25d0*H(1,1)  !(1/4,3/4)
            vnodes(5,2) = 0.75d0*H(2,2)
            vnodes(6,1) = 0.75d0*H(1,1)  !(1/4,3/4)
            vnodes(6,2) = 0.25d0*H(2,2)
          ENDIF
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
      ELSE
        !unrecognized lattice type (already dealt with before)
      ENDIF
      !
      !Orientation of each grain will be random
      !Generate a list of 3*Nnodes random numbers
      CALL GEN_NRANDNUMBERS( 3*Nnodes , randarray )
      !randarray now contains 3*Nnodes real numbers between 0 and 1
      !They are used to generate rotation matrices
      !Multiply them by 2*pi and subtract pi to generate 3 angles alpha, beta and gamma
      randarray(:) = randarray(:)*2.d0*pi - pi
      DO i=1,Nnodes
        P1 = randarray(3*(i-1)+1)
        P2 = randarray(3*(i-1)+3)
        P3 = randarray(3*(i-1)+3)
        !
        vorient(i,:,:) = Hunity(:,:)
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
      !Read position of that grain
      !Note: position may be given with respect to box dimension, e.g. "box/2"
      line = TRIM(ADJUSTL(line(6:)))
      i=SCAN(line,' ')
      temp = line(:i)
      CALL BOX2DBLE( H(1,:),temp,vnodes(Nnodes,1),status )
      IF( status>0 ) GOTO 810
      line = TRIM(ADJUSTL(line(i+1:)))
      i=SCAN(line,' ')
      temp = line(:i)
      CALL BOX2DBLE( H(2,:),temp,vnodes(Nnodes,2),status )
      IF( status>0 ) GOTO 810
      line = TRIM(ADJUSTL(line(i+1:)))
      i=SCAN(line,' ')
      temp = line(:i)
      CALL BOX2DBLE( H(3,:),temp,vnodes(Nnodes,3),status )
      line = TRIM(ADJUSTL(line(i+1:)))
      IF( status>0 ) GOTO 810
      !
      !Read crystallographic orientation of that grain
      !(can be explicitely given as Miller indices, or random)
      IF( line(1:6)=="random" ) THEN
        !Generate 3 random numbers
        CALL GEN_NRANDNUMBERS( 3 , randarray )
        !
        !randarray now contains 3 real numbers between 0 and 1
        !Multiply them by 2*pi and subtract pi to generate 3 angles alpha, beta and gamma
        randarray(:) = randarray(:)*2.d0*pi - pi
        !
        vorient(Nnodes,:,:) = Hunity(:,:)
        !
        IF( twodim==0 .OR. twodim==1 ) THEN
          !Construct the rotation matrix around X
          rotmat(:,:) = 0.d0
          rotmat(1,1) = 1.d0
          rotmat(2,2) = DCOS(randarray(1))
          rotmat(2,3) = -1.d0*DSIN(randarray(1))
          rotmat(3,2) = DSIN(randarray(1))
          rotmat(3,3) = DCOS(randarray(1))
          vorient(Nnodes,:,:) = rotmat(:,:)
        ENDIF
        IF( twodim==0 .OR. twodim==2 ) THEN
          !Construct the rotation matrix around Y
          rotmat(:,:) = 0.d0
          rotmat(2,2) = 1.d0
          rotmat(3,3) = DCOS(randarray(2))
          rotmat(3,1) = -1.d0*DSIN(randarray(2))
          rotmat(1,3) = DSIN(randarray(2))
          rotmat(1,1) = DCOS(randarray(2))
          vorient(Nnodes,:,:) = MATMUL( rotmat(:,:) , vorient(Nnodes,:,:) )
        ENDIF
        IF( twodim==0 .OR. twodim==3 ) THEN
          !Construct the rotation matrix around Z
          rotmat(:,:) = 0.d0
          rotmat(3,3) = 1.d0
          rotmat(1,1) = DCOS(randarray(3))
          rotmat(1,2) = -1.d0*DSIN(randarray(3))
          rotmat(2,1) = DSIN(randarray(3))
          rotmat(2,2) = DCOS(randarray(3))
          vorient(Nnodes,:,:) = MATMUL( rotmat(:,:) , vorient(Nnodes,:,:) )
        ENDIF
        !
      ELSE
        !The user provides 3 angles or 3 Miller indices
        READ(line,*,END=800,ERR=800) or1, or2, or3
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
          !Try to interpret it as Miller indices, if it fails then it is angles
          miller=.TRUE.
          CALL INDEX_MILLER(or1,rotmat,j)
          IF(j>0) miller=.FALSE.
          CALL INDEX_MILLER(or2,rotmat,j)
          IF(j>0) miller=.FALSE.
          CALL INDEX_MILLER(or3,rotmat,j)
          IF(j>0) miller=.FALSE.
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
          !Read the three angles from or1, or2, or3
          j=SCAN(or1,"°")
          IF(j>0) or1(j:j)=" "
          READ(or1,*,END=800,ERR=800) P1
          P1 = DEG2RAD(P1)
          j=SCAN(or2,"°")
          IF(j>0) or2(j:j)=" "
          READ(or2,*,END=800,ERR=800) P2
          P2 = DEG2RAD(P2)
          j=SCAN(or3,"°")
          IF(j>0) or3(j:j)=" "
          READ(or3,*,END=800,ERR=800) P3
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
      !Check that the box was defined
      IF( .NOT.Hset ) THEN
        GOTO 820
      ENDIF
      !Read total number of grains
      READ(line(7:),*,END=800,ERR=800) Nnodes
      IF(Nnodes<1) GOTO 800
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
      !Multiply them by 2*pi and subtract pi to generate 3 angles alpha, beta and gamma
      randarray(3*Nnodes:) = randarray(3*Nnodes:)*2.d0*pi - pi
      DO i=1,Nnodes
        P1 = randarray(3*Nnodes+3*(i-1)+1)
        P2 = randarray(3*Nnodes+3*(i-1)+2)
        P3 = randarray(3*Nnodes+3*(i-1)+3)
        !
        vorient(i,:,:) = Hunity(:,:)
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
    ENDIF
  ENDIF
ENDDO
250 CONTINUE
CLOSE(31)
!
CALL ATOMSK_MSG(4058,(/''/),(/DBLE(Nnodes),DBLE(twodim)/))
!
!If system is pseudo 2-D then place all nodes on the same plane
IF( twodim==1 ) THEN
  vnodes(:,1) = 0.d0
ELSEIF( twodim==2 ) THEN
  vnodes(:,2) = 0.d0
ELSEIF( twodim==3 ) THEN
  vnodes(:,3) = 0.d0
ENDIF
!
!The maximum number of faces of any polyhedron should be 11 in 2-D,
!and 59 in 3-D, for proof see e.g.:
!   http://www.ericharshbarger.org/voronoi.html
!This will be used to limit the number of iterations in the loops on jnode below
IF( twodim > 0 ) THEN
  maxvertex = 20
ELSE
  maxvertex = 65
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
  msg = "Orientation of nodes:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(vnodes,1)
    WRITE(msg,'(i4,a1,3f9.3,a1)') i, '[', (vorient(i,1,j),j=1,3), ']'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,'(i4,a1,3f9.3,a1)') i, '[', (vorient(i,2,j),j=1,3), ']'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,'(i4,a1,3f9.3,a1)') i, '[', (vorient(i,3,j),j=1,3), ']'
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
boxmax = 1.05d0*VECLENGTH( (/ H(1,1) , H(2,2) , H(3,3) /)  )
!
!
!
300 CONTINUE
!Construct a neighbor list of nodes
CALL NEIGHBOR_LIST(H,vnodes,boxmax,vnodesNeighList)
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
ALLOCATE( NPgrains(Nnodes) )
NPgrains(:) = 0
!
!Construct template supercell Pt(:,:)
!By default the template grain fills the whole box
!This template grain will be cut later to construct each grain
IF( VECLENGTH(Huc(1,:)) < VECLENGTH(H(1,:)) .OR.  &
  & VECLENGTH(Huc(1,:)) < VECLENGTH(H(2,:)) .OR.  &
  & VECLENGTH(Huc(1,:)) < VECLENGTH(H(3,:)) .OR.  &
  & VECLENGTH(Huc(2,:)) < VECLENGTH(H(1,:)) .OR.  &
  & VECLENGTH(Huc(2,:)) < VECLENGTH(H(2,:)) .OR.  &
  & VECLENGTH(Huc(2,:)) < VECLENGTH(H(3,:)) .OR.  &
  & VECLENGTH(Huc(3,:)) < VECLENGTH(H(1,:)) .OR.  &
  & VECLENGTH(Huc(3,:)) < VECLENGTH(H(2,:)) .OR.  &
  & VECLENGTH(Huc(3,:)) < VECLENGTH(H(3,:))       ) THEN
  !
  DO i=1,3
    expandmatrix(i) = CEILING( 1.1d0*MAX( VECLENGTH(H(1,:))/Huc(i,i) , &
                    & VECLENGTH(H(2,:))/Huc(i,i) , VECLENGTH(H(3,:))/Huc(i,i) ) )
  ENDDO
  !If the number of grains is small, the template grain may not be large enough
  IF( Nnodes<=6 ) THEN
    DO i=1,3
      expandmatrix(i) = NINT( 1.5d0*DBLE(expandmatrix(i)) )
    ENDDO
  ENDIF
  !
ELSE
  !All dimensions of the seed provided by the user are larger
  !than any dimension of the box => don't expand it
  expandmatrix(:) = 1
ENDIF
!If the system is 2-D, do not expand along the shortest axis
IF( twodim>0 ) THEN
  expandmatrix(twodim) = 1
ENDIF
!Evaluate how many particles the template will contain
m = PRODUCT(expandmatrix(:))*SIZE(Puc,1)
WRITE(msg,'(a25,i18)') "Expected NP for template:", m
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!If m is very large, reduce some value in expandmatrix(:)
IF( m > 1.d8 ) THEN
  IF( Nnodes>8 ) THEN
    !Many nodes in the box => reduce drastically the size of template grain
    expandmatrix(:) = 0.8d0 * expandmatrix(:)
  ELSE
    !There are not many grains => do not reduce too much
    expandmatrix(:) = 0.8d0 * expandmatrix(:)
  ENDIF
ENDIF
!
WRITE(msg,'(a32,3i4)') "Creating template grain, expand:", expandmatrix(:)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
ALLOCATE( Pt( PRODUCT(expandmatrix(:))*SIZE(Puc,1) , 4 ) )
Pt(:,:) = 0.d0
IF( doshells ) THEN
  ALLOCATE( St( SIZE(Pt,1) , 4 ) )
  St(:,:) = 0.d0
ENDIF
IF( doaux ) THEN
  ALLOCATE( AUX_Q( SIZE(Pt,1) , SIZE(AUXuc,2)+1 ) )
ELSE
  ALLOCATE( AUX_Q( SIZE(Pt,1) , 1 ) )
ENDIF
AUX_Q(:,:) = 0.d0
WRITE(msg,'(a47,3i3)') "Creating template grain, expand:", expandmatrix(:)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
qi=0
DO o = 1 , expandmatrix(3)
  DO n = 1 , expandmatrix(2)
    DO m = 1 , expandmatrix(1)
      DO i=1,SIZE(Puc,1)
        qi=qi+1
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
!Construct grains using Voronoi tesselation
NP=0
DO inode=1,Nnodes
  !This is grain # inode
  CALL ATOMSK_MSG(4055,(/''/),(/DBLE(inode)/))
  !
  IF(ALLOCATED(Q)) DEALLOCATE(Q)
  IF(ALLOCATED(T)) DEALLOCATE(T)
  !
  !Compute number of vertices surrounding current node
  Nvertices = 0
  maxdnodes = 0.d0
  !First, check if current node is a neighbor of itself due to periodic boundary conditions
  expandmatrix(:) = 1
  IF( twodim>0 ) THEN
    !System is pseudo 2-D => do not look along the short distance
    expandmatrix(twodim) = 0
  ENDIF
  DO o=-expandmatrix(3),expandmatrix(3)
    DO n=-expandmatrix(2),expandmatrix(2)
      DO m=-expandmatrix(1),expandmatrix(1)
        IF( .NOT.(m==0 .AND. n==0 .AND. o==0) .AND.           &
          & (m.NE.n .OR. n.NE.o .OR. o.NE.m)         ) THEN
          !Position of the periodic image of node #inode
          P1 = vnodes(inode,1) + DBLE(m)*H(1,1) + DBLE(n)*H(2,1) + DBLE(o)*H(3,1)
          P2 = vnodes(inode,2) + DBLE(m)*H(1,2) + DBLE(n)*H(2,2) + DBLE(o)*H(3,2)
          P3 = vnodes(inode,3) + DBLE(m)*H(1,3) + DBLE(n)*H(2,3) + DBLE(o)*H(3,3)
          vector = (/ P1 , P2 , P3 /)
          !This image is a neighbor if distance is smaller than max. box size
          distance = VECLENGTH( vector(:) - vnodes(inode,:) )
          IF( distance>1.d-12 .AND. distance <= boxmax ) THEN
            Nvertices = Nvertices+1
            IF( distance > maxdnodes ) maxdnodes = distance
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  !Then search neighbor list
  IF( ALLOCATED(vnodesNeighList) .AND. SIZE(vnodesNeighList,1)>1 ) THEN
    DO i=1,SIZE(vnodesNeighList,2)  !loop on neighboring nodes
      IF( vnodesNeighList(inode,i).NE.0 ) THEN
        !This node is neighbor of node #inode
        !Check which periodic image(s) are actually neighbor
        DO o=-expandmatrix(3),expandmatrix(3)
          DO n=-expandmatrix(2),expandmatrix(2)
            DO m=-expandmatrix(1),expandmatrix(1)
              !Position of the periodic image of the i-th neighboring node of node #inode
              P1 = vnodes(vnodesNeighList(inode,i),1) + DBLE(m)*H(1,1) + DBLE(n)*H(2,1) + DBLE(o)*H(3,1)
              P2 = vnodes(vnodesNeighList(inode,i),2) + DBLE(m)*H(1,2) + DBLE(n)*H(2,2) + DBLE(o)*H(3,2)
              P3 = vnodes(vnodesNeighList(inode,i),3) + DBLE(m)*H(1,3) + DBLE(n)*H(2,3) + DBLE(o)*H(3,3)
              vector = (/ P1 , P2 , P3 /)
              !This image is a neighbor if distance is smaller than max. box size
              distance = VECLENGTH( vector(:) - vnodes(inode,:) )
              IF( distance>1.d-12 .AND. distance <= boxmax ) THEN
                Nvertices = Nvertices+1
                IF( distance > maxdnodes ) maxdnodes = distance
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !Allocate memory for vertices
  ALLOCATE( vvertex(Nvertices,4) )
  vvertex(:,:) = 0.d0
  !Save the positions of neighbor vertices for current node
  Nvertices = 0
  DO o=-expandmatrix(3),expandmatrix(3)
    DO n=-expandmatrix(2),expandmatrix(2)
      DO m=-expandmatrix(1),expandmatrix(1)
        IF( .NOT.(m==0 .AND. n==0 .AND. o==0) .AND.           &
          & (m.NE.n .OR. n.NE.o .OR. o.NE.m)         ) THEN
          !Position of the periodic image of node #inode
          P1 = vnodes(inode,1) + DBLE(m)*H(1,1) + DBLE(n)*H(2,1) + DBLE(o)*H(3,1)
          P2 = vnodes(inode,2) + DBLE(m)*H(1,2) + DBLE(n)*H(2,2) + DBLE(o)*H(3,2)
          P3 = vnodes(inode,3) + DBLE(m)*H(1,3) + DBLE(n)*H(2,3) + DBLE(o)*H(3,3)
          vector = (/ P1 , P2 , P3 /)
          !This image is a neighbor if distance is smaller than max. box size
          distance = VECLENGTH( vector(:) - vnodes(inode,:) )
          IF( distance>1.d-12 .AND. distance <= boxmax ) THEN
            Nvertices = Nvertices+1
            vvertex(Nvertices,1:3) = vnodes(inode,:) + (vector(:)-vnodes(inode,:))/2.d0
            vvertex(Nvertices,4) = distance
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  IF( ALLOCATED(vnodesNeighList) .AND. SIZE(vnodesNeighList,1)>1 ) THEN
    DO i=1,SIZE( vnodesNeighList,2 )
      IF( vnodesNeighList(inode,i).NE.0 ) THEN
        !This node is neighbor of node #inode
        !Check which periodic image(s) are actually neighbor
        DO o=-expandmatrix(3),expandmatrix(3)
          DO n=-expandmatrix(2),expandmatrix(2)
            DO m=-expandmatrix(1),expandmatrix(1)
              !Position of the periodic image of the i-th neighboring node of node #inode
              P1 = vnodes(vnodesNeighList(inode,i),1) + DBLE(m)*H(1,1) + DBLE(n)*H(2,1) + DBLE(o)*H(3,1)
              P2 = vnodes(vnodesNeighList(inode,i),2) + DBLE(m)*H(1,2) + DBLE(n)*H(2,2) + DBLE(o)*H(3,2)
              P3 = vnodes(vnodesNeighList(inode,i),3) + DBLE(m)*H(1,3) + DBLE(n)*H(2,3) + DBLE(o)*H(3,3)
              vector = (/ P1 , P2 , P3 /)
              !This image is a neighbor if distance is smaller than max. box size
              distance = VECLENGTH( vector(:) - vnodes(inode,:) )
              IF( distance>1.d-12 .AND. distance <= boxmax ) THEN
                Nvertices = Nvertices+1
                vvertex(Nvertices,1:3) = vnodes(inode,:) + (vector(:)-vnodes(inode,:))/2.d0
                vvertex(Nvertices,4) = distance
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !
  WRITE(msg,'(a,i6)') "N neighbors for this grain:", Nvertices
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Sort vertices by increasing distance
  CALL BUBBLESORT(vvertex,4,"up  ")
  !
  !All neighboring vertices will not be used, only the maxvertex first ones
  !If total number of neighboring vertices is greater than maxvertex, then correct maxdnodes
  IF( SIZE(vvertex,1)>maxvertex ) THEN
    maxdnodes = vvertex(maxvertex,4)
  ENDIF
  !
  IF(verbosity==4) THEN
    !Debug messages
    msg = "List of neighboring vertices:"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO i=1,SIZE(vvertex,1)
      WRITE(msg,'(i4,6f12.3)') i, vvertex(i,:)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
    WRITE(msg,*) maxdnodes
    msg = "Max. distance to a node: "//TRIM(ADJUSTL(msg))
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
  !
  !Copy template supercell Pt(:,:) into Q(:,:)
  Ht(:,:) = H(:,:)
  ALLOCATE( Q( SIZE(Pt,1) , SIZE(Pt,2) ) )
  Q(:,:) = 0.d0
  Q(:,:) = Pt(:,:)
  IF(doshells) THEN
    !Copy shells from template into T(:,:)
    ALLOCATE( T( SIZE(Pt,1) , 4 ) )
    T(:,:) = 0.d0
    T(:,:) = St(:,:)
  ENDIF
  !Auxiliary properties: the array AUX_Q(:,:) will be used
  AUX_Q(:,grainID) = DBLE(inode)
  !
  !Rotate this cell to obtain the desired crystallographic orientation
  CALL ORIENT_XYZ(Ht,Q,T,Ht,vorient(inode,:,:),SELECT,C_tensor)
  IF(nerr>0) GOTO 1000
  !
  !Shift oriented supercell so that its center of mass is at the center of the box
  CALL CENTER_XYZ(H,Q,T,0,SELECT)
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
  qi=0 !so far, zero atom in the grain
  DO i=1,SIZE(Q,1)
    !Shift oriented supercell so that its center of mass is at the position of the node
    Q(i,1:3) = Q(i,1:3) - 0.5d0*(/H(:,1)+H(:,2)+H(:,3)/) + GrainCenter(1:3)
    !Determine if this position is inside of the polyhedron
    isinpolyhedron = .TRUE.
    DO jnode = 1 , MIN( SIZE(vvertex,1) , maxvertex )
      !Compute vector between the vertex and current node
      !By definition this vector is normal to the grain boundary
      vnormal(:) = vvertex(jnode,:) - vnodes(inode,:)
      !Compute vector between atom and current node
      vector(:) = Q(i,1:3) - vnodes(inode,:)
      IF( VEC_PLANE(vnormal,VECLENGTH(vnormal),vector) >= -0.1d0 ) THEN
        !Atom is above this plane of cut, hence out of the polyhedron
        !=> exit the loop on jnode
        isinpolyhedron = .FALSE.
        EXIT
      ENDIF
    ENDDO
    !
    IF( isinpolyhedron ) THEN
      !Atom is inside the polyhedron
      qi = qi+1
    ELSE
      !Atom is outside of the polyhedron -> mark it for termination
      Q(i,4) = -1.d0
    ENDIF
    !
  ENDDO
  !
  CALL ATOMSK_MSG(4056,(/''/),(/DBLE(qi)/))
  !
  !Save number of atoms in this grain
  NPgrains(inode) = qi
  !
  IF( verbosity==4 ) THEN
    !Debug: write positions of atoms of current grain into an XYZ file
    WRITE(temp,*) inode
    msg = 'atomsk_grain'//TRIM(ADJUSTL(temp))//'.xyz'
    OPEN(UNIT=36,FILE=msg,STATUS="UNKNOWN",FORM="FORMATTED")
    WRITE(36,*) MIN(SIZE(Q,1),qi)
    msg = '#Debug file for atomsk, mode --voronoi, grain # '//TRIM(ADJUSTL(temp))
    WRITE(36,*) TRIM(msg)
    DO i=1,MIN(SIZE(Q,1),qi)
      CALL ATOMSPECIES(Q(i,4),species)
      WRITE(36,'(a2,2X,3(f16.8,1X))') species, Q(i,1:3)
    ENDDO
    WRITE(36,'(a4)') 'alat'
    WRITE(36,'(a3)') '1.0'
    WRITE(36,'(a9)') 'supercell'
    WRITE(36,'(3f16.6)') H(1,:)
    WRITE(36,'(3f16.6)') H(2,:)
    WRITE(36,'(3f16.6)') H(3,:)
    CLOSE(36)
    WRITE(msg,*) "Wrote atom positions"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
  !
  IF(ALLOCATED(P)) NP = SIZE(P,1)
  IF( qi>0 ) THEN
    WRITE(msg,*) "Old, new NP: ", NP, NP+qi
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !
    !Add content of array Q into final P
    IF(ALLOCATED(newP)) DEALLOCATE(newP)
    ALLOCATE(newP(NP+qi,4))
    newP(:,:) = 0.d0
    IF( ALLOCATED(P) ) THEN
      DO i=1,SIZE(P,1)
        newP(i,:) = P(i,:)
      ENDDO
    ENDIF
    IF(doshells) THEN
      IF(ALLOCATED(newS)) DEALLOCATE(newS)
      ALLOCATE(newS(NP+qi,4))
      newS(:,:) = 0.d0
      IF( ALLOCATED(S) ) THEN
        DO i=1,SIZE(P,1)
          newS(i,:) = S(i,:)
        ENDDO
      ENDIF
    ENDIF
    !
    IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
    ALLOCATE(newAUX(NP+qi,SIZE(AUX_Q,2)))
    newAUX(:,:) = 0.d0
    IF( ALLOCATED(AUX) ) THEN
      DO i=1,SIZE(AUX,1)
        newAUX(i,:) = AUX(i,:)
      ENDDO
    ENDIF
    !
    n=0
    DO i=1,SIZE(Q,1)
      IF( NINT(Q(i,4))>0 ) THEN
        n=n+1
        newP(NP+n,:) = Q(i,:)
        IF(doshells) THEN
          newS(NP+n,:) = T(i,:)
        ENDIF
        newAUX(NP+n,:) = AUX_Q(i,:)
      ENDIF
    ENDDO
    !
    !Save atom positions of all grains into P
    IF(ALLOCATED(P)) DEALLOCATE(P)
    ALLOCATE( P(SIZE(newP,1),4) )
    P(:,:) = newP(:,:)
    DEALLOCATE(newP)
    DEALLOCATE(Q)
    !
    IF(doshells) THEN
      !Save all shell positions into S
      IF(ALLOCATED(S)) DEALLOCATE(S)
      ALLOCATE( S(SIZE(newS,1),4) )
      S(:,:) = newS(:,:)
      DEALLOCATE(newS)
      DEALLOCATE(T)
    ENDIF
    !
    !Save all auxiliary properties into AUX
    IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
    ALLOCATE( AUX( SIZE(newAUX,1) , SIZE(AUX_Q,2) ) )
    AUX(:,:) = newAUX(:,:)
    DEALLOCATE(newAUX)
    !
    NP = NP+qi
    !
  ELSE  !i.e. zero atoms in this grain, which is strange
    !Display a warning message
    nwarn=nwarn+1
    CALL ATOMSK_MSG(4708,(/""/),(/0.d0/))
  ENDIF
  !
  IF(ALLOCATED(vvertex)) DEALLOCATE(vvertex)
  IF(ALLOCATED(Q)) DEALLOCATE(Q)
  !
ENDDO
!
!
!
500 CONTINUE
IF(ALLOCATED(Pt)) DEALLOCATE(Pt)
IF(ALLOCATED(St)) DEALLOCATE(St)
!P now contains positions of all atoms in all the grains
!Apply options to the final system
CALL OPTIONS_AFF(options_array,H,P,S,AUXNAMES,AUX,ORIENT,SELECT)
!
!Write final system to file(s)
CALL WRITE_AFF(prefix,outfileformats,H,P,S,comment,AUXNAMES,AUX)
!
!Find the max. volume occupied by a grain
Vmin = HUGE(1.d0)
Vmax = 0.d0
DO i=1,SIZE(NPgrains)
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
WRITE(msg,*) "Min, max. volume occupied by a grain: ", Vmin, Vmax
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
IF( Vmax > 0.d0 ) THEN
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
!
!
1000 CONTINUE
IF(ALLOCATED(P)) DEALLOCATE(P)
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(S)) DEALLOCATE(S)
IF(ALLOCATED(T)) DEALLOCATE(T)
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
