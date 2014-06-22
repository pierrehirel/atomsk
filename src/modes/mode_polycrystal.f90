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
!* Last modification: P. Hirel - 15 April 2014                                    *
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
CHARACTER(LEN=*),INTENT(IN):: ucfile  !name of file containing unit cell
CHARACTER(LEN=*),INTENT(IN):: vfile   !name of file containing parameters for Voronoi construction
CHARACTER(LEN=*),INTENT(IN):: prefix  !name or prefix for output file (polycrystal)
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE,INTENT(IN):: outfileformats !list of formats to output
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
!
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: distfile   !name of file containing grain size distribution
CHARACTER(LEN=128):: line
CHARACTER(LEN=128):: or1, or2, or3
CHARACTER(LEN=128):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES   !names of auxiliary properties of atoms
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: doshells, doaux !are there shells, auxiliary properties in initial seed?
LOGICAL:: isinpolyhedron  !is atom inside the polyhedron?
LOGICAL:: miller          !are Miller indices given? (if no then angles are given)
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT
INTEGER:: twodim        !=0 if system is 3-D, =1,2,3 if system is thin along x, y, z
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
REAL(dp),DIMENSION(3):: vector    !vector between an atom and a node
REAL(dp),DIMENSION(3):: vnormal   !vector normal to grain boundary
REAL(dp),DIMENSION(3,3):: Huc       !Base vectors of the unit cell
REAL(dp),DIMENSION(3,3):: Huc_orient !Oriented base vectors of the unit cell
REAL(dp),DIMENSION(3,3):: Hunity    !unit matrix
REAL(dp),DIMENSION(3,3):: H         !Base vectors of the final supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystalographic orientation
REAL(dp),DIMENSION(3,3):: rotmat  !rotation matrix
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:),ALLOCATABLE:: randarray   !random numbers
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Puc, Suc  !positions of atoms, shells in unit cell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Puc_orient, Suc_orient  !positions of atoms, shells in unit cell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, S      !positions of atoms, shells in final supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Q, T      !positions of atoms, shells in a grain
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newP, newQ, newS !positions of atoms, shells (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUXuc     !auxiliary properties of atoms in the unit cell
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
!Read initial seed from file (usually, a unit cell)
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
ELSE
  doaux = .FALSE.
ENDIF
!
!Allocate arrays for oriented seed
ALLOCATE(Puc_orient(SIZE(Puc,1),4))
Puc_orient(:,:) = 0.d0
IF(doshells) THEN
  ALLOCATE(Suc_orient(SIZE(Puc,1),4))
  Suc_orient(:,:) = 0.d0
ENDIF
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
DO
  READ(31,'(a)',END=210,ERR=210) line
  line = TRIM(ADJUSTL(line))
  IF( line(1:5)=="node " .OR. line(1:5)=="grain" ) THEN
    Nnodes = Nnodes+1
  ELSEIF( line(1:6)=="random" ) THEN
    !Read total number of grains
    READ(line(7:),*,END=800,ERR=800) Nnodes
  ENDIF
ENDDO
!
210 CONTINUE
!Store final positions of nodes in array vnodes(:,:)
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
    IF( line(1:3)=="box" ) THEN
      !Read size of the final box
      READ(line(4:),*,END=800,ERR=800) P1, P2, P3
      !Set box vectors
      H(:,:) = 0.d0
      H(1,1) = P1
      H(2,2) = P2
      H(3,3) = P3
      !If the user provided a unit cell (and not a supercell), and
      !if the final box is smaller than 2 times the unit cell along one dimension,
      !then consider that it is a 2-D system and use 2-D Voronoi construction
      IF( .NOT.ANY( H(:,:)>20.d0) .OR. SIZE(Puc,1)<100 ) THEN
        !It is a small unit cell
        IF( H(1,1)<1.1d0*VECLENGTH(Huc(1,:)) .OR. H(1,1)<1.1d0*VECLENGTH(Huc(2,:)) .OR. H(1,1)<1.1d0*VECLENGTH(Huc(3,:)) ) THEN
          twodim = 1
        ELSEIF( H(2,2)<1.1d0*VECLENGTH(Huc(1,:)) .OR. H(2,2)<1.1d0*VECLENGTH(Huc(2,:)) .OR. H(2,2)<1.1d0*VECLENGTH(Huc(3,:)) ) THEN
          twodim = 2
        ELSEIF( H(3,3)<1.1d0*VECLENGTH(Huc(1,:)) .OR. H(3,3)<1.1d0*VECLENGTH(Huc(2,:)) .OR. H(3,3)<1.1d0*VECLENGTH(Huc(3,:)) ) THEN
          twodim = 3
        ENDIF
      ENDIF
      WRITE(msg,*) "twodim = ", twodim
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      !
    ELSEIF( line(1:5)=="node " .OR. line(1:5)=="grain" ) THEN
      !Check that the box was defined
      IF( VECLENGTH(H(1,:))==0.d0 .OR. VECLENGTH(H(2,:))==0.d0 .OR. VECLENGTH(H(3,:))==0.d0 ) THEN
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
    ELSEIF( line(1:6)=="random" ) THEN
      !Position and orientations of grains are random
      !Check that the box was defined
      IF( VECLENGTH(H(1,:))==0.d0 .OR. VECLENGTH(H(2,:))==0.d0 .OR. VECLENGTH(H(3,:))==0.d0 ) THEN
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
  maxvertex = 60
ENDIF
!
IF(verbosity==4) THEN
  !Debug messages
  WRITE(msg,*) "Number of nodes:", Nnodes
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  msg = "Positions of nodes:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(vnodes,1)
    WRITE(msg,*) i, vnodes(i,:)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
  msg = "Orientation of nodes:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(vnodes,1)
    WRITE(msg,*) i, (NINT(vorient(i,1,j)),j=1,3)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,*) i, (NINT(vorient(i,2,j)),j=1,3)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,*) i, (NINT(vorient(i,3,j)),j=1,3)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
ENDIF
!
!Check that the user provided enough information
IF( .NOT. ANY( NINT(H)>0 ) ) THEN
  !User did not provide box, or a box with negative or zero size
  nerr=nerr+1
  !**error message*
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
IF(verbosity==4) THEN
  !Debug messages
  msg = "Neighbor list for nodes:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(vnodesNeighList,1)
    WRITE(msg,'(32i4)') i, (vnodesNeighList(i,j),j=1,MIN(SIZE(vnodesNeighList,2),20))
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
ENDIF
!
!Allocate array containing number of atoms in each grain
ALLOCATE( NPgrains(Nnodes) )
NPgrains(:) = 0
!
!Construct grains using Voronoi tesselation
NP=0
DO inode=1,Nnodes
  !This is grain # inode
  CALL ATOMSK_MSG(4055,(/''/),(/DBLE(inode)/))
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
  !Copy original unit cell to oriented unit cell
  Huc_orient(:,:) = Huc(:,:)
  Puc_orient(:,:) = Puc(:,:)
  IF(doshells) Suc_orient(:,:) = Suc(:,:)
  !
  !Rotated unit cell to obtain the desired crystallographic orientation
  CALL ORIENT_XYZ(Huc_orient,Puc_orient,Suc_orient,Hunity,vorient(inode,:,:),SELECT,C_tensor)
  IF(nerr>0) GOTO 1000
  !
  !Shift oriented unit cell so that 1st atom is at the position of the node
  DO i=1,3
    Puc_orient(:,i) = Puc_orient(:,i) + vnodes(inode,i)
  ENDDO
  IF(doshells) THEN
    DO i=1,3
      Suc_orient(:,i) = Suc_orient(:,i) + vnodes(inode,i)
    ENDDO
  ENDIF
  !
  !Estimate by how much the initial unit cell must be expanded
  !(estimate based on maxdnodes = the largest distance between two nodes)
  !Note: the same factor will be used in the 3 directions of space, because we have no
  !     way to know how the seed will be rotated inside the grain. This should not be a problem
  !     for geometries close to cubic, but may become slow if the seed is very elongated
  !     along one or more directions
  IF( maxdnodes<VECLENGTH(Huc(1,:)) .AND. maxdnodes<VECLENGTH(Huc(2,:)) .AND. maxdnodes<VECLENGTH(Huc(3,:)) ) THEN
    !Seed is bigger than the max. distance between two nodes => no need to expand
    expandmatrix(:) = 0
  ELSE
    expandmatrix(:) = CEILING( 0.5d0*maxdnodes / MIN( VECLENGTH(Huc(1,:)) , VECLENGTH(Huc(2,:)) , VECLENGTH(Huc(3,:)) ) )
    !In each direction Huc_orient(:,:), determine how much the seed must be expanded
!     DO i=1,3
!       DO j=1, maxvertex
!         !Find the closest vertex in that direction
!         distance = VECLENGTH( vvertex(j,:) - vnodes(inode,:) )
!       ENDDO
!     ENDDO
  ENDIF
  !
  IF( twodim>0 ) THEN
    !2-D system => don't expand along the short dimension
    expandmatrix(twodim) = 0
  ENDIF
  msg = "Initial seed will be expanded NxNxN with N="//TRIM(ADJUSTL(msg))
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Expand the seed inside its Voronoi polyhedron,
  !and save positions of atoms of current grain into array Q
  IF(ALLOCATED(Q)) DEALLOCATE(Q)
  !Estimate max. number of atoms in this grain
  !Note: the initial seed will be expanded bewteen -expandmatrix and +expandmatrix
  !     in each direction of space, hence the factor 8
  qi = SIZE(Puc_orient,1)*8*MAX(expandmatrix(1),1)*MAX(expandmatrix(2),1)*MAX(expandmatrix(3),1)
  ALLOCATE( Q( qi , 4 ) )
  Q(:,:) = 0.d0
  IF(doshells) THEN
    IF(ALLOCATED(T)) DEALLOCATE(T)
    ALLOCATE( T( qi , 4 ) )
    T(:,:) = 0.d0
  ENDIF
  IF(doaux) THEN
    IF(ALLOCATED(AUX_Q)) DEALLOCATE(AUX_Q)
    ALLOCATE( AUX_Q( qi , SIZE(AUXuc,2) ) )
    AUX_Q(:,:) = 0.d0
  ENDIF
  !
  qi=0
  DO o = -expandmatrix(3) , expandmatrix(3)
    DO n = -expandmatrix(2) , expandmatrix(2)
      DO m = -expandmatrix(1) , expandmatrix(1)
        DO i=1,SIZE(Puc_orient,1)
          !
          !Compute (cartesian) position of the replica of this atom
          P1 = Puc_orient(i,1) + DBLE(m)*Huc_orient(1,1) + DBLE(n)*Huc_orient(2,1) + DBLE(o)*Huc_orient(3,1)
          P2 = Puc_orient(i,2) + DBLE(m)*Huc_orient(1,2) + DBLE(n)*Huc_orient(2,2) + DBLE(o)*Huc_orient(3,2)
          P3 = Puc_orient(i,3) + DBLE(m)*Huc_orient(1,3) + DBLE(n)*Huc_orient(2,3) + DBLE(o)*Huc_orient(3,3)
          !
          !Determine if this position is inside of the polyhedron
          isinpolyhedron = .TRUE.
          DO jnode = 1 , MIN( SIZE(vvertex,1) , maxvertex )
            !Compute vector between the vertex and current node
            !By definition this vector is normal to the grain boundary
            vnormal(:) = vvertex(jnode,:) - vnodes(inode,:)
            !Compute vector between atom replica and current node
            vector(:) = (/P1,P2,P3/) - vnodes(inode,:)
            IF( VEC_PLANE(vnormal,VECLENGTH(vnormal),vector) >= -0.1d0 ) THEN
              !Atom is above this plane of cut, hence out of the polyhedron
              isinpolyhedron = .FALSE.
            ENDIF
          ENDDO
          !
          IF( isinpolyhedron ) THEN
            !This atom is inside of the polyhedron => save it in array Q
            qi = qi+1
            IF( qi>SIZE(Q,1) ) THEN
              !number of atom replica exceeds size of allocated array Q
              !if this ever happens it means that the estimation of qi above is wrong
              !(i.e. programmer is stupid) and should be fixed
              !Meanwhile just exit smoothly to avoid a segfault
              !nerr = nerr+1
              !CALL ATOMSK_MSG(4821,(/""/),(/0.d0/))
              !GOTO 1000
              !Increase size of Q for 1000 more atoms
              IF(ALLOCATED(newQ)) DEALLOCATE(newQ)
              ALLOCATE( newQ( SIZE(Q,1)+1000 , 4 ) )
              DO j=1,SIZE(Q,1)
                newQ(j,:) = Q(j,:)
              ENDDO
              DEALLOCATE(Q)
              ALLOCATE(Q(SIZE(newQ,1),4))
              Q(:,:) = newQ(:,:)
              DEALLOCATE(newQ)
            ENDIF
            Q(qi,1) = P1
            Q(qi,2) = P2
            Q(qi,3) = P3
            Q(qi,4) = Puc_orient(i,4)
            !Also duplicate shells if any
            IF( doshells ) THEN
              T(qi,1) = Suc_orient(i,1) + DBLE(m)*Huc_orient(1,1) + DBLE(n)*Huc_orient(2,1) + DBLE(o)*Huc_orient(3,1)
              T(qi,2) = Suc_orient(i,2) + DBLE(m)*Huc_orient(1,2) + DBLE(n)*Huc_orient(2,2) + DBLE(o)*Huc_orient(3,2)
              T(qi,3) = Suc_orient(i,3) + DBLE(m)*Huc_orient(1,3) + DBLE(n)*Huc_orient(2,3) + DBLE(o)*Huc_orient(3,3)
              T(qi,4) = Suc_orient(i,4)
            ENDIF
            !Duplicated particles will have same auxiliary properties as the originals
            IF(doaux) THEN
              AUX_Q(qi,:) = AUXuc(i,:)
            ENDIF
            !
          ENDIF
          !
        ENDDO
      ENDDO
    ENDDO
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
    DO i=1,qi
      newP(NP+i,:) = Q(i,:)
    ENDDO
    IF(ALLOCATED(P)) DEALLOCATE(P)
    ALLOCATE( P(SIZE(newP,1),4) )
    P(:,:) = newP(:,:)
    DEALLOCATE(newP)
    DEALLOCATE(Q)
    !
    !Same with shells if any
    IF(doshells) THEN
      IF(ALLOCATED(newS)) DEALLOCATE(newS)
      ALLOCATE(newS(NP+qi,4))
      newS(:,:) = 0.d0
      IF( ALLOCATED(S) ) THEN
        DO i=1,SIZE(S,1)
          newS(i,:) = S(i,:)
        ENDDO
      ENDIF
      DO i=1,qi
        newS(NP+i,:) = T(i,:)
      ENDDO
      IF(ALLOCATED(S)) DEALLOCATE(S)
      ALLOCATE( S(SIZE(newS,1),4) )
      S(:,:) = newS(:,:)
      DEALLOCATE(newS)
      DEALLOCATE(T)
    ENDIF
    !
    !Also copy auxiliary properties if any
    IF(doaux) THEN
      IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
      ALLOCATE(newAUX(NP+qi,SIZE(AUXuc,2)))
      newAUX(:,:) = 0.d0
      IF( ALLOCATED(AUX) ) THEN
        DO i=1,SIZE(AUX,1)
          newAUX(i,:) = AUX(i,:)
        ENDDO
      ENDIF
      DO i=1,qi
        newAUX(NP+i,:) = AUX_Q(i,:)
      ENDDO
      IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
      ALLOCATE( AUX( SIZE(newAUX,1) , SIZE(AUXuc,2) ) )
      AUX(:,:) = newAUX(:,:)
      DEALLOCATE(newAUX)
      DEALLOCATE(AUX_Q)
    ENDIF
    !
    NP = NP+qi
    !
  ELSE  !i.e. qi=0
    !Display a warning message
    nwarn=nwarn+1
    CALL ATOMSK_MSG(4708,(/""/),(/0.d0/))
  ENDIF
  !
  IF(ALLOCATED(vvertex)) DEALLOCATE(vvertex)
  !
ENDDO
!
!
!
500 CONTINUE
!P now contains positions of all atoms in all the grains
!Apply options to the final system
CALL OPTIONS_AFF(options_array,H,P,S,AUXNAMES,AUX,ORIENT)
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
