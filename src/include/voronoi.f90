MODULE voronoi
!
!**********************************************************************************
!*  VORONOI                                                                       *
!**********************************************************************************
!* This module contains routines related to Voronoi tesselation.                  *
!* Basics of Delaunay triangulation are explained in Wikipedia:                   *
!*  https://en.wikipedia.org/wiki/Delaunay_triangulation                          *
!**********************************************************************************
!* (C) June 2025                                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 28 Oct. 2025                                     *
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
!* VORONEIGH           for each node, find neighboring vertices                   *
!* TRIANGULATE_2D      Delaunay triangulation in 2-D                              *
!* TRIANGULATE_3D      Delaunay triangulation in 3-D                              *
!**********************************************************************************
!
!
USE comv
USE constants
USE math
USE subroutines
USE sorting
USE neighbors
USE resize
!
IMPLICIT NONE
!
!
CONTAINS
!
!
!********************************************************
!  VORONEIGH
!  Given a list of nodes, find for each node its direct
!  neighboring vertices.
!********************************************************
SUBROUTINE VORONEIGH(H,twodim,vnodes,Nvertices,vvertex,Volumes)
!
IMPLICIT NONE
LOGICAL:: sameplane  !are all nodes nearly in the same plane?
CHARACTER(LEN=4096):: msg
INTEGER:: i, j, k
INTEGER:: inode, jnode
INTEGER:: thvertex, maxvertex   !theoretical/actual max. number of vertices
INTEGER:: Nnodes, Nneighbors, Nvertex
INTEGER:: status
INTEGER,INTENT(IN):: twodim     !=0 if system is 3-D, =1,2,3 if system is thin along x, y, z
INTEGER,DIMENSION(3):: expandmatrix
INTEGER,DIMENSION(:),ALLOCATABLE:: newindex  !list of index after sorting
INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(OUT):: Nvertices !number of vertices for each node
REAL(dp):: distance, radius
REAL(dp):: P1, P2, P3  !temporary position
REAL(dp):: maxdnodes   !maximum distance between 2 nodes
REAL(dp),DIMENSION(3):: vector, center
REAL(dp),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: Volumes
REAL(dp),DIMENSION(3,2):: nodebox   !Bounding box (xmin,xmax,ymin,ymax,zmin,zmax) surrounding a node
REAL(dp),DIMENSION(3,3),INTENT(IN):: H         !Base vectors of the final supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: tempnodes !cartesian coordinate of neighboring nodes
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: vnodes     !cartesian coordinate of each node
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(OUT):: vvertex !for each node, position of neighboring vertices
!
!
WRITE(msg,*) "Entering VORONEIGH"
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
IF(ALLOCATED(Nvertices)) DEALLOCATE(Nvertices)
IF(ALLOCATED(vvertex)) DEALLOCATE(vvertex)
IF(ALLOCATED(tempnodes)) DEALLOCATE(tempnodes)
IF(ALLOCATED(Volumes)) DEALLOCATE(Volumes)
IF( .NOT.ALLOCATED(vnodes).OR.SIZE(vnodes,1)<=0 ) THEN
  RETURN
ENDIF
!
Nnodes = SIZE(vnodes,1)
!
!The maximum number of faces of any polyhedron should be 11 in 2-D,
!and 59 in 3-D, for proof see e.g.:
!   http://www.ericharshbarger.org/voronoi.html
!To that we will add 26 self-neighbors in 3-D
!This will be used to limit the number of iterations in the loops on jnode below
IF( twodim > 0 ) THEN
  thvertex = 40
ELSE
  thvertex = 100
ENDIF
WRITE(msg,*) "Theoretical max. vertices = ", thvertex
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Neighbors will be searched in adjacent cells (+/-1)
expandmatrix(:) = 1
IF( twodim>0 ) THEN
  !System is pseudo 2-D, do not look for replicas along the short distance
  expandmatrix(twodim) = 0
ENDIF
!
!Allocate memory for vertices
!NOTE: at this stage the number of neighbors is expected to be vastly over-estimated.
!      First entries will contain actual neighbors with vvertex(:,:,4) = distance to node (>0).
!      Following entries will have vvertex(:,:,4) = -1
!      This array will be resized later
ALLOCATE( vvertex(Nnodes,thvertex,4) )
vvertex(:,:,:) = 0.d0
vvertex(:,:,4) = -1.d0
ALLOCATE(Nvertices(Nnodes))
Nvertices(:) = 0
ALLOCATE(tempnodes(27*Nnodes,4))
tempnodes(:,:) = 0.d0
ALLOCATE(Volumes(Nnodes))
Volumes(:) = 0.d0
!
maxvertex = 1
!
!For each node, search vertices that are direct neighbors
!$OMP  PARALLEL DO DEFAULT(SHARED) &
!$OMP& PRIVATE(inode,jnode,nodebox,i,j,k,Nvertex,P1,P2,P3,tempnodes,vector,distance,status,msg)
DO inode=1,Nnodes
  !This is node #inode, its position is vnodes(inode,1:3)
  WRITE(msg,'(a,i6,a8)') "   * * *   NODE  #", inode, "   * * *"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  tempnodes(:,:) = 0.d0
  tempnodes(:,4) = -1.d0
  Nvertex = 0
  !
  !Define bounding box centered on node #inode for neighbor search
  P1 = 1.d0
!   IF( twodim==0 ) THEN
!     IF( Nnodes>40 ) THEN
!       P1 = 0.5d0
!     ELSEIF( Nnodes>20 ) THEN
!       P1 = 0.8d0
!     ENDIF
!   ENDIF
  nodebox(1,1) = vnodes(inode,1) - P1*H(1,1) - 0.1d0   !xmin
  nodebox(1,2) = vnodes(inode,1) + P1*H(1,1) + 0.1d0   !xmax
  nodebox(2,1) = vnodes(inode,2) - P1*H(2,2) - 0.1d0   !ymin
  nodebox(2,2) = vnodes(inode,2) + P1*H(2,2) + 0.1d0   !ymax
  nodebox(3,1) = vnodes(inode,3) - P1*H(3,3) - 0.1d0   !zmin
  nodebox(3,2) = vnodes(inode,3) + P1*H(3,3) + 0.1d0   !zmax
  IF( twodim>0 ) THEN
    !System is pseudo-2D, reduce search box along small dimension
    nodebox(twodim,1) = vnodes(inode,twodim) - 0.1d0
    nodebox(twodim,2) = vnodes(inode,twodim) + 0.1d0
  ENDIF
  !
  !First, search for all nodes that are within a box centered on node #inode
  !Loop on all other nodes
  WRITE(msg,'(a,i6)') "Searching for nodes inside bounding box..."
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO jnode=1,Nnodes
    !Parse all periodic replica of node #jnode
    DO i=-expandmatrix(3),expandmatrix(3)
      DO j=-expandmatrix(2),expandmatrix(2)
        DO k=-expandmatrix(1),expandmatrix(1)
          IF( jnode.NE.inode .OR. i.NE.0 .OR. j.NE.0 .OR. k.NE.0 ) THEN
            !Position of the periodic image of node #jnode
            vector(1) = vnodes(jnode,1) + DBLE(k)*H(1,1)
            vector(2) = vnodes(jnode,2) + DBLE(j)*H(2,2)
            vector(3) = vnodes(jnode,3) + DBLE(i)*H(3,3)
            !Node position = middle point between nodes #inode and this replica of #jnode
            !vector(:) = vnodes(inode,:) + ( (/P1,P2,P3/)-vnodes(inode,:) )/2.d0
            !Check if this vertex is inside bounding box
            IF( vector(1)>=nodebox(1,1) .AND. vector(1)<=nodebox(1,2) .AND. &
              & vector(2)>=nodebox(2,1) .AND. vector(2)<=nodebox(2,2) .AND. &
              & vector(3)>=nodebox(3,1) .AND. vector(3)<=nodebox(3,2)       ) THEN
              !This node is inside bounding box: save it
              IF( Nvertex<SIZE(tempnodes,1) ) THEN
                Nvertex = Nvertex+1
                !Compute distance between node #inode and this vertex
                distance = VECLENGTH( vector(:) - vnodes(inode,1:3) )
                !Save vertex position
                tempnodes(Nvertex,1:3) = vector(:)
                tempnodes(Nvertex,4) = distance
              ENDIF
              !
            ENDIF
          ENDIF
        ENDDO  ! end loop on k
      ENDDO  ! j
    ENDDO  ! i
    !
  ENDDO
  !
  WRITE(msg,'(a,i6)') "Found neighboring nodes in bounding box:", Nvertex
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Perform Delaunay triangulation to find which nodes are actual first neighbors
  !and save positions of vertices into vvertex(inode,:,:)
  IF( twodim>0 ) THEN
    !System is pseudo-2D
    WRITE(msg,'(a,i6,a)') "Performing Delaunay triangulation in 2-D (", Nvertex, " neighbors)..."
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    CALL TRIANGULATE_2D(vnodes(inode,1:3),tempnodes(:,1:3),vvertex(inode,:,:),Nvertex,Volumes(inode),status)
    !
  ELSE  !i.e. if twodim==0
    !System is 3-D
    WRITE(msg,'(a,i6,a)') "Performing Delaunay triangulation in 3-D (", Nvertex, " neighbors)..."
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    CALL TRIANGULATE_3D(vnodes(inode,1:3),tempnodes(:,1:3),vvertex(inode,:,:),Nvertex,Volumes(inode),status)
  ENDIF !end if twodim
  !
  !Save final number of vertices
  Nvertices(inode) = Nvertex
  !
  IF(verbosity==4) THEN
    !Debug messages
    WRITE(msg,*) "Done, neighboring vertices remaining:", Nvertices(inode)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,*) "Final list of neighboring vertices after cleanup:"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    msg = "   N       x           y          z        dist.to current node"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,*) inode
    OPEN(UNIT=40,FILE="atomsk_grain"//TRIM(ADJUSTL(msg))//"_vertices.xyz")
    WRITE(40,*) Nvertices(inode)+1
    WRITE(40,*) "# Position of node # "//TRIM(ADJUSTL(msg))//" and its vertices"
    WRITE(40,'(i4,6f12.3)') 2, vnodes(inode,:)
    DO i=1,Nvertices(inode)
      IF( vvertex(inode,i,4)>1.d-6 ) THEN
        WRITE(msg,'(i4,4f12.3)') i, vvertex(inode,i,1:4)
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        WRITE(40,'(i4,3f12.3)') 1, vvertex(inode,i,1:3)
      ENDIF
    ENDDO
    WRITE(msg,*) "Final grain volume:", Volumes(inode)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    CLOSE(40)
  ENDIF
  !
ENDDO !end loop on inode
!$OMP END PARALLEL DO
!
!
!Resize array vvertex to optimize memory usage
maxvertex = MAXVAL(Nvertices(:))
IF( maxvertex < SIZE(vvertex,2) ) THEN
  WRITE(msg,*) "Resizing array 'maxvertex', max. number of vertices: ", maxvertex
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  IF(twodim>0) THEN
    CALL RESIZE_DBLEARRAY3(vvertex,Nnodes,maxvertex,4,status)
  ELSE
    CALL RESIZE_DBLEARRAY3(vvertex,Nnodes,maxvertex,4,status)
  ENDIF
ENDIF
!
!
1000 CONTINUE
IF(ALLOCATED(tempnodes)) DEALLOCATE(tempnodes)
WRITE(msg,*) "Exiting VORONEIGH"
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
END SUBROUTINE VORONEIGH
!
!
!
!********************************************************
!  TRIANGULATE_2D
!  Given a node position and several neighbors, perform
!  Delaunay triangulation to keep only first neighbors.
!  Array node(:) must contain position of central node.
!  Array neighnodes must contain position (x,y,z) of
!  neighboring nodes, as well as their distance to central node.
!  For all sets of 3 points (node + 2 neighbors), their
!  circumcircle is searched. If it contains no other neighbor,
!  then these 2 neighbors are direct neighbors of the node.
!  The subroutine returns the position of neighboring vertices.
!  NOTE: for practical purposes, even for this 2-D version
!  the arrays must be 3-D.
!********************************************************
SUBROUTINE TRIANGULATE_2D(node,neighnodes,neighvert,Nvertex,Area,status)
!
IMPLICIT NONE
REAL(dp),DIMENSION(3),INTENT(IN):: node         !Position of central node
REAL(dp),DIMENSION(:,:),INTENT(IN):: neighnodes !Position of neighboring nodes
REAL(dp),INTENT(OUT):: Area
REAL(dp),DIMENSION(:,:),INTENT(OUT):: neighvert !Final position of vertices
INTEGER,INTENT(INOUT):: Nvertex !Initial number of neighboring nodes/final number of vertices
!
INTEGER:: i, j, m, n
INTEGER:: status, stattmp
REAL(dp):: P1, P2, P3, distance, radius
REAL(dp),DIMENSION(3):: center, normal
REAL(dp),DIMENSION(3,3):: P   !Position of points at surface of sphere
LOGICAL,DIMENSION(SIZE(neighnodes,1)):: keepvertex
!
status = 0
Area = 0.d0
keepvertex(:) = .FALSE.
!
IF( Nvertex>3 ) THEN
  !$OMP  PARALLEL DO DEFAULT(SHARED) &
  !$OMP& PRIVATE(i,j,m,n,P,P1,P2,P3,distance,center,radius,normal,stattmp) REDUCTION(+:Area)
  DO i=1,SIZE(neighnodes,1)-1
    IF( neighnodes(i,4)>1.d-6 ) THEN
      DO j=i+1,SIZE(neighnodes,1)
        IF( neighnodes(j,4)<1.d-6 ) EXIT
        !For this set of points (node + neighbors #i and #j), find their circumcircle
        P(1,:) = node(1:3)
        P(2,:) = neighnodes(i,1:3)
        P(3,:) = neighnodes(j,1:3)
        CALL CIRCUMCIRCLE(P(1:3,1:3),center,radius,normal,stattmp)
        IF( stattmp>0 ) THEN
          !These neighbors are probably aligned due to PBC
          P1 = DOT_PRODUCT( P(2,:)-P(1,:) , P(3,:)-P(1,:) )
          P2 = VECLENGTH( P(2,:)-P(1,:) )
          P3 = VECLENGTH( P(3,:)-P(1,:) )
          IF( P1 > 0.9d0*P2*P3 ) THEN
            !Relative vectors point (almost) in the same direction
            IF( P2<P3 ) THEN
              keepvertex(i) = .TRUE.
            ELSE
              keepvertex(j) = .TRUE.
            ENDIF
          ELSE
            !They point in opposite directions: keep both
            keepvertex(i) = .TRUE.
            keepvertex(j) = .TRUE.
          ENDIF
        ELSE
          !Circumcircle has center and radius
          !Check if there are other nodes inside this circumcircle
          n=0
          DO m=1,SIZE(neighnodes,1)
            IF( neighnodes(m,4)>1.d-6 ) THEN
              IF( m.NE.i .AND. m.NE.j ) THEN
                !Compute distance between node #m and center
                distance = VECLENGTH( neighnodes(m,1:3)-center(:) )
                IF( distance < radius-1.d0 ) THEN
                  !Node #m is inside circumcircle: increment n and exit
                  n=n+1
                  EXIT
                ENDIF
              ENDIF
            ENDIF
          ENDDO
          IF( n==0 ) THEN
            !No other node inside the circumcircle: keep all nodes
            keepvertex(i) = .TRUE.
            keepvertex(j) = .TRUE.
            !Area of triangle = half the area of parallelogram
            Area = Area + 0.5d0*VECLENGTH( CROSS_PRODUCT( P(3,:)-P(1,:) , P(2,:)-P(1,:) ) )
          ENDIF
        ENDIF  !end if stattmp>0
      ENDDO !j
    ENDIF
  ENDDO !i
  !$OMP END PARALLEL DO
  !
  Nvertex = 0
  DO i=1,SIZE(neighnodes,1)
    IF( neighnodes(i,4)<1.d-6 ) EXIT
    IF( keepvertex(i) ) THEN
      Nvertex = Nvertex+1
      neighvert(Nvertex,1:3) = node(:) + ( neighnodes(i,1:3)-node(:) )/2.d0
      neighvert(Nvertex,4) = VECLENGTH( neighvert(Nvertex,1:3) - node(:) )
    ENDIF
  ENDDO
  !
ELSE !i.e. number of neighbors <=3
  IF( Nvertex>0 ) THEN
    !Just compute vertices for all neighbors and save into neighvert
    DO i=1,Nvertex
      neighvert(i,1:3) = node(:) + ( neighnodes(i,1:3)-node(:) )/2.d0
      neighvert(i,4) = VECLENGTH( neighvert(i,1:3) - node(:) )
    ENDDO
  ENDIF
ENDIF
!
END SUBROUTINE TRIANGULATE_2D
!
!
!
!********************************************************
!  TRIANGULATE_3D
!  Given a node position and several neighbors, perform
!  Delaunay triangulation to keep only first neighbors.
!  Array node(:) must contain position of central node.
!  Array neighnodes must contain position (x,y,z) of
!  neighboring nodes, as well as their distance to central node.
!  For all sets of 4 points (node + 3 neighbors), their
!  circumsphere is searched. If it contains no other neighbor,
!  then these 3 neighbors are direct neighbors of the node.
!  The subroutine returns the position of neighboring vertices.
!  Basics of Delaunay triangulation are explained in Wikipedia:
!  https://en.wikipedia.org/wiki/Delaunay_triangulation
!********************************************************
SUBROUTINE TRIANGULATE_3D(node,neighnodes,neighvert,Nvertex,Volume,status)
REAL(dp),DIMENSION(3),INTENT(IN):: node         !Position of central node
REAL(dp),DIMENSION(:,:),INTENT(IN):: neighnodes !Position of neighboring nodes
REAL(dp),INTENT(OUT):: Volume                   !Volume of
REAL(dp),DIMENSION(:,:),INTENT(OUT):: neighvert !Final position of vertices
INTEGER,INTENT(INOUT):: Nvertex !Initial number of neighboring nodes/final number of vertices
!
INTEGER:: i, j, k, m, n
INTEGER:: status, stattmp
REAL(dp):: distance, radius
REAL(dp),DIMENSION(3):: center
REAL(dp),DIMENSION(4,3):: P   !Position of points at surface of sphere
LOGICAL,DIMENSION(SIZE(neighnodes,1)):: keepvertex
!
status = 0
Volume = 0.d0
keepvertex(:) = .FALSE.
!
IF( Nvertex>4 ) THEN
  !$OMP  PARALLEL DO DEFAULT(SHARED) &
  !$OMP& PRIVATE(i,j,k,m,n,P,distance,center,radius,stattmp) REDUCTION(+:Volume)
  DO i=1,SIZE(neighnodes,1)-2
    IF( neighnodes(i,4)>1.d-6 ) THEN
      DO j=i+1,SIZE(neighnodes,1)-1
        IF( neighnodes(j,4)<1.d-6 ) EXIT
        DO k=j+1,SIZE(neighnodes,1)
          IF( neighnodes(k,4)<1.d-6 ) EXIT
          !For this set of points (node #inode + neighbors #i, #j and #k), find their circumsphere
          P(1,:) = node(1:3)
          P(2,:) = neighnodes(i,1:3)
          P(3,:) = neighnodes(j,1:3)
          P(4,:) = neighnodes(k,1:3)
          CALL CIRCUMSPHERE(P,center,radius,stattmp)
          IF( stattmp>0 ) THEN
            !The points may belong to the same plane: compute circumcircle
            !CALL CIRCUMCIRCLE(P(1:3,1:3),center,radius,normal,stattmp)
            !keepvertex(i) = .TRUE.
            !keepvertex(j) = .TRUE.
            !keepvertex(k) = .TRUE.
          ELSE
            !Circumsphere has center and radius
            !Check if there are other nodes inside this circumsphere
            n=0
            DO m=1,SIZE(neighnodes,1)
              IF( neighnodes(m,4)<1.d-6 ) EXIT
              IF( m.NE.i .AND. m.NE.j .AND. m.NE.k ) THEN
                !Compute distance between node #m and center
                distance = VECLENGTH( neighnodes(m,1:3)-center )
                IF( distance < radius-1.d0 ) THEN
                  !Node #m is inside circumsphere: increment n
                  n=n+1
                  EXIT
                ENDIF
              ENDIF
            ENDDO
            IF( n==0 ) THEN
              !No other node inside the circumsphere: keep all nodes
              keepvertex(i) = .TRUE.
              keepvertex(j) = .TRUE.
              keepvertex(k) = .TRUE.
              Volume = Volume + VOLUME_TETRA(P(1,1:3),P(2,1:3),P(3,1:3),P(4,1:3))
            ENDIF
          ENDIF
        ENDDO  !k
      ENDDO  !j
    ENDIF
  ENDDO  !i
  !$OMP END PARALLEL DO
  !
  Nvertex = 0
  DO i=1,SIZE(neighnodes,1)
    IF( neighnodes(i,4)<1.d-6 ) EXIT
    IF( keepvertex(i) ) THEN
      Nvertex = Nvertex+1
      neighvert(Nvertex,1:3) = node(:) + ( neighnodes(i,1:3)-node(:) )/2.d0
      neighvert(Nvertex,4) = VECLENGTH( neighvert(Nvertex,1:3) - node(:) )
    ENDIF
  ENDDO
  !
ELSE !i.e. number of neighbors <=3
  IF( Nvertex>0 ) THEN
    !Just compute vertices for all neighbors and save into neighvert
    DO i=1,Nvertex
      neighvert(i,1:3) = node(:) + ( neighnodes(i,1:3)-node(:) )/2.d0
      neighvert(i,4) = VECLENGTH( neighvert(i,1:3) - node(:) )
    ENDDO
  ENDIF
ENDIF
!
END SUBROUTINE TRIANGULATE_3D
!
!
!
END MODULE voronoi
