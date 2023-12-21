MODULE cmpt_rdf
!
!**********************************************************************************
!*  COMPUTE_RDF                                                                   *
!**********************************************************************************
!* This module computes the radial distribution function (RDF)                    *
!* of a given system. The number of neighbors in a "skin" of radius R and         *
!* of width dR is computed for each atom in the system, then the average          *
!* number of neighbors N is computed. The RDF is then defined as:                 *
!*    RDF(R) = N(R,R+dR) / dR                                                     *
!* When several atom species are present, the partial RDFs are computed.          *
!* E.g. if atoms A and B are present, the distribution of A atoms around          *
!* A atoms will be computed, then B around A, and B around B.                     *
!* If a list of pairs of atoms is provided (in array "pairs"), then RDFs will     *
!* be computed only for these pairs of atoms. Otherwise, a pair list will be      *
!* constructed, and provided as output of the current routine.                    *
!* The results are output in the array "rdf_sys".                                 *
!**********************************************************************************
!* (C) May. 2012 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 21 Dec. 2023                                     *
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
!Load modules
USE atoms
USE comv
USE constants
USE functions
USE messages
USE neighbors
USE files
USE subroutines
!
!
CONTAINS
!
SUBROUTINE COMPUTE_RDF(H,P,rdf_maxR,rdf_dr,NeighList,aentries,pairs,rdf_sys)
!
!Declare variables
IMPLICIT NONE
!Input
REAL(dp),INTENT(IN):: rdf_maxR            !maximum radius for RDF
REAL(dp),INTENT(IN):: rdf_dr              !width of the "skin"
REAL(dp),DIMENSION(3,3),INTENT(IN):: H    !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P  !atom positions
!
CHARACTER(LEN=4096):: msg
INTEGER:: atompair
INTEGER:: i, id, j, k, l, m, n
INTEGER:: Nspecies   !number of different atom species in the system
INTEGER:: Npairs     !number of unique pairs of atoms
INTEGER:: u, umin, umax, v, vmin, vmax, w, wmin, wmax
INTEGER:: progress   !To show calculation progress
INTEGER:: rdf_Nsteps !number of "skins" for RDF
INTEGER,DIMENSION(:,:),ALLOCATABLE:: NeighList  !list of neighbours
REAL(dp):: distance     !distance between 2 atoms
REAL(dp):: rdf_norm     !normalization factor
REAL(dp):: rdf_radius   !radius of the sphere
REAL(dp):: Vsphere, Vskin        !volume of the sphere, skin
REAL(dp):: Vsystem               !volume of the system
REAL(dp),DIMENSION(3):: Vector
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries  !number of atoms of each chem. species
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: pairs     !pairs of species
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: rdf_sys   !result: partial RDFs for input system
!
msg = "ENTERING COMPUTE_RDF..."
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Initialize variables
atompair=0
u=0
v=0
w=0
IF(ALLOCATED(rdf_sys)) DEALLOCATE(rdf_sys)     !partial RDFs for 1 system
!
!
!
100 CONTINUE
!Set number of steps
rdf_Nsteps = NINT(rdf_maxR/rdf_dr) + 2
!
!Determine if we have to look for periodic replica of atoms
!Minimum image convention: look for all replica in the radius rdf_maxR
umin=-1
umax=1
vmin=-1
vmax=1
wmin=-1
wmax=1
IF ( VECLENGTH(H(1,:))<=rdf_maxR ) THEN
  umax = CEILING( rdf_maxR / VECLENGTH(H(1,:)) ) + 1
  umin = -1*umax
ENDIF
IF ( VECLENGTH(H(2,:))<=rdf_maxR ) THEN
  vmax = CEILING( rdf_maxR / VECLENGTH(H(2,:)) ) + 1
  vmin = -1*vmax
ENDIF
IF ( VECLENGTH(H(3,:))<=rdf_maxR ) THEN
  wmax = CEILING( rdf_maxR / VECLENGTH(H(3,:)) ) + 1
  wmin = -1*wmax
ENDIF
!
!Count how many different species exist in the system
IF( .NOT.ALLOCATED(aentries) ) THEN
  CALL FIND_NSP(P(:,4),aentries)
ENDIF
Nspecies = SIZE(aentries,1)
!Determine number of unique pais of atoms
Npairs = Nspecies*(Nspecies+1)/2
!
IF( verbosity>=4 ) THEN
  WRITE(msg,*) "Nspecies = ", Nspecies
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO k=1,SIZE(aentries,1)
    WRITE(msg,*) "    #", NINT(aentries(k,1)), ": ", NINT(aentries(k,2)), " atoms"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
ENDIF
!
IF( .NOT.ALLOCATED(pairs) ) THEN
  !Construct array storing pairs of atoms
  ALLOCATE( pairs(Npairs,2) )
  pairs(:,2) = 0.d0
  m=0
  DO i=1,SIZE(aentries,1)
    DO j=i,SIZE(aentries,1)
      m=m+1
      pairs(m,1) = aentries(i,1)
      pairs(m,2) = aentries(j,1)
    ENDDO
  ENDDO
ENDIF
IF( verbosity>=4 ) THEN
  WRITE(msg,*) "N. pairs = ", SIZE(pairs,1)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO k=1,SIZE(pairs,1)
    WRITE(msg,*) "    #", k, ": ", NINT(pairs(k,1)), NINT(pairs(k,2))
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
ENDIF
!
!If a neighbour list was provided as input, it will be used
!Otherwise, a neighbour list is constructed
IF( .NOT.ALLOCATED(NeighList) .OR. SIZE(NeighList,1)<SIZE(P,1) ) THEN
  IF(ALLOCATED(NeighList)) DEALLOCATE(NeighList)
  !Construct neighbour list
  !NOTE: to compute RDF up to R, we need neighbors up to R+dR
  CALL ATOMSK_MSG(11,(/""/),(/0.d0/))
  CALL NEIGHBOR_LIST(H,P,rdf_maxR+1.1d0*rdf_dr,NeighList)
  CALL ATOMSK_MSG(15,(/""/),(/0.d0/))
ENDIF
!
!
!
200 CONTINUE
!Construct partial RDFs for current system
IF(ALLOCATED(rdf_sys)) DEALLOCATE(rdf_sys)
ALLOCATE( rdf_sys( SIZE(pairs,1) , rdf_Nsteps ) )
rdf_sys(:,:) = 0.d0
!Loop on all atoms
progress=0
CALL ATOMSK_MSG(4052,(/""/),(/0.d0/))
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP& PRIVATE(i,id,j,m,n,u,v,w,atompair,distance,Vector) &
!$OMP& REDUCTION(+:rdf_sys)
DO i=1,SIZE(P,1)
  !
  IF( SIZE(P,1)>20000 ) THEN
    progress = progress+1
    !If there are many atoms, display a fancy progress bar
    CALL ATOMSK_MSG(10,(/""/),(/DBLE(progress),DBLE(SIZE(P,1))/))
  ENDIF
  !
  !Get index of "atompair" for the pair #i-#i
  DO j=1,SIZE(pairs,1)
    IF( DABS(P(i,4)-pairs(j,1)) < 1.d-3 .AND. &
        & DABS(P(i,4)-pairs(j,2)) < 1.d-3     ) THEN
      atompair = j
      EXIT
    ENDIF
  ENDDO
  !
  !Check if replicas of atom #i are neighbors of atom #i
  !NOTE: exclude self image as it would result in distance=0
  DO u=umin,umax
    DO v=vmin,vmax
      DO w=wmin,wmax
        distance = VECLENGTH( DBLE(u)*H(1,:) + DBLE(v)*H(2,:) + DBLE(w)*H(3,:) )
        IF( DABS(distance)>1.d-12 .AND. distance<(rdf_maxR+rdf_dr) ) THEN
          !Add this replica in the appropriate "atompair" and appropriate skin
          !n = index of skin
          n = MAX(1,CEILING(distance/rdf_dr))
          rdf_sys(atompair,n) = rdf_sys(atompair,n) + 1
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  !
  !Parse the neighbour list of atom #i
  m=1
  DO WHILE ( m<=SIZE(NeighList,2)  .AND. NeighList(i,m)>0 )
    !Get the index of the #m-th neighbor of atom #i
    id = NeighList(i,m)
    !Get index of "atompair" for the pair #i-#id
    atompair = 0
    DO j=1,SIZE(pairs,1)
      IF( DABS(P(i,4)-pairs(j,1)) < 1.d-3 .AND. &
        & DABS(P(id,4)-pairs(j,2)) < 1.d-3      ) THEN
        atompair = j
        EXIT
      ENDIF
    ENDDO
    IF( atompair>0 ) THEN
      !Compute distance between atoms #i and #id
      Vector(:) = P(id,1:3) - P(i,1:3)
      DO u=umin,umax
        DO v=vmin,vmax
          DO w=wmin,wmax
            distance = VECLENGTH( Vector(:) &
                    &   + DBLE(u)*H(1,:) + DBLE(v)*H(2,:) + DBLE(w)*H(3,:) )
            IF( distance<(rdf_maxR+rdf_dr) ) THEN
              !Add this neighbour in the appropriate "atompair" and appropriate skin
              !n = index of skin
              n = MAX(1,CEILING(distance/rdf_dr))
              rdf_sys(atompair,n) = rdf_sys(atompair,n) + 1
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDIF !end if atompair>0
    !Go to next neighbour
    m=m+1
  ENDDO
  !
ENDDO  !loop on atoms (i)
!$OMP END PARALLEL DO
!
IF( verbosity>=4 ) THEN
  WRITE(msg,*) "Done, raw neighbour count:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(pairs,1)
    WRITE(msg,'(4X,i3,2X,i3,2X,i9)') NINT(pairs(i,1)), NINT(pairs(i,2)), NINT(SUM(rdf_sys(i,:)))
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
ENDIF
!
!
!
300 CONTINUE
!At this point, rdf_sys(:,:) contains the sum of all neighbours for this system
!It has to be normalized to the atom density
!Compute total volume of current system
CALL VOLUME_PARA(H,Vsystem)
!Normalize each partial RDF for current system
DO i=1,SIZE(pairs,1)
  !Pair of atoms (sp1,sp2): compute average density of second species
  DO j=1,SIZE(aentries,1)
    !Save number of atoms sp1
    IF( DABS(aentries(j,1)-pairs(i,1))<1.d-3 ) THEN
      u = NINT(aentries(j,2))
    ENDIF
    !Save average density of atoms sp2
    IF( DABS(aentries(j,1)-pairs(i,2))<1.d-3 ) THEN
      v = NINT(aentries(j,2))
    ENDIF
  ENDDO
  !
  !Loop on all skins of width rdf_dr
  DO j=1,rdf_Nsteps
    !Radius of current sphere
    rdf_radius = DBLE(j-1) * rdf_dr
    !Volume of sphere of radius R
    Vsphere = (4.d0/3.d0)*pi*(rdf_radius**3)
    !Volume of sphere of radius R+dR
    distance = rdf_radius + rdf_dr
    distance = (4.d0/3.d0)*pi*(distance**3)
    !Volume of the skin = difference between the 2 spheres
    Vskin = distance - Vsphere
    !
    !Compute normalization factor: (Nsp1*Nsp2*Vskin)/V
    rdf_norm = DBLE(u) * DBLE(v) * Vskin /Vsystem
    !
    !Normalize partial RDFs for current system
    rdf_sys(i,j) = rdf_sys(i,j) / rdf_norm
  ENDDO
ENDDO
!
IF( verbosity>=4 ) THEN
  WRITE(msg,*) "Normalized neighbour count:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(pairs,1)
    WRITE(msg,'(4X,i3,2X,i3,2X,i9)') NINT(pairs(i,1)), NINT(pairs(i,2)), NINT(SUM(rdf_sys(i,:)))
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
ENDIF
!
!
!
1000 CONTINUE
msg = "EXITING COMPUTE_RDF..."
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
END SUBROUTINE COMPUTE_RDF
!
!
END MODULE cmpt_rdf
