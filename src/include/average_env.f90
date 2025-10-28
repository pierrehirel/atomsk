MODULE avgenv
!
!**********************************************************************************
!*  AVGENV                                                                        *
!**********************************************************************************
!* This module parses all atoms in the system, detects the different types        *
!* of sites (based on central atom, number of neighbours...)                      *
!* and construct "reference" sites (Pref) by averaging all atoms that have        *
!* the same or similar environments.                                              *
!* If different lattices exist in the system (e.g. fcc and bcc), this module      *
!* should be able to recognize each as a different site/environment.              *
!* Atoms that have an odd number of neighbours (due to free surfaces,             *
!* vacancies, etc.) are not included in the averaging.                            *
!* Averaged environments are returned in array Pref(:,:,:)                        *
!* - first index is the site index [1-Nsites]                                     *
!* - second index is the atom index [1-Natoms];                                   *
!*   #1 = central atom, next ones are neighbours                                  *
!* - third index [1-6] is for atom position ([1-3]=x,y,z), atomic number [4],     *
!*   number of atoms in this type of site [5], and number of neighbours [6]       *
!* So, data for the first type of site is in Pref(1,:,:):                         *
!* - central atom position is Pref(1,1,1:3), atomic number Pref(1,1,4)            *
!* - number of atoms matching this type of site is in Pref(1,1,5)                 *
!* - number of neighbours for this type of site is Nneigh=Pref(1,1,6)             *
!* - rel.positions of neighbours in Pref(1,n,1:4), with n=2 up to Nneigh          *
!* Similarly, data for the second type of site is in Pref(2,:,:), and so on.      *
!**********************************************************************************
!* (C) February 2023 - Pierre Hirel                                               *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 27 Aug. 2025                                     *
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
USE files
USE neighbors
USE sorting
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE AVG_ENV(H,P,NeighList,Pref,NeighFactor,siteindex)
!
IMPLICIT NONE
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: msg, temp
INTEGER:: NeighIndex  !index of a neighbour
INTEGER:: iat, i, j, k, m, n
INTEGER:: Nneighbors
INTEGER:: NNmin=3     !minimum number of neighbours to keep
INTEGER:: Nsites, NneighMax
INTEGER:: progress      !To show calculation progress
INTEGER,DIMENSION(:),ALLOCATABLE:: newindex  !list of index after sorting
INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(OUT):: siteindex !for each atom in P, index of its type of site in Pref
INTEGER,DIMENSION(:,:),ALLOCATABLE:: NeighList !neighbour list can be provided as input, or computed below
REAL(dp):: alpha, theta  !angle between two vectors
REAL(dp),INTENT(IN):: NeighFactor !%of tolerance in the radius for neighbor search
REAL(dp):: radius=8.d0  !R for neighbor search: 8 A should be enough to find some neighbors in any system
REAL(dp):: tempreal
REAL(dp),DIMENSION(9),PARAMETER:: known_angle= &
       & (/30.d0,45.d0,60.d0,90.d0,109.28d0,109.47122d0,120.d0,150.d0,180.d0/)
REAL(dp),DIMENSION(3):: vector                !just a vector
REAL(dp),DIMENSION(3,3),INTENT(IN):: H        !cell vectors
REAL(dp),DIMENSION(:,:),INTENT(IN):: P        !atom positions
REAL(dp),DIMENSION(20,21,6):: Tref            !references for atoms environments (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: PosList !list of positions of neighbors of one atom
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(OUT):: Pref  !references for atoms environments (final)
!
WRITE(msg,*) "  ENTERING  AVG_ENV , NeighFactor = ", NeighFactor
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
Nsites=0
Nneighmax=0
Tref(:,:,:) = 0.d0
IF(ALLOCATED(Pref)) DEALLOCATE(Pref)
!
IF(ALLOCATED(siteindex)) DEALLOCATE(siteindex)
ALLOCATE(siteindex(SIZE(P,1)))
siteindex(:) = 0
!
!If provided array NeighList already contains the neighbour list, use it
!Otherwise, construct the neighbour list now
IF( .NOT.ALLOCATED(NeighList) .OR. SIZE(NeighList,1).NE.SIZE(P,1) ) THEN
  msg = "  NeighList empty: constructing neighbour list ..."
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  CALL NEIGHBOR_LIST(H,P,radius,NeighList)
ENDIF
!
!
100 CONTINUE
!Loop on atoms in the system
msg = "  Parsing atoms ..."
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!!!!!$OMP !PARALLEL DO DEFAULT(SHARED) &
!!!!!$OMP& !PRIVATE(i,iat,j,k,n,alpha,tempreal,PosList,newindex,Nneighbors) &
!!!!!$OMP& !REDUCTION(+:Tref(:,:,5))
DO iat=1,SIZE(P,1)
  !
  progress = progress+1
  !
  IF( SIZE(P,1) > 100000 ) THEN
    !If there are many atoms, display a fancy progress bar
    CALL ATOMSK_MSG(10,(/""/),(/DBLE(progress),DBLE(SIZE(P,1))/))
  ENDIF
  !
  !Using neighbor list, save Cartesian neighbours positions in PosList
  CALL NEIGHBOR_POS(H,P,P(iat,1:3),NeighList(iat,:),ALLOCATED(NeighList),radius,PosList)
  !
  !Sort neighbours by increasing distance
  CALL BUBBLESORT(PosList,4,'up  ',newindex)
  !
  !Get the relative positions of neighbors with respect to atom #iat
  !Keep only neighbors that are within NeighFactor times the distance of third neighbor
  Nneighbors = 0
  DO j=1,MIN(20,SIZE(PosList,1))
    !Shift all neighbours positions so the central atom is at (0,0,0)
    PosList(j,1:3) = PosList(j,1:3) - P(iat,1:3)
    IF( PosList(j,4) <= NeighFactor*PosList(NNmin,4) ) THEN
      !IF( NINT(Psecond(NINT(PosList(j,5)),4))==NINT(Psecond(NINT(PosList(1,5)),4)) ) THEN
        !Keep this neighbor
        Nneighbors = Nneighbors+1
      !ELSE
        !Ignore this neighbor (not same species as first neighbor)
      !ENDIF
    ELSE
      !Wipe out the rest of the list
      PosList(j+1:,:) = 0.d0
      EXIT
    ENDIF
  ENDDO
  !
  !Search in Tref if central atom #iat already matches an environment in Tref
  !If not, create a new type of site in Tref
  i=0
  DO WHILE(i<SIZE(Tref,1))
    i=i+1
    IF(i>SIZE(Tref,1)) THEN
      nerr=nerr+1
      CALL ATOMSK_MSG(4070,(/""/),(/DBLE(i)/))
      EXIT
    ENDIF
    IF( NINT(P(iat,4))==NINT(Tref(i,1,4)) .AND. ABS(NINT(Tref(i,1,6))-Nneighbors)<=1 ) THEN
      !Atom species and number of neighbors match
      !=> atom #iat occupies a site of the type #i
      !Save the index of the site for atom #iat
      siteindex(iat) = i
      !Increment counter of atoms for this type of atom site
      Tref(i,1,5) = Tref(i,1,5) + 1.d0
      EXIT
    ELSE IF( NINT(Tref(i,1,4))==0 .AND. MOD(Nneighbors,2)==0 .AND. MOD(Nneighbors,5)>0 ) THEN
      !No suitable site was found before: if number of neighbors is even,
      !then create a new one in an empty slot (i.e. where Tref(i,1,4)=0)
      !Save species of central atom for this site
      Tref(i,1,4) = P(iat,4)
      !Increment counter of atoms for this type of atom site
      Tref(i,1,5) = Tref(i,1,5) + 1.d0
      !Save number of neighbors for this site
      Tref(i,1,6) = DBLE(Nneighbors)
      !Save the index of the site for atom #iat
      siteindex(iat) = i
      EXIT
    ENDIF
  ENDDO
  !
  !Add atom to average only if number of neighbors is greater than 3, lower than 14, even, and
  !not multiple of 5 (because an atom in any material is not supposed to have 5 or 10 neighbors)
  !i.e. we consider only atoms with 4, 6, 8, 12, 14 neighbors are in "perfect" environments
  IF( Nneighbors>3 .AND. Nneighbors<14 .AND. MOD(Nneighbors,2)==0 .AND. MOD(Nneighbors,5)>0 ) THEN
    !Add neighbors positions to already known positions for averaging
    IF( NINT(Tref(i,1,5))==1 ) THEN
      !No neighbor was found for this site yet
      !Simply store rel. pos. of neighbors of atom #iat in Tref(i,:,:)
      DO k=1,Nneighbors
        NeighIndex = NINT(PosList(k,5))
        Tref(i,k+1,1:3) = PosList(k,1:3)
        Tref(i,k+1,4) = P(NeighIndex,4)
      ENDDO
    ELSE
      !Some atoms were already found in this site
      !For each neighbor in PosList, find its best match in Tref 
      !(i.e. the one that maximizes the dot product) and add its position to it
      DO k=1,Nneighbors
        n=0
        alpha = 0.d0
        j=1
        DO WHILE( j<=SIZE(Tref,2) )
          tempreal = DOT_PRODUCT( PosList(k,1:3) , Tref(i,j,1:3) )
          IF( tempreal > alpha ) THEN
            !This neighbor position is a better match
            n = j
            alpha = tempreal
          ENDIF
          j=j+1
        ENDDO
        IF (n>0) THEN
          !Now n is the index of the best matching neighbor in Tref
          !Add position of atom #n and perform averaging
          Tref(i,n,1:3) = ( Tref(i,n,5)*Tref(i,n,1:3) + PosList(k,1:3) ) / (Tref(i,n,5)+1.d0)
          !Save atomic number of this neighbor (only if it's empty)
          NeighIndex = NINT(PosList(k,5))
          IF( NINT(Tref(i,n,4))==0 ) THEN
            Tref(i,n,4) = P(NeighIndex,4)
          ENDIF
          !Increment number of neighbors at this position
          Tref(i,n,5) = Tref(i,n,5) + 1.d0
        ELSE
          !n is zero, meaning that the dot product was zero for all neighbors
          !This should not happen, but here we are
          !PRINT*, "         ERROR no matching site for neighbor #", k
        ENDIF
        !
      ENDDO
    ENDIF  !end if (Nneighbors...)
    !
  ENDIF
  !
ENDDO !loop on iat
!!!!!$OMP !END PARALLEL DO
!
!
!Get number of different sites and max. number of neighbours
Nsites = MAXVAL(siteindex)
Nneighmax = MAXVAL(Tref(:,:,6))
WRITE(msg,*) "  Sites detected: ", Nsites
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
WRITE(msg,*) "  Max. N.neigh.:  ", Nneighmax
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!Save averaged sites in final array Pref(:,:,:)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
IF( Nsites>0 ) THEN
  ALLOCATE( Pref(Nsites,Nneighmax+1,6) )
  Pref(:,:,:) = 0.d0
  DO i=1,Nsites
    DO j=1,Nneighmax+1
      Pref(i,j,:) = Tref(i,j,:)
    ENDDO
  ENDDO
ELSE
  !No atomic environment found: make sure Pref and siteindex are returned unallocated
  IF(ALLOCATED(Pref)) DEALLOCATE(Pref)
  IF(ALLOCATED(siteindex)) DEALLOCATE(siteindex)
ENDIF
!
!
200 CONTINUE
IF( ALLOCATED(Pref) ) THEN
  !
  !Now Pref contains the averaged relative positions of neighbors for each type of atom site
  !For each central atom, compute angles between neighbors, and if they are close to
  !some known/expected angles, correct them
  !
  DO i=1,SIZE(Pref,1)
    !Loop on all neighbors
    DO j=2,NINT(Pref(i,1,6))-1
      !Loop on all other neighbors
      DO k=j,NINT(Pref(i,1,6))
        n=0
        !Compute angle between neighbors #j and #k
        alpha = ANGVEC( Pref(i,j,1:3)-Pref(i,1,1:3) , Pref(i,k,1:3)-Pref(i,1,1:3) )
        !Check if angle alpha is close to a known value
        DO m=1,SIZE(known_angle)
          IF( DABS(RAD2DEG(alpha)-known_angle(m)) < 0.1d0 ) THEN
            n=m
            EXIT
          ENDIF
        ENDDO
        IF( n>0 ) THEN
          !Correct position of neighbor #k
          vector(:) = Pref(i,j,1:3)-Pref(i,1,1:3)
          tempreal = VECLENGTH( Pref(i,k,1:3)-Pref(i,1,1:3) )

          Pref(i,k,1:3) = Pref(i,1,1:3) + vector(1:3)
        ENDIF
      ENDDO !k
    ENDDO !j
  ENDDO !i
  !
ENDIF !end if allocated(Pref)
!
!
1000 CONTINUE
!Write some debugging information
! IF( verbosity==4 ) THEN
!   IF( ALLOCATED(Pref) .AND. SIZE(Pref,1)>0 ) THEN
!     !For each site, write position of central atom (that should always be (0,0,0))
!     !and positions of neighbors into a XYZ file named "atomsk_site_i.xyz" for visualization
!     msg = " ATOMIC  ENVIRONMENTS  DETECTED:"
!     CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!     DO i=1,MIN(Nsites,SIZE(Pref,1))
!       !PRINT*, "Site #", i, ": ", NINT(Pref(i,1,5)), " atoms in this site"
!       !Check if number of neighbours for site #i is positive
!       IF( NINT(Pref(i,1,6))>0 ) THEN
!         WRITE(msg,*) i
!         OPEN(UNIT=23,FILE="atomsk_site_"//TRIM(ADJUSTL(msg))//".xyz",FORM="FORMATTED")
!         WRITE(23,*) NINT(Pref(i,1,6))+1
!         WRITE(23,*) "# Averaged environment for site #", i
!         CALL ATOMSPECIES(Pref(i,1,4),species)
!         WRITE(msg,*) i
!         WRITE(temp,*) NINT(Pref(i,1,6))
!         WRITE(23,'(a2,3X,3f16.8)') species, Pref(i,1,1:3)
!         WRITE(msg,*) "  Site # "//TRIM(ADJUSTL(msg))//" occupied by "//species// &
!                     & " atom, has "//TRIM(ADJUSTL(temp))//" neighbors:"
!         CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!         DO j=2,NINT(Pref(i,1,6))+1
!           CALL ATOMSPECIES(Pref(i,j,4),species)
!           WRITE(23,'(a2,3X,3f16.8)') species, Pref(i,j,1:3)
!           WRITE(msg,'(8X,a2,3X,3f12.3)') species, Pref(i,j,1:3)
!           CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!         ENDDO
!         CLOSE(23)
!       ENDIF
!     ENDDO
!     OPEN(UNIT=23,FILE="atomsk_siteindex.txt",FORM="FORMATTED")
!     WRITE(23,*) SIZE(P,1)
!     WRITE(23,*) "# Atomsk AVG_ENV: atom id and type of site they occupy"
!     DO i=1,SIZE(siteindex,1)
!       WRITE(23,*) i, siteindex(i)
!     ENDDO
!     CLOSE(23)
!   ELSE
!     msg = " ERROR:  NO ATOM ENVIRONMENT DETECTED!"
!     CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!   ENDIF
! ENDIF
!
msg = "  EXITING  AVG_ENV"
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!
END SUBROUTINE AVG_ENV
!
END MODULE avgenv
