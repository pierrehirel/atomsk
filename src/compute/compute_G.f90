MODULE cmpt_G
!
!**********************************************************************************
!*  COMPUTE_G                                                                     *
!**********************************************************************************
!* This module reads a set of atomic coordinates, and computes the                *
!* lattice correspondance tensor G for each atom.                                 *
!* The algorithm follows the method described in:                                 *
!*    C.S. Hartley, Y. Mishin, Acta Mater. 53 (2005) 1313                         *
!* Equation numbers also refer to this reference.                                 *
!* The names "system1" and "system2" refer to the reference system                *
!* (bulk with no defect) and the system to analyze, respectively.                 *
!* The reference system 1 can be either:                                          *
!* - a full reference system containing the same number of atoms as system2,      *
!*   and indexed in the exact same order;                                         *
!* - a unit cell of the material in system 2;                                     *
!* - NULL, in which case the "reference" environments are built on-the-fly        *
!*   by averaging atomic environments found in the system2.                       *
!* Neighbour lists (NeighList1 and NeighList2) can be provided as input,          *
!* however if they are not they will be constructed by this routine and           *
!* returned as output.                                                            *
!* In the end this routine returns (N=N.atoms,M=N.sites,NN=N.neighbours):         *
!* - neighbour lists (NeighList1 and 2), if they were not provided as input       *
!* - parameters used for the calculation, in array "params"                       *
!* - the lattice correspondance G for each atom, an array of size Nx3x3           *
!* - reference environments, stored in array Pref of size MxNNx6                  *
!* - the site each atom belongs to, store in siteindex(:) of size N               *
!**********************************************************************************
!* (C) July 2023 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 18 July 2023                                     *
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
USE comv
USE subroutines
USE functions
USE messages
USE neighbors
USE avgenv
!
!
!
CONTAINS
!
!
SUBROUTINE COMPUTE_G(H1,P1,H2,P2,params,NeighList1,NeighList2,Pref,siteindex,G)
!
IMPLICIT NONE
!Input variables and arrays
REAL(dp),DIMENSION(3,3),INTENT(IN):: H1, H2 !Box vectors of systems 1 and 2
REAL(dp),DIMENSION(:,:),ALLOCATABLE,TARGET,INTENT(IN):: P1, P2 !atom positions of ref. and deformed system
!
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: conffile
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(4):: user_values !user values of NeighFactor and theta_max
LOGICAL:: fileexists !does the file exist?
LOGICAL:: firstref   !is the first file used as reference system? (default: yes)
LOGICAL:: ucref      !is the reference a unit cell? (default: no)
LOGICAL:: usercutoff,userthetamax !did the user provide these values?
INTEGER:: i, j, k, m, n, iat
INTEGER:: nb_neigh, eps
INTEGER:: Nneighbors
INTEGER:: NNmin   !minimum number of neighbours to keep
INTEGER:: Nsites  !number of atomic environments
INTEGER:: ok
INTEGER:: progress      !To show calculation progress
INTEGER:: Q_matrix_rank, INFO, LWORK
INTEGER,DIMENSION(:),ALLOCATABLE:: Nlist !index of neighbours
INTEGER,DIMENSION(:),ALLOCATABLE:: Tab_PQ    !correspondance table between neighbor lists
INTEGER,DIMENSION(:),ALLOCATABLE:: newindex  !list of index after sorting
INTEGER,DIMENSION(:),ALLOCATABLE:: siteindex !for each atom, index of its type of site in Pref
INTEGER,DIMENSION(:),ALLOCATABLE:: sitescore !"score" of reference sites (only when filefirst is a unit cell)
INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: NeighList1, NeighList2  !list of neighbors for systems 1 and 2
REAL(dp):: alpha, alpha_tmp !angles between two vectors
REAL(dp):: theta_max  !if angle between vectors in Pneigh and Qneigh exceed this value, exclude them
REAL(dp):: NeighFactor !%of tolerance in the cutoff for neighbor search
REAL(dp),PARAMETER:: radius=6.d0 !R for neighbour search (6A should be enough to find neighbours in any structure)
REAL(dp):: cutoff      !cutoff for neighbor search
REAL(dp):: d1n, d2n, angle, tempreal
REAL(dp),DIMENSION(3,3):: test_matrix
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries
REAL(dp),DIMENSION(3,3,3):: A_tensor  !tensor A(IM)
REAL(dp),DIMENSION(:),ALLOCATABLE:: Stemp
REAL(dp),DIMENSION(:),ALLOCATABLE:: work_array
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: IdMat, Q_plus
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S, V_NN
REAL(dp),DIMENSION(:,:),POINTER:: Ppoint !pointer to the system to build reference
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P_neigh, Q_neigh, P_matrix, Q_matrix, P_neigh_tmp, Q_matrix_copy
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: PosList1, PosList2 !list of positions of neighbors of one atom
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE:: Delta_G, Delta_G_matrix, Delta_G_tmp
!
!Output variables and arrays
!Pref(i,j,k) = "reference" atomic environments
!  Reference environment #i contains Pref(i,1,6) atoms (central atom+neighbours)
!  atom j=1 is central atom, j=2 to N are neighbours,
!  for atom #j, k=1-3 are (x,y,z) and k=4 atomic number (=0 if no atom in this site)
!  for atom j=1,
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE:: Pref  !"reference" atomic environments
!G(:,:,:) = lattice correspondance tensor, for each atom G is a 3x3 tensor
!i.e. atom #i has a tensor G(i,1:3,1:3)
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(OUT):: G
!Parameters used for the calculation: cutoff, NNmin, NeighFactor, theta_max
REAL(dp),DIMENSION(4),INTENT(INOUT):: params


!Initialize variables and arrays
Nsites = 0
firstref = .TRUE.
ucref = .FALSE.
usercutoff = .FALSE.
userthetamax = .FALSE.
user_values(:) = ""
IF(ALLOCATED(Pref)) DEALLOCATE(Pref)
IF(ALLOCATED(G)) DEALLOCATE(G)
!
!
!
!!!!  Set up NeighFactor and theta_max  !!!!
!Cutoff for neighbour search: an initial value of 6 A
!should be enough to find some neighbors in any system
!This value will be adjusted later
cutoff = 4.d0
!Work with neighbours that are closer than neighbour #NNmin
!NOTE: first neighbours shell count 8 atoms in bcc, 12 atoms in fcc
!     However it counts only 6 atoms in NaCl structure. This will be adjusted later on
NNmin = 6
!After the neighbor search, neighbors will be sorted by increasing distance
!Neighbor #NNmin will be at distance d
!Only neighbors within a distance d*NeighFactor will be kept and used to compute G
!Default value of 1.2 is expected to be robust for simple lattices (fcc, bcc)
NeighFactor = 1.1d0
!Then, neighbours in reference (Pneigh) and in system (Qneigh) will be compared
!If position vectors form an angle greater than theta_max, they will be excluded from calculation
!Hartley and Mishin found that a value of 27° gives best results in FCC lattices,
!i.e. a little less than half the 60° angle between first neighbours
theta_max = 27.d0   !degrees
!
!Check if user wrote specific values for these parameters in a configuration file
!(UNIX/Linux: "atomsk.conf", Windows: "atomsk.ini"), with the format:
! nye cutoff <value>
! nye nnmin <value>
! nye neighfactor <value>
! nye theta_max <value>
#if defined(WINDOWS)
conffile = "atomsk.ini"
#else
conffile = "atomsk.conf"
#endif
INQUIRE(FILE=conffile,EXIST=fileexists)
IF( fileexists ) THEN
  CALL ATOMSK_MSG(16,(/conffile/),(/0.d0/))
  OPEN(UNIT=31,FILE=conffile,STATUS="OLD",FORM="FORMATTED")
  m = 0
  DO
    n=1
    READ(31,'(a)',ERR=110,END=110) temp
    temp = TRIM(ADJUSTL(temp))
    IF( temp(1:1).NE.'#' ) THEN
      IF( StrDnCase(temp(1:3))=="nye" ) THEN
        temp = TRIM(ADJUSTL(temp(4:)))
        IF( StrDnCase(temp(1:6))=="cutoff" ) THEN
          READ(temp(7:),*,ERR=105,END=105) cutoff
          m = m+1
          WRITE(msg,'(f12.3)') cutoff
          user_values(m) = "cutoff = "//TRIM(ADJUSTL(msg))
          usercutoff = .TRUE.
        ELSEIF( StrDnCase(temp(1:5))=="nnmin" ) THEN
          READ(temp(6:),*,ERR=105,END=105) NNmin
          m = m+1
          WRITE(msg,*) NNmin
          user_values(m) = "NNmin = "//TRIM(ADJUSTL(msg))
        ELSEIF( StrDnCase(temp(1:11))=="neighfactor" ) THEN
          READ(temp(12:),*,ERR=105,END=105) NeighFactor
          m = m+1
          WRITE(msg,'(f12.3)') NeighFactor
          user_values(m) = "NeighFactor = "//TRIM(ADJUSTL(msg))
        ELSEIF( StrDnCase(temp(1:9))=="theta_max" ) THEN
          READ(temp(10:),*,ERR=105,END=105) theta_max
          m = m+1
          WRITE(msg,'(f12.3)') theta_max
          user_values(m) = "theta_max = "//TRIM(ADJUSTL(msg))//"°"
          userthetamax = .TRUE.
        ENDIF
      ENDIF
    ENDIF
    n=0
    105 CONTINUE
    IF( n>0 ) CALL ATOMSK_MSG(808,(/TRIM(temp)/),(/0.d0/))
  ENDDO
  110 CONTINUE
  CLOSE(31)
  CALL ATOMSK_MSG(4075,user_values(:),(/0.d0/))
  !
ELSE
  !User did not define values in a file
  !Default values above fail for complex or distorted systems
  !=> Try to adjust them if system contains more than one atom type
  !Search how many atom types exist in the system2
  CALL FIND_NSP(P2(:,4),aentries)
  !Get number of different elements
  n = SIZE(aentries,1)
  IF( n>1 ) THEN
    !There is more than one element in this system
    !Get max. number of atoms for one element, multiply by 10%
    tempreal = 0.1d0*MAXVAL(aentries(:,2))
    IF( n==2 ) THEN
      !Compare relative concentrations of the two elements
      IF( aentries(1,2)>0.1d0*tempreal .AND. aentries(2,2)>0.1*tempreal ) THEN
        !All elements are present in large number => binary material
        !=> the default NeighFactor will probably be too small, use larger value
        IF( DABS(aentries(1,2)-aentries(2,2)) < tempreal ) THEN
          !Both elements are in comparable concentration
          !Assume rock-salt (NaCl, MgO...), where angle between first neighbours is 90°
          NNmin = 3
          NeighFactor = 1.3d0
          theta_max = 43.d0
        ELSE
          !One element is more concentrated than the other, e.g. SiO2, Al2O3...
          NNmin = 3
          NeighFactor = 1.3d0
          theta_max = 50.d0
        ENDIF
      ELSE
        !Otherwise, one element is present only in small concentration (<10%)
        !=> probably a unitary compound with few impurities, don't change NeighFactor
      ENDIF
      !
    ELSE IF( n>=3 ) THEN
      !Good chances that it is a complex material
      !=> boost NeighFactor
      NNmin = 3
      NeighFactor = 1.35d0   !1.33
      theta_max = 55.d0
    ENDIF
  ENDIF
ENDIF
!
!Convert theta_max from degrees to radians
theta_max = DEG2RAD(theta_max)
!
!
!
100 CONTINUE
IF( ALLOCATED(P1) .AND. SIZE(P1,1)>0 ) THEN
  !A "reference" system was provided
  IF( SIZE(P1,1) < SIZE(P2,1) ) THEN
    !The system1 contains fewer atoms than the system2
    !Consider that the first system is a "unit cell"
    firstref = .FALSE.
    ucref = .TRUE.
  ELSE IF( SIZE(P1,1) == SIZE(P2,1) ) THEN
    !Systems 1 and 2 contain exactly the same number of atoms
    !Assume that system 1 is full reference (perfect crystal),
    !and that atoms have the same IDs in system 2
    firstref = .TRUE.
    ucref = .FALSE.
  ELSE
    CALL ATOMSK_MSG(4810,(/""/),(/DBLE(SIZE(P1,1)) , DBLE(SIZE(P2,1))/))
    nerr=nerr+1
    GOTO 1000
  ENDIF
  !
ELSE  !i.e. P1 is NOT ALLOCATED
  !No "reference": we have only the deformed system (P2) to work with
  firstref = .FALSE.
  ucref = .FALSE.
ENDIF
!
!
!Construct neighbor lists
CALL ATOMSK_MSG(11,(/""/),(/0.d0/))
IF( firstref .OR. ucref ) THEN
  IF( .NOT.ALLOCATED(NeighList1) ) THEN
    !Construct neighbor list for reference system
    CALL NEIGHBOR_LIST(H1,P1,radius,NeighList1)
  ENDIF
  IF( verbosity==4 ) THEN
    !Some debug messages
    WRITE(msg,*) "Size of neighbor list for SYSTEM 1: ", SIZE(NeighList1,1), SIZE(NeighList1,2)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,*) "Sample of neighbor list for SYSTEM 1:"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO i=1,MIN(20,SIZE(NeighList1,1))
      WRITE(msg,'(i5,a1,20i5)') i, "|", NeighList1(i,1:MIN(SIZE(NeighList1,2),20))
      IF( SIZE(NeighList1,1)>20 ) msg = TRIM(ADJUSTL(msg))//" (...)"
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
    IF( i>=20 ) THEN
      WRITE(msg,*) '      (...discontinued...)'
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDIF
  ENDIF
ENDIF
!
!Construct neighbor list for studied system
IF( .NOT.ALLOCATED(NeighList2) ) THEN
  CALL NEIGHBOR_LIST(H2,P2,radius,NeighList2)
ENDIF
IF( verbosity==4 ) THEN
  WRITE(msg,*) "Size of neighbor list for SYSTEM 2: ", SIZE(NeighList2,1), SIZE(NeighList2,2)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "Sample of neighbor list for SYSTEM 2:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,MIN(20,SIZE(NeighList2,1))
    WRITE(msg,'(i5,a1,20i5)') i, "|", NeighList2(i,1:MIN(SIZE(NeighList2,2),20))
    msg = TRIM(ADJUSTL(msg))//' (...)'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
  IF( i>=20 ) THEN
    WRITE(msg,*) '      (...discontinued...)'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
ENDIF
CALL ATOMSK_MSG(15,(/""/),(/0.d0/))
!
!
IF( .NOT.firstref ) THEN
  !No reference file was provided (filefirst was "NULL"),
  !or reference is a unit cell
  IF( ucref ) THEN
    !A unit cell was provided: use it to construct reference environments
    Ppoint => P1
    CALL ATOMSK_MSG(4070,(/""/),(/1.d0/))
    !Each atom in the unit cell is considered as a separate "environment"
    !This may create duplicates, but it is necessary in materials where atoms
    !of same species have different environments
    ALLOCATE( Pref(SIZE(Ppoint,1),21,6) )
    Pref(:,:,:) = 0.d0
    !Allocate an array to save, for each atom in system2, the site it belongs to
    ALLOCATE(siteindex(SIZE(P2,1)))
    siteindex(:) = 0
    !Loop on all atoms in the unit cell
    DO iat=1,SIZE(Ppoint,1)
      !Use neighbor list from unit cell, save it in PosList2
      CALL NEIGHBOR_POS(H1,Ppoint,Ppoint(iat,1:3),NeighList1(iat,:),ALLOCATED(NeighList1),radius,PosList2)
      !
      !Now PosList2(:,:) contains the cartesian positions of all neighbors in the cutoff,
      !their distance to the atom #iat, and their indices.
      !Sort them by increasing distance:
      CALL BUBBLESORT(PosList2,4,'up  ',newindex)
      !
      !Get the relative positions of neighbors with respect to atom #iat
      !Keep only neighbors that are within NeighFactor times the distance of NNmin neighbor
      Nneighbors = 0
      j = MIN( NNmin , SIZE(PosList2,1) )
      IF(.NOT.usercutoff) cutoff = NeighFactor * PosList2(j,4)
      DO j=1,MIN(20,SIZE(PosList2,1))
        !Shift all neighbours positions so the central atom is at (0,0,0)
        PosList2(j,1:3) = PosList2(j,1:3) - Ppoint(iat,1:3)
        IF( PosList2(j,4) <= NeighFactor*PosList2(NNmin,4) ) THEN
          !IF( NINT(P2(NINT(PosList2(j,5)),4))==NINT(P2(NINT(PosList2(1,5)),4)) ) THEN
            !Keep this neighbor
            Nneighbors = Nneighbors+1
          !ELSE
            !Ignore this neighbor (not same species as first neighbor)
          !ENDIF
        ELSE
          !Wipe out the rest of the list
          PosList2(j+1:,:) = 0.d0
          EXIT
        ENDIF
      ENDDO
      !
      !Each atom in the unit cell is considered as a "new" environment
      !This may create duplicates, but it is necessary in materials where atoms
      !of same species may have different environments
      !Save the index of the site for atom #iat
      siteindex(iat) = iat
      !Central atom is at Pref(iat,1,1:3) = (0,0,0)
      !Save species of central atom for this site
      Pref(iat,1,4) = P1(iat,4)
      !Set counter of atoms for this type of atom site
      Pref(iat,1,5) = 1.d0
      !Save number of neighbors for this site
      Pref(iat,1,6) = DBLE(Nneighbors)
      !Simply store rel. pos. of neighbors of atom #iat in Pref(iat,:,:)
      DO k=1,Nneighbors
        Pref(iat,k+1,1:3) = PosList2(k,1:3)
        n = NINT(PosList2(k,5))
        Pref(iat,k+1,4) = P1(n,4)
      ENDDO
    ENDDO
    !
  ELSE
    !No reference provided at all: construct reference environments from P2
    Ppoint => P2
    CALL ATOMSK_MSG(4070,(/""/),(/2.d0/))
    !Build "reference" environments by averaging environments found in system 2
    CALL AVG_ENV(H2,P2,NeighList2,Pref,NeighFactor,siteindex)
  ENDIF
  !
  !
  190 CONTINUE
  !Check that arrays Pref *and* siteindex are allocated
  IF( .NOT.ALLOCATED(Pref) .OR. .NOT.ALLOCATED(siteindex) ) THEN
    !No atomic environment found: cannot compute anything
    CALL ATOMSK_MSG(4830,(/""/),(/0.d0/))
    nerr = nerr+1
    GOTO 1000
  ENDIF
  !
  !Count number of different environments that were found
  IF( ucref ) THEN
    Nsites = SIZE(Pref,1)
  ELSE
    Nsites = MAXVAL(siteindex)
  ENDIF
  IF( Nsites<=0 ) THEN
    !No atomic environment found: cannot compute anything
    CALL ATOMSK_MSG(4830,(/""/),(/0.d0/))
    nerr = nerr+1
    GOTO 1000
  ENDIF
  !
  !Check the occurrence of environments, and remove those that have a very low occurrence
  IF( (.NOT.ucref) .AND. Nsites>0 ) THEN
    !Double loop on all types of sites
    DO i=1,SIZE(Pref,1)-1
      DO j=i+1,SIZE(Pref,1)
        !Compare the atomic number of central atoms in sites #i and #j
        IF( NINT(Pref(i,1,4))>0 .AND. NINT(Pref(i,1,4)) == NINT(Pref(j,1,4)) ) THEN
          !Sites i and j have the same central atom
          !Verify if one of them has a very low occurrence, or has very few neighbours
          IF( Pref(i,1,5) <= MAX(10,NINT(Pref(j,1,5)/10.d0)) .OR. Pref(i,1,6)<3 &
            & .OR. .NOT.ANY(NINT(Pref(i,2:,4))>0)  ) THEN
            !Site of type #i has a very low occurrence compared to site #j
            !Add its occurrence to the site type #j
            Pref(j,1,5) = Pref(j,1,5)+Pref(i,1,5)
            !Delete this type of site from the list
            Pref(i,:,:) = 0.d0
            !Update total number of types of sites
            Nsites = Nsites-1
            !Update siteindex: for all atoms in site #i, change site to type #j
            DO k=1,SIZE(siteindex)
              IF(siteindex(k)==i) siteindex = j
            ENDDO
          ELSEIF( Pref(j,1,5) <= MAX(10,NINT(Pref(i,1,5)/10.d0)) .OR. Pref(j,1,6)<3 &
            & .OR. .NOT.ANY(NINT(Pref(j,2:,4))>0)  ) THEN
            !Site of type #j has a very low occurrence compared to site #i
            !Add its occurrence to the site type #i
            Pref(i,1,5) = Pref(i,1,5)+Pref(j,1,5)
            !Delete it from the list
            Pref(j,:,:) = 0.d0
            !Update total number of types of sites
            Nsites = Nsites-1
            !Update siteindex: for all atoms in site #j, change site to type #i
            DO k=1,SIZE(siteindex)
              IF(siteindex(k)==j) siteindex = i
            ENDDO
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  !
  IF( Nsites<=0 ) THEN
    !No atomic environment found: cannot compute anything
    CALL ATOMSK_MSG(4830,(/""/),(/0.d0/))
    nerr = nerr+1
    GOTO 1000
  ENDIF
  !
  !Display number of environments on screen
  CALL ATOMSK_MSG(4071,(/""/),(/DBLE(Nsites)/))
  !
  !If user did not provide computation parameter, adjust parameters according to environments that were found
  IF( .NOT.usercutoff .OR. .NOT.userthetamax ) THEN
    !For each reference site, adjust cut-off to first neighbour shell,
    !and theta_max to half
    DO i=1,SIZE(Pref,1)
      !Compute distance to first neighbours
      DO j=1,SIZE(Pref,2)

      ENDDO
      !Compute angle between neighbours
    ENDDO
  ENDIF
  !
  !Write some debugging information
  IF( verbosity==4 ) THEN
    IF( ALLOCATED(Pref) .AND. SIZE(Pref,1)>0 ) THEN
      !For each site, write position of central atom (that should always be (0,0,0))
      !and positions of neighbors into a XYZ file named "atomsk_site_i.xyz" for visualization
      msg = " ATOMIC  ENVIRONMENTS  USED  TO  COMPUTE  NYE  TENSOR:"
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      DO i=1,MIN(Nsites,SIZE(Pref,1))
        !PRINT*, "Site #", i, ": ", NINT(Pref(i,1,5)), " atoms in this site"
        !Check if number of neighbours for site #i is positive
        IF( NINT(Pref(i,1,6))>0 ) THEN
          WRITE(msg,*) i
          OPEN(UNIT=23,FILE="atomsk_site_"//TRIM(ADJUSTL(msg))//".xyz",FORM="FORMATTED")
          WRITE(23,*) NINT(Pref(i,1,6))+1
          WRITE(23,*) "# Averaged environment for site #", i
          CALL ATOMSPECIES(Pref(i,1,4),species)
          WRITE(msg,*) i
          WRITE(temp,*) NINT(Pref(i,1,6))
          WRITE(23,'(a2,3X,3f16.8)') species, Pref(i,1,1:3)
          WRITE(msg,*) "  Site # "//TRIM(ADJUSTL(msg))//" occupied by "//species// &
                      & " atom, has "//TRIM(ADJUSTL(temp))//" neighbors:"
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          DO j=2,NINT(Pref(i,1,6))+1
            CALL ATOMSPECIES(Pref(i,j,4),species)
            WRITE(23,'(a2,3X,3f16.8)') species, Pref(i,j,1:3)
            WRITE(msg,'(8X,a2,3X,3f12.3)') species, Pref(i,j,1:3)
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          ENDDO
          CLOSE(23)
        ENDIF
      ENDDO
      OPEN(UNIT=23,FILE="atomsk_siteindex.txt",FORM="FORMATTED")
      WRITE(23,*) SIZE(Ppoint,1)
      WRITE(23,*) "# Atomsk AVG_ENV: atom id and type of site they occupy"
      DO i=1,SIZE(siteindex,1)
        WRITE(23,*) i, siteindex(i)
      ENDDO
      CLOSE(23)
    ENDIF
  ENDIF
  !
  !If a unit cell was provided, then for each atom in system2, find its site type
  IF( ucref ) THEN
    siteindex(:) = 0
    IF(ALLOCATED(sitescore)) DEALLOCATE(sitescore)
    ALLOCATE(sitescore(SIZE(Pref,1)))
    !Loop on all atoms in system 2
    DO i=1,SIZE(P2,1)
      sitescore(:) = 0
      !Detect central atom species
      DO j=1,SIZE(Pref,1)
        IF( NINT(P2(i,4)) == NINT(Pref(j,1,4)) ) THEN
          sitescore(j) = sitescore(j)+1
        ENDIF
      ENDDO
      !If all sites have a score equal to zero we are in trouble
      IF( SUM(sitescore(:))<=0 ) THEN
        nerr=nerr+1
        GOTO 1000
      ENDIF
      !If only one site has a positive score, we are done
      !Otherwise (i.e. if several sites have a non-zero score),
      !then we have to use another criterion to discriminate them
      IF( SUM(sitescore(:))>1 ) THEN
        !Try to disciminate sites by comparing number of neighbours

      ENDIF
      !Save site type in array siteindex(:)
      siteindex(i) = MAXLOC(sitescore(:),DIM=1)
    ENDDO
    IF(ALLOCATED(sitescore)) DEALLOCATE(sitescore)
  ENDIF !end if(ucref)
  !
  NULLIFY(Ppoint)
  !
ENDIF
!
!
!
200 CONTINUE
!**********************************************************************************
!         COMPUTE  LATTICE  CORRESPONDENCE  TENSOR  G  FOR  EACH  ATOM
!**********************************************************************************
CALL ATOMSK_MSG(4062,(/""/),(/0.d0/))
!
ALLOCATE(G(SIZE(P2,1),3,3))
G(:,:,:) = 0.d0
!
!First, loop on all atoms to compute the tensor G for each atom
progress=0
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP& PRIVATE(msg,iat,i,j,k,m,n,ok,Nneighbors,nb_neigh,PosList1,PosList2,P_neigh,P_neigh_tmp,Tab_PQ) &
!$OMP& PRIVATE(IdMat,Q_neigh,Q_plus,P_matrix,Q_matrix,Q_matrix_copy,Stemp,test_matrix,V_NN,newindex) &
!$OMP& PRIVATE(Q_matrix_rank,work_array,LWORK,INFO,tempreal,alpha,alpha_tmp)
DO iat=1,SIZE(P2,1)
  !
  progress = progress+1
  !
  IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
  IF(ALLOCATED(PosList1)) DEALLOCATE(PosList1)
  IF(ALLOCATED(PosList2)) DEALLOCATE(PosList2)
  IF(ALLOCATED(P_neigh)) DEALLOCATE(P_neigh)
  !
  IF( SIZE(P2,1) > 5000 ) THEN
    !If there are many atoms, display a fancy progress bar
    CALL ATOMSK_MSG(10,(/""/),(/DBLE(progress),DBLE(SIZE(P2,1))/))
  ENDIF
  !
  !WRITE(msg,*) iat
  !WRITE(msg,*) '==========   ATOM # '//TRIM(ADJUSTL(msg))
  !CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  IF( firstref ) THEN
    !The system in "filefirst" is used as reference
    !Search for neighbors of atom #iat in the first system
    CALL NEIGHBOR_POS(H1,P1(:,:),P1(iat,1:3),NeighList1(iat,:),ALLOCATED(NeighList1),radius,PosList1)
    !
    IF( SIZE(PosList1,1)>=3 ) THEN
      !Now PosList1(:,:) contains the cartesian positions of all neighbors in the radius,
      !their distance to the atom #iat, and their indices.
      !Sort them by increasing distance:
      CALL BUBBLESORT(PosList1,4,'up  ',newindex)
      !
      !Determine cut-off based on distance to neighbour #NNmin
      j = MIN( NNmin , SIZE(PosList1,1) )
      IF(.NOT.usercutoff) cutoff = NeighFactor*PosList1(j,4)
      !
      !Keep only the first neighbors, save them in V_NN
      !Make sure to keep at least NNmin neighbors: compare distances to that of neighbor #NNmin
      !The neighbor list NeighList1 is cleaned up to keep only neighbors
      !(this will speed up things when the Nye tensor is computed below)
      Nneighbors=0
      DO j=1,MIN(20,SIZE(PosList1,1))
        IF( PosList1(j,4) <= cutoff ) THEN
          !This neighbor is about as close as the NNmin-th neighbor => keep it
          Nneighbors=Nneighbors+1
        ENDIF
      ENDDO
      !Clean neighbor list for this atom
      !NeighList1(iat,:) = 0
      ALLOCATE(V_NN(Nneighbors,3))
      V_NN(:,:) = 0.d0
      Nneighbors=0
      m=0
      DO j=1,MIN(20,SIZE(PosList1,1))
        IF( PosList1(j,4) <= cutoff ) THEN
          !This neighbor is closer than the NNmin-th neighbor => keep it
          Nneighbors=Nneighbors+1
          V_NN(Nneighbors,:) = PosList1(j,1:3)
          !Add the index of this neighbor to the neighbor list NeighList1
!           IF( PosList1(j,5).NE.iat .AND. .NOT.ANY(NeighList1(iat,:)==NINT(PosList1(j,5))) ) THEN
!             m=m+1
!             NeighList1(iat,m) = NINT(PosList1(j,5))
!           ENDIF
        ENDIF
      ENDDO
      !We don't need PosList1 for this atom anymore
      IF(ALLOCATED(PosList1)) DEALLOCATE(PosList1)
    ENDIF
    !
    WRITE(msg,*) "SYSTEM 1: atom #", iat, "(", SIZE(V_NN,1), " neighbors)"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !
    IF( .NOT. ALLOCATED(V_NN) .OR. SIZE(V_NN,1)<3 ) THEN
      !Not enough neighbors to perform calculation => skip
      CALL ATOMSK_MSG(4709,(/""/),(/ DBLE(iat) , DBLE(SIZE(V_NN,1)) /))
      nwarn = nwarn+1
      GOTO 290
      !
    ELSEIF ( SIZE(V_NN,1)>100 ) THEN
      !Atom #iat has more than 100 neighbors => skip calculation
      CALL ATOMSK_MSG(4705,(/""/),(/DBLE(iat)/))
      nwarn = nwarn+1
      GOTO 290
      !
    ELSE
      !
      !Save relative positions of neighbors of atom #iat into P_neigh(:,:)
      ALLOCATE (P_neigh(SIZE(V_NN,1),3))
      P_neigh(:,:)=0.d0
      DO j=1,SIZE(V_NN,1)
        P_neigh(j,:) = V_NN(j,1:3)-P1(iat,1:3)
      ENDDO
      !
    ENDIF
    !
    IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
    !
    !
  ELSE  !i.e. if .NOT.firstref
    !Unit cell was provided as reference, or no reference at all
    !The "reference" for atom #iat is constructed on-the-fly from the environments in Pref(:,:,:)
    IF( siteindex(iat)==0 ) THEN
      !The type of site for this atom was not determined before
      !Determine it now
      IF( ucref ) THEN
        !A unit cell was provided as a reference
        !It means that "reference sites" in Pref were constructed by assuming
        !that each atom in the unit cell belongs to a different "site".
        !This means that for a given atom species, several "references" are possible
        !Parse all reference sites that match species of atom #iat, and decide which
        !one is the best fit. The "score" of each reference site is saved in "sitescore"
        !NB: in the case where two (or more) atoms actually have exactly the same
        !environment in the unit cell, they should have the same score, and it doesn't matter
        !which one is taken as "reference"
        !
        IF( .NOT.ALLOCATED(sitescore) ) THEN
          ALLOCATE( sitescore(SIZE(Pref,1)) )
        ENDIF
        sitescore(:) = 0
        !
        !Determine number of neighbours of atom #iat
        !Array S will contain atom species, and atom count for reference and for analyzed system
        !IF(ALLOCATED(S)) DEALLOCATE(S)
        !ALLOCATE(S(SIZE(Pref,1),3))
        !S(:,:) = 0.d0
        Nneighbors = 0
        DO k=1,SIZE(NeighList2,2)  !loop on neighbours of atom #iat
          !If neighbour index is zero there are no more neighbours => exit loop on k
          IF( NeighList2(iat,k)==0 ) EXIT
          !Otherwise, increment number of neighbours of atom #iat
          Nneighbors = Nneighbors+1
          !Add current neighbour #k to the corresponding species counter in S
!           m = 1
!           DO m=1,SIZE(S,1)
!             IF( NINT(NeighList2(iat,k)) == NINT(S(m,1)) ) THEN
!               S(m,2) = S(m,2) + 1.d0
!               EXIT
!             ELSEIF( NINT(S(m,1))==0 ) THEN
!               S(m,1) = P2(iat,4)
!               S(m,2) = S(m,2) + 1.d0
!             ENDIF
!           ENDDO
        ENDDO  !end k
        DO k=1,MIN(Nsites,SIZE(Pref,1))
          IF( NINT(P2(iat,4))==NINT(Pref(k,1,4)) ) THEN
            !Site #k matches the species of atom #iat
            sitescore(k) = sitescore(k) + 1
            !Check if number of neighbours coincide
!             IF( NINT(Pref(k,1,6)) == Nneighbors ) THEN
!               sitescore(k) = sitescore(k) + 1
!             ENDIF
            !Check if species of neighbours coincide

            !Check if relative positions of neighbours coincide

          ENDIF
        ENDDO
        !Save index of site with maximum score
        siteindex(iat) = MAXLOC(sitescore(:),1)
        !Free memory
        IF(ALLOCATED(S)) DEALLOCATE(S)
      ELSE
        !No reference provided at all
        !It means that "reference sites" in Pref were constructed by averaging
        !sites from system 2, i.e. each atom species has its own reference site
        !=> The type of site can be determined simply from atom species
        DO k=1,MIN(Nsites,SIZE(Pref,1))
          IF( NINT(P2(iat,4))==NINT(Pref(k,1,4)) ) THEN
            siteindex(iat) = k
            EXIT
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    IF( siteindex(iat)>0 ) THEN
      !We know that atom #iat occupies a site of type siteindex(iat) in Pref
      !Get number of neighbors for this type of site
      Nneighbors = NINT(Pref(siteindex(iat),1,6))
      !Save neighbors of atom #iat from Pref into P_neigh
      !(i.e. P_neigh will contain the averaged neighbor positions for this type of site)
      IF( Nneighbors>=3 ) THEN
        ALLOCATE(P_neigh(Nneighbors,3))
        P_neigh(:,:) = 0.d0
        DO i=1,Nneighbors
          P_neigh(i,:) = Pref(siteindex(iat),i+1,1:3)
        ENDDO
      ELSE
        !Not enough neighbors for atom #iat
        CALL ATOMSK_MSG(4709,(/""/),(/ DBLE(iat) , DBLE(Nneighbors) /))
        nwarn = nwarn+1
        GOTO 290
      ENDIF
    ELSE
      PRINT*, "WARNING: no site for atom #", iat
      nwarn=nwarn+1
      GOTO 290
    ENDIF
    !
  ENDIF !end if firstref
  !
  !Now the array P_neigh contains the relative coordinates of neighbors of atom #iat in perfect env.
  !Compute positions of neighbors of atom #iat in second system
  CALL NEIGHBOR_POS(H2,P2,P2(iat,1:3),NeighList2(iat,:),ALLOCATED(NeighList2),radius,PosList2)
  !
  IF( SIZE(PosList2,1)>=NNmin ) THEN
    !Now PosList2(:,:) contains the cartesian positions of all neighbors in the radius,
    !and their distance to the atom #i.
    !Sort them by increasing distance:
    CALL BUBBLESORT(PosList2,4,'up  ',newindex)
    !
    !Determine cut-off based on distance to neighbour #NNmin
    j = MIN( NNmin , SIZE(PosList2,1) )
    IF(.NOT.usercutoff) cutoff = NeighFactor*PosList2(j,4)
    !
    IF( firstref ) THEN
      !The system in "filefirst" is used as reference
      !Keep only the first neighbors, save them in V_NN
      !Make sure to keep at least NNmin neighbors: compare distances to that of neighbor #NNmin
      Nneighbors=0
      DO j=1,SIZE(PosList2,1)
        IF( PosList2(j,4) <= cutoff ) THEN
          !This neighbor is about as close as the 3rd neighbor => keep it
          Nneighbors=Nneighbors+1
        ENDIF
      ENDDO
      !Clean neighbor list for this atom
      !NeighList2(iat,:) = 0
      ALLOCATE(V_NN(Nneighbors,3))
      V_NN(:,:) = 0.d0
      Nneighbors=0
      m=0
      DO j=1,SIZE(PosList2,1)
        IF( PosList2(j,4) <= cutoff ) THEN
          !This neighbor is closer than NNmin-th neighbor => keep it
          Nneighbors=Nneighbors+1
          V_NN(Nneighbors,:) = PosList2(j,1:3)
          !Add the index of this neighbor to the neighbor list NeighList2
!           IF( PosList2(j,5).NE.iat .AND. .NOT.ANY(NeighList2(iat,:)==NINT(PosList2(j,5))) ) THEN
!             m=m+1
!             NeighList2(iat,m) = NINT(PosList2(j,5))
!           ENDIF
        ENDIF
      ENDDO
      !
    ELSE  !i.e. if .NOT.firstref
      !Unit cell was provided as reference, or no reference at all
      !The "reference" for atom #iat is constructed on-the-fly from the environments in Pref(:,:,:)
      Nneighbors=0
      DO j=1,SIZE(PosList2,1)
        IF( PosList2(j,4) <= cutoff ) THEN
          !This neighbor is about as close as the NNmin-th neighbor => keep it
          Nneighbors=Nneighbors+1
        ENDIF
      ENDDO
      !Clean neighbor list for this atom
      !NeighList2(iat,:) = 0
      ALLOCATE(V_NN(Nneighbors,3))
      V_NN(:,:) = 0.d0
      Nneighbors=0
      k=0
      DO j=1,SIZE(PosList2,1)
        IF( PosList2(j,4) <= cutoff ) THEN
          !This neighbor is closer than the NNmin-th neighbor => keep it
          Nneighbors=Nneighbors+1
          V_NN(Nneighbors,:) = PosList2(j,1:3)
          !Add the index of this neighbor to the neighbor list NeighList2
!           IF( PosList2(j,5).NE.iat .AND. .NOT.ANY(NeighList2(iat,:)==NINT(PosList2(j,5))) ) THEN
!             k=k+1
!             NeighList2(iat,k) = NINT(PosList2(j,5))
!           ENDIF
        ENDIF
      ENDDO
    ENDIF  !end if firstref
  ENDIF   !end if size(PosList2)>=3
  !We don't need PosList2 for this atom anymore
  IF(ALLOCATED(PosList2)) DEALLOCATE(PosList2)
  !
  WRITE(msg,*) "SYSTEM 2: atom #", iat, "(", SIZE(V_NN,1), " neighbors)"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !
  IF( .NOT. ALLOCATED(V_NN) .OR. SIZE(V_NN,1)<3 ) THEN
    !Not enough neighbors to perform calculation => skip
    CALL ATOMSK_MSG(4709,(/""/),(/ DBLE(iat) , DBLE(SIZE(V_NN,1)) /))
    nwarn = nwarn+1
    !
  ELSEIF( SIZE(V_NN,1)>100 ) THEN
    CALL ATOMSK_MSG(4705,(/""/),(/DBLE(iat)/))
    nwarn = nwarn+1
    !
  ELSE
    !
    !Save positions of neighbors into the array Q_neigh
    IF(ALLOCATED(Q_neigh)) DEALLOCATE(Q_neigh)
    ALLOCATE (Q_neigh(SIZE(V_NN,1),3))
    Q_neigh(:,:)=0.d0
    DO j=1,SIZE(V_NN,1)
      Q_neigh(j,:) = V_NN(j,1:3) - P2(iat,1:3)
    ENDDO
    IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
    !
    !Build a table of correspondance between indices in P_neigh and those in Q_neigh
    IF(ALLOCATED(Tab_PQ)) DEALLOCATE(Tab_PQ)
    ALLOCATE (Tab_PQ(SIZE(Q_neigh,1)))
    Tab_PQ(:)=0
    DO j=1,SIZE(Q_neigh,1)
      alpha=100
      DO k=1,SIZE(P_neigh,1)
        alpha_tmp=DABS(ANGVEC(Q_neigh(j,1:3),P_neigh(k,1:3)))
        IF( alpha_tmp < alpha ) THEN
          Tab_PQ(j)=k
          alpha=alpha_tmp
        ENDIF
      ENDDO
    ENDDO
    !
    !Matrix P_matrix and Q_matrix
    IF(ALLOCATED(P_neigh_tmp)) DEALLOCATE(P_neigh_tmp)
    ALLOCATE(P_neigh_tmp(SIZE(Q_neigh,1),3))
    P_neigh_tmp(:,:)=0.d0
    DO j=1,SIZE(Q_neigh,1)
      P_neigh_tmp(j,1:3)=P_neigh(Tab_PQ(j),1:3)
    ENDDO
    !
    !At this stage, the tables P and Q are in the same order
    !Compute angles between corresponding vectors, and
    !exclude vectors if the angle is too large
    Tab_PQ(:)=0
    nb_neigh=0
    DO j=1,SIZE(Q_neigh,1)
      alpha=DABS(ANGVEC(Q_neigh(j,1:3),P_neigh_tmp(j,1:3)))
      IF (alpha.gt.theta_max) THEN
        Tab_PQ(j)=0
      ELSE
        ok=0
        !boucle sur les voisins
        DO k=1,SIZE(P_neigh,1)
          alpha_tmp=DABS(ANGVEC(Q_neigh(j,1:3),P_neigh(k,1:3)))
          IF (alpha_tmp.eq.alpha) ok=ok+1
        ENDDO
        IF (ok.eq.1) THEN
          nb_neigh=nb_neigh+1
          Tab_PQ(j)=nb_neigh
        ELSE
          Tab_PQ(j)=0
        ENDIF
      ENDIF
    ENDDO
    !
    !Save final positions of neighbors into P_matrix and Q_matrix
    IF(ALLOCATED(P_matrix)) DEALLOCATE(P_matrix)
    IF(ALLOCATED(Q_matrix)) DEALLOCATE(Q_matrix)
    ALLOCATE(Q_matrix(nb_neigh,3))
    ALLOCATE(P_matrix(nb_neigh,3))
    Q_matrix(:,:)=0.d0
    P_matrix(:,:)=0.d0
    DO j=1,SIZE(Q_neigh,1)
      IF (Tab_PQ(j).NE.0) THEN
        Q_matrix(Tab_PQ(j),:) = Q_neigh(j,:)
        P_matrix(Tab_PQ(j),:) = P_neigh_tmp(j,:)
      ENDIF
    ENDDO
    !
    !
    IF (verbosity==4) THEN
      WRITE(msg,*) '-----Relative positions of neighbors-----'
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      DO j=1,MIN(SIZE(P_neigh,1),SIZE(Q_neigh,1))
        WRITE(msg,'(a9,3f12.4,a6,f12.4)') 'P_neigh: ', P_neigh(j,1:3), " | d =", VECLENGTH(P_neigh(j,1:3))
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDDO
      DO j=1,SIZE(Q_neigh,1)
        WRITE(msg,'(a9,3f12.4,a6,f12.4)') 'Q_neigh: ', Q_neigh(j,1:3), " | d =", VECLENGTH(Q_neigh(j,1:3))
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDDO
      WRITE(msg,*) '-----------------------------------------'
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDIF
    !
    !
    IF( nb_neigh < 3 ) THEN
      !Not enough neighbors to perform calculation => skip
      WRITE(msg,*) 'NOT ENOUGH NEIGHBORS: nb_neigh = ', nb_neigh
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      CALL ATOMSK_MSG(4709,(/""/),(/ DBLE(iat) , DBLE(nb_neigh) /))
      nwarn = nwarn+1
      IF(ALLOCATED(Tab_PQ)) DEALLOCATE(Tab_PQ)
      !
    ELSE
      !
      IF(ALLOCATED(Tab_PQ)) DEALLOCATE(Tab_PQ)
      IF(ALLOCATED(P_neigh)) DEALLOCATE(P_neigh)
      IF(ALLOCATED(Q_neigh)) DEALLOCATE(Q_neigh)
      !
      IF (verbosity==4) THEN
        WRITE(msg,*) '-----Neighbors used for calculation of G----'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        DO j=1,nb_neigh
          WRITE(msg,'(3f12.4,a3,3f12.4)') P_matrix(j,1:3), " | ", Q_matrix(j,1:3)
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        ENDDO
        WRITE(msg,*) '--------------------------------------------'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDIF
      !
      !Compute Q+ = (Q^T * Q)^-1 * Q^T   (Eq.18)
      !This is done thanks to the LAPACK subroutine DGELSS
      ALLOCATE(IdMat(nb_neigh,nb_neigh))
      IdMat(:,:) = 0.d0
      DO j=1,nb_neigh
        IdMat(j,j) = 1.d0
      ENDDO
      !
      ALLOCATE(Q_matrix_copy(nb_neigh,3))
      Q_matrix_copy(:,:) = Q_matrix(:,:)
      !
      LWORK=3*MIN(nb_neigh,3) + MAX( 2*MIN(nb_neigh,3), MAX(nb_neigh,3), nb_neigh )
      IF(ALLOCATED(work_array)) DEALLOCATE(work_array)
      IF(ALLOCATED(Stemp)) DEALLOCATE(Stemp)
      ALLOCATE (work_array(MAX(1,LWORK)))
      ALLOCATE (Stemp(3))
      !
      CALL DGELSS(nb_neigh,3,nb_neigh,Q_matrix,nb_neigh,IdMat,MAX(nb_neigh,3),  &
                  & Stemp,-1.d0,Q_matrix_rank,work_array,LWORK,INFO)
      !
      IF ( INFO.NE.0 ) THEN
        !LAPACK routine returned an error
        CALL ATOMSK_MSG(4710,(/"inverse of Q"/),(/ 0.d0 /))
        nwarn = nwarn+1
        !
      ELSE
        !Keep only the first 3 rows of output IdMat
        IF(ALLOCATED(Q_plus)) DEALLOCATE(Q_plus)
        ALLOCATE(Q_plus(3,nb_neigh))
        Q_plus(:,:)=0.d0
        DO j=1,3
          Q_plus(j,1:nb_neigh) = IdMat(j,1:nb_neigh)
        ENDDO
        !
        !Test if  (Q+ * Qmatrix)  is identity matrix
        test_matrix(:,:) = MATMUL(Q_plus,Q_matrix_copy)
        IF (verbosity==4) THEN
          WRITE(msg,*)  'Q+ * Qmatrix should be id matrix, check:'
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          WRITE(msg,*) test_matrix(1,:)
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          WRITE(msg,*) test_matrix(2,:)
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          WRITE(msg,*) test_matrix(3,:)
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        ENDIF
        ok=0
        DO i=1,3
          DO j=1,3
            tempreal = test_matrix(i,j)
            IF( i==j ) THEN
              tempreal = tempreal-1.d0
            ENDIF
            IF( DABS(tempreal)>1.d-6 ) THEN
              ok = 1
            ENDIF
          ENDDO
        ENDDO
        !
        IF(ALLOCATED(IdMat)) DEALLOCATE(IdMat)
        IF(ALLOCATED(Q_matrix)) DEALLOCATE(Q_matrix)
        IF(ALLOCATED(Q_matrix_copy)) DEALLOCATE(Q_matrix_copy)
        !
        !
        IF( ok.NE.0 ) THEN
          WRITE(msg,*) iat
          CALL ATOMSK_MSG(4717,(/"Q+ * Qmatrix (atom #"//TRIM(ADJUSTL(msg))//")"/),(/0.d0/))
          nwarn = nwarn+1
          !
        ELSE
          !
          !Compute  G = Q+ * P   (Eq.17)
          !Save the final tensor G for atom #iat
          G(iat,:,:) = MATMUL(Q_plus,P_matrix)
        ENDIF
        !
        !
        IF (verbosity==4) THEN
          WRITE(msg,*) 'TENSOR G FOR ATOM # ', iat
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          WRITE(msg,*) "    ", G(iat,1,:)
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          WRITE(msg,*) "    ", G(iat,2,:)
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          WRITE(msg,*) "    ", G(iat,3,:)
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          WRITE(msg,*) '============================================'
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        ENDIF
        !
        IF(ALLOCATED(P_matrix)) DEALLOCATE(P_matrix)
        IF(ALLOCATED(Q_plus)) DEALLOCATE(Q_plus)
        IF(ALLOCATED(work_array)) DEALLOCATE(work_array)
        IF(ALLOCATED(Stemp)) DEALLOCATE(Stemp)
        !
        !
      ENDIF !end IF(INFO.NE.0)
      !
    ENDIF  !end IF(nb_neigh < 3)
    !
  ENDIF  !end IF(.NOT.ALLOCATED(V_NN))... for system 2
  !
  290 CONTINUE
  !
ENDDO
!$OMP END PARALLEL DO
!
!
1000 CONTINUE
!At this point tensor G was computed for all atoms (hopefully), and is returned by this routine
!Save parameters used during calculation in array param(:)
params(1) = cutoff
params(2) = NNmin
params(3) = NeighFactor
params(4) = theta_max
!
!
END SUBROUTINE COMPUTE_G
!
!
END MODULE cmpt_G
