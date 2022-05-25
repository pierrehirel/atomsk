MODULE mode_nye
!
!**********************************************************************************
!*  MODE_NYE                                                                      *
!**********************************************************************************
!* This module reads two sets of atomic coordinates and                           *
!* computes the corresponding Nye tensor.                                         *
!* The algorithm follows the method described in:                                 *
!*    C.S. Hartley, Y. Mishin, Acta Mater. 53 (2005) 1313                         *
!* Equation numbers also refer to this reference.                                 *
!* The names "system1" and "system2" refer to the reference system                *
!* (bulk with no defect) and the system to analyze, respectively.                 *
!* The reference system 1 can be either:                                          *
!* - a full reference system containing the same number of atoms as system2,      *
!*   and indexed in the exact same order;
!* - a unit cell of the material in system 2;                                     *
!* - NULL, in which case the "reference" environments are built on-the-fly        *
!*   by averaging atomic environments found in the system2.                       *
!**********************************************************************************
!* (C) October 2013 - Philippe Carrez                                             *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 10 March 2022                                    *
!**********************************************************************************
!* OUTLINE:                                                                       *
!* 100        Read atom positions systems 1 and 2, construct neighbor lists       *
!* 200        Compute lattice correspondence tensor G for each atom               *
!* 300        Using G, compute Nye tensor for each atom                           *
!* 400        Output final results to file(s)                                     *
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
!Module for reading input files
USE readin
USE writeout
!Module for applying options
USE options
!
!
!
CONTAINS
!
!
SUBROUTINE NYE_TENSOR(filefirst,filesecond,options_array,prefix,outfileformats)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: filefirst, filesecond, prefix
CHARACTER(LEN=2):: species
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE,INTENT(IN):: outfileformats !list of output file formats
CHARACTER(LEN=4096):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
LOGICAL:: firstref  !is the first file used as reference system? (default: yes)
LOGICAL:: lowp      !does an environment have a low probability?
LOGICAL:: ucref     !is the reference a unit cell? (default: no)
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: i, j, k, m, n, iat
INTEGER:: nb_neigh, eps
INTEGER:: Nneighbors
INTEGER:: Nsites  !number of atomic environments
INTEGER:: ok
INTEGER:: progress      !To show calculation progress
INTEGER:: Q_matrix_rank, INFO, LWORK
INTEGER,DIMENSION(:),ALLOCATABLE:: Nlist !index of neighbours
INTEGER,DIMENSION(:),ALLOCATABLE:: Tab_PQ    !correspondance table between neighbor lists
INTEGER,DIMENSION(:),ALLOCATABLE:: newindex  !list of index after sorting
INTEGER,DIMENSION(:),ALLOCATABLE:: siteindex !for each atom, index of its type of site in Pref
INTEGER,DIMENSION(:,:),ALLOCATABLE:: NeighList1, NeighList2  !list of neighbors for systems 1 and 2
REAL(dp):: alpha, alpha_tmp !angles between two vectors
REAL(dp):: beta, gamma      !
REAL(dp),PARAMETER:: maxangvec=27.d0*(2.d0*pi/360.d0)  !if angle between vectors in Pneigh and Qneigh exceed
                                                    !this value, exclude them (Hartley&Mishin used 27°)
REAL(dp):: NeighFactor !%of tolerance in the radius for neighbor search
REAL(dp):: P1, P2, P3  !coordinates of a neighbor
REAL(dp),PARAMETER:: angle_th=pi/5.d0  !threshold to exlude neighbors
REAL(dp),PARAMETER:: radius=8.d0 !R for neighbor search: 8 A should be enough to find some neighbors in any system
REAL(dp):: tempreal
REAL(dp),DIMENSION(3,3):: alpha_tensor, test_matrix
REAL(dp),DIMENSION(3,3):: Huc   !Box vectors of unit cell (unknown, set to 0 here)
REAL(dp),DIMENSION(3,3):: Hfirst, Hsecond !Box vectors of systems 1 and 2
REAL(dp),DIMENSION(3,3):: Hneigh, Hnew    !Base vectors for neighbors
REAL(dp),DIMENSION(3,3):: ORIENT     !crystallographic orientation (not used here but necessary for some calls)
REAL(dp),DIMENSION(3,3):: rot_matrix !a rotation matrix
REAL(dp),DIMENSION(3,3,3):: A_tensor  !tensor A(IM)
REAL(dp),DIMENSION(9,9):: C_tensor    !stiffness tensor (not used here but necessary for some calls)
REAL(dp),DIMENSION(:),ALLOCATABLE:: Stemp
REAL(dp),DIMENSION(:),ALLOCATABLE:: work_array
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: IdMat, Q_plus
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX  !final auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE,TARGET:: Pfirst, Psecond
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S,V_NN
REAL(dp),DIMENSION(:,:),POINTER:: Ppoint !pointer to the system to build reference
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P_neigh, Q_neigh, P_matrix, Q_matrix, P_neigh_tmp, Q_matrix_copy
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: PosList1, PosList2 !list of positions of neighbors of one atom
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE:: Pref  !references for atoms environments
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE:: G, Delta_G, Delta_G_matrix, Delta_G_tmp 


!Initialize variables and arrays
Nsites = 0
firstref = .TRUE.
ucref = .FALSE.
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
NeighFactor = 1.1d0 !distance to 1st neighbor times 1.1 ensures to exclude second neighbor sphere in bcc metals
                    !1.1 is expected to be robust for simple lattices (fcc, bcc) but fails for complex
                    !or distorted systems
Huc(:,:) = 0.d0
ORIENT(:,:) = 0.d0     !crystallographic orientation (not used, set to zero)
 C_tensor(:,:) = 0.d0  !stiffness tensor (not used, set to zero)
!
!
CALL ATOMSK_MSG(4061,(/""/),(/0.d0/))
!
!
!
100 CONTINUE
!**********************************************************************************
!                            READ  INPUT  FILES
!**********************************************************************************
IF( filefirst=="NULL" ) THEN
  !There is no reference file: atoms environment will be generated from the system2
  firstref = .FALSE.
ELSE
  !Read atomic positions from filefirst and store them into Pfirst(:,:)
  CALL READ_AFF(filefirst,Hfirst,Pfirst,S,comment,AUXNAMES,AUX)
  !Get rid of shells and auxiliary properties
  IF (ALLOCATED(S)) DEALLOCATE(S)
  IF (ALLOCATED(comment)) DEALLOCATE(comment)
  IF (ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
  IF (ALLOCATED(AUX)) DEALLOCATE(AUX)
  !Apply options to system 1
  CALL OPTIONS_AFF(options_array,Huc,Hfirst,Pfirst,S,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
ENDIF
!
!Read atomic positions from filesecond and store them into Psecond(:,:)
CALL READ_AFF(filesecond,Hsecond,Psecond,S,comment,AUXNAMES,AUX)
!Get rid of shells and auxiliary properties
IF (ALLOCATED(S)) DEALLOCATE(S)
IF (ALLOCATED(comment)) DEALLOCATE(comment)
IF (ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF (ALLOCATED(AUX)) DEALLOCATE(AUX)
!Apply options to system 2
CALL OPTIONS_AFF(options_array,Huc,Hsecond,Psecond,S,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
!
!
!Search how many atom types exist in the system2
CALL FIND_NSP(Psecond(:,4),aentries)
IF( SIZE(aentries,1)==2 ) THEN
  !Binary material
  !=> the default NeighFactor will probably be too small, use larger value
  NeighFactor = 1.15d0
ELSE IF( SIZE(aentries,1)>=3 ) THEN
  !Good chances that it is a complex material
  !=> boost NeighFactor
  NeighFactor = 1.33d0   !1.33
ENDIF
!
!
IF( firstref ) THEN
  !Check that systems 1 and 2 have the same number of atoms
  IF( SIZE(Pfirst,1) < SIZE(Psecond,1) ) THEN
    !The system1 contains fewer atoms than the system2
    !Consider that the first system is a "unit cell"
    firstref = .FALSE.
    ucref = .TRUE.
  ELSE IF( SIZE(Pfirst,1) .NE. SIZE(Psecond,1) ) THEN
    CALL ATOMSK_MSG(4810,(/""/),(/DBLE(SIZE(Pfirst,1)) , DBLE(SIZE(Psecond,1))/))
    nerr=nerr+1
    GOTO 1000
  ENDIF
ENDIF
!
!
!Construct neighbor lists
CALL ATOMSK_MSG(11,(/""/),(/0.d0/))
IF( firstref .OR. ucref ) THEN
  !Construct neighbor list for reference system
  CALL NEIGHBOR_LIST(Hfirst,Pfirst,radius,NeighList1)
  IF( verbosity==4 ) THEN
    !Some debug messages
    WRITE(msg,*) "Size of neighbor list for SYSTEM 1: ", SIZE(NeighList1,1), SIZE(NeighList1,2)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,*) "Sample of neighbor list for SYSTEM 1:"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO i=1,MIN(20,SIZE(NeighList1,1))
      WRITE(msg,'(i5,a1,20i5)') i, "|", NeighList1(i,1:MIN(SIZE(NeighList1,2),20))
      msg = TRIM(ADJUSTL(msg))//' (...)'
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
CALL NEIGHBOR_LIST(Hsecond,Psecond,radius,NeighList2)
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
  !No reference file was provided (filefirst was "NULL")
  IF( ucref ) THEN
    !A unit cell was provided: use it to construct reference environments
    Ppoint => Pfirst
    CALL ATOMSK_MSG(4070,(/filefirst/),(/1.d0/))
  ELSE
    !No reference provided at all: construct reference environments from Psecond
    Ppoint => Psecond
    CALL ATOMSK_MSG(4070,(/filesecond/),(/2.d0/))
  ENDIF
  !Build a reference by averaging environments found in system 2
  !It is assumed that there are no more than 20 different atom types,
  !each with less than 20 neighbors. For each of these atoms, we store 6 numbers:
  !x, y, z, atomic number, number of atoms in such site, number of neighbors of such site
  ALLOCATE( Pref(20,21,6) )
  Pref(:,:,:) = 0.d0
  !Allocate an array to save, for each atom in system2, the site it belongs to
  ALLOCATE(siteindex(SIZE(Psecond,1)))
  siteindex(:) = 0
  Nsites = 0  !to count number of different atom sites
  !Loop on atoms in reference system: perform average only on the first 1000 atoms should be sufficient
  DO iat=1,SIZE(Ppoint,1)
    !Get positions of neighbors of this atom
    IF( ucref ) THEN
      !Use neighbor list from unit cell, save it in PosList2
      CALL NEIGHBOR_POS(Hfirst,Ppoint,Ppoint(iat,1:3),NeighList1(iat,:),ALLOCATED(NeighList1),radius,PosList2)
    ELSE
      !Use neighbor list from system2, save it in PosList2
      CALL NEIGHBOR_POS(Hsecond,Ppoint,Ppoint(iat,1:3),NeighList2(iat,:),ALLOCATED(NeighList2),radius,PosList2)
    ENDIF
    !
    !Now PosList2(:,:) contains the cartesian positions of all neighbors in the radius,
    !their distance to the atom #iat, and their indices.
    !Sort them by increasing distance:
    CALL BUBBLESORT(PosList2,4,'up  ',newindex)
    !Get the relative positions of neighbors with respect to atom #iat
    !Keep only neighbors that are within NeighFactor times the distance of third neighbor
    Nneighbors = 0
    DO j=1,MIN(20,SIZE(PosList2,1))
      PosList2(j,1:3) = PosList2(j,1:3) - Ppoint(iat,1:3)
      IF( PosList2(j,4) <= NeighFactor*PosList2(3,4) ) THEN
        !IF( NINT(Psecond(NINT(PosList2(j,5)),4))==NINT(Psecond(NINT(PosList2(1,5)),4)) ) THEN
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
    !Find in Pref the index for this type of atom site
    !or create a new type of site
    i=0
    DO WHILE(i<SIZE(Pref,1))
      i=i+1
      IF( NINT(Ppoint(iat,4))==NINT(Pref(i,1,4)) .AND. ABS(NINT(Pref(i,1,6))-Nneighbors)<=1 ) THEN
        !Atom species and number of neighbors match
        !=> atom #iat occupies a site of the type #i
        !Save the index of the site for atom #iat
        siteindex(iat) = i
        !Increment counter of atoms for this type of atom site
        Pref(i,1,5) = Pref(i,1,5) + 1.d0
        EXIT
      ELSE IF( NINT(Pref(i,1,4))==0 .AND. MOD(Nneighbors,2)==0 .AND. MOD(Nneighbors,5)>0 ) THEN
        !No suitable site was found before: if number of neighbors is even,
        !then create a new one in an empty slot (i.e. where Pref(i,1,4)=0)
        !Save species of central atom for this site
        Pref(i,1,4) = Ppoint(iat,4)
        !Increment counter of atoms for this type of atom site
        Pref(i,1,5) = Pref(i,1,5) + 1.d0
        !Save number of neighbors for this site
        Pref(i,1,6) = DBLE(Nneighbors)
        !Save the index of the site for atom #iat
        siteindex(iat) = i
        EXIT
      ENDIF
      IF(i>SIZE(Pref,1)) THEN
        nerr=nerr+1
        CALL ATOMSK_MSG(4070,(/""/),(/DBLE(i)/))
        GOTO 1000
      ENDIF
    ENDDO
    !
    !Add atom to average only if number of neighbors is greater than 3, lower than 14, even, and
    !not multiple of 5 (because an atom in any material is not supposed to have 5 or 10 neighbors)
    !i.e. we consider only atoms with 4, 6, 8, 12, 14 neighbors are in "perfect" environments
    IF( Nneighbors>3 .AND. Nneighbors<=14 .AND. MOD(Nneighbors,2)==0 .AND. MOD(Nneighbors,5)>0 ) THEN
      !Add neighbors positions to already known positions for averaging
      IF( NINT(Pref(i,1,5))==1 ) THEN
        !No neighbor was found for this site yet
        !Simply store rel. pos. of neighbors of atom #iat in Pref(i,:,:)
        DO k=1,Nneighbors
          Pref(i,k+1,1:3) = PosList2(k,1:3)
          Pref(i,k+1,4) = Ppoint(NINT(PosList2(k,5)),4)
        ENDDO
      ELSE
        !Some atoms were already found in this site
        !For each neighbor in PosList, find its best match in Pref 
        !(i.e. the one that maximizes the dot product) and add its position to it
        DO k=1,Nneighbors
          n=0
          alpha = 0.d0
          j=1
          DO WHILE( j<=SIZE(Pref,2) )
            tempreal = DOT_PRODUCT( PosList2(k,1:3) , Pref(i,j,1:3) )
            IF( tempreal > alpha ) THEN
              !This neighbor position is a better match
              n = j
              alpha = tempreal
            ENDIF
            j=j+1
          ENDDO
          IF (n>0) THEN
            !Now n is the index of the best matching neighbor in Pref
            !Add position and perform averaging
            Pref(i,n,1:3) = ( Pref(i,n,5)*Pref(i,n,1:3) + PosList2(k,1:3) ) / (Pref(i,n,5)+1.d0)
            !Save atomic number of this neighbor (only if it's empty)
            IF( NINT(Pref(i,n,4))==0 ) THEN
              Pref(i,n,4) = Ppoint(NINT(PosList2(k,5)),4)
            ENDIF
            !Increment number of neighbors at this position
            Pref(i,n,5) = Pref(i,n,5) + 1.d0
          ELSE
            !n is zero, meaning that the dot product was zero for all neighbors
            !This should not happen, but here we are
            !The site for atom #iat will be decided later (see beginning of section about G tensor)
            !PRINT*, "         ERROR no matching site for neighbor #", k
          ENDIF
          !
        ENDDO
      ENDIF
      !
    ENDIF  !end if MOD(Nneighbors,2)==0
    !
  ENDDO ! end loop on iat
  !
  190 CONTINUE
  Nsites = MAXVAL(siteindex)
  IF( Nsites<=0 ) THEN
    !No atomic environment found: cannot compute anything
    CALL ATOMSK_MSG(4830,(/""/),(/0.d0/))
    nerr = nerr+1
    GOTO 1000
  ENDIF
  !
!   !Nsites atom sites were found
!   !For each site, perform a rotation so that a neighbors's position vector
!   !is aligned with the Cartesian X axis, another is in XY plane,
!   !and a third one pointing out of XY plane.
!   !Then, express all neighbors positions in this new base
!   DO i=1,MIN(Nsites,SIZE(Pref,1))
!     Hneigh(:,:) = 0.d0
!     !First neighbor is used as first base vector
!     Hneigh(1,1:3) = Pref(i,2,1:3)
!     m = 1
!     !Find two other vectors to form a base of 3 linearly independent vectors
!     DO j=3,NINT(Pref(i,1,6))+1  !loop on other neighbors
!       !Compute magnitude of cross product between current vector and first base vector
!       P1 = VECLENGTH( CROSS_PRODUCT( Pref(i,j,1:3) , Hneigh(1,1:3) ) )
!       IF( m>=2 ) THEN
!         !Second base vector was found:
!         !Compute magnitude of cross product between current vector and second base vector
!         P2 = VECLENGTH( CROSS_PRODUCT( Pref(i,j,1:3) , Hneigh(2,1:3) ) )
!       ELSE
!         !Second base vector not found yet: don't enforce this condition
!         P2 = 100.d0
!       ENDIF
!       !If current vector is not aligned with an existing "base vector" (cross product.NE.0),
!       !then save it as an additional base vector
!       IF( m<3 .AND. P1>1.d0 .AND. P2>1.d0 ) THEN
!         m = m+1
!         Hneigh(m,1:3) = Pref(i,j,1:3)
!       ENDIF
!     ENDDO
!     !
!     IF( m==3 ) THEN
!       !Now Hneigh contains the "old base vectors"
!       !Convert all neighbor positions into fractional coordinates
!       CALL CART2FRAC(Pref(i,2:,1:3),Hneigh)
!       !Convert base vectors into "conventional notation"
!       P1 = VECLENGTH(Hneigh(1,:))
!       P2 = VECLENGTH(Hneigh(2,:))
!       P3 = VECLENGTH(Hneigh(3,:))
!       alpha = ANGVEC(Hneigh(2,:),Hneigh(3,:))
!       beta  = ANGVEC(Hneigh(3,:),Hneigh(1,:))
!       gamma = ANGVEC(Hneigh(1,:),Hneigh(2,:))
!       !Then convert them into vectors
!       CALL CONVMAT(P1,P2,P3,alpha,beta,gamma,Hneigh)
!       !Convert all neighbor positions back into Cartesian coordinates
!       CALL FRAC2CART(Pref(i,2:,1:3),Hneigh)
!       !
!     ELSE
!       !m is not equal to 3, meaning that we could not find a suitable base
!       !of 3 vectors that are linearly independent: just leave neighbor list as it is
!     ENDIF
!   ENDDO
  !
  !Check the occurrence of environments, and remove those that have a very low occurrence
  IF( Nsites>0 ) THEN
    !Double loop on all types of sites
    DO i=1,SIZE(Pref,1)-1
      DO j=i+1,SIZE(Pref,1)
        !Compare the atomic number of central atoms in sites #i and #j
        IF( NINT(Pref(i,1,4))>0 .AND. NINT(Pref(i,1,4)) == NINT(Pref(j,1,4)) ) THEN
          !Sites i and j have the same central atom
          !Verify if one of them has a very low occurrence
          IF( Pref(i,1,5) <= MAX(10,NINT(Pref(j,1,5)/10.d0)) ) THEN
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
          ELSEIF( Pref(j,1,5) <= MAX(10,NINT(Pref(i,1,5)/10.d0)) ) THEN
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
  !Now Pref contains the averaged relative positions of neighbors for each type of atom site
  IF( ucref ) THEN
    !Index of first atoms in system2 may not match those of unit cell
    !Wipe out the "siteindex" array, atom sites will be found later
    siteindex(:) = 0
  ENDIF
  !Write some debugging information
  IF( verbosity==4 ) THEN
    !For each site, write position of central atom (that should always be (0,0,0))
    !and positions of neighbors into a XYZ file named "site_i.xyz" for visualization
    msg = " ATOMS ENVIRONMENTS DETECTED:"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO i=1,MIN(Nsites,SIZE(Pref,1))
      !PRINT*, "Site #", i, ": ", NINT(Pref(i,1,5)), " atoms in this site"
      IF( NINT(Pref(i,1,4))>0 ) THEN
        WRITE(msg,*) i
        OPEN(UNIT=23,FILE="atomsk_site_"//TRIM(ADJUSTL(msg))//".xyz",FORM="FORMATTED")
        WRITE(23,*) NINT(Pref(i,1,6))+1
        WRITE(23,*) "# Averaged environment for site #", i
        CALL ATOMSPECIES(Pref(i,1,4),species)
        WRITE(23,'(a2,3X,3f16.8)') species, Pref(i,1,1:3)
        WRITE(msg,*) "  Site #", i, "occupied by "//species//" atom, has ", NINT(Pref(i,1,6)), "neighbors:"
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        DO j=2,NINT(Pref(i,1,6))+1
          IF( NINT(Pref(i,j,4))>0 ) THEN
            CALL ATOMSPECIES(Pref(i,j,4),species)
          ELSE
            species = "H "
          ENDIF
          WRITE(23,'(a2,3X,3f16.8)') species, Pref(i,j,1:3)
          WRITE(msg,'(8X,3f12.3)') Pref(i,j,1:3)
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        ENDDO
        CLOSE(23)
      ENDIF
    ENDDO
    OPEN(UNIT=23,FILE="atomsk_atomsites.xyz",FORM="FORMATTED")
    WRITE(23,*) SIZE(Psecond,1)
    WRITE(23,*) "# Atom index and type of site they occupy"
    DO i=1,SIZE(Psecond,1)
      WRITE(23,*) i, siteindex(i)
    ENDDO
    CLOSE(23)
  ENDIF
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
ALLOCATE(G(SIZE(Psecond,1),3,3))
G(:,:,:) = 0.d0
!
!First, loop on all atoms to compute the tensor G for each atom
progress=0
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP& PRIVATE(msg,iat,i,j,k,m,n,ok,Nneighbors,nb_neigh,PosList1,PosList2,P_neigh,P_neigh_tmp,Tab_PQ) &
!$OMP& PRIVATE(IdMat,Q_neigh,Q_plus,P_matrix,Q_matrix,Q_matrix_copy,Stemp,test_matrix,V_NN,newindex) &
!$OMP& PRIVATE(Q_matrix_rank,work_array,LWORK,INFO,tempreal,alpha,alpha_tmp)
DO iat=1,SIZE(Psecond,1)
  !
  progress = progress+1
  !
  IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
  IF(ALLOCATED(PosList1)) DEALLOCATE(PosList1)
  IF(ALLOCATED(PosList2)) DEALLOCATE(PosList2)
  IF(ALLOCATED(P_neigh)) DEALLOCATE(P_neigh)
  !
  IF( SIZE(Psecond,1) > 5000 ) THEN
    !If there are many atoms, display a fancy progress bar
    CALL ATOMSK_MSG(10,(/""/),(/DBLE(progress),DBLE(SIZE(Psecond,1))/))
  ENDIF
  !
  WRITE(msg,*) iat
  WRITE(msg,*) '==========   ATOM # '//TRIM(ADJUSTL(msg))
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  IF( firstref ) THEN
    !The system in "filefirst" is used as reference
    !Search for neighbors of atom #iat in the first system
    CALL NEIGHBOR_POS(Hfirst,Pfirst(:,:),Pfirst(iat,1:3),NeighList1(iat,:),ALLOCATED(NeighList1),radius,PosList1)
    !
    IF( SIZE(PosList1,1)>=3 ) THEN
      !Now PosList1(:,:) contains the cartesian positions of all neighbors in the radius,
      !their distance to the atom #iat, and their indices.
      !Sort them by increasing distance:
      CALL BUBBLESORT(PosList1,4,'up  ',newindex)
      !
      !Keep only the first neighbors, save them in V_NN
      !Make sure to keep at least 3 neighbors: compare distances to that of the 3rd neighbor
      !The neighbor list NeighList1 is cleaned up to keep only neighbors
      !(this will speed up things when the Nye tensor is computed below)
      Nneighbors=0
      DO j=1,MIN(20,SIZE(PosList1,1))
        IF( PosList1(j,4) <= NeighFactor*PosList1(3,4) ) THEN
          !This neighbor is about as close as the 3rd neighbor => keep it
          Nneighbors=Nneighbors+1
        ENDIF
      ENDDO
      !Clean neighbor list for this atom
      NeighList1(iat,:) = 0
      ALLOCATE(V_NN(Nneighbors,3))
      V_NN(:,:) = 0.d0
      Nneighbors=0
      m=0
      DO j=1,MIN(20,SIZE(PosList1,1))
        IF( PosList1(j,4) <= NeighFactor*PosList1(3,4) ) THEN
          !This neighbor is closer than the 3rd neighbor => keep it
          Nneighbors=Nneighbors+1
          V_NN(Nneighbors,:) = PosList1(j,1:3)
          !Add the index of this neighbor to the neighbor list NeighList1
          IF( PosList1(j,5).NE.iat .AND. .NOT.ANY(NeighList1(iat,:)==NINT(PosList1(j,5))) ) THEN
            m=m+1
            NeighList1(iat,m) = NINT(PosList1(j,5))
          ENDIF
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
        P_neigh(j,:) = V_NN(j,1:3)-Pfirst(iat,1:3)
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
      !Determine it now, based on atom species
      DO k=1,MIN(Nsites,SIZE(Pref,1))
        IF( NINT(Psecond(iat,4))==NINT(Pref(k,1,4)) ) THEN
          siteindex(iat) = k
          EXIT
        ENDIF
      ENDDO
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
  CALL NEIGHBOR_POS(Hsecond,Psecond,Psecond(iat,1:3),NeighList2(iat,:),ALLOCATED(NeighList2),radius,PosList2)
  !
  IF( SIZE(PosList2,1)>=3 ) THEN
    !Now PosList2(:,:) contains the cartesian positions of all neighbors in the radius,
    !and their distance to the atom #i.
    !Sort them by increasing distance:
    CALL BUBBLESORT(PosList2,4,'up  ',newindex)
    !
    IF( firstref ) THEN
      !The system in "filefirst" is used as reference
      !Keep only the first neighbors, save them in V_NN
      !Make sure to keep at least 3 neighbors: compare distances to that of the 3rd neighbor
      Nneighbors=0
      DO j=1,SIZE(PosList2,1)
        IF( PosList2(j,4) <= NeighFactor*PosList2(3,4) ) THEN
          !This neighbor is about as close as the 3rd neighbor => keep it
          Nneighbors=Nneighbors+1
        ENDIF
      ENDDO
      !Clean neighbor list for this atom
      NeighList2(iat,:) = 0
      ALLOCATE(V_NN(Nneighbors,3))
      V_NN(:,:) = 0.d0
      Nneighbors=0
      m=0
      DO j=1,SIZE(PosList2,1)
        IF( PosList2(j,4) <= NeighFactor*PosList2(3,4) ) THEN
          !This neighbor is closer than the 3rd neighbor => keep it
          Nneighbors=Nneighbors+1
          V_NN(Nneighbors,:) = PosList2(j,1:3)
          !Add the index of this neighbor to the neighbor list NeighList1
          IF( PosList2(j,5).NE.iat .AND. .NOT.ANY(NeighList2(iat,:)==NINT(PosList2(j,5))) ) THEN
            m=m+1
            NeighList2(iat,m) = NINT(PosList2(j,5))
          ENDIF
        ENDIF
      ENDDO
      !
    ELSE  !i.e. if .NOT.firstref
      !Unit cell was provided as reference, or no reference at all
      !The "reference" for atom #iat is constructed on-the-fly from the environments in Pref(:,:,:)
      Nneighbors=0
      DO j=1,SIZE(PosList2,1)
        IF( PosList2(j,4) <= NeighFactor*PosList2(3,4) ) THEN
          !This neighbor is about as close as the 3rd neighbor => keep it
          Nneighbors=Nneighbors+1
        ENDIF
      ENDDO
      !Clean neighbor list for this atom
      NeighList2(iat,:) = 0
      ALLOCATE(V_NN(Nneighbors,3))
      V_NN(:,:) = 0.d0
      Nneighbors=0
      k=0
      DO j=1,SIZE(PosList2,1)
        IF( PosList2(j,4) <= NeighFactor*PosList2(3,4) ) THEN
          !This neighbor is closer than the 3rd neighbor => keep it
          Nneighbors=Nneighbors+1
          V_NN(Nneighbors,:) = PosList2(j,1:3)
          !Add the index of this neighbor to the neighbor list NeighList2
          IF( PosList2(j,5).NE.iat .AND. .NOT.ANY(NeighList2(iat,:)==NINT(PosList2(j,5))) ) THEN
            k=k+1
            NeighList2(iat,k) = NINT(PosList2(j,5))
          ENDIF
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
      Q_neigh(j,:) = V_NN(j,1:3) - Psecond(iat,1:3)
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
      IF (alpha.gt.maxangvec) THEN
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
        WRITE(msg,'(a9,3f12.4,a5,f12.4)') 'P_neigh: ', P_neigh(j,1:3), "| d =", VECLENGTH(P_neigh(j,1:3))
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDDO  
      DO j=1,SIZE(Q_neigh,1)
        WRITE(msg,'(a9,3f12.4,a5,f12.4)') 'Q_neigh: ', Q_neigh(j,1:3), "| d =", VECLENGTH(Q_neigh(j,1:3))
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
IF(nerr>0) GOTO 1000
!
IF( verbosity==4 ) THEN
  !DEBUG: Write each component of G tensor into data file (projected in XY plane)
  !If user did not provide output file name or prefix, set a default one
  msg = prefix
  i = 0
  IF( LEN_TRIM(prefix)<=0 ) THEN
    msg = filesecond
    i = SCAN(msg,".",BACK=.TRUE.)
  ENDIF
  IF(i==0) i=LEN_TRIM(msg)+1
  msg = TRIM(ADJUSTL(msg(1:i-1)))//"_G.dat"
  OPEN(UNIT=41,FILE=TRIM(msg),STATUS="UNKNOWN",FORM="FORMATTED")
  WRITE(41,*) "# Lattice distortion tensor G computed with Atomsk"
  DO i=1,SIZE(G,1)
    msg = ""
    DO j=1,3
      WRITE(msg,'(a,3(2X,f12.6))') TRIM(msg), (G(i,j,k),k=1,3)
    ENDDO
    WRITE(41,*) Psecond(i,1), Psecond(i,2), "  "//TRIM(ADJUSTL(msg))
  ENDDO
  CLOSE(41)
ENDIF
!
!Write atom coordinates and per-atom matrix G into a CFG file for visualization
IF( ALLOCATED(AUX) ) DEALLOCATE(AUX)
IF( ALLOCATED(AUXNAMES) ) DEALLOCATE(AUXNAMES)
IF( ALLOCATED(comment) ) DEALLOCATE(comment)
ALLOCATE( AUX( SIZE(Psecond,1) , 9 ) )
AUX(:,:) = 0.d0
AUX(:,1) = G(:,1,1)
AUX(:,2) = G(:,1,2)
AUX(:,3) = G(:,1,3)
AUX(:,4) = G(:,2,1)
AUX(:,5) = G(:,2,2)
AUX(:,6) = G(:,2,3)
AUX(:,7) = G(:,3,1)
AUX(:,8) = G(:,3,2)
AUX(:,9) = G(:,3,3)
ALLOCATE (AUXNAMES(9))
AUXNAMES(1)="G_11"
AUXNAMES(2)="G_12"
AUXNAMES(3)="G_13"
AUXNAMES(4)="G_21"
AUXNAMES(5)="G_22"
AUXNAMES(6)="G_23"
AUXNAMES(7)="G_31"
AUXNAMES(8)="G_32"
AUXNAMES(9)="G_33"
ALLOCATE(comment(1))
comment(1) = "# Per-atom G tensor computed by Atomsk"
!If user did not provide output file name or prefix, set a default one
msg = prefix
IF( LEN_TRIM(prefix)<=0 ) THEN
  msg = filesecond
ENDIF
i = SCAN(msg,".",BACK=.TRUE.)
IF(i==0) i=LEN_TRIM(msg)
msg = TRIM(ADJUSTL(msg(1:i-1)))//"_G.cfg"
CALL WRITE_AFF(msg,outfileformats,Hsecond,Psecond,S,comment,AUXNAMES,AUX)
DEALLOCATE(AUX)
DEALLOCATE(AUXNAMES)
DEALLOCATE(comment)
!
!
!
300 CONTINUE
!**********************************************************************************
!                            COMPUTE  NYE  TENSOR
!**********************************************************************************
!Tensor G is known for all atoms => use it to compute the Nye tensor
CALL ATOMSK_MSG(4063,(/""/),(/0.d0/))
!
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
ALLOCATE( AUX(SIZE(Psecond,1),9) )
AUX(:,:) = 0.d0
!
!New Loop on the atoms
progress=0
DO iat=1,SIZE(Psecond,1)
  !
  progress = progress+1
  !
  !Initialize Nye tensor for atom #iat
  alpha_tensor(:,:) = 0.d0
  !
  IF( SIZE(Psecond,1) > 5000 ) THEN
    !If there are many atoms, display a fancy progress bar
    CALL ATOMSK_MSG(10,(/""/),(/DBLE(progress),DBLE(SIZE(Psecond,1))/))
  ENDIF
  !
  A_tensor(:,:,:) = 0.d0
  IF(ALLOCATED(Nlist)) DEALLOCATE(Nlist)
  IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
  IF(ALLOCATED(P_neigh)) DEALLOCATE(P_neigh)
  IF(ALLOCATED(Q_neigh)) DEALLOCATE(Q_neigh)
  IF(ALLOCATED(Delta_G)) DEALLOCATE(Delta_G)
  !
  !
  IF( firstref ) THEN
    !The system in "filefirst" is used as reference
    !Search neighbors of atom #iat in first system
    CALL NEIGHBOR_POS(Hfirst,Pfirst,Pfirst(iat,1:3),NeighList1(iat,:),ALLOCATED(NeighList1),radius,PosList1)
    !
    IF( SIZE(PosList1,1)>=3 ) THEN
      !Now PosList1(:,:) contains the cartesian positions of all neighbors of atom #iat,
      !their distance to the central atom #iat, and their indices in P(:,:).
      !Sort them by increasing distance:
      CALL BUBBLESORT(PosList1,4,'up  ',newindex)
      !Keep only the first neighbors, save them in V_NN
      !Make sure to keep at least 3 neighbors: compare distances to that of the 3rd neighbor
      Nneighbors=0
      DO j=1,SIZE(PosList1,1)
        IF( PosList1(j,4) <= NeighFactor*PosList1(3,4) ) THEN
          !This neighbor is about as close as the 3rd neighbor => keep it
          Nneighbors=Nneighbors+1
        ENDIF
      ENDDO
      !Save first neighbors positions into V_NN(:,:)
      !Save their index into Nlist(:)
      ALLOCATE(V_NN(Nneighbors,3))
      V_NN(:,:) = 0.d0
      ALLOCATE(Nlist(Nneighbors))
      Nlist(:) = 0
      Nneighbors=0
      DO j=1,SIZE(PosList1,1)
        IF( PosList1(j,4) <= NeighFactor*PosList1(3,4) ) THEN
          !This neighbor is closer than the 3rd neighbor => keep it
          Nneighbors=Nneighbors+1
          V_NN(Nneighbors,:) = PosList1(j,1:3)
          Nlist(Nneighbors) = NINT(PosList1(j,5))
        ENDIF
      ENDDO
    ENDIF
    !
    !
    IF( .NOT. ALLOCATED(V_NN) .OR. SIZE(V_NN,1)<3 ) THEN
      !Not enough neighbors to perform calculation => skip
      CALL ATOMSK_MSG(4709,(/""/),(/ DBLE(iat) , DBLE(SIZE(V_NN,1)) /))
      nwarn = nwarn+1
      GOTO 390
      !
    ELSEIF( SIZE(V_NN,1)>100 ) THEN
      !Atom #iat has more than 100 neighbors => skip calculation
      CALL ATOMSK_MSG(4705,(/""/),(/DBLE(iat)/))
      nwarn = nwarn+1
      GOTO 390
      !
    ELSE
      !
      !Save relative positions of neighbors in P_neigh
      !and compute  Delta_G(IM) = G_neighbor(IM) - G_0(IM)  (Eq.19)
      IF(ALLOCATED(P_neigh)) DEALLOCATE(P_neigh)
      ALLOCATE (P_neigh(SIZE(V_NN,1),3))
      P_neigh(:,:)=0.d0
      ALLOCATE (Delta_G(SIZE(V_NN,1),3,3))
      Delta_G(:,:,:) = 0.d0
      !
      DO j=1,SIZE(V_NN,1)
        P_neigh(j,:) = V_NN(j,1:3) - Pfirst(iat,1:3)
        !
        DO i=1,3
          DO m=1,3
            Delta_G(j,i,m) = G(Nlist(j),i,m) - G(iat,i,m)
          ENDDO
        ENDDO
        !
      ENDDO
      !
    ENDIF
    !
    IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
    !
    !Search neighbors of atom #iat in second system
    CALL NEIGHBOR_POS(Hsecond,Psecond,Psecond(iat,1:3),NeighList2(iat,:),ALLOCATED(NeighList2),radius,PosList2)
    !
    IF( SIZE(PosList2,1)>=3 ) THEN
      !Now PosList2(:,:) contains the cartesian positions of all neighbors in the radius,
      !their distance to the atom #iat, and their indices.
      !Sort them by increasing distance:
      CALL BUBBLESORT(PosList2,4,'up  ',newindex)
      !
      !Keep only the first neighbors, save them in V_NN
      !Make sure to keep at least 3 neighbors: compare distances to that of the 3rd neighbor
      Nneighbors=0
      DO j=1,SIZE(PosList2,1)
        IF( PosList2(j,4) <= NeighFactor*PosList2(3,4) ) THEN
          !This neighbor is about as close as the 3rd neighbor => keep it
          Nneighbors=Nneighbors+1
        ENDIF
      ENDDO
      ALLOCATE(V_NN(Nneighbors,3))
      V_NN(:,:) = 0.d0
      Nneighbors=0
      DO j=1,SIZE(PosList2,1)
        IF( PosList2(j,4) <= NeighFactor*PosList2(3,4) ) THEN
          !This neighbor is closer than the 3rd neighbor => keep it
          Nneighbors=Nneighbors+1
          V_NN(Nneighbors,:) = PosList2(j,1:3)
        ENDIF
      ENDDO
    ENDIF
    !
  ELSE
    !No reference was provided (filefirst was NULL)
    !The atom #iat is in a site of type siteindex(iat)
    !and neighbors of atom #iat are in Pref
    IF( siteindex(iat)>0 ) THEN
      !Get number of neighbors for this type of site
      Nneighbors = NINT(Pref(siteindex(iat),1,6))
      !
      !Save atoms in perfect environment in P_neigh
      ALLOCATE(P_neigh(Nneighbors,3))
      P_neigh(:,:) = 0.d0
      DO j=1,Nneighbors
        P_neigh(j,1:3) = Pref(siteindex(iat),j+1,1:3)
      ENDDO
      !
      !Search neighbors of atom #iat in second system
      CALL NEIGHBOR_POS(Hsecond,Psecond,Psecond(iat,1:3),NeighList2(iat,:),ALLOCATED(NeighList2),radius,PosList2)
      !
      IF( SIZE(PosList2,1)>=3 ) THEN
        !Now PosList2(:,:) contains the cartesian positions of all neighbors of atom #iat,
        !their distance to the central atom #iat, and their indices in Psecond(:,:).
        !Sort them by increasing distance:
        CALL BUBBLESORT(PosList2,4,'up  ',newindex)
        !Keep only the first neighbors, save them in V_NN
        !Make sure to keep at least Nneighbors: compare distances to that of the 3rd neighbor
        nb_neigh=0
        DO j=1,SIZE(PosList2,1)
          IF( PosList2(j,4) <= NeighFactor*PosList2(3,4) ) THEN
            !This neighbor is about as close as the 3rd neighbor => keep it
            nb_neigh=nb_neigh+1
          ENDIF
        ENDDO
        ALLOCATE(V_NN(nb_neigh,5))
        V_NN(:,:) = 0.d0
        ALLOCATE(Nlist(nb_neigh))
        Nlist(:) = 0
        nb_neigh = 0
        DO j=1,SIZE(PosList2,1)
          IF( PosList2(j,4) <= NeighFactor*PosList2(3,4) ) THEN
            !This neighbor is about as close as the 3rd neighbor => keep it
            nb_neigh=nb_neigh+1
            !Save neighbors from second system in V_NN
            V_NN(nb_neigh,1:3) = PosList2(j,1:3)
            !Save index of neighbors in Nlist
            Nlist(nb_neigh) = NINT(PosList2(j,5))
          ENDIF
        ENDDO
        !
        !Compute  Delta_G(IM) = G_neighbor(IM) - G_0(IM)  (Eq.19)
        IF( nb_neigh>=3 ) THEN
          ALLOCATE (Delta_G(nb_neigh,3,3))
          Delta_G(:,:,:) = 0.d0
          DO j=1,nb_neigh
            DO i=1,3
              DO m=1,3
                Delta_G(j,i,m) = G(Nlist(j),i,m) - G(iat,i,m)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          !Not enough neighbors => skip calculation
          GOTO 390
        ENDIF
      ELSE
        !Not enough neighbors => skip calculation
        GOTO 390
      ENDIF
      !
    ELSE
      !No site was found for this atom; at this point, several warnings were already displayed
      !Just skip the calculation for this atom
      GOTO 390
    ENDIF
    !
  ENDIF  !end if firstref
  !
  IF( verbosity==4 ) THEN
    WRITE(msg,*) "atom #", iat, "  (", nb_neigh, " neighbors), Delta_G ="
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO j=1,nb_neigh
      WRITE(msg,'(3(3f9.3,a3))') Delta_G(j,1,1:3), " | ", Delta_G(j,2,1:3), " | ", Delta_G(j,3,1:3)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
  ENDIF
  !
  !
  IF( .NOT. ALLOCATED(V_NN) .OR. SIZE(V_NN,1)<3 ) THEN
    !Not enough neighbors to perform calculation => skip
    CALL ATOMSK_MSG(4709,(/""/),(/ DBLE(iat) , DBLE(SIZE(V_NN,1)) /))
    nwarn = nwarn+1
    !
  ELSEIF( SIZE(V_NN,1)>100 ) THEN
    !Atom #iat has more than 100 neighbors => skip calculation
    CALL ATOMSK_MSG(4705,(/""/),(/DBLE(iat)/))
    nwarn = nwarn+1
    !
  ELSE
    !
    !Save relative positions of neighbors in Q_neigh
    IF(ALLOCATED(Q_neigh)) DEALLOCATE(Q_neigh)
    ALLOCATE (Q_neigh(SIZE(V_NN,1),3))
    Q_neigh(:,:)=0.d0
    DO j=1,SIZE(V_NN,1)
      Q_neigh(j,:) = V_NN(j,1:3) - Psecond(iat,1:3)
    ENDDO
    DEALLOCATE(V_NN)
    !
    !Generate correspondance table between P_neigh and Q_neigh
    ALLOCATE (Tab_PQ(SIZE(Q_neigh,1)))
    Tab_PQ(:)=0
    !
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
    IF(ALLOCATED(P_neigh_tmp)) DEALLOCATE(P_neigh_tmp)
    IF(ALLOCATED(Delta_G_tmp)) DEALLOCATE(Delta_G_tmp)
    ALLOCATE(P_neigh_tmp(SIZE(Q_neigh,1),3))
    ALLOCATE(Delta_G_tmp(SIZE(Q_neigh,1),3,3))
    P_neigh_tmp(:,:)=0.d0
    Delta_G_tmp(:,:,:)=0.d0
    !
    IF( firstref ) THEN
      !A full reference was provided
      DO j=1,SIZE(Q_neigh,1)
        P_neigh_tmp(j,1:3) = P_neigh(Tab_PQ(j),1:3)
        Delta_G_tmp(j,:,:) = Delta_G(Tab_PQ(j),:,:)
      ENDDO
    ELSE
      !No reference was provided
      DO j=1,SIZE(Q_neigh,1)
        P_neigh_tmp(j,1:3) = P_neigh(Tab_PQ(j),1:3)
        Delta_G_tmp(j,:,:) = Delta_G(j,:,:)
      ENDDO
    ENDIF
    DEALLOCATE(Tab_PQ)
    !
    !We won't need Delta_G anymore
    DEALLOCATE(Delta_G)
    !
    !At this stage, the tables P and Q are in the same order
    !Compute angles between corresponding vectors, and
    !exclude vectors if the angle is too large
    ALLOCATE (Tab_PQ(SIZE(Q_neigh,1)))
    Tab_PQ(:)=0
    nb_neigh=0
    DO j=1,SIZE(Q_neigh,1)
      alpha = DABS(ANGVEC(Q_neigh(j,1:3),P_neigh_tmp(j,1:3)))
      IF( alpha > maxangvec ) THEN
        Tab_PQ(j)=0
      ELSE
        ok=0 
        !Loop on neighbors
        DO k=1,SIZE(P_neigh,1)
          alpha_tmp=DABS(ANGVEC(Q_neigh(j,1:3),P_neigh(k,1:3)))
          IF(alpha_tmp==alpha) ok=1  !=ok+1
        ENDDO
        IF( ok==1 ) THEN
          nb_neigh=nb_neigh+1
          Tab_PQ(j)=nb_neigh
        ELSE
          Tab_PQ(j)=0
        ENDIF
      ENDIF
    ENDDO
    !
    !Compute Delta_G(IM) = Q * A(IM)   (Eq.20)
    IF(ALLOCATED(P_matrix)) DEALLOCATE(P_matrix)
    IF(ALLOCATED(Q_matrix)) DEALLOCATE(Q_matrix)
    ALLOCATE(Q_matrix(nb_neigh,3))
    Q_matrix(:,:)=0.d0
    IF(ALLOCATED(Delta_G_matrix)) DEALLOCATE(Delta_G_matrix)
    ALLOCATE(Delta_G_matrix(nb_neigh,3,3))
    Delta_G_matrix(:,:,:)=0.d0
    DO j=1,SIZE(Q_neigh,1)
      IF (Tab_PQ(j).NE.0) THEN
        Q_matrix(Tab_PQ(j),:) = Q_neigh(j,:)
        Delta_G_matrix(Tab_PQ(j),:,:) = Delta_G_tmp(j,:,:)
      ENDIF
    ENDDO
    !
    DEALLOCATE(P_neigh,Q_neigh,P_neigh_tmp,Delta_G_tmp)
    DEALLOCATE(Tab_PQ)
    !
    !
    IF ( nb_neigh < 3 ) THEN
      !Not enough neighbors to perform calculation => skip
      CALL ATOMSK_MSG(4709,(/""/),(/ DBLE(iat) , DBLE(nb_neigh) /))
      nwarn = nwarn+1
      !
    ELSE
      !
      !Compute Q+ = (Q^T * Q)^-1 * Q^T   (Eq.18)
      !This is done thanks to the LAPACK subroutine DGELSS
      IF(ALLOCATED(IdMat)) DEALLOCATE(IdMat)
      ALLOCATE(IdMat(nb_neigh,nb_neigh))
      IdMat(:,:) = 0.d0
      DO j=1,nb_neigh
        IdMat(j,j) = 1.d0
      ENDDO
      !
      LWORK=3*MIN(nb_neigh,3) + MAX( 2*MIN(nb_neigh,3), MAX(nb_neigh,3), nb_neigh )
      ALLOCATE (work_array(MAX(1,LWORK)))
      ALLOCATE (Stemp(3))
      !
      CALL DGELSS(nb_neigh,3,nb_neigh,Q_matrix,nb_neigh,IdMat,MAX(nb_neigh,3),    &
                  & Stemp,-1.d0,Q_matrix_rank,work_array,LWORK,INFO)
      !
      IF (INFO.NE.0) THEN
        CALL ATOMSK_MSG(4710,(/"inverse of Q"/),(/ 0.d0 /))
        nwarn = nwarn+1
        !
      ELSE
        !
        !keep only the first 3 rows of output IdMat
        ALLOCATE(Q_plus(3,nb_neigh))
        Q_plus(:,:)=0.d0
        DO j=1,3
          Q_plus(j,1:nb_neigh) = IdMat(j,1:nb_neigh)
        ENDDO
        DEALLOCATE(IdMat,Q_matrix)
        !
        !Compute A(IM) = Q+ * Delta_G(IM)   (Eq.21)
        DO i=1,3
          DO m=1,3
            A_tensor(i,m,:) = MATMUL(Q_plus,Delta_G_matrix(:,i,m))
          ENDDO
        ENDDO
        !
        DEALLOCATE (Delta_G_matrix,Q_plus,work_array,Stemp)
        !
        !
        !Compute the Nye tensor:  alpha(jk) = -epsilon(imk) * T(ijm)   (Eq.22)
        !where epsilon is the Levi-Civita permutation symbol, and
        !T(ijm) is the tensor of derivatives of G
        !NOTE: Eq.22 is wrong in the published article by Hartley et al.
        DO j=1,3
          DO k=1,3
            alpha_tensor(j,k) = 0.d0
            DO i=1,3
              DO m=1,3
                eps = EPS_LEVI_CIVITA(i,m,k)
                alpha_tensor(j,k) = alpha_tensor(j,k) - DBLE(eps)*A_tensor(i,j,m)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        !
        !Save the values of Nye tensor into AUX(:,:)
        DO i=1,9
          IF (i<=3) THEN
            j=1
            k=i
          ELSEIF (i<=6) THEN
            j=2
            k=i-3
          ELSE
            j=3
            k=i-6
          ENDIF
          AUX(iat,i) = alpha_tensor(j,k)
        ENDDO
        !
        !
      ENDIF  !end IF(INFO.NE.0)
      !
    ENDIF  !end IF(nb_neigh < 3)
    !
  ENDIF  !end IF(.NOT.ALLOCATED(V_NN))... for system 2
    !
  !ENDIF  !end IF(.NOT.ALLOCATED(V_NN))... for system 1
  !
  390 CONTINUE
  !
ENDDO ! End loop on all atoms iat
!
!
!
400 CONTINUE
!Write final results into file(s)
!
!DEALLOCATE (Pfirst)
ALLOCATE (AUXNAMES(9))
!
AUXNAMES(1)="Nye_11"
AUXNAMES(2)="Nye_12"
AUXNAMES(3)="Nye_13"
AUXNAMES(4)="Nye_21"
AUXNAMES(5)="Nye_22"
AUXNAMES(6)="Nye_23"
AUXNAMES(7)="Nye_31"
AUXNAMES(8)="Nye_32"
AUXNAMES(9)="Nye_33"
!
ALLOCATE(comment(1))
 comment(1) = "# Per-atom Nye tensor computed by Atomsk"
!
msg = prefix
!If user did not provide output file name or prefix, set a default one
IF( LEN_TRIM(prefix)<=0 ) THEN
  msg = filesecond
  i = SCAN(msg,".",BACK=.TRUE.)
  IF(i==0) i=LEN_TRIM(msg)
  msg = TRIM(ADJUSTL(msg(1:i-1)))//"_Nye"
  !In addition, if user did not specify an output format, use CFG by default
  IF( .NOT.ALLOCATED(outfileformats) .OR. SIZE(outfileformats)==0 ) THEN
    msg = TRIM(ADJUSTL(msg))//".cfg"
  ENDIF
ENDIF
!
!CALL WRITE_AFF(prefix,outfileformats,Hfirst,Pfirst,S,comment,AUXNAMES,AUX)
CALL WRITE_AFF(msg,outfileformats,Hsecond,Psecond,S,comment,AUXNAMES,AUX)
!
DEALLOCATE (Psecond, AUX, AUXNAMES)
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE NYE_TENSOR
!
!
END MODULE mode_nye
