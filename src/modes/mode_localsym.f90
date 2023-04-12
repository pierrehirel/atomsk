MODULE mode_localsym
!
!**********************************************************************************
!*  MODE_LOCALSYM                                                                 *
!**********************************************************************************
!* This module computes a symmetry parameter for each atom in the system,         *
!* and saves it as an auxiliary parameter in the array AUX.                       *
!* Such a parameter can be useful for visualizing defects.                        *
!* For unary systems the method employed follows the well-established             *
!* centrosymmetry by C.L. Kelchner and co-workers:                                *
!*   C.L. Kelchner et al., Phys. Rev. B 58 (1998) 11085-8                         *
!* Note that this method is already implemented and ready to use in the           *
!* visualization software Atomeye:                                                *
!*   http://li.mit.edu/Archive/Graphics/A/#central_symm_coloring                  *
!*   http://li.mit.edu/Archive/Graphics/A/Doc/CentralSymmetry.pdf                 *
!* In Atomsk, the actual implementation for unary systems (bcc, fcc...)           *
!* follows the algorithm of Ju Li's Atomeye.                                      *
!* For binary and ternary systems a variant of this algorithm is used.            *
!* In addition to centrosymmetry, this mode also computes:                        *
!* - a "tetrahedral parameter", which is naught if atom is in tetrahedral         *
!*   environment, and positive otherwise (useful e.g. in diamond lattice).        *
!* - a "sp2 parameter", which is naught if atom is in a planar sp2 environment    *
!*   (e.g. graphene), and positive otherwise.                                     *
!* In the end, only the smallest of these parameters is saved as                  *
!* a "local symmetry parameter".                                                  *
!**********************************************************************************
!* (C) April 2015 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 16 Jan. 2023                                     *
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
USE subroutines
USE readin
USE options
USE writeout
!
!
CONTAINS
!
!
SUBROUTINE LOCAL_SYM(inputfile,options_array,prefix,outfileformats)
!
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=*),INTENT(IN):: prefix
!
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: msg, temp
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of output file formats
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES, newAUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: c     !column in AUX that will contain the central symmetry parameters c(i)
INTEGER:: i, ipairs, j, k, ktemp
INTEGER:: Mdefault      !most common number of neighbors
INTEGER:: sp1, sp2      !atomic numbers of two atoms
INTEGER,PARAMETER:: Nmax=14 !maximum number of neighbors for any lattice
INTEGER:: Nneigh1           !number of first-neighbors of an atom (even)
INTEGER,DIMENSION(3):: Nneigh  !number of 1st, 2nd, and 3rd neighbors of an atom
INTEGER,DIMENSION(3):: Neighsp !species of 1st, 2nd, and 3rd neighbors of an atom
INTEGER:: Nspecies      !number of species among an atom's neighbours
INTEGER,DIMENSION(:,:),ALLOCATABLE:: NeighList !list of index of neighbors
INTEGER,DIMENSION(:),ALLOCATABLE:: newindex  !list of index after sorting
INTEGER,DIMENSION(:,:),ALLOCATABLE:: pairs  !indexes of pairs of atoms
REAL(dp):: distance, dmin
REAL(dp):: sum_dj       !sum of |d_j|²
REAL(dp),DIMENSION(3,3):: Huc   !Box vectors of unit cell (unknown, set to 0 here)
REAL(dp),DIMENSION(3,3):: H       !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystallographic orientation of the system (mode create)
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:),ALLOCATABLE:: pairs_distances  !distances
REAL(dp),DIMENSION(:),ALLOCATABLE:: three_angles     !angles
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries    !array containing result
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX, newAUX !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, S     !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: PosList
!
!
!Initialize variables
Mdefault = 12
Nspecies = 0
Huc(:,:) = 0.d0
 C_tensor(:,:) = 0.d0
!
!
msg = 'ENTERING LOCAL_SYM...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(4069,(/inputfile/),(/0.d0/))
!
!
100 CONTINUE
! Read the input file
CALL READ_AFF(inputfile,H,P,S,comment,AUXNAMES,AUX)
!
! Apply options if any
CALL OPTIONS_AFF(options_array,Huc,H,P,S,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
IF(nerr>0) GOTO 1000
!
!If auxiliary properties already exist, add a column to save the central symmetry parameter
IF( ALLOCATED(AUX) .AND. SIZE(AUX,2)==SIZE(AUXNAMES) ) THEN
  !If "local_symmetry" already exists in AUX, use its column (i.e. overwrite)
  c = 0
  DO i=1,SIZE(AUXNAMES)
    IF( AUXNAMES(i)=="local_symmetry" ) c=i
  ENDDO
  !Otherwise, create a new column to store local symmetry parameter
  k = SIZE(AUXNAMES)
  IF( c==0 ) THEN
    k = k+1
    c = k
  ENDIF
  !
  IF( k>SIZE(AUXNAMES) ) THEN
    !Resize arrays
    ALLOCATE(newAUXNAMES(k))
    DO i=1,SIZE(AUXNAMES)
      newAUXNAMES(i) = AUXNAMES(i)
    ENDDO
    DEALLOCATE(AUXNAMES)
    ALLOCATE(AUXNAMES(SIZE(newAUXNAMES)))
    AUXNAMES(:) = newAUXNAMES(:)
    DEALLOCATE(newAUXNAMES)
    ALLOCATE(newAUX(SIZE(P,1),k))
    newAUX(:,:) = 0.d0
    DO i=1,SIZE(AUX,1)
      newAUX(i,:) = AUX(i,:)
    ENDDO
    DEALLOCATE(AUX)
    ALLOCATE(AUX(SIZE(newAUX,1),SIZE(newAUX,2)))
    AUX(:,:) = newAUX(:,:)
    DEALLOCATE(newAUX)
  ENDIF
ELSE
  !No auxiliary prop. existed: allocate new arrays
  IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
  IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
  c = 1
  ALLOCATE(AUXNAMES(1))
  ALLOCATE( AUX(SIZE(P,1),1) )
  AUX(:,:) = 0.d0
ENDIF
!
AUXNAMES(c) = "local_symmetry"
!
!
!
200 CONTINUE
! Compute the central symmetry parameter for each atom
!
!Construct neighbor list
CALL ATOMSK_MSG(11,(/""/),(/0.d0/))
CALL VERLET_LIST(H,P,6.d0,NeighList)
CALL ATOMSK_MSG(15,(/""/),(/0.d0/))
!PRINT*, "SIZE NeighList = ", SIZE(NeighList,1), SIZE(NeighList,2)
IF( .NOT.ALLOCATED(NeighList) .OR. SIZE(NeighList)<1 ) THEN
  nerr=nerr+1
  CALL ATOMSK_MSG(1815,(/""/),(/0.d0/))
  GOTO 1000
ENDIF
!
IF( verbosity==4 ) THEN
  !Some debug messages
  IF( ALLOCATED(NeighList) .AND. SIZE(NeighList,1)>0 ) THEN
    WRITE(msg,*) "Sample of neighbor list:"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO i=1,MIN(10,SIZE(NeighList,1))
      WRITE(msg,'(i5,a3,20i5)') i, " | ", NeighList(i,1:MIN(SIZE(NeighList,2),16))
      msg = TRIM(ADJUSTL(msg))//' (...)'
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
    IF( i>=10 ) THEN
      WRITE(msg,*) '      (...discontinued...)'
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDIF
  ENDIF
ENDIF
!
!
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j,k,ktemp,ipairs,Nneigh1,Nneigh,Neighsp,sum_dj,PosList,newindex,sp1,sp2) &
!$OMP& PRIVATE(pairs,pairs_distances,three_angles,dmin,distance,msg,temp)
DO i=1,SIZE(P,1)
  !i is the index of the central atom
  !Initialize variables and arrays for atom #i
  ipairs = 0
  ktemp = 0
  Nspecies = 0
  Nneigh(:) = 0
  Neighsp(:) = 0
  sum_dj = 0.d0
  IF(ALLOCATED(pairs)) DEALLOCATE(pairs)
  IF(ALLOCATED(pairs_distances)) DEALLOCATE(pairs_distances)
  IF(ALLOCATED(three_angles)) DEALLOCATE(three_angles)
  !
  !Get the positions of neighbors of atom #i
  CALL NEIGHBOR_POS(H,P,P(i,1:3),NeighList(i,:),ALLOCATED(NeighList),6.d0,PosList)
  !
  IF( ALLOCATED(PosList) .AND. SIZE(PosList,1)>0 ) THEN
    !
    WRITE(msg,*) "Atom # ", i, " : size of PosList = ", SIZE(PosList,1)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !
    !Sort neighbors by increasing distance
    CALL BUBBLESORT(PosList(:,:),4,'up  ',newindex)
    !
    ! Determine how many different atom species are present in PosList, and their number
    Nspecies = 0
    DO j=1,SIZE(PosList,1)
      
    ENDDO
    WRITE(msg,*) "Nspecies:", Nspecies
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !
    !Count number of 1st, 2nd, and 3rd neighbours
    IF( Nspecies>1 ) THEN
      !If Nspecies>1 ("complex" material), keep only neighbours of same species as 1st neighbour
      k = NINT(PosList(1,5)) !index of first neighbour
      sp1 = NINT(P(k,4))     !species of first neighbour
      Neighsp(1) = P(k,4)
      sp2 = sp1
      j = 1
      DO WHILE( sp2==sp1 .AND. j<=SIZE(PosList,1) )
        !Compute distance between central atom #i and its third neighbour in the list
        k = MIN(1,SIZE(PosList,1))
        distance = PosList(k,4)  !SUM(PosList(1:k,4)) / DBLE(k)
        IF( PosList(j,4) < 1.3d0*distance ) THEN
          !Atom #j's distance to central atom #i is comparable to that of 3rd neighbour
          !=> count atom #j as a first neighbour of atom #i
          Nneigh(1) = Nneigh(1)+1
        ENDIF
        !Save species of j-th neighbour, assume it is the species of all 1st neighbours
        j = j+1
        k = NINT(PosList(j,5))
        sp2 = NINT(P(k,4))  !species of neighbour #j
      ENDDO
      !
    ELSE
      !i.e. Nspecies=1, "simple" material
      DO j=1,MIN(Nmax,SIZE(PosList,1))
        !Compute distance between central atom #i and its third neighbour in the list
        k = MIN(3,SIZE(PosList,1))
        distance = PosList(k,4)  !SUM(PosList(1:k,4)) / DBLE(k)
        IF( PosList(j,4) < 1.1d0*distance ) THEN
          !Atom #j's distance to central atom #i is comparable to that of 3rd neighbour
          !=> count atom #j as a first neighbour of atom #i
          Nneigh(1) = Nneigh(1)+1
        ENDIF
        !Save species of 1st neighbour, assume it is the species of all 1st neighbours
        k = NINT(PosList(1,5))
        Neighsp(1) = P(k,4)
      ENDDO
    ENDIF
    !Now we know that there are Nneigh(1) first neighbours
    !Check that PosList contains other neighbours
    IF( Nneigh(1) < SIZE(PosList,1) ) THEN
      !Count 2nd neighbours
      DO j=Nneigh(1)+1,SIZE(PosList,1)
        !Compute distance between central atom #i and the Nneigh(1)+1 neighbour
        k = MIN(Nneigh(1)+1,SIZE(PosList,1))
        distance = PosList(k,4)
        IF( PosList(j,4) < 1.1d0*distance ) THEN
          !Atom #j is a second neighbour to atom #i
          Nneigh(2) = Nneigh(2)+1
        ENDIF
        !Save species of 2nd neighbour, assume it is the species of all 2nd neighbours
        k = NINT(PosList(Nneigh(1)+1,5))
        Neighsp(2) = P(k,4)
      ENDDO
    ENDIF
    !Now we know that there are Nneigh(1) first neighbours and Nneigh(2) second neighbours
    !Check that PosList still contains other neighbours
    IF( Nneigh(1)+Nneigh(2) < SIZE(PosList,1) ) THEN
      !Count 3rd neighbours
      DO j=Nneigh(1)+Nneigh(2)+1,SIZE(PosList,1)
        !Compute distance between central atom #i and the Nneigh(1)+1 neighbour
        k = MIN(Nneigh(1)+Nneigh(2)+1,SIZE(PosList,1))
        distance = PosList(k,4)
        IF( PosList(j,4) < 1.1d0*distance ) THEN
          !Atom #j is a third neighbour to atom #i
          Nneigh(3) = Nneigh(3)+1
        ENDIF
        !Save species of 3rd neighbour, assume it is the species of all 3rd neighbours
        k = NINT(PosList(Nneigh(1)+Nneigh(2)+1,5))
        Neighsp(3) = P(k,4)
      ENDDO
    ENDIF
    !
    !For the following, keep only first neighbors to compute local symmetry parameter
    !PosList(Nneigh(1)+1:,:) = 0.d0
    !
    !Make the PosList more compact by removing intermediate zeros
!     DO j=2,SIZE(PosList,1)
!       IF( NINT(PosList(j,5))==0 ) THEN
!         DO k=j,SIZE(PosList,1)-1
!           PosList(k,:) = PosList(k+1,:)
!         ENDDO
!       ENDIF
!     ENDDO
    !
    IF( verbosity==4 ) THEN
      !Some debug messages
      CALL ATOMSPECIES(P(i,4),species)
      WRITE(msg,*) i
      WRITE(temp,*) Nneigh(1)
      WRITE(msg,*) "Number of 1st neighbors of atom # "//TRIM(ADJUSTL(msg))//&
                & "("//TRIM(species)//"): "//TRIM(ADJUSTL(temp))
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      WRITE(msg,*) "   j |      x        y        z    |     d_ij    at.number"
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      DO j=1,MIN(20,Nneigh(1))
        WRITE(msg,'(i5,a3,3f9.3,a3,2f9.3)') j, " | ", PosList(j,1:3), " | ", PosList(j,4), P(NINT(PosList(j,5)),4)
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDDO
      IF( j>=20 ) THEN
        WRITE(msg,*) '      (...discontinued...)'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDIF
    ENDIF
    !
    IF( Nneigh(1) >= 2 ) THEN
      !Atom #i has many neighbors, we can work with that
      Nneigh1 = Nneigh(1)
      !
      !We must make pairs of atoms => Nneigh1 has to be an even number
      IF( MOD(Nneigh1,2) .NE. 0 ) THEN
        Nneigh1 = Nneigh1-1
      ENDIF
      !
      !We will need a table containing pairs of atoms
      ALLOCATE( pairs(Nneigh1/2,2) )
      pairs(:,:) = 0
      ALLOCATE( pairs_distances(Nneigh1/2) )
      pairs_distances(:) = 0.d0
      !
      !Find atoms j and k that are opposite
      ipairs=0
      sum_dj = 0.d0
      DO j=1,MIN(Nneigh1,Mdefault)
        !Sum the square of distances of all neighbours
        sum_dj = sum_dj + ( VECLENGTH(PosList(j,1:3)-P(i,1:3)) )**2
        IF( .NOT. ANY(pairs(:,:)==j) ) THEN
          !We have the first member of a new pair
          ipairs = ipairs+1
          pairs(ipairs,1) = j
          !
          !Among the remaining atoms, find the atom #k that minimizes D = |dj + dk|²
          ! *AND* that is of the same species as atom #j
          dmin = 1.d12
          ktemp = 0
          DO k=j+1,MIN(Nneigh1,Mdefault)
            !Atom #k is of same species as atom #j
            IF( .NOT. ANY(pairs(:,:)==k) ) THEN
              !Atom #k was not paired with another atom yet
              distance = ( VECLENGTH( PosList(j,1:3) + PosList(k,1:3) - 2.d0*P(i,1:3) ) )**2
              IF( distance < dmin ) THEN
                dmin = distance
                ktemp = k
              ENDIF
            ENDIF
          ENDDO
          !Now atom #ktemp is the best match to atom #j
          IF( ktemp>0 ) THEN
            pairs(ipairs,2) = ktemp
            pairs_distances(ipairs) = dmin
          ELSE
            !No suitable atom was found to pair atom #j => remove this entry
            pairs(ipairs,:) = 0
            ipairs = ipairs-1
          ENDIF
        ENDIF
      ENDDO
      !
      IF( verbosity==4 ) THEN
        WRITE(msg,*) i
        WRITE(msg,*) "Pairs of atoms around atom # "//TRIM(ADJUSTL(msg))//" :"
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        WRITE(msg,*) "   j    k      |dj+dk|² "
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        DO j=1,SIZE(pairs,1)
          WRITE(msg,'(2i5,f12.3)') pairs(j,1), pairs(j,2), pairs_distances(j)
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        ENDDO
      ENDIF
      !
      !Compute the central symmetry parameter c of atom #i
      AUX(i,c) = SUM(pairs_distances(:)) / (2.d0*sum_dj)
      !
      DEALLOCATE(pairs)
      DEALLOCATE(pairs_distances)
      !
      !Update scores of lattice types
      
      !
      !
      !Check if atom #i has 4 neighbours: maybe it has a tetrahedral environment
      !e.g. in silicon or diamond lattice, or in minerals containing SiO4 tetrahedra
      !Actually do this calculation if Nneigh is 3, 4 or 5
      IF( Nneigh(1)>=3 .AND. Nneigh(1)<=5 ) THEN
        !Compute angles between pairs of neighbours, using atom #i as central atom
        !In a perfect tetrahedron we expect 6 unique pairs, but here Nneigh may be 3, 4 or 5
        !therefore there are  Nneigh*(Nneigh-1)/2  unique pairs
        ipairs = Nneigh(1)*(Nneigh(1)-1)/2
        ALLOCATE( three_angles(ipairs) )
        three_angles(:) = 0.d0
        ipairs=0
        DO j=1,Nneigh(1)-1
          DO k=j+1,Nneigh(1)
            ipairs=ipairs+1
            !Compute angle between position vectors r_ij and r_ik
            dmin = ANGVEC( PosList(j,1:3)-P(i,1:3) , PosList(k,1:3)-P(i,1:3) )
            !Compute the squared difference of cosine, save it into three_angles(:)
            !NOTE: in a perfect tetrahedron the angle is a=acos(-1/3)=109.471°, hence -cos(a)=1/3
            three_angles(ipairs) = ( DCOS(dmin) + 1.d0/3.d0 )**2
          ENDDO
        ENDDO
        !Compute the sum of all angles and normalize
        dmin = SUM(three_angles(:)) / SIZE(three_angles)
        WRITE(msg,*) dmin
        WRITE(msg,*) "Normalized sum of squared dihedral angles: "//TRIM(ADJUSTL(msg))
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        !If this criterion is smaller than centro-symmetry criterion, use it instead
        IF( dmin<AUX(i,c) ) THEN
          AUX(i,c) = dmin
        ENDIF
        DEALLOCATE(three_angles)
      ENDIF
      !
      !
      !Check if atom #i has 3 neighbours of fewer: maybe it is a sp2 ("graphene") environment
      !Actually do this calculation if Nneigh is 3 or 4
      IF( Nneigh(1)<=4 ) THEN
        !Compute angles between pairs of neighbours, using atom #i as central atom
        !In perfect sp2 we expect 3 unique pairs, but here Nneigh may be 3 or 4
        !therefore there are  Nneigh*(Nneigh-1)/2  unique pairs
        ipairs = Nneigh(1)*(Nneigh(1)-1)/2
        ALLOCATE( three_angles(ipairs) )
        three_angles(:) = 0.d0
        ipairs=0
        DO j=1,Nneigh(1)-1
          DO k=j+1,Nneigh(1)
            ipairs=ipairs+1
            !Compute angle between position vectors r_ij and r_ik
            dmin = ANGVEC( PosList(j,1:3)-P(i,1:3) , PosList(k,1:3)-P(i,1:3) )
            !Compute the squared difference of cosine, save it into three_angles(:)
            !NOTE: in perfect sp2 the angle is 120°, hence cos(120°)=-1/2
            three_angles(ipairs) = ( DCOS(dmin) + 0.5d0 )**2
          ENDDO
        ENDDO
        !Compute the sum of all angles and normalize
        dmin = SUM(three_angles(:)) / SIZE(three_angles)
        WRITE(msg,*) dmin
        WRITE(msg,*) "Normalized sum of squared dihedral angles: "//TRIM(ADJUSTL(msg))
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        !If this criterion is smaller than centro-symmetry criterion, use it instead
        IF( dmin<AUX(i,c) ) THEN
          AUX(i,c) = dmin
        ENDIF
        DEALLOCATE(three_angles)
      ENDIF
      !
      !
    ELSEIF( Nneigh(1)==1 ) THEN
      !Atom #i has only one neighbor
      AUX(i,c) = 1.d0
    ELSE
      !Isolated atom with no neighbor: set parameter to zero
      AUX(i,c) = 0.d0
    ENDIF
    !
    WRITE(msg,*) i
    WRITE(msg,*) "c("//TRIM(ADJUSTL(msg))//") = ", AUX(i,c)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !
    IF( AUX(i,c) > 1.d3 .OR. AUX(i,c) < 0.d0 ) THEN
      !Something went wrong in the calculation!
      WRITE(*,*) "/!\ WARNING: auxiliary property out-of-bound for atom: ", i, AUX(i,c)
      AUX(i,c) = 0.d0
    ENDIF
    !
  ENDIF
  !
ENDDO
!$OMP END PARALLEL DO
!
!
!
400 CONTINUE
!Write atom positions and central symmetry parameters into output file(s)
IF( .NOT.ALLOCATED(outfileformats) ) THEN
  ALLOCATE(outfileformats(1))
  outfileformats(1) = "cfg"
ENDIF
CALL WRITE_AFF(prefix,outfileformats,H,P,S,comment,AUXNAMES,AUX)
!
!
GOTO 1000
!
!
!
850 CONTINUE
nerr = nerr+1
GOTO 1000
!
!
!
1000 CONTINUE
IF( ALLOCATED(pairs) ) DEALLOCATE(pairs)
IF( ALLOCATED(pairs_distances) ) DEALLOCATE(pairs_distances)
IF( ALLOCATED(three_angles) ) DEALLOCATE(three_angles)
!
!
END SUBROUTINE LOCAL_SYM
!
!
END MODULE mode_localsym
