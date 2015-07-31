MODULE mode_centrosym
!
!**********************************************************************************
!*  MODE_CENTROSYM                                                                *
!**********************************************************************************
!* This module computes a central symmetry parameter for each atom in the system, *
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
!**********************************************************************************
!* (C) April 2015 - Pierre Hirel                                                  *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 29 July 2015                                     *
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
SUBROUTINE CENTRO_SYM(inputfile,options_array,prefix,outfileformats)
!
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=*),INTENT(IN):: prefix
!
CHARACTER(LEN=2):: species
CHARACTER(LEN=27):: pbar
CHARACTER(LEN=128):: msg
CHARACTER(LEN=4096):: outputfile
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of output file formats
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES, newAUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: c             !column in AUX that will contain the central symmetry parameters c(i)
INTEGER:: i, ipairs, j, k, ktemp
INTEGER:: Mdefault      !most common number of neighbors
INTEGER:: Nneigh        !number of neighbors of an atom
INTEGER:: Nspecies      !number of species in the system
INTEGER:: percent
INTEGER,DIMENSION(:,:),ALLOCATABLE:: NeighList !list of index of neighbors
INTEGER,DIMENSION(:,:),ALLOCATABLE:: pairs  !indexes of pairs of atoms
REAL(dp):: distance, dmin
REAL(dp):: snumber      !atomic number of an atom
REAL(dp):: progress     !progress of the computation
REAL(dp):: d_1N         !distance to first neighbor
REAL(dp):: sum_dj       !sum of |d_j|²
REAL(dp),DIMENSION(3,3):: H       !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystallographic orientation of the system (mode create)
REAL(dp),DIMENSION(:),ALLOCATABLE:: pairs_distances  !distances
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries    !array containing result
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX, newAUX !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, S     !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: PosList
REAL(dp),DIMENSION(:),POINTER:: PropPoint !pointer to the property whose density must be computed
!
!
!Initialize variables
Nspecies = 0
!
!
!CALL ATOMSK_MSG(4066,(/property/),(/DBLE(den_type),Sigma/))
!
!
100 CONTINUE
! Read the input file
CALL READ_AFF(inputfile,H,P,S,comment,AUXNAMES,AUX)
!
! Apply options if any
CALL OPTIONS_AFF(options_array,H,P,S,AUXNAMES,AUX,ORIENT,SELECT)
IF(nerr>0) GOTO 1000
!
!If auxiliary properties already exist, add a column to save the central symmetry parameter
IF( ALLOCATED(AUX) .AND. SIZE(AUX,2)==SIZE(AUXNAMES) ) THEN
  c = SIZE(AUXNAMES) + 1
  ALLOCATE(newAUXNAMES(SIZE(AUXNAMES)+1))
  DO i=1,SIZE(AUXNAMES)
    newAUXNAMES(i) = AUXNAMES(i)
  ENDDO
  DEALLOCATE(AUXNAMES)
  ALLOCATE(AUXNAMES(SIZE(newAUXNAMES)))
  AUXNAMES(:) = newAUXNAMES(:)
  DEALLOCATE(newAUXNAMES)
  ALLOCATE(newAUX(SIZE(P,1),SIZE(AUX,2)+1))
  newAUX(:,:) = 0.d0
  DO i=1,SIZE(AUX,1)
    newAUX(i,:) = AUX(i,:)
  ENDDO
  DEALLOCATE(AUX)
  ALLOCATE(AUX(SIZE(newAUX,1),SIZE(newAUX,2)))
  AUX(:,:) = newAUX(:,:)
  DEALLOCATE(newAUX)
ELSE
  IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
  IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
  c = 1
  ALLOCATE(AUXNAMES(1))
  ALLOCATE( AUX(SIZE(P,1),1) )
  AUX(:,:) = 0.d0
ENDIF
!
AUXNAMES(c) = "central symm."
!
!
!
200 CONTINUE
! Determine how many different atom species are present, and their number
CALL FIND_NSP(P(:,4),aentries)
WRITE(msg,*) "Number of different species:", SIZE(aentries,1)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
IF( SIZE(aentries,1)>1 ) THEN
  !Sort them by increasing number of atoms
  CALL BUBBLESORT(aentries(:,:),2,'up  ')
ENDIF
!
! Determine what type of system it is: unary, binary, ternary, other
IF( SIZE(aentries,1)==1 ) THEN
  !No ambiguity, it is a unary system
  Nspecies = 1
ELSEIF( SIZE(aentries,1)==2 ) THEN
  !It may be a unary system with just a few solute atoms
  !Check if a species is a minority
  IF( aentries(1,2) > 2.9d0*aentries(2,2) .OR. aentries(2,2) > 2.9d0*aentries(1,2) ) THEN
    !One of the species is a minority => consider that it is a unary system
    Nspecies = 1
  ELSE
    !Consider that it is a binary system
    Nspecies = 2
  ENDIF
ELSEIF( SIZE(aentries,1)==3 ) THEN
  !It may be a unary system with just a few solute atoms
  !Check if a species is a minority
  IF( aentries(1,2) > 2.9d0*aentries(2,2) .OR. aentries(2,2) > 2.9d0*aentries(1,2) .OR.  &
    & aentries(1,2) > 2.9d0*aentries(3,2) .OR. aentries(3,2) > 2.9d0*aentries(1,2) .OR.  &
    & aentries(2,2) > 2.9d0*aentries(3,2) .OR. aentries(3,2) > 2.9d0*aentries(2,2)      ) THEN
    !One or two of the species are a minority => find which one
    Nspecies = 3
  ELSE
    !No minority species => it is a ternary material
    Nspecies = 3
  ENDIF
ELSE
  Nspecies = 0
  PRINT*, "X!X ERROR Nspecies is zero"
  GOTO 1000
ENDIF
WRITE(msg,*) "Nspecies:", Nspecies
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!
300 CONTINUE
! Compute the central symmetry parameter for each atom
!
!
!Construct neighbor list
CALL ATOMSK_MSG(11,(/""/),(/0.d0/))
CALL NEIGHBOR_LIST(H,P,8.d0,NeighList)
PRINT*, "SIZE NeighList = ", SIZE(NeighList,1), SIZE(NeighList,2)
!
IF( verbosity==4 ) THEN
  !Some debug messages
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
!
!Determine most common number of neighbors
!   Nneigh=0
!   DO i=1,SIZE(NeighList,1)
!     !Get the positions of neighbors of atom #i
!     CALL NEIGHBOR_POS(H,P,P(i,1:3),NeighList(i,:),6.d0,PosList)
!     !
!     !Sort neighbors by increasing distance
!     CALL BUBBLESORT(PosList,4,'up  ')
!     !
!     !Distance to the closest neighbor
!     d_1N = VECLENGTH( PosList(1,1:3) - P(i,1:3) )
!     !
!     !Count number of first neighbors of atom #i
!     DO j=1,SIZE(NeighList)
!       distance = VECLENGTH( PosList(i,1:3) - P(i,1:3) )
!       IF( distance < d_1N*1.2d0 ) THEN
!         Nneigh=Nneigh+1
!       ENDIF
!     ENDDO
!   ENDDO
!
!Average number of neighbors
!Mdefault = Nneigh / SIZE(P,1)
Mdefault = 12
!We must make pairs of atoms => Mdefault has to be an even number
IF( MOD(Mdefault,2) .NE. 0 ) THEN
  Mdefault = Mdefault-1
ENDIF
WRITE(msg,*) "Most common number of neighbors: ", Mdefault
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!We will need a table containing pairs of atoms
ALLOCATE( pairs(Mdefault/2,2) )
ALLOCATE( pairs_distances(Mdefault/2) )
!
DO i=1,SIZE(P,1)
  !i is the index of the central atom
  !Initialize variables and arrays for atom #i
  ipairs = 0
  Nneigh = 0
  pairs(:,:) = 0
  pairs_distances(:) = 0.d0
  sum_dj = 0.d0
  !
  !Get the positions of neighbors of atom #i
  CALL NEIGHBOR_POS(H,P,P(i,1:3),NeighList(i,:),8.d0,PosList)
  !
  !Sort neighbors by increasing distance
  CALL BUBBLESORT(PosList(:,:),4,'up  ')
  !
  !Select which neighbors to keep: that depends on the type of compound
  SELECT CASE(Nspecies)
  CASE(1)
    !Unary compound => follow the recipe of Kelchner et al.
    !keep only neighbors j that are approx. at the same distance as the four 1st neighbors
    Nneigh = 0
    DO j=1,SIZE(PosList,1)
      IF( PosList(j,4) < 1.4d0*SUM(PosList(1:4,4))/4.d0 ) THEN
        Nneigh = Nneigh+1
        sum_dj = sum_dj + ( VECLENGTH(PosList(j,1:3)-P(i,1:3)) )**2
      ELSE
        PosList(j,:) = 0.d0
      ENDIF
    ENDDO
    !
  CASE DEFAULT
    !Complex compound
    Nneigh = 0
    IF( aentries(1,2)/aentries(2,2) <= 1.d0 .AND. aentries(SIZE(aentries,1),2) > 2.d0*aentries(1,2) ) THEN
      !It is a material of formula ABCx with x>2 (for instance perovskite ABO3)
      !=> Keep only neighbors of different species as the central atom
      !   *and* that are approx. at the same distance as the 1st neighbor
      DO j=1,SIZE(PosList,1)
        !NINT(PosList(:,5)) is the actual index of atom
        IF( NINT(P(i,4))==NINT(aentries(3,1)) ) THEN
          !Central atom i is of type C => keep only neighbors that are not C
          IF( NINT(P(NINT(PosList(j,5)),4)).NE.NINT(aentries(3,1)) .AND.     &
            & PosList(j,4) < 1.25d0*SUM(PosList(1:4,4))/4.d0             ) THEN
            Nneigh = Nneigh+1
            sum_dj = sum_dj + ( VECLENGTH(PosList(j,1:3)-P(i,1:3)) )**2
          ELSE
            PosList(j,:) = 0.d0
          ENDIF
        ELSE
          !Central atom i is A or B => keep only neighbors of type C
          IF( NINT(P(NINT(PosList(j,5)),4))==NINT(aentries(3,1)) .AND.      &
            & PosList(j,4) < 1.4d0*SUM(PosList(1:4,4))/4.d0             ) THEN
            Nneigh = Nneigh+1
            sum_dj = sum_dj + ( VECLENGTH(PosList(j,1:3)-P(i,1:3)) )**2
          ELSE
            PosList(j,:) = 0.d0
          ENDIF
        ENDIF
      ENDDO
    ELSE
      !It is another complex material ABCx with x<2
      DO j=1,SIZE(PosList,1)
        !NINT(PosList(:,5)) is the actual index of atom
        IF( NINT(P(NINT(PosList(j,5)),4)).NE.NINT(P(i,4)) &
          & .AND. PosList(j,4) < 1.4d0*SUM(PosList(1:4,4))/4.d0 ) THEN
          Nneigh = Nneigh+1
          sum_dj = sum_dj + ( VECLENGTH(PosList(j,1:3)-P(i,1:3)) )**2
        ELSE
          PosList(j,:) = 0.d0
        ENDIF
      ENDDO
    ENDIF
    !Make the PosList more compact by removing zeros
    DO j=2,SIZE(PosList,1)
      IF( NINT(PosList(j,5))==0 ) THEN
        DO k=j,SIZE(PosList,1)-1
          PosList(k,:) = PosList(k+1,:)
        ENDDO
      ENDIF
    ENDDO
    !
  !CASE DEFAULT
    !Even more complex compound!!
  END SELECT
  !
  !We must make pairs of atoms => Nneigh has to be an even number
  IF( MOD(Nneigh,2) .NE. 0 ) THEN
    Nneigh = Nneigh-1
  ENDIF
  IF( Nneigh.NE.Mdefault ) THEN
    PRINT*, "Atom # ", i, " has ", Nneigh, "neighbors:" !, (NeighList(i,j),j=1,Nneigh)
  ENDIF
  !
  IF( verbosity==4 ) THEN
    !Some debug messages
    WRITE(msg,*) i
    WRITE(pbar,*) Nneigh
    WRITE(msg,*) "Number of neighbors of atom # "//TRIM(ADJUSTL(msg))//": "//TRIM(ADJUSTL(pbar))
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO j=1,20 !MAX(12,Nneigh)
      WRITE(msg,'(i5,a3,3f9.3,a3,2f9.3)') j, " | ", PosList(j,1:3), " | ", PosList(j,4), P(NINT(PosList(j,5)),4)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
    IF( j>=20 ) THEN
      WRITE(msg,*) '      (...discontinued...)'
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDIF
  ENDIF
  !
  IF( Nneigh >= 2 ) THEN
    !Atom #i has many neighbors, we can work with that
    !Find atoms j and k that are opposite
    DO j=1,MIN(Nneigh,Mdefault)
      IF( .NOT. ANY(pairs(:,:)==j) ) THEN
        !We have the first member of a new pair
        ipairs = ipairs+1
        pairs(ipairs,1) = j
        !
        !Among the remaining atoms, find the one that minimizes D = |dj + dk|²
        Dmin = 1.d12
        DO k=j+1,MIN(Nneigh,Mdefault)
          IF( .NOT. ANY(pairs(:,:)==k) ) THEN
            !Atom #k was not paired with another atom yet
            distance = ( VECLENGTH( PosList(j,1:3) + PosList(k,1:3) - 2.d0*P(i,1:3) ) )**2
            IF( distance < Dmin ) THEN
              Dmin = distance
              ktemp = k
            ENDIF
          ENDIF
        ENDDO
        !Now atom #ktemp is the best match to atom #j
        pairs(ipairs,2) = ktemp
        pairs_distances(ipairs) = Dmin
      ENDIF
    ENDDO
    !
    IF( verbosity==4 ) THEN
      WRITE(msg,*) i
      WRITE(msg,*) "Pairs of atoms around atom # "//TRIM(ADJUSTL(msg))//" :"
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      DO j=1,SIZE(pairs,1)
        WRITE(msg,'(2i5,f9.3)') pairs(j,1), pairs(j,2), pairs_distances(j)
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDDO
    ENDIF
    !
    !Compute the central symmetry parameter c of atom #i
    AUX(i,c) = SUM(pairs_distances(:)) / (2.d0*sum_dj)
    !
  ELSEIF( Nneigh==1 ) THEN
    !Atom #i has only one neighbor
    AUX(i,c) = 1.d0
  ELSE
    !Isolated atom with no neighbor
    AUX(i,c) = 0.d0
  ENDIF
  !
  WRITE(msg,*) i
  WRITE(msg,*) "c("//TRIM(ADJUSTL(msg))//") = ", AUX(i,c)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  IF( AUX(i,c) > 1.d0 .OR. AUX(i,c)<0.d0 ) THEN
    !Something went wrong in the calculation!
    WRITE(*,*) "/!\ WARNING: auxiliary property is out-of-bound: ", AUX(i,c)
    AUX(i,c) = 0.d0
  ENDIF
ENDDO
!
!
!
400 CONTINUE
!Write atom positions and central symmetry parameters into output file(s)
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
!
!
END SUBROUTINE CENTRO_SYM
!
!
END MODULE mode_centrosym
