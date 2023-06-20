MODULE mode_matchid
!
!**********************************************************************************
!*  MODE_MATCHID                                                                  *
!**********************************************************************************
!* This module reads two arrays P1 and P2 containing atom positions,              *
!* and        *
!**********************************************************************************
!* (C) June 2023 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 20 June 2023                                     *
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
USE sorting
!Module for reading input files
USE readin
!Module for applying options
USE options
!Output modules
USE writeout
!
!
CONTAINS
!
!
SUBROUTINE MATCHID_XYZ(filefirst,filesecond,options_array,outputfile,outfileformats)
!
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: filefirst, filesecond
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: outputfile, msg
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES, AUXNAMES1, AUXNAMES2 !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment1, comment2
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: i, j
INTEGER:: NP1    !number of atoms in filefirst
INTEGER:: status
INTEGER,DIMENSION(:),ALLOCATABLE:: paired  !index of atom with which atom #i is paired
REAL(dp),DIMENSION(3,3):: Huc        !Box vectors of unit cell (unknown, set to 0 here)
REAL(dp),DIMENSION(3,3):: H1, H2     !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT     !crystallographic orientation of the system
REAL(dp),DIMENSION(9,9):: C_tensor   !elastic tensor
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX1, AUX2  !auxiliary properties of system 1 and 2
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P1, P2, newP2 !atom positions of system 1 and 2
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S1, S2, newS2 !shell positions of system 1 and 2
!
!
!Initialize variables
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
Huc(:,:) = 0.d0
 C_tensor(:,:) = 0.d0
ORIENT(:,:) = 0.d0
!
!
CALL ATOMSK_MSG(4076,(/''/),(/0.d0/))
!
!
100 CONTINUE
!Read atomic positions from filefirst and store them into P1(:,:)
CALL READ_AFF(filefirst,H1,P1,S1,comment1,AUXNAMES1,AUX1)
!Save N.atoms from system 1
NP1 = SIZE(P1,1)
!
!Read atomic positions from filesecond and store them into P2(:,:)
CALL READ_AFF(filesecond,H2,P2,S2,comment2,AUXNAMES2,AUX2)
!
!
!
300 CONTINUE
CALL ATOMSK_MSG(4077,(/''/),(/0.d0/))
!
!For each atom in P2, find its match in P1
CALL FIND_MATCHING_ID(P2,P1,paired,j)
!
IF( j==0 ) THEN
  !No equivalent pairs of atoms could be found => abort
  CALL ATOMSK_MSG(4835,(/''/),(/0.d0/))
  nerr = nerr+1
  GOTO 1000
ENDIF
!
CALL ATOMSK_MSG(4078,(/''/),(/DBLE(j),DBLE(SIZE(P2,1)-j)/))
!
IF( verbosity>=4 ) THEN
  OPEN(UNIT=31,FILE="atomsk_matchid_pairs.dat",FORM="FORMATTED",STATUS="UNKNOWN")
  WRITE(31,*) "# Atom indices and their replacements"
  DO i=1,SIZE(paired)
    WRITE(31,'(i6,2X,i6)') i, paired(i)
  ENDDO
  CLOSE(31)
ENDIF
!
!We won't need info about the first system any more
DEALLOCATE(P1)
IF(ALLOCATED(S1)) DEALLOCATE(S1)
IF(ALLOCATED(AUX1)) DEALLOCATE(AUX1)
IF(ALLOCATED(AUXNAMES1)) DEALLOCATE(AUXNAMES1)
IF(ALLOCATED(comment1)) DEALLOCATE(comment1)
!
!If P2 is larger than P1, it means it contains "new" atoms
IF( SIZE(P2,1) > NP1 .OR. ANY(paired(:)==0) ) THEN
  CALL ATOMSK_MSG(4079,(/""/),(/0.d0/))
  IF( ALLOCATED(AUX2) ) THEN
    !Aux.prop. already existed: copy their names
    j = SIZE(AUX2,2)+1
    ALLOCATE( AUXNAMES(j) )
    DO i=1,SIZE(AUXNAMES2)
      AUXNAMES(i) = AUXNAMES2(i)
    ENDDO
  ELSE
    !No aux.prop. before: "new_atoms" will be the only aux.prop.
    j = 1
    ALLOCATE( AUXNAMES(j) )
  ENDIF
  IF( SIZE(P2,1) > NP1 ) THEN
    !Atoms without an equivalent are labeled as "new_atom"
    AUXNAMES(j) = "new_atoms"
  ELSE
    !Atoms without an equivalent are labeled as "unpaired_atom"
    AUXNAMES(j) = "unpaired_atoms"
  ENDIF
  IF(ALLOCATED(AUXNAMES2)) DEALLOCATE(AUXNAMES2)
  ALLOCATE(AUXNAMES2(SIZE(AUXNAMES)))
  AUXNAMES2(:) = AUXNAMES(:)
  DEALLOCATE(AUXNAMES)
  !Resize array AUX2 to store the values of the new aux.prop. "new_atoms"
  CALL RESIZE_DBLEARRAY2(AUX2,SIZE(P2,1),j,status)
  IF(status>0) THEN
    CALL ATOMSK_MSG(818,(/"AUX2"/),(/0.d0/))
    nerr = nerr+1
    GOTO 1000
  ENDIF
  !Mark all unpaired atoms as "new_atoms"
  DO i=1,SIZE(paired)
    IF( paired(i)<=0 ) THEN
      AUX2(i,j) = 1.d0
    ENDIF
  ENDDO
ENDIF
!
!Sort atoms (and shells and/or aux.prop) in second system according to index list
CALL ATOMSK_MSG(4080,(/""/),(/0.d0/))
CALL IDSORT(paired,P2)
IF(ALLOCATED(S2)) CALL IDSORT(paired,S2)
IF(ALLOCATED(AUX2)) CALL IDSORT(paired,AUX2)
IF(ALLOCATED(SELECT)) CALL IDSORT_SELECT(paired,SELECT)
CALL ATOMSK_MSG(15,(/""/),(/0.d0/))
!
!
!
400 CONTINUE
!Apply options to final system
CALL OPTIONS_AFF(options_array,Huc,H2,P2,S2,AUXNAMES2,AUX2,ORIENT,SELECT,C_tensor)
!
!Write files
CALL WRITE_AFF(outputfile,outfileformats,H2,P2,S2,comment2,AUXNAMES2,AUX2)
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE MATCHID_XYZ
!
!
END MODULE mode_matchid
