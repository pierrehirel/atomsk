MODULE mode_cpprop
!
!**********************************************************************************
!*  MODE_CPPROP                                                                   *
!**********************************************************************************
!* This module reads two sets of atomic data, and copies auxiliary properties     *
!* from the first into the second.                                                *
!**********************************************************************************
!* (C) March 2023 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 20 March 2023                                    *
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
SUBROUTINE COPY_PROPERTIES(filefirst,filesecond,options_array,prefix,outfileformats)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: filefirst, filesecond, prefix
CHARACTER(LEN=2):: species
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE,INTENT(IN):: outfileformats !list of output file formats
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment1, comment2
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES1, AUXNAMES2 !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list (system 2 only)
INTEGER:: i, j, k, l, m, n
INTEGER:: Nnew, Now, Ntot  !number of "new" and overwritten aux.prop.
REAL(dp):: tempreal
REAL(dp),DIMENSION(3,3):: Huc    !Unit cell vectors (unknown here, set to zero)
REAL(dp),DIMENSION(3,3):: H1, H2 !Box vectors of systems 1 and 2
REAL(dp),DIMENSION(3,3):: ORIENT     !crystallographic orientation (not used here but necessary for some calls)
REAL(dp),DIMENSION(9,9):: C_tensor    !stiffness tensor (not used here but necessary for some calls)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX1, AUX2  !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX         !final auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P1, P2
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S1, S2


!Initialize variables and arrays
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
Nnew=0
Now=0
Ntot=0
C_tensor(:,:) = 0.d0   !stiffness tensor (not used, set to zero)
ORIENT(:,:) = 0.d0     !crystallographic orientation (not used, set to zero)
Huc(:,:) = 0.d0
!
!
CALL ATOMSK_MSG(4073,(/filefirst,filesecond/),(/0.d0/))
!
!
!
100 CONTINUE
!**********************************************************************************
!                            READ  INPUT  FILES
!**********************************************************************************
!Read atomic positions from filefirst and store them into P1(:,:)
CALL READ_AFF(filefirst,H1,P1,S1,comment1,AUXNAMES1,AUX1)
IF(nerr>0 .OR. .NOT.ALLOCATED(P1)) GOTO 1000
!Get rid of atoms and shells
IF (ALLOCATED(P1)) DEALLOCATE(P1)
IF (ALLOCATED(S1)) DEALLOCATE(S1)
!
!If there are no aux.prop. in first system, then nothing to copy => abort
IF( .NOT.ALLOCATED(AUX1) .OR. .NOT.ALLOCATED(AUXNAMES1) ) THEN
  nerr=nerr+1
  CALL ATOMSK_MSG(4834,(/filefirst/),(/0.d0/))
  GOTO 1000
ELSE
  Ntot = SIZE(AUXNAMES1)
ENDIF
!
!Read atomic positions from filesecond and store them into P2(:,:)
CALL READ_AFF(filesecond,H2,P2,S2,comment2,AUXNAMES2,AUX2)
IF(nerr>0 .OR. .NOT.ALLOCATED(P2)) GOTO 1000
!Apply options to system 2
CALL OPTIONS_AFF(options_array,Huc,H2,P2,S2,AUXNAMES1,AUX1,ORIENT,SELECT,C_tensor)
!
!
200 CONTINUE
!Check the size of arrays
IF( SIZE(P2,1) .NE. SIZE(AUX1,1) ) THEN
  !Systems 1 and 2 do not have the same number of atoms: display warning
  nwarn=nwarn+1
  CALL ATOMSK_MSG(4721,(/""/),(/DBLE(SIZE(AUX1,1)),DBLE(SIZE(P2,1))/))
ENDIF
!
!Check if there were already aux.prop. in system 2
IF( ALLOCATED(AUXNAMES2) ) THEN
  !There are aux.prop. in both initial systems
  !Determine how many aux.prop. there will be in final system
  !Assume that they are all different
  Nnew = SIZE(AUXNAMES1) + SIZE(AUXNAMES2)
  !Count duplicates: aux.prop. that exist in AUX1, but not in AUX2
  DO m=1,SIZE(AUXNAMES1)
    DO n=1,SIZE(AUXNAMES2)
      IF( AUXNAMES1(m)==AUXNAMES2(n) ) THEN
        Nnew = Nnew - 1
      ENDIF
    ENDDO
  ENDDO
  !
  !Allocate arrays for final data
  ALLOCATE(AUXNAMES(Nnew))
  AUXNAMES(:) = ""
  ALLOCATE(AUX(SIZE(P2,1),Nnew))
  AUX(:,:) = 0.d0
  !
  !Copy data
  !First, copy all names and aux.prop. from system2
  DO j=1,SIZE(AUXNAMES2)
    AUXNAMES(j) = AUXNAMES2(j)
    DO i=1,SIZE(AUX2,1)
      IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
        AUX(:,j) = AUX2(:,j)
      ENDIF
    ENDDO
  ENDDO
  !
  !Second, copy aux.prop. from system 1 (overwrite previous data if necessary)
  Nnew = SIZE(AUXNAMES2)
  DO m=1,SIZE(AUXNAMES1)
    l=0
    !Colmun #m of AUXNAMES1: determine in which column of AUX to copy it
    DO n=1,SIZE(AUXNAMES2)
      IF( AUXNAMES1(m)==AUXNAMES2(n) ) THEN
        !Aux.prop. exists in both system 1 and system 2
        !Use same column in AUX as in AUX2 (data will be overwritten)
        j = n
        l = 1
        Now = Now+1
      ENDIF
    ENDDO
    IF( l==0 ) THEN
      !This aux.prop. exists in system 1 but not in system 2
      !This will be a new column
      Nnew = Nnew+1
      j = Nnew
    ENDIF
    !Copy data
    !NB: systems 1 and 2 may have different numbers of atoms, hence P2 and AUX1
    !    may not have the same length, use the smallest of the two
    DO i=1,MIN(SIZE(P2,1),SIZE(AUX1,1))
      IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
        AUXNAMES(j) = AUXNAMES1(m)
        AUX(i,j) = AUX1(i,m)
      ENDIF
    ENDDO
  ENDDO
  !
  Ntot = SIZE(AUXNAMES)
  !
ELSE
  !AUXNAMES2 was not allocated (no aux.prop. in system 2)
  !Compare numbers of atoms
  IF( SIZE(P2,1) .NE. SIZE(AUX1,1) ) THEN
    !Systems 1 and 2 do not contain the same number of atoms
    ALLOCATE(AUXNAMES(SIZE(AUXNAMES1)))
    AUXNAMES(:) = AUXNAMES1(:)
    j = MIN(SIZE(P2,1),SIZE(AUX1,1))
    ALLOCATE(AUX(j,SIZE(AUXNAMES1)))
    DO i=1,j
      AUX(i,:) = AUX1(i,:)
    ENDDO
  ENDIF
ENDIF
!
!
!
400 CONTINUE
!Display number of aux.prop. that were copied, that are new, that were overwritten, and total
j = SIZE(AUXNAMES1)
CALL ATOMSK_MSG(4074,(/""/),(/DBLE(j),DBLE(j-Now),DBLE(Now),DBLE(Ntot)/))
!
!Final system is made of atom positions of system 2 + aux.prop. that were copied
!Write final system into file(s)
IF( ALLOCATED(AUXNAMES) .AND. ALLOCATED(AUX) ) THEN
  CALL WRITE_AFF(prefix,outfileformats,H2,P2,S2,comment2,AUXNAMES,AUX)
ELSE
  CALL WRITE_AFF(prefix,outfileformats,H2,P2,S2,comment2,AUXNAMES1,AUX1)
ENDIF
!
!
!
1000 CONTINUE
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
!
!
END SUBROUTINE COPY_PROPERTIES
!
!
END MODULE mode_cpprop
