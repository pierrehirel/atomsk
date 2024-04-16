MODULE selection
!
!**********************************************************************************
!*  SELECTION                                                                     *
!**********************************************************************************
!* This module contains general subroutines dealing                               *
!* with files, manipulation of strings, arrays and so on.                         *
!**********************************************************************************
!* (C) April 2024 - Pierre Hirel (unless specified otherwise, cf each subroutine) *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 16 April 2024                                    *
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
!* IS_SELECTED         determines if an atom is selected or not                   *
!* LIST_SELECTED_ATOMS makes a list of selected atoms (chem.species/number)       *
!**********************************************************************************
!
!
USE comv
!
IMPLICIT NONE
!
!
CONTAINS
!
!
!
!********************************************************
!  IS_SELECTED
!  This function determines if an atom is selected or not.
!  If array SELECT is unallocated, atom is considered
!  selected (function returns .TRUE.).
!  If it is allocated, function returns the value
!  of SELECT(i).
!********************************************************
LOGICAL FUNCTION IS_SELECTED(SELECT,i)
!
IMPLICIT NONE
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT
INTEGER,INTENT(IN):: i
!
IS_SELECTED = .TRUE.
!
IF( ALLOCATED(SELECT) ) THEN
  IF( i>0 .AND. i<=SIZE(SELECT) ) THEN
    IS_SELECTED = SELECT(i)
  ELSE
    IS_SELECTED = .FALSE.
  ENDIF
ELSE
  IS_SELECTED = .TRUE.
ENDIF
!
END FUNCTION IS_SELECTED
!
!
!********************************************************
!  LIST_SELECTED_ATOMS
!  This subroutine generates a table containing
!  chemical species and number of selected atoms.
!********************************************************
!
SUBROUTINE LIST_SELECTED_ATOMS(P,SELECT,selectedlist)
!
IMPLICIT NONE
INTEGER:: i, j
INTEGER:: Nspecies
INTEGER,DIMENSION(20,2):: templist !assume maximum 20 different chemical species
INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: selectedlist
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT
REAL(dp),DIMENSION(:,:),INTENT(IN):: P
!
Nspecies=0
templist(:,:) = 0
IF(ALLOCATED(selectedlist)) DEALLOCATE(selectedlist)
!
DO i=1,SIZE(P,1)
  IF( IS_SELECTED(SELECT,i) ) THEN
    DO j=1,SIZE(templist,1)
      IF( templist(j,1) == NINT(P(i,4)) ) THEN
        templist(j,2) = templist(j,2) + 1
        EXIT
      ELSEIF( templist(j,1) == 0 ) THEN
        templist(j,1) = NINT(P(i,4))
        templist(j,2) = templist(j,2) + 1
        Nspecies = Nspecies + 1
        EXIT
      ENDIF
    ENDDO
  ENDIF
ENDDO
!
ALLOCATE(selectedlist(Nspecies,2))
DO j=1,SIZE(selectedlist,1)
  selectedlist(j,:) = templist(j,:)
ENDDO
!
END SUBROUTINE LIST_SELECTED_ATOMS
!
!
!
END MODULE selection
