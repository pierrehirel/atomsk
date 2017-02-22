MODULE sort
!
!**********************************************************************************
!*  SORT                                                                          *
!**********************************************************************************
!* This module sorts an array of type (species x y z) by increasing               *
!* or decreasing values of one of the columns.                                    *
!**********************************************************************************
!* (C) May 2010 - Pierre Hirel                                                    *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 21 Feb. 2017                                     *
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
USE constants
USE messages
USE files
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE SORT_XYZ(P,S,AUXNAMES,AUX,SELECT,sortcol,sortorder)
!
!
IMPLICIT NONE
CHARACTER(LEN=4):: sortorder  !up or down
CHARACTER(LEN=16):: sortcol    !property to be sorted: x, y, z or s, or any name
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: i, j
INTEGER:: sortnum   !index of column to sort
INTEGER,DIMENSION(:),ALLOCATABLE:: newindex  !list of index after sorting
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX    !auxiliary properties of atoms
REAL(dp),DIMENSION(:,:),ALLOCATABLE,TARGET:: P, S  !positions of atoms, shells
!
!Initialize variables
sortnum = 0
IF( ALLOCATED(newindex) ) DEALLOCATE(newindex) !will be allocated by QUICKSORT or PACKSORT routines
!
WRITE(msg,*) 'Entering SORT_XYZ: ', sortcol, sortorder
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
CALL ATOMSK_MSG(2087,(/sortcol//'   ', sortorder/),(/0.d0/))
!
!
!
100 CONTINUE
!Find out which property must be sorted
SELECT CASE (sortcol)
!
CASE('species','s','S','x','X','y','Y','z','Z')
  !The array P must be sorted
  !Determine the column to be sorted
  IF( sortcol=="x" .OR. sortcol=="X" ) THEN
    sortnum = 1
  ELSEIF( sortcol=="y" .OR. sortcol=="Y" ) THEN
    sortnum = 2
  ELSEIF( sortcol=="z" .OR. sortcol=="Z" ) THEN
    sortnum = 3
  ELSE
    sortnum = 4
  ENDIF
  WRITE(msg,*) 'Array to sort: P ; sortnum = ', sortnum
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  SELECT CASE(sortorder)
  CASE('up','down')
    !Sort P by increasing or decreasing values
    CALL QUICKSORT(P,sortnum,sortorder,newindex)
  CASE('pack')
    !Pack identical values of P together
    CALL PACKSORT(P,sortnum,newindex)
  CASE DEFAULT
  END SELECT
  !
  !newindex(:) now contains the new list of index (after sorting)
  !=> use it to sort other arrays
  !If shells are defined, sort them accordingly
  IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(newindex) ) THEN
    CALL IDSORT(newindex,S)
  ENDIF
  !
  !If auxiliary properties are defined, sort them accordingly
  IF( ALLOCATED(AUX) .AND. SIZE(AUX,1)==SIZE(newindex) ) THEN
    CALL IDSORT(newindex,AUX)
  ENDIF
  !
CASE DEFAULT
  !If it is none of the above, it has to be an auxiliary property
  DO i=1,SIZE(AUXNAMES)
    IF( TRIM(ADJUSTL(AUXNAMES(i)))==TRIM(ADJUSTL(sortcol)) ) THEN
      sortnum = i
    ENDIF
  ENDDO
  WRITE(msg,*) 'Array to sort: AUX ; sortnum = ', sortnum
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  IF(sortnum==0) THEN
    !no such property exist => display a warning and exit this module
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2730,(/TRIM(ADJUSTL(sortcol))/),(/0.d0/))
    GOTO 1000
  ENDIF
  !
  SELECT CASE(sortorder)
  CASE('up','down')
    !Sort AUX by increasing or decreasing values
    CALL QUICKSORT(AUX,sortnum,sortorder,newindex)
  CASE('pack')
    !Pack identical values of AUX together
    CALL PACKSORT(AUX,sortnum,newindex)
  CASE DEFAULT
  END SELECT
  !
  !newindex(:) now contains the new list of index (after sorting)
  !=> use it to sort other arrays
  CALL IDSORT(newindex,P)
  !
  !If shells are defined, sort them accordingly
  IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(newindex) ) THEN
    CALL IDSORT(newindex,S)
  ENDIF
  !
END SELECT
!
!
!
300 CONTINUE
CALL ATOMSK_MSG(2088,(/''/),(/0.d0/))
!
!
!
1000 CONTINUE
IF( ALLOCATED(newindex) ) DEALLOCATE(newindex)
!
!
END SUBROUTINE SORT_XYZ
!
!
!
END MODULE sort
