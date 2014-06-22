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
!* Last modification: P. Hirel - 25 Sept. 2013                                    *
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
SUBROUTINE SORT_XYZ(P,S,AUXNAMES,AUX,sortcol,sortorder)
!
!
IMPLICIT NONE
CHARACTER(LEN=4):: sortorder  !up or down
CHARACTER(LEN=16):: sortcol    !property to be sorted: x, y, z or s, or any name
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
INTEGER:: i, j
INTEGER:: sortnum, sortcorrection
REAL(dp),DIMENSION(:,:),POINTER:: ArrayToSort  !pointer to the array to be sorted
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX    !auxiliary properties of atoms
REAL(dp),DIMENSION(:,:),ALLOCATABLE,TARGET:: P, S  !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,TARGET:: PandALL  !position + shells + aux. properties
!
!Initialize variables
sortcorrection = 0
sortnum = 0
IF(ALLOCATED(PandALL)) DEALLOCATE(PandALL)
ArrayToSort=>P
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
!If auxiliary properties or shells are present, merge everything together in "PandALL"
!This is required to ensure that after the sorting, each atom conserves
!its shells and/or auxiliary properties
IF( ALLOCATED(AUX) .AND. SIZE(AUX(1,:))>0 ) THEN
  !
  IF(ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
    !If auxiliary properties AND shells are present, merge all of them to PandALL
    sortcorrection = 8
    ALLOCATE( PandALL( SIZE(P,1), 8+SIZE(AUX,2) ) )
    PandALL(:,:) = 0.d0
    !Fill PandAUX with S
    DO i=1,SIZE(P,1)
      DO j=1,4
        PandALL(i,4+j) = S(i,j)
      ENDDO
    ENDDO
    !
  ELSE
    !No shells, only auxiliary properties => merge P and AUX
    sortcorrection = 4
    ALLOCATE( PandALL( SIZE(P,1), 4+SIZE(AUX,2) ) )
    PandALL(:,:) = 0.d0
  ENDIF
  !
  !Fill PandALL with P and AUX
  DO i=1,SIZE(P,1)
    DO j=1,4
      PandALL(i,j) = P(i,j)
    ENDDO
    DO j=1,SIZE(AUX,2)
      PandALL(i,sortcorrection+j) = AUX(i,j)
    ENDDO
  ENDDO
  !
  !The array to sort becomes PandALL
  ArrayToSort=>PandALL
  !
  !Find which property must be sorted
  SELECT CASE (sortcol)
  CASE('species','s','S')
    sortnum = 4
  CASE('x','X')
    sortnum = 1
  CASE('y','Y')
    sortnum = 2
  CASE('z','Z')
    sortnum = 3
  CASE DEFAULT
    !If it is none of the above, it has to be an auxiliary property
    DO i=1,SIZE(AUXNAMES)
      IF( TRIM(ADJUSTL(AUXNAMES(i)))==TRIM(ADJUSTL(sortcol)) ) THEN
        sortnum = i
      ENDIF
    ENDDO
    IF(sortnum==0) THEN
      !no such property exist
      nwarn=nwarn+1
      CALL ATOMSK_MSG(2730,(/TRIM(ADJUSTL(sortcol))/),(/0.d0/))
      GOTO 1000
    ELSE
      !correct the column number to match the array PandALL
      sortnum = sortnum + sortcorrection
    ENDIF
  END SELECT
!
ELSEIF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
  !No auxiliary properties, only shells => merge P and S
  sortcorrection = 0
  ALLOCATE( PandALL( SIZE(P,1), SIZE(P,2)+SIZE(S,2) ) )
  DO i=1,SIZE(P,1)
    DO j=1,4
      PandALL(i,j) = P(i,j)
      PandALL(i,4+j) = S(i,j)
    ENDDO
  ENDDO
  !The array to sort becomes PandALL
  ArrayToSort=>PandALL
  !Find which property must be sorted
  SELECT CASE (sortcol)
  CASE('species','s','S')
    sortnum = 4
  CASE('x','X')
    sortnum = 1
  CASE('y','Y')
    sortnum = 2
  CASE('z','Z')
    sortnum = 3
  CASE DEFAULT
    CALL ATOMSK_MSG(2732,(/TRIM(ADJUSTL(sortcol))/),(/0.d0/))
    GOTO 1000
  END SELECT
!
ELSE
  !Simple case: just sort P
  SELECT CASE (sortcol)
  CASE('species','s','S')
    sortnum = 4
  CASE('x','X')
    sortnum = 1
  CASE('y','Y')
    sortnum = 2
  CASE('z','Z')
    sortnum = 3
  CASE DEFAULT
    CALL ATOMSK_MSG(2732,(/TRIM(ADJUSTL(sortcol))/),(/0.d0/))
    GOTO 1000
  END SELECT
ENDIF
!
!
!
200 CONTINUE
IF(SIZE(P(:,1))>500000) THEN
  CALL ATOMSK_MSG(3,(/''/),(/0.d0/))
ENDIF
!
WRITE(msg,*) 'Sorting...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Sort atoms: algorithms are in "subroutines.f90"
SELECT CASE(sortorder)
CASE('up','down')
  CALL BUBBLESORT(ArrayToSort,sortnum,sortorder)
CASE('pack')
  CALL PACKSORT(ArrayToSort,sortnum)
CASE DEFAULT
END SELECT
!
!
!
300 CONTINUE
!If relevant, separate PandALL into P, S and/or AUX
IF( ALLOCATED(AUX) .AND. SIZE(AUX(1,:))>0 ) THEN
  !
  WRITE(msg,*) 'Separating arrays...'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  IF(ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
    !Separate shells and auxiliary properties
    DO i=1,SIZE(P,1)
      DO j=1,4
        P(i,j) = PandALL(i,j)
        S(i,j) = PandALL(i,4+j)
      ENDDO
      DO j=1,SIZE(AUX,2)
        AUX(i,j) = PandALL(i,8+j)
      ENDDO
    ENDDO
    !
  ELSE
    !Separate auxiliary properties
    DO i=1,SIZE(P,1)
      DO j=1,4
        P(i,j) = PandALL(i,j)
      ENDDO
      DO j=1,SIZE(AUX,2)
        AUX(i,j) = PandALL(i,SIZE(P,2)+j)
      ENDDO
    ENDDO
  ENDIF
  !
ELSEIF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
  !
  WRITE(msg,*) 'Separating arrays P and S...'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !Separate cores and shells
  DO i=1,SIZE(P,1)
    DO j=1,4
      P(i,j) = PandALL(i,j)
      S(i,j) = PandALL(i,4+j)
    ENDDO
  ENDDO
ENDIF
!
!
CALL ATOMSK_MSG(2088,(/''/),(/0.d0/))
!
!
!
1000 CONTINUE
IF(ALLOCATED(PandALL)) DEALLOCATE(PandALL)
!
!
END SUBROUTINE SORT_XYZ
!
!
!
END MODULE sort
