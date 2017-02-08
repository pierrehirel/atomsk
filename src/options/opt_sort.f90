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
!* Last modification: P. Hirel - 24 July 2014                                     *
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
INTEGER:: sortnum, sortcorrection
REAL(dp),DIMENSION(:,:),POINTER:: ArrayToSort  !pointer to the array to be sorted
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX    !auxiliary properties of atoms
REAL(dp),DIMENSION(:,:),ALLOCATABLE,TARGET:: P, S  !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,TARGET:: PandALL  !position + shells + aux. properties + SELECT
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
!its shells, auxiliary properties, and whether it is seleted or not
IF( ALLOCATED(AUX) .AND. SIZE(AUX(1,:))>0 ) THEN
  !Auxiliary properties are present
  IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
    !Auxiliary properties AND shells are present
    IF( ALLOCATED(SELECT) .AND. SIZE(SELECT)>0 ) THEN
      !PandALL(:,1:4) = P
      !PandALL(:,5:8) = S
      !PandALL(:,9)   = SELECT  (1.0 if selected, 0 otherwise)
      !PandALL(:,10:) = AUX
      sortcorrection = 9
      ALLOCATE( PandALL( SIZE(P,1), 9+SIZE(AUX,2) ) )
      PandALL(:,:) = 0.d0
      !Fill PandAUX with S and SELECT
      DO i=1,SIZE(P,1)
        DO j=1,4
          PandALL(i,4+j) = S(i,j)
        ENDDO
        IF( SELECT(i) ) THEN
          PandALL(i,9) = 1.d0
        ENDIF
      ENDDO
    ELSE
      !PandALL(:,1:4) = P
      !PandALL(:,5:8) = S
      !PandALL(:,9:)  = AUX
      sortcorrection = 8
      ALLOCATE( PandALL( SIZE(P,1), 8+SIZE(AUX,2) ) )
      PandALL(:,:) = 0.d0
      !Fill PandAUX with S
      DO i=1,SIZE(P,1)
        DO j=1,4
          PandALL(i,4+j) = S(i,j)
        ENDDO
      ENDDO
    ENDIF
    !
  ELSE
    !Auxiliary properties are present, but no shell
    IF( ALLOCATED(SELECT) .AND. SIZE(SELECT)>0 ) THEN
      !PandALL(:,1:4) = P
      !PandALL(:,5)   = SELECT
      !PandALL(:,6:)  = AUX
      sortcorrection = 5
      ALLOCATE( PandALL( SIZE(P,1), 5+SIZE(AUX,2) ) )
      PandALL(:,:) = 0.d0
      IF( SELECT(i) ) THEN
        PandALL(i,5) = 1.d0
      ENDIF
    ELSE
      !PandALL(:,1:4) = P
      !PandALL(:,5:)  = AUX
      sortcorrection = 4
      ALLOCATE( PandALL( SIZE(P,1), 4+SIZE(AUX,2) ) )
      PandALL(:,:) = 0.d0
    ENDIF
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
  !Shells are present, but no auxiliary property
  IF( ALLOCATED(SELECT) .AND. SIZE(SELECT)>0 ) THEN
    !PandALL(:,1:4) = P
    !PandALL(:,5:8) = S
    !PandALL(:,9)   = SELECT
    sortcorrection = 1
    ALLOCATE( PandALL( SIZE(P,1), SIZE(P,2)+SIZE(S,2)+1 ) )
    PandALL(:,:) = 0.d0
    DO i=1,SIZE(P,1)
      DO j=1,4
        PandALL(i,j) = P(i,j)
        PandALL(i,4+j) = S(i,j)
      ENDDO
      IF( SELECT(i) ) THEN
        PandALL(i,9) = 1.d0
      ENDIF
    ENDDO
  ELSE
    !PandALL(:,1:4) = P
    !PandALL(:,5:8) = S
    sortcorrection = 0
    ALLOCATE( PandALL( SIZE(P,1), SIZE(P,2)+SIZE(S,2) ) )
    PandALL(:,:) = 0.d0
    DO i=1,SIZE(P,1)
      DO j=1,4
        PandALL(i,j) = P(i,j)
        PandALL(i,4+j) = S(i,j)
      ENDDO
    ENDDO
  ENDIF
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
  IF( ALLOCATED(SELECT) .AND. SIZE(SELECT)>0 ) THEN
    !PandALL(:,1:4) = P
    !PandALL(:,5)   = SELECT
    sortcorrection = 0
    ALLOCATE( PandALL( SIZE(P,1), SIZE(P,2)+1 ) )
    PandALL(:,:) = 0.d0
    DO i=1,SIZE(P,1)
      DO j=1,4
        PandALL(i,j) = P(i,j)
      ENDDO
      IF( SELECT(i) ) THEN
        PandALL(i,5) = 1.d0
      ENDIF
    ENDDO
    !The array to sort becomes PandALL
    ArrayToSort=>PandALL
  ELSE
    ArrayToSort=>P
  ENDIF
  !
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
IF( .NOT.ASSOCIATED(ArrayToSort) ) THEN
  !For some reason the pointer does not point to anything => stop
  nerr = nerr+1
  GOTO 1000
ENDIF
!
IF(SIZE(P(:,1))>10000000) THEN
  CALL ATOMSK_MSG(3,(/''/),(/0.d0/))
ENDIF
!
!
!Sort atoms: algorithms are in "subroutines.f90"
SELECT CASE(sortorder)
CASE('up','down')
  WRITE(msg,*) 'Calling QUICKSORT...'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !CALL BUBBLESORT(ArrayToSort,sortnum,sortorder)
  CALL QUICKSORT(ArrayToSort,sortnum,sortorder)
CASE('pack')
  WRITE(msg,*) 'Calling PACKSORT...'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  CALL PACKSORT(ArrayToSort,sortnum)
CASE DEFAULT
END SELECT
!
!
!
300 CONTINUE
!If relevant, separate PandALL into P, S, AUX, SELECT
WRITE(msg,*) 'Separating arrays...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
IF( ALLOCATED(AUX) .AND. SIZE(AUX,2)>0 ) THEN
  !Auxiliary properties are present
  IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
    !Auxiliary properties AND shells are present
    IF( ALLOCATED(SELECT) .AND. SIZE(SELECT)>0 ) THEN
      !PandALL(:,1:4) = P
      !PandALL(:,5:8) = S
      !PandALL(:,9)   = SELECT  (1.0 if selected, 0 otherwise)
      !PandALL(:,10:) = AUX
      DO i=1,SIZE(P,1)
        DO j=1,4
          P(i,j) = PandALL(i,j)
          S(i,j) = PandALL(i,4+j)
        ENDDO
        IF( NINT(PandALL(i,9))==1 ) THEN
          SELECT(i) = .TRUE.
        ELSE
          SELECT(i) = .FALSE.
        ENDIF
        DO j=1,SIZE(AUX,2)
          AUX(i,j) = PandALL(i,9+j)
        ENDDO
      ENDDO
    ELSE
      !PandALL(:,1:4) = P
      !PandALL(:,5:8) = S
      !PandALL(:,9:) = AUX
      DO i=1,SIZE(P,1)
        DO j=1,4
          P(i,j) = PandALL(i,j)
          S(i,j) = PandALL(i,4+j)
        ENDDO
        DO j=1,SIZE(AUX,2)
          AUX(i,j) = PandALL(i,8+j)
        ENDDO
      ENDDO
    ENDIF
    !
  ELSE
    !Auxiliary properties are present, but no shell
    IF( ALLOCATED(SELECT) .AND. SIZE(SELECT)>0 ) THEN
      !PandALL(:,1:4) = P
      !PandALL(:,5)   = SELECT
      !PandALL(:,6:)  = AUX
      DO i=1,SIZE(P,1)
        DO j=1,4
          P(i,j) = PandALL(i,j)
        ENDDO
        IF( NINT(PandALL(i,5))==1 ) THEN
          SELECT(i) = .TRUE.
        ELSE
          SELECT(i) = .FALSE.
        ENDIF
        DO j=1,SIZE(AUX,2)
          AUX(i,j) = PandALL(i,4+j)
        ENDDO
      ENDDO
    ELSE
      !PandALL(:,1:4) = P
      !PandALL(:,5:)  = AUX
      DO i=1,SIZE(P,1)
        DO j=1,4
          P(i,j) = PandALL(i,j)
        ENDDO
        DO j=1,SIZE(AUX,2)
          AUX(i,j) = PandALL(i,4+j)
        ENDDO
      ENDDO
    ENDIF
    !
  ENDIF
  !
ELSEIF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
  !Shells are present, but no auxiliary property
  IF( ALLOCATED(SELECT) .AND. SIZE(SELECT)>0 ) THEN
    !PandALL(:,1:4) = P
    !PandALL(:,5:8) = S
    !PandALL(:,9)   = SELECT
    DO i=1,SIZE(P,1)
      DO j=1,4
        P(i,j) = PandALL(i,j)
        S(i,j) = PandALL(i,4+j)
      ENDDO
      IF( NINT(PandALL(i,9))==1 ) THEN
        SELECT(i) = .TRUE.
      ELSE
        SELECT(i) = .FALSE.
      ENDIF
    ENDDO
  ELSE
    !PandALL(:,1:4) = P
    !PandALL(:,5:8) = S
    DO i=1,SIZE(P,1)
      DO j=1,4
        P(i,j) = PandALL(i,j)
        S(i,j) = PandALL(i,4+j)
      ENDDO
    ENDDO
  ENDIF
  !
ELSE
  IF( ALLOCATED(SELECT) .AND. SIZE(SELECT)>0 ) THEN
    !PandALL(:,1:4) = P
    !PandALL(:,5)   = SELECT
    DO i=1,SIZE(P,1)
      DO j=1,4
        P(i,j) = PandALL(i,j)
      ENDDO
      IF( NINT(PandALL(i,5))==1 ) THEN
        SELECT(i) = .TRUE.
      ELSE
        SELECT(i) = .FALSE.
      ENDIF
    ENDDO
  ENDIF
ENDIF
!
!
CALL ATOMSK_MSG(2088,(/''/),(/0.d0/))
!
!
!
1000 CONTINUE
IF(ALLOCATED(PandALL)) DEALLOCATE(PandALL)
NULLIFY(ArrayToSort)
!
!
END SUBROUTINE SORT_XYZ
!
!
!
END MODULE sort
