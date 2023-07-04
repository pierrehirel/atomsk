MODULE resize
!
!**********************************************************************************
!*  RESIZE                                                                        *
!**********************************************************************************
!* This module contains general subroutines to resize arrays.                     *
!**********************************************************************************
!* (C) April 2018 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 21 June 2023                                     *
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
!* RESIZE_DBLEARRAY2   resize a real 2-D array, keeping its content               *
!* RESIZE_INTARRAY2    resize an integer 2-D array, keeping its content           *
!* RESIZE_LOGICAL1     resize a boolean 1-D array, keeping its content            *
!**********************************************************************************
!
!
USE comv
USE messages
!
IMPLICIT NONE
!
!
CONTAINS
!
!
!
!********************************************************
! RESIZE_DBLEARRAY2
! This routine changes the size of the provided Array.
! Initial data is preserved in the new array.
! If the new size is larger, unknown data is set to zero.
! If the new size is smaller, some data is lost.
!********************************************************
SUBROUTINE RESIZE_DBLEARRAY2(Array,L1,L2,status)
!
IMPLICIT NONE
INTEGER:: i, j
INTEGER,OPTIONAL:: status  !Success=0; failure=1
INTEGER,INTENT(IN):: L1, L2  !new sizes of Array
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: Array !the array to resize
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: temp_array !temporary copy of Array
!
IF(PRESENT(status)) status = 0
!
IF( .NOT.ALLOCATED(Array) ) THEN
  !Allocate Array with required size and fill it with zeros
  ALLOCATE(Array(L1,L2),STAT=i)
  IF( i>0 ) THEN
    ! Allocation failed (not enough memory)
    IF(PRESENT(status)) status = 1
    nerr = nerr+1
    CALL ATOMSK_MSG(819,(/''/),(/DBLE(L1)/))
    RETURN
  ENDIF
  Array(:,:) = 0.d0
  !
ELSE
  !Array is already allocated => resize it (only if new sizes are different from previous ones)
  IF( L1>0 .AND. L2>0 .AND. (L1.NE.SIZE(Array,1) .OR. L2.NE.SIZE(Array,2)) ) THEN
    !
    ALLOCATE( temp_array(L1,L2) , STAT=i )
    IF( i>0 ) THEN
      ! Allocation failed (not enough memory)
      IF(PRESENT(status)) status = 1
      nerr = nerr+1
      CALL ATOMSK_MSG(819,(/''/),(/DBLE(SIZE(Array,1)+L1)/))
      RETURN
    ENDIF
    !
    temp_array(:,:) = 0.d0
    DO i=1,MIN(L1,SIZE(Array,1))
      DO j=1,MIN(L2,SIZE(Array,2))
        temp_array(i,j) = Array(i,j)
      ENDDO
    ENDDO
    !
    DEALLOCATE(Array)
    ALLOCATE( Array(L1,L2) , STAT=i )
    IF( i>0 ) THEN
      ! Allocation failed (not enough memory)
      IF(PRESENT(status)) status = 1
      nerr = nerr+1
      CALL ATOMSK_MSG(819,(/''/),(/DBLE(SIZE(Array,1)+L1)/))
      RETURN
    ENDIF
    Array(:,:) = temp_array(:,:)
    !
    DEALLOCATE(temp_array)
    !
  ELSE
    !i.e. if L1<=0 or L2<=0 => problem
    IF(PRESENT(status)) status = 1
  ENDIF
  !
ENDIF
!
END SUBROUTINE RESIZE_DBLEARRAY2
!
!
!
!********************************************************
! RESIZE_INTARRAY2
! This routine changes the size of the provided Array.
! Initial data is preserved in the new array.
! If the new size is larger, unknown data is set to zero.
! If the new size is smaller, some data is lost.
!********************************************************
SUBROUTINE RESIZE_INTARRAY2(Array,L1,L2,status)
!
IMPLICIT NONE
INTEGER:: i, j
INTEGER,OPTIONAL:: status  !Success=0; failure=1
INTEGER,INTENT(IN):: L1, L2  !new sizes of Array
INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: Array !the array to resize
INTEGER,DIMENSION(:,:),ALLOCATABLE:: temp_array !temporary copy of Array
!
IF(PRESENT(status)) status = 0
!
IF( .NOT.ALLOCATED(Array) ) THEN
  !Allocate Array with required size and fill it with zeros
  ALLOCATE(Array(L1,L2),STAT=i)
  IF( i>0 ) THEN
    ! Allocation failed (not enough memory)
    IF(PRESENT(status)) status = 1
    nerr = nerr+1
    CALL ATOMSK_MSG(819,(/''/),(/0.d0/))
    RETURN
  ENDIF
  Array(:,:) = 0
  !
ELSE
  !Array is already allocated => resize it
  IF( L1>0 .AND. L2>0 .AND. (L1.NE.SIZE(Array,1) .OR. L2.NE.SIZE(Array,2)) ) THEN
    !
    ALLOCATE( temp_array(L1,L2) , STAT=i )
    IF( i>0 ) THEN
      ! Allocation failed (not enough memory)
      IF(PRESENT(status)) status = 1
      nerr = nerr+1
      CALL ATOMSK_MSG(819,(/''/),(/0.d0/))
      RETURN
    ENDIF
    !
    temp_array(:,:) = 0
    DO i=1,MIN(L1,SIZE(Array,1))
      DO j=1,MIN(L2,SIZE(Array,2))
        temp_array(i,j) = Array(i,j)
      ENDDO
    ENDDO
    !
    DEALLOCATE(Array)
    ALLOCATE( Array(L1,L2) , STAT=i )
    IF( i>0 ) THEN
      ! Allocation failed (not enough memory)
      IF(PRESENT(status)) status = 1
      nerr = nerr+1
      CALL ATOMSK_MSG(819,(/''/),(/0.d0/))
      RETURN
    ENDIF
    Array(:,:) = temp_array(:,:)
    !
    DEALLOCATE(temp_array)
    !
  ELSE
    !i.e. if L1<=0 or L2<=0 => problem
    IF(PRESENT(status)) status = 1
  ENDIF
  !
ENDIF
!
END SUBROUTINE RESIZE_INTARRAY2
!
!
!
!********************************************************
! RESIZE_LOGICAL1
! This routine changes the size of the provided logical Array.
! Initial data is preserved in the new array.
! If the new size is larger, unknown data is set to .FALSE.
! If the new size is smaller, some data is lost.
!********************************************************
SUBROUTINE RESIZE_LOGICAL1(Array,L1,status)
!
IMPLICIT NONE
INTEGER:: i
INTEGER,OPTIONAL:: status  !Success=0; failure=1
INTEGER,INTENT(IN):: L1    !new size of Array
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: Array !the array to resize
LOGICAL,DIMENSION(:),ALLOCATABLE:: temp_array !temporary copy of Array
!
IF(PRESENT(status)) status = 0
!
IF( .NOT.ALLOCATED(Array) ) THEN
  !Allocate Array with required size and fill it with .FALSE.
  ALLOCATE(Array(L1))
  Array(:) = .FALSE.
  !
ELSE
  !Array is already allocated => resize it
  IF( L1>0 .AND. L1.NE.SIZE(Array) ) THEN
    !
    ALLOCATE( temp_array(L1),STAT=i )
    IF( i>0 ) THEN
      ! Allocation failed (not enough memory)
      IF(PRESENT(status)) status = 1
      nerr = nerr+1
      CALL ATOMSK_MSG(819,(/''/),(/0.d0/))
      RETURN
    ENDIF
    !
    temp_array(:) = .FALSE.
    DO i=1,MIN(L1,SIZE(Array))
      temp_array(i) = Array(i)
    ENDDO
    !
    DEALLOCATE(Array)
    ALLOCATE( Array(L1) ,STAT=i)
    IF( i>0 ) THEN
      ! Allocation failed (not enough memory)
      IF(PRESENT(status)) status = 1
      nerr = nerr+1
      CALL ATOMSK_MSG(819,(/''/),(/0.d0/))
      RETURN
    ENDIF
    Array(:) = temp_array(:)
    !
    DEALLOCATE(temp_array)
    !
  ELSE
    !i.e. if L1<=0 or L2<=0 => problem
    IF(PRESENT(status)) status = 1
  ENDIF
  !
ENDIF
!
END SUBROUTINE RESIZE_LOGICAL1
!
!
!
END MODULE resize
