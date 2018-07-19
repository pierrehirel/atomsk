MODULE sorting
!
!**********************************************************************************
!*  SORTING                                                                       *
!**********************************************************************************
!* This module contains subroutines for sorting arrays.                           *
!**********************************************************************************
!* (C) March 2018 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 05 March 2018                                    *
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
!* BUBBLESORT          sorts an array by increasing or decreasing values          *
!* QUICKSORT           sorts an array by increasing or decreasing values          *
!* PACKSORT            sorts an array by packing identical values together        *
!* IDSORT              sorts an array according to the provided list of index     *
!**********************************************************************************
!
!
USE comv
USE constants
USE functions
USE math
!
IMPLICIT NONE
!
!
CONTAINS
!
!
!********************************************************
! BUBBLESORT
! This subroutine sorts a MxN array by increasing or
! decreasing values of its column 'col'.
! It also returns the list of index after the sorting.
! It uses the so-called "bubble-sort" algorithm.
! IMPORTANT NOTE: this algorithm swaps entire columns
! of an array, up to N*(N-1) times (N=number of lines),
! i.e. it scales as N². For very large arrays, use
! the QUICKSORT algorithm below.
!********************************************************
SUBROUTINE BUBBLESORT(A,col,order,newindex)
!
IMPLICIT NONE
CHARACTER(LEN=4),INTENT(IN):: order  !up or down
LOGICAL:: sorted
INTEGER:: i, j, k
INTEGER,INTENT(IN):: col             !index of column to sort
REAL(dp),DIMENSION(:,:),INTENT(INOUT):: A
REAL(dp),DIMENSION(SIZE(A,2)) :: col_value
INTEGER,DIMENSION(:),ALLOCATABLE:: newindex !list of sorted indexes
!
!
IF( SIZE(A,2)>=1 .AND. col>0 .AND. col<=SIZE(A,2) ) THEN
  col_value(:) = 0.d0
  IF(ALLOCATED(newindex)) DEALLOCATE(newindex)
  ALLOCATE( newindex(SIZE(A,1)) )
  DO i=1,SIZE(newindex)
    newindex(i) = i
  ENDDO
  !
  IF(order=='down') THEN
    DO j=SIZE(A,1),2,-1
      sorted = .TRUE.
      DO i=1,j-1
        !If element i+1 is greater than element i, swap them
        IF( A(i+1,col) > A(i,col) ) THEN
          col_value(:) = A(i,:)
          A(i,:) = A(i+1,:)
          A(i+1,:) = col_value(:)
          !Save new indexes
          k = newindex(i)
          newindex(i) = newindex(i+1)
          newindex(i+1) = k
          !We performed an inversion => list is not sorted
          sorted = .FALSE.
        ENDIF
      ENDDO
      !If no inversion was performed after loop on i, the list is fully sorted
      IF(sorted) EXIT
    ENDDO
    !
  ELSE  !i.e. if order is "up"
    DO j=SIZE(A,1),2,-1
      sorted = .TRUE.
      DO i=1,j-1
        !If element i+1 is smaller than element i, swap them
        IF( A(i+1,col) < A(i,col) ) THEN
          col_value(:) = A(i,:)
          A(i,:) = A(i+1,:)
          A(i+1,:) = col_value(:)
          !Save new indexes
          k = newindex(i)
          newindex(i) = newindex(i+1)
          newindex(i+1) = k
          !We performed an inversion => list is not sorted
          sorted = .FALSE.
        ENDIF
      ENDDO
      !If no inversion was performed after loop on i, the list is fully sorted
      IF(sorted) EXIT
    ENDDO
  ENDIF
  !
ELSE
  PRINT*, "ERROR col = ", col
ENDIF
!
END SUBROUTINE BUBBLESORT
!
!
!********************************************************
! QUICKSORT
! This subroutine sorts a MxN array by increasing or
! decreasing values of its column 'col'.
! It also returns the list of index after the sorting.
! It uses the so-called "quick-sort" algorithm, which
! scales as N*log(N) where N is the number of lines.
! NOTE: this algorithm actually counts three subroutines:
!      QUICKSORT (the one that you should CALL),
!      QSORT (which is recursive), and
!      QS_PARTITION below.
!********************************************************
SUBROUTINE QUICKSORT(A,col,order,newindex)
!
IMPLICIT NONE
CHARACTER(LEN=4),INTENT(IN):: order  !up or down
INTEGER,INTENT(IN):: col             !index of column to sort
INTEGER:: i
INTEGER,DIMENSION(:),ALLOCATABLE:: newindex !list of sorted indexes
REAL(dp),DIMENSION(:,:),INTENT(INOUT):: A
!
IF(ALLOCATED(newindex)) DEALLOCATE(newindex)
ALLOCATE( newindex(SIZE(A,1)) )
DO i=1,SIZE(newindex)
  newindex(i) = i
ENDDO
!
CALL QSORT(A,col,order,newindex)
!
END SUBROUTINE QUICKSORT
!
!********************************************************
RECURSIVE SUBROUTINE QSORT(A,col,order,newindex)
!
IMPLICIT NONE
CHARACTER(LEN=4),INTENT(IN):: order  !up or down
INTEGER,INTENT(IN):: col             !index of column to sort
INTEGER:: iq
INTEGER,DIMENSION(:):: newindex !list of sorted indexes
REAL(dp),DIMENSION(:,:),INTENT(INOUT):: A
!
IF( SIZE(A,1) > 1 ) THEN
  CALL QS_PARTITION(A,col,order,iq,newindex)
  CALL QSORT(A(:iq-1,:),col,order,newindex(:iq-1))
  CALL QSORT(A(iq:,:),col,order,newindex(iq:))
ENDIF
!
END SUBROUTINE QSORT
!
!********************************************************
!
SUBROUTINE QS_PARTITION(A,col,order,marker,newindex)
!
IMPLICIT NONE
CHARACTER(LEN=4),INTENT(IN):: order  !up or down
REAL(dp),INTENT(INOUT),DIMENSION(:,:):: A
INTEGER,INTENT(OUT):: marker
INTEGER,INTENT(IN):: col             !index of column to sort
INTEGER:: i, j, k
INTEGER,DIMENSION(:):: newindex !list of sorted indexes
REAL(dp),DIMENSION(SIZE(A,2)):: tempreal
REAL(dp):: pivot      !value of pivot point
!
pivot = A(1,col)
i = 0
j = SIZE(A,1)+1
!
IF(order=='up') THEN
  DO
    j = j-1
    DO
      IF( A(j,col) <= pivot ) THEN
        EXIT
      ENDIF
      j = j-1
    ENDDO
    !
    i = i+1
    DO
      IF( A(i,col) >= pivot ) THEN
        EXIT
      ENDIF
      i = i+1
    ENDDO
    !
    IF( i < j ) THEN
      ! Exchange A(i,:) and A(j,:)
      tempreal(:) = A(i,:)
      A(i,:) = A(j,:)
      A(j,:) = tempreal(:)
      k = newindex(i)
      newindex(i) = newindex(j)
      newindex(j) = k
    ELSEIF( i==j ) THEN
      marker = i+1
      RETURN
    ELSE
      marker = i
      RETURN
    ENDIF
    !
  ENDDO
  !
ELSE  !i.e. if order != 'up'
  DO
    j = j-1
    DO
      IF( A(j,col) >= pivot ) THEN
        EXIT
      ENDIF
      j = j-1
    ENDDO
    !
    i = i+1
    DO
      IF( A(i,col) <= pivot ) THEN
        EXIT
      ENDIF
      i = i+1
    ENDDO
    !
    IF( i < j ) THEN
      ! Exchange A(i,:) and A(j,:)
      tempreal(:) = A(i,:)
      A(i,:) = A(j,:)
      A(j,:) = tempreal(:)
      k = newindex(i)
      newindex(i) = newindex(j)
      newindex(j) = k
    ELSEIF( i==j ) THEN
      marker = i+1
      RETURN
    ELSE
      marker = i
      RETURN
    ENDIF
    !
  ENDDO
  !
ENDIF  
!
END SUBROUTINE QS_PARTITION
!
!
!********************************************************
! PACKSORT
! This subroutine sorts a MxN array by packing
! identical values together, but without sorting them
! (unlike the BUBBLESORT above). E.g. if values like:
!      2 4 2 3 1 1 2 4 3 3 2 1 4
! are given, this subroutine will pack identical
! values so the result is:
!      2 2 2 2 4 4 4 3 3 3 1 1 1
! It also returns the list of index after the sorting.
!********************************************************
SUBROUTINE PACKSORT(A,col,newindex)
!
IMPLICIT NONE
INTEGER:: i, j, k, l
INTEGER:: col, last
INTEGER,DIMENSION(:),ALLOCATABLE:: newindex !list of sorted indexes
REAL(dp),DIMENSION(:,:),INTENT(INOUT):: A
REAL(dp),DIMENSION(SIZE(A,2)):: Atemp
!
IF(ALLOCATED(newindex)) DEALLOCATE(newindex)
ALLOCATE( newindex(SIZE(A,1)) )
DO i=1,SIZE(newindex)
  newindex(i) = i
ENDDO
!
IF(col>SIZE(A,2)) col = SIZE(A,2)
!
DO i=1,SIZE(A,1)
  last=i
  DO j=i+1,SIZE(A,1)
    IF( A(j,col)==A(i,col) ) THEN
      !Only consider A(j)==A(i)
      IF(j==last+1) THEN
        !If the two values are contiguous just go on
        last=j
      ELSE
        !If the value is further down in the array, pack it with the others
        Atemp(:) = A(j,:)
        l = newindex(j)
        !Shift all previous values in the array
        DO k=j,last+1,-1
          A(k,:) = A(k-1,:)
          newindex(k) = newindex(k-1)
        ENDDO
        last=last+1
        A(last,:) = Atemp(:)
        newindex(last) = l
      ENDIF
    ENDIF
  ENDDO
ENDDO

!
END SUBROUTINE PACKSORT
!
!
!********************************************************
! IDSORT
! This subroutine takes a list of index and an 2-dim.
! array A as input, and re-shuffles array A according
! to the index list.
!********************************************************
SUBROUTINE IDSORT(idlist,A)
!
INTEGER:: i
INTEGER,DIMENSION(:),INTENT(IN):: idlist !list of sorted indexes
REAL(dp),DIMENSION(:,:),INTENT(INOUT):: A
REAL(dp),DIMENSION(SIZE(A,1),SIZE(A,2)):: Atemp
!
DO i=1,SIZE(idlist)
  !Exchange entries in A with index i and idlist(i)
  Atemp(i,:) = A(idlist(i),:)
ENDDO
!
A(:,:) = Atemp(:,:)
!
END SUBROUTINE IDSORT
!
!
!
END MODULE sorting