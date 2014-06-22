MODULE rmprop
!
!**********************************************************************************
!*  RMPROP                                                                        *
!**********************************************************************************
!* This module removes some auxiliary properties.                                 *
!**********************************************************************************
!* (C) May 2011 - Pierre Hirel                                                    *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 15 Jan. 2014                                     *
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
SUBROUTINE RMPROP_XYZ(AUXNAMES,AUX,rmprop_prop,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=32),INTENT(IN):: rmprop_prop
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: newAUXNAMES            !names of auxiliary properties (temp.)
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: c1 !index of column to be removed from AUX and AUXNAMES
INTEGER:: i, j
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX            !auxiliary properties (temporary)
!
!Initialize variables
 c1=0
IF(ALLOCATED(newAUXNAMES)) DEALLOCATE(newAUXNAMES)
IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
!
WRITE(msg,*) 'Entering NOPROP_XYZ: '//TRIM(ADJUSTL(rmprop_prop))
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(2069,(/TRIM(ADJUSTL(rmprop_prop))/),(/0.d0/))
!
IF( .NOT. ALLOCATED(AUXNAMES) ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2729,(/TRIM(ADJUSTL(rmprop_prop))/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
100 CONTINUE
IF( TRIM(ADJUSTL(rmprop_prop))=='all' ) THEN
  IF( .NOT.ALLOCATED(SELECT) .OR. SIZE(SELECT)<=0 ) THEN
    !Remove all auxiliary properties
    DEALLOCATE(AUXNAMES)
    DEALLOCATE(AUX)
  ELSE
    !Some atoms are selected => just set auxiliary properties = 0 for these atoms
    DO i=1,SIZE(AUX,1)
      IF( SELECT(i) ) THEN
        AUX(i,:) = 0.d0
      ENDIF
    ENDDO
  ENDIF
  !
ELSE
  !Remove only auxliary property "rmprop_prop"
  !Find column containing this property
  DO i=1,SIZE(AUXNAMES)
    IF( TRIM(ADJUSTL(AUXNAMES(i)))==TRIM(ADJUSTL(rmprop_prop)) ) THEN
      c1 = i
    ENDIF
  ENDDO
  !
  IF(c1<=0) THEN
    !No property with that name => skip
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2730,(/TRIM(ADJUSTL(rmprop_prop))/),(/0.d0/))
    GOTO 1000
    !
  ELSE
    IF( .NOT.ALLOCATED(SELECT) .OR. SIZE(SELECT)<=0 ) THEN
      !One column must be removed from AUX and AUXNAMES
      IF( SIZE(AUXNAMES)<=1 ) THEN
        !This is the only property => simply deallocate arrays
        DEALLOCATE(AUXNAMES)
        DEALLOCATE(AUX)
        !
      ELSE
        !Remove only the column
        ALLOCATE( newAUXNAMES( SIZE(AUXNAMES)-1 ) )
        newAUXNAMES(:) = ''
        ALLOCATE( newAUX( SIZE(AUX(:,1)), SIZE(AUX(1,:))-1 ) )
        newAUX(:,:) = 0.d0
        !Save values to newAUX
        j=0
        DO i=1,SIZE(AUXNAMES)
          IF(i.NE.c1) THEN
            j=j+1
            newAUXNAMES(j) = AUXNAMES(i)
            newAUX(:,j) = AUX(:,i)
          ENDIF
        ENDDO
        !
        !Replace old AUX by newAUX
        DEALLOCATE(AUX)
        ALLOCATE( AUX( SIZE(newAUX(:,1)), SIZE(newAUX(1,:)) ) )
        AUX = newAUX
        DEALLOCATE(newAUX)
        !Same with newAUXNAMES
        DEALLOCATE(AUXNAMES)
        ALLOCATE( AUXNAMES( SIZE(newAUXNAMES) ) )
        AUXNAMES = newAUXNAMES
        DEALLOCATE(newAUXNAMES)
      ENDIF
      !
    ELSE
      !Some atoms are selected => set property=0 for these atoms
      DO i=1,SIZE(AUX,1)
        IF( SELECT(i) ) THEN
          AUX(i,c1) = 0.d0
        ENDIF
      ENDDO
    ENDIF
    !
  ENDIF  !endif c1==0
ENDIF  !endif rmprop_prop==all
!
!
!
400 CONTINUE
CALL ATOMSK_MSG(2070,(/TRIM(ADJUSTL(rmprop_prop))/),(/0.d0/))
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE RMPROP_XYZ
!
!
!
END MODULE rmprop