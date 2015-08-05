MODULE swap
!
!**********************************************************************************
!*  SWAP                                                                          *
!**********************************************************************************
!* This module reads cartesian coordinates from an array and                      *
!* swaps atoms of given indices.                                                  *
!**********************************************************************************
!* (C) August 2015 - Pierre Hirel                                                 *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 05 Aug. 2015                                     *
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
SUBROUTINE SWAP_XYZ(P,S,AUX,swap_id)
!
!
IMPLICIT NONE
!
CHARACTER(LEN=128):: msg
INTEGER:: i
INTEGER,DIMENSION(2),INTENT(IN):: swap_id  !indices of atoms to swap
REAL(dp),DIMENSION(4):: Vtemp
REAL(dp),DIMENSION(:),ALLOCATABLE:: AUXtemp  !auxiliary properties (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of cores, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX  !auxiliary properties
!
!
!Initialize variables
i = 0
!
!
msg = 'Entering SWAP_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG( 2125, (/""/),(/DBLE(swap_id(1)),DBLE(swap_id(2))/) )
!
!If indices of atoms are out of bounds, abort
DO i=1,2
  IF( swap_id(i)<=0 .OR. swap_id(i)>SIZE(P,1) ) THEN
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2742,(/''/),(/DBLE(swap_id(i))/))
    GOTO 1000
  ENDIF
ENDDO
!
!
!Save position of first atom
Vtemp(:) = P(swap_id(1),:)
!
!Swap atoms
P(swap_id(1),:) = P(swap_id(2),:)
P(swap_id(2),:) = Vtemp(:)
!
!
!Swap shell positions (if any)
IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
  Vtemp(:) = S(swap_id(1),:)
  S(swap_id(1),:) = S(swap_id(2),:)
  S(swap_id(2),:) = Vtemp(:)
ENDIF
!
!Swap auxiliary properties (if any)
IF( ALLOCATED(AUX) .AND. SIZE(AUX,1)==SIZE(P,1) ) THEN
  ALLOCATE( AUXtemp( SIZE(AUX,2) ) )
  AUXtemp(:) = AUX(swap_id(1),:)
  AUX(swap_id(1),:) = AUX(swap_id(2),:)
  AUX(swap_id(2),:) = AUXtemp(:)
  DEALLOCATE(AUXtemp)
ENDIF
!
!
CALL ATOMSK_MSG(2126,(/''/),(/0.d0/))
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
END SUBROUTINE SWAP_XYZ
!
!
END MODULE swap
