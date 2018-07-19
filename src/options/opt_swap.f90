MODULE swap
!
!**********************************************************************************
!*  SWAP                                                                          *
!**********************************************************************************
!* This module reads cartesian coordinates from an array and swaps atoms          *
!* of given indices, or swaps the two given Cartesian axes.                       *
!**********************************************************************************
!* (C) August 2015 - Pierre Hirel                                                 *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 05 April 2016                                    *
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
SUBROUTINE SWAP_XYZ(H,P,S,AUX,swap_id)
!
!
IMPLICIT NONE
!
CHARACTER(LEN=128):: msg
INTEGER:: i
INTEGER,DIMENSION(2):: id   !indices of axis or atoms to swap
CHARACTER(LEN=16),DIMENSION(2),INTENT(IN):: swap_id  !Cartesian axes or indices of atoms to swap
REAL(dp),DIMENSION(4):: Vtemp
REAL(dp),DIMENSION(:),ALLOCATABLE:: AUXtemp  !auxiliary properties (temporary)
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of cores, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX  !auxiliary properties
!
!
!Initialize variables
i = 0
id(:)=0
!
!
msg = 'Entering SWAP_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG( 2125, (/swap_id(1),swap_id(2)/),(/0.d0/) )
!
IF( swap_id(1) == swap_id(2) ) THEN
  !There is nothing to exchange => slip
  nwarn=nwarn+1
  CALL ATOMSK_MSG( 2757, (/""/),(/0.d0/) )
  GOTO 1000
ENDIF
!
SELECT CASE(swap_id(1))
!
CASE('X','x','Y','y','Z','z')
  !Two Cartesian axes must be swapped
  !
  IF( swap_id(1)=="X" .OR. swap_id(1)=="x" ) THEN
    id(1) = 1
  ELSEIF( swap_id(1)=="Y" .OR. swap_id(1)=="y" ) THEN
    id(1) = 2
  ELSE
    id(1) = 3
  ENDIF
  !
  SELECT CASE(swap_id(2))
  CASE('X','x','Y','y','Z','z')
    IF( swap_id(2)=="X" .OR. swap_id(2)=="x" ) THEN
      id(2) = 1
    ELSEIF( swap_id(2)=="Y" .OR. swap_id(2)=="y" ) THEN
      id(2) = 2
    ELSE
      id(2) = 3
    ENDIF
    !
    !Swap the two axes
    Vtemp(1:3) = H(id(1),1:3)
    H(id(1),id(1)) = H(id(2),id(2))
    H(id(1),id(2)) = H(id(2),id(1))
    H(id(2),id(1)) = Vtemp(id(2))
    H(id(2),id(2)) = Vtemp(id(1))
    !
    !Swap coordinates of each atom
    DO i=1,SIZE(P,1)
      Vtemp(:) = P(i,:)
      P(i,id(1)) = Vtemp(id(2))
      P(i,id(2)) = Vtemp(id(1))
    ENDDO
    !
  CASE DEFAULT
    !swap_id2 is not a Cartesian axis => big problem!
    !(this should have been dealt with before and not happen here,
    !however if for some reason it does happen let's ensure a smooth escape)
    nerr=nerr+1
  END SELECT
  !
  !
CASE DEFAULT
  !Two atoms must be swapped
  !Read indices of atoms to swap
  !(at this point swap_id(:) *must* contain integer numbers)
  READ(swap_id(1),*,ERR=800,END=800) id(1)
  READ(swap_id(2),*,ERR=800,END=800) id(2)
  !
  !If indices of atoms are out of bounds, abort
  DO i=1,2
    IF( id(i)<=0 .OR. id(i)>SIZE(P,1) ) THEN
      nwarn=nwarn+1
      CALL ATOMSK_MSG(2742,(/''/),(/DBLE(id(i))/))
      GOTO 1000
    ENDIF
  ENDDO
  !
  !Save position of first atom
  Vtemp(:) = P(id(1),:)
  !
  !Swap atoms
  P(id(1),:) = P(id(2),:)
  P(id(2),:) = Vtemp(:)
  !
  !Swap shell positions (if any)
  IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
    Vtemp(:) = S(id(1),:)
    S(id(1),:) = S(id(2),:)
    S(id(2),:) = Vtemp(:)
  ENDIF
  !
  !Swap auxiliary properties (if any)
  IF( ALLOCATED(AUX) .AND. SIZE(AUX,1)==SIZE(P,1) ) THEN
    ALLOCATE( AUXtemp( SIZE(AUX,2) ) )
    AUXtemp(:) = AUX(id(1),:)
    AUX(id(1),:) = AUX(id(2),:)
    AUX(id(2),:) = AUXtemp(:)
    DEALLOCATE(AUXtemp)
  ENDIF
  !
END SELECT
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
