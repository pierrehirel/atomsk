MODULE rmshells
!
!**********************************************************************************
!*  RMSHELLS                                                                      *
!**********************************************************************************
!* This module removes shells on one type of atoms, or on all atoms.              *
!**********************************************************************************
!* (C) May 2014 - Pierre Hirel                                                    *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 01 March 2017                                    *
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
USE atoms
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
SUBROUTINE RMSHELLS_XYZ(P,S,rmshells_prop,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=2):: species
CHARACTER(LEN=6),INTENT(IN):: rmshells_prop !atom species on which shells should be removed
CHARACTER(LEN=128):: msg
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT  !mask for atom list
INTEGER:: i
INTEGER:: rmshells !number of shells removed
REAL(dp):: snumber !atomic number
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P      !positions of atoms
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: S   !positions of shells
!
species = ''
i = 0
rmshells = 0
snumber = 0.d0
!
WRITE(msg,*) 'Entering RMSHELLS_XYZ: '//TRIM(ADJUSTL(rmshells_prop))
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(2122,(/rmshells_prop/),(/0.d0/))
!
!
!If there are no shells, nothing to be done
IF( .NOT.ALLOCATED(S) ) THEN
  CALL ATOMSK_MSG(2744,(/""/),(/0.d0/))
  !
  !
ELSE
  IF( rmshells_prop(1:3)=="all" .OR. rmshells_prop(1:3)=="any" ) THEN
    !Count how many shells will be removed
    DO i=1,SIZE(S,1)
      IF( S(i,4)>0.5d0 ) THEN
        rmshells = rmshells+1
      ENDIF
    ENDDO
    !Remove all shells
    IF(ALLOCATED(S)) DEALLOCATE(S)
    !
  ELSE
    !Remove shells of the given species
    species = ADJUSTL(rmshells_prop)
    !
    !Check that it is a species
    CALL ATOMNUMBER(species,snumber)
    IF( snumber<0.1d0 ) THEN
      !Unable to recognize atom species => error
      nerr=nerr+1
      CALL ATOMSK_MSG(801,(/rmshells_prop/),(/0.d0/))
      GOTO 1000
    ENDIF
    !
    !Remove shells that belong to that species *and* are selected
    DO i=1,SIZE(P,1)
      IF( NINT(P(i,4))==NINT(snumber) .AND. (.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) ) THEN
        !Remove this shell
        S(i,:) = 0.d0
        rmshells = rmshells+1
      ENDIF
    ENDDO
    !
    !If there is no shell left, de-allocate S
    IF( .NOT. ANY( NINT(S(:,4)).NE.0 ) ) THEN
      DEALLOCATE(S)
    ENDIF
    !
  ENDIF
  !
ENDIF
!
!
!
!
300 CONTINUE
CALL ATOMSK_MSG(2123,(/""/),(/DBLE(rmshells)/))
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
!
END SUBROUTINE RMSHELLS_XYZ
!
END MODULE rmshells