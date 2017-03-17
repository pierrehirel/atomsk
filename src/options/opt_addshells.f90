MODULE addshells
!
!**********************************************************************************
!*  ADDSHELLS                                                                     *
!**********************************************************************************
!* This module creates shells (in the sens of a core-shell model) for some        *
!* or all atoms.                                                                  *
!**********************************************************************************
!* (C) May 2010 - Pierre Hirel                                                    *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 12 June 2014                                     *
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
SUBROUTINE ADDSHELLS_XYZ(P,S,species)
!
!
IMPLICIT NONE
CHARACTER(LEN=2),INTENT(IN):: species
CHARACTER(LEN=128):: temp, msg
INTEGER:: i, Ssize
INTEGER:: Nshells
REAL(dp):: snumber
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: S
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries
!
!Initialize variables
Nshells=0
Ssize = 0
snumber = 0.d0
!
!
msg = 'Entering CSHELLS_XYZ'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
IF(species=='XX') THEN
  temp = 'all'
ELSE
  temp = species
ENDIF
CALL ATOMSK_MSG(2054,(/TRIM(ADJUSTL(temp))/),(/0.d0/))
!
CALL ATOMNUMBER(species,snumber)
IF(species.NE.'XX' .AND. snumber==0.d0) THEN
  nerr=nerr+1
  CALL ATOMSK_MSG(801,(/species/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
100 CONTINUE
!Check if these atoms already have shells
IF( species.NE.'XX' .AND. ALLOCATED(S) ) THEN
  DO i=1,SIZE(S(:,1))
    IF( NINT(S(i,4))==NINT(snumber) ) THEN
      nwarn=nwarn+1
      CALL ATOMSK_MSG(2722,(/species/),(/0.d0/))
      GOTO 1000
    ENDIF
  ENDDO
ENDIF
!
!Find how many atoms belong to the 'species' in P
IF(species=='XX') THEN
  Ssize = SIZE(P(:,1))
  IF(ALLOCATED(S)) DEALLOCATE(S)
ELSE
  CALL FIND_NSP(P(:,4),aentries)
  DO i=1,SIZE(aentries(:,1))
    IF( NINT(aentries(i,1))==NINT(snumber) ) THEN
      Ssize = NINT(aentries(i,2))
    ENDIF
  ENDDO
ENDIF
!
!Allocate S for the shells
IF( Ssize==0 ) THEN
  !No such atom in the system => abort
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2723,(/species/),(/0.d0/))
  GOTO 1000
ELSE
  !
  IF( .NOT.ALLOCATED(S) ) THEN
    !Create the array S for storing shells
    ALLOCATE( S(SIZE(P,1),4) )
    S(:,:)=0.d0
    Ssize=0
  ENDIF
  !
  !Write shell positions in S
  DO i=1,SIZE(P(:,1))
    IF(species=='XX' .OR. NINT(P(i,4))==NINT(snumber) ) THEN
      S(i,1:4) = P(i,1:4)
      Nshells=Nshells+1
    ENDIF
    !Uncomment the following line to shift shells
    !S(Ssize,3)=S(Ssize,3)+0.05d0
  ENDDO
ENDIF
!
!
!
400 CONTINUE
CALL ATOMSK_MSG(2055,(/species/),(/DBLE(Nshells)/))
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
END SUBROUTINE ADDSHELLS_XYZ
!
!
!
END MODULE addshells