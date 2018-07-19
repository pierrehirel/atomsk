MODULE substitute
!
!**********************************************************************************
!*  SUBSTITUTE                                                                    *
!**********************************************************************************
!* This module reads an array of type (species x y z) and substitutes             *
!* one atomic species by another.                                                 *
!**********************************************************************************
!* (C) October 2010 - Pierre Hirel                                                *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 19 Feb. 2014                                     *
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
SUBROUTINE SUBSTITUTE_XYZ(P,S,species1,species2,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=2):: species
CHARACTER(LEN=3),INTENT(IN):: species1  !species1 will be replaced by species2
CHARACTER(LEN=3),INTENT(IN):: species2  !species1 will be replaced by species2
CHARACTER(LEN=3):: sp1, sp2   !sp1=species(1), sp2=species(2)
CHARACTER(LEN=128):: msg
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT  !mask for atom list
INTEGER:: i, Nsub
REAL(dp):: testreal
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of cores, shells
!
!Initialize variables
Nsub = 0
sp1 = species1
sp2 = species2
!
WRITE(msg,*) 'Entering SUBSTITUTE_XYZ: ', sp1, sp2
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!Check if sp1 was entered as a species or atom number
READ(sp1,*,ERR=50) testreal
!No error? Then it was a number => convert it to atom species
CALL ATOMSPECIES(testreal,species)
sp1 = species
!Else if there was an error then sp1 was an atom species, we don't change it
50 CONTINUE
!Check if sp2 was entered as a species or atom number
READ(sp2,*,ERR=51) testreal
!No error? Then it was a number => convert it to atom species
CALL ATOMSPECIES(testreal,species)
sp2 = species
!Else if there was an error then sp2 was an atom species, we don't change it
51 CONTINUE
!
!
CALL ATOMSK_MSG(2089,(/sp1,sp2/),(/0.d0/))
!
IF( sp1==sp2 ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2737,(/sp1,sp2/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!Check that the species exist
CALL ATOMNUMBER(sp1,testreal)
IF(testreal==0.d0) THEN
  CALL ATOMSK_MSG(801,(/sp1/),(/0.d0/))
  GOTO 1000
ENDIF
CALL ATOMNUMBER(sp2,testreal)
IF(testreal==0.d0) THEN
  CALL ATOMSK_MSG(801,(/sp2/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
100 CONTINUE
!Substitute the atoms
DO i=1,SIZE(P,1)
  IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
    CALL ATOMSPECIES(P(i,4),species)
    IF(species==sp1) THEN
      CALL ATOMNUMBER(sp2,P(i,4))
      Nsub = Nsub+1
    ENDIF
    IF( ALLOCATED(S) .AND. SIZE(S,1)>0 ) THEN
      CALL ATOMSPECIES(S(i,4),species)
      IF(species==sp1) THEN
        CALL ATOMNUMBER(sp2,S(i,4))
      ENDIF
    ENDIF
  ENDIF
ENDDO
!
CALL ATOMSK_MSG(2090,(/''/),(/DBLE(Nsub)/))
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE SUBSTITUTE_XYZ
!
!
!
END MODULE substitute