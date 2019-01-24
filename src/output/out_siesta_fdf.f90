MODULE out_siesta_fdf
!
!
!**********************************************************************************
!*  OUT_FDF                                                                       *
!**********************************************************************************
!* This module writes SIESTA FDF format.                                          *
!* The FDF format is described in the SIESTA manual:                              *
!*    https://departments.icmab.es/leem/siesta/Documentation/Manuals/manuals.html *
!**********************************************************************************
!* (C) January 2019 - Pierre Hirel                                                *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 09 Jan. 2019                                     *
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
USE functions
USE messages
USE files
USE subroutines
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE WRITE_FDF(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
INTEGER:: i, j, Nspecies
INTEGER:: Ntypes    !number of atom types
INTEGER:: typecol !index of atom types in AUX
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: atypes
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
typecol = 0
Nspecies = 0
!Find how many different species are in P
CALL FIND_NSP(P(:,4),atypes)
Ntypes = SIZE(atypes(:,1))
!
!
msg = 'entering WRITE_FDF'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Check if some auxiliary properties are present
IF( ALLOCATED(AUXNAMES) ) THEN
  DO i=1,SIZE(AUXNAMES)
    IF( TRIM(ADJUSTL(AUXNAMES(i)))=='type' ) THEN
      typecol = i
      !Find how many different atom types are in AUX
      CALL FIND_NSP(AUX(:,typecol),atypes)
      Ntypes = SIZE(atypes(:,1))
    ENDIF
  ENDDO
ENDIF
!
!
!
OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN',ERR=500)
!
!Write header of FDF file
IF( ALLOCATED(comment) .AND. SIZE(comment)>=1 ) THEN
  temp = TRIM(ADJUSTL(comment(1)))
  IF( temp(1:1)=="#" ) THEN
    temp = ADJUSTL(temp(2:))
  ENDIF
  WRITE(40,*) "SystemName          "//TRIM(ADJUSTL(temp))
ELSE
  WRITE(40,*) "SystemName          Unknown atomic system"
ENDIF
!Determine formula of the compound
msg=''
DO i=1,SIZE(atypes,1)
  CALL ATOMSPECIES(atypes(i,1),species)
  temp=''
  IF( NINT(atypes(i,2))>1 ) THEN
    WRITE(temp,'(i10)') NINT(atypes(i,2))
    temp = TRIM(species)//TRIM(ADJUSTL(temp))
  ELSE
    temp = TRIM(species)
  ENDIF
  msg = TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(temp))
ENDDO
WRITE(40,*) "SystemLabel         "//TRIM(ADJUSTL(msg))
WRITE(40,*) ""
!
!Write number of atoms and species
WRITE(temp,*) SIZE(P,1)
WRITE(40,*) "NumberOfAtoms       "//TRIM(ADJUSTL(temp))
WRITE(temp,*) Ntypes
WRITE(40,*) "NumberOfSpecies     "//TRIM(ADJUSTL(temp))
!
!Write cell vectors
WRITE(40,*) ""
WRITE(40,*) "%block LatticeVectors"
DO i=1,3
  WRITE(40,'(3X,3(f16.9,3X))') H(i,1), H(i,2), H(i,3)
ENDDO
WRITE(40,*) "%endblock LatticeVectors"
!
!
!Write atom types (or "index"), atomic numbers, and chemical symbols
WRITE(40,*) ""
WRITE(40,*) "%block ChemicalSpeciesLabel"
Nspecies=0
DO i=1,SIZE(atypes,1)
  CALL ATOMSPECIES(atypes(i,1),species)
  IF( typecol>0 ) THEN
    !If atom types are in auxiliary properties, use it
    Nspecies = NINT( AUX(i,typecol) )
  ELSE
    !Otherwise replace species by atom types in their order of appearance
    Nspecies=Nspecies+1
  ENDIF
  WRITE(40,'(i3,2X,i3,2X,a2)') Nspecies, NINT(atypes(i,1)), species
ENDDO
WRITE(40,*) "%endblock ChemicalSpeciesLabel"
!
!
!Write atom positions
WRITE(40,*) ""
WRITE(40,*) "AtomicCoordinatesFormat Ang"
WRITE(40,*) "%block AtomicCoordinatesAndAtomicSpecies"
DO i=1,SIZE(P,1)
  IF( typecol>0 ) THEN
    !If atom types are in auxiliary properties, use it
    Nspecies = NINT( AUX(i,typecol) )
  ELSE
    !Otherwise replace species by atom types in their order of appearance
    DO j=1,SIZE(atypes(:,1))
      IF( atypes(j,1)==INT(P(i,4)) ) Nspecies = j
    ENDDO
  ENDIF
  WRITE(40,'(3(f16.9,2X),i3)') P(i,1), P(i,2), P(i,3), Nspecies
ENDDO
WRITE(40,*) "%endblock AtomicCoordinatesAndAtomicSpecies"
!
!
!
!
500 CONTINUE
CLOSE(40)
msg = "FDF"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
1000 CONTINUE
!
END SUBROUTINE WRITE_FDF
!
END MODULE out_siesta_fdf
