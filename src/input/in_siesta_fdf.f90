MODULE in_siesta_fdf
!
!
!**********************************************************************************
!*  IN_FDF                                                                        *
!**********************************************************************************
!* This module reads SIESTA FDF format.                                           *
!* The FDF format is described in the SIESTA manual:                              *
!*    https://departments.icmab.es/leem/siesta/Documentation/Manuals/manuals.html *
!* This module tries to be flexible and to follow SIESTA's recommendations:       *
!* - the file format is case-insensitive                                          *
!* - keywords may or may not contain hyphens (-) or underscores (_), i.e.         *
!*   "LatticeConstant", "Lattice-Constant", "Lattice_Constant" are the same       *
!* In addition, if the keyword "NumberOfAtoms" is missing, this module tries to   *
!* count them.                                                                    *
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
SUBROUTINE READ_FDF(inputfile,H,P,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=16):: unit
CHARACTER(LEN=32):: atPosType !type of atomic positions: Cartesian, ScaledCartesian, etc.
CHARACTER(LEN=4096):: line, msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: comment
INTEGER:: i, j, Nline
INTEGER:: NP, Ntypes    !number of particles, of atom types
INTEGER:: strlength
INTEGER,DIMENSION(20,2):: atypes  !atom types and chemical species
REAL(dp):: alat   !lattice constant
REAL(dp),DIMENSION(3,3),INTENT(OUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: AUX !auxiliary properties
!
!
!Initialize variables
atPosType=""
alat = 1.d0
H(:,:) = 0.d0
IF(ALLOCATED(comment)) DEALLOCATE(comment)
IF(ALLOCATED(P)) DEALLOCATE(P)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
ALLOCATE(comment(2))
 comment(:) = ""
!
!
msg = 'entering READ_FDF'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
!
OPEN(UNIT=30,FILE=inputfile,STATUS='UNKNOWN',ERR=500)
!
REWIND(30)
!
Nline=0  ! init line counter
DO
  Nline = Nline+1
  READ(30,'(a4096)',ERR=200,END=200) line
  line = ADJUSTL(line)
  !
  IF( LEN_TRIM(line)>0 .AND. line(1:1).NE."#" ) THEN
    !
    !Remove dashes (-) and underscores (_),
    !everything that is after a pound character (#)
    strlength = SCAN(line,'-')
    DO WHILE(strlength>0)
      line = TRIM(line(:strlength-1))//ADJUSTL(line(strlength+1:))
      strlength = SCAN(line,'-')
    ENDDO
    strlength = SCAN(line,'_')
    DO WHILE(strlength>0)
      line = TRIM(line(:strlength-1))//ADJUSTL(line(strlength+1:))
      strlength = SCAN(line,'_')
    ENDDO
    strlength = SCAN(line,'#')
    IF( strlength>0 ) THEN
      line = line(1:strlength-1)
    ENDIF
    !
    !Save first keyword as a lowercase string
    READ(line,*) temp
    temp = TRIM(ADJUSTL( StrDnCase(temp) ))
    !
    IF( temp(1:10)=="systemname" ) THEN
      !Save system name as a comment
      temp = ADJUSTL(line(11:))
      IF( temp(1:1)=="#" ) temp = ADJUSTL(temp(2:))
      comment(1) = "# "//TRIM(ADJUSTL(temp))
      !
    ELSEIF( temp(1:11)=="systemlabel" ) THEN
      !Save system label as a comment
      temp = ADJUSTL(line(12:))
      IF( temp(1:1)=="#" ) temp = ADJUSTL(temp(2:))
      comment(2) = "# "//TRIM(ADJUSTL(temp))
      !
    ELSEIF( temp(1:13)=="numberofatoms" ) THEN
      !Read number of atoms
      IF( NP==0 ) READ(line(14:),*,ERR=190,END=190) NP
      !
    ELSEIF( temp(1:15)=="latticeconstant" ) THEN
      !Read lattice constant
      READ(line(16:),*,ERR=190,END=190) alat, unit
      IF( alat<=1.d-12 ) alat = 1.d0
      !
    ELSEIF( temp(1:23)=="atomiccoordinatesformat" ) THEN
      !Read type of atomic coordinates
      READ(line(24:),*,ERR=190,END=190) atPosType
      !
    ELSEIF( temp(1:6)=="%block" ) THEN
      !
      IF( INDEX(line,"LatticeVectors") > 0 ) THEN
        DO i=1,3
          READ(30,*,ERR=190,END=190) H(i,1), H(i,2), H(i,3)
        ENDDO
        !
      ELSEIF( INDEX(line,"ChemicalSpeciesLabel") > 0 ) THEN
        Ntypes=0
        DO WHILE( line(1:9) .NE. "%endblock" )
          READ(30,'(a4096)',ERR=190,END=190) line
          line = TRIM(ADJUSTL(line))
          IF( LEN_TRIM(line)>0 ) THEN
            Ntypes=Ntypes+1
            READ(line,*,ERR=190,END=190) atypes(Ntypes,1), atypes(Ntypes,2)
          ENDIF
        ENDDO
        !If chemical species were not attributed before, do it now
        IF( NP>0 .AND. ALLOCATED(P) ) THEN
          DO i=1,SIZE(P,1)  !loop on all atoms
            IF( NINT(P(i,4)) == 0 ) THEN
              !The species of this atom was not defined
              !Check if it has a type
              IF( ALLOCATED(AUX) .AND. NINT(AUX(i,1))>0 ) THEN
                !Use this type to find the corresponding species
                DO j=1,Ntypes
                  IF( atypes(j,1) == NINT(AUX(NP,1)) ) THEN
                    P(NP,4) = DBLE(atypes(j,2))
                    EXIT
                  ENDIF
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        !
      ELSEIF( INDEX(line,"AtomicCoordinatesAndAtomicSpecies") > 0 ) THEN
        !Verify that such a section was not already read
        IF( .NOT. ALLOCATED(P) ) THEN
          !Read atom coordinates
          IF( NP>0 ) THEN
            !Allocate array to store atom coordinates
            ALLOCATE( P(NP,4) )
            P(:,:) = 0.d0
            !Allocate array to store atom type as auxiliary property
            ALLOCATE( AUX(NP,1) )
            AUX(:,:) = 0.d0
            NP = 0
            READ(30,'(a4096)',ERR=190,END=190) line
            DO WHILE( line(1:9) .NE. "%endblock" )
              line = TRIM(ADJUSTL(line))
              IF( LEN_TRIM(line)>0 ) THEN
                NP = NP+1
                READ(line,*,ERR=190,END=190) P(NP,1), P(NP,2), P(NP,3), j
                !Save atom type as auxiliary property
                AUX(NP,1) = j
                !Determine chemical species
                P(NP,4) = 0.d0
                IF( Ntypes>0 ) THEN
                  !A section ChemicalSpeciesLabel was read before
                  DO j=1,Ntypes
                    IF( atypes(j,1) == NINT(AUX(NP,1)) ) THEN
                      P(NP,4) = DBLE(atypes(j,2))
                      EXIT
                    ENDIF
                  ENDDO
                  IF( P(NP,4)<=0.d0 ) THEN
                    !This atom type was not declared in the section ChemicalSpeciesLabel
                  ENDIF
                ELSE
                  !A section ChemicalSpeciesLabel was not read before
                ENDIF
              ENDIF
              READ(30,'(a4096)',ERR=190,END=190) line
            ENDDO
            !
          ELSE
            !Number of atoms is unknown
            !Count lines before the end of this section
            NP = 0
            READ(30,'(a4096)',ERR=190,END=190) line
            DO WHILE( line(1:9) .NE. "%endblock" )
              line = TRIM(ADJUSTL(line))
              IF( LEN_TRIM(line)>0 ) THEN
                NP = NP+1
              ENDIF
              READ(30,'(a4096)',ERR=190,END=190) line
            ENDDO
            REWIND(30)
            GOTO 190
          ENDIF !end if NP>0
          !
        ENDIF !end if not allocated(P)
        !
      ENDIF
      !
    ENDIF
    !
  ENDIF
  !
  190 CONTINUE
  !
ENDDO
!
!
!
200 CONTINUE
!
!
!
500 CONTINUE
CLOSE(30)
ALLOCATE(AUXNAMES(1))
AUXNAMES(1) = "type"
!Rescale box vectors to account for lattice constant
IF( alat.NE.1.d0 ) THEN
  H(:,:) = alat*H(:,:)
  !Also rescale atom positions if needed
  IF( atPosType=="ScaledCartesian" ) THEN
    !Atom positions were given in units of the lattice constant
    !=> multiply them all by the lattice constant
    P(:,1:3) = alat*P(:,1:3)
  ELSEIF( atPosType=="Fractional" .OR. atPosType=="fractional" .OR. atPosType=="ScaledByLatticeVectors" ) THEN
    !Atom positions were given referred to the lattice vectors
    !=> convert them to Cartesian
    CALL FRAC2CART(P,H)
  ENDIF
ENDIF
!
!
!
1000 CONTINUE
!
END SUBROUTINE READ_FDF
!
END MODULE in_siesta_fdf
