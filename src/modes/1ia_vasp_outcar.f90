MODULE oia_vaspout
!
!**********************************************************************************
!*  1IA_VASPOUT                                                                   *
!**********************************************************************************
!* This module reads VASP output files (OUTCAR), containing                       *
!* several snapshots, and writes each snapshot to a separate file.                *
!* The output produced by VASP is described here:                                 *
!*    https://www.vasp.at/wiki/index.php/OUTCAR                                   *
!**********************************************************************************
!* (C) May 2021 - Pierre Hirel                                                    *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 21 July 2022                                     *
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
!
USE atoms
USE comv
USE constants
USE messages
USE files
USE subroutines
USE options
USE writeout
USE qepw_ibrav
!
!
CONTAINS
!
SUBROUTINE ONEINALL_OUTCAR(inputfile,out_prefix,outfileformats,options_array)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to output
CHARACTER(LEN=128):: msg, test
CHARACTER(LEN=4096):: out_prefix, outputfile, outfileformat
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: fileexists
LOGICAL:: isreduced
LOGICAL:: readforces
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: NtypesPOTCAR, NtypesNIONS  !number of different atom types
INTEGER:: i, j, NP, Nline, Nsnap, snap, strlen
INTEGER:: Nsys  !number of systems converted
INTEGER,DIMENSION(10,2):: atypes !atom "types" (atomic numbers) and their number
REAL(dp):: alat  !lattice constant
REAL(dp):: snumber !atomic number
REAL(dp),DIMENSION(3,3):: Huc   !Box vectors of unit cell (unknown, set to 0 here)
REAL(dp),DIMENSION(3,3):: H, Htemp   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystal orientation
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P,S  !positions of cores, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
outputfile=''
readforces=.FALSE.
NP = 0
Nsnap = 0
Nsys = 0
snap = -1  !so that snap index will start at 0
NtypesPOTCAR = 0
NtypesNIONS = 0
atypes(:,:) = 0
alat = 0.d0
Huc(:,:) = 0.d0
C_tensor(:,:) = 0.d0
ORIENT(:,:) = 0.d0
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
IF(ALLOCATED(P)) DEALLOCATE(P)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
ALLOCATE(comment(1))
!
!
msg = 'entering ONEINALL_VASPOUT'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
!
100 CONTINUE
!Get number of atoms
DO WHILE(NP==0)
  READ(30,'(a128)',ERR=800,END=800) test
  test = TRIM(ADJUSTL(test))
  IF( INDEX(test,"NIONS")>0 ) THEN
    i = INDEX(test,"NIONS")
    test = test(i+1:)
    i = SCAN(test,"=")
    test = test(i+1:)
    READ(test,*,END=800,ERR=800) NP
    EXIT
  ENDIF
ENDDO
!
110 CONTINUE
WRITE(msg,*) 'NIONS = ', NP
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ALLOCATE(P(NP,4))
P(:,:) = 0.d0
!
!Get atom types
!First, try from lines starting with "POTCAR:", for example:
!       POTCAR:   PAW_GGA O 05Jan2001
!       POTCAR:   PAW_PBE Sn_d 06Sep2000
!Assuming that element name is in 3rd position
REWIND(30)
NtypesPOTCAR = 0
DO Nline=1,1000  !restrict reading to the first thousand lines
  READ(30,'(a128)',ERR=800,END=800) test
  test = TRIM(ADJUSTL(test))
  IF( test(1:7)=="POTCAR:" ) THEN
    DO WHILE( test(1:7)=="POTCAR:" )
      test = TRIM(ADJUSTL(test(8:)))
      i = SCAN(test," ")
      test = TRIM(ADJUSTL(test(i+1:)))
      species = test(1:2)
      !Verify that it is a valid atomic species
      CALL ATOMNUMBER(species,snumber)
      IF( snumber==0 ) THEN
        !Try with only first letter
        species(2:2) = " "
        CALL ATOMNUMBER(species,snumber)
      ENDIF
      NtypesPOTCAR = NtypesPOTCAR + 1
      atypes(NtypesPOTCAR,1) = NINT(snumber)
      !Read next line
      READ(30,'(a128)',ERR=800,END=800) test
      test = TRIM(ADJUSTL(test))
    ENDDO
  ELSEIF( test(1:6)=="IBRION" ) THEN
    !Read value of IBRION used for this simulation
    i = SCAN(test,'=')
    IF(i==0) i=6
    test = TRIM(ADJUSTL(test(i+1:)))
    READ(test,*,END=115,ERR=115) j
    IF(j<0) THEN
      !IBRION<0 means that a static calculation was performed
      !(i.e. electronic relaxation, ions are fixed)
      !Therefore the OUTCAR file does not contain any snapshot
      nerr=nerr+1
      CALL ATOMSK_MSG(4833,(/""/),(/DBLE(j)/))
      GOTO 1000
    ENDIF
  ENDIF
  115 CONTINUE
ENDDO
120 CONTINUE
REWIND(30)
!
!If atom types could not be determined, try reading POTCAR file from current folder
IF( NtypesPOTCAR==0 ) THEN
  INQUIRE(FILE='POTCAR',EXIST=fileexists)
  IF(fileexists) THEN
    CALL ATOMSK_MSG(1003,(/''/),(/0.d0/))
    OPEN(UNIT=31,FILE="POTCAR",FORM='FORMATTED',STATUS='OLD')
    DO
      READ(31,'(a128)',END=130,ERR=130) test
      test = TRIM(ADJUSTL(test))
      IF( (test(1:2)=="US" .OR. test(1:3)=="PAW") .AND. INDEX(test,"radial")==0 ) THEN
        strlen = INDEX(test," ")
        test = TRIM(ADJUSTL(test(strlen+1:)))
        species = test(1:2)
        !Verify that it is a valid atomic species
        CALL ATOMNUMBER(species,snumber)
        IF( snumber==0 ) THEN
          !Try with only first letter
          species(2:2) = " "
          CALL ATOMNUMBER(species,snumber)
        ENDIF
        NtypesPOTCAR = NtypesPOTCAR + 1
        atypes(NtypesPOTCAR,1) = NINT(snumber)
      ENDIF
    ENDDO
    CLOSE(31)
  ENDIF
ENDIF
!
130 CONTINUE
WRITE(msg,*) 'NTYPES(POTCAR) = ', NtypesPOTCAR
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
REWIND(30)
!Get the number of each atom type
DO
  READ(30,'(a128)',ERR=800,END=800) test
  test = TRIM(ADJUSTL(test))
  IF( INDEX(test,"ions per type")>0 ) THEN
    i = INDEX(test,"ions per type")
    test = test(i+1:)
    i = SCAN(test,"=")
    test = TRIM(ADJUSTL(test(i+1:)))
    i=0
    DO WHILE(LEN_TRIM(test)>0)
      NtypesNIONS = NtypesNIONS + 1
      READ(test,*,END=130,ERR=130) atypes(NtypesNIONS,2)
      !Truncate string "test"
      j = SCAN(test," ")
      test = TRIM(ADJUSTL(test(j:)))
    ENDDO
    EXIT
  ENDIF
ENDDO
!
140 CONTINUE
REWIND(30)
WRITE(msg,*) 'NTYPES(IONS) = ', NtypesNIONS
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
IF( NtypesNIONS.NE.NtypesPOTCAR ) THEN
  PRINT*, " /!\ WARNING: inconsistent number of types"
ENDIF
WRITE(msg,*) 'ATOM TYPES'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
DO i=1,MAX(NtypesNIONS,NtypesPOTCAR)
  WRITE(msg,*) atypes(i,1), atypes(i,2)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ENDDO
!
!Fill array P with the atom types that were found
NP = 0
DO j=1,NtypesNIONS
  DO i=1,atypes(j,2)
    NP = NP+1
    P(NP,4) = atypes(j,1)
  ENDDO
ENDDO
!
!Get the cell size
DO WHILE( VECLENGTH(H(1,:))<1.d-6 .OR. VECLENGTH(H(2,:))<1.d-6 .OR. VECLENGTH(H(3,:))<1.d-6 )
  READ(30,'(a128)',ERR=800,END=800) test
  test = TRIM(ADJUSTL(test))
  IF( test(1:2)=="A1" ) THEN
    !Remove parenthesis and commas from string "test"
    strlen = SCAN(test,"=")
    test = test(strlen+1:)
    CALL STR_CHAR2SPACE(test,"(),")
    READ(test,*,ERR=800,END=800) H(1,1), H(1,2), H(1,3)
  ELSEIF( test(1:2)=="A2" ) THEN
    !Remove parenthesis and commas from string "test"
    strlen = SCAN(test,"=")
    test = test(strlen+1:)
    CALL STR_CHAR2SPACE(test,"(),")
    READ(test,*,ERR=800,END=800) H(2,1), H(2,2), H(2,3)
  ELSEIF( test(1:2)=="A3" ) THEN
    !Remove parenthesis and commas from string "test"
    strlen = SCAN(test,"=")
    test = test(strlen+1:)
    CALL STR_CHAR2SPACE(test,"(),")
    READ(test,*,ERR=800,END=800) H(3,1), H(3,2), H(3,3)
  ENDIF
ENDDO
WRITE(msg,*) 'CELL'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
DO i=1,3
  WRITE(msg,*) "   ", H(i,1), H(i,2), H(i,3)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ENDDO
!
!
!
200 CONTINUE
REWIND(30)
WRITE(msg,*) 'READING SNAPSHOTS ...'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!Check that arrays are correctly allocated
IF(.NOT.ALLOCATED(P)) GOTO 800
DO  !loop on all snapshots
  !Initialize variables
  isreduced=.FALSE.
  readforces=.FALSE.
  P(:,1:3) = 0.d0
  !
  !
  READ(30,'(a128)',ERR=1000,END=1000) test
  test = TRIM(ADJUSTL(test))
  IF( test(1:22)=="direct lattice vectors" ) THEN
    !read or update cell vectors on the next 3 lines
    DO i=1,3
      READ(30,*,ERR=800,END=800) H(i,1), H(i,2), H(i,3)
    ENDDO
    !
  ELSEIF( test(1:8)=="POSITION" ) THEN
    !New snapshot: read atom positions
    snap = snap+1
    !
    !Set the output file name
    WRITE(msg,*) snap
    outputfile = TRIM(ADJUSTL(out_prefix))//"_"//TRIM(ADJUSTL(msg))
    comment(1) = "# VASP output snapshot # "//TRIM(ADJUSTL(msg))
    WRITE(msg,*) 'snap, outputfile:', snap, TRIM(ADJUSTL(outputfile))
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
    !Check if forces must be read
    IF( INDEX(test,"FORCE")>0 ) THEN
      readforces=.TRUE.
      IF( .NOT.ALLOCATED(AUX) ) THEN
        ALLOCATE(AUX(NP,3))
        ALLOCATE(AUXNAMES(3))
        AUXNAMES(1) = "fx"
        AUXNAMES(2) = "fy"
        AUXNAMES(3) = "fz"
      ENDIF
      AUX(:,1:3) = 0.d0
    ENDIF
    !
    !Ignore next line because it only contains "---"
    READ(30,'(a128)',ERR=500,END=500) test
    !
    !Read all atom positions
    DO i=1,SIZE(P,1)
      READ(30,'(a128)',ERR=500,END=500) test
      IF( readforces ) THEN
        READ(test,*,ERR=500,END=500) P(i,1), P(i,2), P(i,3), AUX(i,1), AUX(i,2), AUX(i,3)
      ELSE
        READ(test,*,ERR=500,END=500) P(i,1), P(i,2), P(i,3)
      ENDIF
    ENDDO
    !
    !Convert positions to cartesian coordinates
    IF(isreduced) THEN
      !Positions are in reduced coordinates
      CALL FRAC2CART(P,H)
    ENDIF
    CALL ATOMSK_MSG(4043,(/''/),(/0.d0/))
    !
    !
    !
    280 CONTINUE
    !H may be modified by options
    !=> pass a copy Htemp instead, to prevent H from being modified at each passage
    Htemp = H
    !
    !Apply options
    CALL OPTIONS_AFF(options_array,Huc,Htemp,P,S,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
    !
    !Output snapshot to one or several format(s)
    outfileformat=''
    CALL WRITE_AFF(outputfile,outfileformats,Htemp,P,S,comment,AUXNAMES,AUX)
    !
    !
    Nsys = Nsys+1
    CALL ATOMSK_MSG(4001,(/''/),(/0.d0/))
    !
  ENDIF
  !
ENDDO
!
!
!
500 CONTINUE
CALL ATOMSK_MSG(4042,(/''/),(/DBLE(Nsnap)/))
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(1801,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
1000 CONTINUE
CLOSE(30)
CALL ATOMSK_MSG(4045,(/''/),(/DBLE(Nsys)/))
!
END SUBROUTINE ONEINALL_OUTCAR
!
!
END MODULE oia_vaspout
