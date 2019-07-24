MODULE in_bopfox
!
!
!**********************************************************************************
!*  IN_BOPFOX                                                                     *
!**********************************************************************************
!* This module reads files used by the BOPfox program, which is described in:     *
!*   T. Hammerschmidt et al., Comput.Phys.Comm. 235 (2019) 221-233                *
!* More information from the Web site:                                            *
!*   http://bopfox.de/                                                            *
!**********************************************************************************
!* (C) November 2016 - Matous Mrovec                                              *
!*     ICAMS                                                                      *
!*     Ruhr-Universitaet Bochum, Germany                                          *
!*     matous.mrovec@icams.rub.de                                                 *
!* Last modification: P. Hirel - 16 July 2019                                     *
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
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE READ_BOPFOX(inputfile,H,P,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=1):: flagx, flagy, flagz  !flags for fixed atoms (T or F)
CHARACTER(LEN=2):: species
CHARACTER(LEN=4):: coord, ptype
CHARACTER(LEN=128):: test, test2
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
LOGICAL:: domag  !read magnetisation?
LOGICAL:: dofix  !read fix flags?
INTEGER:: i, j, NP, NPd, NPinert, NLINES
INTEGER:: magx, magy, magz, fixx, fixy, fixz
REAL(dp):: a0, testreal
REAL(dp):: len1, len2, len3
REAL(dp):: snumber !Mass of atoms
REAL(dp), DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P
REAL(dp), DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
i = 0
NP = 0
NPd = 0
NPinert = 0
dofix = .FALSE.  !by default no flags about fixed atoms
domag = .FALSE.  !by default assume that there is no magnetisation in the file
a0 = 1.d0
len1 = 1.d0
len2 = 1.d0
len3 = 1.d0
H(:,:) = 0.d0
!
!
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
DO
  READ(30,'(a128)',ERR=110,END=110) test
  i=i+1
  !
  test = TRIM(ADJUSTL(test))
  !
  IF( test(1:1).NE."#" .AND. test(1:1).NE."/" ) THEN
    !
    IF( verbosity==4 .AND. i<50 ) THEN
      msg = 'BOPFOX line:  '//TRIM(test)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDIF
    !
    ! Read alat
    IF( test(1:4)=='aLat' .OR. test(1:4)=='alat' .OR. test(1:4)=='ALAT' ) THEN
      j = SCAN(test,"=")
      IF( j<4 ) j=4
      READ(test(j+1:),*,ERR=410,END=410) a0
      !
    ! Read basis vectors
    ELSEIF( test(1:2)=='a1' .OR. test(1:2)=='A1' ) THEN
      j = SCAN(test,"=")
      IF( j<2 ) j=2
      READ(test(j+1:),*,ERR=410,END=410) H(1,1), H(1,2), H(1,3)
    ELSEIF( test(1:2)=='a2' .OR. test(1:2)=='A2' ) THEN
      j = SCAN(test,"=")
      IF( j<2 ) j=2
      READ(test(j+1:),*,ERR=410,END=410) H(2,1), H(2,2), H(2,3)
    ELSEIF( test(1:2)=='a3' .OR. test(1:2)=='A3' ) THEN
      j = SCAN(test,"=")
      IF( j<2 ) j=2
      READ(test(j+1:),*,ERR=410,END=410) H(3,1), H(3,2), H(3,3)
      !
    ! Read coord info
    ELSEIF( test(1:5)=='coord' ) THEN
      j = SCAN(test,"=")
      IF( j<5 ) j=5
      test = TRIM(ADJUSTL(test(j+1:)))
      IF (test(1:6)=="direct") THEN
        coord = 'frac'
        WRITE(msg,*) "Found direct/fractional coordinates"
        CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      ELSE
        coord = 'cart'
        WRITE(msg,*) "Probably cartesian coordinates"
        CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      ENDIF
      !Count number of atoms
      IF( NP==0 ) THEN
        DO
          READ(30,'(a128)',ERR=105,END=105) test
          test = TRIM(ADJUSTL(test))
          IF( test(1:1).NE."#" .AND. test(1:1).NE."/" ) THEN
            !Remove comments at end of line if any
            j = SCAN(test,"/")
            IF( j<=0 ) THEN
              j = SCAN(test,"#")
            ENDIF
            IF( j>0 ) THEN
              test = test(:j-1)
            ENDIF
            !Detect if there are "T" or "F" flags on the line
            !NOTE: do not try to detect an isolated "F",
            !      because it may be fluorine atom at the beginning of line
            j = INDEX(test,"T ")
            IF( j<=0 ) THEN
              j = INDEX(test," T")
            ENDIF
            IF( j<=0 ) THEN
              j = INDEX(test,"F F")
            ENDIF
            IF( j<=0 ) THEN
              j = INDEX(test,"T F")
            ENDIF
            IF( j<=0 ) THEN
              j = INDEX(test,"F T")
            ENDIF
            IF( j>0 ) THEN
              dofix = .TRUE.
            ENDIF
            !Try to interpret the first 2 characters as a chemical symbol
            species = test(1:2)
            CALL ATOMNUMBER(species,snumber)
            IF( NINT(snumber)>0 ) THEN
              !Success: increment number of atoms
              NP = NP+1
            ELSE
              !Failed: we are done counting atoms
              test = ""
              EXIT
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      !A line was read that did not start with an atom species: go back one line
      BACKSPACE(30)
      WRITE(msg,*) 'Counted number of atoms: NP = ', NP
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      !
    ELSEIF( test(1:13)=='magnetisation' .OR. test(1:13)=='magnetization' ) THEN
      j = SCAN(test,"=")
      IF( j<13 ) j=13
      READ(test(j+1:),*,ERR=410,END=410) test2
      test2 = TRIM(ADJUSTL(test2))
      IF( test2(1:4)=="true" ) THEN
        domag = .TRUE.
      ENDIF
      !
    ENDIF
    !
  ENDIF
  !
  105 CONTINUE
  !
ENDDO
!
110 CONTINUE
NLINES=i
i=0
WRITE(msg,*) 'Number of lines in input file : ', NLINES
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
WRITE(msg,*) 'Lattice parameter : ', a0
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
H(1,:) = H(1,:)*a0
H(2,:) = H(2,:)*a0
H(3,:) = H(3,:)*a0
!
IF( dofix ) THEN
  WRITE(msg,*) 'Flags for fixed atoms were detected'
ENDIF
IF( domag ) THEN
  WRITE(msg,*) 'Magnetization section was detected'
ENDIF
!
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!Allocate array for atom positions
IF(NP>0) THEN
  ALLOCATE(P(NP,4))
  P(:,:) = 0.d0
ELSE
  !If NP=0 we are in trouble => error message and exit
  CALL ATOMSK_MSG(804,(/''/),(/0.d0/))
  nerr=nerr+1
  GOTO 400
ENDIF
!
!If magnetisation or fix flags must be read, allocate array AUX
IF( domag .OR. dofix ) THEN
  j = 0
  IF( domag ) THEN
    magx = 1
    magy = 2
    magz = 3
    j=3
  ENDIF
  IF( dofix ) THEN
    fixx = j+1
    fixy = j+2
    fixz = j+3
    j=j+3
  ENDIF
  ALLOCATE(AUX(NP,j))
  AUX(:,:) = 0.d0
  ALLOCATE(AUXNAMES(j))
  j = 0
  IF( domag ) THEN
    AUXNAMES(1) = "magx"
    AUXNAMES(2) = "magy"
    AUXNAMES(3) = "magz"
    j=3
  ENDIF
  IF( dofix ) THEN
    AUXNAMES(j+1) = "fixx"
    AUXNAMES(j+2) = "fixy"
    AUXNAMES(j+3) = "fixz"
    j=j+3
  ENDIF
ENDIF
!
!Read atomic positions and save them to P
REWIND(30)
DO
  READ(30,'(a128)',END=250,ERR=250) test
  test = TRIM(ADJUSTL(test))
  !
  ! Read coord info
  IF( test(1:5)=='coord') THEN
    WRITE(msg,*) 'Reading atom coordinates...'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    i=0
    DO WHILE( i<SIZE(P,1) )
      READ(30,'(a128)',END=220,ERR=400) test
    CALL ATOMSK_MSG(999,(/TRIM(test)/),(/0.d0/))
      IF( test(1:1).NE."#" .AND. test(1:1).NE."/" ) THEN
        i=i+1
        IF( dofix ) THEN
          !Read flags: T=atom is fixed along that direction (=1); F=free to move (=0)
          !NOTE: if some flags are missing, just read the atom positions and skip the flags
          READ(test,*,END=400,ERR=400) species, P(i,1), P(i,2), P(i,3)
          READ(test,*,END=210,ERR=210) species, P(i,1), P(i,2), P(i,3), flagx, flagy, flagz
          IF(flagx=="T") AUX(i,fixx) = 1.d0
          IF(flagy=="T") AUX(i,fixy) = 1.d0
          IF(flagz=="T") AUX(i,fixz) = 1.d0
        ELSE
          READ(test,*,END=400,ERR=400) species, P(i,1), P(i,2), P(i,3)
        ENDIF
        CALL ATOMNUMBER(species,P(i,4))
      ENDIF
      210 CONTINUE
    ENDDO
    220 CONTINUE
    WRITE(msg,*) 'Finished reading atom coordinates'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    BACKSPACE(30)
    !
  ELSEIF( test(1:13)=='magnetisation' .OR. test(1:13)=='magnetization' ) THEN
    WRITE(msg,*) 'Reading magnetization...'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO i=1,SIZE(AUX,1)
      READ(30,'(a128)',END=250,ERR=250) test
      j = INDEX(test,"orbital")
      IF( j>0 ) THEN
        test = ADJUSTL(test(j+8:))
      ENDIF
      READ(test,*,END=400,ERR=400) AUX(i,magx), AUX(i,magy), AUX(i,magz)
    ENDDO
    !
  ENDIF
  240 CONTINUE
ENDDO
!
250 CONTINUE
!If positions were in reduced coordinates, convert them to Cartesian
IF( coord=="frac" ) THEN
  CALL FRAC2CART(P,H)
ENDIF
GOTO 500
!
!
!
400 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
GOTO 500
!
410 CONTINUE
CALL ATOMSK_MSG(1809,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 500
!
420 CONTINUE
CALL ATOMSK_MSG(1810,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 500
!
!
!
500 CONTINUE
CLOSE(30)
!
END SUBROUTINE READ_BOPFOX
!
END MODULE in_bopfox

