MODULE in_vasp_poscar
!
!
!**********************************************************************************
!*  IN_VASP_POSCAR                                                                *
!**********************************************************************************
!* This module reads POSCAR files, as well as CONTCAR files since it              *
!* is the same format. These files are used by VASP.                              *
!* The POSCAR format is described here:                                           *
!*    http://cms.mpi.univie.ac.at/vasp/vasp/POSCAR_file.html                      *
!* Note that in the POSCAR format used by VASP 5.x and greater the atom species   *
!* appear in the file, however in the format used by VASP 4.x and before the      *
!* atom species do not appear in the POSCAR file, only in the associated          *
!* POTCAR file. So, in the latter case this module will check for a POTCAR file   *
!* in the current directory, and if it exists will read atom species from it.     *
!* If no POTCAR file can be found, dummy atom numbers (1, 2, 3...) will be used.  *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 19 March 2014                                    *
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
SUBROUTINE READ_POSCAR(inputfile,H,P,comment,AUXNAMES,AUX)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=1):: fixx, fixy, fixz
CHARACTER(LEN=2):: species
CHARACTER(LEN=4):: coord
CHARACTER(LEN=128):: msg, test
CHARACTER(LEN=2),DIMENSION(20):: speciestxt
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: fileexists
LOGICAL:: foundsp !have we already found species in POSCAR file?
LOGICAL:: seldyn  !use "Selective dynamics"?
INTEGER:: i, j, k
INTEGER:: NP
INTEGER,DIMENSION(20):: NPi, spi
REAL(dp):: a0, tempreal
REAL(dp):: P1, P2, P3
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
foundsp = .FALSE.
seldyn = .FALSE.
i = 0
NP = 0
NPi(:) = 0
spi(:) = 0
a0 = 1.d0
 coord = 'cart'
ALLOCATE(comment(1))
!
!
100 CONTINUE
WRITE(msg,*) 'entering READ_POSCAR'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
READ(30,'(a128)',ERR=400,END=400) comment(1)
READ(30,*,ERR=410,END=410) a0
READ(30,*,ERR=410,END=410) H(1,1), H(1,2), H(1,3)
READ(30,*,ERR=410,END=410) H(2,1), H(2,2), H(2,3)
READ(30,*,ERR=410,END=410) H(3,1), H(3,2), H(3,3)
H(:,:) = a0*H(:,:)
!
!VASP 4.x: next line contains integers = number of atoms for each type
!VASP 5.x: next line contains atom species, and the line after it contains the integers
READ(30,'(a128)',ERR=420,END=420) test
!Try to read an atom species. If not successful, directly attempt to read integers
READ(test,*,ERR=105,END=105) species
CALL ATOMNUMBER(species,tempreal)
IF( NINT(tempreal)==0 ) THEN
  !It was not an atom species => attempt to read integers
  GOTO 105
ELSE
  !It was an atom species => read all species from current line
  DO k=1,20
    READ(test,*,ERR=104,END=104) (speciestxt(j), j=1,k)
  ENDDO
  104 CONTINUE
  DO j=1,k
    CALL ATOMNUMBER(speciestxt(j),tempreal)
    spi(j) = NINT(tempreal)
  ENDDO
  foundsp=.TRUE.
  !Read next line, it should contain the integers
  READ(30,'(a128)',ERR=420,END=420) test
ENDIF
105 CONTINUE
!Read the number of atoms for each species
DO k=1,20
  READ(test,*,ERR=110,END=110) (NPi(j), j=1,k)
ENDDO
!
!
110 CONTINUE
k=k-1
WRITE(msg,*) "Number of different species: ", k
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,*) "NP for each species: ", (NPi(j), j=1,k)
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Check if the next line contains "Selective dynamics"
READ(30,'(a128)',ERR=410,END=410) test
test = TRIM(ADJUSTL(test))
IF( test(1:1)=="S" .OR. test(1:1)=="s" ) THEN
  seldyn = .TRUE.
  WRITE(msg,*) "Selective Dynamics detected"
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ELSE
  BACKSPACE(30)
ENDIF
!
IF( .NOT. foundsp ) THEN
  !If atom species were not found in the POSCAR file, check for a "POTCAR" file
  !If it exists, attempt to read the atom species from it
  INQUIRE(FILE='POTCAR',EXIST=fileexists)
  IF(fileexists) THEN
    CALL ATOMSK_MSG(1003,(/''/),(/0.d0/))
    OPEN(UNIT=31,FILE='POTCAR',FORM='FORMATTED',STATUS='OLD')
    i=0
    DO
      READ(31,'(a128)',ERR=120,END=120) test
      test = ADJUSTL(test)
      IF(test(1:2)=='US') THEN
        i=i+1
        READ(test(3:),*) species
        CALL ATOMNUMBER(species,tempreal)
        spi(i) = NINT(tempreal)
      ENDIF
    ENDDO
  ENDIF
ENDIF
!
120 CONTINUE
CLOSE(31)
!Next line should contain a keyword deciding if the atom positions are in Direct or Cartesian coordinates
!According to VASP documentation, anything starting with C, c, K or k switches to cartesian, anything else means Direct
test=""
DO WHILE( LEN_TRIM(test)==0 .OR. test(1:1)=="#" )
  READ(30,'(a128)',ERR=410,END=410) test
  test = TRIM(ADJUSTL(test))
ENDDO
IF( test(1:1)=="C" .OR. test(1:1)=="c" .OR. test(1:1)=="K" .OR. test(1:1)=="k" ) THEN
  coord = 'cart'
  WRITE(msg,*) "Found cartesian coordinates"
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ELSE
  coord = 'frac'
  WRITE(msg,*) "Found fractional coordinates"
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ENDIF
!
!
!
200 CONTINUE
NP = 0
DO j=1,k
  NP = NP+NPi(j)
ENDDO
ALLOCATE(P(NP,4))
P(:,:) = 0.d0
IF(seldyn) THEN
  ALLOCATE(AUXNAMES(3))
  AUXNAMES(1) = "fixx"
  AUXNAMES(2) = "fixy"
  AUXNAMES(3) = "fixz"
  ALLOCATE(AUX(NP,3))
  AUX(:,:) = 0.d0
ENDIF
!
!Read atomic coordinates
WRITE(msg,*) "Reading coordinates"
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
NP = 0
DO j=1,k
  DO i=1,NPi(j)
    IF(spi(1)==0) THEN
      !If no species read from POTCAR, use j
      P(NP+i,4) = j
    ELSE
      !Else read the species read from POTCAR
      P(NP+i,4) = spi(j)
    ENDIF
    READ(30,'(a128)',END=400,ERR=400) test
    READ(test,*,END=400,ERR=400) P1, P2, P3
    P(NP+i,1) = P1
    P(NP+i,2) = P2
    P(NP+i,3) = P3
    !
    IF(seldyn) THEN
      READ(test,*,END=400,ERR=400) P1, P2, P3, fixx, fixy, fixz
      IF(fixx=="F") AUX(NP+i,1) = 1.d0
      IF(fixy=="F") AUX(NP+i,2) = 1.d0
      IF(fixz=="F") AUX(NP+i,3) = 1.d0
    ENDIF
  ENDDO
  NP = NP+NPi(j)
ENDDO
GOTO 500
!
!
!
400 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(NP+i)/))
nerr = nerr+1
GOTO 1000
!
410 CONTINUE
CALL ATOMSK_MSG(1809,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
420 CONTINUE
CALL ATOMSK_MSG(1810,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
500 CONTINUE
!In case of fractional coordinates, convert them to cartesian
IF(coord=='frac') THEN
  CALL FRAC2CART(P,H)
ENDIF
CLOSE(30)
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE READ_POSCAR
!
END MODULE in_vasp_poscar
