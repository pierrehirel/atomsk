MODULE in_csv 
!
!**********************************************************************************
!*  IN_CSV                                                                        *
!**********************************************************************************
!* This module reads a file in Comma-Separated Values (CSV) format.               *
!* The CSV format has no standard specification. A description is given here:     *
!*    https://tools.ietf.org/html/rfc4180                                         *
!* If the first line contains keywords like species, x, y, z, H, etc.             *
!* (as written by the module "out_csv.f90" for instance), then they are used      *
!* to determine what information is in which column.                              *
!* If the first line does not contain such information, then this module will     *
!* attempt to find real values and consider that they are the atoms positions.    *
!**********************************************************************************
!* (C) Oct. 2018 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 26 Oct. 2018                                     *
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
SUBROUTINE READ_CSV(inputfile,H,P,S,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: fields  !to store values
CHARACTER(LEN=4096):: line, msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: comment
INTEGER:: i, j, k
INTEGER:: q1, q2
INTEGER:: Nlines, Ncol  !Number of lines and columns
INTEGER:: commentcol     !column number for comments
INTEGER:: Hcol, H11, H12, H13, H21, H22, H23, H31, H32, H33 !column numbers for box vectors H
INTEGER:: atcol, scol, xcol, ycol, zcol !column numbers for atom type, species, x,y,z
INTEGER:: Sscol, Sxcol, Sycol, Szcol    !column numbers for shells: species, x, y, z
REAL(dp):: P1, P2, P3
REAL(dp):: tempreal
REAL(dp),DIMENSION(3,3),INTENT(OUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: S
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: AUX !auxiliary properties
!
!
!Initialize variables
Nlines=0
Ncol=0
atcol = 0
Hcol = 0
H11 = 0
H12 = 0
H13 = 0
H21 = 0
H22 = 0
H23 = 0
H31 = 0
H32 = 0
H33 = 0
scol=0
xcol=0
ycol=0
zcol=0
Sscol=0
Sxcol=0
Sycol=0
Szcol=0
 commentcol=0
H(:,:) = 0.d0
!
!
100 CONTINUE
msg = 'entering READ_CSV'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
OPEN(UNIT=30,FILE=inputfile,STATUS='UNKNOWN',ERR=800)
!
! Read first line and get number of fields
! The first line may contain the description of the data (or not)
! Here, if Atomsk recognizes some keywords (such as x, y, z), it will
! use it to store the corresponding data
READ(30,'(a4096)',END=800,ERR=800) line
line = ADJUSTL(line)
Ncol=0
DO WHILE(LEN_TRIM(line)>0)
  Ncol = Ncol+1
  !Check for quotes
  IF( line(1:1)=='"' ) THEN
    !This field starts with a quote => find the next quote
    line = ADJUSTL(line(2:))
    j = SCAN(line,'"')
    !Save field in temp. string
    temp = line(1:j-1)
    !Remove quote
    line(j:j) = " "
  ELSE
    !There are no quotes => find the next comma
    j = SCAN(line,',')
    IF( j==0 ) j=LEN_TRIM(line)+1
    !Save field in temp. string
    temp = TRIM(ADJUSTL(line(1:j-1)))
  ENDIF
  !Check if keywords can be recognized (species, x, y, z...)
  IF( temp(1:7)=="species" .OR. temp(1:7)=="SPECIES" ) THEN
    scol = Ncol
  ELSEIF( scol==0 .AND. ( temp(1:4)=="atom" .OR. temp(1:4)=="ATOM" .OR. &
        & temp(1:4)=="type" .OR. temp(1:4)=="TYPE" ) ) THEN
    scol = Ncol
  ELSEIF ( temp(1:1)=="x" .OR. temp(1:1)=="X" ) THEN
    xcol = Ncol
  ELSEIF ( temp(1:1)=="y" .OR. temp(1:1)=="Y" ) THEN
    ycol = Ncol
  ELSEIF ( temp(1:1)=="z" .OR. temp(1:1)=="Z" ) THEN
    zcol = Ncol
  ELSEIF( temp(1:6)=="shells" .OR. temp(1:6)=="SHELLS" .OR. temp=="S" ) THEN
    Sscol = Ncol
  ELSEIF ( temp(1:2)=="sx" .OR. temp(1:2)=="SX" ) THEN
    Sxcol = Ncol
  ELSEIF ( temp(1:2)=="sy" .OR. temp(1:2)=="SY" ) THEN
    Sycol = Ncol
  ELSEIF ( temp(1:2)=="sz" .OR. temp(1:2)=="SZ" ) THEN
    Szcol = Ncol
  ELSEIF ( temp(1:3)=="H11" ) THEN
    !Box vector components will be in different columns
    H11 = Ncol
  ELSEIF ( temp(1:3)=="H12" ) THEN
    H12 = Ncol
  ELSEIF ( temp(1:3)=="H13" ) THEN
    H13 = Ncol
  ELSEIF ( temp(1:3)=="H21" ) THEN
    H21 = Ncol
  ELSEIF ( temp(1:3)=="H22" ) THEN
    H22 = Ncol
  ELSEIF ( temp(1:3)=="H23" ) THEN
    H23 = Ncol
  ELSEIF ( temp(1:3)=="H31" ) THEN
    H31 = Ncol
  ELSEIF ( temp(1:3)=="H32" ) THEN
    H32 = Ncol
  ELSEIF ( temp(1:3)=="H33" ) THEN
    H33 = Ncol
  ELSEIF ( temp(1:1)=="H" ) THEN
    !All components of box vectors are in this column
    Hcol = Ncol
  ENDIF
  !Erase this field and the comma from the string
  line = ADJUSTL(line(j+1:))
  !If last field is empty, count it anyway
  !IF( LEN_TRIM(line)==0 ) Ncol=Ncol+1
ENDDO
Nlines = Nlines+1
IF( verbosity>=4 ) THEN
  WRITE(msg,*) "Counted ", Ncol, " fields on first line"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "  scol, xcol, ycol, zcol = ", scol, xcol, ycol, zcol
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  IF( Sscol>0 ) THEN
    WRITE(msg,*) "  SHELLS:    s, x, y, z  = ", Sscol, Sxcol, Sycol, Szcol
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
  WRITE(msg,*) "  H = ", Hcol
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!
!
150 CONTINUE
!Count number of lines
DO
  READ(30,'(a4096)',ERR=200,END=200) line
  line = ADJUSTL(line)
  Nlines = Nlines+1
  !Count number of fields in this line
  k=0
  DO WHILE(LEN_TRIM(line)>0)
    k = k+1
    IF( line(1:1)=="," ) line = ADJUSTL(line(2:))
    !Check for quotes
    IF( line(1:1)=='"' ) THEN
      !This field starts with a quote
      !Remove the first quote sign
      line = ADJUSTL(line(2:))
      !Find the next quote
      j = SCAN(line,'"')
      IF( j>1 ) THEN
        temp = TRIM(ADJUSTL(line(1:j-1)))
      ENDIF
      !Find the comma after the quote
      line = ADJUSTL(line(j+1:))
      j = SCAN(line,',')
      IF( j>0 ) THEN
        line = ADJUSTL(line(j:))
      ELSE
        !No more commas
        line = ""
      ENDIF
    ELSE
      !There are no quotes
      !Find the next comma
      j = SCAN(line,',')
      IF( j>1 ) THEN
        temp = TRIM(ADJUSTL(line(1:j-1)))
        line = ADJUSTL(line(j:))
      ELSEIF( j==1 ) THEN
        !This is an empty field
        temp = ""
        line = ADJUSTL(line(j:))
      ELSE
        !No more commas =>
        temp = TRIM(line)
        line = ""
      ENDIF
    ENDIF
    !PRINT*, "field(", k, ") = ", TRIM(temp)
    !
    !If the column for atom species was not identified when reading the first line,
    !try to convert the current text into an atom mass
    IF( scol==0 ) THEN
      IF( LEN_TRIM(temp)<=2 ) THEN
        CALL ATOMMASS(TRIM(temp),tempreal)
        IF( tempreal>=0.99d0 ) THEN
          !An atom species was recognized => assume that atom species are stored in this column
          scol = k
        ENDIF
      ENDIF
    ENDIF
    !
    !If the columns for x,y,z were not identified when reading the first line,
    !try to find columns containing real numbers
    IF( xcol==0 .OR. ycol==0 .OR. zcol==0 ) THEN
      READ(temp,*,END=180,ERR=180) tempreal
      !If we're still here, a real number was successfully read
      IF( xcol==0 ) THEN
        xcol = k
      ELSEIF( ycol==0 ) THEN
        ycol = k
      ELSEIF( zcol==0 ) THEN
        zcol = k
      ENDIF
    ENDIF
    180 CONTINUE
  ENDDO
  !Check that the number of fields is the same as first line
  IF( k.NE.Ncol ) THEN
    CALL ATOMSK_MSG(1708,(/""/),(/DBLE(Nlines),DBLE(Ncol),DBLE(k)/))
    GOTO 190
  ENDIF
  !
  190 CONTINUE
  !
ENDDO
!
!
!
200 CONTINUE
IF( verbosity>=4 ) THEN
  WRITE(msg,*) "Counted ", Nlines, " lines and ", Ncol, " columns"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!
!
!Go back to beginning of file
REWIND(30)
IF( xcol>0 .OR. ycol>0 .OR. zcol>0 ) THEN
  !A keyword was recognized on the first line
  !i.e. the first line contains description and not data
  READ(30,'(a4096)',ERR=800,END=800) line
  Nlines = Nlines-1
ENDIF
!
!
!Allocate arrays to store data
ALLOCATE(fields(Ncol))
ALLOCATE(P(Nlines,4))
P(:,:) = 0.d0
IF( Sscol>0 .AND. (Sxcol>0 .OR. Sycol>0 .OR. Szcol>0) ) THEN
  ALLOCATE(S(Nlines,4))
  S(:,:) = 0.d0
ENDIF
!
!
!Read the data line by line
Nlines=0
DO i=1,SIZE(P,1)
  fields(:) = ""
  READ(30,'(a4096)',ERR=800,END=800) line
  Nlines=Nlines+1
  !
  !Separate the different fields
  k=0
  DO WHILE(LEN_TRIM(line)>0)
    k = k+1
    IF( k<=SIZE(fields) ) THEN
      IF( line(1:1)=="," ) line = ADJUSTL(line(2:))
      !Check for quotes
      IF( line(1:1)=='"' ) THEN
        !This field starts with a quote
        !Remove the first quote sign
        line = ADJUSTL(line(2:))
        !Find the next quote
        j = SCAN(line,'"')
        IF( j>1 ) THEN
          fields(k) = TRIM(ADJUSTL(line(1:j-1)))
        ENDIF
        !Find the comma after the quote
        line = ADJUSTL(line(j+1:))
        j = SCAN(line,',')
        IF( j>0 ) THEN
          line = ADJUSTL(line(j:))
        ELSE
          !No more commas
          line = ""
        ENDIF
      ELSE
        !There are no quotes
        !Find the next comma
        j = SCAN(line,',')
        IF( j>1 ) THEN
          fields(k) = TRIM(ADJUSTL(line(1:j-1)))
          line = ADJUSTL(line(j:))
        ELSEIF( j==1 ) THEN
          !This is an empty field
          fields(k) = ""
          line = ADJUSTL(line(j:))
        ELSE
          !No more commas =>
          fields(k) = TRIM(line)
          line = ""
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  !
  !Read values from the fields
  !Box vector components
  IF( i==1 ) THEN
    IF( H11>0 ) READ(fields(H11),*,ERR=250,END=250) H(1,1)
    IF( H12>0 ) READ(fields(H12),*,ERR=250,END=250) H(1,2)
    IF( H13>0 ) READ(fields(H13),*,ERR=250,END=250) H(1,3)
    IF( H21>0 ) READ(fields(H21),*,ERR=250,END=250) H(2,1)
    IF( H22>0 ) READ(fields(H22),*,ERR=250,END=250) H(2,2)
    IF( H23>0 ) READ(fields(H23),*,ERR=250,END=250) H(2,3)
    IF( H31>0 ) READ(fields(H31),*,ERR=250,END=250) H(3,1)
    IF( H32>0 ) READ(fields(H32),*,ERR=250,END=250) H(3,2)
    IF( H33>0 ) READ(fields(H33),*,ERR=250,END=250) H(3,3)
  ENDIF
  !
  250 CONTINUE
  !
  IF( Hcol>0 .AND. i<=9 ) THEN
    !Read components of box vectors
    !(if one or more are missing, just ignore it)
    IF( i<=3 ) THEN
      READ(fields(Hcol),*,ERR=260,END=260) H(1,i)
    ELSEIF( i<=6 ) THEN
      READ(fields(Hcol),*,ERR=260,END=260) H(2,i-3)
    ELSEIF( i<=9 ) THEN
      READ(fields(Hcol),*,ERR=260,END=260) H(3,i-6)
    ENDIF
  ENDIF
  !
  260 CONTINUE
  !
  !Get atom coordinates and species
  IF( xcol>0 ) THEN
    READ(fields(xcol),*,ERR=261,END=261) P(i,1)
  ENDIF
  261 CONTINUE
  IF( ycol>0 ) THEN
    READ(fields(ycol),*,ERR=262,END=262) P(i,2)
  ENDIF
  262 CONTINUE
  IF( zcol>0 ) THEN
    READ(fields(zcol),*,ERR=263,END=263) P(i,3)
  ENDIF
  263 CONTINUE
  !Consider all atoms as hydrogen by default
  P(i,4) = 1.d0
  IF( scol>0 ) THEN
    species = fields(scol)(1:2)
    CALL ATOMNUMBER(species,tempreal)
    IF( tempreal>0.99d0 ) P(i,4) = tempreal
  ENDIF
  !
  !Get coordinates of shells (if any)
  IF( Sxcol>0 ) THEN
    READ(fields(Sxcol),*,ERR=271,END=271) S(i,1)
  ENDIF
  271 CONTINUE
  IF( Sycol>0 ) THEN
    READ(fields(Sycol),*,ERR=272,END=272) S(i,2)
  ENDIF
  272 CONTINUE
  IF( Szcol>0 ) THEN
    READ(fields(Szcol),*,ERR=273,END=273) S(i,3)
  ENDIF
  273 CONTINUE
  IF( Sscol>0 ) THEN
    S(i,4) = 0.d0  !By default consider that this atom has no shell
    species = fields(Sscol)(1:2)
    CALL ATOMNUMBER(species,tempreal)
    IF( tempreal>0.99d0 ) S(i,4) = tempreal
  ENDIF
  !
  290 CONTINUE
  !
ENDDO
!
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(Nlines)/))
nerr = nerr+1
!
!
!
1000 CONTINUE
CLOSE(40)
!
END SUBROUTINE READ_CSV
!
!
END MODULE in_csv
