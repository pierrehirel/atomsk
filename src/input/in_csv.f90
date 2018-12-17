MODULE in_csv 
!
!**********************************************************************************
!*  IN_CSV                                                                        *
!**********************************************************************************
!* This module reads a file in Comma-Separated Values (CSV) format.               *
!* The CSV format has no standard specification. A description is given here:     *
!*    https://tools.ietf.org/html/rfc4180                                         *
!* Accepted separators are commas (,), semicolons (;), pipe (|),                  *
!* blanck spaces, and tabulation characters. If one of these is present in the    *
!* first line, then it is assumed to be the separator in the whole file.          *
!* In some languages (e.g. French) the comma separates decimals, e.g. the         *
!* number 1.2 is written "1,2". Therefore, if one of the possible separators      *
!* listed above exists in the first line, afterwards all commas are replaced      *
!* by dots to make the numbers readable.                                          *
!* If the first line contains keywords like species, x, y, z, H, etc.             *
!* (as written by the module "out_csv.f90" for instance), then they are used      *
!* to determine what information is in which column.                              *
!* If the first line does not contain such information, then this module will     *
!* attempt to find real values and consider that they are the atoms positions.    *
!* Note that these considerations imply a lot of assumptions, but this routine    *
!* will try its best to convert the data from the CSV file into acceptable        *
!* atomic data. Of course if the CSV file is complete gibberish, it will fail.    *
!**********************************************************************************
!* (C) Oct. 2018 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 17 Dec. 2018                                     *
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
CHARACTER(LEN=1):: fs   !field separator. Default: comma (,)
CHARACTER(LEN=2):: species
CHARACTER(LEN=128),DIMENSION(:,:),ALLOCATABLE:: fields  !to store values
CHARACTER(LEN=4096):: line, msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: comment
LOGICAL:: skipfirst     !Should the 1st line be skipped?
INTEGER:: i, j, k
INTEGER:: q1, q2
INTEGER:: Nlines, Ncol  !Number of lines and columns
INTEGER:: Ncomm, Naux   !Number of comments, of auxiliary properties
INTEGER:: Nreal         !Number of real numbers in a column
INTEGER:: commentcol    !column number for comments
INTEGER:: Hcol, H11, H12, H13, H21, H22, H23, H31, H32, H33 !column numbers for box vectors H
INTEGER:: atcol, scol, xcol, ycol, zcol !column numbers for atom type, species, x,y,z
INTEGER:: Sscol, Sxcol, Sycol, Szcol    !column numbers for shells: species, x, y, z
INTEGER,DIMENSION(100):: auxcol          !Idex of columns containing auxiliary properties (assuming <100)
REAL(dp):: P1, P2, P3
REAL(dp):: tempreal
REAL(dp),DIMENSION(3,3),INTENT(OUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: S
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: AUX !auxiliary properties
!
!
!Initialize variables
skipfirst=.FALSE.
fs = ','  !Default field separator is the comma
Nlines=0
Ncol=0
Ncomm=0
Naux=0
atcol = 0
auxcol(:) = 0
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
IF(ALLOCATED(fields)) DEALLOCATE(fields)
!
!
100 CONTINUE
msg = 'entering READ_CSV'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
OPEN(UNIT=30,FILE=inputfile,STATUS='UNKNOWN',ERR=800)
!
! Read first line and determine what is the separator
READ(30,'(a4096)',END=800,ERR=800) line
line = ADJUSTL(line)
temp = line
i = 1
DO WHILE( i<=LEN_TRIM(temp) )
  IF( temp(i:i)=='"' ) THEN
    !This starts a field => go on until next quote character
    j=i+1
    DO WHILE( temp(j:j).NE.'"' )
      j=j+1
    ENDDO
    i=j
  ELSEIF( temp(i:i)==';' ) THEN
    !Use semicolon as separator
    fs = ';'
  ELSEIF( temp(i:i)=='|' ) THEN
    !Use pipe as separator
    fs = '|'
  ELSEIF( temp(i:i)==' ' ) THEN
    !Use blanck space as separator
    fs = ' '
  ELSEIF( temp(i:i)==ACHAR(9) ) THEN
    !Use tabulation character as separator
    fs = ACHAR(9)
  ENDIF
  i=i+1
ENDDO
!
!
!Go back to beginning of file
REWIND(30)
!Count total number of lines and columns
Nlines = 0
Ncol = 0
DO
  READ(30,'(a4096)',ERR=200,END=110) line
  110 CONTINUE
  line = ADJUSTL(line)
  IF( LEN_TRIM(line)>0 ) THEN
    Nlines = Nlines+1
    !Count number of fields in this line
    k=0
    DO WHILE(LEN_TRIM(line)>0)
      k = k+1
      IF( line(1:1)==fs ) line = ADJUSTL(line(2:))
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
        !Find the separator after the quote
        line = ADJUSTL(line(j+1:))
        j = SCAN(line,fs)
        IF( j>0 ) THEN
          line = ADJUSTL(line(j:))
        ELSE
          !No more separators
          line = ""
        ENDIF
      ELSE
        !There are no quotes
        !Find the next separator
        j = SCAN(line,fs)
        IF( j>1 ) THEN
          temp = TRIM(ADJUSTL(line(1:j-1)))
          line = ADJUSTL(line(j:))
        ELSEIF( j==1 ) THEN
          !This is an empty field
          temp = ""
          line = ADJUSTL(line(j:))
        ELSE
          !No more separator
          temp = TRIM(line)
          line = ""
        ENDIF
      ENDIF
    ENDDO
    !Check that the number of fields is the same as first line
    IF( Nlines>1 .AND. k.NE.Ncol ) THEN
      CALL ATOMSK_MSG(1708,(/""/),(/DBLE(Nlines),DBLE(Ncol),DBLE(k)/))
      GOTO 190
    ENDIF
    !
    190 CONTINUE
    Ncol = MAX(k,Ncol)
    !
  ENDIF
  !
ENDDO
!
!
!
200 CONTINUE
IF( verbosity>=4 ) THEN
  WRITE(msg,*) "Counted ", Nlines, " lines X ", Ncol, " columns"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!
IF( Nlines<=0 .OR. Ncol<=0 ) THEN
  nerr = nerr+1
  GOTO 1000
ENDIF
!
!Allocate the array fields(:,:) to store the data
ALLOCATE( fields(Nlines,Ncol) )
fields(:,:) = ""
!
!Go back to beginning of file
REWIND(30)
!Read all data and store it into the array fields(:,:)
DO i=1,Nlines
  READ(30,'(a4096)',ERR=300,END=300) line
  line = ADJUSTL(line)
  !
  !Read fields in this line
  !(some lines may contain an incorrect number of fields)
  k=0
  DO WHILE(LEN_TRIM(line)>0)
    k = k+1
    temp = ""
    IF( line(1:1)==fs ) line = ADJUSTL(line(2:))
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
      !Find the separator after the quote
      line = ADJUSTL(line(j+1:))
      j = SCAN(line,fs)
      IF( j>0 ) THEN
        line = ADJUSTL(line(j:))
      ELSE
        !No more separators
        line = ""
      ENDIF
    ELSE
      !There are no quotes
      !Find the next separator
      j = SCAN(line,fs)
      IF( j>1 ) THEN
        temp = TRIM(ADJUSTL(line(1:j-1)))
        line = ADJUSTL(line(j:))
      ELSEIF( j==1 ) THEN
        !This is an empty field
        temp = ""
        line = ADJUSTL(line(j:))
      ELSE
        !No more separator
        temp = TRIM(line)
        line = ""
      ENDIF
    ENDIF
    !
    fields(i,k) = TRIM(ADJUSTL(temp))
    !
  ENDDO  !End loop on columns
  !
ENDDO  !End loop on lines
!
!
!
300 CONTINUE
IF( verbosity>=4 ) THEN
  WRITE(msg,*) "Data successfully stored in memory"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "20 first lines and 10 first columns:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1, MIN(20,SIZE(fields,1))
    WRITE(msg,'(10a12)') (fields(i,j), j=1,MIN(10,SIZE(fields,2)))
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
ENDIF
IF( .NOT. ANY(LEN_TRIM(fields)>0) ) THEN
  nerr = nerr+1
  GOTO 1000
ENDIF
!Now, all the data is stored in the form of strings in the array fields(:,:)
!It may contain text, numbers, etc. (let's not assume anything)
!Try to interpret all this data as atomic coordinates
!
!First line: try to detect keywords among the fields (species, x, y, z...)
Ncomm = 0
DO j=1,SIZE(fields,2)
  line = TRIM(ADJUSTL(fields(1,j)))
  !
  IF( line(1:7)=="species" .OR. line(1:7)=="SPECIES" .OR. line(1:7)=="Species" ) THEN
    scol = j
    skipfirst=.TRUE.
  ELSEIF( scol==0 .AND. ( line(1:4)=="atom" .OR. line(1:4)=="ATOM" .OR. &
        & line(1:4)=="type" .OR. line(1:4)=="TYPE" ) ) THEN
    IF( scol==0) scol = j
    skipfirst=.TRUE.
  ELSEIF ( line(1:1)=="x" .OR. line(1:1)=="X" ) THEN
    xcol = j
    !Verify if this column contains something else than real numbers
    DO i=2,MIN(20,SIZE(fields,1))
      IF( .NOT.IS_REAL(fields(i,j)) ) THEN
        xcol=0
        GOTO 340
      ENDIF
    ENDDO
    skipfirst=.TRUE.
  ELSEIF ( line(1:1)=="y" .OR. line(1:1)=="Y" ) THEN
    ycol = j
    !Verify if this column contains something else than real numbers
    DO i=2,MIN(20,SIZE(fields,1))
      IF( .NOT.IS_REAL(fields(i,j)) ) THEN
        ycol=0
        GOTO 340
      ENDIF
    ENDDO
    skipfirst=.TRUE.
  ELSEIF ( line(1:1)=="z" .OR. line(1:1)=="Z" ) THEN
    zcol = j
    !Verify if this column contains something else than real numbers
    DO i=2,MIN(20,SIZE(fields,1))
      IF( .NOT.IS_REAL(fields(i,j)) ) THEN
        zcol=0
        GOTO 340
      ENDIF
    ENDDO
    skipfirst=.TRUE.
  ELSEIF( line(1:6)=="shells" .OR. line(1:6)=="SHELLS" .OR. line(1:2)=="S " ) THEN
    Sscol = j
    skipfirst=.TRUE.
  ELSEIF ( line(1:2)=="sx" .OR. line(1:2)=="SX" ) THEN
    Sxcol = j
    !Verify if this column contains something else than real numbers
    DO i=2,MIN(20,SIZE(fields,1))
      IF( .NOT.IS_REAL(fields(i,j)) ) THEN
        Sxcol=0
        GOTO 340
      ENDIF
    ENDDO
    skipfirst=.TRUE.
  ELSEIF ( line(1:2)=="sy" .OR. line(1:2)=="SY" ) THEN
    Sycol = j
    !Verify if this column contains something else than real numbers
    DO i=2,MIN(20,SIZE(fields,1))
      IF( .NOT.IS_REAL(fields(i,j)) ) THEN
        Sycol=0
        GOTO 340
      ENDIF
    ENDDO
    skipfirst=.TRUE.
  ELSEIF ( line(1:2)=="sz" .OR. line(1:2)=="SZ" ) THEN
    Szcol = j
    !Verify if this column contains something else than real numbers
    DO i=2,MIN(20,SIZE(fields,1))
      IF( .NOT.IS_REAL(fields(i,j)) ) THEN
        Szcol=0
        GOTO 340
      ENDIF
    ENDDO
    skipfirst=.TRUE.
  ELSEIF ( line(1:3)=="H11" ) THEN
    !Box vector components will be in different columns
    H11 = j
    skipfirst=.TRUE.
  ELSEIF ( line(1:3)=="H12" ) THEN
    H12 = j
    skipfirst=.TRUE.
  ELSEIF ( line(1:3)=="H13" ) THEN
    H13 = j
    skipfirst=.TRUE.
  ELSEIF ( line(1:3)=="H21" ) THEN
    H21 = j
    skipfirst=.TRUE.
  ELSEIF ( line(1:3)=="H22" ) THEN
    H22 = j
    skipfirst=.TRUE.
  ELSEIF ( line(1:3)=="H23" ) THEN
    H23 = j
    skipfirst=.TRUE.
  ELSEIF ( line(1:3)=="H31" ) THEN
    H31 = j
    skipfirst=.TRUE.
  ELSEIF ( line(1:3)=="H32" ) THEN
    H32 = j
    skipfirst=.TRUE.
  ELSEIF ( line(1:3)=="H33" ) THEN
    H33 = j
    skipfirst=.TRUE.
  ELSEIF ( line(1:1)=="H" ) THEN
    !All components of box vectors are in this column
    Hcol = j
    skipfirst=.TRUE.
  ELSEIF ( line(1:7)=="comment" ) THEN
    !Consider this column as containing comments
    !NOTE: only one column of comments will be used here
    commentcol = j
    skipfirst=.TRUE.
    !Count how many non-empty fields exist in that column (max.20 comments)
    DO i=2,MIN(21,SIZE(fields,1))
      IF( LEN_TRIM(fields(i,commentcol)) > 0 ) THEN
        Ncomm = Ncomm+1
      ENDIF
    ENDDO
  ELSEIF ( line(1:1)=="#" ) THEN
    !Consider that this column contains comments
    !(only if no column labelled "comments" was found)
    IF( commentcol==0) THEN
      commentcol = j
    ENDIF
    IF( commentcol==j ) THEN
      Ncomm = Ncomm+1
    ENDIF
  ELSEIF( LEN_TRIM(line) > 0 ) THEN
    !This field contains something
    IF( skipfirst ) THEN
      !Other keywords were recognized on the first line
      !Consider that this field is the name of some auxiliary property
      Naux = Naux+1
      IF( Naux<=SIZE(auxcol) ) THEN
        auxcol(Naux) = j
      ENDIF
    ENDIF
  ENDIF
  340 CONTINUE
ENDDO
!
!
!If we don't know where are the coordinates x, y, z,
!try to find them among the data
!NOTE: the first column that contains integers that are successfully converted
!     into atom chemical species is used in that way. Similarly, the first
!     column containing real numbers is marked as containin the X coordinates
!     of atoms, the second column Y, and the third Z
IF( (.NOT.skipfirst) .OR. (scol==0 .AND. ( xcol==0 .OR. ycol==0 .OR. zcol==0 )) ) THEN
  !
  WRITE(msg,*) "Unable to read field names from first line, analyzing data to guess format..."
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  DO j=1,SIZE(fields,2)   ! loop on columns
    !Make sure we don't replace a column that is already known
    IF( j.NE.scol .AND. j.NE.xcol .AND. j.NE.ycol .AND. j.NE.zcol ) THEN
      k = 0
      Nreal = 0
      DO i=MIN(2,SIZE(fields,1)) , MAX(10,SIZE(fields,1))   !loop on lines
        !Try to read a number from this field
        READ(fields(i,j),*,ERR=350,END=350) tempreal
        !A real number was read: try to interpret it
        IF( tempreal>0.d0 .AND. DABS( DBLE(NINT(tempreal)) - tempreal ) < 1.d-12 ) THEN
          !It is actually a positive integer number: check if it is an atomic number
          CALL ATOMSPECIES(tempreal,species)
          IF( species .NE. "XX" ) THEN
            IF(scol==0) scol = j
          ENDIF
        ELSE
          !It is not an integer but a real number
          !=> increase the counter
          k = k+1
        ENDIF
        GOTO 360
        350 CONTINUE
        !Reading a number failed => try to interpret this text as an atom species
        IF( LEN_TRIM(fields(i,j)) <= 2 ) THEN
          CALL ATOMNUMBER(fields(i,j),tempreal)
          IF( tempreal >= 0.9d0 ) THEN
            !It is an atom species: mark this column
            IF(scol==0) scol = j
          ELSE
            !It is a short text we cannot understand: ignore it
          ENDIF
        ELSE
          !String too long to be an atom species: ignore it
        ENDIF
        360 CONTINUE
      ENDDO  !end loop on i
      !
      !Now k is the number of reals that were read
      IF( k >= MIN(10,SIZE(fields,1))-2 ) THEN
        !This column is filled with real numbers
        IF( xcol==0 ) THEN
          xcol = j
        ELSEIF( ycol==0 ) THEN
          ycol = j
        ELSEIF( zcol==0 ) THEN
          zcol = j
        ELSE
          !The columns for X, Y, Z are already known
          !but column #j still contains real numbers
          Nreal = Nreal+1
        ENDIF
      ENDIF
      !
    ENDIF
    !
  ENDDO  !end loop on j
  !
  IF( Nreal >= SIZE(fields,1)-1 ) THEN
    !Colum #j is filled with real numbers that are not atom coordinates
    !Consider that this column contains some auxiliary property
    Naux = Naux+1
    auxcol(Naux) = j
  ENDIF
  !
ENDIF
!
!
!
400 CONTINUE
IF( verbosity>=4 ) THEN
  WRITE(msg,*) "scol, xcol, ycol, zcol: ", scol, xcol, ycol, zcol
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "Number of aux.prop.: ", Naux
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "Number of comments: ", Ncomm
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!
!At this point we should know which column contains what
!If not, display an error and exit
IF( xcol==0 .AND. ycol==0 .AND. zcol==0 ) THEN
  CALL ATOMSK_MSG(1804,(/""/),(/0.d0/))
  nerr = nerr+1
  GOTO 1000
ENDIF
!
!Box vector components (if one fails, just go on)
H(:,:) = 0.d0
IF( H11>0 ) READ(fields(2,H11),*,ERR=450,END=450) H(1,1)
450 CONTINUE
IF( H12>0 ) READ(fields(2,H12),*,ERR=451,END=451) H(1,2)
451 CONTINUE
IF( H13>0 ) READ(fields(2,H13),*,ERR=452,END=452) H(1,3)
452 CONTINUE
IF( H21>0 ) READ(fields(2,H21),*,ERR=453,END=453) H(2,1)
453 CONTINUE
IF( H22>0 ) READ(fields(2,H22),*,ERR=454,END=454) H(2,2)
454 CONTINUE
IF( H23>0 ) READ(fields(2,H23),*,ERR=455,END=455) H(2,3)
455 CONTINUE
IF( H31>0 ) READ(fields(2,H31),*,ERR=456,END=456) H(3,1)
456 CONTINUE
IF( H32>0 ) READ(fields(2,H32),*,ERR=457,END=457) H(3,2)
457 CONTINUE
IF( H33>0 ) READ(fields(2,H33),*,ERR=458,END=458) H(3,3)
458 CONTINUE
!
IF( Hcol>0 .AND. SIZE(fields,1)>=10 ) THEN
  !Read components of box vectors
  !(if one or more are missing, just ignore it)
  DO i=2,4
    READ(fields(i,Hcol),*,ERR=460,END=460) H(1,i-1)
  ENDDO
  DO i=5,7
    READ(fields(i,Hcol),*,ERR=460,END=460) H(2,i-4)
  ENDDO
  DO i=8,10
    READ(fields(i,Hcol),*,ERR=460,END=460) H(3,i-7)
  ENDDO
ENDIF
!
460 CONTINUE
IF( verbosity>=4 ) THEN
  WRITE(msg,*) "READ BOX VECTORS H: "
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "         H(1,:) = ", H(1,1), H(1,2), H(1,3)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "         H(2,:) = ", H(2,1), H(2,2), H(2,3)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "         H(3,:) = ", H(3,1), H(3,2), H(3,3)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!
!Allocate array to store atom positions
k = SIZE(fields,1)
IF( skipfirst ) k=k-1
ALLOCATE(P(k,4))
P(:,:) = 0.d0
IF( verbosity>=4 ) THEN
  WRITE(msg,*) "ATOMS ARRAY SIZE: ", SIZE(P,1), " X ", SIZE(P,2)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!Allocate array for ionic shells (if any)
IF( Sscol>0 .AND. (Sxcol>0 .OR. Sycol>0 .OR. Szcol>0) ) THEN
  ALLOCATE(S(SIZE(P,1),4))
  S(:,:) = 0.d0
  IF( verbosity>=4 ) THEN
    WRITE(msg,*) "SHELLS ARRAY SIZE: ", SIZE(S,1), " X ", SIZE(S,2)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
ENDIF
!Allocate arrays for auxiliary properties (if any)
IF( Naux>0 ) THEN
  ALLOCATE(AUXNAMES(Naux))
  AUXNAMES(:) = ""
  IF( skipfirst ) THEN
    !Get auxiliary properties names from first line
    i=1
    j=0
    DO WHILE( auxcol(i)>0 )
      j=j+1
      AUXNAMES(j) = TRIM(ADJUSTL( fields(1,auxcol(i)) ))
      i=i+1
    ENDDO
  ELSE
    !Fill AUXNAMES with dummy property names
    i=1
    j=0
    DO WHILE( auxcol(i)>0 )
      j=j+1
      WRITE(temp,*) j
      temp = "property_"//TRIM(ADJUSTL(temp))
      AUXNAMES(i) = TRIM(ADJUSTL(temp))
      i=i+1
    ENDDO
  ENDIF
  ALLOCATE( AUX(SIZE(P,1),Naux) )
  AUX(:,:) = 0.d0
  IF( verbosity>=4 ) THEN
    WRITE(msg,*) "AUX ARRAY SIZE: ", SIZE(AUX,1), " X ", SIZE(AUX,2)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,*) "Columns for auxiliary properties:"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    msg = ""
    i=1
    DO WHILE( auxcol(i)>0 )
      WRITE(temp,*) auxcol(i)
      WRITE(msg,*) TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(temp))
      i=i+1
    ENDDO
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    msg = ""
    i=1
    DO WHILE( auxcol(i)>0 )
      WRITE(msg,*) TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(AUXNAMES(i)))
      i=i+1
    ENDDO
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
ENDIF
!Allocate arrays for comments (if any)
IF( Ncomm>0 ) THEN
  ALLOCATE(comment(Ncomm))
  comment(:) = ""
  IF( verbosity>=4 ) THEN
    WRITE(msg,*) "COMMENTS ARRAY SIZE: ", SIZE(comment)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
ENDIF
!
!Read atom data line by line from the array fields(:,:)
Nlines=0
i = 0
Ncomm = 0
DO Ncol=1,SIZE(fields,1)
  !
  IF( skipfirst ) THEN
    IF(Ncol==1) GOTO 490
    i = Ncol-1
  ELSE
    i = Ncol
  ENDIF
  !
  IF( i>0 .AND. i<=SIZE(P,1) ) THEN
    !Get atom coordinates and species
    IF( xcol>0 ) THEN
      READ(fields(Ncol,xcol),*,ERR=461,END=461) P(i,1)
    ENDIF
    461 CONTINUE
    IF( ycol>0 ) THEN
      READ(fields(Ncol,ycol),*,ERR=462,END=462) P(i,2)
    ENDIF
    462 CONTINUE
    IF( zcol>0 ) THEN
      READ(fields(Ncol,zcol),*,ERR=463,END=463) P(i,3)
    ENDIF
    463 CONTINUE
    !
    !Consider all atoms as hydrogen by default
    P(i,4) = 1.d0
    IF( scol>0 ) THEN
      !Try to read a real number
      READ(fields(Ncol,scol),*,END=465,ERR=465) P(i,4)
      GOTO 470   !success: go to next fields
      465 CONTINUE
      !There was an error: this was not a real number
      !Assuming it is a chemical symbol, try to convert it into an atomic number
      species = fields(Ncol,scol)(1:2)
      CALL ATOMNUMBER(species,tempreal)
      IF( tempreal>0.99d0 ) THEN
        !Success: keep this as the atomic number of this atom
        P(i,4) = tempreal
      ENDIF
    ENDIF
    !
    470 CONTINUE
    !Get coordinates of shells (if any)
    IF( ALLOCATED(S) .AND. i<=SIZE(S,1) ) THEN
      IF( Sxcol>0 ) THEN
        READ(fields(Ncol,Sxcol),*,ERR=471,END=471) S(i,1)
      ENDIF
      471 CONTINUE
      IF( Sycol>0 ) THEN
        READ(fields(Ncol,Sycol),*,ERR=472,END=472) S(i,2)
      ENDIF
      472 CONTINUE
      IF( Szcol>0 ) THEN
        READ(fields(Ncol,Szcol),*,ERR=473,END=473) S(i,3)
      ENDIF
      473 CONTINUE
      IF( Sscol>0 ) THEN
        S(i,4) = 0.d0  !By default consider that this atom has no shell
        species = fields(Ncol,Sscol)(1:2)
        CALL ATOMNUMBER(species,tempreal)
        IF( tempreal>0.99d0 ) S(i,4) = tempreal
      ENDIF
    ENDIF
    !
    !Get auxiliary properties (just skip if there are errors)
    IF( ALLOCATED(AUX) .AND. i<=SIZE(AUX,1) ) THEN
      j=1
      Naux=0
      DO WHILE( auxcol(j)>0 .AND. j<=SIZE(AUX,2) )
        Naux=Naux+1
        READ(fields(Ncol,auxcol(j)),*,ERR=475,END=475) AUX(i,Naux)
        j=j+1
      ENDDO
      475 CONTINUE
    ENDIF
    !
    480 CONTINUE
    !
    !Get comments (if any)
    IF( commentcol>0 .AND. ALLOCATED(comment) ) THEN
      IF( LEN_TRIM(fields(Ncol,commentcol)) > 0 ) THEN
        Ncomm = Ncomm+1
        IF( Ncomm <= SIZE(comment) ) THEN
          comment(Ncomm) = fields(Ncol,commentcol)
        ENDIF
      ENDIF
    ENDIF
    !
  ELSE
    GOTO 1000
  ENDIF
  !
  490 CONTINUE
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
IF(ALLOCATED(fields)) DEALLOCATE(fields)
!
END SUBROUTINE READ_CSV
!
!
END MODULE in_csv
