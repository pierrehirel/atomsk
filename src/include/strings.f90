MODULE strings
!
!**********************************************************************************
!*  STRINGS                                                                       *
!**********************************************************************************
!* This module contains functions manipulating strings.                           *
!**********************************************************************************
!* (C) April 2024 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 16 April 2024                                    *
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
!* List of functions in this file:                                                *
!* REAL2TSTR           transforms a real into a string                            *
!* STRUPCASE           converts all letters to upper case in a string             *
!* STRDNCASE           converts all letters to lower case in a string             *
!* CHARLONG2SHRT       convert a long file name to a shorter one                  *
!* STR_REAL            determines if a string contains a number                   *
!* STR_RMSPACE         removes all blank spaces in a string                       *
!* STR_CHAR2SPACE      replaces a character by blank space in a string            *
!* STR2BOOL            transforms a string into a boolean value                   *
!* BOX2DBLE            transforms a string into a real number                     *
!**********************************************************************************
!
!
USE comv
!
!
CONTAINS
!
!
!********************************************************
!  REAL2TSTR
!  This function returns a string containing a real
!  number with no trailing zeros.
!********************************************************
FUNCTION REAL2TSTR(R) RESULT(Rstring)
!
IMPLICIT NONE
CHARACTER(LEN=64):: Rstring
REAL(dp),INTENT(IN):: R
!
Rstring=''
WRITE(Rstring,*) R
Rstring = ADJUSTR(Rstring)
DO WHILE(Rstring(64:64)=='0')
  Rstring = Rstring(1:63)
  Rstring = ADJUSTR(Rstring)
ENDDO
Rstring = TRIM(ADJUSTL(Rstring))
!
END FUNCTION REAL2TSTR
!
!
!********************************************************
!  STRUPCASE
!  This function reads a string of any length
!  and capitalizes all letters.
!********************************************************
FUNCTION StrUpCase (input_string) RESULT (UC_string)
!
IMPLICIT NONE
CHARACTER(*),PARAMETER:: lower_case = 'abcdefghijklmnopqrstuvwxyz'
CHARACTER(*),PARAMETER:: upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
CHARACTER(*),INTENT(IN):: input_string
CHARACTER(LEN(Input_String)):: UC_string !Upper-Case string that is produced
INTEGER:: i, n
!
IF(LEN(input_string)==0) RETURN
UC_string = input_string
! Loop over string elements
DO i=1,LEN(UC_string)
  !Find location of letter in lower case constant string
  n = INDEX( lower_case, UC_string(i:i) )
  !If current substring is a lower case letter, make it upper case
  IF(n>0) THEN
    UC_string(i:i) = upper_case(n:n)
  ENDIF
END DO
!
END FUNCTION StrUpCase
!
!
!********************************************************
!  STRDNCASE
!  This function reads a string of any length
!  and transforms all letters to lower case.
!********************************************************
FUNCTION StrDnCase (input_string) RESULT (lc_string)
!
IMPLICIT NONE
CHARACTER(*),PARAMETER:: lower_case = 'abcdefghijklmnopqrstuvwxyz'
CHARACTER(*),PARAMETER:: upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
CHARACTER(*),INTENT(IN):: input_string
CHARACTER(LEN(Input_String)):: lc_string !Lower-Case string that is produced
INTEGER:: i, n
!
IF(LEN(input_string)==0) RETURN
lc_string = input_string
! Loop over string elements
DO i=1,LEN(lc_string)
  !Find location of letter in upper case constant string
  n = INDEX( upper_case, lc_string(i:i) )
  !If current substring is an upper case letter, make it lower case
  IF(n>0) THEN
    lc_string(i:i) = lower_case(n:n)
  ENDIF
END DO
!
END FUNCTION StrDnCase
!
!
!
!********************************************************
! CHARLONG2SHRT
! This function modifies a long string into a shorter
! one, with dots (...) in the middle. This is typically
! used to display long file names, e.g.:
! /a/very/very/ridiculously/loooong/path/to/my/file
! would be replaced by something like:
! /a/very/.../to/my/file
!********************************************************
FUNCTION CHARLONG2SHRT(longname) RESULT (shortname)
!
IMPLICIT NONE
CHARACTER(*):: longname
CHARACTER(LEN=64):: shortname
CHARACTER(LEN=64):: part1, part2
CHARACTER(LEN=LEN(longname)):: temp
INTEGER:: i, j
!
IF( LEN_TRIM(longname)>64 ) THEN
  !Look for the last path separator
  i = SCAN(longname,pathsep,BACK=.TRUE.)
  IF(i>0) THEN
    !Save file name in part2
    part2 = longname(i:)
    !Look for the separator before that
    temp = longname(1:i-1)
    i = SCAN(temp,pathsep,BACK=.TRUE.)
    IF( i>0 .AND. i>=32 ) THEN
      !Look for the separator before that
      temp = longname(1:i-1)
      j = SCAN(temp,pathsep,BACK=.TRUE.)
      IF( j>0 .AND. j>=32 ) THEN
        !part2 will contain the two last folder names + file name
        part2 = longname(j:)
      ELSE
        !part2 will contain the last folder name + file name
        part2 = longname(i:)
      ENDIF
    ENDIF
    !part1: look for a path separator before the 32-nd character
    part1 = longname(1:32)
    i = SCAN(part1,pathsep,BACK=.TRUE.)
    IF(i>0) THEN
      !Cut the first part at this separator
      part1 = longname(1:i)
    ELSE
      !No separator => part1 will contain the first 28 characters
      part1 = longname(1:28)
    ENDIF
  ELSE
    !no path separator at all, it must be a veeeery long file name...
    !=> part1 will contain the first 28 characters
    !   part2 will contain the last 28 characters
    part1 = longname(1:28)
    i=LEN_TRIM(longname)-28
    part2 = longname(i:)
  ENDIF
  !
  shortname = TRIM(ADJUSTL(part1))//"..."//TRIM(ADJUSTL(part2))
  !
ELSE
  shortname = ADJUSTL(longname)
ENDIF
!
END FUNCTION CHARLONG2SHRT
!
!
!
!********************************************************
!  STR_REAL
!  This function determines if a string can be
!  interpreted as a number (integer or real).
!********************************************************
LOGICAL FUNCTION STR_REAL(string)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: string
REAL(dp):: number
!
STR_REAL = .FALSE.
!
READ(string,*,ERR=100,END=100) number
STR_REAL=.TRUE.
100 CONTINUE
!
END FUNCTION STR_REAL
!
!
!
!********************************************************
! STR_RMSPACE
! This subroutine removes all blank spaces in a string,
!
!********************************************************
SUBROUTINE STR_RMSPACE(string)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(INOUT):: string
INTEGER:: i
!
IF( LEN_TRIM(string) > 1 ) THEN
  DO i=1,LEN_TRIM(string)
    IF( string(i:i)==" " ) THEN
      string(i:) = string(i+1:)
    ENDIF
  ENDDO
ENDIF
!
END SUBROUTINE STR_RMSPACE
!
!
!
!********************************************************
! STR_CHAR2SPACE
! This subroutine parses a string, replacing some
! characters with a space character.
!********************************************************
SUBROUTINE STR_CHAR2SPACE(string,characters)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(INOUT):: string
CHARACTER(LEN=*),INTENT(IN):: characters
INTEGER:: i, j
!
IF( LEN_TRIM(string) > 0 ) THEN
  DO i=1,LEN_TRIM(string)
    DO j=1,LEN_TRIM(characters)
      IF( string(i:i)==characters(j:j) ) THEN
        string(i:i) = " "
      ENDIF
    ENDDO
  ENDDO
ENDIF
!
END SUBROUTINE STR_CHAR2SPACE
!
!
!
!********************************************************
! STR2BOOL
! This subroutine transforms a string into a
! boolean value
!********************************************************
SUBROUTINE STR2BOOL(string,bool)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: string
LOGICAL:: bool
!
SELECT CASE(string)
CASE('n','N','no','NO','No','nO','f','F','false','False','FALSE','.false.','.FALSE.','0','OFF','Off','off')
  bool = .FALSE.
CASE('y','Y','yes','YES','Yes','yES','yEs','yeS','yeah','YEAH','t','T','true','True','TRUE','.true.','.TRUE.','1','ON','On','on')
  bool = .TRUE.
CASE DEFAULT
  WRITE(*,*) 'X!X ERROR: unable to convert this string to a boolean: ', TRIM(string)
  nerr=nerr+1
END SELECT
!
END SUBROUTINE STR2BOOL
!
!
!
!********************************************************
! BOX2DBLE
! This subroutine transforms a string containing one
! of the keywords "BOX" (or "box"), "INF" or "-INF",
! or "%", into a real number.
! The box vector vbox must be provided.
! NOTE: this routine will work only for simple
!      expressions involving only one occurence of
!      "BOX", one operator, and one real number.
!      More complex strings (like "2*BOX-0.2") will
!      NOT be understood!!
!********************************************************
SUBROUTINE BOX2DBLE( vbox,text,number, status )
!
IMPLICIT NONE
CHARACTER(LEN=1):: op
CHARACTER(LEN=128):: text
INTEGER:: posbox, posop  !position of 'BOX' and operator in the string "text"
INTEGER,INTENT(OUT):: status !return status of this routine (0=success; >0=error)
REAL(dp):: areal
REAL(dp):: numbersign
REAL(dp),INTENT(OUT):: number !final real number that is returned
REAL(dp),DIMENSION(3),INTENT(IN):: vbox    !components of box vectors in current direction
!
posbox=0
status = 0
number = 0.d0
numbersign = 1.d0
!
IF( INDEX(text,"-INF")>0 .OR. INDEX(text,"-inf")>0 ) THEN
  number = -1.d0*HUGE(1.d0)
  !
ELSEIF( text=="INF" .OR. INDEX(text,"inf")>0 ) THEN
  number = 1.d0*HUGE(1.d0)
  !
ELSEIF( text=='BOX' .OR. text=='box' ) THEN
  !Total length of the box
  number = MAX(0.d0,MAXVAL(vbox)) + DABS(MIN(0.d0,MINVAL(vbox)))
  !
ELSEIF( text=='-BOX' .OR. text=='-box' ) THEN
  !Opposite of total length of the box
  number = -1.d0 * ( MAX(0.d0,MAXVAL(vbox)) + DABS(MIN(0.d0,MINVAL(vbox))) )
  !
ELSEIF( INDEX(text,'BOX')>0 .OR. INDEX(text,'box')>0) THEN
  !There is an operation, e.g. '0.6*BOX' or 'BOX/2'
  !We have to convert that to a number
  !TODO: replace "BOX" by a number and call a parser (see GNU libmatheval)
  !
  !Check if the text starts with a sign (+ or -)
  !If so, remove it
  IF( text(1:1)=="-" ) THEN
    numbersign = -1.d0
    text = text(2:)
  ELSEIF( text(1:1)=="+" ) THEN
    numbersign = +1.d0
    text = text(2:)
  ENDIF
  !
  text = ADJUSTL(text)
  !Detect the position of the keyword "BOX"
  posbox = INDEX(text,'BOX')
  IF(posbox==0) posbox = INDEX(text,'box')
  !
  !Detect where the operator is: it must be just before or just after "box"
  op = ""
  op = text(posbox+3:posbox+3)
  posop = SCAN(op,'*+-:/')
  IF(posop>0) THEN
    !Save position of operator in text
    posop = posbox+3
  ELSE
    !operator was not after "box" => search before "box"
    op = text(posbox-1:posbox-1)
    posop = SCAN(op,'*+-:/')
    IF(posop>0) THEN
      !Save position of operator in text
      posop = posbox-1
    ELSE
      !no operator => assume it is a multiplication
      op="*"
      IF( posbox==1 ) THEN
        posop = posbox+2
      ELSE
        posop = posbox
      ENDIF
    ENDIF
  ENDIF
  !
  !Evaluate the expression and save it to "number"
  IF( posop<=posbox ) THEN
    !Read what is before 'BOX'
    READ(text(1:posop-1),*,END=800,ERR=800) areal
    !Do the operation
    IF( op=='+' ) THEN
      number = areal + MAX(0.d0,MAXVAL(vbox))
    ELSEIF( op=='-' ) THEN
      number = areal - MAX(0.d0,MAXVAL(vbox))
    ELSEIF( op=='*' ) THEN
      number = areal * SUM(vbox)
    ELSEIF( op==':' .OR. op=='/' ) THEN
      IF( DABS(vbox(1))<1.d-12 .OR. DABS(vbox(2))<1.d-12 .OR. DABS(vbox(3))<1.d-12 ) THEN
        number = HUGE(areal)
      ELSE
        number = areal / SUM(vbox)
      ENDIF
    ENDIF
  ELSE
    !Read what is after the operator
    READ(text(posop+1:),*,END=800,ERR=800) areal
    !Do the operation
    IF( op=='+' ) THEN
      number = MAX(0.d0,MAXVAL(vbox)) + areal
    ELSEIF( op=='-' ) THEN
      number = MAX(0.d0,MAXVAL(vbox)) - areal
    ELSEIF( op=='*' ) THEN
      number = areal * SUM(vbox)
    ELSEIF( op==':' .OR. op=='/' ) THEN
      IF( DABS(areal)<1.d-12 ) THEN
        number = HUGE(vbox)
      ELSE
        number = SUM(vbox) / areal
      ENDIF
    ENDIF
  ENDIF
  !
  !If the text was starting with a sign, apply it
  number = numbersign * number
  !
ELSEIF( INDEX(text,'%')>0 ) THEN
  !This is a number given in percent: read the number and divide it by 100
  posbox=INDEX(text,'%')
  text(posbox:posbox) = " "
  READ( text,*,END=800,ERR=800) areal
  number = SUM(vbox) * areal/100.d0
  !
ELSE
  !easy case: there should just be a number
  READ(text,*,END=800,ERR=800) number
  !
ENDIF
!
RETURN
!
800 CONTINUE
!something failed
status=1
!
END SUBROUTINE BOX2DBLE
!
!
!
END MODULE strings
