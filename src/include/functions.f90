MODULE functions
!
!**********************************************************************************
!*  FUNCTIONS                                                                     *
!**********************************************************************************
!* This module contains functions commonly used by several                        *
!* parts of the ATOMSK program.                                                   *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 09 Nov. 2022                                     *
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
!* DAY_OF_WEEK         finds the day of the week given the date                   *
!* CHARLONG2SHRT       convert a long file name to a shorter one                  *
!* IS_INTEGER          determines if a real number is an integer                  *
!* IS_REAL             determines if a string contains a number                   *
!**********************************************************************************
!
!
USE comv
USE constants
USE math
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
!********************************************************
!  DAY_OF_WEEK
!  This function determines the day of the week
!  (Monday, Tuesday...) from the date.
!********************************************************
FUNCTION DAY_OF_WEEK(day,month,year)
!
IMPLICIT NONE
INTEGER:: Day_of_week, j, k, mm, yy
INTEGER,INTENT(IN) :: day, month, year
!
mm=month
yy=year
IF(mm.le.2) THEN
  mm=mm+12
  yy=yy-1
END IF
j = yy / 100
k = MOD(yy, 100)
Day_of_week = MOD( DBLE(day) + ((DBLE(mm)+1.d0)*26.d0)/10.d0 + DBLE(k) + DBLE(k)/4.d0 + DBLE(j)/4.d0 + 5.d0*DBLE(j), 7.d0)
!
END FUNCTION Day_of_week
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
    IF( i>0 .AND. i>=36 ) THEN
      !Look for the separator before that
      temp = longname(1:i-1)
      j = SCAN(temp,pathsep,BACK=.TRUE.)
      IF( j>0 .AND. j>=36 ) THEN
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
!  IS_INTEGER
!  This function determines if a real number can be
!  interpreted as an integer. For instance 1.d0, 2.d0
!  will return TRUE, but 1.0001d0 will return FALSE.
!********************************************************
LOGICAL FUNCTION IS_INTEGER(number,th)
!
IMPLICIT NONE
REAL(dp),INTENT(IN):: number
REAL(dp),INTENT(IN):: th     !threshold to decide if it is an integer
!
IS_INTEGER = .FALSE.
!
IF( DABS( DBLE(NINT(number)) - number ) < th ) THEN
  IS_INTEGER = .TRUE.
ENDIF
!
END FUNCTION IS_INTEGER
!
!
!
!********************************************************
!  IS_REAL
!  This function determines if a string can be
!  interpreted as a number (integer or real).
!********************************************************
LOGICAL FUNCTION IS_REAL(string)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: string
REAL(dp):: number
!
IS_REAL = .FALSE.
!
READ(string,*,ERR=100,END=100) number
IS_REAL=.TRUE.
100 CONTINUE
!
END FUNCTION IS_REAL
!
!
!
END MODULE functions
