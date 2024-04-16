MODULE exprev
!
!**********************************************************************************
!*  EQUATION                                                                      *
!**********************************************************************************
!* This program interprets a string containing numbers, operations, and           *
!* functions, and evaluates it to produce a unique value (real number).           *
!**********************************************************************************
!* (C) December 2018 - Pierre Hirel                                               *
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
!
USE comv
USE constants
USE messages
USE random
!
!
CONTAINS
!
!
!********************************************************
!  STR_EXP2VAL
!  This subroutine replaces an expression in a string
!  by a value. For example, if the string is:
!    string = "3*smth + 6"
!  and if expr="smth" and value=4, then the string
!  will be modified into:
!    string = "3*4 + 6"
!********************************************************
SUBROUTINE STR_EXP2VAL(string,expr,value,nr)
!
CHARACTER(LEN=*),INTENT(INOUT):: string
CHARACTER(LEN=*),INTENT(IN):: expr
CHARACTER(LEN=4096):: string1, string2, temp
INTEGER:: i1, i2  !start and end position of keyword
INTEGER:: nr      !number of substitutions that were made
REAL(dp),INTENT(IN):: value
!
nr = 0
!
DO WHILE( INDEX(string,expr) > 0 )
  nr = nr+1
  i1=INDEX(string,expr)
  IF( i1>1 ) THEN
    string1 = TRIM(ADJUSTL(string(1:i1-1)))
  ELSE
    string1 = ""
  ENDIF
  i2=i1+LEN_TRIM(expr)
  IF( i2>1 ) THEN
    string2 = TRIM(ADJUSTL(string(i2:)))
  ELSE
    string2 = ""
  ENDIF
  WRITE(temp,*) value
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(string2)
ENDDO

END SUBROUTINE STR_EXP2VAL
!
!
!
!********************************************************
!  EXPREVAL
!  This function reads an input string containing
!  numbers and operations, and evaluates it to produce
!  a real number ("value").
!  This subroutine is recursive and calls itself to
!  evaluate the different parts of the expression.
!  The input string may contain only numbers, operators,
!  and/or functions.
!  The numbers must be real (float) or integers.
!  Complex numbers are not supported.
!  The following operators are supported:
!     +        (addition)
!     -        (subtraction)
!     *        (multiplication)
!     / or :   (division)
!     %        (remainder or modulo)
!     ^ or **  (power)
!     !        (factorial)
!     ()       (parenthesis)
!  Common operator priority is respected.
!  The following functions are supported:
!     cos   acos   exp   sqrt
!     sin   asin   log   int
!     tan   atan   ln    abs
!  The only recognized constant is pi.
!  The input string shall not contain any variables
!  (e.g. "3*x" or "y-2" will not be recognized).
!
!  EXAMPLE:
!  If the input string is "exp(20/6) + cos(pi/3)",
!  then the return value will be 28.531624...
!
!  In case of unmatched opening (or closing) parenthesis,
!  corresponding closing (resp. opening) parenthesis are
!  assumed to be at the end (resp. beginning) of the string.
!  For instance:  "5+6)^2" will be interpreted as "(5+6)^2".
!  "3*(7+8" will be interpreted as "3*(7+8)", etc.
!
!  NOTE:
!  You will probably find this routine "dirty" and
!  sub-optimal, and curse the one who developped it.
!  It has several known flaws and limitations, a few
!  of which are listed below:
!  - It is not accurate. Values are written into a string
!    before being evaluated again, leading to round-off
!    errors that propagate. Also, in some places,
!    values smaller than 10^-15 are rounded off to zero.
!  - It is not optimal. Other faster algorithms exist.
!  - It is not complete: all functions are not supported.
!  - Only real numbers are supported (not complex numbers)
!********************************************************
RECURSIVE FUNCTION EXPREVAL(string,recuri,status) RESULT(eval)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(INOUT):: string
CHARACTER(LEN=1):: ang  !is the angle in degrees or radians (default: radians)
CHARACTER(LEN=7):: operators="+-*/:^%"
CHARACTER(LEN=12):: illegal="&#'|`@°£$øµ§" !list of illegal characters that will lead to an error
CHARACTER(LEN=256):: prefix
CHARACTER(LEN=4096):: string1, string2
CHARACTER(LEN=4096):: temp, temp1, temp2
INTEGER:: i, j, m
INTEGER:: ms  !position of minus sign
INTEGER,INTENT(INOUT):: recuri  !depth of recursion
INTEGER:: p1, p2  !position of opening and closing parenthesis
INTEGER:: status
REAL(dp),PARAMETER:: low_limit=1.d-15   !in trig.functions, values which are smaller than that
                                        !in absolute will be rounded off to zero
REAL(dp):: x, y  !to store numbers
REAL(dp),DIMENSION(:),ALLOCATABLE:: randarray  !random numbers
REAL(dp):: eval
!
!Initialize variables and empty strings
prefix = ""
string1 = ""
string2 = ""
temp = ""
ang="R"    !by default angles are in radians
p1 = 0
p2 = 0
!
!If recursion is too large, exit with an error
IF( recuri>100 ) THEN
  status=recuri
ENDIF
!
IF( recuri<=1 ) status=0
!
IF( status>0 ) RETURN
!
DO i=0,recuri
  prefix = TRIM(prefix)//".."
ENDDO
!
IF(verbosity==4) PRINT*, TRIM(prefix)//"  NEW CALL TO EXPREVAL WITH INPUT STRING: ", TRIM(string)
recuri = recuri+1
!
!In case of empty string, return zero
IF( LEN_TRIM(string)<=0 ) THEN
  eval = 0.d0
  RETURN
ENDIF
!
IF( INDEX(string,"pi")>1 .OR. INDEX(string,"(")>0 .OR. SCAN(string,"[]{}")>0 ) THEN
  !If user used [sqare] or {curly} brackets, replace with regular parenthesis
  i=0
  DO WHILE(i<LEN_TRIM(string))
    i = i+1
    IF( string(i:i)=='[' ) string(i:i)='('
    IF( string(i:i)==']' ) string(i:i)=')'
    IF( string(i:i)=='{' ) string(i:i)='('
    IF( string(i:i)=='}' ) string(i:i)=')'
  ENDDO
  !
  !If user wrote things like "2pi", correct it into "2*pi"
  !Same with functions, e.g. 2cos(pi/3) is corrected into 2*cos(pi/3)
  !Same with parenthesis, e.g. 2(3+1) is corrected into 2*(3+1)
  i=0
  DO WHILE(i<LEN_TRIM(string))
    i = i+1
    m = 0
    IF( string(i:i+1)=="pi" .OR. string(i:i)=='(' .OR.   &
      & string(i:i+2)=="cos" .OR. string(i:i+2)=="sin" .OR.string(i:i+2)=="tan" .OR.    &
      & string(i:i+3)=="acos" .OR. string(i:i+3)=="asin" .OR.string(i:i+3)=="atan" .OR. &
      & string(i:i+4)=="atan2" .OR. string(i:i+3)=="sqrt" .OR.string(i:i+2)=="abs" .OR. &
      & string(i:i+2)=="exp" .OR.string(i:i+2)=="log" .OR.string(i:i+1)=="ln" .OR.      &
      & string(i:i+2)=="int" .OR.string(i:i+2)=="log" .OR.string(i:i+1)=="ln"      ) THEN
      !
      IF( SCAN(string(i-1:i-1),"0123456789")>0 ) THEN
        IF( .NOT. ( i>5 .AND. string(i-5:i-1)=="atan2" ) ) THEN
          !pi or function or parenthesis is multiplied by a number
          string = ADJUSTL(string(1:i-1))//"*"//TRIM(string(i:))
        ENDIF
      ELSEIF( string(i-1:i-1)=='-' ) THEN
        !pi is multiplied by -1
        string = ADJUSTL(string(1:i-2))//"-1.0*"//TRIM(string(i:))
      ENDIF
    ENDIF
  ENDDO
ENDIF
!
!If string contains strange characters, exit with error
IF( SCAN(string,illegal)>0 .OR. SCAN(string,'"')>0 ) THEN
  i = SCAN(string,illegal)
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  ILLEGAL CHARACTER: "//string(i:i)
  status=1
  RETURN
ENDIF
!
!First, try to read a number: if it succeeds, we are done
IF( SCAN(string,operators)==0 .AND. SCAN(string,"+-",BACK=.TRUE.)<=1 ) THEN
  READ(string,*,END=100,ERR=100) eval
  IF( IS_INTEGER(eval,1.d-16) ) eval = DBLE(NINT(eval))
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  NUMBER(1): ", TRIM(string)
  recuri = recuri-1
  RETURN
ENDIF
!
!
100 CONTINUE
!Reading a number failed: try to interpret the expression
eval = 0.d0
!
!Save position of first opening parenthesis (if any)
p1 = SCAN(string,"(")
!Save position of first minus sign
!(that is not the first character, not a negative number, and not part of a mantissa)
ms = SCAN(string(2:),"-")+1
IF( ms>0 ) THEN
  DO WHILE( SCAN(string(ms-1:ms-1),"Ee/*+-:")>0 )
    !The minus sign is part of a mantissa or indicates a negative number: look for another one
    IF( SCAN(string(ms+1:),"-")>0 ) THEN
      !This may be a subtraction sign (will be checked at next loop iteration)
      ms = ms + SCAN(string(ms+1:),"-")
    ELSE
      !No further minus sign => no subtraction in this string
      ms = 0
    ENDIF
  ENDDO
ENDIF
!
!If it is not a function but there is a parenthesis, interpret the parenthesis
IF( p1==0 .AND. SCAN(string,')')>0 ) THEN
  !Something is strange: p1=0 so there is no opening parenthesis,
  !but there is a closing parenthesis: remove it from the string
  i = SCAN(string,')')
  string(i:i) = " "
ENDIF
!
!
IF( INDEX(string,"--")>0 ) THEN
  i = INDEX(string,"--")
  string(i:i) = '+'
  string(i+1:) = string(i+2:)
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( INDEX(string,"+-")>0 ) THEN
  i = INDEX(string,"+-")
  string(i:i) = '-'
  string(i+1:) = string(i+2:)
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( INDEX(string,"-+")>0 ) THEN
  i = INDEX(string,"-+")
  string(i:i) = '-'
  string(i+1:) = string(i+2:)
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( INDEX(string,"**")>0 ) THEN
  i = INDEX(string,"**")
  string(i:i) = '^'
  string(i+1:) = string(i+2:)
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( string=="pi" ) THEN
  eval = pi
  !
ELSE IF( string=="kB" ) THEN
  eval = kB
  !
ELSE IF( string=="hbar" .OR. string=="h_bar" ) THEN
  eval = h_bar
  !
ELSE IF( string=="Navo" .OR. string=="Na" ) THEN
  eval = Navo
  !
ELSE IF( string=="qe" .OR. string=="Qe" ) THEN
  eval = e_charge
  !
ELSE IF( string=="epsilon0" .OR. string=="eps0" .OR. string=="eps_0" ) THEN
  eval = eps_0
  !
ELSE IF( string=="mu0" .OR. string=="mu_0" ) THEN
  eval = mu_0
  !
ELSE IF( string=="rand" .OR. string=="random" ) THEN
  !generate a random number
  CALL GEN_NRANDNUMBERS(1,randarray)
  eval = randarray(1)
  !
!Replace functions by their eval
ELSE IF( INDEX(string,"sqrt(")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED SQUARE ROOT: ", TRIM(string)
  !Save position of the expression "sqrt"
  i = INDEX(string,"sqrt(")
  !Save everything that is before function name into string1
  !(only if function name is not the first thing in string)
  IF( i>1 ) THEN
    string1 = TRIM(ADJUSTL(string(1:i-1)))
  ENDIF
  !Save position of the opening parenthesis
  p1 = i+4
  !Save everything after opening parenthesis into temp
  temp = TRIM(ADJUSTL(string(p1+1:)))
  !Find position of matching closing parenthesis
  m=0
  i=0
  j=1
  p2 = LEN_TRIM(string)
  DO WHILE( m==0 .AND. i<=LEN_TRIM(temp) )
    i=i+1
    IF( temp(i:i)=="(" ) THEN
      j=j+1
    ELSEIF( temp(i:i)==")" ) THEN
      j=j-1
    ENDIF
    IF( j==0 ) THEN
      !i is the position of the matching closing parenthesis
      p2=i
      m=1
      EXIT
    ENDIF
  ENDDO
  !Save everything after closing parenthesis into string2
  string2 = TRIM(ADJUSTL(temp(p2+1:)))
  !Remove from temp everything that is after closing parenthesis
  !NOTE: if no matching closing parenthesis was found, then temp
  !      will contain everything that is after opening parenthesis
  temp = TRIM(ADJUSTL(temp(1:p2-1)))
  !Interpret the expression that is inside the parenthesis, save result in x
  x = EXPREVAL(temp,recuri,status)
  !Eliminate round-off errors
  IF( DABS(x)<low_limit ) x=0.d0
  !Perform the calculation of the square root
  IF( x>=0.d0 ) THEN
    eval = DSQRT(x)
  ELSE
    IF(verbosity==4) PRINT*, " X ! X ERROR: negative eval inside square root"
    status=1
    GOTO 900
  ENDIF
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: sqrt("//TRIM(ADJUSTL(temp))//") = ", eval
  !Write the result into temp
  WRITE(temp,*) eval
  !Re-write the string by concatenating string1, temp, and string2
  string = TRIM(string1)//TRIM(temp)//TRIM(string2)
  !Evaluate the resulting string
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( INDEX(string,"abs(")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED ABSOLUTE VALUE: ", TRIM(string)
  !Save position of the expression "abs"
  i = INDEX(string,"abs(")
  !Save everything that is before function name into string1
  !(only if function name is not the first thing in string)
  IF( i>1 ) THEN
    string1 = TRIM(ADJUSTL(string(1:i-1)))
  ENDIF
  !Save position of the opening parenthesis
  p1 = i+3
  !Save everything after opening parenthesis into temp
  temp = TRIM(ADJUSTL(string(p1+1:)))
  !Find position of matching closing parenthesis
  m=0
  i=0
  j=1
  p2 = LEN_TRIM(string)
  DO WHILE( m==0 .AND. i<=LEN_TRIM(temp) )
    i=i+1
    IF( temp(i:i)=="(" ) THEN
      j=j+1
    ELSEIF( temp(i:i)==")" ) THEN
      j=j-1
    ENDIF
    IF( j==0 ) THEN
      !i is the position of the matching closing parenthesis
      p2=i
      m=1
      EXIT
    ENDIF
  ENDDO
  !Save everything after closing parenthesis into string2
  string2 = TRIM(ADJUSTL(temp(p2+1:)))
  !Remove from temp everything that is after closing parenthesis
  !NOTE: if no matching closing parenthesis was found, then temp
  !      will contain everything that is after opening parenthesis
  temp = temp(1:p2-1)
  !Interpret the expression that is inside the parenthesis, save result in x
  x = EXPREVAL(temp,recuri,status)
  !Perform the calculation of the absolute eval
  eval = DABS(x)
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: |"//TRIM(ADJUSTL(temp))//"| = ", eval
  !Write the result into temp
  WRITE(temp,*) eval
  !Re-write the string by concatenating string1, temp, and string2
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  !Evaluate the resulting string
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( INDEX(string,"int(")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED INT: ", TRIM(string)
  !Save position of the expression "int"
  i = INDEX(string,"int(")
  !Save everything that is before function name into string1
  !(only if function name is not the first thing in string)
  IF( i>1 ) THEN
    string1 = TRIM(ADJUSTL(string(1:i-1)))
  ENDIF
  !Save position of the opening parenthesis
  p1 = i+3
  !Save everything after opening parenthesis into temp
  temp = TRIM(ADJUSTL(string(p1+1:)))
  !Find position of matching closing parenthesis
  m=0
  i=0
  j=1
  p2 = LEN_TRIM(string)
  DO WHILE( m==0 .AND. i<=LEN_TRIM(temp) )
    i=i+1
    IF( temp(i:i)=="(" ) THEN
      j=j+1
    ELSEIF( temp(i:i)==")" ) THEN
      j=j-1
    ENDIF
    IF( j==0 ) THEN
      !i is the position of the matching closing parenthesis
      p2=i
      m=1
      EXIT
    ENDIF
  ENDDO
  !Save everything after closing parenthesis into string2
  string2 = TRIM(ADJUSTL(temp(p2+1:)))
  !Remove from temp everything that is after closing parenthesis
  !NOTE: if no matching closing parenthesis was found, then temp
  !      will contain everything that is after opening parenthesis
  temp = temp(1:p2-1)
  !Interpret the expression that is inside the parenthesis, save result in x
  x = EXPREVAL(temp,recuri,status)
  !Perform the calculation of the nearest integer
  eval = NINT(x)
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: |"//TRIM(ADJUSTL(temp))//"| = ", eval
  !Write the result into temp
  WRITE(temp,*) eval
  !Re-write the string by concatenating string1, temp, and string2
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  !Evaluate the resulting string
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( INDEX(string,"acos(")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED ARCCOSINUS: ", TRIM(string)
  !Save position of the expression "acos"
  i = INDEX(string,"acos(")
  !Save everything that is before function name into string1
  !(only if function name is not the first thing in string)
  IF( i>1 ) THEN
    string1 = TRIM(ADJUSTL(string(1:i-1)))
  ENDIF
  !Save position of the opening parenthesis
  p1 = i+4
  !Save everything after opening parenthesis into temp
  temp = TRIM(ADJUSTL(string(p1+1:)))
  !Find position of matching closing parenthesis
  m=0
  i=0
  j=1
  p2 = LEN_TRIM(string)
  DO WHILE( m==0 .AND. i<=LEN_TRIM(temp) )
    i=i+1
    IF( temp(i:i)=="(" ) THEN
      j=j+1
    ELSEIF( temp(i:i)==")" ) THEN
      j=j-1
    ENDIF
    IF( j==0 ) THEN
      !i is the position of the matching closing parenthesis
      p2=i
      m=1
      EXIT
    ENDIF
  ENDDO
  !Save everything after closing parenthesis into string2
  string2 = TRIM(ADJUSTL(temp(p2+1:)))
  !Remove from temp everything that is after closing parenthesis
  !NOTE: if no matching closing parenthesis was found, then temp
  !      will contain everything that is after opening parenthesis
  temp = temp(1:p2-1)
  !Interpret the expression that is inside the parenthesis, save result in x
  x = EXPREVAL(temp,recuri,status)
  !Perform the calculation of the arccosinus
  eval = DACOS(x)
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: acos("//TRIM(ADJUSTL(temp))//") = ", eval
  !Write the result into temp
  WRITE(temp,*) eval
  !Re-write the string by concatenating string1, temp, and string2
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  !Evaluate the resulting string
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( INDEX(string,"cos(")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED COSINUS: ", TRIM(string)
  !Save position of the expression "cos"
  i = INDEX(string,"cos(")
  !Save everything that is before function name into string1
  !(only if function name is not the first thing in string)
  IF( i>1 ) THEN
    string1 = TRIM(ADJUSTL(string(1:i-1)))
  ENDIF
  !Save position of the opening parenthesis
  p1 = i+3
  !Save everything after opening parenthesis into temp
  temp = TRIM(ADJUSTL(string(p1+1:)))
  !Find position of matching closing parenthesis
  m=0
  i=0
  j=1
  p2 = LEN_TRIM(string)
  DO WHILE( m==0 .AND. i<=LEN_TRIM(temp) )
    i=i+1
    IF( temp(i:i)=="(" ) THEN
      j=j+1
    ELSEIF( temp(i:i)==")" ) THEN
      j=j-1
    ENDIF
    IF( j==0 ) THEN
      !i is the position of the matching closing parenthesis
      p2=i
      m=1
      EXIT
    ENDIF
  ENDDO
  !Save everything after closing parenthesis into string2
  string2 = TRIM(ADJUSTL(temp(p2+1:)))
  !Remove from temp everything that is after closing parenthesis
  !NOTE: if no matching closing parenthesis was found, then temp
  !      will contain everything that is after opening parenthesis
  temp = temp(1:p2-1)
  !Interpret the expression that is inside the parenthesis, save result in x
  x = EXPREVAL(temp,recuri,status)
  !Perform the calculation of the cosinus
  eval = DCOS(x)
  !Eliminate round-off errors
  IF( DABS(eval)<=low_limit ) eval = 0.d0
  IF( DABS(1.d0-DABS(eval))<=low_limit ) eval = eval/DABS(eval)
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: cos("//TRIM(ADJUSTL(temp))//") = ", eval
  !Write the result into temp
  WRITE(temp,*) eval
  !Re-write the string by concatenating string1, temp, and string2
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  !Evaluate the resulting string
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( INDEX(string,"asin(")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED ARCSINUS: ", TRIM(string)
  !Save position of the expression "asin"
  i = INDEX(string,"asin(")
  !Save everything that is before function name into string1
  !(only if function name is not the first thing in string)
  IF( i>1 ) THEN
    string1 = TRIM(ADJUSTL(string(1:i-1)))
  ENDIF
  !Save position of the opening parenthesis
  p1 = i+4
  !Save everything after opening parenthesis into temp
  temp = TRIM(ADJUSTL(string(p1+1:)))
  !Find position of matching closing parenthesis
  m=0
  i=0
  j=1
  p2 = LEN_TRIM(string)
  DO WHILE( m==0 .AND. i<=LEN_TRIM(temp) )
    i=i+1
    IF( temp(i:i)=="(" ) THEN
      j=j+1
    ELSEIF( temp(i:i)==")" ) THEN
      j=j-1
    ENDIF
    IF( j==0 ) THEN
      !i is the position of the matching closing parenthesis
      p2=i
      m=1
      EXIT
    ENDIF
  ENDDO
  !Save everything after closing parenthesis into string2
  string2 = TRIM(ADJUSTL(temp(p2+1:)))
  !Remove from temp everything that is after closing parenthesis
  !NOTE: if no matching closing parenthesis was found, then temp
  !      will contain everything that is after opening parenthesis
  temp = temp(1:p2-1)
  !Interpret the expression that is inside the parenthesis, save result in x
  x = EXPREVAL(temp,recuri,status)
  !Perform the calculation of the arcsinus
  eval = DASIN(x)
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: asin("//TRIM(ADJUSTL(temp))//") = ", eval
  !Write the result into temp
  WRITE(temp,*) eval
  !Re-write the string by concatenating string1, temp, and string2
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  !Evaluate the resulting string
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( INDEX(string,"sin(")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED SINUS: ", TRIM(string)
  !Save position of the expression "sin"
  i = INDEX(string,"sin(")
  !Save everything that is before function name into string1
  !(only if function name is not the first thing in string)
  IF( i>1 ) THEN
    string1 = TRIM(ADJUSTL(string(1:i-1)))
  ENDIF
  !Save position of the opening parenthesis
  p1 = i+3
  !Save everything after opening parenthesis into temp
  temp = TRIM(ADJUSTL(string(p1+1:)))
  !Find position of matching closing parenthesis
  m=0
  i=0
  j=1
  p2 = LEN_TRIM(string)
  DO WHILE( m==0 .AND. i<=LEN_TRIM(temp) )
    i=i+1
    IF( temp(i:i)=="(" ) THEN
      j=j+1
    ELSEIF( temp(i:i)==")" ) THEN
      j=j-1
    ENDIF
    IF( j==0 ) THEN
      !i is the position of the matching closing parenthesis
      p2=i
      m=1
      EXIT
    ENDIF
  ENDDO
  !Save everything after closing parenthesis into string2
  string2 = TRIM(ADJUSTL(temp(p2+1:)))
  !Remove from temp everything that is after closing parenthesis
  !NOTE: if no matching closing parenthesis was found, then temp
  !      will contain everything that is after opening parenthesis
  temp = temp(1:p2-1)
  !Interpret the expression that is inside the parenthesis, save result in x
  x = EXPREVAL(temp,recuri,status)
  !Perform the calculation of the sinus
  eval = DSIN(x)
  !Eliminate round-off errors
  IF( DABS(eval)<=low_limit ) eval = 0.d0
  IF( DABS(1.d0-DABS(eval))<=low_limit ) eval = eval/DABS(eval)
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: sin("//TRIM(ADJUSTL(temp))//") = ", eval
  !Write the result into temp
  WRITE(temp,*) eval
  !Re-write the string by concatenating string1, temp, and string2
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  !Evaluate the resulting string
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( INDEX(string,"atan(")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED ARCTANGENT: ", TRIM(string)
  !Save position of the expression "atan"
  i = INDEX(string,"atan(")
  !Save everything that is before function name into string1
  !(only if function name is not the first thing in string)
  IF( i>1 ) THEN
    string1 = TRIM(ADJUSTL(string(1:i-1)))
  ENDIF
  !Save position of the opening parenthesis
  p1 = i+4
  !Save everything after opening parenthesis into temp
  temp = TRIM(ADJUSTL(string(p1+1:)))
  !Find position of matching closing parenthesis
  m=0
  i=0
  j=1
  p2 = LEN_TRIM(string)
  DO WHILE( m==0 .AND. i<=LEN_TRIM(temp) )
    i=i+1
    IF( temp(i:i)=="(" ) THEN
      j=j+1
    ELSEIF( temp(i:i)==")" ) THEN
      j=j-1
    ENDIF
    IF( j==0 ) THEN
      !i is the position of the matching closing parenthesis
      p2=i
      m=1
      EXIT
    ENDIF
  ENDDO
  !Save everything after closing parenthesis into string2
  string2 = TRIM(ADJUSTL(temp(p2+1:)))
  !Remove from temp everything that is after closing parenthesis
  !NOTE: if no matching closing parenthesis was found, then temp
  !      will contain everything that is after opening parenthesis
  temp = temp(1:p2-1)
  !Interpret the expression that is inside the parenthesis, save result in x
  x = EXPREVAL(temp,recuri,status)
  !Perform the calculation of the arctangent
  eval = DATAN(x)
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: atan("//TRIM(ADJUSTL(temp))//") = ", eval
  !Write the result into temp
  WRITE(temp,*) eval
  !Re-write the string by concatenating string1, temp, and string2
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  !Evaluate the resulting string
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( INDEX(string,"atan2(")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED ARCTANGENT2: ", TRIM(string)
  !Save position of the expression "atan"
  i = INDEX(string,"atan2(")
  !Save everything that is before function name into string1
  !(only if function name is not the first thing in string)
  IF( i>1 ) THEN
    string1 = TRIM(ADJUSTL(string(1:i-1)))
  ENDIF
  !Save position of the opening parenthesis
  p1 = i+5
  !Save everything after opening parenthesis into temp
  temp = TRIM(ADJUSTL(string(p1+1:)))
  !Find position of matching closing parenthesis
  m=0
  i=0
  j=1
  p2 = LEN_TRIM(string)
  DO WHILE( m==0 .AND. i<=LEN_TRIM(temp) )
    i=i+1
    IF( temp(i:i)=="(" ) THEN
      j=j+1
    ELSEIF( temp(i:i)==")" ) THEN
      j=j-1
    ENDIF
    IF( j==0 ) THEN
      !i is the position of the matching closing parenthesis
      p2=i
      m=1
      EXIT
    ENDIF
  ENDDO
  !Save everything after closing parenthesis into string2
  string2 = TRIM(ADJUSTL(temp(p2+1:)))
  !Remove from temp everything that is after closing parenthesis
  !NOTE: if no matching closing parenthesis was found, then temp
  !      will contain everything that is after opening parenthesis
  temp = temp(1:p2-1)
  !Get position of the comma
  i = SCAN(temp,",")
  IF( i>1 ) THEN
    !Save numerator into temp1
    temp1 = temp(1:i-1)
    !Interpret this expression, save result in y
    y = EXPREVAL(temp1,recuri,status)
    !Save denominator into temp2
    temp2 = temp(i+1:)
    !Interpret this expression, save result in x
    x = EXPREVAL(temp2,recuri,status)
    !Perform the calculation of the arctangent
    !IF( DABS(x)>1.d-12 ) THEN
      eval = DATAN2(y,x)
    !ELSE
    !  IF(verbosity==4) PRINT*, " X ! X ERROR: division by zero: ", TRIM(string)
    !  status=10
    !  RETURN
    !ENDIF
  ELSE
    !No division: cannot compute an arctan2
    !Perform the calculation of a regular arctangent
    !Interpret the expression that is inside the parenthesis, save result in x
    x = EXPREVAL(temp,recuri,status)
    !Perform the calculation of the arctangent
    eval = DATAN(x)
  ENDIF
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: atan2("//TRIM(ADJUSTL(temp))//") = ", eval
  !Write the result into temp
  WRITE(temp,*) eval
  !Re-write the string by concatenating string1, temp, and string2
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  !Evaluate the resulting string
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( INDEX(string,"tan(")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED TANGENT: ", TRIM(string)
  !Save position of the expression "tan"
  i = INDEX(string,"tan(")
  !Save everything that is before function name into string1
  !(only if function name is not the first thing in string)
  IF( i>1 ) THEN
    string1 = TRIM(ADJUSTL(string(1:i-1)))
  ENDIF
  !Save position of the opening parenthesis
  p1 = i+3
  !Save everything after opening parenthesis into temp
  temp = TRIM(ADJUSTL(string(p1+1:)))
  !Find position of matching closing parenthesis
  m=0
  i=0
  j=1
  p2 = LEN_TRIM(string)
  DO WHILE( m==0 .AND. i<=LEN_TRIM(temp) )
    i=i+1
    IF( temp(i:i)=="(" ) THEN
      j=j+1
    ELSEIF( temp(i:i)==")" ) THEN
      j=j-1
    ENDIF
    IF( j==0 ) THEN
      !i is the position of the matching closing parenthesis
      p2=i
      m=1
      EXIT
    ENDIF
  ENDDO
  !Save everything after closing parenthesis into string2
  string2 = TRIM(ADJUSTL(temp(p2+1:)))
  !Remove from temp everything that is after closing parenthesis
  !NOTE: if no matching closing parenthesis was found, then temp
  !      will contain everything that is after opening parenthesis
  temp = temp(1:p2-1)
  !Interpret the expression that is inside the parenthesis, save result in x
  x = EXPREVAL(temp,recuri,status)
  !Perform the calculation of the tangent
  eval = DTAN(x)
  !Eliminate round-off errors
  IF( DABS(eval)<=low_limit ) eval = 0.d0
  IF( DABS(1.d0-DABS(eval))<=low_limit ) eval = eval/DABS(eval)
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: tan("//TRIM(ADJUSTL(temp))//") = ", eval
  !Write the result into temp
  WRITE(temp,*) eval
  !Re-write the string by concatenating string1, temp, and string2
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  !Evaluate the resulting string
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( INDEX(string,"exp(")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED EXPONENTIAL: ", TRIM(string)
  !Save position of the expression "exp"
  i = INDEX(string,"exp(")
  !Save everything that is before function name into string1
  !(only if function name is not the first thing in string)
  IF( i>1 ) THEN
    string1 = TRIM(ADJUSTL(string(1:i-1)))
  ENDIF
  !Save position of the opening parenthesis
  p1 = i+3
  !Save everything after opening parenthesis into temp
  temp = TRIM(ADJUSTL(string(p1+1:)))
  !Find position of matching closing parenthesis
  m=0
  i=0
  j=1
  p2 = LEN_TRIM(string)
  DO WHILE( m==0 .AND. i<=LEN_TRIM(temp) )
    i=i+1
    IF( temp(i:i)=="(" ) THEN
      j=j+1
    ELSEIF( temp(i:i)==")" ) THEN
      j=j-1
    ENDIF
    IF( j==0 ) THEN
      !i is the position of the matching closing parenthesis
      p2=i
      m=1
      EXIT
    ENDIF
  ENDDO
  !Save everything after closing parenthesis into string2
  string2 = TRIM(ADJUSTL(temp(p2+1:)))
  !Remove from temp everything that is after closing parenthesis
  !NOTE: if no matching closing parenthesis was found, then temp
  !      will contain everything that is after opening parenthesis
  temp = temp(1:p2-1)
  !Interpret the expression that is inside the parenthesis, save result in x
  x = EXPREVAL(temp,recuri,status)
  !Perform the calculation of the exponential
  eval = DEXP(x)
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: exp("//TRIM(ADJUSTL(temp))//") = ", eval
  !Write the result into temp
  WRITE(temp,*) eval
  !Re-write the string by concatenating string1, temp, and string2
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  !Evaluate the resulting string
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( INDEX(string,"ln(")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED LN: ", TRIM(string)
  !Save position of the expression "exp"
  i = INDEX(string,"ln(")
  !Save everything that is before function name into string1
  !(only if function name is not the first thing in string)
  IF( i>1 ) THEN
    string1 = TRIM(ADJUSTL(string(1:i-1)))
  ENDIF
  !Save position of the opening parenthesis
  p1 = i+2
  !Save everything after opening parenthesis into temp
  temp = TRIM(ADJUSTL(string(p1+1:)))
  !Find position of matching closing parenthesis
  m=0
  i=0
  j=1
  p2 = LEN_TRIM(string)
  string2 = ""
  DO WHILE( m==0 .AND. i<=LEN_TRIM(temp) )
    i=i+1
    IF( temp(i:i)=="(" ) THEN
      j=j+1
    ELSEIF( temp(i:i)==")" ) THEN
      j=j-1
    ENDIF
    IF( j==0 ) THEN
      !i is the position of the matching closing parenthesis
      p2=i
      m=1
      EXIT
    ENDIF
  ENDDO
  !Save everything after closing parenthesis into string2
  string2 = TRIM(ADJUSTL(temp(p2+1:)))
  !Remove from temp everything that is after closing parenthesis
  !NOTE: if no matching closing parenthesis was found, then temp
  !      will contain everything that is after opening parenthesis
  temp = temp(1:p2-1)
  !Interpret the expression that is inside the parenthesis, save result in x
  x = EXPREVAL(temp,recuri,status)
  !Perform the calculation of the natural logarithm
  IF( x>0.d0 ) THEN
    eval = DLOG(x)
  ELSE
    IF(verbosity==4) PRINT*, "X ! X ERROR: cannot compute logarithm of negative number"
    status=1
    GOTO 900
  ENDIF
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: ln("//TRIM(ADJUSTL(temp))//") = ", eval
  !Write the result into temp
  WRITE(temp,*) eval
  !Re-write the string by concatenating string1, temp, and string2
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  !Evaluate the resulting string
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( INDEX(string,"log(")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED LOG: ", TRIM(string)
  !Save position of the expression "exp"
  i = INDEX(string,"log(")
  !Save everything that is before function name into string1
  !(only if function name is not the first thing in string)
  IF( i>1 ) THEN
    string1 = TRIM(ADJUSTL(string(1:i-1)))
  ENDIF
  !Save position of the opening parenthesis
  p1 = i+3
  !Save everything after opening parenthesis into temp
  temp = TRIM(ADJUSTL(string(p1+1:)))
  !Find position of matching closing parenthesis
  m=0
  i=0
  j=1
  p2 = LEN_TRIM(string)
  string2 = ""
  DO WHILE( m==0 .AND. i<=LEN_TRIM(temp) )
    i=i+1
    IF( temp(i:i)=="(" ) THEN
      j=j+1
    ELSEIF( temp(i:i)==")" ) THEN
      j=j-1
    ENDIF
    IF( j==0 ) THEN
      !i is the position of the matching closing parenthesis
      p2=i
      m=1
      EXIT
    ENDIF
  ENDDO
  !Save everything after closing parenthesis into string2
  string2 = TRIM(ADJUSTL(temp(p2+1:)))
  !Remove from temp everything that is after closing parenthesis
  !NOTE: if no matching closing parenthesis was found, then temp
  !      will contain everything that is after opening parenthesis
  temp = temp(1:p2-1)
  !Interpret the expression that is inside the parenthesis, save result in x
  x = EXPREVAL(temp,recuri,status)
  !Perform the calculation of the logarithm
  IF( x>0.d0 ) THEN
    eval = DLOG10(x)
  ELSE
    IF(verbosity==4) PRINT*, "X ! X ERROR: cannot compute logarithm of negative number"
    status=1
    GOTO 900
  ENDIF
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: ln("//TRIM(ADJUSTL(temp))//") = ", eval
  !Write the result into temp
  WRITE(temp,*) eval
  !Re-write the string by concatenating string1, temp, and string2
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  !Evaluate the resulting string
  eval = EXPREVAL(string,recuri,status)
  !
!If it is not a function but there is a parenthesis, interpret the parenthesis
ELSE IF( p1>0 ) THEN
  !An opening parenthesis was detected at the beginning of this routine
  !Save everything that is before this opening parenthesis into string1
  !(only if parenthesis is not the first thing in string)
  IF( p1>1 ) THEN
    string1 = TRIM(ADJUSTL(string(1:p1-1)))
  ENDIF
  !Find position of matching closing parenthesis
  m=0
  i=p1
  j=1
  DO WHILE(m==0 .AND. i<LEN_TRIM(string) )
    i=i+1
    IF( string(i:i)=="(" ) THEN
      j=j+1
    ELSEIF( string(i:i)==")" ) THEN
      j=j-1
    ENDIF
    IF( j==0 ) THEN
      !i is the position of the matching closing parenthesis
      p2=i
      m=1
      EXIT
    ENDIF
  ENDDO
  !Save everything after closing parenthesis into string2
  IF( m>0 ) THEN
    string2 = TRIM(ADJUSTL(string(p2+1:)))
  ELSE
    string2 = ""
    p2=LEN_TRIM(string)+1
  ENDIF
  !Check if this parenthesis has an exponent or factorial attached after it
  IF( string2(1:1)=="^" ) THEN
    temp1 = string(p1+1:p2-1)
    !Get exponent
    i = SCAN(string2,"+*/:")
    j = SCAN(string2(2:),"-")+1
    temp = string(p1+1:p2+i-1)
    IF( i>0 ) THEN
      temp2 = TRIM(ADJUSTL(string2(2:i-1)))
      string2 = TRIM(ADJUSTL(string2(i:)))
    ELSE IF( j>1 ) THEN
      temp2 = string2(2:i-1)
      string2 = TRIM(ADJUSTL(string2(j:)))
    ELSE
      temp2 = string2(2:)
      string2 = TRIM(ADJUSTL(string2(3:)))
    ENDIF
    IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED PARENTHESIS w/ EXPONENT: ", TRIM(temp)
    !Evaluate expression inside parenthesis
    x = EXPREVAL(temp1,recuri,status)
    IF(verbosity==4) PRINT*, "                   ", TRIM(temp2)
    !Evaluate exponent
    y = EXPREVAL(temp2,recuri,status)
    eval = x**y
    IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: ("//TRIM(ADJUSTL(temp1))//&
                   & ") ^ "//TRIM(ADJUSTL(temp2))//"  = ", eval
  ELSEIF( string2(1:1)=="!" ) THEN
    temp1 = string(p1+1:p2-1)
    temp = string(p1:p2+1)
    string2 = TRIM(ADJUSTL(string2(2:)))
    IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED PARENTHESIS w/ FACTORIAL: ", TRIM(temp)
    !Evaluate expression inside parenthesis
    x = EXPREVAL(temp1,recuri,status)
    !Evaluate factorial
    eval = 1.d0
    DO i=1,NINT(x)
      eval = eval * DBLE(i)
    ENDDO
    IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: ("//TRIM(ADJUSTL(temp1))//&
                   & ")! = ", eval
  ELSE
    temp = string(p1+1:p2-1)
    IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED PARENTHESIS: ", TRIM(temp)
    !Evaluate expression inside parenthesis
    eval = EXPREVAL(temp,recuri,status)
    IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: ("//TRIM(ADJUSTL(temp1))//&
                   & ") = ", eval
  ENDIF
  !Write the result into temp
  WRITE(temp,*) eval
  !Re-write the string by concatenating string1, temp, and string2
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  !Evaluate the resulting string
  eval = EXPREVAL(string,recuri,status)
  !
!Replace divisions by their eval
ELSE IF( SCAN(string,"/")>0 .OR. SCAN(string,":")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED DIVISION: ", TRIM(string)
  !Save position of the division operator
  i = SCAN(string,"/")
  IF( i==0 ) THEN
    i = SCAN(string,":")
  ENDIF
  p1 = i
  temp = string(1:i-1)
  j = SCAN(TRIM(temp),operators,BACK=.TRUE.)
  IF( string(j-1:j-1)=="E" .OR. string(j-1:j-1)=="e" ) THEN
    j = SCAN(TRIM(temp(1:j-1)),operators,BACK=.TRUE.)
  ENDIF
  temp1 = TRIM(ADJUSTL(temp(j+1:)))
  string1 = string(1:j)
  IF(verbosity==4) PRINT*, "          ", TRIM(temp1)
  temp = string(i+1:)
  j = SCAN(TRIM(ADJUSTL(temp)),operators)
  IF( temp(1:1)=="-" ) THEN
    j = SCAN(TRIM(ADJUSTL(temp(2:))),operators)
    IF( temp(j:j)=="E" .OR. temp(j:j)=="e" ) THEN
      j = SCAN(TRIM(temp(j+2:)),operators)
    ENDIF
    IF(j>0) THEN
      j=j+1
      temp2 = TRIM(ADJUSTL(temp(:j-1)))
      string2 = string(p1+j:)
    ELSE
      temp2 = TRIM(ADJUSTL(temp))
      string2 = ""
    ENDIF
  ELSE IF(j>1) THEN
    IF( temp(j-1:j-1)=="E" .OR. temp(j-1:j-1)=="e" ) THEN
      j = SCAN(TRIM(temp(j+1:)),operators)
    ENDIF
    IF(j>0) THEN
      temp2 = TRIM(ADJUSTL(temp(:j-1)))
      string2 = string(p1+j:)
    ELSE
      temp2 = TRIM(ADJUSTL(temp))
      string2 = ""
    ENDIF
  ELSE
    temp2 = TRIM(ADJUSTL(temp))
    string2 = ""
  ENDIF
  IF(verbosity==4) PRINT*, "          ", TRIM(temp2)
  x = EXPREVAL(temp1,recuri,status)
  y = EXPREVAL(temp2,recuri,status)
  IF( DABS(y)<low_limit ) THEN
    IF(verbosity==4) PRINT*, " X ! X ERROR: division by zero: ", TRIM(string)
    status=10
    RETURN
  ELSE
    eval = x/y
  ENDIF
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: "//TRIM(ADJUSTL(temp1))//" / "//TRIM(ADJUSTL(temp2))//"  = ", eval
  WRITE(temp,*) eval
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  eval = EXPREVAL(string,recuri,status)
  !
!Replace modulo by their eval
ELSE IF( SCAN(string,"%")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED MODULO: ", TRIM(string)
  !Save position of the modulo operator
  i = SCAN(string,"%")
  p1 = i
  temp = string(1:i-1)
  j = SCAN(TRIM(temp),operators,BACK=.TRUE.)
  IF( string(j-1:j-1)=="E" .OR. string(j-1:j-1)=="e" ) THEN
    j = SCAN(TRIM(temp(1:j-1)),operators,BACK=.TRUE.)
  ENDIF
  temp1 = TRIM(ADJUSTL(temp(j+1:)))
  string1 = string(1:j)
  IF(verbosity==4) PRINT*, "          ", TRIM(temp1)
  temp = string(i+1:)
  j = SCAN(TRIM(ADJUSTL(temp)),operators)
  IF( temp(1:1)=="-" ) THEN
    j = SCAN(TRIM(ADJUSTL(temp(2:))),operators)
    IF( temp(j:j)=="E" .OR. temp(j:j)=="e" ) THEN
      j = SCAN(TRIM(temp(j+2:)),operators)
    ENDIF
    IF(j>0) THEN
      j=j+1
      temp2 = TRIM(ADJUSTL(temp(:j-1)))
      string2 = string(p1+j:)
    ELSE
      temp2 = TRIM(ADJUSTL(temp))
      string2 = ""
    ENDIF
  ELSE IF(j>1) THEN
    IF( temp(j-1:j-1)=="E" .OR. temp(j-1:j-1)=="e" ) THEN
      j = SCAN(TRIM(temp(j+1:)),operators)
    ENDIF
    IF(j>0) THEN
      temp2 = TRIM(ADJUSTL(temp(:j-1)))
      string2 = string(p1+j:)
    ELSE
      temp2 = TRIM(ADJUSTL(temp))
      string2 = ""
    ENDIF
  ELSE
    temp2 = TRIM(ADJUSTL(temp))
    string2 = ""
  ENDIF
  IF(verbosity==4) PRINT*, "          ", TRIM(temp2)
  x = EXPREVAL(temp1,recuri,status)
  y = EXPREVAL(temp2,recuri,status)
  eval = MOD(x,y)
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: "//TRIM(ADJUSTL(temp1))//" % "//TRIM(ADJUSTL(temp2))//"  = ", eval
  WRITE(temp,*) eval
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  eval = EXPREVAL(string,recuri,status)
  !
!Replace exponents by their eval
ELSE IF( SCAN(string,"^")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED EXPONENT: ", TRIM(string)
  i = SCAN(string,"^")
  p1 = i
  IF( string(i-1:i-1)==')' ) THEN
    temp = string(1:i-2)
    !Search for matching opening parenthesis
    m=0
    i=i-1
    j=1
    p2=0
    DO WHILE(m==0 .AND. i>0)
      i=i-1
      IF( string(i:i)=="(" ) THEN
        j=j-1
      ELSEIF( string(i:i)==")" ) THEN
        j=j+1
      ENDIF
      IF( j==0 ) THEN
        !i is the position of the matching opening parenthesis
        p2=i
        m=1
        EXIT
      ENDIF
    ENDDO
    IF( i<1 .OR. m==0 ) THEN
      p2=0
      string1 = ""
    ELSE
      string1 = string(1:i-1)
    ENDIF
    IF( p2<i-2 ) THEN
      temp1 = string(p2+1:i-2)
    ELSE
      temp1 = ""
    ENDIF
    IF(verbosity==4) PRINT*, "          ", TRIM(temp1)
  ELSE
    temp = string(1:i-1)
    j = SCAN(TRIM(temp),"+-*/:^)(",BACK=.TRUE.)
    IF( SCAN(TRIM(temp),"-",BACK=.TRUE.)<=1 ) THEN
      temp1 = TRIM(ADJUSTL(temp(1:)))
    ELSE
      temp1 = TRIM(ADJUSTL(temp(j+1:)))
      DO WHILE( temp(j-1:j-1)=="E" .OR. temp(j-1:j-1)=="e" )
        !Nevermind: it is an exponent
        j = SCAN(temp(:j-1),"-",BACK=.TRUE.)
        temp1 = TRIM(ADJUSTL(temp(1:)))
      ENDDO
    ENDIF
    string1 = string(1:j)
    IF(verbosity==4) PRINT*, "          ", TRIM(temp1)
  ENDIF
  temp = string(p1+1:)
  j = SCAN(TRIM(ADJUSTL(temp)),"*/:^+")
  IF(j>0) THEN
    temp2 = TRIM(ADJUSTL(temp(:j-1)))
    string2 = string(p1+j:)
  ELSE
    temp2 = TRIM(ADJUSTL(temp))
    string2 = ""
  ENDIF
  IF(verbosity==4) PRINT*, "          ", TRIM(temp2)
  x = EXPREVAL(temp1,recuri,status)
  y = EXPREVAL(temp2,recuri,status)
  eval = x**y
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: "//TRIM(ADJUSTL(temp1))//" ^ "//TRIM(ADJUSTL(temp2))//"  = ", eval
  WRITE(temp,*) eval
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  eval = EXPREVAL(string,recuri,status)
  !
ELSE IF( INDEX(string,"!")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED FACTORIAL: ", TRIM(string)
  i = SCAN(string,"!")
  p1 = i
  string1 = ""
  !Save everything after exclamation mark into string2
  string2 = TRIM(ADJUSTL(string(i+1:)))
  !Save everything before exclamation mark into temp
  temp = TRIM(ADJUSTL(string(1:i-1)))
  temp1 = temp
  !Check if exclamation mark is located after a closing parenthesis
  IF( string(i-1:i-1)==')' ) THEN
    !Yes it is: keep into temp only what is before the closing parenthesis
    temp = string(1:i-2)
    !Search for matching opening parenthesis
    m=0
    i=p1
    j=1
    DO WHILE(m==0 .AND. i>0)
      i=i-1
      IF( string(i:i)=="(" ) THEN
        j=j-1
      ELSEIF( string(i:i)==")" ) THEN
        j=j+1
      ENDIF
      IF( j==0 ) THEN
        !i is the position of the matching opening parenthesis
        p1=i
        m=1
        EXIT
      ENDIF
    ENDDO
    IF( i<=1 ) p2=1
    temp1 = string(p2+1:p1-2)
    string1 = TRIM(ADJUSTL(string(1:p1-2)))
    IF(verbosity==4) PRINT*, "          ", TRIM(temp1)
  ELSE
    !There is no closing parenthesis
    !Search for nearest operator before the exclamation mark
    j = SCAN(temp,"+-*:/(",BACK=.TRUE.)
    p1 = j
    IF( j>0 ) THEN
      temp1 = temp(j+1:)
      string1 = TRIM(ADJUSTL(string(1:j-1)))
      temp = ""
    ELSE
      temp1 = TRIM(temp)
    ENDIF
  ENDIF
  !Interpret the expression that is before the exclamation mark, save result in x
  x = EXPREVAL(temp1,recuri,status)
  !Perform the calculation of the factorial
  IF( NINT(x)<0 ) THEN
    IF(verbosity==4) PRINT*, "X ! X ERROR: cannot compute factorial of negative number"
    status=1
    GOTO 900
  ENDIF
  eval = 1.d0
  DO i=1,NINT(x)
    eval = eval * DBLE(i)
  ENDDO
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: "//TRIM(ADJUSTL(temp1))//"! = ", eval
  !Write the result into temp
  WRITE(temp,*) eval
  !Re-write the string by concatenating string1, temp, and string2
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  !Evaluate the resulting string
  eval = EXPREVAL(string,recuri,status)
  !
!Replace products by their eval
ELSE IF( SCAN(string,"*")>0 ) THEN
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED MULTIPLICATION: ", TRIM(string)
  i = SCAN(string,"*")
  p1 = i
  temp = string(1:i-1)
  j = SCAN(TRIM(temp),operators,BACK=.TRUE.)
  IF( string(j-1:j-1)=="E" .OR. string(j-1:j-1)=="e" ) THEN
    j = SCAN(TRIM(temp(1:j-1)),operators,BACK=.TRUE.)
  ENDIF
  temp1 = TRIM(ADJUSTL(temp(j+1:)))
  string1 = string(1:j)
  IF(verbosity==4) PRINT*, "          ", TRIM(temp1)
  temp = string(i+1:)
  j = SCAN(TRIM(ADJUSTL(temp)),operators)
  IF( temp(1:1)=="-" ) THEN
    j = SCAN(TRIM(ADJUSTL(temp(2:))),operators)
    IF( temp(j:j)=="E" .OR. temp(j:j)=="e" ) THEN
      j = SCAN(TRIM(temp(j+2:)),operators)
    ENDIF
    IF(j>0) THEN
      j=j+1
      temp2 = TRIM(ADJUSTL(temp(:j-1)))
      string2 = string(p1+j:)
    ELSE
      temp2 = TRIM(ADJUSTL(temp))
      string2 = ""
    ENDIF
  ELSE IF(j>1) THEN
    IF( temp(j-1:j-1)=="E" .OR. temp(j-1:j-1)=="e" ) THEN
      j = SCAN(TRIM(temp(j+1:)),operators)
    ENDIF
    IF(j>0) THEN
      temp2 = TRIM(ADJUSTL(temp(:j-1)))
      string2 = string(p1+j:)
    ELSE
      temp2 = TRIM(ADJUSTL(temp))
      string2 = ""
    ENDIF
  ELSE
    temp2 = TRIM(ADJUSTL(temp))
    string2 = ""
  ENDIF
  IF(verbosity==4) PRINT*, "          ", TRIM(temp2)
  IF( LEN_TRIM(temp1)>0 ) THEN
    x = EXPREVAL(temp1,recuri,status)
  ELSE
    x = 1.d0
  ENDIF
  IF( LEN_TRIM(temp2)>0 ) THEN
    y = EXPREVAL(temp2,recuri,status)
  ELSE
    y = 1.d0
  ENDIF
  eval = x*y
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: "//TRIM(ADJUSTL(temp1))//" * "//TRIM(ADJUSTL(temp2))//"  = ", eval
  WRITE(temp,*) eval
  string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
  eval = EXPREVAL(string,recuri,status)
  !
!Reject all other characters
ELSE IF( SCAN(string,"ABCDFGHIJKLMNOPQRSTUVWXYZabcdfghijklmnopqrstuvwxyz?,;&éè@àçè|={}")>0 ) THEN
  !IF(verbosity==4) PRINT*, " X ! X ERROR: illegal character in input string: ", TRIM(string)
  status = 1
  GOTO 900
  !
  !
!Replace sums by their eval
ELSE IF( SCAN(string,"+")>0 ) THEN
  i = SCAN(string,"+",BACK=.TRUE.)
  IF( i==1 ) THEN
    !It is a positive number located at the beginning of the string
    !Return its eval
    IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED NUMBER(+): ", TRIM(string)
    string = ADJUSTL(string(2:))
    eval = EXPREVAL(string,recuri,status)
    !
  ELSEIF( string(i-1:i-1)=="E" .OR. string(i-1:i-1)=="e" ) THEN
    !It is a negative number located at the beginning of the string
    !return its eval
    IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED NUMBER(E+): ", TRIM(string)
    READ(string,*,END=900,ERR=900) eval
    !
  ELSE
    IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED ADDITION: ", TRIM(string)
    temp1 = string(:i-1)
    temp2 = string(i+1:)
    x = EXPREVAL(temp1,recuri,status)
    y = EXPREVAL(temp2,recuri,status)
    eval = x+y
    IF(verbosity==4) PRINT*,  TRIM(prefix)//"  RESULT: "//TRIM(ADJUSTL(temp1))//" + "//TRIM(ADJUSTL(temp2))//"  = ", eval
    WRITE(temp,*) eval
!     string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
!     CALL EXPREVAL(string,eval,recuri,status)
  ENDIF
  !
  !
!Replace subtractions by their eval
! ELSE IF( ms>0 ) THEN
!   IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED SUBTRACTION: ", TRIM(string)
!   temp1 = string(:ms-1)
!   temp2 = string(ms+1:)
!   CALL EXPREVAL(temp1,x,recuri,status)
!   CALL EXPREVAL(temp2,y,recuri,status)
!   eval = x-y
!   IF(verbosity==4) PRINT*,  TRIM(prefix)//"  RESULT: "//TRIM(ADJUSTL(temp1))//" - "//TRIM(ADJUSTL(temp2))//"  = ", eval
!   WRITE(temp,*) eval
  !
  !
!Replace subtractions by their eval
ELSE IF( SCAN(string,"-",BACK=.TRUE.)>0 ) THEN
  i = SCAN(string,"-",BACK=.TRUE.)
  IF( i==1 ) THEN
    !It is a negative number located at the beginning of the string
    !return its eval
    IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED NUMBER(-): ", TRIM(string)
    READ(string,*,END=900,ERR=900) eval
    !
  ELSEIF( string(i-1:i-1)=="E" .OR. string(i-1:i-1)=="e" ) THEN
    !It is a negative number located at the beginning of the string
    !return its eval
    IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED NUMBER(E-): ", TRIM(string)
    READ(string,*,END=900,ERR=900) eval
    !
  ELSEIF( SCAN(string(i-1:i-1),"*+:")>0 ) THEN
    !It is a negative number located after another operator
    !Return its eval
    p1 = i
    temp = ADJUSTL(string(i:))
    j = SCAN(TRIM(temp),operators)
    IF(j>0) THEN
      temp2 = TRIM(ADJUSTL(temp(:j-1)))
      string2 = string(p1+j:)
    ELSE
      temp2 = TRIM(ADJUSTL(temp))
      string2 = ""
    ENDIF
    IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED NEGATIVE NUMBER: ", TRIM(temp2)
    READ(temp2,*,END=900,ERR=900) eval
    !
  ELSE
    !It is an actual subtraction
    IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED SUBTRACTION: ", TRIM(string)
    p1 = i
    temp1 = string(:i-1)
    temp2 = string(i+1:)
    x = EXPREVAL(temp1,recuri,status)
    y = EXPREVAL(temp2,recuri,status)
    eval = x-y
    IF(verbosity==4) PRINT*, TRIM(prefix)//"  RESULT: "//TRIM(ADJUSTL(temp1))//" - "//TRIM(ADJUSTL(temp2))//"  = ", eval
    WRITE(temp,*) eval
!     string = TRIM(string1)//TRIM(ADJUSTL(temp))//TRIM(ADJUSTL(string2))
!     CALL EXPREVAL(string,eval,recuri,status)
  ENDIF
  !
ELSE
  IF(verbosity==4) PRINT*, TRIM(prefix)//"  DETECTED NUMBER(F): ", TRIM(string)
  READ(string,*,END=900,ERR=900) eval
ENDIF
!
!IF( DABS(eval)<=low_limit ) eval = 0.d0
IF(verbosity==4) PRINT*, TRIM(prefix)//"  END CALL TO EXPREVAL, RETURN VALUE = ", eval
recuri = recuri-1
RETURN
!
!
900 CONTINUE
!CALL ATOMSK_MSG(2813,(/string/),(/0.d0/))
status = 1
RETURN
!
!
END FUNCTION EXPREVAL
!
!
!
END MODULE exprev
