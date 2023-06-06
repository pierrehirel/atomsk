MODULE subroutines
!
!**********************************************************************************
!*  SUBROUTINES                                                                   *
!**********************************************************************************
!* This module contains general subroutines dealing                               *
!* with files, manipulation of strings, arrays and so on.                         *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel (unless specified otherwise, cf each subroutine)  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 06 June 2023                                     *
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
!* List of subroutines in this module:                                            *
!* CHECKMEM            checks that local computer has enough memory               *
!* CHECK_ARRAY_CONSISTENCY checks that arrays P, S, AUX, AUXNAMES are consistent  *
!* STR_RMSPACE         removes all blank spaces in a string                       *
!* STR_CHAR2SPACE      replaces a character by blank space in a string            *
!* STR2BOOL            transforms a string into a boolean value                   *
!* INT2MONTH           transforms an integer into a month                         *
!* INT2DAY             transforms an integer into a day                           *
!* BOX2DBLE            transforms a string into a real number                     *
!* CHECKNAN            checks an array for NaN (Not a Number) values              *
!* DO_STATS            do some statistics on a 1-dimensional array                *
!* CHECK_ORTHOVEC      checks if two vectors are normal to each other             *
!* CART2FRAC           converts cartesian coordinates to fractional               *
!* FRAC2CART           converts fractional coordinates to cartesian               *
!* COUNT_OUTBOX        given a list of positions, counts how many are out of box  *
!* FIND_NSP            find the number of different entries in an array           *
!* FIND_IF_REDUCED     determines if values of an array are within 0 and 1        *
!* UNWRAP              unwrap atoms that have jumped a boundary                   *
!* LIST_SELECTED_ATOMS makes a list of selected atoms (chem.species/number)       *
!**********************************************************************************
!
!
USE comv
USE constants
USE functions
USE math
!
IMPLICIT NONE
!
!
CONTAINS
!
!
!
!********************************************************
! CHECKMEM
! This subroutine verifies that the given size
! is within the limit of signed integers, and that
! the current computer has enough memory to allocate
! a real array of size Nx4
!********************************************************
SUBROUTINE CHECKMEM(Asize,status)
!
IMPLICIT NONE
REAL(dp),INTENT(IN):: Asize !size of array to allocate
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Atest !test array
INTEGER,INTENT(OUT):: status !=0 if ok, =1 if Asize exceeds integer limit, =2 
                             !1=problem with S; 2=problem with AUX; 3=problem with AUXNAMES
INTEGER:: i
!
status = 0
!
IF( Asize<=0 .OR. Asize>NATOMS_MAX ) THEN
  ! Asize is negative or greater than 2 billions (limit of signed integers)
  ! => It will be impossible to allocate an array that large
  status = 1
ELSE
  !Asize is ok, but does current computer have enough memory?
  !Try to allocate and check if it works
  ALLOCATE( Atest(NINT(Asize),4),STAT=i)
  IF( i>0 ) THEN
    ! Allocation failed (not enough memory)
    status = 2
  ENDIF
  IF(ALLOCATED(Atest)) DEALLOCATE(Atest)
ENDIF
!
END SUBROUTINE CHECKMEM
!
!
!
!********************************************************
! CHECK_ARRAY_CONSISTENCY
! This subroutine verifies that the arrays S (containing
! the positions of shells), AUX (auxiliary properties),
! and AUXNAMES (names of auxiliary properties), have
! sizes that are consistent with the array P (atom positions)
! and with each other.
! Rules: if allocated, S must have same size MxN as P.
!        if AUXNAMES is allocated and has a size L,
!        then AUX must have a size MxL.
!********************************************************
SUBROUTINE CHECK_ARRAY_CONSISTENCY(P,S,AUX,AUXNAMES,status)
!
IMPLICIT NONE
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P !positions of atoms (or ionic cores)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S   !positions of ionic shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
INTEGER,INTENT(OUT):: status !=0 if no problem was detected
                             !1=problem with S; 2=problem with AUX; 3=problem with AUXNAMES
INTEGER:: i
!
status = 0
!
IF(ALLOCATED(S)) THEN
  IF( SIZE(S,1) > SIZE(P,1) ) THEN
    !There cannot be more shells than cores => exit
    status = 1
  ELSEIF( SIZE(S,1)<=0 ) THEN
    !no shell => S should not be allocated
    DEALLOCATE(S)
  ENDIF
ENDIF
IF( ALLOCATED(AUXNAMES) .OR. ALLOCATED(AUX) ) THEN
  !first dimension of AUX must correspond to number of atoms
  IF( SIZE(AUX,1) .NE. SIZE(P,1) ) THEN
    status = 2
  ENDIF
  !second dimension of AUX must correspond to number of auxiliary properties
  IF( SIZE(AUXNAMES) .NE. SIZE(AUX,2) ) THEN
    status = 3
  ENDIF
  !avoid arrays allocated with zero size
  IF( SIZE(AUXNAMES)<=0 .OR. SIZE(AUX,1)<=0 .OR. SIZE(AUX,2)<=0 ) THEN
    IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
    IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
  ENDIF
  !Verify that names of auxiliary properties are left-aligned
  IF( ALLOCATED(AUXNAMES) .AND. SIZE(AUXNAMES)>0 ) THEN
    DO i=1,SIZE(AUXNAMES)
      AUXNAMES(i) = TRIM(ADJUSTL(AUXNAMES(i)))
    ENDDO
  ENDIF
ENDIF
!
END SUBROUTINE CHECK_ARRAY_CONSISTENCY
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
! INT2MONTH
! Converts an integer to strings containing the
! month in long and short forms.
! Integer should be between 1 and 12
! (otherwise "December" is returned).
!********************************************************
SUBROUTINE INT2MONTH(number,month,smonth)
CHARACTER(LEN=16),INTENT(OUT):: month, smonth
INTEGER,INTENT(IN):: number
!
SELECT CASE(number)
CASE(1)
  month = "January"
  smonth = "Jan."
CASE(2)
  month = "February"
  smonth = "Feb."
CASE(3)
  month = "March"
  smonth = "March"
CASE(4)
  month = "April"
  smonth = "Apr."
CASE(5)
  month = "May"
  smonth = "May"
CASE(6)
  month = "June"
  smonth = "June"
CASE(7)
  month = "July"
  smonth = "July"
CASE(8)
  month = "August"
  smonth = "Aug."
CASE(9)
  month = "September"
  smonth = "Sept."
CASE(10)
  month = "October"
  smonth = "Oct."
CASE(11)
  month = "November"
  smonth = "Nov."
CASE(12)
  month = "December"
  smonth = "Dec."
CASE DEFAULT
  month = "Unknown"
  smonth = "N/A"
END SELECT
!
END SUBROUTINE INT2MONTH
!
!
!********************************************************
! INT2DAY
! Converts an integer to strings containing the
! day in long and short forms.
! Integer should be between 1 and 12
! (otherwise "Sunday" is returned).
!********************************************************
SUBROUTINE INT2DAY(number,day,sday)
CHARACTER(LEN=16),INTENT(OUT):: day,sday
INTEGER,INTENT(IN):: number
!
SELECT CASE(number)
CASE(0,7)
  day = "Sunday"
  sday = "Sun."
CASE(1)
  day = "Monday"
  sday = "Mon."
CASE(2)
  day = "Tuesday"
  sday = "Tue."
CASE(3)
  day = "Wednesday"
  sday = "Wed."
CASE(4)
  day = "Thursday"
  sday = "Thu."
CASE(5)
  day = "Friday"
  sday = "Fri."
CASE(6)
  day = "Saturday"
  sday = "Sat."
CASE DEFAULT
  day = "Unknown"
  sday = "N/A"
END SELECT
!
END SUBROUTINE INT2DAY
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
!********************************************************
! CHECKNAN
! This subroutine reads an array and, if any number
! is actually NaN (Not a Number) the index is returned,
! otherwise 0 is returned.
!********************************************************
SUBROUTINE CHECKNAN(A,NaNindex)
!
IMPLICIT NONE
CHARACTER(LEN=32):: temp
INTEGER:: i, j
INTEGER:: NaNindex
REAL(dp),DIMENSION(:,:),INTENT(IN):: A
!
NaNindex = 0
DO i=1,SIZE(A,1)
  DO j=1,SIZE(A,2)
    WRITE(temp,'(f12.8)') A(i,j)
    IF( TRIM(ADJUSTL(temp))=="NaN" .OR. SCAN(temp,"Inf").NE.0 ) NaNindex = i
    !As soon as a NaN is found, no need to parse the rest of the array
    IF( NaNindex.NE.0 ) RETURN
  ENDDO
ENDDO
!
END SUBROUTINE CHECKNAN
!
!
!********************************************************
! DO_STATS
! This subroutine performs some simple statistics
! on the numbers contained in a 1-dimensional array.
!********************************************************
SUBROUTINE DO_STATS(array,mi,M,A,D,S)
!
IMPLICIT NONE
INTEGER:: a_size, i
REAL(dp),INTENT(OUT):: mi, M, A, D, S
REAL(dp),DIMENSION(:),INTENT(IN):: array
!
a_size = SIZE(array)
A = 0.d0
D = 0.d0
S = 0.d0
!
!Find the minimum (mi) and maximum (M)
!values in the array
mi = MINVAL(array)
M  = MAXVAL(array)
!
!Calculate the average (A) of all values in the array
DO i=1,a_size
  A = A+array(i)
ENDDO
A = A/DBLE(a_size)
!
!Calculate average absolute deviation (D)
!and standard deviation (S)
DO i=1,a_size
  D = D + DABS(array(i)-A)
  S = S + (array(i)-A)**2
ENDDO
D = D/DBLE(a_size)
S = DSQRT( S/DBLE(a_size) )
!
END SUBROUTINE DO_STATS
!
!
!********************************************************
! CHECK_ORTHOVEC
! This subroutine checks if two vectors are
! orthogonal and returns .TRUE. if they are,
! .FALSE. otherwise.
!********************************************************
SUBROUTINE CHECK_ORTHOVEC(V1,V2,areortho)
!
IMPLICIT NONE
LOGICAL,INTENT(OUT):: areortho
REAL(dp),DIMENSION(3),INTENT(IN):: V1, V2
!
areortho = .FALSE.
IF(DOT_PRODUCT(V1,V2)==0.d0) THEN
  areortho = .TRUE.
ENDIF
!
END SUBROUTINE CHECK_ORTHOVEC
!
!
!********************************************************
! CART2FRAC
! Converts the 3 first columns of an array A, assumed
! to contain cartesian coordinates, into fractional
! coordinates.
!********************************************************
SUBROUTINE CART2FRAC(A,H)
!
IMPLICIT NONE
INTEGER:: i, j
REAL(dp),DIMENSION(3):: V
REAL(dp),DIMENSION(3,3),INTENT(IN):: H
REAL(dp),DIMENSION(3,3):: G
REAL(dp),DIMENSION(:,:):: A
!
IF( SIZE(A,1)>0 .AND. SIZE(A,2)>=3 ) THEN
  CALL INVMAT(H,G)
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,V)
  DO i=1,SIZE(A,1)
    V(:) = A(i,1:3)
    DO j=1,3
      A(i,j) = V(1)*G(1,j) + V(2)*G(2,j) + V(3)*G(3,j)
    ENDDO
  END DO
  !$OMP END PARALLEL DO
ELSE
  WRITE(*,*) 'X!X ERROR: could not transform to fractional,'
  WRITE(*,*) '          inconsistent array size.'
  nerr = nerr+1
ENDIF
!
END SUBROUTINE CART2FRAC
!
!
!********************************************************
! FRAC2CART
! Converts the 3 first columns of an array A, assumed
! to contain fractional coordinates, into cartesian
! coordinates.
!********************************************************
SUBROUTINE FRAC2CART(A,H)
!
IMPLICIT NONE
INTEGER:: i, j
REAL(dp),DIMENSION(3):: V
REAL(dp),DIMENSION(3,3),INTENT(IN):: H
REAL(dp),DIMENSION(:,:):: A
!
IF( SIZE(A,1)>0 .AND. SIZE(A,2)>=3 ) THEN
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,V)
  DO i=1,SIZE(A,1)
    V(:) = A(i,1:3)
    DO j=1,3
      A(i,j) = V(1)*H(1,j) + V(2)*H(2,j) + V(3)*H(3,j)
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
ELSE
  WRITE(*,*) 'X!X ERROR: could not transform to cartesian,'
  WRITE(*,*) '          inconsistent array size.'
  nerr = nerr+1
ENDIF
!
END SUBROUTINE FRAC2CART
!
!
!********************************************************
! COUNT_OUTBOX
! This subroutine counts how many points of array A(x,y,z)
! are outside of the box.
!********************************************************
SUBROUTINE COUNT_OUTBOX(H,A,Nout)
!
IMPLICIT NONE
INTEGER:: i, Nout
REAL(dp),DIMENSION(3,3),INTENT(IN):: H
REAL(dp),DIMENSION(:,:),INTENT(IN):: A
REAL(dp),DIMENSION(SIZE(A,1),SIZE(A,2)):: Afrac
!
Nout = 0
!
IF( SIZE(A,1)>0 ) THEN
  Afrac(:,:) = A(:,:)
  !
  CALL CART2FRAC(Afrac,H)
  !
  DO i=1,SIZE(Afrac,1)
    IF( Afrac(i,1)<0.d0 .OR. Afrac(i,1)>1.d0 .OR. &
      & Afrac(i,2)<0.d0 .OR. Afrac(i,2)>1.d0 .OR. &
      & Afrac(i,3)<0.d0 .OR. Afrac(i,3)>1.d0     ) THEN
      Nout = Nout+1
    ENDIF
  ENDDO
!
ENDIF
!
END SUBROUTINE COUNT_OUTBOX
!
!
!********************************************************
! FIND_NSP
! This subroutine finds the number of different
! entries in a real column array A. The result is
! output in a real array of the form aentries(entry,Nentry)
! where "entry" is the value of the entries and Nentry
! is the number of times this entry appears in A.
!********************************************************
SUBROUTINE FIND_NSP(A,aentries)
!
IMPLICIT NONE
INTEGER:: i, j
INTEGER:: Ndiff  !Number of different entries
LOGICAL:: foundA
REAL(dp),DIMENSION(:),INTENT(IN):: A   !input array that must be analyzed
REAL(dp),DIMENSION(100,2):: atemp      !temporary array containing results
                                       !This assumes that there are no more than 100 different elements in A
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: aentries !array containing result
!
Ndiff = 0
IF(ALLOCATED(aentries)) DEALLOCATE(aentries)
atemp(:,:) = 0.d0
!
!Loop on A and save values to aentries(:,:)
Ndiff = 0
DO i=1,SIZE(A)
  foundA = .FALSE.
  !
  IF( i>1 ) THEN
    !Check if we can find A in the array aentries(:,1)
    DO j=1,SIZE(atemp,1)
      IF( DABS(A(i)-atemp(j,1))<1.d-12 ) THEN
        foundA = .TRUE.
        atemp(j,2) = atemp(j,2)+1
      ENDIF
    ENDDO
  ENDIF
  !
  !If it was not found then it is a new entry
  IF(.NOT.foundA) THEN
    Ndiff = Ndiff+1
    atemp(Ndiff,1) = A(i)
    atemp(Ndiff,2) = 1
  ENDIF
ENDDO
!
ALLOCATE(aentries(Ndiff,2))
aentries(:,:) = 0.d0
!
DO i=1,Ndiff
  aentries(i,:) = atemp(i,:)
ENDDO
!
END SUBROUTINE FIND_NSP
!
!
!********************************************************
! FIND_IF_REDUCED
! This subroutine reads a real array of dimension NxM,
! reads the columns (maximum up to the third column) and
! tries to determine if values are in reduced coordinates.
! "Reduced coordinates" usually mean that all coordinates
! are between 0 and 1, however it is not always so
! obvious. Sometimes atoms can escape the box and have
! a coordinate smaller than 0 or greater than 1; and they
! can obviously not be wrapped back into the box if we
! are unsure if coordinates are reduced or not. Another
! possibility is that all atoms have been shifted by
! an arbitrary vector; their reduced coordinate can
! then be between X and X+1. A way to recognize such
! extreme cases is to work with the standard deviation
! of coordinates along X, Y and Z. If the standard
! deviation is smaller than 1 in all directions of space
! then coordinates have to be reduced.
! One possibility cannot be solved by this routine: if
! atoms diffuse a lot and are not wrapped, even their
! reduced coordinate can take large values (>1) that
! can be confused with values in angströms.
!********************************************************
SUBROUTINE FIND_IF_REDUCED(H,array,isreduced)
!
IMPLICIT NONE
LOGICAL,INTENT(OUT):: isreduced  !does A contain reduced coordinates?
INTEGER:: i, j
REAL(dp):: avg, minmax, D
REAL(dp):: th
REAL(dp),DIMENSION(3,3),INTENT(IN):: H !box vectors
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: array !atom positions
!
isreduced = .TRUE.
avg = 0.d0
minmax = 0.d0
D = 0.d0
th = 3.d0
!
IF( ANY( DABS(array(:,1:3)) > 0.99d0 ) .OR. ANY( DABS(array(:,1:3)) < 1.d-6 ) ) THEN
  !Some coordinates are not contained between 0 and 1
  !Loop on X, Y, Z
  DO i=1,3
    !
    IF( SIZE(array,1) < 11 ) THEN
      th = 0.4d0*VECLENGTH(H(:,i))
    ENDIF
    !
    !Compute difference between max and min coordinates
    minmax = MAXVAL(array(:,i)) - MINVAL(array(:,i))
    !
    !Compute average of all values
    DO j=1,SIZE(array,1)
      avg = avg + DABS(array(j,i))
    ENDDO
    avg = avg/DBLE(SIZE(array,1))
    !
    !Calculate average absolute deviation (D)
    DO j=1,SIZE(array,1)
      D = D + DABS(array(j,i)-avg)
    ENDDO
    D = D/DBLE(SIZE(array,1))
    !
    !If minmax is much greater than 1, or if D is greater than 1,
    !then coordinates are Cartesian
    IF( avg > 1.5d0 .OR. minmax > th .OR. D > th ) THEN
      isreduced = .FALSE.
    ENDIF
    !
  ENDDO
  !
ENDIF
!
END SUBROUTINE FIND_IF_REDUCED
!
!
!********************************************************
!  UNWRAP
!  This subroutine tries to unwrap a row of atoms.
!  Sometimes atoms cross the boundaries of the supercell,
!  and some programs wraps them back into it.
!  This subroutine attempts to "unwrap" such atoms.
!  It only works with a row of atoms perpendicular to the
!  direction "dir", and PBC are applied only along dir.
!********************************************************
!
SUBROUTINE UNWRAP(atoms_v,dir,vector,threshold,nwrapat)
!
IMPLICIT NONE
INTEGER:: i, j, dir, nwrapat
REAL(dp),INTENT(IN):: vector, threshold
REAL(dp),DIMENSION(:,:),INTENT(INOUT):: atoms_v
!
nwrapat=0
j=0
DO i=2, SIZE(atoms_v(:,1))
  DO WHILE( DABS(atoms_v(i,dir)-atoms_v(i-1,dir)) > threshold )
    DO WHILE( atoms_v(i,dir)-atoms_v(i-1,dir)> threshold )
      atoms_v(i,dir) = atoms_v(i,dir)-vector
    ENDDO
    DO WHILE( atoms_v(i,dir)-atoms_v(i-1,dir)<-threshold )
      atoms_v(i,dir) = atoms_v(i,dir)+vector
    ENDDO
    IF(i.NE.j) THEN
      j=i
      nwrapat=nwrapat+1
    ENDIF
  ENDDO
ENDDO
!
END SUBROUTINE UNWRAP
!
!
!********************************************************
!  LIST_SELECTED_ATOMS
!  This subroutine generates a table containing
!  chemical species and number of selected atoms.
!********************************************************
!
SUBROUTINE LIST_SELECTED_ATOMS(P,SELECT,selectedlist)
!
IMPLICIT NONE
INTEGER:: i, j
INTEGER:: Nspecies
INTEGER,DIMENSION(20,2):: templist !assume maximum 20 different chemical species
INTEGER,DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: selectedlist
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT
REAL(dp),DIMENSION(:,:),INTENT(IN):: P
!
Nspecies=0
templist(:,:) = 0
IF(ALLOCATED(selectedlist)) DEALLOCATE(selectedlist)
!
DO i=1,SIZE(P,1)
  IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
    DO j=1,SIZE(templist,1)
      IF( templist(j,1) == NINT(P(i,4)) ) THEN
        templist(j,2) = templist(j,2) + 1
        EXIT
      ELSEIF( templist(j,1) == 0 ) THEN
        templist(j,1) = NINT(P(i,4))
        templist(j,2) = templist(j,2) + 1
        Nspecies = Nspecies + 1
        EXIT
      ENDIF
    ENDDO
  ENDIF
ENDDO
!
ALLOCATE(selectedlist(Nspecies,2))
DO j=1,SIZE(selectedlist,1)
  selectedlist(j,:) = templist(j,:)
ENDDO
!
END SUBROUTINE LIST_SELECTED_ATOMS
!
!
!
END MODULE subroutines
