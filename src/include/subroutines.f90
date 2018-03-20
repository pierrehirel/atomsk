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
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 06 March 2018                                    *
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
!* RESIZE_ARRAY2       resize a 2-D array, keeping its content                    *
!* CHECK_ARRAY_CONSISTENCY checks that arrays P, S, AUX, AUXNAMES are consistent  *
!* STR2BOOL            transforms a string into a boolean value                   *
!* INT2MONTH           transforms an integer into a month                         *
!* INT2DAY             transforms an integer into a day                           *
!* BOX2DBLE            transforms a string into a real number                     *
!* CHECKNAN            checks an array for NaN (Not a Number) values              *
!* GEN_NRANDNUMBERS    generate a list of N random real numbers in [0.0,1.0]      *
!* GEN_NRANDGAUSS      generate a list of N numbers with gaussian distribution    *
!* DO_STATS            do some statistics on a 1-dimensional array                *
!* CHECK_ORTHOVEC      checks if two vectors are normal to each other             *
!* CART2FRAC           converts cartesian coordinates to fractional               *
!* FRAC2CART           converts fractional coordinates to cartesian               *
!* INDEX_MILLER        find the indices of a plane from a string                  *
!* INDEX_MILLER_HCP    find the indices of a plane from a string (hcp lattices)   *
!* ELAST2TENSOR        converts from Voigt notation to full elastic tensor        *
!* FIND_NSP            find the number of different entries in an array           *
!* FIND_IF_REDUCED     determines if values of an array are within 0 and 1        *
!* UNWRAP              unwrap atoms that have jumped a boundary                   *
!* DERIVATIVE          calculate the derivative of a function                     *
!* CHECK_CTENSOR       checks if an elastic tensor is symmetric                   *
!* COMPFORMULA         extracts a compound formula from atom site lists P and AUX *
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
! RESIZE_DBLEARRAY2
! This routine changes the size of the provided Array.
! Initial data is preserved in the new array.
! If the new size is larger, unknown data is set to zero.
! If the new size is smaller, some data is lost.
!********************************************************
SUBROUTINE RESIZE_DBLEARRAY2(Array,L1,L2,status)
!
IMPLICIT NONE
INTEGER:: i, j
INTEGER,OPTIONAL:: status  !Success=0; failure=1
INTEGER,INTENT(IN):: L1, L2  !new sizes of Array
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: Array !the array to resize
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: temp_array !temporary copy of Array
!
IF(PRESENT(status)) status = 0
!
IF( .NOT.ALLOCATED(Array) ) THEN
  !Allocate Array with required size and fill it with zeros
  ALLOCATE(Array(L1,L2))
  Array(:,:) = 0.d0
  !
ELSE
  !Array is already allocated => resize it
  IF( L1>0 .AND. L2>0 ) THEN
    !
    ALLOCATE( temp_array(L1,L2) )
    temp_array(:,:) = 0.d0
    DO i=1,MIN(L1,SIZE(Array,1))
      DO j=1,MIN(L2,SIZE(Array,2))
        temp_array(i,j) = Array(i,j)
      ENDDO
    ENDDO
    !
    DEALLOCATE(Array)
    ALLOCATE( Array(L1,L2) )
    Array(:,:) = temp_array(:,:)
    !
    DEALLOCATE(temp_array)
    !
  ELSE
    !i.e. if L1<=0 or L2<=0 => problem
    IF(PRESENT(status)) status = 1
  ENDIF
  !
ENDIF
!
END SUBROUTINE RESIZE_DBLEARRAY2
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
      number = MIN(0.d0,MINVAL(vbox)) + areal*( SUM(DABS(vbox)) )
    ELSEIF( op==':' .OR. op=='/' ) THEN
      IF( DABS(vbox(1))<1.d-12 .OR. DABS(vbox(2))<1.d-12 .OR. DABS(vbox(3))<1.d-12 ) THEN
        number = HUGE(areal)
      ELSE
        number = MIN(0.d0,MINVAL(vbox)) + areal/( SUM(DABS(vbox)) )
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
      number = MIN(0.d0,MINVAL(vbox)) + areal*( SUM(DABS(vbox)) )
    ELSEIF( op==':' .OR. op=='/' ) THEN
      IF( DABS(areal)<1.d-12 ) THEN
        number = HUGE(vbox)
      ELSE
        number = MIN(0.d0,MINVAL(vbox)) + ( SUM(DABS(vbox)) )/areal
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
! GEN_NRANDNUMBERS
! This subroutine generates N random real numbers
! between 0 and 1 with a uniform distribution,
! and returns them in the array randarray(:).
! It uses the system clock as a seed,
! which should ensure that the seed is different
! every time this subroutine is called.
! Note: using a KIND=8 integer for system clock returns
!      the clock time in µs or ns, instead of ms for
!      the default KIND=4. This is to make sure that we
!      obtain a different seed if this subroutine
!      is called twice in a row in a very short
!      period of time.
!********************************************************
SUBROUTINE GEN_NRANDNUMBERS(N,randarray)
!
IMPLICIT NONE
INTEGER:: i, k
INTEGER(KIND=8):: clock !system clock
INTEGER,INTENT(IN):: N  !number of random numbers to generate
INTEGER,DIMENSION(:),ALLOCATABLE:: seed !seed for generating random numbers
REAL(dp),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: randarray    !random numbers
!
!Generate N random numbers
CALL RANDOM_SEED(SIZE=k)
IF(ALLOCATED(seed)) DEALLOCATE(seed)
ALLOCATE(seed(k))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 42*(/ (i-1, i=1,k) /)
CALL RANDOM_SEED(PUT=seed)
DEALLOCATE(seed)
IF( ALLOCATED(randarray) ) DEALLOCATE(randarray)
ALLOCATE(randarray(N))
randarray(:) = 0.d0
CALL RANDOM_NUMBER(randarray)
!
END SUBROUTINE GEN_NRANDNUMBERS
!
!
!********************************************************
! GEN_NRANDGAUSS
! This subroutine generates N random real numbers
! between 0 and 1 with a gaussian distribution,
! with a mean of 0 and a variance of 1.
! and returns them in the array randarray(:).
! It uses the system clock as a seed,
! which should ensure that the seed is different
! every time this subroutine is called.
! Note: using a KIND=8 integer for system clock returns
!      the clock time in µs or ns, instead of ms for
!      the default KIND=4. This is to make sure that we
!      obtain a different seed if this subroutine
!      is called twice in a row in a very short
!      period of time.
!********************************************************
SUBROUTINE GEN_NRANDGAUSS(N,randarray)
!
IMPLICIT NONE
INTEGER:: i, k
INTEGER(KIND=8):: clock !system clock
INTEGER,INTENT(IN):: N  !number of random numbers to generate
INTEGER,DIMENSION(:),ALLOCATABLE:: seed !seed for generating random numbers
REAL(dp),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: randarray    !final random numbers
REAL(dp),DIMENSION(:),ALLOCATABLE:: randarraytemp            !random numbers
!
!Generate 2N random numbers with uniform distribution
CALL RANDOM_SEED(SIZE=k)
IF(ALLOCATED(seed)) DEALLOCATE(seed)
ALLOCATE(seed(k))
CALL SYSTEM_CLOCK(COUNT=clock)
seed = clock + 42*(/ (i-1, i=1,k) /)
CALL RANDOM_SEED(PUT=seed)
DEALLOCATE(seed)
IF( ALLOCATED(randarraytemp) ) DEALLOCATE(randarraytemp)
ALLOCATE(randarraytemp(2*N))
randarraytemp(:) = 0.d0
CALL RANDOM_NUMBER(randarraytemp)
!
!Generate the final gaussian distribution
IF( ALLOCATED(randarray) ) DEALLOCATE(randarray)
ALLOCATE(randarray(N))
randarray(:) = 0.d0
DO i=1,N
  randarray(i) = DSQRT( -2.0d0*DLOG(randarraytemp(i)) ) * DCOS( 2.0d0*pi*randarraytemp(N+i) )
ENDDO
!
END SUBROUTINE GEN_NRANDGAUSS
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
INTEGER:: i
REAL(dp):: P1, P2, P3
REAL(dp),DIMENSION(3,3),INTENT(IN):: H
REAL(dp),DIMENSION(3,3):: G
REAL(dp),DIMENSION(:,:):: A
!
IF( SIZE(A(:,1)).NE.0 .AND. SIZE(A(1,:))>=3 ) THEN
    CALL INVMAT(H,G)
    DO i=1,SIZE(A(:,1))
      P1 = A(i,1)
      P2 = A(i,2)
      P3 = A(i,3)
      A(i,1) = P1*G(1,1) + P2*G(2,1) + P3*G(3,1)
      A(i,2) = P1*G(1,2) + P2*G(2,2) + P3*G(3,2)
      A(i,3) = P1*G(1,3) + P2*G(2,3) + P3*G(3,3)
    END DO
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
INTEGER:: i
REAL(dp):: P1, P2, P3
REAL(dp),DIMENSION(3,3),INTENT(IN):: H
REAL(dp),DIMENSION(:,:):: A
!
IF( SIZE(A(:,1)).NE.0 .AND. SIZE(A(1,:))>=3 ) THEN
    DO i=1,SIZE(A(:,1))
      P1 = A(i,1)
      P2 = A(i,2)
      P3 = A(i,3)
      A(i,1) = P1*H(1,1) + P2*H(2,1) + P3*H(3,1)
      A(i,2) = P1*H(1,2) + P2*H(2,2) + P3*H(3,2)
      A(i,3) = P1*H(1,3) + P2*H(2,3) + P3*H(3,3)
    ENDDO
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
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: A
REAL(dp),DIMENSION(SIZE(A,1),SIZE(A,2)):: Afrac
!
Nout = 0
!
IF( ALLOCATED(A) .AND. SIZE(A,1)>0 ) THEN
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
!  INDEX_MILLER
!  This subroutine finds the Miller indices [hkl]
!  from a string. The string can contain signed integers
!  attached together, or separated by underscores,
!  with or without brackets, e.g.:
!         1-10       [1-10]
!        1_-1_0     [1_-1_0]
!  The notation without underscore only allows to use
!  one-digit signed integers (like "1-10"), while the
!  notation allows to use greater integers,
!  like "[12_-15_14]".
!  The "planestring" is converted to a real array
!  containing the Miller indices, e.g. "1-10" will be
!  converted into (1.d0  -1.d0  0.d0).
!********************************************************
!
SUBROUTINE INDEX_MILLER(planestring,planeindices,ifail)
!
IMPLICIT NONE
CHARACTER(LEN=16),INTENT(IN):: planestring
CHARACTER(LEN=16):: temp, temp2
INTEGER:: i, m, mint
INTEGER:: ifail
INTEGER:: strpos
REAL(dp):: msign   !sign of value
REAL(dp),DIMENSION(3):: planeindices
!
ifail=0
planeindices(:) = 0.d0
!
m = 1
mint = 0
msign = 1.d0
temp = planestring
!
!If there are brackets, remove them
IF( temp(1:1)=='[' ) THEN
  temp = ADJUSTL(temp(2:))
ENDIF
strpos = LEN_TRIM(temp)
IF( temp(strpos:strpos)==']' ) THEN
  temp = temp(:strpos-1)
ENDIF
!
IF( SCAN(temp,'_').NE.0 ) THEN
  !values are separated by underscores, e.g. "1_-1_0"
  !read first value
  strpos = SCAN(temp,'_')
  READ(temp(1:strpos-1),*,ERR=100,END=100) planeindices(1)
  temp = temp(strpos+1:)
  !read second value
  strpos = SCAN(temp,'_')
  READ(temp(1:strpos-1),*,ERR=100,END=100) planeindices(2)
  temp = temp(strpos+1:)
  !read third value
  READ(temp,*,ERR=100,END=100) planeindices(3)
ELSE
  !values are given as attached integers, e.g. "1-10"
  !convert them to a vector like [1 -1 0]
  DO i=1, LEN_TRIM(temp)
    READ(temp(i:i),*,ERR=100,END=100) temp2
    IF(temp2=='-') THEN
      msign = -1.d0
    ELSEIF(temp2=='+') THEN
      msign = 1.d0
    ELSEIF(LEN_TRIM(temp2)>0) THEN
      READ(temp2,*,ERR=100,END=100) mint
      planeindices(m) = msign*DBLE(mint)
      m = m+1
      msign = 1.d0
    ENDIF
  ENDDO
ENDIF
RETURN
!
100 CONTINUE
ifail=1
!
END SUBROUTINE INDEX_MILLER
!
!
!********************************************************
!  INDEX_MILLER_HCP
!  This subroutine finds the Miller indices [hkil]
!  from a string. The string can contain signed integers
!  attached together, or separated by underscores,
!  with or without brackets, e.g.:
!         1-100       [1-100]
!        1_-1_00     [1_-1_0_0]
!  The notation without underscore only allows to use
!  one-digit signed integers (like "1-10"), while the
!  notation allows to use greater integers,
!  like "[12_-15_3_14]".
!  The "planestring" is converted to a real array
!  containing only the 3 relevant Miller indices hkl,
! e.g. "1-10" will be converted into (1.d0  -1.d0  0.d0).
!********************************************************
!
SUBROUTINE INDEX_MILLER_HCP(planestring,planeindices,ifail)
!
IMPLICIT NONE
CHARACTER(LEN=16),INTENT(IN):: planestring
CHARACTER(LEN=16):: temp, temp2
INTEGER:: i, m, mint
INTEGER:: ifail !0=success; 1=error while reading string; 2=h+k not equal to -i
INTEGER:: strpos
REAL(dp):: msign   !sign of value
REAL(dp),DIMENSION(3):: planeindices
!
ifail=0
planeindices(:) = 0.d0
!
m = 1
mint = 0
msign = 1.d0
temp = planestring
!
!If there are brackets, remove them
IF( temp(1:1)=='[' ) THEN
  temp = ADJUSTL(temp(2:))
ENDIF
strpos = LEN_TRIM(temp)
IF( temp(strpos:strpos)==']' ) THEN
  temp = temp(:strpos-1)
ENDIF
!
IF( SCAN(temp,'_').NE.0 ) THEN
  !values are separated by underscores, e.g. "1_-1_0_0"
  !read first value
  strpos = SCAN(temp,'_')
  READ(temp(1:strpos-1),*,ERR=100,END=100) planeindices(1)
  temp = temp(strpos+1:)
  !read second value
  strpos = SCAN(temp,'_')
  READ(temp(1:strpos-1),*,ERR=100,END=100) planeindices(2)
  temp = temp(strpos+1:)
  !read third value
  strpos = SCAN(temp,'_')
  READ(temp(1:strpos-1),*,ERR=100,END=100) planeindices(3)
  !Check that -i=h+k
  IF( -1*NINT(planeindices(3))==NINT(planeindices(1))+NINT(planeindices(2)) ) THEN
    !Notation is good, go on
    temp = temp(strpos+1:)
    !read fourth value
    READ(temp,*,ERR=100,END=100) planeindices(3)
  ELSE
    !h+k is not equal to -i => return with error
    ifail = 2
    RETURN
  ENDIF
ELSE
  !values are given as attached integers, e.g. "1-100"
  !convert them to a vector like [1 -1 0]
  strpos=0
  DO i=1, LEN_TRIM(temp)
    READ(temp(i:i),*,ERR=100,END=100) temp2
    IF(temp2=='-') THEN
      msign = -1.d0
    ELSEIF(temp2=='+') THEN
      msign = 1.d0
    ELSEIF(LEN_TRIM(temp2)>0) THEN
      READ(temp2,*,ERR=100,END=100) mint
      planeindices(m) = msign*DBLE(mint)
      strpos=strpos+1
      IF( strpos==3 ) THEN
        !Check that -i=h+k
        IF( -1*NINT(planeindices(3)).NE.NINT(planeindices(1))+NINT(planeindices(2)) ) THEN
          ifail = 2
          RETURN
        ENDIF
      ENDIF
      m = MIN(3,m+1)
      msign = 1.d0
    ENDIF
  ENDDO
ENDIF
RETURN
!
100 CONTINUE
ifail=1
!
END SUBROUTINE INDEX_MILLER_HCP
!
!
!********************************************************
!  ELAST2TENSOR
!  This subroutine sets up the 9x9 elastic tensor,
!  provided the 9 (Voigt notation) elastic constants.
!********************************************************
!
SUBROUTINE ELAST2TENSOR(elcst,eltens)
!
IMPLICIT NONE
REAL(dp),DIMENSION(9),INTENT(IN):: elcst !C11,C22,C33,C23,C31,C12,C44,C55,C66
REAL(dp),DIMENSION(:,:),INTENT(INOUT):: eltens
!
eltens(:,:) = 0.d0
!
eltens(1,1) = elcst(1)
eltens(1,2) = elcst(6)
eltens(1,3) = elcst(5)
eltens(2,1) = elcst(6)
eltens(2,2) = elcst(2)
eltens(2,3) = elcst(4)
eltens(3,1) = elcst(5)
eltens(3,2) = elcst(4)
eltens(3,3) = elcst(3)
eltens(4,4) = elcst(7)
eltens(5,5) = elcst(8)
 eltens(6,6) = elcst(9)
IF( SIZE(eltens,1)==9 .AND. SIZE(eltens,2)==9 ) THEN
  !Fill the rest by symmetry
  eltens(1:3,7:9) = eltens(1:3,4:6)
  eltens(4:6,7:9) = eltens(4:6,4:6)
  eltens(7:9,1:3) = eltens(4:6,1:3)
  eltens(7:9,4:6) = eltens(4:6,4:6)
  eltens(7:9,7:9) = eltens(4:6,4:6)
ENDIF

!
END SUBROUTINE ELAST2TENSOR
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
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries !array containing result
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
! reduced coordinate can take large values.
!********************************************************
SUBROUTINE FIND_IF_REDUCED(array,isreduced)
!
IMPLICIT NONE
LOGICAL,INTENT(OUT):: isreduced  !does A contain reduced coordinates?
INTEGER:: i
REAL(dp):: mi, M, A, D, S  !statistics
REAL(dp),DIMENSION(3):: cm
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: array
!
IF( SIZE(array,1) > 4 ) THEN
  cm(:) = 0.d0
  !We read the columns of array and compute the normal deviation of all values
  !- if array has one, two or three columns we read them
  !- if array has four or more columns we read only the 3 first columns
  DO i=1, MIN( SIZE(array,2),3 )
    !Compute the statistics on coordinates in this column
    CALL DO_STATS(array(:,i),mi,M,A,D,S)
    !Save deviation on atom positions
    cm(i) = S
  ENDDO
  !
  !If deviation is smaller than 1 in all
  !directions, then it is reduced coordinates
  IF( cm(1)<=0.9d0 .AND. cm(2)<=0.9d0 .AND. cm(3)<=0.9d0 ) THEN
    isreduced = .TRUE.
  ELSE
    isreduced = .FALSE.
  ENDIF
  !
ELSE
  !Small number of atoms => the method above may produce wrong results
  !if an atom is placed at the origin. E.g. in a unit cell containing
  !only two atoms and one is at (0,0,0) and the other at (1.3,1.3,1.3),
  !the standard deviation will be smaller than 1, thus wrongly leading
  !to reduced coordinates. Therefore when the number of atoms is small
  !another method is used.
  isreduced = .TRUE.
  DO i=2,SIZE(array,1)
    IF( VECLENGTH( array(i,1:3)-array(1,1:3) ) > 1.25d0 ) THEN
      isreduced = .FALSE.
    ENDIF
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
!  CHECK_CTENSOR
!  This subroutine checks if an elastic tensor is
!  symmetric.
!********************************************************
!
SUBROUTINE CHECK_CTENSOR(C_tensor,status)
!
IMPLICIT NONE
INTEGER:: i, j
INTEGER:: status  !=0 if tensor is OK
                  !=1 if tensor is not symmetric, i.e. C(i,j).NE.C(j,i)
REAL(dp),DIMENSION(9,9),INTENT(IN):: C_tensor
!
status = 0
DO i=2,9
  DO j=1,i-1
    IF( DABS(C_tensor(i,j)-C_tensor(j,i))>1.d-3 ) THEN
      status=1
    ENDIF
  ENDDO
ENDDO
!
END SUBROUTINE CHECK_CTENSOR
!
!
!********************************************************
!  COMPFORMULA
!  This subroutine extracts the compound formula from the
!  given atomic site list P taking the auxiliary list AUX
!  into account if it contains site occupancies.
!  The formula is returned as string.
!  The compound summed atomic mass is returned as well.
!  (C) Juri Barthel 
!*     Gemeinschaftslabor fuer Elektronenmikroskopie
!*     RWTH Aachen (GERMANY)
!*     ju.barthel@fz-juelich.de
!********************************************************
!
SUBROUTINE COMPFORMULA(P,AUXNAMES,AUX,formula,mass)
!
USE atoms
!
IMPLICIT NONE
INTEGER:: i, j, iaux
INTEGER:: occ
INTEGER:: NP, Naux, NS, NA, NF
INTEGER,DIMENSION(ATOMMAXZ+10):: ispecies
REAL(dp):: PO, amass, socc
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
REAL(dp),INTENT(OUT):: mass
CHARACTER(LEN=2):: species
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128), INTENT(OUT):: formula
CHARACTER(LEN=128) :: temp, msg
!
! 
!Initialize variables
mass = 0.0d+0
amass = 0.0d+0
formula = ''
temp = ''
msg = ''
Naux = 0
occ = 0 ! no occupancy information by default
PO = 1.0d+0 ! 1.0 occupancy of each atomic site by default
ispecies = 0 ! preset species hash
NP = SIZE(P,1) ! number of sites
IF (ALLOCATED(AUX).AND.ALLOCATED(AUXNAMES)) THEN ! Get number of auxiliary properties
  Naux = SIZE(AUX,2) ! ? Redundant but maybe better MIN(SIZE(AUX,2),SIZE(AUXNAMES))
ENDIF
IF (Naux>0) THEN
  DO iaux=1, Naux
    IF("occ"==TRIM(AUXNAMES(iaux))) occ = iaux
  ENDDO
END IF
!
!Get different atomic species
!    from list P(:,4) = atomic numbers of all atomic sites
CALL FIND_NSP(P(:,4),aentries)
NS = SIZE(aentries,1) ! number of different species
!Get occupation numbers of the different species
IF (NS>0) THEN
  ispecies = 0
  !Remember the index of each species in aentries
  !     but only for those species which are on atom sites
  DO i=1, NS
    j = NINT(aentries(i,1))
    IF (j>0.AND.aentries(i,2)>0.0d+0) ispecies(j) = i
  ENDDO
  !Erase atom species counts in aentries
  aentries(1:NS,2) = 0.0d+0
  !Sum up the occupancies
  DO i=1, NP
    IF (P(i,4)<1.or.P(i,4)>ATOMMAXZ) CYCLE ! unknown element
    j = ispecies(NINT(P(i,4))) ! species index in aentries
    IF (j<=0) CYCLE ! invalid index
    IF (occ>0) PO = AUX(i,occ) ! set individual occupancy from AUX
    aentries(j,2) = aentries(j,2) + PO ! accumulation
  ENDDO
ENDIF
! Create formula string and sum up atomic masses
DO i=1, NS
  CALL ATOMSPECIES(aentries(i,1),species)
  temp=''
  socc = aentries(i,2)
  IF (socc<=0.d0) CYCLE ! this species is not in the compound
  NA = INT(socc)
  NF = INT((socc - NA)*1.0d+1)
  IF( NA>0 .or. NF>0 ) THEN
    WRITE(temp,'(i10)') NA
    IF( NF>0) THEN
      WRITE(temp,'(a,i1)') TRIM(ADJUSTL(temp))//".", NF
    END IF
    temp = TRIM(species)//TRIM(ADJUSTL(temp))
  ENDIF
  msg = TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(temp))
  CALL ATOMMASS(species,amass)
  mass = mass + amass*socc
ENDDO
formula = msg
!
IF (ALLOCATED(aentries)) DEALLOCATE(aentries)
!
END SUBROUTINE COMPFORMULA
!
!
!
END MODULE subroutines
