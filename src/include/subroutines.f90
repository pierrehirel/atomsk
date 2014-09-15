MODULE subroutines
!
!**********************************************************************************
!*  SUBROUTINES                                                                   *
!**********************************************************************************
!* This module contains general subroutines dealing                               *
!* with files, manipulation of strings, arrays and so on.                         *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel (unless specified otherwise, cf each subroutine)  *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 12 Sept. 2014                                    *
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
!* STR2BOOL            transforms a string into a boolean value                   *
!* INT2MONTH           transforms an integer into a month                         *
!* INT2DAY             transforms an integer into a day                           *
!* BOX2DBLE            transforms a string into a real number                     *
!* PROGBAR             writes a progress bar to a string                          *
!* CHECKNAN            checks an array for NaN (Not a Number) values              *
!* GEN_NRANDNUMBERS    generate a list of N random real numbers in [0.0,1.0]      *
!* GEN_NRANDGAUSS      generate a list of N numbers with gaussian distribution    *
!* BUBBLESORT          sorts an array by increasing or decreasing values          *
!* PACKSORT            sorts an array by packing identical values together        *
!* DO_STATS            do some statistics on a 1-dimensional array                *
!* CHECK_ORTHOVEC      checks if two vectors are normal to each other             *
!* INVMAT              inverts a NxN matrix                                       *
!* CONVMAT             converts conventional vectors into a matrix                *
!* MATCONV             converts matrix into conventional vectors                  *
!* VOLUME_PARA         computes the volume of a parallelepiped                    *
!* CART2FRAC           converts cartesian coordinates to fractional               *
!* FRAC2CART           converts fractional coordinates to cartesian               *
!* INDEX_MILLER        find the indices of a plane from a string                  *
!* ELAST2TENSOR        converts from Voigt notation to full elastic tensor        *
!* FIND_NSP            find the number of different entries in an array           *
!* FIND_IF_REDUCED     determines if values of an array are within 0 and 1        *
!* UNWRAP              unwrap atoms that have jumped a boundary                   *
!* DERIVATIVE          calculate the derivative of a function                     *
!* CHECK_CTENSOR       checks if an elastic tensor is symmetric                   *
!**********************************************************************************
!
!
USE comv
USE constants
USE functions
!
IMPLICIT NONE
!
!
CONTAINS
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
CASE(0)
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
CASE(7)
  day = "Sunday"
  sday = "Sun."
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
  number = MAXVAL(vbox) + DABS(MINVAL(vbox))
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
      !no operator before of after "box" => abandon ship
      nerr=nerr+1
      RETURN
    ENDIF
  ENDIF
  !
  !Evaluate the expression and save it to "number"
  IF( posop<posbox ) THEN
    !Read what is before 'BOX'
    READ(text(1:posop-1),*,END=800,ERR=800) areal
    !Do the operation
    IF( op=='+' ) THEN
      number = areal+MAXVAL(vbox)
    ELSEIF( op=='-' ) THEN
      number = areal-MAXVAL(vbox)
    ELSEIF( op=='*' ) THEN
      number = MINVAL(vbox) + areal*( SUM(DABS(vbox)) )
    ELSEIF( op==':' .OR. op=='/' ) THEN
      IF( ANY(vbox(:)<1.d-12) ) THEN
        number = HUGE(areal)
      ELSE
        number = MINVAL(vbox) + areal/( SUM(DABS(vbox)) )
      ENDIF
    ENDIF
  ELSE
    !Read what is after the operator
    READ(text(posop+1:),*,END=800,ERR=800) areal
    !Do the operation
    IF( op=='+' ) THEN
      number = MAXVAL(vbox)+areal
    ELSEIF( op=='-' ) THEN
      number = MAXVAL(vbox)-areal
    ELSEIF( op=='*' ) THEN
      number = MINVAL(vbox) + areal*( SUM(DABS(vbox)) )
    ELSEIF( op==':' .OR. op=='/' ) THEN
      IF( DABS(areal)<1.d-12 ) THEN
        number = HUGE(vbox)
      ELSE
        number = MINVAL(vbox) + ( SUM(DABS(vbox)) )/areal
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
!********************************************************
! PROGBAR
! This subroutine writes a progress bar in the string
! "pbar", according to the given percentage. The number
! "percent" must be given in % (i.e. between 0 and 100).
!********************************************************
SUBROUTINE PROGBAR(percent,pbar)
!
IMPLICIT NONE
CHARACTER(LEN=22):: pbar
REAL(dp),INTENT(IN):: percent
!
pbar = ""
IF(percent<5.d0) THEN
  pbar = "[>                   ]"
ELSEIF(percent>=5.d0 .AND. percent<10.d0) THEN
  pbar = "[=>                  ]"
ELSEIF(percent>=10.d0 .AND. percent<15.d0) THEN
  pbar = "[==>                 ]"
ELSEIF(percent>=15.d0 .AND. percent<20.d0) THEN
  pbar = "[===>                ]"
ELSEIF(percent>=20.d0 .AND. percent<25.d0) THEN
  pbar = "[====>               ]"
ELSEIF(percent>=25.d0 .AND. percent<30.d0) THEN
  pbar = "[=====>              ]"
ELSEIF(percent>=30.d0 .AND. percent<35.d0) THEN
  pbar = "[======>             ]"
ELSEIF(percent>=35.d0 .AND. percent<40.d0) THEN
  pbar = "[=======>            ]"
ELSEIF(percent>=40.d0 .AND. percent<45.d0) THEN
  pbar = "[========>           ]"
ELSEIF(percent>=45.d0 .AND. percent<50.d0) THEN
  pbar = "[=========>          ]"
ELSEIF(percent>=50.d0 .AND. percent<55.d0) THEN
  pbar = "[==========>         ]"
ELSEIF(percent>=55.d0 .AND. percent<60.d0) THEN
  pbar = "[===========>        ]"
ELSEIF(percent>=60.d0 .AND. percent<65.d0) THEN
  pbar = "[============>       ]"
ELSEIF(percent>=65.d0 .AND. percent<70.d0) THEN
  pbar = "[=============>      ]"
ELSEIF(percent>=70.d0 .AND. percent<75.d0) THEN
  pbar = "[==============>     ]"
ELSEIF(percent>=75.d0 .AND. percent<80.d0) THEN
  pbar = "[===============>    ]"
ELSEIF(percent>=80.d0 .AND. percent<85.d0) THEN
  pbar = "[================>   ]"
ELSEIF(percent>=85.d0 .AND. percent<90.d0) THEN
  pbar = "[=================>  ]"
ELSEIF(percent>=90.d0 .AND. percent<95.d0) THEN
  pbar = "[==================> ]"
ELSEIF(percent>=95.d0 .AND. percent<100.d0) THEN
  pbar = "[===================>]"
ELSEIF(percent>=100.d0) THEN
  pbar = "[====================]"
ENDIF
!
END SUBROUTINE PROGBAR
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
! BUBBLESORT
! This subroutine sorts a MxN array by increasing or
! decreasing values of its column 'col'.
! It uses the so-called "bubble-sort" algorithm.
! IMPORTANT NOTE: this algorithm swaps entire columns
! of an array, up to N*(N-1) times (N=number of columns),
! so you may find it highly unefficient for very
! large arrays.
!********************************************************
SUBROUTINE BUBBLESORT(A,col,order)
!
IMPLICIT NONE
CHARACTER(LEN=4):: order  !up or down
INTEGER:: i, j
INTEGER:: col
REAL(dp),DIMENSION(:,:),INTENT(INOUT):: A
REAL(dp),DIMENSION(SIZE(A,2)) :: col_value
!
 col_value(:) = 0.d0
!
IF(order=='up') THEN
  DO j=1,SIZE(A,1)
    DO i=j+1,SIZE(A,1)
      !If element i is smaller than element j, we swap them
      IF( A(i,col) < A(j,col) ) THEN
        col_value(:) = A(i,:)
        A(i,:) = A(j,:)
        A(j,:) = col_value(:)
      ENDIF
    ENDDO
  ENDDO
!
ELSE
  DO j=1,SIZE(A,1)
    DO i=j+1,SIZE(A,1)
      !If element i is greater than element j, we swap them
      IF( A(i,col) > A(j,col) ) THEN
        col_value(:) = A(i,:)
        A(i,:) = A(j,:)
        A(j,:) = col_value(:)
      ENDIF
    ENDDO
  ENDDO
ENDIF
!
END SUBROUTINE BUBBLESORT
!
!
!********************************************************
! PACKSORT
! This subroutine sorts a MxN array by packing
! identical values together, but without sorting them
! (unlike the BUBBLESORT above). E.g. if values like:
!      2 4 2 3 1 1 2 4 3 3 2 1 4
! are given, this subroutine will pack identical
! values so the result is:
!      2 2 2 2 4 4 4 3 3 3 1 1 1
!********************************************************
SUBROUTINE PACKSORT(A,col)
!
IMPLICIT NONE
INTEGER:: i, j, k
INTEGER:: col, last
REAL(dp),DIMENSION(:,:),INTENT(INOUT):: A
REAL(dp),DIMENSION(SIZE(A(1,:))):: Atemp
!
IF(col>SIZE(A(1,:))) col = SIZE(A(1,:))
!
DO i=1,SIZE(A(:,1))
  last=i
  DO j=i+1,SIZE(A(:,1))
    IF( A(j,col)==A(i,col) ) THEN
      !Only consider A(j)==A(i)
      IF(j==last+1) THEN
        !If the two values are contiguous just go on
        last=j
      ELSE
        !If the value is elsewhere in the array,
        !pack it with the others
        Atemp(:) = A(j,:)
        k=j
        DO k=j,last+1,-1
          A(k,:) = A(k-1,:)
        ENDDO
        last=last+1
        A(last,:) = Atemp(:)
      ENDIF
    ENDIF
  ENDDO
ENDDO

!
END SUBROUTINE PACKSORT
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
! INVMAT
! This subroutine inverts a NxN matrix M
! and outputs the result into the matrix G.
!********************************************************
SUBROUTINE INVMAT(M,G,status)
!
IMPLICIT NONE
REAL(dp),DIMENSION(:,:),INTENT(IN):: M
REAL(dp),DIMENSION(:,:),INTENT(OUT):: G
INTEGER,INTENT(OUT),OPTIONAL:: status
INTEGER:: i
INTEGER:: LWORK !for LAPACK routine DGETRI
INTEGER,DIMENSION(SIZE(M,1)):: IPIV !for LAPACK routine DGETRI
REAL(dp):: det
REAL(dp),DIMENSION(SIZE(M,1)):: WORK !for LAPACK routine DGETRI
!
i=0
!
IF( SIZE(M,1).NE.SIZE(M,2) ) THEN
  !Non-square matrix: cancel
  i=1
  !
ELSE
  IF( SIZE(M,1)==3 .AND. SIZE(M,2)==3 ) THEN
    !3x3 matrix: simple enough, let's do it by hand
    det =   M(1,1)*M(2,2)*M(3,3) - M(1,1)*M(3,2)*M(2,3) &
        & - M(2,1)*M(1,2)*M(3,3) + M(2,1)*M(3,2)*M(1,3) &
        & + M(3,1)*M(1,2)*M(2,3) - M(3,1)*M(2,2)*M(1,3)
    !
    G(1,1) = (M(2,2)*M(3,3) - M(2,3)*M(3,2))/det
    G(2,1) = (M(2,3)*M(3,1) - M(2,1)*M(3,3))/det
    G(3,1) = (M(2,1)*M(3,2) - M(2,2)*M(3,1))/det
    !
    G(1,2) = (M(3,2)*M(1,3) - M(3,3)*M(1,2))/det
    G(2,2) = (M(3,3)*M(1,1) - M(3,1)*M(1,3))/det
    G(3,2) = (M(3,1)*M(1,2) - M(3,2)*M(1,1))/det
    !
    G(1,3) = (M(1,2)*M(2,3) - M(1,3)*M(2,2))/det
    G(2,3) = (M(1,3)*M(2,1) - M(1,1)*M(2,3))/det
    G(3,3) = (M(1,1)*M(2,2) - M(1,2)*M(2,1))/det
    !
  ELSEIF( SIZE(G,1)==SIZE(M,1) .AND. SIZE(G,2)==SIZE(M,2) ) THEN
    !general NxN matrix: call LAPACK routines DGETRF and DGETRI
    G(:,:) = M(:,:)
    CALL DGETRF( SIZE(M,1), SIZE(M,2), G, SIZE(M,1), IPIV, i)
    IF( i==0 ) THEN
      CALL DGETRI( SIZE(M,1), G, SIZE(M,1), IPIV, WORK, SIZE(M,1), i )
    ENDIF
    !
  ELSE
    !non-consistent array sizes: some programmer
    !made a mistake when calling this routine
    i=1
  ENDIF
  !
ENDIF
!
IF(PRESENT(status)) status=i
!
END SUBROUTINE INVMAT
!
!
!********************************************************
!  CONVMAT
!  This subroutine converts conventional vectors
!  defined by a b c alpha beta gamma,
!  into a lower triangular matrix:
!        |  H(1,1)    0       0     |
!   H =  |  H(2,1)  H(2,2)    0     |
!        |  H(3,1)  H(3,2)  H(3,3)  |
!********************************************************
!
SUBROUTINE CONVMAT(a,b,c,alpha,beta,gamma,H)
!
IMPLICIT NONE
REAL(dp),INTENT(IN):: a, b, c, alpha, beta, gamma
REAL(dp),DIMENSION(3,3):: H
!
H(:,:) = 0.d0
H(1,1) = a
H(2,1) = b*DCOS(gamma)
H(2,2) = b*DSIN(gamma)
H(3,1) = c*DCOS(beta)
H(3,2) = c*( DSIN(beta)*( DCOS(alpha)-DCOS(beta)*DCOS(gamma) )/ &
       &      (DSIN(beta)*DSIN(gamma))             )
H(3,3) = c*(DSIN(beta)*                                         &
       &     DSQRT(                                             &
       &           ( DSIN(gamma)**2                             &
       &             -DCOS(beta)**2 - DCOS(alpha)**2            &
       &             +2.d0*DCOS(alpha)*DCOS(beta)*DCOS(gamma)   &
       &           )                                            &
       &          )/(DSIN(beta)*DSIN(gamma))                    &
       &   )
!
!
END SUBROUTINE CONVMAT
!
!
!********************************************************
!  MATCONV
!  This subroutine converts a matrix into conventional
!  vectors defined by a b c alpha beta gamma.
!********************************************************
!
SUBROUTINE MATCONV(H,a,b,c,alpha,beta,gamma)
!
IMPLICIT NONE
REAL(dp):: a, b, c, alpha, beta, gamma
REAL(dp),DIMENSION(3,3),INTENT(IN):: H
!
 a = VECLENGTH(H(1,:))
 b = VECLENGTH(H(2,:))
 c = VECLENGTH(H(3,:))
 alpha = ANGVEC( H(2,:),H(3,:) )
  beta = ANGVEC( H(3,:),H(1,:) )
 gamma = ANGVEC( H(1,:),H(2,:) )
!
!
END SUBROUTINE MATCONV
!
!
!********************************************************
!  VOLUME_PARA
!  This subroutine computes the volume of a parallelepiped
!  defined by three vectors stored in a 3x3 matrix.
!********************************************************
!
SUBROUTINE VOLUME_PARA(Pvec,Volume)
!
IMPLICIT NONE
REAL(dp):: a, b, c, alpha, beta, gamma
REAL(dp):: Volume
REAL(dp), DIMENSION(3,3),INTENT(IN):: Pvec
!
CALL MATCONV(Pvec, a, b, c, alpha, beta, gamma)
!
Volume = a*b*c*                                            &
       & (1.d0 + 2.d0*DCOS(alpha)*DCOS(beta)*DCOS(gamma)   &
       &  -DCOS(alpha)**2 -DCOS(beta)**2 -DCOS(gamma)**2  )
!
!
END SUBROUTINE VOLUME_PARA
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
!  INDEX_MILLER
!  This subroutine finds the Miller indices
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
  !Check if we can find A in the array aentries(:,1)
  DO j=1,SIZE(atemp,1)
    IF( DABS(A(i)-atemp(j,1))<1.d-12 ) THEN
      foundA = .TRUE.
      atemp(j,2) = atemp(j,2)+1
    ENDIF
  ENDDO
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
!  DERIVATIVE
!  This subroutine calculates the centered derivative
!  of an array.
!********************************************************
!
SUBROUTINE DERIVATIVE(func,dfunc)
!
IMPLICIT NONE
INTEGER:: i, funcsize
REAL(dp):: step
REAL(dp),DIMENSION(:,:),INTENT(IN):: func
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: dfunc
!
funcsize = SIZE(func(:,1))
IF(.NOT.ALLOCATED(dfunc)) ALLOCATE(dfunc(funcsize-2,2))
dfunc(:,:)=0.d0
!
DO i=2,funcsize-1
  step = func(i+1,1)-func(i-1,1)
  dfunc(i,1) = func(i,1)
  dfunc(i,2) = (func(i+1,2)-func(i-1,2))/step
ENDDO
!
END SUBROUTINE DERIVATIVE
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
!
END MODULE subroutines
