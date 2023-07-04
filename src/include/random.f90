MODULE random
!
!**********************************************************************************
!*  RANDOM                                                                        *
!**********************************************************************************
!* This module contains procedures generating random numbers.                     *
!**********************************************************************************
!* (C) April 2010 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 03 July 2023                                     *
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
!* GEN_NRANDNUMBERS    generate a list of N random real numbers in [0.0,1.0]      *
!* GEN_NRANDGAUSS      generate a list of N numbers with gaussian distribution    *
!* GEN_NRANDINDEX      generate a list of int. numbers sorted randomly            *
!**********************************************************************************
!
!
USE comv
USE constants
USE math
!
IMPLICIT NONE
!
!
CONTAINS
!
!
!********************************************************
! GEN_NRANDNUMBERS
! This subroutine generates N "random" real numbers
! between 0 and 1 with a uniform distribution,
! and returns them in the array randarray(:).
! On GNU/Linux and BSD or macOS systems, it tries to
! use the file /dev/urandom as a seed, which should
! guaranty to have a "true random" seed.
! On Windows systems, or if the file /dev/urandom can't
! be read, it uses the system clock as a seed,
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
LOGICAL:: fileexists
INTEGER:: i, k
INTEGER(KIND=8):: clock !system clock
INTEGER,INTENT(IN):: N  !number of random numbers to generate
INTEGER,DIMENSION(:),ALLOCATABLE:: seed !seed for generating random numbers
REAL(dp),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: randarray    !random numbers
!
fileexists=.FALSE.
!
CALL RANDOM_SEED(SIZE=k)
IF(ALLOCATED(seed)) DEALLOCATE(seed)
ALLOCATE(seed(k))
seed(:) = 0
!
#if ! defined(WINDOWS)
!Detect if file "/dev/urandom" exists on current system
!(NOTE: we know that this file does not exist on Windows systems, hence the preproc #if instruction)
INQUIRE(FILE="/dev/urandom",EXIST=fileexists)
IF( fileexists ) THEN
  !Use the file /dev/urandom to generate random seed
  OPEN(50,FILE="/dev/urandom",ACCESS="stream",FORM='UNFORMATTED')
  READ(50,ERR=10,END=10) seed
  10 CONTINUE
  CLOSE(50)
  !If seed(:) is still filled with zeros, there was a problem
  IF( .NOT.ANY(seed(:).NE.0) ) fileexists=.FALSE.
ENDIF
#endif
IF( .NOT.fileexists ) THEN
  !Use system clock as seed
  CALL SYSTEM_CLOCK(COUNT=clock)
  seed = clock + 42*(/ (i-1, i=1,k) /)
ENDIF
!
CALL RANDOM_SEED(PUT=seed)
DEALLOCATE(seed)
!Use seed to generate N numbers
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
!********************************************************
SUBROUTINE GEN_NRANDGAUSS(N,randarray)
!
IMPLICIT NONE
INTEGER:: i
INTEGER,INTENT(IN):: N  !number of random numbers to generate
REAL(dp),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: randarray    !final random numbers
REAL(dp),DIMENSION(:),ALLOCATABLE:: randarraytemp            !random numbers
!
!Generate 2N random numbers with uniform distribution
CALL GEN_NRANDNUMBERS(2*N,randarraytemp)
!
!Generate the final gaussian distribution
IF( ALLOCATED(randarray) ) DEALLOCATE(randarray)
ALLOCATE(randarray(N))
randarray(:) = 0.d0
DO i=1,N
  randarray(i) = DSQRT( -2.0d0*DLOG(randarraytemp(i)) ) * DCOS( 2.0d0*pi*randarraytemp(N+i) )
ENDDO
!
IF( ALLOCATED(randarraytemp) ) DEALLOCATE(randarraytemp)
!
END SUBROUTINE GEN_NRANDGAUSS
!
!
!********************************************************
! GEN_NRANDINDEX
! This subroutine generates a list of N atom indices
! in a random order.
!********************************************************
SUBROUTINE GEN_NRANDINDEX(N,idlist)
!
INTEGER:: i, j, k
INTEGER,INTENT(IN):: N  !number of index to generate
INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(OUT):: idlist !random indices
REAL(dp),DIMENSION(:),ALLOCATABLE:: randarray    !random numbers from 0 to 1
!
IF(ALLOCATED(idlist)) DEALLOCATE(idlist)
!
IF( N >= 1 ) THEN
  !Initialize the list with sorted values from 1 to N
  ALLOCATE(idlist(N))
  DO i=1,N
    idlist(i) = i
  ENDDO
  !
  IF( N>1 ) THEN
    !Generate N random numbers
    CALL GEN_NRANDNUMBERS( N , randarray )
    !
    !Use the random list to sort idlist
    DO i=1,N
      !Pick a random element j
      j = FLOOR(N*randarray(i))
      IF(j==0) j=1
      IF(j>N) j=N
      IF( i.NE.j ) THEN
        !Exchange elements i and j
        k = idlist(i)
        idlist(i) = idlist(j)
        idlist(j) = k
      ENDIF
    ENDDO
    !
  ENDIF
  !
ENDIF
!
END SUBROUTINE GEN_NRANDINDEX
!
!
!
END MODULE random
