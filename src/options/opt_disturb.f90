MODULE disturb
!
!**********************************************************************************
!*  DISTURB                                                                       *
!**********************************************************************************
!* This module reads cartesian coordinates from an array and                      *
!* applies a random shift to atom positions.                                      *
!**********************************************************************************
!* (C) July 2013 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 31 May 2023                                      *
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
USE random
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE DISTURB_XYZ(P,S,SELECT,dmode,dmax)
!
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT  !mask for atom list
INTEGER,INTENT(IN):: dmode  !0=user provided values for dmax along x, y, z (default)
                            !Otherwise user provided norm of max. displacement
INTEGER:: i, j
INTEGER:: NP, Nmoved
REAL(dp):: drift  !total displacement along a direction
REAL(dp),DIMENSION(3):: dmax  !maximum displacement of an atom along each cartesian direction
REAL(dp),DIMENSION(:),ALLOCATABLE:: randarray    !random numbers
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S  !positions of cores, shells
!
!
!Initialize variables
i = 0
Nmoved = 0
!
!
WRITE(msg,*) 'Entering DISTURB_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
IF( dmode==0 ) THEN
  CALL ATOMSK_MSG( 2114, (/"xyz"/), (/ dmax(1),dmax(2),dmax(3) /) )
ELSE
  CALL ATOMSK_MSG( 2114, (/''/), (/ dmax(1),0.d0,0.d0 /) )
  !User provided norm of maximal displacement => compute corresponding max. displacements along X, Y, Z
  drift = DSQRT((dmax(1)**2)/3.d0)
  dmax(:) = drift
ENDIF
!
WRITE(msg,*) 'dmax(:) = ', dmax(1), dmax(2), dmax(3)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!If maximum displacement is zero then skip
IF( DABS(dmax(1))<1.d-12 .AND. DABS(dmax(2))<1.d-12 .AND. DABS(dmax(3))<1.d-12 ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2751,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
200 CONTINUE
!Count how many atoms must be disturbed
IF( .NOT.ALLOCATED(SELECT) ) THEN
  NP = SIZE(P,1)
ELSE
  NP=0
  DO i=1,SIZE(SELECT)
    IF(SELECT(i)) NP=NP+1
  ENDDO
ENDIF
WRITE(msg,*) 'NP = ', NP
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Generate 3N random numbers
CALL GEN_NRANDNUMBERS(3*NP,randarray)
!randarray now contains 3N random real numbers uniformly distributed between 0.d0 and 1.d0
!Scale it so that it is between -dxmax and +dxmax
randarray(:) = ( randarray(:) - 0.5d0 ) * 2.d0
!randarray(1:SIZE(P,1)) contains displacements along X
!randarray(SIZE(P,1)+1,2*SIZE(P,1)) contains displacements along Y
!randarray(2*SIZE(P,1)+1,3*SIZE(P,1)) contains displacements along Z
!At this point all values in randarray(:) are positive.
!Adjust them so that total displacement in each direction is zero.
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,drift)
DO i=1,3 !x,y,z
  !Compute global drift along that direction
  drift = SUM( randarray( (i-1)*NP+1 : i*NP) ) * dmax(i)
  !Apply correction
  DO j=1,NP
    randarray((i-1)*NP+j) = randarray((i-1)*NP+j) * dmax(i) - drift/DBLE(SIZE(randarray))
  ENDDO
ENDDO
!$OMP END PARALLEL DO
!
IF( verbosity>=4 ) THEN
  WRITE(msg,*) "randarray (SIZE ", SIZE(randarray), ") ="
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,MIN(SIZE(randarray/3),10)
    WRITE(msg,'(6X,3f12.6,2X)') randarray(i), randarray(SIZE(P,1)+i), randarray(2*SIZE(P,1)+i)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
  IF( SIZE(randarray/3)>10 ) THEN
    msg = "      (...continued...)"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
ENDIF
!
!
!Use random numbers as displacements along X, Y, Z
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) REDUCTION(+:Nmoved)
DO i=1,SIZE(P,1)
  IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
    P(i,1) = P(i,1) + randarray(i)
    P(i,2) = P(i,2) + randarray(NP+i)
    P(i,3) = P(i,3) + randarray(2*NP+i)
    IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
      S(i,1) = S(i,1) + randarray(i)
      S(i,2) = S(i,2) + randarray(NP+i)
      S(i,3) = S(i,3) + randarray(2*NP+i)
    ENDIF
    Nmoved = Nmoved+1
  ENDIF
ENDDO
!$OMP END PARALLEL DO
!
CALL ATOMSK_MSG(2115,(/''/),(/DBLE(Nmoved)/))
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
!
!
!
1000 CONTINUE
IF(ALLOCATED(randarray)) DEALLOCATE(randarray)
!
!
END SUBROUTINE DISTURB_XYZ
!
!
END MODULE disturb
