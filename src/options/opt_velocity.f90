MODULE velocity
!
!**********************************************************************************
!*  VELOCITY                                                                      *
!**********************************************************************************
!* This module reads cartesian coordinates from an array and                      *
!* assigns a random velocity to each atom, such that the distribution             *
!* of velocities correspond to a Maxwell-Boltzmann distribution.                  *
!**********************************************************************************
!* (C) July 2013 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 10 Feb. 2022                                     *
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
USE atoms
USE comv
USE constants
USE messages
USE files
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE VELOCITY_XYZ(P,AUXNAMES,AUX,T,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES, newAUXNAMES !names of auxiliary properties
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT  !mask for atom list
INTEGER:: i, j, k
INTEGER:: NP         !Number of atoms which velocity must be set
INTEGER:: vx, vy, vz !columns in AUX where atom velocities are stored
REAL(dp):: mass          !mass of an atom
REAL(dp):: mi, M, A, D, S  !for statistics
REAL(dp),INTENT(IN):: T  !target temperature (K)
REAL(dp):: tempreal
REAL(dp):: velocitymax  !maximum velocity of an atom
REAL(dp),DIMENSION(:),ALLOCATABLE:: randarray    !random numbers
REAL(dp),DIMENSION(:),ALLOCATABLE:: velocities   !velocity of each atom
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX   !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX              !auxiliary properties (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P  !atom positions (required to get atom mass, but won't be modified)
!
!
!Initialize variables
i = 0
!
!
msg = 'Entering VELOCITY_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG( 2111, (/''/), (/ T /) )
!
!
!If target temperature is zero or negative then skip
IF( T<=0.d0 ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2751,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
200 CONTINUE
!Generate a Maxwell-Boltzmann distribution of velocities
!
!Prepare array AUX
IF( ALLOCATED(AUX) ) THEN
  !AUX is already defined,
  !check if it contains vx, vy, vz already
  vx=0
  vy=0
  vz=0
  DO i=1,SIZE(AUXNAMES(:))
    IF( TRIM(ADJUSTL(AUXNAMES(i))) == "vx" ) THEN
      vx = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i))) == "vy" ) THEN
      vy = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i))) == "vz" ) THEN
      vz = i
    ENDIF
  ENDDO
  !If if doesn't contain that property, expand the arrays
  IF( vx==0 .OR. vy==0 .OR. vz==0 ) THEN
    vx = SIZE(AUX,2)+1
    vy = SIZE(AUX,2)+2
    vz = SIZE(AUX,2)+3
    ALLOCATE( newAUX( SIZE(AUX,1), vz ) )
    newAUX(:,:) = 0.d0
    DO i=1,SIZE(AUX,1)
      DO j=1,SIZE(AUX,2)
        newAUX(i,j) = AUX(i,j)
      ENDDO
    ENDDO
    DEALLOCATE(AUX)
    ALLOCATE( AUX( SIZE(newAUX,1), SIZE(newAUX,2) ) )
    AUX(:,:) = newAUX(:,:)
    DEALLOCATE(newAUX)
    !Same with array AUXNAMES
    ALLOCATE(newAUXNAMES( vz ))
    DO i=1,SIZE(AUXNAMES(:))
      newAUXNAMES(i) = AUXNAMES(i)
    ENDDO
    DEALLOCATE(AUXNAMES)
    ALLOCATE(AUXNAMES(vz))
    AUXNAMES(:) = newAUXNAMES(:)
    DEALLOCATE(newAUXNAMES)
  ENDIF
ELSE
  !AUX is not defined yet => create it
  ALLOCATE( AUX( SIZE(P,1) , 3 ) )
  AUX(:,:) = 0.d0
  ALLOCATE( AUXNAMES(3) )
  vx = 1
  vy = 2
  vz = 3
ENDIF
WRITE(msg,*) "Index vx, vy, vz: ", vx, vy, vz
CALL ATOMSK_MSG(999,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
AUXNAMES(vx) = "vx"
AUXNAMES(vy) = "vy"
AUXNAMES(vz) = "vz"
!
!Count how many atoms are selected
IF( ALLOCATED(SELECT) ) THEN
  NP=0
  DO i=1,SIZE(SELECT)
    IF(SELECT(i)) NP=NP+1
  ENDDO
ELSE
  NP = SIZE(P,1)
ENDIF
!
!Generate 3*NP random numbers
CALL GEN_NRANDGAUSS(3*NP,randarray)
!randarray now contains 3*NP random real numbers with a gaussian distribution
!Save random numbers as auxiliary properties in array AUX
!Multiply them by sqrt(k*T/m) to obtain proper Maxwell-Boltzmann distribution
!(k=Boltzmann constant (J/K); T=target temperature (K); m=mass of atom (a.u.))
!Multiply by 1.d-10/1.d-12 = 1.d2 to obtain A/ps.
j=0
DO i=1,SIZE(AUX,1)
  IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
    j=j+1
    CALL ATOMSPECIES( P(i,4),species )
    CALL ATOMMASS( species,mass )
    IF( NINT(mass)==0 ) mass=1.d0  !avoid division by zero
    mass = mass * mass_amu  !conversion to kg
    AUX(i,vx) = randarray(j)      * DSQRT(kB*T/mass) * 1.d2
    AUX(i,vy) = randarray(NP+j)   * DSQRT(kB*T/mass) * 1.d2
    AUX(i,vz) = randarray(2*NP+j) * DSQRT(kB*T/mass) * 1.d2
  ENDIF
ENDDO
!
!At this point all velocities are positive.
!Apply a correction so that there is no global drift of the system
!(i.e. sum(vx)=0, sum(vy)=0 and sum(vz)=0)
IF(verbosity==4) THEN
  WRITE(msg,*) "Sum of velocities along X, Y, Z: ", SUM(AUX(:,vx)), SUM(AUX(:,vy)), SUM(AUX(:,vz))
  CALL ATOMSK_MSG(999,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
ENDIF
IF( .NOT.ALLOCATED(SELECT) ) THEN
  !No selection defined => rescale velocities of all atoms
  IF( SUM(AUX(:,vx)) .NE. 0.d0 ) THEN
    AUX(i,vx) = AUX(i,vx) - SUM(AUX(:,vx)) / DBLE(SIZE(AUX,1))
  ENDIF
  IF( SUM(AUX(:,vy)) .NE. 0.d0 ) THEN
    AUX(i,vy) = AUX(i,vy) - SUM(AUX(:,vy)) / DBLE(SIZE(AUX,1))
  ENDIF
  IF( SUM(AUX(:,vz)) .NE. 0.d0 ) THEN
    AUX(i,vz) = AUX(i,vz) - SUM(AUX(:,vz)) / DBLE(SIZE(AUX,1))
  ENDIF
  !
ELSE
  !A selection is defined => rescale only velocities of selected atoms
  tempreal = 0.d0
  DO i=1,SIZE(AUX,1)
    IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
      tempreal = tempreal + AUX(i,vx)
    ENDIF
  ENDDO
  tempreal = tempreal / NP
  IF( DABS(tempreal) > 1.d-12 ) THEN
    DO i=1,SIZE(AUX,1)
      IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
        AUX(i,vx) = AUX(i,vx) - tempreal
      ENDIF
    ENDDO
  ENDIF
  !
  tempreal = 0.d0
  DO i=1,SIZE(AUX,1)
    IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
      tempreal = tempreal + AUX(i,vy)
    ENDIF
  ENDDO
  tempreal = tempreal / NP
  IF( DABS(tempreal) > 1.d-12 ) THEN
    DO i=1,SIZE(AUX,1)
      IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
        AUX(i,vy) = AUX(i,vy) - tempreal
      ENDIF
    ENDDO
  ENDIF
  !
  tempreal = 0.d0
  DO i=1,SIZE(AUX,1)
    IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
      tempreal = tempreal + AUX(i,vz)
    ENDIF
  ENDDO
  tempreal = tempreal / NP
  IF( DABS(tempreal) > 1.d-12 ) THEN
    DO i=1,SIZE(AUX,1)
      IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
        AUX(i,vz) = AUX(i,vz) - tempreal
      ENDIF
    ENDDO
  ENDIF
  !
ENDIF
!
!Store velocities in array (this is done for ALL atoms)
ALLOCATE( velocities(SIZE(AUX,1)) )
velocities(:) = 0.d0
DO i=1,SIZE(velocities)
  velocities(i) = DSQRT( AUX(i,vx)**2 + AUX(i,vy)**2 + AUX(i,vz)**2 )
ENDDO
velocitymax = MAXVAL(velocities)
WRITE(msg,*) "Velocity max: ", velocitymax
CALL ATOMSK_MSG(999,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
!
CALL ATOMSK_MSG(2112,(/''/),(/0.d0/))
!
!
!
300 CONTINUE
!Output properties of velocity distribution in a file
msg = "velocity_dist.dat"
IF(.NOT.overw) CALL CHECKFILE(msg,"write")
OPEN(UNIT=35,FILE=msg,FORM="FORMATTED",STATUS="UNKNOWN")
CALL DO_STATS(velocities(:),mi,M,A,D,S)
WRITE(temp,'(f16.3)') T
WRITE(35,'(a)') "# Velocity (A/ps) distribution generated by Atomsk for T = "//TRIM(ADJUSTL(temp))
WRITE(35,'(a26,e12.3)') "# Min. velocity:          ", mi
WRITE(35,'(a26,e12.3)') "# Max. velocity:          ", M
WRITE(35,'(a26,e12.3)') "# Average velocity:       ", A
WRITE(35,'(a26,e12.3)') "# Av. absolute deviation: ", D
WRITE(35,'(a26,e12.3)') "# Standard deviation:     ", S
NP=0
DO j=0,100
  k=0
  DO i=1,SIZE(velocities)
    IF( velocities(i) >  DBLE(j)*velocitymax/100.d0 .AND.                &
      & velocities(i) <= DBLE(j+1)*velocitymax/100.d0    ) THEN
      k=k+1
    ENDIF
  ENDDO
  WRITE(35,'(2e24.4)') DBLE(j)*velocitymax/100.d0, DBLE(k)
  NP=NP+k
ENDDO
CLOSE(35)
CALL ATOMSK_MSG(2113,(/msg/),(/0.d0/))
!
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
END SUBROUTINE VELOCITY_XYZ
!
!
END MODULE velocity
