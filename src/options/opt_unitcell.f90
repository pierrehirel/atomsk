MODULE unitcell
!
!**********************************************************************************
!*  UNITCELL                                                                      *
!**********************************************************************************
!* This module reads an array of type (species x y z) and the crystal             *
!* orientation, and finds the unit cell for that orientation.                     *
!**********************************************************************************
!* (C) April 2011 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@izbs.uni-karlsruhe.de                                         *
!* Last modification: P. Hirel - 12 July 2012                                     *
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
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE UNITCELL_XYZ(H,P,lat_a0,ORIENT)
!
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
INTEGER:: i, m, n, o
INTEGER:: Nuc  !Number of atoms in unit cell
INTEGER:: Rsize
INTEGER,DIMENSION(3):: expandmatrix         !
REAL(dp),DIMENSION(3),INTENT(IN):: lat_a0   !lattice constants a, b, c
REAL(dp),DIMENSION(3):: uv                  !unit vector along each direction
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3),INTENT(IN):: ORIENT !crystal orientation
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, R
!
!Initialize variables
expandmatrix(1)=1
expandmatrix(2)=1
expandmatrix(3)=1
Rsize = SIZE(P(:,1))
uv(:) = 0.d0
IF(ALLOCATED(R)) DEALLOCATE(R)
!
WRITE(msg,*) 'debug --> Entering UNITCELL_XYZ'
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = '>>> Reducing to unit cell...'
CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!
!
100 CONTINUE
!Find the distance between 2 planes along each direction:
!    d_hkl = a0 / sqrt(h^2+k^2+l^2)
DO i=1,3
  uv(i) = 2.d0*lat_a0(i) /                        &
        & DSQRT(ORIENT(i,1)**2.d0 +               &
        &       ORIENT(i,2)**2.d0 +               &
        &       ORIENT(i,3)**2.d0                 &
        &      )
ENDDO
!
!
!
200 CONTINUE
!Shift and expand system so that it fills the unit cell
!1- shift the system so that everything is below (0,0,0)
IF( H(2,1).NE.0.d0 .OR. H(3,1).NE.0.d0 ) THEN
  DO i=1,SIZE(P(:,1))
    P(i,1) = P(i,1)-H(2,1)-H(3,1)
  ENDDO
ENDIF
!
!2- determine the number of times the system must be repeated
DO i=1,3
  IF( H(i,1)<uv(1) ) THEN
    expandmatrix(1)=expandmatrix(1)+1
  ENDIF
ENDDO
DO i=1,3
  IF( H(i,2)<uv(2) ) THEN
    expandmatrix(2)=expandmatrix(2)+1
  ENDIF
ENDDO
DO i=1,3
  IF( H(i,3)<uv(3) ) THEN
    expandmatrix(3)=expandmatrix(3)+1
  ENDIF
ENDDO
!
Rsize = SIZE(P(:,1))*expandmatrix(1)*expandmatrix(2)*expandmatrix(3)
ALLOCATE(R(Rsize,4))
!
!3- duplicate system and save it in R
IF( expandmatrix(1)>1 .OR. expandmatrix(2)>1 .OR. expandmatrix(3)>1 ) THEN
  Nuc=0
  DO i=1,SIZE(P(:,1))
    DO o=1,expandmatrix(3)
      DO n=1,expandmatrix(2)
        DO m=1,expandmatrix(1)
          Nuc = Nuc+1
          R(Nuc,1) = P(i,1)+REAL(m-1)*H(1,1)+REAL(n-1)*H(2,1)+REAL(o-1)*H(3,1)
          R(Nuc,2) = P(i,2)+REAL(m-1)*H(1,2)+REAL(n-1)*H(2,2)+REAL(o-1)*H(3,2)
          R(Nuc,3) = P(i,3)+REAL(m-1)*H(1,3)+REAL(n-1)*H(2,3)+REAL(o-1)*H(3,3)
          R(Nuc,4) = P(i,4)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDIF
!
!
!
300 CONTINUE
DEALLOCATE(P)
!Count atoms that are actually in the box
Nuc=0
DO i=1,SIZE(R(:,1))
  IF( R(i,1)>=0.d0 .AND. R(i,1)<uv(1) .AND.             &
    & R(i,2)>=0.d0 .AND. R(i,2)<uv(2) .AND.             &
    & R(i,3)>=0.d0 .AND. R(i,3)<uv(3)        ) THEN
    Nuc = Nuc+1
  ENDIF
ENDDO
!
!If no atom was found we are in trouble
IF(Nuc==0) THEN
  nerr=nerr+1
  WRITE(msg,*) 'X!X ERROR: no atom found in the unit cell.'
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  GOTO 1000
ENDIF
!
!Save unit cell atoms to P
ALLOCATE(P(Nuc,4))
DO i=1,SIZE(R(:,1))
  IF( R(i,1)>=0.d0 .AND. R(i,1)<uv(1) .AND.             &
    & R(i,2)>=0.d0 .AND. R(i,2)<uv(2) .AND.             &
    & R(i,3)>=0.d0 .AND. R(i,3)<uv(3)        ) THEN
    Nuc = Nuc+1
    P(Nuc,:) = R(i,:)
  ENDIF
ENDDO
!
!
!
400 CONTINUE
!Save new cell vectors to H
H(:,:) = 0.d0
DO i=1,3
  H(i,i) = uv(i)
ENDDO
!
!
msg = '..> System was reduced to unit cell.'
CALL DISPLAY_MSG(verbosity,msg,logfile)
GOTO 1000
!
!
!
800 CONTINUE
WRITE(msg,*) 'X!X ERROR while trying to read atom ', i
CALL DISPLAY_MSG(verbosity,msg,logfile)
nerr = nerr+1
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE UNITCELL_XYZ
!
!
!
END MODULE unitcell