MODULE roll
!
!**********************************************************************************
!*  ROLL                                                                          *
!**********************************************************************************
!* This module reads cartesian coordinates from an array P and "rolls"            *
!* atoms around the given direction. This option is best suited for               *
!* pseudo-2D systems, e.g. to make nanotubes.                                     *
!**********************************************************************************
!* (C) January 2016 - Pierre Hirel                                                *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 22 Jan. 2016                                     *
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
SUBROUTINE ROLL_XYZ(H,P,S,AUXNAMES,AUX,direction,angle,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=1),INTENT(IN):: direction
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: a1, a2, a3
INTEGER:: i
REAL(dp):: radius   !radius of the final tube
REAL(dp):: angle    !total angle of rolled system
REAL(dp):: r, theta !actual position of an atom (cylindre coordinates)
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX  !auxiliary properties of atoms/shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of atoms, shells
!
!Initialize variables
i = 0
a1=0
a2=0
a3=0
radius = 10.d0
!
!
msg = 'Entering ROLL_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(2127,(/direction/),(/angle/))
!
IF( DABS(angle)<=1.d-12) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2734,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!convert the angle into radians
!PRINT*, "angle (degrees) = ", angle
angle = DEG2RAD(angle)
!PRINT*, "angle (radians) = ", angle
!
!Define the axes
!a3 will be the axis of the tube, i.e. only a1 and a2 coordinates will be changed
IF(direction=='x' .OR. direction=='X') THEN
  a1 = 2
  a2 = 3
  a3 = 1
ELSEIF(direction=='y' .OR. direction=='Y') THEN
  a1 = 1
  a2 = 3
  a3 = 2
ELSEIF(direction=='z' .OR. direction=='Z') THEN
  a1 = 1
  a2 = 2
  a3 = 3
ELSE
  nerr = nerr+1
  CALL ATOMSK_MSG(2800,(/direction/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Determine the radius of the tube
!radius = DSQRT( VECLENGTH(H(:,a1)) / pi )
radius = VECLENGTH(H(:,a1)) / angle
!PRINT*, "radius = ", radius
!
!
!
100 CONTINUE
!Now let's roll the sheet into a tube
DO i=1,SIZE(P,1)
  r = radius - P(i,a2)
  theta = angle*P(i,a1)/VECLENGTH(H(:,a1)) - angle
  !PRINT*, i, r, theta*180.d0/pi
  P(i,a1) = r*DCOS(theta)
  P(i,a2) = r*DSIN(theta)
ENDDO
!
!
!Modify supercell vectors to fit the tube
H(a1,:) = 0.d0
H(a1,a1) = 2.d0*radius
H(a2,:) = 0.d0
H(a2,a2) = 2.d0*radius
!
CALL ATOMSK_MSG(2128,(/''/),(/0.d0/))
!
!
!
1000 CONTINUE
!
!
!
END SUBROUTINE ROLL_XYZ
!
!
END MODULE roll