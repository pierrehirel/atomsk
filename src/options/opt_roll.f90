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
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 16 June 2016                                     *
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
SUBROUTINE ROLL_XYZ(H,P,S,AUXNAMES,AUX,direction,angle,axis,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=1),INTENT(IN):: direction, axis
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
a1 = 0
a2 = 0
a3 = 0
radius = 10.d0
!
!
msg = 'Entering ROLL_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(2127,(/direction,axis/),(/angle/))
!
IF( DABS(angle)<=1.d-12) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2734,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!convert the angle into radians
angle = DEG2RAD(angle)
!
!Define the axes
!a3 is the axis of the tube
!a1 is the direction that is rolled around a3
IF(axis=='x' .OR. axis=='X') THEN
  a3 = 1
  IF(direction=='y' .OR. direction=='Y') THEN
    a1 = 2
    a2 = 3
  ELSE !i.e. direction=="Z"
    a1 = 3
    a2 = 2
  ENDIF
ELSEIF(axis=='y' .OR. axis=='Y') THEN
  a3 = 2
  IF(direction=='x' .OR. direction=='X') THEN
    a1 = 1
    a2 = 3
  ELSE !i.e. direction=="Z"
    a1 = 3
    a2 = 1
  ENDIF
ELSEIF(axis=='z' .OR. axis=='Z') THEN
  a3 = 3
  IF(direction=='x' .OR. direction=='X') THEN
    a1 = 1
    a2 = 2
  ELSE !i.e. direction=="Y"
    a1 = 2
    a2 = 1
  ENDIF
ELSE
  nerr = nerr+1
  CALL ATOMSK_MSG(2800,(/axis/),(/0.d0/))
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
  P(i,a1) = r*DCOS(theta)
  P(i,a2) = r*DSIN(theta)
ENDDO
!
!
!Modify supercell vectors along a3 to fit the bent system
!If angle > pi then the box dimension will be increased by 2*radius
H(a2,a2) = H(a2,a2) + radius* ( 1.d0 - DCOS(MIN(angle,pi)) )
!
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