MODULE torsion
!
!**********************************************************************************
!*  TORSION                                                                       *
!**********************************************************************************
!* This module reads cartesian coordinates from an array P, and                   *
!* applies a torsion of the given angle around the given axis.                    *
!**********************************************************************************
!* (C) January 2016 - Pierre Hirel                                                *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
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
SUBROUTINE TORSION_XYZ(H,P,S,AUXNAMES,AUX,direction,angle,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=1),INTENT(IN):: direction
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
!
CHARACTER(LEN=128):: msg
LOGICAL:: velocities, forces
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: a1, a2, a3
INTEGER:: i
INTEGER,DIMENSION(3):: Fxyz, Vxyz !columns of AUX containing forces, velocities
REAL(dp):: angle  !angle of torsion (in degrees)
REAL(dp):: P2, P3
REAL(dp):: rot_angle
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX  !auxiliary properties of atoms/shells
REAL(dp),DIMENSION(3,3):: rot_matrix
!
!Initialize variables
forces=.FALSE.
velocities=.FALSE.
i = 0
Fxyz(:)=0
Vxyz(:)=0
P2 = 0.d0
P3 = 0.d0
rot_matrix(:,:) = 0.d0
DO i=1,3
  rot_matrix(i,i) = 1.d0
ENDDO
!
!
msg = 'Entering TORSION_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(2129,(/direction/),(/angle/))
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
!a1 will be the axis of rotation, i.e. only a2 and a3 coordinates will be changed
IF(direction=='x' .OR. direction=='X') THEN
  a1 = 1
  a2 = 2
  a3 = 3
ELSEIF(direction=='y' .OR. direction=='Y') THEN
  a1 = 2
  a2 = 3
  a3 = 1
ELSEIF(direction=='z' .OR. direction=='Z') THEN
  a1 = 3
  a2 = 1
  a3 = 2
ELSE
  nerr = nerr+1
  CALL ATOMSK_MSG(2800,(/direction/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Parse auxiliary properties, search for known vectors
IF( ALLOCATED(AUXNAMES) .AND. SIZE(AUXNAMES)>0 ) THEN
  DO i=1,SIZE(AUXNAMES)
    IF( AUXNAMES(i)=="fx" ) THEN
      Fxyz(1)=i
    ELSEIF( AUXNAMES(i)=="fy" ) THEN
      Fxyz(2)=i
    ELSEIF( AUXNAMES(i)=="fz" ) THEN
      Fxyz(3)=i
    ELSEIF( AUXNAMES(i)=="vx" ) THEN
      Vxyz(1)=i
    ELSEIF( AUXNAMES(i)=="vy" ) THEN
      Vxyz(2)=i
    ELSEIF( AUXNAMES(i)=="vz" ) THEN
      Vxyz(3)=i
    ENDIF
  ENDDO
  !
  IF( Fxyz(1)>0 .AND. Fxyz(2)>0 .AND. Fxyz(3)>0 ) THEN
    forces = .TRUE.
  ENDIF
  IF( Vxyz(1)>0 .AND. Vxyz(2)>0 .AND. Vxyz(3)>0 ) THEN
    velocities = .TRUE.
  ENDIF
ENDIF
!
!
!
100 CONTINUE
!Apply torsion
DO i=1,SIZE(P,1)
  IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
    !Actual angle of rotation for current atom
    rot_angle = angle*P(i,a1) / VECLENGTH(H(:,a1))
    !
    !set the rotation matrix
    rot_matrix(:,:) = 0.d0
    rot_matrix(a1,a1) = 1.d0
    rot_matrix(a2,a2) = DCOS(rot_angle)
    rot_matrix(a2,a3) = -DSIN(rot_angle)
    rot_matrix(a3,a2) = DSIN(rot_angle)
    rot_matrix(a3,a3) = DCOS(rot_angle)
    !
    !Apply rotation to this atom
    P2 = P(i,a2) - 0.5d0*VECLENGTH(H(:,a2))
    P3 = P(i,a3) - 0.5d0*VECLENGTH(H(:,a3))
    P(i,a2) = P2*rot_matrix(a2,a2) + P3*rot_matrix(a2,a3) + 0.5d0*VECLENGTH(H(:,a2))
    P(i,a3) = P2*rot_matrix(a3,a2) + P3*rot_matrix(a3,a3) + 0.5d0*VECLENGTH(H(:,a3))
    !
    !Same with shell if they exist
    IF( ALLOCATED(S) .AND. SIZE(S,1)>0 ) THEN
      P2 = S(i,a2) - 0.5d0*VECLENGTH(H(:,a2))
      P3 = S(i,a3) - 0.5d0*VECLENGTH(H(:,a3))
      S(i,a2) = P2*rot_matrix(a2,a2) + P3*rot_matrix(a2,a3) + 0.5d0*VECLENGTH(H(:,a2))
      S(i,a3) = P2*rot_matrix(a3,a2) + P3*rot_matrix(a3,a3) + 0.5d0*VECLENGTH(H(:,a3))
    ENDIF
    !Same with forces if they exist
    IF( forces ) THEN
      P2 = AUX(i,Fxyz(a2))
      P3 = AUX(i,Fxyz(a3))
      AUX(i,Fxyz(a2)) = P2*rot_matrix(a2,a2) + P3*rot_matrix(a2,a3)
      AUX(i,Fxyz(a3)) = P2*rot_matrix(a3,a2) + P3*rot_matrix(a3,a3)
    ENDIF
    !Same with velocities if they exist
    IF( velocities ) THEN
      P2 = AUX(i,Vxyz(a2))
      P3 = AUX(i,Vxyz(a3))
      AUX(i,Vxyz(a2)) = P2*rot_matrix(a2,a2) + P3*rot_matrix(a2,a3)
      AUX(i,Vxyz(a3)) = P2*rot_matrix(a3,a2) + P3*rot_matrix(a3,a3)
    ENDIF
  ENDIF
ENDDO
!
CALL ATOMSK_MSG(2130,(/''/),(/0.d0/))
!
!
!
1000 CONTINUE
!
!
!
END SUBROUTINE TORSION_XYZ
!
!
END MODULE torsion