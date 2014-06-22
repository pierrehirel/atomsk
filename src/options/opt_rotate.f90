MODULE rotate
!
!**********************************************************************************
!*  ROTATE                                                                        *
!**********************************************************************************
!* This module reads atomic positions from an array P and rotates                 *
!* the system by a certain angle around the cartesian X, Y or Z axis.             *
!**********************************************************************************
!* (C) May 2010 - Pierre Hirel                                                    *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 10 Feb. 2014                                     *
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
SUBROUTINE ROTATE_XYZ(H,P,S,AUXNAMES,AUX,rot_axis,rot_angle,SELECT,C_tensor)
!
!
IMPLICIT NONE
CHARACTER(LEN=1),INTENT(IN):: rot_axis  ! cartesian x, y or z axis
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
INTEGER:: a1, a2, a3
INTEGER:: i
INTEGER,DIMENSION(3):: Fxyz, Vxyz !columns of AUX containing forces, velocities
LOGICAL:: velocities, forces
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
REAL(dp):: H2, H3
REAL(dp):: rot_angle      !angle in degrees
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: rot_matrix
REAL(dp),DIMENSION(9,9),INTENT(INOUT):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of cores, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX  !auxiliary properties of atoms/shells
!
!Initialize variables
forces=.FALSE.
velocities=.FALSE.
i = 0
Fxyz(:)=0
Vxyz(:)=0
H2 = 0.d0
H3 = 0.d0
rot_matrix(:,:) = 0.d0
DO i=1,3
  rot_matrix(i,i) = 1.d0
ENDDO
!
msg = 'Entering ROTATE_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(2081,(/rot_axis/),(/rot_angle/))
!
!If angle is more than 360°, reduce it
DO WHILE(DABS(rot_angle)>=360.d0)
  IF(rot_angle>=360.d0)  rot_angle = rot_angle-360.d0
  IF(rot_angle<=-360.d0) rot_angle = rot_angle+360.d0
ENDDO
!
IF(rot_angle==0.d0) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2734,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!convert the angle into radians
rot_angle = DEG2RAD(rot_angle)
!
!Define the axes
IF(rot_axis=='x' .OR. rot_axis=='X') THEN
  a1 = 1
  a2 = 2
  a3 = 3
ELSEIF(rot_axis=='y' .OR. rot_axis=='Y') THEN
  a1 = 2
  a2 = 3
  a3 = 1
ELSEIF(rot_axis=='z' .OR. rot_axis=='Z') THEN
  a1 = 3
  a2 = 1
  a3 = 2
ELSE
  CALL ATOMSK_MSG(2800,(/rot_axis/),(/0.d0/))
  nerr = nerr+1
  GOTO 1000
ENDIF
!
!set the rotation matrix
rot_matrix(:,:) = 0.d0
rot_matrix(a1,a1) = 1.d0
rot_matrix(a2,a2) = DCOS(rot_angle)
rot_matrix(a2,a3) = -DSIN(rot_angle)
rot_matrix(a3,a2) = DSIN(rot_angle)
rot_matrix(a3,a3) = DCOS(rot_angle)
!
IF(verbosity==4) THEN
  msg = 'debug --> Rotation matrix:'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,60) '           | ',                              &
              & rot_matrix(1,1), rot_matrix(1,2), rot_matrix(1,3), '|'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,60) ' R_theta = | ',                              &
              & rot_matrix(2,1), rot_matrix(2,2), rot_matrix(2,3), '|'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,60) '           | ',                              &
              & rot_matrix(3,1), rot_matrix(3,2), rot_matrix(3,3), '|'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  60 FORMAT(a13,3(f10.5,2X),a1)
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
!Rotate the atomic positions
!P(a1) is untouched and we rotate only atoms that are selected in SELECT
DO i=1,SIZE(P,1)
  IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
    H2 = P(i,a2)
    H3 = P(i,a3)
    P(i,a2) = H2*rot_matrix(a2,a2) + H3*rot_matrix(a2,a3)
    P(i,a3) = H2*rot_matrix(a3,a2) + H3*rot_matrix(a3,a3)
    !Same with shell if they exist
    IF( ALLOCATED(S) .AND. SIZE(S,1)>0 ) THEN
      H2 = S(i,a2)
      H3 = S(i,a3)
      S(i,a2) = H2*rot_matrix(a2,a2) + H3*rot_matrix(a2,a3)
      S(i,a3) = H2*rot_matrix(a3,a2) + H3*rot_matrix(a3,a3)
    ENDIF
    !Same with forces if they exist
    IF( forces ) THEN
      H2 = AUX(i,Fxyz(a2))
      H3 = AUX(i,Fxyz(a3))
      AUX(i,Fxyz(a2)) = H2*rot_matrix(a2,a2) + H3*rot_matrix(a2,a3)
      AUX(i,Fxyz(a3)) = H2*rot_matrix(a3,a2) + H3*rot_matrix(a3,a3)
    ENDIF
    !Same with velocities if they exist
    IF( velocities ) THEN
      H2 = AUX(i,Vxyz(a2))
      H3 = AUX(i,Vxyz(a3))
      AUX(i,Vxyz(a2)) = H2*rot_matrix(a2,a2) + H3*rot_matrix(a2,a3)
      AUX(i,Vxyz(a3)) = H2*rot_matrix(a3,a2) + H3*rot_matrix(a3,a3)
    ENDIF
  ENDIF
ENDDO
!
IF( .NOT.ALLOCATED(SELECT) ) THEN
  !Rotate the base vectors of the system (only if no selection is defined)
  !The H(:,a1) are left untouched
  H2 = H(a1,a2)
  H3 = H(a1,a3)
  H(a1,a2) = H2*rot_matrix(a2,a2) + H3*rot_matrix(a2,a3)
  H(a1,a3) = H2*rot_matrix(a3,a2) + H3*rot_matrix(a3,a3)
  H2 = H(a2,a2)
  H3 = H(a2,a3)
  H(a2,a2) = H2*rot_matrix(a2,a2) + H3*rot_matrix(a2,a3)
  H(a2,a3) = H2*rot_matrix(a3,a2) + H3*rot_matrix(a3,a3)
  H2 = H(a3,a2)
  H3 = H(a3,a3)
  H(a3,a2) = H2*rot_matrix(a2,a2) + H3*rot_matrix(a2,a3)
  H(a3,a3) = H2*rot_matrix(a3,a2) + H3*rot_matrix(a3,a3)
ENDIF
!
!If elastic tensor is set, rotate it
IF( C_tensor(1,1).NE.0.d0 ) THEN
  C_tensor = ROTELAST( C_tensor, rot_matrix )
  CALL ATOMSK_MSG(2099,(/""/),(/0.d0/))
ENDIF
!
!
300 CONTINUE
CALL ATOMSK_MSG(2082,(/''/),(/0.d0/))
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
!
!
END SUBROUTINE ROTATE_XYZ
!
!
END MODULE rotate
