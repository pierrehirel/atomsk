MODULE rotate
!
!**********************************************************************************
!*  ROTATE                                                                        *
!**********************************************************************************
!* This module reads atomic positions from an array P and rotates                 *
!* the system by a certain angle around the cartesian X, Y or Z axis.             *
!**********************************************************************************
!* (C) May 2010 - Pierre Hirel                                                    *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 20 Sept. 2023                                    *
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
USE crystallography
USE elasticity
USE messages
USE files
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE ROTATE_XYZ(H,P,S,AUXNAMES,AUX,com,rot_axis,rot_angle,ORIENT,SELECT,C_tensor)
!
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: rot_axis  ! Cartesian x, y or z axis, Miller vector, or 3 real numbers
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
INTEGER:: a1, a2, a3
INTEGER,INTENT(IN):: com  !=1 if rotation around the system center of mass, 0 otherwise
INTEGER:: i
INTEGER,DIMENSION(3):: Fxyz, Vxyz !columns of AUX containing forces, velocities
LOGICAL:: velocities, forces
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
REAL(dp):: H1, H2, H3
REAL(dp):: rot_angle      !angle in degrees
REAL(dp):: smass          !mass of an atom
REAL(dp):: totmass        !mass of all (or selected) atoms
REAL(dp):: V1, V2, V3
REAL(dp):: u, v, w, x, z1, z2
REAL(dp),DIMENSION(3):: MILLER       !Miller indices
REAL(dp),DIMENSION(3):: Vcom !position of center of rotation
REAL(dp),DIMENSION(3):: Vrot !vector around which the rotation is made
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3),INTENT(IN):: ORIENT !current crystallographic orientation of the system
REAL(dp),DIMENSION(3,3):: ORIENTN      !normalized ORIENT
REAL(dp),DIMENSION(3,3):: rot_matrix
REAL(dp),DIMENSION(9,9),INTENT(INOUT):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of cores, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX  !auxiliary properties of atoms/shells
!
!Initialize variables
forces=.FALSE.
velocities=.FALSE.
a1=0
a2=0
a3=0
i = 0
Fxyz(:)=0
Vxyz(:)=0
H2 = 0.d0
H3 = 0.d0
MILLER(:) = 0.d0
Vcom(:) = 0.d0
Vrot(:) = 0.d0
rot_matrix(:,:) = 0.d0
DO i=1,3
  rot_matrix(i,i) = 1.d0
ENDDO
!
msg = 'Entering ROTATE_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(2081,(/rot_axis/),(/rot_angle,DBLE(com)/))
!
!If angle is more than 360°, reduce it
DO WHILE(DABS(rot_angle)>=360.d0)
  IF(rot_angle>=360.d0)  rot_angle = rot_angle-360.d0
  IF(rot_angle<=-360.d0) rot_angle = rot_angle+360.d0
ENDDO
!
IF( DABS(rot_angle)<=1.d-12) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2734,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!If rotation around center of mass, compute position of COM
IF( com .NE. 0 ) THEN
  totmass = 0.d0
  Vcom(:) = 0.d0
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,species,smass) REDUCTION(+:Vcom,totmass)
  DO i=1,SIZE(P,1)
    IF( .NOT.(ALLOCATED(SELECT)) .OR. SELECT(i) ) THEN
      CALL ATOMSPECIES(P(i,4),species)
      CALL ATOMMASS(species,smass)
      Vcom(:) = Vcom(:) + smass*P(i,1:3)
      totmass = totmass + smass
    ENDIF
  ENDDO
  !$OMP END PARALLEL DO
  Vcom(:) = Vcom(:) / totmass
ENDIF
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
  !Directions will be simply a1=X, a2=Y, a3=Z
  a1 = 1
  a2 = 2
  a3 = 3
  !
  IF( SCAN(rot_axis,'[]_')>0 ) THEN
    !It may be a Miller direction
    CALL INDEX_MILLER(rot_axis,MILLER,i)
    IF( i>0 ) THEN
      !Miller direction may be given for hexagonal
      CALL INDEX_MILLER_HCP(rot_axis,MILLER,i)
      IF(i==0) THEN
        !Convert [hkil] notation into [uvw]
        CALL HKIL2UVW(MILLER(1),MILLER(2),0.d0,MILLER(3),u,v,w)
        !Set indices in MILLER
        MILLER(:) = u*H(1,:) + v*H(2,:) + w*H(3,:)
      ENDIF
    ENDIF
    IF(i==0) THEN
      !It was a Miller direction
      !Check that it is not [000]
      IF( VECLENGTH(MILLER)<1.d-12 ) THEN
        CALL ATOMSK_MSG(814,(/""/),(/0.d0/))
        nerr=nerr+1
        GOTO 1000
      ENDIF
      !
      !If the system has a defined crystallographic orientation ORIENT,
      !then Vrot(1,:) is defined in that basis
      !=> rotate Vrot(1,:) to express it in cartesian basis
      IF( ANY( NINT(ORIENT(:,:)).NE.0 ) ) THEN
        DO i=1,3
          ORIENTN(i,:) = ORIENT(i,:) / VECLENGTH(ORIENT(i,:))
        ENDDO
        V1 = MILLER(1)
        V2 = MILLER(2)
        V3 = MILLER(3)
        MILLER(1) = ORIENTN(1,1)*V1 + ORIENTN(1,2)*V2 + ORIENTN(1,3)*V3
        MILLER(2) = ORIENTN(2,1)*V1 + ORIENTN(2,2)*V2 + ORIENTN(2,3)*V3
        MILLER(3) = ORIENTN(3,1)*V1 + ORIENTN(3,2)*V2 + ORIENTN(3,3)*V3
      ENDIF
      !
      Vrot(:) = MILLER(:)
    ELSE
      !It was not a Miller vector: impossible to understand, abort
      GOTO 801
    ENDIF
    !
  ELSE
    !Last possibility: rot_axis should contain 3 real numbers
    READ(rot_axis,*,ERR=801,END=801) V1, V2, V3
    Vrot(:) = (/ V1 , V2 , V3 /)
  ENDIF
ENDIF
!
WRITE(msg,'(a6,3f12.3)') 'Vrot: ', Vrot(:)
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
IF( VECLENGTH(Vrot(:))<1.d-6 ) THEN
  !convert the angle into radians
  rot_angle = DEG2RAD(rot_angle)
  !set the rotation matrix
  rot_matrix(:,:) = 0.d0
  rot_matrix(a1,a1) = 1.d0
  rot_matrix(a2,a2) = DCOS(rot_angle)
  rot_matrix(a2,a3) = -DSIN(rot_angle)
  rot_matrix(a3,a2) = DSIN(rot_angle)
  rot_matrix(a3,a3) = DCOS(rot_angle)
  !
ELSE
  !Direction was Miller index or a vector
  !Normalize Vrot
  Vrot(:) = Vrot(:)/VECLENGTH(Vrot(:))
  !
  !Construct rotation matrix
  rot_matrix = ROTMAT_AXIS(Vrot,rot_angle)
ENDIF
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
!Rotate only atoms that are selected in SELECT
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,H1,H2,H3)
DO i=1,SIZE(P,1)
  IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
    H1 = P(i,a1) - Vcom(a1)
    H2 = P(i,a2) - Vcom(a2)
    H3 = P(i,a3) - Vcom(a3)
    P(i,a1) = Vcom(a1) + H1*rot_matrix(a1,a1) + H2*rot_matrix(a1,a2) + H3*rot_matrix(a1,a3)
    P(i,a2) = Vcom(a2) + H1*rot_matrix(a2,a1) + H2*rot_matrix(a2,a2) + H3*rot_matrix(a2,a3)
    P(i,a3) = Vcom(a3) + H1*rot_matrix(a3,a1) + H2*rot_matrix(a3,a2) + H3*rot_matrix(a3,a3)
    !
    !Same with shell if they exist
    IF( ALLOCATED(S) .AND. SIZE(S,1)>0 ) THEN
      H1 = S(i,a1) - Vcom(a1)
      H2 = S(i,a2) - Vcom(a2)
      H3 = S(i,a3) - Vcom(a3)
      S(i,a1) = Vcom(a1) + H1*rot_matrix(a1,a1) + H2*rot_matrix(a1,a2) + H3*rot_matrix(a1,a3)
      S(i,a2) = Vcom(a2) + H1*rot_matrix(a2,a1) + H2*rot_matrix(a2,a2) + H3*rot_matrix(a2,a3)
      S(i,a3) = Vcom(a3) + H1*rot_matrix(a3,a1) + H2*rot_matrix(a3,a2) + H3*rot_matrix(a3,a3)
    ENDIF
    !
    !Same with forces if they exist
    IF( forces ) THEN
      H1 = AUX(i,Fxyz(a1))
      H2 = AUX(i,Fxyz(a2))
      H3 = AUX(i,Fxyz(a3))
      AUX(i,Fxyz(a1)) = H1*rot_matrix(a1,a1) + H2*rot_matrix(a1,a2) + H3*rot_matrix(a1,a3)
      AUX(i,Fxyz(a2)) = H1*rot_matrix(a2,a1) + H2*rot_matrix(a2,a2) + H3*rot_matrix(a2,a3)
      AUX(i,Fxyz(a3)) = H1*rot_matrix(a3,a1) + H2*rot_matrix(a3,a2) + H3*rot_matrix(a3,a3)
    ENDIF
    !
    !Same with velocities if they exist
    IF( velocities ) THEN
      H1 = AUX(i,Vxyz(a1))
      H2 = AUX(i,Vxyz(a2))
      H3 = AUX(i,Vxyz(a3))
      AUX(i,Vxyz(a1)) = H1*rot_matrix(a1,a1) + H2*rot_matrix(a1,a2) + H3*rot_matrix(a1,a3)
      AUX(i,Vxyz(a2)) = H1*rot_matrix(a2,a1) + H2*rot_matrix(a2,a2) + H3*rot_matrix(a2,a3)
      AUX(i,Vxyz(a3)) = H1*rot_matrix(a3,a1) + H2*rot_matrix(a3,a2) + H3*rot_matrix(a3,a3)
    ENDIF
    !
  ENDIF
ENDDO
!$OMP END PARALLEL DO
!
IF( .NOT.ALLOCATED(SELECT) .AND. com==0 ) THEN
  !Rotate the base vectors of the system (only if no selection is defined)
  !The H(:,a1) are left untouched
  H1 = H(a1,a1)
  H2 = H(a1,a2)
  H3 = H(a1,a3)
  H(a1,a1) = H1*rot_matrix(a1,a1) + H2*rot_matrix(a1,a2) + H3*rot_matrix(a1,a3)
  H(a1,a2) = H1*rot_matrix(a2,a1) + H2*rot_matrix(a2,a2) + H3*rot_matrix(a2,a3)
  H(a1,a3) = H1*rot_matrix(a3,a1) + H2*rot_matrix(a3,a2) + H3*rot_matrix(a3,a3)
  H1 = H(a2,a1)
  H2 = H(a2,a2)
  H3 = H(a2,a3)
  H(a2,a1) = H1*rot_matrix(a1,a1) + H2*rot_matrix(a1,a2) + H3*rot_matrix(a1,a3)
  H(a2,a2) = H1*rot_matrix(a2,a1) + H2*rot_matrix(a2,a2) + H3*rot_matrix(a2,a3)
  H(a2,a3) = H1*rot_matrix(a3,a1) + H2*rot_matrix(a3,a2) + H3*rot_matrix(a3,a3)
  H1 = H(a3,a1)
  H2 = H(a3,a2)
  H3 = H(a3,a3)
  H(a3,a1) = H1*rot_matrix(a1,a1) + H2*rot_matrix(a1,a2) + H3*rot_matrix(a1,a3)
  H(a3,a2) = H1*rot_matrix(a2,a1) + H2*rot_matrix(a2,a2) + H3*rot_matrix(a2,a3)
  H(a3,a3) = H1*rot_matrix(a3,a1) + H2*rot_matrix(a3,a2) + H3*rot_matrix(a3,a3)
ENDIF
!
!If elastic tensor is set, rotate it
IF( DABS(C_tensor(1,1))>1.d-6 ) THEN
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
GOTO 1000
!
801 CONTINUE
!Nothing worked: impossible to understand this direction
CALL ATOMSK_MSG(2800,(/rot_axis/),(/0.d0/))
nerr = nerr+1
GOTO 1000
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
