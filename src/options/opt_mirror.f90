MODULE mirror
!
!**********************************************************************************
!*  MIRROR                                                                        *
!**********************************************************************************
!* This module reads atomic coordinates from an array P, and applies a            *
!* mirror transformation. The mirror plane is defined by its direction            *
!* (which can be X, Y, Z or any crystallographic direction) and its distance      *
!* from the origin.                                                               *
!**********************************************************************************
!* (C) May 2014 - Pierre Hirel                                                    *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 01 March 2017                                    *
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
SUBROUTINE MIRROR_XYZ(P,S,AUXNAMES,AUX,mirror_dir,mirror_d,ORIENT,SELECT)
!
IMPLICIT NONE
CHARACTER(LEN=16),INTENT(IN):: mirror_dir   !x, y, z, or crystallographic direction
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: AUXNAMES !names of auxiliary properties
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT   !mask for atom list
INTEGER:: a1
INTEGER:: i, j
INTEGER:: fx, fy, fz, vx, vy, vz  !columns for forces and velocities in AUX
REAL(dp),INTENT(IN):: mirror_d
REAL(dp):: tempreal
REAL(dp):: V1, V2, V3  !vector components
REAL(dp),DIMENSION(3):: vector    !shortest vector between mirror plane and an atom
REAL(dp),DIMENSION(1,3):: Vplane  !vector defining the plane
REAL(dp),DIMENSION(3,3),INTENT(IN):: ORIENT !current crystallographic orientation of the system
REAL(dp),DIMENSION(3,3):: ORIENTN      !normalized ORIENT
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S  !positions of cores, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX   !auxiliary properties of atoms
!
!
!
!Initialize variables
i = 0
fx=0
fy=0
fz=0
vx=0
vy=0
vz=0
!
!
CALL ATOMSK_MSG(2120,(/mirror_dir/),(/mirror_d/))
!
!
!Detect if velocities and/or forces are defined
IF( ALLOCATED(AUXNAMES) ) THEN
  DO i=1,SIZE(AUXNAMES)
    IF( AUXNAMES(i)=="fx" ) fx=i
    IF( AUXNAMES(i)=="fy" ) fy=i
    IF( AUXNAMES(i)=="fz" ) fz=i
    IF( AUXNAMES(i)=="vx" ) vx=i
    IF( AUXNAMES(i)=="vy" ) vy=i
    IF( AUXNAMES(i)=="vz" ) vz=i
  ENDDO
  !
  IF( fx<=0 .OR. fy<=0 .OR. fz<=0 ) THEN
    fx=0
    fy=0
    fz=0
  ENDIF
  IF( vx<=0 .OR. vy<=0 .OR. vz<=0 ) THEN
    vx=0
    vy=0
    vz=0
  ENDIF
ENDIF
!
!
!
100 CONTINUE
!
SELECT CASE(mirror_dir)
CASE("x","X","y","Y","z","Z")
  !mirror_dir is a cartesian direction
  !Define the axes
  IF(mirror_dir=='x' .OR. mirror_dir=='X') THEN
    a1 = 1
  ELSEIF(mirror_dir=='y' .OR. mirror_dir=='Y') THEN
    a1 = 2
  ELSEIF(mirror_dir=='z' .OR. mirror_dir=='Z') THEN
    a1 = 3
  ENDIF
  WRITE(msg,*) 'a1 = ', a1
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !
  DO i=1,SIZE(P,1)
    IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
      !Apply mirror to atom position
      tempreal = P(i,a1) - mirror_d
      P(i,a1) = P(i,a1) - 2.d0*tempreal
      !
      !Same with shells if any
      IF( ALLOCATED(S) ) THEN
        tempreal = S(i,a1) - mirror_d
        S(i,a1) = S(i,a1) - 2.d0*tempreal
      ENDIF
      !
      !Mirror forces and velocities if they are defined
      IF( fx>0 ) THEN
        IF( a1==1 ) THEN
          AUX(i,fx) = -AUX(i,fx)
        ELSEIF( a1==2 ) THEN
          AUX(i,fy) = -AUX(i,fy)
        ELSEIF( a1==3 ) THEN
          AUX(i,fz) = -AUX(i,fz)
        ENDIF
      ENDIF
      !
      IF( vx>0 ) THEN
        IF( a1==1 ) THEN
          AUX(i,vx) = -AUX(i,vx)
        ELSEIF( a1==2 ) THEN
          AUX(i,vy) = -AUX(i,vy)
        ELSEIF( a1==3 ) THEN
          AUX(i,vz) = -AUX(i,vz)
        ENDIF
      ENDIF
      !
    ENDIF
  ENDDO
  !
CASE DEFAULT
  !mirror_dir should contain a crystallographic direction
  !convert it to a vector and save it in Vplane(1,:)
  CALL INDEX_MILLER(mirror_dir,Vplane(1,:),j)
  IF(j>0) GOTO 800
  !
  !If the system has a defined crystallographic orientation ORIENT,
  !then Vplane(1,:) is defined in that basis
  !=> rotate Vplane(1,:) to express it in cartesian basis
  IF( ANY( NINT(ORIENT(:,:)).NE.0 ) ) THEN
    DO i=1,3
      ORIENTN(i,:) = ORIENT(i,:) / VECLENGTH(ORIENT(i,:))
    ENDDO
    V1 = Vplane(1,1)
    V2 = Vplane(1,2)
    V3 = Vplane(1,3)
    Vplane(1,1) = ORIENTN(1,1)*V1 + ORIENTN(1,2)*V2 + ORIENTN(1,3)*V3
    Vplane(1,2) = ORIENTN(2,1)*V1 + ORIENTN(2,2)*V2 + ORIENTN(2,3)*V3
    Vplane(1,3) = ORIENTN(3,1)*V1 + ORIENTN(3,2)*V2 + ORIENTN(3,3)*V3
  ENDIF
  !Normalize Vplane
  Vplane(1,:) = Vplane(1,:)/VECLENGTH(Vplane(1,:))
  WRITE(msg,'(a8,3f12.3)') 'Vplane: ', Vplane(1,:)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !
  DO i=1,SIZE(P,1)
    IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
      !compute shortest vector between atom #i and the mirror plane
      vector = P(i,1:3) - mirror_d*Vplane(1,:)/VECLENGTH(Vplane(1,:))
      P(i,1:3) = P(i,1:3) - 2.d0*vector(:)
      !
      !Same with shells if any
      IF( ALLOCATED(S) ) THEN
        vector = S(i,1:3) - mirror_d*Vplane(1,:)/VECLENGTH(Vplane(1,:))
        S(i,1:3) = S(i,1:3) - 2.d0*vector(:)
      ENDIF
      !
      !Mirror forces and velocities if they are defined
      IF( fx>0 ) THEN
        !Compute force component normal to mirror plane
        V1 = AUX(i,fx)
        V2 = AUX(i,fy)
        V3 = AUX(i,fz)
        vector = DOT_PRODUCT( (/V1,V2,V3/) , Vplane(1,:)/VECLENGTH(Vplane(1,:)) )
        AUX(i,fx) = AUX(i,fx) - 2.d0*vector(1)
        AUX(i,fy) = AUX(i,fy) - 2.d0*vector(2)
        AUX(i,fz) = AUX(i,fz) - 2.d0*vector(3)
      ENDIF
      !
      IF( vx>0 ) THEN
        !Compute velocity component normal to mirror plane
        V1 = AUX(i,vx)
        V2 = AUX(i,vy)
        V3 = AUX(i,vz)
        vector = DOT_PRODUCT( (/V1,V2,V3/) , Vplane(1,:)/VECLENGTH(Vplane(1,:)) )
        AUX(i,vx) = AUX(i,vx) - 2.d0*vector(1)
        AUX(i,vy) = AUX(i,vy) - 2.d0*vector(2)
        AUX(i,vz) = AUX(i,vz) - 2.d0*vector(3)
      ENDIF
      !
    ENDIF
  ENDDO
  !
END SELECT
!
!
CALL ATOMSK_MSG(2121,(/''/),(/0.d0/))
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
!
!
END SUBROUTINE MIRROR_XYZ
!
END MODULE mirror
