MODULE unit
!
!**********************************************************************************
!*  UNIT                                                                          *
!**********************************************************************************
!* This module reads an array of type (species x y z) and converts                *
!* coordinates, atom velocities, or any property, from one unit to another.       *
!**********************************************************************************
!* (C) October 2010 - Pierre Hirel                                                *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
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
SUBROUTINE UNIT_XYZ(H,P,S,AUXNAMES,AUX,u1,u2,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=128),INTENT(IN):: u1, u2
CHARACTER(LEN=16),DIMENSION(4):: units
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128):: property  !name of property that must be rescaled
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT  !mask for atom list
INTEGER:: i, j
INTEGER:: transform  !0=atom coordinates; 1=atom velocities
INTEGER:: vx, vy, vz !columns in AUX where velocities are stored
REAL(dp):: factor, factor2
REAL(dp),DIMENSION(4):: f
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX !auxiliary properties of atoms/shells
!
!Initialize variables
property=''
units(:) = ''
f(:)=1.d0
factor = 1.d0
transform = 0
vx = 0
vy = 0
vz = 0
!
msg = 'Entering UNIT_XYZ, u1, u2 = '//TRIM(u1)//" , "//TRIM(u2)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
i = SCAN(u1,'/')
j = SCAN(u2,'/')
IF( i==0 .AND. j==0 ) THEN
  !Check if u1 or u2 contain a number
  READ(u1,*,ERR=10,END=10) factor
  !u1 contained a number => u2 must be a property name
  transform=-1
  property=ADJUSTL(u2)
  GOTO 50
  10 CONTINUE
  !u1 did not contain a number => try with u2
  READ(u2,*,ERR=20,END=20) factor
  !u2 contained a number => u1 must be a property name
  transform=-1
  property=ADJUSTL(u1)
  GOTO 50
  !
  20 CONTINUE
  !Neither u1 nor u2 are numbers, they must be both strings containing a unit
  !=> atom coordinates will be transformed
  transform = 0
  units(1) = ADJUSTL(u1)
  units(2) = ADJUSTL(u2)
  !
ELSEIF( i>0 .AND. j>0 ) THEN
  !atom velocitites will be transformed
  transform = 1
  !Check that velocities exist in AUX
  IF( ALLOCATED(AUXNAMES) ) THEN
    DO i=1,SIZE(AUXNAMES)
      IF( TRIM(ADJUSTL(AUXNAMES(i)))=="vx" ) THEN
        vx = i
      ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=="vy" ) THEN
        vy = i
      ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=="vz" ) THEN
        vz = i
      ENDIF
    ENDDO
  ENDIF
  IF( vx==0 .AND. vy==0 .AND. vz==0 ) THEN
    !no velocity exist => exit
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2730,(/'velocity'/),(/0.d0/))
    GOTO 1000
  ELSE
    !save units of distance in units(1) and units(2)
    !save units of time in units(3) and units(4)
    i = SCAN(u1,'/')
    units(1) = u1(1:i-1)
    units(3) = u1(i+1:)
    j = SCAN(u2,'/')
    units(2) = u2(1:j-1)
    units(4) = u2(j+1:)
  ENDIF
ELSE
  nerr = nerr+1
  GOTO 1000
ENDIF
!
50 CONTINUE
WRITE(msg,*) 'Entering UNIT_XYZ: ', units(1), units(2), units(3), units(4)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
IF( transform==-1 ) THEN
  msg = u2
  CALL ATOMSK_MSG(2091,(/property/),(/factor/))
ELSE
  IF( transform==1 ) THEN
    msg = "velocities"
  ELSE
    msg = "coordinates"
  ENDIF
  CALL ATOMSK_MSG(2091,(/msg,units(1),units(2),units(3),units(4)/),(/0.d0/))
  IF( TRIM(ADJUSTL(units(1)))==TRIM(ADJUSTL(units(2))) ) THEN
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2738,(/''/),(/0.d0/))
    GOTO 1000
  ENDIF
ENDIF
!
!
!
100 CONTINUE
IF( transform>=0 ) THEN
  !Define the factor for distance units depending on the two units
  DO i=1,2
    SELECT CASE(units(i))
    CASE('km')
      f(i) = 1.d3
    CASE('m')
      f(i) = 1.d0
    CASE('mm')
      f(i) = 1.d-3
    CASE('µm','micron')
      f(i) = 1.d-6
    CASE('nm')
      f(i) = 1.d-9
    CASE('A','Angstrom','angstrom','Angstroms','angstroms')
      f(i) = 1.d-10
    CASE('B','Bohr','bohr','Bohrs','bohrs')
      f(i) = a_bohr  !defined in constants.f90
    CASE('pm')
      f(i) = 1.d-12
    CASE('fm')
      f(i) = 1.d-15
    CASE DEFAULT
      CALL ATOMSK_MSG(803,(/TRIM(units(i))/),(/0.d0/))
      nerr = nerr+1
      GOTO 1000
    END SELECT
  ENDDO
  factor = f(1) / f(2)
  WRITE(msg,*) 'factor: ', factor
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  IF( transform==1 ) THEN
    !Define factor for time units
    DO i=3,4
      SELECT CASE(units(i))
      CASE('h')
        f(i) = 3600.d0
      CASE('min','mn')
        f(i) = 60.d0
      CASE('s')
        f(i) = 1.d0
      CASE('ms')
        f(i) = 1.d-3
      CASE('µs')
        f(i) = 1.d-6
      CASE('ns')
        f(i) = 1.d-9
      CASE('ps')
        f(i) = 1.d-12
      CASE('fs')
        f(i) = 1.d-15
      CASE DEFAULT
        CALL ATOMSK_MSG(803,(/TRIM(units(i))/),(/0.d0/))
        nerr = nerr+1
        GOTO 1000
      END SELECT
    ENDDO
    factor2 = f(3)/f(4)
    WRITE(msg,*) 'factor2: ', factor2
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
  !
ENDIF
!
!
!
200 CONTINUE
IF( transform==-1 ) THEN
  !Multiply the property u1 by the factor u2
  !
  SELECT CASE(property)
  CASE('x','X')
    !Multiply all X coordinates by the given factor
    DO i=1,SIZE(P,1)
      IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
        P(i,1) = factor*P(i,1)
      ENDIF
    ENDDO
  CASE('y','Y')
    !Multiply all Y coordinates by the given factor
    DO i=1,SIZE(P,1)
      IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
        P(i,2) = factor*P(i,2)
      ENDIF
    ENDDO
  CASE('z','Z')
    !Multiply all Z coordinates by the given factor
    DO i=1,SIZE(P,1)
      IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
        P(i,3) = factor*P(i,3)
      ENDIF
    ENDDO
  CASE DEFAULT
    !Search the property number in AUXNAMES
    j=0
    IF( ALLOCATED(AUXNAMES) .AND. SIZE(AUXNAMES)>=1 ) THEN
      DO i=1,SIZE(AUXNAMES)
        IF( TRIM(ADJUSTL(AUXNAMES(i)))==TRIM(ADJUSTL(property)) ) THEN
          j=i
        ENDIF
      ENDDO
    ENDIF
    !
    IF( j>0 ) THEN
      !Multiply all values of this auxiliary property by the factor
      DO i=1,SIZE(AUX,1)
        IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
          AUX(i,j) = factor*AUX(i,j)
        ENDIF
      ENDDO
    ELSE
      !No property with that name
      nwarn=nwarn+1
      CALL ATOMSK_MSG(2730,(/property/),(/0.d0/))
    ENDIF
  END SELECT
  !
  !
ELSEIF( transform==0 ) THEN
  !Convert atom coordinates
  DO i=1,SIZE(P,1)
    IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
      DO j=1,3
        P(i,j) = P(i,j)*factor
      ENDDO
    ENDIF
  ENDDO
  !
  !Same with shells if they exist
  IF( ALLOCATED(S) .AND. SIZE(S,1)>0 ) THEN
    DO i=1,SIZE(S,1)
      IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
        DO j=1,3
          S(i,j) = S(i,j)*factor
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  !
  !Convert base vectors H
  DO i=1,3
    DO j=1,3
      H(i,j) = H(i,j)*factor
    ENDDO
  ENDDO
  !
  !
ELSEIF( transform==1 ) THEN
  !Convert atom velocities
  IF(vx>0 ) THEN
    DO i=1,SIZE(AUX,1)
      IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
        AUX(i,vx) = AUX(i,vx) * factor / factor2
      ENDIF
    ENDDO
  ENDIF
  IF(vy>0 ) THEN
    DO i=1,SIZE(AUX,1)
      IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
        AUX(i,vy) = AUX(i,vy) * factor / factor2
      ENDIF
    ENDDO
  ENDIF
  IF(vz>0 ) THEN
    DO i=1,SIZE(AUX,1)
      IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
        AUX(i,vz) = AUX(i,vz) * factor / factor2
      ENDIF
    ENDDO
  ENDIF
ENDIF
!
!
IF( transform==-1 ) THEN
  msg = property
ELSEIF( transform==1 ) THEN
  msg = "velocities"
ELSE
  msg = "coordinates"
ENDIF
CALL ATOMSK_MSG(2092,(/msg/),(/0.d0/))
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE UNIT_XYZ
!
!
!
END MODULE unit