MODULE roundoff
!
!**********************************************************************************
!*  ROUNDOFF                                                                      *
!**********************************************************************************
!* This module reads an array of type (species x y z) and converts                *
!* coordinates, atom velocities, or any property, from one unit to another.       *
!**********************************************************************************
!* (C) February 2019 - Pierre Hirel                                               *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 28 Feb. 2019                                     *
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
SUBROUTINE ROUNDOFF_XYZ(P,S,AUXNAMES,AUX,prop,vth,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=128),INTENT(IN):: prop   !the property to round off: may be x, y, z, or any auxiliary property
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT  !mask for atom list
INTEGER:: i, j, j1, j2, jprop
INTEGER:: Nval   !number of values that were rounded off
REAL(dp),INTENT(IN):: vth   !threshold value under which a value is rounded off
REAL(dp),DIMENSION(:,:),POINTER:: Ppoint    !pointer to the array containing the property to round off
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT),TARGET:: P, S !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT),TARGET:: AUX !auxiliary properties of atoms/shells
!
!Initialize variables
jprop=0
Nval=0
!
msg = 'Entering ROUNDOFF_XYZ, property = '//TRIM(prop)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(2146,(/prop/),(/vth/))
!
!
!Detect which property must be rounded off
SELECT CASE(prop)
CASE('x','X')
  Ppoint => P
  jprop = 1
CASE('y','Y')
  Ppoint => P
  jprop = 2
CASE('z','Z')
  Ppoint => P
  jprop = 3
CASE('xyz','XYZ')
  Ppoint => P
  jprop = -3
CASE('aux','AUX')
  Ppoint => AUX
  jprop = -1
CASE DEFAULT
  IF( .NOT.ALLOCATED(AUX) .OR. .NOT.ALLOCATED(AUXNAMES) ) THEN
    !No such property defined: skip
    nwarn = nwarn+1
    GOTO 1000
  ENDIF
  Ppoint=>AUX
  DO i=1,SIZE(AUXNAMES)
    IF( AUXNAMES(i) == prop ) THEN
      jprop = i
    ENDIF
  ENDDO
  !
  IF( jprop==0 ) THEN
    !No such property in AUX: skip
    nwarn=nwarn+1
    GOTO 1000
  ENDIF
END SELECT
!
!
!
100 CONTINUE
!Loop on all values
DO i=1,SIZE(Ppoint,1)
  IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
    !
    IF( jprop>0 ) THEN
      !Only one property must be rounded off
      j1 = jprop
      j2 = jprop
    ELSE
      !All coordinates or auxiliary properties must be rounded off
      j1 = 1
      j2 = SIZE(Ppoint,2)
      !If coordinates must be rounded off, don't round off the atom number
      IF(prop=="XYZ") j2=3
    ENDIF
    !
    DO j=j1,j2
      IF( DABS(Ppoint(i,j)-DBLE(NINT(Ppoint(i,j)))) < DABS(vth) ) THEN
        !This value is close to an integer number: round it off
        Ppoint(i,j) = DBLE(NINT(Ppoint(i,j)))
        Nval = Nval+1
      ENDIF
    ENDDO
    !
  ENDIF
  !
ENDDO
CALL ATOMSK_MSG(2147,(/""/),(/DBLE(Nval)/))
!
!
!
1000 CONTINUE
NULLIFY(Ppoint)
!
!
END SUBROUTINE ROUNDOFF_XYZ
!
!
!
END MODULE roundoff
