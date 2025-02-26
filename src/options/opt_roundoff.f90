MODULE roundoff
!
!**********************************************************************************
!*  ROUNDOFF                                                                      *
!**********************************************************************************
!* This module reads an array of type (species x y z) and rounds off the values   *
!* down to the target accuracy.                                                   *
!**********************************************************************************
!* (C) February 2019 - Pierre Hirel                                               *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 26 Feb. 2025                                     *
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
SUBROUTINE ROUNDOFF_XYZ(H,P,S,AUXNAMES,AUX,prop,accuracy,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=16):: sf  !special format
CHARACTER(LEN=128),INTENT(IN):: prop   !the property to round off: may be x, y, z, or any auxiliary property
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT  !mask for atom list
INTEGER:: i, j, j1, j2, jprop
INTEGER:: Nval   !number of values that were rounded off
INTEGER:: expAcc !exponent of accuracy
REAL(dp),INTENT(IN):: accuracy   !target accuracy: can be 0.01, 1e-6, 10, or any real number
REAL(dp),DIMENSION(3,3),INTENT(INOUT),TARGET:: H !box vectors
REAL(dp),DIMENSION(:,:),POINTER:: Ppoint    !pointer to the array containing the property to round off
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT),TARGET:: P, S !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT),TARGET:: AUX !auxiliary properties of atoms/shells
!
!Initialize variables
jprop=0
Nval=0
!
!Get the exponent of the required accuracy (mantissa is ignored)
!Write accuracy into a string in exponent notation
!NOTE: this will write the number 1.0 in the form 0.1E1, 0.1 as 0.1E0, etc.
!     So the exponent will have to be corrected by -1
WRITE(msg,'(e19.4)') accuracy
!Get the number after the exponent
i=SCAN(msg,'Ee')
IF( i>0 ) THEN
  msg=msg(i+1:)
ENDIF
!Read exponent from the string
READ(msg,*) expAcc
!Remove 1 to correct the exponential notation
expAcc = expAcc-1
!Define specific format using the demanded accuracy
IF( expAcc<0 ) THEN
  !Negative power of 10: the format will be f26.n where n is the number of decimals to keep
  WRITE(sf,'("(f26.", i0, ")")') ABS(expAcc)
ELSE
  !Positive power of 10: the format will not be used (we define it to something arbitrary)
  sf="(f26.0)"
ENDIF
!
msg = 'Entering ROUNDOFF_XYZ, property = '//TRIM(prop)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(2147,(/prop/),(/accuracy/))
!
IF( .NOT.ALLOCATED(P) .OR. SIZE(P,1)<=0 ) THEN
  !No atom in system: can not apply option
  GOTO 1000
ENDIF
!
!
!Detect which property must be rounded off
!jprop is the column number (if jprop<0 then all columns will be rounded off)
SELECT CASE(prop)
CASE('H')
  Ppoint => H
  jprop = -3
CASE('H1','Hx')
  Ppoint => H
  jprop = 1
CASE('H2','Hy')
  Ppoint => H
  jprop = 2
CASE('H3','Hz')
  Ppoint => H
  jprop = 3
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
  IF( .NOT.ALLOCATED(AUX) .OR. .NOT.ALLOCATED(AUXNAMES) ) THEN
    !No aux. property defined: skip
    nwarn = nwarn+1
    CALL ATOMSK_MSG(2729,(/""/),(/0.d0/))
    GOTO 1000
  ELSE
    Ppoint => AUX
    jprop = -1
  ENDIF
CASE DEFAULT
  IF( .NOT.ALLOCATED(AUX) .OR. .NOT.ALLOCATED(AUXNAMES) ) THEN
    !No aux. property defined: skip
    nwarn = nwarn+1
    CALL ATOMSK_MSG(2729,(/""/),(/0.d0/))
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
    !No property with that name in AUX: skip
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2730,(/TRIM(ADJUSTL(prop))/),(/0.d0/))
    GOTO 1000
  ENDIF
END SELECT
!
!
!
100 CONTINUE
!Loop on all values
DO i=1,SIZE(Ppoint,1)
  IF( IS_SELECTED(SELECT,i) ) THEN
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
    IF( expAcc>0 ) THEN
      !Positive power of 10
      DO j=j1,j2
        Ppoint(i,j) = DBLE( NINT( Ppoint(i,j) / 10.d0**(expAcc) ) * 10.d0**(expAcc) )
        Nval = Nval+1
      ENDDO
    ELSE
      !Negative power of 10
      DO j=j1,j2
        !Write number into a string with the specified format (sf)
        WRITE(msg,sf) Ppoint(i,j)
        !Read number from the string
        READ(msg,*) Ppoint(i,j)
        Nval = Nval+1
      ENDDO
    ENDIF
    !
  ENDIF
  !
ENDDO
CALL ATOMSK_MSG(2148,(/""/),(/DBLE(Nval)/))
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
