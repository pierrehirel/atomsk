MODULE stress
!
!**********************************************************************************
!*  STRESS                                                                        *
!**********************************************************************************
!* This module deforms the box according to the stress required by the user.      *
!* The elastic tensor C_tensor must be provided.                                  *
!**********************************************************************************
!* (C) October 2014 - Pierre Hirel                                                *
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
SUBROUTINE STRESS_XYZ(H,P,S,stress_in,stress,C_tensor)
!
!
IMPLICIT NONE
CHARACTER(LEN=4096),INTENT(IN):: stress_in  !Voigt stress component, or name of file
REAL(dp),DIMENSION(9,9),INTENT(IN):: C_tensor !elastic tensor
CHARACTER(LEN=128):: msg
INTEGER:: i, j
REAL(dp),INTENT(IN):: stress  !value of applied stress (GPa)
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of atoms, shells
REAL(dp),DIMENSION(6):: veps     !deformation (Voigt notation)
REAL(dp),DIMENSION(3,3):: meps     !deformation tensor
REAL(dp),DIMENSION(6):: vstress  !stress (Voigt notation)
REAL(dp),DIMENSION(3,3):: mstress  !stress tensor
REAL(dp),DIMENSION(6,6):: S_tensor !compliance tensor
!
!
!Initialize variables
i = 0
meps(:,:) = 0.d0
veps(:) = 0.d0
mstress(:,:) = 0.d0
vstress(:) = 0.d0
!
!
WRITE(msg,*) 'Entering STRESS_XYZ: '
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(2124,(/stress_in/),(/stress/))
!
!Check that the elastic tensor C_tensor is defined (i.e. non-zero)
IF( .NOT. ANY( C_tensor(1:3,1:3) > 1.d-12 ) ) THEN
  nerr = nerr+1
  CALL ATOMSK_MSG(2816,(/""/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
100 CONTINUE
SELECT CASE(StrDnCase(stress_in))
CASE('x','xx')
  vstress(1) = stress
CASE('y','yy')
  vstress(2) = stress
CASE('z','zz')
  vstress(3) = stress
CASE('yz','zy')
  vstress(4) = stress
CASE('xz','zx')
  vstress(5) = stress
CASE('xy','yx')
  vstress(6) = stress
CASE('p','P')
  !"stress" is actually the value of isostatic pressure
  !=> define corresponding actual stress
  vstress(1) = -1.d0*stress
  vstress(2) = -1.d0*stress
  vstress(3) = -1.d0*stress
CASE DEFAULT
  !It has to be some file name
  !Check that file exists
  CALL CHECKFILE(stress_in,'read')
  !Read stress tensor from file
  OPEN(31,FILE=stress_in,FORM='FORMATTED')
  READ(31,*) msg
  IF(msg=='Voigt' .OR. msg=="voigt") THEN
    !Read the six Voigt components
    READ(31,*) vstress(1), vstress(2), vstress(3)
    READ(31,*) vstress(4), vstress(5), vstress(6)
  ELSEIF(msg=="stress" .OR. msg=="Stress") THEN
    !Read full stress tensor
    DO i=1,3
      READ(31,*) mstress(i,1), mstress(i,2), mstress(i,3)
    ENDDO
    !Convert that into Voigt components
    vstress(1) = mstress(1,1)
    vstress(2) = mstress(2,2)
    vstress(3) = mstress(3,3)
    vstress(4) = mstress(2,3)
    vstress(5) = mstress(1,3)
    vstress(6) = mstress(1,2)
  ELSE
    !unable to read file
    CLOSE(31)
    nerr=nerr+1
    CALL ATOMSK_MSG(800,(/""/),(/0.d0/))
    GOTO 1000
  ENDIF
  CLOSE(31)
END SELECT
!
!
!
200 CONTINUE
i=0
!Invert the elastic tensor to get the compliance tensor
CALL INVMAT(C_tensor(1:6,1:6),S_tensor,i)
!If i is non-zero then the inversion failed
IF(i.NE.0) THEN
  nerr=nerr+1
  CALL ATOMSK_MSG(2815,(/"C_tensor"/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Compute the Voigt components of deformation
!That's Hooke's law:  [eps] = [S_tensor] * [stress]
veps(:) = 0.d0
DO i=1,6
  DO j=1,6
    veps(i) = veps(i) + S_tensor(i,j) * vstress(j)
  ENDDO
ENDDO
!
!Translate that into the full deformation tensor
!NOTE: we make the deformation tensor a lower triangular matrix,
!     e.g. instead of having meps(2,3) = meps(3,2) = veps(4)/2
!     we use meps(3,2) = veps(4) and meps(2,3)=0.
meps(:,:) = 0.d0
meps(1,1) = 1.d0 + veps(1)
meps(2,2) = 1.d0 + veps(2)
meps(3,3) = 1.d0 + veps(3)
meps(3,2) = veps(4)
meps(3,1) = veps(5)
meps(2,1) = veps(6)
!
IF(verbosity==4) THEN
  msg = 'Deformation tensor:'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,150) '     | ', meps(1,1), meps(1,2), meps(1,3), '|'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,150) ' e = | ', meps(2,1), meps(2,2), meps(2,3), '|'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,150) '     | ', meps(3,1), meps(3,2), meps(3,3), '|'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!
!
!
300 CONTINUE
!Transform atom positions to fractional
CALL CART2FRAC(P,H)
IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
  CALL CART2FRAC(S,H)
ENDIF
!
!
!Apply the strain to the box
H(:,:) = MATMUL ( H(:,:) , meps(:,:) )
!
IF(verbosity==4) THEN
  msg = 'New base vectors:'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,150) '     | ', H(1,1), H(1,2), H(1,3), '|'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,150) ' H = | ', H(2,1), H(2,2), H(2,3), '|'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,150) '     | ', H(3,1), H(3,2), H(3,3), '|'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  150 FORMAT(a17,3(f10.5,2X),a1)
ENDIF
!
!
!Convert back to cartesian
CALL FRAC2CART(P,H)
IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
  CALL FRAC2CART(S,H)
ENDIF
!
!
!
400 CONTINUE
CALL ATOMSK_MSG(2060,(/''/),(/0.d0/))
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
END SUBROUTINE STRESS_XYZ
!
!
END MODULE stress
