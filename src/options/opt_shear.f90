MODULE shear
!
!**********************************************************************************
!*  SHEAR                                                                         *
!**********************************************************************************
!* This module reads atomic positions from an array P and                         *
!* applies simple shear strain.                                                   *
!**********************************************************************************
!* (C) November 2010 - Pierre Hirel                                               *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 25 Sept. 2013                                    *
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
SUBROUTINE SHEAR_XYZ(H,P,S,shear_surf,shear_strain,shear_dir)
!
!
IMPLICIT NONE
CHARACTER(LEN=1):: shear_dir  !direction of applied shear strain (X, Y or Z)
CHARACTER(LEN=1):: shear_surf !normal to the surface to be sheared (X, Y or Z)
CHARACTER(LEN=128):: msg
INTEGER:: a1, a2
INTEGER:: i
REAL(dp):: shear_strain  !amplitude of shear in A
REAL(dp):: sstrain       !shear strain in %
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S
!
!
!Initialize variables
i = 0
a1 = 0
a2 = 0
sstrain=0.d0
!
!
CALL ATOMSK_MSG(2083,(/shear_surf,shear_dir/),(/shear_strain/))
!
IF(shear_strain==0.d0) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2735,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Define the axes
IF(shear_dir=='x' .OR. shear_dir=='X') THEN
  a1 = 1
ELSEIF(shear_dir=='y' .OR. shear_dir=='Y') THEN
  a1 = 2
ELSEIF(shear_dir=='z' .OR. shear_dir=='Z') THEN
  a1 = 3
ELSE
  nerr = nerr+1
  CALL ATOMSK_MSG(2800,(/TRIM(shear_surf)/),(/0.d0/))
  GOTO 1000
ENDIF
!
IF(shear_surf=='x' .OR. shear_surf=='X') THEN
  a2 = 1
ELSEIF(shear_surf=='y' .OR. shear_surf=='Y') THEN
  a2 = 2
ELSEIF(shear_surf=='z' .OR. shear_surf=='Z') THEN
  a2 = 3
ELSE
  nerr = nerr+1
  CALL ATOMSK_MSG(2800,(/TRIM(shear_dir)/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
100 CONTINUE
!Transform atom positions to fractional
CALL CART2FRAC(P,H)
IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
  CALL CART2FRAC(S,H)
ENDIF
!
!Apply the shear to the box
H(a2,a1) = H(a2,a1) + shear_strain
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
  150 FORMAT(a7,3(f10.5,2X),a1)
ENDIF
!
!Convert back to cartesian
CALL FRAC2CART(P,H)
IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
  CALL FRAC2CART(S,H)
ENDIF
!
!Compute shear strain (in %)
sstrain = 100.0*shear_strain/H(a2,a2)
!
!
!
300 CONTINUE
CALL ATOMSK_MSG(2084,(/''/),(/sstrain/))
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
END SUBROUTINE SHEAR_XYZ
!
!
END MODULE shear
