MODULE deform
!
!**********************************************************************************
!*  DEFORM                                                                        *
!**********************************************************************************
!* This module reads cartesian coordinates from an array P and applies            *
!* the specified strain. If a Poisson's ratio is specified, directions normal     *
!* to given axis is deformed accordingly.                                         *
!**********************************************************************************
!* (C) March 2010 - Pierre Hirel                                                  *
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
SUBROUTINE DEFORM_XYZ(H,P,S,def_dir,def_strain,def_poisson,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=2),INTENT(IN):: def_dir  !strain component: x, y, z, xy, xz...
CHARACTER(LEN=128):: msg
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: a1, a2, a3
INTEGER:: i
REAL(dp),INTENT(IN):: def_strain, def_poisson  !applied strain and Poisson's ratio
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of atoms, shells
REAL(dp),DIMENSION(3,3):: meps  !applied strain along X, Y, Z
!
!
!Initialize variables
i = 0
a1 = 0
a2 = 0
a3 = 0
meps(:,:) = 0.d0
!
!
CALL ATOMSK_MSG(2058,(/def_dir/),(/def_strain,def_poisson/))
!
IF(def_strain==0.d0) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2724,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Display a warning if Poisson ratio is not between -1 and +0.5
IF(def_poisson<-1.d0 .OR. def_poisson>0.5d0) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2741,(/''/),(/0.d0/))
ENDIF
!
!
!Define the axes
SELECT CASE(StrDnCase(def_dir))
CASE('x',"xx")
  a1 = 1
  a2 = 2
  a3 = 3
CASE('y',"yy")
  a1 = 2
  a2 = 3
  a3 = 1
CASE('z',"zz")
  a1 = 3
  a2 = 1
  a3 = 2
CASE("xy","xz","yx","yz","zx","zy")
  DO i=1,3
    meps(i,i) = 1.d0
  ENDDO
  CONTINUE
CASE DEFAULT
  nerr = nerr+1
  CALL ATOMSK_MSG(2800,(/def_dir/),(/0.d0/))
  GOTO 1000
END SELECT
!
!Isotropic case:
!calculate the deformation matrix
SELECT CASE(StrDnCase(def_dir))
CASE('x','y','z')
  meps(a1,a1) = 1.d0+def_strain
  meps(a2,a2) = 1.d0-def_poisson*def_strain
  meps(a3,a3) = 1.d0-def_poisson*def_strain
CASE('xy')
  meps(2,1) = def_strain
CASE('xz')
  meps(3,1) = def_strain
CASE('yz')
  meps(3,2) = def_strain
CASE('yx')
  meps(1,2) = def_strain
CASE('zx')
  meps(1,3) = def_strain
CASE('zy')
  meps(2,3) = def_strain
END SELECT
!
!
!
100 CONTINUE
IF( ALLOCATED(P) .AND. SIZE(P,1)>0 ) THEN
  !Transform atom positions to fractional
  CALL CART2FRAC(P,H)
  IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
    CALL CART2FRAC(S,H)
  ENDIF
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
IF( ALLOCATED(P) .AND. SIZE(P,1)>0 ) THEN
  !Convert back to cartesian
  CALL FRAC2CART(P,H)
  IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
    CALL FRAC2CART(S,H)
  ENDIF
ENDIF
!
!
!
300 CONTINUE
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
END SUBROUTINE DEFORM_XYZ
!
!
END MODULE deform
