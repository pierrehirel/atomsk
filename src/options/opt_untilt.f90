MODULE untilt
!
!**********************************************************************************
!*  UNTILT                                                                        *
!**********************************************************************************
!* This module removes non-diagonl components of the cell vectors matrix,         *
!* i.e. it removes the "tilt" of cell vectors. It is possible to remove           *
!* only one given tilt.                                                           *
!**********************************************************************************
!* (C) September 2023 - Pierre Hirel                                              *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 19 Sept. 2023                                    *
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
SUBROUTINE UNTILT_XYZ(H,P,S,untilt_dir)
!
!
IMPLICIT NONE
CHARACTER(LEN=2),INTENT(IN):: untilt_dir  !strain component: xy, xz, yz, yx, zx, zy (or empty)
CHARACTER(LEN=128):: msg
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: i
REAL(dp),DIMENSION(3,3):: Hcp   !copy of base vectors of the supercell
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of atoms, shells
!
!
!Initialize variables
i = 0
!
!
CALL ATOMSK_MSG(2155,(/untilt_dir/),(/0.d0/))
!
!Copy non-diagonal elements of H and verify if they are non-zero
Hcp(:,:) = H(:,:)
DO i=1,3
  Hcp(i,i) = 0.d0
ENDDO
IF( .NOT.ANY(Hcp(:,:).NE.0.d0) ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2767,(/''/),(/0.d0/))
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
!Remove tilt component
SELECT CASE(StrDnCase(untilt_dir))
CASE('xy')
  H(2,1) = 0.d0
CASE('xz')
  H(3,1) = 0.d0
CASE('yz')
  H(3,2) = 0.d0
CASE('yx')
  H(1,2) = 0.d0
CASE('zx')
  H(1,3) = 0.d0
CASE('zy')
  H(2,3) = 0.d0
CASE DEFAULT
  Hcp(:,:) = H(:,:)
  H(:,:) = 0.d0
  DO i=1,3
    H(i,i)= Hcp(i,i)
  ENDDO
END SELECT
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
300 CONTINUE
CALL ATOMSK_MSG(2156,(/''/),(/0.d0/))
GOTO 1000
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE UNTILT_XYZ
!
!
END MODULE untilt
