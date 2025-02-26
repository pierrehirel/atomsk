MODULE carttofrac
!
!**********************************************************************************
!*  CART2FRAC                                                                     *
!**********************************************************************************
!* This module converts atom coordinates from cartesian to reduced coordinates.   *
!**********************************************************************************
!* (C) Nov. 2010 - Pierre Hirel                                                   *
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
!
USE comv
USE constants
USE messages
USE files
USE subroutines
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE CART2FRAC_XYZ(H,P,S)
!
CHARACTER(LEN=128):: msg
LOGICAL:: isreduced
REAL(dp), DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S
!
!
msg = 'entering CART2FRAC_XYZ'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
!
CALL ATOMSK_MSG(2052,(/msg/),(/0.d0/))
!
IF( .NOT.ALLOCATED(P) .OR. SIZE(P,1)<=0 ) THEN
  !No atom in system: can not apply option
  GOTO 1000
ENDIF
!
!
100 CONTINUE
!Check if coordinates are already reduced or not
CALL FIND_IF_REDUCED(H,P,isreduced)
!If it is the case then we skip
IF(isreduced) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2721,(/msg/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
200 CONTINUE
!Convert to reduced coordinates
CALL CART2FRAC(P,H)
IF( ALLOCATED(S) ) THEN
  CALL CART2FRAC(S,H)
ENDIF
!
!
!
300 CONTINUE
CALL ATOMSK_MSG(2053,(/msg/),(/0.d0/))
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE CART2FRAC_XYZ
!
!
END MODULE carttofrac
