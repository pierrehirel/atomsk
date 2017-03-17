MODULE spacegroup
!
!**********************************************************************************
!*  SPACEGROUP                                                                    *
!**********************************************************************************
!* This module reads cartesian coordinates from an array and                      *
!* applies the symmetry operations of the given space group.                      *
!**********************************************************************************
!* (C) January 2016 - Pierre Hirel                                                *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 02 May 2016                                      *
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
USE subroutines
USE symops
!
!
CONTAINS
!
!
SUBROUTINE SPACEGROUP_XYZ(H,P,S,AUXNAMES,AUX,sgroup)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: sgroup  !Hermann-Mauguin symbol or number of space group
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: AUXNAMES !names of auxiliary properties
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX  !auxiliary properties of atoms/shells
!
! Initializations
!
!
WRITE(msg,*) 'Entering SPACEGROUP_XYZ: '//TRIM(ADJUSTL(sgroup))
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(2131,(/TRIM(sgroup)/),(/0.d0/))
!
!
!
100 CONTINUE
!Apply symmetry operations
!(cf. /include/symops.f90)
CALL SG_APPLY_SYMOPS(sgroup,H,P,S,AUXNAMES,AUX)
!
CALL ATOMSK_MSG(2133,(/""/),(/DBLE(SIZE(P,1))/))
!
!
!
1000 CONTINUE
!
!
!
END SUBROUTINE SPACEGROUP_XYZ
!
!
END MODULE spacegroup