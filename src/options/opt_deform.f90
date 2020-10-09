MODULE deform
!
!**********************************************************************************
!*  DEFORM                                                                        *
!**********************************************************************************
!* This module reads cartesian coordinates from an array P and applies            *
!* normal stress according to specified strain and Poisson's ratio.               *
!**********************************************************************************
!* (C) March 2010 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 09 Oct. 2020                                     *
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
CHARACTER(LEN=1),INTENT(IN):: def_dir  !direction of applied strain (X, Y or Z)
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
CALL ATOMSK_MSG(2058,(/def_dir/),(/def_strain/))
IF(def_strain==0.d0) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2724,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
CALL ATOMSK_MSG(2059,(/''/),(/def_poisson/))
!Display a warning if Poisson ratio is not between -1 and +0.5
IF(def_poisson<-1.d0 .OR. def_poisson>0.5d0) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2741,(/''/),(/0.d0/))
ENDIF
!
!
!Define the axes
IF(def_dir=='x' .OR. def_dir=='X') THEN
  a1 = 1
  a2 = 2
  a3 = 3
ELSEIF(def_dir=='y' .OR. def_dir=='Y') THEN
  a1 = 2
  a2 = 3
  a3 = 1
ELSEIF(def_dir=='z' .OR. def_dir=='Z') THEN
  a1 = 3
  a2 = 1
  a3 = 2
ELSE
  nerr = nerr+1
  CALL ATOMSK_MSG(2800,(/def_dir/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Isotropic case:
!calculate the deformation matrix
!def_strain, def_poisson
meps(:,:) = 0.d0
meps(a1,a1) = 1.d0+def_strain
meps(a2,a2) = 1.d0-def_poisson*def_strain
meps(a3,a3) = 1.d0-def_poisson*def_strain
!
!
!
100 CONTINUE
!Apply the strain to the box
IF(.NOT.ALLOCATED(SELECT)) THEN
  H(1,1) = H(1,1)*meps(1,1)
  H(2,2) = H(2,2)*meps(2,2)
  H(3,3) = H(3,3)*meps(3,3)
ENDIF
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
!
200 CONTINUE
!Adjust positions of atoms
DO i=1,SIZE(P(:,1))
  IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
    P(i,1) = P(i,1)*meps(1,1)
    P(i,2) = P(i,2)*meps(2,2)
    P(i,3) = P(i,3)*meps(3,3)
  ENDIF
ENDDO
!
!Same with shells
IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
  DO i=1,SIZE(S(:,1))
    IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
      S(i,1) = S(i,1)*meps(1,1)
      S(i,2) = S(i,2)*meps(2,2)
      S(i,3) = S(i,3)*meps(3,3)
    ENDIF
  ENDDO
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
