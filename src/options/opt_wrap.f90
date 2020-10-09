MODULE wrap
!
!**********************************************************************************
!*  WRAP                                                                          *
!**********************************************************************************
!* This module reads an array of type (species x y z) and wraps atoms             *
!* back into the box defined by the base vectors H.                               *
!* NOTE: if shells are present, they are translated by the same vector            *
!* as their core, so that core-shell vector remain the same.                      *
!**********************************************************************************
!* (C) October 2010 - Pierre Hirel                                                *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 08 Oct. 2020                                     *
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
SUBROUTINE WRAP_XYZ(H,P,S,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
LOGICAL:: doshells !treat shells?
LOGICAL:: isreduced
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT  !mask for atom list
INTEGER:: i, j
INTEGER:: nloop   !number of iterations
INTEGER:: NPwrap  !number of atoms wrapped
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of atoms, shells
!
!Initialize variables
doshells=.FALSE.
NPwrap = 0
!
WRITE(msg,*) 'Entering WRAP_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(2093,(/TRIM(msg)/),(/0.d0/))
!
!Find out if shells are defined
IF( ALLOCATED(S) .AND. SIZE(S,1)>0 ) THEN
  doshells=.TRUE.
ENDIF
!
!
100 CONTINUE
!Determine if coordinates are already reduced or not
CALL FIND_IF_REDUCED(H,P,isreduced)
!
!If not, then transform positions into reduced coordinates
IF(.NOT.isreduced) THEN
  CALL CART2FRAC(P,H)
  IF(doshells) THEN
    CALL CART2FRAC(S,H)
  ENDIF
ENDIF
!
!Second, wrap all coordinates greater than 1
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,nloop) REDUCTION(+:NPwrap)
DO i=1,SIZE(P,1) !loop on all atoms
  IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
    DO j=1,3  !loop on xyz
      nloop=0
      IF(P(i,j)>=1.d0 .OR. P(i,j)<0.d0) NPwrap = NPwrap+1
      DO WHILE( P(i,j)>=1.d0 .AND. nloop<1000)
        P(i,j) = P(i,j)-1.d0
        IF(doshells) THEN
          S(i,j) = S(i,j)-1.d0
        ENDIF
        nloop=nloop+1
      ENDDO
      nloop=0
      DO WHILE( P(i,j)<0.d0  .AND. nloop<1000)
        P(i,j) = P(i,j)+1.d0
        IF(doshells) THEN
          S(i,j) = S(i,j)+1.d0
        ENDIF
        nloop=nloop+1
      ENDDO
    ENDDO
  ENDIF
ENDDO
!$OMP END PARALLEL DO
!
!If coordinates were not originally reduced,
!transform back to cartesian coordinates
IF(.NOT.isreduced) THEN
  CALL FRAC2CART(P,H)
  IF(doshells) THEN
    CALL FRAC2CART(S,H)
  ENDIF
ENDIF
!
!
!
300 CONTINUE
CALL ATOMSK_MSG(2094,(/''/),(/DBLE(NPwrap)/))
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE WRAP_XYZ
!
!
!
END MODULE wrap
