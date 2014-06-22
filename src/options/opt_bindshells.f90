MODULE bindshells
!
!**********************************************************************************
!*  BINDSHELLS                                                                    *
!**********************************************************************************
!* This module re-associates shells (in the sens of a core-shell model) with      *
!* their core.                                                                    *
!**********************************************************************************
!* (C) February 2012 - Pierre Hirel                                               *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 19 Feb. 2014                                     *
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
USE neighbors
USE files
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE BSHELLS_XYZ(H,P,S)
!
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
INTEGER:: i
INTEGER:: Nbound !number of shells that were re-bound
INTEGER,DIMENSION(:),ALLOCATABLE:: Nlist !index of nearest shell
LOGICAL:: exceeds100 !does the number of neighbours exceed 100?
REAL(dp),DIMENSION(4):: Stemp !temporary position of shell
REAL(dp),DIMENSION(3,3),INTENT(IN):: H    !supercell parameters
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P     !positions of cores
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: S  !positions of shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: V_NN  !position of nearest shell
!
!Initialize variables
Nbound=0
!
!
msg = 'Entering BSHELLS_XYZ'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
CALL ATOMSK_MSG(2105,(/''/),(/0.d0/))
!
IF( .NOT.ALLOCATED(S) ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2744,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
100 CONTINUE
DO i=1,SIZE(P,1)
  Stemp(:) = 0.d0
  !Search the shell that is closest to that atom
  CALL FIND_NNN(H,S,P(i,:),1,V_NN,Nlist,exceeds100)
  !V_NN(1,:) = position of closest shell
  !Nlist(1) = index of that shell
  !
  !Check that we found a shell that has a different index than i,
  !and that it is actually a shell
  IF( Nlist(1).NE.i .AND. NINT(S(Nlist(1),4)).NE.0 ) THEN
    !This shell was not associated with this core before
    IF( VECLENGTH( V_NN(1,1:3)-P(i,1:3) ) < 1.5d0 ) THEN
      !If distance is smaller than 1.5 A,then
      !this shell belongs to that core
      IF( NINT(S(i,4)).NE.0 ) THEN
        !Present core already had a shell: save this shell position
        Stemp(:) = S(i,:)
      ENDIF
      !
      !Save closest shell position at current index i
      S(i,1:3) = V_NN(1,1:3)
      S(i,4) = P(i,4)
      !
      IF( NINT(Stemp(4)).NE.0 ) THEN
        !Save previous shell of atom i to index Nlist(1)
        S(Nlist(1),:) = Stemp(:)
        IF( VECLENGTH( S(Nlist(1),1:3)-P(Nlist(1),1:3) ) < 1.5d0 ) THEN
          !This shell belongs to the core at Nlist(1)
          Nbound = Nbound+1
        ENDIF
      ELSE
        !Otherwise, simply delete previous shell position
        S(Nlist(1),:) = 0.d0
      ENDIF
      !
      !One more shell was associated
      Nbound = Nbound+1
    ELSE
      !Otherwise (distance>1.5A) this core has no shell
      S(i,:) = 0.d0
    ENDIF
  ENDIF
ENDDO
!
CALL ATOMSK_MSG(2106,(/''/),(/DBLE(Nbound)/))
!
!
!
1000 CONTINUE
IF(ALLOCATED(Nlist)) DEALLOCATE(Nlist)
IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
!
!
!
END SUBROUTINE BSHELLS_XYZ
!
!
!
END MODULE bindshells