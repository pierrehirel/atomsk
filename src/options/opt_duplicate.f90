MODULE duplicate
!
!**********************************************************************************
!*  DUPLICATE                                                                     *
!**********************************************************************************
!* This module reads atomic coordinates from an array P, and                      *
!* expands it along the directions of the cell vectors.                           *
!**********************************************************************************
!* (C) March 2010 - Pierre Hirel                                                  *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 13 June 2016                                     *
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
CONTAINS
!
!
SUBROUTINE DUPLICATECELL(H,P,S,dupmatrix,SELECT,AUX)
!
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
LOGICAL:: doshells
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: i, newNP
INTEGER:: m, n, o, qi
INTEGER, DIMENSION(3):: dupmatrix
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Q, T
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX !auxiliary properties
!
!Initialize variables
i = 0
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(T)) DEALLOCATE(T)
IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
  doshells=.TRUE.
ELSE
  doshells=.FALSE.
ENDIF
!
CALL ATOMSK_MSG(2066,(/''/),(/ DBLE(dupmatrix(1)), &
     & DBLE(dupmatrix(2)), DBLE(dupmatrix(3)) /))
!
!If expansion is zero along a direction, correct it
DO i=1,3
  IF( dupmatrix(i)==0) THEN
    dupmatrix(i)=1
  ENDIF
ENDDO

!If all dimensions are set to 1 then the system stays the same
IF( dupmatrix(1)==1 .AND. dupmatrix(2)==1 .AND.                 &
  & dupmatrix(3)==1                                ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2728,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
200 CONTINUE
IF( .NOT.ALLOCATED(SELECT) ) THEN
  !All atoms must be duplicated
  newNP = SIZE(P,1)*ABS(dupmatrix(1)*dupmatrix(2)*dupmatrix(3))
ELSE
  !Only selected atoms will be duplicated
  !Count how many atoms are selected in the original system
  qi=0
  DO i=1,SIZE(SELECT)
    IF(SELECT(i)) qi=qi+1
  ENDDO
  newNP = qi*ABS(dupmatrix(1)*dupmatrix(2)*dupmatrix(3))
ENDIF
WRITE(msg,*) "new NP = ", newNP
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ALLOCATE(Q(newNP,4))
IF(ALLOCATED(AUX)) ALLOCATE( newAUX(newNP,SIZE(AUX,2) ) )
CALL ATOMSK_MSG(2067,(/''/),(/DBLE(newNP)/))
!
IF( newNP<=0 ) THEN
  nerr=nerr+1
  GOTO 500
ENDIF
!
IF( doshells ) THEN
  ALLOCATE(T(newNP,4))
ENDIF
!
!
qi = 0
DO o = 0 , dupmatrix(3)-SIGN(1,dupmatrix(3)) , SIGN(1,dupmatrix(3))
  DO n = 0 , dupmatrix(2)-SIGN(1,dupmatrix(2)) , SIGN(1,dupmatrix(2))
    DO m = 0 , dupmatrix(1)-SIGN(1,dupmatrix(1)) , SIGN(1,dupmatrix(1))
      DO i=1,SIZE(P,1)
        IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
          qi = qi+1
          Q(qi,1) = P(i,1)*SIGN(1,dupmatrix(1)) + DBLE(m)*H(1,1) + DBLE(n)*H(2,1) + DBLE(o)*H(3,1)
          Q(qi,2) = P(i,2)*SIGN(1,dupmatrix(2)) + DBLE(m)*H(1,2) + DBLE(n)*H(2,2) + DBLE(o)*H(3,2)
          Q(qi,3) = P(i,3)*SIGN(1,dupmatrix(3)) + DBLE(m)*H(1,3) + DBLE(n)*H(2,3) + DBLE(o)*H(3,3)
          Q(qi,4) = P(i,4)
          !Duplicated particles will have same auxiliary properties as the originals
          IF(ALLOCATED(newAUX)) newAUX(qi,:) = AUX(i,:)
          !Also duplicate shells if any
          IF( doshells ) THEN
            T(qi,1) = S(i,1)*SIGN(1,dupmatrix(1)) + DBLE(m)*H(1,1) + DBLE(n)*H(2,1) + DBLE(o)*H(3,1)
            T(qi,2) = S(i,2)*SIGN(1,dupmatrix(2)) + DBLE(m)*H(1,2) + DBLE(n)*H(2,2) + DBLE(o)*H(3,2)
            T(qi,3) = S(i,3)*SIGN(1,dupmatrix(3)) + DBLE(m)*H(1,3) + DBLE(n)*H(2,3) + DBLE(o)*H(3,3)
            T(qi,4) = S(i,4)
          ENDIF
        ENDIF
      ENDDO
    ENDDO
  ENDDO
ENDDO
!
!Replace old P with the new Q
DEALLOCATE(P)
ALLOCATE(P(SIZE(Q,1),4))
P = Q
DEALLOCATE(Q)
!Replace old AUX by newAUX
IF(ALLOCATED(AUX)) THEN
  DEALLOCATE(AUX)
  ALLOCATE( AUX( SIZE(newAUX,1), SIZE(newAUX,2) ) )
  AUX = newAUX
  IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
ENDIF
!
!Replace old S with the new T
IF( doshells) THEN
  DEALLOCATE(S)
  ALLOCATE(S(SIZE(T(:,1)),4))
  S(:,:) = T(:,:)
  DEALLOCATE(T)
ENDIF
!
!Resize the cell dimensions
H(1,:) = dupmatrix(1)*H(1,:)
H(2,:) = dupmatrix(2)*H(2,:)
H(3,:) = dupmatrix(3)*H(3,:)
!
!
!
500 CONTINUE
CALL ATOMSK_MSG(2068,(/''/),(/DBLE(newNP)/))
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
!
END SUBROUTINE DUPLICATECELL
!
!
END MODULE duplicate