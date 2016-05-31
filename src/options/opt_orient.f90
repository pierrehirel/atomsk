MODULE ORIENT
!
!**********************************************************************************
!*  ORIENT                                                                        *
!**********************************************************************************
!* This module performs a base transformation of an atomic system,                *
!* provided the system is oriented in a starting base Hstart, and the             *
!* user wants it to be oriented in the base Hend.                                 *
!**********************************************************************************
!* (C) July 2010 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 30 May 2016                                      *
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
USE functions
USE messages
USE files
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE ORIENT_XYZ(H,P,S,Hstart,Hend,SELECT,C_tensor)
!
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
LOGICAL,DIMENSION(3):: orthovec  !are vectors orthogonal?
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT  !mask for atom list
INTEGER:: i, j
REAL(dp):: H1, H2, H3
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: Hstart, Hend  !Initial and final orientations
REAL(dp),DIMENSION(3,3):: Gstart   !inverse of Hstart
REAL(dp),DIMENSION(3,3):: rot_matrix  !rotation matrix
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S  !atom positions
!
!
!Initialize variables
orthovec(:) = .FALSE.
i = 0
Gstart(:,:) = 0.d0
!
msg = 'Entering ORIENT_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
!
CALL ATOMSK_MSG(2071,(/''/),(/0.d0/))
!
!
IF(verbosity==4) THEN
  !Some debug messages
  msg = 'Hstart:'
  CALL ATOMSK_MSG(999,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
  DO i=1,3
    WRITE(msg,'(3(f6.2,2X))') (Hstart(i,j), j=1,3)
    CALL ATOMSK_MSG(999,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
  ENDDO
  msg = 'Hend:'
  CALL ATOMSK_MSG(999,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
  DO i=1,3
    WRITE(msg,'(3(f6.2,2X))') (Hend(i,j), j=1,3)
    CALL ATOMSK_MSG(999,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
  ENDDO
  WRITE(msg,*) 'Angles in Hstart, Hend (degrees):'
  CALL ATOMSK_MSG(999,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
  WRITE(msg,*) 'alpha: ', 180.d0*ANGVEC(Hstart(2,:),Hstart(3,:))/pi, &
      &                     180.d0*ANGVEC(Hend(2,:),Hend(3,:))/pi
  CALL ATOMSK_MSG(999,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
  WRITE(msg,*) 'beta:  ', 180.d0*ANGVEC(Hstart(3,:),Hstart(1,:))/pi, &
      &                     180.d0*ANGVEC(Hend(3,:),Hend(1,:))/pi
  CALL ATOMSK_MSG(999,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
  WRITE(msg,*) 'gamma: ', 180.d0*ANGVEC(Hstart(1,:),Hstart(2,:))/pi, &
      &                     180.d0*ANGVEC(Hend(1,:),Hend(2,:))/pi
  CALL ATOMSK_MSG(999,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
ENDIF
!
!
!
100 CONTINUE
!If bases are the same we skip it
!For this we compute the angle between vectors (Hx,H'x), (Hy,H'y), (Hz,H'z)
!If all angles are zero then the system is not rotated at all
IF( DABS( ANGVEC(Hstart(1,:),Hend(1,:)) ) <1.d-6 .AND.          &
  & DABS( ANGVEC(Hstart(2,:),Hend(2,:)) ) <1.d-6 .AND.          &
  & DABS( ANGVEC(Hstart(3,:),Hend(3,:)) ) <1.d-6       ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2731,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Normalize the vectors in each matrix
DO i=1,3
  IF( VECLENGTH(Hstart(i,:)) .NE. 0.d0 ) THEN
    Hstart(i,:) = Hstart(i,:)/VECLENGTH(Hstart(i,:))
  ELSE
    !we have a problem
    nerr = nerr+1
    GOTO 1000
  ENDIF
  !
  IF( VECLENGTH(Hend(i,:)) .NE. 0.d0 ) THEN
    Hend(i,:) = Hend(i,:)/VECLENGTH(Hend(i,:))
  ELSE
    !we have a problem
    nerr = nerr+1
    GOTO 1000
  ENDIF
ENDDO
!
!
!Check that vectors of end matrix form the same angles as vectors in initial matrix
orthovec(:) = .FALSE.
IF( DABS( ANGVEC(Hstart(1,:),Hstart(2,:))-ANGVEC(Hend(1,:),Hend(2,:)) ) < 1.d-6 ) orthovec(1)=.TRUE.
IF( DABS( ANGVEC(Hstart(2,:),Hstart(3,:))-ANGVEC(Hend(2,:),Hend(3,:)) ) < 1.d-6 ) orthovec(2)=.TRUE.
IF( DABS( ANGVEC(Hstart(3,:),Hstart(1,:))-ANGVEC(Hend(3,:),Hend(1,:)) ) < 1.d-6 ) orthovec(3)=.TRUE.
WRITE(msg,*) "orthovec: ", orthovec(:)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
IF( ANY(.NOT.orthovec) ) THEN
  CALL ATOMSK_MSG(2801,(/''/),(/0.d0/))
  nerr = nerr+1
  GOTO 1000
ENDIF
!
!
!
200 CONTINUE
!Perform the transformation:
!First, convert to fractional coordinates
CALL CART2FRAC(P,H)
IF( ALLOCATED(S) ) THEN
  CALL CART2FRAC(S,H)
ENDIF
!
!Second, save the new base vectors
CALL INVMAT(Hstart,Gstart)
DO i=1,3
  H1 = H(i,1)
  H2 = H(i,2)
  H3 = H(i,3)
  H(i,1) = Gstart(1,1)*H1 + Gstart(1,2)*H2 + Gstart(1,3)*H3
  H(i,2) = Gstart(2,1)*H1 + Gstart(2,2)*H2 + Gstart(2,3)*H3
  H(i,3) = Gstart(3,1)*H1 + Gstart(3,2)*H2 + Gstart(3,3)*H3
  H1 = H(i,1)
  H2 = H(i,2)
  H3 = H(i,3)
  H(i,1) = H1*Hend(1,1) + H2*Hend(1,2) + H3*Hend(1,3)
  H(i,2) = H1*Hend(2,1) + H2*Hend(2,2) + H3*Hend(2,3)
  H(i,3) = H1*Hend(3,1) + H2*Hend(3,2) + H3*Hend(3,3)
ENDDO
!
!Third, convert back to cartesian coordinates
CALL FRAC2CART(P,H)
IF( ALLOCATED(S) ) THEN
  CALL FRAC2CART(S,H)
ENDIF
!
!
!If elastic tensor is set, rotate it
IF( C_tensor(1,1).NE.0.d0 ) THEN
  !Set up rotation matrix
  CALL INVMAT(Hstart,Gstart)
  rot_matrix = MATMUL(Gstart,Hend)
  !
  C_tensor = ROTELAST( C_tensor, rot_matrix )
  !
  CALL ATOMSK_MSG(2099,(/""/),(/0.d0/))
ENDIF
!
!
!
400 CONTINUE
CALL ATOMSK_MSG(2072,(/''/),(/0.d0/))
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
END SUBROUTINE ORIENT_XYZ
!
!
END MODULE orient
