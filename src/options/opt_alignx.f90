MODULE alignx
!
!**********************************************************************************
!*  ALIGNX                                                                        *
!**********************************************************************************
!* This module reads atomic positions from an array P and                         *
!* rotates the system so the first vector is aligned with X and the               *
!* second vector is in the XY plane, i.e. so the H matrix is transformed          *
!* into a lower triangular matrix H':                                             *
!*                     |  H'(1,1)     0        0     |                            *
!*                H' = |  H'(2,1)  H'(2,2)     0     |                            *
!*                     |  H'(3,1)  H'(3,2)  H'(3,3)  |                            *
!* Atomic positions are also transformed accordingly.                             *
!**********************************************************************************
!* (C) May 2010 - Pierre Hirel                                                    *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 06 April 2022                                    *
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
USE crystallography
USE functions
USE messages
USE files
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE ALIGN_X(H,P,S,C_tensor)
!
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
INTEGER:: i
REAL(dp):: alpha, beta, gamma !cell angles
REAL(dp):: Hx, Hy, Hz  !Length of H(1,:), H(2,:) and H(3,:)
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: Hstart, Hend !initial/final H, normalized
REAL(dp),DIMENSION(3,3):: Gend !inverse of Hend
REAL(dp),DIMENSION(3,3):: rot_matrix !rotation matrix (only for rotating C_tensor)
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P,S  !positions of atoms, shells
!
!
!Initialize variables
i = 0
alpha = 0.d0
beta  = 0.d0
gamma = 0.d0
Hx = 0.d0
Hy = 0.d0
Hz = 0.d0
Hstart = H  !for now don't normalize anything
!
msg = 'Entering ALIGN_X'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
CALL ATOMSK_MSG(2050,(/''/),(/0.d0/))
!
!
!
100 CONTINUE
IF( H(1,2)==0.d0 .AND. H(1,3)==0.d0 .AND. H(2,3)==0.d0 ) THEN
  !Vectors are already aligned => skip it
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2720,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
200 CONTINUE
!Convert to fractional coordinates
CALL CART2FRAC(P,H)
IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
  CALL CART2FRAC(S,H)
ENDIF
!
!Convert the matrix H into conventional notation
Hx = VECLENGTH(H(1,:))
Hy = VECLENGTH(H(2,:))
Hz = VECLENGTH(H(3,:))
alpha = ANGVEC(H(2,:),H(3,:))
beta  = ANGVEC(H(3,:),H(1,:))
gamma = ANGVEC(H(1,:),H(2,:))
!
!Then convert this conventional notation into lower-triangular matrix H
CALL CONVMAT(Hx,Hy,Hz,alpha,beta,gamma,H)
!
!Convert back to cartesian coordinates
CALL FRAC2CART(P,H)
IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
  CALL FRAC2CART(S,H)
ENDIF
!
!
!If elastic tensor is set, rotate it
IF( C_tensor(1,1).NE.0.d0 ) THEN
  !Normalize initial H (which was saved in Hstart)
  DO i=1,3
    Hstart(i,:) = Hstart(i,:)/VECLENGTH(Hstart(i,:))
  ENDDO
  !Normalize final H -> Hend
  DO i=1,3
    Hend(i,:) = H(i,:)/VECLENGTH(H(i,:))
  ENDDO
  !Compute Gend = inverse of Hend
  CALL INVMAT(Hend,Gend)
  !Define rotation matrix
  rot_matrix = MATMUL(Gend,Hstart)
  !
  C_tensor = ROTELAST( C_tensor, rot_matrix )
  !
  CALL ATOMSK_MSG(2099,(/""/),(/0.d0/))
ENDIF
!
!
CALL ATOMSK_MSG(2051,(/''/),(/0.d0/))
!
!
!
1000 CONTINUE
!
!
!
END SUBROUTINE ALIGN_X
!
!
END MODULE alignx
