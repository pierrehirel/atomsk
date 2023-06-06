MODULE elasticity
!
!**********************************************************************************
!*  ELASTICITY                                                                    *
!**********************************************************************************
!* This module contains routines related to elasticity and elastic tensor.        *
!**********************************************************************************
!* (C) June 2023                                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 06 June 2023                                     *
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
!* List of subroutines in this module:                                            *
!* ELAST2TENSOR        converts from Voigt notation to full elastic tensor        *
!* CHECK_CTENSOR       checks if an elastic tensor is symmetric                   *
!* ELASTINDEX          reduces indices (i,j) into index m for 9x9 matrices        *
!* ELAST2INDEX         convert index m into indices (i,j) for 9x9 matrices        *
!* ROTELAST            rotates a 9x9 matrix                                       *
!**********************************************************************************
!
!
USE comv
USE constants
USE functions
USE math
USE subroutines
!
IMPLICIT NONE
!
!
CONTAINS
!
!
!********************************************************
!  ELAST2TENSOR
!  This subroutine sets up the 9x9 elastic tensor,
!  provided the 9 (Voigt notation) elastic constants.
!********************************************************
!
SUBROUTINE ELAST2TENSOR(elcst,eltens)
!
IMPLICIT NONE
REAL(dp),DIMENSION(9),INTENT(IN):: elcst !C11,C22,C33,C23,C31,C12,C44,C55,C66
REAL(dp),DIMENSION(:,:),INTENT(INOUT):: eltens
!
eltens(:,:) = 0.d0
!
eltens(1,1) = elcst(1)
eltens(1,2) = elcst(6)
eltens(1,3) = elcst(5)
eltens(2,1) = elcst(6)
eltens(2,2) = elcst(2)
eltens(2,3) = elcst(4)
eltens(3,1) = elcst(5)
eltens(3,2) = elcst(4)
eltens(3,3) = elcst(3)
eltens(4,4) = elcst(7)
eltens(5,5) = elcst(8)
 eltens(6,6) = elcst(9)
IF( SIZE(eltens,1)==9 .AND. SIZE(eltens,2)==9 ) THEN
  !Fill the rest by symmetry
  eltens(1:3,7:9) = eltens(1:3,4:6)
  eltens(4:6,7:9) = eltens(4:6,4:6)
  eltens(7:9,1:3) = eltens(4:6,1:3)
  eltens(7:9,4:6) = eltens(4:6,4:6)
  eltens(7:9,7:9) = eltens(4:6,4:6)
ENDIF

!
END SUBROUTINE ELAST2TENSOR
!
!
!********************************************************
!  CHECK_CTENSOR
!  This subroutine checks if an elastic tensor is
!  symmetric.
!********************************************************
!
SUBROUTINE CHECK_CTENSOR(C_tensor,status)
!
IMPLICIT NONE
INTEGER:: i, j
INTEGER:: status  !=0 if tensor is OK
                  !=1 if tensor is not symmetric, i.e. C(i,j).NE.C(j,i)
REAL(dp),DIMENSION(9,9),INTENT(IN):: C_tensor
!
status = 0
DO i=2,9
  DO j=1,i-1
    IF( DABS(C_tensor(i,j)-C_tensor(j,i))>1.d-3 ) THEN
      status=1
    ENDIF
  ENDDO
ENDDO
!
END SUBROUTINE CHECK_CTENSOR
!
!
!********************************************************
!  ELASTINDEX
!  This function converts the indices i and j into
!  single index m.
!********************************************************
FUNCTION ELASTINDEX(i,j) RESULT(m)
!
IMPLICIT NONE
INTEGER, INTENT(IN):: i, j
INTEGER:: m
!
m = 0
IF(i==j) THEN
  m = i
!
ELSEIF(i==1) THEN
  IF(j==2) THEN
    m = 6
  ELSEIF(j==3) THEN
    m = 8
  ENDIF
!
ELSEIF(i==2) THEN
  IF(j==1) THEN
    m = 9
  ELSEIF(j==3) THEN
    m = 4
  ENDIF
!
ELSEIF(i==3) THEN
  IF(j==1) THEN
    m = 5
  ELSEIF(j==2) THEN
    m = 7
  ENDIF
ENDIF
!
END FUNCTION ELASTINDEX
!
!
!********************************************************
!  ELAST2INDEX
!  This function converts the index i into indices m and n.
!  This is for use by the function ROTELAST below.
!********************************************************
FUNCTION ELAST2INDEX(i) RESULT(mn)
!
IMPLICIT NONE
INTEGER,INTENT(IN):: i
INTEGER,DIMENSION(2):: mn !mn(1)=m and mn(2)=n
!
mn(:) = 0
!
SELECT CASE(i)
!
CASE(1,2,3)
  mn(1) = i
  mn(2) = i
CASE(4)
  mn(1) = 2
  mn(2) = 3
CASE(5)
  mn(1) = 3
  mn(2) = 1
CASE(6)
  mn(1) = 1
  mn(2) = 2
CASE(7)
  mn(1) = 3
  mn(2) = 2
CASE(8)
  mn(1) = 1
  mn(2) = 3
CASE(9)
  mn(1) = 2
  mn(2) = 1
CASE DEFAULT
  PRINT*, "Problem with elastindex"
END SELECT
!
RETURN
!
END FUNCTION ELAST2INDEX
!
!
!********************************************************
!  ROTELAST
!  This function applies a rotation defined by a
!  3x3 matrix to a 9x9 tensor. Although this operation
!  is general, it is intended to rotate the elastic tensor.
!********************************************************
FUNCTION ROTELAST(ELTENS,T) RESULT(newELTENS)
!
IMPLICIT NONE
INTEGER:: i, j, k, l, m, n
INTEGER,DIMENSION(2):: mn
REAL(dp), DIMENSION(9,9):: ELTENS    !Elastic tensor before
REAL(dp), DIMENSION(9,9):: newELTENS !Elastic tensor after
REAL(dp), DIMENSION(3,3), INTENT(IN):: T  !Rotation matrix (3x3)
REAL(dp), DIMENSION(9,9):: Q  !Rotation matrix (9x9)
!
!Initialize variables
mn(:) = 0
newELTENS(:,:) = 0.d0
Q(:,:) = 0.d0
!
!Build the 9x9 rotation matrix ROT9
DO i=1,9
  mn = ELAST2INDEX(i)
  m = mn(1)
  n = mn(2)
  DO j=1,9
    mn = ELAST2INDEX(j)
    k = mn(1)
    l = mn(2)
    Q(i,j) = T(k,m)*T(l,n)
  ENDDO
ENDDO
!
!Rotate the elastic tensor:  C' = tQ * C * Q
newELTENS = MATMUL( TRANSPOSE(Q) , MATMUL(ELTENS,Q) )
!
RETURN
!
END FUNCTION ROTELAST
!
!
!
END MODULE elasticity
