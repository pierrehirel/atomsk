MODULE dislocation_aniso
!
!**********************************************************************************
!*  DISLOCATION_ANISO                                                             *
!**********************************************************************************
!* This module contains functions computing the displacements due to              *
!* a straight dislocation in the framework of anisotropic elasticity.             *
!* Formulae can be found for instance in:                                         *
!* J.P. Hirth, J. Lothe, 'Theory of dislocations', 1st ed. (1968), p.426.         *
!* These functions are used by the module "opt_disloc.f90".                       *
!**********************************************************************************
!* (C) January 2018 - Pierre Hirel                                                *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 10 Sept. 2019                                    *
!**********************************************************************************
!* List of subroutines in this module:                                            *
!* ANISO_DISP          applies anisotropic disp. of a disloc. to 1 atom           *
!* ANISO_COEFF         computes coefficients for disp. in anisotropic medium      *
!* ANISO_STRESS        computes stresses due to dislocation                       *
!* List of functions in this module:                                              *
!* ANISO_EFACTOR       computes prelogarithmic energy factor                      *
!* DET_COMPMAT         computes the determinant of a 2x2 complex matrix           *
!* DIAGMUL             multiplies elements of a diagonal of a 9x3 complex array   *
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
USE sorting
USE subroutines
!
!
CONTAINS
!
!
!********************************************************
! ANISO_DISP
! This function calculates the displacements due to a
! dislocation in an anisotropic medium:
!   u(k) = Re{ [-1/(2*i*pi)]*[SUM(n=1,3) A_k(n)*D(n)
!                                  *Ln(x1+P(n)*x2)] }
! The complex coefficients P(n), A_k(n), D(n) must
! be provided as input (see routine ANISO_COEFF below).
!********************************************************
!
FUNCTION ANISO_DISP(i,P,a1,a2,a3,pos1,pos2,A_kn,Dn,Pn) RESULT(disp)
!
IMPLICIT NONE
INTEGER,INTENT(IN):: a1, a2, a3
INTEGER,INTENT(IN):: i            !index of atom to be displaced
INTEGER:: k, n
REAL(dp),DIMENSION(3):: disp
REAL(dp),INTENT(IN):: pos1, pos2  !Position of the dislocation in the plane
                                  !orthogonal to the dislocation line
REAL(dp),DIMENSION(3),INTENT(IN):: P  !Atom position
COMPLEX(dp):: logterm, tempcmplx !the term that is in the neperian log
COMPLEX(dp),PARAMETER:: frac2ipi = DCMPLX(0.d0,-1.d0/(2.d0*pi))  ! -1/(2*i*pi)
COMPLEX(dp),DIMENSION(3),INTENT(IN):: Dn, Pn !anisotropy coefficients D(n), P(n)-
COMPLEX(dp),DIMENSION(3,3),INTENT(IN):: A_kn !-and A_k(n)
!
disp(:) = 0.d0
!
DO k=1,3
  tempcmplx = DCMPLX(0.d0,0.d0)
  !Compute the sum
  DO n=1,3
    logterm = DCMPLX(P(a1)-pos1,0.d0) + Pn(n)*DCMPLX(P(a2)-pos2,0.d0)
    tempcmplx = tempcmplx + frac2ipi*A_kn(k,n)*Dn(n)*LOG(logterm)
  ENDDO
  !We want only the real part
  disp(k) = DBLE(tempcmplx)
ENDDO
!Message if displacement was too large
IF( VECLENGTH(disp(:))>=10.d0 ) CALL ATOMSK_MSG(2727,(/''/),(/DBLE(i)/))
!
END FUNCTION ANISO_DISP
!
!
!********************************************************
! ANISO_COEFF
! This subroutine computes the complex coefficients
! A_k(n), D(n) and P(n) that determine the displacements
! due to a dislocation in an anisotropic medium,
! provided the elastic tensor for that medium and the
! Burgers vector of the dislocation.
! The method is fully explained in:
! J.P. Hirth, J. Lothe, 'Theory of dislocations',
! 1st ed. (1968), p.418.
! The numerotation of equations below follow
! this reference.
! NOTE: the elastic tensor "C_tensor" provided as input
! is assumed to correspond to the current orientation of
! the system, i.e. any necessary rotation must be
! performed BEFORE calling this routine. Furthermore
! it is recommended to provide values in GPa to avoid
! numerical precision issues (addition/multiplication
! of huge and small values).
! NOTE2: as in the aforementioned book, this routine
! assumes that the dislocation line is along the Z axis,
! therefore the k index of coefficients will follow
! the order X (k=1), Y (k=2), Z (k=3).
! If the dislocation line is not along Z, make sure to
! rotate b and C_tensor BEFORE calling this routine.
!********************************************************
!
SUBROUTINE ANISO_COEFF(b,C_tensor,A_kn,Dn,Pn,B_ijk,ifail)
!
IMPLICIT NONE
!Input variables
REAL(dp),DIMENSION(3),INTENT(IN):: b !Burgers vector of the dislocation
REAL(dp),DIMENSION(9,9),INTENT(IN):: C_tensor  !Elastic tensor
!
!Internal variables
CHARACTER(LEN=128):: msg
INTEGER:: i, j, k, l, n, p, q, r
INTEGER,DIMENSION(6):: IPIV !for LAPACK routine DGESV
INTEGER,DIMENSION(:),ALLOCATABLE:: newindex  !list of index after sorting
REAL(dp),PARAMETER:: Ak_threshold=1.d-9 !threshold for normalization of the A_k(n)
                                      !This is because of numerical precision,
                                      !to avoid division by ridiculously small numbers.
                                      !The value 10^-4 assumes C_tensor is in GPa
REAL(dp),DIMENSION(3):: PIm_sign !sign of imaginary part of P(n)
REAL(dp),DIMENSION(6):: WR, WI !real and imaginary parts of eigenvalues of aik_Hess
REAL(dp),DIMENSION(7):: aik_det !factors for sextic equation (13-85)
REAL(dp),DIMENSION(18):: WORK !For LAPACK routine DGEEV
REAL(dp),DIMENSION(6,1):: RHA !Right-Hand Array for Eq.(13-88) and (13-89)
REAL(dp),DIMENSION(6,2):: aik_roots !complex roots of Eq. (13-85) (real/imaginary parts)
REAL(dp),DIMENSION(6,6):: LHA !Left-Hand Array for Eq.(13-88) and (13-89)
REAL(dp),DIMENSION(6,6):: aik_Hess !companion Hesseberg matrix of Eq.(13-85)
REAL(dp),DIMENSION(6,7):: cik_diag !diagonal products of c_ik
REAL(dp),DIMENSION(9,3):: c_ik !contains the c_i1k1, c_i2k1+c_i1k2, and c_i2k2
REAL(dp),DIMENSION(1,6):: VL, VR  !For LAPACK routine DGEEV
COMPLEX(dp):: A_1, A_2, A_3 !subdeterminants of a_ik(n)
COMPLEX(dp):: tempcmplx
COMPLEX(dp),DIMENSION(2,2):: tempmat !temporary 2x2 matrix
COMPLEX(dp),DIMENSION(3,3,3):: a_ik  !the a_ik(n)
!
!Output variables
COMPLEX(dp),DIMENSION(3),INTENT(OUT):: Dn, Pn !the final D(n) and P(n)
COMPLEX(dp),DIMENSION(3,3),INTENT(OUT):: A_kn !the final A_k(n)
COMPLEX(dp),DIMENSION(9,3,3),INTENT(OUT):: B_ijk  !the final B_ijk(n)
INTEGER,INTENT(OUT),OPTIONAL:: ifail !=0 if the routine succeeds
                                     !=1 if roots of Eq.(13-85) cannot be found
                                     !=2 if the A_k(n) cannot be calculated
                                     !=3 if the linear eq. giving D(n) cannot be solved
!
!Initialize variables
ifail = 0 !so far the routine is successful
IPIV(:) = 0
A_kn(:,:) = DCMPLX(0.d0,0.d0)
Dn(:) = DCMPLX(0.d0,0.d0)
Pn(:) = DCMPLX(0.d0,0.d0)
B_ijk(:,:,:) = DCMPLX(0.d0,0.d0)
!
!
IF(verbosity==4) THEN
  !Some debug messages
  !Write elastic tensor to logfile
  msg = 'Provided elastic tensor (GPa):'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,9
    WRITE(msg,'(9(e10.3,2X))') (C_tensor(i,j), j=1,9)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
  !Write Burgers vector to logfile
  WRITE(msg,'(a16,3e10.3)') 'Burgers vector: ', b(1), b(2), b(3)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!
!
!
100 CONTINUE
! 1. Solve the equation:
!      |{a_ik(n)}| = 0     (13-85)
!    where:
!      a_ik(n) = c_i1k1 + (c_i1k2+c_i2k1)*P(n) + c_i2k2*P(n)^2
!    Since a_ik(n) is a 3x3 matrix and we need to calculate its determinant,
!    the roots of an order-6 (or sextic) polynomial must be found.
!
!First, convert the c_ijkl into c_ik
!using the sudoku of elastic tensor:
!       
!                   i or k
!            ||  1  |  2  |  3
!      ======++=====+=====+======
!        1   ||  1  |  9  |  5
!  j   ------++------------------
!  or    2   ||  6  |  2  |  7
!  l   ------++------------------
!        3   ||  8  |  4  |  3
!
!E.g. c_1111 becomes c(1,1), c_2312 becomes c(4,6), etc.
!Then apply symmetry considerations like c(9,1)=c(1,6); c(9,7)=c(4,6) etc.
!This allows to convert the c_i1k1 etc. into the c_ik below.
!You can do it by hand it's a very good exercise.
!
!The c_ik(:,1) are the c_i1k1
!The c_ik(:,2) are the c_i1k2+c_i2k1
!The c_ik(:,3) are the c_i2k2
  c_ik(:,:) = 0.d0
 c_ik(1,1) = C_tensor(1,1)
 c_ik(1,2) = C_tensor(1,6)*2.d0
 c_ik(1,3) = C_tensor(6,6)
  c_ik(2,1) = C_tensor(1,6)
  c_ik(2,2) = C_tensor(1,2)+C_tensor(6,6)
  c_ik(2,3) = C_tensor(2,6)
 c_ik(3,1) = C_tensor(1,5)
 c_ik(3,2) = C_tensor(1,4)+C_tensor(5,6)
 c_ik(3,3) = C_tensor(4,6)
  c_ik(4,1) = C_tensor(1,6)
  c_ik(4,2) = C_tensor(6,6)+C_tensor(1,2)
  c_ik(4,3) = C_tensor(2,6)
 c_ik(5,1) = C_tensor(6,6)
 c_ik(5,2) = C_tensor(2,6)*2.0
 c_ik(5,3) = C_tensor(2,2) 
  c_ik(6,1) = C_tensor(5,6)
  c_ik(6,2) = C_tensor(4,6)+C_tensor(2,5)
  c_ik(6,3) = C_tensor(2,4)
 c_ik(7,1) = C_tensor(1,5)
 c_ik(7,2) = C_tensor(5,6)+C_tensor(1,4)
 c_ik(7,3) = C_tensor(4,6)
  c_ik(8,1) = C_tensor(5,6)
  c_ik(8,2) = C_tensor(2,5)+C_tensor(4,6)
  c_ik(8,3) = C_tensor(2,4)
 c_ik(9,1) = C_tensor(5,5)
 c_ik(9,2) = C_tensor(4,5)*2.d0
 c_ik(9,3) = C_tensor(4,4)
!
IF(verbosity==4) THEN
  !Some debug messages
  msg = "c_ik(:,:) ="
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,9
    WRITE (msg,'(3X,3(e10.3,1X))') (c_ik(i,j), j=1,3)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
ENDIF
!
!Multiply elements in the diagonals of the matrix
!(DIAGMUL function is at end of this module)
 cik_diag(:,:) = 0.d0
 cik_diag(1,:) = DIAGMUL(c_ik,1,5,9)
 cik_diag(2,:) = -1.d0*DIAGMUL(c_ik,1,6,8)
 cik_diag(3,:) = DIAGMUL(c_ik,2,6,7)
 cik_diag(4,:) = -1.d0*DIAGMUL(c_ik,2,4,9)
 cik_diag(5,:) = DIAGMUL(c_ik,3,4,8)
 cik_diag(6,:) = -1.d0*DIAGMUL(c_ik,3,5,7)
!
!Compute the determinants = sum of diagonals
aik_det(:) = 0.d0
DO i=1,7
  aik_det(i) = SUM( cik_diag(:,i) )
ENDDO
!
!Now Eq.(13-85) looks like:
!   aik_det(1)*x^6 + aik_det(2)*x^5 + ... + aik_det(6)*x + aik_det(7) = 0
!Let's find the 6 complex roots!!
!
!Construct the companion matrix (which is a Hessenberg matrix)
aik_Hess(:,:) = 0.d0
DO i=1,5
  aik_Hess(i+1,i) = 1.d0
ENDDO
aik_Hess(1,6) = -1.d0*aik_det(7)/aik_det(1)
aik_Hess(2,6) = -1.d0*aik_det(6)/aik_det(1)
aik_Hess(3,6) = -1.d0*aik_det(5)/aik_det(1)
aik_Hess(4,6) = -1.d0*aik_det(4)/aik_det(1)
aik_Hess(5,6) = -1.d0*aik_det(3)/aik_det(1)
aik_Hess(6,6) = -1.d0*aik_det(2)/aik_det(1)
!
!Call LAPACK routine DGEEV to find the 6 complex eigenvalues of aik_Hess
!( eigenvalues of aik_Hess = roots of Eq.(13-85) )
CALL DGEEV( 'N', 'N', 6, aik_Hess, 6, WR, WI, VL, 1, VR,1, WORK, 18, k )
!If k is different from 0 then it failed => exit
IF(k.NE.0) THEN
  ifail = 1
  GOTO 400
ENDIF
!Otherwise WR and WI are the real and imaginary parts of the solutions.
!Conjugate pairs appear consecutively with WI>0 first
!Save it to the table aik_roots(:,:)
DO i=1,6
  aik_roots(i,1) = WR(i)
  aik_roots(i,2) = WI(i)
ENDDO
!
!The aik_roots(:,:) are now the 6 complex roots
!Sort them by increasing values of real parts
CALL BUBBLESORT(aik_roots,1,'up  ',newindex)
!
IF(verbosity==4) THEN
  !Some debug messages
  msg = "6 complex roots of |{a_ik(n)}| (sorted):"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,6
    WRITE(msg,'(3X,e12.5,a3,e12.5,a2)') aik_roots(i,1), &
         & " + ", aik_roots(i,2), " i"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
ENDIF
!
!Roots come in pairs of complex conjugates, i.e. a+ib and a-ib
!Keep only non-conjugate roots with positive imaginary parts => these are the P(n)
Pn(1) = DCMPLX( aik_roots(1,1), DABS(aik_roots(1,2)) )
Pn(2) = DCMPLX( aik_roots(3,1), DABS(aik_roots(3,2)) )
Pn(3) = DCMPLX( aik_roots(5,1), DABS(aik_roots(5,2)) )
!
IF(verbosity==4) THEN
  !Some debug messages
  msg = "Values of the P(n):"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO n=1,3
    WRITE(msg,'(3X,a2,i1,a4,e12.5,a3,e12.5,a2)') "P(", n, ") = ", &
         & DBLE(Pn(n)), " + ", AIMAG(Pn(n)), " i"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
ENDIF
!
!Now the roots P(n) are known, compute the actual a_ik(n):
!    a_ik(n) = c_i1k1 + (c_i1k2+c_i2k1)*P(n) + c_i2k2*P(n)^2
DO n=1,3
  k=0
  DO i=1,3
    IF(i==2) THEN
      k=3
    ELSEIF(i==3) THEN
      k=6
    ENDIF
    !
    DO j=1,3
      l=j+k
      a_ik(i,j,n) = c_ik(l,1) + c_ik(l,2)*Pn(n) + c_ik(l,3)*Pn(n)**2
    ENDDO
  ENDDO
ENDDO
!
IF(verbosity==4) THEN
  !Some debug messages
  msg = "Values of the a_ik(n):"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO n=1,3
    DO i=1,3
      DO k=1,3
        WRITE (msg,'(3X,a2,i1,i1,a1,i1,a4,e10.3,a3,e10.3,a2)') &
            & "a_", i, k, "(", n, ") = ", &
            & DBLE(a_ik(i,k,n)), " + ", AIMAG(a_ik(i,k,n)), " i"
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDDO
    ENDDO
  ENDDO
ENDIF
!
!
200 CONTINUE
! 2. Determine the A_k(n) by solving the set of equations:
!      a_ik(n).A_k(n) = 0      (13-87)
!
!    The coefficients A_k(n) are given by the subdeterminants
!    of the a_ik(n):
!
!                  | a12(n)  a13(n) |
!      A1(n) = det | a22(n)  a23(n) |
!
!                  | a11(n)  a13(n) |
!      A2(n) = det | a21(n)  a23(n) |
!
!                  | a11(n)  a12(n) |
!      A3(n) = det | a21(n)  a22(n) |
!
!    Then the A_k(n) must be normalized e.g. by dividing all
!    of them by A3(n) so that A3(n)=1. However there are
!    cases where A3(n)=0 which forbids this division. In
!    such cases one must divide by whichever is not zero among
!    A1(n) and A2(n).
!
DO n=1,3
  ! Compute the Ak(n)
  !A1(n)
  tempmat(1,1) = a_ik(1,2,n)
  tempmat(2,1) = a_ik(1,3,n)
  tempmat(1,2) = a_ik(2,2,n)
  tempmat(2,2) = a_ik(2,3,n)
  A_1 = DET_COMPMAT(tempmat)
  !A2(n)
  tempmat(1,1) = a_ik(1,1,n)
  tempmat(1,2) = a_ik(1,3,n)
  tempmat(2,1) = a_ik(2,1,n)
  tempmat(2,2) = a_ik(2,3,n)
  A_2 = DET_COMPMAT(tempmat)
  !A3(n)
  tempmat(1,1) = a_ik(1,1,n)
  tempmat(1,2) = a_ik(1,2,n)
  tempmat(2,1) = a_ik(2,1,n)
  tempmat(2,2) = a_ik(2,2,n)
  A_3 = DET_COMPMAT(tempmat)
  !
  IF(verbosity==4) THEN
    WRITE(msg,'(a24,i1,a2)') "Subdeterminants of a_ik(", n, "):"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,'(a6,i1,a4,e10.3,a3,e10.3,a2)') "   A1(", n, ") = ", &
         & DBLE(A_1), " + ", AIMAG(A_1), " i"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,'(a6,i1,a4,e10.3,a3,e10.3,a2)') "   A2(", n, ") = ", &
         & DBLE(A_2), " + ", AIMAG(A_2), " i"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,'(a6,i1,a4,e10.3,a3,e10.3,a2)') "   A3(", n, ") = ", &
         & DBLE(A_3), " + ", AIMAG(A_3), " i"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
  !
  IF( DABS(DBLE(A_3)) > Ak_threshold ) THEN
    !We are not dividing by zero => normalize the A_k(n) to A3(n)
    msg = "Normalizing by A3(n)..."
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !In this case A3(n)=1
    A_kn(1,n) = A_1/A_3
    A_kn(2,n) = -1.d0*A_2/A_3
    A_kn(3,n) = DCMPLX(1.d0,0.d0)
    !
  ELSEIF( DABS(DBLE(A_1)) > Ak_threshold ) THEN
    !A3(n) is zero => try to divide by A1
    msg = "Normalizing by A1(n)..."
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !In this case A1(n)=1
    A_kn(1,n) = DCMPLX(1.d0,0.d0)
    A_kn(2,n) = -1.d0*A_2/A_1
    A_kn(3,n) = A_3/A_1
    !
  ELSEIF( DABS(DBLE(A_2)) > Ak_threshold ) THEN
    !A3(n) and A1(n) are zero => try to divide by A2
    msg = "Normalizing by A2(n)..."
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !In this case A2(n)=1
    A_kn(1,n) = -1.d0*A_1/A_2
    A_kn(2,n) = DCMPLX(1.d0,0.d0)
    A_kn(3,n) = -1.d0*A_3/A_2
    !
  ELSE
    !All A_k(n) are zero:
    msg = "All Ak(n) are zero, dividing by a12..."
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !
    tempcmplx = DBLE(a_ik(1,2,n))**2 + AIMAG(a_ik(1,2,n))**2
    !
    WRITE(msg,*) "|a12| = ", a_ik(1,2,n)  !tempcmplx
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !
    IF( DBLE(tempcmplx)>1.d-15 ) THEN
      !We can normalize the A_k(n)
      !In that case A1(n)=1 and A3(n)=0
      A_kn(1,n) = DCMPLX(1.d0,0.d0)
      A_kn(3,n) = DCMPLX(0.d0,0.d0)
      !A2(n) = -(a11(n)/a12(n))*A1(n)
      A_kn(2,n) = -1.d0*A_kn(1,n)*a_ik(1,1,n)/a_ik(1,2,n)  !tempcmplx
    ELSE
      !Impossible to have proper Ak(n) => exit
      ifail = 2
      GOTO 400
    ENDIF
  ENDIF
ENDDO  !loop on n
!
!
IF(verbosity==4) THEN
  !Some debug messages
  msg = "Values of the A_k(n):"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO n=1,3
    DO k=1,3
      WRITE (msg,'(3X,a2,i1,a1,i1,a4,e10.3,a3,e10.3,a2)') "A_", k, "(", n, ") = ", &
            & DBLE(A_kn(k,n)), " + ", AIMAG(A_kn(k,n)), " i"
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
  ENDDO
ENDIF
!
!
300 CONTINUE
! 3. Determine the D(n) by solving the set of 6 linear equations:
!      Re{ SUM(n=1,3) +/-A_k(n).D(n) } = b_k          (13-88)
!      Re{ SUM(n=1,3) +/-B_i2k.A_k(n).D(n) } = 0      (13-89)
!    where b_k is the Burgers vector, and +/- signs are used when
!    the imaginary part of P(n) is positive and negative respectively.
!
!Save the sign of imaginary part of the P(n)
DO n=1,3
  IF( AIMAG(Pn(n))<0.d0 ) THEN
    PIm_sign(n) = -1.d0
  ELSE
    PIm_sign(n) = 1.d0
  ENDIF
ENDDO
!
!Compute the B_ijk(n) = c_ijk1 + c_ijk2*P(n)
B_ijk(:,:,:) = DCMPLX(0.d0,0.d0)
DO n=1,3
  DO i=1,3
    DO j=1,3
      DO k=1,3
        p = ELASTINDEX(i,j)
        q = ELASTINDEX(k,1)
        r = ELASTINDEX(k,2)
        B_ijk(p,k,n) = C_tensor(p,q) + C_tensor(p,r)*Pn(n)
      ENDDO
    ENDDO
  ENDDO
ENDDO
!
IF(verbosity==4) THEN
  !Some debug messages
  msg = "Values of the B_ijk(n):"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO n=1,3
    WRITE (msg,'(a4,i1)') "n = ", n
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO i=1,3
      DO j=1,3
        p = ELASTINDEX(i,j)
        DO k=1,3
          WRITE (msg,'(3X,a2,i1,i1,i1,a1,i1,a4,e10.3,a3,e10.3,a2)') "B_", i, j, k,     &
                & "(", n, ") = ", DBLE(B_ijk(p,k,n)), " + ", AIMAG(B_ijk(p,k,n)), " i"
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  msg = "Values of the +/- B_i2k(n)*A_k(n):"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO n=1,3
    DO i=1,3
      p = ELASTINDEX(i,2)    ! p = ij
      tempcmplx = DCMPLX(0.d0,0.d0)
      DO k=1,3
        tempcmplx = tempcmplx + B_ijk(p,k,n)*A_kn(k,n)*PIm_sign(n)
      ENDDO
      WRITE(msg,'(3X,a2,i1,a3,i1,a6,i1,a4,e10.3,a3,e10.3,a2)') "B_", i, &
           & "2k(", n, ")*A_k(", n, ") = ", DBLE(tempcmplx), " + ", -1.d0*AIMAG(tempcmplx), " i"
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
  ENDDO
ENDIF
!
!Write the parameters for the equations (13-88) and (13-89) in arrays:
! LHA = Left-Hand Array, contains the values of the left-hand side of the equations,
!       i.e.  LHA(1:3,1:3) contains the Re{SUM(n=1,3) +/-A_k(n)}
!             LHA(1:3,4:6) contains the Im{SUM(n=1,3) +/-A_k(n)}
!             LHA(4:6,1:3) contains the Re{SUM(n=1,3) +/-B_i2k(n)*A_k(n)}
!             LHA(4:6,4:6) contains the Im{SUM(n=1,3) +/-B_i2k(n)*A_k(n)}
! RHA = Right-Hand Array, contains the values of the right-hand side of the equations,
!       i.e.  RHA(1:3,1) contains the 3 coordinates of the Burgers vector
!             RHA(4:6,1) contains zeros
LHA(:,:) = 0.d0
!Set the +/-A_k(n)
DO n=1,3
  DO k=1,3
    LHA(k,n) = PIm_sign(n)*DBLE(A_kn(k,n))
    LHA(k,n+3) = -1.d0*PIm_sign(n)*AIMAG(A_kn(k,n))
  ENDDO
ENDDO
!Set the +/-B_i2k(n)*A_k(n)
DO n=1,3
  DO i=1,3
    p = ELASTINDEX(i,2)    ! p = ij
    tempcmplx = DCMPLX(0.d0,0.d0)
    DO k=1,3
      tempcmplx = tempcmplx + B_ijk(p,k,n)*A_kn(k,n)*PIm_sign(n)
    ENDDO
    LHA(i+3,n) = DBLE(tempcmplx)
    LHA(i+3,n+3) = -1.d0*AIMAG(tempcmplx)
  ENDDO
ENDDO
!Set the right-hand term
RHA(:,:) = 0.d0
DO i=1,3
  RHA(i,1) = b(i)
ENDDO
!
IF(verbosity==4) THEN
  !Some debug messages
  msg = "Terms for the equations:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO k=1,6
    WRITE (msg,'(6(e10.3,1X),a3,e10.3)') (LHA(k,n),n=1,6), " | ", RHA(k,1)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
  !
  !Some more debug messages
  msg = "Linear equations:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,3
    tempcmplx = DCMPLX( SUM(LHA(i,1:3)), SUM(LHA(i,4:6)) )
    WRITE (msg,'(3X,a5,e10.3,a3,e10.3,a7,i1,a6,e10.3)') &
          & "Re{ (", DBLE(tempcmplx), " + ", AIMAG(tempcmplx), "i) * D(", i, ") } = ", RHA(i,1)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
  DO i=4,6
    tempcmplx = DCMPLX( SUM(LHA(i,1:3)), SUM(LHA(i,4:6)) )
    WRITE (msg,'(3X,a5,e10.3,a3,e10.3,a7,i1,a6,e10.3)') &
          & "Re{ (", DBLE(tempcmplx), " + ", AIMAG(tempcmplx), "i) * D(", i-3, ") } = ", RHA(i,1)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
  msg = "Solving equations..."
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!
!Solving the 6 linear equations
!Note: DGESV is a LAPACK routine for solving linear equations
CALL DGESV(6,1,LHA(:,:),6,IPIV,RHA(:,:),6,k)
!
IF(verbosity==4) THEN
  !Some debug messages
  msg = "Values of RHA returned by DGESV:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO n=1,3
    WRITE (msg,'(3X,2e10.3)') RHA(n,1), RHA(n+3,1)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
ENDIF
!
!If it failed, k!=0
IF(k.NE.0) THEN
  ifail = 3
  GOTO 400
ENDIF
!If successful, k=0 and RHA now contains the complex solutions D(n) of the equations
!
!Save the real and imaginary parts of the solutions to Dn(n)
DO n=1,3
  Dn(n) = DCMPLX( RHA(n,1), RHA(n+3,1) )
ENDDO
!
IF(verbosity==4) THEN
  !Some debug messages
  msg = "Values of D(n):"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO n=1,3
    WRITE (msg,'(3X,a2,i1,a4,e10.3,a3,e10.3,a2)') "D(", n, ") = ", &
          & DBLE(Dn(n)), " + ", AIMAG(Dn(n)), " i"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
ENDIF
!
!
400 CONTINUE
! 4. At this point if ifail==0 then all the coefficients
!    A_k(n), D(n), P(n), and B_ijk(n) are known and are
!    provided as output of this routine.
!
!    If ifail is not zero then the routine failed.
!
!
END SUBROUTINE ANISO_COEFF
!
!
!********************************************************
! ANISO_STRESS
! This function calculates the stresses due to a
! dislocation in an anisotropic medium:
!   sigma_ij = Re{ (-1/2ipi) SUM(n=1,3) B_ijk(n)A_k(n)D(n)/(x1+P(n)x2) }
! The complex coefficients P(n), A_k(n), D(n) must
! be provided as input (see routine ANISO_COEFF above).
!********************************************************
!
SUBROUTINE ANISO_STRESS(P,a1,a2,a3,pos1,pos2,A_kn,B_ijk,Dn,Pn,sigma)
!
INTEGER:: i, j, k, r, n
INTEGER,INTENT(IN):: a1, a2, a3
REAL(dp),INTENT(IN):: pos1, pos2 !Position of the dislocation in the plane
                                 !orthogonal to the dislocation line
REAL(dp),DIMENSION(3),INTENT(IN):: P  !Atom position
REAL(dp),DIMENSION(3,3),INTENT(OUT):: sigma   !dislocation theoretical elastic stresses
COMPLEX(dp):: tempcmplx
COMPLEX(dp),PARAMETER:: frac2ipi = DCMPLX(0.d0,-1.d0/(2.d0*pi))  ! -1/(2*i*pi)
COMPLEX(dp),DIMENSION(3),INTENT(IN):: Dn, Pn    !anisotropy coefficients: D(n), P(n)...
COMPLEX(dp),DIMENSION(3,3),INTENT(IN):: A_kn    !... A_k(n) ...
COMPLEX(dp),DIMENSION(9,3,3),INTENT(IN):: B_ijk !... and B_ijk(n)
!
sigma(:,:) = 0.d0
!
DO i=1,3
  DO j=1,3
    r = ELASTINDEX(i,j)    ! r = ij
    tempcmplx = DCMPLX(0.d0,0.d0)
    !Compute the SUM(n=1,3)
    DO n=1,3
      DO k=1,3
        tempcmplx = tempcmplx + B_ijk(r,k,n)*A_kn(k,n)*Dn(n)/ &
                  & ( (P(a1)-pos1)+Pn(n)*(P(a2)-pos2) )
      ENDDO
    ENDDO
    !Multiply by -1/(2*i*pi)
    tempcmplx = tempcmplx * frac2ipi
    !Save the real part in sigma(i,j)
    sigma(i,j) = DBLE(tempcmplx)
  ENDDO
ENDDO

END SUBROUTINE ANISO_STRESS
!
!
!********************************************************
! ANISO_EFACTOR
! This function calculates the prelogarithmic energy
! factor for a dislocation in an anisotropic medium (13-83):
!   E = Kb²/4pi = (b_i/4pi) * Im{ SUM(n=1,3) B_i2k A_k(n) D(n) }
! The complex coefficients P(n), A_k(n), D(n) must
! be provided as input (see routine ANISO_COEFF above).
!********************************************************
!
FUNCTION ANISO_EFACTOR(b,A_kn,B_ijk,Dn,Pn) RESULT(Efactor)
!
IMPLICIT NONE
INTEGER:: i, k, r, n
REAL(dp),DIMENSION(3),INTENT(IN):: b  !Burgers vector
REAL(dp):: Efactor  ! E = Kb²/4pi
COMPLEX(dp):: tempcmplx
COMPLEX(dp),DIMENSION(3),INTENT(IN):: Dn, Pn    !anisotropy coefficients: D(n), P(n)...
COMPLEX(dp),DIMENSION(3,3),INTENT(IN):: A_kn    !... A_k(n) ...
COMPLEX(dp),DIMENSION(9,3,3),INTENT(IN):: B_ijk !... and B_ijk(n)
!
Efactor = 0.d0
tempcmplx = DCMPLX(0.d0,0.d0)
!
DO i=1,3
  r = ELASTINDEX(i,2)     ! j=2
  DO n=1,3
    DO k=1,3
      tempcmplx = tempcmplx + B_ijk(r,k,n)*A_kn(k,n)*Dn(n)
    ENDDO
  ENDDO
  Efactor = Efactor + b(i) * AIMAG(tempcmplx) / (4.d0*pi)
ENDDO
!
END FUNCTION ANISO_EFACTOR
!
!
!********************************************************
! DET_COMPMAT
! This function computes the determinant of a complex
! 2x2 matrix.
!********************************************************
FUNCTION DET_COMPMAT(M) RESULT(det)
!
IMPLICIT NONE
COMPLEX(dp):: det !determinant
COMPLEX(dp),DIMENSION(2,2):: M !complex 2x2 matrix
!
det = M(1,1)*M(2,2) - M(2,1)*M(1,2)
!
END FUNCTION DET_COMPMAT
      
!********************************************************
! DIAGMUL
! This function multiplies diagonal elements of a
! 9x3 matrix of real numbers. This operation is somewhat
! peculiar to what is done in subroutine ANISO_COEFF.
!********************************************************
FUNCTION DIAGMUL(M,i,j,k) RESULT(DS)
!
IMPLICIT NONE
INTEGER:: i, j, k
REAL(dp),DIMENSION(7):: DS  !results of multiplications
REAL(dp),DIMENSION(9,3):: M !the matrix
!
DS(:) = 0.d0
!
DS(7) = M(i,1)*M(j,1)*M(k,1)
!
DS(6) =  M(i,1)*M(j,1)*M(k,2) &
      & +M(i,1)*M(j,2)*M(k,1) &
      & +M(i,2)*M(j,1)*M(k,1)
!
DS(5) =  M(i,1)*M(j,1)*M(k,3) &
      & +M(i,1)*M(j,2)*M(k,2) &
      & +M(i,1)*M(j,3)*M(k,1) &
      & +M(i,2)*M(j,1)*M(k,2) &
      & +M(i,2)*M(j,2)*M(k,1) &
      & +M(i,3)*M(j,1)*M(k,1)
!
DS(4) =  M(i,1)*M(j,2)*M(k,3) &
      & +M(i,1)*M(j,3)*M(k,2) &
      & +M(i,2)*M(j,2)*M(k,2) &
      & +M(i,2)*M(j,1)*M(k,3) &
      & +M(i,2)*M(j,3)*M(k,1) &
      & +M(i,3)*M(j,1)*M(k,2) &
      & +M(i,3)*M(j,2)*M(k,1)
!
DS(3) =  M(i,1)*M(j,3)*M(k,3) &
      & +M(i,2)*M(j,2)*M(k,3) &
      & +M(i,2)*M(j,3)*M(k,2) &
      & +M(i,3)*M(j,1)*M(k,3) &
      & +M(i,3)*M(j,2)*M(k,2) &
      & +M(i,3)*M(j,3)*M(k,1)
!
DS(2) =  M(i,2)*M(j,3)*M(k,3) &
      & +M(i,3)*M(j,2)*M(k,3) &
      & +M(i,3)*M(j,3)*M(k,2)
!
DS(1) =  M(i,3)*M(j,3)*M(k,3)
!
RETURN
!
END FUNCTION DIAGMUL
!
!
END MODULE dislocation_aniso
