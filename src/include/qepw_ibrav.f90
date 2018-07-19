MODULE qepw_ibrav
!

!**********************************************************************************
!*  QEPW_IBRAV                                                                    *
!**********************************************************************************
!* This subroutine defines the cell vectors H(:,:) according to the Bravais       *
!* lattice index (0<ibrav<14) and the cell parameters celldm(:), defined as:      *
!*      celldm(1) = a               celldm(4) = cos(gamma)                        *
!*      celldm(2) = b/a             celldm(5) = cos(beta)                         *
!*      celldm(3) = c/a             celldm(6) = cos(alpha)                        *
!* where a, b and c are the lengths of the three cell vectors, alpha the angle    *
!* between the second and third cell vectors, i.e. alpha=angle(H(2,:),H(3,:)),    *
!* beta=angle(H(1,:),H(3,:)) and gamma=angle(H(1,:),H(2,:)).                      *
!* This subroutine is particular to Quantum Espresso files, and is called by      *
!* "input/in_qe_pw.f90" and by "modes/1ia_qe_pw.f90".                             *
!* The relations used here can be found in the documentation of QE:               *
!*    http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html        *
!**********************************************************************************
!* (C) January 2013 - Pierre Hirel                                                *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 26 Apr. 2013                                     *
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
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE QE_PW_IBRAV(ibrav,celldm,H)
!
CHARACTER(LEN=128):: msg
INTEGER,INTENT(IN):: ibrav  !index of Bravais lattice
REAL(dp):: tx, ty, tz, u, v
REAL(dp),DIMENSION(6),INTENT(IN):: celldm
REAL(dp),DIMENSION(3,3),INTENT(OUT):: H   !Base vectors of the supercell
!
!
!
msg = "ENTERING QE_PW_IBRAV"
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,*) "ibrav =", ibrav
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,'(a8,6f8.4)') "celldm =", celldm(1:6)
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
H(:,:) = 0.d0
!
SELECT CASE(ibrav)
CASE(0,1)
  !cubic P (simple cubic)
  !v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,1)
  H(1,1) = celldm(1)
  H(2,2) = celldm(1)
  H(3,3) = celldm(1)
CASE(2)
  !cubic F (fcc)
  !v1 = (a/2)(-1,0,1),  v2 = (a/2)(0,1,1), v3 = (a/2)(-1,1,0)
  H(1,1) = -0.5d0*celldm(1)
  H(1,3) = 0.5d0*celldm(1)
  H(2,2) = 0.5d0*celldm(1)
  H(2,3) = 0.5d0*celldm(1)
  H(3,1) = -0.5d0*celldm(1)
  H(3,2) = 0.5d0*celldm(1)
CASE(3)
  !cubic I (bcc)
  !v1 = (a/2)(1,1,1),  v2 = (a/2)(-1,1,1),  v3 = (a/2)(-1,-1,1)
  H(1,:) = 0.5d0*celldm(1)
  H(2,:) = 0.5d0*celldm(1)
  H(2,1) = -0.5d0*celldm(1)
  H(3,:) = -0.5d0*celldm(1)
  H(3,3) = 0.5d0*celldm(1)
CASE(4)
  !Hexagonal and Trigonal P        celldm(3)=c/a
  !v1 = a(1,0,0),  v2 = a(-1/2,sqrt(3)/2,0),  v3 = a(0,0,c/a)
  H(1,1) = celldm(1)
  H(2,1) = -0.5d0*celldm(1)
  H(2,2) = 0.5d0*DSQRT(3.d0)*celldm(1)
  H(3,3) = celldm(1)*celldm(3)
CASE(5,-5)
  tx = DSQRT(0.5d0*(1.d0-celldm(4)))
  ty = DSQRT((1.d0-celldm(4))/6.d0)
  tz = DSQRT((1.d0+2.d0*celldm(4))/3.d0)
  IF(ibrav==5) THEN
    !Trigonal R, 3fold axis c        celldm(4)=cos(alpha)
    !The crystallographic vectors form a three-fold star around
    !the z-axis, the primitive cell is a simple rhombohedron:
    !v1 = a(tx,-ty,tz),   v2 = a(0,2ty,tz),   v3 = a(-tx,-ty,tz)
    !where c=cos(alpha) is the cosine of the angle alpha between
    !any pair of crystallographic vectors, tx, ty, tz are:
    !  tx=sqrt((1-c)/2), ty=sqrt((1-c)/6), tz=sqrt((1+2c)/3)
    H(1,1) = celldm(1)*tx
    H(1,2) = -1.d0*celldm(1)*ty
    H(1,3) = celldm(1)*tz
    H(2,2) = celldm(1)*2.d0*ty
    H(2,3) = celldm(1)*tz
    H(3,1) = -1.d0*celldm(1)*tx
    H(3,2) = -1.d0*celldm(1)*ty
    H(3,3) = celldm(1)*tz
  ELSE  !ibrav=-5
    !Trigonal R, 3fold axis <111>    celldm(4)=cos(alpha)
    !The crystallographic vectors form a three-fold star around
    !<111>. Defining a' = a/sqrt(3) :
    !v1 = a' (u,v,v),   v2 = a' (v,u,v),   v3 = a' (v,v,u)
    !where u and v are defined as
    !  u = tz - 2*sqrt(2)*ty,  v = tz + sqrt(2)*ty
    !and tx, ty, tz as for case ibrav=5
    u = tz - 2.d0*DSQRT(2.d0)*ty
    v = tz + DSQRT(2.d0)*ty
    H(1,1) = (celldm(1)/DSQRT(3.d0))*u
    H(1,2) = (celldm(1)/DSQRT(3.d0))*v
    H(1,3) = (celldm(1)/DSQRT(3.d0))*v
    H(2,1) = (celldm(1)/DSQRT(3.d0))*v
    H(2,2) = (celldm(1)/DSQRT(3.d0))*u
    H(2,3) = (celldm(1)/DSQRT(3.d0))*v
    H(3,1) = (celldm(1)/DSQRT(3.d0))*v
    H(3,2) = (celldm(1)/DSQRT(3.d0))*v
    H(3,3) = (celldm(1)/DSQRT(3.d0))*u
  ENDIF
CASE(6)
  !Tetragonal P (st)               celldm(3)=c/a
  !v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,c/a)
  H(1,1) = celldm(1)
  H(2,2) = celldm(1)
  H(3,3) = celldm(3)
CASE(7)
  !Tetragonal I (bct)              celldm(3)=c/a
  !v1=(a/2)(1,-1,c/a),  v2=(a/2)(1,1,c/a),  v3=(a/2)(-1,-1,c/a)
  H(1,1) = 0.5d0*celldm(1)
  H(1,2) = -0.5d0*celldm(1)
  H(1,3) = 0.5d0*celldm(1)*celldm(3)
  H(2,1) = 0.5d0*celldm(1)
  H(2,2) = 0.5d0*celldm(1)
  H(2,3) = 0.5d0*celldm(1)*celldm(3)
  H(3,1) = -0.5d0*celldm(1)
  H(3,2) = -0.5d0*celldm(1)
  H(3,3) = 0.5d0*celldm(1)*celldm(3)
CASE(8)
  !Orthorhombic P                  celldm(2)=b/a
  !                                celldm(3)=c/a
  !v1 = (a,0,0),  v2 = (0,b,0), v3 = (0,0,c)
  H(1,1) = celldm(1)
  H(2,2) = celldm(1)*celldm(2)
  H(3,3) = celldm(1)*celldm(3)
CASE(9)
  !Orthorhombic base-centered(bco) celldm(2)=b/a
  !                                celldm(3)=c/a
  !v1 = (a/2, b/2,0),  v2 = (-a/2,b/2,0),  v3 = (0,0,c)
  H(1,1) = 0.5d0*celldm(1)
  H(1,2) = 0.5d0*celldm(1)*celldm(2)
  H(2,1) = -0.5d0*celldm(1)
  H(2,2) = 0.5d0*celldm(1)*celldm(2)
  H(3,3) = celldm(1)*celldm(3)
CASE(-9)
  !as 9, alternate description
  !v1 = (a/2,-b/2,0),  v2 = (a/2,-b/2,0),  v3 = (0,0,c)
  H(1,1) = 0.5d0*celldm(1)
  H(1,2) = -0.5d0*celldm(1)*celldm(2)
  H(2,1) = 0.5d0*celldm(1)
  H(2,2) = -0.5d0*celldm(1)*celldm(2)
  H(3,3) = celldm(1)*celldm(3)
CASE(10)
  !Orthorhombic face-centered      celldm(2)=b/a
  !                                celldm(3)=c/a
  !v1 = (a/2,0,c/2),  v2 = (a/2,b/2,0),  v3 = (0,b/2,c/2)
  H(1,1) = 0.5d0*celldm(1)
  H(1,3) = 0.5d0*celldm(1)*celldm(3)
  H(2,1) = 0.5d0*celldm(1)
  H(2,2) = 0.5d0*celldm(1)*celldm(2)
  H(3,2) = 0.5d0*celldm(1)*celldm(2)
  H(3,3) = 0.5d0*celldm(1)*celldm(3)
CASE(11)
  !Orthorhombic body-centered      celldm(2)=b/a
  !                                celldm(3)=c/a
  !v1=(a/2,b/2,c/2),  v2=(-a/2,b/2,c/2),  v3=(-a/2,-b/2,c/2)
  H(1,1) = 0.5d0*celldm(1)
  H(1,2) = 0.5d0*celldm(1)*celldm(2)
  H(1,3) = 0.5d0*celldm(1)*celldm(3)
  H(2,1) = -0.5d0*celldm(1)
  H(2,2) = 0.5d0*celldm(1)*celldm(2)
  H(2,3) = 0.5d0*celldm(1)*celldm(3)
  H(3,1) = -0.5d0*celldm(1)
  H(3,2) = -0.5d0*celldm(1)*celldm(2)
  H(3,3) = 0.5d0*celldm(1)*celldm(3)
CASE(12)
  !Monoclinic P, unique axis c     celldm(2)=b/a
  !                                celldm(3)=c/a,
  !                                celldm(4)=cos(ab)=cos(gamma)
  !v1=(a,0,0), v2=(b*cos(gamma),b*sin(gamma),0),  v3 = (0,0,c)
  !where gamma is the angle between axis a and b.
  H(1,1) = celldm(1)
  H(2,1) = celldm(1)*celldm(2)*celldm(4)
  H(2,2) = celldm(1)*celldm(2)*DSIN(DACOS(celldm(4)))
  H(3,3) = celldm(1)*celldm(3)
CASE(-12)
  !Monoclinic P, unique axis b     celldm(2)=b/a
  !                                celldm(3)=c/a,
  !                                celldm(5)=cos(ac)=cos(beta)
  !v1 = (a,0,0), v2 = (0,b,0), v3 = (c*sin(beta),0,c*cos(beta))
  !where beta is the angle between axis a and c
  H(1,1) = celldm(1)
  H(2,2) = celldm(1)*celldm(2)
  H(3,1) = celldm(1)*celldm(3)*DSIN(DACOS(celldm(5)))
  H(3,3) = celldm(1)*celldm(3)*celldm(5)
CASE(13)
  !Monoclinic base-centered        celldm(2)=b/a
  !                                celldm(3)=c/a,
  !                                celldm(4)=cos(ab)=cos(gamma)
  !v1 = (  a/2,         0,                -c/2),
  !v2 = (b*cos(gamma), b*sin(gamma), 0),
  !v3 = (  a/2,         0,                  c/2),
  !where gamma is the angle between axis a and b
  H(1,1) = 0.5d0*celldm(1)
  H(1,3) = -0.5d0*celldm(1)*celldm(3)
  H(2,1) = celldm(1)*celldm(2)*celldm(4)
  H(2,2) = celldm(1)*celldm(2)*DSIN(DACOS(celldm(4)))
  H(3,1) = 0.5d0*celldm(1)
  H(3,3) = 0.5d0*celldm(1)*celldm(3)
CASE(14)
  !Triclinic                       celldm(2)= b/a,
  !                                celldm(3)= c/a,
  !                                celldm(4)= cos(bc)=cos(alpha)
  !                                celldm(5)= cos(ac)=cos(beta)
  !                                celldm(6)= cos(ab)=cos(gamma)
  !v1 = (a, 0, 0),
  !v2 = (b*cos(gamma), b*sin(gamma), 0)
  !v3 = (c*cos(beta),  c*(cos(alpha)-cos(beta)cos(gamma))/sin(gamma),
  !      c*sqrt( 1 + 2*cos(alpha)cos(beta)cos(gamma)
  !              - cos(alpha)^2-cos(beta)^2-cos(gamma)^2 )/sin(gamma) )
  !where alpha is the angle between axis b and c
  !      beta is the angle between axis a and c
  !      gamma is the angle between axis a and b
  H(1,1) = celldm(1)
  H(2,1) = celldm(1)*celldm(2)*celldm(6)
  H(2,2) = celldm(1)*celldm(2)*DSIN(DACOS(celldm(6)))
  H(3,1) = celldm(1)*celldm(3)*celldm(5)
  H(3,2) = celldm(1)*celldm(3)*( celldm(4) - celldm(5)*celldm(6)/DSIN(DACOS(celldm(6))) )
  H(3,3) = celldm(1)*celldm(3)*DSQRT( 1.d0 + 2.d0*celldm(4)*celldm(5)*celldm(6)       &
          &                            - celldm(4)**2 - celldm(5)**2 - celldm(6)**2    &
          &                          )/DSIN(DACOS(celldm(6)))
CASE DEFAULT
  !error
END SELECT
!
!
WRITE(msg,'(a8,3f8.4)') "H(1,:) =", H(1,:)
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,'(a8,3f8.4)') "H(2,:) =", H(2,:)
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,'(a8,3f8.4)') "H(3,:) =", H(3,:)
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
END SUBROUTINE QE_PW_IBRAV
!
!
END MODULE