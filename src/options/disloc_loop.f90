MODULE dislocation_loop
!
!**********************************************************************************
!*  DISLOCATION_LOOP                                                              *
!**********************************************************************************
!* This module contains functions computing the displacements due to              *
!* a dislocation loop in the framework of isotropic elasticity.                   *
!* These functions were copied from the program Babel developped by               *
!* Emmanuel Clouet, and adapted for use within Atomsk.                            *
!* Babel is available from the following Website:                                 *
!*   http://emmanuel.clouet.free.fr/Programs/Babel/                               *
!* These functions are used by the module "opt_disloc.f90".                       *
!**********************************************************************************
!* (C) October 2017 - Emmanuel Clouet                                             *
!*     Service de Recherches de Métallurgie Physique                              *
!*     SRMP CEA-Saclay, 91191 Gif-sur-Yvette, France                              *
!*     emmanuel.clouet@cea.fr                                                     *
!* Last modification: P. Hirel - 09 June 2022                                     *
!**********************************************************************************
!* List of subroutines in this module:                                            *
!* LOOP_SEGMENTS         builds a list of points that form a circular loop        *
!* LOOP_DISPLACEMENT     compute displacement due to a disloc.loop at a point     *
!* SOLIDANGLE            compute sol.angle used for disp. due to disloc.triangle  *
!* DISLOSEG_DISPLACEMENT_ISO compute disp. due to disloc.segment                  *
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
USE functions
USE subroutines
!
!
CONTAINS
!
!
!********************************************************
! LOOP_SEGMENTS
! Provided the center, radius and normal of a circular loop,
! this routine discretizes the loop into points.
!********************************************************
FUNCTION LOOP_SEGMENTS(center,radius,normal) RESULT(xLoop)
!
IMPLICIT NONE
CHARACTER(LEN=1), INTENT(IN):: normal ! Direction normal to the loop: must be X, Y or Z
REAL(dp),INTENT(IN):: radius   !radius of the dislocation loop
REAL(dp),DIMENSION(3),INTENT(IN):: center
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: xLoop !coordinates of points forming the loop
!
CHARACTER(LEN=128):: msg
INTEGER:: a1, a2, a3  ! indices replacing X and Y
INTEGER:: i
INTEGER:: Npoints  ! number of points
REAL(dp):: perimeter  !perimeter of the loop
REAL(dp):: angle   ! angle for current point
REAL(dp):: theta   ! angle between two points (radians)
!
IF(ALLOCATED(xLoop)) DEALLOCATE(xLoop)
!
! Define directions according to normal
SELECT CASE(normal)
CASE('x','X')
  !Loop in (y,z) plane
  a1 = 2
  a2 = 3
  a3 = 1
CASE('y','Y')
  !Loop in (x,z) plane
  a1 = 3
  a2 = 1
  a3 = 2
CASE DEFAULT
  !Loop in (x,y) plane
  a1 = 1
  a2 = 2
  a3 = 3
END SELECT
!
! Loop perimeter
perimeter = 2.d0*pi*radius
!
! Define number of points forming the loop
! The following criteria are used as a trade-off between
! good accuracy and computational efficiency:
! - each dislocation segment should have a length of 5 angströms;
! - the loop should contain at least 3 points
!    (for very small loops, this will result in segments shorter than 5 A);
! - there should not be more than 100 points
!    (for very large loops, this will result in segments longer than 5 A).
Npoints = MAX( 3 , MIN( NINT(perimeter/5.d0) , 100 ) )
!
! Angle between two consecutive points
theta = 2.d0*pi / DBLE(Npoints)
!
! Allocate xLoop
ALLOCATE( xLoop(Npoints,3) )
xLoop(:,:) = 0.d0
!
! Save the position of each point of the loop
angle = 0.d0
DO i=1,SIZE(xLoop,1)
  xLoop(i,a1) = center(a1) + radius*DCOS(angle)
  xLoop(i,a2) = center(a2) + radius*DSIN(angle)
  xLoop(i,a3) = center(a3)
  ! Increment angle for next point
  angle = angle + theta
ENDDO
!
IF( verbosity==4 ) THEN
  WRITE(msg,'(i3,a42)') Npoints, "  POINTS FOR DISLOCATION LOOP"
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  DO i=1,SIZE(xLoop,1)
    WRITE(msg,'(4X,3f9.3)') xLoop(i,1:3)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  ENDDO
ENDIF
!
END FUNCTION LOOP_SEGMENTS
!
!
!********************************************************
! LOOP_DISPLACEMENT
! Calculate displacement u(:) created by a loop
! at the position R
!********************************************************
FUNCTION LOOP_DISPLACEMENT(R, b, nu, center, xLoop) RESULT(u)
! 
IMPLICIT NONE
REAL(dp),DIMENSION(3),INTENT(IN):: b  ! Burgers vector of the loop
REAL(dp),DIMENSION(3),INTENT(IN):: center  ! center of the loop
REAL(dp),DIMENSION(3),INTENT(IN):: R !atom coordinates
REAL(dp),INTENT(IN):: nu ! Poisson coefficient
REAL(dp),DIMENSION(:,:),INTENT(IN):: xLoop !coordinates of points forming the loop
!
INTEGER :: i
REAL(dp):: omega
REAL(dp),DIMENSION(1:3):: xA, xB, xC
REAL(dp),DIMENSION(1:3):: u !final displacement due to loop

u(:) = 0.d0
omega = 0.d0
xC(:) = center(:) - R(:)
!
! Loop on all segments forming the dislocation loop
DO i=1,SIZE(xLoop,1)
  !Coordinates of point A
  IF( i==1 ) THEN
    xA(:) = xLoop(SIZE(xLoop,1),:) - R(:)
  ELSE
    xA(:) = xLoop(i-1,:) - R(:)
  ENDIF
  ! Coordinates of point B
  xB(:) = xLoop(i,:) - R(:)
  !
  ! Part due to solid angle
  omega = omega + SolidAngle(xA, xB, xC)
  ! Part due to elasticity
  u(:) = u(:) + DisloSeg_displacement_iso(xA, xB, b(:), nu)
ENDDO

! Final total displacement
u(:) = u(:) + b(:)*omega

END FUNCTION LOOP_DISPLACEMENT
!
!
!
!********************************************************
! SOLIDANGLE
! Calculate solid angle (normalized by 4*pi) used to obtain the 
! displacement created by a dislocation triangle loop ABC at the origin
! (field point = origin)
! Ref.: Van Oosterom, A. and Strackee, J., 
!       The Solid Angle of a Plane Triangle, 
!       IEEE Transactions on Biomedical Engineering BME-30, 125 (1983).
!********************************************************
FUNCTION SolidAngle(xA, xB, xC) RESULT(Omega)
!
IMPLICIT NONE
!
! Extremities of the triangle loop
REAL(dp),DIMENSION(3),INTENT(IN) :: xA, xB, xC
!
! Solid angle (normalized by 4*pi)
REAL(dp):: omega
REAL(dp),PARAMETER:: factor=1.d0/(2.d0*pi)
!
REAL(dp) :: rA, rB, rC, numerator, denominator
!
rA = VECLENGTH(xA)
rB = VECLENGTH(xB)
rC = VECLENGTH(xC)
!
numerator = SCALAR_TRIPLE_PRODUCT( xA, xB, xC )
denominator = rA*rB*rC + DOT_PRODUCT( xA, xB )*rC &
            & + DOT_PRODUCT( xB, xC )*rA + DOT_PRODUCT( xC, xA )*rB
!
omega = factor*ATAN2( numerator, denominator )
!
END FUNCTION SolidAngle
!
!
!********************************************************
! DISLOSEG_DISPLACEMENT_ISO
! Calculate displacement created by dislocation segment AB
! once the solid angle part has been removed
! Isotropic elastic calculation with nu Poisson coef.
! Ref.: Eq. (1) in Barnett, D. M. 
!       The Displacement Field of a Triangular Dislocation Loop
!       Philos. Mag. A, 1985, 51, 383-387 
!********************************************************
FUNCTION DisloSeg_displacement_iso(xA, xB, b, nu) RESULT(u)
!
IMPLICIT NONE
REAL(dp),DIMENSION(3),INTENT(IN):: xA, xB ! Extremities of the segment
REAL(dp),DIMENSION(3),INTENT(IN):: b  ! Burgers vector
REAL(dp),INTENT(IN):: nu    ! Poisson coefficient
REAL(dp),DIMENSION(3):: u ! Displacement
REAL(dp):: rA, rB
REAL(dp),DIMENSION(1:3)::  tAB, nAB 
!
rA = VECLENGTH(xA)
rB = VECLENGTH(xB)
!
! Tangent vector
tAB(:) = xB(:) - xA(:)
tAB(:) = tAB(:)/VECLENGTH(tAB)
!
! Normal vector
nAB(:) = CROSS_PRODUCT(xA,xB)
nAB(:) = nAB(:)/VECLENGTH(nAB)
!
u(:) = ( -(1.d0-2.d0*nu)*CROSS_PRODUCT(b, tAB)* &
     &  LOG( (rB + DOT_PRODUCT(xB,tAB))/(rA + DOT_PRODUCT(xA,tAB)) ) &
     &  + DOT_PRODUCT(b,nAB)*CROSS_PRODUCT(xB/rB-xA/rA,nAB) ) &
     &  /(8.d0*pi*(1.d0-nu))
!
END FUNCTION DisloSeg_displacement_iso
!
!
!
END MODULE dislocation_loop
