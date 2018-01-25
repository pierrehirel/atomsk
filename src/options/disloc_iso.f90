MODULE dislocation_iso
!
!**********************************************************************************
!*  DISLOCATION_ISO                                                               *
!**********************************************************************************
!* This module contains functions computing the displacements due to              *
!* a straight dislocation in the framework of isotropic elasticity.               *
!* Formulae can be found for instance in:                                         *
!* J.P. Hirth, J. Lothe, 'Theory of dislocations', 1st ed. (1968), p.75.          *
!* These functions are used by the module "opt_dislocation.f90".                  *
!**********************************************************************************
!* (C) January 2018 - Pierre Hirel                                                *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 25 Jan. 2018                                     *
!**********************************************************************************
!* List of functions in this module:                                              *
!* DISPSCREW           applies isotropic disp. of a screw disloc. to 1 atom       *
!* DISPEDGE            applies isotropic disp. of an edge disloc. to 1 atom       *
!* List of subroutines in this module:                                            *
!* STRESSSCREW         computes stresses due to screw disloc.                     *
!* STRESSEDGE          computes stresses due to edge disloc.                      *
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
!********************************************************
! DISPSCREW
! This function calculates the displacements due to
! a screw dislocation in an isotropic medium.
! Formulae can be found for instance in
! J.P. Hirth, J. Lothe, 'Theory of dislocations',
! 1st ed. (1968), p.59; or in
! D.Hull,D.J.Bacon,'Introduction to dislocations',
! 4th Ed.(2001),p.66.
!********************************************************
!
FUNCTION DISPSCREW(i,P,a1,a2,a3,b,pos1,pos2) RESULT(disp)
!
IMPLICIT NONE
INTEGER,INTENT(IN):: a1, a2, a3
INTEGER,INTENT(IN):: i            !index of atom to be displaced
REAL(dp),INTENT(IN):: b           !norm of Burgers vector
REAL(dp),INTENT(IN):: pos1, pos2  !Position of the dislocation in the plane
                                  !orthogonal to the dislocation line
REAL(dp),DIMENSION(3),INTENT(IN):: P      !Atom position
REAL(dp),DIMENSION(3):: disp  !displacement of the atom
!
disp(:)  = 0.d0
disp(a3) = b/(2.d0*pi)*DATAN2( P(a2)-pos2 , P(a1)-pos1 )
!Message if displacement was too large
IF( VECLENGTH(disp(:)) >= 2.d0*DABS(b) ) CALL ATOMSK_MSG(2727,(/''/),(/DBLE(i)/))
!
END FUNCTION DISPSCREW
!
!
!********************************************************
! DISPEDGE
! This function calculates the displacements due to
! an edge dislocation in an isotropic medium.
!********************************************************
!
FUNCTION DISPEDGE(i,P,a1,a2,b,nu,pos1,pos2) RESULT(disp)
!
IMPLICIT NONE
INTEGER,INTENT(IN):: a1, a2
INTEGER,INTENT(IN):: i           !index of atom to be displaced
REAL(dp),INTENT(IN):: b          !norm of Burgers vector
REAL(dp),INTENT(IN):: nu         !Poisson ratio
REAL(dp),INTENT(IN):: pos1, pos2 !Position of the dislocation in the plane
                                 !orthogonal to the dislocation line
REAL(dp),DIMENSION(3),INTENT(IN):: P  !Atom position
REAL(dp),DIMENSION(3):: disp !edge displacement of the atom
!
disp(:) = 0.d0
disp(a1) = b/(2.d0*pi)*                                     &
          &  (DATAN( (P(a2)-pos2)/(P(a1)-pos1) ) +          &
          &   (P(a1)-pos1)*(P(a2)-pos2)/( 2.d0*(1.d0-nu)*   &
          &   ( (P(a1)-pos1)**2+(P(a2)-pos2)**2) )          &
          &  )
disp(a2) = 0.d0 - b/(2.d0*pi)*                              &
          &  ( (1.d0-2.d0*nu)/(4.d0*(1.d0-nu))*             &
          &   DLOG((P(a1)-pos1)**2+(P(a2)-pos2)**2) +       &
          &   ( (P(a1)-pos1)**2-(P(a2)-pos2)**2 )/          &
          &   ( 4.d0*(1.d0-nu)*                             &
          &     ( (P(a1)-pos1)**2+(P(a2)-pos2)**2 )         &
          &   )                                             &
          &  )
!Message if displacement was too large
IF( VECLENGTH(disp(:)) >= 2.d0*DABS(b) .OR. DABS(disp(a2)) >= 2.d0*DABS(b)) &
  & CALL ATOMSK_MSG(2727,(/''/),(/DBLE(i)/))
!
END FUNCTION DISPEDGE
!
!
!********************************************************
! STRESSSCREW
! This function calculates the displacements due to
! a screw dislocation in an isotropic medium:
!   sigma_11 = sigma_22 = sigma_33 = sigma_12 = 0
!   sigma_13 = (-mu*b/2pi) * y/(x²+y²)
!   sigma_23 = (mu*b/2pi) * x/(x²+y²)
! NOTE: the shear modulus mu is not passed to this
! routine, therefore what is resturned is the stress
! tensor normalized by mu, i.e. sigma/mu.
!********************************************************
!
SUBROUTINE STRESSSCREW(P,a1,a2,b,nu,pos1,pos2,sigma)
!
IMPLICIT NONE
INTEGER,INTENT(IN):: a1, a2
REAL(dp),INTENT(IN):: b          !norm of Burgers vector
REAL(dp),INTENT(IN):: nu         !Poisson ratio
REAL(dp),INTENT(IN):: pos1, pos2 !Position of the dislocation in the plane
                                 !orthogonal to the dislocation line
REAL(dp):: tempreal
REAL(dp),DIMENSION(3),INTENT(IN):: P  !Atom position
REAL(dp),DIMENSION(3,3),INTENT(OUT):: sigma   !dislocation theoretical elastic stresses
!
sigma(:,:) = 0.d0
tempreal = (P(a1)-pos1)**2 + (P(a2)-pos2)**2
!avoid division by zero
IF(DABS(tempreal)>1.d-6) THEN
  sigma(2,3) = -1.d0*b/(2.d0*pi) * (P(a2)-pos2)/tempreal
  sigma(1,3) = b/(2.d0*pi) * (P(a1)-pos1)/tempreal
ENDIF
!
END SUBROUTINE STRESSSCREW
!
!
!********************************************************
! STRESSEDGE
! This function calculates the displacements due to
! an edge dislocation in an isotropic medium:
!   sigma_11 = (-mu*b/2pi(1-nu)) * y*(3x²+y²)/(x²+y²)²
!   sigma_22 = (mu*b/2pi(1-nu)) * y*(x²-y²)/(x²+y²)²
!   sigma_12 = (mu*b/2pi(1-nu)) * x*(x²-y²)/(x²+y²)²
!   sigma_33 = nu*(sigma_11 + sigma_22)
!   sigma_13 = sigma_23 = 0
! NOTE: the shear modulus mu is not passed to this
! routine, therefore what is returned is the stress
! tensor normalized by mu, i.e. sigma/mu.
!********************************************************
!
SUBROUTINE STRESSEDGE(P,a1,a2,b,nu,pos1,pos2,sigma)
!
IMPLICIT NONE
INTEGER,INTENT(IN):: a1, a2
REAL(dp),INTENT(IN):: b          !norm of Burgers vector
REAL(dp),INTENT(IN):: nu         !Poisson ratio
REAL(dp),INTENT(IN):: pos1, pos2 !Position of the dislocation in the plane
                                 !orthogonal to the dislocation line
REAL(dp):: tempreal
REAL(dp),DIMENSION(3),INTENT(IN):: P  !Atom position
REAL(dp),DIMENSION(3,3),INTENT(OUT):: sigma   !dislocation theoretical elastic stresses
!
sigma(:,:) = 0.d0
tempreal = (P(a1)-pos1)**2 + (P(a2)-pos2)**2
!avoid division by zero
IF(DABS(tempreal)>1.d-6) THEN
  sigma(1,1) = -1.d0*b/(2.d0*pi*(1.d0-nu)) *   &
          & (P(a2)-pos2)*(3.d0*(P(a1)-pos1)**2+(P(a2)-pos2)**2) &
          & / (tempreal**2)
  sigma(2,2) = b/(2.d0*pi*(1.d0-nu)) *         &
            & (P(a2)-pos2)*((P(a1)-pos1)**2-(P(a2)-pos2)**2)    &
            & / (tempreal**2)
  sigma(1,2) = b/(2.d0*pi*(1.d0-nu)) *         &
              & (P(a1)-pos1)*((P(a1)-pos1)**2-(P(a2)-pos2)**2)  &
              & / (tempreal**2)
  sigma(3,3) = nu*( sigma(1,1) + sigma(2,2) )
ENDIF
!
END SUBROUTINE STRESSEDGE
!
!
!********************************************************
! ISO_EFACTOR
! This function calculates the prelogarithmic energy
! factor for a dislocation in an isotropic medium (3-87):
!   E = (µb²/4pi) * ( cos²i + sin²i/(1-nu) )
! where i is the angle between the Burgers vector b
! and the dislocation line: screw i=0, edge i=90°.
! NOTE: the shear modulus mu is not passed to this
! routine, therefore what is returned is the energy
! factor normalized by mu, i.e. E/mu
!********************************************************
!
FUNCTION ISO_EFACTOR(b,nu,angle) RESULT(Efactor)
!
IMPLICIT NONE
REAL(dp),INTENT(IN):: b  !Burgers vector
REAL(dp),INTENT(IN):: nu !Poisson ratio
REAL(dp),INTENT(IN):: angle !angle between b and disloc.line (radians)
REAL(dp):: Efactor
!  
Efactor = ( (DCOS(angle))**2 + ((DSIN(angle))**2)/(1.d0-nu)) * b*b / (4.d0*pi)
!
END FUNCTION ISO_EFACTOR
!
!
!
END MODULE dislocation_iso