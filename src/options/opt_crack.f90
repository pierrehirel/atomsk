MODULE crack
!
!**********************************************************************************
!*  CRACK                                                                         *
!**********************************************************************************
!* This module reads atomic positions from an array P and introduces a            *
!* straight crack. The crack tip lies along the direction <crackline>,            *
!* at the position (pos1,pos2) in the plane normal to <crackline>.                *
!* The crack can be mode I (opening), mode II (in-plane shear) or                 *
!* mode III (out-of-plane shear).                                                 *
!* The theoretical (continuum) stresses associated with the crack are computed    *
!* and saved in the array AUX.                                                    *
!**********************************************************************************
!* (C) March 2012 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 16 April 2024                                    *
!**********************************************************************************
!* List of subroutines in this module:                                            *
!* CRACK_XYZ           main subroutine - introduces a crack                       *
!* DISPMODEI           applies isotropic disp. of a mode I crack                  *
!* DISPMODEII          applies isotropic disp. of a mode II crack                 *
!* DISPMODEIII         applies isotropic disp. of a mode III crack                *
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
SUBROUTINE CRACK_XYZ(H,P,S,crackmode,cracktype,crackK,crackline,crackplane,mu,nu,pos1,pos2,SELECT,AUXNAMES,AUX)
!
!
IMPLICIT NONE
!Input variables
CHARACTER(LEN=1),INTENT(IN):: crackline !direction of crack tip line, must be x, y or z
CHARACTER(LEN=1),INTENT(IN):: crackplane !normal to the plane of cut, must be x, y or z
CHARACTER(LEN=3),INTENT(IN):: crackmode !crack mode, must be I, II or III
CHARACTER(LEN=6),INTENT(IN):: cracktype !"stress" or "strain"
REAL(dp),INTENT(IN):: crackK  !stress intensity factor (GPa.A^1/2)
REAL(dp),INTENT(IN):: mu   !shear modulus of the material (GPa)
REAL(dp),INTENT(IN):: nu   !Poisson ratio of the material
!
!Internal variables
CHARACTER(LEN=1):: dir1, dir2, dir3 !directions x, y, z
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES, newAUXNAMES !names of auxiliary prop.
LOGICAL:: doshells !apply displacements also to shells?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: a1, a2, a3 !a3=1,2,3 if crackline=x,y,z respectively
INTEGER:: i
INTEGER:: sig1, sig2, sig3, sig4, sig5, sig6 !for indexing stresses
REAL(dp):: pos1, pos2 !Coordinates of crack in the plane normal to crackline
REAL(dp):: kappa !plane stress=(3-nu)/(1+nu); plane strain=3-4nu
REAL(dp):: R  !distance of atom to the crack
REAL(dp):: theta  !angle between direction a1 and the (crak tip->atom) segment
REAL(dp),DIMENSION(3):: V1, V2 !temporary vectors
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: sigma   !crack theoretical elastic stresses
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S  !Positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Q, T                !Positions of atoms, shells (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX            !auxiliary properties (temporary)
!
!Initialize variables
a1 = 0
a2 = 0
a3 = 0
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(T)) DEALLOCATE(T)
IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
IF(ALLOCATED(newAUXNAMES)) DEALLOCATE(newAUXNAMES)
IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
  doshells=.TRUE.
ELSE
  doshells=.FALSE.
ENDIF
!
!
msg = 'Entering CRACK_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!
100 CONTINUE
!Define the axes: a3=crack line
IF(crackline=='x' .OR. crackline=='X') THEN
  a3 = 1
  dir3 = "x"
ELSEIF(crackline=='y' .OR. crackline=='Y') THEN
  a3 = 2
  dir3 = "y"
ELSEIF(crackline=='z' .OR. crackline=='Z') THEN
  a3 = 3
  dir3 = "z"
ELSE
  CALL ATOMSK_MSG(2800,(/crackline/),(/0.d0/))
  nerr = nerr+1
  GOTO 1000
ENDIF
!a2=normal to plane of cut
!a1=last component
IF(crackplane=='x' .OR. crackplane=='X') THEN
  a2 = 1
  dir2 = "x"
  IF(a3==2) THEN
    a1 = 3
    dir1 = "z"
  ELSEIF(a3==3) THEN
    a1 = 2
    dir1 = "y"
  ELSE
    CALL ATOMSK_MSG(2811,(/""/),(/0.d0/))
    nerr = nerr+1
    GOTO 1000
  ENDIF
ELSEIF(crackplane=='y' .OR. crackplane=='Y') THEN
  a2 = 2
  dir2 = "y"
  IF(a3==1) THEN
    a1 = 3
    dir1 = "z"
  ELSEIF(a3==3) THEN
    a1 = 1
    dir1 = "x"
  ELSE
    CALL ATOMSK_MSG(2811,(/""/),(/0.d0/))
    nerr = nerr+1
    GOTO 1000
  ENDIF
ELSEIF(crackplane=='z' .OR. crackplane=='Z') THEN
  a2 = 3
  dir2 = "z"
  IF(a3==1) THEN
    a1 = 2
    dir1 = "y"
  ELSEIF(a3==2) THEN
    a1 = 1
    dir1 = "x"
  ELSE
    CALL ATOMSK_MSG(2811,(/""/),(/0.d0/))
    nerr = nerr+1
    GOTO 1000
  ENDIF
ELSE
  CALL ATOMSK_MSG(2800,(/crackplane/),(/0.d0/))
  nerr = nerr+1
  GOTO 1000
ENDIF
!
WRITE(msg,*) "a1, a2, a3: ", a1, a2, a3
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Check that the values make sense
IF(pos1<0.d0) pos1 = DABS(pos1)
IF(pos2<0.d0) pos2 = DABS(pos2)
!If pos1 and pos2 are reduced, convert them to cartesian
!IF(pos1>0.d0 .AND. pos1<=1.d0 .AND. pos2>0.d0 .AND. pos2<=1.d0) THEN
!  pos1 = pos1*( SUM(H(:,a1)) )
!  pos2 = pos2*( SUM(H(:,a2)) )
!ENDIF
!
!
!Print messages
CALL ATOMSK_MSG(2108,(/crackmode,crackline//'  '/),(/pos1,pos2/))
!
!
!If the stress intensity factor crackK is zero then skip the whole thing
IF(crackK<=0.d0) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2747,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Check that the directions normal to crack line are not too small
IF( VECLENGTH(H(:,a1))<10.d0 .OR. VECLENGTH(H(:,a2))<10.d0 ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2748,(/''/),(/0.d0/))
ENDIF
!
!Set kappa
IF( cracktype=="stress" ) THEN
  !Plane stress
  kappa = (3.d0-nu) / (1.d0+nu)
ELSE
  !Plane strain
  kappa = 3.d0 - 4.d0*nu
ENDIF
!
!
!
200 CONTINUE
!Introduce displacements due to the crack
IF(crackmode=="I") THEN
  DO i=1,SIZE(P,1)
    IF(IS_SELECTED(SELECT,i)) THEN
      CALL DISPMODEI(i,P(i,:),a1,a2,pos1,pos2,crackK,mu,kappa)
    ENDIF
  ENDDO
  IF(ALLOCATED(S)) THEN
    DO i=1,SIZE(S,1)
      IF(IS_SELECTED(SELECT,i)) THEN
        CALL DISPMODEI(i,S(i,:),a1,a2,pos1,pos2,crackK,mu,kappa)
      ENDIF
    ENDDO
  ENDIF
  !Elongate the box along a2
  H(a2,a2) = H(a2,a2) + crackK/mu * DSQRT(pos2/(2.d0*pi)) * (kappa-1.d0)
  !
ELSEIF(crackmode=="II") THEN
  DO i=1,SIZE(P,1)
    IF(IS_SELECTED(SELECT,i)) THEN
      CALL DISPMODEII(i,P(i,:),a1,a2,pos1,pos2,crackK,mu,kappa)
    ENDIF
  ENDDO
  IF(ALLOCATED(S)) THEN
    DO i=1,SIZE(S,1)
      IF(IS_SELECTED(SELECT,i)) THEN
        CALL DISPMODEII(i,S(i,:),a1,a2,pos1,pos2,crackK,mu,kappa)
      ENDIF
    ENDDO
  ENDIF
  !
ELSEIF(crackmode=="III") THEN
  DO i=1,SIZE(P,1)
    IF(IS_SELECTED(SELECT,i)) THEN
      CALL DISPMODEIII(i,P(i,:),a1,a2,a3,pos1,pos2,crackK,mu)
    ENDIF
  ENDDO
  IF(ALLOCATED(S)) THEN
    DO i=1,SIZE(S,1)
      IF(IS_SELECTED(SELECT,i)) THEN
        CALL DISPMODEIII(i,S(i,:),a1,a2,a3,pos1,pos2,crackK,mu)
      ENDIF
    ENDDO
  ENDIF
  !
ENDIF
!
!Crack was successfully created
CALL ATOMSK_MSG(2109,(/''/),(/0.d0/))
!
!
!
400 CONTINUE
!Compute the theoretical (continuum) crack stresses.
!Formulae are given below.
!Only the Voigt components of the stresses are output
!as auxiliary properties, namely:
! sigma_xx sigma_yy sigma_zz sigma_xy sigma_yz sigma_xz
!If these auxiliary properties already exist, the present
!stresses must be added to the existing ones.
!
sigma(:,:) = 0.d0
!Setup array AUX for crack stresses
sig1=0
sig2=0
sig3=0
sig4=0
sig5=0
sig6=0
IF(.NOT.ALLOCATED(AUX)) THEN
  !No auxiliary property exist => create array
  ALLOCATE(AUXNAMES(6))
  ALLOCATE( AUX( SIZE(P(:,1)),6 ) )
  AUX(:,:) = 0.d0
  sig1=1
  sig2=2
  sig3=3
  sig4=4
  sig5=5
  sig6=6
  !Set names of auxiliary properties
  AUXNAMES(sig1) = "sigma_"//dir1//dir1
  AUXNAMES(sig2) = "sigma_"//dir2//dir2
  AUXNAMES(sig3) = "sigma_"//dir3//dir3
  AUXNAMES(sig4) = "sigma_"//dir2//dir3
  AUXNAMES(sig5) = "sigma_"//dir1//dir3
  AUXNAMES(sig6) = "sigma_"//dir1//dir2
ELSE
  !Auxiliary property already exist
  !Check if it already contains the sigma
  DO i=1,SIZE(AUXNAMES)
    IF( AUXNAMES(i)=="sigma_"//dir1//dir1 ) THEN
      sig1=i
    ELSEIF( AUXNAMES(i)=="sigma_"//dir2//dir2 ) THEN
      sig2=i
    ELSEIF( AUXNAMES(i)=="sigma_"//dir3//dir3) THEN
      sig3=i
    ELSEIF( AUXNAMES(i)=="sigma_"//dir2//dir3 .OR.      &
          & AUXNAMES(i)=="sigma_"//dir3//dir2      ) THEN
      sig4=i
    ELSEIF( AUXNAMES(i)=="sigma_"//dir1//dir3 .OR.      &
          & AUXNAMES(i)=="sigma_"//dir3//dir1      ) THEN
      sig5=i
    ELSEIF( AUXNAMES(i)=="sigma_"//dir1//dir2 .OR.      &
          & AUXNAMES(i)=="sigma_"//dir2//dir1      ) THEN
      sig6=i
    ENDIF
  ENDDO
  !If the sigma do not exist, extend array AUX
  IF(MIN(sig1,sig2,sig3,sig4,sig5,sig6)==0) THEN
    sig1=SIZE(AUXNAMES)+1
    sig2=sig1+1
    sig3=sig2+1
    sig4=sig3+1
    sig5=sig4+1
    sig6=sig5+1
    !Create extended arrays
    ALLOCATE(newAUXNAMES( SIZE(AUXNAMES)+6 ) )
    ALLOCATE(newAUX( SIZE(P(:,1)), SIZE(AUXNAMES)+6 ) )
    newAUX(:,:) = 0.d0
    DO i=1,SIZE(AUXNAMES)
      newAUXNAMES(i) = AUXNAMES(i)
    ENDDO
    DO i=1,SIZE(P(:,1))
      newAUX(i,1:sig1-1) = AUX(i,:)
    ENDDO
    !Replace old arrays with new ones
    DEALLOCATE(AUXNAMES)
    DEALLOCATE(AUX)
    ALLOCATE(AUXNAMES( SIZE(newAUXNAMES) ))
    ALLOCATE(AUX( SIZE(newAUX(:,1)), SIZE(newAUX(1,:)) ))
    AUXNAMES = newAUXNAMES
    AUX = newAUX
    DEALLOCATE(newAUXNAMES)
    DEALLOCATE(newAUX)
    !Set names of auxiliary properties
    AUXNAMES(sig1) = "sigma_"//dir1//dir1
    AUXNAMES(sig2) = "sigma_"//dir2//dir2
    AUXNAMES(sig3) = "sigma_"//dir3//dir3
    AUXNAMES(sig4) = "sigma_"//dir2//dir3
    AUXNAMES(sig5) = "sigma_"//dir1//dir3
    AUXNAMES(sig6) = "sigma_"//dir1//dir2
  ENDIF
ENDIF
!
!Compute stresses for each atom and save them to AUX
IF(crackmode=="I") THEN
  !sigma_11 = K/(sqrt(2piR)) * cos(theta/2) * [1 - sin(theta/2)*sin(3theta/2)]
  !sigma_22 = K/(sqrt(2piR)) * cos(theta/2) * [1 + sin(theta/2)*sin(3theta/2)]
  !sigma_12 = K/(sqrt(2piR)) * cos(theta/2) * sin(theta/2) * cos(3theta/2)
  !sigma_33 = sigma13 = sigma23 = 0
  DO i=1,SIZE(P(:,1))
    V1(:) = 0.d0
    V1(a1) = P(i,a1)-pos1
    V1(a2) = P(i,a2)-pos2
    V2(:) = 0.d0
    V2(a1) = 1.d0
    R = VECLENGTH( V1 )
    theta = ANGVEC( V1,V2 )
    !theta is not signed => correct for that
    IF( V1(a2)<0.d0 ) THEN
      theta = -1.d0*theta
    ENDIF
    !avoid division by zero
    IF(DABS(R)>1.d-6) THEN
      sigma(1,1) = crackK/(DSQRT(2.d0*pi*R))*DCOS(theta/2.d0)*                &
                 & ( 1.d0 - DSIN(theta/2.d0)*DSIN(3.d0*theta/2.d0))
      sigma(2,2) = crackK/(DSQRT(2.d0*pi*R))*DCOS(theta/2.d0)*                &
                 & ( 1.d0 + DSIN(theta/2.d0)*DSIN(3.d0*theta/2.d0))
      sigma(1,2) = crackK/(DSQRT(2.d0*pi*R))*DCOS(theta/2.d0)*                &
                 & DSIN(theta/2.d0)*DCOS(3.d0*theta/2.d0)
      AUX(i,sig1) = AUX(i,sig1) + sigma(1,1)
      AUX(i,sig2) = AUX(i,sig2) + sigma(2,2)
      AUX(i,sig6) = AUX(i,sig6) + sigma(1,2)
      IF( cracktype=="strain" ) THEN
        AUX(i,sig3) = AUX(i,sig3) + nu*sigma(1,1)*sigma(2,2)
      ENDIF
    ENDIF
  ENDDO
  !
ELSEIF(crackmode=="II") THEN
  !sigma_11 = -K/(sqrt(2piR)) * sin(theta/2) * [2 + cos(theta/2)*cos(3theta/2)]
  !sigma_22 = K/(sqrt(2piR)) * sin(theta/2) * cos(theta/2) * cos(3theta/2)
  !sigma_12 = K/(sqrt(2piR)) * cos(theta/2) * [1 - sin(theta/2)*sin(3theta/2)]
  !sigma_33 = sigma13 = sigma23 = 0
  DO i=1,SIZE(P(:,1))
    V1(:) = 0.d0
    V1(a1) = P(i,a1)-pos1
    V1(a2) = P(i,a2)-pos2
    V2(:) = 0.d0
    V2(a1) = 1.d0
    R = VECLENGTH( V1 )
    theta = ANGVEC( V1,V2 )
    !theta is not signed => correct for that
    IF( V1(a2)<0.d0 ) THEN
      theta = -1.d0*theta
    ENDIF
    !avoid division by zero
    IF(DABS(R)>1.d-6) THEN
      sigma(1,1) = -1.d0*crackK/(DSQRT(2.d0*pi*R))*DSIN(theta/2.d0)*          &
                 & ( 2.d0 + DCOS(theta/2.d0)*DCOS(3.d0*theta/2.d0))
      sigma(2,2) = crackK/(DSQRT(2.d0*pi*R))*DSIN(theta/2.d0)*                &
                 & DCOS(theta/2.d0)*DCOS(3.d0*theta/2.d0)
      sigma(1,2) = crackK/(DSQRT(2.d0*pi*R))*DCOS(theta/2.d0)*                &
                 & ( 1.d0 - DSIN(theta/2.d0)*DSIN(3.d0*theta/2.d0))
      AUX(i,sig1) = AUX(i,sig1) + sigma(1,1)
      AUX(i,sig2) = AUX(i,sig2) + sigma(2,2)
      AUX(i,sig6) = AUX(i,sig6) + sigma(1,2)
      IF( cracktype=="strain" ) THEN
        AUX(i,sig3) = AUX(i,sig3) + nu*sigma(1,1)*sigma(2,2)
      ENDIF
    ENDIF
  ENDDO
  !
ELSEIF(crackmode=="III") THEN
  !sigma_23 = K/(sqrt(2piR)) * cos(theta/2)
  !sigma_13 = -K/(sqrt(2piR)) * sin(theta/2)
  !sigma_11 = sigma_22 = sigma_33 = 0
  DO i=1,SIZE(P(:,1))
    V1(:) = 0.d0
    V1(a1) = P(i,a1)-pos1
    V1(a2) = P(i,a2)-pos2
    V2(:) = 0.d0
    V2(a1) = 1.d0
    R = VECLENGTH( V1 )
    theta = ANGVEC( V1,V2 )
    !theta is not signed => correct for that
    IF( V1(a2)<0.d0 ) THEN
      theta = -1.d0*theta
    ENDIF
    !avoid division by zero
    IF(DABS(R)>1.d-6) THEN
      sigma(2,3) = crackK/(DSQRT(2.d0*pi*R))*DCOS(theta/2.d0)
      sigma(1,3) = -1.d0*crackK/(DSQRT(2.d0*pi*R))*DSIN(theta/2.d0)
      AUX(i,sig4) = AUX(i,sig4) + sigma(2,3)
      AUX(i,sig5) = AUX(i,sig5) + sigma(1,3)
    ENDIF
  ENDDO
!
ENDIF
!
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
IF(ALLOCATED(Q)) DEALLOCATE(Q)
!
!
END SUBROUTINE CRACK_XYZ
!
!
!
!********************************************************
! DISPMODEI
! This subroutine calculates the displacements due to
! a mode I crack in an isotropic medium:
! ux = K/2mu * sqrt(R/2pi) * cos(theta/2) * [ kappa - 1 + 2sin²(theta/2) ]
! uy = K/2mu * sqrt(R/2pi) * sin(theta/2) * [ kappa + 1 - 2cos²(theta/2) ]
! Formulae can be found for instance in L.B. Freund,
! "Dynamic fracture mechanics".
!********************************************************
!
SUBROUTINE DISPMODEI(i,P,a1,a2,pos1,pos2,crackK,mu,kappa)
!
IMPLICIT NONE
INTEGER,INTENT(IN):: a1, a2
INTEGER:: i
REAL(dp),INTENT(IN):: crackK, mu, kappa
REAL(dp),INTENT(IN):: pos1, pos2  !Position of the crack in the plane
                                !orthogonal to the crack tip
REAL(dp):: prev1, prev2
REAL(dp):: R  !distance of atom to the crack
REAL(dp):: theta  !angle between direction a1 and the (crak tip->atom) segment
REAL(dp),DIMENSION(3):: V1, V2
REAL(dp),DIMENSION(3),INTENT(INOUT):: P  !Atom position
!
prev1 = P(a1)
prev2 = P(a2)
!
V1(:) = 0.d0
V1(a1) = P(a1)-pos1
V1(a2) = P(a2)-pos2
V2(:) = 0.d0
V2(a1) = 1.d0
R = VECLENGTH( V1 )
theta = ANGVEC( V1,V2 )
!theta is not signed => correct for that
IF( V1(a2)<0.d0 ) THEN
  theta = -1.d0*theta
ENDIF
!
P(a1) = P(a1) + crackK/(2.d0*mu) * DSQRT(R/(2.d0*pi)) * DCOS(theta/2.d0) *  &
      &         ( kappa - 1.d0 + 2.d0*((DSIN(theta/2.d0))**2) )
P(a2) = P(a2) + crackK/(2.d0*mu) * DSQRT(R/(2.d0*pi)) * DSIN(theta/2.d0) *  &
      &         ( kappa + 1.d0 - 2.d0*((DCOS(theta/2.d0))**2) )
!Message if displacement was too large
IF( DABS(prev1-P(a1)) >= 10.d0 .OR. DABS(prev2-P(a2)) >= 10.d0 ) &
                    & CALL ATOMSK_MSG(2727,(/''/),(/DBLE(i)/))
!
END SUBROUTINE DISPMODEI
!
!
!
!********************************************************
! DISPMODEII
! This subroutine calculates the displacements due to
! a mode II crack in an isotropic medium:
! ux = K/2mu * sqrt(R/2pi) * sin(theta/2) * [ kappa + 1 + 2cos²(theta/2) ]
! uy = -K/2mu * sqrt(R/2pi) * cos(theta/2) * [ kappa - 1 - 2sin²(theta/2) ]
! Formulae can be found for instance in L.B. Freund,
! "Dynamic fracture mechanics".
!********************************************************
!
SUBROUTINE DISPMODEII(i,P,a1,a2,pos1,pos2,crackK,mu,kappa)
!
IMPLICIT NONE
INTEGER,INTENT(IN):: a1, a2
INTEGER:: i
REAL(dp),INTENT(IN):: crackK, mu, kappa
REAL(dp),INTENT(IN):: pos1, pos2  !Position of the crack in the plane
                                !orthogonal to the crack tip
REAL(dp):: prev1, prev2
REAL(dp):: R  !distance of atom to the crack
REAL(dp):: theta  !angle between direction a1 and the (crak tip->atom) segment
REAL(dp),DIMENSION(3):: V1, V2
REAL(dp),DIMENSION(3),INTENT(INOUT):: P  !Atom position
!
prev1 = P(a1)
prev2 = P(a2)
!
V1(:) = 0.d0
V1(a1) = P(a1)-pos1
V1(a2) = P(a2)-pos2
V2(:) = 0.d0
V2(a1) = 1.d0
R = VECLENGTH( V1 )
theta = ANGVEC( V1,V2 )
!theta is not signed => correct for that
IF( V1(a2)<0.d0 ) THEN
  theta = -1.d0*theta
ENDIF
!
P(a1) = P(a1) + crackK/(2.d0*mu) * DSQRT(R/(2.d0*pi)) * DSIN(theta/2.d0) *  &
      &         ( kappa + 1.d0 + 2.d0*((DCOS(theta/2.d0))**2) )
P(a2) = P(a2) - crackK/(2.d0*mu) * DSQRT(R/(2.d0*pi)) * DCOS(theta/2.d0) *  &
      &         ( kappa - 1.d0 - 2.d0*((DSIN(theta/2.d0))**2) )
!Message if displacement was too large
IF( DABS(prev1-P(a1)) >= 10.d0 .OR. DABS(prev2-P(a2)) >= 10.d0 ) &
                    & CALL ATOMSK_MSG(2727,(/''/),(/DBLE(i)/))
!
END SUBROUTINE DISPMODEII
!
!
!
!********************************************************
! DISPMODEIII
! This subroutine calculates the displacements due to
! a mode III crack in an isotropic medium:
! uz = K/mu * sqrt(R/2pi) * sin(theta/2)
! Formulae can be found for instance in L.B. Freund,
! "Dynamic fracture mechanics".
!********************************************************
!
SUBROUTINE DISPMODEIII(i,P,a1,a2,a3,pos1,pos2,crackK,mu)
!
IMPLICIT NONE
INTEGER,INTENT(IN):: a1, a2, a3
INTEGER:: i
REAL(dp),INTENT(IN):: crackK, mu
REAL(dp),INTENT(IN):: pos1, pos2  !Position of the crack in the plane
                                !orthogonal to the crack tip
REAL(dp):: prev
REAL(dp):: R  !distance of atom to the crack
REAL(dp):: theta  !angle between direction a1 and the (crak tip->atom) segment
REAL(dp),DIMENSION(3):: V1, V2
REAL(dp),DIMENSION(3),INTENT(INOUT):: P  !Atom position
!
prev = P(a3)
!
V1(:) = 0.d0
V1(a1) = P(a1)-pos1
V1(a2) = P(a2)-pos2
V2(:) = 0.d0
V2(a1) = 1.d0
R = VECLENGTH( V1 )
theta = ANGVEC( V1,V2 )
!theta is not signed => correct for that
IF( V1(a2)<0.d0 ) THEN
  theta = -1.d0*theta
ENDIF
!
P(a3) = P(a3) + crackK/(2.d0*mu) * DSQRT(R/(2.d0*pi)) * DSIN(theta/2.d0)
!Message if displacement was too large
IF( DABS(prev-P(a3)) >= 10.d0 ) CALL ATOMSK_MSG(2727,(/''/),(/DBLE(i)/))
!
END SUBROUTINE DISPMODEIII
!
!
END MODULE crack
