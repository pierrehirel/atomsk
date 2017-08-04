MODULE dislocation
!
!**********************************************************************************
!*  DISLOCATION                                                                   *
!**********************************************************************************
!* This module reads atomic positions from an array P and introduces a            *
!* straight dislocation of Burgers vector "b", with line along "dislocline".      *
!* The position of the dislocation is defined by (pos1,pos2) in the plane normal  *
!* to "dislocline".                                                               *
!* If "C_tensor" only contains zeros, isotropic elasticity is used                *
!* (i.e. subroutines DISPEDGE or DISPSCREW).                                      *
!* If a non-zero elastic tensor "C_tensor" is provided, then the anisotropic      *
!* elasticity is used, by first calculating the coefficients (ANISO_COEFF)        *
!* and then by applying the corresponding displacements (ANISO_DISP).             *
!* Note that the elastic tensor "C_tensor" is assumed to contain the correct      *
!* elastic tensor for the current orientation of the system, i.e. the C'          *
!* and NOT the constants C which correspond to the orientation [100][010][001].   *
!* In other words the elastic tensor must have been correctly oriented            *
!* BEFORE being passed to the present routine.                                    *
!**********************************************************************************
!* (C) May 2010 - Pierre Hirel                                                    *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 15 Feb. 2017                                     *
!**********************************************************************************
!* List of subroutines in this module:                                            *
!* DISLOC_XYZ          main subroutine - introduces a dislocation                 *
!* DISPEDGE            applies isotropic disp. of an edge disloc. to 1 atom       *
!* DISPSCREW           applies isotropic disp. of a screw disloc. to 1 atom       *
!* ANISO_DISP          applies anisotropic disp. of a disloc. to 1 atom           *
!* ANISO_COEFF         computes coefficients for disp. in anisotropic medium      *
!* List of functions in this module:                                              *
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
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE DISLOC_XYZ(H,P,S,disloctype,dislocline,dislocplane,b,nu,pos1,pos2,SELECT,AUXNAMES,AUX,C_tensor)
!
!
IMPLICIT NONE
!Input variables
CHARACTER(LEN=1),INTENT(IN):: dislocline !direction of dislocation line, must be x, y or z
CHARACTER(LEN=1),INTENT(IN):: dislocplane !normal to plane of cut, must be x, y or z
CHARACTER(LEN=5),INTENT(IN):: disloctype !type of dislocation, must be screw, edge or edge2
REAL(dp),INTENT(IN):: nu   !Poissons ratio of the material
                         !(ignored if disloctype=screw, or if aniso=.TRUE.)
REAL(dp),DIMENSION(9,9),INTENT(IN):: C_tensor   !Elastic tensor
REAL(dp),DIMENSION(3),INTENT(IN):: b    !Burgers vector (in Angstroms)
!
!Internal variables
CHARACTER(LEN=1):: dir1, dir2, dir3 !directions x, y, z
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES, newAUXNAMES !names of auxiliary properties
LOGICAL:: aniso !shall anisotropic elasticity be used?
LOGICAL:: doshells !apply displacements also to shells?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: a1, a2, a3 !a3=1,2,3 if dislocline=x,y,z respectively
INTEGER:: i, j, k, n, r, u
INTEGER:: sig1, sig2, sig3, sig4, sig5, sig6 !for indexing stresses
REAL(dp):: bsign !sign of Burgers vector (+1.d0 or -1.d0)
REAL(dp):: pos1, pos2 !Coordinates of dislocation in the plane
REAL(dp):: tempreal
REAL(dp),DIMENSION(3):: disp     !elastic displacement applied to an atom
REAL(dp),DIMENSION(3):: Kfactor  !energy factor (anisotropic elasticity)
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: rot_matrix !rotation matrix
REAL(dp),DIMENSION(3,3):: sigma      !dislocation theoretical elastic stresses
REAL(dp),DIMENSION(9,9):: C_tensor_rot !rotated tensor if dislocline is not along Z
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S   !Positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Q, T                 !Positions of atoms, shells (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX    !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX               !auxiliary properties (temporary)
COMPLEX(dp):: tempcmplx
COMPLEX(dp),PARAMETER:: frac2ipi = DCMPLX(0.d0,-1.d0/(2.d0*pi))  ! -1/(2*i*pi)
COMPLEX(dp),DIMENSION(3):: Dn, Pn    !anisotropy coefficients D(n), P(n)-
COMPLEX(dp),DIMENSION(3,3):: A_kn    !-, A_k(n) -
COMPLEX(dp),DIMENSION(9,3,3):: B_ijk !-, and B_ijk(n)
!
!Initialize variables
a1 = 0
a2 = 0
a3 = 0
k = 0
bsign = 1.d0
aniso = .FALSE.
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
msg = 'Entering DISLOC_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!
100 CONTINUE
!Define the axes: a3=dislocation line
IF(dislocline=='x' .OR. dislocline=='X') THEN
  a3 = 1
  dir3 = "x"
ELSEIF(dislocline=='y' .OR. dislocline=='Y') THEN
  a3 = 2
  dir3 = "y"
ELSEIF(dislocline=='z' .OR. dislocline=='Z') THEN
  a3 = 3
  dir3 = "z"
ELSE
  CALL ATOMSK_MSG(2800,(/dislocline/),(/0.d0/))
  nerr = nerr+1
  GOTO 1000
ENDIF
!
!a2=normal to plane of cut
!a1=last component
IF(dislocplane=='x' .OR. dislocplane=='X') THEN
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
ELSEIF(dislocplane=='y' .OR. dislocplane=='Y') THEN
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
ELSEIF(dislocplane=='z' .OR. dislocplane=='Z') THEN
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
  CALL ATOMSK_MSG(2800,(/dislocplane/),(/0.d0/))
  nerr = nerr+1
  GOTO 1000
ENDIF
!
WRITE(msg,*) "a1, a2, a3: ", a1, a2, a3
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Check if elastic tensor is defined
IF( ANY( C_tensor(1:3,1:3).NE.0.d0 ) ) THEN
  aniso = .TRUE.
  k = 1 !only used by message below
ELSE
  aniso = .FALSE.
  k = 0 !only used by message below
ENDIF
!
!
!Print messages
CALL ATOMSK_MSG(2061,(/disloctype,dislocline//'    '/),(/b(1),b(2),b(3),DBLE(k),pos1,pos2/))
!
!
!If Burgers vector is zero then skip the whole thing
IF(VECLENGTH(b)==0.d0) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2725,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!Check that the directions normal to dislocation line are not too small
IF( VECLENGTH(H(:,a1))<3.d0*VECLENGTH(b(:)) .OR. VECLENGTH(H(:,a2))<3.d0*VECLENGTH(b(:)) ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2726,(/''/),(/0.d0/))
ENDIF
!Check that the position for the dislocation is inside the box
IF( pos1<MINVAL(P(:,a1)) .OR. pos1>MAXVAL(P(:,a1)) .OR. &
  & pos2<MINVAL(P(:,a2)) .OR. pos2>MAXVAL(P(:,a2))     ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2754,(/'dislocation'/),(/0.d0/))
ENDIF
!
!
!
200 CONTINUE
!If anisotropic elasticity is used, solve the equations
!WARNING: at this point it is assumed that the elastic tensor corresponds
!       to the current crystallographic orientation of the system,
!       i.e. that C_tensor contains the C' constants and not the C.
IF(aniso) THEN
  CALL ATOMSK_MSG(2062,(/''/),(/0.d0/))
  !
  !If dislocation line is not along Z then 
  !we have to provide a rotated C_tensor
  IF(a3.NE.3) THEN
    !Set rotation matrix:
    rot_matrix(:,:) = 0.d0
    IF( a3==1 ) THEN
      !dislocline along X => rotate by -90° around Y
      rot_matrix(2,2) = 1.d0
      rot_matrix(1,3) = 1.d0
      rot_matrix(3,1) = -1.d0
    ELSE  !i.e. if a3==2
      !dislocline along Y => rotate by 90° around X
      rot_matrix(1,1) = 1.d0
      rot_matrix(2,3) = 1.d0
      rot_matrix(3,2) = -1.d0
    ENDIF
    !Rotate elastic tensor
    C_tensor_rot = ROTELAST( C_tensor, rot_matrix )
    !
    IF( a2.NE.a3+1 ) THEN
      !Dislocline is along X, but slip plane is not normal to Y, OR
      !dislocline is along Y, but slip plane is not normal to Z
      !=> rotate elastic tensor around the Z axis
      rot_matrix(:,:) = 0.d0
      rot_matrix(1,2) = -1.d0
      rot_matrix(2,1) = 1.d0
      rot_matrix(3,3) = 1.d0
      !Rotate elastic tensor
      C_tensor_rot = ROTELAST( C_tensor, rot_matrix )
    ENDIF
    !
  ELSE
    !Otherwise (i.e. dislocline is along Z) just use C_tensor
    C_tensor_rot = C_tensor
  ENDIF
  !
  !Compute the coefficients A_k(n), D(n) and P(n) for anisotropic displacements
  CALL ANISO_COEFF((/b(a1),b(a2),b(a3)/),C_tensor_rot,A_kn,Dn,Pn,B_ijk,k)
  !
  !If the routine ANISO_COEFF succeeded then k=0
  !Otherwise display an error message and exit
  IF(k.NE.0) THEN
    nerr=nerr+1
    CALL ATOMSK_MSG(2807,(/''/),(/DBLE(k)/))
    GOTO 1000
  ELSE
    !Routine ANISO_COEFF worked: now the coefficients A_kn, Dn, Pn, B_ijk
    !are known, for a dislocation line along Z.
    !If dislocation line is *not* along Z, then 
    !rotate the C_tensor back to its previous orientation
    IF(a3.NE.3) THEN
      !Set rotation matrix:
      rot_matrix(:,:) = 0.d0
      IF( a3==1 ) THEN
        !dislocline along X => rotate by 90° around a1=Y
        rot_matrix(2,2) = 1.d0
        rot_matrix(1,3) = -1.d0
        rot_matrix(3,1) = 1.d0
      ELSE  !i.e. if a3==2
        !dislocline along Y=a3 => rotate by -90° around a2=X
        rot_matrix(1,1) = 1.d0
        rot_matrix(2,3) = -1.d0
        rot_matrix(3,2) = 1.d0
      ENDIF
      !Rotate elastic tensor
      C_tensor_rot = ROTELAST( C_tensor, rot_matrix )
      !
      !In addition, we have to fix the A_kn and B_ijk
      !so that the displacements u_k are applied in the
      !correct direction
      IF( a3==1 ) THEN
        !dislocline along X => swap Z and X
        DO i=1,3
          tempcmplx = A_kn(1,i)
          A_kn(1,i) = A_kn(3,i)
          A_kn(3,i) = tempcmplx
          DO n=1,3
            tempcmplx = B_ijk(i,1,n)
            B_ijk(i,1,n) = B_ijk(i,3,n)
            B_ijk(i,3,n) = tempcmplx
          ENDDO
        ENDDO
      ELSE  !i.e. if a3==2
        !dislocline along Y => swap Z and Y
        DO i=1,3
          tempcmplx = A_kn(2,i)
          A_kn(2,i) = A_kn(3,i)
          A_kn(3,i) = tempcmplx
          DO n=1,3
            tempcmplx = B_ijk(i,2,n)
            B_ijk(i,2,n) = B_ijk(i,3,n)
            B_ijk(i,3,n) = tempcmplx
          ENDDO
        ENDDO
      ENDIF
    ENDIF
    CALL ATOMSK_MSG(2096,(/''/),(/0.d0/))
  ENDIF
ENDIF
!
!
!
300 CONTINUE
!Apply atomic displacements
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  SCREW DISLOCATION  !!!!!!!!!!
IF(disloctype=='screw') THEN
  bsign = b(a3)/DABS(b(a3))
  !
  IF( aniso ) THEN
    !Anisotropic elasticity
    DO i=1,SIZE(P,1)
      IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
        disp = ANISO_DISP(i,P(i,:),a1,a2,a3,pos1,pos2,A_kn,Dn,Pn)
        P(i,1:3) = P(i,1:3) + disp(:)
        !
        !Same if shells exist
        IF( doshells ) THEN
          S(i,1:3) = S(i,1:3) + disp(:)
        ENDIF
      ENDIF
    ENDDO
  ELSE
    !Isotropic elasticity
    DO i=1,SIZE(P,1)
      IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
        disp = DISPSCREW(i,P(i,:),a1,a2,a3,b(a3),pos1,pos2)
        P(i,1:3) = P(i,1:3) + disp(:)
        !
        !Same if shells exist
        IF( doshells ) THEN
          S(i,1:3) = S(i,1:3) + disp(:)
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  EDGE DISLOCATION  !!!!!!!!!!
! Method 1: insert a plane of atoms
! (NP is increased so we have to use a new array Q)
ELSEIF(disloctype=='edge') THEN
  bsign = b(a1)/DABS(b(a1))
  !First, insert a plane of atoms by copying an existing plane
  !Count k = number of atoms that will be inserted
  k = 0
  DO i=1,SIZE(P,1)
    IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
      IF( P(i,a1)>=pos1+DABS(b(a1))            .AND.                  &
        & P(i,a1)<pos1+2.d0*DABS(b(a1))-0.01d0 .AND.                  &
        & (P(i,a2)-pos2)*b(a1)/DABS(b(a1))>=0.d0            ) THEN
        k=k+1
      ENDIF
    ENDIF
  ENDDO
  !
  !Allocate Q to store atom positions of P + inserted atoms
  ALLOCATE( Q( SIZE(P,1)+k,4 ) )
  Q(:,:) = 0.d0
  !If auxiliary properties are present they must be copied for inserted atoms
  IF(ALLOCATED(AUX)) THEN
    ALLOCATE( newAUX( SIZE(P(:,1))+k, SIZE(AUX(1,:)) ) )
    newAUX(:,:) = 0.d0
    DO i=1,SIZE(AUX(:,1))
      newAUX(i,:) = AUX(i,:)
    ENDDO
  ENDIF
  !
  !Also create temporary array T for shells if relevant
  IF( doshells ) THEN
    k = 0
    DO i=1,SIZE(S,1)
      IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
        IF( S(i,a1)>=pos1+DABS(b(a1))            .AND.                  &
          & S(i,a1)<pos1+2.d0*DABS(b(a1))-0.01d0 .AND.                  &
          & (S(i,a2)-pos2)*b(a1)/DABS(b(a1))>=0.d0            ) THEN
          k=k+1
        ENDIF
      ENDIF
    ENDDO
    !
    !Allocate T to store shells positions of S + inserted shells
    ALLOCATE( T( SIZE(S,1)+k,4 ) )
    T(:,:) = 0.d0
  ENDIF
  !
  IF( aniso ) THEN
    !Anisotropic elasticity
    !Store positions of inserted atoms into Q
    k = 0
    DO i=1,SIZE(P,1)
      IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
        IF( P(i,a1)>=pos1+DABS(b(a1))            .AND.                  &
          & P(i,a1)<pos1+2.d0*DABS(b(a1))-0.01d0 .AND.                  &
          & (P(i,a2)-pos2)*b(a1)/DABS(b(a1))>=0.d0            ) THEN
          k=k+1
          n = SIZE(P(:,1))+k
          Q(n,:) = P(i,:)
          !Apply edge displacements
          disp = ANISO_DISP(n,Q(n,:),a1,a2,a3,pos1,pos2,A_kn,Dn,Pn)
          Q(n,1:3) = Q(n,1:3) + disp(:)
          !Along a1, positions are those of P shifted by 3/4 Burgers vector
          Q(n,a1) = P(i,a1) -DABS(b(a1))*0.75d0
          !
          !Apply displacements to shell
          IF( doshells ) THEN
            T(n,:) = S(i,:)
            T(n,1:3) = T(n,1:3) + disp(:)
            T(n,a1) = S(i,a1) -DABS(b(a1))*0.75d0
          ENDIF
          !
          !Copy auxiliary properties for inserted atoms
          IF(ALLOCATED(AUX)) THEN
            newAUX(n,:) = AUX(i,:)
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    !
    CALL ATOMSK_MSG(2063,(/''/),(/DBLE(k)/))
    !
    !Then we apply atomic displacements to all atoms
    !(except the ones we just inserted) and store them into Q
    DO i=1,SIZE(P,1)
      Q(i,:) = P(i,:)
      IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
        IF( P(i,a1)>=pos1 .AND. (P(i,a2)-pos2)*b(a1)/DABS(b(a1))>=0.d0 ) THEN
          Q(i,a1) = P(i,a1)+DABS(b(a1))
          IF( doshells ) THEN
            T(i,a1) = S(i,a1)+DABS(b(a1))
          ENDIF
        ENDIF
        disp = ANISO_DISP(i,Q(i,:),a1,a2,a3,pos1,pos2,A_kn,Dn,Pn)
        Q(i,1:3) = Q(i,1:3) + disp(:)
        !Apply displacements to shell
        IF( doshells ) THEN
          T(i,1:3) = T(i,1:3) + disp(:)
        ENDIF
      ENDIF
    ENDDO
    !
    !
  ELSE
    !Isotropic elasticity
    !Store positions of inserted atoms into Q
    k = 0
    DO i=1,SIZE(P,1)
      IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
        IF( P(i,a1)>=pos1+DABS(b(a1))            .AND.                  &
          & P(i,a1)<pos1+2.d0*DABS(b(a1))-0.01d0 .AND.                  &
          & (P(i,a2)-pos2)*b(a1)/DABS(b(a1))>=0.d0            ) THEN
          k=k+1
          n = SIZE(P,1)+k
          Q(n,:) = P(i,:)
          !Apply edge displacements along a2
          disp = DISPEDGE(n,Q(n,:),a1,a2,b(a1),nu,pos1,pos2)
          Q(n,1:3) = Q(n,1:3) + disp(:)
          !Along a1, positions are those of P shifted by 5/4 Burgers vector
          Q(n,a1) = P(i,a1) -DABS(b(a1))*1.25d0
          !Apply displacements to shells if relevant
          IF( doshells ) THEN
            T(n,:) = S(i,:)
            T(n,1:3) = T(n,1:3) + disp
            !Along a1, positions are those of P shifted by 5/4 Burgers vector
            T(n,a1) = S(i,a1) -DABS(b(a1))*1.25d0
          ENDIF
          !Copy auxiliary properties for inserted atoms
          IF(ALLOCATED(AUX)) THEN
            newAUX(n,:) = AUX(i,:)
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    !
    CALL ATOMSK_MSG(2063,(/''/),(/DBLE(k)/))
    !
    !Then we apply atomic displacements to all atoms
    !(except the ones we just inserted) and store them into Q
    DO i=1,SIZE(P,1)
      Q(i,:) = P(i,:)
      IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
        IF( P(i,a1)>=pos1 ) THEN
          Q(i,a1) = P(i,a1)+DABS(b(a1))/2.d0
        ENDIF
      ENDIF
      disp = DISPEDGE(i,Q(i,:),a1,a2,b(a1),nu,pos1,pos2)
      Q(i,1:3) = Q(i,1:3) + disp(:)
      !Apply displacements to shell
      IF( doshells ) THEN
        T(i,1:3) = T(i,1:3) + disp(:)
      ENDIF
    ENDDO
    !
  ENDIF  !Endif aniso
  !
  !Replace old P with the new Q
  DEALLOCATE(P)
  ALLOCATE(P(SIZE(Q,1),4))
  P(:,:) = Q(:,:)
  DEALLOCATE(Q)
  !Replace old AUX by newAUX
  IF(ALLOCATED(AUX)) THEN
    DEALLOCATE(AUX)
    ALLOCATE( AUX( SIZE(newAUX,1), SIZE(newAUX,2) ) )
    AUX = newAUX
    DEALLOCATE(newAUX)
  ENDIF
  !
  !If shells exist, replace old S with new T
  IF( doshells ) THEN
    DEALLOCATE(S)
    ALLOCATE(S(SIZE(T,1),4))
    S(:,:) = T(:,:)
    DEALLOCATE(T)
  ENDIF
  !
  !Supercell has been expanded along a1: calculate the new length
  H(a1,a1) = H(a1,a1)+DABS(b(a1))/2.d0
  CALL ATOMSK_MSG(2064,(/''/),(/0.d0/))
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  EDGE DISLOCATION  !!!!!!!!!!
! Method 2: inserting a surface step (NP constant)
ELSEIF(disloctype=='edge2') THEN
  bsign = b(a1)/DABS(b(a1))
  !
  IF( aniso ) THEN
    !Anisotropic elasticity
    DO i=1,SIZE(P,1)
      IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
        disp = ANISO_DISP(i,P(i,:),a1,a2,a3,pos1,pos2,A_kn,Dn,Pn)
        P(i,1:3) = P(i,1:3) + disp(:)
        !
        !Same if shells exist
        IF( doshells ) THEN
          S(i,1:3) = S(i,1:3) + disp(:)
        ENDIF
      ENDIF
    ENDDO
    !
  ELSE
    !Isotropic elasticity
    DO i=1,SIZE(P,1)
      IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
        IF( P(i,a1)>pos1 .AND. P(i,a2)*bsign<pos2*bsign ) THEN
          P(i,a1) = P(i,a1)+DABS(b(a1))/2.d0
          IF( doshells ) THEN
            S(i,a1) = S(i,a1)+DABS(b(a1))/2.d0
          ENDIF
        ENDIF
        !
        disp = DISPEDGE(i,P(i,:),a1,a2,b(a1),nu,pos1,pos2)
        P(i,1:3) = P(i,1:3) + disp(:)
        IF( doshells ) THEN
          S(i,1:3) = S(i,1:3) + disp(:)
        ENDIF
        !
        IF( P(i,a1)>pos1 .AND. P(i,a2)*bsign>=pos2*bsign ) THEN
          P(i,a1) = P(i,a1)-DABS(b(a1))/2.d0
          IF( doshells ) THEN
            S(i,a1) = S(i,a1)-DABS(b(a1))/2.d0
          ENDIF
        ENDIF
        !
      ENDIF
    ENDDO
    !
  ENDIF
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!  MIXED DISLOCATION  !!!!!!!!!!
ELSE
  IF( aniso ) THEN
    !Anisotropic elasticity
    DO i=1,SIZE(P,1)
      IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
        disp = ANISO_DISP(i,P(i,:),a1,a2,a3,pos1,pos2,A_kn,Dn,Pn)
        P(i,1:3) = P(i,1:3) + disp(:)
        IF( doshells ) THEN
          S(i,1:3) = S(i,1:3) + disp(:)
        ENDIF
      ENDIF
    ENDDO
    !
  ELSE
    !Isotropic elasticity: cannot build a mixed dislocation, abort
    nerr = nerr+1
    CALL ATOMSK_MSG(2808,(/''/),(/0.d0/))
    GOTO 1000
  ENDIF
!
!
ENDIF
!
!Dislocation was successfully created
CALL ATOMSK_MSG(2065,(/''/),(/0.d0/))
!
!
!
400 CONTINUE
!Compute the theoretical (continuum) dislocation stresses
!and energy factor associated with this dislocation.
!Formula are given below, again see Hirth and Lothe.
!Only the Voigt components of the stresses are output
!as auxiliary properties, namely:
! sigma_xx sigma_yy sigma_zz sigma_xy sigma_yz sigma_xz
!If these auxiliary properties already exist, the present
!stresses must be added to the existing ones.
!
sigma(:,:) = 0.d0
!Setup array AUX for disloc stresses
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
    ALLOCATE(newAUX( SIZE(P,1), SIZE(AUXNAMES)+6 ) )
    newAUX(:,:) = 0.d0
    DO i=1,SIZE(AUXNAMES)
      newAUXNAMES(i) = AUXNAMES(i)
    ENDDO
    DO i=1,SIZE(P,1)
      newAUX(i,1:sig1-1) = AUX(i,:)
    ENDDO
    !Replace old arrays with new ones
    DEALLOCATE(AUXNAMES)
    DEALLOCATE(AUX)
    ALLOCATE(AUXNAMES( SIZE(newAUXNAMES) ))
    ALLOCATE(AUX( SIZE(newAUX,1), SIZE(newAUX,2) ))
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
IF(.NOT.aniso) THEN
  !isotropic case
  !Note: since mu is unknown we use mu=1 here, i.e. what is actually
  !     computed are the sigma(i,j)/mu
  IF(disloctype=="screw") THEN
    !sigma_11 = sigma_22 = sigma_33 = sigma_12 = 0
    !sigma_13 = (-mu*b/2pi) * y/(x²+y²)
    !sigma_23 = (mu*b/2pi) * x/(x²+y²)
    DO u=1,SIZE(P,1)
      tempreal = (P(u,a1)-pos1)**2 + (P(u,a2)-pos2)**2
      !avoid division by zero
      IF(DABS(tempreal)>1.d-6) THEN
        sigma(2,3) = -1.d0*bsign*VECLENGTH(b)/(2.d0*pi) * (P(u,a2)-pos2)/tempreal
        sigma(1,3) = bsign*VECLENGTH(b)/(2.d0*pi) * (P(u,a1)-pos1)/tempreal
        AUX(u,sig4) = AUX(u,sig4) + sigma(2,3)
        AUX(u,sig5) = AUX(u,sig5) + sigma(1,3)
      ENDIF
    ENDDO
  !
  ELSEIF(disloctype(1:4)=="edge") THEN
    !sigma_11 = (-mu*b/2pi(1-nu)) * y*(3x²+y²)/(x²+y²)²
    !sigma_22 = (mu*b/2pi(1-nu)) * y*(x²-y²)/(x²+y²)²
    !sigma_12 = (mu*b/2pi(1-nu)) * x*(x²-y²)/(x²+y²)²
    !sigma_33 = nu*(sigma_11 + sigma_22)
    !sigma_13 = sigma_23 = 0
    DO u=1,SIZE(P,1)
      tempreal = (P(u,a1)-pos1)**2 + (P(u,a2)-pos2)**2
      !avoid division by zero
      IF(DABS(tempreal)>1.d-6) THEN
        sigma(1,1) = -1.d0*bsign*VECLENGTH(b)/(2.d0*pi*(1.d0-nu)) *          &
                 & (P(u,a2)-pos2)*(3.d0*(P(u,a1)-pos1)**2+(P(u,a2)-pos2)**2) &
                 & / (tempreal**2)
        sigma(2,2) = bsign*VECLENGTH(b)/(2.d0*pi*(1.d0-nu)) *                &
                  & (P(u,a2)-pos2)*((P(u,a1)-pos1)**2-(P(u,a2)-pos2)**2)     &
                  & / (tempreal**2)
        sigma(1,2) = bsign*VECLENGTH(b)/(2.d0*pi*(1.d0-nu)) *                &
                    & (P(u,a1)-pos1)*((P(u,a1)-pos1)**2-(P(u,a2)-pos2)**2)   &
                    & / (tempreal**2)
        AUX(u,sig1) = AUX(u,sig1) + sigma(1,1)
        AUX(u,sig2) = AUX(u,sig2) + sigma(2,2)
        AUX(u,sig3) = AUX(u,sig3) + nu*( sigma(1,1) + sigma(2,2) )
        AUX(u,sig6) = AUX(u,sig6) + sigma(1,2)
      ENDIF
    ENDDO
  ENDIF
  !
ELSE
  !anisotropic case
  !Dislocation stresses are given by:
  !sigma_ij = Re{ (-1/2ipi) SUM(n=1,3) B_ijk(n)A_k(n)D(n)/(x1+P(n)x2) }
  DO u=1,SIZE(P,1)
    !Compute stresses
    DO i=1,3
      DO j=1,3
        r = ELASTINDEX(i,j)    ! r = ij
        tempcmplx = DCMPLX(0.d0,0.d0)
        DO n=1,3
          DO k=1,3
            tempcmplx = tempcmplx + B_ijk(r,k,n)*A_kn(k,n)*Dn(n)/ &
                      & ( (P(u,a1)-pos1)+Pn(n)*(P(u,a2)-pos2) )
          ENDDO
        ENDDO
        !Keep only the real part
        sigma(i,j) = DBLE(tempcmplx)
      ENDDO
    ENDDO
    !Save stresses to AUX
    AUX(u,sig1) = AUX(u,sig1) + sigma(1,1)
    AUX(u,sig2) = AUX(u,sig2) + sigma(2,2)
    AUX(u,sig3) = AUX(u,sig3) + sigma(3,3)
    AUX(u,sig4) = AUX(u,sig4) + sigma(2,3)
    AUX(u,sig5) = AUX(u,sig5) + sigma(1,3)
    AUX(u,sig6) = AUX(u,sig6) + sigma(1,2)
  ENDDO !end loop on all atoms
  !
  !Now compute the factor Kfactor = Kb²
  !Energy factor is given by:
  !  Kb² = b_i * Im{ SUM(n=1,3) B_i2k A_k(n) D(n) }
  DO i=1,3
    r = ELASTINDEX(i,2)     ! j=2
    tempcmplx = DCMPLX(0.d0,0.d0)
    DO n=1,3
      DO k=1,3
        tempcmplx = tempcmplx + B_ijk(r,k,n)*A_kn(k,n)*Dn(n)
      ENDDO
    ENDDO
    Kfactor(i) = b(i) * AIMAG(tempcmplx)
  ENDDO
  tempreal = VECLENGTH(Kfactor)
  !CALL ATOMSK_MSG(2101,(/"Kb²"/),(/tempreal/))
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
END SUBROUTINE DISLOC_XYZ
!
!
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
! Formulae can be found for instance in
! J.P. Hirth, J. Lothe, 'Theory of dislocations',
! 1st ed. (1968), p.75
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
! ANISO_DISP
! This function calculates the displacements due to a
! dislocation in an anisotropic medium:
!   u(k) = Re{ [-1/(2*i*pi)]*[SUM(n=1,3) A_k(n)*D(n)
!                                  *Ln(x1+P(n)*x2)] }
! The complex coefficients P(n), A_k(n), D(n) must
! be provided as input (see routine ANISO_COEFF below).
! Formulae can be found for instance in
! J.P. Hirth, J. Lothe, 'Theory of dislocations',
! 1st ed. (1968), p.426.
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
! The numerotation and equation numbers below follow
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
REAL(dp),PARAMETER:: Ak_threshold=1.d-14 !threshold for normalization of the A_k(n)
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
  IF( DBLE(A_3)>Ak_threshold ) THEN
    !We are not dividing by zero => normalize the A_k(n) to A3(n)
    msg = "Normalizing by A3(n)..."
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !In this case A3(n)=1
    A_kn(1,n) = A_1/A_3
    A_kn(2,n) = -1.d0*A_2/A_3
    A_kn(3,n) = DCMPLX(1.d0,0.d0)
    !
  ELSEIF( DBLE(A_1)>Ak_threshold ) THEN
    !A3(n) is zero => try to divide by A1
    msg = "Normalizing by A1(n)..."
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !In this case A1(n)=1
    A_kn(1,n) = DCMPLX(1.d0,0.d0)
    A_kn(2,n) = -1.d0*A_2/A_1
    A_kn(3,n) = A_3/A_1
    !
  ELSEIF( DBLE(A_2)>Ak_threshold ) THEN
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
END MODULE dislocation