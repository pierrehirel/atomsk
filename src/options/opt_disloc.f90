MODULE dislocation
!
!**********************************************************************************
!*  DISLOCATION                                                                   *
!**********************************************************************************
!* This module reads atomic positions from an array P and introduces a            *
!* straight dislocation, or a dislocation loop, of Burgers vector "b".            *
!* FOR A STRAIGHT DISLOCATION:                                                    *
!* dislocation lies along 'dislocline', plane of cut is 'dislocplane'.            *
!* Position of the dislocation is X=pos(1) and Y=pos(2) if dislocline='Z'         *
!* (otherwise, circular permutation applies).                                     *
!* If "C_tensor" only contains zeros, isotropic elasticity is used                *
!* (cf subroutines DISPEDGE or DISPSCREW in module "dislocation_iso").            *
!* If "C_tensor" is non-zero, then anisotropic elasticity is used                 *
!* (cf routines ANISO_COEFF and ANISO_DISP in module "dislocation_aniso").        *
!* FOR A DISLOCATION LOOP:                                                        *
!* 'dislocline' is the direction normal to the loop, and to the plane of cut.     *
!* Center of the loop is placed at X=pos(1), Y=pos(2), Z=pos(3).                  *
!* Radius of the loop is pos(4).                                                  *
!* (cf routines in module "dislocation_loop").                                    *
!**********************************************************************************
!* (C) May 2010 - Pierre Hirel                                                    *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 26 Feb. 2025                                     *
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
USE messages
USE files_msg
USE resize
USE subroutines
USE dislocation_iso
USE dislocation_aniso
USE dislocation_loop
!
!
CONTAINS
!
!
SUBROUTINE DISLOC_XYZ(H,P,S,disloctype,dislocline,dislocplane,b,nu,pos,SELECT,ORIENT,AUXNAMES,AUX,C_tensor)
!
!
IMPLICIT NONE
!Input variables
CHARACTER(LEN=*),INTENT(IN):: dislocline !direction of dislocation line: x, y, z, or Miller index, or file name
CHARACTER(LEN=*),INTENT(IN):: dislocplane !normal to plane of cut: x, y, z, or Miller index
CHARACTER(LEN=*),INTENT(IN):: disloctype !type of dislocation: screw, edge, edge_add, edge_rm, mixed, loop, file
REAL(dp),INTENT(IN):: nu   !Poisson ratio of the material
                           !(not used if disloctype=screw, or if aniso=.TRUE.)
REAL(dp),DIMENSION(3,3),INTENT(IN):: ORIENT   !crystal orientation
REAL(dp),DIMENSION(9,9),INTENT(IN):: C_tensor !Elastic tensor (GPa)
REAL(dp),DIMENSION(3),INTENT(IN):: b    !Burgers vector (angströms)
!
!Internal variables
CHARACTER(LEN=1):: dir1, dir2, dir3 !directions x, y, z
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES, newAUXNAMES !names of auxiliary properties
LOGICAL:: aniso !shall anisotropic elasticity be used?
LOGICAL:: doaux    !treat auxiliary properties?
LOGICAL:: doshells !apply displacements also to shells?
LOGICAL:: rotate  !rotate displacements?
LOGICAL:: straight !construct a straight dislocation line?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT     !mask for atom list
LOGICAL,DIMENSION(:),ALLOCATABLE:: newSELECT  !mask for atom list (temporary)
INTEGER:: a1, a2, a3 !a3=1,2,3 if dislocline=x,y,z respectively
INTEGER:: i, j, k, n, r, u
INTEGER:: Ndisloc  !number of dislocations successfully inserted
INTEGER:: sig1, sig2, sig3, sig4, sig5, sig6 !for indexing stresses
REAL(dp):: bsign !sign of Burgers vector (+1.d0 or -1.d0)
REAL(dp),DIMENSION(5):: pos !if straight: pos(1:1) = position in plane normal to dislocline
                            !if segment: pos(1:2) = position in plane normal to dislocline;
                            !            pos(3:4) = min/max along dislocline
                            !if loop: pos(1:3) = coord. (x,y,z) of center of loop;
                            !         pos(4)   = radius of loop
REAL(dp):: Efactor  !prelogarithmic energy factor
REAL(dp):: mu       !shear modulus (GPa)
REAL(dp):: tempreal
REAL(dp):: V1, V2, V3
REAL(dp),DIMENSION(3):: disp     !elastic displacement applied to an atom
REAL(dp),DIMENSION(1,3):: Vline, Vplane  !vector along disloc. line, normal to cut plane
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENTN    !crystal orientation
REAL(dp),DIMENSION(3,3):: rot_matrix !rotation matrix
REAL(dp),DIMENSION(3,3):: disp_rot_matrix !rotation matrix for displacements
REAL(dp),DIMENSION(3,3):: sigma      !dislocation theoretical elastic stresses
REAL(dp),DIMENSION(9,9):: C_tensor_rot !rotated tensor if dislocline is not along Z
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S   !Positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Q, T                 !Positions of atoms, shells (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX    !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX               !auxiliary properties (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: xLoop !coordinates of points forming the loop
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: disarrayb      !Burgers vectors of dislocation segments
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE:: disarraypos  !coordinates of segments forming dislocations
COMPLEX(dp):: tempcmplx
COMPLEX(dp),DIMENSION(3):: Dn, Pn    !anisotropy coefficients D(n), P(n)-
COMPLEX(dp),DIMENSION(3,3):: A_kn    !-, A_k(n) -
COMPLEX(dp),DIMENSION(9,3,3):: B_ijk !-, and B_ijk(n)
!
!Initialize variables
a1 = 0
a2 = 0
a3 = 0
k = 0
Ndisloc = 0
sig1=0
sig2=0
sig3=0
sig4=0
sig5=0
sig6=0
bsign = 1.d0
mu = 1.d0    !(dummy value)
disp_rot_matrix(:,:) = 0.d0
Vline(:,:) = 0.d0
Vplane(:,:) = 0.d0
aniso = .FALSE.
rotate = .FALSE.
straight = .TRUE.
doaux = .FALSE.
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(T)) DEALLOCATE(T)
IF(ALLOCATED(newSELECT)) DEALLOCATE(newSELECT)
IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
IF(ALLOCATED(newAUXNAMES)) DEALLOCATE(newAUXNAMES)
IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
  doshells=.TRUE.
ELSE
  doshells=.FALSE.
ENDIF
!
!
IF( verbosity>=4 ) THEN
  msg = 'Entering DISLOC_XYZ'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "disloctype:  ", TRIM(disloctype)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,*) "dislocline:  ", TRIM(dislocline)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,*) "dislocplane: ", TRIM(dislocplane)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,'(a6,3f9.3,a1)') "b = ( ", b(1), b(2), b(3), ")"
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,*) "nu = ", nu
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,'(a6,3f9.3,a1)') "position = ( ", pos(1), pos(2), pos(3), ")"
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,*) "Loop radius = ", pos(4)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ENDIF
!
!
!
100 CONTINUE
IF( disloctype=="loop" .OR. disloctype=="array" .OR. disloctype=="file" ) THEN
  k = 0
  straight = .FALSE.
  !
ELSE
  !A straight dislocation must be constructed
  !Define the axes: a3=dislocation line
  Vline(:,:) = 0.d0
  IF(dislocline=='x' .OR. dislocline=='X') THEN
    a3 = 1
    dir3 = "x"
    Vline(1,1) = 1.d0
  ELSEIF(dislocline=='y' .OR. dislocline=='Y') THEN
    a3 = 2
    dir3 = "y"
    Vline(1,2) = 1.d0
  ELSEIF(dislocline=='z' .OR. dislocline=='Z') THEN
    a3 = 3
    dir3 = "z"
    Vline(1,3) = 1.d0
  ELSE
    !It ought to be a vector given by Miller indices
    !Convert "dislocline" into a Cartesian vector and save it in Vline(1,:)
    CALL MILLER2VEC(H,dislocline,ORIENT,Vline(1,:),j)
    !
    !Check return status j (0=success, otherwise there was an error)
    IF( j>0 ) THEN
      IF( j==2 ) THEN
        !The error was because i is not equal to -h-k
        nerr=nerr+1
        CALL ATOMSK_MSG(815,(/dislocline/),(/0.d0/))
        GOTO 1000
      ELSE
        !Other error, unable to convert this string into a proper vector
        CALL ATOMSK_MSG(817,(/TRIM(dislocline)/),(/0.d0/))
        GOTO 1000
      ENDIF
    ENDIF
    !
    IF( VECLENGTH( CROSS_PRODUCT(Vline(1,:),(/0.d0,0.d0,1.d0/)) ) > 0.1d0 ) THEN
      !Vline is not along Z => we will have to rotate displacements accordingly
      rotate = .TRUE.
    ENDIF
    !
  ENDIF
  !
  !a2=normal to plane of cut
  Vplane(:,:) = 0.d0
  !a1=last component
  IF(dislocplane=='x' .OR. dislocplane=='X') THEN
    a2 = 1
    dir2 = "x"
    Vplane(1,1) = 1.d0
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
    Vplane(1,2) = 1.d0
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
    Vplane(1,3) = 1.d0
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
    !It ought to be a vector given by Miller indices
    !Convert "dislocplane" into a Cartesian vector and save it in Vplane(1,:)
    CALL MILLER2VEC(H,dislocplane,ORIENT,Vplane(1,:),j)
    !
    !Check return status j (0=success, otherwise there was an error)
    IF( j>0 ) THEN
      IF( j==2 ) THEN
        !The error was because i is not equal to -h-k
        nerr=nerr+1
        CALL ATOMSK_MSG(815,(/dislocplane/),(/0.d0/))
        GOTO 1000
      ELSE
        !Other error, unable to convert this string into a proper vector
        CALL ATOMSK_MSG(817,(/TRIM(dislocplane)/),(/0.d0/))
        GOTO 1000
      ENDIF
    ENDIF
    !
    IF( VECLENGTH( CROSS_PRODUCT(Vplane(1,:),(/0.d0,1.d0,0.d0/)) ) > 0.1d0 ) THEN
      !Vplane is not along Y => we will have to rotate displacements accordingly
      rotate = .TRUE.
    ENDIF
    !
  ENDIF
  !
  !Debug messages
  IF( verbosity>=4 ) THEN
    WRITE(msg,*) "a1, a2, a3: ", a1, a2, a3
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,*) "Vline = (", Vline(1,:), " )"
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,*) "Vplane = (", Vplane(1,:), " )"
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  ENDIF
  !
  IF( a1==0 .OR. a2==0 .OR. a3==0 ) THEN
  
  ENDIF
  !
  !Check if elastic tensor is defined
  IF( ANY( C_tensor(1:3,1:3).NE.0.d0 ) ) THEN
    aniso = .TRUE.
    k = 1 !only used by message below
  ELSE
    aniso = .FALSE.
    k = 0 !only used by message below
  ENDIF
ENDIF
!
!
!Print messages
CALL ATOMSK_MSG(2061,(/disloctype,dislocline/),(/b(1),b(2),b(3),DBLE(k),pos(1),pos(2),pos(3),pos(4)/))
!
IF( .NOT.ALLOCATED(P) .OR. SIZE(P,1)<=0 ) THEN
  !No atom in system: can not apply option
  GOTO 1000
ENDIF
!
!
!
!
IF( straight ) THEN
  !If Burgers vector is zero then skip the whole thing
  IF( VECLENGTH(b)==0.d0) THEN
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2725,(/''/),(/0.d0/))
    GOTO 1000
  ENDIF
  !
  !Check that the directions normal to dislocation line are not too small
  IF( VECLENGTH(H(:,a1))<3.d0*VECLENGTH(b(:)) .OR. VECLENGTH(H(:,a2))<3.d0*VECLENGTH(b(:)) ) THEN
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2726,(/''/),(/0.d0/))
  ENDIF
  !Check that the position for the dislocation is inside the box
  IF( pos(1)<MINVAL(P(:,a1)) .OR. pos(1)>MAXVAL(P(:,a1)) .OR. &
    & pos(2)<MINVAL(P(:,a2)) .OR. pos(2)>MAXVAL(P(:,a2))     ) THEN
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2754,(/'dislocation'/),(/0.d0/))
  ENDIF
ELSEIF( disloctype=="loop" ) THEN
  IF( DABS(pos(4)) <= 2.d0 ) THEN
    !Dislocation loop radius is extremely small and loop will vanish
    !=> display a warning and skip this option
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2759,(/''/),(/0.d0/))
    GOTO 1000
  ENDIF
ENDIF
!
!
!
200 CONTINUE
IF( straight ) THEN
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
      !rotate the A_kn and B_ijk so that the displacements u_k
      !are applied in the correct direction
      IF(a3.NE.3) THEN
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
  !Prepare arrays AUX and AUXNAMES to store dislocation stresses.
  !Only the Voigt components of the stresses are output
  !as auxiliary properties, namely:
  ! sigma_xx sigma_yy sigma_zz sigma_xy sigma_yz sigma_xz
  !If these auxiliary properties already exist, the present
  !stresses must be added to the existing ones.
  sigma(:,:) = 0.d0
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
ENDIF
!
!
!
300 CONTINUE
IF(disloctype=='screw') THEN
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!  SCREW DISLOCATION  !!!!!!!!!!
  ! Apply atomic displacements of a straight screw dislocation line
  bsign = b(a3)/DABS(b(a3))
  !
  IF( aniso ) THEN
    !Anisotropic elasticity
    DO i=1,SIZE(P,1)
      IF(IS_SELECTED(SELECT,i)) THEN
        disp = ANISO_DISP(i,P(i,:),a1,a2,a3,pos(1),pos(2),A_kn,Dn,Pn)
        P(i,1:3) = P(i,1:3) + disp(:)
        !
        !Same if shells exist
        IF( doshells ) THEN
          S(i,1:3) = S(i,1:3) + disp(:)
        ENDIF
        !
        !Compute stress
        CALL ANISO_STRESS(P(i,1:3),a1,a2,a3,pos(1),pos(2),A_kn,B_ijk,Dn,Pn,sigma)
        !
        !Add stresses due to this disloc. to existing stresses in AUX
        AUX(i,sig1) = AUX(i,sig1) + sigma(1,1)
        AUX(i,sig2) = AUX(i,sig2) + sigma(2,2)
        AUX(i,sig3) = AUX(i,sig3) + sigma(3,3)
        AUX(i,sig4) = AUX(i,sig4) + sigma(2,3)
        AUX(i,sig5) = AUX(i,sig5) + sigma(1,3)
        AUX(i,sig6) = AUX(i,sig6) + sigma(1,2)
      ENDIF
    ENDDO
    !
    !Compute prelogarithmic energy factor
    Efactor = ANISO_EFACTOR(b,A_kn,B_ijk,Dn,Pn)
    CALL ATOMSK_MSG(2101,(/'Kb²/4pi'/),(/Efactor/))
    !
  ELSE
    !Isotropic elasticity
    DO i=1,SIZE(P,1)
      IF(IS_SELECTED(SELECT,i)) THEN
        disp = DISPSCREW(i,P(i,:),a1,a2,a3,b(a3),pos(1),pos(2))
        P(i,1:3) = P(i,1:3) + disp(:)
        !
        !Same if shells exist
        IF( doshells ) THEN
          S(i,1:3) = S(i,1:3) + disp(:)
        ENDIF
        !
        !Compute stress
        CALL STRESSSCREW(P(i,1:3),a1,a2,bsign*VECLENGTH(b),mu,nu,pos(1),pos(2),sigma)
        !
        !Add stresses due to this disloc. to existing stresses in AUX
        AUX(i,sig4) = AUX(i,sig4) + sigma(2,3)
        AUX(i,sig5) = AUX(i,sig5) + sigma(1,3)
      ENDIF
    ENDDO
    !
    !Compute prelogarithmic energy factor
    Efactor = ISO_EFACTOR(VECLENGTH(b),mu,nu,0.d0)
    CALL ATOMSK_MSG(2101,(/'b²/4pi'/),(/Efactor/))
  ENDIF
  !
  Ndisloc = 1
  !
  !
ELSEIF(disloctype=='edge') THEN
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!  EDGE DISLOCATION  !!!!!!!!!!
  ! Method 1: only apply displacements (NP constant)
  ! NOTE: this will introduce a step on one boundary
  bsign = b(a1)/DABS(b(a1))
  !
  IF( aniso ) THEN
    !Anisotropic elasticity
    DO i=1,SIZE(P,1)
      IF(IS_SELECTED(SELECT,i)) THEN
        disp = ANISO_DISP(i,P(i,:),a1,a2,a3,pos(1),pos(2),A_kn,Dn,Pn)
        P(i,1:3) = P(i,1:3) + disp(:)
        !
        !Same if shells exist
        IF( doshells ) THEN
          S(i,1:3) = S(i,1:3) + disp(:)
        ENDIF
        !
        !Compute stress
        CALL ANISO_STRESS(P(i,1:3),a1,a2,a3,pos(1),pos(2),A_kn,B_ijk,Dn,Pn,sigma)
        !
        !Add stresses due to this disloc. to existing stresses in AUX
        AUX(i,sig1) = AUX(i,sig1) + sigma(1,1)
        AUX(i,sig2) = AUX(i,sig2) + sigma(2,2)
        AUX(i,sig3) = AUX(i,sig3) + sigma(3,3)
        AUX(i,sig4) = AUX(i,sig4) + sigma(2,3)
        AUX(i,sig5) = AUX(i,sig5) + sigma(1,3)
        AUX(i,sig6) = AUX(i,sig6) + sigma(1,2)
      ENDIF
    ENDDO
    !
    !Compute prelogarithmic energy factor
    Efactor = ANISO_EFACTOR(b,A_kn,B_ijk,Dn,Pn)
    CALL ATOMSK_MSG(2101,(/'Kb²/4pi'/),(/Efactor/))
    !
    !
  ELSE
    !Isotropic elasticity
    DO i=1,SIZE(P,1)
      IF(IS_SELECTED(SELECT,i)) THEN
        IF( P(i,a1)>pos(1) .AND. (P(i,a2)-pos(2))*bsign<=0.d0 ) THEN
          P(i,a1) = P(i,a1)+DABS(b(a1))/2.d0
          IF( doshells ) THEN
            S(i,a1) = S(i,a1)+DABS(b(a1))/2.d0
          ENDIF
        ENDIF
        !
        disp = DISPEDGE(i,P(i,:),a1,a2,b(a1),nu,pos(1),pos(2))
        P(i,1:3) = P(i,1:3) + disp(:)
        IF( doshells ) THEN
          S(i,1:3) = S(i,1:3) + disp(:)
        ENDIF
        !
        IF( P(i,a1)>pos(1) .AND. (P(i,a2)-pos(2))*bsign>=0.d0 ) THEN
          P(i,a1) = P(i,a1)-DABS(b(a1))/2.d0
          IF( doshells ) THEN
            S(i,a1) = S(i,a1)-DABS(b(a1))/2.d0
          ENDIF
        ENDIF
        !
        !
        !Compute stress
        CALL STRESSEDGE(P(i,1:3),a1,a2,bsign*VECLENGTH(b),mu,nu,pos(1),pos(2),sigma)
        !
        !Add stresses due to this disloc. to existing stresses in AUX
        AUX(i,sig1) = AUX(i,sig1) + sigma(1,1)
        AUX(i,sig2) = AUX(i,sig2) + sigma(2,2)
        AUX(i,sig3) = AUX(i,sig3) + sigma(3,3)
        AUX(i,sig6) = AUX(i,sig6) + sigma(1,2)
      ENDIF
    ENDDO
    !
    !Compute prelogarithmic energy factor
    Efactor = ISO_EFACTOR(VECLENGTH(b),mu,nu,pi/4.d0)
    CALL ATOMSK_MSG(2101,(/'b²/4pi(1-nu)'/),(/Efactor/))
    !
  ENDIF
  !
  Ndisloc = 1
  !
  !
ELSEIF(disloctype=='edge_add' .OR. disloctype=='edge-add') THEN
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!  EDGE DISLOCATION  !!!!!!!!!!
  ! Method 2: insert a plane of atoms
  ! (NP is increased so we have to use a new array Q)
  bsign = b(a1)/DABS(b(a1))
  !First, insert a plane of atoms by copying an existing plane
  !Count k = number of atoms that will be inserted
  k = 0
  DO i=1,SIZE(P,1)
    IF(IS_SELECTED(SELECT,i)) THEN
      IF( P(i,a1)>=pos(1)+DABS(b(a1))            .AND.                  &
        & P(i,a1)<pos(1)+2.d0*DABS(b(a1))-0.01d0 .AND.                  &
        & (P(i,a2)-pos(2))*bsign>=0.d0                ) THEN
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
  !If array SELECT is defined, expand it
  IF( ALLOCATED(SELECT) ) THEN
    CALL RESIZE_LOGICAL1(SELECT,SIZE(P,1)+k,i)
    IF( i>0 ) THEN
      nerr=nerr+1
      CALL ATOMSK_MSG(818,(/"SELECT"/),(/0.d0/))
      GOTO 1000
    ENDIF
  ENDIF
  !
  !Also create temporary array T for shells if relevant
  IF( doshells ) THEN
    k = 0
    DO i=1,SIZE(S,1)
      IF(IS_SELECTED(SELECT,i)) THEN
        IF( S(i,a1)>=pos(1)+DABS(b(a1))            .AND.                  &
          & S(i,a1)<pos(1)+2.d0*DABS(b(a1))-0.01d0 .AND.                  &
          & (S(i,a2)-pos(2))*bsign>=0.d0                ) THEN
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
      IF(IS_SELECTED(SELECT,i)) THEN
        IF( P(i,a1)>=pos(1)+DABS(b(a1))            .AND.                  &
          & P(i,a1)<pos(1)+2.d0*DABS(b(a1))-0.01d0 .AND.                  &
          & (P(i,a2)-pos(2))*bsign>=0.d0                ) THEN
          k=k+1
          n = SIZE(P(:,1))+k
          Q(n,:) = P(i,:)
          !Apply edge displacements
          disp = ANISO_DISP(n,Q(n,:),a1,a2,a3,pos(1),pos(2),A_kn,Dn,Pn)
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
    DO i=1,SIZE(Q,1)
      IF( i <= SIZE(P,1) ) THEN
        Q(i,:) = P(i,:)
        IF(IS_SELECTED(SELECT,i)) THEN
          IF( P(i,a1)>=pos(1) .AND. (P(i,a2)-pos(2))*bsign>=0.d0 ) THEN
            Q(i,a1) = P(i,a1)+DABS(b(a1))
            IF( doshells ) THEN
              T(i,a1) = S(i,a1)+DABS(b(a1))
            ENDIF
          ENDIF
          disp = ANISO_DISP(i,Q(i,:),a1,a2,a3,pos(1),pos(2),A_kn,Dn,Pn)
          Q(i,1:3) = Q(i,1:3) + disp(:)
          !Apply displacements to shell
          IF( doshells ) THEN
            T(i,1:3) = T(i,1:3) + disp(:)
          ENDIF
        ENDIF
        !
        !Compute stress
        CALL ANISO_STRESS(Q(i,1:3),a1,a2,a3,pos(1),pos(2),A_kn,B_ijk,Dn,Pn,sigma)
        !
        !Add stresses due to this disloc. to existing stresses in AUX
        newAUX(i,sig1) = newAUX(i,sig1) + sigma(1,1)
        newAUX(i,sig2) = newAUX(i,sig2) + sigma(2,2)
        newAUX(i,sig3) = newAUX(i,sig3) + sigma(3,3)
        newAUX(i,sig4) = newAUX(i,sig4) + sigma(2,3)
        newAUX(i,sig5) = newAUX(i,sig5) + sigma(1,3)
        newAUX(i,sig6) = newAUX(i,sig6) + sigma(1,2)
      ENDIF
    ENDDO
    !
    !Compute prelogarithmic energy factor
    Efactor = ANISO_EFACTOR(b,A_kn,B_ijk,Dn,Pn)
    CALL ATOMSK_MSG(2101,(/'Kb²/4pi'/),(/Efactor/))
    !
    !
    !
  ELSE
    !Isotropic elasticity
    !Store positions of inserted atoms into Q
    k = 0
    DO i=1,SIZE(P,1)
      IF(IS_SELECTED(SELECT,i)) THEN
        IF( P(i,a1)>=pos(1)+DABS(b(a1))            .AND.                  &
          & P(i,a1)<pos(1)+2.d0*DABS(b(a1))-0.01d0 .AND.                  &
          & (P(i,a2)-pos(2))*bsign>=0.d0                ) THEN
          k=k+1
          n = SIZE(P,1)+k
          Q(n,:) = P(i,:)
          !Apply edge displacements along a2
          disp = DISPEDGE(n,Q(n,:),a1,a2,b(a1),nu,pos(1),pos(2))
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
    DO i=1,SIZE(Q,1)
      IF( i <= SIZE(P,1) ) THEN
        Q(i,:) = P(i,:)
        IF(IS_SELECTED(SELECT,i)) THEN
          IF( P(i,a1)>=pos(1) ) THEN
            Q(i,a1) = P(i,a1)+DABS(b(a1))/2.d0
            IF( doshells ) THEN
              T(i,a1) = S(i,a1)-DABS(b(a1))/2.d0
            ENDIF
          ENDIF
        ENDIF
        disp = DISPEDGE(i,Q(i,:),a1,a2,b(a1),nu,pos(1),pos(2))
        Q(i,1:3) = Q(i,1:3) + disp(:)
        !Apply displacements to shell
        IF( doshells ) THEN
          T(i,1:3) = T(i,1:3) + disp(:)
        ENDIF
      ENDIF
      !
      !Compute stress
      CALL STRESSEDGE(Q(i,1:3),a1,a2,bsign*VECLENGTH(b),mu,nu,pos(1),pos(2),sigma)
      !
      !Add stresses due to this disloc. to existing stresses in AUX
      newAUX(i,sig1) = newAUX(i,sig1) + sigma(1,1)
      newAUX(i,sig2) = newAUX(i,sig2) + sigma(2,2)
      newAUX(i,sig3) = newAUX(i,sig3) + sigma(3,3)
      newAUX(i,sig6) = newAUX(i,sig6) + sigma(1,2)
    ENDDO
    !
    !Compute prelogarithmic energy factor
    Efactor = ISO_EFACTOR(VECLENGTH(b),mu,nu,pi/4.d0)
    CALL ATOMSK_MSG(2101,(/'b²/4pi(1-nu)'/),(/Efactor/))
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
  H(a1,a1) = H(a1,a1) + DABS(b(a1))/2.d0
  CALL ATOMSK_MSG(2064,(/"+b/2"/),(/DBLE(a1)/))
  !
  Ndisloc = 1
  !
  !
ELSEIF(disloctype=='edge_rm' .OR. disloctype=='edge-rm') THEN
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!  EDGE DISLOCATION  !!!!!!!!!!
  ! Method 3: remove a half-plane of atoms
  ! (NP is decreased so we have to use a new array Q)
  bsign = b(a1)/DABS(b(a1))
  !
  !Allocate Q to store atom positions
  !(for now it has the same size as P)
  ALLOCATE( Q( SIZE(P,1),4 ) )
  Q(:,:) = 0.d0
  !
  !If auxiliary properties are present they must be copied
  IF(ALLOCATED(AUX)) THEN
    ALLOCATE( newAUX( SIZE(P,1) , SIZE(AUX,2) ) )
    newAUX(:,:) = 0.d0
  ENDIF
  !Also create temporary array T for shells if relevant
  IF( doshells ) THEN
    !Allocate T to store shells positions of S + inserted shells
    ALLOCATE( T( SIZE(S,1),4 ) )
    T(:,:) = 0.d0
  ENDIF
  !If array SELECT is defined, expand it
  IF( ALLOCATED(SELECT) ) THEN
    ALLOCATE( newSELECT( SIZE(P,1)) )
    newSELECT(:) = .TRUE.
  ENDIF
  !
  !Remove atoms
  !Count j = index of current atom
  j = 0
  !Count k = number of removed atoms
  k = 0
  DO i=1,SIZE(P,1)
    IF(IS_SELECTED(SELECT,i)) THEN
      IF( P(i,a1)>=pos(1) .AND. P(i,a1)<pos(1)+DABS(b(a1))-0.01d0 .AND. &
        & (P(i,a2)-pos(2))*bsign<=0.d0                              ) THEN
        !This atom is removed
        k=k+1
      ELSE
        !This atom stays
        j=j+1
        Q(j,:) = P(i,:)
        !If auxiliary properties exist, copy them
        IF(ALLOCATED(AUX)) THEN
          newAUX(j,:) = AUX(i,:)
        ENDIF
        !If shells are present, copy them
        IF(doshells) THEN
          T(j,:) = S(i,:)
        ENDIF
        !If a selection was defined, copy it
        IF(ALLOCATED(SELECT)) THEN
          newSELECT(j) = SELECT(i)
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  !
  !Get rid of old atom positions P
  DEALLOCATE(P)
  !Copy the j first entries from Q into P
  ALLOCATE( P(j,4) )
  P(:,:) = 0.d0
  IF(ALLOCATED(AUX)) THEN
    DEALLOCATE(AUX)
    ALLOCATE(AUX(j,SIZE(newAUX,2)))
  ENDIF
  DO i=1,SIZE(P,1)
    P(i,:) = Q(i,:)
    IF(ALLOCATED(AUX)) THEN
      AUX(i,:) = newAUX(i,:)
    ENDIF
    IF(doshells) THEN
      S(i,:) = T(i,:)
    ENDIF
    IF(ALLOCATED(SELECT)) THEN
      SELECT(i) = newSELECT(i)
    ENDIF
  ENDDO
  !
  CALL ATOMSK_MSG(2063,(/''/),(/-1.d0*DBLE(k)/))
  !
  !We don't need Q anymore
  DEALLOCATE(Q)
  IF(ALLOCATED(T)) DEALLOCATE(T)
  IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
  IF(ALLOCATED(newSELECT)) DEALLOCATE(newSELECT)
  !
  IF( aniso ) THEN
    !Anisotropic elasticity
    !Apply elastic displacements to all atoms
    DO i=1,SIZE(P,1)
      IF(IS_SELECTED(SELECT,i)) THEN
        IF( P(i,a1)>=pos(1) .AND. (P(i,a2)-pos(2))*bsign<=0.d0 ) THEN
          P(i,a1) = P(i,a1)-DABS(b(a1))
          IF( doshells ) THEN
            S(i,a1) = S(i,a1)-DABS(b(a1))
          ENDIF
        ENDIF
        disp = ANISO_DISP(i,P(i,:),a1,a2,a3,pos(1),pos(2),A_kn,Dn,Pn)
        P(i,1:3) = P(i,1:3) + disp(:)
        !Apply displacements to shell
        IF( doshells ) THEN
          S(i,1:3) = S(i,1:3) + disp(:)
        ENDIF
      ENDIF
      !
      !Compute stress
      CALL ANISO_STRESS(P(i,1:3),a1,a2,a3,pos(1),pos(2),A_kn,B_ijk,Dn,Pn,sigma)
      !
      !Add stresses due to this disloc. to existing stresses in AUX
      AUX(i,sig1) = AUX(i,sig1) + sigma(1,1)
      AUX(i,sig2) = AUX(i,sig2) + sigma(2,2)
      AUX(i,sig3) = AUX(i,sig3) + sigma(3,3)
      AUX(i,sig4) = AUX(i,sig4) + sigma(2,3)
      AUX(i,sig5) = AUX(i,sig5) + sigma(1,3)
      AUX(i,sig6) = AUX(i,sig6) + sigma(1,2)
    ENDDO
    !
    !Compute prelogarithmic energy factor
    Efactor = ANISO_EFACTOR(b,A_kn,B_ijk,Dn,Pn)
    CALL ATOMSK_MSG(2101,(/'Kb²/4pi'/),(/Efactor/))
    !
    !
    !
  ELSE
    !Isotropic elasticity
    !Apply elastic displacements to all atoms
    DO i=1,SIZE(P,1)
      IF(IS_SELECTED(SELECT,i)) THEN
        IF( P(i,a1)>=pos(1)+0.1d0 .AND. (P(i,a2)-pos(2))*bsign<=0.d0 ) THEN
          P(i,a1) = P(i,a1)-DABS(b(a1))/2.d0
          IF( doshells ) THEN
            S(i,a1) = S(i,a1)-DABS(b(a1))/2.d0
          ENDIF
        ENDIF
      ENDIF
      disp = DISPEDGE(i,P(i,:),a1,a2,b(a1),nu,pos(1),pos(2))
      P(i,1:3) = P(i,1:3) + disp(:)
      !Apply displacements to shell
      IF( doshells ) THEN
        S(i,1:3) = S(i,1:3) + disp(:)
      ENDIF
      !
      IF( P(i,a1)>pos(1)+0.1d0 .AND. P(i,a2)*bsign>=bsign*pos(2) ) THEN
        P(i,a1) = P(i,a1)-DABS(b(a1))/2.d0
        IF( doshells ) THEN
          S(i,a1) = S(i,a1)-DABS(b(a1))/2.d0
        ENDIF
      ENDIF
      !
      !Compute stress
      CALL STRESSEDGE(P(i,1:3),a1,a2,bsign*VECLENGTH(b),mu,nu,pos(1),pos(2),sigma)
      !
      !Add stresses due to this disloc. to existing stresses in AUX
      AUX(i,sig1) = AUX(i,sig1) + sigma(1,1)
      AUX(i,sig2) = AUX(i,sig2) + sigma(2,2)
      AUX(i,sig3) = AUX(i,sig3) + sigma(3,3)
      AUX(i,sig6) = AUX(i,sig6) + sigma(1,2)
    ENDDO
    !
    !Compute prelogarithmic energy factor
    Efactor = ISO_EFACTOR(VECLENGTH(b),mu,nu,pi/4.d0)
    CALL ATOMSK_MSG(2101,(/'b²/4pi(1-nu)'/),(/Efactor/))
    !
  ENDIF  !Endif aniso
  !
  !Supercell was shrinked along a1: calculate the new length
  H(a1,a1) = H(a1,a1) - DABS(b(a1))/2.d0
  CALL ATOMSK_MSG(2064,(/"-b/2"/),(/DBLE(a1)/))
  !
  Ndisloc = 1
  !
  !
ELSEIF(disloctype=='mixed') THEN
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!  MIXED DISLOCATION  !!!!!!!!!!
  IF( aniso ) THEN
    !Anisotropic elasticity
    DO i=1,SIZE(P,1)
      IF(IS_SELECTED(SELECT,i)) THEN
        disp = ANISO_DISP(i,P(i,:),a1,a2,a3,pos(1),pos(2),A_kn,Dn,Pn)
        P(i,1:3) = P(i,1:3) + disp(:)
        IF( doshells ) THEN
          S(i,1:3) = S(i,1:3) + disp(:)
        ENDIF
        !
        !Compute stress
        CALL ANISO_STRESS(P(i,1:3),a1,a2,a3,pos(1),pos(2),A_kn,B_ijk,Dn,Pn,sigma)
        !
        !Add stresses due to this disloc. to existing stresses in AUX
        AUX(i,sig1) = AUX(i,sig1) + sigma(1,1)
        AUX(i,sig2) = AUX(i,sig2) + sigma(2,2)
        AUX(i,sig3) = AUX(i,sig3) + sigma(3,3)
        AUX(i,sig4) = AUX(i,sig4) + sigma(2,3)
        AUX(i,sig5) = AUX(i,sig5) + sigma(1,3)
        AUX(i,sig6) = AUX(i,sig6) + sigma(1,2)
      ENDIF
    ENDDO
    !
    !Compute prelogarithmic energy factor
    Efactor = ANISO_EFACTOR(b,A_kn,B_ijk,Dn,Pn)
    CALL ATOMSK_MSG(2101,(/'Kb²/4pi'/),(/Efactor/))
    !
    !
  ELSE
    !Isotropic elasticity: cannot build mixed dislocation
    CALL ATOMSK_MSG(2808,(/""/),(/0.d0/))
    nerr=nerr+1
    GOTO 1000
    !
  ENDIF
  !
  Ndisloc = 1
  !
  !
ELSEIF( disloctype=="loop" ) THEN
  !A dislocation loop must be introduced in the system
  !The direction normal to the loop is dislocline
  !pos(1:3) = position of the center of the loop; pos(4) = radius of the loop
  !
  SELECT CASE(dislocline)
  CASE('x','X')
    !Loop in (y,z) plane
    a1 = 2
    a2 = 3
    a3 = 1
  CASE('y','Y')
    a1 = 3
    a2 = 1
    a3 = 2
  CASE DEFAULT
    a1 = 1
    a2 = 2
    a3 = 3
  END SELECT
  !
  !Discretize loop into segments
  IF( pos(4) < 0.d0 ) THEN
    !"Negative" radius => user wants a square of side pos(4)
    ALLOCATE(xLoop(4,3))
    xLoop(:,:) = 0.d0
    xLoop(1,a1) = pos(a1) - pos(4)
    xLoop(1,a2) = pos(a2) - pos(4)
    xLoop(1,a3) = pos(a3)
     xLoop(2,a1) = pos(a1) + pos(4)
     xLoop(2,a2) = pos(a2) - pos(4)
     xLoop(2,a3) = pos(a3)
    xLoop(3,a1) = pos(a1) + pos(4)
    xLoop(3,a2) = pos(a2) + pos(4)
    xLoop(3,a3) = pos(a3)
     xLoop(4,a1) = pos(a1) - pos(4)
     xLoop(4,a2) = pos(a2) + pos(4)
     xLoop(4,a3) = pos(a3)
  ELSE
    !Positive pos(4) => Generate points that belong to a circle of radius pos(4)
    xLoop = LOOP_SEGMENTS( pos(1:3) , pos(4) , dislocline )
  ENDIF
  !
  IF( verbosity==4 ) THEN
    CALL ATOMSK_MSG(999,(/"Positions of points of the loop:"/),(/0.d0/))
    DO i=1,SIZE(xLoop,1)
      WRITE(msg,'(3(3X,f16.3))') xLoop(i,1), xLoop(i,2), xLoop(i,3)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
    !Write positions of the points in an XYZ file for visualization
    OPEN(UNIT=41,FILE="atomsk_dislocationloop.xyz",FORM="FORMATTED",STATUS="UNKNOWN")
    WRITE(41,*) SIZE(xLoop,1)
    WRITE(41,*) "# Points for dislocation loop generated by Atomsk"
    DO i=1,SIZE(xLoop,1)
      WRITE(41,'(a1,3X,3(3X,f16.3))') "H", xLoop(i,1), xLoop(i,2), xLoop(i,3)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
    WRITE(41,*) "supercell"
    DO i=1,3
      WRITE(41,'(3(3X,f16.3))') H(i,1), H(i,2), H(i,3)
    ENDDO
    CLOSE(41)
  ENDIF
  !
  !Now xLoop should contain the positions of points defining the loop
  !Check if some points are outside of the box
  CALL COUNT_OUTBOX(H,xLoop,i)
  IF(i>0) THEN
    !Some points of the loop are outside of the box => display a warning
    IF( i>=SIZE(xLoop,1) ) THEN
      !ALL points are outside of the box => something is probably very wrong
      CALL ATOMSK_MSG(2754,(/'all'/),(/0.d0/))
    ELSE
      !There are i points that are outside of the box
      CALL ATOMSK_MSG(2754,(/'loop'/),(/DBLE(i)/))
    ENDIF
  ENDIF
  !
  !Save current number of atoms
  u = SIZE(P,1)
  !
  IF(ALLOCATED(AUX) .AND. SIZE(AUX,1)>0) THEN
    doaux=.TRUE.
  ENDIF
  !
  !If Burgers vector has a component normal to the loop plane (=prismatic dislocation loop),
  !duplicate or remove atoms that belong to the loop
  IF( DABS(b(a3))>1.d-2 ) THEN
    IF( b(a3)>0.d0 ) THEN
      !The dislocation loop will separate atomic planes
      !=> new atoms must be inserted inside the loop
      !Get total system volume
      CALL VOLUME_PARA(H,V1)
      !Compute average system density
      tempreal = DBLE( SIZE(P,1)/V1 )
      !Compute area of disloc. loop
      V2 = 4.d0*pos(4)*pos(4)
      !Estimate number of atoms to duplicate = (average density)*(loop area)*(height of Burgers vector)
      !(may be over-estimated)
      n = CEILING( 1.1d0 * tempreal * V2 * VECLENGTH(b(:)) ) + 10
      !New atoms must be inserted inside the loop
      !Allocate memory for new atoms
      ALLOCATE( Q(u+n,4) )
      Q(:,:) = 0.d0
      IF( doaux ) THEN
        ALLOCATE( newAUX(u+n,SIZE(AUX,2)) )
        newAUX(:,:) = 0.d0
        !Copy data from AUX to newAUX
        DO i=1,SIZE(AUX,1)
          newAUX(:,1:SIZE(AUX,2)) = AUX(:,:)
        ENDDO
      ENDIF
      n = u
      !Copy atom positions and duplicated atoms in temporary array Q(:,:)
      DO i=1,SIZE(P,1)
        !Copy this atom into Q(:,:)
        Q(i,:) = P(i,:)
        !Increment r if atom is inside the loop
        r = 0
        IF( P(i,a3)<pos(a3)+b(a3)/2.d0 .AND. P(i,a3)>pos(a3)-b(a3)/2.d0 ) THEN
          IF( pos(4)<0.d0 ) THEN
            !"Negative" loop radius: square loop
            IF( P(i,a1)>pos(a1)-pos(4) .AND. P(i,a1)<pos(a1)+pos(4) .AND. &
              & P(i,a2)>pos(a2)-pos(4) .AND. P(i,a2)<pos(a2)+pos(4)       ) THEN
              r=r+1
            ENDIF
          ELSE
            !Positive loop radius: circular loop
            IF( VECLENGTH( P(i,1:3)-pos(1:3) ) < pos(4) ) THEN
              r=r+1
            ENDIF
          ENDIF
        ENDIF
        IF( r>0 ) THEN
          !Atom #i must be duplicated: add it at the end of atom list
          n = n+1
          Q(n,:) = P(i,:)
          Q(n,a3) = Q(n,a3) + b(a3)/2.d0
          IF(ALLOCATED(AUX) .AND. SIZE(AUX,1)>0) THEN
            newAUX(n,:) = AUX(i,:)
          ENDIF
        ENDIF
      ENDDO
      !Save final atom positions inside P(:,:)
      DEALLOCATE(P)
      ALLOCATE(P(n,4))
      IF(doaux) THEN
        DEALLOCATE(AUX)
        ALLOCATE( AUX(n,SIZE(newAUX,2)) )
      ENDIF
      DO i=1,n
        P(i,:) = Q(i,:)
        IF(doaux) THEN
          AUX(i,:) = newAUX(i,:)
        ENDIF
      ENDDO
      DEALLOCATE(Q)
      !
    ELSE !i.e. if b(a3)>0
      !The dislocation will bring atomic planes closer
      !=> atoms inside the loop must be deleted to avoid overlap
      ALLOCATE(Q(SIZE(P,1),4))
      Q(:,:) = P(:,:)
      n = 0
      !Delete atoms inside dislocation loop
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,r) REDUCTION(+:n)
      DO i=1,SIZE(P,1)
        !Increment r if atom is inside the loop
        r = 0
        IF( P(i,a3)<pos(a3)-b(a3)/2.d0 .AND. P(i,a3)>pos(a3)+b(a3)/2.d0 ) THEN
          IF( pos(4)<0.d0 ) THEN
            !"Negative" loop radius: square loop
            IF( P(i,a1)>pos(a1)-pos(4) .AND. P(i,a1)<pos(a1)+pos(4) .AND. &
              & P(i,a2)>pos(a2)-pos(4) .AND. P(i,a2)<pos(a2)+pos(4)       ) THEN
              r=r+1
            ENDIF
          ELSE
            !Positive loop radius: circular loop
            IF( VECLENGTH( P(i,1:3)-pos(1:3) ) < pos(4) ) THEN
              r=r+1
            ENDIF
          ENDIF
        ENDIF
        !
        IF( r==0 ) THEN
          !Atom will survive: increment number of atoms n
          n=n+1
        ELSE
          !Atom #i must be deleted: set its species to zero
          Q(i,:) = 0.d0
        ENDIF
      ENDDO
      !$OMP END PARALLEL DO
      !Save final atom positions inside P(:,:)
      DEALLOCATE(P)
      ALLOCATE(P(n,4))
      ALLOCATE(newAUX(n,4))
      n = 0
      DO i=1,SIZE(Q,1)
        IF( Q(i,4)>0.1d0 ) THEN
          n = n+1
          P(n,:) = Q(i,:)
          !Copy auxiliary properties (if any)
          IF(doaux) THEN
            newAUX(n,:) = AUX(i,:)
          ENDIF
        ENDIF
      ENDDO
      DEALLOCATE(Q)
      !Save auxiliary properties in AUX
      IF( doaux ) THEN
        DEALLOCATE(AUX)
        ALLOCATE( AUX(SIZE(newAUX,1),SIZE(newAUX,2)) )
        AUX(:,:) = newAUX(:,:)
        DEALLOCATE(newAUX)
      ENDIF
      !
    ENDIF
  ENDIF
  !
  !For each atom, compute its displacement due to the loop
  !(if atoms were duplicated just before, don't displace them)
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,disp)
  DO i=1,MIN(u,SIZE(P,1))
    disp(:) = LOOP_DISPLACEMENT( P(i,1:3) , -1.d0*b, nu, pos(1:3) , xLoop )
    P(i,1:3) = P(i,1:3) + disp(:)
  ENDDO
  !$OMP END PARALLEL DO
  !
  Ndisloc = 1
  !
  !
ELSEIF( disloctype=="file" .OR. disloctype=="array" ) THEN
  !Dislocation segments positions and Burgers vector are read from file "dislocline"
  !This text file must contain the keyword "dislocation" followed by the Burgers vector,
  !then on following lines the postions (x,y,z) through which passes the dislocation
  CALL CHECKFILE(dislocline,"read ")
  !
  OPEN(UNIT=35,FILE=dislocline,FORM="FORMATTED",STATUS="OLD")
  !First loop to count number of dislocation lines and max. number of segments
  i = 0  !counter for lines in the file
  r = 0  !counter for disloc.lines
  u = 0  !counter for disloc. segments
  DO
    READ(35,'(a128)',END=350,ERR=350) msg
    msg = TRIM(ADJUSTL(msg))
    i = i+1
    !Ignore empty lines and lines starting with "#"
    IF( LEN_TRIM(msg)>0 .AND. msg(1:1).NE."#" ) THEN
      IF( msg(1:6)=="disloc" ) THEN
        !Definition of a new dislocation
        IF( INDEX(msg,"box")>0 .OR. INDEX(msg,"BOX")>0 ) GOTO 802
        r = r+1
        u = 0
      ELSE
        !Try to read 3 real numbers
        IF( INDEX(msg,"box")>0 .OR. INDEX(msg,"BOX")>0 ) GOTO 802
        IF( SCAN(msg,"0123456789")==0 ) GOTO 349
        READ(msg,*,END=349,ERR=349) V1, V2, V3
        u = u+1
      ENDIF
    ENDIF
    349 CONTINUE
  ENDDO
  !
  350 CONTINUE
  !Finished counting lines and segments
  !Allocate array for storing Burgers vectors
  ALLOCATE( disarrayb(r,3) )
  disarrayb(:,:) = 0.d0
  !Allocate array for storing positions of points
  ALLOCATE( disarraypos(r,u,3) )
  disarraypos(:,:,1) = -1000.d0
  disarraypos(:,:,2) = -1001.d0
  disarraypos(:,:,3) = -1002.d0
  !
  !Go back to beginning of file and save data into arrays
  REWIND(35)
  i = 0  !counter for lines in the file
  r = 0  !counter for disloc.lines
  u = 0  !counter for disloc. segments
  DO
    READ(35,'(a128)',END=355,ERR=355) msg
    msg = TRIM(ADJUSTL(msg))
    i = i+1
    !Ignore empty lines and lines starting with "#"
    IF( LEN_TRIM(msg)>0 .AND. msg(1:1).NE."#" ) THEN
      IF( msg(1:6)=="disloc" ) THEN
        !Definition of a new dislocation
        IF( INDEX(msg,"box")>0 .OR. INDEX(msg,"BOX")>0 ) GOTO 802
        r = r+1
        u = 0
        j=SCAN(msg," ")
        READ(msg(j+1:),*,END=355,ERR=801) V1, V2, V3
        disarrayb(r,1) = V1
        disarrayb(r,2) = V2
        disarrayb(r,3) = V3
      ELSE
        !Try to read 3 real numbers
        IF( INDEX(msg,"box")>0 .OR. INDEX(msg,"BOX")>0 ) GOTO 802
        IF( SCAN(msg,"0123456789")==0 ) GOTO 354
        READ(msg,*,END=354,ERR=801) V1, V2, V3
        u = u+1
        disarraypos(r,u,1) = V1
        disarraypos(r,u,2) = V2
        disarraypos(r,u,3) = V3
      ENDIF
    ENDIF
    354 CONTINUE
  ENDDO
  !
  355 CONTINUE
  !Finished reading the file
  CLOSE(35)
  !
  !Insert dislocations into the system
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,k,n,r,u,V1,V2,V3,pos,disp) REDUCTION(+:Ndisloc)
  !Loop on all dislocation loops
  DO r=1,SIZE(disarrayb,1)
    IF( verbosity==4 ) THEN
      WRITE(msg,'(a14,i2,a8,3(f9.3,a3))') "Dislocation # ", r, ": b = [ ", &
            & disarrayb(r,1), " , ", disarrayb(r,2), " , ", disarrayb(r,3), " ]."
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      DO u=1,SIZE(disarraypos,2)
        IF( NINT(disarraypos(r,u,1)).NE.-1000 .OR. NINT(disarraypos(r,u,2)).NE.-1001 &
          & .OR. NINT(disarraypos(r,u,3)).NE.-1002 ) THEN
          WRITE(msg,'(6X,3(f9.3,2X))') disarraypos(r,u,1:3)
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        ENDIF
      ENDDO
    ENDIF
    !
    !Check if some points of the loop are out of the box
    CALL COUNT_OUTBOX(H,disarraypos(r,:,1:3),u)
    IF(u>0) THEN
      !Some points of the loop are outside of the box => display a warning
      IF( u>=SIZE(disarraypos(r,u,1:3),1) ) THEN
        !ALL points are outside of the box => something is probably very wrong
        CALL ATOMSK_MSG(2754,(/'all'/),(/0.d0/))
      ELSE
        !There are u points that are outside of the box
        CALL ATOMSK_MSG(2754,(/'loop'/),(/DBLE(u)/))
      ENDIF
    ENDIF
    !
    !Check if initial and final points are not too close
    n = SIZE(disarraypos,2)
    IF( VECLENGTH( disarraypos(r,1,1:3) - disarraypos(r,n,1:3) ) < 1.d-3 ) THEN
      disarraypos(r,n,1) = -1000.d0
      disarraypos(r,n,2) = -1001.d0
      disarraypos(r,n,3) = -1002.d0
    ENDIF
    !
    !Compute position of center of this loop, save it in pos(:)
    k = 0
    pos(:) = 0.d0
    DO u=1,SIZE(disarraypos,2)
      IF( NINT(disarraypos(r,u,1))==-1000 .AND. NINT(disarraypos(r,u,2))==-1001 &
        & .AND. NINT(disarraypos(r,u,3))==-1002 ) THEN
        EXIT
      ELSE
        k = k+1
        pos(1:3) = pos(1:3) + disarraypos(r,u,1:3)
      ENDIF
    ENDDO
    pos(1:3) = pos(1:3) / DBLE(k)
    !
    !Loop on all atoms
    DO i=1,SIZE(P,1)
      !Compute displacements of atom #i due to current dislocation segments
      disp(:) = 0.d0
      disp(:) = disp(:) + LOOP_DISPLACEMENT( P(i,1:3) , disarrayb(r,:), 0.33d0, pos(1:3) , disarraypos(r,1:k,:) )
      !Apply displacement to atom #i
      P(i,1:3) = P(i,1:3) + disp(:)
    ENDDO
    !
    !TO BE IMPLEMENTED: if Burgers vector is not in cut plane, duplicate atoms to fill the void
    !Detect atoms above and below cut plane, duplicate them, translate them by b(:)
    !DO i=1,SIZE(P,1)
    
    !ENDDO
    !
    !Increment dislocation counter
    Ndisloc = Ndisloc + 1
    !
  ENDDO
  !$OMP END PARALLEL DO
  !
  IF(ALLOCATED(disarraypos)) DEALLOCATE(disarraypos)
  !
  !
ELSE
  !Unknown dislocation type
  nerr=nerr+1
  GOTO 1000
ENDIF
!
!Dislocation(s) was successfully inserted: display message on screen
CALL ATOMSK_MSG(2065,(/''/),(/DBLE(Ndisloc)/))
!
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
GOTO 1000
!
801 CONTINUE
CALL ATOMSK_MSG(807,(/''/),(/DBLE(i)/))
nerr = nerr+1
GOTO 1000
!
802 CONTINUE
CALL ATOMSK_MSG(2822,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
1000 CONTINUE
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(disarrayb)) DEALLOCATE(disarrayb)
!
!
END SUBROUTINE DISLOC_XYZ
!
!
!
END MODULE dislocation
