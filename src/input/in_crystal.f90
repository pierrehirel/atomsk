MODULE in_crystal
!
!
!**********************************************************************************
!*  IN_CRYSTAL                                                                    *
!**********************************************************************************
!* This module reads a file in the CRYSTAL format.                                *
!* The CRYSTAL format is officially described here:                               *
!*     http://www.crystal.unito.it/                                               *
!**********************************************************************************
!* (C) May 2019 - Pierre Hirel                                                    *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 31 May 2019                                      *
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
USE resize
USE subroutines
USE symops
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE READ_CRYSTAL(inputfile,H,P,S,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: line
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: symop_list  !list of symmetry operations
LOGICAL:: new
INTEGER:: dimens  !dimensionality of the system: 3-D, 2-D, 1-D or 0-D
INTEGER:: i, j, k, l, m
INTEGER:: lminmax  !min and max values for loops when orienting cubic systems
INTEGER:: NP  !counter for total number of particles
INTEGER:: sgnumber  !space group number
INTEGER:: symconv, symtype, symsetting
REAL(dp):: a, b, c, alpha, beta, gamma !supercell (conventional notation)
REAL(dp),DIMENSION(6):: values   !values read
REAL(dp),DIMENSION(1,4):: tempP  !temporary position
REAL(dp),DIMENSION(3,3),INTENT(OUT):: H    !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT    !crystal orientation (defined with SLABCUT)
REAL(dp),DIMENSION(3,3):: uv        !unit vectors corresponding to new orientation ORIENT(:,:)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: P    !Atom positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S    !Shell positions (dummy array, not used here)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: AUX  !auxiliary properties
!
!
!Initialize variables
dimens = 0
NP = 0
sgnumber = 0
 a = 1.d0
 b = 1.d0
 c = 1.d0
 alpha = 90.d0
 beta = 90.d0
 gamma = 90.d0
H(:,:) = 0.d0
!
msg = 'entering READ_CRYSTAL'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
ALLOCATE(comment(1))
 comment(1) = ""
!
!
!
100 CONTINUE
OPEN(UNIT=30,FILE=inputfile,STATUS='UNKNOWN',ERR=800)
!
i=0
j=0
DO
  i=i+1
  READ(30,'(a128)',ERR=200,END=200) line
  line = ADJUSTL(line)
  !
  IF( i==1 ) THEN
    !First line is a comment
    comment(1) = "# "//TRIM(line)
    !
  ELSEIF( dimens==0 .AND. ( line(1:7)=='CRYSTAL' .OR. line(1:5)=='SLAB ' &
        & .OR. line(1:7)=='POLYMER' .OR. line(1:8)=='MOLECULE' )         ) THEN
    !
    IF( line(1:7)=='CRYSTAL' ) THEN
      dimens = 3
      !Read properties for space group:
      !- symconv = convention for s.g. identification, 0=number, 1=H-M symbol
      !- symtype = type of cell, 0=hexagonal (a,c), 1=rhombohedral (a, alpha)
      !- symsetting = setting for origin of crystal reference frame (ignored here)
      READ(30,*,ERR=190,END=190) symconv, symtype, symsetting
      IF( symsetting>1 ) THEN
        READ(30,'(a128)',ERR=190,END=190) line
      ENDIF
    ELSEIF( line(1:4)=='SLAB' ) THEN
      dimens = 2
    ELSEIF( line(1:7)=='POLYMER' ) THEN
      dimens = 1
    ELSEIF( line(1:8)=='MOLECULE' ) THEN
      dimens = 0
    ENDIF
    !
    !Read space group number or Hermann-Mauguin symbol
    READ(30,'(a128)',ERR=190,END=190) line
    !
    IF( dimens==3 ) THEN
      IF( symconv==0 ) THEN
        !Read space group number
        READ(line,*,ERR=200,END=200) sgnumber
        !Verify that it is a valid space group number; otherwise reset it to zero
        IF( sgnumber <1 .OR. sgnumber > 230 ) THEN
          sgnumber = 0
        ENDIF
      ELSEIF( symconv==1 ) THEN
        !Hermann-Mauguin notation
        !Remove spaces from the string
        DO j=1,LEN_TRIM(line)
          IF(line(j:j)==' ') THEN
            line(j:) = TRIM(line(j+1:))
          ENDIF
        ENDDO
        !Transform H-M symbol into a space group number
        CALL SG_NAMGETNUM(line,sgnumber)
      ELSE
        !Rubbish
      ENDIF
    ENDIF
    !
    !Set default lattice constants
    a = 1.d0
    b = a
    c = a
    alpha = 90.d0
    beta  = 90.d0
    gamma = 90.d0
    IF( dimens>0 ) THEN
      !Read actual lattice constant(s)
      READ(30,'(a128)',ERR=190,END=190) line
      READ(line,*,ERR=800,END=800) values(1)
      a = values(1)
      b = values(1)
      c = values(1)
      !Try to read other lattice parameters
      k = 1
      READ(line,*,ERR=120,END=120) values(1), values(2)
      a = values(1)
      c = values(2)
      k = 2
      READ(line,*,ERR=120,END=120) values(1), values(2), values(3)
      a = values(1)
      b = values(2)
      c = values(3)
      k = 3
      READ(line,*,ERR=120,END=120) values(1), values(2), values(3), values(4)
      a = values(1)
      b = values(2)
      c = values(3)
      gamma = values(4)
      IF( DABS(a-c)<1.d-12 ) THEN
        !Actually read beta angle
        beta = gamma
        gamma = 90.d0
      ELSEIF( DABS(b-c)<1.d-12 ) THEN
        !Actually read alpha angle
        alpha = gamma
        gamma = 90.d0
      ENDIF
      k = 4
      READ(line,*,ERR=120,END=120) values(1), values(2), values(3), values(4), values(5)
      !Not supposed to happen, but accept it as it is
      a = values(1)
      b = values(2)
      c = values(3)
      alpha = values(4)
      beta  = values(5)
      k = 5
      READ(line,*,ERR=120,END=120) values(1), values(2), values(3), values(4), values(5), values(6)
      a = values(1)
      b = values(2)
      c = values(3)
      alpha = values(4)
      beta  = values(5)
      gamma = values(6)
      k = 6
      !
      120 CONTINUE
      IF( dimens<3 ) THEN
        IF( k==3 ) THEN
          !Actually read a, b, and gamma
          gamma = c
          c = a
        ELSEIF( k==2 ) THEN
          !Actually read a and b
          b = c
          c = a
        ENDIF
      ELSEIF( k==2 .AND. symtype==1 ) THEN
        !Rhombohedral lattice: actually read a and gamma
        gamma = c
        c = a
      ENDIF
      !
      !Convert angles into radians
      alpha = DEG2RAD(alpha)
      beta  = DEG2RAD(beta)
      gamma = DEG2RAD(gamma)
      !Convert lattice parameters into box vectors
      CALL CONVMAT(a,b,c,alpha,beta,gamma,H)
      !
    ELSE  !i.e. dimens==0
      !Molecule: a bounding box will be constructed
      !after reading the atom positions (see below)
    ENDIF
    !
    !Read number of atoms to follow
    READ(30,*,ERR=800,END=800) j
    IF(j<=0) GOTO 800
    !
    !Allocate memory to store atoms
    IF( ALLOCATED(P) ) THEN
      ALLOCATE(S(SIZE(P,1)+j,4))
      S(:,:) = 0.d0
      DO k=1,NP
        S(k,:) = P(k,:)
      ENDDO
      DEALLOCATE(P)
      ALLOCATE(P(SIZE(S,1),4))
      P(:,:) = S(:,:)
      DEALLOCATE(S)
    ELSE
      ALLOCATE(P(j,4))
      P(:,:) = 0.d0
    ENDIF
    !
    !Read atom species and positions
    DO k=1,j
      NP = NP+1
      READ(30,*,ERR=190,END=190) P(NP,4), P(NP,1), P(NP,2), P(NP,3)
    ENDDO
    !
    !
    IF( dimens==3 ) THEN
      !3-D system ("CRYSTAL"): all coordinates are fractional, convert them to Cartesian
      CALL FRAC2CART(P,H)
      !
      !Apply symmetry operations (only 3-D periodic crystal supported)
      IF( sgnumber > 1 .AND. sgnumber <= 230 ) THEN
        CALL ATOMSK_MSG(1004,(/""/),(/0.d0/))
        ALLOCATE(symops_trf(sgnumber,SIZE(symop_list)))
        !Initialize symmetry operations
        CALL SYMOPS_INIT()
        !Apply symmetry operations
        WRITE(msg,'(i3)') sgnumber
        CALL SG_APPLY_SYMOPS(msg,H,P,S,AUXNAMES,AUX)  
      ENDIF
      !
    ELSEIF( dimens==2 ) THEN
      !2-D system ("SLAB"): only X and Y are fractional, Z is already Cartesian
      DO i=1,SIZE(P,1)
        a = P(i,1)
        b = P(i,2)
        c = P(i,3)
        P(i,1) = a*H(1,1) + b*H(2,1) + c*H(3,1)
        P(i,2) = a*H(1,2) + b*H(2,2) + c*H(3,2)
      ENDDO
      !
    ELSEIF( dimens==1 ) THEN
      !2-D system ("POLYMER"): only X is fractional, Y and Z are already Cartesian
      DO i=1,SIZE(P,1)
        a = P(i,1)
        b = P(i,2)
        c = P(i,3)
        P(i,1) = a*H(1,1) + b*H(2,1) + c*H(3,1)
      ENDDO
      !
    ELSEIF( dimens==0 ) THEN
      !Molecule: construct a bounding box that encloses all atoms
      H(1,1) = DABS( MAXVAL(P(:,1)) - MINVAL(P(:,1)) )
      H(2,2) = DABS( MAXVAL(P(:,2)) - MINVAL(P(:,2)) )
      H(3,3) = DABS( MAXVAL(P(:,3)) - MINVAL(P(:,3)) )
    ENDIF
    !
    !If all cell vectors are set to 1, it means they are not set => reset them to zero
    IF( DABS(VECLENGTH(H(1,:)))-1.d0<1.d-12 .AND.         &
      & DABS(VECLENGTH(H(2,:)))-1.d0<1.d-12 .AND.         &
      & DABS(VECLENGTH(H(3,:)))-1.d0<1.d-12      ) THEN
      H(:,:) = 0.d0
    ENDIF
    !
    !
  ELSEIF( line(1:5)=='SLAB ' .OR. line(1:7)=='SLABCUT' ) THEN
    !Cut a slab along a specific direction
    IF( .NOT.ALLOCATED(P) ) THEN
      !Cannot orient the crystal if there is not crystal
      nwarn = nwarn+1
    ELSE
      !Read crystal direction (3 integer Miller indices) along the Z axis
      READ(30,'(a128)',ERR=190,END=190) line
      READ(line,*,ERR=190,END=190) ORIENT(3,1), ORIENT(3,2), ORIENT(3,3)
      !Determine best matching crystal axis along X (shortest vector)
      DO k=0,20
        DO j=0,20
          DO i=0,20
            DO l=1,-1,-2
              IF( i.NE.0 .OR. j.NE.0 .OR. k.NE.0 ) THEN
                m = l*i*NINT(ORIENT(3,1)) + j*NINT(ORIENT(3,2)) + k*NINT(ORIENT(3,3))
                IF( m==0 ) THEN
                  ORIENT(1,1) = i*l
                  ORIENT(1,2) = j
                  ORIENT(1,3) = k
                  GOTO 150
                ELSE
                  m = i*NINT(ORIENT(3,1)) + l*j*NINT(ORIENT(3,2)) + k*NINT(ORIENT(3,3))
                  IF( m==0 ) THEN
                    ORIENT(1,1) = i
                    ORIENT(1,2) = j*l
                    ORIENT(1,3) = k
                    GOTO 150
                  ELSE
                    m = i*NINT(ORIENT(3,1)) + j*NINT(ORIENT(3,2)) + l*k*NINT(ORIENT(3,3))
                    IF( m==0 ) THEN
                      ORIENT(1,1) = i
                      ORIENT(1,2) = j
                      ORIENT(1,3) = k*l
                      GOTO 150
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      150 CONTINUE
      !Determine crystal direction along Y (normal to X and Z)
      ORIENT(2,:) = CROSS_PRODUCT( ORIENT(3,:) , ORIENT(1,:) )
      !
      WRITE(msg,*) "SLABCUT: required crystal orientation:"
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      WRITE(msg,'(3f16.6)') ORIENT(1,:)
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      WRITE(msg,'(3f16.6)') ORIENT(2,:)
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      WRITE(msg,'(3f16.6)') ORIENT(3,:)
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      !
      !Rotate crystal
      !
      !Set lminmax = 10 * largest value in ORIENT
      lminmax = 10 * NINT(MAXVAL(DABS(ORIENT(:,:))))
      !
      !The oriented unit cell vectors are defined by the Miller indices
      DO i=1,3
        !Set box vector
        uv(i,:) = ( ORIENT(i,1)*H(1,:) + ORIENT(i,2)*H(2,:) + ORIENT(i,3)*H(3,:) ) !/ x
      ENDDO
      !
      WRITE(msg,*) "Oriented unit cell vectors:"
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      WRITE(msg,'(3f16.6)') uv(1,:)
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      WRITE(msg,'(3f16.6)') uv(2,:)
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      WRITE(msg,'(3f16.6)') uv(3,:)
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      !
      !For each atom in the unit cell H(:,:), keep only periodic replica that are inside the uv(:)
      !and store it in S(:,:)
      WRITE(msg,*) "Duplicating atoms inside oriented unit cell..."
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      !Estimate new number of particles NP by comparing volumes of old and new unit cells
      NP = 1.2d0*CEILING( SIZE(P,1) * DABS( DABS(uv(1,1)*uv(2,2)*uv(3,3)) / &
          & DABS(VECLENGTH(H(1,:))*VECLENGTH(H(2,:))*VECLENGTH(H(3,:))) ) ) + 20
      WRITE(msg,*) "Estimated new number of atoms : ", NP
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      IF(ALLOCATED(S)) DEALLOCATE(S)
      ALLOCATE( S(NP,4) )
      S(:,:) = 0.d0
      !
      !Loop over all replica in a wide range
      NP = 0
      DO i=1,SIZE(P,1)
        DO l=-lminmax,lminmax
          DO k=-lminmax,lminmax
            DO j=-lminmax,lminmax
              !Compute cartesian position of this replica
              tempP(1,1) = P(i,1) + DBLE(j)*H(1,1) + DBLE(k)*H(2,1) + DBLE(l)*H(3,1)
              tempP(1,2) = P(i,2) + DBLE(j)*H(1,2) + DBLE(k)*H(2,2) + DBLE(l)*H(3,2)
              tempP(1,3) = P(i,3) + DBLE(j)*H(1,3) + DBLE(k)*H(2,3) + DBLE(l)*H(3,3)
              tempP(1,4) = P(i,4)
              CALL CART2FRAC(tempP,uv)
              IF( tempP(1,1)>=-1.d-12 .AND. tempP(1,1)<1.d0-1.d-12 .AND.             &
                & tempP(1,2)>=-1.d-12 .AND. tempP(1,2)<1.d0-1.d-12 .AND.             &
                & tempP(1,3)>=-1.d-12 .AND. tempP(1,3)<1.d0-1.d-12       ) THEN
                !This replica is inside the new cell, mark it as new
                new = .TRUE.
                !Verify that its position is different from all previous atoms
                DO m=1,NP
                  IF( DABS( VECLENGTH(tempP(1,1:3)-S(m,1:3)) )<1.d-6 ) THEN
                    new = .FALSE.
                  ENDIF
                ENDDO
                !
                IF( new ) THEN
                  NP = NP+1
                  IF(NP>SIZE(S,1)) THEN
                    !Resize array S
                    CALL RESIZE_DBLEARRAY2(S,SIZE(S,1)+10,SIZE(S,2))
                  ENDIF
                  S(NP,:) = tempP(1,:)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      !
      !Replace old P with the new S
      IF(ALLOCATED(P)) DEALLOCATE(P)
      ALLOCATE(P(NP,4))
      DO i=1,NP
        P(i,:) = S(i,:)
      ENDDO
      IF(ALLOCATED(S)) DEALLOCATE(S)
      WRITE(msg,*) "new NP in oriented cell:", NP
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      !
      !Replace old H with the new (oriented) cell vectors
      !Align 1st box vector with Cartesian X axis
      !Place 2nd box vector in the XY plane
      H(:,:) = 0.d0
      a = VECLENGTH(uv(1,:))
      b = VECLENGTH(uv(2,:))
      c = VECLENGTH(uv(3,:))
      alpha = ANGVEC(uv(2,:),uv(3,:))
      beta  = ANGVEC(uv(3,:),uv(1,:))
      gamma = ANGVEC(uv(1,:),uv(2,:))
      !Then convert this conventional notation into lower-triangular matrix H
      CALL CONVMAT(a,b,c,alpha,beta,gamma,H)
      !Transform atom positions in Cartesian coordinates
      CALL FRAC2CART(P,H)
      !
      !Read label of surface layer (ignored here), and number of slabs
      READ(30,'(a128)',ERR=190,END=190) line
      READ(line,*,ERR=190,END=190) j, k
      !Duplicate slab accordingly along the Z axis
      IF(ALLOCATED(S)) DEALLOCATE(S)
      NP = k*SIZE(P,1)
      ALLOCATE(S(NP,4))
      S(:,:) = 0.d0
      NP = 0
      DO j=0,k-1
        DO i=1,SIZE(P,1)
          NP = NP+1
          S(NP,:) = P(i,:)
          S(NP,3) = P(i,3) + DBLE(j)*H(3,3)
        ENDDO
      ENDDO
      DEALLOCATE(P)
      ALLOCATE(P(SIZE(S,1),4))
      P(:,:) = S(:,:)
      DEALLOCATE(S)
      H(3,:) = k*H(3,:)
    ENDIF
    !
  ELSEIF( line(1:8)=='ATOMSUBS' ) THEN
    !Read number of substitutions to perform
    READ(30,*,ERR=190,END=190) k
    DO j=1,k
      READ(30,*,ERR=190,END=190) l, a
      P(l,4) = a
    ENDDO
    !
  ELSEIF( line(1:8)=='ATOMDISP' ) THEN
    !Read number of displacements to perform
    READ(30,*,ERR=190,END=190) k
    DO j=1,k
      READ(30,*,ERR=190,END=190) l, a, b, c
      P(l,1) = P(l,1) + a
      P(l,2) = P(l,2) + b
      P(l,3) = P(l,3) + c
    ENDDO
    !
  ELSEIF( line(1:8)=='ATOMINSE' ) THEN
    !Read number of atoms to insert
    READ(30,*,ERR=190,END=190) k
    !Increase size of array P
    IF( ALLOCATED(P) ) THEN
      ALLOCATE(S(SIZE(P,1)+j,4))
      S(:,:) = 0.d0
      DO k=1,NP
        S(k,:) = P(k,:)
      ENDDO
      DEALLOCATE(P)
      ALLOCATE(P(SIZE(S,1),4))
      P(:,:) = S(:,:)
      DEALLOCATE(S)
    ELSE
      ALLOCATE(P(j,4))
      P(:,:) = 0.d0
    ENDIF
    !Read positions of new atoms
    DO j=1,k
      NP = NP+1
      READ(30,*,ERR=190,END=190) P(NP,4), P(NP,1), P(NP,2), P(NP,3)
      P(l,4) = a
    ENDDO
    !
  ELSEIF( line(1:8)=='ATOMREMO' ) THEN
    IF( ALLOCATED(P).AND.SIZE(P,1)>1 ) THEN
      !Read number of atoms to remove
      READ(30,*,ERR=190,END=190) k
      !Read labels of atoms to remove
      DO j=1,k
        READ(30,*,ERR=190,END=190) l
        IF(l>0 .AND. l<=SIZE(P,1) ) THEN
          NP = NP-1
          DO m=l,SIZE(P,1)-1
            P(m,:) = P(m+1,:)
          ENDDO
        ELSE
          !Index l is out of range
          nwarn = nwarn+1
          CALL ATOMSK_MSG(2742,(/""/),(/DBLE(l)/))
        ENDIF
      ENDDO
      IF( NP<SIZE(P,1) ) THEN
        !Decrease size of array P
        IF( ALLOCATED(P) ) THEN
          ALLOCATE(S(SIZE(P,1)-k,4))
          S(:,:) = 0.d0
          DO k=1,NP
            S(k,:) = P(k,:)
          ENDDO
          DEALLOCATE(P)
          ALLOCATE(P(SIZE(S,1),4))
          P(:,:) = S(:,:)
          DEALLOCATE(S)
        ENDIF
      ENDIF
      !
    ELSE  !i.e. P is not allocated
      nwarn = nwarn+1
      CALL ATOMSK_MSG(1709,(/"ATOMINSE    ","remove atoms"/),(/0.d0/))
    ENDIF
    !
    !
  ELSEIF( line(1:3)=='END' ) THEN
    !End of geometry section defining the crystal
    EXIT
    !
  ENDIF
  !
  190 CONTINUE
  !
ENDDO
!
!
!
200 CONTINUE
CLOSE(30)
!At this point, if array P(:,:) is not allocated we are in trouble
IF( .NOT.ALLOCATED(P) ) GOTO 800
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(807,(/''/),(/DBLE(i)/))
nerr = nerr+1
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE READ_CRYSTAL
!
END MODULE in_crystal
