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
!* Last modification: P. Hirel - 07 May 2019                                      *
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
INTEGER:: dimens  !dimensionality of the system: 3-D, 2-D, 1-D or 0-D
INTEGER:: i, j, k, l
INTEGER:: NP  !counter for total number of particles
INTEGER:: sgnumber  !space group number
INTEGER:: symconv, symtype, symsetting
REAL(dp):: a, b, c, alpha, beta, gamma !supercell (conventional notation)
REAL(dp),DIMENSION(6):: values   !values read
REAL(dp),DIMENSION(3,3),INTENT(OUT):: H    !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: P    !Atom positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S    !Shell positions (dummy array, not used here)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: AUX  !auxiliary properties
!
!
!Initialize variables
dimens = 3  !default: 3-D periodic
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
  ELSEIF( line(1:7)=='CRYSTAL' ) THEN
    !
    !Read properties for space group:
    !- symconv = convention for s.g. identification, 0=number, 1=H-M symbol
    !- symtype = type of cell, 0=hexagonal (a,c), 1=rhombohedral (a, alpha)
    !- symsetting = setting for origin of crystal reference frame (ignored here)
    READ(30,*,ERR=190,END=190) symconv, symtype, symsetting
    !
    !Read space group number or Hermann-Mauguin symbol
    READ(30,'(a128)',ERR=190,END=190) line
    !
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
    !
    !Set default lattice constants
    a = 4.d0
    b = a
    c = a
    alpha = 90.d0
    beta  = 90.d0
    gamma = 90.d0
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
    IF( k==2 .AND. symtype==1 ) THEN
      !Rhombohedral lattice: actually read a and gamma
      gamma = c
      c = a
    ENDIF
    !Convert angles into radians
    alpha = DEG2RAD(alpha)
    beta  = DEG2RAD(beta)
    gamma = DEG2RAD(gamma)
    !Convert lattice parameters into box vectors
    CALL CONVMAT(a,b,c,alpha,beta,gamma,H)
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
  ELSEIF( line(1:4)=='SLAB' ) THEN
    !dimens = 2
    !
  ELSEIF( line(1:7)=='POLYMER' ) THEN
    dimens = 1
    !
  ELSEIF( line(1:8)=='MOLECULE' ) THEN
    dimens = 0
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
!
IF( dimens==3 ) THEN
  !All coordinates are fractional: convert to cartesian
  CALL FRAC2CART(P,H)
ENDIF
!
!If all cell vectors are set to 1, it means they are not set => reset them to zero
IF( DABS(VECLENGTH(H(1,:)))-1.d0<1.d-12 .AND.         &
  & DABS(VECLENGTH(H(2,:)))-1.d0<1.d-12 .AND.         &
  & DABS(VECLENGTH(H(3,:)))-1.d0<1.d-12      ) THEN
  H(:,:) = 0.d0
ENDIF
!Apply symmetry operations
IF( sgnumber > 1 .AND. sgnumber <= 230 ) THEN
  CALL ATOMSK_MSG(1004,(/""/),(/0.d0/))
  ALLOCATE(symops_trf(sgnumber,SIZE(symop_list)))
  !Initialize symmetry operations
  CALL SYMOPS_INIT()
  !Apply symmetry operations
  WRITE(msg,'(i3)') sgnumber
  CALL SG_APPLY_SYMOPS(msg,H,P,S,AUXNAMES,AUX)  
ENDIF
GOTO 1000
!
!
!
800 CONTINUE
801 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
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
