MODULE in_str
!
!
!**********************************************************************************
!*  IN_STR                                                                        *
!**********************************************************************************
!* This module reads files in the PDFFIT structure format (*.str or *.stru).      *
!* The STR format is described on the Website of PDFFIT:                          *
!*   https://web.pa.msu.edu/cmp/billinge-group/programs/discus/pdffit.html        *
!**********************************************************************************
!* (C) February 2018 - Pierre Hirel                                               *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 12 Feb. 2018                                     *
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
SUBROUTINE READ_STR(inputfile,H,P,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=32):: spcgr  !name of space group
CHARACTER(LEN=128):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: pdffit !if .TRUE. format is pdffit, if .FALSE. format is DISCUS
INTEGER:: i, m, n, o
INTEGER:: NP   ! number of atoms
INTEGER:: Nsym        !number of symmetry operations
INTEGER:: sgroupnum  !space group number
INTEGER:: Nx, Ny, Nz !number of unit cells along X, Y, Z
INTEGER:: strlength
REAL(dp):: a0 !scaling factor
REAL(dp):: a, b, c, alpha, beta, gamma !supercell (conventional notation)
REAL(dp),DIMENSION(3,3):: G, Gi   !Metric tensor and its inverse
REAL(dp),DIMENSION(3,3),INTENT(OUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, Q, S
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
spcgr = ""
pdffit = .FALSE.
i=0
NP=0
Nx = 1
Ny = 1
Nz = 1
Nsym = 0
sgroupnum = 0
a0 = 1.d0  !scaling factor=1 by default
a=0.d0
b=0.d0
c=0.d0
alpha=0.d0
beta=0.d0
gamma=0.d0
IF (ALLOCATED(P)) DEALLOCATE(P)
IF (ALLOCATED(S)) DEALLOCATE(S)
IF (ALLOCATED(AUX)) DEALLOCATE(AUX)
IF (ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF (ALLOCATED(comment)) DEALLOCATE(comment)
ALLOCATE(comment(1))
 comment(1) = ""
!
msg = 'entering READ_STL'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
!
100 CONTINUE
OPEN(UNIT=30,FILE=inputfile,STATUS='UNKNOWN',ERR=800)
!
!Go back to beginning of file and count number of lines
REWIND(30)
i=1  ! init line counter
DO
  READ(30,'(a128)',ERR=207,END=207) temp
  temp = ADJUSTL(temp)
  !
  !Remove commas from the line
  strlength = SCAN(temp,',')
  DO WHILE(strlength>0)
    temp(strlength:strlength) = " "
    strlength = SCAN(temp,',')
  ENDDO
  !
  IF( temp(1:5)=="title" .OR. temp(1:5)=="TITLE" ) THEN
    !Save title as a comment
    comment(1) = "# "//TRIM(ADJUSTL(temp(6:)))
  ELSEIF( temp(1:6)=="format" .OR. temp(1:6)=="FORMAT" ) THEN
    !Read file format
    msg = TRIM(ADJUSTL(temp(7:)))
    IF( msg(1:6)=="pdffit" .OR. msg(1:6)=="PDFFIT" ) THEN
      pdffit = .TRUE.
    ENDIF
  ELSEIF( temp(1:5)=="scale" .OR. temp(1:5)=="SCALE" ) THEN
    !Save scaling factor (will be applied at end of this routine)
    READ(temp(6:),*,ERR=207,END=207) a0
  ELSEIF( temp(1:5)=="spcgr" .OR. temp(1:5)=="SPCGR" ) THEN
    !Save space group name (will be applied at end of this routine)
    temp = ADJUSTL(temp(6:))
    spcgr = temp(1:32)
  ELSEIF( temp(1:4)=="cell" .OR. temp(1:4)=="CELL" ) THEN
    !Read cell parameters (conventional notation, angles in degrees)
    READ(temp(5:),*,ERR=207,END=207) a, b, c, alpha, beta, gamma
    !Convert cell angles from degree to radian
    alpha = DEG2RAD(alpha)
    beta  = DEG2RAD(beta)
    gamma = DEG2RAD(gamma)
    !Save cell vectors into H(:,:)
    CALL CONVMAT(a,b,c,alpha,beta,gamma,H)
  ELSEIF( temp(1:5)=="ncell" .OR. temp(1:5)=="NCELL" ) THEN
    READ(temp(6:),*,ERR=207,END=207) Nx, Ny, Nz, a
    IF( a > NATOMS_MAX ) THEN
      nerr = nerr+1
      CALL ATOMSK_MSG(821,(/""/),(/a/))
      GOTO 1000
    ENDIF
    NP = NINT(a)
  ELSEIF( temp(1:5)=="atoms" .OR. temp(1:5)=="ATOMS" ) THEN
    IF( NP==0 ) THEN
      !Number of atoms is unknown => count them
      NP = 0
      DO
        i = i+1
        READ(30,'(a128)',ERR=110,END=110) temp
        temp = ADJUSTL(temp)
        species = temp(1:2)
        !Convert species to proper casing
        species(1:1) = StrUpCase(species(1:1))
        species(2:2) = StrDnCase(species(2:2))
        !Check if it is a recognizable species
        CALL ATOMNUMBER(species,a)
        IF( a>=0.9 ) THEN
          !This was an atom species => increment NP
          NP = NP+1
        ELSE
          !Atom species missing at beginning of that line
        ENDIF
      ENDDO
    ELSE
      !Number of atoms is known, this was the last keyword
      GOTO 110
    ENDIF
  ENDIF
  !
  i = i+1
  !
ENDDO
!
110 CONTINUE
WRITE(msg,*) 'Counted atoms, NP = ', NP
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
IF (NP<=0) GOTO 800 ! No atoms, -> error
WRITE(msg,*) 'Scaling factor a0 = ', a0
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,*) 'File format is pdffit: ', pdffit
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
! Allocate the atom data P and auxiliary data AUX and AUXNAMES
ALLOCATE(P(NP,4) , STAT=i)
IF( i>0 ) THEN
  ! Allocation failed (not enough memory)
  nerr = nerr+1
  CALL ATOMSK_MSG(819,(/''/),(/DBLE(NP)/))
  GOTO 1000
ENDIF
P = 0.d0
ALLOCATE(AUX(NP,8),AUXNAMES(8))
AUX = 0.d0
AUXNAMES(1) = "occ"
AUXNAMES(2) = "U11"
AUXNAMES(3) = "U22"
AUXNAMES(4) = "U33"
AUXNAMES(5) = "U12"
AUXNAMES(6) = "U13"
AUXNAMES(7) = "U23"
AUXNAMES(8) = "Biso"
!
IF( .NOT.pdffit ) THEN
  !Set metric tensor
  G(1,1) = DOT_PRODUCT( H(1,:) , H(1,:) )
  G(2,2) = DOT_PRODUCT( H(2,:) , H(2,:) )
  G(3,3) = DOT_PRODUCT( H(3,:) , H(3,:) )
  G(1,2) = DOT_PRODUCT( H(1,:) , H(2,:) )
  G(2,1) = G(1,2)
  G(1,3) = DOT_PRODUCT( H(1,:) , H(3,:) )
  G(3,1) = G(1,3)
  G(2,3) = DOT_PRODUCT( H(2,:) , H(3,:) )
  G(3,2) = G(2,3)
  !Invert metric tenro
  CALL INVMAT(G,Gi)
ENDIF
!
!Go back to beginning of file and store data
REWIND(30)
i = 0
DO
  i = i+1 !Line counter
  READ(30,'(a128)',ERR=200,END=200) temp
  temp = ADJUSTL(temp)
  IF( temp(1:5)=="atoms" .OR. temp(1:5)=="ATOMS" ) THEN
    DO n=1,NP
      READ(30,'(a128)',ERR=207,END=207) temp
      temp = ADJUSTL(temp)
      !
      !Remove commas from the line
      strlength = SCAN(temp,',')
      DO WHILE(strlength>0)
        temp(strlength:strlength) = " "
        strlength = SCAN(temp,',')
      ENDDO
      !
      !Convert species to proper casing
      species(1:1) = StrUpCase(temp(1:1))
      species(2:2) = StrDnCase(temp(2:2))
      !Save atomic number
      CALL ATOMNUMBER(species,P(n,4))
      !
      !Read atom position from the rest of the line
      READ(temp(4:),*,ERR=207,END=207) P(n,1), P(n,2), P(n,3)
      !
      IF( pdffit ) THEN
        !Additional info about current atom in the next 5 lines
        !Read occupancy at the end of current line
        READ(temp(4:),*,ERR=207,END=207) a, b, c, AUX(n,1)
        !Read standard deviations (not saved)
        READ(30,'(a128)',ERR=207,END=207) temp
        !Read anisotropic temperature factors U11, U22, U33
        READ(30,*,ERR=207,END=207) AUX(n,2), AUX(n,3), AUX(n,4)
        !Read deviations (not saved)
        READ(30,'(a128)',ERR=207,END=207) temp
        !Read anisotropic temperature factors U12, U13, U23
        READ(30,*,ERR=207,END=207) AUX(n,5), AUX(n,6), AUX(n,7)
        !Read deviations (not saved)
        READ(30,'(a128)',ERR=207,END=207) temp
        !Set "isotropic" temperature factor Biso = 8pi²(U11+U22+U33)/3
        AUX(n,8) = 8.d0*pi*pi*( AUX(n,2) + AUX(n,3) + AUX(n,4) ) / 3.d0
      ELSE
        !DISCUS format
        !Set occupancy of this atom to 1
        AUX(n,1) = 1.d0
        !Read Biso at the end of current line
        READ(temp(4:),*,ERR=207,END=207) a, b, c, AUX(n,8)
        !Convert Biso into "anisotropic" factors: U11 = U22 = U33 = B / (8pi²)
        AUX(n,2) = AUX(n,8) / (8.d0*pi*pi)   ! U11
        AUX(n,3) = AUX(n,8) / (8.d0*pi*pi)   ! U22
        AUX(n,4) = AUX(n,8) / (8.d0*pi*pi)   ! U33
        AUX(n,5) = 0.d0   ! U12
        AUX(n,6) = 0.d0   ! U13
        AUX(n,7) = 0.d0   ! U23
      ENDIF
      !
    ENDDO
    GOTO 200
  ENDIF
ENDDO
!
!
200 CONTINUE
CLOSE(30)
!Convert reduced coordinates to cartesian
CALL FRAC2CART(P,H)
!
!If a scaling factor (different from 1) was provided, apply it
IF( DABS(a0-1.d0)>1.d-12 .AND. a0>1.d-6 ) THEN
  !Rescale the cell
  H(:,:) = a0*H(:,:)
  !Rescale all atom positions
  P(:,1:3) = a0*P(:,1:3)
ENDIF
!
!Apply symmetry operations (if any)
!(cf. /include/symops.f90)
IF( LEN_TRIM(spcgr)>0 ) THEN
  CALL SG_APPLY_SYMOPS(spcgr,H,P,S,AUXNAMES,AUX)
ENDIF
!
!If system must be duplicated, do it
IF( Nx*Ny*Nz>1 ) THEN
  NP = SIZE(P,1)*Nx*Ny*Nz
  ALLOCATE(Q(NP,4))
  NP=0
  DO o=1,Nz
    DO n=1,Ny
      DO m=1,Nx
        DO i=1,SIZE(P,1)
          NP = NP+1
          Q(NP,1) = P(i,1)+REAL(m-1)*H(1,1)+REAL(n-1)*H(2,1)+REAL(o-1)*H(3,1)
          Q(NP,2) = P(i,2)+REAL(m-1)*H(1,2)+REAL(n-1)*H(2,2)+REAL(o-1)*H(3,2)
          Q(NP,3) = P(i,3)+REAL(m-1)*H(1,3)+REAL(n-1)*H(2,3)+REAL(o-1)*H(3,3)
          Q(NP,4) = P(i,4)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !Replace old P by new Q
  DEALLOCATE(P)
  ALLOCATE(P(SIZE(Q,1),4))
  P(:,:) = Q(:,:)
  DEALLOCATE(Q)
  !
  !Also duplicate auxiliary properties
  ALLOCATE( Q( SIZE(P,1) , SIZE(AUXNAMES) ) )
  NP = 0
  DO o=1,Nz
    DO n=1,Ny
      DO m=1,Nx
        DO i=1,SIZE(AUX,1)
          NP = NP+1
          Q(NP,:) = AUX(i,:)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  DEALLOCATE(AUX)
  ALLOCATE(AUX(SIZE(Q,1),SIZE(Q,2)))
  AUX(:,:) = Q(:,:)
  DEALLOCATE(Q)
  !
  !Resize the cell
  H(1,:) = Nx*H(1,:)
  H(2,:) = Ny*H(2,:)
  H(3,:) = Nz*H(3,:)
ENDIF
GOTO 1000
!
207 CONTINUE
CLOSE(30)
CALL ATOMSK_MSG(807,(/''/),(/DBLE(i)/))
nerr = nerr+1
GOTO 1000
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
END SUBROUTINE READ_STR
!
END MODULE in_str
