MODULE in_abinit
!
!
!**********************************************************************************
!*  IN_ABINIT                                                                     *
!**********************************************************************************
!* This module reads a file in the ABINIT format.                                 *
!* The ABINIT documentation is available here:                                    *
!*     https://www.abinit.org/                                                    *
!**********************************************************************************
!* (C) May 2019 - Pierre Hirel                                                    *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 16 April 2024                                    *
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
USE symops
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE READ_ABINIT(inputfile,H,P,S,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=128):: line
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
INTEGER:: i, j, k
INTEGER:: nline  !counter for the line number
INTEGER:: natom  !counter for total number of particles
INTEGER:: ntypat !number of different atom types
INTEGER:: sgnumber  !space group number
INTEGER,DIMENSION(:,:),ALLOCATABLE:: atypes  !atomic numbers and types
REAL(dp):: a, b, c
REAL(dp):: scalecart
REAL(dp),DIMENSION(3):: acell  !factor for cell vectors
REAL(dp),DIMENSION(3):: angdeg !angles of box vectors (degrees)
REAL(dp),DIMENSION(3,3),INTENT(OUT):: H    !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: P    !Atom positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S    !Shell positions (dummy array, not used here)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: AUX  !auxiliary properties
!
!
!Initialize variables
nline = 0
natom = 0
ntypat = 0
sgnumber = 0
scalecart = 1.d0
acell(:) = 1.d0
H(:,:) = 0.d0
DO i=1,3
  H(i,i) = 1.d0
ENDDO
!
msg = 'entering READ_ABINIT'
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
  nline = nline+1
  READ(30,'(a128)',ERR=200,END=200) line
  line = ADJUSTL(line)
  !
  IF( line(1:1)=="#" ) THEN
    !This line is a comment
    comment(1) = "# "//TRIM(line)
    !
  ELSEIF( line(1:5)=='acell' ) THEN
    !Read cell lattice vector scaling
    line = ADJUSTL(line(6:))
    !
    j = INDEX(line,"*")
    IF( j>0 ) THEN
      !This is of the form "3*acell"
      READ(line(j+1:),*,END=800,ERR=800) acell(1)
      acell(2) = acell(1)
      acell(3) = acell(1)
    ELSE
      !Read the three values of acell(:)
      READ(line,*,END=800,ERR=800) acell(1), acell(2), acell(3)
    ENDIF
    !
  ELSEIF( line(1:9)=='scalecart' ) THEN
    !Read cell lattice vector scaling
    line = ADJUSTL(line(10:))
    !
    j = INDEX(line,"*")
    IF( j>0 ) THEN
      !This is of the form "3*scalecart"
      READ(line(j+1:),*,END=800,ERR=800) scalecart
    ELSE
      READ(line,*,END=800,ERR=800) scalecart
    ENDIF
    !
  ELSEIF( line(1:6)=='angdeg' ) THEN
    !Read angles between cell vectors
    READ(line,*,END=800,ERR=800) angdeg(1), angdeg(2), angdeg(3)
    a = VECLENGTH(H(1,:))
    b = VECLENGTH(H(2,:))
    c = VECLENGTH(H(3,:))
    H(:,:) = 0.d0
    IF( angdeg(1)-angdeg(2)<1.d-12 .AND. angdeg(2)-angdeg(3)<1.d-12 .AND. DABS(angdeg(1)-90.d0)>1.d-12 ) THEN
      H(1,1) = a
      H(1,3) = c
       H(2,1) = -1.d0*a/2.d0
       H(2,2) = a*DSQRT(3.d0)/2.d0
       H(2,3) = c
      H(3,1) = -1.d0*a/2.d0
      H(3,2) = -1.d0*a*DSQRT(3.d0)/2.d0
      H(3,3) = c
    ELSE
      H(1,1) = 1.d0
       H(2,1) = a
       H(2,2) = b
      H(3,1) = c
      H(3,2) = c
      H(3,3) = c
    ENDIF
    !
  ELSEIF( line(1:5)=='rprim' ) THEN
    !Read cell lattice vector scaling
    READ(line(6:),*,END=800,ERR=800) H(1,1), H(1,2), H(1,3)
    nline = nline+1
    READ(30,*,END=800,ERR=800) H(2,1), H(2,2), H(2,3)
    nline = nline+1
    READ(30,*,END=800,ERR=800) H(3,1), H(3,2), H(3,3)
    !Apply multiplication factors
    H(1,:) = scalecart * acell(1) * H(1,:)
    H(2,:) = scalecart * acell(2) * H(2,:)
    H(3,:) = scalecart * acell(3)*H(3,:)
    !
  ELSEIF( line(1:5)=='natom' ) THEN
    !Read number of atoms
    IF( .NOT. ALLOCATED(P) ) THEN
      READ(line(6:),*,END=800,ERR=800) natom
      !Allocate array P to store atom positions
      IF( natom>0 ) THEN
        ALLOCATE(P(natom,4))
        P(:,:) = 0.d0
        !Allocate array AUX to save atom types
        ALLOCATE(AUX(natom,1))
        AUX(:,:) = 0.d0
        ALLOCATE(AUXNAMES(1))
        AUXNAMES(1) = "type"
      ELSE
        nerr = nerr+1
        GOTO 1000
      ENDIF
    ELSE
      !ERROR: cannot change number of atoms
      nerr = nerr+1
      GOTO 1000
    ENDIF
    !
  ELSEIF( line(1:6)=='ntypat' ) THEN
    !Read number of different types of atoms
    READ(line(7:),*,END=800,ERR=800) ntypat
    IF( .NOT.ALLOCATED(atypes) .AND. ntypat>0 ) THEN
      ALLOCATE(atypes(ntypat,2))
      atypes(:,:) = 0
    ELSE
    
    ENDIF
    !
  ELSEIF( line(1:5)=='typat' ) THEN
    !Read number of atoms of each type
    READ(line(6:),*,END=800,ERR=800) (AUX(j,1),j=1,natom)
    !
  ELSEIF( line(1:5)=='znucl' ) THEN
    !Read atomic numbers
    READ(line(6:),*,END=800,ERR=800) (atypes(j,1),j=1,ntypat)
    !
  ELSEIF( line(1:5)=='xcart' .OR. line(1:6)=='xangst' .OR. line(1:4)=='xred' ) THEN
    IF( line(1:4)=='xred' ) THEN
      k = 5
    ELSEIF( line(1:4)=='xangst' ) THEN
      k = 7
    ELSE
      k = 6
    ENDIF
    !
    !If array P was not allocated before, we are in trouble
    IF( ALLOCATED(P) .AND. SIZE(P,1)>0 ) THEN
      !Read atom positions
      READ(line(k:),*) P(1,1), P(1,2), P(1,3)
      P(1,4) = atypes(1,1)
      DO i=2,SIZE(P,1)
        nline = nline+1
        READ(30,*,END=800,ERR=800) P(i,1), P(i,2), P(i,3)
        P(i,4) = DBLE(atypes(NINT(AUX(i,1)),1))
      ENDDO
      !
      IF( k==5) THEN
        !Convert coordinates to cartesian
        CALL FRAC2CART(P,H)
      ENDIF
      !
    ELSE
      nerr = nerr+1
      GOTO 1000
    ENDIF
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
!If all cell vectors are set to 1, it means they are not set => reset them to zero
IF( DABS(VECLENGTH(H(1,:)))-1.d0<1.d-12 .AND.         &
  & DABS(VECLENGTH(H(2,:)))-1.d0<1.d-12 .AND.         &
  & DABS(VECLENGTH(H(3,:)))-1.d0<1.d-12      ) THEN
  H(:,:) = 0.d0
ENDIF
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(807,(/''/),(/DBLE(nline)/))
nerr = nerr+1
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE READ_ABINIT
!
END MODULE in_abinit
