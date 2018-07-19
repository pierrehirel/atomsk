MODULE in_jems
!
!
!**********************************************************************************
!*  IN_JEMS                                                                       *
!**********************************************************************************
!* This module reads a crystal structure for JEMS.                                *
!* The JEMS format is used by the TEM simulation software JEMS, and is described  *
!* on the software website:                                                       *
!*     http://www.jems-saas.ch/                                                   *
!**********************************************************************************
!* (C) October 2015 - Pierre Hirel                                                *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 01 March 2017                                    *
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
SUBROUTINE READ_JEMS(inputfile,H,P,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: msg, temp
CHARACTER(LEN=128):: sgroup
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: isreduced
INTEGER:: j
INTEGER:: NP   ! number of atoms
INTEGER:: NSym ! number of symmetry operations
REAL(dp):: a, b, c, alpha, beta, gamma !supercell (conventional notation)
REAL(dp):: P1, P2, P3, dw, occ, absorption
REAL(dp),DIMENSION(3,3),INTENT(OUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: P   !Positions of atoms
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S   !Positions of shells (not used here)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: AUX !auxiliary properties
!
!
!Initialize variables
NP=0
sgroup=""
a=0.d0
b=0.d0
c=0.d0
alpha=0.d0
beta=0.d0
gamma=0.d0
H(:,:) = 0.d0
!
IF (ALLOCATED(P)) DEALLOCATE(P)
IF (ALLOCATED(AUX)) DEALLOCATE(AUX)
IF (ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
!
msg = 'entering READ_JEMS'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
!
100 CONTINUE
OPEN(UNIT=30,FILE=inputfile,STATUS='UNKNOWN',ERR=800)
REWIND(30)
!
!Parse the file a first time to count atoms and number of symmetry operations
NP=0
Nsym=0
DO
  READ(30,'(a128)',ERR=200,END=200) temp
  temp = ADJUSTL(temp)
  IF( temp(1:5)=="atom|" ) THEN
    NP = NP+1
  !
  ELSEIF( temp(1:3)=="rps" ) THEN
    Nsym=Nsym+1
  ENDIF
ENDDO
!
!
!
200 CONTINUE
WRITE(msg,*) 'Counted atoms, NP = ', NP
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
IF (NP<=0) GOTO 800 ! No atoms => error
!
IF( Nsym>0 ) THEN
  IF(ALLOCATED(symops_trf)) DEALLOCATE(symops_trf)
  ALLOCATE(symops_trf(12,Nsym))
  symops_trf(:,:) = 0.d0
ENDIF
!
! Allocate the atom data P and auxiliary data AUX and AUXNAMES
ALLOCATE(P(NP,4))
P(:,:) = 0.d0
ALLOCATE(AUX(NP,3))
AUX(:,:) = 0.d0
ALLOCATE(AUXNAMES(3))
AUXNAMES(1) = 'Debye-Waller'
AUXNAMES(2) = 'occ'
AUXNAMES(3) = 'absorption'
!
!Go back to beginning of file and store data
REWIND(30)
Nsym=0
DO
  READ(30,'(a128)',ERR=300,END=300) temp
  temp = ADJUSTL(temp)
  !Replace all commas by blank spaces
  j=SCAN(temp,",")
  DO WHILE(j>0)
    temp(j:j) = " "
    j=SCAN(temp,",")
  ENDDO
  !Replace vertical bars by blank spaces
  j=SCAN(temp,"|")
  DO WHILE(j>0)
    temp(j:j) = " "
    j=SCAN(temp,"|")
  ENDDO
  !
  IF( temp(1:8)=="HMSymbol" ) THEN
    !This line contains the space group number
    !Replace vertical bars by blank spaces
    j=SCAN(temp,"|")
    DO WHILE(j>0)
      temp(j:j) = " "
      j=SCAN(temp,"|")
    ENDDO
    !Read the space group number
    READ(temp,*,ERR=300,END=300) msg, sgroup
  ENDIF
  !
  IF( temp(1:3)=="rps" ) THEN
    !This line contains a symmetry operation
    Nsym=Nsym+1
    !Replace vertical bars by blank spaces
    j=SCAN(temp,"|")
    DO WHILE(j>0)
      temp(j:j) = " "
      j=SCAN(temp,"|")
    ENDDO
    !Read the space group number
    READ(temp,*,ERR=300,END=300) msg
  ENDIF
  !
  IF( temp(1:8)=="lattice " ) THEN
    READ(temp(9:),*,ERR=800,END=800) j, P1
    IF(j==0) THEN
      a = P1
    ELSEIF(j==1) THEN
      b = P1
    ELSEIF(j==2) THEN
      c = P1
    ELSEIF(j==3) THEN
      alpha = DEG2RAD(P1)
    ELSEIF(j==4) THEN
      beta = DEG2RAD(P1)
    ELSEIF(j==5) THEN
      gamma = DEG2RAD(P1)
    ENDIF
  ELSEIF( temp(1:5)=="atom " ) THEN
    !Read atom position
    !NOTE: in JEMS files the first atom has an index zero, while Atomsk starts counting at 1
    !NOTE2: here we attempt to read atom position, and return an error if it fails (GOTO 800).
    !      However, even if Debye-Waller, occupancies and absorption coefficients are supposed
    !      to be present in a proper JEMS file, Atomsk will not fail if they are not present,
    !      it will just ignore them and go on.
    READ(temp(6:),*,ERR=800,END=800) j, species, msg, P1, P2, P3
    READ(temp(6:),*,ERR=220,END=220) j, species, msg, P1, P2, P3, dw, occ, absorption
    220 CONTINUE
    IF( j+1>0 .AND. j+1<=SIZE(P,1) ) THEN
      P(j+1,1) = P1
      P(j+1,2) = P2
      P(j+1,3) = P3
      CALL ATOMNUMBER(species,P(j+1,4))
      AUX(j+1,1) = dw
      AUX(j+1,2) = occ
      AUX(j+1,3) = absorption
    ELSE
      !Out-of-range atom
      nwarn=nwarn+1
      CALL ATOMSK_MSG(2742,(/''/),(/DBLE(j+1)/))
    ENDIF
  ENDIF
ENDDO
!
!
300 CONTINUE
CLOSE(30)
!Save cell vectors in H(:,:)
CALL CONVMAT(a,b,c,alpha,beta,gamma,H)
!Find out if coordinates are reduced or cartesian
CALL FIND_IF_REDUCED(P,isreduced)
!In case of reduced coordinates, convert them to cartesian
IF(isreduced) THEN
  CALL FRAC2CART(P,H)
ENDIF
!Apply symmetry operations (if any)
!(cf. /include/symops.f90)
IF( LEN_TRIM(sgroup)>0 ) THEN
  CALL SG_APPLY_SYMOPS(sgroup,H,P,S,AUXNAMES,AUX)
ENDIF
GOTO 1000
!
!
800 CONTINUE
801 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(j+1)/))
nerr = nerr+1
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE READ_JEMS
!
END MODULE in_jems
