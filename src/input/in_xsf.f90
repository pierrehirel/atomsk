MODULE in_xsf
!
!
!**********************************************************************************
!*  IN_XSF                                                                        *
!**********************************************************************************
!* This module reads XSF format, initially designed for XCrysDen.                 *
!* The XSF format is described here:                                              *
!*    http://www.xcrysden.org/doc/XSF.html                                        *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 19 March 2014                                    *
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
USE atoms
USE comv
USE constants
USE messages
USE files
USE subroutines
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE READ_XSF(inputfile,H,P,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=7):: xsfformat
CHARACTER(LEN=128):: temp, msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(100):: tempcomment
LOGICAL:: atoms   !true if molecular structure, false if 2-D or 3-D periodic
LOGICAL:: forces  !are forces present in the file?
INTEGER:: i, j, k, l, NP
REAL(dp):: a0, tempreal
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
xsfformat = ''
atoms = .FALSE.
forces = .FALSE.
NP = 0
a0 = 1.d0
!
100 CONTINUE
msg = 'entering READ_XSF'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
i=-1
l=0
DO WHILE(NP==0)
  READ(30,'(a128)',ERR=400,END=400) temp
  temp = TRIM(ADJUSTL(temp))
  IF(temp(1:1)=='#') THEN
    l=l+1
    IF(l<=100) THEN
      tempcomment(l) = TRIM(temp)
    ENDIF
  ELSEIF(temp(1:7)=='PRIMVEC') THEN
    DO j=1,3
      READ(30,*,ERR=410,END=410) (H(j,k), k=1,3)
    ENDDO
  ELSEIF(temp(1:9)=='PRIMCOORD') THEN
    READ(30,*,ERR=420,END=420) NP
  ELSEIF(temp(1:5)=='ATOMS' ) THEN
    atoms = .TRUE.
    !count atoms
    DO
      READ(30,*,ERR=150,END=150) i
      NP = NP+1
    ENDDO
    150 CONTINUE
    !Go back to the first atom
    DO i=1,NP
      BACKSPACE(30)
    ENDDO
  ENDIF
ENDDO
!
200 CONTINUE
WRITE(msg,*) 'NP:', NP
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
ALLOCATE(P(NP,4))
P(:,:) = 0.d0
!
!If comments were found, save them
IF(l>0) THEN
  ALLOCATE(comment(l))
  DO i=1,l
    comment(i) = tempcomment(i)
  ENDDO
ENDIF
!
!Determine if the format is
! (number x y z)  or  (species x y z)
READ(30,*,ERR=430,END=430) temp
temp=TRIM(ADJUSTL(temp))
!If temp(1:2) contains a species name, converting it will give tempreal.NE.0
CALL ATOMNUMBER(temp(1:2),tempreal)
!
IF(tempreal.NE.0.d0) THEN
  !A species name was read: it is "speciesxyz" format
  xsfformat = 'spiexyz'
!
ELSE
  !Otherwise the file should have "massxyz" format
  xsfformat = 'massxyz'
  !Let's confirm it by reading the real number
  !If there is an error then we are screwed -> exit
  READ(temp,*,ERR=430,END=430) tempreal
  !Then convert this real to atom species
  CALL ATOMSPECIES(tempreal,species)
  !And convert this species back to a real
  CALL ATOMNUMBER(species,tempreal)
  !Now the real should be the atomic number
  !If it is still zero then we are in trouble -> exit
  IF(tempreal==0.d0) GOTO 430
ENDIF
!
!
250 CONTINUE
BACKSPACE(30)
i=1
!Read first atom position and count number of columns
READ(30,'(a128)',ERR=400,END=400) temp
IF(xsfformat=='massxyz') THEN
  READ(temp,*,ERR=400,END=400) P(1,4), P(1,1), P(1,2), P(1,3)
ELSE
  READ(temp,*,ERR=400,END=400) species, P(1,1), P(1,2), P(1,3)
  CALL ATOMNUMBER(species,P(1,4))
ENDIF
!Check if forces are present after atom position
READ(temp,*,ERR=260,END=260) msg, msg, msg, msg, msg, msg, msg
forces = .TRUE.
ALLOCATE(AUXNAMES(3))
AUXNAMES(1) = "fx"
AUXNAMES(2) = "fy"
AUXNAMES(3) = "fz"
ALLOCATE(AUX(SIZE(P,1),3))
AUX(:,:) = 0.d0
!
!
260 CONTINUE
!Read all atomic positions
BACKSPACE(30)
DO i=1,NP
  IF(xsfformat=='massxyz') THEN
    IF(forces) THEN
      READ(30,*,ERR=400,END=400) P(i,4), P(i,1), P(i,2), P(i,3), AUX(i,1), AUX(i,2), AUX(i,3)
    ELSE
      READ(30,*,ERR=400,END=400) P(i,4), P(i,1), P(i,2), P(i,3)
    ENDIF
    !
  ELSE
    IF(forces) THEN
      READ(30,*,ERR=400,END=400) species, P(i,1), P(i,2), P(i,3), AUX(i,1), AUX(i,2), AUX(i,3)
    ELSE
      READ(30,*,ERR=400,END=400) species, P(i,1), P(i,2), P(i,3)
    ENDIF
    CALL ATOMNUMBER(species,P(i,4))
  ENDIF
ENDDO
GOTO 500
!
!
!
400 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
GOTO 500
!
410 CONTINUE
CALL ATOMSK_MSG(1809,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 500
!
420 CONTINUE
CALL ATOMSK_MSG(1810,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 500
!
430 CONTINUE
CALL ATOMSK_MSG(1804,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 500
!
!
!
500 CONTINUE
CLOSE(30)
!
END SUBROUTINE READ_XSF
!
END MODULE in_xsf
