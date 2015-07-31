MODULE in_mbpp_coorat
!
!
!**********************************************************************************
!*  IN_MBPP_COORAT                                                                *
!**********************************************************************************
!* This module reads files in the COORAT format, which is used by the             *
!* Mixed-Basis PseudoPotential (MBPP) code by B. Meyer, C. Elsässer,              *
!* F. Lechermann, and M. Fähnle, Max-Plank-Institut für Metalforschung,           *
!* Stuttgart (unpublished).                                                       *
!**********************************************************************************
!* (C) Feb. 2011 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
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
USE comv
USE constants
USE messages
USE files
USE subroutines
!
!
CONTAINS
!
SUBROUTINE READ_COORAT(inputfile,H,P,comment,AUXNAMES,AUX)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=1):: INPchar  !special characters used in INP files
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: msg, Hleft, Hright, test
CHARACTER(LEN=32),DIMENSION(3):: Hstring
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: fileexists
INTEGER:: i, j, strpos
INTEGER:: NP, NPsp
REAL(dp):: a0, atn
REAL(dp):: Vcell, VH !volume of the cell, volume of H
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
i = 0
NP = 0
NPsp = 0
a0 = 1.d0
VH = 0.d0
H(:,:) = 0.d0
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
!
!
100 CONTINUE
msg = 'entering READ_COORAT'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
!Determine the total number of particles
DO
  READ(30,'(a128)',ERR=150,END=150) test
  test = ADJUSTL(test)
  DO i=1,124
    IF(test(i:i+4)=='natom' .OR. test(i:i+4)=='NATOM') THEN
      !Read number of atoms for this species
      strpos = i+SCAN(test(i:128),'=')
      READ(test(strpos:),*,END=420,ERR=420) NPsp
      NP = NP+NPsp
    ENDIF
  ENDDO
ENDDO
!
150 CONTINUE
WRITE(msg,*) 'NP = ', NP
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!Allocate array P
IF(NP.NE.0) THEN
  ALLOCATE(P(NP,4))
  P(:,:) = 0.d0
ELSE
  GOTO 400
ENDIF
!
!Read atom positions and species
REWIND(30)
NP = 0
DO WHILE( NP<SIZE(P(:,1)) )
  READ(30,'(a128)',ERR=400,END=400) test
  test = ADJUSTL(test)
  NPsp=0
  atn=0.d0
  i=0
  DO WHILE(NPsp==0 .OR. atn==0.d0 .OR. i<125)
    i=i+1
    IF(test(i:i+4)=='natom' .OR. test(i:i+4)=='NATOM') THEN
      !Read number of atoms for this species
      strpos = i+SCAN(test(i:),'=')
      READ(test(strpos:),*) NPsp
    ELSEIF(test(i:i+3)=='name' .OR. test(i:i+3)=='NAME') THEN
      !Now read the atom species
      strpos = i+SCAN(test(i:),'=')
      READ(test(strpos:),*) species
      CALL ATOMNUMBER(species,atn)
    ENDIF
  ENDDO
  !
  !Read atom positions for this species
  DO i=1,NPsp
    NP = NP+1
    READ(30,*,ERR=400,END=400) P(NP,1), P(NP,2), P(NP,3)
    P(NP,4) = atn
  ENDDO
ENDDO
CLOSE(30)
!
!
!
200 CONTINUE
!Search for file named "INP" and read supercell
INQUIRE(FILE='INP',EXIST=fileexists)
IF(fileexists) THEN
  CALL ATOMSK_MSG(1002,(/''/),(/0.d0/))
  !
  OPEN(UNIT=31,FILE='INP',FORM='FORMATTED',STATUS='OLD')
  REWIND(31)
  READ(31,'(a64)',ERR=500,END=500) comment
  DO
    READ(31,'(a128)',ERR=500,END=500) test
    test = ADJUSTL(test)
    IF(test(1:4)=='alat' .OR. test(1:3)=='vol') THEN
      IF(test(1:4)=='alat') THEN
        !Read scaling factor
        strpos = SCAN(test,'=')
        READ(test(strpos+1:),*) a0
      ELSE
        !Read cell volume
        strpos = SCAN(test,'=')
        READ(test(strpos+1:),*) Vcell
      ENDIF
      !
      !Read supercell parameters H(i,j)
      DO i=1,3
        READ(31,*,ERR=400,END=400) Hstring(1), Hstring(2), Hstring(3)
        DO j=1,3
          Hstring(j) = TRIM(ADJUSTL( Hstring(j) ))
          !Check if it is just a character
          IF( SCAN(Hstring(j),'ehlnorstxy')>0 ) THEN
            strpos = SCAN(Hstring(j),'*')
            IF(strpos==0) THEN
              !If there is no star then it is a number and we are done
              READ(Hstring(j),*) H(i,j)
              GOTO 210
            !
            ELSE
              !If there is a '*' we must read the numbers on both sides
              READ(Hstring(j)(1:strpos-1),*) Hleft
              READ(Hstring(j)(strpos+1:),*)  Hright
              IF( SCAN(Hstring(j),'ehlnorstxy')>strpos ) THEN
                READ(Hleft,*) H(i,j)
                READ(Hright,*) INPchar
              ELSE
                READ(Hright,*) H(i,j)
                READ(Hleft,*) INPchar
              ENDIF
              H(i,j) = H(i,j)*INP2DBLE(INPchar)
            ENDIF
            !
          ELSE
            !no letter in Hstring => it should be a number
            READ(Hstring(j),*,ERR=400,END=400) H(i,j)
          ENDIF
          210 CONTINUE
        ENDDO  !loop on j
      ENDDO  !loop on i
      !
      !Scale H
      IF(test(1:4)=='alat') THEN
        !just scale all vectors to alat
        H(:,:) = a0*H(:,:)
      ELSE
        !Compute current volume of H
        CALL VOLUME_PARA(H,VH)
        !scale cell vectors
        H(:,:) = ( (Vcell/VH)**(1.d0/3.d0) ) * H(:,:)
      ENDIF
      !
      !Convert to cartesian coordinates
      CALL FRAC2CART(P,H)
    ENDIF
    !
  ENDDO
  !
  CLOSE(31)
  !
ENDIF
GOTO 500
!
!
!
400 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(NP)/))
nerr = nerr+1
GOTO 1000
!
410 CONTINUE
CALL ATOMSK_MSG(1809,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
420 CONTINUE
CALL ATOMSK_MSG(1810,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
430 CONTINUE
CALL ATOMSK_MSG(1804,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
500 CONTINUE
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE READ_COORAT
!
!
!
!********************************************************
! INP2DBLE
! This function transforms some characters used in
! the INP file into a real number.
!********************************************************
FUNCTION INP2DBLE(INPchar) RESULT(value)
CHARACTER(LEN=1):: INPchar
REAL(dp):: value
!
value=0.d0
SELECT CASE(INPchar)
!
CASE('e')
  value = DSQRT(10.d0)
CASE('h')
  value = DSQRT(3.d0)
CASE('l')
  value = DSQRT(5.d0)
CASE('n')
  value = 9.d0
CASE('o')
  value = 1.d0
CASE('r')
  value = DSQRT(2.d0)
CASE('s')
  value = 7.d0
CASE('t')
  value = 3.d0
CASE('x')
  value= DSQRT(8.d0/3.d0)
CASE('y')
  value = DSQRT(6.d0)
CASE DEFAULT
  !impossible to recognize value => error
  nerr=nerr+1
END SELECT
!
END FUNCTION INP2DBLE
!
!
!
END MODULE in_mbpp_coorat
