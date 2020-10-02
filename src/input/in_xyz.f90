MODULE in_xyz
!
!**********************************************************************************
!*  IN_XYZ                                                                        *
!**********************************************************************************
!* This module reads files in the XYZ format.                                     *
!* No standard specification exists, but a widely used format is                  *
!* described for instance at:                                                     *
!*     http://openbabel.org/wiki/XYZ_%28format%29                                 *
!* This module can also read extended XYZ format, where the comment line          *
!* (i.e. the 2nd line) is replaced by a set of keywords/values. The               *
!* extended XYZ format is described for instance at:                              *
!*     http://www.jrkermode.co.uk/quippy/extended_xyz.html                        *
!* Finally, this module can also read special XYZ format, which contains          *
!* additional information at the end of the file:                                 *
!*          NP                                                                    *
!*          #comment line                                                         *
!*          species1 x1 y1 z1                                                     *
!*          species2 x2 y2 z2                                                     *
!*           ...                                                                  *
!*          alat                                                                  *
!*          <a0>                                                                  *
!*          supercell                                                             *
!*          <H(1,1)> <H(1,2)> <H(1,3)>                                            *
!*          <H(2,1)> <H(2,2)> <H(2,3)>                                            *
!*          <H(3,1)> <H(3,2)> <H(3,3)>                                            *
!*          <cartesian/fractional> coordinates                                    *
!*          mass <species1> <smass1>                                              *
!*          mass <species2> <smass2>                                              *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 02 Oct. 2020                                     *
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
!
USE atoms
USE comv
USE constants
USE messages
USE files
USE subroutines
!
!
CONTAINS
!
SUBROUTINE READ_XYZ(inputfile,H,P,comment,AUXNAMES,AUX)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=1):: spchar
CHARACTER(LEN=2):: species
CHARACTER(LEN=7):: xyzformat
CHARACTER(LEN=128):: msg
CHARACTER(LEN=1024):: properties, temp, temp2
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: isreduced, Hset
LOGICAL:: xyz, exyz, sxyz !is it basic XYZ, Extended XYZ, Special XYZ? Cannot be all 3
INTEGER:: auxiliary  !number of auxiliary properties
INTEGER:: i, j, k, l, NP, strlength
REAL(dp):: a0, tempreal
REAL(dp):: a, b, c, alpha, beta, gamma
REAL(dp):: P1, P2, P3
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P  !atom positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
isreduced = .FALSE.
Hset=.FALSE.
xyz=.FALSE.
exyz=.FALSE.
sxyz=.FALSE.
auxiliary = 0
i = 0
j = 0
a0 = 1.d0
xyzformat = 'massxyz'
!
!
msg = 'entering READ_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!
100 CONTINUE
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!First line is always the number of atoms
READ(30,*,ERR=820,END=820) a
IF( a > NATOMS_MAX ) THEN
  nerr = nerr+1
  CALL ATOMSK_MSG(821,(/""/),(/a/))
  GOTO 1000
ENDIF
NP = NINT(a)
IF(NP==0) GOTO 800
ALLOCATE(P(NP,4) , STAT=i)
IF( i>0 ) THEN
  ! Allocation failed (not enough memory)
  nerr = nerr+1
  CALL ATOMSK_MSG(819,(/''/),(/DBLE(NP)/))
  GOTO 1000
ENDIF
!
!In case of special XYZ format, we have to read info at end of file
!This is only possible if the XYZ file contains a single system,
!i.e. if the mode is different from "1-in-all".
msg = 'determine if special XYZ or not'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
DO
  READ(30,*,ERR=145,END=145) temp
  temp = TRIM(ADJUSTL(temp))
  IF(temp=='alat') THEN
    sxyz=.TRUE.
    READ(30,*,ERR=800,END=800) a0
  ELSEIF(temp(1:5)=='super') THEN
    DO i=1,3
      READ(30,*,ERR=810,END=810) (H(i,j), j=1,3)
    ENDDO
    Hset = .TRUE.
  ELSEIF(temp(1:7)=='convent') THEN
    READ(30,*,ERR=810,END=810) a, b, c
    READ(30,*,ERR=810,END=810) alpha, beta, gamma
    alpha = DEG2RAD(alpha)
    beta = DEG2RAD(beta)
    gamma = DEG2RAD(gamma)
    CALL CONVMAT(a,b,c,alpha,beta,gamma,H)
    Hset = .TRUE.
  ELSEIF(temp(1:8)=='property') THEN
    auxiliary = auxiliary+1
  ENDIF
ENDDO
!
145 CONTINUE
IF(Hset) THEN
  H(:,:) = a0*H(:,:)
  msg = 'special XYZ detected'
  sxyz=.TRUE.
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!If auxiliary properties are present, read their name
IF(auxiliary>0) THEN
  WRITE(msg,*) 'Number of auxiliary properties:', auxiliary
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  ALLOCATE( AUXNAMES(auxiliary) )
  AUXNAMES(:) = ''
  ALLOCATE( AUX(NP,auxiliary) )
  AUX(:,:) = 0.d0
  REWIND(30)
  DO
    READ(30,'(a1024)',ERR=200,END=200) temp
    temp = TRIM(ADJUSTL(temp))
    IF(temp(1:8)=='property') THEN
      READ(temp(9:),*,ERR=830,END=830) auxiliary
      READ(temp(9:),*,ERR=830,END=830) i, AUXNAMES(auxiliary)
    ENDIF
  ENDDO
ENDIF
!
!
!
200 CONTINUE
REWIND(30)
!Read the file
READ(30,*,ERR=820,END=820) NP
READ(30,'(a1024)',ERR=830,END=830) properties !Properties if extended XYZ, comment otherwise
!
!Parse this second line for keywords
!If it has any then it is extended XYZ => read properties
DO i=1,1020
  !Parse string "properties"
  temp = properties(i:1024)
  IF(LEN_TRIM(temp)==0) EXIT
  !
  IF(temp(1:9)=='Lattice="') THEN
    msg = 'EXYZ: found Lattice'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    exyz=.TRUE.
    strlength = SCAN(temp(10:),'"')+8
    READ(temp(10:strlength),*,END=250,ERR=250) H(1,1), H(2,1), H(3,1), &
        & H(1,2), H(2,2), H(3,2), H(1,3), H(2,3), H(3,3)
    Hset=.TRUE.
    !
  ELSEIF(temp(1:11)=='Properties=') THEN
    msg = 'EXYZ: found Properties'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    exyz=.TRUE.
    temp = temp(12:)
    !Read auxiliary properties; they should have the format:
    !  Properties=species:S:1:pos:R:3:aux1:R:1:aux2:R:2...
    !each property possessing a name (arbitrary length),
    !a one-letter format specifier (R=real, I=integer, S=string),
    !and the number of columns this property is using,
    !each of these three fields being separated by a colon (:).
    !Note that atomsk can store only numbers (real or integers)
    !as auxiliary properties, the rest will be ignored.
    DO WHILE( LEN_TRIM(temp(1:1)).NE.0 )
      strlength = SCAN(temp(1:),":")
      temp2 = temp(1:strlength-1)
      CALL ATOMSK_MSG(999,(/TRIM(temp2)/),(/0.d0/))
      IF(temp2=="species") THEN
        msg = 'EXYZ: found species'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        !this is no auxiliary property; read the format for species
        temp = temp(strlength+1:)
        READ(temp(1:3),*) temp2
        IF(TRIM(ADJUSTL(temp2)).NE."S:1") THEN
          !There is a problem with the format of species
          nerr=nerr+1
          GOTO 1000
        ELSE
          temp = temp(5:)
        ENDIF
        !
      ELSEIF(temp2=="pos") THEN
        msg = 'EXYZ: found positions'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        !this is no auxiliary property; read the format for positions
        temp = temp(strlength+1:)
        READ(temp(1:3),*) temp2
        IF(TRIM(ADJUSTL(temp2)).NE."R:3") THEN
          !There is a problem with the format of positions
          nerr=nerr+1
          GOTO 1000
        ELSE
          temp = temp(5:)
        ENDIF
        !
      ELSE
        msg = 'EXYZ: auxiliary properties: '//TRIM(temp2)
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        !else there are auxiliary properties
        !First count how many properties are real or integer numbers
        !This is because we can write only real numbers in AUX
        auxiliary=0
        temp2 = temp
        DO k=1, LEN_TRIM(temp2)
          IF( temp2(k:k+2)==":R:" .OR. temp2(k:k+2)==":I:" ) THEN
            !read the number of columns used by this auxiliary property
            READ(temp2(k+3:k+3),*,ERR=800,END=800) j
            auxiliary=auxiliary+j
          ENDIF
        ENDDO
        !
        WRITE(msg,*) 'EXYZ: number of aux. prop.: ', auxiliary
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        !
        IF( auxiliary.NE.0 ) THEN
          !Allocate arrays for auxiliary properties
          ALLOCATE( AUXNAMES(auxiliary) )
          AUXNAMES(:) = ""
          ALLOCATE( AUX(NP,auxiliary) )
          AUX(:,:) = 0.d0
          !
          j=0
          DO WHILE (j<auxiliary)
            j=j+1
            !Read name of this property
            strlength = SCAN(temp,":")
            IF( strlength.NE.0 ) THEN
              temp2 = temp(1:strlength-1)  !name of the property
              msg = 'EXYZ: auxiliary property: '//TRIM(temp2)
              CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
              !Special case for atom charges: transform it to q
              IF( TRIM(temp2(1:6))=="charge" ) temp2="q"
              !Read how many columns it uses
              READ(temp(strlength+3:strlength+3),*,ERR=830,END=830) l
              IF( l==1 ) THEN
                AUXNAMES(j) = TRIM(temp2)
              ELSEIF( l==3 ) THEN
                !We assume it is 3 properties along X, Y and Z
                !Special case for velocities: transform it to vx, vy, vz
                IF( TRIM(temp2)=="vel" ) temp2="v"
                AUXNAMES(j) = TRIM(temp2)//"x"
                AUXNAMES(j+1) = TRIM(temp2)//"y"
                AUXNAMES(j+2) = TRIM(temp2)//"z"
                j=j+3
              ELSE
                DO k=1,l
                  AUXNAMES(j+k) = TRIM(temp2)
                ENDDO
              ENDIF
            ENDIF
            !
            !Remove the property that was just read from "temp" for next iteration
            temp = temp(strlength+5:)
          ENDDO
          !
        ENDIF
        !
        !All properties were read => wipe out temp & get out of the loop
        temp = ""
        EXIT
        !
      ENDIF
    ENDDO !end loop on auxiliary properties
  ENDIF
ENDDO !end parsing the 2nd line
!
!
210 CONTINUE
!If both extended and special XYZ were detected we are in trouble
IF(exyz.AND.sxyz) GOTO 800
!
!If the cell vectors are not set, parse the second line again,
!and attempt to detect 3 consecutive real numbers
!Note: this assumes an orthorhombic box, i.e. only H(1,1), H(2,2) and H(3,3)
IF( .NOT.Hset ) THEN
  temp = properties
  !
  !Values may be separated by special characters, replace them with blank spaces
  DO i=1,7
    IF(i==1) spchar=":"
    IF(i==2) spchar=";"
    IF(i==3) spchar=","
    IF(i==4) spchar="'"
    IF(i==5) spchar='"'
    IF(i==6) spchar='('
    IF(i==7) spchar=')'
    DO WHILE( SCAN(temp,spchar)>0 )
      strlength = SCAN(temp,spchar)
      temp(strlength:strlength) = " "
    ENDDO
  ENDDO
  !
  !Look for 3 consecutive real numbers and save them as cell vectors
  P1=0.d0
  P2=0.d0
  P3=0.d0
  DO WHILE( LEN_TRIM(temp)>0 )
    READ(temp,*,END=212,ERR=212) temp2
    !Something was read, is it a real number?
    READ(temp2,*,END=211,ERR=211) tempreal
    !It was a real number, check that it is positive
    IF( tempreal > 0.d0 ) THEN
      P1 = tempreal
      !Try to read the next 2 real numbers
      strlength = SCAN(temp,' ')
      temp = ADJUSTL(temp(strlength+1:))
      READ(temp,*,END=212,ERR=212) temp2
      !Something was read, is it a real number?
      READ(temp2,*,END=211,ERR=211) tempreal
      !It was a real number, check that it is positive
      IF( tempreal > 0.d0 ) THEN
        P2 = tempreal
        !Try to read the third real number
        strlength = SCAN(temp,' ')
        temp = ADJUSTL(temp(strlength+1:))
        READ(temp,*,END=212,ERR=212) temp2
        !Something was read, is it a real number?
        READ(temp2,*,END=211,ERR=211) tempreal
        !It was a real number, check that it is positive
        IF( tempreal > 0.d0 ) THEN
          P3 = tempreal
          !We have 3 real numbers, no need to continue
          GOTO 212
        ENDIF
      ENDIF
    ENDIF
    !
    211 CONTINUE
    strlength = SCAN(temp,' ')
    temp = ADJUSTL(temp(strlength+1:))
    !
  ENDDO
  !
  212 CONTINUE
  IF( P1>0.d0 .AND. P2>0.d0 .AND. P3>0.d0 ) THEN
    H(1,1) = P1
    H(2,2) = P2
    H(3,3) = P3
    Hset = .TRUE.
  ENDIF
ENDIF
!
220 CONTINUE
IF( .NOT.exyz ) THEN
  ALLOCATE(comment(1))
  comment(1) = properties(1:128)
ENDIF
!
!Determine if the format is (number x y z)  or  (species x y z)
msg = 'determine if format is massxyz or spiexyz'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
READ(30,*,ERR=830,END=830) temp
temp=TRIM(ADJUSTL(temp))
!If temp(1:2) contains a species name, converting it will give tempreal.NE.0
CALL ATOMNUMBER(temp(1:2),tempreal)
!
IF(tempreal.NE.0.d0) THEN
  !A species name was read: it is "speciesxyz" format
  xyzformat = 'spiexyz'
!
ELSE
  !Otherwise the file should have "massxyz" format
  xyzformat = 'massxyz'
  !Let's confirm it by reading the real number
  !If there is an error then we are screwed -> exit
  READ(temp,*,ERR=830,END=830) tempreal
  !Then convert this real to atom species
  CALL ATOMSPECIES(tempreal,species)
  !And convert this species back to a real
  CALL ATOMNUMBER(species,tempreal)
  !Now the real should be the atomic number
  !If it is still zero then we are in trouble -> exit
  IF(tempreal==0.d0) GOTO 830
ENDIF
!
!
250 CONTINUE
WRITE(msg,*) 'format: ', xyzformat
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!Go back one line (because the first line was already read above)
BACKSPACE(UNIT=30)
!Read atomic positions
WRITE(msg,'(a36,i7)') 'read atomic positions; NP=', NP
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
DO i=1,NP
  READ(30,'(a1024)',ERR=800,END=800) temp
  IF(xyzformat=='massxyz') THEN
    READ(temp,*,ERR=800,END=800) P(i,4), P1, P2, P3
  ELSEIF(xyzformat=='spiexyz') THEN
    READ(temp,*,ERR=800,END=800) temp2, P1, P2, P3
    READ(temp2,*) species
    CALL ATOMNUMBER(species,P(i,4))
  ENDIF
  P(i,1) = P1
  P(i,2) = P2
  P(i,3) = P3
  !Read auxiliary properties if any
  IF( ALLOCATED(AUX) ) THEN
    READ(temp,*,ERR=800,END=800) temp2, P1, P2, P3, ( AUX(i,j), j=1,SIZE(AUX,2) )
  ENDIF
ENDDO
!
!
!
300 CONTINUE
!Find out if coordinates are reduced or cartesian
CALL FIND_IF_REDUCED(H,P,isreduced)
!In case of reduced coordinates, convert them to cartesian
IF(isreduced .AND. Hset) THEN
  CALL FRAC2CART(P,H)
ENDIF
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
GOTO 1000
!
810 CONTINUE
CALL ATOMSK_MSG(1809,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
820 CONTINUE
CALL ATOMSK_MSG(1810,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
830 CONTINUE
CALL ATOMSK_MSG(1804,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
1000 CONTINUE
CLOSE(30)
!
END SUBROUTINE READ_XYZ
!
!
END MODULE in_xyz
