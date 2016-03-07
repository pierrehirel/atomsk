MODULE in_cfg
!
!**********************************************************************************
!*  IN_CFG                                                                        *
!**********************************************************************************
!* This modules reads Ju Li's CFG atomic position files, both standard            *
!* and extended. These two CFG formats are described here:                        *
!*    http://mt.seas.upenn.edu/Archive/Graphics/A/#standard_CFG                   *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 02 March 2016                                    *
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
USE comv
USE constants
USE functions
USE messages
USE files
USE subroutines
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE READ_CFG(inputfile,H,P,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=1024):: temp
LOGICAL:: Hset, novelocity
INTEGER:: auxiliary  !number of auxiliary properties
INTEGER:: i, Hread
INTEGER:: Ncomment  !number of comment lines
INTEGER:: NP, NPcount, strlength, strlength2
REAL(dp):: a0, snumber
REAL(dp):: testreal, cfgformat
REAL(dp):: P1, P2, P3
REAL(dp), DIMENSION(3,3),INTENT(OUT):: H   !Base vectors of the supercell
REAL(dp), DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: P   !atom positions
REAL(dp), DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: AUX !auxiliary properties
!
!
!Initialize some variables
Hset = .FALSE.
novelocity = .FALSE.
a0 = 1.d0
auxiliary = 0
Ncomment=0
NP=0
Hread=0
NPcount=0
testreal = 0.d0
H(:,:) = 0.d0
 cfgformat = 0.05d0 ! negative = standard CFG
                    ! positive = extended CFG (assumed to be default)
!
!
!
100 CONTINUE
msg = 'entering READ_CFG'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
!Read the preamble of the file
DO WHILE(.NOT.Hset .OR. NP==0)
  READ(30,'(a128)',ERR=400,END=160) temp
  temp = ADJUSTL(temp)
  IF(temp(1:1).NE.'#') THEN
    !Check if it is a new species
    IF(temp(1:21)=='Number of particles =' .OR. temp(1:20)=='Number of particles=') THEN
        READ(temp(22:128),*,ERR=420,END=420) NP
        ALLOCATE(P(NP,4))
    ELSEIF(temp(1:2)=='A=' .OR. temp(1:3)=='A =') THEN
        READ(temp(4:128),*,ERR=400,END=400) testreal
        a0 = testreal
    ELSEIF(temp(1:7)=='H0(1,1)') THEN
        READ(temp(10:128),*,ERR=410,END=410) testreal
        H(1,1) = testreal
        Hread=Hread+1
    ELSEIF(temp(1:7)=='H0(1,2)') THEN
        READ(temp(10:128),*,ERR=410,END=410) testreal
        H(1,2) = testreal
        Hread=Hread+1
    ELSEIF(temp(1:7)=='H0(1,3)') THEN
        READ(temp(10:128),*,ERR=410,END=410) testreal
        H(1,3) = testreal
        Hread=Hread+1
    ELSEIF(temp(1:7)=='H0(2,1)') THEN
        READ(temp(10:128),*,ERR=410,END=410) testreal
        H(2,1) = testreal
        Hread=Hread+1
    ELSEIF(temp(1:7)=='H0(2,2)') THEN
        READ(temp(10:128),*,ERR=410,END=410) testreal
        H(2,2) = testreal
        Hread=Hread+1
    ELSEIF(temp(1:7)=='H0(2,3)') THEN
        READ(temp(10:128),*,ERR=410,END=410) testreal
        H(2,3) = testreal
        Hread=Hread+1
    ELSEIF(temp(1:7)=='H0(3,1)') THEN
        READ(temp(10:128),*,ERR=410,END=410) testreal
        H(3,1) = testreal
        Hread=Hread+1
    ELSEIF(temp(1:7)=='H0(3,2)') THEN
        READ(temp(10:128),*,ERR=410,END=410) testreal
        H(3,2) = testreal
        Hread=Hread+1
    ELSEIF(temp(1:7)=='H0(3,3)') THEN
        READ(temp(10:128),*,ERR=410,END=410) testreal
        H(3,3) = testreal
        Hread=Hread+1
    ENDIF
    !
    IF(Hread==9) Hset=.TRUE.
  !
  ELSE  !i.e. if the line starts with #
    Ncomment=Ncomment+1
  ENDIF  !Endif #
  !
  150 CONTINUE
  msg = 'CFG header: '//TRIM(temp)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
ENDDO
!
160 CONTINUE
H(:,:) = a0*H(:,:)
!
170 CONTINUE
!Determine if it is standard or extended CFG
!By default, extended CFG format is assumed
!The file is searched for lines with format (mass species x y z)
!If more than two such lines are found, the counter "cfgformat"
!becomes negative and the file is assumed to be standard CFG
DO i=1,500  !Limit the search to the first 500 lines
  READ(30,'(a128)',ERR=200,END=200) temp
  !Check if the line has the format (mass species x y z)
  READ(temp,*,ERR=171,END=171) testreal, species, testreal, testreal, testreal
  !If there was no error it is most likely standard CFG
  !Check it by trying to convert the species
  CALL ATOMNUMBER(species,snumber)
  IF(snumber.NE.0.d0) THEN
    !If conversion works then decrease the counter
    cfgformat = cfgformat-0.1d0
  ENDIF
  !
  171 CONTINUE
  !If the line did not have the format (mass species x y z)
  !then do nothing
ENDDO
!
!
!
200 CONTINUE
IF(Ncomment>0) THEN
  ALLOCATE(comment(Ncomment))
  comment(:)=''
ENDIF
!
IF(cfgformat<0.d0) THEN
  msg = 'standard'
ELSE
  msg = 'extended'
ENDIF
WRITE(msg,*) TRIM(msg)//' CFG: cfgformat = ', cfgformat
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
IF(cfgformat==0.d0) THEN
  nerr=nerr+1
  CALL ATOMSK_MSG(1804,(/''/),(/0.d0/))
  GOTO 500
ENDIF
!
!Go back to the beginning of the file
REWIND(30)
Ncomment = 0
P(:,:) = 0.d0
snumber = 0.d0
msg = 'Finding the starting point for atomic positions...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!Let's make sure to be placed at the beginning of atomic positions
!In both standard and extended CFG it is the first line that
!starts with a number
DO WHILE(snumber==0.d0)
  READ(30,'(a1024)',ERR=400,END=400) temp
  temp = ADJUSTL(temp)
  IF( LEN_TRIM(temp)>0 ) THEN
    !Ignore empty lines
    IF(temp(1:1)=='#') THEN
      Ncomment=Ncomment+1
      IF( Ncomment<=SIZE(comment) ) THEN
        comment(Ncomment) = temp(1:128)
      ENDIF
    ELSEIF(temp(1:13)=='.NO_VELOCITY.') THEN
      novelocity = .TRUE.
    ELSEIF(temp(1:11)=='entry_count') THEN
      strlength = SCAN(temp,'=')
      READ(temp(strlength+1:),*,ERR=400,END=400) auxiliary
      IF(novelocity) THEN
        auxiliary = auxiliary-3
      ELSE
        auxiliary = auxiliary-6
      ENDIF
      WRITE(msg,*) 'Number of auxiliary properties:', auxiliary
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      IF(auxiliary>0) THEN
        ALLOCATE( AUXNAMES(auxiliary) )
        AUXNAMES(:) = ''
        ALLOCATE( AUX( SIZE(P(:,1)), auxiliary ) )
        AUX(:,:) = 0.d0
      ENDIF
    ELSEIF( temp(1:9)=='auxiliary' ) THEN
      IF(.NOT. ALLOCATED(AUX) ) THEN
        !If auxiliary array was not allocated before we are in trouble
        nwarn = nwarn+1
        CALL ATOMSK_MSG(1700,(/TRIM(temp)/),(/0.d0/))
      ELSE
        !Read the name of the auxiliary property
        strlength = SCAN(temp,'[')
        strlength2 = SCAN(temp,']')
        READ(temp(strlength+1:strlength2-1),*) auxiliary
        strlength = SCAN(temp,'=')
        READ(temp(strlength+1:),'(a128)') AUXNAMES(auxiliary+1)
        AUXNAMES(auxiliary+1) = TRIM(ADJUSTL(AUXNAMES(auxiliary+1)))
      ENDIF
    !
    ELSEIF( temp(1:9)=="Transform" .OR. temp(1:3)=="eta" ) THEN
      !ignore that
      !
    ELSE
      !try to read a real
      READ(temp,*,ERR=201,END=201) testreal
      !No error? Then the real was most probably an atom mass.
      !To confirm that, let's try to read the atom species in the
      !second position (standard CFG) or the next line (extended CFG).
      !Note: ideally we should compare the atom mass and species
      !and require that they match. However some atomistic simulation
      !programs write the CFG format in a non-standard way, writing
      !the atomic number or a dummy number instead of the atomic mass.
      !So here we use only the atomic species
      IF(cfgformat<0.d0) THEN
        READ(temp,*,ERR=400,END=400) testreal, species
        CALL ATOMNUMBER(species,snumber)
      ELSE
        READ(30,*,ERR=400,END=400) species
        CALL ATOMNUMBER(species,snumber)
      ENDIF
      !If species was not recognized
      201 CONTINUE
    ENDIF
    !
  ENDIF
ENDDO
!
!Now, read all atomic coordinates
msg = 'Reading atomic positions'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
NPcount = 0
IF(cfgformat<0.d0) THEN
  !Standard CFG format
  !Go back one line
  BACKSPACE(30)
  !Read atom information
  DO WHILE(NPcount<NP)
    !Read the line
    READ(30,'(a1024)',ERR=300,END=300) temp
    temp = TRIM(ADJUSTL(temp))
    !
    IF( LEN_TRIM(temp)>0 .AND. temp(1:1).NE.'#' ) THEN
      !First column = atom mass (ignored)
      !Second column = atom species (converted into atomic number)
      !3rd, 4th, 5th columns = X, Y, Z
      READ(temp,*,ERR=400,END=400) testreal, species, P1, P2, P3
      CALL ATOMNUMBER(species,snumber)
      NPcount=NPcount+1
      P(NPcount,1) = P1
      P(NPcount,2) = P2
      P(NPcount,3) = P3
      P(NPcount,4) = snumber
    ENDIF
  ENDDO
!
ELSE
  !Extended CFG format
  DO WHILE(NPcount<NP)
    !Read the line
    READ(30,'(a1024)',ERR=300,END=240) temp
    240 CONTINUE
    temp = TRIM(ADJUSTL(temp))
    !
    IF( LEN_TRIM(temp)>0 .AND. temp(1:1).NE.'#' ) THEN
      READ(temp,*,ERR=250,END=250) P1, P2, P3
      !
      !No error=we just read the 3 coordinates
      !Load these coordinates in P
      NPcount=NPcount+1
      P(NPcount,1) = P1
      P(NPcount,2) = P2
      P(NPcount,3) = P3
      P(NPcount,4) = snumber
      !Read auxiliary property if any
      IF(ALLOCATED(AUX)) THEN
        READ(temp,*,ERR=400,END=400) P1, P2, P3, &
            & ( AUX(NPcount,auxiliary), auxiliary=1,SIZE(AUX(1,:)) )
      ENDIF
      !Success, jump to the next line of the CFG file
      GOTO 260
      !
      250 CONTINUE
      !If we did not read 3 coordinates then it has to be the
      !lines indicating a new species; let's check that
      READ(temp,*,ERR=400,END=400) testreal   !First line is atomic mass
      READ(30,'(a64)',ERR=400,END=400) temp   !Second line is atomic species
      species = TRIM(ADJUSTL(temp))
      CALL ATOMNUMBER(species,snumber)
    ENDIF
    260 CONTINUE
  ENDDO
ENDIF
!
!
!
300 CONTINUE
IF( NPcount.NE.SIZE(P(:,1)) ) THEN
  nerr=nerr+1
  CALL ATOMSK_MSG(1805,(/''/),(/ DBLE(NPcount), DBLE( SIZE(P(:,1)) ) /))
ELSE
  CLOSE(30)
  !Convert coordinates to cartesian
  CALL FRAC2CART(P,H)
ENDIF
GOTO 500
!
!
!
400 CONTINUE
CLOSE(30)
nerr = nerr+1
CALL ATOMSK_MSG(802,(/''/),(/DBLE(NPcount)/))
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
!
!
500 CONTINUE
!
END SUBROUTINE READ_CFG
!
!
END MODULE in_cfg
