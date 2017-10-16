MODULE in_gulp_gin
!
!**********************************************************************************
!*  IN_GULP_GIN                                                                   *
!**********************************************************************************
!* This module reads GULP input files (.gin) and restart files (.res or .grs).    *
!* This format is described in the GULP manual, which can be found here:          *
!*    https://www.ivec.org/gulp/help/manuals.html                                 *
!* WARNING: this file format is more or less permissive, and people use it        *
!*         in a more or less stupid way, sometimes using incomprehensible         *
!*         atom names, omitting some columns, etc. As a result, this module       *
!*         is more or less a mess, but it tries to read the file no matter what.  *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 16 Oct. 2017                                     *
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
!
CONTAINS
!
SUBROUTINE READ_GIN(inputfile,H,P,S,comment,AUXNAMES,AUX)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!RESTRICTIONS OF THIS ROUTINE:
!1- Symmetry operations are not taken into account. So, if
!   you specify coordinates in a primitive lattice, only
!   atoms in that primitive cell will be read and output.
!   Symmetry is a whole other story to implement, so it was
!   not done in this program.
!
!2- since NEB was implemented in GULP (>3.4), the RES file
!   can contain all atomic positions for all images of NEB.
!   This is not recognized here, and only the positions
!   of the 1st NEB image will be read and converted by
!   this program.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=4):: coord, ptype
CHARACTER(LEN=128):: sgroup  !Hermann-Mauguin symbol or number of space group
CHARACTER(LEN=128):: coord1, coord2, coord3
CHARACTER(LEN=128):: test, test2, test3, temp, temp2
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(100):: tempcomment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
LOGICAL:: chargesC, chargesS  !are charges of cores/shells defined in the "species" section?
LOGICAL:: knowcol      !do we know the number of columns on each line?
LOGICAL:: shellsafter  !are shell positions after all core positions?
LOGICAL:: velocities !are velocities present?
INTEGER:: fixx, fixy, fixz !columns in AUX where the fixes are defined
INTEGER:: q, qs  !columns in AUX where charges of cores/shells are defined
INTEGER:: occ, radius !columns in AUX where where occupancy and radius are defined
INTEGER:: vx, vy, vz !columns in AUX where velocities are defined
INTEGER:: i, j, k
INTEGER:: icoordtemp1, icoordtemp3
INTEGER:: Ncol, icol  !number of columns for each core (or shell), index of volumn
INTEGER:: NP, NS  !number of cores, shells
INTEGER:: sgroupnum
INTEGER:: strlength
REAL(dp):: alpha, beta, gamma
REAL(dp):: Hx, Hy, Hz !Length of H(1,:), H(2,:) and H(3,:)
REAL(dp):: P1, P2, P3
REAL(dp):: snumber !Mass of atoms
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(30,3):: tempq  !temp. array for core/shell charges (max.30 species)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, S !positions of cores, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
sgroup = ''
test=''
tempcomment(:)=''
 coord = 'cart'
 chargesC = .FALSE.
 chargesS = .FALSE.
knowcol=.FALSE.
shellsafter = .FALSE.
velocities = .FALSE.
Ncol = 0
icol = 0
fixx = 0
fixy = 0
fixz = 0
i=0
j=0
occ = 0
q = 0
qs = 0
radius = 0
sgroupnum = 0
vx = 0
vy = 0
vz = 0
H(:,:) = 0.d0
tempq(:,:) = 0.d0
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
!
!
10 CONTINUE
msg = 'entering READ_GIN'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
!
20 CONTINUE
!
!Parse the file to find and read comment, supercell vectors,
!and whether velocities are present
DO
  READ(30,'(a128)',ERR=25,END=25) temp
  temp = TRIM(ADJUSTL(temp))
  IF(temp(1:5)=='title') THEN
    !Read all comment lines (up to 100 are supported here)
    k=0
    DO WHILE( temp(1:3).NE.'end' .AND. k<100)
      k=k+1
      READ(30,'(a128)',ERR=600,END=600) temp
      tempcomment(k) = TRIM(ADJUSTL(temp))
    ENDDO
    !Save comments to array comment(:)
    ALLOCATE(comment(MAX(1,k-1)))
    DO j=1,k-1
      comment(j) = tempcomment(j)
    ENDDO
    !
  ELSEIF(temp(1:3)=='vec') THEN
    msg = 'Found cartesian vectors.'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    READ(30,*,END=820,ERR=820) H(1,1), H(1,2), H(1,3)
    READ(30,*,END=820,ERR=820) H(2,1), H(2,2), H(2,3)
    READ(30,*,END=820,ERR=820) H(3,1), H(3,2), H(3,3)
    !
  ELSEIF(temp(1:3)=='cel') THEN
    msg = 'Found conventional vectors:'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    READ(30,*,END=820,ERR=820) Hx, Hy, Hz,             &
        &                      alpha, beta, gamma
    !
    IF(verbosity==4) THEN
      WRITE(msg,21) '  a = ', Hx, '     alpha = ', alpha
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      WRITE(msg,21) '  b = ', Hy, '      beta = ', beta
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      WRITE(msg,21) '  c = ', Hz, '     gamma = ', gamma
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDIF
    !
    !In GULP angles are in degrees: convert to radians
    alpha = DEG2RAD(alpha)
    beta  = DEG2RAD(beta)
    gamma = DEG2RAD(gamma)
    !Convert basis vectors in cartesian coordinates
    CALL CONVMAT(Hx,Hy,Hz,alpha,beta,gamma,H)
    !
  ELSEIF( temp(1:3)=="vel" ) THEN
    !Velocities are present
    velocities = .TRUE.
    vx=1
    vy=2
    vz=3
    icol = 3
    msg = 'Found velocities.'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !
  ELSEIF( temp(1:7)=="species" ) THEN
    !Read the charges of cores and shells
    !(in case of error just ignore this part)
    !  tempq(:,1) = atomic number
    !  tempq(:,2) = core charge
    !  tempq(:,3) = shell charge (if any)
    READ(temp(8:),*,END=20,ERR=20) k
    IF( k<=SIZE(tempq) ) THEN
      tempq(:,:) = 0.d0
      DO i=1,k
        READ(30,*,END=20,ERR=20) test, test2, test3
        !Determine ion species
        READ(test,*,END=20,ERR=20) species
        CALL ATOMNUMBER(species,snumber)
        IF(NINT(snumber)<=0) GOTO 20
        !Read the type of entry: "core" or "shel"
        test2 = ADJUSTL(test2)
        ptype = test2(1:4)
        !Read the electric charge
        READ(test3,*,END=20,ERR=20) P1
        !Parse tempq and store charge in the appropriate row
        j=0
        DO WHILE ( LEN_TRIM(ptype)>0 )
          j=j+1
          IF( NINT(tempq(j,1)) == NINT(snumber) ) THEN
            !Ion #i is the same species as ion #j
            !Save the new charge into tempq
            IF( ptype=="shel" ) THEN
              tempq(j,3) = P1
              chargesS = .TRUE.
            ELSE
              tempq(j,2) = P1
              chargesC = .TRUE.
            ENDIF
            ptype=""
          ELSEIF( NINT(tempq(j,1)) <= 0 ) THEN
            !Ion #i belongs to a species never seen before
            !=> create a new entry
            tempq(j,1) = snumber
            IF( ptype=="shel" ) THEN
              tempq(j,3) = P1
              chargesS = .TRUE.
            ELSE
              tempq(j,2) = P1
              chargesC = .TRUE.
            ENDIF
            ptype=""
          ENDIF
        ENDDO
      ENDDO
      !
      IF( verbosity==4 ) THEN
        DO i=1,j
          WRITE(msg,*) 'tempq:', NINT(tempq(i,1)), tempq(i,2), tempq(i,3)
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        ENDDO
      ENDIF
      !
    ENDIF
    !
  ELSEIF( temp(1:10)=="spacegroup" ) THEN
    !Read the space group
    READ(30,'(a128)',END=820,ERR=820) sgroup
    !Check if it is an integer or not
    READ(sgroup,*,ERR=22,END=22) sgroupnum
    ! Success => go directly to 23
    GOTO 23
    22 CONTINUE
    ! There was an error => sgroup does not contain an integer
    ! It should contain a Hermann-Mauguin symbol
    sgroup = ADJUSTL(sgroup)
    ! In GULP format, the symbol may upper case letters and blanck spaces
    ! Transform all that into a character chain that is valid for Atomsk routines
    ! 1. Remove spaces if any
    j = SCAN(TRIM(sgroup),' ')
    DO WHILE( j>LEN_TRIM(sgroup) )
      sgroup(j:) = sgroup(j+1:)
      j = SCAN(TRIM(sgroup),' ')
    ENDDO
    ! 2. Replace upper case letters with lower case (except first letter)
    sgroup(2:) = StrDnCase(sgroup(2:))
    !
    msg = "space group name: "//TRIM(sgroup)
    !
    ! Verify if it is a valid space group name
    CALL SG_NAMGETNUM(TRIM(ADJUSTL(sgroup)),sgroupnum)
    IF( sgroupnum <= 0 .OR. sgroupnum > 230 ) THEN
      ! Invalid space group number
      nerr = nerr+1
      CALL ATOMSK_MSG(809,(/TRIM(sgroup)/),(/0.d0/))
      GOTO 1000
    ENDIF
  ENDIF
  23 CONTINUE
ENDDO
21 FORMAT(a6,f9.5,a13,f9.5)
!
!
25 CONTINUE
IF(verbosity==4) THEN
  msg = 'Base vectors:'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,27) '     | ', H(1,1), H(1,2), H(1,3), '|'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,27) ' H = | ', H(2,1), H(2,2), H(2,3), '|'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,27) '     | ', H(3,1), H(3,2), H(3,3), '|'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  27 FORMAT(a7,3(f10.5,2X),a1)
ENDIF
!
!Return to beginning of file, because sometimes
!the vectors are written after the coordinates
REWIND(30)
!
26 CONTINUE
!Find out if coordinates are fractional or cartesian
test = ''
DO WHILE( test(1:3).NE."fra" .AND. test(1:3).NE.'car' )
  READ(30,*,ERR=600,END=600) test
  test = TRIM(ADJUSTL(test))
  IF(test(1:3)=='fra') THEN
    msg = 'Found fractional coordinates.'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    coord='frac'
    GOTO 30
  ELSEIF(test(1:3)=='car') THEN
    msg = 'Found cartesian coordinates.'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    coord='cart'
    GOTO 30
  ENDIF
ENDDO
!
!
!
30 CONTINUE
NP=0
NS=0
!Count the total number of atoms NP
DO
  READ(30,*,ERR=31,END=31), test, test2
  test = ADJUSTL(test)
  IF(test(1:1).NE.'#') THEN
    !In all cases if the first two characters do not correspond
    !to a recognizable species then we are done counting
    species = test(1:2)
    CALL ATOMNUMBER(species,snumber)
    IF( snumber<=0.1d0 ) THEN
      species = species(1:1)
      CALL ATOMNUMBER(species,snumber)
    ENDIF
    IF( DABS(snumber) < 0.1d0 ) GOTO 31
    !
    !If cores and/or shells exist, count them separately
    IF(TRIM(ADJUSTL(test2))=='core') THEN
      NP=NP+1
    ELSEIF(TRIM(ADJUSTL(test2))=='shel') THEN
      NS=NS+1
    ELSE
      !If there is only atom species (i.e. "core" or "shell" are not specified)
      !then consider it a core
      NP=NP+1
    ENDIF
    !If there is no 'shel' or 'core' on the line,
    !then we are done counting
    !IF( TRIM(ADJUSTL(test2)).NE.'core' .AND.                &
    !  & TRIM(ADJUSTL(test2)).NE.'shel'       ) GOTO 31
  ENDIF
ENDDO
!
31 CONTINUE
WRITE(msg,*) 'NP, NS = ', NP, NS
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!Allocate array for atom positions
IF(NP>0) THEN
  ALLOCATE(P(NP,4))
  P(:,:) = 0.d0
ELSE
  !If NP=0 we are in trouble => error message and exit
  CALL ATOMSK_MSG(804,(/''/),(/0.d0/))
  nerr=nerr+1
  GOTO 1000
ENDIF
!
IF(NS>0) THEN
  !Allocate array for shells (same size as P)
  ALLOCATE(S(NP,4))
  S(:,:) = 0.d0
ENDIF
!
REWIND(30)
!Look for beginning of atom positions
test=''
DO WHILE( .NOT.(test(1:3)=='fra' .OR. test(1:3)=='car') )
  READ(30,*,END=800,ERR=800), test
  test = ADJUSTL(test)
ENDDO
!
!
40 CONTINUE
!Read atom positions
species = ''
!
i=0  !Counter for cores
j=0  !Counter for shells
DO WHILE(i<NP .OR. j<NS)
  strlength = 0
  READ(30,'(a128)',ERR=800,END=800) temp
  temp = ADJUSTL(temp)
  !
  !If it is a comment then skip the whole thing and read next line
  IF(temp(1:1).NE.'#' .AND. LEN_TRIM(temp).NE.0) THEN
    READ(temp,*,END=800,ERR=800) test
    test = ADJUSTL(test)
    !Remove trailing comment if any
    strlength = SCAN(test,'#')
    IF( strlength>0 ) THEN
      test(strlength:) = ' '
    ENDIF
    !
    !Read atom species
    !Sometimes in GULP files there are strange "atom names"
    !like "O1", "O2" etc. => Try to determine actual species
    !from the first (or first two) characters
    species=test(1:2)
    CALL ATOMNUMBER(species,snumber)
    IF( snumber==0.d0 ) THEN
      !If we didn't get a proper atom species, try
      !with the first letter only
      species = species(1:1)
      CALL ATOMNUMBER(species,snumber)
      !If that succeeded, we're good.
      !Otherwise if it still didn't work it means that
      !we cannot understand the "atom name".
      !The output will be garbage but I don't care
    ENDIF
    !
    !Determine if it is a core or shell
    READ(temp,*,END=800,ERR=800) test, test2
    IF(TRIM(ADJUSTL(test2))=='core') THEN
      ptype = 'core'
      i=i+1
      !Save the line in "temp2" for later
      strlength = INDEX(temp,'core')+5
    ELSEIF(TRIM(ADJUSTL(test2))=='shel') THEN
      ptype = 'shel'
      IF( i>=NP .AND. j<1 ) THEN
        !All cores were read, but no shell yet, this is the first
        !=> all shells positions are after all core positions
        shellsafter = .TRUE.
      ENDIF
      IF( shellsafter ) THEN
        !This makes it difficult to find which shell
        !corresponds to which core. Here we attempt to find the next
        !core that has the same atom species as this shell.
        !This assumes that the following sequence of shells will
        !match the sequence of cores.
        DO k=j,NP
          IF( DABS(snumber-P(k,4))<0.1d0 ) THEN
            j=k
            EXIT
          ENDIF
        ENDDO
      ELSE
        !Just assume that current shell corresponds to current core
        j=i
      ENDIF
      !Save the line in "temp2" for later
      strlength = INDEX(temp,'shel')+5
    ELSE
      !Nothing specified: it is assumed to be a core
      ptype = 'core'
      i=i+1
      !GOTO 45
      !Save the line in "temp2" for later
      strlength = 3
    ENDIF
    temp2 = temp(strlength:)
    !
    !
    WRITE(msg,*) 'Read '//TRIM(ptype)//' coordinates: ', species, ' ', TRIM(temp2)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !
    IF(coord=='cart') THEN
    !Cartesian coordinates are easy to read
      READ(temp2,*,END=800,ERR=800) P1, P2, P3
      !Save the rest of the line in "temp2" for later
      temp2 = ADJUSTL(temp2)
      DO k=1,3
        READ(temp2,*,END=800,ERR=800) test
        strlength = LEN_TRIM(test)+1
        temp2 = ADJUSTL(temp2(strlength:))
      ENDDO
    !
    ELSEIF(coord=='frac') THEN
      !GULP sometimes outputs reduced coordinates as fractions,
      !eg. 1/3 or 5/6, which Fortran does not consider
      !as REAL... even worse, Fortran stops reading at the
      !slash sign... so we have to deal with that...
      READ(temp2,*,END=800,ERR=800) coord1
      IF( LEN_TRIM(ADJUSTL(coord1))==1 ) THEN
        !Length is 1 because there is a slash: try to read fraction
        strlength = strlength+SCAN(temp2,'0123456789')+2
        temp2 = ADJUSTL(temp2)
        READ(temp2(1:1),*,END=800,ERR=800) icoordtemp1
        READ(temp2(3:3),*,END=800,ERR=800) icoordtemp3
        IF(icoordtemp3==0) GOTO 810
        P1=DBLE(icoordtemp1)/DBLE(icoordtemp3)
        WRITE(msg,*) 'Fraction1:', icoordtemp1, '/', icoordtemp3
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ELSE
        !Easy case: "coord1" is just a real
        READ(coord1,*,END=800,ERR=800) P1
        strlength = strlength+SCAN(temp2,'0123456789')+LEN_TRIM(coord1)
      ENDIF
      !
      !same with Y coordinate
      temp2 = temp(strlength:)
      READ(temp2,*,END=800,ERR=800) coord2
      IF( LEN_TRIM(ADJUSTL(coord2))==1 ) THEN
        !Length is 1 because there is a slash: try to read fraction
        strlength = strlength+SCAN(temp2,'0123456789')+2
        temp2 = ADJUSTL(temp2)
        READ(temp2(1:1),*,END=800,ERR=800) icoordtemp1
        READ(temp2(3:3),*,END=800,ERR=800) icoordtemp3
        IF(icoordtemp3==0) GOTO 810
        P2=DBLE(icoordtemp1)/DBLE(icoordtemp3)
        WRITE(msg,*) 'Fraction2:', icoordtemp1, '/', icoordtemp3
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ELSE
        !Easy case: "coord2" is just a real
        READ(coord2,*,END=800,ERR=800) P2
        strlength = strlength+SCAN(temp2,'0123456789')+LEN_TRIM(coord2)
      ENDIF
      !
      !same with Z coordinate
      temp2 = temp(strlength:)
      READ(temp2,*,END=800,ERR=800) coord3
      IF( LEN_TRIM(ADJUSTL(coord3))==1 ) THEN
        !Length is 1 because there is a slash: try to read fraction
        strlength = strlength+SCAN(temp2,'0123456789')+2
        temp2 = ADJUSTL(temp2)
        READ(temp2(1:1),*,END=800,ERR=800) icoordtemp1
        READ(temp2(3:3),*,END=800,ERR=800) icoordtemp3
        IF(icoordtemp3==0) GOTO 810
        P3=DBLE(icoordtemp1)/DBLE(icoordtemp3)
        WRITE(msg,*) 'Fraction3:', icoordtemp1, '/', icoordtemp3
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ELSE
        !Easy case: "coord3" is just a real
        READ(coord3,*,END=800,ERR=800) P3
        strlength = strlength+SCAN(temp2,'0123456789')+LEN_TRIM(coord3)
      ENDIF
      temp2 = temp(strlength:)
    !
    ENDIF !endif coord=cart
    !
    !Finally, set up atomic species and load coord. to P (cores) or S (shells)
    !Note: at this point, when reading a shell then j=i if the shells
    !always appear just after their core, or if there are all cores positions
    !and then all shells positions in a matching order. If cores and shells
    !are written in random order, then the core-shell association will be wrong.
    !WRITE(msg,*) '       ', P1, P2, P3
    !CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    IF(ptype=='core') THEN
      IF( i>0 .AND. i<=SIZE(P,1) ) THEN
        P(i,4) = snumber
        P(i,1) = P1
        P(i,2) = P2
        P(i,3) = P3
      ELSE
        !out-of-bounds
        nwarn=nwarn+1
        CALL ATOMSK_MSG(2742,(/""/),(/DBLE(i)/))
      ENDIF
    ELSE
      IF( j>0 .AND. j<=SIZE(S,1) ) THEN
        S(j,4) = snumber
        S(j,1) = P1
        S(j,2) = P2
        S(j,3) = P3
      ELSE
        !out-of-bounds
        nwarn=nwarn+1
        CALL ATOMSK_MSG(2742,(/""/),(/DBLE(j)/))
      ENDIF
    ENDIF
    !
    !
    !Attempt to read auxiliary properties that follow on each line
    IF( .NOT.knowcol ) THEN
      !This is the first atom position that we read:
      !read the rest of the line and count how many columns exist after x y z.
      !Note that the columns cannot be counted earlier because of the possible
      !slash signs that may screw everything...
      WRITE(msg,*) 'Columns: ', TRIM(temp2)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      Ncol=0
      test = ADJUSTL(temp2)
      DO
        READ(test,*,ERR=43,END=43) test2
        IF( LEN_TRIM(test2)==0 ) THEN
          EXIT
        ELSE
          Ncol=Ncol+1
          strlength = LEN_TRIM(test2)+1
          test = ADJUSTL(test(strlength:))
        ENDIF
      ENDDO
      43 CONTINUE
      IF( verbosity==4) THEN
        WRITE(msg,*) 'Ncol = ', Ncol
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        WRITE(msg,*) 'charges cores =  ', chargesC
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        WRITE(msg,*) 'charges shells = ', chargesS
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        WRITE(msg,*) 'velocities = ', velocities
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDIF
      !
      !The number of columns fixes the format:
      !(only Ncol=3 is ambiguous, it will be dealt with later)
      !Ncol=0   (no auxiliary property)
      !Ncol=1   q
      !Ncol=2   q occ
      !Ncol=3   q occ radius
      !Ncol=3   fixx fixy fixz
      !Ncol=4   q fixx fixy fixz
      !Ncol=5   q occ fixx fixy fixz
      !Ncol=6   q occ radius fixx fixy fixz
      !
      !Note: atomsk stores the charges of cores as q, and shells as qs.
      !     So, if the number of shells is not zero, allocate a column for qs.
      !
      !At this point icol=3 if velocities were detected, 0 otherwise
      IF( Ncol.NE.0 ) THEN
        icol = icol+1
        IF( Ncol==1 ) THEN
          q = icol
          IF(NS>0) THEN
            icol=icol+1
            qs = icol
          ENDIF
          !
        ELSEIF( Ncol==2 ) THEN
          q = icol
          IF(NS>0) THEN
            icol=icol+1
            qs = icol
          ENDIF
          icol=icol+1
          occ = icol
          !
        ELSEIF( Ncol==3 ) THEN
          READ(temp2,*,ERR=800,END=800) test, test2, test3
          IF( (test=='1'.OR.test=='0') .AND. (test2=='1'.OR.test2=='0')  &
            & .AND. (test3=='1'.OR.test3=='0')                           ) THEN
            !It is the fix flags
            !However, charges of cores/shells may have been defined elsewhere
            !(e.g. in the "species" section)
            IF(chargesC) THEN
              q = icol
              icol=icol+1
            ENDIF
            IF(chargesS .AND. ALLOCATED(S) .AND. SIZE(S,1)>0 ) THEN
              qs = icol
              icol=icol+1
            ENDIF
            fixx = icol
            icol=icol+1
            fixy = icol
            icol=icol+1
            fixz = icol
          ELSE
            !It is the charge, occupancy, radius
            q = icol
            IF(NS>0) THEN
              icol=icol+1
              qs = icol
            ENDIF
            icol=icol+1
            occ = icol
            icol=icol+1
            radius = icol
          ENDIF
          !
        ELSEIF( Ncol==4 ) THEN
          q = icol
          IF(NS>0) THEN
            icol=icol+1
            qs = icol
          ENDIF
          icol=icol+1
          fixx = icol
          icol=icol+1
          fixy = icol
          icol=icol+1
          fixz = icol
          !
        ELSEIF( Ncol==5 ) THEN
          q = icol
          IF(NS>0) THEN
            icol=icol+1
            qs = icol
          ENDIF
          icol=icol+1
          occ = icol
          icol=icol+1
          fixx = icol
          icol=icol+1
          fixy = icol
          icol=icol+1
          fixz = icol
          !
        ELSEIF( Ncol==6 ) THEN
          q = icol
          IF(NS>0) THEN
            icol=icol+1
            qs = icol
          ENDIF
          icol=icol+1
          occ = icol
          icol=icol+1
          radius = icol
          icol=icol+1
          fixx = icol
          icol=icol+1
          fixy = icol
          icol=icol+1
          fixz = icol
          !
        ENDIF
        !
        !Allocate arrays: number of auxiliary properties = icol
        ALLOCATE( AUX(SIZE(P,1),icol) )
        AUX(:,:) = 0.d0
        ALLOCATE( AUXNAMES(icol) )
        !Set names of auxiliary properties
        IF(occ>0 .AND. occ<=SIZE(AUX,1)) AUXNAMES(occ) = 'occupancy'
        IF(radius>0 .AND. radius<=SIZE(AUX,1)) AUXNAMES(radius) = 'bsradius'
        IF(fixx>0 .AND. fixx<=SIZE(AUX,1)) AUXNAMES(fixx) = 'fixx'
        IF(fixy>0 .AND. fixy<=SIZE(AUX,1)) AUXNAMES(fixy) = 'fixy'
        IF(fixz>0 .AND. fixz<=SIZE(AUX,1)) AUXNAMES(fixz) = 'fixz'
        !
        IF( chargesC .OR. chargesS ) THEN
          !Charges of cores or shells were defined before (in section "species")
          !Define them for all 
          DO k=1,SIZE(AUX,1)
            
          ENDDO
        ENDIF
        !
        WRITE(msg,*) 'SIZE AUX = ', SIZE(AUX,1), SIZE(AUX,2)
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        !
      ELSE
        !No additional auxiliary property, each line only contains x y z
        !BUT charges and velocities may have been defined elsewhere
        !=> allocate arrays for that
        k=0
        IF( velocities ) THEN
          k=k+3
          vx=1
          vy=2
          vz=3
        ENDIF
        IF( chargesC ) THEN
          k=k+1
          q=k
        ENDIF
        IF( chargesS .AND. ALLOCATED(S) .AND. SIZE(S,1)>0 ) THEN
          k=k+1
          qs=k
        ENDIF
        IF(k>0) THEN
          ALLOCATE( AUX(SIZE(P,1),k) )
          AUX(:,:) = 0.d0
          ALLOCATE( AUXNAMES(k) )
        ENDIF
      ENDIF !endif Ncol.NE.0
      !
      !Set names of auxiliary properties
      IF( velocities ) THEN
        AUXNAMES(vx) = 'vx'
        AUXNAMES(vy) = 'vy'
        AUXNAMES(vz) = 'vz'
      ENDIF
      IF( q>0 .AND. q<=SIZE(AUXNAMES) ) THEN
        AUXNAMES(q) = 'q'
      ENDIF
      IF( qs>0 .AND. qs<=SIZE(AUXNAMES) ) THEN
        AUXNAMES(qs) = 'qs'
      ENDIF
      !
      !We know the number of columns: go back one line
      knowcol = .TRUE.
      BACKSPACE(30)
      i=i-1
      WRITE(msg,*) 'Ncol = ', Ncol
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      !
    ELSE
      !The number of columns is known: read properties
      !
      !If charges were defined before, set it for this atom
      !(note that it will be overriden later if q.NE.0)
      IF( chargesC .AND. q>0 ) THEN
        DO k=1,SIZE(tempq,1)
          IF( NINT(tempq(k,1))==NINT(P(i,4)) ) THEN
            AUX(i,q) = tempq(k,2)
            EXIT
          ENDIF
        ENDDO
      ENDIF
      IF( chargesS .AND. qs>0 ) THEN
        DO k=1,SIZE(tempq,1)
          IF( NINT(tempq(k,1))==NINT(P(i,4)) ) THEN
            AUX(i,qs) = tempq(k,3)
            EXIT
          ENDIF
        ENDDO
      ENDIF
      !
      !Read auxiliary properties depending on the number of columns in the file
      !Note: for GULP, a fix=1 means that the particle is mobile,
      !     while for atomsk fix=1 means that it is fixed, so this
      !     needs to be converted
      IF( Ncol>0 ) THEN
        IF( Ncol==1 ) THEN
          !Only charge
          IF(ptype=='shel' .AND. qs>0) THEN
            READ(temp2,*,END=800,ERR=800) AUX(i,qs)
          ELSE
            READ(temp2,*,END=800,ERR=800) AUX(i,q)
          ENDIF
        ELSEIF( Ncol==2 ) THEN
          !charge occ
          IF(ptype=='shel') THEN
            READ(temp2,*,END=800,ERR=800) AUX(i,qs), AUX(i,occ)
          ELSE
            READ(temp2,*,END=800,ERR=800) AUX(i,q), AUX(i,occ)
          ENDIF
        ELSEIF( Ncol==3 ) THEN
          !ambiguous case
          IF( fixx>0 .AND. fixy>0 .AND. fixz>0 ) THEN
            !Fixed coordinates fixx fixy fixz
            READ(temp2,*,END=800,ERR=800) AUX(i,fixx), AUX(i,fixy), AUX(i,fixz)
            !Replace 0 by 1 and vice-versa
            AUX(i,fixx) = DABS(AUX(i,fixx)-1.d0)
            AUX(i,fixy) = DABS(AUX(i,fixy)-1.d0)
            AUX(i,fixz) = DABS(AUX(i,fixz)-1.d0)
          ELSE
            !charge occ radius
            IF(ptype=='shel') THEN
              READ(temp2,*,END=800,ERR=800) AUX(i,qs), AUX(i,occ), AUX(i,radius)
            ELSE
              READ(temp2,*,END=800,ERR=800) AUX(i,q), AUX(i,occ), AUX(i,radius)
            ENDIF
          ENDIF
        ELSEIF( Ncol==4 ) THEN
          !charge fixx fixy fixz
          IF(ptype=='shel') THEN
            READ(temp2,*,END=800,ERR=800) AUX(i,qs), AUX(i,fixx),               &
                                        & AUX(i,fixy), AUX(i,fixz)
          ELSE
            READ(temp2,*,END=800,ERR=800) AUX(i,q), AUX(i,fixx),                &
                                        & AUX(i,fixy), AUX(i,fixz)
          ENDIF
          !Replace 0 by 1 and vice-versa
          AUX(i,fixx) = DABS(AUX(i,fixx)-1.d0)
          AUX(i,fixy) = DABS(AUX(i,fixy)-1.d0)
          AUX(i,fixz) = DABS(AUX(i,fixz)-1.d0)
        ELSEIF( Ncol==5 ) THEN
          !charge occ fixx fixy fixz
          IF(ptype=='shel') THEN
            READ(temp2,*,END=800,ERR=800) AUX(i,qs), AUX(i,occ),                &
                                        & AUX(i,fixx), AUX(i,fixy), AUX(i,fixz)
          ELSE
            READ(temp2,*,END=800,ERR=800) AUX(i,q), AUX(i,occ),                 &
                                        & AUX(i,fixx), AUX(i,fixy), AUX(i,fixz)
          ENDIF
          AUX(i,fixx) = DABS(AUX(i,fixx)-1.d0)
          AUX(i,fixy) = DABS(AUX(i,fixy)-1.d0)
          AUX(i,fixz) = DABS(AUX(i,fixz)-1.d0)
        ELSEIF( Ncol==6 ) THEN
          !charge occ radius fixx fixy fixz
          IF(ptype=='shel') THEN
            READ(temp2,*,END=800,ERR=800) AUX(i,qs), AUX(i,occ), AUX(i,radius), &
                                        & AUX(i,fixx), AUX(i,fixy), AUX(i,fixz)
          ELSE
            READ(temp2,*,END=800,ERR=800) AUX(i,q), AUX(i,occ), AUX(i,radius),  &
                                        & AUX(i,fixx), AUX(i,fixy), AUX(i,fixz)
          ENDIF
          AUX(i,fixx) = DABS(AUX(i,fixx)-1.d0)
          AUX(i,fixy) = DABS(AUX(i,fixy)-1.d0)
          AUX(i,fixz) = DABS(AUX(i,fixz)-1.d0)
        ENDIF
        !Any other number of columns: don't attempt to read aux. properties
      ENDIF
      !
    ENDIF
    !
  ENDIF
  !
  45 CONTINUE
ENDDO
!
!
!
500 CONTINUE
!Convert fractional coordinates to cartesian
IF(coord=='frac') THEN
  msg = 'converting to cartesian coordinates'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  CALL FRAC2CART(P,H)
  IF(ALLOCATED(S)) THEN
    CALL FRAC2CART(S,H)
  ENDIF
ENDIF
!
!Continue reading, look for velocities
IF( velocities .AND. vx>0 .AND. vy>0 .AND. vz>0 ) THEN
  msg = 'looking for velocities...'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !Atom velocities are present: read them
  !Note: velocities may not be defined for all atoms, and will
  !     remain zero if not defined in this section.
  !     They have the format  "i vx vy vz".
  DO
    READ(30,'(a128)',ERR=550,END=550) temp
    temp = ADJUSTL(temp)
    IF( temp(1:3)=="vel" ) THEN
      msg = 'Found velocities'
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      DO
        READ(30,'(a128)',ERR=550,END=550) temp
        temp = ADJUSTL(temp)
        READ(temp,*,ERR=550,END=550) i, P1, P2, P3
        IF( i>0 .AND. i<=SIZE(AUX,1) ) THEN
          AUX(i,vx) = P1
          AUX(i,vy) = P2
          AUX(i,vz) = P3
        ELSE
          !out-of-bounds
          nwarn=nwarn+1
          CALL ATOMSK_MSG(2742,(/""/),(/DBLE(i)/))
        ENDIF
      ENDDO
    ENDIF
  ENDDO
ENDIF
!
550 CONTINUE
IF( sgroupnum.NE.0 ) THEN
  !A space group number was read
  IF( sgroupnum > 0 .AND. sgroupnum <=230 ) THEN
    ! Apply symmetry operations
    CALL SG_APPLY_SYMOPS(sgroup,H,P,S,AUXNAMES,AUX)
  ELSE
    ! Invalid space group number
    nerr = nerr+1
    CALL ATOMSK_MSG(809,(/TRIM(sgroup)/),(/0.d0/))
  ENDIF
ENDIF
GOTO 1000
!
!
!
600 CONTINUE
CALL ATOMSK_MSG(1804,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i+j)/))
nerr = nerr+1
GOTO 1000
!
810 CONTINUE
CALL ATOMSK_MSG(1806,(/''/),(/DBLE(i+j)/))
nerr = nerr+1
GOTO 1000
!
820 CONTINUE
CALL ATOMSK_MSG(1809,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
830 CONTINUE
CALL ATOMSK_MSG(1810,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
1000 CONTINUE
CLOSE(30)
!
END SUBROUTINE READ_GIN
!
!
END MODULE in_gulp_gin
