MODULE in_lmp_data
!
!
!**********************************************************************************
!*  IN_LMP_DATA                                                                   *
!**********************************************************************************
!* This module reads a data file for LAMMPS (lmp format).                         *
!* WARNING: it assumes the X, Y and Z coordinates are the three last numbers of   *
!* each line, so this will not work for atom_style = dipole, ellipsoid, hybrid.   *
!**********************************************************************************
!* (C) November 2010 - Pierre Hirel                                               *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 16 Oct. 2018                                     *
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
USE functions
USE messages
USE files
USE resize
USE subroutines
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE READ_LMP_DATA(inputfile,H,P,S,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=16) datatype  !type of data file: "atom" or "charge" or "molecule" or...
CHARACTER(LEN=128):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(100):: tempcomment
LOGICAL:: allidzero  !are all atom id equal to zero?
LOGICAL:: dobonds    !are there bonds in the file?
LOGICAL:: doneAtoms  !are atom positions read already?
LOGICAL:: shells     !are some particles shells? (in the sense of ionic core-shell model)
LOGICAL:: flags      !are there replica flags at the end of each line?
LOGICAL:: molecule   !is there a column "moleculeID" before the column "atom-type"?
LOGICAL:: velocities !are velocities in the file?
INTEGER:: flagx, flagy, flagz !replica flags at end of lines (optional)
INTEGER:: i, id, j, k, Ncol
INTEGER:: Naux     !number of auxiliary properties
INTEGER:: Ncomment !number of comment lines
INTEGER:: NP, Nspieces !number of particles, of atomic spieces
INTEGER:: Ncores, Nshells    !number of cores, shells (core-shell model)
INTEGER:: strlength
INTEGER:: molID, q, vx, vy, vz !columns for molecule ID, electric charge, velocities in AUX
INTEGER,DIMENSION(:,:),ALLOCATABLE:: Masses   !Relation between atom type, mass, core/shell
INTEGER,DIMENSION(:,:),ALLOCATABLE:: BONDS   !(1) atom 1, (2) atom 2, (3) bond type
REAL(dp):: alpha, beta, gamma
REAL(dp):: a, b, c
REAL(dp):: atomtype
REAL(dp):: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz !triclinic box parameters
REAL(dp),DIMENSION(3):: vector  !temporary vector
REAL(dp),DIMENSION(20):: column
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Ptemp  !positions of particles (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P   !positions of atoms
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S   !positions of shells (if any)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
tempcomment(:)=''
datatype = ''
allidzero = .TRUE.
dobonds = .FALSE.
doneAtoms = .FALSE.
flags = .FALSE.
molecule = .FALSE.
shells = .FALSE.
atomtype = 1
i = 0
Naux = 0
Ncomment = 0
Nspieces = 0
Ncores = 0
Nshells = 0
strlength = 0
Ncol = 0
molID = 0
q = 0
vx = 0
vy = 0
vz = 0
 column(:) = 0.d0
xlo = 0.d0
xhi = 0.d0
ylo = 0.d0
yhi = 0.d0
zlo = 0.d0
zhi = 0.d0
xy = 0.d0
xz = 0.d0
yz = 0.d0
H(:,:) = 0.d0
IF(ALLOCATED(Masses)) DEALLOCATE(Masses)
IF(ALLOCATED(Ptemp)) DEALLOCATE(Ptemp)
IF(ALLOCATED(P)) DEALLOCATE(P)
IF(ALLOCATED(S)) DEALLOCATE(S)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
IF(ALLOCATED(BONDS)) DEALLOCATE(BONDS)
!
!
100 CONTINUE
msg = 'entering READ_LMP_DATA'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
OPEN(UNIT=30,FILE=inputfile,STATUS='UNKNOWN',ERR=500)
REWIND(30)
!
!Parse the file, detect keywords
DO
  READ(30,'(a128)',ERR=500,END=500) temp
  temp = ADJUSTL(temp)
  strlength = LEN_TRIM(temp)
  IF(temp(1:1)=='#') THEN
    Ncomment = Ncomment+1
    IF(Ncomment<=SIZE(tempcomment)) THEN
      tempcomment(Ncomment) = TRIM(ADJUSTL(temp))
    ENDIF
  ELSEIF(temp(strlength-4:)=='atoms') THEN
    READ(temp,*,ERR=820,END=820) NP
    ALLOCATE(P(NP,4))
  ELSEIF(temp(strlength-9:)=='atom types') THEN
    READ(temp,*,ERR=820,END=820) Nspieces
  ELSEIF(temp(strlength-4:)=='bonds') THEN
    dobonds = .TRUE.
    molecule = .TRUE.
    READ(temp,*,ERR=820,END=820) i
    ALLOCATE(BONDS(i,3))
  ELSEIF(temp(strlength-9:)=='bond types') THEN
    dobonds = .TRUE.
    molecule = .TRUE.
  ELSEIF(temp(strlength-6:)=='xlo xhi') THEN
    READ(temp,*,ERR=810,END=810) xlo, xhi
    xhi = xhi-xlo
  ELSEIF(temp(strlength-6:)=='ylo yhi') THEN
    READ(temp,*,ERR=810,END=810) ylo, yhi
    yhi = yhi-ylo
  ELSEIF(temp(strlength-6:)=='zlo zhi') THEN
    READ(temp,*,ERR=810,END=810) zlo, zhi
    zhi = zhi-zlo
  ELSEIF(temp(strlength-7:)=='xy xz yz') THEN
    READ(temp,*,ERR=810,END=810) xy, xz, yz
    !
  ELSEIF(temp(1:6)=='Masses') THEN
    !The particles masses follow
    !=> read them, and store them as atomic numbers
    !The array stores [ type | atomic number | 0=atom or core; 1=shell ]
    !There must be one blank line after the keyword "Masses"
    READ(30,*)
    ALLOCATE(Masses(100,3))  !allow for maximum 100 different atom types
    Masses(:,:) = 0
    k=0
    Ncores = 0
    Nshells = 0
    DO WHILE( k<SIZE(Masses,1) )
      READ(30,'(a128)',ERR=120,END=120) temp
      IF( LEN_TRIM(temp)>0 ) THEN
        READ(temp,*,ERR=120,END=120) i, a
        k=k+1
        Masses(k,1) = i       !Atom type is saved in Masses(:,1)
        CALL ATOMMASSSPECIES(a,species)
        CALL ATOMNUMBER(species,b)
        Masses(k,2) = NINT(b) !Atomic number is saved in Masses(:,2)
        Masses(k,3) = 0       !by default the particle is an atom
        IF( INDEX(temp,'shell')>0 .OR. INDEX(temp,'Shell')>0 .OR. INDEX(temp,'SHELL')>0 ) THEN
          !This type of particle is a shell => mark it as such in Masses(:,3)
          Masses(k,3) = 1
          shells = .TRUE.
          Nshells = Nshells+1
          IF( k>1 ) THEN
            !Try to detect to which core it is associated
            !Parse the previous Masses
            Masses(k,2) = 0
            DO i=1,k-1
              CALL ATOMSPECIES(DBLE(Masses(i,2)),species)
              CALL ATOMMASS(species,c)
              IF( DABS( c - 10.d0*a )<=0.9d0 ) THEN
                Masses(k,2) = Masses(i,2)
              ENDIF
            ENDDO
          ENDIF
          !
          IF( Masses(k,2)==0 ) THEN
            !Failed to find a core with 10*mass of shell
            !=> assume that cores and shells are listed in same order
            !  and use atomic number of N-th core
            Ncores = 0
            DO i=1,k-1
              IF( Masses(k,3)==0 ) THEN
                Ncores = Ncores+1
                IF( Ncores==Nshells ) THEN
                  Masses(k,2) = Masses(i,2)
                  EXIT
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ELSE
          Ncores = Ncores+1
        ENDIF
      ENDIF
    ENDDO
    120 CONTINUE
    BACKSPACE(30)
    IF( verbosity>=4 ) THEN
      WRITE(msg,*) "Masses: Type    Atomic number    0=core/1=shell"
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      DO i=1,k
        WRITE(msg,*) Masses(i,:)
        CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      ENDDO
    ENDIF
    IF( doneAtoms ) THEN
      !Atom positions were read before the Masses were read
      !=> parse them again and save the correct atomic number in P(:,4)
      !NOTE: array AUX should be allocated and contain atom type in its first column
      IF( ALLOCATED(P) .AND. SIZE(P,1)>0 .AND. ALLOCATED(AUX) ) THEN
        DO i=1,SIZE(P,1)
          DO j=1,SIZE(Masses,1)
            IF( NINT(AUX(i,1)) == Masses(j,1) ) THEN
              P(i,4) = Masses(j,2)
              EXIT
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF
    !
  ELSEIF(temp(1:5)=='Atoms') THEN
    !Atom positions follow
    WRITE(msg,*) 'Reading Atoms...'
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !The type of data file may follow, e.g. "Atoms  # charge"
    IF( LEN_TRIM(temp(6:))>0 ) THEN
      strlength = SCAN(temp,"#")
      IF( strlength>=6 ) THEN
        temp = TRIM(ADJUSTL(temp(strlength+1:)))
        SELECT CASE(temp)
        CASE("angle","atomic","body","bond","charge","dipole","dpd","edpd","mdpd","tdpd",   &
            & "electron","ellipsoid","full","line","meso","molecular","peri","smd","sphere", &
            & "template","tri","wavepacket","hybrid")
          datatype = temp(1:16)
        CASE DEFAULT
          !Unknown data type or garbage => ignore it and proceed
          datatype = ""
        END SELECT
      ENDIF
    ENDIF
    msg = "Data type read from Atoms line: "//TRIM(datatype)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
    !There must be one blank line after the keyword "Atoms"
    READ(30,*,ERR=800,END=800)
    !
    !Convert supercell parameters
    a = xhi
    b = DSQRT( yhi**2 + xy**2 )
    c = DSQRT( zhi**2 + xz**2 + yz**2 )
    alpha = DACOS( (xy*xz+yhi*yz)/(b*c) )
    beta = DACOS(xz/c)
    gamma = DACOS(xy/b)
    CALL CONVMAT(a,b,c,alpha,beta,gamma,H)
    !
    !Set number of columns that follow atom-ID according to data type
    !Also set the number of auxiliary properties to allocate AUX and AUXNAMES
    Ncol = 0
    Naux = 1 !always at least 1 auxiliary property: the atom type
    SELECT CASE(datatype)
    CASE("atomic","atom")
      Ncol = 4
    CASE("angle","bond","dpd","mdpd","molecular")
      Ncol = 5
    CASE("charge")
      Ncol = 5
      Naux = Naux+1  !charge
    CASE("body","edpd","ellipsoid","full","peri","sphere")
      Ncol = 6
    CASE("electron","line","meso","template","tri")
      Ncol = 7
    CASE("dipole")
      Ncol = 8
      Naux = Naux+1  !charge
      q = Naux
    CASE("smd")
      Ncol = 9
    CASE("wavepacket")
      Ncol = 10
    CASE DEFAULT
      !Unknown data type of garbage => ignore it and proceed
      Ncol = 0
    END SELECT
    WRITE(msg,*) "Number of columns in Atoms section: ", Ncol
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
    !At this point, if P is not allocated then we are in trouble
    IF(.NOT.ALLOCATED(P)) GOTO 820
    !
     IF( shells ) THEN
       !Some particles are ionic shells => allocate memory
       ALLOCATE( S(SIZE(P,1),4) )
       S(:,:) = 0.d0
     ENDIF
    !
    !Allocate array AUX to store auxiliary properties
    ALLOCATE( AUX(SIZE(P,1),Naux) )
    AUX(:,:) = 0.d0
    ALLOCATE(AUXNAMES(Naux))
    AUXNAMES(1) = "type"
    IF( Naux>=2 ) THEN
      q = 2
      AUXNAMES(2) = "q"
    ENDIF
    !
    !Read atom coordinates
    WRITE(msg,*) 'Reading atom coordinates...'
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,*) 'SIZE P | AUX: ', SIZE(P,1), " | ", SIZE(AUX,1), SIZE(AUX,2)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    i=0
    id=0
    k=0
    allidzero = .TRUE.
    Ncores = 0   !counter for ionic cores
    Nshells = 0 !counter for ionic shells
    DO i=1,SIZE(P,1)
      shells = .FALSE.  !by default particle #i is not a shell
      READ(30,'(a128)',ERR=800,END=800) temp
      temp = ADJUSTL(temp)
      !Remove trailing comment if any
      strlength = SCAN(temp,'#')
      IF( strlength>0 ) THEN
        temp(strlength:) = ' '
      ENDIF
      flagx=0
      flagy=0
      flagz=0
      !
      IF( i==1 .AND. Ncol==0 .AND. LEN_TRIM(datatype)==0 ) THEN
        !The number of columns is unknown, determine it
        WRITE(msg,*) 'Trying to determine number of columns...'
        CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
        !Try to read numbers until there is an error
        DO j=1,20
          READ(temp,*,END=210,ERR=210) (column(k), k=1,j)
          Ncol = Ncol+1
        ENDDO
        210 CONTINUE
        WRITE(msg,*) 'Number of columns counted from 1st line of data: ', Ncol
        CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
        Ncol = Ncol-1 !first column is always atom-ID => remove it from count
        IF( molecule ) THEN
          !Don't count the "moleculeID" column
          Ncol = Ncol-1
        ENDIF
        !
        !We know Ncol, allocate arrays
        Naux=1  !first aux. for storing atom type
        IF(molecule) THEN
          Naux = Naux+1
          molID = Naux
        ENDIF
        IF( Ncol==4 ) THEN
          !Only one possibility: data type is "atomic"
          datatype = "atomic"
        ELSEIF( Ncol==5 ) THEN
          !Assume "atom_style charge" = atom-ID atom-type q x y z
          !Column #3 contains charge
          datatype = "charge"
          Naux = Naux+1
          q = Naux
        ELSEIF( Ncol==6 ) THEN
          !Assume "atom_style full" = atom-ID molecule-ID atom-type q x y z
          !Column #3 contains charge
          datatype = "full"
        ENDIF
        !
        IF( .NOT. ALLOCATED(AUX) ) THEN
          ALLOCATE( AUX(SIZE(P,1),Naux) )
          AUX(:,:) = 0.d0
          ALLOCATE( AUXNAMES(Naux) )
          AUXNAMES(1) = 'type'
          IF(q>0) THEN
            AUXNAMES(2) = 'q'
          ENDIF
        ENDIF
        !
        !Now read the coordinates for that atom
        READ(temp,*,ERR=800,END=800) id, (column(j), j=1,Ncol)
        !
        IF( id==0 ) THEN
          id = 1
        ELSE
          !All atom id are not zero
          allidzero = .FALSE.
        ENDIF
        !
        IF( Ncol>=6 .AND. datatype.NE."full" ) THEN
          !If there are at least 6 columns, determine if
          !the last three columns contain replica flags (check for 3 integers)
          flagx = NINT(column(Ncol-2))
          flagy = NINT(column(Ncol-1))
          flagz = NINT(column(Ncol))
          IF( DABS(column(Ncol-2) - DBLE(flagx)) < 1.d-12 .AND.      &
            & DABS(column(Ncol-1) - DBLE(flagy)) < 1.d-12 .AND.      &
            & DABS(column(Ncol)   - DBLE(flagz)) < 1.d-12       ) THEN
            !The 3 last columns are indeed 3 integers
            flags = .TRUE.
            WRITE(msg,*) 'Detected replica flags on first line: ', flagx, flagy, flagz
            CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
          ENDIF
        ENDIF
        220 CONTINUE
        !
      ELSE
        !The number of columns Ncol is known => read values
        !First value must always be the atom-ID
        READ(temp,*,ERR=800,END=800) k, (column(j), j=1,Ncol)
        !
        !Check that the particle id (here read as "k") is within the bounds
        IF( k==0 ) THEN
          IF( allidzero ) THEN
            !So far, all atoms encountered have an id equal to zero
            !Just increment the id automatically
            id = id+1
          ELSE
            !All atom id are not zero (at least one previous atom had a non-zero id)
            !=> this particular atom is out of bounds (warning will be displayed later)
            id = k
          ENDIF
        ELSE
          !k is not zero => this atom id is not zero
          id = k
          allidzero = .FALSE.
        ENDIF
        !
        !Read values of flags
        IF( flags ) THEN
          flagx = NINT(column(Ncol-2))
          flagy = NINT(column(Ncol-1))
          flagz = NINT(column(Ncol))
        ENDIF
      ENDIF
      !
      IF( id<=0 .OR. id>SIZE(P,1) ) THEN
        !Atom out of bounds, cannot write out of array P, skip this atom
        nwarn=nwarn+1
        CALL ATOMSK_MSG(2742,(/""/),(/DBLE(id)/))
        !
      ELSE
        !Save atom type as an auxiliary property
        IF ( datatype=="full" ) THEN
          atomtype = column(Ncol-4)
        ELSEIF ( datatype=="angle" .OR. datatype=="bond" .OR. datatype=="molecular" .OR. datatype=="template" ) THEN
          atomtype = column(Ncol-3)
        ELSE
          atomtype = column(1)
        ENDIF
        !Save atom species
        IF( ALLOCATED(Masses) .AND. SIZE(Masses,1)>0 ) THEN
          !Look for the mass of this type of atom
          j=0
          DO WHILE( j<=SIZE(Masses,1) )
            j = j+1
            IF( Masses(j,1) == NINT(atomtype) ) THEN
              IF( Masses(j,3)>=1 ) THEN
                !This atom type is an ionic shell
                Nshells = Nshells+1
                shells = .TRUE.
                S(Nshells,4) = Masses(j,2)
                EXIT
              ELSE
                !This atom type is an ionic core
                Ncores = Ncores+1
                P(Ncores,4) = Masses(j,2)
                EXIT
              ENDIF
            ENDIF
          ENDDO
          !Make sure that the atomic number is not zero
          IF( shells ) THEN
            IF( S(Nshells,4)<1.d-12 ) THEN
              !Mass of this shell unknown => use atom type = atom species
              S(Nshells,4) = DBLE(atomtype)
            ENDIF
          ELSE
            IF( P(Ncores,4)<1.d-12 ) THEN
              !Mass of this atom unknown => use atom type = atom species
              P(Ncores,4) = DBLE(atomtype)
            ENDIF
          ENDIF
        ELSE
          !No Masses section in the file => use atom type = atom species
          P(id,4) = DBLE(atomtype)
        ENDIF
        !
        !Set coordinates for this particle
        IF ( datatype=="dipole" ) THEN
          !Style  "dipole": x y z are in columns 4 5 6
          vector(1) = column(3)
          vector(2) = column(4)
          vector(3) = column(5)
        ELSEIF ( datatype=="tdpd" .OR. datatype=="hybrid" ) THEN
          !Styles "tdpd" and "hybrid": x y z are in columns 3 4 5
          vector(1) = column(2)
          vector(2) = column(3)
          vector(3) = column(4)
        ELSE
          !All other styles: x y z are in the last 3 columns
          !(except if three integer flags are in the last 3 columns)
          IF( flags ) THEN
            !Replica flags are in the last 3 columns
            !Atom positions are in the 3 columns before that
            !Use the flags to shift atoms by box vectors appropriately
            vector(1) = column(Ncol-5) + DBLE(flagx)*H(1,1) + DBLE(flagy)*H(2,1) + DBLE(flagz)*H(3,1)
            vector(2) = column(Ncol-4) + DBLE(flagx)*H(1,2) + DBLE(flagy)*H(2,2) + DBLE(flagz)*H(3,2)
            vector(3) = column(Ncol-3) + DBLE(flagx)*H(1,3) + DBLE(flagy)*H(2,3) + DBLE(flagz)*H(3,3)
          ELSE
            !Atom positions are in the last 3 columns
            vector(1) = column(Ncol-2)
            vector(2) = column(Ncol-1)
            vector(3) = column(Ncol)
          ENDIF
        ENDIF
        !
        !Save particle position in P or in S depending on its type
        IF( shells ) THEN
          S(Nshells,1:3) = vector(:)
          shells = .FALSE. !we don't know if next particle is a shell
        ELSEIF( Ncores>0 ) THEN
          P(Ncores,1:3) = vector(:)
          AUX(Ncores,1) = DBLE(atomtype)
        ELSE
          P(id,1:3) = vector(:)
          AUX(id,1) = DBLE(atomtype)
        ENDIF
        !
!         IF(molecule) THEN
!           AUX(id,molID) = molID
!         ENDIF
!         IF(q>0) THEN
!           AUX(id,q) = column(1)
!         ENDIF
      ENDIF
      !
    ENDDO  !end loop on particles positions
    !
    !If cores/shells were detected, resize arrays P, S and AUX
    IF( Ncores>0 .AND. Nshells>0 ) THEN
      CALL RESIZE_DBLEARRAY2(P,Ncores,4,i)
      IF( i>0 ) THEN
        nerr=nerr+1
        CALL ATOMSK_MSG(818,(/"P"/),(/0.d0/))
        GOTO 1000
      ENDIF
      CALL RESIZE_DBLEARRAY2(S,Ncores,4,i)
      IF( i>0 ) THEN
        nerr=nerr+1
        CALL ATOMSK_MSG(818,(/"S"/),(/0.d0/))
        GOTO 1000
      ENDIF
      CALL RESIZE_DBLEARRAY2(AUX,Ncores,SIZE(AUX,2),i)
      IF( i>0 ) THEN
        nerr=nerr+1
        CALL ATOMSK_MSG(818,(/"AUX"/),(/0.d0/))
        GOTO 1000
      ENDIF
    ENDIF
    !
    doneAtoms = .TRUE.
    !
    !
  ELSEIF(temp(1:10)=='Velocities') THEN
    WRITE(msg,*) 'Reading Velocities...'
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    IF( vx==0 .AND. vy==0 .AND. vz==0 ) THEN
      vx = SIZE(AUXNAMES)+1
      vy = vx+1
      vz = vy+1
      CALL RESIZE_DBLEARRAY2(AUX,SIZE(P,1),SIZE(AUX,2)+3,j)
      IF( j>0 ) THEN
        nerr=nerr+1
        CALL ATOMSK_MSG(818,(/"AUX"/),(/0.d0/))
        GOTO 1000
      ENDIF
      !
      IF( ALLOCATED(AUXNAMES) ) DEALLOCATE(AUXNAMES)
      ALLOCATE( AUXNAMES(SIZE(AUX,2)) )
      AUXNAMES(1) = 'type'
      IF(q>0) THEN
        AUXNAMES(q) = 'q'
      ENDIF
      AUXNAMES(vx) = 'vx'
      AUXNAMES(vy) = 'vy'
      AUXNAMES(vz) = 'vz'
    ENDIF
    READ(30,'(a128)',ERR=800,END=800) temp !first line must be empty
    DO
      !Each line has format  "id vx vy vz"
      !Note: velocities may not be defined for all atoms,
      !i.e. "id" may not span from 1 up to NP. In addition
      !the id of atoms may be in arbitrary order
      READ(30,*,ERR=500,END=500) id, (column(j), j=1,3)
      IF( id>0 .AND. id<=SIZE(AUX,1) ) THEN
        AUX(id,vx) = column(1)
        AUX(id,vy) = column(2)
        AUX(id,vz) = column(3)
      ELSE
        !id is out-of-bounds
      ENDIF
    ENDDO
    !
!   ELSEIF(temp(1:5)=='Bonds') THEN
!     IF( .NOT.ALLOCATED(BONDS) ) THEN
!       !Skip
!     ELSE
!       WRITE(msg,*) 'Reading Bonds...'
!       CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!       READ(30,'(a128)',ERR=800,END=800) temp !first line must be empty
!       DO
!         !In LAMMPS each line has format "id type atom1 atom2"
!         !Atomsk saves it with th format "atom1 atom2 type" in array BONDS
!         READ(30,*,ERR=500,END=500) id, (column(j), j=1,3)
!         IF( id>0 .AND. id<=SIZE(BONDS) ) THEN
!           BONDS(id,1) = NINT(column(2))
!           BONDS(id,2) = NINT(column(3))
!           BONDS(id,3) = NINT(column(1))
!         ELSE
!           !id is out-of-bounds
!         ENDIF
!       ENDDO
!     ENDIF
    !
  ENDIF
  !
  350 CONTINUE
  !
ENDDO
!
!
!
500 CONTINUE
CLOSE(30)
!
!If comments were found, save them
IF(Ncomment>0) THEN
  ALLOCATE(comment(Ncomment))
  DO i=1,Ncomment
    comment(i) = tempcomment(i)
  ENDDO
ENDIF
!
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
!Unable to read box vectors
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
IF(ALLOCATED(Masses)) DEALLOCATE(Masses)
!
END SUBROUTINE READ_LMP_DATA
!
END MODULE in_lmp_data
