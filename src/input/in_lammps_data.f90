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
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 21 April 2017                                    *
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
SUBROUTINE READ_LMP_DATA(inputfile,H,P,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=128):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(100):: tempcomment
LOGICAL:: allidzero  !are all atom id equal to zero?
LOGICAL:: dobonds    !are there bonds in the file?
LOGICAL:: flags      !are there replica flags at the end of each line?
LOGICAL:: molecule   !is there a column "moleculeID" before the column "atom-type"?
LOGICAL:: velocities !are velocities in the file?
INTEGER:: flagx, flagy, flagz !replica flags at end of lines (optional)
INTEGER:: i, id, j, k, Ncol
INTEGER:: Naux  !number of auxiliary properties
INTEGER:: NP, Nspieces !number of particles, of atomic spieces
INTEGER:: strlength
INTEGER:: molID, q, vx, vy, vz !columns for molecule ID, electric charge, velocities in AUX
INTEGER,DIMENSION(:,:),ALLOCATABLE:: BONDS   !(1) atom 1, (2) atom 2, (3) bond type
REAL(dp):: alpha, beta, gamma
REAL(dp):: a, b, c
REAL(dp):: atomtype
REAL(dp):: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz !triclinic box parameters
REAL(dp),DIMENSION(20):: column
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
tempcomment(:)=''
allidzero = .TRUE.
dobonds = .FALSE.
flags = .FALSE.
molecule = .FALSE.
velocities = .FALSE.
i = 0
Naux = 0
Nspieces = 0
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
IF(ALLOCATED(P)) DEALLOCATE(P)
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
!Determine if velocities are present in the file
DO
  READ(30,'(a128)',ERR=110,END=110) temp
  temp = ADJUSTL(temp)
  strlength = LEN_TRIM(temp)
  IF(temp(1:10)=='Velocities') THEN
    velocities = .TRUE.
    EXIT
  ENDIF
ENDDO
!
110 CONTINUE
WRITE(msg,*) 'Velocities = ', velocities
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
REWIND(30)
!
!Read header of the file
k=0
DO
  READ(30,'(a128)',ERR=810,END=810) temp
  temp = ADJUSTL(temp)
  strlength = LEN_TRIM(temp)
  IF(temp(1:1)=='#') THEN
    k=k+1
    IF(k<100) THEN
      tempcomment(k) = TRIM(ADJUSTL(temp))
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
  ELSEIF(temp(strlength-6:)=='ylo yhi') THEN
    READ(temp,*,ERR=810,END=810) ylo, yhi
  ELSEIF(temp(strlength-6:)=='zlo zhi') THEN
    READ(temp,*,ERR=810,END=810) zlo, zhi
  ELSEIF(temp(strlength-7:)=='xy xz yz') THEN
    READ(temp,*,ERR=810,END=810) xy, xz, yz
  ELSEIF(temp(1:5)=='Atoms') THEN
    READ(30,*)
    GOTO 200
  ENDIF
ENDDO
!
!
!
200 CONTINUE
!At this point, if P is not allocated then we are in trouble
IF(.NOT.ALLOCATED(P)) GOTO 820
!
!If comments were found, save them
IF(k>0) THEN
  ALLOCATE(comment(k))
  DO i=1,k
    comment(i) = tempcomment(i)
  ENDDO
ENDIF
!
!Convert supercell parameters
xhi = xhi-xlo
yhi = yhi-ylo
zhi = zhi-zlo
 a = xhi
 b = DSQRT( yhi**2 + xy**2 )
 c = DSQRT( zhi**2 + xz**2 + yz**2 )
alpha = DACOS( (xy*xz+yhi*yz)/(b*c) )
beta = DACOS(xz/c)
gamma = DACOS(xy/b)
CALL CONVMAT(a,b,c,alpha,beta,gamma,H)
!
!Read atom coordinates
WRITE(msg,*) 'Reading atom coordinates...'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
i=0
id=0
k=0
allidzero = .TRUE.
DO i=1,NP
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
  IF(i==1) THEN
    !The number of columns is unknown, determine it
    !Try to read numbers until there is an error
    DO j=1,20
      READ(temp,*,END=210,ERR=210) (column(k), k=1,j)
      Ncol = Ncol+1
    ENDDO
    210 CONTINUE
    Ncol = Ncol-2 !number of columns besides index and type
    IF( molecule ) THEN
      !Don't count the "moleculeID" column
      Ncol = Ncol-1
    ENDIF
    !
    !We know Ncol, allocate arrays
    Naux=1  !first aux. for storing atom type
    IF( velocities ) THEN
      vx=Naux+1
      vy=vx+1
      vz=vy+1
      Naux=Naux+3
    ENDIF
    IF(molecule) THEN
      Naux = Naux+1
      molID = Naux
    ENDIF
    IF( Ncol==4 ) THEN
      !Assume "atom_style charge" = id type q x y z
      !Column #3 contains charge
      Naux = Naux+1
      q=Naux
    ENDIF
    ALLOCATE( AUX(SIZE(P,1),Naux) )
    AUX(:,:) = 0.d0
    ALLOCATE( AUXNAMES(Naux) )
    IF( velocities ) THEN
      AUXNAMES(vx) = 'vx'
      AUXNAMES(vy) = 'vy'
      AUXNAMES(vz) = 'vz'
    ENDIF
    IF(q>0) THEN
      AUXNAMES(q) = 'q'
    ENDIF
    AUXNAMES(1) = 'type'
    !
    !Now read the coordinates for that atom
    IF( molecule ) THEN
      READ(temp,*,ERR=800,END=800) id, molID, atomtype, (column(j), j=1,Ncol)
    ELSE
      READ(temp,*,ERR=800,END=800) id, atomtype, (column(j), j=1,Ncol)
    ENDIF
    !
    IF( id==0 ) THEN
      id = 1
    ELSE
      !All atom id are not zero
      allidzero = .FALSE.
    ENDIF
    !
    IF( Ncol>=6 ) THEN
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
    !The number of columns Ncol is known => read coordinates
    !"atomtype" is the second column for styles "atomic", "charge"
    !which are the only ones we consider here
    IF( molecule ) THEN
      READ(temp,*,ERR=800,END=800) k, molID, atomtype, (column(j), j=1,Ncol)
    ELSE
      READ(temp,*,ERR=800,END=800) k, atomtype, (column(j), j=1,Ncol)
    ENDIF
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
  !
  IF( id<=0 .OR. id>SIZE(P,1) ) THEN
    !Atom out of bounds, cannot write out of array P, skip this atom
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2742,(/""/),(/DBLE(id)/))
    !
  ELSE
    !Set coordinates and spieces for this particle:
    !Coordinates are always the last three columns
    !except for styles "dipole" and "hybrid" which we ignore here
    IF( flags ) THEN
      !Replica flags are in the last 3 columns
      !Atom positions are in the 3 columns before that
      !Use the flags to shift atoms by box vectors appropriately
      P(id,1) = column(Ncol-5) + DBLE(flagx)*H(1,1) + DBLE(flagy)*H(2,1) + DBLE(flagz)*H(3,1)
      P(id,2) = column(Ncol-4) + DBLE(flagx)*H(1,2) + DBLE(flagy)*H(2,2) + DBLE(flagz)*H(3,2)
      P(id,3) = column(Ncol-3) + DBLE(flagx)*H(1,3) + DBLE(flagy)*H(2,3) + DBLE(flagz)*H(3,3)
    ELSE
      !Atom positions are in the last 3 columns
      P(id,1) = column(Ncol-2)
      P(id,2) = column(Ncol-1)
      P(id,3) = column(Ncol)
    ENDIF
    P(id,4) = DBLE(atomtype)
    AUX(id,1) = DBLE(atomtype)
    IF(molecule) THEN
      AUX(id,molID) = molID
    ENDIF
    IF(q>0) THEN
      AUX(id,q) = column(1)
    ENDIF
  ENDIF
  !
ENDDO
!
!
!
300 CONTINUE
!Try to read bonds and velocities
!(in this file format the section 'Velocities'
! always comes after the 'Atoms' section)
DO
  READ(30,'(a128)',ERR=500,END=500) temp
  temp = ADJUSTL(temp)
  strlength = LEN_TRIM(temp)
  IF(temp(1:10)=='Velocities') THEN
    WRITE(msg,*) 'Reading Velocities...'
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
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
  ELSEIF(temp(1:5)=='Bonds') THEN
    IF( .NOT.ALLOCATED(BONDS) ) THEN
    
    ELSE
      WRITE(msg,*) 'Reading Bonds...'
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      READ(30,'(a128)',ERR=800,END=800) temp !first line must be empty
      DO
        !In LAMMPS each line has format "id type atom1 atom2"
        !Atomsk saves it with th format "atom1 atom2 type" in array BONDS
        READ(30,*,ERR=500,END=500) id, (column(j), j=1,3)
        IF( id>0 .AND. id<=SIZE(BONDS) ) THEN
          BONDS(id,1) = NINT(column(2))
          BONDS(id,2) = NINT(column(3))
          BONDS(id,3) = NINT(column(1))
        ELSE
          !id is out-of-bounds
        ENDIF
      ENDDO
    ENDIF
    !
  ENDIF
  !
ENDDO
!
!
!
500 CONTINUE
CLOSE(30)
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
!
END SUBROUTINE READ_LMP_DATA
!
END MODULE in_lmp_data