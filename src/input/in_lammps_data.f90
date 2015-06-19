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
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 17 June 2014                                     *
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
LOGICAL:: velocities !are velocities in the file?
INTEGER:: i, id, j, k, Ncol
INTEGER:: Naux  !number of auxiliary properties
INTEGER:: NP, Nspieces !number of particles, of atomic spieces
INTEGER:: strlength
INTEGER:: q, vx, vy, vz !columns for electric charge, velocities
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
velocities = .FALSE.
i = 0
q = 0
Naux = 0
Nspieces = 0
strlength = 0
Ncol = 0
vx = 0
vy = 0
vz = 0
 column(:) = 0.d0
xy = 0.d0
xz = 0.d0
yz = 0.d0
H(:,:) = 0.d0
IF(ALLOCATED(P)) DEALLOCATE(P)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
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
  ELSEIF(temp(strlength-6:)=='xlo xhi') THEN
    READ(temp,*,ERR=810,END=810) xlo, xhi
  ELSEIF(temp(strlength-6:)=='ylo yhi') THEN
    READ(temp,*,ERR=810,END=810) ylo, yhi
  ELSEIF(temp(strlength-6:)=='zlo zhi') THEN
    READ(temp,*,ERR=810,END=810) zlo, zhi
  ELSEIF(temp(strlength-7:)=='xy xz yz') THEN
    READ(temp,*,ERR=810,END=810) xy, xz, yz
  ELSEIF(temp(strlength-4:)=='Atoms') THEN
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
DO i=1,NP
  READ(30,'(a128)',ERR=800,END=800) temp
  temp = ADJUSTL(temp)
  !Remove trailing comment if any
  strlength = SCAN(temp,'#')
  IF( strlength>0 ) THEN
    temp(strlength:) = ' '
  ENDIF
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
    !
    !We know Ncol, allocate arrays
    Naux=1  !first aux. for storing atom type
    IF( velocities ) THEN
      vx=Naux+1
      vy=vx+1
      vz=vy+1
      Naux=Naux+3
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
    READ(temp,*,ERR=800,END=800) id, atomtype, (column(j), j=1,Ncol)
  ELSE
    !The number of columns Ncol is known => read coordinates
    !"atomtype" is the second column for styles "atomic", "charge"
    !which are the only ones we consider here
    READ(temp,*,ERR=800,END=800) id, atomtype, (column(j), j=1,Ncol)
  ENDIF
  !
  !Check that the particle id is within the bounds
  IF(id>NP .OR. id<=0) GOTO 800
  !
  !Set coordinates and spieces for this particle:
  !Coordinates are always the last three columns
  !except for styles "dipole" and "hybrid" which we ignore here
  P(id,1) = column(Ncol-2)
  P(id,2) = column(Ncol-1)
  P(id,3) = column(Ncol)
  P(id,4) = DBLE(atomtype)
  AUX(id,1) = DBLE(atomtype)
  IF(q>0) THEN
    AUX(id,q) = column(1)
  ENDIF
ENDDO
!
!
!
300 CONTINUE
!Try to read velocities
!(in this file format the section 'Velocities'
! always comes after the 'Atoms' section)
WRITE(msg,*) 'Reading Velocities...'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
DO
  READ(30,'(a128)',ERR=500,END=500) temp
  temp = ADJUSTL(temp)
  strlength = LEN_TRIM(temp)
  IF(temp(1:10)=='Velocities') THEN
    READ(30,'(a128)',ERR=800,END=800) temp !first line must be empty
    DO
      !Each line has format  "id vx vy vz"
      !Note: velocities may not be defined for all atoms,
      !i.e. "id" may not span from 1 up to NP. In addition
      !the id of atoms may be in arbitrary order
      READ(30,*,ERR=500,END=500) id, (column(j), j=1,3)
      AUX(id,vx) = column(1)
      AUX(id,vy) = column(2)
      AUX(id,vz) = column(3)
    ENDDO
    !We are done
    GOTO 500
  ENDIF
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
