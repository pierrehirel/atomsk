MODULE in_dlp_cfg
!
!**********************************************************************************
!*  IN_DLP_CFG                                                                    *
!**********************************************************************************
!* This module reads DL_POLY configuration files.                                 *
!* This format is described in the DL_POLY manual:                                *
!*   http://www.cse.scitech.ac.uk/ccg/software/DL_POLY/MANUALS/USRMAN4.01.pdf     *
!**********************************************************************************
!* (C) May 2011 - Pierre Hirel                                                    *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
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
SUBROUTINE READ_DLP_CFG(inputfile,H,P,S,comment,AUXNAMES,AUX)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: temp
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
INTEGER:: i, j
INTEGER:: levcfg, imcon, megatm
INTEGER:: Natoms, Nshells
REAL(dp):: tempreal
REAL(dp), DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, S !Atom and Shells positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: tempP, tempS !temporary positions
REAL(dp), DIMENSION(:,:),ALLOCATABLE:: AUX, tempAUX !auxiliary properties
!
!
!Initialize variables
i = 0
Natoms = 0
Nshells = 0
H(:,:) = 0.d0
IF(ALLOCATED(tempP)) DEALLOCATE(tempP)
IF(ALLOCATED(tempS)) DEALLOCATE(tempS)
IF(ALLOCATED(tempAUX)) DEALLOCATE(tempAUX)
!
!
!
100 CONTINUE
msg = 'entering READ_DLP_CFG'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
!1st line = comment
ALLOCATE(comment(1))
READ(30,'(a128)',END=810,ERR=810) comment(1)
!
!2nd line = levcfg, imcon, megatm
!levcfg: 0=positions only; 1= pos.+velocities; 2=pos.+vel.+forces
!imcon: 0=no periodicity; 1=cubic BC; 2=orthorhombic BC; 3= parallelepipedic BC; 
!megatm = total number of particles
READ(30,*,END=820,ERR=820) levcfg, imcon, megatm
IF(levcfg==1) THEN
  !Velocities follow atom positions
  ALLOCATE(AUXNAMES(3))
  AUXNAMES(1) = 'vx'
  AUXNAMES(2) = 'vy'
  AUXNAMES(3) = 'vz'
  ALLOCATE( tempAUX( megatm, 3 ) )
  tempAUX(:,:) = 0.d0
  !
ELSEIF(levcfg==2) THEN
  !Velocities + forces follow atom positions
  ALLOCATE(AUXNAMES(6))
  AUXNAMES(1) = 'vx'
  AUXNAMES(2) = 'vy'
  AUXNAMES(3) = 'vz'
  AUXNAMES(4) = 'fx'
  AUXNAMES(5) = 'fy'
  AUXNAMES(6) = 'fz'
  ALLOCATE( tempAUX( megatm, 6 ) )
  tempAUX(:,:) = 0.d0
  !
ELSEIF(levcfg>2) THEN
  nerr=nerr+1
  CALL ATOMSK_MSG(1808,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
IF(megatm<=0) THEN
  nerr=nerr+1
  CALL ATOMSK_MSG(804,(/''/),(/0.d0/))
  GOTO 1000
ELSE
  !Allocate both tempP and tempS to megatm
  !They should be both smaller that megatm,
  !the size will be adjusted later
  ALLOCATE( tempP(megatm,4) )
  tempP(:,:) = 0.d0
  ALLOCATE( tempS(megatm,4) )
  tempS(:,:) = 0.d0
ENDIF
!
!In case of HISTORY file, an additional line starting with "timestep" lies here
READ(30,'(a128)',END=810,ERR=810) temp
temp = ADJUSTL(temp)
IF(temp(1:8).NE."timestep") THEN
  !Otherwise just go back one line
  BACKSPACE(30)
ENDIF
!
!Supercell parameters
DO i=1,3
  READ(30,*,END=810,ERR=810) (H(i,j), j=1,3)
ENDDO
!
!Count atoms and shells
i = 0
Natoms=0
Nshells=0
DO i=1,megatm
  !First entry is not necessarily the atom species but "atom name"
  !meaning basically that it could be anything...
  !Here we assume that meaningful names are used, i.e. that
  !the atom species is in the first or two first letters.
  !We also assume that the shells are specified by the suffix "_s",
  !and that shells are written in the same order as atom cores
  !(either alternating core, shell, core, shell... or all
  !core positions and then all shells positions)
  READ(30,'(a72)',END=800,ERR=800) temp
  !
  !Read atom species
  temp=ADJUSTL(temp)
  species=temp(1:2)
  CALL ATOMNUMBER(species,tempreal)
  IF( tempreal==0.d0 ) THEN
    !If we didn't get a proper atom species, try
    !with the first letter only
    species = species(1:1)
    CALL ATOMNUMBER(species,tempreal)
    !If that succeeded, we're good.
    !Otherwise if it still didn't work it means that
    !we cannot understand the "atom name".
    !The output will be garbage but I don't care
  ENDIF
  !
  !Detect if it is atom or shell
  IF( INDEX(temp,"_s")==0 ) THEN
    !it's an atom
    Natoms=Natoms+1
    !Read atom position
    READ(30,*,END=800,ERR=800) (tempP(Natoms,j),j=1,3)
    !Set atom species
    tempP(Natoms,4) = tempreal
    !Read atom auxiliary properties
    IF(levcfg>=1) THEN
      !Line containing velocities
      READ(30,*,END=800,ERR=800) (tempAUX(Natoms,j),j=1,3)
    ENDIF
    IF(levcfg==2) THEN
      !Line containing forces
      READ(30,*,END=800,ERR=800) (tempAUX(Natoms,j),j=4,6)
    ENDIF
    ! 
  ELSE
    !it's a shell
    Nshells=Nshells+1
    !Check that index is not zero
    IF(Natoms==0) Natoms=Natoms+1
    !Read shell position
    READ(30,*,END=800,ERR=800) (tempS(Natoms,j),j=1,3)
    !Set shell species
    tempS(Natoms,4) = tempreal
    !Pretend to read auxiliary properties but don't save them
    IF(levcfg>=1) THEN
      !Line containing velocities
      READ(30,*,END=800,ERR=800) temp
    ENDIF
    IF(levcfg==2) THEN
      !Line containing forces
      READ(30,*,END=800,ERR=800) temp
    ENDIF
  ENDIF
  !
ENDDO
!
!
!Allocate arrays for atoms (and shells if any)
ALLOCATE(P(Natoms,4))
P(:,:) = 0.d0
IF(Nshells.NE.0) THEN
  !Note: array for shells has same size as P
  ALLOCATE(S(Natoms,4))
  S(:,:) = 0.d0
ENDIF
!
!
!Copy atoms positions to the final P (and shells to S if any)
i = 0
DO i=1,Natoms
  P(i,:) = tempP(i,:)
  IF(Nshells.NE.0) THEN
    IF( DABS(tempS(i,4))>0.1d0 ) THEN
      S(i,:) = tempS(i,:)
    ENDIF
  ENDIF
ENDDO
!
!Same with auxiliary properties if any
IF(ALLOCATED(tempAUX)) THEN
  ALLOCATE( AUX( Natoms,SIZE(tempAUX(1,:)) ) )
  DO i=1,Natoms
    AUX(i,:) = tempAUX(i,:)
  ENDDO
ENDIF
!
!
!Deallocate temporary arrays
IF(ALLOCATED(tempP)) DEALLOCATE(tempP)
IF(ALLOCATED(tempS)) DEALLOCATE(tempS)
IF(ALLOCATED(tempAUX)) DEALLOCATE(tempAUX)
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
CALL ATOMSK_MSG(1809,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
820 CONTINUE
CALL ATOMSK_MSG(1810,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
1000 CONTINUE
CLOSE(30)
!
!
END SUBROUTINE READ_DLP_CFG
!
!
END MODULE in_dlp_cfg
