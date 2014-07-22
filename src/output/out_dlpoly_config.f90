MODULE out_dlp_cfg
!
!**********************************************************************************
!*  OUT_DLP_CFG                                                                   *
!**********************************************************************************
!* This module reads in an arrays containing atomic positions, and                *
!* writes a file in the DL_POLY CONFIG format.                                    *
!* This file format is described in the manual of DL_POLY:                        *
!*   http://www.cse.scitech.ac.uk/ccg/software/DL_POLY/MANUALS/USRMAN4.01.pdf     *
!**********************************************************************************
!* (C) May 2011 - Pierre Hirel                                                    *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 04 April 2014                                    *
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
SUBROUTINE WRITE_DLP_CFG(H,P,comment,AUXNAMES,AUX,outputfile,S_in)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
INTEGER:: fx, fy, fz, vx, vy, vz !Columns where Forces/Velocities are stored in AUX
INTEGER:: i, j
INTEGER:: levcfg, megatm  !level of CONFIG file, total number of particles
INTEGER:: Natoms, Nshells !indices for atoms and shells
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P          !positions of ionic cores
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN),OPTIONAL:: S_in !positions of ionic shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S !positions of ionic shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
IF(ALLOCATED(S)) DEALLOCATE(S)
IF(PRESENT(S_in) .AND. ALLOCATED(S_in) .AND. SIZE(S_in,1)>0 ) THEN
  ALLOCATE( S(SIZE(S_in,1),SIZE(S_in,2)) )
  S = S_in
ENDIF
vx = 0
vy = 0
vz = 0
fx = 0
fy = 0
fz = 0
levcfg = 0
Nshells = 0
IF(ALLOCATED(S)) THEN
  DO i=1, SIZE(S(:,1))
    IF( DABS(S(i,4))>0.1d0 ) Nshells=Nshells+1
  ENDDO
ENDIF
!
!
!
100 CONTINUE
msg = 'entering WRITE_DLP_CFG'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Check if velocities or forces are present as auxiliary properties
IF(ALLOCATED(AUXNAMES)) THEN
  DO i=1,SIZE(AUXNAMES)
    IF( AUXNAMES(i)=='fx' ) THEN
      fx = i
    ELSEIF( AUXNAMES(i)=='fy' ) THEN
      fy = i
    ELSEIF( AUXNAMES(i)=='fz' ) THEN
      fz = i
    ELSEIF( AUXNAMES(i)=='vx' ) THEN
      vx = i
    ELSEIF( AUXNAMES(i)=='vy' ) THEN
      vy = i
    ELSEIF( AUXNAMES(i)=='vz' ) THEN
      vz = i
    ENDIF
  ENDDO
  !
  IF(vx>0 .AND. vy>0 .AND. vz>0) THEN
    IF(fx>0 .AND. fy>0 .AND. fz>0) THEN
      levcfg = 2
    ELSE
      levcfg = 1
    ENDIF
  ENDIF
ENDIF
!
!
!
200 CONTINUE
OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN',ERR=1000)
!1st line = comment (limited to 72 characters)
WRITE(40,'(a)') TRIM( comment(1)(1:72) )
!
!levcfg: 0=positions only; 1= pos.+velocities; 2=pos.+vel.+forces
!imcon: 0=no periodicity; 1=cubic BC; 2=orthorhombic BC; 3= parallelepipedic BC
!megatm = total number of particles
megatm = SIZE(P(:,1))+Nshells
WRITE(40,'(3i10)') levcfg, 3, megatm
!
!Supercell parameters
DO i=1,3
  WRITE(40,'(3f20.8)') (H(i,j), j=1,3)
ENDDO
!
!Write atoms and shells positions
Natoms=0
Nshells=0
DO WHILE( Natoms+Nshells < megatm )
  Natoms = Natoms+1
  CALL ATOMSPECIES(P(Natoms,4),species)
  WRITE(temp,'(a2,3X,i8)') species, Natoms+Nshells
  WRITE(40,'(a)') TRIM(ADJUSTL(temp))
  WRITE(40,'(3f20.8)') (P(Natoms,j), j=1,3)
  IF(levcfg>=1) THEN
    WRITE(40,'(3f20.8)') AUX(Natoms,vx), AUX(Natoms,vy), AUX(Natoms,vz)
  ENDIF
  IF(levcfg==2) THEN
    WRITE(40,'(3f20.8)') AUX(Natoms,fx), AUX(Natoms,fy), AUX(Natoms,fz)
  ENDIF
  !
  !Now check if this atom has a shell
  IF( ALLOCATED(S) ) THEN
    IF( DABS(S(Natoms,4))>0.1d0 ) THEN
      !This shell corresponds to that atom: write its position
      WRITE(temp,'(a2,a2,3X,i8)') TRIM(species), "_s ", Natoms+Nshells+1
      WRITE(40,'(a)') TRIM(ADJUSTL(temp))
      WRITE(40,'(3f20.8)') (S(Natoms,j), j=1,3)
      !If levcfg>0, just write zeros for velocities and forces
      IF(levcfg>=1) THEN
        WRITE(40,'(3f20.8)') zero, zero, zero
      ENDIF
      IF(levcfg==2) THEN
        WRITE(40,'(3f20.8)') zero, zero, zero
      ENDIF
      Nshells = Nshells+1
    ENDIF
  ENDIF
ENDDO
!
!
CLOSE(40)
!
msg = "DL_POLY CONFIG"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
1000 CONTINUE
!
!
END SUBROUTINE WRITE_DLP_CFG
!
END MODULE out_dlp_cfg
