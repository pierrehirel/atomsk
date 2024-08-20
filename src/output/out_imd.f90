MODULE out_imd
!
!**********************************************************************************
!*  OUT_IMD                                                                       *
!**********************************************************************************
!* This module reads in an arrays containing atomic positions, and                *
!* writes a file in the IMD format.                                               *
!* This file format is described for instance on that page:                       *
!*   http://www.itap.physik.uni-stuttgart.de/~imd/userguide/config.html           *
!**********************************************************************************
!* (C) Feb. 2011 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 20 Aug. 2024                                     *
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
SUBROUTINE WRITE_IMD(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: velocities
INTEGER:: i, j, Nspecies
INTEGER:: vx, vy, vz, typecol, masscol !index of velocities, charges, types, mass in AUX
REAL(dp):: smass
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: atypes
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
velocities = .FALSE.
vx = 0
vy = 0
vz = 0
masscol = 0
typecol = 0
!
msg = 'entering WRITE_IMD'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Check if velocities are in auxiliary properties
IF( ALLOCATED(AUXNAMES) ) THEN
  DO i=1,SIZE(AUXNAMES)
    IF( TRIM(ADJUSTL(AUXNAMES(i)))=='vx' ) THEN
      vx = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='vy' ) THEN
      vy = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='vz' ) THEN
      vz = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='mass' ) THEN
      masscol = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='type' ) THEN
      typecol = i
    ENDIF
  ENDDO
  !
  IF( vx.NE.0 .AND. vy.NE.0 .AND. vz.NE.0 ) velocities = .TRUE.
ENDIF
!
!
100 CONTINUE
IF(ofu.NE.6) THEN
  OPEN(UNIT=ofu,FILE=outputfile,STATUS='UNKNOWN',ERR=1000)
ENDIF
!Write header
!1st line:  #F f n t m c v d = Format of IMD file
!  f: A for ASCII, B or b for binary big endian, L or l for binary little endian
!  n: number of columns for the atom number (0 or 1)
!  t: number of columns for the atom type (0 or 1)
!  m: number of columns for the mass (0 or 1)
!  c: number of columns for the atom coordinates (2 or 3)
!  v: number of columns for the atom velocities (0, 2 or 3)
!  d: number of columns for the atom data, e.g. Epot (0, 1,...)
!2nd line: #C describes the format of each line
IF(velocities) THEN
  WRITE(ofu,'(a16)') '#F A 1 1 1 3 3 0'
  WRITE(ofu,'(a34)') '#C number type mass x y z vx vy vz'
ELSE
  WRITE(ofu,'(a16)') '#F A 1 1 1 3 0 0'
  WRITE(ofu,'(a25)') '#C number type mass x y z'
ENDIF
!Supercell parameters
WRITE(ofu,'(a3,3(f12.6,1X))') '#X ', H(1,1), H(1,2), H(1,3)
WRITE(ofu,'(a3,3(f12.6,1X))') '#Y ', H(2,1), H(2,2), H(2,3)
WRITE(ofu,'(a3,3(f12.6,1X))') '#Z ', H(3,1), H(3,2), H(3,3)
!Comment
IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
  WRITE(ofu,'(a)') '#'//TRIM(comment(1))
ENDIF
!Label indicating end of header
WRITE(ofu,'(a2)') '#E'
!
!Determine how many different species are present
CALL FIND_NSP(P(:,4),atypes)
!
!Write atom positions
DO i=1,SIZE(P(:,1))
  IF(typecol.NE.0) THEN
    !use the defined types
    !Note: in AUX the "types" start at 1, but in
    !     IMD we have to make it start at 0
    Nspecies = NINT(AUX(i,typecol))-1
  ELSE
    !Replace atomic number by atom type
    DO j=1,SIZE(atypes(:,1))
      IF( atypes(j,1)==INT(P(i,4)) ) Nspecies = j-1
    ENDDO
  ENDIF
  !
  IF(masscol>0) THEN
    !Atom mass is defined as an auxiliary property
    smass = AUX(i,masscol)
  ELSE
    !Compute atom mass
    CALL ATOMSPECIES(P(i,4),species)
    CALL ATOMMASS(species,smass)
  ENDIF
  !
  IF(velocities) THEN
    WRITE(temp,150) i, Nspecies, smass, P(i,1), P(i,2), P(i,3), &
                  & AUX(i,vx), AUX(i,vy), AUX(i,vz)
  ELSE
    WRITE(temp,150) i, Nspecies, smass, P(i,1), P(i,2), P(i,3)
  ENDIF
  !Write line to file
  WRITE(ofu,'(a)') TRIM(temp)
ENDDO
150 FORMAT(i8,1X,i2,1X,7(f14.8,1X))
!
!
!
200 CONTINUE
IF(ofu.NE.6) THEN
  CLOSE(ofu)
ENDIF
!
!
1000 CONTINUE
!
!
END SUBROUTINE WRITE_IMD
!
END MODULE out_imd
