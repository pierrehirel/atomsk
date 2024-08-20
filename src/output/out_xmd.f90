MODULE out_xmd
!
!**********************************************************************************
!*  OUT_XMD                                                                       *
!**********************************************************************************
!* This module reads in an arrays containing atomic positions, and                *
!* writes a file in the XMD format.                                               *
!* This file format is described for instance on that page:                       *
!*   http://xmd.sourceforge.net/doc/manual/xmd.html                               *
!**********************************************************************************
!* (C) Nov. 2013 - Pierre Hirel                                                   *
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
SUBROUTINE WRITE_XMD(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: isreduced, velocities
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
msg = 'entering WRITE_XMD'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Determine if atom positions are in reduced coordinates
CALL FIND_IF_REDUCED(H,P,isreduced)
WRITE(msg,*) 'isreduced:', isreduced
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Determine how many different species are present
CALL FIND_NSP(P(:,4),atypes)
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
!Write comment(s)
IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
  DO i=1,SIZE(comment)
    WRITE(ofu,'(a)') '#'//TRIM(comment(1))
  ENDDO
ENDIF
!
!
!If box is not orthogonal, warn the user
IF( DABS(H(1,2))>1.d-6 .OR. DABS(H(1,3))>1.d-6 .OR.  &
  & DABS(H(2,1))>1.d-6 .OR. DABS(H(2,3))>1.d-6 .OR.  &
  & DABS(H(3,1))>1.d-6 .OR. DABS(H(3,2))>1.d-6      ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(3709,(/""/),(/0.d0/))
ENDIF
!
!Write box size
IF( isreduced ) THEN
  !All atom positions will be scaled to box size
  WRITE(ofu,'(a10,3(f12.6,1X))') 'BOX SCALE ', H(1,1), H(2,2), H(3,3)
ELSE
  WRITE(ofu,'(a4,3(f12.6,1X))') 'BOX ', H(1,1), H(2,2), H(3,3)
ENDIF
!
!
!Write atom positions and velocities
WRITE(msg,*) SIZE(P,1)
IF( velocities ) THEN
  WRITE(ofu,'(a)') "POSVEL "//TRIM(ADJUSTL(msg))
ELSE
  WRITE(ofu,'(a)') "POSITION "//TRIM(ADJUSTL(msg))
ENDIF
!
DO i=1,SIZE(P,1)
  IF(typecol.NE.0) THEN
    !use the defined types
    Nspecies = NINT(AUX(i,typecol))
  ELSE
    !Replace atomic number by atom type
    DO j=1,SIZE(atypes,1)
      IF( atypes(j,1)==INT(P(i,4)) ) Nspecies = j
    ENDDO
  ENDIF
  !
  IF(velocities) THEN
    WRITE(temp,150) Nspecies, P(i,1), P(i,2), P(i,3), &
                  & AUX(i,vx), AUX(i,vy), AUX(i,vz)
  ELSE
    WRITE(temp,150) Nspecies, P(i,1), P(i,2), P(i,3)
  ENDIF
  !Write line to file
  WRITE(ofu,'(a)') TRIM(ADJUSTL(temp))
ENDDO
150 FORMAT(i2,1X,6(f14.8,1X))
!
!
!Write atom masses and species
DO i=1,SIZE(atypes,1)
  WRITE(ofu,'(a12,i2)') "SELECT TYPE ", i
  CALL ATOMSPECIES(atypes(i,1),species)
  IF(masscol>0) THEN
    !Atom mass is defined as an auxiliary property
    !Find an atom of this type and use its mass
    smass=0.d0
    j=0
    DO WHILE(smass<=0.d0)
      j=j+1
      IF( NINT(P(j,4)) == NINT(atypes(i,1)) ) THEN
        smass = AUX(j,masscol)
      ENDIF
    ENDDO
  ELSE
    !Compute atom mass
    CALL ATOMMASS(species,smass)
  ENDIF
  WRITE(ofu,'(a5,f12.3)') "MASS ", smass
  WRITE(ofu,'(a9,i2,1X,a2)') "TYPENAME ", i, species
ENDDO
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
END SUBROUTINE WRITE_XMD
!
END MODULE out_xmd
