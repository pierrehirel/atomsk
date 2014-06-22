MODULE out_siesta_xv
!
!
!**********************************************************************************
!*  OUT_XV                                                                        *
!**********************************************************************************
!* This module writes SIESTA XV format.                                           *
!* The XV format is described in the SIESTA manual:                               *
!*    http://www.icmab.es/dmmis/leem/siesta/                                      *
!**********************************************************************************
!* (C) Jan. 2012 - Eva Marie Kalivoda                                             *
!*     Fraunhofer Institute für Werkstoffmechanik IWM                             *
!*     Wöhlerstr. 11, 79108 Freiburg im Breisgau                                  *
!* Last modification: P. Hirel - 26 March 2014                                    *
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
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE WRITE_XV(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: vel !are there velocities or not?
INTEGER:: vx, vy, vz !index for velocities in AUX
INTEGER:: i, j, Nspecies
INTEGER:: Ntypes    !number of atom types
INTEGER:: typecol !index of atom types in AUX
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: atypes
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
typecol = 0
Nspecies = 0
vx = 0
vy = 0
vz = 0
vel=.FALSE.
!Find how many different species are in P
CALL FIND_NSP(P(:,4),atypes)
Ntypes = SIZE(atypes(:,1))
!
!
msg = 'entering WRITE_XV'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Check if some auxiliary properties are present
IF( ALLOCATED(AUXNAMES) ) THEN
  DO i=1,SIZE(AUXNAMES)
    IF( TRIM(ADJUSTL(AUXNAMES(i)))=='type' ) THEN
      typecol = i
      !Find how many different species are in AUX
      CALL FIND_NSP(AUX(:,typecol),atypes)
      Ntypes = SIZE(atypes(:,1))
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='vx' ) THEN
      vx = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='vy' ) THEN
      vy = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='vz' ) THEN
      vz = i
    ENDIF
  ENDDO
ENDIF
!
IF(vx.NE.0 .AND. vy.NE.0 .AND. vz.NE.0) vel=.TRUE.
!
!
!
OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN',ERR=500)
!
!Write header of XV file
WRITE(40,105) H(1,1), H(1,2), H(1,3), zero, zero, zero
WRITE(40,105) H(2,1), H(2,2), H(2,3), zero, zero, zero
WRITE(40,105) H(3,1), H(3,2), H(3,3), zero, zero, zero
WRITE(40,'(3X,i9)') SIZE(P(:,1))
105 FORMAT(3X,3(2X,f16.9),3X,3(2X,f16.9))
!
!Write atom positions
DO i=1,SIZE(P(:,1))
  IF( typecol.NE.0 ) THEN
    !If atom types are in auxiliary properties, use it
    Nspecies = NINT( AUX(i,typecol) )
  ELSE
    !Otherwise replace species by atom types in their order of appearance
    DO j=1,SIZE(atypes(:,1))
      IF( atypes(j,1)==INT(P(i,4)) ) Nspecies = j
    ENDDO
  ENDIF
  IF(vel) THEN
    WRITE(40,110) Nspecies, INT(P(i,4)), P(i,1), P(i,2), P(i,3), &
                & AUX(i,vx), AUX(i,vy), AUX(i,vz)
  ELSE
    WRITE(40,110) Nspecies, INT(P(i,4)), P(i,1), P(i,2), P(i,3), &
                & zero, zero, zero
  ENDIF
ENDDO
110 FORMAT(i3,3X,i3,3(2X,f16.9),3X,3(2X,f16.9))
!
!
!
!
500 CONTINUE
CLOSE(40)
msg = "XV"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
1000 CONTINUE
!
END SUBROUTINE WRITE_XV
!
END MODULE out_siesta_xv
