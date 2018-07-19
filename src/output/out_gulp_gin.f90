MODULE out_gulp_gin
!
!
!**********************************************************************************
!*  OUT_GULP_GIN                                                                  *
!**********************************************************************************
!* This module writes a GULP input file (.gin) draft.                             *
!* Only atomic positions and cell parameters are written, so the rest of          *
!* the GIN file (atomic potential, simulation parameters, etc.) has to be         *
!* completed by the user.                                                         *
!* This format is described in the GULP manual, which can be found here:          *
!*    https://www.ivec.org/gulp/help/manuals.html                                 *
!**********************************************************************************
!* (C) March 2010 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
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
USE atoms
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
SUBROUTINE WRITE_GIN(H,P,comment,AUXNAMES,AUX,outputfile,S_in)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
INTEGER:: i, NP
LOGICAL:: isreduced
INTEGER:: fixx, fixy, fixz !columns in AUX where the fixes are defined
INTEGER:: q, qs  !columns in AUX where charges of cores/shells are defined
INTEGER:: occ, radius !columns in AUX where where occupancy and radius are defined
INTEGER:: vx, vy, vz !columns in AUX where the velocities are defined
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
fixx=0
fixy=0
fixz=0
vx = 0
vy = 0
vz = 0
q = 0
qs = 0
occ = 0
radius = 0
NP = SIZE(P,1)
!NP will be total number of particles (cores+shells)
!=> count how many atoms actually have shells
IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
  DO i=1, SIZE(S(:,1))
    IF( DABS(S(i,4))>0.1d0 ) NP=NP+1
  ENDDO
ENDIF
!
!
100 CONTINUE
msg = 'entering WRITE_GIN'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Check if some atoms are fixed
IF( ALLOCATED(AUX) .AND. SIZE(AUX,1).NE.0 ) THEN
  DO i=1,SIZE(AUXNAMES)
    IF(AUXNAMES(i)=="vx") vx=i
    IF(AUXNAMES(i)=="vy") vy=i
    IF(AUXNAMES(i)=="vz") vz=i
    IF(AUXNAMES(i)=="fixx") fixx=i
    IF(AUXNAMES(i)=="fixy") fixy=i
    IF(AUXNAMES(i)=="fixz") fixz=i
    IF(AUXNAMES(i)=="q") q=i             !charges for cores
    IF(AUXNAMES(i)=="qs") qs=i           !charges for shells
    IF(AUXNAMES(i)=="occupancy") occ=i   !site occupancy
    IF(AUXNAMES(i)=="bsradius") radius=i !radius for breathing shell
  ENDDO
ENDIF
!
!Check if coordinates are reduced or cartesian
CALL FIND_IF_REDUCED(P,isreduced)
!
!
!
200 CONTINUE
OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN')
!
!Write header of GIN file
WRITE(40,'(a10)') '# Keywords'
WRITE(40,'(a4)') 'opti'
WRITE(40,*) ''
WRITE(40,'(a5)') 'title'
DO i=1,SIZE(comment)
  WRITE(40,*) ' '//TRIM(ADJUSTL(comment(i)))
ENDDO
WRITE(40,'(a3)') 'end'
WRITE(40,*) ''
WRITE(40,*) '### insert simulation parameters here ###'
WRITE(40,*) ''
WRITE(40,'(a7)') 'vectors'
WRITE(40,210) '  ', H(1,1), H(1,2), H(1,3)
WRITE(40,210) '  ', H(2,1), H(2,2), H(2,3)
WRITE(40,210) '  ', H(3,1), H(3,2), H(3,3)
IF( fixx.NE.0 .OR. fixy.NE.0 .OR. fixz.NE.0 ) THEN
  WRITE(40,'(a15)') '   1 1 1  1 1 1'
ENDIF
WRITE(40,*) ''
210 FORMAT(a2,3(f16.8))
!
!Write fractional or cartesian
IF(isreduced) THEN
  WRITE(40,'(a11,1X,i8)') 'fractional ', NP
ELSE
  WRITE(40,'(a10,1X,i8)') 'cartesian ', NP
ENDIF
!
!Write positions of cores (and shells if any)
!Format of each line:  at.no. x y z charge occupancy radius fixx fixy fixz
DO i=1,SIZE(P,1)
  CALL ATOMSPECIES(P(i,4),species)
  WRITE(msg,221) species, ' core ', P(i,1), P(i,2), P(i,3)
  !
  IF( q.NE.0 ) THEN
    !If charges are defined, write them
    WRITE(temp,'(f16.8)') AUX(i,q)
    msg = TRIM(msg)//' '//TRIM(ADJUSTL(temp))
    !
    IF( occ.NE.0 ) THEN
      !If charge is specified, site occupancy may also be
      WRITE(temp,'(f4.1)') AUX(i,occ)
      msg = TRIM(msg)//' '//TRIM(ADJUSTL(temp))
      IF( radius.NE.0 ) THEN
        !If charge and occupancy are given, radius of breathing shell may also be
        !(breathing radius will be defined for the shell, not the core,
        ! so just write a zero for the core in this column)
        WRITE(temp,'(f4.1)') zero
        msg = TRIM(msg)//' '//TRIM(ADJUSTL(temp))
      ENDIF
    ENDIF
  ENDIF
  !
  IF( fixx.NE.0 .OR. fixy.NE.0 .OR. fixz.NE.0 ) THEN
    !Caution: internally if AUX(fix)==1 then atom is fixed
    !while in GULP the flag "1" means that atom is mobile
    IF( AUX(i,fixx)>0.5d0 ) THEN
      msg = TRIM(msg)//" 0"
    ELSE
      msg = TRIM(msg)//" 1"
    ENDIF
    IF( AUX(i,fixy)>0.5d0 ) THEN
      msg = TRIM(msg)//" 0"
    ELSE
      msg = TRIM(msg)//" 1"
    ENDIF
    IF( AUX(i,fixz)>0.5d0 ) THEN
      msg = TRIM(msg)//" 0"
    ELSE
      msg = TRIM(msg)//" 1"
    ENDIF
  ELSE
    msg = TRIM(msg)
  ENDIF
  WRITE(40,'(a)') TRIM(msg)
  !
  !Write positions of corresponding shell (if any)
  !Note: shells are never fixed, line always ends with "1 1 1"
  IF( ALLOCATED(S) .AND. SIZE(S,1)>0 ) THEN
    IF( DABS(S(i,4))>0.1d0 ) THEN
      CALL ATOMSPECIES(S(i,4),species)
      WRITE(msg,221) species, ' shel ', S(i,1), S(i,2), S(i,3)
      !
      IF( q.NE.0 ) THEN
        !Charges must be written
        IF( qs.NE.0 ) THEN
          WRITE(temp,'(f16.8)') AUX(i,qs)
        ELSE
          !If qs (shell charges) are not defined, just write zero
          WRITE(temp,'(f16.8)') zero
        ENDIF
        msg = TRIM(msg)//' '//TRIM(ADJUSTL(temp))
        IF( occ.NE.0 ) THEN
          !If charge is specified, site occupancy may also be
          !Note: shell has same occupancy as its core
          WRITE(temp,'(f4.1)') AUX(i,occ)
          msg = TRIM(msg)//' '//TRIM(ADJUSTL(temp))
          IF( radius.NE.0 ) THEN
            !If charge and occupancy are given,
            !radius of breathing shell may also be
            WRITE(temp,'(f4.1)') AUX(i,radius)
            msg = TRIM(msg)//' '//TRIM(ADJUSTL(temp))
          ENDIF
        ENDIF
      ENDIF
      !
      IF( fixx.NE.0 .OR. fixy.NE.0 .OR. fixz.NE.0 ) THEN
        !Write fix flags
        !Note: shells are never fixed
        msg = TRIM(msg)//" 1 1 1"
      ENDIF
      WRITE(40,'(a)') TRIM(msg)
    ENDIF
  ENDIF
ENDDO
!
221 FORMAT(a2,a6,3(f16.8,1X))
!
!
250 CONTINUE
!Write velocities if they are defined
IF( vx.NE.0 .AND. vy.NE.0 .AND. vz.NE.0 ) THEN
  WRITE(40,'(a18)') 'velocities angs/ps'
  DO i=1,SIZE(P,1)
    WRITE(40,'(i8,3(1X,f16.8))') i, AUX(i,vx), AUX(i,vy), AUX(i,vz)
  ENDDO
ENDIF
!
!
!
500 CONTINUE
CLOSE(40)
msg = "GIN"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
1000 CONTINUE
!
END SUBROUTINE WRITE_GIN
!
END MODULE out_gulp_gin
