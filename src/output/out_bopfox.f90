MODULE out_bopfox
!
!
!**********************************************************************************
!*  OUT_BOPFOX                                                                    *
!**********************************************************************************
!* This module writes files used by the BOPfox program, which is described in:    *
!*   T. Hammerschmidt et al., Comput.Phys.Comm. 235 (2019) 221-233                *
!* More information from the Web site:                                            *
!*   http://bopfox.de/                                                            *
!**********************************************************************************
!* (C) November 2016 - Matous Mrovec                                              *
!*     ICAMS                                                                      *
!*     Ruhr-Universitaet Bochum, Germany                                          *
!*     matous.mrovec@icams.rub.de                                                 *
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
SUBROUTINE WRITE_BOPFOX(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=2):: species
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: isreduced
INTEGER:: fixx, fixy, fixz
INTEGER:: magx, magy, magz
INTEGER:: i, NPd, NPinert
REAL(dp):: len1, len2, len3
REAL(dp):: P1, P2, P3
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
fixx=0
fixy=0
fixz=0
magx=0
magy=0
magz=0
NPd = 0
NPinert = 0
len1 = 1.d0
len2 = 1.d0
len3 = 1.d0
!
!
100 CONTINUE
IF(ofu.NE.6) THEN
  OPEN(UNIT=ofu,FILE=outputfile,FORM='FORMATTED',STATUS='UNKNOWN')
ENDIF
!
!Check if coordinates are already reduced or not
CALL FIND_IF_REDUCED(H,P,isreduced)
WRITE(msg,*) 'isreduced:', isreduced
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Check if some atoms are fixed
IF(ALLOCATED(AUX)) THEN
  DO i=1,SIZE(AUXNAMES)
    IF(AUXNAMES(i)=="fixx") fixx=i
    IF(AUXNAMES(i)=="fixy") fixy=i
    IF(AUXNAMES(i)=="fixz") fixz=i
    IF(AUXNAMES(i)=="magx") magx=i
    IF(AUXNAMES(i)=="magy") magy=i
    IF(AUXNAMES(i)=="magz") magz=i
  ENDDO
ENDIF
!
!Count the number of fixed atoms
IF( fixx.NE.0 .AND. fixy.NE.0 .AND. fixz.NE.0 ) THEN
  DO i=1,SIZE(P(:,1))
    IF( AUX(i,fixx)>0.1d0 .OR. AUX(i,fixy)>0.1d0 .OR. AUX(i,fixz)>0.1d0 ) THEN
      NPinert = NPinert+1
    ENDIF
  ENDDO
ELSE
  NPinert = 0
ENDIF
!
!Set the number of mobile atoms
NPd = SIZE(P(:,1))-NPinert
!
!
!
200 CONTINUE
201 FORMAT(a2,2X,3(f16.8,2X))
WRITE(ofu,'(a10)') "aLat = 1.0"
WRITE(ofu,102) 'a1 = ', H(1,1), H(1,2), H(1,3)
WRITE(ofu,102) 'a2 = ', H(2,1), H(2,2), H(2,3)
WRITE(ofu,102) 'a3 = ', H(3,1), H(3,2), H(3,3)
102 FORMAT(a5,2X,3(f16.8,2X))
!
WRITE(msg,*) "# Natoms = ", SIZE(P,1)
WRITE(ofu,'(a)') TRIM(ADJUSTL(msg))
!
!Write atom coordinates
IF(isreduced) THEN
  WRITE(ofu,'(a14)') 'coord = Direct'
ELSE
  WRITE(ofu,'(a17)') 'coord = Cartesian'
ENDIF
!
!Write atomic coordinates
DO i=1,SIZE(P,1)
  CALL ATOMSPECIES(P(i,4),species)
  WRITE(msg,201) species, P(i,1), P(i,2), P(i,3)
  !
  !If some atoms are fixed, add flags at the end of line
  temp = ""
  !Note: internally if AUX(i,fix)==1 then atom is fixed, if it is 0 it is mobile
  !In BOPfox, T indicates that atom is fixed along a direction, and F that it is free to move
  IF( fixx>0 .OR. fixy>0 .OR. fixz>0 ) THEN
    IF( fixx>0 ) THEN
      IF( AUX(i,fixx)>0.5d0 ) THEN
        temp = "T"
      ELSE
        temp = "F"
      ENDIF
    ELSE
      temp = "F"
    ENDIF
    IF( fixy>0 ) THEN
      IF( AUX(i,fixy)>0.5d0 ) THEN
        temp = TRIM(temp)//" T"
      ELSE
        temp = TRIM(temp)//" F"
      ENDIF
    ELSE
      temp = TRIM(temp)//" F"
    ENDIF
    IF( fixz>0 ) THEN
      IF( AUX(i,fixz)>0.5d0 ) THEN
        temp = TRIM(temp)//" T"
      ELSE
        temp = TRIM(temp)//" F"
      ENDIF
    ELSE
      temp = TRIM(temp)//" F"
    ENDIF
    temp = "   "//TRIM(ADJUSTL(temp))
  ENDIF
  !
  WRITE(ofu,'(a)') TRIM(ADJUSTL(msg))//TRIM(temp)
ENDDO
!
!Write magnetisation section (if values are defined in AUX)
IF( magx>0 .AND. magy>0 .AND. magz>0 ) THEN
  WRITE(ofu,'(a20)') 'magnetisation = true'
  DO i=1,SIZE(AUX,1)
    WRITE(ofu,'(3(f16.6,2X))') AUX(i,magx), AUX(i,magy), AUX(i,magz)
  ENDDO
ENDIF
GOTO 1000
!
!
!
400 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
!
!
!
1000 CONTINUE
IF(ofu.NE.6) THEN
  CLOSE(ofu)
ENDIF
!
END SUBROUTINE WRITE_BOPFOX
!
END MODULE out_bopfox
