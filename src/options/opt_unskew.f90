MODULE unskew
!
!**********************************************************************************
!*  UNSKEW                                                                        *
!**********************************************************************************
!* This module reads a supercell and "unskews" it.                                *
!* The "skew parameters" of the supercell are the non-diagonal elements.          *
!* They represent the tilt (or inclination) of the supercell vectors              *
!* along X, Y and Z. If a skew parameter is larger than its main component        *
!* (i.e. if |H(i,j)| > H(j,j) for i.NE.j) then this module will "unskew" it,      *
!* i.e. it will add +/-H(j,j) until it becomes smaller than H(j,j).               *
!**********************************************************************************
!* (C) February 2012 - Pierre Hirel                                               *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 25 Sept. 2013                                    *
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
!
CONTAINS
!
!
SUBROUTINE UNSKEW_XYZ(H)
!
!
IMPLICIT NONE
CHARACTER(LEN=2):: skew
CHARACTER(LEN=128):: msg
INTEGER:: i, j, iloop
INTEGER:: Nunskewed !number of 
REAL(dp):: tiltbefore
REAL(dp),DIMENSION(6):: tilt
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
!
!Initialize variables
Nunskewed = 0
tilt(:) = 0.d0
!
WRITE(msg,*) 'Entering UNSKEW_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(2103,(/TRIM(msg)/),(/0.d0/))
!
!If all skew parameters are zero (rectangular box),
!then there is nothing to do => exit
IF( H(2,1)==0.d0 .AND. H(3,1)==0.d0 .AND. H(3,2)==0.d0 .AND. &
  & H(1,2)==0.d0 .AND. H(1,3)==0.d0 .AND. H(2,3)==0.d0       ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2743,(/TRIM(msg)/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
100 CONTINUE
tilt(1) = H(2,1)
tilt(2) = H(3,1)
tilt(3) = H(3,2)
tilt(4) = H(1,2)
tilt(5) = H(1,3)
tilt(6) = H(2,3)
DO i=1,6
  IF(i==1) THEN
    j = 1
    skew = 'xy'
  ELSEIF(i==2) THEN
    j = 1
    skew = 'xz'
  ELSEIF(i==3) THEN
    j = 2
    skew = 'yz'
  ELSEIF(i==4) THEN
    j = 2
    skew = 'yx'
  ELSEIF(i==5) THEN
    j = 3
    skew = 'zx'
  ELSEIF(i==6) THEN
    j = 3
    skew = 'zy'
  ENDIF
  !
  IF( DABS(tilt(i))>0.5d0*H(j,j) ) THEN
    tiltbefore = tilt(i)
    iloop=0
    DO WHILE( tilt(i)>0.5d0*H(j,j) )
      tilt(i) = tilt(i)-H(j,j)
      iloop=iloop+1
      IF(iloop>100) EXIT
    ENDDO
    !
    DO WHILE( tilt(i)<=-0.5d0*H(j,j) )
      tilt(i) = tilt(i)+H(j,j)
      iloop=iloop+1
      IF(iloop>200) EXIT
    ENDDO
    !
    IF(iloop>100) THEN
      nwarn=nwarn+1
      CALL ATOMSK_MSG(3706,(/skew/),(/0.d0/))
      tilt(i) = tiltbefore
    ELSE
      CALL ATOMSK_MSG(3004,(/skew/),(/0.d0/))
      Nunskewed=Nunskewed+1
    ENDIF
    !
  ENDIF
ENDDO
!
CALL ATOMSK_MSG(2104,(/''/),(/DBLE(Nunskewed)/))
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE UNSKEW_XYZ
!
!
!
END MODULE unskew