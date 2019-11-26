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
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
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
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
!
!Initialize variables
Nunskewed = 0
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
DO i=1,3
  skew(2:2)='x'
  IF(i==2) skew(2:2)='y'
  IF(i==3) skew(2:2)='z'
  !
  DO j=3,1,-1
    !Don't consider diagonal elements
    IF(i.NE.j) THEN
      skew(1:1)='x'
      IF(j==2) skew(1:1)='y'
      IF(j==3) skew(1:1)='z'
      !
      !Unskew tilt H(i,j)
      tiltbefore = H(i,j)
      iloop=0
      !If tilt is too large, remove the matching box vector
      DO WHILE( H(i,j)>0.5d0*H(j,j) )
        H(i,:) = H(i,:) - H(j,:)
        iloop=iloop+1
        IF(iloop>100) EXIT
      ENDDO
      !If tilt is too negative, add the matching box vector
      DO WHILE( H(i,j)<-0.5d0*H(j,j) )
        H(i,:) = H(i,:) + H(j,:)
        iloop=iloop+1
        IF(iloop>100) EXIT
      ENDDO
      !Check that the loops did not go crazy
      IF(iloop>100) THEN
        !After 100 loops no solution was found
        !Display a warning and restore initial value
        nwarn=nwarn+1
        CALL ATOMSK_MSG(3706,(/skew/),(/0.d0/))
        H(i,j) = tiltbefore
      ELSEIF(iloop>0) THEN
        !This tilt was corrected: display message
        CALL ATOMSK_MSG(3004,(/skew/),(/0.d0/))
        Nunskewed=Nunskewed+1
      ENDIF
    ENDIF
  ENDDO
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
