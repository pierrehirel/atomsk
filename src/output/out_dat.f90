MODULE out_dat
!
!**********************************************************************************
!*  OUT_DAT                                                                       *
!**********************************************************************************
!* This module writes a file in "data file" format.                               *
!* This format has no standard specification, it just contains columns of data.   *
!* Here it will only contain the atomic number and atom position,                 *
!* followed by auxiliary properties if any.                                       *
!* Box dimensions and positions of shells (if any) are ignored.                   *
!* This type of data file can be visualized e.g. with Gnuplot.                    *
!**********************************************************************************
!* (C) April 2022 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 16 April 2024                                    *
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
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE WRITE_DAT(H,P,S,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=1):: fs=' '  !field separator. Default is a blank space. Change it here if you want something else
CHARACTER(LEN=128):: v1,v2,v3,v4  !to store values
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
INTEGER:: i, j
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: S
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
!
!
100 CONTINUE
msg = 'entering WRITE_DAT'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
IF(ofu.NE.6) THEN
  OPEN(UNIT=ofu,FILE=outputfile,STATUS='UNKNOWN',ERR=250)
ENDIF
!
!
! Write coordinates in DAT file
DO i=1,SIZE(P,1)
  !Get position and atomic number of current atom
  WRITE(v1,'(f16.8)') P(i,1)
  WRITE(v2,'(f16.8)') P(i,2)
  WRITE(v3,'(f16.8)') P(i,3)
  WRITE(v4,'(f16.8)') P(i,4)
  !Write that into the line
  msg = TRIM(ADJUSTL(v4))//fs//TRIM(ADJUSTL(v1))//fs//TRIM(ADJUSTL(v2))//fs//TRIM(ADJUSTL(v3))
  !
  IF( ALLOCATED(AUXNAMES) .AND. ALLOCATED(AUX) .AND. SIZE(AUX,1)==SIZE(P,1) .AND. SIZE(AUX,2)>0 ) THEN
    DO j=1,SIZE(AUX,2)
      WRITE(v1,'(f16.8)') AUX(i,j)
      !Write that into the line
      msg = TRIM(ADJUSTL(msg))//fs//TRIM(ADJUSTL(v1))
    ENDDO
  ENDIF
  !
  !Write the line to the file
  WRITE(ofu,'(a)') TRIM(ADJUSTL(msg))
END DO
GOTO 300
!
250 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
300 CONTINUE
msg = "DAT"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
1000 CONTINUE
IF(ofu.NE.6) THEN
  CLOSE(ofu)
ENDIF
!
END SUBROUTINE WRITE_DAT
!
!
END MODULE out_dat
