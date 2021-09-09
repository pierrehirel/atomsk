MODULE out_xsf
!
!
!**********************************************************************************
!*  OUT_XSF                                                                       *
!**********************************************************************************
!* This module writes XSF format, initially designed for XCrysDen.                *
!* The XSF format is described here:                                              *
!*    http://www.xcrysden.org/doc/XSF.html                                        *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 31 May 2021                                      *
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
SUBROUTINE WRITE_XSF(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: forces !are there forces or not?
INTEGER:: fx, fy, fz !index for forces in AUX
INTEGER:: i
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
forces = .FALSE.
fx = 0
fy = 0
fz = 0
!
msg = 'entering WRITE_XSF'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!Check if forces are present in auxiliary properties
IF(ALLOCATED(AUXNAMES)) THEN
  DO i=1,SIZE(AUXNAMES)
    IF( TRIM(ADJUSTL(AUXNAMES(i)))=='fx') THEN
      fx = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='fy') THEN
      fy = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='fz') THEN
      fz = i
    ENDIF
  ENDDO
  IF( fx.NE.0 .AND. fy.NE.0 .AND. fz.NE.0 ) THEN
    forces = .TRUE.
  ENDIF
ENDIF
WRITE(temp,*) 'forces: ', forces
CALL ATOMSK_MSG(999,(/temp/),(/0.d0/))
!
!
!
100 CONTINUE
IF(ofu.NE.6) THEN
  OPEN(UNIT=ofu,FILE=outputfile,STATUS='UNKNOWN',ERR=500)
ENDIF
!
!Write header of XSF file
WRITE(ofu,*) TRIM(ADJUSTL(comment(1)))
WRITE(ofu,*) "CRYSTAL"
WRITE(ofu,*) "PRIMVEC"
WRITE(ofu,105) H(1,1), H(1,2), H(1,3)
WRITE(ofu,105) H(2,1), H(2,2), H(2,3)
WRITE(ofu,105) H(3,1), H(3,2), H(3,3)
WRITE(ofu,*) "CONVVEC"
WRITE(ofu,105) H(1,1), H(1,2), H(1,3)
WRITE(ofu,105) H(2,1), H(2,2), H(2,3)
WRITE(ofu,105) H(3,1), H(3,2), H(3,3)
WRITE(ofu,*) "PRIMCOORD"
WRITE(ofu,'(i9,a2)') SIZE(P(:,1)), ' 1'
105 FORMAT(3(f16.8,2X))
!
!Write atom positions
DO i=1,SIZE(P(:,1))
  WRITE(temp,110) INT(P(i,4)), P(i,1), P(i,2), P(i,3)
  !If forces are present, add them to the line
  IF(forces) THEN
    WRITE(msg,'(3(e16.8,2X))') AUX(i,fx), AUX(i,fy), AUX(i,fz)
    temp = TRIM(ADJUSTL(temp))//'  '//TRIM(ADJUSTL(msg))
  ENDIF
  !Write line to the file
  WRITE(ofu,'(a)') TRIM(ADJUSTL(temp))
ENDDO
110 FORMAT(i5,2X,3(f16.8,2X))
!
!
!
500 CONTINUE
IF(ofu.NE.6) THEN
  CLOSE(40)
ENDIF
msg = "XSF"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
1000 CONTINUE
!
END SUBROUTINE WRITE_XSF
!
END MODULE out_xsf
