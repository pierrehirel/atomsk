MODULE out_str
!
!**********************************************************************************
!*  OUT_STR                                                                       *
!**********************************************************************************
!* This module write files in the PDFFIT structure format (*.str or *.stru).      *
!* The STR format is described on the Website of PDFFIT:                          *
!*   https://web.pa.msu.edu/cmp/billinge-group/programs/discus/pdffit.html        *
!**********************************************************************************
!* (C) February 2018 - Pierre Hirel                                               *
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
SUBROUTINE WRITE_STR(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=12):: test
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: isreduced
INTEGER:: i, j
INTEGER:: occ, U11, U22, U33, U12, U13, U23, Biso !index of these properties in AUX
REAL(dp):: a, b, c, alpha, beta, gamma !supercell (conventional notation)
REAL(dp):: P1, P2, P3
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: G
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
isreduced = .FALSE.
G(:,:) = 0.d0
i=1
occ = 0
U11 = 0
U22 = 0
U33 = 0
U12 = 0
U13 = 0
U23 = 0
!
100 CONTINUE
msg = 'entering WRITE_STR'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Check if coordinates are already reduced or not
CALL FIND_IF_REDUCED(H,P,isreduced)
WRITE(msg,*) 'isreduced:', isreduced
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
IF( .NOT.isreduced ) THEN
  !Calculate the inverse of matrix H
  msg = 'inverting matrix H'
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  CALL INVMAT(H,G)
ENDIF
!
!Check if properties are present
IF( ALLOCATED(AUXNAMES) .AND. SIZE(AUXNAMES)>=1 ) THEN
  temp = ""
  DO i=1,SIZE(AUXNAMES)
    IF( TRIM(ADJUSTL(AUXNAMES(i)))=='occ') THEN
      occ = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='Biso') THEN
      Biso = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='U11') THEN
      U11 = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='U22') THEN
      U22 = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='U33') THEN
      U33 = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='U12') THEN
      U12 = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='U13') THEN
      U13 = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='U23') THEN
      U23 = i
    ENDIF
  ENDDO
ENDIF
!
IF(ofu.NE.6) THEN
  OPEN(UNIT=ofu,FILE=outputfile,STATUS='UNKNOWN',ERR=250)
ENDIF
!
! Write header of STR file
CONTINUE
!Remove leading # from comment
temp = ADJUSTL(comment(1))
IF( temp(1:1)=="#" ) THEN
  temp(1:1) = " "
ENDIF
WRITE(ofu,'(a)') "title  "//TRIM(ADJUSTL(temp))
WRITE(ofu,'(a)') "format pdffit"
WRITE(ofu,'(a)') "scale  1.0000"
WRITE(ofu,'(a)') "spcgr  P1"
!Convert cell vectors (conventional notation, size in angströms, angles in degrees)
CALL MATCONV(H,a,b,c,alpha,beta,gamma)
WRITE(temp,'(6(f12.6,a3))') a, ",  ", b, ",  ", c, ",  ", &
                         & RAD2DEG(alpha), ",  ", RAD2DEG(beta), ",  ", RAD2DEG(gamma)
WRITE(ofu,'(a)') "cell   "//TRIM(ADJUSTL(temp))
!
WRITE(ofu,'(a)') "dcell  0.000000,  0.000000,  0.000000,  0.000000,  0.000000,  0.000000"
!
!Write number of replicas (here always 1x1x1) and number of atoms
WRITE(temp,*) SIZE(P,1)
WRITE(ofu,'(a)') "ncell  1,1,1,"//TRIM(ADJUSTL(temp))
!
!Write atom positions
WRITE(ofu,'(a)') "atoms"
DO i=1,SIZE(P,1)
  !Get atom species
  CALL ATOMSPECIES(P(i,4),species)
  !STR format uses all-capital letters
  species(1:1) = StrUpCase(species(1:1))
  species(2:2) = StrUpCase(species(2:2))
  !
  !Coordinates
  IF(isreduced) THEN
    WRITE(temp,'(a2,3X,3(f12.8,2X))') species, P(i,1), P(i,2), P(i,3)
  ELSE
    P1 = P(i,1)
    P2 = P(i,2)
    P3 = P(i,3)
    WRITE(temp,'(a2,2X,3(f12.8,6X))')  species, &
                              &  P1*G(1,1) + P2*G(2,1) + P3*G(3,1),     &
                              &  P1*G(1,2) + P2*G(2,2) + P3*G(3,2),     &
                              &  P1*G(1,3) + P2*G(2,3) + P3*G(3,3)
  ENDIF
  !
  !Append the auxiliray properties if defined
  IF( occ>0 ) THEN
    WRITE(msg,'(f9.4)') AUX(i,occ)
    temp = TRIM(ADJUSTL(temp))//'       '//TRIM(ADJUSTL(msg))
  ELSE
    !Default occupancy = 1
    temp = TRIM(ADJUSTL(temp))//'       1.0000'
  ENDIF
  !
  !Write the line to the file
  WRITE(ofu,'(a)') TRIM(ADJUSTL(temp))
  !
  !Next line contains standard deviations on position and occupancy: here we just write zeros
  WRITE(ofu,'(6X,a)') "0.00000000        0.00000000        0.00000000       0.0000"
  !
  !Next line contains anisotropy factors U11, U22, U33
  temp = ""
  IF( U11>0 ) THEN
    WRITE(temp,'(f12.8)') AUX(i,U11)
  ELSE
    WRITE(temp,'(a10)') "0.00000000"
  ENDIF
  IF( U22>0 ) THEN
    WRITE(msg,'(f12.8)') AUX(i,U22)
  ELSE
    WRITE(msg,'(a10)') "0.00000000"
  ENDIF
  temp = TRIM(ADJUSTL(temp))//"        "//TRIM(ADJUSTL(msg))
  IF( U33>0 ) THEN
    WRITE(msg,'(f12.8)') AUX(i,U33)
  ELSE
    WRITE(msg,'(a10)') "0.00000000"
  ENDIF
  temp = TRIM(ADJUSTL(temp))//"        "//TRIM(ADJUSTL(msg))
  !Write this line to the file
  WRITE(ofu,'(6X,a)') TRIM(ADJUSTL(temp))
  !
  !Next line contains standard deviations on U11,U22,U33: here we just write zeros
  WRITE(ofu,'(6X,a)') "0.00000000        0.00000000        0.00000000"
  !
  !Next line contains anisotropy factors U12, U13, U23
  temp = ""
  IF( U12>0 ) THEN
    WRITE(temp,'(f12.8)') AUX(i,U12)
  ELSE
    WRITE(temp,'(a10)') "0.00000000"
  ENDIF
  IF( U13>0 ) THEN
    WRITE(msg,'(f12.8)') AUX(i,U13)
  ELSE
    WRITE(msg,'(a10)') "0.00000000"
  ENDIF
  temp = TRIM(ADJUSTL(temp))//"        "//TRIM(ADJUSTL(msg))
  IF( U23>0 ) THEN
    WRITE(msg,'(f12.8)') AUX(i,U23)
  ELSE
    WRITE(msg,'(a10)') "0.00000000"
  ENDIF
  temp = TRIM(ADJUSTL(temp))//"        "//TRIM(ADJUSTL(msg))
  !Write this line to the file
  WRITE(ofu,'(6X,a)') TRIM(ADJUSTL(temp))
  !
  !Next line contains standard deviations on U12,U13,U23: here we just write zeros
  WRITE(ofu,'(6X,a)') "0.00000000        0.00000000        0.00000000"
  !
ENDDO
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
msg = "STR"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
1000 CONTINUE
IF(ofu.NE.6) THEN
  CLOSE(40)
ENDIF
!
END SUBROUTINE WRITE_STR
!
!
END MODULE out_str
