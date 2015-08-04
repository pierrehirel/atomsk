MODULE out_cel
!
!
!**********************************************************************************
!*  OUT_CEL                                                                       *
!**********************************************************************************
!* This module writes a Dr. Probe super-cell structure definition file (CEL).     *
!* The CEL format is used by the TEM simulation software Dr. Probe and described  *
!* on the software website:                                                       *
!*     http://www.er-c.org/barthel/drprobe/celfile.html                           *
!**********************************************************************************
!* (C) July 2015 - Juri Barthel                                                   *
!*     Gemeinschaftslabor fuer Elektronenmikroskopie                              *
!*     RWTH Aachen (GERMANY)                                                      *
!*     ju.barthel@fz-juelich.de                                                   *
!* Last modification: P. Hirel - 04 Aug. 2015                                     *
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
SUBROUTINE WRITE_CEL(H,P,comment,AUXNAMES,AUX,outputfile)
!
REAL(dp),PARAMETER :: r2d = 180.d0/pi
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=5):: zone
CHARACTER(LEN=8):: date
CHARACTER(LEN=10):: time
CHARACTER(LEN=16):: month, smonth
CHARACTER(LEN=4096):: msg, tempt, temp, tempocc, tempbiso
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
INTEGER:: i, iaux
INTEGER:: occ
INTEGER:: biso
INTEGER:: Naux, NP
INTEGER,DIMENSION(8):: values
REAL(dp):: a, b, c, alpha, beta, gamma !supercell (conventional notation)
REAL(dp):: P1, P2, P3, PO, PB
REAL(dp):: smass_tot !mass of the compound
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: G   !Invert of H
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
CALL INVMAT(H,G)
Naux = 0
IF (ALLOCATED(AUX).AND.ALLOCATED(AUXNAMES)) THEN ! Get number of auxiliary properties
  Naux = SIZE(AUX,2) ! ? Redundant but maybe better MIN(SIZE(AUX,2),SIZE(AUXNAMES))
ENDIF
occ = 0 ! no occupancy information by default
biso = 0 ! no thermal vibration information by default
PO = 1.0d+0 ! 1.0 occupancy of each atomic site by default
PB = 0.0d+0 ! 0.0 thermal vibration parameter (Biso) by default
NP = SIZE(P,1)
IF (Naux>0) THEN
  DO iaux=1, Naux
    SELECT CASE(TRIM(AUXNAMES(iaux)))
    CASE("occ")
      occ = iaux
    CASE("biso")
      biso = iaux
    ENDSELECT
  ENDDO
END IF
!
!Warnings in case of missing occupancy and thermal vibration parameters
IF (occ<=0) THEN
  nwarn = nwarn+1
  CALL ATOMSK_MSG(3711,(/''/),(/0.d0/))
ENDIF
IF (biso<=0) THEN
  nwarn = nwarn+1
  CALL ATOMSK_MSG(3712,(/''/),(/0.d0/))
ENDIF
!
!
msg = 'entering WRITE_CEL'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
!
100 CONTINUE
OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN',ERR=800)
!
!
!Prepare a time stamp
CALL DATE_AND_TIME(DATE, TIME, ZONE, VALUES)
CALL INT2MONTH(VALUES(2),month,smonth)
WRITE(msg,'(i2,a2,i4)') VALUES(3), ", ", VALUES(1)
WRITE(tempt,*) TRIM(smonth)//" "//TRIM(ADJUSTL(msg))
!
!
!Find number of atoms for each species
!Determine formula and mass of the compound
!This is a modified version, which includes the optional
!occupancy parameter of an atomic site in the formula
!and mass calculation.
CALL COMPFORMULA(P,AUXNAMES,AUX,msg,smass_tot)
!
!Write a comment in the first line
WRITE(40,'(a)') comment(1)
!WRITE(40,'(a)') "# Draft CEL file of "//TRIM(ADJUSTL(msg))// &
!     &          " generated by atomsk - "//TRIM(ADJUSTL(tempt))
!
!Write cell vectors (conventional notation, size in nm, angles in degrees)
CALL MATCONV(H,a,b,c,alpha,beta,gamma)
WRITE(temp,'(6(f8.4,2X))') a*0.1,b*0.1,c*0.1,alpha*r2d,beta*r2d,gamma*r2d
WRITE(40,'(a)') " 0  "//TRIM(ADJUSTL(temp))
!
!
!Write atom site list
DO i=1,SIZE(P,1)
  ! Get species string
  CALL ATOMSPECIES(P(i,4),species)
  ! Prepare string with fractional x,y,z atom position in the cell
  P1 = P(i,1)
  P2 = P(i,2)
  P3 = P(i,3)
  WRITE(temp,'(3(f9.6,1X))')  P1*G(1,1) + P2*G(2,1) + P3*G(3,1),     &
                           &  P1*G(1,2) + P2*G(2,2) + P3*G(3,2),     &
                           &  P1*G(1,3) + P2*G(2,3) + P3*G(3,3)
  ! add strings of auxiliary properties to temp
  IF (occ>0) PO = AUX(i,occ)
  WRITE(tempocc,'(f9.6)') PO
  IF (biso>0) PB = AUX(i,biso)*0.01 ! from A**2 to nm**2
  WRITE(tempbiso,'(f9.6)') PB
  ! combine the strings and add 3 not used entries
  WRITE(temp,*) species//'  '//TRIM(ADJUSTL(temp))//'  '//           &
              & TRIM(ADJUSTL(tempocc))//'  '//                       &
              & TRIM(ADJUSTL(tempbiso))//                            &
              & '  0.000000  0.000000  0.000000'
  ! write to file
  WRITE(40,'(a)') TRIM(ADJUSTL(temp))
ENDDO
!
!Write end line
WRITE(40,'(a1)') "*"
!
!
!
200 CONTINUE
CLOSE(40)
msg = "CEL"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
GOTO 1000
!
!
!
800 CONTINUE
nerr=nerr+1
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE WRITE_CEL
!
END MODULE out_cel
