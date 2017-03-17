MODULE out_jems
!
!
!**********************************************************************************
!*  OUT_JEMS                                                                      *
!**********************************************************************************
!* This module writes a crystal structure for JEMS.                               *
!* The JEMS format is used by the TEM simulation software JEMS, and is described  *
!* on the software website:                                                       *
!*     http://www.jems-saas.ch/                                                   *
!**********************************************************************************
!* (C) October 2015 - Pierre Hirel                                                *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 01 March 2017                                    *
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
SUBROUTINE WRITE_JEMS(H,P,comment,AUXNAMES,AUX,outputfile)
!
REAL(dp),PARAMETER :: r2d = 180.d0/pi
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: temp1, temp2, temp3
CHARACTER(LEN=1096):: msg, tempi, temp, tempabs, tempocc, tempdw
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
INTEGER:: i, iaux
INTEGER:: absorption, occ, dw
INTEGER:: NP
REAL(dp):: a, b, c, alpha, beta, gamma !supercell (conventional notation)
REAL(dp):: P1, P2, P3, PA, PO, PB
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: G   !Invert of H
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
absorption = 0  ! Assume that no absorption data is available
occ = 0         ! Assume that no occupancy information is available
dw = 0          ! Assume that no thermal vibration information is available
PA = 0.03d0     ! Default absorption coefficient
PO = 1.0d0      ! Default occupancy of each atomic site
PB = 0.0d0      ! Default Debye-Waller parameter
NP = SIZE(P,1)  ! Number of atoms in the system
!
!Look for absorption, occupancy, Debye-Waller param. in the auxiliary properties
IF( ALLOCATED(AUXNAMES) .AND. SIZE(AUXNAMES)>0 ) THEN
  DO iaux=1,SIZE(AUXNAMES)
    SELECT CASE(TRIM(AUXNAMES(iaux)))
    CASE("occ")
      occ = iaux
    CASE("Debye-Waller")
      dw = iaux
    CASE("absorption")
      absorption = iaux
    ENDSELECT
  ENDDO
END IF
!
!Warnings in case of missing occupancy or thermal vibration parameters
IF (occ<=0) THEN
  nwarn = nwarn+1
  CALL ATOMSK_MSG(3711,(/''/),(/0.d0/))
ENDIF
IF (dw<=0) THEN
  nwarn = nwarn+1
  CALL ATOMSK_MSG(3712,(/''/),(/0.d0/))
ENDIF
!
!
msg = 'entering WRITE_JEMS'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
!
100 CONTINUE
OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN',ERR=800)
!
! Write header of JEMS file
CALL GETCWD(msg)
WRITE(40,'(a)') "file|"//TRIM(ADJUSTL(msg))//pathsep//TRIM(ADJUSTL(outputfile))
WRITE(40,'(a16)') "system|triclinic"
WRITE(40,'(a21)') "HMSymbol|1|1|0|0| P 1"
WRITE(40,'(a16)') "rps|0| x , y , z"
!
! Write cell vectors (conventional notation, angles in degrees)
CALL MATCONV(H,a,b,c,alpha,beta,gamma)
WRITE(temp,'(f16.4)') a
WRITE(40,'(a)') "lattice|0|"//TRIM(ADJUSTL(temp))
WRITE(temp,'(f16.4)') b
WRITE(40,'(a)') "lattice|1|"//TRIM(ADJUSTL(temp))
WRITE(temp,'(f16.4)') c
WRITE(40,'(a)') "lattice|2|"//TRIM(ADJUSTL(temp))
WRITE(temp,'(f16.4)') alpha*r2d
WRITE(40,'(a)') "lattice|3|"//TRIM(ADJUSTL(temp))
WRITE(temp,'(f16.4)') beta*r2d
WRITE(40,'(a)') "lattice|4|"//TRIM(ADJUSTL(temp))
WRITE(temp,'(f16.4)') gamma*r2d
WRITE(40,'(a)') "lattice|5|"//TRIM(ADJUSTL(temp))
!
!Invert matrix of cell vectors
CALL INVMAT(H,G)
!
!Write atom site list
DO i=1,SIZE(P,1)
  WRITE(tempi,*) i-1  !JEMS starts counting at zero, not 1
  ! Get species string
  CALL ATOMSPECIES(P(i,4),species)
  ! Prepare string with fractional x,y,z atom position in the cell
  P1 = P(i,1)
  P2 = P(i,2)
  P3 = P(i,3)
  WRITE(temp1,'(f16.4)')  P1*G(1,1) + P2*G(2,1) + P3*G(3,1)
  WRITE(temp2,'(f16.4)')  P1*G(1,2) + P2*G(2,2) + P3*G(3,2)
  WRITE(temp3,'(f16.4)')  P1*G(1,3) + P2*G(2,3) + P3*G(3,3)
  ! add strings of auxiliary properties to temp
  IF (occ>0) PO = AUX(i,occ)
  WRITE(tempocc,'(f9.6)') PO
  IF (dw>0) PB = AUX(i,dw)
  WRITE(tempdw,'(f9.6)') PB
  IF (absorption>0) PA = AUX(i,dw)
  WRITE(tempabs,'(f9.6)') PA
  !
  !Format:  atom|i|XX,a,x,y,z,Debye-Waller,occupancy,absorption
  !Example: atom|0|Ga,a,0.0000,0.0000,0.0000,0.0050,1.0000,0.0520
  WRITE(40,'(a)') "atom|"//TRIM(ADJUSTL(tempi))//"|"//TRIM(species)//",a,"//         &
                & TRIM(ADJUSTL(temp1))//","//TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))// &
                & ","//TRIM(ADJUSTL(tempdw))//","//TRIM(ADJUSTL(tempocc))//","//TRIM(ADJUSTL(tempabs))
ENDDO
!
!
!
200 CONTINUE
CLOSE(40)
msg = "JEMS"
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
END SUBROUTINE WRITE_JEMS
!
END MODULE out_jems
