MODULE out_dd
!
!**********************************************************************************
!*  OUT_DD                                                                        *
!**********************************************************************************
!* This module reads in two arrays containing atomic positions, and               *
!* writes a file in the ddplot (.dd) and PLT formats.                             *
!* These file formats are described for instance at:                              *
!*     http://groger.ipm.cz/download/ddplot/ddplot.html#formats                   *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 06 April 2020                                    *
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
SUBROUTINE WRITE_DD(H,Pfirst,Psecond,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=4096):: msg, temp
INTEGER:: i
!REAL(dp):: alat
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: Pfirst, Psecond
!
!
!Initialize variables
!
!
!If arrays do not contain the same number of atoms, exit
IF( SIZE(Pfirst,1) .NE. SIZE(Psecond,1) ) THEN
  CALL ATOMSK_MSG(806,(/outputfile/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Display a message if number of atoms is too big for ddplot
IF(SIZE(Pfirst,1)>6000) THEN
  nwarn = nwarn+1
  CALL ATOMSK_MSG(3707,(/''/),(/DBLE(SIZE(Pfirst(:,1)))/))
ENDIF
!
!
!
100 CONTINUE
!CALL NAME_OUTFILE(outputfile,outputfile,"dd   ")
OPEN(UNIT=40,FILE=outputfile,FORM='FORMATTED',STATUS='UNKNOWN')
!Write header of DD file
WRITE(40,'(a4)') 'CSYS'
WRITE(40,'(a5)') '1 0 0'
WRITE(40,'(a5)') '0 1 0'
WRITE(40,'(a5)') '0 0 1'
WRITE(40,'(a6)') 'PERIOD'
WRITE(40,'(3(f10.6,1X))') H(1,1), H(2,2), H(3,3)
WRITE(40,'(a12)') 'DISLO_CENTER'
WRITE(40,*) '0.000 0.000'
WRITE(40,'(a9)') 'NUM_UNREL'
WRITE(temp,'(i16)') SIZE(Pfirst,1)
WRITE(40,*) TRIM(ADJUSTL(temp))
WRITE(40,'(a10)') 'COOR_UNREL'
!
!Write positions of ideal system
DO i=1,SIZE(Pfirst,1)
  WRITE(40,110) Pfirst(i,1), Pfirst(i,2), Pfirst(i,3), Pfirst(i,4)
ENDDO
!
!Write second section of DD file
WRITE(40,'(a7)') 'NUM_REL'
WRITE(temp,'(i16)') SIZE(Psecond,1)
WRITE(40,*) TRIM(ADJUSTL(temp))
WRITE(40,'(a8)') 'COOR_REL'
!
!Write positions of relaxed system
DO i=1,SIZE(Psecond,1)
  WRITE(40,110) Psecond(i,1), Psecond(i,2), Psecond(i,3)
ENDDO
110 FORMAT(3(f16.8,1X),i3)
!
CLOSE(40)
msg = "ddplot"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
300 CONTINUE
!Write PLT file
!This is not portable so it was commented out
!
! alat = 3.905d0
! CALL NAME_OUTFILE(outputfile,outputfile,"plt")
! OPEN(UNIT=42,FILE=outputfile,FORM='FORMATTED',STATUS='UNKNOWN')
! !First we need to convert to fractional coordinates
! !CALL CART2FRAC(Pfirst,H)
! !CALL CART2FRAC(Psecond,H)
! !
! !Number of atoms in relaxed system
! WRITE(42,'(i5)') SIZE(Pfirst(:,1))
! !Z-coordinates of atoms in relaxed system
! DO i=1,SIZE(Psecond,1)
!   WRITE(42,110) Psecond(i,3)/alat
! ENDDO
! !Atom positions in relaxed system
! DO i=1,SIZE(Psecond(:,1))
!   WRITE(42,110) Psecond(i,1)/alat, Psecond(i,2)/alat, Psecond(i,4)
! ENDDO
! !
! WRITE(42,'(i5)') SIZE(Pfirst(:,1))
! !Z-coordinates of atoms in ideal system
! DO i=1,SIZE(Psecond,1)
!   WRITE(42,110) Pfirst(i,3)
! ENDDO
! !Atom positions in ideal system
! DO i=1,SIZE(Pfirst(:,1))
!   WRITE(42,110) Pfirst(i,1)/alat, Pfirst(i,2)/alat
! ENDDO
! !
! WRITE(42,'(a1)') "0"  !unused, don't ask me why
! !
! WRITE(42,'(f10.6,1X)') H(1,1)
! WRITE(42,'(f10.6,1X)') H(2,2)
! WRITE(42,'(f10.6,1X)') H(3,3)
! !
! !WRITE(42,*) "CORE {0.0;0.0}"   !(x,y) position of the screw dislocation
! !
! CLOSE(42)
! msg = "plt"
! CALL ATOMSK_MSG(3002,(/msg,outputfile/),(/0.d0/))
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE WRITE_DD
!
END MODULE out_dd
