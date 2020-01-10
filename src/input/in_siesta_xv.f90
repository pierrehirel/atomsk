MODULE in_siesta_xv
!
!
!**********************************************************************************
!*  IN_SIESA_XV                                                                   *
!**********************************************************************************
!* This module reads XV format, initially designed for SIESTA.                    *
!**********************************************************************************
!* This module reads SIESTA XV format.                                            *
!* The XV format is described in the SIESTA manual:                               *
!*    http://www.icmab.es/dmmis/leem/siesta/                                      *
!**********************************************************************************
!* (C) Jan. 2012 - Eva Marie Kalivoda                                             *
!*     Fraunhofer Institute für Werkstoffmechanik IWM                             *
!*     Wöhlerstr. 11, 79108 Freiburg im Breisgau                                  *
!* Last modification: P. Hirel - 08 Jan. 2020                                     *
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
SUBROUTINE READ_XV(inputfile,H,P,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
INTEGER:: i, NP
REAL(dp):: a
REAL(dp), DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P
REAL(dp), DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
NP = 0
!
100 CONTINUE
msg = 'entering READ_XV'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
READ(30,*,ERR=400,END=400) H(1,1), H(1,2), H(1,3)
READ(30,*,ERR=400,END=400) H(2,1), H(2,2), H(2,3)
READ(30,*,ERR=400,END=400) H(3,1), H(3,2), H(3,3)
READ(30,*,ERR=400,END=400) a
IF( a > NATOMS_MAX ) THEN
  nerr = nerr+1
  CALL ATOMSK_MSG(821,(/""/),(/a/))
  GOTO 1000
ENDIF
NP = NINT(a)
!
ALLOCATE(P(NP,4) , STAT=i)
IF( i>0 ) THEN
  ! Allocation failed (not enough memory)
  nerr = nerr+1
  CALL ATOMSK_MSG(819,(/''/),(/DBLE(NP)/))
  GOTO 1000
ENDIF
P(:,:) = 0.d0
!
ALLOCATE(AUX(NP,4))
AUX(:,:)=0.d0
ALLOCATE(AUXNAMES(4))
AUXNAMES(1)="type"
AUXNAMES(2)="vx"
AUXNAMES(3)="vy"
AUXNAMES(4)="vz"
!
DO i=1,NP
  READ(30,*,ERR=400,END=400) AUX(i,1), P(i,4), P(i,1), P(i,2), P(i,3), &
                           & AUX(i,2), AUX(i,3), AUX(i,4)
ENDDO
GOTO 1000
!
400 CONTINUE
CALL ATOMSK_MSG(1801,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
1000 CONTINUE
CLOSE(30)
!
END SUBROUTINE READ_XV
!
END MODULE in_siesta_xv
