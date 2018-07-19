MODULE out_moldy
!
!
!**********************************************************************************
!*  OUT_MOLDY                                                                     *
!**********************************************************************************
!* This module writes a system file for MOLDY (by default called "system.in").    *
!* The format is described for instance at:                                       *
!*     https://www.wiki.ed.ac.uk/display/ComputerSim/MOLDY+Information            *
!**********************************************************************************
!* (C) Oct. 2011 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 26 March 2014                                    *
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
SUBROUTINE WRITE_MOLDY(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: isreduced
INTEGER:: i, j
REAL(dp):: P1, P2, P3
REAL(dp):: smass  !atom mass
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: G   !Inverse of H
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P !atom positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
species=''
G(:,:) = 0.d0
isreduced = .FALSE.
!
msg = 'entering WRITE_MOLDY'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
!
100 CONTINUE
!Check if coordinates are reduced or cartesian
CALL FIND_IF_REDUCED(P,isreduced)
CALL INVMAT(H,G)
!
OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN',ERR=500)
!First line: number of particles
WRITE(msg,*) SIZE(P(:,1))
WRITE(40,*) TRIM(ADJUSTL(msg))
!
!Second line: number of replicas along X, Y, Z
WRITE(40,'(a5)') "1 1 1"
!
!Three next lines: supercell vectors
DO i=1,3
  WRITE(40,'(3(f16.8,1X))') (H(i,j),j=1,3)
ENDDO
!
!
!
200 CONTINUE
!Write atomic positions
DO i=1,SIZE(P,1)
  CALL ATOMSPECIES(P(i,4),species)
  CALL ATOMMASS(species,smass)
  IF(isreduced) THEN
    WRITE(temp,220) (P(i,j),j=1,3), NINT(P(i,4)), smass
  ELSE
    P1 = P(i,1)
    P2 = P(i,2)
    P3 = P(i,3)
    WRITE(temp,220)  P1*G(1,1) + P2*G(2,1) + P3*G(3,1),     &
                  &  P1*G(1,2) + P2*G(2,2) + P3*G(3,2),     &
                  &  P1*G(1,3) + P2*G(2,3) + P3*G(3,3),     &
                  &  NINT(P(i,4)), smass
  ENDIF
  WRITE(40,'(a)') TRIM(ADJUSTL(temp))
ENDDO
220 FORMAT(3(f16.8,1X),i3,1X,f6.2)
!
!
!
500 CONTINUE
CLOSE(40)
msg = "MOLDY"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
1000 CONTINUE
!
END SUBROUTINE WRITE_MOLDY
!
END MODULE out_moldy
