MODULE out_mbpp_coorat
!
!
!**********************************************************************************
!*  OUT_MBPP_COORAT                                                               *
!**********************************************************************************
!* This module writes files in the COORAT format, which is used by the            *
!* Mixed-Basis PseudoPotential (MBPP) code by B. Meyer, C. Elsässer,              *
!* F. Lechermann, and M. Fähnle, Max-Plank-Institut für Metalforschung,           *
!* Stuttgart (unpublished).                                                       *
!**********************************************************************************
!* (C) March 2011 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 02 Oct. 2020                                     *
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
SUBROUTINE WRITE_COORAT(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: isreduced
INTEGER:: i, j
REAL(dp):: P1, P2, P3
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: G
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
CALL INVMAT(H,G)
!
100 CONTINUE
msg = 'entering WRITE_COORAT'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Find number of species
CALL FIND_NSP(P(:,4),aentries)
!
!Check if coordinates are already reduced or not
CALL FIND_IF_REDUCED(H,P,isreduced)
!
OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN',ERR=1000)
!
DO j=1,SIZE(aentries(:,1))
  !Write the line concerning this new species
  CALL ATOMSPECIES(aentries(j,1),species)
  WRITE(msg,*) NINT(aentries(j,2))
  WRITE(40,*) 'name='//TRIM(species)//' '//'natom='//TRIM(ADJUSTL(msg))
  !
  !Write atoms positions (must be reduced coordinates)
  DO i=1,SIZE(P(:,1))
    IF( P(i,4) == aentries(j,1) ) THEN
      IF(isreduced) THEN
        WRITE(temp,'(3(f12.8,2X))') P(i,1), P(i,2), P(i,3)
      ELSE
        P1 = P(i,1)
        P2 = P(i,2)
        P3 = P(i,3)
        WRITE(40,110)  P1*G(1,1) + P2*G(2,1) + P3*G(3,1),     &
                       P1*G(1,2) + P2*G(2,2) + P3*G(3,2),     &
                       P1*G(1,3) + P2*G(2,3) + P3*G(3,3)
      ENDIF
    ENDIF
  ENDDO
ENDDO
110 FORMAT(3(f16.8,2X))
CLOSE(40)
msg = "COORAT"
CALL ATOMSK_MSG(3002,(/msg,outputfile/),(/0.d0/))
!
!
!
200 CONTINUE
!Write supercell parameters to INP file
msg = 'INP'
CALL CHECKFILE(msg,'writ')
OPEN(UNIT=41,FILE=msg,STATUS='UNKNOWN',ERR=1000)
WRITE(41,'(a)') TRIM(comment(1))
WRITE(41,'(a18)') '#define VERSION1.0'
WRITE(41,*) ''
WRITE(41,'(a2)') '10'
WRITE(msg,*) SIZE(aentries(:,1))
WRITE(temp,*) MAXVAL(aentries(:,2))
WRITE(msg,*) ' ntype='//TRIM(ADJUSTL(msg))//' natomax='// &
     &       TRIM(ADJUSTL(temp))//' source=file'
WRITE(41,'(a)') TRIM(msg)
WRITE(41,'(a17)') ' alat= 1.00000000'
DO i=1,3
  WRITE(41,'(3(f16.8,2X))') H(i,1), H(i,2), H(i,3)
ENDDO
CLOSE(41)
msg = "INP"
CALL ATOMSK_MSG(3002,(/msg,msg/),(/0.d0/))
!
!
!
1000 CONTINUE
!
END SUBROUTINE WRITE_COORAT
!
END MODULE out_mbpp_coorat
