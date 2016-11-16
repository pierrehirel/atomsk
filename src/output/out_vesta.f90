MODULE out_vesta
!
!
!**********************************************************************************
!*  OUT_VESTA                                                                     *
!**********************************************************************************
!* This module writes VESTA format, designed for visualization with VESTA.        *
!* The VESTA format is described here:                                            *
!*    http://jp-minerals.org/vesta/en/doc.html                                    *
!**********************************************************************************
!* (C) March 2015 - Pierre Hirel                                                  *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 28 Oct. 2016                                     *
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
SUBROUTINE WRITE_VESTA(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: vectors !are there vectors to display or not?
LOGICAL:: isreduced  !are coordinates already reduced?
INTEGER:: vx, vy, vz, occ !index for vector and occupancies in AUX
INTEGER:: i, j, k
REAL(dp):: a, b, c, alpha, beta, gamma   !supercell (conventional notation)
REAL(dp):: P1, P2, P3
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: G              !Inverse of H
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
vectors = .FALSE.
isreduced = .FALSE.
G(:,:) = 0.d0
vx=0
vy=0
vz=0
occ=0
!
msg = 'entering WRITE_VESTA'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Check if coordinates are already reduced or not
CALL FIND_IF_REDUCED(P,isreduced)
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
!Check if a vector quantity is present in auxiliary properties
IF(ALLOCATED(AUXNAMES)) THEN
  temp = ""
  DO i=1,SIZE(AUXNAMES)
    IF( TRIM(ADJUSTL(AUXNAMES(i)))=='occ') THEN
      occ = i
    ENDIF
    !
    IF( .NOT.vectors ) THEN
      k = SCAN(AUXNAMES(i),"x")
      IF( k>0 ) THEN
        vx = i
        temp = AUXNAMES(i)
        temp = ADJUSTL(temp(1:k-1))
        DO j=i+1,SIZE(AUXNAMES)
          IF( TRIM(ADJUSTL(AUXNAMES(j)))==TRIM(temp)//'y' ) THEN
            vy = j
          ELSEIF( TRIM(ADJUSTL(AUXNAMES(j)))==TRIM(temp)//'z') THEN
            vz = j
          ENDIF
        ENDDO
      ENDIF
      !
      IF( vx>0 .AND. vy>0 .AND. vz>0 ) THEN
        vectors = .TRUE.
      ELSE
        vectors = .FALSE.
        vx=0
        vy=0
        vz=0
      ENDIF
    ENDIF
  ENDDO
ENDIF
WRITE(temp,*) 'vectors found: ', vectors
CALL ATOMSK_MSG(999,(/temp/),(/0.d0/))
!
!
!
100 CONTINUE
OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN',ERR=500)
!
!Write header of VESTA file
WRITE(40,'(a27)') "#VESTA_FORMAT_VERSION 3.0.0"
WRITE(40,*) ""
WRITE(40,'(a7)') "CRYSTAL"
WRITE(40,*) ""
WRITE(40,'(a5)') "TITLE"
WRITE(40,'(a)') TRIM(comment(1))
WRITE(40,*) ""
WRITE(40,'(a5)') "GROUP"
WRITE(40,'(a7)') "1 1 P 1"
WRITE(40,'(a5)') "SYMOP"
WRITE(40,'(a)') " 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1"
WRITE(40,'(a)') " -1.0 -1.0 -1.0  0 0 0  0 0 0  0 0 0"
WRITE(40,'(a5)') "TRANM"
WRITE(40,'(a)') " 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1"
WRITE(40,'(a7)') "LTRANSL"
WRITE(40,'(a)') " -1"
WRITE(40,'(a)') " 0.000000  0.000000  0.000000  0.000000  0.000000  0.000000"
WRITE(40,'(a7)') "LORIENT"
WRITE(40,'(a)') " -1   0   0   0   0"
WRITE(40,'(a)') " 1.000000  0.000000  0.000000  1.000000  0.000000  0.000000"
WRITE(40,'(a)') " 0.000000  0.000000  1.000000  0.000000  0.000000  1.000000"
WRITE(40,'(a7)') "LMATRIX"
WRITE(40,'(a)') " 1.000000  0.000000  0.000000  0.000000"
WRITE(40,'(a)') " 0.000000  1.000000  0.000000  0.000000"
WRITE(40,'(a)') " 0.000000  0.000000  1.000000  0.000000"
WRITE(40,'(a)') " 0.000000  0.000000  0.000000  1.000000"
WRITE(40,'(a)') " 0.000000  0.000000  0.000000"
!
!
!Write cell data
CALL MATCONV(H,a,b,c,alpha,beta,gamma)
WRITE(40,'(a5)') "CELLP"
WRITE(40,'(X,6(f12.6,3X))') a, b, c, alpha*180.d0/pi, beta*180.d0/pi, gamma*180.d0/pi
WRITE(40,'(a)') "  0.000000   0.000000   0.000000   0.000000   0.000000   0.000000"
!
!
!Write atom positions
a=1.d0
WRITE(40,'(a5)') "STRUC"
DO i=1,SIZE(P,1)
  CALL ATOMSPECIES(P(i,4),species)
  !
  IF(isreduced) THEN
    WRITE(temp,'(3(f12.8,2X))') P(i,1), P(i,2), P(i,3)
  ELSE
    P1 = P(i,1)
    P2 = P(i,2)
    P3 = P(i,3)
    WRITE(temp,'(3(f12.8,2X))')  P1*G(1,1) + P2*G(2,1) + P3*G(3,1),     &
                              &  P1*G(1,2) + P2*G(2,2) + P3*G(3,2),     &
                              &  P1*G(1,3) + P2*G(2,3) + P3*G(3,3)
  ENDIF
  !
  !Write line to the file
  IF( occ>0 ) THEN
    a = AUX(i,occ)
  ENDIF
  WRITE(40,110) i, species, species, a, TRIM(ADJUSTL(temp))//"    1a       1"
  WRITE(40,*) "                            0.000000   0.000000   0.000000  0.00"
ENDDO
WRITE(40,'(a)') "  0 0 0 0 0 0 0"
110 FORMAT(i6,2X,a2,2X,a2,2X,f12.8,2X,a)
!
!
IF( vectors ) THEN
  !Write the vector coordinates
  WRITE(40,'(a5)') "VECTR"
  DO i=1,SIZE(AUX,1)
    WRITE(40,'(i7,3(1X,f12.8),a3)') i, AUX(i,vx), AUX(i,vy), AUX(i,vz), "  0"
    WRITE(40,*) i, " 0    0    0    0"
    WRITE(40,*) " 0 0 0 0 0"
  ENDDO
  WRITE(40,'(a)') " 0 0 0 0 0"
  !Write vector properties
  WRITE(40,'(a5)') "VECTT"
  DO i=1,SIZE(AUX,1)
    WRITE(40,'(i7,a21)') i, "  0.500 255   0   0 0"
  ENDDO
  WRITE(40,'(a)') " 0 0 0 0 0"
ENDIF
!
!
!
500 CONTINUE
CLOSE(40)
msg = "VESTA"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
1000 CONTINUE
!
END SUBROUTINE WRITE_VESTA
!
END MODULE out_vesta
