MODULE out_cfg 
!
!**********************************************************************************
!*  OUT_CFG                                                                       *
!**********************************************************************************
!* This module writes Ju Li's extended CFG format.                                *
!* The CFG format is described here:                                              *
!*    http://mt.seas.upenn.edu/Archive/Graphics/A/#standard_CFG                   *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
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
SUBROUTINE WRITE_CFG(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=12):: test
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: isreduced
INTEGER:: i, j
REAL(dp):: smass, snumber
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
!
100 CONTINUE
msg = 'entering WRITE_CFG'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN',ERR=250)
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
! Write header of CFG file
CONTINUE
msg = 'writing header of CFG file: '//TRIM(outputfile)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
WRITE(temp,*) SIZE(P(:,1))
WRITE(40,'(a)') 'Number of particles = '//TRIM(ADJUSTL(temp))
DO i=1,SIZE(comment)
  WRITE(40,'(a)') TRIM(ADJUSTL(comment(i)))
ENDDO
WRITE(40,'(a45)') 'A = 1.000000000 Angstrom (basic length-scale)'
WRITE(40,102) 'H0(1,1) = ', H(1,1)
WRITE(40,102) 'H0(1,2) = ', H(1,2)
WRITE(40,102) 'H0(1,3) = ', H(1,3)
WRITE(40,102) 'H0(2,1) = ', H(2,1)
WRITE(40,102) 'H0(2,2) = ', H(2,2)
WRITE(40,102) 'H0(2,3) = ', H(2,3)
WRITE(40,102) 'H0(3,1) = ', H(3,1)
WRITE(40,102) 'H0(3,2) = ', H(3,2)
WRITE(40,102) 'H0(3,3) = ', H(3,3)
WRITE(40,'(a13)') '.NO_VELOCITY.'
102 FORMAT(a10,f16.8)
!Check if auxiliary properties are present
IF( ALLOCATED(AUX) .AND. SIZE(AUXNAMES)>0 ) THEN
  !Atomeye supports up to 32 auxiliary properties
  !If more exist, display a warning
  IF( SIZE(AUXNAMES)>32 ) THEN
    nwarn=nwarn+1
    CALL ATOMSK_MSG(3708,(/msg/),(/0.d0/))
  ENDIF
  !Write the number of entries for each line
  WRITE(40,'(a14,i3)') 'entry_count = ', MIN( SIZE(AUXNAMES(:)),32 )+3
  !Name of each auxiliary property
  DO i=1, MIN( SIZE(AUXNAMES(:)),32 )
    WRITE(temp,'(i3)') i-1
    WRITE(40,'(a)') 'auxiliary['//TRIM(ADJUSTL(temp))//'] = ' &
         &          //TRIM(ADJUSTL(AUXNAMES(i)))
  ENDDO
ELSE
  WRITE(40,'(a15)') 'entry_count = 3'
ENDIF
!
!
200 CONTINUE
! Write coordinates in CFG file
WRITE(msg,'(a32,i6)') 'write coordinates, NP=', SIZE(P,1)
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
test = ''
smass = 0.d0
snumber = 0.d0
species = ''
!
DO i=1,SIZE(P,1)
  !If it is a new species we have to write it in the CFG file
  IF(snumber.NE.P(i,4)) THEN
    !This atom is different from the previous one
    snumber = P(i,4)
    IF(snumber.NE.0.d0) THEN
      !Determine mass and species of this atom
      CALL ATOMSPECIES(snumber,species)
      CALL ATOMMASS(species,smass)
      !Atomic mass
      WRITE(40,'(f9.4)') smass
      !Atomic species
      WRITE(40,'(a2)') species
    ELSE
      !Set dummy mass and species
      WRITE(40,'(f9.4)') 0.d0
      WRITE(40,'(a2)') "XX"
    ENDIF
  ENDIF
  !
  !Coordinates
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
  !Append the auxiliray properties if any
  IF(ALLOCATED(AUX)) THEN
    DO j=1, MIN( SIZE(AUX(1,:)),32 )
      WRITE(msg,'(e12.5)')  AUX(i,j)
      temp = TRIM(ADJUSTL(temp))//'  '//TRIM(ADJUSTL(msg))
    ENDDO
  ENDIF
  !
  !Write the line to the file
  WRITE(40,'(a)') TRIM(ADJUSTL(temp))
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
msg = "CFG"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
1000 CONTINUE
CLOSE(40)
!
END SUBROUTINE WRITE_CFG
!
!
END MODULE out_cfg
