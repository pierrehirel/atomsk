MODULE in_cel
!
!
!**********************************************************************************
!*  IN_CEL                                                                        *
!**********************************************************************************
!* This module reads a Dr. Probe super-cell structure definition file (CEL).      *
!* The CEL format is used by the TEM simulation software Dr. Probe and described  *
!* on the software website:                                                       *
!*     http://www.er-c.org/barthel/drprobe/celfile.html                           *
!**********************************************************************************
!* (C) July 2015 - Juri Barthel                                                   *
!*     Gemeinschaftslabor fuer Elektronenmikroskopie                              *
!*     RWTH Aachen (GERMANY)                                                      *
!*     ju.barthel@fz-juelich.de                                                   *
!* Last modification: P. Hirel - 16 April 2024                                    *
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
SUBROUTINE READ_CEL(inputfile,H,P,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: isreduced
INTEGER:: i
INTEGER:: Naux !number of auxiliary properties
INTEGER:: Ncol, NP ! number of columns and number of atoms
INTEGER:: occ  !index of occupancies in AUXNAMES
INTEGER:: biso !index of Biso in AUXNAMES
REAL(dp):: a, b, c, alpha, beta, gamma !supercell (conventional notation)
REAL(dp):: sbiso
REAL(dp),DIMENSION(3,3),INTENT(OUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
Naux=0
NP=0
Ncol = 9
occ=0
biso=0
a=0.d0
b=0.d0
c=0.d0
alpha=0.d0
beta=0.d0
gamma=0.d0
!
msg = 'entering READ_CEL'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
!
100 CONTINUE
OPEN(UNIT=30,FILE=inputfile,STATUS='UNKNOWN',ERR=800)
!
!Go back to beginning of file and count number of lines
REWIND(30)
i=1  ! init line counter
!We try to read the whole file here, this is also a check
!for a corrupt input file.
!A proper CEL file should end with a '*' character in the
!last line.
DO
  READ(30,'(a128)',ERR=207,END=207) temp ! o.O ! Error or end of file -> corrupt cel file
  temp = ADJUSTL(temp)
  !CALL ATOMSK_MSG(999,(/temp/),(/0.d0/))
  !
  IF( temp(1:1)=='*' ) EXIT ! find the line which terminates the cel file
  !
  ! count number of lines before the end
  i = i+1
  !
ENDDO
!
NP = i-3 ! number of atom sites = number of lines - number of lines for header and termination
WRITE(msg,*) 'Counted atoms, NP = ', NP
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
IF (NP<=0) GOTO 800 ! No atoms, -> error
!
! allocate the atom data P and auxiliary data AUX and AUXNAMES
IF (ALLOCATED(P)) DEALLOCATE(P)
ALLOCATE(P(NP,4))
IF (ALLOCATED(AUX)) DEALLOCATE(AUX)
IF (ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
Naux = 2
WRITE(msg,*) 'Naux = ', Naux
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ALLOCATE(AUX(NP,Naux),AUXNAMES(Naux))
! init the new arrays
P = 0.d0
AUX = 0.d0
occ = 1
biso = 2
AUXNAMES(occ) = 'occ'
AUXNAMES(biso) = 'Debye-Waller'
!
!Go back to beginning of file and store data
REWIND(30)
! read the first header line
READ(30,'(a128)',ERR=200,END=200) temp ! this is a comment and can be ignored
! read the 2nd header line
READ(30,*,ERR=200,END=200) i, a, b, c, alpha, beta, gamma ! cell dimensions (nm by default in CEL files)
! convert the cell dimensions from nm to A
!a = 1.0d+1*a
!b = 1.0d+1*b
!c = 1.0d+1*c
! convert the cell angles from degree to radian
alpha = DEG2RAD(alpha)
beta  = DEG2RAD(beta)
gamma = DEG2RAD(gamma)
! read the atom site list
DO i=1, NP
  ! line format: <Symbol string> <f_x> <f_y> <f_z> <occ> <biso> <unused> <unused> <unused>
  ! using unformatted reading
  READ(30,*,ERR=200,END=200) temp, P(i,1), P(i,2), P(i,3), &
                           & AUX(i,occ), sbiso ! ignore the rest of the line
  ! interpret the symbol and set the atomic number in P(i,4)
  temp = ADJUSTL(temp)
  species = temp(1:2)
  CALL ATOMNUMBER(species,P(i,4))
  IF( P(i,4)<=0.d0 ) THEN
    species = temp(1:1)
    CALL ATOMNUMBER(species,P(i,4))
  ENDIF
ENDDO
!
!
200 CONTINUE
CLOSE(30)
!Save cell vectors in H(:,:)
CALL CONVMAT(a,b,c,alpha,beta,gamma,H)
!Find out if coordinates are reduced or cartesian
CALL FIND_IF_REDUCED(H,P,isreduced)
!In case of reduced coordinates, convert them to cartesian
IF(isreduced) THEN
  CALL FRAC2CART(P,H)
ENDIF
GOTO 1000
!
207 CONTINUE
CLOSE(30)
CALL ATOMSK_MSG(807,(/''/),(/DBLE(i)/))
nerr = nerr+1
GOTO 1000
!
!
800 CONTINUE
801 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE READ_CEL
!
END MODULE in_cel
