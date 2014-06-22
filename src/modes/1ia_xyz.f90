MODULE oia_xyz
!
!**********************************************************************************
!*  1IA_XYZ                                                                       *
!**********************************************************************************
!* This module reads files in the XYZ format that contain several snapshots.      *
!* No standard specification exists, but a widely used format is                  *
!* described for instance at:                                                     *
!*     http://openbabel.org/wiki/XYZ_%28format%29                                 *
!* This module can also read extended XYZ format, where the comment line          *
!* (i.e. the 2nd line) is replaced by a set of keywords/values. The               *
!* extended XYZ format is described for instance at:                              *
!*     http://www.jrkermode.co.uk/quippy/extended_xyz.html                        *
!* Note that this module can NOT read special XYZ format.                         *
!**********************************************************************************
!* (C) April 2011 - Pierre Hirel                                                  *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 04 June 2014                                     *
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
USE atoms
USE comv
USE constants
USE messages
USE files
USE subroutines
USE options
USE writeout
!
!
CONTAINS
!
SUBROUTINE ONEINALL_XYZ(inputfile,out_prefix,outfileformats,options_array)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to output
CHARACTER(LEN=7):: xyzformat
CHARACTER(LEN=128):: msg, temp
CHARACTER(LEN=4096):: out_prefix, outputfile, outfileformat
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: isreduced, Hset
INTEGER:: i, NP, strlength
INTEGER:: Nsys  !number of systems converted
INTEGER:: snap  !snapshot number
REAL(dp):: a0, testreal
REAL(dp):: snumber
REAL(dp):: P1, P2, P3
REAL(dp), DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystal orientation
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P,S  !positions of cores, shells
REAL(dp), DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
outputfile=''
Nsys = 0
snap = 0
ORIENT(:,:) = 0.d0
ALLOCATE(comment(1))
 comment(1)=''
!
!
msg = 'entering ONEINALL_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
!
!
100 CONTINUE
DO
  !Initialize variables
  isreduced = .FALSE.
  Hset=.FALSE.
  H(:,:) = 0.d0
  i=0
  snap = snap+1
  a0 = 1.d0
  xyzformat = 'massxyz'
  IF(ALLOCATED(P)) DEALLOCATE(P)
  !
  !
  !Set the output file name
  WRITE(msg,*) snap
  outputfile = TRIM(ADJUSTL(out_prefix))//TRIM(ADJUSTL(msg))
  WRITE(msg,*) 'snap, outputfile:', snap, TRIM(ADJUSTL(outputfile))
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  CALL ATOMSK_MSG(4041,(/''/),(/DBLE(snap)/))
  !
  !Read the file (this is copied from in_xyz.f90)
  READ(30,*,ERR=800,END=500) NP
  IF(NP==0) GOTO 800
  ALLOCATE(P(NP,4))
  READ(30,'(a128)',ERR=800,END=800) comment(1)   !Comment
  !Check if the comment contains the supercell parameters
  DO i=1,120
    temp = ADJUSTL(comment(1)(i:128))
    IF(temp(1:9)=='Lattice="') THEN
      strlength = SCAN(temp(10:128),'"')+8
      READ(temp(10:strlength),*,ERR=250) H(1,1), H(2,1), H(3,1), &
	  & H(1,2), H(2,2), H(3,2), H(1,3), H(2,3), H(3,3)
      Hset=.TRUE.
      comment(1)=''
    !ELSEIF(temp(1:11)=='Properties=') THEN
    !  READ(temp,*)
    ENDIF
  ENDDO
  !
  !
  210 CONTINUE
  !Determine if the format is (number x y z)  or  (species x y z)
  msg = 'determine if format is massxyz or spiexyz'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  READ(30,*,ERR=800,END=800) temp
  READ(temp,*,ERR=211,END=800) testreal
  !If it is a real, then the file has 'massxyz' format
  xyzformat = 'massxyz'
  GOTO 250
  211 CONTINUE
  !If it was not a real, then it has to be 'speciesxyz' format
  xyzformat = 'spiexyz'
  !Let's confirm that: check that the string does correspond to
  !a known species
  READ(temp,*) species
  CALL ATOMNUMBER(species,snumber)
  WRITE(msg,*) 'spiexyz: species, snumber = ', species, snumber
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  IF(snumber==0.d0) THEN
    i = -1
    GOTO 800
  ENDIF
  !
  !
  250 CONTINUE
  !Go back one line (because the first line was already read above)
  BACKSPACE(UNIT=30)
  !Read atomic positions
  WRITE(msg,'(a36,i7)') 'read atomic positions; NP=', NP
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,NP
    IF(xyzformat=='massxyz') THEN
      READ(30,*,ERR=800,END=800) P(i,4), P1, P2, P3
    ELSEIF(xyzformat=='spiexyz') THEN
      READ(30,*,ERR=800,END=800) temp, P1, P2, P3
      READ(temp,*) species
      CALL ATOMNUMBER(species,P(i,4))
    ENDIF
    P(i,1) = P1
    P(i,2) = P2
    P(i,3) = P3
  ENDDO
  !
  CALL ATOMSK_MSG(4043,(/''/),(/0.d0/))
  !
  !
  !
  300 CONTINUE
  !Find out if coordinates are reduced or cartesian
  CALL FIND_IF_REDUCED(P,isreduced)
  !In case of reduced coordinates, convert them to cartesian
  IF(isreduced .AND. Hset) THEN
    CALL FRAC2CART(P,H)
  ENDIF
  !
  !
  !
  400 CONTINUE
  !Apply options
  CALL OPTIONS_AFF(options_array,H,P,S,AUXNAMES,AUX,ORIENT)
  !
  !Output snapshot to one or several format(s)
  outfileformat=''
  CALL WRITE_AFF(outputfile,outfileformats,H,P,S,comment,AUXNAMES,AUX)
  !
  !
  Nsys=Nsys+1
  CALL ATOMSK_MSG(4001,(/''/),(/0.d0/))
ENDDO
!
!
!
500 CONTINUE
CALL ATOMSK_MSG(4042,(/''/),(/DBLE(snap)/))
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(1801,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
1000 CONTINUE
CLOSE(30)
CALL ATOMSK_MSG(4045,(/''/),(/DBLE(Nsys)/))
!
END SUBROUTINE ONEINALL_XYZ
!
!
END MODULE oia_xyz
