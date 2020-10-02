MODULE oia_xsf
!
!**********************************************************************************
!*  1IA_XSF                                                                       *
!**********************************************************************************
!* This module reads files in the XSF format that contain several snapshots       *
!* (so-called "animated XSF" format).                                             *
!* The XSF format is described here:                                              *
!*    http://www.xcrysden.org/doc/XSF.html                                        *
!**********************************************************************************
!* (C) April 2011 - Pierre Hirel                                                  *
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
SUBROUTINE ONEINALL_XSF(inputfile,out_prefix,outfileformats,options_array)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to output
CHARACTER(LEN=7):: xsfformat
CHARACTER(LEN=128):: msg, test
CHARACTER(LEN=4096):: out_prefix, outputfile, outfileformat
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: atoms, isreduced
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: i, j, NP, Nsnap, snap
INTEGER:: Nsys  !number of systems converted
REAL(dp):: a0, testreal
REAL(dp):: snumber
REAL(dp),DIMENSION(3,3):: Huc   !Box vectors of unit cell (unknown, set to 0 here)
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystal orientation
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P,S  !positions of cores, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
outputfile=''
xsfformat = 'massxyz'
Nsys = 0
snap = 0
C_tensor(:,:) = 0.d0
ORIENT(:,:) = 0.d0
IF(ALLOCATED(P)) DEALLOCATE(P)
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
!
!
msg = 'entering ONEINALL_XSF'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
!Read number of snapshots
READ(30,'(a128)') test
test=ADJUSTL(test)
!If the first line is not 'ANIMSTEPS nstep' then the format is wrong
IF(test(1:9).NE.'ANIMSTEPS') THEN
  CALL ATOMSK_MSG(4812,(/''/),(/0.d0/))
  nerr=nerr+1
  GOTO 1000
ENDIF
READ(test(10:),*,ERR=800,END=800) Nsnap
!
!
!
100 CONTINUE
DO j=1,Nsnap
  !Initialize variables
  NP=0  
  snap = snap+1
  a0 = 1.d0
  !
  !
  !Set the output file name
  WRITE(msg,*) snap
  outputfile = TRIM(ADJUSTL(out_prefix))//TRIM(ADJUSTL(msg))
  WRITE(msg,*) 'snap, outputfile:', snap, TRIM(ADJUSTL(outputfile))
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !
  CALL ATOMSK_MSG(4041,(/''/),(/DBLE(snap)/))
  !
  !Read the file
  DO WHILE(NP==0)
    READ(30,'(a128)',ERR=400,END=400) test
    test = TRIM(ADJUSTL(test))
    IF(test(1:7)=='PRIMVEC') THEN
      READ(30,*,ERR=400,END=400) H(1,1), H(1,2), H(1,3)
      READ(30,*,ERR=400,END=400) H(2,1), H(2,2), H(2,3)
      READ(30,*,ERR=400,END=400) H(3,1), H(3,2), H(3,3)
    ELSEIF(test(1:9)=='PRIMCOORD') THEN
      READ(30,*,ERR=400,END=400) NP
    ELSEIF(test(1:5)=='ATOMS' ) THEN
      atoms = .TRUE.
      !count atoms
      DO
        READ(30,*,ERR=200,END=200) i
        NP = NP+1
      ENDDO
      !Go back to the first atom
      DO i=1,NP
        BACKSPACE(30)
      ENDDO
    ENDIF
  ENDDO
  !
  !
  200 CONTINUE
  WRITE(msg,*) 'NP:', NP
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !
  IF(.NOT.ALLOCATED(P)) THEN
    ALLOCATE(P(NP,4))
  ELSE
    IF( NP.NE.SIZE(P(:,1)) ) THEN
      GOTO 800
    ENDIF
  ENDIF
  P(:,:) = 0.d0
  !
  !We have to determine if the format is
  ! (number x y z)  or  (species x y z)
  READ(30,'(a128)',ERR=400,END=400) test
  test = ADJUSTL(test)
  READ(test,*,ERR=220,END=400) testreal
  !If it is a real, then the file has 'massxyz' format
  xsfformat = 'massxyz'
  GOTO 250
  220 CONTINUE
  !If it was not a real, then it has to be 'speciesxyz' format
  xsfformat = 'spiexyz'
  !Let's confirm that: check that the string does correspond to
  !a known species
  READ(test,*) species
  CALL ATOMNUMBER(species,snumber)
  WRITE(msg,*) 'spiexyz: species, snumber = ', species, snumber
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  IF(snumber==0.d0) THEN
    i = -1
    GOTO 400
  ENDIF
  !
  !
  250 CONTINUE
  BACKSPACE(30)
  !Read atomic positions
  DO i=1,NP
    IF(xsfformat=='massxyz') THEN
      READ(30,*,ERR=400,END=400) P(i,4), P(i,1), P(i,2), P(i,3)
    ELSE
      READ(30,*,ERR=400,END=400) species, P(i,1), P(i,2), P(i,3)
      CALL ATOMNUMBER(species,P(i,4))
    ENDIF
  ENDDO
  !
  CALL ATOMSK_MSG(4043,(/''/),(/0.d0/))
  !
  !
  !
  300 CONTINUE
  !Find out if coordinates are reduced or cartesian
  CALL FIND_IF_REDUCED(H,P,isreduced)
  !In case of reduced coordinates, convert them to cartesian
  IF(isreduced) THEN
    CALL FRAC2CART(P,H)
  ENDIF
  !
  !
  !
  400 CONTINUE
  !Apply options
  CALL OPTIONS_AFF(options_array,Huc,H,P,S,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
  !
  !Output snapshot to one or several format(s)
  outfileformat=''
  CALL WRITE_AFF(outputfile,outfileformats,H,P,S,comment,AUXNAMES,AUX)
  !
  !
  Nsys = Nsys+1
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
END SUBROUTINE ONEINALL_XSF
!
!
END MODULE oia_xsf
