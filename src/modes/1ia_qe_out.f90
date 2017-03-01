MODULE oia_qeout
!
!**********************************************************************************
!*  1IA_QEOUT                                                                     *
!**********************************************************************************
!* This module reads Quantum Espresso PW output files (*.out), containing         *
!* several snapshots, and writes each snapshot to a separate file.                *
!* The output produced by Quantum Espresso is described here:                     *
!*    http://www.quantum-espresso.org/wp-content/uploads/Doc/pw_user_guide/       *
!* Note on units: by default in the Quantum Espresso PW file format,              *
!* all quantities whose dimensions are not explicitly specified are in            *
!* RYDBERG ATOMIC UNITS, in particular cell dimensions and atom positions         *
!* are in units of Bohr radius. If cell dimensions and atom coordinates are all   *
!* expressed in the same units (e.g. all in Bohr or all in anstroms), then        *
!* atomsk conserves this unit (i.e. coordinates are NOT converted to angstroms).  *
!* But, if no unit is specified for the cell (meaning it is in Bohr), while the   *
!* atom coordinates are specified in angstroms, then the cell dimensions          *
!* will be converted to angstroms for consistency.                                *
!**********************************************************************************
!* (C) September 2012 - Pierre Hirel                                              *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
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
!
USE atoms
USE comv
USE constants
USE messages
USE files
USE subroutines
USE options
USE writeout
USE qepw_ibrav
!
!
CONTAINS
!
SUBROUTINE ONEINALL_QEOUT(inputfile,out_prefix,outfileformats,options_array)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=1):: cell_units, atpos_units !units of cell vectors, atom positions: A (angstroms) or B (Bohrs)
CHARACTER(LEN=2):: species
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to output
CHARACTER(LEN=128):: msg, test
CHARACTER(LEN=4096):: out_prefix, outputfile, outfileformat
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: isreduced
LOGICAL:: cell_defined
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: ibrav  !index of Bravais lattice
INTEGER:: i, NP, Nsnap, snap, strlength, strlength2
INTEGER:: Nsys  !number of systems converted
REAL(dp):: alat  !lattice constant
REAL(dp),DIMENSION(6):: celldm
REAL(dp),DIMENSION(3,3):: H, Htemp   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystal orientation
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P,S  !positions of cores, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
outputfile=''
ibrav = 0
NP = 0
Nsnap = 0
Nsys = 0
snap = -1  !so that snap index will start at 0
alat = 0.d0
ORIENT(:,:) = 0.d0
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
IF(ALLOCATED(P)) DEALLOCATE(P)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
ALLOCATE(comment(1))
!
!
msg = 'entering ONEINALL_QEOUT'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
!
!
100 CONTINUE
DO  !loop on all snapshots
  !Initialize variables
  IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
  IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
  cell_defined = .FALSE.
  isreduced=.FALSE.
  atpos_units='B'
  cell_units='B'
  snap = snap+1
  !
  !
  !Set the output file name
  WRITE(msg,*) snap
  outputfile = TRIM(ADJUSTL(out_prefix))//"_"//TRIM(ADJUSTL(msg))
  comment(1) = "# Quantum Espresso PWSCF output snapshot # "//TRIM(ADJUSTL(msg))
  WRITE(msg,*) 'snap, outputfile:', snap, TRIM(ADJUSTL(outputfile))
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !
  !Find the number of atoms in each snapshot
  IF(NP==0) THEN
    DO WHILE(NP==0)
      READ(30,'(a128)',ERR=800,END=800) test
      test = TRIM(ADJUSTL(test))
      IF(test(1:15)=='number of atoms') THEN
        strlength = SCAN(test,'=')
        READ(test(strlength+1:),*,ERR=800,END=800) NP
      ENDIF
    ENDDO
    !
    WRITE(msg,*) 'NP=', NP
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
    ALLOCATE(P(NP,4))
    !
    !Other information (like alat) may appear before NP
    !=> go back to start of file
    REWIND(30)
    !
  ENDIF
  P(:,:) = 0.d0
  !
  !
  200 CONTINUE
  !Read file until cell vectors or atom positions are found
  DO
    READ(30,'(a128)',ERR=500,END=500) test
    test = TRIM(ADJUSTL(test))
    !
    IF( test(1:21)=="bravais-lattice index" ) THEN
      !Read the index of Bravais lattice ibrav
      strlength = SCAN(test,"=")
      READ(test(strlength+1:),*,ERR=800,END=800) ibrav
      !
    ELSEIF( test(1:17)=="lattice parameter" ) THEN
      !Here is the lattice constant, with the format:
      !  lattice parameter (alat)  =       6.6604  a.u.
      strlength = SCAN(test,"=")
      test = TRIM(ADJUSTL(test(strlength+1:)))
      READ(test,*,ERR=800,END=800) alat
      !Try to read units, go on if it fails
      READ(test,*,ERR=210,END=210) alat, msg
      210 CONTINUE
      msg=ADJUSTL(msg)
      IF( msg(1:4)=="a.u." ) THEN
        cell_units = 'B'
      ENDIF
      cell_defined = .TRUE.
      !
    ELSEIF( test(1:12)=="crystal axes" ) THEN
      !After this line come the three cell vectors, with the format:
      !  a(1) = (   1.000000   0.000000   0.000000 )
      !  a(2) = (   0.000000   1.024404   0.000000 )
      !  a(3) = (   0.000000   0.000000   1.445291 )
      DO i=1,3
        READ(30,'(a128)',ERR=800,END=800) test
        test = TRIM(ADJUSTL(test))
        strlength = SCAN(test,"=")
        test = test(strlength+1:)
        strlength = SCAN(test,"(")
        strlength2 = SCAN(test,")")
        READ(test(strlength+1:strlength2-1),*,ERR=800,END=800) H(i,1), H(i,2), H(i,3)
      ENDDO
      IF(alat>0.d0) THEN
        H(:,:) = alat*H(:,:)
      ENDIF
      cell_defined = .TRUE.
      !
    ELSEIF( test(1:9)=="celldm(1)" ) THEN
      !Here are the cell parameters, with the format:
      !  celldm(1)=   8.783574  celldm(2)=   0.000000  celldm(3)=   0.000000
      !  celldm(4)=   0.000000  celldm(5)=   0.000000  celldm(6)=   0.000000
      DO i=1,3
        strlength = SCAN(test,"=")
        test = test(strlength+1:)
        strlength2 = SCAN(test,"c")
        IF(strlength2==0) strlength2 = LEN(test)+1
        READ(test(1:strlength2-1),*,ERR=800,END=800) celldm(i)
      ENDDO
      READ(30,'(a128)',ERR=500,END=500) test
      test = TRIM(ADJUSTL(test))
      DO i=1,3
        strlength = SCAN(test,"=")
        test = test(strlength+1:)
        strlength2 = SCAN(test,"c")
        IF(strlength2==0) strlength2 = LEN(test)
        READ(test(1:strlength2-1),*,ERR=800,END=800) celldm(3+i)
      ENDDO
      !
      IF(alat<=1.d-12 .AND. celldm(1)>0.d0) THEN
        alat = celldm(1)
      ENDIF
      IF(ibrav.NE.0) THEN
        !Define the cell vectors H(:,:) according to the value of ibrav and the celldm(:)
        CALL QE_PW_IBRAV(ibrav,celldm(:),H)
      ENDIF
      !
      !No unit was specified => H(:,:) is in Bohrs
      cell_units = 'B'
      cell_defined = .TRUE.
      !
    ELSEIF( test(1:15)=="CELL_PARAMETERS" ) THEN
      !In case of variable cell simulations ("vc-relax" or "vc-md")
      !the cell parameters are given for each snapshot
      IF( INDEX(test(16:),"alat").NE.0 ) THEN
        !Read lattice constant alat
        strlength = SCAN(test,"=")
        strlength2 = SCAN(test,")")
        READ(test(strlength+1:strlength2-1),*,ERR=220,END=220) alat
        !No unit was specified => it is in Bohrs
        cell_units = 'B'
      ENDIF
      !
      220 CONTINUE
      WRITE(msg,*) 'alat=', alat
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      !Read supercell vectors
      DO i=1,3
        READ(30,*,ERR=800,END=800) H(i,1), H(i,2), H(i,3)
      ENDDO
      H(:,:) = alat*H(:,:)
      cell_defined = .TRUE.
      !
    ELSEIF( test(1:16)=="ATOMIC_POSITIONS" ) THEN
      !Atom positions follow, with the format:
      !  H       4.5000000   1.5000000   6.5000000
      CALL ATOMSK_MSG(4041,(/''/),(/DBLE(snap)/))
      !
      !Check if positions are in reduced coordinates
      test = TRIM(ADJUSTL(test(17:)))
      IF( INDEX(test,"crystal")>0 ) THEN
        isreduced = .TRUE.
      ELSEIF( INDEX(test,"angst")>0 ) THEN
        atpos_units = "A"
      ELSEIF( INDEX(test,"bohr")>0 ) THEN
        atpos_units = "B"
      ELSEIF( LEN_TRIM(test)==0 ) THEN
        atpos_units = "z"
      ENDIF
      !
      !Read atomic positions
      DO i=1,SIZE(P,1)
        READ(30,*,ERR=800,END=800) species, P(i,1), P(i,2), P(i,3)
        CALL ATOMNUMBER(species,P(i,4))
      ENDDO
      !
      !Finished reading this snapshot: exit the loop
      Nsnap = Nsnap+1
      GOTO 300
      !
    ELSEIF( test(1:22)=="Forces acting on atoms" ) THEN
      !There should be one empty line
      READ(30,*,ERR=800,END=800) test
      !If the line was not empty, go back
      IF( LEN_TRIM(test).NE.0 ) THEN
        BACKSPACE(30)
      ENDIF
      !
      !Allocate arrays
      ALLOCATE( AUXNAMES(3) )
      AUXNAMES(1) = "fx"
      AUXNAMES(2) = "fy"
      AUXNAMES(3) = "fz"
      ALLOCATE( AUX(SIZE(P,1),3) )
      AUX(:,:) = 0.d0
      !
      !Read the forces
      DO i=1,SIZE(AUX,1)
        READ(30,'(a128)',ERR=800,END=800) test
        strlength = SCAN(test,"=")
        READ(test(strlength+1:),*,ERR=800,END=800) AUX(i,1), AUX(i,2), AUX(i,3)
      ENDDO
      !
    ENDIF
    !
  ENDDO
  !
  !
  !
  300 CONTINUE
  WRITE(msg,*) "Units of cell, atom positions: ", cell_units, atpos_units
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !
  IF( cell_defined ) THEN
    !A new cell was defined
    IF(cell_units=="B" .AND. atpos_units=="A") THEN
      !Atom coordinates are in angstroms while cell vectors are in Bohrs
      !=> convert cell into angströms for consistency
      nwarn=nwarn+1
      CALL ATOMSK_MSG(1706,(/msg/),(/0.d0/))
      H(:,:) = 1.d10*a_bohr * H(:,:)
      
      ! => convert atom positions into Bohrs for consistency
!       nwarn=nwarn+1
!       CALL ATOMSK_MSG(1706,(/msg/),(/0.d0/))
!       DO i=1,SIZE(P,1)
!         P(i,1:3) = P(i,1:3) / (1.d10*a_bohr)
!       ENDDO
    ENDIF
  ENDIF
  !
  !Convert positions to cartesian coordinates
  IF(isreduced) THEN
    !Positions are in reduced coordinates
    CALL FRAC2CART(P,H)
  ELSEIF(atpos_units=="z") THEN
    !Positions are in units of alat=celldm(1)
    P(:,1:3) = alat*P(:,1:3)
  ENDIF
  CALL ATOMSK_MSG(4043,(/''/),(/0.d0/))
  !
  !
  !
  400 CONTINUE
  !H may be modified by options
  !=> pass a copy Htemp instead, to prevent H from being modified at each passage
  Htemp = H
  !
  !Apply options
  CALL OPTIONS_AFF(options_array,Htemp,P,S,AUXNAMES,AUX,ORIENT,SELECT)
  !
  !Output snapshot to one or several format(s)
  outfileformat=''
  CALL WRITE_AFF(outputfile,outfileformats,Htemp,P,S,comment,AUXNAMES,AUX)
  !
  !
  Nsys = Nsys+1
  CALL ATOMSK_MSG(4001,(/''/),(/0.d0/))
ENDDO
!
!
!
500 CONTINUE
CALL ATOMSK_MSG(4042,(/''/),(/DBLE(Nsnap)/))
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
END SUBROUTINE ONEINALL_QEOUT
!
!
END MODULE oia_qeout
