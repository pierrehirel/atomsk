MODULE out_qe_pw
!
!
!**********************************************************************************
!*  OUT_QE_PW                                                                     *
!**********************************************************************************
!* This module writes pw input files for Quantum Espresso.                        *
!* This file format is described here:                                            *
!*    http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html        *
!**********************************************************************************
!* (C) June 2012 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 04 Dec. 2015                                     *
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
SUBROUTINE WRITE_QEPW(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: pseudo_dir
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: fileexists !does file exist?
LOGICAL:: isreduced  !are positions in reduced coordinates?
LOGICAL:: pseudo_dir_exists !does the pseudo_dir exist?
INTEGER:: i, j
INTEGER:: fx, fy, fz       !position of forces (x,y,z) in AUX
INTEGER:: fixx, fixy, fixz !position of flags for fixed atoms in AUX
REAL(dp):: smass  !mass of atoms
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
pseudo_dir = ""
pseudo_dir_exists = .FALSE.
fx=0
fy=0
fz=0
fixx=0
fixy=0
fixz=0
!
WRITE(msg,*) 'entering WRITE_QEPW'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Find number of species
CALL FIND_NSP(P(:,4),aentries)
!
!Check if coordinates are reduced
CALL FIND_IF_REDUCED(P,isreduced)
WRITE(msg,*) 'isreduced:', isreduced
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Check if some atoms are fixed
IF( ALLOCATED(AUXNAMES) .AND. SIZE(AUXNAMES)>0 ) THEN
  DO i=1,SIZE(AUXNAMES)
    IF(TRIM(ADJUSTL(AUXNAMES(i)))=="fx") THEN
      fx=i
    ELSEIF(TRIM(ADJUSTL(AUXNAMES(i)))=="fy") THEN
      fy=i
    ELSEIF(TRIM(ADJUSTL(AUXNAMES(i)))=="fz") THEN
      fz=i
    ELSEIF(TRIM(ADJUSTL(AUXNAMES(i)))=="fixx") THEN
      fixx=i
    ELSEIF(TRIM(ADJUSTL(AUXNAMES(i)))=="fixy") THEN
      fixy=i
    ELSEIF(TRIM(ADJUSTL(AUXNAMES(i)))=="fixz") THEN
      fixz=i
    ENDIF
  ENDDO
ENDIF
!
!
!
100 CONTINUE
OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN',ERR=500)
!
!Write control section
WRITE(40,'(a8)') "&CONTROL"
WRITE(40,'(a)') "  title = '"//TRIM(comment(1))//"'"
!pseudo_dir: value of the $ESPRESSO_PSEUDO environment variable if set;
!            '$HOME/espresso/pseudo/' otherwise
CALL GET_ENVIRONMENT_VARIABLE('ESPRESSO_PSEUDO',pseudo_dir)
IF( LEN_TRIM(pseudo_dir)<=0 ) THEN
  CALL GET_ENVIRONMENT_VARIABLE('HOME',msg)
  pseudo_dir = TRIM(ADJUSTL(msg))//"/espresso/pseudo/"
ENDIF
!Verify if the pseudo_dir actually exists
temp = TRIM(ADJUSTL(pseudo_dir))//".atomsk.tmp"
OPEN(UNIT=49,FILE=temp,FORM="FORMATTED",STATUS="UNKNOWN",ERR=110,IOSTAT=j)
110 CONTINUE
CLOSE(49,STATUS='DELETE')
IF( j==0 ) THEN
  !Directory exists
  pseudo_dir_exists = .TRUE.
ELSE
  !There was an error when trying to open a file in that directory
  !=> directory does not exist
  pseudo_dir_exists = .FALSE.
ENDIF
!
WRITE(40,'(a)') "  pseudo_dir = '"//TRIM(ADJUSTL(pseudo_dir))//"'"
WRITE(40,'(a)') "  calculation = 'scf'"
WRITE(40,'(a1)') "/"
!
!Write system section
WRITE(40,*) ""
WRITE(40,'(a7)') "&SYSTEM"
WRITE(msg,*) SIZE(P,1)
WRITE(40,'(a)') "  nat= "//TRIM(ADJUSTL(msg))
WRITE(msg,*) SIZE(aentries,1)
WRITE(40,'(a)') "  ntyp= "//TRIM(ADJUSTL(msg))
WRITE(40,'(a)') "  ibrav= 0"
WRITE(40,'(a)') "  ecutwfc= 20.0"
WRITE(40,'(a1)') "/"
!
!Write electrons section
WRITE(40,*) ""
WRITE(40,'(a10)') "&ELECTRONS"
WRITE(40,'(a19)') "  mixing_beta = 0.7"
WRITE(40,'(a20)') "  conv_thr =  1.0d-8"
WRITE(40,'(a1)') "/"
!
!Write empty "ions" and "cell" sections
WRITE(40,*) ""
WRITE(40,'(a5)') "&IONS"
WRITE(40,'(a1)') "/"
WRITE(40,*) ""
WRITE(40,'(a5)') "&CELL"
WRITE(40,'(a1)') "/"
!
!Write mass of species
!Note: the user will have to append the pseudopotential file names
WRITE(40,*) ""
WRITE(40,'(a14)') "ATOMIC_SPECIES"
DO i=1,SIZE(aentries,1)
  fileexists = .FALSE.
  CALL ATOMSPECIES(aentries(i,1),species)
  !
  IF( pseudo_dir_exists ) THEN
    !Look for files starting with element name in pseudo_dir
    msg = TRIM(ADJUSTL(pseudo_dir))//TRIM(ADJUSTL(species))//".*"
    CALL SYSTEM(system_ls//" "//TRIM(msg)//" > .atomsk.tmp.out_qe_pw")
    INQUIRE(FILE=".atomsk.tmp.out_qe_pw",EXIST=fileexists)
    msg = ""
    temp = ""
    IF( fileexists ) THEN
      OPEN(UNIT=50,FILE=".atomsk.tmp.out_qe_pw",FORM="FORMATTED",STATUS="UNKNOWN")
      READ(50,'(a4096)',ERR=150,END=150) temp
      150 CONTINUE
      CLOSE(50,STATUS='DELETE')
    ENDIF
    !
    IF( LEN_TRIM(temp)>0 ) THEN
      !Verify that file exists
      INQUIRE(FILE=temp,EXIST=fileexists)
      IF( fileexists ) THEN
        j=SCAN(temp,"/",BACK=.TRUE.)
        msg = TRIM(ADJUSTL(temp(j+1:)))
      ENDIF
    ELSE
      fileexists = .FALSE.
    ENDIF
  ENDIF
  !
  !If no suitable file was found, use a dummy one
  IF( .NOT.fileexists ) THEN
    msg = TRIM(species)//".fixme.upf"
  ENDIF
  !
  CALL ATOMMASS(species,smass)
  WRITE(40,'(a2,2X,f6.3,2X,a)') species, smass, TRIM(msg)
ENDDO
!
!Write cell parameters
WRITE(40,*) ""
WRITE(40,'(a24)') "CELL_PARAMETERS angstrom"
WRITE(40,201) H(1,1), H(1,2), H(1,3)
WRITE(40,201) H(2,1), H(2,2), H(2,3)
WRITE(40,201) H(3,1), H(3,2), H(3,3)
201 FORMAT(3(f16.8,2X))
!
!Write atom coordinates
WRITE(40,*) ""
msg = "ATOMIC_POSITIONS"
IF( isreduced ) THEN
  msg = TRIM(ADJUSTL(msg))//" crystal"
ELSE
  msg = TRIM(ADJUSTL(msg))//" angstrom"
ENDIF
WRITE(40,'(a)') TRIM(msg)
DO i=1,SIZE(P,1)
  CALL ATOMSPECIES(P(i,4),species)
  WRITE(msg,210) species, P(i,1), P(i,2), P(i,3)
  IF( fixx>0 .OR. fixy>0 .OR. fixz>0 ) THEN
    !Caution: internally if AUX(fix)==1 then atom is fixed, but
    !  in Quantum Espresso the forces on atoms are multiplied by these numbers,
    !  as a result flag "0" means that atom is fixed, "1" that it's mobile
    !Note: even if only one fix is defined (e.g. fixx>0 but fixy=fixz=0), the
    !  three flags must appear, so in undefined directions just write a 1.
    IF( fixx>0 ) THEN
      IF( AUX(i,fixx)>0.5d0 ) THEN
        msg = TRIM(msg)//" 0"
      ELSE
        msg = TRIM(msg)//" 1"
      ENDIF
    ELSE
      msg = TRIM(msg)//" 1"
    ENDIF
    IF( fixy>0 ) THEN
      IF( AUX(i,fixy)>0.5d0 ) THEN
        msg = TRIM(msg)//" 0"
      ELSE
        msg = TRIM(msg)//" 1"
      ENDIF
    ELSE
      msg = TRIM(msg)//" 1"
    ENDIF
    IF( fixz>0 ) THEN
      IF( AUX(i,fixz)>0.5d0 ) THEN
        msg = TRIM(msg)//" 0"
      ELSE
        msg = TRIM(msg)//" 1"
      ENDIF
    ELSE
      msg = TRIM(msg)//" 1"
    ENDIF
  ENDIF
  WRITE(40,'(a)') TRIM(msg)
ENDDO
210 FORMAT(a2,2X,3(f16.8,2X))
!
!Write k points section
WRITE(40,*) ""
WRITE(40,'(a18)') "K_POINTS automatic"
WRITE(40,'(a)') "2 2 2  0 0 0"
!
!Write forces on atoms
IF( fx>0 .AND. fy>0 .AND. fz>0 ) THEN
  WRITE(40,*) ""
  msg = "ATOMIC_FORCES"
  WRITE(40,'(a)') TRIM(msg)
  DO i=1,SIZE(P,1)
    CALL ATOMSPECIES(P(i,4),species)
    WRITE(msg,210) species, AUX(i,fx), AUX(i,fy), AUX(i,fz)
    WRITE(40,'(a)') TRIM(msg)
  ENDDO
ENDIF
GOTO 500
!
!
!
500 CONTINUE
CLOSE(40)
msg = "PW"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
1000 CONTINUE
!
!
!
END SUBROUTINE WRITE_QEPW
!
END MODULE out_qe_pw
