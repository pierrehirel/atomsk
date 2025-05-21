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
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 20 May 2025                                      *
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
CHARACTER(LEN=4096):: filename, msg, temp
CHARACTER(LEN=4096):: tmpfile
CHARACTER(LEN=4096),DIMENSION(10):: ppfiles
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: cubic
LOGICAL:: fileexists !does file exist?
LOGICAL:: isreduced  !are positions in reduced coordinates?
LOGICAL:: pseudo_dir_exists !does the pseudo_dir exist?
INTEGER:: i, j, k
INTEGER:: fx, fy, fz       !position of forces (x,y,z) in AUX
INTEGER:: fixx, fixy, fixz !position of flags for fixed atoms in AUX
REAL(dp):: ecutwfc, ecutrho !minimum cutoff read from pseudopotential files
REAL(dp):: smass  !mass of atoms
REAL(dp):: tempreal
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
pseudo_dir = ""
ppfiles(:) = ""
tmpfile = ".atomsk.tmp.out_qe_pw"
cubic = .FALSE.
pseudo_dir_exists = .FALSE.
fx=0
fy=0
fz=0
fixx=0
fixy=0
fixz=0
!Default values for ecutwfc and ecutrho
!If PP files are found, values will be read from them
ecutwfc = 25.d0
ecutrho = 100.d0
!
WRITE(msg,*) 'entering WRITE_QEPW'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Determine if cell is cubic
IF( DABS(H(1,1)-H(2,2))<1.d-12 .AND. DABS(H(1,1)-H(3,3))<1.d-12 .AND. &
  & ORTHOVEC(H(1,:),H(2,:)) .AND. ORTHOVEC(H(2,:),H(3,:))             ) THEN
  cubic = .TRUE.
ENDIF
!
!Find number of species
CALL FIND_NSP(P(:,4),aentries)
!
!Check if coordinates are reduced
CALL FIND_IF_REDUCED(H,P,isreduced)
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
    ELSEIF(TRIM(ADJUSTL(AUXNAMES(i)))=="freeze_x") THEN
      fixx=i
    ELSEIF(TRIM(ADJUSTL(AUXNAMES(i)))=="freeze_y") THEN
      fixy=i
    ELSEIF(TRIM(ADJUSTL(AUXNAMES(i)))=="freeze_z") THEN
      fixz=i
    ENDIF
  ENDDO
ENDIF
!
!pseudo_dir: value of the $ESPRESSO_PSEUDO environment variable if set;
!            '$HOME/espresso/pseudo/' otherwise
CALL GET_ENVIRONMENT_VARIABLE('ESPRESSO_PSEUDO',pseudo_dir)
IF( LEN_TRIM(pseudo_dir)<=0 ) THEN
  CALL GET_ENVIRONMENT_VARIABLE('HOME',msg)
  pseudo_dir = TRIM(ADJUSTL(msg))//"/espresso/pseudo/"
ENDIF
!Last character should be a slash
j=LEN_TRIM(pseudo_dir)
IF( pseudo_dir(j:j) .NE. "/" ) THEN
  pseudo_dir = TRIM(pseudo_dir)//"/"
ENDIF
!Verify if the pseudo_dir actually exists
temp = TRIM(ADJUSTL(pseudo_dir))//".atomsk.tmp"
OPEN(UNIT=49,FILE=temp,FORM="FORMATTED",STATUS="UNKNOWN",ERR=60,IOSTAT=j)
60 CONTINUE
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
DO i=1,SIZE(aentries,1)
  fileexists = .FALSE.
  CALL ATOMSPECIES(aentries(i,1),species)
  !
  IF( pseudo_dir_exists ) THEN
    !Look for files starting with element name in pseudo_dir
    !Save file name pattern into msg, something like "/home/user/espresso/Al*.*"
    msg = TRIM(ADJUSTL(pseudo_dir))//TRIM(ADJUSTL(species))//".*"
    !List all files matching this pattern, and save the list in temporary file tmpfile
    !This will execute something like "ls /home/user/espresso/Al*.* >.atomsk.tmp.out_qe_pw > /dev/null 2>&1"
    CALL SYSTEM(system_ls//" "//TRIM(msg)//" > "//tmpfile//" "//pathnull)
    !Verify if tmpfile exists
    INQUIRE(FILE=tmpfile,EXIST=fileexists)
    msg = ""
    temp = ""
    !If it does, read the first file name from it, save it to filename
    IF( fileexists ) THEN
      OPEN(UNIT=50,FILE=tmpfile,FORM="FORMATTED",STATUS="UNKNOWN")
      READ(50,'(a4096)',ERR=70,END=70) filename
      70 CONTINUE
      CLOSE(50,STATUS='DELETE')
    ENDIF
    !
    IF( LEN_TRIM(filename)>0 ) THEN
      !Verify that this file exists
      INQUIRE(FILE=filename,EXIST=fileexists)
      IF( fileexists ) THEN
        !This pseudopotential file exists: save its name in ppfiles()
        j=SCAN(filename,"/",BACK=.TRUE.)
        ppfiles(i) = TRIM(ADJUSTL(filename(j+1:)))
        !Read min. values of ecutwcf and ecutrho from PP file
        OPEN(UNIT=51,FILE=filename,FORM="FORMATTED",STATUS="UNKNOWN")
        j=0
        DO WHILE(j<2)
          READ(51,'(a)',END=90,ERR=90) msg
          msg = TRIM(ADJUSTL(msg))
          k = INDEX(msg,"minimum cutoff for wavefunctions")
          IF( k>0 ) THEN
            k = SCAN(msg,':',BACK=.TRUE.)
            msg = ADJUSTL(msg(k+1:))
            READ(msg,*,ERR=90,END=90) tempreal
            IF( tempreal>ecutwfc ) THEN
              ecutwfc = tempreal
              ecutrho = 4.d0*ecutwfc
            ENDIF
            j=j+1
          ENDIF
          k = INDEX(msg,"minimum cutoff for charge density")
          IF( k>0 ) THEN
            k = SCAN(msg,':',BACK=.TRUE.)
            msg = ADJUSTL(msg(k+1:))
            READ(msg,*,ERR=90,END=90) tempreal
            IF( tempreal>ecutrho ) THEN
              ecutrho = tempreal
            ENDIF
            j=j+1
          ENDIF
        ENDDO
        90 CONTINUE
        CLOSE(51)
      ENDIF
    ELSE
      fileexists = .FALSE.
    ENDIF
  ENDIF
  !
  !If no suitable file was found, write a dummy file name into the PW file
  IF( .NOT.fileexists ) THEN
    ppfiles(i) = TRIM(species)//".fixme.upf"
  ENDIF
ENDDO
!
!
!
100 CONTINUE
IF(ofu.NE.6) THEN
  OPEN(UNIT=ofu,FILE=outputfile,STATUS='UNKNOWN',ERR=500)
ENDIF
!
!Write control section
WRITE(ofu,'(a8)') "&CONTROL"
IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
  WRITE(ofu,'(a)') "  title = '"//TRIM(comment(1))//"'"
ELSE
  WRITE(ofu,'(a)') "  title = ''"
ENDIF
!
WRITE(ofu,'(a)') "  pseudo_dir = '"//TRIM(ADJUSTL(pseudo_dir))//"'"
WRITE(ofu,'(a)') "  calculation = 'scf'"
WRITE(ofu,'(a)') "  restart_mode = 'from_scratch'"
WRITE(ofu,'(a1)') "/"
!
!Write system section
WRITE(ofu,*) ""
WRITE(ofu,'(a7)') "&SYSTEM"
WRITE(msg,*) SIZE(P,1)
WRITE(ofu,'(a)') "  nat = "//TRIM(ADJUSTL(msg))
WRITE(msg,*) SIZE(aentries,1)
WRITE(ofu,'(a)') "  ntyp = "//TRIM(ADJUSTL(msg))
IF( cubic ) THEN
  WRITE(ofu,'(a)') "  ibrav = 1"
  !celldm() must be written in units of Bohr radius (atomic units)
  !Here we convert H(1,1) assuming it is in angströms
  WRITE(temp,'(f15.6)') H(1,1) / (a_bohr*1.d10)
  WRITE(ofu,'(a)') "  celldm(1) = "//TRIM(ADJUSTL(temp))
ELSE
  WRITE(ofu,'(a)') "  ibrav = 0"
ENDIF
!Values of ecutwfc and ecutrho were read from PP files if they were found (see above)
!Otherwise default values are used
WRITE(temp,'(i5)') CEILING(ecutwfc)
WRITE(ofu,'(a)') "  ecutwfc = "//TRIM(ADJUSTL(temp))
WRITE(temp,'(i5)') CEILING(ecutrho)
WRITE(ofu,'(a)') "  ecutrho = "//TRIM(ADJUSTL(temp))
WRITE(ofu,'(a1)') "/"
!
!Write electrons section
WRITE(ofu,*) ""
WRITE(ofu,'(a10)') "&ELECTRONS"
WRITE(ofu,'(a19)') "  mixing_beta = 0.7"
WRITE(ofu,'(a20)') "  conv_thr =  1.0d-8"
WRITE(ofu,'(a1)') "/"
!
!Write empty "ions" and "cell" sections
WRITE(ofu,*) ""
WRITE(ofu,'(a5)') "&IONS"
WRITE(ofu,'(a1)') "/"
WRITE(ofu,*) ""
WRITE(ofu,'(a5)') "&CELL"
WRITE(ofu,'(a1)') "/"
!
!Write mass of species
!Note: the user will have to append the pseudopotential file names
WRITE(ofu,*) ""
WRITE(ofu,'(a14)') "ATOMIC_SPECIES"
DO i=1,SIZE(aentries,1)
  CALL ATOMSPECIES(aentries(i,1),species)
  CALL ATOMMASS(species,smass)
  WRITE(ofu,'(a2,2X,f9.3,2X,a)') species, smass, TRIM(ppfiles(i))
ENDDO
!
IF( .NOT.cubic ) THEN
  !Write cell parameters
  WRITE(ofu,*) ""
  WRITE(ofu,'(a24)') "CELL_PARAMETERS angstrom"
  WRITE(ofu,201) H(1,1), H(1,2), H(1,3)
  WRITE(ofu,201) H(2,1), H(2,2), H(2,3)
  WRITE(ofu,201) H(3,1), H(3,2), H(3,3)
  201 FORMAT(3(f16.8,2X))
ENDIF
!
!Write atom coordinates
WRITE(ofu,*) ""
msg = "ATOMIC_POSITIONS"
IF( isreduced ) THEN
  msg = TRIM(ADJUSTL(msg))//" crystal"
ELSE
  msg = TRIM(ADJUSTL(msg))//" angstrom"
ENDIF
WRITE(ofu,'(a)') TRIM(msg)
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
  WRITE(ofu,'(a)') TRIM(msg)
ENDDO
210 FORMAT(a2,2X,3(f16.8,2X))
!
!Write k points section
WRITE(ofu,*) ""
WRITE(ofu,'(a18)') "K_POINTS automatic"
WRITE(ofu,'(a)') "2 2 2  0 0 0"
!
!Write forces on atoms
IF( fx>0 .AND. fy>0 .AND. fz>0 ) THEN
  WRITE(ofu,*) ""
  msg = "ATOMIC_FORCES"
  WRITE(ofu,'(a)') TRIM(msg)
  DO i=1,SIZE(P,1)
    CALL ATOMSPECIES(P(i,4),species)
    WRITE(msg,210) species, AUX(i,fx), AUX(i,fy), AUX(i,fz)
    WRITE(ofu,'(a)') TRIM(msg)
  ENDDO
ENDIF
GOTO 500
!
!
!
500 CONTINUE
IF(ofu.NE.6) THEN
  CLOSE(ofu)
ENDIF
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
