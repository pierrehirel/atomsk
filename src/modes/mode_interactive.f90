MODULE mode_interactive
!
!**********************************************************************************
!*  MODE_INTERACTIVE                                                              *
!**********************************************************************************
!* This module offers a command-line interpreter (CLI) where the user can type    *
!* commands that make calls to various routines of atomsk.                        *
!**********************************************************************************
!* (C) May 2010 - Pierre Hirel                                                    *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 04 July 2014                                     *
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
USE guess_form
USE read_cla
USE mode_create
USE mode_normal
!
CONTAINS
!
!
SUBROUTINE INTERACT()
!
!
IMPLICIT NONE
CHARACTER(LEN=1):: answer
CHARACTER(LEN=2):: species !atom species
CHARACTER(LEN=2),DIMENSION(20):: create_species !chemical species of atoms (mode create)
CHARACTER(LEN=5):: outfileformat
CHARACTER(LEN=10):: create_struc  !lattice type (mode create)
CHARACTER(LEN=12):: mode     !mode in which the program runs
CHARACTER(LEN=16):: helpsection
CHARACTER(LEN=128):: msg, test
CHARACTER(LEN=128):: username
CHARACTER(LEN=4096):: inputfile, outputfile, prefix
CHARACTER(LEN=4096):: instruction, command   !instruction given by the user
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats
CHARACTER(LEN=16),DIMENSION(45):: optnames !names of options
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: mode_param  !parameters for some special modes
CHARACTER(LEN=4096),DIMENSION(:),ALLOCATABLE:: cla  !command-line arguments
CHARACTER(LEN=4096),DIMENSION(5):: pfiles !pfiles(1)=file1
                                          !pfiles(2)=file2
                                          !pfiles(3)=filefirst
                                          !pfiles(4)=filesecond
                                          !pfiles(5)=listfile
INTEGER:: i, j
INTEGER,DIMENSION(2):: NT_mn
LOGICAL:: WrittenToFile  !was the system written to a file?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
REAL(dp),DIMENSION(3):: create_a0    !the lattice constants (mode create)
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystal orientation
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries !array containing atomic number, N atoms
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P     !atomic positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S     !shell positions (is any)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX  !auxiliary properties of atoms/shells

!
!Initialize variables
WrittenToFile = .FALSE.
IF(ALLOCATED(comment)) DEALLOCATE(comment)
IF(ALLOCATED(P)) DEALLOCATE(P)
IF(ALLOCATED(S)) DEALLOCATE(S)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
IF(ALLOCATED(options_array)) DEALLOCATE(options_array) !no option in this mode
IF(ALLOCATED(outfileformats)) DEALLOCATE(outfileformats)
H(:,:) = 0.d0
ORIENT(:,:) = 0.d0
optnames(:) = (/ "add-atoms       ", "addatoms        ", "add-shells      ", "addshells       ", &
            &    "alignx          ", "bind-shells     ", "bs              ", "center          ", &
            &    "crack           ", "cut             ", "deform          ", "def             ", &
            &    "dislocation     ", "disloc          ", "disturb         ", "expand          ", &
            &    "e               ", "fix             ", "fractional      ", "frac            ", &
            &    "mirror          ", "orient          ", "properties      ", "prop            ", &
            &    "rebox           ", "remove-atom     ", "rmatom          ", "remove-doubles  ", &
            &    "rmd             ", "remove-property ", "rmprop          ", "remove-shells   ", &
            &    "rmshells        ", "rotate          ", "rot             ", "select          ", &
            &    "shear           ", "shift           ", "sort            ", "substitute      ", &
            &    "sub             ", "unit            ", "unskew          ", "velocity        ", &
            &    "wrap            "                                                              &
            &/)
!
!Get user name: this is environment-dependent
username=""
#if defined(WINDOWS)
CALL GET_ENVIRONMENT_VARIABLE('USERNAME',username)
#else
CALL GET_ENVIRONMENT_VARIABLE('USER',username)
#endif
!
msg = 'Entering INTERACT'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!Print a nice message
CALL DISPLAY_HEADER()
CALL ATOMSK_MSG(4021,(/''/),(/0.d0/))
CALL DATE_MSG()
!
!
!
100 CONTINUE
CALL ATOMSK_MSG(4022,(/''/),(/0.d0/))
!
DO
  instruction = ""
  nerr=0  !ignore errors
  !
  !The Atomsk prompt!
  WRITE(*,'(a8)',ADVANCE='NO') "atomsk> "
  !
  !Read the instruction given by the user
  READ(*,'(a4096)',ERR=500,END=500) instruction
  instruction = ADJUSTL(instruction)
  !
  !Try to interpret the instruction and act accordingly
  IF( LEN_TRIM(instruction) > 0 ) THEN
    !
    READ(instruction,*) command
    command = ADJUSTL(command)
    !
    SELECT CASE(command)
    !
    !Some speciel commands
    CASE("clear")
      !Wipe out everything from memory
      WrittenToFile = .FALSE.
      H(:,:) = 0.d0
      ORIENT(:,:) = 0.d0
      IF(ALLOCATED(comment)) DEALLOCATE(comment)
      IF(ALLOCATED(P)) DEALLOCATE(P)
      IF(ALLOCATED(S)) DEALLOCATE(S)
      IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
      IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
      IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
      IF(ALLOCATED(options_array)) DEALLOCATE(options_array) !no option in this mode
      IF(ALLOCATED(outfileformats)) DEALLOCATE(outfileformats)
      !
    CASE("ls","dir")
      CALL SYSTEM(system_ls)
      !
    CASE("exit","quit","bye")
      IF( ALLOCATED(P) .AND. .NOT.WrittenToFile ) THEN
        !User may have forgotten to write file => Display a warning
        CALL ATOMSK_MSG(4713,(/''/),(/0.d0/))
        READ(*,*) answer
        IF( answer==langyes .OR. answer==langBigYes ) THEN
          GOTO 500
        ENDIF
      ELSE
        GOTO 500
      ENDIF
      !
    CASE("help")
      helpsection = TRIM(ADJUSTL(instruction(6:)))
      IF( helpsection=="options" .OR. helpsection=="modes" ) THEN
        CALL DISPLAY_HELP(helpsection)
      ELSE
        CALL ATOMSK_MSG(4023,(/''/),(/0.d0/))
      ENDIF
      !
    CASE("memory","mem")
      !Display some info about array sizes
      i=0
      IF( ALLOCATED(P) ) THEN
        i=i+1
        IF( ALLOCATED(S) ) THEN
          WRITE(*,'(a17, i9)') " N ionic cores:  ", SIZE(P,1)
          WRITE(*,'(a17, i9)') " N ionic shells: ", SIZE(S,1)
          j = SIZE(P,1)+SIZE(S,1)
        ELSE
          j = SIZE(P,1)
        ENDIF
        WRITE(*,'(a17, i9)') " N particles:    ", j
        !Write atoms species and their number
        CALL FIND_NSP(P(:,4),aentries)
        IF( SIZE(aentries,1)>0 ) THEN
          msg = "Species:"
          DO j=1,SIZE(aentries,1)
            CALL ATOMSPECIES(aentries(j,1) , species)
            msg = TRIM(msg)//" "//species
            WRITE(test,*) NINT(aentries(j,2))
            msg = TRIM(msg)//"("//TRIM(ADJUSTL(test))//")"
          ENDDO
          WRITE(*,*) TRIM(msg)
        ENDIF
      ENDIF
      IF( ALLOCATED(AUX) ) THEN
        i=i+1
        msg = AUXNAMES(1)
        DO j=2,SIZE(AUXNAMES)
          msg = TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(AUXNAMES(j)))
        ENDDO
        WRITE(*,*) " Auxiliary properties: ", msg
      ENDIF
      IF(i==0) THEN
        WRITE(*,*) " Nothing in memory"
      ENDIF
      !
    CASE("version")
      CALL DISPLAY_COPYRIGHT()
      !
    !Commands to read and write files
    CASE("read")
      !Read a file
      inputfile = TRIM(ADJUSTL(instruction(6:)))
      CALL READ_AFF(inputfile,H,P,S,comment,AUXNAMES,AUX)
      !
    CASE("write")
      !Write system to a file
      prefix = TRIM(ADJUSTL(instruction(7:)))
      CALL WRITE_AFF(prefix,outfileformats,H,P,S,comment,AUXNAMES,AUX)
      WrittenToFile = .TRUE.
      !
    !Command for creating an atomic system
    CASE("create")
      200 CONTINUE
      create_struc = ""
      create_species(:) = ""
      create_a0(:) = 0.d0
      NT_mn(:) = 0
      IF( LEN_TRIM(instruction(7:)) > 0 ) THEN
        !Read parameters, store them into array cla(:)
        READ(instruction(7:),*,ERR=400,END=400)  
        !Call command-line parameters interpreter
        CALL GET_CLA(cla,mode,options_array,outfileformats,pfiles,mode_param)
      ELSE
        !user just typed "create" with no parameter => ask for the parameters
        CALL ATOMSK_MSG(4200,(/''/),(/0.d0/))
        READ(*,*,END=400,ERR=400) create_struc
        IF(create_struc=="q") GOTO 400
        CALL ATOMSK_MSG(4201,(/'a0'/),(/0.d0/))
        READ(*,'(a128)',END=400,ERR=400) test
        IF(test=="q") GOTO 400
        READ(test,*,END=400,ERR=400) create_a0(1)
        create_a0(2) = create_a0(1)
        create_a0(3) = create_a0(1)
        IF( create_struc=="nt" .OR. create_struc=="nanotube" ) THEN
          CALL ATOMSK_MSG(4203,(/'a0'/),(/0.d0/))
          READ(*,'(a128)',END=400,ERR=400) test
          IF(test=="q") GOTO 400
          READ(test,*,END=400,ERR=400) NT_mn(1), NT_mn(2)
        ELSEIF( create_struc=="hcp" .OR. create_struc=="graphite" ) THEN
          CALL ATOMSK_MSG(4201,(/'c0'/),(/0.d0/))
          READ(*,*,END=400,ERR=400) create_a0(3)
        ENDIF
        CALL ATOMSK_MSG(4202,(/''/),(/0.d0/))
        READ(*,'(a128)',END=400,ERR=400) test
        IF(test=="q") GOTO 400
        READ(test,*,END=400,ERR=400) create_species(1)
        test = ADJUSTL(test(3:))
        READ(test,*,END=250,ERR=250) create_species(2)
        test = ADJUSTL(test(3:))
        READ(test,*,END=250,ERR=250) create_species(3)
      ENDIF
      250 CONTINUE
      PRINT*, create_species(:)
      !Run the mode
      CALL CREATE_CELL(create_a0,create_struc,create_species,NT_mn,ORIENT,options_array,outputfile,outfileformats,.FALSE.,H,P)
      !
    !Misc. commands
    CASE("Hello","hello","Hi","hi","bonjour","Bonjour")
      command(1:1) = StrUpCase(command(1:1))
      WRITE(*,*) TRIM(ADJUSTL(command))//", "//TRIM(ADJUSTL(username))//"."
      !
    CASE DEFAULT
      IF( ANY( command == optnames(:) ) ) THEN
        !This is an option: check if an atomic system exists in memory
        IF( ALLOCATED(P) .AND. SIZE(P)>0 ) THEN
          !Apply the option
          IF(ALLOCATED(options_array)) DEALLOCATE(options_array)
          ALLOCATE(options_array(1))
          options_array(1) = "-"//TRIM(ADJUSTL(instruction))
          CALL OPTIONS_AFF(options_array,H,P,S,AUXNAMES,AUX,ORIENT,SELECT)
          DEALLOCATE(options_array)
        ELSE
          !The user wants to apply an option, but no system in memory
          CALL ATOMSK_MSG(2814,(/''/),(/0.d0/))
        ENDIF
        !
      ELSE
        WRITE(*,*) "Unknown command: "//TRIM(ADJUSTL(command))
      ENDIF
      !
    END SELECT
    !
  ENDIF
  !
  400 CONTINUE
  !
ENDDO
!
!
!
500 CONTINUE
IF(ALLOCATED(P)) DEALLOCATE(P)
IF(ALLOCATED(S)) DEALLOCATE(S)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
!
!
!
1000 CONTINUE
!
!
!
END SUBROUTINE INTERACT
!
!
END MODULE mode_interactive
