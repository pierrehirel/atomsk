MODULE mode_interactive
!
!**********************************************************************************
!*  MODE_INTERACTIVE                                                              *
!**********************************************************************************
!* This module offers a command-line interpreter (CLI) where the user can type    *
!* commands that make calls to various routines of Atomsk.                        *
!**********************************************************************************
!* (C) May 2010 - Pierre Hirel                                                    *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 05 Nov. 2025                                     *
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
USE crystallography
USE elasticity
USE messages
USE files
USE random
USE subroutines
USE guess_form
USE read_cla
USE modes
!
CONTAINS
!
!
SUBROUTINE INTERACT()
!
!
IMPLICIT NONE
CHARACTER(LEN=1):: answer
CHARACTER(LEN=2):: species, species2 !atom species
CHARACTER(LEN=2),DIMENSION(20):: create_species !chemical species of atoms (mode create)
CHARACTER(LEN=5):: zone
CHARACTER(LEN=8):: date
CHARACTER(LEN=10):: time
CHARACTER(LEN=10):: create_struc  !lattice type (mode create)
CHARACTER(LEN=12):: mode     !mode in which the program runs
CHARACTER(LEN=16):: helpsection
CHARACTER(LEN=32),DIMENSION(3):: create_Miller  !Miller indices for mode "create"
CHARACTER(LEN=32),DIMENSION(9):: criteria  !stability criteria not satisfied by Cij tensor
CHARACTER(LEN=128):: cwd, msg
CHARACTER(LEN=128):: prompt
CHARACTER(LEN=128):: test, temp
CHARACTER(LEN=128):: question, solution
CHARACTER(LEN=128):: username
CHARACTER(LEN=4096):: inputfile, outputfile, prefix
CHARACTER(LEN=4096):: instruction, command   !instruction given by the user
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats
CHARACTER(LEN=16),DIMENSION(58),PARAMETER:: optnames = & !names of options
  & (/ "add-atoms       ", "addatoms        ", "add-shells      ", "addshells       ", &
  &    "alignx          ", "bind-shells     ", "bs              ", "box             ", &
  &    "cell            ", "change_cell     ", "change_box      ", "center          ", &
  &    "crack           ", "cut             ", "deform          ", "def             ", &
  &    "dislocation     ", "disloc          ", "disturb         ", "duplicate       ", &
  &    "dup             ", "fix             ", "freeze          ", "fractional      ", &
  &    "frac            ", "mirror          ", "orient          ", "orthogonal-cell ", &
  &    "orthocell       ", "properties      ", "prop            ", "rebox           ", &
  &    "reduce-cell     ", "remove-atom     ", "rmatom          ", "remove-doubles  ", &
  &    "rmd             ", "remove-property ", "rmprop          ", "remove-shells   ", &
  &    "roll            ", "rmshells        ", "rotate          ", "rot             ", &
  &    "roundoff        ", "round-off       ", "select          ", "shear           ", &
  &    "shift           ", "sort            ", "substitute      ", "sub             ", &
  &    "swap            ", "torsion         ", "unit            ", "unskew          ", &
  &    "velocity        ", "wrap            "  &
  &/)
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
CHARACTER(LEN=32),DIMENSION(:),ALLOCATABLE:: var_names   !names of user-defined variables
CHARACTER(LEN=32),DIMENSION(:),ALLOCATABLE:: array_names !names of user-defined arrays
INTEGER:: a1, a2, a3
INTEGER:: i, j, k, status
INTEGER:: try, maxtries
INTEGER,DIMENSION(2):: NT_mn
INTEGER,DIMENSION(8):: timeval !values for DATE_AND_TIME function
LOGICAL:: cubic !is the lattice cubic?
LOGICAL:: exists !does file or directory exist?
LOGICAL:: WrittenToFile  !was the system written to a file?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
REAL(dp):: smass, snumber  !atomic mass, atomic number
REAL(dp):: C11, C22, C33, C12, C13, C23, C44, C55, C66  !elastic constants
REAL(dp):: x, y, z  !coordinates of an atom
REAL(dp),DIMENSION(3):: create_a0    !the lattice constants (mode create)
REAL(dp),DIMENSION(3,3):: Huc !Base vectors of the unit cell
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystal orientation
REAL(dp),DIMENSION(3,3):: rot_matrix  !rotation matrix
REAL(dp),DIMENSION(9,9):: C_tensor, S_tensor   !elastic and compliance tensor
REAL(dp),DIMENSION(:),ALLOCATABLE:: randarray  !random numbers
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries !array containing atomic number, N atoms
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, Ptemp !atomic positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S        !shell positions (is any)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX      !auxiliary properties of atoms/shells
REAL(dp),DIMENSION(:),ALLOCATABLE:: variables  !variables that user can define and work with
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE:: arrays !arrays that user can define and work with
!
!
!Initialize global variables
!In interactive mode verbosity cannot be 0 or 2
!If it is the case, reset it to 1
IF(verbosity==0 .OR. verbosity==2) verbosity=1
!
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
IF(ALLOCATED(variables)) DEALLOCATE(variables)
IF(ALLOCATED(var_names)) DEALLOCATE(var_names)
IF(ALLOCATED(arrays)) DEALLOCATE(arrays)
IF(ALLOCATED(array_names)) DEALLOCATE(array_names)
maxtries=5
C11 = 0.d0
C22 = 0.d0
C33 = 0.d0
C12 = 0.d0
C13 = 0.d0
C23 = 0.d0
C44 = 0.d0
C55 = 0.d0
C66 = 0.d0
H(:,:) = 0.d0
Huc(:,:) = 0.d0
 C_tensor(:,:) = 0.d0
 S_tensor(:,:) = 0.d0
ORIENT(:,:) = 0.d0
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
!Information message
CALL ATOMSK_MSG(4021,(/''/),(/0.d0/))
!
!
!
100 CONTINUE
CALL ATOMSK_MSG(4022,(/''/),(/0.d0/))
!
DO
  instruction = ""
  command = ""
  nerr=0  !do not exit program upon errors
  i = 0
  j = 0
  x=0.d0
  IF(ALLOCATED(randarray)) DEALLOCATE(randarray)
  !
  !Get current working directory
  CALL getcwd(cwd)
  IF( cwd=="/home/"//username ) THEN
    cwd="~"
  ELSE
    !Keep only what is after the last separator
    k = INDEX(cwd,pathsep,BACK=.TRUE.)
    cwd = TRIM(ADJUSTL(cwd(k+1:)))
  ENDIF
  !
  !Define the Atomsk prompt (contains the path to current directory)
  WRITE(prompt,*) TRIM(ADJUSTL(username))//"@atomsk:"//TRIM(ADJUSTL(cwd))//"> "
  !Set colour for the prompt
  prompt = COLOUR_MSG(prompt,colourprompt)
  prompt = TRIM(prompt)//ACHAR(27)//"[0m"
  !Display the Atomsk prompt!
  !k = LEN_TRIM(prompt)
  WRITE(*,'(a)',ADVANCE='NO') TRIM(prompt)//" "
  !
  !Read the instruction given by the user
  READ(*,'(a4096)',ERR=500,END=500) instruction
  instruction = ADJUSTL(instruction)
  !
  !Try to interpret the instruction and act accordingly
  IF( LEN_TRIM(instruction) > 0 ) THEN
    !
    !Remove equal sign if any
    i = SCAN(instruction,"=")
    IF(i>0) THEN
      instruction(i:i) = " "
    ENDIF
    j = SCAN(instruction," ")
    command = TRIM(ADJUSTL(instruction(1:j)))
    !
    IF( command(1:2)=="#!" ) THEN
      IF( SCAN(command,"atomsk")>0 ) THEN
        !We are reading a script file
      ENDIF
    ENDIF
    !
    !Ignore lines starting with "#" and empty lines
    IF( command(1:1).NE."#" ) THEN
      !Interpret the command
      SELECT CASE(command)
      !
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!         COMMANDS FOR INFO
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CASE("license","licence")
        CALL DISPLAY_LICENSE()
        !
      CASE("version","ver")
        CALL DISPLAY_COPYRIGHT()
      !
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!         COMMANDS FOR PROGRAM BEHAVIOR
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CASE("ignore","ig")
        IF( ignore ) THEN
          WRITE(*,*) " ignore = FALSE"
          ignore = .FALSE.
        ELSE
          WRITE(*,*) " ignore = TRUE"
          ignore = .TRUE.
        ENDIF
      !
      CASE("overwrite","ow")
        IF( overw ) THEN
          WRITE(*,*) " overw = FALSE"
          overw = .FALSE.
        ELSE
          WRITE(*,*) " overw = TRUE"
          overw = .TRUE.
        ENDIF
        !
      CASE("language","lang")
        j = SCAN(instruction," ")
        temp = TRIM(instruction(j:))
        IF( LEN_TRIM(temp)>0 ) THEN
          IF(temp=="fr") THEN
            lang = "fr"
          ELSEIF(temp=="de") THEN
            lang = "de"
          ELSE
            lang = "en"
          ENDIF
        ENDIF
        WRITE(*,*) " language = "//TRIM(lang)
      !
      CASE("verbosity","v")
        j = SCAN(instruction," ")
        temp = TRIM(instruction(j:))
        IF( LEN_TRIM(temp)>0 ) THEN
          READ(temp,*,ERR=400,END=400) i
          verbosity = MIN(4,ABS(i))
        ENDIF
        WRITE(*,'(a13,i1)') " verbosity = ", verbosity
      !
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!         SPECIAL COMMANDS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CASE("cd","CD")
        temp = TRIM(ADJUSTL( instruction(3:) ))
        IF( LEN_TRIM(temp) <= 0 .OR. temp=="~" ) THEN
          !Go to root directory
#if defined(WINDOWS)
          CALL CHDIR("C:")
#else
          CALL CHDIR("/home/"//username)
#endif
        ELSEIF( INDEX(temp,'!§;,?:#&{}<>=')==0 ) THEN
          !Test if directory exists
          INQUIRE( FILE="."//pathsep//TRIM(temp)//pathsep//"." , EXIST=exists )
          IF( exists ) THEN
            CALL CHDIR(temp)
          ELSE
            !Directory does not exist => display error message
            CALL ATOMSK_MSG(813,(/''/),(/0.d0/))
          ENDIF
        ELSE
          !Strange characters in dir name => display error message
          CALL ATOMSK_MSG(813,(/''/),(/0.d0/))
        ENDIF
        !
      CASE("ls","dir")
        CALL SYSTEM(system_ls)
        !
      CASE("print")
        !Print information about the system on the screen
        temp = TRIM(ADJUSTL( instruction(6:) ))
        IF( temp=="atoms" ) THEN
          IF( ALLOCATED(P) .AND. SIZE(P,1)>0 ) THEN
            WRITE(temp,*) SIZE(P,1)
            WRITE(*,*) "  "//TRIM(ADJUSTL(temp))//" atoms"
            DO i=1,MIN(20,SIZE(P,1))
              CALL ATOMSPECIES(P(i,4),species)
              WRITE(*,'(2X,a2,2X,3(f12.6,2X))') species, P(i,1), P(i,2), P(i,3)
            ENDDO
            IF( SIZE(P,1)>20 ) THEN
              WRITE(*,*) "  (...discontinued...)"
            ENDIF
          ENDIF
          !
        ELSEIF( temp=="box" ) THEN
          DO i=1,3
            WRITE(*,'(a4,i1,a3,3f12.6)') "   H", i, " = ", H(i,1), H(i,2), H(i,3)
          ENDDO
          !
        ELSEIF( temp=="Cij" ) THEN
          IF( ANY(C_tensor(:,:).NE.0.d0) ) THEN
            !Print elastic tensor
            WRITE(*,*) "  Current elastic tensor Cij (GPa):"
            DO i=1,9
              WRITE(*,'(2X,9(f12.6,2X))') (C_tensor(i,j) , j=1,9)
            ENDDO
          ELSE
            WRITE(*,*) "Cij not set"
          ENDIF
          !
        ELSEIF( temp=="Sij" ) THEN
          IF( ANY(S_tensor(:,:).NE.0.d0) ) THEN
            !Print elastic tensor
            WRITE(*,*) "  Current compliance tensor Sij (GPa):"
            DO i=1,9
              WRITE(*,'(2X,9(f12.6,2X))') (S_tensor(i,j) , j=1,9)
            ENDDO
          ELSE
            WRITE(*,*) "Sij not set"
          ENDIF
          !
        ELSE
          !Default: print box vectors and atom positions
          IF( ALLOCATED(P) .AND. SIZE(P,1)>0 ) THEN
            DO i=1,3
              WRITE(*,'(a4,i1,a3,3f12.6)') "   H", i, " = ", H(i,1), H(i,2), H(i,3)
            ENDDO
            WRITE(temp,*) SIZE(P,1)
            WRITE(*,*) "  "//TRIM(ADJUSTL(temp))//" atoms"
            DO i=1,MIN(20,SIZE(P,1))
              CALL ATOMSPECIES(P(i,4),species)
              WRITE(*,'(2X,a2,2X,3(f12.6,2X))') species, P(i,1), P(i,2), P(i,3)
            ENDDO
            IF( SIZE(P,1)>20 ) THEN
              WRITE(*,*) "  (...discontinued...)"
            ENDIF
          ENDIF
        ENDIF
        !
      CASE("pwd","PWD")
        CALL SYSTEM("pwd")
        !
      CASE("time","TIME","date","DATE","day","DAY","month","MONTH","year","YEAR")
        CALL DATE_AND_TIME(DATE, TIME, ZONE, TIMEVAL)
        WRITE(*,*) DATE(1:4)//"-"//DATE(5:6)//"-"//DATE(7:8)//" "//TIME(1:2)//":"//TIME(3:4)
        !
      CASE("whoami")
        !Print user name
        WRITE(*,*) "  "//TRIM(ADJUSTL(username))
        !
      CASE("clear")
        IF( ALLOCATED(P) .AND. .NOT.WrittenToFile ) THEN
          !User may have forgotten to write file => Display a warning
          CALL ATOMSK_MSG(4713,(/'erase it'/),(/0.d0/))
          READ(*,*) answer
          IF( .NOT. (answer==langyes .OR. answer==langBigYes) ) THEN
            GOTO 400
          ENDIF
        ENDIF
        !Wipe out everything from memory
        nerr=0
        nwarn=0
        WrittenToFile = .FALSE.
        C11 = 0.d0
        C22 = 0.d0
        C33 = 0.d0
        C12 = 0.d0
        C13 = 0.d0
        C23 = 0.d0
        C44 = 0.d0
        C55 = 0.d0
        C66 = 0.d0
        H(:,:) = 0.d0
        ORIENT(:,:) = 0.d0
        C_tensor(:,:) = 0.d0
        S_tensor(:,:) = 0.d0
        IF(ALLOCATED(comment)) DEALLOCATE(comment)
        IF(ALLOCATED(P)) DEALLOCATE(P)
        IF(ALLOCATED(S)) DEALLOCATE(S)
        IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
        IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
        IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
        IF(ALLOCATED(options_array)) DEALLOCATE(options_array) !no option in this mode
        IF(ALLOCATED(outfileformats)) DEALLOCATE(outfileformats)
        !
      CASE("exit","quit","bye")
        IF( ALLOCATED(P) .AND. .NOT.WrittenToFile ) THEN
          !User may have forgotten to write file => Display a warning
          CALL ATOMSK_MSG(4713,(/'quit'/),(/0.d0/))
          READ(*,*) answer
          IF( answer==langyes .OR. answer==langBigYes ) THEN
            GOTO 500
          ENDIF
        ELSE
          GOTO 500
        ENDIF
        !
      CASE("help","aide","hilfe","wtf","WTF")
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
        IF( ANY(H(:,:)>1.d-12) ) THEN
          i=i+1
          WRITE(*,'(a14,3f9.3,a2)') "            | ", H(1,1), H(1,2), H(1,3), " |"
          WRITE(*,'(a14,3f9.3,a2)') "     cell = | ", H(2,1), H(2,2), H(2,3), " |"
          WRITE(*,'(a14,3f9.3,a2)') "            | ", H(3,1), H(3,2), H(3,3), " |"
        ENDIF
        IF( ALLOCATED(P) ) THEN
          i=i+1
          IF( ALLOCATED(S) ) THEN
            WRITE(*,'(a18, i9)') "  N ionic cores:  ", SIZE(P,1)
            WRITE(*,'(a18, i9)') "  N ionic shells: ", SIZE(S,1)
            j = SIZE(P,1)+SIZE(S,1)
          ELSE
            j = SIZE(P,1)
          ENDIF
          WRITE(*,'(a18, i9)') "  N particles:    ", j
          !Write atoms species and their number
          CALL FIND_NSP(P(:,4),aentries)
          IF( SIZE(aentries,1)>0 ) THEN
            msg = "  Species:"
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
          WRITE(*,*) " Auxiliary properties: "//TRIM(ADJUSTL(msg))
        ENDIF
        IF( ALLOCATED(SELECT) ) THEN
          i=i+1
          j=0
          DO k=1,SIZE(SELECT)
            IF( SELECT(k) ) j=j+1
          ENDDO
          WRITE(*,*) " Atoms selected: ", j
        ENDIF
        IF(i==0) THEN
          WRITE(*,*) " Nothing in memory"
        ENDIF
        !
      !Commands to read and write files
      CASE("read")
        !Verify that there is no previous system in memory
        answer = langyes
        IF( ALLOCATED(P) .AND. .NOT.WrittenToFile ) THEN
          !User may have forgotten to write file => Display a warning
          CALL ATOMSK_MSG(4713,(/'erase it'/),(/0.d0/))
          READ(*,*) answer
        ENDIF
        IF( answer==langyes .OR. answer==langBigYes ) THEN
          !Read the file
          inputfile = TRIM(ADJUSTL(instruction(6:)))
          CALL READ_AFF(inputfile,H,P,S,comment,AUXNAMES,AUX)
        ENDIF
        !
      CASE("write")
        !Write system to a file
        prefix = TRIM(ADJUSTL(instruction(7:)))
        CALL WRITE_AFF(prefix,outfileformats,H,P,S,comment,AUXNAMES,AUX)
        WrittenToFile = .TRUE.
        !
      CASE("box","H","cell")
        !User wants to define an orthogonal box
        i = SCAN(instruction," ")
        IF( LEN_TRIM(instruction(i+1:))>0 ) THEN
          READ(instruction(i+1:),*,END=400,ERR=400) x, y, z
        ELSE
          WRITE(*,'(a14)',ADVANCE="NO") "     H(1,1) = "
          READ(*,*,END=400,ERR=400) x
          WRITE(*,'(a14)',ADVANCE="NO") "     H(2,2) = "
          READ(*,*,END=400,ERR=400) y
          WRITE(*,'(a14)',ADVANCE="NO") "     H(3,3) = "
          READ(*,*,END=400,ERR=400) z
        ENDIF
        H(:,:) = 0.d0
        H(1,1) = x
        H(2,2) = y
        H(3,3) = z
        WRITE(*,'(a14,3f9.3,a2)') "            | ", H(1,1), H(1,2), H(1,3), " |"
        WRITE(*,'(a14,3f9.3,a2)') "     cell = | ", H(2,1), H(2,2), H(2,3), " |"
        WRITE(*,'(a14,3f9.3,a2)') "            | ", H(3,1), H(3,2), H(3,3), " |"
        !
      CASE("atom")
        !User adds an atom
        READ(instruction(5:),*,END=400,ERR=400) species, x, y, z
        IF( .NOT. ALLOCATED(P) ) THEN
          ALLOCATE(P(1,4))
          k=1
        ELSE
          IF(ALLOCATED(Ptemp)) DEALLOCATE(Ptemp)
          k=SIZE(P,1)
          ALLOCATE(Ptemp(k,4))
          Ptemp(:,:) = P(:,:)
          DEALLOCATE(P)
          ALLOCATE(P(k+1,4))
          DO i=1,k
            P(i,:) = Ptemp(i,:)
          ENDDO
          DEALLOCATE(Ptemp)
          k=k+1
        ENDIF
        P(k,1) = x
        P(k,2) = y
        P(k,3) = z
        CALL ATOMNUMBER(species,P(k,4))
        WRITE(*,'(2X,i3,2X,3(f12.6,2X))') NINT(P(k,4)), P(k,1), P(k,2), P(k,3)
        !
      CASE("C11")
        command = TRIM(ADJUSTL(instruction(4:)))
        IF( LEN_TRIM(command)>0 ) THEN
          READ(command,*,END=400,ERR=400) C11
        ENDIF
        WRITE(temp,'(f12.6)') C11
        WRITE(*,*) "  C11 = "//TRIM(ADJUSTL(temp))//" GPa"
        IF( DABS(C22)<1.d-12 ) THEN
          C22 = C11
          WRITE(temp,'(f12.6)') C22
          WRITE(*,*) "  C22 = "//TRIM(ADJUSTL(temp))//" GPa"
        ENDIF
        IF( DABS(C33)<1.d-12 ) THEN
          C33 = C22
          WRITE(temp,'(f12.6)') C33
          WRITE(*,*) "  C33 = "//TRIM(ADJUSTL(temp))//" GPa"
        ENDIF
        !
      CASE("C22")
        command = TRIM(ADJUSTL(instruction(4:)))
        IF( LEN_TRIM(command)>0 ) THEN
          IF( command=="C11" ) THEN
            C22 = C11
          ELSE
            READ(command,*,END=400,ERR=400) C22
          ENDIF
        ENDIF
        WRITE(temp,'(f12.6)') C22
        WRITE(*,*) "  C22 = "//TRIM(ADJUSTL(temp))//" GPa"
        IF( DABS(C33)<1.d-12 ) THEN
          C33 = C22
          WRITE(temp,'(f12.6)') C33
          WRITE(*,*) "  C33 = "//TRIM(ADJUSTL(temp))//" GPa"
        ENDIF
        !
      CASE("C33")
        command = TRIM(ADJUSTL(instruction(4:)))
        IF( LEN_TRIM(command)>0 ) THEN
          IF( command=="C11" ) THEN
            C33 = C11
          ELSEIF( command=="C22" ) THEN
            C33 = C22
          ELSE
            READ(command,*,END=400,ERR=400) C33
          ENDIF
        ENDIF
        WRITE(temp,'(f12.6)') C33
        WRITE(*,*) "  C33 = "//TRIM(ADJUSTL(temp))//" GPa"
        !
      CASE("C12","C21")
        command = TRIM(ADJUSTL(instruction(4:)))
        IF( LEN_TRIM(command)>0 ) THEN
          READ(command,*,END=400,ERR=400) C12
        ENDIF
        WRITE(temp,'(f12.6)') C12
        WRITE(*,*) "  C12 = "//TRIM(ADJUSTL(temp))//" GPa"
        IF( DABS(C13)<1.d-12 ) THEN
          C13 = C12
          WRITE(temp,'(f12.6)') C13
          WRITE(*,*) "  C13 = "//TRIM(ADJUSTL(temp))//" GPa"
        ENDIF
        IF( DABS(C23)<1.d-12 ) THEN
          C23 = C12
          WRITE(temp,'(f12.6)') C23
          WRITE(*,*) "  C23 = "//TRIM(ADJUSTL(temp))//" GPa"
        ENDIF
        !
      CASE("C13","C31")
        command = TRIM(ADJUSTL(instruction(4:)))
        IF( LEN_TRIM(command)>0 ) THEN
          READ(command,*,END=400,ERR=400) C13
        ENDIF
        WRITE(temp,'(f12.6)') C13
        WRITE(*,*) "  C13 = "//TRIM(ADJUSTL(temp))//" GPa"
        !
      CASE("C23","C32")
        command = TRIM(ADJUSTL(instruction(4:)))
        IF( LEN_TRIM(command)>0 ) THEN
          READ(command,*,END=400,ERR=400) C23
        ENDIF
        WRITE(temp,'(f12.6)') C23
        WRITE(*,*) "  C23 = "//TRIM(ADJUSTL(temp))//" GPa"
        !
      CASE("C44")
        command = TRIM(ADJUSTL(instruction(4:)))
        IF( LEN_TRIM(command)>0 ) THEN
          READ(command,*,END=400,ERR=400) C44
        ENDIF
        WRITE(temp,'(f12.6)') C44
        WRITE(*,*) "  C44 = "//TRIM(ADJUSTL(temp))//" GPa"
        IF( DABS(C55)<1.d-12 ) THEN
          C55 = C44
          WRITE(temp,'(f12.6)') C55
          WRITE(*,*) "  C55 = "//TRIM(ADJUSTL(temp))//" GPa"
        ENDIF
        IF( DABS(C66)<1.d-12 ) THEN
          C66 = C55
          WRITE(temp,'(f12.6)') C66
          WRITE(*,*) "  C66 = "//TRIM(ADJUSTL(temp))//" GPa"
        ENDIF
        !
      CASE("C55")
        command = TRIM(ADJUSTL(instruction(4:)))
        IF( LEN_TRIM(command)>0 ) THEN
          IF( command=="C44" ) THEN
            C55 = C44
          ELSE
            READ(command,*,END=400,ERR=400) C55
          ENDIF
        ENDIF
        WRITE(temp,'(f12.6)') C55
        WRITE(*,*) "  C55 = "//TRIM(ADJUSTL(temp))//" GPa"
        IF( DABS(C66)<1.d-12 ) THEN
          C66 = C55
          WRITE(temp,'(f12.6)') C66
          WRITE(*,*) "  C66 = "//TRIM(ADJUSTL(temp))//" GPa"
        ENDIF
        !
      CASE("C66")
        command = TRIM(ADJUSTL(instruction(4:)))
        IF( LEN_TRIM(command)>0 ) THEN
          IF( command=="C44" ) THEN
            C66 = C44
          ELSEIF( command=="C55" ) THEN
            C66 = C55
          ELSE
            READ(command,*,END=400,ERR=400) C66
          ENDIF
        ENDIF
        WRITE(temp,'(f12.6)') C66
        WRITE(*,*) "  C66 = "//TRIM(ADJUSTL(temp))//" GPa"
        !
      CASE("Cij","Ctensor")
        IF( .NOT. ANY( (/C11,C22,C33,C12,C13,C23,C44,C55,C66/)>0.d0 ) ) THEN
          WRITE(*,*) "  Must set C11, C22, C33, C12, C13, C23, C44, C55, C66, before constructing Cij tensor."
        ELSE
          !Convert elastic constants into a proper 9x9 elastic tensor
          CALL ELAST2TENSOR( (/C11,C22,C33,C23,C13,C12,C44,C55,C66/) , C_tensor)
          !Print elastic tensor
          WRITE(*,*) "  Current elastic tensor Cij (GPa):"
          DO i=1,9
            WRITE(*,'(2X,9(f12.6,2X))') (C_tensor(i,j) , j=1,9)
          ENDDO
          !Check tensor for stability criteria
          CALL CTENSOR_STABILITY(C_tensor,criteria)
          IF( ANY(LEN_TRIM(criteria)>0) ) THEN
            CALL ATOMSK_MSG(2762,criteria,(/0.d0/))
          ENDIF
          !Invert elastic tensor to get compliance tensor
          CALL INVMAT(C_tensor(1:6,1:6),S_tensor(1:6,1:6),i)
          !If i is non-zero then the inversion failed
          IF(i==0) THEN
            S_tensor(7:9,1:3) = S_tensor(4:6,1:3)
            S_tensor(7:9,4:6) = S_tensor(4:6,4:6)
            S_tensor(1:9,7:9) = S_tensor(1:9,4:6)
          ELSE
            CALL ATOMSK_MSG(2815,(/"C_tensor"/),(/0.d0/))
          ENDIF
        ENDIF
        !
      CASE("Sij","Stensor","compliance")
          IF( .NOT.ANY(S_tensor(:,:).NE.0.d0) ) THEN
            IF( .NOT. ANY( (/C11,C22,C33,C12,C13,C23,C44,C55,C66/)>0.d0 ) ) THEN
              WRITE(*,*) "  Must set C11, C22, C33, C12, C13, C23, C44, C55, C66, before constructing Cij tensor."
            ELSE
              !Convert elastic constants into a proper 9x9 elastic tensor
              CALL ELAST2TENSOR( (/C11,C22,C33,C23,C13,C12,C44,C55,C66/) , C_tensor)
              !Print elastic tensor
              WRITE(*,*) "  Current elastic tensor Cij (GPa):"
              DO i=1,9
                WRITE(*,'(2X,9(f12.6,2X))') (C_tensor(i,j) , j=1,9)
              ENDDO
              !Check tensor for stability criteria
              CALL CTENSOR_STABILITY(C_tensor,criteria)
              IF( ANY(LEN_TRIM(criteria)>0) ) THEN
                CALL ATOMSK_MSG(2762,criteria,(/0.d0/))
              ENDIF
              !Invert elastic tensor to get compliance tensor
              CALL INVMAT(C_tensor(1:6,1:6),S_tensor(1:6,1:6),i)
              !If i is non-zero then the inversion failed
              IF(i==0) THEN
                S_tensor(7:9,1:3) = S_tensor(4:6,1:3)
                S_tensor(7:9,4:6) = S_tensor(4:6,4:6)
                S_tensor(1:9,7:9) = S_tensor(1:9,4:6)
              ELSE
                CALL ATOMSK_MSG(2815,(/"C_tensor"/),(/0.d0/))
              ENDIF
            ENDIF
          ENDIF
          !
          IF( .NOT.ANY(S_tensor(:,:).NE.0.d0) ) THEN
            WRITE(*,*) "  S_tensor not set"
          ELSE
            !Print compliance tensor
            WRITE(*,*) "  Current compliance tensor Sij (GPa):"
            DO i=1,9
              WRITE(*,'(2X,9(f12.6,2X))') (S_tensor(i,j) , j=1,9)
            ENDDO
          ENDIF
        !
      CASE("variable","var")
        !User wants to define a new variable
        IF( .NOT.ALLOCATED(variables) ) THEN
          k=1
          IF(ALLOCATED(var_names)) DEALLOCATE(var_names)
          ALLOCATE(variables(1))
          variables(1) = 0.d0
          ALLOCATE(var_names(1))
          var_names(1) = " "
        ELSE
          k=SIZE(variables)+1
          CALL RESIZE_CHAR1(var_names,k)
          CALL RESIZE_DBLEARRAY1(variables,k)
        ENDIF
        i = SCAN(instruction," ")
        READ(instruction(i+1:),*,END=400,ERR=400) var_names(k), variables(k)

      !
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!         COMMANDS FOR MODES
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CASE("create")
        IF( ALLOCATED(P) .AND. .NOT.WrittenToFile ) THEN
          !User may have forgotten to write file => Display a warning
          CALL ATOMSK_MSG(4713,(/'continue'/),(/0.d0/))
          READ(*,*) answer
          IF( answer==langyes .OR. answer==langBigYes ) THEN
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
          ELSE
            GOTO 400
          ENDIF
        ENDIF
        create_struc = ""
        create_species(:) = ""
        create_a0(:) = 0.d0
        NT_mn(:) = 0
        IF( LEN_TRIM(instruction(7:)) > 0 ) THEN
          !Read parameters, store them into array cla(:)
          IF(ALLOCATED(cla)) DEALLOCATE(cla)
          ALLOCATE(cla(10))
          cla(:) = ""
          cla(1) = "--create"
          DO i=1,10
            READ(instruction(7:),*,ERR=205,END=205) (cla(j+1), j=1,i)
          ENDDO
          !
          205 CONTINUE
          !
          !Call command-line parameters interpreter
          CALL GET_CLA(cla,mode,options_array,outfileformats,pfiles,mode_param)
          !
          !Get parameters from mode_param
          READ(mode_param(1),*,END=400,ERR=400) create_struc
          !Get the lattice constant(s)
          READ(mode_param(2),*,END=400,ERR=400) create_a0(1)
          create_a0(2) = create_a0(1)
          create_a0(3) = create_a0(1)
          i=2
          IF(create_struc=='graphite' .OR. create_struc=='hcp') THEN
            i=i+1
            READ(mode_param(i),*,END=400,ERR=400) create_a0(3)
          ELSEIF(create_struc=='nanotube' .OR. create_struc=='NT' .OR. create_struc=='nt') THEN
            i=i+1
            READ(mode_param(i),*,END=400,ERR=400) NT_mn(1)
            i=i+1
            READ(mode_param(i),*,END=400,ERR=400) NT_mn(2)
          ENDIF
          !Get the atomic species for the structure
          !There must always be at least one species
          i=i+1
          READ(mode_param(i),*,END=400,ERR=400) create_species(1)
          !Check that it is a correct atom type
          CALL ATOMNUMBER(create_species(1),smass)
          IF(smass==0.d0) THEN
            CALL ATOMSK_MSG(801,(/TRIM(create_species(1))/),(/0.d0/))
            GOTO 400
          ENDIF
          !Check if other atoms are specified
          DO j=2,20
            i=i+1
            READ(mode_param(i),*,END=210,ERR=210) temp
            IF(LEN_TRIM(temp)<=2 .AND. temp(1:1).NE.'-' ) THEN
              READ(temp,*,END=210,ERR=210) species
              CALL ATOMNUMBER(species,smass)
              IF(smass==0.d0) THEN
                i=i-1
                GOTO 210
              ELSE
                create_species(j) = species
              ENDIF
            ELSE
              i=i-1
            ENDIF
          ENDDO
          210 CONTINUE
          !Check if a crystallographic orientation was given
          i=i+1
          READ(mode_param(i),*,END=220,ERR=220) temp
          IF(temp=="orient") THEN
            DO j=1,3
              i=i+1
              READ(mode_param(i),*,END=220,ERR=220) temp
              CALL INDEX_MILLER(temp,ORIENT(j,:),k)
              IF(k>0) GOTO 400
            ENDDO
          ENDIF
          220 CONTINUE
          IF(ALLOCATED(mode_param)) DEALLOCATE(mode_param)
          !
          !
        ELSE
          !user just typed "create" with no parameter => ask for the parameters
          !Ask for the lattice type
          CALL ATOMSK_MSG(4200,(/''/),(/0.d0/))
          READ(*,*,END=400,ERR=400) create_struc
          IF(create_struc=="q") GOTO 400
          !
          SELECT CASE(create_struc)
          CASE('sc','SC','fcc','FCC','bcc','BCC','diamond','dia','zincblende','zb','ZB','perovskite','per','rocksalt','rs','RS')
            cubic = .TRUE.
          CASE DEFAULT
            cubic = .FALSE.
          END SELECT
          !
          !Ask for lattice parameter(s)
          CALL ATOMSK_MSG(4201,(/'a0'/),(/0.d0/))
          READ(*,'(a128)',END=400,ERR=400) test
          IF(test=="q") GOTO 400
          READ(test,*,END=400,ERR=400) create_a0(1)
          create_a0(2) = create_a0(1)
          create_a0(3) = create_a0(1)
          k=3
          IF( create_struc=="nt" .OR. create_struc=="nanotube" ) THEN
            CALL ATOMSK_MSG(4203,(/'a0'/),(/0.d0/))
            READ(*,'(a128)',END=400,ERR=400) test
            IF(test=="q") GOTO 400
            READ(test,*,END=400,ERR=400) NT_mn(1), NT_mn(2)
          ELSE
            IF( create_struc=="hcp" .OR. create_struc=="graphite" ) THEN
              CALL ATOMSK_MSG(4201,(/'c0'/),(/0.d0/))
              READ(*,*,END=400,ERR=400) create_a0(3)
              k=5
            ENDIF
          ENDIF
          CALL ATOMSK_MSG(4202,(/''/),(/0.d0/))
          READ(*,'(a128)',END=230,ERR=230) test
          IF(test=="q") GOTO 400
          READ(test,*,END=230,ERR=230) create_species(1)
          test = ADJUSTL(test(3:))
          IF(LEN_TRIM(test)>0) THEN
            READ(test,*,END=230,ERR=230) create_species(2)
            test = ADJUSTL(test(3:))
            IF(LEN_TRIM(test)>0) THEN
              READ(test,*,END=230,ERR=230) create_species(3)
            ENDIF
          ENDIF
          230 CONTINUE
          !
          !If the lattice is cubic, ask for the orientation
          CALL ATOMSK_MSG(4204,(/''/),(/0.d0/))
          IF(ALLOCATED(mode_param)) DEALLOCATE(mode_param)
          ALLOCATE(mode_param(3))
          mode_param(:) = ""
          READ(*,*,END=400,ERR=400) mode_param(1), mode_param(2), mode_param(3)
          i=0
          DO j=1,3
            i=i+1
            READ(mode_param(i),*,END=400,ERR=400) create_Miller
          ENDDO
          !
        ENDIF
        !
        250 CONTINUE
        !Run the mode create (but don't write any output file)
        CALL CREATE_CELL(create_a0,create_struc,create_species,NT_mn,create_Miller, &
                        & options_array,outputfile,outfileformats,.FALSE.,H,P)
        !
        !
        !
        !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!         MISC.  COMMANDS
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CASE("Hello","hello","Hi","hi","bonjour","Bonjour")
        command(1:1) = StrUpCase(command(1:1))
        WRITE(*,*) TRIM(ADJUSTL(command))//", "//TRIM(ADJUSTL(username))//"."
        !
      CASE("atomsk","Atomsk","ATOMSK")
        !CALL DISPLAY_HEADER()
        CALL ATOMSK_MSG(4825,(/""/),(/0.d0/))
        !
      CASE("all-in-one","average","centrosymmetry","density","difference","electric-dipoles","edm", &
          &"electric-polarization", "interpolate","list","merge","nye","polycrystal","rdf","unwrap")
        CALL ATOMSK_MSG(4826,(/""/),(/0.d0/))
        !
      CASE("random","rand")
        !generate and display N random numbers between 0 and 1000
        j=SCAN(instruction," ")
        msg = TRIM(ADJUSTL(instruction(j+1:)))
        IF( LEN_TRIM(msg)>0 ) THEN
          READ(msg,*,ERR=400,END=400) a1
          IF( a1<=0 ) a1=1
        ELSE
          a1 = 1
        ENDIF
        CALL GEN_NRANDNUMBERS(a1,randarray)
        DO i=1,SIZE(randarray)
          WRITE(*,'(a2,f12.9)',ADVANCE="NO") "  ", randarray(i)
          IF(MOD(i,10)==0) WRITE(*,'(a1)',ADVANCE="YES") " "
        ENDDO
        WRITE(*,'(a1)',ADVANCE="YES") " "
        !
      CASE("colour","color")
        !change default colour for text messages in terminal
        j=SCAN(instruction," ")
        msg = TRIM(ADJUSTL(instruction(j+1:)))
        IF( LEN_TRIM(msg)>0 ) READ(msg,*,ERR=400,END=400) colourdef
        !
      CASE("1337")
        WRITE(*,*) "U R 1337 :-)"
        !
      CASE("motd")
        !Message of the day
        CALL DATE_MSG()
        !
      CASE("quizz","quiz","play","game")
        !The user wants to play: ask a question about a random atom
        CALL ATOMSK_MSG(4300,(/""/),(/DBLE(maxtries)/))
        try=0
        CALL GEN_NRANDNUMBERS(2,randarray)
        !randarray contains 2 random numbers between 0 and 1
        !First number decides the form of the question
        IF( randarray(1)<=0.5 ) THEN
          !Given an atom species, user must guess its atomic number
          !Use second random number to pick up an atom species
          snumber = DBLE( NINT( randarray(2)*100.d0 ) )
          CALL ATOMSPECIES(snumber,species)
          question = species
          WRITE(solution,*) NINT(snumber)
          solution = ADJUSTL(solution)
          i=0
          DO WHILE( i.NE.NINT(snumber) .AND. try<maxtries )
            !Ask for the atomic number of that atom
            CALL ATOMSK_MSG(4301,(/question/),(/1.d0,DBLE(try),DBLE(i)-snumber/))
            READ(*,*) temp
            !Try to read a number from that string
            READ(temp,*,ERR=300,END=300) i
            300 CONTINUE
            IF(i.NE.snumber) try=try+1
          ENDDO
        ELSEIF( randarray(1)<=0.75 ) THEN
          !Given an atomic number, user must guess its species
          !Use second random number to pick up an atom species
          snumber = DBLE( NINT( randarray(2)*100.d0 ) )
          CALL ATOMSPECIES(snumber,species)
          WRITE(question,*) NINT(snumber)
          question = ADJUSTL(question)
          solution = species
          temp=""
          DO WHILE( temp.NE.solution .AND. try<maxtries )
            !Ask for the species of that atom
            CALL ATOMSK_MSG(4301,(/question/),(/2.d0,DBLE(try)/))
            READ(*,*) temp
            IF(temp.NE.solution) try=try+1
          ENDDO
        ELSE
          !User must guess the next species in periodic table
          !Use second random number to pick up an atom species
          snumber = DBLE( NINT( randarray(2)*100.d0 ) )
          CALL ATOMSPECIES(snumber,species)
          CALL ATOMSPECIES(snumber+1.d0,species2)
          question = species
          solution = species2
          temp=""
          DO WHILE( temp.NE.solution .AND. try<maxtries )
            !Ask for the next atom in the periodic table
            CALL ATOMSK_MSG(4301,(/question/),(/3.d0,DBLE(try)/))
            READ(*,*) temp
            IF(temp.NE.solution) try=try+1
          ENDDO
        ENDIF
        !
        CALL ATOMSK_MSG(4302,(/solution/),(/DBLE(try)/DBLE(maxtries)/))
        !
      CASE("progress-bar")
        !This command is just to test the display of a progress bar
        j=SCAN(instruction," ")
        msg = TRIM(ADJUSTL(instruction(j+1:)))
        IF( LEN_TRIM(msg)>0 ) THEN
          READ(msg,*,ERR=400,END=400) progressbar
        ELSE
          DO i=1,100
            CALL DISPLAY_PROGBAR(DBLE(i),100.d0)
            CALL SYSTEM("sleep 0.1") !wait 0.1 second
          ENDDO
        ENDIF
        !
        !
      CASE DEFAULT
        IF( command=="rotate" .AND. .NOT.ALLOCATED(P) .AND. ANY(C_tensor(:,:).NE.0.d0) ) THEN
          !User asks for rotation, no atomic system is defined, elastic tensor is defined
          !Read axis of rotation and angle of rotation
          READ(instruction(7:),*,ERR=400,END=400) answer, smass
          IF(answer=='x' .OR. answer=='X') THEN
            a1 = 1
            a2 = 2
            a3 = 3
          ELSEIF(answer=='y' .OR. answer=='Y') THEN
            a1 = 2
            a2 = 3
            a3 = 1
          ELSEIF(answer=='z' .OR. answer=='Z') THEN
            a1 = 3
            a2 = 1
            a3 = 2
          ELSE
            !Directions will be simply a1=X, a2=Y, a3=Z
            a1 = 1
            a2 = 2
            a3 = 3
          ENDIF
          !convert the angle into radians
          smass = DEG2RAD(smass)
          !set the rotation matrix
          rot_matrix(:,:) = 0.d0
          rot_matrix(a1,a1) = 1.d0
          rot_matrix(a2,a2) = DCOS(smass)
          rot_matrix(a2,a3) = -DSIN(smass)
          rot_matrix(a3,a2) = DSIN(smass)
          rot_matrix(a3,a3) = DCOS(smass)
          !Rotate elastic tensor
          C_tensor = ROTELAST( C_tensor, rot_matrix )
          !Print elastic tensor
          WRITE(*,*) "Current elastic tensor Cij (GPa):"
          DO i=1,9
            WRITE(*,'(9(f12.6,2X))') (C_tensor(i,j) , j=1,9)
          ENDDO
          !Check tensor for stability criteria
          CALL CTENSOR_STABILITY(C_tensor,criteria)
          IF( ANY(LEN_TRIM(criteria)>0) ) THEN
            CALL ATOMSK_MSG(2762,criteria,(/0.d0/))
          ENDIF
          !
        ELSEIF( ANY( command == optnames(:) ) ) THEN
          !Apply the option
          IF(ALLOCATED(options_array)) DEALLOCATE(options_array)
          ALLOCATE(options_array(1))
          options_array(1) = "-"//TRIM(ADJUSTL(instruction))
          CALL OPTIONS_AFF(options_array,Huc,H,P,S,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
          DEALLOCATE(options_array)
          !System was changed: consider that modified version was not written to a file
          WrittenToFile = .FALSE.
          !
        ELSEIF( SCAN(instruction,"0123456789.+-*/^%")>0 .OR. INDEX(instruction,"pi")>0 &
              & .OR. INDEX(instruction,"rand")>0 ) THEN
          !Maybe it is a formula: try to interpret it
          i=1
          status=0
          command = instruction
          x = EXPREVAL(instruction,i,status)
          IF( status==0 ) THEN
            !No error: display result
            IF( IS_INTEGER(x,1.d-64) ) THEN
              !It is an integer, or close to it: display an integer
              WRITE(msg,*) NINT(x)
            ELSE
              !It is a real number
              IF( x>=1.d4 .OR. x<=1.d-4 ) THEN
                WRITE(msg,'(e21.12)') x
              ELSE
                WRITE(msg,'(f21.12)') x
              ENDIF
            ENDIF
            WRITE(*,*) "          "//TRIM(ADJUSTL(msg))
          ELSE
            IF(status==10) THEN
              !There was a division by zero
              CALL ATOMSK_MSG(816,(/""/),(/0.d0/))
            ELSE
              !Unable to convert that string into a number
              CALL ATOMSK_MSG(2813,(/instruction/),(/0.d0/))
            ENDIF
          ENDIF
          !
        ELSE
          !Unknown command
          CALL ATOMSK_MSG(4824,(/command/),(/0.d0/))
          !
        ENDIF
        !
      END SELECT
      !
    ENDIF
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
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
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
