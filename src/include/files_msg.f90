MODULE files_msg
!
!**********************************************************************************
!*  FILES_MSG                                                                     *
!**********************************************************************************
!* This module contains subroutines checking the status of files                  *
!* and interacting with user.                                                     *
!**********************************************************************************
!* (C) Aug. 2024 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 20 Aug. 2024                                     *
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
!* List of subroutines in this module:                                            *
!* CHECKFILE           checks if a file exists, or prompt user what to do         *
!* FIND_INOUT          among 2files, determine or ask which is input or output    *
!**********************************************************************************
!
!
USE comv
USE constants
USE messages
!
IMPLICIT NONE
!
!
CONTAINS
!
!
!********************************************************
!  CHECKFILE
!  This subroutine checks if a file exists. Next
!  action depends on the file status:
!  - if the file has to be read but doesn't exist, then
!    prompt the user for another file to read;
!  - if the file has to be written but already exists,
!    then prompt the user whether to overwrite the file
!    or chose another name.
!********************************************************
!
SUBROUTINE CHECKFILE(filename,fstatus)
!
IMPLICIT NONE
CHARACTER(LEN=*):: filename
CHARACTER(LEN=1):: test
CHARACTER(LEN=*):: fstatus   !file status; 'read' or 'write'
CHARACTER(LEN=4096):: temp
LOGICAL:: fileexists
!
INQUIRE(FILE=filename,EXIST=fileexists)
IF(fstatus(1:1)=='w' .OR. fstatus(1:1)=='W') THEN   !if we want to write the file
  IF(.NOT.overw) THEN
    DO WHILE(fileexists)
      CALL ATOMSK_MSG(4,(/filename/),(/0.d0/))
      test = ""
      READ(*,*,END=100,ERR=100) test
      100 CONTINUE
      IF( LEN_TRIM(test)<=0 ) test=langyes
      IF(test==langyes) THEN
        CALL ATOMSK_MSG(5,(/filename/),(/0.d0/))
        !Delete the file if it exists
        !INQUIRE(FILE=filename,OPENED=fileisopen)
        !IF(.NOT.fileisopen) OPEN(UNIT=35,FILE=filename)
        !CLOSE(UNIT=35,STATUS='DELETE')
        fileexists = .FALSE.
      ELSEIF(test==langBigYes) THEN
        CALL ATOMSK_MSG(6,(/''/),(/0.d0/))
        fileexists = .FALSE.
        overw = .TRUE.
      ELSE
        CALL ATOMSK_MSG(7,(/''/),(/0.d0/))
        READ(*,*) filename
        INQUIRE(FILE=filename,EXIST=fileexists)
      ENDIF
    ENDDO
  ENDIF
!
ELSE  !if we want to read the file
  DO WHILE(.NOT.fileexists)
    CALL ATOMSK_MSG(8,(/filename/),(/0.d0/))
    READ(*,*) temp
    IF(TRIM(ADJUSTL(temp))==system_ls) THEN
      CALL SYSTEM(system_ls)
      CALL ATOMSK_MSG(9,(/''/),(/0.d0/))
      READ(*,*) filename
    ELSE
      filename = TRIM(ADJUSTL(temp))
    ENDIF
    INQUIRE(FILE=filename,EXIST=fileexists)
  ENDDO
!
ENDIF
!
END SUBROUTINE CHECKFILE
!
!
!********************************************************
! FIND_INOUT
! This subroutine, provided two file names, determines
! which file(s) exist. The existing file will become
! the "inputfile", and the non-existing one will become
! the "outputfile". If both files exist (or do not exist)
! then the user is prompted which file to use as
! "inputfile".
!********************************************************
SUBROUTINE FIND_INOUT(inputfile,outputfile)
!
IMPLICIT NONE
CHARACTER(LEN=1):: answer
CHARACTER(LEN=4096):: test
CHARACTER(LEN=4096),INTENT(INOUT):: inputfile,outputfile
LOGICAL:: fileexists
!
test = ''
!
IF(inputfile.NE.'') THEN
    INQUIRE(FILE=inputfile,EXIST=fileexists)
    IF(fileexists) THEN
      IF(outputfile.NE.'') THEN
        INQUIRE(FILE=outputfile,EXIST=fileexists)
        IF(fileexists) THEN
          !If both inputfile and outputfile already exist, we're in trouble:
          !we have to ask the user what to do
          answer=""
          DO WHILE( answer.NE."1" .AND. answer.NE."2" )
            CALL ATOMSK_MSG(700,(/inputfile,outputfile/),(/0.d0/))
            READ(*,*) answer
          ENDDO
          IF(answer=="1") THEN
            inputfile = TRIM(inputfile)
            outputfile = TRIM(outputfile)
          ELSEIF(answer=="2") THEN
            test = inputfile
            inputfile = outputfile
            outputfile = test
          ENDIF
        ELSE
          !if outputfile doesn't exist then inputfile is the input file
          inputfile = TRIM(inputfile)
          outputfile = TRIM(outputfile)
        ENDIF
      ENDIF
    !
    ELSE !i.e. if inputfile doesn't exist
      IF(outputfile.NE.'') THEN
        INQUIRE(FILE=outputfile,EXIST=fileexists)
        IF(fileexists) THEN
          !If outputfile exists then it is the input file
          test = inputfile
          inputfile = outputfile
          outputfile = test
        ELSE
          !If none of the files specified by the user exist, we're in trouble:
          !we have to ask the user what to do
          CALL ATOMSK_MSG(701,(/inputfile,outputfile/),(/0.d0/))
          READ(*,*) inputfile
          outputfile = ''
        ENDIF
      ENDIF
    ENDIF
  ELSE !i.e. if file1=='' (which also implies that file2=='')
    CALL ATOMSK_MSG(702,(/inputfile,outputfile/),(/0.d0/))
    READ(*,*) inputfile
    outputfile = ''
  ENDIF
!
END SUBROUTINE FIND_INOUT
!
!
!
END MODULE files_msg
