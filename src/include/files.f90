MODULE files
!
!**********************************************************************************
!*  FILES                                                                         *
!**********************************************************************************
!* This module contains subroutines performing special                            *
!* transformations on files for ATOMSK.                                           *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 15 Oct. 2014                                     *
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
!* CHECHFILE           checks if a file exists and behaves accordingly            *
!* SET_OUTPUT          sets one or all output to TRUE or FALSE                    *
!* FIND_INOUT          among 2files, determine which will be input or output      *
!* NAME_OUTFILE        names an output file based on the name of an input file    *
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
IF(fstatus=='writ') THEN   !if we want to write the file
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
! SET_OUTPUT
! This subroutine adds or removes output format(s)
! from the array "outfileformats".
!********************************************************
SUBROUTINE SET_OUTPUT(outfileformats,extension,l_value)
!
IMPLICIT NONE
CHARACTER(LEN=5):: extension !file extension to add to array
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats, tempout
LOGICAL:: l_value  !TRUE if extension must be added,
                   !FALSE if it must be removed
INTEGER:: i,j
INTEGER:: sizeini
!
IF(ALLOCATED(tempout)) DEALLOCATE(tempout)
sizeini = 0
IF(ALLOCATED(outfileformats)) THEN
  sizeini = SIZE(outfileformats)
ENDIF
!
IF(l_value) THEN
  !the file format must be added to the array
  IF(ALLOCATED(outfileformats)) THEN
    ALLOCATE(tempout(sizeini+1))
    DO i=1,SIZE(outfileformats)
      tempout(i) = outfileformats(i)
    ENDDO
    tempout(SIZE(outfileformats)+1) = extension
    DEALLOCATE(outfileformats)
    ALLOCATE(outfileformats(SIZE(tempout)))
    outfileformats = tempout
    DEALLOCATE(tempout)
  ELSE
    ALLOCATE(outfileformats(1))
    outfileformats(1) = extension
  ENDIF
  !
ELSE
  !one or all file format(s) must be removed from the array
  IF(extension=='all  ' ) THEN
    IF(ALLOCATED(outfileformats)) DEALLOCATE(outfileformats)
    RETURN
  ELSE
    IF(ALLOCATED(outfileformats) .AND. SIZE(outfileformats)>0) THEN
      !First, check if the file format actually exists in array
      j=0
      DO i=1,SIZE(outfileformats)
        IF( outfileformats(i)==extension ) j=1
      ENDDO
      !
      IF(j>0) THEN
        IF(sizeini==1) THEN
          !the extension exists in the array, and is the only entry: wipe out array
          DEALLOCATE(outfileformats)
          RETURN
        ELSE
          !the extension exists in the array => remove it
          ALLOCATE(tempout(SIZE(outfileformats)-1))
          j=0
          DO i=1,SIZE(outfileformats)
            IF( outfileformats(i).NE.extension ) THEN
              j=j+1
              tempout(j) = outfileformats(i)
            ENDIF
          ENDDO
          DEALLOCATE(outfileformats)
          ALLOCATE(outfileformats(SIZE(tempout)))
          outfileformats = tempout
          DEALLOCATE(tempout)
        ENDIF
        !
      ELSE
        !this file extension was not in the array => quit
        RETURN
      ENDIF
      !
    ELSE
      !array "outfileformats" is unallocated or has zero size, nothing to remove
      IF(ALLOCATED(outfileformats)) DEALLOCATE(outfileformats)
      RETURN
    ENDIF
  ENDIF
  !
ENDIF
!
!
END SUBROUTINE SET_OUTPUT
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
CHARACTER(LEN=128):: msg
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
!********************************************************
! NAME_OUTFILE
! This subroutine finds a name for an output file
! based on the name of an input file and the extension
! for the output file.
! For instance if an input file 'test' is to be converted
! to .def, the output file will be 'test.def'.
! But if the name of the input file is 'test.abc',
! the name of the output file must be 'test.def',
! and not 'test.abc.def'.
!********************************************************
SUBROUTINE NAME_OUTFILE(inputfile,outputfile,outfileformat)
!
IMPLICIT NONE
CHARACTER(LEN=5),INTENT(IN):: outfileformat
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=*),INTENT(OUT):: outputfile
CHARACTER(LEN=4096):: test, test2
INTEGER:: strlength, strlength2
!
test = TRIM(ADJUSTL(inputfile))
!Find where the extension of inputfile starts
!(we assume here that it has to be after the last path separator, if any)
strlength = SCAN(test,'.',BACK=.TRUE.)
strlength2 = SCAN(test,pathsep,BACK=.TRUE.)
IF(strlength2>=LEN_TRIM(test)) strlength2=0
IF(strlength>strlength2) test=test(1:strlength-1)
!
!suffix must not have a dot
test2 = TRIM(ADJUSTL(outfileformat))
IF(test2(1:1)==".") test2 = test2(2:)
!
!Set the name of outputfile
outputfile = TRIM(ADJUSTL(test))//'.'//TRIM(ADJUSTL(test2))
!
END SUBROUTINE NAME_OUTFILE
!
!
!
END MODULE files
