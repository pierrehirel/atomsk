MODULE mode_list
!
!**********************************************************************************
!*  MODE_LIST                                                                     *
!**********************************************************************************
!* This module converts a list of files to one or several other file(s).          *
!**********************************************************************************
!* (C) October 2013 - Pierre Hirel                                                *
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
!Load modules
USE comv
USE constants
USE functions
USE messages
USE files
USE subroutines
!SPECIAL CASE: this mode depends on "mode_normal.f90"
USE mode_normal
!
!
CONTAINS
!
SUBROUTINE LIST_XYZ(filelist,options_array,outputfile,outfileformats)
!
!Declare variables
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: filelist
CHARACTER(LEN=5):: curr_outfileformat !output file format for one file
CHARACTER(LEN=5):: list_outfileformat !output file format that applies to the whole list
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of output file formats
CHARACTER(LEN=128):: msg
CHARACTER(LEN=4096):: inputfile, outputfile, prefix
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
CHARACTER(LEN=4096),DIMENSION(:),ALLOCATABLE:: filearray !list of files to convert
LOGICAL:: fileexists  !does the file already exist?
LOGICAL:: ignorefile  !should this file be ignored?
INTEGER:: i, j
INTEGER:: Nfiles, Nignored !number of files converted, ignored
INTEGER:: strlength
!
!Initialize
Nfiles = 0
Nignored = 0
list_outfileformat = ""
!
!
!
msg = 'ENTERING CONVERT_AFF...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Initialize variables
!
CALL CHECKFILE(filelist,'read')
CALL ATOMSK_MSG(4001,(/''/),(/0.d0/))
!
!
!
100 CONTINUE
!First, count how many files exist in the list
OPEN(UNIT=50,FILE=filelist,STATUS='OLD',FORM='FORMATTED')
Nfiles=0
DO
  !Read the whole line
  READ(50,'(a4096)',ERR=120,END=120) inputfile
  IF(LEN_TRIM(inputfile).NE.0 .AND. inputfile(1:1).NE.'#') THEN
    IF( inputfile(1:3)=='all' ) THEN
      !'all <format>' means that all the following files have to be converted
      !into the <format> by default
      READ(inputfile(4:20),*) list_outfileformat
      CALL SET_OUTPUT(outfileformats,list_outfileformat,.TRUE.)
      CALL ATOMSK_MSG(4002,(/TRIM(list_outfileformat)/),(/0.d0/))
      msg = 'list_outfileformat: '//TRIM(list_outfileformat)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ELSE
      Nfiles = Nfiles+1
    ENDIF
  ENDIF
ENDDO
!
120 CONTINUE
IF(Nfiles<=0) THEN
  GOTO 335
ENDIF
!We know how many files have to be converted
!Save file names into array "filearray"
ALLOCATE( filearray(Nfiles) )
filearray(:) = ''
!
REWIND(50)
Nfiles=0
DO
  !Read the whole line
  READ(50,'(a4096)',ERR=140,END=140) inputfile
  IF(LEN_TRIM(inputfile).NE.0 .AND. inputfile(1:1).NE.'#') THEN
    Nfiles = Nfiles+1
    filearray(Nfiles) = ADJUSTL(inputfile)
  ENDIF
ENDDO
!
140 CONTINUE
CLOSE(50)
!
!
!
200 CONTINUE
IF( ignore ) THEN
  !"ignore" = if the output file exists, it must not be overwritten
  !Let's determine right now which files can be ignored and removed from the list
  !That will speed things up by a lot if many files already exist
  DO i=1,SIZE(filearray)
    !
    !Ignore the file unless we find at least one output file that doesn't exist
    ignorefile = .TRUE.
    !
    IF( filearray(i)(1:3).NE.'all' ) THEN
      !On a line there should be an input file name (mandatory)
      !and a format (optional), separated by a space
      strlength = SCAN(filearray(i)," ")
      inputfile = filearray(i)(1:strlength)
      !Check if the input file exists
      INQUIRE(FILE=inputfile,EXIST=fileexists)
      !
      IF( fileexists ) THEN
        !Check if the user specified a format for that peculiar file
        strlength = LEN_TRIM(inputfile)
        curr_outfileformat = filearray(i)(strlength+1:128)
        IF(LEN_TRIM(curr_outfileformat).NE.0) THEN
          CALL SET_OUTPUT(outfileformats,curr_outfileformat,.TRUE.)
        ENDIF
        !
        !Give a prefix to the output file
        strlength = SCAN(inputfile,'.',BACK=.TRUE.)
        IF(strlength.NE.0) THEN
          prefix = inputfile(1:strlength-1)
        ELSE
          prefix = inputfile
        ENDIF
        !
        !Check if output file(s) all exist
        DO j=1,SIZE(outfileformats)
          IF( LEN_TRIM(outfileformats(j)) > 0 ) THEN
            CALL NAME_OUTFILE(prefix,outputfile,outfileformats(j))
            INQUIRE(FILE=outputfile,EXIST=fileexists)
            IF( .NOT.fileexists ) THEN
              !This output file doesn't exist => this file will not be ignored
              ignorefile = .FALSE.
            ENDIF
          ENDIF
        ENDDO
        !
        IF(ignorefile) THEN
          !Remove this input file from the list
          filearray(i) = ""
          Nignored = Nignored+1
        ENDIF
        !
        !The format "curr_outfileformat" was only for the current file => disable it
        IF(LEN_TRIM(curr_outfileformat).NE.0) THEN
          CALL SET_OUTPUT(outfileformats,curr_outfileformat,.FALSE.)
        ENDIF
        !
      ENDIF
    ENDIF
    !
  ENDDO
ENDIF
!
!
!
300 CONTINUE
!Convert all files in the list
Nfiles=0
DO i=1,SIZE(filearray)
  outputfile = ""
  !Reset the file-specific output format
  curr_outfileformat = ""
  !If an output format applies to all files, (re)activate it
  IF(LEN_TRIM(list_outfileformat).NE.0) CALL SET_OUTPUT(outfileformats,list_outfileformat,.TRUE.)
  !
  !
  IF( LEN_TRIM(filearray(i)).NE.0 ) THEN
    !
    !On a line there should be an input file name (mandatory)
    !and a format (optional), separated by a space
    strlength = SCAN(filearray(i)," ")
    inputfile = filearray(i)(1:strlength)
    !Check if the input file exists
    INQUIRE(FILE=inputfile,EXIST=fileexists)
    IF(fileexists) THEN
      !Check if the user specified a format for that peculiar file
      strlength = LEN_TRIM(inputfile)
      curr_outfileformat = filearray(i)(strlength+1:128)
      IF(LEN_TRIM(curr_outfileformat).NE.0) THEN
        CALL SET_OUTPUT(outfileformats,curr_outfileformat,.TRUE.)
      ENDIF
      !
      !Give a prefix to the output file
      strlength = SCAN(inputfile,'.',BACK=.TRUE.)
      IF(strlength.NE.0) THEN
        outputfile = inputfile(1:strlength-1)
      ELSE
        outputfile = inputfile
      ENDIF
      !
      msg = 'inputfile, outputfile: '//TRIM(inputfile)//', '//TRIM(outputfile)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      !
      !Now we know the input file and the format to convert it to
      !=> convert the file
      CALL CONVERT_AFF(inputfile,options_array,outputfile,outfileformats)
      IF(nerr>0) GOTO 1000
      !
      !Increment the counter of converted files
      Nfiles = Nfiles+1
      !
      !The format "curr_outfileformat" was only for the current file => disable it
      IF(LEN_TRIM(curr_outfileformat).NE.0) THEN
        CALL SET_OUTPUT(outfileformats,curr_outfileformat,.FALSE.)
      ENDIF
      !
    ELSE
      !File does not exist => skip to next file
      nwarn=nwarn+1
      CALL ATOMSK_MSG(4700,(/TRIM(inputfile)/),(/0.d0/))
    ENDIF
    !
    !Message "file was converted"
    CALL ATOMSK_MSG(4001,(/''/),(/0.d0/))
    !
  ENDIF
  !
  !
ENDDO
335 CONTINUE
!
CALL ATOMSK_MSG(4045,(/''/),(/DBLE(Nfiles),DBLE(Nignored)/))
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  TERMINATE MODULE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
1000 CONTINUE
msg = 'LEAVING MODE_LIST...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
IF(ALLOCATED(filearray)) DEALLOCATE(filearray)
!
END SUBROUTINE LIST_XYZ
!
!
END MODULE mode_list
