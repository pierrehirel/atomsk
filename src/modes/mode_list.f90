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
!* Last modification: P. Hirel - 22 Oct. 2013                                     *
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
CHARACTER(LEN=3):: infileformat
CHARACTER(LEN=5):: curr_outfileformat !output file format for one file
CHARACTER(LEN=5):: list_outfileformat !output file format that applies to the whole list
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of output file formats
CHARACTER(LEN=128):: temp
CHARACTER(LEN=128):: msg
CHARACTER(LEN=4096):: inputfile, outputfile
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
LOGICAL:: fileexists
INTEGER:: Nfiles
INTEGER:: strlength
REAL(dp),DIMENSION(3,3):: ORIENT  !crystal orientation
!
!Initialize
Nfiles = 0
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
OPEN(UNIT=50,FILE=filelist,STATUS='OLD',FORM='FORMATTED')
DO
  outputfile = ""
  !Reset the file-specific output format
  curr_outfileformat = ""
  !If an output format applies to all files, (re)activate it
  IF(LEN_TRIM(list_outfileformat).NE.0) CALL SET_OUTPUT(outfileformats,list_outfileformat,.TRUE.)
  !
  !Read the whole line
  READ(50,'(a128)',ERR=235,END=235) temp
  temp = ADJUSTL(temp)
  !
  IF(LEN_TRIM(temp).NE.0 .AND. temp(1:1).NE.'#') THEN
    IF(temp(1:3)=='all') THEN
      !'all <format>' means that all the following files have to be converted
      !into the <format> by default
      READ(temp(4:20),*) list_outfileformat
      CALL SET_OUTPUT(outfileformats,list_outfileformat,.TRUE.)
      CALL ATOMSK_MSG(4002,(/TRIM(list_outfileformat)/),(/0.d0/))
      msg = 'list_outfileformat: '//TRIM(list_outfileformat)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      !
    ELSE
      !On a line there should be an input file name (mandatory)
      !and a format (optional), separated by a space
      strlength = SCAN(temp," ")
      inputfile = temp(1:strlength)
      !Check if the input file exists
      INQUIRE(FILE=inputfile,EXIST=fileexists)
      IF(fileexists) THEN
        !Check if the user specified a format for that peculiar file
        strlength = LEN_TRIM(inputfile)
        READ(temp(strlength+1:128),*,ERR=231,END=231) curr_outfileformat
        231 CONTINUE
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
        !
        !Increment the counter of converted files
        Nfiles = Nfiles+1
        !
        !The format "outfileformat" was only for the current file => disable it
        IF(LEN_TRIM(curr_outfileformat).NE.0) THEN
          CALL SET_OUTPUT(outfileformats,curr_outfileformat,.FALSE.)
        ENDIF
        !
      ELSE
        !File does not exist => skip to next file
        nwarn=nwarn+1
        CALL ATOMSK_MSG(4700,(/TRIM(inputfile)/),(/0.d0/))
      ENDIF
    ENDIF
    CALL ATOMSK_MSG(4001,(/''/),(/0.d0/))
  ENDIF
ENDDO
235 CONTINUE
CLOSE(50)
!
CALL ATOMSK_MSG(4045,(/''/),(/DBLE(Nfiles)/))
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  TERMINATE MODULE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
1000 CONTINUE
msg = 'LEAVING MODE_LIST...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
END SUBROUTINE LIST_XYZ
!
!
END MODULE mode_list
