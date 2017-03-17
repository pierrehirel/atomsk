MODULE mode_normal
!
!**********************************************************************************
!*  MODE_NORMAL                                                                   *
!**********************************************************************************
!* This module converts a single file to one or several other file(s).            *
!* It does that in three steps:                                                   *
!* 1- read the input file and load information into memory                        *
!* 2- apply options to the system (deform, expand...)                             *
!* 3- write information to file(s) in one or several format(s)                    *
!*                                                                                *
!* WHAT MUST BE PROVIDED:                                                         *
!* inputfile                                                                      *
!*           string containing the name of the input file to read.                *
!* options_array                                                                  *
!*           if allocated, string array of unknown dimension containing the       *
!*           names and parameters of options to be applied to the system.         *
!*           If not allocated, then no option will be applied.                    *
!* outputfile                                                                     *
!*           string containing the name for output file(s). The extensions        *
!*           (e.g. ".cfg", ".xyz", etc.) will be appened to this prefix           *
!*           accordingly for each written file.                                   *
!* outfileformats                                                                 *
!*           string array of unknown dimension, containing the extensions         *
!*           (e.g. "cfg", "xyz", etc.) corresponding to file formats that         *
!*           will be written. If not allocated then no file will be written.      *
!*                                                                                *
!* WHAT IS RETURNED BY THIS ROUTINE:                                              *
!* no variable nor array, only files are written on the disk.                     *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 16 Dec. 2015                                     *
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
USE guess_form
USE deterH
!
!Module for reading input files
USE readin
!
!Module for options
USE options
!
!Module for writing output files
USE writeout
!
!
CONTAINS
!
SUBROUTINE CONVERT_AFF(inputfile,options_array,outputfile,outfileformats)
!
!Declare variables
IMPLICIT NONE
CHARACTER(LEN=3):: infileformat
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats
CHARACTER(LEN=128):: temp
CHARACTER(LEN=128):: msg
CHARACTER(LEN=4096):: inputfile, outputfile
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment, tempcomment
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: strlength
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystal orientation
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, S !positions of cores, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!
msg = 'ENTERING CONVERT_AFF...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Initialize variables
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
H(:,:) = 0.d0
ORIENT(:,:) = 0.d0
infileformat=''
strlength=0
ALLOCATE(tempcomment(1))
!
!If no output file was specified, prompt the user
IF(LEN_TRIM(outputfile)==0) THEN
  strlength = SCAN(inputfile,'.',BACK=.TRUE.)
  outputfile = inputfile(1:strlength)
ENDIF
!
!If the interactive mode runs, then we must use the
!'outfileformat' defined previously
!Else, the output file formats are defined by the
!global logical variables (output_xyz, output_cfg, etc.)
!so we have to reset the 'outfileformat' variable
!IF(mode.NE.'interactive') outfileformat = ''
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  READ INPUT FILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
100 CONTINUE
CALL READ_AFF(inputfile,H,P,S,comment,AUXNAMES,AUX)
!
IF(verbosity==4) THEN
  tempcomment(1) = '#Debug file for atomsk, from input file: '//TRIM(inputfile)
  temp = 'atomsk.xyz'
  CALL WRITE_XYZ(H,P,tempcomment,AUXNAMES,AUX,temp,'sxyz ')
ENDIF
IF(nerr>0 .OR. .NOT.ALLOCATED(P)) GOTO 1000
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  DETERMINE H
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
300 CONTINUE
!!Determine if the cell vectors were found or not
IF( VECLENGTH(H(1,:))==0.d0 .OR. VECLENGTH(H(2,:))==0.d0 .OR. VECLENGTH(H(3,:))==0.d0 ) THEN
  !!If basis vectors H were not found in the input file, then compute them
  ! (this does not work perfectly yet, see determine_H.f90)
  CALL ATOMSK_MSG(4701,(/''/),(/0.d0/))
  !
  CALL DETERMINE_H(H,P)
  !
ENDIF
!
IF(verbosity==4) THEN
  msg = 'Base vectors:'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,310) '     | ', H(1,1), H(1,2), H(1,3), '|'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,310) ' H = | ', H(2,1), H(2,2), H(2,3), '|'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,310) '     | ', H(3,1), H(3,2), H(3,3), '|'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  310 FORMAT(a7,3(f10.5,2X),a1)
ENDIF
!
IF(nerr>=1 .OR. .NOT.ALLOCATED(P)) GOTO 810
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  APPLY OPTIONS TO THE SYSTEM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
400 CONTINUE
msg = 'APPLYING OPTIONS TO THE SYSTEM:'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
IF( ALLOCATED(options_array) ) THEN
  CALL OPTIONS_AFF(options_array,H,P,S,AUXNAMES,AUX,ORIENT,SELECT)
ENDIF
!
450 CONTINUE
IF(nerr>0 .OR. .NOT.ALLOCATED(P)) GOTO 820
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  WRITE OUTPUT FILE(S)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
500 CONTINUE
CALL WRITE_AFF(outputfile,outfileformats,H,P,S,comment,AUXNAMES,AUX)
!
IF(verbosity==4) THEN
  tempcomment(1) = '#Debug file for atomsk'
  temp = 'atomsk.xyz'
  CALL WRITE_XYZ(H,P,tempcomment,AUXNAMES,AUX,temp,'sxyz ')
ENDIF
!
IF(nerr>0) THEN
  GOTO 830
ELSE
  GOTO 1000
ENDIF
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  ERROR MESSAGES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
800 CONTINUE
GOTO 1000
!
810 CONTINUE
CALL ATOMSK_MSG(2803,(/''/),(/0.d0/))
GOTO 1000
!
820 CONTINUE
CALL ATOMSK_MSG(2804,(/''/),(/0.d0/))
GOTO 1000
!
830 CONTINUE
CALL ATOMSK_MSG(3802,(/''/),(/0.d0/))
GOTO 1000
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  TERMINATE MODULE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
1000 CONTINUE
IF(ALLOCATED(P)) DEALLOCATE(P)
IF(ALLOCATED(S)) DEALLOCATE(S)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
IF(ALLOCATED(comment)) DEALLOCATE(comment)
IF(ALLOCATED(tempcomment)) DEALLOCATE(tempcomment)
msg = 'LEAVING CONVERT_AFF...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
END SUBROUTINE CONVERT_AFF
!
!
END MODULE mode_normal
