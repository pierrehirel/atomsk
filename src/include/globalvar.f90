MODULE comv
!
!**********************************************************************************
!*  COMV                                                                          *
!**********************************************************************************
!* This module contains the global variables used by Atomsk.                      *
!* Global variables should be limited to the strict minimum, and NEVER            *
!* include variables or arrays about atomic systems (atom positions etc.).        *
!**********************************************************************************
!* (C) March 2010 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel (see date in the variable "version" below)         *
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
CHARACTER(LEN=24),PARAMETER:: version = 'master-2024-04-16'
!
!**********************************
!*  DATA TYPES / PRECISION
!**********************************
INTEGER,PARAMETER:: il = SELECTED_INT_KIND(9)        !integers up to 10^9
INTEGER,PARAMETER:: dp = SELECTED_REAL_KIND(15,307)  !reals with 64-bits precision
INTEGER(il),PARAMETER:: NATOMS_MAX = HUGE(INT(0,il)) !maximum number of atoms that Atomsk can handle
!
!**********************************
!*  ENVIRONMENT-DEPENDENT VARIABLES
!**********************************
CHARACTER(LEN=3):: lang !language in which the program will run (should be 2 letters)
CHARACTER(LEN=1):: langyes, langBigYes, langno !one-letter shortcuts for "yes" and "no", e.g. "y", "n"
#if defined(WINDOWS)
  CHARACTER(LEN=1),PARAMETER:: pathsep='\'     !path separator for Windows
  CHARACTER(LEN=3),PARAMETER:: system_ls="dir" !command to list current directory
  CHARACTER(LEN=16),PARAMETER:: pathnull="2>nul"  !redirection to NULL for Windows
  LOGICAL:: colourtext=.FALSE.  !by default, don't colour text in Windows
#else
  CHARACTER(LEN=1),PARAMETER:: pathsep='/'     !path separator for UNIX/Linux
  CHARACTER(LEN=3),PARAMETER:: system_ls="ls " !command to list current directory
  CHARACTER(LEN=16),PARAMETER:: pathnull="2>/dev/null"  !redirection to NULL for UNIX/Linux
  LOGICAL:: colourtext=.FALSE.   !by default, don't colour text in UNIX/Linux environments
#endif
!
!**********************************
!*  PROGRAM BEHAVIOR
!**********************************
CHARACTER(LEN=16):: neighsearch !algorithm for neighbor search (default: auto)
CHARACTER(LEN=128):: logfile !name of logfile for the program
LOGICAL:: overw, ignore      !automatically overwrite/ignore existing files?
INTEGER:: nwarn, nerr        !number of warnings/errors encountered during run
INTEGER:: verbosity          !level of verbosity of the program
!
!**********************************
!*  DISPLAY STYLES
!**********************************
CHARACTER(LEN=32):: colourdef="none"               !default colour for all messages
CHARACTER(LEN=32):: colourerr="red bold blink"     !default colour for ERROR
CHARACTER(LEN=32):: colourwarn="yellow bold blink" !default colour for WARNING
CHARACTER(LEN=32):: progressbar="linear"           !default style for progress bar
!
!**********************************
!*  INPUT/OUTPUT
!**********************************
!The following array contains a list of all formats (input+output) supported by Atomsk.
!First column gives file format extension or name,
!second column tells if format can be read by Atomsk ("yes" or "no"),
!third column if it can be written by Atomsk.
!Additional formats may be added here, don't forget to change array size.
!Note that each entry must be *exactly* 5 characters long (add spaces if necessary)
CHARACTER(LEN=5),DIMENSION(33,3),PARAMETER:: flist = RESHAPE( (/ &
  & "atsk ","yes  ","yes  ", &
  & "abin ","yes  ","yes  ", &
  & "bop  ","yes  ","yes  ", &
  & "bx   ","yes  ","yes  ", &
  & "cfg  ","yes  ","yes  ", &
  & "cel  ","yes  ","yes  ", &
  & "cif  ","yes  ","yes  ", &
  & "coo  ","yes  ","yes  ", &
  & "csv  ","yes  ","yes  ", &
  & "d12  ","yes  ","yes  ", &
  & "dat  ","no   ","yes  ", &
  & "dd   ","no   ","yes  ", &
  & "dlp  ","yes  ","yes  ", &
  & "fdf  ","yes  ","yes  ", &
  & "gin  ","yes  ","yes  ", &
  & "imd  ","yes  ","yes  ", &
  & "jems ","yes  ","yes  ", &
  & "lmc  ","yes  ","no   ", &
  & "lmp  ","yes  ","yes  ", &
  & "mol  ","yes  ","yes  ", &
  & "pdb  ","yes  ","yes  ", &
  & "pos  ","yes  ","yes  ", &
  & "pw   ","yes  ","yes  ", &
  & "pwo  ","yes  ","no   ", &
  & "str  ","yes  ","yes  ", &
  & "vout ","yes  ","no   ", &
  & "vesta","yes  ","yes  ", &
  & "xmd  ","yes  ","yes  ", &
  & "xsf  ","yes  ","yes  ", &
  & "xv   ","yes  ","yes  ", &
  & "xyz  ","yes  ","yes  ", &
  & "exyz ","yes  ","yes  ", &
  & "sxyz ","yes  ","yes  "  &
  & /) , SHAPE(flist) , ORDER=(/2,1/) )
!
!Unit for output files. Default is 40, can be changed to 6 for output to stdout
INTEGER:: ofu=40
!
END MODULE comv
