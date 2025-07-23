MODULE readconf
!
!**********************************************************************************
!*  READCONF                                                                      *
!**********************************************************************************
!* This module reads a configuration file containing default parameters           *
!* for the atomsk program. The absolute path to the file is passed through        *
!* the 'conffile' variable.                                                       *
!* The config files are read in the following order (see atomsk.f90):             *
!* - /etc/atomsk.conf (global config for all users)                               *
!* - /home/user/.config/atomsk.conf (user personal config)                        *
!* - current working directory (./atomsk.conf)                                    *
!* Any parameter found in a file overrides previous definitions.                  *
!**********************************************************************************
!* (C) Dec. 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 23 July 2025                                     *
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
USE strings
USE files
USE messages
!
!
!
CONTAINS
!
SUBROUTINE READ_CONF(conffile)
!
!
!Declare variables
IMPLICIT NONE
CHARACTER(LEN=5):: ext
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to write
CHARACTER(LEN=128):: conffile, msg, temp, keyword
INTEGER:: i, sp, Nthreads
!
!
Nthreads = 0
sp = 0
!
!
msg = "reading conffile = "//TRIM(conffile)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
OPEN(UNIT=19,FILE=conffile,FORM='FORMATTED',ERR=800)
!
DO
  READ(19,'(a128)',END=1000,ERR=800) temp
  temp = TRIM(ADJUSTL(temp))
  !
  IF( LEN_TRIM(temp)>0 .AND. temp(1:1).NE.'#' ) THEN
    !
    !Detect position of next space character
    sp = SCAN(temp," ")
    IF(sp<=0) sp=1
    !
    !Read first word ("keyword") from line
    READ(temp,*) keyword
    !Convert keyword into lower-case characters
    keyword = TRIM(ADJUSTL( StrDnCase(keyword) ))
    !
    SELECT CASE(keyword)
    CASE("verbosity")
      !read level of verbosity
      READ(temp(sp:),*,END=800,ERR=800) verbosity
    !
    CASE("format")
      !read formats that will always be activated for output
      READ(temp(sp:),*,END=800,ERR=800) ext
      CALL SET_OUTPUT(outfileformats,ext,.TRUE.)
    !
    CASE("ow","overw")
      !if true, file will always be overwritten without prompt
      READ(temp(sp:),*,END=800,ERR=800) temp
      CALL STR2BOOL(temp,overw)
    !
    CASE("ig","ignore")
      !if true, existing file will always be ignored without prompt
      READ(temp(sp:),*,END=800,ERR=800) temp
      CALL STR2BOOL(temp,ignore)
    !
    CASE("lang","language")
      !read language for messages
      READ(temp(sp:),*,END=800,ERR=800) msg
      msg = ADJUSTL(msg)
      IF( INDEX(msg,"fr")>0 .OR. INDEX(msg,"FR")>0 ) THEN
        lang="fr"
      ELSEIF( INDEX(msg,"de")>0 .OR. INDEX(msg,"DE")>0 ) THEN
        lang="de"
      ELSE
        lang="en"
      ENDIF
    !
    CASE("nthreads","threads")
      !read max. number of OpenMP threads to use
      READ(temp(9:),*,END=800,ERR=800) Nthreads
#if defined(OPENMP)
      IF( Nthreads>0 ) THEN
        CALL OMP_SET_NUM_THREADS(Nthreads)
      ENDIF
#else
      Nthreads=0
      CALL ATOMSK_MSG(751,(/conffile/),(/0.d0/))
#endif
    !
    CASE("neigh","neighlist","neighbor","neighbour","neigh_search","neighbor_search","neighbour_search")
      !user wants to use a specific search algorithm
      i = SCAN(temp," ")
      READ(temp(i+1:),*,END=800,ERR=800) neighsearch
    !
    CASE("colour","color")
      !if true, display some texts in colour
      IF( LEN_TRIM(temp(7:))>0 ) THEN
        READ(temp(7:),*,END=800,ERR=800) temp
        CALL STR2BOOL(temp,colourtext)
      ENDIF
    !
    CASE("colour_default","color_default")
      !change default text colour
      IF( colourtext ) THEN
        colourdef = TRIM(temp(15:))
      ENDIF
    !
    CASE("colour_warning","color_warning")
      !change warning text colour
      IF( colourtext ) THEN
        colourwarn = TRIM(temp(15:))
      ENDIF
    !
    CASE("colour_error","color_error")
      !change error text colour
      IF( colourtext ) THEN
        colourerr = TRIM(temp(13:))
      ENDIF
    !
    CASE("colour_prompt","color_prompt")
      !change colour of prompt in interactive mode
      IF( colourtext ) THEN
        colourprompt = TRIM(temp(13:))
      ENDIF
    !
    CASE("progressbar")
      !change style of progress bar: linear, dots, rotate, inflate, bounce
      i = SCAN(temp," ")
      READ(temp(i+1:),*,END=800,ERR=800) progressbar
    !
    !
    !Special keywords that will be treated by modes/options are ignored here
    CASE("density","nye","orthocell")
      CONTINUE
    !
    !
    !Unrecognized keywords will trigger a warning
    CASE DEFAULT
      nwarn = nwarn+1
      CALL ATOMSK_MSG(1702,(/temp,conffile/),(/0.d0/))
    !
    END SELECT
    !
  ENDIF
ENDDO
!
CLOSE(19)
GOTO 1000
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  ERROR MESSAGES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
800 CONTINUE
CALL ATOMSK_MSG(1703,(/conffile/),(/0.d0/))
nwarn = nwarn+1
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  TERMINATE MODULE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
1000 CONTINUE
!
!
END SUBROUTINE READ_CONF
!
!
END MODULE readconf
