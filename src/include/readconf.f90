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
!* Last modification: P. Hirel - 06 June 2023                                     *
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
CHARACTER(LEN=128):: conffile, msg, temp
INTEGER:: i, Nthreads
!
!
Nthreads = 0
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
    IF( StrDnCase(temp(1:9))=='verbosity' ) THEN
      !read level of verbosity
      READ(temp(10:),*,END=800,ERR=800) verbosity
    !
    ELSEIF( StrDnCase(temp(1:6))=='format') THEN
      !read formats that will always be activated for output
      READ(temp(7:),*,END=800,ERR=800) ext
      CALL SET_OUTPUT(outfileformats,ext,.TRUE.)
    !
    ELSEIF( StrDnCase(temp(1:5))=='overw'  ) THEN
      !if true, file will always be overwritten without prompt
      READ(temp(6:),*,END=800,ERR=800) temp
      CALL STR2BOOL(temp,overw)
    !
    ELSEIF( StrDnCase(temp(1:6))=='ignore' ) THEN
      !if true, existing file will always be ignored without prompt
      READ(temp(7:),*,END=800,ERR=800) temp
      CALL STR2BOOL(temp,ignore)
    !
    ELSEIF( StrDnCase(temp(1:4))=='lang' ) THEN
      !read language for messages
      READ(temp(5:),*,END=800,ERR=800) msg
      msg = ADJUSTL(msg)
      IF( INDEX(msg,"fr")>0 .OR. INDEX(msg,"FR")>0 ) THEN
        lang="fr"
      ELSEIF( INDEX(msg,"de")>0 .OR. INDEX(msg,"DE")>0 ) THEN
        lang="de"
      ELSE
        lang="en"
      ENDIF
    !
    ELSEIF( StrDnCase(temp(1:8))=='nthreads' ) THEN
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
    ELSEIF( INDEX(StrDnCase(temp),"neigh") == 1 ) THEN
      !user wants to use a specific search algorithm
      i = SCAN(temp," ")
      READ(temp(i+1:),*,END=800,ERR=800) neighsearch
    !
    ELSEIF( StrDnCase(temp(1:7))=='colour ' ) THEN
      !if true, display some texts in colour
      IF( LEN_TRIM(temp(7:))>0 ) THEN
        READ(temp(7:),*,END=800,ERR=800) temp
        CALL STR2BOOL(temp,colourtext)
      ENDIF
    !
    ELSEIF( StrDnCase(temp(1:14))=='colour_default' ) THEN
      !change default text colour
      IF( colourtext ) THEN
        colourdef = TRIM(temp(15:))
      ENDIF
    !
    ELSEIF( StrDnCase(temp(1:14))=='colour_warning' ) THEN
      !change warning text colour
      IF( colourtext ) THEN
        colourwarn = TRIM(temp(15:))
      ENDIF
    !
    ELSEIF( StrDnCase(temp(1:12))=='colour_error' ) THEN
      !change error text colour
      IF( colourtext ) THEN
        colourerr = TRIM(temp(13:))
      ENDIF
    !
    !
    !Special keywords that will be treated by modes are ignored here
    ELSEIF( StrDnCase(temp(1:3))=='nye' ) THEN
      CONTINUE
    !
    !
    !Unrecognized keywords will trigger a warning
    ELSE
      nwarn = nwarn+1
      CALL ATOMSK_MSG(1702,(/temp,conffile/),(/0.d0/))
    ENDIF
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
