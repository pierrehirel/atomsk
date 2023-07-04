MODULE display_messages
!
!**********************************************************************************
!*  DISPLAY_MESSAGES                                                              *
!**********************************************************************************
!* This module contains routines to display messages for the Atomsk program.      *
!**********************************************************************************
!* (C) June 2011 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 27 June 2023                                     *
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
!
USE comv
!
!
CONTAINS
!
!
!********************************************************
! DISPLAY_MSG
! This subroutine displays a message depending on the
! level of verbosity
!********************************************************
SUBROUTINE DISPLAY_MSG(verb,msg,logfile)
!
IMPLICIT NONE
CHARACTER(LEN=64):: logfile     !The log file to write messages into
CHARACTER(LEN=128):: msg, msg2  !The message itself
INTEGER:: verb   !Level of verbosity:
                 !0=silent, no message is ever displayed nor written in log file
                 !1=messages are written on the screen only
                 !2=messages are written in log file only
                 !3=messages are written in log+screen
                 !4=additional debugging messages (starting with 'debug')
                 !  are written in logfile
                 !The latter value (4) should be used only for
                 !programming/debugging purposes
!
msg2 = TRIM(ADJUSTL(msg))
!Write message in log file if verbosity>=2 or if debug message
IF( (verb>=2.AND.msg2(1:5).NE.'debug') .OR. verb==4 ) THEN
  OPEN(UNIT=20,FILE=logfile,FORM='FORMATTED',POSITION='APPEND')
  WRITE(20,*) TRIM(msg)
  CLOSE(20)
ENDIF
!Display message on screen if verbosity==1,3,4, or if error/warning message
IF( (verb==1 .OR. verb>=3).AND.msg2(1:5).NE.'debug' .OR. &
  & msg2(1:3)=='/!\' .OR. msg2(1:3)=='X!X' ) THEN
  WRITE(*,*) COLOUR_MSG(TRIM(msg),colourdef)
ENDIF
!
END SUBROUTINE DISPLAY_MSG
!
!
!
!********************************************************
! DISPLAY_COPYRIGHT
! This subroutine displays the copyright of atomsk
!********************************************************
SUBROUTINE DISPLAY_COPYRIGHT()
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
!
msg = "(C) P. Hirel 2010 - Version "//TRIM(ADJUSTL(version))
CALL DISPLAY_MSG(verbosity,msg,logfile)
!
END SUBROUTINE DISPLAY_COPYRIGHT
!
!
!********************************************************
! DISPLAY_HEADER
! This subroutine displays the "welcome message"
! of atomsk
!********************************************************
SUBROUTINE DISPLAY_HEADER()
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
!
!Print a nice message
msg = " ___________________________________________________"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "|              ___________                          |"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "|     o---o    "//COLOUR_MSG("A T O M S K","bold")
msg = TRIM(msg)//"                          |"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "|    o---o|    Version "//TRIM(ADJUSTL(version))
msg(53:53) = "|"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "|    |   |o    (C) 2010 Pierre Hirel                |"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "|    o---o     https://atomsk.univ-lille.fr         |"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "|___________________________________________________|"
CALL DISPLAY_MSG(verbosity,msg,logfile)
!
END SUBROUTINE DISPLAY_HEADER
!
!
!********************************************************
! DISPLAY_PROGBAR
! This subroutine displays a fancy progress bar.
! Available styles are:
!     linear   dots    rotate  inflate
!     jump     wave    snail
! Default style is "linear", it can be customized
! in config. file "atomsk.conf".
!********************************************************
SUBROUTINE DISPLAY_PROGBAR(r1,r2)
!
IMPLICIT NONE
CHARACTER(LEN=64):: bar, pc
CHARACTER(LEN=64):: msg, temp
INTEGER:: i, j
REAL(dp),INTENT(IN):: r1,r2
REAL(dp):: tempreal
!
bar = ""
pc = ""
msg = ""
temp = ""
tempreal = 100.d0*r1/r2 !percentage of progress
WRITE(pc,'(i3)') NINT(tempreal)
!
SELECT CASE(progressbar)
CASE("linear")
  !Display a linear progress bar (default)
  DO j=1,NINT(tempreal/5.d0)
    bar = TRIM(ADJUSTL(bar))//"="
  ENDDO
  bar = "["//TRIM(ADJUSTL(bar))//">"
  bar(22:22) = "]"
  msg = TRIM(ADJUSTL(bar))//" "//TRIM(ADJUSTL(pc))//"%"
  !
CASE("dots")
  !Display dots
  DO j=1,NINT(tempreal/5.d0)
    bar = TRIM(ADJUSTL(bar))//"."
  ENDDO
  msg = TRIM(ADJUSTL(bar))//"  ["//TRIM(ADJUSTL(pc))//"%]"
  !
CASE("rotate")
  !Imitate rotation
  j = MOD(NINT(tempreal),4)
  IF( j==1 ) THEN
    msg = "[\] "//TRIM(ADJUSTL(pc))//"%"
  ELSEIF( j==2 ) THEN
    msg = "[|] "//TRIM(ADJUSTL(pc))//"%"
  ELSEIF( j==3 ) THEN
    msg = "[/] "//TRIM(ADJUSTL(pc))//"%"
  ELSE
    msg = "[-] "//TRIM(ADJUSTL(pc))//"%"
  ENDIF
  !
CASE("clock")
  !Imitate a clock
  j = MOD(NINT(tempreal),6)
  IF( j==1 ) THEN
    msg = "(` ) "//TRIM(ADJUSTL(pc))//"%"
  ELSEIF( j==2 ) THEN
    msg = "( ´) "//TRIM(ADJUSTL(pc))//"%"
  ELSEIF( j==3 ) THEN
    msg = "( -) "//TRIM(ADJUSTL(pc))//"%"
  ELSEIF( j==4 ) THEN
    msg = "( .) "//TRIM(ADJUSTL(pc))//"%"
  ELSEIF( j==5 ) THEN
    msg = "(, ) "//TRIM(ADJUSTL(pc))//"%"
  ELSE
    msg = "(- ) "//TRIM(ADJUSTL(pc))//"%"
  ENDIF
  !
CASE("barrier")
  !Imitate a dot crossing a barrier
  j = MOD(NINT(tempreal),4)
  IF( j==1 ) THEN
    msg = "[ : ] "//TRIM(ADJUSTL(pc))//"%"
  ELSEIF( j==2 ) THEN
    msg = "[ ..] "//TRIM(ADJUSTL(pc))//"%"
  ELSEIF( j==3 ) THEN
    msg = "[ : ] "//TRIM(ADJUSTL(pc))//"%"
  ELSE
    msg = "[.. ] "//TRIM(ADJUSTL(pc))//"%"
  ENDIF
  !
CASE("inflate")
  !Display round characters of different sizes
  j = MOD(NINT(tempreal),4)
  IF( j==1 ) THEN
    msg = "[o] "//TRIM(ADJUSTL(pc))//"%"
  ELSEIF( j==2 ) THEN
    msg = "[O] "//TRIM(ADJUSTL(pc))//"%"
  ELSEIF( j==3 ) THEN
    msg = "[o] "//TRIM(ADJUSTL(pc))//"%"
  ELSE
    msg = "[.] "//TRIM(ADJUSTL(pc))//"%"
  ENDIF
  !
CASE("jump")
  !Imitate a jumping dot
  j = MOD(NINT(tempreal),4)
  IF( j==1 ) THEN
    msg = "[:] "//TRIM(ADJUSTL(pc))//"%"
  ELSEIF( j==2 ) THEN
    msg = "['] "//TRIM(ADJUSTL(pc))//"%"
  ELSEIF( j==3 ) THEN
    msg = "[:] "//TRIM(ADJUSTL(pc))//"%"
  ELSE
    msg = "[.] "//TRIM(ADJUSTL(pc))//"%"
  ENDIF
  !
CASE("wave")
  !Imitate a wave
  j = MOD(NINT(tempreal),9)+1
  temp = "..::!!!::...::!!!::...::!!!::...::!!!::.."
  bar = temp(j+1:j+20)
  msg = TRIM(ADJUSTL(bar))//"  ["//TRIM(ADJUSTL(pc))//"%]"
  !
CASE("bounce")
  !Display a bar bouncing from left to right
  j = 2*MIN( MOD(NINT(tempreal),20)+1 , 20-MOD(NINT(tempreal),20)-1 ) + 1
  bar(j:j+2) = "-=-"
  bar(1:1) = "["
  bar(22:23) = "] "
  msg = TRIM(ADJUSTL(bar))//" "//TRIM(ADJUSTL(pc))//"%"
  !
CASE("snail")
  !Display a snail
  j = NINT(tempreal/5.d0)
  bar(j:j+4) = ' _@" '
  msg = TRIM(bar)//"  ["//TRIM(ADJUSTL(pc))//"%]"
  !
CASE("pacman")
  !Display Pacman
  j = NINT(tempreal/5.d0)
  bar = " *   *   *   *   *   *"
  bar(1:j) = ""
  IF( bar(j+2:j+2)=="*" ) THEN
    bar(j:j+2) = ' (='
  ELSE
    bar(j:j+2) = ' (<'
  ENDIF
  msg = bar(1:22)//"  ["//TRIM(ADJUSTL(pc))//"%]"
  !
CASE("face")
  !Imitate a human face
  IF( NINT(tempreal)>90 ) THEN
    msg = "(^o^)  "//TRIM(ADJUSTL(pc))//"%"
  ELSE
    j = MOD(NINT(tempreal),5)
    IF( j==0 ) THEN
      msg = "(-.-)  "//TRIM(ADJUSTL(pc))//"%"
    ELSE
      msg = "(°.°)  "//TRIM(ADJUSTL(pc))//"%"
    ENDIF
  ENDIF
  !
CASE DEFAULT
  !Default style: only percentage is displayed
  msg = "["//TRIM(ADJUSTL(pc))//"%]"
END SELECT
!
!Display the progress bar
IF( NINT(r1) >= NINT(r2) .OR. NINT(tempreal)==100 ) THEN
  WRITE(*,'(a)',ADVANCE="NO") CHAR(13)
ELSE
  WRITE(*,'(a)',ADVANCE="NO") CHAR(13)//"     "//TRIM(msg)
ENDIF
!
END SUBROUTINE DISPLAY_PROGBAR
!
!
!********************************************************
! COLOUR_MSG
! This function colours a text. It is possible to
! combine a colour and a typesetting, e.g. "red bold".
! WARNING: this is expected to work in Linux/bash,
! but may not in other environments (macOS, Windows...)
!********************************************************
FUNCTION COLOUR_MSG(intxt,colour) RESULT(outtxt)
!
CHARACTER(LEN=*):: intxt
CHARACTER(LEN=LEN_TRIM(intxt)+16):: outtxt
CHARACTER(LEN=*):: colour
CHARACTER(LEN=16):: code
!
outtxt = intxt
code="[0"
!
!Apply colour only if global variable "colourtext" is true
IF( colourtext ) THEN
  !First, set colour / only first colour in the string "colour" is used
  IF( INDEX(colour,"black")>0 ) THEN
    code="[095"
  ELSEIF( INDEX(colour,"red")>0 ) THEN
    code="[031"
  ELSEIF( INDEX(colour,"green")>0 ) THEN
    code="[032"
  ELSEIF( INDEX(colour,"yellow")>0 ) THEN
    code="[033"
  ELSEIF( INDEX(colour,"blue")>0 ) THEN
    code="[034"
  ELSEIF( INDEX(colour,"magenta")>0 ) THEN
    code="[035"
  ELSEIF( INDEX(colour,"cyan")>0 ) THEN
    code="[036"
  ELSEIF( INDEX(colour,"grey")>0 .OR. INDEX(colour,"gray")>0 ) THEN
    code="[090"
  ELSEIF( INDEX(colour,"white")>0 ) THEN
    code="[097"
  ENDIF
  !
  !Add typesetting, several typesettings are possible
  IF( INDEX(colour,"bold")>0 .OR. INDEX(colour,"bright")>0 ) code=TRIM(code)//";1"
  IF( INDEX(colour,"italic")>0 ) code=TRIM(code)//";3"
  IF( INDEX(colour,"underline")>0 ) code=TRIM(code)//";4"
  IF( INDEX(colour,"blink")>0 ) code=TRIM(code)//";5"
  !
  IF( LEN_TRIM(code)>2 ) THEN
    code=TRIM(code)//"m"
    outtxt = ACHAR(27)//TRIM(ADJUSTL(code))//intxt//ACHAR(27)//"[0m"
  ENDIF
ENDIF
!
END FUNCTION COLOUR_MSG
!
!
!
END MODULE display_messages
