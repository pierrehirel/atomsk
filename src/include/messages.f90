MODULE messages
!
!**********************************************************************************
!*  MESSAGES                                                                      *
!**********************************************************************************
!* This module deals with messages displayed by the ATOMSK program.               *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 27 Oct. 2014                                     *
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
!* DISPLAY_COPYRIGHT   displays the copyright of the program                      *
!* DISPLAY_HEADER      displays the headline of the program                       *
!* DISPLAY_LICENSE     displays the license of the program                        *
!* DISPLAY_HELP        displays the help of the program                           *
!* ATOMSK_MSG          all messages used by atomsk                                *
!* DATE_MSG            displays a nice message according to the date              *
!**********************************************************************************
!
!
USE comv
USE display_messages
!Localization modules
USE messages_en  !English (default)
USE messages_fr  !French/francais
USE messages_de  !German/Deutsch
!--- add other language modules here ---
!
IMPLICIT NONE
!
!
CONTAINS
!
!
!********************************************************
! DISPLAY_LICENSE
! This subroutine displays the license message
! of atomsk
!********************************************************
SUBROUTINE DISPLAY_LICENSE()
!
SELECT CASE(lang)
!
CASE("de")
  CALL DISPLAY_LICENSE_DE()
CASE("fr")
  CALL DISPLAY_LICENSE_FR()
CASE DEFAULT
  CALL DISPLAY_LICENSE_EN()
END SELECT
!
END SUBROUTINE DISPLAY_LICENSE
!
!
!********************************************************
! DISPLAY_HELP
! This subroutine displays all or part of
! the help for atomsk
!********************************************************
SUBROUTINE DISPLAY_HELP(helpsection)
!
IMPLICIT NONE
CHARACTER(LEN=16):: helpsection
!
SELECT CASE(lang)
!
CASE("fr")
  CALL DISPLAY_HELP_FR(helpsection)
CASE("de")
  CALL DISPLAY_HELP_DE(helpsection)
CASE DEFAULT
  CALL DISPLAY_HELP_EN(helpsection)
END SELECT
!
!
END SUBROUTINE DISPLAY_HELP
!
!
!
!********************************************************
! ATOMSK_MSG
! This routine displays the message corresponding to
! the given index "imsg". The character array "strings"
! and the real array "reals" may be unallocated, or
! contain values relevent to the message.
!********************************************************
SUBROUTINE ATOMSK_MSG(imsg,strings,reals)
!
IMPLICIT NONE
CHARACTER(LEN=*),DIMENSION(:):: strings !Character strings that may be part of the message
INTEGER,INTENT(IN):: imsg  !index of message to display
INTEGER:: i
REAL(dp),DIMENSION(:),INTENT(IN):: reals  !real numbers that may be part of the message
!
!Strings: reduce the length if it is too long
!(this will mainly apply to file names)
DO i=1,SIZE(strings)
  IF( SCAN(strings(i),"/") > 0 ) THEN
    strings(i) = CHARLONG2SHRT(strings(i))
  ENDIF
ENDDO
!
SELECT CASE(lang)
!
CASE("de")
  langBigYes = "J"
  langyes = "j"
  langno = "n"
  CALL ATOMSK_MSG_DE(imsg,strings,reals)
CASE("fr")
  langBigYes = "O"
  langyes = "o"
  langno = "n"
  CALL ATOMSK_MSG_FR(imsg,strings,reals)
CASE DEFAULT
  langBigYes = "Y"
  langyes = "y"
  langno = "n"
  CALL ATOMSK_MSG_EN(imsg,strings,reals)
END SELECT
!
END SUBROUTINE
!
!
!
!********************************************************
! DATE_MSG
! Displays a nice message on certain dates.
!********************************************************
SUBROUTINE DATE_MSG()
!
SELECT CASE(lang)
!
CASE("de")
  CALL DATE_MSG_DE()
CASE("fr")
  CALL DATE_MSG_FR()
CASE DEFAULT
  CALL DATE_MSG_EN()
END SELECT
!
END SUBROUTINE DATE_MSG
!
!
!
!********************************************************
! CREATE_DATE
! Write the username and.
!********************************************************
SUBROUTINE CREATE_DATE(VALUES,username,msg)
!
IMPLICIT NONE
CHARACTER(LEN=128),INTENT(IN):: username
INTEGER,DIMENSION(8),INTENT(IN):: VALUES
CHARACTER(LEN=128),INTENT(OUT):: msg
!
SELECT CASE(lang)
!
CASE("de")
  CALL ATOMSK_CREATE_DATE_DE(VALUES,username,msg)
CASE("fr")
  CALL ATOMSK_CREATE_DATE_FR(VALUES,username,msg)
CASE DEFAULT
  CALL ATOMSK_CREATE_DATE(VALUES,username,msg)
END SELECT
!
END SUBROUTINE CREATE_DATE
!
!
!
END MODULE messages
