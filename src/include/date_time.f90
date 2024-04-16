MODULE date_time
!
!**********************************************************************************
!*  DATE_TIME                                                                     *
!**********************************************************************************
!* This module contains functions manipulating date and time.                     *
!**********************************************************************************
!* (C) April 2024 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 16 April 2024                                    *
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
!* List of functions in this file:                                                *
!* DAY_OF_WEEK         finds the day of the week given the date                   *
!* INT2MONTH           transforms an integer into a month                         *
!* INT2DAY             transforms an integer into a day                           *
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
!  DAY_OF_WEEK
!  This function determines the day of the week
!  (Monday, Tuesday...) from the date.
!********************************************************
FUNCTION DAY_OF_WEEK(day,month,year)
!
IMPLICIT NONE
INTEGER:: Day_of_week, j, k, mm, yy
INTEGER,INTENT(IN) :: day, month, year
!
mm=month
yy=year
IF(mm.le.2) THEN
  mm=mm+12
  yy=yy-1
END IF
j = yy / 100
k = MOD(yy, 100)
Day_of_week = MOD( DBLE(day) + ((DBLE(mm)+1.d0)*26.d0)/10.d0 + DBLE(k) + DBLE(k)/4.d0 + DBLE(j)/4.d0 + 5.d0*DBLE(j), 7.d0)
!
END FUNCTION Day_of_week
!
!
!
!********************************************************
! INT2MONTH
! Converts an integer to strings containing the
! month in long and short forms.
! Integer should be between 1 and 12
! (otherwise "December" is returned).
!********************************************************
SUBROUTINE INT2MONTH(number,month,smonth)
CHARACTER(LEN=16),INTENT(OUT):: month, smonth
INTEGER,INTENT(IN):: number
!
SELECT CASE(number)
CASE(1)
  month = "January"
  smonth = "Jan."
CASE(2)
  month = "February"
  smonth = "Feb."
CASE(3)
  month = "March"
  smonth = "March"
CASE(4)
  month = "April"
  smonth = "Apr."
CASE(5)
  month = "May"
  smonth = "May"
CASE(6)
  month = "June"
  smonth = "June"
CASE(7)
  month = "July"
  smonth = "July"
CASE(8)
  month = "August"
  smonth = "Aug."
CASE(9)
  month = "September"
  smonth = "Sept."
CASE(10)
  month = "October"
  smonth = "Oct."
CASE(11)
  month = "November"
  smonth = "Nov."
CASE(12)
  month = "December"
  smonth = "Dec."
CASE DEFAULT
  month = "Unknown"
  smonth = "N/A"
END SELECT
!
END SUBROUTINE INT2MONTH
!
!
!********************************************************
! INT2DAY
! Converts an integer to strings containing the
! day in long and short forms.
! Integer should be between 1 and 12
! (otherwise "Sunday" is returned).
!********************************************************
SUBROUTINE INT2DAY(number,day,sday)
CHARACTER(LEN=16),INTENT(OUT):: day,sday
INTEGER,INTENT(IN):: number
!
SELECT CASE(number)
CASE(0,7)
  day = "Sunday"
  sday = "Sun."
CASE(1)
  day = "Monday"
  sday = "Mon."
CASE(2)
  day = "Tuesday"
  sday = "Tue."
CASE(3)
  day = "Wednesday"
  sday = "Wed."
CASE(4)
  day = "Thursday"
  sday = "Thu."
CASE(5)
  day = "Friday"
  sday = "Fri."
CASE(6)
  day = "Saturday"
  sday = "Sat."
CASE DEFAULT
  day = "Unknown"
  sday = "N/A"
END SELECT
!
END SUBROUTINE INT2DAY
!
!
!
END MODULE date_time
