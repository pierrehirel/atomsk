MODULE out_csv 
!
!**********************************************************************************
!*  OUT_CSV                                                                       *
!**********************************************************************************
!* This module writes a file in Comma-Separated Values (CSV) format.              *
!* The CSV format has no standard specification. A description is given here:     *
!*    https://tools.ietf.org/html/rfc4180                                         *
!**********************************************************************************
!* (C) Oct. 2018 - Pierre Hirel                                                   *
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
!
!
USE comv
USE constants
USE messages
USE files
USE subroutines
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE WRITE_CSV(H,P,S,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=1):: fs=','  !field separator. Default is a comma (,). Change it here if you want something else
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: v1,v2,v3  !to store values
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
INTEGER:: i, j, k
INTEGER:: Ncol  !Number of columns to write
INTEGER:: Hcol, commentcol !column numbers for box vectors H and comments
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: S
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
Hcol = 0
 commentcol=0
!
!
100 CONTINUE
msg = 'entering WRITE_CSV'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
IF(ofu.NE.6) THEN
  OPEN(UNIT=ofu,FILE=outputfile,STATUS='UNKNOWN',ERR=250)
ENDIF
!
! First line is header to indicate the fields names
msg = 'species'//fs//'x'//fs//'y'//fs//'z'
Ncol = 4
IF( ALLOCATED(S) ) THEN
  msg = TRIM(ADJUSTL(msg))//fs//'S'//fs//'sx'//fs//'sy'//fs//'sz'
  Ncol=Ncol+4
ENDIF
IF( ALLOCATED(AUXNAMES) .AND. SIZE(AUXNAMES)>0 ) THEN
  DO i=1,SIZE(AUXNAMES)
    msg = TRIM(ADJUSTL(msg))//fs//TRIM(ADJUSTL(AUXNAMES(i)))
    Ncol=Ncol+1
  ENDDO
ENDIF
!
!Add columns for the box vectors
IF( SIZE(P,1)>=9 ) THEN
  !Store all values H(:,:) in one column
  msg = TRIM(ADJUSTL(msg))//fs//'H'
  Ncol=Ncol+1
  Hcol=-1*Ncol
ELSE
  !Use 9 columns to store box vectors
  msg = TRIM(ADJUSTL(msg))//fs//'H11'//fs//'H12'//fs//'H13'//fs// &
      &   'H21'//fs//'H22'//fs//'H23'//fs//'H31'//fs//'H32'//fs//'H33'
  Hcol=Ncol+1
  Ncol=Ncol+9
ENDIF
!
!Add one column for comments (if necessary)
!NOTE: if there are more comment lines than atoms, some comments will be lost
IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
  msg = TRIM(ADJUSTL(msg))//fs//'comment'
  Ncol=Ncol+1
ENDIF
!Write the line to the file
WRITE(ofu,'(a)') TRIM(ADJUSTL(msg))
!
!
! Write coordinates in CSV file
DO i=1,SIZE(P,1)
  !
  !Get species of current atom
  CALL ATOMSPECIES(P(i,4),species)
  !Get position of current atom
  WRITE(v1,'(f16.8)') P(i,1)
  WRITE(v2,'(f16.8)') P(i,2)
  WRITE(v3,'(f16.8)') P(i,3)
  !Write that into the line
  msg = TRIM(species)//fs//TRIM(ADJUSTL(v1))//fs//TRIM(ADJUSTL(v2))//fs//TRIM(ADJUSTL(v3))
  !
  IF( ALLOCATED(S) ) THEN
    !Get species of current atom
    CALL ATOMSPECIES(S(i,4),species)
    !Get position of current shell
    WRITE(v1,'(f16.8)') S(i,1)
    WRITE(v2,'(f16.8)') S(i,2)
    WRITE(v3,'(f16.8)') S(i,3)
    !Write that into the line
    msg = TRIM(ADJUSTL(msg))//fs//TRIM(ADJUSTL(species))//fs// &
        & TRIM(ADJUSTL(v1))//fs//TRIM(ADJUSTL(v2))//fs//TRIM(ADJUSTL(v3))
  ENDIF
  !
  IF( ALLOCATED(AUX) ) THEN
    DO j=1,SIZE(AUX,2)
      WRITE(v1,'(f16.8)') AUX(i,j)
      !Write that into the line
      msg = TRIM(ADJUSTL(msg))//fs//TRIM(ADJUSTL(v1))
    ENDDO
  ENDIF
  !
  IF( Hcol<0 ) THEN
    !All values of H are stored in same column
    IF( i<=3 ) THEN
      !Write components of first box vector H(1,:)
      WRITE(v1,'(f16.8)') H(1,i)
      msg = TRIM(ADJUSTL(msg))//fs//TRIM(ADJUSTL(v1))
    ELSEIF( i<=6 ) THEN
      !Write components of second box vector H(2,:)
      WRITE(v1,'(f16.8)') H(2,i-3)
      msg = TRIM(ADJUSTL(msg))//fs//TRIM(ADJUSTL(v1))
    ELSEIF( i<=9 ) THEN
      !Write components of third box vector H(3,:)
      WRITE(v1,'(f16.8)') H(3,i-6)
      msg = TRIM(ADJUSTL(msg))//fs//TRIM(ADJUSTL(v1))
    ELSE
      !No more values to write: leave a blanck field
      msg = TRIM(ADJUSTL(msg))//fs
    ENDIF
  ELSE
    !All values of H(:,:) are stored in first line
    IF( i==1 ) THEN
      !This is the first line: write values of H
      WRITE(v1,'(f16.8)') H(1,1)
      WRITE(v2,'(f16.8)') H(1,2)
      WRITE(v3,'(f16.8)') H(1,3)
      msg = TRIM(msg)//fs//TRIM(ADJUSTL(v1))//fs//TRIM(ADJUSTL(v2))//fs//TRIM(ADJUSTL(v3))
      WRITE(v1,'(f16.8)') H(2,1)
      WRITE(v2,'(f16.8)') H(2,2)
      WRITE(v3,'(f16.8)') H(2,3)
      msg = TRIM(msg)//fs//TRIM(ADJUSTL(v1))//fs//TRIM(ADJUSTL(v2))//fs//TRIM(ADJUSTL(v3))
      WRITE(v1,'(f16.8)') H(3,1)
      WRITE(v2,'(f16.8)') H(3,2)
      WRITE(v3,'(f16.8)') H(3,3)
      msg = TRIM(msg)//fs//TRIM(ADJUSTL(v1))//fs//TRIM(ADJUSTL(v2))//fs//TRIM(ADJUSTL(v3))
    ELSE
      !This is another line: just write empty fields
      DO j=1,9
        msg = TRIM(ADJUSTL(msg))//fs
      ENDDO
    ENDIF
  ENDIF
  !
  IF( SIZE(comment)>0 ) THEN
    !Store comments in corresponding column
    IF( i<=SIZE(comment) ) THEN
      k = SCAN(comment(i),fs)
      IF( k>0 ) THEN
        !There is one or several commas in this comment => use quotes
        msg = TRIM(msg)//fs//'"'//TRIM(ADJUSTL(comment(i)))//'"'
      ELSE
        !No commas => it is safe to just write the comment
        msg = TRIM(msg)//fs//TRIM(ADJUSTL(comment(i)))
      ENDIF
    ELSE
      msg = TRIM(ADJUSTL(msg))//fs
    ENDIF
  ENDIF
  !
  !Write the line to the file
  WRITE(ofu,'(a)') TRIM(ADJUSTL(msg))
END DO
GOTO 300
!
250 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
300 CONTINUE
msg = "CSV"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
1000 CONTINUE
IF(ofu.NE.6) THEN
  CLOSE(ofu)
ENDIF
!
END SUBROUTINE WRITE_CSV
!
!
END MODULE out_csv
