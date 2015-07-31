MODULE in_imd
!
!**********************************************************************************
!*  IN_IMD                                                                        *
!**********************************************************************************
!* This module reads IMD configuration files.                                     *
!* This format is described on this page:                                         *
!*   http://www.itap.physik.uni-stuttgart.de/~imd/userguide/config.html           *
!**********************************************************************************
!* (C) Feb. 2011 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 30 July 2015                                     *
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
!
CONTAINS
!
SUBROUTINE READ_IMD(inputfile,H,P,comment,AUXNAMES,AUX)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=128):: temp, temp2, temp3
CHARACTER(LEN=128):: msg
CHARACTER(LEN=1),DIMENSION(7):: F
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
LOGICAL:: arereduced  !are the coordinates reduced?
LOGICAL:: velocities  !are there velocities?
INTEGER:: i, id, itemp, j, NP
INTEGER:: Ncol   !number of columns for each atom
INTEGER:: typecol, xcol, ycol, zcol, idcol !column for species, X, Y, Z, id
INTEGER:: masscol !column for mass of atoms
INTEGER:: vx, vy, vz !columns for velocities
REAL(dp),DIMENSION(30):: acol  !all colmuns for each atom
                        !(we assume there are no more than 30 columns)
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
velocities = .FALSE.
i = 0
Ncol = 0
typecol = 0
idcol = 0
masscol = 0
xcol = 0
ycol = 0
zcol = 0
vx = 0
vy = 0
vz = 0
H(:,:) = 0.d0
acol(:) = 0.d0
ALLOCATE(comment(1))
 comment(1)=''
!
!
100 CONTINUE
msg = 'entering READ_IMD'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
!
140 CONTINUE
!Read the header
!Note: in principle this section is at the beginning of the file
!     (as suggests its name "header"), however IMD works just fine
!     if it is at the end. Therefore, here it is assumed that the
!     "header" can be at the beginning or at the end of the file
DO WHILE(temp.NE.'#E')
  READ(30,'(a128)',ERR=800,END=1000) temp
  temp = TRIM(ADJUSTL(temp))
  WRITE(msg,*) 'Header: ', TRIM(temp)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  IF(temp(1:2)=='#F') THEN
    !Format line
    READ(temp(3:128),*,END=800,ERR=800) (F(i), i=1,7)
    !Check if file is ASCII format
    IF(F(1).NE.'A') THEN
      nerr = nerr+1
      CALL ATOMSK_MSG(1807,(/''/),(/0.d0/))
      GOTO 1000
    ENDIF
    !Check if velocities are present
    IF(F(6)=='3') THEN
      velocities = .TRUE.
    ENDIF
  !
  ELSEIF(temp(1:2)=='#C') THEN
    !Contents line: find in which column each property is stored
    itemp = 2
    temp2 = TRIM(ADJUSTL(temp))
    DO i=1,SIZE(acol)
      temp2 = TRIM(ADJUSTL(temp2(itemp+1:128)))
      READ(temp2,*,END=150,ERR=150) temp3
      !No error? Then read what this column contains
      temp3 = TRIM(ADJUSTL(temp3))
      IF(temp3(1:4)=='type') THEN
        typecol = i
      ELSEIF(temp3(1:1)=='x' .AND. xcol==0) THEN
        xcol = i
      ELSEIF(temp3(1:1)=='y' .AND. ycol==0) THEN
        ycol = i
      ELSEIF(temp3(1:1)=='z' .AND. zcol==0) THEN
        zcol = i
      ELSEIF(temp3(1:4)=='mass' .AND. masscol==0) THEN
        masscol = i
      ELSEIF(temp3(1:6)=='number' .AND. idcol==0) THEN
        idcol = i
      ELSEIF(temp3(1:2)=='vx' .AND. vx==0) THEN
        vx = i
      ELSEIF(temp3(1:2)=='vy' .AND. vy==0) THEN
        vy = i
      ELSEIF(temp3(1:2)=='vz' .AND. vz==0) THEN
        vz = i
      ENDIF
      !Increase the number of columns that the file actually contains
      Ncol = Ncol+1
      itemp = LEN_TRIM(temp3)
    ENDDO
    !
    150 CONTINUE
    WRITE(msg,*) 'Ncol = ', Ncol
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,*) 'Col X, Y, Z, type: ', xcol, ycol, zcol, typecol
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !If some columns were not detected we are in trouble
    IF(xcol==0 .OR. ycol==0 .OR. zcol==0 .OR. typecol==0) THEN
      GOTO 800
    ENDIF
  !
  ELSEIF(temp(1:2)=='#X') THEN
    READ(temp(3:128),*,ERR=810,END=810) H(1,1), H(1,2), H(1,3)
  ELSEIF(temp(1:2)=='#Y') THEN
    READ(temp(3:128),*,ERR=810,END=810) H(2,1), H(2,2), H(2,3)
  ELSEIF(temp(1:2)=='#Z') THEN
    READ(temp(3:128),*,ERR=810,END=810) H(3,1), H(3,2), H(3,3)
  !
  ELSEIF(comment(1)=='') THEN
    READ(temp,'(a128)') comment(1)
  ENDIF
ENDDO
!
!
!
200 CONTINUE
!Find the beginning of atom positions
!Note: this is necessary because the header may appear before or after the atom positions
REWIND(30)
temp=''
i=-1
DO
  210 CONTINUE
  READ(30,'(a128)',ERR=800,END=1000) temp
  temp = TRIM(ADJUSTL(temp))
  !Try to read an integer, if it fails then go to next line
  READ(temp(1:1),*,ERR=210,END=210) i
  !Else, an integer was read: go back one line and start reading atom positions
  IF(i>=0) THEN
    BACKSPACE(30)
    GOTO 220
  ENDIF
ENDDO
220 CONTINUE
!Determine the number of atoms
NP=0
DO
  READ(30,'(a128)',END=250,ERR=250) temp3
  READ(temp3,*,END=250,ERR=250) (acol(j), j=1,Ncol)
  NP=NP+1
ENDDO
!
250 CONTINUE
WRITE(msg,*) 'NP = ', NP
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ALLOCATE(P(NP,4))
P(:,:) = 0.d0
!
!If velocities are present, allocate AUX
WRITE(msg,*) 'velocities = ', velocities
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
IF(velocities) THEN
  IF(masscol>0) THEN
    ALLOCATE( AUXNAMES(5) )
    ALLOCATE( AUX(NP,4) )
    AUXNAMES(5) = 'mass'
  ELSE
    ALLOCATE( AUXNAMES(4) )
    ALLOCATE( AUX(NP,4) )
  ENDIF
  AUXNAMES(2) = 'vx'
  AUXNAMES(3) = 'vy'
  AUXNAMES(4) = 'vz'
  AUX(:,:) = 0.d0
ELSE
  IF(masscol>0) THEN
    ALLOCATE( AUXNAMES(2) )
    ALLOCATE( AUX(NP,2) )
    AUXNAMES(2) = 'mass'
  ELSE
    ALLOCATE( AUXNAMES(1) )
    ALLOCATE( AUX(NP,1) )
  ENDIF
  AUX(:,:) = 0.d0
ENDIF
AUXNAMES(1) = 'type'
!
!
!
300 CONTINUE
!Go back to the beginning of atom positions
REWIND(30)
temp=''
i=-1
DO
  310 CONTINUE
  READ(30,'(a128)',ERR=800,END=1000) temp
  temp = TRIM(ADJUSTL(temp))
  !Try to read an integer, if it fails then go to next line
  READ(temp(1:1),*,ERR=310,END=310) i
  !Else, an integer was read: go back one line and start reading atom positions
  IF(i>=0) THEN
    BACKSPACE(30)
    GOTO 320
  ENDIF
ENDDO
!
320 CONTINUE
!Read each line and store positions and species to P
i=0
DO i=1,SIZE(P,1)
  READ(30,'(a128)',ERR=800,END=800) temp3
  READ(temp3,*,ERR=800,END=800) (acol(j), j=1,Ncol)
  id = NINT(acol(idcol))
  IF( id<=0 .OR. id>SIZE(P,1) ) THEN
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2742,(/""/),(/DBLE(id)/))
  ELSE
    P(id,1) = acol(xcol)
    P(id,2) = acol(ycol)
    IF(zcol.NE.0) P(id,3) = acol(zcol)
    !Types in IMD start at 0 => add 1
    P(id,4) = acol(typecol)+1
    AUX(id,1) = P(id,4)
    !
    IF(velocities) THEN
      !Read velocities
      AUX(id,2) = acol(vx)
      AUX(id,3) = acol(vy)
      AUX(id,4) = acol(vz)
      !Read mass of atom
      IF(masscol>0) THEN
        AUX(id,5) = acol(masscol)
      ENDIF
    ELSE
      !Read mass of atom
      IF(masscol>0) THEN
        AUX(id,2) = acol(masscol)
      ENDIF
    ENDIF
  ENDIF
  !
ENDDO
!
!Check if it is fractional coordinates
CALL FIND_IF_REDUCED(P,arereduced)
!If it is the case then convert to cartesian
IF(arereduced) THEN
  CALL FRAC2CART(P,H)
ENDIF
GOTO 1000
!

!
800 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
GOTO 1000
!
810 CONTINUE
CALL ATOMSK_MSG(1809,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
820 CONTINUE
CALL ATOMSK_MSG(1810,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
1000 CONTINUE
CLOSE(30)
!
!
END SUBROUTINE READ_IMD
!
!
END MODULE in_imd
