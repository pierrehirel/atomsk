MODULE in_xmd
!
!**********************************************************************************
!*  IN_XMD                                                                        *
!**********************************************************************************
!* This module reads XMD configuration files.                                     *
!* This format is described on this page:                                         *
!*   http://xmd.sourceforge.net/doc/manual/xmd.html                               *
!**********************************************************************************
!* (C) Nov. 2013 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 08 Jan. 2020                                     *
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
SUBROUTINE READ_XMD(inputfile,H,P,comment,AUXNAMES,AUX)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: temp, temp2
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
LOGICAL:: arereduced  !are the coordinates reduced?
LOGICAL:: velocities  !are there velocities?
INTEGER:: i, j, k, NP
INTEGER:: Ncol   !number of columns for each atom
INTEGER:: strlength
INTEGER:: typecol, vx, vy, vz !columns for atom types, velocities
INTEGER,DIMENSION(30,2):: atypes !relationship between atom types and atom species
                                 !(assuming there are no more than 30 different species in the system)
REAL(dp):: snumber  !atomic number
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
arereduced = .FALSE.
velocities = .FALSE.
i = 0
Ncol = 0
typecol = 0
vx = 0
vy = 0
vz = 0
atypes(:,:) = 0
H(:,:) = 0.d0
ALLOCATE(comment(1))
 comment(1)=''
!
!
100 CONTINUE
msg = 'entering READ_XMD'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
!Parse the file a first time to count the number of auxiliary properties
Ncol=1 !at least 1 column to store atom types
DO
  READ(30,'(a128)',ERR=200,END=200) temp
  temp = TRIM(ADJUSTL(temp))
  IF( temp(1:6)=="POSVEL" ) THEN
    velocities = .TRUE.
    Ncol = Ncol+3
  ENDIF
ENDDO
!
!
!
200 CONTINUE
ALLOCATE(AUXNAMES(Ncol))
typecol=1
AUXNAMES(typecol) = "type"
IF( velocities ) THEN
  vx = typecol+1
  vy=vx+1
  vz=vy+1
  AUXNAMES(vx) = "vx"
  AUXNAMES(vy) = "vy"
  AUXNAMES(vz) = "vz"
ENDIF
!
REWIND(30)
DO
  READ(30,'(a128)',ERR=800,END=300) temp
  temp = TRIM(ADJUSTL(temp))
  !
  IF( temp(1:1).NE.'#' ) THEN
    !
    IF( temp(1:3)=="BOX" ) THEN
      H(:,:) = 0.d0
      temp2 = ADJUSTL(temp(4:))
      IF( temp2(1:5)=="SCALE" ) THEN
        arereduced = .TRUE.
        READ(temp2(6:),*,ERR=800,END=800) H(1,1), H(2,2), H(3,3)
      ELSE
        arereduced = .FALSE.
        READ(temp2,*,ERR=800,END=800) H(1,1), H(2,2), H(3,3)
      ENDIF
      !
      !
    ELSEIF( temp(1:8)=="POSITION" .OR. temp(1:8)=="PARTICLE" .OR. temp(1:6)=="POSVEL" ) THEN
      IF( temp(1:8)=="POSITION" .OR. temp(1:8)=="PARTICLE" ) THEN
        READ(temp(9:),*,ERR=800,END=800) snumber
        velocities = .FALSE.
      ELSE
        READ(temp(7:),*,ERR=800,END=800) snumber
        velocities = .TRUE.
      ENDIF
      ! 
      IF( snumber > NATOMS_MAX ) THEN
        nerr = nerr+1
        CALL ATOMSK_MSG(821,(/""/),(/snumber/))
        GOTO 1000
      ENDIF
      NP = NINT(snumber)
      !
      ALLOCATE(P(NP,4) , STAT=i)
      IF( i>0 ) THEN
        ! Allocation failed (not enough memory)
        nerr = nerr+1
        CALL ATOMSK_MSG(819,(/''/),(/DBLE(NP)/))
        GOTO 1000
      ENDIF
      P(:,:) = 0.d0
      ALLOCATE(AUX(NP,SIZE(AUXNAMES)))
      AUX(:,:) = 0.d0
      !
      IF( velocities ) THEN
        DO i=1,NP
          READ(30,*,ERR=800,END=800) AUX(i,typecol), P(i,1), P(i,2), P(i,3), AUX(i,vx), AUX(i,vy), AUX(i,vz)
          P(i,4) = AUX(i,typecol)
        ENDDO
      ELSE
        DO i=1,NP
          READ(30,*,ERR=800,END=800) AUX(i,typecol), P(i,1), P(i,2), P(i,3)
          P(i,4) = AUX(i,typecol)
        ENDDO
      ENDIF
      !
      !If coordinates are reduced then convert to cartesian
      IF(arereduced) THEN
        CALL FRAC2CART(P,H)
      ENDIF
      !
      !
    ELSEIF( temp(1:9)=="TYPENAME " ) THEN
      !Save atom types in atypes(:,1) and atom species in atypes(:,2)
      !Note: at this point atom positions and types may not have been read yet,
      !     and P may not be allocated, therefore atomic numbers cannot be saved in P(:,4).
      temp2 = ADJUSTL(temp(10:))
      DO WHILE( LEN_TRIM(temp2).NE.0 )
        READ(temp2,*,ERR=800,END=800) k
        strlength = SCAN(temp2," ")
        temp2 = ADJUSTL(temp2(strlength:))
        READ(temp2,*,ERR=800,END=800) species
        CALL ATOMNUMBER(species,snumber)
        DO j=1,SIZE(atypes,1)
          IF( atypes(j,1)==0 .OR. atypes(j,1)==k ) THEN
            atypes(j,1) = k
            atypes(j,2) = NINT(snumber)
            EXIT
          ENDIF
        ENDDO
        strlength = SCAN(temp2," ")
        temp2 = ADJUSTL(temp2(strlength:))
      ENDDO
      !
      !
    ENDIF
    !
  ENDIF
ENDDO
!
!
!
300 CONTINUE
!If atom species were given in the input file, set them
IF( atypes(1,1).NE.0 ) THEN
  DO i=1,SIZE(P,1)
    DO j=1,SIZE(atypes,1)
      IF( NINT(AUX(i,typecol)) == atypes(j,1) ) THEN
        P(i,4) = DBLE(atypes(j,2))
      ENDIF
    ENDDO
  ENDDO
ENDIF
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
!
!
!
1000 CONTINUE
CLOSE(30)
!
!
END SUBROUTINE READ_XMD
!
!
END MODULE in_xmd
