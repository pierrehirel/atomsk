MODULE freeze
!
!**********************************************************************************
!*  FREEZE                                                                        *
!**********************************************************************************
!* This module allows to freeze some atoms along a direction, and save this       *
!* state as an auxiliary property that can be used in some output formats,        *
!* typically to prevent some atoms to relax or exclude them from dynamics.        *
!* A value of 1.d0 means that the atom is fixed along X, Y and Z respectively;    *
!* a value of 0.d0 means that the atom is not fixed.                              *
!**********************************************************************************
!* (C) July 2011 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 14 Jan. 2025                                     *
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
USE comv
USE constants
USE crystallography
USE messages
USE files
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE FREEZE_XYZ(H,P,AUXNAMES,AUX,freeze_coord,freeze_side,freeze_dist,freeze_normal,ORIENT,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=5):: freeze_side  !above or below
CHARACTER(LEN=5),INTENT(IN)::  freeze_coord   !coordinate to freeze (x,y,z,xy,xz,yz,all)
CHARACTER(LEN=16),INTENT(IN):: freeze_normal  !x, y, z, or crystallographic direction
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: newAUXNAMES !names of auxiliary properties
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT  !mask for atom list
INTEGER:: a1
INTEGER:: freezex, freezey, freezez !columns of AUX in which fixes are stored
INTEGER:: i, j
INTEGER:: NPfrozen !number of atoms frozen
REAL(dp),INTENT(IN):: freeze_dist
REAL(dp):: tempreal
REAL(dp):: V1, V2, V3  !vector components
REAL(dp),DIMENSION(1,3):: Vplane  !crystallographic vector defining the plane
REAL(dp),DIMENSION(3,3),INTENT(IN):: H      !box vectors
REAL(dp),DIMENSION(3,3),INTENT(IN):: ORIENT !current crystallographic orientation of the system
REAL(dp),DIMENSION(3,3):: ORIENTN      !normalized ORIENT
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P !atom positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX !auxiliary properties
!
!Initialize variables
freezex=0
freezey=0
freezez=0
NPfrozen=0
IF(ALLOCATED(newAUXNAMES)) DEALLOCATE(newAUXNAMES)
IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
!
WRITE(msg,*) "Entering FREEZE_XYZ: ", freeze_coord, freeze_side, " ", freeze_normal, freeze_dist
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
IF( freeze_side.NE."above" .AND. freeze_side.NE."below" ) THEN
  IF( ALLOCATED(SELECT) ) THEN
    freeze_side='selec'
  ELSE
    freeze_side=''
  ENDIF
ENDIF
CALL ATOMSK_MSG(2097,(/freeze_coord//"     ",freeze_side//"     ",freeze_normal/),(/freeze_dist/))
!
!
!
100 CONTINUE
!Set up AUX array
IF( .NOT.ALLOCATED(AUX) ) THEN
  !No auxiliary property exist => create a new array
  ALLOCATE(AUXNAMES(3))
  ALLOCATE( AUX( SIZE(P(:,1)),3 ) )
  freezex=1
  freezey=2
  freezez=3
  AUX(:,:) = 0.d0  !Default: all atoms are NOT frozen (i.e. they are "mobile")
  !
ELSE
  !Some auxiliary properties already exist
  !
  !Check if properties named  already exist in AUXNAMES
  !If they do exist, use their indices
  DO i=1,SIZE(AUXNAMES)
    IF( AUXNAMES(i)(1:3)=="fix" .OR. AUXNAMES(i)(1:6)=="freeze" ) THEN
      j = LEN_TRIM(AUXNAMES(i))+1
      IF( SCAN(AUXNAMES(i)(j:),'x')>0 .OR. SCAN(AUXNAMES(i)(j:),'X')>0 ) freezex = i
      IF( SCAN(AUXNAMES(i)(j:),'y')>0 .OR. SCAN(AUXNAMES(i)(j:),'Y')>0 ) freezex = i
      IF( SCAN(AUXNAMES(i)(j:),'z')>0 .OR. SCAN(AUXNAMES(i)(j:),'Z')>0 ) freezex = i
    ENDIF
  ENDDO
  !
  !If "freeze" coordinates do not already exist, create them
  !NOTE: if one or more is missing we ignore previous proeprties and create new ones
  IF( freezex==0 .OR. freezey==0 .OR. freezez==0 ) THEN
    ALLOCATE( newAUXNAMES( SIZE(AUXNAMES)+3 ) )
    ALLOCATE( newAUX( SIZE(P(:,1)), SIZE(AUXNAMES)+3 ) )
    newAUX(:,:) = 0.d0
    freezex=SIZE(AUXNAMES)+1
    freezey=freezex+1
    freezez=freezey+1
    !
    !Set the names
    DO i=1,SIZE(AUXNAMES)
      newAUXNAMES(i) = AUXNAMES(i)
    ENDDO
    !
    DO i=1,SIZE(P(:,1))
      newAUX(i,1:freezex-1) = AUX(i,:)
    ENDDO
    !
    !Replace old AUXNAMES by newAUXNAMES
    DEALLOCATE(AUXNAMES)
    ALLOCATE( AUXNAMES( SIZE(newAUXNAMES) ) )
    AUXNAMES = newAUXNAMES
    DEALLOCATE(newAUXNAMES)
    !Same with AUX
    DEALLOCATE(AUX)
    ALLOCATE( AUX( SIZE(newAUX(:,1)), SIZE(newAUX(1,:)) ) )
    AUX = newAUX
    DEALLOCATE(newAUX)
  ENDIF
ENDIF
!
!Set names of auxiliary properties
AUXNAMES(freezex) = "freeze_x"
AUXNAMES(freezey) = "freeze_y"
AUXNAMES(freezez) = "freeze_z"
!
!
!
200 CONTINUE
WRITE(msg,*) "freeze_side = ", freeze_side
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
IF( freeze_side=="above" .OR. freeze_side=="below" ) THEN
  WRITE(msg,*) "FREEZE: ", freeze_normal, " ", freeze_coord
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  SELECT CASE(StrDnCase(freeze_normal))
  CASE('x','y','z')
    !Define the normal to the plane
    IF( StrDnCase(freeze_normal)=='x' ) THEN
      a1 = 1
    ELSEIF( StrDnCase(freeze_normal)=='y' ) THEN
      a1 = 2
    ELSEIF( StrDnCase(freeze_normal)=='z' ) THEN
      a1 = 3
    ELSE
      nerr = nerr+1
      CALL ATOMSK_MSG(2800,(/freeze_normal/),(/0.d0/))
      GOTO 1000
    ENDIF
    !
    !Freeze atoms that are above/below <freeze_dist> along <freeze_normal>
    ! *or* atoms that are selected inside the region
    DO i=1,SIZE(P,1)
      IF(IS_SELECTED(SELECT,i)) THEN
        !Atom #i is selected: determine if it is on the correct side of the plane
        IF( (freeze_side=='above'  .AND.  P(i,a1)>freeze_dist) .OR.        &
          & (freeze_side=='below'  .AND.  P(i,a1)<freeze_dist)      ) THEN
          !Atom #i is on the frozen side: freeze it
          SELECT CASE(StrDnCase(freeze_coord))
          CASE('x')
            AUX(i,freezex) = 1.d0
          CASE('y')
            AUX(i,freezey) = 1.d0
          CASE('z')
            AUX(i,freezez) = 1.d0
          CASE("xy","yx")
            AUX(i,freezex) = 1.d0
            AUX(i,freezey) = 1.d0
          CASE("xz","zx")
            AUX(i,freezex) = 1.d0
            AUX(i,freezez) = 1.d0
          CASE("yz","zy")
            AUX(i,freezey) = 1.d0
            AUX(i,freezez) = 1.d0
          CASE("all","xyz","xzy","yxz","yzx","zxy","zyx")
            AUX(i,freezex) = 1.d0
            AUX(i,freezey) = 1.d0
            AUX(i,freezez) = 1.d0
          END SELECT
          NPfrozen = NPfrozen+1
        ENDIF
      ENDIF
    ENDDO
    !
  CASE DEFAULT
    !freeze_normal should contain a crystallograhic direction
    WRITE(msg,*) 'Looking for a crystal direction... '
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !
    !Convert "freeze_normal" into a Cartesian vector and save it in Vplane(1,:)
    CALL MILLER2VEC(H,freeze_normal,ORIENT,Vplane(1,:),j)
    !
    !Check return status j (0=success, otherwise there was an error)
    IF( j>0 ) THEN
      IF( j==2 ) THEN
        !The error was because i is not equal to -h-k
        nerr=nerr+1
        CALL ATOMSK_MSG(815,(/freeze_normal/),(/0.d0/))
        GOTO 1000
      ELSE
        !Other error, unable to convert this string into a proper vector
        CALL ATOMSK_MSG(817,(/TRIM(freeze_normal)/),(/0.d0/))
        GOTO 1000
      ENDIF
    ENDIF
    !
    DO i=1,SIZE(P,1)
      IF(IS_SELECTED(SELECT,i)) THEN
        !Atom #i is selected: determine if it is on the correct side of the plane
        tempreal = VEC_PLANE( Vplane(1,:) , freeze_dist , P(i,1:3) )
        IF( (freeze_side=='above' .AND. tempreal>0.d0) .OR.        &
          & (freeze_side=='below' .AND. tempreal<0.d0)       ) THEN
          !Atom #i is on the frozen side: freeze it
          SELECT CASE(StrDnCase(freeze_coord))
          CASE('x')
            AUX(i,freezex) = 1.d0
          CASE('y')
            AUX(i,freezey) = 1.d0
          CASE('z')
            AUX(i,freezez) = 1.d0
          CASE("xy","yx")
            AUX(i,freezex) = 1.d0
            AUX(i,freezey) = 1.d0
          CASE("xz","zx")
            AUX(i,freezex) = 1.d0
            AUX(i,freezez) = 1.d0
          CASE("yz","zy")
            AUX(i,freezey) = 1.d0
            AUX(i,freezez) = 1.d0
          CASE("all","xyz","xzy","yxz","yzx","zxy","zyx")
            AUX(i,freezex) = 1.d0
            AUX(i,freezey) = 1.d0
            AUX(i,freezez) = 1.d0
          END SELECT
          NPfrozen = NPfrozen+1
        ENDIF
      ENDIF
    ENDDO
    !
  END SELECT
  !
  !
ELSE
  !Freeze all atoms, or selected atoms
  WRITE(msg,*) "FREEZE: all atoms"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(P,1)
    IF(IS_SELECTED(SELECT,i)) THEN
      !Atom #i is selected: freeze it
      SELECT CASE(StrDnCase(freeze_coord))
      CASE('x')
        AUX(i,freezex) = 1.d0
      CASE('y')
        AUX(i,freezey) = 1.d0
      CASE('z')
        AUX(i,freezez) = 1.d0
      CASE("xy","yx")
        AUX(i,freezex) = 1.d0
        AUX(i,freezey) = 1.d0
      CASE("xz","zx")
        AUX(i,freezex) = 1.d0
        AUX(i,freezez) = 1.d0
      CASE("yz","zy")
        AUX(i,freezey) = 1.d0
        AUX(i,freezez) = 1.d0
      CASE("all","xyz","xzy","yxz","yzx","zxy","zyx")
        AUX(i,freezex) = 1.d0
        AUX(i,freezey) = 1.d0
        AUX(i,freezez) = 1.d0
      END SELECT
      NPfrozen = NPfrozen+1
    ENDIF
  ENDDO
ENDIF
!
!
!
400 CONTINUE
CALL ATOMSK_MSG(2098,(/""/),(/DBLE(NPfrozen)/))
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
!
!
END SUBROUTINE FREEZE_XYZ
!
!
!
END MODULE freeze
