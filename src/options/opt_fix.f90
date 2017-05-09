MODULE fix
!
!**********************************************************************************
!*  FIX                                                                           *
!**********************************************************************************
!* This module allows to fix some atoms along a direction, and save this          *
!* state in auxiliary properties with names fixx, fixy and fixz. A value          *
!* of 1.d0 means that the atom is fixed along X, Y and Z respectively; a          *
!* value of 0.d0 means that the atom is not fixed.                                *
!**********************************************************************************
!* (C) July 2011 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 09 May 2017                                      *
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
USE messages
USE files
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE FIX_XYZ(P,AUXNAMES,AUX,fixaxis,fix_dir,fixdistance,fixdir,ORIENT,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=5):: fix_dir  !above or below
CHARACTER(LEN=5),INTENT(IN):: fixaxis  !direction along which atoms are fixed (x,y,z,all)
CHARACTER(LEN=16),INTENT(IN):: fixdir   !x, y, z, or crystallographic direction
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: newAUXNAMES !names of auxiliary properties
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT  !mask for atom list
INTEGER:: a1
INTEGER:: fixx, fixy, fixz !columns of AUX in which fixes are stored
INTEGER:: i, j
INTEGER:: NPfixed !number of atoms fixed
REAL(dp),INTENT(IN):: fixdistance
REAL(dp):: tempreal
REAL(dp):: V1, V2, V3  !vector components
REAL(dp),DIMENSION(1,3):: Vplane  !crystallographic vector defining the plane
REAL(dp),DIMENSION(3,3),INTENT(IN):: ORIENT !current crystallographic orientation of the system
REAL(dp),DIMENSION(3,3):: ORIENTN      !normalized ORIENT
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P !atom positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX !auxiliary properties
!
!Initialize variables
fixx=0
fixy=0
fixz=0
NPfixed=0
IF(ALLOCATED(newAUXNAMES)) DEALLOCATE(newAUXNAMES)
IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
!
WRITE(msg,*) 'Entering FIX_XYZ: ', fixaxis, fix_dir, fixdir
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
IF( fix_dir.NE.'above' .AND. fix_dir.NE.'below' ) THEN
  IF( ALLOCATED(SELECT) ) THEN
    fix_dir='selec'
  ELSE
    fix_dir=''
  ENDIF
ENDIF
CALL ATOMSK_MSG(2097,(/fixaxis,fix_dir,fixdir//'    '/),(/fixdistance/))
!
!
!
100 CONTINUE
!Set up AUX array
IF( .NOT.ALLOCATED(AUX) ) THEN
  !No auxiliary property exist => create a new array
  ALLOCATE(AUXNAMES(3))
  ALLOCATE( AUX( SIZE(P(:,1)),3 ) )
  fixx=1
  fixy=2
  fixz=3
  AUX(:,:) = 0.d0  !Default: all atoms are NOT fixed
  !
ELSE
  !Some auxiliary properties already exist
  !
  !Check if fixx, fixy and fixz already exist in AUXNAMES
  !If they do exist, use their indices
  DO i=1,SIZE(AUXNAMES)
    IF(AUXNAMES(i)=="fixx") fixx = i
    IF(AUXNAMES(i)=="fixy") fixy = i
    IF(AUXNAMES(i)=="fixz") fixz = i
  ENDDO
  !
  !If fixx, fixy and fixz do not already exist, create them
  IF( fixx==0 .OR. fixy==0 .OR. fixz==0 ) THEN
    ALLOCATE( newAUXNAMES( SIZE(AUXNAMES)+3 ) )
    ALLOCATE( newAUX( SIZE(P(:,1)), SIZE(AUXNAMES)+3 ) )
    newAUX(:,:) = 0.d0
    fixx=SIZE(AUXNAMES)+1
    fixy=fixx+1
    fixz=fixy+1
    !
    !Set the names
    DO i=1,SIZE(AUXNAMES)
      newAUXNAMES(i) = AUXNAMES(i)
    ENDDO
    !
    DO i=1,SIZE(P(:,1))
      newAUX(i,1:fixx-1) = AUX(i,:)
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
AUXNAMES(fixx) = "fixx"
AUXNAMES(fixy) = "fixy"
AUXNAMES(fixz) = "fixz"
!
!
!
200 CONTINUE
WRITE(msg,*) 'fix_dir = ', fix_dir
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
IF( fix_dir=='above' .OR. fix_dir=='below' ) THEN
  WRITE(msg,*) 'fixdir = ', fixdir
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  SELECT CASE(fixdir)
  CASE("x","X","y","Y","z","Z")
    !Define the axes
    IF(fixdir=='x' .OR. fixdir=='X') THEN
      a1 = 1
    ELSEIF(fixdir=='y' .OR. fixdir=='Y') THEN
      a1 = 2
    ELSEIF(fixdir=='z' .OR. fixdir=='Z') THEN
      a1 = 3
    ELSE
      nerr = nerr+1
      CALL ATOMSK_MSG(2800,(/fixdir/),(/0.d0/))
      GOTO 1000
    ENDIF
    !
    !
    !Fix atoms that are above/below <fixdistance> along <fixdir>
    ! *or* atoms that are selected inside the region
    DO i=1,SIZE(P,1)
      IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
        IF( (fix_dir=='above' .AND. P(i,a1)>fixdistance) .OR.        &
          & (fix_dir=='below' .AND. P(i,a1)<fixdistance)      ) THEN
          !If atom is in the region of interest, fix it
          SELECT CASE(fixaxis)
          CASE("x","X")
            AUX(i,fixx) = 1.d0
          CASE("y","Y")
            AUX(i,fixy) = 1.d0
          CASE("z","Z")
            AUX(i,fixz) = 1.d0
          CASE("all")
            AUX(i,fixx) = 1.d0
            AUX(i,fixy) = 1.d0
            AUX(i,fixz) = 1.d0
          END SELECT
          NPfixed = NPfixed+1
        ENDIF
      ENDIF
    ENDDO
    !
  CASE DEFAULT
    !cutdir should contain a crystallograhic direction
    !convert it to a vector and save it in Vplane(1,:)
    CALL INDEX_MILLER(fixdir,Vplane(1,:),j)
    IF(j>0) GOTO 800
    !
    !If the system has a defined crystallographic orientation ORIENT,
    !then Vplane(1,:) is defined in that basis
    !=> rotate Vplane(1,:) to express it in cartesian basis
    IF( ANY( NINT(ORIENT(:,:)).NE.0 ) ) THEN
      DO i=1,3
        ORIENTN(i,:) = ORIENT(i,:) / VECLENGTH(ORIENT(i,:))
      ENDDO
      V1 = Vplane(1,1)
      V2 = Vplane(1,2)
      V3 = Vplane(1,3)
      Vplane(1,1) = ORIENTN(1,1)*V1 + ORIENTN(1,2)*V2 + ORIENTN(1,3)*V3
      Vplane(1,2) = ORIENTN(2,1)*V1 + ORIENTN(2,2)*V2 + ORIENTN(2,3)*V3
      Vplane(1,3) = ORIENTN(3,1)*V1 + ORIENTN(3,2)*V2 + ORIENTN(3,3)*V3
    ENDIF
    !Normalize Vplane
    Vplane(1,:) = Vplane(1,:)/VECLENGTH(Vplane(1,:))
    !
    DO i=1,SIZE(P,1)
      !determine if atom is above or below the plane
      tempreal = VEC_PLANE( Vplane(1,:) , fixdistance , P(i,1:3) )
      IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
        IF( (fix_dir=='above' .AND. tempreal>0.d0) .OR.        &
          & (fix_dir=='below' .AND. tempreal<0.d0)       ) THEN
          !If atom is in the region of interest, fix it
          SELECT CASE(fixaxis)
          CASE("x","X")
            AUX(i,fixx) = 1.d0
          CASE("y","Y")
            AUX(i,fixy) = 1.d0
          CASE("z","Z")
            AUX(i,fixz) = 1.d0
          CASE("all")
            AUX(i,fixx) = 1.d0
            AUX(i,fixy) = 1.d0
            AUX(i,fixz) = 1.d0
          END SELECT
          NPfixed = NPfixed+1
        ENDIF
      ENDIF
    ENDDO
    !
  END SELECT
  !
  !
ELSE
  !Fix all atoms, or selected atom
  DO i=1,SIZE(P,1)
    IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
      !If atom is in the region of interest, fix it
      SELECT CASE(fixaxis)
      CASE("x","X")
        AUX(i,fixx) = 1.d0
      CASE("y","Y")
        AUX(i,fixy) = 1.d0
      CASE("z","Z")
        AUX(i,fixz) = 1.d0
      CASE("all")
        AUX(i,fixx) = 1.d0
        AUX(i,fixy) = 1.d0
        AUX(i,fixz) = 1.d0
      END SELECT
      NPfixed = NPfixed+1
    ENDIF
  ENDDO
ENDIF
!
!
!
400 CONTINUE
CALL ATOMSK_MSG(2098,(/""/),(/DBLE(NPfixed)/))
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
END SUBROUTINE FIX_XYZ
!
!
!
END MODULE fix