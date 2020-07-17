MODULE cell
!
!**********************************************************************************
!*  CELL                                                                          *
!**********************************************************************************
!* This module modifies the cell vectors (but does NOT modify atom positions).    *
!* The cell may be elongated ("add") or shortened ("rm") along one or all         *
!* directions, or given a tilt. This module can also try to provide a             *
!* rectangular bounding box, however this is not very accurate and may generate   *
!* a box unsuitable for some systems.                                             *
!**********************************************************************************
!* (C) July 2020 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 16 July 2020                                     *
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
USE deterH
USE messages
USE files
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE CELL_XYZ(H,P,cellop,celllength,celldir)
!
!
IMPLICIT NONE
CHARACTER(LEN=5),INTENT(IN):: celldir  !direction in which the cell is modified: x,y,z,xy,yx,xz,zx,yz,zy,xyz,all
CHARACTER(LEN=5),INTENT(IN):: cellop   !operation to perform on the cell: add,rm,set,auto (or "rebox")
CHARACTER(LEN=128):: msg
INTEGER:: a1, a2, a3
INTEGER:: i
INTEGER:: Nmodified  !number of box vectors that were modified
INTEGER,DIMENSION(3):: Hmodified !index of modified vectors
REAL(dp):: alpha !angle between a cell vector and a Cartesian axis
REAL(dp):: dH  !change of cell vector (+celllength for "add" or -celllength for "rm")
REAL(dp):: vl  !length of a cell vector
REAL(dp),DIMENSION(3):: cartvec  !coordinates of a Cartesian vector
REAL(dp),INTENT(IN):: celllength  !value by which the cell will be modified
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P  !positions of atoms (not modified by this module)
!
!
!Initialize variables
i = 0
a1 = 0
a2 = 0
a3 = 0
Nmodified = 0
Hmodified(:) = 0
dH = 0.d0
 cartvec(:) = 0.d0
!
!
100 CONTINUE
msg = 'Entering CELL_XYZ...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(2151,(/cellop,celldir/),(/celllength/))
!
!Verify that celllength is not zero
IF( celllength==0.d0 ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2765,(/""/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Set direction(s) to modify
SELECT CASE(celldir)
CASE("H1",'x','X')
  Hmodified(1) = 1
  Hmodified(2) = 0
  Hmodified(3) = 0
  cartvec(1) = 1.d0
CASE("H2",'y','Y')
  Hmodified(1) = 0
  Hmodified(2) = 1
  Hmodified(3) = 0
  cartvec(2) = 1.d0
CASE("H3",'z','Z')
  Hmodified(1) = 0
  Hmodified(2) = 0
  Hmodified(3) = 1
  cartvec(3) = 1.d0
CASE("xyz","XYZ","all","ALL")
  Hmodified(1) = 1
  Hmodified(2) = 1
  Hmodified(3) = 1
CASE("xy","XY")
  a1 = 2
  a2 = 1
CASE("xz","XZ")
  a1 = 3
  a2 = 1
CASE("yx","YX")
  a1 = 1
  a2 = 2
CASE("yz","YZ")
  a1 = 3
  a2 = 2
CASE("zx","ZX")
  a1 = 1
  a2 = 3
CASE("zy","ZY")
  a1 = 2
  a2 = 3
CASE DEFAULT
  Hmodified(:) = 0
END SELECT
!
!Set dH according to operation to perform (add or rm)
IF( cellop=="add" .OR. cellop=="ADD" ) THEN
  dH = celllength
ELSEIF( cellop=="rm" .OR. cellop=="RM" ) THEN
  dH = -1.d0 * celllength
ENDIF
!
!
!
200 CONTINUE
!Modify the box
IF( cellop=="add" .OR. cellop=="rm" .OR. cellop=="set" ) THEN
  SELECT CASE(celldir)
  CASE("H1","H2","H3")
    !Modify vector length along given directions
    DO i=1,3
      IF( Hmodified(i)==1 ) THEN
        !Compute length of cell vector
        vl = VECLENGTH(H(i,:))
        IF( cellop=="set" ) THEN
          !Resize vector to specified length
          H(i,:) = celllength*H(i,:) / vl
        ELSE
          !Add or remove given length to cell vector
          H(i,:) = (vl + dH)*H(i,:) / vl
        ENDIF
        Nmodified = Nmodified+1
      ENDIF
    ENDDO
    !
  CASE('x','X','y','Y','z','Z',"xyz","XYZ","all","ALL")
    !Expand box along Cartesian axis
    !Modify vector length along given directions
    DO i=1,3
      IF( Hmodified(i)==1 ) THEN
        !Determine which cell vector has longest component along Cartesian direction i
        IF( H(1,i)>H(2,i) .AND. H(1,i)>H(3,i) ) THEN
          !H1 will be modified
          a1 = 1
        ELSEIF( H(2,i)>H(1,i) .AND. H(2,i)>H(3,i) ) THEN
          !H2 will be modified
          a1 = 2
        ELSE
          !H3 will be modified
          a1 = 3
        ENDIF
        !Compute length of this cell vector
        vl = VECLENGTH(H(a1,:))
        !Compute angle between this cell vector and given Cartesian axis
        alpha = ANGVEC( H(a1,:) , cartvec )
        IF( cellop=="set" ) THEN
          !Resize vector to specified length
          H(a1,:) = celllength*H(a1,:) / (vl*DSIN(alpha))
        ELSE
          !Add or remove given length to cell vector
          H(a1,:) = ( vl + (dH/DCOS(alpha)) ) * H(a1,:) / vl
        ENDIF
        Nmodified = Nmodified+1
      ENDIF
    ENDDO
    !
  CASE("xy","XY","yx","YX","xz","XZ","zx","ZX","yz","YZ","zy","ZY")
    !Change box tilt
    IF( cellop=="set" ) THEN
      !Set this component to specified value
      H(a1,a2) = celllength
    ELSE
      !Add or remove given value to this component
      H(a1,a2) = H(a1,a2) + dH
    ENDIF
    Hmodified(a1) = 1
    !
  END SELECT
  !
ENDIF
!
!
!
300 CONTINUE
CALL ATOMSK_MSG(2152,(/''/),(/0.d0/))
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
END SUBROUTINE CELL_XYZ
!
!
END MODULE cell
