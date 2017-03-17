MODULE shift
!
!**********************************************************************************
!*  SHIFT                                                                         *
!**********************************************************************************
!* This module reads cartesian coordinates from an array and                      *
!* shifts atoms that are above (or below) a given plane.                          *
!**********************************************************************************
!* (C) March 2010 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 04 Aug. 2015                                     *
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
SUBROUTINE SHIFT_XYZ(P,S,shift_dir,shift_dist,shift_axis,shift_tau1,shift_tau2,shift_tau3,ORIENT,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=5):: shift_dir    !'above' or 'below' (or empty)
CHARACTER(LEN=16):: shift_axis  !x, y, z, or crystallographic direction
CHARACTER(LEN=128):: msg
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT  !mask for atom list
INTEGER:: a1
INTEGER:: shiftdir       !-1=below; +1=above
INTEGER:: i, j, NPshifted
REAL(dp):: shift_dist       !distance of the "cut plane" to the origin
REAL(dp),INTENT(IN):: shift_tau1, shift_tau2, shift_tau3  !shift vector
REAL(dp):: tempreal
REAL(dp):: V1, V2, V3  !vector components
REAL(dp),DIMENSION(1,3):: Vplane  !crystallographic vector defining the plane
REAL(dp),DIMENSION(3,3),INTENT(IN):: ORIENT !current crystallographic orientation of the system
REAL(dp),DIMENSION(3,3):: ORIENTN      !normalized ORIENT
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of cores, shells
!
!
!Initialize variables
i = 0
shiftdir = 0
NPshifted = 0    !Number of shifted atoms
!
!
msg = 'Entering SHIFT_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
IF(shift_dir=='above') THEN
  shiftdir = 1
ELSEIF(shift_dir=='below') THEN
  shiftdir = -1
ELSE
  IF( ALLOCATED(SELECT) ) THEN
    shift_dir = "selec"
  ELSE
    shift_dir = ""
  ENDIF
ENDIF
!
CALL ATOMSK_MSG( 2085, (/ shift_dir, shift_axis//'    ' /), &
                & (/shift_dist, shift_tau1, shift_tau2, shift_tau3/) )
!
WRITE(msg,*) 'a1: ', a1
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!If shift is zero then skip
IF(shift_tau1==0.d0 .AND. shift_tau2==0.d0 .AND. shift_tau3==0.d0) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2736,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
200 CONTINUE
IF( shift_dir=='above' .OR. shift_dir=='below' ) THEN
  !Shift atoms above or below a plane
  SELECT CASE(shift_axis)
  CASE("x","X","y","Y","z","Z")
    !Shift atoms above or below the shift_dist along a cartesian axis
    !Define the axes
    IF(shift_axis=='x' .OR. shift_axis=='X') THEN
      a1 = 1
    ELSEIF(shift_axis=='y' .OR. shift_axis=='Y') THEN
      a1 = 2
    ELSEIF(shift_axis=='z' .OR. shift_axis=='Z') THEN
      a1 = 3
    ELSE
      nwarn=nwarn+1
      CALL ATOMSK_MSG(2800,(/TRIM(shift_dir)/),(/0.d0/))
      GOTO 1000
    ENDIF
    !
    !Shift atoms (or ion cores)
    DO i=1, SIZE(P(:,1))
      IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
        IF( shiftdir>0 .AND. P(i,a1)>shift_dist  .OR.                       &
          & shiftdir<0 .AND. P(i,a1)<shift_dist          ) THEN
          NPshifted = NPshifted+1
          P(i,1) = P(i,1) + shift_tau1
          P(i,2) = P(i,2) + shift_tau2
          P(i,3) = P(i,3) + shift_tau3
          !
          !If shells exist, do the same
          IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
            S(i,1) = S(i,1) + shift_tau1
            S(i,2) = S(i,2) + shift_tau2
            S(i,3) = S(i,3) + shift_tau3
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    !
    !
  CASE DEFAULT
    !cutdir should contain a crystallograhic direction
    !convert it to a vector and save it in Vplane(1,:)
    CALL INDEX_MILLER(shift_axis,Vplane(1,:),j)
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
    !Shift atoms (or ion cores)
    DO i=1, SIZE(P(:,1))
      IF(.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) THEN
        !determine if atom position is above or below the plane
        tempreal = VEC_PLANE( Vplane(1,:) , shift_dist , P(i,1:3) )
        !
        IF( shiftdir>0 .AND. tempreal>0.d0  .OR.                       &
          & shiftdir<0 .AND. tempreal<0.d0           ) THEN
          NPshifted = NPshifted+1
          P(i,1) = P(i,1) + shift_tau1
          P(i,2) = P(i,2) + shift_tau2
          P(i,3) = P(i,3) + shift_tau3
          !
          !If shells exist, do the same
          IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
            S(i,1) = S(i,1) + shift_tau1
            S(i,2) = S(i,2) + shift_tau2
            S(i,3) = S(i,3) + shift_tau3
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    !
    !
  END SELECT
  !
  !
ELSE
  !Shift all atoms (or selected atoms)
  DO i=1, SIZE(P,1)
    IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
      NPshifted = NPshifted+1
      P(i,1) = P(i,1) + shift_tau1
      P(i,2) = P(i,2) + shift_tau2
      P(i,3) = P(i,3) + shift_tau3
      !If shells exist, do the same
      IF( ALLOCATED(S) .AND. SIZE(S,1).NE.0 ) THEN
        S(i,1) = S(i,1) + shift_tau1
        S(i,2) = S(i,2) + shift_tau2
        S(i,3) = S(i,3) + shift_tau3
      ENDIF
    ENDIF
  ENDDO
ENDIF
!
!
!
300 CONTINUE
CALL ATOMSK_MSG(2086,(/''/),(/DBLE(NPshifted)/))
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
END SUBROUTINE SHIFT_XYZ
!
!
END MODULE shift
