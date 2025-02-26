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
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 26 Feb. 2025                                     *
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
SUBROUTINE SHIFT_XYZ(H,P,S,shift_dir,shift_dist,shift_axis,shift_tau1,shift_tau2,shift_tau3,ORIENT,SELECT)
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
REAL(dp),DIMENSION(3,3),INTENT(IN):: H      !box vectors
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
IF( .NOT.ALLOCATED(P) .OR. SIZE(P,1)<=0 ) THEN
  !No atom in system: can not apply option
  GOTO 1000
ENDIF
!
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
    DO i=1,SIZE(P,1)
      IF(IS_SELECTED(SELECT,i)) THEN
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
    !shift_axis should contain a crystallograhic direction
    !Convert "shift_axis" into a Cartesian vector and save it in Vplane(1,:)
    CALL MILLER2VEC(H,shift_axis,ORIENT,Vplane(1,:),j)
    !
    !Check return status j (0=success, otherwise there was an error)
    IF( j>0 ) THEN
      IF( j==2 ) THEN
        !The error was because i is not equal to -h-k
        nerr=nerr+1
        CALL ATOMSK_MSG(815,(/shift_axis/),(/0.d0/))
        GOTO 1000
      ELSE
        !Other error, unable to convert this string into a proper vector
        CALL ATOMSK_MSG(817,(/TRIM(shift_axis)/),(/0.d0/))
        GOTO 1000
      ENDIF
    ENDIF
    !
    !Shift atoms (or ion cores)
    DO i=1, SIZE(P,1)
      IF(IS_SELECTED(SELECT,i)) THEN
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
  DO i=1,SIZE(P,1)
    IF( IS_SELECTED(SELECT,i) ) THEN
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
