MODULE swap
!
!**********************************************************************************
!*  SWAP                                                                          *
!**********************************************************************************
!* This module reads cartesian coordinates from an array and swaps atoms          *
!* of given indices, or swaps the two given Cartesian axes, or swaps two          *
!* given auxiliary properties.                                                    *
!**********************************************************************************
!* (C) August 2015 - Pierre Hirel                                                 *
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
USE messages
USE files
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE SWAP_XYZ(H,P,S,AUXNAMES,AUX,swap_id,SELECT)
!
!
IMPLICIT NONE
!
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: AUXNAMES !names of auxiliary properties
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT  !mask for atom list
INTEGER:: atype, type1, type2
INTEGER:: i
INTEGER:: Nswap  !number of atoms that were swapped
INTEGER,DIMENSION(2):: id   !indices of axis or atoms to swap
CHARACTER(LEN=16),DIMENSION(2),INTENT(IN):: swap_id  !Cartesian axes or indices of atoms to swap
REAL(dp),DIMENSION(2):: swap_sp  !atomic numbers to swap
REAL(dp),DIMENSION(4):: Vtemp
REAL(dp),DIMENSION(:),ALLOCATABLE:: AUXtemp  !auxiliary properties (temporary)
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of cores, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX  !auxiliary properties
!
!
!Initialize variables
atype = 0
type1=0
type2=0
i = 0
Nswap = -1
id(:)=0
!
!
msg = 'Entering SWAP_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Try to read integer numbers from swap_id
!If it succeeds then it means that two atoms must be swapped
!If it fails it does not matter
READ(swap_id(1),*,ERR=100,END=100) id(1)
READ(swap_id(2),*,ERR=100,END=100) id(2)
!
!
100 CONTINUE
SELECT CASE(swap_id(1))
!
CASE('X','x','Y','y','Z','z')
  !Two Cartesian axes must be swapped
  CALL ATOMSK_MSG( 2125, (/swap_id(1),swap_id(2)/),(/0.d0/) )
  !
  IF( swap_id(1) == swap_id(2) ) THEN
    !There is nothing to exchange => skip
    nwarn=nwarn+1
    CALL ATOMSK_MSG( 2757, (/""/),(/0.d0/) )
    GOTO 1000
  ENDIF
  !
  IF( swap_id(1)=="X" .OR. swap_id(1)=="x" ) THEN
    id(1) = 1
  ELSEIF( swap_id(1)=="Y" .OR. swap_id(1)=="y" ) THEN
    id(1) = 2
  ELSE
    id(1) = 3
  ENDIF
  !
  SELECT CASE(swap_id(2))
  CASE('X','x','Y','y','Z','z')
    IF( swap_id(2)=="X" .OR. swap_id(2)=="x" ) THEN
      id(2) = 1
    ELSEIF( swap_id(2)=="Y" .OR. swap_id(2)=="y" ) THEN
      id(2) = 2
    ELSE
      id(2) = 3
    ENDIF
    !
    !Swap the two axes
    Vtemp(1:3) = H(id(1),1:3)
    H(id(1),id(1)) = H(id(2),id(2))
    H(id(1),id(2)) = H(id(2),id(1))
    H(id(2),id(1)) = Vtemp(id(2))
    H(id(2),id(2)) = Vtemp(id(1))
    !
    IF( ALLOCATED(P) .AND. SIZE(P,1)>0 ) THEN
      !Swap coordinates of each atom
      DO i=1,SIZE(P,1)
        Vtemp(:) = P(i,:)
        P(i,id(1)) = Vtemp(id(2))
        P(i,id(2)) = Vtemp(id(1))
        !Swap shell positions (if any)
        IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
          Vtemp(:) = S(i,:)
          S(i,id(1)) = Vtemp(id(2))
          S(i,id(2)) = Vtemp(id(1))
        ENDIF
      ENDDO
    ENDIF
    Nswap=2
    !
  CASE DEFAULT
    !swap_id2 is not a Cartesian axis => big problem!
    !(this should have been dealt with before and not happen here,
    !however if for some reason it does happen let's ensure a smooth escape)
    nerr=nerr+1
  END SELECT
  !
  !
CASE DEFAULT
  IF( id(1)>0 .AND. id(2)>0 ) THEN
    !Two atoms must be swapped
    CALL ATOMSK_MSG( 2125, (/swap_id(1),swap_id(2)/),(/1.d0/) )
    !
    IF( .NOT.ALLOCATED(P) .OR. SIZE(P,1)<=0 ) THEN
      !No atom in system: can not apply option
      GOTO 1000
    ENDIF
    !
    IF( swap_id(1) == swap_id(2) ) THEN
      !There is nothing to exchange => skip
      nwarn=nwarn+1
      CALL ATOMSK_MSG( 2757, (/""/),(/0.d0/) )
      GOTO 1000
    ENDIF
    !Read indices of atoms to swap
    !
    !If indices of atoms are out of bounds, abort
    DO i=1,2
      IF( id(i)<=0 .OR. id(i)>SIZE(P,1) ) THEN
        nwarn=nwarn+1
        CALL ATOMSK_MSG(2742,(/''/),(/DBLE(id(i))/))
        GOTO 1000
      ENDIF
    ENDDO
    !
    !Save position of first atom
    Vtemp(:) = P(id(1),:)
    !
    !Swap atoms
    P(id(1),:) = P(id(2),:)
    P(id(2),:) = Vtemp(:)
    !
    !Swap shell positions (if any)
    IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
      Vtemp(:) = S(id(1),:)
      S(id(1),:) = S(id(2),:)
      S(id(2),:) = Vtemp(:)
    ENDIF
    !
    !Swap auxiliary properties (if any)
    IF( ALLOCATED(AUX) .AND. SIZE(AUX,1)==SIZE(P,1) ) THEN
      ALLOCATE( AUXtemp( SIZE(AUX,2) ) )
      AUXtemp(:) = AUX(id(1),:)
      AUX(id(1),:) = AUX(id(2),:)
      AUX(id(2),:) = AUXtemp(:)
      DEALLOCATE(AUXtemp)
    ENDIF
    Nswap=2
    !
  ELSE
    !The swap_id(:) may contain atom species to swap, *or* the names of auxiliary properties
    Nswap = 0
    CALL ATOMNUMBER(swap_id(1),swap_sp(1))
    CALL ATOMNUMBER(swap_id(2),swap_sp(2))
!     IF( swap_sp(1)<0.1d0 ) THEN
!       CALL ATOMSK_MSG(801,(/swap_id(1)/),(/0.d0/))
!       nerr = nerr+1
!       GOTO 1000
!     ELSEIF( swap_sp(2)<0.1d0 ) THEN
!       CALL ATOMSK_MSG(801,(/swap_id(2)/),(/0.d0/))
!       nerr = nerr+1
!       GOTO 1000
!     ENDIF
    !
    !
    IF( swap_sp(1)>0.1d0 .AND. swap_sp(2)>0.1d0 ) THEN
      !Two atomic species must be swapped
      CALL ATOMSK_MSG( 2125, (/swap_id(1),swap_id(2)/),(/2.d0/) )
      !
      IF( .NOT.ALLOCATED(P) .OR. SIZE(P,1)<=0 ) THEN
        !No atom in system: can not apply option
        GOTO 1000
      ENDIF
      !
      IF( swap_id(1) == swap_id(2) ) THEN
        !There is nothing to exchange => skip
        nwarn=nwarn+1
        CALL ATOMSK_MSG( 2757, (/""/),(/0.d0/) )
        GOTO 1000
      ENDIF
      !If "type" is defined as auxiliary property, also exchange the two atom types
      IF( ALLOCATED(AUXNAMES) .AND. SIZE(AUXNAMES)>0 ) THEN
        !Determine if a column contains "type"
        atype = 0
        DO i=1,SIZE(AUXNAMES)
          IF( AUXNAMES(i)=="type" ) THEN
            atype = i
          ENDIF
        ENDDO
        !
        !If "type" was found, determine the type of the atoms to swap
        IF( atype>0 ) THEN
          type1 = 0
          type2 = 0
          i = 0
          DO WHILE( type1==0 .OR. type2==0 )
            i=i+1
            IF( i>SIZE(P,1) ) EXIT
            IF( NINT(P(i,4))==NINT(swap_sp(1)) ) THEN
              type1 = NINT(AUX(i,atype))
            ELSEIF( NINT(P(i,4))==NINT(swap_sp(2)) ) THEN
              type2 = NINT(AUX(i,atype))
            ENDIF
          ENDDO
        ENDIF
      ENDIF
      !Swap species of types 1 and 2
      DO i=1,SIZE(P,1)
        IF( IS_SELECTED(SELECT,i) ) THEN
          IF( NINT(P(i,4))==NINT(swap_sp(1)) ) THEN
            P(i,4) = swap_sp(2)
            Nswap = Nswap+1
            !Swap shell positions (if any)
            IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
              S(i,4) = swap_sp(2)
            ENDIF
            !Swap the "type" of atoms (if any)
            IF( atype>0 ) THEN
              AUX(i,atype) = DBLE(type2)
            ENDIF
          ELSEIF( NINT(P(i,4))==NINT(swap_sp(2)) ) THEN
            P(i,4) = swap_sp(1)
            Nswap = Nswap+1
            !Swap shell positions (if any)
            IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
              S(i,4) = swap_sp(1)
            ENDIF
            !Swap the "type" of atoms (if any)
            IF( atype>0 ) THEN
              AUX(i,atype) = DBLE(type1)
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      !
    ELSE
      !The swap_id did not contain atom species
      !Then it must be the names of auxiliary properties
      CALL ATOMSK_MSG( 2125, (/swap_id(1),swap_id(2)/),(/3.d0/) )
      IF( swap_id(1) == swap_id(2) ) THEN
        !There is nothing to exchange => skip
        nwarn=nwarn+1
        CALL ATOMSK_MSG( 2757, (/""/),(/0.d0/) )
        GOTO 1000
      ENDIF
      !Locate those aux.prop. in AUXNAMES
      type1 = 0
      type2 = 0
      DO i=1,SIZE(AUXNAMES)
        IF( AUXNAMES(i)==swap_id(1) ) THEN
          type1 = i
        ENDIF
        IF( AUXNAMES(i)==swap_id(2) ) THEN
          type2 = i
        ENDIF
      ENDDO
      !
      !Check that indices type1 and type2 are positive
      !NB: it was already tested at the beginning that swap_id(1) != swap_id(2),
      !    so at this point type1 and type2 *must* be different
      IF( type1>0 .AND. type2>0 ) THEN
        !Swap the two aux.prop
        !Save data from column #type1 in AUXtemp
        ALLOCATE(AUXtemp(SIZE(AUX,1)))
        DO i=1,SIZE(AUX,1)
          IF( IS_SELECTED(SELECT,i) ) THEN
            AUXtemp(i) = AUX(i,type1)
          ENDIF
        ENDDO
        !Copy data from colum #type2 into column #type1
        DO i=1,SIZE(AUX,1)
          IF( IS_SELECTED(SELECT,i) ) THEN
            AUX(i,type1) = AUX(i,type2)
          ENDIF
        ENDDO
        !Copy data from AUXtemp into column #type2
        DO i=1,SIZE(AUX,1)
          IF( IS_SELECTED(SELECT,i) ) THEN
            AUX(i,type2) = AUXtemp(i)
            Nswap = Nswap+1
          ENDIF
        ENDDO
        DEALLOCATE(AUXtemp)
        !
        !Swap the names of aux.prop. (only if no selection was defined)
        IF( IS_SELECTED(SELECT,i) ) THEN
          msg = TRIM(ADJUSTL(AUXNAMES(type1)))
          AUXNAMES(type1) = TRIM(ADJUSTL(AUXNAMES(type2)))
          AUXNAMES(type2) = TRIM(ADJUSTL(msg))
        ENDIF
        !
      ELSE
        !No such aux.prop. in AUXNAMES
        nwarn=nwarn+1
        IF( type1==0 ) THEN
          CALL ATOMSK_MSG(2730,(/swap_id(1)/),(/0.d0/))
        ENDIF
        IF( type2==0 ) THEN
          CALL ATOMSK_MSG(2730,(/swap_id(2)/),(/0.d0/))
        ENDIF
        GOTO 1000
      ENDIF
    ENDIF
    !
  ENDIF
  !
END SELECT
!
!
CALL ATOMSK_MSG(2126,(/''/),(/DBLE(Nswap)/))
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
END SUBROUTINE SWAP_XYZ
!
!
END MODULE swap
