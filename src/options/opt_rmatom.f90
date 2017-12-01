MODULE rmatom
!
!**********************************************************************************
!*  RMATOM                                                                        *
!**********************************************************************************
!* This module removes either one atom given its index, or all atoms of           *
!* a given species, or all selected atoms.                                        *
!**********************************************************************************
!* (C) Dec 2011 - Pierre Hirel                                                    *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 30 Nov. 2017                                     *
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
USE atoms
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
SUBROUTINE RMATOM_XYZ(P,S,AUX,rmatom_prop,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=2):: species
CHARACTER(LEN=32),INTENT(IN):: rmatom_prop
CHARACTER(LEN=128):: msg
LOGICAL:: clear_sel  !clear selection at the end of option?
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: SELECT   !mask for atom list
LOGICAL,DIMENSION(:),ALLOCATABLE:: newSELECT              !mask for atom list (temporary)
INTEGER:: atomindex !index of atom that must be removed
INTEGER:: method    !1=selected atoms; 2=atom species; 3=atom index; 0=error
INTEGER:: NP !number of particles
INTEGER:: i
INTEGER:: rmatoms !number of atoms removed
REAL(dp):: snumber !atomic number
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S  !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newP, newS          !positions of atoms, shells (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX   !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX              !auxiliary properties (temporary)
!
species = ''
 clear_sel = .FALSE.
atomindex = 0
i = 0
method=0
NP = 0
rmatoms = 0
snumber = 0.d0
IF(ALLOCATED(newP)) DEALLOCATE(newP)
IF(ALLOCATED(newS)) DEALLOCATE(newS)
IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
IF(ALLOCATED(newSELECT)) DEALLOCATE(newSELECT)
!
WRITE(msg,*) 'Entering RMATOM_XYZ: '//TRIM(ADJUSTL(rmatom_prop))
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!
100 CONTINUE
!Read "rmatom_prop" to know what to remove
IF( rmatom_prop(1:3)=="sel" .OR. rmatom_prop(1:3)=="SEL" ) THEN
  method=1
  CALL ATOMSK_MSG(2102,(/"SEL"/),(/0.d0/))
ELSE
  !Check if rmatom_prop contains an atom species
  IF( LEN_TRIM(rmatom_prop)<=2 ) THEN
    species = TRIM(ADJUSTL(rmatom_prop))
    CALL ATOMNUMBER(species,snumber)
  ENDIF
  !
  IF( NINT(snumber)>0 ) THEN
    method=2
  ELSE
    species = ''
    !rmatom_prop does not contain an atom spcies
    !=> it should be an atom index, check for that
    READ(rmatom_prop,*,END=120,ERR=120) atomindex
    method=3
  ENDIF
  !
  120 CONTINUE
  CALL ATOMSK_MSG(2102,(/species/),(/DBLE(atomindex)/))
ENDIF
!
!
!
200 CONTINUE
IF( method==1 ) THEN
  !Remove all selected atoms
  NP=0
  IF( ALLOCATED(SELECT) ) THEN
    DO i=1,SIZE(SELECT)
      IF( .NOT.SELECT(i) ) NP=NP+1
    ENDDO
    !
    IF(NP>0) THEN
      ALLOCATE( newP(NP,4) )
      newP(:,:) = 0.d0
      IF( ALLOCATED(S) ) THEN
        ALLOCATE( newS(NP,4) )
        newS(:,:) = 0.d0
      ENDIF
      IF( ALLOCATED(AUX) ) THEN
        ALLOCATE( newAUX(NP,SIZE(AUX,2)) )
        newAUX(:,:) = 0.d0
      ENDIF
      !
      NP = 0
      DO i=1,SIZE(P,1)
        IF( .NOT.SELECT(i) ) THEN
          !This atom is not selected => it lives
          NP=NP+1
          newP(NP,:) = P(i,:)
          IF( ALLOCATED(S) ) newS(NP,:) = S(i,:)
          IF( ALLOCATED(AUX) ) newAUX(NP,:) = AUX(i,:)
        ELSE
          !This atom is selected => it will be removed
          rmatoms = rmatoms+1
        ENDIF
      ENDDO
      !
    ELSE
      !All atoms must be removed
      rmatoms = SIZE(P,1)
    ENDIF
    !
    !Since all selected atoms were removed, clear selection
    clear_sel = .TRUE.
    !
  ELSE
    !No selection is defined => do not remove any atom
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2752,(/species/),(/DBLE(atomindex)/))
  ENDIF
  !
  !
ELSEIF( method==2 ) THEN
  !Remove all atoms of the given species
  CALL ATOMSPECIES(snumber,species)
  !
  !Allocate arrays
  !For now, use the same size as the original arrays because
  !we don't know how many atoms will eventually be removed
  ALLOCATE( newP( SIZE(P,1), 4 ) )
  newP(:,:) = 0.d0
  IF(ALLOCATED(S)) THEN
    ALLOCATE( newS( SIZE(S,1), 4 ) )
    newS(:,:) = 0.d0
  ENDIF
  IF(ALLOCATED(AUX)) THEN
    ALLOCATE( newAUX( SIZE(AUX,1), SIZE(AUX,2) ) )
    newAUX(:,:) = 0.d0
  ENDIF
  IF(ALLOCATED(SELECT)) THEN
    ALLOCATE( newSELECT( SIZE(P,1) ) )
    newSELECT(:) = .FALSE.
  ENDIF
  !
  !Remove all atoms that belong to that species *and* are selected
  DO i=1,SIZE(P,1)
    IF( NINT(P(i,4))==NINT(snumber) .AND. (.NOT.ALLOCATED(SELECT) .OR. SELECT(i)) ) THEN
      !This atom dies
      rmatoms=rmatoms+1
      !
    ELSE
      !This atom lives
      NP=NP+1
      newP(NP,:) = P(i,:)
      IF(ALLOCATED(S)) THEN
        newS(NP,:) = S(i,:)
      ENDIF
      IF(ALLOCATED(AUX)) THEN
        newAUX(NP,:) = AUX(i,:)
      ENDIF
      IF(ALLOCATED(SELECT)) THEN
        newSELECT(NP) = SELECT(i)
      ENDIF
    ENDIF
    !
  ENDDO
  !
  !
ELSEIF( method==3 ) THEN
  !Remove atom with the given index
  IF( atomindex>SIZE(P,1) .OR. atomindex<=0 ) THEN
    !Atom index is out-of-bounds => warn and skip
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2742,(/""/),(/DBLE(atomindex)/))
    !
  ELSE
    !Allocate arrays
    ALLOCATE( newP( SIZE(P,1)-1, 4 ) )
    newP(:,:) = 0.d0
    IF(ALLOCATED(S)) THEN
      ALLOCATE( newS( SIZE(S,1)-1, 4 ) )
      newS(:,:) = 0.d0
    ENDIF
    IF(ALLOCATED(AUX)) THEN
      ALLOCATE( newAUX( SIZE(AUX,1)-1, SIZE(AUX,2) ) )
    ENDIF
    IF(ALLOCATED(SELECT)) THEN
      ALLOCATE( newSELECT( SIZE(newP,1) ) )
      newSELECT(:) = .FALSE.
    ENDIF
    !
    !Remove only atom with that index (but only if it is selected)
    IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(atomindex) ) THEN
      !Atom that must be removed is in the selection => proceed
      DO i=1,SIZE(P,1)
        IF(i==atomindex) THEN
          !This atom dies
          rmatoms=rmatoms+1
          CALL ATOMSPECIES(P(i,4),species)
        ELSE
          !This atom lives
          NP=NP+1
          newP(NP,:) = P(i,:)
          IF(ALLOCATED(S)) THEN
            newS(NP,:) = S(i,:)
          ENDIF
          IF(ALLOCATED(AUX)) THEN
            newAUX(NP,:) = AUX(i,:)
          ENDIF
          IF(ALLOCATED(SELECT)) THEN
            newSELECT(NP) = SELECT(i)
          ENDIF
        ENDIF
      ENDDO
    ENDIF
    !
  ENDIF
  !
  !
ELSE
  !Unable to determine which atoms to remove
  nerr=nerr+1
  CALL ATOMSK_MSG(2812,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
300 CONTINUE
!Replace old P by newP; same with S, AUX and SELECT if necessary
IF( rmatoms>=SIZE(P,1) ) THEN
  !All atoms were removed! Just de-allocate arrays
  DEALLOCATE(P)
  IF(ALLOCATED(S)) DEALLOCATE(S)
  IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
  !
ELSEIF( rmatoms>0 ) THEN
  DEALLOCATE(P)
  IF(ALLOCATED(S)) DEALLOCATE(S)
  IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
  IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
  !
  ALLOCATE(P(NP,4))
  P(:,:) = 0.d0
  IF(ALLOCATED(newS)) THEN
    ALLOCATE(S(NP,4))
    S(:,:) = 0.d0
  ENDIF
  IF(ALLOCATED(newAUX)) THEN
    ALLOCATE( AUX( NP,SIZE(newAUX,2) ) )
    AUX(:,:) = 0.d0
  ENDIF
  IF(ALLOCATED(newSELECT)) THEN
    ALLOCATE( SELECT(NP) )
    SELECT(:) = .FALSE.
  ENDIF
  !
  DO i=1,NP
    P(i,:) = newP(i,:)
    IF(ALLOCATED(S)) S(i,:) = newS(i,:)
    IF(ALLOCATED(AUX)) AUX(i,:) = newAUX(i,:)
    IF(ALLOCATED(SELECT)) SELECT(i) = newSELECT(i)
  ENDDO
ENDIF
IF(ALLOCATED(newP)) DEALLOCATE(newP)
IF(ALLOCATED(newS)) DEALLOCATE(newS)
IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
!
!
CALL ATOMSK_MSG(2080,(/species/),(/DBLE(rmatoms),DBLE(SIZE(P,1))/))
IF( clear_sel ) THEN
  !All selected atoms were removed => clear the selection
  nwarn=nwarn+1
  CALL ATOMSK_MSG(2750,(/''/),(/0.d0/))
  IF( ALLOCATED(SELECT) ) DEALLOCATE(SELECT)
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
!
!
!
END SUBROUTINE RMATOM_XYZ
!
END MODULE rmatom