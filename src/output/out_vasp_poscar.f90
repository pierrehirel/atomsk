MODULE out_vasp_poscar
!
!
!**********************************************************************************
!*  OUT_VASP_POSCAR                                                               *
!**********************************************************************************
!* This module writes POSCAR files, as well as CONTCAR files since it             *
!* is the same format. These files are used by VASP.                              *
!* The POSCAR format is described here:                                           *
!*    http://cms.mpi.univie.ac.at/vasp/vasp/POSCAR_file.html                      *
!* Beware that atomic species are not specified in the POSCAR format,             *
!* hence that information is lost when converting from any format to POSCAR.      *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 29 July 2016                                     *
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
USE functions
USE messages
USE files
USE subroutines
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE WRITE_POSCAR(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=1):: answer
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: fixedatoms  !are some atoms fixed?
LOGICAL:: isreduced
INTEGER:: i, j, last
INTEGER:: fixx, fixy, fixz  !are atoms fixed along X, Y, Z? (0=no, 1=yes)
LOGICAL:: contiguous   !are atom with same species contiguous? must they be packed?
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: atypes
REAL(dp),DIMENSION(:,:),POINTER:: Ppoint  !pointer to P or Q
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN),TARGET:: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,TARGET:: Q  !copy of P (if needed)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN),TARGET:: AUX !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE,TARGET:: AUX2 !auxiliary properties (temporary)
REAL(dp),DIMENSION(:,:),POINTER:: AUXpoint  !pointer to AUX or AUX2
!
!
!Initialize variables
fixx=0
fixy=0
fixz=0
 contiguous = .TRUE.
 fixedatoms = .FALSE.
!
WRITE(msg,*) 'entering WRITE_POSCAR'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Check if some atoms are fixed
IF( ALLOCATED(AUXNAMES) .AND. SIZE(AUXNAMES)>0 ) THEN
  DO i=1,SIZE(AUXNAMES)
    IF(AUXNAMES(i)=="fixx") fixx=i
    IF(AUXNAMES(i)=="fixy") fixy=i
    IF(AUXNAMES(i)=="fixz") fixz=i
  ENDDO
ENDIF
IF( fixx>0 .AND. fixy>0 .AND. fixz>0 ) THEN
  fixedatoms = .TRUE.
ENDIF
!
!Check if coordinates are reduced
CALL FIND_IF_REDUCED(P,isreduced)
!
!
!
100 CONTINUE
!First make sure that atoms of same species are all contiguous
DO i=1,SIZE(P,1)
  last=i
  DO j=i+1,SIZE(P,1)
    IF( P(j,4)==P(i,4) ) THEN
      IF(j==last+1) THEN
        last=j
      ELSE
        contiguous=.FALSE.
        GOTO 120
      ENDIF
    ENDIF
  ENDDO
ENDDO
!
120 CONTINUE
IF(.NOT.contiguous) THEN
  !Display a message
  nwarn=nwarn+1
  CALL ATOMSK_MSG(3703,(/msg,outputfile/),(/0.d0/))
  READ(*,*) answer
  IF(answer.NE.langyes .AND. answer.NE.langBigYes) THEN
    contiguous = .TRUE.
  ENDIF
ENDIF
!
IF(.NOT.contiguous) THEN
  !If species are not contiguous then rearrange them
  !Note that the arrays P and AUX must remained untouched
  !Therefore the columns fixx, fixy, fixz of AUX are copied into Q for sorting
  IF( fixedatoms ) THEN
    ALLOCATE( Q( SIZE(P,1), 7 ) )
    Q(:,1:4) = P(:,1:4)
    Q(:,5) = AUX(:,fixx)
    Q(:,6) = AUX(:,fixy)
    Q(:,7) = AUX(:,fixz)
  ELSE
    ALLOCATE( Q( SIZE(P,1), SIZE(P,2) ) )
    Q(:,:) = P(:,:)
  ENDIF
  !
  !Sort atoms (and aux. prop.) to pack atoms with the same species
  CALL PACKSORT(Q,4)
  Ppoint=>Q
  !
  IF( ALLOCATED(AUX) .AND. SIZE(AUX,1)>0 ) THEN
    !Save which atoms are fixed in AUX2
    ALLOCATE( AUX2( SIZE(AUX,1) , 3 ) )
    AUX2(:,1) = Q(:,5)
    AUX2(:,2) = Q(:,6)
    AUX2(:,3) = Q(:,7)
    fixx = 1
    fixy = 2
    fixz = 3
    AUXpoint=>AUX2
  ENDIF
  !
  !Count the number of atoms for each atomic species
  CALL FIND_NSP(Ppoint(:,4),atypes)
  !
  CALL ATOMSK_MSG(3003,(/''/),(/atypes(:,1)/))
  !
ELSE
  !If everything is fine, just use P in the following
  Ppoint=>P
  AUXpoint=>AUX
  !Count the number of atoms for each atomic species
  CALL FIND_NSP(Ppoint(:,4),atypes)
ENDIF
!
!
!
200 CONTINUE
WRITE(msg,*) "Number of different species: ", SIZE(atypes(:,1))
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN',ERR=500)
!
!Write header of POSCAR file
WRITE(40,*) TRIM(ADJUSTL(comment(1)))
WRITE(40,'(a8)') '1.000000'
WRITE(40,201) H(1,1), H(1,2), H(1,3)
WRITE(40,201) H(2,1), H(2,2), H(2,3)
WRITE(40,201) H(3,1), H(3,2), H(3,3)
201 FORMAT(3(f16.8,2X))
!
!!!  VASP 5.x: write the element symbol for each species
!!!  Uncomment the following lines to enable this feature
!!!  NOTE: this makes the POSCAR file incompatible with VASP 4.x !!!!!
!temp=""
!DO i=1,SIZE(atypes,1)
!  CALL ATOMSPECIES(atypes(i,1),species)
!  temp = TRIM(ADJUSTL(temp))//"  "//species
!ENDDO
!WRITE(40,*) TRIM(ADJUSTL(temp))
!!!  END OF VASP 5.x
!
!Write the number of atoms for each species
WRITE(40,'(20(i6,2X))') ( NINT(atypes(j,2)), j=1,SIZE(atypes(:,1)) )
!
IF( fixx.NE.0 .OR. fixy.NE.0 .OR. fixz.NE.0 ) THEN
  WRITE(40,'(a18)') 'Selective dynamics'
ENDIF
!
!Write atom coordinates
IF(isreduced) THEN
  WRITE(40,'(a6)') 'Direct'
ELSE
  WRITE(40,'(a9)') 'Cartesian'
ENDIF
DO i=1,SIZE(Ppoint,1)
  WRITE(msg,210) Ppoint(i,1), Ppoint(i,2), Ppoint(i,3)
  IF( fixedatoms ) THEN
    !Caution: internally if AUX(fix)==1 then atom is fixed.
    !On the contrary in VASP the flag "F" means that atom is fixed, "T" that it's mobile
    IF( AUXpoint(i,fixx)>0.5d0 ) THEN
      msg = TRIM(msg)//" F"
    ELSE
      msg = TRIM(msg)//" T"
    ENDIF
    IF( AUXpoint(i,fixy)>0.5d0 ) THEN
      msg = TRIM(msg)//" F"
    ELSE
      msg = TRIM(msg)//" T"
    ENDIF
    IF( AUXpoint(i,fixz)>0.5d0 ) THEN
      msg = TRIM(msg)//" F"
    ELSE
      msg = TRIM(msg)//" T"
    ENDIF
  ENDIF
  WRITE(40,'(a)') TRIM(msg)
ENDDO
210 FORMAT(3(f16.8,2X))
GOTO 500
!
!
!
500 CONTINUE
CLOSE(40)
msg = "POSCAR"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
1000 CONTINUE
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(AUX2)) DEALLOCATE(AUX2)
NULLIFY(Ppoint)
NULLIFY(AUXpoint)
!
!
!
END SUBROUTINE WRITE_POSCAR
!
END MODULE out_vasp_poscar
