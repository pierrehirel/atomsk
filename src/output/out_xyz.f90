MODULE out_xyz
!
!
!**********************************************************************************
!*  OUT_XYZ                                                                       *
!**********************************************************************************
!* This module writes an XYZ file. It can also write an "extended"                *
!* XYZ format, i.e. the comment line (2nd line) is replaced by a series           *
!* of keywords/values, for instance:                                              *
!*  NP                                                                            *
!*  Lattice="H11 H21 H31 H12 H22 H32 H13 H23 H33" Properties=species:S:1:pos:R:3  *
!*  species1 x1 y1 z1                                                             *
!*  species2 x2 y2 z2                                                             *
!*  ...                                                                           *
!* This module can also write a "special" XYZ format, which instead adds          *
!* information at the end of the file:                                            *
!*  NP                                                                            *
!*  #comment line                                                                 *
!*  species1 x1 y1 z1                                                             *
!*  species2 x2 y2 z2                                                             *
!*   ...                                                                          *
!*  alat                                                                          *
!*  <a0>                                                                          *
!*  supercell                                                                     *
!*  <H(1,1)> <H(1,2)> <H(1,3)>                                                    *
!*  <H(2,1)> <H(2,2)> <H(2,3)>                                                    *
!*  <H(3,1)> <H(3,2)> <H(3,3)>                                                    *
!*  <cartesian/fractional> coordinates                                            *
!*  mass <species1> <smass1>                                                      *
!*  mass <species2> <smass2>                                                      *
!*   ...                                                                          *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 02 June 2022                                     *
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
SUBROUTINE WRITE_XYZ(H,P,comment,AUXNAMES,AUX,outputfile,xyzformatin)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=5),INTENT(IN),OPTIONAL:: xyzformatin !'exyz' for Extended XYZ, 'sxyz' for Special XYZ, otherwise normal XYZ
CHARACTER(LEN=5):: xyzformat
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: isreduced
INTEGER:: i, j
REAL(dp):: smass
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
IF(PRESENT(xyzformatin)) THEN
  xyzformat = xyzformatin
ELSE
  xyzformat = 'xyz'
ENDIF
isreduced = .FALSE.
IF(ALLOCATED(aentries)) DEALLOCATE(aentries)
!
msg = 'entering WRITE_XYZ'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
IF( .NOT.ALLOCATED(P) .OR. SIZE(P,1)==0 ) THEN
  nerr = nerr+1
  CALL ATOMSK_MSG(3800,(/TRIM(msg)/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
100 CONTINUE
!Check if coordinates are reduced or cartesian
CALL FIND_IF_REDUCED(H,P,isreduced)
!
IF(ofu.NE.6) THEN
  OPEN(UNIT=ofu,FILE=outputfile,STATUS='UNKNOWN',ERR=500)
ENDIF
!First line: number of particles
WRITE(msg,*) SIZE(P(:,1))
WRITE(ofu,'(a)') TRIM(ADJUSTL(msg))
!
!Second line: keywords or comment
IF(xyzformat=='exyz' .OR. xyzformat=='EXYZ') THEN
  !For extended XYZ, replace comment by keywords
  msg=''
  temp=''
  DO i=1,3
    DO j=1,3
      IF(H(i,j)==0.d0) THEN
        WRITE(temp,*) '0.0'
      ELSE
        WRITE(temp,'(f16.8)') H(i,j)
      ENDIF
      WRITE(msg,*) TRIM(ADJUSTL(msg))//' '//TRIM(ADJUSTL(temp))
    ENDDO
  ENDDO
  temp = 'Lattice="'//TRIM(ADJUSTL(msg))//'" '//'Properties=species:S:1:pos:R:3'
  !If auxiliary properties are present concatenate their names here
  IF( ALLOCATED(AUXNAMES) ) THEN
    DO i=1,SIZE(AUXNAMES)
      temp = TRIM(temp)//':'//TRIM(ADJUSTL(AUXNAMES(i)))//':R:1'
    ENDDO
  ENDIF
  WRITE(ofu,'(a)') TRIM(temp)
ELSE
  !Otherwise just write the comment
  IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
    WRITE(ofu,'(a)') TRIM(comment(1))
  ELSE
    WRITE(ofu,'(a)') "#"
  ENDIF
ENDIF
!
!Write atomic positions
DO i=1,SIZE(P(:,1))
  CALL ATOMSPECIES(P(i,4),species)
  WRITE(temp,*) species
  temp = TRIM(ADJUSTL(temp))
  DO j=1,3
    WRITE(msg,'(f16.8)')  P(i,j)
    temp = TRIM(ADJUSTL(temp))//'  '//TRIM(ADJUSTL(msg))
  ENDDO
  !If auxiliary properties are present, concatenate them at end of line
  !This is done only for extended and special XYZ
  IF( ALLOCATED(AUXNAMES) .AND.                                                                  &
      ( xyzformat=='sxyz' .OR. xyzformat=='SXYZ' .OR. xyzformat=='exyz' .OR. xyzformat=='EXYZ' ) &
    ) THEN
    DO j=1,SIZE(AUX(1,:))
      WRITE(msg,'(e12.5)')  AUX(i,j)
      temp = TRIM(ADJUSTL(temp))//'  '//TRIM(ADJUSTL(msg))
    ENDDO
  ENDIF
  WRITE(ofu,'(a)') TRIM(ADJUSTL(temp))
ENDDO
!
!For special XYZ, write footer of XYZ file
IF(xyzformat=='sxyz' .OR. xyzformat=='SXYZ') THEN
  WRITE(ofu,'(a4)') 'alat  1.0'
  WRITE(ofu,'(a3)') '1.0'
  WRITE(ofu,'(a9)') 'supercell'
  WRITE(ofu,'(3(f12.6,1X))') H(1,1), H(1,2), H(1,3)
  WRITE(ofu,'(3(f12.6,1X))') H(2,1), H(2,2), H(2,3)
  WRITE(ofu,'(3(f12.6,1X))') H(3,1), H(3,2), H(3,3)
  !
  !Masses of the atoms
  CALL FIND_NSP(P(:,4),aentries)
  DO i=1,SIZE(aentries(:,1))
    CALL ATOMSPECIES(DBLE(aentries(i,1)),species)
    CALL ATOMMASS(species,smass)
    WRITE(msg,'(f12.3)') smass
    WRITE(msg,*) 'mass '//species//' '//TRIM(ADJUSTL(msg))
    WRITE(ofu,'(a)') TRIM(ADJUSTL(msg))
  ENDDO
  !
  !Names of auxiliary properties
  IF( ALLOCATED(AUXNAMES) ) THEN
    DO i=1,SIZE(AUXNAMES)
      WRITE(temp,*) i
      temp = 'property '//TRIM(ADJUSTL(temp))//' '//AUXNAMES(i)
      WRITE(ofu,'(a)') TRIM(ADJUSTL(temp))
    ENDDO
  ENDIF
  !
  IF(isreduced) THEN
    WRITE(ofu,'(a19)') 'reduced coordinates'
  ELSE
    WRITE(ofu,'(a21)') 'cartesian coordinates'
  ENDIF
ENDIF
!
!
!
500 CONTINUE
IF(ofu.NE.6) THEN
  CLOSE(ofu)
ENDIF
IF(xyzformat=='exyz' .OR. xyzformat=='EXYZ') THEN
  msg = "extended XYZ"
ELSEIF(xyzformat=='sxyz' .OR. xyzformat=='SXYZ') THEN
  msg = "special XYZ"
ELSE
  msg = "XYZ"
ENDIF
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
1000 CONTINUE
IF(ALLOCATED(aentries)) DEALLOCATE(aentries)
!
END SUBROUTINE WRITE_XYZ
!
END MODULE out_xyz
