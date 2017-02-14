MODULE out_pdb
!
!
!**********************************************************************************
!*  OUT_PDB                                                                       *
!**********************************************************************************
!* This module writes a file in the Protein Data Bank (PDB) format.               *
!* The produced file is only a draft though, and does not fully comply            *
!* with the PDB standard, so proceed with care.                                   *
!* The PDB format is officially described here:                                   *
!*     http://www.wwpdb.org/docs.html                                             *
!**********************************************************************************
!* (C) Oct. 2012 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 14 Feb. 2017                                     *
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
SUBROUTINE WRITE_PDB(H,P,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=1):: atom_altLoc, atom_chainID, atom_iCode
CHARACTER(LEN=2):: atom_element, atom_charge
CHARACTER(LEN=1):: atom_resName
CHARACTER(LEN=4):: atom_name
CHARACTER(LEN=80):: pdbline  !exactly the length of a line
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: isreduced
INTEGER:: atom_serial, atom_resSeq
INTEGER:: i, j
INTEGER:: q, occ  !position of properties in AUX: atom charge (q), occupancy
REAL(dp):: a, b, c, alpha, beta, gamma !supercell (conventional notation)
REAL(dp):: atom_occupancy, atom_tempFactor
REAL(dp):: smass
REAL(dp),DIMENSION(3):: TN, UN
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIGXN, SCALEN   !Origin and scale matrices
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
isreduced = .FALSE.
occ=0
q=0
TN(:) = 0.d0
UN(:) = 0.d0
ORIGXN(:,:)=0.d0
SCALEN(:,:) = 0.d0
DO i=1,3
  ORIGXN(i,i) = 1.d0
  SCALEN(i,i) = 1.d0
ENDDO
IF(ALLOCATED(aentries)) DEALLOCATE(aentries)
!
msg = 'entering WRITE_PDB'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
IF( .NOT.ALLOCATED(P) .OR. SIZE(P,1)==0 ) THEN
  nerr = nerr+1
  CALL ATOMSK_MSG(3800,(/TRIM(msg)/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
IF( ALLOCATED(AUXNAMES) ) THEN
  DO i=1,SIZE(AUXNAMES)
    IF( TRIM(ADJUSTL(AUXNAMES(i)))=='occupancy' ) THEN
      occ = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='q' ) THEN
      q = i
    ENDIF
  ENDDO
ENDIF
!
!
!
100 CONTINUE
!Check if coordinates are reduced or cartesian
CALL FIND_IF_REDUCED(P,isreduced)
!
OPEN(UNIT=40,FILE=outputfile,STATUS='UNKNOWN',ERR=500)
!
!  ===  TITLE SECTION  ===
!Lines for this section may be stored in the array comment(:).
!*BUT* atomsk has added a hash sign (#) at the beginning of each comment(:),
!so we have to look for characters 2:7.
!If a keyword cannot be found in comment(:) then produce dummy lines.
!
!Header (remove leading # from comment)
pdbline = 'HEADER   '//comment(1)(2:73)
IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
  DO i=1,SIZE(comment)
    IF( comment(i)(2:7)=='HEADER' ) THEN
      pdbline = comment(i)(2:80)
      EXIT
    ENDIF
  ENDDO
ENDIF
WRITE(40,'(a80)') pdbline
!
!Title (remove leading # from comment)
pdbline = 'TITLE    '//comment(1)(2:73)
IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
  DO i=1,SIZE(comment)
    IF( comment(i)(2:7)=='TITLE' ) THEN
      pdbline = comment(i)(2:80)
      EXIT
    ENDIF
  ENDDO
ENDIF
WRITE(40,'(a80)') pdbline
!
!Caveat: this is to warn that this file does not fully comply to the PDB standard
WRITE(pdbline,'(a6)') 'CAVEAT      DRAFT FILE PRODUCED WITH ATOMSK'
WRITE(40,'(a80)') pdbline
!
!Compound
pdbline=''
IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
  DO i=1,SIZE(comment)
    IF( comment(i)(2:7)=='COMPND' ) THEN
      pdbline = comment(i)(2:80)
      EXIT
    ENDIF
  ENDDO
ENDIF
IF( LEN_TRIM(pdbline)==0 ) THEN
  !Find how many atoms of each species exist, and write a formula
  temp = ''
  WRITE(pdbline,'(a)') 'COMPND '//temp(1:73)
ENDIF
WRITE(40,'(a80)') pdbline
!
!Source
pdbline=''
IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
  DO i=1,SIZE(comment)
    IF( comment(i)(2:7)=='SOURCE' ) THEN
      pdbline = comment(i)(2:80)
      EXIT
    ENDIF
  ENDDO
ENDIF
IF( LEN_TRIM(pdbline)==0 ) THEN
  WRITE(pdbline,'(a)') 'SOURCE '
ENDIF
WRITE(40,'(a80)') pdbline
!
!Keywords
pdbline=''
IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
  DO i=1,SIZE(comment)
    IF( comment(i)(2:7)=='KEYWDS' ) THEN
      pdbline = comment(i)(2:80)
      EXIT
    ENDIF
  ENDDO
ENDIF
IF( LEN_TRIM(pdbline)==0 ) THEN
  WRITE(pdbline,'(a)') 'KEYWDS '
ENDIF
WRITE(40,'(a80)') pdbline
!
!Experimental data
pdbline=''
IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
  DO i=1,SIZE(comment)
    IF( comment(i)(2:7)=='EXPDTA' ) THEN
      pdbline = comment(i)(2:80)
      EXIT
    ENDIF
  ENDDO
ENDIF
IF( LEN_TRIM(pdbline)==0 ) THEN
  WRITE(pdbline,'(a)') 'EXPDTA '
ENDIF
WRITE(40,'(a80)') pdbline
!
!Author
pdbline=''
IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
  DO i=1,SIZE(comment)
    IF( comment(i)(2:7)=='AUTHOR' ) THEN
      pdbline = comment(i)(2:80)
      EXIT
    ENDIF
  ENDDO
ENDIF
IF( LEN_TRIM(pdbline)==0 ) THEN
  WRITE(pdbline,'(a)') 'AUTHOR '
ENDIF
WRITE(40,'(a80)') pdbline
!
!Data revision
pdbline=''
IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
  DO i=1,SIZE(comment)
    IF( comment(i)(2:7)=='REVDAT' ) THEN
      pdbline = comment(i)(2:80)
      EXIT
    ENDIF
  ENDDO
ENDIF
IF( LEN_TRIM(pdbline)==0 ) THEN
  WRITE(pdbline,'(a)') 'REVDAT '
ENDIF
WRITE(40,'(a80)') pdbline
!
!
!  ===  REMARKS SECTION  ===
pdbline=''
IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
  DO i=1,SIZE(comment)
    IF( comment(i)(2:7)=='REMARK' ) THEN
      pdbline = comment(i)(2:80)
      WRITE(40,'(a80)') pdbline
    ENDIF
  ENDDO
ENDIF
IF( LEN_TRIM(pdbline)==0 ) THEN
  WRITE(pdbline,'(a)') 'REMARK '
  WRITE(40,'(a80)') pdbline
ENDIF
!
!
!  ===  PRIMARY STRUCTURE  ===
WRITE(40,'(a6)') 'SEQRES'
!
!
!  ===  CRYSTALLOGRAPHIC  ===
!Supercell parameters
CALL MATCONV(H,a,b,c,alpha,beta,gamma)

WRITE(pdbline,601) 'CRYST1', a, b, c, RAD2DEG(alpha), RAD2DEG(beta), RAD2DEG(gamma), 'P 1       ', 1
WRITE(40,'(a80)') pdbline
!
!
!  ===  COORDINATE TRANSFORMATION  ===
WRITE(pdbline,160) 'ORIGX1', (ORIGXN(1,j),j=1,3), TN(1)
WRITE(40,'(a80)') pdbline
WRITE(pdbline,160) 'ORIGX2', (ORIGXN(2,j),j=1,3), TN(1)
WRITE(40,'(a80)') pdbline
WRITE(pdbline,160) 'ORIGX3', (ORIGXN(3,j),j=1,3), TN(1)
WRITE(40,'(a80)') pdbline
WRITE(pdbline,160) 'SCALE1', (SCALEN(1,j),j=1,3), UN(1)
WRITE(40,'(a80)') pdbline
WRITE(pdbline,160) 'SCALE2', (SCALEN(2,j),j=1,3), UN(1)
WRITE(40,'(a80)') pdbline
WRITE(pdbline,160) 'SCALE3', (SCALEN(3,j),j=1,3), UN(1)
WRITE(40,'(a80)') pdbline
160 FORMAT(a6,5X,3f10.6,6X,f10.5)
!
!
!  ===  COORDINATE  ===
!Write atomic positions
DO i=1,SIZE(P,1)
  CALL ATOMSPECIES(P(i,4),atom_element)
  !PDB uses only upper-case symbols...
  atom_element = STRUPCASE(atom_element)
  !
  atom_name = ADJUSTL(atom_element)
  atom_altLoc = ''    !alternate location indicator
  atom_resName = ''   !residue name
  atom_chainID = ''   !chain identifier
  atom_resSeq = 0     !residue sequence number
  atom_iCode = ''     !code for insertion of residue
  !
  !Set atom occupancy
  IF(occ>0) THEN
    atom_occupancy = AUX(i,occ)
  ELSE
    atom_occupancy = 1.d0
  ENDIF
  !
  !Set atom temperature factor
  atom_tempFactor = 0.d0
  !
  !Set atom charge
  atom_charge = "  "
  IF(q>0) THEN
    IF( NINT(AUX(i,q)).NE.0 ) THEN
      WRITE(atom_charge(1:1),'(i1)') NINT(AUX(i,q))
      IF( AUX(i,q)<0.d0 ) THEN
        atom_charge(2:2) = '-'
      ELSE
        atom_charge(2:2) = '+'
      ENDIF
    ENDIF
  ENDIF
  !
  !Set atom element name
  !it must be right-justified
  atom_element=ADJUSTR(atom_element)
  !
  !Write line to file
  pdbline(1:6) = "ATOM  "
  WRITE(pdbline(7:11),'(i5)') i
  pdbline(13:16) = TRIM(atom_name)
  pdbline(17:17) = atom_altLoc
  pdbline(18:20) = TRIM(atom_resName)
  pdbline(22:22) = atom_chainID
  WRITE(pdbline(23:26),'(i4)') atom_resSeq
  pdbline(27:27) = atom_iCode
  WRITE(pdbline(31:38),'(f8.3)') P(i,1)
  WRITE(pdbline(39:46),'(f8.3)') P(i,2)
  WRITE(pdbline(47:54),'(f8.3)') P(i,3)
  WRITE(pdbline(55:60),'(f6.2)') atom_occupancy
  WRITE(pdbline(61:66),'(f6.2)') atom_tempFactor
  pdbline(77:78) = ADJUSTR(atom_element)
  WRITE(pdbline(79:80),'(a2)') atom_charge
  WRITE(40,'(a80)') pdbline
ENDDO
WRITE(pdbline,'(a6)') 'TER   '
WRITE(40,'(a80)') pdbline
!
!
!  ===  BOOKKEEPING  ===
WRITE(pdbline,'(a6)') 'MASTER'
WRITE(40,'(a80)') pdbline
WRITE(pdbline,'(a3)') 'END'
WRITE(40,'(a80)') pdbline
!
!
!
500 CONTINUE
CLOSE(40)
msg = "PDB"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
!
!
!
!The standard line formats of PDB files
!ATOM
600 FORMAT(a6,i5,1X,a4,a1,a3,1X,a1,i4,a1,3X,3f8.3,2f6.2,2X,2a2)
!CRYST1
601 FORMAT(a6,3f9.3,3f7.2,1X,a10,i3)
!
!
!
1000 CONTINUE
IF(ALLOCATED(aentries)) DEALLOCATE(aentries)
!
END SUBROUTINE WRITE_PDB
!
END MODULE out_pdb
