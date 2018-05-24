MODULE in_qe_pw
!
!
!**********************************************************************************
!*  IN_QE_PW                                                                      *
!**********************************************************************************
!* This module reads pw input files for Quantum Espresso.                         *
!* This file format is described here:                                            *
!*    http://www.quantum-espresso.org/wp-content/uploads/Doc/INPUT_PW.html        *
!**********************************************************************************
!* (C) June 2012 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 16 May 2018                                      *
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
USE qepw_ibrav
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE READ_QEPW(inputfile,H,P,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=1):: cell_units, atpos_units !units of cell vectors, atom positions: A (angstroms) or B (Bohrs)
CHARACTER(LEN=2):: species
CHARACTER(LEN=16):: section  !section of the PW file
CHARACTER(LEN=128):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: fields
LOGICAL:: convnot !are conventional vectors a,b,c,alpha,beta,gamma defined?
LOGICAL:: celldm_defined !are the celldm(:) defined?
INTEGER:: fx, fy, fz       !position of forces (x,y,z) in AUX
INTEGER:: fixx, fixy, fixz !position of flags for fixed atoms in AUX
INTEGER:: i
INTEGER:: ibrav  !index of Bravais lattice
INTEGER:: Naux   !number of auxiliary properties
INTEGER:: NP     !total number of particles in the system
INTEGER:: strlength
REAL(dp):: a, b, c, cosalpha, cosbeta, cosgamma
REAL(dp):: tempreal
REAL(dp),DIMENSION(6):: celldm
REAL(dp),DIMENSION(3,3),INTENT(OUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: AUX !auxiliary properties
!
!
!Initialize variables
atpos_units='B'
 cell_units='B'
section = ''
 convnot=.FALSE.
 celldm_defined=.FALSE.
ibrav=0
fx=0
fy=0
fz=0
fixx=0
fixy=0
fixz=0
Naux=0
NP=0
i=0
a = 0.d0
b = 0.d0
 c = 0.d0
 cosalpha = 0.d0   !default: alphe=beta=gamma=90°
 cosbeta = 0.d0
 cosgamma = 0.d0
 celldm(:)=0.d0
H(:,:) = 0.d0
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
!
WRITE(msg,*) 'entering READ_QEPW'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
!
100 CONTINUE
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
!Parse the file for number of atoms and cell parameters
DO
  IF(ALLOCATED(fields)) DEALLOCATE(fields)
  READ(30,'(a128)',END=150,ERR=150) temp
  temp=ADJUSTL(temp)
  CALL ATOMSK_MSG(999,(/temp/),(/0.d0/))
  !
  IF( temp(1:1).NE."!" .AND. temp(1:1).NE."#" ) THEN
    !Check if there are several keywords (they should be separated by commas)
    strlength=SCAN(temp,',')
    IF(strlength<=0) THEN
      !only one keyword/value on this line
      ALLOCATE(fields(1))
      fields(1) = temp
    ELSE
      !several keywords/values on this line, separate them
      i=0
      ALLOCATE(fields(20))  !assuming there are no more than 20 values on a line
      fields(:)=''
      DO WHILE(LEN_TRIM(temp)>0)
        i=i+1
        IF(i>20) GOTO 800
        strlength = SCAN(temp,",")
        IF(strlength==0) EXIT
        fields(i) = TRIM(ADJUSTL(temp(:strlength-1)))
        temp = temp(strlength+1:)
      ENDDO
      !
      IF( section=="atomic_positions" .AND. fixx==0 .AND. fixy==0 .AND. fixz==0 ) THEN
        IF( i==7 ) THEN
          !This section contains 7 fields:  X x y z if_pos(1) if_pos(2) if_pos(3)
          !Set indexes to store fixed coordinates into AUX later
          fixx=Naux+1
          fixy=Naux+2
          fixz=Naux+3
          Naux=Naux+3
        ENDIF
      ENDIF
    ENDIF
    !
    DO i=1,SIZE(fields)
      IF(LEN_TRIM(fields(i))>0) THEN
        CALL ATOMSK_MSG(999,(/fields(i)/),(/0.d0/))
        !
        IF( fields(i)(1:5) == '&cell' .OR. fields(i)(1:5) == '&CELL' ) THEN
          section = "cell"
        ELSEIF( fields(i)(1:8) == '&control' .OR. fields(i)(1:8) == '&CONTROL' ) THEN
          section = "control"
        ELSEIF( fields(i)(1:10) == '&electrons' .OR. fields(i)(1:10) == '&ELECTRONS' ) THEN
          section = "electrons"
        ELSEIF( fields(i)(1:5) == '&ions' .OR. fields(i)(1:5) == '&IONS' ) THEN
          section = "ions"
        ELSEIF( fields(i)(1:8) == '&system' .OR. fields(i)(1:8) == '&SYSTEM' ) THEN
          section = "system"
        ELSEIF( fields(i)(1:2) == '/' ) THEN
          section = ""
        ELSEIF( fields(i)(1:13) == 'ATOMIC_FORCES' .OR. fields(i)(1:13) == 'atomic_forces' ) THEN
          section = "atomic_forces"
          fx=1
          fy=2
          fz=3
          Naux=Naux+3
        ELSEIF( temp(1:16)=='ATOMIC_POSITIONS' .OR. temp(1:16)=='atomic_positions' ) THEN
          section = "atomic_positions"
          !
        ELSEIF( section=='system' ) THEN
          !
          IF( fields(i)(1:3)=='nat' .OR. fields(i)(1:3)=='Nat' .OR. fields(i)(1:3)=='NAT' ) THEN
            !Read number of atoms NP
            strlength = SCAN(fields(i),'=')
            READ(fields(i)(strlength+1:),*,END=800,ERR=800) NP
          ELSEIF( fields(i)(1:5)=='ibrav' .OR. fields(i)(1:5)=='IBRAV' ) THEN
            !Defines the type of lattice: read it
            strlength = SCAN(fields(i),'=')
            READ(fields(i)(strlength+1:),*,END=800,ERR=800) ibrav
          ELSEIF( fields(i)(1:9)=='celldm(1)' .OR. fields(i)(1:9)=='CELLDM(1)' ) THEN
            !Read cell parameter alat
            strlength = SCAN(fields(i),'=')
            READ(fields(i)(strlength+1:),*,END=800,ERR=800) celldm(1)
          ELSEIF( fields(i)(1:9)=='celldm(2)' .OR. fields(i)(1:9)=='CELLDM(2)' ) THEN
            !Read cell parameter alat
            strlength = SCAN(fields(i),'=')
            READ(fields(i)(strlength+1:),*,END=800,ERR=800) celldm(2)
          ELSEIF( fields(i)(1:9)=='celldm(3)' .OR. fields(i)(1:9)=='CELLDM(3)' ) THEN
            !Read cell parameter alat
            strlength = SCAN(fields(i),'=')
            READ(fields(i)(strlength+1:),*,END=800,ERR=800) celldm(3)
          ELSEIF( fields(i)(1:9)=='celldm(4)' .OR. fields(i)(1:9)=='CELLDM(4)' ) THEN
            !Read cell parameter alat
            strlength = SCAN(fields(i),'=')
            READ(fields(i)(strlength+1:),*,END=800,ERR=800) celldm(4)
          ELSEIF( fields(i)(1:9)=='celldm(5)' .OR. fields(i)(1:9)=='CELLDM(5)' ) THEN
            !Read cell parameter alat
            strlength = SCAN(fields(i),'=')
            READ(fields(i)(strlength+1:),*,END=800,ERR=800) celldm(5)
          ELSEIF( fields(i)(1:9)=='celldm(6)' .OR. fields(i)(1:9)=='CELLDM(6)' ) THEN
            !Read cell parameter alat
            strlength = SCAN(fields(i),'=')
            READ(fields(i)(strlength+1:),*,END=800,ERR=800) celldm(6)
          ELSEIF( fields(i)(1:5)=='cosAB' .OR. fields(i)(1:5)=='cosab' ) THEN
            !Read cell angle gamma
            strlength = SCAN(fields(i),'=')
            READ(fields(i)(strlength+1:),*,END=800,ERR=800) cosgamma
          ELSEIF( fields(i)(1:5)=='cosAC' .OR. fields(i)(1:5)=='cosac' ) THEN
            !Read cell angle beta
            strlength = SCAN(fields(i),'=')
            READ(fields(i)(strlength+1:),*,END=800,ERR=800) cosbeta
          ELSEIF( fields(i)(1:5)=='cosBC' .OR. fields(i)(1:5)=='cosbc' ) THEN
            !Read cell angle alpha
            strlength = SCAN(fields(i),'=')
            READ(fields(i)(strlength+1:),*,END=800,ERR=800) cosalpha
          ELSEIF( fields(i)(1:5)=='angle' .OR. fields(i)(1:5)=='block' .OR. &
                & fields(i)(1:3)=='ass'  ) THEN
            !ignore these keywords
          ELSEIF( fields(i)(1:1)=='A' .OR. fields(i)(1:1)=='a' ) THEN
            !Read cell parameter A
            strlength = SCAN(fields(i),'=')
            IF( strlength>0 ) READ(fields(i)(strlength+1:),*,END=800,ERR=800) a
          ELSEIF( fields(i)(1:1)=='B' .OR. fields(i)(1:1)=='b' ) THEN
            !Read cell parameter B
            strlength = SCAN(fields(i),'=')
            IF( strlength>0 ) READ(fields(i)(strlength+1:),*,END=800,ERR=800) b
          ELSEIF( fields(i)(1:1)=='C' .OR. fields(i)(1:1)=='c' ) THEN
            !Read cell parameter C
            strlength = SCAN(fields(i),'=')
            IF( strlength>0 ) READ(fields(i)(strlength+1:),*,END=800,ERR=800) c
          ENDIF
          !
        ENDIF
      ENDIF
      !
    ENDDO !loop on i
    !
  ENDIF !end if temp=="!" .OR. temp=="#"
  !
ENDDO
!
150 CONTINUE
!At this point if we did not find the number of atoms we are in trouble:
!count them "manually"
IF( NP==0 ) THEN
  REWIND(30)
  NP=0
  DO
    READ(30,'(a128)',END=200,ERR=200) temp
    IF( temp(1:16)=='ATOMIC_POSITIONS' .OR. temp(1:16)=='atomic_positions' ) THEN
      DO
        !Attempt to read species + 3 real numbers until it fails
        READ(30,*,END=200,ERR=200) species, tempreal,tempreal,tempreal
        NP=NP+1
      ENDDO
    ENDIF
  ENDDO
ENDIF
!
!
!
200 CONTINUE
WRITE(msg,*) 'NP:', NP
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
REWIND(30)
!
IF( NP.NE.0 ) THEN
  ALLOCATE(P(NP,4))
  IF( Naux>0 ) THEN
    ALLOCATE(AUXNAMES(Naux))
    ALLOCATE(AUX(SIZE(P,1),Naux))
    AUX(:,:) = 0.d0
    IF( fx>0 .AND. fy>0 .AND. fz>0 ) THEN
      AUXNAMES(fx) = "fx"
      AUXNAMES(fy) = "fy"
      AUXNAMES(fz) = "fz"
    ENDIF
    IF( fixx>0 .AND. fixy>0 .AND. fixz>0 ) THEN
      AUXNAMES(fixx) = "fixx"
      AUXNAMES(fixy) = "fixy"
      AUXNAMES(fixz) = "fixz"
    ENDIF
  ENDIF
ELSE
  nerr=nerr+1
  GOTO 1000
ENDIF
!
!
!Check whether conventional notation (a,b,c,alpha,beta,gamma) or the celldm(:) were defined.
!According to Quantum Espresso documentation, either one OR the other must be used, not both.
!Here, if both were defined then output a warning, use the celldm(:) and ignore the conventional notation.
IF( NINT(a).NE.0 .AND. NINT(b).NE.0 .AND. NINT(c).NE.0 ) THEN
  convnot = .TRUE.
ENDIF
IF( NINT(SUM(celldm(2:6))).NE.0 ) THEN
  celldm_defined = .TRUE.
ENDIF
IF( convnot .AND. celldm_defined ) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(1705,(/""/),(/0.d0/))
  convnot = .FALSE.
ENDIF
!
!
!Define the cell vectors according to the value of ibrav
!Note: if ibrav==0 then the cell vectors will be defined later,
!     in the section CELL_PARAMETERS
WRITE(msg,*) "ibrav =", ibrav
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
IF(ibrav.NE.0) THEN
  IF(convnot) THEN
    !vectors were defined in conventional notation a,b,c,alpha,beta,gamma
    ! => use them to define the values of the celldm(:)
    celldm(1) = a
    celldm(2) = b/a
    celldm(3) = c/a
    celldm(4) = cosgamma
    celldm(5) = cosbeta
    IF(ibrav==14) THEN
      celldm(4) = cosalpha
      celldm(5) = cosbeta
      celldm(6) = cosgamma
    ENDIF
  ENDIF
  !Define the cell vectors H(:,:) according to the value of ibrav and the celldm(:)
  CALL QE_PW_IBRAV(ibrav,celldm(:),H)
ENDIF
!
!
DO
  READ(30,'(a128)',END=1000,ERR=800) temp
  IF( temp(1:15)=='CELL_PARAMETERS' .OR. temp(1:15)=='cell_parameters' ) THEN
    !Note: as in Quantum Espresso, this section is ignored if ibrav is not zero
    IF(ibrav==0) THEN
      !Read cell parameters
      READ(30,*,END=800,ERR=800) H(1,1), H(1,2), H(1,3)
      READ(30,*,END=800,ERR=800) H(2,1), H(2,2), H(2,3)
      READ(30,*,END=800,ERR=800) H(3,1), H(3,2), H(3,3)
      msg = ADJUSTL(temp(16:))
      IF( LEN_TRIM(msg)<=0 .OR. INDEX(msg,'alat')>0 .OR. INDEX(msg,'ALAT')>0 ) THEN     
        !Cell dimensions are in alat units (this is the default if nothing is specified)
        !If celldm(1) was not specified, the H(:,:) are in Bohrs
        !Otherwise the H(:,:) must be normalized to alat=celldm(1)
        IF( celldm(1).NE.0.d0 ) THEN
          H(:,:) = celldm(1)*H(:,:)
        ENDIF
        cell_units='B' 
      ELSEIF( msg(1:4)=="bohr" .OR. msg(1:4)=="Bohr" ) THEN
        !Cell dimensions are already in Bohrs
        cell_units='B'
      ELSEIF( msg(1:3)=='ang' ) THEN
        !Cell dimensions are in Angströms
        cell_units='A'
      ENDIF
    ENDIF
    !
  ELSEIF( temp(1:16)=='ATOMIC_POSITIONS' .OR. temp(1:16)=='atomic_positions' ) THEN
    msg = ADJUSTL(temp(17:))
    DO i=1,SIZE(P,1)
      IF( fixx>0 .AND. fixy>0 .AND. fixz>0 ) THEN
        READ(30,*,END=800,ERR=800) species, P(i,1), P(i,2), P(i,3), AUX(i,fixx), AUX(i,fixy), AUX(i,fixz)
      ELSE
        READ(30,*,END=800,ERR=800) species, P(i,1), P(i,2), P(i,3)
      ENDIF
      CALL ATOMNUMBER(species,P(i,4))
    ENDDO
    !If coordinates were reduced or in lattice coordinates, convert them
    IF( LEN_TRIM(msg)<=0 .OR. INDEX(msg,'alat')>0 .OR. INDEX(msg,'ALAT')>0 ) THEN
      !Atom coordinates are in alat units (this is the default if nothing is specified)
      P(:,1:3) = celldm(1)*P(:,1:3)
      atpos_units=cell_units
    ELSEIF( msg(1:4)=="bohr" .OR. msg(1:4)=="Bohr" ) THEN
      !Atom positions are in Bohrs
      atpos_units='B'
    ELSEIF( msg(1:7)=='crystal' .OR. msg(1:7)=='CRYSTAL' ) THEN
      CALL FRAC2CART(P,H)
    ELSEIF( msg(1:3)=='ang' ) THEN
      !Atom positions are in Angströms
      atpos_units='A'
    ENDIF
    !
    !
  ELSEIF( temp(1:13)=='ATOMIC_FORCES' .OR. temp(1:16)=='atomic_forces' ) THEN
    DO i=1,SIZE(AUX,1)
      READ(30,*,END=800,ERR=800) species, AUX(i,fx), AUX(i,fy), AUX(i,fz)
    ENDDO
    !
  ENDIF
ENDDO
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
CLOSE(30)
!
IF( cell_units=="B" .AND. atpos_units=="A" ) THEN
  !Atom coordinates are in angstroms while cell vectors are in Bohrs
  ! => convert cell vectors into angströms for consistency
  nwarn=nwarn+1
  CALL ATOMSK_MSG(1706,(/msg/),(/0.d0/))
  H(:,:) = H(:,:) * 1.d10*a_bohr
ENDIF
!
!
!
END SUBROUTINE READ_QEPW
!
END MODULE in_qe_pw