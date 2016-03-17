MODULE guess_form
!
!**********************************************************************************
!*  GUESS_FORM                                                                    *
!**********************************************************************************
!* This module tries to determine the format of a file.                           *
!* If the file has an extension then the file format will be determined           *
!* from the extension *and* from the file content.                                *
!* If the file has no extension then the file format will be determined           *
!* only from the file content.                                                    *
!* The purpose of this module is to recognize the format of a file                *
!* EVEN IF IT HAS NO EXTENSION, OR IF ITS EXTENSION IS ERRONEOUS.                 *
!* (e.g. a CFG file with extension ".xyz").                                       *
!* For that purpose, if the file exists then its content is ALWAYS parsed         *
!* and searched for specific keywords.                                            *
!* Those keywords are used to recognize the file format.                          *
!* If the file doesn't exist or has to be written, then only the                  *
!* extension is taken into account.                                               *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 17 March 2016                                    *
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
!
USE comv
USE messages
USE files
USE subroutines
!
!
CONTAINS

SUBROUTINE GUESS_FORMAT(inputfile,infileformat,fstatus)
!
IMPLICIT NONE
CHARACTER(LEN=4):: fstatus   !file status; 'read' or 'writ'
CHARACTER(LEN=5):: extension, infileformat
CHARACTER(LEN=128):: msg
CHARACTER(LEN=1024):: test
CHARACTER(LEN=1024):: filename
CHARACTER(LEN=4096),INTENT(IN):: inputfile
LOGICAL:: fileexists
LOGICAL:: fileisopened
INTEGER:: strlength, i
INTEGER:: NP
REAL(dp):: isatsk
REAL(dp):: isbop, iscfg, iscel, iscif, iscml, iscoorat, isdd, isdlp, isgin, isimd
REAL(dp):: isjems, islmp, islmpc, ismoldy, ispdb, isposcar, isqepw, isqeout, isvesta
REAL(dp):: isxsf, isxv,isxmd, isxyz, isexyz, issxyz
REAL(dp):: likely
REAL(dp):: testreal
REAL(dp):: certainty
!
!
!Initialize variables
 certainty = 0.9d0
infileformat = ''
extension = ''
test=''
i = 1
testreal = 0.d0
isbop = 0.d0     !Bond-Order Potential format
iscel = 0.d0     !CEL Super-cell file for Dr. Probe
iscfg = 0.d0     !Atomeye CFG format
iscif = 0.d0     !Crystallographic Information File
iscml = 0.d0     !Chemical Markup Language
iscoorat = 0.d0  !Mixed-Basis PseudoPotential format
isdd = 0.d0      !ddplot format
isdlp = 0.d0     !DL_POLY format
isgin = 0.d0     !GULP input file format
isimd = 0.d0     !IMD format
isjems = 0.d0    !JEMS format
islmp = 0.d0     !LAMMPS data file
islmpc = 0.d0    !LAMMPS custom dump file
ismoldy=0.d0     !MOLDY file
isatsk=0.d0      !Binary atomsk format
ispdb=0.d0       !Protein Data Bank format
isposcar = 0.d0  !VASP POSCAR format
isqepw = 0.d0    !Quantum Espresso PWscf format
isqeout = 0.d0   !Quantum Espresso PWscf output format
isvesta = 0.d0   !VESTA format
isxmd = 0.d0     !XMD format
isxsf = 0.d0     !xCrySDen format
isxv = 0.d0      !SIESTA XV format
isxyz = 0.d0     !XYZ format
isexyz = 0.d0    !XYZ format (extended)
issxyz = 0.d0    !XYZ format (special)
! -- initialize other formats here --
!
msg = 'Entering GUESS_FORMAT...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
INQUIRE(FILE=inputfile,EXIST=fileexists)
IF(.NOT.fileexists .OR. fstatus=='writ') THEN
  !File does not exist => detection will be based only on file name or extension
  certainty = 0.2d0
ELSE
  !File exists => detection will be based on file name, extension, and content
  certainty = 0.8d0
ENDIF
!
!
100 CONTINUE
!Get the input file name, remove path if necessary
!(e.g. if input file is /path/to/file.ext then only keep file.ext)
strlength=SCAN(inputfile,pathsep,BACK=.TRUE.)
IF(strlength.NE.0) THEN
  filename = inputfile(strlength+1:)
ELSE
  filename = inputfile
ENDIF
!
msg = 'Check extension of the file: '//TRIM(filename)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
strlength = SCAN(filename,'.',BACK=.TRUE.)
IF( strlength > 0 ) THEN
  extension = filename(strlength+1:)
  msg = 'Extension: '//extension
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !The extension is considered as a strong indication of the file format:
  !increase corresponding counter by a lot
  SELECT CASE(extension)
  CASE('atsk','ATSK')
    isatsk = isatsk+0.6d0
  CASE('cfg','CFG')
    iscfg = iscfg+0.6d0
    iscoorat = iscoorat-0.4d0
  CASE('cel','CEL')
    iscel = iscel+0.6d0
  CASE('cif','CIF')
    iscif = iscif+0.6d0
  CASE('cml','CML')
    iscml = iscml+0.6d0
  CASE('DD','dd')
    isdd = isdd+0.6d0
  CASE('gin','GIN','res','RES','grs','GRS')
    isgin = isgin+0.6d0
  CASE('imd','IMD')
    isimd = isimd+0.6d0
  CASE('jems','JEMS')
    isjems = isjems+0.6d0
  CASE('lmc','LMC')
    islmpc = islmpc+0.6d0
  CASE('lmp','LMP')
    islmp = islmp+0.6d0
  CASE('mol','MOL')
    ismoldy = ismoldy+0.15d0
  CASE('pdb','PDB')
    ispdb = ispdb+0.6d0
  CASE('pw','PW')
    isqepw = isqepw+0.15d0
  CASE('vesta')
    isvesta = isvesta+0.6d0
  CASE('xmd','XMD')
    isxmd = isxmd+0.6d0
  CASE('xsf','XSF')
    isxsf = isxsf+0.6d0
    iscoorat = iscoorat-0.4d0
  CASE('XV','xv')
    isxv = isxv+0.6d0
  CASE('xyz','XYZ')
    isxyz = isxyz+0.6d0
    iscoorat = iscoorat-0.4d0
    iscfg = iscfg-0.4d0
  CASE('exyz','EXYZ')
    isexyz = isexyz+0.6d0
    iscoorat = iscoorat-0.4d0
    iscfg = iscfg-0.4d0
  CASE('sxyz','SXYZ')
    issxyz = issxyz+0.6d0
    iscoorat = iscoorat-0.4d0
    iscfg = iscfg-0.4d0
  CASE DEFAULT
    CONTINUE
  END SELECT
  !Special case for MOLDY files
  IF( TRIM(ADJUSTL(filename))=='system.in' .OR. &
    & TRIM(ADJUSTL(filename))=='system.out'     ) THEN
    ismoldy = ismoldy+0.15d0
  ENDIF
  !
ELSE
  !Special case for POSCAR/CONTCAR files
  IF( filename=='POSCAR' .OR. filename=='CONTCAR' ) THEN
    isposcar = certainty+0.4d0
  ENDIF
  !Special case for COORAT files
  IF( filename=='COORAT' ) THEN
    iscoorat = certainty+0.4d0
  ENDIF
  IF( fstatus=="read" .AND. INDEX(filename,'COORAT').NE.0 ) THEN
    iscoorat = iscoorat+0.2d0
  ENDIF
  !Special case for DL_POLY CONFIG files
  IF( filename=='CONFIG' .OR. filename=='REVCON' .OR. &
    & filename=='CFGMIN' .OR. filename=='HISTORY'     ) THEN
    isdlp = certainty+0.4d0
  ENDIF
  IF( fstatus=="read" .AND.                    &
      & INDEX(filename,'CONFIG').NE.0    .OR.  &
      & INDEX(filename,'REVCON').NE.0    .OR.  &
      & INDEX(filename,'CFGMIN').NE.0    .OR.  &
      & INDEX(filename,'HISTORY').NE.0         ) THEN
    isdlp = isdlp+0.2d0
  ENDIF
ENDIF
!
!If we want to write the file we don't care about its content
IF(.NOT.fileexists .OR. fstatus=='writ') THEN
  GOTO 300
ENDIF
!
!
200 CONTINUE
IF(fstatus.NE.'writ') THEN
  OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='UNKNOWN',ERR=800)
  REWIND(30)
ENDIF
msg = 'Parsing file content...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
IF(fileexists) THEN
  !Try to confirm the format by reading the contents of the file
  !We limit it to the first 50 lines to avoid loosing too much time
  DO i=1,50
    READ(30,'(a1024)',END=300,ERR=300) test
    test = TRIM(ADJUSTL(test))
    strlength = LEN_TRIM(test)
    !
    IF(verbosity==4) THEN
      IF(i<50) THEN
        msg = 'GUESS_FORMAT line:  '//TRIM(test)
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ELSEIF(i==50) THEN
        msg = 'GUESS_FORMAT line:  ... discontinued ...'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDIF
    ENDIF
    !
    !Search for patterns corresponding to BOP format
    IF(test(1:4)=='ITER') THEN
      isbop = isbop+0.2d0
    ELSEIF(test(1:4)=='A   ') THEN
      isbop = isbop+0.1d0
    ELSEIF(test(1:4)=='LEN ') THEN
      isbop = isbop+0.1d0
    ELSEIF(test(1:6)=='LATPAR') THEN
      isbop = isbop+0.2d0
    ELSEIF(test(1:4)=='ND  ') THEN
      isbop = isbop+0.2d0
    ELSEIF(test(1:3)=='D  ') THEN
      isbop = isbop+0.2d0
    ELSEIF(test(1:6)=='NINERT') THEN
      isbop = isbop+0.4d0
    ELSEIF(test(1:6)=='DINERT') THEN
      isbop = isbop+0.4d0
    ELSEIF(test(1:3)=='UNRLD') THEN
      isbop = isbop+0.4d0
    !
    !Search for patterns corresponding to CEL format
    ELSEIF(TRIM(ADJUSTL(test))=='*') THEN
      iscel = iscel+0.5d0 ! EOF signal (short CEL file)
    !
    !Search for patterns corresponding to CFG format
    ELSEIF(test(1:19)=='Number of particles') THEN
      IF( i<=2 ) THEN
        iscfg = iscfg+0.6d0
      ENDIF
    ELSEIF(test(1:3)=='A =' .OR. test(1:2)=='A=') THEN
      iscfg = iscfg+0.2d0
    ELSEIF(test(1:6)=='H(1,1)' .OR. test(1:6)=='H(1,2)' .OR. &
          &test(1:6)=='H(2,1)' .OR. test(1:6)=='H(2,2)' .OR. &
          &test(1:6)=='H(3,1)' .OR. test(1:6)=='H(3,2)' .OR. &
          &test(1:6)=='H(3,3)') THEN
      iscfg = iscfg+0.1d0
    ELSEIF(test(1:13)=='.NO_VELOCITY.') THEN
      iscfg = iscfg+0.6d0
    ELSEIF(test(1:11)=='entry_count') THEN
      iscfg = iscfg+0.6d0
    ELSEIF(test(1:9)=='auxiliary') THEN
      iscfg = iscfg+0.2d0
    !
    !Search for patterns corresponding to CIF format
    ELSEIF( test(1:7)=='_audit_' ) THEN
      iscif = iscif+0.3d0
    ELSEIF( test(1:6)=='_cell_' ) THEN
      iscif = iscif+0.3d0
    ELSEIF( test(1:5)=='loop_' ) THEN
      iscif = iscif+0.3d0
    ELSEIF( test(1:6)=='_atom_' ) THEN
      iscif = iscif+0.3d0
    ELSEIF( test(1:8)=='_diffrn_' ) THEN
      iscif = iscif+0.3d0
    ELSEIF( test(1:7)=='_exptl_' ) THEN
      iscif = iscif+0.3d0
    !
    !Search for patterns corresponding to CML format
    !ELSEIF(test(1,1)=='<') THEN
    !  iscml = iscml+0.1d0
    !
    !Search for patterns corresponding to COORAT format
    ELSEIF(test(1:6)=='natom=') THEN
      iscoorat = iscoorat+0.4d0
    !
    !Search for patterns corresponding to GIN format
    ELSEIF(test(1:5)=='title') THEN
      isgin = isgin+0.2d0
    ELSEIF(test(1:4)=='cell') THEN
      isgin = isgin+0.2d0
    ELSEIF(test(1:4)=='vect') THEN
      isgin = isgin+0.1d0
    ELSEIF(test(1:3)=='car') THEN
      isgin = isgin+0.1d0
    ELSEIF(test(1:4)=='buck') THEN
      isgin = isgin+0.3d0
    ELSEIF(test(1:7)=='species') THEN
      isgin = isgin+0.2d0
    ELSEIF(test(1:6)=='spring') THEN
      isgin = isgin+0.2d0
    ELSEIF(test(1:3)=='fra') THEN
      isgin = isgin+0.1d0
    ELSEIF(test(1:7)=='keyword') THEN
      isgin = isgin+0.3d0
    ELSEIF(test(1:4)=='cutp') THEN
      isgin = isgin+0.2d0
    ELSEIF(test(4:7)=='core') THEN
      isgin = isgin+0.3d0
    ELSEIF(test(4:7)=='shel') THEN
      isgin = isgin+0.3d0
    ELSEIF(test(1:5)=='total') THEN
      isgin = isgin+0.1d0
    ELSEIF(test(1:10)=='dump every') THEN
      isgin = isgin+0.2d0
    ELSEIF(test(1:4)=='output') THEN
      isgin = isgin+0.1d0
    !
    !Search for patterns corresponding to IMD format
    ELSEIF(test(1:3)=='#F ') THEN
      isimd = isimd+0.4d0
    ELSEIF(test(1:3)=='#C ') THEN
      isimd = isimd+0.4d0
    ELSEIF(test(1:3)=='#X ') THEN
      isimd = isimd+0.4d0
    ELSEIF(test(1:3)=='#Y ') THEN
      isimd = isimd+0.4d0
    ELSEIF(test(1:3)=='#Z ') THEN
      isimd = isimd+0.4d0
    ELSEIF(test(1:3)=='#E ') THEN
      isimd = isimd+0.4d0
    !
    !Search for patterns corresponding to JEMS format
    ELSEIF(test(1:5)=='file|') THEN
      isjems = isjems+0.5d0
    ELSEIF(test(1:7)=='system|') THEN
      isjems = isjems+0.5d0
    ELSEIF(test(1:9)=='HMSymbol|') THEN
      isjems = isjems+0.5d0
    ELSEIF(test(1:4)=='rps|') THEN
      isjems = isjems+0.6d0
    ELSEIF(test(1:8)=='lattice|') THEN
      isjems = isjems+0.3d0
    ELSEIF(test(1:5)=='atom|') THEN
      isjems = isjems+0.3d0
    ELSEIF(test(1:3)=='aff|') THEN
      isjems = isjems+0.3d0
    !
    !Search for patterns corresponding to LAMMPS data format
    ELSEIF(strlength-4==INDEX(test,'atoms').AND.strlength>4) THEN
      islmp = islmp+0.2d0
    ELSEIF(strlength-4==INDEX(test,'Atoms').AND.strlength>4) THEN
      islmp = islmp+0.1d0 ! Note: 'Atoms' and 'atoms' have different weights?
    ELSEIF(strlength-9==INDEX(test,'atom types').AND.strlength>9) THEN
      islmp = islmp+0.4d0
    ELSEIF(strlength-6==INDEX(test,'xlo xhi').AND.strlength>6) THEN
      islmp = islmp+0.4d0
    ELSEIF(strlength-6==INDEX(test,'ylo yhi').AND.strlength>6) THEN
      islmp = islmp+0.4d0
    ELSEIF(strlength-6==INDEX(test,'zlo zhi').AND.strlength>6) THEN
      islmp = islmp+0.4d0
    ELSEIF(strlength-7==INDEX(test,'xy xz yz').AND.strlength>7) THEN
      islmp = islmp+0.4d0
    !
    !Search for patterns corresponding to LAMMPS custom format
    ELSEIF(test(1:5)=='ITEM:') THEN
      islmpc = islmpc+0.4d0
      isxv = isxv-0.4d0
    !
    !Search for patterns corresponding to Protein Data Bank (PDB) format
    ELSEIF(test(1:6)=='HEADER') THEN
      ispdb = ispdb+0.1d0
    ELSEIF(test(1:6)=='TITLE') THEN
      ispdb = ispdb+0.1d0
    ELSEIF(test(1:6)=='AUTHOR') THEN
      ispdb = ispdb+0.1d0
    ELSEIF(test(1:6)=='REMARK') THEN
      ispdb = ispdb+0.1d0
    ELSEIF(test(1:6)=='ATOM  ') THEN
      ispdb = ispdb+0.1d0
    !
    !Search for patterns corresponding to POSCAR format (VASP)
    ELSEIF(test(1:18)=='Selective dynamics') THEN
      isposcar = isposcar+0.6d0
    ELSEIF(test(1:7)=='Direct ') THEN
      isposcar = isposcar+0.4d0
    ELSEIF(test(1:9)=='Cartesian') THEN
      isposcar = isposcar+0.3d0
      isdlp = isdlp-0.3d0
    !
    !Search for patterns corresponding to Quantum Espresso PW format
    ELSEIF(test(1:8)=='&CONTROL' .OR. test(1:8)=='&control') THEN
      isqepw = isqepw+0.4d0
    ELSEIF(test(1:7)=='&SYSTEM' .OR. test(1:7)=='&system') THEN
      isqepw = isqepw+0.4d0
    ELSEIF(test(1:9)=='&ELECTRON' .OR. test(1:9)=='&electron') THEN
      isqepw = isqepw+0.4d0
    ELSEIF(test(1:5)=='&IONS' .OR. test(1:5)=='&ions') THEN
      isqepw = isqepw+0.4d0
    ELSEIF(test(1:5)=='&CELL' .OR. test(1:5)=='&cell') THEN
      isqepw = isqepw+0.4d0
    ELSEIF(test(1:14)=='ATOMIC_SPECIES' .OR. test(1:14)=='atomic_species') THEN
      isqepw = isqepw+0.3d0
    ELSEIF(test(1:16)=='ATOMIC_POSITIONS' .OR. test(1:16)=='atomic_positions') THEN
      isqepw = isqepw+0.2d0
      isqeout = isqeout+0.2d0
    ELSEIF(test(1:8)=='K_POINTS' .OR. test(1:8)=='k_points') THEN
      isqepw = isqepw+0.1d0
    ELSEIF(test(1:15)=='CELL_PARAMETERS' .OR. test(1:15)=='cell_parameters') THEN
      isqepw = isqepw+0.2d0
      isqeout = isqeout+0.2d0
    ELSEIF(test(1:11)=='CONSTRAINTS' .OR. test(1:11)=='constraints') THEN
      isqepw = isqepw+0.2d0
    ELSEIF(test(1:10)=='OCCUPATION' .OR. test(1:10)=='occupation') THEN
      isqepw = isqepw+0.2d0
    !
    !Search for patterns corresponding to Quantum Espresso PW output format
    ELSEIF(test(1:13)=='Program PWSCF') THEN
      isqeout = isqeout+1.d0
    !
    !Search for patterns corresponding to VESTA format
    ELSEIF(test(1:13)=='#VESTA_FORMAT') THEN
      isvesta = isvesta+0.6d0
    ELSEIF(test(1:5)=='TRANM') THEN
      isvesta = isvesta+0.1d0
    ELSEIF(test(1:7)=='LTRANSL') THEN
      isvesta = isvesta+0.1d0
    ELSEIF(test(1:7)=='LORIENT') THEN
      isvesta = isvesta+0.1d0
    ELSEIF(test(1:7)=='LMATRIX') THEN
      isvesta = isvesta+0.1d0
    ELSEIF(test(1:7)=='LCELLP') THEN
      isvesta = isvesta+0.1d0
    !
    !Search for patterns corresponding to XMD format
    ELSEIF(test(1:9)=='POSITION ') THEN
      isxmd = isxmd+0.2d0
    ELSEIF(test(1:7)=='POSVEL ') THEN
      isxmd = isxmd+0.6d0
    !
    !Search for patterns corresponding to XSF format
    ELSEIF(test=='PRIMCOORD') THEN
      isxsf = isxsf+0.4d0
    ELSEIF(test=='PRIMVEC') THEN
      isxsf = isxsf+0.4d0
    ELSEIF(test=='CONVVEC') THEN
      isxsf = isxsf+0.4d0
    ELSEIF(test=='CRYSTAL') THEN
      isxsf = isxsf+0.4d0
    ELSEIF(test=='ATOMS') THEN
      isxsf = isxsf+0.4d0
    ELSEIF(test=='SLAB') THEN
      isxsf = isxsf+0.4d0
    !
    !Search for patterns corresponding to special XYZ format
    ELSEIF(test(1:4)=='alat') THEN
      isxyz = isxyz+0.2d0
      isdlp = isdlp-0.3d0
      ismoldy = ismoldy-0.3d0
    ELSEIF(test(1:5)=='super') THEN
      isxyz = isxyz+0.2d0
      isdlp = isdlp-0.3d0
      ismoldy = ismoldy-0.3d0
    !
    !None of the above worked => this format has no recognizable keyword
    ! => it can be XYZ, MOLDY, XV...
    ELSE
      IF( i==1 ) THEN
        !Try to read the string "Atomsk binary file"
        IF( INDEX(test,'Atomsk binary file')>0 ) THEN
          isatsk = isatsk+1.d0
        ELSE
          !Try to read an integer
          !If successful, it can be XYZ or MOLDY format
          READ(test,*,ERR=220,END=220) NP
          isxyz = isxyz+0.3d0
          ismoldy = ismoldy+0.3d0
          isdlp = isdlp-0.3d0
          isposcar = isposcar-0.3d0
        ENDIF
        !
      ELSEIF( i==2 ) THEN
        !Try to read the keywords "Lattice" and "Properties" from 2nd line
        !If successful, it is extended XYZ format
        IF( INDEX(test,'Lattice=').NE.0 ) isxyz = isxyz+0.6d0
        IF( INDEX(test,'Properties=').NE.0 ) isxyz = isxyz+0.6d0
        !Otherwise, try to read 3 integers
        !If successful, it could be a MOLDY file or a DL_POLY file
        READ(test,*,ERR=220,END=220) NP, NP, NP
        ismoldy = ismoldy + 0.4d0
        isdlp = isdlp + 0.4d0
        !
      ENDIF
      !
      220 CONTINUE
      IF( i==1 .OR. i==2 .OR. i==3 ) THEN
        !Try to read 6 reals
        !If successful, it could be SIESTA XV format
        READ(test,*,ERR=230,END=230) testreal, testreal, testreal, &
                                   & testreal, testreal, testreal
        isxv = isxv+0.2d0
      ENDIF
      !
      230 CONTINUE
      IF( i==3 .OR. i==4 .OR. i==5 ) THEN
        !Try to read 3 reals
        !If successful for i=3,4,5 then it can be MOLDY or DL_POLY file
        READ(test,*,ERR=240,END=240) testreal, testreal, testreal
        IF(i==3) THEN
          ismoldy = ismoldy+0.1d0
          isdlp = isdlp+0.1d0
        ELSE
          ismoldy = ismoldy+0.15d0
          isdlp = isdlp+0.15d0
        ENDIF
      ENDIF
      !
      240 CONTINUE
      IF( i>2 .AND. test(1:1).NE.'#' ) THEN
        !XYZ format has lines  "N x y z"
        READ(test,*,ERR=241,END=241) NP, testreal, testreal, testreal
        isxyz = isxyz+0.01d0
        isdlp = isdlp-0.1d0
        ismoldy = ismoldy-0.1d0
        isposcar = isposcar-0.1d0
        GOTO 249
        241 CONTINUE
        !MOLDY format has lines  "x y z N"
        READ(test,*,ERR=242,END=242) testreal, testreal, testreal, NP
        ismoldy = ismoldy+0.01d0
        isdlp = isdlp-0.1d0
        isposcar = isposcar-0.1d0
        GOTO 249
        242 CONTINUE
        !SIESTA XV format has lines  "N Z x y z fx fy fz"
        READ(test,*,ERR=243,END=243) NP, NP, testreal, testreal, testreal, &
                                    & testreal, testreal, testreal
        isxv = isxv+0.05d0
        isdlp = isdlp-0.1d0
        ismoldy = ismoldy-0.1d0
        isxyz = isxyz-0.1d0
        isposcar = isposcar-0.1d0
        243 CONTINUE
        GOTO 249
        !Atomeye, DL_POLY CONFIG (or REVCON), and VASP POSCAR formats have lines "x y z"
        READ(test,*,ERR=249,END=249) testreal, testreal, testreal
        isposcar = isposcar+0.01d0
        iscfg = iscfg+0.01d0
        isxyz = isxyz-0.1d0
        isdlp = isdlp+0.01d0
        ismoldy = ismoldy-0.1d0
        isxv = isxv-0.1d0
        249 CONTINUE
      ENDIF
    !
    ! -- Add other formats here --
    ENDIF
    !
    290 CONTINUE
  !
  ENDDO  !End loop on reading lines in the file
  !
  CLOSE(30)
  !
ENDIF   !If fileexists
!
!
!
300 CONTINUE
!Find the best score
likely = MAX(isbop,iscfg,iscel,iscif,iscml,iscoorat,isdd,isdlp,isgin,isimd,isjems,islmp, &
       &     islmpc,ismoldy,isatsk,ispdb,isposcar,isqepw,isqeout,isvesta,isxmd,isxsf, &
       &     isxv,isxyz,isexyz,issxyz)
!
IF( verbosity==4 ) THEN
  WRITE(msg,*) 'Scores of file formats: '
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   ATSK ', isatsk
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   BOP ', isbop
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   CFG ', iscfg
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   CEL ', iscel
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   CIF ', iscif
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   CML ', iscml
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   COORAT ', iscoorat
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   DDPLOT ', isdd
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   DLPOLY ', isdlp
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   GIN ', isgin
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   IMD ', isimd
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   JEMS', isjems
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   LMP ', islmp
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   LMC ', islmpc
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   MOLDY ', ismoldy
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   PDB ', ispdb
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   POSCAR ', isposcar
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   QEPW ', isqepw
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   QEOUT ', isqeout
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   VESTA ', isvesta
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   XMD ', isxmd
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   XSF ', isxsf
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   XV ', isxv
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   XYZ ', isxyz
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   EXYZ ', isexyz
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) '   SXYZ ', issxyz
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) 'MAX SCORE: ', likely
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!
!
IF(likely<certainty/3.d0) THEN
  !If the best score is very low then it is probably an unknown format
  infileformat = 'xxx'
!
ELSE
  !Otherwise we have detected a known file format
  IF(isatsk==likely) THEN
    infileformat = 'atsk'
  ELSEIF(isbop==likely) THEN
    infileformat = 'bop'
  ELSEIF(iscfg==likely) THEN
    infileformat = 'cfg'
  ELSEIF(iscel==likely) THEN
    infileformat = 'cel'
  ELSEIF(iscif==likely) THEN
    infileformat = 'cif'
  ELSEIF(iscoorat==likely) THEN
    infileformat = 'coo'
  ELSEIF(iscml==likely) THEN
    infileformat = 'cml'
  ELSEIF(isdd==likely) THEN
    infileformat = 'dd'
  ELSEIF(isdlp==likely) THEN
    infileformat = 'dlp'
  ELSEIF(isgin==likely) THEN
    infileformat = 'gin'
  ELSEIF(isimd==likely) THEN
    infileformat = 'imd'
  ELSEIF(isjems==likely) THEN
    infileformat = 'jems'
  ELSEIF(islmpc==likely) THEN
    infileformat = 'lmc'
  ELSEIF(islmp==likely) THEN
    infileformat = 'lmp'
  ELSEIF(ismoldy==likely) THEN
    infileformat = 'mol'
  ELSEIF(ispdb==likely) THEN
    infileformat = 'pdb'
  ELSEIF(isposcar==likely) THEN
    infileformat = 'pos'
  ELSEIF(isqepw==likely) THEN
    infileformat = 'pw'
  ELSEIF(isqeout==likely) THEN
    infileformat = 'pwo'
  ELSEIF(isvesta==likely) THEN
    infileformat = 'vesta'
  ELSEIF(isxmd==likely) THEN
    infileformat = 'xmd'
  ELSEIF(isxyz==likely) THEN
    infileformat = 'xyz'
  ELSEIF(isexyz==likely) THEN
    infileformat = 'exyz'
  ELSEIF(issxyz==likely) THEN
    infileformat = 'sxyz'
  ELSEIF(isxv==likely) THEN
    infileformat = 'xv'
  ELSEIF(isxsf==likely) THEN
    infileformat = 'xsf'
  ! -- add other formats here --
  !
  ELSE
    infileformat = 'xxx'
  ENDIF
  !
  !If the score is lower than the threshold, display a warning
  IF(verbosity==4 .AND. likely<certainty) THEN
    nwarn = nwarn+1
    CALL ATOMSK_MSG(1701,(/TRIM(infileformat)/),(/0.d0/))
  ENDIF
!
ENDIF
!
!
!
500 CONTINUE
msg = 'Detected file format: '//TRIM(infileformat)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
GOTO 1000
!
!
!
800 CONTINUE
nerr=nerr+1
CALL ATOMSK_MSG(1801,(/TRIM(inputfile)/),(/0.d0/))
!
!
!
1000 CONTINUE
!Make sure the file is closed
INQUIRE(UNIT=30,OPENED=fileisopened)
IF(fileisopened) CLOSE(30)
!
!
END SUBROUTINE GUESS_FORMAT
!
!
!
END MODULE guess_form
