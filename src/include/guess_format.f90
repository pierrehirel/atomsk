MODULE guess_form
!
!**********************************************************************************
!*  GUESS_FORM                                                                    *
!**********************************************************************************
!* This module tries to determine the format of a file.                           *
!* List of supported formats: see flist in "globalvar.f90".                       *
!* The purpose of this module is to recognize the format of a file                *
!* EVEN IF IT HAS NO EXTENSION, OR IF ITS EXTENSION IS ERRONEOUS.                 *
!* (e.g. a CFG file with extension ".xyz").                                       *
!* If the file doesn't exist or has to be written, then only its name or          *
!* extension is taken into account.                                               *
!* Otherwise if the file exists, then its content is ALWAYS parsed and            *
!* searched for specific keywords. Then:                                          *
!* - if the file has an extension then the format is guessed from both            *
!*   the extension *and* the file content;                                        *
!* - if the file has no extension then the format is guessed from file content.   *
!* The guessing is based on a "score": the higher the score, the more likely      *
!* the file format. Scores are saved in an integer array "fscore" that has        *
!* the same size as the array "flist" (see "globalvar.f90").                      *
!* The routine "SET_SCORE" (at the end of this module) is used to modify scores.  *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 22 Feb. 2024                                     *
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
CHARACTER:: c
CHARACTER(LEN=4):: fstatus   !file status; 'read' or 'writ'
CHARACTER(LEN=5):: extension, infileformat
CHARACTER(LEN=8):: nline     !line number
CHARACTER(LEN=128):: msg
CHARACTER(LEN=1024):: test, test2
CHARACTER(LEN=1024):: filename
CHARACTER(LEN=4096),INTENT(IN):: inputfile
LOGICAL:: fileexists
LOGICAL:: fileisopened
LOGICAL:: formatted  !is the file formatted?
INTEGER:: certainty
INTEGER:: likely
INTEGER:: strlength, i, j, k
INTEGER:: Ncol
INTEGER:: NP
INTEGER,DIMENSION(SIZE(flist,1)):: fscore
REAL(dp):: testreal
!
!
!Initialize variables
 certainty = 10
infileformat = "     "
extension = "     "
test=''
i = 1
testreal = 0.d0
fscore(:) = 0  !set all scores to zero
formatted = .TRUE.
!
msg = 'Entering GUESS_FORMAT...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
INQUIRE(FILE=inputfile,EXIST=fileexists)
IF(.NOT.fileexists .OR. fstatus=='writ') THEN
  !File does not exist => detection will be based only on file name or extension
  certainty = 2
ELSE
  !File exists => detection will be based on file name, extension, and content
  certainty = 8
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
  !increase corresponding score by a lot, and decrease other scores
  IF( TRIM(ADJUSTL(filename))=='system.in' .OR. &
    & TRIM(ADJUSTL(filename))=='system.out'     ) THEN
    !Special case for MOLDY files
    CALL SET_SCORE(fscore,extension,4,1)
  ELSE
    CALL SET_SCORE(fscore,extension,10,2)
  ENDIF
  !
ENDIF
!
IF( .NOT.ANY(fscore(:)>=10) ) THEN
  !Either there was no extension, or extension was not a recognized format
  !Let us try to guess format from file name
  !Special case for VASP POSCAR/CONTCAR/CHG files
  IF( filename=='POSCAR' .OR. filename=='CONTCAR' .OR. filename=='CHG' .OR. filename=='CHGCAR' ) THEN
    CALL SET_SCORE(fscore,"pos  ",10,1)
  ENDIF
  !Special case for VASP OUTCAR files
  IF( filename=='OUTCAR' ) THEN
    CALL SET_SCORE(fscore,"vout  ",10,1)
  ENDIF
  !Special case for VASP POSCAR/CONTCAR/CHG files
  IF( fstatus=="read" .AND. ( INDEX(filename,'POSCAR')>0 .OR. INDEX(filename,'CONTCAR')>0 .OR. &
    &  INDEX(filename,'CHG')>0 ) ) THEN
    CALL SET_SCORE(fscore,"pos  ",5,1)
  ENDIF
  !Special case for COORAT files
  IF( filename=='COORAT' ) THEN
    CALL SET_SCORE(fscore,"coo  ",10,1)
  ENDIF
  IF( fstatus=="read" .AND. INDEX(filename,'COORAT')>0 ) THEN
    CALL SET_SCORE(fscore,"coo  ",4,1)
  ENDIF
  !Special case for DL_POLY CONFIG files
  IF( filename=='CONFIG' .OR. filename=='REVCON' .OR. &
    & filename=='CFGMIN' .OR. filename=='HISTORY'     ) THEN
    CALL SET_SCORE(fscore,"dlp  ",10,1)
  ENDIF
  IF( fstatus=="read" .AND.                  &
      & ( INDEX(filename,'CONFIG')>0  .OR.   &
      &   INDEX(filename,'REVCON')>0  .OR.   &
      &   INDEX(filename,'CFGMIN')>0  .OR.   &
      &   INDEX(filename,'HISTORY')>0     )   ) THEN
    CALL SET_SCORE(fscore,"dlp  ",4,1)
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
  OPEN(UNIT=25,FILE=inputfile,FORM='FORMATTED',STATUS='UNKNOWN',ERR=800)
  REWIND(25)
ENDIF
msg = 'Parsing file content...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
IF(fileexists) THEN
!
  !Try to detect non-ASCII characters to check if file is binary
  j=0
  k=1  !counter for lines. Reading the first 5 lines should be sufficient
  DO WHILE( j==0 .AND. formatted .AND. k<=5 )
    READ(25,'(a1)',IOSTAT=j) c
    formatted = formatted .AND. ( IACHAR(c)<=127 )
  ENDDO
  IF( .NOT.formatted ) THEN
    !File is in binary format: not supported for now
    infileformat = "xxx"
    GOTO 1000
  ENDIF
  !
  !Go back to beginning of file
  REWIND(25)
  !
  !Try to confirm the format by reading the contents of the file
  !We limit it to the first 50 lines to avoid loosing too much time
  DO i=1,50
    READ(25,'(a1024)',END=300,ERR=300) test
    test = TRIM(ADJUSTL(test))
    strlength = LEN_TRIM(test)
    !
    IF(verbosity==4) THEN
      IF(i<50) THEN
        WRITE(nline,'(i8)') i
        msg = 'GUESS_FORMAT line '//TRIM(ADJUSTL(nline))//':  '//TRIM(test)
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ELSEIF(i==50) THEN
        msg = 'GUESS_FORMAT         ... discontinued ...'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDIF
    ENDIF
    !
    IF( test(1:2)=="#!" ) THEN
      !Shebang: determine what comes after
      IF( INDEX(test(3:),"sh")>0 ) THEN
        !It is a Linux shell script
      ELSEIF( INDEX(test(3:),"atomsk")>0 ) THEN
        !It is an Atomsk script
      ENDIF
    ELSEIF( test(1:1)=='#' ) THEN
      !It is just a comment, skip line
    !
    !Search for patterns corresponding to BOP format
    ELSEIF(test(1:4)=='ITER') THEN
      CALL SET_SCORE(fscore,"bop  ",2,1)
    ELSEIF(test(1:4)=='A   ') THEN
      CALL SET_SCORE(fscore,"bop  ",2,1)
    ELSEIF(test(1:4)=='LEN ') THEN
      CALL SET_SCORE(fscore,"bop  ",4,1)
    ELSEIF(test(1:6)=='LATPAR') THEN
      CALL SET_SCORE(fscore,"bop  ",4,1)
    ELSEIF(test(1:4)=='ND  ') THEN
      CALL SET_SCORE(fscore,"bop  ",4,1)
    ELSEIF(test(1:3)=='D  ') THEN
      CALL SET_SCORE(fscore,"bop  ",4,1)
    ELSEIF(test(1:6)=='NINERT') THEN
      CALL SET_SCORE(fscore,"bop  ",10,1)
    ELSEIF(test(1:6)=='DINERT') THEN
      CALL SET_SCORE(fscore,"bop  ",10,1)
    ELSEIF(test(1:3)=='UNRLD') THEN
      CALL SET_SCORE(fscore,"bop  ",10,1)
    !
    !Search for patterns corresponding to BOPfox format
    ELSEIF(test(1:6)=='aLat =') THEN
      CALL SET_SCORE(fscore,"bx   ",4,1)
    ELSEIF(test(1:4)=='a1 =') THEN
      CALL SET_SCORE(fscore,"bx   ",4,1)
    ELSEIF(test(1:4)=='a2 =') THEN
      CALL SET_SCORE(fscore,"bx   ",4,1)
    ELSEIF(test(1:4)=='a3 =') THEN
      CALL SET_SCORE(fscore,"bx   ",4,1)
    ELSEIF(test(1:7)=='coord =') THEN
      CALL SET_SCORE(fscore,"bx   ",4,1)
    ELSEIF(test(1:15)=='magnetisation =') THEN
      CALL SET_SCORE(fscore,"bx   ",6,2)
    !
    !Search for patterns corresponding to ABINIT format
    ELSEIF(test(1:5)=='acell') THEN
      CALL SET_SCORE(fscore,"abin ",6,1)
    ELSEIF(test(1:5)=='natom') THEN
      CALL SET_SCORE(fscore,"abin ",4,1)
    ELSEIF(test(1:6)=='ntypat') THEN
      CALL SET_SCORE(fscore,"abin ",10,1)
    ELSEIF(test(1:5)=='typat') THEN
      CALL SET_SCORE(fscore,"abin ",6,1)
    ELSEIF(test(1:5)=='rprim') THEN
      CALL SET_SCORE(fscore,"abin ",10,1)
    ELSEIF(test(1:4)=='xred') THEN
      CALL SET_SCORE(fscore,"abin ",10,1)
    ELSEIF(test(1:5)=='xcart') THEN
      CALL SET_SCORE(fscore,"abin ",10,1)
    ELSEIF(test(1:5)=='znucl') THEN
      CALL SET_SCORE(fscore,"abin ",10,1)
    !
    !Search for patterns corresponding to CEL format
    ELSEIF(TRIM(ADJUSTL(test))=='*') THEN ! EOF signal (short CEL file)
      CALL SET_SCORE(fscore,"cel  ",5,1)
    !
    !Search for patterns corresponding to CFG format
    ELSEIF(test(1:19)=='Number of particles') THEN
      IF( i<=2 ) THEN
        CALL SET_SCORE(fscore,"cfg  ",10,1)
      ELSE
        CALL SET_SCORE(fscore,"cfg  ",6,1)
      ENDIF
    ELSEIF(test(1:3)=='A =' .OR. test(1:2)=='A=') THEN
      CALL SET_SCORE(fscore,"cfg  ",4,1)
    ELSEIF(test(1:6)=='H(1,1)' .OR. test(1:6)=='H(1,2)' .OR. &
          &test(1:6)=='H(2,1)' .OR. test(1:6)=='H(2,2)' .OR. &
          &test(1:6)=='H(3,1)' .OR. test(1:6)=='H(3,2)' .OR. &
          &test(1:6)=='H(3,3)') THEN
      CALL SET_SCORE(fscore,"cfg  ",4,1)
    ELSEIF(test(1:13)=='.NO_VELOCITY.') THEN
      CALL SET_SCORE(fscore,"cfg  ",10,1)
    ELSEIF(test(1:11)=='entry_count') THEN
      CALL SET_SCORE(fscore,"cfg  ",10,1)
    ELSEIF(test(1:9)=='auxiliary') THEN
      CALL SET_SCORE(fscore,"cfg  ",4,1)
    !
    !Search for patterns corresponding to CIF format
    ELSEIF( test(1:7)=='_audit_' ) THEN
      CALL SET_SCORE(fscore,"cif  ",10,1)
    ELSEIF( test(1:6)=='_cell_' ) THEN
      CALL SET_SCORE(fscore,"cif  ",6,1)
    ELSEIF( test(1:5)=='loop_' ) THEN
      CALL SET_SCORE(fscore,"cif  ",10,1)
    ELSEIF( test(1:6)=='_atom_' ) THEN
      CALL SET_SCORE(fscore,"cif  ",6,1)
    ELSEIF( test(1:8)=='_diffrn_' ) THEN
      CALL SET_SCORE(fscore,"cif  ",10,1)
    ELSEIF( test(1:7)=='_exptl_' ) THEN
      CALL SET_SCORE(fscore,"cif  ",10,1)
    !
    !Search for patterns corresponding to CML format
    !ELSEIF(test(1,1)=='<') THEN
    !  iscml = iscml+0.1d0
    !
    !Search for patterns corresponding to COORAT format
    ELSEIF(test(1:6)=='natom=') THEN
      CALL SET_SCORE(fscore,"coo  ",6,1)
    !
    !Search for patterns corresponding to CRYSTAL format (*.d12)
    !NOTE: keywords "CRYSTAL" and "SLAB" can also correspond
    !      to XSF file format (see below)
    ELSEIF(test(1:8)=='ATOMSUBS') THEN
      CALL SET_SCORE(fscore,"d12  ",10,1)
    ELSEIF(test(1:8)=='MOLECULE') THEN
      CALL SET_SCORE(fscore,"d12  ",4,1)
    ELSEIF(test(1:7)=='POLYMER') THEN
      CALL SET_SCORE(fscore,"d12  ",4,1)
    ELSEIF(test(1:8)=='ATOMDISP') THEN
      CALL SET_SCORE(fscore,"d12  ",4,1)
    ELSEIF(test(1:3)=='END') THEN
      CALL SET_SCORE(fscore,"d12  ",1,1)
    ELSEIF(test(1:7)=='OPTGEOM') THEN
      CALL SET_SCORE(fscore,"d12  ",6,1)
    ELSEIF(test(1:8)=='FULLOPTG') THEN
      CALL SET_SCORE(fscore,"d12  ",10,1)
    ELSEIF(test(1:8)=='CELLONLY') THEN
      CALL SET_SCORE(fscore,"d12  ",4,1)
    ELSEIF(test(1:6)=='SHRINK') THEN
      CALL SET_SCORE(fscore,"d12  ",10,1)
    ELSEIF(test(1:7)=='FMIXING') THEN
      CALL SET_SCORE(fscore,"d12  ",10,1)
    !
    !Search for patterns corresponding to SIESTA FDF format
    ELSEIF(test(1:10)=='SystemName') THEN
      CALL SET_SCORE(fscore,"fdf  ",4,1)
    ELSEIF(test(1:11)=='SystemLabel') THEN
      CALL SET_SCORE(fscore,"fdf  ",6,1)
    ELSEIF(test(1:13)=='NumberOfAtoms') THEN
      CALL SET_SCORE(fscore,"fdf  ",6,1)
    ELSEIF(test(1:15)=='NumberOfSpecies') THEN
      CALL SET_SCORE(fscore,"fdf  ",6,1)
    ELSEIF(test(1:23)=='AtomicCoordinatesFormat') THEN
      CALL SET_SCORE(fscore,"fdf  ",10,1)
    ELSEIF(test(1:7)=='%block ') THEN
      CALL SET_SCORE(fscore,"fdf  ",10,1)
    ELSEIF(test(1:10)=='%endblock ') THEN
      CALL SET_SCORE(fscore,"fdf  ",10,1)
    !
    !Search for patterns corresponding to GIN format
    ELSEIF(test(1:5)=='title') THEN
      CALL SET_SCORE(fscore,"gin  ",4,1)
      CALL SET_SCORE(fscore,"str  ",4,1)
    ELSEIF(test(1:4)=='cell') THEN
      CALL SET_SCORE(fscore,"gin  ",4,1)
      CALL SET_SCORE(fscore,"str  ",4,1)
    ELSEIF(test(1:4)=='vect') THEN
      CALL SET_SCORE(fscore,"gin  ",2,1)
    ELSEIF(test(1:3)=='car') THEN
      CALL SET_SCORE(fscore,"gin  ",2,1)
    ELSEIF(test(1:4)=='buck') THEN
      CALL SET_SCORE(fscore,"gin  ",4,1)
    ELSEIF(test(1:7)=='species') THEN
      CALL SET_SCORE(fscore,"gin  ",2,1)
    ELSEIF(test(1:6)=='spring') THEN
      CALL SET_SCORE(fscore,"gin  ",4,1)
    ELSEIF(test(1:3)=='fra') THEN
      CALL SET_SCORE(fscore,"gin  ",2,1)
    ELSEIF(test(1:7)=='keyword') THEN
      CALL SET_SCORE(fscore,"gin  ",2,1)
    ELSEIF(test(1:4)=='cutp') THEN
      CALL SET_SCORE(fscore,"gin  ",2,1)
    ELSEIF(test(4:7)=='core') THEN
      CALL SET_SCORE(fscore,"gin  ",2,1)
    ELSEIF(test(4:7)=='shel') THEN
      CALL SET_SCORE(fscore,"gin  ",6,1)
    ELSEIF(test(1:5)=='total') THEN
      CALL SET_SCORE(fscore,"gin  ",2,1)
    ELSEIF(test(1:10)=='dump every') THEN
      CALL SET_SCORE(fscore,"gin  ",2,1)
    ELSEIF(test(1:4)=='output') THEN
      CALL SET_SCORE(fscore,"gin  ",2,1)
    !
    !Search for patterns corresponding to IMD format
    ELSEIF(test(1:3)=='#F ') THEN
      CALL SET_SCORE(fscore,"imd  ",4,1)
    ELSEIF(test(1:3)=='#C ') THEN
      CALL SET_SCORE(fscore,"imd  ",4,1)
    ELSEIF(test(1:3)=='#X ') THEN
      CALL SET_SCORE(fscore,"imd  ",4,1)
    ELSEIF(test(1:3)=='#Y ') THEN
      CALL SET_SCORE(fscore,"imd  ",4,1)
    ELSEIF(test(1:3)=='#Z ') THEN
      CALL SET_SCORE(fscore,"imd  ",4,1)
    ELSEIF(test(1:3)=='#E ') THEN
      CALL SET_SCORE(fscore,"imd  ",4,1)
    !
    !Search for patterns corresponding to JEMS format
    ELSEIF(test(1:5)=='file|') THEN
      CALL SET_SCORE(fscore,"jems ",10,1)
    ELSEIF(test(1:7)=='system|') THEN
      CALL SET_SCORE(fscore,"jems ",10,1)
    ELSEIF(test(1:9)=='HMSymbol|') THEN
      CALL SET_SCORE(fscore,"jems ",10,1)
    ELSEIF(test(1:4)=='rps|') THEN
      CALL SET_SCORE(fscore,"jems ",10,1)
    ELSEIF(test(1:8)=='lattice|') THEN
      CALL SET_SCORE(fscore,"jems ",10,1)
    ELSEIF(test(1:5)=='atom|') THEN
      CALL SET_SCORE(fscore,"jems ",10,1)
    ELSEIF(test(1:3)=='aff|') THEN
      CALL SET_SCORE(fscore,"jems ",10,1)
    !
    !Search for patterns corresponding to LAMMPS data format
    ELSEIF(strlength-4==INDEX(test,'atoms').AND.strlength>4) THEN
      CALL SET_SCORE(fscore,"lmp  ",3,1)
    ELSEIF(strlength-4==INDEX(test,'Atoms').AND.strlength>4) THEN
      CALL SET_SCORE(fscore,"lmp  ",4,1)
    ELSEIF(strlength-9==INDEX(test,'atom types').AND.strlength>9) THEN
      CALL SET_SCORE(fscore,"lmp  ",10,2)
    ELSEIF(strlength-6==INDEX(test,'xlo xhi').AND.strlength>6) THEN
      CALL SET_SCORE(fscore,"lmp  ",10,2)
    ELSEIF(strlength-6==INDEX(test,'ylo yhi').AND.strlength>6) THEN
      CALL SET_SCORE(fscore,"lmp  ",10,2)
    ELSEIF(strlength-6==INDEX(test,'zlo zhi').AND.strlength>6) THEN
      CALL SET_SCORE(fscore,"lmp  ",10,2)
    ELSEIF(strlength-7==INDEX(test,'xy xz yz').AND.strlength>7) THEN
      CALL SET_SCORE(fscore,"lmp  ",10,1)
    !
    !Search for patterns corresponding to LAMMPS custom format
    ELSEIF(test(1:5)=='ITEM:') THEN
      CALL SET_SCORE(fscore,"lmc  ",20,5)
    !
    !Search for patterns corresponding to Protein Data Bank (PDB) format
    ELSEIF(test(1:6)=='HEADER') THEN
      CALL SET_SCORE(fscore,"pdb  ",4,1)
    ELSEIF(test(1:6)=='TITLE ') THEN
      CALL SET_SCORE(fscore,"pdb  ",4,1)
    ELSEIF(test(1:6)=='AUTHOR') THEN
      CALL SET_SCORE(fscore,"pdb  ",4,1)
    ELSEIF(test(1:6)=='REMARK') THEN
      CALL SET_SCORE(fscore,"pdb  ",4,1)
    ELSEIF(test(1:6)=='ATOM  ') THEN
      CALL SET_SCORE(fscore,"pdb  ",4,1)
    !
    !Search for patterns corresponding to POSCAR format (VASP)
    ELSEIF(test(1:18)=='Selective dynamics') THEN
      CALL SET_SCORE(fscore,"pos  ",10,1)
    ELSEIF(test(1:7)=='Direct ') THEN
      CALL SET_SCORE(fscore,"pos  ",10,1)
    ELSEIF(test(1:9)=='Cartesian') THEN
      CALL SET_SCORE(fscore,"pos  ",10,1)
    !
    !Search for patterns corresponding to OUTCAR format (VASP)
    ELSEIF(i==1 .AND. test(1:4)=='vasp') THEN
      CALL SET_SCORE(fscore,"vout ",10,1)
    ELSEIF(test(1:6)=='INCAR:') THEN
      CALL SET_SCORE(fscore,"vout ",6,1)
    ELSEIF(test(1:7)=='POTCAR:') THEN
      CALL SET_SCORE(fscore,"vout ",6,1)
    ELSEIF(test(1:8)=='KPOINTS:') THEN
      CALL SET_SCORE(fscore,"vout ",6,1)
    ELSEIF(test(1:6)=='distr:') THEN
      CALL SET_SCORE(fscore,"vout ",10,1)
    ELSEIF(test(1:7)=='distrk:') THEN
      CALL SET_SCORE(fscore,"vout ",10,1)
    !
    !Search for patterns corresponding to Quantum Espresso PW format
    ELSEIF(test(1:8)=='&CONTROL' .OR. test(1:8)=='&control') THEN
      CALL SET_SCORE(fscore,"pw   ",10,1)
    ELSEIF(test(1:7)=='&SYSTEM' .OR. test(1:7)=='&system') THEN
      CALL SET_SCORE(fscore,"pw   ",10,1)
    ELSEIF(test(1:9)=='&ELECTRON' .OR. test(1:9)=='&electron') THEN
      CALL SET_SCORE(fscore,"pw   ",10,1)
    ELSEIF(test(1:5)=='&IONS' .OR. test(1:5)=='&ions') THEN
      CALL SET_SCORE(fscore,"pw   ",10,1)
    ELSEIF(test(1:5)=='&CELL' .OR. test(1:5)=='&cell') THEN
      CALL SET_SCORE(fscore,"pw   ",10,1)
    ELSEIF(test(1:14)=='ATOMIC_SPECIES' .OR. test(1:14)=='atomic_species') THEN
      CALL SET_SCORE(fscore,"pw   ",10,1)
    ELSEIF(test(1:16)=='ATOMIC_POSITIONS' .OR. test(1:16)=='atomic_positions') THEN
      CALL SET_SCORE(fscore,"pw   ",3,1)
      CALL SET_SCORE(fscore,"pwo  ",3,1)
    ELSEIF(test(1:8)=='K_POINTS' .OR. test(1:8)=='k_points') THEN
      CALL SET_SCORE(fscore,"pw   ",4,1)
    ELSEIF(test(1:15)=='CELL_PARAMETERS' .OR. test(1:15)=='cell_parameters') THEN
      CALL SET_SCORE(fscore,"pw   ",4,1)
      CALL SET_SCORE(fscore,"pwo  ",4,1)
    ELSEIF(test(1:11)=='CONSTRAINTS' .OR. test(1:11)=='constraints') THEN
      CALL SET_SCORE(fscore,"pw   ",6,1)
    ELSEIF(test(1:10)=='OCCUPATION' .OR. test(1:10)=='occupation') THEN
      CALL SET_SCORE(fscore,"pw   ",4,1)
    !
    !Search for patterns corresponding to Quantum Espresso PW output format
    ELSEIF(test(1:13)=='Program PWSCF') THEN
      CALL SET_SCORE(fscore,"pwo  ",10,1)
    !
    !Search for patterns corresponding to PDFFIT structure file format (*.str)
    ELSEIF(test(1:13)=='dcell') THEN
      CALL SET_SCORE(fscore,"str  ",6,1)
    ELSEIF(test(1:13)=='spcgr') THEN
      CALL SET_SCORE(fscore,"str  ",10,1)
    !
    !Search for patterns corresponding to VESTA format
    ELSEIF(test(1:13)=='#VESTA_FORMAT') THEN
      CALL SET_SCORE(fscore,"vesta",10,1)
    ELSEIF(test(1:5)=='TRANM') THEN
      CALL SET_SCORE(fscore,"vesta",10,1)
    ELSEIF(test(1:7)=='LTRANSL') THEN
      CALL SET_SCORE(fscore,"vesta",10,1)
    ELSEIF(test(1:7)=='LORIENT') THEN
      CALL SET_SCORE(fscore,"vesta",10,1)
    ELSEIF(test(1:7)=='LMATRIX') THEN
      CALL SET_SCORE(fscore,"vesta",10,1)
    ELSEIF(test(1:7)=='LCELLP') THEN
      CALL SET_SCORE(fscore,"vesta",10,1)
    !
    !Search for patterns corresponding to XMD format
    ELSEIF(test(1:9)=='POSITION ') THEN
      CALL SET_SCORE(fscore,"xmd  ",3,1)
    ELSEIF(test(1:7)=='POSVEL ') THEN
      CALL SET_SCORE(fscore,"xmd  ",10,1)
    !
    !Search for patterns corresponding to XSF format
    ELSEIF(test=='PRIMCOORD') THEN
      CALL SET_SCORE(fscore,"xsf  ",10,1)
    ELSEIF(test=='PRIMVEC') THEN
      CALL SET_SCORE(fscore,"xsf  ",10,1)
    ELSEIF(test=='CONVVEC') THEN
      CALL SET_SCORE(fscore,"xsf  ",10,1)
    ELSEIF(test=='CRYSTAL') THEN
      CALL SET_SCORE(fscore,"xsf  ",6,1)
      CALL SET_SCORE(fscore,"d12  ",6,1)
    ELSEIF(test=='ATOMS') THEN
      CALL SET_SCORE(fscore,"xsf  ",6,1)
    ELSEIF(test=='SLAB') THEN
      CALL SET_SCORE(fscore,"xsf  ",6,1)
      CALL SET_SCORE(fscore,"d12  ",6,1)
    !
    !Search for patterns corresponding to special XYZ format
    ELSEIF(test(1:4)=='alat') THEN
      CALL SET_SCORE(fscore,"sxyz ",6,1)
    ELSEIF(test(1:5)=='super') THEN
      CALL SET_SCORE(fscore,"sxyz ",6,1)
    !
    !None of the above worked => this format has no recognizable keyword
    ! => it can be XYZ, MOLDY, XV...
    ELSE
      IF( i==1 ) THEN
        !Try to read the string "Atomsk binary file"
        IF( INDEX(test,'Atomsk binary file')>0 ) THEN
          CALL SET_SCORE(fscore,"atsk ",10,2)
        ELSE
          !Try to read an integer
          !If successful, it can be XYZ or MOLDY format
          READ(test,*,ERR=220,END=220) NP
          CALL SET_SCORE(fscore,"xyz  ",2,1)
          CALL SET_SCORE(fscore,"mol  ",2,1)
        ENDIF
        !
      ELSEIF( i==2 ) THEN
        !Try to read the keywords "Lattice" and "Properties" from 2nd line
        !If successful, it is extended XYZ format
        IF( INDEX(test,'Lattice=').NE.0 ) CALL SET_SCORE(fscore,"xyz  ",6,1)
        IF( INDEX(test,'Properties=').NE.0 ) CALL SET_SCORE(fscore,"xyz  ",6,1)
        !Otherwise, try to read 3 integers
        !If successful, it could be a MOLDY file or a DL_POLY file
        READ(test,*,ERR=220,END=220) NP, NP, NP
        CALL SET_SCORE(fscore,"mol  ",6,1)
        CALL SET_SCORE(fscore,"dlp  ",6,1)
        !
      ENDIF
      !
      220 CONTINUE
      IF( i==1 .OR. i==2 .OR. i==3 ) THEN
        !Try to read 6 reals
        !If successful, it could be SIESTA XV format
        READ(test,*,ERR=230,END=230) testreal, testreal, testreal, &
                                   & testreal, testreal, testreal
        CALL SET_SCORE(fscore,"xv   ",6,0)
      ENDIF
      !
      230 CONTINUE
      IF( i==3 .OR. i==4 .OR. i==5 ) THEN
        !Try to read 3 reals
        !If successful for i=3,4,5 then it can be MOLDY or DL_POLY file
        READ(test,*,ERR=240,END=240) testreal, testreal, testreal
        CALL SET_SCORE(fscore,"dlp  ",2,0)
        CALL SET_SCORE(fscore,"mol  ",2,0)
      ENDIF
      !
      240 CONTINUE
      IF( i>2 .AND. test(1:1).NE.'#' ) THEN
        !Determine number of columns of data
        Ncol = 0
        test2 = ADJUSTL(test)
        DO WHILE(LEN_TRIM(test2)>0)
          READ(test2,*,ERR=241,END=241) msg
          PRINT*, TRIM(msg)
          j = LEN_TRIM(msg) + 1
          test2 = TRIM(ADJUSTL(test2(j:)))
          Ncol=Ncol+1
        ENDDO
        241 CONTINUE
        !
        IF( Ncol==3 ) THEN
          !Atomeye, DL_POLY CONFIG (or REVCON), and VASP POSCAR formats have lines "x y z"
          READ(test,*,ERR=260,END=260) testreal, testreal, testreal
          CALL SET_SCORE(fscore,"cfg  ",2,0)
          CALL SET_SCORE(fscore,"dlp  ",2,0)
          CALL SET_SCORE(fscore,"pos  ",2,0)
        ELSEIF( Ncol==4 ) THEN
          !XYZ format has lines  "N x y z"
          j=1
          READ(test,*,ERR=243,END=243) NP, testreal, testreal, testreal
          CALL SET_SCORE(fscore,"lmp  ",1,0)
          CALL SET_SCORE(fscore,"xyz  ",1,0)
          GOTO 260
          243 CONTINUE
          j=2
          !MOLDY format has lines  "x y z N"
          READ(test,*,ERR=244,END=244) testreal, testreal, testreal, NP
          CALL SET_SCORE(fscore,"mol  ",1,0)
          GOTO 260
          244 CONTINUE
        ELSEIF( Ncol==8 ) THEN
          !SIESTA XV format has lines  "N Z x y z fx fy fz"
          READ(test,*,ERR=260,END=260) NP, NP, testreal, testreal, testreal, &
                                      & testreal, testreal, testreal
          CALL SET_SCORE(fscore,"xv   ",1,0)
        ENDIF
        !
        260 CONTINUE
      ENDIF
    !
    ! -- Add other formats here --
    ENDIF
    !
    290 CONTINUE
  !
  ENDDO  !End loop on reading lines in the file
  !
  CLOSE(25)
  !
ENDIF   !If fileexists
!
!
!
300 CONTINUE
!Find the best score
likely = MAXVAL( fscore(:) )
!
IF( verbosity==4 ) THEN
  WRITE(msg,*) 'Scores of file formats: '
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(flist,1)
    WRITE(msg,*) "      ", flist(i,1), "   ", fscore(i)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
ENDIF
!
!
IF( likely<NINT(DBLE(certainty)/3.d0) ) THEN
  !If the best score is very low then it is probably an unknown format
  infileformat = 'xxx'
!
ELSE
  !Otherwise we have detected a known file format
  !Get position of highest score
  i = MAXLOC( fscore(:) , DIM=1 )
  !Save the name of matching file format
  infileformat = flist(i,1)
  !
  !If the score is lower than the threshold, display a warning
  IF( verbosity==4 .AND. likely<certainty ) THEN
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
INQUIRE(UNIT=25,OPENED=fileisopened)
IF(fileisopened) CLOSE(25)
!
!
END SUBROUTINE GUESS_FORMAT
!
!
!
!********************************************************
! SET_SCORE
! Increases the score of the given file format (fformat)
! by incr, and decreases the score of all other
! file formats by decr.
!********************************************************
SUBROUTINE SET_SCORE(fscore,fformat,incr,decr)
!
IMPLICIT NONE
CHARACTER(LEN=5),INTENT(IN):: fformat
INTEGER:: i
INTEGER,INTENT(IN):: incr !increment given file format by this much
INTEGER,INTENT(IN):: decr !decrement other file formats by this much
INTEGER,DIMENSION(:),INTENT(INOUT):: fscore
!
DO i=1,SIZE(flist,1)
  IF(flist(i,1)==fformat) THEN
    fscore(i) = fscore(i)+ABS(incr)
  ELSE
    IF(fscore(i)>0) fscore(i) = fscore(i)-ABS(decr)
  ENDIF
ENDDO
!
END SUBROUTINE SET_SCORE
!
!
!
END MODULE guess_form
