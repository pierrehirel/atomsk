MODULE writeout
!
!**********************************************************************************
!*  WRITEOUT                                                                      *
!**********************************************************************************
!* This module writes atomic data to one or several files.                        *
!*                                                                                *
!* WHAT MUST BE PROVIDED:                                                         *
!* prefix    string containing the name for output file(s). The extensions        *
!*           (e.g. ".cfg", ".xyz", etc.) will be appened to this prefix           *
!*           accordingly for each written file.                                   *
!* outfileformats                                                                 *
!*           string array of unknown dimension, containing the extensions         *
!*           (e.g. "cfg", "xyz", etc.) corresponding to file formats that         *
!*           will be written. If not allocated then the user will be prompted.    *
!* H         3 x 3 real array containing vectors of the supercell. H(1,:) is      *
!*           the first vector, H(2,:) the second, H(3,:) the third.               *
!* P         N x 4 real array containing atom positions. N=number of atoms,       *
!*           each line contains the coordinates (X,Y,Z) and atomic number.        *
!* S         if allocated, N x 4 real array containing positions of shells        *
!*           (for ionic core/shell model), and number of shells N must be         *
!*           identical to N in P. If S is not allocated, then no shell exist.     *
!* comment   an array of strings that can contain anything, or be unallocated.    *
!* AUX       if allocated, N x M real array containing values of auxiliary        *
!*           properties of atoms, e.g. their velocity, forces, etc. The first     *
!*           dimension N must be equal to the number of atoms N in P, the second  *
!*           dimension M must be equal to the number of auxiliary properties.     *
!*           If AUX is not allocated then no auxiliary property exists, and       *
!*           AUXNAMES must also be non-allocated.                                 *
!* AUXNAMES  if allocated, string array of dimension M, M being the number of     *
!*           existing auxiliary properties (must be equal to M in AUX). If        *
!*           AUXNAMES is not allocated, then AUX must also be non-allocated.      *
!*                                                                                *
!* WHAT IS RETURNED BY THIS ROUTINE:                                              *
!* no variable nor array, only files are written on the disk.                     *
!**********************************************************************************
!* (C) October 2010 - Pierre Hirel                                                *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 10 July 2019                                     *
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
!Load modules
USE comv
USE constants
USE functions
USE guess_form
USE messages
USE files
USE subroutines
!
!Modules managing output files
USE out_atsk
USE out_abinit
USE out_bop
USE out_bopfox
USE out_cfg
USE out_cel
USE out_cif
USE out_crystal
USE out_csv
USE out_dlp_cfg
USE out_gulp_gin
USE out_imd
USE out_jems
USE out_lammps_data
USE out_mbpp_coorat
USE out_moldy
USE out_pdb
USE out_qe_pw
USE out_siesta_fdf
USE out_siesta_xv
USE out_str
USE out_vasp_poscar
USE out_vesta
USE out_xmd
USE out_xsf
USE out_xyz
! -- please add other formats in alphabetical order --
!
!
!
CONTAINS
!
SUBROUTINE WRITE_AFF(prefix,outfileformats,H,P,S,comment,AUXNAMES,AUX)
!
!Declare variables
IMPLICIT NONE
CHARACTER(LEN=5):: newformat
CHARACTER(LEN=5):: zone
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to write
CHARACTER(LEN=8):: date
CHARACTER(LEN=10):: time
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128):: username
CHARACTER(LEN=*):: prefix
CHARACTER(LEN=4096):: outputfile
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN),OPTIONAL:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,OPTIONAL:: comment
LOGICAL:: fileexists
INTEGER:: i, j, k
INTEGER:: q, qs, itypes  !position of charges, shell charges, atom types in AUX
INTEGER,DIMENSION(8):: values
REAL(dp):: Q_total  !total electric charge of the cell
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P          !positions of atoms (or ionic cores)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN),OPTIONAL:: S !positions of ionic shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: atypes    !atom types and their number
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!
username=""
q = 0
qs = 0
Q_total = 0.d0
!
IF(verbosity==4) THEN
  msg = 'Entering WRITE_AFF'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  msg = 'prefix: '//TRIM(prefix)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) 'Size P: ', SIZE(P,1), SIZE(P,2)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  IF(ALLOCATED(S)) THEN
    WRITE(msg,*) 'Size S: ', SIZE(S,1), SIZE(S,2)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
  IF(ALLOCATED(AUXNAMES)) THEN
    WRITE(msg,*) 'Size AUXNAMES, AUX: ', SIZE(AUXNAMES), SIZE(AUX,1)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
  IF(ALLOCATED(outfileformats)) THEN
    WRITE(msg,*) 'Activated output file formats: '
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    msg = '  '
    DO i=1, SIZE(outfileformats)
      msg = TRIM(msg)//' '//TRIM(outfileformats(i))
    ENDDO
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
ENDIF
!
!
100 CONTINUE
!If user only specified "NULL" instead of an output file name, don't write any file
IF( prefix=="NULL" .AND. (.NOT.ALLOCATED(outfileformats) .OR. SIZE(outfileformats)<=1) ) THEN
  CALL ATOMSK_MSG(3007,(/TRIM("")/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Check that array P is allocated, if not output error message
IF( .NOT.ALLOCATED(P) .OR. SIZE(P,1)==0 ) THEN
  nerr = nerr+1
  CALL ATOMSK_MSG(3800,(/TRIM(msg)/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Parse array P and check for NaN
!                   THIS IS VERY TIME CONSUMING!!!
!Therefore this is done only for relatively small systems.
!When building bigger systems, it is assumed that the user has tested
!his scripts on smaller systems already, and knows what he is doing
IF( SIZE(P,1)<30000 ) THEN
  WRITE(msg,*) 'Looking for NaN in array P...'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  CALL CHECKNAN(P,i)
  IF( i.NE.0 ) THEN
    nerr=nerr+1
    CALL ATOMSK_MSG(3803,(/""/),(/DBLE(i)/))
    GOTO 1000
  ENDIF
  !Also check auxiliary properties
  IF( ALLOCATED(AUX) ) THEN
    WRITE(msg,*) 'Looking for NaN in array AUX...'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    CALL CHECKNAN(AUX,i)
    IF( i.NE.0 ) THEN
      nerr=nerr+1
      CALL ATOMSK_MSG(3804,(/""/),(/DBLE(i)/))
      GOTO 1000
    ENDIF
  ENDIF
ENDIF
!
!Check sizes of the arrays
CALL CHECK_ARRAY_CONSISTENCY(P,S,AUX,AUXNAMES,i)
IF( i.NE.0 ) THEN
  IF( i==1 ) THEN
    msg = 'S'
  ELSEIF( i==2 ) THEN
    msg = 'AUX'
  ELSEIF( i==3 ) THEN
    msg = 'AUXNAMES'
  ENDIF
  nerr=nerr+1
  CALL ATOMSK_MSG(1802,(/msg/),(/0.d0/))
ENDIF
!
!
!If no comment is defined, define a generic one
!Note: Intel Fortran Compiler seems to be unable to interpret the following "IF" when
!     all conditions are in the same parenthesis, so use a trick by setting a value for i
i=0
IF( .NOT.ALLOCATED(comment) .OR. SIZE(comment)<=0 ) THEN
  i=1
ELSE
  IF( LEN_TRIM(comment(1))==0 .OR. comment(1)(1:17)=="# File generated " ) THEN
    i=1
  ENDIF
ENDIF
IF(i==1) THEN
  WRITE(msg,*) 'Generating comment...'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !Get user name: this is environment-dependent
#if defined(WINDOWS)
  CALL GET_ENVIRONMENT_VARIABLE('USERNAME',username)
#else
  CALL GET_ENVIRONMENT_VARIABLE('USER',username)
#endif
  !
  IF( ALLOCATED(comment) .AND. SIZE(comment)<=0 ) DEALLOCATE(comment)
  IF(.NOT.ALLOCATED(comment)) ALLOCATE(comment(1))
  CALL DATE_AND_TIME(DATE, TIME, ZONE, VALUES)
  !Generate compound formula
  !CALL COMPFORMULA(P,AUXNAMES,AUX,formula,smass_tot)
  !Generate message for the comment
  CALL CREATE_DATE(VALUES,username,comment(1))
ENDIF
!
!Make sure that each comment line starts with a hash sign (#)
IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
  DO i=1,SIZE(comment)
    msg = TRIM(ADJUSTL(comment(i)))
    IF(msg(1:1).NE.'#') comment(i) = '# '//TRIM(comment(i))
  ENDDO
ENDIF
!
!
IF(LEN_TRIM(prefix)<=0) THEN
  CALL ATOMSK_MSG(3700,(/TRIM(msg)/),(/0.d0/))
  READ(*,*) prefix
ENDIF
!
!Check if the file name has a recognizable extension
IF( SCAN(prefix,'.',BACK=.TRUE.)>0 ) THEN
  CALL GUESS_FORMAT(prefix,newformat,'writ')
  IF(newformat.NE.'xxx') THEN
    !A file format was recognized => add it to the list of outfileformats(:)
    CALL SET_OUTPUT(outfileformats,newformat,.TRUE.)
  ENDIF
ENDIF
!
!
150 CONTINUE
!If no output format is active then ask the user
IF( .NOT.ALLOCATED(outfileformats) .OR. SIZE(outfileformats)==0 ) THEN
  nwarn = nwarn+1
  CALL ATOMSK_MSG(3701,(/TRIM(prefix)/),(/0.d0/))
  READ(*,*) newformat
  !check that this format can be recognized
  IF(.NOT.ANY(listofformats==newformat)) THEN
    GOTO 150
  ELSE
    CALL SET_OUTPUT(outfileformats,newformat,.TRUE.)
  ENDIF
ENDIF
!
!Check that the array outfileformats does not contain twice the same output format
!If so, just replace duplicates by blanks
IF(SIZE(outfileformats)>1) THEN
  outfileformats(:) = ADJUSTL(outfileformats(:))
  DO i=1,SIZE(outfileformats)-1
    DO j=i+1,SIZE(outfileformats)
      IF(outfileformats(j) == outfileformats(i) ) THEN
        outfileformats(j) = ''
      ELSEIF( outfileformats(i)=="sxyz" ) THEN
        !Special case for xyz format: only one *.xyz file can be output, give extended XYZ (exyz) highest priority
        IF( outfileformats(j)=="xyz" .OR. outfileformats(j)=="sxyz" ) THEN
          outfileformats(j) = ''
        ELSEIF( outfileformats(j)=="exyz" ) THEN
          outfileformats(i) = ''
        ENDIF
      ENDIF
    ENDDO
  ENDDO
ENDIF
!
!Check for the total electric charge and atom types
IF( ALLOCATED(AUXNAMES) .AND. SIZE(AUXNAMES)>0 ) THEN
  DO i=1,SIZE(AUXNAMES)
    msg = TRIM(ADJUSTL(AUXNAMES(i)))
    IF( msg(1:2)=='q ' ) THEN
      q = i
    ELSEIF( msg(1:3)=='qs ' ) THEN
      qs = i
    ELSEIF( msg(1:5)=='type ' ) THEN
      itypes = i
    ENDIF
  ENDDO
  !
  IF( .NOT.ALLOCATED(S) .OR. SIZE(S,1).NE.SIZE(AUX,1) ) THEN
    qs = 0
  ENDIF
  !
  Q_total = 0.d0
  IF( q>0 ) THEN
    DO i=1,SIZE(AUX,1)
      Q_total = Q_total + AUX(i,q)
      IF( qs>0 ) THEN
        Q_total = Q_total + AUX(i,qs)
      ENDIF
    ENDDO
  ENDIF
  !
  IF( DABS(Q_total) > 1.d-6 ) THEN
    !Total electric charge is non-zero => display a warning
    nwarn = nwarn+1
    CALL ATOMSK_MSG(3717,(/""/),(/Q_total/))
  ENDIF
  !
  !If output is activated for file formats that depend on atom "type",
  !then check against atoms that have different species but same type
  !(again, only for relatively small systems)
  IF( SIZE(P,1)<30000 ) THEN
    IF( itypes>0 .AND. itypes<=SIZE(AUX,2) ) THEN
      fileexists = .FALSE.
      !Check if user wants to write to a file format that doesn't support shells
      DO i=1,SIZE(outfileformats)
        IF( LEN_TRIM(outfileformats(i)) > 0 ) THEN
          SELECT CASE( outfileformats(i) )
          CASE('abinit','in','imd','IMD','fdf','FDF','siesta','lmp','LMP','lammps','LAMMPS','xmd','XMD')
            !Those file formats require atom "type"
            fileexists = .TRUE.
          CASE DEFAULT
            !Other file formats don't
          END SELECT
        ENDIF
      ENDDO
      IF(fileexists) THEN
        WRITE(msg,*) 'Checking against shared atom types...'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        k = 0
        DO i=1,SIZE(P,1)-1
          DO j=SIZE(P,1),i+1,-1
            IF( NINT(P(i,4)).NE.NINT(P(j,4)) .AND. NINT(AUX(i,itypes))==NINT(AUX(j,itypes)) ) THEN
              k = k+1
              EXIT
            ENDIF
          ENDDO
        ENDDO
        !
        IF( k>0 ) THEN
          !The same "type" was given to atoms of different species
          !Display a warning
          nwarn = nwarn+1
          CALL ATOMSK_MSG(3718,(/""/),(/0.d0/))
        ENDIF
      ENDIF
    ENDIF
  ENDIF
  !
ENDIF
!
!
!
200 CONTINUE
!!Then, write output file(s)
! by calling the module corresponding to the output format
! Output to several formats is possible
IF( ALLOCATED(S) .AND. SIZE(S,1)>0 ) THEN
  !Count actual number of shells
  j = 0
  DO i=1,SIZE(S,1)
    IF( NINT(S(i,4))>0 ) THEN
      j = j+1
    ENDIF
  ENDDO
  CALL ATOMSK_MSG(3000,(/""/),(/DBLE(SIZE(P,1)),DBLE(j)/))
ELSE
  CALL ATOMSK_MSG(3000,(/""/),(/DBLE(SIZE(P,1)),0.d0/))
ENDIF
!
!If shells exist (ionic core-shell model), then we must make sure
!that all the current output file formats support it. Otherwise, the output file
!will contain only the positions of the cores, but the shells will be lost.
!Atomsk will go with it, but a warning will be displayed
IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
  fileexists = .FALSE.
  !Check if user wants to write to a file format that doesn't support shells
  DO i=1,SIZE(outfileformats)
    IF( LEN_TRIM(outfileformats(i)) > 0 ) THEN
      SELECT CASE( outfileformats(i) )
      CASE('atsk','ATSK','csv','CSV','d12','dlp','DLP','gin','GIN','lmp','LMP')
        !Those file formats do support shells
      CASE DEFAULT
        !Other file formats don't
        fileexists = .TRUE.
      END SELECT
    ENDIF
  ENDDO
  !
  IF( fileexists ) THEN
    !User wants to write shell positions, but some output file formats don't support it
    !Display a warning
    nwarn = nwarn+1
    CALL ATOMSK_MSG(3716,(/""/),(/0.d0/))
  ENDIF
ENDIF
!
!If AUX contains data about partial occupancies, then we must make sure
!that all the current output file formats support it. Otherwise, the output file
!will contain atoms at the exact same position, and no information about their occupancy,
!which may not be desirable. Atomsk will go with it, but a warning will be displayed
IF( ALLOCATED(AUXNAMES) .AND. SIZE(AUXNAMES)>0 ) THEN
  j = 0
  DO i=1,SIZE(AUXNAMES)
    IF( AUXNAMES(i) == "occ" ) THEN
      !Occupancies are present
      j = i
    ENDIF
  ENDDO
  !
  IF( j>0 ) THEN
    fileexists = .FALSE.
    !Check if user wants to write to a file format that doesn't support partial occupancies
    DO i=1,SIZE(outfileformats)
      IF( LEN_TRIM(outfileformats(i)) > 0 ) THEN
        SELECT CASE( outfileformats(i) )
        CASE('atsk','ATSK','cel','CEL','cfg','CFG','cif','CIF','csv','CSV','d12','gin','GIN','jems','pdb','str','vesta')
          !Those file formats do support partial occupancies
        CASE DEFAULT
          !Other file formats don't
          fileexists = .TRUE.
        END SELECT
      ENDIF
    ENDDO
    !
    IF( fileexists ) THEN
      !User wants to write info about partial occupancies, but some output file formats don't support it
      !This is not a problem if all occupancies equal 1, but it is a problem otherwise
      !=> Check if some values of partial occupancies are different from 1
      fileexists = .FALSE.
      DO i=1,SIZE(AUX,1)
        IF( DABS( AUX(i,j) - 1.d0 ) > 1.d-9 ) THEN
          !This value is different from 1
          fileexists = .TRUE.
        ENDIF
      ENDDO
      !
      IF( fileexists ) THEN
        !Some occupancies are different from 1 => potential problem, display a warning
        nwarn = nwarn+1
        CALL ATOMSK_MSG(3715,(/""/),(/0.d0/))
      ENDIF
      !
    ENDIF
    !
  ENDIF
ENDIF
!
IF(SIZE(P(:,1))>10000000) THEN
  CALL ATOMSK_MSG(3,(/''/),(/0.d0/))
ENDIF
!
!
DO i=1,SIZE(outfileformats)
  SELECT CASE(outfileformats(i))
  !
  CASE('atsk','ATSK')
    CALL NAME_OUTFILE(prefix,outputfile,'atsk ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_ATSK(H,P,comment,AUXNAMES,AUX,outputfile,S)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('abin')
    CALL NAME_OUTFILE(prefix,outputfile,'in   ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_ABINIT(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('bop','BOP')
    CALL NAME_OUTFILE(prefix,outputfile,'bop  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_BOP(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('bx','BX')
    CALL NAME_OUTFILE(prefix,outputfile,'bx   ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_BOPFOX(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('cfg','CFG')
    CALL NAME_OUTFILE(prefix,outputfile,'cfg  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_CFG(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('cel','CEL')
    CALL NAME_OUTFILE(prefix,outputfile,'cel  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_CEL(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('cif','CIF')
    CALL NAME_OUTFILE(prefix,outputfile,'cif  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_CIF(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('coo','COO')
    outputfile = 'COORAT'
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_COORAT(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('csv','CSV')
    CALL NAME_OUTFILE(prefix,outputfile,'csv  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_CSV(H,P,S,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('d12')
    CALL NAME_OUTFILE(prefix,outputfile,'d12  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_CRYSTAL(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('dd','DD')
    !cannot write ddplot format here, only in mode "--ddplot"
    nwarn = nwarn+1
    CALL ATOMSK_MSG(3702,(/TRIM(outputfile)/),(/0.d0/))
  !
  CASE('dlp','DLP')
    outputfile = 'CONFIG'
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_DLP_CFG(H,P,comment,AUXNAMES,AUX,outputfile,S)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('fdf','FDF')
    CALL NAME_OUTFILE(prefix,outputfile,'fdf  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_FDF(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('gin','GIN')
    CALL NAME_OUTFILE(prefix,outputfile,'gin  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_GIN(H,P,comment,AUXNAMES,AUX,outputfile,S)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('imd','IMD')
    CALL NAME_OUTFILE(prefix,outputfile,'imd  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_IMD(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('jems','JEMS')
    CALL NAME_OUTFILE(prefix,outputfile,'txt  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_JEMS(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('lmp','LMP')
    CALL NAME_OUTFILE(prefix,outputfile,'lmp  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_LMP_DATA(H,P,S,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('mol','MOL')
    CALL NAME_OUTFILE(prefix,outputfile,'mol  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_MOLDY(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('pdb','PDB')
    CALL NAME_OUTFILE(prefix,outputfile,'pdb  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_PDB(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('pos','POS')
    outputfile = 'POSCAR'
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_POSCAR(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('pw','PW')
    CALL NAME_OUTFILE(prefix,outputfile,'pw   ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_QEPW(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('str','STR','stru','STRU')
    CALL NAME_OUTFILE(prefix,outputfile,'str  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_STR(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('vesta','VESTA')
    CALL NAME_OUTFILE(prefix,outputfile,'vesta')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_VESTA(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('xmd','XMD')
    CALL NAME_OUTFILE(prefix,outputfile,'xmd  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_XMD(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('xsf','XSF')
    CALL NAME_OUTFILE(prefix,outputfile,'xsf  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_XSF(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('xv','XV')
    CALL NAME_OUTFILE(prefix,outputfile,'XV   ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_XV(H,P,comment,AUXNAMES,AUX,outputfile)
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  CASE('xyz','XYZ','exyz','EXYZ','sxyz','SXYZ')
    CALL NAME_OUTFILE(prefix,outputfile,'xyz  ')
    INQUIRE(FILE=outputfile,EXIST=fileexists)
    IF( (fileexists .AND. .NOT.ignore) .OR. .NOT.fileexists) THEN
      IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
      CALL WRITE_XYZ(H,P,comment,AUXNAMES,AUX,outputfile,outfileformats(i))
    ELSE
      CALL ATOMSK_MSG(3001,(/TRIM(outputfile)/),(/0.d0/))
    ENDIF
  !
  ! -- please add other formats in alphabetical order --
  !
  CASE('')
    !empty entry: just ignore it
    CONTINUE
  CASE DEFAULT
    !all other cases: unknown format
    CALL ATOMSK_MSG(3710,(/TRIM(outfileformats(i))/),(/0.d0/))
  !
  END SELECT
  !
  IF(nerr>0) GOTO 800
  !
ENDDO
!
!
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(3801,(/TRIM(outputfile)/),(/0.d0/))
!
!
!
1000 CONTINUE
!
END SUBROUTINE WRITE_AFF
!
!
END MODULE
