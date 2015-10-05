MODULE modes
!
!
!**********************************************************************************
!*  MODES                                                                         *
!**********************************************************************************
!* This module deals with the different modes that atomsk can run.                *
!* Modes can be divided in two categories:                                        *
!* (1) Modes that only deal with file formats conversion or operations on         *
!*     atom coordinates; that includes the normal mode, list mode, ddplot mode,   *
!*     mode unwrap, etc.                                                          *
!* (2) Modes that analyze the system or compute some properties; that includes    *
!*     the mode "electric dipole moments", etc.                                   *
!* The modes of type (2) appear after the label 5000. They will be available      *
!* in the final executable only if the code is compiled with the flag -DMODES.    *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 22 Sept. 2015                                    *
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
!Modules for input and options
USE readin
USE options
USE deterH
!Modules for the different modes
USE mode_normal
USE mode_list
USE out_dd
USE mode_merge
USE mode_create
USE mode_unwrap
USE mode_average
!Modules for mode 1-in-all
USE oia_dlp_history
USE oia_qeout
USE oia_xsf
USE oia_xyz
!Module for mode all-in-one
USE aio
!Particular modules that compute something
USE mode_centrosym
USE mode_density
USE mode_difference
USE edm
USE mode_epola
USE mode_interpolate
USE mode_nye
USE mode_polycrystal
USE mode_rdf
!
!
CONTAINS
!
!
SUBROUTINE RUN_MODE(mode,options_array,outfileformats,pfiles,mode_param)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: mode     !mode in which the program runs
CHARACTER(LEN=*),DIMENSION(5),INTENT(IN):: pfiles !pfiles(1)=file1
                                                  !pfiles(2)=file2
                                                  !pfiles(3)=filefirst
                                                  !pfiles(4)=filesecond
                                                  !pfiles(5)=listfile
CHARACTER(LEN=1):: axis !an axis: x, y or z
CHARACTER(LEN=1):: answer !"y" or "n"
CHARACTER(LEN=2):: Pspecies  !species forming polyhedra (mode EDM)
CHARACTER(LEN=2):: species
CHARACTER(LEN=5):: outfileformat
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to output
CHARACTER(LEN=10):: create_struc
CHARACTER(LEN=128):: msg, test, temp
CHARACTER(LEN=128):: property !name of the property whose density will be calculated (mode DENSITY)
CHARACTER(LEN=4096):: file1, file2  !should be inputfile, ouputfile (or vice-versa)
CHARACTER(LEN=4096):: listfile      !used by modes 'filelist' and 'rdf'
CHARACTER(LEN=4096):: filefirst, filesecond !used by modes 'ddplot' and 'merge'
CHARACTER(LEN=4096):: inputfile, outputfile
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: commentfirst, commentsecond !used by modes 'ddplot' and 'merge'
CHARACTER(LEN=4096),DIMENSION(:),ALLOCATABLE:: merge_files  !files for merge mode
CHARACTER(LEN=2),DIMENSION(20):: create_species
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: mode_param  !parameters for some special modes
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: i, ioptions, j, k, m
INTEGER:: Nimages   !number of images (mode --interpolate)
INTEGER:: strlength
INTEGER,DIMENSION(2):: NT_mn
REAL(dp):: NNN, rdf_maxR, rdf_dr, smass
REAL(dp),DIMENSION(3):: create_a0
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystallographic orientation of the system (mode create)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, Pfirst, Psecond, S
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties of atoms
!
!
!Initialize variables
inputfile=''
outputfile=''
axis = ''
NT_mn(:) = 0
file1 = pfiles(1)
file2 = pfiles(2)
filefirst = pfiles(3)
filesecond = pfiles(4)
listfile = pfiles(5)
 create_species(:) = ''
ORIENT(:,:) = 0.d0
!Setup the files correctly for some modes
IF(mode=='list') THEN
  file1=listfile
ENDIF
IF(ALLOCATED(S)) DEALLOCATE(S)
!
!
IF(verbosity==4) THEN
  msg = 'file1, file2:'//TRIM(file1)//', '//TRIM(file2)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  IF(ALLOCATED(mode_param)) THEN
    msg = 'MODE_PARAM ARRAY:'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    i=1
    DO WHILE( LEN_TRIM(mode_param(i)).NE.0 .AND. i<=SIZE(mode_param) )
      msg = mode_param(i)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      i=i+1
    ENDDO
  ENDIF
  !
  IF(ALLOCATED(options_array)) THEN
    msg = 'OPTIONS ARRAY:'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ioptions = 0
    DO i=1,SIZE(options_array)
      ioptions = ioptions+1
      msg = TRIM(options_array(ioptions))
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
  ENDIF
  !
  msg = 'running mode: '//TRIM(mode)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!
!
!
100 CONTINUE
SELECT CASE(mode)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  MODE NORMAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('normal')
  msg = 'Mode normal'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !Normal mode: we just convert a single file to another format.
  !
  !If inputfile is not known already, and if file1 and/or file2
  ! were specified as command-line parameters, then figure out which
  ! one is the input file and which one is the output file
  IF(inputfile=='') THEN
    file1 = TRIM(ADJUSTL(file1))
    file2 = TRIM(ADJUSTL(file2))
    msg = 'file1, file2: '//TRIM(file1)//', '//TRIM(file2)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !
    CALL FIND_INOUT(file1,file2)
    !
    inputfile = TRIM(ADJUSTL(file1))
    outputfile = TRIM(ADJUSTL(file2))
  ENDIF
  !
  msg = 'inputfile, outputfile: '//TRIM(inputfile)//', '//TRIM(outputfile)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !If we have an output file name, activate its extension
  IF(LEN_TRIM(outputfile).NE.0) THEN
    CALL GUESS_FORMAT(outputfile,outfileformat,'writ')
    IF(outfileformat.NE.'xxx') CALL SET_OUTPUT(outfileformats,outfileformat,.TRUE.)
  ELSE
    j=SCAN(inputfile,pathsep,BACK=.TRUE.)
    strlength = SCAN(inputfile,'.',BACK=.TRUE.)
    IF(strlength>j) THEN
      outputfile = inputfile(1:strlength-1)
    ELSE
      outputfile = inputfile
    ENDIF
  ENDIF
  msg = 'outputfile: '//TRIM(outputfile)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !We can now convert the inputfile into outputfile
  !The "outfileformat" parameter here is dummy, since the variables
  !for output (output_xyz, output_cfg, etc.) are passed to the
  !subroutine CONVERT_AFF
  CALL CONVERT_AFF(inputfile,options_array,outputfile,outfileformats)
  msg = 'LEAVING MODE_NORMAL...'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
200 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  MODE LIST
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('list')
  msg = 'Mode list'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  listfile = file1
  !List mode:
  !in this mode the program has to read a file containing
  !a list of files to convert.
  CALL ATOMSK_MSG(4000,(/TRIM(listfile)/),(/0.d0/))
  !
  CALL LIST_XYZ(listfile,options_array,outputfile,outfileformats)
!
!
300 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  MODE DDPLOT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('ddplot')
  CALL SET_OUTPUT(outfileformats,'all  ',.FALSE.)
  msg = 'Mode ddplot'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !DDplot mode:
  !in this mode the program has to read 2 files (no matter their format)
  !and create a file for ddplot.
  CALL ATOMSK_MSG(4010,(/''/),(/0.d0/))
  !Check if the files exit
  CALL CHECKFILE(filefirst,'read')
  CALL CHECKFILE(filesecond,'read')
  outputfile = TRIM(ADJUSTL(file1))
  IF(outputfile=='') CALL NAME_OUTFILE(filesecond,outputfile,'dd   ')
  IF(.NOT. overw) CALL CHECKFILE(outputfile,'writ')
  overw = .TRUE.
  !
  !Read the two input files and apply options if any
  CALL READ_AFF(filefirst,H,Pfirst,S,commentfirst,AUXNAMES,AUX)
  IF(nerr>=1) GOTO 10000
  CALL OPTIONS_AFF(options_array,H,Pfirst,S,AUXNAMES,AUX,ORIENT,SELECT)
  CALL READ_AFF(filesecond,H,Psecond,S,commentsecond,AUXNAMES,AUX)
  IF(nerr>=1) GOTO 10000
  CALL OPTIONS_AFF(options_array,H,Psecond,S,AUXNAMES,AUX,ORIENT,SELECT)
  !
  !!Determine if the cell vectors were found or not
  IF( VECLENGTH(H(1,:))==0.d0 .OR. VECLENGTH(H(2,:))==0.d0 &
    & .OR. VECLENGTH(H(3,:))==0.d0 ) THEN
    !!If basis vectors H were not found in the input file, then compute them
    ! (this does not work perfectly yet, see determine_H.f90)
    CALL ATOMSK_MSG(4701,(/''/),(/0.d0/))
    !
    CALL DETERMINE_H(H,Pfirst)
    !
  ENDIF
  !
  !Use the two XYZ files to produce the DD file
  CALL WRITE_DD(H,Pfirst,Psecond,outputfile)
!
!
!
400 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  MODE MERGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('merge')
  !Merge mode:
  !in this mode the program reads m files (no matter their format)
  !and merges them into one single file
  !
  msg = 'MERGE MODE:'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Read array mode_param(:)
  i=1
  READ(mode_param(i),*,ERR=7000,END=7000) temp
  IF( temp=='x' .OR. temp=='y' .OR. temp=='z' .OR. &
    & temp=='X' .OR. temp=='Y' .OR. temp=='Z'      ) THEN
    READ(temp,*,ERR=7000,END=7000) axis
    i=i+1
    READ(mode_param(i),*,ERR=7000,END=7000) m
  ELSE
    READ(temp,*,ERR=7000,END=7000) m
  ENDIF
  !
  IF(ALLOCATED(merge_files)) DEALLOCATE(merge_files)
  ALLOCATE(merge_files(m))
  !
  !Read the m file names that follow
  DO j=1,SIZE(merge_files)
    i=i+1
    READ(mode_param(i),'(a128)',ERR=7000,END=7000) temp
    merge_files(j) = temp
    !Check if that input file exists
    CALL CHECKFILE(merge_files(j),'read')
    CALL ATOMSK_MSG(999,(/merge_files(j)/),(/0.d0/))
  ENDDO
  !
  !The outputfile is file1
  i=i+1
  outputfile = file1
  CALL GUESS_FORMAT(outputfile,outfileformat,'writ')
  IF(outfileformat.NE.'xxx') CALL SET_OUTPUT(outfileformats,outfileformat,.TRUE.)
  !
  msg = 'outputfile: '//TRIM(outputfile)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Merge all files
  CALL MERGE_XYZ(merge_files(:),axis,options_array,outputfile,outfileformats)
!
!
!
500 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  MODE CREATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('create')
  !Create mode:
  !in this mode the program creates a structure
  !and writes it to one or several format(s)
  !Read parameters from array 'mode_param'
  READ(mode_param(1),*,END=7000,ERR=7000) create_struc
  !Get the lattice constant(s)
  READ(mode_param(2),*,END=7000,ERR=7000) create_a0(1)
  create_a0(2) = create_a0(1)
  create_a0(3) = create_a0(1)
  i=2
  IF( create_struc=='graphite' .OR. create_struc=='hcp' .OR.            &
    & create_struc=='wurtzite' .OR. create_struc=='wz'      ) THEN
    i=i+1
    READ(mode_param(i),*,END=7000,ERR=7000) create_a0(3)
  ELSEIF(create_struc=='nanotube' .OR. create_struc=='NT' .OR. create_struc=='nt') THEN
    i=i+1
    READ(mode_param(i),*,END=7000,ERR=7000) NT_mn(1)
    i=i+1
    READ(mode_param(i),*,END=7000,ERR=7000) NT_mn(2)
  ENDIF
  !Get the atomic species for the structure
  !There must always be at least one species
  i=i+1
  READ(mode_param(i),*,END=7000,ERR=7000) create_species(1)
  !Check that it is a correct atom type
  CALL ATOMNUMBER(create_species(1),smass)
  IF(smass==0.d0) THEN
    CALL ATOMSK_MSG(801,(/TRIM(create_species(1))/),(/0.d0/))
    GOTO 8000
  ENDIF
  !Check if other atoms are specified
  DO j=2,20
    i=i+1
    READ(mode_param(i),*,END=510,ERR=510) temp
    IF(LEN_TRIM(temp)<=2 .AND. temp(1:1).NE.'-' ) THEN
      READ(temp,*,END=510,ERR=510) species
      CALL ATOMNUMBER(species,smass)
      IF(smass==0.d0) THEN
        i=i-1
        GOTO 510
      ELSE
        create_species(j) = species
      ENDIF
    ELSE
      i=i-1
    ENDIF
  ENDDO
  510 CONTINUE
  !Check if a crystallographic orientation was given
  i=i+1
  READ(mode_param(i),*,END=520,ERR=520) temp
  IF(temp=="orient") THEN
    DO j=1,3
      i=i+1
      READ(mode_param(i),*,END=520,ERR=520) temp
      CALL INDEX_MILLER(temp,ORIENT(j,:),k)
      IF(k>0) GOTO 7000
    ENDDO
  ENDIF
  520 CONTINUE
  IF(ALLOCATED(mode_param)) DEALLOCATE(mode_param)
  !
  !Make sure we have an output file name
  outputfile = TRIM(ADJUSTL(file1))
  IF(LEN_TRIM(outputfile)==0) outputfile = TRIM(ADJUSTL(filefirst))
  IF(LEN_TRIM(outputfile)==0) THEN
    IF(create_struc=='L12') THEN
      outputfile = TRIM(create_species(1))//'3'//TRIM(create_species(2))
    ELSE
      DO i=1,SIZE(create_species)
        outputfile = TRIM(outputfile)//TRIM(create_species(i))
      ENDDO
      IF(create_struc=='per' .OR. create_struc=='perovskite') THEN
        outputfile = TRIM(outputfile)//'3'
      ENDIF
    ENDIF
  ENDIF
  msg = 'CREATE mode: outputfile: '//TRIM(outputfile)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !Activate the extension for this output file
  IF(LEN_TRIM(outputfile).NE.0) THEN
    CALL GUESS_FORMAT(outputfile,outfileformat,'writ')
    IF(outfileformat.NE.'xxx') CALL SET_OUTPUT(outfileformats,outfileformat,.TRUE.)
  ENDIF
  IF(nerr>0) GOTO 10000
  !
  !Create the unit cell
  CALL CREATE_CELL(create_a0,create_struc,create_species,NT_mn,ORIENT,options_array,outputfile,outfileformats,.TRUE.,H,P)
  IF(nerr>0) GOTO 10000
!
!
!
600 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  MODE ALL-IN-ONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('ai1')
  !All-in-one mode:
  !in this mode the program reads many files (no matter their format)
  !and concatenates them into one single file
  msg = 'ALL-IN-ONE mode: file: '//TRIM(file1)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Check that list file exists
  CALL CHECKFILE(file1,'read')
  !
  !If we have an output file name, activate its extension
  IF(LEN_TRIM(file2).NE.0) THEN
    CALL GUESS_FORMAT(file2,outfileformat,'writ')
    IF(outfileformat.NE.'xxx') CALL SET_OUTPUT(outfileformats,outfileformat,.TRUE.)
  ELSE
    strlength = SCAN(file1,'.',BACK=.TRUE.)
    IF(strlength.NE.0) THEN
      file2 = file1(1:strlength-1)
    ELSE
      file2 = file1
    ENDIF
  ENDIF
  msg = 'ALL-IN-ONE mode: files: '//TRIM(file1)//', '//TRIM(file2)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Call all-in-one module
  CALL ALLINONE(file1,file2,outfileformats,options_array)
!
!
!
700 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  MODE ONE-IN-ALL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('1ia')
  !One-in-all mode:
  !in this mode the program reads one file containing several
  !"snapshots" and writes each of them to a separate file
  !
  CALL ATOMSK_MSG(4020,(/TRIM(file1)/),(/0.d0/))
  !
  !Check if input file exists
  CALL CHECKFILE(file1,'read')
  !Guess its format
  CALL GUESS_FORMAT(file1,outfileformat,'read')
  IF(outfileformat=='xxx') THEN
    CALL ATOMSK_MSG(800,(/''/),(/0.d0/))
    nerr=nerr+1
    GOTO 10000
  ENDIF
  !
  !If we have an output file name, activate its extension
  IF(LEN_TRIM(file2).NE.0) THEN
    CALL GUESS_FORMAT(file2,outfileformat,'writ')
    IF(outfileformat.NE.'xxx') CALL SET_OUTPUT(outfileformats,outfileformat,.TRUE.)
  ELSE
    strlength = SCAN(file1,'.',BACK=.TRUE.)
    IF(strlength.NE.0) THEN
      file2 = file1(1:strlength-1)
    ELSE
      file2 = file1
    ENDIF
  ENDIF
  msg = 'ONE-IN-ALL mode: files: '//TRIM(file1)//', '//TRIM(file2)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Call corresponding 1-in-all module
  IF(outfileformat=='dlp') THEN
    CALL ONEINALL_DLP_HISTORY(file1,file2,outfileformats,options_array)
  ELSEIF(outfileformat=='pwo') THEN
    CALL ONEINALL_QEOUT(file1,file2,outfileformats,options_array)
  ELSEIF(outfileformat=='xsf') THEN
    CALL ONEINALL_XSF(file1,file2,outfileformats,options_array)
  ELSEIF(outfileformat=='xyz') THEN
    CALL ONEINALL_XYZ(file1,file2,outfileformats,options_array)
  ! --- add other formats here ---
  ELSE
    CALL ATOMSK_MSG(800,(/TRIM(outfileformat)/),(/0.d0/))
    nerr=nerr+1
    GOTO 10000
  ENDIF
!
!
!
800 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  MODE UNWRAP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('unwrap')
  msg = 'UNWRAP mode: file1,file2: '//TRIM(filefirst)//', '//TRIM(filesecond)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Unwrap mode:
  !in this mode the program reads two files
  !and attempts to unwrap coordinates of the second file
  !
  outputfile = TRIM(ADJUSTL(file1))
  !If we have an output file name, activate its extension
  IF(LEN_TRIM(outputfile).NE.0) THEN
    CALL GUESS_FORMAT(outputfile,outfileformat,'writ')
    IF(outfileformat.NE.'xxx') CALL SET_OUTPUT(outfileformats,outfileformat,.TRUE.)
  ELSE
    j=SCAN(inputfile,pathsep,BACK=.TRUE.)
    strlength = SCAN(filesecond,'.',BACK=.TRUE.)
    IF(strlength>j) THEN
      outputfile = filesecond(1:strlength-1)
    ELSE
      outputfile = filesecond
    ENDIF
  ENDIF
  msg = 'outputfile: '//TRIM(outputfile)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Unwrap positions in Psecond using Pfirst as reference
  CALL UNWRAP_XYZ(filefirst,filesecond,options_array,outputfile,outfileformats)
!
!
!
900 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  MODE POLYCRYSTAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('polycrystal')
  msg = 'POLYCRYSTAL mode: file1,file2: '//TRIM(file1)//', '//TRIM(file2)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  msg = 'POLYCRYSTAL mode: parameters: '//TRIM(filefirst)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Polycrystal mode:
  !in this mode the program reads one file containing atom positions,
  !one file containing parameters for the Voronoi construction,
  !and generates a polycrystal
  !
  outputfile = TRIM(ADJUSTL(file2))
  !If we have an output file name, activate its extension
  IF(LEN_TRIM(outputfile).NE.0) THEN
    CALL GUESS_FORMAT(outputfile,outfileformat,'writ')
    IF(outfileformat.NE.'xxx') CALL SET_OUTPUT(outfileformats,outfileformat,.TRUE.)
  ELSE
    j=SCAN(inputfile,pathsep,BACK=.TRUE.)
    strlength = SCAN(file1,'.',BACK=.TRUE.)
    IF(strlength>j) THEN
      outputfile = file1(1:strlength-1)
    ELSE
      outputfile = file1
    ENDIF
  ENDIF
  msg = 'outputfile: '//TRIM(outputfile)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  CALL POLYCRYS(file1,filefirst,options_array,file2,outfileformats)
!
!
!
1000 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  MODE INTERPOLATE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('interpolate')
  msg = 'INTERPOLATE mode: filefirst,filesecond: '//TRIM(filefirst)//', '//TRIM(filesecond)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Interpolate mode:
  !in this mode the program reads two files containing atom positions,
  !and builds a chain of N configurations (or images) between them
  !by a simple interpolation of atom positions
  !
  READ(mode_param(1),*,END=7000,ERR=7000) Nimages
  !
  outputfile = TRIM(ADJUSTL(file1))
  !If we have an output file name, activate its extension
  IF(LEN_TRIM(outputfile).NE.0) THEN
    CALL GUESS_FORMAT(outputfile,outfileformat,'writ')
    IF(outfileformat.NE.'xxx') CALL SET_OUTPUT(outfileformats,outfileformat,.TRUE.)
  ELSE
    j=SCAN(inputfile,pathsep,BACK=.TRUE.)
    strlength = SCAN(filefirst,'.',BACK=.TRUE.)
    IF(strlength>j) THEN
      outputfile = filefirst(1:strlength-1)
    ELSE
      outputfile = filefirst
    ENDIF
  ENDIF
  msg = 'outputfile: '//TRIM(outputfile)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  CALL INTERPOLATE_XYZ(filefirst,filesecond,Nimages,outputfile,outfileformats,options_array)
!
!
!
1100 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  MODE AVERAGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('average')
  msg = 'AVERAGE mode: listfile: '//TRIM(filefirst)//', '//TRIM(filesecond)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Average mode:
  !in this mode the program reads N files containing atom positions,
  !and averages the atom positions
  !
  outputfile = TRIM(ADJUSTL(file1))
  !If we have an output file name, activate its extension
  IF(LEN_TRIM(outputfile).NE.0) THEN
    CALL GUESS_FORMAT(outputfile,outfileformat,'writ')
    IF(outfileformat.NE.'xxx') CALL SET_OUTPUT(outfileformats,outfileformat,.TRUE.)
  ELSE
    j=SCAN(inputfile,pathsep,BACK=.TRUE.)
    strlength = SCAN(filefirst,'.',BACK=.TRUE.)
    IF(strlength>j) THEN
      outputfile = filefirst(1:strlength-1)
    ELSE
      outputfile = filefirst
    ENDIF
  ENDIF
  msg = 'outputfile: '//TRIM(outputfile)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  CALL AVERAGE_XYZ(listfile,outputfile,outfileformats,options_array)
!
!
!
!! === ADD OTHER MODES HERE ===
!
!
!
!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  !!!   Modes below perform some particular computation    !!!
!  !!!                                                      !!!
!
5000 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  MODE DIFFERENCE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('diff')
  CALL SET_OUTPUT(outfileformats,'all  ',.FALSE.)
  msg = 'DIFF mode: file1,file2: '//TRIM(file1)//', '//TRIM(file2)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Difference mode:
  !in this mode the program reads two files (no matter their format)
  !converts them to extended XYZ, and calculates the difference
  !in atomic positions
  CALL CHECKFILE(file1,'read')
  CALL CHECKFILE(file2,'read')
  !
  !Read first file
  CALL READ_AFF(file1,H,Pfirst,S,commentfirst,AUXNAMES,AUX)
  IF(nerr>=1) GOTO 10000
  CALL OPTIONS_AFF(options_array,H,Pfirst,S,AUXNAMES,AUX,ORIENT,SELECT)
  !Read second file
  CALL READ_AFF(file2,H,Psecond,S,commentsecond,AUXNAMES,AUX)
  IF(nerr>=1) GOTO 10000
  CALL OPTIONS_AFF(options_array,H,Psecond,S,AUXNAMES,AUX,ORIENT,SELECT)
  !
  !!Determine if the cell vectors were found or not
  IF( VECLENGTH(H(1,:))==0.d0 .OR. VECLENGTH(H(2,:))==0.d0 &
    & .OR. VECLENGTH(H(3,:))==0.d0 ) THEN
    !!If basis vectors H were not found in the input file, then compute them
    ! (this does not work perfectly yet, see determine_H.f90)
    CALL ATOMSK_MSG(4701,(/''/),(/0.d0/))
    !
    CALL DETERMINE_H(H,Pfirst)
    !
  ENDIF
  !
  !Set the name for output file
  strlength = SCAN(file2,'.',BACK=.TRUE.)
  IF(strlength.NE.0) THEN
    test = file2(1:strlength-1)
  ELSE
    test = file2
  ENDIF
  !
  !Calculate the difference between the two sets of positions
  CALL DIFF_XYZ(H,Pfirst,Psecond,test)
!
!
!
5100 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  MODE ELECTRIC DIPOLE MOMENTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('edm')
  !EDM mode:
  !in this mode the program reads one file
  !and computes individual electric dipole moments for
  !one or several kinds of polyhedra
  !Read parameters from 'mode_param'
  READ(mode_param(1),*,END=7000,ERR=7000) NNN !Number of nearest neighbours
  READ(mode_param(2),*,END=7000,ERR=7000) Pspecies !Atoms forming polyhedra
  IF(ALLOCATED(mode_param)) DEALLOCATE(mode_param)
  !
  outputfile = TRIM(ADJUSTL(file1))
  IF(outputfile=='') outputfile = TRIM(ADJUSTL(filefirst))
  IF(LEN_TRIM(outputfile).NE.0) THEN
    CALL GUESS_FORMAT(outputfile,outfileformat,'writ')
    IF(outfileformat.NE.'xxx') CALL SET_OUTPUT(outfileformats,outfileformat,.TRUE.)
  ENDIF
  msg = 'EDM mode file: '//TRIM(file1)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Check if input file exists
  CALL CHECKFILE(file1,'read')
  !Read the file
  CALL READ_AFF(file1,H,Pfirst,S,commentfirst,AUXNAMES,AUX)
  IF(nerr>0) GOTO 10000
  !
  !!Determine if the cell vectors were found or not
  IF( VECLENGTH(H(1,:))==0.d0 .OR. VECLENGTH(H(2,:))==0.d0 &
    & .OR. VECLENGTH(H(3,:))==0.d0 ) THEN
    !!If basis vectors H were not found in the input file, then compute them
    ! (this does not work perfectly yet, see determine_H.f90)
    CALL ATOMSK_MSG(4701,(/''/),(/0.d0/))
    !
    CALL DETERMINE_H(H,Pfirst)
    !
  ENDIF
  !
  !Check if user asked for wrapping of atoms
  CALL CHECK_OPTION_WRAP(options_array)
  !
  !Set the name for output file
  strlength = SCAN(file1,'.',BACK=.TRUE.)
  IF(strlength.NE.0) THEN
    test = file1(1:strlength-1)
  ELSE
    test = file1
  ENDIF
  !
  !Apply options
  CALL OPTIONS_AFF(options_array,H,Pfirst,S,AUXNAMES,AUX,ORIENT,SELECT)
  IF(nerr>0) GOTO 10000
  !
  !Compute polarization
  CALL E_DIPOLES(H,Pfirst,S,AUXNAMES,AUX,NNN,Pspecies,commentfirst,test)
!
!
!
5200 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  MODE ELECTRONIC POLARIZATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('PE')
  !Electronic polarization mode:
  !in this mode the program reads one file
  !and computes electronic polarization of core-shell ions
  !
  outputfile = TRIM(ADJUSTL(file1))
  IF(outputfile=='') outputfile = TRIM(ADJUSTL(filefirst))
  IF(LEN_TRIM(outputfile).NE.0) THEN
    CALL GUESS_FORMAT(outputfile,outfileformat,'writ')
    IF(outfileformat.NE.'xxx') CALL SET_OUTPUT(outfileformats,outfileformat,.TRUE.)
  ENDIF
  msg = 'ELECTRONIC POLARIZATION mode file: '//TRIM(file1)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Check if input file exists
  CALL CHECKFILE(file1,'read')
  !Read the file
  CALL READ_AFF(file1,H,Pfirst,S,commentfirst,AUXNAMES,AUX)
  IF(nerr>0) GOTO 10000
  !
  !!Determine if the cell vectors were found or not
  IF( VECLENGTH(H(1,:))==0.d0 .OR. VECLENGTH(H(2,:))==0.d0 &
    & .OR. VECLENGTH(H(3,:))==0.d0 ) THEN
    !!If basis vectors H were not found in the input file, then compute them
    ! (this does not work perfectly yet, see determine_H.f90)
    CALL ATOMSK_MSG(4701,(/''/),(/0.d0/))
    !
    CALL DETERMINE_H(H,Pfirst)
    !
  ENDIF
  !
  !Set the name for output file
  strlength = SCAN(file1,'.',BACK=.TRUE.)
  IF(strlength.NE.0) THEN
    test = file1(1:strlength-1)
  ELSE
    test = file1
  ENDIF
  !
  !Apply options
  CALL OPTIONS_AFF(options_array,H,Pfirst,S,AUXNAMES,AUX,ORIENT,SELECT)
  IF(nerr>0) GOTO 10000
  !
  !Compute polarization
  CALL E_POLA(H,Pfirst,S,AUXNAMES,AUX,comment,test)
!
!
!
5300 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  MODE RADIAL DISTRIBUTION FUNCTION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('rdf')
  !Radial Distribution Function (RDF) mode:
  !in this mode the program reads a list of files,
  !computes the space-averaged RDF for each of them,
  !and then makes a time-average
  !
  READ(mode_param(1),*,END=7000,ERR=7000) rdf_maxR !max. radius R
  READ(mode_param(2),*,END=7000,ERR=7000) rdf_dr   !skin radius dR
  !
  !Check if user asked for wrapping of atoms
  CALL CHECK_OPTION_WRAP(options_array)
  !
  !Compute polarization
  CALL RDF_XYZ(listfile,rdf_maxR,rdf_dr,options_array)
!
!
5400 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  MODE NYE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('nye')
  !in this mode the code reads two files and computes the Nye tensor
  !
  !Check if user asked for wrapping of atoms
  CALL CHECK_OPTION_WRAP(options_array)
  !
  CALL NYE_TENSOR(filefirst,filesecond,options_array,file1,outfileformats)
!
!
!
5500 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  DENSITY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('density')
  !DENSITY mode:
  !in this mode the program reads one file
  !and computes the density of a given property
  msg = 'DENSITY mode file: '//TRIM(file1)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !Read parameters from 'mode_param'
  READ(mode_param(1),*,END=7000,ERR=7000) property !name of the property whose density will be calculated
  READ(mode_param(2),*,END=7000,ERR=7000) j        !type of density to output: 1, 2 or 3 (for 1-D, 2-D or 3-D)
  IF( j==1 .OR. j==2 ) THEN
    READ(mode_param(3),*,END=7000,ERR=7000) axis
  ENDIF
  READ(mode_param(4),*,END=7000,ERR=7000) NNN      !Sigma = square of variance
  IF(ALLOCATED(mode_param)) DEALLOCATE(mode_param)
  !
  !Set the name for output file
  strlength = SCAN(file1,'.',BACK=.TRUE.)
  IF(strlength.NE.0) THEN
    test = file1(1:strlength-1)
  ELSE
    test = file1
  ENDIF
  !
  !Check if user asked for wrapping of atoms
  CALL CHECK_OPTION_WRAP(options_array)
  !
  !Compute polarization
  CALL DENSITY_XYZ(file1,options_array,property,j,axis,NNN,test)
!
!
!
5600 CONTINUE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!  CENTRAL SYMMETRY PARAMETER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('cs')
  !CENTRAL SYMMETRY PARAMETER mode:
  !in this mode the program reads one file
  !and computes the central symmetry parameter for each atom
  msg = 'CENTRAL SYMMETRY PARAMETER mode file: '//TRIM(filefirst)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  outputfile = TRIM(ADJUSTL(file1))
  !If we have an output file name, activate its extension
  IF(LEN_TRIM(outputfile).NE.0) THEN
    CALL GUESS_FORMAT(outputfile,outfileformat,'writ')
    IF(outfileformat.NE.'xxx') CALL SET_OUTPUT(outfileformats,outfileformat,.TRUE.)
  ELSE
    j=SCAN(inputfile,pathsep,BACK=.TRUE.)
    strlength = SCAN(filefirst,'.',BACK=.TRUE.)
    IF(strlength>j) THEN
      outputfile = filefirst(1:strlength-1)
    ELSE
      outputfile = filefirst
    ENDIF
  ENDIF
  msg = 'outputfile: '//TRIM(outputfile)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Check if user asked for wrapping of atoms
  CALL CHECK_OPTION_WRAP(options_array)
  !
  !Compute central symmetry parameter
  CALL CENTRO_SYM(filefirst,options_array,outputfile,outfileformats)
!
!
!
CASE DEFAULT
  GOTO 7000
END SELECT  !END OF MODES
!
GOTO 10000
!
!
!
7000 CONTINUE
CALL ATOMSK_MSG(4800,(/TRIM(mode)/),(/0.d0/))
temp = TRIM(ADJUSTL(mode))
CALL DISPLAY_HELP(temp)
!
!
!
8000 CONTINUE
nerr=nerr+1
!
!
!
10000 CONTINUE
IF(ALLOCATED(Pfirst)) DEALLOCATE(Pfirst)
IF(ALLOCATED(Psecond)) DEALLOCATE(Psecond)
!
!
END SUBROUTINE RUN_MODE
!
!
!
SUBROUTINE CHECK_OPTION_WRAP(options_array)
!
IMPLICIT NONE
CHARACTER(LEN=1):: answer !"y" or "n"
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: options_array !options and their parameters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_temp  !options and their parameters (temporary)
INTEGER:: i
!
!Check if user asked for wrapping of atoms
IF( .NOT.ALLOCATED(options_array) ) THEN
  !User did not ask to wrap atoms => warn and ask
  nwarn=nwarn+1
  CALL ATOMSK_MSG(4712,(/""/),(/0.d0/))
  READ(*,*) answer
  IF( answer==langyes .OR. answer==langBigYes ) THEN
    ALLOCATE(options_array(1))
    options_array(1) = "-wrap"
  ENDIF
  !
ELSEIF( .NOT.ANY(options_array=="-wrap") ) THEN
  !User did not ask to wrap atoms => warn and ask
  nwarn=nwarn+1
  CALL ATOMSK_MSG(4712,(/""/),(/0.d0/))
  READ(*,*) answer
  IF( answer==langyes .OR. answer==langBigYes ) THEN
    !Add option "-wrap" to options_array (at the beginning)
    ALLOCATE( options_temp(SIZE(options_array)+1) )
    options_temp(:) = ""
    DO i=1,SIZE(options_array)
      options_temp(i+1) = options_array(i)
    ENDDO
    DEALLOCATE(options_array)
    ALLOCATE(options_array(SIZE(options_temp)))
    options_array(:) = options_temp(:)
    DEALLOCATE(options_temp)
    options_array(1) = "-wrap"
  ENDIF
ENDIF
!
END SUBROUTINE CHECK_OPTION_WRAP
!
!
!
END MODULE modes
