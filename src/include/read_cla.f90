MODULE read_cla
!
!**********************************************************************************
!*  READ_CLA                                                                      *
!**********************************************************************************
!* This module reads command-line arguments for ATOMSK.                           *
!**********************************************************************************
!* (C) March 2011 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 04 April 2018                                    *
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
SUBROUTINE GET_CLA(cla,mode,options_array,outfileformats,pfiles,mode_param)
!
IMPLICIT NONE
CHARACTER(LEN=2):: species
CHARACTER(LEN=5),DIMENSION(100):: tempout !temporary list of formats to write
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: outfileformats !list of formats to write
CHARACTER(LEN=12),INTENT(OUT):: mode     !mode in which the program runs
CHARACTER(LEN=16):: region_geom  !geometry of the region: "box" or "sphere"
CHARACTER(LEN=4096):: clarg      !one command-line argument
CHARACTER(LEN=4096):: msg, temp, temp2, temp3, temp4
CHARACTER(LEN=4096),DIMENSION(5):: pfiles !pfiles(1)=file1
                                          !pfiles(2)=file2
                                          !pfiles(3)=filefirst
                                          !pfiles(4)=filesecond
                                          !pfiles(5)=listfile
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: options_array !options and their parameters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: mode_param  !parameters for some special modes
CHARACTER(LEN=4096),DIMENSION(:),ALLOCATABLE,INTENT(IN):: cla  !command-line arguments
INTEGER:: i, ioptions, j, m
INTEGER:: Nout !number of output formats
REAL(dp):: smass, tempreal
!
!
!Initialize variables
pfiles(:)=''
i=0
ioptions = 0
Nout = 0
IF(ALLOCATED(mode_param)) DEALLOCATE(mode_param)
IF(ALLOCATED(outfileformats)) DEALLOCATE(outfileformats)
!
IF( ALLOCATED(options_array) .AND. SIZE(options_array)>=1 ) THEN
  options_array(:) = ""
ENDIF
!
!
IF(verbosity==4) THEN
  !Debug messages
  msg = "Entering GET_CLA"
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  DO i=1,SIZE(cla)
    WRITE(msg,'(a18,i3,3X,a2,a)') 'Command-line arg #',        &
                            & i, ': ', TRIM(ADJUSTL(cla(i)))
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  ENDDO
ENDIF
!
IF( .NOT.ALLOCATED(cla) .OR. SIZE(cla)<=0 ) THEN
  GOTO 1000
ENDIF
!
!
i=0
DO WHILE(i<SIZE(cla))
  i=i+1
  !Initialize variables
  clarg = ''
  msg = ''
  temp = ''
  tempreal = 0.d0
  !
  clarg = TRIM(ADJUSTL( cla(i) ))
  IF(LEN_TRIM(clarg)==0) EXIT
  !
  !
  !Deal with modes
  IF(clarg=='--all-in-one' .OR. clarg=='-AI1') THEN
    mode = 'ai1'
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(1)
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(2)
    !
  ELSEIF(clarg=='--average') THEN
    mode = 'average'
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(5)
    IF(LEN_TRIM(pfiles(5))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  ELSEIF(clarg=='--central-symmetry' .OR. clarg=="--cs") THEN
    mode = 'cs'
    i=i+1
    READ(cla(i),'(a4096)',END=130,ERR=130) pfiles(3) !File containing atom positions
    IF(LEN_TRIM(pfiles(3))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  ELSEIF(clarg=='--create' .OR. clarg=='-C') THEN
    IF(ALLOCATED(mode_param)) GOTO 150
    ALLOCATE(mode_param(20))
    mode_param(:) = ''
    mode='create'
    !Get the type of structure
    i=i+1
    m=1
    READ(cla(i),*,END=130,ERR=130) temp
    !Make the mode_param(1) variable straight
    !in case the user mistyped 'diamand' or 'pervoskite'
    SELECT CASE(temp(1:2))
    CASE('di','Di','DI')
      temp = 'diamond'
    CASE('pe','Pe','PE')
      temp = 'perovskite'
    CASE('zi','Zi','zb','ZB')
      temp = 'zincblende'
    CASE('na','NA','nt','NT')
      temp = 'nanotube'
    CASE('wu','Wu','WU','wz','Wz','WZ')
      temp = 'wurtzite'
    END SELECT
    mode_param(m) = TRIM(ADJUSTL(temp))
    !Get the lattice constant a
    i=i+1
    m=m+1
    READ(cla(i),*,END=130,ERR=130) mode_param(m)
    IF( mode_param(1)=='graphite' .OR. mode_param(1)=='hcp' .OR.            &
      & mode_param(1)=='wurtzite' .OR. mode_param(1)=='wz' .OR. mode_param(1)=='c14'      ) THEN
      !Get lattice constant c
      i=i+1
      m=m+1
      READ(cla(i),*,END=130,ERR=130) mode_param(m)
    ELSEIF(mode_param(1)=='nanotube' .OR. mode_param(1)=='NT' .OR. mode_param(1)=='nt') THEN
      !Get chiral indices (m,n) of the nanotube
      i=i+1
      m=m+1
      READ(cla(i),*,END=130,ERR=130) mode_param(m)
      i=i+1
      m=m+1
      READ(cla(i),*,END=130,ERR=130) mode_param(m)
    ENDIF
    !Get the atomic species for the structure
    DO
      i=i+1
      temp = ADJUSTL(cla(i))
      IF( LEN_TRIM(temp)>0 ) THEN
        IF(LEN_TRIM(temp)<=2 .AND. temp(1:1).NE.'-' ) THEN
          !it may be an atom species => try it
          m=m+1
          READ(temp,*,ERR=130,END=130) species
          CALL ATOMNUMBER(species,smass)
          IF(smass==0.d0) THEN
            !nope, not an atom species
            i=i-1
            GOTO 110
          ELSE
            !save this atom species and try the next argument
            WRITE(mode_param(m),*) species
          ENDIF
        ELSE
          !It is not an atom species => go back and exit the loop
          i=i-1
          EXIT
        ENDIF
      ELSE
        !Empty entry => exit
        EXIT
      ENDIF
    ENDDO
    !Look for crystallographic orientation (if any)
    i=i+1
    temp = ADJUSTL(cla(i))
    IF( temp=="orient" ) THEN
      SELECT CASE(mode_param(1))
      CASE('sc','SC','fcc','FCC','L12','L1_2','bcc','BCC','CsCl','diamond','dia','zincblende','zb','ZB','B3', &
          & 'perovskite','per','rocksalt','rs','RS','B1','fluorite','fluorine')
        !Cubic lattice => read crystal orientation
        m=m+1
        mode_param(m) = "orient"
        !Get Miller indices [hkl] of the crystal orientation along X, Y, Z
        i=i+1
        m=m+1
        READ(cla(i),*,END=130,ERR=130) mode_param(m)
        i=i+1
        m=m+1
        READ(cla(i),*,END=130,ERR=130) mode_param(m)
        i=i+1
        m=m+1
        READ(cla(i),*,END=130,ERR=130) mode_param(m)
      CASE('hcp','HCP','wurtzite','wz','WZ','graphite','c14')
        !Hexagonal lattice => read crystal orientation
        m=m+1
        mode_param(m) = "orient"
        !Get Miller indices [hkil] of the crystal orientation
        i=i+1
        m=m+1
        READ(cla(i),*,END=130,ERR=130) mode_param(m)
        i=i+1
        m=m+1
        READ(cla(i),*,END=130,ERR=130) mode_param(m)
        i=i+1
        m=m+1
        READ(cla(i),*,END=130,ERR=130) mode_param(m)
      CASE DEFAULT
        !Lattice is not cubic: crystal orientation not implemented yet
        !Display error message and exit
        CALL ATOMSK_MSG(4827,(/""/),(/0.d0/))
        nerr = nerr+1
        GOTO 1000
      END SELECT
    ELSE
      i=i-1
    ENDIF
    IF(LEN_TRIM(mode_param(1))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  ELSEIF(clarg=='--density') THEN
    IF(ALLOCATED(mode_param)) GOTO 150
    ALLOCATE(mode_param(4))
    mode_param(:) = " "
    mode = 'density'
    i=i+1
    READ(cla(i),'(a4096)',END=130,ERR=130) pfiles(1) !File containing atom positions
    IF(LEN_TRIM(pfiles(1))==0 .OR. i>SIZE(cla)) GOTO 130
    i=i+1
    READ(cla(i),*,END=130,ERR=130) mode_param(1)  !property whose density will be calculated
    i=i+1
    READ(cla(i),*,END=130,ERR=130) temp  !den_type: 1, 2 or 3 (for 1-D, 2-D or 3-D density)
    SELECT CASE(temp)
    CASE("1","1d","1D","1-d","1-D")
      j=1
    CASE("2","2d","2D","2-d","2-D")
      j=2
    CASE("3","3d","3D","3-d","3-D")
      j=3
    CASE DEFAULT
      !It is neither 1-D, 2-D nor 3-D => don't know what it is, abort
      j=-1
      nerr = nerr+1
      GOTO 1000
    END SELECT
    WRITE(mode_param(2),*) j
    IF( j==1 .OR. j==2 ) THEN
      i=i+1
      READ(cla(i),*,END=130,ERR=130) mode_param(3)  !axis: x, y or z
    ENDIF
    i=i+1
    READ(cla(i),*,END=130,ERR=130) tempreal  !Sigma = square of variance
    WRITE(mode_param(4),*) tempreal
    IF( i>SIZE(cla) ) GOTO 130
    !
  ELSEIF(clarg=='--difference' .OR. clarg=='--diff' .OR. clarg=='-D') THEN
    mode = 'diff'
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(1)
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(2)
    IF(LEN_TRIM(pfiles(1))==0 .OR. LEN_TRIM(pfiles(2))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  ELSEIF(clarg=='--ddplot') THEN
    mode = 'ddplot'
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(3)
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(4)
    IF(LEN_TRIM(pfiles(3))==0 .OR. LEN_TRIM(pfiles(4))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  ELSEIF(clarg=='--edm') THEN
    IF(ALLOCATED(mode_param)) GOTO 150
    ALLOCATE(mode_param(20))
    mode_param(:) = ''
    mode = 'edm'
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(1) !File containing atom positions
    i=i+1
    READ(cla(i),*,END=130,ERR=130) species  !Atom species forming polyhedra
    WRITE(mode_param(2),*) species
    CALL ATOMNUMBER(species,tempreal)
    i=i+1
    READ(cla(i),*,END=130,ERR=130) tempreal  !NNN
    WRITE(mode_param(1),*) tempreal
    IF(LEN_TRIM(pfiles(1))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  ELSEIF(clarg=='--electronic-polarization' .OR. clarg=='-PE') THEN
    mode = 'PE'
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(1) !File containing atom positions
    IF(LEN_TRIM(pfiles(1))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  ELSEIF(clarg=='--interpolate') THEN
    IF(ALLOCATED(mode_param)) GOTO 150
    ALLOCATE(mode_param(1))
    mode = 'interpolate'
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(3) !Atom positions of first configuration
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(4) !Atom positions of final configuration
    i=i+1
    READ(cla(i),*,END=130,ERR=130) m  !number of images to interpolate
    WRITE(mode_param(1),*) m
    IF(LEN_TRIM(pfiles(3))==0 .OR. LEN_TRIM(pfiles(4))==0 .OR. m<=0 .OR. i>SIZE(cla)) GOTO 130
    !
  ELSEIF(clarg=='--list' .OR. clarg=='-L') THEN
    mode = 'list'
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(5)
    IF(LEN_TRIM(pfiles(5))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  ELSEIF(clarg=='--merge' .OR. clarg=='-M') THEN
    IF(ALLOCATED(mode_param)) GOTO 150
    ALLOCATE(mode_param(100))
    mode_param(:) = ''
    mode = 'merge'
    !Store all parameters and file names in mode_param(:)
    j=0
    i=i+1
    READ(cla(i),*,ERR=130,END=130) temp
    temp = TRIM(ADJUSTL(temp))
    IF( temp=='x' .OR. temp=='y' .OR. temp=='z' .OR. &
      & temp=='X' .OR. temp=='Y' .OR. temp=='Z'      ) THEN
      j=j+1
      mode_param(j) = temp
      i=i+1
      READ(cla(i),*,ERR=130,END=130) m
      j=j+1
      WRITE(mode_param(j),*) m
      m=m+2
    ELSE
      READ(temp,*,ERR=130,END=130) m
      j=j+1
      WRITE(mode_param(j),*) m
      m=m+1
    ENDIF
    !Read the m file names that follow
    DO WHILE(j<m)
      i=i+1
      READ(cla(i),'(a128)',ERR=130,END=130) temp
      j=j+1
      mode_param(j) = TRIM(ADJUSTL(temp))
    ENDDO
    !Read name of the output file
    i=i+1
    READ(cla(i),*,ERR=130,END=130) pfiles(1)
    IF(LEN_TRIM(pfiles(1))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  ELSEIF(clarg=='--normal') THEN
    !nothing special to do, this is normal mode
    !
  ELSEIF(clarg=='--nye') THEN
    mode = 'nye'
    i=i+1
    READ(cla(i),'(a)',END=130,ERR=130) pfiles(3)
    i=i+1
    READ(cla(i),'(a)',END=130,ERR=130) pfiles(4)
    !
  ELSEIF(clarg=='--one-in-all' .OR. clarg=='-1IA') THEN
    mode = '1ia'
    i=i+1
    READ(cla(i),'(a128)',ERR=130,END=130) pfiles(1)
    !
  ELSEIF(clarg=='--polycrystal') THEN
    mode = 'polycrystal'
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(1) !File containing atom positions
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(3) !File containing parameters for Voronoi construction
    IF(LEN_TRIM(pfiles(1))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  ELSEIF(clarg=='--rdf') THEN
    IF(ALLOCATED(mode_param)) GOTO 150
    ALLOCATE(mode_param(2))
    mode = 'rdf'
    i=i+1
    READ(cla(i),'(a4096)',END=130,ERR=130) pfiles(5) !File list
    IF(LEN_TRIM(pfiles(5))==0 .OR. i>SIZE(cla)) GOTO 130
    i=i+1
    READ(cla(i),*,END=130,ERR=130) mode_param(1)  !R
    i=i+1
    READ(cla(i),*,END=130,ERR=130) mode_param(2)  !dR
    IF( i>SIZE(cla) ) GOTO 130
    !
  ELSEIF(clarg=='--unwrap') THEN
    IF(ALLOCATED(mode_param)) GOTO 150
    mode = 'unwrap'
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(3) !"Reference"
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(4) !"Configuration" to unwrap
    IF(LEN_TRIM(pfiles(3))==0 .OR. LEN_TRIM(pfiles(4))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  ELSEIF(clarg(1:2)=='--') THEN
    !if it starts with "--" we assume it is a wrong mode entered by the user
    nerr=nerr+1
    CALL ATOMSK_MSG(4813,(/TRIM(clarg)/),(/0.d0/))
    GOTO 1000
  !
  !
  !Deal with output formats
  ELSEIF(clarg=='atsk' .OR. clarg=='atomsk' .OR. clarg=='ATSK' .OR. clarg=='ATOMSK' .OR. clarg=='Atomsk') THEN
    Nout = Nout+1
    tempout(Nout) = 'atsk'
  ELSEIF(clarg=='bop' .OR. clarg=="BOP") THEN
    Nout = Nout+1
    tempout(Nout) = 'bop'
  ELSEIF( clarg=='cfg' .OR. clarg=='CFG' .OR. clarg=='atomeye' .OR. clarg=='Atomeye' .OR. &
        & clarg=='ATOMEYE' .OR. clarg=='QSTEM' .OR. clarg=='qstem' .OR. clarg=='Qstem' ) THEN
    Nout = Nout+1
    tempout(Nout) = 'cfg'
  ELSEIF(clarg=='cel' .OR. clarg=='CEL' .OR. clarg=='drprobe' .OR. clarg=='DrProbe') THEN
    Nout = Nout+1
    tempout(Nout) = 'cel'
  ELSEIF(clarg=='cif' .OR. clarg=='CIF') THEN
    Nout = Nout+1
    tempout(Nout) = 'cif'
  ELSEIF(clarg=='cml' .OR. clarg=='CML') THEN
    Nout = Nout+1
    tempout(Nout) = 'cml'
  ELSEIF(clarg=='dlp' .OR. clarg=='DLP' .OR. clarg=='dlpoly' .OR. clarg=='DLPOLY' .OR.  &
        & clarg=="dl_poly" .OR. clarg=="DL_POLY" ) THEN
    Nout = Nout+1
    tempout(Nout) = 'dlp'
  ELSEIF(clarg=='config' .OR. clarg=='CONFIG') THEN
    IF(LEN_TRIM(pfiles(1)).NE.0) THEN
      IF(LEN_TRIM(pfiles(2))==0) pfiles(2) = 'CONFIG'
      Nout = Nout+1
      tempout(Nout) = 'dlp'
    ELSE
      pfiles(1) = 'CONFIG'
    ENDIF
  ELSEIF(clarg=='coo'.OR. clarg=='COO' .OR. clarg=='mbpp' .OR. clarg=='MBPP') THEN
    Nout = Nout+1
    tempout(Nout) = 'coo'
  ELSEIF(clarg=='coorat' .OR. clarg=='COORAT' .OR. clarg=='mbpp' .OR. clarg=='MBPP') THEN
    !if a file name was defined before, 
    IF(LEN_TRIM(pfiles(1)).NE.0) THEN
      IF(LEN_TRIM(pfiles(2))==0) pfiles(2) = 'COORAT'
      Nout = Nout+1
      tempout(Nout) = 'coo'
    ELSE
      pfiles(1) = 'COORAT'
    ENDIF
  ELSEIF(clarg=='dd' .OR. clarg=='DD' .OR. clarg=='ddplot') THEN
    Nout = Nout+1
    tempout(Nout) = 'dd'
  ELSEIF(clarg=='exyz' .OR. clarg=='EXYZ') THEN
    Nout = Nout+1
    tempout(Nout) = 'exyz'
  ELSEIF(clarg=='gin' .OR. clarg=='GIN' .OR. clarg=='gulp' .OR. clarg=='GULP') THEN
    Nout = Nout+1
    tempout(Nout) = 'gin'
  ELSEIF(clarg=='imd' .OR. clarg=='IMD') THEN
    Nout = Nout+1
    tempout(Nout) = 'imd'
  ELSEIF(clarg=='jems' .OR. clarg=='JEMS' .OR. clarg=='Jems') THEN
    Nout = Nout+1
    tempout(Nout) = 'jems'
  ELSEIF(clarg=='lammps' .OR. clarg=='LAMMPS' .OR. clarg=='lmp' .OR. clarg=='LMP') THEN
    Nout = Nout+1
    tempout(Nout) = 'lmp'
  ELSEIF( clarg=='mol' .OR. clarg=='MOL' .OR. clarg=='moldy' &
        & .OR. clarg=='mold' .OR. clarg=='MOLDY'                 ) THEN
    Nout = Nout+1
    tempout(Nout) = 'mol'
  ELSEIF(clarg=='pdb' .OR. clarg=='PDB') THEN
    Nout = Nout+1
    tempout(Nout) = 'pdb'
  ELSEIF(clarg=='pos' .OR. clarg=='POS' .OR. clarg=='vasp' .OR. clarg=='VASP') THEN
    Nout = Nout+1
    tempout(Nout) = 'pos'
  ELSEIF(clarg=='poscar' .OR. clarg=='POSCAR') THEN
    IF(LEN_TRIM(pfiles(1)).NE.0) THEN
      IF(LEN_TRIM(pfiles(2))==0) pfiles(2) = 'POSCAR'
      Nout = Nout+1
      tempout(Nout) = 'pos'
    ELSE
      IF(LEN_TRIM(pfiles(2))==0) pfiles(1) = 'POSCAR'
    ENDIF
  ELSEIF( clarg=='pw' .OR. clarg=='PW' .OR. clarg=='pwscf' .OR. clarg=='PWscf' .OR. &
        & clarg=='PWSCF' .OR. clarg=='qe' .OR. clarg=='QE' ) THEN
    Nout = Nout+1
    tempout(Nout) = 'pw'
  ELSEIF(clarg=='vesta' .OR. clarg=='VESTA') THEN
    Nout = Nout+1
    tempout(Nout) = 'vesta'
  ELSEIF(clarg=='xmd' .OR. clarg=='XMD') THEN
    Nout = Nout+1
    tempout(Nout) = 'xmd'
  ELSEIF(clarg=='xsf' .OR. clarg=='XSF' .OR. clarg=='xcrysden' .OR. clarg=='XCrySDen' &
        & .OR. clarg=='XCrysDen'.OR. clarg=='XCRYSDEN') THEN
    Nout = Nout+1
    tempout(Nout) = 'xsf'
  ELSEIF(clarg=='sxyz' .OR. clarg=='SXYZ') THEN
    Nout = Nout+1
    tempout(Nout) = 'sxyz'
  ELSEIF(clarg=='xv' .OR. clarg=='XV' .OR. clarg=='siesta' .OR. clarg=='Siesta' .OR. clarg=='SIESTA') THEN
    Nout = Nout+1
    tempout(Nout) = 'xv'
  ELSEIF(clarg=='xyz' .OR. clarg=='XYZ') THEN
    Nout = Nout+1
    tempout(Nout) = 'xyz'
  ! -- please add new formats in alphabetical order --
  !
  !
  !Deal with special options
  !For each option we store the option parameters in the string array "options_array"
  ELSEIF(clarg=='-add-atom' .OR. clarg=='-add-atoms' .OR. clarg=='-addatom' .OR. clarg=='-addatoms') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Read the species of atom to add
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !Read keyword: must be "at", "near" or "random"
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    IF( .NOT.( temp(1:2)=="at" .OR. temp(1:2)=="AT" .OR. temp(1:1)=="@" .OR.  &
      &        temp(1:4)=="near" .OR. temp(1:4)=="NEAR" .OR.                  &
      &        temp(1:6)=="random" .OR. temp(1:6)=="RANDOM" .OR.              &
      &        temp(1:8)=="relative" .OR. temp(1:3)=="rel"                    &
      &      ) ) THEN
      GOTO 400
    ENDIF
    temp2 = temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !Read parameters
    SELECT CASE(temp2)
    CASE("at","AT","@")
      !Read coordinate x of atom to add
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !Read coordinate y of atom to add
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !Read coordinate z of atom to add
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !
    CASE("relative","rel")
      !Read index of atom near which the new atom must be added
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      READ(temp,*,END=120,ERR=120) tempreal
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !Read coordinate x of atom to add
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !Read coordinate y of atom to add
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !Read coordinate z of atom to add
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !
    CASE("near","NEAR")
      !Read index of atom near which the new atom must be added
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      READ(temp,*,END=120,ERR=120) tempreal
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !
    CASE("random")
      !Read how many atoms must be inserted randomly in the system
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      READ(temp,*,END=120,ERR=120) tempreal
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !
    END SELECT
  !
  ELSEIF( clarg=='-add-shells' .OR. clarg=='-addshells' .OR. clarg=='-as' .OR. &
        & clarg=='-create-shells' .OR. clarg=='-cs') THEN
    IF( clarg=='-create-shells' .OR. clarg=='-cs') THEN
      nwarn=nwarn+1
      CALL ATOMSK_MSG(2799,(/'-create-shells','-add-shells   '/),(/0.d0/))
    ENDIF
    ioptions = ioptions+1
    options_array(ioptions) = '-add-shells'
    !read atom species for which the shells must be created
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  ELSEIF(clarg=='-alignx') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
  !
  ELSEIF(clarg=='-bind-shells' .OR. clarg=='-bs') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
  !
  ELSEIF(clarg=='-crack') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Read crack mode (I, II or III)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF(temp.NE.'I' .AND. temp.NE.'II' .AND. temp.NE.'III') GOTO 120
    !Read type of crack displacements: "stress" or "strain"
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF(temp.NE.'stress' .AND. temp.NE.'strain') GOTO 120
    !Read stress intensity factor K
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    READ(temp,*,END=120,ERR=120) tempreal
    !Read first coordinate of crack tip (pos1)
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
      & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
      & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
      & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') GOTO 120
    !Read normal to plane of cut (X, Y or Z)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
      & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') GOTO 120
    !Read shear modulus (mu)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    READ(temp,*,END=120,ERR=120) tempreal
    !Read Poisson ratio (mu)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    READ(temp,*,END=120,ERR=120) tempreal
    !if slashes are present (user wants to perform a division), replace them by a colon (:)
    j = SCAN(options_array(ioptions),'/')
    DO WHILE(j>0)
      options_array(ioptions)(j:j) = ':'
      j = SCAN(options_array(ioptions),'/')
    ENDDO
  !
  ELSEIF(clarg=='-cell') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  ELSEIF(clarg=='-center') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    IF(temp=="com") WRITE(temp,'(a1)') "0"
    READ(temp,*,END=120,ERR=120) j
    IF(j<0) WRITE(temp,'(a1)') "0"
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  ELSEIF(clarg=='-cut') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read 'above' or 'below'
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF( temp(1:5).NE.'above' .AND. temp(1:5).NE.'below' ) GOTO 120
    !read cut distance
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
      & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read cut direction
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  ELSEIF(clarg=='-deform' .OR. clarg=='-def') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read deformation direction (x, y or z)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
      & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') GOTO 120
    !read deformation
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    !if user entered the symbol "%" then remove it and convert
    m=SCAN(temp,"%")
    IF(m>0) THEN
      temp(m:m)=" "
      READ(temp,*,END=120,ERR=120) tempreal
      WRITE(temp,'(f16.3)') tempreal/100.d0
    ENDIF
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    READ(temp,*,END=120,ERR=120) tempreal
    !read Poisson ratio
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    READ(temp,*,END=120,ERR=120) tempreal
  !
  ELSEIF(clarg=='-dislocation' .OR. clarg=='-disloc') THEN
    m=0
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read first coordinate of disloc.
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF( TRIM(temp)=='loop' ) THEN
      !Construct a dislocation loop
      !read coordinates (x,y,z) of dislocation loop
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      !Read direction normal to loop plane
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
        & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') GOTO 120
      !Read dislocation loop radius
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      READ(temp,*,END=120,ERR=120) tempreal
      options_array(ioptions) = TRIM(options_array(ioptions))//" "//TRIM(temp)
      !Read the Burgers vector
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//" "//TRIM(temp)
      !If we just read Miller indices, finish here
      !Otherwise, we just read a real number equal to bx => read by and bz
      IF( SCAN(temp,'[')==0 .AND. SCAN(temp,']')==0 .AND. SCAN(temp,'_')==0 ) THEN
        !Check that bx was indeed a real number
        READ(temp,*,END=120,ERR=120) tempreal
        !Read the coordinates by and bz of Burgers vector
        i=i+1
        READ(cla(i),*,END=400,ERR=400) temp
        READ(temp,*,END=120,ERR=120) tempreal
        options_array(ioptions) = TRIM(options_array(ioptions))//" "//TRIM(temp)
        i=i+1
        READ(cla(i),*,END=400,ERR=400) temp
        READ(temp,*,END=120,ERR=120) tempreal
        options_array(ioptions) = TRIM(options_array(ioptions))//" "//TRIM(temp)
        !Read Poisson ratio
        i=i+1
        READ(cla(i),*,END=400,ERR=400) temp
        READ(temp,*,END=120,ERR=120) tempreal
        options_array(ioptions) = TRIM(options_array(ioptions))//" "//TRIM(temp)
      ENDIF
      !
    ELSE !i.e. if temp.NE.'loop'
      !Construct a straight dislocation line
      !=> temp should contain a real number corresponding to p1 = pos(1)
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      !READ(temp,*,END=120,ERR=120) tempreal
      !read second coordinate of disloc.
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
          IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
            & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      !READ(temp,*,END=120,ERR=120) tempreal
      !read the character of dislocation (edge, edge2, screw or mixed)
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      temp = TRIM(ADJUSTL(temp))
      temp2 = temp !important for reading the Burgers vector afterwards
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      IF( TRIM(temp).NE.'screw' .AND. TRIM(temp(1:4)).NE.'edge' &
        & .AND. TRIM(temp).NE.'mixed' ) GOTO 120
      IF( temp(1:4)=="edge" ) m=1
      !read the dislocation line direction (x, y or z)
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      temp = TRIM(ADJUSTL(temp))
      temp3 = temp !important for reading the Burgers vector afterwards
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
        & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') GOTO 120
      !read the normal to cut plane (x, y or z)
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      temp = TRIM(ADJUSTL(temp))
      temp4 = temp !important for reading the Burgers vector afterwards
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
        & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') GOTO 120
      !read the Burgers vector
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      READ(temp,*,END=120,ERR=120) tempreal
      IF(temp2=="mixed") THEN
        !three components are given: the first one was b(1)
        options_array(ioptions) = TRIM(options_array(ioptions))//" "//TRIM(temp)
        !Read b(2)
        i=i+1
        READ(cla(i),*,END=400,ERR=400) temp
        READ(temp,*,END=120,ERR=120) tempreal
        options_array(ioptions) = TRIM(options_array(ioptions))//" "//TRIM(temp)
        !Read b(3)
        i=i+1
        READ(cla(i),*,END=400,ERR=400) temp
        READ(temp,*,END=120,ERR=120) tempreal
        options_array(ioptions) = TRIM(options_array(ioptions))//" "//TRIM(temp)
      ELSE
        !Only one component is given
        options_array(ioptions) = TRIM(options_array(ioptions))//" "//TRIM(temp)
      ENDIF
      !read Poisson's ratio (only if it is of edge character)
      IF( m==1 ) THEN
        i=i+1
        READ(cla(i),*,END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        READ(temp,*,END=120,ERR=120) tempreal
      ENDIF
    ENDIF
    !scan the final option line
    !if slashes are present (user wants to perform a division), replace them by a colon (:)
    j = SCAN(options_array(ioptions),'/')
    DO WHILE(j>0)
      options_array(ioptions)(j:j) = ':'
      j = SCAN(options_array(ioptions),'/')
    ENDDO
  !
  ELSEIF(clarg=='-disturb') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read the maximum displacement along X
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    READ(temp,*,END=120,ERR=120) tempreal
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !try to read the maximum displacement along Y
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp2
    !Verify that it is a real number
    READ(temp2,*,END=101,ERR=101) tempreal
    !No error? It was the max. displacement along Y, save it
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp2)
    !Proceed with max. disp. along Z
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp2
    READ(temp2,*,END=101,ERR=101) tempreal
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp2)
    !We are done reading the 3 real numbers
    GOTO 110
    101 CONTINUE
    !There was an error while trying to real max. displacement along Y
    !=> user gave only one value, that will apply in all directions X, Y and Z
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)//' '//TRIM(temp)
    !Go one step back
    i=i-1
  !
  ELSEIF(clarg=='-duplicate' .OR. clarg=='-dup' .OR. clarg=='-replicate' ) THEN
    ioptions = ioptions+1
    options_array(ioptions) = '-duplicate'
    DO m=1,3
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      READ(temp,*,END=120,ERR=120) tempreal
    ENDDO
  !
  ELSEIF(clarg=='-fix' .OR. clarg=='-freeze') THEN
    ioptions = ioptions+1
    options_array(ioptions) = "-fix"
    !read fixaxis: must be x, y, z, or all
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    IF( temp(1:3)=="xyz" .OR. temp(1:3)=="XYZ" ) temp = "all"
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
      & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z' .AND.  &
      & temp(1:3).NE.'all' ) GOTO 120
    !read 'above' or 'below'
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    IF(temp(1:5)=='above' .OR. temp(1:5)=='below') THEN
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !read fix distance
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      !READ(temp,*,END=120,ERR=120) tempreal
      !read fix direction (x, y or z)
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      temp = TRIM(ADJUSTL(temp))
      IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
        & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z' .AND.  &
        & temp(1:3).NE.'all' ) GOTO 120
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    ELSE
      i=i-1
    ENDIF
  !
  ELSEIF(clarg=='-frac' .OR. clarg=='-fractional' .OR. clarg=='-reduce' .OR. clarg=='-reduced') THEN
    ioptions = ioptions+1
    options_array(ioptions) = '-fractional'
  !
  ELSEIF(clarg=='-mirror') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read distance between mirror plane and origin
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
      & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
    !read normal to mirror plane (x, y, z, or Miller indices)
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  ELSEIF(clarg=='-orient') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read the 1st vector of ancient base
    i=i+1
    READ(cla(i),'(a128)',END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read the 2nd vector of ancient base
    i=i+1
    READ(cla(i),'(a128)',END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read the 3rd vector of ancient base
    i=i+1
    READ(cla(i),'(a128)',END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read the 1st vector of new base
    i=i+1
    READ(cla(i),'(a128)',END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read the 2nd vector of new base
    i=i+1
    READ(cla(i),'(a128)',END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read the 3rd vector of new base
    i=i+1
    READ(cla(i),'(a128)',END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  ELSEIF(clarg=='-orthogonal-cell' .OR. clarg=='-orthorhombic-cell' .OR. &
        & clarg=='-orthogonal-box' .OR. clarg=='-orthorhombic-box'  .OR. &
        & clarg=='-orthocell' .OR. clarg=='-orthobox' ) THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
  !
  ELSEIF(clarg=='-prop'.OR.clarg=='-properties'.OR.clarg=="property") THEN
    ioptions = ioptions+1
    options_array(ioptions) = "-prop"
    !read the name of property file
    i=i+1
    READ(cla(i),'(a128)',END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  ELSEIF(clarg=='-rebox') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
  !
  ELSEIF(clarg=='-remove-atom'  .OR. clarg=='-rmatom'  .OR. &
        &clarg=='-remove-atoms' .OR. clarg=='-rmatoms'      ) THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read the atom species or index that will be removed
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  ELSEIF(clarg=='-remove-doubles' .OR. clarg=='-remove-double' .OR. clarg=='-rmd') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Read the max. distance
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    READ(temp,*,END=120,ERR=120) tempreal
  !
  ELSEIF(clarg=='-remove-property' .OR. clarg=='-rmprop') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read the name of the property to be removed
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  ELSEIF( clarg=='-remove-shells' .OR. clarg=='-remove-shell' .OR. &
        & clarg=='-rmshell' .OR. clarg=='-rmshells' ) THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read the atom species on which shells will be removed
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  ELSEIF(clarg=='-roll' .OR. clarg=="-bend") THEN
    ioptions = ioptions+1
    options_array(ioptions) = "-roll"
    !read the direction (x, y or z)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
      & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') THEN
      !No axis detected => error
      GOTO 120
    ENDIF
    !read the angle (degrees)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp2
    j=SCAN(temp2,"°")
    IF(j>0) THEN
      temp2 = temp2(1:j-1)
    ENDIF
    READ(temp2,*,END=120,ERR=120) tempreal
    !read the axis of rolling (x, y or z)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp3
    IF( temp3(1:1).NE.'x' .AND. temp3(1:1).NE.'y' .AND. temp3(1:1).NE.'z' .AND.  &
      & temp3(1:1).NE.'X' .AND. temp3(1:1).NE.'Y' .AND. temp3(1:1).NE.'Z') THEN
      !No axis detected => error
      GOTO 120
    ENDIF
    !temp and temp3 must *not* be the same direction
    IF( TRIM(ADJUSTL(temp))==TRIM(ADJUSTL(temp3)) ) THEN
      GOTO 120
    ENDIF
    options_array(ioptions) = TRIM(options_array(ioptions))//' '// &
                            & TRIM(temp)//' '//TRIM(temp2)//' '//TRIM(temp3)
  !
  ELSEIF(clarg=='-rotate' .OR. clarg=='-rot') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Check if first
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    IF( temp(1:3)=="com" ) THEN
      options_array(ioptions) = TRIM(options_array(ioptions))//' com '
      !read the axis of rotation (x, y or z) and the angle of rotation
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
    ENDIF
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp2
    temp = TRIM(ADJUSTL(temp))
    temp2 = TRIM(ADJUSTL(temp2))
    j=SCAN(temp,"°")
    IF(j>0) THEN
      temp = temp(1:j-1)
    ENDIF
    j=SCAN(temp2,"°")
    IF(j>0) THEN
      temp2 = temp2(1:j-1)
    ENDIF
    !User may have entered "-rotate axis angle" or "-rotate angle axis"
    !detect which is which and save "-rotate axis angle" in options_array
    IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
      & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') THEN
      IF( temp2(1:1).NE.'x' .AND. temp2(1:1).NE.'y' .AND. temp2(1:1).NE.'z' .AND.  &
      & temp2(1:1).NE.'X' .AND. temp2(1:1).NE.'Y' .AND. temp2(1:1).NE.'Z') THEN
        !No axis detected => error
        GOTO 120
      ELSE
        !Axis is in temp2
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp2)//' '//TRIM(temp)
        READ(temp,*,END=120,ERR=120) tempreal
      ENDIF
      !
    ELSE
      !Axis is in temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)//' '//TRIM(temp2)
      READ(temp2,*,END=120,ERR=120) tempreal
    ENDIF
  !
  ELSEIF(clarg=='-select') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read first keyword
    i=i+1
    temp = ADJUSTL(cla(i))
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !
    IF( temp=='all' .OR. temp=='any' .OR. temp=='invert' .OR. temp=='none' ) THEN
      !No other parameter to read
      CONTINUE
    ELSEIF( temp=='above' .OR. temp=='below' ) THEN
      !read distance normal to plane of selection
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !read axis normal to plane of selection: can be x, y, z, or a crystallographic vector
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      temp = TRIM(ADJUSTL(temp))
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    ELSEIF( temp=='in' .OR. temp=='out' ) THEN
      !A geometric region is defined
      !read the region geometry: 'cell', 'box', 'sphere', 'cylinder'
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !IF( temp.NE.'box' .AND. temp.NE.'sphere' .AND. temp.NE.'cylinder' .AND. &
      !  & temp.NE.'prism' .AND. temp.NE.'torus' ) GOTO 120
      region_geom = temp(1:16)
      !Next parameters depend on the geometry
      IF(region_geom=='cell') THEN
        !No other parameters
        CONTINUE
      ELSEIF(region_geom=='box') THEN
        !Read the coordinates of the first corner of the box
        !X
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !Y
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !Z
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !read the coordinates of the last corner of the box
        !X
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !Y
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !Z
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
      ELSEIF(region_geom=='sphere') THEN
        !Read the coordinates of the center of the sphere
        !X
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !Y
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !Z
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !read the radius of the sphere
        i=i+1
        READ(cla(i),*,END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        READ(temp,*,END=120,ERR=120) tempreal
      ELSEIF(region_geom=='cylinder') THEN
        !read the axis of the cylinder
        i=i+1
        READ(cla(i),*,END=400,ERR=400) temp
        temp = TRIM(ADJUSTL(temp))
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
          & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') GOTO 120
        !Read the 2 coordinates of the center of the cylinder
        !in the plane normal to its axis
        !X
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !Y
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !read the radius of the cylinder
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        READ(temp,*,END=120,ERR=120) tempreal
      ELSEIF(region_geom=='torus') THEN
        !read the axis (normal to the torus plane / to the base of pyramid)
        i=i+1
        READ(cla(i),*,END=400,ERR=400) temp
        temp = TRIM(ADJUSTL(temp))
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
          & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') GOTO 120
        !read the 3 coordinates of the center of the torus / of the base of pyramid
        !X
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !Y
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !Z
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !read main radius of torus / side of pyramid base
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 ) GOTO 120
        !read secondary radius of torus / height of pyramid
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 ) GOTO 120
      ELSE
        !Unrecognized shape => display error
        GOTO 120
      ENDIF
      !
    ELSEIF( temp=='grid' ) THEN
      !Read the name of the file containing grid elements
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !
    ELSEIF( temp=='list' ) THEN
      !Read the name of the file atom indices will be read from
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !
    ELSEIF( temp=='prop' .OR. temp=='property' ) THEN
      !Atoms must be selected according to a property
      !Read name of the property
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !Read criterion (will be interpreted later)
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !
    ELSEIF( temp=='stl' .OR. temp=='STL' ) THEN
      !Atoms inside a 3-D model must be selected
      !Read name of the STL file
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      temp = TRIM(ADJUSTL(temp))
      IF( temp(1:6)=="center" ) THEN
        i=i+1
        READ(cla(i),*,END=400,ERR=400) temp
        options_array(ioptions) = '-select stl center '//TRIM(temp)
      ELSE
        options_array(ioptions) = '-select stl '//TRIM(temp)
      ENDIF
      !
    ELSEIF( temp=='random' .OR. temp=='rand' ) THEN
      !A number N of atoms of given species must be selected at random
      !read number of atoms N
      !read the atom species or index that will be removed
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !read atom species; can be "all" or "any"
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      IF(LEN_TRIM(temp)>3) GOTO 120
      !
    ELSEIF( i+2<=SIZE(cla) .AND.                                                           &
          &  ( cla(i+2)=="neighbors" .OR. cla(i+2)=="neighbours" .OR. cla(i+2)=="neigh" )  &
          & ) THEN
      options_array(ioptions) = '-select neigh '
      !Read number of neighbors, or cutoff radius for neighbor search
      READ(cla(i),'(a)',END=400,ERR=400) temp
      READ(temp,*,END=120,ERR=120) tempreal  !check that it is a number
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !Read species of neighbor(s) (can be "all" or "any")
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !skip keyword "neighbor"
      i=i+1
      !Read index of atom whose neighbors must be found
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      READ(temp,*,END=120,ERR=120) tempreal  !check that it is a number
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !
    ELSEIF( LEN_TRIM(temp)<=2 .AND. temp(1:2).NE.'in' ) THEN
      !It should be an atom species: try to recognize it
      CALL ATOMNUMBER(temp(1:2),tempreal)
      IF( NINT(tempreal)==0 ) THEN
        !it is not an atom species, then it can be an integer,
        !or a range of integer, or a list of integers separated by a comma
        !This will be dealt with inside the option (see "opt_select.f90")
        !READ(temp,*,ERR=120,END=120) m
      ENDIF
      !
    ENDIF
    !if slashes are present (user wants to perform a division), replace them by a colon (:)
    j = SCAN(options_array(ioptions),'/')
    DO WHILE(j>0)
      options_array(ioptions)(j:j) = ':'
      j = SCAN(options_array(ioptions),'/')
    ENDDO
  !
  ELSEIF(clarg=='-separate' .OR. clarg=='-sep') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Read the max. distance
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    READ(temp,*,END=120,ERR=120) tempreal
    !Read the separation distance
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    READ(temp,*,END=120,ERR=120) tempreal
  !
  ELSEIF(clarg=='-shear') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read the normal to the sheared surface (x, y or z)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
      & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') GOTO 120
    !read the tilt vector (in A) or shear strain (in %)
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !make sure it contains a real number
    m = SCAN(temp,'%')
    IF( m>0 ) THEN
      temp2 = temp(1:m-1)
      READ(temp2,'(a)',END=120,ERR=120) tempreal
    ELSE
      m = SCAN(temp,'A')
      IF(m>0) temp(m:m)=" "
      READ(temp,'(a)',END=120,ERR=120) tempreal
    ENDIF
    !read the shear direction (x, y or z)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
      & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') GOTO 120
  !
  ELSEIF(clarg=='-shift') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read part to shift ('above' or 'below', or a number)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    IF( temp(1:5)=='above' .OR. temp(1:5)=='below' ) THEN
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !read distance of plane of shift
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !read axis perpendicular to the plane of shift (x, y, z, or crystallographic direction)
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      temp = TRIM(ADJUSTL(temp))
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    ELSE
      i=i-1
    ENDIF
    !read tau1
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
      & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read tau2
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
      & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read tau3
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
      & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !if slashes are present (user wants to perform a division), replace them by a colon (:)
    j = SCAN(options_array(ioptions),'/')
    DO WHILE(j>0)
      options_array(ioptions)(j:j) = ':'
      j = SCAN(options_array(ioptions),'/')
    ENDDO
  !
  ELSEIF(clarg=='-sort') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Read which property must be sorted
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !Read how to sort it
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF( temp.NE.'up' .AND. temp.NE.'down' .AND. temp.NE.'pack') GOTO 120
  !
  ELSEIF(clarg=='-spacegroup' .OR. clarg=='-space-group' .OR. clarg=='-sgroup' &
        &.OR. clarg=='-sg') THEN
    ioptions = ioptions+1
    options_array(ioptions) = '-spacegroup'
    !Read space group number or name
    i=i+1
    READ(cla(i),'(a128)',END=400,ERR=400) temp
    temp = ADJUSTL(temp)
    !If first letter is lower case, make it upper case
    temp(1:1) = StrUpCase(temp(1:1))
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  ELSEIF(clarg=='-stress') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Read stress component
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = ADJUSTL(temp)
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    SELECT CASE(temp)
    CASE('x','X','xx','XX','y','Y','yy','YY','z','Z','zz','ZZ', &
        & 'xy','XY','yx','YX','zx','ZX','xz','XZ','zy','ZY','yz','YZ','p','P')
      !Read value of stress
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      READ(temp,'(a)',END=120,ERR=120) tempreal
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    END SELECT
  !
  ELSEIF(clarg=='-substitute' .OR. clarg=='-sub') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF(LEN_TRIM(temp)>3) GOTO 120
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF(LEN_TRIM(temp)>3) GOTO 120
  !
  ELSEIF(clarg=='-swap') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Read first value to swap (must be X,Y,Z or an integer)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    m = 0
    IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
      & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') THEN
      !It *must* be an integer
      READ(temp,*,END=400,ERR=400) m
    ENDIF
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !Read second value to swap
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    IF( m.NE.0 ) THEN
      !the second value must also be an integer
      READ(temp,*,END=400,ERR=400) m
    ELSE
      !the second value must be X,Y or Z
      IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
        & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') THEN
        GOTO 120
      ENDIF
    ENDIF
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  ELSEIF(clarg=='-torsion') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read the axis of torsion (x, y or z) and the angle of rotation
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp2
    temp = TRIM(ADJUSTL(temp))
    temp2 = TRIM(ADJUSTL(temp2))
    j=SCAN(temp,"°")
    IF(j>0) THEN
      temp = temp(1:j-1)
    ENDIF
    j=SCAN(temp2,"°")
    IF(j>0) THEN
      temp2 = temp2(1:j-1)
    ENDIF
    !User may have entered "-torsion axis angle" or "-torsion angle axis"
    !detect which is which and save "-torsion axis angle" in options_array
    IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
      & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') THEN
      IF( temp2(1:1).NE.'x' .AND. temp2(1:1).NE.'y' .AND. temp2(1:1).NE.'z' .AND.  &
      & temp2(1:1).NE.'X' .AND. temp2(1:1).NE.'Y' .AND. temp2(1:1).NE.'Z') THEN
        !No axis detected => error
        GOTO 120
      ELSE
        !Axis is in temp2
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp2)//' '//TRIM(temp)
        READ(temp,*,END=120,ERR=120) tempreal
      ENDIF
      !
    ELSE
      !Axis is in temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)//' '//TRIM(temp2)
      READ(temp2,*,END=120,ERR=120) tempreal
    ENDIF
  !
  ELSEIF(clarg=='-unit' .OR. clarg=='-units' .OR. clarg=='-u') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = cla(i)
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = cla(i)
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  ELSEIF(clarg=='-unskew') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
  !
  ELSEIF(clarg=='-velocity') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Read target temperature for Maxwell-Boltzmann distribution
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    READ(temp,'(a)',END=120,ERR=120) tempreal
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  ELSEIF(clarg=='-wrap') THEN
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
  !
  ! -- add other options in alphabetical order --
  !
  ELSEIF(clarg(1:2)=='--') THEN
    !if it starts with "--" we assume it is a wrong mode entered by the user
    nerr=nerr+1
    CALL ATOMSK_MSG(4813,(/TRIM(clarg)/),(/0.d0/))
    GOTO 1000
  !
  ELSEIF(clarg(1:1)=='-') THEN
    !if it starts with "-" we assume it is a wrong option entered by the user
    nerr=nerr+1
    CALL ATOMSK_MSG(2805,(/TRIM(clarg)/),(/0.d0/))
    GOTO 1000
  !
  !If it is none of the above, we assume it is some file name
  ELSE !IF( LEN_TRIM(pfiles(1))==0 .OR. LEN_TRIM(pfiles(2))==0 ) THEN
    IF(pfiles(1)=='') THEN
      WRITE(pfiles(1),*) TRIM(clarg)
    ELSEIF(pfiles(2)=='') THEN
      WRITE(pfiles(2),*) TRIM(clarg)
    ELSE
      !More than two file names => display warning
      nwarn=nwarn+1
      CALL ATOMSK_MSG(704,(/TRIM(clarg)/),(/0.d0/))
    ENDIF
  !
  !
  !And the rest we do not understand => display a warning
!   ELSE
!     nwarn = nwarn+1
!     CALL ATOMSK_MSG(703,(/TRIM(clarg)/),(/0.d0/))
  !
  ENDIF
  !
  110 CONTINUE
  IF(i>SIZE(cla)) GOTO 120
ENDDO
!
!Save output formats to array outfileformats
IF(Nout>0) THEN
  ALLOCATE(outfileformats(Nout))
  DO i=1,Nout
    outfileformats(i) = tempout(i)
  ENDDO
ENDIF
GOTO 1000
!
120 CONTINUE
READ(options_array(ioptions),*) temp
CALL ATOMSK_MSG(2806,(/TRIM(temp)/),(/0.d0/))
CALL DISPLAY_HELP(temp)
GOTO 400
!
130 CONTINUE
temp = TRIM(ADJUSTL(mode))
CALL ATOMSK_MSG(4800,(/TRIM(temp)/),(/0.d0/))
CALL DISPLAY_HELP(temp)
GOTO 400
!
150 CONTINUE
CALL ATOMSK_MSG(4814,(/''/),(/0.d0/))
GOTO 400
!
!
!
400 CONTINUE
nerr=nerr+1
!
!
!
1000 CONTINUE
!Check size of array options_array(:)
!if it is empty or no option was detected in the command-line parameters, de-allocate it
IF( ALLOCATED(options_array) ) THEN
  IF( ioptions<=0 ) THEN
    DEALLOCATE(options_array)
  ELSEIF( LEN_TRIM(options_array(1))<1 ) THEN
    DEALLOCATE(options_array)
  ELSEIF( SIZE(options_array)<1 ) THEN
    DEALLOCATE(options_array)
  ENDIF
ENDIF
!
!
END SUBROUTINE GET_CLA
!
END MODULE read_cla
