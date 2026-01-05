MODULE read_cla
!
!**********************************************************************************
!*  READ_CLA                                                                      *
!**********************************************************************************
!* This module reads command-line arguments for Atomsk.                           *
!**********************************************************************************
!* (C) March 2011 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 15 Dec. 2025                                     *
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
USE strings
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
CHARACTER(LEN=16):: region_geom  !geometry of the region: "box" or "sphere" or "cylinder" or... etc.
CHARACTER(LEN=16):: select_mul   !keyword for multiple selection ("add" or "rm" or "intersect" etc.)
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
INTEGER:: i, ioptions, j, k, m
INTEGER:: Nout !number of output formats
REAL(dp):: smass, tempreal
!
!
!Initialize variables
 clarg=""
msg=""
temp=""
temp2=""
temp3=""
temp4=""
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
  clarg = TRIM(ADJUSTL(cla(i)))
  IF(LEN_TRIM(clarg)==0) EXIT
  !
  !
  !Deal with modes
  SELECT CASE(StrDnCase(clarg))
  CASE("--gather","--all-in-one","--ai1")
    mode = "gather"
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(1)
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(2)
    IF(i>SIZE(cla)) GOTO 130
    !
  CASE("--average","--avg")
    mode = "average"
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(5)
    IF(LEN_TRIM(pfiles(5))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  CASE( "--local-symmetry","--localsymmetry","--central-symmetry","--centralsymmetry", &
      & "--centro-symmetry","--centrosymmetry","--symmetry","--sym","--cs","--ls" )
    mode = "cs"
    i=i+1
    READ(cla(i),'(a4096)',END=130,ERR=130) pfiles(3) !File containing atom positions
    IF(LEN_TRIM(pfiles(3))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  CASE("--copy-properties","--cpprop")
    mode = "cpprop"
    i=i+1
    READ(cla(i),'(a)',END=130,ERR=130) pfiles(3)
    i=i+1
    READ(cla(i),'(a)',END=130,ERR=130) pfiles(4)
    IF(LEN_TRIM(pfiles(3))==0 .OR. LEN_TRIM(pfiles(4))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  CASE("--create")
    IF(ALLOCATED(mode_param)) GOTO 150
    ALLOCATE(mode_param(20))
    mode_param(:) = ''
    mode="create"
    !Get the type of structure
    i=i+1
    m=1
    READ(cla(i),*,END=130,ERR=130) temp
    mode_param(m) = TRIM(ADJUSTL(temp))
    !Get the lattice constant a
    i=i+1
    m=m+1
    READ(cla(i),*,END=130,ERR=130) mode_param(m)
    !Detect if a second lattice constant (c) or chiral indices must be read
    SELECT CASE(StrDnCase(temp))
    CASE("graphite","hcp","wurtzite","wz","c14","c36","l10","l1_0","limo2", &
        & "st","bct","fct")
      !Get lattice constant c
      i=i+1
      m=m+1
      READ(cla(i),*,END=130,ERR=130) mode_param(m)
    CASE("nanotube","nt")
      !Get chiral indices (m,n) of the nanotube
      i=i+1
      m=m+1
      READ(cla(i),*,END=130,ERR=130) mode_param(m)
      i=i+1
      m=m+1
      READ(cla(i),*,END=130,ERR=130) mode_param(m)
    END SELECT
    !Get the atomic species for the structure
    DO WHILE( i+1<=SIZE(cla) )
      i=i+1
      temp=""
      IF(i<=SIZE(cla)) temp = ADJUSTL(cla(i))
      IF( LEN_TRIM(temp)>0 ) THEN
!         IF( SCAN(temp,'0123456789')>0 ) THEN
!           !It is a number: user probably entered too many lattice constants
!           GOTO 130
!         ELSE
        IF(LEN_TRIM(temp)>2 .OR. temp(1:1)=='-' ) THEN
          !It is not an atom species => go back and exit the loop
          i=i-1
          EXIT
        ELSE
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
        ENDIF
      ELSE
        !Empty entry => exit
        EXIT
      ENDIF
    ENDDO
    !Look for crystallographic orientation (if any)
    IF( i+1 <= SIZE(cla)-3 ) THEN
      i=i+1
      temp = ADJUSTL(cla(i))
      IF( temp=="orient" ) THEN
        m=m+1
        mode_param(m) = "orient"
        !Get Miller indices [hkl] or [hkil] of the crystal orientation along X, Y, Z
        i=i+1
        m=m+1
        READ(cla(i),*,END=130,ERR=130) mode_param(m)
        i=i+1
        m=m+1
        READ(cla(i),*,END=130,ERR=130) mode_param(m)
        i=i+1
        m=m+1
        READ(cla(i),*,END=130,ERR=130) mode_param(m)
      ELSE
        i=i-1
      ENDIF
    ENDIF
    IF(LEN_TRIM(mode_param(1))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  CASE("--density")
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
  CASE("--difference","--diff")
    mode = 'diff'
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(1)
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(2)
    IF(LEN_TRIM(pfiles(1))==0 .OR. LEN_TRIM(pfiles(2))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  CASE("--ddplot")
    mode = "ddplot"
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(3)
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(4)
    IF(LEN_TRIM(pfiles(3))==0 .OR. LEN_TRIM(pfiles(4))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  CASE("--edm","--electric-dipole-moments")
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
  CASE("--electronic-polarization","--pe")
    mode = 'PE'
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(1) !File containing atom positions
    IF(LEN_TRIM(pfiles(1))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  CASE("--interpolate")
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
  CASE("--list","--l")
    mode = 'list'
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(5)
    IF(LEN_TRIM(pfiles(5))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  CASE("--match-id","--matchid")
    mode = "match-id"
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(3)
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(4)
    IF(LEN_TRIM(pfiles(3))==0 .OR. LEN_TRIM(pfiles(4))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  CASE("--merge","--m","--stack")
    IF(ALLOCATED(mode_param)) GOTO 150
    ALLOCATE(mode_param(100))  !assuming user won't merge more than 100 files...
    mode_param(:) = ''
    mode = "merge"
    !Store all parameters and file names in mode_param(:)
    j=0
    i=i+1
    m=0
    DO WHILE(m<=0)
      READ(cla(i),*,ERR=130,END=130) temp
      IF( StrDnCase(temp)=="stack" ) THEN
        j = j+1
        mode_param(j) = TRIM(ADJUSTL(temp))  !"stack"
        !Read direction of stacking
        i=i+1
        READ(cla(i),*,ERR=130,END=130) temp
        IF( StrUpCase(temp)=='X' .OR. StrUpCase(temp)=='Y' .OR. StrUpCase(temp)=='Z' ) THEN
          j = j+1
          mode_param(j) = TRIM(ADJUSTL(temp))
        ENDIF
        i=i+1
      ELSEIF( StrDnCase(temp)=="rescale" .OR. StrDnCase(temp)=="scale" .OR. &
            & StrDnCase(temp)=="match" ) THEN
        j = j+1
        mode_param(j) = TRIM(ADJUSTL(temp))  !"rescale"
        !Read direction along which systems will be rescaled to match 1st system's size
        i=i+1
        READ(cla(i),*,ERR=130,END=130) temp
        !IF( StrUpCase(temp)=='X' .OR. StrUpCase(temp)=='Y' .OR. StrUpCase(temp)=='Z' ) THEN
          j = j+1
          mode_param(j) = TRIM(ADJUSTL(temp))
        !ENDIF
        i=i+1
      ELSEIF( StrDnCase(temp)=="x" .OR. StrDnCase(temp)=="y" .OR. StrDnCase(temp)=="z" ) THEN
        !User directly gave direction without keyword: assume it is the stacking direction
        j = j+1
        mode_param(j) = "stack "//TRIM(ADJUSTL(temp))
      ELSE
        !If none of the above, must be the number of file names to follow
        READ(temp,*,ERR=130,END=130) m
        IF(m<=0) GOTO 130
        j=j+1
        WRITE(mode_param(j),*) m
        EXIT
      ENDIF
    ENDDO
    !Read the m file names that follow
    DO k=1,m
      i=i+1
      READ(cla(i),'(a128)',ERR=130,END=130) temp
      j=j+1
      mode_param(j) = TRIM(ADJUSTL(temp))
    ENDDO
    !Read name of the output file
    i=i+1
    READ(cla(i),'(a128)',ERR=130,END=130) pfiles(1)
    IF(LEN_TRIM(pfiles(1))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  CASE("--normal","convert")
    !nothing special to do, this is normal mode
    !
  CASE("--nye")
    mode = "nye"
    i=i+1
    READ(cla(i),'(a)',END=130,ERR=130) pfiles(3)
    i=i+1
    READ(cla(i),'(a)',END=130,ERR=130) pfiles(4)
    IF(LEN_TRIM(pfiles(3))==0 .OR. LEN_TRIM(pfiles(4))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  CASE("--unfold","--one-in-all","-1ia")
    mode = "unfold"
    i=i+1
    READ(cla(i),'(a128)',ERR=130,END=130) pfiles(1)
    IF(LEN_TRIM(pfiles(1))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  CASE("--polycrystal","--polyx")
    mode = "polycrystal"
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(1) !File containing atom positions
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(3) !File containing parameters for Voronoi construction
    IF(LEN_TRIM(pfiles(1))==0 .OR. i>SIZE(cla)) GOTO 130
    !
  CASE("--rdf","--radial-distribution-function")
    IF(ALLOCATED(mode_param)) GOTO 150
    ALLOCATE(mode_param(2))
    mode = "rdf"
    i=i+1
    READ(cla(i),'(a4096)',END=130,ERR=130) pfiles(5) !File list
    IF(LEN_TRIM(pfiles(5))==0 .OR. i>SIZE(cla)) GOTO 130
    i=i+1
    READ(cla(i),*,END=130,ERR=130) mode_param(1)  !R
    i=i+1
    READ(cla(i),*,END=130,ERR=130) mode_param(2)  !dR
    IF( i>SIZE(cla) ) GOTO 130
    !
  CASE("--unwrap")
    IF(ALLOCATED(mode_param)) GOTO 150
    mode = 'unwrap'
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(3) !"Reference"
    i=i+1
    READ(cla(i),'(a128)',END=130,ERR=130) pfiles(4) !"Configuration" to unwrap
    IF(LEN_TRIM(pfiles(3))==0 .OR. LEN_TRIM(pfiles(4))==0 .OR. i>SIZE(cla)) GOTO 130
  !
  !
  !
  !
  !Deal with options
  !For each option we store the option parameters in the string array "options_array"
  CASE("-add-atom","-add-atoms","-addatom","-addatoms")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Read the species of atom to add
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    IF( LEN_TRIM(temp)<=2 ) THEN
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    ELSE
      GOTO 120
    ENDIF
    !Read keyword: must be "at", "near" or "random"
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(StrDnCase(temp)))
    IF( .NOT.( temp(1:2)=="at" .OR. temp(1:1)=="@" .OR. temp(1:4)=="near" .OR. &
      &        temp(1:6)=="random" .OR. temp(1:3)=="rel"  ) ) THEN
      GOTO 400
    ENDIF
    temp2 = temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !Read parameters
    SELECT CASE(StrDnCase(temp2))
    CASE("at","@")
      !Read coordinate x of atom to add
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !Read coordinate y of atom to add
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !Read coordinate z of atom to add
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
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
      IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !Read coordinate y of atom to add
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !Read coordinate z of atom to add
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !
    CASE("near")
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
  CASE("-add-shells","-addshells")
    ioptions = ioptions+1
    options_array(ioptions) = '-add-shells'
    !read atom species for which the shells must be created
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = ADJUSTL(temp)
    IF( LEN_TRIM(temp)<=2 .OR. temp=="all" ) THEN
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    ELSE
      GOTO 120
    ENDIF
  !
  CASE("-alignx")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
  !
  CASE("-bind-shells","-bs")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
  !
  CASE("-crack")
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
      IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
    IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
      & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !Read second coordinate of crack tip (pos2)
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
      IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
    IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
      & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !Read direction along crack tip (X, Y or Z)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
      IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
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
  CASE("-cell","-box","-change-cell","-change_cell","-change-box","-change_box")
    ioptions = ioptions+1
    options_array(ioptions) = "-cell"
    !Read direction
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(StrDnCase(temp)))
    SELECT CASE(temp)
    CASE( "h1","h2","h3",'x','y','z',"xx","yy","zz", &
        & "xy","xz","yx","yz","zx","zy","xyz","all"  )
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    CASE DEFAULT
      GOTO 120
    END SELECT
    !Read operation to perform (add, rm, set)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(StrDnCase(temp)))
    IF( temp.NE."add" .AND. temp.NE."rm" .AND. temp.NE."set" ) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !Read value (real number)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    IF( SCAN(temp,'0123456789')==0 ) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  CASE("-center")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    IF(StrDnCase(temp)=="com") WRITE(temp,'(a1)') "0"
    READ(temp,*,END=120,ERR=120) j
    IF(j<0) WRITE(temp,'(a1)') "0"
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  CASE("-cut")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read 'above' or 'below'
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(StrDnCase(temp)))
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF( temp(1:5).NE.'above' .AND. temp(1:5).NE.'below' ) GOTO 120
    !read cut distance
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
    IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
      & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read cut direction
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  CASE("-deform","-def")
    j=0
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read deformation direction (x, y or z)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
      & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') GOTO 120
    IF( StrDnCase(temp(1:2))=="x " .OR. StrDnCase(temp(1:2))=="y " .OR. StrDnCase(temp(1:2))=="z " ) &
      & j=1
    !read deformation
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    IF( j==0 .AND. (StrDnCase(temp(1:7))=="untilt " .OR. StrDnCase(temp(1:8))=="unshear ") ) THEN
      !Tilt
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    ELSE
      !if user entered the symbol "%" then remove it and convert
      m=SCAN(temp,"%")
      IF(m>0) THEN
        temp(m:m)=" "
        READ(temp,*,END=120,ERR=120) tempreal
        WRITE(temp,'(f16.6)') tempreal/100.d0
      ELSE
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
      ENDIF
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      READ(temp,*,END=120,ERR=120) tempreal
    ENDIF
    IF( j==1 ) THEN
      !attempt reading Poisson ratio
      i=i+1
      READ(cla(i),'(a)',END=110,ERR=110) temp
      READ(temp,*,END=102,ERR=102) tempreal
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      GOTO 110
      102 CONTINUE
      i=i-1
    ENDIF
  !
  CASE("-denoise")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
      READ(temp,*,END=120,ERR=120) tempreal
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  CASE("-dislocation","-disloc","-dislo")
    m=0
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read first coordinate of disloc.
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF( TRIM(temp)=='loop' ) THEN
      !Construct a dislocation loop
      !read coordinates (x,y,z) of dislocation loop
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
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
      IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
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
    ELSEIF( TRIM(temp)=="array" .OR. TRIM(temp)=="file" ) THEN
      !Last possiblity: read dislocation coordinates from a file
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//" "//TRIM(temp)
      !Read Poisson ratio
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      READ(temp,*,END=120,ERR=120) tempreal
      options_array(ioptions) = TRIM(options_array(ioptions))//" "//TRIM(temp)
      !
    ELSE
      !This should be parameters for constructing a straight dislocation line
      !=> temp should contain a real number corresponding to p1 = pos(1)
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      !READ(temp,*,END=120,ERR=120) tempreal
      !read second coordinate of disloc.
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      !READ(temp,*,END=120,ERR=120) tempreal
      !read the character of dislocation: edge, edge_add, edge_rm, screw, or mixed
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      temp = TRIM(ADJUSTL(temp))
      temp2 = temp !important for reading the Burgers vector afterwards
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      IF( TRIM(temp).NE.'screw' .AND. TRIM(temp).NE.'edge' .AND. TRIM(temp).NE.'edge_add' &
        & .AND. TRIM(temp).NE.'edge-add' .AND. TRIM(temp).NE.'edge_rm' .AND. TRIM(temp).NE.'edge-rm' &
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
      IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
      READ(temp,*,END=120,ERR=120) tempreal
      IF(temp2=="mixed") THEN
        !three components are given: the first one was b(1)
        options_array(ioptions) = TRIM(options_array(ioptions))//" "//TRIM(temp)
        !Read b(2)
        i=i+1
        READ(cla(i),*,END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
        READ(temp,*,END=120,ERR=120) tempreal
        options_array(ioptions) = TRIM(options_array(ioptions))//" "//TRIM(temp)
        !Read b(3)
        i=i+1
        READ(cla(i),*,END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
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
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
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
  CASE("-disturb")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read the maximum displacement along X
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    READ(temp,*,END=120,ERR=120) tempreal
    tempreal = DABS(tempreal)
    WRITE(temp,*) tempreal
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(ADJUSTL(temp))
    !try to read the maximum displacement along Y
    i=i+1
    IF( i<=SIZE(cla) ) THEN
      READ(cla(i),*,END=101,ERR=101) temp2
      !Verify that it is a real number
      READ(temp2,*,END=101,ERR=101) tempreal
      !No error? It was the max. displacement along Y, save it
      tempreal = DABS(tempreal)
      WRITE(temp2,*) tempreal
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(ADJUSTL(temp2))
      !Proceed with max. disp. along Z
      !At this point user provided 2 values, there MUST be a third, otherwise it is an error
      i=i+1
      IF(i>SIZE(cla)) GOTO 120
      READ(cla(i),*,END=120,ERR=120) temp2
      READ(temp2,*,END=120,ERR=120) tempreal
      tempreal = DABS(tempreal)
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(ADJUSTL(temp2))
      !We are done reading the 3 real numbers
      GOTO 110
    ENDIF
    101 CONTINUE
    !There was an error while trying to real max. displacement along Y
    !=> user gave only one value which is the maximum displacement
    !Go one step back
    i=i-1
  !
  CASE("-duplicate","-dup","-replicate")
    ioptions = ioptions+1
    options_array(ioptions) = "-duplicate"
    DO m=1,3
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      READ(temp,*,END=120,ERR=120) tempreal
    ENDDO
  !
  CASE("-freeze","-fix")
    IF(clarg=="-fix") THEN
      CALL ATOMSK_MSG(2799,(/"-fix   ","-freeze"/),(/0.d0/))
    ENDIF
    ioptions = ioptions+1
    options_array(ioptions) = "-freeze"
    !read coordinate(s) to freeze: must be x, y, z, xy, xz, yz, xyz, or all
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    SELECT CASE(StrDnCase(temp))
    CASE('x','y','z',"xy","xz","yx","yz","zx","zy","xyz","xzy","yxz","yzx","zxy","zyx","all")
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    CASE DEFAULT
      GOTO 120
    END SELECT
    !read 'above' or 'below'
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    IF(temp(1:5)=='above' .OR. temp(1:5)=='below') THEN
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !read distance to plane
      i=i+1
      READ(cla(i),'(a)',END=400,ERR=400) temp
      IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
        & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
      !READ(temp,*,END=120,ERR=120) tempreal
      !read direction normal to plane (x, y or z)
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
  CASE("-frac","-fractional","-reduce","-reduced")
    ioptions = ioptions+1
    options_array(ioptions) = '-fractional'
  !
  CASE("-mirror")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read distance between mirror plane and origin
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
      & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
    !read normal to mirror plane (x, y, z, or Miller indices)
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  CASE("-orient")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read the 1st vector of ancient base
    i=i+1
    IF(i>SIZE(cla)) GOTO 120
    READ(cla(i),'(a128)',END=400,ERR=400) temp
    IF( temp(1:1)=='-' .OR. SCAN(temp,'0123456789')==0 ) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read the 2nd vector of ancient base
    i=i+1
    IF(i>SIZE(cla)) GOTO 120
    READ(cla(i),'(a128)',END=400,ERR=400) temp
    IF( temp(1:1)=='-' .OR. SCAN(temp,'0123456789')==0 ) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read the 3rd vector of ancient base
    i=i+1
    IF(i>SIZE(cla)) GOTO 120
    READ(cla(i),'(a128)',END=400,ERR=400) temp
    IF( temp(1:1)=='-' .OR. SCAN(temp,'0123456789')==0 ) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read the 1st vector of new base
    i=i+1
    IF(i>SIZE(cla)) GOTO 120
    READ(cla(i),'(a128)',END=400,ERR=400) temp
    IF( temp(1:1)=='-' .OR. SCAN(temp,'0123456789')==0 ) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read the 2nd vector of new base
    i=i+1
    IF(i>SIZE(cla)) GOTO 120
    READ(cla(i),'(a128)',END=400,ERR=400) temp
    IF( temp(1:1)=='-' .OR. SCAN(temp,'0123456789')==0 ) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read the 3rd vector of new base
    i=i+1
    IF(i>SIZE(cla)) GOTO 120
    READ(cla(i),'(a128)',END=400,ERR=400) temp
    IF( temp(1:1)=='-' .OR. SCAN(temp,'0123456789')==0 ) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  CASE("-orthogonal-cell","-orthorhombic-cell","-orthogonal-box","-orthorhombic-box", &
      & "-orthocell","-orthobox")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
  !
  CASE("-prop","-properties","property")
    ioptions = ioptions+1
    options_array(ioptions) = "-prop"
    !read the name of property file
    i=i+1
    READ(cla(i),'(a128)',END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  CASE("-rebox","-shrink-wrap","-shrinkwrap")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
  !
  CASE("-reduce-cell","-reducecell","-reduce-box","-reducebox")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Read direction along which the cell must be reduced (optional)
    IF( i<SIZE(cla) ) THEN
      i=i+1
      READ(cla(i),'(a128)',END=400,ERR=400) temp
      IF( temp=="x".OR.temp=="X" .OR. temp=="y".OR.temp=="Y" .OR. temp=="z".OR.temp=="Z" &
        &  .OR. temp=="p".OR.temp=="P") THEN
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      ELSE
        i=i-1
      ENDIF
    ENDIF
  !
  CASE("-remove-atom","-rmatom","-remove-atoms","-rmatoms","-delete-atoms","-delete_atoms")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read the atom species or index that will be removed
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  CASE("-remove-doubles","-remove-double","-rmd")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Read the max. distance
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    READ(temp,*,END=120,ERR=120) tempreal
  !
  CASE("-remove-property","-rmprop")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read the name of the property to be removed
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  CASE("-remove-shells","-remove-shell","-rmshell","-rmshells")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read the atom species on which shells will be removed
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    IF( LEN_TRIM(temp)<=2 .OR. temp=="all" ) THEN
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    ELSE
      GOTO 120
    ENDIF
  !
  CASE("-roll","-bend")
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
  CASE("-rotate","-rot")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Check if first
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    IF( temp(1:3)=="com" ) THEN
      options_array(ioptions) = TRIM(options_array(ioptions))//' com '
      !read the axis of rotation (x, y, z, or Miller vector, or 3 real numbers)
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      temp = TRIM(ADJUSTL(temp))
    ENDIF
    !
    IF( temp=='x' .OR. temp=='X' .OR. temp=='y' .OR. temp=='Y' .OR. temp=='z' .OR.temp=='Z') THEN
      !Read the rotation angle
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp2
      !Verify that it is a real number
      j=SCAN(temp2,'°')
      IF(j>0) temp2=temp2(:j-1)
      READ(temp2,*,END=120,ERR=120) tempreal
      !Save it into options_array
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)//' '//TRIM(temp2)
    ELSEIF( (SCAN(temp,'[')>0 .OR. SCAN(temp,']')>0 .OR. SCAN(temp,'_')>0) .AND. SCAN(temp,'.')==0 ) THEN
      !It should be a Miller vector
      !Read the rotation angle
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp2
      !Verify that it is a real number
      j=SCAN(temp2,'°')
      IF(j>0) temp2=temp2(:j-1)
      READ(temp2,*,END=120,ERR=120) tempreal
      !Save it into options_array
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)//' '//TRIM(temp2)
    ELSEIF( SCAN(StrDnCase(temp),'abcdefghijklmnopqrstuvwxyz')>0 ) THEN
      !It is probably the name of a file containing a rotation matrix
      !Save it into options_array
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    ELSE
      !There should be 3 real numbers, the first one is already in temp
      !Verify that it is a real number
      READ(temp,*,END=120,ERR=120) tempreal
      !Read the second real number
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp2
      !Verify that it is a real number
      READ(temp2,*,END=120,ERR=120) tempreal
      !Read the third real number
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp3
      !Verify that it is a real number
      READ(temp3,*,END=120,ERR=120) tempreal
      !Read the rotation angle
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp4
      !Verify that it is a real number
      j=SCAN(temp4,'°')
      IF(j>0) temp4=temp4(:j-1)
      READ(temp4,*,END=120,ERR=120) tempreal
      !Save it into options_array
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)//' '//&
                              TRIM(temp2)//' '//TRIM(temp3)//' '//TRIM(temp4)
    ENDIF
  !
  CASE("-roundoff","-round-off")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read property that will be rounded off
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read threshold value
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  CASE("-select")
    select_mul = ""
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !read first keyword
    i=i+1
    temp = ADJUSTL(StrDnCase(cla(i)))
    !
    IF( temp=="add" .OR. temp=="union" .OR. temp=="rm" .OR. temp=="remove" .OR.     &
      & temp=="subtract" .OR. temp=="intersect" .OR. temp=="xor" .OR. temp=="among" ) THEN
      select_mul = TRIM(ADJUSTL(temp))
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      i=i+1
      temp = ADJUSTL(StrDnCase(cla(i)))
    ENDIF
    !
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
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
      region_geom = StrDnCase(temp(1:16))
      !Next parameters depend on the geometry
      IF(region_geom=='cell') THEN
        !No other parameters
        CONTINUE
      ELSEIF(region_geom=='box' .OR. region_geom=='block') THEN
        !Read the coordinates of the first corner of the box
        !X
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !Y
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !Z
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !read the coordinates of the last corner of the box
        !X
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !Y
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !Z
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
      ELSEIF(region_geom=='sphere') THEN
        !Read the coordinates of the center of the sphere
        !X
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !Y
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !Z
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !read the radius of the sphere
        i=i+1
        READ(cla(i),*,END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
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
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !READ(temp,*,END=120,ERR=120) tempreal
        !Y
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
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
      ELSEIF(region_geom=='cone') THEN
        !read the axis of the cone
        i=i+1
        READ(cla(i),*,END=400,ERR=400) temp
        temp = TRIM(ADJUSTL(temp))
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( temp(1:1).NE.'x' .AND. temp(1:1).NE.'y' .AND. temp(1:1).NE.'z' .AND.  &
          & temp(1:1).NE.'X' .AND. temp(1:1).NE.'Y' .AND. temp(1:1).NE.'Z') GOTO 120
        !read the 3 coordinates of the tip of the cone
        !X
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !Y
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !Z
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !read opening angle of the cone (in degrees)
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 ) GOTO 120
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
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !Y
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
        options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
        IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
          & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
        !Z
        i=i+1
        READ(cla(i),'(a)',END=400,ERR=400) temp
        IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
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
      temp2 = ""
      !Detect keywords "center", "scale" or "rescale", "fill", concatenate them into one string
      DO WHILE( StrDnCase(temp(1:6))=="center" .OR. StrDnCase(temp(1:7))=="rescale" .OR.  &
              & StrDnCase(temp(1:5))=="scale"  .OR. StrDnCase(temp(1:4))=="fill"          )
        temp2 = TRIM(temp2)//TRIM(ADJUSTL(temp))
        i=i+1
        READ(cla(i),*,END=400,ERR=400) temp
        temp = TRIM(ADJUSTL(temp))
      ENDDO
      !Copy the keywords
      options_array(ioptions) = TRIM(ADJUSTL(options_array(ioptions)))//" "//TRIM(temp2)//" "//TRIM(temp)
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
    ELSEIF( i+1<=SIZE(cla) .AND. ( cla(i+1)=="modulo" .OR. cla(i+1)=="mod" ) ) THEN
      options_array(ioptions) = "-select "//TRIM(select_mul)//" modulo "
      !Read index
      READ(cla(i),'(a)',END=400,ERR=400) temp
      READ(temp,*,END=120,ERR=120) tempreal  !check that it is a number
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      !Read modulo
      i=i+2
      READ(cla(i),'(a)',END=400,ERR=400) temp
      READ(temp,*,END=120,ERR=120) tempreal  !check that it is a number
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      IF( NINT(tempreal)==0 ) THEN
        !Display error message and abort (division by zero)
        CALL ATOMSK_MSG(2821,(/""/),(/0.d0/))
        GOTO 120
      ENDIF
      !
    ELSEIF( i+2<=SIZE(cla) .AND.                                                           &
          &  ( cla(i+2)=="neighbors" .OR. cla(i+2)=="neighbours" .OR. cla(i+2)=="neigh"    &
          &    .OR. cla(i+2)=="neighbor" .OR. cla(i+2)=="neighbour" )             ) THEN
      options_array(ioptions) = "-select "//TRIM(select_mul)//" neigh "
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
        !it is not an atom species, then it can be an integer
        !This will be dealt with inside the option (see "opt_select.f90")
        !READ(temp,*,ERR=120,END=120) m
      ENDIF
      !
    ELSE
      !It may be an integer, or a range of integer, or a list of integers separated by a comma
      !e.g. 3,6,13 or 4:12 => verify that
      DO j=1,LEN_TRIM(temp)
        IF( SCAN(temp(j:j),' 0123456789:,') == 0 ) THEN
          !Illegal character in this string
          GOTO 120
        ENDIF
      ENDDO
      !
    ENDIF
    !if slashes are present (user wants to perform a division), replace them by a colon (:)
    j = SCAN(options_array(ioptions),'/')
    DO WHILE(j>0)
      options_array(ioptions)(j:j) = ':'
      j = SCAN(options_array(ioptions),'/')
    ENDDO
  !
  CASE("-separate","-sep")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Read the max. distance
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    READ(temp,*,END=120,ERR=120) tempreal
    !Read the translation distance/vector
    !Read first real number
    i=i+1
    READ(cla(i),*,END=110,ERR=110) temp
    temp = ADJUSTL(temp)
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !Try to read a second real number
!     i=i+1
!     READ(cla(i),*,END=110,ERR=110,IOSTAT=j) temp
!     READ(temp,*,END=110,ERR=110,IOSTAT=j) tempreal
!     IF( j==0 ) THEN
!       !Succeeded reading 2nd number: there MUST be a third
!       i=i+1
!       READ(cla(i),*,END=400,ERR=400) temp
!       READ(temp,*,END=120,ERR=120) tempreal
!     ELSE
!       i=i-1
!     ENDIF
  !
  CASE("-shift","-move","-translate","-translation")
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
    IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
    IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
      & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read tau2
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
    IF( SCAN(temp,'0123456789')==0 .AND. INDEX(temp,'INF')==0 .AND. &
      & INDEX(temp,'box')==0 .AND. INDEX(temp,'BOX')==0) GOTO 120
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !read tau3
    i=i+1
    READ(cla(i),'(a)',END=400,ERR=400) temp
    IF( SCAN(temp,"*")>0 ) temp = "'"//TRIM(ADJUSTL(temp))//"'"
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
  CASE("-sort")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Read which property must be sorted
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    IF( temp.NE."random" .AND. temp.NE."reverse" ) THEN
      !Read how to sort it
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
      SELECT CASE(temp)
      CASE("up","UP","down","DOWN","pack","PACK")
        !These keywords are valid
        temp = StrDnCase(temp)
        CONTINUE
      CASE DEFAULT
        !Invalid keyword: exit with error
        GOTO 120
      END SELECT
    ENDIF
  !
  CASE("-spacegroup","-space-group","-sgroup","-sg")
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
  CASE("-stress")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Read stress component
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = ADJUSTL(temp)
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    SELECT CASE(StrDnCase(temp))
    CASE('x','xx','y','yy','z','zz','xy','yx','zx','xz','zy','yz','p')
      !Read value of stress
      i=i+1
      READ(cla(i),*,END=400,ERR=400) temp
      READ(temp,'(a)',END=120,ERR=120) tempreal
      options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    END SELECT
  !
  CASE("-substitute","-sub")
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
  CASE("-swap")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Read first value to swap
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    temp = TRIM(ADJUSTL(temp))
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
    !Read second value to swap
    i=i+1
    READ(cla(i),*,END=120,ERR=120) temp
    temp = TRIM(ADJUSTL(temp))
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  CASE("-torsion")
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
  CASE("-unit","-units","-u")
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
  CASE("-unskew")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
  !
  CASE("-velocity")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
    !Read target temperature for Maxwell-Boltzmann distribution
    i=i+1
    READ(cla(i),*,END=400,ERR=400) temp
    READ(temp,'(a)',END=120,ERR=120) tempreal
    options_array(ioptions) = TRIM(options_array(ioptions))//' '//TRIM(temp)
  !
  CASE("-wrap")
    ioptions = ioptions+1
    options_array(ioptions) = TRIM(clarg)
  !
  ! -- add other options in alphabetical order --
  !
  CASE DEFAULT
    !Deal with output formats
    !First, check if argument is one of the formats in flist(:) (see "globalvar.f90")
    IF( ANY(flist(:,1)==clarg) ) THEN
      Nout = Nout+1
      tempout(Nout) = StrDnCase(clarg)
    !Second, accept some keywords and save them as format
    ELSEIF(clarg=="atomsk") THEN
      Nout = Nout+1
      tempout(Nout) = 'atsk'
    ELSEIF( clarg=='abinit' ) THEN
      Nout = Nout+1
      tempout(Nout) = 'abin'
    ELSEIF( clarg=="bopfox" ) THEN
      Nout = Nout+1
      tempout(Nout) = 'bx'
    ELSEIF( clarg=='atomeye' .OR. clarg=='qstem' ) THEN
      Nout = Nout+1
      tempout(Nout) = 'cfg'
    ELSEIF( clarg=='drprobe' ) THEN
      Nout = Nout+1
      tempout(Nout) = 'cel'
    ELSEIF( clarg=='dlpoly' ) THEN
      Nout = Nout+1
      tempout(Nout) = 'dlp'
    ELSEIF( clarg=='config' ) THEN
      IF(LEN_TRIM(pfiles(1)).NE.0) THEN
        IF(LEN_TRIM(pfiles(2))==0) pfiles(2) = 'CONFIG'
        Nout = Nout+1
        tempout(Nout) = 'dlp'
      ELSE
        pfiles(1) = 'CONFIG'
      ENDIF
    ELSEIF( clarg=='mbpp' ) THEN
      Nout = Nout+1
      tempout(Nout) = 'coo'
    ELSEIF( clarg=='coorat' ) THEN
      !if a file name was defined before,
      IF(LEN_TRIM(pfiles(1)).NE.0) THEN
        IF(LEN_TRIM(pfiles(2))==0) pfiles(2) = 'COORAT'
        Nout = Nout+1
        tempout(Nout) = 'coo'
      ELSE
        pfiles(1) = 'COORAT'
      ENDIF
    ELSEIF( clarg=='crystal' ) THEN
      Nout = Nout+1
      tempout(Nout) = 'd12'
    ELSEIF( clarg=='ddplot' ) THEN
      Nout = Nout+1
      tempout(Nout) = 'dd'
    ELSEIF( clarg=='siesta' ) THEN
      Nout = Nout+1
      tempout(Nout) = 'fdf'
    ELSEIF( clarg=='gulp' ) THEN
      Nout = Nout+1
      tempout(Nout) = 'gin'
    ELSEIF( clarg=='lammps' ) THEN
      Nout = Nout+1
      tempout(Nout) = 'lmp'
    ELSEIF( clarg=='moldy' ) THEN
      Nout = Nout+1
      tempout(Nout) = 'mol'
    ELSEIF( clarg=='vasp' ) THEN
      Nout = Nout+1
      tempout(Nout) = 'pos'
    ELSEIF( clarg=='poscar' ) THEN
      IF(LEN_TRIM(pfiles(1)).NE.0) THEN
        IF(LEN_TRIM(pfiles(2))==0) pfiles(2) = 'POSCAR'
        Nout = Nout+1
        tempout(Nout) = 'pos'
      ELSE
        IF(LEN_TRIM(pfiles(2))==0) pfiles(1) = 'POSCAR'
      ENDIF
    ELSEIF( clarg=='pwscf' .OR. clarg=='qe' ) THEN
      Nout = Nout+1
      tempout(Nout) = 'pw'
    ELSEIF( clarg=='xcrysden' ) THEN
      Nout = Nout+1
      tempout(Nout) = 'xsf'
    ELSEIF( clarg=='stru' ) THEN
      Nout = Nout+1
      tempout(Nout) = 'stru'
    ! -- please add new formats in alphabetical order --
    !
    ELSEIF(clarg(1:2)=='--') THEN
      !if it starts with "--" we assume it is a wrong mode entered by the user
      nerr=nerr+1
      CALL ATOMSK_MSG(4813,(/TRIM(clarg)/),(/0.d0/))
      GOTO 1000
    !
    ELSEIF(clarg=='-') THEN
      !output to stdout
      ofu=6
      verbosity=0  !disable all other messages
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
    ENDIF
  !
  !
  !And the rest we do not understand => display a warning
!   ELSE
!     nwarn = nwarn+1
!     CALL ATOMSK_MSG(703,(/TRIM(clarg)/),(/0.d0/))
  !
  END SELECT  !clarg
  !
  GOTO 110
  !
  109 CONTINUE
  i=i-1
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
IF( verbosity>=4 ) THEN
  IF( ALLOCATED(options_array) ) THEN
    WRITE(msg,*) "SIZE OF ARRAY options: ", SIZE(options_array)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    DO i=1,SIZE(options_array)
      WRITE(msg,*) TRIM(options_array(i))
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    ENDDO
  ELSE
    msg = "ARRAY options NOT ALLOCATED"
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  ENDIF
ENDIF
!
!
END SUBROUTINE GET_CLA
!
END MODULE read_cla
