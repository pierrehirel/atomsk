MODULE options
!
!**********************************************************************************
!*  OPTIONS                                                                       *
!**********************************************************************************
!* This module applies one or several transformations to the system.              *
!*                                                                                *
!* WHAT MUST BE PROVIDED:                                                         *
!* options_array                                                                  *
!*           if allocated, string array of unknown dimension containing the       *
!*           names and parameters of options to be applied to the system.         *
!*           If not allocated, then no option will be applied.                    *
!* H         3 x 3 real array containing vectors of the supercell. H(1,:) is      *
!*           the first vector, H(2,:) the second, H(3,:) the third.               *
!* P         N x 4 real array containing atom positions. N=number of atoms,       *
!*           each line contains the coordinates (X,Y,Z) and atomic number.        *
!* S         if allocated, N x 4 real array containing positions of shells        *
!*           (for ionic core/shell model), and number of shells N must be         *
!*           identical to N in P. If S is not allocated, then no shell exist.     *
!* AUXNAMES  if allocated, string array of dimension M, M being the number of     *
!*           existing auxiliary properties (must be equal to M in AUX). If        *
!*           AUXNAMES is not allocated, then AUX must also be non-allocated.      *
!* AUX       if allocated, N x M real array containing values of auxiliary        *
!*           properties of atoms, e.g. their velocity, forces, etc. The first     *
!*           dimension N must be equal to the number of atoms N in P, the second  *
!*           dimension M must be equal to the number of auxiliary properties.     *
!*           If AUX is not allocated then no auxiliary property exists, and       *
!*           AUXNAMES must also be non-allocated.                                 *
!*                                                                                *
!* WHAT IS RETURNED BY THIS ROUTINE:                                              *
!* modified arrays H, P, S, AUXNAMES and AUX.                                     *
!*                                                                                *
!**********************************************************************************
!* (C) October 2010 - Pierre Hirel                                                *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 18 Sept. 2014                                    *
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
USE atoms
USE comv
USE constants
USE functions
USE messages
USE files
USE subroutines
USE guess_form
USE deterH
!
!Modules managing options
USE addatom
USE addshells
USE alignx
USE bindshells
USE carttofrac
USE center
USE crack
USE cut_cell
USE deform
USE dislocation
USE disturb
USE expand
USE fix
USE mirror
USE orient
USE properties
USE rmatom
USE remdoubles
USE rmprop
USE rmshells
USE rotate
USE select
USE shear
USE shift
USE sort
USE substitute
USE unit
USE unskew
USE velocity
USE wrap
! -- add other options in alphabetical order --
!
!
!
CONTAINS
!
SUBROUTINE OPTIONS_AFF(options_array,H,P,S,AUXNAMES,AUX,ORIENT,SELECT)
!
!Declare variables
IMPLICIT NONE
CHARACTER(LEN=2):: species
CHARACTER(LEN=32):: optionname
CHARACTER(LEN=128):: temp, msg
CHARACTER(LEN=128),DIMENSION(10):: treal !text containing a real number and maybe a word
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
INTEGER:: i, ioptions, j
INTEGER:: status
INTEGER:: strlength
REAL(dp):: tempreal !temporary real number
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: HS  !Copy of H for shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P  !atomic positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S  !shell positions (is any)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX, AUXdummy !auxiliary properties of atoms/shells
!
!Variables relative to Option: add-atom
CHARACTER(LEN=2):: addatom_species !species of atom(s) to add
CHARACTER(LEN=8):: addatom_type    !"at" or "near" or "random"
REAL(dp),DIMENSION(3):: addatom_prop  !properties of atom(s) to add
                                      !if addatom_type=="at", position x,y,z of new atom
                                      !if addatom_type=="near", index of atom
                                      !if addatom-type=="random", number of atoms to add
!
!Variables relative to Option: crack
CHARACTER(LEN=1):: crackline !direction of crack tip line, must be x, y or z
CHARACTER(LEN=1):: crackplane !normal to the plane of cut, must be x, y or z
CHARACTER(LEN=3):: crackmode !crack mode, must be I, II or III
CHARACTER(LEN=11):: cracktype !planestress or planestrain
REAL(dp):: crackK  !stress intensity factor
REAL(dp):: mu   !shear modulus of the material
!
!Variables for option: center
INTEGER:: center_atom  !index of atom that must be centered; if 0, center the center of mass
!
!Variables for option: create shells
CHARACTER(LEN=2):: cs_species
!
!Variables relative to Option: cut
CHARACTER(LEN=5):: cut_dir   !above or below
CHARACTER(LEN=16):: cutdir   !x, y, z, or crystallographic direction
REAL(dp):: cutdistance
!
!Variables relative to Option: deform
CHARACTER(LEN=1):: def_dir       !direction of applied strain (X, Y or Z)
REAL(dp):: def_strain, def_poisson  !applied strain and Poisson's ratio
!
!Variables relative to Option: dislocation
CHARACTER(LEN=1):: dislocline     !(x, y or z)
CHARACTER(LEN=1):: dislocplane    !(x, y or z)
CHARACTER(LEN=5):: disloctype     !edge or screw
REAL(dp):: nu, pos1, pos2
REAL(dp),DIMENSION(3):: b !Burgers vector
!
!Variables relative to Option: disturb
REAL(dp):: dist_dxmax  !maximum translation vector
!
!Variables relative to Option: expand
INTEGER, DIMENSION(3):: expandmatrix  !number of times the system must be
                                      !repeated in each direction of space
!
!Variables relative to Option: fix
CHARACTER(LEN=5):: fix_dir, fixaxis
CHARACTER(LEN=16):: fixdir   !x, y, z, or crystallographic direction
REAL(dp):: fixdistance
!
!Variables relative to Option: mirror
CHARACTER(LEN=16):: mirror_dir   !x, y, z, or crystallographic direction
REAL(dp):: mirror_d  !distance between mirror plane and origin of coordinates
!
!Variables relative to Option: orient
CHARACTER(LEN=16):: oldvec1, oldvec2, oldvec3, newvec1, newvec2, newvec3
REAL(dp),DIMENSION(3,3):: Hstart, Hend
!
!Variables relative to Option: properties
CHARACTER(LEN=128):: propfile  !file containing the properties
REAL(dp),DIMENSION(3):: lat_a0   !lattice constants a, b, c
REAL(dp),DIMENSION(3,3):: ORIENT  !crystalographic orientation
REAL(dp),DIMENSION(9,9):: C_tensor, C_tensor_dummy !elastic tensor (and dummy)
!
!Variables relative to Option: rmatom
CHARACTER(LEN=32):: rmatom_prop  !what atom(s) to erase
!
!Variables relative to Option: remove doubles
REAL(dp):: rmd_radius            !in angstroms
!
!Variables relative to Option: rmprop
CHARACTER(LEN=32):: rmprop_prop  !property that must be erased
!
!Variables relative to Option: rmshells
CHARACTER(LEN=6):: rmshells_prop  !species on which shells are removed
!
!Variables relative to Option: rotate
CHARACTER(LEN=1):: rot_axis     !(x, y or z)
REAL(dp):: rot_angle              !in degrees
!
!Variables relative to Option: select
CHARACTER(LEN=16):: region_dir   !x, y, z
CHARACTER(LEN=8):: region_side   !'in' or 'out' or 'all'
CHARACTER(LEN=16):: region_geom  !geometry of the region: "box" or "sphere"
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
REAL(dp),DIMENSION(3):: region_1 !First corner for'box', or center of sphere
REAL(dp),DIMENSION(3):: region_2 !Last corner for'box', or radius of sphere
!
!Variables relative to Option: shear
CHARACTER(LEN=1):: shear_dir  !direction of applied shear strain (X, Y or Z)
CHARACTER(LEN=1):: shear_surf !normal to the surface to be sheared (X, Y or Z)
REAL(dp):: shear_strain  !applied shear strain in %
!
!Variables relative to Option: shift
CHARACTER(LEN=5):: shift_dir       !'above' or 'below'
CHARACTER(LEN=16):: shift_axis     !x, y, z or crystallographic direction
REAL(dp):: shift_dist       !distance of the "cut plane" to the origin
REAL(dp):: shift_tau1, shift_tau2, shift_tau3  !shift vector
!
!Variables relative to Option: sort
CHARACTER(LEN=4):: sortorder
CHARACTER(LEN=16):: sortcol
!
!Variables relative to Option: substitute
CHARACTER(LEN=3):: sp1, sp2
!
!Variables relative to Option: unit
CHARACTER(LEN=16):: unit1, unit2
!
!Variables relative to Option: velocity
REAL(dp):: vel_T  !target temperature for Maxwell-Boltzmann distribution
!
! -- add other options variables in alphabetical order --
!
j=0
region_1(:) = 0.d0
region_2(:) = 0.d0
!
!
msg = 'ENTERING OPTIONS'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!Initialize variables
lat_a0(:) = 1.d0
!If orientation of the system is unknown, default orientation is assumed to be [100] [010] [001]
IF( NINT(VECLENGTH(ORIENT(1,:))) == 0 .OR.          &
  & NINT(VECLENGTH(ORIENT(2,:))) == 0 .OR.          &
  & NINT(VECLENGTH(ORIENT(3,:))) == 0       ) THEN
  ORIENT(:,:) = 0.d0
  DO i=1,3
    ORIENT(i,i) = 1.d0
  ENDDO
ENDIF
!By default the elastic tensor is not set
 C_tensor(:,:) = 0.d0
 C_tensor_dummy(:,:) = 0.d0
!
!If there is no option then exit this module
IF(.NOT.ALLOCATED(options_array)) GOTO 1000
IF( verbosity==4 ) THEN
  WRITE(msg,*) 'options_array size: ', SIZE(options_array)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  IF( SIZE(options_array)>0 ) THEN
    DO i=1,SIZE(options_array)
      WRITE(msg,*) 'options_array(', i, ') : ', TRIM(ADJUSTL(options_array(i)))
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
  ENDIF
ENDIF
IF(SIZE(options_array)<=0) GOTO 1000
!
IF( .NOT.ALLOCATED(P) .OR. SIZE(P)<=0 ) THEN
  !The user wants to apply an option, but no system exists in memory
  !=> error
  nerr=nerr+1
  CALL ATOMSK_MSG(2814,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
100 CONTINUE
!Loop on all options
!NOTE: options must be applied to shells and auxiliary properties at the same time
!     as they are applied to atoms! All modules that change the number
!     of atoms or change their order in the array P must also add, remove or
!     re-order the shells in S and/or auxiliary properties in AUX accordingly.
!
DO ioptions=1,SIZE(options_array)
  msg = 'option: '//TRIM(options_array(ioptions))
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  READ(options_array(ioptions),*,END=1000,ERR=1000) optionname
  optionname = TRIM(ADJUSTL(optionname))
  IF(optionname=='') EXIT
  HS=H
  !Make sure that the dummy arrays are empty or deallocated
  IF(ALLOCATED(AUXdummy)) DEALLOCATE(AUXdummy)
  C_tensor_dummy(:,:) = 0.d0
  !
  SELECT CASE(optionname)
  !
  CASE('')
    !empty entry => just ignore it
    CONTINUE
    !
  CASE('-add-atom','add-atoms','addatom','addatoms')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, addatom_species, &
        & addatom_type
    IF( addatom_type=="at" .OR. addatom_type=="AT" .OR. addatom_type=="@" ) THEN
      !Read coordinates x, y, z of atom to add
      !Can be real numbers or fractional coordinates with the keyword "box"
      READ(options_array(ioptions),*,END=800,ERR=800) optionname, addatom_species, &
        & addatom_type, treal(1), treal(2), treal(3)
      CALL BOX2DBLE( H(:,1) , treal(1) , addatom_prop(1) , status )
      IF(status>0) THEN
        temp = treal(1)
        GOTO 810
      ENDIF
      CALL BOX2DBLE( H(:,2) , treal(2) , addatom_prop(2) , status )
      IF(status>0) THEN
        temp = treal(2)
        GOTO 810
      ENDIF
      CALL BOX2DBLE( H(:,3) , treal(3) , addatom_prop(3) , status )
      IF(status>0) THEN
        temp = treal(3)
        GOTO 810
      ENDIF
    ELSEIF( addatom_type=="near" .OR. addatom_type=="NEAR" .OR.      &
          & addatom_type=="random" .OR. addatom_type=="RANDOM" .OR.  &
          & addatom_type=="rand" .OR. addatom_type=="RAND"           ) THEN
      !Read index of atom near which the new atom must be added,
      !or the number of atoms to insert randomly
      READ(options_array(ioptions),*,END=800,ERR=800) optionname, addatom_species, &
        & addatom_type, addatom_prop(1)
    ELSE
      nerr=nerr+1
      GOTO 1000
    ENDIF
    CALL ADDATOM_XYZ(H,P,S,AUX,addatom_species,addatom_type,addatom_prop,SELECT)
  !
  CASE('-add-shells','-addshells','-as')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, temp
    IF(temp=='all') THEN
      cs_species = 'XX'
    ELSE
      READ(temp,*,ERR=800,END=800) cs_species
    ENDIF
    CALL ADDSHELLS_XYZ(P,S,cs_species)
  !
  CASE('-alignx')
    CALL ALIGN_X(H,P,S,C_tensor)
  !
  CASE('-bind-shells','-bs')
    CALL BSHELLS_XYZ(H,P,S)
  !
  CASE('-center')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, center_atom
    CALL CENTER_XYZ(H,P,S,center_atom,SELECT)
  !
  CASE('-crack')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, crackmode, &
        & cracktype, crackK, treal(1), treal(2), crackline, crackplane, mu, nu
    CALL BOX2DBLE( H(:,1) , treal(1) , pos1 , status )
    IF(status>0) THEN
      temp = treal(1)
      GOTO 810
    ENDIF
    CALL BOX2DBLE( H(:,2) , treal(2) , pos2 , status )
    IF(status>0) THEN
      temp = treal(2)
      GOTO 810
    ENDIF
    CALL CRACK_XYZ(H,P,S,crackmode,cracktype,crackK,crackline,crackplane,mu,nu,pos1,pos2,SELECT,AUXNAMES,AUX)
  !
  CASE('-cut')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
        & cut_dir, treal(1), cutdir
    !Check if numbers contain a keyword like "BOX" or "INF"
    i=1
    SELECT CASE(cutdir)
    CASE('x','X')
      i=1
    CASE('y','Y')
      i=2
    CASE('z','Z')
      i=3
    END SELECT
    CALL BOX2DBLE( H(:,i) , treal(1) , cutdistance , status )
    IF(status>0) THEN
      temp = treal(1)
      GOTO 810
    ENDIF
    CALL CUTCELL(P,S,AUX,cut_dir,cutdistance,cutdir,ORIENT,SELECT)
  !
  CASE('-def', '-deform')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
        & def_dir, treal(1), def_poisson
    !treal may contain the shear strain expressed in percent, e.g. "3%"
    SELECT CASE(def_dir)
    CASE('x','X')
      i=1
    CASE('y','Y')
      i=2
    CASE('z','Z')
      i=3
    END SELECT
    CALL BOX2DBLE( H(:,i) , treal(1) , def_strain , status )
    IF(status>0) THEN
      temp = treal(1)
      GOTO 810
    ENDIF
    CALL DEFORM_XYZ(H,P,S,def_dir,def_strain,def_poisson,SELECT)
  !
  CASE('-disloc', '-dislocation')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
        & treal(9), treal(10), disloctype
    IF( disloctype=="mixed" ) THEN
      READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
          & treal(9), treal(10), disloctype, dislocline, dislocplane, b(1), b(2), b(3), nu
    ELSE
      READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
          & treal(9), treal(10), disloctype, dislocline, dislocplane, tempreal, nu
      !Determine the complete Burgers vector b(:)
      b(:) = 0.d0
      IF(disloctype=="screw") THEN
        !only one component is given, the one along dislocline
        SELECT CASE(dislocline)
        CASE("x","X")
          b(1) = tempreal
        CASE("y","Y")
          b(2) = tempreal
        CASE("z","Z")
          b(3) = tempreal
        END SELECT
      ELSEIF(disloctype=="edge") THEN
        !only one component is given, normal to dislocline and plane of cut
        SELECT CASE(dislocline)
        CASE("x","X")
          SELECT CASE(dislocplane)
          CASE('z','Z')
            b(2) = tempreal
          CASE('y','Y')
            b(3) = tempreal
          END SELECT
        CASE("y","Y")
          SELECT CASE(dislocplane)
          CASE('x','X')
            b(3) = tempreal
          CASE('z','Z')
            b(1) = tempreal
          END SELECT
        CASE("z","Z")
          SELECT CASE(dislocplane)
          CASE('x','X')
            b(2) = tempreal
          CASE('y','Y')
            b(1) = tempreal
          END SELECT
        END SELECT
      ENDIF
    ENDIF
    !Check if numbers contain a keyword like "BOX" or "INF"
    SELECT CASE(dislocplane)
    CASE('x','X')
      j=1
    CASE('y','Y')
      j=2
    CASE('z','Z')
      j=3
    END SELECT
    SELECT CASE(dislocline)
    CASE('x','X')
      treal(2) = treal(9)
      treal(3) = treal(10)
      IF(j==2) THEN
        i=3
      ELSEIF(j==3) THEN
        i=2
      ENDIF
    CASE('y','Y')
      treal(3) = treal(9)
      treal(1) = treal(10)
      IF(j==1) THEN
        i=3
      ELSEIF(j==3) THEN
        i=1
      ENDIF
    CASE('z','Z')
      treal(1) = treal(9)
      treal(2) = treal(10)
      IF(j==1) THEN
        i=2
      ELSEIF(j==2) THEN
        i=1
      ENDIF
    END SELECT
    CALL BOX2DBLE( H(:,i) , treal(i) , pos1 , status )
    IF(status>0) THEN
      temp = treal(i)
      GOTO 810
    ENDIF
    CALL BOX2DBLE( H(:,j) , treal(j) , pos2 , status )
    IF(status>0) THEN
      temp = treal(j)
      GOTO 810
    ENDIF
    CALL DISLOC_XYZ(H,P,S,disloctype,dislocline,dislocplane,b,nu,pos1,pos2,SELECT,AUXNAMES,AUX,C_tensor)
  !
  CASE('-disturb')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, dist_dxmax
    CALL DISTURB_XYZ(dist_dxmax,P,SELECT)
  !
  CASE('-e', '-expand', '-duplicate', '-dup')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
        & expandmatrix(1), expandmatrix(2), expandmatrix(3)
    !If expansion is zero along a direction, correct it
    DO i=1,3
      IF(expandmatrix(i)<=0) expandmatrix(i)=1
    ENDDO
    CALL EXPANDCELL(H,P,S,expandmatrix,SELECT,AUX)
  !
  CASE('-fix', '-freeze')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
        & fixaxis, fix_dir, treal(1), fixdir
    !Check if numbers contain a keyword like "BOX" or "INF"
    SELECT CASE(fixdir)
    CASE('x','X')
      i=1
    CASE('y','Y')
      i=2
    CASE('z','Z')
      i=3
    END SELECT
    CALL BOX2DBLE( H(:,i) , treal(1) , fixdistance , status )
    IF(status>0) THEN
      temp = treal(1)
      GOTO 810
    ENDIF
    CALL FIX_XYZ(P,AUXNAMES,AUX,fixaxis,fix_dir,fixdistance,fixdir,ORIENT,SELECT)
  !
  CASE('-frac','-fractional')
    CALL CART2FRAC_XYZ(H,P,S)
  !
  CASE('-mirror')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
        & treal(1), mirror_dir
    !Check if numbers contain a keyword like "BOX" or "INF"
    i=1
    SELECT CASE(mirror_dir)
    CASE('x','X')
      i=1
    CASE('y','Y')
      i=2
    CASE('z','Z')
      i=3
    END SELECT
    CALL BOX2DBLE( H(:,i) , treal(1) , mirror_d , status )
    IF(status>0) THEN
      temp = treal(1)
      GOTO 810
    ENDIF
    CALL MIRROR_XYZ(P,S,AUXNAMES,AUX,mirror_dir,mirror_d,ORIENT,SELECT)
  !
  CASE('-orient')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, oldvec1, oldvec2,  &
        & oldvec3, newvec1, newvec2, newvec3
    CALL INDEX_MILLER(oldvec1,Hstart(1,:),j)
    IF(j>0) GOTO 800
    CALL INDEX_MILLER(oldvec2,Hstart(2,:),j)
    IF(j>0) GOTO 800
    CALL INDEX_MILLER(oldvec3,Hstart(3,:),j)
    IF(j>0) GOTO 800
    CALL INDEX_MILLER(newvec1,Hend(1,:),j)
    IF(j>0) GOTO 800
    CALL INDEX_MILLER(newvec2,Hend(2,:),j)
    IF(j>0) GOTO 800
    CALL INDEX_MILLER(newvec3,Hend(3,:),j)
    IF(j>0) GOTO 800
    CALL ORIENT_XYZ(H,P,S,Hstart,Hend,SELECT,C_tensor)
  !
  CASE('-prop','-properties','-property')
    READ(options_array(ioptions),'(a128)',END=800,ERR=800) temp
    READ(temp(6:),'(a128)') propfile
    propfile = TRIM(ADJUSTL(propfile))
    CALL READ_PROPERTIES(propfile,H,P,ORIENT,C_tensor,AUXNAMES,AUX)
  !
  CASE('-rebox')
    CALL DETERMINE_H(H,P)
  !
  CASE('-remove-atom', '-rmatom')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, rmatom_prop
    CALL RMATOM_XYZ(P,S,AUX,rmatom_prop,SELECT)
  !
  CASE('-remove-doubles', '-rmd')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
        & rmd_radius
    CALL REMDOUBLES_XYZ(H,P,S,AUX,rmd_radius,SELECT)
  !
  CASE('-remove-property', '-rmprop')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, rmprop_prop
    CALL RMPROP_XYZ(AUXNAMES,AUX,rmprop_prop,SELECT)
  !
  CASE('-remove-shells','-remove-shell','-rmshells','-rmshell')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, rmshells_prop
    rmshells_prop = ADJUSTL(rmshells_prop)
    CALL RMSHELLS_XYZ(P,S,rmshells_prop,SELECT)
  !
  CASE('-rot', '-rotate')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
        & rot_axis, rot_angle
    CALL ROTATE_XYZ(H,P,S,AUXNAMES,AUX,rot_axis,rot_angle,SELECT,C_tensor)
  !
  CASE('-select')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, region_side
    IF( region_side=="above" .OR. region_side=="below" ) THEN
      READ(options_array(ioptions),*,END=800,ERR=800) optionname, region_side, treal(1), region_dir
      i=1
      SELECT CASE(region_dir)
      CASE('x','X')
        i=1
      CASE('y','Y')
        i=2
      CASE('z','Z')
        i=3
      CASE DEFAULT
        !nerr=nerr+1
        !GOTO 1000
      END SELECT
      CALL BOX2DBLE( H(:,i) , treal(1) , region_1(1) , status )
      IF(status>0) THEN
        temp = treal(1)
        GOTO 810
      ENDIF
    ELSEIF( region_side=="in" .OR. region_side=="out" ) THEN
      READ(options_array(ioptions),*,END=800,ERR=800) optionname, region_side, region_geom
      IF(region_geom=='box') THEN
        READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
            & region_side, region_geom, treal(1), treal(2), treal(3), &
            & treal(4), treal(5), treal(6)
        !Check if numbers contain a keyword like "BOX" or "INF"
        DO i=1,3
          CALL BOX2DBLE( H(:,i) , treal(i) , region_1(i) , status )
          IF(status>0) THEN
            temp = treal(i)
            GOTO 810
          ENDIF
        ENDDO
        DO i=1,3
          CALL BOX2DBLE( H(:,i) , treal(i+3) , region_2(i) , status )
          IF(status>0) THEN
            temp = treal(i+3)
            GOTO 810
          ENDIF
        ENDDO
      ELSEIF(region_geom=='sphere') THEN
        READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
            & region_side, region_geom, treal(1), treal(2), treal(3), treal(4)
        !Check if numbers contain a keyword like "BOX" or "INF"
        DO i=1,3
          CALL BOX2DBLE( H(:,i) , treal(i) , region_1(i) , status )
          IF(status>0) THEN
            temp = treal(i)
            GOTO 810
          ENDIF
        ENDDO
        !Radius of sphere must be in Angstroms
        READ(treal(4),*,END=800,ERR=800) region_2(1)
        region_2(2) = 0.d0
        region_2(3) = 0.d0
      ELSEIF(region_geom=='cylinder') THEN
        READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
            & region_side, region_geom, region_dir, treal(1), treal(2), treal(4)
        !Check if numbers contain a keyword like "BOX" or "INF"
        SELECT CASE(region_dir)
        CASE('x','X')
          CALL BOX2DBLE( H(:,2) , treal(1) , region_1(1) , status )
          IF(status>0) THEN
            temp = treal(1)
            GOTO 810
          ENDIF
          CALL BOX2DBLE( H(:,3) , treal(2) , region_1(2) , status )
          IF(status>0) THEN
            temp = treal(2)
            GOTO 810
          ENDIF
        CASE('y','Y')
          CALL BOX2DBLE( H(:,3) , treal(1) , region_1(1) , status )
          IF(status>0) THEN
            temp = treal(1)
            GOTO 810
          ENDIF
          CALL BOX2DBLE( H(:,1) , treal(2) , region_1(2) , status )
          IF(status>0) THEN
            temp = treal(2)
            GOTO 810
          ENDIF
        CASE('z','Z')
          CALL BOX2DBLE( H(:,1) , treal(1) , region_1(1) , status )
          IF(status>0) THEN
            temp = treal(1)
            GOTO 810
          ENDIF
          CALL BOX2DBLE( H(:,2) , treal(2) , region_1(2) , status )
          IF(status>0) THEN
            temp = treal(2)
            GOTO 810
          ENDIF
        CASE DEFAULT
          nerr=nerr+1
          GOTO 1000
        END SELECT
        !Radius of cylinder must be in Angstroms
        READ(treal(4),*,END=800,ERR=800) region_2(1)
        region_1(3) = 0.d0
        region_2(2) = 0.d0
        region_2(3) = 0.d0
      ENDIF
    ELSEIF( region_side=="prop" ) THEN
      !store property name in "region_geom"
      READ(options_array(ioptions),*,END=800,ERR=800) optionname, region_side, region_geom, temp
      !store value(s) of the property into 
      temp = ADJUSTL(temp)
      strlength = SCAN(temp,":")
      IF( strlength>0 ) THEN
        !user gives a range => read min and max values (can be -INF or +INF)
        treal(1) = temp(1:strlength-1)
        CALL BOX2DBLE( H(:,1) , treal(1) , region_1(1) , status )
        treal(2) = temp(strlength+1:)
        CALL BOX2DBLE( H(:,1) , treal(2) , region_1(2) , status )
        region_2(1) = 10.d0
        !make sure region_1(1) is smaller than region_1(2)
        IF( region_1(1) > region_1(2) ) THEN
          tempreal = region_1(2)
          region_1(2) = region_1(1)
          region_1(1) = tempreal
        ENDIF
      ELSE
        !user gives only one value => save it to region_1(1)
        READ(temp,*,END=800,ERR=800) region_1(1)
        region_2(1) = 0.d0
      ENDIF
    ELSEIF( region_side=="random" ) THEN
      !store number of atoms N to region_1(1) and species to region_geom
      READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
          & region_side, treal(1), region_geom
      j = INDEX(treal(1),"%")
      IF( j>0 ) THEN
        !treal(1) contains a percentage
        !=> append the symbol "%" to region_side
        region_side = "random%"
        !Read the number from treal(1)
        treal(1)(j:j) = " "
        READ(treal(1),*,END=800,ERR=800) region_1(1)
        region_1(1) = region_1(1) / 100.d0
      ELSE
        READ(treal(1),*,END=800,ERR=800) region_1(1)
      ENDIF
      !Read the number in treal(1) and save it into region_1(1)
      !CALL BOX2DBLE( (/1.d0,0.d0,0.d0/) , treal(1) , region_1(1) , status )
      !IF(status>0) THEN
      !  temp = treal(1)
      !  GOTO 810
      !ENDIF
    ELSEIF( region_side=="neigh" ) THEN
      !Store number of neighbors, or cutoff radius for neighbor search, into region_1(1)
      !Store species of neighbors (can be "all" or "any") in region_geom
      !Store index of atoms whose neighbors must be searched in region_1(2)
      READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
          & region_side, region_1(1), region_geom, region_1(2)
    ENDIF
    IF(nerr.NE.0) GOTO 1000
    CALL SELECT_XYZ(H,P,AUXNAMES,AUX,region_side,region_geom,region_dir,region_1,region_2,ORIENT,SELECT)
  !
  CASE('-shear')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
        & shear_surf, treal(1), shear_dir
    !treal may contain the shear strain expressed in percent, e.g. "3%"
    SELECT CASE(shear_surf)
    CASE('x','X')
      i=1
    CASE('y','Y')
      i=2
    CASE('z','Z')
      i=3
    END SELECT
    CALL BOX2DBLE( H(:,i) , treal(1) , shear_strain , status )
    IF(status>0) then
      temp = treal(1)
      goto 810
    endif
    CALL SHEAR_XYZ(H,P,S,shear_surf,shear_strain,shear_dir)
  !
  CASE('-shift')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, &
        & shift_dir, treal(4), shift_axis, treal(1), treal(2), treal(3)
    !Check if numbers contain a keyword like "BOX" or "INF"
    i=1
    SELECT CASE(shift_axis)
    CASE('x','X')
      i=1
    CASE('y','Y')
      i=2
    CASE('z','Z')
      i=3
    END SELECT
    CALL BOX2DBLE( H(:,i) , treal(4) , shift_dist , status )
    IF(status>0) THEN
      temp = treal(4)
      GOTO 810
    ENDIF
    CALL BOX2DBLE( H(:,1) , treal(1) , shift_tau1 , status )
    IF(status>0) THEN
      temp = treal(1)
      GOTO 810
    ENDIF
    CALL BOX2DBLE( H(:,2) , treal(2) , shift_tau2 , status )
    IF(status>0) THEN
      temp = treal(2)
      GOTO 810
    ENDIF
    CALL BOX2DBLE( H(:,3) , treal(3) , shift_tau3 , status )
    IF(status>0) THEN
      temp = treal(3)
      GOTO 810
    ENDIF
    CALL SHIFT_XYZ(P,S,shift_dir,shift_dist,shift_axis,shift_tau1,shift_tau2,shift_tau3,ORIENT,SELECT)
  !
  CASE('-sort')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, sortcol, sortorder
    CALL SORT_XYZ(P,S,AUXNAMES,AUX,SELECT,sortcol,sortorder)
  !
  CASE('-substitute', '-sub')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, sp1, sp2
    CALL SUBSTITUTE_XYZ(P,S,sp1,sp2,SELECT)
  !
  CASE('-unit', '-u')
    READ(options_array(ioptions),'(a128)',END=800,ERR=800) temp  !optionname, unit1, unit2
    temp = ADJUSTL(temp)
    strlength = SCAN(temp,' ')
    temp = temp(strlength:)
    temp = ADJUSTL(temp)
    strlength = SCAN(temp,' ')
    unit1 = temp(1:strlength)
    temp = temp(strlength:)
    unit2 = ADJUSTL(temp)
    CALL UNIT_XYZ(H,P,S,AUXNAMES,AUX,unit1,unit2)
  !
  CASE('-unskew')
    CALL UNSKEW_XYZ(H)
  !
  CASE('-velocity')
    READ(options_array(ioptions),*,END=800,ERR=800) optionname, vel_T
    CALL VELOCITY_XYZ(P,AUXNAMES,AUX,vel_T,SELECT)
  !
  CASE('-wrap')
    CALL WRAP_XYZ(H,P,S,SELECT)
  !
  ! -- please add other options in alphabetical order --
  !
  CASE DEFAULT
    !Could not understand option => display a warning
    nwarn = nwarn+1
    CALL ATOMSK_MSG(2700,(/options_array(ioptions)/),(/0.d0/))
  !
  END SELECT
  !
  !
  IF(verbosity==4) THEN
    !DEBUGGING
    !Write atoms positions to "atomsk.xyz"
    temp = 'atomsk.xyz'
    OPEN(UNIT=36,FILE=temp,STATUS="UNKNOWN",FORM="FORMATTED")
    WRITE(36,*) SIZE(P(:,1))
    msg = '#Debug file for atomsk, after option: '//TRIM(optionname)
    WRITE(36,*) TRIM(msg)
    DO i=1, SIZE(P(:,1))
      CALL ATOMSPECIES(P(i,4),species)
      WRITE(36,'(a2,2X,3(f16.8,1X))') species, P(i,1:3)
    ENDDO
    WRITE(36,'(a4)') 'alat'
    WRITE(36,'(a3)') '1.0'
    WRITE(36,'(a9)') 'supercell'
    WRITE(36,'(3f16.6)') H(1,:)
    WRITE(36,'(3f16.6)') H(2,:)
    WRITE(36,'(3f16.6)') H(3,:)
    CLOSE(36)
    !
    !Write elastic tensor to logfile
    msg = 'Elastic tensor (GPa):'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO i=1,9
      WRITE(msg,'(9(e10.3,2X))') (C_tensor(i,j), j=1,9)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
    !
    IF( ALLOCATED(AUXNAMES) ) THEN
      WRITE(msg,*) 'Size AUXNAMES, AUX:', SIZE(AUXNAMES), ',', &
                 & SIZE(AUX(:,1)), SIZE(AUX(1,:))
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      !Write auxiliary properties to logfile
      msg = 'Auxiliary properties (1-4):'
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      WRITE(msg,*) (TRIM(ADJUSTL(AUXNAMES(i)))//'  ', i=1,MIN(SIZE(AUXNAMES),6) )
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      DO i=1,MIN( 20,SIZE(AUX(:,1)) )
        WRITE(msg,'(6(f10.3,2X))') (AUX(i,j), j=1,MIN(SIZE(AUX(1,:)),6))
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDDO
      IF(i>=20) THEN
        WRITE(msg,*) '  (...discontinued...)'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDIF
    ENDIF
    !
  ENDIF !end v==4
  !
  !
  !In case of error break loop and exit
  IF(nerr.NE.0) RETURN
  !
  !
  !If shells exist, check that array S has a smaller size than P
  IF( ALLOCATED(S) ) THEN
    IF( SIZE(S,1) .NE. SIZE(P,1) ) THEN
      nerr = nerr+1
      CALL ATOMSK_MSG(2810,(/""/),(/0.d0/))
    ENDIF
  ENDIF
  !
  !If a selection exists, check that at least one atom is selected, otherwise clear selection
  IF( ALLOCATED(SELECT) ) THEN
    IF( SIZE(SELECT)<=0 .OR. .NOT.ANY(SELECT(:)) ) THEN
      nwarn=nwarn+1
      DEALLOCATE(SELECT)
      CALL ATOMSK_MSG(2750,(/""/),(/0.d0/))
    ENDIF
  ENDIF
  !
  !If an elastic tensor was defined, check it
  IF( ANY( C_tensor(1:3,1:3).NE.0.d0 ) ) THEN
    !Check its shape
    CALL CHECK_CTENSOR(C_tensor,i)
    IF(i.NE.0) THEN
      nwarn=nwarn+1
      CALL ATOMSK_MSG(2740,(/""/),(/0.d0/))
    ENDIF
    !Also check that there are no "NaN" values
    CALL CHECKNAN(C_tensor,i)
    IF(i.NE.0) THEN
      nerr=nerr+1
      CALL ATOMSK_MSG(2809,(/""/),(/0.d0/))
      GOTO 1000
    ENDIF
  ENDIF
  !
  !
!
ENDDO
GOTO 1000
!
!
!
800 CONTINUE
nerr=nerr+1
temp = options_array(ioptions)
CALL ATOMSK_MSG(2806,(/TRIM(temp)/),(/0.d0/))
CALL DISPLAY_HELP(temp)
GOTO 1000
!
810 CONTINUE
CALL ATOMSK_MSG(2813,(/TRIM(temp)/),(/0.d0/))
GOTO 1000
!
!
!
1000 CONTINUE
msg = 'EXITING OPTIONS'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
END SUBROUTINE OPTIONS_AFF
!
!
END MODULE
