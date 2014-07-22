MODULE readin
!
!**********************************************************************************
!*  READIN                                                                        *
!**********************************************************************************
!* This module reads a single file and loads atomic data to memory.               *
!*                                                                                *
!* First the input file format is detected by the "guess_format" subroutine,      *
!* and then the appropriate module (in_*) is called for reading the file.         *
!*                                                                                *
!* WHAT MUST BE PROVIDED:                                                         *
!* inputfile                                                                      *
!*           string containing the name of the input file to read.                *
!*                                                                                *
!*                                                                                *
!* WHAT IS RETURNED BY THIS ROUTINE:                                              *
!* H         3 x 3 real array containing vectors of the supercell. H(1,:) is      *
!*           the first vector, H(2,:) the second, H(3,:) the third.               *
!* P         N x 4 real array containing atom positions. N=number of atoms,       *
!*           each line contains the coordinates (X,Y,Z) and atomic number.        *
!* S         if allocated, N x 4 real array containing positions of shells        *
!*           (for ionic core/shell model), and number of shells N must be         *
!*           identical to N in P. If S is not allocated, then no shell exist.     *
!* comment   an array of strings that can contain anything.                       *
!* AUXNAMES  if allocated, string array of dimension M, M being the number of     *
!*           existing auxiliary properties (must be equal to M in AUX). If        *
!*           AUXNAMES is not allocated, then AUX must also be non-allocated.      *
!* AUX       if allocated, N x M real array containing values of auxiliary        *
!*           properties of atoms, e.g. their velocity, forces, etc. The first     *
!*           dimension N must be equal to the number of atoms N in P, the second  *
!*           dimension M must be equal to the number of auxiliary properties.     *
!*           If AUX is not allocated then no auxiliary property exists, and       *
!*           AUXNAMES must also be non-allocated.                                 *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 22 July 2014                                     *
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
USE messages
USE files
USE subroutines
USE guess_form
!
!Modules managing input files
USE in_atsk
USE in_bop
USE in_cfg
USE in_cif
USE in_dlp_cfg
USE in_gulp_gin
USE in_imd
USE in_lmp_c
USE in_lmp_data
USE in_mbpp_coorat
USE in_moldy
USE in_pdb
USE in_qe_pw
USE in_siesta_xv
USE in_vasp_poscar
USE in_xmd
USE in_xyz
USE in_xsf
! -- add other formats here --
!
!
!
CONTAINS
!
SUBROUTINE READ_AFF(inputfile,H,P,S,comment,AUXNAMES,AUX)
!
!
!
!Declare variables
IMPLICIT NONE
CHARACTER(LEN=*):: inputfile
CHARACTER(LEN=5):: infileformat
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(OUT),OPTIONAL:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(OUT),OPTIONAL:: comment
CHARACTER(LEN=1024):: msg
INTEGER:: i, strlength
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: P !positions of atoms (or ionic cores)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT),OPTIONAL:: S !positions of ionic shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT),OPTIONAL:: AUX !auxiliary properties
!
!
!Initialize variables
H(:,:) = 0.d0
IF(ALLOCATED(comment)) DEALLOCATE(comment)
IF(ALLOCATED(P)) DEALLOCATE(P)
IF(ALLOCATED(S)) DEALLOCATE(S)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
infileformat=''
strlength=0
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  CHECK THE FILES, IDENTIFY FORMATS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
100 CONTINUE
msg = 'READ_AFF'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
msg = 'Inputfile: '//TRIM(inputfile)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!If no input file was specified, or if it does not exist, then prompt the user
msg = 'Check if file exists: '//TRIM(inputfile)
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
CALL CHECKFILE(inputfile,'read')
!
120 CONTINUE
CALL ATOMSK_MSG(1000,(/TRIM(ADJUSTL(inputfile))/),(/0.d0/))
!
!!Open input file
!!and find out what is its format
!Note: the format of the input file is not guessed only
!     from its extension, but also from its contents; see
!     the routine GUESS_FORMAT for more details.
CALL GUESS_FORMAT(inputfile,infileformat,'read')
IF(infileformat=='xxx') THEN
  CALL ATOMSK_MSG(1800,(/''/),(/0.d0/))
  nerr=nerr+1
  GOTO 1000
ENDIF
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  READ INPUT FILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
200 CONTINUE
msg = 'Calling READ module...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!!Then, read input file
! by calling the module corresponding to the input format
SELECT CASE(infileformat)
CASE('atsk')
  CALL READ_ATSK(inputfile,H,P,S,comment,AUXNAMES,AUX)
CASE('bop')
  CALL READ_BOP(inputfile,H,P,comment,AUXNAMES,AUX)
CASE('cfg')
  CALL READ_CFG(inputfile,H,P,comment,AUXNAMES,AUX)
CASE('cif')
  CALL READ_CIF(inputfile,H,P,comment,AUXNAMES,AUX)
CASE('coo')
  CALL READ_COORAT(inputfile,H,P,comment,AUXNAMES,AUX)
CASE('dlp')
  CALL READ_DLP_CFG(inputfile,H,P,S,comment,AUXNAMES,AUX)
CASE('gin')
  CALL READ_GIN(inputfile,H,P,S,comment,AUXNAMES,AUX)
CASE('imd')
  CALL READ_IMD(inputfile,H,P,comment,AUXNAMES,AUX)
CASE('lmc')
  CALL READ_LMP_CUSTOM(inputfile,H,P,comment,AUXNAMES,AUX)
CASE('lmp')
  CALL READ_LMP_DATA(inputfile,H,P,comment,AUXNAMES,AUX)
CASE('mol')
  CALL READ_MOLDY(inputfile,H,P,comment,AUXNAMES,AUX)
CASE('pdb')
  CALL READ_PDB(inputfile,H,P,comment,AUXNAMES,AUX)
CASE('pos')
  CALL READ_POSCAR(inputfile,H,P,comment,AUXNAMES,AUX)
CASE('pw')
  CALL READ_QEPW(inputfile,H,P,comment,AUXNAMES,AUX)
CASE('xmd')
  CALL READ_XMD(inputfile,H,P,comment,AUXNAMES,AUX)
CASE('xsf')
  CALL READ_XSF(inputfile,H,P,comment,AUXNAMES,AUX)
CASE('xv','XV')
  CALL READ_XV(inputfile,H,P,comment,AUXNAMES,AUX)
CASE('xyz')
  CALL READ_XYZ(inputfile,H,P,comment,AUXNAMES,AUX)
!
! -- please add other formats in alphabetical order --
!
CASE DEFAULT
  GOTO 800
END SELECT
IF(nerr>=1) GOTO 800
!
IF(.NOT.ALLOCATED(P)) nerr=nerr+1
IF(nerr>=1) GOTO 800
!
!Check sizes of the arrays
IF(ALLOCATED(S)) THEN
  IF( SIZE(S,1) > SIZE(P,1) ) THEN
    !There cannot be more shells than cores => exit
    nerr=nerr+1
    CALL ATOMSK_MSG(1802,(/'S'/),(/0.d0/))
  ELSEIF( SIZE(S,1)<=0 ) THEN
    !no shell => S should not be allocated
    DEALLOCATE(S)
  ENDIF
ENDIF
IF( ALLOCATED(AUXNAMES) .OR. ALLOCATED(AUX) ) THEN
  !first dimension of AUX must correspond to number of atoms
  IF( SIZE(AUX,1) .NE. SIZE(P,1) ) THEN
    nerr=nerr+1
    CALL ATOMSK_MSG(1803,(/'AUX'/),(/0.d0/))
  ENDIF
  !second dimension of AUX must correspond to number of auxiliary properties
  IF( SIZE(AUXNAMES) .NE. SIZE(AUX,2) ) THEN
    nerr=nerr+1
    CALL ATOMSK_MSG(1802,(/'AUXNAMES'/),(/0.d0/))
  ENDIF
  !avoid arrays allocated with zero size
  IF( SIZE(AUXNAMES)<=0 .OR. SIZE(AUX,1)<=0 .OR. SIZE(AUX,2)<=0 ) THEN
    IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
    IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
  ENDIF
ENDIF
!
IF(verbosity==4) THEN
  IF(ALLOCATED(AUX)) THEN
    WRITE(msg,'(a5,12(a9,1X))') "AUXNAMES: ", (TRIM(AUXNAMES(i))//' ', i=1,SIZE(AUXNAMES(:)) )
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO strlength=1,MIN(20,SIZE(AUX,1))
      WRITE(msg,'(a5,20e10.3)') "AUX: ", (AUX(strlength,i), i=1,SIZE(AUX(1,:)) )
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
    IF(SIZE(AUX,1)>20) THEN
      WRITE(msg,*) "AUX:    ... discontinued ..."
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDIF
  ENDIF
ENDIF
CALL ATOMSK_MSG(1001,(/''/),(/0.d0/))
GOTO 1000
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  ERROR MESSAGES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
800 CONTINUE
CALL ATOMSK_MSG(1801,(/TRIM(inputfile)/),(/0.d0/))
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  TERMINATE MODULE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
1000 CONTINUE
!
!
END SUBROUTINE READ_AFF
!
!
END MODULE readin
