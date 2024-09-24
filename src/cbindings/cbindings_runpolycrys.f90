MODULE cbindings_runpolycrys
!
!**********************************************************************************
!*  RUN POLYCRYS                                                                  *
!**********************************************************************************
!* This module defines the cbindings for the RUN_OPTIONS subroutine  .            *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 16 Dec. 2015                                     *
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
!Module for options
USE ISO_C_BINDING
USE mode_polycrystal
!
!Module for writing output files
USE writeout
!
REAL(dp),POINTER:: PPOINT(:,:)  !atomic positions
!
CONTAINS
!
SUBROUTINE RUN_POLYCRYS(c_ucfile,UCFILESIZE,c_vfile,VFILESIZE,&
                       &NUMOPTIONS,c_options_array,&
                       &NUMATOMS,c_POUT,c_H) BIND(C)
!
IMPLICIT NONE
!Declare C variables
INTEGER(C_INT)::UCFILESIZE
CHARACTER(KIND=C_CHAR),DIMENSION(UCFILESIZE):: c_ucfile
INTEGER(C_INT)::VFILESIZE
CHARACTER(KIND=C_CHAR),DIMENSION(VFILESIZE):: c_vfile
INTEGER(C_INT):: NUMOPTIONS
CHARACTER(KIND=C_CHAR),DIMENSION(*):: c_options_array
INTEGER(C_INT):: NUMATOMS
TYPE(C_PTR),TARGET:: c_POUT
REAL(C_DOUBLE),DIMENSION(3,3):: c_H
!Declare Fortran variables
CHARACTER(LEN=UCFILESIZE):: ucfile  !name of file containing seed (usually a unit cell, but can be anything)
CHARACTER(LEN=VFILESIZE):: vfile   !name of file containing parameters for Voronoi construction
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P  !atomic positions
CHARACTER(LEN=4096):: prefix
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to output
LOGICAL:: wof !write output file?
INTEGER:: i
INTEGER:: j
!Transfer from c-style to fortran-style variables
ucfile=TRANSFER(c_ucfile, ucfile)
vfile=TRANSFER(c_vfile, vfile)
ALLOCATE(options_array(NUMOPTIONS))
DO i=1,NUMOPTIONS
  !Set up options array from passed in options
  options_array(i)=TRANSFER(c_options_array(128*j+1:128*(j+1)),options_array(i))
  j=j+1
END DO
DO i=1,3
  H(i,:)=TRANSFER(c_H(i,:),H(i,:))
END DO
wof=.FALSE.
!Run creation
CALL POLYCRYS(ucfile,vfile,options_array,prefix,outfileformats,wof,H,P)
!Save to c-vars
NUMATOMS=SIZE(P,1)
ALLOCATE(PPOINT(NUMATOMS,4))
DO i=1,NUMATOMS
  PPOINT(i,:)=P(i,:)
END DO
c_POUT=C_LOC(PPOINT)
DO i=1,3
  c_H(i,:)=TRANSFER(H(i,:),c_H(i,:))
END DO
!
END SUBROUTINE RUN_POLYCRYS
!
SUBROUTINE DEALLOCATE_POLYCRYS(POUT) BIND(C)
IMPLICIT NONE
TYPE(C_PTR)::POUT
!
DEALLOCATE(PPOINT)
POUT=C_NULL_PTR
END SUBROUTINE DEALLOCATE_POLYCRYS
!
END MODULE cbindings_runpolycrys
