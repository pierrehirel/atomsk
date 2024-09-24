MODULE cbindings_runcreate
!
!**********************************************************************************
!*  RUN CREATE                                                                    *
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
USE mode_create
!
!Module for writing output files
USE writeout
!
REAL(dp),POINTER:: PPOINT(:,:)  !atomic positions
!
CONTAINS
!
SUBROUTINE RUN_CREATE(c_create_a0,c_create_struc,c_create_species,c_NT_mn,&
                     &NUMOPTIONS,c_options_array,&
                     &NUMATOMS,c_POUT,&
                     &c_H,c_ORIENT) BIND(C)
!
IMPLICIT NONE
!Declare C variables
REAL(C_DOUBLE),DIMENSION(3):: c_create_a0
CHARACTER(KIND=C_CHAR),DIMENSION(10):: c_create_struc
CHARACTER(KIND=C_CHAR),DIMENSION(40):: c_create_species
INTEGER(C_INT),DIMENSION(2):: c_NT_mn
INTEGER(C_INT):: NUMOPTIONS
CHARACTER(KIND=C_CHAR),DIMENSION(*):: c_options_array
INTEGER(C_INT):: NUMATOMS
TYPE(C_PTR),TARGET:: c_POUT
REAL(C_DOUBLE),DIMENSION(3,3):: c_H
REAL(C_DOUBLE),DIMENSION(3,3):: c_ORIENT
!Declare Fortran variables
REAL(dp),DIMENSION(3):: create_a0               !the lattice constants
CHARACTER(LEN=10):: create_struc  !structure to create (fcc, bcc...)
CHARACTER(LEN=2),DIMENSION(20):: create_species !chemical species of atoms
INTEGER,DIMENSION(2):: NT_mn
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P  !atomic positions
REAL(dp),DIMENSION(3,3):: ORIENT    !crystalographic orientation
CHARACTER(LEN=4096):: outputfile
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to output
LOGICAL:: wof !write output file?
INTEGER:: i
INTEGER:: j
!Transfer from c-style to fortran-style variables
create_a0=TRANSFER(c_create_a0, create_a0)
create_struc=TRANSFER(c_create_struc, create_struc)
j=0
DO i=1,20
  create_species(i)=TRANSFER(c_create_species(2*j+1:2*(j+1)), create_species(i))
  j=j+1
END DO
NT_mn=TRANSFER(c_NT_mn, NT_mn)
ALLOCATE(options_array(NUMOPTIONS))
DO i=1,NUMOPTIONS
  !Set up options array from passed in options
  options_array(i)=TRANSFER(c_options_array(128*j+1:128*(j+1)),options_array(i))
  j=j+1
END DO
DO i=1,3
  H(i,:)=TRANSFER(c_H(i,:),H(i,:))
  ORIENT(i,:)=TRANSFER(c_ORIENT(i,:),ORIENT(i,:))
END DO
wof=.FALSE.
!Run creation
CALL CREATE_CELL(create_a0,create_struc,create_species,NT_mn,ORIENT,options_array,outputfile,outfileformats,wof,H,P)
!Save to c-vars
NUMATOMS=SIZE(P,1)
ALLOCATE(PPOINT(NUMATOMS,4))
DO i=1,NUMATOMS
  PPOINT(i,:)=P(i,:)
END DO
c_POUT=C_LOC(PPOINT)
DO i=1,3
  c_H(i,:)=TRANSFER(H(i,:),c_H(i,:))
  c_ORIENT(i,:)=TRANSFER(ORIENT(i,:),c_ORIENT(i,:))
END DO
!
END SUBROUTINE RUN_CREATE
!
SUBROUTINE DEALLOCATE_CREATE(POUT) BIND(C)
IMPLICIT NONE
TYPE(C_PTR)::POUT
!
DEALLOCATE(PPOINT)
POUT=C_NULL_PTR
END SUBROUTINE DEALLOCATE_CREATE
!
END MODULE cbindings_runcreate
