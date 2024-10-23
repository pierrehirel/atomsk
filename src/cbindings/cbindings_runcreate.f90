MODULE cbindings_runcreate
!
!**********************************************************************************
!*  RUN CREATE                                                                    *
!**********************************************************************************
!* This module defines the cbindings for the CREATE_CELL subroutine.              *
!* Additionally, the DEALLOCATE_CREATE routine deallocates memory used by c       *
!* pointers.                                                                      *
!*                                                                                *
!* WHAT MUST BE PROVIDED:                                                         *
!* c_create_a0                                                                    *
!*           a vector of C_DOUBLEs that defines the lattice constants.            *
!* c_create_struc                                                                 *
!*           a vector of C_CHARs of length 10 that defines the lattice type (e.g. *
!*           'fcc', 'bcc', etc.) and any unneeded space should be set to ' '.     *
!* c_create species                                                               *
!*           a vector of C_CHARs of length 40 that defines the atomic species in  *
!*           the unit cell. Each species should occupy two C_CHARs and any        *
!*           unneeded space should be set to ' '.                                 *
!* c_NT_mn   a vector of C_INTs of length 2 that define the nanotube.             *
!* c_create_Miller                                                                *
!*           a vector of C_CHARs of length 96 that define the Miller indices of   *
!*           the system. Each Miller index should occupy 32 C_CHARs and any       *
!*           unneeded space should be set to ' '.                                 *
!* NUMOPTIONS                                                                     *
!*           a C_INT that defines the number of options to be applied to the      *
!*           system.                                                              *
!* c_options_array                                                                *
!*           a vector of C_CHARs of length 128*NUMOPTIONS that contains the       *
!*           options to be applied. Each set of 128 characters defines one option *
!*           and any unneeded space should be set to ' '.                         *
!* NUMATOMS  a C_INT that defines the number of atoms in the system.              *
!* c_POUT    a C_PTR that will point to the updated coordinates and atomic        *
!*           numbers of the system after the function is run.                     *
!* c_H       an array of contiguous C_DOUBLEs of size (3,3) that define the base  *
!*           vectors of the supercell (often this is not necessary and can be     *
!*           zeros).                                                              *
!*                                                                                *
!* WHAT IS RETURNED BY THIS ROUTINE:                                              *
!* modified C_PTR c_POUT                                                          *
!* modified array c_H                                                             *
!*                                                                                *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: J. Meziere - 25 Sep. 2024                                   *
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
SUBROUTINE RUN_CREATE(c_create_a0,c_create_struc,c_create_species,c_NT_mn,c_create_Miller,&
                     &NUMOPTIONS,c_options_array,&
                     &NUMATOMS,c_POUT,&
                     &c_H) BIND(C)
!
IMPLICIT NONE
!Declare C variables
REAL(C_DOUBLE),DIMENSION(3):: c_create_a0
CHARACTER(KIND=C_CHAR),DIMENSION(10):: c_create_struc
CHARACTER(KIND=C_CHAR),DIMENSION(40):: c_create_species
INTEGER(C_INT),DIMENSION(2):: c_NT_mn
CHARACTER(KIND=C_CHAR),DIMENSION(96):: c_create_Miller
INTEGER(C_INT):: NUMOPTIONS
CHARACTER(KIND=C_CHAR),DIMENSION(*):: c_options_array
INTEGER(C_INT):: NUMATOMS
TYPE(C_PTR),TARGET:: c_POUT
REAL(C_DOUBLE),DIMENSION(3,3):: c_H
!Declare Fortran variables
REAL(dp),DIMENSION(3):: create_a0               !the lattice constants
CHARACTER(LEN=10):: create_struc  !structure to create (fcc, bcc...)
CHARACTER(LEN=2),DIMENSION(20):: create_species !chemical species of atoms
INTEGER,DIMENSION(2):: NT_mn
CHARACTER(LEN=32),DIMENSION(3):: create_Miller    !Miller vectors along X, Y, Z
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P  !atomic positions
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
j=0
DO i=1,3
  create_Miller(i)=TRANSFER(c_create_Miller(32*j+1:32*(j+1)), create_Miller(i))
  j=j+1
END DO
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
CALL CREATE_CELL(create_a0,create_struc,create_species,NT_mn,create_Miller,options_array,outputfile,outfileformats,wof,H,P)
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
