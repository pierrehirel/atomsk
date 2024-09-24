MODULE cbindings_runwrite
!
!**********************************************************************************
!*  RUN WRITE                                                                     *
!**********************************************************************************
!* This module defines the cbindings for the RUN_WRITE_AFF subroutine             *
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
USE writeout
!
CONTAINS
!
SUBROUTINE RUN_WRITE(NUMATOMS,c_P,NUMSHELLS,c_S,&
                    &NUMAUXNAMES,c_AUXNAMES,c_AUX,&
                    &c_H,c_prefix,PREFIXSIZE,c_fileformats,NUMFORMATS) BIND(C)
!
IMPLICIT NONE
!Declare C variables
INTEGER(C_INT):: NUMATOMS
REAL(C_DOUBLE),DIMENSION(NUMATOMS,4):: c_P
INTEGER(C_INT):: NUMSHELLS
REAL(C_DOUBLE),DIMENSION(NUMSHELLS,4):: c_S
INTEGER(C_INT):: NUMAUXNAMES
CHARACTER(KIND=C_CHAR),DIMENSION(*):: c_AUXNAMES
REAL(C_DOUBLE),DIMENSION(NUMATOMS,NUMAUXNAMES):: c_AUX
REAL(C_DOUBLE),DIMENSION(3,3):: c_H
INTEGER(C_INT)::PREFIXSIZE
CHARACTER(KIND=C_CHAR),DIMENSION(PREFIXSIZE):: c_prefix
INTEGER(C_INT):: NUMFORMATS
CHARACTER(KIND=C_CHAR),DIMENSION(*):: c_fileformats
!Declare Fortran variables
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P  !atomic positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S  !shell positions (if any)
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties of atoms/shells
CHARACTER(LEN=PREFIXSIZE):: prefix
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: fileformats
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
INTEGER:: i
INTEGER:: j
!
!Transfer from c-style to fortran-style variables
ALLOCATE(P(NUMATOMS,4))
ALLOCATE(AUX(NUMATOMS,NUMAUXNAMES))
DO i=1,NUMATOMS
  P(i,:)=TRANSFER(c_P(i,:),P(i,:))
  AUX(i,:)=TRANSFER(c_AUX(i,:),AUX(i,:))
END DO
ALLOCATE(S(NUMSHELLS,4))
DO i=1,NUMSHELLS
  S(i,:)=TRANSFER(c_S(i,:),S(i,:))
END DO
ALLOCATE(AUXNAMES(NUMAUXNAMES))
j=0
DO i=1,NUMAUXNAMES
  AUXNAMES(i)=TRANSFER(c_AUXNAMES(128*j+1:128*(j+1)),AUXNAMES(i))
  j=j+1
END DO
DO i=1,3
  H(i,:)=TRANSFER(c_H(i,:),H(i,:))
END DO
prefix = TRANSFER(c_prefix(:),prefix)
j=0
ALLOCATE(fileformats(NUMFORMATS))
DO i=1,NUMFORMATS
  fileformats(i)=TRANSFER(c_fileformats(5*j+1:5*(j+1)),fileformats(i))
  j=j+1
END DO
!Run write
CALL WRITE_AFF(prefix,fileformats,H,P,S,comment,AUXNAMES,AUX)
!
END SUBROUTINE RUN_WRITE
!
END MODULE cbindings_runwrite
