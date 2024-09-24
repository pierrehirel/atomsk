MODULE cbindings_runoptions
!
!**********************************************************************************
!*  RUN OPTIONS                                                                   *
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
USE options
!
!Module for writing output files
USE writeout
!
REAL(dp),POINTER:: PPOINT(:,:)  !atomic positions
REAL(dp),POINTER:: SPOINT(:,:)  !shell positions (if any)
REAL(dp),POINTER:: AUXPOINT(:,:) !auxiliary properties of atoms/shells
LOGICAL,POINTER:: SELECTPOINT(:)  !mask for atom list
!
CONTAINS
!
SUBROUTINE RUN_OPTIONS(NUMOPTIONS,c_options_array,&
                      &NUMATOMS,c_PIN,c_POUT,&
                      &NUMSHELLS,c_SIN,c_SOUT,&
                      &NUMAUXNAMES,c_AUXNAMES,c_AUXIN,c_AUXOUT,&
                      &c_Huc,c_H,c_ORIENT,c_SELECTIN,c_SELECTOUT,c_C_tensor) BIND(C)
!
IMPLICIT NONE
!Declare C variables
INTEGER(C_INT):: NUMOPTIONS
CHARACTER(KIND=C_CHAR),DIMENSION(*):: c_options_array
INTEGER(C_INT):: NUMATOMS
REAL(C_DOUBLE),DIMENSION(NUMATOMS,4):: c_PIN
TYPE(C_PTR),TARGET:: c_POUT
INTEGER(C_INT):: NUMSHELLS
REAL(C_DOUBLE),DIMENSION(NUMSHELLS,4):: c_SIN
TYPE(C_PTR):: c_SOUT
INTEGER(C_INT):: NUMAUXNAMES
CHARACTER(KIND=C_CHAR),DIMENSION(*):: c_AUXNAMES
REAL(C_DOUBLE),DIMENSION(NUMATOMS,NUMAUXNAMES):: c_AUXIN
TYPE(C_PTR):: c_AUXOUT
REAL(C_DOUBLE),DIMENSION(3,3):: c_Huc
REAL(C_DOUBLE),DIMENSION(3,3):: c_H
REAL(C_DOUBLE),DIMENSION(3,3):: c_ORIENT
INTEGER(C_INT),DIMENSION(NUMATOMS):: c_SELECTIN
TYPE(C_PTR):: c_SELECTOUT
REAL(C_DOUBLE),DIMENSION(9,9):: c_C_tensor
!Declare Fortran variables
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
REAL(dp),DIMENSION(3,3):: Huc !Base vectors of the unit cell
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P  !atomic positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S  !shell positions (if any)
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties of atoms/shells
REAL(dp),DIMENSION(3,3):: ORIENT    !crystalographic orientation
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
REAL(dp),DIMENSION(9,9):: C_tensor  !stiffness tensor
INTEGER:: i
INTEGER:: j
!
!Transfer from c-style to fortran-style variables
ALLOCATE(P(NUMATOMS,4))
ALLOCATE(AUX(NUMATOMS,NUMAUXNAMES))
ALLOCATE(SELECT(NUMATOMS))
DO i=1,NUMATOMS
  P(i,:)=TRANSFER(c_PIN(i,:),P(i,:))
  AUX(i,:)=TRANSFER(c_AUXIN(i,:),AUX(i,:))
  SELECT(i)=TRANSFER(c_SELECTIN(i),SELECT(i))
END DO
ALLOCATE(S(NUMSHELLS,4))
DO i=1,NUMSHELLS
  S(i,:)=TRANSFER(c_SIN(i,:),S(i,:))
END DO
ALLOCATE(AUXNAMES(NUMAUXNAMES))
j=0
DO i=1,NUMAUXNAMES
  AUXNAMES(i)=TRANSFER(c_AUXNAMES(128*j+1:128*(j+1)),AUXNAMES(i))
  j=j+1
END DO
DO i=1,3
  Huc(i,:)=TRANSFER(c_Huc(i,:),Huc(i,:))
  H(i,:)=TRANSFER(c_H(i,:),H(i,:))
  ORIENT(i,:)=TRANSFER(c_ORIENT(i,:),ORIENT(i,:))
END DO
DO i=1,9
  C_tensor(i,:)=TRANSFER(c_C_tensor(i,:),C_tensor(i,:))
END DO
j=0
ALLOCATE(options_array(NUMOPTIONS))
DO i=1,NUMOPTIONS
  !Set up options array from passed in options
  options_array(i)=TRANSFER(c_options_array(128*j+1:128*(j+1)),options_array(i))
  j=j+1
END DO
!Run all options
CALL OPTIONS_AFF(options_array,Huc,H,P,S,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
!Save to c-style pointers
NUMATOMS=SIZE(P,1)
ALLOCATE(PPOINT(NUMATOMS,4))
ALLOCATE(AUXPOINT(NUMATOMS,NUMAUXNAMES))
ALLOCATE(SELECTPOINT(NUMATOMS))
DO i=1,NUMATOMS
  PPOINT(i,:)=P(i,:)
  AUXPOINT(i,:)=AUX(i,:)
  SELECTPOINT(i)=SELECT(i)
END DO
c_POUT=C_LOC(PPOINT)
c_AUXOUT=C_LOC(AUXPOINT)
c_SELECTOUT=C_LOC(SELECTPOINT)
NUMSHELLS=SIZE(S,1)
ALLOCATE(SPOINT(NUMSHELLS,4))
DO i=1,NUMSHELLS
  SPOINT(i,:)=S(i,:)
END DO
c_SOUT=C_LOC(SPOINT)
DO i=1,3
  c_Huc(i,:)=TRANSFER(Huc(i,:),c_Huc(i,:))
  c_H(i,:)=TRANSFER(H(i,:),c_H(i,:))
  c_ORIENT(i,:)=TRANSFER(ORIENT(i,:),c_ORIENT(i,:))
END DO
DO i=1,9
  c_C_tensor(i,:)=TRANSFER(C_tensor(i,:),c_C_tensor(i,:))
END DO
!
END SUBROUTINE RUN_OPTIONS
!
SUBROUTINE DEALLOCATE_OPTIONS(POUT, SOUT, AUXOUT, SELECTOUT) BIND(C)
IMPLICIT NONE
TYPE(C_PTR)::POUT
TYPE(C_PTR)::SOUT
TYPE(C_PTR)::AUXOUT
TYPE(C_PTR)::SELECTOUT
!
DEALLOCATE(PPOINT)
DEALLOCATE(SPOINT)
DEALLOCATE(AUXPOINT)
DEALLOCATE(SELECTPOINT)
POUT=C_NULL_PTR
SOUT=C_NULL_PTR
AUXOUT=C_NULL_PTR
SELECTOUT=C_NULL_PTR
END SUBROUTINE DEALLOCATE_OPTIONS
!
END MODULE cbindings_runoptions
