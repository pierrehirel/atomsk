MODULE mode_create
!
!**********************************************************************************
!*  MODE_CREATE                                                                   *
!**********************************************************************************
!* This module creates an atomic structure based on some keywords.                *
!* Unit cell vectors are return in the array H,                                   *
!* atom positions and atomic numbers in P.                                        *
!**********************************************************************************
!* (C) March 2010 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 24 June 2020                                     *
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
USE resize
USE out_xyz
USE options
USE writeout
!
CONTAINS
!
!
SUBROUTINE CREATE_CELL(create_a0,create_struc,create_species,NT_mn,ORIENT,options_array,outputfile,outfileformats,wof,H,P)
!
!
IMPLICIT NONE
!Input parameters
CHARACTER(LEN=2),DIMENSION(20),INTENT(IN):: create_species !chemical species of atoms
REAL(dp),DIMENSION(3),INTENT(IN):: create_a0               !the lattice constants
REAL(dp),DIMENSION(3,3):: ORIENT  !crystallographic orientation of the cell
!
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to output
CHARACTER(LEN=10):: create_struc  !structure to create (fcc, bcc...)
CHARACTER(LEN=32):: NT_type       !type of nanotube (zig-zag or armchair or chiral)
CHARACTER(LEN=128):: msg, temp
CHARACTER(LEN=4096):: outputfile
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
LOGICAL:: cubic, hexagonal, oriented   !is the system cubic? is it hcp? should it be oriented?
LOGICAL:: new
LOGICAL,INTENT(IN):: wof !write output file?
LOGICAL,DIMENSION(3):: orthovec  !are vectors orthogonal?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: i, j, k, l, m
INTEGER:: lminmax  !min and max values for loops when orienting cubic systems
INTEGER:: NP
INTEGER:: d, NT_NP, nspecies, r, t
INTEGER,DIMENSION(2):: NT_mn
REAL(dp):: H1, H2, H3
REAL(dp):: l_perp, l_para, NT_radius, vol, vol_cell, x, y, z1, z2  !for nanotubes
REAL(dp):: u, v, w
REAL(dp),DIMENSION(3):: a_perp, a_para, b_perp, b_para, coord, Hrecip !for nanotubes
REAL(dp),DIMENSION(1,4):: tempP  !temporary position
REAL(dp),DIMENSION(3,3):: Huc !Base vectors of the unit cell
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ips, uv     !interplanar spacing, unit vectors corresponding to new orientation ORIENT(:,:)
REAL(dp),DIMENSION(3,3):: ORIENTN     !normalized ORIENT
REAL(dp),DIMENSION(3,3):: rot_matrix, rot_matrix2  !rotation matrices
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, S  !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Q     !positions of atoms (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!Initialize variables
 cubic = .FALSE.
 hexagonal = .FALSE.
 oriented = .FALSE.
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
nspecies = 0
Huc(:,:) = 0.d0
H(:,:) = 0.d0
 C_tensor(:,:) = 0.d0
ips(:,:) = 0.d0
uv(:,:) = 0.d0
IF(ALLOCATED(P)) DEALLOCATE(P)
ALLOCATE(comment(1))
 comment(1)=''
!
!
CALL ATOMSK_MSG(4027,(/''/),(/0.d0/))
!
WRITE(msg,*) "lattice type: "//create_struc
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,'(a19,3f12.3)') "lattice constants: ", create_a0(:)
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
!
100 CONTINUE
!Make the create_struc variable straight
!in case the user mistyped 'diamand' or 'pervoskite'
SELECT CASE(create_struc(1:2))
CASE('di','Di','DI')
  create_struc = 'diamond'
CASE('pe','Pe','PE')
  create_struc = 'perovskite'
CASE('zi','Zi','zb','ZB')
  create_struc = 'zincblende'
CASE('na','NA','nt','NT')
  create_struc = 'nanotube'
CASE('wu','Wu','WU','wz','Wz','WZ')
  create_struc = 'wurtzite'
END SELECT
!
!Determine the number of species
nspecies=0
DO i=1,SIZE(create_species)
  IF(create_species(i).NE.'') nspecies=nspecies+1
ENDDO
!
!
!
200 CONTINUE
!Define base vectors H(:,:) and atom positions P(:,:) for the given lattice
SELECT CASE(create_struc)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  CUBIC LATTICES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('sc')
  cubic = .TRUE.
  IF(nspecies.NE.1) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0 /))
    GOTO 810
  ENDIF
  ALLOCATE(P(1,4))
  !Set up atom positions
  P(:,:) = 0.d0
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,2) = create_a0(2)
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up the messages
  WRITE(comment(1),*) 'Simple cubic '//TRIM(ADJUSTL((create_species(1))))
!
!
CASE('bcc','CsCl')
  cubic = .TRUE.
  IF(nspecies.NE.1 .AND. nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0,2.d0 /))
    GOTO 810
  ENDIF
  ALLOCATE(P(2,4))
  !Set up atom positions
  P(:,:) = 0.d0
  P(2,1) = 0.5d0
  P(2,2) = 0.5d0
  P(2,3) = 0.5d0
  P(:,1) = create_a0(1)*P(:,1)
  P(:,2) = create_a0(2)*P(:,2)
  P(:,3) = create_a0(3)*P(:,3)
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  IF(nspecies==2) THEN
    CALL ATOMNUMBER(create_species(2),P(2,4))
  ELSE
    P(2,4) = P(1,4)
  ENDIF
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,2) = create_a0(2)
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up the messages
  WRITE(comment(1),*) TRIM(create_species(1))
  IF(nspecies==1) THEN
    comment(1) = 'Bcc '//TRIM(ADJUSTL(comment(1)))
  ELSE
    comment(1) = 'Bcc '//TRIM(ADJUSTL(comment(1)))//TRIM(ADJUSTL(create_species(2)))//' alloy'
  ENDIF
!
!
CASE('fcc')
  cubic = .TRUE.
  IF(nspecies.NE.1 .AND. nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0,2.d0 /))
    GOTO 810
  ENDIF
  ALLOCATE(P(4,4))
  !Set up atom positions
  P(:,:) = 0.d0
  P(2,1) = 0.5d0
  P(2,2) = 0.5d0
   P(3,2) = 0.5d0
   P(3,3) = 0.5d0
  P(4,1) = 0.5d0
  P(4,3) = 0.5d0
  P(:,1) = create_a0(1)*P(:,1)
  P(:,2) = create_a0(2)*P(:,2)
  P(:,3) = create_a0(3)*P(:,3)
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  P(2,4) = P(1,4)
  IF(nspecies==2) THEN
    CALL ATOMNUMBER(create_species(2),P(3,4))
  ELSE
    P(3,4) = P(1,4)
  ENDIF
  P(4,4) = P(3,4)
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,2) = create_a0(2)
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up the messages
  WRITE(comment(1),*) TRIM(create_species(1))
  IF(nspecies==1) THEN
    comment(1) = 'Fcc '//TRIM(ADJUSTL(comment(1)))
  ELSE
    comment(1) = 'Fcc '//TRIM(ADJUSTL(comment(1)))//TRIM(ADJUSTL(create_species(2)))//' alloy'
  ENDIF
!
!
CASE('L12','L1_2')
  cubic = .TRUE.
  IF(nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0,2.d0 /))
    GOTO 810
  ENDIF
  ALLOCATE(P(4,4))
  !Set up atom positions
  P(:,:) = 0.d0
  P(1,1) = 0.5d0
  P(1,2) = 0.5d0
   P(2,2) = 0.5d0
   P(2,3) = 0.5d0
  P(3,1) = 0.5d0
  P(3,3) = 0.5d0
  P(:,1) = create_a0(1)*P(:,1)
  P(:,2) = create_a0(2)*P(:,2)
  P(:,3) = create_a0(3)*P(:,3)
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  P(2,4) = P(1,4)
  P(3,4) = P(1,4)
  CALL ATOMNUMBER(create_species(2),P(4,4))
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,2) = create_a0(2)
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up the messages
  WRITE(comment(1),*) TRIM(create_species(1))//"3"//TRIM(create_species(2))
  comment(1) = 'L12 '//TRIM(ADJUSTL(comment(1)))
!
!
CASE('dia','diamond','zincblende','zc')
  cubic = .TRUE.
  IF(nspecies.NE.1 .AND. nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0,2.d0 /))
    GOTO 810
  ENDIF
  ALLOCATE(P(8,4))
  !Set up atom positions
  P(:,:) = 0.d0
  P(2,1) = 0.5d0
  P(2,2) = 0.5d0
   P(3,2) = 0.5d0
   P(3,3) = 0.5d0
  P(4,1) = 0.5d0
  P(4,3) = 0.5d0
   P(5,1) = 0.25d0
   P(5,2) = 0.25d0
   P(5,3) = 0.25d0
  P(6,1) = 0.75d0
  P(6,2) = 0.75d0
  P(6,3) = 0.25d0
   P(7,1) = 0.75d0
   P(7,2) = 0.25d0
   P(7,3) = 0.75d0
  P(8,1) = 0.25d0
  P(8,2) = 0.75d0
  P(8,3) = 0.75d0
  P(:,1) = create_a0(1)*P(:,1)
  P(:,2) = create_a0(2)*P(:,2)
  P(:,3) = create_a0(3)*P(:,3)
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  IF(nspecies==2) THEN
    CALL ATOMNUMBER(create_species(2),P(5,4))
  ELSE
    P(5,4) = P(1,4)
  ENDIF
  DO i=2,4
    P(i,4) = P(1,4)
  ENDDO
  DO i=5,8
    P(i,4) = P(5,4)
  ENDDO
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,2) = create_a0(2)
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up the messages
  WRITE(comment(1),*) TRIM(create_species(1))
  IF(nspecies==1) THEN
    comment(1) = TRIM(ADJUSTL(comment(1)))//' with diamond structure'
  ELSE
    comment(1) = TRIM(ADJUSTL(comment(1)))//TRIM(ADJUSTL(create_species(2)))//' with zincblende structure'
  ENDIF
!
!
CASE('rocksalt','rs','B1')
  cubic = .TRUE.
  IF(nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 2.d0 /))
    GOTO 810
  ENDIF
  ALLOCATE(P(8,4))
  !Set up atom positions
  P(:,:) = 0.d0
  P(2,1) = 0.5d0
  P(2,2) = 0.5d0
   P(3,2) = 0.5d0
   P(3,3) = 0.5d0
  P(4,1) = 0.5d0
  P(4,3) = 0.5d0
   P(5,1) = 0.5d0
  P(6,2) = 0.5d0
   P(7,3) = 0.5d0
  P(8,1:3) = 0.5d0
  P(:,1) = create_a0(1)*P(:,1)
  P(:,2) = create_a0(2)*P(:,2)
  P(:,3) = create_a0(3)*P(:,3)
  !Set up atom species
  DO i=1,4
    CALL ATOMNUMBER(create_species(1),P(i,4))
    CALL ATOMNUMBER(create_species(2),P(i+4,4))
  ENDDO
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,2) = create_a0(2)
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up the messages
  WRITE(comment(1),*) TRIM(create_species(1))
  comment(1) = 'Rocksalt '//TRIM(ADJUSTL(comment(1)))//TRIM(create_species(2))
!
!
CASE('fluorite','fluorine')
  cubic = .TRUE.
  IF(nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0,2.d0 /))
    GOTO 810
  ENDIF
  ALLOCATE(P(12,4))
  !Set up atom positions
  ! 4 for the fcc classical lattice
  P(:,:) = 0.d0
  P(1,1) = 0.5d0
  P(1,2) = 0.5d0
  P(2,2) = 0.5d0
  P(2,3) = 0.5d0
  P(3,1) = 0.5d0
  P(3,3) = 0.5d0
  ! the 4 atoms at z=0.25
  P(5,1) = 0.25d0
  P(5,2) = 0.25d0
  P(5,3) = 0.25d0
  P(6,1) = 0.75d0
  P(6,2) = 0.25d0
  P(6,3) = 0.25d0
  P(7,1) = 0.25d0
  P(7,2) = 0.75d0
  P(7,3) = 0.25d0
  P(8,1) = 0.75d0
  P(8,2) = 0.75d0
  P(8,3) = 0.25d0
  ! the 4 atoms at z=0.75
  P(9,1) = 0.25d0
  P(9,2) = 0.25d0
  P(9,3) = 0.75d0
  P(10,1) = 0.75d0
  P(10,2) = 0.25d0
  P(10,3) = 0.75d0
  P(11,1) = 0.25d0
  P(11,2) = 0.75d0
  P(11,3) = 0.75d0
  P(12,1) = 0.75d0
  P(12,2) = 0.75d0
  P(12,3) = 0.75d0
  P(:,1) = create_a0(1)*P(:,1)
  P(:,2) = create_a0(2)*P(:,2)
  P(:,3) = create_a0(3)*P(:,3)
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  P(2,4) = P(1,4)
  P(3,4) = P(1,4)
  P(4,4) = P(1,4)
  CALL ATOMNUMBER(create_species(2),P(5,4))
  P(6,4) = P(5,4)
  P(7,4) = P(5,4)
  P(8,4) = P(5,4)
  P(9,4) = P(5,4)
  P(10,4) = P(5,4)
  P(11,4) = P(5,4)
  P(12,4) = P(5,4)
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,2) = create_a0(2)
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up the messages
  WRITE(comment(1),*) TRIM(create_species(1))//TRIM(create_species(2))//"2"
  comment(1) = "Fluorite "//TRIM(ADJUSTL(comment(1)))
!
!
CASE('c15','C15')
  cubic = .TRUE.
  IF(nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0,2.d0 /))
    GOTO 810
  ENDIF
  ALLOCATE(P(24,4))
  !Set up atom positions
  P(:,:) = 0.d0
  !Consider the prototype Cu2Mg
  !Cu atoms
  P(1,1:3) = 0.5d0
    P(2,1) = 0.5d0
    P(2,2) = 0.75d0
    P(2,3) = 0.75d0
  P(3,1) = 0.75d0
  P(3,2) = 0.5d0
  P(3,3) = 0.75d0
    P(4,1) = 0.75d0
    P(4,2) = 0.75d0
    P(4,3) = 0.5d0
  P(5,1) = 0.25d0
  P(5,2) = 0.5d0
  P(5,3) = 0.25d0
    P(6,1) = 0.5d0
    P(6,2) = 0.25d0
    P(6,3) = 0.25d0
  P(7,1) = 0.25d0
  P(7,2) = 0.25d0
  P(7,3) = 0.5d0
    P(8,1) = 0.25d0
    P(8,2) = 0.75d0
    P(8,3) = 0.d0
  P(9,1) = 0.5d0
  P(9,2) = 0.d0
  P(9,3) = 0.d0
    P(10,1) = 0.25d0
    P(10,2) = 0.d0
    P(10,3) = 0.75d0
  P(11,1) = 0.d0
  P(11,2) = 0.5d0
  P(11,3) = 0.d0
    P(12,1) = 0.75d0
    P(12,2) = 0.25d0
    P(12,3) = 0.d0
  P(13,1) = 0.d0
  P(13,2) = 0.25d0
  P(13,3) = 0.75d0
    P(14,1) = 0.d0
    P(14,2) = 0.75d0
    P(14,3) = 0.25d0
  P(15,1) = 0.75d0
  P(15,2) = 0.d0
  P(15,3) = 0.25d0
    P(16,1) = 0.d0
    P(16,2) = 0.d0
    P(16,3) = 0.5d0
  !Mg atoms are at positions of the diamond lattice
   P(17,1) = 0.125d0
   P(17,2) = 0.125d0
   P(17,3) = 0.125d0
  P(18,1) = 0.875d0
  P(18,2) = 0.875d0
  P(18,3) = 0.875d0
   P(19,1) = 0.875d0
   P(19,2) = 0.375d0
   P(19,3) = 0.375d0
  P(20,1) = 0.375d0
  P(20,2) = 0.875d0
  P(20,3) = 0.375d0
   P(21,1) = 0.625d0
   P(21,2) = 0.125d0
   P(21,3) = 0.625d0
  P(22,1) = 0.375d0
  P(22,2) = 0.375d0
  P(22,3) = 0.875d0
   P(23,1) = 0.125d0
   P(23,2) = 0.625d0
   P(23,3) = 0.625d0
  P(24,1) = 0.625d0
  P(24,2) = 0.625d0
  P(24,3) = 0.125d0
  P(:,1) = create_a0(1)*P(:,1)
  P(:,2) = create_a0(2)*P(:,2)
  P(:,3) = create_a0(3)*P(:,3)
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  DO i=2,16
    P(i,4) = P(1,4)
  ENDDO
  CALL ATOMNUMBER(create_species(2),P(17,4))
  DO i=18,24
    P(i,4) = P(17,4)
  ENDDO
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,2) = create_a0(2)
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up the messages
  WRITE(comment(1),*) TRIM(create_species(1))//"2"//TRIM(create_species(2))
  comment(1) = TRIM(ADJUSTL(comment(1)))//' with C15 Laves structure'
!
!
CASE('per','perovskite')
  cubic = .TRUE.
  IF(nspecies.NE.3) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 3.d0 /))
    GOTO 810
  ENDIF
  ALLOCATE(P(5,4))
  !Set up atom positions
  P(:,:) = 0.d0
  P(1,:) = 0.5d0
  P(3,1) = 0.5d0
   P(4,2) = 0.5d0
  P(5,3) = 0.5d0
  P(:,1) = create_a0(1)*P(:,1)
  P(:,2) = create_a0(2)*P(:,2)
  P(:,3) = create_a0(3)*P(:,3)
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  CALL ATOMNUMBER(create_species(2),P(2,4))
  CALL ATOMNUMBER(create_species(3),P(3,4))
  P(4,4) = P(3,4)
  P(5,4) = P(3,4)
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,2) = create_a0(2)
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up the messages
  temp = TRIM(create_species(1))//TRIM(create_species(2))//TRIM(create_species(3))//'3'
  WRITE(comment(1),*) 'Cubic perovskite '//TRIM(temp)
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  TETRAGONAL LATTICES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('st')
  IF(nspecies.NE.1) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0 /))
    GOTO 810
  ENDIF
  ALLOCATE(P(1,4))
  !Set up atom positions
  P(:,:) = 0.d0
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,2) = create_a0(2)
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up the messages
  WRITE(comment(1),*) 'Simple tetragonal '//TRIM(ADJUSTL((create_species(1))))
!
!
CASE('bct')
  IF(nspecies.NE.1 .AND. nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0,2.d0 /))
    GOTO 810
  ENDIF
  ALLOCATE(P(2,4))
  !Set up atom positions
  P(:,:) = 0.d0
  P(2,1) = 0.5d0
  P(2,2) = 0.5d0
  P(2,3) = 0.5d0
  P(:,1) = create_a0(1)*P(:,1)
  P(:,2) = create_a0(2)*P(:,2)
  P(:,3) = create_a0(3)*P(:,3)
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  IF(nspecies==2) THEN
    CALL ATOMNUMBER(create_species(2),P(2,4))
  ELSE
    P(2,4) = P(1,4)
  ENDIF
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,2) = create_a0(2)
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up the messages
  WRITE(comment(1),*) TRIM(create_species(1))
  IF(nspecies==1) THEN
    comment(1) = 'Tetragonal body-centered'//TRIM(ADJUSTL(comment(1)))
  ELSE
    comment(1) = 'Tetragonal body-centered'//TRIM(ADJUSTL(comment(1)))//TRIM(ADJUSTL(create_species(2)))//' alloy'
  ENDIF
!
!
CASE('fct','L10','L1_0')
  IF(nspecies.NE.1 .AND. nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0,2.d0 /))
    GOTO 810
  ENDIF
  ALLOCATE(P(4,4))
  !Set up atom positions
  P(:,:) = 0.d0
  P(2,1) = 0.5d0
  P(2,2) = 0.5d0
   P(3,2) = 0.5d0
   P(3,3) = 0.5d0
  P(4,1) = 0.5d0
  P(4,3) = 0.5d0
  P(:,1) = create_a0(1)*P(:,1)
  P(:,2) = create_a0(2)*P(:,2)
  P(:,3) = create_a0(3)*P(:,3)
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  P(:,4) = P(1,4)
  IF(nspecies==2) THEN
    CALL ATOMNUMBER(create_species(2),P(3,4))
  ELSE
    P(3,4) = P(1,4)
  ENDIF
  P(4,4) = P(3,4)
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,2) = create_a0(2)
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up the messages
  WRITE(comment(1),*) TRIM(create_species(1))
  IF(nspecies==1) THEN
    comment(1) = 'Tetragonal face-centered '//TRIM(ADJUSTL(comment(1)))
  ELSE
    comment(1) = 'Tetragonal face-centered'//TRIM(ADJUSTL(comment(1)))//TRIM(ADJUSTL(create_species(2)))//' alloy'
  ENDIF
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!  HEXAGONAL LATTICES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('hcp')
  hexagonal = .TRUE.
  IF(nspecies.NE.1 .AND. nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0,2.d0 /))
    GOTO 810
  ENDIF
  ALLOCATE(P(2,4))
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,1) = create_a0(2)*DCOS(DEG2RAD(120.d0))
  H(2,2) = create_a0(2)*DSIN(DEG2RAD(120.d0))
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up atom positions
  P(:,:) = 0.d0
  x = 1.d0/3.d0
  y = 2.d0/3.d0
  z1 = 0.5d0
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  IF(nspecies==2) THEN
    CALL ATOMNUMBER(create_species(2),P(2,4))
  ELSE
    P(2,4) = P(1,4)
  ENDIF
  !Transform atom positions to cartesian
  P(2,1) = x*H(1,1) + y*H(2,1)
  P(2,2) = y*H(2,2)
  P(2,3) = z1*H(3,3)
  !Set up the messages
  WRITE(comment(1),*) TRIM(create_species(1))
  IF(nspecies==1) THEN
    comment(1) = 'Hcp '//TRIM(ADJUSTL(comment(1)))
  ELSE
    comment(1) = 'Hcp '//TRIM(ADJUSTL(comment(1)))//TRIM(ADJUSTL(create_species(2)))//' alloy'
  ENDIF
!
CASE('wurtzite','wz')
  hexagonal = .TRUE.
  IF(nspecies.NE.1 .AND. nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0,2.d0 /))
    GOTO 810
  ENDIF
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,1) = create_a0(2)*DCOS(DEG2RAD(120.d0))
  H(2,2) = create_a0(2)*DSIN(DEG2RAD(120.d0))
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up atom positions
  x = 1.d0/3.d0
  y = 2.d0/3.d0
  z1 = 3.d0/8.d0
  ALLOCATE(P(4,4))
  P(:,:) = 0.d0
  P(1,1) = x*H(1,1) + y*H(2,1)
  P(1,2) = y*H(2,2)
  P(1,3) = 0.d0
   P(2,1) = y*H(1,1) + x*H(2,1)
   P(2,2) = x*H(2,2)
   P(2,3) = 0.5d0*H(3,3)
  P(3,1) = x*H(1,1) + y*H(2,1)
  P(3,2) = y*H(2,2)
  P(3,3) = z1*H(3,3)
   z1 = 7.d0/8.d0
   P(4,1) = y*H(1,1) + x*H(2,1)
   P(4,2) = x*H(2,2)
   P(4,3) = z1*H(3,3)
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  P(2:4,4) = P(1,4)
  IF(nspecies==2) THEN
    CALL ATOMNUMBER(create_species(2),P(3,4))
    P(4,4) = P(3,4)
    WRITE(comment(1),*) TRIM(create_species(1))//TRIM(create_species(2))
  ELSE
    WRITE(comment(1),*) TRIM(create_species(1))
  ENDIF
  comment(1) = TRIM(ADJUSTL(comment(1)))//' with wurtzite structure'
!
!
CASE('graphite')
  hexagonal = .TRUE.
  IF(nspecies.NE.1 .AND. nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0,2.d0 /))
    GOTO 810
  ENDIF
  !Set up the unit cell
  H(1,1) = create_a0(1)
    H(2,1) = create_a0(2)*DCOS(DEG2RAD(120.d0))
    H(2,2) = create_a0(2)*DSIN(DEG2RAD(120.d0))
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up atom positions
  ALLOCATE(P(4,4))
  P(:,:) = 0.d0
   P(2,3) = 0.5d0*H(3,3)
  x = 1.d0/3.d0
  y = 2.d0/3.d0
  P(3,:) = x*H(1,:) + y*H(2,:)
   x = 2.d0/3.d0
   y = 1.d0/3.d0
   P(4,:) = x*H(1,:) + y*H(2,:) + 0.5*H(3,:)
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  IF(nspecies==2) THEN
    CALL ATOMNUMBER(create_species(2),P(2,4))
  ELSE
    P(2,4) = P(1,4)
  ENDIF
  P(3,4) = P(1,4)
  P(4,4) = P(2,4)
  !Set up the messages
  temp = TRIM(create_species(1))//TRIM(create_species(2))
  WRITE(comment(1),*) TRIM(temp)//' with hexagonal graphite structure'
!
!
CASE('c14','C14')
  hexagonal = .TRUE.
  IF(nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0,2.d0 /))
    GOTO 810
  ENDIF
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,1) = create_a0(2)*DCOS(DEG2RAD(120.d0))
  H(2,2) = create_a0(2)*DSIN(DEG2RAD(120.d0))
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up atom positions
  !Consider the prototype MgZn2
  ALLOCATE(P(12,4))
  P(:,:) = 0.d0
  !Positions of Mg
  P(1,2) = 0.667d0
  P(1,3) = 0.063d0
  P(2,2) = 0.667d0
  P(2,3) = 0.437d0
  P(3,1) = 0.500d0
  P(3,2) = 0.333d0
  P(3,3) = 0.937d0
  P(4,1) = 0.500d0
  P(4,2) = 0.333d0
  P(4,3) = 0.563d0
  !Positions of Zn
  P(6,3) = 0.500d0
  P(7,2) = 0.338d0
  P(7,3) = 0.750d0
  P(8,1) = -0.247d0
  P(8,2) = 0.831d0
  P(8,3) = 0.750d0
  P(9,1) = 0.253d0
  P(9,2) = 0.169d0
  P(9,3) = 0.250d0
  P(10,1) = 0.247d0
  P(10,2) = 0.831d0
  P(10,3) = 0.750d0
  P(11,1) = 0.747d0
  P(11,2) = 0.167d0
  P(11,3) = 0.250d0
  P(12,1) = 0.500d0
  P(12,2) = 0.662d0
  P(12,3) = 0.250d0
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  P(2,4) = P(1,4)
  P(3,4) = P(1,4)
  P(4,4) = P(1,4)
  CALL ATOMNUMBER(create_species(2),P(5,4))
  P(6,4) = P(5,4)
  P(7,4) = P(5,4)
  P(8,4) = P(5,4)
  P(9,4) = P(5,4)
  P(10,4) = P(5,4)
  P(11,4) = P(5,4)
  P(12,4) = P(5,4)
  !Transform atom positions to cartesian
  P(:,1) = H(1,1)*P(:,1)
  P(:,2) = H(2,2)*P(:,2)
  P(:,3) = H(3,3)*P(:,3)
  !Set up the messages
  WRITE(comment(1),*) TRIM(create_species(1))//TRIM(create_species(2))//"2"
  comment(1) = TRIM(ADJUSTL(comment(1)))//' with C14 Laves structure'
!
!
CASE('c36','C36')
  hexagonal = .TRUE.
  IF(nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0,2.d0 /))
    GOTO 810
  ENDIF
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,1) = create_a0(2)*DCOS(DEG2RAD(120.d0))
  H(2,2) = create_a0(2)*DSIN(DEG2RAD(120.d0))
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up atom positions
  !Consider the prototype MgNi2
  ALLOCATE(P(24,4))
  P(:,:) = 0.d0
  !Positions of Mg
  P(1,3) = 0.406d0
  P(2,3) = 0.594d0
  P(3,3) = 0.906d0
  P(4,3) = 0.094d0
  P(5,2) = 0.667d0
  P(5,3) = 0.656d0
  P(6,1) = 0.500d0
  P(6,2) = 0.333d0
  P(6,3) = 0.344d0
  P(7,1) = 0.500d0
  P(7,2) = 0.333d0
  P(7,3) = 0.156d0
  P(8,2) = 0.667d0
  P(8,3) = 0.844d0
  !Positions of Ni
  P(9,2) = 0.667d0
  P(9,3) = 0.375d0
  P(10,1) = 0.500d0
  P(10,2) = 0.333d0
  P(10,3) = 0.625d0
  P(11,1) = 0.500d0
  P(11,2) = 0.333d0
  P(11,3) = 0.875d0
  P(12,2) = 0.667d0
  P(12,3) = 0.125d0
  P(13,1) = 0.250d0
  P(13,2) = 0.500d0
  P(13,3) = 0.500d0
  P(14,1) = 0.500d0
  P(14,3) = 0.500d0
  P(15,1) = -0.250d0
  P(15,2) = 0.500d0
  P(15,3) = 0.500d0
  P(16,1) = 0.250d0
  P(16,2) = 0.500d0
  P(17,1) = 0.500d0
  P(18,1) = -0.250d0
  P(18,2) = 0.500d0
  P(19,1) = -0.254d0
  P(19,2) = 0.836d0
  P(19,3) = 0.250d0
  P(20,2) = 0.329d0
  P(20,3) = 0.250d0
  P(21,1) = 0.254d0
  P(21,2) = 0.836d0
  P(21,3) = 0.250d0
  P(22,1) = 0.754d0
  P(22,2) = 0.164d0
  P(22,3) = 0.750d0
  P(23,1) = 0.500d0
  P(23,2) = 0.671d0
  P(23,3) = 0.750d0
  P(24,1) = 0.246d0
  P(24,2) = 0.164d0
  P(24,3) = 0.750d0
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  P(2,4) = P(1,4)
  P(3,4) = P(1,4)
  P(4,4) = P(1,4)
  P(5,4) = P(1,4)
  P(6,4) = P(1,4)
  P(7,4) = P(1,4)
  P(8,4) = P(1,4)
  CALL ATOMNUMBER(create_species(2),P(9,4))
  P(10,4) = P(9,4)
  P(11,4) = P(9,4)
  P(12,4) = P(9,4)
  P(13,4) = P(9,4)
  P(14,4) = P(9,4)  
  P(15,4) = P(9,4)
  P(16,4) = P(9,4)
  P(17,4) = P(9,4)
  P(18,4) = P(9,4)
  P(21,4) = P(9,4)
  P(20,4) = P(9,4)
  P(21,4) = P(9,4)
  P(22,4) = P(9,4)
  P(23,4) = P(9,4)
  P(24,4) = P(9,4)
  !Transform atom positions to cartesian
  P(:,1) = H(1,1)*P(:,1)
  P(:,2) = H(2,2)*P(:,2)
  P(:,3) = H(3,3)*P(:,3)
  !Set up the messages
  WRITE(comment(1),*) TRIM(create_species(1))//TRIM(create_species(2))//"2"
  comment(1) = TRIM(ADJUSTL(comment(1)))//' with C36 Laves structure'
!
!
CASE('limo2','LiMO2')
  hexagonal = .TRUE.
  IF(nspecies.NE.3) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 3.d0 /))
    GOTO 810
  ENDIF
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,1) = create_a0(2)*DCOS(DEG2RAD(120.d0))
  H(2,2) = create_a0(2)*DSIN(DEG2RAD(120.d0))
  H(3,3) = create_a0(3)
  Huc(:,:) = H(:,:)
  !Set up atom positions
  !Consider the prototype LiMO2
  ALLOCATE(P(12,4))
  P(:,:) = 0.d0
  x = 1.d0/3.d0
  y = 2.d0/3.d0
  !Positions of Li
  !   0     0     0
  ! 2/3   1/3   1/3
  ! 1/3   2/3   2/3
  P(2,1) = y*H(1,1) + x*H(2,1)  
  P(2,2) = x*H(2,2)
  P(2,3) = x*H(3,3)
  P(3,1) = x*H(1,1) + y*H(2,1)
  P(3,2) = y*H(2,2)
  P(3,3) = y*H(3,3)
  !Positions of M
  !   0     0     1/5
  ! 2/3   1/3   0.833
  ! 1/3   2/3   0.167
  P(4,3) = 0.50d0*H(3,3)
  P(5,1) = y*H(1,1) + x*H(2,1)  
  P(5,2) = x*H(2,2)
  P(5,3) = 0.833d0*H(3,3) 
  P(6,1) = x*H(1,1) + y*H(2,1)
  P(6,2) = y*H(2,2)
  P(6,3) = 0.167d0*H(3,3)
  !Positions of O
  !   0     0   0.231
  !   0     0   0.769
  ! 2/3   1/3   0.564
  ! 2/3   1/3   0.102
  ! 1/3   2/3   0.898
  ! 1/3   2/3   0.436
  P(7,3) = 0.231d0*H(3,3)
  P(8,3) = 0.769d0*H(3,3)
  P(9,1) = y*H(1,1) + x*H(2,1)
  P(9,2) = x*H(2,2)
  P(9,3) = 0.564d0*H(3,3)
  P(10,1) = y*H(1,1) + x*H(2,1)
  P(10,2) = x*H(2,2)
  P(10,3) = 0.102d0*H(3,3)
  P(11,1) = x*H(1,1) + y*H(2,1)
  P(11,2) = y*H(2,2)
  P(11,3) = 0.898d0*H(3,3)
  P(12,1) = x*H(1,1) + y*H(2,1)
  P(12,2) = y*H(2,2)
  P(12,3) = 0.436d0*H(3,3)
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),P(1,4))
  P(2,4) = P(1,4)
  P(3,4) = P(1,4)
  CALL ATOMNUMBER(create_species(2),P(4,4))
  P(5,4) = P(4,4)
  P(6,4) = P(4,4)
  CALL ATOMNUMBER(create_species(3),P(7,4))
  P(8,4) = P(7,4)
  P(9,4) = P(7,4)
  P(10,4) = P(7,4)
  P(11,4) = P(7,4)
  P(12,4) = P(7,4)
  !Set up the messages
  temp = TRIM(create_species(1))//TRIM(create_species(2))//TRIM(create_species(3))//'2'
  WRITE(comment(1),*) TRIM(temp)//' with hexagonal structure'
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  OTHER STRUCTURES  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE('nanotube','NT','nt')
  IF(nspecies.NE.1 .AND. nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0,2.d0 /))
    GOTO 810
  ENDIF
  IF(NT_mn(1)<=0 .AND. NT_mn(2)<=0) THEN
    nerr = nerr+1
    CALL ATOMSK_MSG(4802,(/''/),(/0.d0/))
    GOTO 1000
  ENDIF
  !Make sure the first index is the smallest
  IF(NT_mn(1)>NT_mn(2)) THEN
    r = NT_mn(2)
    NT_mn(2) = NT_mn(1)
    NT_mn(1) = r
  ENDIF
  !Set type of nanotube
  WRITE(NT_type,*) NT_mn(1)
  WRITE(temp,*) NT_mn(2)
  NT_type = '('//TRIM(ADJUSTL(NT_type))//','//TRIM(ADJUSTL(temp))//')'
  IF(NT_mn(1)==NT_mn(2)) THEN
    NT_type = 'armchair '//TRIM(NT_type)
  ELSEIF( NT_mn(1)==0 .OR. NT_mn(2)==0 ) THEN
    NT_type = 'zigzag '//TRIM(NT_type)
  ELSE
    NT_type = 'chiral '//TRIM(NT_type)
  ENDIF
  !
  ! Definition of lattice vectors for the graphene-like sheet
  H(1,1) = DSQRT(3.d0)*create_a0(1)/2.d0
  H(1,2) = create_a0(2)/2.d0
   H(2,1) = DSQRT(3.d0)*create_a0(1)/2.d0
   H(2,2) = -create_a0(2)/2.d0
  H(3,3) = 1.d0  !Arbitrary
  Huc(:,:) = H(:,:)
  !Positions of the 2 atoms of the graphene-like sheet
  ALLOCATE(Q(2,4))
  Q(:,:) = 0.d0
  Q(2,1) = create_a0(1)*DSQRT(3.d0)/3.d0
  !Set up atom species
  CALL ATOMNUMBER(create_species(1),Q(1,4))
  IF(nspecies==2) THEN
    CALL ATOMNUMBER(create_species(2),Q(2,4))
  ELSE
    Q(2,4) = Q(1,4)
  ENDIF
  !
  Hrecip = CROSS_PRODUCT(H(2,:), H(3,:))
  vol_cell = DABS(DOT_PRODUCT(H(1,:), Hrecip))
  !
  !Calculate coordinates of nanotube basis vectors
  d = GCD(2*NT_mn(2)+NT_mn(1), 2*NT_mn(1)+NT_mn(2))
  a_perp(:) = DBLE(NT_mn(2))*H(1,:) + DBLE(NT_mn(1))*H(2,:)
  a_para(:) = (( 2.d0*DBLE(NT_mn(1))+DBLE(NT_mn(2)) )/DBLE(d))*H(1,:) -    &
            & (( 2.d0*DBLE(NT_mn(2))+DBLE(NT_mn(1)) )/DBLE(d))*H(2,:)
  !l_perp, l_para = sides of the rectangle to cut from the sheet
  l_perp = DSQRT(DOT_PRODUCT(a_perp,a_perp))
  l_para = DSQRT(DOT_PRODUCT(a_para,a_para))
  !Compute reciprocal base vectors
  b_perp(:) = CROSS_PRODUCT(a_para(:),H(3,:))
  b_para(:) = CROSS_PRODUCT(H(3,:),a_perp(:))
  !Compute volume of unit cell for nanotube
  vol = DABS(DOT_PRODUCT(a_perp,b_perp))
  b_perp(:) = 2.d0*pi*b_perp(:)/vol
  b_para(:) = 2.d0*pi*b_para(:)/vol
  !
  !Calculate radius of nanotube
  NT_radius = l_perp/(2.d0*pi)
  !Calculate number of atoms for this nanotube
  NT_NP = FLOOR(2.d0*(vol/vol_cell)+0.1d0)
  WRITE(msg,*) '2*(vol/vol_cell) = ', 2.d0*(vol/vol_cell)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) 'Theoretical number of atoms: ', NT_NP
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  ALLOCATE(P(NT_NP,4))
  P(:,:) = 0.d0
  !
  !Compute positions of atoms that will form the nanotube
  NT_NP=0
  DO r=-1000, 1000
    DO t=-1000, 1000
      DO i=1,2 !Loop on the 2 atoms of the unit cell
        !Set coordinates of atom
        DO j=1,3
          coord(j) = Q(i,j) + DBLE(r)*H(1,j) + DBLE(t)*H(2,j)
        ENDDO
        x = DOT_PRODUCT(coord, a_perp)/l_perp
        y = DOT_PRODUCT(coord, a_para)/l_para
        !Check if coordinates are inside the cell
        IF( x>=0.d0 .AND. x<=l_perp .AND. y>=0.d0 .AND. y<=l_para ) THEN
          new=.TRUE.
          !Check if current atom is actually a new one or a periodic replica
          IF(NT_NP>0) THEN
            DO j=1, NT_NP
              z1 = DOT_PRODUCT(coord(:)-P(j,1:3), b_perp)/(2.d0*pi)
              z2 = DOT_PRODUCT(coord(:)-P(j,1:3), b_para)/(2.d0*pi)
              IF( DABS(z1-DBLE(NINT(z1)))<1.d-12 .AND.                   &
                & DABS(z2-DBLE(NINT(z2)))<1.d-12       ) THEN
                new=.FALSE.
                EXIT
              ENDIF
            ENDDO
          ENDIF
          !Not a replica? Then save it to P
          IF(new) THEN
            NT_NP = NT_NP+1
            IF( NT_NP>SIZE(P(:,1)) ) GOTO 820
            P(NT_NP,1:3) = coord(:)
            P(NT_NP,4) = Q(i,4)
          END IF
        END IF
      END DO ! end loop over i
    END DO ! end loop over t
  END DO !end loop over r
  !
  DEALLOCATE(Q)
  !
  IF(verbosity==4) THEN
    WRITE(msg,*) 'Number of atoms found: ', NT_NP
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    temp = 'atomsk_NTsheet.xyz'
    CALL WRITE_XYZ(H,P,comment,AUXNAMES,AUX,temp)
  ENDIF
  IF( NT_NP.NE.SIZE(P(:,1)) ) GOTO 820
  !
  !Now let's roll the nanotube
  DO i=1,NT_NP
    DO j=1,3
      coord(j) = P(i,j)
    ENDDO
    x = DOT_PRODUCT(coord(:),a_perp(:)) / (l_perp*NT_radius)
    P(i,1) = NT_radius*DCOS(x)
    P(i,2) = NT_radius*DSIN(x)
    P(i,3) = DOT_PRODUCT(coord(:),a_para) / l_para
  ENDDO
  !
  !Modify supercell vectors to fit the nanotube
  H(:,:) = 0.d0
  H(1,1) = 4.d0*NT_radius
   H(2,1) = 0.5d0*H(1,1)
   H(2,2) = 0.5d0*H(1,1)*DSQRT(3.d0)
  H(3,3) = l_para
  !
  !Set up the messages
  temp = TRIM(create_species(1))//TRIM(create_species(2))
  WRITE(comment(1),*) TRIM(temp)//' nanotube of type '//NT_type
  WRITE(temp,'(f16.3)') NT_radius
  comment(1) = TRIM(ADJUSTL(comment(1)))//'; radius: '//TRIM(ADJUSTL(temp))//' A.'
!
!
! -- add other structures according to the groups defined above (cubic/hexagonal/etc.) --
!
CASE DEFAULT
  !Structure not recogized
  nerr = nerr+1
  CALL ATOMSK_MSG(4805,(/create_struc/),(/0.d0/))
  GOTO 1000
!
END SELECT
!
!
!
300 CONTINUE
IF( verbosity==4 ) THEN
  IF( cubic ) THEN
    WRITE(msg,*) "System is cubic:", cubic
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  ENDIF
  IF( hexagonal ) THEN
    WRITE(msg,*) "System is hexagonal:", hexagonal
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  ENDIF
ENDIF
!
!Orient the system according to ORIENT(:,:) - CUBIC SYSTEMS ONLY!
IF( cubic ) THEN
  IF( VECLENGTH(ORIENT(1,:))>1.d-12 .OR. VECLENGTH(ORIENT(2,:))>1.d-12 &
    & .OR. VECLENGTH(ORIENT(3,:))>1.d-12 ) THEN
    !
    !Check that no vector in ORIENT is [000]
    IF( VECLENGTH(ORIENT(1,:))<1.d-12 .OR. VECLENGTH(ORIENT(2,:))<1.d-12 &
    & .OR. VECLENGTH(ORIENT(3,:))<1.d-12 ) THEN
      CALL ATOMSK_MSG(814,(/""/),(/0.d0/))
      nerr=nerr+1
      GOTO 1000
    ENDIF
    !
    !Reduce Miller indices
    !For example direction [333] will be replaced by [111],
    !direction [840] replaced by [210], and so on
    DO i=1,3
      u = ORIENT(i,1)
      v = ORIENT(i,2)
      w = ORIENT(i,3)
      IF( DABS(u)>0.1d0 .AND. NINT(DABS(v))>0.1d0 ) THEN
        z1 = GCD( NINT(DABS(u)) , NINT(DABS(v)) )
      ELSE
        z1 = MAX(DABS(u),DABS(v))
      ENDIF
      IF( DABS(u)>0.1d0 .AND. NINT(DABS(w))>0.1d0 ) THEN
        z2 = GCD( NINT(DABS(u)) , NINT(DABS(w)) )
      ELSE
        z2 = MAX(DABS(u),DABS(w))
      ENDIF
      IF( DABS(z1)>0.1d0 .AND. NINT(z2)>0.1d0 ) THEN
        x = GCD( NINT(DABS(z1)),NINT(DABS(z2)) )
      ELSE  !i.e. z1==0 or z2==0
        x = MAX( DABS(z1) , DABS(z2) )
      ENDIF
      IF( DABS(x)<0.1d0 ) x=1.d0  !avoid division by zero
      !Set box vector
      ORIENT(i,1) = u / x
      ORIENT(i,2) = v / x
      ORIENT(i,3) = w / x
    ENDDO
    !
    !Generate text with target orientation
    oriented = .TRUE.
    comment(1) = TRIM(comment(1))//' oriented'
    DO i=1,3
      IF(i==1) comment(1) = TRIM(comment(1))//' X=['
      IF(i==2) comment(1) = TRIM(comment(1))//' Y=['
      IF(i==3) comment(1) = TRIM(comment(1))//' Z=['
      DO j=1,3
        m = NINT(ORIENT(i,j))
        WRITE(msg,*) m
        IF( j>1 .AND. ANY(ABS(NINT(ORIENT(1,:)))>=10) ) THEN
          comment(1) = TRIM(comment(1))//"_"//TRIM(ADJUSTL(msg))
        ELSE
          comment(1) = TRIM(comment(1))//TRIM(ADJUSTL(msg))
        ENDIF
        IF(j==3) comment(1) = TRIM(comment(1))//']'
      ENDDO
    ENDDO
    !
    WRITE(msg,*) "orienting the system:"
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,'(3f16.6)') ORIENT(1,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,'(3f16.6)') ORIENT(2,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,'(3f16.6)') ORIENT(3,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
    !Set lminmax = 10 * largest value in ORIENT
    lminmax = 10 * NINT(MAXVAL(DABS(ORIENT(:,:))))
    !
    !Check that vectors of ORIENT are orthogonal
    orthovec(:) = .FALSE.
    IF( DABS( ANGVEC(ORIENT(1,:),ORIENT(2,:))-(pi/2.d0) ) < 1.d-6 ) THEN
      orthovec(1)=.TRUE.
    ELSE
      IF( DABS( ANGVEC(ORIENT(3,:),ORIENT(1,:))-(pi/2.d0) ) < 1.d-6 ) THEN
        !Vectors 1 and 3 are orthogonal, but vector 2 is not
        !Compute what vector 2 should be, and suggest it to the user
        a_perp(:) = CROSS_PRODUCT( ORIENT(3,:) , ORIENT(1,:) )
        msg = '['
        DO i=1,3
          WRITE(temp,*) NINT(a_perp(i))
          msg = TRIM(ADJUSTL(msg))//TRIM(ADJUSTL(temp))
        ENDDO
        msg = TRIM(ADJUSTL(msg))//']'
        temp = 'Y'
        CALL ATOMSK_MSG(4819,(/msg,temp/),(/0.d0/))
      ELSEIF( DABS( ANGVEC(ORIENT(2,:),ORIENT(3,:))-(pi/2.d0) ) < 1.d-6 ) THEN
        !Vectors 2 and 3 are orthogonal, but vector 1 is not
        !Compute what vector 1 should be, and suggest it to the user
        a_perp(:) = CROSS_PRODUCT( ORIENT(2,:) , ORIENT(3,:) )
        msg = '['
        DO i=1,3
          WRITE(temp,*) NINT(a_perp(i))
          msg = TRIM(ADJUSTL(msg))//TRIM(ADJUSTL(temp))
        ENDDO
        msg = TRIM(ADJUSTL(msg))//']'
        temp = 'X'
        CALL ATOMSK_MSG(4819,(/msg,temp/),(/0.d0/))
      ELSE
        !None of the 3 vectors are orthogonal
        CALL ATOMSK_MSG(4819,(/''/),(/0.d0/))
      ENDIF
      nerr = nerr+1
      GOTO 1000
    ENDIF
    IF( DABS( ANGVEC(ORIENT(2,:),ORIENT(3,:))-(pi/2.d0) ) < 1.d-6 ) THEN
      orthovec(2)=.TRUE.
    ELSE
      IF( DABS( ANGVEC(ORIENT(1,:),ORIENT(2,:))-(pi/2.d0) ) < 1.d-6 ) THEN
        !Vectors 1 and 2 are orthogonal, but vector 3 is not
        !Compute what vector 3 should be, and suggest it to the user
        a_perp(:) = CROSS_PRODUCT( ORIENT(1,:) , ORIENT(2,:) )
        msg = '['
        DO i=1,3
          WRITE(temp,*) NINT(a_perp(i))
          msg = TRIM(ADJUSTL(msg))//TRIM(ADJUSTL(temp))
        ENDDO
        msg = TRIM(ADJUSTL(msg))//']'
        temp = 'Z'
        CALL ATOMSK_MSG(4819,(/msg,temp/),(/0.d0/))
      ELSEIF( DABS( ANGVEC(ORIENT(3,:),ORIENT(1,:))-(pi/2.d0) ) < 1.d-6 ) THEN
        !Vectors 1 and 3 are orthogonal, but vector 2 is not
        !Compute what vector 2 should be, and suggest it to the user
        a_perp(:) = CROSS_PRODUCT( ORIENT(3,:) , ORIENT(1,:) )
        msg = '['
        DO i=1,3
          WRITE(temp,*) NINT(a_perp(i))
          msg = TRIM(ADJUSTL(msg))//TRIM(ADJUSTL(temp))
        ENDDO
        msg = TRIM(ADJUSTL(msg))//']'
        temp = 'Y'
        CALL ATOMSK_MSG(4819,(/msg,temp/),(/0.d0/))
      ELSE
        !None of the 3 vectors are orthogonal
        CALL ATOMSK_MSG(4819,(/''/),(/0.d0/))
      ENDIF
      nerr = nerr+1
      GOTO 1000
    ENDIF
    IF( DABS( ANGVEC(ORIENT(3,:),ORIENT(1,:))-(pi/2.d0) ) < 1.d-6 ) THEN
      orthovec(3)=.TRUE.
    ELSE
      IF( DABS( ANGVEC(ORIENT(1,:),ORIENT(2,:))-(pi/2.d0) ) < 1.d-6 ) THEN
        !Vectors 1 and 2 are orthogonal, but vector 3 is not
        !Compute what vector 3 should be, and suggest it to the user
        a_perp(:) = CROSS_PRODUCT( ORIENT(1,:) , ORIENT(2,:) )
        msg = '['
        DO i=1,3
          WRITE(temp,*) NINT(a_perp(i))
          msg = TRIM(ADJUSTL(msg))//TRIM(ADJUSTL(temp))
        ENDDO
        msg = TRIM(ADJUSTL(msg))//']'
        temp = 'Z'
        CALL ATOMSK_MSG(4819,(/msg,temp/),(/0.d0/))
      ELSEIF( DABS( ANGVEC(ORIENT(2,:),ORIENT(3,:))-(pi/2.d0) ) < 1.d-6 ) THEN
        !Vectors 2 and 3 are orthogonal, but vector 1 is not
        !Compute what vector 1 should be, and suggest it to the user
        a_perp(:) = CROSS_PRODUCT( ORIENT(2,:) , ORIENT(3,:) )
        msg = '['
        DO i=1,3
          WRITE(temp,*) NINT(a_perp(i))
          msg = TRIM(ADJUSTL(msg))//TRIM(ADJUSTL(temp))
        ENDDO
        msg = TRIM(ADJUSTL(msg))//']'
        temp = 'X'
        CALL ATOMSK_MSG(4819,(/msg,temp/),(/0.d0/))
      ELSE
        !None of the 3 vectors are orthogonal
        CALL ATOMSK_MSG(4819,(/''/),(/0.d0/))
      ENDIF
      nerr = nerr+1
      GOTO 1000
    ENDIF
    !
    !Currently, unit cell is oriented so that H(1,:)//[100], H(2,:)//[010], H(3,:)//[001]
    !=> express the vectors H(:,:) and atom positions P(:,:) in the new base ORIENT
    !Normalize the vectors in ORIENT, save them to ORIENTN
    DO i=1,3
      ORIENTN(i,:) = ORIENT(i,:)/VECLENGTH(ORIENT(i,:))
    ENDDO
    !Convert atom positions to reduced coordinates
    CALL CART2FRAC(P,H)
    !Rotate vectors in H(:,:)
    DO i=1,3
      H1 = H(i,1)
      H2 = H(i,2)
      H3 = H(i,3)
      H(i,1) = H1*ORIENTN(1,1) + H2*ORIENTN(1,2) + H3*ORIENTN(1,3)
      H(i,2) = H1*ORIENTN(2,1) + H2*ORIENTN(2,2) + H3*ORIENTN(2,3)
      H(i,3) = H1*ORIENTN(3,1) + H2*ORIENTN(3,2) + H3*ORIENTN(3,3)
    ENDDO
    !Convert atom positions back to cartesian coordinates
    CALL FRAC2CART(P,H)
    WRITE(msg,*) "oriented H(:,:):"
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,'(3f16.6)') H(1,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,'(3f16.6)') H(2,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,'(3f16.6)') H(3,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
    !The oriented unit cell vectors are defined by the Miller indices
    uv(:,:) = 0.d0
    uv(1,1) = ORIENT(1,1)*H(1,1) + ORIENT(1,2)*H(2,1) + ORIENT(1,3)*H(3,1)
    uv(2,2) = ORIENT(2,1)*H(1,2) + ORIENT(2,2)*H(2,2) + ORIENT(2,3)*H(3,2)
    uv(3,3) = ORIENT(3,1)*H(1,3) + ORIENT(3,2)*H(2,3) + ORIENT(3,3)*H(3,3)
    WRITE(msg,*) "Oriented unit cell vectors:"
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,'(3f16.6)') uv(1,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,'(3f16.6)') uv(2,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,'(3f16.6)') uv(3,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
    !Correct cell length for some special orientations to find the minimal repetition unit
    SELECT CASE(create_struc)
    CASE('bcc','BCC')
      !in bcc structure, shortest period along a <hkl> direction is actually 1/2<hkl> if h, k, l are all odd
      !Examples: [111] is replaced by 1/2[111], [531] becomes 1/2[531], etc.
      DO i=1,3
        IF( MOD(NINT(DABS(ORIENT(i,1))),2).NE.0 .AND. MOD(NINT(DABS(ORIENT(i,2))),2).NE.0  &
          & .AND. MOD(NINT(DABS(ORIENT(i,3))),2).NE.0 ) THEN
          uv(i,i) = uv(i,i)/2.d0
        ENDIF
      ENDDO
      !
    CASE('fcc','FCC','diamond','dia','zincblende','zb','rocksalt','rs')
      !in fcc, diamond/zb, and rocksalt structures,
      !shortest period along a <hkl> direction is actually 1/2<hkl> if h and k are odd and l is even
      !Examples: [110] is replaced by 1/2[110], [112] becomes 1/2[112], etc.
      DO i=1,3
        IF(  MOD(NINT(DABS(ORIENT(i,1))),2).NE.0 .AND. MOD(NINT(DABS(ORIENT(i,2))),2).NE.0  &
          &  .AND. MOD(NINT(DABS(ORIENT(i,3))),2)==0 .OR.                                   &
          &  MOD(NINT(DABS(ORIENT(i,1))),2).NE.0 .AND. MOD(NINT(DABS(ORIENT(i,3))),2).NE.0  &
          &  .AND. MOD(NINT(DABS(ORIENT(i,2))),2)==0 .OR.                                   &
          &  MOD(NINT(DABS(ORIENT(i,2))),2).NE.0 .AND. MOD(NINT(DABS(ORIENT(i,3))),2).NE.0  &
          &  .AND. MOD(NINT(DABS(ORIENT(i,1))),2)==0   ) THEN
          uv(i,i) = uv(i,i)/2.d0
        ENDIF
      ENDDO
      !
    CASE DEFAULT
      !
    END SELECT
    !
    !For each atom in the unit cell H(:,:), keep only periodic replica that are inside the uv(:)
    !and store it in Q(:,:)
    WRITE(msg,*) "Duplicating atoms inside oriented unit cell..."
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !Estimate new number of particles NP by comparing volumes of old and new unit cells
    NP = 1.2d0*CEILING( SIZE(P,1) * DABS( DABS(uv(1,1)*uv(2,2)*uv(3,3)) / &
        & DABS(VECLENGTH(H(1,:))*VECLENGTH(H(2,:))*VECLENGTH(H(3,:))) ) ) + 20
    WRITE(msg,*) "Estimated new number of atoms : ", NP
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    IF(ALLOCATED(Q)) DEALLOCATE(Q)
    ALLOCATE( Q(NP,4) )
    Q(:,:) = 0.d0
    !
    !Loop over all replica in a wide range
    NP = 0
    DO i=1,SIZE(P,1)
      DO l=-lminmax,lminmax
        DO k=-lminmax,lminmax
          DO j=-lminmax,lminmax
            !Compute cartesian position of this replica
            tempP(1,1) = P(i,1) + DBLE(j)*H(1,1) + DBLE(k)*H(2,1) + DBLE(l)*H(3,1)
            tempP(1,2) = P(i,2) + DBLE(j)*H(1,2) + DBLE(k)*H(2,2) + DBLE(l)*H(3,2)
            tempP(1,3) = P(i,3) + DBLE(j)*H(1,3) + DBLE(k)*H(2,3) + DBLE(l)*H(3,3)
            tempP(1,4) = P(i,4)
            IF( tempP(1,1)>=-1.d-12 .AND. tempP(1,1)<uv(1,1)-1.d-12 .AND.             &
              & tempP(1,2)>=-1.d-12 .AND. tempP(1,2)<uv(2,2)-1.d-12 .AND.             &
              & tempP(1,3)>=-1.d-12 .AND. tempP(1,3)<uv(3,3)-1.d-12       ) THEN
              !This replica is inside the new cell, mark it as new
              new = .TRUE.
              !Verify that its position is different from all previous atoms
              DO m=1,NP
                IF( DABS( VECLENGTH(tempP(1,1:3)-Q(m,1:3)) )<1.d-6 ) THEN
                  new = .FALSE.
                ENDIF
              ENDDO
              !
              IF( new ) THEN
                NP = NP+1
                IF(NP>SIZE(Q,1)) THEN
                  nerr = nerr+1
                  CALL ATOMSK_MSG(4821,(/""/),(/DBLE(NP),DBLE(SIZE(Q,1))/))
                  GOTO 1000
                ENDIF
                Q(NP,:) = tempP(1,:)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    !Replace old P with the new Q
    IF(ALLOCATED(P)) DEALLOCATE(P)
    ALLOCATE(P(NP,4))
    DO i=1,NP
      P(i,:) = Q(i,:)
    ENDDO
    IF(ALLOCATED(Q)) DEALLOCATE(Q)
    WRITE(msg,*) "new NP in oriented cell:", NP
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
    !Replace old H with the new (oriented) cell vectors
    H(:,:) = 0.d0
    DO i=1,3
      H(i,i) = uv(i,i)
    ENDDO
    !
  ELSE
    comment(1) = TRIM(comment(1))//" oriented X=[100] Y=[010] Z=[001]"
    oriented = .FALSE. !set to .FALSE. to avoid duplicating atoms later
  ENDIF  !end if oriented
    !
    !
ELSEIF( hexagonal ) THEN
  IF( VECLENGTH(ORIENT(1,:))>1.d-12  .OR. VECLENGTH(ORIENT(2,:))>1.d-12 &
    & .OR. VECLENGTH(ORIENT(3,:))>1.d-12 ) THEN
    !
    !Check that no vector in ORIENT is [000]
    IF( VECLENGTH(ORIENT(1,:))<1.d-12 .OR. VECLENGTH(ORIENT(2,:))<1.d-12 &
    & .OR. VECLENGTH(ORIENT(3,:))<1.d-12 ) THEN
      CALL ATOMSK_MSG(814,(/""/),(/0.d0/))
      nerr=nerr+1
      GOTO 1000
    ENDIF
    !
    !Check that no box vector is a linear combination of the other two
    tempP(1,:) = SCALAR_TRIPLE_PRODUCT( ORIENT(1,:) , ORIENT(2,:) , ORIENT(3,:) )
    IF( VECLENGTH(tempP(1,:))<1.d-12 ) THEN
      nerr = nerr+1
      CALL ATOMSK_MSG(4829,(/""/),(/0.d0/))
      GOTO 1000
    ENDIF
    !
    oriented = .TRUE.
    comment(1) = TRIM(comment(1))//' with box vectors'
    DO i=1,3
      IF(i==1) comment(1) = TRIM(comment(1))//' H1=['
      IF(i==2) comment(1) = TRIM(comment(1))//' H2=['
      IF(i==3) comment(1) = TRIM(comment(1))//' H3=['
      DO j=1,3
        m = NINT(ORIENT(i,j))
        WRITE(msg,*) m
        IF( j>1 .AND. ( ANY(ABS(NINT(ORIENT(i,:)))>=10) .OR. ABS(NINT(ORIENT(i,1)+ORIENT(i,2)))>=10 ) ) THEN
          comment(1) = TRIM(comment(1))//"_"//TRIM(ADJUSTL(msg))
        ELSE
          comment(1) = TRIM(comment(1))//TRIM(ADJUSTL(msg))
        ENDIF
        IF( j==2 ) THEN
          m = -1*NINT(ORIENT(i,1)) - NINT(ORIENT(i,2))
          WRITE(msg,*) m
          IF( j>1 .AND. ( ANY(ABS(NINT(ORIENT(i,:)))>=10) .OR. ABS(NINT(ORIENT(i,1)+ORIENT(i,2)))>=10 ) ) THEN
            comment(1) = TRIM(comment(1))//"_"//TRIM(ADJUSTL(msg))
          ELSE
            comment(1) = TRIM(comment(1))//TRIM(ADJUSTL(msg))
          ENDIF
        ENDIF
        IF(j==3) comment(1) = TRIM(comment(1))//']'
      ENDDO
    ENDDO
    !
    WRITE(msg,*) "orienting the system:"
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,'(3f16.6)') ORIENT(1,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,'(3f16.6)') ORIENT(2,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,'(3f16.6)') ORIENT(3,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
    !Set lminmax = 10 * largest value in ORIENT
    lminmax = 10 * NINT(MAXVAL(DABS(ORIENT(:,:))))
    !
    !The oriented unit cell vectors are defined by the Miller indices
    DO i=1,3
      !Convert [hkil] notation into [uvw]
      u = 2.d0*ORIENT(i,1) + ORIENT(i,2)
      v = ORIENT(i,1) + 2.d0*ORIENT(i,2)
      w = ORIENT(i,3)
      !Check for common divisor
      IF( DABS(u)>0.1d0 .AND. NINT(DABS(v))>0.1d0 ) THEN
        z1 = GCD( NINT(DABS(u)) , NINT(DABS(v)) )
      ELSE
        z1 = MAX(DABS(u),DABS(v))
      ENDIF
      IF( DABS(u)>0.1d0 .AND. NINT(DABS(w))>0.1d0 ) THEN
        z2 = GCD( NINT(DABS(u)) , NINT(DABS(w)) )
      ELSE
        z2 = MAX(DABS(u),DABS(w))
      ENDIF
      IF( DABS(z1)>0.1d0 .AND. NINT(z2)>0.1d0 ) THEN
        x = GCD( NINT(DABS(z1)),NINT(DABS(z2)) )
      ELSE  !i.e. z1==0 or z2==0
        x = MAX( DABS(z1) , DABS(z2) )
      ENDIF
      IF( DABS(x)<0.1d0 ) x=1.d0  !avoid division by zero
      !Set box vector
      uv(i,:) = ( u*H(1,:) + v*H(2,:) + w*H(3,:) ) / x
    ENDDO
    !
    WRITE(msg,*) "Oriented unit cell vectors:"
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,'(3f16.6)') uv(1,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,'(3f16.6)') uv(2,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,'(3f16.6)') uv(3,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
    !For each atom in the unit cell H(:,:), keep only periodic replica that are inside the uv(:)
    !and store it in Q(:,:)
    WRITE(msg,*) "Duplicating atoms inside oriented unit cell..."
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !Estimate new number of particles NP by comparing volumes of old and new unit cells
    NP = 1.2d0*CEILING( SIZE(P,1) * DABS( DABS(uv(1,1)*uv(2,2)*uv(3,3)) / &
        & DABS(VECLENGTH(H(1,:))*VECLENGTH(H(2,:))*VECLENGTH(H(3,:))) ) ) + 20
    WRITE(msg,*) "Estimated new number of atoms : ", NP
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    IF(ALLOCATED(Q)) DEALLOCATE(Q)
    ALLOCATE( Q(NP,4) )
    Q(:,:) = 0.d0
    !
    !Loop over all replica in a wide range
    NP = 0
    DO i=1,SIZE(P,1)
      DO l=-lminmax,lminmax
        DO k=-lminmax,lminmax
          DO j=-lminmax,lminmax
            !Compute cartesian position of this replica
            tempP(1,1) = P(i,1) + DBLE(j)*H(1,1) + DBLE(k)*H(2,1) + DBLE(l)*H(3,1)
            tempP(1,2) = P(i,2) + DBLE(j)*H(1,2) + DBLE(k)*H(2,2) + DBLE(l)*H(3,2)
            tempP(1,3) = P(i,3) + DBLE(j)*H(1,3) + DBLE(k)*H(2,3) + DBLE(l)*H(3,3)
            tempP(1,4) = P(i,4)
            CALL CART2FRAC(tempP,uv)
            IF( tempP(1,1)>=-1.d-12 .AND. tempP(1,1)<1.d0-1.d-12 .AND.             &
              & tempP(1,2)>=-1.d-12 .AND. tempP(1,2)<1.d0-1.d-12 .AND.             &
              & tempP(1,3)>=-1.d-12 .AND. tempP(1,3)<1.d0-1.d-12       ) THEN
              !This replica is inside the new cell, mark it as new
              new = .TRUE.
              !Verify that its position is different from all previous atoms
              DO m=1,NP
                IF( DABS( VECLENGTH(tempP(1,1:3)-Q(m,1:3)) )<1.d-6 ) THEN
                  new = .FALSE.
                ENDIF
              ENDDO
              !
              IF( new ) THEN
                NP = NP+1
                IF(NP>SIZE(Q,1)) THEN
                  !Resize array Q
                  CALL RESIZE_DBLEARRAY2(Q,SIZE(Q,1)+10,SIZE(Q,2))
                ENDIF
                Q(NP,:) = tempP(1,:)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    !Replace old P with the new Q
    IF(ALLOCATED(P)) DEALLOCATE(P)
    ALLOCATE(P(NP,4))
    DO i=1,NP
      P(i,:) = Q(i,:)
    ENDDO
    IF(ALLOCATED(Q)) DEALLOCATE(Q)
    WRITE(msg,*) "new NP in oriented cell:", NP
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
    !Replace old H with the new (oriented) cell vectors
    !Align 1st box vector with Cartesian X axis
    !Place 2nd box vector in the XY plane
    H(:,:) = 0.d0
    H1 = VECLENGTH(uv(1,:))
    H2 = VECLENGTH(uv(2,:))
    H3 = VECLENGTH(uv(3,:))
    x = ANGVEC(uv(2,:),uv(3,:))
    y  = ANGVEC(uv(3,:),uv(1,:))
    z1 = ANGVEC(uv(1,:),uv(2,:))
    !Then convert this conventional notation into lower-triangular matrix H
    CALL CONVMAT(H1,H2,H3,x,y,z1,H)
    !Transform atom positions in Cartesian coordinates
    CALL FRAC2CART(P,H)
    !
  ELSE
    comment(1) = TRIM(comment(1))//" with box vectors H1=[2-1-10], H2=[-12-10], H3=[0001]"
    oriented = .FALSE. !set to .FALSE. to avoid duplicating atoms later
  ENDIF  !end if oriented
  !
  !
ELSE
  !The lattice is not cubic nor hexagonal
  IF( VECLENGTH(ORIENT(1,:)).NE.0.d0 .AND. VECLENGTH(ORIENT(2,:)).NE.0.d0 &
    & .AND. VECLENGTH(ORIENT(3,:)).NE.0.d0 ) THEN
    !The user asked for a crystal orientation, but the lattice is not cubic
    !=> write error message and exit
    CALL ATOMSK_MSG(4827,(/""/),(/0.d0/))
    nerr = nerr+1
    GOTO 1000
  ENDIF
  !
ENDIF !end if cubic or hexagonal
!
!
WRITE(msg,*) "final cell vectors:"
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,'(3f16.6)') H(1,:)
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,'(3f16.6)') H(2,:)
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,'(3f16.6)') H(3,:)
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
 comment(1) = TRIM(comment(1))//'.'
!
CALL ATOMSK_MSG(4028,(/TRIM(comment(1))/),(/0.d0/))
!
 comment(1) = '# '//TRIM(comment(1))
CALL ATOMSK_MSG(4029,(/''/),(/0.d0/))
!
!
!
400 CONTINUE
!Check that array P is allocated and number of atoms is not zero, otherwise display error message
IF( .NOT.ALLOCATED(P) .OR. SIZE(P,1)==0 ) THEN
  nerr = nerr+1
  CALL ATOMSK_MSG(3800,(/TRIM(msg)/),(/0.d0/))
  GOTO 1000
ENDIF
!Apply options to the created system
CALL OPTIONS_AFF(options_array,Huc,H,P,S,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
IF(nerr>0) GOTO 1000
!
!
!
500 CONTINUE
IF(wof) THEN
  !Write output file(s)
  msg = 'outputfile: '//TRIM(outputfile)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  CALL WRITE_AFF(outputfile,outfileformats,H,P,S,comment,AUXNAMES,AUX)
ENDIF
!
GOTO 1000
!
!
!
810 CONTINUE
nerr = nerr+1
GOTO 1000
!
820 CONTINUE
nerr = nerr+1
CALL ATOMSK_MSG(4803,(/''/),(/DBLE(SIZE(P,1)), DBLE(NT_NP) /))
GOTO 1000
!
!
!
1000 CONTINUE
!
!
!
END SUBROUTINE CREATE_CELL
!
!
END MODULE mode_create
