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
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 13 Jan. 2014                                     *
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
REAL(dp),DIMENSION(3,3),INTENT(IN):: ORIENT  !crystallographic orientation of the cell
!
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to output
CHARACTER(LEN=10):: create_struc  !structure to create (fcc, bcc...)
CHARACTER(LEN=32):: NT_type       !type of nanotube (zig-zag or armchair or chiral)
CHARACTER(LEN=128):: msg, temp
CHARACTER(LEN=4096):: outputfile
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
LOGICAL:: cubic     !is the system cubic?
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
REAL(dp),DIMENSION(3):: a_perp, a_para, b_perp, b_para, coord, Hrecip !for nanotubes
REAL(dp),DIMENSION(4):: tempP  !temporary position
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ips, uv     !interplanar spacing, unit vectors corresponding to new orientation ORIENT(:,:)
REAL(dp),DIMENSION(3,3):: ORIENTN     !normalized ORIENT
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, S  !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Q     !positions of atoms (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!Initialize variables
 cubic = .FALSE.
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
nspecies = 0
H(:,:) = 0.d0
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
END SELECT
!
!Determine if structure is cubic or not
SELECT CASE(create_struc)
CASE('sc','SC','fcc','FCC','L12','bcc','BCC','diamond','dia','zincblende','zb','ZB','perovskite','per','rocksalt','rs','RS')
  cubic = .TRUE.
CASE DEFAULT
  cubic = .FALSE.
END SELECT
WRITE(msg,*) "cubic:", cubic
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
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
CASE('sc')
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
  !Set up the messages
  WRITE(comment(1),*) 'Simple cubic '//TRIM(ADJUSTL((create_species(1))))
!
!
CASE('bcc')
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
  !Set up the messages
  WRITE(comment(1),*) TRIM(create_species(1))
  IF(nspecies==1) THEN
    comment(1) = 'Fcc '//TRIM(ADJUSTL(comment(1)))
  ELSE
    comment(1) = 'Fcc '//TRIM(ADJUSTL(comment(1)))//TRIM(ADJUSTL(create_species(2)))//' alloy'
  ENDIF
!
!
CASE('L12')
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
  !Set up the messages
  WRITE(comment(1),*) TRIM(create_species(1))//"3"//TRIM(create_species(2))
  comment(1) = 'L12 '//TRIM(ADJUSTL(comment(1)))
!
!
CASE('hcp')
  IF(nspecies.NE.1 .AND. nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0,2.d0 /))
    GOTO 810
  ENDIF
  ALLOCATE(P(2,4))
  !Set up the unit cell
  H(1,1) = create_a0(1)
  H(2,1) = create_a0(2)*DCOS(DEG2RAD(60.d0))
  H(2,2) = create_a0(2)*DSIN(DEG2RAD(60.d0))
  H(3,3) = create_a0(3)
  !Set up atom positions
  P(:,:) = 0.d0
  x = 1.d0/3.d0
  y = 1.d0/3.d0
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
!
CASE('dia','diamond','zincblende','zc')
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
  !Set up the messages
  WRITE(comment(1),*) TRIM(create_species(1))
  IF(nspecies==1) THEN
    comment(1) = TRIM(ADJUSTL(comment(1)))//' with diamond structure'
  ELSE
    comment(1) = TRIM(ADJUSTL(comment(1)))//TRIM(ADJUSTL(create_species(2)))//' with zincblende structure'
  ENDIF
!
!
CASE('rocksalt','rs')
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
  !Set up the messages
  WRITE(comment(1),*) TRIM(create_species(1))
  comment(1) = 'Rocksalt '//TRIM(ADJUSTL(comment(1)))//TRIM(create_species(2))
!
!
CASE('per','perovskite')
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
  !Set up the messages
  temp = TRIM(create_species(1))//TRIM(create_species(2))//TRIM(create_species(3))//'3'
  WRITE(comment(1),*) 'Cubic perovskite '//TRIM(temp)
!
!
CASE('graphite')
  IF(nspecies.NE.1 .AND. nspecies.NE.2) THEN
    CALL ATOMSK_MSG(4804,(/''/),(/ 1.d0,2.d0 /))
    GOTO 810
  ENDIF
  ALLOCATE(P(4,4))
  !Set up atom positions
  P(:,:) = 0.d0
   P(2,1) = 1.d0
   P(2,2) = 1.d0/DSQRT(3.d0)
  P(3,3) = 0.5d0
   P(4,1) = 0.5d0
   P(4,2) = P(2,2)/2.d0
   P(4,3) = 0.5d0
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
  P(3,4) = P(1,4)
  P(4,4) = P(2,4)
  !Set up the unit cell
  H(1,1) = create_a0(1)
   H(2,1) = create_a0(1)/2.d0
   H(2,2) = DSQRT(3.d0)*create_a0(2)/2.d0
  H(3,3) = create_a0(3)
  !Set up the messages
  temp = TRIM(create_species(1))//TRIM(create_species(2))
  WRITE(comment(1),*) TRIM(temp)//' with hexagonal graphite structure'
!
!
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
! -- add other structures here --
!
CASE DEFAULT
  nerr = nerr+1
  CALL ATOMSK_MSG(4805,(/create_struc/),(/0.d0/))
  GOTO 1000
!
END SELECT
!
!
!
300 CONTINUE
!Orient the system so that the cartesian axes are X=ORIENT(1,:); Y=ORIENT(2,:); Z=ORIENT(3,:)
IF( cubic .AND. VECLENGTH(ORIENT(1,:)).NE.0.d0 .AND. VECLENGTH(ORIENT(2,:)).NE.0.d0 .AND. VECLENGTH(ORIENT(3,:)).NE.0.d0 ) THEN
  comment(1) = TRIM(comment(1))//' oriented'
  DO i=1,3
    IF(i==1) comment(1) = TRIM(comment(1))//' X=['
    IF(i==2) comment(1) = TRIM(comment(1))//' Y=['
    IF(i==3) comment(1) = TRIM(comment(1))//' Z=['
    DO j=1,3
      WRITE(msg,*) NINT(ORIENT(i,j))
      comment(1) = TRIM(comment(1))//TRIM(ADJUSTL(msg))
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
  IF( ANGVEC(ORIENT(1,:),ORIENT(2,:))-90.d0<1.d-6 ) orthovec(1)=.TRUE.
  IF( ANGVEC(ORIENT(2,:),ORIENT(3,:))-90.d0<1.d-6 ) orthovec(2)=.TRUE.
  IF( ANGVEC(ORIENT(3,:),ORIENT(1,:))-90.d0<1.d-6 ) orthovec(3)=.TRUE.
  WRITE(msg,*) "orthovec: ", orthovec(:)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  IF( ANY(.NOT.orthovec) ) THEN
    CALL ATOMSK_MSG(4819,(/''/),(/0.d0/))
    nerr = nerr+1
    GOTO 1000
  ENDIF
  !
  !Compute the interplanar spacings that correspond to the new orientation
  !For a cubic system the interplanar spacing along a Miller axis hkl is defined by:
  !   d_hkl = a² / sqrt( h² + k² + l² )
  ips(:,:) = 0.d0
  DO i=1,3
    ips(i,i) = create_a0(1) / DSQRT( ORIENT(i,1)**2 + ORIENT(i,2)**2 + ORIENT(i,3)**2 )
  ENDDO
  WRITE(msg,*) "interplanar spacings along X, Y, Z:"
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,'(3f16.6)') VECLENGTH(ips(1,:)), VECLENGTH(ips(2,:)), VECLENGTH(ips(3,:))
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !
  !Unit cell is oriented so that H(1,:)//[100], H(2,:)//[010], H(3,:)//[001]
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
  !The unit cell is *not* defined by the interplanar spacing, but by a repetition vector
  !=> find the shortest repetition vectors
  !   they must be a multiple of the ips(:) and a linear combination of the H(:,:)
  uv(:,:) = HUGE(1.d0)
  DO i=1,3
    !multiply ips(i,i) by m and check it against a linear combination of H(:,:)
    DO m=1,10
      DO l=-lminmax,lminmax
        DO k=-lminmax,lminmax
          DO j=-lminmax,lminmax
            !Compute the difference between current m*ips(i,i) and current linear combination of H(:,:)
            x = VECLENGTH( DBLE(m)*ips(i,:) + DBLE(j)*H(1,:) + DBLE(k)*H(2,:) + DBLE(l)*H(3,:) )
            IF( x < 1.d-6 .AND. DBLE(m)*ips(i,i)<uv(i,i) ) THEN
              !We have found a suitable vector, and it is shorter than the previous one
              !=> save it to uv(i,i)
              uv(i,i) = DBLE(m)*ips(i,i)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    320 CONTINUE
  ENDDO
  !
  !Correct for some special orientations to find the minimal repetition unit
  SELECT CASE(create_struc)
  CASE('bcc','BCC')
    !in bcc structure, period along a <111> direction is actually 1/2<111>
    DO i=1,3
      IF( DABS(ORIENT(i,1))==DABS(ORIENT(i,2)) .AND. DABS(ORIENT(i,1))==DABS(ORIENT(i,3)) ) THEN
        uv(i,i) = uv(i,i)/2.d0
      ENDIF
    ENDDO
    !
  CASE('fcc','FCC','diamond','dia','zincblende','zb','rocksalt','rs')
    !in fcc, diamond/zb, and rocksalt structures,
    !period along a <110> direction is actually 1/2<110>
    !period along a <112> direction is actually 1/2<112>
    DO i=1,3
      IF( DABS(ORIENT(i,1))==DABS(ORIENT(i,2)) .AND. DABS(ORIENT(i,3))==0.d0 .OR.           &
        & DABS(ORIENT(i,1))==DABS(ORIENT(i,3)) .AND. DABS(ORIENT(i,2))==0.d0 .OR.           &
        & DABS(ORIENT(i,2))==DABS(ORIENT(i,3)) .AND. DABS(ORIENT(i,1))==0.d0       ) THEN
        uv(i,i) = uv(i,i)/2.d0
      ELSEIF( DABS(ORIENT(i,1))==DABS(ORIENT(i,2)) .AND. DABS(ORIENT(i,3))==2.d0*DABS(ORIENT(i,1)) .OR.           &
            & DABS(ORIENT(i,1))==DABS(ORIENT(i,3)) .AND. DABS(ORIENT(i,2))==2.d0*DABS(ORIENT(i,1)) .OR.           &
            & DABS(ORIENT(i,2))==DABS(ORIENT(i,3)) .AND. DABS(ORIENT(i,1))==2.d0*DABS(ORIENT(i,2))       ) THEN
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
  IF(ALLOCATED(Q)) DEALLOCATE(Q)
  ALLOCATE( Q(1000,4) )    !assuming there are less than 1,000 atoms in the unit cell...
  Q(:,:) = 0.d0
  NP=0
  !
  !Loop over all replica in a wide range
  DO i=1,SIZE(P,1)
    DO l=-lminmax,lminmax
      DO k=-lminmax,lminmax
        DO j=-lminmax,lminmax
          !Compute cartesian position of this replica
          tempP(1) = P(i,1) + DBLE(j)*H(1,1) + DBLE(k)*H(2,1) + DBLE(l)*H(3,1)
          tempP(2) = P(i,2) + DBLE(j)*H(1,2) + DBLE(k)*H(2,2) + DBLE(l)*H(3,2)
          tempP(3) = P(i,3) + DBLE(j)*H(1,3) + DBLE(k)*H(2,3) + DBLE(l)*H(3,3)
          tempP(4) = P(i,4)
          IF( tempP(1)>=0.d0 .AND. tempP(1)<uv(1,1)-1.d-16 .AND.             &
            & tempP(2)>=0.d0 .AND. tempP(2)<uv(2,2)-1.d-16 .AND.             &
            & tempP(3)>=0.d0 .AND. tempP(3)<uv(3,3)-1.d-16       ) THEN
            !This replica is inside the new cell => keep it
            NP = NP+1
            IF(NP>SIZE(Q,1)) THEN
              nerr = nerr+1
              GOTO 1000
            ENDIF
            Q(NP,:) = tempP(:)
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
  WRITE(msg,*) "new cell vectors:"
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,'(3f16.6)') H(1,:)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,'(3f16.6)') H(2,:)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,'(3f16.6)') H(3,:)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !
ELSEIF( cubic ) THEN
  comment(1) = TRIM(comment(1))//" oriented X=[100], Y=[010], Z=[001]"
ENDIF
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
!Apply options to the created system
CALL OPTIONS_AFF(options_array,H,P,S,AUXNAMES,AUX,ORIENT,SELECT)
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
GOTO 900
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
900 CONTINUE
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
