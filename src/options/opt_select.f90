MODULE select
!
!**********************************************************************************
!*  SELECT                                                                        *
!**********************************************************************************
!* This module selects a region of the system, by returning a mask array          *
!* called SELECT: selected atoms have a .TRUE. value, others have                 *
!* a .FALSE. value.                                                               *
!**********************************************************************************
!* (C) November 2010 - Pierre Hirel                                               *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 28 June 2018                                     *
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
USE comv
USE constants
USE messages
USE neighbors
USE files
USE functions
USE in_stl
USE subroutines
USE sorting
!
!
CONTAINS
!
!
SUBROUTINE SELECT_XYZ(H,P,AUXNAMES,AUX,select_multiple,region_side,region_geom,region_dir,region_1,region_2,ORIENT,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=2):: species
CHARACTER(LEN=3):: rand_sp      !species of atoms to select randomly
CHARACTER(LEN=16):: select_multiple  !empty or 'add' or 'rm' or 'xor' or 'among'
CHARACTER(LEN=16):: region_dir  !x, y, z, or crystallographic direction
CHARACTER(LEN=128):: msg
CHARACTER(LEN=4096):: temp
CHARACTER(LEN=128):: region_side  !'in' or 'out' or 'all' or 'inv' or 'neigh' or 'list' or 'grid' or 'stl'
CHARACTER(LEN=4096):: region_geom !geometry of the region: "box" or "sphere". If "neighbors" then
                                  !neighbors of an atom must be searched. If region_side=="list" then
                                  !region_geom contains the name of the file.
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
LOGICAL:: add, rm, intersect, xor, among !add, remove, intersect, XOR, among from previous selection?
LOGICAL:: selectmul  !modify (add/rm/intersect/xor) previous selection with current criteria?
LOGICAL:: exceeds100 !are there more than 100 neighboring atoms?
LOGICAL:: isreduced  !are positions in reduced coordinates?
LOGICAL:: keep  !keep atom?
LOGICAL:: toselect !atomindices(:) contains atoms to select? (if FALSE, atoms to un-select)
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: SELECT !the final selection array
LOGICAL,DIMENSION(:),ALLOCATABLE:: prevSELECT  !copy of previous selection
INTEGER:: a1, a2, a3
INTEGER:: atomrank, atomrank2  !rank of atom to be selected
INTEGER:: gridformat !format of the file (for select grid)
INTEGER:: i, j, k, l, m, n
INTEGER:: line
INTEGER:: Ngrid    !number of elements in the grid
INTEGER:: Nperline !number of grid element per line
INTEGER:: Nselect  !number of atoms selected
INTEGER:: Nadded, Nrm  !number of atoms added or removed from selection
INTEGER:: rand_N   !number of atoms to select randomly
INTEGER:: sp_N     !number of atoms of the given species that exist in P
INTEGER,DIMENSION(:),ALLOCATABLE:: atomindices !indices of atom(s) that must be selected
INTEGER,DIMENSION(:),ALLOCATABLE:: newindex    !list of sorted indexes
INTEGER,DIMENSION(:),ALLOCATABLE:: Nlist !index of neighbours
INTEGER,DIMENSION(:,:),ALLOCATABLE:: NeighList  !the neighbor list
REAL(dp):: distance  !distance of an atom to the center of the sphere
REAL(dp):: dz        !increment along Z
REAL(dp):: ta, tf, tt, tu, tv !used to detect ray-triangle intersections
REAL(dp):: snumber   !atomic number
REAL(dp):: tempreal
REAL(dp):: V1, V2, V3  !vector components
REAL(dp):: xmin, xmax, ymin, ymax, zmin, zmax !box parameters
REAL(dp):: txmin, txmax, tymin, tymax, tzmax  !min/max X,Y,Z coordinates of a triangle
REAL(dp),DIMENSION(3):: e1, e2, td, th, ts, tq !vectors used to detect ray-triangle intersections
REAL(dp),DIMENSION(3):: region_1 !First corner for'box', or center of sphere
REAL(dp),DIMENSION(3):: region_2 !Last corner for'box', or radius of sphere
REAL(dp),DIMENSION(3):: T1, T2, T3 !positions of vertices of a triangle
REAL(dp),DIMENSION(1,3):: Vplane  !crystallographic vector defining the plane
REAL(dp),DIMENSION(3,3),INTENT(IN):: H      !supercell vectors
REAL(dp),DIMENSION(3,3),INTENT(IN):: ORIENT !current crystallographic orientation of the system
REAL(dp),DIMENSION(3,3):: ORIENTN      !normalized ORIENT
REAL(dp),DIMENSION(:),ALLOCATABLE:: randarray    !random numbers
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries   !species, Natoms of this species
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX        !auxiliary properties of atoms/shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: GRID       !positions of elements in a finite-element grid
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: facenormals !vectors normal to faces (prism)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: triangles  !triangles defining a 3-D surface
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Q     !positions of atoms of a given species
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: PosList !positions of neighbors
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: V_NN    ! positions of 1st nearest neighbors
!
!
!Initialize variables
add = .FALSE.
rm = .FALSE.
intersect = .FALSE.
xor = .FALSE.
among = .FALSE.
selectmul = .FALSE.
a1 = 1
a2 = 1
a3 = 1
atomrank=-1
atomrank2=-1
i = 0
Nselect = 0
Nadded = 0
Nrm = 0
snumber = 0.d0
IF(ALLOCATED(aentries)) DEALLOCATE(aentries)
IF(ALLOCATED(atomindices)) DEALLOCATE(atomindices)
IF(ALLOCATED(facenormals)) DEALLOCATE(facenormals)
IF(ALLOCATED(GRID)) DEALLOCATE(GRID)
IF(ALLOCATED(prevSELECT)) DEALLOCATE(prevSELECT)
!
!
msg = 'Entering SELECT_XYZ'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
SELECT CASE(region_side)
!
CASE('all','any','none','above','below','invert','inv','in','list','out','prop','property','random','rand','random%','rand%','stl')
  CONTINUE
  !
CASE DEFAULT
  !Check if region_side contains atom index (indices)
  j=SCAN(region_side,":")
  k=SCAN(region_side,",")
  IF( j>0 .OR. k>0 ) THEN
    !region_side contains several indices, or a range of atom indices
    IF( k>0 ) THEN
      region_geom = "list "
    ELSE
      !j>0, i.e. there are only two integers separated by a colon
      region_geom = "range"
    ENDIF
    !Replace all commas by blanck spaces
    DO WHILE(k>0)
      region_side(k:k) = " "
      k=SCAN(region_side,",")
    ENDDO
    !
    IF(verbosity==4) THEN
      msg = "Detected a "//TRIM(ADJUSTL(region_geom))//" of index"
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      msg = 'List: '//TRIM(region_side)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDIF
    !Determine how many atoms must be selected
    Nselect=0
    msg = ADJUSTL(region_side)
    DO WHILE( LEN_TRIM(msg)>0 )
      READ(msg,*) temp
      j=SCAN(temp,":")
      IF( j>0 ) THEN
        !Read range of numbers:
        !First number
        IF( LEN_TRIM(temp(1:j-1)) <=0 ) THEN
          !Nothing before the colon: consider it is atom #1
          atomrank = 1
        ELSE
          READ(temp(1:j-1),*,ERR=800,END=800) atomrank
        ENDIF
        region_1(1) = atomrank
        !Second number
        IF( LEN_TRIM(temp(j+1:)) <=0 ) THEN
          !Nothing before the colon: consider it is last atom
          atomrank2 = SIZE(P,1)
        ELSE
          READ(temp(j+1:),*,ERR=800,END=800) atomrank2
        ENDIF
        !Check that atomrank2 > atomrank
        IF( atomrank2<atomrank ) THEN
          k = atomrank
          atomrank = atomrank2
          atomrank2 = k
        ENDIF
        region_1(2) = atomrank2
        !All atoms between atomrank and atomrank2 will have to be selected
        Nselect = Nselect + (atomrank2-atomrank+1)
      ELSE
        !Read number
        READ(temp,*) atomrank
        Nselect = Nselect+1
      ENDIF
      j=SCAN(msg," ")
      msg = ADJUSTL(msg(j:))
    ENDDO
    !
    !Indices of atoms that must be selected will be saved in atomindices
    ALLOCATE(atomindices(Nselect))
    atomindices(:) = 0
    !
    !Store each number and range into atomindices(:)
    Nselect=0
    msg = ADJUSTL(region_side)
    DO WHILE( LEN_TRIM(msg)>0 )
      READ(msg,*) temp
      j=SCAN(temp,":")
      IF( j>0 ) THEN
        !Read range of numbers:
        !First number
        IF( LEN_TRIM(temp(1:j-1)) <=0 ) THEN
          !Nothing before the colon: consider it is atom #1
          atomrank = 1
        ELSE
          READ(temp(1:j-1),*,ERR=800,END=800) atomrank
        ENDIF
        !Second number
        IF( LEN_TRIM(temp(j+1:)) <=0 ) THEN
          !Nothing before the colon: consider it is last atom
          atomrank2 = SIZE(P,1)
        ELSE
          READ(temp(j+1:),*,ERR=800,END=800) atomrank2
        ENDIF
        !Check that atomrank2 > atomrank
        IF( atomrank2<atomrank ) THEN
          k = atomrank
          atomrank = atomrank2
          atomrank2 = k
        ENDIF
        DO i=atomrank,atomrank2
          Nselect = Nselect+1
          atomindices(Nselect) = i
        ENDDO
      ELSE
        !Read number
        READ(temp,*) atomrank
        Nselect = Nselect+1
        atomindices(Nselect) = atomrank
      ENDIF
      j=SCAN(msg," ")
      msg = ADJUSTL(msg(j:))
    ENDDO
    !
    !Set region_side to "index" for later
    region_side = "index"
    !
  ELSE
    !try to read an integer
    READ(region_side,*,ERR=10,END=10) atomrank
    !It does contain an integer => save it to atomindices(1)
    ALLOCATE(atomindices(1))
    atomindices(:) = atomrank
    !modify region_side
    region_side = "index"
    !Save index into region_1(:) for the message
    region_1(1) = DBLE(atomrank)
  ENDIF
  !
END SELECT
!
!
10 CONTINUE
IF( verbosity==4 ) THEN
  WRITE(msg,*) "region_side = ", TRIM(region_side)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "region_geom = ", TRIM(ADJUSTL(region_geom))
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "region_dir = ", TRIM(region_dir)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "region_1 = ", region_1(:)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "region_2 = ", region_2(:)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!
!
!A message to the user of what this option will do
CALL ATOMSK_MSG(2077, (/ region_side//'    ',                 &
     &                 region_geom, region_dir//'       ',    &
     &                 select_multiple//'    '            /), &
     & (/ region_1(1), region_1(2), region_1(3),              &
     &    region_2(1), region_2(2), region_2(3) /))
!
selectmul = .TRUE.
SELECT CASE(select_multiple)
CASE("add")
  add = .TRUE.
CASE("rm")
  rm = .TRUE.
CASE("intersect")
  intersect = .TRUE.
CASE("xor")
  xor = .TRUE.
CASE("among")
  among = .TRUE.
  selectmul = .FALSE.
CASE DEFAULT
  selectmul = .FALSE.
END SELECT
!
IF( .NOT. ALLOCATED(SELECT) ) THEN
  IF( add .OR. rm .OR. intersect .OR. xor .OR. among ) THEN
    !User asked to add or remove atoms from previous selection,
    !but no previous selection exists => display warning and
    !proceed with a "normal" selection
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2761,(/select_multiple/),(/0.d0/))
    selectmul = .FALSE.
    add = .FALSE.
    rm = .FALSE.
    intersect = .FALSE.
    xor = .FALSE.
  ENDIF
ENDIF
!
IF( (selectmul.OR.among) .AND. ALLOCATED(SELECT) ) THEN
  !Copy previous selection to a new array
  ALLOCATE(prevSELECT(SIZE(SELECT)))
  prevSELECT(:) = SELECT(:)
  DEALLOCATE(SELECT)
  !Prepare new array
  ALLOCATE(SELECT(SIZE(P,1)))
  SELECT(:) = .FALSE.
ELSE
  !This option will work with array SELECT directly
ENDIF
!
IF( verbosity==4 ) THEN
  IF( add ) THEN
    WRITE(msg,*) "new selection will be ADDED to previous one"
  ELSEIF( rm ) THEN
    WRITE(msg,*) "new selection will be REMOVED from previous one"
  ELSEIF( intersect ) THEN
    WRITE(msg,*) "new selection will be INTERSECTED with previous one"
  ELSEIF( xor ) THEN
    WRITE(msg,*) "new selection will be XOR with previous one"
  ELSEIF( xor ) THEN
    WRITE(msg,*) "new selection will be chosen AMONG previous one"
  ENDIF
ENDIF
!
!
!
100 CONTINUE
!The content of region_side decides what must be done
SELECT CASE(region_side)
!
CASE('all','any','none')
  !The selection must be cleared (this is equivalent to selecting all atoms)
  IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
  Nselect = SIZE(P,1)
  selectmul = .FALSE.
  !
  !
CASE('invert','inv')
  !Previous selection must be inverted
  selectmul = .FALSE.
  IF( .NOT.ALLOCATED(SELECT) .OR. SIZE(SELECT)<=0 ) THEN
    !Previously all atoms were selected => inverting would mean selecting no atoms
    !However selecting no atoms doesn't make sense, so we do nothing (SELECT remains unallocated)
    Nselect = SIZE(P,1)
  ELSE
    !Atoms were selected before => invert selection
    DO i=1,SIZE(SELECT)
      SELECT(i) = .NOT.SELECT(i)
      IF(SELECT(i)) THEN
        Nselect = Nselect+1
      ENDIF
    ENDDO
  ENDIF
  !
  !
CASE('index')
  !All atoms must be un-selected, but one
  IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
  Nselect=0
  IF( ALLOCATED(atomindices) ) THEN
    ALLOCATE( SELECT(SIZE(P,1)) )
    SELECT(:) = .FALSE.
    DO i=1,SIZE(atomindices)
      IF( atomindices(i)>0 .AND. atomindices(i)<=SIZE(P,1) ) THEN
        !WARNING: an index can appear twice in the list of atomindices(:)
        ! do not double-count selected atoms
        IF( .NOT.SELECT(atomindices(i)) ) THEN
          SELECT(atomindices(i)) = .TRUE.
          Nselect = Nselect+1
        ELSE
          nwarn = nwarn+1
          CALL ATOMSK_MSG(2753,(/""/),(/DBLE(atomindices(i))/))
        ENDIF
      ELSE
        !Atom index is out-of-bounds
        nwarn = nwarn+1
        CALL ATOMSK_MSG(2742,(/""/),(/DBLE(atomindices(i))/))
      ENDIF
    ENDDO
  ENDIF
  !
  !
CASE('above','below')
  !Select atoms that are above or below a distance along an axis
  !region_dir=axis normal to plane; region_1(1)=distance along that axis
  IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
  ALLOCATE( SELECT( SIZE(P,1) ) )
  SELECT(:) = .TRUE.
  !
  SELECT CASE(region_dir)
  CASE("x","X","y","Y","z","Z")
    !region_dir contains a cartesian direction
    !Define the axes: a3 is the direction of the main axis
    SELECT CASE(region_dir)
    CASE("x","X")
      a1 = 1
    CASE("y","Y")
      a1 = 2
    CASE("z","Z")
      a1 = 3
    CASE DEFAULT
      CALL ATOMSK_MSG(2800,(/region_dir/),(/0.d0/))
      nerr = nerr+1
      GOTO 1000
    END SELECT
    WRITE(msg,*) 'a1 = ', a1
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    WRITE(msg,*) 'region_1(1) = ', region_1(1)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
    IF(region_side=="above") THEN
      !Select only atoms that are above the plane
      DO i=1,SIZE(P,1)
        IF( P(i,a1)>region_1(1) ) THEN
          !atom is above the plane
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ELSE
          !atom is below the plane
          SELECT(i) = .FALSE.
        ENDIF
      ENDDO
    ELSE
      !Select only atoms that are below the plane
      DO i=1,SIZE(P,1)
        IF( P(i,a1)<region_1(1) ) THEN
          !atom is below the plane
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ELSE
          !atom is above the plane
          SELECT(i) = .FALSE.
        ENDIF
      ENDDO
    ENDIF
    !
  CASE DEFAULT
    !region_dir should contain a crystallograhic direction
    !convert it to a vector and save it in Vplane(1,:)
    CALL INDEX_MILLER(region_dir,Vplane(1,:),j)
    IF(j>0) GOTO 800
    !
    !If the system has a defined crystallographic orientation ORIENT,
    !then Vplane(1,:) is defined in that basis
    !=> rotate Vplane(1,:) to express it in cartesian basis
    IF( ANY( NINT(ORIENT(:,:)).NE.0 ) ) THEN
      DO i=1,3
        ORIENTN(i,:) = ORIENT(i,:) / VECLENGTH(ORIENT(i,:))
      ENDDO
      V1 = Vplane(1,1)
      V2 = Vplane(1,2)
      V3 = Vplane(1,3)
      Vplane(1,1) = ORIENTN(1,1)*V1 + ORIENTN(1,2)*V2 + ORIENTN(1,3)*V3
      Vplane(1,2) = ORIENTN(2,1)*V1 + ORIENTN(2,2)*V2 + ORIENTN(2,3)*V3
      Vplane(1,3) = ORIENTN(3,1)*V1 + ORIENTN(3,2)*V2 + ORIENTN(3,3)*V3
    ENDIF
    !Normalize Vplane
    Vplane(1,:) = Vplane(1,:)/VECLENGTH(Vplane(1,:))
    WRITE(msg,'(a8,3f12.3)') 'Vplane: ', Vplane(1,:)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
    !For each atom, determine if it is above or below the plane and select it
    IF(region_side=="above") THEN
      !Select only atoms that are above the plane
      DO i=1,SIZE(P,1)
        !determine if atom position is above or below the plane
        tempreal = VEC_PLANE( Vplane(1,:) , region_1(1) , P(i,1:3) )
        IF(tempreal>0.d0) THEN
          !atom is above the plane
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ELSE
          !atom is below or in the plane
          SELECT(i) = .FALSE.
        ENDIF
      ENDDO
      !
    ELSE
      !Select only atoms that are below the plane
      DO i=1,SIZE(P,1)
        !determine if atom position is above or below the plane
        tempreal = VEC_PLANE( Vplane(1,:) , region_1(1) , P(i,1:3) )
        IF(tempreal<0.d0) THEN
          !atom is below the plane
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ELSE
          !atom is above or in the plane
          SELECT(i) = .FALSE.
        ENDIF
      ENDDO
      !
    ENDIF
    !
  END SELECT
  !
  !
  !
CASE('in','out')
  !Select atoms inside or outside of a region
  IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
  ALLOCATE( SELECT( SIZE(P,1) ) )
  SELECT(:) = .TRUE.
  !
  SELECT CASE(region_geom)
  !
  CASE('cell')
    !Select all atoms that are inside or outside of the simulation cell
    !Convert atom positions into reduced/fractional coordinates
    !Check if coordinates are already reduced or not
    CALL FIND_IF_REDUCED(P,isreduced)
    !If not reduced, then reduce them
    IF( .NOT.isreduced ) THEN
      CALL CART2FRAC(P,H)
    ENDIF
    !
    !Loop on all atoms
    DO i=1,SIZE(P,1)
      IF( P(i,1)>=0.d0 .AND. P(i,1)<1.d0 .AND. &
        & P(i,2)>=0.d0 .AND. P(i,2)<1.d0 .AND. &
        & P(i,3)>=0.d0 .AND. P(i,3)<1.d0     ) THEN
        !If atom is inside the cell...
        IF(region_side=='in') THEN
          !...and if we want to select the inside of the cell, set to true
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ELSE
          !...and if we want to select the outside of the cell, set to false
          SELECT(i) = .FALSE.
        ENDIF
      !
      ELSE
        !If atom is NOT inside the cell...
        IF(region_side=='in') THEN
          !...and if we want to select the inside of the cell, set to false
          SELECT(i) = .FALSE.
        ELSE
          !...and if we want to select the outside of the cell, set to true
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ENDIF
      ENDIF
    ENDDO
    !
    !Convert atom positions back to Cartesian coordinates
    IF( .NOT.isreduced ) THEN
      CALL FRAC2CART(P,H)
    ENDIF
    !
  CASE('box')
    !Select atoms that are inside/outside a rectangular box
    !region_1(:) = position of first corner of the box
    !region_2(:) = position of last corner of the box
    xmin = MIN( region_1(1),region_2(1) )
    xmax = MAX( region_1(1),region_2(1) )
    ymin = MIN( region_1(2),region_2(2) )
    ymax = MAX( region_1(2),region_2(2) )
    zmin = MIN( region_1(3),region_2(3) )
    zmax = MAX( region_1(3),region_2(3) )
    !
    DO i=1,SIZE(P(:,1))
      IF( P(i,1)>xmin .AND. P(i,1)<xmax .AND.                &
        & P(i,2)>ymin .AND. P(i,2)<ymax .AND.                &
        & P(i,3)>zmin .AND. P(i,3)<zmax       ) THEN
        !If atom is inside the box...
        IF(region_side=='in') THEN
          !...and if we want to select the inside of the box, set to true
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ELSE
          !...and if we want to select the outside of the box, set to false
          SELECT(i) = .FALSE.
        ENDIF
      !
      ELSE
        !If atom is NOT inside the box...
        IF(region_side=='in') THEN
          !...and if we want to select the inside of the box, set to false
          SELECT(i) = .FALSE.
        ELSE
          !...and if we want to select the outside of the box, set to true
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ENDIF
      ENDIF
    ENDDO
    !
    !
  CASE('sphere')
    !region_1(:) = center of the sphere
    !region_2(1) = radius of the sphere
    DO i=1,SIZE(P(:,1))
      distance = VECLENGTH( P(i,1:3)-region_1(1:3) )
      IF( distance<region_2(1) ) THEN
        !If atom is inside the sphere...
        IF(region_side=='in') THEN
          !...and if we want to select the inside of the sphere, set to true
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ELSE
          !...and if we want to select the outside of the sphere, set to false
          SELECT(i) = .FALSE.
        ENDIF
      ELSE
        !If atom is NOT inside the sphere...
        IF(region_side=='in') THEN
          !...and if we want to select the inside of the sphere, set to false
          SELECT(i) = .FALSE.
        ELSE
          !...and if we want to select the outside of the sphere, set to true
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ENDIF
      ENDIF
    ENDDO
    !
    !
  CASE('cylinder')
    !Define the axes: a3 is the direction of the main axis
    SELECT CASE(region_dir)
    CASE("x","X")
      a1 = 2
      a2 = 3
      a3 = 1
    CASE("y","Y")
      a1 = 3
      a2 = 1
      a3 = 2
    CASE("z","Z")
      a1 = 1
      a2 = 2
      a3 = 3
    CASE DEFAULT
      CALL ATOMSK_MSG(2800,(/region_dir/),(/0.d0/))
      nerr = nerr+1
      GOTO 1000
    END SELECT
    !region_1(1:2) = center of the cylinder (region_1(3) is ignored here)
    !region_2(1) = radius of the cylinder
    DO i=1,SIZE(P(:,1))
      distance = VECLENGTH( (/P(i,a1)-region_1(1), P(i,a2)-region_1(2),0.d0/) )
      IF( distance<region_2(1) ) THEN
        !If atom is inside the cylinder...
        IF(region_side=='in') THEN
          !...and if we want to select the inside of the cylinder, set to true
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ELSE
          !...and if we want to select the outside of the cylinder, set to false
          SELECT(i) = .FALSE.
        ENDIF
      ELSE
        !If atom is NOT inside the sphere...
        IF(region_side=='in') THEN
          !...and if we want to select the inside of the cylinder, set to false
          SELECT(i) = .FALSE.
        ELSE
          !...and if we want to select the outside of the cylinder, set to true
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ENDIF
      ENDIF
    ENDDO
    !
    !
  CASE('cone')
    !Define the axes: a3 is the direction normal to the base of the cone
    SELECT CASE(region_dir)
    CASE("x","X")
      a1 = 2
      a2 = 3
      a3 = 1
    CASE("y","Y")
      a1 = 3
      a2 = 1
      a3 = 2
    CASE("z","Z")
      a1 = 1
      a2 = 2
      a3 = 3
    CASE DEFAULT
      CALL ATOMSK_MSG(2800,(/region_dir/),(/0.d0/))
      nerr = nerr+1
      GOTO 1000
    END SELECT
    !region_1(1:3) = position of the tip of the cone
    !region_2(1) = angle between direction a3 and surface of the cone (degrees)
    !Convert region_2(1) into radians
    DO WHILE( region_2(1) > 90.d0 )
      region_2(1) = region_2(1) - 180.d0
    ENDDO
    DO WHILE( region_2(1) < -90.d0 )
      region_2(1) = region_2(1) + 180.d0
    ENDDO
    region_2(1) = DEG2RAD(region_2(1))
    !
    !Select atoms that are in/out of the cone
    DO i=1,SIZE(P,1)
      !Compute position of the torus in the direction of the atom
      distance = ( P(i,a1)-region_1(a1) )**2 + ( P(i,a2)-region_1(a2) )**2 &
               & - ( (P(i,a3)-region_1(a3))**2 )*((DTAN(region_2(1)))**2)
      IF( P(i,a3)>region_1(a3) .AND. distance < 0.d0 ) THEN
        !The atom is inside the cone
        IF(region_side=='in') THEN
          !...and we want to select the inside of the cone => set to true
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ELSE
          !...and we want to select the outside of the cone => set to false
          SELECT(i) = .FALSE.
        ENDIF
      ELSE
        !The atom is outside the cone
        IF(region_side=='in') THEN
          !...and we want to select the inside of the cone => set to false
          SELECT(i) = .FALSE.
        ELSE
          !...and we want to select the outside of the cone => set to true
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ENDIF
      ENDIF
    ENDDO
    !
    !
  CASE('torus')
    !Define the axes: a3 is the direction normal to the torus plane
    SELECT CASE(region_dir)
    CASE("x","X")
      a1 = 2
      a2 = 3
      a3 = 1
    CASE("y","Y")
      a1 = 3
      a2 = 1
      a3 = 2
    CASE("z","Z")
      a1 = 1
      a2 = 2
      a3 = 3
    CASE DEFAULT
      CALL ATOMSK_MSG(2800,(/region_dir/),(/0.d0/))
      nerr = nerr+1
      GOTO 1000
    END SELECT
    !region_1(1:3) = position of the center of the torus
    !region_2(1) = main radius of torus
    !region_2(2) = secondary radius of torus
    !
    !Select atoms that are in/out of the torus
    DO i=1,SIZE(P,1)
      !Compute position of the torus in the direction of the atom
      Vplane(1,a1) = P(i,a1) - region_1(a1)
      Vplane(1,a2) = P(i,a2) - region_1(a2)
      Vplane(1,a3) = 0.d0
      IF( VECLENGTH(Vplane(1,:)) > 1.d-6 ) THEN
        Vplane(1,:) = region_1(1:3) + Vplane(1,:) * region_2(1) / VECLENGTH(Vplane(1,:))
      ENDIF
      Vplane(1,a3) = region_1(a3)
      distance = VECLENGTH( Vplane(1,:) - P(i,1:3) )
      IF( distance < region_2(2) ) THEN
        !The atom is inside the torus
        IF(region_side=='in') THEN
          !...and we want to select the inside of the torus => set to true
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ELSE
          !...and we want to select the outside of the torus => set to false
          SELECT(i) = .FALSE.
        ENDIF
      ELSE
        !The atom is outside the torus
        IF(region_side=='in') THEN
          !...and we want to select the inside of the torus => set to false
          SELECT(i) = .FALSE.
        ELSE
          !...and we want to select the outside of the torus => set to true
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ENDIF
      ENDIF
    ENDDO
    !
  END SELECT
  !
  !
CASE('prop','property')
  !Select atoms according to the given property
  IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
  IF( .NOT. ALLOCATED(AUXNAMES) .OR. .NOT.ALLOCATED(AUX) ) THEN
    !No auxiliary property defined => abort
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2729,(/""/),(/0.d0/))
    !
  ELSE
    !Search the index j for the given property
    i=0
    j=0
    DO WHILE(j==0)
      i=i+1
      IF( TRIM(ADJUSTL(AUXNAMES(i))) == TRIM(ADJUSTL(region_geom)) ) THEN
        j=i
      ENDIF
    ENDDO
    !
    IF( j<=0 ) THEN
      !No such property is defined => abort
      nwarn=nwarn+1
      CALL ATOMSK_MSG(2730,(/region_geom/),(/0.d0/))
      !
    ELSE
      !The given property does exist
      ALLOCATE( SELECT( SIZE(P,1) ) )
      SELECT(:) = .FALSE.
      !
      IF( region_dir(1:3)=="min" ) THEN
        !Select atom(s) whose property takes the smallest value
        !First loop to determine the smallest value of that property
        tempreal = 1.d12
        DO i=1,SIZE(AUX,1)
          IF( among .AND. prevSELECT(i) ) THEN
            IF( AUX(i,j)<tempreal ) THEN
              tempreal = AUX(i,j)
            ENDIF
          ENDIF
        ENDDO
        !Second loop to select atoms that have this value (relative error +/-10^-6)
        DO i=1,SIZE(AUX,1)
          IF( DABS((AUX(i,j)-tempreal)/tempreal)<=1.d-6 ) THEN
            SELECT(i) = .TRUE.
            Nselect = Nselect+1
          ENDIF
        ENDDO
        !
      ELSEIF( region_dir(1:3)=="max" ) THEN
        !Select atom(s) whose property takes the largest value
        !First loop to determine the smallest value of that property
        tempreal = -1.d12
        DO i=1,SIZE(AUX,1)
          IF( among .AND. prevSELECT(i) ) THEN
            IF( AUX(i,j)>tempreal ) THEN
              tempreal = AUX(i,j)
            ENDIF
          ENDIF
        ENDDO
        !Second loop to select atoms that have this value (relative error +/-10^-6)
        DO i=1,SIZE(AUX,1)
          IF( DABS((AUX(i,j)-tempreal)/tempreal)<=1.d-6 ) THEN
            SELECT(i) = .TRUE.
            Nselect = Nselect+1
          ENDIF
        ENDDO
        !
      ELSEIF( region_2(1)>2.d0 ) THEN
        !User gave a range of values
        !Select atoms whose property is in the given range
        DO i=1,SIZE(AUX,1)
          IF( AUX(i,j)>=region_1(1) .AND. AUX(i,j)<=region_1(2) ) THEN
            SELECT(i) = .TRUE.
            Nselect = Nselect+1
          ENDIF
        ENDDO
        !
      ELSE
        !All other cases: user gave only one value (and region_2(1)==0)
        DO i=1,SIZE(AUX,1)
          IF( DABS(AUX(i,j)-region_1(1)) < 1.d-12 ) THEN
            SELECT(i) = .TRUE.
            Nselect = Nselect+1
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDIF
  !
  !
CASE('random','rand','random%','rand%')
  !Select N atoms of the given species at random
  IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
  !
  IF( region_1(1)<=0.d0 ) THEN
    !The user gave a zero or negative number of atoms to select => warning and exit
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2746,(/rand_sp/),(/DBLE(rand_N)/))
    GOTO 1000
  ENDIF
  !
  ALLOCATE( SELECT( SIZE(P,1) ) )
  SELECT(:) = .FALSE.
  snumber = 0.d0
  rand_N  = NINT(region_1(1))
  rand_sp = ADJUSTL(region_geom)
  !Check if rand_sp contains an atom species
  IF( LEN_TRIM(rand_sp)<=2 ) THEN
    species = TRIM(ADJUSTL(rand_sp))
    CALL ATOMNUMBER(species,snumber)
  ENDIF
  !
  IF( NINT(snumber)==0 ) THEN
    species = ''
    !rand_sp does not contain an atom species
  ENDIF
  !
  i=0
  sp_N = 0
  IF(snumber>0.1d0) THEN
    !rand_sp contains a species: only atoms of that species must be selected
    !Count how many atoms of the given species exist in P
    DO i=1,SIZE(P,1)
      IF( NINT(P(i,4))==NINT(snumber) ) THEN
        sp_N = sp_N+1
      ENDIF
    ENDDO
  ELSE
    !rand_sp is "all" or "any": atoms can be selected regardless of their species
    sp_N = SIZE(P,1)
  ENDIF
  !
  !If region_side contains "%" then the user wants to select a given percentage of atoms
  !NOTE: the number was already divided by 100, don't do it again
  IF( INDEX(region_side,"%") > 0 ) THEN
    rand_N = NINT( region_1(1)*DBLE(sp_N) )
  ENDIF
  !
  WRITE(msg,*) 'sp_N, rand_N = ', sp_N, rand_N
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  IF( sp_N==0 ) THEN
    !no such atom in the system => exit
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2723,(/rand_sp/),(/0.d0/))
    !
  ELSEIF( sp_N < rand_N ) THEN
    !The user asked to select more atoms than there are in the system.
    !Warn about that and adjust the numbers
    nwarn=nwarn+1
    CALL ATOMSK_MSG(2745,(/''/),(/DBLE(sp_N)/))
    rand_N = sp_N
    !This becomes simple, there is no randomness possible
    IF(snumber>0.1d0) THEN
      !All atoms of the given species must be selected
      SELECT(:) = .FALSE.
      DO i=1,SIZE(P,1)
        IF( NINT(P(i,4)) == NINT(snumber) ) THEN
          SELECT(i) = .TRUE.
        ENDIF
      ENDDO
    ELSE
      !All atoms must be selected
      IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
    ENDIF
    Nselect = sp_N
    !
  ELSE
    !At last, this is what this part was written for:
    !rand_N atoms must be selected at random.
    !
    !If rand_N < 0.5*sp_N then generate a random list of atoms to select.
    !Otherwise it is faster to generate a list of atoms to be un-selected.
    IF( DBLE(rand_N) <= 0.5d0*DBLE(sp_N) ) THEN
      toselect = .TRUE.
    ELSE
      toselect = .FALSE.
    ENDIF
    !
    !Generate random numbers
    IF(toselect) THEN
      CALL GEN_NRANDNUMBERS( rand_N , randarray )
    ELSE
      CALL GEN_NRANDNUMBERS( sp_N-rand_N , randarray )
    ENDIF
    !
    !randarray now contains random real numbers between 0.d0 and 1.d0
    !Sort them by increasing values (this is bubble sort algorithm)
    DO j=1,SIZE(randarray)
      DO i=j+1,SIZE(randarray)
        !If element i is smaller than element j, swap them
        IF( randarray(i) < randarray(j) ) THEN
          tempreal = randarray(i)
          randarray(i) = randarray(j)
          randarray(j) = tempreal
        ENDIF
      ENDDO
    ENDDO
    IF(verbosity==4) THEN
      WRITE(msg,*) 'randarray (sorted):', SIZE(randarray), ' entries'
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      DO i=1,SIZE(randarray)
        WRITE(msg,*) '     ', randarray(i)
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDDO
    ENDIF
    !
    !Use it to generate the indices of the atoms that will be (un-)selected
    ALLOCATE( atomindices(SIZE(randarray)) )
    atomindices(:) = 0
    DO i=1,SIZE(randarray)
      !
      atomrank = NINT(randarray(i)*sp_N)
      !
      IF( snumber>0.1d0 ) THEN
        !Look for the atomrank-th atom of the given <species>
        !and set atomindices(i) to its index
        k=0
        DO j=1,SIZE(P,1)
          IF( NINT(P(j,4))==NINT(snumber) ) THEN
            k=k+1
            IF(k==atomrank) THEN
              atomindices(i) = j
              EXIT
            ENDIF
          ENDIF
        ENDDO
        !
      ELSE
        !Don't bother with atom species: look for the atomrank-th atom
        !and set atomindices(i) to its index
        k=0
        DO j=1,SIZE(P,1)
          k=k+1
          IF(k==atomrank) THEN
            atomindices(i) = j
            EXIT
          ENDIF
        ENDDO
      ENDIF
      !
      !Refuse illegal indices
      IF( atomindices(i)==0 ) atomindices(i) = 1
      !
      !NOTE: There is a probability that the same index appears twice
      !     in atomindices(:). Of course this probability is expected to be small
      !     if rand_N << sp_N, but still the possibility should not be overlooked.
      !     As a result we have to check for duplicate indices.
      !     If the current index already exists in atomindices(:),
      !     then increase it by 1 until it is different from all other indices
      IF(i>1) THEN
        161 CONTINUE
        DO j=1,i-1  !loop on all values in atomindices(:)
          IF ( atomindices(i)==atomindices(j) ) THEN
            WRITE(msg,*) 'duplicate index in atomindices: ', &
                      & atomindices(i), atomindices(j)
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
            !
            atomindices(i) = atomindices(i)+1
            !
            IF( atomrank>=sp_N ) THEN
              atomindices(i) = 1
            ENDIF
            IF( snumber>0.1d0 ) THEN
              !Increase index until finding the next atom of the given species
              DO WHILE( NINT(P(atomindices(i),4)).NE.NINT(snumber) )
                atomindices(i) = atomindices(i)+1
                IF( atomindices(i)>SIZE(P,1) ) THEN
                  atomindices(i) = 1
                ENDIF
              ENDDO
            ENDIF
            !
            WRITE(msg,*) '                     new index: ', atomindices(i)
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
            !Index of atomindices(i) has changed:
            !go back to checking if it exists elsewhere in atomindices(:)
            GOTO 161
            !
          ENDIF
          !
        ENDDO !atomindices(:)
      ENDIF
    ENDDO !i
    !
    IF(verbosity==4) THEN
      WRITE(msg,*) 'indices of random atoms:'
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      DO i=1,SIZE(atomindices)
        WRITE(msg,*) '     ', atomindices(i)
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDDO
    ENDIF
    !
    !Save the selection into SELECT
    IF(toselect) THEN
      !atomindices(:) contains the indices of atoms to select
      SELECT(:) = .FALSE.
      DO i=1,SIZE(atomindices)
        SELECT(atomindices(i)) = .TRUE.
        Nselect=Nselect+1
      ENDDO
    ELSE
      !atomindices(:) contains the indices of atoms to un-select
      Nselect = sp_N
      SELECT(:) = .FALSE.
      IF( snumber<0.1d0 ) THEN
        !Select all atoms
        SELECT(:) = .TRUE.
      ELSE
        !Select all atoms of the given species
        DO i=1,SIZE(P,1)
          IF( NINT(P(i,4))==NINT(snumber) ) THEN
            SELECT(i) = .TRUE.
          ENDIF
        ENDDO
      ENDIF
      DO i=1,SIZE(atomindices)
        !Un-select atoms with the given indices
        SELECT(atomindices(i)) = .FALSE.
        Nselect=Nselect-1
      ENDDO
    ENDIF
    !
  ENDIF
  !
  !
CASE("neigh","neighbors")
  !The neighbors of the atom with the given index must be searched
  Nselect=0
  exceeds100 = .FALSE.
  IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
  ALLOCATE( SELECT( SIZE(P,1) ) )
  SELECT(:) = .FALSE.
  !region_1(1) = number of neighbors, or cutoff radius for neighbor search
  !region_geom = species of neighboring atoms
  !region_1(2) = index of central atom (=atom whose neighbors must be found)
  !
  atomrank = NINT( region_1(2) )
  IF( atomrank>SIZE(P,1) .OR. atomrank<=0 ) THEN
    !index provided by user is out-of-bound => display error and abandon ship
    nerr = nerr+1
    CALL ATOMSK_MSG(1811,(/""/),(/DBLE(atomrank)/))
    Nselect = 0
    GOTO 1000
  ENDIF
  !
  snumber=-1.d0
  !try to recognize if region_geom contains a recognizable species
  species = TRIM(ADJUSTL(region_geom))
  CALL ATOMNUMBER(species,snumber)
  !
  !Find how many different species exist in the system
  CALL FIND_NSP(P(:,4),aentries)
  !
  IF( snumber>0.1d0 ) THEN
    !Only atoms of a given species must be selected
    !Check if the desired atom species exists in this system
    k=0
    DO i=1,SIZE(aentries,1)
      IF( DABS(aentries(i,1)-snumber)<1.d-12 ) THEN
        k=1  !found it
      ENDIF
    ENDDO
    !
    IF( k==0 ) THEN
      !The species that should be selected does not exist in this system
      !=> display warning and exit
      CALL ATOMSPECIES(snumber,species)
      CALL ATOMSK_MSG(2723,(/species/),(/0.d0/))
      GOTO 1000
    ENDIF
    !
    IF( SIZE(aentries,1)==1 ) THEN
      !There is only one atom species and it must be selected => ignore species
      snumber=-1.d0
      !
    ELSEIF( SIZE(aentries,1)>1 ) THEN
      !There are several atom species
      !Save positions of atoms of the species that must be selected
      DO i=1,SIZE(aentries,1)
        IF( .NOT.ALLOCATED(Q) .AND. DABS(aentries(i,1)-snumber)<1.d-12 ) THEN
          !aentries(i,1) is the species that must be selected
          !aentries(i,2) is the number of atoms of that species
          ALLOCATE( Q( NINT(aentries(i,2)),4) )
          Q(:,:) = 0.d0
          !Save positions of atoms of given species into Q
          k=0
          DO j=1,SIZE(P,1)
            IF( DABS(P(j,4)-snumber)<1.d-12 ) THEN
              k=k+1
              Q(k,:) = P(j,:)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      !
    ENDIF
  ENDIF
  !
  !Search for the neighbors
  !NOTE: the central atom must not be selected (only its neighbors are selected)
  IF( region_1(1)==0.d0 ) THEN
    !Search for the first nearest neighbors (their number or distance are unknown)
    IF( snumber>0.1d0 .AND. ALLOCATED(Q) .AND. SIZE(Q,1)>0 ) THEN
      CALL FIND_1NN(H,Q,P(atomrank,:),V_NN,Nlist,exceeds100)
      IF(ALLOCATED(Q)) DEALLOCATE(Q)
    ELSE
      CALL FIND_1NN(H,P,P(atomrank,:),V_NN,Nlist,exceeds100)
    ENDIF
    !
    IF( exceeds100 ) THEN
      !This atom has more than 100 neighbors => display a warning
      nwarn=nwarn+1
      CALL ATOMSK_MSG(4705,(/''/),(/DBLE(atomrank)/))
    ENDIF
    !
    !The indices of neighbors are stored in Nlist => select these atoms
    IF( SIZE(Nlist)>0 ) THEN
      IF( snumber>0.1d0 ) THEN
        !Count atoms of given species, select those who appear in Nlist(:)
        k=0
        DO i=1,SIZE(P,1)
          IF( DABS(P(i,4)-snumber)<1.d-12 ) THEN
            !This atom is of the desired species => increase counter
            k=k+1
            !Search if this index is in Nlist(:)
            DO j=1,SIZE(Nlist)
              IF( k==Nlist(j) .AND. k.NE.atomrank ) THEN
                !Atom #i is in the list => select it
                SELECT(i) = .TRUE.
                Nselect = Nselect+1
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        !
      ELSE
        DO i=1,SIZE(Nlist)
          IF( Nlist(i).NE.atomrank ) THEN
            SELECT(Nlist(i)) = .TRUE.
            Nselect = Nselect+1
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    !
  ELSEIF( region_1(1)>0.d0 .AND. IS_INTEGER(region_1(1)) ) THEN
    !Search for the N closest neighbors
    k=NINT(region_1(1))
    IF( snumber<0.1d0 .OR. DABS(P(atomrank,4)-snumber)<1.d-12 ) THEN
      !Look for 1 more neighbor because the central atom will be found as its own neighbor
      !(this will be corrected for afterwards)
      k = k+1
    ENDIF
    !
    IF( snumber>0.1d0 .AND. ALLOCATED(Q) .AND. SIZE(Q,1)>0 ) THEN
      CALL FIND_NNN(H,Q,P(atomrank,:),k,V_NN,Nlist,exceeds100)
      IF(ALLOCATED(Q)) DEALLOCATE(Q)
    ELSE
      CALL FIND_NNN(H,P,P(atomrank,:),k,V_NN,Nlist,exceeds100)
    ENDIF
    !
    IF( exceeds100 ) THEN
      !This atom has more than 100 neighbors => display a warning
      nwarn=nwarn+1
      CALL ATOMSK_MSG(4705,(/''/),(/DBLE(atomrank)/))
    ENDIF
    !
    !The indices of neighbors are stored in Nlist => select these atoms
    IF( SIZE(Nlist)>0 ) THEN
      IF( snumber>0.1d0 ) THEN
        !Count atoms of given species, select those who appear in Nlist(:)
        k=0
        DO i=1,SIZE(P,1)
          IF( DABS(P(i,4)-snumber)<1.d-12 ) THEN
            !This atom is of the desired species => increase counter
            k=k+1
            IF( i.NE.atomrank ) THEN
              !Do not select central atom
              !Search if this index is in Nlist(:)
              DO j=1,SIZE(Nlist)
                IF( k==Nlist(j) ) THEN
                  !Atom #i is in the list => select it
                  SELECT(i) = .TRUE.
                  Nselect = Nselect+1
                ENDIF
              ENDDO
            ENDIF
          ENDIF
        ENDDO
        !
      ELSE
        !Easy case: select atoms with indices in Nlist(:)
        DO i=1,SIZE(Nlist)
          IF( Nlist(i).NE.atomrank ) THEN
            SELECT(Nlist(i)) = .TRUE.
            Nselect = Nselect+1
          ENDIF
        ENDDO
      ENDIF
    ENDIF
    !
    !
  ELSE
    !Search for all neighbors in the radius R=DABS(region_1(1))
    !This is similar to "-select in sphere"
    !P(atomrank,:) = center of the sphere
    !region_1(1) = radius of the sphere
    DO i=1,SIZE(P,1)
      distance = VECLENGTH( P(i,1:3) - P(atomrank,1:3) )
      IF( distance < DABS(region_1(1)) ) THEN
        !Atom is inside the sphere
        IF( snumber<0.1d0 .OR. DABS(P(i,4)-snumber)<1.d-12 ) THEN
          !We don't care about the atom species, or atom is of the right species => select it
          IF( i.NE.atomrank ) THEN
            SELECT(i) = .TRUE.
            Nselect = Nselect+1
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDIF
  !
  !
CASE("list")
  !A list of atom indices is read from the file "region_geom"
  IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
  Nselect=0
  !Check that the file actually exists
  CALL CHECKFILE(region_geom,"read")
  OPEN(UNIT=30,FILE=region_geom,FORM='FORMATTED')
  DO
    READ(30,'(a128)',ERR=281,END=281) temp
    temp = ADJUSTL(temp)
    IF( LEN_TRIM(temp)>0 .AND. temp(1:1).NE."#" ) THEN
      READ(temp,*,ERR=280,END=280) i
      IF( i>0 .AND. i<=SIZE(P,1) ) THEN
        IF( .NOT. ALLOCATED(SELECT) ) THEN
          ALLOCATE( SELECT(SIZE(P,1)) )
          SELECT(:) = .FALSE.
        ENDIF
        !Select atom #i
        SELECT(i) = .TRUE.
        Nselect = Nselect+1
      ELSE
        !Atom index is out-of-bounds
        nwarn = nwarn+1
        CALL ATOMSK_MSG(2742,(/""/),(/DBLE(i)/))
      ENDIF
    ENDIF
    280 CONTINUE
  ENDDO
  281 CONTINUE
  CLOSE(30)
  !
  !
CASE('grid')
  !Elements of a grid are read from the file "region_geom"
  !This file is a text file that can have one of three formats:
  !
  !(1) It may contain one line with three integers "NX NY NZ",
  !    indicating the number of grid elements along each direction
  !    (in that case all elements have the same size),
  !    followed by lines of 0 and 1 indicating if atoms inside each element
  !    are selected (1) or not (0).
  !
  !(2) It may contain several lines with values "x y z 1|0".
  !    The (x,y,z) are actually the positions of nodes and are used
  !    to define Voronoi polyhedra. The last value (1 or 0) indicates
  !    if atoms inside a grid element are selected (1) or not (0)
  !
  !(3) It may contain rows of 0 and 1 that represent a 2-D grid.
  !    The values may be separated by blank spaces or contiguous.
  !
  gridformat = 2  !by default, assume format "x y z 1|0"
  IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
  Nselect=0
  !Check that the file actually exists
  CALL CHECKFILE(region_geom,"read")
  OPEN(UNIT=30,FILE=region_geom,FORM='FORMATTED')
  !
  !Parse file a first time to count total number of elements in the grid
  Ngrid=0
  Nperline=0
  k=0
  line = 1
  DO
    READ(30,'(a4096)',ERR=294,END=294) temp
    Nperline=LEN_TRIM(temp)
    temp = ADJUSTL(temp)
    k=k+1
    IF( k==1 ) THEN
      !First line of data => determine the format of the file
      !Check if line contains only numbers
      IF( temp(1:6)=="FINITE" .OR. temp(1:6)=="Finite" .OR. temp(1:6)=="finite" &
        & .OR. temp(1:2)=="FE" .OR. temp(1:2)=="fe" ) THEN
        !Finite-element grid: lines of "xyz 0|1" should follow
        gridformat = 2
      ELSEIF( SCAN(temp,"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz?,;:!§%*$£=+-^&~'#{}[]()<>/") == 0 ) THEN
        !Try to read three integers
        READ(temp,*,ERR=292,END=292) a1, a2, a3
        !Success => we have NX=a1, NY=a2, NZ=a3
        gridformat=1
        IF( a1*a2*a3 <= 0 ) THEN
          !cannot have a grid with zero dimension
          !=> maybe those 3 integers are part of a 2-D grid
          GOTO 292
        ELSE
          Ngrid = a1*a2*a3
          GOTO 294
        ENDIF
        !
        292 CONTINUE
        !Failed to read 3 integers, and failed to read 3 real numbers
        !=> assume 2-D grid made of random characters
        gridformat = 3
        !Determine max. length of a line and number of lines
        a2=1  !one line was read already
        DO
          READ(30,'(a4096)',ERR=293,END=293) temp
          a2 = a2+1
          IF( LEN_TRIM(temp) > Nperline ) THEN
            Nperline = LEN_TRIM(temp)
          ENDIF
        ENDDO
        293 CONTINUE
        !Now we know the number of lines, and the number of elements per line
        a1 = Nperline
        a3 = 1
        Ngrid = a1*a2
      ENDIF
      !
    ELSE
      !Not first line of data => read line assuming same format
      IF( gridformat==2 ) THEN
        !Verify that 3 real numbers + 1 integer can be read
        READ(temp,*,ERR=801,END=801) V1, V2, V3, i
        Ngrid=Ngrid+1
      ELSEIF( gridformat==3 ) THEN
        Ngrid=Ngrid+Nperline
      ENDIF
    ENDIF
    !
  ENDDO
  294 CONTINUE
  IF( verbosity==4 ) THEN
    IF( gridformat==1 ) THEN
      WRITE(msg,*) "Format of the grid file: (1) NX NY NZ + lines of 0 and 1"
    ELSEIF( gridformat==2 ) THEN
      WRITE(msg,*) "Format of the grid file: (2) lines of  x y z 1|0"
    ELSE
      WRITE(msg,*) "Format of the grid file: (3) lines of ASCII characters forming a 2-D pattern"
    ENDIF
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,*) "Number of grid elements: ", Ngrid
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
  !
  IF( Ngrid>0 ) THEN
    ALLOCATE(GRID(Ngrid,4))
    GRID(:,:) = 0.d0
    !
    !Now store all positions of grid elements into the array GRID
    !Also, count the number of grid elements that must be selected
    Ngrid=0
    Nselect=0
    REWIND(30)
    !
    IF( gridformat==1 ) THEN
      !First line contains NX, NY, NZ (it was already read, so there should not be any error)
      READ(30,*,ERR=801,END=801) temp
      !Read lines of 0 and 1
      !Each new value corresponds to a new grid element
      !Note that in that case, positions of grid elements are not used
      Ngrid=0 !grid element counter
      DO
        READ(30,'(a4096)',ERR=295,END=295) region_geom
        region_geom = ADJUSTL(region_geom)
        IF( LEN_TRIM(region_geom)>0 .AND. region_geom(1:1).NE."#" ) THEN
          !Strip all blank spaces
          j=SCAN(TRIM(region_geom)," ")
          DO WHILE( j>0 )
            region_geom = ADJUSTL( region_geom(:j-1)//ADJUSTL(region_geom(j+1:)) )
            j=SCAN(TRIM(region_geom)," ")
          ENDDO
          !Read values (should be only 0 or 1)
          DO i=1,LEN_TRIM(region_geom)
            Ngrid=Ngrid+1
            IF( Ngrid<=SIZE(GRID,1) ) THEN
              READ(region_geom(i:i),'(i1)',ERR=800,END=800) j
              IF( j==0 ) THEN
                GRID(Ngrid,4) = 0.d0
              ELSE
                GRID(Ngrid,4) = 1.d0
              ENDIF
            ENDIF
          ENDDO !i
        ENDIF
      ENDDO
      GOTO 295
      !
    ELSEIF( gridformat==2 ) THEN
      !First line contains "FE" or similar (it was already read, so there should not be any error)
      READ(30,*,ERR=801,END=801) temp
      !Read position and status (0 or 1) of each grid element
      DO
        READ(30,'(a4096)',ERR=295,END=295) temp
        temp = ADJUSTL(temp)
        IF( LEN_TRIM(temp)>0 .AND. temp(1:1).NE."#" ) THEN
          Ngrid = Ngrid+1
          READ(temp,*,ERR=801,END=801) GRID(Ngrid,1), GRID(Ngrid,2), GRID(Ngrid,3), GRID(Ngrid,4)
          IF( GRID(Ngrid,4)>0.1d0 ) THEN
            Nselect=Nselect+1
          ENDIF
        ENDIF
      ENDDO
      !
    ELSEIF( gridformat==3 ) THEN
      !Lines may contain arbitrary ASCII characters
      !Blank spaces and 0 count as "0", any other character as "1"
      WRITE(msg,*) "Number of elements per line: ", Nperline
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      Ngrid = 0
      !Read lines, one at a time
      DO
        READ(30,'(a4096)',ERR=295,END=295) temp
        !Parse line, one character at a time
        !Note: it is assumed that all lines contain Nperline characters
        !     Shorter lines will be completed with zeros if necessary
        DO i=1,Nperline
          Ngrid = Ngrid+1
          IF( Ngrid<=SIZE(GRID,1) ) THEN
            IF( temp(i:i)==" " .OR. temp(i:i)=="0" ) THEN
              GRID(Ngrid,4) = 0.d0
            ELSE
              GRID(Ngrid,4) = 1.d0
            ENDIF
          ELSE
            EXIT
          ENDIF
        ENDDO
      ENDDO
    ENDIF
    295 CONTINUE
    CLOSE(30)
    !
    IF( verbosity==4 ) THEN
      WRITE(msg,*) "Finite-element grid:"
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      DO i=1,SIZE(GRID,1)
        WRITE(msg,'(4f9.3,i3)') GRID(i,1), GRID(i,2), GRID(i,3), GRID(i,4)
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDDO
    ENDIF
    !
    ALLOCATE(SELECT(SIZE(P,1)))
    SELECT(:) = .FALSE.
    Nselect=0
    !
    IF( gridformat == 1 .OR. gridformat==3 ) THEN
      !All elements have the same size, define dimensions of an element
      !(reminder: a1=NX, a2=NY, a3=NZ)
      region_1(1) = H(1,1)/DBLE(a1)
      region_1(2) = H(2,2)/DBLE(a2)
      region_1(3) = H(3,3)/DBLE(a3)
      !
      DO i=1,SIZE(P,1)
        !Parse all grid elements
        Ngrid=0 !counter for grid elements
        DO n=1,a3  !outer loop along Z
          DO m=1,a2
            DO l=1,a1  !inner loop along X
              Ngrid=Ngrid+1
              !Find if atom #i is in this grid element
              IF( P(i,1)>=DBLE(l-1)*region_1(1) .AND. P(i,1)<DBLE(l)*region_1(1) .AND. &
                & P(i,2)>=DBLE(m-1)*region_1(2) .AND. P(i,2)<DBLE(m)*region_1(2) .AND. &
                & P(i,3)>=DBLE(n-1)*region_1(3) .AND. P(i,3)<DBLE(n)*region_1(3)      ) THEN
                !Atom #i is inside that grid element
                !Now check if atoms inside current grid elements must be selected
                IF( NINT(GRID(Ngrid,4)) > 0 ) THEN
                  SELECT(i) = .TRUE.
                  Nselect=Nselect+1
                  !Stop parsing grid elements, jump to next atom
                  GOTO 297
                ENDIF
              ENDIF
            ENDDO !l
          ENDDO !m
        ENDDO !n
        297 CONTINUE
      ENDDO
      !
    ELSE  !i.e. gridformat==2
      !Check if coordinates of GRID elements are already reduced or not
      CALL FIND_IF_REDUCED(GRID,isreduced)
      !If reduced, convert them to Cartesian
      IF( isreduced ) THEN
        CALL FRAC2CART(GRID,H)
      ENDIF
      !
      !Generate a neighbor list for GRID elements
      !Guess radius for neighbor search, do not let it fall below 8 A
      distance = 2.d0*MAX(VECLENGTH(H(1,:)),VECLENGTH(H(2,:)),VECLENGTH(H(3,:))) / DSQRT(DBLE(Ngrid))
      CALL NEIGHBOR_LIST(H,GRID,distance,NeighList)
      !
      IF( verbosity==4 ) THEN
        DO i=1,SIZE(NeighList,1)
          WRITE(msg,*) "Grid neighbors:"
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          WRITE(msg,'(i3,a3,10(1X,f9.3))') i, " | ", (NeighList(i,j),j=1,10)
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        ENDDO
      ENDIF
      !
      IF( Nselect < SIZE(GRID,1)/2 ) THEN
        !Unselect all atoms, and parse only grid elements that must be selected
        SELECT(:) = .FALSE.
        Nselect=0
        !
        DO i=1,SIZE(P,1)
          DO Ngrid=1,SIZE(GRID,1)
            IF( DABS(GRID(Ngrid,4)) > 0.1d0 ) THEN
              !Find neighbor elements of element #Ngrid
              CALL NEIGHBOR_POS(H,GRID,GRID(Ngrid,1:3),NeighList(Ngrid,:),ALLOCATED(NeighList),distance,PosList)
              !
              IF( ALLOCATED(PosList) .AND. SIZE(PosList,1)>0 ) THEN
                !PosList(:,1:3) positions of neighbors
                !PosList(:,4)   distances of neighbors
                !PosList(:,5)   index of neighbor (=index of element in GRID)
                !
                !Check if atom #i is inside element #Ngrid
                keep=.TRUE.
                DO k=1,SIZE(PosList,1)  !Loop on all neighbors of element #Ngrid
                  !Compute vector between the neighboring element #k and element #Ngrid
                  region_1(:) = PosList(k,1:3) - GRID(Ngrid,1:3)
                  !Compute vector between atom and element #Ngrid
                  region_2(:) = P(i,1:3) - GRID(Ngrid,1:3)
                  IF( VEC_PLANE(region_1,VECLENGTH(region_1),region_2) > -0.001d0 ) THEN
                    !Atom is above this plane of cut, hence out of the grid element #Ngrid
                    !=> un-select it and exit the loop on k
                    keep = .FALSE.
                    EXIT
                  ENDIF
                ENDDO
                !
                IF( keep ) THEN
                  !Atom #i is inside grid element #Ngrid => select it
                  SELECT(i) = .TRUE.
                  Nselect=Nselect+1
                  !An atom can only be inside one grid element => exit loop on Ngrid
                  EXIT
                ENDIF
                !
              ELSE
                !PRINT*, "WARNING: grid element #", Ngrid, " has no neighbor"
              ENDIF
              !
            ENDIF
            !
          ENDDO !loop on Ngrid
        ENDDO  !loop on i
      ELSE
        !Select all atoms, and parse only grid elements that must be un-selected
        SELECT(:) = .TRUE.
        Nselect=SIZE(P,1)
        !
        DO i=1,SIZE(P,1)
          DO Ngrid=1,SIZE(GRID,1)
            IF( DABS(GRID(Ngrid,4)) < 0.1d0 ) THEN
              !Find neighbor elements of element #Ngrid
              CALL NEIGHBOR_POS(H,GRID,GRID(Ngrid,1:3),NeighList(Ngrid,:),ALLOCATED(NeighList),distance,PosList)
              !
              IF( ALLOCATED(PosList) .AND. SIZE(PosList,1)>0 ) THEN
                !PosList(:,1:3) positions of neighbors
                !PosList(:,4)   distances of neighbors
                !PosList(:,5)   index of neighbor (=index of element in GRID)
                !
                !Check if atom #i is inside element #Ngrid
                keep=.TRUE.
                DO k=1,SIZE(PosList,1)  !Loop on all neighbors of element #Ngrid
                  !Compute vector between the neighboring element #k and element #Ngrid
                  region_1(:) = PosList(k,1:3) - GRID(Ngrid,1:3)
                  !Compute vector between atom and element #Ngrid
                  region_2(:) = P(i,1:3) - GRID(Ngrid,1:3)
                  IF( VEC_PLANE(region_1,VECLENGTH(region_1),region_2) > -0.001d0 ) THEN
                    !Atom is above this plane of cut, hence out of the grid element #Ngrid
                    !=> un-select it and exit the loop on k
                    keep = .FALSE.
                    EXIT
                  ENDIF
                ENDDO
                !
                IF( keep ) THEN
                  !Atom #i is inside grid element #Ngrid => unselect it
                  SELECT(i) = .FALSE.
                  Nselect=Nselect-1
                  !An atom can only be inside one grid element => exit loop on Ngrid
                  EXIT
                ENDIF
                !
              ELSE
                !PRINT*, "WARNING: grid element #", Ngrid, " has no neighbor"
              ENDIF
              !
            ENDIF
            !
          ENDDO !loop on Ngrid
        ENDDO  !loop on i
        !
      ENDIF
      !
    ENDIF
    !
  ENDIF
  !
  !
CASE('stl','STL')
  !Select atoms inside polygons (triangles) defined in a STL file
  !NOTE: in order to determine if an atom at coordinates (x,y,z) is inside
  !     or outside of the 3-D shape, the program follows a line starting at
  !     (x,y,0), progressively increasing Z, and incrementing a variable n
  !     every time it crosses the surface of a triangle. In the end, if
  !     n is even then the atom is inside the 3-D shape, otherwise it is outside.
  !
  IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
  !
  !Define increment to parse along Z
  !NOTE: if two triangles of the 3D surface are separated by less than
  !     this distance, the detection will fail
  dz = 10.d0  !MAX ( 0.1d0 , 0.01d0*VECLENGTH(H(:,3)) )
  !
  !Define
  m = NINT( VECLENGTH(H(:,3)) / dz ) + 1
  !
  !Read triangles from STL file
  CALL READ_STL(region_geom,triangles)
  !
  IF( ALLOCATED(triangles) .AND. SIZE(triangles,1) > 3 ) THEN
    CALL ATOMSK_MSG(2142,(/""/),(/DBLE(SIZE(triangles,1))/))
    !
    !The array "triangles" should contain all triangles defining the 3-D object
    !It is a N x 12 array (where N is the number of triangles):
    !  triangles(:,1:3)     vector normal to the triangle surface
    !  triangles(:,4:6)     position of 1st vertex of triangle
    !  triangles(:,7:9)     position of 2nd vertex of triangle
    !  triangles(:,10:12)   position of 3rd vertex of triangle
    !
    !Sort triangles by increasing Z coordinates
    !NOTE: if the file respects the STL standard, it is already the case.
    !     This is just a precaution to make sure it is the case.
    !NOTE2: here we sort only the Z coordinates of the 1st vertices
    CALL QUICKSORT( triangles(:,:) , 6 , 'up  ' , newindex )
    IF(ALLOCATED(newindex)) DEALLOCATE(newindex)
    !
    !Determine bounding box of 3-D model
    xmin = MIN( MINVAL(triangles(:,4)) , MINVAL(triangles(:,7)) , MINVAL(triangles(:,10)) )
    xmax = MAX( MAXVAL(triangles(:,4)) , MAXVAL(triangles(:,7)) , MAXVAL(triangles(:,10)) )
    ymin = MIN( MINVAL(triangles(:,5)) , MINVAL(triangles(:,8)) , MINVAL(triangles(:,11)) )
    ymax = MAX( MAXVAL(triangles(:,5)) , MAXVAL(triangles(:,8)) , MAXVAL(triangles(:,11)) )
    zmin = MIN( MINVAL(triangles(:,6)) , MINVAL(triangles(:,9)) , MINVAL(triangles(:,12)) )
    zmax = MAX( MAXVAL(triangles(:,6)) , MAXVAL(triangles(:,9)) , MAXVAL(triangles(:,12)) )
    WRITE(msg,'(a14,6f6.2)') "BOUNDING BOX: ", xmin, xmax, ymin, ymax, zmin, zmax
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !
    !Rescale the 3-D model to the dimensions of current simulation box
    !conserving the proportions of the 3-D model
    !tempreal is the rescaling factor
    tempreal = MIN( VECLENGTH(H(:,1)) / (xmax-xmin) ,  &
                  & VECLENGTH(H(:,2)) / (ymax-ymin) ,  &
                  & VECLENGTH(H(:,3)) / (zmax-zmin)  )
    triangles(:,4)  = triangles(:,4)  * tempreal
    triangles(:,7)  = triangles(:,7)  * tempreal
    triangles(:,10) = triangles(:,10) * tempreal
    triangles(:,5)  = triangles(:,5)  * tempreal
    triangles(:,8)  = triangles(:,8)  * tempreal
    triangles(:,11) = triangles(:,11) * tempreal
    triangles(:,6)  = triangles(:,6)  * tempreal
    triangles(:,9)  = triangles(:,9)  * tempreal
    triangles(:,12) = triangles(:,12) * tempreal
    !
    !Update bounding box of 3-D model
    xmin = MIN( MINVAL(triangles(:,4)) , MINVAL(triangles(:,7)) , MINVAL(triangles(:,10)) )
    xmax = MAX( MAXVAL(triangles(:,4)) , MAXVAL(triangles(:,7)) , MAXVAL(triangles(:,10)) )
    ymin = MIN( MINVAL(triangles(:,5)) , MINVAL(triangles(:,8)) , MINVAL(triangles(:,11)) )
    ymax = MAX( MAXVAL(triangles(:,5)) , MAXVAL(triangles(:,8)) , MAXVAL(triangles(:,11)) )
    zmin = MIN( MINVAL(triangles(:,6)) , MINVAL(triangles(:,9)) , MINVAL(triangles(:,12)) )
    zmax = MAX( MAXVAL(triangles(:,6)) , MAXVAL(triangles(:,9)) , MAXVAL(triangles(:,12)) )
    WRITE(msg,'(a22,6f6.2)') "UPDATED BOUNDING BOX: ", xmin, xmax, ymin, ymax, zmin, zmax
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !
    IF( region_dir=="center" ) THEN
      !Translate 3-D model so that it is at center of simulation box
      tempreal = xmin + (xmax-xmin)/2.d0 - MAXVAL(H(1,:))/2.d0
      triangles(:,4)  = triangles(:,4) - tempreal
      triangles(:,7)  = triangles(:,7) - tempreal
      triangles(:,10) = triangles(:,10) - tempreal
      tempreal = ymin + (ymax-ymin)/2.d0 - MAXVAL(H(2,:))/2.d0
      triangles(:,5)  = triangles(:,5) - tempreal
      triangles(:,8)  = triangles(:,8) - tempreal
      triangles(:,11) = triangles(:,11) - tempreal
      tempreal = zmin + (zmax-zmin)/2.d0 - MAXVAL(H(3,:))/2.d0
      triangles(:,6)  = triangles(:,6) - tempreal
      triangles(:,9)  = triangles(:,9) - tempreal
      triangles(:,12) = triangles(:,12) - tempreal
      !
      !Update bounding box of 3-D model once again
      xmin = MIN( MINVAL(triangles(:,4)) , MINVAL(triangles(:,7)) , MINVAL(triangles(:,10)) )
      xmax = MAX( MAXVAL(triangles(:,4)) , MAXVAL(triangles(:,7)) , MAXVAL(triangles(:,10)) )
      ymin = MIN( MINVAL(triangles(:,5)) , MINVAL(triangles(:,8)) , MINVAL(triangles(:,11)) )
      ymax = MAX( MAXVAL(triangles(:,5)) , MAXVAL(triangles(:,8)) , MAXVAL(triangles(:,11)) )
      zmin = MIN( MINVAL(triangles(:,6)) , MINVAL(triangles(:,9)) , MINVAL(triangles(:,12)) )
      zmax = MAX( MAXVAL(triangles(:,6)) , MAXVAL(triangles(:,9)) , MAXVAL(triangles(:,12)) )
      WRITE(msg,'(a22,6f6.2)') "UPDATED BOUNDING BOX: ", xmin, xmax, ymin, ymax, zmin, zmax
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDIF
    !
    ALLOCATE( SELECT(SIZE(P,1)) )
    SELECT(:) = .FALSE.
    !
    !Loop on all atoms
    td = (/0,0,1/)  !vector of ray
    DO i=1,SIZE(P,1)
      !
      IF( SIZE(P,1)>500000 .OR. SIZE(triangles,1)>100000 ) THEN
        !If there are many atoms and triangles, display a fancy progress bar
        CALL ATOMSK_MSG(10,(/""/),(/DBLE(i),DBLE(SIZE(P,1))/))
      ENDIF
      !
      IF( P(i,1)<xmin .OR. P(i,1)>xmax .OR. &
        & P(i,2)<ymin .OR. P(i,2)>ymax .OR. &
        & P(i,3)<zmin .OR. P(i,3)>zmax      ) THEN
        !Atom is outside of the bounding box, hence outside of the 3-D shape
        !SELECT(i) = .FALSE.
      ELSE
        !Loop on Z coordinates to check how many vertices are crossed
        !At beginning atom is assumed to be outside the 3-D shape
        keep = .FALSE.
        !
        !Loop on all triangles
        k = 1
        DO WHILE( k <= SIZE(triangles,1) )
          txmin = MIN( triangles(k,4) , triangles(k,7) , triangles(k,10) ) - 0.1d0
          txmax = MAX( triangles(k,4) , triangles(k,7) , triangles(k,10) ) + 0.1d0
          tymin = MIN( triangles(k,5) , triangles(k,8) , triangles(k,11) ) - 0.1d0
          tymax = MAX( triangles(k,5) , triangles(k,8) , triangles(k,11) ) + 0.1d0
          tzmax = MAX( triangles(k,6) , triangles(k,9) , triangles(k,12) ) + 0.1d0
          IF( P(i,3)<tzmax .AND. P(i,1)>txmin .AND. P(i,1)<txmax .AND. &
            & P(i,2)>tymin .AND. P(i,2)<tymax   ) THEN
            !Current triangle appears to be above (or close) to current atom
            !Determine if ray [001] passing through atom #i intersects the triangle #k
            !Note: this part follows the algorithm proposed in this Web site:
            !  http://www.lighthouse3d.com/tutorials/maths/ray-triangle-intersection/
            e1(:) = triangles(k,7:9) - triangles(k,4:6)
            e2(:) = triangles(k,10:12) - triangles(k,4:6)
            !
            th = CROSS_PRODUCT(td,e2)
            ta = DOT_PRODUCT(e1,th)
            !
            IF( DABS(ta) > 1.d-12 ) THEN
              tf = 1.d0/ta
              ts(:) = P(i,1:3) - triangles(k,4:6)
              tu = tf * DOT_PRODUCT(ts,th)
              !
              IF( tu >= 0.d0 .AND. tu <= 1.d0 ) THEN
                !Ray seems to intersect triangle
                tq = CROSS_PRODUCT(ts,e1)
                tv = tf * DOT_PRODUCT(td,tq)
                !
                IF( tv >= 0.d0 .AND. tu+tv <= 1.d0 ) THEN
                  !There is a line or ray intersection
                  ! at this stage we can compute t to find out where
                  ! the intersection point is on the line
                  tt = tf * DOT_PRODUCT(e2,tq)
                  !
                  IF( tt > 1.d-12 ) THEN
                    !The ray intersects the triangle
                    !=> invert status of atom
                    keep = .NOT.keep
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
            !
          ENDIF
          !
!           IF( triangles(k,6)>P(i,3) ) THEN
!             n=0
!           ENDIF
          k = k+1
          !
        ENDDO !end loop on k
        !
        IF( keep ) THEN
          !This atom is inside the 3-D shape
          !PRINT*, "Keep atom #", i
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ENDIF
      ENDIF
      !
    ENDDO !end loop on i
    !
    !Free memory
    IF(ALLOCATED(triangles)) DEALLOCATE(triangles)
    !
  ELSE
    !Array "triangles" is not allocated or has zero size
    nerr = nerr+1
    CALL ATOMSK_MSG(2818,(/''/),(/0.d0/))
    GOTO 1000
  ENDIF
  !
  !
  !
CASE DEFAULT
  IF( LEN_TRIM(region_side)<=2 ) THEN
    !Last case, it should be an atom species: select all atoms of the given species
    IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
    !
    CALL ATOMNUMBER(region_side(1:2),snumber)
    IF( NINT(snumber)==0 ) THEN
      !could not recognize an atom species => error
      nerr=nerr+1
      GOTO 1000
    ELSE
      !Select all atoms of the given species
      ALLOCATE( SELECT( SIZE(P,1) ) )
      SELECT(:) = .FALSE.
      DO i=1,SIZE(P,1)
        IF( NINT(P(i,4))==NINT(snumber) ) THEN
          SELECT(i) = .TRUE.
          Nselect = Nselect+1
        ENDIF
      ENDDO
    ENDIF
  ENDIF
  !
  !
END SELECT
!
!
!
300 CONTINUE
IF( selectmul ) THEN
  !Final array SELECT will be a combination of previous selection and current one
  IF( add ) THEN
    !New atoms must be added to selection
    Nadded = 0
    Nselect = 0
    DO i=1,SIZE(SELECT)
      SELECT(i) = prevSELECT(i) .OR. SELECT(i)
      IF( SELECT(i) ) THEN
        Nselect = Nselect+1
        IF( .NOT.prevSELECT(i) ) THEN
          !Atom was not selected before => count it as added to the selection
          Nadded = Nadded+1
        ENDIF
      ENDIF
    ENDDO
    !
  ELSEIF( rm ) THEN
    !Some atoms must be removed from previous selection
    Nrm = 0
    Nselect = 0
    DO i=1,SIZE(SELECT)
      SELECT(i) = prevSELECT(i) .AND. (.NOT.SELECT(i))
      IF( SELECT(i) ) THEN
        Nselect = Nselect+1
      ELSE
        IF( prevSELECT(i) ) THEN
          !Atom was selected before => count it as removed from the selection
          Nrm = Nrm+1
        ENDIF
      ENDIF
    ENDDO
    !
  ELSEIF( intersect ) THEN
    !Only atoms that belong to both SELECT arrays are kept
    Nselect = 0
    DO i=1,SIZE(SELECT)
      SELECT(i) = prevSELECT(i) .AND. SELECT(i)
      IF( SELECT(i) ) THEN
        !Atom is now selected
        Nselect = Nselect+1
        IF( .NOT.prevSELECT(i) ) THEN
          !Atom was not selected before
          Nadded = Nadded+1
        ENDIF
      ELSE
        !Atom is now un-selected
        IF( prevSELECT(i) ) THEN
          !Atom was selected before
          Nrm = Nrm+1
        ENDIF
      ENDIF
    ENDDO
    !
  ELSEIF( xor ) THEN
    !Only atoms that belong to one SELECT array or the other, but not both, are kept
    Nselect = 0
    DO i=1,SIZE(SELECT)
      SELECT(i) = SELECT(i).NEQV.prevSELECT(i)
      IF( SELECT(i) ) THEN
        !Atom is now selected
        Nselect = Nselect+1
        IF( .NOT.prevSELECT(i) ) THEN
          !Atom was not selected before
          Nadded = Nadded+1
        ENDIF
      ELSE
        !Atom is now un-selected
        IF( prevSELECT(i) ) THEN
          !Atom was selected before
          Nrm = Nrm+1
        ENDIF
      ENDIF
    ENDDO
  ENDIF
ENDIF
CALL ATOMSK_MSG(2078,(/''/),(/DBLE(Nselect),DBLE(Nadded),DBLE(Nrm)/))
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
801 CONTINUE
CALL ATOMSK_MSG(807,(/''/),(/DBLE(line)/))
nerr = nerr+1
!
!
!
1000 CONTINUE
!
IF(ALLOCATED(atomindices)) DEALLOCATE(atomindices)
IF(ALLOCATED(aentries)) DEALLOCATE(aentries)
IF(ALLOCATED(prevSELECT)) DEALLOCATE(prevSELECT)
! IF( Nselect<=0 ) THEN
!   IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
! ENDIF
! IF( Nselect==SIZE(P,1) ) THEN
!   IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
! ENDIF
!
!
!
END SUBROUTINE SELECT_XYZ
!
!
END MODULE select
