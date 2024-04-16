MODULE crystallography
!
!**********************************************************************************
!*  CRYSTALLOGRAPHY                                                               *
!**********************************************************************************
!* This module contains routines related to crystallography, like                 *
!* manipulating Miller indexes.                                                   *
!**********************************************************************************
!* (C) April 2022                                                                 *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 16 April 2024                                    *
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
!* List of subroutines in this module:                                            *
!* INDEX_MILLER        find the indices of a plane from a string                  *
!* INDEX_MILLER_HCP    find the indices of a plane from a string (hcp lattices)   *
!* HKIL2UVW            converts Miller-Bravais [hkil] into Miller [uvw]           *
!* UVW2HKIL            converts Miller [uvw] into Miller-Bravais [hkil]           *
!* MILLER2VEC          translates Miller indices into Cartesian vector            *
!* COMPFORMULA         extracts a compound formula from atom site lists P and AUX *
!**********************************************************************************
!
!
USE comv
USE constants
USE math
USE subroutines
!
IMPLICIT NONE
!
!
CONTAINS
!
!
!********************************************************
!  INDEX_MILLER
!  This subroutine finds the Miller indices [hkl]
!  from a string. The string can contain signed integers
!  attached together, or separated by underscores,
!  with or without brackets, e.g.:
!         1-10       [1-10]
!        1_-1_0     [1_-1_0]
!  The notation without underscore only allows to use
!  one-digit signed integers (like "1-10"), while the
!  notation allows to use greater integers,
!  like "[12_-15_14]".
!  The "planestring" is converted to a real array
!  containing the Miller indices, e.g. "1-10" will be
!  converted into (1.d0  -1.d0  0.d0).
!********************************************************
!
SUBROUTINE INDEX_MILLER(planestring,planeindices,ifail)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: planestring
CHARACTER(LEN=LEN_TRIM(planestring)):: temp, temp2
INTEGER:: i, m, mint
INTEGER,INTENT(OUT):: ifail !0=success; 1=error while reading string
INTEGER:: strpos
REAL(dp):: msign   !sign of value
REAL(dp),DIMENSION(3),INTENT(OUT):: planeindices
!
ifail=0
planeindices(:) = 0.d0
!
m = 1
mint = 0
msign = 1.d0
temp = planestring
!
!Detect if there are dots, commas or other forbidden characters
IF( SCAN(planestring,".,:/!")>0 ) THEN
  ifail = 2
  RETURN
ENDIF
!
!If there are brackets, remove them
strpos = SCAN(temp,'[')
IF( strpos>0 ) THEN
  temp = ADJUSTL(temp(strpos+1:))
ENDIF
strpos = SCAN(temp,']')
IF( strpos>0 ) THEN
  temp = ADJUSTL(temp(:strpos-1))
ENDIF
!
IF( SCAN(temp,'_').NE.0 ) THEN
  !values are separated by underscores, e.g. "1_-1_0"
  !read first value
  strpos = SCAN(temp,'_')
  READ(temp(1:strpos-1),*,ERR=100,END=100) planeindices(1)
  temp = temp(strpos+1:)
  !read second value
  strpos = SCAN(temp,'_')
  READ(temp(1:strpos-1),*,ERR=100,END=100) planeindices(2)
  temp = temp(strpos+1:)
  !read third value
  READ(temp,*,ERR=100,END=100) planeindices(3)
ELSE
  !values are given as attached integers, e.g. "1-10"
  !convert them to a vector like [1 -1 0]
  DO i=1, LEN_TRIM(temp)
    IF(m>3) GOTO 100
    READ(temp(i:i),*,ERR=100,END=100) temp2
    IF(temp2=='-') THEN
      msign = -1.d0
    ELSEIF(temp2=='+') THEN
      msign = 1.d0
    ELSEIF(LEN_TRIM(temp2)>0) THEN
      READ(temp2,*,ERR=100,END=100) mint
      planeindices(m) = msign*DBLE(mint)
      m = m+1
      msign = 1.d0
    ENDIF
  ENDDO
ENDIF
RETURN
!
100 CONTINUE
ifail=1
!
END SUBROUTINE INDEX_MILLER
!
!
!********************************************************
!  INDEX_MILLER_HCP
!  This subroutine finds the Miller indices [hkil]
!  from a string. The string can contain signed integers
!  attached together, or separated by underscores,
!  with or without brackets, e.g.:
!         1-100       [1-100]
!        1_-1_00     [1_-1_0_0]
!  The notation without underscore only allows to use
!  one-digit signed integers (like "1-10"), while the
!  notation allows to use greater integers,
!  like "[12_-15_3_14]".
!  The "planestring" is converted to a real array
!  containing only the 3 relevant Miller indices hkl,
!  e.g. "1-10" will be converted into (1.d0  -1.d0  0.d0).
!********************************************************
!
SUBROUTINE INDEX_MILLER_HCP(planestring,planeindices,ifail)
!
IMPLICIT NONE
CHARACTER(LEN=16),INTENT(IN):: planestring
CHARACTER(LEN=16):: temp, temp2
INTEGER:: i, m, mint
INTEGER,INTENT(OUT):: ifail !0=success; 1=error while reading string; 2=h+k not equal to -i
INTEGER:: strpos
REAL(dp):: msign   !sign of value
REAL(dp),DIMENSION(3),INTENT(OUT):: planeindices
!
ifail=0
planeindices(:) = 0.d0
!
m = 1
mint = 0
msign = 1.d0
temp = planestring
!
!Detect if there are dots, commas or other forbidden characters
IF( SCAN(planestring,".,:/!")>0 ) THEN
  ifail = 2
  RETURN
ENDIF
!
!If there are brackets, remove them
strpos = SCAN(temp,'[')
IF( strpos>0 ) THEN
  temp = ADJUSTL(temp(strpos+1:))
ENDIF
strpos = SCAN(temp,']')
IF( strpos>0 ) THEN
  temp = ADJUSTL(temp(:strpos-1))
ENDIF
!
IF( SCAN(temp,'_').NE.0 ) THEN
  !values are separated by underscores, e.g. "1_-1_0_0"
  !read first value
  strpos = SCAN(temp,'_')
  READ(temp(1:strpos-1),*,ERR=100,END=100) planeindices(1)
  temp = temp(strpos+1:)
  !read second value
  strpos = SCAN(temp,'_')
  READ(temp(1:strpos-1),*,ERR=100,END=100) planeindices(2)
  temp = temp(strpos+1:)
  !read third value (save it in mint)
  strpos = SCAN(temp,'_')
  READ(temp(1:strpos-1),*,ERR=100,END=100) mint
  temp = temp(strpos+1:)
  !read fourth value
  READ(temp,*,ERR=100,END=100) planeindices(3)
  !Check that -i=h+k
  IF( -1*mint==NINT(planeindices(1))+NINT(planeindices(2)) ) THEN
    !Notation is good, go on
    RETURN
  ELSE
    !h+k is not equal to -i => return with error
    ifail = 2
    RETURN
  ENDIF
ELSE
  !values are given as attached integers, e.g. "1-100"
  !convert them to a vector like [1 -1 0]
  strpos=0
  DO i=1, LEN_TRIM(temp)
    IF(strpos>4) GOTO 100
    READ(temp(i:i),*,ERR=100,END=100) temp2
    IF(temp2=='-') THEN
      msign = -1.d0
    ELSEIF(temp2=='+') THEN
      msign = 1.d0
    ELSEIF(LEN_TRIM(temp2)>0) THEN
      READ(temp2,*,ERR=100,END=100) mint
      planeindices(m) = msign*DBLE(mint)
      strpos=strpos+1
      IF( strpos==3 ) THEN
        !Check that -i=h+k
        IF( -1*NINT(planeindices(3)).NE.NINT(planeindices(1))+NINT(planeindices(2)) ) THEN
          ifail = 2
        ENDIF
      ENDIF
      m = MIN(3,m+1)
      msign = 1.d0
    ENDIF
  ENDDO
ENDIF
IF(strpos.NE.4) ifail=1
RETURN
!
100 CONTINUE
ifail=1
!
END SUBROUTINE INDEX_MILLER_HCP
!
!
!********************************************************
!  HKIL2UVW
!  This subroutine converts Miller-Bravais [hkil] notation
!  into Miller [uvw] notation.
!  This routine returns the minimum integer indices,
!  i.e. it attempts to divide u,v,w by their common divisor.
!  NOTE: index i is passed to this routine but not used,
!  therefore its value can be set to anything.
!********************************************************
!
SUBROUTINE HKIL2UVW(h,k,i,l,u,v,w)
!
REAL(dp):: a, b, c
REAL(dp),INTENT(IN):: h, k, i, l
REAL(dp),INTENT(OUT):: u, v, w
!
c = 1.d0
!
!Convert
u = 2.d0*h + k
v = h + 2.d0*k
w = l
!Check for common divisor
IF( DABS(u)>0.1d0 .AND. NINT(DABS(v))>0.1d0 ) THEN
  a = GCD( NINT(DABS(u)) , NINT(DABS(v)) )
ELSE
  a = MAX(DABS(u),DABS(v))
ENDIF
IF( DABS(u)>0.1d0 .AND. NINT(DABS(w))>0.1d0 ) THEN
  b = GCD( NINT(DABS(u)) , NINT(DABS(w)) )
ELSE
  b = MAX(DABS(u),DABS(w))
ENDIF
IF( DABS(a)>0.1d0 .AND. NINT(b)>0.1d0 ) THEN
  c = GCD( NINT(DABS(a)),NINT(DABS(b)) )
ELSE  !i.e. a==0 or b==0
  c = MAX( DABS(a) , DABS(b) )
ENDIF
IF( DABS(c)<0.1d0 ) c=1.d0  !avoid division by zero
!Update uvw
u = u / c
v = v / c
w = w / c
!
!
END SUBROUTINE HKIL2UVW
!
!
!********************************************************
!  UVW2HKIL
!  This subroutine converts Miller notation [uvw]
!  into Miller-Bravais notation [hkil].
!  This routine returns the minimum integer indices,
!  i.e. it attempts to divide u,v,w by their common divisor.
!********************************************************
!
SUBROUTINE UVW2HKIL(u,v,w,h,k,i,l)
!
REAL(dp):: a, b, c
REAL(dp),INTENT(IN):: u, v, w
REAL(dp),INTENT(OUT):: h, k, i, l
!
a = 1.d0
b = 1.d0
c = 1.d0
!
!Convert
h = (2.d0*u-v)/3.d0
k = (2.d0*v-u)/3.d0
i = -1.d0*(h+k)
l = l
!
!
END SUBROUTINE UVW2HKIL
!
!
!********************************************************
!  MILLER2VEC
!  This subroutine converts a string containing
!  Miller indices [hkl] or [hkil], into a
!  Cartesian vector (x y z). It includes the appropriate
!  accounting for crystal orientation. The final Cartesian
!  vector is normalized, i.e. it is a unit vector.
!********************************************************
!
SUBROUTINE MILLER2VEC(H,dir,ORIENT,vector,ifail)
!
IMPLICIT NONE
CHARACTER(LEN=16),INTENT(IN):: dir   !crystallographic direction, e.g. [110] or [2-110]
INTEGER:: i
INTEGER,INTENT(OUT):: ifail !0=success; 1=error while reading string; 2=h+k not equal to -i
REAL(dp):: u, v, w, x, z1, z2
REAL(dp):: V1, V2, V3
REAL(dp),DIMENSION(3),INTENT(OUT):: vector
REAL(dp),DIMENSION(3):: MILLER       !Miller indices
REAL(dp),DIMENSION(3,3),INTENT(IN):: H !box vectors
REAL(dp),DIMENSION(3,3),INTENT(IN):: ORIENT !current crystallographic orientation of the system
REAL(dp),DIMENSION(3,3):: ORIENTN !normalized ORIENT
!
ifail=0
vector(:) = 0.d0
MILLER(:) = 0.d0
ORIENTN(:,:) = 0.d0
!
!Try reading Miller indices [hkl]
CALL INDEX_MILLER(dir,MILLER,ifail)
IF(ifail==0) THEN
  !Reading Miller indices [hkl] was successful
  !
  !If the system has a defined crystallographic orientation ORIENT,
  !then Miller indices are defined in that basis
  !=> rotate Vplane(1,:) to express it in the basis of H(:,:)
  IF( ANY( NINT(ORIENT(:,:)).NE.0 ) ) THEN
    !Normalize orientation vectors
    DO i=1,3
      ORIENTN(i,:) = ORIENT(i,:) / VECLENGTH(ORIENT(i,:))
    ENDDO
    V1 = MILLER(1)
    V2 = MILLER(2)
    V3 = MILLER(3)
    MILLER(1) = ORIENTN(1,1)*V1 + ORIENTN(1,2)*V2 + ORIENTN(1,3)*V3
    MILLER(2) = ORIENTN(2,1)*V1 + ORIENTN(2,2)*V2 + ORIENTN(2,3)*V3
    MILLER(3) = ORIENTN(3,1)*V1 + ORIENTN(3,2)*V2 + ORIENTN(3,3)*V3
  ENDIF
  !
ELSE
  !Try to read [hkil] Miller indices
  CALL INDEX_MILLER_HCP(dir,MILLER,ifail)
  IF( ifail==0 ) THEN
    !Reading Miller indices [hkil] was successful
    !
    !If the system has a defined crystallographic orientation ORIENT,
    !then Miller indices are defined in that basis
    !=> apply rotation to express it in the basis of H(:,:)
    IF( ANY( NINT(ORIENT(:,:)).NE.0 ) ) THEN
      !Normalize orientation vectors
      DO i=1,3
        ORIENTN(i,:) = ORIENT(i,:) / VECLENGTH(ORIENT(i,:))
      ENDDO
      V1 = MILLER(1)
      V2 = MILLER(2)
      V3 = MILLER(3)
      MILLER(1) = ORIENTN(1,1)*V1 + ORIENTN(1,2)*V2 + ORIENTN(1,3)*V3
      MILLER(2) = ORIENTN(2,1)*V1 + ORIENTN(2,2)*V2 + ORIENTN(2,3)*V3
      MILLER(3) = ORIENTN(3,1)*V1 + ORIENTN(3,2)*V2 + ORIENTN(3,3)*V3
    ENDIF
    !
    !Convert [hkil] notation into [uvw]
    u = 2.d0*MILLER(1) + MILLER(2)
    v = MILLER(1) + 2.d0*MILLER(2)
    w = MILLER(3)
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
    !Save final indices into MILLER
    MILLER(1) = u / x
    MILLER(2) = v / x
    MILLER(3) = w / x
    !
  ENDIF
  !
ENDIF
!
IF( ifail==0 ) THEN
  !Use Miller indices to define Cartesian vector
  vector(:) = MILLER(1)*H(1,:) + MILLER(2)*H(2,:) + MILLER(3)*H(3,:)
  !Normalize vector
  vector(:) = vector(:) / VECLENGTH(vector(:))
ENDIF
!
END SUBROUTINE MILLER2VEC
!
!
!********************************************************
!  COMPFORMULA
!  This subroutine extracts the compound formula from the
!  given atomic site list P taking the auxiliary list AUX
!  into account if it contains site occupancies.
!  The formula is returned as string.
!  The compound summed atomic mass is returned as well.
!  (C) Juri Barthel 
!*     Gemeinschaftslabor fuer Elektronenmikroskopie
!*     RWTH Aachen (GERMANY)
!*     ju.barthel@fz-juelich.de
!********************************************************
!
SUBROUTINE COMPFORMULA(P,AUXNAMES,AUX,formula,mass)
!
USE atoms
!
IMPLICIT NONE
INTEGER:: i, j, iaux
INTEGER:: occ
INTEGER:: NP, Naux, NS, NA, NF
INTEGER,DIMENSION(ATOMMAXZ+10):: ispecies
REAL(dp):: PO, amass, socc
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
REAL(dp),INTENT(OUT):: mass
CHARACTER(LEN=2):: species
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128), INTENT(OUT):: formula
CHARACTER(LEN=128) :: temp, msg
!
! 
!Initialize variables
mass = 0.0d+0
amass = 0.0d+0
formula = ''
temp = ''
msg = ''
Naux = 0
occ = 0 ! no occupancy information by default
PO = 1.0d+0 ! 1.0 occupancy of each atomic site by default
ispecies = 0 ! preset species hash
NP = SIZE(P,1) ! number of sites
IF (ALLOCATED(AUX).AND.ALLOCATED(AUXNAMES)) THEN ! Get number of auxiliary properties
  Naux = SIZE(AUX,2) ! ? Redundant but maybe better MIN(SIZE(AUX,2),SIZE(AUXNAMES))
ENDIF
IF (Naux>0) THEN
  DO iaux=1, Naux
    IF("occ"==TRIM(AUXNAMES(iaux))) occ = iaux
  ENDDO
END IF
!
!Get different atomic species
!    from list P(:,4) = atomic numbers of all atomic sites
CALL FIND_NSP(P(:,4),aentries)
NS = SIZE(aentries,1) ! number of different species
!Get occupation numbers of the different species
IF (NS>0) THEN
  ispecies = 0
  !Remember the index of each species in aentries
  !     but only for those species which are on atom sites
  DO i=1, NS
    j = NINT(aentries(i,1))
    IF (j>0.AND.aentries(i,2)>0.0d+0) ispecies(j) = i
  ENDDO
  !Erase atom species counts in aentries
  aentries(1:NS,2) = 0.0d+0
  !Sum up the occupancies
  DO i=1, NP
    IF (P(i,4)<1.or.P(i,4)>ATOMMAXZ) CYCLE ! unknown element
    j = ispecies(NINT(P(i,4))) ! species index in aentries
    IF (j<=0) CYCLE ! invalid index
    IF (occ>0) PO = AUX(i,occ) ! set individual occupancy from AUX
    aentries(j,2) = aentries(j,2) + PO ! accumulation
  ENDDO
ENDIF
! Create formula string and sum up atomic masses
DO i=1, NS
  CALL ATOMSPECIES(aentries(i,1),species)
  temp=''
  socc = aentries(i,2)
  IF (socc<=0.d0) CYCLE ! this species is not in the compound
  NA = INT(socc)
  NF = INT((socc - NA)*1.0d+1)
  IF( NA>0 .or. NF>0 ) THEN
    WRITE(temp,'(i10)') NA
    IF( NF>0) THEN
      WRITE(temp,'(a,i1)') TRIM(ADJUSTL(temp))//".", NF
    END IF
    temp = TRIM(species)//TRIM(ADJUSTL(temp))
  ENDIF
  msg = TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(temp))
  CALL ATOMMASS(species,amass)
  mass = mass + amass*socc
ENDDO
formula = msg
!
IF (ALLOCATED(aentries)) DEALLOCATE(aentries)
!
END SUBROUTINE COMPFORMULA
!
!
!
END MODULE crystallography
