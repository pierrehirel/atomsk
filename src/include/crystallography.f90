MODULE crystallography
!
!**********************************************************************************
!*  SUBROUTINES                                                                   *
!**********************************************************************************
!* This module contains routines related to crystallography, like                 *
!* manipulating Miller indexes or elastic tensor.                                 *
!**********************************************************************************
!* (C) April 2022                                                                 *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 06 April 2022                                    *
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
!* ELAST2TENSOR        converts from Voigt notation to full elastic tensor        *
!* CHECK_CTENSOR       checks if an elastic tensor is symmetric                   *
!* COMPFORMULA         extracts a compound formula from atom site lists P and AUX *
!* ELASTINDEX          reduces indices (i,j) into index m for 9x9 matrices        *
!* ELAST2INDEX         convert index m into indices (i,j) for 9x9 matrices        *
!* ROTELAST            rotates a 9x9 matrix                                       *
!**********************************************************************************
!
!
USE comv
USE constants
USE functions
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
INTEGER:: ifail
INTEGER:: strpos
REAL(dp):: msign   !sign of value
REAL(dp),DIMENSION(3):: planeindices
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
IF( temp(1:1)=='[' ) THEN
  temp = ADJUSTL(temp(2:))
ENDIF
strpos = LEN_TRIM(temp)
IF( temp(strpos:strpos)==']' ) THEN
  temp = temp(:strpos-1)
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
! e.g. "1-10" will be converted into (1.d0  -1.d0  0.d0).
!********************************************************
!
SUBROUTINE INDEX_MILLER_HCP(planestring,planeindices,ifail)
!
IMPLICIT NONE
CHARACTER(LEN=16),INTENT(IN):: planestring
CHARACTER(LEN=16):: temp, temp2
INTEGER:: i, m, mint
INTEGER:: ifail !0=success; 1=error while reading string; 2=h+k not equal to -i
INTEGER:: strpos
REAL(dp):: msign   !sign of value
REAL(dp),DIMENSION(3):: planeindices
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
IF( temp(1:1)=='[' ) THEN
  temp = ADJUSTL(temp(2:))
ENDIF
strpos = LEN_TRIM(temp)
IF( temp(strpos:strpos)==']' ) THEN
  temp = temp(:strpos-1)
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
  !read third value
  strpos = SCAN(temp,'_')
  READ(temp(1:strpos-1),*,ERR=100,END=100) planeindices(3)
  !Check that -i=h+k
  IF( -1*NINT(planeindices(3))==NINT(planeindices(1))+NINT(planeindices(2)) ) THEN
    !Notation is good, go on
    temp = temp(strpos+1:)
    !read fourth value
    READ(temp,*,ERR=100,END=100) planeindices(3)
    strpos=4
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
          RETURN
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
!  ELAST2TENSOR
!  This subroutine sets up the 9x9 elastic tensor,
!  provided the 9 (Voigt notation) elastic constants.
!********************************************************
!
SUBROUTINE ELAST2TENSOR(elcst,eltens)
!
IMPLICIT NONE
REAL(dp),DIMENSION(9),INTENT(IN):: elcst !C11,C22,C33,C23,C31,C12,C44,C55,C66
REAL(dp),DIMENSION(:,:),INTENT(INOUT):: eltens
!
eltens(:,:) = 0.d0
!
eltens(1,1) = elcst(1)
eltens(1,2) = elcst(6)
eltens(1,3) = elcst(5)
eltens(2,1) = elcst(6)
eltens(2,2) = elcst(2)
eltens(2,3) = elcst(4)
eltens(3,1) = elcst(5)
eltens(3,2) = elcst(4)
eltens(3,3) = elcst(3)
eltens(4,4) = elcst(7)
eltens(5,5) = elcst(8)
 eltens(6,6) = elcst(9)
IF( SIZE(eltens,1)==9 .AND. SIZE(eltens,2)==9 ) THEN
  !Fill the rest by symmetry
  eltens(1:3,7:9) = eltens(1:3,4:6)
  eltens(4:6,7:9) = eltens(4:6,4:6)
  eltens(7:9,1:3) = eltens(4:6,1:3)
  eltens(7:9,4:6) = eltens(4:6,4:6)
  eltens(7:9,7:9) = eltens(4:6,4:6)
ENDIF

!
END SUBROUTINE ELAST2TENSOR
!
!
!********************************************************
!  CHECK_CTENSOR
!  This subroutine checks if an elastic tensor is
!  symmetric.
!********************************************************
!
SUBROUTINE CHECK_CTENSOR(C_tensor,status)
!
IMPLICIT NONE
INTEGER:: i, j
INTEGER:: status  !=0 if tensor is OK
                  !=1 if tensor is not symmetric, i.e. C(i,j).NE.C(j,i)
REAL(dp),DIMENSION(9,9),INTENT(IN):: C_tensor
!
status = 0
DO i=2,9
  DO j=1,i-1
    IF( DABS(C_tensor(i,j)-C_tensor(j,i))>1.d-3 ) THEN
      status=1
    ENDIF
  ENDDO
ENDDO
!
END SUBROUTINE CHECK_CTENSOR
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
!********************************************************
!  ELASTINDEX
!  This function converts the indices i and j into
!  single index m.
!********************************************************
FUNCTION ELASTINDEX(i,j) RESULT(m)
!
IMPLICIT NONE
INTEGER, INTENT(IN):: i, j
INTEGER:: m
!
m = 0
IF(i==j) THEN
  m = i
!
ELSEIF(i==1) THEN
  IF(j==2) THEN
    m = 6
  ELSEIF(j==3) THEN
    m = 8
  ENDIF
!
ELSEIF(i==2) THEN
  IF(j==1) THEN
    m = 9
  ELSEIF(j==3) THEN
    m = 4
  ENDIF
!
ELSEIF(i==3) THEN
  IF(j==1) THEN
    m = 5
  ELSEIF(j==2) THEN
    m = 7
  ENDIF
ENDIF
!
END FUNCTION ELASTINDEX
!
!
!********************************************************
!  ELAST2INDEX
!  This function converts the index i into indices m and n.
!  This is for use by the function ROTELAST below.
!********************************************************
FUNCTION ELAST2INDEX(i) RESULT(mn)
!
IMPLICIT NONE
INTEGER,INTENT(IN):: i
INTEGER,DIMENSION(2):: mn !mn(1)=m and mn(2)=n
!
mn(:) = 0
!
SELECT CASE(i)
!
CASE(1,2,3)
  mn(1) = i
  mn(2) = i
CASE(4)
  mn(1) = 2
  mn(2) = 3
CASE(5)
  mn(1) = 3
  mn(2) = 1
CASE(6)
  mn(1) = 1
  mn(2) = 2
CASE(7)
  mn(1) = 3
  mn(2) = 2
CASE(8)
  mn(1) = 1
  mn(2) = 3
CASE(9)
  mn(1) = 2
  mn(2) = 1
CASE DEFAULT
  PRINT*, "Problem with elastindex"
END SELECT
!
RETURN
!
END FUNCTION ELAST2INDEX
!
!
!********************************************************
!  ROTELAST
!  This function applies a rotation defined by a
!  3x3 matrix to a 9x9 tensor. Although this operation
!  is general, it is intended to rotate the elastic tensor.
!********************************************************
FUNCTION ROTELAST(ELTENS,T) RESULT(newELTENS)
!
IMPLICIT NONE
INTEGER:: i, j, k, l, m, n
INTEGER,DIMENSION(2):: mn
REAL(dp), DIMENSION(9,9):: ELTENS    !Elastic tensor before
REAL(dp), DIMENSION(9,9):: newELTENS !Elastic tensor after
REAL(dp), DIMENSION(3,3), INTENT(IN):: T  !Rotation matrix (3x3)
REAL(dp), DIMENSION(9,9):: Q  !Rotation matrix (9x9)
!
!Initialize variables
mn(:) = 0
newELTENS(:,:) = 0.d0
Q(:,:) = 0.d0
!
!Build the 9x9 rotation matrix ROT9
DO i=1,9
  mn = ELAST2INDEX(i)
  m = mn(1)
  n = mn(2)
  DO j=1,9
    mn = ELAST2INDEX(j)
    k = mn(1)
    l = mn(2)
    Q(i,j) = T(k,m)*T(l,n)
  ENDDO
ENDDO
!
!Rotate the elastic tensor:  C' = tQ * C * Q
newELTENS = MATMUL( TRANSPOSE(Q) , MATMUL(ELTENS,Q) )
!
RETURN
!
END FUNCTION ROTELAST
!
!
!
END MODULE crystallography
