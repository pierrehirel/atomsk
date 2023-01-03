MODULE mode_density
!
!**********************************************************************************
!*  MODE_DENSITY                                                                  *
!**********************************************************************************
!* This module reads an input file containing atom positions and computes         *
!* the "density" of a given property. The property must have a defined value for  *
!* each atom. Instead of having zero-dimensional atoms with a given property,     *
!* atoms are replaced by Gaussian functions so as to smear the property.          *
!* The "density" at any point in space is defined as the sum of all Gaussian      *
!* functions. Depending on the value of <den_type>, the density will be           *
!* computed along a given <axis> (if den_type=1), or in the plane normal          *
!* to the <axis> (if den_type=2), or the in the whole volume.                     *
!**********************************************************************************
!* (C) February 2015 - Pierre Hirel                                               *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 13 Dec. 2022                                     *
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
USE readin
USE options
USE writeout
!
!
CONTAINS
!
!
SUBROUTINE DENSITY_XYZ(inputfile,options_array,property,den_type,axis,Sigma,prefix)
!
!
IMPLICIT NONE
CHARACTER(LEN=1),INTENT(IN):: axis  !density will be computed along this axis (if den_type=1)
                                    !or in the plane normal to this axis (if den_type=2)
                                    !or in the volume (if den_type=3) (default)
CHARACTER(LEN=128),INTENT(IN):: property !name of the property whose density will be calculated
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=*),INTENT(IN):: prefix
INTEGER,INTENT(IN):: den_type       !type of density to output: 1, 2 or 3 (for 1-D, 2-D or 3-D)
REAL(dp),INTENT(IN):: Sigma         !square root of variance
!
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: msg
CHARACTER(LEN=4096):: outputfile
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of output file formats
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
LOGICAL:: peak, dip  !is it a peak/dip in the density?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: a1, a2, a3
INTEGER:: progress      !To show calculation progress
INTEGER:: i, j, k, l, m, n, o
INTEGER:: jmin, jmax, kmin, kmax, lmin, lmax
INTEGER:: Ninter, Nvac  !Number of detected interstitials, vacancies
INTEGER:: Nx, Ny, Nz    !number of points in the grid
REAL(dp):: A
REAL(dp):: cutoff       !cut-off for Gaussian functions (to speed up calculation)
REAL(dp):: distance
REAL(dp):: dx, dy,dz    !line/grid step along X, Y
REAL(dp):: prefactor
REAL(dp):: snumber      !atomic number of an atom
REAL(dp):: tempreal
REAL(dp):: x, y, z      !coordinates along the <axis> / in the grid
REAL(dp),DIMENSION(3,3):: Huc   !Box vectors of unit cell (unknown, set to 0 here)
REAL(dp),DIMENSION(3,3):: H       !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystallographic orientation of the system
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:),ALLOCATABLE:: DenGrid1     !density grid (1-D)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: DenGrid2   !density grid (2-D)
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE:: DenGrid3 !density grid (3-D)
REAL(dp),DIMENSION(:),ALLOCATABLE,TARGET:: DUMMY_PROP  !dummy property used for element density
REAL(dp),DIMENSION(:,:),ALLOCATABLE,TARGET:: AUX      !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE,TARGET:: P, S     !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Pinter, Pintertemp !positions of interstitials
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Pvac, Pvactemp     !positions of vacancies
REAL(dp),DIMENSION(:),POINTER:: PropPoint !pointer to the property whose density must be computed
!
!
!Initialize variables
a1=1
a2=2
a3=3
 cutoff=2.5d0*Sigma  !cut-off for Gaussian functions (to speed up calculation)
                     !rationalization: FWHM of Gaussian is about 2.36*Sigma
                     !therefore a radius of 2.5*Sigma should do it
dx = 0.2d0
dy = 0.2d0
dz = 0.2d0
Huc(:,:) = 0.d0
 C_tensor(:,:) = 0.d0
snumber = 0.d0
IF(ALLOCATED(DenGrid1)) DEALLOCATE(DenGrid1)
IF(ALLOCATED(DenGrid2)) DEALLOCATE(DenGrid2)
IF(ALLOCATED(DenGrid3)) DEALLOCATE(DenGrid3)
!
!
CALL ATOMSK_MSG(4066,(/property/),(/DBLE(den_type),Sigma/))
!
!
100 CONTINUE
! Read the input file
CALL READ_AFF(inputfile,H,P,S,comment,AUXNAMES,AUX)
!
! Apply options if any
CALL OPTIONS_AFF(options_array,Huc,H,P,S,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
IF(nerr>0) GOTO 1000
!
!
!
200 CONTINUE
!Define properties for the simulation
! Define axes
IF( den_type==1 ) THEN
  !density will be computed along the direction a1
  SELECT CASE(axis)
  CASE("x","X")
    a1=1
    a2=2
    a3=3
  CASE("y","Y")
    a1=2
    a2=3
    a3=1
  CASE("z","Z")
    a1=3
    a2=1
    a3=2
  END SELECT
ELSEIF( den_type==2 ) THEN
  !density will be computed in the plane (a1,a2) normal to a3
  SELECT CASE(axis)
  CASE("x","X")
    a1=2
    a2=3
    a3=1
  CASE("y","Y")
    a1=3
    a2=1
    a3=2
  CASE("z","Z")
    a1=1
    a2=2
    a3=3
  END SELECT
ELSE
  a1=1
  a2=2
  a3=3
ENDIF
!
! Define the discrete grid
! If Nx, Ny or Nz is too large (i.e. if the box H(:,:) is very wide)
! then set Nx, Ny and/or Nz to a smaller value, and increase the step dx, dy and/or dz
IF( den_type==1 ) THEN
  ! By default use a step dx=0.1 A
  dx = 0.1d0
  !Compute Nx
  Nx = NINT( ( MAXVAL(P(:,a1))-MINVAL(P(:,a1)) ) / dx )
  IF( Nx<200 ) THEN
    Nx = 200
  ELSEIF( Nx>10000 ) THEN
    Nx = 10000
  ENDIF
  dx = H(a1,a1) / DBLE(Nx)
  !
ELSEIF( den_type==2 ) THEN
  ! By default use a step dx=dy=0.1 A
  dx = 0.1d0
  dy = 0.1d0
  !Compute Nx and Ny
  Nx = NINT( ( MAXVAL(P(:,a1))-MINVAL(P(:,a1)) ) / dx )
  Ny = NINT( ( MAXVAL(P(:,a2))-MINVAL(P(:,a2)) ) / dy )
  IF( Nx<20 ) THEN
    Nx = 20
  ELSEIF( Nx>100 ) THEN
    Nx = 100
  ENDIF
  IF( Ny<20 ) THEN
    Ny = 20
  ELSEIF( Ny>100 ) THEN
    Ny = 100
  ENDIF
  dx = H(a1,a1) / DBLE(Nx)
  dy = H(a2,a2) / DBLE(Ny)
  !
ELSE
  ! By default use a step dx=dy=2.0 A
  dx = 2.d0
  dy = 2.d0
  dz = 2.d0
  !Compute Nx, Ny, Nz
  Nx = NINT( ( MAXVAL(P(:,a1))-MINVAL(P(:,a1)) ) / dx )
  Ny = NINT( ( MAXVAL(P(:,a2))-MINVAL(P(:,a2)) ) / dy )
  Nz = NINT( ( MAXVAL(P(:,a3))-MINVAL(P(:,a3)) ) / dz )
  IF( Nx<40 .OR. Nx>100 ) THEN
    Nx = 100
    dx = H(a1,a1) / DBLE(Nx)
  ENDIF
  IF( Ny<40 .OR. Ny>100 ) THEN
    Ny = 100
    dy = H(a2,a2) / DBLE(Ny)
  ENDIF
  IF( Nz<40 .OR. Nz>100 ) THEN
    Nz = 100
    dz = H(a3,a3) / DBLE(Nz)
  ENDIF
ENDIF
WRITE(msg,'(a13,3(i6,1X))') 'Nx, Ny, Nz = ', Nx, Ny, Nz
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
! Set PropPoint to the property whose density must be computed
NULLIFY(PropPoint)
IF( property=="mass" ) THEN
  !Check if masses are defined in array AUX
  k=0
  IF( ALLOCATED(AUXNAMES) ) THEN
    DO i=1,SIZE(AUXNAMES)
      IF( AUXNAMES(i)=="mass" ) THEN
        k=i
      ENDIF
    ENDDO
  ENDIF
  !
  IF( k>0 ) THEN
    !Masses are defined in column #k of array AUX
    PropPoint => AUX(:,k)
  ELSE
    !Masses are not defined in AUX
    !Fetch standard (NIST) atomic masses and store them inside array DUMMY_PROP
    ALLOCATE(DUMMY_PROP(SIZE(P,1)))
    DUMMY_PROP(:) = 0.d0
    DO i=1,SIZE(P,1)
      CALL ATOMSPECIES(P(i,4),species)
      CALL ATOMMASS(species,DUMMY_PROP(i))
    ENDDO
    PropPoint => DUMMY_PROP(:)
  ENDIF
ELSE
  !Check if "property" is an atom species
  IF( LEN_TRIM(property) <= 2 ) THEN
    species = ADJUSTL(property)
    CALL ATOMNUMBER(species,snumber)
    IF( snumber.NE.0.d0 ) THEN
      !It is indeed an atomic number => compute the density of these atoms only
      !For that we create a dummy array that will act as a mask:
      !it will have a value of 1 for atoms that match, 0 otherwise
      ALLOCATE(DUMMY_PROP(SIZE(P,1)))
      DUMMY_PROP(:) = 0.d0
      DO i=1,SIZE(P,1)
        IF( NINT(P(i,4)) == NINT(snumber) ) THEN
          DUMMY_PROP(:) = 1.d0
        ENDIF
      ENDDO
      PropPoint => DUMMY_PROP(:)
    ELSE
      !Look for the "property" in the auxiliary properties
      IF( ALLOCATED(AUXNAMES) ) THEN
        DO i=1,SIZE(AUXNAMES)
          IF( AUXNAMES(i) == property ) THEN
            PropPoint => AUX(:,i)
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDIF
  !
  !At this point, if pointer is not associated with anything then we are in trouble
  IF( snumber<=0 .AND. .NOT. ASSOCIATED(PropPoint) ) THEN
    CALL ATOMSK_MSG(2817,(/TRIM(property)/),(/0.d0/))
    nerr = nerr+1
    GOTO 1000
  ENDIF
ENDIF
!
!
!
300 CONTINUE
! Compute density
CALL ATOMSK_MSG(4067,(/property/),(/DBLE(den_type),Sigma/))
!
!
IF( den_type==1 ) THEN    !!!!!!!   1-D DENSITY   !!!!!!!
  !density profile will be computed along the direction a1
  ALLOCATE(DenGrid1(Nx))
  DenGrid1(:) = 0.d0
  A = 1.d0 / DSQRT(2.d0*pi*Sigma**2)
  !
  progress=0  !to count progress
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,m,x) &
  !$OMP& REDUCTION(+:DenGrid1)
  DO i=1,SIZE(P,1)
    progress = progress+1
    IF( SIZE(P,1)>100000 ) THEN
      !If there are many atoms, display a fancy progress bar
      CALL ATOMSK_MSG(10,(/""/),(/DBLE(progress),DBLE(SIZE(P,1))/))
    ENDIF
    !
    DO j=1,Nx
      x = DBLE(j)*dx
      !
      prefactor = PropPoint(i) * A
      !
      IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
        distance = x - P(i,a1)
        DenGrid1(j) = DenGrid1(j) + prefactor * DEXP( -(distance**2) / (2.d0*Sigma**2) )
        !
        !Add the contribution of periodic images
        DO m=-1,1
          IF( m.NE.0 ) THEN
            distance = DABS( x - (P(i,a1)+DBLE(m)*H(a1,a1)) )
            DenGrid1(j) = DenGrid1(j) + prefactor * DEXP( -(distance**2) / (2.d0*Sigma**2) )
          ENDIF
        ENDDO
        !
      ENDIF
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  !
  !
ELSEIF( den_type==2 ) THEN   !!!!!!!   2-D DENSITY   !!!!!!!
  !density will be computed in the plane (a1,a2) normal to a3
  ALLOCATE(DenGrid2(Nx,Ny))
  DenGrid2(:,:) = 0.d0
  A = 1.d0 / DSQRT(2.d0*pi*Sigma**2)
  !
  !
  progress=0  !to count progress
  tempreal=0.d0
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,m,n,x,y) &
  !$OMP& REDUCTION(+:tempreal,DenGrid2)
  DO i=1,SIZE(P,1)
    progress = progress+1
    IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
      !
      IF( SIZE(P,1)>10000 ) THEN
        !If there are many atoms, display a fancy progress bar
        CALL ATOMSK_MSG(10,(/""/),(/DBLE(progress),DBLE(SIZE(P,1))/))
      ENDIF
      !
      prefactor = PropPoint(i) * A
      DO j=1,Nx
        x = DBLE(j)*dx
        DO k=1,Ny
          tempreal = 0.d0
          y = DBLE(k)*dy
          !
          distance = DSQRT( (x-P(i,a1))**2 + (y-P(i,a2))**2 )
          tempreal = tempreal + prefactor*DEXP( -distance**2 / (2.d0*Sigma**2) )
          !
          !Add the contribution of its periodic images
          DO m=-1,1
            DO n=-1,1
              IF( m.NE.0 .OR. n.NE.0 ) THEN
                distance = DSQRT( (x-(P(i,a1)+DBLE(m)*H(a1,a1)))**2 + (y-(P(i,a2)+DBLE(n)*H(a2,a2)))**2 )
                tempreal = tempreal + prefactor*DEXP( -distance**2 / (2.d0*Sigma**2) )
              ENDIF
            ENDDO
          ENDDO
          !
          DenGrid2(j,k) = DenGrid2(j,k) + tempreal
          !
        ENDDO !loop on k
      ENDDO !loop on j
      !
    ENDIF
  ENDDO !loop on i
  !$OMP END PARALLEL DO
  !
  !Integrate values in given area
  z = 0.d0
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,k) &
  !$OMP& REDUCTION(+:z)
  DO j=1,Nx
    x = DBLE(j)*dx
    DO k=1,Ny
      y = DBLE(k)*dy
      z = z + DenGrid2(j,k)*dx*dy
    ENDDO
  ENDDO
  !$OMP END PARALLEL DO
  !
  !
ELSE                       !!!!!!!   3-D DENSITY   !!!!!!!
  !density will be computed in the volume
  ALLOCATE(DenGrid3(Nx,Ny,Nz))
  DenGrid3(:,:,:) = 0.d0
  A = 1.d0 / ( (DSQRT(2.d0*pi*Sigma**2))**3 )
  !
  progress=0  !to count progress
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,l,m,n,o,x,y,z) &
  !$OMP& REDUCTION(+:DenGrid3)
  DO i=1,SIZE(P,1)
    progress = progress+1
    IF( .NOT.ALLOCATED(SELECT) .OR. SELECT(i) ) THEN
      IF( SIZE(P,1)>1000 ) THEN
        !If there are many atoms, display a fancy progress bar
        CALL ATOMSK_MSG(10,(/""/),(/DBLE(progress),DBLE(SIZE(P,1))/))
      ENDIF
      !
      prefactor = PropPoint(i) * A
      !
      !Define jmin and jmax, considering that the Gaussian function vanishes for a radius > cutoff
      jmin=MAX( 1  , FLOOR(   (P(i,a1)-cutoff)/dx ) )
      jmax=MIN( Nx , CEILING( (P(i,a1)+cutoff)/dx ) )
      kmin=MAX( 1  , FLOOR(   (P(i,a2)-cutoff)/dy ) )
      kmax=MIN( Ny , CEILING( (P(i,a2)+cutoff)/dy ) )
      lmin=MAX( 1  , FLOOR(   (P(i,a3)-cutoff)/dz ) )
      lmax=MIN( Nz , CEILING( (P(i,a3)+cutoff)/dz ) )
      !
      DO j=1,Nx
        x = DBLE(j)*dx
        DO k=1,Ny
          y = DBLE(k)*dy
          DO l=1,Nz
            z = DBLE(l)*dz
            distance = DSQRT( (x-P(i,a1))**2 + (y-P(i,a2))**2 + (z-P(i,a3))**2 )
            DenGrid3(j,k,l) = DenGrid3(j,k,l) + prefactor * DEXP( -distance**2 / (2.d0*Sigma**2) )
            !
            !Add the contribution of periodic images
            DO m=-1,1
              DO n=-1,1
                DO o=-1,1
                  IF( m.NE.0 .OR. n.NE.0 .OR. o.NE.0 ) THEN
                    distance = DSQRT( (x-(P(i,a1)+DBLE(m)*H(a1,a1)))**2 + &
                             &        (y-(P(i,a2)+DBLE(n)*H(a2,a2)))**2 + &
                             &        (z-(P(i,a3)+DBLE(o)*H(a3,a3)))**2   )
                    DenGrid3(j,k,l) = DenGrid3(j,k,l) + prefactor*DEXP( -distance**2 / (2.d0*Sigma**2) )
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
            !
          ENDDO
        ENDDO
      ENDDO
      !
    ENDIF
  ENDDO !loop on i
  !$OMP END PARALLEL DO
  !
  !
  WRITE(msg,*) 'Detecting peaks and dips in the density...'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ! Now we know the value of density in the whole volume
  ! => compute average density
  A = SUM( DenGrid3(:,:,:) ) / DBLE(Nx*Ny*Nz)
  !
  ALLOCATE( Pintertemp(SIZE(P,1),4) )
  Pintertemp(:,:) = 0.d0
  ALLOCATE( Pvactemp(SIZE(P,1),4) )
  Pvactemp(:,:) = 0.d0
  ! Detect positions (x,y,z) where the density is significantly higher than average
  Ninter=0
  Nvac=0
  DO j=1,Nx
    x = DBLE(j)*dx
    DO k=1,Ny
      y = DBLE(k)*dy
      DO l=1,Nz
        z = DBLE(l)*dz
        !
        peak = .TRUE.
        dip  = .TRUE.
        DO m=MAX(1-j,-1),MIN(Nx-j,1)
          DO n=MAX(1-k,-1),MIN(Ny-k,1)
            DO o=MAX(1-l,-1),MIN(Nz-l,1)
              IF( DenGrid3(j,k,l) < DenGrid3(j+m,k+n,l+o) ) THEN
                peak = .FALSE.
              ENDIF
              IF( DenGrid3(j,k,l) > DenGrid3(j+m,k+n,l+o) ) THEN
                dip = .FALSE.
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        !
        IF( peak ) THEN
          !We detected a peak in the density
          WRITE(msg,'(a9,3f9.3)') 'Peak at: ', x, y, z
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          Ninter = Ninter+1
          IF( Ninter>SIZE(Pintertemp,1) ) PRINT*, "oversize Ninter = ", Ninter
          Pintertemp(Ninter,1) = x
          Pintertemp(Ninter,2) = y
          Pintertemp(Ninter,3) = z
          Pintertemp(Ninter,4) = 1.d0
        ENDIF
        IF( dip ) THEN
          !We detected a dip in the density
          WRITE(msg,'(a9,3f9.3)') 'Dip at:  ', x, y, z
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          Nvac = Nvac+1
          IF( Nvac>SIZE(Pvactemp,1) ) PRINT*, "oversize Nvac = ", Nvac
          Pvactemp(Nvac,1) = x
          Pvactemp(Nvac,2) = y
          Pvactemp(Nvac,3) = z
          Pvactemp(Nvac,4) = 1.d0
        ENDIF
        !
      ENDDO
    ENDDO
  ENDDO
  !
  ! Store final positions of intersitials and vacancies
  IF( Ninter > 0 ) THEN
    ALLOCATE( Pinter(Ninter,4) )
    Pinter(:,:) = Pintertemp(1:Ninter,:)
    DEALLOCATE(Pintertemp)
  ENDIF
  IF( Nvac > 0 ) THEN
    ALLOCATE( Pvac(Nvac,4) )
    Pvac(:,:) = Pvactemp(1:Nvac,:)
    DEALLOCATE(Pvactemp)
  ENDIF
  !
ENDIF
!
!
!
400 CONTINUE
CALL ATOMSK_MSG(4068,(/''/),(/0.d0/))
!
! Output results to data file
outputfile = TRIM(ADJUSTL(prefix))//"_"//TRIM(ADJUSTL(property))//"_density.dat"
!
IF(.NOT.overw) CALL CHECKFILE(outputfile,'writ')
OPEN(UNIT=30,FILE=outputfile,STATUS="UNKNOWN",FORM="FORMATTED")
!
IF( den_type==1 ) THEN
  WRITE(30,*) "# Linear density of "//TRIM(ADJUSTL(property))//" along "//axis//" computed with Atomsk"
  WRITE(30,*) "# Can be visualized with Gnuplot"
  WRITE(30,*) "# Step: ", dx
  WRITE(30,*) "#        x            density"
  DO i=1,Nx
    x = DBLE(i)*dx
    WRITE(30,'(2f16.8,1X)') x, DenGrid1(i)
  ENDDO
  !
ELSEIF( den_type==2 ) THEN
  WRITE(30,*) "# 2-D density plot of "//TRIM(ADJUSTL(property))//" computed with Atomsk"
  WRITE(30,*) "# Can be visualized with Gnuplot"
  WRITE(30,*) "# Step: ", dx, dy
  WRITE(30,*) "#        x               y            density"
  DO j=1,Nx
    x = DBLE(j)*dx
    DO k=1,Ny
      y = DBLE(k)*dy
      WRITE(30,'(3f16.8,1X)') x, y, DenGrid2(j,k)
    ENDDO
    WRITE(30,*) ""
  ENDDO
!
ELSE
  WRITE(30,*) "# 3-D density plot of "//TRIM(ADJUSTL(property))//" computed with Atomsk"
  WRITE(30,*) "# Step: ", dx, dy, dz
  WRITE(30,*) "#        x               y             z             density"
  DO j=1,Nx
    x = DBLE(j)*dx
    DO k=1,Ny
      y = DBLE(k)*dy
      DO l=1,Nz
        z = DBLE(l)*dz
        WRITE(30,'(4f16.8,1X)') x, y, z, DenGrid3(j,k,l)
      ENDDO
    ENDDO
    WRITE(30,*) ""
  ENDDO
ENDIF
CLOSE(30)
CALL ATOMSK_MSG(4039,(/TRIM(outputfile)/),(/0.d0/))
!
!
!Clear atomic data
IF(ALLOCATED(P)) DEALLOCATE(P)
IF(ALLOCATED(S)) DEALLOCATE(S)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
!
!
IF( den_type==3 ) THEN
  ! Write positions of interstitials / vacancies if any were detected
  ALLOCATE(outfileformats(1))
  outfileformats(1) = "xsf"
  IF( Ninter > 0 ) THEN
    outputfile = TRIM(ADJUSTL(prefix))//"_"//TRIM(ADJUSTL(property))//"_interstitials"
    CALL WRITE_AFF(outputfile,outfileformats,H,Pinter,S,comment,AUXNAMES,AUX)
  ENDIF
  IF( Nvac > 0 ) THEN
    outputfile = TRIM(ADJUSTL(prefix))//"_"//TRIM(ADJUSTL(property))//"_vacancies"
    CALL WRITE_AFF(outputfile,outfileformats,H,Pvac,S,comment,AUXNAMES,AUX)
  ENDIF
ENDIF
!
!
GOTO 1000
!
!
!
850 CONTINUE
nerr = nerr+1
GOTO 1000
!
!
!
1000 CONTINUE
IF(ALLOCATED(DUMMY_PROP)) DEALLOCATE(DUMMY_PROP)
IF(ALLOCATED(DenGrid1)) DEALLOCATE(DenGrid1)
IF(ALLOCATED(DenGrid2)) DEALLOCATE(DenGrid2)
IF(ALLOCATED(DenGrid3)) DEALLOCATE(DenGrid3)
!
!
END SUBROUTINE DENSITY_XYZ
!
!
END MODULE mode_density
