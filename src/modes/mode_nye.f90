MODULE mode_nye
!
!**********************************************************************************
!*  MODE_NYE                                                                      *
!**********************************************************************************
!* This module reads two sets of atomic coordinates and                           *
!* computes the corresponding Nye tensor.                                         *
!* The algorithm follows the method described in:                                 *
!*    C.S. Hartley, Y. Mishin, Acta Mater. 53 (2005) 1313                         *
!* Equation numbers also refer to this reference.                                 *
!* The names "system1" and "system2" refer to the reference system                *
!* (bulk with no defect) and the system to analyze, respectively.                 *
!* The reference system 1 can be either:                                          *
!* - a full reference system containing the same number of atoms as system2,      *
!*   and indexed in the exact same order;
!* - a unit cell of the material in system 2;                                     *
!* - NULL, in which case the "reference" environments are built on-the-fly        *
!*   by averaging atomic environments found in the system2.                       *
!**********************************************************************************
!* (C) October 2013 - Philippe Carrez                                             *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 28 Feb. 2025                                     *
!**********************************************************************************
!* OUTLINE:                                                                       *
!* 100        Read atom positions systems 1 and 2, construct neighbor lists       *
!* 200        Compute lattice correspondence tensor G for each atom               *
!* 300        Using G, compute Nye tensor for each atom                           *
!* 400        Output final results to file(s)                                     *
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
!General modules and routines
USE comv
USE subroutines
USE math
USE messages
USE neighbors
USE avgenv
!Modules for computing things
USE cmpt_G
!Module for reading input files
USE readin
USE writeout
!Module for applying options
USE options
!
!
!
CONTAINS
!
!
SUBROUTINE NYE_TENSOR(filefirst,filesecond,options_array,prefix,outfileformats)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: filefirst, filesecond, prefix
CHARACTER(LEN=2):: species
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE,INTENT(IN):: outfileformats !list of output file formats
CHARACTER(LEN=4096):: conffile
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
LOGICAL:: fileexists !does the file exist?
LOGICAL:: firstref   !is the first file used as reference system? (default: yes)
LOGICAL:: ucref      !is the reference a unit cell? (default: no)
LOGICAL:: usercutoff !did the user provide a fixed cutoff?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: i, j, k, m, n, iat
INTEGER:: nb_neigh, eps
INTEGER:: Nneighbors
INTEGER:: NNmin   !minimum number of neighbours to keep
INTEGER:: Nsites  !number of atomic environments
INTEGER:: ok
INTEGER:: progress      !To show calculation progress
INTEGER:: Q_matrix_rank, INFO, LWORK
INTEGER,DIMENSION(:),ALLOCATABLE:: Nlist !index of neighbours
INTEGER,DIMENSION(:),ALLOCATABLE:: Tab_PQ    !correspondance table between neighbor lists
INTEGER,DIMENSION(:),ALLOCATABLE:: newindex  !list of index after sorting
INTEGER,DIMENSION(:),ALLOCATABLE:: siteindex !for each atom, index of its type of site in Pref
INTEGER,DIMENSION(:),ALLOCATABLE:: sitescore !"score" of reference sites (only when filefirst is a unit cell)
INTEGER,DIMENSION(:,:),ALLOCATABLE:: NeighList1, NeighList2  !list of neighbors for systems 1 and 2
REAL(dp):: alpha, alpha_tmp !angles between two vectors
REAL(dp):: theta_max  !if angle between vectors in Pneigh and Qneigh exceed this value, exclude them
REAL(dp):: NeighFactor !%of tolerance in the cutoff for neighbor search
REAL(dp),PARAMETER:: radius=6.d0 !R for neighbour search (6A should be enough to find neighbours in any structure)
REAL(dp):: cutoff      !cutoff for neighbor search
REAL(dp):: tempreal
!Parameters used for the calculation: cutoff, NNmin, NeighFactor, theta_max
REAL(dp),DIMENSION(4):: params
REAL(dp),DIMENSION(3,3):: alpha_tensor, test_matrix
REAL(dp),DIMENSION(3,3):: Huc   !Box vectors of unit cell (unknown, set to 0 here)
REAL(dp),DIMENSION(3,3):: H1, H2 !Box vectors of systems 1 and 2
REAL(dp),DIMENSION(3,3):: ORIENT      !crystallographic orientation (not used here but necessary for some calls)
REAL(dp),DIMENSION(3,3,3):: A_tensor  !tensor A(IM)
REAL(dp),DIMENSION(9,9):: C_tensor    !stiffness tensor (not used here but necessary for some calls)
REAL(dp),DIMENSION(:),ALLOCATABLE:: Stemp
REAL(dp),DIMENSION(:),ALLOCATABLE:: work_array
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: IdMat, Q_plus
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX  !final auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE,TARGET:: P1, P2
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S,V_NN
REAL(dp),DIMENSION(:,:),POINTER:: Ppoint !pointer to the system to build reference
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P_neigh, Q_neigh, P_matrix, Q_matrix, P_neigh_tmp, Q_matrix_copy
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: PosList1, PosList2 !list of positions of neighbors of one atom
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE:: Pref  !references for atoms environments
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE:: G, Delta_G, Delta_G_matrix, Delta_G_tmp 


!Initialize variables and arrays
Nsites = 0
firstref = .TRUE.
ucref = .FALSE.
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
Nnmin=3   !minimum number of neighbours to keep (may change depending on composition, see below)
Huc(:,:) = 0.d0
ORIENT(:,:) = 0.d0     !crystallographic orientation (not used, set to zero)
 C_tensor(:,:) = 0.d0  !stiffness tensor (not used, set to zero)
!
!
CALL ATOMSK_MSG(4061,(/""/),(/0.d0/))
!
!
!
100 CONTINUE
!**********************************************************************************
!                            READ  INPUT  FILES
!**********************************************************************************
IF( StrUpCase(filefirst)=="NULL" ) THEN
  !There is no reference file: atoms environment will be generated from the system2
  firstref = .FALSE.
ELSE
  !Read atomic positions from filefirst and store them into P1(:,:)
  CALL READ_AFF(filefirst,H1,P1,S,comment,AUXNAMES,AUX)
  IF(nerr>0 .OR. .NOT.ALLOCATED(P1)) GOTO 1000
  !Get rid of shells and auxiliary properties
  IF (ALLOCATED(S)) DEALLOCATE(S)
  IF (ALLOCATED(comment)) DEALLOCATE(comment)
  IF (ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
  IF (ALLOCATED(AUX)) DEALLOCATE(AUX)
  !Apply options to system 1
  CALL OPTIONS_AFF(options_array,Huc,H1,P1,S,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
ENDIF
!
!Read atomic positions from filesecond and store them into P2(:,:)
CALL READ_AFF(filesecond,H2,P2,S,comment,AUXNAMES,AUX)
IF(nerr>0 .OR. .NOT.ALLOCATED(P2)) GOTO 1000
!Get rid of shells and auxiliary properties
IF (ALLOCATED(S)) DEALLOCATE(S)
IF (ALLOCATED(comment)) DEALLOCATE(comment)
IF (ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF (ALLOCATED(AUX)) DEALLOCATE(AUX)
!Apply options to system 2
CALL OPTIONS_AFF(options_array,Huc,H2,P2,S,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
!
!
IF( ALLOCATED(P1) .AND. SIZE(P1,1)>0 ) THEN
  !A "reference" system was provided
  IF( SIZE(P1,1) < SIZE(P2,1) ) THEN
    !The system1 contains fewer atoms than the system2
    !Consider that the first system is a "unit cell"
    firstref = .FALSE.
    ucref = .TRUE.
  ELSE IF( SIZE(P1,1) == SIZE(P2,1) ) THEN
    !Systems 1 and 2 contain exactly the same number of atoms
    !Assume that system 1 is full reference (perfect crystal),
    !and that atoms have the same IDs in system 2
    firstref = .TRUE.
    ucref = .FALSE.
  ELSE
    !System 1 contains MORE atoms than system 2
    CALL ATOMSK_MSG(4810,(/""/),(/DBLE(SIZE(P1,1)) , DBLE(SIZE(P2,1))/))
    nerr=nerr+1
    GOTO 1000
  ENDIF
  !
ELSE  !i.e. P1 is NOT ALLOCATED
  !No "reference": we have only the deformed system (P2) to work with
  firstref = .FALSE.
  ucref = .FALSE.
ENDIF
!
!
!
200 CONTINUE
!**********************************************************************************
!         COMPUTE  LATTICE  CORRESPONDENCE  TENSOR  G  FOR  EACH  ATOM
!**********************************************************************************
CALL ATOMSK_MSG(4062,(/""/),(/0.d0/))
!
!The following routine performs the computation of G (see /src/compute/compute_G.f90)
CALL COMPUTE_G(H1,P1,H2,P2,params,NeighList1,NeighList2,Pref,siteindex,G)
!
IF(nerr>0) GOTO 1000
!
!Get parameters used during calculation of G
cutoff = params(1)
NNmin = params(2)
NeighFactor = params(3)
theta_max = params(4)
!
!
!
!**********************************************************************************
!         WRITE  G  AND  DERIVED  TENSORS  INTO  FILES
!**********************************************************************************
!Write atom coordinates and per-atom matrix G into a CFG file for visualization
IF( ALLOCATED(AUX) ) DEALLOCATE(AUX)
IF( ALLOCATED(AUXNAMES) ) DEALLOCATE(AUXNAMES)
IF( ALLOCATED(comment) ) DEALLOCATE(comment)
ALLOCATE( AUX( SIZE(P2,1) , 9 ) )
AUX(:,:) = 0.d0
AUX(:,1) = G(:,1,1)
AUX(:,2) = G(:,1,2)
AUX(:,3) = G(:,1,3)
AUX(:,4) = G(:,2,1)
AUX(:,5) = G(:,2,2)
AUX(:,6) = G(:,2,3)
AUX(:,7) = G(:,3,1)
AUX(:,8) = G(:,3,2)
AUX(:,9) = G(:,3,3)
ALLOCATE (AUXNAMES(9))
AUXNAMES(1)="G_11"
AUXNAMES(2)="G_12"
AUXNAMES(3)="G_13"
AUXNAMES(4)="G_21"
AUXNAMES(5)="G_22"
AUXNAMES(6)="G_23"
AUXNAMES(7)="G_31"
AUXNAMES(8)="G_32"
AUXNAMES(9)="G_33"
ALLOCATE(comment(1))
comment(1) = "# Per-atom G tensor computed by Atomsk"
!If user did not provide output file name or prefix, set a default one
msg = prefix
IF( LEN_TRIM(prefix)<=0 ) THEN
  msg = filesecond
ENDIF
i = SCAN(msg,".",BACK=.TRUE.)
IF(i==0) i=LEN_TRIM(msg)+1
msg = TRIM(ADJUSTL(msg(1:i-1)))//"_G.cfg"
CALL WRITE_AFF(msg,outfileformats,H2,P2,S,comment,AUXNAMES,AUX)
DEALLOCATE(AUX)
DEALLOCATE(AUXNAMES)
DEALLOCATE(comment)
!
!Write atom coordinates and strain tensor into a CFG file for visualization
ALLOCATE( AUX( SIZE(P2,1) , 12 ) )
!Components of strain tensor
AUX(:,:) = 0.d0
AUX(:,1) = 1.d0 - G(:,1,1)
AUX(:,2) = 1.d0 - G(:,2,2)
AUX(:,3) = 1.d0 - G(:,3,3)
AUX(:,4) = -1.d0 * ( G(:,1,2) + G(:,2,1) ) / 2.d0
AUX(:,5) = -1.d0 * ( G(:,1,3) + G(:,3,1) ) / 2.d0
AUX(:,6) = -1.d0 * ( G(:,2,3) + G(:,3,2) ) / 2.d0
!Strain invariants I1, I2, I3
AUX(:,7) = AUX(:,1) + AUX(:,2) + AUX(:,3)
AUX(:,8) = AUX(:,1)*AUX(:,2) + AUX(:,1)*AUX(:,3) + AUX(:,2)*AUX(:,3) - AUX(:,4)**2 - AUX(:,5)**2 - AUX(:,6)**2
DO i=1,SIZE(AUX,1)
  AUX(i,9) = MATDET(G(i,1:3,1:3))  ! I3 = det(G)
ENDDO
AUX(:,11) = AUX(:,7)**2 - AUX(:,8)
AUX(:,12) = (2.d0*AUX(:,7)**3)/27.d0 - AUX(:,7)*AUX(:,8)/3.d0 + AUX(:,9)
ALLOCATE (AUXNAMES(12))
AUXNAMES(1)="strain_11"
AUXNAMES(2)="strain_22"
AUXNAMES(3)="strain_33"
AUXNAMES(4)="strain_12"
AUXNAMES(5)="strain_13"
AUXNAMES(6)="strain_23"
AUXNAMES(7)="I1"
AUXNAMES(8)="I2"
AUXNAMES(9)="I3"
AUXNAMES(10)="J1"
AUXNAMES(11)="J2"
AUXNAMES(12)="J3"
ALLOCATE(comment(1))
comment(1) = "# Per-atom strain tensor computed by Atomsk"
!If user did not provide output file name or prefix, set a default one
msg = prefix
IF( LEN_TRIM(prefix)<=0 ) THEN
  msg = filesecond
ENDIF
i = SCAN(msg,".",BACK=.TRUE.)
IF(i==0) i=LEN_TRIM(msg)+1
msg = TRIM(ADJUSTL(msg(1:i-1)))//"_strain.cfg"
CALL WRITE_AFF(msg,outfileformats,H2,P2,S,comment,AUXNAMES,AUX)
DEALLOCATE(AUX)
DEALLOCATE(AUXNAMES)
DEALLOCATE(comment)
!
!Write atom coordinates and rotation tensor into a CFG file for visualization
ALLOCATE( AUX( SIZE(P2,1) , 6 ) )
AUX(:,:) = 0.d0
AUX(:,4) = ( G(:,2,1) - G(:,1,2) ) / 2.d0
AUX(:,5) = ( G(:,3,1) - G(:,1,3) ) / 2.d0
AUX(:,6) = ( G(:,3,2) - G(:,2,3) ) / 2.d0
ALLOCATE (AUXNAMES(6))
AUXNAMES(1)="rot_11"
AUXNAMES(2)="rot_22"
AUXNAMES(3)="rot_33"
AUXNAMES(4)="rot_12"
AUXNAMES(5)="rot_13"
AUXNAMES(6)="rot_23"
ALLOCATE(comment(1))
comment(1) = "# Per-atom rotation tensor computed by Atomsk"
!If user did not provide output file name or prefix, set a default one
msg = prefix
IF( LEN_TRIM(prefix)<=0 ) THEN
  msg = filesecond
ENDIF
i = SCAN(msg,".",BACK=.TRUE.)
IF(i==0) i=LEN_TRIM(msg)+1
msg = TRIM(ADJUSTL(msg(1:i-1)))//"_rot.cfg"
CALL WRITE_AFF(msg,outfileformats,H2,P2,S,comment,AUXNAMES,AUX)
DEALLOCATE(AUX)
DEALLOCATE(AUXNAMES)
DEALLOCATE(comment)
!
!
!
300 CONTINUE
!**********************************************************************************
!                            COMPUTE  NYE  TENSOR
!**********************************************************************************
!Tensor G is known for all atoms => use it to compute the Nye tensor
CALL ATOMSK_MSG(4063,(/""/),(/0.d0/))
!
IF( verbosity>=4 ) THEN
  WRITE(msg,*) "     cutoff = ", cutoff
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "      NNmin = ", NNmin
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "NeighFactor = ", NeighFactor
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) " theta_max  = ", theta_max
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
!
!Prepare array AUX to store Nye tensor
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
ALLOCATE( AUX(SIZE(P2,1),9) )
AUX(:,:) = 0.d0
!
!New Loop on the atoms
progress=0
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP& PRIVATE(msg,iat,i,j,k,m,n,ok,Nneighbors,nb_neigh,cutoff,PosList1,PosList2,P_neigh,P_neigh_tmp,Tab_PQ) &
!$OMP& PRIVATE(IdMat,Q_neigh,Q_plus,P_matrix,Q_matrix,Q_matrix_copy,Stemp,test_matrix,V_NN,newindex) &
!$OMP& PRIVATE(Q_matrix_rank,work_array,LWORK,INFO,tempreal,alpha,alpha_tmp) &
!$OMP& PRIVATE(Nlist,alpha_tensor,A_tensor,eps,Delta_G,Delta_G_matrix,Delta_G_tmp)
DO iat=1,SIZE(P2,1)
  !
  progress = progress+1
  !
  !Initialize Nye tensor for atom #iat
  alpha_tensor(:,:) = 0.d0
  !
  IF( SIZE(P2,1) > 5000 ) THEN
    !If there are many atoms, display a fancy progress bar
    CALL ATOMSK_MSG(10,(/""/),(/DBLE(progress),DBLE(SIZE(P2,1))/))
  ENDIF
  !
  A_tensor(:,:,:) = 0.d0
  IF(ALLOCATED(Nlist)) DEALLOCATE(Nlist)
  IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
  IF(ALLOCATED(P_neigh)) DEALLOCATE(P_neigh)
  IF(ALLOCATED(Q_neigh)) DEALLOCATE(Q_neigh)
  IF(ALLOCATED(Delta_G)) DEALLOCATE(Delta_G)
  !
  !
  IF( firstref ) THEN
    !The system in "filefirst" is used as reference
    !Search neighbors of atom #iat in first system
    CALL NEIGHBOR_POS(H1,P1,P1(iat,1:3),NeighList1(iat,:),ALLOCATED(NeighList1),radius,PosList1)
    !
    IF( SIZE(PosList1,1)>=3 ) THEN
      !Now PosList1(:,:) contains the cartesian positions of all neighbors of atom #iat,
      !their distance to the central atom #iat, and their indices in P(:,:).
      !Sort them by increasing distance:
      CALL BUBBLESORT(PosList1,4,'up  ',newindex)
      !Keep only the first neighbors, save them in V_NN
      !Make sure to keep at least NNmin neighbors: compare distances to that of neighbor #NNmin
      Nneighbors=0
      DO j=1,SIZE(PosList1,1)
        !IF( PosList1(j,4) <= NeighFactor*PosList1(NNmin,4) ) THEN
        !IF( PosList1(j,4) <= cutoff ) THEN
        !IF( PosList1(j,5) > 0 ) THEN
          !This neighbor is about as close as the 3rd neighbor => keep it
          Nneighbors=Nneighbors+1
        !ENDIF
      ENDDO
      !Save first neighbors positions into V_NN(:,:)
      !Save their index into Nlist(:)
      ALLOCATE(V_NN(Nneighbors,3))
      V_NN(:,:) = 0.d0
      ALLOCATE(Nlist(Nneighbors))
      Nlist(:) = 0
      Nneighbors=0
      DO j=1,SIZE(PosList1,1)
        !IF( PosList1(j,4) <= NeighFactor*PosList1(NNmin,4) ) THEN
        !IF( PosList1(j,4) <= cutoff ) THEN
        !IF( PosList1(j,5) > 0 ) THEN
          !This neighbor is closer than the NNmin-th neighbor => keep it
          Nneighbors=Nneighbors+1
          V_NN(Nneighbors,:) = PosList1(j,1:3)
          Nlist(Nneighbors) = NINT(PosList1(j,5))
        !ENDIF
      ENDDO
    ENDIF
    !
    !
    IF( .NOT. ALLOCATED(V_NN) .OR. SIZE(V_NN,1)<3 ) THEN
      !Not enough neighbors to perform calculation => skip
      CALL ATOMSK_MSG(4709,(/""/),(/ DBLE(iat) , DBLE(SIZE(V_NN,1)) /))
      nwarn = nwarn+1
      GOTO 390
      !
    ELSEIF( SIZE(V_NN,1)>100 ) THEN
      !Atom #iat has more than 100 neighbors => skip calculation
      CALL ATOMSK_MSG(4705,(/""/),(/DBLE(iat)/))
      nwarn = nwarn+1
      GOTO 390
      !
    ELSE
      !
      !Save relative positions of neighbors in P_neigh
      !and compute  Delta_G(IM) = G_neighbor(IM) - G_0(IM)  (Eq.19)
      IF(ALLOCATED(P_neigh)) DEALLOCATE(P_neigh)
      ALLOCATE (P_neigh(SIZE(V_NN,1),3))
      P_neigh(:,:)=0.d0
      ALLOCATE (Delta_G(SIZE(V_NN,1),3,3))
      Delta_G(:,:,:) = 0.d0
      !
      DO j=1,SIZE(V_NN,1)
        P_neigh(j,:) = V_NN(j,1:3) - P1(iat,1:3)
        !
        DO i=1,3
          DO m=1,3
            Delta_G(j,i,m) = G(Nlist(j),i,m) - G(iat,i,m)
          ENDDO
        ENDDO
        !
      ENDDO
      !
    ENDIF
    !
    IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
    !
    !Search neighbors of atom #iat in second system
    CALL NEIGHBOR_POS(H2,P2,P2(iat,1:3),NeighList2(iat,:),ALLOCATED(NeighList2),radius,PosList2)
    !
    IF( SIZE(PosList2,1)>=3 ) THEN
      !Now PosList2(:,:) contains the cartesian positions of all neighbors in the cutoff,
      !their distance to the atom #iat, and their indices.
      !Sort them by increasing distance:
      CALL BUBBLESORT(PosList2,4,'up  ',newindex)
      !
      !Keep only the first neighbors, save them in V_NN
      !Make sure to keep at least NNmin neighbors: compare distances to that of neighbor #NNmin
      Nneighbors=0
      DO j=1,SIZE(PosList2,1)
        !IF( PosList2(j,4) <= NeighFactor*PosList2(NNmin,4) ) THEN
        !IF( PosList2(j,4) <= cutoff ) THEN
        !IF( PosList2(j,5) > 0 ) THEN
          !This neighbor is about as close as the NNmin-th neighbor => keep it
          Nneighbors=Nneighbors+1
        !ENDIF
      ENDDO
      ALLOCATE(V_NN(Nneighbors,3))
      V_NN(:,:) = 0.d0
      Nneighbors=0
      DO j=1,SIZE(PosList2,1)
        !IF( PosList2(j,4) <= NeighFactor*PosList2(NNmin,4) ) THEN
        !IF( PosList2(j,4) <= cutoff ) THEN
        !IF( PosList1(j,5) > 0 ) THEN
          !This neighbor is closer than the NNmin-th neighbor => keep it
          Nneighbors=Nneighbors+1
          V_NN(Nneighbors,:) = PosList2(j,1:3)
        !ENDIF
      ENDDO
    ENDIF
    !
  ELSE
    !No reference was provided (filefirst was NULL)
    !The atom #iat is in a site of type siteindex(iat)
    !and neighbors of atom #iat are in Pref
    IF( siteindex(iat)>0 ) THEN
      !Get number of neighbors for this type of site
      Nneighbors = NINT(Pref(siteindex(iat),1,6))
      !
      !Save atoms in perfect environment in P_neigh
      ALLOCATE(P_neigh(Nneighbors,3))
      P_neigh(:,:) = 0.d0
      DO j=1,Nneighbors
        P_neigh(j,1:3) = Pref(siteindex(iat),j+1,1:3)
      ENDDO
      !
      !Search neighbors of atom #iat in second system
      CALL NEIGHBOR_POS(H2,P2,P2(iat,1:3),NeighList2(iat,:),ALLOCATED(NeighList2),radius,PosList2)
      !
      IF( SIZE(PosList2,1)>=NNmin ) THEN
        !Now PosList2(:,:) contains the cartesian positions of all neighbors of atom #iat,
        !their distance to the central atom #iat, and their indices in P2(:,:).
        !Sort them by increasing distance:
        CALL BUBBLESORT(PosList2,4,'up  ',newindex)
        !Keep only the first neighbors, save them in V_NN
        !Make sure to keep at least NNmin neighbors: compare distances to that of neighbor #NNmin
        nb_neigh=0
        DO j=1,SIZE(PosList2,1)
          !IF( PosList2(j,4) <= NeighFactor*PosList2(NNmin,4) ) THEN
          !IF( PosList2(j,4) <= cutoff ) THEN
            !This neighbor is about as close as the NNmin-th neighbor => keep it
            nb_neigh=nb_neigh+1
          !ENDIF
        ENDDO
        ALLOCATE(V_NN(nb_neigh,5))
        V_NN(:,:) = 0.d0
        ALLOCATE(Nlist(nb_neigh))
        Nlist(:) = 0
        nb_neigh = 0
        DO j=1,SIZE(PosList2,1)
          !IF( PosList2(j,4) <= NeighFactor*PosList2(NNmin,4) ) THEN
          !IF( PosList2(j,4) <= cutoff ) THEN
            !This neighbor is about as close as the NNmin-th neighbor => keep it
            nb_neigh=nb_neigh+1
            !Save neighbors from second system in V_NN
            V_NN(nb_neigh,1:3) = PosList2(j,1:3)
            !Save index of neighbors in Nlist
            Nlist(nb_neigh) = NINT(PosList2(j,5))
          !ENDIF
        ENDDO
        !
        !Compute  Delta_G(IM) = G_neighbor(IM) - G_0(IM)  (Eq.19)
        IF( nb_neigh>=3 ) THEN
          ALLOCATE (Delta_G(nb_neigh,3,3))
          Delta_G(:,:,:) = 0.d0
          DO j=1,nb_neigh
            DO i=1,3
              DO m=1,3
                Delta_G(j,i,m) = G(Nlist(j),i,m) - G(iat,i,m)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          !Not enough neighbors => skip calculation
          GOTO 390
        ENDIF
      ELSE
        !Not enough neighbors => skip calculation
        GOTO 390
      ENDIF
      !
    ELSE
      !No site was found for this atom; at this point, several warnings were already displayed
      !Just skip the calculation for this atom
      GOTO 390
    ENDIF
    !
  ENDIF  !end if firstref
  !
  IF( verbosity==4 ) THEN
    WRITE(msg,*) "atom #", iat, "  (", nb_neigh, " neighbors), Delta_G ="
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO j=1,nb_neigh
      WRITE(msg,'(3(3f9.3,a3))') Delta_G(j,1,1:3), " | ", Delta_G(j,2,1:3), " | ", Delta_G(j,3,1:3)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
  ENDIF
  !
  !
  IF( .NOT. ALLOCATED(V_NN) .OR. SIZE(V_NN,1)<3 ) THEN
    !Not enough neighbors to perform calculation => skip
    CALL ATOMSK_MSG(4709,(/""/),(/ DBLE(iat) , DBLE(SIZE(V_NN,1)) /))
    nwarn = nwarn+1
    !
  ELSEIF( SIZE(V_NN,1)>100 ) THEN
    !Atom #iat has more than 100 neighbors => skip calculation
    CALL ATOMSK_MSG(4705,(/""/),(/DBLE(iat)/))
    nwarn = nwarn+1
    !
  ELSE
    !
    !Save relative positions of neighbors in Q_neigh
    IF(ALLOCATED(Q_neigh)) DEALLOCATE(Q_neigh)
    ALLOCATE (Q_neigh(SIZE(V_NN,1),3))
    Q_neigh(:,:)=0.d0
    DO j=1,SIZE(V_NN,1)
      Q_neigh(j,:) = V_NN(j,1:3) - P2(iat,1:3)
    ENDDO
    DEALLOCATE(V_NN)
    !
    !Generate correspondance table between P_neigh and Q_neigh
    ALLOCATE (Tab_PQ(SIZE(Q_neigh,1)))
    Tab_PQ(:)=0
    !
    DO j=1,SIZE(Q_neigh,1)
      alpha=100
      DO k=1,SIZE(P_neigh,1)
        alpha_tmp=DABS(ANGVEC(Q_neigh(j,1:3),P_neigh(k,1:3)))
        IF( alpha_tmp < alpha ) THEN
          Tab_PQ(j)=k
          alpha=alpha_tmp
        ENDIF
      ENDDO
    ENDDO
    !
    IF(ALLOCATED(P_neigh_tmp)) DEALLOCATE(P_neigh_tmp)
    IF(ALLOCATED(Delta_G_tmp)) DEALLOCATE(Delta_G_tmp)
    ALLOCATE(P_neigh_tmp(SIZE(Q_neigh,1),3))
    ALLOCATE(Delta_G_tmp(SIZE(Q_neigh,1),3,3))
    P_neigh_tmp(:,:)=0.d0
    Delta_G_tmp(:,:,:)=0.d0
    !
    IF( firstref ) THEN
      !A full reference was provided
      DO j=1,SIZE(Q_neigh,1)
        P_neigh_tmp(j,1:3) = P_neigh(Tab_PQ(j),1:3)
        Delta_G_tmp(j,:,:) = Delta_G(Tab_PQ(j),:,:)
      ENDDO
    ELSE
      !No reference was provided
      DO j=1,SIZE(Q_neigh,1)
      !DO j=1,SIZE(Tab_PQ)
        P_neigh_tmp(j,1:3) = P_neigh(Tab_PQ(j),1:3)
        Delta_G_tmp(j,:,:) = Delta_G(j,:,:)
      ENDDO
    ENDIF
    DEALLOCATE(Tab_PQ)
    !
    !We won't need Delta_G anymore
    DEALLOCATE(Delta_G)
    !
    !At this stage, the tables P and Q are in the same order
    !Compute angles between corresponding vectors, and
    !exclude vectors if the angle is too large
    ALLOCATE (Tab_PQ(SIZE(Q_neigh,1)))
    Tab_PQ(:)=0
    nb_neigh=0
    DO j=1,SIZE(Q_neigh,1)
      alpha = DABS(ANGVEC(Q_neigh(j,1:3),P_neigh_tmp(j,1:3)))
      IF( alpha > theta_max ) THEN
        Tab_PQ(j)=0
      ELSE
        ok=0 
        !Loop on neighbors
        DO k=1,SIZE(P_neigh,1)
          alpha_tmp=DABS(ANGVEC(Q_neigh(j,1:3),P_neigh(k,1:3)))
          IF( DABS(alpha_tmp-alpha)<1.d-12 ) THEN
            ok=1  !=ok+1
            EXIT
          ENDIF
        ENDDO
        IF( ok==1 ) THEN
          nb_neigh=nb_neigh+1
          Tab_PQ(j)=nb_neigh
        ELSE
          Tab_PQ(j)=0
        ENDIF
      ENDIF
    ENDDO
    !
    !Copy Delta_G into Delta_G_matrix with correct P/Q indices
    IF(ALLOCATED(P_matrix)) DEALLOCATE(P_matrix)
    IF(ALLOCATED(Q_matrix)) DEALLOCATE(Q_matrix)
    ALLOCATE(Q_matrix(nb_neigh,3))
    Q_matrix(:,:)=0.d0
    IF(ALLOCATED(Delta_G_matrix)) DEALLOCATE(Delta_G_matrix)
    ALLOCATE(Delta_G_matrix(nb_neigh,3,3))
    Delta_G_matrix(:,:,:)=0.d0
    DO j=1,SIZE(Q_neigh,1)
      IF( Tab_PQ(j) > 0 ) THEN
        Q_matrix(Tab_PQ(j),:) = Q_neigh(j,:)
        Delta_G_matrix(Tab_PQ(j),:,:) = Delta_G_tmp(j,:,:)
      ENDIF
    ENDDO
    !
    DEALLOCATE(P_neigh,Q_neigh,P_neigh_tmp,Delta_G_tmp)
    DEALLOCATE(Tab_PQ)
    !
    !
    IF ( nb_neigh < 3 ) THEN
      !Not enough neighbors to perform calculation => skip
      CALL ATOMSK_MSG(4709,(/""/),(/ DBLE(iat) , DBLE(nb_neigh) /))
      nwarn = nwarn+1
      !
    ELSE
      !
      !Compute Q+ = (Q^T * Q)^-1 * Q^T   (Eq.18)
      !This is done thanks to the LAPACK subroutine DGELSS
      IF(ALLOCATED(IdMat)) DEALLOCATE(IdMat)
      ALLOCATE(IdMat(nb_neigh,nb_neigh))
      IdMat(:,:) = 0.d0
      DO j=1,nb_neigh
        IdMat(j,j) = 1.d0
      ENDDO
      !
      LWORK=3*MIN(nb_neigh,3) + MAX( 2*MIN(nb_neigh,3), MAX(nb_neigh,3), nb_neigh )
      ALLOCATE (work_array(MAX(1,LWORK)))
      ALLOCATE (Stemp(3))
      !
      CALL DGELSS(nb_neigh,3,nb_neigh,Q_matrix,nb_neigh,IdMat,MAX(nb_neigh,3),    &
                  & Stemp,-1.d0,Q_matrix_rank,work_array,LWORK,INFO)
      !
      IF (INFO.NE.0) THEN
        CALL ATOMSK_MSG(4710,(/"inverse of Q"/),(/ 0.d0 /))
        nwarn = nwarn+1
        !
      ELSE
        !
        !keep only the first 3 rows of output IdMat
        ALLOCATE(Q_plus(3,nb_neigh))
        Q_plus(:,:)=0.d0
        DO j=1,3
          Q_plus(j,1:nb_neigh) = IdMat(j,1:nb_neigh)
        ENDDO
        DEALLOCATE(IdMat,Q_matrix)
        !
        !Compute A(IM) = Q+ * Delta_G(IM)   (Eq.21)
        DO i=1,3
          DO m=1,3
            A_tensor(i,m,:) = MATMUL( Q_plus(1:3,:) , Delta_G_matrix(:,i,m) )
          ENDDO
        ENDDO
        !
        DEALLOCATE (Delta_G_matrix,Q_plus,work_array,Stemp)
        !
        !
        !Compute the Nye tensor:  alpha(jk) = epsilon(imk) * T(ijm)   (Eq.22)
        !where epsilon is the Levi-Civita permutation symbol, and
        !T(ijm) is the tensor of derivatives of G
        !NOTE: Eq.22 is wrong in the published article by Hartley et al.
        DO j=1,3
          DO k=1,3
            alpha_tensor(j,k) = 0.d0
            DO i=1,3
              DO m=1,3
                eps = EPS_LEVI_CIVITA(i,m,k)
                alpha_tensor(j,k) = alpha_tensor(j,k) + DBLE(eps)*A_tensor(i,j,m)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        !
        !Save the values of Nye tensor into AUX(:,:)
        DO i=1,9
          IF (i<=3) THEN
            j=1
            k=i
          ELSEIF (i<=6) THEN
            j=2
            k=i-3
          ELSE
            j=3
            k=i-6
          ENDIF
          AUX(iat,i) = alpha_tensor(j,k)
        ENDDO
        !
        !
      ENDIF  !end IF(INFO.NE.0)
      !
    ENDIF  !end IF(nb_neigh < 3)
    !
  ENDIF  !end IF(.NOT.ALLOCATED(V_NN))... for system 2
    !
  !ENDIF  !end IF(.NOT.ALLOCATED(V_NN))... for system 1
  !
  390 CONTINUE
  !
ENDDO ! End loop on all atoms iat
!$OMP END PARALLEL DO
!
!
!
400 CONTINUE
!Won't need P1 any longer
IF(ALLOCATED(P1)) DEALLOCATE (P1)
!
!Write final results into file(s)
ALLOCATE (AUXNAMES(9))
AUXNAMES(1)="Nye_11"
AUXNAMES(2)="Nye_12"
AUXNAMES(3)="Nye_13"
AUXNAMES(4)="Nye_21"
AUXNAMES(5)="Nye_22"
AUXNAMES(6)="Nye_23"
AUXNAMES(7)="Nye_31"
AUXNAMES(8)="Nye_32"
AUXNAMES(9)="Nye_33"
!
ALLOCATE(comment(1))
 comment(1) = "# Per-atom Nye tensor computed by Atomsk"
!
msg = prefix
!If user did not provide output file name or prefix, set a default one
IF( LEN_TRIM(prefix)<=0 ) THEN
  msg = filesecond
  i = SCAN(msg,".",BACK=.TRUE.)
  IF(i==0) i=LEN_TRIM(msg)
  msg = TRIM(ADJUSTL(msg(1:i-1)))//"_Nye"
  !In addition, if user did not specify an output format, use CFG by default
  IF( .NOT.ALLOCATED(outfileformats) .OR. SIZE(outfileformats)==0 ) THEN
    msg = TRIM(ADJUSTL(msg))//".cfg"
  ENDIF
ENDIF
!
!CALL WRITE_AFF(prefix,outfileformats,H1,P1,S,comment,AUXNAMES,AUX)
CALL WRITE_AFF(msg,outfileformats,H2,P2,S,comment,AUXNAMES,AUX)
!
DEALLOCATE (P2, AUX, AUXNAMES)
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE NYE_TENSOR
!
!
END MODULE mode_nye
