MODULE mode_nye
!
!**********************************************************************************
!*  MODE_NYE                                                                      *
!**********************************************************************************
!* This module reads two sets of atomic coordinates and                           *
!* calculate the corresponding nye tensor.                                        *
!* The algorithm follows the method described in:                                 *
!*    C.S. Hartley, Y. Mishin, Acta Mater. 53 (2005) 1313                         *
!* Equation numbers also refer to this reference.                                 *
!**********************************************************************************
!* (C) October 2013 - Philippe Carrez                                             *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     philippe.carrez@univ-lille1.fr                                             *
!* Last modification: P. Hirel - 28 May 2014                                      *
!**********************************************************************************
!* OUTLINE:                                                                       *
!* 100        Read atom positions systems 1 and 2, construct neighbor lists       *
!* 200        Compute tensor G for each atom                                      *
!* 300        Using G, compute gradient of G and Nye tensor for each atom         *
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
USE comv
USE subroutines
USE functions
USE messages
USE neighbors
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
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE,INTENT(IN):: outfileformats !list of output file formats
CHARACTER(LEN=22):: pbar
CHARACTER(LEN=128):: temp, temp2
CHARACTER(LEN=4096):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE :: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE :: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
INTEGER:: i, j, k, m, iat
INTEGER:: nb_neigh, eps
INTEGER:: Nneighbors
INTEGER:: ok
INTEGER:: Q_matrix_rank, INFO, LWORK
INTEGER,DIMENSION(:),ALLOCATABLE:: Nlist !index of neighbours
INTEGER,DIMENSION(:),ALLOCATABLE:: TAB_PQ
INTEGER,DIMENSION(:,:),ALLOCATABLE:: NeighList1, NeighList2  !list of neighbors for systems 1 and 2
REAL(dp):: alpha, alpha_tmp !angles between two vectors
REAL(dp):: NeighFactor !%of tolerance in the radius for neighbor search
REAL(dp):: radius
REAL(dp):: tempreal
REAL(dp),DIMENSION(3,3):: alpha_tensor, test_matrix 
REAL(dp),DIMENSION(3,3):: Hfirst,Hsecond !Box vectors of systems 1 and 2
REAL(dp),DIMENSION(3,3):: ORIENT   !crystallographic orientation of the system
REAL(dp),DIMENSION(3,3,3):: A_tensor  !tensor A(IM)
REAL(dp),DIMENSION(:),ALLOCATABLE:: Stemp
REAL(dp),DIMENSION(:),ALLOCATABLE:: work_array
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Id_matrix, Q_plus
REAL(dp),DIMENSION(:,:),ALLOCATABLE :: AUX  !final auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Pfirst,Psecond,S,V_NN
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P_neigh, Q_neigh, P_matrix, Q_matrix, P_neigh_tmp, Q_matrix_copy
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: PosList1, PosList2 !list of positions of neighbors of one atom
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE:: G, Delta_G, Delta_G_matrix, Delta_G_tmp 


!Initialize variables and arrays
NeighFactor = 1.1d0 !distance to 1st neighbor times 1.1 ensures to exclude second neighbor sphere in bcc metals
                    !1.1 is expected to be robust for simple lattices (fcc, bcc) but fails for complex
                    !or distorted systems
radius = 8.d0  !radius for neighbor search: 8 A should be enough to find some neighbors in any system
ORIENT(:,:) = 0.d0
!
!
CALL ATOMSK_MSG(4061,(/""/),(/0.d0/))
!
!
!
100 CONTINUE
!Read atomic positions from filefirst and store them into Pfirst(:,:)
CALL READ_AFF(filefirst,Hfirst,Pfirst,S,comment,AUXNAMES,AUX)
IF (ALLOCATED(S)) DEALLOCATE(S)
IF (ALLOCATED(comment)) DEALLOCATE(comment)
IF (ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF (ALLOCATED(AUX)) DEALLOCATE(AUX)
!Apply options to system 1
CALL OPTIONS_AFF(options_array,Hfirst,Pfirst,S,AUXNAMES,AUX,ORIENT)
!
!Read atomic positions from filesecond and store them into Psecond(:,:)
CALL READ_AFF(filesecond,Hsecond,Psecond,S,comment,AUXNAMES,AUX)
IF (ALLOCATED(S)) DEALLOCATE(S)
IF (ALLOCATED(comment)) DEALLOCATE(comment)
IF (ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF (ALLOCATED(AUX)) DEALLOCATE(AUX)
!Apply options to system 2
CALL OPTIONS_AFF(options_array,Hsecond,Psecond,S,AUXNAMES,AUX,ORIENT)
!
!
!Check that systems 1 and 2 have the same number of atoms
IF( SIZE(Pfirst,1) .NE. SIZE(Psecond,1) ) THEN
  CALL ATOMSK_MSG(4810,(/""/),(/DBLE(SIZE(Pfirst,1)) , DBLE(SIZE(Psecond,1))/))
  nerr=nerr+1
  GOTO 1000
ENDIF
!
!
!Search how many atom types exist in the system
CALL FIND_NSP(Pfirst(:,4),aentries)
IF( SIZE(aentries,1)>=3 ) THEN
  !Good chances that it is a complex material
  !=> the default NeighFactor will probably be too small, use larger value
  NeighFactor = 1.33d0
ENDIF
!
!
!Construct neighbor lists for systems 1 and 2
CALL ATOMSK_MSG(11,(/""/),(/0.d0/))
CALL NEIGHBOR_LIST(Hfirst,Pfirst,radius,NeighList1)
CALL NEIGHBOR_LIST(Hsecond,Psecond,radius,NeighList2)
!
IF( verbosity==4 ) THEN
  !Some debug messages
  WRITE(msg,*) "Size of neighbor list for SYSTEM 1: ", SIZE(NeighList1,1), SIZE(NeighList1,2)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "Sample of neighbor list for SYSTEM 1:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,MIN(10,SIZE(NeighList1,1))
    WRITE(msg,'(20i5)') i, NeighList1(i,1:MIN(SIZE(NeighList1,2),16))
    msg = TRIM(ADJUSTL(msg))//' (...)'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
  IF( i>=10 ) THEN
    WRITE(msg,*) '      (...discontinued...)'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
  WRITE(msg,*) "Size of neighbor list for SYSTEM 2: ", SIZE(NeighList2,1), SIZE(NeighList2,2)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "Sample of neighbor list for SYSTEM 2:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,MIN(10,SIZE(NeighList2,1))
    WRITE(msg,'(20i5)') i, NeighList2(i,1:MIN(SIZE(NeighList2,2),16))
    msg = TRIM(ADJUSTL(msg))//' (...)'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
  IF( i>=10 ) THEN
    WRITE(msg,*) '      (...discontinued...)'
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDIF
ENDIF
!
!
!
200 CONTINUE
CALL ATOMSK_MSG(4062,(/""/),(/0.d0/))
!
ALLOCATE(G(SIZE(Pfirst,1),3,3))
G(:,:,:) = 0.d0
!
!First, loop on all atoms to compute the tensor G for each atom
DO iat=1,SIZE(Pfirst,1)
  !
  IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
  !
  IF( SIZE(Pfirst,1) > 5000 ) THEN
    !If there are many atoms, display a fancy progress bar
    tempreal = 100.d0*DBLE(iat)/SIZE(Pfirst,1)
    CALL PROGBAR(tempreal,pbar)
    WRITE(temp,'(i3)') NINT(tempreal)
    WRITE(temp,*) "     "//TRIM(ADJUSTL(temp))//"% "//pbar
    CALL ATOMSK_MSG(10,(/TRIM(temp)/),(/0.d0/))
    IF(tempreal>=100.d0) WRITE(*,*) ''
  ENDIF
  !
  !Search for neighbors of atom #iat in the first system
  CALL NEIGHBOR_POS(Hfirst,Pfirst(:,:),Pfirst(iat,1:3),NeighList1(iat,:),radius,PosList1)
  !
  IF( SIZE(PosList1,1)>=3 ) THEN
    !Now PosList1(:,:) contains the cartesian positions of all neighbors in the radius,
    !their distance to the atom #iat, and their indices.
    !Sort them by increasing distance:
    CALL BUBBLESORT(PosList1,4,'up  ')
    !
    !Keep only the first neighbors, save them in V_NN
    !Make sure to keep at least 3 neighbors: compare distances to that of the 3rd neighbor
    Nneighbors=0
    DO j=1,SIZE(PosList1,1)
      IF( PosList1(j,4) <= NeighFactor*PosList1(3,4) ) THEN
        !This neighbor is about as close as the 3rd neighbor => keep it
        Nneighbors=Nneighbors+1
      ENDIF
    ENDDO
    ALLOCATE(V_NN(Nneighbors,3))
    V_NN(:,:) = 0.d0
    Nneighbors=0
    DO j=1,SIZE(PosList1,1)
      IF( PosList1(j,4) <= NeighFactor*PosList1(3,4) ) THEN
        !This neighbor is closer than the 3rd neighbor => keep it
        Nneighbors=Nneighbors+1
        V_NN(Nneighbors,:) = PosList1(j,1:3)
      ENDIF
    ENDDO
  ENDIF
  !
  WRITE(msg,*) 'SYSTEM 1: atom #, number of neighbors: ', iat, SIZE(V_NN,1)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  IF( .NOT. ALLOCATED(V_NN) .OR. SIZE(V_NN,1)<3 ) THEN
    !Not enough neighbors to perform calculation => skip
    CALL ATOMSK_MSG(4709,(/""/),(/ DBLE(iat) , DBLE(SIZE(V_NN,1)) /))
    nwarn = nwarn+1
    !
  ELSEIF ( SIZE(V_NN,1)>100 ) THEN
    !Atom #iat has more than 100 neighbors => skip calculation
    CALL ATOMSK_MSG(4705,(/""/),(/DBLE(iat)/))
    nwarn = nwarn+1
    !
  ELSE
    !
    !Save relative positions of neighbors of atom #iat into P_neigh(:,:)
    IF(ALLOCATED(P_neigh)) DEALLOCATE(P_neigh)
    ALLOCATE (P_neigh(SIZE(V_NN,1),3))
    P_neigh(:,:)=0.d0
    DO j=1,SIZE(V_NN,1)
      P_neigh(j,:) = V_NN(j,1:3)-Pfirst(iat,1:3)
    ENDDO
    !
    DEALLOCATE(V_NN)
    !
    !Now search neighbors of atom #iat in second system
    CALL NEIGHBOR_POS(Hsecond,Psecond,Psecond(iat,1:3),NeighList2(iat,:),radius,PosList2)
    IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
    !
    IF( SIZE(PosList2,1)>=3 ) THEN
      !Now PosList2(:,:) contains the cartesian positions of all neighbors in the radius,
      !and their distance to tha atom #i.
      !Sort them by increasing distance:
      CALL BUBBLESORT(PosList2,4,'up  ')
      !Keep only the first neighbors, save them in V_NN
      !Make sure to keep at least 3 neighbors: compare distances to that of the 3rd neighbor
      Nneighbors=0
      DO j=1,SIZE(PosList2,1)
        IF( PosList2(j,4) <= NeighFactor*PosList2(3,4) ) THEN
          !This neighbor is about as close as the 3rd neighbor => keep it
          Nneighbors=Nneighbors+1
        ENDIF
      ENDDO
      ALLOCATE(V_NN(Nneighbors,3))
      V_NN(:,:) = 0.d0
      Nneighbors=0
      DO j=1,SIZE(PosList2,1)
        IF( PosList2(j,4) <= NeighFactor*PosList2(3,4) ) THEN
          !This neighbor is closer than the 3rd neighbor => keep it
          Nneighbors=Nneighbors+1
          V_NN(Nneighbors,:) = PosList2(j,1:3)
        ENDIF
      ENDDO
    ENDIF
    !
    WRITE(msg,*) 'SYSTEM 2: atom #, number of neighbors: ', iat, SIZE(V_NN,1)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    !
    !
    IF( .NOT. ALLOCATED(V_NN) .OR. SIZE(V_NN,1)<3 ) THEN
      !Not enough neighbors to perform calculation => skip
      CALL ATOMSK_MSG(4709,(/""/),(/ DBLE(iat) , DBLE(SIZE(V_NN,1)) /))
      nwarn = nwarn+1
      !
    ELSEIF( SIZE(V_NN,1)>100 ) THEN
      CALL ATOMSK_MSG(4705,(/""/),(/DBLE(iat)/))
      nwarn = nwarn+1
      !
    ELSE
      !
      !Save positions of neighbors into the array Q_neigh
      IF(ALLOCATED(Q_neigh)) DEALLOCATE(Q_neigh)
      ALLOCATE (Q_neigh(SIZE(V_NN,1),3))
      Q_neigh(:,:)=0.d0
      DO j=1,SIZE(V_NN,1)
        Q_neigh(j,:) = V_NN(j,1:3) - Psecond(iat,1:3)
      ENDDO
      DEALLOCATE(V_NN)
      !
      IF (verbosity==4) THEN
        WRITE(msg,*) '-----Relative positions of neighbors-----'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        DO j=1,SIZE(P_neigh,1)
          WRITE(msg,'(a9,4f12.4)') 'P_neigh: ', P_neigh(j,1:3), VECLENGTH(P_neigh(j,1:3))
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        ENDDO  
        DO j=1,SIZE(Q_neigh,1)
          WRITE(msg,'(a9,4f12.4)') 'Q_neigh: ', Q_neigh(j,1:3), VECLENGTH(Q_neigh(j,1:3))
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        ENDDO  
        WRITE(msg,*) '-----------------------------------------'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDIF
      !
      !Build a table of correspondance between indices in P_neigh and those in Q_neigh
      ALLOCATE (Tab_PQ(SIZE(Q_neigh,1)))
      Tab_PQ(:)=0
      DO j=1,SIZE(Q_neigh,1)
        alpha=100
        DO k=1,SIZE(P_neigh,1)
          alpha_tmp=DABS(ANGVEC(Q_neigh(j,1:3),P_neigh(k,1:3)))
          IF (alpha_tmp.le.alpha) THEN 
            Tab_PQ(j)=k
            alpha=alpha_tmp
          ENDIF
        ENDDO
      ENDDO
      !
      !Matrix P_matrix and Q_matrix
      ALLOCATE(P_neigh_tmp(SIZE(Q_neigh,1),3))
      P_neigh_tmp(:,:)=0.d0
      DO j=1,SIZE(Q_neigh,1)
        P_neigh_tmp(j,1:3)=P_neigh(Tab_PQ(j),1:3)
      ENDDO
      !
      !At this stage, the tables P and Q are in the same order
      !Compute angles between corresponding vectors, and
      !exclude vectors if the angle is too large
      Tab_PQ(:)=0
      nb_neigh=0
      DO j=1,SIZE(Q_neigh,1)
        alpha=DABS(ANGVEC(Q_neigh(j,1:3),P_neigh_tmp(j,1:3)))
        IF (alpha.gt.0.6d0) THEN
          Tab_PQ(j)=0 
        ELSE
          ok=0 
          !boucle sur les voisins
          DO k=1,SIZE(P_neigh,1)
            alpha_tmp=DABS(ANGVEC(Q_neigh(j,1:3),P_neigh(k,1:3)))
            IF (alpha_tmp.eq.alpha) ok=ok+1
          ENDDO
          IF (ok.eq.1) THEN
            nb_neigh=nb_neigh+1
            Tab_PQ(j)=nb_neigh
          ELSE
            Tab_PQ(j)=0
          ENDIF
        ENDIF
      ENDDO
      !
      !Save final positions of neighbors into P_matrix and Q_matrix
      IF(ALLOCATED(P_matrix)) DEALLOCATE(P_matrix)
      IF(ALLOCATED(Q_matrix)) DEALLOCATE(Q_matrix)
      ALLOCATE(Q_matrix(nb_neigh,3))
      ALLOCATE(P_matrix(nb_neigh,3))
      Q_matrix(:,:)=0.d0
      P_matrix(:,:)=0.d0
      DO j=1,SIZE(Q_neigh,1)
        IF (Tab_PQ(j).ne.0) THEN
          Q_matrix(Tab_PQ(j),:) = Q_neigh(j,:)
          P_matrix(Tab_PQ(j),:) = P_neigh_tmp(j,:)
        ENDIF
      ENDDO
      !
      DEALLOCATE(P_neigh,Q_neigh,P_neigh_tmp)
      DEALLOCATE(Tab_PQ)
      !
      IF (verbosity==4) THEN
        WRITE(msg,*) '-----Neighbors used for calculation----'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        DO j=1,nb_neigh
          WRITE(msg,'(a9,4f12.4)') 'P_matrix:', P_matrix(j,1:3) 
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        ENDDO
        DO j=1,nb_neigh
          WRITE(msg,'(a9,4f12.4)') 'Q_matrix:', Q_matrix(j,1:3) 
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        ENDDO  
        WRITE(msg,*) '----------------------------------------'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDIF
      !
      !
      IF (nb_neigh < 3) THEN
        !Not enough neighbors to perform calculation => skip
        CALL ATOMSK_MSG(4709,(/""/),(/ DBLE(iat) , DBLE(nb_neigh) /))
        nwarn = nwarn+1
        !
      ELSE
        !
        !Compute Q+ = (Q^T * Q)^-1 * Q^T   (Eq.18)
        !This is done thanks to the LAPACK subroutine DGELSS
        ALLOCATE(Id_matrix(nb_neigh,nb_neigh))
        Id_matrix(:,:)=0.d0
        DO j=1,nb_neigh
          Id_matrix(j,j)=1.d0
        ENDDO
        !
        ALLOCATE(Q_matrix_copy(nb_neigh,3))
        Q_matrix_copy(:,:) = Q_matrix(:,:)
        !
        LWORK=3*min(nb_neigh,3) + max( 2*min(nb_neigh,3), max(nb_neigh,3), nb_neigh )
        ALLOCATE (work_array(max(1,LWORK)))
        ALLOCATE (Stemp(3))
        !
        CALL DGELSS(nb_neigh,3,nb_neigh,Q_matrix,nb_neigh,Id_matrix,MAX(nb_neigh,3),  &
                    & Stemp,-1.d0,Q_matrix_rank,work_array,LWORK,INFO)
        !
        IF ( INFO.NE.0 ) THEN
          !LAPACK routine returned an error
          CALL ATOMSK_MSG(4710,(/"inverse of Q"/),(/ 0.d0 /))
          nwarn = nwarn+1
          !
        ELSE
          !Keep only the first 3 rows of output I_matrix
          ALLOCATE(Q_plus(3,nb_neigh))
          Q_plus(:,:)=0.d0
          DO j=1,3
            Q_plus(j,1:nb_neigh) = Id_matrix(j,1:nb_neigh)
          ENDDO
          !
          !Test if  (Q+ * Qmatrix)  is identity matrix
          test_matrix(:,:) = MATMUL(Q_plus,Q_matrix_copy)
          IF (verbosity==4) THEN
            WRITE(msg,*)  'Q+ * Qmatrix should be id matrix, check:'
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
            WRITE(msg,*) test_matrix(1,:)
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
            WRITE(msg,*) test_matrix(2,:)
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
            WRITE(msg,*) test_matrix(3,:)
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          ENDIF
          ok=0
          DO i=1,3
            DO j=1,3
              tempreal = test_matrix(i,j)
              IF( i==j ) THEN
                tempreal = tempreal-1.d0
              ENDIF
              IF( DABS(tempreal)>1.d-6 ) THEN
                ok = 1
              ENDIF
            ENDDO
          ENDDO
          !
          DEALLOCATE(Id_matrix,Q_matrix)
          DEALLOCATE(Q_matrix_copy)
          !
          !
          IF( ok.NE.0 ) THEN
            nwarn=nwarn+1
            PRINT*, "NOT IDENTITY MATRIX!!!"
            !
          ELSE
            !
            !Compute  G = Q+ * P   (Eq.17)
            !Save the final tensor G for atom #iat
            G(iat,:,:) = MATMUL(Q_plus,P_matrix)
          ENDIF
          !
          !
          IF (verbosity==4) THEN
            WRITE(msg,*) 'TENSOR G FOR ATOM # ', iat
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
            WRITE(msg,*) "    ", G(iat,1,:)
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
            WRITE(msg,*) "    ", G(iat,2,:)
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
            WRITE(msg,*) "    ", G(iat,3,:)
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          ENDIF
          !
          DEALLOCATE (P_matrix,Q_plus,work_array,Stemp)
          !
          !
        ENDIF !end IF(INFO.NE.0)
        !
      ENDIF  !end IF(nb_neigh < 3)
      !
    ENDIF  !end IF(.NOT.ALLOCATED(V_NN))... for system 2
    !
  ENDIF  !end IF(.NOT.ALLOCATED(V_NN))... for system 1
  !
  !
ENDDO
!
IF( verbosity==4 ) THEN
  !Write atom coordinates and per-atom matrix G into a file
  ALLOCATE( AUX( SIZE(Pfirst,1) , 9 ) )
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
  !
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
  !
  ALLOCATE(comment(1))
  comment(1) = "# Per-atom G matrix computed by atomsk"
  !
  msg = "atomsk_G.cfg"
  CALL WRITE_AFF(msg,outfileformats,Hfirst,Pfirst,S,comment,AUXNAMES,AUX)
  !
  DEALLOCATE(AUX)
  DEALLOCATE(AUXNAMES)
  DEALLOCATE(comment)
  !
ENDIF
!
!
!
300 CONTINUE
!Tensor G is known for all atoms => use it to compute the Nye tensor
CALL ATOMSK_MSG(4063,(/""/),(/0.d0/))
!
ALLOCATE( AUX(SIZE(Psecond,1),9) )
!
!New Loop on the atom
DO iat=1,SIZE(Pfirst,1)
  !
  IF( SIZE(Pfirst,1) > 5000 ) THEN
    !If there are many atoms, display a fancy progress bar
    tempreal = 100.d0*DBLE(iat)/SIZE(Pfirst,1)
    CALL PROGBAR(tempreal,pbar)
    WRITE(temp,'(i3)') NINT(tempreal)
    WRITE(temp,*) "     "//TRIM(ADJUSTL(temp))//"% "//pbar
    CALL ATOMSK_MSG(10,(/TRIM(temp)/),(/0.d0/))
    IF(tempreal>=100.d0) WRITE(*,*) ''
  ENDIF
  !
  A_tensor(:,:,:)=0.d0
  !
  !Search neighbors of atom #iat in first system
  CALL NEIGHBOR_POS(Hfirst,Pfirst,Pfirst(iat,1:3),NeighList1(iat,:),radius,PosList1)
  IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
  !
  IF( SIZE(PosList1,1)>=3 ) THEN
    !Now PosList1(:,:) contains the cartesian positions of all neighbors in the radius,
    !their distance to the atom #iat, and their indices.
    !Sort them by increasing distance:
    CALL BUBBLESORT(PosList1,4,'up  ')
    !Keep only the first neighbors, save them in V_NN
    !Make sure to keep at least 3 neighbors: compare distances to that of the 3rd neighbor
    Nneighbors=0
    DO j=1,SIZE(PosList1,1)
      IF( PosList1(j,4) <= NeighFactor*PosList1(3,4) ) THEN
        !This neighbor is about as close as the 3rd neighbor => keep it
        Nneighbors=Nneighbors+1
      ENDIF
    ENDDO
    ALLOCATE(V_NN(Nneighbors,3))
    V_NN(:,:) = 0.d0
    IF(ALLOCATED(Nlist)) DEALLOCATE(Nlist)
    ALLOCATE(Nlist(Nneighbors))
    Nlist(:) = 0
    Nneighbors=0
    DO j=1,SIZE(PosList1,1)
      IF( PosList1(j,4) <= NeighFactor*PosList1(3,4) ) THEN
        !This neighbor is closer than the 3rd neighbor => keep it
        Nneighbors=Nneighbors+1
        V_NN(Nneighbors,:) = PosList1(j,1:3)
        Nlist(Nneighbors) = NINT(PosList1(j,5))
        IF( Nlist(Nneighbors)==0 ) THEN
          !It means that atom #iat is a neighbor of itself (because of PBC)
          Nlist(Nneighbors) = iat
        ENDIF
      ENDIF
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
    !Save relative positions of neighbors in P_neigh
    !and compute  Delta_G(IM) = G_neighbor(IM) - G_0(IM)  (Eq.19)
    IF(ALLOCATED(P_neigh)) DEALLOCATE(P_neigh)
    ALLOCATE (P_neigh(SIZE(V_NN,1),3))
    P_neigh(:,:)=0.d0
    ALLOCATE (Delta_G(SIZE(V_NN,1),3,3))
    Delta_G(:,:,:) = 0.d0
    !
    DO j=1,SIZE(V_NN,1)
      P_neigh(j,:) = V_NN(j,1:3) - Pfirst(iat,1:3)
      !
      DO i=1,3
        DO m=1,3
          Delta_G(j,i,m) = G(Nlist(j),i,m) - G(iat,i,m)
        ENDDO
      ENDDO
      !
    ENDDO
    !
    !
    !Search neighbors of atom #iat in second system
    CALL NEIGHBOR_POS(Hsecond,Psecond,Psecond(iat,1:3),NeighList2(iat,:),radius,PosList2)
    IF(ALLOCATED(V_NN)) DEALLOCATE(V_NN)
    !
    IF( SIZE(PosList2,1)>=3 ) THEN
      !Now PosList2(:,:) contains the cartesian positions of all neighbors in the radius,
      !their distance to the atom #iat, and their indices.
      !Sort them by increasing distance:
      CALL BUBBLESORT(PosList2,4,'up  ')
      !Keep only the first neighbors, save them in V_NN
      !Make sure to keep at least 3 neighbors: compare distances to that of the 3rd neighbor
      Nneighbors=0
      DO j=1,SIZE(PosList2,1)
        IF( PosList2(j,4) <= NeighFactor*PosList2(3,4) ) THEN
          !This neighbor is about as close as the 3rd neighbor => keep it
          Nneighbors=Nneighbors+1
        ENDIF
      ENDDO
      ALLOCATE(V_NN(Nneighbors,3))
      V_NN(:,:) = 0.d0
      Nneighbors=0
      DO j=1,SIZE(PosList2,1)
        IF( PosList2(j,4) <= NeighFactor*PosList2(3,4) ) THEN
          !This neighbor is closer than the 3rd neighbor => keep it
          Nneighbors=Nneighbors+1
          V_NN(Nneighbors,:) = PosList2(j,1:3)
        ENDIF
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
        Q_neigh(j,:) = V_NN(j,1:3) - Psecond(iat,1:3)
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
          IF (alpha_tmp.le.alpha) THEN 
            Tab_PQ(j)=k
            alpha=alpha_tmp
          ENDIF
        ENDDO
      ENDDO
      !
      ALLOCATE(P_neigh_tmp(SIZE(Q_neigh,1),3))
      ALLOCATE(Delta_G_tmp(SIZE(Q_neigh,1),3,3))
      P_neigh_tmp(:,:)=0.d0
      Delta_G_tmp(:,:,:)=0.d0
      !
      DO j=1,SIZE(Q_neigh,1)
        P_neigh_tmp(j,1:3) = P_neigh(Tab_PQ(j),1:3)
        Delta_G_tmp(j,:,:) = Delta_G(Tab_PQ(j),:,:)
      ENDDO
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
        IF (alpha.gt.0.6d0) THEN
          Tab_PQ(j)=0 
        ELSE
          ok=0 
          !Loop on neighbors
          DO k=1,SIZE(P_neigh,1)
            alpha_tmp=DABS(ANGVEC(Q_neigh(j,1:3),P_neigh(k,1:3)))
            IF (alpha_tmp==alpha) ok=ok+1
          ENDDO
          IF (ok==1) THEN
            nb_neigh=nb_neigh+1
            Tab_PQ(j)=nb_neigh
          ELSE
            Tab_PQ(j)=0
          ENDIF
        ENDIF
      ENDDO
      !
      !Compute Delta_G(IM) = Q * A(IM)   (Eq.20)
      IF(ALLOCATED(P_matrix)) DEALLOCATE(P_matrix)
      IF(ALLOCATED(Q_matrix)) DEALLOCATE(Q_matrix)
      ALLOCATE(Q_matrix(nb_neigh,3))
      Q_matrix(:,:)=0.d0
      IF(ALLOCATED(Delta_G_matrix)) DEALLOCATE(Delta_G_matrix)
      ALLOCATE(Delta_G_matrix(nb_neigh,3,3))
      Delta_G_matrix(:,:,:)=0.d0
      DO j=1,SIZE(Q_neigh,1)
        IF (Tab_PQ(j).NE.0) THEN
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
        IF(ALLOCATED(Id_matrix)) DEALLOCATE(Id_matrix)
        ALLOCATE(Id_matrix(nb_neigh,nb_neigh))
        Id_matrix(:,:)=0.d0
        DO j=1,nb_neigh
          Id_matrix(j,j)=1.d0
        ENDDO
        !
        LWORK=3*min(nb_neigh,3) + max( 2*min(nb_neigh,3), max(nb_neigh,3), nb_neigh )
        ALLOCATE (work_array(max(1,LWORK)))
        ALLOCATE (Stemp(3))
        !
        CALL DGELSS(nb_neigh,3,nb_neigh,Q_matrix,nb_neigh,Id_matrix,MAX(nb_neigh,3),    &
                    & Stemp,-1.d0,Q_matrix_rank,work_array,LWORK,INFO)
        !
        IF (INFO.NE.0) THEN
          CALL ATOMSK_MSG(4710,(/"inverse of Q"/),(/ 0.d0 /))
          nwarn = nwarn+1
          !
        ELSE
          !
          !keep only the first 3 rows of output I_matrix
          ALLOCATE(Q_plus(3,nb_neigh))
          Q_plus(:,:)=0.d0
          DO j=1,3
            Q_plus(j,1:nb_neigh)=Id_matrix(j,1:nb_neigh)
          ENDDO
          DEALLOCATE(Id_matrix,Q_matrix)
          !
          !Compute A(IM) = Q+ * Delta_G(IM)   (Eq.21)
          DO i=1,3
            DO m=1,3
              A_tensor(i,m,:) = MATMUL(Q_plus,Delta_G_matrix(:,i,m))
            ENDDO
          ENDDO
          !
          DEALLOCATE (Delta_G_matrix,Q_plus,work_array,Stemp)
          !
          !
          !Compute the Nye tensor:  alpha(jk) = -epsilon(jim) * T(imk)   (Eq.22)
          DO j=1,3
            DO k=1,3
              alpha_tensor(j,k)=0.d0
              DO i=1,3
                DO m=1,3
                  eps = EPS_LEVI_CIVITA(j,i,m)
                  alpha_tensor(j,k) = alpha_tensor(j,k) - eps*A_tensor(i,m,k)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          !
          !Save the values of Nye tensor into AUX(:,:)
          DO i=1,9
            IF (i.le.3) THEN
              j=1
              k=i
            ELSEIF (i.le.6) THEN
              j=2
              k=i-3
            ELSE
              j=3
              k=i-6
            ENDIF
            AUX(iat,i)=alpha_tensor(j,k)
          ENDDO
          !
          !
        ENDIF  !end IF(INFO.NE.0)
        !
      ENDIF  !end IF(nb_neigh < 3)
      !
    ENDIF  !end IF(.NOT.ALLOCATED(V_NN))... for system 2
    !
  ENDIF  !end IF(.NOT.ALLOCATED(V_NN))... for system 1
  !
  !
ENDDO ! End loop on all atoms at_nb
!
!
!
400 CONTINUE
!Write final results into file(s)
!DEALLOCATE (Pfirst)
ALLOCATE (AUXNAMES(9))
!
AUXNAMES(1)="alpha_11"
AUXNAMES(2)="alpha_12"
AUXNAMES(3)="alpha_13"
AUXNAMES(4)="alpha_21"
AUXNAMES(5)="alpha_22"
AUXNAMES(6)="alpha_23"
AUXNAMES(7)="alpha_31"
AUXNAMES(8)="alpha_32"
AUXNAMES(9)="alpha_33"
!
ALLOCATE(comment(1))
comment(1) = "# Per-atom G matrix computed by atomsk"
!
CALL WRITE_AFF(prefix,outfileformats,Hfirst,Pfirst,S,comment,AUXNAMES,AUX)
!CALL WRITE_AFF(prefix,outfileformats,Hsecond,Psecond,S,comment,AUXNAMES,AUX)
!
DEALLOCATE (Psecond, AUX, AUXNAMES)
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
