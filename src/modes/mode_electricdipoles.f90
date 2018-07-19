MODULE edm
!
!**********************************************************************************
!* MODE_EDM                                                                       *
!**********************************************************************************
!* This mode is meant to compute electric dipole moments in ionic systems.        *
!* The user has to provide the ionic species "Pspecies" forming the polyhedra,    *
!* as well as the coordination number (NNN) between these two species.            *
!* For instance if sp1 is expected to be in tetrahedra formed by sp2,             *
!* then NNN=4. For octahedra one would have NNN=6, etc. If the user               *
!* provides NNN<0 then all neighbours within a distance of DABS(NNN) will         *
!* be considered as part of the polyhedron. If NNN=0 is provided then the         *
!* program will attempt to automatically find the nearest neighbours.             *
!* (for details see FIND_NNN, FIND_NNR and FIND_1NN in subroutines.f90).          *
!* For each polyhedron the electric dipole moment is computed as:                 *
!*       P(r) =  Q1(r1-rcm) + SUM [ Q2 (ri-rcm) ]                                 *
!* where rcm is the position of the center of mass of the polyhedron, and         *
!* the SUM runs over all atoms of type 2 forming the polyhedron.                  *
!**********************************************************************************
!* (C) Jan. 2011 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 01 March 2017                                    *
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
!
!Load modules
USE atoms
USE comv        !global variables
USE constants   !math and physics constants
USE functions   !functions used by the program
USE messages
USE neighbors
USE subroutines !subroutines for this program
USE out_cfg     !For writing CFG file
USE out_xsf     !For writing XSF file
!
!
CONTAINS
!
SUBROUTINE E_DIPOLES(H,P,S,AUXNAMES,AUX,NNN,Pspecies,comment,outputfile)
!
!Declare variables
IMPLICIT NONE
CHARACTER(LEN=2),INTENT(IN):: Pspecies  !atom species forming the polyhedra (type #2)
CHARACTER(LEN=2):: species
CHARACTER(LEN=16):: ion1, ion2
CHARACTER(LEN=128):: msg, temp
CHARACTER(LEN=128):: outputfile, pxsf, pcfg, pnorm, pstat  !Output file names
CHARACTER(LEN=128):: poladat, polaxsf, ptot                !Output file names
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: CFGEDMNAME !name of aux.propr. for CFG file
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: fakeAUXNAMES !fake names for polarization vectors
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: autodetect  !autodetect polyhedra?
LOGICAL:: exceeds100 !does the number of neighbours exceed 100?
LOGICAL:: shells !are ions composed of core+shell?
INTEGER:: i, j, k
INTEGER:: qcol, qscol  !column containing charges in AUX array
INTEGER:: N1, N2   !number of particles #1 and #2 in the system
INTEGER,DIMENSION(:),ALLOCATABLE:: Nlist !index of neighbours (not used here)
REAL(dp):: mi, M, A, D, T  !for statistics
REAL(dp):: mass1, mass2 !mass of atoms of type #1 and #2
REAL(dp):: Q1, Q2 !charge of current species#1, of Pspecies
REAL(dp):: sp1number, sp2number  !atomic number of atoms of type #1, type #2 (Pspecies)
REAL(dp):: totcharge  !total charge of the system
REAL(dp):: Volume  !volume of the supercell
REAL(dp),INTENT(IN):: NNN     !number of particles #2 that are neighbours of particle 1
REAL(dp),DIMENSION(3):: center  !point of observation
REAL(dp),DIMENSION(3):: totpola  !total polarization
REAL(dp),DIMENSION(:),ALLOCATABLE:: normpola  !norm of polarization
REAL(dp), DIMENSION(3,3):: H, G !Base vectors of the supercell, inverse
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: CFGEDMAUX
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: P, S  !positions of ions, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P2  !positions of ions of type #2
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: pola  !polarization vectors
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: V_NN  !positions of 1st NN
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX  !auxiliary properties (should contain charges)
!
!Initialize variables
autodetect = .FALSE.
IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
  shells = .TRUE.
ELSE
  shells = .FALSE.
ENDIF
 qcol = 0
 qscol = 0
N1 = 0
N2 = 0
Q1 = 0.d0
Q2 = 0.d0
totcharge = 0.d0
totpola(:) = 0.d0
 msg = ''
ptot = TRIM(outputfile)//'_edm_Ptot'//'.txt'
temp = TRIM(outputfile)//'_edm'
polaxsf = TRIM(temp)//'_all.xsf'
poladat = TRIM(temp)//'_all.dat'
pxsf = TRIM(temp)//'.xsf'
pcfg = TRIM(temp)//'.cfg'
pnorm = TRIM(temp)//'.norm'
pstat = TRIM(temp)//'.stat'
IF(ALLOCATED(comment)) DEALLOCATE(comment)
ALLOCATE(comment(1))
 comment(1)=''
!
!
!
CALL ATOMSK_MSG(4032,(/''/),(/0.d0/))
!
!Invert matrix H
CALL INVMAT(H,G)
!
!Run some checks
!
!Find the column(s) containing charges in AUXNAMES
IF(ALLOCATED(AUX)) THEN
  i=0
  DO i=1,SIZE(AUXNAMES)
    IF( TRIM(ADJUSTL(AUXNAMES(i)))=='q' ) THEN
      qcol = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='qs' ) THEN
      qscol = i
    ENDIF
  ENDDO
ENDIF
WRITE(msg,*) "columns for q, qs: ", qcol, qscol
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!If no electric charges are specified, abort
IF(.NOT.ALLOCATED(AUX) .OR. qcol==0) THEN
  nerr=nerr+1
  CALL ATOMSK_MSG(4807,(/msg/),(/0.d0/))
  GOTO 1000
ENDIF
!
!Check consistency of shells
IF( qscol>0 .AND. .NOT.shells ) THEN
  !If shell charges are defined but there are no shells, don't use shells
  qscol = 0
ELSEIF( shells.AND.qscol==0 ) THEN
  !If shells are present but have zero charge, abort
  nerr=nerr+1
  CALL ATOMSK_MSG(4816,(/''/),(/DBLE(qcol),DBLE(qscol)/))
  GOTO 1000
ENDIF
!
!Set mass of atoms of type #2
CALL ATOMMASS(Pspecies,mass2)
!Set atomic number of atoms of type #2
CALL ATOMNUMBER(Pspecies,sp2number)
!Determine Q2 = charge of "type 2" (Pspecies) ions
i=0
DO WHILE(Q2==0.d0)
  i=i+1
  IF( DABS(P(i,4)-sp2number)<1.d-9 ) THEN
    Q2 = AUX(i,qcol)
    IF(qscol>0 .AND. NINT(S(i,4)).NE.0) Q2 = Q2+AUX(i,qscol)
    !...but if charge is zero we are in trouble
    !and we must exit the infinite loop
    IF(Q2==0.d0) THEN
      nerr=nerr+1
      CALL ATOMSK_MSG(4808,(/TRIM(Pspecies)/),(/0.d0/))
      GOTO 1000
    ENDIF
  ENDIF
ENDDO
WRITE(msg,*) "Q2 = ", Q2
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!Count atoms of type #2 in the system
CALL FIND_NSP(P(:,4),aentries)
DO i=1,SIZE(aentries,1)
  IF( DABS(aentries(i,1)-sp2number)<1.d-9 ) THEN
    N2 = NINT(aentries(i,2))
  ENDIF
ENDDO
WRITE(msg,*) "N2 = ", N2
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!If N2 is zero we are in trouble
IF(N2==0) THEN
  nerr = nerr+1
  CALL ATOMSK_MSG(4809,(/TRIM(ion2)/),(/0.d0/))
  GOTO 1000
ENDIF
!
!
!
100 CONTINUE
!Compute total electric polarization of the cell
totcharge = 0.d0
totpola(:) = 0.d0
DO i=1,SIZE(P,1)
  !Add the charge of current atom
  totcharge = totcharge + AUX(i,qcol)
  IF(qscol>0 .AND. NINT(S(i,4)).NE.0) totcharge = totcharge + AUX(i,qscol)
  !
  !Add contribution to total polarization
  totpola(:) = totpola(:) + AUX(i,qcol)*P(i,1:3)
  IF(qscol>0 .AND. NINT(S(i,4)).NE.0) totpola(:) = totpola(:) + AUX(i,qscol)*S(i,1:3)
ENDDO
CALL VOLUME_PARA(H,Volume)
totpola(:) = totpola(:)/Volume
!
!Write it to a file
IF(.NOT.overw) CALL CHECKFILE(ptot,'writ')
OPEN(UNIT=40,FILE=ptot)
WRITE(40,*) 'Total charge of the cell: ', totcharge
WRITE(40,*) 'Total polarization vector:', totpola(1), totpola(2), totpola(3)
!
CALL ATOMSK_MSG(4033,(/TRIM(ptot)/),(/0.d0/))
!
IF(totcharge>1.d-3) THEN
  nwarn=nwarn+1
  CALL ATOMSK_MSG(4702,(/''/),(/0.d0/))
  WRITE(40,*) ''
  WRITE(40,*) TRIM(msg)
ENDIF
CLOSE(40)
!
!
!
200 CONTINUE
!Compute individual electric dipole moments
!
!Store positions of atoms of type #2 in a separate array
IF(ALLOCATED(P2)) DEALLOCATE(P2)
ALLOCATE( P2(N2,4) )
N2 = 0
DO i=1,SIZE(P,1)
  IF( DABS(P(i,4)-sp2number)<1.d-9 ) THEN
    N2 = N2+1
    P2(N2,:) = P(i,:)
  ENDIF
ENDDO
!
!
!Allocate array to store polarizations
IF(ALLOCATED(pola)) DEALLOCATE(pola)
ALLOCATE( pola(SIZE(P,1),3) )
pola(:,:) = 0.d0
!Array to store the norms (current species only)
IF(ALLOCATED(normpola)) DEALLOCATE(normpola)
ALLOCATE( normpola(SIZE(P,1)) )
normpola(:) = 0.d0
!
totpola(:) = 0.d0
!
!Loop on all types of ions
DO i=1,SIZE(aentries,1)
  !
  !We consider only ions different from Pspecies
  IF( aentries(i,1).NE.sp2number ) THEN
    !Initialize variables
    Q1 = 0.d0
    !
    !Count atoms of type #1
    N1 = NINT(aentries(i,2))
    IF(N1==0) THEN
      nwarn = nwarn+1
      CALL ATOMSK_MSG(2723,(/TRIM(ion1)/),(/0.d0/))
      GOTO 250
    ENDIF
    !
    !Set mass of atoms of type #1
    CALL ATOMSPECIES(aentries(i,1),species)
    CALL ATOMMASS(species,mass1)
    !
    !Determine Q1 = charge of ions of type #1
    CALL ATOMNUMBER(species,sp1number)
    j=0
    DO WHILE(Q1==0.d0)
      j=j+1
      IF( DABS(P(j,4)-sp1number)<1.d-9 ) THEN
        Q1 = AUX(j,qcol)
        IF(qscol>0 .AND. NINT(S(j,4))>0) Q1 = Q1 + AUX(j,qscol)
        !...but if charge is zero we are in trouble
        !and we must exit the infinite loop
        IF(Q1==0.d0) THEN
          EXIT
        ENDIF
      ENDIF
    ENDDO
    !
    !Setup the text messages
    CALL ATOMSPECIES(aentries(i,1),ion1)
    CALL ATOMSK_MSG(4034,(/TRIM(ion1),TRIM(Pspecies)/),(/NNN/))
    IF( N1>5000 ) THEN
      CALL ATOMSK_MSG(3,(/''/),(/0.d0/))
    ENDIF
    !
    IF(Q1==0.d0) THEN
      nwarn = nwarn+1
      CALL ATOMSK_MSG(4703,(/TRIM(ion1),TRIM(ion2)/),(/NNN/))
      !
    ELSE
      !Name individual files
      temp = TRIM(outputfile)//'_edm_'//TRIM(species)//TRIM(Pspecies)
      pxsf = TRIM(temp)//'.xsf'
      pnorm = TRIM(temp)//'.norm'
      !
      !Write header of XSF file
      IF(.NOT.overw) CALL CHECKFILE(pxsf,'writ')
      OPEN(UNIT=40,FILE=pxsf)
      WRITE(40,*) "# Electric dipole moments for "               &
                & //TRIM(species)//TRIM(Pspecies)//" polyhedra"  &
                & //" computed with atomsk"
      WRITE(40,*) "CRYSTAL"
      WRITE(40,*) "PRIMVEC"
      WRITE(40,210) H(1,1), H(1,2), H(1,3)
      WRITE(40,210) H(2,1), H(2,2), H(2,3)
      WRITE(40,210) H(3,1), H(3,2), H(3,3)
      WRITE(40,*) "CONVVEC"
      WRITE(40,210) H(1,1), H(1,2), H(1,3)
      WRITE(40,210) H(2,1), H(2,2), H(2,3)
      WRITE(40,210) H(3,1), H(3,2), H(3,3)
      WRITE(40,*) "PRIMCOORD"
      WRITE(40,'(i9,a2)') N1, ' 1'
      210 FORMAT(3(f16.8,2X))
      !
      !Norms of polarizations
      IF(.NOT.overw) CALL CHECKFILE(pnorm,'writ')
      OPEN(UNIT=41,FILE=pnorm)
      CALL ATOMSPECIES(aentries(i,1),species)
      !
      !Compute the ionic polarization vectors
      N1=0
      !Loop on all atoms
      DO j=1,SIZE(P,1)
        !
        !We want only atoms of type 1
        IF( DABS(aentries(i,1)-P(j,4))<1.d-9 ) THEN
          N1=N1+1 !counter for atoms of type #1
          !
          IF( aentries(i,2)>5000 ) THEN
            !If there are many atoms, display a fancy progress bar
            CALL ATOMSK_MSG(10,(/""/),(/DBLE(N1),aentries(i,2)/))
          ENDIF
          !
          !Find nearest neighbours of current atom j
          !Search is performed only among P2, i.e. only for atoms of type #2
          IF(NNN>0.d0) THEN
            !Search N nearest neighbours
            CALL FIND_NNN(H,P2,P(j,:),NINT(NNN),V_NN,Nlist,exceeds100)
          ELSEIF(NNN<0.d0) THEN
            !Search all neighbours in the radius DABS(NNN)
            CALL FIND_NNRdR(H,P2,P(j,:),0.d0,DABS(NNN),V_NN,Nlist,exceeds100)
          ELSE
            !Automatically find the first nearest neighbours
            CALL FIND_1NN(H,P2,P(j,:),V_NN,Nlist,exceeds100)
          ENDIF
          !
          IF( verbosity==4 ) THEN
            WRITE(msg,*) j
            WRITE(msg,*) 'Atom #'//TRIM(ADJUSTL(msg)), P(j,1), P(j,2), P(j,3)
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
            WRITE(msg,*) SIZE(V_NN(:,1))
            WRITE(msg,*) 'Number of neighbours: '//TRIM(ADJUSTL(msg))
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
            DO k=1,SIZE(V_NN(:,1))
              WRITE(msg,'(i2,3f9.3)') k, V_NN(k,1), V_NN(k,2), V_NN(k,3)
              CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
            ENDDO
          ENDIF
          !
          IF(exceeds100) THEN
            nwarn=nwarn+1
            CALL ATOMSK_MSG(4705,(/''/),(/DBLE(j)/))
            !
          ELSEIF(SIZE(V_NN(:,1)).NE.0) THEN
            !Position of the center of mass of current polyhedron
            center(:) = mass1*P(j,1:3)
            DO k=1,SIZE(V_NN,1)
              center(:) = center(:) + mass2*V_NN(k,1:3)
            ENDDO
            center(:) = center(:) / ( mass1 + mass2*SIZE(V_NN,1) )
            !
            !Compute the polarization vector
            pola(j,1:3) = Q1*(P(j,1:3)-center(1:3))
            DO k=1,SIZE(V_NN,1)
              pola(j,1:3) = pola(j,1:3) + Q2*(V_NN(k,1:3)-center(:))
            ENDDO
            !
            !Compute norm of current pola. vector
            normpola(j) = VECLENGTH(pola(j,:))
            !
            !Compute total polarization
            totpola(:) = totpola(:)+pola(j,:)
            !
            !Output
            !XSF file
            WRITE(40,'(a2,6f16.8)') species, (P(j,k),k=1,3), (pola(j,k),k=1,3)
            !norm file
            WRITE(41,*) P(j,1), P(j,2), P(j,3), normpola(j)
            !
          ELSE
            !SIZE(V_NN) is zero => no neighbour was found
            nwarn=nwarn+1
            CALL ATOMSK_MSG(4706,(/TRIM(Pspecies)/),(/DBLE(j)/))
          ENDIF
        ENDIF
      ENDDO
      !
      CLOSE(40)
      CALL ATOMSK_MSG(4035,(/pxsf/),(/0.d0/))
      CALL ATOMSK_MSG(4036,(/pnorm/),(/0.d0/))
      !
      !
    ENDIF  !End IF(Q1==0)
    !
    !
    250 CONTINUE
    !
    !
  ENDIF
ENDDO  !End of big loop on all species
!
!
CALL ATOMSK_MSG(4044,(/""/),(/0.d0/))
!
!Write all polarization vectors to DAT file
IF(.NOT.overw) CALL CHECKFILE(poladat,'writ')
OPEN(UNIT=40,FILE=poladat)
DO i=1,SIZE(P,1)
  WRITE(40,'(6f16.8)') (P(i,j),j=1,3), (1.d0*pola(i,j),j=1,3)
ENDDO
CLOSE(40)
CALL ATOMSK_MSG(4037,(/poladat/),(/0.d0/))
!
!Write all polarization vectors to XSF file
!Note: we use fake names for auxiliary properties to make the
!     routine WRITE_XSF believe that pola vectors are forces,
!     otherwise it will not write it in the XSF file
 comment(1) = '# Electric dipole moments computed with atomsk'
ALLOCATE(fakeAUXNAMES(3))
fakeAUXNAMES(1)="fx"
fakeAUXNAMES(2)="fy"
fakeAUXNAMES(3)="fz"
IF(.NOT.overw) CALL CHECKFILE(polaxsf,'writ')
CALL WRITE_XSF(H,P,comment,fakeAUXNAMES,pola,polaxsf)
DEALLOCATE(fakeAUXNAMES)
!
!Write CFG file with norm as auxiliary property
pcfg = TRIM(outputfile)//'_edm_all.cfg'
ALLOCATE(CFGEDMNAME(4))
 CFGEDMNAME(1) = 'edm X (e.A)'
 CFGEDMNAME(2) = 'edm Y (e.A)'
 CFGEDMNAME(3) = 'edm Z (e.A)'
 CFGEDMNAME(4) = 'edm norm (e.A)'
ALLOCATE( CFGEDMAUX(SIZE(P(:,1)),4) )
 CFGEDMAUX(:,1:3) = pola(:,1:3)
 CFGEDMAUX(:,4) = normpola(:)
IF(.NOT.overw) CALL CHECKFILE(pcfg,'writ')
CALL WRITE_CFG(H,P,comment,CFGEDMNAME,CFGEDMAUX,pcfg)
DEALLOCATE(CFGEDMNAME)
DEALLOCATE(CFGEDMAUX)
!
! "normpola" contains norms of EDM for all atoms
!Compute statistics and output it
pstat = TRIM(outputfile)//'_edm_stat.txt'
CALL DO_STATS(normpola,mi,M,A,D,T)
IF(.NOT.overw) CALL CHECKFILE(pstat,'writ')
OPEN(UNIT=40,FILE=pstat)
WRITE(40,*) 'Values in e.A (1.602 10^-29 C.m)'
WRITE(40,*) ''
WRITE(40,*) '        Min. value: ', mi
WRITE(40,*) '        Max. value: ', M
WRITE(40,*) '           Average: ', A
WRITE(40,*) 'Av. abs. deviation: ', D
WRITE(40,*) 'Standard deviation: ', T
WRITE(40,*) ''
mi = SUM( pola(:,1) )
WRITE(40,*) '       Sum along X: ', mi
mi = SUM( pola(:,2) )
WRITE(40,*) '       Sum along Y: ', mi
mi = SUM( pola(:,3) )
WRITE(40,*) '       Sum along Z: ', mi
mi = SUM( normpola(:) )
WRITE(40,*) '  Sum of all norms: ', mi
CLOSE(40)
CALL ATOMSK_MSG(4037,(/pstat/),(/0.d0/))
!
!
!
1000 CONTINUE
IF(ALLOCATED(P2)) DEALLOCATE(P2)
IF(ALLOCATED(pola)) DEALLOCATE(pola)
IF(ALLOCATED(normpola)) DEALLOCATE(normpola)
!
!
!
END SUBROUTINE E_DIPOLES
!
!
END MODULE