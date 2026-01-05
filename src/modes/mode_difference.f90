MODULE mode_difference
!
!**********************************************************************************
!*  MODE_DIFFERENCE                                                               *
!**********************************************************************************
!* This module reads two arrays P1 and P2 containing atom positions               *
!* and computes the difference between them, taking P1 as the reference:          *
!*     difference = P2 - P1                                                       *
!* Then it outputs the displacement vectors and some statistics to files.         *
!* If the two arrays P1 and P2 have the same size, then it is assumed that        *
!* they contain atoms listed in the same order.                                   *
!* If the two arrays P1 and P2 have different sizes, then an attempt is made to   *
!* match each atom of P1 with an atom in P2. If a match is found, differences     *
!* (i.e. displacements, aux.prop.) are computed. If no match is found, then       *
!* atom is marked as "new" (i.e. it exists in P1 but not in P2).                  *
!* Systems P1 and P2 *should* contain similar systems (e.g. two snaphots from MD) *
!* for the computation to make sense, although this is not checked here.          *
!* If auxiliary properties exist in both system, the difference between matching  *
!* auxiliary properties is also computed. Auxiliary properties that exist in      *
!* one system but not both are ignored and will not appear in output files.       *
!**********************************************************************************
!* (C) March 2010 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 18 Nov. 2025                                     *
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
USE sorting
!Module for reading input files
USE readin
!Module for applying options
USE options
!Output modules
USE writeout
USE out_xyz
!
!
CONTAINS
!
!
SUBROUTINE DIFF_XYZ(filefirst,filesecond,options_array,outputfile)
!
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: filefirst, filesecond
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: outputfile, msg
CHARACTER(LEN=128):: diffcfgfile
CHARACTER(LEN=128):: histdatfile, stattxtfile
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to write
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES, AUXNAMES1, AUXNAMES2 !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment, comment2
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: i, j, k, l
INTEGER:: histcount, histcumul !counters for histogram
INTEGER:: Naux  !number of auxiliary properties that exist in both systems
INTEGER:: status
INTEGER,DIMENSION(:),ALLOCATABLE:: datnum  !difference in atomic number
INTEGER,DIMENSION(:),ALLOCATABLE:: paired  !index of atom with which atom #i is paired
INTEGER,DIMENSION(100,2):: tabAUX !correspondance table for aux.prop. (assuming <100 properties)
REAL(dp):: histstep  !step for histogram in Angstroms
REAL(dp):: histdown, histup !boundaries for histogram
REAL(dp):: stat_min, stat_M, stat_A, stat_D, stat_S !statistics on displacements
REAL(dp),DIMENSION(3,3):: Huc        !Box vectors of unit cell (unknown, set to 0 here)
REAL(dp),DIMENSION(3,3):: H1, H2     !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT     !crystallographic orientation of the system
REAL(dp),DIMENSION(9,9):: C_tensor   !elastic tensor
REAL(dp),DIMENSION(:),ALLOCATABLE:: spi_table !table containing atomic numbers
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX, AUX1, AUX2  !auxiliary properties of system 1 and 2
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: atypes
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P1, P2 !atom positions of system 1 and 2
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S      !shell positions (ignored here)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Pboth  !atom positions of both systems
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: stattable !table containing norms of displacements
!
!
!Initialize variables
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
Naux=0
tabAUX(:,:) = 0
histstep = 1.d0
Huc(:,:) = 0.d0
 C_tensor(:,:) = 0.d0
ORIENT(:,:) = 0.d0
IF(ALLOCATED(Pboth)) DEALLOCATE(Pboth)
ALLOCATE(comment(1))
 comment(1) = "# Differences between '"//TRIM(filefirst)//"' and '"//TRIM(filesecond)//"'"
!Set file names
diffcfgfile = TRIM(ADJUSTL(outputfile))//'_diff.cfg'
histdatfile = TRIM(ADJUSTL(outputfile))//'_hist.dat'
stattxtfile = TRIM(ADJUSTL(outputfile))//'_stat.txt'
!
!
CALL ATOMSK_MSG(4038,(/''/),(/0.d0/))
!
!
100 CONTINUE
!Read atomic positions from filefirst and store them into P1(:,:)
CALL READ_AFF(filefirst,H1,P1,S,comment2,AUXNAMES1,AUX1)
IF(nerr>0) GOTO 1000
!Remove shells
IF(ALLOCATED(S)) DEALLOCATE(S)
!Apply options to system 1
CALL OPTIONS_AFF(options_array,Huc,H1,P1,S,AUXNAMES1,AUX1,ORIENT,SELECT,C_tensor)
!
!Read atomic positions from filesecond and store them into P2(:,:)
CALL READ_AFF(filesecond,H2,P2,S,comment2,AUXNAMES2,AUX2)
IF(nerr>0) GOTO 1000
!Remove shells
IF(ALLOCATED(S)) DEALLOCATE(S)
!Apply options to system 2
CALL OPTIONS_AFF(options_array,Huc,H2,P2,S,AUXNAMES2,AUX2,ORIENT,SELECT,C_tensor)
!
!Construct a correspondance table between auxiliary properties of systems 1 and 2.
!Naux will be the number of aux.prop. that exist in both systems,
!and for which we can compute a difference.
!NOTE: aux.prop. that exist in only one of the two systems will be ignored
Naux=0
IF( ALLOCATED(AUXNAMES1) .AND. ALLOCATED(AUXNAMES2) .AND.    &
  & ALLOCATED(AUX1) .AND. ALLOCATED(AUX2) .AND.              &
  & SIZE(AUX1,1)==SIZE(P1,1) .AND. SIZE(AUX2,1)==SIZE(P2,1) ) THEN
  !
  DO i=1,SIZE(AUXNAMES1)
    DO j=1,SIZE(AUXNAMES2)
      IF( AUXNAMES1(i)==AUXNAMES2(j) ) THEN
        !Auxiliary property #i in AUX1 has same name as aux.prop. #j in AUX2
        Naux = Naux+1
        tabAUX(Naux,1) = i
        tabAUX(Naux,2) = j
      ENDIF
    ENDDO
  ENDDO
  !
ENDIF
!
!
!
300 CONTINUE
!Verify if both systems contain the same number of atoms
IF( SIZE(P1,1) == SIZE(P2,1) ) THEN
  !Both systems contain the same number of atoms
  !IMPORTANT NOTE: it is assumed here that atoms are sorted in the same order in both arrays
  !
  !Prepare array AUX to store differences
  ALLOCATE(AUXNAMES(Naux+4))
  IF( Naux>0 ) THEN
    DO j=1,Naux
      AUXNAMES(j) = "diff_"//TRIM(ADJUSTL(AUXNAMES1(tabAUX(j,1))))
    ENDDO
  ENDIF
  AUXNAMES(Naux+1) = 'dx'
  AUXNAMES(Naux+2) = 'dy'
  AUXNAMES(Naux+3) = 'dz'
  AUXNAMES(Naux+4) = 'dtot'
  ALLOCATE(AUX(SIZE(P1,1),Naux+4))
  AUX(:,:) = 0.d0
  ALLOCATE( datnum(SIZE(P1,1)) )
  datnum(:) = 0
  !
  !Loop on all atoms
  CALL ATOMSK_MSG(4038,(/''/),(/0.d0/))
  DO i=1,SIZE(P1,1)
    IF( Naux>0 ) THEN
      !Compute differences in auxiliary properties
      DO j=1,Naux
        AUX(i,j) = AUX2(i,tabAUX(j,2)) - AUX1(i,tabAUX(j,1))
      ENDDO
    ENDIF
    !Compute displacements
    !Note: unwrapped cartesian coordinates are assumed here
    AUX(i,Naux+1) = P2(i,1) - P1(i,1)
    AUX(i,Naux+2) = P2(i,2) - P1(i,2)
    AUX(i,Naux+3) = P2(i,3) - P1(i,3)
    AUX(i,Naux+4) = VECLENGTH(AUX(i,Naux+1:Naux+3))
    !Compute difference in atomic number
    datnum = NINT( P2(i,4) - P1(i,4) )
  ENDDO
  !
  !We won't need initial auxiliary properties any more
  IF(ALLOCATED(AUXNAMES1)) DEALLOCATE(AUXNAMES1)
  IF(ALLOCATED(AUXNAMES2)) DEALLOCATE(AUXNAMES2)
  IF(ALLOCATED(AUX1)) DEALLOCATE(AUX1)
  IF(ALLOCATED(AUX2)) DEALLOCATE(AUX2)
  !
  !If there is any difference in atomic number, add it as a new auxiliary property
  IF( ANY(datnum.NE.0) ) THEN
    Naux = Naux+1
    ALLOCATE( AUXNAMES1(Naux) )
    IF( Naux>0 ) THEN
      DO j=1,Naux-1
        AUXNAMES1(j) = AUXNAMES(j)
      ENDDO
    ENDIF
    AUXNAMES1(Naux) = "diff_Z"
    DEALLOCATE(AUXNAMES)
    ALLOCATE(AUXNAMES(SIZE(AUXNAMES1)))
    AUXNAMES(:) = AUXNAMES1(:)
    DEALLOCATE(AUXNAMES)
    CALL RESIZE_DBLEARRAY2(AUX,SIZE(AUX,1),Naux,status)
    IF(status>0) THEN
      CALL ATOMSK_MSG(818,(/"AUX"/),(/0.d0/))
      nerr = nerr+1
      GOTO 1000
    ENDIF
    DO j=1,SIZE(AUX,1)
      AUX(j,Naux) = DBLE(datnum(j))
    ENDDO
  ENDIF
  !
ELSE
  !The two systems DO NOT have same number of atoms
  !Attempt to match each atom in P1 with an atom in P2,
  !and to compute differences (P2 - P1).
  !If P1 is smaller, then each atom should have an equivalent in P2
  !Otherwise (i.e.P1 larger), atoms that don't have an equivalent will be
  !declared as "deleted"
  !
  !Assign pointers to smallest and largest system
  IF( SIZE(P1,1) < SIZE(P2,1) ) THEN
    l = 4
  ELSE
    l = 5   !5th column to store "deleted" atoms value
  ENDIF
  !
  !Prepare array AUX to store differences
  ALLOCATE(AUXNAMES(Naux+l))
  DO j=1,Naux
    AUXNAMES(j) = "diff_"//TRIM(ADJUSTL(AUXNAMES1(tabAUX(j,1))))
  ENDDO
  AUXNAMES(Naux+1) = 'dx'
  AUXNAMES(Naux+2) = 'dy'
  AUXNAMES(Naux+3) = 'dz'
  AUXNAMES(Naux+4) = 'dtot'
  IF(l==5) AUXNAMES(Naux+5) = 'deleted_atoms'
  ALLOCATE(AUX(SIZE(P1,1),Naux+l))
  AUX(:,:) = 0.d0
  !
  !Find matching indices between P1 and P2
  CALL ATOMSK_MSG(4077,(/""/),(/0.d0/))
  CALL FIND_MATCHING_ID(P1,P2,paired,j)
  !
  IF( j==0 ) THEN
    !No equivalent pairs of atoms could be found => abort
    CALL ATOMSK_MSG(4835,(/''/),(/0.d0/))
    nerr = nerr+1
    GOTO 1000
  ENDIF
  !
  CALL ATOMSK_MSG(4078,(/""/),(/DBLE(j),DBLE(SIZE(P1,1)-j)/))
  !NOTE: since P1 and P2 have different sizes, paired(:) may contain zeros
  !
  !Compute differences
  CALL ATOMSK_MSG(4038,(/''/),(/0.d0/))
  DO i=1,SIZE(paired)
    j = paired(i)
    !Check that atom #i was paired
    IF( j<=0 ) THEN
      !Atom #i was NOT paired: can not compute differences, mark it "deleted"
      IF( l==5 ) AUX(i,Naux+5) = 1.d0
      !
    ELSEIF( j<=SIZE(P2,1) ) THEN
      !Atom #i in P1 was paired with atom #j in P2
      !Compute displacements
      AUX(i,Naux+1) = P2(j,1) - P1(i,1)
      AUX(i,Naux+2) = P2(j,2) - P1(i,2)
      AUX(i,Naux+3) = P2(j,3) - P1(i,3)
      AUX(i,Naux+4) = VECLENGTH(AUX(i,Naux+1:Naux+3))
      !Compute differences in auxiliary properties, if any
      IF( Naux>0 ) THEN
        DO j=1,Naux
          AUX(i,j) = AUX2(j,tabAUX(j,2)) - AUX1(i,tabAUX(j,1))
        ENDDO
      ENDIF
    ENDIF
  ENDDO !i
  !
  IF(ALLOCATED(paired)) DEALLOCATE(paired)
  !
ENDIF !end if systems have same number of atoms
!
!
!
400 CONTINUE
!We won't need info about second system any more
IF(ALLOCATED(P2)) DEALLOCATE(P2)
IF(ALLOCATED(AUXNAMES2)) DEALLOCATE(AUXNAMES2)
IF(ALLOCATED(AUX2)) DEALLOCATE(AUX2)
!
!Write output files
!
!CFG file containing atom positions of 1st system + displacements
ALLOCATE(outfileformats(1))
outfileformats(1) = "cfg"
CALL WRITE_AFF(diffcfgfile,outfileformats,H1,P1,S,comment,AUXNAMES,AUX)
!
!
!
440 CONTINUE
!Store norm of displacements in a table
!for later statistics calculations (see label 450)
ALLOCATE(stattable(SIZE(P1,1),2))
DO i=1,SIZE(P1,1)
  stattable(i,1) = P1(i,4)
  stattable(i,2) = AUX(i,Naux+4)
ENDDO
!Don't need AUX any more => free memory
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
!
!Write statistics of displacements in a file
IF(.NOT.overw) CALL CHECKFILE(stattxtfile,'writ')
OPEN(UNIT=39,FILE=stattxtfile,FORM='FORMATTED',STATUS='UNKNOWN')
WRITE(39,*) 'STATISTICS ON THE DISPLACEMENTS'
WRITE(39,*) 'All values are in angstroms'
WRITE(39,*) TRIM(comment(1))
WRITE(39,*) ''
WRITE(39,*) 'Sp. = atomic Spieces (TO=Total statistics on all atoms)'
WRITE(39,*) 'm = Minimum displacement'
WRITE(39,*) 'M = Maximum displacement'
WRITE(39,*) 'A = Average displacement'
WRITE(39,*) 'D = Average absolute deviation'
WRITE(39,*) 'S = Standard deviation'
WRITE(39,*) ''
WRITE(39,*) 'Sp. |     m      |     M      |     A      |     D      |     S'
WRITE(39,*) '----+------------+------------+------------+------------+------------'
!
!Statistics on individual species
!Find total number of particles for each species
CALL FIND_NSP(P1(:,4),atypes)
!
!For each species we calculate the statistics on the norm of displacements
DO i=1,SIZE(atypes,1)
  WRITE(msg,*) 'species, NP: ', atypes(i,1), atypes(i,1)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Create a table to store norm of displacements
  ALLOCATE( spi_table(INT(atypes(i,2))) )
  spi_table(:) = 0.d0
  !For all atoms of this species, place norm of displ. in the table
  k = 0
  DO l=1,SIZE(P1,1)
    IF( stattable(l,1) == atypes(i,1) ) THEN
      k = k+1
      spi_table(k) = stattable(l,2)
    ENDIF
  ENDDO
  CALL ATOMSPECIES(DBLE(atypes(i,1)),species)
  !Calculate statistics on the norm of displ.
  CALL DO_STATS( spi_table(:),stat_min,stat_M,stat_A,stat_D,stat_S )
  !Write the stats in the file
  WRITE(39,460) species, '|', stat_min, '|', stat_M, '|', &
       &         stat_A, '|',stat_D, '|', stat_S
  !Don't forget to deallocate the table for next loop
  DEALLOCATE(spi_table)
ENDDO
!
!Total statistics on all atoms
CALL DO_STATS( stattable(:,2),stat_min,stat_M,stat_A,stat_D,stat_S )
WRITE(39,*) '----+------------+------------+------------+------------+------------'
WRITE(39,460) 'TO', '|', stat_min, '|', stat_M, '|', stat_A, '|', stat_D, '|', stat_S
460 FORMAT(1X,a2,2X,a1,1X,5(f10.4,1X,a1,1X))
!
CLOSE(39)
CALL ATOMSK_MSG(4039,(/TRIM(stattxtfile)/),(/0.d0/))
!
!
450 CONTINUE
!Text file containing a histogram:
!Natoms that have a displacement between histdown and histup
!Use a step that will display at least 20 intervals,
!but that is maximum 1 Angströms
histstep = MIN( 1.d0 , stat_M/20.d0 )
!
IF(.NOT.overw) CALL CHECKFILE(histdatfile,'writ')
OPEN(UNIT=35,FILE=histdatfile,FORM='FORMATTED',STATUS='UNKNOWN')
histdown=0.d0
histup=histstep
histcumul=0
DO WHILE( histup < stat_M )
  histcount = 0
  !Parse all norms of displacements
  DO i=1,SIZE(stattable,1)
    IF( stattable(i,2)>=histdown .AND. stattable(i,2)<histup ) THEN
      !If displacement is in the interval, increase the counters
      histcount = histcount+1
      histcumul = histcumul+1
    ENDIF
  ENDDO
  !Write the number of atoms that have this displacement
  WRITE(35,'(f12.3,2i12)') histdown, histcount, histcumul
  !Go to the next interval
  histup = histup+histstep
  histdown = histdown+histstep
ENDDO
CLOSE(35)
CALL ATOMSK_MSG(4039,(/TRIM(histdatfile)/),(/0.d0/))
!
!
!
500 CONTINUE
CALL ATOMSK_MSG(4040,(/''/),(/0.d0/))
GOTO 1000
!
!
!
850 CONTINUE
CALL ATOMSK_MSG(4811,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
1000 CONTINUE
IF(ALLOCATED(stattable)) DEALLOCATE(stattable)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
!
!
END SUBROUTINE DIFF_XYZ
!
!
END MODULE mode_difference
