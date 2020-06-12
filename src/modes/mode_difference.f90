MODULE mode_difference
!
!**********************************************************************************
!*  MODE_DIFFERENCE                                                               *
!**********************************************************************************
!* This module reads two arrays Pfirst and Psecond containing atom positions      *
!* and computes the difference between them, taking Pfirst as the reference:      *
!*     difference = Psecond - Pfirst                                              *
!* Then it outputs the displacement vectors and some statistics to files.         *
!* The two arrays Pfirst and Psecond must have the same size. They should         *
!* also contain similar systems (e.g. two snaphots from MD) for the               *
!* computation to make sense, although this is not checked here.                  *
!* If auxiliary properties exist in both system, the difference between matching  *
!* auxiliary properties is also computed. Auxiliary properties that exist in      *
!* one system but not both are ignored and will not appear in output files.       *
!**********************************************************************************
!* (C) March 2010 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 25 May 2020                                      *
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
!Module for reading input files
USE readin
!Module for applying options
USE options
!Output modules
USE out_cfg
USE out_xyz
USE out_xsf
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
CHARACTER(LEN=128):: diffcfgfile, diffxyzfile, diffnormfile, diffxsffile
CHARACTER(LEN=128):: bothxsffile, histdatfile, stattxtfile
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES, AUXNAMES1, AUXNAMES2 !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: i, j, k, l
INTEGER:: histcount, histcumul !counters for histogram
INTEGER:: Naux  !number of auxiliary properties that exist in both systems
INTEGER,DIMENSION(100,2):: tabAUX !correspondance table for aux.prop. (assuming <100 properties)
REAL(dp):: histstep  !step for histogram in Angstroms
REAL(dp):: histdown, histup !boundaries for histogram
REAL(dp):: stat_min, stat_M, stat_A, stat_D, stat_S !statistics on displacements
REAL(dp),DIMENSION(3,3):: Huc   !Box vectors of unit cell (unknown, set to 0 here)
REAL(dp),DIMENSION(3,3):: H, H1, H2   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT     !crystallographic orientation of the system
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:),ALLOCATABLE:: spi_table !table containing atomic numbers
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX, AUX1, AUX2  !auxiliary properties of system 1 and 2
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: atypes
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Pfirst, Psecond !atom positions of system 1 and 2
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S   !shell positions (ignored here)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Pboth !atom positions of both systems
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
 comment(1) = '# Differences in atomic positions'
!Set file names
diffcfgfile = TRIM(ADJUSTL(outputfile))//'_diff.cfg'
diffxyzfile = TRIM(ADJUSTL(outputfile))//'_diff.xyz'
diffnormfile = TRIM(ADJUSTL(outputfile))//'_norm.dat'
diffxsffile = TRIM(ADJUSTL(outputfile))//'_diff.xsf'
bothxsffile = TRIM(ADJUSTL(outputfile))//'_both.xsf'
histdatfile = TRIM(ADJUSTL(outputfile))//'_hist.dat'
stattxtfile = TRIM(ADJUSTL(outputfile))//'_stat.txt'
!
!
CALL ATOMSK_MSG(4038,(/''/),(/0.d0/))
!
!
100 CONTINUE
!Read atomic positions from filefirst and store them into Pfirst(:,:)
CALL READ_AFF(filefirst,H1,Pfirst,S,comment,AUXNAMES1,AUX1)
!Remove shells
IF(ALLOCATED(S)) DEALLOCATE(S)
!Apply options to system 1
CALL OPTIONS_AFF(options_array,Huc,H1,Pfirst,S,AUXNAMES1,AUX1,ORIENT,SELECT,C_tensor)
!
!Read atomic positions from filesecond and store them into Psecond(:,:)
CALL READ_AFF(filesecond,H2,Psecond,S,comment,AUXNAMES2,AUX2)
!Remove shells
IF(ALLOCATED(S)) DEALLOCATE(S)
!Apply options to system 2
CALL OPTIONS_AFF(options_array,Huc,H2,Psecond,S,AUXNAMES2,AUX2,ORIENT,SELECT,C_tensor)
!
!Verify that both systems contain the same number of atoms
IF( SIZE(Pfirst,1) .NE. SIZE(Psecond,1) ) THEN
  CALL ATOMSK_MSG(4810,(/''/), &
       & (/ DBLE(SIZE(Pfirst(:,1))), DBLE(SIZE(Psecond(:,1))) /))
  nerr = nerr+1
  GOTO 1000
ENDIF
!
!Construct a correspondance table between auxiliary properties of systems 1 and 2
Naux=0
IF( ALLOCATED(AUXNAMES1) .AND. ALLOCATED(AUXNAMES2) .AND.                  &
  & ALLOCATED(AUX1) .AND. ALLOCATED(AUX2) .AND. SIZE(AUX1,1)==SIZE(AUX2,1) ) THEN
  !
  DO i=1,SIZE(AUXNAMES1)
    DO j=1,SIZE(AUXNAMES2)
      IF( AUXNAMES1(i)==AUXNAMES2(j) ) THEN
        !Auxiliary property #i in AUX1 has same name as au.prop. #j in AUX2
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
!Compute differences and save them in AUX
ALLOCATE(AUXNAMES(Naux+4))
DO j=1,Naux
  AUXNAMES(j) = "diff_"//TRIM(ADJUSTL(AUXNAMES1(tabAUX(j,1))))
ENDDO
AUXNAMES(Naux+1) = 'dx'
AUXNAMES(Naux+2) = 'dy'
AUXNAMES(Naux+3) = 'dz'
AUXNAMES(Naux+4) = 'dtot'
ALLOCATE(AUX(SIZE(Pfirst,1),Naux+4))
AUX(:,:) = 0.d0
DO i=1,SIZE(Pfirst,1)
  IF( Naux>0 ) THEN
    !Compute differences in auxiliary properties
    DO j=1,Naux
      AUX(i,j) = AUX2(i,tabAUX(j,2)) - AUX1(i,tabAUX(j,1))
    ENDDO
  ENDIF
  !Compute displacements
  !Note: unwrapped cartesian coordinates are assumed here
  AUX(i,Naux+1) = Psecond(i,1) - Pfirst(i,1)
  AUX(i,Naux+2) = Psecond(i,2) - Pfirst(i,2)
  AUX(i,Naux+3) = Psecond(i,3) - Pfirst(i,3)
  AUX(i,Naux+4) = VECLENGTH(AUX(i,Naux+1:Naux+3))
ENDDO
!
!
!
400 CONTINUE
!Check if files already exist
IF(.NOT.overw) CALL CHECKFILE(diffxyzfile,'writ')
IF(.NOT.overw) CALL CHECKFILE(diffnormfile,'writ')
IF(.NOT.overw) CALL CHECKFILE(diffxsffile,'writ')
IF(.NOT.overw) CALL CHECKFILE(bothxsffile,'writ')
IF(.NOT.overw) CALL CHECKFILE(stattxtfile,'writ')
!
!Write files
!
!XYZ file containing atom positions of 1st system + displacements
CALL WRITE_XYZ(H1,Pfirst,comment,AUXNAMES,AUX,diffxyzfile,'exyz ')
!
!CFG file containing atom positions of 1st system + displacements
CALL WRITE_CFG(H1,Pfirst,comment,AUXNAMES,AUX,diffcfgfile)
!
!XSF file containing both systems
ALLOCATE( Pboth( SIZE(Pfirst,1)+SIZE(Psecond,1), 4 ) )
DO i=1,SIZE(Pfirst,1)
  Pboth(i,:) = Pfirst(i,:)
ENDDO
DO i=1,SIZE(Psecond,1)
  Pboth(SIZE(Pfirst,1)+i,:) = Psecond(i,:)
ENDDO
CALL WRITE_XSF(H1,Pboth,comment,AUXNAMES,AUX,bothxsffile)
DEALLOCATE(Pboth)
!
!XSF file containing atom positions of 1st system + displacements
!Change name to trick module into writing displacements to XSF
AUXNAMES(Naux+1) = 'fx'
AUXNAMES(Naux+2) = 'fy'
AUXNAMES(Naux+3) = 'fz'
CALL WRITE_XSF(H1,Pfirst,comment,AUXNAMES,AUX,diffxsffile)
!
!Text file containing the norms of displacements
ALLOCATE(stattable(SIZE(Pfirst,1),2))
OPEN(UNIT=36,FILE=diffnormfile,FORM='FORMATTED',STATUS='UNKNOWN')
DO i=1,SIZE(Pfirst,1)
  CALL ATOMSPECIES(Pfirst(i,4),species)
  WRITE(36,416) i, AUX(i,Naux+4), species
  !
  !Store norm of displacements in a table
  !for later statistics calculations (see label 450)
  stattable(i,1) = Pfirst(i,4)
  stattable(i,2) = AUX(i,Naux+4)
ENDDO
416 FORMAT(i6,1X,f12.6,1X,a3)
CLOSE(36)
CALL ATOMSK_MSG(4039,(/TRIM(diffnormfile)/),(/0.d0/))
!
!
430 CONTINUE
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
!
!
440 CONTINUE
!Calculate the statistics of displacements
OPEN(UNIT=39,FILE=stattxtfile,FORM='FORMATTED',STATUS='UNKNOWN')
WRITE(39,*) 'STATISTICS ON THE DISPLACEMENTS'
WRITE(39,*) 'All values are in angstroms'
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
CALL FIND_NSP(Pfirst(:,4),atypes)
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
  DO l=1,SIZE(Pfirst,1)
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
