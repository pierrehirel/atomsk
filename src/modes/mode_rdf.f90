MODULE mode_rdf
!
!**********************************************************************************
!*  MODE_RDF                                                                      *
!**********************************************************************************
!* This module computes the radial distribution function (RDF)                    *
!* of a given system. The number of neighbors in a "skin" of radius R and         *
!* of width dR is computed for each atom in the system, then the average          *
!* number of neighbors N is computed. The RDF is then defined as:                 *
!*    RDF(R) = N(R,R+dR) / dR                                                     *
!* When several atom species are present, the partial RDFs are computed.          *
!* E.g. if atoms A and B are present, the distribution of A atoms around          *
!* A atoms will be computed, then B around A, and B around B.                     *
!* The results are output in special files.                                       *
!**********************************************************************************
!* (C) May. 2012 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 19 Dec. 2023                                     *
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
!Load modules
USE atoms
USE comv
USE constants
USE functions
USE messages
USE neighbors
USE files
USE subroutines
!Module for reading input files
USE readin
!Module for applying options
USE options
!Modules for computing things
USE cmpt_rdf
!
!
CONTAINS
!
SUBROUTINE RDF_XYZ(listfile,rdf_maxR,rdf_dr,options_array)
!
!Declare variables
IMPLICIT NONE
!Input
CHARACTER(LEN=*),INTENT(IN):: listfile    !file containing the names of files to analyze
REAL(dp),INTENT(IN):: rdf_maxR            !maximum radius for RDF
REAL(dp),INTENT(IN):: rdf_dr              !width of the "skin"
!
CHARACTER(LEN=2):: sp1, sp2               !species of atoms type 1, type 2
CHARACTER(LEN=5):: infileformat
CHARACTER(LEN=128):: rdfdat               !Output file names
CHARACTER(LEN=4096):: msg
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties (not used)
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=4096),DIMENSION(:),ALLOCATABLE:: inputfiles    !names of files to analyze
LOGICAL:: fileexists !does the file exist?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: atompair
INTEGER:: i, id, j, k, l, m, n
INTEGER:: u, umin, umax, v, vmin, vmax, w, wmin, wmax
INTEGER:: fid, Nfiles !file index, number of files analyzed
INTEGER:: progress   !To show calculation progress
INTEGER,DIMENSION(:,:),ALLOCATABLE:: NeighList  !list of neighbours
REAL(dp):: nl, nt, ntot !to count atoms of each type
REAL(dp),DIMENSION(3):: Vector
REAL(dp),DIMENSION(3,3):: Huc    !Base vectors of unit cell (unknown, set to 0 here)
REAL(dp),DIMENSION(3,3):: H      !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystal orientation
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX          !auxiliary properties of atoms (not used)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P            !atom positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Pref         !"reference" atom positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S            !shell positions (not used)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: pairs        !pairs of species of 1st system, current system
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: rdf_temp     !temporary partial RDFs for 1 system
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: rdf_partial  !final values of the partial RDFs (space- and time-averaged)
REAL(dp),DIMENSION(:),ALLOCATABLE:: rdf_total      !final values of the total RDF (space- and time-averaged)
!
msg = "ENTERING RDF_XYZ..."
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
Huc(:,:) = 0.d0
 C_tensor(:,:) = 0.d0
!
CALL ATOMSK_MSG(4051,(/""/),(/rdf_maxR,rdf_dr/))
!
!Initialize variables
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
IF(ALLOCATED(inputfiles)) DEALLOCATE(inputfiles)
atompair=0
u=0
v=0
w=0
ORIENT(:,:) = 0.d0
IF(ALLOCATED(aentries)) DEALLOCATE(aentries)
IF(ALLOCATED(Pref)) DEALLOCATE(Pref)
IF(ALLOCATED(pairs)) DEALLOCATE(pairs)
IF(ALLOCATED(rdf_temp)) DEALLOCATE(rdf_temp)        !partial RDFs for 1 system
IF(ALLOCATED(rdf_partial)) DEALLOCATE(rdf_partial)  !partial RDFs
IF(ALLOCATED(rdf_total)) DEALLOCATE(rdf_total)      !total RDF
!
!
!
100 CONTINUE
!Check that input file exists
CALL CHECKFILE(listfile,'read')
!
!Try to determine the format of input file
CALL GUESS_FORMAT(listfile,infileformat,"read")
IF( infileformat=="xxx" ) THEN
  !Unknown file format => assume that it is a list file
  !Open file containing the names of all files to analyze
  OPEN(UNIT=50,FILE=listfile,STATUS='OLD',FORM='FORMATTED')
  REWIND(50)
  !Count how many file names are in the list
  i=0
  DO
    READ(50,*,ERR=110,END=110) msg
    msg = TRIM(ADJUSTL(msg))
    IF(msg(1:1).NE.'#') i=i+1
  ENDDO
  110 CONTINUE
  IF( i<=0 ) GOTO 1000
  ALLOCATE(inputfiles(i))
  inputfiles(:) = ""
  REWIND(50)
  i=0
  DO
    READ(50,*,ERR=120,END=120) msg
    msg = TRIM(ADJUSTL(msg))
    IF(msg(1:1).NE.'#') THEN
      i=i+1
      inputfiles(i) = TRIM(msg)
    ENDIF
  ENDDO
  120 CONTINUE
  CLOSE(50)
ELSE
  !File format was recognized => we have only this file to analyze
  ALLOCATE(inputfiles(1))
  inputfiles(1) = listfile
ENDIF
!
!
!
200 CONTINUE
!Compute the RDFs, space-average for each file, and time-average over all files
!Note: if several atom species exist the partial RDFs must be computed,
!     i.e. for all couples of species (k,l) with l>=k.
Nfiles=0
DO fid=1,SIZE(inputfiles) !Loop on files
  !Initialization
  H(:,:) = 0.d0
  !
  !Check if file actually exists
  INQUIRE(FILE=inputfiles(fid),EXIST=fileexists)
  !
  IF(fileexists) THEN
    !Read data: cell size, atom positions...
    CALL READ_AFF(inputfiles(fid),H,P,S,comment,AUXNAMES,AUX)
    !
    !We won't use shell positions nor auxiliary properties: free memory
    IF(ALLOCATED(S)) DEALLOCATE(S)
    IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
    IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
    !
    !Apply options to the system
    CALL OPTIONS_AFF(options_array,Huc,H,P,S,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
    IF(nerr>0) GOTO 1000
    !
    !If this is the first file, then save atom positions as "reference"
    IF( Nfiles==0 ) THEN
      ALLOCATE( Pref(SIZE(P,1),SIZE(P,2)) )
      Pref(:,:) = P(:,:)
    ELSE
      !Otherwise, compute the difference in atom positions
      m=0
      DO i=1,SIZE(Pref,1)
        nl = VECLENGTH( P(i,2:4)-Pref(i,2:4) )
        IF( nl>rdf_dr ) THEN
          m=m+1
          EXIT
        ENDIF
      ENDDO
      !If max.disp. is smaller than skin dr, keep neighbour list and use it for next snapshot
      IF( m.NE.0 ) THEN
        !Otherwise, if a displacement exceeds threshold value, destroy neighbour list
        !This will force its reconstruction by "COMPUTE_RDF"
        IF(ALLOCATED(NeighList)) DEALLOCATE(NeighList)
        !Current system also becomes "reference"
        Pref(:,:) = P(:,:)
      ENDIF
    ENDIF
    !
    !
    !Construct partial RDFs for current system
    !NOTE: for the first system, array "pairs" is unallocated and will be built by COMPUTE_RDF
    !For the following systems, "pairs" is provided as an input array to COMPUTE_RDF
    CALL COMPUTE_RDF(H,P,rdf_maxR,rdf_dr,NeighList,aentries,pairs,rdf_temp)
    !
    !We are done with this system: free memory
    IF(ALLOCATED(P)) DEALLOCATE(P)
    !
    !If this is the first file, then allocate array "rdf_partial"
    IF( Nfiles==0 ) THEN
      IF(ALLOCATED(rdf_partial)) DEALLOCATE(rdf_partial)
      ALLOCATE( rdf_partial(SIZE(rdf_temp,1),SIZE(rdf_temp,2)) )
      rdf_partial(:,:) = 0.d0
    ENDIF
    !
    !At this point, rdf_temp(:,:) contains the partial RDFs of current system
    !Add partial RDFs of current system into rdf_partial
    !(NOTE: it will be time-averaged later, see label 300)
    rdf_partial(:,:) = rdf_partial(:,:) + rdf_temp(:,:)
    !
    !We don't need rdf_temp any more
    IF(ALLOCATED(rdf_temp)) DEALLOCATE(rdf_temp)
    !
    !File was successfully analyzed
    Nfiles = Nfiles+1
    !
  ELSE
    !Input file doesn't exist: output a warning and go to the next file
    nwarn=nwarn+1
    CALL ATOMSK_MSG(4700,(/TRIM(inputfiles(fid))/),(/0.d0/))
    !
  ENDIF   !end if fileexists
  !
  290 CONTINUE
  !
ENDDO  !loop on m files
!
!
!
300 CONTINUE
CLOSE(50)
IF(ALLOCATED(inputfiles)) DEALLOCATE(inputfiles)
IF(ALLOCATED(NeighList)) DEALLOCATE(NeighList)
IF(ALLOCATED(Pref)) DEALLOCATE(Pref)
!
IF(Nfiles<=0) THEN
  !no file was analyzed => exit
  nerr=nerr+1
  GOTO 1000
ENDIF
!
msg = "Loop on files terminated, averaging RDFs over time..."
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!rdf_partial contains the partial RDFs for each atom pair, summed over all analyzed files:
!rdf_partial( atom_pair , R , g(R) )
!=> divide g(R) by the number of files to make the time-average
rdf_partial(:,:) = rdf_partial(:,:) / DBLE(Nfiles)
!
!Compute the total RDF
!= sum of all partial RDFs weighted with number of atoms
ALLOCATE( rdf_total(SIZE(rdf_partial,2)) )
rdf_total(:) = 0.d0
ntot = 0.d0
DO i=1,SIZE(rdf_partial,1)
  nt = 1.d0
  nl = 0.d0
  !Loop over the two atom types of the pair
  DO k=1,2
    !Search for atoms of type pairs(i,k) in the table aentries
    DO j=1,SIZE(aentries,1)
      IF( DABS(aentries(j,1)-pairs(i,k)) < 1.d-3 ) THEN
        !Found a match
        IF(k==1) THEN
          !First atom of the pair: save number of atoms of type k=1
          nt = aentries(j,2)
          nl = aentries(j,1)
        ELSE !i.e. if k==2
          !Second atom of the pair
          nt = nt * aentries(j,2)
          !If species are different, account for factor 2
          IF( DABS(aentries(j,1)-nl)>0.1d0 ) nt = 2.d0*nt
        ENDIF
        GOTO 310
      ENDIF
    ENDDO
    310 CONTINUE
  ENDDO
  rdf_total(:) = rdf_total(:) + nt*rdf_partial(i,:)
  ntot = ntot + nt
ENDDO
!Normalize total RDF
rdf_total(:) = rdf_total(:) / ntot
!
!
!
400 CONTINUE
!Knowing the RDF, compute the structure factor S
!=Fourier transform of the total correlation function

!
!
!
500 CONTINUE
!Output each partial RDF to files
IF(SIZE(aentries,1)>1) THEN
  atompair=0
  DO k=1,SIZE(aentries,1)
    DO l=k,SIZE(aentries,1)
      atompair=atompair+1
      !
      CALL ATOMSPECIES(aentries(k,1),sp1)
      CALL ATOMSPECIES(aentries(l,1),sp2)
      rdfdat = 'rdf_'//TRIM(sp1)//TRIM(sp2)//'.dat'
      IF(.NOT.overw) CALL CHECKFILE(rdfdat,'writ')
      OPEN(UNIT=40,FILE=rdfdat,FORM="FORMATTED",STATUS="UNKNOWN")
      msg = "# "//TRIM(sp1)//"-"//TRIM(sp2)//" partial RDF computed with Atomsk"
      WRITE(40,'(a)') TRIM(ADJUSTL(msg))
      !
      DO i=1,SIZE(rdf_partial,2)-1
        WRITE(40,'(2f24.8)') DBLE(i-1)*rdf_dr, rdf_partial(atompair,i)
      ENDDO
      !
      CLOSE(40)
      CALL ATOMSK_MSG(4039,(/rdfdat/),(/0.d0/))
    ENDDO
  ENDDO
ENDIF
!
!Output total RDF
rdfdat = 'rdf_total.dat'
IF(.NOT.overw) CALL CHECKFILE(rdfdat,'writ')
OPEN(UNIT=40,FILE=rdfdat,FORM="FORMATTED",STATUS="UNKNOWN")
msg = "# Total RDF computed with Atomsk"
WRITE(40,'(a)') TRIM(ADJUSTL(msg))
DO i=1,SIZE(rdf_total,1)-1
  WRITE(40,'(2f24.8)') DBLE(i-1)*rdf_dr, rdf_total(i)
ENDDO
CLOSE(40)
CALL ATOMSK_MSG(4039,(/rdfdat/),(/0.d0/))
!
!
!
CALL ATOMSK_MSG(4053,(/""/),(/DBLE(Nfiles)/))
GOTO 1000
!
!
!
1000 CONTINUE
IF(ALLOCATED(rdf_partial)) DEALLOCATE(rdf_partial)
IF(ALLOCATED(rdf_total)) DEALLOCATE(rdf_total)
!
!
END SUBROUTINE RDF_XYZ
!
!
END MODULE mode_rdf
