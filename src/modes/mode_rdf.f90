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
!* Last modification: P. Hirel - 14 June 2022                                     *
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
USE readin
USE options
!
!
CONTAINS
!
SUBROUTINE RDF_XYZ(listfile,rdf_maxR,rdf_dr,options_array)
!
!Declare variables
IMPLICIT NONE
!Input
CHARACTER(LEN=*),INTENT(IN):: listfile  !file containing the names of files to analyze
REAL(dp),INTENT(IN):: rdf_maxR            !maximum radius for RDF
REAL(dp),INTENT(IN):: rdf_dr              !width of the "skin"
!
CHARACTER(LEN=2):: sp1, sp2               !species of atoms type 1, type 2
CHARACTER(LEN=128):: msg
CHARACTER(LEN=128):: rdfdat               !Output file names
CHARACTER(LEN=4096):: inputfile           !name of a file to analyze
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties (not used)
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: fileexists !does the file exist?
LOGICAL:: doanalysis
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: atompair
INTEGER:: i, id, j, k, l, m, n
INTEGER:: u, umin, umax, v, vmin, vmax, w, wmin, wmax
INTEGER:: Nfiles     !number of files analyzed
INTEGER:: Nspecies   !number of different atom species in the system
INTEGER:: progress   !To show calculation progress
INTEGER:: rdf_Nsteps !number of "skins" for RDF
INTEGER(KIND=SELECTED_INT_KIND(15)):: nt, ntot !to count atoms of each type
INTEGER,DIMENSION(:,:),ALLOCATABLE:: NeighList  !list of neighbours
REAL(dp):: distance     !distance between 2 atoms
REAL(dp):: rdf_norm     !normalization factor
REAL(dp):: rdf_radius   !radius of the sphere
REAL(dp):: Vsphere, Vskin        !volume of the sphere, skin
REAL(dp):: Vsystem               !volume of the system
REAL(dp),DIMENSION(3):: Vector
REAL(dp),DIMENSION(3,3):: Huc    !Base vectors of unit cell (unknown, set to 0 here)
REAL(dp),DIMENSION(3,3):: H      !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystal orientation
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries, aentries2
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX          !auxiliary properties of atoms (not used)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P            !atom positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S            !shell positions (not used)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: pairs        !pairs of species
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: rdf_temp     !temporary partial RDFs for 1 system
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: rdf_partial  !final values of the partial RDFs (space- and time-averaged)
REAL(dp),DIMENSION(:),ALLOCATABLE:: rdf_total      !final values of the total RDF (space- and time-averaged)
!
msg = 'ENTERING RDF_XYZ...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
Huc(:,:) = 0.d0
 C_tensor(:,:) = 0.d0
!
CALL ATOMSK_MSG(4051,(/""/),(/rdf_maxR,rdf_dr/))
!
!Initialize variables
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
atompair=0
u=0
v=0
w=0
ORIENT(:,:) = 0.d0
IF(ALLOCATED(pairs)) DEALLOCATE(pairs)
IF(ALLOCATED(rdf_partial)) DEALLOCATE(rdf_temp)     !partial RDFs for 1 system
IF(ALLOCATED(rdf_partial)) DEALLOCATE(rdf_partial)  !partial RDFs
IF(ALLOCATED(rdf_total)) DEALLOCATE(rdf_total)      !total RDF
!
!
!
100 CONTINUE
!Set number of steps
rdf_Nsteps = NINT(rdf_maxR/rdf_dr) + 2
!Open file containing the names of all files to analyze
CALL CHECKFILE(listfile,'read')
OPEN(UNIT=50,FILE=listfile,STATUS='OLD',FORM='FORMATTED')
REWIND(50)
!
!
!
200 CONTINUE
!Compute the RDFs, space-average for each file, and time-average over all files
!Note: if several atom species exist the partial RDFs must be computed,
!     i.e. for all couples of species (k,l) with l>=k.
Nfiles=0
DO  !Loop on files
  !Initialization
  H(:,:) = 0.d0
  Vsystem = 0.d0
  IF(ALLOCATED(NeighList)) DEALLOCATE(NeighList)
  !
  !Read the name of the file containing the system to analyze
  READ(50,'(a128)',END=300,ERR=300) inputfile
  inputfile = ADJUSTL(inputfile)
  !
  IF( LEN_TRIM(inputfile)>0 .AND. inputfile(1:1).NE.'#' ) THEN
    !Check if file actually exists
    INQUIRE(FILE=inputfile,EXIST=fileexists)
    !
    IF(fileexists) THEN
      doanalysis = .TRUE.
      !Read data: cell size, atom positions...
      CALL READ_AFF(inputfile,H,P,S,comment,AUXNAMES,AUX)
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
      !Determine if we have to look for periodic replica of atoms
      !Minimum image convention: look for all replica in the radius rdf_maxR
      umin=-1
      umax=1
      vmin=-1
      vmax=1
      wmin=-1
      wmax=1
      IF ( VECLENGTH(H(1,:))<=rdf_maxR ) THEN
        umax = CEILING( rdf_maxR / VECLENGTH(H(1,:)) ) + 1
        umin = -1*umax
      ENDIF
      IF ( VECLENGTH(H(2,:))<=rdf_maxR ) THEN
        vmax = CEILING( rdf_maxR / VECLENGTH(H(2,:)) ) + 1
        vmin = -1*vmax
      ENDIF
      IF ( VECLENGTH(H(3,:))<=rdf_maxR ) THEN
        wmax = CEILING( rdf_maxR / VECLENGTH(H(3,:)) ) + 1
        wmin = -1*wmax
      ENDIF
      !
      !If it is the first system, initialize arrays and variables
      IF(Nfiles==0) THEN
        !Count how many different species exist in the system
        CALL FIND_NSP(P(:,4),aentries)
        Nspecies = SIZE(aentries,1)
        IF( verbosity>=4 ) THEN
          WRITE(msg,*) "Nspecies = ", Nspecies
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          DO k=1,SIZE(aentries,1)
            WRITE(msg,*) "    #", NINT(aentries(k,1)), ": ", NINT(aentries(k,2)), " atoms"
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          ENDDO
        ENDIF
        !
        !Construct array storing pairs of atoms
        k = Nspecies*(Nspecies+1)/2
        ALLOCATE( pairs(k,2) )
        pairs(:,2) = 0.d0
        m=0
        DO i=1,SIZE(aentries,1)
          DO j=i,SIZE(aentries,1)
            m=m+1
            pairs(m,1) = aentries(i,1)
            pairs(m,2) = aentries(j,1)
          ENDDO
        ENDDO
        IF( verbosity>=4 ) THEN
          WRITE(msg,*) "N. pairs = ", SIZE(pairs,1)
          CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          DO k=1,SIZE(pairs,1)
            WRITE(msg,*) "    #", k, ": ", NINT(pairs(k,1)), NINT(pairs(k,2))
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
          ENDDO
        ENDIF
        !
        !Allocate array to store the partial RDFs
        ALLOCATE( rdf_partial( SIZE(pairs,1) , rdf_Nsteps ) )
        rdf_partial(:,:) = 0.d0
        !
      ELSE
        !This is not the first file
        !Count how many different species exist in the system
        CALL FIND_NSP(P(:,4),aentries2)
        !Compare with the entries of the first system
        m = 0
        DO i=1,SIZE(aentries,1)
          n=0
          DO j=1,SIZE(aentries2,1)
            IF( DABS( aentries(i,1)-aentries2(j,1) ) < 1.d-3 ) THEN
              n=n+1
              m=m+1
              EXIT
            ENDIF
          ENDDO
          IF( n==0 ) THEN
            !Atoms of type aentries(i,1) existed in first system, but not in current system
            nwarn=nwarn+1
            CALL ATOMSPECIES(aentries(i,1),sp1)
            CALL ATOMSK_MSG(4719,(/sp1/),(/DBLE(Nfiles+1)/))
          ENDIF
        ENDDO
        IF( m==0 ) THEN
          !No atom species of the first system was found in current system
          !Something is very wrong (probably the user's fault?)
          !Anyway, no point in analyzing current system
          GOTO 290
        ENDIF
      ENDIF
      !
      !Construct neighbour list
      !NOTE: to compute RDF up to R, we need neighbors up to R+dR
      CALL ATOMSK_MSG(11,(/""/),(/0.d0/))
      CALL NEIGHBOR_LIST(H,P,rdf_maxR+1.1d0*rdf_dr,NeighList)
      CALL ATOMSK_MSG(15,(/""/),(/0.d0/))
      !
      !
      !Construct partial RDFs for current system
      CALL ATOMSK_MSG(4052,(/""/),(/0.d0/))
      IF(ALLOCATED(rdf_temp)) DEALLOCATE(rdf_temp)
      ALLOCATE( rdf_temp( SIZE(rdf_partial,1) , rdf_Nsteps ) )
      rdf_temp(:,:) = 0.d0
      !Loop on all atoms
      progress=0
      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP& PRIVATE(i,id,j,m,n,u,v,w,atompair,distance,Vector) &
      !$OMP& REDUCTION(+:rdf_temp)
      DO i=1,SIZE(P,1)
        !
        IF( SIZE(P,1)>20000 ) THEN
          progress = progress+1
          !If there are many atoms, display a fancy progress bar
          CALL ATOMSK_MSG(10,(/""/),(/DBLE(progress),DBLE(SIZE(P,1))/))
        ENDIF
        !
        !Get index of "atompair" for the pair #i-#i
        DO j=1,SIZE(pairs,1)
          IF( DABS(P(i,4)-pairs(j,1)) < 1.d-3 .AND. &
              & DABS(P(i,4)-pairs(j,2)) < 1.d-3     ) THEN
            atompair = j
            EXIT
          ENDIF
        ENDDO
        !
        !Check if replicas of atom #i are neighbors of atom #i
        !NOTE: exclude self image as it would result in distance=0
        DO u=umin,umax
          DO v=vmin,vmax
            DO w=wmin,wmax
              distance = VECLENGTH( DBLE(u)*H(1,:) + DBLE(v)*H(2,:) + DBLE(w)*H(3,:) )
              IF( DABS(distance)>1.d-12 .AND. distance<(rdf_maxR+rdf_dr) ) THEN
                !Add this replica in the appropriate "atompair" and appropriate skin
                !n = index of skin
                n = MAX(1,CEILING(distance/rdf_dr))
                rdf_temp(atompair,n) = rdf_temp(atompair,n) + 1
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        !
        !Parse the neighbour list of atom #i
        m=1
        DO WHILE ( m<=SIZE(NeighList,2)  .AND. NeighList(i,m)>0 )
          !Get the index of the #m-th neighbor of atom #i
          id = NeighList(i,m)
          !Get index of "atompair" for the pair #i-#j
          atompair = 0
          DO j=1,SIZE(pairs,1)
            IF( DABS(P(i,4)-pairs(j,1)) < 1.d-3 .AND. &
              & DABS(P(id,4)-pairs(j,2)) < 1.d-3      ) THEN
              atompair = j
              EXIT
            ENDIF
          ENDDO
          IF( atompair>0 ) THEN
            !Compute distance between atoms #i and #id
            Vector(:) = P(id,1:3) - P(i,1:3)
            DO u=umin,umax
              DO v=vmin,vmax
                DO w=wmin,wmax
                  distance = VECLENGTH( Vector(:) &
                          &   + DBLE(u)*H(1,:) + DBLE(v)*H(2,:) + DBLE(w)*H(3,:) )
                  IF( distance<(rdf_maxR+rdf_dr) ) THEN
                    !Add this neighbour in the appropriate &atompair" and appropriate skin
                    !n = index of skin
                    n = MAX(1,CEILING(distance/rdf_dr))
                    rdf_temp(atompair,n) = rdf_temp(atompair,n) + 1
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
          ENDIF !end if atompair>0
          !Go to next neighbour
          m=m+1
        ENDDO
        !
      ENDDO  !loop on atoms (i)
      !$OMP END PARALLEL DO
      !
      !We are done with this system: free memory
      IF(ALLOCATED(P)) DEALLOCATE(P)
      IF(ALLOCATED(NeighList)) DEALLOCATE(NeighList)
      !
      !
      !At this point, rdf_temp(:,:) contains the sum of all neighbours for this system
      !Compute total volume of current system
      CALL VOLUME_PARA(H,Vsystem)
      !Normalize partial RDFs for current system
      DO i=1,SIZE(pairs,1)
        !Pair of atoms (sp1,sp2): compute average density of second species
        DO j=1,SIZE(aentries,1)
          !Save number of atoms sp1
          IF( DABS(aentries(j,1)-pairs(i,1))<1.d-3 ) THEN
            u = NINT(aentries(j,2))
          ENDIF
          !Save average density of atoms sp2
          IF( DABS(aentries(j,1)-pairs(i,2))<1.d-3 ) THEN
            v = NINT(aentries(j,2))
            !average_dens = aentries(j,2) / Vsystem
          ENDIF
        ENDDO
        !
        w = u*v
        !
        !Loop on all skins of width rdf_dr
        DO j=1,rdf_Nsteps
          !Radius of current sphere
          rdf_radius = DBLE(j-1) * rdf_dr
          !Volume of sphere of radius R
          Vsphere = (4.d0/3.d0)*pi*(rdf_radius**3)
          !Volume of sphere of radius R+dR
          distance = rdf_radius + rdf_dr
          distance = (4.d0/3.d0)*pi*(distance**3)
          !Volume of the skin = difference between the 2 spheres
          Vskin = distance - Vsphere
          !
          !Compute normalization factor: (Nsp1*Nsp2*Vskin)/V
          rdf_norm = DBLE(w) * Vskin /Vsystem
          !
          !Normalize partial RDFs for current system
          rdf_temp(i,j) = rdf_temp(i,j) / rdf_norm
        ENDDO
      ENDDO
      !
      !Add partial RDFs of current system into rdf_partial
      !(NOTE: it will be time-averaged later, see label 300)
      rdf_partial(:,:) = rdf_partial(:,:) + rdf_temp(:,:)
      !We don't need rdf_temp any more
      IF(ALLOCATED(rdf_temp)) DEALLOCATE(rdf_temp)
      !
      !File was successfully analyzed
      Nfiles = Nfiles+1
      !
    ELSE
      !Input file doesn't exist: output a warning and go to the next file
      nwarn=nwarn+1
      CALL ATOMSK_MSG(4700,(/TRIM(inputfile)/),(/0.d0/))
      !
    ENDIF   !end if fileexists
    !
  ENDIF   !end if temp.NE."#"
  !
  290 CONTINUE
  !
ENDDO  !loop on m files
!
!
!
300 CONTINUE
CLOSE(50)
!
IF(Nfiles<=0) THEN
  !no file was analyzed => exit
  nerr=nerr+1
  GOTO 1000
ENDIF
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
ntot = 0
DO i=1,SIZE(rdf_partial,1)
  nt = 1
  l = 0
  !Loop over the two atom types of the pair
  DO k=1,2
    !Search for atoms of type pairs(i,k) in the table aentries
    DO j=1,SIZE(aentries,1)
      IF( DABS(aentries(j,1)-pairs(i,k)) < 1.d-3 ) THEN
        !Found a match
        IF(k==1) THEN
          !First atom of the pair: save number of atoms of type k=1
          nt = NINT(aentries(j,2))
          l=NINT(aentries(j,1))
        ELSE !i.e. if k==2
          !Second atom of the pair
          nt = nt * NINT(aentries(j,2))
          !If species are different, account for factor 2
          IF( NINT(aentries(j,1)).NE.l ) nt = 2*nt
        ENDIF
        GOTO 310
      ENDIF
    ENDDO
    310 CONTINUE
  ENDDO
  rdf_total(:) = rdf_total(:) + DBLE(nt)*rdf_partial(i,:)
  ntot = ntot + nt
ENDDO
!Normalize total RDF
rdf_total(:) = rdf_total(:) / DBLE(ntot)
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
!Output partial RDFs to files
IF(Nspecies>1) THEN
  atompair=0
  DO k=1,Nspecies
    DO l=k,Nspecies
      atompair=atompair+1
      !
      CALL ATOMSPECIES(aentries(k,1),sp1)
      CALL ATOMSPECIES(aentries(l,1),sp2)
      rdfdat = 'rdf_'//TRIM(sp1)//TRIM(sp2)//'.dat'
      IF(.NOT.overw) CALL CHECKFILE(rdfdat,'writ')
      OPEN(UNIT=40,FILE=rdfdat)
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
OPEN(UNIT=40,FILE=rdfdat)
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
