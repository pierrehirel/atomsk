MODULE out_lammps_data
!
!
!**********************************************************************************
!*  OUT_LAMMPS_DATA                                                               *
!**********************************************************************************
!* This module writes a data file for LAMMPS, that can be called                  *
!* with the  'read_data <filename>' in the LAMMPS input file, see:                *
!*     http://lammps.sandia.gov/doc/read_data.html                                *
!* Support for bonds is limited to adiabatic core-shell model, described          *
!* in LAMMPS documentation section How-to 26:                                     *
!*     http://lammps.sandia.gov/doc/Section_howto.html#howto-26                   *
!**********************************************************************************
!* (C) June 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 08 March 2023                                    *
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
USE comv
USE constants
USE functions
USE messages
USE files
USE subroutines
USE resize
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE WRITE_LMP_DATA(H,P,S,comment,AUXNAMES,AUX,outputfile)
!
CHARACTER(LEN=*),INTENT(IN):: outputfile
CHARACTER(LEN=2):: answer, skew
CHARACTER(LEN=2):: species
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(IN):: comment
LOGICAL:: charges, shells, velocities !are atom charges, shells, velocities, defined?
INTEGER:: atomID
INTEGER:: bondtype   !bond type (core-shell model)
INTEGER:: i, iloop, j, l, m
INTEGER:: molID      !position of the molecule-ID in AUX (=0 if not defined)
INTEGER:: moleculeID !actual ID of a molecule
INTEGER:: NPdigits   !number of digits in the number of particles
INTEGER:: Nspecies
INTEGER:: Ntypes     !number of particles types (atoms, or cores+shells)
INTEGER:: Nbond      !bond number (core-shell model)
INTEGER:: Nshells    !number of shells (core-shell model)
INTEGER:: Nshelltypes !number of shell types (core-shell model)
INTEGER:: vx, vy, vz, q, typecol !index of atom (or core) velocities, charges, types in AUX
INTEGER:: qs, Stypecol           !index of shell charges, types in AUX
REAL(dp):: alpha, beta, gamma, Kx, Ky, Kz !for conventional notation
REAL(dp):: Qcore, Qshell !electric charge of core, shell
REAL(dp):: smass         !mass of an atom
REAL(dp):: Smassratio=0.1d0  !ratio (mass of shell)/(mass of core)
REAL(dp):: tiltbefore
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: K   !copy of H
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: typemass  !array for storing at.number and mass of each atom type
REAL(dp),DIMENSION(:,:),POINTER:: Ppoint  !pointer to P or R
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN),TARGET:: P !atom positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE,TARGET:: R  !copy of P (if necessary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: S !positions of shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: atypes    !atom types and their number
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentries  !atom species and their number
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: aentriesS !atom species and their number (shells)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(IN):: AUX !auxiliary properties
!
!
!Initialize variables
 charges = .FALSE.
shells = .FALSE.
velocities = .FALSE.
atomID=0
i=0
j=0
q = 0
qs = 0
typecol = 0
Stypecol = 0
vx = 0
vy = 0
vz = 0
Nbond = 0
Nshells = 0
Nshelltypes = 0
Nspecies = 0
K=H
IF(ALLOCATED(atypes)) DEALLOCATE(atypes)
IF(ALLOCATED(aentries)) DEALLOCATE(aentries)
IF(ALLOCATED(aentriesS)) DEALLOCATE(aentriesS)
IF(ALLOCATED(R)) DEALLOCATE(R)
Ppoint=>P
!
!
msg = 'entering WRITE_LMP_DATA'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
!Find how many different species are in P
CALL FIND_NSP(P(:,4),aentries)
Ntypes = SIZE(aentries,1)
IF(verbosity==4) THEN
  msg = "aentries:   at.number |  occurrence"
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  DO i=1,SIZE(aentries,1)
    WRITE(msg,'(10X,f9.3,a5,i5)') aentries(i,1), "   | ", NINT(aentries(i,2))
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  ENDDO
ENDIF
!
!Count number of digits composing the number of particles
WRITE(temp,*) SIZE(P,1)
NPdigits = LEN_TRIM(ADJUSTL(temp))
!
!If shells are present (ionic core-shell model), count number of shells and number of bonds
IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
  !count how many species have shells
  CALL FIND_NSP( S(:,4) , aentriesS )
  Nshells = 0
  Nshelltypes = 0
  !remove zeros
  DO i=1,SIZE(aentriesS,1)
    IF( aentriesS(i,1) < 0.1d0 ) THEN
      DO j=i,SIZE(aentriesS,1)-1
        aentriesS(j,:) = aentriesS(j+1,:)
        aentriesS(j+1,:) = 0.d0
      ENDDO
    ENDIF
    IF( aentriesS(i,1) > 0.1d0 ) THEN
      Nshelltypes = Nshelltypes+1
      Nshells = Nshells + NINT(aentriesS(i,2))
    ENDIF
  ENDDO
  CALL RESIZE_DBLEARRAY2(aentriesS,Nshelltypes,SIZE(aentriesS,2))
  !Update total number of atom types
  Ntypes = Ntypes + Nshelltypes
  IF(verbosity==4) THEN
    msg = "aentriesS:  at.number |  occurrence"
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    DO i=1,SIZE(aentriesS,1)
      WRITE(msg,'(10X,f9.3,a5,i5)') aentriesS(i,1), "   | ", NINT(aentriesS(i,2))
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    ENDDO
  ENDIF
ENDIF
!
WRITE(msg,*) "Total number of particles   NP = ", SIZE(P,1)
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,*) "Number of shells       Nshells = ", Nshells
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,*) "      Number of particle types = ", Ntypes
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,*) "         Number of shell types = ", Nshelltypes
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Check if auxiliary properties relevant to LAMMPS data files are present
IF( ALLOCATED(AUXNAMES) ) THEN
  DO i=1,SIZE(AUXNAMES)
    IF( TRIM(ADJUSTL(AUXNAMES(i)))=='vx' ) THEN
      vx = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='vy' ) THEN
      vy = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='vz' ) THEN
      vz = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='molID' ) THEN
      molID = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='q' ) THEN
      q = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='qs' ) THEN
      qs = i
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='type' ) THEN
      typecol = i
      IF( Nshells==0 ) THEN
        !Could not determine number of atom types before
        !Count how many different atom types are in AUX
        !CALL FIND_NSP(AUX(:,typecol),atypes)
        Ntypes = MAXVAL(AUX(:,typecol))
        !Verify that atom types are all greater than zero
        IF( ANY(AUX(:,typecol)<0.9d0) ) THEN
          nwarn=nwarn+1
          CALL ATOMSK_MSG(3714,(/""/),(/0.d0/))
        ENDIF
      ENDIF
    ELSEIF( TRIM(ADJUSTL(AUXNAMES(i)))=='Stype' ) THEN
      Stypecol = i
    ENDIF
  ENDDO
  !
  IF( vx>0 .AND. vy>0 .AND. vz>0 ) velocities = .TRUE.
  IF( q>0 ) charges = .TRUE.
  !
  WRITE(msg,*) "AUX columns: vx vy vz ", vx, vy, vz, "; type ", typecol
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ENDIF
!
!
IF( typecol>0 ) THEN
  !A column for atom types is defined
  !For each atom type, find its atomic number and mass and store them in array "typemass"
  ALLOCATE( typemass(Ntypes,2) )
  typemass(:,:) = 0.d0
  
  DO j=1,Ntypes-Nshelltypes
    i = 0
    DO WHILE( i<SIZE(P,1) .AND. typemass(j,1)<0.1d0 )
      i = i+1
      IF( NINT(AUX(i,typecol)) == j ) THEN
        typemass(j,1) = P(i,4)
        CALL ATOMSPECIES(P(i,4),species)
        CALL ATOMMASS(species , typemass(j,2))
      ENDIF
    ENDDO
  ENDDO
  !
  !Then, append shells if any
  IF( Nshelltypes>0 .AND. ALLOCATED(S) ) THEN
    i = 0
    IF( Stypecol>0 ) THEN
      !Use types as defined in AUX
      msg = "Using shell types defined in AUX"
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      DO WHILE( i<SIZE(S,1) .AND. j<=Nshelltypes )
        i=i+1
        IF( NINT(AUX(i,Stypecol)) == j+1 ) THEN
          IF( S(i,4)>0.9d0 ) THEN
            j = j+1
            typemass(j,1) = S(i,4)
            CALL ATOMSPECIES(S(i,4),species)
            CALL ATOMMASS(species , typemass(j,2))
            typemass(j,2) = typemass(j,2)*Smassratio
          ENDIF
        ENDIF
      ENDDO
    ELSE
      !Generate new "types" for shells
      msg = "Generating new particle types for shells..."
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      !Try to put shells in same order as cores
      !At this point iloop = number of types of core
      iloop = Ntypes-Nshelltypes
      !Parse typemass(:,:) from i=1 to number of core types (iloop)
      i=0
      j=0
      DO WHILE( i<iloop .AND. j<=Nshelltypes )
        i = i+1
        !Check if atoms of type #i have a shell
        DO m=1,SIZE(aentriesS,1)  !Loop on all shell types
          IF( DABS(typemass(i,1)-aentriesS(m,1)) < 0.1d0 ) THEN
            !Shells of type #m have the same species as cores of type #i
            j = j+1
            typemass(iloop+j,:) = typemass(i,:)
            EXIT
          ENDIF
        ENDDO
      ENDDO
    ENDIF
  ENDIF
ENDIF
!
!
!
100 CONTINUE
IF( verbosity>=4 ) THEN
  IF( ALLOCATED(typemass) ) THEN
    msg = "Particle Type | At.Number |  Mass"
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    DO i=1,SIZE(typemass,1)
      WRITE(msg,'(8X,i3, 2f12.3)') i, typemass(i,1), typemass(i,2)
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    ENDDO
  ENDIF
ENDIF
!
IF(ofu.NE.6) THEN
  OPEN(UNIT=ofu,FILE=outputfile,STATUS='UNKNOWN',ERR=500)
ENDIF
!
!Write header of data file
IF( ALLOCATED(comment) .AND. SIZE(comment)>0 ) THEN
  WRITE(ofu,*) TRIM(ADJUSTL(comment(1)))
  WRITE(ofu,*) ' '
ENDIF
WRITE(ofu,*) SIZE(P,1)+Nshells, ' atoms'
IF( Nshells>0 ) THEN
  WRITE(ofu,*) Nshells, ' bonds'
ENDIF
!
!Determine how many different species are present
WRITE(ofu,*) Ntypes, ' atom types'
IF( Nshelltypes>0 ) THEN
  WRITE(ofu,*) Nshelltypes, ' bond types'
ENDIF
WRITE(ofu,*) ''
!
!
!Check that supercell vectors form a lower rectangular matrix
IF( DABS(H(1,2))>1.d-12 .OR. DABS(H(1,3))>1.d-12 .OR. DABS(H(2,3))>1.d-12 ) THEN
  !Prompt user whether to rotate the system
  nwarn = nwarn+1
  CALL ATOMSK_MSG(3704,(/''/),(/0.d0/))
  READ(*,*) answer
  IF(answer==langyes .OR. answer==langBigYes) THEN
    !The following does the same as option "-alignx"
    !Copy P to R
    ALLOCATE( R( SIZE(P,1), SIZE(P,2) ) )
    R=P
    !
    !Convert to fractional coordinates
    CALL CART2FRAC(R,K)
    !
    !Convert the matrix K into conventional notation
    Kx = VECLENGTH(K(1,:))
    Ky = VECLENGTH(K(2,:))
    Kz = VECLENGTH(K(3,:))
    alpha = ANGVEC(K(2,:),K(3,:))
    beta  = ANGVEC(K(3,:),K(1,:))
    gamma = ANGVEC(K(1,:),K(2,:))
    !
    !Then convert this conventional notation into lower-triangular matrix K
    CALL CONVMAT(Kx,Ky,Kz,alpha,beta,gamma,K)
    !
    !Convert back to cartesian coordinates
    CALL FRAC2CART(R,K)
    !Modify pointer
    Ppoint=>R
    CALL ATOMSK_MSG(2051,(/''/),(/0.d0/))
    !
  ELSE
    !If user does not want to rotate the system,
    !display a final warning just to be sure
    CALL ATOMSK_MSG(3709,(/''/),(/0.d0/))
  ENDIF
ENDIF
!
!
!Write supercell data
IF( VECLENGTH(K(1,:)) > 1.d-12 ) THEN
  WRITE(ofu,160) MIN(K(1,1), zero), MAX(K(1,1), zero), '  xlo xhi'
ENDIF
IF( VECLENGTH(K(2,:)) > 1.d-12 ) THEN
  WRITE(ofu,160) MIN(K(2,2), zero), MAX(K(2,2), zero), '  ylo yhi'
ENDIF
IF( VECLENGTH(K(3,:)) > 1.d-12 ) THEN
  WRITE(ofu,160) MIN(K(3,3), zero), MAX(K(3,3), zero), '  zlo zhi'
ENDIF
!
!LAMMPS requires that skew parameters are less than half the box
!length in each direction. If it is not the case, warn the user, and
!propose to "unskew" the box
IF( DABS(K(2,1))>1.d-12 .OR. DABS(K(3,1))>1.d-12 .OR. DABS(K(3,2))>1.d-12 ) THEN
  DO i=2,3
    IF(i==2) skew(2:2)='y'
    IF(i==3) skew(2:2)='z'
    !
    DO j=i-1,1,-1
      !Don't consider diagonal elements
      IF(i.NE.j) THEN
        IF(j==1) skew(1:1)='x'
        IF(j==2) skew(1:1)='y'
        !
        !Check if tilt is too large
        IF( DABS(K(i,j))>0.5d0*K(j,j)+1.d-12 ) THEN
          nwarn=nwarn+1
          !Tilt is too large: ask if it should be reduced
          CALL ATOMSK_MSG(3705,(/skew/),(/0.d0/))
          READ(*,*) answer
          IF(answer==langyes .OR. answer==langBigYes) THEN
            !Unskew tilt K(i,j)
            tiltbefore = K(i,j)
            iloop=0
            !If tilt is too large, remove the matching box vector
            DO WHILE( K(i,j)>0.5d0*K(j,j)+1.d-12 )
              K(i,:) = K(i,:) - K(j,:)
              iloop=iloop+1
              IF(iloop>100) EXIT
            ENDDO
            !If tilt is too negative, add the matching box vector
            DO WHILE( K(i,j)<-0.5d0*K(j,j)-1.d-12 )
              K(i,:) = K(i,:) + K(j,:)
              iloop=iloop+1
              IF(iloop>100) EXIT
            ENDDO
            !Check that the loops did not go crazy
            IF(iloop>100) THEN
              !After 100 loops no solution was found
              !Display a warning and restore initial value
              nwarn=nwarn+1
              CALL ATOMSK_MSG(3706,(/skew/),(/0.d0/))
              K(i,j) = tiltbefore
            ELSEIF(iloop>0) THEN
              !This tilt was corrected: display message
              CALL ATOMSK_MSG(3004,(/skew/),(/0.d0/))
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  WRITE(ofu,161) K(2,1), K(3,1), K(3,2), ' xy xz yz'
ENDIF
WRITE(ofu,*) ''
160 FORMAT(f20.12,1X,f20.12,a9)
161 FORMAT(3(f20.12,1X),a9)
!
!
!
200 CONTINUE
!Write the Masses section
WRITE(ofu,'(a6)') "Masses"
WRITE(ofu,*) ''
IF( Nshells>0 .AND. ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
  !In case of core/shell model, write mass ratio
  !Write mass of cores (and chemical species as comment)
  IF( typecol>0 ) THEN
    !Atom types are defined in auxiliary properties
    !and were saved in array "typemass"
    !=> write that information
    DO i=1,Ntypes-Nshelltypes
      CALL ATOMSPECIES(typemass(i,1),species)
      WRITE(temp,'(f16.8)') typemass(i,2)*(1.d0-Smassratio)
      WRITE(msg,*) i, "  ", TRIM(ADJUSTL(temp))
      msg(40:) = "# "//species//" core"
      WRITE(ofu,*) TRIM(msg)
    ENDDO
    DO i=Ntypes-Nshelltypes+1,Ntypes
      CALL ATOMSPECIES(typemass(i,1),species)
      WRITE(temp,'(f16.8)') typemass(i,2)*Smassratio
      WRITE(msg,*) i, "  ", TRIM(ADJUSTL(temp))
      msg(40:) = "# "//species//" shell"
      WRITE(ofu,*) TRIM(msg)
    ENDDO
  ELSE
    !Atom types are undefined
    !=> set them by order of appearence in array "aentries"
    DO i=1,SIZE(aentries,1)
      !Determine the mass of this type of atom
      CALL ATOMSPECIES(aentries(i,1),species)
      CALL ATOMMASS(species,smass)
      WRITE(temp,'(f16.8)') smass*(1.d0-Smassratio)
      WRITE(msg,*) i, "  ", TRIM(ADJUSTL(temp))
      msg(40:) = "# "//species//" core"
      WRITE(ofu,*) TRIM(msg)
    ENDDO
    !Same with shells
    j = SIZE(aentries,1)
    DO i=1,SIZE(aentriesS,1)
      IF( aentriesS(i,1)>0.1d0 ) THEN
        j=j+1
        CALL ATOMSPECIES(aentriesS(i,1),species)
        CALL ATOMMASS(species,smass)
        WRITE(temp,'(f16.8)') smass*Smassratio
        WRITE(msg,*) j, "  ", TRIM(ADJUSTL(temp))
        msg(40:) = "# "//species//" shell"
        WRITE(ofu,*) TRIM(msg)
      ENDIF
    ENDDO
    !
  ENDIF
  !
ELSE
  !Write mass of each atom type (and chemical species as comment)
  IF( typecol>0 ) THEN
    !Atom types are defined in auxiliary properties
    !and were saved in array "typemass"
    !=> write that
    DO i=1,SIZE(typemass,1)
      CALL ATOMSPECIES(typemass(i,1),species)
      WRITE(temp,'(f16.8)') typemass(i,2)
      WRITE(msg,*) i, "  ", TRIM(ADJUSTL(temp))
      msg(40:) = "# "//species
      WRITE(ofu,*) TRIM(msg)
    ENDDO
  ELSE
    !Atom types are undefined
    !=> set them by order of appearence in array "aentries"
    DO i=1,SIZE(aentries,1)
      !Determine the mass of this type of atom
      CALL ATOMSPECIES(aentries(i,1),species)
      CALL ATOMMASS(species,smass)
      WRITE(temp,'(f16.8)') smass
      !Write atom type and mass to file
      WRITE(msg,*) i, "  ", TRIM(ADJUSTL(temp))
      msg(40:) = "# "//species
      WRITE(ofu,*) TRIM(msg)
    ENDDO
  ENDIF
ENDIF
WRITE(ofu,*) ''
!
!
molID = 0
Nspecies = 0
!Write atom positions
IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
  WRITE(ofu,'(a12)') 'Atoms # full'
ELSEIF (q>0 ) THEN
  WRITE(ofu,'(a14)') 'Atoms # charge'
ELSE
  WRITE(ofu,'(a14)') 'Atoms # atomic'
ENDIF
WRITE(ofu,*) ''
!
IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
  !Ion shells are present (core-shell model)
  !=> the format is described in Section How-to 26
  !   http://lammps.sandia.gov/doc/Section_howto.html#howto-26
  !
  !(1) Write positions of cores (atoms) and shells
  !    with format "full": atom-ID molecule-ID atom-type q x y z
  atomID=0
  moleculeID = 0
  DO i=1,SIZE(Ppoint,1)
    atomID=atomID+1
    !If "molecule-ID" is defined in AUX, it is used,
    !otherwise the moleculeID is incremented at each new atom-shell pair.
    IF( molID>0 ) THEN
      moleculeID = NINT(AUX(i,molID))
    ELSE
      moleculeID = moleculeID+1
    ENDIF
    IF( q>0 ) THEN
      Qcore = AUX(i,q)
    ELSE
      Qcore = 0.d0
    ENDIF
    IF( typecol>0 ) THEN
      !Atom types are in auxiliary properties, use it
      Nspecies = NINT(AUX(i,typecol))
    ELSE
      !Replace species by atom types in their order of appearance
      DO iloop=1,SIZE(aentries,1)
        IF( aentries(iloop,1)==NINT(Ppoint(i,4)) ) Nspecies = iloop
      ENDDO
    ENDIF
    WRITE(ofu,212) atomID, moleculeID, Nspecies, Qcore, Ppoint(i,1), Ppoint(i,2), Ppoint(i,3)
    !
    IF( NINT(S(i,4)) == NINT(P(i,4)) ) THEN
      !Write shell info immediately after its core
      !NOTE: i is the index of core in P(:,:) AND the index of its shell in S(:,:)
      !     But i cannot be used as "atom-ID" in LAMMPS data file!
      !     => use a separate variable "atomID"
      atomID=atomID+1
      IF( qs>0 ) THEN
        Qshell = AUX(i,qs)
      ELSE
        Qshell = 0.d0
      ENDIF
      !Determine type of this shell
      IF( Stypecol>0 ) THEN
        !Shell types are in auxiliary properties, use it
        Nspecies = AUX(i,Stypecol)
      ELSE
        !Shell types are undefined
        !=> set them by order of appearence in array "aentriesS"
        DO iloop=1,SIZE(aentriesS,1)
          IF( NINT(S(i,4))==NINT(aentriesS(iloop,1)) ) THEN
            Nspecies = Ntypes - Nshelltypes + iloop
          ENDIF
        ENDDO
      ENDIF
      WRITE(ofu,212) atomID, moleculeID, Nspecies, Qshell, S(i,1), S(i,2), S(i,3)
    ENDIF
  ENDDO
  !
  !(2) Write bond pairs with format "bondID bondtype atom1 atom2"
  !    ID   = bond number (1-Nbonds)
  !    type = bond type (1-Nbondtype)
  !    atom1,atom2 = IDs of 1st,2nd atoms in bond
  WRITE(ofu,*) ''
  WRITE(ofu,'(a5)') 'Bonds'
  WRITE(ofu,*) ''
  j=0
  DO i=1,SIZE(S,1)
    j=j+1
    IF( S(i,4)>0.1d0 ) THEN
      !Atom #i has a shell
      !Determine the type of that bond
      bondtype = 0
      DO l=1,SIZE(aentriesS,1)
        IF( aentriesS(l,1)>0.1d0 .AND. NINT(aentriesS(l,1)) == NINT(S(i,4)) ) THEN
          bondtype = l
          EXIT
        ENDIF
      ENDDO
      !Ionic core #j is bonded with ionic shell #j+1
      Nbond = Nbond+1 ! = ID (bond number)
      WRITE(ofu,*) Nbond, bondtype, j, j+1
      !
      j=j+1
    ENDIF
  ENDDO
  !
ELSE
  !No ionic shells, only atoms are present
  IF( Ntypes>1 ) THEN
    IF( charges ) THEN
      !Atom charges are defined
      IF( typecol>0 ) THEN
        !Atom types are in auxiliary properties, use it
        DO i=1,SIZE(Ppoint,1)
          WRITE(ofu,210) i, NINT(AUX(i,typecol)), AUX(i,q), Ppoint(i,1), Ppoint(i,2), Ppoint(i,3)
        ENDDO
        !
      ELSE
        !Replace species by atom types in their order of appearance
        DO i=1,SIZE(Ppoint,1)
          DO j=1,SIZE(aentries,1)
            IF( NINT(aentries(j,1))==NINT(Ppoint(i,4)) ) Nspecies = j
          ENDDO
          WRITE(ofu,210) i, Nspecies, AUX(i,q), Ppoint(i,1), Ppoint(i,2), Ppoint(i,3)
        ENDDO
      ENDIF
      !
    ELSE
      !No atom charge defined
      IF( typecol>0 ) THEN
        !Atom types are in auxiliary properties, use it
        DO i=1,SIZE(Ppoint,1)
          WRITE(ofu,211) i, NINT(AUX(i,typecol)), Ppoint(i,1), Ppoint(i,2), Ppoint(i,3)
        ENDDO
        !
      ELSE
        !Replace species by atom types in their order of appearance
        DO i=1,SIZE(Ppoint,1)
          DO j=1,SIZE(aentries,1)
            IF( NINT(aentries(j,1))==NINT(Ppoint(i,4)) ) Nspecies = j
          ENDDO
          WRITE(ofu,211) i, Nspecies, Ppoint(i,1), Ppoint(i,2), Ppoint(i,3)
        ENDDO
      ENDIF
    ENDIF
    !
  ELSE
    Nspecies = 1 !only one type of atom
    !
    IF( charges ) THEN
      DO i=1,SIZE(Ppoint,1)
        !Format when atom charges are defined
        WRITE(ofu,210) i, Nspecies, AUX(i,q), Ppoint(i,1), Ppoint(i,2), Ppoint(i,3)
      ENDDO
      !
    ELSE
      DO i=1,SIZE(Ppoint,1)
        !Format when no atom charge is defined
        WRITE(ofu,211) i, Nspecies, Ppoint(i,1), Ppoint(i,2), Ppoint(i,3)
      ENDDO
    ENDIF
    !
  ENDIF
  !
ENDIF
210 FORMAT(i10,2X,i3,2X,f9.6,2X,3(f20.12,1X))
211 FORMAT(i10,2X,i3,2X,3(f20.12,1X))
212 FORMAT(i10,2X,i10,2X,i3,2X,f9.6,2X,3(f20.12,1X))
!
!Write velocities
IF( velocities ) THEN
  WRITE(ofu,*) ''
  WRITE(ofu,'(a10)') 'Velocities'
  WRITE(ofu,*) ''
  DO i=1,SIZE(Ppoint(:,1))
    WRITE(ofu,'(i10,2X,3(e16.8,1X))') i, AUX(i,vx), AUX(i,vy), AUX(i,vz)
  ENDDO
ENDIF
!
!Write shell positions (if any)
!IF( ALLOCATED(S) .AND. SIZE(S,1)>0 ) THEN
!  WRITE(ofu,*) ''
!  WRITE(ofu,'(a10)') 'Shells'
!  WRITE(ofu,*) ''
!  DO i=1,SIZE(S,1)
!    WRITE(ofu,'(i8,2X,3(f16.8,1X))') NINT(S(i,4)), (S(i,j), j=1,3)
!  ENDDO
!ENDIF
!
!
!
500 CONTINUE
IF(ofu.NE.6) THEN
  CLOSE(ofu)
ENDIF
msg = "LAMMPS data"
temp = outputfile
CALL ATOMSK_MSG(3002,(/msg,temp/),(/0.d0/))
GOTO 1000
!
!
!
800 CONTINUE
!
!
!
1000 CONTINUE
IF(ALLOCATED(atypes)) DEALLOCATE(atypes)
IF(ALLOCATED(aentries)) DEALLOCATE(aentries)
IF(ALLOCATED(R)) DEALLOCATE(R)
!
END SUBROUTINE WRITE_LMP_DATA
!
END MODULE out_lammps_data
