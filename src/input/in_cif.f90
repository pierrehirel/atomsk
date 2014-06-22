MODULE in_cif
!
!
!**********************************************************************************
!*  IN_CIF                                                                        *
!**********************************************************************************
!* This module reads a Crystallographic Information Files (CIF).                  *
!* The CIF format is officially described in the following articles:              *
!*     S.R.Hall, F.H. Allen, I.D. Brown, Acta Cryst. A 47 (1991) 655–685          *
!*     I.D. Brown, B. McMahon, Acta Cryst. B 58 (2002) 317–24                     *
!* and also on the Web site of the International Union of Crystallography:        *
!*     http://www.iucr.org/resources/cif/                                         *
!**********************************************************************************
!* (C) Oct. 2012 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 19 March 2014                                    *
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
USE constants
USE functions
USE messages
USE files
USE subroutines
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE READ_CIF(inputfile,H,P,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: msg, temp, temp2
CHARACTER(LEN=32),DIMENSION(32):: columns !contents of columns
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: atpos_done, isreduced
INTEGER:: at_occ
INTEGER:: at_sp, at_x, at_y, at_z !position of atom species and coordinates in the line
INTEGER:: i, j
INTEGER:: Naux !number of auxiliary properties
INTEGER:: Ncol, NP, sp_NP
INTEGER:: occ  !index of occupancies in AUXNAMES
INTEGER:: strlength
REAL(dp):: a, b, c, alpha, beta, gamma !supercell (conventional notation)
REAL(dp):: snumber
REAL(dp),DIMENSION(3,3),INTENT(OUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
atpos_done=.FALSE.
at_sp=0
at_x=0
at_y=0
at_z=0
Naux=0
NP=0
occ=0
sp_NP=9
a=0.d0
b=0.d0
 c=0.d0
alpha=0.d0
beta=0.d0
gamma=0.d0
!
msg = 'entering READ_CIF'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!
!
100 CONTINUE
OPEN(UNIT=30,FILE=inputfile,STATUS='UNKNOWN',ERR=800)
!
!Go back to beginning of file and store data
REWIND(30)
i=0
DO
  READ(30,'(a128)',ERR=200,END=200) temp
  temp = ADJUSTL(temp)
  !CALL ATOMSK_MSG(999,(/temp/),(/0.d0/))
  !
  IF( temp(1:14)=='_cell_length_a' ) THEN
    !Read cell length a. It may be followed by precision in brackets
    strlength = SCAN(temp,'(')
    IF(strlength<=0) strlength=LEN(temp)
    READ(temp(15:strlength-1),*,ERR=800,END=800) a
  !
  ELSEIF( temp(1:14)=='_cell_length_b' ) THEN
    !Read cell length b. It may be followed by precision in brackets
    strlength = SCAN(temp,'(')
    IF(strlength<=0) strlength=LEN(temp)
    READ(temp(15:strlength-1),*,ERR=800,END=800) b
  !
  ELSEIF( temp(1:14)=='_cell_length_c' ) THEN
    !Read cell length c. It may be followed by precision in brackets
    strlength = SCAN(temp,'(')
    IF(strlength<=0) strlength=LEN(temp)
    READ(temp(15:strlength-1),*,ERR=800,END=800) c
  !
  ELSEIF( temp(1:17)=='_cell_angle_alpha' ) THEN
    !Read cell angle alpha. It may be followed by precision in brackets
    strlength = SCAN(temp,'(')
    IF(strlength<=0) strlength=LEN(temp)
    READ(temp(18:strlength-1),*,ERR=800,END=800) alpha
    alpha = DEG2RAD(alpha)
  !
  ELSEIF( temp(1:16)=='_cell_angle_beta' ) THEN
    !Read cell angle beta. It may be followed by precision in brackets
    strlength = SCAN(temp,'(')
    IF(strlength<=0) strlength=LEN(temp)
    READ(temp(17:strlength-1),*,ERR=800,END=800) beta
    beta = DEG2RAD(beta)
  !
  ELSEIF( temp(1:17)=='_cell_angle_gamma' ) THEN
    !Read cell angle gamma. It may be followed by precision in brackets
    strlength = SCAN(temp,'(')
    IF(strlength<=0) strlength=LEN(temp)
    READ(temp(18:strlength-1),*,ERR=800,END=800) gamma
    gamma = DEG2RAD(gamma)
  !
  ELSEIF( temp(1:21)=='_chemical_formula_sum' ) THEN
    IF(.NOT.ALLOCATED(P)) THEN
      !Read the number of atoms from the formula
      !The formula should be element symbols followed by numbers
      !and separated by blank spaces, e.g. 'C4 H16 O'.
      !
      !Remove the keyword "_chemical_formula_sum"
      temp = TRIM(ADJUSTL(temp(22:)))
      !Remove the quotes if any
      strlength=SCAN(temp,"'")
      DO WHILE(strlength>0)
        temp(strlength:strlength)=" "
        strlength=SCAN(temp,"'")
      ENDDO
      temp = ADJUSTL(temp)
      !
      WRITE(msg,*) 'Found chemical formula: ', TRIM(temp)
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      !
      !READ(temp,*) temp2
      !temp = ADJUSTL(temp(LEN_TRIM(temp2)+1:))
      !
      !Parse the line to find atom element names
      NP=0
      DO WHILE(LEN_TRIM(temp)>0)
        sp_NP=0
        species=temp(1:2)
        CALL ATOMNUMBER(species,snumber)
        IF( NINT(snumber)>0 ) THEN
          !This is an atom species: read how many of them there are
          temp2=ADJUSTL(temp(3:))
          strlength=VERIFY(temp2,'0123456789')
          IF(strlength<=1) THEN
            !There is only one such atom
            sp_NP=1
          ELSE
            !Read how many such atoms exist
            READ(temp2,*,ERR=800,END=800) sp_NP
          ENDIF
          WRITE(msg,*) "Species, sp_NP: ", NINT(snumber), sp_NP
          CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
        ELSE
          !Try with only the first letter
          species=temp(1:1)
          CALL ATOMNUMBER(species,snumber)
          IF( NINT(snumber)>0 ) THEN
            !This is an atom species: read how many of them there are
            temp2=ADJUSTL(temp(2:))
            strlength=VERIFY(temp2,'0123456789')
            IF(strlength<=1) THEN
              !There is only one such atom
              sp_NP=1
            ELSE
              !Read how many such atoms exist
              READ(temp2,*,ERR=800,END=800) sp_NP
            ENDIF
            WRITE(msg,*) "Species, sp_NP: ", NINT(snumber), sp_NP
            CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
          ENDIF
        ENDIF
        !
        IF(sp_NP>0) THEN
          NP = NP+sp_NP
        ENDIF
        !
        strlength=SCAN(temp," ")
        temp = ADJUSTL(temp(strlength:))
      ENDDO
      !
      WRITE(msg,*) 'NP = ', NP
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      !
      !Allocate P
      IF(NP>0) THEN
        ALLOCATE(P(NP,4))
        P(:,:) = 0.d0
      ELSE
        !Zero atom in the formula? Something is wrong
        GOTO 800
      ENDIF
    ENDIF
  !
  ELSEIF( temp(1:5)=="loop_" ) THEN
    IF( .NOT.atpos_done ) THEN
      !Check if it is a loop over atom positions
      READ(30,'(a128)',ERR=200,END=200) temp
      temp = ADJUSTL(temp)
      IF( temp(1:5)=="_atom" ) THEN
        !It is a loop over atoms: continue reading the header,
        !count number of columns (Ncol) until reaching the actual atom positions
        at_sp=0
        Naux=0
        Ncol=0
        DO WHILE ( temp(1:1)=="_" )
          Ncol=Ncol+1
          IF( temp=="_atom_site_type_symbol" ) THEN
            at_sp = Ncol
          ELSEIF( temp=="_atom_site_label" ) THEN
            IF(at_sp==0) at_sp = Ncol
          ELSEIF( temp=="_atom_site_fract_x" ) THEN
            at_x = Ncol
          ELSEIF( temp=="_atom_site_fract_y" ) THEN
            at_y = Ncol
          ELSEIF( temp=="_atom_site_fract_z" ) THEN
            at_z = Ncol
          ELSEIF( temp=="_atom_site_occupancy" ) THEN
            at_occ = Ncol
            Naux=Naux+1
          ENDIF
          READ(30,'(a128)',ERR=200,END=200) temp
          temp = ADJUSTL(temp)
        ENDDO
        !
        IF( at_sp.NE.0 .AND. at_x.NE.0 .AND. at_y.NE.0 .AND. at_z.NE.0 ) THEN
          !This section indeed contains atom positions
          !Go back one line
          BACKSPACE(30)
          !
          WRITE(msg,*) 'Ncol = ', Ncol
          CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
          WRITE(msg,*) 'Column sp x y z: ', at_sp, at_x, at_y, at_z
          CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
          !
          IF( .NOT.ALLOCATED(P) ) THEN
            !If number of atoms was not declared before, then P is not allocated
            !=> count atoms
            NP=0
            DO
              !Read lines 
              READ(30,*,ERR=170,END=170) (columns(j),j=1,Ncol)
              !Get rid of parethesis
              strlength = SCAN(columns(at_x),'(')
              IF(strlength>0) columns(at_x) = columns(at_x)(1:strlength-1)
              strlength = SCAN(columns(at_y),'(')
              IF(strlength>0) columns(at_y) = columns(at_y)(1:strlength-1)
              strlength = SCAN(columns(at_z),'(')
              IF(strlength>0) columns(at_z) = columns(at_z)(1:strlength-1)
              !
              !Make sure that column #at_sp contains an atom species
              READ(columns(at_sp),*,ERR=170,END=170) species
              CALL ATOMNUMBER(species,snumber)
              IF( NINT(snumber)>0 ) THEN
                !This line indeed contains information about an atom
                NP = NP+1
              ENDIF
            ENDDO
            !
            170 CONTINUE
            WRITE(msg,*) 'Counted atoms, NP = ', NP
            CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
            IF(NP>0) THEN
              ALLOCATE(P(NP,4))
              P(:,:) = 0.d0
              REWIND(30)
            ELSE
              !Zero atom found? Something is wrong
              GOTO 800
            ENDIF
            !
          ELSE
            !P is allocated => we can proceed with atom positions and auxiliary properties
            WRITE(msg,*) 'Naux = ', Naux
            CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
            IF(Naux>0) THEN
              ALLOCATE( AUX(SIZE(P,1),Naux) )
              AUX(:,:) = 0.d0
              ALLOCATE( AUXNAMES(Naux) )
              Naux=1
              IF( at_occ>0 ) THEN
                occ=Naux
                AUXNAMES(occ) = "occ"
              ENDIF
            ENDIF
            !
            DO i=1,SIZE(P,1)
              READ(30,*,ERR=200,END=200) (columns(j),j=1,Ncol)
              !Get rid of parethesis
              strlength = SCAN(columns(at_x),'(')
              IF(strlength>0) columns(at_x) = columns(at_x)(1:strlength-1)
              strlength = SCAN(columns(at_y),'(')
              IF(strlength>0) columns(at_y) = columns(at_y)(1:strlength-1)
              strlength = SCAN(columns(at_z),'(')
              IF(strlength>0) columns(at_z) = columns(at_z)(1:strlength-1)
              !
              !Store data to arrays
              READ(columns(at_x),*,ERR=800,END=800) P(i,1)
              READ(columns(at_y),*,ERR=800,END=800) P(i,2)
              READ(columns(at_z),*,ERR=800,END=800) P(i,3)
              READ(columns(at_sp),*,ERR=800,END=800) temp
              temp = ADJUSTL(temp)
              species = temp(1:2)
              CALL ATOMNUMBER(species,P(i,4))
              IF( P(i,4)<=0.d0 ) THEN
                species = temp(1:1)
                CALL ATOMNUMBER(species,P(i,4))
              ENDIF
              !
              !Store auxiliary properties if any
              IF(ALLOCATED(AUX)) THEN
                Naux=1
                IF(occ>0) THEN
                  !Get rid of parethesis
                  strlength = SCAN(columns(at_occ),'(')
                  IF(strlength>0) columns(at_occ) = columns(at_occ)(1:strlength-1)
                  READ(columns(at_occ),*,ERR=800,END=800) AUX(i,occ)
                ENDIF
              ENDIF
              !
            ENDDO  !loop on atoms
            atpos_done=.TRUE.
            !
          ENDIF !end if ALLOCATED(P)
          !
        ENDIF
        !
      ELSEIF( temp(1:5)=="_symmetry" ) THEN
        !Warn user that symmetry operations are not taken into account
        nwarn=nwarn+1
        CALL ATOMSK_MSG(1704,(/""/),(/0.d0/))
      ENDIF !if "_atom"
      !
    ENDIF !if not atpos_done
    !
  ENDIF
ENDDO
!
!
!
200 CONTINUE
CLOSE(30)
IF(.NOT.atpos_done) GOTO 800
!Save cell vectors in H(:,:)
CALL CONVMAT(a,b,c,alpha,beta,gamma,H)
!Find out if coordinates are reduced or cartesian
CALL FIND_IF_REDUCED(P,isreduced)
!In case of reduced coordinates, convert them to cartesian
IF(isreduced) THEN
  CALL FRAC2CART(P,H)
ENDIF
GOTO 1000
!
!
!
800 CONTINUE
801 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE READ_CIF
!
END MODULE in_cif
