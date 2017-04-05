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
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 27 March 2017                                    *
!**********************************************************************************
!* Note on how Biso and Usio parameters are handled (by J. Barthel)               *
!*     The data is stored in Biso form, thus Uiso input is translated here to     *
!*     Biso = 8 * Pi**2 * Usio (in units of angstrom squared).                    *
!*     The auxiliary name is set to "biso".                                       *
!*     The decision about whether the input is Biso or Usio is made depending on  *
!*     the presence of the CIF label                                              *
!*         _atom_site_B_iso_or_equiv   -> Biso                                    *
!*     or  _atom_site_U_iso_or_equiv   -> Uiso.                                   *
!*     Sometimes CIF files define the parameter type by an extra label            *
!*         _atom_site_thermal_displace_type                                       *
!*     Since this label is not always present, I decided to rely on the two       *
!*     other labels as described above.                                           *
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
USE symops
!
IMPLICIT NONE
!
!
CONTAINS
!
!
SUBROUTINE READ_CIF(inputfile,H,P,comment,AUXNAMES,AUX)
!
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=64):: j_name, j_volume, j_page1, j_page2, j_year  !article information
CHARACTER(LEN=128):: chemical_name, chemical_formula
CHARACTER(LEN=128):: msg, temp, temp2
CHARACTER(LEN=32),DIMENSION(32):: columns !contents of columns
CHARACTER(LEN=1),DIMENSION(:),ALLOCATABLE:: Wyckoff_letters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: symop_list  !list of symmetry operations
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(OUT):: comment
LOGICAL:: atpos_done, symops_done, isreduced
INTEGER:: at_occ     !atom site occupancy
INTEGER:: at_biso
INTEGER:: at_uiso
INTEGER:: at_Wyckoff
INTEGER:: at_sp, at_x, at_y, at_z !position of atom species and coordinates in the line
INTEGER:: sgnumber  !space group number
INTEGER:: sy_pxyz !position of symmetry operations string in the line
INTEGER:: i, j, k, m, n
INTEGER:: Naux !number of auxiliary properties
INTEGER:: Ncol, NP, sp_NP, q
INTEGER:: Nsym !number of symmetry operations
INTEGER:: occ  !index of occupancies in AUXNAMES
INTEGER:: biso !index of Biso in AUXNAMES
INTEGER:: strlength
REAL(dp):: a, b, c, alpha, beta, gamma !supercell (conventional notation)
REAL(dp):: snumber, sbiso
REAL(dp):: tempreal
REAL(dp),DIMENSION(3,3),INTENT(OUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: P   !positions of atoms
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S               !positions of shells (remains unallocated in this module)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: AUX !auxiliary properties
!
!
!Initialize variables
 chemical_name = ""
 chemical_formula = ""
j_name = ""
j_volume = ""
j_page1 = ""
j_page2 = ""
j_year = ""
atpos_done=.FALSE.
symops_done=.FALSE.
at_sp=0
at_x=0
at_y=0
at_z=0
at_occ=0
at_biso=0
at_uiso=0
at_Wyckoff=0
sgnumber=0
sy_pxyz=0
Naux=0
Ncol=0
NP=0
Nsym=0
occ=0
biso=0
sp_NP=9
q=0
a=0.d0
b=0.d0
c=0.d0
alpha=0.d0
beta=0.d0
gamma=0.d0
IF (ALLOCATED(S)) DEALLOCATE(S)
IF (ALLOCATED(symops_trf)) DEALLOCATE(symops_trf)
IF (ALLOCATED(Wyckoff_letters)) DEALLOCATE(Wyckoff_letters)
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
DO !Main reading loop
  READ(30,'(a128)',ERR=200,END=200) temp
  temp = ADJUSTL(temp)
  !CALL ATOMSK_MSG(999,(/temp/),(/0.d0/))
  !
  !Parse cif content for main keywords
  !
  IF( temp(1:18)=='_journal_name_full' ) THEN
    !Read journal name
    j_name = ADJUSTL(temp(19:))
    !Remove quotes if any
    strlength = SCAN(j_name,"'")
    DO WHILE (strlength>0)
      j_name(strlength:strlength) = " "
      strlength = SCAN(j_name,"'")
    ENDDO
    strlength = SCAN(j_name,'"')
    DO WHILE (strlength>0)
      j_name(strlength:strlength) = " "
      strlength = SCAN(j_name,'"')
    ENDDO
  !
  ELSEIF( temp(1:15)=='_journal_volume' ) THEN
    !Read journal volume
    j_volume = ADJUSTL(temp(16:))
  !
  ELSEIF( temp(1:19)=='_journal_page_first' ) THEN
    !Read journal page
    j_page1 = ADJUSTL(temp(20:))
  !
  ELSEIF( temp(1:18)=='_journal_page_last' ) THEN
    !Read journal page
    j_page2 = ADJUSTL(temp(19:))
  !
  ELSEIF( temp(1:13)=='_journal_year' ) THEN
    !Read journal page
    j_year = ADJUSTL(temp(14:))
  !
  ELSEIF( temp(1:14)=='_cell_length_a' ) THEN
    !Read cell length a. It may be followed by precision in brackets
    IF( NINT(a)==0 ) THEN
      strlength = SCAN(temp,'(')
      IF(strlength<=0) strlength=LEN(temp)
      READ(temp(15:strlength-1),*,ERR=800,END=800) a
    ENDIF
  !
  ELSEIF( temp(1:14)=='_cell_length_b' ) THEN
    !Read cell length b. It may be followed by precision in brackets
    IF( NINT(b)==0 ) THEN
      strlength = SCAN(temp,'(')
      IF(strlength<=0) strlength=LEN(temp)
      READ(temp(15:strlength-1),*,ERR=800,END=800) b
    ENDIF
  !
  ELSEIF( temp(1:14)=='_cell_length_c' ) THEN
    !Read cell length c. It may be followed by precision in brackets
    IF( NINT(c)==0 ) THEN
      strlength = SCAN(temp,'(')
      IF(strlength<=0) strlength=LEN(temp)
      READ(temp(15:strlength-1),*,ERR=800,END=800) c
    ENDIF
  !
  ELSEIF( temp(1:17)=='_cell_angle_alpha' ) THEN
    !Read cell angle alpha. It may be followed by precision in brackets
    IF( NINT(alpha)==0 ) THEN
      strlength = SCAN(temp,'(')
      IF(strlength<=0) strlength=LEN(temp)
      READ(temp(18:strlength-1),*,ERR=800,END=800) alpha
      alpha = DEG2RAD(alpha)
    ENDIF
  !
  ELSEIF( temp(1:16)=='_cell_angle_beta' ) THEN
    !Read cell angle beta. It may be followed by precision in brackets
    IF( NINT(beta)==0 ) THEN
      strlength = SCAN(temp,'(')
      IF(strlength<=0) strlength=LEN(temp)
      READ(temp(17:strlength-1),*,ERR=800,END=800) beta
      beta = DEG2RAD(beta)
    ENDIF
  !
  ELSEIF( temp(1:17)=='_cell_angle_gamma' ) THEN
    !Read cell angle gamma. It may be followed by precision in brackets
    IF( NINT(gamma)==0 ) THEN
      strlength = SCAN(temp,'(')
      IF(strlength<=0) strlength=LEN(temp)
      READ(temp(18:strlength-1),*,ERR=800,END=800) gamma
      gamma = DEG2RAD(gamma)
    ENDIF
  !
  ELSEIF( temp(1:25)=='_chemical_name_systematic' ) THEN
    IF( LEN_TRIM(chemical_name) == 0 ) THEN
      temp = TRIM(ADJUSTL(temp(26:)))
      !Remove the quotes if any
      strlength=SCAN(temp,"'")
      DO WHILE(strlength>0)
        temp(strlength:strlength)=" "
        strlength=SCAN(temp,"'")
      ENDDO
      temp = ADJUSTL(temp)
      chemical_name = TRIM(temp)
    ENDIF
  !
  ELSEIF( temp(1:22)=='_chemical_name_mineral' ) THEN
    temp = TRIM(ADJUSTL(temp(23:)))
    !Remove the quotes if any
    strlength=SCAN(temp,"'")
    DO WHILE(strlength>0)
      temp(strlength:strlength)=" "
      strlength=SCAN(temp,"'")
    ENDDO
    temp = ADJUSTL(temp)
    chemical_name = TRIM(temp)
  !
  ELSEIF( temp(1:21)=='_chemical_formula_sum' ) THEN
    !The formula should be element symbols followed by numbers
    !and separated by blank spaces, e.g. 'C4 H16 O'.
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
    chemical_formula = TRIM(temp)
    !
    WRITE(msg,*) 'Found chemical formula: ', TRIM(temp)
    CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
    !
    !
  ELSEIF( temp(1:27)=='_symmetry_Int_Tables_number' ) THEN
    !Read space group number
    READ(temp(28:),*,ERR=190,END=190) sgnumber
    !Verify that it is a valid space group number; otherwise reset it to zero
    IF( sgnumber <1 .OR. sgnumber > 230 ) THEN
      sgnumber = 0
    ENDIF
    !
    !
  ELSEIF( temp(1:22)=='_space_group_IT_number' ) THEN
    !Read space group number
    READ(temp(23:),*,ERR=190,END=190) sgnumber
    !Verify that it is a valid space group number; otherwise reset it to zero
    IF( sgnumber <1 .OR. sgnumber > 230 ) THEN
      sgnumber = 0
    ENDIF
    !
    !
  ELSEIF( temp(1:30)=='_symmetry_space_group_name_H-M' ) THEN
    !Read Hermann-Mauguin symbol (only if space group number was not already found)
    IF( sgnumber==0 ) THEN
      temp2 = ADJUSTL(temp(31:))
      !Remove quotes if any
      strlength=SCAN(temp2,"'")
      DO WHILE(strlength>0)
        temp2(strlength:strlength) = " "
        strlength=SCAN(temp2,"'")
      ENDDO
      temp2 = ADJUSTL(temp2)
      !Transform H-M symbol into a space group number
      CALL SG_NAMGETNUM(temp2,sgnumber)
    ENDIF
    !
    !
  ELSEIF( temp(1:5)=="loop_" ) THEN
    !
    !Check if it is a loop over atom positions or symmetry operations
    READ(30,'(a128)',ERR=200,END=200) temp
    temp = ADJUSTL(temp)
    IF( temp(1:10)=="_atom_site" .AND. .NOT.atpos_done) THEN
      !It is a loop over atoms and we still need to read it:
      !continue reading the header,
      !count number of columns (Ncol) until reaching the actual atom positions
      at_sp=0
      Naux=0
      Ncol=0
      DO WHILE ( temp(1:1)=="_" ) ! begin "_atom" keyword loop
        Ncol=Ncol+1
        IF( temp=="_atom_site_type_symbol" ) THEN ! begin "_atom" keyword parser
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
        ELSEIF( temp=="_atom_site_Wyckoff_symbol" ) THEN
          at_Wyckoff = Ncol
        ELSEIF( temp=="_atom_site_B_iso_or_equiv" ) THEN
          ! a well defined CIF file contains either biso or uiso,
          ! not both at the same time.
          at_biso = Ncol
          ! mark the auxiliary input as biso
          at_uiso = 0
          Naux=Naux+1
        ELSEIF( temp=="_atom_site_U_iso_or_equiv" ) THEN
          at_biso = Ncol
          ! mark the auxiliary input as uiso
          at_uiso = 1
          Naux=Naux+1
        ENDIF ! end of "_atom" keyword parser
        !Read next keyword of the loop
        READ(30,'(a128)',ERR=200,END=200) temp
        temp = ADJUSTL(temp)
        !
      ENDDO ! end of "_atom" keyword loop
      !
      !Check if this "_atom_" loop contains all required data
      IF( at_sp.NE.0 .AND. at_x.NE.0 .AND. at_y.NE.0 .AND. at_z.NE.0 ) THEN
        !This section indeed contains atom species and positions
        !Check if this string contains the charge of the atom
        strlength = SCAN(temp,'+')
        IF( strlength>0 ) THEN
          q = 1
        ELSE
          strlength = SCAN(temp,'-')
          IF( strlength>0 ) THEN
            q = 1
          ENDIF
        ENDIF
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
            READ(30,'(a128)',ERR=170,END=170) temp
            READ(temp,*,ERR=170,END=170) (columns(j),j=1,Ncol)
            !Get rid of parenthesis
            strlength = SCAN(columns(at_x),'(')
            IF(strlength>0) columns(at_x) = columns(at_x)(1:strlength-1)
            strlength = SCAN(columns(at_y),'(')
            IF(strlength>0) columns(at_y) = columns(at_y)(1:strlength-1)
            strlength = SCAN(columns(at_z),'(')
            IF(strlength>0) columns(at_z) = columns(at_z)(1:strlength-1)
            !
            !Make sure that column #at_sp contains an atom species
            READ(columns(at_sp),*,ERR=170,END=170) temp
            temp = ADJUSTL(temp)
            species = temp(1:2)
            CALL ATOMNUMBER(species,snumber)
            IF( NINT(snumber) <= 0 ) THEN
              species = temp(1:1)
              CALL ATOMNUMBER(species,snumber)
            ENDIF
            IF( NINT(snumber)>0 ) THEN
              !This line indeed contains information about an atom
              NP = NP+1
            ELSE
              GOTO 170
            ENDIF
          ENDDO
          !
          170 CONTINUE
          WRITE(msg,*) 'Counted atoms, NP = ', NP
          CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
          IF(NP>0) THEN
            ALLOCATE(P(NP,4))
            P(:,:) = 0.d0
            IF( at_Wyckoff > 0 ) THEN
              WRITE(msg,*) 'Detected Wyckoff letter in column # ', at_Wyckoff
              CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
              IF( ALLOCATED(Wyckoff_letters) ) DEALLOCATE(Wyckoff_letters)
              ALLOCATE(Wyckoff_letters(NP))
              Wyckoff_letters(:) = " "
            ENDIF
            REWIND(30) !Set file pointer back to start.
                       ! -> reads everything again, except that
                       !    this time P is already allocated
          ELSE
            !Zero atom found? Something is wrong
            GOTO 800
          ENDIF
          !
        ELSE
          !P is allocated => we can proceed with atom positions and auxiliary properties
          Naux = 0
          IF( at_occ>0 ) THEN
            Naux = Naux+1
            occ = Naux
          ENDIF
          IF( at_biso>0 ) THEN
            Naux = Naux+1
            biso = Naux
          END IF
          IF( q>0 ) THEN
            Naux = Naux+1
            q = Naux
          ENDIF
          WRITE(msg,*) 'Naux = ', Naux
          CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
          IF(Naux>0) THEN
            ALLOCATE( AUX(SIZE(P,1),Naux) )
            AUX(:,:) = 0.d0
            ALLOCATE( AUXNAMES(Naux) )
            AUXNAMES(:) = ""
            IF( at_occ>0 ) THEN
              AUXNAMES(occ) = "occ"
            ENDIF
            IF( at_biso>0 ) THEN
              AUXNAMES(biso) = "biso"
            END IF
            IF( q>0 ) THEN
              AUXNAMES(q) = "q"
            END IF
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
            IF( at_Wyckoff > 0 ) THEN
              READ(columns(at_Wyckoff),*,ERR=175,END=175) Wyckoff_letters(i)
            ENDIF
            175 CONTINUE
            !Store atomic species
            species = temp(1:2)
            CALL ATOMNUMBER(species,P(i,4))
            IF( P(i,4)<=1.d-12 ) THEN
              !Failed to read an atom species from the first two characters
              !Try with only the first character
              species = temp(1:1)
              CALL ATOMNUMBER(species,P(i,4))
              !At this point, if neither the first not the two first characters contain
              !a recognizable chemical species, something is wrong and P(i,4) will be zero
              IF( q>0 .AND. q<=SIZE(AUX,1) ) THEN
                !Try to read the charge of this atom
                !(note: if this fails the atom will still have a charge equal to zero)
                temp = ADJUSTL(temp(2:))
                strlength = SCAN(temp,'+')
                IF( strlength==1 ) THEN
                  AUX(i,q) = 1.d0
                ELSEIF( strlength>1 ) THEN
                  READ(temp(1:strlength-1),*,ERR=176,END=176) tempreal
                  AUX(i,q) = tempreal
                ELSE
                  strlength = SCAN(temp,'-')
                  IF( strlength==1 ) THEN
                    AUX(i,q) = -1.d0
                  ELSEIF( strlength>1 ) THEN
                    READ(temp(1:strlength-1),*,ERR=176,END=176) tempreal
                    AUX(i,q) = -1.d0*tempreal
                  ENDIF
                ENDIF
              ENDIF
            ELSE !i.e. P(i,4)>0
              IF( q>0 .AND. q<=SIZE(AUX,1) ) THEN
                !Try to read the charge of this atom
                !(note: if this fails the atom will still have a charge equal to zero)
                temp = ADJUSTL(temp(3:))
                strlength = SCAN(temp,'+')
                IF( strlength==1 ) THEN
                  AUX(i,q) = 1.d0
                ELSEIF( strlength>1 ) THEN
                  READ(temp(1:strlength-1),*,ERR=176,END=176) tempreal
                  AUX(i,q) = tempreal
                ELSE
                  strlength = SCAN(temp,'-')
                  IF( strlength==1 ) THEN
                    AUX(i,q) = -1.d0
                  ELSEIF( strlength>1 ) THEN
                    READ(temp(1:strlength-1),*,ERR=176,END=176) tempreal
                    AUX(i,q) = -1.d0*tempreal
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
            176 CONTINUE
            !
            !
            !Store auxiliary properties if any
            IF(ALLOCATED(AUX)) THEN
              ! Naux is never used in this scope, prev. version: Naux=1
              IF(occ>0) THEN
                !Get rid of parenthesis
                strlength = SCAN(columns(at_occ),'(')
                IF(strlength>0) columns(at_occ) = columns(at_occ)(1:strlength-1)
                READ(columns(at_occ),*,ERR=800,END=800) AUX(i,occ)
              ENDIF
              IF(biso>0) THEN
                !Get rid of parethesis
                strlength = SCAN(columns(at_biso),'(')
                IF(strlength>0) columns(at_biso) = columns(at_biso)(1:strlength-1)
                strlength = SCAN(columns(at_biso),'?')
                IF(strlength>0) THEN
                  !No biso specified, set it to zero by default
                  AUX(i,biso) = 0.d0
                ELSE
                  READ(columns(at_biso),*,ERR=800,END=800) sbiso ! read to a temp variable
                  IF(at_uiso==1) THEN
                    ! translate from Uiso to Biso
                    AUX(i,biso) = sbiso * 8.d0 * pi*pi   ! * 0.789568352d+2
                  ELSE
                    ! transfer Biso directly
                    AUX(i,biso) = sbiso
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
            !
          ENDDO  !loop on atoms
          atpos_done=.TRUE. ! status: atomic positions read
          !
        ENDIF !end if ALLOCATED(P)
        !
      ENDIF ! atomic positions found.
      !
    ELSEIF( (temp(1:9)=="_symmetry" .OR. temp(1:12)=="_space_group") .AND. .NOT.symops_done) THEN
      !This is a loop over symmetry data and we still need to read it.
      ! - continue reading the header
      ! - count number of columns (Ncol) until reaching the actual data
      !
      sy_pxyz=0 ! reset the column index for symmetry operation strings
      Ncol=0 ! reset tne number of columns in this loop
      DO WHILE ( temp(1:1)=="_" ) ! header reading loop until data line found
        Ncol=Ncol+1 ! increase the number of columns
        ! Check for symmetry operation columns.
        ! The first keyword "_symmetry_equiv_pos_as_xyz" is the old and soon obsolete form.
        ! The second keyword "_space_group_symop_operation_xyz" is the new form.
        ! See: http://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Isymmetry_equiv_pos_as_xyz.html
        IF ( temp=="_symmetry_equiv_pos_as_xyz" .OR. &
           & temp=="_space_group_symop_operation_xyz" ) THEN
          sy_pxyz = Ncol ! Save the column number !
        ENDIF
        READ(30,'(a128)',ERR=200,END=200) temp ! Read the next line !
        temp = ADJUSTL(temp)
      ENDDO ! symmetry header reading
      !
      IF (sy_pxyz.NE.0) THEN ! This section contains data that we want to read.
        !
        BACKSPACE(30) ! Go back one line !
        !
        ! debug out
        WRITE(msg,*) 'Symmetry data table:'
        CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
        WRITE(msg,*) 'Ncol = ', Ncol
        CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
        WRITE(msg,*) 'Column of pos_xyz: ', sy_pxyz
        CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
        !
        !
        IF (.NOT.ALLOCATED(symops_trf)) THEN 
          WRITE(msg,*) 'symops_trf unallocated, counting number of symmetry operations'
          CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
          ! The array symsops_trf is not declared, thus the number
          ! of symmetry operations is unknown.
          ! Determine the number of symmetry operations now!
          !
          Nsym=0
          columns=""
          temp=""
          !
          ! Read in columns until an error occurs or until a new header
          ! line is found.
          DO
            ! Read lines !
            READ(30,'(a128)',ERR=180,END=180) temp
            !NOTE: a CIF file may contain symmetry operations within quotes,
            ! and with spaces after commas, e.g. '-y, x-y, z'
            j = SCAN(temp,"'")   !position of (eventual) first quote sign
            DO WHILE (j>0)
              temp(j:j) = " "
              temp = ADJUSTL(temp)
              j = SCAN(temp,"'")  !position of the second quote sign
              !Remove spaces inside the quotes
              k = SCAN(temp(1:j)," ")
              DO WHILE(k>0)
                temp = temp(1:k-1)//TRIM(temp(k+1:))
                j = SCAN(temp,"'")
                k = SCAN(temp(1:j)," ")
              ENDDO
              j = SCAN(temp,"'")
            ENDDO
            j = SCAN(temp,'"')   !position of (eventual) first quote sign
            DO WHILE (j>0)
              temp(j:j) = " "
              temp = ADJUSTL(temp)
              j = SCAN(temp,'"')  !position of the second quote sign
              !Remove spaces inside the quotes
              k = SCAN(temp(1:j)," ")
              DO WHILE(k>0)
                temp = temp(1:k-1)//TRIM(temp(k+1:))
                j = SCAN(temp,'"')
                k = SCAN(temp(1:j)," ")
              ENDDO
              j = SCAN(temp,'"')
            ENDDO
            !NOTE: a CIF file may contain slash characters (e.g. in fractions like "1/2")
            temp = ADJUSTL(temp)
            j = SCAN(temp,"/")
            DO WHILE (j>0)
              READ(temp(j-1:j-1),*,ERR=180,END=180) m
              READ(temp(j+1:j+1),*,ERR=180,END=180) n
              IF( m==1 ) THEN
                IF( n==2 ) THEN
                  temp2 = "0.5"
                ELSEIF( n==3 ) THEN
                  temp2 = "0.333"
                ELSE
                  WRITE(temp2,'(f6.3)') DBLE(m)/DBLE(n)
                ENDIF
              ELSEIF( m==3 ) THEN
                IF( n==2 ) THEN
                  temp2 = "1.5"
                ELSE
                  WRITE(temp2,'(f6.3)') DBLE(m)/DBLE(n)
                ENDIF
              ELSE
                WRITE(temp2,'(f6.3)') DBLE(m)/DBLE(n)
              ENDIF
              temp = temp(:j-2)//TRIM(ADJUSTL(temp2))//temp(j+2:)
              j = SCAN(temp,"/")
            ENDDO
            !Read the content of each column
            !READ(temp,*,ERR=178,END=178) (columns(j),j=1,Ncol)
            DO j=1,Ncol
              k = SCAN(temp," ")
              columns(j) = ADJUSTL(temp(1:k-1))
              temp = ADJUSTL(temp(k:))
            ENDDO
            ! Extract symmetry operation string !
            temp = columns(sy_pxyz)
            temp = ADJUSTL(temp)
            ! Check if this is really a symmetry operation string !
            CALL SYMOPS_CHECK_STR(TRIM(temp),k)
            ! If not, exit here !
            IF(k==0) GOTO 180
            ! Increase count of symmetry operations !
            Nsym = Nsym+1
            ! debug out
            WRITE(temp,*) Nsym 
            WRITE(msg,*) 'Symmetry operation #'//TRIM(ADJUSTL(temp))// &
                       & ': '//TRIM(ADJUSTL(columns(sy_pxyz)))
            CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
            !
          ENDDO
          !
          180 CONTINUE ! finalize counting of symmetry operations
          ! debug out
          WRITE(msg,*) 'Counted symmetry operations, Nsym = ', Nsym
          CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
          !
          IF(NSym>0) THEN
            ! There are symmetry operations to load.
            ! Allocate the transformation array in the symmetry operation module!
            ALLOCATE(symops_trf(symops_nltrf,Nsym))
            symops_trf(:,:) = 0.d0
            REWIND(30) !Set file pointer back to start.
                       ! -> reads everything again, except that
                       !    this time symops_trf is already allocated.
          ENDIF
          !
        ELSE !
          !
          WRITE(msg,*) 'symops_trf allocated, converting strings into symmetry operations'
          CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
          ! Number of symmetry operations Nsym is known and symops_trf is allocated.
          ! Clear the transformation table and read the data in now!
          !
          call SYMOPS_INIT() ! initialize transformation list
          columns=""
          temp=""
          !
          ! Read Nsym lines and transform the symmetry strings to a transformation.
          DO i=1, Nsym
            ! Read line
            READ(30,'(a128)',ERR=800,END=800) temp
            temp = ADJUSTL(temp)
            !Remove quotes
            j = SCAN(temp,"'")   !position of (eventual) first quote sign
            DO WHILE (j>0)
              temp(j:j) = " "
              temp = ADJUSTL(temp)
              j = SCAN(temp,"'")  !position of the second quote sign
              !Remove spaces inside the quotes
              k = SCAN(temp(1:j)," ")
              DO WHILE(k>0)
                temp = temp(1:k-1)//TRIM(temp(k+1:))
                j = SCAN(temp,"'")
                k = SCAN(temp(1:j)," ")
              ENDDO
              j = SCAN(temp,"'")
            ENDDO
            j = SCAN(temp,'"')   !position of (eventual) first quote sign
            DO WHILE (j>0)
              temp(j:j) = " "
              temp = ADJUSTL(temp)
              j = SCAN(temp,'"')  !position of the second quote sign
              !Remove spaces inside the quotes
              k = SCAN(temp(1:j)," ")
              DO WHILE(k>0)
                temp = temp(1:k-1)//TRIM(temp(k+1:))
                j = SCAN(temp,'"')
                k = SCAN(temp(1:j)," ")
              ENDDO
              j = SCAN(temp,'"')
            ENDDO
            !Replace divisions by actual numerical values
            j = SCAN(temp,"/")
            DO WHILE (j>0)
              READ(temp(j-1:j-1),*,ERR=800,END=800) m
              READ(temp(j+1:j+1),*,ERR=800,END=800) n
              IF( m==1 ) THEN
                IF( n==2 ) THEN
                  temp2 = "0.5"
                ELSEIF( n==3 ) THEN
                  temp2 = "0.333"
                ELSEIF( n==4 ) THEN
                  temp2 = "0.25"
                ELSE
                  WRITE(temp2,'(f9.3)') DBLE(m)/DBLE(n)
                ENDIF
              ELSEIF( m==3 ) THEN
                IF( n==2 ) THEN
                  temp2 = "1.5"
                ELSEIF( n==4 ) THEN
                  temp2 = "0.75"
                ELSE
                  WRITE(temp2,'(f9.3)') DBLE(m)/DBLE(n)
                ENDIF
              ELSE
                WRITE(temp2,'(f9.3)') DBLE(m)/DBLE(n)
              ENDIF
              temp = temp(:j-2)//TRIM(ADJUSTL(temp2))//temp(j+2:)
              j = SCAN(temp,"/")
            ENDDO
            !READ(temp,'(a)',ERR=800,END=800) (columns(j),j=1,Ncol)
            DO j=1,Ncol
              k = SCAN(temp," ")
              columns(j) = ADJUSTL(temp(1:k-1))
              temp = ADJUSTL(temp(k:))
            ENDDO
            !READ(30,*,ERR=800,END=800) (columns(j),j=1,Ncol)
            ! Get symmetry operation string
            temp = ADJUSTL(columns(sy_pxyz))
            ! Translate the string to a transformation in the
            ! symmetry operation module
            CALL SYMOPS_SET_STR(TRIM(temp),i,k)
            IF (k==0) THEN ! report interpretation error
              CALL ATOMSK_MSG(812,(/TRIM(temp)/),(/0.d0/))
              nerr = nerr+1
            END IF
            !
          ENDDO
          !
        ENDIF ! symops allocation state
        !
      ENDIF ! symmetry data reading
      !
    ENDIF !if "_atom_site" elseif "_symmetry" elseif "atom_type"
    !
  ENDIF ! end of main keyword parser "_loop"
  !
  190 CONTINUE
  !
ENDDO ! end of main reading loop
!
!
!
!
200 CONTINUE
CLOSE(30)
IF( NINT(a)==0 .AND. NINT(b)==0 .AND. NINT(c)==0 ) THEN
  !No cell parameters were found
  CALL ATOMSK_MSG(1809,(/""/),(/0.d0/))
  GOTO 1000
ENDIF
IF( .NOT.atpos_done .OR. .NOT.ALLOCATED(P) ) THEN
  !No atoms in this file
  CALL ATOMSK_MSG(1810,(/""/),(/0.d0/))
  GOTO 1000
ENDIF
!Save cell vectors in H(:,:)
CALL CONVMAT(a,b,c,alpha,beta,gamma,H)
!Find out if coordinates are reduced or cartesian
CALL FIND_IF_REDUCED(P,isreduced)
!In case of reduced coordinates, convert them to cartesian
IF(isreduced) THEN
  CALL FRAC2CART(P,H)
ENDIF
!
!Apply symmetry operations
IF ( Nsym>0 .AND. ALLOCATED(symops_trf) ) THEN
  !Symmetry operations were detected after "_symmetry_equiv_pos_as_xyz"
  !Apply them
  CALL ATOMSK_MSG(1004,(/""/),(/0.d0/))
  CALL SYMOPS_APPLY(H,P,S,AUXNAMES,AUX,0.5d0,i)
  DEALLOCATE(symops_trf)
  !
ELSE
  !No symmetry operation was defined, but a group symmetry was
  !Apply symmetry operations of that symmetry group
  IF( sgnumber > 1 .AND. sgnumber <= 230 ) THEN
    ALLOCATE(symops_trf(sgnumber,SIZE(symop_list)))
    !Initialize symmetry operations
    CALL SYMOPS_INIT()
    !Apply symmetry operations
    WRITE(msg,'(i3)') sgnumber
    CALL SG_APPLY_SYMOPS(msg,H,P,S,AUXNAMES,AUX)  
  ENDIF
  !
ENDIF
!
!Generate comments based on information gathered in the CIF file
IF (ALLOCATED(comment)) DEALLOCATE(comment)
IF( LEN_TRIM(j_name)>0 ) THEN
  !Save journal name and references in second line of comment(:)
  ALLOCATE(comment(2))
  comment(2) = TRIM(j_name)//" "//TRIM(j_volume)//" ("//TRIM(j_year)//") "//TRIM(j_page1)//"-"//TRIM(j_page2)
ELSE
  ALLOCATE(comment(1))
ENDIF
!In first line of comment(:), create a nice comment based on chemical name and formula
 comment(1) = ""
IF( LEN_TRIM(chemical_name) > 0 ) THEN
  comment(1) = TRIM(chemical_name)//" -"
ENDIF
 comment(1) = TRIM(comment(1))//" "//TRIM(chemical_formula)
 comment(1) = ADJUSTL(comment(1))
!
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
!
!
END MODULE in_cif