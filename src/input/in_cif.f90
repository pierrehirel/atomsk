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
!* Last modification: P. Hirel - 05 Oct. 2015                                     *
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
CHARACTER(LEN=128):: msg, temp, temp2
CHARACTER(LEN=32),DIMENSION(32):: columns !contents of columns
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: atpos_done, symops_done, isreduced
INTEGER:: at_occ
INTEGER:: at_biso
INTEGER:: at_uiso
INTEGER:: at_sp, at_x, at_y, at_z !position of atom species and coordinates in the line
INTEGER:: sy_pxyz !position of symmetry operations string in the line
INTEGER:: i, j, iaux, k
INTEGER:: Naux !number of auxiliary properties
INTEGER:: Ncol, NP, sp_NP
INTEGER:: Nsym !number of symmetry operations
INTEGER:: occ  !index of occupancies in AUXNAMES
INTEGER:: biso !index of Biso in AUXNAMES
INTEGER:: strlength
REAL(dp):: a, b, c, alpha, beta, gamma !supercell (conventional notation)
REAL(dp):: snumber, sbiso
REAL(dp),DIMENSION(3,3),INTENT(OUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
atpos_done=.FALSE.
symops_done=.FALSE.
at_sp=0
at_x=0
at_y=0
at_z=0
at_occ=0
at_biso=0
at_uiso=0
sy_pxyz=0
Naux=0
NP=0
Nsym=0
occ=0
biso=0
sp_NP=9
a=0.d0
b=0.d0
c=0.d0
alpha=0.d0
beta=0.d0
gamma=0.d0
IF (ALLOCATED(symops_trf)) DEALLOCATE(symops_trf)
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
    !
    ! Note: In order to support symmetry operations, where the number
    !       of listed atomic sites may differ from the chemical formula,
    !       I have removed (commented out) the code, which determines
    !       the number of atoms to read from the chemical formula.
    !       The new code will still read and report the chemical
    !       formula, but will not determine NP and will not allocate P
    !       at this point.
    !       Instead, the number of listed atoms is always determined
    !       from the "_atom" loops below, which will be read twice. On
    !       the first run, it will determine NP and allocate P, while on
    !       the second run, it will read the atomic site data into P.
    !       At the end, the symmetry operations are applied.
    !       As a further consistency check, we could move this code
    !       to the end of the reading routine, in order to check whether
    !       the composition obtained after reading the "_atom" list and
    !       applying the symmetry operations corresponds to the chemical
    !       formula. However, I have not done this.
    !
    !       J. Barthel, 2015-07-31
    !
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
    !
  ELSEIF( temp(1:5)=="loop_" ) THEN
    !
    !Check if it is a loop over atom positions or symmetry operations
    READ(30,'(a128)',ERR=200,END=200) temp
    temp = ADJUSTL(temp)
    IF( temp(1:5)=="_atom" .AND. .NOT.atpos_done) THEN
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
            READ(columns(at_sp),*,ERR=170,END=170) temp
            temp = ADJUSTL(temp)
            species = temp(1:2)
            CALL ATOMNUMBER(species,snumber)
            IF( snumber<=0.d0 ) THEN
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
          WRITE(msg,*) 'Naux = ', Naux
          CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
          IF(Naux>0) THEN
            ALLOCATE( AUX(SIZE(P,1),Naux) )
            AUX(:,:) = 0.d0
            ALLOCATE( AUXNAMES(Naux) )
            ! keep Naux in the new version, prev. version: Naux=1
            iaux = 1
            IF( at_occ>0 ) THEN
              occ = iaux
              iaux = iaux+1
              AUXNAMES(occ) = "occ"
            ENDIF
            IF( at_biso>0 ) THEN
              biso = iaux
              iaux = iaux+1
              AUXNAMES(biso) = "biso"
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
              ! Naux is never used in this scope, prev. version: Naux=1
              IF(occ>0) THEN
                !Get rid of parethesis
                strlength = SCAN(columns(at_occ),'(')
                IF(strlength>0) columns(at_occ) = columns(at_occ)(1:strlength-1)
                READ(columns(at_occ),*,ERR=800,END=800) AUX(i,occ)
              ENDIF
              IF(biso>0) THEN
                !Get rid of parethesis
                strlength = SCAN(columns(at_biso),'(')
                IF(strlength>0) columns(at_biso) = columns(at_biso)(1:strlength-1)
                READ(columns(at_biso),*,ERR=800,END=800) sbiso ! read to a temp variable
                IF(at_uiso==1) THEN
                  ! translate from Uiso to Biso
                  AUX(i,biso) = sbiso * 0.789568352d+2 ! * 8 * Pi**2
                ELSE
                  ! transfer Biso directly
                  AUX(i,biso) = sbiso
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
    ELSEIF( temp(1:9)=="_symmetry" .AND. .NOT.symops_done) THEN
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
            READ(30,*,ERR=180,END=180) (columns(j),j=1,Ncol)
            ! Extract symmetry operation string !
            temp = ADJUSTL(columns(sy_pxyz))
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
            symops_trf = 0.d0
            REWIND(30) !Set file pointer back to start.
                       ! -> reads everything again, except that
                       !    this time symops_trf is already allocated.
          ENDIF
          !
        ELSE !
          !
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
            READ(30,*,ERR=800,END=800) (columns(j),j=1,Ncol)
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
    ENDIF !if "_atom" elseif "_symmetry"
    !
  ENDIF ! end of main keyword parser
  !
ENDDO ! end of main reading loop
!
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
!Apply symmetry operations
IF (Nsym>0.AND.ALLOCATED(symops_trf)) THEN
  CALL SYMOPS_APPLY(H,P,AUXNAMES,AUX,0.5d0,i)
  NP=SIZE(P,1) ! update number of atom sites
ENDIF
!
IF (ALLOCATED(symops_trf)) DEALLOCATE(symops_trf)
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
