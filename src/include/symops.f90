MODULE symops
!
!
!**********************************************************************************
!*  SYMOPS                                                                        *
!**********************************************************************************
!* This module handles symmetry operations of linear type.                        *
!* A symmetry operation transforms a 3D position P into a new position P' with    *
!*     P' = M.P + S,                                                              *
!* where M is a 3x3 matrix and S is a 3D shift vector.                            *
!* Internally a symmetry operation is defined by a 12*REAL(dp) array, where       *
!* the items 1-3 define S and the items 4-12 define M, with 4-6 for the x',       *
!* 7-9 for the y', and 10-12 for the z' component.                                *
!**********************************************************************************
!* (C) July 2015 - Juri Barthel                                                   *
!*     Gemeinschaftslabor fuer Elektronenmikroskopie                              *
!*     RWTH Aachen (GERMANY)                                                      *
!*     ju.barthel@fz-juelich.de                                                   *
!* Last modification: P. Hirel - 22 March 2017                                    *
!**********************************************************************************
!* Symmetry operation handling and parsing, added by J. Barthel, July 2015        *
!* SYMOPS_INIT    initializes the array symops_trf to identity                    *
!*                operations. Since we leave the allocation of the array to       *
!*                other routines, and keep the array for public access accross    *
!*                the program, this routine is a pure initialization routine.     *
!*                Note:  The allocation state of symops_trf will be checked,      *
!*                but not altered. No error message will occur if not allocated.  *
!* SYMOPS_APPLY   applies the currently defined symmetry operations               *
!*                from the array "symops_trf" to the passed arrays P and AUX.     *
!*                P(:,1:3) are assumed to be cartesian atom coordinates.          *
!*                Note: P and AUX may be re-allocated at the end of this routine  *
!*                since the number of atoms may have increased. Update size       *
!*                variables in the calling routine after calling SYMOPS_APPLY!    *
!* SYMOPS_CHECK_STR checks if a string defines a symmetry operation               *
!* SYMOPS_SET_STR translates a string symmetry operation to numbers defining      *
!*                linear transformations and stores these numbers in symops       *
!* SYMOPS_PARSE_STR_LINTRF translates a 1D operation string to linear             *
!*                transformation parameters by parsing the string.                *
!* SYMOPS_SET_SGNAME sets symmetry operations for a spacegroup given by the       *
!*                space group Hermann-Mauguin symbol.                             *
!* SYMOPS_SET_SGNUM sets symmetry operations for a spacegroup given by the        *
!*                space group number (1 - 230).                                   *
!* SG_APPLY_SYMOPS applies the symmetry operations of the given space group.      *
!* Remark: SYMOPS_PARSE_STR_LINTRF is a recursive routine.                        *
!*                Use the respective compiler option to enable recursice routines *
!*                (ifort: /recursive)                                             *
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
USE spacegroups
!
IMPLICIT NONE
!
!
! MODULE DATA
!
CHARACTER(LEN=3),PARAMETER,PUBLIC:: symops_chanstr = 'xyz'
INTEGER,PARAMETER,PUBLIC:: symops_nltrf = 12
REAL(dp),DIMENSION(:,:),ALLOCATABLE,PUBLIC:: symops_trf ! symmetry operations list
!
!
!
! MODULE ROUTINES
!
CONTAINS
!
!
SUBROUTINE SYMOPS_INIT()
!
INTEGER:: n, i, j
!
! Initialization
IF(.NOT.ALLOCATED(symops_trf)) RETURN ! nothing to do, exit
j=SIZE(symops_trf,1) ! Get number of row entries
n=SIZE(symops_trf,2) ! Get number of transformations
!
! Check number of row entries for consistency:
IF(j.NE.symops_nltrf) RETURN ! inconsistent row size, exit
                             ! Leave initialization to external routine.
IF(n<=0) RETURN ! Invalid number of rows, exit
                ! Leave initialization to external routine.
symops_trf = 0.d0 ! set all entries to zero
DO i=1, n ! loop through rows and set M to Identity
  symops_trf( 4,i) = 1.d0
  symops_trf( 8,i) = 1.d0
  symops_trf(12,i) = 1.d0
ENDDO
!
END SUBROUTINE SYMOPS_INIT
!
!
!
!
SUBROUTINE SYMOPS_APPLY(H,P,S,AUXNAMES,AUX,dmindist,nchk)
!
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE::P, S, AUX
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
REAL(DP),INTENT(IN) :: dmindist ! minimum atom distance
INTEGER,INTENT(OUT) :: nchk ! success code: 0=failure, 1=success
!
CHARACTER(LEN=128) :: temp, msg ! strings
INTEGER::i,j,isym,occ,idup,j3
INTEGER::n1,n2
INTEGER::MP,NP,NS,Naux,Nsym
REAL(dp)::a,b,c,alpha,beta,gamma
REAL(dp)::trf(symops_nltrf),rtmp,rfdif(3),rdif(3),rdist,p1(3),p2(3),ds
REAL(dp)::socc1,socc2,soccs
INTEGER,DIMENSION(:),ALLOCATABLE::SUSE ! intermediate array for new data
REAL(dp),DIMENSION(:,:),ALLOCATABLE::S1,S2,R1,R2,SAUX1,SAUX2 ! intermediate arrays for new data 
REAL(dp),DIMENSION(3,3):: G ! inverse of H
!
WRITE(msg,*) 'entering SYMOPS_APPLY, Nsym = ', Nsym
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Initial checks.
nchk=0
NP=0
NS=0
IF (.NOT.ALLOCATED(P)) RETURN ! Just exit, since there are no atoms.
IF (.NOT.ALLOCATED(symops_trf)) RETURN ! Just exit, since no symmetries are defined.
!
!Initialization
CALL MATCONV(H,a,b,c,alpha,beta,gamma)
CALL INVMAT(H,G)
Nsym=SIZE(symops_trf,2)
NP=SIZE(P,1)
MP=SIZE(P,2)
IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
  NS=SIZE(S,1)
ENDIF
!
IF( verbosity==4 ) THEN
  WRITE(msg,*) NP
  msg = '- input number of atoms: '//TRIM(ADJUSTL(msg))
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  WRITE(msg,*) NS
  msg = '- input number of shells: '//TRIM(ADJUSTL(msg))
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ENDIF
!
IF( Nsym<=0 .OR. NP<=0 ) THEN
  !No atom or no symmetry operation to apply => exit
  RETURN
ENDIF
!
Naux=0
idup=0
IF (ALLOCATED(AUX)) Naux=SIZE(AUX,1) ! Get number of auxiliary data
IF (Naux>0 .AND. Naux.NE.NP) THEN ! P and AUX should have the same length
  CALL ATOMSK_MSG(1802,(/'AUX'/),(/0.d0/))
  nerr = nerr+1
  RETURN
ENDIF
IF (Naux>0) Naux=SIZE(AUX,2)
occ=0
IF (Naux>0) THEN
  DO i=1, Naux
    IF(TRIM(ADJUSTL(AUXNAMES(i)))=='occ') THEN
      occ=i
      EXIT
    ENDIF
  ENDDO
ENDIF
socc1=1.d0
socc2=1.d0
soccs=0.d0
n1=NP
n2=2*n1
!Prepare initial arrays
ALLOCATE(S1(n1,MP),S2(n2,MP),SUSE(n2))
S1=P ! copy P to S1
CALL CART2FRAC(S1,H) ! transform to fractional coordinates
S1(1:n1,1:3)=MODULO(S1(1:n1,1:3),1.0) ! wrap into [0,1[
S2=0.d0 ! reset S2
S2(1:n1,1:MP)=S1(1:n1,1:MP) ! fill first half of S2 with S1

IF( NS>0 ) THEN
  !Deal with shells (in sense of ionic core-shell model)
  ALLOCATE(R1(n1,MP),R2(n2,MP))
  R1=S
  CALL CART2FRAC(R1,H) ! transform to fractional coordinates
  R1(1:n1,1:3)=MODULO(R1(1:n1,1:3),1.0) ! wrap into [0,1[
  R2=0.d0
  R2(1:n1,1:MP)=R1(1:n1,1:MP) ! fill first half of R2 with R1
ENDIF

SUSE=0 ! reset SUSE
SUSE(1:n1)=1 ! mark the copied content of P to be used
IF (Naux>0) THEN ! there are auciliaries ... 
  ALLOCATE(SAUX1(n1,Naux),SAUX2(n2,Naux)) ! ... allocate arrays to handle them
  SAUX1=AUX ! copy AUX to SAUX1
  SAUX2=0.d0 ! reset SAUX2
  SAUX2(1:n1,1:Naux)=SAUX1(1:n1,1:Naux) ! fill the first half of SAUX2 with SAUX1
ENDIF
!
DO isym=1,Nsym ! apply each symmetry op. ...
  !Get current transformation parameters
  trf(:) = symops_trf(:,isym)
  WRITE(msg,'(a9,i3,a3,96e9.3)') "Sym.op. #", isym, " : ", trf(:)
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !Apply this transformation from S1 to S2
  DO i=1,n1 ! ... to each of the current positions
    ! get current atomic position in fractional coordinates
    p1(1:3) = MODULO( S1(i,1:3), 1.0 ) ! wrapped into the unit cell
    DO j=1,3
      j3=3*(j-1)
      ! get new fractional coordinate from symmetry operation
      rtmp = trf(j)+trf(4+j3)*p1(1)+trf(5+j3)*p1(2)+trf(6+j3)*p1(3)
      ! wrap it periodically into the cell
      S2(i+n1,j) = MODULO( rtmp, 1.d0 )
    ENDDO
    IF (MP>3) S2(i+n1,4:MP) = S1(i,4:MP) ! copy the rest of the data row (at 4 we have the atomic numbers)
    IF (NS>0) THEN
      R2(i+n1,4:MP) = R1(i,4:MP)
    ENDIF
    IF (Naux>0) THEN
      SAUX2(i+n1,1:Naux) = SAUX1(i,1:Naux) ! ... and copy the aux data
    ENDIF
    !
    ! Now that we have a new atomic site, we should check if this atom is
    ! a duplicate, which is absolutely identical to one of the old atoms.
    ! Identical means S2(i,1:4)==S2(j,1:4), including the atomic number.
    ! We will not allow duplicates here and leave them unmarked (0) in the
    !  array "SUSE". We will mark all unique new positions (1) in "SUSE".
    !
    ! Note on partial site occupancies:
    ! The above check allows partial site occupancies by different atomic
    ! species, since their atomic number is different. It doesn't allow a
    ! site sharing by two identical species, regardless of their occupancy
    ! factor. This would anyway be a very strange structure definition.
    !
    idup=0 ! reset the duplicate indicator
    p1 = S2(i+n1,1:3) ! get the new symmetry position
    DO j=1,n1+i-1 ! look through the S2 data up to the new position
      IF (SUSE(j)==0) CYCLE ! skip any unmarked previous position.
      p2 = S2(j,1:3) ! get fractional coordinate (should already be in [0,1[ )
      rfdif = MODULO(p2-p1,1.0) ! fractional distance vector (modulo unit cell)
      ! Note: We are looking for the shortest fractional distance
      !       under periodic boundary conditions. This may require to look
      !       for shorter vectors across the cell boundary:
      IF (rfdif(1)>0.5d0) rfdif(1) = rfdif(1)-1.0d0
      IF (rfdif(2)>0.5d0) rfdif(2) = rfdif(2)-1.0d0
      IF (rfdif(3)>0.5d0) rfdif(3) = rfdif(3)-1.0d0
      ! calculate cartesian distance vector
      rdif(1) = SUM( H(1:3,1)*rfdif )
      rdif(2) = SUM( H(1:3,2)*rfdif )
      rdif(3) = SUM( H(1:3,3)*rfdif )
      rdist = VECLENGTH(rdif) ! get distance in A
      ! check if this is a "duplicate"
      ! ( distance to other atom < min. distance and atomic numbers equal )
      IF (rdist<=dmindist .AND. DABS(S2(j,4)-S2(i,4))<0.5d0) THEN
        idup = idup + 1 ! raise the duplicate count
        EXIT ! Exit the loop over other atoms, we found one duplicate, that's enough.
      ENDIF
    ENDDO
    ! handle non-duplicates
    IF (idup==0) SUSE(i+n1)=1 ! mark new position as unique
    !
  ENDDO ! (positions)
  !
  ! There is a new set of positions in S2.
  ! We transfer it now back to S1 and reinitiate everything.
  ! Update number of atoms n1
  n1 = SUM(SUSE)
  IF (ALLOCATED(S1)) DEALLOCATE(S1)
  IF ( NS>0 .AND. ALLOCATED(R1) ) DEALLOCATE(R1)
  IF (Naux>0 .AND. ALLOCATED(SAUX1)) DEALLOCATE(SAUX1)
  ALLOCATE(S1(1:n1,1:MP))
  IF (NS>0) ALLOCATE(R1(1:n1,1:MP))
  IF (Naux>0) ALLOCATE(SAUX1(1:n1,1:Naux))
  ! Transfer the usable data from S2 to the new S1
  j = 0
  DO i=1, n2
    IF (SUSE(i)==0) CYCLE ! skip unused atoms
    j = j + 1 ! raise atom list index
    S1(j,1:MP) = S2(i,1:MP) ! positions and species
    IF (NS>0) THEN
      R1(j,1:MP) = R2(i,1:MP) ! positions and species
    ENDIF
    IF (Naux>0) THEN ! auxiliaries
      SAUX1(j,1:Naux) = SAUX2(i,1:Naux)
    ENDIF
  ENDDO
  ! Update the secondary array set
  n2 = 2*n1 ! possible number of new atoms
  IF (ALLOCATED(S2)) DEALLOCATE(S2)
  IF (ALLOCATED(SUSE)) DEALLOCATE(SUSE)
  IF ( NS>0 .AND. ALLOCATED(R2) ) DEALLOCATE(R2)
  IF (Naux>0 .AND. ALLOCATED(SAUX2)) DEALLOCATE(SAUX2)
  ALLOCATE(S2(1:n2,1:MP),SUSE(1:n2)) ! new S2 and SUSE
  S2 = 0.d0
  SUSE = 0
  S2(1:n1,1:MP) = S1(1:n1,1:MP) ! copy S1
  SUSE(1:n1) = 1 ! mark S1 as usable
  IF (NS>0) THEN
    ALLOCATE(R2(1:n2,1:MP))
    R2 = 0.d0
    R2(1:n1,1:MP) = R1(1:n1,1:MP) ! copy R1
  ENDIF
  IF (Naux>0) THEN
    ALLOCATE(SAUX2(1:n2,1:Naux)) ! new SAUX2
    SAUX2 = 0.d0
    SAUX2(1:n1,1:Naux) = SAUX1(1:n1,1:Naux)! copy SAUX1
  ENDIF
  !
ENDDO ! (symmetry operations)
!
! Finally update P and AUX
DEALLOCATE(P)
IF (Naux>0 .AND. ALLOCATED(AUX)) DEALLOCATE(AUX)
NP = n1 ! set new number of atoms
ALLOCATE(P(1:NP,1:MP))
P = S1 ! set new atomic site data
CALL FRAC2CART(P,H) ! transform P to cartesian coordinates
IF (NS>0) THEN
  IF(ALLOCATED(S)) DEALLOCATE(S)
  NS = n1 ! set new number of atoms
  ALLOCATE(S(1:NP,1:MP))
  S = R1 ! set new atomic site data
  CALL FRAC2CART(S,H) ! transform P to cartesian coordinates
ENDIF
IF (Naux>0) THEN
  ALLOCATE(AUX(1:NP,1:Naux))
  AUX = SAUX1 ! set new auxiliary data
ENDIF
! Clean up the heap.
IF (ALLOCATED(S1)) DEALLOCATE(S1)
IF (ALLOCATED(S2)) DEALLOCATE(S2)
IF (ALLOCATED(R1)) DEALLOCATE(R1)
IF (ALLOCATED(R2)) DEALLOCATE(R2)
IF (ALLOCATED(SUSE)) DEALLOCATE(SUSE)
IF (Naux>0 .AND. ALLOCATED(SAUX2)) DEALLOCATE(SAUX2)
WRITE(msg,*) NP, NS
msg = '- output number of atoms, shells: '//TRIM(ADJUSTL(msg))
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
nchk=1
!
END SUBROUTINE SYMOPS_APPLY
!
!
!
!
SUBROUTINE SYMOPS_CHECK_STR(instr,nchk)
!
INTEGER,PARAMETER:: numtests = 2 ! NUMBER OF CHECKS: INCREASE IF YOU IMPLEMENT MORE!
CHARACTER(LEN=*), INTENT(IN) :: instr ! input string
INTEGER, INTENT(OUT) :: nchk ! output check result:
                             !   0: no symmetry operation
                             !   1: symmetry operation
INTEGER:: i, l, c1, c2, ichk ! iterators, helpers, and subcheck results
CHARACTER(LEN=128) :: temp, msg, symopchars ! strings
!
! Initialize
nchk=0
ichk=0
i=0
c1=0
c2=0
temp=ADJUSTL(instr)
l=LEN_TRIM(temp)
symopchars='xyz1234567890+-*/.'
!
msg = 'entering SYMOPS_CHECK_STR'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
msg = "checking '"//TRIM(temp)//"' for CIF symmetry operation style"
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Check #1: A well defined operation contains the right sides
!          of 3 equations, separated by extacly 2 "," characters
c1=INDEX(TRIM(temp(1:l)),",",BACK=.FALSE.)
c2=INDEX(TRIM(temp(1:l)),",",BACK=.TRUE.)
IF (c1>0 .AND. c2>0 .AND. c2>c1) ichk = ichk + 1
IF (ichk==0) THEN
  msg = "- failed the 2-comma test."
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  RETURN
END IF
!Remove the commas.
temp(c1:c1) = " "
temp(c2:c2) = " "
!
!Check #2: After removing the 2 commas, the string should contain
!          only the characters in symopchars.
!Replace all the occurrences of symopchars in temp by " ":
c1=0
DO i=1, l
  c1=INDEX(TRIM(symopchars),temp(i:i))
  IF (c1>0) temp(i:i) = " "
END DO
c2=LEN_TRIM(temp(1:l))
IF (c2>0) THEN
  msg = "- failed the restricted characters test."
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  RETURN
ENDIF
ichk = ichk+1 ! 2nd check passed.
!
IF (ichk==numtests) THEN
  nchk=1 ! success, this is a symmetry operation
  temp=ADJUSTL(instr)
  msg = "- '"//TRIM(temp)// &
      & "' passed all tests for a CIF symmetry op."
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
ENDIF
!
END SUBROUTINE SYMOPS_CHECK_STR
!
!
!
!
SUBROUTINE SYMOPS_SET_STR(instr,irow,nchk)
! This routine transforms instr to linear transformation parameters
! of the form M.(x,y,z)+S, where S is a 3D shift vector and and M
! is a 3x3 tranformation matrix.
! The values of S and M are stored in sequence (S,M) in the array
! symops_trf(:,irow) of the module symops (link 'symops.f90'!).
! symops_trf must be preallocated and initialized.
! 'instr' should be of the form 'x-y, y+1/2, -z/2', which would result
! in S = (/ 0, 0.5, 0 /) and
! and M = (/ (/ 1, -1, 0 /), (/ 0, 1, 0 /), (/ 0, 0, -0.5 /) /).
! The input parameter 'irow' determines the row of the data storage.
! The output parameter nchk signalizes the success or failur of the
! translation in this routine.
!
CHARACTER(LEN=*),INTENT(IN):: instr ! input string
INTEGER,INTENT(IN):: irow ! input row index in the transformation table
INTEGER, INTENT(OUT):: nchk ! output success code: 0: failure, 1: success
!
CHARACTER(LEN=128) :: temp, msg, substr(3) ! strings
INTEGER:: i, j, k, l ! iterators
INTEGER:: c1, c2 ! positions of substring separation
REAL(dp):: trf(symops_nltrf) ! local transformation
!
! Initialize
trf = 0.d0
trf(4)  = 1.d0
trf(8)  = 1.d0
trf(12) = 1.d0
substr=""
 c1=0
 c2=0
i=0
j=0
k=0
temp=ADJUSTL(instr)
l=LEN_TRIM(temp)
nchk=0
IF (.NOT.ALLOCATED(symops_trf)) GOTO 850
IF (SIZE(symops_trf,1).NE.symops_nltrf) GOTO 850
IF (SIZE(symops_trf,2)<irow) GOTO 850
! The minimum valid string is 'x,y,z'.
IF (l<5) GOTO 801 ! warn & report invalid CIF symmetry operation string.
!
msg = 'entering SYMOPS_SET_STR: '//TRIM(instr)
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
! Determine the two comma positions
 c1=INDEX(temp(1:l),",",BACK=.FALSE.)
 c2=INDEX(temp(1:l),",",BACK=.TRUE.)
! There should be two commas in the string. The first comma should be at
! a position >1, such that the first substring fits in before it. The
! second comma should come at least 2 characters beyond the first comma,
! such that the second substring fits in between. The second comma
! shouldn't be at the end of the string, such that there is still space
! for the third substring.
IF (c1<2 .OR. c2-c1<2 .or. c2>=l) GOTO 801
! Extract the 3 substrings
substr(1) = TRIM(ADJUSTL(temp(1   :c1-1)))
substr(2) = TRIM(ADJUSTL(temp(c1+1:c2-1)))
substr(3) = TRIM(ADJUSTL(temp(c2+1:l   )))
! Parsing the 3 substrings
temp = TRIM(substr(1))
call SYMOPS_PARSE_STR_LINTRF(temp,trf(1),trf(4:6),k)
IF (k==0) GOTO 801
temp = TRIM(substr(2))
call SYMOPS_PARSE_STR_LINTRF(temp,trf(2),trf(7:9),k)
IF (k==0) GOTO 801
temp = TRIM(substr(3))
call SYMOPS_PARSE_STR_LINTRF(temp,trf(3),trf(10:12),k)
IF (k==0) GOTO 801
!
! Store the transformation in the symops module
symops_trf(:,irow) = trf(:)
nchk=1 ! success
! debug out
msg = "- found transformation for '"//TRIM(ADJUSTL(instr))//"':"
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,901) trf(1), trf(4), trf(5), trf(6)
msg = "  x' = "//TRIM(ADJUSTL(msg))
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,901) trf(2), trf(7), trf(8), trf(9)
msg = "  y' = "//TRIM(ADJUSTL(msg))
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
WRITE(msg,901) trf(3), trf(10), trf(11), trf(12)
msg = "  z' = "//TRIM(ADJUSTL(msg))
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
901 FORMAT(F8.3," + ",F8.3,"*x + ",F8.3,"*y + ",F8.3,"*z")

GOTO 1000 ! skip warnings / errors.
!
801 CONTINUE
CALL ATOMSK_MSG(1707,(/TRIM(temp)/),(/0.d0/))
nwarn = nwarn+1
GOTO 1000
850 CONTINUE
!Unable to store transformation in symops module => error message
CALL ATOMSK_MSG(1707,(/msg/),(/0.d0/))
!
1000 CONTINUE
!
END SUBROUTINE SYMOPS_SET_STR
!
!
!
!
RECURSIVE SUBROUTINE SYMOPS_PARSE_STR_LINTRF(instr,shift,slope,nchk)
!
CHARACTER(LEN=*),INTENT(IN):: instr        !the string to parse
REAL(dp),INTENT(OUT):: shift               !Shift along X, Y, Z
REAL(dp),DIMENSION(3),INTENT(OUT):: slope  !Multiplication factor along X, Y, Z
INTEGER,INTENT(OUT)::nchk
!
INTEGER:: i,j,l
CHARACTER(LEN=128):: temp,msg
INTEGER:: ichan, k1, k2
REAL(dp):: rtmp, sh1, sh2
REAL(dp),DIMENSION(3):: sl1, sl2
!
nchk=0
shift=0.d0
slope=0.d0
temp=ADJUSTL(instr)
l=LEN_TRIM(temp)
!
msg = 'entering SYMOPS_PARSE_STR_LINTRF'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
msg = "- parsing '"//temp(1:l)//"'..."
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
IF (l<=2) GOTO 400 ! TRIVIAL CASE: numbers or characters, no operation
!
IF(temp(1:1)=="+") THEN
  temp = temp(2:)
ENDIF
! Check for basic operations: the operation character should be beyond
! position 1 of the string "temp". Store the position in "j".
j=INDEX(temp(1:l),'+')
IF (j>1) GOTO 500
j=INDEX(temp(1:l),'-')
IF (j>1) GOTO 520
j=INDEX(temp(1:l),'-',BACK=.TRUE.) !necessary if first character is "-", e.g. "-z-0.5"
IF (j>1) GOTO 520
j=INDEX(temp(1:l),'*')
IF (j>1) GOTO 540
j=INDEX(temp(1:l),'/')
IF (j>1) GOTO 560
! No operation character found. We assume now, that temp is a trivial
! number or character, which would mean something like "0.5" or "-10"
! since temp is longer than 2.
!
400 CONTINUE ! TRIVIAL CASE: single numbers or single characters, no operation
!                          ! Solve it here !
msg = "- assuming TRIVIAL CASE: "//TRIM(temp(1:l))
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
! Determine the output channel: shift=0, x=1, y=2, z=3
ichan=0
j=0
DO i=1,3
  j=INDEX(temp(1:l),symops_chanstr(i:i))
  IF(j>0) THEN
    ichan=i ! got the character channel
    temp(j:j)='1' ! replace the character by a number
    EXIT ! done, exit this loop
  ENDIF
END DO
! "temp" should now contain just a trivial number
! Translate to a numerical value
READ(temp,*,ERR=801,END=801) rtmp
! Put numerical value into channel
IF(ichan==0) THEN
  shift = rtmp
ELSE
  slope(ichan) = rtmp
ENDIF
! Success of trivial case
GOTO 800
!
500 CONTINUE ! ADDITION CASE       : operation '+' with two substrings
!            ! Bifurcate here !
msg = "- assuming ADDITION: "//TRIM(temp(1:l))
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
call SYMOPS_PARSE_STR_LINTRF(ADJUSTL(temp(1:j-1)),sh1,sl1,k1)
call SYMOPS_PARSE_STR_LINTRF(ADJUSTL(temp(j+1:l)),sh2,sl2,k2)
IF (k1==0 .OR. k2==0) GOTO 802
! Combine
shift = sh1+sh2
slope = sl1+sl2
GOTO 800
!
520 CONTINUE ! SUBTRACTION CASE    : operation '-' with two substrings
!            ! Bifurcate here !
msg = "- assuming SUBTRACTION: "//TRIM(temp(1:l))
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
call SYMOPS_PARSE_STR_LINTRF(ADJUSTL(temp(1:j-1)),sh1,sl1,k1)
call SYMOPS_PARSE_STR_LINTRF(ADJUSTL(temp(j+1:l)),sh2,sl2,k2)
IF (k1==0 .OR. k2==0) GOTO 802
! Combine
shift = sh1-sh2
slope = sl1-sl2
GOTO 800
!
540 CONTINUE ! MULTIPLICATION CASE : operation '*' with two substrings
!            ! Bifurcate here !
msg = "- assuming MULTIPLICATION: "//TRIM(temp(1:l))
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
call SYMOPS_PARSE_STR_LINTRF(ADJUSTL(temp(1:j-1)),sh1,sl1,k1)
call SYMOPS_PARSE_STR_LINTRF(ADJUSTL(temp(j+1:l)),sh2,sl2,k2)
IF (k1==0 .OR. k2==0) GOTO 802
! Combine: Possible scenarios
! 1) S1=A,M1=0, S2=0,M2=X -> S=0, M=A*X -> S=S1*S2, M=S1*M2+S2*M1
! 2) S1=0,M1=X, S2=A,M2=0 -> S=0, M=X*A -> S=S1*S2, M=S1*M2+S2*M1
! 3) S1=A,M1=0, S2=B,M2=0 -> S=A*B, M=0 -> S=S1*S2, M=S1*M2+S2*M1
! - no scenario where M1/=0 and M2/=0, this would be a non-linear form
! - no scenario like A*(B+X), since brackets are not supported.
shift = sh1*sh2
slope = sh1*sl2+sh2*sl1
GOTO 800
!
560 CONTINUE ! DIVISION CASE       : operation '/' with two substrings
!            ! Bifurcate here !
msg = "- assuming DIVISION: "//TRIM(temp(1:l))
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
call SYMOPS_PARSE_STR_LINTRF(ADJUSTL(temp(1:j-1)),sh1,sl1,k1)
call SYMOPS_PARSE_STR_LINTRF(ADJUSTL(temp(j+1:l)),sh2,sl2,k2)
IF (k1==0 .OR. k2==0) GOTO 802
! Combine: Possible scenarios
! 1) S1=0,M1=X, S2=A,M2=0 -> S=0, M=X/A -> S=S1/S2, M=M1/S2
! 2) S1=B,M1=0, S2=A,M2=0 -> S=B/A, M=0 -> S=S1/S2, M=M1/S2
! - requires S2/=0 and M2==0 ! Only division by pure numbers
! - no scenario with A/M2 or M1/M2, this would be a non-linear form
IF (sh2==0.d0) GOTO 802 ! error, division by zero
IF ( sl2(1)/=0.d0 .OR. sl2(2)/=0.d0 .OR. sl2(3)/=0.d0 ) GOTO 802 ! error, division by variable
shift = sh1/sh2
slope = sl1/sh2
GOTO 800
!
!
800 CONTINUE ! success
nchk = 1
! debug out
WRITE(msg,'(F8.3," + ",F8.3,"*x + ",F8.3,"*y + ",F8.3,"*z")') &
     &    shift, slope(1), slope(2), slope(3)
CALL ATOMSK_MSG(999,(/"Op.string"//msg/),(/0.d0/))
msg = "- transformation: "//TRIM(ADJUSTL(msg))
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
GOTO 1000
801 CONTINUE ! error, failed to convert to number
msg = "- FAILED TO CONVERT: "//TRIM(temp(1:l))
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
CALL ATOMSK_MSG(808,(/temp(1:l)/),(/0.d0/))
nerr = nerr+1
GOTO 1000
802 CONTINUE ! error, failed to convert operation
msg = "- FAILED TO CONVERT: "//TRIM(temp(1:l))
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
CALL ATOMSK_MSG(1813,(/temp(1:l)/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
1000 CONTINUE
!
END SUBROUTINE SYMOPS_PARSE_STR_LINTRF
!
!
!
!
SUBROUTINE SYMOPS_SET_SGNAME(sgname,nchk)
!
CHARACTER(LEN=*),INTENT(IN)::sgname ! Hermann-Mauguin symbol identifying the space group
INTEGER,INTENT(OUT)::nchk ! output success code: 0: failure, 1: success
!
INTEGER::nsgnum,i
CHARACTER(LEN=128)::temp,msg
!
!Initialization
nchk=0
nsgnum=0
i=0
msg = 'entering SYMOPS_SET_SGNAME'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Get the space group number. This call will initialize the spacegroup
!module if it isn't initialized yet.
temp=ADJUSTL(sgname)
CALL SG_NAMGETNUM(TRIM(temp),nsgnum)
IF (nsgnum<1 .OR. nsgnum>sg_nummax) GOTO 801 ! invalid space group name
WRITE(msg,*) nsgnum
msg = "Identified space group '"//TRIM(temp)//"' as number "// &
    & TRIM(ADJUSTL(msg))//"."
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Apply set the symmetry operations based on the SG number
CALL SYMOPS_SET_SGNUM(nsgnum,i)
IF (i.NE.0) THEN
  nchk=0
  GOTO 1000 ! pipe error through, error messages were posted by SYMOPS_SET_SGNUM
ENDIF
!
800 CONTINUE ! success
nchk=1
GOTO 1000
!
!Error handling
801 CONTINUE ! error, invalid space group name
CALL ATOMSK_MSG(809,(/TRIM(temp)/),(/0.d0/))
nerr = nerr+1
nchk = 0
GOTO 1000
!
!Routine exit.
1000 CONTINUE
!
END SUBROUTINE SYMOPS_SET_SGNAME
!
!
!
!
SUBROUTINE SYMOPS_SET_SGNUM(nsgnumber,nchk)
!
INTEGER,INTENT(IN)::nsgnumber ! space group number (1 - 230)
INTEGER,INTENT(OUT)::nchk ! output success code: 0: failure, 1: success
!
!
INTEGER::nsgnum,nsymnum,ichk,i,nfail
CHARACTER(LEN=128)::temp,msg
!
!Initialization
nchk=0
ichk=0
nsgnum=nsgnumber
nsymnum=0
nfail=0
msg = 'entering SYMOPS_SET_SGNUM'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
IF (nsgnum<1 .OR. nsgnum>sg_nummax) GOTO 801 ! invalid space group number
!Get number of symmetry operations for the space group
CALL SG_NUMGETSYMNUM(nsgnum,nsymnum,ichk)
IF (ichk.NE.1) THEN
  IF (ichk==-1) GOTO 802
  IF (ichk==-2) GOTO 801
ENDIF
IF (nsymnum.LE.0) GOTO 802
!Prepare the symmetry data array
IF (ALLOCATED(symops_trf)) DEALLOCATE(symops_trf)
ALLOCATE(symops_trf(symops_nltrf,nsymnum))
!Initialize symmetry operations
CALL SYMOPS_INIT()
!
!Set the operations from the space group data
DO i=1, nsymnum
  !Get the symmetry operation string
  temp=''
  CALL SG_NUMGETSYMOP(nsgnum,i,temp(1:sg_soplen),ichk)
  IF (ichk.LE.0) GOTO 802
  !Set the operation parameters from the string
  CALL SYMOPS_SET_STR(TRIM(temp),i,ichk)
  IF (ichk.LE.0) THEN ! report interpretation error
                      ! do not stop here, go on
    CALL ATOMSK_MSG(812,(/TRIM(temp)/),(/0.d0/))
    nerr = nerr+1
    nfail = nfail+1
  ENDIF
ENDDO
!
!
800 CONTINUE 
IF (nfail==0) nchk=1 ! success
GOTO 1000
!
!Error handling
801 CONTINUE ! error, invalid space group number
CALL ATOMSK_MSG(810,(/''/),(/DBLE(nsgnum)/))
nerr = nerr+1
nchk = 0
GOTO 1000
802 CONTINUE ! error, failed to access space group data
CALL ATOMSK_MSG(811,(/''/),(/DBLE(nsgnum)/))
nerr = nerr+1
nchk = 0
GOTO 1000
803 CONTINUE ! error, failed to access space group data
CALL ATOMSK_MSG(812,(/TRIM(temp)/),(/0.d0/))
nerr = nerr+1
nchk = 0
GOTO 1000
!
!Routine exit.
1000 CONTINUE
!
END SUBROUTINE SYMOPS_SET_SGNUM
!
!
!
SUBROUTINE SG_APPLY_SYMOPS(sgroup,H,P,S,AUXNAMES,AUX)
! Applies the symmetry operations of the given space group
! The string "sgroup" may contain a group number, or a
! Hermann-Mauguin symbol appropriately formatted for routine "SG_NAMGETNUM"
CHARACTER(LEN=*),INTENT(IN):: sgroup  !Hermann-Mauguin symbol or number of space group
CHARACTER(LEN=128):: msg
CHARACTER(LEN=sg_soplen),DIMENSION(:),ALLOCATABLE:: strsymops
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: AUXNAMES !names of auxiliary properties
INTEGER:: i, k
INTEGER:: nsymnum, nchk
INTEGER:: sgroupnum  !space group number
REAL(dp),DIMENSION(3,3),INTENT(IN):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P, S !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX  !auxiliary properties of atoms/shells
!
IF (ALLOCATED(symops_trf)) DEALLOCATE(symops_trf)
!
! Determine if "sgroup" contains an integer number
READ(sgroup,*,ERR=10,END=10) sgroupnum
! Success => go directly to 200
GOTO 20
10 CONTINUE
! There was an error => sgroup does not contain an integer
! It should contain a Hermann-Mauguin symbol
! Try to read it
CALL SG_NAMGETNUM(TRIM(ADJUSTL(sgroup)),sgroupnum)
CALL ATOMSK_MSG(2132,(/""/),(/DBLE(sgroupnum)/))
!
20 CONTINUE
! Now, sgroupnum is the space group number
WRITE(msg,*) 'space group number = ', sgroupnum
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
IF( sgroupnum <= 0 .OR. sgroupnum > 230 ) THEN
  ! Invalid space group number
  nerr = nerr+1
  CALL ATOMSK_MSG(809,(/TRIM(sgroup)/),(/0.d0/))
  RETURN
ENDIF
!
! Determine the number of symmetry operations
CALL SG_NUMGETSYMNUM(sgroupnum,nsymnum,nchk)
WRITE(msg,*) 'number of symmetry operations = ', nsymnum
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
! Determine symmetry operations
CALL SG_NUMGETSYMOPS(sgroupnum,strsymops,nsymnum,nchk)
!
IF( verbosity==4 ) THEN
  !Print symmetry operations
  WRITE(msg,*) "Symmetry operations (strsymops):"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(strsymops)
    WRITE(msg,'(i3,2X,a32)') i, TRIM( ADJUSTL(strsymops(i)) )
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
ENDIF
!
! Allocate the transformation list "symops_trf"
ALLOCATE(symops_trf(symops_nltrf,nsymnum))
!
! Initialize transformation list (=fill it with zeros)
CALL SYMOPS_INIT()
!
! Fill transformation list
IF( nsymnum>0 ) THEN
  DO i=1,SIZE(strsymops)
    CALL SYMOPS_SET_STR(strsymops(i),i,k)
  ENDDO
ELSE
  CALL ATOMSK_MSG(2758,(/""/),(/0.d0/))
ENDIF
!
IF( verbosity==4 ) THEN
  !Print symmetry operations
  WRITE(msg,*) "Contents of symops_trf:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  IF( ALLOCATED(symops_trf) ) THEN
    DO i=1,SIZE(symops_trf,1)
      k = MIN(12,SIZE(symops_trf,2))
      WRITE(msg,'(i3,2X,20(f6.2,1X))') i, symops_trf(i,1:k)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDDO
  ENDIF
ENDIF
!
! Apply symmetry operations (stored in symops_trf) to the current system
IF( nsymnum>0 .AND. ALLOCATED(symops_trf) ) THEN
  CALL SYMOPS_APPLY(H,P,S,AUXNAMES,AUX,0.5d0,i)
ELSE
  CALL ATOMSK_MSG(2758,(/""/),(/0.d0/))
ENDIF
!
IF (ALLOCATED(symops_trf)) DEALLOCATE(symops_trf)
!
END SUBROUTINE SG_APPLY_SYMOPS
!
!
!
END MODULE symops
