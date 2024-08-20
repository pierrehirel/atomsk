MODULE properties
!
!**********************************************************************************
!*  PROPERTIES                                                                    *
!**********************************************************************************
!* This module reads some system properties from a file.                          *
!* We distinguish between "per-atom properties" which are saved in the            *
!* array AUX as auxiliary properties, and "system-wide properties" which are      *
!* saved in different arrays or strings.                                          *
!**********************************************************************************
!* (C) March 2011 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 24 Aug. 2024                                     *
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
USE crystallography
USE elasticity
USE files_msg
USE exprev
USE messages
USE subroutines
!
!
CONTAINS
!
!
SUBROUTINE READ_PROPERTIES(propfile,H,P,S,ORIENT,C_tensor,AUXNAMES,AUX,SELECT)
!
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: propfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=16):: miller  !Miller indices
CHARACTER(LEN=4096):: msg, msg2
CHARACTER(LEN=4096):: temp, temp2
CHARACTER(LEN=4096),DIMENSION(3):: func_uxyz, func_ui   !functions for displacements
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES, newAUXNAMES !names of auxiliary properties
LOGICAL:: chargeshell !are the charges for shells defined?
LOGICAL:: readprop    !was the property read?
LOGICAL,DIMENSION(3):: areortho
LOGICAL,DIMENSION(:),ALLOCATABLE,INTENT(IN):: SELECT  !mask for atom list
INTEGER:: auxcol, qcol, qscol  !index of columns in AUX
INTEGER:: i, j, k
INTEGER:: Nprop  !number of per-atom properties declared in the file
INTEGER:: status, strlength
INTEGER:: vx, vy, vz !columns in AUX where atom velocities are stored
REAL(dp):: a, b, c, alpha, beta, gamma
REAL(dp):: aniA, aniH !anisotropy factors
REAL(dp):: snumber
REAL(dp):: tempreal, tempreal2, tempreal3
REAL(dp):: u, v, w
REAL(dp),DIMENSION(9):: Voigt     !Elastic constants (Voigt notation)
REAL(dp),DIMENSION(3,3):: rot_matrix  !rotation matrix
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: ORIENT  !crystalographic orientation
REAL(dp),DIMENSION(6,6):: S_tensor  !compliance tensor (6x6 tensor for inversion)
REAL(dp),DIMENSION(9,9),INTENT(INOUT):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(100,3):: tempprop  !temp. array for per-atom properties (max.100 species)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: AUX  !auxiliary properties
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: newAUX             !auxiliary properties (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: P  !atom positions
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(INOUT):: S  !shell positions (if any)
!
!Initialize variables
areortho(:) = .TRUE.
 chargeshell = .FALSE.
auxcol=0
qcol = 0
qscol = 0
status = 0
!
!
WRITE(msg,*) 'Entering READ_PROPERTIES: '
CALL ATOMSK_MSG(999,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
!
CALL ATOMSK_MSG(2073,(/TRIM(propfile)/),(/0.d0/))
!
!Check that file exists
CALL CHECKFILE(propfile,'read')
!
!
!
100 CONTINUE
!Read properties from the file
OPEN(UNIT=35,FILE=propfile,FORM='FORMATTED',STATUS='OLD')
REWIND(35)
DO
  !Initialize variables
  IF(ALLOCATED(newAUXNAMES)) DEALLOCATE(newAUXNAMES)
  IF(ALLOCATED(newAUX)) DEALLOCATE(newAUX)
  Nprop=0
  tempprop(:,:) = 0.d0
  !
  READ(35,'(a128)',END=200,ERR=200) temp
  temp = ADJUSTL(temp)
  CALL ATOMSK_MSG(999,(/TRIM(ADJUSTL(temp))/),(/0.d0/))
  !
  !! SYSTEM-WIDE PROPERTIES !!
  IF(LEN_TRIM(temp).NE.0 .AND. temp(1:1).NE.'#') THEN
    IF(temp(1:9)=='supercell') THEN
      msg = 'supercell vectors'
      READ(35,*,END=800,ERR=800) H(1,1), H(1,2), H(1,3)
      READ(35,*,END=800,ERR=800) H(2,1), H(2,2), H(2,3)
      READ(35,*,END=800,ERR=800) H(3,1), H(3,2), H(3,3)
      CALL ATOMSK_MSG(2074,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
    !
    ELSEIF(temp(1:12)=='conventional') THEN
      msg = 'conventional vectors'
      READ(35,*,END=800,ERR=800) a, b, c
      READ(35,*,END=800,ERR=800) alpha, beta, gamma
      alpha = DEG2RAD(alpha)
      beta = DEG2RAD(beta)
      gamma = DEG2RAD(gamma)
      CALL CONVMAT(a,b,c,alpha,beta,gamma,H)
      CALL ATOMSK_MSG(2074,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
    !
    ELSEIF(temp(1:11)=='orientation') THEN
      msg = 'system orientation'
      READ(35,*,END=800,ERR=800) miller
      WRITE(temp2,*) TRIM(miller)
      CALL ATOMSK_MSG(999,(/TRIM(temp2)/),(/0.d0/))
      CALL INDEX_MILLER(miller,ORIENT(1,:),j)
      IF(j>0) THEN !Failed to read [hkl]: try to read [hkil] indices
        DO i=1,3
          CALL INDEX_MILLER_HCP(miller,ORIENT(i,:),j)
          IF( j>0 ) THEN
            IF( j==2 ) THEN
              !The error was because i is not equal to -h-k
              nerr=nerr+1
              CALL ATOMSK_MSG(815,(/miller/),(/0.d0/))
              GOTO 1000
            ELSE
              !Other error, unable to convert this string into a proper vector
              nerr = nerr+1
              CALL ATOMSK_MSG(817,(/TRIM(miller)/),(/0.d0/))
              GOTO 1000
            ENDIF
          ENDIF
          !Convert [hkil] into [uvw]
          CALL HKIL2UVW(ORIENT(i,1),ORIENT(i,2),0.d0,ORIENT(i,3),u,v,w)
          !Set box vectors in ORIENT
          ORIENT(i,1) = u
          ORIENT(i,2) = v
          ORIENT(i,3) = w
          !read next line
          IF(i<3) THEN
            READ(35,*,END=800,ERR=800) miller
            WRITE(temp2,*) TRIM(miller)
            CALL ATOMSK_MSG(999,(/TRIM(temp2)/),(/0.d0/))
          ENDIF
        ENDDO
        !
      ELSE
        !Succeeded reading [hkl]: read the [hkl] along the other two directions
        READ(35,*,END=800,ERR=800) miller
        WRITE(temp2,*) TRIM(miller)
        CALL ATOMSK_MSG(999,(/TRIM(temp2)/),(/0.d0/))
        CALL INDEX_MILLER(miller,ORIENT(2,:),j)
        IF(j>0) GOTO 800
        READ(35,*,END=800,ERR=800) miller
        WRITE(temp2,*) TRIM(miller)
        CALL ATOMSK_MSG(999,(/TRIM(temp2)/),(/0.d0/))
        CALL INDEX_MILLER(miller,ORIENT(3,:),j)
        IF(j>0) GOTO 800
      ENDIF
      !
      IF(verbosity==4) THEN
        !Debug messages
        DO i=1,3
          WRITE(temp2,'(3(f10.3,1X))') (ORIENT(i,j), j=1,3)
          CALL ATOMSK_MSG(999,(/TRIM(temp2)/),(/0.d0/))
        ENDDO
      ENDIF
      CALL ATOMSK_MSG(2074,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
      !
      !If vectors are not orthonormal, display a warning
      areortho(1) = ORTHOVEC( ORIENT(1,:), ORIENT(2,:) )
      areortho(2) = ORTHOVEC( ORIENT(2,:), ORIENT(3,:)  )
      areortho(3) = ORTHOVEC( ORIENT(3,:), ORIENT(1,:) )
      IF( ANY(.NOT.areortho(:)) ) THEN
        nwarn=nwarn+1
        CALL ATOMSK_MSG(2739,(/""/),(/0.d0/))
      ENDIF
      !
      !If the elastic tensor is defined, rotate it
      IF( C_tensor(1,1).NE.0.d0 ) THEN
        !Define rotation matrix
        DO i=1,3
          rot_matrix(i,:) = ORIENT(i,:)/VECLENGTH(ORIENT(i,:))
        ENDDO
        C_tensor = ROTELAST( C_tensor, rot_matrix )
        CALL ATOMSK_MSG(2099,(/""/),(/0.d0/))
      ENDIF
    !
    ELSEIF(temp(1:7)=='elastic') THEN
      msg = 'elastic tensor'
      C_tensor(:,:) = 0.d0
      !Check if it is Voigt notation
      temp2 = TRIM(ADJUSTL(temp(8:)))
      IF( TRIM(temp2)=="Voigt" .OR. TRIM(temp2)=="voigt" ) THEN
        READ(35,*,END=800,ERR=800) Voigt(1), Voigt(2), Voigt(3)  !C11, C22, C33
        READ(35,*,END=800,ERR=800) Voigt(4), Voigt(5), Voigt(6)  !C23, C31, C12
        READ(35,*,END=800,ERR=800) Voigt(7), Voigt(8), Voigt(9)  !C44, C55, C66
        CALL ELAST2TENSOR(Voigt,C_tensor)
      ELSE
        !Read the 6x6 tensor
        DO i=1,6
          READ(35,*,END=800,ERR=800) (C_tensor(i,j), j=1,6)
        ENDDO
        !Completion of the 9x9 tensor
        C_tensor(1:3,7:9) = C_tensor(1:3,4:6)
        C_tensor(4:6,7:9) = C_tensor(4:6,4:6)
        C_tensor(7:9,1:3) = C_tensor(4:6,1:3)
        C_tensor(7:9,4:6) = C_tensor(4:6,4:6)
        C_tensor(7:9,7:9) = C_tensor(4:6,4:6)
      ENDIF
      IF(verbosity==4) THEN
        temp2 = "Elastic tensor (GPa):"
        CALL ATOMSK_MSG(999,(/TRIM(temp2)/),(/0.d0/))
        DO i=1,9
          WRITE(temp2,'(9(e10.3,2X))') (C_tensor(i,j), j=1,9)
          CALL ATOMSK_MSG(999,(/TRIM(temp2)/),(/0.d0/))
        ENDDO
      ENDIF
      CALL ATOMSK_MSG(2074,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
      !Check tensor for stability criteria
      CALL CTENSOR_STABILITY(C_tensor)
      !Print the anisotropy ratio (A) and anisotropy factor (H)
      aniA = 2.d0*C_tensor(4,4)/(C_tensor(1,1)-C_tensor(1,2))
      aniH = 2.d0*C_tensor(4,4) + C_tensor(1,2) - C_tensor(1,1)
      CALL ATOMSK_MSG(2100,(/""/),(/aniA,aniH/))
    !
    ELSEIF(temp(1:10)=='compliance') THEN
      msg = 'compliance tensor'
      !Read compliances, save them in "S_tensor", and invert it to get C_tensor
      !NOTE: the compliance tensor is a 9x9 tensor, however as such it
      !     is not inversible because all of its lines are not linearly independant
      !     (e.g. S(1,7)=S(1,4), S(1,8)=S(1,5) and so on). In order to invert it,
      !     it is reduced to a 6x6 matrix, which is inverted to get the 6x6 elastic
      !     tensor. Then it is completed to get the 9x9 elastic tensor.
      C_tensor(:,:) = 0.d0
      S_tensor(:,:) = 0.d0
      !Check if it is Voigt notation or 6x6 tensor
      temp2 = TRIM(ADJUSTL(temp(11:)))
      IF( TRIM(temp2)=="Voigt" .OR. TRIM(temp2)=="voigt" ) THEN
        !Voigt notation
        READ(35,*,END=800,ERR=800) Voigt(1), Voigt(2), Voigt(3)
        READ(35,*,END=800,ERR=800) Voigt(4), Voigt(5), Voigt(6)
        READ(35,*,END=800,ERR=800) Voigt(7), Voigt(8), Voigt(9)
        CALL ELAST2TENSOR(Voigt,S_tensor)
      ELSE
        !User gives the 6x6 compliance matrix
        DO i=1,6
          READ(35,*,END=800,ERR=800) (S_tensor(i,j), j=1,6)
        ENDDO
      ENDIF
      IF(verbosity==4) THEN
        temp2 = "Compliance tensor (1/GPa):"
        CALL ATOMSK_MSG(999,(/TRIM(temp2)/),(/0.d0/))
        DO i=1,6
          WRITE(temp2,'(6(e10.3,2X))') (S_tensor(i,j), j=1,6)
          CALL ATOMSK_MSG(999,(/TRIM(temp2)/),(/0.d0/))
        ENDDO
      ENDIF
      !Invert the compliance tensor to get the elastic tensor
      CALL INVMAT(S_tensor,C_tensor(1:6,1:6),i)
      !If i is non-zero then the inversion failed
      IF(i.NE.0) THEN
        nerr=nerr+1
        CALL ATOMSK_MSG(2815,(/"C_tensor"/),(/0.d0/))
        GOTO 800
      ENDIF
      !Otherwise (i=0) the C_tensor(:,:) now contains the 6x6 elastic tensor
      !Completion of the 9x9 elastic tensor
      C_tensor(1:3,7:9) = C_tensor(1:3,4:6)
      C_tensor(4:6,7:9) = C_tensor(4:6,4:6)
      C_tensor(7:9,1:3) = C_tensor(4:6,1:3)
      C_tensor(7:9,4:6) = C_tensor(4:6,4:6)
      C_tensor(7:9,7:9) = C_tensor(4:6,4:6)
      IF(verbosity==4) THEN
        temp2 = "Elastic tensor (GPa):"
        CALL ATOMSK_MSG(999,(/TRIM(temp2)/),(/0.d0/))
        DO i=1,9
          WRITE(temp2,'(9(e10.3,2X))') (C_tensor(i,j), j=1,9)
          CALL ATOMSK_MSG(999,(/TRIM(temp2)/),(/0.d0/))
        ENDDO
      ENDIF
      CALL ATOMSK_MSG(2074,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
      !Check tensor for stability criteria
      CALL CTENSOR_STABILITY(C_tensor)
      !Print the anisotropy ratio (A) and anisotropy factor (H)
      aniA = 2.d0*C_tensor(4,4)/(C_tensor(1,1)-C_tensor(1,2))
      aniH = 2.d0*C_tensor(4,4) + C_tensor(1,2) - C_tensor(1,1)
      CALL ATOMSK_MSG(2100,(/""/),(/aniA,aniH/))
    !
    ! -- add other system-wide properties here --
    !
    !
    !! PER-ATOM PROPERTIES !!
    ELSEIF(temp(1:12)=='displacement' .OR. temp(1:4)=='disp') THEN
      IF( INDEX(temp,"function")>0 ) THEN
        msg = 'atom displacements functions'
        func_uxyz(:) = ""
        func_ui(:) = ""
        status = 0
        DO j=1,3
          READ(35,'(a128)',END=171,ERR=171) msg2
          msg2 = ADJUSTL(msg2)
          IF( msg2(1:2)=="ux" .OR. msg(1:2)=="Ux" .OR. msg2(1:2)=="UX" ) THEN
            k = SCAN(msg2,"=")
            IF(k==0) k=2
            func_uxyz(1) = ADJUSTL(msg2(k+1:))
          ELSEIF( msg2(1:2)=="uy" .OR. msg(1:2)=="Uy" .OR. msg2(1:2)=="UY" ) THEN
            k = SCAN(msg2,"=")
            IF(k==0) k=2
            func_uxyz(2) = ADJUSTL(msg2(k+1:))
          ELSEIF( msg2(1:2)=="uz" .OR. msg(1:2)=="Uz" .OR. msg2(1:2)=="UZ" ) THEN
            k = SCAN(msg2,"=")
            IF(k==0) k=2
            func_uxyz(3) = ADJUSTL(msg2(k+1:))
          ENDIF
        ENDDO
        171 CONTINUE
        !For each function, replace Hx, Hy, Hz by the appropriate box length
        DO j=1,3
          tempreal = MAX(0.d0,MAXVAL(H(:,1))) + DABS(MIN(0.d0,MINVAL(H(:,1))))
          CALL STR_EXP2VAL(func_uxyz(j),"Hx",tempreal,strlength)
          tempreal = MAX(0.d0,MAXVAL(H(:,2))) + DABS(MIN(0.d0,MINVAL(H(:,2))))
          CALL STR_EXP2VAL(func_uxyz(j),"Hy",tempreal,strlength)
          tempreal = MAX(0.d0,MAXVAL(H(:,3))) + DABS(MIN(0.d0,MINVAL(H(:,3))))
          CALL STR_EXP2VAL(func_uxyz(j),"Hz",tempreal,strlength)
        ENDDO
        CALL ATOMSK_MSG(999,(/"ux = "//TRIM(func_uxyz(1))/),(/0.d0/))
        CALL ATOMSK_MSG(999,(/"uy = "//TRIM(func_uxyz(2))/),(/0.d0/))
        CALL ATOMSK_MSG(999,(/"uz = "//TRIM(func_uxyz(3))/),(/0.d0/))
        !Apply displacements immediately to all atoms
        CALL ATOMSK_MSG(2146,(/""/),(/0.d0/))
        !$OMP PARALLEL DO DEFAULT(SHARED) &
        !$OMP& PRIVATE(i,j,func_ui,strlength,status,tempreal,tempreal2,tempreal3)
        DO i=1,SIZE(P,1)
          IF( IS_SELECTED(SELECT,i) ) THEN
            !In each function string, replace the  variables x, y, z
            !by their actual values for atom i
            status = 0
            strlength = 0
            tempreal = 0.d0
            tempreal2 = 0.d0
            tempreal3 = 0.d0
            DO j=1,3
              func_ui(j) = func_uxyz(j)
              CALL STR_EXP2VAL(func_ui(j),"box",H(j,j),strlength)
              CALL STR_EXP2VAL(func_ui(j),"BOX",H(j,j),strlength)
              CALL STR_EXP2VAL(func_ui(j),"x",P(i,1),strlength)
              CALL STR_EXP2VAL(func_ui(j),"X",P(i,1),strlength)
              CALL STR_EXP2VAL(func_ui(j),"y",P(i,2),strlength)
              CALL STR_EXP2VAL(func_ui(j),"Y",P(i,2),strlength)
              CALL STR_EXP2VAL(func_ui(j),"z",P(i,3),strlength)
              CALL STR_EXP2VAL(func_ui(j),"Z",P(i,3),strlength)
              !Evaluate the function with those values
              IF( LEN_TRIM(func_ui(j))>0 ) THEN
                strlength=0
                status=0
                tempreal3 = EXPREVAL(func_ui(j),strlength,status)
              ELSE
                tempreal3 = 0.d0
              ENDIF
              IF(status==10) THEN
                !There was a division by zero
                CALL ATOMSK_MSG(816,(/""/),(/0.d0/))
                !Just consider zero displacement
                tempreal3 = 0.d0
              ELSEIF(status>0) THEN
                !Unable to convert that string into a number
                CALL ATOMSK_MSG(2813,(/func_ui(j)/),(/0.d0/))
                tempreal3 = 0.d0
                nerr = nerr+1
                IF( nerr>50 ) THEN
                  !Too many errors, abort
                  EXIT
                ENDIF
              ENDIF
              IF(j==1) tempreal = tempreal3
              IF(j==2) tempreal2 = tempreal3
            ENDDO
            !Check that displacement vector is not too large
            IF( VECLENGTH( (/tempreal,tempreal2,tempreal3/) ) > 100.d0 ) THEN
              nwarn = nwarn+1
              CALL ATOMSK_MSG(2727,(/""/),(/DBLE(i)/))
            ENDIF
            !Apply displacements
            P(i,1) = P(i,1) + tempreal
            P(i,2) = P(i,2) + tempreal2
            P(i,3) = P(i,3) + tempreal3
          ENDIF  !end if SELECT
        ENDDO
        !$OMP END PARALLEL DO
        IF(nerr>0) GOTO 1000
        !
      ELSE
        msg = 'atom displacements'
        !Displacement vectors of atoms (and shells) follow
        CALL ATOMSK_MSG(2146,(/""/),(/0.d0/))
        !Read the displacement for each atom and apply it immediately
        !Note: the values may not be given for all atoms
        DO
          j=0
          READ(35,'(a128)',END=173,ERR=173) msg2
          IF( LEN_TRIM(msg2)>0 ) THEN
            !Read displacement of atom (or ionic core)
            READ(msg2,*,END=173,ERR=173) i, a, b, c
            !Read displacement of ionic shell (optional)
            READ(msg2,*,END=172,ERR=172) i, a, b, c, tempreal, tempreal2, tempreal3
            j=1
            172 CONTINUE
            IF( i>0 .AND. i<=SIZE(P,1) ) THEN
              IF( IS_SELECTED(SELECT,i) ) THEN
                P(i,1) = P(i,1) + a
                P(i,2) = P(i,2) + b
                P(i,3) = P(i,3) + c
                IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
                  IF( j>0 ) THEN
                    !Displacements were explicitely given for the shell => use that
                    S(i,1) = S(i,1) + tempreal
                    S(i,2) = S(i,2) + tempreal2
                    S(i,3) = S(i,3) + tempreal3
                  ELSE
                    !User did not specify any displacement for the shell
                    !=> use same displacement vector as for the core
                    S(i,1) = S(i,1) + a
                    S(i,2) = S(i,2) + b
                    S(i,3) = S(i,3) + c
                  ENDIF  ! end if j>0
                ENDIF  !end if S
              ENDIF  !end if SELECT
            ELSE
              !given index is out-of-bounds
              nwarn=nwarn+1
              CALL ATOMSK_MSG(2742,(/""/),(/DBLE(i)/))
            ENDIF
          ELSE
            EXIT
          ENDIF
        ENDDO
      ENDIF
      173 CONTINUE
      CALL ATOMSK_MSG(2074,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
      !
    ELSEIF(temp(1:6)=='charge') THEN
      msg = 'atom charges'
      !Read charges and store it in "tempprop"
      DO
        READ(35,'(a128)',END=136,ERR=135) temp2
        READ(temp2,*,END=136,ERR=135) species, tempreal
        IF(LEN_TRIM(species)==0) GOTO 800
        !
        WRITE(msg2,*) TRIM(species), tempreal
        CALL ATOMSK_MSG(999,(/TRIM(ADJUSTL(msg2))/),(/0.d0/))
        !
        CALL ATOMNUMBER(species,tempreal2)
        IF(tempreal2>1.d-9) THEN
          Nprop=Nprop+1
          tempprop(Nprop,1) = tempreal2
          tempprop(Nprop,2) = tempreal
        ELSE
          GOTO 135
        ENDIF
        !
        !Look if the charge of shell is also defined
        tempreal2=0.d0
        READ(temp2,*,END=131,ERR=131) species, tempreal, tempreal2
        IF( tempreal2.NE.0.d0 ) THEN
          chargeshell = .TRUE.
          tempprop(Nprop,3) = tempreal2
        ENDIF
        !
        131 CONTINUE
      ENDDO
      135 CONTINUE
      !Last line did not contain charge => go back two lines
      BACKSPACE(35)
      !Last line was an error => go back one line
      136 CONTINUE
      BACKSPACE(35)
      !The charges must be saved as an auxiliary property in the array "AUX"
      IF(ALLOCATED(AUX)) THEN
        !AUX is already defined,
        !check if it contains q and qs already
        qcol=0
        DO i=1,SIZE(AUXNAMES)
          IF( AUXNAMES(i)=="q" ) THEN
            qcol=i
            EXIT
          ELSEIF( AUXNAMES(i)=="qs" ) THEN
            qscol=i
            EXIT
          ENDIF
        ENDDO
        !If if doesn't contain q, expand the arrays
        !Also expand it if we need to store the charges of shells
        IF(qcol==0) THEN
          qcol = SIZE(AUX,2)+1
          IF(chargeshell) THEN
            qscol = SIZE(AUX,2)+2
            ALLOCATE( newAUX( SIZE(AUX,1), qscol ) )
          ELSE
            ALLOCATE( newAUX( SIZE(AUX,1), qcol ) )
          ENDIF
          DO i=1,SIZE(AUX,1)
            DO j=1,SIZE(AUX,2)
              newAUX(i,j) = AUX(i,j)
            ENDDO
          ENDDO
          newAUX(:,qcol) = 0.d0
          IF(chargeshell) THEN
            newAUX(:,qscol) = 0.d0
          ENDIF
          DEALLOCATE(AUX)
          ALLOCATE( AUX( SIZE(newAUX,1), SIZE(newAUX,2) ) )
          AUX(:,:) = newAUX(:,:)
          DEALLOCATE(newAUX)
          !Same with array AUXNAMES
          IF(chargeshell) THEN
            ALLOCATE(newAUXNAMES( qscol ))
          ELSE
            ALLOCATE(newAUXNAMES( qcol ))
          ENDIF
          DO i=1,SIZE(AUXNAMES)
            newAUXNAMES(i) = AUXNAMES(i)
          ENDDO
          DEALLOCATE(AUXNAMES)
          ALLOCATE(AUXNAMES( SIZE(newAUXNAMES) ))
          AUXNAMES(:) = newAUXNAMES(:)
          DEALLOCATE(newAUXNAMES)
          !
        ELSEIF( qscol==0 .AND. chargeshell ) THEN
          qscol = SIZE(AUX,2)+1
          ALLOCATE( newAUX( SIZE(AUX,1), qscol ) )
          DO i=1,SIZE(AUX,1)
            DO j=1,SIZE(AUX,2)
              newAUX(i,j) = AUX(i,j)
            ENDDO
          ENDDO
          newAUX(:,qscol) = 0.d0
          DEALLOCATE(AUX)
          ALLOCATE( AUX( SIZE(newAUX,1), SIZE(newAUX,2) ) )
          AUX(:,:) = newAUX(:,:)
          DEALLOCATE(newAUX)
          !Same with array AUXNAMES
            ALLOCATE(newAUXNAMES( qscol ))
          DO i=1,SIZE(AUXNAMES(:))
            newAUXNAMES(i) = AUXNAMES(i)
          ENDDO
          DEALLOCATE(AUXNAMES)
          ALLOCATE(AUXNAMES( SIZE(newAUXNAMES(:)) ))
          AUXNAMES(:) = newAUXNAMES(:)
          DEALLOCATE(newAUXNAMES)
        ENDIF
        !
      ELSE
        !No auxiliary property exist => allocate AUX and store charges inside
        IF(chargeshell) THEN
          ALLOCATE( AUX(SIZE(P,1),2) )
          ALLOCATE( AUXNAMES(2) )
          qcol = 1
          qscol = 2
        ELSE
          ALLOCATE( AUX(SIZE(P,1),1) )
          ALLOCATE( AUXNAMES(1) )
          qcol = 1
        ENDIF
        AUX(:,:) = 0.d0
      ENDIF
      !Set names of auxiliary properties
      AUXNAMES(qcol) = 'q'
      IF(chargeshell) THEN
        AUXNAMES(qscol) = 'qs'
      ENDIF
      !Save charges in AUX
      DO i=1,SIZE(AUX,1)
        DO j=1,SIZE(tempprop,1)
          IF( DABS( tempprop(j,1)-P(i,4) )<1.d-9 ) THEN
            AUX(i,qcol) = tempprop(j,2)
            IF(chargeshell) THEN
              AUX(i,qscol) = tempprop(j,3)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      CALL ATOMSK_MSG(2074,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
    !
    ELSEIF(temp(1:4)=='type') THEN
      msg = 'atom types'
      !Read pairs of (species,type) and store it in "tempprop"
      DO
        READ(35,*,END=141,ERR=140) species, tempreal
        IF(LEN_TRIM(species)==0) EXIT
        !
        WRITE(temp2,*) TRIM(species), tempreal
        CALL ATOMSK_MSG(999,(/TRIM(ADJUSTL(temp2))/),(/0.d0/))
        !
        CALL ATOMNUMBER(species,tempreal2)
        IF(tempreal2>1.d-9) THEN
          Nprop=Nprop+1
          tempprop(Nprop,1) = tempreal2
          tempprop(Nprop,2) = tempreal
        ELSE
          GOTO 140
        ENDIF
      ENDDO
      140 CONTINUE
      !Last line did not contain type => go back two lines
      BACKSPACE(35)
      141 CONTINUE
      BACKSPACE(35)
      !The types must be saved as an auxiliary property in the array "AUX"
      IF(ALLOCATED(AUX)) THEN
        !AUX is already defined,
        !check if it contains type already
        auxcol=0
        DO i=1,SIZE(AUXNAMES(:))
          IF( AUXNAMES(i)=="type" ) THEN
            auxcol=i
            EXIT
          ENDIF
        ENDDO
        !If if doesn't contain type, expand the arrays
        IF(auxcol==0) THEN
          auxcol = SIZE(AUX,2)+1
          ALLOCATE( newAUX( SIZE(AUX,1), auxcol ) )
          DO i=1,SIZE(AUX,1)
            DO j=1,SIZE(AUX,2)
              newAUX(i,j) = AUX(i,j)
            ENDDO
          ENDDO
          newAUX(:,auxcol) = 0.d0
          DEALLOCATE(AUX)
          ALLOCATE( AUX( SIZE(newAUX,1), auxcol ) )
          AUX(:,:) = newAUX(:,:)
          DEALLOCATE(newAUX)
          !Same with array AUXNAMES
          ALLOCATE(newAUXNAMES( auxcol ))
          DO i=1,SIZE(AUXNAMES(:))
            newAUXNAMES(i) = AUXNAMES(i)
          ENDDO
          DEALLOCATE(AUXNAMES)
          ALLOCATE(AUXNAMES(auxcol))
          AUXNAMES(:) = newAUXNAMES(:)
          DEALLOCATE(newAUXNAMES)
        ENDIF
        !
      ELSE
        !No auxiliary property exist => allocate AUX and store types inside
        ALLOCATE( AUX(SIZE(P,1),1) )
        AUX(:,:) = 0.d0
        ALLOCATE( AUXNAMES(1) )
        auxcol = 1
      ENDIF
      !Set 'type' as name of auxiliary property
      AUXNAMES(auxcol) = 'type'
      !Save atom types in AUX
      DO i=1,SIZE(P,1)
        DO j=1,SIZE(tempprop,1)
          IF( DABS( tempprop(j,1)-P(i,4) )<1.d-9 ) THEN
            AUX(i,auxcol) = tempprop(j,2)
          ENDIF
        ENDDO
      ENDDO
      CALL ATOMSK_MSG(2074,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
    !
    ELSEIF(temp(1:8)=='velocity') THEN
      msg = 'atom velocities'
      temp2 = TRIM(ADJUSTL(temp(9:)))
      msg2=""
      READ(temp2,*,ERR=150,END=150) msg2
      150 CONTINUE
      !Read velocities and save them in the array "AUX"
      IF(ALLOCATED(AUX)) THEN
        !AUX is already defined,
        !check if it contains vx, vy, vz already
        vx=0
        vy=0
        vz=0
        DO i=1,SIZE(AUXNAMES(:))
          IF( TRIM(ADJUSTL(AUXNAMES(i))) == "vx" ) THEN
            vx = i
          ELSEIF( TRIM(ADJUSTL(AUXNAMES(i))) == "vy" ) THEN
            vy = i
          ELSEIF( TRIM(ADJUSTL(AUXNAMES(i))) == "vz" ) THEN
            vz = i
          ENDIF
        ENDDO
        !If if doesn't contain that property, expand the arrays
        IF( vx==0 .OR. vy==0 .OR. vz==0 ) THEN
          vx = SIZE(AUX,2)+1
          vy = SIZE(AUX,2)+2
          vz = SIZE(AUX,2)+3
          ALLOCATE( newAUX( SIZE(AUX,1), vz ) )
          DO i=1,SIZE(AUX,1)
            DO j=1,SIZE(AUX,2)
              newAUX(i,j) = AUX(i,j)
            ENDDO
          ENDDO
          newAUX(:,auxcol) = 0.d0
          DEALLOCATE(AUX)
          ALLOCATE( AUX( SIZE(newAUX,1), SIZE(newAUX,2) ) )
          AUX(:,:) = newAUX(:,:)
          DEALLOCATE(newAUX)
          !Same with array AUXNAMES
          ALLOCATE(newAUXNAMES( vz ))
          DO i=1,SIZE(AUXNAMES(:))
            newAUXNAMES(i) = AUXNAMES(i)
          ENDDO
          DEALLOCATE(AUXNAMES)
          ALLOCATE(AUXNAMES(auxcol))
          AUXNAMES(:) = newAUXNAMES(:)
          DEALLOCATE(newAUXNAMES)
        ENDIF
        !
      ELSE
        !No auxiliary property exist at all => allocate AUX and store velocities inside
        ALLOCATE( AUX(SIZE(P,1),3) )
        AUX(:,:) = 0.d0
        ALLOCATE( AUXNAMES(3) )
        vx = 1
        vy = 2
        vz = 3
      ENDIF
      !Set the name of that auxiliary property
      AUXNAMES(vx) = "vx"
      AUXNAMES(vy) = "vy"
      AUXNAMES(vz) = "vz"
      !Read the velocity for each atom and save it
      !Note: the values may not be given for all atoms
      DO
        readprop=.FALSE.
        READ(35,'(a128)',END=153,ERR=153) msg2
        READ(msg2,*,END=152,ERR=152) i, tempreal, tempreal2, tempreal3
        IF( i>0 .AND. i<=SIZE(AUX,1) ) THEN
          AUX(i,vx) = tempreal
          AUX(i,vy) = tempreal2
          AUX(i,vz) = tempreal3
          readprop=.TRUE.
        ELSE
          !given index is out-of-bounds
          nwarn=nwarn+1
          CALL ATOMSK_MSG(2742,(/""/),(/DBLE(i)/))
        ENDIF
        !
        152 CONTINUE
        IF(.NOT.readprop) THEN
          CALL ATOMSK_MSG(2749,(/AUXNAMES(auxcol)/),(/DBLE(i)/))
          nwarn=nwarn+1
        ENDIF
      ENDDO
      153 CONTINUE
      CALL ATOMSK_MSG(2074,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
      !
    ELSEIF(temp(1:9)=='auxiliary') THEN
      !A list of values for atoms is defined
      !Read the name of the auxiliary property
      msg = TRIM(ADJUSTL(temp(10:)))
      IF( LEN_TRIM(msg)==0 ) msg = "undefined"
      !This auxiliary property must be saved in the array "AUX"
      IF(ALLOCATED(AUX)) THEN
        !AUX is already defined,
        !check if it contains that property already
        auxcol=0
        DO i=1,SIZE(AUXNAMES(:))
          IF( TRIM(ADJUSTL(AUXNAMES(i))) == TRIM(ADJUSTL(msg)) ) THEN
            auxcol=i
            EXIT
          ENDIF
        ENDDO
        !If if doesn't contain that property, expand the arrays
        IF(auxcol==0) THEN
          auxcol = SIZE(AUX,2)+1
          ALLOCATE( newAUX( SIZE(AUX,1), auxcol ) )
          DO i=1,SIZE(AUX,1)
            DO j=1,SIZE(AUX,2)
              newAUX(i,j) = AUX(i,j)
            ENDDO
          ENDDO
          newAUX(:,auxcol) = 0.d0
          DEALLOCATE(AUX)
          ALLOCATE( AUX( SIZE(newAUX,1), auxcol ) )
          AUX(:,:) = newAUX(:,:)
          DEALLOCATE(newAUX)
          !Same with array AUXNAMES
          ALLOCATE(newAUXNAMES( auxcol ))
          DO i=1,SIZE(AUXNAMES(:))
            newAUXNAMES(i) = AUXNAMES(i)
          ENDDO
          DEALLOCATE(AUXNAMES)
          ALLOCATE(AUXNAMES(auxcol))
          AUXNAMES(:) = newAUXNAMES(:)
          DEALLOCATE(newAUXNAMES)
        ENDIF
        !
      ELSE
        !No auxiliary property exist => allocate AUX and store types inside
        ALLOCATE( AUX(SIZE(P,1),1) )
        AUX(:,:) = 0.d0
        ALLOCATE( AUXNAMES(1) )
        auxcol = 1
      ENDIF
      !Set the name of that auxiliary property
      AUXNAMES(auxcol) = TRIM(msg)
      !Read the property for each atom and save it
      !Note: the values may not be given for all atoms
      DO
        readprop=.FALSE.
        READ(35,'(a128)',END=163,ERR=163) msg2
        msg2 = ADJUSTL(msg2)
        IF( LEN_TRIM(msg2)>0 ) THEN
          READ(msg2,*,END=162,ERR=162) i, tempreal
          ! Succeeded reading an integer number i
          IF( i>0 .AND. i<=SIZE(AUX,1) ) THEN
            AUX(i,auxcol) = tempreal
            readprop=.TRUE.
          ELSE
            !given index is out-of-bounds
            nwarn=nwarn+1
            CALL ATOMSK_MSG(2742,(/""/),(/DBLE(i)/))
          ENDIF
          !
          162 CONTINUE
          ! Failed reading an integer number i
          ! => maybe it is an atom species
          IF(.NOT.readprop) THEN
            species = msg2(1:2)
            CALL ATOMNUMBER(species,snumber)
            IF( snumber>1.d-6 ) THEN
              !We have an atomic number
              !Read the value of the property
              READ(msg2,*,END=162,ERR=162) species, tempreal
              !Give this value of the property to all atoms of this species
              DO i=1,SIZE(P,1)
                IF( NINT(P(i,4))==NINT(snumber) ) THEN
                  AUX(i,auxcol) = tempreal
                ENDIF
              ENDDO
            ELSE
              !Cannot make sense out of this line
              CALL ATOMSK_MSG(2749,(/AUXNAMES(auxcol)/),(/DBLE(i)/))
              nwarn=nwarn+1
            ENDIF
          ENDIF
        ELSE
          EXIT
        ENDIF
      ENDDO
      163 CONTINUE
      CALL ATOMSK_MSG(2074,(/AUXNAMES(auxcol)/),(/0.d0/))
    !
    ! -- add other per-atom properties here --
    !
    ELSE
      nwarn=nwarn+1
      CALL ATOMSK_MSG(2732,(/TRIM(ADJUSTL(temp))/),(/0.d0/))
      GOTO 190
    ENDIF
    !
  ENDIF
  !
  190 CONTINUE
ENDDO
!
!
!
200 CONTINUE
CLOSE(35)
CALL ATOMSK_MSG(2075,(/''/),(/0.d0/))
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(2802,(/TRIM(ADJUSTL(msg))/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE READ_PROPERTIES
!
!
!
SUBROUTINE CTENSOR_STABILITY(C_tensor)
!
CHARACTER(LEN=32):: temp, msg
INTEGER:: i, j
REAL(dp),DIMENSION(9,9),INTENT(IN):: C_tensor  !elastic tensor
!
IF( C_tensor(1,1) < DABS(C_tensor(1,2)) ) THEN
  msg = "C11 > |C12|"
  CALL ATOMSK_MSG(2762,(/msg/),(/0.d0/))
ENDIF
DO i=1,6
  IF( C_tensor(i,i)<=0.d0 ) THEN
    WRITE(temp,'(2i1)') i, i
    msg = "C"//TRIM(ADJUSTL(temp))//" > 0"
    CALL ATOMSK_MSG(2762,(/msg/),(/0.d0/))
  ENDIF
  DO j=i,6
    IF( j>i ) THEN
      IF( C_tensor(i,j)**2 > C_tensor(i,i)*C_tensor(j,j) ) THEN
        WRITE(temp,'(2i1)') i, i
        msg = "C"//TRIM(ADJUSTL(temp))
        WRITE(temp,'(2i1)') j, j
        msg = TRIM(ADJUSTL(msg))//" * C"//TRIM(ADJUSTL(temp))
        WRITE(temp,'(2i1)') i, j
        msg = "(C"//TRIM(ADJUSTL(temp))//")² < "//TRIM(ADJUSTL(msg))
        CALL ATOMSK_MSG(2762,(/msg/),(/0.d0/))
      ENDIF
    ENDIF
  ENDDO
ENDDO
!
END SUBROUTINE CTENSOR_STABILITY
!
!
!
END MODULE properties
