MODULE mode_average
!
!**********************************************************************************
!*  MODE_AVERAGE                                                                  *
!**********************************************************************************
!* This module reads many files containing atom positions, and produces           *
!* an output file containing an average of atom positions.                        *
!**********************************************************************************
!* (C) Dec. 2013 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 16 April 2024                                    *
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
USE comv
USE constants
USE messages
USE files
USE subroutines
USE readin
USE options
USE writeout
!
!
CONTAINS
!
SUBROUTINE AVERAGE_XYZ(listfile,outputfile,outfileformats,options_array)
!
!Declare variables
IMPLICIT NONE
!Input
CHARACTER(LEN=*),INTENT(IN):: listfile  !file containing the names of files to analyze
!
CHARACTER(LEN=128):: msg, temp
CHARACTER(LEN=4096):: inputfile       !name of a file to analyze
CHARACTER(LEN=*):: outputfile      !name of the final output file
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to output
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES      !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMEStemp  !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment, commenttemp
LOGICAL:: auxok
LOGICAL:: fileexists !does the file exist?
LOGICAL:: first      !is this the first file?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: i, j
INTEGER:: Nfiles   !number of files
REAL(dp),DIMENSION(3,3):: Huc   !Box vectors of unit cell (unknown, set to 0 here)
REAL(dp),DIMENSION(3,3):: H       !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: Htemp   !Base vectors of the supercell (temporary)
REAL(dp),DIMENSION(3,3):: ORIENT  !crystal orientation
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX          !auxiliary properties of atoms
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUXtemp      !auxiliary properties of atoms (temporary)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, S         !positions of atoms, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Ptemp, Stemp !positions of atoms, shells (temporary)
!
msg = 'ENTERING AVERAGE_XYZ...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
CALL ATOMSK_MSG(4064,(/listfile/),(/0.d0/))
!
!Initialize variables
H(:,:) = 0.d0
Huc(:,:) = 0.d0
 C_tensor(:,:) = 0.d0
ORIENT(:,:) = 0.d0
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
!
!
!
100 CONTINUE
OPEN(UNIT=50,FILE=listfile,STATUS='OLD',FORM='FORMATTED')
REWIND(50)
!
!
!
200 CONTINUE
!Read each file and store the sum of positions for each atom
!Then the final total sum will be divided by the number of files
Nfiles=0
first = .TRUE.
DO
  !
  !Read the name of the file
  READ(50,'(a128)',END=300,ERR=300) inputfile
  inputfile = ADJUSTL(inputfile)
  !
  IF( temp(1:1).NE.'#' ) THEN
    !Check if file actually exists
    INQUIRE(FILE=inputfile,EXIST=fileexists)
    !
    IF(fileexists) THEN
      !Read atom positions
      CALL READ_AFF(inputfile,Htemp,Ptemp,Stemp,commenttemp,AUXNAMEStemp,AUXtemp)
      IF(nerr>0 .OR. .NOT.ALLOCATED(P)) GOTO 1000
      !
      !Convert atom positions to reduced coordinates
      CALL CART2FRAC(Ptemp,Htemp)
      IF( ALLOCATED(Stemp) .AND. SIZE(Stemp,1)==SIZE(P,1) ) THEN
        CALL CART2FRAC(Stemp,Htemp)
      ENDIF
      !
      IF( first ) THEN
        !First file => allocate arrays
        Nfiles = Nfiles+1
        H(:,:) = Htemp(:,:)
        ALLOCATE(P(SIZE(Ptemp,1),4))
        P(:,:) = Ptemp(:,:)
        IF( ALLOCATED(Stemp) .AND. SIZE(Stemp,1)==SIZE(P,1) ) THEN
          ALLOCATE( S(SIZE(P,1),4) )
          S(:,:) = Stemp(:,:)
        ENDIF
        IF( ALLOCATED(AUXtemp) .AND. SIZE(AUXtemp,1)==SIZE(P,1) ) THEN
          ALLOCATE( AUX(SIZE(AUXtemp,1),SIZE(AUXtemp,2)) )
          AUX(:,:) = AUXtemp(:,:)
          ALLOCATE( AUXNAMES(SIZE(AUX,2)) )
          AUXNAMES(:) = AUXNAMEStemp(:)
        ENDIF
        IF( ALLOCATED(commenttemp) .AND. SIZE(commenttemp)>0 ) THEN
          ALLOCATE( comment(SIZE(commenttemp)) )
          comment(:) = commenttemp(:)
        ELSE
          ALLOCATE( comment(1) )
          comment(1) = ""
        ENDIF
        first = .FALSE.
        !
      ELSE
        !Check that this file has the same number of atoms as the first file
        IF( SIZE(Ptemp,1) == SIZE(P,1) ) THEN
          !Add these positions to the previous ones
          Nfiles = Nfiles+1
          H(:,:) = H(:,:) + Htemp(:,:)
          P(:,1:3) = P(:,1:3) + Ptemp(:,1:3)
          !
          IF( ALLOCATED(S) .AND. SIZE(S,1)==SIZE(P,1) ) THEN
            IF( ALLOCATED(Stemp) .AND. SIZE(Stemp,1)==SIZE(S,1) ) THEN
              !Current files has shells
              S(:,1:3) = S(:,1:3) + Stemp(:,1:3)
            ELSE
              !Current file does not have shells => ignore shells for the average
              IF(ALLOCATED(S)) DEALLOCATE(S)
            ENDIF
            !
          ELSE
            IF(ALLOCATED(S)) DEALLOCATE(S)
          ENDIF
          !
          IF( ALLOCATED(AUX) .AND. SIZE(AUX,1)==SIZE(P,1) .AND. SIZE(AUX,2)==SIZE(AUXNAMES) ) THEN
            IF( ALLOCATED(AUXtemp) .AND. SIZE(AUXtemp,1)==SIZE(AUX,1) .AND. SIZE(AUXNAMEStemp)==SIZE(AUXNAMES) ) THEN
              !Current files has auxiliary properties
              !Make sure they are the same as the previous file(s)
              auxok = .TRUE.
              DO i=1,SIZE(AUXNAMEStemp)
                IF( AUXNAMEStemp(i) .NE. AUXNAMES(i) ) THEN
                  !This property does not match => remove all auxiliary properties
                  auxok = .FALSE.
                ENDIF
              ENDDO
              !
              IF( auxok ) THEN
                !General case: auxiliary properties are also averaged blindly
                !However we know that for some of them it doesn't make sense:
                !atom type
                DO j=1,SIZE(AUX,2)
                  IF( AUXNAMES(j).NE."type" ) THEN
                    AUX(:,j) = AUX(:,j) + AUXtemp(:,j)
                  ENDIF
                ENDDO
                !
              ELSE
                IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
                IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
              ENDIF
              !
            ELSE
              !Current files does not have auxiliary properties => ignore aux.prop. for average
              IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
              IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
            ENDIF
            !
          ENDIF
          !
        ELSE
          !This system has a different number of atoms, warn and skip to next file
          nwarn=nwarn+1
          CALL ATOMSK_MSG(4711,(/TRIM(inputfile)/),(/0.d0/))
        ENDIF
        !
      ENDIF
      !
      !Free arrays
      IF(ALLOCATED(Ptemp)) DEALLOCATE(Ptemp)
      IF(ALLOCATED(Stemp)) DEALLOCATE(Stemp)
      IF(ALLOCATED(AUXtemp)) DEALLOCATE(AUXtemp)
      IF(ALLOCATED(AUXNAMEStemp)) DEALLOCATE(AUXNAMEStemp)
      !
    ELSE
      !File does not exist => display a warning and skip to next file
      nwarn=nwarn+1
      CALL ATOMSK_MSG(4700,(/TRIM(inputfile)/),(/0.d0/))
      !
    ENDIF
    !
  ENDIF !end if temp.NE."#"
        
  !
ENDDO  !loop on m files
!
!
!
300 CONTINUE
CLOSE(50)
!
!Divide positions by number of files
IF(Nfiles>1) THEN
  H(:,:) = H(:,:) / DBLE(Nfiles)
  P(:,1:3) = P(:,1:3) / DBLE(Nfiles)
  IF( ALLOCATED(S) .AND. SIZE(S,1)>0 ) THEN
    S(:,1:3) = S(:,1:3) / DBLE(Nfiles)
  ENDIF
  IF( ALLOCATED(AUX) .AND. SIZE(AUX,1)>0 ) THEN
    DO j=1,SIZE(AUX,2)
      IF( AUXNAMES(j).NE."type" ) THEN
        AUX(:,:) = AUX(:,:) / DBLE(Nfiles)
      ENDIF
    ENDDO
  ENDIF
  CALL ATOMSK_MSG(4065,(/""/),(/DBLE(Nfiles)/))
  !
ELSE
  !no file was analyzed => exit
  CALL ATOMSK_MSG(4822,(/""/),(/0.d0/))
  nerr=nerr+1
  GOTO 1000
ENDIF
!
!
!
400 CONTINUE
!Atom positions are in reduced coordinates => convert to cartesian coordinates
CALL FRAC2CART(P,H)
IF( ALLOCATED(S) .AND. SIZE(S,1)>0 ) THEN
  CALL FRAC2CART(S,H)
ENDIF
!
!Apply options to the final system
CALL OPTIONS_AFF(options_array,Huc,H,P,S,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
IF(nerr>0) GOTO 1000
!
!
!
500 CONTINUE
WRITE(temp,*) Nfiles
 comment(1) = "# Positions averaged over "//TRIM(ADJUSTL(temp))//" files"
!Output final system to file(s)
CALL WRITE_AFF(outputfile,outfileformats,H,P,S,comment,AUXNAMES,AUX)
!
!
!
1000 CONTINUE
!
!
END SUBROUTINE AVERAGE_XYZ
!
!
END MODULE mode_average
