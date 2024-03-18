MODULE mode_interpolate
!
!**********************************************************************************
!*  MODE_INTERPOLATE                                                              *
!**********************************************************************************
!* This module reads two files containing atom positions of similar systems,      *
!* and generates N new configurations by interpolating atom positions.            *
!**********************************************************************************
!* (C) July 2013 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 18 March 2024                                    *
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
USE messages
USE files
USE subroutines
!Modules to read and write files
USE readin
USE writeout
!Module to apply options
USE options
!
!
CONTAINS
!
!
SUBROUTINE INTERPOLATE_XYZ(file1,file2,Nimages,prefix,outfileformats,options_array)
!
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: file1, file2
CHARACTER(LEN=*),INTENT(IN):: prefix
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats
CHARACTER(LEN=4096):: msg, outputfile
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES1, AUXNAMES2 !names of auxiliary properties of ini.&final images
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMESimg          !names of auxiliary properties of an image
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment, comment1, comment2
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
LOGICAL:: doaux    !take auxiliary properties into account?
LOGICAL:: doshells !also take shells into account?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: i, j, k
INTEGER,INTENT(IN):: Nimages  !number of configurations to interpolate (that excludes initial & final images)
INTEGER,DIMENSION(:,:),ALLOCATABLE:: AUXid !index of matching aux.prop. in AUX1, AUX2
REAL(dp),DIMENSION(3,3):: Huc   !Box vectors of unit cell (unknown, set to 0 here)
REAL(dp),DIMENSION(3,3):: H1, H2   !Base vectors of the supercells of initial, final images
REAL(dp),DIMENSION(3,3):: Himg     !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystal orientation
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX1, AUX2 !auxiliary properties of initial, final images
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P1, P2  !positions of initial, final images
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S1, S2  !positions of initial, final shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Pimg, Simg !positions of atoms, shells in an image
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUXimg     !auxiliary properties of an image
!
!
!Initialize variables
doaux = .FALSE.
doshells = .FALSE.
IF(ALLOCATED(Pimg)) DEALLOCATE(Pimg)
IF(ALLOCATED(Simg)) DEALLOCATE(Simg)
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
IF(ALLOCATED(AUXid)) DEALLOCATE(AUXid)
C_tensor(:,:) = 0.d0
ORIENT(:,:) = 0.d0
!
!
CALL ATOMSK_MSG(4059,(/''/),(/0.d0/))
!
!
100 CONTINUE
!Read files
CALL READ_AFF(file1,H1,P1,S1,comment1,AUXNAMES1,AUX1)
CALL READ_AFF(file2,H2,P2,S2,comment2,AUXNAMES2,AUX2)
!
!Check if the number of atoms coincide
IF( SIZE(P1,1) .NE. SIZE(P2,1) ) THEN
  RETURN
ENDIF
!
!If both configurations have shells, their positions will be interpolated
IF( ALLOCATED(S1) .AND. ALLOCATED(S2) ) THEN
  IF( SIZE(S1,1) == SIZE(S2,1) ) THEN
    doshells = .TRUE.
  ENDIF
ENDIF
!
IF( .NOT. doshells ) THEN
  IF(ALLOCATED(S1)) DEALLOCATE(S1)
  IF(ALLOCATED(S2)) DEALLOCATE(S2)
ENDIF
!
!Deal with auxiliary properties
IF( ALLOCATED(AUXNAMES1) .AND. ALLOCATED(AUXNAMES2) ) THEN
  !Auxiliary properties are defined in both systems
  !Construct a table matching identical properties from both systems
  ALLOCATE( AUXid( MIN(SIZE(AUXNAMES1),SIZE(AUXNAMES2)) , 2 ) )
  AUXid(:,:) = 0
  k=0
  DO i=1,SIZE(AUXNAMES1)
    DO j=1,SIZE(AUXNAMES2)
      IF( AUXNAMES1(i) == AUXNAMES2(j) ) THEN
        k = k+1
        AUXid(k,1) = i
        AUXid(k,2) = j
      ENDIF
    ENDDO
  ENDDO
  ALLOCATE( AUXNAMESimg(k) )
  AUXNAMESimg(:) = ""
  k=0
  DO i=1,SIZE(AUXNAMES1)
    DO j=1,SIZE(AUXNAMES2)
      IF( AUXNAMES1(i) == AUXNAMES2(j) ) THEN
        k = k+1
        AUXNAMESimg(k) = TRIM(AUXNAMES1(i))
      ENDIF
    ENDDO
  ENDDO
  IF( k>0 ) THEN
    !Matching aux.prop. were found, and their names are now in AUXNAMESimg(:)
    doaux = .TRUE.
    IF(verbosity>=4) THEN
      WRITE(msg,*) "Found ", k, " matching properties:"
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      DO j=1,SIZE(AUXNAMESimg)
        WRITE(msg,*) j, "  ", TRIM(AUXNAMESimg(j))
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDDO
    ENDIF
  ELSE
    !Found zero matching aux.prop. between the two systems: ignore aux.prop.
    doaux = .FALSE.
    IF(ALLOCATED(AUXNAMESimg)) DEALLOCATE(AUXNAMESimg)
  ENDIF
ENDIF
!
IF(.NOT.doaux) THEN
  !Ignore auxiliary properties
  IF(ALLOCATED(AUX1)) DEALLOCATE(AUX1)
  IF(ALLOCATED(AUX2)) DEALLOCATE(AUX2)
  IF(ALLOCATED(AUXimg)) DEALLOCATE(AUXimg)
  IF(ALLOCATED(AUXNAMES1)) DEALLOCATE(AUXNAMES1)
  IF(ALLOCATED(AUXNAMES2)) DEALLOCATE(AUXNAMES2)
  IF(ALLOCATED(AUXNAMESimg)) DEALLOCATE(AUXNAMESimg)
ENDIF
!
!
!
200 CONTINUE
!Loop on all interpolated images
DO i=0,Nimages+1
  CALL ATOMSK_MSG(4060,(/''/),(/DBLE(i)/))
  !
  IF( i==0 .OR. i==Nimages+1 ) THEN
    !Initial or final configuration
    ALLOCATE( Pimg( SIZE(P1,1) , 4) )
    Pimg(:,:) = 0.d0
    !Copy positions of relevant configuration
    IF( i==0 ) THEN
      Himg = H1
      Pimg = P1
    ELSE
      Himg = H2
      Pimg = P2
    ENDIF
    !Same with positions of shells, if any
    IF( doshells ) THEN
      ALLOCATE( Simg(SIZE(P1,1) , 4) )
      IF( i==0 ) THEN
        Simg = S1
      ELSE
        Simg = S2
      ENDIF
    ENDIF
    !Same with aux.prop, if any
    IF( doaux ) THEN
      ALLOCATE( AUXimg(SIZE(P1,1) , SIZE(AUXNAMESimg) ) )
      IF( i==0 ) THEN
        DO j=1,SIZE(AUXNAMESimg)
          DO k=1,SIZE(AUXNAMES1)
            IF( AUXNAMESimg(j)==AUXNAMES1(k) ) THEN
              AUXimg(:,j) = AUX1(:,k)
            ENDIF
          ENDDO
        ENDDO
      ELSE
        DO j=1,SIZE(AUXNAMESimg)
          DO k=1,SIZE(AUXNAMES2)
            IF( AUXNAMESimg(j)==AUXNAMES2(k) ) THEN
              AUXimg(:,j) = AUX2(:,k)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF
    !
  ELSE
    !Compute interpolated box vectors
    DO j=1,3
      DO k=1,3
        Himg(j,k) = H1(j,k) + (H2(j,k)-H1(j,k))*DBLE(i)/DBLE(Nimages)
      ENDDO
    ENDDO
    !
    !Compute interpolated positions, save them into Pimg
    ALLOCATE( Pimg( SIZE(P1,1) , 4) )
    Pimg(:,:) = 0.d0
    DO j=1,SIZE(P1,1)
      Pimg(j,1:3) = P1(j,1:3) + (P2(j,1:3)-P1(j,1:3))*DBLE(i)/DBLE(Nimages)
      Pimg(j,4) = P1(j,4)
    ENDDO
    !
    !Same with positions of shells, if any
    IF( doshells ) THEN
      ALLOCATE( Simg(SIZE(P1,1) , 4) )
      Simg(:,:) = 0.d0
      DO j=1,SIZE(S1,1)
        Simg(j,1:3) = S1(j,1:3) + (S2(j,1:3)-S1(j,1:3))*DBLE(i)/DBLE(Nimages)
        Simg(j,4) = S1(j,4)
      ENDDO
    ENDIF
    !
    !Interpolate matching aux.prop., save them into AUXimg
    IF( doaux ) THEN
      ALLOCATE( AUXimg(SIZE(P1,1),SIZE(AUXNAMESimg) ) )
      AUXimg(:,:) = 0.d0
      DO j=1,SIZE(AUXimg,2)
        AUXimg(:,j) = AUX1(:,AUXid(j,1)) + (AUX2(:,AUXid(j,2))-AUX1(:,AUXid(j,1)))*DBLE(i)/DBLE(Nimages)
      ENDDO
    ENDIF
    !
  ENDIF
  !
  !Apply options to current image
  CALL OPTIONS_AFF(options_array,Huc,Himg,Pimg,Simg,AUXNAMESimg,AUXimg,ORIENT,SELECT,C_tensor)
  !
  !Write current image to file(s)
  ALLOCATE(comment(1))
  WRITE(comment(1),*) i
  outputfile = TRIM(ADJUSTL(prefix))//"_img"//TRIM(ADJUSTL(comment(1)))
  comment(1) = "Image #"//TRIM(ADJUSTL(comment(1)))
  CALL WRITE_AFF(outputfile,outfileformats,Himg,Pimg,Simg,comment,AUXNAMESimg,AUXimg)
  !
  !Free memory
  IF(ALLOCATED(comment)) DEALLOCATE(comment)
  IF(ALLOCATED(Pimg)) DEALLOCATE(Pimg)
  IF(ALLOCATED(Simg)) DEALLOCATE(Simg)
  IF(ALLOCATED(AUXimg)) DEALLOCATE(AUXimg)
ENDDO
!
!
!
1000 CONTINUE
IF(ALLOCATED(P1)) DEALLOCATE(P1)
IF(ALLOCATED(S1)) DEALLOCATE(S1)
IF(ALLOCATED(P2)) DEALLOCATE(P2)
IF(ALLOCATED(S2)) DEALLOCATE(S2)
IF(ALLOCATED(AUX1)) DEALLOCATE(AUX1)
IF(ALLOCATED(AUX2)) DEALLOCATE(AUX2)
IF(ALLOCATED(AUXid)) DEALLOCATE(AUXid)
IF(ALLOCATED(AUXNAMESimg)) DEALLOCATE(AUXNAMESimg)
!
!
END SUBROUTINE INTERPOLATE_XYZ
!
!
END MODULE mode_interpolate
