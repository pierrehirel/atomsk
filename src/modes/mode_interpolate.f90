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
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 08 Feb. 2018                                     *
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
CHARACTER(LEN=4096):: outputfile
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES1, AUXNAMES2 !names of auxiliary properties of ini.&final images
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMESimg          !names of auxiliary properties of an image
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment, comment1, comment2
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
LOGICAL:: doshells !also take shells into account?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: i, j, k
INTEGER,INTENT(IN):: Nimages  !number of configurations to interpolate (that excludes initial & final images)
REAL(dp),DIMENSION(3,3):: Huc   !Box vectors of unit cell (unknown, set to 0 here)
REAL(dp),DIMENSION(3,3):: H1, H2   !Base vectors of the supercells of initial, final images
REAL(dp),DIMENSION(3,3):: Himg     !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystal orientation
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX1, AUX2 !auxiliary properties of initial, final images
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P1, P2  !positions of initial, final images
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S1, S2  !positions of initial, final shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Pimg, Simg !positions of atoms, shells in an image
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUXimg     !auxiliary properties of an image
!
!
!Initialize variables
doshells = .FALSE.
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
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
!If one configuration has shells but not the other, ignore shells
IF( ALLOCATED(S1) .OR. ALLOCATED(S2) ) THEN
  IF( SIZE(S1,1) == SIZE(S2,1) ) THEN
    doshells = .TRUE.
  ELSE
    doshells = .FALSE.
    IF(ALLOCATED(S1)) DEALLOCATE(S1)
    IF(ALLOCATED(S2)) DEALLOCATE(S2)
  ENDIF
ENDIF
!
!Ignore auxiliary properties
IF(ALLOCATED(AUX1)) DEALLOCATE(AUX1)
IF(ALLOCATED(AUX2)) DEALLOCATE(AUX2)
IF(ALLOCATED(AUXNAMES1)) DEALLOCATE(AUXNAMES1)
IF(ALLOCATED(AUXNAMES2)) DEALLOCATE(AUXNAMES2)
!
!
!
200 CONTINUE
!Loop on all interpolated images
DO i=1,Nimages
  CALL ATOMSK_MSG(4060,(/''/),(/DBLE(i)/))
  !
  !Compute interpolated box vectors
  DO j=1,3
    DO k=1,3
      Himg(j,k) = H1(j,k) + (H2(j,k)-H1(j,k))*DBLE(i)/DBLE(Nimages)
    ENDDO
  ENDDO
  !
  !Compute interpolated positions, save them into Pimg
  IF(ALLOCATED(Pimg)) DEALLOCATE(Pimg)
  ALLOCATE( Pimg( SIZE(P1,1) , 4) )
  Pimg(:,:) = 0.d0
  DO j=1,SIZE(P1,1)
    Pimg(j,1:3) = P1(j,1:3) + (P2(j,1:3)-P1(j,1:3))*DBLE(i)/DBLE(Nimages)
    Pimg(j,4) = P1(j,4)
  ENDDO
  !
  !Same with positions of shells, if any
  IF( doshells ) THEN
    IF(ALLOCATED(Simg)) DEALLOCATE(Simg)
    ALLOCATE( Simg(SIZE(P1,1) , 4) )
    Simg(:,:) = 0.d0
    DO j=1,SIZE(S1,1)
      Simg(j,1:3) = S1(j,1:3) + (S2(j,1:3)-S1(j,1:3))*DBLE(i)/DBLE(Nimages)
      Simg(j,4) = S1(j,4)
    ENDDO
  ENDIF
  !
  !Apply options to current image
  CALL OPTIONS_AFF(options_array,Huc,Himg,Pimg,Simg,AUXNAMESimg,AUXimg,ORIENT,SELECT)
  !
  !Write current image to file(s)
  IF(ALLOCATED(comment)) DEALLOCATE(comment)
  ALLOCATE(comment(1))
  WRITE(comment(1),*) i
  outputfile = TRIM(ADJUSTL(prefix))//"_img"//TRIM(ADJUSTL(comment(1)))
  comment(1) = "Image #"//TRIM(ADJUSTL(comment(1)))
  CALL WRITE_AFF(outputfile,outfileformats,Himg,Pimg,Simg,comment,AUXNAMESimg,AUXimg)
ENDDO
!
!
!
1000 CONTINUE

!
!
END SUBROUTINE INTERPOLATE_XYZ
!
!
END MODULE mode_interpolate
