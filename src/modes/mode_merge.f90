MODULE mode_merge
!
!**********************************************************************************
!*  MODE_MERGE                                                                    *
!**********************************************************************************
!* This module reads two or more sets of atomic coordinates and                   *
!* merges them, i.e. unites all atomic positions in one single array.             *
!**********************************************************************************
!* (C) March 2010 - Pierre Hirel                                                  *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
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
USE subroutines
USE deterh
USE messages
USE files
USE readin
USE writeout
USE options
!
!
CONTAINS
!
!
SUBROUTINE MERGE_XYZ(merge_files,merge_dir,options_array,outputfile,outfileformats)
!
!
IMPLICIT NONE
CHARACTER(LEN=1),INTENT(IN):: merge_dir  !if systems must be concatenated another along X, Y or Z
CHARACTER(LEN=*),DIMENSION(:),INTENT(IN):: merge_files
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats
CHARACTER(LEN=128):: msg
CHARACTER(LEN=4096):: outputfile
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES, currAUXNAMES, tempAUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment, tempcomment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array
LOGICAL:: auxexists  !does current auxiliary property already exist in AUX?
LOGICAL:: cat  !must the files be concatenated?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: a1   !direction along which the files must be concatenated
INTEGER:: auxcol
INTEGER:: i, j, k
INTEGER:: Nfiles  !number of files merged
REAL(dp),DIMENSION(3,3):: H, Htemp      !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT        !Crystallographic orientation of the system
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, S !Coordinates read from one file (atoms/shells)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Q, T !Final coordinates of atoms/shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: R, U !Temporary arrays storing coordinates of atoms/shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX, currAUX, tempAUX  !auxiliary properties of atoms
!
!
!Initialize variables
a1 = 0
Nfiles = 0
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(T)) DEALLOCATE(T)
IF(ALLOCATED(R)) DEALLOCATE(R)
IF(ALLOCATED(U)) DEALLOCATE(U)
ALLOCATE(comment(1))
 comment(1)=''
ORIENT(:,:) = 0.d0
!
!
100 CONTINUE
IF( SIZE(merge_files)==0 ) THEN
  nerr = nerr+1
  CALL ATOMSK_MSG(4806,(/''/),(/0.d0/))
  GOTO 1000
ENDIF
!
IF( verbosity==4 ) THEN
  msg = 'Entering MERGE_XYZ, list of files:'
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  DO i=1,SIZE(merge_files)
    msg = TRIM(ADJUSTL(merge_files(i)))
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  ENDDO
ENDIF
!
IF(merge_dir=='x' .OR. merge_dir=='X') THEN
  cat = .TRUE.
  a1 = 1
ELSEIF(merge_dir=='y' .OR. merge_dir=='Y') THEN
  cat = .TRUE.
  a1 = 2
ELSEIF(merge_dir=='z' .OR. merge_dir=='Z') THEN
  cat = .TRUE.
  a1 = 3
ELSE
  cat = .FALSE.
ENDIF
!
IF(cat) THEN
  j=1
ELSE
  j=-1
ENDIF
CALL ATOMSK_MSG(4030,(/merge_dir/),(/DBLE(j)/))
!
!
!
200 CONTINUE
!Read files
DO i=1,SIZE(merge_files)
  !Read the current file
  CALL READ_AFF(merge_files(i),Htemp,P,S,tempcomment,currAUXNAMES,currAUX)
  IF(nerr>0) GOTO 800
  !
  IF( verbosity==4 ) THEN
    !Print some debug messages
    msg = 'Size of arrays for file:'//TRIM(ADJUSTL(merge_files(i)))
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,*) "P: ", SIZE(P,1), SIZE(P,2)
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    IF( ALLOCATED(S) ) THEN
      WRITE(msg,*) "S: ", SIZE(S,1), SIZE(S,2)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDIF
    IF( ALLOCATED(currAUX) ) THEN
      WRITE(msg,*) "currAUXNAMES: ", SIZE(currAUXNAMES)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      WRITE(msg,*) "currAUX: ", SIZE(AUX,1), SIZE(AUX,2)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDIF
  ENDIF
  !
  IF(i==1) THEN
    !First file: allocate arrays and fill them
    ALLOCATE( Q( SIZE(P,1), SIZE(P,2) ) )
    Q = P
    IF(ALLOCATED(S)) THEN
      ALLOCATE( T( SIZE(S,1), SIZE(S,2) ) )
      T = S
    ENDIF
    H = Htemp
    IF( ALLOCATED(tempcomment) ) THEN
      comment(1) = tempcomment(1)
    ELSE
      comment(1) = ''
    ENDIF
    IF( ALLOCATED(currAUXNAMES) .AND. SIZE(currAUXNAMES)>0 ) THEN
      ALLOCATE( AUXNAMES(SIZE(currAUXNAMES)) )
      AUXNAMES(:) = currAUXNAMES(:)
      ALLOCATE( AUX( SIZE(currAUX,1),SIZE(currAUX,2) ) )
      AUX(:,:) = currAUX(:,:)
    ENDIF
    !
  ELSE
    !Next files: merge them with the previous system H, Q
    !Create array R to store arrays Q and P
    IF(ALLOCATED(R)) DEALLOCATE(R)
    ALLOCATE( R( SIZE(Q,1)+SIZE(P,1), SIZE(Q,2) ) )
    IF(ALLOCATED(S)) THEN
      IF(ALLOCATED(U)) DEALLOCATE(U)
      ALLOCATE( U( SIZE(T,1)+SIZE(S,1), SIZE(S,2) ) )
    ENDIF
    !
    IF(cat) THEN
      !Concatenate the files
      !The coordinates of P and S must be shifted by H(a1,:)
      DO j=1,SIZE(P,1)
        P(j,:) = P(j,:) + H(a1,:)
      ENDDO
      IF(ALLOCATED(S)) THEN
        DO j=1,SIZE(S,1)
          S(j,:) = S(j,:) + H(a1,:)
        ENDDO
      ENDIF
      !
      !The box must be extended along a1
      H(a1,:) = H(a1,:) + Htemp(a1,:)
    ENDIF
    !
    !Copy atom coordinates in R
    DO j=1,SIZE(Q,1)
      R(j,:) = Q(j,:)
    ENDDO
    DO j=1,SIZE(P,1)
      R(SIZE(Q,1)+j,:) = P(j,:)
    ENDDO
    !Copy shells in U
    IF(ALLOCATED(S)) THEN
      DO j=1,SIZE(T,1)
        U(j,:) = T(j,:)
      ENDDO
      DO j=1,SIZE(S,1)
        U(SIZE(T,1)+j,:) = S(j,:)
      ENDDO
    ENDIF
    !
    !Replace old Q with new R
    DEALLOCATE(Q)
    ALLOCATE( Q( SIZE(R,1), SIZE(R,2) ) )
    Q = R
    !Shells: replace old T with new U
    IF(ALLOCATED(S)) THEN
      DEALLOCATE(T)
      ALLOCATE( T( SIZE(U,1), SIZE(U,2) ) )
      T = U
      DEALLOCATE(U)
      DEALLOCATE(S)
    ENDIF
    DEALLOCATE(R)
    !we don't need P anymore
    DEALLOCATE(P)
    !
    !If AUX exists, extend it to fit the number of particles in Q
    !For now all auxiliary properties of new atoms are assigned zero values,
    !that will be corrected later if these properties are defined in currAUX
    IF( ALLOCATED(AUXNAMES) ) THEN
      ALLOCATE( tempAUX( SIZE(Q,1),SIZE(AUX,2) ) )
      tempAUX(:,:) = 0.d0
      DO j=1,SIZE(AUX,1)
        tempAUX(j,:) = AUX(j,:)
      ENDDO
      DEALLOCATE(AUX)
      ALLOCATE( AUX( SIZE(tempAUX,1),SIZE(tempAUX,2) ) )
      AUX(:,:) = tempAUX(:,:)
      DEALLOCATE(tempAUX)
    ENDIF
    !
    !If auxiliary properties were read from current file, merge them with existing ones
    IF( ALLOCATED(currAUXNAMES) .AND. ALLOCATED(currAUX) ) THEN
      !
      IF( ALLOCATED(AUXNAMES) .AND. SIZE(AUXNAMES)>0 ) THEN
        !Some auxiliary properties already existed
        !Check if some property names already exist in AUXNAMES
        DO j=1,SIZE(currAUXNAMES)
          auxexists=.FALSE.
          auxcol=0
          k=0
          DO WHILE ( k<=SIZE(AUXNAMES) .AND. auxcol==0 )
            k=k+1
            IF( TRIM(ADJUSTL(currAUXNAMES(j)))==TRIM(ADJUSTL(AUXNAMES(k))) ) THEN
              auxexists=.TRUE.
              auxcol=k
            ENDIF
          ENDDO
          !
          IF( auxexists ) THEN
            WRITE(msg,*) 'Property already exists: ', TRIM(ADJUSTL(currAUXNAMES(j)))
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
            !This property already has a column in AUX => complete it
            DO k=1,SIZE(currAUX,1)
              AUX(SIZE(AUX,1)-SIZE(currAUX,1)+k,auxcol) = currAUX(k,auxcol)
            ENDDO
          ELSE
            WRITE(msg,*) 'New property: ', TRIM(ADJUSTL(currAUXNAMES(j)))
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
            !This property does not exist in AUX yet => create a new column for it
            ALLOCATE(tempAUXNAMES(SIZE(AUXNAMES)+1))
            DO k=1,SIZE(AUXNAMES)
              tempAUXNAMES(k) = AUXNAMES(k)
            ENDDO
            tempAUXNAMES(SIZE(AUXNAMES)+1) = currAUXNAMES(j)
            DEALLOCATE(AUXNAMES)
            ALLOCATE( AUXNAMES( SIZE(tempAUXNAMES) ) )
            AUXNAMES(:) = tempAUXNAMES(:)
            DEALLOCATE(tempAUXNAMES)
            !
            ALLOCATE( tempAUX( SIZE(AUX,1) , SIZE(AUXNAMES) ) )
            tempAUX(:,:) = 0.d0
            DO k=1,SIZE(AUX,2)
              tempAUX(:,k) = AUX(:,k)
            ENDDO
            DO k=1,SIZE(currAUX,1)
              tempAUX(SIZE(AUX,1)-SIZE(currAUX,1)+k,SIZE(AUXNAMES)) = currAUX(k,j)
            ENDDO
            DEALLOCATE(AUX)
            ALLOCATE( AUX( SIZE(tempAUX,1),SIZE(tempAUX,2) ) )
            AUX(:,:) = tempAUX(:,:)
            DEALLOCATE(tempAUX)
          ENDIF
        ENDDO
        !
      ELSE
        IF( ALLOCATED(AUXNAMES) ) DEALLOCATE(AUXNAMES)
        IF( ALLOCATED(AUX) ) DEALLOCATE(AUX)
        WRITE(msg,*) 'Allocating AUX'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        !No auxiliary property existed before, the current ones are the first ones
        !=> define arrays and fill them
        ALLOCATE( AUXNAMES(SIZE(currAUXNAMES)) )
        AUXNAMES(:) = currAUXNAMES(:)
        ALLOCATE( AUX( SIZE(Q,1),SIZE(currAUX,2) ) )
        AUX(:,:) = 0.d0
        DO k=1,SIZE(currAUX,1)
          AUX(SIZE(AUX,1)-SIZE(currAUX,1)+k,:) = currAUX(k,:)
        ENDDO
      ENDIF
    ENDIF
    !
  ENDIF !Endif i==1
  !
  Nfiles = Nfiles+1
ENDDO
!
!
!
400 CONTINUE
CALL ATOMSK_MSG(4031,(/''/),(/DBLE(Nfiles)/))
!
!Apply options to the final system
Htemp(:,:) = 0.d0
CALL OPTIONS_AFF(options_array,Htemp,H,Q,T,AUXNAMES,AUX,ORIENT,SELECT)
!
!Write final system into file(s)
CALL WRITE_AFF(outputfile,outfileformats,H,Q,T,comment,AUXNAMES,AUX)
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(1801,(/TRIM(merge_files(i))/),(/0.d0/))
!
!
!
1000 CONTINUE
IF(ALLOCATED(R)) DEALLOCATE(R)
!
!
END SUBROUTINE MERGE_XYZ
!
!
END MODULE mode_merge
