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
!* Last modification: P. Hirel - 26 March 2025                                    *
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
SUBROUTINE MERGE_XYZ(merge_files,merge_stack,merge_scale,options_array,outputfile,outfileformats)
!
!
IMPLICIT NONE
CHARACTER(LEN=1),INTENT(IN):: merge_stack  !if systems must be concatenated another along X, Y or Z
CHARACTER(LEN=3),INTENT(IN):: merge_scale  !if systems must match dimension along X, Y, Z, XY, YZ, XZ, XYZ
CHARACTER(LEN=*),DIMENSION(:),INTENT(IN):: merge_files
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats
CHARACTER(LEN=4096):: msg
CHARACTER(LEN=4096):: outputfile
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES, currAUXNAMES, tempAUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment, currcomment, tempcomment
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array
LOGICAL:: auxexists  !does current auxiliary property already exist in AUX?
LOGICAL:: stack !must systems be stacked?
LOGICAL:: match !must systems be rescaled to match size of 1st system?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: a1   !direction along which the files must be concatenated
INTEGER:: a2, a3 !if system lengths must be rescaled to match along other directions
INTEGER:: auxcol
INTEGER:: i, j, k
INTEGER:: NP       !number of atoms in current system
INTEGER:: Nfiles  !number of files merged
INTEGER:: oldsize
INTEGER:: sysID !column in AUX for storing system ID
REAL(dp):: def  !deformation factor to apply
REAL(dp):: disp !displacement vector
REAL(dp),DIMENSION(3,3):: H, Htemp      !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT        !Crystallographic orientation of the system
REAL(dp),DIMENSION(9,9):: C_tensor  !elastic tensor
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P, S !Coordinates read from one file (atoms/shells)
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Q, T !Final coordinates of atoms/shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: R, U !Temporary arrays storing coordinates of atoms/shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX, currAUX, tempAUX  !auxiliary properties of atoms
!
!
!Initialize variables
a1 = 0
a2 = 0
a3 = 0
stack = .FALSE.
match = .FALSE.
Nfiles = 0
sysID = 0
IF(ALLOCATED(SELECT)) DEALLOCATE(SELECT)
IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
IF(ALLOCATED(Q)) DEALLOCATE(Q)
IF(ALLOCATED(T)) DEALLOCATE(T)
IF(ALLOCATED(R)) DEALLOCATE(R)
IF(ALLOCATED(U)) DEALLOCATE(U)
IF(ALLOCATED(comment)) DEALLOCATE(comment)
def = 0.d0
ORIENT(:,:) = 0.d0
 C_tensor(:,:) = 0.d0
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
SELECT CASE( StrUpCase(merge_stack) )
CASE('X')
  stack = .TRUE.
  a1 = 1
CASE('Y')
  stack = .TRUE.
  a1 = 2
CASE('Z')
  stack = .TRUE.
  a1 = 3
END SELECT
!
IF(stack) THEN
  j=1
ELSE
  j=-1
ENDIF
!
SELECT CASE( StrUpCase(merge_scale) )
CASE('X')
  a2 = 1
CASE('Y')
  a2 = 2
CASE('Z')
  a2 = 3
CASE('XY','YX')
  a2 = 1
  a3 = 2
CASE('XZ')
  a2 = 1
  a3 = 3
CASE('YZ','ZY')
  a2 = 1
  a3 = 3
CASE('XYZ')
  a2 = 4
END SELECT
!
IF( a2>0 .OR. a3>0 ) THEN
  match = .TRUE.
ENDIF
!
CALL ATOMSK_MSG(4030,(/merge_stack,merge_scale/),(/DBLE(j)/))
!
!
!
200 CONTINUE
!Read files
DO i=1,SIZE(merge_files)
  !
  WRITE(msg,*) "MODE MERGE: reading file #", i, " : ", TRIM(ADJUSTL(merge_files(i)))
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  !
  !Read the current file
  CALL READ_AFF(merge_files(i),Htemp,P,S,currcomment,currAUXNAMES,currAUX)
  IF(nerr>0) GOTO 800
  !
  IF( verbosity==4 ) THEN
    !Print some debug messages
    msg = 'Size of arrays for file: '//TRIM(ADJUSTL(merge_files(i)))
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
      WRITE(msg,*) "currAUX: ", SIZE(currAUX,1), SIZE(currAUX,2)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDIF
  ENDIF
  !
  IF(i==1) THEN
    !First file: allocate arrays and fill them
    IF(ALLOCATED(Q)) DEALLOCATE(Q)
    ALLOCATE( Q( SIZE(P,1), SIZE(P,2) ) )
    Q = P
    IF(ALLOCATED(S)) THEN
      IF(ALLOCATED(T)) DEALLOCATE(T)
      ALLOCATE( T( SIZE(S,1), SIZE(S,2) ) )
      T = S
    ENDIF
    H = Htemp
    IF( ALLOCATED(currcomment) ) THEN
      IF(ALLOCATED(comment)) DEALLOCATE(comment)
      ALLOCATE(comment(SIZE(currcomment)))
      comment(:) = currcomment(:)
    ENDIF
    !Create array AUX to save auxiliary properties
    !We add one column to store the system ID ("sysID")
    IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
    IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
    IF( ALLOCATED(currAUXNAMES) .AND. SIZE(currAUXNAMES)>0 ) THEN
      ALLOCATE( AUXNAMES(SIZE(currAUXNAMES)+1) )
      sysID = SIZE(AUXNAMES)
      ALLOCATE( AUX( SIZE(currAUX,1),SIZE(currAUX,2)+1 ) )
      AUX(:,:) = 0.d0
      DO j=1,SIZE(currAUXNAMES)
        AUXNAMES(j) = currAUXNAMES(j)
      ENDDO
      DO j=1,SIZE(currAUX,1)
        DO k=1,SIZE(currAUX,2)
          AUX(j,k) = currAUX(j,k)
        ENDDO
      ENDDO
    ELSE
      ALLOCATE( AUXNAMES(1) )
      sysID = 1
      ALLOCATE( AUX( SIZE(P,1) , 1 ) )
      AUX(:,:) = 0.d0
    ENDIF
    AUXNAMES(sysID) = "sysID"
    AUX(:,sysID) = 1.d0
    !
  ELSE
    !Next files: merge them with the previous system H, Q
    NP = SIZE(P,1)  !number of atoms in current system
    !Resize array R to store arrays Q and P
    CALL RESIZE_DBLEARRAY2( R , SIZE(Q,1)+NP , SIZE(Q,2) , k )
    IF( k.NE.0 ) THEN
      CALL ATOMSK_MSG(818,(/"R"/),(/0.d0/))
      nerr = nerr+1
      GOTO 1000
    ENDIF
    !
    IF( ALLOCATED(S) ) THEN
      CALL RESIZE_DBLEARRAY2( U , SIZE(T,1)+SIZE(S,1) , SIZE(S,2) , k )
      IF( k.NE.0 ) THEN
        CALL ATOMSK_MSG(818,(/"U"/),(/0.d0/))
        nerr = nerr+1
        GOTO 1000
      ENDIF
    ENDIF
    !
    IF( match ) THEN
      !Match the size of current system along a2, and a3 if relevant
      WRITE(msg,*)  "MODE MERGE: matching system size along dimensions ", a2, a3
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      DO k=1,3 !loop on X, Y, Z
        IF( a2==k .OR. a3==k .OR. a2>=4 ) THEN
          !Compute deformation along this direction
          def = ( H(k,k)-Htemp(k,k) ) / Htemp(k,k)
          DO j=1,SIZE(P,1)  !loop on all atoms
            !Compute displacement
            disp = P(j,k) * def
            !Apply displacement
            P(j,k) = P(j,k) + disp
            IF(ALLOCATED(S)) THEN
              !Apply same displacement to shell
              S(j,k) = S(j,k) + disp
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    !
    IF(stack) THEN
      !Stack current system on top of previous ones along direction a1
      WRITE(msg,*)  "MODE MERGE: stacking systems along dimension ", a1
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      !The coordinates of P and S must be shifted by H(a1,:)
      DO j=1,NP
        P(j,:) = P(j,:) + H(a1,:)
      ENDDO
      IF(ALLOCATED(S)) THEN
        DO j=1,SIZE(S,1)
          S(j,:) = S(j,:) + H(a1,:)
        ENDDO
      ENDIF
      !
      !Final box must be extended along a1
      H(a1,:) = H(a1,:) + Htemp(a1,:)
      !
      WRITE(msg,*)  "MODE MERGE: new system size ", H(a1,:)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    ENDIF
    !
    !Copy atom coordinates in R
    WRITE(msg,*)  "MODE MERGE: merging Q and P into R"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    DO j=1,SIZE(Q,1)
      R(j,:) = Q(j,:)
    ENDDO
    DO j=1,NP
      R(SIZE(Q,1)+j,:) = P(j,:)
    ENDDO
    !Copy shells in U
    IF(ALLOCATED(S)) THEN
      WRITE(msg,*)  "MODE MERGE: merging T and S into U"
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      DO j=1,SIZE(T,1)
        U(j,:) = T(j,:)
      ENDDO
      DO j=1,SIZE(S,1)
        U(SIZE(T,1)+j,:) = S(j,:)
      ENDDO
    ENDIF
    !
    !Replace old Q with new R
    WRITE(msg,*)  "MODE MERGE: replacing old Q with new R"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    CALL RESIZE_DBLEARRAY2( Q , SIZE(R,1) , SIZE(R,2) , k )
    IF( k.NE.0 ) THEN
      CALL ATOMSK_MSG(818,(/"Q"/),(/0.d0/))
      nerr = nerr+1
      GOTO 1000
    ENDIF
    !IF(ALLOCATED(Q)) DEALLOCATE(Q)
    !ALLOCATE( Q( SIZE(R,1), SIZE(R,2) ) )
    Q = R
    !Shells: replace old T with new U
    IF(ALLOCATED(S)) THEN
      WRITE(msg,*)  "MODE MERGE: replacing old T with new U"
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      CALL RESIZE_DBLEARRAY2( T , SIZE(U,1) , SIZE(U,2) , k )
      IF( k.NE.0 ) THEN
        CALL ATOMSK_MSG(818,(/"T"/),(/0.d0/))
        nerr = nerr+1
        GOTO 1000
      ENDIF
      !IF(ALLOCATED(T)) DEALLOCATE(T)
      !ALLOCATE( T( SIZE(U,1), SIZE(U,2) ) )
      T = U
      IF(ALLOCATED(U)) DEALLOCATE(U)
      IF(ALLOCATED(S)) DEALLOCATE(S)
    ENDIF
    IF(ALLOCATED(R)) DEALLOCATE(R)
    !we don't need P anymore
    IF(ALLOCATED(P)) DEALLOCATE(P)
    !
    !If comments exist, append them to the "comment" array
    IF( ALLOCATED(currcomment) .AND. SIZE(currcomment)>0 ) THEN
      IF(ALLOCATED(tempcomment)) DEALLOCATE(tempcomment)
      IF(ALLOCATED(comment)) THEN
        k = SIZE(comment)
      ELSE
        k = 0
      ENDIF
      ALLOCATE(tempcomment(k+SIZE(currcomment)))
      IF(ALLOCATED(comment)) THEN
        DO j=1,k
          tempcomment(j) = comment(j)
        ENDDO
      ENDIF
      DO j=1,SIZE(currcomment)
        tempcomment(k+j) = currcomment(j)
      ENDDO
      IF(ALLOCATED(comment)) DEALLOCATE(comment)
      ALLOCATE(comment(SIZE(tempcomment)))
      comment(:) = tempcomment(:)
      DEALLOCATE(tempcomment)
    ENDIF
    !
    !Extend array AUX to fit the number of particles in Q
    !NOTE: at this point AUX *must have been allocated* while reading first system,
    !      so we know that there is at least 1 aux.prop. (sysID).
    !For now all auxiliary properties of new atoms are assigned zero values.
    !If current system #i has aux.prop. (in currAUX), they will be copied later
    oldsize = SIZE(AUX,1)
    WRITE(msg,*)  "MODE MERGE: resizing AUX"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    CALL RESIZE_DBLEARRAY2( AUX , SIZE(Q,1) , SIZE(AUXNAMES) , k )
    IF( k.NE.0 ) THEN
      CALL ATOMSK_MSG(818,(/"AUX"/),(/0.d0/))
      nerr = nerr+1
      GOTO 1000
    ENDIF
    !
    !Add current system ID
    DO j=1,NP
      AUX(oldsize+j,sysID) = DBLE(i)
    ENDDO
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
          IF( auxexists .AND. auxcol>0 ) THEN
            WRITE(msg,*) 'Property already exists: ', TRIM(ADJUSTL(currAUXNAMES(j)))
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
            !This property already has a column in AUX => complete it
            DO k=1,SIZE(currAUX,1)
              AUX(oldsize+k,auxcol) = currAUX(k,j)
            ENDDO
          ELSE
            WRITE(msg,*) 'New property: ', TRIM(ADJUSTL(currAUXNAMES(j)))
            CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
            !This property does not exist in AUX yet => create a new column for it
            auxcol = SIZE(AUXNAMES)+1
            ALLOCATE(tempAUXNAMES(auxcol))
            DO k=1,SIZE(AUXNAMES)
              tempAUXNAMES(k) = AUXNAMES(k)
            ENDDO
            tempAUXNAMES(auxcol) = currAUXNAMES(j)
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
              tempAUX(oldsize+k,auxcol) = currAUX(k,j)
            ENDDO
            DEALLOCATE(AUX)
            ALLOCATE( AUX( SIZE(tempAUX,1),SIZE(tempAUX,2) ) )
            AUX(:,:) = tempAUX(:,:)
            DEALLOCATE(tempAUX)
          ENDIF
        ENDDO
        !
      ELSE
        WRITE(msg,*) 'First system with auxiliary properties: allocating AUX'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        !No auxiliary property existed before, the current ones are the first ones
        !=> define arrays and fill them
        IF( ALLOCATED(AUXNAMES) ) DEALLOCATE(AUXNAMES)
        IF( ALLOCATED(AUX) ) DEALLOCATE(AUX)
        ALLOCATE( AUXNAMES(SIZE(currAUXNAMES)) )
        AUXNAMES(:) = currAUXNAMES(:)
        ALLOCATE( AUX( SIZE(Q,1),SIZE(currAUX,2) ) )
        AUX(:,:) = 0.d0
        DO k=1,SIZE(currAUX,1)
          AUX(SIZE(AUX,1)-SIZE(currAUX,1)+k,:) = currAUX(k,:)
        ENDDO
      ENDIF  !endif AUXNAMES
      !
    ENDIF !endif currAUX
    !
  ENDIF !Endif i==1
  !
  IF( verbosity==4 ) THEN
    WRITE(msg,*) "Finished treating file #", i, " , '", TRIM(ADJUSTL(merge_files(i))), "'"
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
    WRITE(msg,*) "   =  =  =  =  =  =  =  =  =  ="
    CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!     DO j=1,SIZE(AUX,1)
!       WRITE(msg,*) AUX(j,:)
!       CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!     ENDDO
  ENDIF
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
CALL OPTIONS_AFF(options_array,Htemp,H,Q,T,AUXNAMES,AUX,ORIENT,SELECT,C_tensor)
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
