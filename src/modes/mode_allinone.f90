MODULE aio
!
!**********************************************************************************
!*  AI1_XYZ                                                                       *
!**********************************************************************************
!* This module reads many files, and writes them to a single file in              *
!* the XYZ format, or to the animated XSF format.                                 *
!* No standard specification exists, but a widely used format is                  *
!* described for instance at:                                                     *
!*     http://openbabel.org/wiki/XYZ_%28format%29                                 *
!* This module can also write to extended XYZ format, where the comment line      *
!* (i.e. the 2nd line) is replaced by a set of keywords/values. The               *
!* extended XYZ format is described for instance at:                              *
!*     http://www.jrkermode.co.uk/quippy/extended_xyz.html                        *
!* Note that this module can NOT write special XYZ format.                        *
!* The animated XSF format is described here (note that this module can only      *
!* write to the "Variable-cell animated XSF" format):                             *
!*    http://www.xcrysden.org/doc/XSF.html#__toc__10                              *
!**********************************************************************************
!* (C) December 2011 - Pierre Hirel                                               *
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
!
USE atoms
USE comv
USE constants
USE messages
USE files
USE subroutines
USE options
USE readin
!
!
CONTAINS
!
SUBROUTINE ALLINONE(listfile,outputfile,outfileformats,options_array)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: listfile, outputfile
CHARACTER(LEN=2):: species
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE,INTENT(IN):: outfileformats !list of formats to write
CHARACTER(LEN=128):: outfile_xyz, outfile_xsf
CHARACTER(LEN=128):: msg
CHARACTER(LEN=1024):: temp
CHARACTER(LEN=4096):: inputfile
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
LOGICAL:: fileexists
LOGICAL:: output_xyz, output_exyz,output_xsf
LOGICAL:: xsf_forces  !write forces in XSF format?
INTEGER:: fx, fy, fz
INTEGER:: i, j, k
INTEGER:: snap, totsnap  !current/total snapshot
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystal orientation
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX_FORCES !forces on atoms
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P,S  !positions of cores, shells
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!Initialize variables
outfile_xyz=""
output_xyz=.FALSE.
output_exyz=.FALSE.
output_xsf=.FALSE.
xsf_forces = .FALSE.
fx=0
fy=0
fz=0
snap = 0
totsnap = 0
ORIENT(:,:) = 0.d0
ALLOCATE(comment(1))
 comment(1)=''
!
!
msg = 'entering MODE ALL-IN-ONE'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
msg = "XYZ"
CALL ATOMSK_MSG(4050,(/listfile,msg,outputfile/),(/0.d0/))
!
DO i=1,SIZE(outfileformats)
  IF(outfileformats(i)=='xyz')  output_xyz = .TRUE.
  IF(outfileformats(i)=='exyz') output_exyz = .TRUE.
  IF(outfileformats(i)=='xsf')  output_xsf = .TRUE.
ENDDO
!
!
!
100 CONTINUE
OPEN(UNIT=35,FILE=listfile,FORM='FORMATTED',STATUS='OLD')
REWIND(35)
!Parse listfile to count the total number of snapshots
!Don't count lines starting with #, and ignore files that don't exist
DO
  READ(35,*,END=110,ERR=110) inputfile
  inputfile = ADJUSTL(inputfile)
  INQUIRE(FILE=inputfile,EXIST=fileexists)
  IF( inputfile(1:1).NE."#" .AND. fileexists ) THEN
    totsnap=totsnap+1
  ENDIF
ENDDO
110 CONTINUE
REWIND(35)
WRITE(msg,*) 'totsnap = ', totsnap
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!Check that number of snapshots is not zero
IF(totsnap<=0) THEN
  CALL ATOMSK_MSG(4818,(/listfile/),(/0.d0/))
  nerr = nerr+1
  GOTO 1000
ENDIF
!
!Construct file name(s) according to output format(s)
IF(output_xyz .OR. output_exyz) THEN
  !Name XYZ file
  CALL NAME_OUTFILE(outputfile,outfile_xyz,'xyz  ')
  IF(.NOT. overw) CALL CHECKFILE(outfile_xyz,'writ')
  OPEN(UNIT=40,FILE=outfile_xyz,FORM='FORMATTED',STATUS='UNKNOWN')
ENDIF
IF(output_xsf) THEN
  !Name XSF file
  CALL NAME_OUTFILE(outputfile,outfile_xsf,'xsf  ')
  IF(.NOT. overw) CALL CHECKFILE(outfile_xsf,'writ')
  OPEN(UNIT=41,FILE=outfile_xsf,FORM='FORMATTED',STATUS='UNKNOWN')
  !1st line = total number of snapshots
  WRITE(temp,*) "ANIMSTEPS ", totsnap
  WRITE(41,'(a)') TRIM(ADJUSTL(temp))
ENDIF
!IF(.FALSE.) THEN
!  CALL ATOMSK_MSG(4818,(/""/),(/0.d0/))
!  nerr=nerr+1
!  GOTO 1000
!ENDIF
!
!
!
200 CONTINUE
!Loop on all files to convert
DO i=1,totsnap
  !Initialize variables and arrays
  inputfile = ""
  H(:,:) = 0.d0
  IF(ALLOCATED(P)) DEALLOCATE(P)
  IF(ALLOCATED(S)) DEALLOCATE(S)
  IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
  IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
  IF(ALLOCATED(AUX_FORCES)) DEALLOCATE(AUX_FORCES)
  !
  !Read name of file to read
  READ(35,'(a128)',END=300,ERR=300) inputfile
  inputfile = ADJUSTL(inputfile)
  !
  !Ignore lines starting with #
  IF( inputfile(1:1).NE.'#' ) THEN
    !Check if file exists
    INQUIRE(FILE=inputfile,EXIST=fileexists)
    !
    IF(fileexists) THEN
      snap = snap+1
      WRITE(msg,*) 'Opening snapshot ', snap, TRIM(inputfile)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      !
      !Read file
      CALL READ_AFF(inputfile,H,P,S,comment,AUXNAMES,AUX)
      IF(nerr>0) GOTO 300
      !
      !Apply options if any
      CALL OPTIONS_AFF(options_array,H,P,S,AUXNAMES,AUX,ORIENT)
      IF(nerr>0) GOTO 300
      !
      !
      IF(output_xyz .OR. output_exyz) THEN
        WRITE(msg,*) 'Writing snapshot to XYZ file: ', TRIM(outfile_xyz)
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        !Write current snapshot to XYZ file
        !First line = number of particles
        WRITE(40,*) SIZE(P,1)
        !2nd line = comment or system properties
        IF(output_exyz) THEN
          !In extended xyz format the comment is replaced by properties
          msg=''
          temp=''
          DO j=1,3
            DO k=1,3
              IF(H(k,j)==0.d0) THEN
                WRITE(temp,*) '0.0'
              ELSE
                WRITE(temp,'(f16.8)') H(k,j)
              ENDIF
              WRITE(msg,*) TRIM(ADJUSTL(msg))//' '//TRIM(ADJUSTL(temp))
            ENDDO
          ENDDO
          temp = 'Lattice="'//TRIM(ADJUSTL(msg))//'" '//'Properties=species:S:1:pos:R:3'
          !If auxiliary properties are present concatenate their names here
          IF( ALLOCATED(AUXNAMES) ) THEN
            DO j=1,SIZE(AUXNAMES)
              temp = TRIM(temp)//':'//TRIM(ADJUSTL(AUXNAMES(j)))//':R:1'
            ENDDO
          ENDIF
        ELSE
          temp = comment(1)
        ENDIF
        WRITE(40,*) TRIM(ADJUSTL(temp))
        !Following lines = atom coordinates
        DO j=1,SIZE(P,1)
          temp=''
          CALL ATOMSPECIES(P(j,4),species)
          WRITE(temp,*) species
          temp = TRIM(ADJUSTL(temp))
          DO k=1,3
            WRITE(msg,'(f16.8)')  P(j,k)
            temp = TRIM(ADJUSTL(temp))//'  '//TRIM(ADJUSTL(msg))
          ENDDO
          !If auxiliary properties are present, concatenate them at end of line
          !This is done only for extended XYZ
          IF( ALLOCATED(AUXNAMES) .AND. (output_exyz) ) THEN
            DO k=1,SIZE(AUX,2)
              WRITE(msg,'(e12.5)')  AUX(j,k)
              temp = TRIM(ADJUSTL(temp))//'  '//TRIM(ADJUSTL(msg))
            ENDDO
          ENDIF
          WRITE(40,'(a)') TRIM(ADJUSTL(temp))
        ENDDO
        WRITE(msg,*) 'Snapshot written to XYZ file'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDIF
      !
      !
      IF(output_xsf) THEN
        WRITE(msg,*) 'Writing snapshot to animated XSF file: ', TRIM(outfile_xsf)
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
        !Write current snapshot to animated XSF file
        !
        !First of all, deal with forces
        IF( snap==1 ) THEN
          !if forces are not here in the first snapshot, ignore them for the rest
          IF( ALLOCATED(AUXNAMES) .AND. SIZE(AUXNAMES).NE.0 ) THEN
            fx=0
            fy=0
            fz=0
            DO j=1,SIZE(AUXNAMES)
              IF( AUXNAMES(j)=="fx" ) THEN
                fx=j
              ELSEIF( AUXNAMES(j)=="fy" ) THEN
                fy=j
              ELSEIF( AUXNAMES(j)=="fz" ) THEN
                fz=j
              ENDIF
            ENDDO
            !
            IF( fx.NE.0 .OR. fy.NE.0 .OR. fz.NE.0 ) THEN
              !The force is strong in AUX => save to AUX_FORCES
              xsf_forces = .TRUE.
              IF(ALLOCATED(AUX_FORCES)) DEALLOCATE(AUX_FORCES)
              ALLOCATE( AUX_FORCES( SIZE(P,1),3 ) )
              AUX_FORCES(:,:) = 0.d0
              IF(fx.NE.0) AUX_FORCES(:,1) = AUX(:,fx)
              IF(fy.NE.0) AUX_FORCES(:,2) = AUX(:,fy)
              IF(fz.NE.0) AUX_FORCES(:,3) = AUX(:,fz)
            ELSE
              !No force in the first file of the list
              !=> forces will not be written at all for any snapshot
              xsf_forces = .FALSE.
            ENDIF
            !
          ELSE
            xsf_forces = .FALSE.
          ENDIF
          !
        ELSE !i.e. this is not the first snapshot
          IF( xsf_forces ) THEN
            !Forces must be written to XSF file => make sure AUX_FORCES is properly set
            !
            !By default set all forces to zero
            IF(ALLOCATED(AUX_FORCES)) DEALLOCATE(AUX_FORCES)
            ALLOCATE( AUX_FORCES( SIZE(P,1),3 ) )
            AUX_FORCES(:,:) = 0.d0
            !
            IF(ALLOCATED(AUX)) THEN
              !some auxiliary properties were read from the present file
              !=> check if it contains forces
              fx=0
              fy=0
              fz=0
              DO j=1,SIZE(AUXNAMES)
                IF( AUXNAMES(j)=="fx" ) THEN
                  fx=j
                ELSEIF( AUXNAMES(j)=="fy" ) THEN
                  fy=j
                ELSEIF( AUXNAMES(j)=="fz" ) THEN
                  fz=j
                ENDIF
              ENDDO
              !
              IF( fx.NE.0 .OR. fy.NE.0 .OR. fz.NE.0 ) THEN
                !Some forces were found => save them in AUX_FORCES
                IF(fx.NE.0) AUX_FORCES(:,1) = AUX(:,fx)
                IF(fy.NE.0) AUX_FORCES(:,2) = AUX(:,fy)
                IF(fz.NE.0) AUX_FORCES(:,3) = AUX(:,fz)
              ENDIF
            ENDIF
            !
          ENDIF  !endif xsf_forces
          !
        ENDIF  !end snap
        !
        !Write current snapshot to animated XSF file
        WRITE(41,*) "CRYSTAL"
        WRITE(41,*) "PRIMVEC", snap
        WRITE(41,220) H(1,1), H(1,2), H(1,3)
        WRITE(41,220) H(2,1), H(2,2), H(2,3)
        WRITE(41,220) H(3,1), H(3,2), H(3,3)
        WRITE(41,*) "CONVVEC", snap
        WRITE(41,220) H(1,1), H(1,2), H(1,3)
        WRITE(41,220) H(2,1), H(2,2), H(2,3)
        WRITE(41,220) H(3,1), H(3,2), H(3,3)
        WRITE(41,*) "PRIMCOORD", snap
        WRITE(41,'(i9,a2)') SIZE(P,1), ' 1'
        220 FORMAT(3(f16.8,2X))
        !Write atom positions
        DO j=1,SIZE(P,1)
          IF( xsf_forces ) THEN
            WRITE(temp,221) INT(P(j,4)), P(j,1), P(j,2), P(j,3), (AUX_FORCES(j,k), k=1,3)
          ELSE
            WRITE(temp,222) INT(P(j,4)), P(j,1), P(j,2), P(j,3)
          ENDIF
          !Write line to the file
          WRITE(41,*) TRIM(ADJUSTL(temp))
        ENDDO
        221 FORMAT(i5,2X,6(f16.8,2X))
        222 FORMAT(i5,2X,3(f16.8,2X))
        WRITE(msg,*) 'Snapshot written to animated XSF file'
        CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      ENDIF
      !
      !-- add other formats here --
      !
      !
    ELSE
      !File does not exist => skip to next file
      nwarn = nwarn+1
      CALL ATOMSK_MSG(4700,(/TRIM(inputfile)/),(/0.d0/))
    ENDIF
    !
    CALL ATOMSK_MSG(4001,(/''/),(/0.d0/))
    !
  ENDIF  !endif #
ENDDO
!
!
!
300 CONTINUE
CLOSE(35)
CALL ATOMSK_MSG(4042,(/''/),(/DBLE(snap)/))
IF(output_xyz .OR. output_exyz) THEN
  CLOSE(40)
  msg = "XYZ"
  CALL ATOMSK_MSG(3002,(/msg,outfile_xyz/),(/0.d0/))
ENDIF
IF(output_xsf) THEN
  CLOSE(41)
  msg = "XSF"
  CALL ATOMSK_MSG(3002,(/msg,outfile_xsf/),(/0.d0/))
ENDIF
!
!
!
1000 CONTINUE
!
!
!
END SUBROUTINE ALLINONE
!
END MODULE