MODULE oia_lmc
!
!**********************************************************************************
!*  1IA_LMC                                                                       *
!**********************************************************************************
!* This module reads files in the LAMMPS custom dump format, that contains        *
!* several snapshots.                                                             *
!* The LAMMPS custom dump format is described here:                               *
!*    http://lammps.sandia.gov/doc/dump.html                                      *
!**********************************************************************************
!* (C) January 2019 - Pierre Hirel                                                *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 28 Jan. 2019                                     *
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
USE writeout
!
!
CONTAINS
!
SUBROUTINE ONEINALL_LMC(inputfile,out_prefix,outfileformats,options_array)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=128):: msg
CHARACTER(LEN=4096):: out_prefix, outputfile, outfileformat
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to output
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMES !names of auxiliary properties
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment
CHARACTER(LEN=128),DIMENSION(30):: acol  !all columns for each atom (no more than 30 columns)
CHARACTER(LEN=128),DIMENSION(30):: namecol  !names of columns (assuming no more than 30 columns)
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: options_array !options and their parameters
CHARACTER(LEN=1024):: temp, temp2, temp3
LOGICAL:: finished
LOGICAL,DIMENSION(3):: reduced  !are the coordinates reduced?
LOGICAL,DIMENSION(:),ALLOCATABLE:: SELECT  !mask for atom list
INTEGER:: i, id, itemp2, j, NP
INTEGER:: Ncol   !number of columns for each atom
INTEGER:: Nsnap  !counter for snapshots
INTEGER:: sizeaux  !size of array for auxiliary properties
INTEGER:: typecol, xcol, ycol, zcol, idcol, elecol !column for species, X, Y, Z, id, element
INTEGER:: vx, vy, vz, fx, fy, fz, q !columns for velocities, forces
REAL(dp):: a, b, c, alpha, beta, gamma
REAL(dp):: testreal
REAL(dp):: lx, ly, lz
REAL(dp):: xlo_bound, xhi_bound, ylo_bound, yhi_bound, zlo_bound, zhi_bound
REAL(dp):: xy, xz, yz  !for triclinic boxes
REAL(dp),DIMENSION(3,3):: H   !Base vectors of the supercell
REAL(dp),DIMENSION(3,3):: ORIENT  !crystal orientation
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: P
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: S
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUX !auxiliary properties
!
!
!Initialize variables
outputfile=''
Nsnap = 0
H(:,:) = 0.d0
ORIENT(:,:) = 0.d0
!
!
100 CONTINUE
msg = 'debug --> entering READ_LMP_CUSTOM'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
OPEN(UNIT=30,FILE=inputfile,FORM='FORMATTED',STATUS='OLD')
REWIND(30)
!
!
140 CONTINUE
DO
  !Parse the file and search for keywords
  READ(30,'(a1024)',ERR=1000,END=1000) temp
  temp = TRIM(ADJUSTL(temp))
  IF(temp(1:5)=='ITEM:') THEN
    temp2 = TRIM(ADJUSTL(temp(6:)))
    IF(temp2(1:8)=='TIMESTEP') THEN
      !This is a new snapshot
      Nsnap = Nsnap+1
      !Initialize variables
      !NOTE: box vectors H are not initialized at each snapshot
      !     so that, by default, vectors of the first snapshot are used
      acol(:) = ""
      namecol(:) = ""
      i = 0
      Ncol = 0
      elecol = 0
      idcol = 0
      typecol = 0
      xcol = 0
      ycol = 0
      zcol = 0
      xlo_bound = 0.d0
      ylo_bound = 0.d0
      zlo_bound = 0.d0
      vx = 0
      vy = 0
      vz = 0
      fx = 0
      fy = 0
      fz = 0
      q = 0
      finished = .FALSE.
      reduced(:) = .FALSE.
      IF(ALLOCATED(P)) DEALLOCATE(P)
      IF(ALLOCATED(comment)) DEALLOCATE(comment)
      IF(ALLOCATED(AUX)) DEALLOCATE(AUX)
      IF(ALLOCATED(AUXNAMES)) DEALLOCATE(AUXNAMES)
      !Read timestep and use it to generate a comment
      READ(30,*,END=820,ERR=820) id
      WRITE(msg,*) id
      ALLOCATE(comment(1))
      comment(1) = "# LAMMPS timestep "//TRIM(ADJUSTL(msg))
      !Set the output file name
      WRITE(msg,*) Nsnap
      outputfile = TRIM(ADJUSTL(out_prefix))//TRIM(ADJUSTL(msg))
      WRITE(msg,*) 'snap, outputfile:', Nsnap, TRIM(ADJUSTL(outputfile))
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      !
      CALL ATOMSK_MSG(4041,(/''/),(/DBLE(Nsnap)/))
      !
    ELSEIF(temp2(1:15)=='NUMBER OF ATOMS') THEN
      READ(30,*,END=820,ERR=820) NP
      ALLOCATE(P(NP,4))
      P(:,:) = 0.d0
    ELSEIF(temp2(1:10)=='BOX BOUNDS') THEN
      H(:,:) = 0.d0
      IF( temp2(12:19)=='xy xz yz' ) THEN
        !If the box is triclinic we read the bounding box and tilt factors
        READ(30,*,END=810,ERR=810) xlo_bound, xhi_bound, xy
        READ(30,*,END=810,ERR=810) ylo_bound, yhi_bound, xz
        READ(30,*,END=810,ERR=810) zlo_bound, zhi_bound, yz
      ELSE
        !If the box is not triclinic we set the tilt factors to zero
        READ(30,*,END=810,ERR=810) xlo_bound, xhi_bound
        READ(30,*,END=810,ERR=810) ylo_bound, yhi_bound
        READ(30,*,END=810,ERR=810) zlo_bound, zhi_bound
        xy=0.d0
        xz=0.d0
        yz=0.d0
      ENDIF
      !Transform LAMMPS parameters and tilt factors into a proper H matrix
      lx = xhi_bound-MAX(0.d0,xy,xz,xy+xz) - xlo_bound+MIN(0.d0,xy,xz,xy+xz)
      ly = yhi_bound-MAX(0.d0,yz) - ylo_bound+MIN(0.d0,yz)
      lz = zhi_bound-zlo_bound
      a = lx
      b = DSQRT( ly**2 + xy**2 )
      c = DSQRT( lz**2 + xz**2 + yz**2 )
      alpha = DACOS( (xy*xz+ly*yz)/(b*c) )
      beta = DACOS(xz/c)
      gamma = DACOS(xy/b)
      CALL CONVMAT(a,b,c,alpha,beta,gamma,H)
    ELSEIF(temp2(1:5)=='ATOMS') THEN
      !At this point if P is not allocated we are in trouble
      IF(.NOT.ALLOCATED(P)) GOTO 820
      !
      itemp2 = 5
      IF( LEN_TRIM(temp2(itemp2+1:))<=0 ) THEN
        !No column specified: use defaults
        Ncol = 5
        idcol = 1
        typecol = 2
        xcol = 3
        ycol = 4
        zcol = 5
        reduced(:) = .TRUE.
        !
      ELSE
        !Columns are specified, interpret them
        DO i=1,SIZE(acol)
          temp2 = TRIM(ADJUSTL(temp2(itemp2+1:)))
          READ(temp2,*,END=150,ERR=150) temp3
          !No error? Then read what this column contains
          temp3 = TRIM(ADJUSTL(temp3))
          IF(temp3=='type') THEN
            typecol = i
          ELSEIF(temp3(1:1)=='x' .AND. xcol==0) THEN
            xcol = i
            IF(temp3(2:2)=='s') reduced(1) = .TRUE.
          ELSEIF(temp3(1:1)=='y' .AND. ycol==0) THEN
            ycol = i
            IF(temp3(2:2)=='s') reduced(2) = .TRUE.
          ELSEIF(temp3(1:1)=='z' .AND. zcol==0) THEN
            zcol = i
            IF(temp3(2:2)=='s') reduced(3) = .TRUE.
          ELSEIF(temp3(1:2)=='vx' .AND. vx==0) THEN
            vx = i
          ELSEIF(temp3(1:2)=='vy' .AND. vy==0) THEN
            vy = i
          ELSEIF(temp3(1:2)=='vz' .AND. vz==0) THEN
            vz = i
          ELSEIF(temp3(1:2)=='fx' .AND. fx==0) THEN
            fx = i
          ELSEIF(temp3(1:2)=='fy' .AND. fy==0) THEN
            fy = i
          ELSEIF(temp3(1:2)=='fz' .AND. fz==0) THEN
            fz = i
          ELSEIF(temp3(1:2)=='id' .AND. idcol==0) THEN
            idcol = i
          ELSEIF(temp3(1:2)=='q ' .AND. q==0) THEN
            q = i
          ELSEIF(temp3(1:7)=='element' .AND. elecol==0) THEN
            elecol = i
          ELSE
            namecol(i) = TRIM(temp3)
          ENDIF
          !Increase the number of columns that the file actually contains
          Ncol = Ncol+1
          itemp2 = LEN_TRIM(temp3)
        ENDDO
        !
      ENDIF
      !
      150 CONTINUE
      IF(verbosity==4) THEN
        WRITE(msg,*) 'debug --> Ncol = ', Ncol
        CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
        WRITE(msg,*) 'debug --> Col type, element, X, Y, Z: ', typecol, elecol, xcol, ycol, zcol
        CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
        WRITE(msg,*) 'debug --> Col Vxyz, Fxyz, q: ', vx, vy, vz, fx, fy, fz, q
        CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      ENDIF
      !If some columns were not detected we are in trouble
      IF(xcol==0 .AND. ycol==0 .AND. zcol==0 ) THEN
        GOTO 810
      ENDIF
      !Allocate auxiliary array if auxiliary properties were detected
      sizeaux = 0
      IF( typecol>0 ) sizeaux = sizeaux+1
      IF( vx>0 ) sizeaux = sizeaux+1
      IF( vy>0 ) sizeaux = sizeaux+1
      IF( vz>0 ) sizeaux = sizeaux+1
      IF( fx>0 ) sizeaux = sizeaux+1
      IF( fy>0 ) sizeaux = sizeaux+1
      IF( fz>0 ) sizeaux = sizeaux+1
      IF( q>0 ) sizeaux = sizeaux+1
      !If other properties are present, their names are in namecol(:)
      !=> count them
      DO i=1,Ncol
        IF( LEN_TRIM(namecol(i))>0 ) sizeaux = sizeaux+1
      ENDDO
      !
      IF(sizeaux>0) THEN
        ALLOCATE( AUXNAMES(sizeaux) )
        ALLOCATE( AUX( SIZE(P,1), sizeaux ) )
        AUX(:,:) = 0.d0
        !Set the names of auxiliary properties
        sizeaux=0
        IF( typecol.NE.0 ) THEN
          sizeaux = sizeaux+1
          AUXNAMES(sizeaux) = 'type'
        ENDIF
        IF( vx>0 ) THEN
          sizeaux = sizeaux+1
          AUXNAMES(sizeaux) = 'vx'
        ENDIF
        IF( vy>0 ) THEN
          sizeaux = sizeaux+1
          AUXNAMES(sizeaux) = 'vy'
        ENDIF
        IF( vz>0 ) THEN
          sizeaux = sizeaux+1
          AUXNAMES(sizeaux) = 'vz'
        ENDIF
        IF( fx>0 ) THEN
          sizeaux = sizeaux+1
          AUXNAMES(sizeaux) = 'fx'
        ENDIF
        IF( fy>0 ) THEN
          sizeaux = sizeaux+1
          AUXNAMES(sizeaux) = 'fy'
        ENDIF
        IF( fz>0 ) THEN
          sizeaux = sizeaux+1
          AUXNAMES(sizeaux) = 'fz'
        ENDIF
        IF( q>0 ) THEN
          sizeaux = sizeaux+1
          AUXNAMES(sizeaux) = 'q'
        ENDIF
        DO i=1,Ncol
          IF( LEN_TRIM(namecol(i))>0 ) THEN
            sizeaux = sizeaux+1
            AUXNAMES(sizeaux) = namecol(i)
          ENDIF
        ENDDO
      ENDIF
      WRITE(msg,*) 'debug --> sizeaux, SIZE(AUX) = ', sizeaux, SIZE(AUX,2)
      CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
      !
      !
      !Read each line and store positions and species to P
      DO i=1,SIZE(P,1)
        READ(30,'(a1024)',END=840,ERR=840) temp3
        READ(temp3,*,END=840,ERR=840) (acol(j), j=1,Ncol)
        IF(idcol>0) THEN
          READ(acol(idcol),*,END=840,ERR=840) id
        ELSE
          id = i
        ENDIF
        IF( id<=0 .OR. id>SIZE(P,1) ) THEN
          nwarn=nwarn+1
          CALL ATOMSK_MSG(2742,(/""/),(/DBLE(id)/))
        ELSE
          IF(xcol>0) READ(acol(xcol),*,END=840,ERR=840) a
          IF(ycol>0) READ(acol(ycol),*,END=840,ERR=840) b
          IF(zcol>0) READ(acol(zcol),*,END=840,ERR=840) c
          !If reduced, convert to Cartesian coordinates
          !Also, the box was shifted so that the lower bounds are at (0,0,0)
          !=> shift atoms to keep their relative positions to the box
          IF( reduced(1) ) THEN
            P(id,1) = a*H(1,1) + b*H(2,1) + c*H(3,1) - xlo_bound
          ELSE
            P(id,1) = a - xlo_bound
          ENDIF
          IF( reduced(2) ) THEN
            P(id,2) = a*H(1,2) + b*H(2,2) + c*H(3,2) - ylo_bound
          ELSE
            P(id,2) = b - ylo_bound
          ENDIF
          IF( reduced(3) ) THEN
            P(id,3) = a*H(1,3) + b*H(2,3) + c*H(3,3) - zlo_bound
          ELSE
            P(id,3) = c - zlo_bound
          ENDIF
          !
          !Read element if present, otherwise type = element
          IF( elecol>0 ) THEN
            !Could be actual element name or a number
            !Try to pass the string, see if it is an element name
            CALL ATOMNUMBER(acol(elecol),testreal)
            IF( DABS(testreal)<1.d-12 ) THEN
              !String was not an element name => try a real
              READ( acol(elecol),*,END=170,ERR=170 ) testreal
              !It succeeded => store it in P(id,4)
              GOTO 180
              170 CONTINUE
              !It failed, it was not a real number => ignore it, use atom type
              READ(acol(typecol),*,END=840,ERR=840) testreal
            ENDIF
            !
          ELSEIF( typecol>0 ) THEN
            !No 'element' column => use atom type as atomic number
            READ(acol(typecol),*,END=840,ERR=840) testreal
            !
          ELSE
            !No atom type or element in dump file
            !=> we cannot guess elements, all atoms will be hydrogen
            testreal = 1.d0
          ENDIF
          180 CONTINUE
          P(id,4) = testreal
          !
          !Read auxiliary properties (if any)
          IF(ALLOCATED(AUX)) THEN
            sizeaux = 0
            IF( typecol>0 ) THEN
              sizeaux = 1
              READ(acol(typecol),*,END=840,ERR=840) a
              AUX(id,sizeaux) = a
            ENDIF
            IF( vx>0 ) THEN
              sizeaux = sizeaux+1
              READ(acol(vx),*,END=840,ERR=840) a
              AUX(id,sizeaux) = a
            ENDIF
            IF( vy>0 ) THEN
              sizeaux = sizeaux+1
              READ(acol(vy),*,END=840,ERR=840) a
              AUX(id,sizeaux) = a
            ENDIF
            IF( vz>0 ) THEN
              sizeaux = sizeaux+1
              READ(acol(vz),*,END=840,ERR=840) a
              AUX(id,sizeaux) = a
            ENDIF
            IF( fx>0 ) THEN
              sizeaux = sizeaux+1
              READ(acol(fx),*,END=840,ERR=840) a
              AUX(id,sizeaux+1) = a
            ENDIF
            IF( fy>0 ) THEN
              sizeaux = sizeaux+1
              READ(acol(fy),*,END=840,ERR=840) a
              AUX(id,sizeaux+1) = a
            ENDIF
            IF( fz>0 ) THEN
              sizeaux = sizeaux+1
              READ(acol(fz),*,END=840,ERR=840) a
              AUX(id,sizeaux+1) = a
            ENDIF
            IF( q>0 ) THEN
              READ(acol(q),*,END=840,ERR=840) a
              AUX(id,sizeaux+1) = a
              sizeaux = sizeaux+1
            ENDIF
            DO j=1,Ncol
              IF( LEN_TRIM(namecol(j))>0 ) THEN
                READ(acol(j),*,END=840,ERR=840) a
                sizeaux = sizeaux+1
                IF( sizeaux <= SIZE(AUX,2) ) THEN
                  AUX(id,sizeaux) = a
                ENDIF
              ENDIF
            ENDDO
          ENDIF
          !
        ENDIF
        !
      ENDDO !end loop on atom positions
      finished = .TRUE.
      !
    ENDIF
  ENDIF
  !
  IF( finished ) THEN
    !Apply options
    CALL OPTIONS_AFF(options_array,H,H,P,S,AUXNAMES,AUX,ORIENT,SELECT)
    !
    !Output snapshot to one or several format(s)
    outfileformat=''
    CALL WRITE_AFF(outputfile,outfileformats,H,P,S,comment,AUXNAMES,AUX)
    !
    CALL ATOMSK_MSG(4001,(/''/),(/0.d0/))
    finished = .FALSE.
  ENDIF
  !
ENDDO !end loop on snapshots
!
!
GOTO 1000
!
!
!
800 CONTINUE
CALL ATOMSK_MSG(802,(/''/),(/DBLE(i)/))
nerr = nerr+1
GOTO 1000
!
810 CONTINUE
CALL ATOMSK_MSG(1809,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
820 CONTINUE
CALL ATOMSK_MSG(1810,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
830 CONTINUE
CALL ATOMSK_MSG(1804,(/''/),(/0.d0/))
nerr = nerr+1
GOTO 1000
!
840 CONTINUE
CALL ATOMSK_MSG(1812,(/''/),(/DBLE(i)/))
nerr = nerr+1
GOTO 1000
!
!
!
1000 CONTINUE
CLOSE(30)
CALL ATOMSK_MSG(4045,(/''/),(/DBLE(Nsnap)/))
!
!
END SUBROUTINE ONEINALL_LMC
!
END MODULE oia_lmc
