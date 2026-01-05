MODULE polyx_readparam
!
!**********************************************************************************
!*  POLYX_READ_PARAM                                                              *
!**********************************************************************************
!* This module reads parameters from a file concerning the construction of        *
!* a polycrystal. It is called by subroutine "POLYCRYS"                           *
!* (see modes/mode_polycrystal.f90)                                               *
!**********************************************************************************
!* (C) June 2025 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 05 Jan. 2026                                     *
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
USE math
USE crystallography
USE messages
USE files
USE files_msg
USE random
USE subroutines
!
CONTAINS
!
!
SUBROUTINE POLYX_READ_PARAM(vfile,Huc,Puc,H,Nnodes,vnodes,vorient,twodim,clearance,outparam,status)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: vfile   !name of file containing parameters for Voronoi construction
CHARACTER(LEN=128):: lattice  !if grains are organized according to a lattice
CHARACTER(LEN=128):: or1, or2, or3
CHARACTER(LEN=4096):: line
CHARACTER(LEN=4096):: msg, temp
LOGICAL:: Hset            !are the box vectors H(:,:) defined?
LOGICAL:: miller          !are Miller indices given? (if no then angles are given)
LOGICAL,DIMENSION(4),INTENT(OUT):: outparam  !Are keywords (1) "node", (2) "lattice", (3) "random" used?
                                             !(4) Must parameters be saved in a text file?
INTEGER:: i, j, m, n, o
INTEGER:: linenumber    !line number
INTEGER:: status        !=0 if successful, >0 otherwise
INTEGER:: twodim        !=0 if system is 3-D, =1,2,3 if system is thin along x, y, z
INTEGER,INTENT(INOUT):: Nnodes      !number of nodes
REAL(dp):: distance    !distance between two points
REAL(dp):: grad1, grad2!gradient values
REAL(dp):: P1, P2, P3  !temporary position
REAL(dp):: Volume, Vmin !min. volume occupied by a grain
REAL(dp),INTENT(INOUT):: clearance   !clear atoms that close to the GB
REAL(dp),DIMENSION(3):: vector    !vector between an atom and a node
REAL(dp),DIMENSION(3,3),INTENT(IN):: Huc        !Base vectors of the seed
REAL(dp),DIMENSION(3,3),INTENT(INOUT):: H       !Base vectors of the final supercell
REAL(dp),DIMENSION(3,3):: rotmat  !rotation matrix
REAL(dp),DIMENSION(:),ALLOCATABLE:: randarray   !random numbers
REAL(dp),DIMENSION(:,:),INTENT(IN):: Puc       !positions of atoms in the seed
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: vnodes    !cartesian coordinate of each node
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(OUT):: vorient !rotation matrix for each node
!
!Initialize variables
Hset = .FALSE.
miller = .FALSE.
outparam(:) = .FALSE.
linenumber = 0
status = 0
twodim = 0
H(:,:) = 0.d0
IF(ALLOCATED(vnodes)) DEALLOCATE(vnodes)
IF(ALLOCATED(vorient)) DEALLOCATE(vorient)
!
!
msg = 'ENTERING POLYX_READ_PARAM...'
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!Read the file containing parameters for Voronoi construction
CALL CHECKFILE(vfile,'read')
OPEN(UNIT=31,FILE=vfile)
!Parse the file a first time to count number of nodes
Nnodes = 0
Hset=.FALSE.
linenumber = 0
DO
  READ(31,'(a)',END=210,ERR=210) line
  line = TRIM(ADJUSTL(line))
  linenumber = linenumber+1
  !
  !Ignore empty lines and lines starting with #
  IF( line(1:1).NE."#" .AND. LEN_TRIM(line)>0 ) THEN
    IF( StrDnCase(line(1:3))=="box" ) THEN
      !Read size of the final box
      READ(line(4:),*,END=830,ERR=830) P1, P2, P3
      !Set final box vectors
      H(:,:) = 0.d0
      H(1,1) = DABS(P1)
      H(2,2) = DABS(P2)
      H(3,3) = DABS(P3)
      WRITE(msg,'(a11,3(f9.3,a3))') " Read box: ", H(1,1), " x ", H(2,2), " x ",  H(3,3)
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      !
      IF( H(1,1)<1.1d0*Huc(1,1) .OR. H(2,2)<1.1d0*Huc(2,2) .OR. H(3,3)<1.1d0*Huc(3,3) ) THEN
        !The user asked for a final cell with at least one small dimension (smaller than seed)
        !=> Look along which dimension the final box is "small"
        twodim=0  !=1 if cell is small along X; =2 along Y; =3 along Z
        m=0  !counter for how many dimensions are small (only one is allowed to be small)
        DO i=1,3
          IF( i==1 ) msg="X"
          IF( i==2 ) msg="Y"
          IF( i==3 ) msg="Z"
          IF( VECLENGTH(H(i,:)) < 1.1d0*VECLENGTH(Huc(i,:)) ) THEN
            !The final box is zero along this dimension
            m=m+1
            IF( m>=2 ) THEN
              !Final box is small in many directions => error
              nerr=nerr+1
              CALL ATOMSK_MSG(4828,(/""/),(/0.d0/))
              GOTO 1000
            ELSE
              !The final box is small in only one dimension => pseudo-2D system
              IF( VECLENGTH(H(i,:)) < 1.d-12 ) THEN
                !Cell size was actually zero along that dimension => warn that it will be resized
                nwarn=nwarn+1
                CALL ATOMSK_MSG(4714,(/TRIM(msg)/),(/Huc(i,i)/))
              ENDIF
              twodim = i
              !Make sure that the final box dimension matches the seed dimension in that direction
              H(i,i) = Huc(i,i)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
      WRITE(msg,*) "            twodim = ", twodim
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      Hset=.TRUE.
      !
      Volume = H(1,1) * H(2,2) * H(3,3)
      IF( twodim==0 .AND. Volume < 8000.d0 ) THEN
        !Final cell is very small: display warning
        nwarn = nwarn+1
        CALL ATOMSK_MSG(4719,(/""/),(/Volume/))
      ENDIF
      !
      !Check that there will be enough memory for final system
      Vmin = VOLUME_PARA(Huc)  ! Volume of seed
      Volume = VOLUME_PARA(H)  ! Volume of final box
      P1 = SIZE(Puc,1) * Volume/Vmin  !estimate of number of atoms in final box
      !Check if number of atoms (P1) is ok
      CALL CHECKMEM(P1,i)
      IF( i>0 ) THEN
        ! N.atoms too large, unable to allocate
        nerr = nerr+1
        CALL ATOMSK_MSG(819,(/''/),(/P1/))
        GOTO 1000
      ENDIF
      !
    ELSEIF( StrDnCase(line(1:5))=="node " .OR. StrDnCase(line(1:5))=="grain" ) THEN
      !Check that the box was defined
      IF( .NOT.Hset ) THEN
        GOTO 820
      ENDIF
      Nnodes = Nnodes+1
      outparam(1) = .TRUE.
      !
    ELSEIF( StrDnCase(line(1:7))=="lattice" ) THEN
      !Check that the box was defined
      IF( .NOT.Hset ) THEN
        GOTO 820
      ENDIF
      outparam(2) = .TRUE.
      lattice = TRIM(ADJUSTL(line(8:)))
      !Set the number of nodes according to the lattice type
      !Beware of pseudo-2D systems
      IF( StrDnCase(lattice)=="sc" ) THEN
        Nnodes = 1
      ELSEIF( StrDnCase(lattice)=="bcc" ) THEN
        Nnodes = 2
      ELSEIF( StrDnCase(lattice)=="fcc" ) THEN
        IF(twodim>0) THEN
          Nnodes = 3
        ELSE
          Nnodes = 4
        ENDIF
      ELSEIF( StrDnCase(lattice)=="diamond" .OR. StrDnCase(lattice)=="dia"  ) THEN
        IF(twodim>0) THEN
          Nnodes = 6
        ELSE
          Nnodes = 8
        ENDIF
      ELSEIF( StrDnCase(lattice)=="hcp" ) THEN
        IF(twodim>0) THEN
          Nnodes = 2
        ELSE
          Nnodes = 4
        ENDIF
      ELSE
        !Unrecognized lattice type => abort
        nerr=nerr+1
        GOTO 1000
      ENDIF
      !
    ELSEIF( StrDnCase(line(1:6))=="random" ) THEN
      !Check that the box was defined
      IF( .NOT.Hset ) THEN
        GOTO 820
      ENDIF
      !Read total number of grains
      READ(line(7:),*,END=800,ERR=800) Nnodes
      outparam(3) = .TRUE.
      WRITE(msg,*) " Read random ", Nnodes
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      !
    ELSEIF( StrDnCase(line(1:9))=="clearance" .OR. StrDnCase(line(1:5))=="clear" ) THEN
      !User wants to change default clearance around GB plane
      READ(line,*,END=800,ERR=800) temp, P1
      !Make sure clearance is a positive number (or zero)
      clearance = DABS(P1)
      WRITE(msg,*) " Read clearance ", clearance
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      !
    ELSEIF( StrDnCase(line(1:8))=="gradient" ) THEN
      !User wants to create a gradient along a direction between P1 and P2
      READ(line(9:),*,END=800,ERR=800) or1
      WRITE(msg,*) " Read gradient ", or1
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      !
!     ELSEIF( line(1:7)=="protect" ) THEN
!       !User asked to preserve unit cell
!       protectuc = .TRUE.
      !
    ELSE
      !Unknown command => display warning
      nwarn=nwarn+1
      temp(1:64) = CHARLONG2SHRT(vfile)
      CALL ATOMSK_MSG(1702,(/line,temp/),(/0.d0/))
    ENDIF
    !
  ENDIF
ENDDO
!
210 CONTINUE
!Keywords "lattice", "node" and "random" are mutually exclusive
!If two of them appear in the param file, display an error message and exit
IF( outparam(1) .AND. outparam(2) ) THEN
  CALL ATOMSK_MSG(4823,(/"node   ","random "/),(/0.d0/))
  nerr=nerr+1
  GOTO 1000
ELSEIF( outparam(1) .AND. outparam(3) ) THEN
  CALL ATOMSK_MSG(4823,(/"node   ","lattice"/),(/0.d0/))
  nerr=nerr+1
  GOTO 1000
ELSEIF( outparam(2) .AND. outparam(3) ) THEN
  CALL ATOMSK_MSG(4823,(/"random ","lattice"/),(/0.d0/))
  nerr=nerr+1
  GOTO 1000
ENDIF
!
!Final positions of nodes will be stored in array vnodes(:,:)
!Final rotation matrix for each grain will be stored in array vorient(:,:)
ALLOCATE(vnodes(Nnodes,3))
vnodes(:,:) = 0.d0
ALLOCATE(vorient(Nnodes,3,3))
vorient(:,:,:) = 0.d0
!
REWIND(31)
linenumber = 0
Nnodes = 0
DO
  READ(31,'(a)',END=250,ERR=250) line
  line = TRIM(ADJUSTL(line))
  linenumber = linenumber+1
  !Ignore empty lines and lines starting with #
  IF( line(1:1).NE."#" .AND. LEN_TRIM(line)>0 ) THEN
    !
    IF( StrDnCase(line(1:3))=="box" ) THEN
      !This part was already dealt with before
      !
      !
    ELSEIF( StrDnCase(line(1:7))=="lattice" ) THEN
      !Nodes will be placed according to a pattern
      !Patterns are defined by "lattice", can be fcc, bcc, etc.
      !The crystallographic orientation of the grains will be random
      !
      !Check that the box was defined
      IF( .NOT.Hset ) THEN
        GOTO 820
      ENDIF
      !
      lattice = TRIM(ADJUSTL(line(8:)))
      IF( StrDnCase(lattice)=="sc" ) THEN
        Nnodes = 1
        vnodes(1,:) = 0.d0
      ELSEIF( StrDnCase(lattice)=="bcc" ) THEN
        Nnodes = 2
        vnodes(1,:) = 0.d0          !(0,0,0)
        vnodes(2,1) = 0.5d0*H(1,1)  !(1/2,1/2,1/2)
        vnodes(2,2) = 0.5d0*H(2,2)
        vnodes(2,3) = 0.5d0*H(3,3)
        IF(twodim>0) THEN
          vnodes(2,twodim) = 0.d0
        ENDIF
      ELSEIF( StrDnCase(lattice)=="fcc" ) THEN
        IF(twodim>0) THEN
          !System is 2-D => define only 3 nodes
          Nnodes = 3
          vnodes(1,:) = 0.d0          !(0,0)
          IF(twodim==1) THEN
            vnodes(2,2) = 0.d0*H(2,2)  !(0,0)
            vnodes(2,3) = 0.5d0*H(3,3)
            vnodes(3,2) = 0.5d0*H(2,2) !(1/2,0)
            vnodes(3,3) = 0.d0*H(3,3)
          ELSEIF(twodim==2) THEN
            vnodes(2,1) = 0.d0*H(1,1)  !(0,0)
            vnodes(2,3) = 0.5d0*H(3,3)
            vnodes(3,1) = 0.5d0*H(1,1) !(1/2,0)
            vnodes(3,3) = 0.d0*H(3,3)
          ELSEIF(twodim==3) THEN
            vnodes(2,1) = 0.d0*H(1,1)  !(0,0)
            vnodes(2,2) = 0.5d0*H(2,2)
            vnodes(3,1) = 0.5d0*H(1,1) !(1/2,0)
            vnodes(3,2) = 0.d0*H(2,2)
          ENDIF
          vnodes(:,twodim) = 0.5d0
        ELSE
          !System is 3-D => define fcc lattice
          Nnodes = 4
          vnodes(:,:) = 0.d0          !(0,0,0)
          vnodes(2,1) = 0.5d0*H(1,1)  !(1/2,1/2,0)
          vnodes(2,2) = 0.5d0*H(2,2)
          vnodes(3,1) = 0.5d0*H(1,1)  !(1/2,0,1/2)
          vnodes(3,3) = 0.5d0*H(3,3)
          vnodes(4,2) = 0.5d0*H(2,2)  !(0,1/2,1/2)
          vnodes(4,3) = 0.5d0*H(3,3)
        ENDIF
      ELSEIF( StrDnCase(lattice)=="diamond" .OR. StrDnCase(lattice)=="dia" ) THEN
        IF(twodim>0) THEN
          !System is 2-D => define only 6 nodes
          Nnodes = 6
          vnodes(:,:) = 0.d0
          vnodes(2,:) = 0.5d0*H(1,1)  !(1/2,1/2)
          vnodes(2,2) = 0.5d0*H(2,2)
          vnodes(2,3) = 0.5d0*H(3,3)
          vnodes(3,1) = 0.25d0*H(1,1) !(1/4,1/4)
          vnodes(3,2) = 0.25d0*H(2,2)
          vnodes(3,3) = 0.25d0*H(3,3)
          vnodes(4,1) = 0.75d0*H(1,1) !(3/4,3/4)
          vnodes(4,2) = 0.75d0*H(2,2)
          vnodes(4,3) = 0.75d0*H(3,3)
          IF(twodim==1) THEN
            vnodes(5,2) = 0.25d0*H(2,2)  !(1/4,3/4)
            vnodes(5,3) = 0.75d0*H(3,3)
            vnodes(6,2) = 0.75d0*H(2,2)  !(3/4,1/4)
            vnodes(6,3) = 0.25d0*H(3,3)
          ELSEIF(twodim==2) THEN
            vnodes(5,1) = 0.25d0*H(1,1)  !(1/4,3/4)
            vnodes(5,3) = 0.75d0*H(3,3)
            vnodes(6,1) = 0.75d0*H(1,1)  !(3/4,1/4)
            vnodes(6,3) = 0.25d0*H(3,3)
          ELSE   !i.e. twodim==3
            vnodes(5,1) = 0.25d0*H(1,1)  !(1/4,3/4)
            vnodes(5,2) = 0.75d0*H(2,2)
            vnodes(6,1) = 0.75d0*H(1,1)  !(3/4,1/4)
            vnodes(6,2) = 0.25d0*H(2,2)
          ENDIF
          vnodes(:,twodim) = 0.5d0
        ELSE
          !System is 3-D => define diamond lattice
          Nnodes = 8
          vnodes(:,:) = 0.d0
          vnodes(2,1) = 0.5d0*H(1,1)  !(1/2,1/2,0)
          vnodes(2,2) = 0.5d0*H(2,2)
          vnodes(3,1) = 0.5d0*H(1,1)  !(1/2,0,1/2)
          vnodes(3,3) = 0.5d0*H(3,3)
          vnodes(4,2) = 0.5d0*H(2,2)  !(0,1/2,1/2)
          vnodes(4,3) = 0.5d0*H(3,3)
          vnodes(5,1) = 0.25d0*H(1,1) !(1/4,1/4,1/4)
          vnodes(5,2) = 0.25d0*H(2,2)
          vnodes(5,3) = 0.25d0*H(3,3)
          vnodes(6,1) = 0.75d0*H(1,1) !(3/4,3/4,1/4)
          vnodes(6,2) = 0.75d0*H(2,2)
          vnodes(6,3) = 0.25d0*H(3,3)
          vnodes(7,1) = 0.75d0*H(1,1) !(3/4,1/4,3/4)
          vnodes(7,2) = 0.25d0*H(2,2)
          vnodes(7,3) = 0.75d0*H(3,3)
          vnodes(8,1) = 0.25d0*H(1,1) !(1/4,3/4,3/4)
          vnodes(8,2) = 0.75d0*H(2,2)
          vnodes(8,3) = 0.75d0*H(3,3)
        ENDIF
      ELSEIF( StrDnCase(lattice)=="hcp" ) THEN
        IF(twodim>0) THEN
          !System is 2-D => define only 2 nodes
          Nnodes = 2
          vnodes(:,:) = 0.d0          !(0,0,0)
          vnodes(2,1) = 0.5d0*H(1,1)  !(1/2,1/2,1/2)
          vnodes(2,2) = 0.5d0*H(2,2)
          vnodes(2,3) = 0.5d0*H(3,3)
          vnodes(2,twodim) = 0.d0
        ELSE
          !System is 3-D => define hcp lattice
          Nnodes = 4
          vnodes(:,:) = 0.d0          !(0,0,0)
          vnodes(2,2) = H(2,2)/3.d0   !(0,1/3,1/2)
          vnodes(2,3) = 0.5d0*H(3,3)
          vnodes(3,1) = 0.5d0*H(1,1)  !(1/2,1/2,0)
          vnodes(3,2) = 0.5d0*H(2,2)
          vnodes(4,1) = 0.5d0*H(1,1)  !(1/2,5/6,1/2)
          vnodes(4,2) = 5.d0*H(2,2)/6.d0
          vnodes(4,3) = 0.5d0*H(3,3)
        ENDIF
      ELSE
        !unrecognized lattice type (already dealt with before)
      ENDIF
      !
      !Orientation of each grain will be random
      !Generate a list of 3*Nnodes random numbers
      CALL GEN_NRANDNUMBERS( 3*Nnodes , randarray )
      !randarray now contains 3*Nnodes real numbers between 0 and 1
      !They are used to generate rotation matrices
      IF( twodim>0 ) THEN
        !Only one random number will be used => Multiply all of them by 2*pi and subtract pi
        randarray(:) = randarray(:)*2.d0*pi - pi
      ELSE
        DO i=1,Nnodes
          m = 3*(i-1) + 1
          n = 3*(i-1) + 2
          o = 3*(i-1) + 3
          randarray(m) = randarray(m)*2.d0*pi - pi
          randarray(n) = DACOS(2.d0*randarray(n) - 1.d0)
          randarray(o) = randarray(o)*2.d0*pi - pi
        ENDDO
      ENDIF
      DO i=1,Nnodes
        P1 = randarray(3*(i-1)+1)
        P2 = randarray(3*(i-1)+2)
        P3 = randarray(3*(i-1)+3)
        !
        vorient(i,:,:) = Id_Matrix(:,:)  !unity matrix
        IF( twodim==0 .OR. twodim==3 ) THEN
          !Construct the rotation matrix around Z
          rotmat(:,:) = 0.d0
          rotmat(3,3) = 1.d0
          rotmat(1,1) = DCOS(P3)
          rotmat(1,2) = -1.d0*DSIN(P3)
          rotmat(2,1) = DSIN(P3)
          rotmat(2,2) = DCOS(P3)
          vorient(i,:,:) = rotmat(:,:)
        ENDIF
        IF( twodim==0 .OR. twodim==2 ) THEN
          !Construct the rotation matrix around Y
          rotmat(:,:) = 0.d0
          rotmat(2,2) = 1.d0
          rotmat(3,3) = DCOS(P2)
          rotmat(3,1) = -1.d0*DSIN(P2)
          rotmat(1,3) = DSIN(P2)
          rotmat(1,1) = DCOS(P2)
          vorient(i,:,:) = MATMUL( rotmat(:,:) , vorient(i,:,:) )
        ENDIF
        IF( twodim==0 .OR. twodim==1 ) THEN
          !Construct the rotation matrix around X
          rotmat(:,:) = 0.d0
          rotmat(1,1) = 1.d0
          rotmat(2,2) = DCOS(P1)
          rotmat(2,3) = -1.d0*DSIN(P1)
          rotmat(3,2) = DSIN(P1)
          rotmat(3,3) = DCOS(P1)
          vorient(i,:,:) = MATMUL( rotmat(:,:) , vorient(i,:,:) )
        ENDIF
      ENDDO
      !
      GOTO 250
      !
      !
    ELSEIF( StrDnCase(line(1:5))=="node " .OR. StrDnCase(line(1:5))=="grain" ) THEN
      !Check that the box was defined
      IF( .NOT.Hset ) THEN
        GOTO 820
      ENDIF
      Nnodes = Nnodes+1
      !
      IF( .NOT.ALLOCATED(randarray) ) THEN
        !Generate N random numbers, just in case one or more nodes have random orientation
        CALL GEN_NRANDNUMBERS( 3*SIZE(vnodes,1) , randarray )
        !randarray now contains N real numbers between 0 and 1
        !Multiply them by 2*pi and subtract pi to generate 3 angles alpha, beta and gamma
        IF( twodim>0 ) THEN
          randarray(:) = randarray(:)*2.d0*pi - pi
        ELSE
          DO i=1,Nnodes
            m = 3*(i-1) + 1
            n = 3*(i-1) + 2
            randarray(m) = randarray(m)*2.d0*pi - pi
            randarray(n) = randarray(n)*2.d0*pi - pi
          ENDDO
        ENDIF
      ENDIF
      !
      !Read position of that grain
      !Note: position may be given with respect to box dimension, e.g. "box/2"
      line = TRIM(ADJUSTL(line(6:)))
      i=SCAN(line,' ')
      temp = TRIM(ADJUSTL(line(:i)))
      CALL BOX2DBLE( H(1,:),temp,vnodes(Nnodes,1),status )
      IF( status>0 ) GOTO 830
      line = TRIM(ADJUSTL(line(i+1:)))
      i=SCAN(line,' ')
      temp = TRIM(ADJUSTL(line(:i)))
      CALL BOX2DBLE( H(2,:),temp,vnodes(Nnodes,2),status )
      IF( status>0 ) GOTO 830
      line = TRIM(ADJUSTL(line(i+1:)))
      i=SCAN(line,' ')
      temp = TRIM(ADJUSTL(line(:i)))
      CALL BOX2DBLE( H(3,:),temp,vnodes(Nnodes,3),status )
      line = TRIM(ADJUSTL(line(i+1:)))
      IF( status>0 ) GOTO 830
      !
      !Read crystallographic orientation of that grain
      !(can be explicitely given as Miller indices, or random)
      IF( line(1:6)=="random" ) THEN
        !Generated random parameters will be written into a file later
        outparam(4) = .TRUE.
        !
        IF( twodim>0 ) THEN
          !Pick a random angle from randarray
          P1 = randarray(3*(Nnodes-1)+1)
          !Get indices
          m = twodim
          n = twodim+1
          IF(n>3) n=n-3
          o = twodim+2
          IF(o>3) o=o-3
          !Construct the rotation matrix around short axis
          rotmat(:,:) = 0.d0
          rotmat(m,m) = 1.d0
          rotmat(n,n) = DCOS(P1)
          rotmat(n,o) = -1.d0*DSIN(P1)
          rotmat(o,n) = DSIN(P1)
          rotmat(o,o) = DCOS(P1)
          !Save rotation matrix in vorient
          vorient(Nnodes,:,:) = rotmat(:,:)
        ELSE
          !Generate a random rotation matrix for this grain, using random numbers generated before
          !Method is from "Fast random rotation matrices", James Arvo, Cornell University
          !NOTE: the distribution of orientations will be completely random only at the
          !     condition that all grain orientations are specified as "random".
          !     If specific angles or Miller indices are given explicitely by the user
          !     for some grains, then the distribution will not be completely random.
          P1 = randarray(3*(Nnodes-1)+1)
          P2 = randarray(3*(Nnodes-1)+2)
          P3 = randarray(3*(Nnodes-1)+3)
          !Compute vector V
          vector(1) = DSQRT(P3)*DCOS(P1)
          vector(2) = DSQRT(P3)*DSIN(P1)
          vector(3) = DSQRT(1.d0-P3)
          !Compute matrix R
          rotmat(:,:) = 0.d0
          rotmat(1,1) = DCOS(P2)
          rotmat(1,2) = DSIN(P2)
          rotmat(2,1) = -1.d0*DSIN(P2)
          rotmat(2,2) = DCOS(P2)
          rotmat(3,3) = 1.d0
          !Compute final rotation matrix:  M = ( 2*V^T*V - I ) * R
          vorient(Nnodes,:,:) = MATMUL( 2.d0*VECMAT(vector,vector) - Id_Matrix , rotmat )
        ENDIF
        !
      ELSE
        !The user provides 3 angles or 3 Miller indices
        READ(line,*,END=830,ERR=830) or1, or2, or3
        miller=.TRUE.
        IF( SCAN(or1,"[")>0 .OR. SCAN(or2,"[")>0 .OR. SCAN(or3,"[")>0 .OR.          &
          & SCAN(or1,"]")>0 .OR. SCAN(or2,"]")>0 .OR. SCAN(or3,"]")>0 .OR.          &
          & SCAN(or1,"_")>0 .OR. SCAN(or2,"_")>0 .OR. SCAN(or3,"_")>0 ) THEN
          !No ambiguity: it should be Miller indices (given in brackets and/or with underscores)
          miller=.TRUE.
        ELSEIF( SCAN(or1,"°")>0 .OR. SCAN(or2,"°")>0 .OR. SCAN(or3,"°")>0 .OR.      &
              & SCAN(or1,".")>0 .OR. SCAN(or2,".")>0 .OR. SCAN(or3,".")>0 ) THEN
          !No ambiguity: it is angles
          miller=.FALSE.
        ELSE
          !Ambiguous data: the user entered something like "110", is it an angle or Miller indices?
          IF( LEN_TRIM(or1)>=3 .AND. LEN_TRIM(or2)>=3 .AND. LEN_TRIM(or3)>=3 ) THEN
            !Try to interpret it as Miller indices, if it fails then it is angles
            miller=.TRUE.
            CALL INDEX_MILLER(or1,rotmat,j)
            IF(j>0) miller=.FALSE.
            CALL INDEX_MILLER(or2,rotmat,j)
            IF(j>0) miller=.FALSE.
            CALL INDEX_MILLER(or3,rotmat,j)
            IF(j>0) miller=.FALSE.
          ELSE
            !or1, or2 and/or or3 contain only 2 digits => consider they are angles
            miller=.FALSE.
          ENDIF
        ENDIF
        !
        IF( miller ) THEN
          !Read and interpret the Miller indices,
          !save the rotation matrix in vorient(Nnodes,:,:)
          CALL INDEX_MILLER(or1,rotmat(1,:),j)
          IF(j>0) GOTO 800
          CALL INDEX_MILLER(or2,rotmat(2,:),j)
          IF(j>0) GOTO 800
          CALL INDEX_MILLER(or3,rotmat(3,:),j)
          IF(j>0) GOTO 800
          CALL MILLER2ROTMAT(Id_Matrix,rotmat,vorient(Nnodes,1:3,1:3),j)
        ELSE
          !Read and interpret the angles,
          !save the rotation matrix in vorient(Nnodes,:,:)
          !Rotation will be done in the order X, Y, Z
          !Read the first angle from or1
          j=SCAN(or1,"°")
          IF(j>0) or1(j:j)=" "
          READ(or1,*,END=830,ERR=830) vector(1)
          !Read the second angle from or2
          j=SCAN(or2,"°")
          IF(j>0) or2(j:j)=" "
          READ(or2,*,END=830,ERR=830) vector(2)
          !Read the third angle from or3
          j=SCAN(or3,"°")
          IF(j>0) or3(j:j)=" "
          READ(or3,*,END=830,ERR=830) vector(3)
          !Make sure angles are between -180° and +180°, convert into radians
          DO j=1,3
            DO WHILE( vector(j)<=-180.d0 )
              vector(j) = vector(j)+360.d0
            ENDDO
            DO WHILE( vector(j)>180.d0 )
              vector(j) = vector(j)-360.d0
            ENDDO
            !Convert into radians
            vector(j) = DEG2RAD(vector(j))
          ENDDO
          !Construct rotation matrix
          CALL EULER2MAT_ZYX(vector(1),vector(2),vector(3),vorient(Nnodes,:,:))
          IF( .NOT.IS_ROTMAT(vorient(Nnodes,:,:)) ) THEN
            WRITE(msg,*) Nnodes
            WRITE(msg,*) "Node # "//TRIM(ADJUSTL(msg))//": matrix constructed from angles"
            CALL ATOMSK_MSG(2823,(/msg/),(/0.d0/))
          ENDIF
        ENDIF
      ENDIF
      !
      !
    ELSEIF( StrDnCase(line(1:6))=="random" ) THEN
      !Position and orientations of grains are random
      !Generated random parameters will be written into a file later
      outparam(4) = .TRUE.
      !Check that the box was defined
      IF( .NOT.Hset ) THEN
        GOTO 820
      ENDIF
      !Read total number of grains
      READ(line(7:),*,END=830,ERR=830) Nnodes
      IF(Nnodes<1) THEN
        CALL ATOMSK_MSG(4831,(/vfile/),(/0.d0/))
        nerr=nerr+1
        GOTO 1000
      ENDIF
      !Generate a list of 6*Nnodes random numbers
      CALL GEN_NRANDNUMBERS( 6*Nnodes , randarray )
      !randarray now contains 6*Nnodes real numbers between 0 and 1
      !The 3*Nnodes first random numbers are used to generate positions
      DO i=1,Nnodes
        vnodes(i,1) = randarray(3*(i-1)+1) * H(1,1)
        vnodes(i,2) = randarray(3*(i-1)+2) * H(2,2)
        vnodes(i,3) = randarray(3*(i-1)+3) * H(3,3)
      ENDDO
      !
      !The last 3*Nnodes random numbers are used to generate rotation matrices
      !Modify them to generate angles
      IF( twodim>0 ) THEN
        !Grains are rotated only around one axis
        !Multiply all random numbers by 2*pi to generate random angles
        randarray(3*Nnodes:) = randarray(3*Nnodes:)*2.d0*pi - pi
        DO i=1,Nnodes
          !Pick a random angle
          P1 = randarray(3*Nnodes+3*(i-1)+1)
          !Get indices
          m = twodim
          n = twodim+1
          IF(n>3) n=n-3
          o = twodim+2
          IF(o>3) o=o-3
          !Construct the rotation matrix around short axis
          rotmat(:,:) = 0.d0
          rotmat(m,m) = 1.d0
          rotmat(n,n) = DCOS(P1)
          rotmat(n,o) = -1.d0*DSIN(P1)
          rotmat(o,n) = DSIN(P1)
          rotmat(o,o) = DCOS(P1)
          !Save rotation matrix in vorient
          vorient(i,:,:) = rotmat(:,:)
        ENDDO
      ELSE
        !Generate a random rotation matrix for each grain
        !Method is from "Fast random rotation matrices", James Arvo, Cornell University
        DO i=1,Nnodes
          !Generate two random angles; third random number will be used below
          m = 3*Nnodes + 3*(i-1) + 1
          n = 3*Nnodes + 3*(i-1) + 2
          randarray(m) = randarray(m)*2.d0*pi - pi
          randarray(n) = randarray(n)*2.d0*pi - pi
          !Compute vector V
          P1 = randarray(3*Nnodes+3*(i-1)+3)
          vector(1) = DSQRT(P1)*DCOS(randarray(3*Nnodes+3*(i-1)+1))
          vector(2) = DSQRT(P1)*DSIN(randarray(3*Nnodes+3*(i-1)+1))
          vector(3) = DSQRT(1.d0-P1)
          !Compute matrix R
          rotmat(:,:) = 0.d0
          rotmat(1,1) = DCOS(randarray(3*Nnodes+3*(i-1)+2))
          rotmat(1,2) = DSIN(randarray(3*Nnodes+3*(i-1)+2))
          rotmat(2,1) = -1.d0*DSIN(randarray(3*Nnodes+3*(i-1)+2))
          rotmat(2,2) = DCOS(randarray(3*Nnodes+3*(i-1)+2))
          rotmat(3,3) = 1.d0
          !Compute final rotation matrix:  M = ( 2*V^T*V - I ) * R
          vorient(i,:,:) = MATMUL( 2.d0*VECMAT(vector,vector) - Id_Matrix , rotmat )
        ENDDO
      ENDIF
      !
      !Check for additional keyword after
      READ(line(7:),*,END=830,ERR=830) temp  !that's Nnodes
      temp = ADJUSTL(temp)
      m = LEN_TRIM(temp)
      temp = TRIM(ADJUSTL(line(7+m+1:)))
      IF( LEN_TRIM(temp)>0 ) THEN
        IF( temp=="linear" ) THEN
          !Adjust node positions to generate linear gradient along given direction
          DO i=1,Nnodes
            vnodes(i,1) = vnodes(i,1)**2
          ENDDO
        ENDIF
      ENDIF
      !
      !
      IF( verbosity==4 ) THEN
        !Write angles into a file for visualization/debug purposes
        !Rotation matrices are applied to vector [100]
        !File is in XYZ format to be visualized with VESTA, gnuplot or other softwares
        OPEN(UNIT=43,FILE="atomsk_angles.xyz",STATUS="UNKNOWN")
        WRITE(43,*) Nnodes
        WRITE(43,*) "# Distribution of random angles generated by Atomsk"
        DO i=1,Nnodes
          vector(1) = vorient(i,1,1)
          vector(2) = vorient(i,2,1)
          vector(3) = vorient(i,3,1)
          WRITE(43,'(a3,3(f16.3,2X))') "1  ", vector(1), vector(2), vector(3)
        ENDDO
        CLOSE(43)
      ENDIF
      !
    ELSEIF( StrDnCase(line(1:8))=="gradient" ) THEN
      !User wants to create a gradient along a direction between P1 and P2
      IF( ALLOCATED(vnodes) .AND. SIZE(vnodes,1)>1 ) THEN
        READ(line(9:),*,END=800,ERR=800) or1
        IF( StrDnCase(or1)=="linear" ) THEN
          !Read direction of linear gradient and min/max values
          READ(line(9:),*,END=800,ERR=800) or1, or2, P1, grad1, P2, grad2
          IF( StrDnCase(or2)=='x' ) THEN
            j=1
          ELSEIF( StrDnCase(or2)=='y' ) THEN
            j=2
          ELSEIF( StrDnCase(or2)=='y' ) THEN
            j=3
          ENDIF
          !Compute gradient and normalize it
          !distance = (grad2-grad1)/(P2-P1)
          !Modify node positions (keep positions outside of gradient zone)
          DO i=1,SIZE(vnodes,1)
            IF( vnodes(i,j)>=P1 .AND. vnodes(i,j)<=P2 ) THEN
              vnodes(i,j) = vnodes(i,j)*(grad2-grad1)/(P2-P1)
            ENDIF
          ENDDO
        ELSEIF( StrDnCase(or1)=="radial" ) THEN
          !Modify node positions (keep positions outside of gradient zone)
          vector(:) = 0.5d0*(H(1,:)+H(2,:)+H(3,:))
          DO i=1,SIZE(vnodes,1)
            distance = VECLENGTH( vnodes(i,1:3) - vector(1:3) )
            vnodes(i,1:3) = (vnodes(i,1:3) - vector(:)) + distance*vnodes(i,1:3)
          ENDDO
        ENDIF
        !
      ELSE  !i.e. if vnodes unallocated or has only 1 node

      ENDIF
      !
      GOTO 250
      !
      !
    ENDIF
  ENDIF
ENDDO
250 CONTINUE
CLOSE(31)
!
!
300 CONTINUE
!Check that number of nodes is not zero
IF(Nnodes<1) THEN
  CALL ATOMSK_MSG(4831,(/vfile/),(/0.d0/))
  nerr=nerr+1
  GOTO 1000
ENDIF
!
!Make sure that all nodes are inside the final box H(:,:)
CALL CART2FRAC(vnodes,H)
DO i=1,SIZE(vnodes,1) !loop on all nodes
  m=0
  DO j=1,3  !loop on xyz
    DO WHILE( vnodes(i,j)>=1.d0 )
      vnodes(i,j) = vnodes(i,j)-1.d0
      m=m+1
    ENDDO
    DO WHILE( vnodes(i,j)<0.d0 )
      vnodes(i,j) = vnodes(i,j)+1.d0
      m=m+1
    ENDDO
  ENDDO
  IF(m>0) THEN
    !This node was wrapped: display message
    P1 = vnodes(i,1)*H(1,1)
    P2 = vnodes(i,2)*H(2,2)
    P3 = vnodes(i,3)*H(3,3)
    nwarn=nwarn+1
    CALL ATOMSK_MSG(4716,(/''/),(/P1,P2,P3,DBLE(i)/))
  ENDIF
ENDDO
CALL FRAC2CART(vnodes,H)
!
!Check that nodes are not too close to one another
DO i=1,SIZE(vnodes,1)-1 !loop on all nodes
  DO j=i+1,SIZE(vnodes,1)  !loop on all nodes
    !Compute distance between the two nodes #i and #j
    distance = VECLENGTH( vnodes(i,1:3)-vnodes(j,1:3) )
    IF( distance<clearance+0.1d0 ) THEN
      !Nodes #i and #j are extremely close
      IF( outparam(3) ) THEN
        !User asked for random positions: proximity was caused by random generator
        !Correct it silently (user will probably not complain, since he asked for randomness anyway)
        vnodes(j,1:3) = vnodes(j,1:3) + (clearance+0.1d0)*(vnodes(j,1:3)-vnodes(i,1:3))
      ELSE
        !These positions were defined by user: display error and exit
        !(do not try to automatically fix something that was explicitely defined by user)
        nerr = nerr+1
        CALL ATOMSK_MSG(4832,(/''/),(/DBLE(i),DBLE(j)/))
        GOTO 1000
      ENDIF
    ELSEIF( distance<2.d0 ) THEN
      !Nodes #i and #j are very close to one another: display warning
      nwarn = nwarn+1
      CALL ATOMSK_MSG(4718,(/''/),(/DBLE(i),DBLE(j),distance/))
    ENDIF
  ENDDO
ENDDO
GOTO 1000
!
!
!
800 CONTINUE
status = 1
nerr=nerr+1
GOTO 1000
!
810 CONTINUE
status = 2
nerr=nerr+1
CALL ATOMSK_MSG(2813,(/TRIM(temp)/),(/0.d0/))
GOTO 1000
!
820 CONTINUE
status = 3
nerr=nerr+1
CALL ATOMSK_MSG(4820,(/TRIM(temp)/),(/0.d0/))
GOTO 1000
!
830 CONTINUE
status = 4
nerr=nerr+1
CALL ATOMSK_MSG(1801,(/TRIM(vfile)/),(/DBLE(linenumber)/))
GOTO 1000
!
!
!
1000 CONTINUE
IF( verbosity==4 ) THEN
  WRITE(msg,*) " Final parameters:"
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,'(a8,3(f9.3,a3))') "    Box: ", H(1,1), " x ", H(2,2), " x ",  H(3,3)
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "   Nnodes = ", Nnodes
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
  WRITE(msg,*) "   Clearance = ", clearance
  CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
ENDIF
WRITE(msg,*) "EXITING POLYX_READ_PARAM (", status, ")"
CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
!
!
!
END SUBROUTINE POLYX_READ_PARAM
!
!
END MODULE polyx_readparam
