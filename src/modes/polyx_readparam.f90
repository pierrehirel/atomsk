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
!* Last modification: P. Hirel - 24 July 2026                                     *
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
USE readin
USE subroutines
USE writeout
!
CONTAINS
!
!
SUBROUTINE POLYX_READ_PARAM(vfile,Hs,Ps,Ss,AUXNAMESs,AUXs,NPs,idseed,vnodes,vorient,H,clearance,status)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: vfile   !name of file containing parameters for Voronoi construction
CHARACTER(LEN=128):: lattice  !if grains are organized according to a lattice
CHARACTER(LEN=128):: or1, or2, or3
CHARACTER(LEN=4096):: line, prefix
CHARACTER(LEN=4096):: msg, temp
CHARACTER(LEN=4096):: outparamfile  !file where grain parameters are written (if some parameters equal "random")
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats !list of formats to output
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: comment     !comments of a seed
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMESt   !names of aux.prop. of a seed
LOGICAL:: doshells        !are shells present?
LOGICAL:: Hset            !are the box vectors H(:,:) defined?
LOGICAL:: miller          !are Miller indices given? (if no then angles are given)
LOGICAL:: seedasfile      !was seed already defined in a file?
LOGICAL,DIMENSION(4):: outparam  !Are keywords (1) "node", (2) "lattice", (3) "random" used?
                                             !(4) Must parameters be saved in a text file?
INTEGER:: i, j, m, n, o
INTEGER:: linenumber    !line number
INTEGER:: seed_id       !index of current seed
INTEGER:: status        !=0 if successful, >0 otherwise
INTEGER:: twodim        !=0 if system is 3-D, =1,2,3 if system is thin along x, y, z
INTEGER:: Nnodes      !number of nodes
REAL(dp):: distance    !distance between two points
REAL(dp):: grad1, grad2!gradient values
REAL(dp):: P1, P2, P3  !temporary position
REAL(dp):: Volume, Vmin !min. volume occupied by a grain
REAL(dp),INTENT(INOUT):: clearance   !clear atoms that close to the GB
REAL(dp),DIMENSION(3):: vector    !vector between an atom and a node
REAL(dp),DIMENSION(3,3):: Ht      !cell vectors of a seed
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: Pt, St !positions of atoms/shells of a seed
REAL(dp),DIMENSION(:,:),ALLOCATABLE:: AUXt   !auxiliary properties of a seed
REAL(dp),DIMENSION(3,3):: rotmat  !rotation matrix
REAL(dp),DIMENSION(:),ALLOCATABLE:: randarray   !random numbers
!
CHARACTER(LEN=128),DIMENSION(:),ALLOCATABLE:: AUXNAMESs       !Names of aux.prop. of all seeds
CHARACTER(LEN=4096),DIMENSION(:),ALLOCATABLE:: seedfiles      !Names of seed file for each grain
INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: NPs          !Number of atoms for each seed
INTEGER,DIMENSION(:),ALLOCATABLE,INTENT(INOUT):: idseed       !Seed index for each node
REAL(dp),DIMENSION(3,3),INTENT(OUT):: H                       !Base vectors of the final supercell
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT):: Hs      !Base vectors of all seeds
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT):: Ps, Ss  !positions of atoms/shells in seeds
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(INOUT):: AUXs    !aux.prop. in seeds
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: vnodes      !cartesian coordinate of each node
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(OUT):: vorient   !rotation matrix for each node
!
!Initialize variables
outparamfile = TRIM_EXT(vfile)//"_param.txt"
doshells = .FALSE.
Hset = .FALSE.
miller = .FALSE.
outparam(:) = .FALSE.
linenumber = 0
seed_id = 0
status = 0
twodim = 0
H(:,:) = 0.d0
IF(ALLOCATED(AUXNAMESs)) DEALLOCATE(AUXNAMESs)
IF(ALLOCATED(seedfiles)) DEALLOCATE(seedfiles)
IF(ALLOCATED(vnodes)) DEALLOCATE(vnodes)
IF(ALLOCATED(vorient)) DEALLOCATE(vorient)
!
!Check if a seed was already defined as input in arrays Hs, Ps, Ss
seedasfile = .FALSE.
IF( ALLOCATED(Hs) .AND. ALLOCATED(Ps) ) THEN
  IF( ANY(DABS(Hs(:,:,:))>1.d-3) .AND. SIZE(Ps,1)>0 .AND. SIZE(Ps,2)>0 .AND. ANY(DABS(Ps(:,:,:))>1.d-3) ) THEN
    !A seed was already provided as a "seed file"
    seedasfile = .TRUE.
  ELSE
    !Something is wrong with Ps: deallocate everything
    DEALLOCATE(Ps)
    IF(ALLOCATED(Hs)) DEALLOCATE(Hs)
    IF(ALLOCATED(Ss)) DEALLOCATE(Ss)
    IF(ALLOCATED(AUXs)) DEALLOCATE(AUXs)
    IF(ALLOCATED(AUXNAMESs)) DEALLOCATE(AUXNAMESs)
  ENDIF
ELSE
  !Ps is not allocated: make sure the other arrays are unallocated too
  IF(ALLOCATED(Hs)) DEALLOCATE(Hs)
  IF(ALLOCATED(Ps)) DEALLOCATE(Ps)
  IF(ALLOCATED(Ss)) DEALLOCATE(Ss)
  IF(ALLOCATED(AUXs)) DEALLOCATE(AUXs)
  IF(ALLOCATED(AUXNAMESs)) DEALLOCATE(AUXNAMESs)
ENDIF
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
      Hset=.TRUE.
      !
      IF(ALLOCATED(Hs)) THEN
        twodim = CHECK2DIM(Hs,H)
      ENDIF
      !
      WRITE(msg,*) "            twodim = ", twodim
      CALL ATOMSK_MSG(999,(/TRIM(msg)/),(/0.d0/))
      !
      IF( twodim>0 ) THEN
        IF( DABS(H(twodim,twodim))<1.d-3 ) THEN
          !Set the dimension of final cell to that of the seed along short axis
          H(twodim,twodim) = MAXVAL(Hs(:,twodim,twodim))
        ENDIF
      ENDIF
      !
      Volume = H(1,1) * H(2,2) * H(3,3)
      IF( twodim==0 .AND. Volume < 8000.d0 ) THEN
        !Final cell is very small: display warning
        nwarn = nwarn+1
        CALL ATOMSK_MSG(4719,(/""/),(/Volume/))
      ENDIF
      !
    ELSEIF( StrDnCase(line(1:5))=="seed " ) THEN
      !Name of input file containing the seed (there can be only one)
      IF( ALLOCATED(Ps) ) THEN
        !seed was already defined before
        PRINT*, "Warning: only one seed can be defined"
      ELSE
        !Read data from file
        temp = TRIM(ADJUSTL(line(5:)))
        CALL READ_AFF(temp,Ht,Pt,St,comment,AUXNAMESt,AUXt)
        IF(nerr>0 ) GOTO 1000
        IF( .NOT.ALLOCATED(Pt) .OR. SIZE(Pt,1)<1 ) THEN
          CALL ATOMSK_MSG(804,(/''/),(/0.d0/))
          nerr=nerr+1
          GOTO 1000
        ENDIF
        IF(ALLOCATED(comment)) DEALLOCATE(comment)
        !
        IF(ALLOCATED(idseed)) DEALLOCATE(idseed)
        ALLOCATE(idseed(1))
        idseed(1) = 1
        IF(ALLOCATED(NPs)) DEALLOCATE(NPs)
        ALLOCATE(NPs(1))
        NPs(1) = SIZE(Pt,1)
        !
        !Copy this data to seeds arrays: Hs, Ps, Ss, AUXs
        !NOTE: in this case the first dimension of these arrays will be =1
        CALL RESIZE_DBLEARRAY3(Hs,1,3,3,status)
        IF( status>0 ) GOTO 830
        Hs(1,:,:) = Ht(:,:)
        !Copy atom positions
        n = MAX( SIZE(Ps,2) , SIZE(Pt,1) )
        CALL RESIZE_DBLEARRAY3(Ps,1,n,4,status)
        IF( status>0 ) GOTO 830
        m = SIZE(Pt,1)
        Ps(1,1:m,:) = Pt(1:m,:)
        !Copy shells
        IF( ALLOCATED(St) .AND. SIZE(St,1)==SIZE(Pt,1) ) THEN
          CALL RESIZE_DBLEARRAY3(Ss,1,n,4,status)
          IF( status>0 ) GOTO 830
          IF( ALLOCATED(St) ) THEN
            m = SIZE(St,1)
            Ss(1,1:m,:) = St(1:m,:)
          ENDIF
        ENDIF
        !Copy auxiliary properties
        IF( ALLOCATED(AUXt) .AND. SIZE(AUXt,1)==SIZE(Pt,1) ) THEN
          ALLOCATE(AUXNAMESs(SIZE(AUXNAMESt)))
          AUXNAMESs(:) = AUXNAMESt(:)
          CALL RESIZE_DBLEARRAY3(AUXs,1,n,SIZE(AUXt,2),status)
          IF( status>0 ) GOTO 830
          IF( ALLOCATED(AUXt) ) THEN
            m = SIZE(AUXt,1)
            AUXs(1,1:m,:) = AUXt(1:m,:)
          ENDIF
        ENDIF
        !Free temporary arrays
        IF(ALLOCATED(Pt)) DEALLOCATE(Pt)
        IF(ALLOCATED(St)) DEALLOCATE(St)
        IF(ALLOCATED(AUXt)) DEALLOCATE(AUXt)
        IF(ALLOCATED(AUXNAMESt)) DEALLOCATE(AUXNAMESt)
        !
        twodim = CHECK2DIM(Hs,H)
        !
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
!Check if system is pseudo-2D
twodim = CHECK2DIM(Hs,H)
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
          ELSEIF( StrDnCase(or2)=='z' ) THEN
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
          !Read sigma value for Gaussian dstribution
          READ(line(9:),*,END=800,ERR=800) or1, P2
          !Modify node positions (keep positions outside of gradient zone)
          P1 = 1.d0 / DSQRT(2.d0*pi*P2**2)
          vector(:) = 0.5d0*(H(1,:)+H(2,:)+H(3,:))
          IF(twodim>0) vector(twodim) = 0.5d0*H(twodim,twodim)
          DO i=1,SIZE(vnodes,1)
            distance = VECLENGTH( vnodes(i,1:3) - vector(1:3) )
            vnodes(i,1:3) = vnodes(i,1:3) * DEXP( -distance**2 / (2.d0*P2**2) )
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
!Check one last time if system is pseudo-2D
twodim = CHECK2DIM(Hs,H)
!
IF( twodim > 0 ) THEN
  !System is pseudo 2-D => place all nodes on the same plane
  !in the middle of the cell along the short direction
  vnodes(:,twodim) = 0.5d0*H(twodim,twodim)
ENDIF
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
!
IF( outparam(4) .AND. ofu.NE.6 ) THEN
  !Write positions and orientations in a parameter file
  OPEN(41,FILE=outparamfile,STATUS="UNKNOWN")
  WRITE(41,'(a62)') "# Random positions and rotations of grains generated by Atomsk"
  WRITE(41,'(a66)') "# This parameter file can be used to generate the same polycrystal"
  WRITE(41,'(a4,3f16.6)') "box ", H(1,1), H(2,2), H(3,3)
  DO i=1,SIZE(vorient,1)
    !Compute corresponding rotation vectors and write them into parameter file
    CALL MAT2EULER_ZYX(vorient(i,:,:),P1,P2,P3)
    WRITE(41,'(a5,6f16.6)') "node ", vnodes(i,:), RAD2DEG(P1), RAD2DEG(P2), RAD2DEG(P3)
  ENDDO
  CLOSE(41)
  !Random parameters were written into a file => inform user
  temp = "parameter"
  CALL ATOMSK_MSG(3002,(/outparamfile,temp/),(/0.d0/))
ENDIF
!
IF( ofu.NE.6 ) THEN
  !Write positions of nodes into a file
  ALLOCATE(comment(1))
  comment(1) = "# Positions of nodes"
  ALLOCATE(outfileformats(1))
  outfileformats(1) = "xsf"
  ALLOCATE(Pt(SIZE(vnodes,1),4))
  DO i=1,SIZE(vnodes,1)
    Pt(i,1:3) = vnodes(i,1:3)
    Pt(i,4) = 1.d0  !nodes are given the atomic number of hydrogen
  ENDDO
  CALL NAME_OUTFILE(prefix,temp,"     ")
  i=SCAN(temp,".",BACK=.TRUE.)  !remove trailing dot
  temp = temp(1:i-1)
  outparamfile = TRIM(ADJUSTL(temp))//"_nodes.xsf"
  CALL WRITE_AFF(outparamfile,outfileformats,H,Pt,St,comment,AUXNAMESt,AUXt)
  IF(ALLOCATED(Pt)) DEALLOCATE(Pt)
  IF(ALLOCATED(comment)) DEALLOCATE(comment)
ENDIF
!
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
IF(ALLOCATED(seedfiles)) DEALLOCATE(seedfiles)
!
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
!********************************************************
! CHECK2DIM
! This subroutine checks if the final supercell H(:,:)
! can be considered as "2-dimensional" (or pseudo-2D), by
! comparing its length with existing see vectors.
!********************************************************
FUNCTION CHECK2DIM(Hs,H) RESULT(twodim)
!
IMPLICIT NONE
INTEGER:: i
INTEGER:: twodim  !=0 if system is 3-D, =1,2,3 if system is thin along x, y, z
REAL(dp):: dmin
REAL(dp),DIMENSION(3,3),INTENT(IN):: H     !Base vectors of the final supercell
REAL(dp),DIMENSION(:,:,:),ALLOCATABLE,INTENT(IN):: Hs  !Base vectors of all seeds
!
twodim = 0
!
!Check if final supercell H is set
IF( ANY(DABS(H(:,:))>1.d-6) ) THEN
  !Check if at least one seed cell is defined
  IF( ALLOCATED(Hs) ) THEN
    IF( SIZE(Hs,1)>=1 .AND. SIZE(Hs,2)>=3 .AND. SIZE(Hs,3)>=3 .AND. ANY(DABS(Hs(:,:,:))>1.d-6) ) THEN
      !Compare each dimension of H with Hs
      DO i=1,3
        IF( H(i,i)<=MAXVAL(Hs(:,i,i)) ) THEN
          twodim = i
        ENDIF
      ENDDO
    ENDIF
    !
  ELSE
    !Hs unallocated, no seed is defined yet
    IF( H(1,1)<6.d0 ) THEN
      twodim = 1
    ELSEIF( H(2,2)<6.d0 ) THEN
      twodim = 2
    ELSEIF( H(3,3)<6.d0 ) THEN
      twodim = 3
    ENDIF
  ENDIF
ENDIF
!
END FUNCTION CHECK2DIM
!
!
END MODULE polyx_readparam
