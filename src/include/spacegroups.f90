MODULE spacegroups
!
!**********************************************************************************
!*  SPACEGROUPS                                                                   *
!**********************************************************************************
!* This module provides space group symmetry data, together with the spacegroup   *
!* number, name and Patterson space group number.                                 *
!* The data was taken from <http://cci.lbl.gov/cctbx/explore_symmetry.html>       *
!* and edited by Dr. J. Barthel, RWTH-Aachen University, Aachen, Germany          *
!**********************************************************************************
!* (C) July 2015 - Juri Barthel                                                   *
!*     Gemeinschaftslabor fuer Elektronenmikroskopie                              *
!*     RWTH Aachen (GERMANY)                                                      *
!*     ju.barthel@fz-juelich.de                                                   *
!* Last modification: J. Barthel - 04 Aug. 2015                                   *
!**********************************************************************************
!* Notes on the model and on how to use it.                                       *
!* The space group data can be accessed from its number 1 ... 230 and from its    *
!* name, which is the Hermann-Mauguin symbol, e.g. P1, F4-3m, I4/mcm, and so on.  *
!* The data is stored in the allocatable arrays:                                  *
!* sg_name, sg_patn, sg_symnum, and sg_symmetry.                                  *
!* We use allocatable arrays here in order to not overload the stack.             *
!* On instantiation, the module data is not allocated by default.                 *
!* Call the routine SG_INIT to allocate and initialize the data arrays.           *
!* The total ammount of allocated heap is about 1.4MB.                            *
!* Calling any of the data access routines (see below) with an uninitialized      *
!* module will cause an auto-initialization.                                      *
!* Call the routine SG_UNINIT to deallocate the data arrays.                      *
!* Call the routine SG_ISREADY to probe the initialization state of the module.   *
!* Data access routines can be used to retrieve the space group data. In order to *
!* to create a fail-safe module, the actual arrays cannot be accessed directly.   *
!* The module should never cause access violations. Each data access routine will *
!* return a success-code, such that the calling routines may identify problems.   *
!* None of the routines writes messages to an I/O unit.                           *
!**********************************************************************************
!* NOTE ON COMPILATION (P. Hirel, Jan.2016):                                      *
!* Compiling this module with gfortran 4.4 and the flag "-ftree-pre" results in   *
!* a huge workload of RAM, and eventually a compiler internal error.              *
!* I was unable to determine which routine caused that. The flag "-free-pre"      *
!* is activated for "-O2" and higher, therefore to solve this problem             *
!* please compile with "-O1" or with "-fno-tree-pre".                             *
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
!
IMPLICIT NONE
!
SAVE
!
!* general module routines
PUBLIC::SG_ISREADY          ! returns the initialization state (0=not initialized)
PUBLIC::SG_INIT             ! initializes the module, allocates memory
PUBLIC::SG_UNINIT           ! de-initializes the module, deallocates memory
!* module data access routines
!* - routines retrieving data based on the space group number (index 1 ... 230)
PUBLIC::SG_NUMGETNAME       ! returns the name of a space group (Hermann-Mauguin symbol)
PUBLIC::SG_NUMGETPATN       ! returns the Patterson space group number
PUBLIC::SG_NUMGETSYMNUM     ! returns the number of symmetry operations
PUBLIC::SG_NUMGETSYMOP      ! returns a specific symmetry operation string
PUBLIC::SG_NUMGETSYMOPS     ! returns all symmetry operation strings
!* - routines retrieving data based on the space group name (Hermann-Mauguin symbol)
PUBLIC::SG_NAMGETNUM        ! returns the space group number (index 1 ... 230)
PUBLIC::SG_NAMGETPATN       ! returns the Patterson space group number
PUBLIC::SG_NAMGETSYMNUM     ! returns the number of symmetry operations
PUBLIC::SG_NAMGETSYMOP      ! returns a specific symmetry operation string
PUBLIC::SG_NAMGETSYMOPS     ! returns all symmetry operation strings
!
!* module parameter declarations
INTEGER,PUBLIC,PARAMETER::sg_nummax=230     ! maximum number of space groups
!                                           ! = table size after initialization
INTEGER,PRIVATE,PARAMETER::sg_opmax=192     ! maximum number of symmetry operations
!                                           ! = table size for the operator strings
INTEGER,PUBLIC,PARAMETER::sg_soplen=32      ! declaration string length
!
!* module data array declarations (all private)
CHARACTER(len=sg_soplen),PRIVATE,ALLOCATABLE::sg_name(:)        ! space group names
CHARACTER(len=sg_soplen),PRIVATE,ALLOCATABLE::sg_symmetry(:,:)  ! space group symmetry operation strings
INTEGER,PRIVATE,ALLOCATABLE::sg_patn(:)     ! Patterson space group number
INTEGER,PRIVATE,ALLOCATABLE::sg_symnum(:)   ! Number of space group symmetry operations
!
!
CONTAINS
!
!
!
SUBROUTINE SG_NUMGETNAME(nsgnum,strname,nchk)
!
! Returns the space group H-M symbol "strname" for the
! space group number "nsgnum".
! In case of an error "strname" returns '' and
! an error code is returned by "nchk".
! Possible error codes:
! - nchk = -1: a failure in the module initialization
! - nhck = -2: invalid space group number (out of range 1 - 230)
! In case of success, nchk return 1.
! The name string variable "strname" MUST BE of the form
!     CHARACTER(len=sg_soplen)::strname
!
IMPLICIT NONE
!
INTEGER,INTENT(IN)::nsgnum
CHARACTER(len=sg_soplen),INTENT(OUT)::strname
INTEGER,INTENT(OUT)::nchk
!
INTEGER::state
!
state=0
nchk=0
strname=''
IF (nsgnum<1 .OR. nsgnum>sg_nummax) GOTO 802 ! invalid sg index
!
CALL SG_ISREADY(state)
IF (0==state) CALL SG_INIT(state)
IF (0==state) GOTO 801 ! initialization error
!
strname=sg_name(nsgnum)
nchk=1
GOTO 1000
!
! error handling
801 CONTINUE ! module initialization problem
nchk = -1
GOTO 900
802 CONTINUE ! invalid space group name
nchk = -2
GOTO 900
!
! error clean up, go here only in case of a routine error!
900 CONTINUE
strname=''
!
! finish
1000 CONTINUE
!
END SUBROUTINE SG_NUMGETNAME
!
!
!
!
SUBROUTINE SG_NUMGETPATN(nsgnum,npatn,nchk)
!
! Returns the Patterson space group number "npatn" for the
! space group number "nsgnum".
! In case of an error "npatn" returns 0 and
! an error code is returned by "nchk".
! Possible error codes:
! - nchk = -1: a failure in the module initialization
! - nhck = -2: invalid space group number (out of range 1 - 230)
! In case of success, nchk return 1.
!
IMPLICIT NONE
!
INTEGER,INTENT(IN)::nsgnum
INTEGER,INTENT(OUT)::npatn,nchk
!
INTEGER::state
!
state=0
nchk=0
npatn=0
IF (nsgnum<1 .OR. nsgnum>sg_nummax) GOTO 802 ! invalid sg index
!
CALL SG_ISREADY(state)
IF (0==state) CALL SG_INIT(state)
IF (0==state) GOTO 801 ! initialization error
!
npatn=sg_patn(nsgnum)
nchk=1
GOTO 1000
!
! error handling
801 CONTINUE ! module initialization problem
nchk = -1
GOTO 900
802 CONTINUE ! invalid space group name
nchk = -2
GOTO 900
!
! error clean up, go here only in case of a routine error!
900 CONTINUE
npatn=0
!
! finish
1000 CONTINUE
!
END SUBROUTINE SG_NUMGETPATN
!
!
!
!
SUBROUTINE SG_NUMGETSYMNUM(nsgnum,nsymnum,nchk)
!
! Returns the number of symmetry operations "nsymnum" for the
! space group number "nsgnum".
! In case of an error "nsymnum" returns 0 and
! an error code is returned by "nchk".
! Possible error codes:
! - nchk = -1: a failure in the module initialization
! - nhck = -2: invalid space group number (out of range 1 - 230)
! In case of success, nchk return 1.
!
IMPLICIT NONE
!
INTEGER,INTENT(IN)::nsgnum
INTEGER,INTENT(OUT)::nsymnum,nchk
!
INTEGER::state
!
state=0
nchk=0
nsymnum=0
IF (nsgnum<1 .OR. nsgnum>sg_nummax) GOTO 802 ! invalid sg index
!
CALL SG_ISREADY(state)
IF (0==state) CALL SG_INIT(state)
IF (0==state) GOTO 801 ! initialization error
!
nsymnum=sg_symnum(nsgnum)
nchk=1
GOTO 1000
!
! error handling
801 CONTINUE ! module initialization problem
nchk = -1
GOTO 900
802 CONTINUE ! invalid space group name
nchk = -2
GOTO 900
!
! error clean up, go here only in case of a routine error!
900 CONTINUE
nsymnum=0
!
! finish
1000 CONTINUE
!
END SUBROUTINE SG_NUMGETSYMNUM
!
!
!
!
SUBROUTINE SG_NUMGETSYMOP(nsgnum,isymop,strsymop,nchk)
!
! Returns the symmetry operation string number "isymop" for the
! space group number "nsgnum" in the string "strsymop".
! The variable "nchk" returns <=0 in case of a routine error
! or 1 in case of success.
! Possible causes for a routine errors:
! - nchk = -1: a failure in the module initialization
! - nhck = -2: an invalid space group number (out of range 1 - 230)
! - nhck = -3: an invalid symmetry operation index
! The string "strsymop" MUST BE of the form
!     CHARACTER(len=sg_soplen)::strsymop
!
IMPLICIT NONE
!
INTEGER,INTENT(IN)::nsgnum,isymop
CHARACTER(len=sg_soplen),INTENT(OUT)::strsymop
INTEGER,INTENT(OUT)::nchk
!
INTEGER::state
INTEGER::nsymnum
!
nchk=0
nsymnum=0
strsymop=''
IF (nsgnum<1 .OR. nsgnum>sg_nummax) GOTO 802 ! false index error
!
CALL SG_ISREADY(state)
IF (0==state) CALL SG_INIT(state)
IF (0==state) GOTO 801 ! initialization error
!
! space group identified, get all the data
nsymnum = sg_symnum(nsgnum)
IF (isymop<1 .OR. isymop>nsymnum) GOTO 803 ! invalid operation index.
! all is ok now.
strsymop = sg_symmetry(isymop,nsgnum) ! get the string.
nchk=1 ! success!
GOTO 1000
!
! error handling
801 CONTINUE ! module initialization problem
nchk = -1
GOTO 900
802 CONTINUE ! invalid space group name
nchk = -2
GOTO 900
803 CONTINUE ! invalid index of symmetry operation
nchk = -3
GOTO 900
!
! error clean up, go here only in case of a routine error!
900 CONTINUE
strsymop=''
!
! finish
1000 CONTINUE
!
END SUBROUTINE SG_NUMGETSYMOP
!
!
!
!
SUBROUTINE SG_NUMGETSYMOPS(nsgnum,strsymops,nsymnum,nchk)
!
! Returns the list of all symmetry operation strings for the
! space group number "nsgnum" in the array "strsymops".
! The number of symmetry operation strings is returned in
! the variable nsymnum, which should be identical to
! SIZE(strsymops) on the routine exit in case of success.
! The variable "nchk" returns <=0 in case of a routine error
! or 1 in case of success.
! Possible causes for a routine errors:
! - nchk =  0: an unknown error
! - nchk = -1: a failure in the module initialization
! - nhck = -2: an invalid space group number (out of range 1 - 230)
! The array "strsymops" MUST BE an allocatable array of the
! form CHARACTER(len=sg_soplen),ALLOCATABLE::strsymops(:)
! The array will be reallocated by the routine. Or returns
! as unallocated in case of a failure.
!
IMPLICIT NONE
!
INTEGER,INTENT(IN)::nsgnum
CHARACTER(len=sg_soplen),DIMENSION(:),ALLOCATABLE::strsymops
INTEGER,INTENT(OUT)::nsymnum,nchk
!
INTEGER::nalloc
INTEGER::state
!
nchk=0
nalloc=0
nsymnum=0
IF (ALLOCATED(strsymops)) DEALLOCATE(strsymops,stat=nalloc)
IF (nsgnum<1 .OR. nsgnum>sg_nummax) GOTO 802 ! false index error
!
CALL SG_ISREADY(state)
IF (0==state) CALL SG_INIT(state)
IF (0==state) GOTO 801 ! initialization error
!
! space group identified, get all the data
nsymnum = sg_symnum(nsgnum)
ALLOCATE( strsymops(1:nsymnum), stat=nalloc )
IF (nalloc.NE.0) GOTO 800 ! This should not happen.
strsymops='' ! init.
strsymops(1:nsymnum) = sg_symmetry(1:nsymnum,nsgnum) ! copy.
nchk=1 ! got it, success!
GOTO 1000
!
! error handling
800 CONTINUE ! unknown error
nchk = 0
GOTO 900
801 CONTINUE ! module initialization problem
nchk = -1
GOTO 900
802 CONTINUE ! invalid space group name
nchk = -2
GOTO 900
!
! error clean up, go here only in case of a routine error!
900 CONTINUE
IF (ALLOCATED(strsymops)) DEALLOCATE(strsymops,stat=nalloc)
nsymnum=0
!
! finish
1000 CONTINUE
!
END SUBROUTINE SG_NUMGETSYMOPS
!
!
!
!
SUBROUTINE SG_NAMGETNUM(strname,nsgnum)
!
! Returns the space group number "nsgnum" for the
! space group H-M symbol "strname".
! In case of an error "nsgnum" returns 0.
!
IMPLICIT NONE
!
CHARACTER(len=*),INTENT(IN)::strname
INTEGER,INTENT(OUT)::nsgnum
!
INTEGER::state, i
CHARACTER(len=2*sg_soplen)::temp
!
state=0
nsgnum=0
temp=ADJUSTL(strname)
IF (LEN_TRIM(temp)<2) RETURN ! false name error
!
CALL SG_ISREADY(state)
IF (0==state) CALL SG_INIT(state)
IF (0==state) RETURN ! initialization error
!
! Get the index of the space group from its name.
! The name must be written correctly. Case sensitive!
DO i=1, sg_nummax
  IF (TRIM(sg_name(i))==TRIM(temp)) THEN
    nsgnum=i
    EXIT
  ENDIF
ENDDO
!
END SUBROUTINE SG_NAMGETNUM
!
!
!
!
SUBROUTINE SG_NAMGETPATN(strname,npatn,nchk)
!
! Returns the Patterson space group number "npatn" for the
! space group H-M symbol "strname".
! The variable "nchk" returns <=0 in case of a routine error
! or 1 in case of success.
! Possible causes for a routine errors:
! - nchk = -1: a failure in the module initialization
! - nhck = -2: a wrong H-M symbol
! In case of an error "npatn" returns 0.
!
IMPLICIT NONE
!
CHARACTER(len=*),INTENT(IN)::strname
INTEGER,INTENT(OUT)::npatn,nchk
!
INTEGER::state
INTEGER::nsgnum
!
nchk=0
npatn=0
nsgnum=0
IF (LEN_TRIM(strname)<2) GOTO 802 ! false name error
!
CALL SG_ISREADY(state)
IF (0==state) CALL SG_INIT(state)
IF (0==state) GOTO 801 ! initialization error
!
! Get the index of the space group from its name.
CALL SG_NAMGETNUM(strname,nsgnum)
IF (nsgnum==0) GOTO 802 ! false name error
!
! space group identified, get all the data
npatn = sg_patn(nsgnum)
nchk=1 ! success!
GOTO 1000
!
! error handling
801 CONTINUE ! module initialization problem
nchk = -1
GOTO 900
802 CONTINUE ! invalid space group name
nchk = -2
GOTO 900
!
! error clean up, go here only in case of a routine error!
900 CONTINUE
npatn=0
!
! finish
1000 CONTINUE
!
END SUBROUTINE SG_NAMGETPATN
!
!
!
!
SUBROUTINE SG_NAMGETSYMNUM(strname,nsymnum,nchk)
!
! Returns the number of symmetry operations "nsymnum" for the
! space group H-M symbol "strname".
! The variable "nchk" returns <=0 in case of a routine error
! or 1 in case of success.
! Possible causes for a routine errors:
! - nchk = -1: a failure in the module initialization
! - nhck = -2: a wrong H-M symbol
! In case of an error "nsymnum" returns 0.
!
IMPLICIT NONE
!
CHARACTER(len=*),INTENT(IN)::strname
INTEGER,INTENT(OUT)::nsymnum,nchk
!
INTEGER::state
INTEGER::nsgnum
!
nchk=0
nsymnum=0
nsgnum=0
IF (LEN_TRIM(strname)<2) GOTO 802 ! false name error
!
CALL SG_ISREADY(state)
IF (0==state) CALL SG_INIT(state)
IF (0==state) GOTO 801 ! initialization error
!
! Get the index of the space group from its name.
CALL SG_NAMGETNUM(strname,nsgnum)
IF (nsgnum==0) GOTO 802 ! false name error
!
! space group identified, get all the data
nsymnum = sg_symnum(nsgnum)
nchk=1 ! success!
GOTO 1000
!
! error handling
801 CONTINUE ! module initialization problem
nchk = -1
GOTO 900
802 CONTINUE ! invalid space group name
nchk = -2
GOTO 900
!
! error clean up, go here only in case of a routine error!
900 CONTINUE
nsymnum=0
!
! finish
1000 CONTINUE
!
END SUBROUTINE SG_NAMGETSYMNUM
!
!
!
!
SUBROUTINE SG_NAMGETSYMOP(strname,isymop,strsymop,nchk)
!
! Returns the symmetry operation string number "isymop" for the
! space group H-M symbol "strname" in the string "strsymop".
! The variable "nchk" returns <=0 in case of a routine error
! or 1 in case of success.
! Possible causes for a routine errors:
! - nchk = -1: a failure in the module initialization
! - nhck = -2: a wrong H-M symbol
! - nhck = -3: an invalid symmetry operation index
! The string "strsymop" MUST BE of the form
!     CHARACTER(len=sg_soplen)::strsymop
!
IMPLICIT NONE
!
CHARACTER(len=*),INTENT(IN)::strname
INTEGER,INTENT(IN)::isymop
CHARACTER(len=sg_soplen),DIMENSION(:),INTENT(OUT)::strsymop
INTEGER,INTENT(OUT)::nchk
!
INTEGER::state
INTEGER::nsgnum
INTEGER::nsymnum
!
nchk=0
nsymnum=0
nsgnum=0
strsymop=''
IF (LEN_TRIM(strname)<2) GOTO 802 ! false name error
!
CALL SG_ISREADY(state)
IF (0==state) CALL SG_INIT(state)
IF (0==state) GOTO 801 ! initialization error
!
! Get the index of the space group from its name.
CALL SG_NAMGETNUM(strname,nsgnum)
IF (nsgnum==0) GOTO 802 ! false name error
!
! space group identified, get all the data
nsymnum = sg_symnum(nsgnum)
IF (isymop<1 .OR. isymop>nsymnum) GOTO 803 ! invalid operation index.
! all is ok now.
strsymop = sg_symmetry(isymop,nsgnum) ! get the string.
nchk=1 ! success!
GOTO 1000
!
! error handling
801 CONTINUE ! module initialization problem
nchk = -1
GOTO 900
802 CONTINUE ! invalid space group name
nchk = -2
GOTO 900
803 CONTINUE ! invalid index of symmetry operation
nchk = -3
GOTO 900
!
! error clean up, go here only in case of a routine error!
900 CONTINUE
strsymop=''
!
! finish
1000 CONTINUE
!
END SUBROUTINE SG_NAMGETSYMOP
!
!
!
!
SUBROUTINE SG_NAMGETSYMOPS(strname,strsymops,nsymnum,nchk)
!
! Returns the list of all symmetry operation strings for the
! space group H-M symbol "strname" in the array "strsymops".
! The number of symmetry operation strings is returned in
! the variable nsymnum, which should be identical to
! SIZE(strsymops) on the routine exit in case of success.
! The variable "nchk" returns <=0 in case of a routine error
! or 1 in case of success.
! Possible causes for a routine errors:
! - nchk =  0: an unknown error
! - nchk = -1: a failure in the module initialization
! - nhck = -2: a wrong H-M symbol
! The array "strsymops" MUST BE an allocatable array of the
! form CHARACTER(len=sg_soplen),ALLOCATABLE::strsymops(:)
! The array will be reallocated by the routine. Or returns
! as unallocated in case of a failure.
!
IMPLICIT NONE
!
CHARACTER(len=*),INTENT(IN)::strname
CHARACTER(len=sg_soplen),DIMENSION(:),ALLOCATABLE::strsymops
INTEGER,INTENT(OUT)::nsymnum,nchk
!
INTEGER::nalloc
INTEGER::state
INTEGER::nsgnum
!
nchk=0
nalloc=0
nsymnum=0
nsgnum=0
IF (ALLOCATED(strsymops)) DEALLOCATE(strsymops,stat=nalloc)
IF (LEN_TRIM(strname)<2) GOTO 802 ! false name error
!
CALL SG_ISREADY(state)
IF (0==state) CALL SG_INIT(state)
IF (0==state) GOTO 801 ! initialization error
!
! Get the index of the space group from its name.
CALL SG_NAMGETNUM(strname,nsgnum)
IF (nsgnum==0) GOTO 802 ! false name error
!
! space group identified, get all the data
nsymnum = sg_symnum(nsgnum)
ALLOCATE( strsymops(1:nsymnum), stat=nalloc )
IF (nalloc.NE.0) GOTO 800 ! Again, this should not happen.
strsymops='' ! init.
strsymops(1:nsymnum) = sg_symmetry(1:nsymnum,nsgnum) ! copy.
nchk=1 ! got it, success!
GOTO 1000
!
! error handling
800 CONTINUE ! unknown error
nchk = 0
GOTO 900
801 CONTINUE ! module initialization problem
nchk = -1
GOTO 900
802 CONTINUE ! invalid space group name
nchk = -2
GOTO 900
!
! error clean up, go here only in case of a routine error!
900 CONTINUE
IF (ALLOCATED(strsymops)) DEALLOCATE(strsymops,stat=nalloc)
nsymnum=0
!
! finish
1000 CONTINUE
!
END SUBROUTINE SG_NAMGETSYMOPS
!
!
!
!
SUBROUTINE SG_ISREADY(state)
!
! Reports the module initialization status.
! state == 0: module is not initialized.
! state == 1: module is initialized.
!
IMPLICIT NONE
!
INTEGER,INTENT(OUT)::state
! Initialize with "initialized" by default.
state=1
! Check if all arrays for their allocation state. If one allocation
! is missing, switch to "not initialized" ...
IF(.NOT.ALLOCATED(sg_name)) state=0
IF(.NOT.ALLOCATED(sg_symmetry)) state=0
IF(.NOT.ALLOCATED(sg_patn)) state=0
IF(.NOT.ALLOCATED(sg_symnum)) state=0
!
END SUBROUTINE SG_ISREADY
!
!
!
!
SUBROUTINE SG_UNINIT()
!
! De-initializes the module by de-allocating the data arrays.
! 
IMPLICIT NONE
!
INTEGER::nalloc ! You may use nalloc to post de-allocation error messages.
!               ! .oO( Although this usually never happens. )
!
nalloc=0
IF(ALLOCATED(sg_name)) DEALLOCATE(sg_name,stat=nalloc)
IF(ALLOCATED(sg_symmetry)) DEALLOCATE(sg_symmetry,stat=nalloc)
IF(ALLOCATED(sg_patn)) DEALLOCATE(sg_patn,stat=nalloc)
IF(ALLOCATED(sg_symnum)) DEALLOCATE(sg_symnum,stat=nalloc)
!
END SUBROUTINE SG_UNINIT
!
!
!
!
!
SUBROUTINE SG_INIT(state)
!
! Initializes the module and reports the final initialization status.
! state == 0: module is not initialized = Failure of the routine.
! state == 1: module is initialized     = Success of the routine.
!
IMPLICIT NONE
!
INTEGER,INTENT(OUT)::state
INTEGER::nalloc ! You may use nalloc to post allocation error messages.
!               ! .oO( Although this usually never happens nowadays. )
!
! Get the current initialization status
CALL SG_UNINIT() ! Deinitialize the module, deallocate all arrays
state=0
nalloc=0
!
! Allocations
ALLOCATE(sg_name(1:sg_nummax), stat=nalloc)
IF (nalloc.NE.0) GOTO 800
ALLOCATE(sg_symmetry(1:sg_opmax,1:sg_nummax), stat=nalloc)
IF (nalloc.NE.0) GOTO 800
ALLOCATE(sg_patn(1:sg_nummax), stat=nalloc)
IF (nalloc.NE.0) GOTO 800
ALLOCATE(sg_symnum(1:sg_nummax), stat=nalloc)
IF (nalloc.NE.0) GOTO 800
!
! Pre-Initializations
sg_name = ''
sg_symmetry = ''
sg_patn = 0
sg_symnum = 0
!
! Data initializations (There's a very long list now.)
sg_name(1) = 'P1'
sg_patn(1) = 2
sg_symnum(1) = 1
sg_symmetry(1, 1) = 'x,y,z'
sg_name(2) = 'P-1'
sg_patn(2) = 2
sg_symnum(2) = 2
sg_symmetry(1, 2) = 'x,y,z'
sg_symmetry(2, 2) = '-x,-y,-z'
sg_name(3) = 'P2'
sg_patn(3) = 10
sg_symnum(3) = 2
sg_symmetry(1, 3) = 'x,y,z'
sg_symmetry(2, 3) = '-x,y,-z'
sg_name(4) = 'P2(1)'
sg_patn(4) = 10
sg_symnum(4) = 2
sg_symmetry(1, 4) = 'x,y,z'
sg_symmetry(2, 4) = '-x,y+1/2,-z'
sg_name(5) = 'C2'
sg_patn(5) = 12
sg_symnum(5) = 4
sg_symmetry(1, 5) = 'x,y,z'
sg_symmetry(2, 5) = '-x,y,-z'
sg_symmetry(3, 5) = '1/2+x,1/2+y,z'
sg_symmetry(4, 5) = '1/2-x,1/2+y,-z'
sg_name(6) = 'Pm'
sg_patn(6) = 10
sg_symnum(6) = 2
sg_symmetry(1, 6) = 'x,y,z'
sg_symmetry(2, 6) = 'x,-y,z'
sg_name(7) = 'Pc'
sg_patn(7) = 10
sg_symnum(7) = 2
sg_symmetry(1, 7) = 'x,y,z'
sg_symmetry(2, 7) = 'x,-y,1/2+z'
sg_name(8) = 'Cm'
sg_patn(8) = 12
sg_symnum(8) = 4
sg_symmetry(1, 8) = 'x,y,z'
sg_symmetry(2, 8) = 'x,-y,z'
sg_symmetry(3, 8) = '1/2+x,1/2+y,z'
sg_symmetry(4, 8) = '1/2+x,1/2-y,z'
sg_name(9) = 'Cc'
sg_patn(9) = 12
sg_symnum(9) = 4
sg_symmetry(1, 9) = 'x,y,z'
sg_symmetry(2, 9) = 'x,-y,1/2+z'
sg_symmetry(3, 9) = '1/2+x,1/2+y,z'
sg_symmetry(4, 9) = '1/2+x,1/2-y,1/2+z'
sg_name(10) = 'P2/m'
sg_patn(10) = 10
sg_symnum(10) = 4
sg_symmetry(1, 10) = 'x,y,z'
sg_symmetry(2, 10) = 'x,-y,z'
sg_symmetry(3, 10) = '-x,y,-z'
sg_symmetry(4, 10) = '-x,-y,-z'
sg_name(11) = 'P2(1)/m'
sg_patn(11) = 10
sg_symnum(11) = 4
sg_symmetry(1, 11) = 'x,y,z'
sg_symmetry(2, 11) = '-x,1/2+y,-z'
sg_symmetry(3, 11) = '-x,-y,-z'
sg_symmetry(4, 11) = 'x,1/2-y,z'
sg_name(12) = 'C2/m'
sg_patn(12) = 12
sg_symnum(12) = 8
sg_symmetry(1, 12) = 'x,y,z'
sg_symmetry(2, 12) = 'x,-y,z'
sg_symmetry(3, 12) = '-x,y,-z'
sg_symmetry(4, 12) = '-x,-y,-z'
sg_symmetry(5, 12) = '1/2+x,1/2+y,z'
sg_symmetry(6, 12) = '1/2+x,1/2-y,z'
sg_symmetry(7, 12) = '1/2-x,1/2+y,-z'
sg_symmetry(8, 12) = '1/2-x,1/2-y,-z'
sg_name(13) = 'P2/c'
sg_patn(13) = 10
sg_symnum(13) = 4
sg_symmetry(1, 13) = 'x,y,z'
sg_symmetry(2, 13) = '-x,y,1/2-z'
sg_symmetry(3, 13) = '-x,-y,-z'
sg_symmetry(4, 13) = 'x,-y,1/2+z'
sg_name(14) = 'P2(1)/c'
sg_patn(14) = 10
sg_symnum(14) = 4
sg_symmetry(1, 14) = 'x,y,z'
sg_symmetry(2, 14) = '-x,-y,-z'
sg_symmetry(3, 14) = '-x,1/2+y,1/2-z'
sg_symmetry(4, 14) = 'x,1/2-y,1/2+z'
sg_name(15) = 'C2/c'
sg_patn(15) = 12
sg_symnum(15) = 8
sg_symmetry(1, 15) = 'x,y,z'
sg_symmetry(2, 15) = '-x,y,1/2-z'
sg_symmetry(3, 15) = '-x,-y,-z'
sg_symmetry(4, 15) = 'x,-y,1/2+z'
sg_symmetry(5, 15) = '1/2+x,1/2+y,z'
sg_symmetry(6, 15) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(7, 15) = '1/2-x,1/2-y,-z'
sg_symmetry(8, 15) = '1/2+x,1/2-y,1/2+z'
sg_name(16) = 'P222'
sg_patn(16) = 47
sg_symnum(16) = 4
sg_symmetry(1, 16) = 'x,y,z'
sg_symmetry(2, 16) = '-x,-y,z'
sg_symmetry(3, 16) = '-x,y,-z'
sg_symmetry(4, 16) = 'x,-y,-z'
sg_name(17) = 'P222(1)'
sg_patn(17) = 47
sg_symnum(17) = 4
sg_symmetry(1, 17) = 'x,y,z'
sg_symmetry(2, 17) = '-x,-y,1/2+z'
sg_symmetry(3, 17) = '-x,y,1/2-z'
sg_symmetry(4, 17) = 'x,-y,-z'
sg_name(18) = 'P2(1)2(1)2'
sg_patn(18) = 47
sg_symnum(18) = 4
sg_symmetry(1, 18) = 'x,y,z'
sg_symmetry(2, 18) = '-x,-y,z'
sg_symmetry(3, 18) = '1/2-x,1/2+y,-z'
sg_symmetry(4, 18) = '1/2+x,1/2-y,-z'
sg_name(19) = 'P2(1)2(1)2(1)'
sg_patn(19) = 47
sg_symnum(19) = 4
sg_symmetry(1, 19) = 'x,y,z'
sg_symmetry(2, 19) = '1/2-x,-y,1/2+z'
sg_symmetry(3, 19) = '-x,1/2+y,1/2-z'
sg_symmetry(4, 19) = '1/2+x,1/2-y,-z'
sg_name(20) = 'C222(1)'
sg_patn(20) = 65
sg_symnum(20) = 8
sg_symmetry(1, 20) = 'x,y,z'
sg_symmetry(2, 20) = '-x,-y,1/2+z'
sg_symmetry(3, 20) = '-x,y,1/2-z'
sg_symmetry(4, 20) = 'x,-y,-z'
sg_symmetry(5, 20) = '1/2+x,1/2+y,z'
sg_symmetry(6, 20) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(7, 20) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(8, 20) = '1/2+x,1/2-y,-z'
sg_name(21) = 'C222'
sg_patn(21) = 65
sg_symnum(21) = 8
sg_symmetry(1, 21) = 'x,y,z'
sg_symmetry(2, 21) = '-x,-y,z'
sg_symmetry(3, 21) = '-x,y,-z'
sg_symmetry(4, 21) = 'x,-y,-z'
sg_symmetry(5, 21) = '1/2+x,1/2+y,z'
sg_symmetry(6, 21) = '1/2-x,1/2-y,z'
sg_symmetry(7, 21) = '1/2-x,1/2+y,-z'
sg_symmetry(8, 21) = '1/2+x,1/2-y,-z'
sg_name(22) = 'F222'
sg_patn(22) = 69
sg_symnum(22) = 16
sg_symmetry(1, 22) = 'x,y,z'
sg_symmetry(2, 22) = '-x,-y,z'
sg_symmetry(3, 22) = '-x,y,-z'
sg_symmetry(4, 22) = 'x,-y,-z'
sg_symmetry(5, 22) = 'x,1/2+y,1/2+z'
sg_symmetry(6, 22) = '-x,1/2-y,1/2+z'
sg_symmetry(7, 22) = '-x,1/2+y,1/2-z'
sg_symmetry(8, 22) = 'x,1/2-y,1/2-z'
sg_symmetry(9, 22) = '1/2+x,y,1/2+z'
sg_symmetry(10, 22) = '1/2-x,-y,1/2+z'
sg_symmetry(11, 22) = '1/2-x,y,1/2-z'
sg_symmetry(12, 22) = '1/2+x,-y,1/2-z'
sg_symmetry(13, 22) = '1/2+x,1/2+y,z'
sg_symmetry(14, 22) = '1/2-x,1/2-y,z'
sg_symmetry(15, 22) = '1/2-x,1/2+y,-z'
sg_symmetry(16, 22) = '1/2+x,1/2-y,-z'
sg_name(23) = 'I222'
sg_patn(23) = 71
sg_symnum(23) = 8
sg_symmetry(1, 23) = 'x,y,z'
sg_symmetry(2, 23) = '-x,-y,z'
sg_symmetry(3, 23) = 'x,-y,-z'
sg_symmetry(4, 23) = '-x,y,-z'
sg_symmetry(5, 23) = 'x+1/2,y+1/2,z+1/2'
sg_symmetry(6, 23) = '-x+1/2,-y+1/2,z+1/2'
sg_symmetry(7, 23) = 'x+1/2,-y+1/2,-z+1/2'
sg_symmetry(8, 23) = '-x+1/2,y+1/2,-z+1/2'
sg_name(24) = 'I2(1)2(1)2(1)'
sg_patn(24) = 71
sg_symnum(24) = 8
sg_symmetry(1, 24) = 'x,y,z'
sg_symmetry(2, 24) = '1/2-x,-y,1/2+z'
sg_symmetry(3, 24) = '-x,1/2+y,1/2-z'
sg_symmetry(4, 24) = '1/2+x,1/2-y,-z'
sg_symmetry(5, 24) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(6, 24) = '-x,1/2-y,z'
sg_symmetry(7, 24) = '1/2-x,y,-z'
sg_symmetry(8, 24) = 'x,-y,1/2-z'
sg_name(25) = 'Pmm2'
sg_patn(25) = 47
sg_symnum(25) = 4
sg_symmetry(1, 25) = 'x,y,z'
sg_symmetry(2, 25) = '-x,-y,z'
sg_symmetry(3, 25) = 'x,-y,z'
sg_symmetry(4, 25) = '-x,y,z'
sg_name(26) = 'Pmc2(1)'
sg_patn(26) = 47
sg_symnum(26) = 4
sg_symmetry(1, 26) = 'x,y,z'
sg_symmetry(2, 26) = '-x,-y,1/2+z'
sg_symmetry(3, 26) = 'x,-y,1/2+z'
sg_symmetry(4, 26) = '-x,y,z'
sg_name(27) = 'Pcc2'
sg_patn(27) = 47
sg_symnum(27) = 4
sg_symmetry(1, 27) = 'x,y,z'
sg_symmetry(2, 27) = '-x,-y,z'
sg_symmetry(3, 27) = 'x,-y,1/2+z'
sg_symmetry(4, 27) = '-x,y,1/2+z'
sg_name(28) = 'Pma2'
sg_patn(28) = 47
sg_symnum(28) = 4
sg_symmetry(1, 28) = 'x,y,z'
sg_symmetry(2, 28) = '-x,-y,z'
sg_symmetry(3, 28) = '1/2+x,-y,z'
sg_symmetry(4, 28) = '1/2-x,y,z'
sg_name(29) = 'Pca2(1)'
sg_patn(29) = 47
sg_symnum(29) = 4
sg_symmetry(1, 29) = 'x,y,z'
sg_symmetry(2, 29) = '-x,-y,1/2+z'
sg_symmetry(3, 29) = '1/2+x,-y,z'
sg_symmetry(4, 29) = '1/2-x,y,1/2+z'
sg_name(30) = 'Pnc2'
sg_patn(30) = 47
sg_symnum(30) = 4
sg_symmetry(1, 30) = 'x,y,z'
sg_symmetry(2, 30) = '-x,-y,z'
sg_symmetry(3, 30) = 'x,1/2-y,1/2+z'
sg_symmetry(4, 30) = '-x,1/2+y,1/2+z'
sg_name(31) = 'Pmn2(1)'
sg_patn(31) = 47
sg_symnum(31) = 4
sg_symmetry(1, 31) = 'x,y,z'
sg_symmetry(2, 31) = '1/2-x,-y,1/2+z'
sg_symmetry(3, 31) = '1/2+x,-y,1/2+z'
sg_symmetry(4, 31) = '-x,y,z'
sg_name(32) = 'Pba2'
sg_patn(32) = 47
sg_symnum(32) = 4
sg_symmetry(1, 32) = 'x,y,z'
sg_symmetry(2, 32) = '-x,-y,z'
sg_symmetry(3, 32) = '1/2+x,1/2-y,z'
sg_symmetry(4, 32) = '1/2-x,1/2+y,z'
sg_name(33) = 'Pna2(1)'
sg_patn(33) = 47
sg_symnum(33) = 4
sg_symmetry(1, 33) = 'x,y,z'
sg_symmetry(2, 33) = '-x,-y,1/2+z'
sg_symmetry(3, 33) = '1/2+x,1/2-y,z'
sg_symmetry(4, 33) = '1/2-x,1/2+y,1/2+z'
sg_name(34) = 'Pnn2'
sg_patn(34) = 47
sg_symnum(34) = 4
sg_symmetry(1, 34) = 'x,y,z'
sg_symmetry(2, 34) = '-x,-y,z'
sg_symmetry(3, 34) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(4, 34) = '1/2-x,1/2+y,1/2+z'
sg_name(35) = 'Cmm2'
sg_patn(35) = 65
sg_symnum(35) = 8
sg_symmetry(1, 35) = 'x,y,z'
sg_symmetry(2, 35) = '-x,-y,z'
sg_symmetry(3, 35) = 'x,-y,z'
sg_symmetry(4, 35) = '-x,y,z'
sg_symmetry(5, 35) = '1/2+x,1/2+y,z'
sg_symmetry(6, 35) = '1/2-x,1/2-y,z'
sg_symmetry(7, 35) = '1/2+x,1/2-y,z'
sg_symmetry(8, 35) = '1/2-x,1/2+y,z'
sg_name(36) = 'Cmc2(1)'
sg_patn(36) = 65
sg_symnum(36) = 8
sg_symmetry(1, 36) = 'x,y,z'
sg_symmetry(2, 36) = '-x,-y,1/2+z'
sg_symmetry(3, 36) = 'x,-y,1/2+z'
sg_symmetry(4, 36) = '-x,y,z'
sg_symmetry(5, 36) = '1/2+x,1/2+y,z'
sg_symmetry(6, 36) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(7, 36) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(8, 36) = '1/2-x,1/2+y,z'
sg_name(37) = 'Ccc2'
sg_patn(37) = 65
sg_symnum(37) = 8
sg_symmetry(1, 37) = 'x,y,z'
sg_symmetry(2, 37) = '-x,-y,z'
sg_symmetry(3, 37) = 'x,-y,1/2+z'
sg_symmetry(4, 37) = '-x,y,1/2+z'
sg_symmetry(5, 37) = '1/2+x,1/2+y,z'
sg_symmetry(6, 37) = '1/2-x,1/2-y,z'
sg_symmetry(7, 37) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(8, 37) = '1/2-x,1/2+y,1/2+z'
sg_name(38) = 'Amm2'
sg_patn(38) = 65
sg_symnum(38) = 8
sg_symmetry(1, 38) = 'x,y,z'
sg_symmetry(2, 38) = '-x,-y,z'
sg_symmetry(3, 38) = 'x,-y,z'
sg_symmetry(4, 38) = '-x,y,z'
sg_symmetry(5, 38) = 'x,1/2+y,1/2+z'
sg_symmetry(6, 38) = '-x,1/2-y,1/2+z'
sg_symmetry(7, 38) = 'x,1/2-y,1/2+z'
sg_symmetry(8, 38) = '-x,1/2+y,1/2+z'
sg_name(39) = 'Abm2'
sg_patn(39) = 65
sg_symnum(39) = 8
sg_symmetry(1, 39) = 'x,y,z'
sg_symmetry(2, 39) = '-x,-y,z'
sg_symmetry(3, 39) = 'x,1/2-y,z'
sg_symmetry(4, 39) = '-x,1/2+y,z'
sg_symmetry(5, 39) = 'x,1/2+y,1/2+z'
sg_symmetry(6, 39) = '-x,1/2-y,1/2+z'
sg_symmetry(7, 39) = 'x,-y,1/2+z'
sg_symmetry(8, 39) = '-x,y,1/2+z'
sg_name(40) = 'Ama2'
sg_patn(40) = 65
sg_symnum(40) = 8
sg_symmetry(1, 40) = 'x,y,z'
sg_symmetry(2, 40) = '-x,-y,z'
sg_symmetry(3, 40) = '1/2+x,-y,z'
sg_symmetry(4, 40) = '1/2-x,y,z'
sg_symmetry(5, 40) = 'x,1/2+y,1/2+z'
sg_symmetry(6, 40) = '-x,1/2-y,1/2+z'
sg_symmetry(7, 40) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(8, 40) = '1/2-x,1/2+y,1/2+z'
sg_name(41) = 'Aba2'
sg_patn(41) = 65
sg_symnum(41) = 8
sg_symmetry(1, 41) = 'x,y,z'
sg_symmetry(2, 41) = '-x,-y,z'
sg_symmetry(3, 41) = '1/2+x,1/2-y,z'
sg_symmetry(4, 41) = '1/2-x,1/2+y,z'
sg_symmetry(5, 41) = 'x,1/2+y,1/2+z'
sg_symmetry(6, 41) = '-x,1/2-y,1/2+z'
sg_symmetry(7, 41) = '1/2+x,-y,1/2+z'
sg_symmetry(8, 41) = '1/2-x,y,1/2+z'
sg_name(42) = 'Fmm2'
sg_patn(42) = 69
sg_symnum(42) = 16
sg_symmetry(1, 42) = 'x,y,z'
sg_symmetry(2, 42) = '-x,-y,z'
sg_symmetry(3, 42) = 'x,-y,z'
sg_symmetry(4, 42) = '-x,y,z'
sg_symmetry(5, 42) = 'x,1/2+y,1/2+z'
sg_symmetry(6, 42) = '-x,1/2-y,1/2+z'
sg_symmetry(7, 42) = 'x,1/2-y,1/2+z'
sg_symmetry(8, 42) = '-x,1/2+y,1/2+z'
sg_symmetry(9, 42) = '1/2+x,y,1/2+z'
sg_symmetry(10, 42) = '1/2-x,-y,1/2+z'
sg_symmetry(11, 42) = '1/2+x,-y,1/2+z'
sg_symmetry(12, 42) = '1/2-x,y,1/2+z'
sg_symmetry(13, 42) = '1/2+x,1/2+y,z'
sg_symmetry(14, 42) = '1/2-x,1/2-y,z'
sg_symmetry(15, 42) = '1/2+x,1/2-y,z'
sg_symmetry(16, 42) = '1/2-x,1/2+y,z'
sg_name(43) = 'Fdd2'
sg_patn(43) = 69
sg_symnum(43) = 16
sg_symmetry(1, 43) = 'x,y,z'
sg_symmetry(2, 43) = '-x,-y,z'
sg_symmetry(3, 43) = '1/4+x,1/4-y,1/4+z'
sg_symmetry(4, 43) = '1/4-x,1/4+y,1/4+z'
sg_symmetry(5, 43) = 'x,1/2+y,1/2+z'
sg_symmetry(6, 43) = '-x,1/2-y,1/2+z'
sg_symmetry(7, 43) = '1/4+x,3/4-y,3/4+z'
sg_symmetry(8, 43) = '1/4-x,3/4+y,3/4+z'
sg_symmetry(9, 43) = '1/2+x,y,1/2+z'
sg_symmetry(10, 43) = '1/2-x,-y,1/2+z'
sg_symmetry(11, 43) = '3/4+x,1/4-y,3/4+z'
sg_symmetry(12, 43) = '3/4-x,1/4+y,3/4+z'
sg_symmetry(13, 43) = '1/2+x,1/2+y,z'
sg_symmetry(14, 43) = '1/2-x,1/2-y,z'
sg_symmetry(15, 43) = '3/4+x,3/4-y,1/4+z'
sg_symmetry(16, 43) = '3/4-x,3/4+y,1/4+z'
sg_name(44) = 'Imm2'
sg_patn(44) = 71
sg_symnum(44) = 8
sg_symmetry(1, 44) = 'x,y,z'
sg_symmetry(2, 44) = '-x,-y,z'
sg_symmetry(3, 44) = 'x,-y,z'
sg_symmetry(4, 44) = '-x,y,z'
sg_symmetry(5, 44) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(6, 44) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(7, 44) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(8, 44) = '1/2-x,1/2+y,1/2+z'
sg_name(45) = 'Iba2'
sg_patn(45) = 71
sg_symnum(45) = 8
sg_symmetry(1, 45) = 'x,y,z'
sg_symmetry(2, 45) = '-x,-y,z'
sg_symmetry(3, 45) = '1/2+x,1/2-y,z'
sg_symmetry(4, 45) = '1/2-x,1/2+y,z'
sg_symmetry(5, 45) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(6, 45) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(7, 45) = 'x,-y,1/2+z'
sg_symmetry(8, 45) = '-x,y,1/2+z'
sg_name(46) = 'Ima2'
sg_patn(46) = 71
sg_symnum(46) = 8
sg_symmetry(1, 46) = 'x,y,z'
sg_symmetry(2, 46) = '-x,-y,z'
sg_symmetry(3, 46) = '1/2+x,-y,z'
sg_symmetry(4, 46) = '1/2-x,y,z'
sg_symmetry(5, 46) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(6, 46) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(7, 46) = 'x,1/2-y,1/2+z'
sg_symmetry(8, 46) = '-x,1/2+y,1/2+z'
sg_name(47) = 'Pmmm'
sg_patn(47) = 47
sg_symnum(47) = 8
sg_symmetry(1, 47) = 'x,y,z'
sg_symmetry(2, 47) = '-x,-y,z'
sg_symmetry(3, 47) = '-x,y,-z'
sg_symmetry(4, 47) = 'x,-y,-z'
sg_symmetry(5, 47) = '-x,-y,-z'
sg_symmetry(6, 47) = 'x,y,-z'
sg_symmetry(7, 47) = 'x,-y,z'
sg_symmetry(8, 47) = '-x,y,z'
sg_name(48) = 'Pnnn'
sg_patn(48) = 47
sg_symnum(48) = 8
sg_symmetry(1, 48) = 'x,y,z'
sg_symmetry(2, 48) = '-x,-y,z'
sg_symmetry(3, 48) = '-x,y,-z'
sg_symmetry(4, 48) = 'x,-y,-z'
sg_symmetry(5, 48) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(6, 48) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(7, 48) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(8, 48) = '1/2-x,1/2+y,1/2+z'
sg_name(49) = 'Pccm'
sg_patn(49) = 47
sg_symnum(49) = 8
sg_symmetry(1, 49) = 'x,y,z'
sg_symmetry(2, 49) = '-x,-y,z'
sg_symmetry(3, 49) = '-x,y,1/2-z'
sg_symmetry(4, 49) = 'x,-y,1/2-z'
sg_symmetry(5, 49) = '-x,-y,-z'
sg_symmetry(6, 49) = 'x,y,-z'
sg_symmetry(7, 49) = 'x,-y,1/2+z'
sg_symmetry(8, 49) = '-x,y,1/2+z'
sg_name(50) = 'Pban'
sg_patn(50) = 47
sg_symnum(50) = 8
sg_symmetry(1, 50) = 'x,y,z'
sg_symmetry(2, 50) = '-x,-y,z'
sg_symmetry(3, 50) = '-x,y,-z'
sg_symmetry(4, 50) = 'x,-y,-z'
sg_symmetry(5, 50) = '1/2-x,1/2-y,-z'
sg_symmetry(6, 50) = '1/2+x,1/2+y,-z'
sg_symmetry(7, 50) = '1/2+x,1/2-y,z'
sg_symmetry(8, 50) = '1/2-x,1/2+y,z'
sg_name(51) = 'Pmma'
sg_patn(51) = 47
sg_symnum(51) = 8
sg_symmetry(1, 51) = 'x,y,z'
sg_symmetry(2, 51) = '1/2-x,-y,z'
sg_symmetry(3, 51) = '-x,y,-z'
sg_symmetry(4, 51) = '1/2+x,-y,-z'
sg_symmetry(5, 51) = '-x,-y,-z'
sg_symmetry(6, 51) = '1/2+x,y,-z'
sg_symmetry(7, 51) = 'x,-y,z'
sg_symmetry(8, 51) = '1/2-x,y,z'
sg_name(52) = 'Pnna'
sg_patn(52) = 47
sg_symnum(52) = 8
sg_symmetry(1, 52) = 'x,y,z'
sg_symmetry(2, 52) = '1/2-x,-y,z'
sg_symmetry(3, 52) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(4, 52) = 'x,1/2-y,1/2-z'
sg_symmetry(5, 52) = '-x,-y,-z'
sg_symmetry(6, 52) = '1/2+x,y,-z'
sg_symmetry(7, 52) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(8, 52) = '-x,1/2+y,1/2+z'
sg_name(53) = 'Pmna'
sg_patn(53) = 47
sg_symnum(53) = 8
sg_symmetry(1, 53) = 'x,y,z'
sg_symmetry(2, 53) = '1/2-x,-y,1/2+z'
sg_symmetry(3, 53) = '1/2-x,y,1/2-z'
sg_symmetry(4, 53) = 'x,-y,-z'
sg_symmetry(5, 53) = '-x,-y,-z'
sg_symmetry(6, 53) = '1/2+x,y,1/2-z'
sg_symmetry(7, 53) = '1/2+x,-y,1/2+z'
sg_symmetry(8, 53) = '-x,y,z'
sg_name(54) = 'Pcca'
sg_patn(54) = 47
sg_symnum(54) = 8
sg_symmetry(1, 54) = 'x,y,z'
sg_symmetry(2, 54) = '1/2-x,-y,z'
sg_symmetry(3, 54) = '-x,y,1/2-z'
sg_symmetry(4, 54) = '1/2+x,-y,1/2-z'
sg_symmetry(5, 54) = '-x,-y,-z'
sg_symmetry(6, 54) = '1/2+x,y,-z'
sg_symmetry(7, 54) = 'x,-y,1/2+z'
sg_symmetry(8, 54) = '1/2-x,y,1/2+z'
sg_name(55) = 'Pbam'
sg_patn(55) = 47
sg_symnum(55) = 8
sg_symmetry(1, 55) = 'x,y,z'
sg_symmetry(2, 55) = '-x,-y,z'
sg_symmetry(3, 55) = '1/2-x,1/2+y,-z'
sg_symmetry(4, 55) = '1/2+x,1/2-y,-z'
sg_symmetry(5, 55) = '-x,-y,-z'
sg_symmetry(6, 55) = 'x,y,-z'
sg_symmetry(7, 55) = '1/2+x,1/2-y,z'
sg_symmetry(8, 55) = '1/2-x,1/2+y,z'
sg_name(56) = 'Pccn'
sg_patn(56) = 47
sg_symnum(56) = 8
sg_symmetry(1, 56) = 'x,y,z'
sg_symmetry(2, 56) = '1/2-x,1/2-y,z'
sg_symmetry(3, 56) = '-x,1/2+y,1/2-z'
sg_symmetry(4, 56) = '1/2+x,-y,1/2-z'
sg_symmetry(5, 56) = '-x,-y,-z'
sg_symmetry(6, 56) = '1/2+x,1/2+y,-z'
sg_symmetry(7, 56) = 'x,1/2-y,1/2+z'
sg_symmetry(8, 56) = '1/2-x,y,1/2+z'
sg_name(57) = 'Pbcm'
sg_patn(57) = 47
sg_symnum(57) = 8
sg_symmetry(1, 57) = 'x,y,z'
sg_symmetry(2, 57) = '-x,-y,1/2+z'
sg_symmetry(3, 57) = '-x,1/2+y,1/2-z'
sg_symmetry(4, 57) = 'x,1/2-y,-z'
sg_symmetry(5, 57) = '-x,-y,-z'
sg_symmetry(6, 57) = 'x,y,1/2-z'
sg_symmetry(7, 57) = 'x,1/2-y,1/2+z'
sg_symmetry(8, 57) = '-x,1/2+y,z'
sg_name(58) = 'Pnnm'
sg_patn(58) = 47
sg_symnum(58) = 8
sg_symmetry(1, 58) = 'x,y,z'
sg_symmetry(2, 58) = '-x,-y,z'
sg_symmetry(3, 58) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(4, 58) = '1/2+x,1/2-y,1/2-z'
sg_symmetry(5, 58) = '-x,-y,-z'
sg_symmetry(6, 58) = 'x,y,-z'
sg_symmetry(7, 58) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(8, 58) = '1/2-x,1/2+y,1/2+z'
sg_name(59) = 'Pmmn'
sg_patn(59) = 47
sg_symnum(59) = 8
sg_symmetry(1, 59) = 'x,y,z'
sg_symmetry(2, 59) = '-x,-y,z'
sg_symmetry(3, 59) = '1/2-x,y+1/2,-z'
sg_symmetry(4, 59) = 'x+1/2,1/2-y,-z'
sg_symmetry(5, 59) = '1/2-x,1/2-y,-z'
sg_symmetry(6, 59) = 'x+1/2,y+1/2,-z'
sg_symmetry(7, 59) = 'x,-y,z'
sg_symmetry(8, 59) = '-x,y,z'
sg_name(60) = 'Pbcn'
sg_patn(60) = 47
sg_symnum(60) = 8
sg_symmetry(1, 60) = 'x,y,z'
sg_symmetry(2, 60) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(3, 60) = '-x,y,1/2-z'
sg_symmetry(4, 60) = '1/2+x,1/2-y,-z'
sg_symmetry(5, 60) = '-x,-y,-z'
sg_symmetry(6, 60) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(7, 60) = 'x,-y,1/2+z'
sg_symmetry(8, 60) = '1/2-x,1/2+y,z'
sg_name(61) = 'Pbca'
sg_patn(61) = 47
sg_symnum(61) = 8
sg_symmetry(1, 61) = 'x,y,z'
sg_symmetry(2, 61) = '1/2-x,-y,1/2+z'
sg_symmetry(3, 61) = '-x,1/2+y,1/2-z'
sg_symmetry(4, 61) = '1/2+x,1/2-y,-z'
sg_symmetry(5, 61) = '-x,-y,-z'
sg_symmetry(6, 61) = '1/2+x,y,1/2-z'
sg_symmetry(7, 61) = 'x,1/2-y,1/2+z'
sg_symmetry(8, 61) = '1/2-x,1/2+y,z'
sg_name(62) = 'Pnma'
sg_patn(62) = 47
sg_symnum(62) = 8
sg_symmetry(1, 62) = 'x,y,z'
sg_symmetry(2, 62) = '-x+1/2,-y,z+1/2'
sg_symmetry(3, 62) = '-x,y+1/2,-z'
sg_symmetry(4, 62) = 'x+1/2,-y+1/2,-z+1/2'
sg_symmetry(5, 62) = '-x,-y,-z'
sg_symmetry(6, 62) = 'x+1/2,y,-z+1/2'
sg_symmetry(7, 62) = 'x,-y+1/2,z'
sg_symmetry(8, 62) = '-x+1/2,y+1/2,z+1/2'
sg_name(63) = 'Cmcm'
sg_patn(63) = 65
sg_symnum(63) = 16
sg_symmetry(1, 63) = 'x,y,z'
sg_symmetry(2, 63) = '-x,-y,1/2+z'
sg_symmetry(3, 63) = '-x,y,1/2-z'
sg_symmetry(4, 63) = 'x,-y,-z'
sg_symmetry(5, 63) = '-x,-y,-z'
sg_symmetry(6, 63) = 'x,y,1/2-z'
sg_symmetry(7, 63) = 'x,-y,1/2+z'
sg_symmetry(8, 63) = '-x,y,z'
sg_symmetry(9, 63) = '1/2+x,1/2+y,z'
sg_symmetry(10, 63) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(11, 63) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(12, 63) = '1/2+x,1/2-y,-z'
sg_symmetry(13, 63) = '1/2-x,1/2-y,-z'
sg_symmetry(14, 63) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(15, 63) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(16, 63) = '1/2-x,1/2+y,z'
sg_name(64) = 'Cmca'
sg_patn(64) = 65
sg_symnum(64) = 16
sg_symmetry(1, 64) = 'x,y,z'
sg_symmetry(2, 64) = '-x,1/2-y,1/2+z'
sg_symmetry(3, 64) = '-x,1/2+y,1/2-z'
sg_symmetry(4, 64) = 'x,-y,-z'
sg_symmetry(5, 64) = '-x,-y,-z'
sg_symmetry(6, 64) = 'x,1/2+y,1/2-z'
sg_symmetry(7, 64) = 'x,1/2-y,1/2+z'
sg_symmetry(8, 64) = '-x,y,z'
sg_symmetry(9, 64) = '1/2+x,1/2+y,z'
sg_symmetry(10, 64) = '1/2-x,-y,1/2+z'
sg_symmetry(11, 64) = '1/2-x,y,1/2-z'
sg_symmetry(12, 64) = '1/2+x,1/2-y,-z'
sg_symmetry(13, 64) = '1/2-x,1/2-y,-z'
sg_symmetry(14, 64) = '1/2+x,y,1/2-z'
sg_symmetry(15, 64) = '1/2+x,-y,1/2+z'
sg_symmetry(16, 64) = '1/2-x,1/2+y,z'
sg_name(65) = 'Cmmm'
sg_patn(65) = 65
sg_symnum(65) = 16
sg_symmetry(1, 65) = 'x,y,z'
sg_symmetry(2, 65) = '-x,-y,z'
sg_symmetry(3, 65) = '-x,y,-z'
sg_symmetry(4, 65) = 'x,-y,-z'
sg_symmetry(5, 65) = '-x,-y,-z'
sg_symmetry(6, 65) = 'x,y,-z'
sg_symmetry(7, 65) = 'x,-y,z'
sg_symmetry(8, 65) = '-x,y,z'
sg_symmetry(9, 65) = '1/2+x,1/2+y,z'
sg_symmetry(10, 65) = '1/2-x,1/2-y,z'
sg_symmetry(11, 65) = '1/2-x,1/2+y,-z'
sg_symmetry(12, 65) = '1/2+x,1/2-y,-z'
sg_symmetry(13, 65) = '1/2-x,1/2-y,-z'
sg_symmetry(14, 65) = '1/2+x,1/2+y,-z'
sg_symmetry(15, 65) = '1/2+x,1/2-y,z'
sg_symmetry(16, 65) = '1/2-x,1/2+y,z'
sg_name(66) = 'Cccm'
sg_patn(66) = 65
sg_symnum(66) = 16
sg_symmetry(1, 66) = 'x,y,z'
sg_symmetry(2, 66) = '-x,-y,z'
sg_symmetry(3, 66) = '-x,y,1/2-z'
sg_symmetry(4, 66) = 'x,-y,1/2-z'
sg_symmetry(5, 66) = '-x,-y,-z'
sg_symmetry(6, 66) = 'x,y,-z'
sg_symmetry(7, 66) = 'x,-y,1/2+z'
sg_symmetry(8, 66) = '-x,y,1/2+z'
sg_symmetry(9, 66) = '1/2+x,1/2+y,z'
sg_symmetry(10, 66) = '1/2-x,1/2-y,z'
sg_symmetry(11, 66) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(12, 66) = '1/2+x,1/2-y,1/2-z'
sg_symmetry(13, 66) = '1/2-x,1/2-y,-z'
sg_symmetry(14, 66) = '1/2+x,1/2+y,-z'
sg_symmetry(15, 66) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(16, 66) = '1/2-x,1/2+y,1/2+z'
sg_name(67) = 'Cmma'
sg_patn(67) = 65
sg_symnum(67) = 16
sg_symmetry(1, 67) = 'x,y,z'
sg_symmetry(2, 67) = '-x,1/2-y,z'
sg_symmetry(3, 67) = '-x,1/2+y,-z'
sg_symmetry(4, 67) = 'x,-y,-z'
sg_symmetry(5, 67) = '-x,-y,-z'
sg_symmetry(6, 67) = 'x,1/2+y,-z'
sg_symmetry(7, 67) = 'x,1/2-y,z'
sg_symmetry(8, 67) = '-x,y,z'
sg_symmetry(9, 67) = '1/2+x,1/2+y,z'
sg_symmetry(10, 67) = '1/2-x,-y,z'
sg_symmetry(11, 67) = '1/2-x,y,-z'
sg_symmetry(12, 67) = '1/2+x,1/2-y,-z'
sg_symmetry(13, 67) = '1/2-x,1/2-y,-z'
sg_symmetry(14, 67) = '1/2+x,y,-z'
sg_symmetry(15, 67) = '1/2+x,-y,z'
sg_symmetry(16, 67) = '1/2-x,1/2+y,z'
sg_name(68) = 'Ccca'
sg_patn(68) = 65
sg_symnum(68) = 16
sg_symmetry(1, 68) = 'x,y,z'
sg_symmetry(2, 68) = '1/2-x,1/2-y,z'
sg_symmetry(3, 68) = '-x,y,-z'
sg_symmetry(4, 68) = '1/2+x,1/2-y,-z'
sg_symmetry(5, 68) = '-x,1/2-y,1/2-z'
sg_symmetry(6, 68) = '1/2+x,y,1/2-z'
sg_symmetry(7, 68) = 'x,1/2-y,1/2+z'
sg_symmetry(8, 68) = '1/2-x,y,1/2+z'
sg_symmetry(9, 68) = '1/2+x,1/2+y,z'
sg_symmetry(10, 68) = '-x,-y,z'
sg_symmetry(11, 68) = '1/2-x,1/2+y,-z'
sg_symmetry(12, 68) = 'x,-y,-z'
sg_symmetry(13, 68) = '1/2-x,-y,1/2-z'
sg_symmetry(14, 68) = 'x,1/2+y,1/2-z'
sg_symmetry(15, 68) = '1/2+x,-y,1/2+z'
sg_symmetry(16, 68) = '-x,1/2+y,1/2+z'
sg_name(69) = 'Fmmm'
sg_patn(69) = 69
sg_symnum(69) = 32
sg_symmetry(1, 69) = 'x,y,z'
sg_symmetry(2, 69) = '-x,-y,z'
sg_symmetry(3, 69) = '-x,y,-z'
sg_symmetry(4, 69) = 'x,-y,-z'
sg_symmetry(5, 69) = '-x,-y,-z'
sg_symmetry(6, 69) = 'x,y,-z'
sg_symmetry(7, 69) = 'x,-y,z'
sg_symmetry(8, 69) = '-x,y,z'
sg_symmetry(9, 69) = 'x,1/2+y,1/2+z'
sg_symmetry(10, 69) = '-x,1/2-y,1/2+z'
sg_symmetry(11, 69) = '-x,1/2+y,1/2-z'
sg_symmetry(12, 69) = 'x,1/2-y,1/2-z'
sg_symmetry(13, 69) = '-x,1/2-y,1/2-z'
sg_symmetry(14, 69) = 'x,1/2+y,1/2-z'
sg_symmetry(15, 69) = 'x,1/2-y,1/2+z'
sg_symmetry(16, 69) = '-x,1/2+y,1/2+z'
sg_symmetry(17, 69) = '1/2+x,y,1/2+z'
sg_symmetry(18, 69) = '1/2-x,-y,1/2+z'
sg_symmetry(19, 69) = '1/2-x,y,1/2-z'
sg_symmetry(20, 69) = '1/2+x,-y,1/2-z'
sg_symmetry(21, 69) = '1/2-x,-y,1/2-z'
sg_symmetry(22, 69) = '1/2+x,y,1/2-z'
sg_symmetry(23, 69) = '1/2+x,-y,1/2+z'
sg_symmetry(24, 69) = '1/2-x,y,1/2+z'
sg_symmetry(25, 69) = '1/2+x,1/2+y,z'
sg_symmetry(26, 69) = '1/2-x,1/2-y,z'
sg_symmetry(27, 69) = '1/2-x,1/2+y,-z'
sg_symmetry(28, 69) = '1/2+x,1/2-y,-z'
sg_symmetry(29, 69) = '1/2-x,1/2-y,-z'
sg_symmetry(30, 69) = '1/2+x,1/2+y,-z'
sg_symmetry(31, 69) = '1/2+x,1/2-y,z'
sg_symmetry(32, 69) = '1/2-x,1/2+y,z'
sg_name(70) = 'Fddd'
sg_patn(70) = 69
sg_symnum(70) = 32
sg_symmetry(1, 70) = 'x,y,z'
sg_symmetry(2, 70) = '-x,-y,z'
sg_symmetry(3, 70) = '-x,y,-z'
sg_symmetry(4, 70) = 'x,-y,-z'
sg_symmetry(5, 70) = '1/4-x,1/4-y,1/4-z'
sg_symmetry(6, 70) = '1/4+x,1/4+y,1/4-z'
sg_symmetry(7, 70) = '1/4+x,1/4-y,1/4+z'
sg_symmetry(8, 70) = '1/4-x,1/4+y,1/4+z'
sg_symmetry(9, 70) = 'x,1/2+y,1/2+z'
sg_symmetry(10, 70) = '-x,1/2-y,1/2+z'
sg_symmetry(11, 70) = '-x,1/2+y,1/2-z'
sg_symmetry(12, 70) = 'x,1/2-y,1/2-z'
sg_symmetry(13, 70) = '1/4-x,3/4-y,3/4-z'
sg_symmetry(14, 70) = '1/4+x,3/4+y,3/4-z'
sg_symmetry(15, 70) = '1/4+x,3/4-y,3/4+z'
sg_symmetry(16, 70) = '1/4-x,3/4+y,3/4+z'
sg_symmetry(17, 70) = '1/2+x,y,1/2+z'
sg_symmetry(18, 70) = '1/2-x,-y,1/2+z'
sg_symmetry(19, 70) = '1/2-x,y,1/2-z'
sg_symmetry(20, 70) = '1/2+x,-y,1/2-z'
sg_symmetry(21, 70) = '3/4-x,1/4-y,3/4-z'
sg_symmetry(22, 70) = '3/4+x,1/4+y,3/4-z'
sg_symmetry(23, 70) = '3/4+x,1/4-y,3/4+z'
sg_symmetry(24, 70) = '3/4-x,1/4+y,3/4+z'
sg_symmetry(25, 70) = '1/2+x,1/2+y,z'
sg_symmetry(26, 70) = '1/2-x,1/2-y,z'
sg_symmetry(27, 70) = '1/2-x,1/2+y,-z'
sg_symmetry(28, 70) = '1/2+x,1/2-y,-z'
sg_symmetry(29, 70) = '3/4-x,3/4-y,1/4-z'
sg_symmetry(30, 70) = '3/4+x,3/4+y,1/4-z'
sg_symmetry(31, 70) = '3/4+x,3/4-y,1/4+z'
sg_symmetry(32, 70) = '3/4-x,3/4+y,1/4+z'
sg_name(71) = 'Immm'
sg_patn(71) = 71
sg_symnum(71) = 16
sg_symmetry(1, 71) = 'x,y,z'
sg_symmetry(2, 71) = '-x,-y,z'
sg_symmetry(3, 71) = '-x,y,-z'
sg_symmetry(4, 71) = 'x,-y,-z'
sg_symmetry(5, 71) = '-x,-y,-z'
sg_symmetry(6, 71) = 'x,y,-z'
sg_symmetry(7, 71) = 'x,-y,z'
sg_symmetry(8, 71) = '-x,y,z'
sg_symmetry(9, 71) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(10, 71) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(11, 71) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(12, 71) = '1/2+x,1/2-y,1/2-z'
sg_symmetry(13, 71) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(14, 71) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(15, 71) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(16, 71) = '1/2-x,1/2+y,1/2+z'
sg_name(72) = 'Ibam'
sg_patn(72) = 71
sg_symnum(72) = 16
sg_symmetry(1, 72) = 'x,y,z'
sg_symmetry(2, 72) = '-x,-y,z'
sg_symmetry(3, 72) = '1/2-x,1/2+y,-z'
sg_symmetry(4, 72) = '1/2+x,1/2-y,-z'
sg_symmetry(5, 72) = '-x,-y,-z'
sg_symmetry(6, 72) = 'x,y,-z'
sg_symmetry(7, 72) = '1/2+x,1/2-y,z'
sg_symmetry(8, 72) = '1/2-x,1/2+y,z'
sg_symmetry(9, 72) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(10, 72) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(11, 72) = '-x,y,1/2-z'
sg_symmetry(12, 72) = 'x,-y,1/2-z'
sg_symmetry(13, 72) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(14, 72) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(15, 72) = 'x,-y,1/2+z'
sg_symmetry(16, 72) = '-x,y,1/2+z'
sg_name(73) = 'Ibca'
sg_patn(73) = 71
sg_symnum(73) = 16
sg_symmetry(1, 73) = 'x,y,z'
sg_symmetry(2, 73) = '1/2-x,-y,1/2+z'
sg_symmetry(3, 73) = '-x,1/2+y,1/2-z'
sg_symmetry(4, 73) = '1/2+x,1/2-y,-z'
sg_symmetry(5, 73) = '-x,-y,-z'
sg_symmetry(6, 73) = '1/2+x,y,1/2-z'
sg_symmetry(7, 73) = 'x,1/2-y,1/2+z'
sg_symmetry(8, 73) = '1/2-x,1/2+y,z'
sg_symmetry(9, 73) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(10, 73) = '-x,1/2-y,z'
sg_symmetry(11, 73) = '1/2-x,y,-z'
sg_symmetry(12, 73) = 'x,-y,1/2-z'
sg_symmetry(13, 73) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(14, 73) = 'x,1/2+y,-z'
sg_symmetry(15, 73) = '1/2+x,-y,z'
sg_symmetry(16, 73) = '-x,y,1/2+z'
sg_name(74) = 'Imma'
sg_patn(74) = 71
sg_symnum(74) = 16
sg_symmetry(1, 74) = 'x,y,z'
sg_symmetry(2, 74) = '-x,1/2-y,z'
sg_symmetry(3, 74) = '-x,1/2+y,-z'
sg_symmetry(4, 74) = 'x,-y,-z'
sg_symmetry(5, 74) = '-x,-y,-z'
sg_symmetry(6, 74) = 'x,1/2+y,-z'
sg_symmetry(7, 74) = 'x,1/2-y,z'
sg_symmetry(8, 74) = '-x,y,z'
sg_symmetry(9, 74) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(10, 74) = '1/2-x,-y,1/2+z'
sg_symmetry(11, 74) = '1/2-x,y,1/2-z'
sg_symmetry(12, 74) = '1/2+x,1/2-y,1/2-z'
sg_symmetry(13, 74) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(14, 74) = '1/2+x,y,1/2-z'
sg_symmetry(15, 74) = '1/2+x,-y,1/2+z'
sg_symmetry(16, 74) = '1/2-x,1/2+y,1/2+z'
sg_name(75) = 'P4'
sg_patn(75) = 83
sg_symnum(75) = 4
sg_symmetry(1, 75) = 'x,y,z'
sg_symmetry(2, 75) = '-x,-y,z'
sg_symmetry(3, 75) = '-y,x,z'
sg_symmetry(4, 75) = 'y,-x,z'
sg_name(76) = 'P4(1)'
sg_patn(76) = 83
sg_symnum(76) = 4
sg_symmetry(1, 76) = 'x,y,z'
sg_symmetry(2, 76) = '-x,-y,1/2+z'
sg_symmetry(3, 76) = '-y,x,1/4+z'
sg_symmetry(4, 76) = 'y,-x,3/4+z'
sg_name(77) = 'P4(2)'
sg_patn(77) = 83
sg_symnum(77) = 4
sg_symmetry(1, 77) = 'x,y,z'
sg_symmetry(2, 77) = '-x,-y,z'
sg_symmetry(3, 77) = '-y,x,1/2+z'
sg_symmetry(4, 77) = 'y,-x,1/2+z'
sg_name(78) = 'P4(3)'
sg_patn(78) = 83
sg_symnum(78) = 4
sg_symmetry(1, 78) = 'x,y,z'
sg_symmetry(2, 78) = '-x,-y,1/2+z'
sg_symmetry(3, 78) = '-y,x,3/4+z'
sg_symmetry(4, 78) = 'y,-x,1/4+z'
sg_name(79) = 'I4'
sg_patn(79) = 87
sg_symnum(79) = 8
sg_symmetry(1, 79) = 'x,y,z'
sg_symmetry(2, 79) = '-x,-y,z'
sg_symmetry(3, 79) = '-y,x,z'
sg_symmetry(4, 79) = 'y,-x,z'
sg_symmetry(5, 79) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(6, 79) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(7, 79) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(8, 79) = '1/2+y,1/2-x,1/2+z'
sg_name(80) = 'I4(1)'
sg_patn(80) = 87
sg_symnum(80) = 8
sg_symmetry(1, 80) = 'x,y,z'
sg_symmetry(2, 80) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(3, 80) = '-y,1/2+x,1/4+z'
sg_symmetry(4, 80) = '1/2+y,-x,3/4+z'
sg_symmetry(5, 80) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(6, 80) = '-x,-y,z'
sg_symmetry(7, 80) = '1/2-y,x,3/4+z'
sg_symmetry(8, 80) = 'y,1/2-x,1/4+z'
sg_name(81) = 'P-4'
sg_patn(81) = 83
sg_symnum(81) = 4
sg_symmetry(1, 81) = 'x,y,z'
sg_symmetry(2, 81) = '-x,-y,z'
sg_symmetry(3, 81) = 'y,-x,-z'
sg_symmetry(4, 81) = '-y,x,-z'
sg_name(82) = 'I-4'
sg_patn(82) = 87
sg_symnum(82) = 8
sg_symmetry(1, 82) = 'x,y,z'
sg_symmetry(2, 82) = '-x,-y,z'
sg_symmetry(3, 82) = 'y,-x,-z'
sg_symmetry(4, 82) = '-y,x,-z'
sg_symmetry(5, 82) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(6, 82) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(7, 82) = '1/2+y,1/2-x,1/2-z'
sg_symmetry(8, 82) = '1/2-y,1/2+x,1/2-z'
sg_name(83) = 'P4/m'
sg_patn(83) = 83
sg_symnum(83) = 8
sg_symmetry(1, 83) = 'x,y,z'
sg_symmetry(2, 83) = '-x,-y,z'
sg_symmetry(3, 83) = '-y,x,z'
sg_symmetry(4, 83) = 'y,-x,z'
sg_symmetry(5, 83) = '-x,-y,-z'
sg_symmetry(6, 83) = 'x,y,-z'
sg_symmetry(7, 83) = 'y,-x,-z'
sg_symmetry(8, 83) = '-y,x,-z'
sg_name(84) = 'P4(2)/m'
sg_patn(84) = 83
sg_symnum(84) = 8
sg_symmetry(1, 84) = 'x,y,z'
sg_symmetry(2, 84) = '-x,-y,z'
sg_symmetry(3, 84) = '-y,x,1/2+z'
sg_symmetry(4, 84) = 'y,-x,1/2+z'
sg_symmetry(5, 84) = '-x,-y,-z'
sg_symmetry(6, 84) = 'x,y,-z'
sg_symmetry(7, 84) = 'y,-x,1/2-z'
sg_symmetry(8, 84) = '-y,x,1/2-z'
sg_name(85) = 'P4/n'
sg_patn(85) = 83
sg_symnum(85) = 8
sg_symmetry(1, 85) = 'x,y,z'
sg_symmetry(2, 85) = '-x,-y,z'
sg_symmetry(3, 85) = '1/2-y,1/2+x,z'
sg_symmetry(4, 85) = '1/2+y,1/2-x,z'
sg_symmetry(5, 85) = '1/2-x,1/2-y,-z'
sg_symmetry(6, 85) = '1/2+x,1/2+y,-z'
sg_symmetry(7, 85) = 'y,-x,-z'
sg_symmetry(8, 85) = '-y,x,-z'
sg_name(86) = 'P4(2)/n'
sg_patn(86) = 83
sg_symnum(86) = 8
sg_symmetry(1, 86) = 'x,y,z'
sg_symmetry(2, 86) = '-x,-y,z'
sg_symmetry(3, 86) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(4, 86) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(5, 86) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(6, 86) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(7, 86) = 'y,-x,-z'
sg_symmetry(8, 86) = '-y,x,-z'
sg_name(87) = 'I4/m'
sg_patn(87) = 87
sg_symnum(87) = 16
sg_symmetry(1, 87) = 'x,y,z'
sg_symmetry(2, 87) = '-x,-y,z'
sg_symmetry(3, 87) = '-y,x,z'
sg_symmetry(4, 87) = 'y,-x,z'
sg_symmetry(5, 87) = '-x,-y,-z'
sg_symmetry(6, 87) = 'x,y,-z'
sg_symmetry(7, 87) = 'y,-x,-z'
sg_symmetry(8, 87) = '-y,x,-z'
sg_symmetry(9, 87) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(10, 87) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(11, 87) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(12, 87) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(13, 87) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(14, 87) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(15, 87) = '1/2+y,1/2-x,1/2-z'
sg_symmetry(16, 87) = '1/2-y,1/2+x,1/2-z'
sg_name(88) = 'I4(1)/a'
sg_patn(88) = 87
sg_symnum(88) = 16
sg_symmetry(1, 88) = 'x,y,z'
sg_symmetry(2, 88) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(3, 88) = '-y,1/2+x,1/4+z'
sg_symmetry(4, 88) = '1/2+y,-x,3/4+z'
sg_symmetry(5, 88) = '-x,1/2-y,1/4-z'
sg_symmetry(6, 88) = '1/2+x,y,3/4-z'
sg_symmetry(7, 88) = 'y,-x,-z'
sg_symmetry(8, 88) = '1/2-y,1/2+x,1/2-z'
sg_symmetry(9, 88) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(10, 88) = '-x,-y,z'
sg_symmetry(11, 88) = '1/2-y,x,3/4+z'
sg_symmetry(12, 88) = 'y,1/2-x,1/4+z'
sg_symmetry(13, 88) = '1/2-x,-y,3/4-z'
sg_symmetry(14, 88) = 'x,1/2+y,1/4-z'
sg_symmetry(15, 88) = '1/2+y,1/2-x,1/2-z'
sg_symmetry(16, 88) = '-y,x,-z'
sg_name(89) = 'P422'
sg_patn(89) = 123
sg_symnum(89) = 8
sg_symmetry(1, 89) = 'x,y,z'
sg_symmetry(2, 89) = '-x,-y,z'
sg_symmetry(3, 89) = '-y,x,z'
sg_symmetry(4, 89) = 'y,-x,z'
sg_symmetry(5, 89) = '-x,y,-z'
sg_symmetry(6, 89) = 'x,-y,-z'
sg_symmetry(7, 89) = 'y,x,-z'
sg_symmetry(8, 89) = '-y,-x,-z'
sg_name(90) = 'P42(1)2'
sg_patn(90) = 123
sg_symnum(90) = 8
sg_symmetry(1, 90) = 'x,y,z'
sg_symmetry(2, 90) = '-x,-y,z'
sg_symmetry(3, 90) = '1/2-y,1/2+x,z'
sg_symmetry(4, 90) = '1/2+y,1/2-x,z'
sg_symmetry(5, 90) = '1/2-x,1/2+y,-z'
sg_symmetry(6, 90) = '1/2+x,1/2-y,-z'
sg_symmetry(7, 90) = 'y,x,-z'
sg_symmetry(8, 90) = '-y,-x,-z'
sg_name(91) = 'P4(1)22'
sg_patn(91) = 123
sg_symnum(91) = 8
sg_symmetry(1, 91) = 'x,y,z'
sg_symmetry(2, 91) = '-x,-y,1/2+z'
sg_symmetry(3, 91) = '-y,x,1/4+z'
sg_symmetry(4, 91) = 'y,-x,3/4+z'
sg_symmetry(5, 91) = '-x,y,-z'
sg_symmetry(6, 91) = 'x,-y,1/2-z'
sg_symmetry(7, 91) = 'y,x,3/4-z'
sg_symmetry(8, 91) = '-y,-x,1/4-z'
sg_name(92) = 'P4(1)2(1)2'
sg_patn(92) = 123
sg_symnum(92) = 8
sg_symmetry(1, 92) = 'x,y,z'
sg_symmetry(2, 92) = '-x,-y,1/2+z'
sg_symmetry(3, 92) = '1/2-y,1/2+x,1/4+z'
sg_symmetry(4, 92) = '1/2+y,1/2-x,3/4+z'
sg_symmetry(5, 92) = '1/2-x,1/2+y,1/4-z'
sg_symmetry(6, 92) = '1/2+x,1/2-y,3/4-z'
sg_symmetry(7, 92) = 'y,x,-z'
sg_symmetry(8, 92) = '-y,-x,1/2-z'
sg_name(93) = 'P4(2)22'
sg_patn(93) = 123
sg_symnum(93) = 8
sg_symmetry(1, 93) = 'x,y,z'
sg_symmetry(2, 93) = '-x,-y,z'
sg_symmetry(3, 93) = '-y,x,1/2+z'
sg_symmetry(4, 93) = 'y,-x,1/2+z'
sg_symmetry(5, 93) = '-x,y,-z'
sg_symmetry(6, 93) = 'x,-y,-z'
sg_symmetry(7, 93) = 'y,x,1/2-z'
sg_symmetry(8, 93) = '-y,-x,1/2-z'
sg_name(94) = 'P4(2)2(1)2'
sg_patn(94) = 123
sg_symnum(94) = 8
sg_symmetry(1, 94) = 'x,y,z'
sg_symmetry(2, 94) = '-x,-y,z'
sg_symmetry(3, 94) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(4, 94) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(5, 94) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(6, 94) = '1/2+x,1/2-y,1/2-z'
sg_symmetry(7, 94) = 'y,x,-z'
sg_symmetry(8, 94) = '-y,-x,-z'
sg_name(95) = 'P4(3)22'
sg_patn(95) = 123
sg_symnum(95) = 8
sg_symmetry(1, 95) = 'x,y,z'
sg_symmetry(2, 95) = '-x,-y,1/2+z'
sg_symmetry(3, 95) = '-y,x,3/4+z'
sg_symmetry(4, 95) = 'y,-x,1/4+z'
sg_symmetry(5, 95) = '-x,y,-z'
sg_symmetry(6, 95) = 'x,-y,1/2-z'
sg_symmetry(7, 95) = 'y,x,1/4-z'
sg_symmetry(8, 95) = '-y,-x,3/4-z'
sg_name(96) = 'P4(3)2(1)2'
sg_patn(96) = 123
sg_symnum(96) = 8
sg_symmetry(1, 96) = 'x,y,z'
sg_symmetry(2, 96) = '-x,-y,1/2+z'
sg_symmetry(3, 96) = '1/2-y,1/2+x,3/4+z'
sg_symmetry(4, 96) = '1/2+y,1/2-x,1/4+z'
sg_symmetry(5, 96) = '1/2-x,1/2+y,3/4-z'
sg_symmetry(6, 96) = '1/2+x,1/2-y,1/4-z'
sg_symmetry(7, 96) = 'y,x,-z'
sg_symmetry(8, 96) = '-y,-x,1/2-z'
sg_name(97) = 'I422'
sg_patn(97) = 139
sg_symnum(97) = 16
sg_symmetry(1, 97) = 'x,y,z'
sg_symmetry(2, 97) = '-x,-y,z'
sg_symmetry(3, 97) = '-y,x,z'
sg_symmetry(4, 97) = 'y,-x,z'
sg_symmetry(5, 97) = '-x,y,-z'
sg_symmetry(6, 97) = 'x,-y,-z'
sg_symmetry(7, 97) = 'y,x,-z'
sg_symmetry(8, 97) = '-y,-x,-z'
sg_symmetry(9, 97) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(10, 97) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(11, 97) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(12, 97) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(13, 97) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(14, 97) = '1/2+x,1/2-y,1/2-z'
sg_symmetry(15, 97) = '1/2+y,1/2+x,1/2-z'
sg_symmetry(16, 97) = '1/2-y,1/2-x,1/2-z'
sg_name(98) = 'I4(1)22'
sg_patn(98) = 139
sg_symnum(98) = 16
sg_symmetry(1, 98) = 'x,y,z'
sg_symmetry(2, 98) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(3, 98) = '-y,1/2+x,1/4+z'
sg_symmetry(4, 98) = '1/2+y,-x,3/4+z'
sg_symmetry(5, 98) = '1/2-x,y,3/4-z'
sg_symmetry(6, 98) = 'x,1/2-y,1/4-z'
sg_symmetry(7, 98) = '1/2+y,1/2+x,1/2-z'
sg_symmetry(8, 98) = '-y,-x,-z'
sg_symmetry(9, 98) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(10, 98) = '-x,-y,z'
sg_symmetry(11, 98) = '1/2-y,x,3/4+z'
sg_symmetry(12, 98) = 'y,1/2-x,1/4+z'
sg_symmetry(13, 98) = '-x,1/2+y,1/4-z'
sg_symmetry(14, 98) = '1/2+x,-y,3/4-z'
sg_symmetry(15, 98) = 'y,x,-z'
sg_symmetry(16, 98) = '1/2-y,1/2-x,1/2-z'
sg_name(99) = 'P4mm'
sg_patn(99) = 123
sg_symnum(99) = 8
sg_symmetry(1, 99) = 'x,y,z'
sg_symmetry(2, 99) = '-x,-y,z'
sg_symmetry(3, 99) = '-y,x,z'
sg_symmetry(4, 99) = 'y,-x,z'
sg_symmetry(5, 99) = 'x,-y,z'
sg_symmetry(6, 99) = '-x,y,z'
sg_symmetry(7, 99) = '-y,-x,z'
sg_symmetry(8, 99) = 'y,x,z'
sg_name(100) = 'P4bm'
sg_patn(100) = 123
sg_symnum(100) = 8
sg_symmetry(1, 100) = 'x,y,z'
sg_symmetry(2, 100) = '-x,-y,z'
sg_symmetry(3, 100) = '-y,x,z'
sg_symmetry(4, 100) = 'y,-x,z'
sg_symmetry(5, 100) = '1/2+x,1/2-y,z'
sg_symmetry(6, 100) = '1/2-x,1/2+y,z'
sg_symmetry(7, 100) = '1/2-y,1/2-x,z'
sg_symmetry(8, 100) = '1/2+y,1/2+x,z'
sg_name(101) = 'P4(2)cm'
sg_patn(101) = 123
sg_symnum(101) = 8
sg_symmetry(1, 101) = 'x,y,z'
sg_symmetry(2, 101) = '-x,-y,z'
sg_symmetry(3, 101) = '-y,x,1/2+z'
sg_symmetry(4, 101) = 'y,-x,1/2+z'
sg_symmetry(5, 101) = 'x,-y,1/2+z'
sg_symmetry(6, 101) = '-x,y,1/2+z'
sg_symmetry(7, 101) = '-y,-x,z'
sg_symmetry(8, 101) = 'y,x,z'
sg_name(102) = 'P4(2)nm'
sg_patn(102) = 123
sg_symnum(102) = 8
sg_symmetry(1, 102) = 'x,y,z'
sg_symmetry(2, 102) = '-x,-y,z'
sg_symmetry(3, 102) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(4, 102) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(5, 102) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(6, 102) = '1/2-x,1/2+y,1/2+z'
sg_symmetry(7, 102) = '-y,-x,z'
sg_symmetry(8, 102) = 'y,x,z'
sg_name(103) = 'P4cc'
sg_patn(103) = 123
sg_symnum(103) = 8
sg_symmetry(1, 103) = 'x,y,z'
sg_symmetry(2, 103) = '-x,-y,z'
sg_symmetry(3, 103) = '-y,x,z'
sg_symmetry(4, 103) = 'y,-x,z'
sg_symmetry(5, 103) = 'x,-y,1/2+z'
sg_symmetry(6, 103) = '-x,y,1/2+z'
sg_symmetry(7, 103) = '-y,-x,1/2+z'
sg_symmetry(8, 103) = 'y,x,1/2+z'
sg_name(104) = 'P4nc'
sg_patn(104) = 123
sg_symnum(104) = 8
sg_symmetry(1, 104) = 'x,y,z'
sg_symmetry(2, 104) = '-x,-y,z'
sg_symmetry(3, 104) = '-y,x,z'
sg_symmetry(4, 104) = 'y,-x,z'
sg_symmetry(5, 104) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(6, 104) = '1/2-x,1/2+y,1/2+z'
sg_symmetry(7, 104) = '1/2-y,1/2-x,1/2+z'
sg_symmetry(8, 104) = '1/2+y,1/2+x,1/2+z'
sg_name(105) = 'P4(2)mc'
sg_patn(105) = 123
sg_symnum(105) = 8
sg_symmetry(1, 105) = 'x,y,z'
sg_symmetry(2, 105) = '-x,-y,z'
sg_symmetry(3, 105) = '-y,x,1/2+z'
sg_symmetry(4, 105) = 'y,-x,1/2+z'
sg_symmetry(5, 105) = 'x,-y,z'
sg_symmetry(6, 105) = '-x,y,z'
sg_symmetry(7, 105) = '-y,-x,1/2+z'
sg_symmetry(8, 105) = 'y,x,1/2+z'
sg_name(106) = 'P4(2)bc'
sg_patn(106) = 123
sg_symnum(106) = 8
sg_symmetry(1, 106) = 'x,y,z'
sg_symmetry(2, 106) = '-x,-y,z'
sg_symmetry(3, 106) = '-y,x,1/2+z'
sg_symmetry(4, 106) = 'y,-x,1/2+z'
sg_symmetry(5, 106) = '1/2+x,1/2-y,z'
sg_symmetry(6, 106) = '1/2-x,1/2+y,z'
sg_symmetry(7, 106) = '1/2-y,1/2-x,1/2+z'
sg_symmetry(8, 106) = '1/2+y,1/2+x,1/2+z'
sg_name(107) = 'I4mm'
sg_patn(107) = 139
sg_symnum(107) = 16
sg_symmetry(1, 107) = 'x,y,z'
sg_symmetry(2, 107) = '-x,-y,z'
sg_symmetry(3, 107) = '-y,x,z'
sg_symmetry(4, 107) = 'y,-x,z'
sg_symmetry(5, 107) = 'x,-y,z'
sg_symmetry(6, 107) = '-x,y,z'
sg_symmetry(7, 107) = '-y,-x,z'
sg_symmetry(8, 107) = 'y,x,z'
sg_symmetry(9, 107) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(10, 107) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(11, 107) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(12, 107) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(13, 107) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(14, 107) = '1/2-x,1/2+y,1/2+z'
sg_symmetry(15, 107) = '1/2-y,1/2-x,1/2+z'
sg_symmetry(16, 107) = '1/2+y,1/2+x,1/2+z'
sg_name(108) = 'I4cm'
sg_patn(108) = 139
sg_symnum(108) = 16
sg_symmetry(1, 108) = 'x,y,z'
sg_symmetry(2, 108) = '-x,-y,z'
sg_symmetry(3, 108) = '-y,x,z'
sg_symmetry(4, 108) = 'y,-x,z'
sg_symmetry(5, 108) = 'x,-y,1/2+z'
sg_symmetry(6, 108) = '-x,y,1/2+z'
sg_symmetry(7, 108) = '-y,-x,1/2+z'
sg_symmetry(8, 108) = 'y,x,1/2+z'
sg_symmetry(9, 108) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(10, 108) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(11, 108) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(12, 108) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(13, 108) = '1/2+x,1/2-y,z'
sg_symmetry(14, 108) = '1/2-x,1/2+y,z'
sg_symmetry(15, 108) = '1/2-y,1/2-x,z'
sg_symmetry(16, 108) = '1/2+y,1/2+x,z'
sg_name(109) = 'I4(1)md'
sg_patn(109) = 139
sg_symnum(109) = 16
sg_symmetry(1, 109) = 'x,y,z'
sg_symmetry(2, 109) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(3, 109) = '-y,1/2+x,1/4+z'
sg_symmetry(4, 109) = '1/2+y,-x,3/4+z'
sg_symmetry(5, 109) = 'x,-y,z'
sg_symmetry(6, 109) = '1/2-x,1/2+y,1/2+z'
sg_symmetry(7, 109) = '-y,1/2-x,1/4+z'
sg_symmetry(8, 109) = '1/2+y,x,3/4+z'
sg_symmetry(9, 109) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(10, 109) = '-x,-y,z'
sg_symmetry(11, 109) = '1/2-y,x,3/4+z'
sg_symmetry(12, 109) = 'y,1/2-x,1/4+z'
sg_symmetry(13, 109) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(14, 109) = '-x,y,z'
sg_symmetry(15, 109) = '1/2-y,-x,3/4+z'
sg_symmetry(16, 109) = 'y,1/2+x,1/4+z'
sg_name(110) = 'I4(1)cd'
sg_patn(110) = 139
sg_symnum(110) = 16
sg_symmetry(1, 110) = 'x,y,z'
sg_symmetry(2, 110) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(3, 110) = '-y,1/2+x,1/4+z'
sg_symmetry(4, 110) = '1/2+y,-x,3/4+z'
sg_symmetry(5, 110) = 'x,-y,1/2+z'
sg_symmetry(6, 110) = '1/2-x,1/2+y,z'
sg_symmetry(7, 110) = '-y,1/2-x,3/4+z'
sg_symmetry(8, 110) = '1/2+y,x,1/4+z'
sg_symmetry(9, 110) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(10, 110) = '-x,-y,z'
sg_symmetry(11, 110) = '1/2-y,x,3/4+z'
sg_symmetry(12, 110) = 'y,1/2-x,1/4+z'
sg_symmetry(13, 110) = '1/2+x,1/2-y,z'
sg_symmetry(14, 110) = '-x,y,1/2+z'
sg_symmetry(15, 110) = '1/2-y,-x,1/4+z'
sg_symmetry(16, 110) = 'y,1/2+x,3/4+z'
sg_name(111) = 'P-42m'
sg_patn(111) = 123
sg_symnum(111) = 8
sg_symmetry(1, 111) = 'x,y,z'
sg_symmetry(2, 111) = '-x,-y,z'
sg_symmetry(3, 111) = '-y,x,-z'
sg_symmetry(4, 111) = 'y,-x,-z'
sg_symmetry(5, 111) = '-x,y,-z'
sg_symmetry(6, 111) = 'x,-y,-z'
sg_symmetry(7, 111) = '-y,-x,z'
sg_symmetry(8, 111) = 'y,x,z'
sg_name(112) = 'P-42c'
sg_patn(112) = 123
sg_symnum(112) = 8
sg_symmetry(1, 112) = 'x,y,z'
sg_symmetry(2, 112) = '-x,-y,z'
sg_symmetry(3, 112) = '-y,x,-z'
sg_symmetry(4, 112) = 'y,-x,-z'
sg_symmetry(5, 112) = '-x,y,1/2-z'
sg_symmetry(6, 112) = 'x,-y,1/2-z'
sg_symmetry(7, 112) = '-y,-x,1/2+z'
sg_symmetry(8, 112) = 'y,x,1/2+z'
sg_name(113) = 'P-42(1)m'
sg_patn(113) = 123
sg_symnum(113) = 8
sg_symmetry(1, 113) = 'x,y,z'
sg_symmetry(2, 113) = '-x,-y,z'
sg_symmetry(3, 113) = '-y,x,-z'
sg_symmetry(4, 113) = 'y,-x,-z'
sg_symmetry(5, 113) = '1/2-x,1/2+y,-z'
sg_symmetry(6, 113) = '1/2+x,1/2-y,-z'
sg_symmetry(7, 113) = '1/2-y,1/2-x,z'
sg_symmetry(8, 113) = '1/2+y,1/2+x,z'
sg_name(114) = 'P-42(1)c'
sg_patn(114) = 123
sg_symnum(114) = 8
sg_symmetry(1, 114) = 'x,y,z'
sg_symmetry(2, 114) = '-x,-y,z'
sg_symmetry(3, 114) = '-y,x,-z'
sg_symmetry(4, 114) = 'y,-x,-z'
sg_symmetry(5, 114) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(6, 114) = '1/2+x,1/2-y,1/2-z'
sg_symmetry(7, 114) = '1/2-y,1/2-x,1/2+z'
sg_symmetry(8, 114) = '1/2+y,1/2+x,1/2+z'
sg_name(115) = 'P-4m2'
sg_patn(115) = 123
sg_symnum(115) = 8
sg_symmetry(1, 115) = 'x,y,z'
sg_symmetry(2, 115) = '-x,-y,z'
sg_symmetry(3, 115) = 'y,-x,-z'
sg_symmetry(4, 115) = '-y,x,-z'
sg_symmetry(5, 115) = 'x,-y,z'
sg_symmetry(6, 115) = '-x,y,z'
sg_symmetry(7, 115) = 'y,x,-z'
sg_symmetry(8, 115) = '-y,-x,-z'
sg_name(116) = 'P-4c2'
sg_patn(116) = 123
sg_symnum(116) = 8
sg_symmetry(1, 116) = 'x,y,z'
sg_symmetry(2, 116) = '-x,-y,z'
sg_symmetry(3, 116) = '-y,x,-z'
sg_symmetry(4, 116) = 'y,-x,-z'
sg_symmetry(5, 116) = 'x,-y,1/2+z'
sg_symmetry(6, 116) = '-x,y,1/2+z'
sg_symmetry(7, 116) = 'y,x,1/2-z'
sg_symmetry(8, 116) = '-y,-x,1/2-z'
sg_name(117) = 'P-4b2'
sg_patn(117) = 123
sg_symnum(117) = 8
sg_symmetry(1, 117) = 'x,y,z'
sg_symmetry(2, 117) = '-x,-y,z'
sg_symmetry(3, 117) = '-y,x,-z'
sg_symmetry(4, 117) = 'y,-x,-z'
sg_symmetry(5, 117) = '1/2+x,1/2-y,z'
sg_symmetry(6, 117) = '1/2-x,1/2+y,z'
sg_symmetry(7, 117) = '1/2+y,1/2+x,-z'
sg_symmetry(8, 117) = '1/2-y,1/2-x,-z'
sg_name(118) = 'P-4n2'
sg_patn(118) = 123
sg_symnum(118) = 8
sg_symmetry(1, 118) = 'x,y,z'
sg_symmetry(2, 118) = '-x,-y,z'
sg_symmetry(3, 118) = '-y,x,-z'
sg_symmetry(4, 118) = 'y,-x,-z'
sg_symmetry(5, 118) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(6, 118) = '1/2-x,1/2+y,1/2+z'
sg_symmetry(7, 118) = '1/2+y,1/2+x,1/2-z'
sg_symmetry(8, 118) = '1/2-y,1/2-x,1/2-z'
sg_name(119) = 'I-4m2'
sg_patn(119) = 139
sg_symnum(119) = 16
sg_symmetry(1, 119) = 'x,y,z'
sg_symmetry(2, 119) = '-x,-y,z'
sg_symmetry(3, 119) = '-y,x,-z'
sg_symmetry(4, 119) = 'y,-x,-z'
sg_symmetry(5, 119) = 'x,-y,z'
sg_symmetry(6, 119) = '-x,y,z'
sg_symmetry(7, 119) = 'y,x,-z'
sg_symmetry(8, 119) = '-y,-x,-z'
sg_symmetry(9, 119) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(10, 119) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(11, 119) = '1/2-y,1/2+x,1/2-z'
sg_symmetry(12, 119) = '1/2+y,1/2-x,1/2-z'
sg_symmetry(13, 119) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(14, 119) = '1/2-x,1/2+y,1/2+z'
sg_symmetry(15, 119) = '1/2+y,1/2+x,1/2-z'
sg_symmetry(16, 119) = '1/2-y,1/2-x,1/2-z'
sg_name(120) = 'I-4c2'
sg_patn(120) = 139
sg_symnum(120) = 16
sg_symmetry(1, 120) = 'x,y,z'
sg_symmetry(2, 120) = '-x,-y,z'
sg_symmetry(3, 120) = 'y,-x,-z'
sg_symmetry(4, 120) = '-y,x,-z'
sg_symmetry(5, 120) = 'x,-y,1/2+z'
sg_symmetry(6, 120) = '-x,y,1/2+z'
sg_symmetry(7, 120) = 'y,x,1/2-z'
sg_symmetry(8, 120) = '-y,-x,1/2-z'
sg_symmetry(9, 120) = 'x+1/2,y+1/2,z+1/2'
sg_symmetry(10, 120) = '-x+1/2,-y+1/2,z+1/2'
sg_symmetry(11, 120) = 'y+1/2,-x+1/2,-z+1/2'
sg_symmetry(12, 120) = '-y+1/2,x+1/2,-z+1/2'
sg_symmetry(13, 120) = 'x+1/2,-y+1/2,z'
sg_symmetry(14, 120) = '-x+1/2,y+1/2,z'
sg_symmetry(15, 120) = 'y+1/2,x+1/2,-z'
sg_symmetry(16, 120) = '-y+1/2,-x+1/2,-z'
sg_name(121) = 'I-42m'
sg_patn(121) = 139
sg_symnum(121) = 16
sg_symmetry(1, 121) = 'x,y,z'
sg_symmetry(2, 121) = '-x,-y,z'
sg_symmetry(3, 121) = '-y,x,-z'
sg_symmetry(4, 121) = 'y,-x,-z'
sg_symmetry(5, 121) = '-x,y,-z'
sg_symmetry(6, 121) = 'x,-y,-z'
sg_symmetry(7, 121) = '-y,-x,z'
sg_symmetry(8, 121) = 'y,x,z'
sg_symmetry(9, 121) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(10, 121) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(11, 121) = '1/2-y,1/2+x,1/2-z'
sg_symmetry(12, 121) = '1/2+y,1/2-x,1/2-z'
sg_symmetry(13, 121) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(14, 121) = '1/2+x,1/2-y,1/2-z'
sg_symmetry(15, 121) = '1/2-y,1/2-x,1/2+z'
sg_symmetry(16, 121) = '1/2+y,1/2+x,1/2+z'
sg_name(122) = 'I-42d'
sg_patn(122) = 139
sg_symnum(122) = 16
sg_symmetry(1, 122) = 'x,y,z'
sg_symmetry(2, 122) = '-x,-y,z'
sg_symmetry(3, 122) = '-y,x,-z'
sg_symmetry(4, 122) = 'y,-x,-z'
sg_symmetry(5, 122) = '1/2-x,y,3/4-z'
sg_symmetry(6, 122) = '1/2+x,-y,3/4-z'
sg_symmetry(7, 122) = '1/2-y,-x,3/4+z'
sg_symmetry(8, 122) = '1/2+y,x,3/4+z'
sg_symmetry(9, 122) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(10, 122) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(11, 122) = '1/2-y,1/2+x,1/2-z'
sg_symmetry(12, 122) = '1/2+y,1/2-x,1/2-z'
sg_symmetry(13, 122) = '-x,1/2+y,1/4-z'
sg_symmetry(14, 122) = 'x,1/2-y,1/4-z'
sg_symmetry(15, 122) = '-y,1/2-x,1/4+z'
sg_symmetry(16, 122) = 'y,1/2+x,1/4+z'
sg_name(123) = 'P4/mmm'
sg_patn(123) = 123
sg_symnum(123) = 16
sg_symmetry(1, 123) = 'x,y,z'
sg_symmetry(2, 123) = '-x,-y,z'
sg_symmetry(3, 123) = '-y,x,z'
sg_symmetry(4, 123) = 'y,-x,z'
sg_symmetry(5, 123) = '-x,y,-z'
sg_symmetry(6, 123) = 'x,-y,-z'
sg_symmetry(7, 123) = 'y,x,-z'
sg_symmetry(8, 123) = '-y,-x,-z'
sg_symmetry(9, 123) = '-x,-y,-z'
sg_symmetry(10, 123) = 'x,y,-z'
sg_symmetry(11, 123) = 'y,-x,-z'
sg_symmetry(12, 123) = '-y,x,-z'
sg_symmetry(13, 123) = 'x,-y,z'
sg_symmetry(14, 123) = '-x,y,z'
sg_symmetry(15, 123) = '-y,-x,z'
sg_symmetry(16, 123) = 'y,x,z'
sg_name(124) = 'P4/mcc'
sg_patn(124) = 123
sg_symnum(124) = 16
sg_symmetry(1, 124) = 'x,y,z'
sg_symmetry(2, 124) = '-x,-y,z'
sg_symmetry(3, 124) = '-y,x,z'
sg_symmetry(4, 124) = 'y,-x,z'
sg_symmetry(5, 124) = '-x,y,1/2-z'
sg_symmetry(6, 124) = 'x,-y,1/2-z'
sg_symmetry(7, 124) = 'y,x,1/2-z'
sg_symmetry(8, 124) = '-y,-x,1/2-z'
sg_symmetry(9, 124) = '-x,-y,-z'
sg_symmetry(10, 124) = 'x,y,-z'
sg_symmetry(11, 124) = 'y,-x,-z'
sg_symmetry(12, 124) = '-y,x,-z'
sg_symmetry(13, 124) = 'x,-y,1/2+z'
sg_symmetry(14, 124) = '-x,y,1/2+z'
sg_symmetry(15, 124) = '-y,-x,1/2+z'
sg_symmetry(16, 124) = 'y,x,1/2+z'
sg_name(125) = 'P4/nbm'
sg_patn(125) = 123
sg_symnum(125) = 16
sg_symmetry(1, 125) = 'x,y,z'
sg_symmetry(2, 125) = '-x,-y,z'
sg_symmetry(3, 125) = '-y,x,z'
sg_symmetry(4, 125) = 'y,-x,z'
sg_symmetry(5, 125) = '-x,y,-z'
sg_symmetry(6, 125) = 'x,-y,-z'
sg_symmetry(7, 125) = 'y,x,-z'
sg_symmetry(8, 125) = '-y,-x,-z'
sg_symmetry(9, 125) = '1/2-x,1/2-y,-z'
sg_symmetry(10, 125) = '1/2+x,1/2+y,-z'
sg_symmetry(11, 125) = '1/2+y,1/2-x,-z'
sg_symmetry(12, 125) = '1/2-y,1/2+x,-z'
sg_symmetry(13, 125) = '1/2+x,1/2-y,z'
sg_symmetry(14, 125) = '1/2-x,1/2+y,z'
sg_symmetry(15, 125) = '1/2-y,1/2-x,z'
sg_symmetry(16, 125) = '1/2+y,1/2+x,z'
sg_name(126) = 'P4/nnc'
sg_patn(126) = 123
sg_symnum(126) = 16
sg_symmetry(1, 126) = 'x,y,z'
sg_symmetry(2, 126) = '-x,-y,z'
sg_symmetry(3, 126) = '-y,x,z'
sg_symmetry(4, 126) = 'y,-x,z'
sg_symmetry(5, 126) = '-x,y,-z'
sg_symmetry(6, 126) = 'x,-y,-z'
sg_symmetry(7, 126) = 'y,x,-z'
sg_symmetry(8, 126) = '-y,-x,-z'
sg_symmetry(9, 126) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(10, 126) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(11, 126) = '1/2+y,1/2-x,1/2-z'
sg_symmetry(12, 126) = '1/2-y,1/2+x,1/2-z'
sg_symmetry(13, 126) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(14, 126) = '1/2-x,1/2+y,1/2+z'
sg_symmetry(15, 126) = '1/2-y,1/2-x,1/2+z'
sg_symmetry(16, 126) = '1/2+y,1/2+x,1/2+z'
sg_name(127) = 'P4/mbm'
sg_patn(127) = 123
sg_symnum(127) = 16
sg_symmetry(1, 127) = 'x,y,z'
sg_symmetry(2, 127) = '-x,-y,z'
sg_symmetry(3, 127) = '-y,x,z'
sg_symmetry(4, 127) = 'y,-x,z'
sg_symmetry(5, 127) = '1/2-x,1/2+y,-z'
sg_symmetry(6, 127) = '1/2+x,1/2-y,-z'
sg_symmetry(7, 127) = '1/2+y,1/2+x,-z'
sg_symmetry(8, 127) = '1/2-y,1/2-x,-z'
sg_symmetry(9, 127) = '-x,-y,-z'
sg_symmetry(10, 127) = 'x,y,-z'
sg_symmetry(11, 127) = 'y,-x,-z'
sg_symmetry(12, 127) = '-y,x,-z'
sg_symmetry(13, 127) = '1/2+x,1/2-y,z'
sg_symmetry(14, 127) = '1/2-x,1/2+y,z'
sg_symmetry(15, 127) = '1/2-y,1/2-x,z'
sg_symmetry(16, 127) = '1/2+y,1/2+x,z'
sg_name(128) = 'P4/mnc'
sg_patn(128) = 123
sg_symnum(128) = 16
sg_symmetry(1, 128) = 'x,y,z'
sg_symmetry(2, 128) = '-x,-y,z'
sg_symmetry(3, 128) = '-y,x,z'
sg_symmetry(4, 128) = 'y,-x,z'
sg_symmetry(5, 128) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(6, 128) = '1/2+x,1/2-y,1/2-z'
sg_symmetry(7, 128) = '1/2+y,1/2+x,1/2-z'
sg_symmetry(8, 128) = '1/2-y,1/2-x,1/2-z'
sg_symmetry(9, 128) = '-x,-y,-z'
sg_symmetry(10, 128) = 'x,y,-z'
sg_symmetry(11, 128) = 'y,-x,-z'
sg_symmetry(12, 128) = '-y,x,-z'
sg_symmetry(13, 128) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(14, 128) = '1/2-x,1/2+y,1/2+z'
sg_symmetry(15, 128) = '1/2-y,1/2-x,1/2+z'
sg_symmetry(16, 128) = '1/2+y,1/2+x,1/2+z'
sg_name(129) = 'P4/nmm'
sg_patn(129) = 123
sg_symnum(129) = 16
sg_symmetry(1, 129) = 'x,y,z'
sg_symmetry(2, 129) = '-x,-y,z'
sg_symmetry(3, 129) = '1/2-y,1/2+x,z'
sg_symmetry(4, 129) = '1/2+y,1/2-x,z'
sg_symmetry(5, 129) = '1/2-x,1/2+y,-z'
sg_symmetry(6, 129) = '1/2+x,1/2-y,-z'
sg_symmetry(7, 129) = 'y,x,-z'
sg_symmetry(8, 129) = '-y,-x,-z'
sg_symmetry(9, 129) = '1/2-x,1/2-y,-z'
sg_symmetry(10, 129) = '1/2+x,1/2+y,-z'
sg_symmetry(11, 129) = 'y,-x,-z'
sg_symmetry(12, 129) = '-y,x,-z'
sg_symmetry(13, 129) = 'x,-y,z'
sg_symmetry(14, 129) = '-x,y,z'
sg_symmetry(15, 129) = '1/2-y,1/2-x,z'
sg_symmetry(16, 129) = '1/2+y,1/2+x,z'
sg_name(130) = 'P4/ncc'
sg_patn(130) = 123
sg_symnum(130) = 16
sg_symmetry(1, 130) = 'x,y,z'
sg_symmetry(2, 130) = '-x,-y,z'
sg_symmetry(3, 130) = '1/2-y,1/2+x,z'
sg_symmetry(4, 130) = '1/2+y,1/2-x,z'
sg_symmetry(5, 130) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(6, 130) = '1/2+x,1/2-y,1/2-z'
sg_symmetry(7, 130) = 'y,x,1/2-z'
sg_symmetry(8, 130) = '-y,-x,1/2-z'
sg_symmetry(9, 130) = '1/2-x,1/2-y,-z'
sg_symmetry(10, 130) = '1/2+x,1/2+y,-z'
sg_symmetry(11, 130) = 'y,-x,-z'
sg_symmetry(12, 130) = '-y,x,-z'
sg_symmetry(13, 130) = 'x,-y,1/2+z'
sg_symmetry(14, 130) = '-x,y,1/2+z'
sg_symmetry(15, 130) = '1/2-y,1/2-x,1/2+z'
sg_symmetry(16, 130) = '1/2+y,1/2+x,1/2+z'
sg_name(131) = 'P4(2)/mmc'
sg_patn(131) = 123
sg_symnum(131) = 16
sg_symmetry(1, 131) = 'x,y,z'
sg_symmetry(2, 131) = '-x,-y,z'
sg_symmetry(3, 131) = '-y,x,1/2+z'
sg_symmetry(4, 131) = 'y,-x,1/2+z'
sg_symmetry(5, 131) = '-x,y,-z'
sg_symmetry(6, 131) = 'x,-y,-z'
sg_symmetry(7, 131) = 'y,x,1/2-z'
sg_symmetry(8, 131) = '-y,-x,1/2-z'
sg_symmetry(9, 131) = '-x,-y,-z'
sg_symmetry(10, 131) = 'x,y,-z'
sg_symmetry(11, 131) = 'y,-x,1/2-z'
sg_symmetry(12, 131) = '-y,x,1/2-z'
sg_symmetry(13, 131) = 'x,-y,z'
sg_symmetry(14, 131) = '-x,y,z'
sg_symmetry(15, 131) = '-y,-x,1/2+z'
sg_symmetry(16, 131) = 'y,x,1/2+z'
sg_name(132) = 'P4(2)/mcm'
sg_patn(132) = 123
sg_symnum(132) = 16
sg_symmetry(1, 132) = 'x,y,z'
sg_symmetry(2, 132) = '-x,-y,z'
sg_symmetry(3, 132) = '-y,x,1/2+z'
sg_symmetry(4, 132) = 'y,-x,1/2+z'
sg_symmetry(5, 132) = '-x,y,1/2-z'
sg_symmetry(6, 132) = 'x,-y,1/2-z'
sg_symmetry(7, 132) = 'y,x,-z'
sg_symmetry(8, 132) = '-y,-x,-z'
sg_symmetry(9, 132) = '-x,-y,-z'
sg_symmetry(10, 132) = 'x,y,-z'
sg_symmetry(11, 132) = 'y,-x,1/2-z'
sg_symmetry(12, 132) = '-y,x,1/2-z'
sg_symmetry(13, 132) = 'x,-y,1/2+z'
sg_symmetry(14, 132) = '-x,y,1/2+z'
sg_symmetry(15, 132) = '-y,-x,z'
sg_symmetry(16, 132) = 'y,x,z'
sg_name(133) = 'P4(2)/nbc'
sg_patn(133) = 123
sg_symnum(133) = 16
sg_symmetry(1, 133) = 'x,y,z'
sg_symmetry(2, 133) = '-x,-y,z'
sg_symmetry(3, 133) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(4, 133) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(5, 133) = '-x,y,1/2-z'
sg_symmetry(6, 133) = 'x,-y,1/2-z'
sg_symmetry(7, 133) = '1/2+y,1/2+x,-z'
sg_symmetry(8, 133) = '1/2-y,1/2-x,-z'
sg_symmetry(9, 133) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(10, 133) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(11, 133) = 'y,-x,-z'
sg_symmetry(12, 133) = '-y,x,-z'
sg_symmetry(13, 133) = '1/2+x,1/2-y,z'
sg_symmetry(14, 133) = '1/2-x,1/2+y,z'
sg_symmetry(15, 133) = '-y,-x,1/2+z'
sg_symmetry(16, 133) = 'y,x,1/2+z'
sg_name(134) = 'P4(2)/nnm'
sg_patn(134) = 123
sg_symnum(134) = 16
sg_symmetry(1, 134) = 'x,y,z'
sg_symmetry(2, 134) = '-x,-y,z'
sg_symmetry(3, 134) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(4, 134) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(5, 134) = '-x,y,-z'
sg_symmetry(6, 134) = 'x,-y,-z'
sg_symmetry(7, 134) = '1/2+y,1/2+x,1/2-z'
sg_symmetry(8, 134) = '1/2-y,1/2-x,1/2-z'
sg_symmetry(9, 134) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(10, 134) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(11, 134) = 'y,-x,-z'
sg_symmetry(12, 134) = '-y,x,-z'
sg_symmetry(13, 134) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(14, 134) = '1/2-x,1/2+y,1/2+z'
sg_symmetry(15, 134) = '-y,-x,z'
sg_symmetry(16, 134) = 'y,x,z'
sg_name(135) = 'P4(2)/mbc'
sg_patn(135) = 123
sg_symnum(135) = 16
sg_symmetry(1, 135) = 'x,y,z'
sg_symmetry(2, 135) = '-x,-y,z'
sg_symmetry(3, 135) = '-y,x,1/2+z'
sg_symmetry(4, 135) = 'y,-x,1/2+z'
sg_symmetry(5, 135) = '1/2-x,1/2+y,-z'
sg_symmetry(6, 135) = '1/2+x,1/2-y,-z'
sg_symmetry(7, 135) = '1/2+y,1/2+x,1/2-z'
sg_symmetry(8, 135) = '1/2-y,1/2-x,1/2-z'
sg_symmetry(9, 135) = '-x,-y,-z'
sg_symmetry(10, 135) = 'x,y,-z'
sg_symmetry(11, 135) = 'y,-x,1/2-z'
sg_symmetry(12, 135) = '-y,x,1/2-z'
sg_symmetry(13, 135) = '1/2+x,1/2-y,z'
sg_symmetry(14, 135) = '1/2-x,1/2+y,z'
sg_symmetry(15, 135) = '1/2-y,1/2-x,1/2+z'
sg_symmetry(16, 135) = '1/2+y,1/2+x,1/2+z'
sg_name(136) = 'P4(2)/mnm'
sg_patn(136) = 123
sg_symnum(136) = 16
sg_symmetry(1, 136) = 'x,y,z'
sg_symmetry(2, 136) = '-x,-y,z'
sg_symmetry(3, 136) = '1/2-y,x+1/2,z+1/2'
sg_symmetry(4, 136) = 'y+1/2,1/2-x,z+1/2'
sg_symmetry(5, 136) = '1/2-x,y+1/2,1/2-z'
sg_symmetry(6, 136) = 'x+1/2,1/2-y,1/2-z'
sg_symmetry(7, 136) = 'y,x,-z'
sg_symmetry(8, 136) = '-y,-x,-z'
sg_symmetry(9, 136) = '-x,-y,-z'
sg_symmetry(10, 136) = 'x,y,-z'
sg_symmetry(11, 136) = 'y+1/2,1/2-x,1/2-z'
sg_symmetry(12, 136) = '1/2-y,x+1/2,1/2-z'
sg_symmetry(13, 136) = 'x+1/2,1/2-y,z+1/2'
sg_symmetry(14, 136) = '1/2-x,y+1/2,z+1/2'
sg_symmetry(15, 136) = '-y,-x,z'
sg_symmetry(16, 136) = 'y,x,z'
sg_name(137) = 'P4(2)/nmc'
sg_patn(137) = 123
sg_symnum(137) = 16
sg_symmetry(1, 137) = 'x,y,z'
sg_symmetry(2, 137) = '-x,-y,z'
sg_symmetry(3, 137) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(4, 137) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(5, 137) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(6, 137) = '1/2+x,1/2-y,1/2-z'
sg_symmetry(7, 137) = 'y,x,-z'
sg_symmetry(8, 137) = '-y,-x,-z'
sg_symmetry(9, 137) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(10, 137) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(11, 137) = 'y,-x,-z'
sg_symmetry(12, 137) = '-y,x,-z'
sg_symmetry(13, 137) = 'x,-y,z'
sg_symmetry(14, 137) = '-x,y,z'
sg_symmetry(15, 137) = '1/2-y,1/2-x,1/2+z'
sg_symmetry(16, 137) = '1/2+y,1/2+x,1/2+z'
sg_name(138) = 'P4(2)/ncm'
sg_patn(138) = 123
sg_symnum(138) = 16
sg_symmetry(1, 138) = 'x,y,z'
sg_symmetry(2, 138) = '-x,-y,z'
sg_symmetry(3, 138) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(4, 138) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(5, 138) = '1/2-x,1/2+y,-z'
sg_symmetry(6, 138) = '1/2+x,1/2-y,-z'
sg_symmetry(7, 138) = 'y,x,1/2-z'
sg_symmetry(8, 138) = '-y,-x,1/2-z'
sg_symmetry(9, 138) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(10, 138) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(11, 138) = 'y,-x,-z'
sg_symmetry(12, 138) = '-y,x,-z'
sg_symmetry(13, 138) = 'x,-y,1/2+z'
sg_symmetry(14, 138) = '-x,y,1/2+z'
sg_symmetry(15, 138) = '1/2-y,1/2-x,z'
sg_symmetry(16, 138) = '1/2+y,1/2+x,z'
sg_name(139) = 'I4/mmm'
sg_patn(139) = 139
sg_symnum(139) = 32
sg_symmetry(1, 139) = 'x,y,z'
sg_symmetry(2, 139) = '-x,-y,z'
sg_symmetry(3, 139) = '-y,x,z'
sg_symmetry(4, 139) = 'y,-x,z'
sg_symmetry(5, 139) = '-x,y,-z'
sg_symmetry(6, 139) = 'x,-y,-z'
sg_symmetry(7, 139) = 'y,x,-z'
sg_symmetry(8, 139) = '-y,-x,-z'
sg_symmetry(9, 139) = '-x,-y,-z'
sg_symmetry(10, 139) = 'x,y,-z'
sg_symmetry(11, 139) = 'y,-x,-z'
sg_symmetry(12, 139) = '-y,x,-z'
sg_symmetry(13, 139) = 'x,-y,z'
sg_symmetry(14, 139) = '-x,y,z'
sg_symmetry(15, 139) = '-y,-x,z'
sg_symmetry(16, 139) = 'y,x,z'
sg_symmetry(17, 139) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(18, 139) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(19, 139) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(20, 139) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(21, 139) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(22, 139) = '1/2+x,1/2-y,1/2-z'
sg_symmetry(23, 139) = '1/2+y,1/2+x,1/2-z'
sg_symmetry(24, 139) = '1/2-y,1/2-x,1/2-z'
sg_symmetry(25, 139) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(26, 139) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(27, 139) = '1/2+y,1/2-x,1/2-z'
sg_symmetry(28, 139) = '1/2-y,1/2+x,1/2-z'
sg_symmetry(29, 139) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(30, 139) = '1/2-x,1/2+y,1/2+z'
sg_symmetry(31, 139) = '1/2-y,1/2-x,1/2+z'
sg_symmetry(32, 139) = '1/2+y,1/2+x,1/2+z'
sg_name(140) = 'I4/mcm'
sg_patn(140) = 139
sg_symnum(140) = 32
sg_symmetry(1, 140) = 'x,y,z'
sg_symmetry(2, 140) = '-x,-y,z'
sg_symmetry(3, 140) = '-y,x,z'
sg_symmetry(4, 140) = 'y,-x,z'
sg_symmetry(5, 140) = '-x,y,1/2-z'
sg_symmetry(6, 140) = 'x,-y,1/2-z'
sg_symmetry(7, 140) = 'y,x,1/2-z'
sg_symmetry(8, 140) = '-y,-x,1/2-z'
sg_symmetry(9, 140) = '-x,-y,-z'
sg_symmetry(10, 140) = 'x,y,-z'
sg_symmetry(11, 140) = 'y,-x,-z'
sg_symmetry(12, 140) = '-y,x,-z'
sg_symmetry(13, 140) = 'x,-y,1/2+z'
sg_symmetry(14, 140) = '-x,y,1/2+z'
sg_symmetry(15, 140) = '-y,-x,1/2+z'
sg_symmetry(16, 140) = 'y,x,1/2+z'
sg_symmetry(17, 140) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(18, 140) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(19, 140) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(20, 140) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(21, 140) = '1/2-x,1/2+y,-z'
sg_symmetry(22, 140) = '1/2+x,1/2-y,-z'
sg_symmetry(23, 140) = '1/2+y,1/2+x,-z'
sg_symmetry(24, 140) = '1/2-y,1/2-x,-z'
sg_symmetry(25, 140) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(26, 140) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(27, 140) = '1/2+y,1/2-x,1/2-z'
sg_symmetry(28, 140) = '1/2-y,1/2+x,1/2-z'
sg_symmetry(29, 140) = '1/2+x,1/2-y,z'
sg_symmetry(30, 140) = '1/2-x,1/2+y,z'
sg_symmetry(31, 140) = '1/2-y,1/2-x,z'
sg_symmetry(32, 140) = '1/2+y,1/2+x,z'
sg_name(141) = 'I4(1)/amd'
sg_patn(141) = 139
sg_symnum(141) = 32
sg_symmetry(1, 141) = 'x,y,z'
sg_symmetry(2, 141) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(3, 141) = '-y,1/2+x,1/4+z'
sg_symmetry(4, 141) = '1/2+y,-x,3/4+z'
sg_symmetry(5, 141) = '1/2-x,y,3/4-z'
sg_symmetry(6, 141) = 'x,1/2-y,1/4-z'
sg_symmetry(7, 141) = '1/2+y,1/2+x,1/2-z'
sg_symmetry(8, 141) = '-y,-x,-z'
sg_symmetry(9, 141) = '-x,1/2-y,1/4-z'
sg_symmetry(10, 141) = '1/2+x,y,3/4-z'
sg_symmetry(11, 141) = 'y,-x,-z'
sg_symmetry(12, 141) = '1/2-y,1/2+x,1/2-z'
sg_symmetry(13, 141) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(14, 141) = '-x,y,z'
sg_symmetry(15, 141) = '1/2-y,-x,3/4+z'
sg_symmetry(16, 141) = 'y,1/2+x,1/4+z'
sg_symmetry(17, 141) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(18, 141) = '-x,-y,z'
sg_symmetry(19, 141) = '1/2-y,x,3/4+z'
sg_symmetry(20, 141) = 'y,1/2-x,1/4+z'
sg_symmetry(21, 141) = '-x,1/2+y,1/4-z'
sg_symmetry(22, 141) = '1/2+x,-y,3/4-z'
sg_symmetry(23, 141) = 'y,x,-z'
sg_symmetry(24, 141) = '1/2-y,1/2-x,1/2-z'
sg_symmetry(25, 141) = '1/2-x,-y,3/4-z'
sg_symmetry(26, 141) = 'x,1/2+y,1/4-z'
sg_symmetry(27, 141) = '1/2+y,1/2-x,1/2-z'
sg_symmetry(28, 141) = '-y,x,-z'
sg_symmetry(29, 141) = 'x,-y,z'
sg_symmetry(30, 141) = '1/2-x,1/2+y,1/2+z'
sg_symmetry(31, 141) = '-y,1/2-x,1/4+z'
sg_symmetry(32, 141) = '1/2+y,x,3/4+z'
sg_name(142) = 'I4(1)/acd'
sg_patn(142) = 139
sg_symnum(142) = 32
sg_symmetry(1, 142) = 'x,y,z'
sg_symmetry(2, 142) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(3, 142) = '-y,1/2+x,1/4+z'
sg_symmetry(4, 142) = '1/2+y,-x,3/4+z'
sg_symmetry(5, 142) = '1/2-x,y,1/4-z'
sg_symmetry(6, 142) = 'x,1/2-y,3/4-z'
sg_symmetry(7, 142) = '1/2+y,1/2+x,-z'
sg_symmetry(8, 142) = '-y,-x,1/2-z'
sg_symmetry(9, 142) = '-x,1/2-y,1/4-z'
sg_symmetry(10, 142) = '1/2+x,y,3/4-z'
sg_symmetry(11, 142) = 'y,-x,-z'
sg_symmetry(12, 142) = '1/2-y,1/2+x,1/2-z'
sg_symmetry(13, 142) = '1/2+x,1/2-y,z'
sg_symmetry(14, 142) = '-x,y,1/2+z'
sg_symmetry(15, 142) = '1/2-y,-x,1/4+z'
sg_symmetry(16, 142) = 'y,1/2+x,3/4+z'
sg_symmetry(17, 142) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(18, 142) = '-x,-y,z'
sg_symmetry(19, 142) = '1/2-y,x,3/4+z'
sg_symmetry(20, 142) = 'y,1/2-x,1/4+z'
sg_symmetry(21, 142) = '-x,1/2+y,3/4-z'
sg_symmetry(22, 142) = '1/2+x,-y,1/4-z'
sg_symmetry(23, 142) = 'y,x,1/2-z'
sg_symmetry(24, 142) = '1/2-y,1/2-x,-z'
sg_symmetry(25, 142) = '1/2-x,-y,3/4-z'
sg_symmetry(26, 142) = 'x,1/2+y,1/4-z'
sg_symmetry(27, 142) = '1/2+y,1/2-x,1/2-z'
sg_symmetry(28, 142) = '-y,x,-z'
sg_symmetry(29, 142) = 'x,-y,1/2+z'
sg_symmetry(30, 142) = '1/2-x,1/2+y,z'
sg_symmetry(31, 142) = '-y,1/2-x,3/4+z'
sg_symmetry(32, 142) = '1/2+y,x,1/4+z'
sg_name(143) = 'P3'
sg_patn(143) = 147
sg_symnum(143) = 3
sg_symmetry(1, 143) = 'x,y,z'
sg_symmetry(2, 143) = '-y,x-y,z'
sg_symmetry(3, 143) = 'y-x,-x,z'
sg_name(144) = 'P3(1)'
sg_patn(144) = 147
sg_symnum(144) = 3
sg_symmetry(1, 144) = 'x,y,z'
sg_symmetry(2, 144) = '-y,x-y,z+1/3'
sg_symmetry(3, 144) = 'y-x,-x,z+2/3'
sg_name(145) = 'P3(2)'
sg_patn(145) = 147
sg_symnum(145) = 3
sg_symmetry(1, 145) = 'x,y,z'
sg_symmetry(2, 145) = '-y,x-y,z+2/3'
sg_symmetry(3, 145) = 'y-x,-x,z+1/3'
sg_name(146) = 'R3'
sg_patn(146) = 148
sg_symnum(146) = 9
sg_symmetry(1, 146) = 'x,y,z'
sg_symmetry(2, 146) = '-y,x-y,z'
sg_symmetry(3, 146) = 'y-x,-x,z'
sg_symmetry(4, 146) = 'x+2/3,y+1/3,z+1/3'
sg_symmetry(5, 146) = '-y+2/3,x-y+1/3,z+1/3'
sg_symmetry(6, 146) = 'y-x+2/3,-x+1/3,z+1/3'
sg_symmetry(7, 146) = 'x+1/3,y+2/3,z+2/3'
sg_symmetry(8, 146) = '-y+1/3,x-y+2/3,z+2/3'
sg_symmetry(9, 146) = 'y-x+1/3,-x+2/3,z+2/3'
sg_name(147) = 'P-3'
sg_patn(147) = 147
sg_symnum(147) = 6
sg_symmetry(1, 147) = 'x,y,z'
sg_symmetry(2, 147) = '-y,x-y,z'
sg_symmetry(3, 147) = 'y-x,-x,z'
sg_symmetry(4, 147) = '-x,-y,-z'
sg_symmetry(5, 147) = 'y,y-x,-z'
sg_symmetry(6, 147) = 'x-y,x,-z'
sg_name(148) = 'R-3'
sg_patn(148) = 148
sg_symnum(148) = 18
sg_symmetry(1, 148) = 'x,y,z'
sg_symmetry(2, 148) = '-y,x-y,z'
sg_symmetry(3, 148) = 'y-x,-x,z'
sg_symmetry(4, 148) = '-x,-y,-z'
sg_symmetry(5, 148) = 'y,y-x,-z'
sg_symmetry(6, 148) = 'x-y,x,-z'
sg_symmetry(7, 148) = '2/3+x,1/3+y,1/3+z'
sg_symmetry(8, 148) = '2/3-y,1/3+x-y,1/3+z'
sg_symmetry(9, 148) = '2/3+y-x,1/3-x,1/3+z'
sg_symmetry(10, 148) = '2/3-x,1/3-y,1/3-z'
sg_symmetry(11, 148) = '2/3+y,1/3+y-x,1/3-z'
sg_symmetry(12, 148) = '2/3+x-y,1/3+x,1/3-z'
sg_symmetry(13, 148) = '1/3+x,2/3+y,2/3+z'
sg_symmetry(14, 148) = '1/3-y,2/3+x-y,2/3+z'
sg_symmetry(15, 148) = '1/3+y-x,2/3-x,2/3+z'
sg_symmetry(16, 148) = '1/3-x,2/3-y,2/3-z'
sg_symmetry(17, 148) = '1/3+y,2/3+y-x,2/3-z'
sg_symmetry(18, 148) = '1/3+x-y,2/3+x,2/3-z'
sg_name(149) = 'P312'
sg_patn(149) = 162
sg_symnum(149) = 6
sg_symmetry(1, 149) = 'x,y,z'
sg_symmetry(2, 149) = '-y,x-y,z'
sg_symmetry(3, 149) = 'y-x,-x,z'
sg_symmetry(4, 149) = '-y,-x,-z'
sg_symmetry(5, 149) = 'y-x,y,-z'
sg_symmetry(6, 149) = 'x,x-y,-z'
sg_name(150) = 'P321'
sg_patn(150) = 164
sg_symnum(150) = 6
sg_symmetry(1, 150) = 'x,y,z'
sg_symmetry(2, 150) = '-y,x-y,z'
sg_symmetry(3, 150) = 'y-x,-x,z'
sg_symmetry(4, 150) = 'y,x,-z'
sg_symmetry(5, 150) = 'x-y,-y,-z'
sg_symmetry(6, 150) = '-x,y-x,-z'
sg_name(151) = 'P3(1)12'
sg_patn(151) = 162
sg_symnum(151) = 6
sg_symmetry(1, 151) = 'x,y,z'
sg_symmetry(2, 151) = '-y,x-y,1/3+z'
sg_symmetry(3, 151) = 'y-x,-x,2/3+z'
sg_symmetry(4, 151) = '-y,-x,2/3-z'
sg_symmetry(5, 151) = 'y-x,y,1/3-z'
sg_symmetry(6, 151) = 'x,x-y,-z'
sg_name(152) = 'P3(1)21'
sg_patn(152) = 164
sg_symnum(152) = 6
sg_symmetry(1, 152) = 'x,y,z'
sg_symmetry(2, 152) = '-y,x-y,z+1/3'
sg_symmetry(3, 152) = 'y-x,-x,z+2/3'
sg_symmetry(4, 152) = 'y,x,-z'
sg_symmetry(5, 152) = 'x-y,-y,2/3-z'
sg_symmetry(6, 152) = '-x,y-x,1/3-z'
sg_name(153) = 'P3(2)12'
sg_patn(153) = 162
sg_symnum(153) = 6
sg_symmetry(1, 153) = 'x,y,z'
sg_symmetry(2, 153) = '-y,x-y,2/3+z'
sg_symmetry(3, 153) = 'y-x,-x,1/3+z'
sg_symmetry(4, 153) = '-y,-x,1/3-z'
sg_symmetry(5, 153) = 'y-x,y,2/3-z'
sg_symmetry(6, 153) = 'x,x-y,-z'
sg_name(154) = 'P3(2)21'
sg_patn(154) = 164
sg_symnum(154) = 6
sg_symmetry(1, 154) = 'x,y,z'
sg_symmetry(2, 154) = '-y,x-y,z+2/3'
sg_symmetry(3, 154) = 'y-x,-x,z+1/3'
sg_symmetry(4, 154) = 'y,x,-z'
sg_symmetry(5, 154) = 'x-y,-y,1/3-z'
sg_symmetry(6, 154) = '-x,y-x,2/3-z'
sg_name(155) = 'R32'
sg_patn(155) = 166
sg_symnum(155) = 18
sg_symmetry(1, 155) = 'x,y,z'
sg_symmetry(2, 155) = '-y,x-y,z'
sg_symmetry(3, 155) = 'y-x,-x,z'
sg_symmetry(4, 155) = 'y,x,-z'
sg_symmetry(5, 155) = 'x-y,-y,-z'
sg_symmetry(6, 155) = '-x,y-x,-z'
sg_symmetry(7, 155) = '2/3+x,1/3+y,1/3+z'
sg_symmetry(8, 155) = '2/3-y,1/3+x-y,1/3+z'
sg_symmetry(9, 155) = '2/3+y-x,1/3-x,1/3+z'
sg_symmetry(10, 155) = '2/3+y,1/3+x,1/3-z'
sg_symmetry(11, 155) = '2/3+x-y,1/3-y,1/3-z'
sg_symmetry(12, 155) = '2/3-x,1/3+y-x,1/3-z'
sg_symmetry(13, 155) = '1/3+x,2/3+y,2/3+z'
sg_symmetry(14, 155) = '1/3-y,2/3+x-y,2/3+z'
sg_symmetry(15, 155) = '1/3+y-x,2/3-x,2/3+z'
sg_symmetry(16, 155) = '1/3+y,2/3+x,2/3-z'
sg_symmetry(17, 155) = '1/3+x-y,2/3-y,2/3-z'
sg_symmetry(18, 155) = '1/3-x,2/3+y-x,2/3-z'
sg_name(156) = 'P3m1'
sg_patn(156) = 164
sg_symnum(156) = 6
sg_symmetry(1, 156) = 'x,y,z'
sg_symmetry(2, 156) = '-y,x-y,z'
sg_symmetry(3, 156) = 'y-x,-x,z'
sg_symmetry(4, 156) = '-y,-x,z'
sg_symmetry(5, 156) = 'y-x,y,z'
sg_symmetry(6, 156) = 'x,x-y,z'
sg_name(157) = 'P31m'
sg_patn(157) = 162
sg_symnum(157) = 6
sg_symmetry(1, 157) = 'x,y,z'
sg_symmetry(2, 157) = '-y,x-y,z'
sg_symmetry(3, 157) = 'y-x,-x,z'
sg_symmetry(4, 157) = 'y,x,z'
sg_symmetry(5, 157) = 'x-y,-y,z'
sg_symmetry(6, 157) = '-x,y-x,z'
sg_name(158) = 'P3c1'
sg_patn(158) = 164
sg_symnum(158) = 6
sg_symmetry(1, 158) = 'x,y,z'
sg_symmetry(2, 158) = '-y,x-y,z'
sg_symmetry(3, 158) = 'y-x,-x,z'
sg_symmetry(4, 158) = '-y,-x,1/2+z'
sg_symmetry(5, 158) = 'y-x,y,1/2+z'
sg_symmetry(6, 158) = 'x,x-y,1/2+z'
sg_name(159) = 'P31c'
sg_patn(159) = 162
sg_symnum(159) = 6
sg_symmetry(1, 159) = 'x,y,z'
sg_symmetry(2, 159) = '-y,x-y,z'
sg_symmetry(3, 159) = 'y-x,-x,z'
sg_symmetry(4, 159) = 'y,x,1/2+z'
sg_symmetry(5, 159) = 'x-y,-y,1/2+z'
sg_symmetry(6, 159) = '-x,y-x,1/2+z'
sg_name(160) = 'R3m'
sg_patn(160) = 166
sg_symnum(160) = 18
sg_symmetry(1, 160) = 'x,y,z'
sg_symmetry(2, 160) = '-y,x-y,z'
sg_symmetry(3, 160) = 'y-x,-x,z'
sg_symmetry(4, 160) = '-y,-x,z'
sg_symmetry(5, 160) = 'y-x,y,z'
sg_symmetry(6, 160) = 'x,x-y,z'
sg_symmetry(7, 160) = '2/3+x,1/3+y,1/3+z'
sg_symmetry(8, 160) = '2/3-y,1/3+x-y,1/3+z'
sg_symmetry(9, 160) = '2/3+y-x,1/3-x,1/3+z'
sg_symmetry(10, 160) = '2/3-y,1/3-x,1/3+z'
sg_symmetry(11, 160) = '2/3+y-x,1/3+y,1/3+z'
sg_symmetry(12, 160) = '2/3+x,1/3+x-y,1/3+z'
sg_symmetry(13, 160) = '1/3+x,2/3+y,2/3+z'
sg_symmetry(14, 160) = '1/3-y,2/3+x-y,2/3+z'
sg_symmetry(15, 160) = '1/3+y-x,2/3-x,2/3+z'
sg_symmetry(16, 160) = '1/3-y,2/3-x,2/3+z'
sg_symmetry(17, 160) = '1/3+y-x,2/3+y,2/3+z'
sg_symmetry(18, 160) = '1/3+x,2/3+x-y,2/3+z'
sg_name(161) = 'R3c'
sg_patn(161) = 166
sg_symnum(161) = 18
sg_symmetry(1, 161) = 'x,y,z'
sg_symmetry(2, 161) = '-y,x-y,z'
sg_symmetry(3, 161) = 'y-x,-x,z'
sg_symmetry(4, 161) = '-y,-x,1/2+z'
sg_symmetry(5, 161) = 'y-x,y,1/2+z'
sg_symmetry(6, 161) = 'x,x-y,1/2+z'
sg_symmetry(7, 161) = '2/3+x,1/3+y,1/3+z'
sg_symmetry(8, 161) = '2/3-y,1/3+x-y,1/3+z'
sg_symmetry(9, 161) = '2/3+y-x,1/3-x,1/3+z'
sg_symmetry(10, 161) = '2/3-y,1/3-x,5/6+z'
sg_symmetry(11, 161) = '2/3+y-x,1/3+y,5/6+z'
sg_symmetry(12, 161) = '2/3+x,1/3+x-y,5/6+z'
sg_symmetry(13, 161) = '1/3+x,2/3+y,2/3+z'
sg_symmetry(14, 161) = '1/3-y,2/3+x-y,2/3+z'
sg_symmetry(15, 161) = '1/3+y-x,2/3-x,2/3+z'
sg_symmetry(16, 161) = '1/3-y,2/3-x,1/6+z'
sg_symmetry(17, 161) = '1/3+y-x,2/3+y,1/6+z'
sg_symmetry(18, 161) = '1/3+x,2/3+x-y,1/6+z'
sg_name(162) = 'P-31m'
sg_patn(162) = 162
sg_symnum(162) = 12
sg_symmetry(1, 162) = 'x,y,z'
sg_symmetry(2, 162) = '-y,x-y,z'
sg_symmetry(3, 162) = 'y-x,-x,z'
sg_symmetry(4, 162) = '-y,-x,-z'
sg_symmetry(5, 162) = 'y-x,y,-z'
sg_symmetry(6, 162) = 'x,x-y,-z'
sg_symmetry(7, 162) = '-x,-y,-z'
sg_symmetry(8, 162) = 'y,y-x,-z'
sg_symmetry(9, 162) = 'x-y,x,-z'
sg_symmetry(10, 162) = 'y,x,z'
sg_symmetry(11, 162) = 'x-y,-y,z'
sg_symmetry(12, 162) = '-x,y-x,z'
sg_name(163) = 'P-31c'
sg_patn(163) = 162
sg_symnum(163) = 12
sg_symmetry(1, 163) = 'x,y,z'
sg_symmetry(2, 163) = '-y,x-y,z'
sg_symmetry(3, 163) = 'y-x,-x,z'
sg_symmetry(4, 163) = '-y,-x,1/2-z'
sg_symmetry(5, 163) = 'y-x,y,1/2-z'
sg_symmetry(6, 163) = 'x,x-y,1/2-z'
sg_symmetry(7, 163) = '-x,-y,-z'
sg_symmetry(8, 163) = 'y,y-x,-z'
sg_symmetry(9, 163) = 'x-y,x,-z'
sg_symmetry(10, 163) = 'y,x,1/2+z'
sg_symmetry(11, 163) = 'x-y,-y,1/2+z'
sg_symmetry(12, 163) = '-x,y-x,1/2+z'
sg_name(164) = 'P-3m1'
sg_patn(164) = 164
sg_symnum(164) = 12
sg_symmetry(1, 164) = 'x,y,z'
sg_symmetry(2, 164) = '-y,x-y,z'
sg_symmetry(3, 164) = 'y-x,-x,z'
sg_symmetry(4, 164) = 'y,x,-z'
sg_symmetry(5, 164) = 'x-y,-y,-z'
sg_symmetry(6, 164) = '-x,y-x,-z'
sg_symmetry(7, 164) = '-x,-y,-z'
sg_symmetry(8, 164) = 'y,y-x,-z'
sg_symmetry(9, 164) = 'x-y,x,-z'
sg_symmetry(10, 164) = '-y,-x,z'
sg_symmetry(11, 164) = 'y-x,y,z'
sg_symmetry(12, 164) = 'x,x-y,z'
sg_name(165) = 'P-3c1'
sg_patn(165) = 164
sg_symnum(165) = 12
sg_symmetry(1, 165) = 'x,y,z'
sg_symmetry(2, 165) = '-y,x-y,z'
sg_symmetry(3, 165) = 'y-x,-x,z'
sg_symmetry(4, 165) = 'y,x,1/2-z'
sg_symmetry(5, 165) = 'x-y,-y,1/2-z'
sg_symmetry(6, 165) = '-x,y-x,1/2-z'
sg_symmetry(7, 165) = '-x,-y,-z'
sg_symmetry(8, 165) = 'y,y-x,-z'
sg_symmetry(9, 165) = 'x-y,x,-z'
sg_symmetry(10, 165) = '-y,-x,1/2+z'
sg_symmetry(11, 165) = 'y-x,y,1/2+z'
sg_symmetry(12, 165) = 'x,x-y,1/2+z'
sg_name(166) = 'R-3m'
sg_patn(166) = 166
sg_symnum(166) = 36
sg_symmetry(1, 166) = 'x,y,z'
sg_symmetry(2, 166) = '-y,x-y,z'
sg_symmetry(3, 166) = 'y-x,-x,z'
sg_symmetry(4, 166) = 'y,x,-z'
sg_symmetry(5, 166) = 'x-y,-y,-z'
sg_symmetry(6, 166) = '-x,y-x,-z'
sg_symmetry(7, 166) = '-x,-y,-z'
sg_symmetry(8, 166) = 'y,y-x,-z'
sg_symmetry(9, 166) = 'x-y,x,-z'
sg_symmetry(10, 166) = '-y,-x,z'
sg_symmetry(11, 166) = 'y-x,y,z'
sg_symmetry(12, 166) = 'x,x-y,z'
sg_symmetry(13, 166) = '2/3+x,1/3+y,1/3+z'
sg_symmetry(14, 166) = '2/3-y,1/3+x-y,1/3+z'
sg_symmetry(15, 166) = '2/3+y-x,1/3-x,1/3+z'
sg_symmetry(16, 166) = '2/3+y,1/3+x,1/3-z'
sg_symmetry(17, 166) = '2/3+x-y,1/3-y,1/3-z'
sg_symmetry(18, 166) = '2/3-x,1/3+y-x,1/3-z'
sg_symmetry(19, 166) = '2/3-x,1/3-y,1/3-z'
sg_symmetry(20, 166) = '2/3+y,1/3+y-x,1/3-z'
sg_symmetry(21, 166) = '2/3+x-y,1/3+x,1/3-z'
sg_symmetry(22, 166) = '2/3-y,1/3-x,1/3+z'
sg_symmetry(23, 166) = '2/3+y-x,1/3+y,1/3+z'
sg_symmetry(24, 166) = '2/3+x,1/3+x-y,1/3+z'
sg_symmetry(25, 166) = '1/3+x,2/3+y,2/3+z'
sg_symmetry(26, 166) = '1/3-y,2/3+x-y,2/3+z'
sg_symmetry(27, 166) = '1/3+y-x,2/3-x,2/3+z'
sg_symmetry(28, 166) = '1/3+y,2/3+x,2/3-z'
sg_symmetry(29, 166) = '1/3+x-y,2/3-y,2/3-z'
sg_symmetry(30, 166) = '1/3-x,2/3+y-x,2/3-z'
sg_symmetry(31, 166) = '1/3-x,2/3-y,2/3-z'
sg_symmetry(32, 166) = '1/3+y,2/3+y-x,2/3-z'
sg_symmetry(33, 166) = '1/3x-y,2/3+x,2/3-z'
sg_symmetry(34, 166) = '1/3-y,2/3-x,2/3+z'
sg_symmetry(35, 166) = '1/3+y-x,2/3+y,2/3+z'
sg_symmetry(36, 166) = '1/3+x,2/3+x-y,2/3+z'
sg_name(167) = 'R-3c'
sg_patn(167) = 166
sg_symnum(167) = 36
sg_symmetry(1, 167) = 'x,y,z'
sg_symmetry(2, 167) = '-y,x-y,z'
sg_symmetry(3, 167) = 'y-x,-x,z'
sg_symmetry(4, 167) = 'y,x,1/2-z'
sg_symmetry(5, 167) = 'x-y,-y,1/2-z'
sg_symmetry(6, 167) = '-x,y-x,1/2-z'
sg_symmetry(7, 167) = '-x,-y,-z'
sg_symmetry(8, 167) = 'y,y-x,-z'
sg_symmetry(9, 167) = 'x-y,x,-z'
sg_symmetry(10, 167) = '-y,-x,1/2+z'
sg_symmetry(11, 167) = 'y-x,y,1/2+z'
sg_symmetry(12, 167) = 'x,x-y,1/2+z'
sg_symmetry(13, 167) = '2/3+x,1/3+y,1/3+z'
sg_symmetry(14, 167) = '2/3-y,1/3+x-y,1/3+z'
sg_symmetry(15, 167) = '2/3+y-x,1/3-x,1/3+z'
sg_symmetry(16, 167) = '2/3+y,1/3+x,5/6-z'
sg_symmetry(17, 167) = '2/3+x-y,1/3-y,5/6-z'
sg_symmetry(18, 167) = '2/3-x,1/3+y-x,5/6-z'
sg_symmetry(19, 167) = '2/3-x,1/3-y,1/3-z'
sg_symmetry(20, 167) = '2/3+y,1/3+y-x,1/3-z'
sg_symmetry(21, 167) = '2/3+x-y,1/3+x,1/3-z'
sg_symmetry(22, 167) = '2/3-y,1/3-x,5/6+z'
sg_symmetry(23, 167) = '2/3+y-x,1/3+y,5/6+z'
sg_symmetry(24, 167) = '2/3+x,1/3+x-y,5/6+z'
sg_symmetry(25, 167) = '1/3+x,2/3+y,2/3+z'
sg_symmetry(26, 167) = '1/3-y,2/3+x-y,2/3+z'
sg_symmetry(27, 167) = '1/3+y-x,2/3-x,2/3+z'
sg_symmetry(28, 167) = '1/3+y,2/3+x,1/6-z'
sg_symmetry(29, 167) = '1/3+x-y,2/3-y,1/6-z'
sg_symmetry(30, 167) = '1/3-x,2/3+y-x,1/6-z'
sg_symmetry(31, 167) = '1/3-x,2/3-y,2/3-z'
sg_symmetry(32, 167) = '1/3+y,2/3+y-x,2/3-z'
sg_symmetry(33, 167) = '1/3+x-y,2/3+x,2/3-z'
sg_symmetry(34, 167) = '1/3-y,2/3-x,1/6+z'
sg_symmetry(35, 167) = '1/3+y-x,2/3+y,1/6+z'
sg_symmetry(36, 167) = '1/3+x,2/3+x-y,1/6+z'
sg_name(168) = 'P6'
sg_patn(168) = 175
sg_symnum(168) = 6
sg_symmetry(1, 168) = 'x,y,z'
sg_symmetry(2, 168) = '-y,x-y,z'
sg_symmetry(3, 168) = 'y-x,-x,z'
sg_symmetry(4, 168) = '-x,-y,z'
sg_symmetry(5, 168) = 'y,y-x,z'
sg_symmetry(6, 168) = 'x-y,x,z'
sg_name(169) = 'P6(1)'
sg_patn(169) = 175
sg_symnum(169) = 6
sg_symmetry(1, 169) = 'x,y,z'
sg_symmetry(2, 169) = '-y,x-y,z+1/3'
sg_symmetry(3, 169) = 'y-x,-x,z+2/3'
sg_symmetry(4, 169) = '-x,-y,z+1/2'
sg_symmetry(5, 169) = 'y,y-x,z+5/6'
sg_symmetry(6, 169) = 'x-y,x,z+1/6'
sg_name(170) = 'P6(5)'
sg_patn(170) = 175
sg_symnum(170) = 6
sg_symmetry(1, 170) = 'x,y,z'
sg_symmetry(2, 170) = '-y,x-y,z+2/3'
sg_symmetry(3, 170) = 'y-x,-x,z+1/3'
sg_symmetry(4, 170) = '-x,-y,z+1/2'
sg_symmetry(5, 170) = 'y,y-x,z+1/6'
sg_symmetry(6, 170) = 'x-y,x,z+5/6'
sg_name(171) = 'P6(2)'
sg_patn(171) = 175
sg_symnum(171) = 6
sg_symmetry(1, 171) = 'x,y,z'
sg_symmetry(2, 171) = '-y,x-y,2/3+z'
sg_symmetry(3, 171) = 'y-x,-x,1/3+z'
sg_symmetry(4, 171) = '-x,-y,z'
sg_symmetry(5, 171) = 'y,y-x,2/3+z'
sg_symmetry(6, 171) = 'x-y,x,1/3+z'
sg_name(172) = 'P6(4)'
sg_patn(172) = 175
sg_symnum(172) = 6
sg_symmetry(1, 172) = 'x,y,z'
sg_symmetry(2, 172) = '-y,x-y,1/3+z'
sg_symmetry(3, 172) = 'y-x,-x,2/3+z'
sg_symmetry(4, 172) = '-x,-y,z'
sg_symmetry(5, 172) = 'y,y-x,1/3+z'
sg_symmetry(6, 172) = 'x-y,x,2/3+z'
sg_name(173) = 'P6(3)'
sg_patn(173) = 175
sg_symnum(173) = 6
sg_symmetry(1, 173) = 'x,y,z'
sg_symmetry(2, 173) = '-y,x-y,z'
sg_symmetry(3, 173) = 'y-x,-x,z'
sg_symmetry(4, 173) = '-x,-y,1/2+z'
sg_symmetry(5, 173) = 'y,y-x,1/2+z'
sg_symmetry(6, 173) = 'x-y,x,1/2+z'
sg_name(174) = 'P-6'
sg_patn(174) = 175
sg_symnum(174) = 6
sg_symmetry(1, 174) = 'x,y,z'
sg_symmetry(2, 174) = '-y,x-y,z'
sg_symmetry(3, 174) = 'y-x,-x,z'
sg_symmetry(4, 174) = 'x,y,-z'
sg_symmetry(5, 174) = '-y,x-y,-z'
sg_symmetry(6, 174) = 'y-x,-x,-z'
sg_name(175) = 'P6/m'
sg_patn(175) = 175
sg_symnum(175) = 12
sg_symmetry(1, 175) = 'x,y,z'
sg_symmetry(2, 175) = '-y,x-y,z'
sg_symmetry(3, 175) = 'y-x,-x,z'
sg_symmetry(4, 175) = '-x,-y,z'
sg_symmetry(5, 175) = 'y,y-x,z'
sg_symmetry(6, 175) = 'x-y,x,z'
sg_symmetry(7, 175) = '-x,-y,-z'
sg_symmetry(8, 175) = 'y,y-x,-z'
sg_symmetry(9, 175) = 'x-y,x,-z'
sg_symmetry(10, 175) = 'x,y,-z'
sg_symmetry(11, 175) = '-y,x-y,-z'
sg_symmetry(12, 175) = 'y-x,-x,-z'
sg_name(176) = 'P6(3)/m'
sg_patn(176) = 175
sg_symnum(176) = 12
sg_symmetry(1, 176) = 'x,y,z'
sg_symmetry(2, 176) = '-y,x-y,z'
sg_symmetry(3, 176) = 'y-x,-x,z'
sg_symmetry(4, 176) = '-x,-y,1/2+z'
sg_symmetry(5, 176) = 'y,y-x,1/2+z'
sg_symmetry(6, 176) = 'x-y,x,1/2+z'
sg_symmetry(7, 176) = '-x,-y,-z'
sg_symmetry(8, 176) = 'y,y-x,-z'
sg_symmetry(9, 176) = 'x-y,x,-z'
sg_symmetry(10, 176) = 'x,y,1/2-z'
sg_symmetry(11, 176) = '-y,x-y,1/2-z'
sg_symmetry(12, 176) = 'y-x,-x,1/2-z'
sg_name(177) = 'P622'
sg_patn(177) = 191
sg_symnum(177) = 12
sg_symmetry(1, 177) = 'x,y,z'
sg_symmetry(2, 177) = '-y,x-y,z'
sg_symmetry(3, 177) = 'y-x,-x,z'
sg_symmetry(4, 177) = '-x,-y,z'
sg_symmetry(5, 177) = 'y,y-x,z'
sg_symmetry(6, 177) = 'x-y,x,z'
sg_symmetry(7, 177) = 'y,x,-z'
sg_symmetry(8, 177) = 'x-y,-y,-z'
sg_symmetry(9, 177) = '-x,y-x,-z'
sg_symmetry(10, 177) = '-y,-x,-z'
sg_symmetry(11, 177) = 'y-x,y,-z'
sg_symmetry(12, 177) = 'x,x-y,-z'
sg_name(178) = 'P6(1)22'
sg_patn(178) = 191
sg_symnum(178) = 12
sg_symmetry(1, 178) = 'x,y,z'
sg_symmetry(2, 178) = '-y,x-y,1/3+z'
sg_symmetry(3, 178) = 'y-x,-x,2/3+z'
sg_symmetry(4, 178) = '-x,-y,1/2+z'
sg_symmetry(5, 178) = 'y,y-x,5/6+z'
sg_symmetry(6, 178) = 'x-y,x,1/6+z'
sg_symmetry(7, 178) = 'y,x,1/3-z'
sg_symmetry(8, 178) = 'x-y,-y,-z'
sg_symmetry(9, 178) = '-x,y-x,2/3-z'
sg_symmetry(10, 178) = '-y,-x,5/6-z'
sg_symmetry(11, 178) = 'y-x,y,1/2-z'
sg_symmetry(12, 178) = 'x,x-y,1/6-z'
sg_name(179) = 'P6(5)22'
sg_patn(179) = 191
sg_symnum(179) = 12
sg_symmetry(1, 179) = 'x,y,z'
sg_symmetry(2, 179) = '-y,x-y,2/3+z'
sg_symmetry(3, 179) = 'y-x,-x,1/3+z'
sg_symmetry(4, 179) = '-x,-y,1/2+z'
sg_symmetry(5, 179) = 'y,y-x,1/6+z'
sg_symmetry(6, 179) = 'x-y,x,5/6+z'
sg_symmetry(7, 179) = 'y,x,2/3-z'
sg_symmetry(8, 179) = 'x-y,-y,-z'
sg_symmetry(9, 179) = '-x,y-x,1/3-z'
sg_symmetry(10, 179) = '-y,-x,1/6-z'
sg_symmetry(11, 179) = 'y-x,y,1/2-z'
sg_symmetry(12, 179) = 'x,x-y,5/6-z'
sg_name(180) = 'P6(2)22'
sg_patn(180) = 191
sg_symnum(180) = 12
sg_symmetry(1, 180) = 'x,y,z'
sg_symmetry(2, 180) = '-y,x-y,2/3+z'
sg_symmetry(3, 180) = 'y-x,-x,1/3+z'
sg_symmetry(4, 180) = '-x,-y,z'
sg_symmetry(5, 180) = 'y,y-x,2/3+z'
sg_symmetry(6, 180) = 'x-y,x,1/3+z'
sg_symmetry(7, 180) = 'y,x,2/3-z'
sg_symmetry(8, 180) = 'x-y,-y,-z'
sg_symmetry(9, 180) = '-x,y-x,1/3-z'
sg_symmetry(10, 180) = '-y,-x,2/3-z'
sg_symmetry(11, 180) = 'y-x,y,-z'
sg_symmetry(12, 180) = 'x,x-y,1/3-z'
sg_name(181) = 'P6(4)22'
sg_patn(181) = 191
sg_symnum(181) = 12
sg_symmetry(1, 181) = 'x,y,z'
sg_symmetry(2, 181) = '-y,x-y,1/3+z'
sg_symmetry(3, 181) = 'y-x,-x,2/3+z'
sg_symmetry(4, 181) = '-x,-y,z'
sg_symmetry(5, 181) = 'y,y-x,1/3+z'
sg_symmetry(6, 181) = 'x-y,x,2/3+z'
sg_symmetry(7, 181) = 'y,x,1/3-z'
sg_symmetry(8, 181) = 'x-y,-y,-z'
sg_symmetry(9, 181) = '-x,y-x,2/3-z'
sg_symmetry(10, 181) = '-y,-x,1/3-z'
sg_symmetry(11, 181) = 'y-x,y,-z'
sg_symmetry(12, 181) = 'x,x-y,2/3-z'
sg_name(182) = 'P6(3)22'
sg_patn(182) = 191
sg_symnum(182) = 12
sg_symmetry(1, 182) = 'x,y,z'
sg_symmetry(2, 182) = '-y,x-y,z'
sg_symmetry(3, 182) = 'y-x,-x,z'
sg_symmetry(4, 182) = '-x,-y,1/2+z'
sg_symmetry(5, 182) = 'y,y-x,1/2+z'
sg_symmetry(6, 182) = 'x-y,x,1/2+z'
sg_symmetry(7, 182) = 'y,x,-z'
sg_symmetry(8, 182) = 'x-y,-y,-z'
sg_symmetry(9, 182) = '-x,y-x,-z'
sg_symmetry(10, 182) = '-y,-x,1/2-z'
sg_symmetry(11, 182) = 'y-x,y,1/2-z'
sg_symmetry(12, 182) = 'x,x-y,1/2-z'
sg_name(183) = 'P6mm'
sg_patn(183) = 191
sg_symnum(183) = 12
sg_symmetry(1, 183) = 'x,y,z'
sg_symmetry(2, 183) = '-y,x-y,z'
sg_symmetry(3, 183) = 'y-x,-x,z'
sg_symmetry(4, 183) = '-x,-y,z'
sg_symmetry(5, 183) = 'y,y-x,z'
sg_symmetry(6, 183) = 'x-y,x,z'
sg_symmetry(7, 183) = '-y,-x,z'
sg_symmetry(8, 183) = 'y-x,y,z'
sg_symmetry(9, 183) = 'x,x-y,z'
sg_symmetry(10, 183) = 'y,x,z'
sg_symmetry(11, 183) = 'x-y,-y,z'
sg_symmetry(12, 183) = '-x,y-x,z'
sg_name(184) = 'P6cc'
sg_patn(184) = 191
sg_symnum(184) = 12
sg_symmetry(1, 184) = 'x,y,z'
sg_symmetry(2, 184) = '-y,x-y,z'
sg_symmetry(3, 184) = 'y-x,-x,z'
sg_symmetry(4, 184) = '-x,-y,z'
sg_symmetry(5, 184) = 'y,y-x,z'
sg_symmetry(6, 184) = 'x-y,x,z'
sg_symmetry(7, 184) = '-y,-x,1/2+z'
sg_symmetry(8, 184) = 'y-x,y,1/2+z'
sg_symmetry(9, 184) = 'x,x-y,1/2+z'
sg_symmetry(10, 184) = 'y,x,1/2+z'
sg_symmetry(11, 184) = 'x-y,-y,1/2+z'
sg_symmetry(12, 184) = '-x,y-x,1/2+z'
sg_name(185) = 'P6(3)cm'
sg_patn(185) = 191
sg_symnum(185) = 12
sg_symmetry(1, 185) = 'x,y,z'
sg_symmetry(2, 185) = '-y,x-y,z'
sg_symmetry(3, 185) = 'y-x,-x,z'
sg_symmetry(4, 185) = '-x,-y,1/2+z'
sg_symmetry(5, 185) = 'y,y-x,1/2+z'
sg_symmetry(6, 185) = 'x-y,x,1/2+z'
sg_symmetry(7, 185) = '-y,-x,1/2+z'
sg_symmetry(8, 185) = 'y-x,y,1/2+z'
sg_symmetry(9, 185) = 'x,x-y,1/2+z'
sg_symmetry(10, 185) = 'y,x,z'
sg_symmetry(11, 185) = 'x-y,-y,z'
sg_symmetry(12, 185) = '-x,y-x,z'
sg_name(186) = 'P6(3)mc'
sg_patn(186) = 191
sg_symnum(186) = 12
sg_symmetry(1, 186) = 'x,y,z'
sg_symmetry(2, 186) = '-y,x-y,z'
sg_symmetry(3, 186) = 'y-x,-x,z'
sg_symmetry(4, 186) = '-x,-y,1/2+z'
sg_symmetry(5, 186) = 'y,y-x,1/2+z'
sg_symmetry(6, 186) = 'x-y,x,1/2+z'
sg_symmetry(7, 186) = '-y,-x,z'
sg_symmetry(8, 186) = 'y-x,y,z'
sg_symmetry(9, 186) = 'x,x-y,z'
sg_symmetry(10, 186) = 'y,x,1/2+z'
sg_symmetry(11, 186) = 'x-y,-y,1/2+z'
sg_symmetry(12, 186) = '-x,y-x,1/2+z'
sg_name(187) = 'P-6m2'
sg_patn(187) = 191
sg_symnum(187) = 12
sg_symmetry(1, 187) = 'x,y,z'
sg_symmetry(2, 187) = '-y,x-y,z'
sg_symmetry(3, 187) = 'y-x,-x,z'
sg_symmetry(4, 187) = 'x,y,-z'
sg_symmetry(5, 187) = '-y,x-y,-z'
sg_symmetry(6, 187) = 'y-x,-x,-z'
sg_symmetry(7, 187) = '-y,-x,z'
sg_symmetry(8, 187) = 'y-x,y,z'
sg_symmetry(9, 187) = 'x,x-y,z'
sg_symmetry(10, 187) = '-y,-x,-z'
sg_symmetry(11, 187) = 'y-x,y,-z'
sg_symmetry(12, 187) = 'x,x-y,-z'
sg_name(188) = 'P-6c2'
sg_patn(188) = 191
sg_symnum(188) = 12
sg_symmetry(1, 188) = 'x,y,z'
sg_symmetry(2, 188) = '-y,x-y,z'
sg_symmetry(3, 188) = 'y-x,-x,z'
sg_symmetry(4, 188) = 'x,y,1/2-z'
sg_symmetry(5, 188) = '-y,x-y,1/2-z'
sg_symmetry(6, 188) = 'y-x,-x,1/2-z'
sg_symmetry(7, 188) = '-y,-x,1/2+z'
sg_symmetry(8, 188) = 'y-x,y,1/2+z'
sg_symmetry(9, 188) = 'x,x-y,1/2+z'
sg_symmetry(10, 188) = '-y,-x,-z'
sg_symmetry(11, 188) = 'y-x,y,-z'
sg_symmetry(12, 188) = 'x,x-y,-z'
sg_name(189) = 'P-62m'
sg_patn(189) = 191
sg_symnum(189) = 12
sg_symmetry(1, 189) = 'x,y,z'
sg_symmetry(2, 189) = '-y,x-y,z'
sg_symmetry(3, 189) = 'y-x,-x,z'
sg_symmetry(4, 189) = 'x,y,-z'
sg_symmetry(5, 189) = '-y,x-y,-z'
sg_symmetry(6, 189) = 'y-x,-x,-z'
sg_symmetry(7, 189) = 'y,x,-z'
sg_symmetry(8, 189) = 'x-y,-y,-z'
sg_symmetry(9, 189) = '-x,y-x,-z'
sg_symmetry(10, 189) = 'y,x,z'
sg_symmetry(11, 189) = 'x-y,-y,z'
sg_symmetry(12, 189) = '-x,y-x,z'
sg_name(190) = 'P-62c'
sg_patn(190) = 191
sg_symnum(190) = 12
sg_symmetry(1, 190) = 'x,y,z'
sg_symmetry(2, 190) = '-y,x-y,z'
sg_symmetry(3, 190) = 'y-x,-x,z'
sg_symmetry(4, 190) = 'x,y,1/2-z'
sg_symmetry(5, 190) = '-y,x-y,1/2-z'
sg_symmetry(6, 190) = 'y-x,-x,1/2-z'
sg_symmetry(7, 190) = 'y,x,-z'
sg_symmetry(8, 190) = 'x-y,-y,-z'
sg_symmetry(9, 190) = '-x,y-x,-z'
sg_symmetry(10, 190) = 'y,x,1/2+z'
sg_symmetry(11, 190) = 'x-y,-y,1/2+z'
sg_symmetry(12, 190) = '-x,y-x,1/2+z'
sg_name(191) = 'P6/mmm'
sg_patn(191) = 191
sg_symnum(191) = 24
sg_symmetry(1, 191) = 'x,y,z'
sg_symmetry(2, 191) = '-y,x-y,z'
sg_symmetry(3, 191) = 'y-x,-x,z'
sg_symmetry(4, 191) = '-x,-y,z'
sg_symmetry(5, 191) = 'y,y-x,z'
sg_symmetry(6, 191) = 'x-y,x,z'
sg_symmetry(7, 191) = 'y,x,-z'
sg_symmetry(8, 191) = 'x-y,-y,-z'
sg_symmetry(9, 191) = '-x,y-x,-z'
sg_symmetry(10, 191) = '-y,-x,-z'
sg_symmetry(11, 191) = 'y-x,y,-z'
sg_symmetry(12, 191) = 'x,x-y,-z'
sg_symmetry(13, 191) = '-x,-y,-z'
sg_symmetry(14, 191) = 'y,y-x,-z'
sg_symmetry(15, 191) = 'x-y,x,-z'
sg_symmetry(16, 191) = 'x,y,-z'
sg_symmetry(17, 191) = 'y-x,-x,-z'
sg_symmetry(18, 191) = '-y,x-y,-z'
sg_symmetry(19, 191) = '-y,-x,z'
sg_symmetry(20, 191) = 'y-x,y,z'
sg_symmetry(21, 191) = 'x,x-y,z'
sg_symmetry(22, 191) = 'y,x,z'
sg_symmetry(23, 191) = 'x-y,-y,z'
sg_symmetry(24, 191) = '-x,y-x,z'
sg_name(192) = 'P6/mcc'
sg_patn(192) = 191
sg_symnum(192) = 24
sg_symmetry(1, 192) = 'x,y,z'
sg_symmetry(2, 192) = '-y,x-y,z'
sg_symmetry(3, 192) = 'y-x,-x,z'
sg_symmetry(4, 192) = '-x,-y,z'
sg_symmetry(5, 192) = 'y,y-x,z'
sg_symmetry(6, 192) = 'x-y,x,z'
sg_symmetry(7, 192) = 'y,x,1/2-z'
sg_symmetry(8, 192) = 'x-y,-y,1/2-z'
sg_symmetry(9, 192) = '-x,y-x,1/2-z'
sg_symmetry(10, 192) = '-y,-x,1/2-z'
sg_symmetry(11, 192) = 'y-x,y,1/2-z'
sg_symmetry(12, 192) = 'x,x-y,1/2-z'
sg_symmetry(13, 192) = '-x,-y,-z'
sg_symmetry(14, 192) = 'y,y-x,-z'
sg_symmetry(15, 192) = 'x-y,x,-z'
sg_symmetry(16, 192) = 'x,y,-z'
sg_symmetry(17, 192) = 'y-x,-x,-z'
sg_symmetry(18, 192) = '-y,x-y,-z'
sg_symmetry(19, 192) = '-y,-x,1/2+z'
sg_symmetry(20, 192) = 'y-x,y,1/2+z'
sg_symmetry(21, 192) = 'x,x-y,1/2+z'
sg_symmetry(22, 192) = 'y,x,1/2+z'
sg_symmetry(23, 192) = 'x-y,-y,1/2+z'
sg_symmetry(24, 192) = '-x,y-x,1/2+z'
sg_name(193) = 'P6(3)/mcm'
sg_patn(193) = 191
sg_symnum(193) = 24
sg_symmetry(1, 193) = 'x,y,z'
sg_symmetry(2, 193) = '-y,x-y,z'
sg_symmetry(3, 193) = 'y-x,-x,z'
sg_symmetry(4, 193) = '-x,-y,1/2+z'
sg_symmetry(5, 193) = 'y,y-x,1/2+z'
sg_symmetry(6, 193) = 'x-y,x,1/2+z'
sg_symmetry(7, 193) = 'y,x,1/2-z'
sg_symmetry(8, 193) = 'x-y,-y,1/2-z'
sg_symmetry(9, 193) = '-x,y-x,1/2-z'
sg_symmetry(10, 193) = '-y,-x,-z'
sg_symmetry(11, 193) = 'y-x,y,-z'
sg_symmetry(12, 193) = 'x,x-y,-z'
sg_symmetry(13, 193) = '-x,-y,-z'
sg_symmetry(14, 193) = 'y,y-x,-z'
sg_symmetry(15, 193) = 'x-y,x,-z'
sg_symmetry(16, 193) = 'x,y,1/2-z'
sg_symmetry(17, 193) = 'y-x,-x,1/2-z'
sg_symmetry(18, 193) = '-y,x-y,1/2-z'
sg_symmetry(19, 193) = '-y,-x,1/2+z'
sg_symmetry(20, 193) = 'y-x,y,1/2+z'
sg_symmetry(21, 193) = 'x,x-y,1/2+z'
sg_symmetry(22, 193) = 'y,x,z'
sg_symmetry(23, 193) = 'x-y,-y,z'
sg_symmetry(24, 193) = '-x,y-x,z'
sg_name(194) = 'P6(3)/mmc'
sg_patn(194) = 191
sg_symnum(194) = 24
sg_symmetry(1, 194) = 'x,y,z'
sg_symmetry(2, 194) = '-y,x-y,z'
sg_symmetry(3, 194) = 'y-x,-x,z'
sg_symmetry(4, 194) = '-x,-y,1/2+z'
sg_symmetry(5, 194) = 'y,y-x,1/2+z'
sg_symmetry(6, 194) = 'x-y,x,1/2+z'
sg_symmetry(7, 194) = 'y,x,-z'
sg_symmetry(8, 194) = 'x-y,-y,-z'
sg_symmetry(9, 194) = '-x,y-x,-z'
sg_symmetry(10, 194) = '-y,-x,1/2-z'
sg_symmetry(11, 194) = 'y-x,y,1/2-z'
sg_symmetry(12, 194) = 'x,x-y,1/2-z'
sg_symmetry(13, 194) = '-x,-y,-z'
sg_symmetry(14, 194) = 'y,y-x,-z'
sg_symmetry(15, 194) = 'x-y,x,-z'
sg_symmetry(16, 194) = 'x,y,1/2-z'
sg_symmetry(17, 194) = 'y-x,-x,1/2-z'
sg_symmetry(18, 194) = '-y,x-y,1/2-z'
sg_symmetry(19, 194) = '-y,-x,z'
sg_symmetry(20, 194) = 'y-x,y,z'
sg_symmetry(21, 194) = 'x,x-y,z'
sg_symmetry(22, 194) = 'y,x,1/2+z'
sg_symmetry(23, 194) = 'x-y,-y,1/2+z'
sg_symmetry(24, 194) = '-x,y-x,1/2+z'
sg_name(195) = 'P23'
sg_patn(195) = 200
sg_symnum(195) = 12
sg_symmetry(1, 195) = 'x,y,z'
sg_symmetry(2, 195) = '-x,-y,z'
sg_symmetry(3, 195) = '-x,y,-z'
sg_symmetry(4, 195) = 'x,-y,-z'
sg_symmetry(5, 195) = 'z,x,y'
sg_symmetry(6, 195) = 'z,-x,-y'
sg_symmetry(7, 195) = '-z,-x,y'
sg_symmetry(8, 195) = '-z,x,-y'
sg_symmetry(9, 195) = 'y,z,x'
sg_symmetry(10, 195) = '-y,z,-x'
sg_symmetry(11, 195) = 'y,-z,-x'
sg_symmetry(12, 195) = '-y,-z,x'
sg_name(196) = 'F23'
sg_patn(196) = 202
sg_symnum(196) = 48
sg_symmetry(1, 196) = 'x,y,z'
sg_symmetry(2, 196) = '-x,-y,z'
sg_symmetry(3, 196) = '-x,y,-z'
sg_symmetry(4, 196) = 'x,-y,-z'
sg_symmetry(5, 196) = 'z,x,y'
sg_symmetry(6, 196) = 'z,-x,-y'
sg_symmetry(7, 196) = '-z,-x,y'
sg_symmetry(8, 196) = '-z,x,-y'
sg_symmetry(9, 196) = 'y,z,x'
sg_symmetry(10, 196) = '-y,z,-x'
sg_symmetry(11, 196) = 'y,-z,-x'
sg_symmetry(12, 196) = '-y,-z,x'
sg_symmetry(13, 196) = 'x,1/2+y,1/2+z'
sg_symmetry(14, 196) = '-x,1/2-y,1/2+z'
sg_symmetry(15, 196) = '-x,1/2+y,1/2-z'
sg_symmetry(16, 196) = 'x,1/2-y,1/2-z'
sg_symmetry(17, 196) = 'z,1/2+x,1/2+y'
sg_symmetry(18, 196) = 'z,1/2-x,1/2-y'
sg_symmetry(19, 196) = '-z,1/2-x,1/2+y'
sg_symmetry(20, 196) = '-z,1/2+x,1/2-y'
sg_symmetry(21, 196) = 'y,1/2+z,1/2+x'
sg_symmetry(22, 196) = '-y,1/2+z,1/2-x'
sg_symmetry(23, 196) = 'y,1/2-z,1/2-x'
sg_symmetry(24, 196) = '-y,1/2-z,1/2+x'
sg_symmetry(25, 196) = '1/2+x,y,1/2+z'
sg_symmetry(26, 196) = '1/2-x,-y,1/2+z'
sg_symmetry(27, 196) = '1/2-x,y,1/2-z'
sg_symmetry(28, 196) = '1/2+x,-y,1/2-z'
sg_symmetry(29, 196) = '1/2+z,x,1/2+y'
sg_symmetry(30, 196) = '1/2+z,-x,1/2-y'
sg_symmetry(31, 196) = '1/2-z,-x,1/2+y'
sg_symmetry(32, 196) = '1/2-z,x,1/2-y'
sg_symmetry(33, 196) = '1/2+y,z,1/2+x'
sg_symmetry(34, 196) = '1/2-y,z,1/2-x'
sg_symmetry(35, 196) = '1/2+y,-z,1/2-x'
sg_symmetry(36, 196) = '1/2-y,-z,1/2+x'
sg_symmetry(37, 196) = '1/2+x,1/2+y,z'
sg_symmetry(38, 196) = '1/2-x,1/2-y,z'
sg_symmetry(39, 196) = '1/2-x,1/2+y,-z'
sg_symmetry(40, 196) = '1/2+x,1/2-y,-z'
sg_symmetry(41, 196) = '1/2+z,1/2+x,y'
sg_symmetry(42, 196) = '1/2+z,1/2-x,-y'
sg_symmetry(43, 196) = '1/2-z,1/2-x,y'
sg_symmetry(44, 196) = '1/2-z,1/2+x,-y'
sg_symmetry(45, 196) = '1/2+y,1/2+z,x'
sg_symmetry(46, 196) = '1/2-y,1/2+z,-x'
sg_symmetry(47, 196) = '1/2+y,1/2-z,-x'
sg_symmetry(48, 196) = '1/2-y,1/2-z,x'
sg_name(197) = 'I23'
sg_patn(197) = 204
sg_symnum(197) = 24
sg_symmetry(1, 197) = 'x,y,z'
sg_symmetry(2, 197) = '-x,-y,z'
sg_symmetry(3, 197) = '-x,y,-z'
sg_symmetry(4, 197) = 'x,-y,-z'
sg_symmetry(5, 197) = 'z,x,y'
sg_symmetry(6, 197) = 'z,-x,-y'
sg_symmetry(7, 197) = '-z,-x,y'
sg_symmetry(8, 197) = '-z,x,-y'
sg_symmetry(9, 197) = 'y,z,x'
sg_symmetry(10, 197) = '-y,z,-x'
sg_symmetry(11, 197) = 'y,-z,-x'
sg_symmetry(12, 197) = '-y,-z,x'
sg_symmetry(13, 197) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(14, 197) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(15, 197) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(16, 197) = '1/2+x,1/2-y,1/2-z'
sg_symmetry(17, 197) = '1/2+z,1/2+x,1/2+y'
sg_symmetry(18, 197) = '1/2+z,1/2-x,1/2-y'
sg_symmetry(19, 197) = '1/2-z,1/2-x,1/2+y'
sg_symmetry(20, 197) = '1/2-z,1/2+x,1/2-y'
sg_symmetry(21, 197) = '1/2+y,1/2+z,1/2+x'
sg_symmetry(22, 197) = '1/2-y,1/2+z,1/2-x'
sg_symmetry(23, 197) = '1/2+y,1/2-z,1/2-x'
sg_symmetry(24, 197) = '1/2-y,1/2-z,1/2+x'
sg_name(198) = 'P2(1)3'
sg_patn(198) = 200
sg_symnum(198) = 12
sg_symmetry(1, 198) = 'x,y,z'
sg_symmetry(2, 198) = '1/2-x,-y,1/2+z'
sg_symmetry(3, 198) = '-x,1/2+y,1/2-z'
sg_symmetry(4, 198) = '1/2+x,1/2-y,-z'
sg_symmetry(5, 198) = 'z,x,y'
sg_symmetry(6, 198) = '1/2+z,1/2-x,-y'
sg_symmetry(7, 198) = '1/2-z,-x,1/2+y'
sg_symmetry(8, 198) = '-z,1/2+x,1/2-y'
sg_symmetry(9, 198) = 'y,z,x'
sg_symmetry(10, 198) = '-y,1/2+z,1/2-x'
sg_symmetry(11, 198) = '1/2+y,1/2-z,-x'
sg_symmetry(12, 198) = '1/2-y,-z,1/2+x'
sg_name(199) = 'I2(1)3'
sg_patn(199) = 204
sg_symnum(199) = 24
sg_symmetry(1, 199) = 'x,y,z'
sg_symmetry(2, 199) = '1/2-x,-y,1/2+z'
sg_symmetry(3, 199) = '-x,1/2+y,1/2-z'
sg_symmetry(4, 199) = '1/2+x,1/2-y,-z'
sg_symmetry(5, 199) = 'z,x,y'
sg_symmetry(6, 199) = '1/2+z,1/2-x,-y'
sg_symmetry(7, 199) = '1/2-z,-x,1/2+y'
sg_symmetry(8, 199) = '-z,1/2+x,1/2-y'
sg_symmetry(9, 199) = 'y,z,x'
sg_symmetry(10, 199) = '-y,1/2+z,1/2-x'
sg_symmetry(11, 199) = '1/2+y,1/2-z,-x'
sg_symmetry(12, 199) = '1/2-y,-z,1/2+x'
sg_symmetry(13, 199) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(14, 199) = '-x,1/2-y,z'
sg_symmetry(15, 199) = '1/2-x,y,-z'
sg_symmetry(16, 199) = 'x,-y,1/2-z'
sg_symmetry(17, 199) = '1/2+z,1/2+x,1/2+y'
sg_symmetry(18, 199) = 'z,-x,1/2-y'
sg_symmetry(19, 199) = '-z,1/2-x,y'
sg_symmetry(20, 199) = '1/2-z,x,-y'
sg_symmetry(21, 199) = '1/2+y,1/2+z,1/2+x'
sg_symmetry(22, 199) = '1/2-y,z,-x'
sg_symmetry(23, 199) = 'y,-z,1/2-x'
sg_symmetry(24, 199) = '-y,1/2-z,x'
sg_name(200) = 'Pm-3'
sg_patn(200) = 200
sg_symnum(200) = 24
sg_symmetry(1, 200) = 'x,y,z'
sg_symmetry(2, 200) = '-x,-y,z'
sg_symmetry(3, 200) = '-x,y,-z'
sg_symmetry(4, 200) = 'x,-y,-z'
sg_symmetry(5, 200) = 'z,x,y'
sg_symmetry(6, 200) = 'z,-x,-y'
sg_symmetry(7, 200) = '-z,-x,y'
sg_symmetry(8, 200) = '-z,x,-y'
sg_symmetry(9, 200) = 'y,z,x'
sg_symmetry(10, 200) = '-y,z,-x'
sg_symmetry(11, 200) = 'y,-z,-x'
sg_symmetry(12, 200) = '-y,-z,x'
sg_symmetry(13, 200) = '-x,-y,-z'
sg_symmetry(14, 200) = 'x,y,-z'
sg_symmetry(15, 200) = 'x,-y,z'
sg_symmetry(16, 200) = '-x,y,z'
sg_symmetry(17, 200) = '-z,-x,-y'
sg_symmetry(18, 200) = '-z,x,y'
sg_symmetry(19, 200) = 'z,x,-y'
sg_symmetry(20, 200) = 'z,-x,y'
sg_symmetry(21, 200) = '-y,-z,-x'
sg_symmetry(22, 200) = 'y,-z,x'
sg_symmetry(23, 200) = '-y,z,x'
sg_symmetry(24, 200) = 'y,z,-x'
sg_name(201) = 'Pn-3'
sg_patn(201) = 200
sg_symnum(201) = 24
sg_symmetry(1, 201) = 'x,y,z'
sg_symmetry(2, 201) = '-x,-y,z'
sg_symmetry(3, 201) = '-x,y,-z'
sg_symmetry(4, 201) = 'x,-y,-z'
sg_symmetry(5, 201) = 'z,x,y'
sg_symmetry(6, 201) = 'z,-x,-y'
sg_symmetry(7, 201) = '-z,-x,y'
sg_symmetry(8, 201) = '-z,x,-y'
sg_symmetry(9, 201) = 'y,z,x'
sg_symmetry(10, 201) = '-y,z,-x'
sg_symmetry(11, 201) = 'y,-z,-x'
sg_symmetry(12, 201) = '-y,-z,x'
sg_symmetry(13, 201) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(14, 201) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(15, 201) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(16, 201) = '1/2-x,1/2+y,1/2+z'
sg_symmetry(17, 201) = '1/2-z,1/2-x,1/2-y'
sg_symmetry(18, 201) = '1/2-z,1/2+x,1/2+y'
sg_symmetry(19, 201) = '1/2+z,1/2+x,1/2-y'
sg_symmetry(20, 201) = '1/2+z,1/2-x,1/2+y'
sg_symmetry(21, 201) = '1/2-y,1/2-z,1/2-x'
sg_symmetry(22, 201) = '1/2+y,1/2-z,1/2+x'
sg_symmetry(23, 201) = '1/2-y,1/2+z,1/2+x'
sg_symmetry(24, 201) = '1/2+y,1/2+z,1/2-x'
sg_name(202) = 'Fm-3'
sg_patn(202) = 202
sg_symnum(202) = 96
sg_symmetry(1, 202) = 'x,y,z'
sg_symmetry(2, 202) = '-x,-y,z'
sg_symmetry(3, 202) = '-x,y,-z'
sg_symmetry(4, 202) = 'x,-y,-z'
sg_symmetry(5, 202) = 'z,x,y'
sg_symmetry(6, 202) = 'z,-x,-y'
sg_symmetry(7, 202) = '-z,-x,y'
sg_symmetry(8, 202) = '-z,x,-y'
sg_symmetry(9, 202) = 'y,z,x'
sg_symmetry(10, 202) = '-y,z,-x'
sg_symmetry(11, 202) = 'y,-z,-x'
sg_symmetry(12, 202) = '-y,-z,x'
sg_symmetry(13, 202) = '-x,-y,-z'
sg_symmetry(14, 202) = 'x,y,-z'
sg_symmetry(15, 202) = 'x,-y,z'
sg_symmetry(16, 202) = '-x,y,z'
sg_symmetry(17, 202) = '-z,-x,-y'
sg_symmetry(18, 202) = '-z,x,y'
sg_symmetry(19, 202) = 'z,x,-y'
sg_symmetry(20, 202) = 'z,-x,y'
sg_symmetry(21, 202) = '-y,-z,-x'
sg_symmetry(22, 202) = 'y,-z,x'
sg_symmetry(23, 202) = '-y,z,x'
sg_symmetry(24, 202) = 'y,z,-x'
sg_symmetry(25, 202) = 'x,1/2+y,1/2+z'
sg_symmetry(26, 202) = '-x,1/2-y,1/2+z'
sg_symmetry(27, 202) = '-x,1/2+y,1/2-z'
sg_symmetry(28, 202) = 'x,1/2-y,1/2-z'
sg_symmetry(29, 202) = 'z,1/2+x,1/2+y'
sg_symmetry(30, 202) = 'z,1/2-x,1/2-y'
sg_symmetry(31, 202) = '-z,1/2-x,1/2+y'
sg_symmetry(32, 202) = '-z,1/2+x,1/2-y'
sg_symmetry(33, 202) = 'y,1/2+z,1/2+x'
sg_symmetry(34, 202) = '-y,1/2+z,1/2-x'
sg_symmetry(35, 202) = 'y,1/2-z,1/2-x'
sg_symmetry(36, 202) = '-y,1/2-z,1/2+x'
sg_symmetry(37, 202) = '-x,1/2-y,1/2-z'
sg_symmetry(38, 202) = 'x,1/2+y,1/2-z'
sg_symmetry(39, 202) = 'x,1/2-y,1/2+z'
sg_symmetry(40, 202) = '-x,1/2+y,1/2+z'
sg_symmetry(41, 202) = '-z,1/2-x,1/2-y'
sg_symmetry(42, 202) = '-z,1/2+x,1/2+y'
sg_symmetry(43, 202) = 'z,1/2+x,1/2-y'
sg_symmetry(44, 202) = 'z,1/2-x,1/2+y'
sg_symmetry(45, 202) = '-y,1/2-z,1/2-x'
sg_symmetry(46, 202) = 'y,1/2-z,1/2+x'
sg_symmetry(47, 202) = '-y,1/2+z,1/2+x'
sg_symmetry(48, 202) = 'y,1/2+z,1/2-x'
sg_symmetry(49, 202) = '1/2+x,y,1/2+z'
sg_symmetry(50, 202) = '1/2-x,-y,1/2+z'
sg_symmetry(51, 202) = '1/2-x,y,1/2-z'
sg_symmetry(52, 202) = '1/2+x,-y,1/2-z'
sg_symmetry(53, 202) = '1/2+z,x,1/2+y'
sg_symmetry(54, 202) = '1/2+z,-x,1/2-y'
sg_symmetry(55, 202) = '1/2-z,-x,1/2+y'
sg_symmetry(56, 202) = '1/2-z,x,1/2-y'
sg_symmetry(57, 202) = '1/2+y,z,1/2+x'
sg_symmetry(58, 202) = '1/2-y,z,1/2-x'
sg_symmetry(59, 202) = '1/2+y,-z,1/2-x'
sg_symmetry(60, 202) = '1/2-y,-z,1/2+x'
sg_symmetry(61, 202) = '1/2-x,-y,1/2-z'
sg_symmetry(62, 202) = '1/2+x,y,1/2-z'
sg_symmetry(63, 202) = '1/2+x,-y,1/2+z'
sg_symmetry(64, 202) = '1/2-x,y,1/2+z'
sg_symmetry(65, 202) = '1/2-z,-x,1/2-y'
sg_symmetry(66, 202) = '1/2-z,x,1/2+y'
sg_symmetry(67, 202) = '1/2+z,x,1/2-y'
sg_symmetry(68, 202) = '1/2+z,-x,1/2+y'
sg_symmetry(69, 202) = '1/2-y,-z,1/2-x'
sg_symmetry(70, 202) = '1/2+y,-z,1/2+x'
sg_symmetry(71, 202) = '1/2-y,z,1/2+x'
sg_symmetry(72, 202) = '1/2+y,z,1/2-x'
sg_symmetry(73, 202) = '1/2+x,1/2+y,z'
sg_symmetry(74, 202) = '1/2-x,1/2-y,z'
sg_symmetry(75, 202) = '1/2-x,1/2+y,-z'
sg_symmetry(76, 202) = '1/2+x,1/2-y,-z'
sg_symmetry(77, 202) = '1/2+z,1/2+x,y'
sg_symmetry(78, 202) = '1/2+z,1/2-x,-y'
sg_symmetry(79, 202) = '1/2-z,1/2-x,y'
sg_symmetry(80, 202) = '1/2-z,1/2+x,-y'
sg_symmetry(81, 202) = '1/2+y,1/2+z,x'
sg_symmetry(82, 202) = '1/2-y,1/2+z,-x'
sg_symmetry(83, 202) = '1/2+y,1/2-z,-x'
sg_symmetry(84, 202) = '1/2-y,1/2-z,x'
sg_symmetry(85, 202) = '1/2-x,1/2-y,-z'
sg_symmetry(86, 202) = '1/2+x,1/2+y,-z'
sg_symmetry(87, 202) = '1/2+x,1/2-y,z'
sg_symmetry(88, 202) = '1/2-x,1/2+y,z'
sg_symmetry(89, 202) = '1/2-z,1/2-x,-y'
sg_symmetry(90, 202) = '1/2-z,1/2+x,y'
sg_symmetry(91, 202) = '1/2+z,1/2+x,-y'
sg_symmetry(92, 202) = '1/2+z,1/2-x,y'
sg_symmetry(93, 202) = '1/2-y,1/2-z,-x'
sg_symmetry(94, 202) = '1/2+y,1/2-z,x'
sg_symmetry(95, 202) = '1/2-y,1/2+z,x'
sg_symmetry(96, 202) = '1/2+y,1/2+z,-x'
sg_name(203) = 'Fd-3'
sg_patn(203) = 202
sg_symnum(203) = 96
sg_symmetry(1, 203) = 'x,y,z'
sg_symmetry(2, 203) = '-x,-y,z'
sg_symmetry(3, 203) = '-x,y,-z'
sg_symmetry(4, 203) = 'x,-y,-z'
sg_symmetry(5, 203) = 'z,x,y'
sg_symmetry(6, 203) = 'z,-x,-y'
sg_symmetry(7, 203) = '-z,-x,y'
sg_symmetry(8, 203) = '-z,x,-y'
sg_symmetry(9, 203) = 'y,z,x'
sg_symmetry(10, 203) = '-y,z,-x'
sg_symmetry(11, 203) = 'y,-z,-x'
sg_symmetry(12, 203) = '-y,-z,x'
sg_symmetry(13, 203) = '1/4-x,1/4-y,1/4-z'
sg_symmetry(14, 203) = '1/4+x,1/4+y,1/4-z'
sg_symmetry(15, 203) = '1/4+x,1/4-y,1/4+z'
sg_symmetry(16, 203) = '1/4-x,1/4+y,1/4+z'
sg_symmetry(17, 203) = '1/4-z,1/4-x,1/4-y'
sg_symmetry(18, 203) = '1/4-z,1/4+x,1/4+y'
sg_symmetry(19, 203) = '1/4+z,1/4+x,1/4-y'
sg_symmetry(20, 203) = '1/4+z,1/4-x,1/4+y'
sg_symmetry(21, 203) = '1/4-y,1/4-z,1/4-x'
sg_symmetry(22, 203) = '1/4+y,1/4-z,1/4+x'
sg_symmetry(23, 203) = '1/4-y,1/4+z,1/4+x'
sg_symmetry(24, 203) = '1/4+y,1/4+z,1/4-x'
sg_symmetry(25, 203) = 'x,1/2+y,1/2+z'
sg_symmetry(26, 203) = '-x,1/2-y,1/2+z'
sg_symmetry(27, 203) = '-x,1/2+y,1/2-z'
sg_symmetry(28, 203) = 'x,1/2-y,1/2-z'
sg_symmetry(29, 203) = 'z,1/2+x,1/2+y'
sg_symmetry(30, 203) = 'z,1/2-x,1/2-y'
sg_symmetry(31, 203) = '-z,1/2-x,1/2+y'
sg_symmetry(32, 203) = '-z,1/2+x,1/2-y'
sg_symmetry(33, 203) = 'y,1/2+z,1/2+x'
sg_symmetry(34, 203) = '-y,1/2+z,1/2-x'
sg_symmetry(35, 203) = 'y,1/2-z,1/2-x'
sg_symmetry(36, 203) = '-y,1/2-z,1/2+x'
sg_symmetry(37, 203) = '1/4-x,3/4-y,3/4-z'
sg_symmetry(38, 203) = '1/4+x,3/4+y,3/4-z'
sg_symmetry(39, 203) = '1/4+x,3/4-y,3/4+z'
sg_symmetry(40, 203) = '1/4-x,3/4+y,3/4+z'
sg_symmetry(41, 203) = '1/4-z,3/4-x,3/4-y'
sg_symmetry(42, 203) = '1/4-z,3/4+x,3/4+y'
sg_symmetry(43, 203) = '1/4+z,3/4+x,3/4-y'
sg_symmetry(44, 203) = '1/4+z,3/4-x,3/4+y'
sg_symmetry(45, 203) = '1/4-y,3/4-z,3/4-x'
sg_symmetry(46, 203) = '1/4+y,3/4-z,3/4+x'
sg_symmetry(47, 203) = '1/4-y,3/4+z,3/4+x'
sg_symmetry(48, 203) = '1/4+y,3/4+z,3/4-x'
sg_symmetry(49, 203) = '1/2+x,y,1/2+z'
sg_symmetry(50, 203) = '1/2-x,-y,1/2+z'
sg_symmetry(51, 203) = '1/2-x,y,1/2-z'
sg_symmetry(52, 203) = '1/2+x,-y,1/2-z'
sg_symmetry(53, 203) = '1/2+z,x,1/2+y'
sg_symmetry(54, 203) = '1/2+z,-x,1/2-y'
sg_symmetry(55, 203) = '1/2-z,-x,1/2+y'
sg_symmetry(56, 203) = '1/2-z,x,1/2-y'
sg_symmetry(57, 203) = '1/2+y,z,1/2+x'
sg_symmetry(58, 203) = '1/2-y,z,1/2-x'
sg_symmetry(59, 203) = '1/2+y,-z,1/2-x'
sg_symmetry(60, 203) = '1/2-y,-z,1/2+x'
sg_symmetry(61, 203) = '3/4-x,1/4-y,3/4-z'
sg_symmetry(62, 203) = '3/4+x,1/4+y,3/4-z'
sg_symmetry(63, 203) = '3/4+x,1/4-y,3/4+z'
sg_symmetry(64, 203) = '3/4-x,1/4+y,3/4+z'
sg_symmetry(65, 203) = '3/4-z,1/4-x,3/4-y'
sg_symmetry(66, 203) = '3/4-z,1/4+x,3/4+y'
sg_symmetry(67, 203) = '3/4+z,1/4+x,3/4-y'
sg_symmetry(68, 203) = '3/4+z,1/4-x,3/4+y'
sg_symmetry(69, 203) = '3/4-y,1/4-z,3/4-x'
sg_symmetry(70, 203) = '3/4+y,1/4-z,3/4+x'
sg_symmetry(71, 203) = '3/4-y,1/4+z,3/4+x'
sg_symmetry(72, 203) = '3/4+y,1/4+z,3/4-x'
sg_symmetry(73, 203) = '1/2+x,1/2+y,z'
sg_symmetry(74, 203) = '1/2-x,1/2-y,z'
sg_symmetry(75, 203) = '1/2-x,1/2+y,-z'
sg_symmetry(76, 203) = '1/2+x,1/2-y,-z'
sg_symmetry(77, 203) = '1/2+z,1/2+x,y'
sg_symmetry(78, 203) = '1/2+z,1/2-x,-y'
sg_symmetry(79, 203) = '1/2-z,1/2-x,y'
sg_symmetry(80, 203) = '1/2-z,1/2+x,-y'
sg_symmetry(81, 203) = '1/2+y,1/2+z,x'
sg_symmetry(82, 203) = '1/2-y,1/2+z,-x'
sg_symmetry(83, 203) = '1/2+y,1/2-z,-x'
sg_symmetry(84, 203) = '1/2-y,1/2-z,x'
sg_symmetry(85, 203) = '3/4-x,3/4-y,1/4-z'
sg_symmetry(86, 203) = '3/4+x,3/4+y,1/4-z'
sg_symmetry(87, 203) = '3/4+x,3/4-y,z+1/4'
sg_symmetry(88, 203) = '3/4-x,3/4+y,z+1/4'
sg_symmetry(89, 203) = '3/4-z,3/4-x,1/4-y'
sg_symmetry(90, 203) = '3/4-z,3/4+x,1/4+y'
sg_symmetry(91, 203) = '3/4+z,3/4+x,1/4-y'
sg_symmetry(92, 203) = '3/4+z,3/4-x,1/4+y'
sg_symmetry(93, 203) = '3/4-y,3/4-z,1/4-x'
sg_symmetry(94, 203) = '3/4+y,3/4-z,1/4+x'
sg_symmetry(95, 203) = '3/4-y,3/4+z,1/4+x'
sg_symmetry(96, 203) = '3/4+y,3/4+z,1/4-x'
sg_name(204) = 'Im-3'
sg_patn(204) = 204
sg_symnum(204) = 48
sg_symmetry(1, 204) = 'x,y,z'
sg_symmetry(2, 204) = '-x,-y,z'
sg_symmetry(3, 204) = '-x,y,-z'
sg_symmetry(4, 204) = 'x,-y,-z'
sg_symmetry(5, 204) = 'z,x,y'
sg_symmetry(6, 204) = 'z,-x,-y'
sg_symmetry(7, 204) = '-z,-x,y'
sg_symmetry(8, 204) = '-z,x,-y'
sg_symmetry(9, 204) = 'y,z,x'
sg_symmetry(10, 204) = '-y,z,-x'
sg_symmetry(11, 204) = 'y,-z,-x'
sg_symmetry(12, 204) = '-y,-z,x'
sg_symmetry(13, 204) = '-x,-y,-z'
sg_symmetry(14, 204) = 'x,y,-z'
sg_symmetry(15, 204) = 'x,-y,z'
sg_symmetry(16, 204) = '-x,y,z'
sg_symmetry(17, 204) = '-z,-x,-y'
sg_symmetry(18, 204) = '-z,x,y'
sg_symmetry(19, 204) = 'z,x,-y'
sg_symmetry(20, 204) = 'z,-x,y'
sg_symmetry(21, 204) = '-y,-z,-x'
sg_symmetry(22, 204) = 'y,-z,x'
sg_symmetry(23, 204) = '-y,z,x'
sg_symmetry(24, 204) = 'y,z,-x'
sg_symmetry(25, 204) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(26, 204) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(27, 204) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(28, 204) = '1/2+x,1/2-y,1/2-z'
sg_symmetry(29, 204) = '1/2+z,1/2+x,1/2+y'
sg_symmetry(30, 204) = '1/2+z,1/2-x,1/2-y'
sg_symmetry(31, 204) = '1/2-z,1/2-x,1/2+y'
sg_symmetry(32, 204) = '1/2-z,1/2+x,1/2-y'
sg_symmetry(33, 204) = '1/2+y,1/2+z,1/2+x'
sg_symmetry(34, 204) = '1/2-y,1/2+z,1/2-x'
sg_symmetry(35, 204) = '1/2+y,1/2-z,1/2-x'
sg_symmetry(36, 204) = '1/2-y,1/2-z,1/2+x'
sg_symmetry(37, 204) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(38, 204) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(39, 204) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(40, 204) = '1/2-x,1/2+y,1/2+z'
sg_symmetry(41, 204) = '1/2-z,1/2-x,1/2-y'
sg_symmetry(42, 204) = '1/2-z,1/2+x,1/2+y'
sg_symmetry(43, 204) = '1/2+z,1/2+x,1/2-y'
sg_symmetry(44, 204) = '1/2+z,1/2-x,1/2+y'
sg_symmetry(45, 204) = '1/2-y,1/2-z,1/2-x'
sg_symmetry(46, 204) = '1/2+y,1/2-z,1/2+x'
sg_symmetry(47, 204) = '1/2-y,1/2+z,1/2+x'
sg_symmetry(48, 204) = '1/2+y,1/2+z,1/2-x'
sg_name(205) = 'Pa-3'
sg_patn(205) = 200
sg_symnum(205) = 24
sg_symmetry(1, 205) = 'x,y,z'
sg_symmetry(2, 205) = '1/2-x,-y,1/2+z'
sg_symmetry(3, 205) = '-x,1/2+y,1/2-z'
sg_symmetry(4, 205) = '1/2+x,1/2-y,-z'
sg_symmetry(5, 205) = 'z,x,y'
sg_symmetry(6, 205) = '1/2+z,1/2-x,-y'
sg_symmetry(7, 205) = '1/2-z,-x,1/2+y'
sg_symmetry(8, 205) = '-z,1/2+x,1/2-y'
sg_symmetry(9, 205) = 'y,z,x'
sg_symmetry(10, 205) = '-y,1/2+z,1/2-x'
sg_symmetry(11, 205) = '1/2+y,1/2-z,-x'
sg_symmetry(12, 205) = '1/2-y,-z,1/2+x'
sg_symmetry(13, 205) = '-x,-y,-z'
sg_symmetry(14, 205) = '1/2+x,y,1/2-z'
sg_symmetry(15, 205) = 'x,1/2-y,1/2+z'
sg_symmetry(16, 205) = '1/2-x,1/2+y,z'
sg_symmetry(17, 205) = '-z,-x,-y'
sg_symmetry(18, 205) = '1/2-z,1/2+x,y'
sg_symmetry(19, 205) = '1/2+z,x,1/2-y'
sg_symmetry(20, 205) = 'z,1/2-x,1/2+y'
sg_symmetry(21, 205) = '-y,-z,-x'
sg_symmetry(22, 205) = 'y,1/2-z,1/2+x'
sg_symmetry(23, 205) = '1/2-y,1/2+z,x'
sg_symmetry(24, 205) = '1/2+y,z,1/2-x'
sg_name(206) = 'Ia-3'
sg_patn(206) = 204
sg_symnum(206) = 48
sg_symmetry(1, 206) = 'x,y,z'
sg_symmetry(2, 206) = '1/2-x,-y,1/2+z'
sg_symmetry(3, 206) = '-x,1/2+y,1/2-z'
sg_symmetry(4, 206) = '1/2+x,1/2-y,-z'
sg_symmetry(5, 206) = 'z,x,y'
sg_symmetry(6, 206) = '1/2+z,1/2-x,-y'
sg_symmetry(7, 206) = '1/2-z,-x,1/2+y'
sg_symmetry(8, 206) = '-z,1/2+x,1/2-y'
sg_symmetry(9, 206) = 'y,z,x'
sg_symmetry(10, 206) = '-y,1/2+z,1/2-x'
sg_symmetry(11, 206) = '1/2+y,1/2-z,-x'
sg_symmetry(12, 206) = '1/2-y,-z,1/2+x'
sg_symmetry(13, 206) = '-x,-y,-z'
sg_symmetry(14, 206) = '1/2+x,y,1/2-z'
sg_symmetry(15, 206) = 'x,1/2-y,1/2+z'
sg_symmetry(16, 206) = '1/2-x,1/2+y,z'
sg_symmetry(17, 206) = '-z,-x,-y'
sg_symmetry(18, 206) = '1/2-z,1/2+x,y'
sg_symmetry(19, 206) = '1/2+z,x,1/2-y'
sg_symmetry(20, 206) = 'z,1/2-x,1/2+y'
sg_symmetry(21, 206) = '-y,-z,-x'
sg_symmetry(22, 206) = 'y,1/2-z,1/2+x'
sg_symmetry(23, 206) = '1/2-y,1/2+z,x'
sg_symmetry(24, 206) = '1/2+y,z,1/2-x'
sg_symmetry(25, 206) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(26, 206) = '-x,1/2-y,z'
sg_symmetry(27, 206) = '1/2-x,+y,-z'
sg_symmetry(28, 206) = 'x,-y,1/2-z'
sg_symmetry(29, 206) = '1/2+z,1/2+x,1/2+y'
sg_symmetry(30, 206) = 'z,-x,1/2-y'
sg_symmetry(31, 206) = '-z,1/2-x,y'
sg_symmetry(32, 206) = '1/2-z,x,-y'
sg_symmetry(33, 206) = '1/2+y,1/2+z,1/2+x'
sg_symmetry(34, 206) = '1/2-y,z,-x'
sg_symmetry(35, 206) = 'y,-z,1/2-x'
sg_symmetry(36, 206) = '-y,1/2-z,x'
sg_symmetry(37, 206) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(38, 206) = 'x,1/2+y,-z'
sg_symmetry(39, 206) = '1/2+x,-y,z'
sg_symmetry(40, 206) = '-x,y,1/2+z'
sg_symmetry(41, 206) = '1/2-z,1/2-x,1/2-y'
sg_symmetry(42, 206) = '-z,x,1/2+y'
sg_symmetry(43, 206) = 'z,1/2+x,-y'
sg_symmetry(44, 206) = '1/2+z,-x,y'
sg_symmetry(45, 206) = '1/2-y,1/2-z,1/2-x'
sg_symmetry(46, 206) = '1/2+y,-z,x'
sg_symmetry(47, 206) = '-y,z,1/2+x'
sg_symmetry(48, 206) = 'y,1/2+z,-x'
sg_name(207) = 'P432'
sg_patn(207) = 221
sg_symnum(207) = 24
sg_symmetry(1, 207) = 'x,y,z'
sg_symmetry(2, 207) = '-x,-y,z'
sg_symmetry(3, 207) = '-x,y,-z'
sg_symmetry(4, 207) = 'x,-y,-z'
sg_symmetry(5, 207) = 'z,x,y'
sg_symmetry(6, 207) = 'z,-x,-y'
sg_symmetry(7, 207) = '-z,-x,y'
sg_symmetry(8, 207) = '-z,x,-y'
sg_symmetry(9, 207) = 'y,z,x'
sg_symmetry(10, 207) = '-y,z,-x'
sg_symmetry(11, 207) = 'y,-z,-x'
sg_symmetry(12, 207) = '-y,-z,x'
sg_symmetry(13, 207) = 'y,x,-z'
sg_symmetry(14, 207) = '-y,-x,-z'
sg_symmetry(15, 207) = 'y,-x,z'
sg_symmetry(16, 207) = '-y,x,z'
sg_symmetry(17, 207) = 'x,z,-y'
sg_symmetry(18, 207) = '-x,z,y'
sg_symmetry(19, 207) = '-x,-z,-y'
sg_symmetry(20, 207) = 'x,-z,y'
sg_symmetry(21, 207) = 'z,y,-x'
sg_symmetry(22, 207) = 'z,-y,x'
sg_symmetry(23, 207) = '-z,y,x'
sg_symmetry(24, 207) = '-z,-y,-x'
sg_name(208) = 'P4(2)32'
sg_patn(208) = 221
sg_symnum(208) = 24
sg_symmetry(1, 208) = 'x,y,z'
sg_symmetry(2, 208) = '-x,-y,z'
sg_symmetry(3, 208) = '-x,y,-z'
sg_symmetry(4, 208) = 'x,-y,-z'
sg_symmetry(5, 208) = 'z,x,y'
sg_symmetry(6, 208) = 'z,-x,-y'
sg_symmetry(7, 208) = '-z,-x,y'
sg_symmetry(8, 208) = '-z,x,-y'
sg_symmetry(9, 208) = 'y,z,x'
sg_symmetry(10, 208) = '-y,z,-x'
sg_symmetry(11, 208) = 'y,-z,-x'
sg_symmetry(12, 208) = '-y,-z,x'
sg_symmetry(13, 208) = '1/2+y,1/2+x,1/2-z'
sg_symmetry(14, 208) = '1/2-y,1/2-x,1/2-z'
sg_symmetry(15, 208) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(16, 208) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(17, 208) = '1/2+x,1/2+z,1/2-y'
sg_symmetry(18, 208) = '1/2-x,1/2+z,1/2+y'
sg_symmetry(19, 208) = '1/2-x,1/2-z,1/2-y'
sg_symmetry(20, 208) = '1/2+x,1/2-z,1/2+y'
sg_symmetry(21, 208) = '1/2+z,1/2+y,1/2-x'
sg_symmetry(22, 208) = '1/2+z,1/2-y,1/2+x'
sg_symmetry(23, 208) = '1/2-z,1/2+y,1/2+x'
sg_symmetry(24, 208) = '1/2-z,1/2-y,1/2-x'
sg_name(209) = 'F432'
sg_patn(209) = 225
sg_symnum(209) = 96
sg_symmetry(1, 209) = 'x,y,z'
sg_symmetry(2, 209) = '-x,-y,z'
sg_symmetry(3, 209) = '-x,y,-z'
sg_symmetry(4, 209) = 'x,-y,-z'
sg_symmetry(5, 209) = 'z,x,y'
sg_symmetry(6, 209) = 'z,-x,-y'
sg_symmetry(7, 209) = '-z,-x,y'
sg_symmetry(8, 209) = '-z,x,-y'
sg_symmetry(9, 209) = 'y,z,x'
sg_symmetry(10, 209) = '-y,z,-x'
sg_symmetry(11, 209) = 'y,-z,-x'
sg_symmetry(12, 209) = '-y,-z,x'
sg_symmetry(13, 209) = 'y,x,-z'
sg_symmetry(14, 209) = '-y,-x,-z'
sg_symmetry(15, 209) = 'y,-x,z'
sg_symmetry(16, 209) = '-y,x,z'
sg_symmetry(17, 209) = 'x,z,-y'
sg_symmetry(18, 209) = '-x,z,y'
sg_symmetry(19, 209) = '-x,-z,-y'
sg_symmetry(20, 209) = 'x,-z,y'
sg_symmetry(21, 209) = 'z,y,-x'
sg_symmetry(22, 209) = 'z,-y,x'
sg_symmetry(23, 209) = '-z,y,x'
sg_symmetry(24, 209) = '-z,-y,-x'
sg_symmetry(25, 209) = 'x,1/2+y,1/2+z'
sg_symmetry(26, 209) = '-x,1/2-y,1/2+z'
sg_symmetry(27, 209) = '-x,1/2+y,1/2-z'
sg_symmetry(28, 209) = 'x,1/2-y,1/2-z'
sg_symmetry(29, 209) = 'z,1/2+x,1/2+y'
sg_symmetry(30, 209) = 'z,1/2-x,1/2-y'
sg_symmetry(31, 209) = '-z,1/2-x,1/2+y'
sg_symmetry(32, 209) = '-z,1/2+x,1/2-y'
sg_symmetry(33, 209) = 'y,1/2+z,1/2+x'
sg_symmetry(34, 209) = '-y,1/2+z,1/2-x'
sg_symmetry(35, 209) = 'y,1/2-z,1/2-x'
sg_symmetry(36, 209) = '-y,1/2-z,1/2+x'
sg_symmetry(37, 209) = 'y,1/2+x,1/2-z'
sg_symmetry(38, 209) = '-y,1/2-x,1/2-z'
sg_symmetry(39, 209) = 'y,1/2-x,1/2+z'
sg_symmetry(40, 209) = '-y,1/2+x,1/2+z'
sg_symmetry(41, 209) = 'x,1/2+z,1/2-y'
sg_symmetry(42, 209) = '-x,1/2+z,1/2+y'
sg_symmetry(43, 209) = '-x,1/2-z,1/2-y'
sg_symmetry(44, 209) = 'x,1/2-z,1/2+y'
sg_symmetry(45, 209) = 'z,1/2+y,1/2-x'
sg_symmetry(46, 209) = 'z,1/2-y,1/2+x'
sg_symmetry(47, 209) = '-z,1/2+y,1/2+x'
sg_symmetry(48, 209) = '-z,1/2-y,1/2-x'
sg_symmetry(49, 209) = '1/2+x,y,1/2+z'
sg_symmetry(50, 209) = '1/2-x,-y,1/2+z'
sg_symmetry(51, 209) = '1/2-x,y,1/2-z'
sg_symmetry(52, 209) = '1/2+x,-y,1/2-z'
sg_symmetry(53, 209) = '1/2+z,x,1/2+y'
sg_symmetry(54, 209) = '1/2+z,-x,1/2-y'
sg_symmetry(55, 209) = '1/2-z,-x,1/2+y'
sg_symmetry(56, 209) = '1/2-z,x,1/2-y'
sg_symmetry(57, 209) = '1/2+y,z,1/2+x'
sg_symmetry(58, 209) = '1/2-y,z,1/2-x'
sg_symmetry(59, 209) = '1/2+y,-z,1/2-x'
sg_symmetry(60, 209) = '1/2-y,-z,1/2+x'
sg_symmetry(61, 209) = '1/2+y,x,1/2-z'
sg_symmetry(62, 209) = '1/2-y,-x,1/2-z'
sg_symmetry(63, 209) = '1/2+y,-x,1/2+z'
sg_symmetry(64, 209) = '1/2-y,x,1/2+z'
sg_symmetry(65, 209) = '1/2+x,z,1/2-y'
sg_symmetry(66, 209) = '1/2-x,z,1/2+y'
sg_symmetry(67, 209) = '1/2-x,-z,1/2-y'
sg_symmetry(68, 209) = '1/2+x,-z,1/2+y'
sg_symmetry(69, 209) = '1/2+z,y,1/2-x'
sg_symmetry(70, 209) = '1/2+z,-y,1/2+x'
sg_symmetry(71, 209) = '1/2-z,y,1/2+x'
sg_symmetry(72, 209) = '1/2-z,-y,1/2-x'
sg_symmetry(73, 209) = '1/2+x,1/2+y,z'
sg_symmetry(74, 209) = '1/2-x,1/2-y,z'
sg_symmetry(75, 209) = '1/2-x,1/2+y,-z'
sg_symmetry(76, 209) = '1/2+x,1/2-y,-z'
sg_symmetry(77, 209) = '1/2+z,1/2+x,y'
sg_symmetry(78, 209) = '1/2+z,1/2-x,-y'
sg_symmetry(79, 209) = '1/2-z,1/2-x,y'
sg_symmetry(80, 209) = '1/2-z,1/2+x,-y'
sg_symmetry(81, 209) = '1/2+y,1/2+z,x'
sg_symmetry(82, 209) = '1/2-y,1/2+z,-x'
sg_symmetry(83, 209) = '1/2+y,1/2-z,-x'
sg_symmetry(84, 209) = '1/2-y,1/2-z,x'
sg_symmetry(85, 209) = '1/2+y,1/2+x,-z'
sg_symmetry(86, 209) = '1/2-y,1/2-x,-z'
sg_symmetry(87, 209) = '1/2+y,1/2-x,z'
sg_symmetry(88, 209) = '1/2-y,1/2+x,z'
sg_symmetry(89, 209) = '1/2+x,1/2+z,-y'
sg_symmetry(90, 209) = '1/2-x,1/2+z,y'
sg_symmetry(91, 209) = '1/2-x,1/2-z,-y'
sg_symmetry(92, 209) = '1/2+x,1/2-z,y'
sg_symmetry(93, 209) = '1/2+z,1/2+y,-x'
sg_symmetry(94, 209) = '1/2+z,1/2-y,x'
sg_symmetry(95, 209) = '1/2-z,1/2+y,x'
sg_symmetry(96, 209) = '1/2-z,1/2-y,-x'
sg_name(210) = 'F4(1)32'
sg_patn(210) = 225
sg_symnum(210) = 96
sg_symmetry(1, 210) = 'x,y,z'
sg_symmetry(2, 210) = '-x,1/2-y,1/2+z'
sg_symmetry(3, 210) = '1/2-x,1/2+y,-z'
sg_symmetry(4, 210) = '1/2+x,-y,1/2-z'
sg_symmetry(5, 210) = 'z,x,y'
sg_symmetry(6, 210) = '1/2+z,-x,1/2-y'
sg_symmetry(7, 210) = '-z,1/2-x,1/2+y'
sg_symmetry(8, 210) = '1/2-z,1/2+x,-y'
sg_symmetry(9, 210) = 'y,z,x'
sg_symmetry(10, 210) = '1/2-y,1/2+z,-x'
sg_symmetry(11, 210) = '1/2+y,-z,1/2-x'
sg_symmetry(12, 210) = '-y,1/2-z,1/2+x'
sg_symmetry(13, 210) = '3/4+y,1/4+x,3/4-z'
sg_symmetry(14, 210) = '1/4-y,1/4-x,1/4-z'
sg_symmetry(15, 210) = '1/4+y,3/4-x,3/4+z'
sg_symmetry(16, 210) = '3/4-y,3/4+x,1/4+z'
sg_symmetry(17, 210) = '3/4+x,1/4+z,3/4-y'
sg_symmetry(18, 210) = '3/4-x,3/4+z,1/4+y'
sg_symmetry(19, 210) = '1/4-x,1/4-z,1/4-y'
sg_symmetry(20, 210) = '1/4+x,3/4-z,3/4+y'
sg_symmetry(21, 210) = '3/4+z,1/4+y,3/4-x'
sg_symmetry(22, 210) = '1/4+z,3/4-y,3/4+x'
sg_symmetry(23, 210) = '3/4-z,3/4+y,1/4+x'
sg_symmetry(24, 210) = '1/4-z,1/4-y,1/4-x'
sg_symmetry(25, 210) = 'x,1/2+y,1/2+z'
sg_symmetry(26, 210) = '-x,-y,z'
sg_symmetry(27, 210) = '1/2-x,y,1/2-z'
sg_symmetry(28, 210) = '1/2+x,1/2-y,-z'
sg_symmetry(29, 210) = 'z,1/2+x,1/2+y'
sg_symmetry(30, 210) = '1/2+z,1/2-x,-y'
sg_symmetry(31, 210) = '-z,-x,y'
sg_symmetry(32, 210) = '1/2-z,x,1/2-y'
sg_symmetry(33, 210) = 'y,1/2+z,1/2+x'
sg_symmetry(34, 210) = '1/2-y,z,1/2-x'
sg_symmetry(35, 210) = '1/2+y,1/2-z,-x'
sg_symmetry(36, 210) = '-y,-z,x'
sg_symmetry(37, 210) = '3/4+y,3/4+x,1/4-z'
sg_symmetry(38, 210) = '1/4-y,3/4-x,3/4-z'
sg_symmetry(39, 210) = '1/4+y,1/4-x,1/4+z'
sg_symmetry(40, 210) = '3/4-y,1/4+x,3/4+z'
sg_symmetry(41, 210) = '3/4+x,3/4+z,1/4-y'
sg_symmetry(42, 210) = '3/4-x,1/4+z,3/4+y'
sg_symmetry(43, 210) = '1/4-x,3/4-z,3/4-y'
sg_symmetry(44, 210) = '1/4+x,1/4-z,1/4+y'
sg_symmetry(45, 210) = '3/4+z,3/4+y,1/4-x'
sg_symmetry(46, 210) = '1/4+z,1/4-y,1/4+x'
sg_symmetry(47, 210) = '3/4-z,1/4+y,3/4+x'
sg_symmetry(48, 210) = '1/4-z,3/4-y,3/4-x'
sg_symmetry(49, 210) = '1/2+x,y,1/2+z'
sg_symmetry(50, 210) = '1/2-x,1/2-y,z'
sg_symmetry(51, 210) = '-x,1/2+y,1/2-z'
sg_symmetry(52, 210) = 'x,-y,-z'
sg_symmetry(53, 210) = '1/2+z,x,1/2+y'
sg_symmetry(54, 210) = 'z,-x,-y'
sg_symmetry(55, 210) = '1/2-z,1/2-x,y'
sg_symmetry(56, 210) = '-z,1/2+x,1/2-y'
sg_symmetry(57, 210) = '1/2+y,z,1/2+x'
sg_symmetry(58, 210) = '-y,1/2+z,1/2-x'
sg_symmetry(59, 210) = 'y,-z,-x'
sg_symmetry(60, 210) = '1/2-y,1/2-z,x'
sg_symmetry(61, 210) = '1/4+y,1/4+x,1/4-z'
sg_symmetry(62, 210) = '3/4-y,1/4-x,3/4-z'
sg_symmetry(63, 210) = '3/4+y,3/4-x,1/4+z'
sg_symmetry(64, 210) = '1/4-y,3/4+x,3/4+z'
sg_symmetry(65, 210) = '1/4+x,1/4+z,1/4-y'
sg_symmetry(66, 210) = '1/4-x,3/4+z,3/4+y'
sg_symmetry(67, 210) = '3/4-x,1/4-z,3/4-y'
sg_symmetry(68, 210) = '3/4+x,3/4-z,1/4+y'
sg_symmetry(69, 210) = '1/4+z,1/4+y,1/4-x'
sg_symmetry(70, 210) = '3/4+z,3/4-y,1/4+x'
sg_symmetry(71, 210) = '1/4-z,3/4+y,3/4+x'
sg_symmetry(72, 210) = '3/4-z,1/4-y,3/4-x'
sg_symmetry(73, 210) = '1/2+x,1/2+y,z'
sg_symmetry(74, 210) = '1/2-x,-y,1/2+z'
sg_symmetry(75, 210) = '-x,y,-z'
sg_symmetry(76, 210) = 'x,1/2-y,1/2-z'
sg_symmetry(77, 210) = '1/2+z,1/2+x,y'
sg_symmetry(78, 210) = 'z,1/2-x,1/2-y'
sg_symmetry(79, 210) = '1/2-z,-x,1/2+y'
sg_symmetry(80, 210) = '-z,x,-y'
sg_symmetry(81, 210) = '1/2+y,1/2+z,x'
sg_symmetry(82, 210) = '-y,z,-x'
sg_symmetry(83, 210) = 'y,1/2-z,1/2-x'
sg_symmetry(84, 210) = '1/2-y,-z,1/2+x'
sg_symmetry(85, 210) = '1/4+y,3/4+x,3/4-z'
sg_symmetry(86, 210) = '3/4-y,3/4-x,1/4-z'
sg_symmetry(87, 210) = '3/4+y,1/4-x,3/4+z'
sg_symmetry(88, 210) = '1/4-y,1/4+x,1/4+z'
sg_symmetry(89, 210) = '1/4+x,3/4+z,3/4-y'
sg_symmetry(90, 210) = '1/4-x,1/4+z,1/4+y'
sg_symmetry(91, 210) = '3/4-x,3/4-z,1/4-y'
sg_symmetry(92, 210) = '3/4+x,1/4-z,3/4+y'
sg_symmetry(93, 210) = '1/4+z,3/4+y,3/4-x'
sg_symmetry(94, 210) = '3/4+z,1/4-y,3/4+x'
sg_symmetry(95, 210) = '1/4-z,1/4+y,1/4+x'
sg_symmetry(96, 210) = '3/4-z,3/4-y,1/4-x'
sg_name(211) = 'I432'
sg_patn(211) = 229
sg_symnum(211) = 48
sg_symmetry(1, 211) = 'x,y,z'
sg_symmetry(2, 211) = '-x,-y,z'
sg_symmetry(3, 211) = '-x,y,-z'
sg_symmetry(4, 211) = 'x,-y,-z'
sg_symmetry(5, 211) = 'z,x,y'
sg_symmetry(6, 211) = 'z,-x,-y'
sg_symmetry(7, 211) = '-z,-x,y'
sg_symmetry(8, 211) = '-z,x,-y'
sg_symmetry(9, 211) = 'y,z,x'
sg_symmetry(10, 211) = '-y,z,-x'
sg_symmetry(11, 211) = 'y,-z,-x'
sg_symmetry(12, 211) = '-y,-z,x'
sg_symmetry(13, 211) = 'y,x,-z'
sg_symmetry(14, 211) = '-y,-x,-z'
sg_symmetry(15, 211) = 'y,-x,z'
sg_symmetry(16, 211) = '-y,x,z'
sg_symmetry(17, 211) = 'x,z,-y'
sg_symmetry(18, 211) = '-x,z,y'
sg_symmetry(19, 211) = '-x,-z,-y'
sg_symmetry(20, 211) = 'x,-z,y'
sg_symmetry(21, 211) = 'z,y,-x'
sg_symmetry(22, 211) = 'z,-y,x'
sg_symmetry(23, 211) = '-z,y,x'
sg_symmetry(24, 211) = '-z,-y,-x'
sg_symmetry(25, 211) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(26, 211) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(27, 211) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(28, 211) = '1/2+x,1/2-y,1/2-z'
sg_symmetry(29, 211) = '1/2+z,1/2+x,1/2+y'
sg_symmetry(30, 211) = '1/2+z,1/2-x,1/2-y'
sg_symmetry(31, 211) = '1/2-z,1/2-x,1/2+y'
sg_symmetry(32, 211) = '1/2-z,1/2+x,1/2-y'
sg_symmetry(33, 211) = '1/2+y,1/2+z,1/2+x'
sg_symmetry(34, 211) = '1/2-y,1/2+z,1/2-x'
sg_symmetry(35, 211) = '1/2+y,1/2-z,1/2-x'
sg_symmetry(36, 211) = '1/2-y,1/2-z,1/2+x'
sg_symmetry(37, 211) = '1/2+y,1/2+x,1/2-z'
sg_symmetry(38, 211) = '1/2-y,1/2-x,1/2-z'
sg_symmetry(39, 211) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(40, 211) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(41, 211) = '1/2+x,1/2+z,1/2-y'
sg_symmetry(42, 211) = '1/2-x,1/2+z,1/2+y'
sg_symmetry(43, 211) = '1/2-x,1/2-z,1/2-y'
sg_symmetry(44, 211) = '1/2+x,1/2-z,1/2+y'
sg_symmetry(45, 211) = '1/2+z,1/2+y,1/2-x'
sg_symmetry(46, 211) = '1/2+z,1/2-y,1/2+x'
sg_symmetry(47, 211) = '1/2-z,1/2+y,1/2+x'
sg_symmetry(48, 211) = '1/2-z,1/2-y,1/2-x'
sg_name(212) = 'P4(3)32'
sg_patn(212) = 221
sg_symnum(212) = 24
sg_symmetry(1, 212) = 'x,y,z'
sg_symmetry(2, 212) = '1/2-x,-y,1/2+z'
sg_symmetry(3, 212) = '-x,1/2+y,1/2-z'
sg_symmetry(4, 212) = '1/2+x,1/2-y,-z'
sg_symmetry(5, 212) = 'z,x,y'
sg_symmetry(6, 212) = '1/2+z,1/2-x,-y'
sg_symmetry(7, 212) = '1/2-z,-x,1/2+y'
sg_symmetry(8, 212) = '-z,1/2+x,1/2-y'
sg_symmetry(9, 212) = 'y,z,x'
sg_symmetry(10, 212) = '-y,1/2+z,1/2-x'
sg_symmetry(11, 212) = '1/2+y,1/2-z,-x'
sg_symmetry(12, 212) = '1/2-y,-z,1/2+x'
sg_symmetry(13, 212) = '1/4+y,3/4+x,3/4-z'
sg_symmetry(14, 212) = '1/4-y,1/4-x,1/4-z'
sg_symmetry(15, 212) = '3/4+y,3/4-x,1/4+z'
sg_symmetry(16, 212) = '3/4-y,1/4+x,3/4+z'
sg_symmetry(17, 212) = '1/4+x,3/4+z,3/4-y'
sg_symmetry(18, 212) = '3/4-x,1/4+z,3/4+y'
sg_symmetry(19, 212) = '1/4-x,1/4-z,1/4-y'
sg_symmetry(20, 212) = '3/4+x,3/4-z,1/4+y'
sg_symmetry(21, 212) = '1/4+z,3/4+y,3/4-x'
sg_symmetry(22, 212) = '3/4+z,3/4-y,1/4+x'
sg_symmetry(23, 212) = '3/4-z,1/4+y,3/4+x'
sg_symmetry(24, 212) = '1/4-z,1/4-y,1/4-x'
sg_name(213) = 'P4(1)32'
sg_patn(213) = 221
sg_symnum(213) = 24
sg_symmetry(1, 213) = 'x,y,z'
sg_symmetry(2, 213) = '1/2-x,-y,1/2+z'
sg_symmetry(3, 213) = '-x,1/2+y,1/2-z'
sg_symmetry(4, 213) = '1/2+x,1/2-y,-z'
sg_symmetry(5, 213) = 'z,x,y'
sg_symmetry(6, 213) = '1/2+z,1/2-x,-y'
sg_symmetry(7, 213) = '1/2-z,-x,1/2+y'
sg_symmetry(8, 213) = '-z,1/2+x,1/2-y'
sg_symmetry(9, 213) = 'y,z,x'
sg_symmetry(10, 213) = '-y,1/2+z,1/2-x'
sg_symmetry(11, 213) = '1/2+y,1/2-z,-x'
sg_symmetry(12, 213) = '1/2-y,-z,1/2+x'
sg_symmetry(13, 213) = '3/4+y,1/4+x,1/4-z'
sg_symmetry(14, 213) = '3/4-y,3/4-x,3/4-z'
sg_symmetry(15, 213) = '1/4+y,1/4-x,3/4+z'
sg_symmetry(16, 213) = '1/4-y,3/4+x,1/4+z'
sg_symmetry(17, 213) = '3/4+x,1/4+z,1/4-y'
sg_symmetry(18, 213) = '1/4-x,3/4+z,1/4+y'
sg_symmetry(19, 213) = '3/4-x,3/4-z,3/4-y'
sg_symmetry(20, 213) = '1/4+x,1/4-z,3/4+y'
sg_symmetry(21, 213) = '3/4+z,1/4+y,1/4-x'
sg_symmetry(22, 213) = '1/4+z,1/4-y,3/4+x'
sg_symmetry(23, 213) = '1/4-z,3/4+y,1/4+x'
sg_symmetry(24, 213) = '3/4-z,3/4-y,3/4-x'
sg_name(214) = 'I4(1)32'
sg_patn(214) = 229
sg_symnum(214) = 48
sg_symmetry(1, 214) = 'x,y,z'
sg_symmetry(2, 214) = '1/2-x,-y,1/2+z'
sg_symmetry(3, 214) = '-x,1/2+y,1/2-z'
sg_symmetry(4, 214) = '1/2+x,1/2-y,-z'
sg_symmetry(5, 214) = 'z,x,y'
sg_symmetry(6, 214) = '1/2+z,1/2-x,-y'
sg_symmetry(7, 214) = '1/2-z,-x,1/2+y'
sg_symmetry(8, 214) = '-z,1/2+x,1/2-y'
sg_symmetry(9, 214) = 'y,z,x'
sg_symmetry(10, 214) = '-y,1/2+z,1/2-x'
sg_symmetry(11, 214) = '1/2+y,1/2-z,-x'
sg_symmetry(12, 214) = '1/2-y,-z,1/2+x'
sg_symmetry(13, 214) = '3/4+y,1/4+x,1/4-z'
sg_symmetry(14, 214) = '3/4-y,3/4-x,3/4-z'
sg_symmetry(15, 214) = '1/4+y,1/4-x,3/4+z'
sg_symmetry(16, 214) = '1/4-y,3/4+x,1/4+z'
sg_symmetry(17, 214) = '3/4+x,1/4+z,1/4-y'
sg_symmetry(18, 214) = '1/4-x,3/4+z,1/4+y'
sg_symmetry(19, 214) = '3/4-x,3/4-z,3/4-y'
sg_symmetry(20, 214) = '1/4+x,1/4-z,3/4+y'
sg_symmetry(21, 214) = '3/4+z,1/4+y,1/4-x'
sg_symmetry(22, 214) = '1/4+z,1/4-y,3/4+x'
sg_symmetry(23, 214) = '1/4-z,3/4+y,1/4+x'
sg_symmetry(24, 214) = '3/4-z,3/4-y,3/4-x'
sg_symmetry(25, 214) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(26, 214) = '-x,1/2-y,z'
sg_symmetry(27, 214) = '1/2-x,y,-z'
sg_symmetry(28, 214) = 'x,-y,1/2-z'
sg_symmetry(29, 214) = '1/2+z,1/2+x,1/2+y'
sg_symmetry(30, 214) = 'z,-x,1/2-y'
sg_symmetry(31, 214) = '-z,1/2-x,y'
sg_symmetry(32, 214) = '1/2-z,x,-y'
sg_symmetry(33, 214) = '1/2+y,1/2+z,1/2+x'
sg_symmetry(34, 214) = '1/2-y,z,-x'
sg_symmetry(35, 214) = 'y,-z,1/2-x'
sg_symmetry(36, 214) = '-y,1/2-z,x'
sg_symmetry(37, 214) = '1/4+y,3/4+x,3/4-z'
sg_symmetry(38, 214) = '1/4-y,1/4-x,1/4-z'
sg_symmetry(39, 214) = '3/4+y,3/4-x,1/4+z'
sg_symmetry(40, 214) = '3/4-y,1/4+x,3/4+z'
sg_symmetry(41, 214) = '1/4+x,3/4+z,3/4-y'
sg_symmetry(42, 214) = '3/4-x,1/4+z,3/4+y'
sg_symmetry(43, 214) = '1/4-x,1/4-z,1/4-y'
sg_symmetry(44, 214) = '3/4+x,3/4-z,1/4+y'
sg_symmetry(45, 214) = '1/4+z,3/4+y,3/4-x'
sg_symmetry(46, 214) = '3/4+z,3/4-y,1/4+x'
sg_symmetry(47, 214) = '3/4-z,1/4+y,3/4+x'
sg_symmetry(48, 214) = '1/4-z,1/4-y,1/4-x'
sg_name(215) = 'P-43m'
sg_patn(215) = 221
sg_symnum(215) = 24
sg_symmetry(1, 215) = 'x,y,z'
sg_symmetry(2, 215) = '-x,-y,z'
sg_symmetry(3, 215) = '-x,y,-z'
sg_symmetry(4, 215) = 'x,-y,-z'
sg_symmetry(5, 215) = 'z,x,y'
sg_symmetry(6, 215) = 'z,-x,-y'
sg_symmetry(7, 215) = '-z,-x,y'
sg_symmetry(8, 215) = '-z,x,-y'
sg_symmetry(9, 215) = 'y,z,x'
sg_symmetry(10, 215) = '-y,z,-x'
sg_symmetry(11, 215) = 'y,-z,-x'
sg_symmetry(12, 215) = '-y,-z,x'
sg_symmetry(13, 215) = 'y,x,z'
sg_symmetry(14, 215) = '-y,-x,z'
sg_symmetry(15, 215) = 'y,-x,-z'
sg_symmetry(16, 215) = '-y,x,-z'
sg_symmetry(17, 215) = 'x,z,y'
sg_symmetry(18, 215) = '-x,z,-y'
sg_symmetry(19, 215) = '-x,-z,y'
sg_symmetry(20, 215) = 'x,-z,-y'
sg_symmetry(21, 215) = 'z,y,x'
sg_symmetry(22, 215) = 'z,-y,-x'
sg_symmetry(23, 215) = '-z,y,-x'
sg_symmetry(24, 215) = '-z,-y,x'
sg_name(216) = 'F4-3m'
sg_patn(216) = 225
sg_symnum(216) = 96
sg_symmetry(1, 216) = 'x,y,z'
sg_symmetry(2, 216) = '-x,-y,z'
sg_symmetry(3, 216) = '-x,y,-z'
sg_symmetry(4, 216) = 'x,-y,-z'
sg_symmetry(5, 216) = 'z,x,y'
sg_symmetry(6, 216) = 'z,-x,-y'
sg_symmetry(7, 216) = '-z,-x,y'
sg_symmetry(8, 216) = '-z,x,-y'
sg_symmetry(9, 216) = 'y,z,x'
sg_symmetry(10, 216) = '-y,z,-x'
sg_symmetry(11, 216) = 'y,-z,-x'
sg_symmetry(12, 216) = '-y,-z,x'
sg_symmetry(13, 216) = 'y,x,z'
sg_symmetry(14, 216) = '-y,-x,z'
sg_symmetry(15, 216) = 'y,-x,-z'
sg_symmetry(16, 216) = '-y,x,-z'
sg_symmetry(17, 216) = 'x,z,y'
sg_symmetry(18, 216) = '-x,z,-y'
sg_symmetry(19, 216) = '-x,-z,y'
sg_symmetry(20, 216) = 'x,-z,-y'
sg_symmetry(21, 216) = 'z,y,x'
sg_symmetry(22, 216) = 'z,-y,-x'
sg_symmetry(23, 216) = '-z,y,-x'
sg_symmetry(24, 216) = '-z,-y,x'
sg_symmetry(25, 216) = 'x,1/2+y,1/2+z'
sg_symmetry(26, 216) = '-x,1/2-y,1/2+z'
sg_symmetry(27, 216) = '-x,1/2+y,1/2-z'
sg_symmetry(28, 216) = 'x,1/2-y,1/2-z'
sg_symmetry(29, 216) = 'z,1/2+x,1/2+y'
sg_symmetry(30, 216) = 'z,1/2-x,1/2-y'
sg_symmetry(31, 216) = '-z,1/2-x,1/2+y'
sg_symmetry(32, 216) = '-z,1/2+x,1/2-y'
sg_symmetry(33, 216) = 'y,1/2+z,1/2+x'
sg_symmetry(34, 216) = '-y,1/2+z,1/2-x'
sg_symmetry(35, 216) = 'y,1/2-z,1/2-x'
sg_symmetry(36, 216) = '-y,1/2-z,1/2+x'
sg_symmetry(37, 216) = 'y,1/2+x,1/2+z'
sg_symmetry(38, 216) = '-y,1/2-x,1/2+z'
sg_symmetry(39, 216) = 'y,1/2-x,1/2-z'
sg_symmetry(40, 216) = '-y,1/2+x,1/2-z'
sg_symmetry(41, 216) = 'x,1/2+z,1/2+y'
sg_symmetry(42, 216) = '-x,1/2+z,1/2-y'
sg_symmetry(43, 216) = '-x,1/2-z,1/2+y'
sg_symmetry(44, 216) = 'x,1/2-z,1/2-y'
sg_symmetry(45, 216) = 'z,1/2+y,1/2+x'
sg_symmetry(46, 216) = 'z,1/2-y,1/2-x'
sg_symmetry(47, 216) = '-z,1/2+y,1/2-x'
sg_symmetry(48, 216) = '-z,1/2-y,1/2+x'
sg_symmetry(49, 216) = '1/2+x,y,1/2+z'
sg_symmetry(50, 216) = '1/2-x,-y,1/2+z'
sg_symmetry(51, 216) = '1/2-x,y,1/2-z'
sg_symmetry(52, 216) = '1/2+x,-y,1/2-z'
sg_symmetry(53, 216) = '1/2+z,x,1/2+y'
sg_symmetry(54, 216) = '1/2+z,-x,1/2-y'
sg_symmetry(55, 216) = '1/2-z,-x,1/2+y'
sg_symmetry(56, 216) = '1/2-z,x,1/2-y'
sg_symmetry(57, 216) = '1/2+y,z,1/2+x'
sg_symmetry(58, 216) = '1/2-y,z,1/2-x'
sg_symmetry(59, 216) = '1/2+y,-z,1/2-x'
sg_symmetry(60, 216) = '1/2-y,-z,1/2+x'
sg_symmetry(61, 216) = '1/2+y,x,1/2+z'
sg_symmetry(62, 216) = '1/2-y,-x,1/2+z'
sg_symmetry(63, 216) = '1/2+y,-x,1/2-z'
sg_symmetry(64, 216) = '1/2-y,x,1/2-z'
sg_symmetry(65, 216) = '1/2+x,z,1/2+y'
sg_symmetry(66, 216) = '1/2-x,z,1/2-y'
sg_symmetry(67, 216) = '1/2-x,-z,1/2+y'
sg_symmetry(68, 216) = '1/2+x,-z,1/2-y'
sg_symmetry(69, 216) = '1/2+z,y,1/2+x'
sg_symmetry(70, 216) = '1/2+z,-y,1/2-x'
sg_symmetry(71, 216) = '1/2-z,y,1/2-x'
sg_symmetry(72, 216) = '1/2-z,-y,1/2+x'
sg_symmetry(73, 216) = '1/2+x,1/2+y,z'
sg_symmetry(74, 216) = '1/2-x,1/2-y,z'
sg_symmetry(75, 216) = '1/2-x,1/2+y,-z'
sg_symmetry(76, 216) = '1/2+x,1/2-y,-z'
sg_symmetry(77, 216) = '1/2+z,1/2+x,y'
sg_symmetry(78, 216) = '1/2+z,1/2-x,-y'
sg_symmetry(79, 216) = '1/2-z,1/2-x,y'
sg_symmetry(80, 216) = '1/2-z,1/2+x,-y'
sg_symmetry(81, 216) = '1/2+y,1/2+z,x'
sg_symmetry(82, 216) = '1/2-y,1/2+z,-x'
sg_symmetry(83, 216) = '1/2+y,1/2-z,-x'
sg_symmetry(84, 216) = '1/2-y,1/2-z,x'
sg_symmetry(85, 216) = '1/2+y,1/2+x,z'
sg_symmetry(86, 216) = '1/2-y,1/2-x,z'
sg_symmetry(87, 216) = '1/2+y,1/2-x,-z'
sg_symmetry(88, 216) = '1/2-y,1/2+x,-z'
sg_symmetry(89, 216) = '1/2+x,1/2+z,y'
sg_symmetry(90, 216) = '1/2-x,1/2+z,-y'
sg_symmetry(91, 216) = '1/2-x,1/2-z,y'
sg_symmetry(92, 216) = '1/2+x,1/2-z,-y'
sg_symmetry(93, 216) = '1/2+z,1/2+y,x'
sg_symmetry(94, 216) = '1/2+z,1/2-y,-x'
sg_symmetry(95, 216) = '1/2-z,1/2+y,-x'
sg_symmetry(96, 216) = '1/2-z,1/2-y,x'
sg_name(217) = 'I-43m'
sg_patn(217) = 229
sg_symnum(217) = 48
sg_symmetry(1, 217) = 'x,y,z'
sg_symmetry(2, 217) = '-x,-y,z'
sg_symmetry(3, 217) = '-x,y,-z'
sg_symmetry(4, 217) = 'x,-y,-z'
sg_symmetry(5, 217) = 'z,x,y'
sg_symmetry(6, 217) = 'z,-x,-y'
sg_symmetry(7, 217) = '-z,-x,y'
sg_symmetry(8, 217) = '-z,x,-y'
sg_symmetry(9, 217) = 'y,z,x'
sg_symmetry(10, 217) = '-y,z,-x'
sg_symmetry(11, 217) = 'y,-z,-x'
sg_symmetry(12, 217) = '-y,-z,x'
sg_symmetry(13, 217) = 'y,x,z'
sg_symmetry(14, 217) = '-y,-x,z'
sg_symmetry(15, 217) = 'y,-x,-z'
sg_symmetry(16, 217) = '-y,x,-z'
sg_symmetry(17, 217) = 'x,z,y'
sg_symmetry(18, 217) = '-x,z,-y'
sg_symmetry(19, 217) = '-x,-z,y'
sg_symmetry(20, 217) = 'x,-z,-y'
sg_symmetry(21, 217) = 'z,y,x'
sg_symmetry(22, 217) = 'z,-y,-x'
sg_symmetry(23, 217) = '-z,y,-x'
sg_symmetry(24, 217) = '-z,-y,x'
sg_symmetry(25, 217) = '1/2+x,1/2+y,1/2+z'
sg_symmetry(26, 217) = '1/2-x,1/2-y,1/2+z'
sg_symmetry(27, 217) = '1/2-x,1/2+y,1/2-z'
sg_symmetry(28, 217) = '1/2+x,1/2-y,1/2-z'
sg_symmetry(29, 217) = '1/2+z,1/2+x,1/2+y'
sg_symmetry(30, 217) = '1/2+z,1/2-x,1/2-y'
sg_symmetry(31, 217) = '1/2-z,1/2-x,1/2+y'
sg_symmetry(32, 217) = '1/2-z,1/2+x,1/2-y'
sg_symmetry(33, 217) = '1/2+y,1/2+z,1/2+x'
sg_symmetry(34, 217) = '1/2-y,1/2+z,1/2-x'
sg_symmetry(35, 217) = '1/2+y,1/2-z,1/2-x'
sg_symmetry(36, 217) = '1/2-y,1/2-z,1/2+x'
sg_symmetry(37, 217) = '1/2+y,1/2+x,1/2+z'
sg_symmetry(38, 217) = '1/2-y,1/2-x,1/2+z'
sg_symmetry(39, 217) = '1/2+y,1/2-x,1/2-z'
sg_symmetry(40, 217) = '1/2-y,1/2+x,1/2-z'
sg_symmetry(41, 217) = '1/2+x,1/2+z,1/2+y'
sg_symmetry(42, 217) = '1/2-x,1/2+z,1/2-y'
sg_symmetry(43, 217) = '1/2-x,1/2-z,1/2+y'
sg_symmetry(44, 217) = '1/2+x,1/2-z,1/2-y'
sg_symmetry(45, 217) = '1/2+z,1/2+y,1/2+x'
sg_symmetry(46, 217) = '1/2+z,1/2-y,1/2-x'
sg_symmetry(47, 217) = '1/2-z,1/2+y,1/2-x'
sg_symmetry(48, 217) = '1/2-z,1/2-y,1/2+x'
sg_name(218) = 'P-43n'
sg_patn(218) = 221
sg_symnum(218) = 24
sg_symmetry(1, 218) = 'x,y,z'
sg_symmetry(2, 218) = '-x,-y,z'
sg_symmetry(3, 218) = '-x,y,-z'
sg_symmetry(4, 218) = 'x,-y,-z'
sg_symmetry(5, 218) = 'z,x,y'
sg_symmetry(6, 218) = 'z,-x,-y'
sg_symmetry(7, 218) = '-z,-x,y'
sg_symmetry(8, 218) = '-z,x,-y'
sg_symmetry(9, 218) = 'y,z,x'
sg_symmetry(10, 218) = '-y,z,-x'
sg_symmetry(11, 218) = 'y,-z,-x'
sg_symmetry(12, 218) = '-y,-z,x'
sg_symmetry(13, 218) = '1/2+y,1/2+x,1/2+z'
sg_symmetry(14, 218) = '1/2-y,1/2-x,1/2+z'
sg_symmetry(15, 218) = '1/2+y,1/2-x,1/2-z'
sg_symmetry(16, 218) = '1/2-y,1/2+x,1/2-z'
sg_symmetry(17, 218) = '1/2+x,1/2+z,1/2+y'
sg_symmetry(18, 218) = '1/2-x,1/2+z,1/2-y'
sg_symmetry(19, 218) = '1/2-x,1/2-z,1/2+y'
sg_symmetry(20, 218) = '1/2+x,1/2-z,1/2-y'
sg_symmetry(21, 218) = '1/2+z,1/2+y,1/2+x'
sg_symmetry(22, 218) = '1/2+z,1/2-y,1/2-x'
sg_symmetry(23, 218) = '1/2-z,1/2+y,1/2-x'
sg_symmetry(24, 218) = '1/2-z,1/2-y,1/2+x'
sg_name(219) = 'F-43c'
sg_patn(219) = 225
sg_symnum(219) = 96
sg_symmetry(1, 219) = 'x,y,z'
sg_symmetry(2, 219) = '-x,-y,z'
sg_symmetry(3, 219) = '-x,y,-z'
sg_symmetry(4, 219) = 'x,-y,-z'
sg_symmetry(5, 219) = 'z,x,y'
sg_symmetry(6, 219) = 'z,-x,-y'
sg_symmetry(7, 219) = '-z,-x,y'
sg_symmetry(8, 219) = '-z,x,-y'
sg_symmetry(9, 219) = 'y,z,x'
sg_symmetry(10, 219) = '-y,z,-x'
sg_symmetry(11, 219) = 'y,-z,-x'
sg_symmetry(12, 219) = '-y,-z,x'
sg_symmetry(13, 219) = '1/2+y,1/2+x,1/2+z'
sg_symmetry(14, 219) = '1/2-y,1/2-x,1/2+z'
sg_symmetry(15, 219) = '1/2+y,1/2-x,1/2-z'
sg_symmetry(16, 219) = '1/2-y,1/2+x,1/2-z'
sg_symmetry(17, 219) = '1/2+x,1/2+z,1/2+y'
sg_symmetry(18, 219) = '1/2-x,1/2+z,1/2-y'
sg_symmetry(19, 219) = '1/2-x,1/2-z,1/2+y'
sg_symmetry(20, 219) = '1/2+x,1/2-z,1/2-y'
sg_symmetry(21, 219) = '1/2+z,1/2+y,1/2+x'
sg_symmetry(22, 219) = '1/2+z,1/2-y,1/2-x'
sg_symmetry(23, 219) = '1/2-z,1/2+y,1/2-x'
sg_symmetry(24, 219) = '1/2-z,1/2-y,1/2+x'
sg_symmetry(25, 219) = 'x,y+1/2,z+1/2'
sg_symmetry(26, 219) = '-x,-y+1/2,z+1/2'
sg_symmetry(27, 219) = '-x,y+1/2,-z+1/2'
sg_symmetry(28, 219) = 'x,-y+1/2,-z+1/2'
sg_symmetry(29, 219) = 'z,x+1/2,y+1/2'
sg_symmetry(30, 219) = 'z,-x+1/2,-y+1/2'
sg_symmetry(31, 219) = '-z,-x+1/2,y+1/2'
sg_symmetry(32, 219) = '-z,x+1/2,-y+1/2'
sg_symmetry(33, 219) = 'y,z+1/2,x+1/2'
sg_symmetry(34, 219) = '-y,z+1/2,-x+1/2'
sg_symmetry(35, 219) = 'y,-z+1/2,-x+1/2)'
sg_symmetry(36, 219) = '-y,-z+1/2,x+1/2'
sg_symmetry(37, 219) = '1/2+y,x,z'
sg_symmetry(38, 219) = '1/2-y,-x,z'
sg_symmetry(39, 219) = '1/2+y,-x,-z'
sg_symmetry(40, 219) = '1/2-y,+x,-z'
sg_symmetry(41, 219) = '1/2+x,+z,y'
sg_symmetry(42, 219) = '1/2-x,+z,-y'
sg_symmetry(43, 219) = '1/2-x,-z,y'
sg_symmetry(44, 219) = '1/2+x,-z,-y'
sg_symmetry(45, 219) = '1/2+z,+y,x'
sg_symmetry(46, 219) = '1/2+z,-y,-x'
sg_symmetry(47, 219) = '1/2-z,+y,-x'
sg_symmetry(48, 219) = '1/2-z,-y,x'
sg_symmetry(49, 219) = 'x+1/2,y,z+1/2'
sg_symmetry(50, 219) = '-x+1/2,-y,z+1/2'
sg_symmetry(51, 219) = '-x+1/2,y,-z+1/2'
sg_symmetry(52, 219) = 'x+1/2,-y,-z+1/2'
sg_symmetry(53, 219) = 'z+1/2,x,y+1/2'
sg_symmetry(54, 219) = 'z+1/2,-x,-y+1/2'
sg_symmetry(55, 219) = '-z+1/2,-x,y+1/2'
sg_symmetry(56, 219) = '-z+1/2,x,-y+1/2'
sg_symmetry(57, 219) = 'y+1/2,z,x+1/2'
sg_symmetry(58, 219) = '-y+1/2,z,-x+1/2'
sg_symmetry(59, 219) = 'y+1/2,-z,-x+1/2'
sg_symmetry(60, 219) = '-y+1/2,-z,x+1/2'
sg_symmetry(61, 219) = 'y,1/2+x,z'
sg_symmetry(62, 219) = '-y,1/2-x,z'
sg_symmetry(63, 219) = 'y,1/2-x,-z'
sg_symmetry(64, 219) = '-y,1/2+x,-z'
sg_symmetry(65, 219) = 'x,1/2+z,y'
sg_symmetry(66, 219) = '(-x,1/2+z,-y'
sg_symmetry(67, 219) = '-x,1/2-z,y'
sg_symmetry(68, 219) = 'x,1/2-z,-y'
sg_symmetry(69, 219) = 'z,1/2+y,x'
sg_symmetry(70, 219) = 'z,1/2-y,-x'
sg_symmetry(71, 219) = '-z,1/2+y,-x'
sg_symmetry(72, 219) = '-z,1/2-y,x'
sg_symmetry(73, 219) = 'x+1/2,y+1/2,z'
sg_symmetry(74, 219) = '-x+1/2,-y+1/2,z'
sg_symmetry(75, 219) = '-x+1/2,y+1/2,-z'
sg_symmetry(76, 219) = 'x+1/2,-y+1/2,-z'
sg_symmetry(77, 219) = 'z+1/2,x+1/2,y'
sg_symmetry(78, 219) = 'z+1/2,-x+1/2,-y'
sg_symmetry(79, 219) = '-z+1/2,-x+1/2,y'
sg_symmetry(80, 219) = '-z+1/2,x+1/2,-y'
sg_symmetry(81, 219) = 'y+1/2,z+1/2,x'
sg_symmetry(82, 219) = '-y+1/2,z+1/2,-x'
sg_symmetry(83, 219) = 'y+1/2,-z+1/2,-x'
sg_symmetry(84, 219) = '-y+1/2,-z+1/2,x'
sg_symmetry(85, 219) = 'y,x,1/2+z'
sg_symmetry(86, 219) = '-y,-x,1/2+z'
sg_symmetry(87, 219) = 'y,-x,1/2-z'
sg_symmetry(88, 219) = '-y,x,1/2-z'
sg_symmetry(89, 219) = 'x,z,1/2+y'
sg_symmetry(90, 219) = '-x,z,1/2-y'
sg_symmetry(91, 219) = '-x,-z,1/2+y'
sg_symmetry(92, 219) = 'x,-z,1/2-y'
sg_symmetry(93, 219) = 'z,y,1/2+x'
sg_symmetry(94, 219) = 'z,-y,1/2-x'
sg_symmetry(95, 219) = '-z,y,1/2-x'
sg_symmetry(96, 219) = '-z,-y,1/2+x'
sg_name(220) = 'I-43d'
sg_patn(220) = 229
sg_symnum(220) = 48
sg_symmetry(1, 220) = 'x,y,z'
sg_symmetry(2, 220) = '1/2-x,-y,1/2+z'
sg_symmetry(3, 220) = '-x,1/2+y,1/2-z'
sg_symmetry(4, 220) = '1/2+x,1/2-y,-z'
sg_symmetry(5, 220) = 'z,x,y'
sg_symmetry(6, 220) = '1/2+z,1/2-x,-y'
sg_symmetry(7, 220) = '1/2-z,-x,1/2+y'
sg_symmetry(8, 220) = '-z,1/2+x,1/2-y'
sg_symmetry(9, 220) = 'y,z,x'
sg_symmetry(10, 220) = '-y,1/2+z,1/2-x'
sg_symmetry(11, 220) = '1/2+y,1/2-z,-x'
sg_symmetry(12, 220) = '1/2-y,-z,1/2+x'
sg_symmetry(13, 220) = '1/4+y,1/4+x,1/4+z'
sg_symmetry(14, 220) = '1/4-y,3/4-x,3/4+z'
sg_symmetry(15, 220) = '3/4+y,1/4-x,3/4-z'
sg_symmetry(16, 220) = '3/4-y,3/4+x,1/4-z'
sg_symmetry(17, 220) = '1/4+x,1/4+z,1/4+y'
sg_symmetry(18, 220) = '3/4-x,3/4+z,1/4-y'
sg_symmetry(19, 220) = '1/4-x,3/4-z,3/4+y'
sg_symmetry(20, 220) = '3/4+x,1/4-z,3/4-y'
sg_symmetry(21, 220) = '1/4+z,1/4+y,1/4+x'
sg_symmetry(22, 220) = '3/4+z,1/4-y,3/4-x'
sg_symmetry(23, 220) = '3/4-z,3/4+y,1/4-x'
sg_symmetry(24, 220) = '1/4-z,3/4-y,3/4+x'
sg_symmetry(25, 220) = 'x+1/2,y+1/2,z+1/2'
sg_symmetry(26, 220) = '-x,-y+1/2,z'
sg_symmetry(27, 220) = '-x+1/2,y,-z'
sg_symmetry(28, 220) = 'x,-y,-z+1/2'
sg_symmetry(29, 220) = 'z+1/2,x+1/2,y+1/2'
sg_symmetry(30, 220) = 'z,-x,-y+1/2'
sg_symmetry(31, 220) = '-z,-x+1/2,y'
sg_symmetry(32, 220) = '-z+1/2,x,-y'
sg_symmetry(33, 220) = 'y+1/2,z+1/2,x+1/2'
sg_symmetry(34, 220) = '-y+1/2,z,-x'
sg_symmetry(35, 220) = 'y,-z,-x+1/2'
sg_symmetry(36, 220) = '-y,-z+1/2,x'
sg_symmetry(37, 220) = '3/4+y,3/4+x,3/4+z'
sg_symmetry(38, 220) = '3/4-y,1/4-x,1/4+z'
sg_symmetry(39, 220) = '1/4+y,3/4-x,1/4-z'
sg_symmetry(40, 220) = '1/4-y,1/4+x,3/4-z'
sg_symmetry(41, 220) = '3/4+x,3/4+z,3/4+y'
sg_symmetry(42, 220) = '1/4-x,1/4+z,3/4-y'
sg_symmetry(43, 220) = '3/4-x,1/4-z,1/4+y'
sg_symmetry(44, 220) = '1/4+x,3/4-z,1/4-y'
sg_symmetry(45, 220) = '3/4+z,3/4+y,3/4+x'
sg_symmetry(46, 220) = '1/4+z,3/4-y,1/4-x'
sg_symmetry(47, 220) = '1/4-z,1/4+y,3/4-x'
sg_symmetry(48, 220) = '3/4-z,1/4-y,1/4+x'
sg_name(221) = 'Pm-3m'
sg_patn(221) = 221
sg_symnum(221) = 48
sg_symmetry(1, 221) = 'x,y,z'
sg_symmetry(2, 221) = '-x,-y,z'
sg_symmetry(3, 221) = '-x,y,-z'
sg_symmetry(4, 221) = 'x,-y,-z'
sg_symmetry(5, 221) = 'z,x,y'
sg_symmetry(6, 221) = 'z,-x,-y'
sg_symmetry(7, 221) = '-z,-x,y'
sg_symmetry(8, 221) = '-z,x,-y'
sg_symmetry(9, 221) = 'y,z,x'
sg_symmetry(10, 221) = '-y,z,-x'
sg_symmetry(11, 221) = 'y,-z,-x'
sg_symmetry(12, 221) = '-y,-z,x'
sg_symmetry(13, 221) = 'y,x,-z'
sg_symmetry(14, 221) = '-y,-x,-z'
sg_symmetry(15, 221) = 'y,-x,z'
sg_symmetry(16, 221) = '-y,x,z'
sg_symmetry(17, 221) = 'x,z,-y'
sg_symmetry(18, 221) = '-x,z,y'
sg_symmetry(19, 221) = '-x,-z,-y'
sg_symmetry(20, 221) = 'x,-z,y'
sg_symmetry(21, 221) = 'z,y,-x'
sg_symmetry(22, 221) = 'z,-y,x'
sg_symmetry(23, 221) = '-z,y,x'
sg_symmetry(24, 221) = '-z,-y,-x'
sg_symmetry(25, 221) = '-x,-y,-z'
sg_symmetry(26, 221) = 'x,y,-z'
sg_symmetry(27, 221) = 'x,-y,z'
sg_symmetry(28, 221) = '-x,y,z'
sg_symmetry(29, 221) = '-z,-x,-y'
sg_symmetry(30, 221) = '-z,x,y'
sg_symmetry(31, 221) = 'z,x,-y'
sg_symmetry(32, 221) = 'z,-x,y'
sg_symmetry(33, 221) = '-y,-z,-x'
sg_symmetry(34, 221) = 'y,-z,x'
sg_symmetry(35, 221) = '-y,z,x'
sg_symmetry(36, 221) = 'y,z,-x'
sg_symmetry(37, 221) = '-y,-x,z'
sg_symmetry(38, 221) = 'y,x,z'
sg_symmetry(39, 221) = '-y,x,-z'
sg_symmetry(40, 221) = 'y,-x,-z'
sg_symmetry(41, 221) = '-x,-z,y'
sg_symmetry(42, 221) = 'x,-z,-y'
sg_symmetry(43, 221) = 'x,z,y'
sg_symmetry(44, 221) = '-x,z,-y'
sg_symmetry(45, 221) = '-z,-y,x'
sg_symmetry(46, 221) = '-z,y,-x'
sg_symmetry(47, 221) = 'z,-y,-x'
sg_symmetry(48, 221) = 'z,y,x'
sg_name(222) = 'Pn-3n'
sg_patn(222) = 221
sg_symnum(222) = 48
sg_symmetry(1, 222) = 'x,y,z'
sg_symmetry(2, 222) = '-x,-y,z'
sg_symmetry(3, 222) = '-x,y,-z'
sg_symmetry(4, 222) = 'x,-y,-z'
sg_symmetry(5, 222) = 'z,x,y'
sg_symmetry(6, 222) = 'z,-x,-y'
sg_symmetry(7, 222) = '-z,-x,y'
sg_symmetry(8, 222) = '-z,x,-y'
sg_symmetry(9, 222) = 'y,z,x'
sg_symmetry(10, 222) = '-y,z,-x'
sg_symmetry(11, 222) = 'y,-z,-x'
sg_symmetry(12, 222) = '-y,-z,x'
sg_symmetry(13, 222) = 'y,x,-z'
sg_symmetry(14, 222) = '-y,-x,-z'
sg_symmetry(15, 222) = 'y,-x,z'
sg_symmetry(16, 222) = '-y,x,z'
sg_symmetry(17, 222) = 'x,z,-y'
sg_symmetry(18, 222) = '-x,z,y'
sg_symmetry(19, 222) = '-x,-z,-y'
sg_symmetry(20, 222) = 'x,-z,y'
sg_symmetry(21, 222) = 'z,y,-x'
sg_symmetry(22, 222) = 'z,-y,x'
sg_symmetry(23, 222) = '-z,y,x'
sg_symmetry(24, 222) = '-z,-y,-x'
sg_symmetry(25, 222) = '-x,-y,-z'
sg_symmetry(26, 222) = 'x,y,-z'
sg_symmetry(27, 222) = 'x,-y,z'
sg_symmetry(28, 222) = '-x,y,z'
sg_symmetry(29, 222) = '-z,-x,-y'
sg_symmetry(30, 222) = '-z,x,y'
sg_symmetry(31, 222) = 'z,x,-y'
sg_symmetry(32, 222) = 'z,-x,y'
sg_symmetry(33, 222) = '-y,-z,-x'
sg_symmetry(34, 222) = 'y,-z,x'
sg_symmetry(35, 222) = '-y,z,x'
sg_symmetry(36, 222) = 'y,z,-x'
sg_symmetry(37, 222) = '-y,-x,z'
sg_symmetry(38, 222) = 'y,x,z'
sg_symmetry(39, 222) = '-y,x,-z'
sg_symmetry(40, 222) = 'y,-x,-z'
sg_symmetry(41, 222) = '-x,-z,y'
sg_symmetry(42, 222) = 'x,-z,-y'
sg_symmetry(43, 222) = 'x,z,y'
sg_symmetry(44, 222) = '-x,z,-y'
sg_symmetry(45, 222) = '-z,-y,x'
sg_symmetry(46, 222) = '-z,y,-x'
sg_symmetry(47, 222) = 'z,-y,-x'
sg_symmetry(48, 222) = 'z,y,x'
sg_name(223) = 'Pm-3n'
sg_patn(223) = 221
sg_symnum(223) = 48
sg_symmetry(1, 223) = 'x,y,z'
sg_symmetry(2, 223) = '-x,-y,z'
sg_symmetry(3, 223) = '-x,y,-z'
sg_symmetry(4, 223) = 'x,-y,-z'
sg_symmetry(5, 223) = 'z,x,y'
sg_symmetry(6, 223) = 'z,-x,-y'
sg_symmetry(7, 223) = '-z,-x,y'
sg_symmetry(8, 223) = '-z,x,-y'
sg_symmetry(9, 223) = 'y,z,x'
sg_symmetry(10, 223) = '-y,z,-x'
sg_symmetry(11, 223) = 'y,-z,-x'
sg_symmetry(12, 223) = '-y,-z,x'
sg_symmetry(13, 223) = '1/2+y,1/2+x,1/2-z'
sg_symmetry(14, 223) = '1/2-y,1/2-x,1/2-z'
sg_symmetry(15, 223) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(16, 223) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(17, 223) = '1/2+x,1/2+z,1/2-y'
sg_symmetry(18, 223) = '1/2-x,1/2+z,1/2+y'
sg_symmetry(19, 223) = '1/2-x,1/2-z,1/2-y'
sg_symmetry(20, 223) = '1/2+x,1/2-z,1/2+y'
sg_symmetry(21, 223) = '1/2+z,1/2+y,1/2-x'
sg_symmetry(22, 223) = '1/2+z,1/2-y,1/2+x'
sg_symmetry(23, 223) = '1/2-z,1/2+y,1/2+x'
sg_symmetry(24, 223) = '1/2-z,1/2-y,1/2-x'
sg_symmetry(25, 223) = '-x,-y,-z'
sg_symmetry(26, 223) = 'x,y,-z'
sg_symmetry(27, 223) = 'x,-y,z'
sg_symmetry(28, 223) = '-x,y,z'
sg_symmetry(29, 223) = '-z,-x,-y'
sg_symmetry(30, 223) = '-z,x,y'
sg_symmetry(31, 223) = 'z,x,-y'
sg_symmetry(32, 223) = 'z,-x,y'
sg_symmetry(33, 223) = '-y,-z,-x'
sg_symmetry(34, 223) = 'y,-z,x'
sg_symmetry(35, 223) = '-y,z,x'
sg_symmetry(36, 223) = 'y,z,-x'
sg_symmetry(37, 223) = '1/2-y,1/2-x,1/2+z'
sg_symmetry(38, 223) = '1/2+y,1/2+x,1/2+z'
sg_symmetry(39, 223) = '1/2-y,1/2+x,1/2-z'
sg_symmetry(40, 223) = '1/2+y,1/2-x,1/2-z'
sg_symmetry(41, 223) = '1/2-x,1/2-z,1/2+y'
sg_symmetry(42, 223) = '1/2+x,1/2-z,1/2-y'
sg_symmetry(43, 223) = '1/2+x,1/2+z,1/2+y'
sg_symmetry(44, 223) = '1/2-x,1/2+z,1/2-y'
sg_symmetry(45, 223) = '1/2-z,1/2-y,1/2+x'
sg_symmetry(46, 223) = '1/2-z,1/2+y,1/2-x'
sg_symmetry(47, 223) = '1/2+z,1/2-y,1/2-x'
sg_symmetry(48, 223) = '1/2+z,1/2+y,1/2+x'
sg_name(224) = 'Pn-3m'
sg_patn(224) = 221
sg_symnum(224) = 48
sg_symmetry(1, 224) = 'x,y,z'
sg_symmetry(2, 224) = '-x,-y,z'
sg_symmetry(3, 224) = '-x,y,-z'
sg_symmetry(4, 224) = 'x,-y,-z'
sg_symmetry(5, 224) = 'z,x,y'
sg_symmetry(6, 224) = 'z,-x,-y'
sg_symmetry(7, 224) = '-z,-x,y'
sg_symmetry(8, 224) = '-z,x,-y'
sg_symmetry(9, 224) = 'y,z,x'
sg_symmetry(10, 224) = '-y,z,-x'
sg_symmetry(11, 224) = 'y,-z,-x'
sg_symmetry(12, 224) = '-y,-z,x'
sg_symmetry(13, 224) = '1/2+y,1/2+x,1/2-z'
sg_symmetry(14, 224) = '1/2-y,1/2-x,1/2-z'
sg_symmetry(15, 224) = '1/2+y,1/2-x,1/2+z'
sg_symmetry(16, 224) = '1/2-y,1/2+x,1/2+z'
sg_symmetry(17, 224) = '1/2+x,1/2+z,1/2-y'
sg_symmetry(18, 224) = '1/2-x,1/2+z,1/2+y'
sg_symmetry(19, 224) = '1/2-x,1/2-z,1/2-y'
sg_symmetry(20, 224) = '1/2+x,1/2-z,1/2+y'
sg_symmetry(21, 224) = '1/2+z,1/2+y,1/2-x'
sg_symmetry(22, 224) = '1/2+z,1/2-y,1/2+x'
sg_symmetry(23, 224) = '1/2-z,1/2+y,1/2+x'
sg_symmetry(24, 224) = '1/2-z,1/2-y,1/2-x'
sg_symmetry(25, 224) = '1/2-x,1/2-y,1/2-z'
sg_symmetry(26, 224) = '1/2+x,1/2+y,1/2-z'
sg_symmetry(27, 224) = '1/2+x,1/2-y,1/2+z'
sg_symmetry(28, 224) = '1/2-x,1/2+y,1/2+z'
sg_symmetry(29, 224) = '1/2-z,1/2-x,1/2-y'
sg_symmetry(30, 224) = '1/2-z,1/2+x,1/2+y'
sg_symmetry(31, 224) = '1/2+z,1/2+x,1/2-y'
sg_symmetry(32, 224) = '1/2+z,1/2-x,1/2+y'
sg_symmetry(33, 224) = '1/2-y,1/2-z,1/2-x'
sg_symmetry(34, 224) = '1/2+y,1/2-z,1/2+x'
sg_symmetry(35, 224) = '1/2-y,1/2+z,1/2+x'
sg_symmetry(36, 224) = '1/2+y,1/2+z,1/2-x'
sg_symmetry(37, 224) = '-y,-x,z'
sg_symmetry(38, 224) = 'y,x,z'
sg_symmetry(39, 224) = '-y,x,-z'
sg_symmetry(40, 224) = 'y,-x,-z'
sg_symmetry(41, 224) = '-x,-z,y'
sg_symmetry(42, 224) = 'x,-z,-y'
sg_symmetry(43, 224) = 'x,z,y'
sg_symmetry(44, 224) = '-x,z,-y'
sg_symmetry(45, 224) = '-z,-y,x'
sg_symmetry(46, 224) = '-z,y,-x'
sg_symmetry(47, 224) = 'z,-y,-x'
sg_symmetry(48, 224) = 'z,y,x'
sg_name(225) = 'Fm-3m'
sg_patn(225) = 225
sg_symnum(225) = 192
sg_symmetry(1, 225) = 'x,y,z'
sg_symmetry(2, 225) = '-x,-y,z'
sg_symmetry(3, 225) = '-x,y,-z'
sg_symmetry(4, 225) = 'x,-y,-z'
sg_symmetry(5, 225) = 'z,x,y'
sg_symmetry(6, 225) = 'z,-x,-y'
sg_symmetry(7, 225) = '-z,-x,y'
sg_symmetry(8, 225) = '-z,x,-y'
sg_symmetry(9, 225) = 'y,z,x'
sg_symmetry(10, 225) = '-y,z,-x'
sg_symmetry(11, 225) = 'y,-z,-x'
sg_symmetry(12, 225) = '-y,-z,x'
sg_symmetry(13, 225) = 'y,x,-z'
sg_symmetry(14, 225) = '-y,-x,-z'
sg_symmetry(15, 225) = 'y,-x,z'
sg_symmetry(16, 225) = '-y,x,z'
sg_symmetry(17, 225) = 'x,z,-y'
sg_symmetry(18, 225) = '-x,z,y'
sg_symmetry(19, 225) = '-x,-z,-y'
sg_symmetry(20, 225) = 'x,-z,y'
sg_symmetry(21, 225) = 'z,y,-x'
sg_symmetry(22, 225) = 'z,-y,x'
sg_symmetry(23, 225) = '-z,y,x'
sg_symmetry(24, 225) = '-z,-y,-x'
sg_symmetry(25, 225) = '-x,-y,-z'
sg_symmetry(26, 225) = 'x,y,-z'
sg_symmetry(27, 225) = 'x,-y,z'
sg_symmetry(28, 225) = '-x,y,z'
sg_symmetry(29, 225) = '-z,-x,-y'
sg_symmetry(30, 225) = '-z,x,y'
sg_symmetry(31, 225) = 'z,x,-y'
sg_symmetry(32, 225) = 'z,-x,y'
sg_symmetry(33, 225) = '-y,-z,-x'
sg_symmetry(34, 225) = 'y,-z,x'
sg_symmetry(35, 225) = '-y,z,x'
sg_symmetry(36, 225) = 'y,z,-x'
sg_symmetry(37, 225) = '-y,-x,z'
sg_symmetry(38, 225) = 'y,x,z'
sg_symmetry(39, 225) = '-y,x,-z'
sg_symmetry(40, 225) = 'y,-x,-z'
sg_symmetry(41, 225) = '-x,-z,y'
sg_symmetry(42, 225) = 'x,-z,-y'
sg_symmetry(43, 225) = 'x,z,y'
sg_symmetry(44, 225) = '-x,z,-y'
sg_symmetry(45, 225) = '-z,-y,x'
sg_symmetry(46, 225) = '-z,y,-x'
sg_symmetry(47, 225) = 'z,-y,-x'
sg_symmetry(48, 225) = 'z,y,x'
sg_symmetry(49, 225) = 'x,1/2+y,1/2+z'
sg_symmetry(50, 225) = '-x,1/2-y,1/2+z'
sg_symmetry(51, 225) = '-x,1/2+y,1/2-z'
sg_symmetry(52, 225) = 'x,1/2-y,1/2-z'
sg_symmetry(53, 225) = 'z,1/2+x,1/2+y'
sg_symmetry(54, 225) = 'z,1/2-x,1/2-y'
sg_symmetry(55, 225) = '-z,1/2-x,1/2+y'
sg_symmetry(56, 225) = '-z,1/2+x,1/2-y'
sg_symmetry(57, 225) = 'y,1/2+z,1/2+x'
sg_symmetry(58, 225) = '-y,1/2+z,1/2-x'
sg_symmetry(59, 225) = 'y,1/2-z,1/2-x'
sg_symmetry(60, 225) = '-y,1/2-z,1/2+x'
sg_symmetry(61, 225) = 'y,1/2+x,1/2-z'
sg_symmetry(62, 225) = '-y,1/2-x,1/2-z'
sg_symmetry(63, 225) = 'y,1/2-x,1/2+z'
sg_symmetry(64, 225) = '-y,1/2+x,1/2+z'
sg_symmetry(65, 225) = 'x,1/2+z,1/2-y'
sg_symmetry(66, 225) = '-x,1/2+z,1/2+y'
sg_symmetry(67, 225) = '-x,1/2-z,1/2-y'
sg_symmetry(68, 225) = 'x,1/2-z,1/2+y'
sg_symmetry(69, 225) = 'z,1/2+y,1/2-x'
sg_symmetry(70, 225) = 'z,1/2-y,1/2+x'
sg_symmetry(71, 225) = '-z,1/2+y,1/2+x'
sg_symmetry(72, 225) = '-z,1/2-y,1/2-x'
sg_symmetry(73, 225) = '-x,1/2-y,1/2-z'
sg_symmetry(74, 225) = 'y,-x+1/2,-z+1/2'
sg_symmetry(75, 225) = 'x,y+1/2,-z+1/2'
sg_symmetry(76, 225) = '-y,x+1/2,-z+1/2'
sg_symmetry(77, 225) = '-x,y+1/2,z+1/2'
sg_symmetry(78, 225) = '-y,-x+1/2,z+1/2'
sg_symmetry(79, 225) = 'x,-y+1/2,z+1/2'
sg_symmetry(80, 225) = 'y,x+1/2,z+1/2'
sg_symmetry(81, 225) = '-z,-x+1/2,-y+1/2'
sg_symmetry(82, 225) = 'x,-z+1/2,-y+1/2'
sg_symmetry(83, 225) = 'z,x+1/2,-y+1/2'
sg_symmetry(84, 225) = '-x,z+1/2,-y+1/2'
sg_symmetry(85, 225) = '-z,x+1/2,y+1/2'
sg_symmetry(86, 225) = '-x,-z+1/2,y+1/2'
sg_symmetry(87, 225) = 'z,-x+1/2,y+1/2'
sg_symmetry(88, 225) = 'x,z+1/2,y+1/2'
sg_symmetry(89, 225) = '-y,-z+1/2,-x+1/2'
sg_symmetry(90, 225) = '-y,z+1/2,x+1/2'
sg_symmetry(91, 225) = '-z,-y+1/2,x+1/2'
sg_symmetry(92, 225) = 'y,-z+1/2,x+1/2'
sg_symmetry(93, 225) = 'z,y+1/2,x+1/2'
sg_symmetry(94, 225) = 'y,z+1/2,-x+1/2'
sg_symmetry(95, 225) = '-z,y+1/2,-x+1/2'
sg_symmetry(96, 225) = 'z,-y+1/2,-x+1/2'
sg_symmetry(97, 225) = 'x+1/2,y,z+1/2'
sg_symmetry(98, 225) = '-y+1/2,x,z+1/2'
sg_symmetry(99, 225) = '-x+1/2,-y,z+1/2'
sg_symmetry(100, 225) = 'y+1/2,-x,z+1/2'
sg_symmetry(101, 225) = 'x+1/2,-y,-z+1/2'
sg_symmetry(102, 225) = 'y+1/2,x,-z+1/2'
sg_symmetry(103, 225) = '-x+1/2,y,-z+1/2'
sg_symmetry(104, 225) = '-y+1/2,-x,-z+1/2'
sg_symmetry(105, 225) = 'z+1/2,x,y+1/2'
sg_symmetry(106, 225) = '-x+1/2,z,y+1/2'
sg_symmetry(107, 225) = '-z+1/2,-x,y+1/2'
sg_symmetry(108, 225) = 'x+1/2,-z,y+1/2'
sg_symmetry(109, 225) = 'z+1/2,-x,-y+1/2'
sg_symmetry(110, 225) = 'x+1/2,z,-y+1/2'
sg_symmetry(111, 225) = '-z+1/2,x,-y+1/2'
sg_symmetry(112, 225) = '-x+1/2,-z,-y+1/2'
sg_symmetry(113, 225) = 'y+1/2,z,x+1/2'
sg_symmetry(114, 225) = 'y+1/2,-z,-x+1/2'
sg_symmetry(115, 225) = 'z+1/2,y,-x+1/2'
sg_symmetry(116, 225) = '-y+1/2,z,-x+1/2'
sg_symmetry(117, 225) = '-z+1/2,-y,-x+1/2'
sg_symmetry(118, 225) = '-y+1/2,-z,x+1/2'
sg_symmetry(119, 225) = 'z+1/2,-y,x+1/2'
sg_symmetry(120, 225) = '-z+1/2,y,x+1/2'
sg_symmetry(121, 225) = '-x+1/2,-y,-z+1/2'
sg_symmetry(122, 225) = 'y+1/2,-x,-z+1/2'
sg_symmetry(123, 225) = 'x+1/2,y,-z+1/2'
sg_symmetry(124, 225) = '-y+1/2,x,-z+1/2'
sg_symmetry(125, 225) = '-x+1/2,y,z+1/2'
sg_symmetry(126, 225) = '-y+1/2,-x,z+1/2'
sg_symmetry(127, 225) = 'x+1/2,-y,z+1/2'
sg_symmetry(128, 225) = 'y+1/2,x,z+1/2'
sg_symmetry(129, 225) = '-z+1/2,-x,-y+1/2'
sg_symmetry(130, 225) = 'x+1/2,-z,-y+1/2'
sg_symmetry(131, 225) = 'z+1/2,x,-y+1/2'
sg_symmetry(132, 225) = '-x+1/2,z,-y+1/2'
sg_symmetry(133, 225) = '-z+1/2,x,y+1/2'
sg_symmetry(134, 225) = '-x+1/2,-z,y+1/2'
sg_symmetry(135, 225) = 'z+1/2,-x,y+1/2'
sg_symmetry(136, 225) = 'x+1/2,z,y+1/2'
sg_symmetry(137, 225) = '-y+1/2,-z,-x+1/2'
sg_symmetry(138, 225) = '-y+1/2,z,x+1/2'
sg_symmetry(139, 225) = '-z+1/2,-y,x+1/2'
sg_symmetry(140, 225) = 'y+1/2,-z,x+1/2'
sg_symmetry(141, 225) = 'z+1/2,y,x+1/2'
sg_symmetry(142, 225) = 'y+1/2,z,-x+1/2'
sg_symmetry(143, 225) = '-z+1/2,y,-x+1/2'
sg_symmetry(144, 225) = 'z+1/2,-y,-x+1/2'
sg_symmetry(145, 225) = 'x+1/2,y+1/2,z'
sg_symmetry(146, 225) = '-y+1/2,x+1/2,z'
sg_symmetry(147, 225) = '-x+1/2,-y+1/2,z'
sg_symmetry(148, 225) = 'y+1/2,-x+1/2,z'
sg_symmetry(149, 225) = 'x+1/2,-y+1/2,-z'
sg_symmetry(150, 225) = 'y+1/2,x+1/2,-z'
sg_symmetry(151, 225) = '-x+1/2,y+1/2,-z'
sg_symmetry(152, 225) = '-y+1/2,-x+1/2,-z'
sg_symmetry(153, 225) = 'z+1/2,x+1/2,y'
sg_symmetry(154, 225) = '-x+1/2,z+1/2,y'
sg_symmetry(155, 225) = '-z+1/2,-x+1/2,y'
sg_symmetry(156, 225) = 'x+1/2,-z+1/2,y'
sg_symmetry(157, 225) = 'z+1/2,-x+1/2,-y'
sg_symmetry(158, 225) = 'x+1/2,z+1/2,-y'
sg_symmetry(159, 225) = '-z+1/2,x+1/2,-y'
sg_symmetry(160, 225) = '-x+1/2,-z+1/2,-y'
sg_symmetry(161, 225) = 'y+1/2,z+1/2,x'
sg_symmetry(162, 225) = 'y+1/2,-z+1/2,-x'
sg_symmetry(163, 225) = 'z+1/2,y+1/2,-x'
sg_symmetry(164, 225) = '-y+1/2,z+1/2,-x'
sg_symmetry(165, 225) = '-z+1/2,-y+1/2,-x'
sg_symmetry(166, 225) = '-y+1/2,-z+1/2,x'
sg_symmetry(167, 225) = 'z+1/2,-y+1/2,x'
sg_symmetry(168, 225) = '-z+1/2,y+1/2,x'
sg_symmetry(169, 225) = '-x+1/2,-y+1/2,-z'
sg_symmetry(170, 225) = 'y+1/2,-x+1/2,-z'
sg_symmetry(171, 225) = 'x+1/2,y+1/2,-z'
sg_symmetry(172, 225) = '-y+1/2,x+1/2,-z'
sg_symmetry(173, 225) = '-x+1/2,y+1/2,z'
sg_symmetry(174, 225) = '-y+1/2,-x+1/2,z'
sg_symmetry(175, 225) = 'x+1/2,-y+1/2,z'
sg_symmetry(176, 225) = 'y+1/2,x+1/2,z'
sg_symmetry(177, 225) = '-z+1/2,-x+1/2,-y'
sg_symmetry(178, 225) = 'x+1/2,-z+1/2,-y'
sg_symmetry(179, 225) = 'z+1/2,x+1/2,-y'
sg_symmetry(180, 225) = '-x+1/2,z+1/2,-y'
sg_symmetry(181, 225) = '-z+1/2,x+1/2,y'
sg_symmetry(182, 225) = '-x+1/2,-z+1/2,y'
sg_symmetry(183, 225) = 'z+1/2,-x+1/2,y'
sg_symmetry(184, 225) = 'x+1/2,z+1/2,y'
sg_symmetry(185, 225) = '-y+1/2,-z+1/2,-x'
sg_symmetry(186, 225) = '-y+1/2,z+1/2,x'
sg_symmetry(187, 225) = '-z+1/2,-y+1/2,x'
sg_symmetry(188, 225) = 'y+1/2,-z+1/2,x'
sg_symmetry(189, 225) = 'z+1/2,y+1/2,x'
sg_symmetry(190, 225) = 'y+1/2,z+1/2,-x'
sg_symmetry(191, 225) = '-z+1/2,y+1/2,-x'
sg_symmetry(192, 225) = 'z+1/2,-y+1/2,-x'
sg_name(226) = 'Fm-3c'
sg_patn(226) = 226
sg_symnum(226) = 192
sg_symmetry(1, 226) = 'x,y,z'
sg_symmetry(2, 226) = '-y+1/2,x,z'
sg_symmetry(3, 226) = '-x+1/2,-y+1/2,z'
sg_symmetry(4, 226) = 'y,-x+1/2,z'
sg_symmetry(5, 226) = 'x,-y,-z'
sg_symmetry(6, 226) = 'y+1/2,x,-z'
sg_symmetry(7, 226) = '-x+1/2,y+1/2,-z'
sg_symmetry(8, 226) = '-y,-x+1/2,-z'
sg_symmetry(9, 226) = 'z,x,y'
sg_symmetry(10, 226) = '-x+1/2,z,y'
sg_symmetry(11, 226) = '-z+1/2,-x+1/2,y'
sg_symmetry(12, 226) = 'x,-z+1/2,y'
sg_symmetry(13, 226) = 'z,-x,-y'
sg_symmetry(14, 226) = 'x+1/2,z,-y'
sg_symmetry(15, 226) = '-z+1/2,x+1/2,-y'
sg_symmetry(16, 226) = '-x,-z+1/2,-y'
sg_symmetry(17, 226) = 'y,z,x'
sg_symmetry(18, 226) = 'y,-z+1/2,-x+1/2'
sg_symmetry(19, 226) = 'z,y,-x+1/2'
sg_symmetry(20, 226) = '-y+1/2,z,-x+1/2'
sg_symmetry(21, 226) = '-z+1/2,-y,-x'
sg_symmetry(22, 226) = '-y,-z,x'
sg_symmetry(23, 226) = 'z,-y,x+1/2'
sg_symmetry(24, 226) = '-z+1/2,y+1/2,x+1/2'
sg_symmetry(25, 226) = '-x,-y,-z'
sg_symmetry(26, 226) = 'y-1/2,-x,-z'
sg_symmetry(27, 226) = 'x-1/2,y-1/2,-z'
sg_symmetry(28, 226) = '-y,x-1/2,-z'
sg_symmetry(29, 226) = '-x,y,z'
sg_symmetry(30, 226) = '-y-1/2,-x,z'
sg_symmetry(31, 226) = 'x-1/2,-y-1/2,z'
sg_symmetry(32, 226) = 'y,x-1/2,z'
sg_symmetry(33, 226) = '-z,-x,-y'
sg_symmetry(34, 226) = 'x-1/2,-z,-y'
sg_symmetry(35, 226) = 'z-1/2,x-1/2,-y'
sg_symmetry(36, 226) = '-x,z-1/2,-y'
sg_symmetry(37, 226) = '-z,x,y'
sg_symmetry(38, 226) = '-x-1/2,-z,y'
sg_symmetry(39, 226) = 'z-1/2,-x-1/2,y'
sg_symmetry(40, 226) = 'x,z-1/2,y'
sg_symmetry(41, 226) = '-y,-z,-x'
sg_symmetry(42, 226) = '-y,z-1/2,x-1/2'
sg_symmetry(43, 226) = '-z,-y,x-1/2'
sg_symmetry(44, 226) = 'y-1/2,-z,x-1/2'
sg_symmetry(45, 226) = 'z-1/2,y,x'
sg_symmetry(46, 226) = 'y,z,-x'
sg_symmetry(47, 226) = '-z,y,-x-1/2'
sg_symmetry(48, 226) = 'z-1/2,-y-1/2,-x-1/2'
sg_symmetry(49, 226) = 'x,y+1/2,z+1/2'
sg_symmetry(50, 226) = '-y+1/2,x+1/2,z+1/2'
sg_symmetry(51, 226) = '-x+1/2,-y+1,z+1/2'
sg_symmetry(52, 226) = 'y,-x+1,z+1/2'
sg_symmetry(53, 226) = 'x,-y+1/2,-z+1/2'
sg_symmetry(54, 226) = 'y+1/2,x+1/2,-z+1/2'
sg_symmetry(55, 226) = '-x+1/2,y+1,-z+1/2'
sg_symmetry(56, 226) = '-y,-x+1,-z+1/2'
sg_symmetry(57, 226) = 'z,x+1/2,y+1/2'
sg_symmetry(58, 226) = '-x+1/2,z+1/2,y+1/2'
sg_symmetry(59, 226) = '-z+1/2,-x+1,y+1/2'
sg_symmetry(60, 226) = 'x,-z+1,y+1/2'
sg_symmetry(61, 226) = 'z,-x+1/2,-y+1/2'
sg_symmetry(62, 226) = 'x+1/2,z+1/2,-y+1/2'
sg_symmetry(63, 226) = '-z+1/2,x+1,-y+1/2'
sg_symmetry(64, 226) = '-x,-z+1,-y+1/2'
sg_symmetry(65, 226) = 'y,z+1/2,x+1/2'
sg_symmetry(66, 226) = 'y,-z+1,-x+1'
sg_symmetry(67, 226) = 'z,y+1/2,-x+1'
sg_symmetry(68, 226) = '-y+1/2,z+1/2,-x+1'
sg_symmetry(69, 226) = '-z+1/2,-y+1/2,-x+1/2'
sg_symmetry(70, 226) = '-y,-z+1/2,x+1/2'
sg_symmetry(71, 226) = 'z,-y+1/2,x+1'
sg_symmetry(72, 226) = '-z+1/2,y+1,x+1'
sg_symmetry(73, 226) = '-x,-y+1/2,-z+1/2'
sg_symmetry(74, 226) = 'y-1/2,-x+1/2,-z+1/2'
sg_symmetry(75, 226) = 'x-1/2,y,-z+1/2'
sg_symmetry(76, 226) = '-y,x,-z+1/2'
sg_symmetry(77, 226) = '-x,y+1/2,z+1/2'
sg_symmetry(78, 226) = '-y-1/2,-x+1/2,z+1/2'
sg_symmetry(79, 226) = 'x-1/2,-y,z+1/2'
sg_symmetry(80, 226) = 'y,x,z+1/2'
sg_symmetry(81, 226) = '-z,-x+1/2,-y+1/2'
sg_symmetry(82, 226) = 'x-1/2,-z+1/2,-y+1/2'
sg_symmetry(83, 226) = 'z-1/2,x,-y+1/2'
sg_symmetry(84, 226) = '-x,z,-y+1/2'
sg_symmetry(85, 226) = '-z,x+1/2,y+1/2'
sg_symmetry(86, 226) = '-x-1/2,-z+1/2,y+1/2'
sg_symmetry(87, 226) = 'z-1/2,-x,y+1/2'
sg_symmetry(88, 226) = 'x,z,y+1/2'
sg_symmetry(89, 226) = '-y,-z+1/2,-x+1/2'
sg_symmetry(90, 226) = '-y,z,x'
sg_symmetry(91, 226) = '-z,-y+1/2,x'
sg_symmetry(92, 226) = 'y-1/2,-z+1/2,x'
sg_symmetry(93, 226) = 'z-1/2,y+1/2,x+1/2'
sg_symmetry(94, 226) = 'y,z+1/2,-x+1/2'
sg_symmetry(95, 226) = '-z,y+1/2,-x'
sg_symmetry(96, 226) = 'z-1/2,-y,-x'
sg_symmetry(97, 226) = 'x+1/2,y,z+1/2'
sg_symmetry(98, 226) = '-y+1,x,z+1/2'
sg_symmetry(99, 226) = '-x+1,-y+1/2,z+1/2'
sg_symmetry(100, 226) = 'y+1/2,-x+1/2,z+1/2'
sg_symmetry(101, 226) = 'x+1/2,-y,-z+1/2'
sg_symmetry(102, 226) = 'y+1,x,-z+1/2'
sg_symmetry(103, 226) = '-x+1,y+1/2,-z+1/2'
sg_symmetry(104, 226) = '-y+1/2,-x+1/2,-z+1/2'
sg_symmetry(105, 226) = 'z+1/2,x,y+1/2'
sg_symmetry(106, 226) = '-x+1,z,y+1/2'
sg_symmetry(107, 226) = '-z+1,-x+1/2,y+1/2'
sg_symmetry(108, 226) = 'x+1/2,-z+1/2,y+1/2'
sg_symmetry(109, 226) = 'z+1/2,-x,-y+1/2'
sg_symmetry(110, 226) = 'x+1,z,-y+1/2'
sg_symmetry(111, 226) = '-z+1,x+1/2,-y+1/2'
sg_symmetry(112, 226) = '-x+1/2,-z+1/2,-y+1/2'
sg_symmetry(113, 226) = 'y+1/2,z,x+1/2'
sg_symmetry(114, 226) = 'y+1/2,-z+1/2,-x+1'
sg_symmetry(115, 226) = 'z+1/2,y,-x+1'
sg_symmetry(116, 226) = '-y+1,z,-x+1'
sg_symmetry(117, 226) = '-z+1,-y,-x+1/2'
sg_symmetry(118, 226) = '-y+1/2,-z,x+1/2'
sg_symmetry(119, 226) = 'z+1/2,-y,x+1'
sg_symmetry(120, 226) = '-z+1,y+1/2,x+1'
sg_symmetry(121, 226) = '-x+1/2,-y,-z+1/2'
sg_symmetry(122, 226) = 'y,-x,-z+1/2'
sg_symmetry(123, 226) = 'x,y-1/2,-z+1/2'
sg_symmetry(124, 226) = '-y+1/2,x-1/2,-z+1/2'
sg_symmetry(125, 226) = '-x+1/2,y,z+1/2'
sg_symmetry(126, 226) = '-y,-x,z+1/2'
sg_symmetry(127, 226) = 'x,-y-1/2,z+1/2'
sg_symmetry(128, 226) = 'y+1/2,x-1/2,z+1/2'
sg_symmetry(129, 226) = '-z+1/2,-x,-y+1/2'
sg_symmetry(130, 226) = 'x,-z,-y+1/2'
sg_symmetry(131, 226) = 'z,x-1/2,-y+1/2'
sg_symmetry(132, 226) = '-x+1/2,z-1/2,-y+1/2'
sg_symmetry(133, 226) = '-z+1/2,x,y+1/2'
sg_symmetry(134, 226) = '-x,-z,y+1/2'
sg_symmetry(135, 226) = 'z,-x-1/2,y+1/2'
sg_symmetry(136, 226) = 'x+1/2,z-1/2,y+1/2'
sg_symmetry(137, 226) = '-y+1/2,-z,-x+1/2'
sg_symmetry(138, 226) = '-y+1/2,z-1/2,x'
sg_symmetry(139, 226) = '-z+1/2,-y,x'
sg_symmetry(140, 226) = 'y,-z,x'
sg_symmetry(141, 226) = 'z,y,x+1/2'
sg_symmetry(142, 226) = 'y+1/2,z,-x+1/2'
sg_symmetry(143, 226) = '-z+1/2,y,-x'
sg_symmetry(144, 226) = 'z,-y-1/2,-x'
sg_symmetry(145, 226) = 'x+1/2,y+1/2,z'
sg_symmetry(146, 226) = '-y+1,x+1/2,z'
sg_symmetry(147, 226) = '-x+1,-y+1,z'
sg_symmetry(148, 226) = 'y+1/2,-x+1,z'
sg_symmetry(149, 226) = 'x+1/2,-y+1/2,-z'
sg_symmetry(150, 226) = 'y+1,x+1/2,-z'
sg_symmetry(151, 226) = '-x+1,y+1,-z'
sg_symmetry(152, 226) = '-y+1/2,-x+1,-z'
sg_symmetry(153, 226) = 'z+1/2,x+1/2,y'
sg_symmetry(154, 226) = '-x+1,z+1/2,y'
sg_symmetry(155, 226) = '-z+1,-x+1,y'
sg_symmetry(156, 226) = 'x+1/2,-z+1,y'
sg_symmetry(157, 226) = 'z+1/2,-x+1/2,-y'
sg_symmetry(158, 226) = 'x+1,z+1/2,-y'
sg_symmetry(159, 226) = '-z+1,x+1,-y'
sg_symmetry(160, 226) = '-x+1/2,-z+1,-y'
sg_symmetry(161, 226) = 'y+1/2,z+1/2,x'
sg_symmetry(162, 226) = 'y+1/2,-z+1,-x+1/2'
sg_symmetry(163, 226) = 'z+1/2,y+1/2,-x+1/2'
sg_symmetry(164, 226) = '-y+1,z+1/2,-x+1/2'
sg_symmetry(165, 226) = '-z+1,-y+1/2,-x'
sg_symmetry(166, 226) = '-y+1/2,-z+1/2,x'
sg_symmetry(167, 226) = 'z+1/2,-y+1/2,x+1/2'
sg_symmetry(168, 226) = '-z+1,y+1,x+1/2'
sg_symmetry(169, 226) = '-x+1/2,-y+1/2,-z'
sg_symmetry(170, 226) = 'y,-x+1/2,-z'
sg_symmetry(171, 226) = 'x,y,-z'
sg_symmetry(172, 226) = '-y+1/2,x,-z'
sg_symmetry(173, 226) = '-x+1/2,y+1/2,z'
sg_symmetry(174, 226) = '-y,-x+1/2,z'
sg_symmetry(175, 226) = 'x,-y,z'
sg_symmetry(176, 226) = 'y+1/2,x,z'
sg_symmetry(177, 226) = '-z+1/2,-x+1/2,-y'
sg_symmetry(178, 226) = 'x,-z+1/2,-y'
sg_symmetry(179, 226) = 'z,x,-y'
sg_symmetry(180, 226) = '-x+1/2,z,-y'
sg_symmetry(181, 226) = '-z+1/2,x+1/2,y'
sg_symmetry(182, 226) = '-x,-z+1/2,y'
sg_symmetry(183, 226) = 'z,-x,y'
sg_symmetry(184, 226) = 'x+1/2,z,y'
sg_symmetry(185, 226) = '-y+1/2,-z+1/2,-x'
sg_symmetry(186, 226) = '-y+1/2,z,x-1/2'
sg_symmetry(187, 226) = '-z+1/2,-y+1/2,x-1/2'
sg_symmetry(188, 226) = 'y,-z+1/2,x-1/2'
sg_symmetry(189, 226) = 'z,y+1/2,x'
sg_symmetry(190, 226) = 'y+1/2,z+1/2,-x'
sg_symmetry(191, 226) = '-z+1/2,y+1/2,-x-1/2'
sg_symmetry(192, 226) = 'z,-y,-x-1/2'
sg_name(227) = 'Fd-3m'
sg_patn(227) = 227
sg_symnum(227) = 192
sg_symmetry(1, 227) = 'x,y,z'
sg_symmetry(2, 227) = '-y,x+1/4,z+1/4'
sg_symmetry(3, 227) = '-x+3/4,-y+1/4,z+1/2'
sg_symmetry(4, 227) = 'y+3/4,-x,z+3/4'
sg_symmetry(5, 227) = 'x,-y+1/4,-z+1/4'
sg_symmetry(6, 227) = 'y+3/4,x+1/4,-z+1/2'
sg_symmetry(7, 227) = '-x+3/4,y,-z+3/4'
sg_symmetry(8, 227) = '-y,-x,-z'
sg_symmetry(9, 227) = 'z,x,y'
sg_symmetry(10, 227) = '-x,z+1/4,y+1/4'
sg_symmetry(11, 227) = '-z+3/4,-x+1/4,y+1/2'
sg_symmetry(12, 227) = 'x+3/4,-z,y+3/4'
sg_symmetry(13, 227) = 'z,-x+1/4,-y+1/4'
sg_symmetry(14, 227) = 'x+3/4,z+1/4,-y+1/2'
sg_symmetry(15, 227) = '-z+3/4,x,-y+3/4'
sg_symmetry(16, 227) = '-x,-z,-y'
sg_symmetry(17, 227) = 'y,z,x'
sg_symmetry(18, 227) = 'y+1/2,-z+3/4,-x+1/4'
sg_symmetry(19, 227) = 'z+1/4,y+3/4,-x+1/2'
sg_symmetry(20, 227) = '-y+1/4,z+1/2,-x+3/4'
sg_symmetry(21, 227) = '-z,-y+1/2,-x+1/2'
sg_symmetry(22, 227) = '-y+1/4,-z+1/4,x'
sg_symmetry(23, 227) = 'z+1/4,-y,x+1/4'
sg_symmetry(24, 227) = '-z+1/2,y+1/4,x+3/4'
sg_symmetry(25, 227) = '-x,-y,-z'
sg_symmetry(26, 227) = 'y,-x-1/4,-z-1/4'
sg_symmetry(27, 227) = 'x-3/4,y-1/4,-z-1/2'
sg_symmetry(28, 227) = '-y-3/4,x,-z-3/4'
sg_symmetry(29, 227) = '-x,y-1/4,z-1/4'
sg_symmetry(30, 227) = '-y-3/4,-x-1/4,z-1/2'
sg_symmetry(31, 227) = 'x-3/4,-y,z-3/4'
sg_symmetry(32, 227) = 'y,x,z'
sg_symmetry(33, 227) = '-z,-x,-y'
sg_symmetry(34, 227) = 'x,-z-1/4,-y-1/4'
sg_symmetry(35, 227) = 'z-3/4,x-1/4,-y-1/2'
sg_symmetry(36, 227) = '-x-3/4,z,-y-3/4'
sg_symmetry(37, 227) = '-z,x-1/4,y-1/4'
sg_symmetry(38, 227) = '-x-3/4,-z-1/4,y-1/2'
sg_symmetry(39, 227) = 'z-3/4,-x,y-3/4'
sg_symmetry(40, 227) = 'x,z,y'
sg_symmetry(41, 227) = '-y,-z,-x'
sg_symmetry(42, 227) = '-y-1/2,z-3/4,x-1/4'
sg_symmetry(43, 227) = '-z-1/4,-y-3/4,x-1/2'
sg_symmetry(44, 227) = 'y-1/4,-z-1/2,x-3/4'
sg_symmetry(45, 227) = 'z,y-1/2,x-1/2'
sg_symmetry(46, 227) = 'y-1/4,z-1/4,-x'
sg_symmetry(47, 227) = '-z-1/4,y,-x-1/4'
sg_symmetry(48, 227) = 'z-1/2,-y-1/4,-x-3/4'
sg_symmetry(49, 227) = 'x,y+1/2,z+1/2'
sg_symmetry(50, 227) = '-y,x+3/4,z+3/4'
sg_symmetry(51, 227) = '-x+3/4,-y+3/4,z+1'
sg_symmetry(52, 227) = 'y+3/4,-x+1/2,z+5/4'
sg_symmetry(53, 227) = 'x,-y+3/4,-z+3/4'
sg_symmetry(54, 227) = 'y+3/4,x+3/4,-z+1'
sg_symmetry(55, 227) = '-x+3/4,y+1/2,-z+5/4'
sg_symmetry(56, 227) = '-y,-x+1/2,-z+1/2'
sg_symmetry(57, 227) = 'z,x+1/2,y+1/2'
sg_symmetry(58, 227) = '-x,z+3/4,y+3/4'
sg_symmetry(59, 227) = '-z+3/4,-x+3/4,y+1'
sg_symmetry(60, 227) = 'x+3/4,-z+1/2,y+5/4'
sg_symmetry(61, 227) = 'z,-x+3/4,-y+3/4'
sg_symmetry(62, 227) = 'x+3/4,z+3/4,-y+1'
sg_symmetry(63, 227) = '-z+3/4,x+1/2,-y+5/4'
sg_symmetry(64, 227) = '-x,-z+1/2,-y+1/2'
sg_symmetry(65, 227) = 'y,z+1/2,x+1/2'
sg_symmetry(66, 227) = 'y+1/2,-z+5/4,-x+3/4'
sg_symmetry(67, 227) = 'z+1/4,y+5/4,-x+1'
sg_symmetry(68, 227) = '-y+1/4,z+1,-x+5/4'
sg_symmetry(69, 227) = '-z,-y+1,-x+1'
sg_symmetry(70, 227) = '-y+1/4,-z+3/4,x+1/2'
sg_symmetry(71, 227) = 'z+1/4,-y+1/2,x+3/4'
sg_symmetry(72, 227) = '-z+1/2,y+3/4,x+5/4'
sg_symmetry(73, 227) = '-x,-y+1/2,-z+1/2'
sg_symmetry(74, 227) = 'y,-x+1/4,-z+1/4'
sg_symmetry(75, 227) = 'x-3/4,y+1/4,-z'
sg_symmetry(76, 227) = '-y-3/4,x+1/2,-z-1/4'
sg_symmetry(77, 227) = '-x,y+1/4,z+1/4'
sg_symmetry(78, 227) = '-y-3/4,-x+1/4,z'
sg_symmetry(79, 227) = 'x-3/4,-y+1/2,z-1/4'
sg_symmetry(80, 227) = 'y,x+1/2,z+1/2'
sg_symmetry(81, 227) = '-z,-x+1/2,-y+1/2'
sg_symmetry(82, 227) = 'x,-z+1/4,-y+1/4'
sg_symmetry(83, 227) = 'z-3/4,x+1/4,-y'
sg_symmetry(84, 227) = '-x-3/4,z+1/2,-y-1/4'
sg_symmetry(85, 227) = '-z,x+1/4,y+1/4'
sg_symmetry(86, 227) = '-x-3/4,-z+1/4,y'
sg_symmetry(87, 227) = 'z-3/4,-x+1/2,y-1/4'
sg_symmetry(88, 227) = 'x,z+1/2,y+1/2'
sg_symmetry(89, 227) = '-y,-z+1/2,-x+1/2'
sg_symmetry(90, 227) = '-y-1/2,z-1/4,x+1/4'
sg_symmetry(91, 227) = '-z-1/4,-y-1/4,x'
sg_symmetry(92, 227) = 'y-1/4,-z,x-1/4'
sg_symmetry(93, 227) = 'z,y,x'
sg_symmetry(94, 227) = 'y-1/4,z+1/4,-x+1/2'
sg_symmetry(95, 227) = '-z-1/4,y+1/2,-x+1/4'
sg_symmetry(96, 227) = 'z-1/2,-y+1/4,-x-1/4'
sg_symmetry(97, 227) = 'x+1/2,y,z+1/2'
sg_symmetry(98, 227) = '-y+1/2,x+1/4,z+3/4'
sg_symmetry(99, 227) = '-x+5/4,-y+1/4,z+1'
sg_symmetry(100, 227) = 'y+5/4,-x,z+5/4'
sg_symmetry(101, 227) = 'x+1/2,-y+1/4,-z+3/4'
sg_symmetry(102, 227) = 'y+5/4,x+1/4,-z+1'
sg_symmetry(103, 227) = '-x+5/4,y,-z+5/4'
sg_symmetry(104, 227) = '-y+1/2,-x,-z+1/2'
sg_symmetry(105, 227) = 'z+1/2,x,y+1/2'
sg_symmetry(106, 227) = '-x+1/2,z+1/4,y+3/4'
sg_symmetry(107, 227) = '-z+5/4,-x+1/4,y+1'
sg_symmetry(108, 227) = 'x+5/4,-z,y+5/4'
sg_symmetry(109, 227) = 'z+1/2,-x+1/4,-y+3/4'
sg_symmetry(110, 227) = 'x+5/4,z+1/4,-y+1'
sg_symmetry(111, 227) = '-z+5/4,x,-y+5/4'
sg_symmetry(112, 227) = '-x+1/2,-z,-y+1/2'
sg_symmetry(113, 227) = 'y+1/2,z,x+1/2'
sg_symmetry(114, 227) = 'y+1,-z+3/4,-x+3/4'
sg_symmetry(115, 227) = 'z+3/4,y+3/4,-x+1'
sg_symmetry(116, 227) = '-y+3/4,z+1/2,-x+5/4'
sg_symmetry(117, 227) = '-z+1/2,-y+1/2,-x+1'
sg_symmetry(118, 227) = '-y+3/4,-z+1/4,x+1/2'
sg_symmetry(119, 227) = 'z+3/4,-y,x+3/4'
sg_symmetry(120, 227) = '-z+1,y+1/4,x+5/4'
sg_symmetry(121, 227) = '-x+1/2,-y,-z+1/2'
sg_symmetry(122, 227) = 'y+1/2,-x-1/4,-z+1/4'
sg_symmetry(123, 227) = 'x-1/4,y-1/4,-z'
sg_symmetry(124, 227) = '-y-1/4,x,-z-1/4'
sg_symmetry(125, 227) = '-x+1/2,y-1/4,z+1/4'
sg_symmetry(126, 227) = '-y-1/4,-x-1/4,z'
sg_symmetry(127, 227) = 'x-1/4,-y,z-1/4'
sg_symmetry(128, 227) = 'y+1/2,x,z+1/2'
sg_symmetry(129, 227) = '-z+1/2,-x,-y+1/2'
sg_symmetry(130, 227) = 'x+1/2,-z-1/4,-y+1/4'
sg_symmetry(131, 227) = 'z-1/4,x-1/4,-y'
sg_symmetry(132, 227) = '-x-1/4,z,-y-1/4'
sg_symmetry(133, 227) = '-z+1/2,x-1/4,y+1/4'
sg_symmetry(134, 227) = '-x-1/4,-z-1/4,y'
sg_symmetry(135, 227) = 'z-1/4,-x,y-1/4'
sg_symmetry(136, 227) = 'x+1/2,z,y+1/2'
sg_symmetry(137, 227) = '-y+1/2,-z,-x+1/2'
sg_symmetry(138, 227) = '-y,z-3/4,x+1/4'
sg_symmetry(139, 227) = '-z+1/4,-y-3/4,x'
sg_symmetry(140, 227) = 'y+1/4,-z-1/2,x-1/4'
sg_symmetry(141, 227) = 'z+1/2,y-1/2,x'
sg_symmetry(142, 227) = 'y+1/4,z-1/4,-x+1/2'
sg_symmetry(143, 227) = '-z+1/4,y,-x+1/4'
sg_symmetry(144, 227) = 'z,-y-1/4,-x-1/4'
sg_symmetry(145, 227) = 'x+1/2,y+1/2,z'
sg_symmetry(146, 227) = '-y+1/2,x+3/4,z+1/4'
sg_symmetry(147, 227) = '-x+5/4,-y+3/4,z+1/2'
sg_symmetry(148, 227) = 'y+5/4,-x+1/2,z+3/4'
sg_symmetry(149, 227) = 'x+1/2,-y+3/4,-z+1/4'
sg_symmetry(150, 227) = 'y+5/4,x+3/4,-z+1/2'
sg_symmetry(151, 227) = '-x+5/4,y+1/2,-z+3/4'
sg_symmetry(152, 227) = '-y+1/2,-x+1/2,-z'
sg_symmetry(153, 227) = 'z+1/2,x+1/2,y'
sg_symmetry(154, 227) = '-x+1/2,z+3/4,y+1/4'
sg_symmetry(155, 227) = '-z+5/4,-x+3/4,y+1/2'
sg_symmetry(156, 227) = 'x+5/4,-z+1/2,y+3/4'
sg_symmetry(157, 227) = 'z+1/2,-x+3/4,-y+1/4'
sg_symmetry(158, 227) = 'x+5/4,z+3/4,-y+1/2'
sg_symmetry(159, 227) = '-z+5/4,x+1/2,-y+3/4'
sg_symmetry(160, 227) = '-x+1/2,-z+1/2,-y'
sg_symmetry(161, 227) = 'y+1/2,z+1/2,x'
sg_symmetry(162, 227) = 'y+1,-z+5/4,-x+1/4'
sg_symmetry(163, 227) = 'z+3/4,y+5/4,-x+1/2'
sg_symmetry(164, 227) = '-y+3/4,z+1,-x+3/4'
sg_symmetry(165, 227) = '-z+1/2,-y+1,-x+1/2'
sg_symmetry(166, 227) = '-y+3/4,-z+3/4,x'
sg_symmetry(167, 227) = 'z+3/4,-y+1/2,x+1/4'
sg_symmetry(168, 227) = '-z+1,y+3/4,x+3/4'
sg_symmetry(169, 227) = '-x+1/2,-y+1/2,-z'
sg_symmetry(170, 227) = 'y+1/2,-x+1/4,-z-1/4'
sg_symmetry(171, 227) = 'x-1/4,y+1/4,-z-1/2'
sg_symmetry(172, 227) = '-y-1/4,x+1/2,-z-3/4'
sg_symmetry(173, 227) = '-x+1/2,y+1/4,z-1/4'
sg_symmetry(174, 227) = '-y-1/4,-x+1/4,z-1/2'
sg_symmetry(175, 227) = 'x-1/4,-y+1/2,z-3/4'
sg_symmetry(176, 227) = 'y+1/2,x+1/2,z'
sg_symmetry(177, 227) = '-z+1/2,-x+1/2,-y'
sg_symmetry(178, 227) = 'x+1/2,-z+1/4,-y-1/4'
sg_symmetry(179, 227) = 'z-1/4,x+1/4,-y-1/2'
sg_symmetry(180, 227) = '-x-1/4,z+1/2,-y-3/4'
sg_symmetry(181, 227) = '-z+1/2,x+1/4,y-1/4'
sg_symmetry(182, 227) = '-x-1/4,-z+1/4,y-1/2'
sg_symmetry(183, 227) = 'z-1/4,-x+1/2,y-3/4'
sg_symmetry(184, 227) = 'x+1/2,z+1/2,y'
sg_symmetry(185, 227) = '-y+1/2,-z+1/2,-x'
sg_symmetry(186, 227) = '-y,z-1/4,x-1/4'
sg_symmetry(187, 227) = '-z+1/4,-y-1/4,x-1/2'
sg_symmetry(188, 227) = 'y+1/4,-z,x-3/4'
sg_symmetry(189, 227) = 'z+1/2,y,x-1/2'
sg_symmetry(190, 227) = 'y+1/4,z+1/4,-x'
sg_symmetry(191, 227) = '-z+1/4,y+1/2,-x-1/4'
sg_symmetry(192, 227) = 'z,-y+1/4,-x-3/4'
sg_name(228) = 'Fd-3c'
sg_patn(228) = 228
sg_symnum(228) = 192
sg_symmetry(1, 228) = 'x,y,z'
sg_symmetry(2, 228) = '-y+1/2,x+1/4,z+1/4'
sg_symmetry(3, 228) = '-x+1/4,-y+3/4,z+1/2'
sg_symmetry(4, 228) = 'y+3/4,-x+1/2,z+3/4'
sg_symmetry(5, 228) = 'x,-y+1/4,-z+1/4'
sg_symmetry(6, 228) = 'y+1/4,x+1/4,-z+1/2'
sg_symmetry(7, 228) = '-x+1/4,y+1/2,-z+3/4'
sg_symmetry(8, 228) = '-y,-x+1/2,-z'
sg_symmetry(9, 228) = 'z,x,y'
sg_symmetry(10, 228) = '-x+1/2,z+1/4,y+1/4'
sg_symmetry(11, 228) = '-z+1/4,-x+3/4,y+1/2'
sg_symmetry(12, 228) = 'x+3/4,-z+1/2,y+3/4'
sg_symmetry(13, 228) = 'z,-x+1/4,-y+1/4'
sg_symmetry(14, 228) = 'x+1/4,z+1/4,-y+1/2'
sg_symmetry(15, 228) = '-z+1/4,x+1/2,-y+3/4'
sg_symmetry(16, 228) = '-x,-z+1/2,-y'
sg_symmetry(17, 228) = 'y,z,x'
sg_symmetry(18, 228) = 'y+1/2,-z+1/4,-x+3/4'
sg_symmetry(19, 228) = 'z+1/4,y+3/4,-x'
sg_symmetry(20, 228) = '-y+3/4,z+1/2,-x+1/4'
sg_symmetry(21, 228) = '-z+1/2,-y+1/2,-x+1/2'
sg_symmetry(22, 228) = '-y+1/4,-z+1/4,x'
sg_symmetry(23, 228) = 'z+1/4,-y,x+3/4'
sg_symmetry(24, 228) = '-z,y+3/4,x+1/4'
sg_symmetry(25, 228) = '-x,-y,-z'
sg_symmetry(26, 228) = 'y-1/2,-x-1/4,-z-1/4'
sg_symmetry(27, 228) = 'x-1/4,y-3/4,-z-1/2'
sg_symmetry(28, 228) = '-y-3/4,x-1/2,-z-3/4'
sg_symmetry(29, 228) = '-x,y-1/4,z-1/4'
sg_symmetry(30, 228) = '-y-1/4,-x-1/4,z-1/2'
sg_symmetry(31, 228) = 'x-1/4,-y-1/2,z-3/4'
sg_symmetry(32, 228) = 'y,x-1/2,z'
sg_symmetry(33, 228) = '-z,-x,-y'
sg_symmetry(34, 228) = 'x-1/2,-z-1/4,-y-1/4'
sg_symmetry(35, 228) = 'z-1/4,x-3/4,-y-1/2'
sg_symmetry(36, 228) = '-x-3/4,z-1/2,-y-3/4'
sg_symmetry(37, 228) = '-z,x-1/4,y-1/4'
sg_symmetry(38, 228) = '-x-1/4,-z-1/4,y-1/2'
sg_symmetry(39, 228) = 'z-1/4,-x-1/2,y-3/4'
sg_symmetry(40, 228) = 'x,z-1/2,y'
sg_symmetry(41, 228) = '-y,-z,-x'
sg_symmetry(42, 228) = '-y-1/2,z-1/4,x-3/4'
sg_symmetry(43, 228) = '-z-1/4,-y-3/4,x'
sg_symmetry(44, 228) = 'y-3/4,-z-1/2,x-1/4'
sg_symmetry(45, 228) = 'z-1/2,y-1/2,x-1/2'
sg_symmetry(46, 228) = 'y-1/4,z-1/4,-x'
sg_symmetry(47, 228) = '-z-1/4,y,-x-3/4'
sg_symmetry(48, 228) = 'z,-y-3/4,-x-1/4'
sg_symmetry(49, 228) = 'x,y+1/2,z+1/2'
sg_symmetry(50, 228) = '-y+1/2,x+3/4,z+3/4'
sg_symmetry(51, 228) = '-x+1/4,-y+5/4,z+1'
sg_symmetry(52, 228) = 'y+3/4,-x+1,z+5/4'
sg_symmetry(53, 228) = 'x,-y+3/4,-z+3/4'
sg_symmetry(54, 228) = 'y+1/4,x+3/4,-z+1'
sg_symmetry(55, 228) = '-x+1/4,y+1,-z+5/4'
sg_symmetry(56, 228) = '-y,-x+1,-z+1/2'
sg_symmetry(57, 228) = 'z,x+1/2,y+1/2'
sg_symmetry(58, 228) = '-x+1/2,z+3/4,y+3/4'
sg_symmetry(59, 228) = '-z+1/4,-x+5/4,y+1'
sg_symmetry(60, 228) = 'x+3/4,-z+1,y+5/4'
sg_symmetry(61, 228) = 'z,-x+3/4,-y+3/4'
sg_symmetry(62, 228) = 'x+1/4,z+3/4,-y+1'
sg_symmetry(63, 228) = '-z+1/4,x+1,-y+5/4'
sg_symmetry(64, 228) = '-x,-z+1,-y+1/2'
sg_symmetry(65, 228) = 'y,z+1/2,x+1/2'
sg_symmetry(66, 228) = 'y+1/2,-z+3/4,-x+5/4'
sg_symmetry(67, 228) = 'z+1/4,y+5/4,-x+1/2'
sg_symmetry(68, 228) = '-y+3/4,z+1,-x+3/4'
sg_symmetry(69, 228) = '-z+1/2,-y+1,-x+1'
sg_symmetry(70, 228) = '-y+1/4,-z+3/4,x+1/2'
sg_symmetry(71, 228) = 'z+1/4,-y+1/2,x+5/4'
sg_symmetry(72, 228) = '-z,y+5/4,x+3/4'
sg_symmetry(73, 228) = '-x,-y+1/2,-z+1/2'
sg_symmetry(74, 228) = 'y-1/2,-x+1/4,-z+1/4'
sg_symmetry(75, 228) = 'x-1/4,y-1/4,-z'
sg_symmetry(76, 228) = '-y-3/4,x,-z-1/4'
sg_symmetry(77, 228) = '-x,y+1/4,z+1/4'
sg_symmetry(78, 228) = '-y-1/4,-x+1/4,z'
sg_symmetry(79, 228) = 'x-1/4,-y,z-1/4'
sg_symmetry(80, 228) = 'y,x,z+1/2'
sg_symmetry(81, 228) = '-z,-x+1/2,-y+1/2'
sg_symmetry(82, 228) = 'x-1/2,-z+1/4,-y+1/4'
sg_symmetry(83, 228) = 'z-1/4,x-1/4,-y'
sg_symmetry(84, 228) = '-x-3/4,z,-y-1/4'
sg_symmetry(85, 228) = '-z,x+1/4,y+1/4'
sg_symmetry(86, 228) = '-x-1/4,-z+1/4,y'
sg_symmetry(87, 228) = 'z-1/4,-x,y-1/4'
sg_symmetry(88, 228) = 'x,z,y+1/2'
sg_symmetry(89, 228) = '-y,-z+1/2,-x+1/2'
sg_symmetry(90, 228) = '-y-1/2,z+1/4,x-1/4'
sg_symmetry(91, 228) = '-z-1/4,-y-1/4,x+1/2'
sg_symmetry(92, 228) = 'y-3/4,-z,x+1/4'
sg_symmetry(93, 228) = 'z-1/2,y,x'
sg_symmetry(94, 228) = 'y-1/4,z+1/4,-x+1/2'
sg_symmetry(95, 228) = '-z-1/4,y+1/2,-x-1/4'
sg_symmetry(96, 228) = 'z,-y-1/4,-x+1/4'
sg_symmetry(97, 228) = 'x+1/2,y,z+1/2'
sg_symmetry(98, 228) = '-y+1,x+1/4,z+3/4'
sg_symmetry(99, 228) = '-x+3/4,-y+3/4,z+1'
sg_symmetry(100, 228) = 'y+5/4,-x+1/2,z+5/4'
sg_symmetry(101, 228) = 'x+1/2,-y+1/4,-z+3/4'
sg_symmetry(102, 228) = 'y+3/4,x+1/4,-z+1'
sg_symmetry(103, 228) = '-x+3/4,y+1/2,-z+5/4'
sg_symmetry(104, 228) = '-y+1/2,-x+1/2,-z+1/2'
sg_symmetry(105, 228) = 'z+1/2,x,y+1/2'
sg_symmetry(106, 228) = '-x+1,z+1/4,y+3/4'
sg_symmetry(107, 228) = '-z+3/4,-x+3/4,y+1'
sg_symmetry(108, 228) = 'x+5/4,-z+1/2,y+5/4'
sg_symmetry(109, 228) = 'z+1/2,-x+1/4,-y+3/4'
sg_symmetry(110, 228) = 'x+3/4,z+1/4,-y+1'
sg_symmetry(111, 228) = '-z+3/4,x+1/2,-y+5/4'
sg_symmetry(112, 228) = '-x+1/2,-z+1/2,-y+1/2'
sg_symmetry(113, 228) = 'y+1/2,z,x+1/2'
sg_symmetry(114, 228) = 'y+1,-z+1/4,-x+5/4'
sg_symmetry(115, 228) = 'z+3/4,y+3/4,-x+1/2'
sg_symmetry(116, 228) = '-y+5/4,z+1/2,-x+3/4'
sg_symmetry(117, 228) = '-z+1,-y+1/2,-x+1'
sg_symmetry(118, 228) = '-y+3/4,-z+1/4,x+1/2'
sg_symmetry(119, 228) = 'z+3/4,-y,x+5/4'
sg_symmetry(120, 228) = '-z+1/2,y+3/4,x+3/4'
sg_symmetry(121, 228) = '-x+1/2,-y,-z+1/2'
sg_symmetry(122, 228) = 'y,-x-1/4,-z+1/4'
sg_symmetry(123, 228) = 'x+1/4,y-3/4,-z'
sg_symmetry(124, 228) = '-y-1/4,x-1/2,-z-1/4'
sg_symmetry(125, 228) = '-x+1/2,y-1/4,z+1/4'
sg_symmetry(126, 228) = '-y+1/4,-x-1/4,z'
sg_symmetry(127, 228) = 'x+1/4,-y-1/2,z-1/4'
sg_symmetry(128, 228) = 'y+1/2,x-1/2,z+1/2'
sg_symmetry(129, 228) = '-z+1/2,-x,-y+1/2'
sg_symmetry(130, 228) = 'x,-z-1/4,-y+1/4'
sg_symmetry(131, 228) = 'z+1/4,x-3/4,-y'
sg_symmetry(132, 228) = '-x-1/4,z-1/2,-y-1/4'
sg_symmetry(133, 228) = '-z+1/2,x-1/4,y+1/4'
sg_symmetry(134, 228) = '-x+1/4,-z-1/4,y'
sg_symmetry(135, 228) = 'z+1/4,-x-1/2,y-1/4'
sg_symmetry(136, 228) = 'x+1/2,z-1/2,y+1/2'
sg_symmetry(137, 228) = '-y+1/2,-z,-x+1/2'
sg_symmetry(138, 228) = '-y,z-1/4,x-1/4'
sg_symmetry(139, 228) = '-z+1/4,-y-3/4,x+1/2'
sg_symmetry(140, 228) = 'y-1/4,-z-1/2,x+1/4'
sg_symmetry(141, 228) = 'z,y-1/2,x'
sg_symmetry(142, 228) = 'y+1/4,z-1/4,-x+1/2'
sg_symmetry(143, 228) = '-z+1/4,y,-x-1/4'
sg_symmetry(144, 228) = 'z+1/2,-y-3/4,-x+1/4'
sg_symmetry(145, 228) = 'x+1/2,y+1/2,z'
sg_symmetry(146, 228) = '-y+1,x+3/4,z+1/4'
sg_symmetry(147, 228) = '-x+3/4,-y+5/4,z+1/2'
sg_symmetry(148, 228) = 'y+5/4,-x+1,z+3/4'
sg_symmetry(149, 228) = 'x+1/2,-y+3/4,-z+1/4'
sg_symmetry(150, 228) = 'y+3/4,x+3/4,-z+1/2'
sg_symmetry(151, 228) = '-x+3/4,y+1,-z+3/4'
sg_symmetry(152, 228) = '-y+1/2,-x+1,-z'
sg_symmetry(153, 228) = 'z+1/2,x+1/2,y'
sg_symmetry(154, 228) = '-x+1,z+3/4,y+1/4'
sg_symmetry(155, 228) = '-z+3/4,-x+5/4,y+1/2'
sg_symmetry(156, 228) = 'x+5/4,-z+1,y+3/4'
sg_symmetry(157, 228) = 'z+1/2,-x+3/4,-y+1/4'
sg_symmetry(158, 228) = 'x+3/4,z+3/4,-y+1/2'
sg_symmetry(159, 228) = '-z+3/4,x+1,-y+3/4'
sg_symmetry(160, 228) = '-x+1/2,-z+1,-y'
sg_symmetry(161, 228) = 'y+1/2,z+1/2,x'
sg_symmetry(162, 228) = 'y+1,-z+3/4,-x+3/4'
sg_symmetry(163, 228) = 'z+3/4,y+5/4,-x'
sg_symmetry(164, 228) = '-y+5/4,z+1,-x+1/4'
sg_symmetry(165, 228) = '-z+1,-y+1,-x+1/2'
sg_symmetry(166, 228) = '-y+3/4,-z+3/4,x'
sg_symmetry(167, 228) = 'z+3/4,-y+1/2,x+3/4'
sg_symmetry(168, 228) = '-z+1/2,y+5/4,x+1/4'
sg_symmetry(169, 228) = '-x+1/2,-y+1/2,-z'
sg_symmetry(170, 228) = 'y,-x+1/4,-z-1/4'
sg_symmetry(171, 228) = 'x+1/4,y-1/4,-z-1/2'
sg_symmetry(172, 228) = '-y-1/4,x,-z-3/4'
sg_symmetry(173, 228) = '-x+1/2,y+1/4,z-1/4'
sg_symmetry(174, 228) = '-y+1/4,-x+1/4,z-1/2'
sg_symmetry(175, 228) = 'x+1/4,-y,z-3/4'
sg_symmetry(176, 228) = 'y+1/2,x,z'
sg_symmetry(177, 228) = '-z+1/2,-x+1/2,-y'
sg_symmetry(178, 228) = 'x,-z+1/4,-y-1/4'
sg_symmetry(179, 228) = 'z+1/4,x-1/4,-y-1/2'
sg_symmetry(180, 228) = '-x-1/4,z,-y-3/4'
sg_symmetry(181, 228) = '-z+1/2,x+1/4,y-1/4'
sg_symmetry(182, 228) = '-x+1/4,-z+1/4,y-1/2'
sg_symmetry(183, 228) = 'z+1/4,-x,y-3/4'
sg_symmetry(184, 228) = 'x+1/2,z,y'
sg_symmetry(185, 228) = '-y+1/2,-z+1/2,-x'
sg_symmetry(186, 228) = '-y,z+1/4,x-3/4'
sg_symmetry(187, 228) = '-z+1/4,-y-1/4,x'
sg_symmetry(188, 228) = 'y-1/4,-z,x-1/4'
sg_symmetry(189, 228) = 'z,y,x-1/2'
sg_symmetry(190, 228) = 'y+1/4,z+1/4,-x'
sg_symmetry(191, 228) = '-z+1/4,y+1/2,-x-3/4'
sg_symmetry(192, 228) = 'z+1/2,-y-1/4,-x-1/4'
sg_name(229) = 'Im-3m'
sg_patn(229) = 229
sg_symnum(229) = 96
sg_symmetry(1, 229) = 'x,y,z'
sg_symmetry(2, 229) = '-y,x,z'
sg_symmetry(3, 229) = '-x,-y,z'
sg_symmetry(4, 229) = 'y,-x,z'
sg_symmetry(5, 229) = 'x,-y,-z'
sg_symmetry(6, 229) = 'y,x,-z'
sg_symmetry(7, 229) = '-x,y,-z'
sg_symmetry(8, 229) = '-y,-x,-z'
sg_symmetry(9, 229) = 'z,x,y'
sg_symmetry(10, 229) = '-x,z,y'
sg_symmetry(11, 229) = '-z,-x,y'
sg_symmetry(12, 229) = 'x,-z,y'
sg_symmetry(13, 229) = 'z,-x,-y'
sg_symmetry(14, 229) = 'x,z,-y'
sg_symmetry(15, 229) = '-z,x,-y'
sg_symmetry(16, 229) = '-x,-z,-y'
sg_symmetry(17, 229) = 'y,z,x'
sg_symmetry(18, 229) = 'y,-z,-x'
sg_symmetry(19, 229) = 'z,y,-x'
sg_symmetry(20, 229) = '-y,z,-x'
sg_symmetry(21, 229) = '-z,-y,-x'
sg_symmetry(22, 229) = '-y,-z,x'
sg_symmetry(23, 229) = 'z,-y,x'
sg_symmetry(24, 229) = '-z,y,x'
sg_symmetry(25, 229) = '-x,-y,-z'
sg_symmetry(26, 229) = 'y,-x,-z'
sg_symmetry(27, 229) = 'x,y,-z'
sg_symmetry(28, 229) = '-y,x,-z'
sg_symmetry(29, 229) = '-x,y,z'
sg_symmetry(30, 229) = '-y,-x,z'
sg_symmetry(31, 229) = 'x,-y,z'
sg_symmetry(32, 229) = 'y,x,z'
sg_symmetry(33, 229) = '-z,-x,-y'
sg_symmetry(34, 229) = 'x,-z,-y'
sg_symmetry(35, 229) = 'z,x,-y'
sg_symmetry(36, 229) = '-x,z,-y'
sg_symmetry(37, 229) = '-z,x,y'
sg_symmetry(38, 229) = '-x,-z,y'
sg_symmetry(39, 229) = 'z,-x,y'
sg_symmetry(40, 229) = 'x,z,y'
sg_symmetry(41, 229) = '-y,-z,-x'
sg_symmetry(42, 229) = '-y,z,x'
sg_symmetry(43, 229) = '-z,-y,x'
sg_symmetry(44, 229) = 'y,-z,x'
sg_symmetry(45, 229) = 'z,y,x'
sg_symmetry(46, 229) = 'y,z,-x'
sg_symmetry(47, 229) = '-z,y,-x'
sg_symmetry(48, 229) = 'z,-y,-x'
sg_symmetry(49, 229) = 'x+1/2,y+1/2,z+1/2'
sg_symmetry(50, 229) = '-y+1/2,x+1/2,z+1/2'
sg_symmetry(51, 229) = '-x+1/2,-y+1/2,z+1/2'
sg_symmetry(52, 229) = 'y+1/2,-x+1/2,z+1/2'
sg_symmetry(53, 229) = 'x+1/2,-y+1/2,-z+1/2'
sg_symmetry(54, 229) = 'y+1/2,x+1/2,-z+1/2'
sg_symmetry(55, 229) = '-x+1/2,y+1/2,-z+1/2'
sg_symmetry(56, 229) = '-y+1/2,-x+1/2,-z+1/2'
sg_symmetry(57, 229) = 'z+1/2,x+1/2,y+1/2'
sg_symmetry(58, 229) = '-x+1/2,z+1/2,y+1/2'
sg_symmetry(59, 229) = '-z+1/2,-x+1/2,y+1/2'
sg_symmetry(60, 229) = 'x+1/2,-z+1/2,y+1/2'
sg_symmetry(61, 229) = 'z+1/2,-x+1/2,-y+1/2'
sg_symmetry(62, 229) = 'x+1/2,z+1/2,-y+1/2'
sg_symmetry(63, 229) = '-z+1/2,x+1/2,-y+1/2'
sg_symmetry(64, 229) = '-x+1/2,-z+1/2,-y+1/2'
sg_symmetry(65, 229) = 'y+1/2,z+1/2,x+1/2'
sg_symmetry(66, 229) = 'y+1/2,-z+1/2,-x+1/2'
sg_symmetry(67, 229) = 'z+1/2,y+1/2,-x+1/2'
sg_symmetry(68, 229) = '-y+1/2,z+1/2,-x+1/2'
sg_symmetry(69, 229) = '-z+1/2,-y+1/2,-x+1/2'
sg_symmetry(70, 229) = '-y+1/2,-z+1/2,x+1/2'
sg_symmetry(71, 229) = 'z+1/2,-y+1/2,x+1/2'
sg_symmetry(72, 229) = '-z+1/2,y+1/2,x+1/2'
sg_symmetry(73, 229) = '-x+1/2,-y+1/2,-z+1/2'
sg_symmetry(74, 229) = 'y+1/2,-x+1/2,-z+1/2'
sg_symmetry(75, 229) = 'x+1/2,y+1/2,-z+1/2'
sg_symmetry(76, 229) = '-y+1/2,x+1/2,-z+1/2'
sg_symmetry(77, 229) = '-x+1/2,y+1/2,z+1/2'
sg_symmetry(78, 229) = '-y+1/2,-x+1/2,z+1/2'
sg_symmetry(79, 229) = 'x+1/2,-y+1/2,z+1/2'
sg_symmetry(80, 229) = 'y+1/2,x+1/2,z+1/2'
sg_symmetry(81, 229) = '-z+1/2,-x+1/2,-y+1/2'
sg_symmetry(82, 229) = 'x+1/2,-z+1/2,-y+1/2'
sg_symmetry(83, 229) = 'z+1/2,x+1/2,-y+1/2'
sg_symmetry(84, 229) = '-x+1/2,z+1/2,-y+1/2'
sg_symmetry(85, 229) = '-z+1/2,x+1/2,y+1/2'
sg_symmetry(86, 229) = '-x+1/2,-z+1/2,y+1/2'
sg_symmetry(87, 229) = 'z+1/2,-x+1/2,y+1/2'
sg_symmetry(88, 229) = 'x+1/2,z+1/2,y+1/2'
sg_symmetry(89, 229) = '-y+1/2,-z+1/2,-x+1/2'
sg_symmetry(90, 229) = '-y+1/2,z+1/2,x+1/2'
sg_symmetry(91, 229) = '-z+1/2,-y+1/2,x+1/2'
sg_symmetry(92, 229) = 'y+1/2,-z+1/2,x+1/2'
sg_symmetry(93, 229) = 'z+1/2,y+1/2,x+1/2'
sg_symmetry(94, 229) = 'y+1/2,z+1/2,-x+1/2'
sg_symmetry(95, 229) = '-z+1/2,y+1/2,-x+1/2'
sg_symmetry(96, 229) = 'z+1/2,-y+1/2,-x+1/2'
sg_name(230) = 'Ia-3d'
sg_patn(230) = 230
sg_symnum(230) = 96
sg_symmetry(1, 230) = 'x,y,z'
sg_symmetry(2, 230) = '-y+1/4,x+3/4,z+1/4'
sg_symmetry(3, 230) = '-x+1/2,-y,z+1/2'
sg_symmetry(4, 230) = 'y+1/4,-x+1/4,z+3/4'
sg_symmetry(5, 230) = 'x,-y,-z+1/2'
sg_symmetry(6, 230) = 'y+1/4,x+3/4,-z+3/4'
sg_symmetry(7, 230) = '-x+1/2,y,-z'
sg_symmetry(8, 230) = '-y+1/4,-x+1/4,-z+1/4'
sg_symmetry(9, 230) = 'z,x,y'
sg_symmetry(10, 230) = '-x+1/4,z+3/4,y+1/4'
sg_symmetry(11, 230) = '-z+1/2,-x,y+1/2'
sg_symmetry(12, 230) = 'x+1/4,-z+1/4,y+3/4'
sg_symmetry(13, 230) = 'z,-x,-y+1/2'
sg_symmetry(14, 230) = 'x+1/4,z+3/4,-y+3/4'
sg_symmetry(15, 230) = '-z+1/2,x,-y'
sg_symmetry(16, 230) = '-x+1/4,-z+1/4,-y+1/4'
sg_symmetry(17, 230) = 'y,z,x'
sg_symmetry(18, 230) = 'y+1/2,-z+1/2,-x'
sg_symmetry(19, 230) = 'z+3/4,y+1/4,-x+1/4'
sg_symmetry(20, 230) = '-y,z+1/2,-x+1/2'
sg_symmetry(21, 230) = '-z+1/4,-y+1/4,-x+1/4'
sg_symmetry(22, 230) = '-y+1/2,-z,x+1/2'
sg_symmetry(23, 230) = 'z+3/4,-y+3/4,x+1/4'
sg_symmetry(24, 230) = '-z+3/4,y+1/4,x+3/4'
sg_symmetry(25, 230) = '-x,-y,-z'
sg_symmetry(26, 230) = 'y-1/4,-x-3/4,-z-1/4'
sg_symmetry(27, 230) = 'x-1/2,y,-z-1/2'
sg_symmetry(28, 230) = '-y-1/4,x-1/4,-z-3/4'
sg_symmetry(29, 230) = '-x,y,z-1/2'
sg_symmetry(30, 230) = '-y-1/4,-x-3/4,z-3/4'
sg_symmetry(31, 230) = 'x-1/2,-y,z'
sg_symmetry(32, 230) = 'y-1/4,x-1/4,z-1/4'
sg_symmetry(33, 230) = '-z,-x,-y'
sg_symmetry(34, 230) = 'x-1/4,-z-3/4,-y-1/4'
sg_symmetry(35, 230) = 'z-1/2,x,-y-1/2'
sg_symmetry(36, 230) = '-x-1/4,z-1/4,-y-3/4'
sg_symmetry(37, 230) = '-z,x,y-1/2'
sg_symmetry(38, 230) = '-x-1/4,-z-3/4,y-3/4'
sg_symmetry(39, 230) = 'z-1/2,-x,y'
sg_symmetry(40, 230) = 'x-1/4,z-1/4,y-1/4'
sg_symmetry(41, 230) = '-y,-z,-x'
sg_symmetry(42, 230) = '-y-1/2,z-1/2,x'
sg_symmetry(43, 230) = '-z-3/4,-y-1/4,x-1/4'
sg_symmetry(44, 230) = 'y,-z-1/2,x-1/2'
sg_symmetry(45, 230) = 'z-1/4,y-1/4,x-1/4'
sg_symmetry(46, 230) = 'y-1/2,z,-x-1/2'
sg_symmetry(47, 230) = '-z-3/4,y-3/4,-x-1/4'
sg_symmetry(48, 230) = 'z-3/4,-y-1/4,-x-3/4'
sg_symmetry(49, 230) = 'x+1/2,y+1/2,z+1/2'
sg_symmetry(50, 230) = '-y+3/4,x+5/4,z+3/4'
sg_symmetry(51, 230) = '-x+1,-y+1/2,z+1'
sg_symmetry(52, 230) = 'y+3/4,-x+3/4,z+5/4'
sg_symmetry(53, 230) = 'x+1/2,-y+1/2,-z+1'
sg_symmetry(54, 230) = 'y+3/4,x+5/4,-z+5/4'
sg_symmetry(55, 230) = '-x+1,y+1/2,-z+1/2'
sg_symmetry(56, 230) = '-y+3/4,-x+3/4,-z+3/4'
sg_symmetry(57, 230) = 'z+1/2,x+1/2,y+1/2'
sg_symmetry(58, 230) = '-x+3/4,z+5/4,y+3/4'
sg_symmetry(59, 230) = '-z+1,-x+1/2,y+1'
sg_symmetry(60, 230) = 'x+3/4,-z+3/4,y+5/4'
sg_symmetry(61, 230) = 'z+1/2,-x+1/2,-y+1'
sg_symmetry(62, 230) = 'x+3/4,z+5/4,-y+5/4'
sg_symmetry(63, 230) = '-z+1,x+1/2,-y+1/2'
sg_symmetry(64, 230) = '-x+3/4,-z+3/4,-y+3/4'
sg_symmetry(65, 230) = 'y+1/2,z+1/2,x+1/2'
sg_symmetry(66, 230) = 'y+1,-z+1,-x+1/2'
sg_symmetry(67, 230) = 'z+5/4,y+3/4,-x+3/4'
sg_symmetry(68, 230) = '-y+1/2,z+1,-x+1'
sg_symmetry(69, 230) = '-z+3/4,-y+3/4,-x+3/4'
sg_symmetry(70, 230) = '-y+1,-z+1/2,x+1'
sg_symmetry(71, 230) = 'z+5/4,-y+5/4,x+3/4'
sg_symmetry(72, 230) = '-z+5/4,y+3/4,x+5/4'
sg_symmetry(73, 230) = '-x+1/2,-y+1/2,-z+1/2'
sg_symmetry(74, 230) = 'y+1/4,-x-1/4,-z+1/4'
sg_symmetry(75, 230) = 'x,y+1/2,-z'
sg_symmetry(76, 230) = '-y+1/4,x+1/4,-z-1/4'
sg_symmetry(77, 230) = '-x+1/2,y+1/2,z'
sg_symmetry(78, 230) = '-y+1/4,-x-1/4,z-1/4'
sg_symmetry(79, 230) = 'x,-y+1/2,z+1/2'
sg_symmetry(80, 230) = 'y+1/4,x+1/4,z+1/4'
sg_symmetry(81, 230) = '-z+1/2,-x+1/2,-y+1/2'
sg_symmetry(82, 230) = 'x+1/4,-z-1/4,-y+1/4'
sg_symmetry(83, 230) = 'z,x+1/2,-y'
sg_symmetry(84, 230) = '-x+1/4,z+1/4,-y-1/4'
sg_symmetry(85, 230) = '-z+1/2,x+1/2,y'
sg_symmetry(86, 230) = '-x+1/4,-z-1/4,y-1/4'
sg_symmetry(87, 230) = 'z,-x+1/2,y+1/2'
sg_symmetry(88, 230) = 'x+1/4,z+1/4,y+1/4'
sg_symmetry(89, 230) = '-y+1/2,-z+1/2,-x+1/2'
sg_symmetry(90, 230) = '-y,z,x+1/2'
sg_symmetry(91, 230) = '-z-1/4,-y+1/4,x+1/4'
sg_symmetry(92, 230) = 'y+1/2,-z,x'
sg_symmetry(93, 230) = 'z+1/4,y+1/4,x+1/4'
sg_symmetry(94, 230) = 'y,z+1/2,-x'
sg_symmetry(95, 230) = '-z-1/4,y-1/4,-x+1/4'
sg_symmetry(96, 230) = 'z-1/4,-y+1/4,-x-1/4'
!
state=1
GOTO 1000
!
! error handling (for the future mayhaps)
800 CONTINUE
RETURN
!
! finish.
1000 CONTINUE
!
END SUBROUTINE SG_INIT
!
!
!
!
END MODULE spacegroups