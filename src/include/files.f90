MODULE files
!
!**********************************************************************************
!*  FILES                                                                         *
!**********************************************************************************
!* This module contains subroutines performing special                            *
!* transformations on files for ATOMSK.                                           *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 27 Aug. 2024                                     *
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
!* List of subroutines in this module:                                            *
!* FREEUNIT            returns an unused I/O UNIT number                          *
!* SET_OUTPUT          sets one or all output to TRUE or FALSE                    *
!* NAME_OUTFILE        names an output file based on the name of an input file    *
!**********************************************************************************
!
!
USE comv
USE constants
!
IMPLICIT NONE
!
!
CONTAINS
!
!
!********************************************************
!  FREEUNIT
!  This function returns a number (between 11 and 99)
!  corresponding to a I/O UNIT that is not used.
!  This number can then safely be used in an OPEN
!  statement, e.g.:
!    fu = FREEUNIT()
!    OPEN(UNIT=fu,FORM="FORMATTED",STATUS="UNKNOWN")
!********************************************************
FUNCTION FREEUNIT(imin,imax) RESULT(iunit)
!
IMPLICIT NONE
INTEGER:: i, iunit, istart, iend
INTEGER,OPTIONAL:: imin
INTEGER,OPTIONAL:: imax
LOGICAL:: isopen
!
iunit = 0
IF( PRESENT(imin) ) THEN
  istart = imin
ELSE
  istart = 11
ENDIF
IF( PRESENT(imax) ) THEN
  iend = imax
ELSE
  iend = 99
ENDIF
!
i=istart
DO WHILE( i<iend .AND. iunit==0 )
  !
  i=i+1
  INQUIRE (UNIT=i,OPENED=isopen)
  !
  IF ( .NOT. isopen ) THEN
    iunit = i
  ENDIF
ENDDO
!
END FUNCTION FREEUNIT
!
!
!********************************************************
! SET_OUTPUT
! This subroutine adds or removes output format(s)
! from the array "outfileformats".
!********************************************************
SUBROUTINE SET_OUTPUT(outfileformats,extension,l_value)
!
IMPLICIT NONE
CHARACTER(LEN=5):: extension !file extension to add to array
CHARACTER(LEN=5),DIMENSION(:),ALLOCATABLE:: outfileformats, tempout
LOGICAL:: l_value  !TRUE if extension must be added,
                   !FALSE if it must be removed
INTEGER:: i,j
INTEGER:: sizeini
!
IF(ALLOCATED(tempout)) DEALLOCATE(tempout)
sizeini = 0
IF(ALLOCATED(outfileformats)) THEN
  sizeini = SIZE(outfileformats)
ENDIF
!
IF(l_value) THEN
  !the file format must be added to the array
  IF(ALLOCATED(outfileformats)) THEN
    ALLOCATE(tempout(sizeini+1))
    DO i=1,SIZE(outfileformats)
      tempout(i) = outfileformats(i)
    ENDDO
    tempout(SIZE(outfileformats)+1) = extension
    DEALLOCATE(outfileformats)
    ALLOCATE(outfileformats(SIZE(tempout)))
    outfileformats = tempout
    DEALLOCATE(tempout)
  ELSE
    ALLOCATE(outfileformats(1))
    outfileformats(1) = extension
  ENDIF
  !
ELSE
  !one or all file format(s) must be removed from the array
  IF(extension=='all  ' ) THEN
    IF(ALLOCATED(outfileformats)) DEALLOCATE(outfileformats)
    RETURN
  ELSE
    IF(ALLOCATED(outfileformats) .AND. SIZE(outfileformats)>0) THEN
      !First, check if the file format actually exists in array
      j=0
      DO i=1,SIZE(outfileformats)
        IF( outfileformats(i)==extension ) j=1
      ENDDO
      !
      IF(j>0) THEN
        IF(sizeini==1) THEN
          !the extension exists in the array, and is the only entry: wipe out array
          DEALLOCATE(outfileformats)
          RETURN
        ELSE
          !the extension exists in the array => remove it
          ALLOCATE(tempout(SIZE(outfileformats)-1))
          j=0
          DO i=1,SIZE(outfileformats)
            IF( outfileformats(i).NE.extension ) THEN
              j=j+1
              tempout(j) = outfileformats(i)
            ENDIF
          ENDDO
          DEALLOCATE(outfileformats)
          ALLOCATE(outfileformats(SIZE(tempout)))
          outfileformats = tempout
          DEALLOCATE(tempout)
        ENDIF
        !
      ELSE
        !this file extension was not in the array => quit
        RETURN
      ENDIF
      !
    ELSE
      !array "outfileformats" is unallocated or has zero size, nothing to remove
      IF(ALLOCATED(outfileformats)) DEALLOCATE(outfileformats)
      RETURN
    ENDIF
  ENDIF
  !
ENDIF
!
!
END SUBROUTINE SET_OUTPUT
!
!
!********************************************************
! NAME_OUTFILE
! This subroutine finds a name for an output file
! based on the name of an input file and the extension
! for the output file.
! For instance if an input file 'test' is to be converted
! to .def, the output file will be 'test.def'.
! But if the name of the input file is 'test.abc',
! the name of the output file must be 'test.def',
! and not 'test.abc.def'.
!********************************************************
SUBROUTINE NAME_OUTFILE(inputfile,outputfile,outfileformat)
!
IMPLICIT NONE
CHARACTER(LEN=5),INTENT(IN):: outfileformat
CHARACTER(LEN=*),INTENT(IN):: inputfile
CHARACTER(LEN=*),INTENT(OUT):: outputfile
CHARACTER(LEN=5):: ext  !extension of input file
CHARACTER(LEN=4096):: test, test2
INTEGER:: strlength, strlength2
!
IF( ofu==6 ) THEN
  outputfile=""
  !
ELSE
  test = TRIM(ADJUSTL(inputfile))
  !Find where the extension of inputfile starts
  !(we assume here that it has to be after the last path separator, if any)
  strlength = SCAN(test,'.',BACK=.TRUE.)
  strlength2 = SCAN(test,pathsep,BACK=.TRUE.)
  IF(strlength2>=LEN_TRIM(test)) strlength2=0
  IF(strlength>strlength2) THEN
    !Get extension of input file
    ext = TRIM(ADJUSTL(test(strlength+1:)))
    !Verify that ext is one of the known extensions
    IF( ANY( flist(:,1)==ext ) ) THEN
      test = test(1:strlength-1)
    ENDIF
  ENDIF
  !
  !suffix must not have a dot
  test2 = TRIM(ADJUSTL(outfileformat))
  IF(test2(1:1)==".") test2 = test2(2:)
  !
  !Set the name of outputfile
  outputfile = TRIM(ADJUSTL(test))//'.'//TRIM(ADJUSTL(test2))
ENDIF
!
END SUBROUTINE NAME_OUTFILE
!
!
!********************************************************
! FILE_SIZE
! This function determines the size of a file (in bytes)
! and returns a string in human-friendly format,
! e.g. "40.2k" or "3.6G". If file does not exist,
! a blank string is returned.
! Note that only the unit multiple (k, M, G, T) is written,
! as the symbol for "bytes" depends on the language
! (e.g. kilo-byte is "kB" in English, and "ko" in French)
!********************************************************
FUNCTION FILE_SIZE(filename) RESULT(filesize)
!
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: filename
CHARACTER(LEN=16):: filesize, temp
LOGICAL:: fileexists
INTEGER:: sizebytes
REAL(dp):: tempreal
!
filesize=""
!
INQUIRE(FILE=filename,EXIST=fileexists,SIZE=sizebytes)
!
IF( fileexists ) THEN
  !File exists, convert size from bytes to
  tempreal = sizebytes
  IF( tempreal>1024d0**4 ) THEN
    WRITE(temp,'(f12.1)') tempreal/1024d0**4
    filesize = TRIM(ADJUSTL(temp))//" T"
  ELSEIF( tempreal>1024d0**3 ) THEN
    WRITE(temp,'(f12.1)') tempreal/1024d0**3
    filesize = TRIM(ADJUSTL(temp))//" G"
  ELSEIF( tempreal>1024d0**2 ) THEN
    WRITE(temp,'(f12.1)') tempreal/1024d0**2
    filesize = TRIM(ADJUSTL(temp))//" M"
  ELSEIF( tempreal>1024 ) THEN
    WRITE(temp,'(f12.1)') tempreal/1024
    filesize = TRIM(ADJUSTL(temp))//" k"
  ELSEIF( sizebytes>0 ) THEN
    WRITE(temp,*) sizebytes
    filesize = TRIM(ADJUSTL(temp))
  ENDIF
ENDIF
!
END FUNCTION FILE_SIZE
!
!
!
END MODULE files
