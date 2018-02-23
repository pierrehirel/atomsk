 MODULE in_stl
!
!
!**********************************************************************************
!*  IN_STL                                                                        *
!**********************************************************************************
!* This module reads a file in the STL file format. STL stands for                *
!* "STereoLithography", or "Standard Triangle Language", or                       *
!* "Standard Tessellation Language". It is a standard file format widely used     *
!* in 3-D printing. The STL file format is documented on Wikipedia:               *
!*    https://en.wikipedia.org/wiki/STL_(file_format)                             *
!* or, alternatively, in the following Web sites:                                 *
!*    https://www.3dsystems.com/quickparts/learning-center/what-is-stl-file       *
!*    https://all3dp.com/what-is-stl-file-format-extension-3d-printing/           *
!* NOTE: this format is not a file format of atom positions. For this reason,     *
!* this source file is placed in the folder "include" of Atomsk, and not          *
!* in the folder "input".                                                         *
!**********************************************************************************
!* (C) July 2017 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 22 Feb. 2018                                     *
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
USE files
USE messages
!
IMPLICIT NONE
!
!
CONTAINS
!
SUBROUTINE READ_STL(stlfile,triangles)
!
!Declarations
CHARACTER(LEN=*),INTENT(IN):: stlfile
CHARACTER(LEN=1):: c
CHARACTER(LEN=80):: header
CHARACTER(LEN=128):: line, msg
LOGICAL:: formatted  !is the file formatted?
INTEGER(KIND=2):: ABC  !Attribute Byte Count (binary format)
INTEGER(KIND=4):: N
INTEGER:: i, j, k
INTEGER:: Ntriangles  !number of triangles declared in STL file
REAL(KIND=4):: tempreal
REAL(dp),DIMENSION(:,:),ALLOCATABLE,INTENT(OUT):: triangles !normal vector, and positions of vertices
!
!Initializations
formatted = .FALSE.
Ntriangles = 0
!
!
msg = 'entering READ_STL'
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
!Detect if file exists
CALL CHECKFILE(stlfile,'read')
!
!Detect if file is binary or not
OPEN(UNIT=32,FILE=stlfile,STATUS="OLD",FORM="FORMATTED",ERR=800)
!Attempt to read the keyword "solid"
READ(32,'(a128)',END=800,ERR=800) line
line = ADJUSTL(line)
IF( line(1:5) == "solid" ) THEN
  !99% chances that the file is formatted
  formatted = .TRUE.
ENDIF
!Try to detect non-ASCII characters to check if file is binary
!READ(32,'(a128)',END=800,ERR=800) line
j=0
k=1  !counter for lines. Reading the first 5 lines should be sufficient
DO WHILE( j==0 .AND. formatted .AND. k<=5 )
  READ(32,'(a1)',IOSTAT=j) c
  formatted = formatted .AND. ( IACHAR(c)<=127 )
  k=k+1
ENDDO
!
!
!
200 CONTINUE
CLOSE(32)
WRITE(msg,*) "STL File is formatted: ", formatted
CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
!
IF( .NOT.formatted ) THEN
  !File is in binary format
  !File format specification is as follows:
  !     header (80 characters)
  !     Ntriangles (4-bytes unsigned integer)
  !then a loop over all triangles:
  !     triangle information (twelve 32-bits floating point numbers)
  OPEN(UNIT=32,FILE=stlfile,STATUS="OLD",FORM="UNFORMATTED",ACCESS="STREAM",ERR=800)
  REWIND(32)
  !
  !Read header (not used here)
  READ(32,ERR=800) header
  !
  !Read number of triangles
  READ(32,ERR=800) N
  Ntriangles = N
  WRITE(msg,*) "Ntriangles = ", Ntriangles
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !
  !Allocate the "triangles" array
  ALLOCATE( triangles(Ntriangles,12) )
  triangles(:,:) = 0.d0
  !
  DO i=1,Ntriangles
    !For each triangle, read 12 real numbers
    DO j=1,12
      READ(32,ERR=800) tempreal
      triangles(i,j) = tempreal
    ENDDO
    !Read attribute byte count
    READ(32,ERR=800) ABC
  ENDDO
  !
  CLOSE(32)
  !
  !
ELSE
  !File is formatted/ASCII
  !File format specification is as follows:
  !     solid <name>
  !       facet normal nx ny nz
  !         outer loop
  !           vertex v1x v1y v1z
  !           vertex v2x v2y v2z
  !           vertex v3x v3y v3z
  !         endloop
  !       endfacet
  !       facet normal nx ny nz
  !         ...
  !       endfacet
  !     endsolid <name>
  OPEN(UNIT=32,FILE=stlfile,STATUS="OLD",FORM="FORMATTED",ERR=800)
  !
  READ(32,'(a128)',END=800,ERR=800) line
  line = ADJUSTL(line)
  IF( line(1:5) .NE. "solid" ) THEN
    PRINT*, "WARNING: mandatory keyword 'solid' is missing."
  ENDIF
  !
  !Parse the file a first time to count triangles
  Ntriangles = 0
  DO
    READ(32,'(a128)',END=250,ERR=250) line
    line = ADJUSTL(line)
    IF( line(1:12) == "facet normal" ) THEN
      Ntriangles = Ntriangles + 1
    ENDIF
  ENDDO
  250 CONTINUE
  WRITE(msg,*) "Ntriangles = ", Ntriangles
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  !
  !Allocate the "triangles" array
  ALLOCATE( triangles(Ntriangles,12) )
  triangles(:,:) = 0.d0
  !
  !Go back to beginning of file
  REWIND(32)
  !
  !Re-read file and store data into array "triangles"
  line = ""
  Ntriangles = 0
  DO WHILE ( line(1:8) .NE. "endsolid" )
    READ(32,'(a128)',END=800,ERR=800) line
    line = ADJUSTL(line)
    IF( line(1:12) == "facet normal" ) THEN
      !Read normal to triangle facet
      Ntriangles = Ntriangles + 1
      READ(line(13:),*,END=800,ERR=800) triangles(Ntriangles,1:3)
    ELSEIF( line(1:10) == "outer loop" ) THEN
      !Read the vertex positions
      DO i=1,3
        j = 3*i + 1
        READ(32,'(a128)',END=800,ERR=800) line
        line = ADJUSTL(line)
        READ(line(7:),*,END=800,ERR=800) triangles(Ntriangles,j:j+2)
      ENDDO
    ENDIF
  ENDDO
  !
  CLOSE(32)
  !
ENDIF
!
!
!
500 CONTINUE
GOTO 1000
!
!
!
800 CONTINUE
nerr=nerr+1
!
!
!
1000 CONTINUE
!PRINT*, "STL file was read successfully"
!
IF( verbosity==4 ) THEN
  msg = 'leaving READ_STL'
  CALL ATOMSK_MSG(999,(/msg/),(/0.d0/))
  OPEN(UNIT=33,FILE="atomsk.stli",STATUS="UNKNOWN",FORM="FORMATTED")
  WRITE(33,*) "# Triangle information (first 20 entries) read from STL file: "//TRIM(ADJUSTL(stlfile))
  IF( Ntriangles>0 ) THEN
    DO i=1,MIN(20,SIZE(triangles,1))
      WRITE(33,'(3f6.2,a3,9f9.3)') triangles(i,1:3), " | ", triangles(i,4:)
    ENDDO
  ENDIF
  CLOSE(33)
ENDIF
!
!
END SUBROUTINE READ_STL
!
END MODULE in_stl