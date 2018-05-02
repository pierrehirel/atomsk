MODULE messages_EN
!
!**********************************************************************************
!*  MESSAGES_EN                                                                   *
!**********************************************************************************
!* This module contains the ENGLISH (default)                                     *
!* version of the messages displayed by the Atomsk program.                       *
!**********************************************************************************
!* (C) June 2011 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 19 April 2018                                    *
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
!* DISPLAY_LICENSE_EN  displays the license of the program                        *
!* DISPLAY_HELP_EN     displays the help of the program                           *
!* ATOMSK_MSG_EN       all messages used by Atomsk                                *
!* DATE_MSG_EN         displays a nice message according to the date              *
!**********************************************************************************
!
!
USE atoms
USE comv
USE constants
USE functions
USE subroutines
USE display_messages
!
IMPLICIT NONE
!
!
CONTAINS
!
!
!********************************************************
! DISPLAY_LICENSE
! This subroutine displays the license message
! of Atomsk
!********************************************************
SUBROUTINE DISPLAY_LICENSE_EN()
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
!
msg = ""
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "Atomsk - A tool for manipulating and converting atomic data files."
CALL DISPLAY_MSG(verbosity,msg,logfile)
!
CALL DISPLAY_COPYRIGHT()
!
msg = ""
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "This program is free software: you can redistribute it and/or modify"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "it under the terms of the GNU General Public License as published by"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "the Free Software Foundation, either version 3 of the License, or"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "(at your option) any later version."
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = ""
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "This program is distributed in the hope that it will be useful,"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "but WITHOUT ANY WARRANTY; without even the implied warranty of"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "GNU General Public License for more details."
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = ""
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "You should have received a copy of the GNU General Public License"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "along with this program.  If not, see <http://www.gnu.org/licenses/>."
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = ""
CALL DISPLAY_MSG(verbosity,msg,logfile)
!
END SUBROUTINE DISPLAY_LICENSE_EN
!
!
!********************************************************
! DISPLAY_HELP
! This subroutine displays all or part of
! the help for Atomsk
!********************************************************
SUBROUTINE DISPLAY_HELP_EN(helpsection)
!
IMPLICIT NONE
CHARACTER(LEN=16):: helpsection
!
IF(TRIM(helpsection)=="") helpsection="general"
!
!
IF(helpsection=="general" .OR. helpsection=="") THEN
WRITE(*,*) ">>> USAGE:"
WRITE(*,*) "       atomsk [--<mode>] <inputfile> [<outputfile>] [<formats>] [options...]"
WRITE(*,*) "    where [] are optional parameters, and <> must be replaced"
WRITE(*,*) "    by specific values or strings."
WRITE(*,*) ""
WRITE(*,*) ">>> EXAMPLE:"
WRITE(*,*) "..> Convert 'file.xsf' to CFG format:"
WRITE(*,*) "       atomsk file.xsf cfg"
WRITE(*,*) ""
WRITE(*,*) ">>> FOR MORE HELP:"
WRITE(*,*) "..> On modes:   atomsk --help modes"
WRITE(*,*) "..> On options: atomsk --help options"
WRITE(*,*) "..> On formats: atomsk --help formats"
ENDIF
!
IF(helpsection=="modes") THEN
  WRITE(*,*) ">>> MODES:"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="interactive") THEN
  WRITE(*,*) "..> Interactive mode:"
  WRITE(*,*) "          atomsk"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="normal") THEN
  WRITE(*,*) "..> Normal mode:"
  WRITE(*,*) "          atomsk <inputfile> [<outputfile>] [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="list") THEN
  WRITE(*,*) "..> List mode:"
  WRITE(*,*) "          atomsk --list <listfile> [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="ai1") THEN
  WRITE(*,*) "..> All-in-one mode:"
  WRITE(*,*) "          atomsk --all-in-one <listfile> <outputfile> [options] [<formats>]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="1ia") THEN
  WRITE(*,*) "..> One-in-all mode:"
  WRITE(*,*) "          atomsk --one-in-all <inputfile> [<outputfile>] [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="create") THEN
  WRITE(*,*) "..> Create mode:"
  WRITE(*,*) "          atomsk --create <structure> <a0> <species> <outputfile> [orient hkl hkl hkl] [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="ddplot") THEN
  WRITE(*,*) "..> DDplot mode:"
  WRITE(*,*) "          atomsk --ddplot <file1> <file2> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="merge") THEN
  WRITE(*,*) "..> Merge mode:"
  WRITE(*,*) "          atomsk --merge [<x|y|z>] <Nfiles> <file1>...<fileN> <outputfile> [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="unwrap") THEN
  WRITE(*,*) "..> Unwrap mode:"
  WRITE(*,*) "          atomsk --unwrap <reference> <system> [<outputfile>] [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="density") THEN
  WRITE(*,*) "..> Density mode:"
  WRITE(*,*) "          atomsk --density <file> <property> <dimension> <axis> <sigma> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="difference") THEN
  WRITE(*,*) "..> List mode:"
  WRITE(*,*) "          atomsk --difference <file1> <file2> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="edm") THEN
  WRITE(*,*) "..> Mode electric dipole moments:"
  WRITE(*,*) "          atomsk --edm <system> <Pspecies> <NNN> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="PE") THEN
  WRITE(*,*) "..> Mode electronic polarization:"
  WRITE(*,*) "          atomsk --electronic-polarization <system> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="rdf") THEN
  WRITE(*,*) "..> Mode radial distribution function:"
  WRITE(*,*) "          atomsk --rdf <listfile> <Rmax> <dR> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="polycrystal") THEN
  WRITE(*,*) "..> Mode polycrystal:"
  WRITE(*,*) "          atomsk --polycrystal <unitcell> <parameters> <outputfile> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="average") THEN
  WRITE(*,*) "..> Mode average:"
  WRITE(*,*) "          atomsk --average <listfile> <outputfile> [options] [<formats>]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="nye") THEN
  WRITE(*,*) "..> Mode Nye:"
  WRITE(*,*) "          atomsk --nye <reference> <dislocation> <outputfile> [options] [<formats>]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="interpolate") THEN
  WRITE(*,*) "..> Mode Interpolate:"
  WRITE(*,*) "          atomsk --interpolate <file1> <file2> <N> [options] [<formats>]"
ENDIF
!
IF(helpsection=="options") THEN
  WRITE(*,*) ">>> OPTIONS (distances=Angströms, angles=degrees):"
  WRITE(*,*) ""
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-add-atom" .OR. helpsection=="-add-atoms" .OR. &
  &helpsection=="-addatom" .OR. helpsection=="-addatoms" ) THEN
  WRITE(*,*) "..> Add new atoms to the system:"
  WRITE(*,*) "          -add-atom <species> at <x> <y> <z>"
  WRITE(*,*) "          -add-atom <species> relative <index> <x> <y> <z>"
  WRITE(*,*) "          -add-atom <species> near <index>"
  WRITE(*,*) "          -add-atom <species> random <N>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-add-shells".OR. helpsection=="-as" &
  & .OR. helpsection=="-create-shells".OR. helpsection=="-cs") THEN
  WRITE(*,*) "..> Create shells for some or all atoms:"
  WRITE(*,*) "          -add-shells <all|species>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-alignx") THEN
  WRITE(*,*) "..> Align the first supercell vector with the X axis:"
  WRITE(*,*) "          -alignx"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-bind-shells" .OR. helpsection=="-bs") THEN
  WRITE(*,*) "..> Re-assign ionic shells to their respective cores:"
  WRITE(*,*) "          -bind-shells"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-center") THEN
  WRITE(*,*) "..> Place an atom or the center of mass of the system at the center of the box:"
  WRITE(*,*) "          -center <index>"
  WRITE(*,*) "          -center com"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-crack") THEN
  WRITE(*,*) "..> Insert a crack in the system:"
  WRITE(*,*) "          -crack <I|II|III> <stress|strain> <K> <pos1> <pos2> "//&
           &            "<crackline> <crackplane> <μ> <ν>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-cut") THEN
  WRITE(*,*) "..> Cut part of the system:"
  WRITE(*,*) "          -cut <above|below> <cutdistance> <x|y|z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-deform" .OR. helpsection=="-def") THEN
  WRITE(*,*) "..> Apply uniaxial stress or strain:"
  WRITE(*,*) "          -def <x|y|z> <strain> <Poissons ratio>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-dislocation" .OR. helpsection=="-disloc") THEN
  WRITE(*,*) "..> Insert a dislocation in the system:"
  WRITE(*,*) "          -disloc <pos1> <pos2> screw <x|y|z> <x|y|z> <b>"
  WRITE(*,*) "          -disloc <pos1> <pos2> <edge|edge2> <x|y|z> <x|y|z> <b> <ν>"
  WRITE(*,*) "          -disloc <pos1> <pos2> mixed <x|y|z> <x|y|z> <b1> <b2> <b3>"
  WRITE(*,*) "          -disloc loop <x> <y> <z> <x|y|z> <radius> <bx> <by> <bz> <nu>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-disturb" ) THEN
  WRITE(*,*) "..> Apply a random displacement to atom positions:"
  WRITE(*,*) "          -disturb <dmax>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-duplicate" .OR. helpsection=="-dup") THEN
  WRITE(*,*) "..> Duplicate the system in the 3 directions of space:"
  WRITE(*,*) "          -duplicate <Nx> <Ny> <Nz>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-fix") THEN
  WRITE(*,*) "..> Fix some atoms:"
  WRITE(*,*) "          -fix <x|y|z>"
  WRITE(*,*) "          -fix <x|y|z> <above|below> <distance> <x|y|z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-fractional" .OR. helpsection=="-frac") THEN
  WRITE(*,*) "..> Convert coordinates to fractional:"
  WRITE(*,*) "          -frac"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-mirror") THEN
  WRITE(*,*) "..> Apply a mirror plane:"
  WRITE(*,*) "          -mirror <d> <normal>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-options") THEN
  WRITE(*,*) "..> Apply a list of options read from a file:"
  WRITE(*,*) "          -options <file>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-orient") THEN
  WRITE(*,*) "..> Change the crystallographic orientation of the system:"
  WRITE(*,*) "          -orient <Hx> <Hy> <Hz> <H'x> <H'y> <H'z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-prop") THEN
  WRITE(*,*) "..> Read system properties:"
  WRITE(*,*) "          -prop <file>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-rebox") THEN
  WRITE(*,*) "..> (Re-)calculate vectors of the bounding box:"
  WRITE(*,*) "          -rebox"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-remove-atom" .OR. helpsection=="-rmatom") THEN
  WRITE(*,*) "..> Remove one atom given its index, or all atoms of a species:"
  WRITE(*,*) "          -rmatom <index|species>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-remove-doubles" .OR. helpsection=="-rmd") THEN
  WRITE(*,*) "..> Remove duplicate atoms:"
  WRITE(*,*) "          -rmd <distance>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-remove-property" .OR. helpsection=="-rmprop") THEN
  WRITE(*,*) "..> Remove one or all auxiliary properties:"
  WRITE(*,*) "          -rmprop <property>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-remove-shell" .OR. helpsection=="-rmshell" .OR. &
  &helpsection=="-remove-shells" .OR. helpsection=="-rmshells" ) THEN
  WRITE(*,*) "..> Remove shells from atoms of given species, or on all atoms :"
  WRITE(*,*) "          -rmshells <species|all>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-roll" .OR. helpsection=="-bend") THEN
  WRITE(*,*) "..> Roll the system around an axis:"
  WRITE(*,*) "          -roll <x|y|z> <angle> <x|y|z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-rotate" .OR. helpsection=="-rot") THEN
  WRITE(*,*) "..> Rotate the system around an axis:"
  WRITE(*,*) "          -rotate [com] <x|y|z> <angle>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-select") THEN
  WRITE(*,*) "..> Select atoms according to a criterion:"
  WRITE(*,*) "          -select all"
  WRITE(*,*) "          -select invert"
  WRITE(*,*) "          -select <species>"
  WRITE(*,*) "          -select <index>"
  WRITE(*,*) "          -select <above|below> <d> <normal>"
  WRITE(*,*) "          -select <in|out> <box|sphere|cylinder|torus> [<axis>] <x1> <y1> <z1> <x2> [<y2> <z2>]]"
  WRITE(*,*) "          -select prop <prop> <value>"
  WRITE(*,*) "          -select random <N> <species>"
  WRITE(*,*) "          -select <NNN> <species> neighbors <index>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-separate") THEN
  WRITE(*,*) "..> Separate atoms that are too close:"
  WRITE(*,*) "          -separate <distance> <shift>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-shear") THEN
  WRITE(*,*) "..> Apply simple shear strain to the system:"
  WRITE(*,*) "          -shear <x|y|z> <amplitude> <x|y|z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-shift") THEN
  WRITE(*,*) "..> Shift part of the system:"
  WRITE(*,*) "          -shift <tauX> <tauY> <tauZ>"
  WRITE(*,*) "          -shift <above|below> <distance> <x|y|z> <tauX> <tauY> <tauZ>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-sort") THEN
  WRITE(*,*) "..> Sort atoms:"
  WRITE(*,*) "          -sort <s|x|y|z> <up|down|pack>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-spacegroup") THEN
  WRITE(*,*) "..> Apply symmetry operations of a given space group:"
  WRITE(*,*) "          -spacegroup <group>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-stress") THEN
  WRITE(*,*) "..> Apply stress:"
  WRITE(*,*) "          -stress <xx|yy|zz|xy|xz|yz|P> <value(GPa)>"
  WRITE(*,*) "          -stress <file>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-substitute" .OR. helpsection=="-sub") THEN
  WRITE(*,*) "..> Substitute atoms of species sp1 by species sp2:"
  WRITE(*,*) "          -sub <sp1> <sp2>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-swap") THEN
  WRITE(*,*) "..> Swap two atoms or Cartesian axis:"
  WRITE(*,*) "          -swap <id1> <id2>"
  WRITE(*,*) "          -swap <x|y|z> <x|y|z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-torsion") THEN
  WRITE(*,*) "..> Apply torsion around an axis:"
  WRITE(*,*) "          -torsion <x|y|z> <angle>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-unit" .OR. helpsection=="-u") THEN
  WRITE(*,*) "..> Convert coordinates to another unit:"
  WRITE(*,*) "          -u <unit1> <unit2>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-unskew") THEN
  WRITE(*,*) "..> Reduce the skew of the supercell:"
  WRITE(*,*) "          -unskew"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-velocity") THEN
  WRITE(*,*) "..> Assign velocities to atoms according to a Maxwell-Boltzmann distribution:"
  WRITE(*,*) "          -velocity <Temperature>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-wrap") THEN
  WRITE(*,*) "..> Wrap atoms back into the box:"
  WRITE(*,*) "          -wrap"
ENDIF
!
! -- add other options in alphabetical order --
!
!
IF(helpsection=="formats") THEN
WRITE(*,*) ">>> FORMATS:"
WRITE(*,*) "    Formats must be specified as written in the first column."
WRITE(*,*) "    Atomsk can convert from any 'yes' in the INPUT column"
WRITE(*,*) "    to any 'yes' in the OUTPUT column:"
WRITE(*,*) "                            |  INPUT | OUTPUT"
WRITE(*,*) "    ------------------------+--------+--------"
WRITE(*,*) "    atsk (Atomsk format)    |   yes  |  yes"
WRITE(*,*) "    bop (Bond-Order format) |   yes  |  yes"
WRITE(*,*) "    cel (Dr. Probe/EMS)     |   yes  |  yes"
WRITE(*,*) "    coo (COORAT/MBPP)       |   yes  |  yes"
WRITE(*,*) "    cfg (Atomeye)           |   yes  |  yes"
WRITE(*,*) "    cif (Cryst.Info.File)   |   yes  |  yes"
WRITE(*,*) "    dd  (ddplot)            |    no  | yes (1)"
WRITE(*,*) "    dlp (DL_POLY CONFIG)    |   yes  |  yes"
WRITE(*,*) "    gin (GULP input)        |   yes  |  yes"
WRITE(*,*) "    imd (IMD input)         |   yes  |  yes"
WRITE(*,*) "    jems (JEMS input)       |   yes  |  yes"
WRITE(*,*) "    lmc (LAMMPS output)     |   yes  |   no"
WRITE(*,*) "    lmp (LAMMPS data)       |   yes  |  yes"
WRITE(*,*) "    mol (MOLDY format)      |   yes  |  yes"
WRITE(*,*) "    pdb (Protein Data Bank) |   yes  |  yes"
WRITE(*,*) "    pos (POSCAR/VASP)       |   yes  |  yes"
WRITE(*,*) "    pw (Quantum Espresso)   |   yes  |  yes"
WRITE(*,*) "    pwout (QE output file)  |  yes(2)|   no"
WRITE(*,*) "    str (PDFFIT)            |   yes  |  yes"
WRITE(*,*) "    vesta (VESTA file)      |   yes  |  yes"
WRITE(*,*) "    xmd (XMD file)          |   yes  |  yes"
WRITE(*,*) "    xsf (XCrySDen)          |   yes  |  yes"
WRITE(*,*) "    xv (SIESTA format)      |   yes  |  yes"
WRITE(*,*) "    xyz/exyz/sxyz           |   yes  |  yes"
WRITE(*,*) "        (1) Mode ddplot only."
WRITE(*,*) "        (2) Mode one-in-all only."
ENDIF
!
WRITE(*,*) ""
WRITE(*,*) ">>> Look at the /doc folder provided with the program"
WRITE(*,*) "    or go to: http://atomsk.univ-lille1.fr/"
WRITE(*,*) ""
!
!
END SUBROUTINE DISPLAY_HELP_EN
!
!
!
!********************************************************
! ATOMSK_CREATE_DATE
! This routine 
!********************************************************
SUBROUTINE ATOMSK_CREATE_DATE(VALUES,username,msg)
!
IMPLICIT NONE
CHARACTER(LEN=128),INTENT(IN):: username
INTEGER,DIMENSION(8),INTENT(IN):: VALUES
CHARACTER(LEN=128),INTENT(OUT):: msg
!
WRITE(msg,'(i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)')  &
  & VALUES(1), "-", VALUES(2),"-", VALUES(3)," ", VALUES(5), ":", VALUES(6), ":", VALUES(7)
!
msg = 'File generated with Atomsk by '//TRIM(ADJUSTL(username))//' on '//TRIM(ADJUSTL(msg))
!
END SUBROUTINE ATOMSK_CREATE_DATE
!
!
!
!********************************************************
! ATOMSK_MSG
! This routine displays the message corresponding to
! a given index. These messages include warning and
! error messages.
!
! Strings and real numbers can be passed to this routine
! as arrays of dimension 1, for example call:
!      ATOMSK_MSG(1, (/"string"/) ,(/1.d0, 2.d0/) )
! If no string or number must be passed, simply use
! empty array, e.g.:
!      ATOMSK_MSG(1, (/""/) ,(/0.d0/) )
! In each error message below, a comment specifies
! which parameters must be passed as strings and/or
! real numbers for the message to be displayed properly.
!
! Indices of the messages are defined as follows:
!         0- 999: general messages
!      1000-1999: messages for input modules
!      2000-2999: messages for options modules
!      3000-3999: messages for output modules
!      4000-4999: messages for modes
! In each of them:
!      xxx0-x499: various messages
!      x700-x799: warning messages
!      x800-x899: error messages
!      x900-x999: debug messages
!********************************************************
SUBROUTINE ATOMSK_MSG_EN(imsg,strings,reals)
!
IMPLICIT NONE
CHARACTER(LEN=2):: species
CHARACTER(LEN=128):: msg  !The message to be displayed
CHARACTER(LEN=128):: temp, temp2, temp3, temp4
CHARACTER(LEN=*),DIMENSION(:):: strings !Character strings that may be part of the message
INTEGER:: i, j
INTEGER,INTENT(IN):: imsg  !index of message to display
REAL(dp):: tempreal
REAL(dp),DIMENSION(:),INTENT(IN):: reals  !real numbers that may be part of the message
!
!
SELECT CASE(imsg)
!
!**************************
! 1- 999: general messages
CASE(1)
  !Message at program termination
  !nerr and nwarn are global variables
  !reals(1) = total time
  !reals(2) = CPU time
  msg = "\o/ Program terminated successfully!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(temp,"(f30.3)") reals(1)
  WRITE(temp2,"(f30.3)") reals(2)
  WRITE(msg,*) "   Total time: "//TRIM(ADJUSTL(temp))//         &
           & " s.; CPU time: "//TRIM(ADJUSTL(temp2))//" s."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  IF( nwarn>0 .OR. nerr>0 ) THEN
    !In case of warning or error, display a big box
    msg = " _________________________________________"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    IF( nwarn>0 ) THEN
      WRITE(temp,*) nwarn
      temp = ADJUSTL(temp)
      msg = "|  /!\ WARNINGS: "//TRIM(temp)
      msg = msg(1:42)//"|"
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ENDIF
    IF( nerr>0 ) THEN
      WRITE(temp,*) nerr
      temp = ADJUSTL(temp)
      msg = "|  X!X ERRORS:   "//TRIM(temp)
      msg = msg(1:42)//"|"
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ENDIF
    msg = "|_________________________________________|"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
  !
CASE(2)
  !Message giving the elasped time
  !reals(1) = total time
  !reals(2) = CPU time
  WRITE(temp,"(f30.3)") reals(1)
  WRITE(temp2,"(f30.3)") reals(2)
  WRITE(msg,*) "   Total time: "//TRIM(ADJUSTL(temp))//         &
           & " s.; CPU time: "//TRIM(ADJUSTL(temp2))//" s."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3)
  msg = "..> This may take some time..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4)
  !strings(1) = name of file that already exists
  msg = "<?> This file already exists: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Do you want to overwrite it ("//langyes//"/"//langno//") ("&
      & //langBigYes//"=overwrite all)?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5)
  !strings(1) = name of file
  msg = "..> OK, I will overwrite "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(6)
  msg = "..> OK, I will overwrite all files from now on."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(7)
  msg = "<?> Enter a name for the file to write:"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(8)
  !strings(1) = name of file
  IF( LEN_TRIM(strings(1))>0 ) THEN
    msg = "<?> This file does not exist: "//TRIM(strings(1))
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
  msg = "    Please provide the name of an existing file:"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    (type '"//TRIM(system_ls)//"' for a list of files of current directory)"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(9)
  msg = "<?> Enter the name of an existing file:"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(10)
  !reals(1) = index of current element
  !reals(2) = total number of elements
  !SPECIAL: this writes a message on the screen without advancing, so
  ! it is restricted to verbosity levels that display something on screen
  IF( verbosity==1 .OR. verbosity>=3 ) THEN
    temp = ""
    tempreal = 100.d0*reals(1)/reals(2) !percentage of progress
    WRITE(temp2,'(i3)') NINT(tempreal)
    IF( tempreal >=50.d0 ) THEN
      IF(tempreal>=100.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [====================]"
      ELSEIF(tempreal>=95.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [===================>]"
      ELSEIF(tempreal>=90.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [==================> ]"
      ELSEIF(tempreal>=85.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [=================>  ]"
      ELSEIF(tempreal>=80.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [================>   ]"
      ELSEIF(tempreal>=75.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [===============>    ]"
      ELSEIF(tempreal>=70.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [==============>     ]"
      ELSEIF(tempreal>=65.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [=============>      ]"
      ELSEIF(tempreal>=60.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [============>       ]"
      ELSEIF(tempreal>=55.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [===========>        ]"
      ELSE
        temp = TRIM(ADJUSTL(temp2))//"% [==========>         ]"
      ENDIF
    ELSE  
      IF(tempreal<5.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [>                   ]"
      ELSEIF(tempreal<10.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [=>                  ]"
      ELSEIF(tempreal<15.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [==>                 ]"
      ELSEIF(tempreal<20.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [===>                ]"
      ELSEIF(tempreal<25.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [====>               ]"
      ELSEIF(tempreal<30.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [=====>              ]"
      ELSEIF(tempreal<35.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [======>             ]"
      ELSEIF(tempreal<40.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [=======>            ]"
      ELSEIF(tempreal<45.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [========>           ]"
      ELSEIF(tempreal<50.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//"% [=========>          ]"
      ELSE
        temp = TRIM(ADJUSTL(temp2))//"% [==========>         ]"
      ENDIF
    ENDIF
    !Display the progress bar
    WRITE(*,'(a)',ADVANCE="NO") CHAR(13)//" "//TRIM(temp)
    IF( NINT(reals(1)) == NINT(reals(2)) ) THEN
      WRITE(*,*) ""
    ENDIF
  ENDIF
  !
CASE(11)
  msg = ">>> Constructing neighbor list..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(12)
  msg = ">>> If you use Atomsk in your work, please cite the following article:"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "    Pierre Hirel, Comput. Phys. Comm. 197 (2015) 212"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(14)
  msg = "<i> Selected atoms were removed: the selection previously defined was cleared."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
! 700- 799: WARNING MESSAGES
CASE(700)
  !strings(1) = first file name
  !strings(2) = second file name
  msg = "/!\ Both files exist: "//TRIM(strings(1))//" and "//TRIM(strings(2))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "<?> Please indicate which one you want to use as input file:"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = '    1- '//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = '    2- '//TRIM(strings(2))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(701)
  !strings(1) = first file name
  !strings(2) = second file name
  msg = "/!\ None of these files exist: "//TRIM(strings(1))//" nor "//TRIM(strings(2))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "<?> Please provide the name of an existing file:"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(702)
  msg = "/!\ No input file was specified."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "<?> Please provide the name of an existing file:"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(703)
  !strings(1) = name of command line argument
  msg = "/!\ WARNING: Unrecognized command-line argument: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(704)
  !strings(1) = name of file that will be ignored
  msg = "/!\ WARNING: too many file names passed as arguments."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "            The following file will be ignored: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(750)
  msg = ""
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "/!\ WARNING: YOU SHOULD NOT RUN THIS PROGRAM AS ROOT!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  temp=""
  DO WHILE(temp.NE."ok")
    msg = "    Enter 'ok' to continue anyway, or Ctrl+C to quit."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    READ(*,*) temp
  ENDDO
  msg = ""
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
! 800- 899: ERROR MESSAGES
CASE(800)
  msg = "X!X ERROR: unrecognized file format."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(801)
  !strings(1) = atom type
  msg = "X!X ERROR: unrecognized atom species: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(802)
  !reals(1) = index of atom that caused the error
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = "X!X ERROR while trying to read atom #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(803)
  !strings(1) = unit
  msg = "X!X ERROR: unknown unit: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(804)
  msg = "X!X ERROR: number of atoms is zero, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(806)
  msg = "X!X ERROR: the systems do not have the same number of atoms!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(807)
  !reals(1) = index of line in file that caused the error on read
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = "X!X ERROR while trying to read line #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(808)
  !strings(1) = string where conversion failed
  msg = "X!X ERROR: failed to convert '"//TRIM(strings(1))// &
      & "' to numerical value."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(809)
  !strings(1) = space group H-M symbol which couldn't be identified
  msg = "X!X ERROR: failed to identify the space group '"// &
      & TRIM(strings(1))//"'."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(810)
  !reals(1) = space group number
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = "X!X ERROR: invalid space group: "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(811)
  !reals(1) = space group number
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = "X!X ERROR: no access to data of space group "// &
      & TRIM(ADJUSTL(msg))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(812)
  !strings(1) = non-conformal symmetry operation string
  msg = "X!X ERROR: failed to interpret the symmetry operations '"// &
      & TRIM(strings(1))//"'."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(813)
  msg = "X!X ERROR: no such file or directory"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(814)
  msg = "X!X ERROR: Miller vector cannot be [000]."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(815)
  msg = "X!X ERROR: in hexagonal lattices, Miller indices [hkil] must have h+k+i=0."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(816)
  msg = "X!X ERROR: division by zero."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(817)
  !strings(1) = a string that could not be interpreted
  msg = "X!X ERROR: unable to interpret this string as Miller indices: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(818)
  !strings(1) = name of the array that was supposed to be resized
  msg = "X!X ERROR: there was an error while attempting to resize the array "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
!
! 900- 999: DEBUG MESSAGES
CASE(999)
  !strings(1) = debug message itself
  msg = "debug: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!
!
!**************************
! 1000-1999: MESSAGES FOR INPUT
CASE(1000)
  !strings(1) = file name
  msg = ">>> Opening the input file: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(1001)
  !reals(1) = number of atoms
  !reals(2) = number of shells
  WRITE(temp,*) NINT(reals(1))
  IF( NINT(reals(2))>0 ) THEN
    WRITE(temp2,*) NINT(reals(2))
    msg = "..> Input file was read successfully ("//TRIM(ADJUSTL(temp))//" cores + "// &
        & TRIM(ADJUSTL(temp2))//" shells)."
  ELSE
    msg = "..> Input file was read successfully ("//TRIM(ADJUSTL(temp))//" atoms)."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(1002)
  msg = "..> Found INP file, reading supercell from it..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(1003)
  msg = "..> Found POTCAR file, reading atom species from it..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(1004)
  msg = "..> Applying symmetry operations..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!1700-1799: WARNING MESSAGES
CASE(1700)
  !strings(1) = auxiliary property that cannot be loaded
  msg = "/!\ WARNING: cannot load auxiliary property: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1701)
  !strings(1) = input file format
  msg = "/!\ WARNING: the most probable file format was found to be "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    However the format could not be determined with certainty."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1702)
  !strings(1) = name of parameter
  !strings(2) = name of custom config file
  msg = "/!\ WARNING: unknown parameter '"//TRIM(strings(1))//&
      & "' in "//TRIM(strings(2))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1703)
  !strings(1) = name of personal config file
  msg = "/!\ WARNING: errors were encountered while reading configuration file "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Some personal parameters may not have been properly set."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1704)
  !strings(1) = name of personal config file
  msg = "/!\ WARNING: symmetry operations are not taken into account."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1705)
  msg = "/!\ WARNING: both celldm(:) and conventional notation were found."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            The celldm(:) will be used, and conventional notation ignored."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1706)
  msg = "/!\ WARNING: cell dimensions are in Bohrs, while atom positions are in angströms."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Cell vectors will be converted to angströms for consistency."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1707) ! invalid symmetry operation string input.
  !strings(1) = failed symmetry operation string
  msg = "/!\ WARNING: invalid symmetry operation string '"// &
      & TRIM(strings(1))//"', skipping..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1799)
  msg = "/!\ WARNING: the data file had an unknown format. Atomsk tried to extract"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            atomic data from it, but it may be wrong. Tread carefully!"
  CALL DISPLAY_MSG(1,msg,logfile)
!
!1800-1899: ERROR MESSAGES
CASE(1800)
  msg = "X!X ERROR: I could not guess the format of this input file!"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Maybe it is a format unsupported by Atomsk yet?"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Anyway it appears I cannot help, I will skip this file."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1801)
  !strings(1) = file name
  !reals(1) = line number
  msg = "X!X ERROR: there were errors while reading the file: " &
      & //TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
  IF( NINT(reals(1))>0 ) THEN
    WRITE(temp,*) NINT(reals(1))
    msg = "          Error appears to be at line # "//TRIM(ADJUSTL(temp))
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(1802)
  !strings(1) = bad array
  msg = "X!X ERROR: inconsistent array size in "//TRIM(strings(1))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1803)
  msg = "X!X ERROR: size of auxiliary properties "// &
             & "is not consistent with number of atoms."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1804)
  msg = "X!X ERROR: unknown format."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1805)
  !reals(1) = number of particles read
  !reals(2) = number of particles declared
  msg = "X!X ERROR: number of atoms read differs from the number"
  CALL DISPLAY_MSG(1,msg,logfile)
  WRITE(temp,*) NINT(reals(1))
  WRITE(msg,*) NINT(reals(2))
  msg = "           of atoms declared: "// &
      & TRIM(ADJUSTL(temp))//"/"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1806)
  !reals(1) = index of atom
  WRITE(msg,*) NINT(reals(1))
  msg = "X!X ERROR: denominator is zero in coordinate of atom #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1807)
  msg = "X!X ERROR: file is not in ASCII format, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1808)
  msg = "X!X ERROR: levcfg cannot be greater than 2."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1809)
  msg = "X!X ERROR: unable to read the supercell parameters."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1810)
  msg = "X!X ERROR: unable to read the number of atoms."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1811)
  !reals(1) = index of atom
  WRITE(msg,*) NINT(reals(1))
  IF( reals(1)<0.1d0 ) THEN
    msg = "X!X ERROR: atom index #"//TRIM(ADJUSTL(msg))//" is negative."
  ELSE
    msg = "X!X ERROR: atom index #"//TRIM(ADJUSTL(msg))//" is greater than the number of atoms."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1812)
  !reals(1) = index of atom
  WRITE(msg,*) NINT(reals(1))
  msg = "X!X ERROR: unable to read properties of atom #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1813)
  !strings(1) = string where conversion failed
  msg = "X!X ERROR: failed to parse symmetry operation in '"// &
      & TRIM(strings(1))//"'."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1814)
  !strings(1) = name of file format
  !strings(2) = name of mode to use to read this file format
  IF( LEN_TRIM(strings(2))>0 ) THEN
    msg = "X!X ERROR: files in "//TRIM(ADJUSTL(strings(1)))//  &
        & " format can only be read with the mode "//TRIM(ADJUSTL(strings(2)))//"."
  ELSE
    msg = "X!X ERROR: files in "//TRIM(ADJUSTL(strings(1)))//" format cannot be read."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
!
!
!
!**************************
! 2000-2999: MESSAGES FOR OPTIONS
!
CASE(2050)
  msg = ">>> Aligning cell vectors..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2051)
  msg = "..> Cell vectors were successfully aligned."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2052)
  msg = ">>> Converting to fractional coordinates..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2053)
  msg = "..> Coordinates were reduced."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2054)
  !strings(1) = atom species
  msg = ">>> Creating shells for "//TRIM(strings(1))//" atoms."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2055)
  !reals(1) = number of shells
  IF( NINT(reals(1))==0 ) THEN
    msg = "..> No shell was added."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 shell was added to the system."
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" shells were added to the system."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2056)
  !string(1) = cut_dir, i.e. "above" or "below"
  !string(2) = cutdir, i.e. x, y or z
  !reals(1) = cutdistance in angstroms
  IF( DABS(reals(1))<1.d12 ) THEN
    WRITE(msg,"(3f16.3)") reals(1)
  ELSEIF( reals(1)<-1.d12 ) THEN
    WRITE(msg,"(a4)") "-INF"
  ELSEIF( reals(1)>1.d12 ) THEN
    WRITE(msg,"(a4)") "+INF"
  ENDIF
  msg = ">>> Cutting the system "//TRIM(strings(1))//" "// &
    & TRIM(ADJUSTL(msg))//" A along "//TRIM(strings(2))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2057)
  !reals(1) = NPcut, number of deleted atoms
  !reals(2) = number of atoms left
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 atom was wiped out"
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" atoms were wiped out"
  ENDIF
  WRITE(temp,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(msg))//", "//TRIM(ADJUSTL(temp))//" atoms left."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2058)
  !reals(1) = deformation
  !strings(1) = direction of deformation: x, y or z
  WRITE(msg,"(f16.3)") reals(1)*100.d0
  msg = ">>> Deforming the system by "//TRIM(ADJUSTL(msg))//"% along "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2059)
  !reals(1) = Poisson coefficient
  IF(reals(1)==0.d0) THEN
    msg = "..> Uniaxial strain."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSE
    WRITE(msg,"(f16.3)") reals(1)
    msg = "..> Uniaxial stress, Poisson ratio: "//TRIM(ADJUSTL(msg))
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
CASE(2060)
  msg = "..> System was successfully deformed."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2061)
  !strings(1) = disloctype: screw, edge, edge2, mixed, loop
  !strings(2) = direction of dislocline: x, y or z
  !reals(1) = X component of Burgers vector
  !reals(2) = Y component of Burgers vector
  !reals(3) = Z component of Burgers vector
  !reals(4) = positive if C_tensor is defined, zero otherwise
  !reals(5) = pos1, position of dislocation along first axis
  !reals(6) = pos2, position of dislocation along second axis
  !reals(7) = pos3 (only for loops)
  !reals(8) = loop radius
  temp = TRIM(ADJUSTL(strings(1)))
  IF(TRIM(temp)=="screw") THEN
    msg = ">>> Inserting a screw dislocation with line along"
  ELSEIF(temp(1:4)=="edge") THEN
    msg = ">>> Inserting an edge dislocation with line along"
  ELSEIF(temp(1:5)=="mixed") THEN
    msg = ">>> Inserting a mixed dislocation with line along"
  ELSEIF(temp(1:4)=="loop") THEN
    msg = ">>> Inserting a dislocation loop in a plane normal to"
  ENDIF
  msg = TRIM(msg)//' '//TRIM(strings(2))//","
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  !
  IF( reals(4)>0.1d0 ) THEN
    msg = "    using anisotropic elasticity,"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
  !
  IF(TRIM(strings(1))=="edge") THEN
    WRITE(msg,"(a34)") "    by inserting a plane of atoms,"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF(TRIM(strings(1))=="edge2") THEN
    WRITE(msg,"(a41)") "    conserving the total number of atoms,"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
  !
  WRITE(msg,"(f16.3)") reals(1)
  WRITE(temp,"(f16.3)") reals(2)
  WRITE(temp2,"(f16.3)") reals(3)
  msg = "["//TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(temp))//" "//TRIM(ADJUSTL(temp2))//"]"
  WRITE(temp,"(f16.3)") reals(5)
  WRITE(temp2,"(f16.3)") reals(6)
  IF( TRIM(ADJUSTL(strings(1)))=="loop" ) THEN
    WRITE(temp3,"(f16.3)") reals(7)
    IF( reals(8)>0.d0 ) THEN
      WRITE(temp4,"(f16.3)") reals(8)
      msg = "    Center ("//TRIM(ADJUSTL(temp))//","//TRIM(ADJUSTL(temp2))//","// &
          & TRIM(ADJUSTL(temp3))//"); Radius "//TRIM(ADJUSTL(temp4))//" A; b="//TRIM(ADJUSTL(msg))
    ELSE
      WRITE(temp4,"(f16.3)") DABS(reals(8))
      msg = "    Center ("//TRIM(ADJUSTL(temp))//","//TRIM(ADJUSTL(temp2))//","// &
          & TRIM(ADJUSTL(temp3))//"); Side "//TRIM(ADJUSTL(temp4))//" A; b="//TRIM(ADJUSTL(msg))
    ENDIF
  ELSE
    msg = "    b="//TRIM(ADJUSTL(msg))//" at ("// &
        & TRIM(ADJUSTL(temp))//","//TRIM(ADJUSTL(temp2))//")"
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2062)
  msg = "..> Searching the solutions to the anisotropic elasticity equations..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2063)
  !reals(1) = number of inserted atoms
  WRITE(msg,*) NINT(reals(1))
  msg = "..> "//TRIM(ADJUSTL(msg))//" atoms were inserted."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2064)
  msg = "..> Supercell was expanded."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2065)
  msg = "..> Dislocation was successfully created."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2066)
  !reals(1) = number of repetitions along X
  !reals(2) = number of repetitions along Y
  !reals(3) = number of repetitions along Z
  WRITE(temp,*) NINT( reals(1) )
  WRITE(temp2,*) NINT( reals(2) )
  WRITE(msg,*) NINT( reals(3) )
  msg = ">>> Duplicating the system: "//TRIM(ADJUSTL(temp))//" x "// &
    & TRIM(ADJUSTL(temp2))//" x "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2067)
  !reals(1) = new number of particles
  WRITE(msg,*) NINT( reals(1) )
  msg = "..> New number of particles: "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2068)
  !reals(1) = new number of atoms
  WRITE(temp,*) NINT( reals(1) )
  msg = "..> System was successfully duplicated ("//TRIM(ADJUSTL(temp))//" atoms)."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2069)
  !strings(1) = "all" or property to be removed
  IF(strings(1)=="all") THEN
    WRITE(msg,'(a)') ">>> Removing all auxiliary properties."
  ELSE
    msg = ">>> Removing auxiliary property: "//TRIM(strings(1))
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2070)
  !strings(1) = "all" or property to be removed
  IF(strings(1)=="all") THEN
    msg = "..> Auxiliary properties were successfully removed."
  ELSE
    msg = "..> Auxiliary property was successfully removed."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2071)
  msg = ">>> Orienting the system..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2072)
  msg = "..> System was successfully oriented."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2073)
  !strings(1) = file containing properties
  msg = ">>> Reading system properties from "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2074)
  !strings(1) = property that was read
  msg = "..> Property '"//TRIM(ADJUSTL(strings(1)))//"' was read."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2075)
  msg = "..> Finished reading system properties."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2076)
  !msg = ">>> Selecting all atoms."
  !CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2077)
  !strings(1) = side of selected region: in or out
  !strings(2) = region geometry: sphere or box or cylinder
  !strings(3) = axis of cylinder
  !reals(1) = region_1(1)
  !reals(2) = region_1(2)
  !reals(3) = region_1(3)
  !reals(4) = region_2(1)
  !reals(5) = region_2(2)
  !reals(6) = region_2(3)
  IF( strings(1)=="all" ) THEN
    msg = ">>> Selecting all atoms."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="index" ) THEN
    IF( strings(2)=="list " ) THEN
      msg = ">>> Selecting a list of atoms."
    ELSEIF( strings(2)=="range" ) THEN
      WRITE(temp,*) NINT(reals(1))
      WRITE(temp2,*) NINT(reals(2))
      msg = ">>> Selecting atoms #"//TRIM(ADJUSTL(temp))//" to "//TRIM(ADJUSTL(temp2))//"."
    ELSE
      WRITE(temp,*) NINT(reals(1))
      msg = ">>> Selecting atom #"//TRIM(ADJUSTL(temp))//"."
    ENDIF
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="invert" ) THEN
    msg = ">>> Inverting the selection of atoms..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="above" .OR. strings(1)=="below" ) THEN
    IF( DABS(reals(1))<1.d12 ) THEN
      WRITE(temp,"(3f16.3)") reals(1)
    ELSEIF( reals(1)<-1.d12 ) THEN
      WRITE(temp,"(a4)") "-INF"
    ELSEIF( reals(1)>1.d12 ) THEN
      WRITE(temp,"(a4)") "+INF"
    ENDIF
    msg = ">>> Selecting atoms "//TRIM(strings(1))//" "//TRIM(ADJUSTL(temp))// &
        & " A along the "//TRIM(strings(3))//" axis."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="in" .OR. strings(1)=="out" ) THEN
    IF(strings(1)=="in") THEN
      temp = "inside the"
    ELSE
      temp = "outside of the"
    ENDIF
    IF(TRIM(strings(2))=="cell") THEN
      temp = TRIM(ADJUSTL(temp))//" cell"
    ENDIF
    msg = ">>> Selecting atoms "//TRIM(temp)//" "//TRIM(strings(2))//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    IF(TRIM(strings(2))=="box") THEN
      msg = "..> Box bounds: ("
      DO i=1,3
        IF( DABS(reals(i))<1.d12 ) THEN
          WRITE(temp,"(3f16.3)") reals(i)
        ELSEIF( reals(i)<-1.d12 ) THEN
          WRITE(temp,"(a4)") "-INF"
        ELSEIF( reals(i)>1.d12 ) THEN
          WRITE(temp,"(a4)") "+INF"
        ENDIF
        IF(i==1.OR.i==2) temp = TRIM(ADJUSTL(temp))//","
        msg = TRIM(msg)//TRIM(ADJUSTL(temp))
      ENDDO
      msg = TRIM(msg)//"); ("
      DO i=4,6
        IF( DABS(reals(i))<1.d12 ) THEN
          WRITE(temp,"(3f16.3)") reals(i)
        ELSEIF( reals(i)<-1.d12 ) THEN
          WRITE(temp,"(a4)") "-INF"
        ELSEIF( reals(i)>1.d12 ) THEN
          WRITE(temp,"(a4)") "+INF"
        ENDIF
        IF(i==4.OR.i==5) temp = TRIM(ADJUSTL(temp))//","
        msg = TRIM(msg)//TRIM(ADJUSTL(temp))
      ENDDO
      msg = TRIM(msg)//")."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ELSEIF(TRIM(strings(2))=="sphere") THEN
      WRITE(temp,"(f16.3)") reals(1)
      msg = "..> Center: ("//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(2)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(3)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))//")"
      WRITE(temp,"(f16.3)") reals(4)
      msg = TRIM(msg)//"; Radius: "//TRIM(ADJUSTL(temp))//" A."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ELSEIF(TRIM(strings(2))=="cylinder") THEN
      msg = "..> Axis along "//TRIM(strings(3))//"."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
      WRITE(temp,"(f16.3)") reals(1)
      msg = "..> Center: ("//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(2)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))//")"
      WRITE(temp,"(f16.3)") reals(4)
      msg = TRIM(msg)//"; Radius: "//TRIM(ADJUSTL(temp))//" A."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ELSEIF(TRIM(strings(2))=="torus") THEN
      WRITE(temp,"(f16.3)") reals(1)
      msg = "..> Axis along "//TRIM(strings(3))//". Center: ("//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(2)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(3)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))//")."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
      WRITE(temp,"(f16.3)") reals(4)
      WRITE(temp2,"(f16.3)") reals(5)
      msg = "..> Main radius: "//TRIM(ADJUSTL(temp))//" A; secondary radius: "//TRIM(ADJUSTL(temp2))//" A."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ENDIF
  ELSEIF( strings(1)=="prop" ) THEN
    msg = ">>> Selecting atoms with "//TRIM(strings(2))
    WRITE(temp,"(f16.3)") reals(1)
    IF( reals(4)>2.d0 ) THEN
      WRITE(temp2,"(f16.3)") reals(2)
      msg = TRIM(ADJUSTL(msg))//" between "//TRIM(ADJUSTL(temp))//" and "//TRIM(ADJUSTL(temp2))//"."
    ELSE
      msg = TRIM(ADJUSTL(msg))//" equal to "//TRIM(ADJUSTL(temp))//"."
    ENDIF
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="random" .OR. strings(1)=="rand" ) THEN
    WRITE(temp,*) NINT(reals(1))
    msg = ">>> Randomly selecting "//TRIM(ADJUSTL(temp))
    IF( NINT(reals(1))<=1 ) THEN
      msg = TRIM(msg)//" atom"
    ELSE
      msg = TRIM(msg)//" atoms"
    ENDIF
    IF( strings(2).NE."any" .AND. strings(2).NE."all" ) THEN
      msg = TRIM(msg)//" of "//TRIM(strings(2))
    ENDIF
    msg = TRIM(msg)//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="random%" ) THEN
    WRITE(temp,'(f16.3)') reals(1)*100.d0
    msg = ">>> Randomly selecting "//TRIM(ADJUSTL(temp))//"% of atoms"
    IF( strings(2).NE."any" .AND. strings(2).NE."all" ) THEN
      msg = TRIM(msg)//" of "//TRIM(strings(2))
    ENDIF
    msg = TRIM(msg)//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="neigh" ) THEN
    !strings(2) = species of neighbors
    !reals(1) = number of neighbors or cutoff radius for neighbor search
    !reals(2) = species of central atom (i.e. atom whose neighbors are searched for)
    IF( strings(2)=="all" .OR. strings(2)=="any" ) THEN
      temp2 = " "
    ELSE
      temp2 = " "//TRIM(ADJUSTL(strings(2)))
    ENDIF
    !
    IF( reals(1)==0.d0 ) THEN
      msg = "first nearest"//TRIM(temp2)//" neighbors"
    ELSEIF( reals(1)>0.d0 .AND. DBLE(NINT(reals(1)))-reals(1)<1.d-12 ) THEN
      WRITE(temp,*) NINT(reals(1))
      IF( NINT(reals(1))==1 ) THEN
        msg = "the first"//TRIM(temp2)//" neighbor"
      ELSE
        msg = "the "//TRIM(ADJUSTL(temp))//" nearest"//TRIM(temp2)//" neighbors"
      ENDIF
    ELSE
      WRITE(temp,'(f16.3)') DABS(reals(1))
      msg = TRIM(ADJUSTL(temp2))//" neighbors within a radius of "//TRIM(ADJUSTL(temp))//" A"
    ENDIF
    WRITE(temp,*) NINT(reals(2))
    msg = ">>> Selecting "//TRIM(ADJUSTL(msg))//" of atom #"//TRIM(ADJUSTL(temp))//"..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="stl" ) THEN
    msg = ">>> Selecting atoms inside the 3-D shape from the STL file: "//TRIM(ADJUSTL(strings(2)))
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSE
    !Last case: strings(1) should be an atom species
    msg = ">>> Selecting all "//TRIM(ADJUSTL(strings(1)))//" atoms..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
CASE(2078)
  !reals(1) = number of atoms that were selected
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 atom was selected."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(temp,*) NINT( reals(1) )
    msg = "..> "//TRIM(ADJUSTL(temp))//" atoms were selected."
  ELSE
    msg = "..> No atom was selected."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2079)
  !reals(1) = cutoff radius for removing atoms
  WRITE(msg,"(f16.3)") reals(1)
  msg = ">>> Removing atoms that are closer than "//TRIM(ADJUSTL(msg))//" A."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2080)
  !strings(1) = species of removed atom(s) (may be empty or contain anything)
  !reals(1) = number of atoms removed
  !reals(2) = number of atoms left
  species = strings(1)
  CALL ATOMNUMBER(species,tempreal)
  IF(NINT(reals(1))==0) THEN
    msg = "..> No atom was removed"
  ELSEIF(NINT(reals(1))==1) THEN
    IF(tempreal>0.5d0) THEN
      msg = "..> 1 atom of "//TRIM(ADJUSTL(species))//" was removed"
    ELSE
      msg = "..> 1 atom was removed"
    ENDIF
  ELSE
    WRITE(msg,*) NINT(reals(1))
    IF(tempreal>0.5d0) THEN
      msg = "..> "//TRIM(ADJUSTL(msg))//" atoms of "//TRIM(ADJUSTL(species))//" were removed"
    ELSE
      msg = "..> "//TRIM(ADJUSTL(msg))//" atoms were removed"
    ENDIF
  ENDIF
  WRITE(temp,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(msg))//", "//TRIM(ADJUSTL(temp))//" atoms left."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2081)
  !strings(1) = rotation axis: x, y or z
  !reals(1) = rotation angle in degrees
  WRITE(msg,"(f16.2)") reals(1)
  msg = ">>> Rotating the system by "//TRIM(ADJUSTL(msg)) &
      & //"° around "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2082)
  msg = "..> System was successfully rotated."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2083)
  !strings(1) = axis that is tilted to produce shear: x, y or z
  !strings(2) = direction of tilt: x, y or z
  !reals(1) = magnitude of tilt in Angstroms
  WRITE(temp,"(f16.3)") reals(1)
  msg = ">>> Tilting the "//TRIM(strings(1))//" axis by "// &
      & TRIM(ADJUSTL(temp))//" A along "//TRIM(strings(2))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2084)
  !reals(1) = shear strain in %
  msg = "..> System was successfully sheared."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,"(f16.3)") reals(1)
  msg = "..> Applied shear strain = "//TRIM(ADJUSTL(msg))//" %."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2085)
  !strings(1) = shift direction: above or below (or empty)
  !strings(2) = shift axis
  !reals(1) = distance above axis
  !reals(2), reals(3), reals(4) = shift vector
  WRITE(temp,"(f16.3)") reals(2)
  WRITE(temp2,"(f16.3)") reals(3)
  WRITE(temp3,"(f16.3)") reals(4)
  IF( strings(1)(1:5)=='above' .OR. strings(1)(1:5)=='below' ) THEN
    IF( DABS(reals(1))<1.d12 ) THEN
      WRITE(msg,"(3f16.3)") reals(1)
    ELSEIF( reals(1)<-1.d12 ) THEN
      WRITE(msg,"(a4)") "-INF"
    ELSEIF( reals(1)>1.d12 ) THEN
      WRITE(msg,"(a4)") "+INF"
    ENDIF
    msg = ">>> Shifting atoms "//TRIM(strings(1))//" "//TRIM(ADJUSTL(msg))//  &
          & "A along "//TRIM(strings(2))//" by ("//TRIM(ADJUSTL(temp))//","//     &
          & TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))//")"
  ELSEIF( strings(1)=="select" ) THEN
    msg = ">>> Shifting selected atoms by ("//TRIM(ADJUSTL(temp))//","//     &
          & TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))//")"
  ELSE
    msg = ">>> Shifting all atoms by ("//TRIM(ADJUSTL(temp))//","//     &
          & TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))//")"
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2086)
  !reals(1) = number of atoms that were shifted
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 atom was shifted."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" atoms were shifted."
  ELSE
    msg = "..> No atom was shifted."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2087)
  !strings(1) = property to be sorted
  !strings(2) = sort order: up, down, pack
  IF(strings(1)=="s") THEN
    temp = "atomic weight"
  ELSEIF( strings(1)=="x" .OR. strings(1)=="y" .OR. strings(1)=="z" .OR. &
        & strings(1)=="X" .OR. strings(1)=="Y" .OR. strings(1)=="Z"  ) THEN
    temp = TRIM(strings(1))//" coordinate"
  ELSE
    temp = TRIM(strings(1))
  ENDIF
  IF(strings(2)=="up") temp = "increasing "//TRIM(temp)
  IF(strings(2)=="down") temp = "decreasing "//TRIM(temp)
  IF(strings(2)=="up" .OR. strings(2)=="down") THEN
    msg = ">>> Sorting atoms by "//TRIM(ADJUSTL(temp))//"."
  ELSE
    msg = ">>> Packing atoms by "//TRIM(ADJUSTL(temp))//"."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2088)
  msg = "..> Atoms were successfully sorted."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2089)
  !strings(1) = first atomic species
  !strings(2) = second atomic species
  msg = ">>> Substituting "//TRIM(strings(1))//" atoms with "//TRIM(strings(2))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2090)
  !reals(1) = number of substituted atoms
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 atom was substituted."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" atoms were substituted."
  ELSE
    msg = "..> No atom was substituted."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2091)
  !strings(1) = what to convert / property name
  !strings(2) = first unit of distance
  !strings(3) = second unit of distance
  !strings(4) = first unit of time
  !strings(5) = second unit of time
  !reals(1) = 0 or factor
  IF( strings(1)=="velocities" .OR. strings(1)=="coordinates" ) THEN
    IF( LEN_TRIM(strings(4))>0 ) THEN
      temp = TRIM(strings(2))//"/"//TRIM(strings(4))
    ELSE
      temp = TRIM(strings(2))
    ENDIF
    IF( LEN_TRIM(strings(5))>0 ) THEN
      temp2 = TRIM(strings(3))//"/"//TRIM(strings(5))
    ELSE
      temp2 = TRIM(strings(3))
    ENDIF
    msg = ">>> Converting "//TRIM(strings(1))//" from "//TRIM(temp)//&
        & " to "//TRIM(temp2)//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSE
    WRITE(temp,'(f16.3)') reals(1)
    msg = ">>> Rescaling "//TRIM(strings(1))//" by a factor "//TRIM(ADJUSTL(temp))//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
CASE(2092)
  !strings(1) = what was converted
  IF( strings(1)=="velocities" .OR. strings(1)=="coordinates" ) THEN
    msg = "..> "//TRIM(ADJUSTL(strings(1)))//" were converted."
  ELSE
    msg = "..> "//TRIM(ADJUSTL(strings(1)))//" were rescaled."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2093)
  msg = ">>> Wrapping atoms..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2094)
  !reals(1) = number of atoms wrapped
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 atom was wrapped."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" atoms were wrapped."
  ELSE
    msg = "..> No atom was wrapped."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2095)
  msg = "..> Cell vectors were computed."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2096)
  msg = "..> The solutions were found."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2097)
  !strings(1) = axis along which coordinate is fixed: x,y,z,all
  !strings(2) = "above" or "below" or "select" or nothing
  !strings(3) = axis 
  !reals(1) = fix distance in angstroms
  IF( DABS(reals(1))<1.d12 ) THEN
    WRITE(msg,"(3f16.3)") reals(1)
  ELSEIF( reals(1)<-1.d12 ) THEN
    WRITE(msg,"(a4)") "-INF"
  ELSEIF( reals(1)>1.d12 ) THEN
    WRITE(msg,"(a4)") "+INF"
  ENDIF
  IF( strings(2)=='above' .OR. strings(2)=='below' ) THEN
    msg = ">>> Fixing "//TRIM(strings(1))//" coordinate of atoms "//&
        & TRIM(strings(2))//" "//TRIM(ADJUSTL(msg))//"A along "//TRIM(strings(3))//"."
  ELSEIF( strings(2)=='selec' ) THEN
    msg = ">>> Fixing "//TRIM(strings(1))//" coordinate of selected atoms."
  ELSE
    msg = ">>> Fixing "//TRIM(strings(1))//" coordinate of all atoms."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2098)
  !reals(1) = number of atoms that were fixed
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 atom was fixed."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" atoms were fixed."
  ELSE
    msg = "..> No atom was fixed."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2099)
  msg = "..> Elastic tensor was rotated."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2100)
  !reals(1) = anisotropy ratio A = 2*C44 / (C11-C12)
  !reals(2) = anisotropy factor H = 2*C44 + C12 - C11
  WRITE(msg,'(f16.3)') reals(1)
  msg = "..> Anisotropy ratio: A = "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,'(e16.3)') reals(2)
  msg = "..> Anisotropy factor: H = "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2101)
  !strings(1) = formula of energy factor, e.g. "Kb²/4pi"
  !reals(1) = energy factor
  msg = "..> Dislocation stresses were computed."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,'(f24.8)') reals(1)
  msg = "..> Prelogarithmic energy factor: "//TRIM(ADJUSTL(strings(1)))//&
      & " = "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2102)
  !strings(1) = empty or atom species to remove or "SEL"
  !reals(1) = 0 or atom index to remove
  IF( strings(1)=="SEL" ) THEN
    msg = ">>> Removing all selected atoms from the system."
  ELSE
    IF( NINT(reals(1)).NE.0 ) THEN
      WRITE(msg,*) NINT(reals(1))
      msg = ">>> Removing atom #"//TRIM(ADJUSTL(msg))//"."
    ELSE
      msg = ">>> Removing all "//TRIM(strings(1))//" atoms from the system."
    ENDIF
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2103)
  msg = ">>> Reducing the tilt of the supercell vectors..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2104)
  !reals(1) = number of components that were unskewed
  IF( NINT(reals(1))==0 ) THEN
    msg = "..> No supercell vector was unskewed."
  ELSE
    msg = "..> Supercell was unskewed."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2105)
  msg = ">>> Binding ionic shells with their respective cores..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2106)
  !reals(1) = number of shells that were re-assigned
  IF( NINT(reals(1))<=0 ) THEN
    msg = "..> No shell was rebound."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 shell was rebound to its core."
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" shells were rebound to their core."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2107)
  !strings(1) = atom species to be removed
  !reals(1) = number of atoms to be randomly removed
  WRITE(msg,*) NINT(reals(1))
  IF( strings(1)=="any" .OR. strings(1)=="all" ) THEN
    msg = ">>> Randomly removing "//TRIM(ADJUSTL(msg))//" atoms from the system."
  ELSE
    msg = ">>> Randomly removing "//TRIM(ADJUSTL(msg))//" "//TRIM(strings(1))//&
        & " atoms from the system."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2108)
  !strings(1) = crack mode (I, II or III)
  !strings(2) = direction of crack line (X, Y or Z)
  !reals(1) = pos1
  !reals(2) = pos2
  WRITE(msg,"(f16.3)") reals(1)
  WRITE(temp,"(f16.3)") reals(2)
  msg = ">>> Inserting a mode "//TRIM(ADJUSTL(strings(1)))//" crack along "// &
      & TRIM(ADJUSTL(strings(2)))//" at ("//TRIM(ADJUSTL(msg))//","//         &
      & TRIM(ADJUSTL(temp))//")."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2109)
  msg = "..> Crack was successfully created."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2110)
  !msg = ">>> Inverting the selection of atoms..."
  !CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2111)
  !reals(1) = target temperature
  WRITE(msg,'(f24.3)') reals(1)
  msg = ">>> Assigning velocities according to a Maxwell-Boltzmann distribution at "// &
      & TRIM(ADJUSTL(msg))//" K."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2112)
  msg = "..> Atom velocities were successfuly set up."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2113)
  !strings(1) = file name
  msg = "..> Distribution of velocities was written to the file: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2114)
  !reals(1) = max. displacement of an atom along X
  !reals(2) = max. displacement of an atom along Y
  !reals(3) = max. displacement of an atom along Z
  msg = ">>> Applying a perturbation to atom positions,"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  IF( DABS(reals(1)-reals(2))<1.d-12 .AND. DABS(reals(1)-reals(3))<1.d-12 ) THEN
    !All values are equal
    WRITE(msg,'(f24.3)') reals(1)
    msg = "..> maximum magnitude: "//TRIM(ADJUSTL(msg))//" A."
  ELSE
    !Values are different along each direction
    WRITE(temp,'(f24.3)') reals(1)
    WRITE(temp2,'(f24.3)') reals(2)
    WRITE(temp3,'(f24.3)') reals(3)
    msg = "..> maximum magnitude: dx="//TRIM(ADJUSTL(temp))//", dy=" &
        & //TRIM(ADJUSTL(temp2))//", dz="//TRIM(ADJUSTL(temp3))//"."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2115)
  msg = "..> Atom positions were disturbed."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2116)
  !strings(1) = species of added atom(s)
  !strings(2) = "at" or "near" or "random"
  !reals(1) = position x of added atom, or index of atom, or number of atoms added
  !reals(2) = position y
  !reals(3) = position z
  SELECT CASE(strings(2))
  CASE("at","AT","@")
    WRITE(temp,"(f16.3)") reals(1)
    WRITE(temp2,"(f16.3)") reals(2)
    WRITE(temp3,"(f16.3)") reals(3)
    msg = ">>> Adding an atom of "//TRIM(ADJUSTL(strings(1)))//" at ("//TRIM(ADJUSTL(temp))//","//     &
        & TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))//")."
  CASE("near","NEAR")
    WRITE(temp,*) NINT(reals(1))
    msg = ">>> Adding an atom of "//TRIM(ADJUSTL(strings(1)))//" near atom #"//TRIM(ADJUSTL(temp))
  CASE("random","RANDOM","rand","RAND")
    WRITE(temp,*) NINT(reals(1))
    msg = ">>> Adding "//TRIM(ADJUSTL(temp))//" atoms of "//TRIM(ADJUSTL(strings(1)))//" at random positions..."
  END SELECT
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2117)
  !reals(1) = number of atoms added
  !reals(2) = number of atoms in the system
  IF(NINT(reals(1))==0) THEN
    msg = "..> No atom was added"
  ELSEIF(NINT(reals(1))==1) THEN
    msg = "..> 1 atom was added"
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" atoms were added"
  ENDIF
  WRITE(temp,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(msg))//", "//TRIM(ADJUSTL(temp))//" atoms in the system."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2118)
  !reals(1) = index of atom to center; if <=0 then center of mass
  IF( NINT(reals(1))<=0 ) THEN
    msg = ">>> Placing the center of mass at the center of the box..."
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = ">>> Shifting system so that atom #"//TRIM(ADJUSTL(msg))//" is at the center of the box..."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2119)
  !reals(1) = X component of shift vector
  !reals(2) = Y component of shift vector
  !reals(3) = Z component of shift vector
  WRITE(temp,'(f16.3)') reals(1)
  WRITE(temp2,'(f16.3)') reals(2)
  WRITE(temp3,'(f16.3)') reals(3)
  msg = "..> System was re-centered, shift vector: ("//TRIM(ADJUSTL(temp))// &
      & ","//TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))//")."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2120)
  !strings(1) = direction normal to mirror plane
  !reals(1) = distance between mirror plane and origin of coordinates
  WRITE(temp,'(f16.3)') reals(1)
  msg = ">>> Applying a mirror plane at "//TRIM(ADJUSTL(temp))//" A along "//TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2121)
  msg = "..> System was successfully mirrored."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2122)
  !strings(1) = species on which shells are removed
  msg = ">>> Removing shells on "//TRIM(ADJUSTL(strings(1)))//" ions..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2123)
  !reals(1) = number of shells removed
  WRITE(temp,*) NINT(reals(1))
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 shell was removed."
  ELSE
    msg = "..> "//TRIM(ADJUSTL(temp))//" shells were removed."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2124)
  !strings(1) = stress component or name of file
  !reals(1) = value of stress
  WRITE(temp,'(f16.3)') reals(1)
  SELECT CASE(strings(1))
  CASE('x','X','xx','XX','y','Y','yy','YY','z','Z','zz','ZZ','xy','XY','yx','YX','zx','ZX','xz','XZ','zy','ZY','yz','YZ')
    temp2 = ADJUSTL(strings(1))
    IF( LEN_TRIM(temp2)<=1 ) THEN
      temp2 = TRIM(ADJUSTL(temp2))//TRIM(ADJUSTL(temp2))
    ENDIF
    msg = ">>> Applying a stress σ_"//TRIM(ADJUSTL(temp2))//" = "//TRIM(ADJUSTL(temp))//" GPa."
  CASE('p','P')
    msg = ">>> Applying an isostatic pressure of "//TRIM(ADJUSTL(temp))//" GPa."
  CASE DEFAULT
    msg = ">>> Applying a stress state as read from "//TRIM(ADJUSTL(temp))//"."
  END SELECT
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2125)
  !strings(1) = Cartesian axis, or integer
  !strings(2) = same type as strings(1)
  SELECT CASE(strings(1))
  CASE('x','X','y','Y','z','Z')
    msg = ">>> Swapping Cartesian axes "//TRIM(ADJUSTL(strings(1)))//" and "//TRIM(ADJUSTL(strings(2)))//"."
  CASE DEFAULT
    msg = ">>> Swapping atoms #"//TRIM(ADJUSTL(strings(1)))//" and #"//TRIM(ADJUSTL(strings(2)))//"."
  END SELECT
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2126)
  msg = "..> Swap successful."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2127)
  !strings(1) = rolled direction: x, y or z
  !strings(2) = roll axis: x, y or z
  !reals(1) = roll angle in degrees
  WRITE(msg,"(f16.2)") reals(1)
  msg = ">>> Rolling the "//TRIM(ADJUSTL(strings(1)))//" direction by " // &
      & TRIM(ADJUSTL(msg))//"° around the "//TRIM(strings(2))//" axis."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2128)
  msg = "..> System was successfully rolled."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2129)
  !strings(1) = torsion axis: x, y or z
  !reals(1) = torsion angle in degrees
  WRITE(msg,"(f16.2)") reals(1)
  msg = ">>> Applying a torsion of "//TRIM(ADJUSTL(msg)) &
      & //"° around "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2130)
  msg = "..> Torsion was successfully applied."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2131)
  !strings(1) = space group name or number
  msg = ">>> Applying symmetry operations of space group: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2132)
  !reals(1) = space group number
  WRITE(temp,"(i5)") NINT(reals(1))
  msg = "..> Space group number: "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2133)
  !reals(1) = new number of atoms
  WRITE(temp,"(i16)") NINT(reals(1))
  msg = "..> Symmetry operations were successfully applied, new number of atoms: "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2140)
  !reals(1) = min. separation distance
  !reals(2) = separate atoms by this amount
  WRITE(temp,"(f16.2)") reals(1)
  WRITE(temp2,"(f16.2)") reals(2)
  msg = ">>> Separating atoms closer than "//TRIM(ADJUSTL(temp))//" A by a distance of "//TRIM(ADJUSTL(temp2))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2141)
  !reals(1) = number of pairs of atoms that were separated
  IF( reals(1) < 1.d-3 ) THEN
    msg = "..> No atoms were separated."
  ELSEIF( NINT(reals(1)) == 1 ) THEN
    msg = "..> One pair of atoms was separated."
  ELSE
    WRITE(temp,"(i16)") NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(temp))//" pairs of atoms were separated."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2142)
  !reals(1) = number of triangles in STL file
  WRITE(temp,*) NINT(reals(1))
  msg = "..> STL file was read successfully ("//TRIM(ADJUSTL(temp))//" triangles)."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2143)
  msg = ">>> Converting system into an orthorhombic cell..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2144)
  WRITE(temp,*) NINT(reals(1))
  msg = "..> Cell is now orthorhombic ("//TRIM(ADJUSTL(temp))//" atoms)."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!2700-2799: WARNING MESSAGES
CASE(2700)
  !strings(1) = option name
  msg = "/!\ WARNING: could not understand this option: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Try `atomsk --help options` for a summary of options."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2720)
  msg = "/!\ WARNING: axis is already aligned, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2721)
  msg = "/!\ WARNING: coordinates are already reduced, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2722)
  !strings(1) = atom species
  msg = "/!\ WARNING: "//TRIM(strings(1))//" atoms already have shells."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2723)
  !strings(1) = atomic species
  msg = "/!\ WARNING: there is no "//TRIM(strings(1))//" atom in the system, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2724)
  msg = "/!\ WARNING: deformation is zero, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2725)
  msg = "/!\ WARNING: Burgers vector is zero, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2726)
  msg = "/!\ WARNING: supercell is very small in one direction normal to the"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    dislocation line! Are you sure you know what you are doing?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2727)
  !reals(1) = index of atom with large displacement
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ WARNING: displacement for atom "//TRIM(ADJUSTL(msg))// " is large."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2728)
  msg = "/!\ WARNING: expansion factors are all equal to 1, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2729)
  msg = "/!\ WARNING: no auxiliary property is defined, skipping..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2730)
  !string(1) = property
  msg = "/!\ WARNING: No "//TRIM(ADJUSTL(strings(1)))//" found in auxiliary properties, skipping..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2731)
  msg = "/!\ WARNING: Hstart = Hend, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2732)
  !strings(1) = name of unknown property
  msg = "/!\ WARNING: unknown property: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2733)
  msg = "/!\ WARNING: specified radius is negative, no atom will be removed."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2734)
  msg = "/!\ WARNING: rotation angle is modulo zero, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2735)
  msg = "/!\ WARNING: shear is zero, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2736)
  msg = "/!\ WARNING: shift vector is naught, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2737)
  msg = "/!\ WARNING: species are the same, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2738)
  msg = "/!\ WARNING: units are the same, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2739)
  msg = "/!\ WARNING: base vectors are not orthonormal."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2740)
  msg = "/!\ WARNING: elastic tensor is not symmetric!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2741)
  msg = "/!\ WARNING: Poisson ratio is out of range [-1 , 0.5]."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2742)
  !reals(1) = atom index
  WRITE(temp,*) NINT(reals(1))
  msg = "/!\ WARNING: atom index #"//TRIM(ADJUSTL(temp))//" is out of bounds, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2743)
  msg = "/!\ WARNING: skew parameters are zero, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2744)
  msg = "/!\ WARNING: there is no ionic shell in the system, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2745)
  !reals(1) = number of atoms that will actually be removed
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ WARNING: cannot select so many atoms, only "// &
      & TRIM(ADJUSTL(msg))//" will be removed."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2746)
  msg = "/!\ WARNING: nothing to be done, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2747)
  msg = "/!\ WARNING: stress intensity factor K is zero, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2748)
  msg = "/!\ WARNING: supercell is very small in one direction normal to the"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    crack line! Are you sure you know what you are doing?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2749)
  !strings(1) = name of property
  !reals(1) = atom index
  WRITE(temp,*) NINT(reals(1))
  msg = "/!\ WARNING: could not assign the value of the property '"//TRIM(ADJUSTL(strings(1)))// &
      & "' to atom #"//TRIM(ADJUSTL(temp))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2750)
  msg = "/!\ WARNING: a selection was defined but it does not contain any atom anymore."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "             Selection was cleared, all atoms are now selected."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2751)
  msg = "/!\ WARNING: target temperature is zero, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2752)
  msg = "/!\ WARNING: no selection is defined, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2753)
  !reals(1) = atom index
  WRITE(temp,*) NINT(reals(1))
  msg = "/!\ WARNING: atom #"//TRIM(ADJUSTL(temp))//" was already selected, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2754)
  !strings(1) = type of object (e.g. "dislocation" or "crack" or "loop" or "all"
  !             ("all" means that all points of a loop are out of the box)
  !reals(1) = number of points that are out of the box
  IF( strings(1)=="loop" ) THEN
    WRITE(temp,*) NINT(reals(1))
    msg = "/!\ WARNING: "//TRIM(ADJUSTL(temp))//" points of the loop are out of the box."
    CALL DISPLAY_MSG(1,msg,logfile)
  ELSEIF( strings(1)=="all" ) THEN
    msg = "/!\ WARNING: all points of the loop are out of the box!"
    CALL DISPLAY_MSG(1,msg,logfile)
    msg = "    Are you sure you know what you are doing?"
    CALL DISPLAY_MSG(1,msg,logfile)
  ELSE
    msg = "/!\ WARNING: the position of the "//TRIM(ADJUSTL(strings(1)))//" is out of the box."
    CALL DISPLAY_MSG(1,msg,logfile)
    msg = "    Are you sure you know what you are doing?"
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(2755)
  msg = "/!\ WARNING: the factor is zero, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2756)
  !strings(1) = direction
  msg = "/!\ WARNING: supercell is quite large along the "//TRIM(ADJUSTL(strings(1)))//" direction!"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Are you sure you know what you are doing?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2757)
  msg = "/!\ WARNING: the indices are the same, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2758)
  msg = "/!\ WARNING: no operation to apply, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2759)
  msg = "/!\ WARNING: dislocation loop radius is too small, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2760)
  msg = "/!\ WARNING: cell is already orthorhombic, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
  !
CASE(2799)
  !strings(1) = name of obsolete option
  !strings(2) = name of new option
  msg = "/!\ WARNING: option "//TRIM(ADJUSTL(strings(1)))//" is obsolete and will be removed."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Please use option "//TRIM(ADJUSTL(strings(2)))//" instead."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!2800-2899: ERROR MESSAGES
CASE(2800)
  !string(1) = axis (if we are here it"s because it is different from x, y or z)
  msg = "X!X ERROR: unknown axis: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2801)
  WRITE(msg,*) "X!X ERROR: the base Hend is not a rotation of Hstart."
  CALL DISPLAY_MSG(1,msg,logfile)
  WRITE(msg,*) "    Check that angles are equal in Hstart and Hend."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2802)
  !strings(1) = property that was not read properly
  msg = "X!X ERROR while reading "//TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2803)
  msg = "X!X ERROR: there were errors while trying to determine supercell vectors."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2804)
  msg = "X!X ERROR: there were errors while applying options."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2805)
  !strings(1) = name of unknown option
  msg = "X!X ERROR: unknown option: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2806)
  !strings(1) = name of option
  msg = "X!X ERROR: non-conform statement in option: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2807)
  !reals(1) = 1 if roots of Eq.(13-85) cannot be found
  !         = 2 if the A_k(n) cannot be calculated
  !         = 3 if the linear equations defining D(n) cannot be solved
  msg = "X!X ERROR:"
  IF(NINT(reals(1))==1) THEN
    msg = TRIM(msg)//" unable to determine the P(n), aborting."
  ELSEIF(NINT(reals(1))==2) THEN
    msg = TRIM(msg)//" unable to determine the A_k(n), aborting."
  ELSEIF(NINT(reals(1))==3) THEN
    msg = TRIM(msg)//" unable to determine the D(n), aborting."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2808)
  msg = "X!X ERROR: cannot build a mixed dislocation with isotropic elasticity."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2809)
  msg = "X!X ERROR: the elastic tensor contains NaN values, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2810)
  msg = "X!X ERROR: inconsistent array size for shells."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2811)
  msg = "X!X ERROR: dislocline and dislocplane must be normal to each other, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2812)
  msg = "X!X ERROR: unable to determine what atom(s) to remove, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2813)
  !strings(1) = string that could not be converted to a number
  msg = "X!X ERROR: unable to convert this string into a number: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2814)
  msg = "X!X ERROR: there is no atomic system to apply this option to."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2815)
  !strings(1) = name of matrix to invert
  msg = "X!X ERROR: unable to invert the matrix "//strings(1)//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2816)
  msg = "X!X ERROR: the elastic tensor is not defined, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2817)
  msg = "X!X ERROR: the property '"//TRIM(ADJUSTL(strings(1)))//"' is not defined, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2818)
  msg = "X!X ERROR: there was an error while reading the STL file, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2819)
  msg = "X!X ERROR: unable to find an orthogonal cell from initial cell vectors."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!
!
!**************************
! 3000-3999: MESSAGES FOR OUTPUT
CASE(3000)
  !reals(1) = number of atoms
  !reals(2) = number of shells
  WRITE(temp,*) NINT(reals(1))
  IF( NINT(reals(2))>0 ) THEN
    WRITE(temp2,*) NINT(reals(2))
    msg = ">>> Writing output file(s) ("//TRIM(ADJUSTL(temp))//" cores + "// &
        & TRIM(ADJUSTL(temp2))//" shells):"
  ELSE
    msg = ">>> Writing output file(s) ("//TRIM(ADJUSTL(temp))//" atoms):"
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3001)
  !strings(1) = name of file
  msg = "..> This file already exists and was not converted again: "// &
      & TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3002)
  !strings(1) = type of file, e.g. "XSF" or "CFG"
  !strings(2) = name of file
  msg = "..> Successfully wrote "//TRIM(strings(1))//" file: " &
      & //TRIM(strings(2))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3003)
  !reals(:) = list of atomic numbers
  msg = "..> Atom species were packed:"
  DO i=1,SIZE(reals(:))
    CALL ATOMSPECIES(reals(i),species)
    msg = TRIM(msg)//" "//TRIM(species)//", "
  ENDDO
  j=LEN_TRIM(msg)
  msg = TRIM(msg(:j-1))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "    Check that this is consistent with the POTCAR file."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3004)
  !strings(1) = skew parameter
  msg = "..> Box skew "//TRIM(strings(1))//" was reduced."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3005)
  msg = "..> OK, leaving skew parameter "//TRIM(strings(1))//" untouched."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3006)
  !strings(1) = name of ddplot file
  msg = ">>> Building the ddplot file: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!3700-3799: WARNING MESSAGES
CASE(3700)
  msg = "/!\ WARNING: no output file name was specified, please provide one:"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3701)
  !strings(1) = name of output file
  msg = "/!\ WARNING: I could not guess the format for the output file: " &
      & //TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "<?> Please specify a format for the output file:"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = '    ('//TRIM(listofformats(1))
  DO i=2,SIZE(listofformats)
    msg = TRIM(msg)//','//TRIM(listofformats(i))
  ENDDO
  msg = TRIM(msg)//')'
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3702)
  msg = "/!\ WARNING: ddplot format is available only when using DDPLOT mode."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            No dd file will be output, please refer to the documentation."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3703)
  msg = "/!\ WARNING: atom species are not contiguous. Do you want to pack them? ("&
      & //langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    (this will affect only the POSCAR file)"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3704)
  msg = "/!\ WARNING: supercell does not form a lower-triangular matrix, which is"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    required by LAMMPS. Do you want to re-align the system? (" &
      & //langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    (This will affect only the LAMMPS output file)"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3705)
  !strings(1) = skew parameter
  msg = "/!\ WARNING: triclinic box skew "//TRIM(strings(1))//" is too large."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Running LAMMPS with such a supercell may result in the same"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    error message. Do you want to reduce the skew? ("&
      & //langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    (This will affect only the LAMMPS output file)"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3706)
  !strings(1) = skew parameter
  msg = "/!\ WARNING: unable to reduce box skew "//TRIM(strings(1))//", aborting..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3707)
  !reals(1) = number of atoms
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ WARNING: number of particles if very large: "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    If ddplot cannot display so many atoms, then use Atomsk"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    with the -cut option to reduce the number of atoms."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3708)
  msg = "/!\ WARNING: only the first 32 auxiliary properties will be "// &
      & "written to the CFG file."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3709)
  msg = "/!\ WARNING: some supercell parameters cannot be written and"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    this data file WILL BE WRONG!!!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3710)
  !strings(1) = file format
  msg = "/!\ WARNING: unknown format '"//TRIM(strings(1))//"', skipping..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3711) ! missing occupancy data
  msg = "/!\ WARNING: occupancy data is missing, occupancy of all atoms will be set to 1."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3712) ! missing thermal vibration data
  msg = "/!\ WARNING: Debye-Waller factors are missing,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            thermal vibration will be set to zero for all atoms."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3713) ! missing absorption data
  msg = "/!\ WARNING: absorption factors are missing, they will be set to 0.03 for all atoms."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3714)
  msg = "/!\ WARNING: some atoms have an invalid type."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3715)
  msg = "/!\ WARNING: input data contains partial occupancies, which are"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            not supported by some output format(s)."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Some atoms may overlap in the output file(s), which is not physical."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3716)
  msg = "/!\ WARNING: data contains ionic shells, which are"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            not supported by some output format(s)."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Shells will be lost in some output file(s)."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3717)
  !reals(1) = total charge
  WRITE(temp,'(f9.3)') reals(1)
  msg = "/!\ WARNING: cell has a non-zero electric charge: Q_total = "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(1,msg,logfile)
!
!3800-3899: ERROR MESSAGES
CASE(3800)
  msg = "X!X ERROR: no atom position to write, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3801)
  !strings(1) = name of file
  msg = "X!X ERROR: there were errors while writing the file: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3802)
  msg = "X!X ERROR: there were errors while writing files."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3803)
  !reals(1) = index of atom that has a NaN coordinate
  WRITE(msg,*) NINT(reals(1))
  msg = "X!X ERROR: atom #"//TRIM(ADJUSTL(msg))//" has NaN coordinate, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3804)
  !reals(1) = index of atom that has a NaN coordinate
  WRITE(msg,*) NINT(reals(1))
  msg = "X!X ERROR: auxiliary property of atom #"//&
      & TRIM(ADJUSTL(msg))//" is NaN, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!
!
!**************************
! 4000-4999: MESSAGES FOR MODES
CASE(4000)
  !strings(1) = name of file for mode "list"
  msg = ">>> Reading the file list from: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4001)
  !Just a separator
  msg = "         --===============--"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4002)
  !strings(1) = output format activated, e.g. "xyz" or "cfg"
  msg = ">>> Converting all following files to: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4010)
  msg = ">>> Using ddplot mode."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4020)
  !strings(1) = file containing several snapshots
  msg = ">>> Unfolding "//TRIM(strings(1))//" to many files..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4021)
  msg = ">>> Atomsk is a free, Open Source software."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "    To learn more, enter 'license'."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4022)
  msg = ">>> Atomsk command-line interpreter:"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = '..> Type "help" for a summary of commands.'
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = ""
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4023)
  msg = "Atomsk is currently running in INTERACTIVE MODE."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "Only the following commands are available:"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "help              Display this help"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "cd                Change working directory"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = system_ls//"               List files in current directory"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "pwd               Print current working directory"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "memory            Summary of what is in memory"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "create            Create an atomic system"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "read <file>       Read the <file> and load its content in memory"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "write <file>      Write current system into <file>"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "clear             Clear memory (destroy atomic system)"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "quit              Exit Atomsk"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "OPTIONS: the options of Atomsk can be used in this command-line interpreter,"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "        enter 'help options' to display the available options."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "        In interactive mode options must be called without the leading minus sign (-)."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "MODES: modes cannot be used in this command-line interpreter."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4024)
  msg = "<?> To which format do you want to convert it?"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = '    ('//TRIM(listofformats(1))
  DO i=2,SIZE(listofformats)
    msg = TRIM(msg)//','//TRIM(listofformats(i))
  ENDDO
  msg = TRIM(msg)//')'
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4025)
  msg = ">>> Oops, I could not convert your file, I am truly sorry..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "    I hope we are still friends! :-)"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4026)
  msg = ">>> The file was successfully converted."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4027)
  msg = ">>> Creating system:"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4028)
  !strings(1) = description of system that is created
  msg = "..> "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4029)
  msg = "..> System was successfully created."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4030)
  !strings(1) = merge direction
  !reals(1) >0 if systems are concatenated, <0 otherwise
  IF( reals(1)>0.d0 ) THEN
    msg = ">>> Concatenating the systems along "//TRIM(strings(1))//"..."
  ELSE
    msg = ">>> Merging the systems..."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4031)
  !reals(1) = number of files merged
  WRITE(msg,*) NINT(reals(1))
  msg = "..> "//TRIM(ADJUSTL(msg))//" systems were merged."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4032)
  msg = ">>> Computing electric dipole moments:"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4033)
  !strings(1) = name of file where total polarization is written
  msg = "..> Wrote the total polarization: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4034)
  !strings(1) = first atomic species
  !strings(2) = second atomic species
  !reals(1) = coordination
  IF(INT(reals(1))==4) THEN
    msg = " tetrahedra."
  ELSEIF(INT(reals(1))==6) THEN
    msg = " octahedra."
  ELSEIF(INT(reals(1))==8) THEN
    msg = " hexahedra."
  ELSE
    msg = " polyhedra."
  ENDIF
  msg = ">>> Computing the moments of "//TRIM(strings(1))// &
      & "-"//TRIM(strings(2))//TRIM(msg)
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  IF(reals(1)<0.d0) THEN
    WRITE(msg,"(f16.2)") DABS(reals(1))
    msg = "..> Using neighbors in a radius of "//TRIM(ADJUSTL(msg))// &
        & " A around "//TRIM(strings(1))//" ions."
  ELSEIF(reals(1)==0.d0) THEN
    msg = "..> Trying to find nearest neighbors automatically."
  ELSE
    WRITE(msg,*) INT(reals(1))
    msg = "..> Using the "//TRIM(ADJUSTL(msg))//" first neighbors of "// &
        & TRIM(strings(1))//" ions."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4035)
  !strings(1) = file name
  msg = "..> Successfully wrote the XSF file: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4036)
  !strings(1) = file name
  msg = "..> Successfully wrote the norms: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4037)
  !strings(1) = file name
  msg = "..> Successfully wrote the statistics: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4038)
  msg = ">>> Computing the difference..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4039)
  !strings(1) = name of file
  msg = "..> Successfully wrote the file: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4040)
  msg = "..> Differencial displacements were computed."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4041)
  !reals(1) = snapshots number
  WRITE(msg,*) NINT(reals(1))
  msg = ">>> Reading snapshot #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4042)
  !reals(1) = number of snapshots read
  WRITE(msg,*) NINT(reals(1))
  IF( NINT(reals(1)) <= 0 ) THEN
    msg = ">>> No snapshot was converted."
  ELSEIF( NINT(reals(1)) == 1 ) THEN
    msg = ">>> 1 snapshot was successfully converted."
  ELSE
    msg = ">>> "//TRIM(ADJUSTL(msg))//" snapshots were successfully converted."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4043)
  msg = "..> Snapshot was read successfully."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4044)
  msg = ">>> Writing properties for all atoms:"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4045)
  !reals(1) = number of files converted
  !reals(1) = number of files ignored
  WRITE(msg,*) NINT(reals(1))
  IF( NINT(reals(1))==0 ) THEN
    msg = ">>> No file was converted"
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = ">>> 1 file was converted"
  ELSE
    msg = ">>> "//TRIM(ADJUSTL(msg))//" files were converted"
  ENDIF
  IF( SIZE(reals)>1 .AND. NINT(reals(2))>0 ) THEN
    WRITE(temp,*) NINT(reals(2))
    IF( NINT(reals(2))==1 ) THEN
      msg = TRIM(ADJUSTL(msg))//", 1 file was ignored"
    ELSE
      msg = TRIM(ADJUSTL(msg))//", "//TRIM(ADJUSTL(temp))//" files were ignored"
    ENDIF
  ENDIF
  msg = TRIM(ADJUSTL(msg))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4046)
  msg = ">>> Unwrapping atoms..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4047)
  !reals(1) = number of atoms unwrapped
  IF( NINT(reals(1))==0 ) THEN
    msg = "..> No atom was unwrapped."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 atom was unwrapped."
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" atoms were unwrapped."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4048)
  msg = ">>> Computing the electronic polarization of all ions in the system..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4049)
  msg = "..> Electronic polarization was computed."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4050)
  !strings(1) = file containing the names of files to include
  msg = ">>> Folding all files listed in "//TRIM(strings(1))// &
      & " into one file..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4051)
  !reals(1) = width of the skin (A)
  msg = ">>> Computing the radial distribution function,"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,'(f16.3)') reals(1)
  msg = "    using a skin width of "//TRIM(ADJUSTL(msg))//" A."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4052)
  !strings(1) = species 1
  !strings(2) = species 2
  msg = "..> Computing RDF of "//TRIM(strings(2))//" atoms around "// &
      & TRIM(strings(1))//" atoms..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4053)
  msg = "..> RDFs were computed."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4054)
  msg = ">>> Constructing a polycrystal using Voronoi method."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4055)
  !reals(1) = index of the grain
  WRITE(msg,*) NINT(reals(1))
  msg = ">>> Generating grain #"//TRIM(ADJUSTL(msg))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4056)
  !reals(1) = number of atoms in the grain
  WRITE(msg,*) NINT(reals(1))
  msg = "..> Number of atoms in this grain: "//TRIM(ADJUSTL(msg))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4057)
  !strings(1) = name of input file
  msg = ">>> Reading parameters for Voronoi construction from: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4058)
  !reals(1) = number of grains
  !reals(2) = 0 if 3-D, 1,2,3 if thin along x, y, z
  msg = "..> File was successfully read."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,*) NINT(reals(1))
  msg = "..> Number of grains to be constructed: "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  IF( NINT(reals(2))==0 ) THEN
    msg = "..> Using a 3-D Voronoi tesselation."
  ELSE
    IF( NINT(reals(2))==1 ) THEN
      msg = "x"
    ELSEIF( NINT(reals(2))==2 ) THEN
      msg = "y"
    ELSEIF( NINT(reals(2))==3 ) THEN
      msg = "z"
    ENDIF
    msg = "..> Using a 2-D Voronoi tesselation, rotation axis: "//TRIM(ADJUSTL(msg))
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4059)
  msg = ">>> Constructing a chain of interpolated configurations."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4060)
  !reals(1) = index of the image
  WRITE(msg,*) NINT(reals(1))
  msg = ">>> Constructing configuration #"//TRIM(ADJUSTL(msg))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4061)
  msg = ">>> Computation of the Nye tensor."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4062)
  msg = ">>> Computing per-atom G matrix..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4063)
  msg = ">>> Computing per-atom Nye tensor..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4064)
  !strings(1) = name of the file containing file list
  msg = ">>> Averaging atom positions for the list of files: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4065)
  !reals(1) = number of files
  WRITE(temp,*) NINT(reals(1))
  msg = "..> Atom positions were averaged over "//TRIM(ADJUSTL(temp))//" configurations."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4066)
  msg = ">>> Running in mode density."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4067)
  !strings(1) = property
  !reals(1) = dimension (1, 2 or 3)
  WRITE(temp,*) NINT(reals(1))
  msg = ">>> Computing the "//TRIM(ADJUSTL(temp))//"-D density of "//TRIM(ADJUSTL(strings(1)))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4068)
  WRITE(msg,*) NINT(reals(1))
  msg = "..> Density was successfully computed."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4069)
  !strings(1) = name of file
  msg = ">>> Computing the central symmetry parameter for: "//TRIM(ADJUSTL(strings(1)))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4200)
  WRITE(*,*) " (type q to cancel)"
  WRITE(*,'(a39)',ADVANCE='NO') " Lattice type (sc,bcc,fcc,dia,rs,per): "
CASE(4201)
  !strings(1) = name of lattice parameter (a, b, c, a0...)
  WRITE(*,'(a28)',ADVANCE='NO') " Lattice parameter "//TRIM(ADJUSTL(strings(1)))//" (Å): "
CASE(4202)
  WRITE(*,'(a15)',ADVANCE='NO') " Atom species: "
CASE(4203)
  WRITE(*,'(a22)',ADVANCE='NO') " Chiral indices (n,m): "
CASE(4204)
  WRITE(*,'(a21)',ADVANCE='NO') " Lattice orientation: "
CASE(4300)
  !reals(1) = number of max tries
  msg = "   ***   ATOMSK QUIZZ   ***"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(temp,*) NINT(reals(1))
  msg = " Answer the following question, you have "//TRIM(ADJUSTL(temp))//" tries."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4301)
  !strings(1) = answer
  !reals(1) = type of quizz (1=guess atomic number; 2=guess atom species; 3=guess next atom)
  !reals(2) = number of tries
  !reals(3) = user's answer - actual solution (quizz type 1 only)
  msg=""
  IF( NINT(reals(2))>0 ) THEN
    msg = " WRONG!"
  ENDIF
  IF( NINT(reals(1))==1 ) THEN
    IF( NINT(reals(2))>0 ) THEN
      IF( reals(3) > 0.d0 ) THEN
        msg = TRIM(msg)//" It's less."
      ELSEIF( reals(3) < 0.d0 ) THEN
        msg = TRIM(msg)//" It's more."
      ENDIF
    ENDIF
    msg = TRIM(msg)//" What is the atomic number of the atom: "//TRIM(ADJUSTL(strings(1)))//"?"
  ELSEIF( NINT(reals(1))==2 ) THEN
    msg = TRIM(msg)//" Which atom has the atomic number "//TRIM(ADJUSTL(strings(1)))//"?"
  ELSE
    msg = TRIM(msg)//" Which atom comes after "//TRIM(ADJUSTL(strings(1)))//" in the periodic table?"
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4302)
  !reals(1) = try/maxtry
  IF( reals(1)>=1.d0 ) THEN
    msg = " GAME OVER!"
  ELSE
    msg = " CORRECT!"
  ENDIF
  msg = TRIM(msg)//" The answer was "//TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!4700-4799: WARNING messages
CASE(4700)
  !strings(1) = name of file that does not exist
  msg = "/!\ WARNING: This file does not exist: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Skipping..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4701)
  msg = "/!\ I did not find the cell vectors in the input file."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    I will try to guess them, but this is very inaccurate."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    I advise that you use the option -prop to setup the supercell vectors."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4702)
  msg = "/!\ WARNING: cell is charged, total polarization may be wrong!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4703)
  msg = "/!\ WARNING: species has zero charge, skipping..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4705)
  !reals(1) = index of atom that has too many neighbors
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ WARNING: number of neighbors exceeds 100 for atom #" &
      & //TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4706)
  !strings(1) = atomic species of absent neighbors
  !reals(1) = index of atom that has no neighbour
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ WARNING: I did not find any "//TRIM(strings(1))// &
      &       " neighbor for atom #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4707)
  !reals(1) = index of atom that has a shell with zero charge
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ WARNING: the shell of ion #"//TRIM(ADJUSTL(msg)) &
      & //" has a zero charge."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4708)
  msg = "/!\ WARNING: this grain contains zero atom."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4709)
  !reals(1) = index of atom
  !reals(2) = number of neighbors
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = "/!\ WARNING: insufficient number of neighbors ("//TRIM(ADJUSTL(temp2))// &
      &              ") for atom #"//TRIM(ADJUSTL(temp))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4710)
  !strings(1) = name of the matrix
  msg = "/!\ WARNING: could not compute the matrix "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4711)
  !strings(1) = name of the file
  msg = "/!\ WARNING: this file has a different number of atoms and will not be treated: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4712)
  msg = "/!\ WARNING: it is recommended to run this mode with the option '-wrap',"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            otherwise atoms that are out of the box may lead to wrong results."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Do you wish to wrap them now? ("//langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4713)
  !strings(1) = action to follow if user says "yes"
  IF(LEN_TRIM(strings(1))<=0) THEN
    strings(1) = "proceed"
  ENDIF
  msg = "/!\ WARNING: apparently you did not save the system into a file."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Are you sure you want to "//TRIM(ADJUSTL(strings(1)))//"? ("//langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4714)
  !strings(1) = direction along which cell is small
  !reals(1) = new cell size along that direction
  WRITE(temp,'(f16.3)') reals(1)
  msg = "/!\ WARNING: final cell has a small dimension along " &
      & //TRIM(ADJUSTL(strings(1)))//", setting it to "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(1,msg,logfile)
!
!4800-4899: ERROR MESSAGES
CASE(4800)
  !strings(1) = mode
  msg = "X!X ERROR: non-conform statement in mode: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4801)
  !strings(1) = file format (e.g. "xyz", "cfg"...)
  msg = "X!X ERROR: this file format is not yet supported in mode 1-in-all:" &
      & //TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4802)
  msg = "X!X ERROR: both chiral indices cannot be zero, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4803)
  !reals(1) = theoretical number of atoms in nanotube
  !reals(2) = number found by Atomsk
  WRITE(msg,*) "X!X ERROR: inconsistent number of atoms in nanotube."
  CALL DISPLAY_MSG(1,msg,logfile)
  WRITE(msg,*) NINT(reals(1))
  WRITE(temp,*) NINT(reals(2))
  msg = "          Theory: "//TRIM(ADJUSTL(msg))//"; Found: "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4804)
  !reals(1) = number of species required
  !reals(2) = alternative number of species required (optional)
  !reals(3) = alternative number of species required (optional)
  WRITE(temp,*) NINT(reals(1))
  IF(SIZE(reals)>1) THEN
    WRITE(temp2,*) NINT(reals(2))
    IF(SIZE(reals)>2) THEN
      WRITE(temp3,*) NINT(reals(3))
      temp = TRIM(ADJUSTL(temp))//", "//TRIM(ADJUSTL(temp2))//" or "//TRIM(ADJUSTL(temp3))
    ELSE
      temp = TRIM(ADJUSTL(temp))//" or "//TRIM(ADJUSTL(temp2))
    ENDIF
  ENDIF
  msg = "X!X ERROR: this structure requires "//TRIM(ADJUSTL(temp))//" atom species."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4805)
  !strings(1) = structure type
  msg = "X!X ERROR: unknown structure: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4806)
  msg = "X!X ERROR: no file to merge, aborting"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4807)
  msg = "X!X ERROR: all charges are zero, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Charges must be specified with the option -properties."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4808)
  !strings(1) = atomic species that has zero charge
  msg = "X!X ERROR: "//TRIM(strings(1))//" ions cannot have zero charge."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4809)
  !strings(1) = atomic species
  msg = "X!X ERROR: no "//TRIM(strings(1))//" ion in the system."
  msg = TRIM(ADJUSTL(msg))
CASE(4810)
  !reals(1) = number of particles in first system
  !reals(2) = number of particles in second system
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = "X!X ERROR: number of particles of the two systems differ: "// &
      & TRIM(ADJUSTL(temp))//"/"//TRIM(ADJUSTL(temp2))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "           Program will EXIT."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4811)
  msg = "X!X ERROR: the two systems do not have the same base vectors."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4812)
  msg = "X!X ERROR: file does not appear to have animated XSF format, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4813)
  !strings(1) = name of unknown mode
  msg = "X!X ERROR: unknown mode: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4814)
  msg = "X!X ERROR: only one mode can be used at a time, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4815)
  msg = "X!X ERROR: file does not appear to have DL_POLY HISTORY format, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4816)
  !reals(1) = index for core charges
  !reals(2) = index for shell charges
  IF( reals(1)<reals(2) ) THEN
    msg = "X!X ERROR: the charges of ionic cores are not defined, aborting."
  ELSE
    msg = "X!X ERROR: the charges of ionic shells are not defined, aborting."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4817)
  msg = "X!X ERROR: there is no ionic shell in the system, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4818)
  !strings(1) = file format (e.g. "xyz", "cfg"...)
  msg = "X!X ERROR: the file "//TRIM(strings(1))//" seems to be empty,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          or does not contain any valid file name."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4819)
  !strings(1) = suggestion for a vector
  !strings(2) = direction of suggested vector (X,Y or Z)
  msg = "X!X ERROR: base vectors are not orthogonal."
  CALL DISPLAY_MSG(1,msg,logfile)
  IF( LEN_TRIM(strings(1))>0 ) THEN
    msg = "    Suggested vector along "//TRIM(ADJUSTL(strings(2)))//": "//TRIM(ADJUSTL(strings(1)))
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(4820)
  msg = "X!X ERROR: supercell dimensions were not defined, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4821)
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = "X!X ERROR: number of atoms ("//TRIM(ADJUSTL(temp))// &
      & ") exceeds size of allocated array ("//TRIM(ADJUSTL(temp2))//")."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4822)
  msg = "X!X ERROR: no file to be treated, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4823)
  !strings(1) = keyword 1
  !strings(2) = keyword 2
  msg = "X!X ERROR: keywords '"//TRIM(ADJUSTL(strings(1)))//"' and '" &
      &  //TRIM(ADJUSTL(strings(2)))//"' are mutually exclusive, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4824)
  !strings(1) = name of unknown command
  msg = "X!X ERROR: unknown command: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Type 'help' for a list of available commands."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4825)
  msg = "X!X ERROR: Atomsk cannot run within itself!"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          You have entered 'atomsk' without any arguments, or you have clicked the"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          executable from a menu, therefore Atomsk is currently running in interactive mode"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          where only a limited subset of commands are available, please refer to the documentation."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          In order to use Atomsk with command-line arguments, you must first run"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          a command shell, and then type in your command."
  CALL DISPLAY_MSG(1,msg,logfile)
#if defined(WINDOWS)
  msg = "          In a Microsoft Windows system, follow these steps:"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          1. Open the Windows menu, go to Accessories, and run Windows Power Shell."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          2. Run atomsk.exe with the arguments, for instance:"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = '             & "C:\Program Files\Atomsk\atomsk.exe" initial.xsf final.cfg'
  CALL DISPLAY_MSG(1,msg,logfile)
#endif
CASE(4826)
  msg = "X!X ERROR: this mode is not available in interactive mode, please refer to the documentation."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4827)
  msg = "X!X ERROR: unable to create lattice with specified orientation."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4828)
  msg = "X!X ERROR: at least two cell dimensions are too small, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4829)
  msg = "X!X ERROR: box vectors are not linearly independent, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4900)
  msg = "X!X ERROR: only one mode can be used at a time."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!
!
!Special case for mode PREFERENCES: many questions can be added there
CASE(5000)
  !strings(1) = user's configuration file
  msg = ">>> The answers to the following questions will set up your preferences"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "..> for Atomsk, which will be saved in "//TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5001)
  msg = "<?> Do you want Atomsk to always overwrite files by default?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5002)
  msg = "<?> Do you want Atomsk to always ignore existing files by default?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5003)
  msg = "<?> Enter a file format you always want to write to ('no' to skip):"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5004)
  msg = "<?> Enter the default language (NO to skip):"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5005)
  msg = "<?> Enter the default level of verbosity (default=1):"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    0=silent; no message is ever displayed nor written in log file"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    1=all messages are written on the screen only (default)."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    2=all messages are written in log file only."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    3=messages are written in log+screen."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    4=debug, additional messages are written in log file."
  CALL DISPLAY_MSG(1,msg,logfile)
!
CASE(5700)
  !strings(1) = user's configuration file
  msg = "/!\ WARNING: it is necessary to overwite the existing file " &
      & //TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Are you sure you want to proceed? (yes/no)"
  CALL DISPLAY_MSG(1,msg,logfile)
!
!
!
CASE DEFAULT
  WRITE(msg,*) imsg
  msg = "@@@ UNKNOWN MESSAGE #"//TRIM(ADJUSTL(msg))//" (this is a bug) @@@"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  !
  !Other languages: in this section "CASE DEFAULT", comment the message above
  !                and uncomment the line below
  !CALL ATOMSK_MSG_EN(imsg,strings,reals)
END SELECT
!
END SUBROUTINE ATOMSK_MSG_EN
!
!
!
!********************************************************
! DATE_MSG
! Displays a nice message on certain dates.
!********************************************************
SUBROUTINE DATE_MSG_EN()
!
IMPLICIT NONE
CHARACTER(LEN=5):: zone
CHARACTER(LEN=8):: date
CHARACTER(LEN=10):: time
CHARACTER(LEN=128):: msg
INTEGER,DIMENSION(8):: values
REAL(dp),DIMENSION(:),ALLOCATABLE:: randarray    !random numbers
!
!Get the current date
CALL DATE_AND_TIME(DATE, TIME, ZONE, VALUES)
!VALUES(1): 	The year
!VALUES(2): 	The month
!VALUES(3): 	The day of the month
!VALUES(4): 	Time difference with UTC in minutes
!VALUES(5): 	The hour of the day
!VALUES(6): 	The minutes of the hour
!VALUES(7): 	The seconds of the minute
!VALUES(8): 	The milliseconds of the second 
!
!New year
IF(values(2)==1 .AND. values(3)<=10) THEN
  WRITE(msg,*) values(1)
  msg = TRIM(ADJUSTL(msg))
  msg = "*** HAPPY NEW YEAR "//TRIM(msg)//"!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!14 March (3/14): Pi day
ELSEIF(values(2)==3 .AND. values(3)==14) THEN
  WRITE(msg,'(a8,f51.48)') "    π = ", pi
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!31 March
ELSEIF(values(2)==3 .AND. values(3)==31) THEN
  msg = "  WORLD BACKUP DAY - REMEMBER TO BACKUP YOUR DATA!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!1 April
ELSEIF(values(2)==4 .AND. values(3)==1) THEN
  msg = "              \,-^--._"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "              /`--;-` "
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!1 May
ELSEIF(values(2)==5 .AND. values(3)==1) THEN
  msg = "*** You working? Workers Day is a HOLIDAY! :-)"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!4 May: Star Wars day ("May the fourth/May the force be with you")
ELSEIF(values(2)==5 .AND. values(3)==4) THEN
  msg = "               .=."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "              '==c|"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "              [)-+|"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "              //'_|"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "             /]==;\"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!14 July
ELSEIF(values(2)==14 .AND. values(3)==07) THEN
  msg = "*** Today is France National Day!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!8 August: international cat day
ELSEIF(values(2)==8 .AND. values(3)==8) THEN
  msg = "             /'---'\"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "            ( o   o )"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "            = ._Y_. ="
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!31 October: Halloween
ELSEIF(values(2)==10 .AND. values(3)==31) THEN
  CALL GEN_NRANDNUMBERS(1,randarray)
  IF( randarray(1)<0.3333d0 ) THEN
    msg = "                ,"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "           ,---'---,"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "          /  0   0  \"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "         |     A     |"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "          \ ^~~^~~^ /"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( randarray(1)<0.6666d0 ) THEN
    msg = "               _/\_"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "               ((°>"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "          _    /^|"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "          =>--/__|m--"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "               ^^"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSE
    msg = "_______________"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = " \ \/_|_\/ / /("
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "\ \/__|__\/ / )"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = " \/___|___\/  ("
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = " /    |    \  )"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "              (_"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "            _\(_)/_"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "             /(o)\"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
!
!Christmas
ELSEIF(values(2)==12 .AND. values(3)>=20 .AND. values(3)<=25) THEN
  msg = "                <>"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "                /\"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "               /  \"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "               / °\"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "              /°   \"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "             /   ° `\"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "             / `   ,\"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "            /___° ___\"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "                ||"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = " * * *   MERRY CHRISTMAS!   * * *"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!If it's late at night or in the week end
ELSEIF( values(5)>=20 .OR. values(5)<=6 ) THEN
  msg = "*** Working out of office hours? You should sleep sometimes. :-)"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!If it's lunch time
ELSEIF(values(5)==12 .AND. values(6)>=30) THEN
  msg = "*** I am hungry! :-p"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!13:37, it's leet time :-)
ELSEIF(values(5)==13 .AND. values(6)==37) THEN
  msg = "*** 17'5 1337 71|V|3!!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
ENDIF
!
END SUBROUTINE DATE_MSG_EN
!
!
!
END MODULE messages_en
