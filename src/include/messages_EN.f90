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
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 27 Oct. 2023                                     *
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
USE random
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
IF(helpsection=="modes" .OR. helpsection=="gather"  .OR. helpsection=="ai1" .OR. helpsection=="all-in-one") THEN
  WRITE(*,*) "..> Gather mode:"
  WRITE(*,*) "          atomsk --gather <listfile> <outputfile> [options] [<formats>]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="unfold" .OR. helpsection=="1ia" .OR. helpsection=="one-in-all") THEN
  WRITE(*,*) "..> Unfold mode:"
  WRITE(*,*) "          atomsk --unfold <inputfile> [<outputfile>] [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="create") THEN
  WRITE(*,*) "..> Create mode:"
  WRITE(*,*) "          atomsk --create <structure> <a> [<c>] <species> <outputfile> [orient hkl hkl hkl] [<formats>] [options]"
  WRITE(*,*) "                            <structure> | N.lattice cst. | N.at.sp."
  WRITE(*,*) "                           -------------+----------------+----------"
  WRITE(*,*) "             CUBIC                  sc  |       1        |     1"
  WRITE(*,*) "             LATTICES              bcc  |       1        |   1 or 2"
  WRITE(*,*) "                                   fcc  |       1        |   1 or 2"
  WRITE(*,*) "                               diamond  |       1        |   1 or 2"
  WRITE(*,*) "                                  L1_2  |       1        |     2"
  WRITE(*,*) "                              fluorite  |       1        |     2"
  WRITE(*,*) "                             rock-salt  |       1        |     2"
  WRITE(*,*) "                            perovskite  |       1        |     3"
  WRITE(*,*) "                                   A15  |       1        |     2"
  WRITE(*,*) "                                   C15  |       1        |     2"
  WRITE(*,*) "                           -------------+----------------+----------"
  WRITE(*,*) "             TETRAGONAL             st  |  2 (a and c)   |   1 or 2"
  WRITE(*,*) "             LATTICES              bct  |  2 (a and c)   |   1 or 2"
  WRITE(*,*) "                                   fct  |  2 (a and c)   |   1 or 2"
  WRITE(*,*) "                           -------------+----------------+----------"
  WRITE(*,*) "             HEXAGONAL             hcp  |  2 (a and c)   |   1 or 2"
  WRITE(*,*) "             LATTICES         wurtzite  |  2 (a and c)   |     2"
  WRITE(*,*) "                              graphite  |  2 (a and c)   |   1 or 2"
  WRITE(*,*) "                                    BN  |  2 (a and c)   |     2"
  WRITE(*,*) "                                   C14  |  2 (a and c)   |     2"
  WRITE(*,*) "                                   C36  |  2 (a and c)   |     2"
  WRITE(*,*) "          atomsk --create nanotube <a> <m> <n> <sp1> [<sp2>] [options] <outputfile> [<formats>]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="cpprop") THEN
  WRITE(*,*) "..> Mode copy properties:"
  WRITE(*,*) "          atomsk --copy-properties <file1> <file2> [options] <outputfile> [<formats>]"
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
  WRITE(*,*) "          atomsk --density <file> <property> <1d|2d|3d> <x|y|z> <sigma> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="difference") THEN
  WRITE(*,*) "..> Difference mode:"
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
  WRITE(*,*) "          atomsk --nye <reference> <defective> <outputfile> [options] [<formats>]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="interpolate") THEN
  WRITE(*,*) "..> Mode Interpolate:"
  WRITE(*,*) "          atomsk --interpolate <file1> <file2> <N> [options] [<formats>]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="cs") THEN
  WRITE(*,*) "..> Mode Local Symmetry:"
  WRITE(*,*) "          atomsk --local-symmetry <file> [options] [<formats>]"
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
IF(helpsection=="options" .OR. helpsection=="-cell") THEN
  WRITE(*,*) "..> Modify the cell vectors:"
  WRITE(*,*) "          -cell <add|rm|set> <length>  <H1|H2|H3|x|y|z|xy|xz|yx|yz|zx|zy|xyz>"
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
  WRITE(*,*) "..> Apply uniaxial or shear strain:"
  WRITE(*,*) "          -def <x|y|z> <strain> [<Poissons ratio>]"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-dislocation" .OR. helpsection=="-disloc") THEN
  WRITE(*,*) "..> Insert a dislocation in the system:"
  WRITE(*,*) "          -disloc <pos1> <pos2> screw <x|y|z> <x|y|z> <b>"
  WRITE(*,*) "          -disloc <pos1> <pos2> <edge|edge_add|edge_rm> <x|y|z> <x|y|z> <b> <ν>"
  WRITE(*,*) "          -disloc <pos1> <pos2> mixed <x|y|z> <x|y|z> <b1> <b2> <b3>"
  WRITE(*,*) "          -disloc loop <x> <y> <z> <x|y|z> <radius> <bx> <by> <bz> <nu>"
  WRITE(*,*) "          -disloc file <file> <nu>"
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
IF(helpsection=="options" .OR. helpsection=="-reduce-cell") THEN
  WRITE(*,*) "..> Reduce the system size by taking advantage of its symmetries:"
  WRITE(*,*) "          -reduce-cell"
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
  WRITE(*,*) "          -rotate [com] [hkl] <angle>"
  WRITE(*,*) "          -rotate [com] vx vy vz <angle>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-roundoff" .OR. helpsection=="-round-off") THEN
  WRITE(*,*) "..> Round off atom coordinates or the values of a property:"
  WRITE(*,*) "          -roundoff <property> <threshold>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-select") THEN
  WRITE(*,*) "..> Select atoms according to a criterion:"
  WRITE(*,*) "          -select all"
  WRITE(*,*) "          -select invert"
  WRITE(*,*) "          -select <species>"
  WRITE(*,*) "          -select <index>"
  WRITE(*,*) "          -select <above|below> <d> <normal>"
  WRITE(*,*) "          -select <in|out> <box|sphere|cylinder|cone|torus> [<axis>] <x1> <y1> [<z1>] [<x2> <y2> <z2>] [R [r]]"
  WRITE(*,*) "          -select prop <prop> <value>"
  WRITE(*,*) "          -select random <N> <species>"
  WRITE(*,*) "          -select <NNN> <species> neighbors <index>"
  WRITE(*,*) "          -select <i> modulo <j>"
  WRITE(*,*) "          -select [add|rm|intersect|xor] <any of the above>"
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
  WRITE(*,*) "          -sort <s|x|y|z> <up|down|pack|reverse>"
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
  WRITE(*,*) "..> Swap two atoms or two atom species or two Cartesian axis:"
  WRITE(*,*) "          -swap <id1> <id2>"
  WRITE(*,*) "          -swap <sp1> <sp2>"
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
WRITE(*,*) "                            |  INPUT  | OUTPUT"
WRITE(*,*) "    ------------------------+---------+--------"
WRITE(*,*) "    atsk (Atomsk format)    |   yes   |  yes"
WRITE(*,*) "    bop (Bond-Order format) |   yes   |  yes"
WRITE(*,*) "    cel (Dr. Probe/EMS)     |   yes   |  yes"
WRITE(*,*) "    coo (COORAT/MBPP)       |   yes   |  yes"
WRITE(*,*) "    cfg (Atomeye)           |   yes   |  yes"
WRITE(*,*) "    cif (Cryst.Info.File)   |   yes   |  yes"
WRITE(*,*) "    csv (Comma-Sep.Values)  |   yes   |  yes"
WRITE(*,*) "    d12 (CRYSTAL)           |   yes   |  yes"
WRITE(*,*) "    dat                     |   yes   |  yes"
WRITE(*,*) "    dd  (ddplot)            |    no   | yes (1)"
WRITE(*,*) "    dlp (DL_POLY CONFIG)    |   yes   |  yes"
WRITE(*,*) "    fdf (SIESTA format)     |   yes   |  yes"
WRITE(*,*) "    gin (GULP input)        |   yes   |  yes"
WRITE(*,*) "    imd (IMD input)         |   yes   |  yes"
WRITE(*,*) "    in (ABINIT input)       |   yes   |  yes"
WRITE(*,*) "    jems (JEMS input)       |   yes   |  yes"
WRITE(*,*) "    lmc (LAMMPS output)     |   yes   |   no"
WRITE(*,*) "    lmp (LAMMPS data)       |   yes   |  yes"
WRITE(*,*) "    mol (MOLDY format)      |   yes   |  yes"
WRITE(*,*) "    OUTCAR (POSCAR/VASP)    |  yes(2) |   no"
WRITE(*,*) "    pdb (Protein Data Bank) |   yes   |  yes"
WRITE(*,*) "    POSCAR (POSCAR/VASP)    |   yes   |  yes"
WRITE(*,*) "    pw (Quantum Espresso)   |   yes   |  yes"
WRITE(*,*) "    pwout (QE output file)  |  yes(2) |   no"
WRITE(*,*) "    str (PDFFIT)            |   yes   |  yes"
WRITE(*,*) "    vesta (VESTA file)      |   yes   |  yes"
WRITE(*,*) "    xmd (XMD file)          |   yes   |  yes"
WRITE(*,*) "    xsf (XCrySDen)          |   yes   |  yes"
WRITE(*,*) "    xv (SIESTA format)      |   yes   |  yes"
WRITE(*,*) "    xyz/exyz/sxyz           |   yes   |  yes"
WRITE(*,*) "        (1) Mode ddplot only."
WRITE(*,*) "        (2) Mode unfold only."
ENDIF
!
WRITE(*,*) ""
WRITE(*,*) ">>> Look at the /doc/ folder provided with the program"
WRITE(*,*) "    or go to: https://atomsk.univ-lille.fr/"
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
CHARACTER(LEN=32):: errmsg, warnmsg
CHARACTER(LEN=256):: msg  !The message to be displayed
CHARACTER(LEN=256):: temp, temp2, temp3, temp4
CHARACTER(LEN=*),DIMENSION(:):: strings !Character strings that may be part of the message
INTEGER:: i, j
INTEGER,INTENT(IN):: imsg  !index of message to display
REAL(dp):: tempreal
REAL(dp),DIMENSION(:),INTENT(IN):: reals  !real numbers that may be part of the message
!
!Set colours for error and warning headers
errmsg = COLOUR_MSG("X!X ERROR:",colourerr)
warnmsg = COLOUR_MSG("/!\ WARNING:",colourwarn)
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
    msg = " ___________________________________________________"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    IF( nwarn>0 ) THEN
      WRITE(temp,*) nwarn
      temp = ADJUSTL(temp)
      msg = COLOUR_MSG("/!\ WARNINGS: ",colourwarn)
      msg = "|  "//TRIM(ADJUSTL(msg))//" "//TRIM(temp)
      msg(53:53)="|"
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ENDIF
    IF( nerr>0 ) THEN
      WRITE(temp,*) nerr
      temp = ADJUSTL(temp)
      msg = COLOUR_MSG("X!X ERRORS: ",colourerr)
      msg = "|  "//TRIM(ADJUSTL(msg))//"   "//TRIM(temp)
      msg = TRIM(msg)//"                                 |"
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ENDIF
    msg = "|___________________________________________________|"
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
  CALL DISPLAY_MSG(1,msg,logfile,'NO')
CASE(5)
  !strings(1) = name of file
  msg = "..> OK, I will overwrite "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(6)
  msg = "..> OK, I will overwrite all files from now on."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(7)
  msg = "<?> Enter a name for the file to write:"
  CALL DISPLAY_MSG(1,msg,logfile,'NO')
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
  CALL DISPLAY_MSG(1,msg,logfile,'NO')
CASE(10)
  !reals(1) = index of current element
  !reals(2) = total number of elements
  !SPECIAL: this writes a message on the screen without advancing, so
  ! it is restricted to verbosity levels that display something on screen
  ! (see include/display_messages.f90)
  IF( verbosity==1 .OR. verbosity>=3 ) THEN
    IF(progressbar.NE."none") CALL DISPLAY_PROGBAR(reals(1),reals(2))
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
CASE(15)
  msg = "..> Done."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(16)
  !strings(1) = name of file
  msg = ">>> Reading user configuration file: "//TRIM(ADJUSTL(strings(1)))
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
  CALL DISPLAY_MSG(1,msg,logfile,'NO')
CASE(702)
  msg = "/!\ No input file was specified."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "<?> Please provide the name of an existing file:"
  CALL DISPLAY_MSG(1,msg,logfile,'NO')
CASE(703)
  !strings(1) = name of command line argument
  msg = TRIM(ADJUSTL(warnmsg))//" Unrecognized command-line argument: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(704)
  !strings(1) = name of file that will be ignored
  msg = TRIM(ADJUSTL(warnmsg))//" too many file names passed as arguments."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "            The following file will be ignored: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(750)
  msg = ""
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = TRIM(ADJUSTL(warnmsg))//" YOU SHOULD NOT RUN THIS PROGRAM AS ROOT!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  temp=""
  DO WHILE(temp.NE."ok")
    msg = "    Enter 'ok' to continue anyway, or Ctrl+C to quit."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    READ(*,*) temp
  ENDDO
  msg = ""
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(751)
  !strings(1) = origin of "nthreads": command-line or config. file
  msg = TRIM(ADJUSTL(warnmsg))//" ignoring directive 'nthreads' from "//TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "    This version of Atomsk was compiled without OpenMP support."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
! 800- 899: ERROR MESSAGES
CASE(800)
  msg = TRIM(ADJUSTL(errmsg))//" unrecognized file format."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(801)
  !strings(1) = atom type
  msg = TRIM(ADJUSTL(errmsg))//" unrecognized atom species: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(802)
  !reals(1) = index of atom that caused the error
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" while trying to read atom #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(803)
  !strings(1) = unit
  msg = TRIM(ADJUSTL(errmsg))//" unknown unit: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(804)
  msg = TRIM(ADJUSTL(errmsg))//" number of atoms is zero, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(806)
  msg = TRIM(ADJUSTL(errmsg))//" the systems do not have the same number of atoms!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(807)
  !reals(1) = index of line in file that caused the error on read
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" incorrect data format. Error is apparently on line "//TRIM(ADJUSTL(msg))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(808)
  !strings(1) = string where conversion failed
  msg = TRIM(ADJUSTL(errmsg))//" failed to convert '"//TRIM(strings(1))// &
      & "' to numerical value."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(809)
  !strings(1) = space group H-M symbol which couldn't be identified
  msg = TRIM(ADJUSTL(errmsg))//" failed to identify the space group '"// &
      & TRIM(strings(1))//"'."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(810)
  !reals(1) = space group number
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" invalid space group: "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(811)
  !reals(1) = space group number
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" no access to data of space group "// &
      & TRIM(ADJUSTL(msg))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(812)
  !strings(1) = non-conformal symmetry operation string
  msg = TRIM(ADJUSTL(errmsg))//" failed to interpret the symmetry operations '"// &
      & TRIM(strings(1))//"'."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(813)
  msg = TRIM(ADJUSTL(errmsg))//" no such file or directory"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(814)
  msg = TRIM(ADJUSTL(errmsg))//" Miller vector cannot be [000]."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(815)
  !strings(1) = Miller indices provided by the user (string)
  !reals(1:3) = values of Miller indices
  msg = TRIM(ADJUSTL(errmsg))//" provided Miller indices do not satisfy h+k+i=0: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
  IF( ANY(DABS(reals(1:3))>0.1d0) ) THEN
    WRITE(temp,*) NINT(reals(1))
    WRITE(temp2,*) NINT(reals(2))
    WRITE(temp3,*) NINT(-1.d0*(reals(1)+reals(2)))
    WRITE(temp4,*) NINT(reals(3))
    IF( ANY(DABS(reals(1:3))>9.9d0) ) THEN
      msg = "           Suggested Miller indices: ["//TRIM(ADJUSTL(temp))//"_"// &
          & TRIM(ADJUSTL(temp2))//"_"//TRIM(ADJUSTL(temp3))//"_"//TRIM(ADJUSTL(temp4))//"]."
    ELSE
      msg = "           Suggested Miller indices: ["//TRIM(ADJUSTL(temp))// &
          & TRIM(ADJUSTL(temp2))//TRIM(ADJUSTL(temp3))//TRIM(ADJUSTL(temp4))//"]."
    ENDIF
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(816)
  msg = TRIM(ADJUSTL(errmsg))//" division by zero."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(817)
  !strings(1) = a string that could not be interpreted
  msg = TRIM(ADJUSTL(errmsg))//" unable to interpret this string as Miller indices: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(818)
  !strings(1) = name of the array that was supposed to be resized
  msg = TRIM(ADJUSTL(errmsg))//" there was an error while attempting to resize the array "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(819)
  !reals(1) = estimated number of particles
  msg = TRIM(ADJUSTL(errmsg))//" memory (RAM) is insufficient to perform this operation."
  CALL DISPLAY_MSG(1,msg,logfile)
  IF( reals(1)>0.d0 ) THEN
    tempreal = reals(1)*32.d0
    IF( tempreal>1.d12 ) THEN
      WRITE(temp,'(f32.1)') tempreal/1.d12
      temp = TRIM(ADJUSTL(temp))//" TB"
    ELSEIF( tempreal>1.d9 ) THEN
      WRITE(temp,'(f12.1)') tempreal/1.d9
      temp = TRIM(ADJUSTL(temp))//" GB"
    ELSEIF( tempreal>1.d6 ) THEN
      WRITE(temp,'(f12.1)') tempreal/1.d6
      temp = TRIM(ADJUSTL(temp))//" MB"
    ELSEIF( tempreal>1.d3 ) THEN
      WRITE(temp,'(f12.1)') tempreal/1.d3
      temp = TRIM(ADJUSTL(temp))//" kB"
    ELSE
      WRITE(temp,'(f12.1)') tempreal
      temp = TRIM(ADJUSTL(temp))//" bytes"
    ENDIF
    msg = "          (estimated required memory: "//TRIM(ADJUSTL(temp))//")."
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
  msg = "          Try on a computer with more memory."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(820)
  msg = TRIM(ADJUSTL(errmsg))//" cannot mix up [hkil] and [uvw] Miller notations."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(821)
  !reals(1) = estimated number of atoms
  WRITE(temp,'(f18.0)') reals(1)
  msg = TRIM(ADJUSTL(errmsg))//" this operation results in a very large number of atoms: "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(1,msg,logfile)
  WRITE(temp,*) NATOMS_MAX
  msg = "          Maximum number of atoms that Atomsk can handle: "//TRIM(ADJUSTL(temp))
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
  msg = TRIM(ADJUSTL(warnmsg))//" cannot load auxiliary property: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1701)
  !strings(1) = input file format
  msg = TRIM(ADJUSTL(warnmsg))//" the most probable file format was found to be "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    However the format could not be determined with certainty."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1702)
  !strings(1) = name of parameter
  !strings(2) = name of custom config file
  msg = TRIM(ADJUSTL(warnmsg))//" unknown parameter '"//TRIM(strings(1))//&
      & "' in "//TRIM(strings(2))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1703)
  !strings(1) = name of personal config file
  msg = TRIM(ADJUSTL(warnmsg))//" errors were encountered while reading configuration file "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Some personal parameters may not have been properly set."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1704)
  !strings(1) = name of personal config file
  msg = TRIM(ADJUSTL(warnmsg))//" symmetry operations are not taken into account."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1705)
  msg = TRIM(ADJUSTL(warnmsg))//" both celldm(:) and conventional notation were found."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            The celldm(:) will be used, and conventional notation ignored."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1706)
  msg = TRIM(ADJUSTL(warnmsg))//" cell dimensions are in Bohrs, while atom positions are in angströms."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Cell vectors will be converted to angströms for consistency."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1707) ! invalid symmetry operation string input.
  !strings(1) = failed symmetry operation string
  msg = TRIM(ADJUSTL(warnmsg))//" invalid symmetry operation string '"// &
      & TRIM(strings(1))//"', skipping..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1708)
  !reals(1) = line number
  !reals(2) = expected number of fields
  !reals(3) = actual number of fields
  IF( NINT(reals(2)).NE.NINT(reals(3)) ) THEN
    WRITE(temp,*) NINT(reals(1))
    WRITE(temp2,*) NINT(reals(2))
    WRITE(temp3,*) NINT(reals(3))
    msg = TRIM(ADJUSTL(warnmsg))//" on line "//TRIM(ADJUSTL(temp))//", expected "//TRIM(ADJUSTL(temp2))// &
        & " fields, found "//TRIM(ADJUSTL(temp3))//"."
    CALL DISPLAY_MSG(1,msg,logfile)
    IF( NINT(reals(2))>NINT(reals(3)) ) THEN
      !Actual number of fields is smaller than expected
      msg = "          Missing data will be replaced by naught values."
    ELSE
      !Actual number of fields is larger than expected
      msg = "          Extra data will be ignored."
    ENDIF
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(1709)
  !strings(1) = keyword
  !strings(2) = type of transformation that cannot be performed
  msg = TRIM(ADJUSTL(warnmsg))//" from keyword '"//TRIM(ADJUSTL(strings(1)))//"':"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Cannot "//TRIM(ADJUSTL(strings(2)))//" because no atom position was read yet."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1799)
  msg = TRIM(ADJUSTL(warnmsg))//" the data file had an unknown format. Atomsk tried to extract"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            atomic data from it, but it may be wrong. Tread carefully!"
  CALL DISPLAY_MSG(1,msg,logfile)
!
!1800-1899: ERROR MESSAGES
CASE(1800)
  msg = TRIM(ADJUSTL(errmsg))//" I could not guess the format of this input file!"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Maybe it is a format unsupported by Atomsk yet?"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Anyway it appears I cannot help, I will skip this file."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1801)
  !strings(1) = file name
  !reals(1) = line number
  msg = TRIM(ADJUSTL(errmsg))//" there were errors while reading the file: " &
      & //TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
  IF( NINT(reals(1))>0 ) THEN
    WRITE(temp,*) NINT(reals(1))
    msg = "          Error appears to be at line # "//TRIM(ADJUSTL(temp))
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(1802)
  !strings(1) = bad array
  msg = TRIM(ADJUSTL(errmsg))//" inconsistent array size in "//TRIM(strings(1))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1803)
  msg = TRIM(ADJUSTL(errmsg))//" size of auxiliary properties "// &
             & "is not consistent with number of atoms."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1804)
  msg = TRIM(ADJUSTL(errmsg))//" unknown or unsupported format."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1805)
  !reals(1) = number of particles read
  !reals(2) = number of particles declared
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(errmsg))//" number of atoms read ("//TRIM(ADJUSTL(temp))//") differs from the number"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          of atoms declared ("//TRIM(ADJUSTL(temp2))//")."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1806)
  !reals(1) = index of atom
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" denominator is zero in coordinate of atom #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1807)
  msg = TRIM(ADJUSTL(errmsg))//" file is not in ASCII format, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1808)
  msg = TRIM(ADJUSTL(errmsg))//" levcfg cannot be greater than 2."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1809)
  msg = TRIM(ADJUSTL(errmsg))//" unable to read the supercell parameters."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1810)
  msg = TRIM(ADJUSTL(errmsg))//" unable to read the number of atoms."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1811)
  !reals(1) = index of atom
  WRITE(msg,*) NINT(reals(1))
  IF( reals(1)<0.1d0 ) THEN
    msg = TRIM(ADJUSTL(errmsg))//" atom index #"//TRIM(ADJUSTL(msg))//" is negative."
  ELSE
    msg = TRIM(ADJUSTL(errmsg))//" atom index #"//TRIM(ADJUSTL(msg))//" is greater than the number of atoms."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1812)
  !reals(1) = index of atom
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" unable to read properties of atom #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1813)
  !strings(1) = string where conversion failed
  msg = TRIM(ADJUSTL(errmsg))//" failed to parse symmetry operation in '"// &
      & TRIM(strings(1))//"'."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1814)
  !strings(1) = name of file format
  !strings(2) = name of mode to use to read this file format
  !strings(3) = name of input file
  IF( LEN_TRIM(strings(2))>0 ) THEN
    msg = TRIM(ADJUSTL(errmsg))//" files in "//TRIM(ADJUSTL(strings(1)))//  &
        & " format can only be read with the mode '--"//TRIM(ADJUSTL(strings(2)))//"'."
    CALL DISPLAY_MSG(1,msg,logfile)
    msg = "    Usage:"
    CALL DISPLAY_MSG(1,msg,logfile)
    IF( LEN_TRIM(strings(3))>0 ) THEN
      msg = TRIM(ADJUSTL(strings(3)))
    ELSE
      msg = "<inputfile>"
    ENDIF
    msg = "      atomsk --one-in-all "//TRIM(ADJUSTL(msg))//" <format> [<options>]"
    CALL DISPLAY_MSG(1,msg,logfile)
  ELSE
    msg = TRIM(ADJUSTL(errmsg))//" files in "//TRIM(ADJUSTL(strings(1)))//" format cannot be read."
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(1815)
  msg = TRIM(ADJUSTL(errmsg))//" neighbor list is empty, probably because atoms are too far from one another."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1816)
  msg = TRIM(ADJUSTL(errmsg))//" input file is in unsupported binary format."
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
  !reals(2) = Poisson coefficient
  !strings(1) = direction of deformation: x, y or z
  WRITE(temp,"(f16.3)") reals(1)*100.d0
  SELECT CASE(StrDnCase(strings(1)))
  CASE('x','y','z')
    msg = ">>> Deforming the system by "//TRIM(ADJUSTL(temp))//"% along "//TRIM(strings(1))
    IF(DABS(reals(2))>1.d-16) THEN
      WRITE(temp2,"(f16.3)") reals(2)
      msg = TRIM(ADJUSTL(msg))//", Poisson ratio: "//TRIM(ADJUSTL(temp2))
    ENDIF
  CASE DEFAULT
    msg = ">>> Shearing the system by "//TRIM(ADJUSTL(temp))//"% along "//TRIM(strings(1))
  END SELECT
  msg = TRIM(ADJUSTL(msg))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2059)
  !
CASE(2060)
  msg = "..> System was successfully deformed."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2061)
  !strings(1) = disloctype: screw, edge, edge_add, edge_rm, mixed, loop
  !strings(2) = direction of dislocline: x, y or z, or Miller index, or file name
  !reals(1) = X component of Burgers vector
  !reals(2) = Y component of Burgers vector
  !reals(3) = Z component of Burgers vector
  !reals(4) = positive if C_tensor is defined, zero otherwise
  !reals(5) = pos1, position of dislocation along first axis
  !reals(6) = pos2, position of dislocation along second axis
  !reals(7) = pos3 (only for loops)
  !reals(8) = loop radius
  j=0
  temp = TRIM(ADJUSTL(strings(1)))
  IF(temp(1:4)=="file" .OR. temp(1:5)=="array") THEN
    msg = ">>> Inserting dislocations defined in the file: "//TRIM(ADJUSTL(strings(2)))
  ELSE
    IF(TRIM(temp)=="screw") THEN
      msg = ">>> Inserting a screw dislocation with line along"
    ELSEIF(temp(1:4)=="edge") THEN
      msg = ">>> Inserting an edge dislocation with line along"
    ELSEIF(temp(1:5)=="mixed") THEN
      msg = ">>> Inserting a mixed dislocation with line along"
    ELSEIF(temp(1:4)=="loop") THEN
      msg = ">>> Inserting a dislocation loop in a plane normal to"
      SELECT CASE(strings(2))
      CASE("x","X")
        IF( reals(1)>0.1d0 ) THEN
          j=-1
        ELSEIF( reals(1)<0.1d0 ) THEN
          j=1
        ENDIF
      CASE("y","Y")
        IF( reals(2)>0.1d0 ) THEN
          j=-1
        ELSEIF( reals(2)<0.1d0 ) THEN
          j=1
        ENDIF
      CASE("z","Z")
        IF( reals(3)>0.1d0 ) THEN
          j=-1
        ELSEIF( reals(3)<0.1d0 ) THEN
          j=1
        ENDIF
      END SELECT
    ENDIF
    msg = TRIM(msg)//' '//TRIM(strings(2))//","
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  !
  IF(temp(1:4).NE."file" .AND. temp(1:5).NE."array") THEN
    IF( reals(4)>0.1d0 ) THEN
      msg = "    using anisotropic elasticity,"
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ENDIF
    !
    IF( TRIM(strings(1))=="edge_add" .OR. j>0 ) THEN
      WRITE(msg,"(a34)") "    by inserting a plane of atoms,"
    ELSEIF( TRIM(strings(1))=="edge_rm" .OR. j<0 ) THEN
      WRITE(msg,"(a33)") "    by removing a plane of atoms,"
    ELSEIF( TRIM(strings(1))=="edge" .OR. TRIM(strings(1))=="screw" .OR. TRIM(strings(1))=="mixed" .OR. j==0 ) THEN
      WRITE(msg,"(a41)") "    conserving the total number of atoms,"
    ENDIF
    CALL DISPLAY_MSG(verbosity,msg,logfile)
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
  ENDIF
CASE(2062)
  msg = "..> Searching the solutions to the anisotropic elasticity equations..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2063)
  !reals(1) = number of inserted/removed atoms
  IF( NINT(reals(1)) < 0 ) THEN
    WRITE(msg,*) NINT(ABS(reals(1)))
    msg = "..> "//TRIM(ADJUSTL(msg))//" atoms were removed."
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" atoms were inserted."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2064)
  !reals(1) = axis of elongation (1=X, 2=Y, 3=Z)
  !strings(1) = magnitude of elongation
  IF( NINT(reals(1))==1 ) THEN
    temp = "X"
  ELSEIF( NINT(reals(1))==2 ) THEN
    temp = "Y"
  ELSE
    temp = "Z"
  ENDIF
  msg = "..> Cell length was changed by "//TRIM(ADJUSTL(strings(1)))//" along "//TRIM(ADJUSTL(temp))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2065)
  !reals(1) = number of dislocations inserted (default 1)
  IF( NINT(reals(1)) <= 0 ) THEN
    msg = "..> No dislocation was inserted."
  ELSEIF( NINT(reals(1)) == 1 ) THEN
    msg = "..> One dislocation was successfully inserted."
  ELSE
    WRITE(temp,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(temp))//" dislocations were successfully inserted."
  ENDIF
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
  msg = ">>> Rotating the system to change the crystal orientation..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2072)
  !strings(1:3) = Miller indices for new orientation
  msg = "..> System was successfully rotated, new orientation: "// &
      & TRIM(ADJUSTL(strings(1)))//" "//TRIM(ADJUSTL(strings(2)))//" "//TRIM(ADJUSTL(strings(3)))
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
  msg = "..> Selecting atoms..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2077)
  !strings(1) = region_side: in or out
  !strings(2) = region_geom: sphere or box or cylinder
  !strings(3) = region_dir: axis of cylinder, or "min" or "max"
  !strings(4) = select_multiple: empty or "add" or "rm" or "intersect" or "xor" or "among"
  !reals(1) = region_1(1)
  !reals(2) = region_1(2)
  !reals(3) = region_1(3)
  !reals(4) = region_2(1)
  !reals(5) = region_2(2)
  !reals(6) = region_2(3)
  IF( strings(4)=="rm" ) THEN
    msg = ">>> Un-selecting"
  ELSE
    msg = ">>> Selecting"
  ENDIF
  IF( strings(1)=="all" .OR. strings(1)=="none" ) THEN
    msg = TRIM(ADJUSTL(msg))//" all atoms."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="invert" ) THEN
    msg = ">>> Inverting the selection of atoms..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="index" ) THEN
    IF( strings(2)=="list " ) THEN
      msg = TRIM(ADJUSTL(msg))//" a list of atoms."
    ELSEIF( strings(2)=="range" ) THEN
      WRITE(temp,*) NINT(reals(1))
      WRITE(temp2,*) NINT(reals(2))
      msg = TRIM(ADJUSTL(msg))//" atoms #"//TRIM(ADJUSTL(temp))//" to "//TRIM(ADJUSTL(temp2))//"."
    ELSE
      WRITE(temp,*) NINT(reals(1))
      msg = TRIM(ADJUSTL(msg))//" atom #"//TRIM(ADJUSTL(temp))//"."
    ENDIF
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="above" .OR. strings(1)=="below" ) THEN
    IF( DABS(reals(1))<1.d12 ) THEN
      WRITE(temp,"(3f16.3)") reals(1)
    ELSEIF( reals(1)<-1.d12 ) THEN
      WRITE(temp,"(a4)") "-INF"
    ELSEIF( reals(1)>1.d12 ) THEN
      WRITE(temp,"(a4)") "+INF"
    ENDIF
    msg = TRIM(ADJUSTL(msg))//" atoms "//TRIM(strings(1))//" "//TRIM(ADJUSTL(temp))// &
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
    msg = TRIM(ADJUSTL(msg))//" atoms "//TRIM(temp)//" "//TRIM(strings(2))//"."
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
    ELSEIF(TRIM(strings(2))=="cone") THEN
      WRITE(temp,"(f16.3)") reals(1)
      msg = "..> Axis along "//TRIM(strings(3))//". Tip at: ("//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(2)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(3)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))//"). Opening angle: "
      WRITE(temp,"(f16.3)") reals(4)
      msg = TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(temp))//"°."
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
    msg = TRIM(ADJUSTL(msg))//" atoms with "//TRIM(strings(2))
    WRITE(temp,"(f16.3)") reals(1)
    IF( strings(3)=="min" .OR. strings(3)=="max" ) THEN
      msg = TRIM(ADJUSTL(msg))//" "//TRIM(strings(3))//"imum."
    ELSEIF( reals(4)>2.d0 ) THEN
      WRITE(temp2,"(f16.3)") reals(2)
      msg = TRIM(ADJUSTL(msg))//" between "//TRIM(ADJUSTL(temp))//" and "//TRIM(ADJUSTL(temp2))
    ELSE
      msg = TRIM(ADJUSTL(msg))//" equal to "//TRIM(ADJUSTL(temp))
    ENDIF
    IF(strings(4)=="among") msg = TRIM(ADJUSTL(msg))//" among selected atoms"
    msg = TRIM(ADJUSTL(msg))//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="random" .OR. strings(1)=="rand" ) THEN
    WRITE(temp,*) NINT(reals(1))
    msg = TRIM(ADJUSTL(msg))//" randomly "//TRIM(ADJUSTL(temp))
    IF( NINT(reals(1))<=1 ) THEN
      msg = TRIM(msg)//" atom"
    ELSE
      msg = TRIM(msg)//" atoms"
    ENDIF
    IF( strings(2).NE."any" .AND. strings(2).NE."all" ) THEN
      msg = TRIM(msg)//" of "//TRIM(strings(2))
    ENDIF
    IF(strings(4)=="among") msg = TRIM(ADJUSTL(msg))//" among selected atoms"
    msg = TRIM(msg)//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="random%" ) THEN
    WRITE(temp,'(f16.3)') reals(1)*100.d0
    msg = TRIM(ADJUSTL(msg))//" randomly "//TRIM(ADJUSTL(temp))//"% of atoms"
    IF( strings(2).NE."any" .AND. strings(2).NE."all" ) THEN
      msg = TRIM(msg)//" of "//TRIM(strings(2))
    ENDIF
    IF(strings(4)=="among") msg = TRIM(ADJUSTL(msg))//" among selected atoms"
    msg = TRIM(msg)//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="modulo" .OR. strings(1)=="mod" ) THEN
    !reals(1) = index of an atom
    !reals(2) = modulo
    WRITE(temp,*) NINT(reals(1))
    WRITE(temp2,*) NINT(reals(2))
    msg = TRIM(ADJUSTL(msg))//" atoms with an index equal to "//TRIM(ADJUSTL(temp))//" modulo "//TRIM(ADJUSTL(temp2))//"."
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
      msg = TRIM(ADJUSTL(msg))//" first nearest"//TRIM(temp2)//" neighbors"
    ELSEIF( reals(1)>0.d0 .AND. DBLE(NINT(reals(1)))-reals(1)<1.d-12 ) THEN
      WRITE(temp,*) NINT(reals(1))
      IF( NINT(reals(1))==1 ) THEN
        msg = TRIM(ADJUSTL(msg))//" the first"//TRIM(temp2)//" neighbor"
      ELSE
        msg = TRIM(ADJUSTL(msg))//" the "//TRIM(ADJUSTL(temp))//" nearest"//TRIM(temp2)//" neighbors"
      ENDIF
    ELSE
      WRITE(temp,'(f16.3)') DABS(reals(1))
      msg = TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(temp2))//" neighbors within a radius of "//TRIM(ADJUSTL(temp))//" A"
    ENDIF
    WRITE(temp,*) NINT(reals(2))
    msg = TRIM(ADJUSTL(msg))//" of atom #"//TRIM(ADJUSTL(temp))//"..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="stl" ) THEN
    msg = TRIM(ADJUSTL(msg))//" atoms inside the 3-D shape from the STL file: "//TRIM(ADJUSTL(strings(2)))
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSE
    !Last case: strings(1) should be an atom species
    msg = TRIM(ADJUSTL(msg))//" all "//TRIM(ADJUSTL(strings(1)))//" atoms..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
CASE(2078)
  !reals(1) = number of atoms that were selected
  !reals(2) = number of atoms that were ADDED to the selection
  !reals(3) = number of atoms that were REMOVED from the selection
  !strings(1) = list of atoms selected
  IF( LEN_TRIM(strings(1))>0 ) THEN
    temp3 = " ("//TRIM(ADJUSTL(strings(1)))//")"
  ELSE
    temp3 = ""
  ENDIF
  IF( NINT(reals(2))>0 .OR. NINT(reals(3))>0 ) THEN
    msg = "..>"
    IF( NINT(reals(2))>0 ) THEN
      IF( NINT(reals(2))==1 ) THEN
        msg = TRIM(ADJUSTL(msg))//" 1 atom was added to the selection."
      ELSEIF( NINT(reals(2))>1 ) THEN
        WRITE(temp,*) NINT( reals(2) )
        msg = TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(temp))//" atoms were added to the selection."
      ENDIF
    ENDIF
    IF( NINT(reals(3))>0 ) THEN
      IF( NINT(reals(2))>0 ) THEN
        msg = TRIM(ADJUSTL(msg))//","
      ENDIF
      IF( NINT(reals(3))==1 ) THEN
        msg = TRIM(ADJUSTL(msg))//" 1 atom was removed from the selection."
      ELSEIF( NINT(reals(3))>1 ) THEN
        WRITE(temp,*) NINT( reals(3) )
        msg = TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(temp))//" atoms were removed from the selection."
      ENDIF
    ENDIF
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 atom"//TRIM(temp3)//" is now selected."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(temp,*) NINT( reals(1) )
    msg = "..> "//TRIM(ADJUSTL(temp))//" atoms are now selected"//TRIM(temp3)//"."
  ELSE
    msg = "..> No atom is selected."
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
  msg = ">>> Rotating the system by "//TRIM(ADJUSTL(msg))//"° "
  IF( SCAN(strings(1),'xXyYzZ[]')>0 ) THEN
    msg = TRIM(ADJUSTL(msg))//" around the "//TRIM(strings(1))//" axis."
  ELSE
    msg = TRIM(ADJUSTL(msg))//" around the vector ("//TRIM(strings(1))//")."
  ENDIF
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
  ELSEIF( strings(1)=="selec" ) THEN
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
  !strings(2) = sort order: up, down, pack, random, reverse
  IF( strings(2)=="random" ) THEN
    msg = ">>> Shuffling atom list in random order..."
  ELSEIF( strings(2)=="reverse" ) THEN
    msg = ">>> Reversing atom list..."
  ELSE
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
  !reals(1) = number of box vectors that were modified
  !reals(2) = index of first box vector that was modified
  !reals(3) = index of second box vector that was modified
  IF( NINT(reals(2))==1 ) THEN
    temp = "First"
  ELSEIF( NINT(reals(2))==2 ) THEN
    temp = "Second"
  ELSEIF( NINT(reals(2))==3 ) THEN
    temp = "Third"
  ENDIF
  IF( NINT(reals(1))==0 ) THEN
    msg = "..> Cell vectors were not modified."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> "//TRIM(ADJUSTL(temp))//" cell vector was modified."
  ELSEIF( NINT(reals(1))==2 ) THEN
    IF( NINT(reals(3))==2 ) THEN
      temp2 = "second"
    ELSEIF( NINT(reals(3))==3 ) THEN
      temp2 = "third"
    ENDIF
    msg = "..> "//TRIM(ADJUSTL(temp))//" and "//TRIM(ADJUSTL(temp2))//" cell vectors were modified."
  ELSE
    msg = "..> Cell vectors were modified."
  ENDIF
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
  msg = "..> Anisotropy ratio:  A = 2*C44 / (C11-C12) = "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,'(e16.3)') reals(2)
  msg = "..> Anisotropy factor: H = 2*C44 + C12 - C11 = "//TRIM(ADJUSTL(msg))
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
  msg = "..> Atom velocities were successfully set up."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2113)
  !strings(1) = file name
  msg = "..> Distribution of velocities was written to the file: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2114)
  !strings(1) = "xyz" if user provided max.disp. along X,Y,Z
  !reals(1) = max. displacement of an atom along X
  !reals(2) = max. displacement of an atom along Y
  !reals(3) = max. displacement of an atom along Z
  msg = ">>> Applying a random perturbation to atom positions,"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  IF( strings(1)=="xyz" ) THEN
    !User provided different values along each direction
    WRITE(temp,'(f24.3)') reals(1)
    WRITE(temp2,'(f24.3)') reals(2)
    WRITE(temp3,'(f24.3)') reals(3)
    msg = "..> maximum magnitude: dx="//TRIM(ADJUSTL(temp))//", dy=" &
        & //TRIM(ADJUSTL(temp2))//", dz="//TRIM(ADJUSTL(temp3))//"."
  ELSE
    !User provided norm of max.disp.
    WRITE(msg,'(f24.3)') reals(1)
    msg = "..> maximum magnitude: "//TRIM(ADJUSTL(msg))//" A."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2115)
  !reals(1) = number of atoms that were displaced
  WRITE(temp,*) NINT(reals(1))
  msg = "..> "//TRIM(ADJUSTL(temp))//" atoms were disturbed."
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
  CASE("relative","rel")
    WRITE(temp,"(f16.3)") reals(2)
    WRITE(temp2,"(f16.3)") reals(3)
    WRITE(temp3,"(f16.3)") reals(4)
    WRITE(temp4,*) NINT(reals(1))
    msg = ">>> Adding an atom of "//TRIM(ADJUSTL(strings(1)))//" at ("//TRIM(ADJUSTL(temp))//","//     &
        & TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))//") relative to atom #"//TRIM(ADJUSTL(temp4))//"."
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
  !reals(1) = type of item to swap: 0=Cartesian axes, 1=atom id, 2= atom species, 3=aux.prop.
  !strings(1) = Cartesian axis, or integer
  !strings(2) = same type as strings(1)
  SELECT CASE(NINT(reals(1)))
  CASE(0)
    msg = ">>> Swapping Cartesian axes "//TRIM(ADJUSTL(strings(1)))//" and "//TRIM(ADJUSTL(strings(2)))//"."
  CASE(1)
    msg = ">>> Swapping atoms #"//TRIM(ADJUSTL(strings(1)))//" and #"//TRIM(ADJUSTL(strings(2)))//"."
  CASE(2)
    msg = ">>> Swapping atoms of "//TRIM(ADJUSTL(strings(1)))//" and "//TRIM(ADJUSTL(strings(2)))//"."
  CASE(3)
    msg = ">>> Swapping auxiliary properties '"//TRIM(ADJUSTL(strings(1)))//"' and '"//TRIM(ADJUSTL(strings(2)))//"'."
  END SELECT
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2126)
  !reals(1) = number of values that were swapped
  IF( NINT(reals(1))<=0 ) THEN
    msg = "..> No such atom exist in the system."
  ELSE
    msg = "..> Swap was successful."
  ENDIF
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
  msg = "..> A suitable cell was found, now filling it with atoms..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(temp,*) NINT(reals(1))
  msg = "    (estimated new number of atoms: "//TRIM(ADJUSTL(temp))//")"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2145)
  WRITE(temp,*) NINT(reals(1))
  msg = "..> Cell is now orthorhombic ("//TRIM(ADJUSTL(temp))//" atoms)."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2146)
  msg = "..> Applying displacements to atoms..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2147)
  !strings(1) = name of property that is rounded off
  msg = ">>> Rounding off the values of"
  IF( strings(1)=="AUX" ) THEN
    msg = TRIM(ADJUSTL(msg))//" all auxiliary properties..."
  ELSEIF( strings(1)=="XYZ" ) THEN
    msg = TRIM(ADJUSTL(msg))//" all atoms coordinates..."
  ELSEIF( strings(1)=="X" .OR. strings(1)=="x" ) THEN
    msg = TRIM(ADJUSTL(msg))//" all X coordinates..."
  ELSEIF( strings(1)=="Y" .OR. strings(1)=="y" ) THEN
    msg = TRIM(ADJUSTL(msg))//" all Y coordinates..."
  ELSEIF( strings(1)=="Z" .OR. strings(1)=="z" ) THEN
    msg = TRIM(ADJUSTL(msg))//" all Z coordinates..."
  ELSE
    msg = TRIM(ADJUSTL(msg))//" the property '"//TRIM(ADJUSTL(strings(1)))//"'..."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2148)
  !reals(1) = number of values that were rounded off
  WRITE(temp,*) NINT(reals(1))
  msg = "..> Done, "//TRIM(ADJUSTL(temp))//" values were rounded off."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2149)
  !strings(1) = direction to reduce
  IF( strings(1)=='p' .OR. strings(1)=='P' ) THEN
    msg = ">>> Reducing system into primitive cell..."
  ELSE
    msg = ">>> Reducing system size while preserving periodicity..."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2150)
  !reals(1) = 1 if cell was reduced along X, 0 otherwise
  !reals(2) = 1 if cell was reduced along Y, 0 otherwise
  !reals(3) = 1 if cell was reduced along Z, 0 otherwise
  !reals(4) = new number of atoms
  j=0
  IF( ANY( NINT(reals(1:3))==0 ) ) THEN
    temp = ""
    IF( NINT(reals(1))==1 ) THEN
      temp = "first"
      j=j+1
    ENDIF
    IF( NINT(reals(2))==1 ) THEN
      IF( NINT(reals(1))==1 ) temp = TRIM(ADJUSTL(temp))//" and"
      temp = TRIM(ADJUSTL(temp))//" second"
      j=j+1
    ENDIF
    IF( NINT(reals(3))==1 ) THEN
      IF( NINT(reals(1))==1 ) temp = TRIM(ADJUSTL(temp))//" and"
      temp = TRIM(ADJUSTL(temp))//" third"
      j=j+1
    ENDIF
    IF(j==1) THEN
      temp2 = "cell vector was reduced"
    ELSE
      temp2 = "cell vectors were reduced"
    ENDIF
    msg = "..> The "//TRIM(ADJUSTL(temp))//" "//TRIM(ADJUSTL(temp2))
  ELSE
    msg = "..> The cell was reduced"
  ENDIF
  WRITE(temp,*) NINT(reals(4))
  msg = TRIM(ADJUSTL(msg))//" ("//TRIM(ADJUSTL(temp))//" atoms left)."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2151)
  !strings(1) = operation performed on cell vector (add, rm, set)
  !strings(2) = component of cell vector
  !reals(1) = distance added (or removed)
  WRITE(temp2,'(f16.3)') reals(1)
  IF( strings(2)=="H1" ) THEN
    temp3 = "the first box vector"
  ELSEIF( strings(2)=="H2" ) THEN
    temp3 = "the second box vector"
  ELSEIF( strings(2)=="H3" ) THEN
    temp3 = "the third box vector"
  ELSEIF( strings(2)=="x" .OR. strings(2)=="X" ) THEN
    temp3 = "the X axis"
  ELSEIF( strings(2)=="y" .OR. strings(2)=="Y" ) THEN
    temp3 = "the Y axis"
  ELSEIF( strings(2)=="z" .OR. strings(2)=="Z" ) THEN
    temp3 = "the Z axis"
  ELSEIF( strings(2)=="xy" .OR. strings(2)=="XY" ) THEN
    temp3 = "the XY tilt"
  ELSEIF( strings(2)=="xz" .OR. strings(2)=="XZ" ) THEN
    temp3 = "the XZ tilt"
  ELSEIF( strings(2)=="yx" .OR. strings(2)=="YX" ) THEN
    temp3 = "the YX tilt"
  ELSEIF( strings(2)=="yz" .OR. strings(2)=="YZ" ) THEN
    temp3 = "the YZ tilt"
  ELSEIF( strings(2)=="zx" .OR. strings(2)=="ZX" ) THEN
    temp3 = "the ZX tilt"
  ELSEIF( strings(2)=="zy" .OR. strings(2)=="ZY" ) THEN
    temp3 = "the ZY tilt"
  ELSEIF( strings(2)=="xyz" .OR. strings(2)=="XYZ" ) THEN
    temp3 = "all box vectors"
  ENDIF
  IF( strings(1)=="add" ) THEN
    msg = ">>> Adding "//TRIM(ADJUSTL(temp2))//" A to "//TRIM(ADJUSTL(temp3))//"..."
  ELSEIF( strings(1)=="rm" ) THEN
    msg = ">>> Removing "//TRIM(ADJUSTL(temp2))//" A to "//TRIM(ADJUSTL(temp3))//"..."
  ELSE
    msg = ">>> Setting "//TRIM(ADJUSTL(temp3))//" to "//TRIM(ADJUSTL(temp2))//" A..."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2152)
  msg = "..> Cell vector was modified."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2153)
  msg = ">>> Attempting to re-adjust the cell vectors automatically..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2154)
  !reals(1) = number of shells detected
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = "..> Detected "//TRIM(ADJUSTL(temp))//" cores and "//TRIM(ADJUSTL(temp2))//" shells."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2155)
  !strings(1) = component to be un-tilted
  IF( LEN_TRIM(strings(1))==0 ) THEN
    msg = ">>> Un-tilting the box..."
  ELSE
    msg = ">>> Removing the tilt component "//TRIM(ADJUSTL(strings(1)))//"..."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2156)
  msg = "..> Box was untilted."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2600)
  !strings(1) = first option
  !strings(2) = second option
  msg = "<!> INFO: for better performance, use option '"//TRIM(ADJUSTL(strings(2)))// &
      & "' before option '"//TRIM(ADJUSTL(strings(1)))//"'."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!2700-2799: WARNING MESSAGES
CASE(2700)
  !strings(1) = option name
  msg = TRIM(ADJUSTL(warnmsg))//" could not understand this option: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Try `atomsk --help options` for a summary of options."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2720)
  msg = TRIM(ADJUSTL(warnmsg))//" axis is already aligned, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2721)
  msg = TRIM(ADJUSTL(warnmsg))//" coordinates are already reduced, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2722)
  !strings(1) = atom species
  msg = TRIM(ADJUSTL(warnmsg))//" "//TRIM(strings(1))//" atoms already have shells."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2723)
  !strings(1) = atomic species
  msg = TRIM(ADJUSTL(warnmsg))//" there is no "//TRIM(strings(1))//" atom in the system, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2724)
  msg = TRIM(ADJUSTL(warnmsg))//" deformation is zero, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2725)
  msg = TRIM(ADJUSTL(warnmsg))//" Burgers vector is zero, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2726)
  msg = TRIM(ADJUSTL(warnmsg))//" supercell is very small in one direction normal to the"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    dislocation line! Are you sure you know what you are doing?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2727)
  !reals(1) = index of atom with large displacement
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" displacement for atom "//TRIM(ADJUSTL(msg))// " is large."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2728)
  msg = TRIM(ADJUSTL(warnmsg))//" expansion factors are all equal to 1, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2729)
  msg = TRIM(ADJUSTL(warnmsg))//" no auxiliary property is defined, skipping..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2730)
  !string(1) = property
  msg = TRIM(ADJUSTL(warnmsg))//" No "//TRIM(ADJUSTL(strings(1)))//" found in auxiliary properties, skipping..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2731)
  msg = TRIM(ADJUSTL(warnmsg))//" Hstart = Hend, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2732)
  !strings(1) = name of unknown property
  msg = TRIM(ADJUSTL(warnmsg))//" unknown property: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2733)
  msg = TRIM(ADJUSTL(warnmsg))//" specified radius is negative, no atom will be removed."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2734)
  msg = TRIM(ADJUSTL(warnmsg))//" rotation angle is zero (modulo 2π), skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2735)
  msg = TRIM(ADJUSTL(warnmsg))//" shear is zero, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2736)
  msg = TRIM(ADJUSTL(warnmsg))//" shift vector is naught, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2737)
  msg = TRIM(ADJUSTL(warnmsg))//" species are the same, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2738)
  msg = TRIM(ADJUSTL(warnmsg))//" units are the same, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2739)
  msg = TRIM(ADJUSTL(warnmsg))//" base vectors are not orthonormal."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2740)
  msg = TRIM(ADJUSTL(warnmsg))//" elastic tensor is not symmetric!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2741)
  msg = TRIM(ADJUSTL(warnmsg))//" Poisson ratio is out of range [-1 , 0.5]."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2742)
  !reals(1) = atom index
  WRITE(temp,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" atom index #"//TRIM(ADJUSTL(temp))//" is out of bounds, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2743)
  msg = TRIM(ADJUSTL(warnmsg))//" skew parameters are zero, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2744)
  msg = TRIM(ADJUSTL(warnmsg))//" there is no ionic shell in the system, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2745)
  !reals(1) = number of atoms that will actually be removed
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" cannot select so many atoms, only "// &
      & TRIM(ADJUSTL(msg))//" will be removed."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2746)
  msg = TRIM(ADJUSTL(warnmsg))//" nothing to be done, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2747)
  msg = TRIM(ADJUSTL(warnmsg))//" stress intensity factor K is zero, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2748)
  msg = TRIM(ADJUSTL(warnmsg))//" supercell is very small in one direction normal to the"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    crack line! Are you sure you know what you are doing?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2749)
  !strings(1) = name of property
  !reals(1) = atom index
  WRITE(temp,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" could not assign the value of the property '"//TRIM(ADJUSTL(strings(1)))// &
      & "' to atom #"//TRIM(ADJUSTL(temp))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2750)
  msg = TRIM(ADJUSTL(warnmsg))//" a selection was defined but it does not contain any atom anymore."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "             Selection was cleared, all atoms are now selected."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2751)
  msg = TRIM(ADJUSTL(warnmsg))//" target temperature is zero, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2752)
  msg = TRIM(ADJUSTL(warnmsg))//" no selection is defined, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2753)
  !reals(1) = atom index
  !reals(2) = 0 (selected) or 1 (un-selected)
  WRITE(temp,*) NINT(reals(1))
  IF( NINT(reals(2))>0 ) THEN
    msg = TRIM(ADJUSTL(warnmsg))//" atom #"//TRIM(ADJUSTL(temp))//" was already unselected, skipping."
  ELSE
    msg = TRIM(ADJUSTL(warnmsg))//" atom #"//TRIM(ADJUSTL(temp))//" was already selected, skipping."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2754)
  !strings(1) = type of object (e.g. "dislocation" or "crack" or "loop" or "all"
  !             ("all" means that all points of a loop are out of the box)
  !reals(1) = number of points that are out of the box
  IF( strings(1)=="loop" ) THEN
    WRITE(temp,*) NINT(reals(1))
    msg = TRIM(ADJUSTL(warnmsg))//" "//TRIM(ADJUSTL(temp))//" points of the loop are out of the box."
    CALL DISPLAY_MSG(1,msg,logfile)
  ELSEIF( strings(1)=="all" ) THEN
    msg = TRIM(ADJUSTL(warnmsg))//" all points of the loop are out of the box!"
    CALL DISPLAY_MSG(1,msg,logfile)
    msg = "    Are you sure you know what you are doing?"
    CALL DISPLAY_MSG(1,msg,logfile)
  ELSE
    msg = TRIM(ADJUSTL(warnmsg))//" the position of the "//TRIM(ADJUSTL(strings(1)))//" is out of the box."
    CALL DISPLAY_MSG(1,msg,logfile)
    msg = "    Are you sure you know what you are doing?"
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(2755)
  msg = TRIM(ADJUSTL(warnmsg))//" the factor is zero, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2756)
  !strings(1) = direction
  msg = TRIM(ADJUSTL(warnmsg))//" supercell is quite large along the "//TRIM(ADJUSTL(strings(1)))//" direction!"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Are you sure you know what you are doing?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2757)
  msg = TRIM(ADJUSTL(warnmsg))//" the indices are the same, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2758)
  msg = TRIM(ADJUSTL(warnmsg))//" no operation to apply, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2759)
  msg = TRIM(ADJUSTL(warnmsg))//" dislocation loop radius is too small, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2760)
  msg = TRIM(ADJUSTL(warnmsg))//" cell is already orthorhombic, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2761)
  !strings(1) = "add" or "rm" or "intersect" or "xor"
  msg = TRIM(ADJUSTL(warnmsg))//" cannot modify selection because no selection was previously defined."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2762)
  !strings(1) = elastic tensor stability criterion
  msg = TRIM(ADJUSTL(warnmsg))//" elastic tensor does not comply to stability criterion: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2763)
  msg = TRIM(ADJUSTL(warnmsg))//" modulo is equal to 1, selecting all atoms."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2764)
  msg = TRIM(ADJUSTL(warnmsg))//" no shorter perdiodic vector was found, system remains the same."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2765)
  msg = TRIM(ADJUSTL(warnmsg))//" distance d is zero, skipping..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2766)
  msg = TRIM(ADJUSTL(warnmsg))//" Poisson ratio cannot be used for shear strain, its value will be ignored."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2767)
  msg = TRIM(ADJUSTL(warnmsg))//" box is not tilted, skipping."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2768)
  msg = TRIM(ADJUSTL(warnmsg))//" there is no atom left in the system."
  CALL DISPLAY_MSG(1,msg,logfile)
  !
CASE(2799)
  !strings(1) = name of obsolete option
  !strings(2) = name of new option
  msg = TRIM(ADJUSTL(warnmsg))//" option '"//TRIM(ADJUSTL(strings(1)))//"' is deprecated and will be removed."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Please use option '"//TRIM(ADJUSTL(strings(2)))//"' instead."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!2800-2899: ERROR MESSAGES
CASE(2800)
  !string(1) = axis (if we are here it"s because it is different from x, y or z)
  msg = TRIM(ADJUSTL(errmsg))//" unknown axis: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2801)
  msg = TRIM(ADJUSTL(errmsg))//" the base Hend is not a rotation of Hstart."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Check that angles are equal in Hstart and Hend."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2802)
  !strings(1) = property that was not read properly
  msg = TRIM(ADJUSTL(errmsg))//" while reading "//TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2803)
  msg = TRIM(ADJUSTL(errmsg))//" there were errors while trying to determine supercell vectors."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2804)
  msg = TRIM(ADJUSTL(errmsg))//" there were errors while applying options."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2805)
  !strings(1) = name of unknown option
  msg = TRIM(ADJUSTL(errmsg))//" unknown option: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2806)
  !strings(1) = name of option
  msg = TRIM(ADJUSTL(errmsg))//" non-conform statement in option: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2807)
  !reals(1) = 1 if roots of Eq.(13-85) cannot be found
  !         = 2 if the A_k(n) cannot be calculated
  !         = 3 if the linear equations defining D(n) cannot be solved
  msg = TRIM(ADJUSTL(errmsg))
  IF(NINT(reals(1))==1) THEN
    msg = TRIM(msg)//" unable to determine the P(n), aborting."
  ELSEIF(NINT(reals(1))==2) THEN
    msg = TRIM(msg)//" unable to determine the A_k(n), aborting."
  ELSEIF(NINT(reals(1))==3) THEN
    msg = TRIM(msg)//" unable to determine the D(n), aborting."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2808)
  msg = TRIM(ADJUSTL(errmsg))//" cannot build a mixed dislocation with isotropic elasticity."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          You may either combine an edge and a screw dislocations,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          or define the elastic tensor to use anisotropic elasticity."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2809)
  msg = TRIM(ADJUSTL(errmsg))//" the elastic tensor contains NaN values, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2810)
  msg = TRIM(ADJUSTL(errmsg))//" inconsistent array size for shells."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2811)
  msg = TRIM(ADJUSTL(errmsg))//" dislocline and dislocplane must be normal to each other, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2812)
  msg = TRIM(ADJUSTL(errmsg))//" unable to determine what atom(s) to remove, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2813)
  !strings(1) = string that could not be converted to a number
  msg = TRIM(ADJUSTL(errmsg))//" unable to convert this string into a number: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2814)
  msg = TRIM(ADJUSTL(errmsg))//" there is no atomic system to apply this option to."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2815)
  !strings(1) = name of matrix to invert
  msg = TRIM(ADJUSTL(errmsg))//" unable to invert the matrix "//strings(1)//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2816)
  msg = TRIM(ADJUSTL(errmsg))//" the elastic tensor is not defined, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2817)
  msg = TRIM(ADJUSTL(errmsg))//" the property '"//TRIM(ADJUSTL(strings(1)))//"' is not defined, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2818)
  msg = TRIM(ADJUSTL(errmsg))//" there was an error while reading the STL file, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2819)
  msg = TRIM(ADJUSTL(errmsg))//" unable to find an orthogonal cell from initial cell vectors."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2820)
  !reals(1) = estimated number of atoms
  msg = TRIM(ADJUSTL(errmsg))//" unable to fill new cell with atoms."
  CALL DISPLAY_MSG(1,msg,logfile)
  WRITE(temp,*) NINT(reals(1))
  msg = "          Estimation of number of atoms required: "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2821)
  msg = TRIM(ADJUSTL(errmsg))//" modulo cannot be naught (division by zero)."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2822)
  msg = TRIM(ADJUSTL(errmsg))//" this option does not accept operations with 'box'."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          Please specify distances in angströms."
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
CASE(3007)
  msg = ">>> Output is NULL, no file will be written."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!3700-3799: WARNING MESSAGES
CASE(3700)
  msg = TRIM(ADJUSTL(warnmsg))//" no output file name was specified, please provide one:"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3701)
  !strings(1) = name of output file
  msg = TRIM(ADJUSTL(warnmsg))//" I could not guess the format for the output file: " &
      & //TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "<?> Please specify a format for the output file among the following:"
  CALL DISPLAY_MSG(1,msg,logfile)
  j=0
  DO i=1,CEILING(SIZE(flist,1)/10.0)
    msg = ""
    DO j=1,10
      IF( 10*(i-1)+j <= SIZE(flist,1) ) THEN
        IF( flist(10*(i-1)+j,3)=="yes  " ) THEN
          msg = TRIM(msg)//' '//flist(10*(i-1)+j,1)
        ENDIF
      ENDIF
    ENDDO
    msg = "        "//TRIM(ADJUSTL(msg))
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDDO
CASE(3702)
  msg = TRIM(ADJUSTL(warnmsg))//" ddplot format is available only when using DDPLOT mode."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            No dd file will be output, please refer to the documentation."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3703)
  msg = TRIM(ADJUSTL(warnmsg))//" atom species are not contiguous. Do you want to pack them? ("&
      & //langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    (this will affect only the POSCAR file)"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3704)
  msg = TRIM(ADJUSTL(warnmsg))//" supercell does not form a lower-triangular matrix, which is"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    required by LAMMPS. Do you want to re-align the system? (" &
      & //langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    (This will affect only the LAMMPS output file)"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3705)
  !strings(1) = skew parameter
  msg = TRIM(ADJUSTL(warnmsg))//" triclinic box skew "//TRIM(strings(1))//" is too large."
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
  msg = TRIM(ADJUSTL(warnmsg))//" unable to reduce box skew "//TRIM(strings(1))//", aborting..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3707)
  !reals(1) = number of atoms
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" number of particles if very large: "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    If ddplot cannot display so many atoms, then use Atomsk"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    with the -cut option to reduce the number of atoms."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3708)
  msg = TRIM(ADJUSTL(warnmsg))//" only the first 32 auxiliary properties will be "// &
      & "written to the CFG file."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3709)
  msg = TRIM(ADJUSTL(warnmsg))//" some supercell parameters cannot be written and"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    this data file WILL BE WRONG!!!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3710)
  !strings(1) = file format
  msg = TRIM(ADJUSTL(warnmsg))//" unknown format '"//TRIM(strings(1))//"', skipping..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3711) ! missing occupancy data
  msg = TRIM(ADJUSTL(warnmsg))//" occupancy data is missing, occupancy of all atoms will be set to 1."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3712) ! missing thermal vibration data
  msg = TRIM(ADJUSTL(warnmsg))//" Debye-Waller factors are missing,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            thermal vibration will be set to zero for all atoms."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3713) ! missing absorption data
  msg = TRIM(ADJUSTL(warnmsg))//" absorption factors are missing, they will be set to 0.03 for all atoms."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3714)
  !reals(1) = number of atoms with a zero "type"
  !reals(2) = index of first atom with a zero "type"
  WRITE(temp,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(warnmsg))
  IF( NINT(reals(1))==1 ) THEN
    msg = TRIM(ADJUSTL(msg))//" atom #"//TRIM(ADJUSTL(temp))//" has an invalid 'type'."
  ELSE
    msg = TRIM(ADJUSTL(msg))//" some atoms have an invalid 'type' (first is atom #" &
        & //TRIM(ADJUSTL(temp))//")."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            You may use the option '-remove-property type' to remove atom types,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            or the option '-properties' to set atom types manually."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3715)
  msg = TRIM(ADJUSTL(warnmsg))//" input data contains partial occupancies, which are"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            not supported by some output format(s)."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Some atoms may overlap in the output file(s), which is not physical."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3716)
  msg = TRIM(ADJUSTL(warnmsg))//" data contains ionic shells, which are"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            not supported by some output format(s)."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Shells will be lost in some output file(s)."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3717)
  !reals(1) = total charge
  WRITE(temp,'(f9.3)') reals(1)
  msg = TRIM(ADJUSTL(warnmsg))//" cell has a non-zero electric charge: Q_total = "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3718)
  msg = TRIM(ADJUSTL(warnmsg))//" some atoms of different species have the same 'type'."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            You may use the option '-remove-property type' to remove atom types,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            or the option '-properties' to set atom types manually."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3719)
  !reals(1) = atom "type"
  WRITE(temp,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" new atoms were assigned the 'type' "//TRIM(ADJUSTL(temp))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!3800-3899: ERROR MESSAGES
CASE(3800)
  msg = TRIM(ADJUSTL(errmsg))//" no atom position to write, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3801)
  !strings(1) = name of file
  msg = TRIM(ADJUSTL(errmsg))//" there were errors while writing the file: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3802)
  msg = TRIM(ADJUSTL(errmsg))//" there were errors while writing files."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3803)
  !reals(1) = index of atom that has a NaN coordinate
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" atom #"//TRIM(ADJUSTL(msg))//" has NaN coordinate, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3804)
  !reals(1) = index of atom that has a NaN coordinate
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" auxiliary property of atom #"//&
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
  msg = "help                       Display this help"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "cd                         Change working directory"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = system_ls
  msg(28:) = "List files in current directory"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "pwd                        Print current working directory"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "print                      Print current box vectors and atom positions"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "memory                     Summary of what is in memory"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "create                     Create an atomic system"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "read <file>                Read the <file> and load its content in memory"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "write <file>               Write current system into <file>"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "box <H11> <H22> <H33>      Define dimensions of orthogonal box"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "atom <sp> <x> <y> <z>      Add a new atom in the system with given species and coordinates"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "C11 <value>                Set value of elastic constant (C11,C22,C33,C12,C13,C23,C44,C55,C66)"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "Cij                        Construct elastic tensor based on previous values, and print it"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "clear                      Clear memory (destroy atomic system)"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "quit                       Exit Atomsk"
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
  msg = '    ('//TRIM(flist(1,1))
  DO i=2,SIZE(flist,1)
    msg = TRIM(msg)//','//TRIM(flist(i,1))
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
  !reals(1) >0 if systems are stacked along a direction, <0 otherwise
  IF( reals(1)>0.d0 ) THEN
    msg = ">>> Stacking the systems along "//TRIM(strings(1))//"..."
  ELSE
    msg = ">>> Merging the systems in the same box..."
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
  !reals(1) = maximum distance
  !reals(2) = width of the skin (A)
  msg = ">>> Computing the radial distribution function,"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(temp,'(f16.3)') reals(1)
  WRITE(temp2,'(f16.3)') reals(2)
  msg = "    up to "//TRIM(ADJUSTL(temp))//" A, using a skin width of "//TRIM(ADJUSTL(temp2))//" A."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4052)
  msg = "..> Computing RDFs..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4053)
  !reals(1) = number of files analyzed
  WRITE(temp,*) NINT(reals(1))
  IF( NINT(reals(1))<=0 ) THEN
    msg = "..> No file was analyzed."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> RDF was computed from 1 file."
  ELSE
    msg = "..> RDF was averaged over "//TRIM(ADJUSTL(temp))//" files."
  ENDIF
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
  msg = ">>> Computing per-atom lattice correspondence tensor G..."
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
  msg = ">>> Computing the local symmetry parameter for: "//TRIM(ADJUSTL(strings(1)))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4070)
  !strings(1) = name of file
  !reals(2) = 1 (reference is a unit cell) or 2 (no reference)
  IF( NINT(reals(1))==1 ) THEN
    msg = ">>> Unit cell provided as a reference : "//TRIM(ADJUSTL(strings(1)))//"..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSE
    temp=""
    IF( LEN_TRIM(strings(1))>0 ) temp=": "//TRIM(ADJUSTL(strings(1)))
    msg = ">>> No reference provided. Generating reference atomic environments"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "    by averaging sites from system to be analyzed"//TRIM(ADJUSTL(temp))//"..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
CASE(4071)
  !reals(1) = number of different atomic environments found
  WRITE(temp,*) NINT(reals(1))
  IF( NINT(reals(1))>1 ) THEN
    msg = "..> Done, found "//TRIM(ADJUSTL(temp))//" different atomic environments."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> Done, found 1 atomic environment."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4072)
  !reals(1) = number of files analyzed
  WRITE(temp,*) NINT(reals(1))
  IF( NINT(reals(1))<=0 ) THEN
    msg = "..> No file was analyzed."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 file was analyzed."
  ELSE
    msg = "..> "//TRIM(ADJUSTL(temp))//" files were analyzed."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4073)
  !strings(1) = name of 1st file
  !strings(2) = name of 2nd file
  msg = ">>> Copying auxiliary properties from "//TRIM(ADJUSTL(strings(1))) &
      & //" to "//TRIM(ADJUSTL(strings(2)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4074)
  !reals(1) = number of aux.prop.copied
  !reals(2) = number of new aux.prop.
  !reals(3) = number of overwritten aux.prop.
  !reals(4) = new total number of aux.prop.
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  WRITE(temp3,*) NINT(reals(3))
  msg = ">>> "//TRIM(ADJUSTL(temp))//" properties were copied: "//TRIM(ADJUSTL(temp2))//" new, " &
      & //TRIM(ADJUSTL(temp3))//" overwritten."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(temp,*) NINT(reals(4))
  msg = "..> System now contains "//TRIM(ADJUSTL(temp))//" auxiliary properties."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4075)
  !strings(:) = string containing user values
  msg = "..> Forcing user-defined values: "
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  DO i=1,SIZE(strings)
    IF( LEN_TRIM(strings(i))>0 ) THEN
      msg = "    "//TRIM(ADJUSTL(strings(i)))
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ENDIF
  ENDDO
CASE(4076)
  msg = ">>> Matching atom IDs of two systems..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4077)
  msg = ">>> Searching for equivalent atoms between the two systems..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4078)
  !reals(1) = number of atom pairs matched
  !reals(2) = number of unpaired atoms
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = "..> Done, "//TRIM(ADJUSTL(temp))//" equivalences found"
  IF( NINT(reals(2))>0 ) THEN
    msg = TRIM(ADJUSTL(msg))//", "//TRIM(ADJUSTL(temp2))//" unpaired."
  ELSE
    msg = TRIM(ADJUSTL(msg))//"."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4079)
  msg = "..> Marking new atoms..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4080)
  msg = "..> Sorting atoms of the second system..."
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
  msg = TRIM(ADJUSTL(warnmsg))//" This file does not exist: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Skipping..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4701)
  msg = TRIM(ADJUSTL(warnmsg))//" I did not find the cell vectors in the input file."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    A new bounding box will be generated, but it may not be accurate."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    I advise using option '-prop' or '-cell' to setup the cell vectors."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4702)
  msg = TRIM(ADJUSTL(warnmsg))//" cell is charged, total polarization may be wrong!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4703)
  msg = TRIM(ADJUSTL(warnmsg))//" species has zero charge, skipping..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4705)
  !reals(1) = index of atom that has too many neighbors
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" number of neighbors exceeds 100 for atom #" &
      & //TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4706)
  !strings(1) = atomic species of absent neighbors
  !strings(2) = atomic species of central atom
  !reals(1) = index of atom that has no neighbour
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" I did not find any "//TRIM(strings(1))// &
      &       " neighbor for atom #"//TRIM(ADJUSTL(msg))//" ("//TRIM(strings(2))//")"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4707)
  !reals(1) = index of atom that has a shell with zero charge
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" the shell of ion #"//TRIM(ADJUSTL(msg)) &
      & //" has a zero charge."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4708)
  msg = TRIM(ADJUSTL(warnmsg))//" this grain contains zero atom."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4709)
  !reals(1) = index of atom
  !reals(2) = number of neighbors
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(warnmsg))//" insufficient number of neighbors ("//TRIM(ADJUSTL(temp2))// &
      &              ") for atom #"//TRIM(ADJUSTL(temp))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4710)
  !strings(1) = name of the matrix
  msg = TRIM(ADJUSTL(warnmsg))//" could not compute the matrix "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4711)
  !strings(1) = name of the file
  msg = TRIM(ADJUSTL(warnmsg))//" this file has a different number of atoms and will not be treated: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4712)
  msg = TRIM(ADJUSTL(warnmsg))//" it is recommended to run this mode with the option '-wrap',"
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
  msg = TRIM(ADJUSTL(warnmsg))//" apparently you did not save the system into a file."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Are you sure you want to "//TRIM(ADJUSTL(strings(1)))//"? ("//langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4714)
  !strings(1) = direction along which cell is small
  !reals(1) = new cell size along that direction
  WRITE(temp,'(f16.3)') reals(1)
  msg = TRIM(ADJUSTL(warnmsg))//" final cell has a small dimension along " &
      & //TRIM(ADJUSTL(strings(1)))//", setting it to "//TRIM(ADJUSTL(temp))//" Å."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4715)
  !reals(1) = index of rotation axis
  IF( NINT(reals(1))==1 ) THEN
    temp = "X"
  ELSEIF( NINT(reals(1))==2 ) THEN
    temp = "Y"
  ELSE
    temp = "Z"
  ENDIF
  msg = TRIM(ADJUSTL(warnmsg))//" a 2-D polycrystal is to be constructed normal to the "//TRIM(temp)//" axis,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            but you specified rotation angles around the other Cartesian axes."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Are you sure that you know what you are doing?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4716)
  !reals(1) = X coordinate of node
  !reals(2) = Y coordinate of node
  !reals(3) = Z coordinate of node
  !reals(4) = index of node
  WRITE(temp,'(f12.2)') reals(1)
  WRITE(temp2,'(f12.2)') reals(2)
  WRITE(temp3,'(f12.2)') reals(3)
  WRITE(temp4,*) NINT(reals(4))
  msg = TRIM(ADJUSTL(warnmsg))//" node #"//TRIM(ADJUSTL(temp4))//" was out of bounds and was wrapped, new position: (" &
      & //TRIM(ADJUSTL(temp))//","//TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))//")."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4717)
  !strings(1) = name of matrix
  msg = TRIM(ADJUSTL(warnmsg))//" "//TRIM(ADJUSTL(strings(1)))//" is not an identity matrix."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4718)
  !reals(1) = index of first node
  !reals(2) = index of 2nd node
  !reals(3) = distance between the two nodes
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  WRITE(temp3,'(f18.1)') reals(3)
  msg = TRIM(ADJUSTL(warnmsg))//" nodes #"//TRIM(ADJUSTL(temp))//" and #"//TRIM(ADJUSTL(temp2))//&
      & " are very close to one another (d = "//TRIM(ADJUSTL(temp3))//" A)."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Are you sure that you know what you are doing?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4719)
  !reals(1) = volume of final cell (A^3)
  WRITE(temp,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" final cell is very small ("//TRIM(ADJUSTL(temp))//" A^3)."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Are you sure that you know what you are doing?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4720)
  !strings(1) = chemical symbol
  !reals(1) = file number
  WRITE(temp,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" No "//TRIM(ADJUSTL(strings(1)))// &
      & " atom found in system N. "//TRIM(ADJUSTL(temp))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4721)
  !reals(1) = N. atoms in 1st system
  !reals(2) = N. atoms in 2nd system
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  WRITE(temp3,*) MIN(NINT(reals(1)),NINT(reals(2)))
  msg = TRIM(ADJUSTL(warnmsg))//" number of atoms differ ("//TRIM(ADJUSTL(temp))// &
      & " vs "//TRIM(ADJUSTL(temp2))//")."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Properties will be copied only for the first "//TRIM(ADJUSTL(temp3))//" atoms."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!4800-4899: ERROR MESSAGES FOR MODES
CASE(4800)
  !strings(1) = mode
  msg = TRIM(ADJUSTL(errmsg))//" non-conform statement in mode: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4801)
  !strings(1) = file format (e.g. "xyz", "cfg"...)
  msg = TRIM(ADJUSTL(errmsg))//' this file format is not yet supported in mode "--unfold":' &
      & //TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4802)
  msg = TRIM(ADJUSTL(errmsg))//" both chiral indices cannot be zero, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4803)
  !reals(1) = theoretical number of atoms in nanotube
  !reals(2) = number found by Atomsk
  WRITE(msg,*) TRIM(ADJUSTL(errmsg))//" inconsistent number of atoms in nanotube."
  CALL DISPLAY_MSG(1,msg,logfile)
  WRITE(msg,*) NINT(reals(1))
  WRITE(temp,*) NINT(reals(2))
  msg = "          Theory: "//TRIM(ADJUSTL(msg))//"; Found: "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4804)
  !strings(1) = name of lattice
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
  IF( LEN_TRIM(strings(1))>0 ) THEN
    temp2 = "the "//TRIM(ADJUSTL(strings(1)))//" lattice"
  ELSE
    temp2 = "this lattice"
  ENDIF
  msg = TRIM(ADJUSTL(errmsg))//TRIM(ADJUSTL(temp2))//" requires "//TRIM(ADJUSTL(temp))//" atom species."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4805)
  !strings(1) = structure type
  msg = TRIM(ADJUSTL(errmsg))//" unknown structure: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4806)
  msg = TRIM(ADJUSTL(errmsg))//" no file to merge, aborting"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4807)
  msg = TRIM(ADJUSTL(errmsg))//" all charges are zero, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Charges must be specified with the option -properties."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4808)
  !strings(1) = atomic species that has zero charge
  msg = TRIM(ADJUSTL(errmsg))//" "//TRIM(strings(1))//" ions cannot have zero charge."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4809)
  !strings(1) = atomic species
  msg = TRIM(ADJUSTL(errmsg))//" no "//TRIM(strings(1))//" ion in the system."
  msg = TRIM(ADJUSTL(msg))
CASE(4810)
  !reals(1) = number of particles in first system
  !reals(2) = number of particles in second system
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(errmsg))//" number of particles of the two systems differ: "// &
      & TRIM(ADJUSTL(temp))//"/"//TRIM(ADJUSTL(temp2))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "           Program will EXIT."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4811)
  msg = TRIM(ADJUSTL(errmsg))//" the two systems do not have the same base vectors."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4812)
  msg = TRIM(ADJUSTL(errmsg))//" file does not appear to have animated XSF format, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4813)
  !strings(1) = name of unknown mode
  msg = TRIM(ADJUSTL(errmsg))//" unknown mode: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  !Try to guess what mode the user wanted to use
  temp=""
  IF( INDEX(strings(1),"av")>0 ) THEN
    temp = "--average"
  ELSEIF( INDEX(strings(1),"create")>0 .OR. INDEX(strings(1),"make")>0 ) THEN
    temp = "--create"
  ELSEIF( INDEX(strings(1),"dif")>0 ) THEN
    temp = "--difference"
  ELSEIF( INDEX(strings(1),"dplo")>0 ) THEN
    temp = "--ddplot"
  ELSEIF( INDEX(strings(1),"nie")>0 .OR. INDEX(strings(1),"Nie")>0 ) THEN
    temp = "--nye"
  ELSEIF( INDEX(strings(1),"poly")>0 ) THEN
    temp = "--polycrystal"
  ELSEIF( INDEX(strings(1),"sym")>0 ) THEN
    temp = "--local-symmetry"
  ENDIF
  IF( LEN_TRIM(temp)>0 ) THEN
    msg = "<?> Did you mean to use the mode '"//TRIM(ADJUSTL(temp))//"'?"
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(4814)
  msg = TRIM(ADJUSTL(errmsg))//" only one mode can be used at a time, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4815)
  msg = TRIM(ADJUSTL(errmsg))//" file does not appear to have DL_POLY HISTORY format, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4816)
  !reals(1) = index for core charges
  !reals(2) = index for shell charges
  IF( reals(1)<reals(2) ) THEN
    msg = TRIM(ADJUSTL(errmsg))//" the charges of ionic cores are not defined, aborting."
  ELSE
    msg = TRIM(ADJUSTL(errmsg))//" the charges of ionic shells are not defined, aborting."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4817)
  msg = TRIM(ADJUSTL(errmsg))//" there is no ionic shell in the system, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4818)
  !strings(1) = file format (e.g. "xyz", "cfg"...)
  msg = TRIM(ADJUSTL(errmsg))//" the file "//TRIM(strings(1))//" seems to be empty,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          or does not contain any valid file name."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4819)
  !strings(1) = suggestion for a vector
  !strings(2) = direction of suggested vector (X,Y or Z)
  msg = TRIM(ADJUSTL(errmsg))//" base vectors are not orthogonal."
  CALL DISPLAY_MSG(1,msg,logfile)
  IF( LEN_TRIM(strings(1))>0 ) THEN
    msg = "    Suggested vector along "//TRIM(ADJUSTL(strings(2)))//": "//TRIM(ADJUSTL(strings(1)))
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(4820)
  msg = TRIM(ADJUSTL(errmsg))//" supercell dimensions were not defined, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4821)
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(errmsg))//" number of atoms ("//TRIM(ADJUSTL(temp))// &
      & ") exceeds size of allocated array ("//TRIM(ADJUSTL(temp2))//")."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4822)
  msg = TRIM(ADJUSTL(errmsg))//" no file to be treated, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4823)
  !strings(1) = keyword 1
  !strings(2) = keyword 2
  msg = TRIM(ADJUSTL(errmsg))//" keywords '"//TRIM(ADJUSTL(strings(1)))//"' and '" &
      &  //TRIM(ADJUSTL(strings(2)))//"' are mutually exclusive, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4824)
  !strings(1) = name of unknown command
  msg = TRIM(ADJUSTL(errmsg))//" unknown command: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Type 'help' for a list of available commands."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4825)
  msg = TRIM(ADJUSTL(errmsg))//" Atomsk cannot run within itself!"
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
  msg = TRIM(ADJUSTL(errmsg))//" this mode is not available in interactive mode, please refer to the documentation."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4827)
  msg = TRIM(ADJUSTL(errmsg))//" unable to create lattice with specified orientation."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4828)
  msg = TRIM(ADJUSTL(errmsg))//" at least two cell dimensions are too small, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4829)
  msg = TRIM(ADJUSTL(errmsg))//" box vectors are not linearly independent, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4830)
  msg = TRIM(ADJUSTL(errmsg))//" unable to construct reference environment, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4831)
  msg = TRIM(ADJUSTL(errmsg))//" no node defined in the parameter file '"//TRIM(ADJUSTL(strings(1)))//"'."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4832)
  !reals(1) = index of first node
  !reals(2) = index of 2nd node
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(errmsg))//" nodes #"//TRIM(ADJUSTL(temp))//" and #"//TRIM(ADJUSTL(temp2))//&
      & " are at the same position, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4833)
  !reals(1) = value of IBRION
  WRITE(temp,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" apparently VASP was run with IBRION = "//TRIM(ADJUSTL(temp))//","
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          meaning that this OUTCAR file does not contain any atomic configuration."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4834)
  !strings(1) = name of file
  msg = TRIM(ADJUSTL(errmsg))//" the file "//TRIM(ADJUSTL(strings(1)))// &
      & " does not contain any auxiliary property, aborting."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4835)
  msg = TRIM(ADJUSTL(errmsg))//" No equivalence could be found between atoms of system 1 and those of system 2."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          Are you sure the two systems are of the same sort?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4900)
  msg = TRIM(ADJUSTL(errmsg))//" only one mode can be used at a time."
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
  msg = TRIM(ADJUSTL(warnmsg))//" it is necessary to overwite the existing file " &
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
!VALUES(1): The year
!VALUES(2): The month
!VALUES(3): The day of the month
!VALUES(4): Time difference with UTC in minutes
!VALUES(5): The hour of the day
!VALUES(6): The minutes of the hour
!VALUES(7): The seconds of the minute
!VALUES(8): The milliseconds of the second
!
!
!If it's late at night or in the week end
IF( values(5)>=20 .OR. values(5)<=6 ) THEN
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
!New year
ELSEIF(values(2)==1 .AND. values(3)<=10) THEN
  WRITE(msg,*) values(1)
  msg = TRIM(ADJUSTL(msg))
  msg = "*** HAPPY NEW YEAR "//TRIM(msg)//"!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!March 14 (3/14): Pi day
ELSEIF(values(2)==3 .AND. values(3)==14) THEN
  WRITE(msg,'(a10,f19.16,a2)') "    ( π = ", pi, " )"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!March 31
ELSEIF(values(2)==3 .AND. values(3)==31) THEN
  msg = "  WORLD BACKUP DAY - REMEMBER TO BACKUP YOUR DATA!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!April 1
ELSEIF(values(2)==4 .AND. values(3)==1) THEN
  msg = "              \,-^--._"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "              /`--;-` "
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!April 27: Morse day
ELSEIF(values(2)==4 .AND. values(3)==27) THEN
  !"HAPPY MORSE CODE DAY"
  msg = ".... .- .--. .--. -.-- / -- --- .-. ... . / -.-. --- -.. . / -.. .- -.--"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!May 1
ELSEIF(values(2)==5 .AND. values(3)==1) THEN
  msg = "*** You working? Workers Day is a HOLIDAY! :-)"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!May 4: Star Wars day ("May the fourth/May the force be with you")
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
!July 14
ELSEIF(values(2)==7 .AND. values(3)==14) THEN
  msg = "*** Today is France National Day!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!August 8: international cat day
ELSEIF(values(2)==8 .AND. values(3)==8) THEN
  CALL GEN_NRANDNUMBERS(1,randarray)
  IF( randarray(1)<0.25d0 ) THEN
    msg = "             /'---'\"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "            ( o   o )"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "            = ._Y_. ="
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( randarray(1)<0.5d0 ) THEN
    msg = "   _                ___       _.--."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "   \`.|\..----...-'`   `-._.-'_.-'`"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "   /  ' `         ,       __.--'"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "   )/' _/     \   `-_,   /"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "   `-'"//'" `"\_  ,_.-;_.-\_ '//"'"//',"'
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "       _.-'_./   {_.'   ; /"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "      {_.-``-'         {_/"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( randarray(1)<0.75d0 ) THEN
    msg = "        |\___/|"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "        ) . . ("
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "       =\  :  /="
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "         )===(    _"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "        /     \  //"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "        |  |  | (("
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "       / \ | / \ ))"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "       \  |||  ///"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "        \_M_M_/_/"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSE
    msg = "     _._     _,-'""`-._"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "    (,-.`._,'(       |\`-/|"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "        `-.-' \ )-`( , o o)"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "              `-    \`_`"//'"'//"'-"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
!
!October 31: Halloween
ELSEIF(values(2)==10 .AND. values(3)==31) THEN
  CALL GEN_NRANDNUMBERS(1,randarray)
  IF( randarray(1)<0.25d0 ) THEN
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
  ELSEIF( randarray(1)<0.5d0 ) THEN
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
  ELSEIF( randarray(1)<0.75d0 ) THEN
    msg = "        |\___/|"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "        ) . . ("
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "       =\  :  /="
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "         )===(    _"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "        /     \  //"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "        |  |  | (("
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "       / \ | / \ ))"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "       \  |||  ///"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "        \_M_M_/_/"
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
!Around December 25: Christmas
ELSEIF(values(2)==12 .AND. values(3)>=24 .AND. values(3)<=25) THEN
  msg = " * * *   MERRY CHRISTMAS!   * * *"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
ENDIF
!
END SUBROUTINE DATE_MSG_EN
!
!
!
END MODULE messages_en
