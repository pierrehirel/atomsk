MODULE messages_DE
!
!**********************************************************************************
!*  MESSAGES_DE                                                                   *
!**********************************************************************************
!* This module contains the GERMAN                                                *
!* version of the messages displayed by the Atomsk program.                       *
!**********************************************************************************
!* (C) July 2015 - Juri Barthel                                                   *
!*     Gemeinschaftslabor fuer Elektronenmikroskopie                              *
!*     RWTH Aachen (GERMANY)                                                      *
!*     ju.barthel@fz-juelich.de                                                   *
!* Last modification: P. Hirel - 06 March 2023                                    *
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
!* DISPLAY_LICENSE_DE  displays the license of the program                        *
!* DISPLAY_HELP_DE     displays the help of the program                           *
!* ATOMSK_MSG_DE       all messages used by atomsk                                *
!* DATE_MSG_DE         displays a nice message according to the date              *
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
USE messages_en
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
! of atomsk
!********************************************************
SUBROUTINE DISPLAY_LICENSE_DE()
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
!
msg = ""
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "Atomsk - Ein Programm fuer den Umgang mit atomaren Systemen."
CALL DISPLAY_MSG(verbosity,msg,logfile)
!
CALL DISPLAY_COPYRIGHT()
!
msg = ""
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "Dieses Programm ist freie Software. Sie koennen es unter den Bedingungen"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "der GNU General Public License, wie von der Free Software Foundation veroeffentlicht,"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "weitergeben und/oder modifizieren, entweder gemaeß Version 3 der Lizenz oder"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "(nach Ihrer Option) jeder spaeteren Version."
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = ""
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "Die Veroeffentlichung dieses Programms erfolgt in der Hoffnung, daß es Ihnen"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "von Nutzen sein wird, aber OHNE IRGENDEINE GARANTIE, sogar ohne die implizite"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "Garantie der MARKTREIFE oder der VERWENDBARKEIT FUeR EINEN BESTIMMTEN ZWECK."
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "Details finden Sie in der GNU General Public License."
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = ""
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "Sie sollten ein Exemplar der GNU General Public License zusammen mit"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "diesem Programm erhalten haben. Falls nicht, siehe <http://www.gnu.org/licenses/>."
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = ""
CALL DISPLAY_MSG(verbosity,msg,logfile)
!
END SUBROUTINE DISPLAY_LICENSE_DE
!
!
!********************************************************
! DISPLAY_HELP
! This subroutine displays all or part of
! the help for atomsk
!********************************************************
SUBROUTINE DISPLAY_HELP_DE(helpsection)
!
IMPLICIT NONE
CHARACTER(LEN=16):: helpsection
!
IF(TRIM(helpsection)=="") helpsection="general"
!
!
IF(helpsection=="general" .OR. helpsection=="") THEN
WRITE(*,*) ">>> VERWENDUNG:"
WRITE(*,*) "       atomsk [--<mode>] <inputfile> [<outputfile>] [<formats>] [options...]"
WRITE(*,*) "    wobei [] optionale Parameter sind, and <> durch Werte oder Zeichenketten"
WRITE(*,*) "    ersetzt werden muss."
WRITE(*,*) ""
WRITE(*,*) ">>> BEISPIEL:"
WRITE(*,*) "..> Konvertiert 'file.xsf' ins CFG Format:"
WRITE(*,*) "       atomsk file.xsf cfg"
WRITE(*,*) ""
WRITE(*,*) ">>> FUeR WEITERE HILFE:"
WRITE(*,*) "..> Ueber Modi:     atomsk --help modes"
WRITE(*,*) "..> Ueber Optionen: atomsk --help options"
WRITE(*,*) "..> Ueber Formate:  atomsk --help formats"
ENDIF
!
IF(helpsection=="modes") THEN
  WRITE(*,*) ">>> MODI:"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="interactive") THEN
  WRITE(*,*) "..> Interaktiver Modus:"
  WRITE(*,*) "          atomsk"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="normal") THEN
  WRITE(*,*) "..> Normaler Modus:"
  WRITE(*,*) "          atomsk <inputfile> [<outputfile>] [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="list") THEN
  WRITE(*,*) "..> Listen Modus:"
  WRITE(*,*) "          atomsk -L <listfile> [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="gather" .OR. helpsection=="all-in-one" .OR. helpsection=="ai1") THEN
  WRITE(*,*) "..> Gather Modus:"
  WRITE(*,*) "          atomsk --gather <listfile> <outputfile> [options] [<formats>]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="unfold"  .OR. helpsection=="1ia" .OR. helpsection=="one-in-all") THEN
  WRITE(*,*) "..> Unfold Modus:"
  WRITE(*,*) "          atomsk --unfold <inputfile> [<outputfile>] [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="create") THEN
  WRITE(*,*) "..> Erstellen Modus:"
  WRITE(*,*) "          atomsk --create <structure> <a> [<c>] <species> <outputfile> [orient hkl hkl hkl] [<formats>] [options]"
  WRITE(*,*) "                            <structure> | N.lattice cst. | N.at.sp."
  WRITE(*,*) "                           -------------+----------------+----------"
  WRITE(*,*) "             KUBISCHE               sc  |       1        |     1"
  WRITE(*,*) "             GITTER                bcc  |       1        | 1 oder 2"
  WRITE(*,*) "                                   fcc  |       1        | 1 oder 2"
  WRITE(*,*) "                               diamond  |       1        | 1 oder 2"
  WRITE(*,*) "                                  L1_2  |       1        |     2"
  WRITE(*,*) "                              fluorite  |       1        |     2"
  WRITE(*,*) "                             rock-salt  |       1        |     2"
  WRITE(*,*) "                            perovskite  |       1        |     3"
  WRITE(*,*) "                                   C15  |       1        |     2"
  WRITE(*,*) "                           -------------+----------------+----------"
  WRITE(*,*) "             TETRAGONALE            st  |  2 (a und c)   | 1 oder 2"
  WRITE(*,*) "             GITTER                bct  |  2 (a und c)   | 1 oder 2"
  WRITE(*,*) "                                   fct  |  2 (a und c)   | 1 oder 2"
  WRITE(*,*) "                           -------------+----------------+----------"
  WRITE(*,*) "             HEXAGONALE            hcp  |  2 (a und c)   | 1 oder 2"
  WRITE(*,*) "             GITTER           wurtzite  |  2 (a und c)   |     2"
  WRITE(*,*) "                              graphite  |  2 (a und c)   | 1 oder 2"
  WRITE(*,*) "                                   C14  |  2 (a und c)   |     2"
  WRITE(*,*) "                                   C36  |  2 (a und c)   |     2"
  WRITE(*,*) "          atomsk --create nanotube <a> <m> <n> <sp1> [<sp2>] [options] <outputfile> [<formats>]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="cpprop") THEN
  WRITE(*,*) "..> Mode copy properties:"
  WRITE(*,*) "          atomsk --copy-properties <file1> <file2> [options] <outputfile> [<formats>]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="ddplot") THEN
  WRITE(*,*) "..> DDplot Modus:"
  WRITE(*,*) "          atomsk --ddplot <file1> <file2> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="merge") THEN
  WRITE(*,*) "..> Zusammenfassen Modus:"
  WRITE(*,*) "          atomsk -M [<x|y|z>] <Nfiles> <file1>...<fileN> <outputfile> [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="unwrap") THEN
  WRITE(*,*) "..> Entfalten Modus:"
  WRITE(*,*) "          atomsk --unwrap <reference> <system> [<outputfile>] [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="density") THEN
  WRITE(*,*) "..> Density Modus:"
  WRITE(*,*) "          atomsk --density <file> <property> <1d|2d|3d> <x|y|z> <sigma> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="difference") THEN
  WRITE(*,*) "..> Difference Modus:"
  WRITE(*,*) "          atomsk --difference <file1> <file2> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="edm") THEN
  WRITE(*,*) "..> Modus zu elektrischen Dipolmomenten:"
  WRITE(*,*) "          atomsk --edm <system> <Pspecies> <NNN> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="PE") THEN
  WRITE(*,*) "..> Modus zur elektrischen Polarisation:"
  WRITE(*,*) "          atomsk -PE <system> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="rdf") THEN
  WRITE(*,*) "..> Modus zur radialen Verteilungsfunktion:"
  WRITE(*,*) "          atomsk --rdf <listfile> <Rmax> <dR> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="polycrystal") THEN
  WRITE(*,*) "..> Modus Polykristall:"
  WRITE(*,*) "          atomsk --polycrystal <unitcell> <parameters> <outputfile> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="average") THEN
  WRITE(*,*) "..> Modus Mittelwert:"
  WRITE(*,*) "          atomsk --average <listfile> <outputfile> [options] [<formats>]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="nye") THEN
  WRITE(*,*) "..> Modus Nye:"
  WRITE(*,*) "          atomsk --nye <reference> <defective> <outputfile> [options] [<formats>]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="interpolate") THEN
  WRITE(*,*) "..> Modus Interpolate:"
  WRITE(*,*) "          atomsk --interpolate <file1> <file2> <N> [options] [<formats>]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="cs") THEN
  WRITE(*,*) "..> Modus Lokale Symmetrie:"
  WRITE(*,*) "          atomsk --local-symmetry <file> [options] [<formats>]"
ENDIF
!
IF(helpsection=="options") THEN
  WRITE(*,*) ">>> OPTIONEN (Abstaende=Angstroems, Winkel=Grad):"
  WRITE(*,*) ""
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-add-atom" .OR. helpsection=="-add-atoms" .OR. &
  &helpsection=="-addatom" .OR. helpsection=="-addatoms" ) THEN
  WRITE(*,*) "..> Fügen Sie dem System neue Atome hinzu:"
  WRITE(*,*) "          -add-atom <species> at <x> <y> <z>"
  WRITE(*,*) "          -add-atom <species> relative <index> <x> <y> <z>"
  WRITE(*,*) "          -add-atom <species> near <index>"
  WRITE(*,*) "          -add-atom <species> random <N>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-add-shells".OR. helpsection=="-as" &
  & .OR. helpsection=="-create-shells".OR. helpsection=="-cs") THEN
  WRITE(*,*) "..> Erstellen Sie Shells für einige oder alle Atome:"
  WRITE(*,*) "          -add-shells <all|species>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-alignx") THEN
  WRITE(*,*) "..> Lege den ersten Zellenvektor auf die X Achse:"
  WRITE(*,*) "          -alignx"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-bind-shells" .OR. helpsection=="-bs") THEN
  WRITE(*,*) "..> Weise Ionenhuellen den entsprechenden Kernen zu:"
  WRITE(*,*) "          -bind-shells"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-cell") THEN
  WRITE(*,*) "..> Zellvektoren ändern:"
  WRITE(*,*) "          -cell <add|rm|set> <length>  <H1|H2|H3|x|y|z|xy|xz|yx|yz|zx|zy|xyz>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-center") THEN
  WRITE(*,*) "..> Legen Sie ein Atom oder der Masseschwerpunkt des Systems in der Mitte der Box:"
  WRITE(*,*) "          -center <index>"
  WRITE(*,*) "          -center com"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-crack") THEN
  WRITE(*,*) "..> Fuege einen Riss ins System ein:"
  WRITE(*,*) "          -crack <I|II|III> <stress|strain> <K> <pos1> <pos2> "//&
           &            "<crackline> <crackplane> <μ> <ν>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-create-shells".OR. helpsection=="-cs") THEN
  WRITE(*,*) "..> Erzeuge Schalen fuer einige oder alle Atome:"
  WRITE(*,*) "          -cs <all|species>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-cut") THEN
  WRITE(*,*) "..> Schneide einen Teil des Systems ab:"
  WRITE(*,*) "          -cut <above|below> <cutdistance> <x|y|z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-deform" .OR. helpsection=="-def") THEN
  WRITE(*,*) "..> Fuehre uniaxiale Spannung oder Verzerrung ein:"
  WRITE(*,*) "          -def <x|y|z> <strain> <Poissons ratio>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-dislocation" .OR. helpsection=="-disloc") THEN
  WRITE(*,*) "..> Fuege eine Versetzung ins System ein:"
  WRITE(*,*) "          -disloc <pos1> <pos2> screw <x|y|z> <x|y|z> <b>"
  WRITE(*,*) "          -disloc <pos1> <pos2> <edge|edge_add|edge_rm> <x|y|z> <x|y|z> <b> <ν>"
  WRITE(*,*) "          -disloc <pos1> <pos2> mixed <x|y|z> <x|y|z> <b1> <b2> <b3>"
  WRITE(*,*) "          -disloc loop <x> <y> <z> <x|y|z> <radius> <bx> <by> <bz> <nu>"
  WRITE(*,*) "          -disloc file <Datei> <nu>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-duplicate" .OR. helpsection=="-dup") THEN
  WRITE(*,*) "..> Vervielfache das System entlang der drei Raumrichtungen:"
  WRITE(*,*) "          -duplicate <Nx> <Ny> <Nz>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-fix") THEN
  WRITE(*,*) "..> Fixiere einige Atome:"
  WRITE(*,*) "          -fix <x|y|z> <above|below> <distance> <x|y|z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-fractional" .OR. helpsection=="-frac") THEN
  WRITE(*,*) "..> Konvertiere Koordinaten zu fraktionalen Koordinaten:"
  WRITE(*,*) "          -frac"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-mirror") THEN
  WRITE(*,*) "..> Spiegelung an einer Ebene:"
  WRITE(*,*) "          -mirror <d> <normal>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-options") THEN
  WRITE(*,*) "..> Wende eine Liste von Optionen aus einer Datei an:"
  WRITE(*,*) "          -options <file>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-orient") THEN
  WRITE(*,*) "..> Aendere die kristallographische Orientierung des Systems:"
  WRITE(*,*) "          -orient <Hx> <Hy> <Hz> <H'x> <H'y> <H'z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-prop") THEN
  WRITE(*,*) "..> Lese Systemeigenschaften ein:"
  WRITE(*,*) "          -prop <file>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-rebox") THEN
  WRITE(*,*) "..> (Erneute) Berechung der Vektoren des Begrenzungsrahmens:"
  WRITE(*,*) "          -rebox"
!
IF(helpsection=="options" .OR. helpsection=="-reduce-cell") THEN
  WRITE(*,*) "..> Reduzieren Sie das System in eine Einheitszelle:"
  WRITE(*,*) "          -reduce-cell"
ENDIF
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-remove-atom" .OR. helpsection=="-rmatom") THEN
  WRITE(*,*) "..> Entferne ein Atom nach Index, oder alle Atome einer Spezies:"
  WRITE(*,*) "          -rmatom <index|species>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-remove-doubles" .OR. helpsection=="-rmd") THEN
  WRITE(*,*) "..> Entferne Atomduplikate oder nahe Atome:"
  WRITE(*,*) "          -rmd <distance>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-remove-property" .OR. helpsection=="-rmprop") THEN
  WRITE(*,*) "..> Entferne eine oder alle Hilfseigenschaften:"
  WRITE(*,*) "          -rmprop <property>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-roll" .OR. helpsection=="-bend") THEN
  WRITE(*,*) "..> Rollen das System um eine Achse:"
  WRITE(*,*) "          -roll <x|y|z> <angle> <x|y|z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-rotate" .OR. helpsection=="-rot") THEN
  WRITE(*,*) "..> Rotiere das System um eine Achse:"
  WRITE(*,*) "          -rotate [com] <x|y|z> <angle>"
  WRITE(*,*) "          -rotate [com] [hkl] <angle>"
  WRITE(*,*) "          -rotate [com] vx vy vz <angle>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-roundoff" .OR. helpsection=="-round-off") THEN
  WRITE(*,*) "..> Runden Atomkoordinaten oder die Werte einer Eigenschaft ab:"
  WRITE(*,*) "          -roundoff <property> <threshold>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-select") THEN
  WRITE(*,*) "..> Waehle Atome nach verschiedenen Kriterien aus:"
  WRITE(*,*) "          -select all"
  WRITE(*,*) "          -select invert"
  WRITE(*,*) "          -select <species>"
  WRITE(*,*) "          -select <index>"
  WRITE(*,*) "          -select <above|below> <d> <normal>"
  WRITE(*,*) "          -select <in|out> <box|sphere|cylinder|cone|torus> [<Achse>] <x1> <y1> <z1> <x2> [<y2> <z2>] [alpha]]"
  WRITE(*,*) "          -select prop <prop> <value>"
  WRITE(*,*) "          -select random <N> <species>"
  WRITE(*,*) "          -select <NNN> <species> neighbors <index>"
  WRITE(*,*) "          -select <i> modulo <j>"
  WRITE(*,*) "          -select [add|rm|intersect|xor] <any of the above>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-separate") THEN
  WRITE(*,*) "..> Separate Atome, die zu nahe sind:"
  WRITE(*,*) "          -separate <distance> <shift>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-shear") THEN
  WRITE(*,*) "..> Wende einfache Scherung oder Verzerrung auf das System an:"
  WRITE(*,*) "          -shear <x|y|z> <amplitude> <x|y|z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-shift") THEN
  WRITE(*,*) "..> Verschiebe einen Teil des Systems:"
  WRITE(*,*) "          -shift <above|below> <distance> <x|y|z> <tauX> <tauY> <tauZ>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-sort") THEN
  WRITE(*,*) "..> Sortiere Atome:"
  WRITE(*,*) "          -sort <s|x|y|z> <up|down|pack|reverse>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-stress") THEN
  WRITE(*,*) "..> Wende eine mechanische Spannung an:"
  WRITE(*,*) "          -stress <xx|yy|zz|xy|xz|yz|P> <value(GPa)>"
  WRITE(*,*) "          -stress <file>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-substitute" .OR. helpsection=="-sub") THEN
  WRITE(*,*) "..> Ersetze Atomsorte sp1 durch Sorte sp2:"
  WRITE(*,*) "          -sub <sp1> <sp2>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-swap") THEN
  WRITE(*,*) "..> Vertauschen zwei Atomen:"
  WRITE(*,*) "          -swap <id1> <id2>"
  WRITE(*,*) "          -swap <sp1> <sp2>"
  WRITE(*,*) "          -swap <x|y|z> <x|y|z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-torsion") THEN
  WRITE(*,*) "..> Torsion um eine Achse auftragen:"
  WRITE(*,*) "          -torsion <x|y|z> <angle>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-unit" .OR. helpsection=="-u") THEN
  WRITE(*,*) "..> Konvertiere Koordinaten in eine andere Einheit:"
  WRITE(*,*) "          -u <unit1> <unit2>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-unskew") THEN
  WRITE(*,*) "..> Reduziere Schraegheit der Superzelle:"
  WRITE(*,*) "          -unskew"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-velocity") THEN
  WRITE(*,*) "..> Setze Geschwindigkeiten der Atome entsprechend der Maxwell-Boltzmann Verteilung:"
  WRITE(*,*) "          -velocity <Temperature>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-wrap") THEN
  WRITE(*,*) "..> Setze Atome periodisch zurueck in die Superzelle:"
  WRITE(*,*) "          -wrap"
ENDIF
!
! -- add other options in alphabetical order --
!
!
IF(helpsection=="formats") THEN
WRITE(*,*) ">>> FORMATE:"
WRITE(*,*) "    Formate muessen so angegeben werden wie in der ersten Spalte."
WRITE(*,*) "    Atomsk konvertiert von jedem 'ja' in der INPUT Spalte"
WRITE(*,*) "    zu jedem 'ja' in der OUTPUT Spalte:"
WRITE(*,*) "                            |  INPUT  | OUTPUT"
WRITE(*,*) "    ------------------------+---------+--------"
WRITE(*,*) "    atsk (Atomsk format)    |  ja     |  ja "
WRITE(*,*) "    bop (Bond-Order format) |  ja     |  ja "
WRITE(*,*) "    cel (Dr. Probe/EMS)     |  ja     |  ja "
WRITE(*,*) "    coo (COORAT/MBPP)       |  ja     |  ja "
WRITE(*,*) "    cfg (Atomeye)           |  ja     |  ja "
WRITE(*,*) "    cif (Cryst.Info.File)   |  ja     |  ja "
WRITE(*,*) "    d12 (CRYSTAL)           |  ja     |  ja"
WRITE(*,*) "    dat                     |  ja     |  ja"
WRITE(*,*) "    dd  (ddplot)            | nein    | ja (1)"
WRITE(*,*) "    dlp (DL_POLY CONFIG)    |  ja     |  ja "
WRITE(*,*) "    fdf (SIESTA format)     |  ja     |  ja "
WRITE(*,*) "    gin (GULP input)        |  ja     |  ja "
WRITE(*,*) "    imd (IMD input)         |  ja     |  ja "
WRITE(*,*) "    in (ABINIT input)       |  ja     |  ja"
WRITE(*,*) "    jems (JEMS input)       |  ja     |  ja "
WRITE(*,*) "    lmc (LAMMPS output)     |  ja     | nein"
WRITE(*,*) "    lmp (LAMMPS data)       |  ja     |  ja "
WRITE(*,*) "    mol (MOLDY format)      |  ja     |  ja "
WRITE(*,*) "    OUTCAR (POSCAR/VASP)    | ja(2)   | nein"
WRITE(*,*) "    pdb (Protein Data Bank) |  ja     |  ja "
WRITE(*,*) "    POSCAR (POSCAR/VASP)    |  ja     |  ja "
WRITE(*,*) "    pw (Quantum Espresso)   |  ja     |  ja "
WRITE(*,*) "    pwout (QE output file)  | ja (2)  | nein"
WRITE(*,*) "    str (PDFFIT)            |  ja     |  ja"
WRITE(*,*) "    vesta (VESTA file)      |  ja     |  ja"
WRITE(*,*) "    xmd (XMD file)          |  ja     |  ja "
WRITE(*,*) "    xsf (XCrySDen)          |  ja     |  ja "
WRITE(*,*) "    xv (SIESTA format)      |  ja     |  ja "
WRITE(*,*) "    xyz/exyz/sxyz           |  ja     |  ja "
WRITE(*,*) "        (1) nur im ddplot Modus."
WRITE(*,*) "        (2) nur im unfold Modus."
ENDIF
!
WRITE(*,*) ""
WRITE(*,*) ">>> Schau im mitgelieferten /doc/ Ordner nach oder "
WRITE(*,*) "    besuche: https://atomsk.univ-lille.fr/"
WRITE(*,*) ""
!
!
END SUBROUTINE DISPLAY_HELP_DE
!
!
!
!********************************************************
! ATOMSK_CREATE_DATE
! This routine 
!********************************************************
SUBROUTINE ATOMSK_CREATE_DATE_DE(VALUES,username,msg)
!
IMPLICIT NONE
CHARACTER(LEN=128),INTENT(IN):: username
INTEGER,DIMENSION(8),INTENT(IN):: VALUES
CHARACTER(LEN=128),INTENT(OUT):: msg
!
WRITE(msg,'(i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') &
  & VALUES(1), "-", VALUES(2),"-", VALUES(3)," ", VALUES(5), ":", VALUES(6), ":", VALUES(7)
!
msg = 'Datei mit Atomsk von '//TRIM(ADJUSTL(username))//' am '//TRIM(ADJUSTL(msg))//" generiert."
!
END SUBROUTINE ATOMSK_CREATE_DATE_DE
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
SUBROUTINE ATOMSK_MSG_DE(imsg,strings,reals)
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
errmsg = COLOUR_MSG("X!X FEHLER:",colourerr)
warnmsg = COLOUR_MSG("/!\ WARNUNG:",colourwarn)
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
  msg = "\o/ Programm erfolgreich beendet!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(temp,"(f30.3)") reals(1)
  WRITE(temp2,"(f30.3)") reals(2)
  WRITE(msg,*) "   Gesamtzeit: "//TRIM(ADJUSTL(temp))//         &
           & " s.; CPU zeit: "//TRIM(ADJUSTL(temp2))//" s."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  IF( nwarn>0 .OR. nerr>0 ) THEN
    !In case of warning or error, display a big box
    msg = " ___________________________________________________"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    IF( nwarn>0 ) THEN
      WRITE(temp,*) nwarn
      temp = ADJUSTL(temp)
      msg = "|  /!\ WARNUNGEN: "//TRIM(temp)
      msg = msg(1:52)//"|"
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ENDIF
    IF( nerr>0 ) THEN
      WRITE(temp,*) nerr
      temp = ADJUSTL(temp)
      msg = COLOUR_MSG("X!X FEHLER: ",colourerr)
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
  WRITE(msg,*) "   Gesamtzeit: "//TRIM(ADJUSTL(temp))//         &
           & " s.; CPU zeit: "//TRIM(ADJUSTL(temp2))//" s."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3)
  msg = "..> Dies kann einige Zeit in Anspruch nehmen..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4)
  !strings(1) = name of file that already exists
  msg = "<?> Diese Datei ist bereits vorhanden: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Moechten Sie sie ueberschreiben? ("//langyes//"/"//langno//") ("&
      & //langBigYes//"=alles ueberschreiben)?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5)
  !strings(1) = name of file
  msg = "..> OK, ich will "//TRIM(strings(1))//" ueberschreiben."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(6)
  msg = "..> OK, zukuenftig will ich alle Dateien ueberschreiben."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(7)
  msg = "<?> Geben Sie einen Namen fuer die zu schreibende Datei an:"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(8)
  !strings(1) = name of file
  IF( LEN_TRIM(strings(1))>0 ) THEN
    msg = "<?> Diese Datei existiert nicht: "//TRIM(strings(1))
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
  msg = "    Bitte geben Sie den Namen einer vorhandenen Datei an:"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    (Geben Sie '"//TRIM(system_ls)//"' fuer eine Liste der Dateien im aktuellen Verzeichnis)"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(9)
  msg = "<?> Geben Sie den Namen einer vorhandenen Datei an:"
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
    IF( NINT(reals(1)) >= NINT(reals(2)) ) THEN
      WRITE(*,'(a)',ADVANCE="YES") ""
    ENDIF
  ENDIF
CASE(11)
  msg = ">>> Konstruiere Nachbarliste..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(12)
  msg = ">>> Wenn Sie Atomsk in Ihrer Arbeit nutzen, bitte zietieren Sie den folgenden Artikel:"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "    Pierre Hirel, Comput. Phys. Comm. 197 (2015) 212"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(14)
  msg = "<i> Ausgewählte Atome wurden entfernt: Die zuvor definierte Auswahl wurde gelöscht."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(15)
  msg = "..> Fertig."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(16)
  !strings(1) = name of file
  msg = ">>> Benutzerkonfigurationsdatei lesen: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
! 700- 799: WARNUNG MESSAGES
CASE(700)
  !strings(1) = first file name
  !strings(2) = second file name
  msg = "/!\ Beide Datei vorhanden: "//TRIM(strings(1))//" und "//TRIM(strings(2))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "<?> Bitte geben Sie an, welche als Eingabedatei verwendet werden soll:"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = '    1- '//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = '    2- '//TRIM(strings(2))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(701)
  !strings(1) = first file name
  !strings(2) = second file name
  msg = "/!\ Keine dieser Dateien ist vorhanden: "//TRIM(strings(1))//" oder "//TRIM(strings(2))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "<?> Bitte geben Sie den Namen einer vorhandenen Datei an:"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(702)
  msg = "/!\ Es wurde Keine Eingabedatei angegeben."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "<?> Bitte geben Sie den Namen einer vorhandenen Datei an:"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(703)
  !strings(1) = name of command line argument
  msg = TRIM(ADJUSTL(warnmsg))//" Nicht erkannte Befehlszeilenargument: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(704)
  !strings(1) = name of file that will be ignored
  msg = TRIM(ADJUSTL(warnmsg))//" zu viele Dateinamen wurden als Argumente übergeben."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "            Die folgende Datei wird ignoriert: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(750)
  msg = ""
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = TRIM(ADJUSTL(warnmsg))//" SIE SOLLTEN DIESES PROGRAM NICHT ALS ROOT EVALUIEREN!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  temp=""
  DO WHILE(temp.NE."ok")
    msg = "    Zum fortfahren drücken Sie 'ok', oder Ctrl+C zum abschließen."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    READ(*,*) temp
  ENDDO
  msg = ""
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(751)
  !strings(1) = origin of "nthreads": command-line or config. file
  IF(strings(1)=="command line") THEN
    temp = "der Kommandozeile"
  ELSE
    temp = "der Datei "//TRIM(ADJUSTL(strings(1)))
  ENDIF
  msg = TRIM(ADJUSTL(warnmsg))//" Ignorieren der Direktive 'nthreads' von "//TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "    Diese Version von Atomsk wurde ohne OpenMP-Unterstützung kompiliert."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
! 800- 899: FEHLER MESSAGES
CASE(800)
  msg = TRIM(ADJUSTL(errmsg))//" Unbekanntes Dateiformat."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(801)
  !strings(1) = atom type
  msg = TRIM(ADJUSTL(errmsg))//" Unerkannte Atomsorte: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(802)
  !reals(1) = index of atom that caused the error
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = "X!X FEHLER beim Versuch, Atom #"//TRIM(ADJUSTL(msg))//" zu lesen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(803)
  !strings(1) = unit
  msg = TRIM(ADJUSTL(errmsg))//" unbekannte Einheit: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(804)
  msg = TRIM(ADJUSTL(errmsg))//" Anzahl der Atome ist Null, Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(806)
  msg = TRIM(ADJUSTL(errmsg))//" die Systeme weisen nicht die gleiche Anzahl von Atomen auf!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(807)
  !reals(1) = index of line in file that caused the error on read
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" falsches Datenformat. Fehler ist anscheinend in Zeile "//TRIM(ADJUSTL(msg))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(808)
  !strings(1) = string where conversion failed
  msg = "X!X FEHLER beim Konvertieren von '"//TRIM(strings(1))// &
      & "' in einen numerischen Wert."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(809)
  !strings(1) = space group H-M symbol which couldn't be identified
  msg = TRIM(ADJUSTL(errmsg))//" unbekannte Raumgruppe: '"//TRIM(strings(1))//"'."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(810)
  !reals(1) = invalid space group number
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" ungueltige Raumgruppe: "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(811)
  !reals(1) = space group number
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = "X!X FEHLER beim Zugriff auf Daten der Raumgruppe "// &
      & TRIM(ADJUSTL(msg))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(812)
  !strings(1) = non-conformal symmetry operation string
  msg = "X!X FEHLER beim Interpretieren der Symmetrieoperationen '"// &
      & TRIM(strings(1))//"'."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(813)
  msg = TRIM(ADJUSTL(errmsg))//" Datei oder Verzeichnis existiert nicht"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(814)
  msg = TRIM(ADJUSTL(errmsg))//" Miller Vektor kann nicht [000] sein."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(815)
  !strings(1) = Miller indices provided by the user (string)
  !reals(1:3) = values of Miller indices
  msg = TRIM(ADJUSTL(errmsg))//" angegebenen Miller-Indizes erfüllen h+k+i=0 nicht: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
  IF( ANY(DABS(reals(1:3))>0.1d0) ) THEN
    WRITE(temp,*) NINT(reals(1))
    WRITE(temp2,*) NINT(reals(2))
    WRITE(temp3,*) NINT(-1.d0*(reals(1)+reals(2)))
    WRITE(temp4,*) NINT(reals(3))
    IF( ANY(DABS(reals(1:3))>9.9d0) ) THEN
      msg = "           Vorgeschlagene Miller-Indizes: ["//TRIM(ADJUSTL(temp))//"_"// &
          & TRIM(ADJUSTL(temp2))//"_"//TRIM(ADJUSTL(temp3))//"_"//TRIM(ADJUSTL(temp4))//"]."
    ELSE
      msg = "           Vorgeschlagene Miller-Indizes: ["//TRIM(ADJUSTL(temp))// &
          & TRIM(ADJUSTL(temp2))//TRIM(ADJUSTL(temp3))//TRIM(ADJUSTL(temp4))//"]."
    ENDIF
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(816)
  msg = TRIM(ADJUSTL(errmsg))//" Durch Null teilen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(817)
  !strings(1) = a string that could not be interpreted
  msg = TRIM(ADJUSTL(errmsg))//" Diese Zeichenfolge kann nicht als Miller-Index interpretiert werden: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(818)
  !strings(1) = name of the array that was supposed to be resized
  msg = TRIM(ADJUSTL(errmsg))//" Beim Versuch die Größe des Arrays zu ändern, ist ein Fehler aufgetreten: "&
      & //TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(819)
  !reals(1) = estimated number of particles
  msg = TRIM(ADJUSTL(errmsg))//" Speicher (RAM) reicht nicht aus, um diese Operation auszuführen."
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
    msg = "          (Schätzung des erforderlichen Speichers: "//TRIM(ADJUSTL(temp))//")."
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
  msg = "          Versuchen Sie einen Computer mit mehr Speicher."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(820)
  msg = TRIM(ADJUSTL(errmsg))//" kann nicht [hkl] und [hkil] Miller-Notationen mischen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(821)
  !reals(1) = estimated number of atoms
  WRITE(temp,'(f18.0)') reals(1)
  msg = "X!X ERROR: Diese Operation führt zu einer sehr großen Anzahl von Atomen: "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(1,msg,logfile)
  WRITE(temp,*) NATOMS_MAX
  msg = "          Anzahl der Atome, die Atomsk verarbeiten kann: "//TRIM(ADJUSTL(temp))
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
  msg = ">>> Oeffnen der Eingabedatei: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(1001)
  !reals(1) = number of atoms
  !reals(2) = number of shells
  WRITE(temp,*) NINT(reals(1))
  IF( NINT(reals(2))>0 ) THEN
    WRITE(temp2,*) NINT(reals(2))
    msg = "..> Eingabedatei wurde erfolgreich gelesen ("//TRIM(ADJUSTL(temp))//" cores + "// &
        & TRIM(ADJUSTL(temp2))//" shells)."
  ELSE
    msg = "..> Eingabedatei wurde erfolgreich gelesen ("//TRIM(ADJUSTL(temp))//" Atomen)."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(1002)
  msg = "..> INP Datei erkannt, entnehme Daten der Superzelle..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(1003)
  msg = "..> POTCAR Datei erkannt, entnehme Atomsorten..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(1004)
  msg = "..> Applying symmetry operations..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!1700-1799: WARNUNG MESSAGES
CASE(1700)
  !strings(1) = auxiliary property that cannot be loaded
  msg = TRIM(ADJUSTL(warnmsg))//" kann Hilfseigenschaften nicht lesen: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1701)
  !strings(1) = input file format
  msg = TRIM(ADJUSTL(warnmsg))//" das wahrscheinlichste Dateiformat ist "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Geringe Sicherheit bei der Bestimmung des Formats."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1702)
  !strings(1) = name of parameter
  !strings(2) = name of custom config file
  msg = TRIM(ADJUSTL(warnmsg))//" Unbekannter Parameter '"//TRIM(strings(1))//&
      & "' in "//TRIM(strings(2))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1703)
  !strings(1) = name of personal config file
  msg = TRIM(ADJUSTL(warnmsg))//" Fehler beim Lesen der Konfigurationsdatei "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "             Einige Einstellungen sind moeglicherweise nicht gesetzt worden."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1704)
  !strings(1) = name of personal config file
  msg = TRIM(ADJUSTL(warnmsg))//" Symmetrie-Opertionen wurden nicht beruecksichtigt."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1705)
  msg = TRIM(ADJUSTL(warnmsg))//" Sowohl celldm(:) als auch die konventionelle Notation liegen vor."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "             celldm(:) wird verwendet, konventionelle Notation wird ignoriert."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1706)
  msg = TRIM(ADJUSTL(warnmsg))//" Zell-Dimensionen in Bohrs, aber atomare Positionen in Angstroems."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "             Zellenvektoren werden in Angstroems konvertiert."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1707) ! invalid symmetry operation string input.
  !strings(1) = failed symmetry operation string
  msg = TRIM(ADJUSTL(warnmsg))//" unzulaessige Zeichenkette fuer "// &
      & "Symmetrieoperation '"//TRIM(strings(1))//"'. Ueberspringe..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1708)
  !reals(1) = line number
  !reals(2) = expected number of fields
  !reals(3) = actual number of fields
  IF( NINT(reals(2)).NE.NINT(reals(3)) ) THEN
    WRITE(temp,*) NINT(reals(1))
    WRITE(temp2,*) NINT(reals(2))
    WRITE(temp3,*) NINT(reals(3))
    msg = TRIM(ADJUSTL(warnmsg))//" in Zeile "//TRIM(ADJUSTL(temp))//", erwartete "//TRIM(ADJUSTL(temp2))// &
        & " Felder, "//TRIM(ADJUSTL(temp3))//" gefunden."
    CALL DISPLAY_MSG(1,msg,logfile)
    IF( NINT(reals(2))>NINT(reals(3)) ) THEN
      !Actual number of fields is smaller than expected
      msg = "          Fehlende Daten werden durch Nullwerte ersetzt."
    ELSE
      !Actual number of fields is larger than expected
      msg = "          Zusätzliche Daten werden ignoriert."
    ENDIF
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(1709)
  !strings(1) = keyword
  !strings(2) = type of transformation that cannot be performed
  msg = TRIM(ADJUSTL(warnmsg))//" vom Schlüsselwort '"//TRIM(ADJUSTL(strings(1)))//"':"
  CALL DISPLAY_MSG(1,msg,logfile)
  IF( strings(2)=="remove atoms" ) THEN
    msg = "            Atome können nicht entfernt werden, da noch keine Atomposition gelesen wurde."
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(1799)
  msg = TRIM(ADJUSTL(warnmsg))//" Die Datendatei hatte ein unbekanntes Format. Atomsk hat versucht,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            atomare Daten daraus zu extrahieren, aber es ist möglicherweise falsch. Vorsichtig auftreten!"
  CALL DISPLAY_MSG(1,msg,logfile)
!
!1800-1899: FEHLER MESSAGES
CASE(1800)
  msg = TRIM(ADJUSTL(errmsg))//" Format der Eingabedatei konnte nicht bestimmt werden."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Moeglicherweise wird dieses Format von atomsk noch nicht unterstuetzt."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Die Datei wird uebersprungen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1801)
  !strings(1) = file name
  !reals(1) = line number
  msg = "X!X FEHLER beim Lesen aus der Datei: " &
      & //TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
  IF( NINT(reals(1))>0 ) THEN
    WRITE(temp,*) NINT(reals(1))
    msg = "          Error appears to be at line # "//TRIM(ADJUSTL(temp))
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(1802)
  !strings(1) = bad array
  msg = TRIM(ADJUSTL(errmsg))//" inkonsistente Feldgroesse in "//TRIM(strings(1))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1803)
  msg = TRIM(ADJUSTL(errmsg))//" Datenmenge der Hilfseigenschaften "// &
             & "stimmt nicht mit der Anzahl der Atome ueberein."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1804)
  msg = TRIM(ADJUSTL(errmsg))//" Unbekanntes oder nicht unterstütztes Format."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1805)
  !reals(1) = number of particles read
  !reals(2) = number of particles declared
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(errmsg))//" Anzahl der eingelesenen Atome ("//TRIM(ADJUSTL(temp))//") stimmt nicht mit"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            der festgelegten Anzahl ueberein ("//TRIM(ADJUSTL(temp2))//")."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1806)
  !reals(1) = index of atom
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" Nenner ist Null fuer die Koordinate von Atom #"// &
      & TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1807)
  msg = TRIM(ADJUSTL(errmsg))//" Datei ist nicht im ASCII Format, Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1808)
  msg = TRIM(ADJUSTL(errmsg))//" levcfg darf nicht groesser als 2 sein."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1809)
  msg = TRIM(ADJUSTL(errmsg))//" Kann Parameter der Superzelle nicht einlesen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1810)
  msg = TRIM(ADJUSTL(errmsg))//" Kann Anzahl der Atome nicht einlesen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1811)
  !reals(1) = index of atom
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" Atom index #"//TRIM(ADJUSTL(msg))// &
      & " ist groesser als Anzahl der Atome."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1812)
  !reals(1) = index of atom
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" Kann Eigenschaften von Atom #"//TRIM(ADJUSTL(msg))&
      & //" nicht einlesen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1813)
  !strings(1) = string where conversion failed
  msg = "X!X FEHLER beim Analysieren der Operation in '"// &
      & TRIM(strings(1))//"'."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1814)
  !strings(1) = name of file format
  !strings(2) = name of mode to use to read this file format
  IF( LEN_TRIM(strings(2))>0 ) THEN
    msg = TRIM(ADJUSTL(errmsg))//" Datei im Format "//TRIM(ADJUSTL(strings(1)))//  &
        & " könnten nur im Modus "//TRIM(ADJUSTL(strings(2)))//" gelesen werden."
    CALL DISPLAY_MSG(1,msg,logfile)
    msg = "    Verwendung:"
    CALL DISPLAY_MSG(1,msg,logfile)
    IF( LEN_TRIM(strings(3))>0 ) THEN
      msg = TRIM(ADJUSTL(strings(3)))
    ELSE
      msg = "<inputfile>"
    ENDIF
    msg = "      atomsk --one-in-all "//TRIM(ADJUSTL(msg))//" <format> [<options>]"
    CALL DISPLAY_MSG(1,msg,logfile)
  ELSE
    msg = TRIM(ADJUSTL(errmsg))//" Datei im Format "//TRIM(ADJUSTL(strings(1)))//" könnten nicht gelesen werden."
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(1815)
  msg = TRIM(ADJUSTL(errmsg))//" die Nachbarliste ist leer, wahrscheinlich weil Atome zu weit voneinander entfernt sind."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1816)
  msg = TRIM(ADJUSTL(errmsg))//" Eingabedatei hat ein nicht unterstütztes Binärformat."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!
!
!**************************
! 2000-2999: MESSAGES FOR OPTIONS
!
CASE(2050)
  msg = ">>> Richte Zellenvektoren aus..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2051)
  msg = "..> Zellenvektoren erfolgreich ausgerichtet."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2052)
  msg = ">>> Konvertiere in fraktionale Koordinaten..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2053)
  msg = "..> Koordinaten wurden reduziert."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2054)
  !strings(1) = atom species
  msg = ">>> Erstelle Schalen fuer "//TRIM(strings(1))//" Atome."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2055)
  !reals(1) = number of shells
  IF( NINT(reals(1))==0 ) THEN
    msg = "..> Keine Schalen hinzugefuegt."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 Schale wurde dem System hinzugefuegt."
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))// &
        & " Schalen wurden dem System hinzugefuegt."
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
  msg = ">>> Schneide "//TRIM(strings(1))//" "// &
    & TRIM(ADJUSTL(msg))//" A entlang "//TRIM(strings(2))//" ab."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2057)
  !reals(1) = NPcut, number of deleted atoms
  !reals(2) = number of atoms left
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 Atom wurde geloescht."
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" Atome wurden geloescht."
  ENDIF
  WRITE(temp,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(msg))//", "//TRIM(ADJUSTL(temp))//" Atome uebrig."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2058)
  !reals(1) = deformation
  !strings(1) = direction of deformation: x, y or z
  WRITE(msg,"(f16.3)") reals(1)*100.d0
  msg = ">>> Deformiere das System um "//TRIM(ADJUSTL(msg))// &
      & "% entlang "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2059)
  !reals(1) = Poisson coefficient
  IF(reals(1)==0.d0) THEN
    msg = "..> Uniaxiale Verzerrung."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSE
    WRITE(msg,"(f16.3)") reals(1)
    msg = "..> Uniaxiale Spannung, Poisson Verhaeltnis: "//TRIM(ADJUSTL(msg))
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
CASE(2060)
  msg = "..> Das System wurde erfolgreich deformiert."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2061)
  !strings(1) = disloctype: screw, edge, edge_add, edge_rm
  !strings(2) = direction of dislocline: x, y or z
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
  temp = TRIM(ADJUSTL(strings(1)))
  IF(temp(1:4)=="file" .OR. temp(1:5)=="array") THEN
    msg = ">>> In der Datei definierte Versetzungen einfügen: "//TRIM(ADJUSTL(strings(2)))
  ELSE
    IF(TRIM(temp)=="screw") THEN
      msg = ">>> Fuege eine Schraubenversetzung ein, entlang der Linie"
    ELSEIF(temp(1:4)=="edge") THEN
      msg = ">>> Fuege eine Stufenversetzung ein, entlang der Linie"
    ELSEIF(temp(1:5)=="mixed") THEN
      msg = ">>> Fuege einer gemischten Versetzung mit der Linie entlang"
    ELSEIF(temp(1:4)=="loop") THEN
      msg = ">>> Fuege einer Versetzungsschleife in eine Ebene normal zu"
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
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
  !
  IF(temp(1:4).NE."file" .AND. temp(1:5).NE."array") THEN
    !
    IF( reals(4)>0.1d0 ) THEN
      msg = "    wende anisotrope Elastizitaet an,"
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ENDIF
    !
    IF( TRIM(strings(1))=="edge_add" .OR. j>0 ) THEN
      WRITE(msg,"(a34)") "    durch Einfuegen einer atomaren Ebene,"
    ELSEIF( TRIM(strings(1))=="edge_rm" .OR. j<0 ) THEN
      WRITE(msg,"(a34)") "    durch Entfernen einer atomaren Ebene,"
    ELSEIF(TRIM(strings(1))=="edge" .OR. TRIM(strings(1))=="screw" .OR. TRIM(strings(1))=="mixed" .OR. j==0) THEN
      WRITE(msg,"(a41)") "    Erhaltung der Anzahl der Atome,"
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
        msg = "    Zentrum ("//TRIM(ADJUSTL(temp))//","//TRIM(ADJUSTL(temp2))//","// &
            & TRIM(ADJUSTL(temp3))//"); Radius "//TRIM(ADJUSTL(temp4))//" A; b="//TRIM(ADJUSTL(msg))
      ELSE
        WRITE(temp4,"(f16.3)") DABS(reals(8))
        msg = "    Zentrum ("//TRIM(ADJUSTL(temp))//","//TRIM(ADJUSTL(temp2))//","// &
            & TRIM(ADJUSTL(temp3))//"); Side "//TRIM(ADJUSTL(temp4))//" A; b="//TRIM(ADJUSTL(msg))
      ENDIF
    ELSE
      msg = "    b="//TRIM(ADJUSTL(msg))//" at ("// &
          & TRIM(ADJUSTL(temp))//","//TRIM(ADJUSTL(temp2))//")"
    ENDIF
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
CASE(2062)
  msg = "..> Bestimme Loesungen der anisotropen Elastizitaetsgleichungen..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2063)
  !reals(1) = number of inserted atoms
  IF( NINT(reals(1)) < 0 ) THEN
    WRITE(msg,*) NINT(ABS(reals(1)))
    msg = "..> "//TRIM(ADJUSTL(msg))//" Atome wurden entfernt."
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" Atome wurden eingefuegt."
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
  msg = "..> Superzelle wurde  erlang "//TRIM(ADJUSTL(temp))//" um "//TRIM(ADJUSTL(strings(1)))//" verändert."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2065)
  !reals(1) = number of dislocations inserted (default 1)
  IF( NINT(reals(1)) <= 0 ) THEN
    msg = "..> Kein Versetzung wurde erzeugt."
  ELSEIF( NINT(reals(1)) == 1 ) THEN
    msg = "..> Ein Versetzung wurde erfolgreich erzeugt."
  ELSE
    WRITE(temp,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(temp))//" Versetzungen wurden erfolgreich erzeugt."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2066)
  !reals(1) = number of repetitions along X
  !reals(2) = number of repetitions along Y
  !reals(3) = number of repetitions along Z
  WRITE(temp,*) NINT( reals(1) )
  WRITE(temp2,*) NINT( reals(2) )
  WRITE(msg,*) NINT( reals(3) )
  msg = ">>> Vervielfaeltige das System um: "//TRIM(ADJUSTL(temp))//" x "// &
    & TRIM(ADJUSTL(temp2))//" x "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2067)
  !reals(1) = new number of particles
  WRITE(msg,*) NINT( reals(1) )
  msg = "..> Neue Teilchenzahl: "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2068)
  !reals(1) = new number of atoms
  WRITE(temp,*) NINT( reals(1) )
  msg = "..> Das System wurde erfolgreich vervielfaeltigt ("//TRIM(ADJUSTL(temp))//" Atome)."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2069)
  !strings(1) = "all" or property to be removed
  IF(strings(1)=="all") THEN
    WRITE(msg,'(a)') ">>> Entferne alle Hilfseigenschaften."
  ELSE
    msg = ">>> Entferne Hilfseigenschaft: "//TRIM(strings(1))
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2070)
  !strings(1) = "all" or property to be removed
  IF(strings(1)=="all") THEN
    msg = "..> Hilfseigenschaften wurde erfolgreich entfernt."
  ELSE
    msg = "..> Hilfseigenschaft wurde erfolgreich entfernt."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2071)
  msg = ">>> Drehen des Systems, um die Kristallorientierung zu ändern..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2072)
  msg = "..> System erfolgreich gedreht. Neuausrichtung:"// &
      & TRIM(ADJUSTL(strings(1)))//" "//TRIM(ADJUSTL(strings(2)))//" "//TRIM(ADJUSTL(strings(3)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2073)
  !strings(1) = file containing properties
  msg = ">>> Lese Systemeigenschaften aus "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2074)
  !strings(1) = property that was read
  msg = "..> Eigenschaft '"//TRIM(ADJUSTL(strings(1)))//"' wurde eingelesen."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2075)
  msg = "..> Einlesen der Systemeigenschaften abgeschlossen."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2076)
  !msg = ">>> Waehle alle Atome aus."
  !CALL DISPLAY_MSG(verbosity,msg,logfile)
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
    msg = ">>> Abwaehle"
  ELSE
    msg = ">>> Waehle"
  ENDIF
  IF( strings(1)=="all" .OR. strings(1)=="none" ) THEN
    msg = ">>> Waehle alle Atome aus."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="invert" ) THEN
    msg = ">>> Invertiere die Auswahl der Atome..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="index" ) THEN
    IF( strings(2)=="list " ) THEN
      msg = TRIM(ADJUSTL(msg))//" eine Liste von Atomen aus."
    ELSEIF( strings(2)=="range" ) THEN
      WRITE(temp,*) NINT(reals(1))
      WRITE(temp2,*) NINT(reals(2))
      msg = TRIM(ADJUSTL(msg))//" Atom #"//TRIM(ADJUSTL(temp))//" bis "//TRIM(ADJUSTL(temp2))//" aus."
    ELSE
      WRITE(temp,*) NINT(reals(1))
      msg = TRIM(ADJUSTL(msg))//" Atom #"//TRIM(ADJUSTL(temp))//" aus."
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
    msg = TRIM(ADJUSTL(msg))//" Atome "//TRIM(strings(1))//" "//TRIM(ADJUSTL(temp))// &
        & " A entlang der "//TRIM(strings(3))//" Achse aus."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="in" .OR. strings(1)=="out" ) THEN
    IF(strings(1)=="in") THEN
      temp = "innerhalb der"
    ELSE
      temp = "ausserhalb der"
    ENDIF
    msg = TRIM(ADJUSTL(msg))//" die Atome "//TRIM(temp)//" "//TRIM(strings(2))//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    IF(TRIM(strings(2))=="box") THEN
      msg = "..> Zellenbegrenzung: ("
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
      msg = "..> Mittelpunkt: ("//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(2)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(3)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))//")"
      WRITE(temp,"(f16.3)") reals(4)
      msg = TRIM(msg)//"; Radius: "//TRIM(ADJUSTL(temp))//" A."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ELSEIF(TRIM(strings(2))=="cylinder") THEN
      msg = "..> Achse entlang "//TRIM(strings(3))//"."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
      WRITE(temp,"(f16.3)") reals(1)
      msg = "..> Mittelpunkt: ("//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(2)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))//")"
      WRITE(temp,"(f16.3)") reals(4)
      msg = TRIM(msg)//"; Radius: "//TRIM(ADJUSTL(temp))//" A."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ELSEIF(TRIM(strings(2))=="cone") THEN
      WRITE(temp,"(f16.3)") reals(1)
      msg = "..> Achse entlang "//TRIM(strings(3))//". Tipp bei: ("//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(2)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(3)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))//"). Öffnungswinkel: "
      WRITE(temp,"(f16.3)") reals(4)
      msg = TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(temp))//"°."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ELSEIF(TRIM(strings(2))=="torus") THEN
      WRITE(temp,"(f16.3)") reals(1)
      msg = "..> Achse entlang "//TRIM(strings(3))//". Mittelpunkt: ("//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(2)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(3)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))//")."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
      WRITE(temp,"(f16.3)") reals(4)
      WRITE(temp2,"(f16.3)") reals(5)
      msg = "..> Main Radius: "//TRIM(ADJUSTL(temp))//" A; Secündar Radius: "//TRIM(ADJUSTL(temp2))//" A."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ENDIF
  ELSEIF( strings(1)=="prop" ) THEN
    msg = TRIM(ADJUSTL(msg))//" Atome mit "//TRIM(strings(2))
    WRITE(temp,"(f16.3)") reals(1)
    IF( strings(3)=="min" .OR. strings(3)=="max" ) THEN
      msg = TRIM(ADJUSTL(msg))//" "//TRIM(strings(3))//"imum."
    ELSEIF( reals(4)>2.d0 ) THEN
      WRITE(temp2,"(f16.3)") reals(2)
      msg = TRIM(ADJUSTL(msg))//" zwischen "//TRIM(ADJUSTL(temp))//" und "//TRIM(ADJUSTL(temp2))
    ELSE
      msg = TRIM(ADJUSTL(msg))//" gleich "//TRIM(ADJUSTL(temp))
    ENDIF
    IF(strings(4)=="among") msg = TRIM(ADJUSTL(msg))//" unter ausgewählten Atomen"
    msg = TRIM(ADJUSTL(msg))//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="random" .OR. strings(1)=="rand" ) THEN
    WRITE(temp,*) NINT(reals(1))
    msg = ">>> Zufaellige Auswahl von "//TRIM(ADJUSTL(temp))//" Atomen"
    IF( strings(2).NE."any" .AND. strings(2).NE."all" ) THEN
      msg = TRIM(msg)//" von "//TRIM(strings(2))
    ENDIF
    IF(strings(4)=="among") msg = TRIM(ADJUSTL(msg))//" unter ausgewählten Atomen"
    msg = TRIM(msg)//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="random%" ) THEN
    WRITE(temp,'(f16.3)') reals(1)*100.d0
    msg = ">>> Zufaellige Auswahl von "//TRIM(ADJUSTL(temp))//"% der Atome"
    IF( strings(2).NE."any" .AND. strings(2).NE."all" ) THEN
      msg = TRIM(msg)//" von "//TRIM(strings(2))
    ENDIF
    IF(strings(4)=="among") msg = TRIM(ADJUSTL(msg))//" unter ausgewählten Atomen"
    msg = TRIM(msg)//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="modulo" .OR. strings(1)=="mod" ) THEN
    !reals(1) = index of an atom
    !reals(2) = modulo
    WRITE(temp,*) NINT(reals(1))
    WRITE(temp2,*) NINT(reals(2))
    msg = TRIM(ADJUSTL(msg))//" Atomen mit einem Index von "//TRIM(ADJUSTL(temp))//" Modulo "//TRIM(ADJUSTL(temp2))//"."
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
      msg = TRIM(ADJUSTL(msg))//" erste naechste"//TRIM(temp2)//" Nachbarn"
    ELSEIF( reals(1)>0.d0 .AND. DBLE(NINT(reals(1)))-reals(1)<1.d-12 ) THEN
      WRITE(temp,*) NINT(reals(1))
      IF( NINT(reals(1))==1 ) THEN
        msg = TRIM(ADJUSTL(msg))//" der erste"//TRIM(temp2)//" Nachbar"
      ELSE
        msg = TRIM(ADJUSTL(msg))//" die "//TRIM(ADJUSTL(temp))//" naehsten"//TRIM(temp2)//" Nachbarn"
      ENDIF
    ELSE
      WRITE(temp,'(f16.3)') DABS(reals(1))
      msg = TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(temp2))//" Nachbarn innerhalb eines Radius von "//TRIM(ADJUSTL(temp))//" A"
    ENDIF
    WRITE(temp,*) NINT(reals(2))
    msg = TRIM(ADJUSTL(msg))//" von Atom #"//TRIM(ADJUSTL(temp))//"..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="stl" ) THEN
    msg = TRIM(ADJUSTL(msg))//" atoms inside the 3-D shape from the STL file: "//TRIM(ADJUSTL(strings(2)))
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSE
    !Last case: strings(1) should be an atom species
    msg = TRIM(ADJUSTL(msg))//" alle "//TRIM(ADJUSTL(strings(1)))//" Atome..."
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
        msg = TRIM(ADJUSTL(msg))//" 1 Atome wurde zur Auswahl hinzugefügt"
      ELSEIF( NINT(reals(2))>1 ) THEN
        WRITE(temp,*) NINT( reals(2) )
        msg = TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(temp))//" Atome wurden zur Auswahl hinzugefügt"
      ENDIF
    ENDIF
    IF( NINT(reals(3))>0 ) THEN
      IF( NINT(reals(2))>0 ) THEN
        msg = TRIM(ADJUSTL(msg))//","
      ENDIF
      IF( NINT(reals(3))==1 ) THEN
        msg = TRIM(ADJUSTL(msg))//" 1 Atom wurde zur Auswahl entfernt."
      ELSEIF( NINT(reals(3))>1 ) THEN
        WRITE(temp,*) NINT( reals(3) )
        msg = TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(temp))//" Atome wurden zur Auswahl entfernt."
      ENDIF
    ENDIF
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 Atom"//TRIM(temp3)//" ist ausgewaehlt."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(temp,*) NINT( reals(1) )
    msg = "..> "//TRIM(ADJUSTL(temp))//" Atome sind ausgewaehlt"//TRIM(temp3)//"."
  ELSE
    msg = "..> Kein Atom ist ausgewaehlt."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2079)
  !reals(1) = cutoff radius for removing atoms
  WRITE(msg,"(f16.3)") reals(1)
  msg = ">>> Entferne Atome die naeher als "//TRIM(ADJUSTL(msg))//" A sind."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2080)
  !strings(1) = species of removed atom(s) (may be empty or contain anything)
  !reals(1) = number of atoms removed
  !reals(2) = number of atoms left
  species = strings(1)
  CALL ATOMNUMBER(species,tempreal)
  IF(NINT(reals(1))==0) THEN
    msg = "..> Kein Atom wurde entfernt"
  ELSEIF(NINT(reals(1))==1) THEN
    IF(tempreal>0.5d0) THEN
      msg = "..> 1 Atom von "//TRIM(ADJUSTL(species))//" wurde entfernt"
    ELSE
      msg = "..> 1 Atom wurde entfernt"
    ENDIF
  ELSE
    WRITE(msg,*) NINT(reals(1))
    IF(tempreal>0.5d0) THEN
      msg = "..> "//TRIM(ADJUSTL(msg))//" Atome von "// &
          & TRIM(ADJUSTL(species))//" wurden entfernt"
    ELSE
      msg = "..> "//TRIM(ADJUSTL(msg))//" Atome wurden entfernt"
    ENDIF
  ENDIF
  WRITE(temp,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(msg))//", "//TRIM(ADJUSTL(temp))//" Atome uebrig."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2081)
  !strings(1) = rotation axis: x, y or z
  !reals(1) = rotation angle in degrees
  WRITE(msg,"(f16.2)") reals(1)
  msg = ">>> Rotatiere das System um "//TRIM(ADJUSTL(msg))//"°"
  IF( SCAN(strings(1),'xXyYzZ[]')>0 ) THEN
    msg = TRIM(msg)//" um die Achse "//TRIM(strings(1))//"."
  ELSE
    msg = TRIM(msg)//" um die Vektor ("//TRIM(strings(1))//")."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2082)
  msg = "..> System erfolgreich rotiert."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2083)
  !strings(1) = axis that is tilted to produce shear: x, y or z
  !strings(2) = direction of tilt: x, y or z
  !reals(1) = magnitude of tilt in Angstroms
  WRITE(temp,"(f16.3)") reals(1)
  msg = ">>> Neige die "//TRIM(strings(1))//" Achse um "// &
      & TRIM(ADJUSTL(temp))//" A entlang "//TRIM(strings(2))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2084)
  !reals(1) = shear strain in %
  msg = "..> System erfolgreich geschert."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,"(f16.3)") reals(1)
  msg = "..> Wende Schaerspannung = "//TRIM(ADJUSTL(msg))//" % an."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2085)
  !strings(1) = shift direction: above or below
  !strings(2) = shift axis
  !reals(1) = distance above axis
  !reals(2), reals(3), reals(4) = shift vector
  WRITE(temp,"(f16.3)") reals(2)
  WRITE(temp2,"(f16.3)") reals(3)
  WRITE(temp3,"(f16.3)") reals(4)
  IF( DABS(reals(1))<1.d12 ) THEN
    WRITE(msg,"(3f16.3)") reals(1)
  ELSEIF( reals(1)<-1.d12 ) THEN
    WRITE(msg,"(a4)") "-INF"
  ELSEIF( reals(1)>1.d12 ) THEN
    WRITE(msg,"(a4)") "+INF"
  ENDIF
  msg = ">>> Verschiebe Atome "//TRIM(strings(1))//" "//TRIM(ADJUSTL(msg))//  &
        & "A entlang "//TRIM(strings(2))//" by ("//TRIM(ADJUSTL(temp))//","//     &
        & TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))//")"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2086)
  !reals(1) = number of atoms that were shifted
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 Atom wurde verschoben."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" Atome wurden verschoben."
  ELSE
    msg = "..> Kein Atom wurde verschoben."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2087)
  !strings(1) = property to be sorted
  !strings(2) = sort order: up, down, pack, random
  IF( strings(2)=="random" ) THEN
    msg = ">>> Mischen der Atomliste in zufälliger Reihenfolge..."
  ELSEIF( strings(2)=="reverse" ) THEN
    msg = ">>> Atomliste umkehren..."
  ELSE
    IF(strings(1)=="s") THEN
      temp = "atomare Masse"
    ELSEIF( strings(1)=="x" .OR. strings(1)=="y" .OR. strings(1)=="z" .OR. &
          & strings(1)=="X" .OR. strings(1)=="Y" .OR. strings(1)=="Z"  ) THEN
      temp = TRIM(strings(1))//" Koordinate"
    ELSE
      temp = TRIM(strings(1))
    ENDIF
    IF(strings(2)=="up") temp = "anwachsend "//TRIM(temp)
    IF(strings(2)=="down") temp = "abfallend "//TRIM(temp)
    IF(strings(2)=="up" .OR. strings(2)=="down") THEN
      msg = ">>> Sortiere Atome nach "//TRIM(ADJUSTL(temp))//"."
    ELSE
      msg = ">>> Fasse Atome nach "//TRIM(ADJUSTL(temp))//" zusammen."
    ENDIF
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2088)
  msg = "..> Atome erfolgreich sortiert."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2089)
  !strings(1) = first atomic species
  !strings(2) = second atomic species
  msg = ">>> Ersetze "//TRIM(strings(1))//" Atoms durch "//TRIM(strings(2))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2090)
  !reals(1) = number of substituted atoms
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 Atom ersetzt."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" Atome ersetzt."
  ELSE
    msg = "..> Kein Atom ersetzt."
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
    msg = ">>> Konvertiere "//TRIM(strings(1))//" von "//TRIM(temp)//&
        & " nach "//TRIM(temp2)//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSE
    WRITE(temp,'(f16.3)') reals(1)
    msg = ">>> Neuskalierung die Eigenschaft "//TRIM(strings(1))//" um den Faktor "//TRIM(ADJUSTL(temp))//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
CASE(2092)
  !strings(1) = what was converted
  IF( strings(1)=="velocities" .OR. strings(1)=="coordinates" ) THEN
    msg = "..> "//TRIM(ADJUSTL(strings(1)))//" konvertiert."
  ELSE
    msg = "..> "//TRIM(ADJUSTL(strings(1)))//" neu skaliert."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2093)
  msg = ">>> Periodisches umbrechen der Atome..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2094)
  !reals(1) = number of atoms wrapped
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 Atom wurde periodisch verschoben."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" Atome wurden periodisch verschoben."
  ELSE
    msg = "..> Kein Atom wurde periodisch verschoben."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2095)
  !reals(1) = number of box vectors that were modified
  !reals(2) = index of first box vector that was modified
  !reals(3) = index of second box vector that was modified
  IF( NINT(reals(2))==1 ) THEN
    temp = "Der erste"
  ELSEIF( NINT(reals(2))==2 ) THEN
    temp = "Der zweite"
  ELSEIF( NINT(reals(2))==3 ) THEN
    temp = "Der dritte"
  ENDIF
  IF( NINT(reals(1))==0 ) THEN
    msg = "..> Zellvektoren wurden nicht modifiziert."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> "//TRIM(ADJUSTL(temp))//" Zellvektoren wurden modifiziert."
  ELSEIF( NINT(reals(1))==2 ) THEN
    IF( NINT(reals(3))==2 ) THEN
      temp2 = "zweite"
    ELSEIF( NINT(reals(3))==3 ) THEN
      temp2 = "dritte"
    ENDIF
    msg = "..> "//TRIM(ADJUSTL(temp))//" und "//TRIM(ADJUSTL(temp2))//" Zellvektoren wurden modifiziert."
  ELSE
    msg = "..> Zellvektoren wurden modifiziert."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2096)
  msg = "..> Die Loesungen wurden gefunden."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2097)
  !strings(1) = axis along which coordinate is fixed: x,y,z,all
  !strings(2) = "above" or "below"
  !strings(3) = axis 
  !reals(1) = fix distance in angstroms
  IF( DABS(reals(1))<1.d12 ) THEN
    WRITE(msg,"(3f16.3)") reals(1)
  ELSEIF( reals(1)<-1.d12 ) THEN
    WRITE(msg,"(a4)") "-INF"
  ELSEIF( reals(1)>1.d12 ) THEN
    WRITE(msg,"(a4)") "+INF"
  ENDIF
  msg = ">>> Fixiere "//TRIM(strings(1))//" Atomkoordinaten "//&
      & TRIM(strings(2))//" "//TRIM(ADJUSTL(msg))//"A along "//TRIM(strings(3))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2098)
  !reals(1) = number of atoms that were fixed
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 Atom fixiert."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" Atome fixiert."
  ELSE
    msg = "..> Kein Atom fixiert."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2099)
  msg = "..> Elastischer Tensor wurde rotiert."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2100)
  !reals(1) = anisotropy ratio A = 2*C44 / (C11-C12)
  !reals(2) = anisotropy factor H = 2*C44 + C12 - C11
  WRITE(msg,'(f16.3)') reals(1)
  msg = "..> Anisotropieverhaeltnis: A = 2*C44 / (C11-C12) = "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,'(e16.3)') reals(2)
  msg = "..> Anisotropiefaktor:      H = 2*C44 + C12 - C11 = "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2101)
  !strings(1) = formula of energy factor, e.g. "Kb²"
  !reals(1) = energy factor
  msg = "..> Versetzungsspannungen wurden berechnet."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,'(f24.8)') reals(1)
  msg = "..> Energie Faktor "//TRIM(ADJUSTL(strings(1)))//&
      & " = "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2102)
  !strings(1) = empty or atom species to remove or "SEL"
  !reals(1) = 0 or atom index to remove
  IF( strings(1)=="SEL" ) THEN
    msg = ">>> Entferne die ausgewaehlten Atome aus dem System."
  ELSE
    IF( NINT(reals(1)).NE.0 ) THEN
      WRITE(msg,*) NINT(reals(1))
      msg = ">>> Entferne Atom #"//TRIM(ADJUSTL(msg))//"."
    ELSE
      msg = ">>> Entferne alle "//TRIM(strings(1))//" Atome aus dem System."
    ENDIF
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2103)
  msg = ">>> Reduziere die Neigung der Superzellen Vektoren..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2104)
  !reals(1) = number of components that were unskewed
  IF( NINT(reals(1))==0 ) THEN
    msg = "..> Kein Superzellen Vektor wurde veraendert."
  ELSE
    msg = "..> Neigung wurde aus der Superzelle entfernt."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2105)
  msg = ">>> Binde ionische Schalen and die entsprechenden Kerne..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2106)
  !reals(1) = number of shells that were re-assigned
  IF( NINT(reals(1))<=0 ) THEN
    msg = "..> Keine Schale wurde neu gebunden."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 Schale wurde an ihren Kern gebunden."
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" Schalen wurden an ihre Kerne gebunden."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2107)
  !strings(1) = atom species to be removed
  !reals(1) = number of atoms to be randomly removed
  WRITE(msg,*) NINT(reals(1))
  IF( strings(1)=="any" .OR. strings(1)=="all" ) THEN
    msg = ">>> Entferne "//TRIM(ADJUSTL(msg))//" zufaellige Atome aus dem System."
  ELSE
    msg = ">>> Entferne "//TRIM(ADJUSTL(msg))//" zufaellige "//TRIM(strings(1))//&
        & " Atome aus dem System."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2108)
  !strings(1) = crack mode (I, II or III)
  !strings(2) = direction of crack line (X, Y or Z)
  !reals(1) = pos1
  !reals(2) = pos2
  WRITE(msg,"(f16.3)") reals(1)
  WRITE(temp,"(f16.3)") reals(2)
  msg = ">>> Erzeuge einen Typ "//TRIM(ADJUSTL(strings(1)))//" Riss entlang "// &
      & TRIM(ADJUSTL(strings(2)))//" bei ("//TRIM(ADJUSTL(msg))//","//         &
      & TRIM(ADJUSTL(temp))//")."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2109)
  msg = "..> Riss erfolgreich erzeugt."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2110)
  !msg = ">>> Invertiere die Atomauswahl..."
  !CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2111)
  !reals(1) = target temperature
  WRITE(msg,'(f24.3)') reals(1)
  msg = ">>> Weise Geschwindigkeiten entsprechend einer Maxwell-Boltzmann Verteilung bei "// &
      & TRIM(ADJUSTL(msg))//" K zu."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2112)
  msg = "..> Gesweindigkeiten der Atome erfolgrich gesetzt."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2113)
  !strings(1) = file name
  msg = "..> Geschwindigkeitsverteilung gespeichert in der Datei: "// &
      & TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2114)
  !strings(1) = "xyz" if user provided max.disp. along X,Y,Z
  !reals(1) = max. displacement of an atom along X
  !reals(2) = max. displacement of an atom along Y
  !reals(3) = max. displacement of an atom along Z
  msg = ">>> Zufaellige Verschiebung der Atome,"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  IF( strings(1)=="xyz" ) THEN
    !User provided different values along each direction
    WRITE(temp,'(f24.3)') reals(1)
    WRITE(temp2,'(f24.3)') reals(2)
    WRITE(temp3,'(f24.3)') reals(3)
    msg = "..> max. Distanz: dx="//TRIM(ADJUSTL(temp))//", dy=" &
        & //TRIM(ADJUSTL(temp2))//", dz="//TRIM(ADJUSTL(temp3))//"."
  ELSE
    !User provided norm of max.disp.
    WRITE(msg,'(f24.3)') reals(1)
    msg = "..> maximum Distanz: "//TRIM(ADJUSTL(msg))//" A."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2115)
  !reals(1) = number of atoms that were displaced
  WRITE(temp,*) NINT(reals(1))
  msg = "..> "//TRIM(ADJUSTL(temp))//" Atome wurden aus ihrer Gleichgewichtslage verschoben."
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
    msg = ">>> Fuege ein "//TRIM(ADJUSTL(strings(1)))//" Atom bei ("// &
        & TRIM(ADJUSTL(temp))//","//     &
        & TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))//") ein."
  CASE("relative","rel")
    WRITE(temp,"(f16.3)") reals(2)
    WRITE(temp2,"(f16.3)") reals(3)
    WRITE(temp3,"(f16.3)") reals(4)
    WRITE(temp4,*) NINT(reals(1))
    msg = ">>> Fuege ein "//TRIM(ADJUSTL(strings(1)))//" Atom bei ("//TRIM(ADJUSTL(temp))//","//     &
        & TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))//") relativ zum Atom #"//TRIM(ADJUSTL(temp4))//"."
  CASE("near","NEAR")
    WRITE(temp,*) NINT(reals(1))
    msg = ">>> Fuege ein "//TRIM(ADJUSTL(strings(1)))// &
        & " Atom in der Naehe von Atom #"//TRIM(ADJUSTL(temp))//" ein."
  CASE("random","RANDOM","rand","RAND")
    WRITE(temp,*) NINT(reals(1))
    msg = ">>> Fuege "//TRIM(ADJUSTL(temp))//" "// &
        & TRIM(ADJUSTL(strings(1)))// &
        & " Atome an zufaelligen Positionen ein..."
  END SELECT
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2117)
  !reals(1) = number of atoms added
  !reals(2) = number of atoms in the system
  IF(NINT(reals(1))==0) THEN
    msg = "..> Kein Atom hinzugefuegt"
  ELSEIF(NINT(reals(1))==1) THEN
    msg = "..> 1 Atom hinzugefuegt"
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" Atome hinzugefuegt"
  ENDIF
  WRITE(temp,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(msg))//", "//TRIM(ADJUSTL(temp))// &
      & " Atome im System."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2118)
  !reals(1) = index of atom to center; if <=0 then center of mass
  IF( NINT(reals(1))<=0 ) THEN
    msg = ">>> Bringe das Massezentrum in die Mitte der Zelle..."
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = ">>> Verschiebe das System so, dass Atom #"// &
        & TRIM(ADJUSTL(msg))//" in der Mitte der Zelle liegt..."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2119)
  !reals(1) = X component of shift vector
  !reals(2) = Y component of shift vector
  !reals(3) = Z component of shift vector
  WRITE(temp,'(f16.3)') reals(1)
  WRITE(temp2,'(f16.3)') reals(2)
  WRITE(temp3,'(f16.3)') reals(3)
  msg = "..> System wurde zentriert, Verschiebungsvektor: ("// &
      & TRIM(ADJUSTL(temp))//","//TRIM(ADJUSTL(temp2))//","// &
      & TRIM(ADJUSTL(temp3))//")."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2120)
  !strings(1) = direction normal to mirror plane
  !reals(1) = distance between mirror plane and origin of coordinates
  WRITE(temp,'(f16.3)') reals(1)
  msg = ">>> Spiegelung an der Ebene bei "//TRIM(ADJUSTL(temp))// &
      & " A entlang "//TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2121)
  msg = "..> System erfolgreich gespiegelt."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2122)
  !strings(1) = species on which shells are removed
  msg = ">>> Entferne Schalen der "//TRIM(ADJUSTL(strings(1)))//" Ionen..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2123)
  !reals(1) = number of shells removed
  WRITE(temp,*) NINT(reals(1))
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 Schale wurde entfernt."
  ELSE
    msg = "..> "//TRIM(ADJUSTL(temp))//" Schalen wurden entfernt."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2124)
  !strings(1) = stress component or name of file
  !reals(1) = value of stress
  WRITE(temp,'(f16.3)') reals(1)
  SELECT CASE(strings(1))
  CASE('x','X','xx','XX','y','Y','yy','YY','z','Z','zz','ZZ', &
     & 'xy','XY','yx','YX','zx','ZX','xz','XZ','zy','ZY','yz','YZ')
    temp2 = ADJUSTL(strings(1))
    IF( LEN_TRIM(temp2)<=1 ) THEN
      temp2 = TRIM(ADJUSTL(temp2))//TRIM(ADJUSTL(temp2))
    ENDIF
    msg = ">>> Wende Spannung von σ_"//TRIM(ADJUSTL(temp2))//" = "// &
        & TRIM(ADJUSTL(temp))//" GPa an."
  CASE('p','P')
    msg = ">>> Wende isostatischen Druck von "//TRIM(ADJUSTL(temp))// &
        & " GPa an."
  CASE DEFAULT
    msg = ">>> Wende den Spannungszustand aus der Datei "// &
        & TRIM(ADJUSTL(temp))//" an."
  END SELECT
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2125)
  !reals(1) = type of item to swap: 0=Cartesian axes, 1=atom id, 2= atom species, 3=aux.prop.
  !strings(1) = Cartesian axis, or integer
  !strings(2) = same type as strings(1)
  SELECT CASE(NINT(reals(1)))
  CASE(0)
    msg = ">>> Tauschen die Achsen "//TRIM(ADJUSTL(strings(1)))//" und "//TRIM(ADJUSTL(strings(2)))//"."
  CASE(1)
    msg = ">>> Tauschen die Atomen #"//TRIM(ADJUSTL(strings(1)))//" und #"//TRIM(ADJUSTL(strings(2)))//"."
  CASE(2)
    msg = ">>> Tauschen die "//TRIM(ADJUSTL(strings(1)))//" und "//TRIM(ADJUSTL(strings(2)))//" Atomen."
  CASE(3)
    msg = ">>> Tauschen Nebeneigenschaften '"//TRIM(ADJUSTL(strings(1)))//"' und '"//TRIM(ADJUSTL(strings(2)))//"'."
  END SELECT
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2126)
  !reals(1) = number of atoms that were swapped
  IF( NINT(reals(1))==0 ) THEN
    msg = "..> Kein solches Atom wurde im System gefunden."
  ELSEIF( NINT(reals(1))>0 ) THEN
    WRITE(temp,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(temp))//" Atomarten wurden geändert."
  ELSE
    msg = "..> Atomen wurden vertauschen."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2127)
  !strings(1) = roll axis: x, y or z
  !reals(1) = roll angle in degrees
  WRITE(msg,"(f16.2)") reals(1)
  msg = ">>> Rollen die "//TRIM(ADJUSTL(strings(1)))//" Richtung um " // &
      & TRIM(ADJUSTL(msg))//"° um die "//TRIM(strings(2))//" Achse."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2128)
  msg = "..> System wurde erfolgreich gerollt."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2129)
  !strings(1) = torsion axis: x, y or z
  !reals(1) = torsion angle in degrees
  WRITE(msg,"(f16.2)") reals(1)
  msg = ">>> Torsion von "//TRIM(ADJUSTL(msg)) &
      & //"° um die "//TRIM(strings(1))//" Achse anwenden..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2130)
  msg = "..> Torsion wurde erfolgreich angewendet."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2131)
  !strings(1) = space group name or number
  msg = ">>> Anwenden von Symmetrieoperationen der Raumgruppe: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2132)
  !reals(1) = space group number
  WRITE(temp,"(i5)") NINT(reals(1))
  msg = "..> Raumgruppe nummer: "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2133)
  !reals(1) = new number of atoms
  WRITE(temp,"(i16)") NINT(reals(1))
  msg = "..> Symmetrieoperationen wurden erfolgreich angewendet, neue Anzahl von Atomen: "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2140)
  !reals(1) = min. separation distance
  !reals(2) = separate atoms by this amount
  WRITE(temp,"(f16.2)") reals(1)
  WRITE(temp2,"(f16.2)") reals(2)
  msg = ">>> Trennen von Atomen näher als "//TRIM(ADJUSTL(temp))//" A  in einer Entfernung von "//TRIM(ADJUSTL(temp2))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2141)
  !reals(1) = number of pairs of atoms that were separated
  IF( reals(1) < 1.d-3 ) THEN
    msg = "..> Keine Atome wurden getrennt."
  ELSEIF( NINT(reals(1)) == 1 ) THEN
    msg = "..> Ein Paar Atome wurde getrennt."
  ELSE
    WRITE(temp,"(i16)") NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(temp))//" Paare von Atomen wurden getrennt."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2142)
  !reals(1) = number of triangles in STL file
  WRITE(temp,*) NINT(reals(1))
  msg = "..> STL Datei wurde erfolgreich gelesen ("//TRIM(ADJUSTL(temp))//" Dreiecke)."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2143)
  msg = ">>> Konvertieren des Systems in eine orthorhombische Zelle..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2144)
  msg = "..> Es wurde eine geeignete Zelle gefunden, die nun mit Atomen gefüllt wurde..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(temp,*) NINT(reals(1))
  msg = "    (geschätzte neue Anzahl von Atomen: "//TRIM(ADJUSTL(temp))//")"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2145)
  WRITE(temp,*) NINT(reals(1))
  msg = "..> Die Zelle ist jetzt orthorhombisch ("//TRIM(ADJUSTL(temp))//" Atome)."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2146)
  msg = "..> Verschiebung von Atomen anwenden..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2147)
  !strings(1) = name of property that is rounded off
  msg = ">>> Rundung der Werte"
  IF( strings(1)=="AUX" ) THEN
    msg = TRIM(ADJUSTL(msg))//" der Hilfseigenschaften..."
  ELSEIF( strings(1)=="XYZ" ) THEN
    msg = TRIM(ADJUSTL(msg))//" der Atomkoordinaten..."
  ELSEIF( strings(1)=="X" .OR. strings(1)=="x" ) THEN
    msg = TRIM(ADJUSTL(msg))//" von X-koordinaten..."
  ELSEIF( strings(1)=="Y" .OR. strings(1)=="y" ) THEN
    msg = TRIM(ADJUSTL(msg))//" von Y-koordinaten..."
  ELSEIF( strings(1)=="Z" .OR. strings(1)=="z" ) THEN
    msg = TRIM(ADJUSTL(msg))//" von Z-koordinaten..."
  ELSE
    msg = TRIM(ADJUSTL(msg))//" der Eigenschaft '"//TRIM(ADJUSTL(strings(1)))//"'..."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2148)
  !reals(1) = number of values that were rounded off
  WRITE(temp,*) NINT(reals(1))
  msg = "..> Fertig, "//TRIM(ADJUSTL(temp))//" Werte wurden gerundet."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2149)
  !strings(1) = direction to reduce
  IF( strings(1)=='p' .OR. strings(1)=='P' ) THEN
    msg = ">>> Reduzierung der Systems in eine primitive Zelle..."
  ELSE
    msg = ">>> Reduzierung der Systemgröße unter Beibehaltung der Periodizität..."
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
      temp = "erste"
      j=j+1
    ENDIF
    IF( NINT(reals(2))==1 ) THEN
      IF( NINT(reals(1))==1 ) temp = TRIM(ADJUSTL(temp))//" und der"
      temp = TRIM(ADJUSTL(temp))//" zweite"
      j=j+1
    ENDIF
    IF( NINT(reals(3))==1 ) THEN
      IF( NINT(reals(1))==1 ) temp = TRIM(ADJUSTL(temp))//" und der"
      temp = TRIM(ADJUSTL(temp))//" dritte"
      j=j+1
    ENDIF
    IF(j==1) THEN
      temp2 = "Zellvektor wurde reduziert"
    ELSE
      temp2 = "Zellvektoren wurden reduziert"
    ENDIF
    msg = "..> Der "//TRIM(ADJUSTL(temp))//" "//TRIM(ADJUSTL(temp2))
  ELSE
    msg = "..> Die Zelle wurde reduziert"
  ENDIF
  WRITE(temp,*) NINT(reals(4))
  msg = TRIM(ADJUSTL(msg))//" ("//TRIM(ADJUSTL(temp))//" Atomen)."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2151)
  !strings(1) = operation performed on cell vector (add, rm, set)
  !strings(2) = component of cell vector
  !reals(1) = distance added (or removed)
  WRITE(temp2,'(f16.3)') reals(1)
  IF( strings(2)=="H1" ) THEN
    temp3 = "zum ersten Zellvektor"
  ELSEIF( strings(2)=="H2" ) THEN
    temp3 = "zum zweiten Zellvektor"
  ELSEIF( strings(2)=="H3" ) THEN
    temp3 = "zum dritten Zellvektor"
  ELSEIF( strings(2)=="x" .OR. strings(2)=="X" ) THEN
    temp3 = "zur X-Achse"
  ELSEIF( strings(2)=="y" .OR. strings(2)=="Y" ) THEN
    temp3 = "zur Y-Achse"
  ELSEIF( strings(2)=="z" .OR. strings(2)=="Z" ) THEN
    temp3 = "zur Z-Achse"
  ELSEIF( strings(2)=="xy" .OR. strings(2)=="XY" ) THEN
    temp3 = "zur XY-Neigung"
  ELSEIF( strings(2)=="xz" .OR. strings(2)=="XZ" ) THEN
    temp3 = "zur XZ-Neigung"
  ELSEIF( strings(2)=="yx" .OR. strings(2)=="YX" ) THEN
    temp3 = "zur YX-Neigung"
  ELSEIF( strings(2)=="yz" .OR. strings(2)=="YZ" ) THEN
    temp3 = "zur YZ-Neigung"
  ELSEIF( strings(2)=="zx" .OR. strings(2)=="ZX" ) THEN
    temp3 = "zur ZX-Neigung"
  ELSEIF( strings(2)=="zy" .OR. strings(2)=="ZY" ) THEN
    temp3 = "zur ZY-Neigung"
  ELSEIF( strings(2)=="xyz" .OR. strings(2)=="XYZ" ) THEN
    temp3 = "zu allen Zellvektoren"
  ENDIF
  IF( strings(1)=="add" ) THEN
    msg = ">>> Hinzufügen von "//TRIM(ADJUSTL(temp2))//" A to "//TRIM(ADJUSTL(temp3))//"..."
  ELSEIF( strings(1)=="rm" ) THEN
    msg = ">>> Entfernen von "//TRIM(ADJUSTL(temp2))//" A to "//TRIM(ADJUSTL(temp3))//"..."
  ELSE
    msg = ">>> Setting "//TRIM(ADJUSTL(temp3))//" to "//TRIM(ADJUSTL(temp2))//" A..."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2152)
  msg = "..> Der Zellvektor wurde modifiziert."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2153)
  msg = ">>> Versuch, die Zellvektoren automatisch neu einzustellen..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2154)
  !reals(1) = number of shells detected
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = "..> Erkannte "//TRIM(ADJUSTL(temp))//" Kerne und "//TRIM(ADJUSTL(temp2))//" Schalen."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2600)
  !strings(1) = first option
  !strings(2) = second option
  msg = "<!> INFO: für eine bessere Leistung, verwenden Sie Option '"//TRIM(ADJUSTL(strings(2)))// &
      & "' vor Option '"//TRIM(ADJUSTL(strings(1)))//"'."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!2700-2799: WARNUNG MESSAGES
CASE(2700)
  !strings(1) = option name
  msg = TRIM(ADJUSTL(warnmsg))//" folgende Anweisung konnte nicht interpretiert"// &
      & " werden: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Siehe `atomsk --help options` fuer eine Uebersicht"// &
      & " aller Optionen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2720)
  msg = TRIM(ADJUSTL(warnmsg))//" Achse ist bereits ausgerichtet. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2721)
  msg = TRIM(ADJUSTL(warnmsg))//" Koordinaten sind bereits reduziert. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2722)
  !strings(1) = atom species
  msg = TRIM(ADJUSTL(warnmsg))//" "//TRIM(strings(1))//" Atome haben bereits Schalen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2723)
  !strings(1) = atomic species
  msg = TRIM(ADJUSTL(warnmsg))//" kein "//TRIM(strings(1))//" Atom im System. "// &
      & "Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2724)
  msg = TRIM(ADJUSTL(warnmsg))//" Deformation ist Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2725)
  msg = TRIM(ADJUSTL(warnmsg))//" Burgers vector ist Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2726)
  msg = TRIM(ADJUSTL(warnmsg))//" Superzelle hat eine sehr geringe Ausdehung"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    senkrecht zur Versetzungslinie. Das macht wenig Sinn."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Bitte ueberpruefen!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2727)
  !reals(1) = index of atom with large displacement
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" Grosse Verschiebung fuer Atom #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2728)
  msg = TRIM(ADJUSTL(warnmsg))//" Vergroesserungsfaktoren sind alle 1. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2729)
  msg = TRIM(ADJUSTL(warnmsg))//" Keine Hilfseigenchaft definiert. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2730)
  !string(1) = property
  msg = TRIM(ADJUSTL(warnmsg))//" Hilfseigenschaft "//TRIM(ADJUSTL(strings(1)))// &
      & " nicht gefunden. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2731)
  msg = TRIM(ADJUSTL(warnmsg))//" Hstart = Hend. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2732)
  !strings(1) = name of unknown property
  msg = TRIM(ADJUSTL(warnmsg))//" Unbekannte Eigenschaft: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2733)
  msg = TRIM(ADJUSTL(warnmsg))//" Angegebener Radius ist negativ. "// &
      & "Kein Atom wird verschoben."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2734)
  msg = TRIM(ADJUSTL(warnmsg))//" Rotationswinkel Modulo 2*Pi ist Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2735)
  msg = TRIM(ADJUSTL(warnmsg))//" Scherung ist Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2736)
  msg = TRIM(ADJUSTL(warnmsg))//" Verschiebevektor ist Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2737)
  msg = TRIM(ADJUSTL(warnmsg))//" Atomsorten sind identisch. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2738)
  msg = TRIM(ADJUSTL(warnmsg))//" Einheiten sind identisch. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2739)
  msg = TRIM(ADJUSTL(warnmsg))//" Basisvektoren sind nicht orthonormal."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2740)
  msg = TRIM(ADJUSTL(warnmsg))//" Elastischer Tensor ist nicht symmetrisch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2741)
  msg = TRIM(ADJUSTL(warnmsg))//" Poisson Verhaeltnis ist ausserhalb des"// &
      & " Bereichs [-1 , 0.5]."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2742)
  msg = TRIM(ADJUSTL(warnmsg))//" Ungueltiger Atomindex. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2743)
  msg = TRIM(ADJUSTL(warnmsg))//" Verzerrungsparameter sind Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2744)
  msg = TRIM(ADJUSTL(warnmsg))//" Keine ionischen Schalen im System. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2745)
  !reals(1) = number of atoms that will actually be removed
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" kann so viele Atome nicht auswaehlen."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Es werden nur "//TRIM(ADJUSTL(msg))//" entfernt."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2746)
  msg = TRIM(ADJUSTL(warnmsg))//" Nichts mehr zu tun. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2747)
  msg = TRIM(ADJUSTL(warnmsg))//" Spannungsintensitaet K ist Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2748)
  msg = TRIM(ADJUSTL(warnmsg))//" Superzelle ist sehr klein senkrecht zur "
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Risslinie. Das macht wenig Sinn. Bitte ueberpruefen!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2749)
  !strings(1) = name of property
  !reals(1) = atom index
  WRITE(temp,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" Konnte Wert der Eingenschaft '"// &
      & TRIM(ADJUSTL(strings(1)))//"' nicht dem Atom #"// &
      & TRIM(ADJUSTL(temp))//" zuweisen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2750)
  msg = TRIM(ADJUSTL(warnmsg))//" Eine gesetzte Auswahl enthaehlt keine Atome mehr."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Die Auswahl wurde geloescht. "// &
      & "Alle Atome sind jetzt ausgewaehlt."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2751)
  msg = TRIM(ADJUSTL(warnmsg))//" Zieltemperatur ist Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2752)
  msg = TRIM(ADJUSTL(warnmsg))//" Keine Auswahl definiert. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2753)
  !reals(1) = atom index
  WRITE(temp,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" atom #"//TRIM(ADJUSTL(temp))//"  wurde bereits ausgewählt. Ueberspringe."
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
    msg = TRIM(ADJUSTL(warnmsg))//" Das "//TRIM(ADJUSTL(strings(1)))//" wird aus der Box platziert."
    CALL DISPLAY_MSG(1,msg,logfile)
    msg = "    Das macht wenig Sinn. Bitte ueberpruefen!"
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(2755)
  msg = TRIM(ADJUSTL(warnmsg))//" Faktor ist Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2756)
  !strings(1) = direction
  msg = TRIM(ADJUSTL(warnmsg))//" supercell is quite large along the "//TRIM(ADJUSTL(strings(1)))//" direction!"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Are you sure you know what you are doing?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2757)
  msg = TRIM(ADJUSTL(warnmsg))//" die Indizes sind gleich. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2758)
  msg = TRIM(ADJUSTL(warnmsg))//" keine Operation zum Anwendenno operation to apply. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2759)
  msg = TRIM(ADJUSTL(warnmsg))//" Schleifenradius ist zu klein. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2760)
  msg = TRIM(ADJUSTL(warnmsg))//" Zelle ist bereits orthorhombisch. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2761)
  !strings(1) = "add" or "rm" or "intersect" or "xor"
  msg = TRIM(ADJUSTL(warnmsg))//" Auswahl kann nicht geändert werden, da zuvor keine Auswahl definiert wurde."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2762)
  !strings(1) = elastic tensor stability criterion
  msg = TRIM(ADJUSTL(warnmsg))//" das Kriterium der elastischen Tensorstabilität ist nicht erfüllt: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2763)
  msg = TRIM(ADJUSTL(warnmsg))//" Modulo ist gleich 1 und wählt alle Atome aus."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2764)
  msg = TRIM(ADJUSTL(warnmsg))//" es wurde kein kürzerer periodischer Vektor gefunden, das System bleibt das gleiche."
  CALL DISPLAY_MSG(1,msg,logfile)
  !
CASE(2799)
  !strings(1) = name of obsolete option
  !strings(2) = name of new option
  msg = TRIM(ADJUSTL(warnmsg))//" Option "//TRIM(ADJUSTL(strings(1)))// &
      & " ist ueberholt und wird bald entfernt."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Bitte benutzen Sie stattdessen die Option "// &
      & TRIM(ADJUSTL(strings(2)))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!2800-2899: FEHLER MESSAGES
CASE(2800)
  !string(1) = axis (if we are here it's because it is different from x, y or z)
  msg = TRIM(ADJUSTL(errmsg))//" Unbekannte Achse: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2801)
  msg = TRIM(ADJUSTL(errmsg))//" Die Basis Hend kann nicht auf"// &
             & " Hstart gedreht werden."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Ueberpruefen Sie bitte, ob die Winkel in"// &
             & " Hstart and Hend identisch sind."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2802)
  !strings(1) = property that was not read properly
  msg = "X!X FEHLER beim Lesen von "//TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2803)
  msg = "X!X FEHLER beim Bestimmen der Basisvektoren der Superzellen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2804)
  msg = "X!X FEHLER beim Ausfuehren von Anweisungen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2805)
  !strings(1) = name of unknown option
  msg = TRIM(ADJUSTL(errmsg))//" Unbekannte Anweisung: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2806)
  !strings(1) = name of option
  msg = TRIM(ADJUSTL(errmsg))//" Unzulaessige Anweisungsform fuer "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2807)
  !reals(1) = 1 if roots of Eq.(13-85) cannot be found
  !         = 2 if the A_k(n) cannot be calculated
  !         = 3 if the linear equations defining D(n) cannot be solved
  msg = TRIM(ADJUSTL(errmsg))//""
  IF(NINT(reals(1))==1) THEN
    msg = TRIM(msg)//" P(n) kann nicht bestimmt werden. Abbruch."
  ELSEIF(NINT(reals(1))==2) THEN
    msg = TRIM(msg)//" A_k(n) kann nicht bestimmt werden. Abbruch."
  ELSEIF(NINT(reals(1))==3) THEN
    msg = TRIM(msg)//" D(n) kann nicht bestimmt werden. Abbruch."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2808)
  msg = TRIM(ADJUSTL(errmsg))//" Kann keine gemischte Versetzung mit isotroper"// &
      & " Elastizitaet bilden."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          Sie können eine Kanten- und eine Schraubenversetzung kombinieren,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          oder den elastischen Tensor definieren, um die anisotrope Elastizität zu verwenden.."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2809)
  msg = TRIM(ADJUSTL(errmsg))//" Elastizitaetstensor enthaelt NaN Werte. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2810)
  msg = TRIM(ADJUSTL(errmsg))//" Inkonsistenz in den Datenlisten fuer Schalen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2811)
  msg = TRIM(ADJUSTL(errmsg))//" Versetzungsline und Versetzungsebene muessen"// &
      & " senkrecht aufeinander stehe. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2812)
  msg = TRIM(ADJUSTL(errmsg))//" Kann zu loeschende Atome nicht festlegen. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2813)
  !strings(1) = string that could not be converted to a number
  msg = TRIM(ADJUSTL(errmsg))//" Kann folgende Zeichenkette nicht in eine Zahl"// &
      & " umwandeln: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2814)
  msg = TRIM(ADJUSTL(errmsg))//" Kein atomares System, um diese Anweisung anzuwenden."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2815)
  !strings(1) = name of matrix to invert
  msg = TRIM(ADJUSTL(errmsg))//" Kann die Matrix "//strings(1)//" nicht invertieren."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2816)
  msg = TRIM(ADJUSTL(errmsg))//" Elastizitaetstensor nicht definiert. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2817)
  msg = TRIM(ADJUSTL(errmsg))//" Die Eigenschaft '"//TRIM(ADJUSTL(strings(1)))//"' ist nicht definiert. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2818)
  msg = TRIM(ADJUSTL(errmsg))//" Beim Lesen der STL-Datei ist ein Fehler aufgetreten. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2819)
  msg = TRIM(ADJUSTL(errmsg))//" nicht in der Lage, eine orthogonale Zelle von anfänglichen Zellenvektoren zu finden."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2820)
  !reals(1) = estimated number of atoms
  msg = TRIM(ADJUSTL(errmsg))//" Neue Zelle kann nicht mit Atomen gefüllt werden."
  CALL DISPLAY_MSG(1,msg,logfile)
  WRITE(temp,*) NINT(reals(1))
  msg = "           Schätzung der Anzahl der erforderlichen Atome: "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2821)
  msg = TRIM(ADJUSTL(errmsg))//" Modulo kann nicht umsonst sein (Division durch Null)."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2822)
  msg = TRIM(ADJUSTL(errmsg))//" diese Option akzeptiert keine Operationen mit 'box'."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "           Bitte Entfernungen in Angström angeben."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!
!
!**************************
! 3000-3999: MESSAGES FOR OUTPUT
CASE(3000)
  !reals(2) = number of shells
  WRITE(temp,*) NINT(reals(1))
  IF( NINT(reals(2))>0 ) THEN
    WRITE(temp2,*) NINT(reals(2))
    msg = ">>> Schreibe Ausgabedatei(en) ("//TRIM(ADJUSTL(temp))//" cores + "// &
        & TRIM(ADJUSTL(temp2))//" shells):"
  ELSE
    msg = ">>> Schreibe Ausgabedatei(en) ("//TRIM(ADJUSTL(temp))//" Atomen):"
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3001)
  !strings(1) = name of file
  msg = "..> Diese Datei existiert bereits und wurde nicht erneut konvertiert: "// &
      & TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3002)
  !strings(1) = type of file, e.g. "XSF" or "CFG"
  !strings(2) = name of file
  msg = "..> Datei "//TRIM(strings(1))//" erfolgreich geschrieben: " &
      & //TRIM(strings(2))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3003)
  !reals(:) = list of atomic numbers
  msg = "..> Liste der Atomsorten:"
  DO i=1,SIZE(reals(:))
    CALL ATOMSPECIES(reals(i),species)
    msg = TRIM(msg)//" "//TRIM(species)//", "
  ENDDO
  j=LEN_TRIM(msg)
  msg = TRIM(msg(:j-1))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "    Ueberpruefe ob das konsistent mit der POTCAR Datei ist."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3004)
  !strings(1) = skew parameter
  msg = "..> Zellen Scherung "//TRIM(strings(1))//" verringert."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3005)
  msg = "..> OK, lasse Scherparameter "//TRIM(strings(1))//" unveraendert."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3006)
  !strings(1) = name of ddplot file
  msg = ">>> Erstelle ddplot Datei: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3007)
  msg = ">>> Das Ausgabeformat ist NULL, es wird keine Datei geschrieben."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!3700-3799: WARNUNG MESSAGES
CASE(3700)
  msg = TRIM(ADJUSTL(warnmsg))//" kein Name fuer Ausgabedatei vergeben, bitte eingeben:"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3701)
  !strings(1) = name of output file
  msg = TRIM(ADJUSTL(warnmsg))//" Format der Ausgabedatei konnte nicht erkannt werden: " &
      & //TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "<?> Welches Format der Ausgabedatei soll ich verwenden:"
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
  msg = TRIM(ADJUSTL(warnmsg))//" ddplot Format ist nur im DDPLOT Modus verfuegbar."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    dd Datei wird nicht geschrieben, siehe Dokumentation."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3703)
  msg = TRIM(ADJUSTL(warnmsg))//" Atomsorten sind nicht zusammengefasst."// &
      & " Soll ich das jetzt machen? ("//langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    (Das hat nur Einfluss auf die POSCAR Datei)"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3704)
  msg = TRIM(ADJUSTL(warnmsg))//" Die Zelle spannt keine untere Dreiecksmatrix auf."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    LAMMPS erwartet das aber."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Soll die Zelle jetzt neu ausgerichtet werden? (" &
      & //langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    (Das hat nur Einfluss auf die LAMMPS Datei)"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3705)
  !strings(1) = skew parameter
  msg = TRIM(ADJUSTL(warnmsg))//" Triklinischer Zellscherung "// &
      & TRIM(strings(1))//" ist zu gross."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    LAMMPS wird mit dieser Superzelle eine aehnliche"// &
      & " Fehlermeldung geben."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Moechten sie die Scherung jetzt reduzieren? ("&
      & //langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    (Diese Aktion hat nur Einfluss auf die LAMMPS Ausgabedatei)"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3706)
  !strings(1) = skew parameter
  msg = TRIM(ADJUSTL(warnmsg))//" Kann Zellscherung von "//TRIM(strings(1))// &
      & " nicht reduzieren. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3707)
  !reals(1) = number of atoms
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" Atomanzahl ist sehr gross: "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Reduzieren die sie Anzahl mit der -cut Anweisung"// &
      & " von Atomsk,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    falls ddplot nicht so viele Atome darstellen kann."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3708)
  msg = TRIM(ADJUSTL(warnmsg))//" Nur die ersten 32 Hilfseigenschaften werden"// &
      & " in die CFG Datei geschrieben."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3709)
  msg = TRIM(ADJUSTL(warnmsg))//" Einige Parameter der Superzelle koennen nicht"// &
      & " geschrieben werden."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Die Ausgabedatei wird Fehler enhalten!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3710)
  !strings(1) = file format
  msg = TRIM(ADJUSTL(warnmsg))//" Unbekanntes Format '"//TRIM(strings(1))// &
      & "'. Ueberspringe "
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3711) ! missing occupancy data
  msg = TRIM(ADJUSTL(warnmsg))//" Besetzungsdaten fehlen."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "             Verwende Standard Besetzungen von 1 fuer alle Atome."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3712) ! missing occupancy data for cel file output
  msg = TRIM(ADJUSTL(warnmsg))//" Daten zur thermischen Vibration fehlen."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "             Setze die Biso Werte aller Atome auf Null."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3713) ! missing absorption data
  msg = TRIM(ADJUSTL(warnmsg))//" absorption factors are missing, they will be set to 0.03 for all atoms."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3714)
  msg = TRIM(ADJUSTL(warnmsg))//" einige Atome haben einen ungültigen 'type'."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Sie können die Option '-remove-property type' verwenden, um Atomtypen zu entfernen,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            oder die Option '-properties', um Atomtypen manuell festzulegen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3715)
  msg = TRIM(ADJUSTL(warnmsg))//" Eingabedaten enthalten Besetzungen, die nicht"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            von einigen Ausgabeformaten unterstützt werden."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Einige Atome können sich in der Ausgabedatei überlappen, die nicht physisch ist."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3716)
  msg = TRIM(ADJUSTL(warnmsg))//" Daten enthalten Shells, die nicht"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            von einigen Ausgabeformaten unterstützt werden."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Shells werden in einer Ausgabedatei verloren gehen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3717)
  !reals(1) = total charge
  WRITE(temp,'(f9.3)') reals(1)
  msg = TRIM(ADJUSTL(warnmsg))//" Zelle hat eine von Null verschiedene elektrische Ladung: Q_total = "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3718)
  msg = TRIM(ADJUSTL(warnmsg))//" einige Atome verschiedener Arten haben den gleichen 'type'."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Sie können die Option '-remove-property type' verwenden, um Atomtypen zu entfernen,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            oder die Option '-properties', um Atomtypen manuell festzulegen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3719)
  !reals(1) = atom "type"
  WRITE(temp,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" neuen Atomen wurde der Typ "//TRIM(ADJUSTL(temp))//" zugewiesen."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!3800-3899: FEHLER MESSAGES
CASE(3800)
  msg = TRIM(ADJUSTL(errmsg))//" Keine Atomposition zu schreiben. Breche ab."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3801)
  !strings(1) = name of file
  msg = "X!X FEHLER waehrend des Schreibens der Datei: "// &
      & TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3802)
  msg = "X!X FEHLER waehrend des Schreibens von Dateien."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3803)
  !reals(1) = index of atom that has a NaN coordinate
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" Atom #"//TRIM(ADJUSTL(msg))// &
      & " hat eine nicht-numerische (NaN) Koordinate. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3804)
  !reals(1) = index of atom that has a NaN coordinate
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" Wert der Hilfseigenschaft von Atom #"//&
      & TRIM(ADJUSTL(msg))//" ist nicht-numerisch (NaN). Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!
!
!**************************
! 4000-4999: MESSAGES FOR MODES
CASE(4000)
  !strings(1) = name of file for mode "list"
  msg = ">>> Lese Dateiliste aus: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4001)
  !Just a separator
  msg = "         --===============--"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4002)
  !strings(1) = output format activated, e.g. "xyz" or "cfg"
  msg = ">>> Konvertiere alle folgenden Dateien zu: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4010)
  msg = ">>> Benutze ddplot Modus."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4020)
  !strings(1) = file containing several snapshots
  msg = ">>> Entpacke "//TRIM(strings(1))//" in mehrere Dateien..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4021)
  msg = ">>> Atomsk ist eine freie Open Source software."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "    Für weitere Informationen, geben Sie 'license' ein."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4022)
  msg = ">>> Atomsk command-line Interpreter:"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = '..> Schreibe "help" fuer eine Uebersicht der Kommandos.'
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = ""
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4023)
  msg = "Verfuegbare Kommandos:"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "help                       Zeigt diese Hilfe"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "cd                         Ändern das aktuelle Arbeitsverzeichnis"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = system_ls
  msg(28:) = "Zeigt Dateien im aktuellen Verzeichnis"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "pwd                        Zeigt das aktuelle Arbeitsverzeichnis an"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "print                      Zeigt die aktuelle Boxvektoren und Atompositionen"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "memory                     Zeigt Speicheruebersicht"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "create                     Erstellt ein atomares System"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "read <file>                Datei <file> in den Speicher einlesen"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "write <file>               Schreibt aktuelles System in die Datei <file>"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "box <H11> <H22> <H33>      Definiert die Abmessungen der orthogonalen Box"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "atom <sp> <x> <y> <z>      Fügen dem System ein neues Atom mit bestimmten Arten und Koordinaten hinzu"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "C11 <value>                Einstellt der elastischen Konstante (C11,C22,C33,C12,C13,C23,C44,C55,C66)"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "Cij                        Erstellt einen Steifheitstensor basierend auf vorherigen Werten"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "clear                      Loescht den Speicher (zerstoert atomares System)"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "quit                       Beendet atomsk"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "OPTIONS: Atomsk Anweisungen fuer den Command-Line Interpreter,"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "        gebe 'help options' ein um verfuegbare Anweisungen"// &
      & " aufzulisten."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "        Im interaktiven Modus muessen Anweisungen mit einem"// &
      & " (-) Symbol beginnen."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "MODI: Modi koennen nicht benutzen werden in den Command-Line Interpreter."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4024)
  msg = "<?> In welches Format soll konvertiert werden?"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = '    ('//TRIM(flist(1,1))
  DO i=2,SIZE(flist,1)
    msg = TRIM(msg)//','//TRIM(flist(i,1))
  ENDDO
  msg = TRIM(msg)//')'
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4025)
  msg = ">>> Oops, tut mir leid, ich kann diese Datei nicht"// &
      & " konvertieren..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "    Ich hoffe Du magst mich trotzdem noch! :-)"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4026)
  msg = ">>> Die Datei wurde erfolgreich konvertiert."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4027)
  msg = ">>> Erstelle System:"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4028)
  !strings(1) = description of system that is created
  msg = "..> "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4029)
  msg = "..> System erfolgreich erstellt."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4030)
  !strings(1) = merge direction
  !reals(1) >0 if systems are concatenated, <0 otherwise
  IF( reals(1)>0.d0 ) THEN
    msg = ">>> Stapeln der Systeme entlang "//TRIM(strings(1))//"..."
  ELSE
    msg = ">>> Verbinde die Systeme in der gleichen Box..."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4031)
  !reals(1) = number of files merged
  WRITE(msg,*) NINT(reals(1))
  msg = "..> "//TRIM(ADJUSTL(msg))//" Systeme verbunden."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4032)
  msg = ">>> Berechne elektrische Dipolmomente:"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4033)
  !strings(1) = name of file where total polarization is written
  msg = "..> Totale Polarisation gespeichert in: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4034)
  !strings(1) = first atomic species
  !strings(2) = second atomic species
  !reals(1) = coordination
  IF(INT(reals(1))==4) THEN
    msg = " Tetraeder."
  ELSEIF(INT(reals(1))==6) THEN
    msg = " Oktaeder."
  ELSEIF(INT(reals(1))==8) THEN
    msg = " Hexaeder."
  ELSE
    msg = " Polyeder."
  ENDIF
  msg = ">>> Berechne Momente von "//TRIM(strings(1))// &
      & "-"//TRIM(strings(2))//TRIM(msg)
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  IF(reals(1)<0.d0) THEN
    WRITE(msg,"(f16.2)") DABS(reals(1))
    msg = "..> Beruecksichtige Nachbarn im Radius von "// &
        & TRIM(ADJUSTL(msg))//" A um "//TRIM(strings(1))//" Ionen."
  ELSEIF(reals(1)==0.d0) THEN
    msg = "..> Versuche naechste Nachbarn automatisch zu finden."
  ELSE
    WRITE(msg,*) INT(reals(1))
    msg = "..> Beruecksichtige die "//TRIM(ADJUSTL(msg))// &
        & " ersten Nachbarn der "//TRIM(strings(1))//" Ionen."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4035)
  !strings(1) = file name
  msg = "..> XSF Datei erfolgreich gespeichert: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4036)
  !strings(1) = file name
  msg = "..> Normierung erfolgreich gespeichert: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4037)
  !strings(1) = file name
  msg = "..> Statistik erfolgreich gespeichert: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4038)
  msg = ">>> Berechne Differenzen..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4039)
  !strings(1) = name of file
  msg = "..> Datei erfolgreich geschrieben: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4040)
  msg = "..> Verschiebungen wurden berechnet."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4041)
  !reals(1) = snapshots number
  WRITE(msg,*) NINT(reals(1))
  msg = ">>> Lese Schnappschuss #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4042)
  !reals(1) = number of snapshots read
  WRITE(msg,*) NINT(reals(1))
  IF( NINT(reals(1)) <= 0 ) THEN
    msg = ">>> Kein Schnappschuss konvertiert."
  ELSEIF( NINT(reals(1)) == 1 ) THEN
    msg = ">>> 1 Schnappschuss erfolgreich konvertiert."
  ELSE
    msg = ">>> "//TRIM(ADJUSTL(msg))// &
        & " Schnappschuesse erfolgreich konvertiert."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4043)
  msg = "..> Schnappschuss erfolgreich geladen."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4044)
  msg = ">>> Schreibe Eigenschaften aller Atome:"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4045)
  !reals(1) = number of files converted
  !reals(1) = number of files ignored
  WRITE(msg,*) NINT(reals(1))
  IF( NINT(reals(1))==0 ) THEN
    msg = ">>> Keine wurde Datei konvertiert"
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = ">>> 1 Datei wurde konvertiert"
  ELSE
    msg = ">>> "//TRIM(ADJUSTL(msg))//" Dateien wurden konvertiert"
  ENDIF
  IF( SIZE(reals)>1 .AND. NINT(reals(2))>0 ) THEN
    WRITE(temp,*) NINT(reals(2))
    IF( NINT(reals(2))==1 ) THEN
      msg = TRIM(ADJUSTL(msg))//", 1 Datei wurde ignoriert"
    ELSE
      msg = TRIM(ADJUSTL(msg))//", "//TRIM(ADJUSTL(temp))// &
          & " Dateien wurden ignoriert"
    ENDIF
  ENDIF
  msg = TRIM(ADJUSTL(msg))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4046)
  msg = ">>> Verschiebe Atome periodisch nach Referenz..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4047)
  !reals(1) = number of atoms unwrapped
  IF( NINT(reals(1))==0 ) THEN
    msg = "..> Kein Atom wurde verschoben."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 wurde verschoben."
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" Atome wurden verschoben."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4048)
  msg = ">>> Berechne elektrische Polarisation aller Ionen im System..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4049)
  msg = "..> Elektrische Polarisation wurde berechnet."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4050)
  !strings(1) = file containing the names of files to include
  msg = ">>> Fuege alle Dateien aus "//TRIM(strings(1))// &
      & " zu eine Datei hinzu..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4051)
  !reals(1) = maximum distance
  !reals(2) = width of the skin (A)
  msg = ">>> Berechne radiale Verteilungsfunktion (RDF),"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(temp,'(f16.3)') reals(1)
  WRITE(temp2,'(f16.3)') reals(2)
  msg = "    bis zu "//TRIM(ADJUSTL(temp))//" A, innerhalb einer Schale der Breite "//TRIM(ADJUSTL(temp2))//" A."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4052)
  msg = "..> Berechne RDF..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4053)
  !reals(1) = number of files analyzed
  WRITE(temp,*) NINT(reals(1))
  IF( NINT(reals(1))<=0 ) THEN
    msg = "..> Es wurde keine Datei analysiert."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> RDFs werden aus 1 Datei berechnet."
  ELSE
    msg = "..> RDFs werden über "//TRIM(ADJUSTL(temp))//" Dateien gemittelt."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4054)
  msg = ">>> Konstruiere Polykristall mit der Voronoi Methode."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4055)
  !reals(1) = index of the grain
  WRITE(msg,*) NINT(reals(1))
  msg = ">>> Erstelle Korn #"//TRIM(ADJUSTL(msg))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4056)
  !reals(1) = number of atoms in the grain
  WRITE(msg,*) NINT(reals(1))
  msg = "..> Anzahl der Atome in diesem Korn: "//TRIM(ADJUSTL(msg))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4057)
  !strings(1) = name of input file
  msg = ">>> Lese Parameter fuer Voronoi-Konstruktion von: "// &
      & TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4058)
  !reals(1) = number of grains
  !reals(2) = 0 if 3-D, 1,2,3 if thin along x, y, z
  msg = "..> Datei erfolgreich eingelesen."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,*) NINT(reals(1))
  msg = "..> Anzahl der zu erstellenden Koerner: "// &
      & TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  IF( NINT(reals(2))==0 ) THEN
    msg = "..> Verwende eine 3-D Voronoi Parkettierung."
  ELSE
    IF( NINT(reals(2))==1 ) THEN
      msg = "x"
    ELSEIF( NINT(reals(2))==2 ) THEN
      msg = "y"
    ELSEIF( NINT(reals(2))==3 ) THEN
      msg = "z"
    ENDIF
    msg = "..> Benutze eine 2-D Voronoi Parkettierung"// &
        & ", Rotationsachse: "//TRIM(ADJUSTL(msg))
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4059)
  msg = ">>> Konstruiere eine Kette interpolierter Konfigurationen."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4060)
  !reals(1) = index of the image
  WRITE(msg,*) NINT(reals(1))
  msg = ">>> Konstruiere Konfiguration #"//TRIM(ADJUSTL(msg))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4061)
  msg = ">>> Berechne Nye Tensor."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4062)
  msg = ">>> Berechne atomare G Matrix..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4063)
  msg = ">>> Berechne atomaren Nye Tensor..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4064)
  !strings(1) = name of the file containing file list
  msg = ">>> Mittle Atompositionen fuer eine Dateiliste: "// &
      & TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4065)
  !reals(1) = number of files
  WRITE(temp,*) NINT(reals(1))
  msg = ">>> Atompositionen wurden ueber "//TRIM(ADJUSTL(temp))// &
      & " Konfigurationen gemittelt."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4066)
  msg = ">>> Laufen in der Modus density."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4067)
  !strings(1) = property
  !reals(1) = dimension (1, 2 or 3)
  WRITE(temp,*) NINT(reals(1))
  msg = ">>> Berechnen der "//TRIM(ADJUSTL(temp))//"-D Dichte der "//TRIM(ADJUSTL(strings(1)))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4068)
  WRITE(msg,*) NINT(reals(1))
  msg = "..> Dichte wurde erfolgreich berechnet."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4069)
  !strings(1) = name of file
  msg = ">>> Berechnen des lokalen Symmetrieparameters für die Datei: "//TRIM(ADJUSTL(strings(1)))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4070)
  !strings(1) = name of file
  !reals(2) = 1 (reference is a unit cell) or 2 (no reference)
  IF( NINT(reals(1))==1 ) THEN
    msg = ">>> Elementarzelle als Referenz : "//TRIM(ADJUSTL(strings(1)))//"..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSE
    msg = ">>> Kein Bezugssystem vorhanden. Erzeugen atomarer Referenzumgebungen durch Mittelung"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    msg = "    von Standorten aus dem bereitgestellten System: "//TRIM(ADJUSTL(strings(1)))//"..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4071)
  !reals(1) = number of different atomic environments found
  IF( NINT(reals(1))>1 ) THEN
    msg = "..> Fertig, "//TRIM(ADJUSTL(temp))//" verschiedene atomare Umgebungen gefunden."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> Fertig, 1 atomare Umgebung gefunden."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4072)
  !reals(1) = number of files analyzed
  WRITE(temp,*) NINT(reals(1))
  IF( NINT(reals(1))<=0 ) THEN
    msg = "..> Es wurde keine Datei analysiert."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 Datei wurde analysiert."
  ELSE
    msg = "..> "//TRIM(ADJUSTL(temp))//" Dateien wurden analysiert."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4073)
  !strings(1) = name of 1st file
  !strings(2) = name of 2nd file
  msg = ">>> Hilfseigenschaften von "//TRIM(ADJUSTL(strings(1))) &
      & //" nach "//TRIM(ADJUSTL(strings(2)))//" kopieren."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4074)
  !reals(1) = number of aux.prop.copied
  !reals(2) = number of new aux.prop.
  !reals(3) = number of overwritten aux.prop.
  !reals(4) = new total number of aux.prop.
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  WRITE(temp3,*) NINT(reals(3))
  msg = ">>> "//TRIM(ADJUSTL(temp))//" Eigenschaften wurden kopiert: "//TRIM(ADJUSTL(temp2))//" neu, " &
      & //TRIM(ADJUSTL(temp3))//" überschrieben."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(temp,*) NINT(reals(4))
  msg = "..> Das System enthält jetzt "//TRIM(ADJUSTL(temp))//" Hilfseigenschaften."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4075)
  !strings(1) = string containing user values
  msg = "..> Benutzerdefinierte Werte erzwingen: "
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  DO i=1,SIZE(strings)
    IF( LEN_TRIM(strings(i))>0 ) THEN
      msg = "    "//TRIM(ADJUSTL(strings(i)))
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ENDIF
  ENDDO
CASE(4200)
  WRITE(*,*) " (Gib 'q' ein um abzubrechen)"
  WRITE(*,'(a39)',ADVANCE='NO') " Gittertyp (sc,bcc,fcc,dia,rs,per): "
CASE(4201)
  !strings(1) = name of lattice parameter (a, b, c, a0...)
  WRITE(*,'(a28)',ADVANCE='NO') " Gitterparameter "// &
      & TRIM(ADJUSTL(strings(1)))//" (Å): "
CASE(4202)
  WRITE(*,'(a15)',ADVANCE='NO') " Atom Spezies: "
CASE(4203)
  WRITE(*,'(a22)',ADVANCE='NO') " Chirale Indizes (n,m): "
CASE(4204)
  WRITE(*,'(a21)',ADVANCE='NO') " Gitterorientierung: "
CASE(4300)
  !reals(1) = number of max tries
  msg = "   ***   ATOMSK QUIZZ   ***"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(temp,*) NINT(reals(1))
  msg = " Beantworte die folgende Frage! Du hast "// &
      & TRIM(ADJUSTL(temp))//" Versuche."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4301)
  !strings(1)
  !reals(1) = 
  !reals(2) = 
  IF( NINT(reals(2))>0 ) THEN
    msg = " FALSCH!"
  ELSE
    msg = ""
  ENDIF
  IF( NINT(reals(1))==1 ) THEN
    msg = TRIM(msg)//" Was ist die Ordnungszahl von "// &
        & TRIM(ADJUSTL(strings(1)))//"?"
  ELSEIF( NINT(reals(1))==2 ) THEN
    msg = TRIM(msg)//" Welches Element hat die Ordnungszahl "// &
        & TRIM(ADJUSTL(strings(1)))//"?"
  ELSE
    msg = TRIM(msg)//" Welches Element kommt nach "// &
        & TRIM(ADJUSTL(strings(1)))//" im Periodensystem?"
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4302)
  !reals(1) = try/maxtry
  IF( reals(1)>=1.d0 ) THEN
    msg = " GAME OVER!"
  ELSE
    msg = " RICHTIG!"
  ENDIF
  msg = TRIM(msg)//" Die richtige Antwort wäre: "//TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!4700-4799: WARNUNG messages
CASE(4700)
  !strings(1) = name of file that does not exist
  msg = TRIM(ADJUSTL(warnmsg))//" Datei nicht gefunden: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Ueberspringe..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4701)
  msg = TRIM(ADJUSTL(warnmsg))//" Konnte Basisvektoren der Superzelle in der"// &
      & " Eingabedatei nicht finden."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Versuche Basisvektoren zu erraten."// &
      & " Vorsicht, das kann sehr ungenau sein."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Es wird empfohlen die Basisvektoren der Superzelle"// &
      & " mit der Anweisungn '-prop' oder '-cell' zu definieren."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4702)
  msg = TRIM(ADJUSTL(warnmsg))//" Zelle ist geladen."// &
      & " Die totale Polarisation kann falsch sein!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4703)
  msg = TRIM(ADJUSTL(warnmsg))//" Spezies hat keine Ladung, ueberspringe..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4705)
  !reals(1) = index of atom that has too many neighbors
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" Anzahl der Nachbar ist groesser als 100"// &
      & " fuer Atom #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4706)
  !strings(1) = atomic species of absent neighbors
  !strings(2) = atomic species of central atom
  !reals(1) = index of atom that has no neighbour
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" Konnte keinen "//TRIM(strings(1))// &
      &       " Nachbarn fuer #"//TRIM(ADJUSTL(msg))//" ("//TRIM(strings(2))//") finden."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4707)
  !reals(1) = index of atom that has a shell with zero charge
  WRITE(msg,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" Die Schale von Ion #"//TRIM(ADJUSTL(msg)) &
      & //" hat keine Ladung."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4708)
  msg = TRIM(ADJUSTL(warnmsg))//" Dieses Korn enthaelt kein Atom."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4709)
  !reals(1) = index of atom
  !reals(2) = number of neighbors
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(warnmsg))//" Unzureichende Anzahl von Nachbarn ("// &
      & TRIM(ADJUSTL(temp2))//") fuer Atom #"//TRIM(ADJUSTL(temp))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4710)
  !strings(1) = name of the matrix
  msg = TRIM(ADJUSTL(warnmsg))//" Konnte die Matrix "//TRIM(ADJUSTL(strings(1)))// &
      & " nicht berechnen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4711)
  !strings(1) = name of the file
  msg = TRIM(ADJUSTL(warnmsg))//" Abweichende Anzahl von Atomen in der Datei: "// &
      & TRIM(ADJUSTL(strings(1)))//". Datei wird nicht verwendet."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4712)
  msg = TRIM(ADJUSTL(warnmsg))//" Einige Atome liegen ausserhalb der Zelle."// &
      & " Das kann zu falschen Ergebnissen fuehren."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Atome jetzt periodisch in die Zelle schieben? ("//&
      & langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4713)
  !strings(1) = action to follow if user says "yes"
  IF(LEN_TRIM(strings(1))<=0) THEN
    strings(1) = "proceed"
  ENDIF
  msg = TRIM(ADJUSTL(warnmsg))//" Offenbar wurde das System noch nicht in einer"// &
      & " Datei gespeichert."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Sind Sie sicher, dass Sie "// &
      & TRIM(ADJUSTL(strings(1)))//" wollen? ("// &
      & langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4714)
  !strings(1) = direction along which cell is small
  !reals(1) = new cell size along that direction
  msg = TRIM(ADJUSTL(warnmsg))//" Die letzte Zelle hat eine kleine Dimension entlang " &
      & //TRIM(ADJUSTL(strings(1)))//", einstellung auf  "//TRIM(ADJUSTL(temp))//" Å."
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
  msg = TRIM(ADJUSTL(warnmsg))//" Ein 2-D-Polykristall soll normal zur "//TRIM(temp)//"-Achse konstruiert werden,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            aber Sie haben Rotationswinkel um die anderen kartesischen Achsen angegeben."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Bitte ueberpruefen!"
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
  msg = TRIM(ADJUSTL(warnmsg))//" Knot #"//TRIM(ADJUSTL(temp4))//&
      & " war außerhalb der Grenzen und wurde zurück in die Schachtel gelegt, neue Position: (" &
      & //TRIM(ADJUSTL(temp))//","//TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))//")."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4717)
  !strings(1) = name of matrix
  msg = TRIM(ADJUSTL(warnmsg))//" "//TRIM(ADJUSTL(strings(1)))//" ist keine Identitätsmatrix."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4718)
  !reals(1) = index of first node
  !reals(2) = index of 2nd node
  !reals(3) = distance between the two nodes
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  WRITE(temp3,'(f18.1)') reals(3)
  msg = TRIM(ADJUSTL(warnmsg))//" Knoten #"//TRIM(ADJUSTL(temp))//" und #"//TRIM(ADJUSTL(temp2))//&
      & " liegen sehr nah beieinander  (d = "//TRIM(ADJUSTL(temp3))//" A)."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Bitte ueberpruefen!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4719)
  !reals(1) = volume of final cell (A^3)
  WRITE(temp,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" letzte Zelle ist sehr klein ("//TRIM(ADJUSTL(temp))//" A^3)."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Bitte ueberpruefen!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4720)
  !strings(1) = chemical symbol
  !reals(1) = file number
  WRITE(temp,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(warnmsg))//" Kein "//TRIM(ADJUSTL(strings(1)))// &
      & " Atom in System N. "//TRIM(ADJUSTL(temp))//" gefunden."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4721)
  !reals(1) = N. atoms in 1st system
  !reals(2) = N. atoms in 2nd system
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  WRITE(temp3,*) MIN(NINT(reals(1)),NINT(reals(2)))
  msg = TRIM(ADJUSTL(warnmsg))//" Anzahl der Atome unterscheiden ("//TRIM(ADJUSTL(temp))// &
      & " vs "//TRIM(ADJUSTL(temp2))//")."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Eigenschaften werden nur für die ersten 4 "//TRIM(ADJUSTL(temp3))//" Atome kopiert."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!4800-4899: FEHLER MESSAGES
CASE(4800)
  !strings(1) = mode
  msg = TRIM(ADJUSTL(errmsg))//" Unzulaessiges Anweisung im Modus: "// &
      & TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4801)
  !strings(1) = file format (e.g. "xyz", "cfg"...)
  msg = TRIM(ADJUSTL(errmsg))//" Dieses Dateiformat wird noch nicht im 1-in-all"// &
      & " Modus unterstuetzt: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4802)
  msg = TRIM(ADJUSTL(errmsg))//" Beide chiralen Indizes sollten nicht Null sein."//&
      & " Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4803)
  !reals(1) = theoretical number of atoms in nanotube
  !reals(2) = number found by Atomsk
  WRITE(msg,*) TRIM(ADJUSTL(errmsg))//" Unerwartete Anzahl von Atomen in der"// &
             & " Nanoroehre."
  CALL DISPLAY_MSG(1,msg,logfile)
  WRITE(msg,*) NINT(reals(1))
  WRITE(temp,*) NINT(reals(2))
  msg = "          Theorie: "//TRIM(ADJUSTL(msg))//"; Vorgefunden: "// &
      & TRIM(ADJUSTL(temp))
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
      temp = TRIM(ADJUSTL(temp))//", "//TRIM(ADJUSTL(temp2))// &
           & " oder "//TRIM(ADJUSTL(temp3))
    ELSE
      temp = TRIM(ADJUSTL(temp))//" or "//TRIM(ADJUSTL(temp2))
    ENDIF
  ENDIF
  msg = TRIM(ADJUSTL(errmsg))//" Diese Struktur benoetigt "//TRIM(ADJUSTL(temp))// &
      & " Elemente."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4805)
  !strings(1) = structure type
  msg = TRIM(ADJUSTL(errmsg))//" Unbekannte Struktur: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4806)
  msg = TRIM(ADJUSTL(errmsg))//" Keine Datei zum Zusammenfuegen. Abbruch"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4807)
  msg = TRIM(ADJUSTL(errmsg))//" Alle Ladungen sind Null. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Ladungen muessen mit der Anweisung -properties"// &
      & " gesetzt werden."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4808)
  !strings(1) = atomic species that has zero charge
  msg = TRIM(ADJUSTL(errmsg))//" "//TRIM(strings(1))//" Ionen haben keine Ladung."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4809)
  !strings(1) = atomic species
  msg = TRIM(ADJUSTL(errmsg))//" Kein "//TRIM(strings(1))//" Ion im System."
  msg = TRIM(ADJUSTL(msg))
CASE(4810)
  !reals(1) = number of particles in first system
  !reals(2) = number of particles in second system
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(errmsg))//" Unterschiedliche Partikelzahl in beiden"// &
      & " Systemen: "//TRIM(ADJUSTL(temp))//"/"//TRIM(ADJUSTL(temp2))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "           Programm wird BEENDET."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4811)
  msg = TRIM(ADJUSTL(errmsg))//" Die beiden Systeme haben unterschiedliche"// &
      & " Basisvektoren."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4812)
  msg = TRIM(ADJUSTL(errmsg))//" Die Datei scheint kein animated XSF"// &
      & " zu enthalten. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4813)
  !strings(1) = name of unknown mode
  msg = TRIM(ADJUSTL(errmsg))//" Unbekannter Modus: "//TRIM(strings(1))
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
    msg = "<?> Wollten Sie den Modus '"//TRIM(ADJUSTL(temp))//"' verwenden?"
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(4814)
  msg = TRIM(ADJUSTL(errmsg))//" es kann immer nur ein Modus aktiv sein. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4815)
  msg = TRIM(ADJUSTL(errmsg))//" Datei scheint nicht im DL_POLY HISTORY Format"// &
      & " zu sein. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4816)
  !reals(1) = index for core charges
  !reals(2) = index for shell charges
  IF( reals(1)<reals(2) ) THEN
    msg = TRIM(ADJUSTL(errmsg))//" Kernladungen der Ionen nicht definiert. Abbruch."
  ELSE
    msg = TRIM(ADJUSTL(errmsg))//" Ladung der ionischen Schalen sind nicht"// &
        & " definiert. Abbruch."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4817)
  msg = TRIM(ADJUSTL(errmsg))//" Es gibt keine ionische Schale im System. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4818)
  !strings(1) = file format (e.g. "xyz", "cfg"...)
  msg = TRIM(ADJUSTL(errmsg))//" Die Dateiliste in "//TRIM(strings(1))// &
      & " scheint leer zu sein."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4819)
  !strings(1) = suggestion for a vector
  !strings(2) = direction of suggested vector (X,Y or Z)
  msg = TRIM(ADJUSTL(errmsg))//" Basisvektoren sind nicht orthogonal."
  CALL DISPLAY_MSG(1,msg,logfile)
  IF( LEN_TRIM(strings(1))>0 ) THEN
    msg = "    Suggested vector along "//TRIM(ADJUSTL(strings(2)))//": "//TRIM(ADJUSTL(strings(1)))
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(4820)
  msg = TRIM(ADJUSTL(errmsg))//" Dimension der Superzelle wurden nicht"// &
      & " definiert. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4821)
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(errmsg))//" Atomanzahl ("//TRIM(ADJUSTL(temp))// &
      & ") zu gross fuer bereitgestellten Speicher ("//TRIM(ADJUSTL(temp2))//")."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4822)
  msg = TRIM(ADJUSTL(errmsg))//" Keine Datei zu bearbeiten. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4823)
  !strings(1) = keyword 1
  !strings(2) = keyword 2
  msg = TRIM(ADJUSTL(errmsg))//" Die Schluesselwoerter '"// &
      & TRIM(ADJUSTL(strings(1)))//"' and '"// &
      & TRIM(ADJUSTL(strings(2)))//"' schliessen sich gegenseitig"// &
      & " aus. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4824)
  !strings(1) = name of unknown command
  msg = TRIM(ADJUSTL(errmsg))//" Unbekanntes Kommando: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4825)
  msg = TRIM(ADJUSTL(errmsg))//" Atomsk kann nicht in sich laufen!"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          Sie haben 'Atomsk' ohne Argumente eingegeben, oder Sie haben aus einem Menü "
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "           auf die ausführbare Datei geklickt, daher läuft Atomsk derzeit im interaktiven Modus,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          wo nur eine begrenzte Untermenge von Befehlen verfügbar ist. Bitte beachten Sie die Dokumentation."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          Um Atomsk mit Kommandozeilenargumenten zu verwenden, musst du zuerst"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          eine Kommando-Shell ausführen und dann deinen Befehl eingeben."
  CALL DISPLAY_MSG(1,msg,logfile)
#if defined(WINDOWS)
  msg = "          Gehen Sie in einem Microsoft Windows-System folgendermaßen vor:"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          1. Öffnen Sie das Windows-Menü, gehen Sie zu Zubehör und führen Sie Windows Powershell aus."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          2. Führen Sie atomsk.exe mit den Argumenten aus, zum Beispiel:"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = '             & "C:\Program Files\Atomsk\atomsk.exe" initial.xsf final.cfg'
  CALL DISPLAY_MSG(1,msg,logfile)
#endif
CASE(4826)
  msg = TRIM(ADJUSTL(errmsg))//" Dieser Modus ist im interaktiven Modus nicht verfügbar. Bitte beachten Sie die Dokumentation."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4827)
  msg = TRIM(ADJUSTL(errmsg))//" Unmöglich, Gitter mit vorgegebener Orientierung zu schaffen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4828)
  msg = TRIM(ADJUSTL(errmsg))//" Mindestens zwei Zellendimensionen sind zu klein. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4829)
  msg = TRIM(ADJUSTL(errmsg))//" Boxvektoren sind nicht linear unabhängig. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4830)
  msg = TRIM(ADJUSTL(errmsg))//" Referenzumgebung kann nicht erstellt werden."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4831)
  msg = TRIM(ADJUSTL(errmsg))//" In der Parameterdatei ist kein Knoten definiert: '"//TRIM(ADJUSTL(strings(1)))//"'."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4832)
  !reals(1) = index of first node
  !reals(2) = index of 2nd node
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(errmsg))//" Knoten #"//TRIM(ADJUSTL(temp))//" und #"//TRIM(ADJUSTL(temp2))//&
      & " sind an der gleichen Position. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4833)
  !reals(1) = value of IBRION
  WRITE(temp,*) NINT(reals(1))
  msg = TRIM(ADJUSTL(errmsg))//" anscheinend wurde VASP mit IBRION = "//TRIM(ADJUSTL(temp))//"  ausgeführt,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          daher enthält diese OUTCAR-Datei keine atomare Konfiguration."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4834)
  !strings(1) = name of file
  msg = TRIM(ADJUSTL(errmsg))//" die Datei "//TRIM(ADJUSTL(strings(1)))// &
      & " enthält keine Hilfseigenschaft, Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4900)
  msg = TRIM(ADJUSTL(errmsg))//" Es kann immer nur ein Modus verwendet werden."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!
!
!Special case for mode PREFERENCES: many questions can be added there
CASE(5000)
  !strings(1) = user's configuration file
  msg = ">>> Anhand der folgenden Fragen werden die atomsk"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "..> Voreinstelungen gesetzt und in "// &
      & TRIM(ADJUSTL(strings(1)))//" gespeichert."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5001)
  msg = "<?> Soll atomsk Dateien standardmaessig immer ueberschreiben?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5002)
  msg = "<?> Soll atomsk bestehende Dateien standardmaessig immer"// &
      & " ignorieren?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5003)
  msg = "<?> Standard Ausgabeformat ('nein' fuer keines):"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5004)
  msg = "<?> Standardsprache ('nein' fuer keine):"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5005)
  msg = "<?> Standard Ausgabestufen fuer Nachrichten:"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    0=stumm; keine Ausgabe von Nachrichten."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    1=auf den Bildschirm (standard)."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    2=in die Logdatei."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    3=auf den Bildschirm und in die Logdatei."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    4=debug; zusaetzliche Nachrichten in die Logdatei."
  CALL DISPLAY_MSG(1,msg,logfile)
!
CASE(5700)
  !strings(1) = user's configuration file
  msg = TRIM(ADJUSTL(warnmsg))//" es ist notwendig die bestehende Datei "// &
      & "zu ueberschreiben: "//TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Sind Sie sicher, dass Sie fortfahren wollen? (ja/nein)"
  CALL DISPLAY_MSG(1,msg,logfile)
!
!
!
CASE DEFAULT
  CALL ATOMSK_MSG_EN(imsg,strings,reals)
END SELECT
!
END SUBROUTINE ATOMSK_MSG_DE
!
!
!
!********************************************************
! DATE_MSG
! Displays a nice message on certain dates.
!********************************************************
SUBROUTINE DATE_MSG_DE()
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
  msg = "*** Arbeiten außerhalb der Buerozeiten? Sie sollten mal schlafen. :-)"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!If it's lunch time
ELSEIF(values(5)==12 .AND. values(6)>=30) THEN
  msg = "*** Ich habe Hunger! :-p"
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
  msg = "*** GLÜCKLICHES NEUES JAHR "//TRIM(msg)//"!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!March 14 (3/14): Pi day
ELSEIF(values(2)==3 .AND. values(3)==14) THEN
  WRITE(msg,'(a10,f19.16,a2)') "    ( π = ", pi, " )"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!March 31
ELSEIF(values(2)==3 .AND. values(3)==31) THEN
  msg = "  WORLD BACKUP DAY - Denken Sie daran, Ihre Daten zu sichern!"
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
  msg = "*** Sie arbeiten? Tag der Arbeit ist ein FEIERTAG! :-)"
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
  msg = " * * *   FROHE WEIHNACHTEN!   * * *"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
ENDIF
!
END SUBROUTINE DATE_MSG_DE
!
!
!
END MODULE messages_de
