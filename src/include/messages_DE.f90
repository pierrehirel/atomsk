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
!* Last modification: P. Hirel - 09 Feb. 2016                                     *
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
IF(helpsection=="modes" .OR. helpsection=="ai1") THEN
  WRITE(*,*) "..> Alle-in-Einen Modus:"
  WRITE(*,*) "          atomsk -AI1 <listfile> <outputfile> [options] [<formats>]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="1ia") THEN
  WRITE(*,*) "..> Einer-in-Alle Modus:"
  WRITE(*,*) "          atomsk -1IA <inputfile> [<outputfile>] [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="create") THEN
  WRITE(*,*) "..> Erstellen Modus:"
  WRITE(*,*) "          atomsk -C <structure> <a0> <species> <outputfile> [<formats>] [options]"
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
!
IF(helpsection=="options") THEN
  WRITE(*,*) ">>> OPTIONEN (Abstaende=Angstroems, Winkel=Grad):"
  WRITE(*,*) ""
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
  WRITE(*,*) "          -disloc <pos1> <pos2> <screw|edge|edge2> <x|y|z> <x|y|z> <b> <ν>"
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
  WRITE(*,*) "..> Roll the system around an axis:"
  WRITE(*,*) "          -roll <x|y|z> <angle> <x|y|z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-rotate" .OR. helpsection=="-rot") THEN
  WRITE(*,*) "..> Rotiere das System um eine Achse:"
  WRITE(*,*) "          -rot <x|y|z> <angle>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-select") THEN
  WRITE(*,*) "..> Waehle Atome nach verschiedenen Kriterien aus:"
  WRITE(*,*) "          -select all"
  WRITE(*,*) "          -select invert"
  WRITE(*,*) "          -select <species>"
  WRITE(*,*) "          -select <index>"
  WRITE(*,*) "          -select <above|below> <d> <normal>"
  WRITE(*,*) "          -select <in|out> <box|sphere|cylinder> <x1> <y1> <z1> <x2> [<y2> <z2>]]"
  WRITE(*,*) "          -select prop <prop> <value>"
  WRITE(*,*) "          -select random <N> <species>"
  WRITE(*,*) "          -select <NNN> <species> neighbors <index>"
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
  WRITE(*,*) "          -sort <s|x|y|z> <up|down|pack>"
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
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-torsion") THEN
  WRITE(*,*) "..> Apply torsion around an axis:"
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
WRITE(*,*) "    atomsk konvertiert von jedem 'ja' in der INPUT Spalte"
WRITE(*,*) "    zu jedem 'ja' in der OUTPUT Spalte:"
WRITE(*,*) "                            |  INPUT | OUTPUT"
WRITE(*,*) "    ------------------------+--------+--------"
WRITE(*,*) "    atsk (atomsk format)    |  ja    |  ja "
WRITE(*,*) "    bop (Bond-Order format) |  ja    |  ja "
WRITE(*,*) "    cel (Dr. Probe/EMS)     |  ja    |  ja "
WRITE(*,*) "    coo (COORAT/MBPP)       |  ja    |  ja "
WRITE(*,*) "    cfg (Atomeye)           |  ja    |  ja "
WRITE(*,*) "    cif (Cryst.Info.File)   |  ja    |  ja "
WRITE(*,*) "    dd  (ddplot)            |  nein  |  ja (1)"
WRITE(*,*) "    dlp (DL_POLY CONFIG)    |  ja    |  ja "
WRITE(*,*) "    gin (GULP input)        |  ja    |  ja "
WRITE(*,*) "    imd (IMD input)         |  ja    |  ja "
WRITE(*,*) "    jems (JEMS input)       |  ja    |  ja "
WRITE(*,*) "    lmc (LAMMPS output)     |  ja    |  nein"
WRITE(*,*) "    lmp (LAMMPS data)       |  ja    |  ja "
WRITE(*,*) "    mol (MOLDY format)      |  ja    |  ja "
WRITE(*,*) "    pdb (Protein Data Bank) |  ja    |  ja "
WRITE(*,*) "    pos (POSCAR/VASP)       |  ja    |  ja "
WRITE(*,*) "    pw (Quantum Espresso)   |  ja    |  ja "
WRITE(*,*) "    pwout (QE output file)  |  ja (2)|  nein"
WRITE(*,*) "    xmd (XMD file)          |  ja    |  ja "
WRITE(*,*) "    xsf (XCrySDen)          |  ja    |  ja "
WRITE(*,*) "    xv (SIESTA format)      |  ja    |  ja "
WRITE(*,*) "    xyz/exyz/sxyz           |  ja    |  ja "
WRITE(*,*) "        (1) nur im ddplot Modus."
WRITE(*,*) "        (2) nur im ein-in-alle Modus."
ENDIF
!
WRITE(*,*) ""
WRITE(*,*) ">>> Schau im mitgelieferten /doc Ordner nach oder "
WRITE(*,*) "    besuche: http://atomsk.univ-lille1.fr/"
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
SUBROUTINE ATOMSK_CREATE_DATE_DE(VALUES,formula,username,msg)
!
IMPLICIT NONE
CHARACTER(LEN=128),INTENT(IN):: username
INTEGER,DIMENSION(8),INTENT(IN):: VALUES
CHARACTER(LEN=128),INTENT(IN):: formula
CHARACTER(LEN=128),INTENT(OUT):: msg
!
WRITE(msg,'(i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') &
  & VALUES(1), "-", VALUES(2),"-", VALUES(3)," ", VALUES(5), ":", VALUES(6), ":", VALUES(7)
!
msg = TRIM(ADJUSTL(formula))//' - Datei mit Atomsk von '//TRIM(ADJUSTL(username))//' am '//TRIM(ADJUSTL(msg))//" generiert."
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
CHARACTER(LEN=128):: msg  !The message to be displayed
CHARACTER(LEN=128):: temp, temp2, temp3
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
  msg = "\o/ Programm erfolgreich beendet!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(temp,*) nwarn
  temp = ADJUSTL(temp)
  WRITE(temp2,*) nerr
  temp2 = ADJUSTL(temp2)
  WRITE(msg,*) "   Warnungen: "//TRIM(temp)//" ; Fehler: "//TRIM(temp2)
  CALL DISPLAY_MSG(verbosity,msg,logfile)
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
    IF( NINT(reals(1)) == NINT(reals(2)) ) THEN
      WRITE(*,*) ""
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
  msg = "/!\ WARNUNG: Nicht erkannte Befehlszeilenargument: "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(750)
  msg = ""
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "/!\ WARNUNG: SIE SOLLTEN DIESES PROGRAM NICHT ALS ROOT EVALUIEREN!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  temp=""
  DO WHILE(temp.NE."ok")
    msg = "    Zum fortfahren drücken Sie 'ok', oder Ctrl+C zum abschließen."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    READ(*,*) temp
  ENDDO
  msg = ""
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
! 800- 899: FEHLER MESSAGES
CASE(800)
  msg = "X!X FEHLER: Unbekanntes Dateiformat."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(801)
  !strings(1) = atom type
  msg = "X!X FEHLER: Unerkannte Atomsorte: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(802)
  !reals(1) = index of atom that caused the error
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = "X!X FEHLER beim Versuch, Atom #"//TRIM(ADJUSTL(msg))//" zu lesen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(803)
  !strings(1) = unit
  msg = "X!X FEHLER: unbekannte Einheit: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(804)
  msg = "X!X FEHLER: Anzahl der Atome ist Null, Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(806)
  msg = "X!X FEHLER: die Systeme weisen nicht die gleiche Anzahl von Atomen auf!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(807)
  !reals(1) = index of line in file that caused the error on read
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = "X!X FEHLER beim Versuch, Zeile #"//TRIM(ADJUSTL(msg))//" zu lesen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(808)
  !strings(1) = string where conversion failed
  msg = "X!X FEHLER beim Konvertieren von '"//TRIM(strings(1))// &
      & "' in einen numerischen Wert."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(809)
  !strings(1) = space group H-M symbol which couldn't be identified
  msg = "X!X FEHLER: unbekannte Raumgruppe: '"//TRIM(strings(1))//"'."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(810)
  !reals(1) = invalid space group number
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = "X!X FEHLER: ungueltige Raumgruppe: "//TRIM(ADJUSTL(msg))
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
  WRITE(temp,*) NINT(reals(1))
  msg = "..> Eingabedatei wurde erfolgreich gelesen ("//TRIM(ADJUSTL(temp))//" Atomen)."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(1002)
  msg = "..> INP Datei erkannt, entnehme Daten der Superzelle..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(1003)
  msg = "..> POTCAR Datei erkannt, entnehme Atomsorten..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!1700-1799: WARNUNG MESSAGES
CASE(1700)
  !strings(1) = auxiliary property that cannot be loaded
  msg = "/!\ WARNUNG: kann Hilfseigenschaften nicht lesen: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1701)
  !strings(1) = input file format
  msg = "/!\ WARNUNG: das wahrscheinlichste Dateiformat ist "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Geringe Sicherheit bei der Bestimmung des Formats."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1702)
  !strings(1) = name of parameter
  !strings(2) = name of custom config file
  msg = "/!\ WARNUNG: Unbekannter Parameter '"//TRIM(strings(1))//&
      & "' in "//TRIM(strings(2))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1703)
  !strings(1) = name of personal config file
  msg = "/!\ WARNUNG: Fehler beim Lesen der Konfigurationsdatei "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "             Einige Einstellungen sind moeglicherweise nicht gesetzt worden."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1704)
  !strings(1) = name of personal config file
  msg = "/!\ WARNUNG: Symmetrie-Opertionen wurden nicht beruecksichtigt."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1705)
  msg = "/!\ WARNUNG: Sowohl celldm(:) als auch die konventionelle Notation liegen vor."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "             celldm(:) wird verwendet, konventionelle Notation wird ignoriert."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1706)
  msg = "/!\ WARNUNG: Zell-Dimensionen in Bohrs, aber atomare Positionen in Angstroems."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "             Atomare Positionen werden in Bohr konvertiert."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1707) ! invalid symmetry operation string input.
  !strings(1) = failed symmetry operation string
  msg = "/!\ WARNUNG: unzulaessige Zeichenkette fuer "// &
      & "Symmetrieoperation '"//TRIM(strings(1))//"'. Ueberspringe..."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!1800-1899: FEHLER MESSAGES
CASE(1800)
  msg = "X!X FEHLER: Format der Eingabedatei konnte nicht bestimmt werden."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Moeglicherweise wird dieses Format von atomsk noch nicht unterstuetzt."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Die Datei wird uebersprungen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1801)
  !strings(1) = file name
  msg = "X!X FEHLER beim Lesen aus der Datei: " &
      & //TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1802)
  !strings(1) = bad array
  msg = "X!X FEHLER: inkonsistente Feldgroesse in "//TRIM(strings(1))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1803)
  msg = "X!X FEHLER: Datenmenge der Hilfseigenschaften "// &
             & "stimmt nicht mit der Anzahl der Atome ueberein."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1804)
  msg = "X!X FEHLER: Unbekanntes Format."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1805)
  !reals(1) = number of particles read
  !reals(2) = number of particles declared
  msg = "X!X FEHLER: Anzahl der eingelesenen Atome stimmt nicht mit"
  CALL DISPLAY_MSG(1,msg,logfile)
  WRITE(temp,*) NINT(reals(1))
  WRITE(msg,*) NINT(reals(2))
  msg = "            der festgelegten Anzahl ueberein: "// &
      & TRIM(ADJUSTL(temp))//"/"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1806)
  !reals(1) = index of atom
  WRITE(msg,*) NINT(reals(1))
  msg = "X!X FEHLER: Nenner ist Null fuer die Koordinate von Atom #"// &
      & TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1807)
  msg = "X!X FEHLER: Datei ist nicht im ASCII Format, Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1808)
  msg = "X!X FEHLER: levcfg darf nicht groesser als 2 sein."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1809)
  msg = "X!X FEHLER: Kann Parameter der Superzelle nicht einlesen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1810)
  msg = "X!X FEHLER: Kann Anzahl der Atome nicht einlesen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1811)
  !reals(1) = index of atom
  WRITE(msg,*) NINT(reals(1))
  msg = "X!X FEHLER: Atom index #"//TRIM(ADJUSTL(msg))// &
      & " ist groesser als Anzahl der Atome."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1812)
  !reals(1) = index of atom
  WRITE(msg,*) NINT(reals(1))
  msg = "X!X FEHLER: Kann Eigenschaften von Atom #"//TRIM(ADJUSTL(msg))&
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
    msg = "X!X FEHLER: Datei im Format "//TRIM(ADJUSTL(strings(1)))//  &
        & " könnten nur im Modus "//TRIM(ADJUSTL(strings(2)))//" gelesen werden."
  ELSE
    msg = "X!X FEHLER: Datei im Format "//TRIM(ADJUSTL(strings(1)))//" könnten nicht gelesen werden."
  ENDIF
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
  !strings(1) = disloctype: screw, edge, edge2
  !strings(2) = direction of dislocline: x, y or z
  !reals(1) = X component of Burgers vector
  !reals(2) = Y component of Burgers vector
  !reals(3) = Z component of Burgers vector
  !reals(4) = positive if C_tensor is defined, zero otherwise
  !reals(5) = pos1, position of dislocation along first axis
  !reals(6) = pos2, position of dislocation along second axis
  temp = TRIM(ADJUSTL(strings(1)))
  IF(TRIM(temp)=="screw") THEN
    msg = ">>> Fuege eine Schraubenversetzung ein, entlang der Linie"
  ELSEIF(temp(1:4)=="edge") THEN
    msg = ">>> Fuege eine Stufenversetzung ein, entlang der Linie"
  ENDIF
  msg = TRIM(msg)//' '//TRIM(strings(2))//","
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  !
  IF( reals(4)>0.1d0 ) THEN
    msg = "    wende anisotrope Elastizitaet an,"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
  !
  IF(TRIM(strings(1))=="edge") THEN
    WRITE(msg,"(a34)") "    durch Einfuegen einer atomaren Ebene,"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF(TRIM(strings(1))=="edge2") THEN
    WRITE(msg,"(a41)") "    Erhaltung der Anzahl der Atome,"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
  !
  WRITE(msg,"(f16.3)") reals(1)
  WRITE(temp,"(f16.3)") reals(2)
  WRITE(temp2,"(f16.3)") reals(3)
  msg = "["//TRIM(ADJUSTL(msg))//" "//TRIM(ADJUSTL(temp))//" "//TRIM(ADJUSTL(temp2))//"]"
  WRITE(temp,"(f16.3)") reals(5)
  WRITE(temp2,"(f16.3)") reals(6)
  msg = "    b="//TRIM(ADJUSTL(msg))//" at ("// &
    & TRIM(ADJUSTL(temp))//","//TRIM(ADJUSTL(temp2))//")"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2062)
  msg = "..> Bestimme Loesungen der anisotropen Elastizitaetsgleichungen..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2063)
  !reals(1) = number of inserted atoms
  WRITE(msg,*) NINT(reals(1))
  msg = "..> "//TRIM(ADJUSTL(msg))//" Atome wurden eingefuegt."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2064)
  msg = "..> Superzelle wurde erweitert."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2065)
  msg = "..> Versetzung wurde erfolgreich erzeugt."
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
  msg = "..> Das System wurde erfolgreich vervielfaeltigt."
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
  msg = ">>> Richte das System aus..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2072)
  msg = "..> System erfolgreich ausgerichtet."
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
    msg = ">>> Waehle alle Atome aus."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="index" ) THEN
    IF( strings(2)=="list " ) THEN
      msg = ">>> Waehle eine Liste von Atomen aus."
    ELSEIF( strings(2)=="range" ) THEN
      WRITE(temp,*) NINT(reals(1))
      WRITE(temp2,*) NINT(reals(2))
      msg = ">>> Waehle Atom #"//TRIM(ADJUSTL(temp))//" bis "//TRIM(ADJUSTL(temp2))//" aus."
    ELSE
      WRITE(temp,*) NINT(reals(1))
      msg = ">>> Waehle Atom #"//TRIM(ADJUSTL(temp))//" aus."
    ENDIF
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="invert" ) THEN
    msg = ">>> Invertiere die Auswahl der Atome..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="above" .OR. strings(1)=="below" ) THEN
    IF( DABS(reals(1))<1.d12 ) THEN
      WRITE(temp,"(3f16.3)") reals(1)
    ELSEIF( reals(1)<-1.d12 ) THEN
      WRITE(temp,"(a4)") "-INF"
    ELSEIF( reals(1)>1.d12 ) THEN
      WRITE(temp,"(a4)") "+INF"
    ENDIF
    msg = ">>> Waehle Atome "//TRIM(strings(1))//" "//TRIM(ADJUSTL(temp))// &
        & " A entlang der "//TRIM(strings(3))//" Achse aus."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="in" .OR. strings(1)=="out" ) THEN
    IF(strings(1)=="in") THEN
      temp = "innerhalb der"
    ELSE
      temp = "ausserhalb der"
    ENDIF
    msg = ">>> Waehle die Atome "//TRIM(temp)//" "//TRIM(strings(2))//"."
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
      msg = "..> Center: ("//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(2)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))//")"
      WRITE(temp,"(f16.3)") reals(4)
      msg = TRIM(msg)//"; Radius: "//TRIM(ADJUSTL(temp))//" A."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ENDIF
  ELSEIF( strings(1)=="prop" ) THEN
    msg = ">>> Waehle Atome mit "//TRIM(strings(2))
    WRITE(temp,"(f16.3)") reals(1)
    IF( reals(4)>2.d0 ) THEN
      WRITE(temp2,"(f16.3)") reals(2)
      msg = TRIM(ADJUSTL(msg))//" zwischen "//TRIM(ADJUSTL(temp))//" und "//TRIM(ADJUSTL(temp2))//"."
    ELSE
      msg = TRIM(ADJUSTL(msg))//" gleich "//TRIM(ADJUSTL(temp))//"."
    ENDIF
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="random" .OR. strings(1)=="rand" ) THEN
    WRITE(temp,*) NINT(reals(1))
    msg = ">>> Zufaellige Auswahl von "//TRIM(ADJUSTL(temp))//" Atomen"
    IF( strings(2).NE."any" .AND. strings(2).NE."all" ) THEN
      msg = TRIM(msg)//" von "//TRIM(strings(2))
    ENDIF
    msg = TRIM(msg)//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="random%" ) THEN
    WRITE(temp,'(f16.3)') reals(1)*100.d0
    msg = ">>> Zufaellige Auswahl von "//TRIM(ADJUSTL(temp))//"% der Atome"
    IF( strings(2).NE."any" .AND. strings(2).NE."all" ) THEN
      msg = TRIM(msg)//" von "//TRIM(strings(2))
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
      msg = "erste naechste"//TRIM(temp2)//" Nachbarn"
    ELSEIF( reals(1)>0.d0 .AND. DBLE(NINT(reals(1)))-reals(1)<1.d-12 ) THEN
      WRITE(temp,*) NINT(reals(1))
      IF( NINT(reals(1))==1 ) THEN
        msg = "der erste"//TRIM(temp2)//" Nachbar"
      ELSE
        msg = "die "//TRIM(ADJUSTL(temp))//" naehsten"//TRIM(temp2)//" Nachbarn"
      ENDIF
    ELSE
      WRITE(temp,'(f16.3)') DABS(reals(1))
      msg = TRIM(ADJUSTL(temp2))//" Nachbarn innerhalb eines Radius von "//TRIM(ADJUSTL(temp))//" A"
    ENDIF
    WRITE(temp,*) NINT(reals(2))
    msg = ">>> Waehle "//TRIM(ADJUSTL(msg))//" von Atom #"//TRIM(ADJUSTL(temp))//"..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSE
    !Last case: strings(1) should be an atom species
    msg = ">>> Waehle alle "//TRIM(ADJUSTL(strings(1)))//" Atome..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
CASE(2078)
  !reals(1) = number of atoms that were selected
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 Atom wurde ausgewaehlt."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(temp,*) NINT( reals(1) )
    msg = "..> "//TRIM(ADJUSTL(temp))//" Atome wurden ausgewaehlt."
  ELSE
    msg = "..> Kein Atom wurde ausgewaehlt."
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
  msg = ">>> Rotatiere das System um "//TRIM(ADJUSTL(msg)) &
      & //"° um die Achse "//TRIM(strings(1))
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
  !strings(2) = sort order: up, down, pack
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
  msg = "..> Zellenvektoren wurden berechnet."
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
  msg = "..> Anisotropieverhaeltnis: A = "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,'(e16.3)') reals(2)
  msg = "..> Anisotropiefaktor     : H = "//TRIM(ADJUSTL(msg))
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
  !reals(1) = max. displacement of an atom along a cartesian direction
  WRITE(msg,'(f24.3)') reals(1)
  msg = ">>> Zufaellige Verschiebung der Atome, max. Distanz"// &
      & TRIM(ADJUSTL(msg))//" A."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2115)
  msg = "..> Atome wurden aus ihrer Gleichgewichtslage verschoben."
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
  !reals(1) = index of first atom
  !reals(2) = index of second atom
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = ">>> Tauschen die Atomen #"//TRIM(ADJUSTL(temp))//" und #"//TRIM(ADJUSTL(temp2))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2126)
  msg = "..> Atomen wurden vertauschen."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2127)
  !strings(1) = roll axis: x, y or z
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
!
!2700-2799: WARNUNG MESSAGES
CASE(2700)
  !strings(1) = option name
  msg = "/!\ WARNUNG: folgende Anweisung konnte nicht interpretiert"// &
      & " werden: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Siehe `atomsk -help options` fuer eine Uebersicht"// &
      & " aller Optionen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2720)
  msg = "/!\ WARNUNG: Achse ist bereits ausgerichtet. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2721)
  msg = "/!\ WARNUNG: Koordinaten sind bereits reduziert. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2722)
  !strings(1) = atom species
  msg = "/!\ WARNUNG: "//TRIM(strings(1))//" Atome haben bereits Schalen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2723)
  !strings(1) = atomic species
  msg = "/!\ WARNUNG: kein "//TRIM(strings(1))//" Atom im System. "// &
      & "Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2724)
  msg = "/!\ WARNUNG: Deformation ist Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2725)
  msg = "/!\ WARNUNG: Burgers vector ist Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2726)
  msg = "/!\ WARNUNG: Superzelle hat eine sehr geringe Ausdehung"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    senkrecht zur Versetzungslinie. Das macht wenig Sinn."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Bitte ueberpruefen!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2727)
  !reals(1) = index of atom with large displacement
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ WARNUNG: Grosse Verschiebung fuer Atom #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2728)
  msg = "/!\ WARNUNG: Vergroesserungsfaktoren sind alle 1. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2729)
  msg = "/!\ WARNUNG: Keine Hilfseigenchaft definiert. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2730)
  !string(1) = property
  msg = "/!\ WARNUNG: Hilfseigenschaft "//TRIM(ADJUSTL(strings(1)))// &
      & " nicht gefunden. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2731)
  msg = "/!\ WARNUNG: Hstart = Hend. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2732)
  !strings(1) = name of unknown property
  msg = "/!\ WARNUNG: Unbekannte Eigenschaft: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2733)
  msg = "/!\ WARNUNG: Angegebener Radius ist negativ. "// &
      & "Kein Atom wird verschoben."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2734)
  msg = "/!\ WARNUNG: Rotationswinkel Modulo 2*Pi ist Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2735)
  msg = "/!\ WARNUNG: Scherung ist Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2736)
  msg = "/!\ WARNUNG: Verschiebevektor ist Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2737)
  msg = "/!\ WARNUNG: Atomsorten sind identisch. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2738)
  msg = "/!\ WARNUNG: Einheiten sind identisch. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2739)
  msg = "/!\ WARNUNG: Basisvektoren sind nicht orthonormal."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2740)
  msg = "/!\ WARNUNG: Elastischer Tensor ist nicht symmetrisch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2741)
  msg = "/!\ WARNUNG: Poisson Verhaeltnis ist ausserhalb des"// &
      & " Bereichs [-1 , 0.5]."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2742)
  msg = "/!\ WARNUNG: Ungueltiger Atomindex. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2743)
  msg = "/!\ WARNUNG: Verzerrungsparameter sind Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2744)
  msg = "/!\ WARNUNG: Keine ionischen Schalen im System. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2745)
  !reals(1) = number of atoms that will actually be removed
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ WARNUNG: kann so viele Atome nicht auswaehlen."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Es werden nur "//TRIM(ADJUSTL(msg))//" entfernt."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2746)
  msg = "/!\ WARNUNG: Nichts mehr zu tun. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2747)
  msg = "/!\ WARNUNG: Spannungsintensitaet K ist Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2748)
  msg = "/!\ WARNUNG: Superzelle ist sehr klein senkrecht zur "
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Risslinie. Das macht wenig Sinn. Bitte ueberpruefen!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2749)
  !strings(1) = name of property
  !reals(1) = atom index
  WRITE(temp,*) NINT(reals(1))
  msg = "/!\ WARNUNG: Konnte Wert der Eingenschaft '"// &
      & TRIM(ADJUSTL(strings(1)))//"' nicht dem Atom #"// &
      & TRIM(ADJUSTL(temp))//" zuweisen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2750)
  msg = "/!\ WARNUNG: Eine gesetzte Auswahl enthaehlt keine Atome mehr."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Die Auswahl wurde geloescht. "// &
      & "Alle Atome sind jetzt ausgewaehlt."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2751)
  msg = "/!\ WARNUNG: Zieltemperatur ist Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2752)
  msg = "/!\ WARNUNG: Keine Auswahl definiert. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2753)
  !reals(1) = atom index
  WRITE(temp,*) NINT(reals(1))
  msg = "/!\ WARNUNG: atom #"//TRIM(ADJUSTL(temp))//"  wurde bereits ausgewählt. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2754)
  !strings(1) = type of object (e.g. "dislocation" or "crack")
  msg = "/!\ WARNUNG: Das "//TRIM(ADJUSTL(strings(1)))//" wird aus der Box platziert."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Das macht wenig Sinn. Bitte ueberpruefen!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2755)
  msg = "/!\ WARNUNG: Faktor ist Null. Ueberspringe."
  CALL DISPLAY_MSG(1,msg,logfile)
  !
CASE(2799)
  !strings(1) = name of obsolete option
  !strings(2) = name of new option
  msg = "/!\ WARNUNG: Option "//TRIM(ADJUSTL(strings(1)))// &
      & " ist ueberholt und wird bald entfernt."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Bitte benutzen Sie stattdessen die Option "// &
      & TRIM(ADJUSTL(strings(2)))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!2800-2899: FEHLER MESSAGES
CASE(2800)
  !string(1) = axis (if we are here it's because it is different from x, y or z)
  msg = "X!X FEHLER: Unbekannte Achse: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2801)
  WRITE(msg,*) "X!X FEHLER: Die Basis Hend kann nicht auf"// &
             & " Hstart gedreht werden."
  CALL DISPLAY_MSG(1,msg,logfile)
  WRITE(msg,*) "    Ueberpruefen Sie bitte, ob die Winkel in"// &
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
  msg = "X!X FEHLER: Unbekannte Anweisung: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2806)
  !strings(1) = name of option
  msg = "X!X FEHLER: Unzulaessige Anweisungsform fuer "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2807)
  !reals(1) = 1 if roots of Eq.(13-85) cannot be found
  !         = 2 if the A_k(n) cannot be calculated
  !         = 3 if the linear equations defining D(n) cannot be solved
  msg = "X!X FEHLER:"
  IF(NINT(reals(1))==1) THEN
    msg = TRIM(msg)//" P(n) kann nicht bestimmt werden. Abbruch."
  ELSEIF(NINT(reals(1))==2) THEN
    msg = TRIM(msg)//" A_k(n) kann nicht bestimmt werden. Abbruch."
  ELSEIF(NINT(reals(1))==3) THEN
    msg = TRIM(msg)//" D(n) kann nicht bestimmt werden. Abbruch."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2808)
  msg = "X!X FEHLER: Kann keine gemischte Versetzung mit isotroper"// &
      & " Elastizitaet bilden."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2809)
  msg = "X!X FEHLER: Elastizitaetstensor enthaelt NaN Werte. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2810)
  msg = "X!X FEHLER: Inkonsistenz in den Datenlisten fuer Schalen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2811)
  msg = "X!X FEHLER: Versetzungsline und Versetzungsebene muessen"// &
      & " senkrecht aufeinander stehe. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2812)
  msg = "X!X FEHLER: Kann zu loeschende Atome nicht festlegen. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2813)
  !strings(1) = string that could not be converted to a number
  msg = "X!X FEHLER: Kann folgende Zeichenkette nicht in eine Zahl"// &
      & " umwandeln: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2814)
  msg = "X!X FEHLER: Kein atomares System, um diese Anweisung anzuwenden."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2815)
  !strings(1) = name of matrix to invert
  msg = "X!X FEHLER: Kann die Matrix "//strings(1)//" nicht invertieren."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2816)
  msg = "X!X FEHLER: Elastizitaetstensor nicht definiert. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!
!
!**************************
! 3000-3999: MESSAGES FOR OUTPUT
CASE(3000)
  !reals(1) = number of atoms
  WRITE(temp,*) NINT(reals(1))
  msg = ">>> Schreibe Ausgabedatei(en) ("//TRIM(ADJUSTL(temp))//" Atomen):"
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
!
!3700-3799: WARNUNG MESSAGES
CASE(3700)
  msg = "/!\ WARNUNG: kein Name fuer Ausgabedatei vergeben, bitte eingeben:"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3701)
  !strings(1) = name of output file
  msg = "/!\ WARNUNG: Format der Ausgabedatei konnte nicht erkannt werden: " &
      & //TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "<?> Welches Format der Ausgabedatei soll ich verwenden:"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = '    ('//TRIM(listofformats(1))
  DO i=2,SIZE(listofformats)
    msg = TRIM(msg)//','//TRIM(listofformats(i))
  ENDDO
  msg = TRIM(msg)//')'
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3702)
  msg = "/!\ WARNUNG: ddplot Format ist nur im DDPLOT Modus verfuegbar."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    dd Datei wird nicht geschrieben, siehe Dokumentation."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3703)
  msg = "/!\ WARNUNG: Atomsorten sind nicht zusammengefasst."// &
      & " Soll ich das jetzt machen? ("//langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    (Das hat nur Einfluss auf die POSCAR Datei)"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3704)
  msg = "/!\ WARNUNG: Die Zelle spannt keine untere Dreiecksmatrix auf."
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
  msg = "/!\ WARNUNG: Triklinischer Zellscherung "// &
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
  msg = "/!\ WARNUNG: Kann Zellscherung von "//TRIM(strings(1))// &
      & " nicht reduzieren. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3707)
  !reals(1) = number of atoms
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ WARNUNG: Atomanzahl ist sehr gross: "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Reduzieren die sie Anzahl mit der -cut Anweisung"// &
      & " von atomsk,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    falls ddplot nicht so viele Atome darstellen kann."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3708)
  msg = "/!\ WARNUNG: Nur die ersten 32 Hilfseigenschaften werden"// &
      & " in die CFG Datei geschrieben."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3709)
  msg = "/!\ WARNUNG: Einige Parameter der Superzelle koennen nicht"// &
      & " geschrieben werden."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Die Ausgabedatei wird Fehler enhalten!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3710)
  !strings(1) = file format
  msg = "/!\ WARNUNG: Unbekanntes Format '"//TRIM(strings(1))// &
      & "'. Ueberspringe "
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3711) ! missing occupancy data
  msg = "/!\ WARNUNG: Besetzungsdaten fehlen."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "             Verwende Standard Besetzungen von 1 fuer alle Atome."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3712) ! missing occupancy data for cel file output
  msg = "/!\ WARNUNG: Daten zur thermischen Vibration fehlen."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "             Setze die Biso Werte aller Atome auf Null."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3713) ! missing absorption data
  msg = "/!\ WARNUNG: absorption factors are missing, they will be set to 0.03 for all atoms."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!3800-3899: FEHLER MESSAGES
CASE(3800)
  msg = "X!X FEHLER: Keine Atomposition zu schreiben. Breche ab."
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
  msg = "X!X FEHLER: Atom #"//TRIM(ADJUSTL(msg))// &
      & " hat eine nicht-numerische (NaN) Koordinate. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3804)
  !reals(1) = index of atom that has a NaN coordinate
  WRITE(msg,*) NINT(reals(1))
  msg = "X!X FEHLER: Wert der Hilfseigenschaft von Atom #"//&
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
  msg = ">>> Atomsk ist eine freie Open Source software"// &
      & ", siehe 'atomsk --license'."
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
  msg = "help              Zeigt diese Hilfe"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = system_ls//"               Zeigt Dateien im aktuellen"// &
      & " Verzeichnis"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "memory            Zeigt Speicheruebersicht"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "create            Erstellt ein atomares System"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "read <file>       Datei <file> in den Speicher einlesen"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "write <file>      Schreibt aktuelles System in die Datei"// &
      & " <file>"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "clear             Loescht den Speicher"// &
      & " (zerstoert atomares System)"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "quit              Beendet atomsk"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "OPTIONS: atomsk Anweisungen fuer den Command-Line Interpreter,"
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
  msg = '    ('//TRIM(listofformats(1))
  DO i=2,SIZE(listofformats)
    msg = TRIM(msg)//','//TRIM(listofformats(i))
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
    msg = ">>> Verbinde die Systeme entlang "//TRIM(strings(1))//"..."
  ELSE
    msg = ">>> Verbinde die Systeme..."
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
  !reals(1) = width of the skin (A)
  msg = ">>> Berechne radiale Verteilungsfunktion (RDF),"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,'(f16.3)') reals(1)
  msg = "    innerhalb einer Schale der Breite "//TRIM(ADJUSTL(msg))//" A."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4052)
  !strings(1) = species 1
  !strings(2) = species 2
  msg = "..> Berechne RDF von "//TRIM(strings(2))//" Atomen um "// &
      & TRIM(strings(1))//" Atome..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4053)
  msg = "..> RDFs werden berechnet."
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
  msg = "/!\ WARNUNG: Datei nicht gefunden: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Ueberspringe..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4701)
  msg = "/!\ Konnte Basisvektoren der Superzelle in der"// &
      & " Eingabedatei nicht finden."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Versuche Basisvektoren zu erraten."// &
      & " Vorsicht, das kann sehr ungenau sein."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Es wird empfohlen die Basisvektoren der Superzelle"// &
      & " mit der Anweisungn -prop zu definieren."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4702)
  msg = "/!\ WARNUNG: Zelle ist geladen."// &
      & " Die totale Polarisation kann falsch sein!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4703)
  msg = "/!\ WARNUNG: Spezies hat keine Ladung, ueberspringe..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4705)
  !reals(1) = index of atom that has too many neighbors
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ WARNUNG: Anzahl der Nachbar ist groesser als 100"// &
      & " fuer Atom #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4706)
  !strings(1) = atomic species of absent neighbors
  !reals(1) = index of atom that has no neighbour
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ WARNUNG: Konnte keinen "//TRIM(strings(1))// &
      &       " Nachbarn fuer #"//TRIM(ADJUSTL(msg))//" finden."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4707)
  !reals(1) = index of atom that has a shell with zero charge
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ WARNUNG: Die Schale von Ion #"//TRIM(ADJUSTL(msg)) &
      & //" hat keine Ladung."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4708)
  msg = "/!\ WARNUNG: Dieses Korn enthaelt kein Atom."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4709)
  !reals(1) = index of atom
  !reals(2) = number of neighbors
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = "/!\ WARNUNG: Unzureichende Anzahl von Nachbarn ("// &
      & TRIM(ADJUSTL(temp2))//") fuer Atom #"//TRIM(ADJUSTL(temp))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4710)
  !strings(1) = name of the matrix
  msg = "/!\ WARNUNG: Konnte die Matrix "//TRIM(ADJUSTL(strings(1)))// &
      & " nicht berechnen."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4711)
  !strings(1) = name of the file
  msg = "/!\ WARNUNG: Abweichende Anzahl von Atomen in der Datei: "// &
      & TRIM(ADJUSTL(strings(1)))//". Datei wird nicht verwendet."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4712)
  msg = "/!\ WARNUNG: Einige Atome liegen ausserhalb der Zelle."// &
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
  msg = "/!\ WARNUNG: Offenbar wurde das System noch nicht in einer"// &
      & " Datei gespeichert."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Sind Sie sicher, dass Sie "// &
      & TRIM(ADJUSTL(strings(1)))//" wollen? ("// &
      & langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
!
!4800-4899: FEHLER MESSAGES
CASE(4800)
  !strings(1) = mode
  msg = "X!X FEHLER: Unzulaessiges Anweisung im Modus: "// &
      & TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4801)
  !strings(1) = file format (e.g. "xyz", "cfg"...)
  msg = "X!X FEHLER: Dieses Dateiformat wird noch nicht im 1-in-all"// &
      & " Modus unterstuetzt: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4802)
  msg = "X!X FEHLER: Beide chiralen Indizes sollten nicht Null sein."//&
      & " Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4803)
  !reals(1) = theoretical number of atoms in nanotube
  !reals(2) = number found by atomsk
  WRITE(msg,*) "X!X FEHLER: Unerwartete Anzahl von Atomen in der"// &
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
  msg = "X!X FEHLER: Diese Struktur benoetigt "//TRIM(ADJUSTL(temp))// &
      & " Elemente."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4805)
  !strings(1) = structure type
  msg = "X!X FEHLER: Unbekannte Struktur: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4806)
  msg = "X!X FEHLER: Keine Datei zum Zusammenfuegen. Abbruch"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4807)
  msg = "X!X FEHLER: Alle Ladungen sind Null. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Ladungen muessen mit der Anweisung -properties"// &
      & " gesetzt werden."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4808)
  !strings(1) = atomic species that has zero charge
  msg = "X!X FEHLER: "//TRIM(strings(1))//" Ionen haben keine Ladung."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4809)
  !strings(1) = atomic species
  msg = "X!X FEHLER: Kein "//TRIM(strings(1))//" Ion im System."
  msg = TRIM(ADJUSTL(msg))
CASE(4810)
  !reals(1) = number of particles in first system
  !reals(2) = number of particles in second system
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = "X!X FEHLER: Unterschiedliche Partikelzahl in beiden"// &
      & " Systemen: "//TRIM(ADJUSTL(temp))//"/"//TRIM(ADJUSTL(temp2))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "           Programm wird BEENDET."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4811)
  msg = "X!X FEHLER: Die beiden Systeme haben unterschiedliche"// &
      & " Basisvektoren."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4812)
  msg = "X!X FEHLER: Die Datei scheint kein animated XSF"// &
      & " zu enthalten. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4813)
  !strings(1) = name of unknown mode
  msg = "X!X FEHLER: Unbekannter Modus: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4814)
  msg = "X!X FEHLER: es kann immer nur ein Modus aktiv sein. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4815)
  msg = "X!X FEHLER: Datei scheint nicht im DL_POLY HISTORY Format"// &
      & " zu sein. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4816)
  !reals(1) = index for core charges
  !reals(2) = index for shell charges
  IF( reals(1)<reals(2) ) THEN
    msg = "X!X FEHLER: Kernladungen der Ionen nicht definiert. Abbruch."
  ELSE
    msg = "X!X FEHLER: Ladung der ionischen Schalen sind nicht"// &
        & " definiert. Abbruch."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4817)
  msg = "X!X FEHLER: Es gibt keine ionische Schale im System. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4818)
  !strings(1) = file format (e.g. "xyz", "cfg"...)
  msg = "X!X FEHLER: Die Dateiliste in "//TRIM(strings(1))// &
      & " scheint leer zu sein."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4819)
  msg = "X!X FEHLER: Basisvektoren sind nicht orthogonal. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4820)
  msg = "X!X FEHLER: Dimension der Superzelle wurden nicht"// &
      & " definiert. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4821)
  msg = "X!X FEHLER: Atomanzahl zu gross fuer bereitgestellten Speicher."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4822)
  msg = "X!X FEHLER: Keine Datei zu bearbeiten. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4823)
  !strings(1) = keyword 1
  !strings(2) = keyword 2
  msg = "X!X FEHLER: Die Schluesselwoerter '"// &
      & TRIM(ADJUSTL(strings(1)))//"' and '"// &
      & TRIM(ADJUSTL(strings(2)))//"' schliessen sich gegenseitig"// &
      & " aus. Abbruch."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4824)
  !strings(1) = name of unknown command
  msg = "X!X FEHLER: Unbekanntes Kommando: "//TRIM(ADJUSTL(strings(1)))
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
  msg = "/!\ WARNUNG: es ist notwendig die bestehende Datei "// &
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
  msg = "*** GLUeCKLICHES NEUES JAHR "//TRIM(msg)//"! ***"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!14 March (3/14): Pi day
ELSEIF(values(2)==3 .AND. values(3)==14) THEN
  msg = "               π"
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
  msg = "*** Sie arbeiten? Tag der Arbeit ist ein FEIERTAG! :-)"
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
  msg = " * * *  FROHE WEIHNACHTEN  !  * * *"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!If it's late at night or in the week end
ELSEIF( values(5)>=20 .OR. values(5)<=6 ) THEN
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
ENDIF
!
END SUBROUTINE DATE_MSG_DE
!
!
!
END MODULE messages_de
