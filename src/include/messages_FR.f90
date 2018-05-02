MODULE messages_FR
!
!**********************************************************************************
!*  MESSAGES_FR                                                                   *
!**********************************************************************************
!* This module contains the FRENCH                                                *
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
!* DISPLAY_LICENSE_FR  affiche la licence du programme                            *
!* DISPLAY_HELP_FR     affiche l'aide du programme                                *
!* ATOMSK_MSG_FR       tous les messages utilisés par le programme                *
!* DATE_MSG_FR         affiche un message en fonction de la date et l'heure       *
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
! of Atomsk
!********************************************************
SUBROUTINE DISPLAY_LICENSE_FR()
!
IMPLICIT NONE
CHARACTER(LEN=128):: msg
!
msg = ""
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "Atomsk - Un outil pour manipuler et convertir les fichiers de données atomiques."
CALL DISPLAY_MSG(verbosity,msg,logfile)
!
CALL DISPLAY_COPYRIGHT()
!
msg = ""
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "Ce programme est un logiciel libre : vous pouvez le redistribuer et/ou le modifier"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "selon les termes de la Licence Publique Générale GNU publiée par la"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "Free Software Foundation, soit la version 3 de la licence, soit"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "(à votre choix) n'importe quelle version supérieure."
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = ""
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "Ce programme est distribué dans l'espoir qu'il sera utile,"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "mais SANS AUCUNE GARANTIE ; sans même les garanties implicites de"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "QUALITÉ MARCHANDE ni d'APTITUDE POUR UN BUT PARTICULIER. Référez-vous"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "à la Licence Publique Générale GNU pour plus de détails."
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = ""
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "Vous devriez avoir reçu une copie de la Licence Publique Générale GNU"
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = "avec ce programme. Si ce n'est pas le cas, rendez-vous sur <http://www.gnu.org/licenses/>."
CALL DISPLAY_MSG(verbosity,msg,logfile)
msg = ""
CALL DISPLAY_MSG(verbosity,msg,logfile)
!
END SUBROUTINE DISPLAY_LICENSE_FR
!
!
!********************************************************
! DISPLAY_HELP
! This subroutine displays all or part of
! the help for Atomsk
!********************************************************
SUBROUTINE DISPLAY_HELP_FR(helpsection)
!
IMPLICIT NONE
CHARACTER(LEN=16):: helpsection
!
IF(TRIM(helpsection)=="") helpsection="general"
!
!
IF(helpsection=="general" .OR. helpsection=="") THEN
WRITE(*,*) ">>> UTILISATION :"
WRITE(*,*) "       atomsk [--<mode>] <inputfile> [<outputfile>] [<formats>] [options...]"
WRITE(*,*) "    où [] sont des paramètres optionnels, et <> doivent être remplacés"
WRITE(*,*) "    par des valeurs ou caractères spécifiques."
WRITE(*,*) ""
WRITE(*,*) ">>> EXEMPLE :"
WRITE(*,*) "..> Convertir 'file.xsf' au format CFG :"
WRITE(*,*) "       atomsk file.xsf cfg"
WRITE(*,*) ""
WRITE(*,*) ">>> POUR PLUS D'AIDE:"
WRITE(*,*) "..> sur les modes:   atomsk --help modes"
WRITE(*,*) "..> sur les options: atomsk --help options"
WRITE(*,*) "..> sur les formats: atomsk --help formats"
ENDIF
!
IF(helpsection=="modes") THEN
  WRITE(*,*) ">>> MODES:"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="interactive") THEN
  WRITE(*,*) "..> Mode interactif :"
  WRITE(*,*) "          atomsk"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="normal") THEN
  WRITE(*,*) "..> Mode normal :"
  WRITE(*,*) "          atomsk <inputfile> [<outputfile>] [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="list") THEN
  WRITE(*,*) "..> Mode liste:"
  WRITE(*,*) "          atomsk --list <listfile> [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="ai1") THEN
  WRITE(*,*) "..> Mode tout-en-un:"
  WRITE(*,*) "          atomsk --all-in-one <listfile> <outputfile> [options] [<formats>]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="1ia") THEN
  WRITE(*,*) "..> Mode un-en-tout :"
  WRITE(*,*) "          atomsk --one-in-all <inputfile> [<outputfile>] [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="create") THEN
  WRITE(*,*) "..> Mode création :"
  WRITE(*,*) "          atomsk --create <structure> <a0> <species> <outputfile> [orient hkl hkl hkl] [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="ddplot") THEN
  WRITE(*,*) "..> Mode DDplot :"
  WRITE(*,*) "          atomsk --ddplot <file1> <file2> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="merge") THEN
  WRITE(*,*) "..> Mode fusion :"
  WRITE(*,*) "          atomsk --merge [<x|y|z>] <Nfiles> <file1>...<fileN> <outputfile> [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="unwrap") THEN
  WRITE(*,*) "..> Mode déballage:"
  WRITE(*,*) "          atomsk --unwrap <reference> <system> [<outputfile>] [<formats>] [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="density") THEN
  WRITE(*,*) "..> Mode densité :"
  WRITE(*,*) "          atomsk --density <fichier> <propriété> <dimension> <axe> <sigma> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="difference") THEN
  WRITE(*,*) "..> List mode:"
  WRITE(*,*) "          atomsk --difference <file1> <file2> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="edm") THEN
  WRITE(*,*) "..> Mode moments dipolaires électriques :"
  WRITE(*,*) "          atomsk --edm <system> <Pspecies> <NNN> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="PE") THEN
  WRITE(*,*) "..> Mode polarisation électronique :"
  WRITE(*,*) "          atomsk --electronic-polarization <system> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="rdf") THEN
  WRITE(*,*) "..> Mode fonction de distribution radiale :"
  WRITE(*,*) "          atomsk --rdf <listfile> <Rmax> <dR> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="polycrystal") THEN
  WRITE(*,*) "..> Mode polycristal :"
  WRITE(*,*) "          atomsk --polycrystal <cellule> <paramètres> <outputfile> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="average") THEN
  WRITE(*,*) "..> Mode moyenne :"
  WRITE(*,*) "          atomsk --average <listfile> <outputfile> [options]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="nye") THEN
  WRITE(*,*) "..> Mode tenseur de Nye:"
  WRITE(*,*) "          atomsk --nye <reference> <dislocation> <outputfile> [options] [<formats>]"
ENDIF
IF(helpsection=="modes" .OR. helpsection=="interpolate") THEN
  WRITE(*,*) "..> Mode Interpolate:"
  WRITE(*,*) "          atomsk --interpolate <fichier1> <fichier2> <N> [options] [<formats>]"
ENDIF
!
IF(helpsection=="options") THEN
  WRITE(*,*) ">>> OPTIONS (distances=Angströms, angles=degrés):"
  WRITE(*,*) ""
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-add-atom" .OR. helpsection=="-add-atoms" .OR. &
  &helpsection=="-addatom" .OR. helpsection=="-addatoms" ) THEN
  WRITE(*,*) "..> Ajouter de nouveaux atomes au système :"
  WRITE(*,*) "          -add-atom <espèce> at <x> <y> <z>"
  WRITE(*,*) "          -add-atom <espèce> relative <indice> <x> <y> <z>"
  WRITE(*,*) "          -add-atom <espèce> near <indice>"
  WRITE(*,*) "          -add-atom <espèce> random <N>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-add-shells".OR. helpsection=="-as" &
  & .OR. helpsection=="-create-shells".OR. helpsection=="-cs") THEN
  WRITE(*,*) "..> Ajouter des coquilles à une ou toutes les espèces :"
  WRITE(*,*) "          -add-shells <all|species>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-alignx") THEN
  WRITE(*,*) "..> Aligner le premier vecteur avec l'axe X :"
  WRITE(*,*) "          -alignx"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-bind-shells" .OR. helpsection=="-bs") THEN
  WRITE(*,*) "..> Réassocier les coquilles avec leurs cœurs :"
  WRITE(*,*) "          -bind-shells"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-center") THEN
  WRITE(*,*) "..> Placer un atome ou le centre de masse du système au centre de la boîte :"
  WRITE(*,*) "          -center <index>"
  WRITE(*,*) "          -center com"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-crack") THEN
  WRITE(*,*) "..> Insérer une fracture dans le système :"
  WRITE(*,*) "          -crack <I/II/III> <stress/strain> <K> <pos1> <pos2> "//&
           &            "<crackline> <crackplane> <μ> <ν>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-cut") THEN
  WRITE(*,*) "..> Couper une partie du système :"
  WRITE(*,*) "          -cut <above|below> <cutdistance> <x|y|z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-deform" .OR. helpsection=="-def") THEN
  WRITE(*,*) "..> Appliquer une déformation ou contrainte uniaxiale :"
  WRITE(*,*) "          -def <x|y|z> <strain> <Poissons ratio>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-dislocation" .OR. helpsection=="-disloc") THEN
  WRITE(*,*) "..> Insérer une dislocation dans le système :"
  WRITE(*,*) "          -disloc <pos1> <pos2> screw <x|y|z> <x|y|z> <b>"
  WRITE(*,*) "          -disloc <pos1> <pos2> <edge|edge2> <x|y|z> <x|y|z> <b> <ν>"
  WRITE(*,*) "          -disloc <pos1> <pos2> mixed <x|y|z> <x|y|z> <b1> <b2> <b3>"
  WRITE(*,*) "          -disloc loop <x> <y> <z> <x|y|z> <rayon> <bx> <by> <bz> <nu>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-disturb" ) THEN
  WRITE(*,*) "..> Déplacer aléatoirement les atomes :"
  WRITE(*,*) "          -disturb <dmax>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-duplicate" .OR. helpsection=="-dup") THEN
  WRITE(*,*) "..> Dupliquer le système dans les 3 directions de l'espace :"
  WRITE(*,*) "          -duplicate <Nx> <Ny> <Nz>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-fix") THEN
  WRITE(*,*) "..> Fixer des atomes:"
  WRITE(*,*) "          -fix <x|y|z> <above|below> <distance> <x|y|z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-fractional" .OR. helpsection=="-frac") THEN
  WRITE(*,*) "..> Convert coordinates to fractional :"
  WRITE(*,*) "          -frac"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-mirror") THEN
  WRITE(*,*) "..> Appliquer un plan miroir :"
  WRITE(*,*) "          -mirror <d> <normal>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-options") THEN
  WRITE(*,*) "..> Appliquer une liste d'options lue depuis un fichier :"
  WRITE(*,*) "          -options <fichier>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-orient") THEN
  WRITE(*,*) "..> Changer l'orientation cristallographique du système:"
  WRITE(*,*) "          -orient <Hx> <Hy> <Hz> <H'x> <H'y> <H'z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-prop") THEN
  WRITE(*,*) "..> Lire les propriétés du système :"
  WRITE(*,*) "          -prop <fichier>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-rebox") THEN
  WRITE(*,*) "..> (Re-)calculer les vecteurs de la boîte :"
  WRITE(*,*) "          -rebox"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-remove-atom" .OR. helpsection=="-rmatom") THEN
  WRITE(*,*) "..> Supprimer un atome d'un indice donné, ou tous les atomes d'une espèce :"
  WRITE(*,*) "          -rmatom <index|species>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-remove-doubles" .OR. helpsection=="-rmd") THEN
  WRITE(*,*) "..> Supprimer les atomes en double :"
  WRITE(*,*) "          -rmd <distance>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-remove-property" .OR. helpsection=="-rmprop") THEN
  WRITE(*,*) "..> Supprimer une ou toutes les propriétés auxiliaires :"
  WRITE(*,*) "          -rmprop <property>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-remove-shell" .OR. helpsection=="-rmshell" .OR. &
  &helpsection=="-remove-shells" .OR. helpsection=="-rmshells" ) THEN
  WRITE(*,*) "..> Supprimer les coquilles sur une espèce chimique, ou sur tous les atomes :"
  WRITE(*,*) "          -rmshells <espèce|all>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-roll" .OR. helpsection=="-bend") THEN
  WRITE(*,*) "..> Enrouler le système autour d'un axe:"
  WRITE(*,*) "          -roll <x|y|z> <angle> <x|y|z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-rotate" .OR. helpsection=="-rot") THEN
  WRITE(*,*) "..> Tourner le système autour d'un axe :"
  WRITE(*,*) "          -rotate [com] <x|y|z> <angle>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-select") THEN
  WRITE(*,*) "..> Selectionner les atomes selon un critère :"
  WRITE(*,*) "          -select all"
  WRITE(*,*) "          -select invert"
  WRITE(*,*) "          -select <species>"
  WRITE(*,*) "          -select <index>"
  WRITE(*,*) "          -select <above|below> <d> <normal>"
  WRITE(*,*) "          -select <in|out> <box|sphere|cylinder|torus> [<axe>] <x1> <y1> <z1> <x2> [<y2> <z2>]]"
  WRITE(*,*) "          -select prop <prop> <value>"
  WRITE(*,*) "          -select random <N> <species>"
  WRITE(*,*) "          -select <NNN> <species> neighbors <index>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-separate") THEN
  WRITE(*,*) "..> Séparer les atomes qui sont trop proches:"
  WRITE(*,*) "          -separate <distance> <décalage>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-shear") THEN
  WRITE(*,*) "..> Appliquer un cisaillement au système :"
  WRITE(*,*) "          -shear <x|y|z> <amplitude> <x|y|z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-shift") THEN
  WRITE(*,*) "..> Translater une partie du système :"
  WRITE(*,*) "          -shift <above|below> <distance> <x|y|z> <tauX> <tauY> <tauZ>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-sort") THEN
  WRITE(*,*) "..> Ranger les atomes :"
  WRITE(*,*) "          -sort <s|x|y|z> <up|down|pack>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-spacegroup") THEN
  WRITE(*,*) "..> Appliquer les opérations de symétrie d'un groupe d'espace donné:"
  WRITE(*,*) "          -spacegroup <groupe>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-stress") THEN
  WRITE(*,*) "..> Appliquer une contrainte :"
  WRITE(*,*) "          -stress <xx|yy|zz|xy|xz|yz|P> <valeur(GPa)>"
  WRITE(*,*) "          -stress <file>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-substitute" .OR. helpsection=="-sub") THEN
  WRITE(*,*) "..> Substituer les atomes sp1 par des atomes sp2 :"
  WRITE(*,*) "          -sub <sp1> <sp2>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-swap") THEN
  WRITE(*,*) "..> Échanger les indices de deux atomes ou des axes cartésiens :"
  WRITE(*,*) "          -swap <id1> <id2>"
  WRITE(*,*) "          -swap <x|y|z> <x|y|z>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-torsion") THEN
  WRITE(*,*) "..> Appliquer une torsion autour d'un axe :"
  WRITE(*,*) "          -torsion <x|y|z> <angle>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-unit" .OR. helpsection=="-u") THEN
  WRITE(*,*) "..> Changer l'unité de distance :"
  WRITE(*,*) "          -u <unit1> <unit2>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-unskew") THEN
  WRITE(*,*) "..> Réduire l'inclinaison de la boîte :"
  WRITE(*,*) "          -unskew"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-velocity") THEN
  WRITE(*,*) "..> Donner une vitesse aléatoire aux atomes selon une distribution de Maxwell-Boltzmann :"
  WRITE(*,*) "          -velocity <T>"
ENDIF
!
IF(helpsection=="options" .OR. helpsection=="-wrap") THEN
  WRITE(*,*) "..> Ramener les atomes dans la boîte :"
  WRITE(*,*) "          -wrap"
ENDIF
!
! -- add other options in alphabetical order --
!
!
IF(helpsection=="formats") THEN
WRITE(*,*) ">>> FORMATS :"
WRITE(*,*) "    Les formats doivent être spécifiés tels que dans la première colonne."
WRITE(*,*) "    Atomusk peut convertir de n'importe quel 'oui' dans la colonne ENTRÉE"
WRITE(*,*) "    vers n'importe quel 'oui' dans la colonne SORTIE :"
WRITE(*,*) "                            | ENTRÉE | SORTIE"
WRITE(*,*) "    ------------------------+--------+--------"
WRITE(*,*) "    atsk (format Atomsk)    |   oui  |  oui"
WRITE(*,*) "    bop                     |   oui  |  oui"
WRITE(*,*) "    cel (Dr. Probe/EMS)     |   oui  |  oui"
WRITE(*,*) "    coo (COORAT/MBPP)       |   oui  |  oui"
WRITE(*,*) "    cfg (Atomeye)           |   oui  |  oui"
WRITE(*,*) "    cif (Cryst.Info.File)   |   oui  |  oui"
WRITE(*,*) "    dd  (ddplot)            |   non  | oui (1)"
WRITE(*,*) "    dlp (DL_POLY CONFIG)    |   oui  |  oui"
WRITE(*,*) "    gin (fichier GULP)      |   oui  |  oui"
WRITE(*,*) "    imd (fichier IMD)       |   oui  |  oui"
WRITE(*,*) "    jems (fichier JEMS)     |   oui  |  oui"
WRITE(*,*) "    lmc (LAMMPS output)     |   oui  |  non"
WRITE(*,*) "    lmp (LAMMPS data)       |   oui  |  oui"
WRITE(*,*) "    mol (fichier MOLDY)     |   oui  |  oui"
WRITE(*,*) "    pdb (Protein Data Bank) |   oui  |  oui"
WRITE(*,*) "    pos (POSCAR/VASP)       |   oui  |  oui"
WRITE(*,*) "    pw (Quantum Espresso)   |   oui  |  oui"
WRITE(*,*) "    pwout (sortie QE)       |  oui(2)|  non"
WRITE(*,*) "    str (PDFFIT)            |   oui  |  oui"
WRITE(*,*) "    vesta (fichier VESTA)   |   oui  |  oui"
WRITE(*,*) "    xsf (XCrySDen)          |   oui  |  oui"
WRITE(*,*) "    xv (fichier SIESTA)     |   oui  |  oui"
WRITE(*,*) "    xyz/exyz/sxyz           |   oui  |  oui"
WRITE(*,*) "        (1) Mode ddplot seulement."
WRITE(*,*) "        (2) Mode one-in-all seulement."
ENDIF
!
WRITE(*,*) ""
WRITE(*,*) ">>> Référez-vous à la documentation dans le dossier /doc fourni avec le programme"
WRITE(*,*) "    ou à cette adresse : http://atomsk.univ-lille1.fr/fr/"
WRITE(*,*) ""
!
!
END SUBROUTINE DISPLAY_HELP_FR
!
!
!
!********************************************************
! ATOMSK_CREATE_DATE
! This routine 
!********************************************************
SUBROUTINE ATOMSK_CREATE_DATE_FR(VALUES,username,msg)
!
IMPLICIT NONE
CHARACTER(LEN=128),INTENT(IN):: username
INTEGER,DIMENSION(8),INTENT(IN):: VALUES
CHARACTER(LEN=128),INTENT(OUT):: msg
!
WRITE(msg,'(i4,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2,a1,i2.2)') &
  & VALUES(1), "-", VALUES(2),"-", VALUES(3)," ", VALUES(5), ":", VALUES(6), ":", VALUES(7)
!
msg = 'Fichier généré avec Atomsk par '//TRIM(ADJUSTL(username))//' le '//TRIM(ADJUSTL(msg))
!
END SUBROUTINE ATOMSK_CREATE_DATE_FR
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
SUBROUTINE ATOMSK_MSG_FR(imsg,strings,reals)
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
  msg = "\o/ Le programme s'est terminé avec succès !"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(temp,"(f30.3)") reals(1)
  WRITE(temp2,"(f30.3)") reals(2)
  WRITE(msg,*) "   Temps total : "//TRIM(ADJUSTL(temp))//         &
           & " s.; temps CPU : "//TRIM(ADJUSTL(temp2))//" s."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  IF( nwarn>0 .OR. nerr>0 ) THEN
    !In case of warning or error, display a big box
    msg = " _________________________________________"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    IF( nwarn>0 ) THEN
      WRITE(temp,*) nwarn
      temp = ADJUSTL(temp)
      msg = "|  /!\ ALERTES : "//TRIM(temp)
      msg = msg(1:42)//"|"
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ENDIF
    IF( nerr>0 ) THEN
      WRITE(temp,*) nerr
      temp = ADJUSTL(temp)
      msg = "|  X!X ERREURS : "//TRIM(temp)
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
  WRITE(msg,*) "   Temps total : "//TRIM(ADJUSTL(temp))//         &
           & " s.; temps CPU : "//TRIM(ADJUSTL(temp2))//" s."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3)
  msg = "..> Ceci peut prendre du temps..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4)
  !strings(1) = name of file that already exists
  msg = "<?> Ce fichier existe déjà : "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Voulez-vous l'écraser ("//langyes//"/"//langno//") ("&
      & //langBigYes//"=tout écraser)?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5)
  !strings(1) = name of file
  msg = "..> Écrasement du fichier : "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(6)
  msg = "..> OK j'écraserai tous les fichiers."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(7)
  msg = "<?> Entrez le nom du fichier à écrire:"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(8)
  !strings(1) = name of file
  IF( LEN_TRIM(strings(1))>0 ) THEN
    msg = "<?> Ce fichier n'existe pas : "//TRIM(strings(1))
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
  msg = "    Veuillez fournir le nom d'un fichier existant :"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    (tapez '"//TRIM(system_ls)//"' pour une liste des fichier du dossier courant)"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(9)
  msg = "<?> Entrez le nom d'un fichier existant :"
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
        temp = TRIM(ADJUSTL(temp2))//" % [====================]"
      ELSEIF(tempreal>=95.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [===================>]"
      ELSEIF(tempreal>=90.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [==================> ]"
      ELSEIF(tempreal>=85.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [=================>  ]"
      ELSEIF(tempreal>=80.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [================>   ]"
      ELSEIF(tempreal>=75.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [===============>    ]"
      ELSEIF(tempreal>=70.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [==============>     ]"
      ELSEIF(tempreal>=65.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [=============>      ]"
      ELSEIF(tempreal>=60.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [============>       ]"
      ELSEIF(tempreal>=55.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [===========>        ]"
      ELSE
        temp = TRIM(ADJUSTL(temp2))//" % [==========>         ]"
      ENDIF
    ELSE  
      IF(tempreal<5.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [>                   ]"
      ELSEIF(tempreal<10.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [=>                  ]"
      ELSEIF(tempreal<15.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [==>                 ]"
      ELSEIF(tempreal<20.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [===>                ]"
      ELSEIF(tempreal<25.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [====>               ]"
      ELSEIF(tempreal<30.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [=====>              ]"
      ELSEIF(tempreal<35.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [======>             ]"
      ELSEIF(tempreal<40.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [=======>            ]"
      ELSEIF(tempreal<45.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [========>           ]"
      ELSEIF(tempreal<50.d0) THEN
        temp = TRIM(ADJUSTL(temp2))//" % [=========>          ]"
      ELSE
        temp = TRIM(ADJUSTL(temp2))//" % [==========>         ]"
      ENDIF
    ENDIF
    !Display the progress bar
    WRITE(*,'(a)',ADVANCE="NO") CHAR(13)//" "//TRIM(temp)
    IF( NINT(reals(1)) == NINT(reals(2)) ) THEN
      WRITE(*,*) ""
    ENDIF
  ENDIF
CASE(11)
  msg = ">>> Construction de la liste de voisins..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(12)
  msg = ">>> Si vous utilisez Atomsk dans vos travaux, merci de citer l'article suivant :"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "    Pierre Hirel, Comput. Phys. Comm. 197 (2015) 212"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(14)
  msg = "<i> Les atomes sélectionnés ont été supprimés : aucune sélection n'est plus définie."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
! 700- 799: WARNING MESSAGES
CASE(700)
  !strings(1) = first file name
  !strings(2) = second file name
  msg = "/!\ Ces deux fichiers existent: "//TRIM(strings(1))//" et "//TRIM(strings(2))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "<?> Veuillez indiquer lequel vous souhaitez utiliser en entrée :"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = '    1- '//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = '    2- '//TRIM(strings(2))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(701)
  !strings(1) = first file name
  !strings(2) = second file name
  msg = "/!\ Aucun de ces fichiers n'existe : "//TRIM(strings(1))//" ni "//TRIM(strings(2))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "<?> Veuillez fournir le nom d'un fichier existant :"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(702)
  msg = "/!\ Aucun fichier d'entrée spécifié."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "<?> Veuillez fournir le nom d'un fichier existant :"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(703)
  !strings(1) = name of command line argument
  msg = "/!\ ALERTE : argument inconnu dans la ligne de commande : "//TRIM(strings(1))
  nwarn = nwarn+1
CASE(704)
  !strings(1) = name of file that will be ignored
  msg = "/!\ ALERTE : trop de noms de fichiers parmi les arguments."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "            Le fichier suivant sera ignoré : "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(750)
  msg = ""
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "/!\ ALERTE : VOUS NE DEVRIEZ PAS EXÉCUTER CE PROGRAMME EN TANT QUE ROOT !"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  temp=""
  DO WHILE(temp.NE."ok")
    msg = "    Entrez 'ok' pour continuer malgré tout, ou Ctrl+C pour quitter."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    READ(*,*) temp
  ENDDO
  msg = ""
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
! 800- 899: ERROR MESSAGES
CASE(800)
  msg = "X!X ERREUR : format de fichier inconnu."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(801)
  !strings(1) = atom type
  msg = "X!X ERREUR : espèce atomique inconnue : "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(802)
  !reals(1) = index of atom that caused the error
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = "X!X ERREUR en tentant de lire l'atome #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(803)
  !strings(1) = unit
  msg = "X!X ERREUR : unité inconnue: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(804)
  msg = "X!X ERREUR : le nombre d'atomes est zéro, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(806)
  msg = "X!X ERREUR : les systèmes n'ont pas le même nombre d'atomes !"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(807)
  !reals(1) = index of line in file that caused the error on read
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = "X!X ERREUR en essayant de lire la ligne #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(808)
  !strings(1) = string where conversion failed
  msg = "X!X ERREUR : impossible de convertir l'expression '"//TRIM(strings(1))// &
      & "' en valeur numérique."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(809)
  !strings(1) = space group H-M symbol which couldn't be identified
  msg = "X!X ERREUR : groupe d'espace non reconnu : '"// &
      & TRIM(strings(1))//"'."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(810)
  !reals(1) = space group number
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = "X!X ERREUR : groupe d'espace invalide : "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(811)
  !reals(1) = space group number
  WRITE(msg,"(i18)") NINT(reals(1))
  msg = "X!X ERREUR : impossible d'accéder aux données du groupe d'espace : "// &
      & TRIM(ADJUSTL(msg))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(812)
  !strings(1) = non-conformal symmetry operation string
  msg = "X!X ERREUR : impossible d'appliquer les opérations de symétrie : '"// &
      & TRIM(strings(1))//"'."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(813)
  msg = "X!X ERREUR : aucun fichier ou dossier de ce type"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(814)
  msg = "X!X ERREUR : un vecteur de Miller ne peut pas être [000]."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(815)
  msg = "X!X ERREUR : dans les mailles hexagonales, les indices de Miller [hkil] doivent respecter h+k+i=0."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(816)
  msg = "X!X ERREUR : division par zero."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(817)
  !strings(1) = a string that could not be interpreted
  msg = "X!X ERREUR : impossible de convertir cette chaîne en indices de Miller : "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(818)
  !strings(1) = name of the array that was supposed to be resized
  msg = "X!X ERREUR : une erreur est survenue en tentant de redimensionner le tableau "//TRIM(ADJUSTL(strings(1)))
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
  msg = ">>> Ouverture du fichier : "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(1001)
  !reals(1) = number of atoms
  !reals(2) = number of shells
  WRITE(temp,*) NINT(reals(1))
  IF( NINT(reals(2))>0 ) THEN
    WRITE(temp2,*) NINT(reals(2))
    msg = "..> Le fichier a bien été lu ("//TRIM(ADJUSTL(temp))//" cœurs + "// &
        & TRIM(ADJUSTL(temp2))//" coquilles)."
  ELSE
    msg = "..> Le fichier a bien été lu ("//TRIM(ADJUSTL(temp))//" atomes)."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(1002)
  msg = "..> Fichier INP trouvé, lecture des vecteurs..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(1003)
  msg = "..> Fichier POTCAR trouvé, lecture des espèces atomiques..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(1004)
  msg = "..> Application des opérations de symétrie..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!1700-1799: WARNING MESSAGES
CASE(1700)
  !strings(1) = auxiliary property that cannot be loaded
  msg = "/!\ ALERTE : impossible de lire la propriété auxiliaire : "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1701)
  !strings(1) = input file format
  msg = "/!\ ALERTE : le format de fichier le plus probable est "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Cependant il n'a pas pu être déterminé de façon certaine."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1702)
  !strings(1) = name of parameter
  !strings(2) = name of custom config file
  msg = "/!\ ALERTE : paramètre inconnu '"//TRIM(ADJUSTL(strings(1)))//&
      & "' dans "//TRIM(strings(2))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1703)
  msg = "/!\ ALERTE : des erreurs sont survenues lors de la lecture du fichier "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Certains paramètres personnels n'ont pas pu être chargés."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1704)
  !strings(1) = name of personal config file
  msg = "/!\ ALERTE : les opérations de symétrie ne sont pas prises en charge."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1705)
  msg = "/!\ ALERTE : les deux notations celldm(:) et conventionnelle ont été trouvées."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Les celldm(:) seront utilisés, et la notation conventionnelle ignorée."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1706)
  msg = "/!\ ALERTE : les paramètres de boîte sont en Bohrs, alors que les positions atomiques sont en angströms."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Les vecteurs de boîte seront converties en angströms pour rétablir la cohérence."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1707) ! invalid symmetry operation string input.
  !strings(1) = failed symmetry operation string
  msg = "/!\ ALERTE : opération de symétrie invalide : '"// &
      & TRIM(strings(1))//"', ignorée..."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!1800-1899: ERROR MESSAGES
CASE(1800)
  msg = "X!X ERREUR: je n'ai pas pu déterminer le format de ce fichier !"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    C'est peut-être un format non supporté par Atomsk ?"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    En tout cas je ne peux apparemment rien faire, j'abandonne."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1801)
  !strings(1) = file name
  !reals(1) = line number
  msg = "X!X ERREUR : il y a eu des erreurs en lisant le fichier : " &
      & //TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
  IF( NINT(reals(1))>0 ) THEN
    WRITE(temp,*) NINT(reals(1))
    msg = "          L'erreur semble se trouver à la ligne # "//TRIM(ADJUSTL(temp))
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(1802)
  !strings(1) = bad array
  msg = "X!X ERREUR : problème de dimension dans la matrice "//TRIM(strings(1))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1803)
  msg = "X!X ERREUR : la dimension des propriétés auxiliaires "// &
             & "n'est pas cohérente avec le nombre d'atomes."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1804)
  msg = "X!X ERREUR : format inconnu."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1805)
  !reals(1) = number of particles read
  !reals(2) = number of particles declared
  msg = "X!X ERREUR : le nombre d'atomes trouvés est différent"
  CALL DISPLAY_MSG(1,msg,logfile)
  WRITE(temp,*) NINT(reals(1))
  WRITE(msg,*) NINT(reals(2))
  msg = "           du nombre d'atomes déclarés : "// &
      & TRIM(ADJUSTL(temp))//"/"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1806)
  !reals(1) = index of atom
  WRITE(msg,*) NINT(reals(1))
  msg = "X!X ERREUR : division par zéro dans la position de l'atome #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1807)
  msg = "X!X ERREUR : le fichier n'est pas au format ASCII."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1808)
  msg = "X!X ERREUR : levcfg ne peut pas être plus grand que 2."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1809)
  msg = "X!X ERREUR : impossible de lire les vecteurs de la boîte."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1810)
  msg = "X!X ERREUR : impossible de lire le nombre d'atomes."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1811)
  !reals(1) = index of atom
  WRITE(msg,*) NINT(reals(1))
  IF( reals(1)<0.1d0 ) THEN
    msg = "X!X ERREUR : l'indice de l'atome #"//TRIM(ADJUSTL(msg))//" est négatif."
  ELSE
    msg = "X!X ERREUR : l'indice de l'atome #"//TRIM(ADJUSTL(msg))//" est plus grand que le nombre d'atomes."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1812)
  !reals(1) = index of atom
  WRITE(msg,*) NINT(reals(1))
  msg = "X!X ERREUR : impossible de lire les propriétés de l'atome #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1813)
  !strings(1) = string where conversion failed
  msg = "X!X ERREUR : échec lors de la lecture des opérations de symétrie dans '"// &
      & TRIM(strings(1))//"'."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(1814)
  !strings(1) = name of file format
  !strings(2) = name of mode to use to read this file format
  IF( LEN_TRIM(strings(2))>0 ) THEN
    msg = "X!X ERREUR : les fichiers au format "//TRIM(ADJUSTL(strings(1)))//  &
        & " ne peuvent être lus que dans le mode "//TRIM(ADJUSTL(strings(2)))//"."
  ELSE
    msg = "X!X ERREUR : les fichiers au format "//TRIM(ADJUSTL(strings(1)))//" ne peuvent pas être lus."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
!
!
!
!**************************
! 2000-2999: MESSAGES FOR OPTIONS
!
CASE(2050)
  msg = ">>> Alignement des vecteurs de boîte..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2051)
  msg = "..> Les vecteurs ont bien été alignés."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2052)
  msg = ">>> Conversion en coordonnées réduites..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2053)
  msg = "..> Les coordonnées ont été réduites."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2054)
  !strings(1) = atom species
  IF(strings(1)=="all") THEN
    temp = "tous les atomes"
  ELSE
    temp = "les atomes de "//TRIM(strings(1))
  ENDIF
  msg = ">>> Création de coquilles pour "//TRIM(ADJUSTL(temp))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2055)
  !reals(1) = number of shells
  IF( NINT(reals(1))==0 ) THEN
    msg = "..> Aucune coquille n'a été ajoutée."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 coquille a été ajoutée au système."
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" coquilles ont été ajoutées au système."
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
  IF(strings(1)=="above") THEN
    temp = "au-dessus de "
  ELSE
    temp = "en-dessous de "
  ENDIF
  msg = ">>> Troncation du système "//TRIM(temp)//" "// &
    & TRIM(ADJUSTL(msg))//" A suivant "//TRIM(strings(2))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2057)
  !reals(1) = NPcut, number of deleted atoms
  !reals(2) = number of atoms left
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 atome a été supprimé"
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" atomes ont été supprimés"
  ENDIF
  WRITE(temp,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(msg))//", "//TRIM(ADJUSTL(temp))//" atomes restants."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2058)
  !reals(1) = deformation in %
  !strings(1) = direction of deformation: x, y or z
  WRITE(msg,"(f16.3)") reals(1)*100.d0
  msg = ">>> Déformation du système de "//TRIM(ADJUSTL(msg))//"% suivant "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2059)
  !reals(1) = Poisson coefficient
  IF(reals(1)==0.d0) THEN
    msg = "..> Déformation uniaxiale."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSE
    WRITE(msg,"(f16.3)") reals(1)
    msg = "..> Contrainte uniaxiale, coefficient de Poisson : "//TRIM(ADJUSTL(msg))
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
CASE(2060)
  msg = "..> Le système a bien été déformé."
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
  !reals(7) = pos3 (only for loops)
  !reals(8) = loop radius
  temp = TRIM(ADJUSTL(strings(1)))
  IF(TRIM(temp)=="screw") THEN
    msg = ">>> Insertion d'une dislocation vis suivant"
  ELSEIF(temp(1:4)=="edge") THEN
    msg = ">>> Insertion d'une dislocation coin suivant"
  ELSEIF(temp(1:5)=="mixed") THEN
    msg = ">>> Insertion d'une dislocation mixte suivant"
  ELSEIF(temp(1:4)=="loop") THEN
    msg = ">>> Insertion d'une boucle de dislocation dans un plan normal à"
  ENDIF
  msg = TRIM(msg)//' '//TRIM(strings(2))//","
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  !
  IF( reals(4)>0.1d0 ) THEN
    msg = "    en utilisant l'élasticité anisotrope,"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
  !
  IF(TRIM(strings(1))=="edge") THEN
    WRITE(msg,"(a34)") "    en insérant un plan d'atomes,"
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF(TRIM(strings(1))=="edge2") THEN
    WRITE(msg,"(a43)") "    en conservant le nombre total d'atomes,"
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
      msg = "    Centre ("//TRIM(ADJUSTL(temp))//","//TRIM(ADJUSTL(temp2))//","// &
          & TRIM(ADJUSTL(temp3))//") ; Rayon "//TRIM(ADJUSTL(temp4))//" A ; b="//TRIM(ADJUSTL(msg))
    ELSE
      WRITE(temp4,"(f16.3)") DABS(reals(8))
      msg = "    Centre ("//TRIM(ADJUSTL(temp))//","//TRIM(ADJUSTL(temp2))//","// &
          & TRIM(ADJUSTL(temp3))//") ; Côté "//TRIM(ADJUSTL(temp4))//" A ; b="//TRIM(ADJUSTL(msg))
    ENDIF
  ELSE
    msg = "    b="//TRIM(ADJUSTL(msg))//" à ("// &
        & TRIM(ADJUSTL(temp))//","//TRIM(ADJUSTL(temp2))//")"
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2062)
  msg = "..> Calcul des solutions aux équations anisotropes..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2063)
  !reals(1) = number of inserted atoms
  WRITE(msg,*) NINT(reals(1))
  msg = "..> "//TRIM(ADJUSTL(msg))//" atomes ont été introduits."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2064)
  msg = "..> La boîte a été agrandie."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2065)
  msg = "..> La dislocation a bien été créée."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2066)
  !reals(1) = number of repetitions along X
  !reals(2) = number of repetitions along Y
  !reals(3) = number of repetitions along Z
  WRITE(temp,*) NINT( reals(1) )
  WRITE(temp2,*) NINT( reals(2) )
  WRITE(msg,*) NINT( reals(3) )
  msg = ">>> Duplication du système : "//TRIM(ADJUSTL(temp))//" x "// &
    & TRIM(ADJUSTL(temp2))//" x "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2067)
  !reals(1) = new number of particles
  WRITE(msg,*) NINT( reals(1) )
  msg = "..> Nouveau nombre de particules : "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2068)
  !reals(1) = new number of atoms
  WRITE(temp,*) NINT( reals(1) )
  msg = "..> Le système a bien été dupliqué ("//TRIM(ADJUSTL(temp))//" atomes)."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2069)
  !strings(1) = "all" or property to be removed
  IF(strings(1)=="all") THEN
    WRITE(msg,"(a)") ">>> Suppression de toutes les propriétés auxiliaires."
  ELSE
    msg = ">>> Suppression de la propriété auxiliaire : "//TRIM(strings(1))
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2070)
  !strings(1) = "all" or property to be removed
  IF(strings(1)=="all") THEN
    msg = "..> Les propriétés auxiliaires ont bien été supprimées."
  ELSE
    msg = "..> La propriété auxiliaire a bien été supprimée."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2071)
  msg = ">>> Orientation du système..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2072)
  msg = "..> L'orientation du système a bien été changée."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2073)
  !strings(1) = file containing properties
  msg = ">>> Lecture des propriétés du système depuis le fichier "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2074)
  !strings(1) = property that was read
  IF(strings(1)=="supercell vectors" .OR. strings(1)=="conventional vectors") THEN
    msg = "Les vecteurs de boîte ont bien été lus."
  ELSEIF(strings(1)=="elastic tensor") THEN
    msg = "Le tenseur élastique a bien été lu."
  ELSEIF(strings(1)=="compliance tensor") THEN
    msg = "Le tenseur de complaisance a bien été lu."
  ELSEIF(strings(1)=="atom charges") THEN
    msg = "Les charges des atomes ont bien été lues."
  ELSEIF(strings(1)=="atom types") THEN
    msg = "Les types des atomes ont bien été lues."
  ELSEIF(strings(1)=="system orientation") THEN
    msg = "L'orientation du système a bien été lue."
  ELSE
    msg = "La propriété '"//TRIM(ADJUSTL(strings(1)))//"' a bien été lue."
  ENDIF
  msg = "..> "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2075)
  msg = "..> Fin de la lecture des propriétés."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2076)
  !msg = ">>> Sélection de tous les atomes."
  !CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2077)
  !strings(1) = side of selected region: in or out
  !strings(2) = region geometry: sphere or box
  !reals(1) = region_1(1)
  !reals(2) = region_1(2)
  !reals(3) = region_1(3)
  !reals(4) = region_2(1)
  !reals(5) = region_2(2)
  !reals(6) = region_2(3)
  temp = ""
  temp2 = ""
  IF( strings(1)=="all" ) THEN
    msg = ">>> Selection de tous les atomes."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="index" ) THEN
    IF( strings(2)=="list " ) THEN
      msg = ">>> Selection d'une liste d'atomes."
    ELSEIF( strings(2)=="range" ) THEN
      WRITE(temp,*) NINT(reals(1))
      WRITE(temp2,*) NINT(reals(2))
      msg = ">>> Selection des atomes #"//TRIM(ADJUSTL(temp))//" à "//TRIM(ADJUSTL(temp2))//"."
    ELSE
      WRITE(temp,*) NINT(reals(1))
      msg = ">>> Selection de l'atome #"//TRIM(ADJUSTL(temp))//"."
    ENDIF
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="invert" ) THEN
    msg = ">>> Inversion de la sélection..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="above" .OR. strings(1)=="below" ) THEN
    IF(strings(1)=="above") THEN
      temp2 = "au-dessus"
    ELSE
      temp2 = "en-dessous"
    ENDIF
    IF( DABS(reals(1))<1.d12 ) THEN
      WRITE(temp,"(3f16.3)") reals(1)
    ELSEIF( reals(1)<-1.d12 ) THEN
      WRITE(temp,"(a4)") "-INF"
    ELSEIF( reals(1)>1.d12 ) THEN
      WRITE(temp,"(a4)") "+INF"
    ENDIF
    msg = ">>> Sélection des atomes "//TRIM(temp2)//" de "//TRIM(ADJUSTL(temp))// &
        & " A suivant l'axe "//TRIM(strings(3))//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="in" .OR. strings(1)=="out" ) THEN
    IF(strings(1)=="in") THEN
      temp = "à l'intérieur"
    ELSE
      temp = "à l'extérieur"
    ENDIF
    IF(strings(2)=="sphere") THEN
      temp2 = "de la sphère"
    ELSEIF(strings(2)=="box" .OR. strings(2)=="cell") THEN
      temp2 = "de la boîte"
    ELSEIF(strings(2)=="cylinder") THEN
      temp2 = "du cylindre"
    ENDIF
    msg = ">>> Selection des atomes "//TRIM(temp)//" "//TRIM(temp2)//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
    IF(TRIM(strings(2))=="box") THEN
      msg = "..> Limites de la boîte : ("
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
      msg = "..> Centre : ("//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(2)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(3)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))//")"
      WRITE(temp,"(f16.3)") reals(4)
      msg = TRIM(msg)//" ; Rayon : "//TRIM(ADJUSTL(temp))//" A."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ELSEIF(TRIM(strings(2))=="cylinder") THEN
      msg = "..> Axe suivant "//TRIM(strings(3))//"."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
      WRITE(temp,"(f16.3)") reals(1)
      msg = "..> Centre : ("//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(2)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))//")"
      WRITE(temp,"(f16.3)") reals(4)
      msg = TRIM(msg)//" ; Rayon : "//TRIM(ADJUSTL(temp))//" A."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ELSEIF(TRIM(strings(2))=="torus") THEN
      WRITE(temp,"(f16.3)") reals(1)
      msg = "..> Axe suivant "//TRIM(strings(3))//". Centre : ("//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(2)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))
      WRITE(temp,"(f16.3)") reals(3)
      msg = TRIM(msg)//","//TRIM(ADJUSTL(temp))//")."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
      WRITE(temp,"(f16.3)") reals(4)
      WRITE(temp2,"(f16.3)") reals(5)
      msg = "..> Rayon principal : "//TRIM(ADJUSTL(temp))//" A ; rayon secondaire : "//TRIM(ADJUSTL(temp2))//" A."
      CALL DISPLAY_MSG(verbosity,msg,logfile)
    ENDIF
  ELSEIF( strings(1)=="prop" ) THEN
    msg = ">>> Sélection des atomes avec "//TRIM(strings(2))
    WRITE(temp,"(f16.3)") reals(1)
    IF( reals(4)>2.d0 ) THEN
      WRITE(temp2,"(f16.3)") reals(2)
      msg = TRIM(ADJUSTL(msg))//" compris entre "//TRIM(ADJUSTL(temp))//" et "//TRIM(ADJUSTL(temp2))//"."
    ELSE
      msg = TRIM(ADJUSTL(msg))//" égal à "//TRIM(ADJUSTL(temp))//"."
    ENDIF
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="random" .OR. strings(1)=="rand" ) THEN
    WRITE(temp,*) NINT(reals(1))
    msg = ">>> Sélection aléatoire de "//TRIM(ADJUSTL(temp))
    IF( NINT(reals(1))<=1 ) THEN
      msg = TRIM(msg)//" atome"
    ELSE
      msg = TRIM(msg)//" atomes"
    ENDIF
    IF( strings(2).NE."any" .AND. strings(2).NE."all" ) THEN
      msg = TRIM(msg)//" de "//TRIM(strings(2))
    ENDIF
    msg = TRIM(msg)//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="random%" ) THEN
    WRITE(temp,'(f16.3)') reals(1)*100.d0
    msg = ">>> Sélection aléatoire de "//TRIM(ADJUSTL(temp))//"% des atomes"
    IF( strings(2).NE."any" .AND. strings(2).NE."all" ) THEN
      msg = TRIM(msg)//" de "//TRIM(strings(2))
    ENDIF
    msg = TRIM(msg)//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="neigh" ) THEN
    !strings(2) = species of neighbors
    !reals(1) = number of neighbors or cutoff radius for neighbor search
    !reals(2) = species of central atom (i.e. atom whose neighbors are searched for)
    IF( strings(2)=="all" .OR. strings(2)=="any" ) THEN
      temp2 = ""
    ELSE
      temp2 = " "//TRIM(ADJUSTL(strings(2)))
    ENDIF
    !
    IF( reals(1)==0.d0 ) THEN
      msg = " des premiers"//TRIM(temp2)//" voisins"
    ELSEIF( reals(1)>0.d0 .AND. DBLE(NINT(reals(1)))-reals(1)<1.d-12 ) THEN
      WRITE(temp,*) NINT(reals(1))
      IF( NINT(reals(1))==1 ) THEN
        msg = "du premier"//TRIM(temp2)//" voisin"
      ELSE
        msg = "des "//TRIM(ADJUSTL(temp))//" premiers"//TRIM(temp2)//" voisins"
      ENDIF
    ELSE
      WRITE(temp,'(f16.3)') DABS(reals(1))
      msg = "des"//TRIM(temp2)//" voisins dans un rayon de "//TRIM(ADJUSTL(temp))//" A"
    ENDIF
    WRITE(temp,*) NINT(reals(2))
    msg = ">>> Sélection "//TRIM(ADJUSTL(msg))//" de l'atome #"//TRIM(ADJUSTL(temp))//"..."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSEIF( strings(1)=="stl" ) THEN
    msg = ">>> Sélection des atomes dans l'objet 3-D défini dans le fichier STL : "//TRIM(ADJUSTL(strings(2)))
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSE
    !Last case: strings(1) should be an atom species
    msg = ">>> Selection de tous les atomes de "//TRIM(ADJUSTL(strings(1)))//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
CASE(2078)
  !reals(1) = number of atoms that were selected
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 atome a été sélectionné."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(temp,*) NINT( reals(1) )
    msg = "..> "//TRIM(ADJUSTL(temp))//" atomes ont été sélectionnés."
  ELSE
    msg = "..> Aucun atome n'a été sélectionné."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2079)
  !reals(1) = cutoff radius for removing atoms
  WRITE(msg,"(f16.3)") reals(1)
  msg = ">>> Suppression des atomes proches de moins de "//TRIM(ADJUSTL(msg))//" A."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2080)
  !strings(1) = species of removed atom(s) (may be empty or contain anything)
  !reals(1) = number of atoms removed
  !reals(2) = number of atoms left
  species = strings(1)
  CALL ATOMNUMBER(species,tempreal)
  IF(NINT(reals(1))<=0) THEN
    msg = "..> Aucun atome n'a été supprimé"
  ELSEIF(NINT(reals(1))==1) THEN
    IF(tempreal>0.5d0) THEN
      msg = "..> 1 atome de "//TRIM(ADJUSTL(species))//" a été supprimé"
    ELSE
      msg = "..> 1 atome a été supprimé"
    ENDIF
  ELSE
    WRITE(msg,*) NINT(reals(1))
    IF(tempreal>0.5d0) THEN
      msg = "..> "//TRIM(ADJUSTL(msg))//" atomes de "//TRIM(ADJUSTL(species))//" ont été supprimés"
    ELSE
      msg = "..> "//TRIM(ADJUSTL(msg))//" atomes ont été supprimés"
    ENDIF
  ENDIF
  WRITE(temp,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(msg))//", "//TRIM(ADJUSTL(temp))//" atomes restants."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2081)
  !strings(1) = rotation axis: x, y or z
  !reals(1) = rotation angle in degrees
  WRITE(msg,"(f16.2)") reals(1)
  msg = ">>> Rotation du système de "//TRIM(ADJUSTL(msg)) &
      & //"° autour de l'axe "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2082)
  msg = "..> Le système a bien été tourné."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2083)
  !strings(1) = axis that is tilted to produce shear: x, y or z
  !strings(2) = direction of tilt: x, y or z
  !reals(1) = magnitude of tilt in Angstroms
  WRITE(temp,"(f16.3)") reals(1)
  msg = ">>> Inclinaison de l'axe "//TRIM(strings(1))//" de "// &
      & TRIM(ADJUSTL(temp))//" A suivant "//TRIM(strings(2))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2084)
  !reals(1) = shear strain in %
  msg = "..> Le système a bien été cisaillé."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,"(f16.3)") reals(1)
  msg = "..> Déformation de cisaillement appliquée = "//TRIM(ADJUSTL(msg))//" %."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2085)
  !strings(1) = shift direction: above or below
  !strings(2) = shift axis
  !reals(1) = distance above axis
  !reals(2), reals(3), reals(4) = shift vector
  IF( DABS(reals(1))<1.d12 ) THEN
    WRITE(msg,"(3f16.3)") reals(1)
  ELSEIF( reals(1)<-1.d12 ) THEN
    WRITE(msg,"(a4)") "-INF"
  ELSEIF( reals(1)>1.d12 ) THEN
    WRITE(msg,"(a4)") "+INF"
  ENDIF
  WRITE(temp,"(f16.3)") reals(2)
  WRITE(temp2,"(f16.3)") reals(3)
  WRITE(temp3,"(f16.3)") reals(4)
  IF( strings(1)(1:5)=='above' .OR. strings(1)(1:5)=='below' ) THEN
    IF(strings(1)=="above") THEN
      msg = "au-dessus de "//TRIM(ADJUSTL(msg))
    ELSE
      msg = "en-dessous de "//TRIM(ADJUSTL(msg))
    ENDIF
    msg = ">>> Translation des atomes "//TRIM(ADJUSTL(msg))//"A suivant "  &
          & //TRIM(strings(2))//" de ("//TRIM(ADJUSTL(temp))//","//        &
          & TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))//")"
  ELSEIF( strings(1)=="select" ) THEN
    msg = ">>> Translation des atomes sélectionnés de ("//TRIM(ADJUSTL(temp))//","//     &
          & TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))//")"
  ELSE
    msg = ">>> Translation de tous les atomes de ("//TRIM(ADJUSTL(temp))//","//     &
          & TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))//")"
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2086)
  !reals(1) = number of atoms that were shifted
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 atome a été déplacé."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" atomes ont été déplacés."
  ELSE
    msg = "..> Aucun atome n'a été déplacé."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2087)
  !strings(1) = property to be sorted
  !strings(2) = sort order: up, down, pack
  IF( strings(1)=="s" .OR. strings(1)=="species" ) THEN
    temp = "numéros atomiques"
    IF(strings(2)=="up") temp = TRIM(temp)//" croissants"
    IF(strings(2)=="down") temp = TRIM(temp)//" décroissants"
  ELSEIF( strings(1)=="x" .OR. strings(1)=="y" .OR. strings(1)=="z" .OR. &
        & strings(1)=="X" .OR. strings(1)=="Y" .OR. strings(1)=="Z"  ) THEN
    temp = " coordonnées "//TRIM(strings(1))
    IF(strings(2)=="up") temp = TRIM(temp)//" croissantes"
    IF(strings(2)=="down") temp = TRIM(temp)//" décroissantes"
  ELSE
    temp = TRIM(strings(1))
  ENDIF
  IF(strings(2)=="up" .OR. strings(2)=="down") THEN
    msg = ">>> Rangement des atomes par "//TRIM(ADJUSTL(temp))//"."
  ELSE
    msg = ">>> Entassement des atomes par "//TRIM(ADJUSTL(temp))//"."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2088)
  msg = "..> Les atomes ont bien été rangés."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2089)
  !strings(1) = first atomic species
  !strings(2) = second atomic species
  msg = ">>> Substitution des atomes de "//TRIM(strings(1))//&
      & " par des atomes de "//TRIM(strings(2))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2090)
  !reals(1) = number of substituted atoms
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 atome a été substitué."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" atomes ont été substitués."
  ELSE
    msg = "..> Aucun atome n'a été substitué."
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
    IF( strings(1)=="velocities" ) THEN
      temp = "vitesses"
    ELSE
      temp = "coordonnées"
    ENDIF
    msg = ">>> Conversion des "//TRIM(strings(2))//" de "//TRIM(temp)//&
        & " vers "//TRIM(temp2)//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ELSE
    WRITE(temp,'(f16.3)') reals(1)
    msg = ">>> Multiplication des valeurs de "//TRIM(strings(1))//" par un facteur "//TRIM(ADJUSTL(temp))//"."
    CALL DISPLAY_MSG(verbosity,msg,logfile)
  ENDIF
CASE(2092)
  !strings(1) = what was converted
  IF( strings(1)=="velocities" .OR. strings(1)=="coordinates" ) THEN
    IF( strings(1)=="velocities" ) THEN
      temp = "vitesses"
    ELSE
      temp = "coordonnées"
    ENDIF
    msg = "..> Les "//TRIM(ADJUSTL(temp))//" ont été converties."
  ELSE
    msg = "..> "//TRIM(ADJUSTL(strings(1)))//" ont été redimensionnées."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2093)
  msg = ">>> Collecte des atomes à replacer dans la boîte..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2094)
  !reals(1) = number of atoms wrapped
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 atome a été replacé dans la boîte."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" atomes ont été replacés dans la boîte."
  ELSE
    msg = "..> Aucun atome n'a été déplacé."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2095)
  msg = "..> Les vecteurs de la boîte ont été ré-estimés."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2096)
  msg = "..> Les solutions ont bien été trouvées."
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
  IF( strings(2)=='above' .OR. strings(2)=='below' ) THEN
    IF(strings(1)=="all") THEN
      temp = "des coordonnées"
    ELSE
      temp = "des coordonnées "//TRIM(strings(1))
    ENDIF
    IF(strings(2)=="above") THEN
      temp2 = "au dessus de"
    ELSE
      temp2 = "en dessous de"
    ENDIF
    msg = ">>> Fixation "//TRIM(ADJUSTL(temp))//" des atomes "// &
        & TRIM(ADJUSTL(temp2))//" "//TRIM(ADJUSTL(msg))//"A suivant "//TRIM(strings(3))//"."
  ELSEIF( strings(2)=='selec' ) THEN
    msg = ">>> Fixation "//TRIM(strings(1))//" des atomes sélectionnés."
  ELSE
    msg = ">>> Fixation "//TRIM(strings(1))//" de tous les atomes."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2098)
  !reals(1) = number of atoms that were fixed
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 atome a été fixé."
  ELSEIF( NINT(reals(1))>1 ) THEN
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" atomes ont été fixés."
  ELSE
    msg = "..> Aucun atome n'a été fixé."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2099)
  msg = "..> Le tenseur élastique a été tourné."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2100)
  !reals(1) = anisotropy ratio A = 2*C44 / (C11-C12)
  !reals(2) = anisotropy factor H = 2*C44 + C12 - C11
  WRITE(msg,'(e16.3)') reals(1)
  msg = "..> Ratio d'anisotropie : A = "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,'(e16.3)') reals(2)
  msg = "..> Facteur d'anisotropie : H = "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2101)
  !strings(1) = formula of energy factor, e.g. "Kb²"
  !reals(1) = energy factor
  msg = "..> Les contraintes dues à la dislocation ont été calculées."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,'(f24.8)') reals(1)
  msg = "..> Facteur d'énergie pré-logarithmique : "//TRIM(ADJUSTL(strings(1)))//&
      & " = "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2102)
  !strings(1) = empty or atom species to remove or "SEL"
  !reals(1) = 0 or atom index to remove
  IF( strings(1)=="SEL" ) THEN
    msg = ">>> Suppression de tous les atomes sélectionnés."
  ELSE
    IF( NINT(reals(1)).NE.0 ) THEN
      WRITE(msg,*) NINT(reals(1))
      msg = ">>> Suppression de l'atome #"//TRIM(ADJUSTL(msg))//"."
    ELSE
      msg = ">>> Suppression de tous les atomes de "//TRIM(strings(1))//"."
    ENDIF
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2103)
  msg = ">>> Réduction de l'inclinaison des vecteurs de boîte..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2104)
  !reals(1) = number of components that were unskewed
  IF( NINT(reals(1))==0 ) THEN
    msg = "..> Aucun vecteur de la boîte n'a été redressé."
  ELSE
    msg = "..> La boîte a été redressée."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2105)
  msg = ">>> Ré-association des coquilles avec leurs cœurs..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2106)
  !reals(1) = number of shells that were re-assigned
  IF( NINT(reals(1))==0 ) THEN
    msg = "..> Aucune coquille n'a été ré-attribuée."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 coquille a été ré-attribuée à son cœur."
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" coquilles ont été ré-attribuées à leurs cœurs."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2107)
  !strings(1) = atom species to be removed
  !reals(1) = number of atoms to be randomly removed
  WRITE(msg,*) NINT(reals(1))
  IF( strings(1)=="any" .OR. strings(1)=="all" ) THEN
    msg = ">>> Suppression aléatoire de "//TRIM(ADJUSTL(msg))//" atomes..."
  ELSE
    msg = ">>> Suppression aléatoire de "//TRIM(ADJUSTL(msg))//" atomes de " &
        & //TRIM(strings(1))//"..."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2108)
  !strings(1) = crack mode (I, II or III)
  !strings(2) = direction of crack line (X, Y or Z)
  !reals(1) = pos1
  !reals(2) = pos2
  WRITE(msg,"(f16.3)") reals(1)
  WRITE(temp,"(f16.3)") reals(2)
  msg = ">>> Création d'une fracture de mode "//TRIM(ADJUSTL(strings(1)))//  &
      & " suivant "//TRIM(ADJUSTL(strings(2)))//" à ("//TRIM(ADJUSTL(msg))// &
      & ","//TRIM(ADJUSTL(temp))//")."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2109)
  msg = "..> La fracture a été créée."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2110)
  msg = ">>> Inversion de la zone de sélection..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2111)
  !reals(1) = target temperature
  WRITE(msg,'(f24.3)') reals(1)
  msg = ">>> Attribution d'une distribution de Maxwell-Boltzmann des vitesses à "// &
      & TRIM(ADJUSTL(msg))//" K."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2112)
  msg = "..> Les vitesses des atomes ont bien été attribuées."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2113)
  !strings(1) = file name
  msg = "..> La distribution de vitesses a été écrite dans le fichier : "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2114)
  !reals(1) = max. displacement of an atom along X
  !reals(2) = max. displacement of an atom along Y
  !reals(3) = max. displacement of an atom along Z
  msg = ">>> Perturbation des positions des atomes,"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  IF( DABS(reals(1)-reals(2))<1.d-12 .AND. DABS(reals(1)-reals(3))<1.d-12 ) THEN
    !All values are equal
    WRITE(msg,'(f24.3)') reals(1)
    msg = "..> déplacement maximal : "//TRIM(ADJUSTL(msg))//" A."
  ELSE
    !Values are different along each direction
    WRITE(temp,'(f24.3)') reals(1)
    WRITE(temp2,'(f24.3)') reals(2)
    WRITE(temp3,'(f24.3)') reals(3)
    msg = "..> déplacement maximal : dx="//TRIM(ADJUSTL(temp))//", dy=" &
        & //TRIM(ADJUSTL(temp2))//", dz="//TRIM(ADJUSTL(temp3))//"."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2115)
  msg = "..> Les atomes ont été déplacés."
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
    msg = ">>> Ajout d'un atome de "//TRIM(ADJUSTL(strings(1)))//" en ("//TRIM(ADJUSTL(temp))//","//     &
        & TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))//")."
  CASE("near","NEAR")
    WRITE(temp,*) NINT(reals(1))
    msg = ">>> Ajout d'un atome de "//TRIM(ADJUSTL(strings(1)))//" près de l'atome #"//TRIM(ADJUSTL(temp))
  CASE("random","RANDOM","rand","RAND")
    WRITE(temp,*) NINT(reals(1))
    msg = ">>> Ajout de "//TRIM(ADJUSTL(temp))//" atomes de "//TRIM(ADJUSTL(strings(1)))//" à des positions aléatoires..."
  END SELECT
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2117)
  !reals(1) = number of atoms added
  !reals(2) = number of atoms in the system
  IF(NINT(reals(1))==0) THEN
    msg = "..> Aucun atome n'a été ajouté"
  ELSEIF(NINT(reals(1))==1) THEN
    msg = "..> 1 atome a été ajouté"
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" atomes ont été ajoutés"
  ENDIF
  WRITE(temp,*) NINT(reals(2))
  msg = TRIM(ADJUSTL(msg))//", "//TRIM(ADJUSTL(temp))//" atomes dans le système."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2118)
  !reals(1) = index of atom to center; if <=0 then center of mass
  IF( NINT(reals(1))<=0 ) THEN
    msg = ">>> Translation du centre de masse vers le centre de la boîte..."
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = ">>> Translation du système pour que l'atome #"//TRIM(ADJUSTL(msg))//" soit au centre de la boîte..."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2119)
  !reals(1) = X component of shift vector
  !reals(2) = Y component of shift vector
  !reals(3) = Z component of shift vector
  WRITE(temp,'(f16.3)') reals(1)
  WRITE(temp2,'(f16.3)') reals(2)
  WRITE(temp3,'(f16.3)') reals(3)
  msg = "..> Le système a été re-centré, vecteur déplacement : ("//TRIM(ADJUSTL(temp))// &
      & ","//TRIM(ADJUSTL(temp2))//","//TRIM(ADJUSTL(temp3))//")."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2120)
  !strings(1) = direction normal to mirror plane
  !reals(1) = distance between mirror plane and origin of coordinates
  WRITE(temp,'(f16.3)') reals(1)
  msg = ">>> Application d'un plan miroir à "//TRIM(ADJUSTL(temp))//" A suivant "//TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2121)
  msg = "..> Le système a bien été transformé."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2122)
  !strings(1) = species on which shells are removed
  msg = ">>> Suppression des coquilles sur les ions de "//TRIM(ADJUSTL(strings(1)))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2123)
  !reals(1) = number of shells removed
  WRITE(temp,*) NINT(reals(1))
  IF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 coquille a été supprimée."
  ELSE
    msg = "..> "//TRIM(ADJUSTL(temp))//" coquilles ont été supprimées."
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
    msg = ">>> Application d'une contrainte σ_"//TRIM(ADJUSTL(temp2))//" = "//TRIM(ADJUSTL(temp))//" GPa."
  CASE('p','P')
    msg = ">>> Application d'une pression hydrostatique de "//TRIM(ADJUSTL(temp))//" GPa."
  CASE DEFAULT
    msg = ">>> Application de l'état de contrainte défini dans "//TRIM(ADJUSTL(temp))//"."
  END SELECT
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2125)
  !strings(1) = Cartesian axis, or integer
  !strings(2) = same type as strings(1)
  SELECT CASE(strings(1))
  CASE('x','X','y','Y','z','Z')
    msg = ">>> Échange des axes cartésiens "//TRIM(ADJUSTL(strings(1)))//" et "//TRIM(ADJUSTL(strings(2)))//"."
  CASE DEFAULT
    msg = ">>> Échange des atomes #"//TRIM(ADJUSTL(strings(1)))//" et #"//TRIM(ADJUSTL(strings(2)))//"."
  END SELECT
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2126)
  msg = "..> Échange réussi."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2127)
  !strings(1) = roll axis: x, y or z
  !reals(1) = roll angle in degrees
  WRITE(msg,"(f16.2)") reals(1)
  msg = ">>> Enroulement de la direction "//TRIM(ADJUSTL(strings(1)))//" d'un angle de " // &
      & TRIM(ADJUSTL(msg))//"° autour de l'axe "//TRIM(strings(2))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2128)
  msg = "..> Le système a bien été enroulé."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2129)
  !strings(1) = torsion axis: x, y or z
  !reals(1) = torsion angle in degrees
  WRITE(msg,"(f16.2)") reals(1)
  msg = ">>> Torsion du système de "//TRIM(ADJUSTL(msg)) &
      & //"° autour de "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2130)
  msg = "..> La torsion a bien été appliquée."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2131)
  !strings(1) = space group name or number
  msg = ">>> Application des opérations de symétrie du groupe d'espace : "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2132)
  !reals(1) = space group number
  WRITE(temp,"(i5)") NINT(reals(1))
  msg = "..> Numéro du groupe d'espace : "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2133)
  !reals(1) = new number of atoms
  WRITE(temp,"(i16)") NINT(reals(1))
  msg = "..> Les opérations de symétrie ont bien été appliquées, nouveau nombre d'atomes : "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2140)
  !reals(1) = min. separation distance
  !reals(2) = separate atoms by this amount
  WRITE(temp,"(f16.2)") reals(1)
  WRITE(temp2,"(f16.2)") reals(2)
  msg = ">>> Séparation des atomes proches de moins de "//TRIM(ADJUSTL(temp))//" A d'une distance de "//TRIM(ADJUSTL(temp2))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2141)
  !reals(1) = number of pairs of atoms that were separated
  IF( reals(1) < 1.d-3 ) THEN
    msg = "..> Aucun atome n'a été séparé."
  ELSEIF( NINT(reals(1)) == 1 ) THEN
    msg = "..> Une paire d'atomes a été séparée."
  ELSE
    WRITE(temp,"(i16)") NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(temp))//" paires d'atomes ont été séparées."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2142)
  !reals(1) = number of triangles in STL file
  WRITE(temp,*) NINT(reals(1))
  msg = "..> Le fichier STL a bien été lu ("//TRIM(ADJUSTL(temp))//" triangles)."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2143)
  msg = ">>> Conversion du système en une boîte orthorhombique..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(2144)
  WRITE(temp,*) NINT(reals(1))
  msg = "..> La boîte est désormais orthorhombique ("//TRIM(ADJUSTL(temp))//" atomes)."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!2700-2799: WARNING MESSAGES
CASE(2700)
  !strings(1) = option name
  msg = "/!\ ALERTE : impossible de comprendre cette option : "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Essayez 'atomsk --help options' pour un résumé des options."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2720)
  msg = "/!\ ALERTE : l'axe est déjà aligné, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2721)
  msg = "/!\ ALERTE : les coordonnées sont déjà réduites, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2722)
  !strings(1) = atom species
  msg = "/!\ ALERTE : les atomes de "//TRIM(strings(1))//" ont déjà des coquilles."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2723)
  !strings(1) = atomic species
  msg = "/!\ ALERTE : il n'y a pas d'atome de "//TRIM(strings(1))//" dans le système, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2724)
  msg = "/!\ ALERTE : la déformation est nulle, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2725)
  msg = "/!\ ALERTE : le vecteur de Burgers est nul, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2726)
  msg = "/!\ ALERTE : la boîte est très petite dans une dimension normale à la ligne"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    de dislocation ! Êtes-vous sûr de savoir ce que vous faites ?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2727)
  !reals(1) = index of atom with large displacement
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ ALERTE : le déplacement de l'atome "//TRIM(ADJUSTL(msg))// " est grand."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2728)
  msg = "/!\ ALERTE : tous les facteurs sont égaux à 1, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2729)
  msg = "/!\ ALERTE : aucune propriété auxiliaire n'est définie, abandon"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2730)
  !string(1) = property
  msg = "/!\ ALERTE: pas de "//TRIM(ADJUSTL(strings(1)))//&
      & " dans les propriétés auxiliaires, abandon..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2731)
  msg = "/!\ ALERTE : Hstart = Hend, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2732)
  !strings(1) = name of unknown property
  msg = "/!\ ALERTE : propriété inconnue : "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2733)
  msg = "/!\ ALERTE : le rayon spécifié est négatif, aucun atome ne sera supprimé."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2734)
  msg = "/!\ ALERTE : l'angle de rotation est modulo zéro, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2735)
  msg = "/!\ ALERTE : cisaillement nul, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2736)
  msg = "/!\ ALERTE : le vecteur de translation est nul, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2737)
  msg = "/!\ ALERTE : les espèces sont les mêmes, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2738)
  msg = "/!\ ALERTE : les unités sont les mêmes, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2739)
  msg = "/!\ ALERTE : les vecteurs de base ne sont pas orthonormaux."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2740)
  msg = "/!\ ALERTE : le tenseur élastique n'est pas symétrique !"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2741)
  msg = "/!\ ALERTE : le coefficient de Poisson est en dehors des limites [-1 , 0.5]."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2742)
  !reals(1) = atom index
  WRITE(temp,*) NINT(reals(1))
  msg = "/!\ ALERTE : l'indice donné #"//TRIM(ADJUSTL(temp))//" est en dehors des limites."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2743)
  msg = "/!\ ALERTE : les vecteurs de boîte ne sont pas inclinés, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2744)
  msg = "/!\ ALERTE : aucune coquille dans le système, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2745)
  !reals(1) = number of atoms that will actually be removed
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ ALERTE : impossible de sélectionner autant d'atomes, seulement "// &
      & TRIM(ADJUSTL(msg))//" seront sélectionnés."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2746)
  msg = "/!\ ALERTE : il n'y a rien à faire, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2747)
  msg = "/!\ ALERTE : le facteur K est nul, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2748)
  msg = "/!\ ALERTE : la boîte est très petite dans une dimension normale à la ligne"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    de fracture ! Êtes-vous sûr de savoir ce que vous faites ?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2749)
  !strings(1) = name of property
  !reals(1) = atom index
  WRITE(temp,*) NINT(reals(1))
  msg = "/!\ ALERTE : impossible d'assigner la valeur de la propriété "//TRIM(ADJUSTL(strings(1)))// &
      & " à l'atome #"//TRIM(ADJUSTL(temp))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2750)
  msg = "/!\ ALERTE : une sélection était définie mais elle ne contient plus aucun atome."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "             La sélection a été annulée, tous les atomes sont désormais sélectionnés."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2751)
  msg = "/!\ ALERTE : la temperature cible est nulle, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2752)
  msg = "/!\ ALERTE : aucune sélection n'est définie, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2753)
  !reals(1) = atom index
  WRITE(temp,*) NINT(reals(1))
  msg = "/!\ ALERTE : l'atome #"//TRIM(ADJUSTL(temp))//" est déjà sélectionné."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2754)
  !strings(1) = type of object (e.g. "dislocation" or "crack" or "loop" or "all"
  !             ("all" means that all points of a loop are out of the box)
  !reals(1) = number of points that are out of the box
  IF( strings(1)=="loop" ) THEN
    WRITE(temp,*) NINT(reals(1))
    msg = "/!\ ALERTE : "//TRIM(ADJUSTL(temp))//" points de la boucle sont en dehors de la boîte."
    CALL DISPLAY_MSG(1,msg,logfile)
  ELSEIF( strings(1)=="all" ) THEN
    msg = "/!\ ALERTE : tous les points de la boucle sont en dehors de la boîte !"
    CALL DISPLAY_MSG(1,msg,logfile)
    msg = "    Êtes-vous sûr de savoir ce que vous faites ?"
    CALL DISPLAY_MSG(1,msg,logfile)
  ELSE
    msg = "/!\ ALERTE : la position de la "//TRIM(ADJUSTL(strings(1)))//" est à l'extérieur de la boîte."
    CALL DISPLAY_MSG(1,msg,logfile)
    msg = "    Êtes-vous sûr de savoir ce que vous faites ?"
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(2755)
  msg = "/!\ ALERTE : le facteur est égal à zéro, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2756)
  !strings(1) = direction
  msg = "/!\ ALERTE : la boîte semble très large suivant la direction "//TRIM(ADJUSTL(strings(1)))//" !"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Êtes-vous sûr de savoir ce que vous faites ?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2757)
  msg = "/!\ ALERTE : les indices sont identiques, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2758)
  msg = "/!\ ALERTE : aucune opération à appliquer, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2759)
  msg = "/!\ ALERTE : le rayon de la boucle de dislocation est trop petit, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2760)
  msg = "/!\ ALERTE : la boîte est déjà orthorhombique, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
  !
CASE(2799)
  !strings(1) = name of obsolete option
  !strings(2) = name of new option
  msg = "/!\ ALERTE : l'option "//TRIM(ADJUSTL(strings(1)))//" est obsolète et sera supprimée."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Veuillez employer l'option "//TRIM(ADJUSTL(strings(2)))//" à la place."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!2800-2899: ERROR MESSAGES
CASE(2800)
  !string(1) = axis (if we are here it"s because it is different from x, y or z)
  msg = "X!X ERREUR : axe inconnu: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2801)
  WRITE(msg,*) "X!X ERREUR : la base Hend n'est pas une rotation de Hstart."
  CALL DISPLAY_MSG(1,msg,logfile)
  WRITE(msg,*) "    Vérifiez que les angles sont égaux dans Hstart et Hend."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2802)
  !strings(1) = property that was not read properly
  msg = "X!X ERREUR en lisant "//TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2803)
  msg = "X!X ERREUR : il y a eu une erreur en tentant de recalculer les vecteurs de boîte."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2804)
  msg = "X!X ERREUR : il y a eu des erreurs en appliquant des options."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2805)
  !strings(1) = name of unknown option
  msg = "X!X ERREUR: option inconnue : "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2806)
  !strings(1) = name of option
  msg = "X!X ERREUR : syntaxe incorrecte dans l'option : "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2807)
  !reals(1) = 1 if roots of Eq.(13-85) cannot be found
  !         = 2 if the A_k(n) cannot be calculated
  !         = 3 if the linear equations defining D(n) cannot be solved
  msg = "X!X ERREUR :"
  IF(NINT(reals(1))==1) THEN
    msg = TRIM(msg)//" impossible de déterminer les P(n), abandon."
  ELSEIF(NINT(reals(1))==2) THEN
    msg = TRIM(msg)//" impossible de déterminer les A_k(n), abandon."
  ELSEIF(NINT(reals(1))==3) THEN
    msg = TRIM(msg)//" impossible de déterminer les D(n), abandon."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2808)
  msg = "X!X ERREUR : impossible de construire une dislocation mixte en élasticité isotrope."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2809)
  msg = "X!X ERREUR : le tenseur élastique contient des valeurs NaN, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2810)
  msg = "X!X ERREUR : taille incohérente du tableau pour les coquilles."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2811)
  msg = "X!X ERREUR : dislocline et dislocplane doivent être perpendiculaires, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2812)
  msg = "X!X ERREUR : impossible de déterminer quel(s) atome(s) supprimer, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2813)
  !strings(1) = string that could not be converted to a number
  msg = "X!X ERREUR : impossible de convertir ces caractères en nombre : "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2814)
  msg = "X!X ERREUR : il n'existe aucun système auquel appliquer cette option."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2815)
  !strings(1) = name of matrix to invert
  msg = "X!X ERREUR : impossible d'inverser la matrice "//strings(1)//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2816)
  msg = "X!X ERREUR : le tenseur élastique n'est pas défini, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2817)
  msg = "X!X ERREUR : la propriété '"//TRIM(ADJUSTL(strings(1)))//"' n'est pas définie, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2818)
  msg = "X!X ERREUR : une erreur s'est produite lors de la lecture du fichier STL, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(2819)
  msg = "X!X ERREUR : impossible de trouver une boîte orthorhombique à partir des vecteurs de boîte initiaux."
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
    msg = ">>> Écriture des fichiers ("//TRIM(ADJUSTL(temp))//" cœurs + "// &
        & TRIM(ADJUSTL(temp2))//" coquilles):"
  ELSE
    msg = ">>> Écriture des fichiers ("//TRIM(ADJUSTL(temp))//" atomes):"
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3001)
  !strings(1) = name of file
  msg = "..> Ce fichier existait déjà et n'a pas été reconverti : "// &
      & TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3002)
  !strings(1) = type of file, e.g. "XSF" or "CFG"
  !strings(2) = name of file
  msg = "..> Le fichier "//TRIM(strings(1))//" a bien été écrit : " &
      & //TRIM(strings(2))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3003)
  !reals(:) = list of atomic numbers
  msg = "..> Les espèces atomiques ont été tassées :"
  DO i=1,SIZE(reals(:))
    CALL ATOMSPECIES(reals(i),species)
    msg = TRIM(msg)//" "//TRIM(species)//", "
  ENDDO
  j=LEN_TRIM(msg)
  msg = TRIM(msg(:j-1))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "    Vérifiez que cet ordre est cohérent avec le fichier POTCAR."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3004)
  !strings(1) = skew parameter
  msg = "..> L'inclinaison "//TRIM(strings(1))//" a été réduite."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3005)
  msg = "..> OK, l'inclinaison "//TRIM(strings(1))//" ne sera pas changée."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(3006)
  !strings(1) = name of ddplot file
  msg = ">>> Écriture du fichier ddplot : "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!3700-3799: WARNING MESSAGES
CASE(3700)
  msg = "/!\ ALERTE : aucun nom de fichier de sortie n'a été spécifié, veuillez en fournir un :"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3701)
  !strings(1) = name of output file
  msg = "/!\ ALERTE : je n'ai pas pu déterminer le format de ce fichier de sortie: " &
      & //TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "<?> Veuillez indiquer un format pour ce fichier :"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = '    ('//TRIM(listofformats(1))
  DO i=2,SIZE(listofformats)
    msg = TRIM(msg)//','//TRIM(listofformats(i))
  ENDDO
  msg = TRIM(msg)//')'
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3702)
  msg = "/!\ ALERTE : l'écriture au format ddplot n'est disponible qu'en mode DDplot."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Aucun fichier ddplot ne sera écrit, référez-vous à la documentation."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3703)
  msg = "/!\ ALERTE : les espèces atomiques ne sont pas contigus. Voulez-vous les tasser ? ("&
      & //langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    (ceci n'affectera que le fichier POSCAR)"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3704)
  msg = "/!\ ALERTE : les vecteurs de la boîte ne forment pas une matrice triangle, ce qui"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    est requis par LAMMPS. Voulez-vous ré-aligner le système ? (" &
      & //langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    (ceci n'affectera que le fichier LAMMPS)"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3705)
  !strings(1) = skew parameter
  msg = "/!\ ALERTE: l'inclinaison "//TRIM(strings(1))//" est trop grande."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Exécuter LAMMPS avec de tels vecteurs de boîte engendrera le même"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    message d'erreur. Voulez-vous réduire l'inclinaison ? ("&
      & //langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    (ceci n'affectera que le fichier LAMMPS)"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3706)
  !strings(1) = skew parameter
  msg = "/!\ ALERTE : impossible de réduire l'inclinaison "//TRIM(strings(1))//", abandon..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3707)
  !reals(1) = number of atoms
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ ALERTE : le nombre de particules est très grand : "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Si ddplot ne peut pas afficher autant d'atomes, vous pouvez utiliser Atomsk"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    avec l'option -cut pour réduire le nombre d'atomes."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3708)
  msg = "/!\ ALERTE : seules les 32 premières propriétés auxiliaires seront "// &
      & "écrites dans le fichier CFG."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3709)
  msg = "/!\ ALERTE: certaines composantes des vecteurs de boîte ne peuvent"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    pas être écrits et ce fichier SERA INCORRECT !!!"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3710)
  !strings(1) = file format
  msg = "/!\ ALERTE : format '"//TRIM(strings(1))//"' inconnu, abandon..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3711) ! missing occupancy data
  msg = "/!\ ALERTE : l'occupation des sites est manquante, l'occupation sera fixée à 1 pour tous les atomes."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3712) ! missing thermal vibration data
  msg = "/!\ ALERTE : les facteurs de Debye-Waller sont manquants,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            la vibration thermique sera fixée à zéro pour tous les atomes."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3713) ! missing absorption data
  msg = "/!\ ALERTE : les facteurs d'absorption sont manquants, ils seront fixés à 0.03 pour tous les atomes."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3714)
  msg = "/!\ ALERTE : certains atomes ont un type non valide."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3715)
  msg = "/!\ ALERTE : les données d'entrée contiennent des occupations partielles,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            que certains formats de sortie ne supportent pas."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Certains atomes risquent de se chevaucher dans les fichiers de sortie."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3716)
  msg = "/!\ ALERTE : les données contiennent les positions de coquilles (shells),"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            que certains formats de sortie ne supportent pas."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Les coquilles seront perdues dans certains fichiers de sortie."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3717)
  !reals(1) = total charge
  WRITE(temp,'(f9.3)') reals(1)
  msg = "/!\ ALERTE : la boîte a une charge électrique non nulle : Q_total = "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(1,msg,logfile)
!
!3800-3899: ERROR MESSAGES
CASE(3800)
  msg = "X!X ERREUR : aucun atome dans le système, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3801)
  !strings(1) = name of file
  msg = "X!X ERREUR : il y a eu des erreurs en tentant d'écrire le fichier: "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3802)
  msg = "X!X ERREUR : il y a eu des erreurs en tentant d'écrire les fichiers."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3803)
  !reals(1) = index of atom that has a NaN coordinate
  WRITE(msg,*) NINT(reals(1))
  msg = "X!X ERREUR : l'atome #"//TRIM(ADJUSTL(msg))//" a une coordonnée NaN, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(3804)
  !reals(1) = index of atom that has a NaN coordinate
  WRITE(msg,*) NINT(reals(1))
  msg = "X!X ERREUR : une propriété auxiliaire de l'atome #"// &
      & TRIM(ADJUSTL(msg))//" est NaN, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!
!
!**************************
! 4000-4999: MESSAGES FOR MODES
CASE(4000)
  !strings(1) = name of file for mode "list"
  msg = ">>> Lecture de la liste de fichiers depuis : "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4001)
  !Just a separator
  msg = "         --===============--"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4002)
  !strings(1) = output format activated, e.g. "xyz" or "cfg"
  msg = ">>> Conversion de tous les fichiers suivants au format "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4010)
  msg = ">>> Utilisation du mode DDplot."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4020)
  msg = ">>> Décomposition de "//TRIM(strings(1))//" en plusieurs fichiers..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4021)
  msg = ">>> Atomsk est un logiciel libre et Open Source."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "    Pour plus d'informations, entrez 'license'."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4022)
  msg = ">>> Interpréteur de ligne de commande de Atomsk :"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = '..> Tapez "help" pour connaître les commandes disponibles.'
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = ""
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4023)
  msg = "Atomsk s'exécute en ce moment en MODE INTERACTIF."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "Seules les commandes suivantes sont disponibles :"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "help               Affiche cette aide"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "cd                 Changer de dossier"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = system_ls//"                Affiche la liste des fichiers du dossier courant"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "pwd                Affiche le dossier courant"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "memory             Résumé du contenu de la mémoire"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "create             Créer un système atomique"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "read <fichier>     Lit le <fichier> et charge son contenu en mémoire"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "write <fichier>    Ecrit le système courant dans le <fichier>"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "clear              Efface la mémoire (détruit le système atomique)"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "quit               Quitte Atomsk"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "OPTIONS : les options de Atomsk peuvent être utilisées,"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "         entrez 'help options' pour afficher les options disponibles."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "         En mode interactif les options doivent être appelées sans le signe moins (-)."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "MODES : les modes ne peuvent pas être utilisés dans cet interpréteur."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4024)
  msg = "<?> Vers quel format souhaitez-vous le convertir ?"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = '    ('//TRIM(listofformats(1))
  DO i=2,SIZE(listofformats)
    msg = TRIM(msg)//','//TRIM(listofformats(i))
  ENDDO
  msg = TRIM(msg)//')'
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4025)
  msg = ">>> Oups, je n'ai pas pu convertir votre fichier, désolé..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "    J'espère que nous sommes toujours amis ! :-)"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4026)
  msg = ">>> Le fichier a bien été converti."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4027)
  msg = ">>> Création du système :"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4028)
  !strings(1) = description of system that is created
  msg = "..> "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4029)
  msg = "..> Le système a bien été créé."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4030)
  !strings(1) = merge direction
  !reals(1) >0 if systems are concatenated, <0 otherwise
  IF( reals(1)>0.d0 ) THEN
    msg = ">>> Collage des systèmes suivant "//TRIM(strings(1))//"..."
  ELSE
    msg = ">>> Fusion des systèmes..."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4031)
  !reals(1) = number of files merged
  WRITE(msg,*) NINT(reals(1))
  msg = "..> "//TRIM(ADJUSTL(msg))//" systèmes ont été fusionnés."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4032)
  msg = ">>> Calcul des moments électriques :"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4033)
  !strings(1) = name of file where total polarization is written
  msg = "..> La polarisation totale a été écrite : "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4034)
  !strings(1) = first atomic species
  !strings(2) = second atomic species
  !reals(1) = coordination
  IF(INT(reals(1))==4) THEN
    msg = " tetraèdres"
  ELSEIF(INT(reals(1))==6) THEN
    msg = " octaèdres"
  ELSEIF(INT(reals(1))==8) THEN
    msg = " hexaèdres"
  ELSE
    msg = " polyèdres"
  ENDIF
  msg = ">>> Calcul des moments des"//TRIM(msg)//" "//TRIM(strings(1))// &
      & "-"//TRIM(strings(2))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  IF(reals(1)<0.d0) THEN
    WRITE(msg,"(f16.2)") DABS(reals(1))
    msg = "..> En utilisant les voisins dans un rayon de "//TRIM(ADJUSTL(msg))// &
        & " A autour des ions "//TRIM(strings(1))//"."
  ELSEIF(reals(1)==0.d0) THEN
    msg = "..> En essayant de trouver les premiers voisins automatiquement."
  ELSE
    WRITE(msg,*) INT(reals(1))
    msg = "..> En utilisant les "//TRIM(ADJUSTL(msg))//" premiers voisins des ions "// &
        & TRIM(strings(1))//"."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4035)
  !strings(1) = file name
  msg = "..> Le fichier XSF a bien été écrit : "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4036)
  !strings(1) = file name
  msg = "..> Les normes ont bien été écrites : "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4037)
  !strings(1) = file name
  msg = "..> Les statistiques ont bien été écrites : "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4038)
  msg = ">>> Calcul des différences..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4039)
  !strings(1) = name of file
  msg = "..> Ce fichier a bien été écrit : "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4040)
  msg = "..> Les déplacements ont bien été calculés."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4041)
  !reals(1) = snapshots number
  WRITE(msg,*) NINT(reals(1))
  msg = ">>> Lecture de l'image #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4042)
  !reals(1) = number of snapshots read
  WRITE(msg,*) NINT(reals(1))
  IF( NINT(reals(1)) <= 0 ) THEN
    msg = ">>> Aucune image n'a été convertie."
  ELSEIF( NINT(reals(1)) == 1 ) THEN
    msg = ">>> 1 image a bien été convertie."
  ELSE
    msg = ">>> "//TRIM(ADJUSTL(msg))//" images ont bien été converties."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4043)
  msg = "..> L'image a bien été lue."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4044)
  msg = ">>> Écriture des propriétés pour tous les atomes :"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4045)
  !reals(1) = number of files converted
  WRITE(msg,*) NINT(reals(1))
  IF( NINT(reals(1))==0 ) THEN
    msg = ">>> Aucun fichier n'a été converti"
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = ">>> 1 fichier a été converti"
  ELSE
    msg = ">>> "//TRIM(ADJUSTL(msg))//" fichiers ont été convertis"
  ENDIF
  IF( SIZE(reals)>1 .AND. NINT(reals(2))>0 ) THEN
    WRITE(temp,*) NINT(reals(2))
    IF( NINT(reals(2))==1 ) THEN
      msg = TRIM(ADJUSTL(msg))//", 1 fichier a été ignoré"
    ELSE
      msg = TRIM(ADJUSTL(msg))//", "//TRIM(ADJUSTL(temp))//" fichiers ont été ignorés"
    ENDIF
  ENDIF
  msg = TRIM(ADJUSTL(msg))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4046)
  msg = ">>> Déballage des atomes..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4047)
  !reals(1) = number of atoms unwrapped
  IF( NINT(reals(1))==0 ) THEN
    msg = "..> Aucun atome n'a été déballé."
  ELSEIF( NINT(reals(1))==1 ) THEN
    msg = "..> 1 atome a été déballé."
  ELSE
    WRITE(msg,*) NINT(reals(1))
    msg = "..> "//TRIM(ADJUSTL(msg))//" atomes ont été déballés."
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4048)
  msg = ">>> Calcul de la polarisation électronique des ions..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4049)
  msg = "..> La polarisation électronique a été calculée."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4050)
  !strings(1) = file containing the names of files to include
  msg = ">>> Écriture de tous les fichiers listés dans "//TRIM(strings(1))// &
      & " vers un seul fichier..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4051)
  !reals(1) = width of the skin (A)
  msg = ">>> Calcul de la fonction de distribution radiale,"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,'(f16.3)') reals(1)
  msg = "    en utilisant un pas de "//TRIM(ADJUSTL(msg))//" A."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4052)
  !strings(1) = species 1
  !strings(2) = species 2
  msg = "..> Calcul de la FDR des atomes de "//TRIM(strings(2))// &
      & " autour des atomes de "//TRIM(strings(1))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4053)
  msg = "..> Les FDR ont été calculées."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4054)
  msg = ">>> Construction d'un polycrystal avec la méthode de Voronoï."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4055)
  !reals(1) = index of the grain
  WRITE(msg,*) NINT(reals(1))
  msg = ">>> Construction du grain #"//TRIM(ADJUSTL(msg))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4056)
  !reals(1) = number of atoms in the grain
  WRITE(msg,*) NINT(reals(1))
  msg = "..> Nombre d'atomes dans ce grain: "//TRIM(ADJUSTL(msg))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4057)
  !strings(1) = name of input file
  msg = ">>> Lecture des paramètres pour la construction de Voronoï depuis: "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4058)
  !reals(1) = number of grains
  !reals(2) = 0 if 3-D, 1,2,3 if thin along x, y, z
  msg = "..> Le fichier a bien été lu."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(msg,*) NINT(reals(1))
  msg = "..> Nombre de grains à construire : "//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  IF( NINT(reals(2))==0 ) THEN
    msg = "..> Tesselation de Voronoi en 3-D."
  ELSE
    IF( NINT(reals(2))==1 ) THEN
      msg = "x"
    ELSEIF( NINT(reals(2))==2 ) THEN
      msg = "y"
    ELSEIF( NINT(reals(2))==3 ) THEN
      msg = "z"
    ENDIF
    msg = "..> Tesselation de Voronoi en 2-D, axe de rotation : "//TRIM(ADJUSTL(msg))
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4059)
  msg = ">>> Construction d'une chaîne de configurations par interpolation."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4060)
  !reals(1) = index of the image
  WRITE(msg,*) NINT(reals(1))
  msg = ">>> Construction de la configuration #"//TRIM(ADJUSTL(msg))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4061)
  msg = ">>> Calcul du tenseur de Nye."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4062)
  msg = ">>> Calcul de la matrice G sur chaque atome..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4063)
  msg = ">>> Calcul du tenseur de Nye sur chaque atome..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4064)
  !strings(1) = name of the file containing file list
  msg = ">>> Moyenne des positions atomiques pour la liste de fichiers : "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4065)
  !reals(1) = number of files
  WRITE(temp,*) NINT(reals(1))
  msg = ">>> Les positions atomiques ont été moyennées sur "//TRIM(ADJUSTL(temp))//" configurations."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4066)
  msg = ">>> Exécution du mode densité."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4067)
  !strings(1) = property
  !reals(1) = dimension (1, 2 or 3)
  WRITE(temp,*) NINT(reals(1))
  msg = ">>> Calcul de la densité "//TRIM(ADJUSTL(temp))//"-D de "//TRIM(ADJUSTL(strings(1)))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4068)
  WRITE(msg,*) NINT(reals(1))
  msg = "..> La densité a bien été calculée."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4069)
  !strings(1) = name of file
  msg = ">>> Calcul du paramètre de symétrie centrale pour : "//TRIM(ADJUSTL(strings(1)))//"..."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4200)
  WRITE(*,*) " (tapez q pour quitter)"
  WRITE(*,'(a45)',ADVANCE='NO') " Réseau cristallin (sc,bcc,fcc,dia,rs,per) : "
CASE(4201)
  !strings(1) = name of lattice parameter (a, b, c, a0...)
  WRITE(*,'(a30)',ADVANCE='NO') " Paramètre de maille "//TRIM(ADJUSTL(strings(1)))//" (Å) : "
CASE(4202)
  WRITE(*,'(a21)',ADVANCE='NO') " Espèces chimiques : "
CASE(4203)
  WRITE(*,'(a24)',ADVANCE='NO') " Indices chiraux (n,m) : "
CASE(4204)
  WRITE(*,'(a33)',ADVANCE='NO') " Orientation cristallographique : "
CASE(4300)
  !reals(1) = number of max tries
  msg = "   ***   ATOMSK QUIZZ   ***"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  WRITE(temp,*) NINT(reals(1))
  msg = " Répondez à la question suivante, vous avez "//TRIM(ADJUSTL(temp))//" essais."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4301)
  !strings(1) = answer
  !reals(1) = type of quizz (1=guess atomic number; 2=guess atom species; 3=guess next atom)
  !reals(2) = number of tries
  !reals(3) = user's answer - actual solution (quizz type 1 only)
  msg=""
  IF( NINT(reals(2))>0 ) THEN
    msg = " FAUX !"
  ENDIF
  IF( NINT(reals(1))==1 ) THEN
    IF( NINT(reals(2))>0 ) THEN
      IF( reals(3) > 0.d0 ) THEN
        msg = TRIM(msg)//" C'est moins."
      ELSEIF( reals(3) < 0.d0 ) THEN
        msg = TRIM(msg)//" C'est plus."
      ENDIF
    ENDIF
    msg = TRIM(msg)//" Quel est le numéro atomique de l'atome : "//TRIM(ADJUSTL(strings(1)))//" ?"
  ELSEIF( NINT(reals(1))==2 ) THEN
    msg = TRIM(msg)//" Quel atome a pour numéro atomique "//TRIM(ADJUSTL(strings(1)))//" ?"
  ELSEIF( NINT(reals(1))==3 ) THEN
    msg = TRIM(msg)//" Quel atome vient après "//TRIM(ADJUSTL(strings(1)))//" dans le tableau périodique ?"
  ENDIF
  CALL DISPLAY_MSG(verbosity,msg,logfile)
CASE(4302)
  !reals(1) = try/maxtry
  IF( reals(1)>=1.d0 ) THEN
    msg = " VOUS AVEZ PERDU !"
  ELSE
    msg = " CORRECT !"
  ENDIF
  msg = TRIM(msg)//" La réponse était "//TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  
!
!4700-4799: WARNING messages
CASE(4700)
  !strings(1) = name of file that does not exist
  msg = "/!\ ALERTE : ce fichier n'existe pas : "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Abandon..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4701)
  msg = "/!\ Je n'ai pas trouvé les vecteurs de la boîte dans le fichier d'entrée."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Je vais tenter de les deviner mais cela peut être imprécis."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Je vous recommande d'utiliser l'option -prop pour définir ces vecteurs."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4702)
  msg = "/!\ ALERTE : le système est chargé, la polarisation totale peut être fausse !"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4703)
  msg = "/!\ ALERTE : cet élément a une charge nulle, abandon..."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4705)
  !reals(1) = index of atom that has too many neighbours
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ ALERTE : le nombre de voisins est supérieur à 100 pour l'atome #" &
      & //TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4706)
  !strings(1) = atomic species of absent neighbours
  !reals(1) = index of atom that has no neighbour
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ ALERTE : je n'ai pas trouvé d'atome "//TRIM(strings(1))// &
      &       " voisin de l'atome #"//TRIM(ADJUSTL(msg))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4707)
  !reals(1) = index of atom that has a shell with zero charge
  WRITE(msg,*) NINT(reals(1))
  msg = "/!\ ALERTE : la coquille de l'ion #"//TRIM(ADJUSTL(msg)) &
      & //" a une charge nulle."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4708)
  msg = "/!\ ALERTE : ce grain ne contient aucun atome."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4709)
  !reals(1) = index of atom
  !reals(2) = number of neighbors
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = "/!\ ALERTE : nombre de voisins ("//TRIM(ADJUSTL(temp2))// &
      &              ") insuffisant pour l'atome #"//TRIM(ADJUSTL(temp))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4710)
  !strings(1) = name of the matrix
  msg = "/!\ ALERTE : impossible de calculer la matrice "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4711)
  !strings(1) = name of the file
  msg = "/!\ ALERTE : ce fichier contient un nombre différent d'atomes et ne sera pas traité : "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4712)
  msg = "/!\ ALERTE : il est recommandé d'exécuter ce mode avec l'option '-wrap',"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            sans quoi les atomes qui sont en dehors de la boîte pourraient fausser les résultats."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Voulez-vous les remettre dans la boîte maintenant ? ("//langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4713)
  !strings(1) = action to follow if user says "yes"
  temp = strings(1)
  IF(LEN_TRIM(strings(1))<=0) THEN
    temp = "continuer"
  ENDIF
  IF(strings(1)=="quit") THEN
    temp = "quitter"
  ELSEIF(strings(1)=="continue") THEN
    temp = "continuer"
  ELSEIF(strings(1)=="erase it") THEN
    temp = "l'effacer"
  ENDIF
  msg = "/!\ ALERTE : apparemment vous n'avez pas enregistré votre système dans un fichier."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            Êtes-vous sûr de vouloir "//TRIM(ADJUSTL(temp))//" ? ("//langyes//"/"//langno//")"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4714)
  !strings(1) = direction along which cell is small
  !reals(1) = new cell size along that direction
  WRITE(temp,'(f16.3)') reals(1)
  msg = "/!\ ALERTE : la boîte finale a une petite dimension suivant "&
      & //TRIM(ADJUSTL(strings(1)))//", ajustement à "//TRIM(ADJUSTL(temp))
  CALL DISPLAY_MSG(1,msg,logfile)
!
!4800-4899: ERROR MESSAGES
CASE(4800)
  !strings(1) = mode
  msg = "X!X ERREUR : erreur de syntaxe dans le mode : "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4801)
  !strings(1) = file format (e.g. "xyz", "cfg"
  msg = "X!X ERREUR : ce format de fichier n'est pas encore supporté dans le mode un-en-tout :" &
      & //TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4802)
  msg = "X!X ERREUR : les deux indices chiraux ne peuvent être nuls, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4803)
  !reals(1) = theoretical number of atoms in nanotube
  !reals(2) = number found by Atomsk
  WRITE(msg,*) "X!X ERREUR : incohérence dans le nombre d'atomes du nanotube."
  CALL DISPLAY_MSG(1,msg,logfile)
  WRITE(msg,*) NINT(reals(1))
  WRITE(temp,*) NINT(reals(2))
  msg = "          Théorie : "//TRIM(ADJUSTL(msg))//"; Trouvés : "//TRIM(ADJUSTL(temp))
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
      temp = TRIM(ADJUSTL(temp))//", "//TRIM(ADJUSTL(temp2))//" ou "//TRIM(ADJUSTL(temp3))
    ELSE
      temp = TRIM(ADJUSTL(temp))//" ou "//TRIM(ADJUSTL(temp2))
    ENDIF
  ENDIF
  msg = "X!X ERREUR : cette structure requiert "//TRIM(ADJUSTL(temp))//" espèces chimiques."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4805)
  !strings(1) = structure type
  msg = "X!X ERREUR : structure inconnue : "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4806)
  msg = "X!X ERREUR : aucun fichier à joindre, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4807)
  msg = "X!X ERREUR : Toutes les charges sont nulles, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Les charges électriques doivent être spécifiées avec l'option -properties."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4808)
  !strings(1) = atomic species that has zero charge
  msg = "X!X ERREUR: les ions "//TRIM(strings(1))//" ne peuvent pas avoir une charge nulle."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4809)
  !strings(1) = atomic species
  msg = "X!X ERREUR : aucun ion "//TRIM(strings(1))//" dans le système."
  msg = TRIM(ADJUSTL(msg))
CASE(4810)
  !reals(1) = number of particles in first system
  !reals(2) = number of particles in second system
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = "X!X ERREUR : le nombre de particules n'est pas le même dans les deux systèmes: "// &
      & TRIM(ADJUSTL(temp))//"/"//TRIM(ADJUSTL(temp2))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4811)
  msg = "X!X ERREUR : les deux systèmes n'ont pas les mêmes vecteurs de base."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4812)
  msg = "X!X ERREUR : le fichier ne semble pas être au format XSF animé, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4813)
  !strings(1) = name of unknown mode
  msg = "X!X ERREUR : mode inconnu : "//TRIM(strings(1))
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4814)
  msg = "X!X ERREUR : un seul mode peut être utilisé à la fois, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4815)
  msg = "X!X ERREUR : le fichier ne semble pas être au format DL_POLY HISTORY, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4816)
  !reals(1) = index for core charges
  !reals(2) = index for shell charges
  IF( reals(1)<reals(2) ) THEN
    msg = "X!X ERREUR : les charges des cœurs des ions ne sont pas définies, abandon."
  ELSE
    msg = "X!X ERREUR : les charges des coquilles ne sont pas définies, abandon."
  ENDIF
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4817)
  msg = "X!X ERREUR : il n'y a pas de coquille dans le système, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4818)
  !strings(1) = file format (e.g. "xyz", "cfg"...)
  msg = "X!X ERREUR : le fichier "//TRIM(strings(1))//" semble vide,"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "            ou il ne contient aucun nom de fichier valide."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4819)
  !strings(1) = suggestion for a vector
  !strings(2) = direction of suggested vector (X,Y or Z)
  msg = "X!X ERREUR : les vecteurs de base ne sont pas orthogonaux."
  CALL DISPLAY_MSG(1,msg,logfile)
  IF( LEN_TRIM(strings(1))>0 ) THEN
    msg = "    Vecteur suggéré suivant "//TRIM(ADJUSTL(strings(2)))//" : "//TRIM(ADJUSTL(strings(1)))
    CALL DISPLAY_MSG(1,msg,logfile)
  ENDIF
CASE(4820)
  msg = "X!X ERREUR : les dimensions de la supercellule n'ont pas été définis, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4821)
  WRITE(temp,*) NINT(reals(1))
  WRITE(temp2,*) NINT(reals(2))
  msg = "X!X ERREUR : le nombre d'atomes ("//TRIM(ADJUSTL(temp))// &
      & ") dépasse la taille du tableau alloué ("//TRIM(ADJUSTL(temp2))//")."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4822)
  msg = "X!X ERREUR : aucun fichier à traiter, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4823)
  !strings(1) = keyword 1
  !strings(2) = keyword 2
  msg = "X!X ERREUR : les mots-clés '"//TRIM(ADJUSTL(strings(1)))//"' et '" &
      &  //TRIM(ADJUSTL(strings(2)))//"' sont mutuellement exclusifs, abandon."
CASE(4824)
  !strings(1) = name of unknown command
  msg = "X!X ERREUR : commande inconnue : "//TRIM(ADJUSTL(strings(1)))
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Tapez 'help' pour afficher la liste des commandes disponibles."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4825)
  msg = "X!X ERREUR : Atomsk ne peut pas s'exécuter à l'intérieur de lui-même !"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          Vous avez entré 'atomsk' sans argument, ou vous avez cliqué sur l'exécutable"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          depuis un menu, et donc Atomsk s'est exécuté en mode interactif, où seul un"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          nombre limité de commandes sont disponibles, veuillez vous référer à la documentation."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          Pour utiliser Atomsk avec des arguments, vous devez d'abord lancer une"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          invite de commandes, et ensuite y taper votre commande."
  CALL DISPLAY_MSG(1,msg,logfile)
#if defined(WINDOWS)
  msg = "          Sur un système Microsoft Windows, suivez les étapes suivantes :"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          1. Ouvrez le menu Windows, allez dans Accessoires, et lancez Windows Power Shell."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "          2. Exécutez atomsk.exe avec les arguments, par exemple :"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = '             & "C:\Program Files\Atomsk\atomsk.exe" initial.xsf final.cfg'
  CALL DISPLAY_MSG(1,msg,logfile)
#endif
CASE(4826)
  msg = "X!X ERREUR : ce mode n'est pas disponible en mode interactif, veuillez vous référer à la documentation."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4827)
  msg = "X!X ERREUR : impossible de créer le système avec l'orientation donnée."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4828)
  msg = "X!X ERREUR : au moins deux dimensions de la boîte sont trop petites, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4829)
  msg = "X!X ERREUR : les vecteurs de boîte ne sont pas linéairement indépendants, abandon."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(4900)
  msg = "X!X ERREUR : un seul mode peut être utilisé à la fois."
  CALL DISPLAY_MSG(1,msg,logfile)
!
!
!
!Special case for mode PREFERENCES: many questions can be added there
CASE(5000)
  !strings(1) = user's configuration file
  msg = ">>> Les réponses aux questions suivantes déterminent vos préférences"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "..> pour Atomsk, qui seront sauvegardées dans "//TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5001)
  msg = "<?> Voulez-vous que Atomsk écrase les fichiers existants par défaut ?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5002)
  msg = "<?> Voulez-vous que Atomsk ignore les fichiers existants par défaut ?"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5003)
  msg = "<?> Entrez un format de fichier de sortie qui sera toujours activé ('non' pour ignorer) :"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5004)
  msg = "<?> Entrez la langue par défaut ('non' pour ignorer):"
  CALL DISPLAY_MSG(1,msg,logfile)
CASE(5005)
  msg = "<?> Entrez le niveau de verbosité par défaut :"
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    0=silencieux ; aucun message à l'écran ni dans le fichier log."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    1=les messages apparaissent seulement à l'écran (défaut)."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    2=les messages sont seulement écrits dans le fichier log."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    3=les messages sont écrits à l'écran et dans le fichier log."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    4=debug, des messages supplémentaires sont affichés dans le fichier log."
  CALL DISPLAY_MSG(1,msg,logfile)
!
CASE(5700)
  !strings(1) = user's configuration file
  msg = "/!\ ALERTE : il est nécessaire d'écraser le fichier existant " &
      & //TRIM(ADJUSTL(strings(1)))//"."
  CALL DISPLAY_MSG(1,msg,logfile)
  msg = "    Êtes-vous sûr de vouloir poursuivre ? (oui/non)"
  CALL DISPLAY_MSG(1,msg,logfile)
!
!
!
CASE DEFAULT
  CALL ATOMSK_MSG_EN(imsg,strings,reals)
END SELECT
!
END SUBROUTINE ATOMSK_MSG_FR
!
!
!
!********************************************************
! DATE_MSG
! Displays a nice message on certain dates.
!********************************************************
SUBROUTINE DATE_MSG_FR()
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
!Nouvel an
IF(values(2)==1 .AND. values(3)<=10) THEN
  WRITE(msg,*) values(1)
  msg = TRIM(ADJUSTL(msg))
  msg = "*** BONNE ANNÉE "//TRIM(msg)//"!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!14 Mars (3/14) : jour du nombre pi
ELSEIF(values(2)==3 .AND. values(3)==14) THEN
  WRITE(msg,'(a8,f51.48)') "    π = ", pi
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!31 March
ELSEIF(values(2)==3 .AND. values(3)==31) THEN
  msg = "  JOURNÉE MONDIALE DE LA SAUVEGARDE - N'oubliez pas de sauvegarder vos données !"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!1 Avril
ELSEIF(values(2)==4 .AND. values(3)==1) THEN
  msg = "              \,-^--._"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
  msg = "              /`--;-` "
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!1er Mai
ELSEIF(values(2)==5 .AND. values(3)==1) THEN
  msg = "*** Le 1er mai est un jour férié ;-)"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!4 Mai : Star Wars day ("May the fourth/May the force be with you")
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
!14 Juillet
ELSEIF(values(2)==14 .AND. values(3)==07) THEN
  msg = "*** Je prendrais bien une Bastille pour la toux."
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
!Noël
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
  msg = "   * * *   JOYEUX NOËL !   * * *"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!S'il se fait tard ou que c'est le week-end
ELSEIF( values(5)>=20 .OR. values(5)<=6 ) THEN
  msg = "*** Vous travaillez tard ? Il faut dormir parfois :-)"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!À l'heure de manger
ELSEIF(values(5)==12 .AND. values(6)>=30) THEN
  msg = "*** J'ai faim ! :-p"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
!13:37, it's leet time :-)
ELSEIF(values(5)==13 .AND. values(6)==37) THEN
  msg = "*** 17'5 1337 71|V|3!!"
  CALL DISPLAY_MSG(verbosity,msg,logfile)
!
ENDIF
!
END SUBROUTINE DATE_MSG_FR
!
!
!
END MODULE messages_fr
