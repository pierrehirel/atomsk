MODULE atoms
!
!**********************************************************************************
!*  ATOMS                                                                         *
!**********************************************************************************
!* This module contains subroutines concerning atoms.                             *
!**********************************************************************************
!* (C) Feb. 2014 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 17 Jan. 2025                                     *
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
!* ATOMNUMBER          gives the atomic number of an atom provided its symbol     *
!* ATOMMASS            gives the mass of an atom provided its symbol              *
!* ATOMSPECIES         gives an atom symbol provided its atomic number            *
!* ATOMMASSSPECIES     gives an atom symbol provided its mass                     *
!**********************************************************************************
!
!
USE comv
USE strings
!
PUBLIC:: Element
INTEGER,PARAMETER,PUBLIC:: ATOMMAXZ = 118  !maximum value of atomic number
TYPE:: Element
  CHARACTER(LEN=32):: name   !element name ("hydrogen", "helium", etc.)
  !CHARACTER(LEN=32):: eltype !element type ("alkali metal", "noble gaz", etc.)
  CHARACTER(LEN=2):: sp      !chemical species ("H", "He", etc.)
  INTEGER:: atnumber         !atomic number
  REAL(dp):: atmass          !atomic mass
END TYPE Element
!
!
!
CONTAINS
!
!
!
SUBROUTINE INIT_ELEMENTS()
!< Factory for creating element types
TYPE(Element),DIMENSION(ATOMMAXZ):: Atom
!
! n=1
Atom(1) = Element(name="hydrogen", sp="H", atmass=1.00794d0, atnumber=1)
!       Atom = Element(name="deuterium", sp="D", atmass=2.014101778d0, atnumber=1)
!       Atom = Element(name="tritium", sp="T", atmass=3.0160492675d0, atnumber=1)
Atom(2) = Element(name="helium", sp="He", atmass=4.002602d0, atnumber=2)
!
! n=2
Atom(3) = Element(name="lithium", sp="Li", atmass=6.941d0, atnumber=3)
Atom(4) = Element(name="beryllium", sp="Be", atmass=9.012182d0, atnumber=4)
Atom(5) = Element(name="boron", sp="B", atmass=10.811d0, atnumber=5)
Atom(6) = Element(name="carbon", sp="C", atmass=12.0107d0, atnumber=6)
Atom(7) = Element(name="nitrogen", sp="N", atmass=14.0067d0, atnumber=7)
Atom(8) = Element(name="oxygen", sp="O", atmass=15.9994d0, atnumber=8)
Atom(9) = Element(name="fluorine", sp="F", atmass=18.9984032d0, atnumber=9)
Atom(10) = Element(name="neon", sp="Ne", atmass=20.1797d0, atnumber=10)
!
! n=3
Atom(11) = Element(name="sodium", sp="Na", atmass=22.98977d0, atnumber=11)
Atom(12) = Element(name="magnesium", sp="Mg", atmass=24.305d0, atnumber=12)
Atom(13) = Element(name="aluminum", sp="Al", atmass=26.981538d0, atnumber=13)
Atom(14) = Element(name="silicon", sp="Si", atmass=28.0855d0, atnumber=14)
Atom(15) = Element(name="phosphorus", sp="P", atmass=30.973761d0, atnumber=15)
Atom(16) = Element(name="sulfur", sp="S", atmass=32.065d0, atnumber=16)
Atom(17) = Element(name="chlorine", sp="Cl", atmass=35.453d0, atnumber=17)
Atom(18) = Element(name="argon", sp="Ar", atmass=39.948d0, atnumber=18)
!
! n=4
Atom(19) = Element(name="potassium", sp="K", atmass=39.0983d0, atnumber=19)
Atom(20) = Element(name="calcium", sp="Ca", atmass=40.078d0, atnumber=20)
Atom(21) = Element(name="scandium", sp="Sc", atmass=44.95591d0, atnumber=21)
Atom(22) = Element(name="titanium", sp="Ti", atmass=47.867d0, atnumber=22)
Atom(23) = Element(name="vanadium", sp="V", atmass=50.9415d0, atnumber=23)
Atom(24) = Element(name="chromium", sp="Cr", atmass=51.9961d0, atnumber=24)
Atom(25) = Element(name="manganese", sp="Mn", atmass=54.938049d0, atnumber=25)
Atom(26) = Element(name="iron", sp="Fe", atmass=55.845d0, atnumber=26)
Atom(27) = Element(name="cobalt", sp="Co", atmass=58.9332d0, atnumber=27)
Atom(28) = Element(name="nickel", sp="Ni", atmass=58.6934d0, atnumber=28)
Atom(29) = Element(name="copper", sp="Cu", atmass=63.546d0, atnumber=29)
Atom(30) = Element(name="zinc", sp="Zn", atmass=65.409d0, atnumber=30)
Atom(31) = Element(name="gallium", sp="Ga", atmass=69.723d0, atnumber=31)
Atom(32) = Element(name="germanium", sp="Ge", atmass=72.64d0, atnumber=32)
Atom(33) = Element(name="arsenic", sp="As", atmass=74.9216d0, atnumber=33)
Atom(34) = Element(name="selenium", sp="Se", atmass=78.96d0, atnumber=34)
Atom(35) = Element(name="bromine", sp="Br", atmass=79.904d0, atnumber=35)
Atom(36) = Element(name="krypton", sp="Kr", atmass=83.798d0, atnumber=36)
!
! n=5
Atom(37) = Element(name="rubidium", sp="Rb", atmass=85.4678d0, atnumber=37)
Atom(38) = Element(name="strontium", sp="Sr", atmass=87.62d0, atnumber=38)
Atom(39) = Element(name="yttrium", sp="Y", atmass=88.90585d0, atnumber=39)
Atom(40) = Element(name="zirconium", sp="Zr", atmass=91.224d0, atnumber=40)
Atom(41) = Element(name="niobium", sp="Nb", atmass=92.90638d0, atnumber=41)
Atom(42) = Element(name="molybdenum", sp="Mo", atmass=95.94d0, atnumber=42)
Atom(43) = Element(name="technetium", sp="Tc", atmass=98d0, atnumber=43)
Atom(44) = Element(name="ruthenium", sp="Ru", atmass=101.07d0, atnumber=44)
Atom(45) = Element(name="rhodium", sp="Rh", atmass=102.9055d0, atnumber=45)
Atom(46) = Element(name="palladium", sp="Pd", atmass=106.42d0, atnumber=46)
Atom(47) = Element(name="silver", sp="Ag", atmass=107.8682d0, atnumber=47)
Atom(48) = Element(name="cadmium", sp="Cd", atmass=112.411d0, atnumber=48)
Atom(49) = Element(name="indium", sp="In", atmass=114.818d0, atnumber=49)
Atom(50) = Element(name="tin", sp="Sn", atmass=118.71d0, atnumber=50)
Atom(51) = Element(name="antimony", sp="Sb", atmass=121.76d0, atnumber=51)
Atom(52) = Element(name="tellurium", sp="Te", atmass=127.6d0, atnumber=52)
Atom(53) = Element(name="iodine", sp="I", atmass=126.90447d0, atnumber=53)
Atom(54) = Element(name="xenon", sp="Xe", atmass=131.293d0, atnumber=54)
!
! n=6
Atom(55) = Element(name="cesium", sp="Cs", atmass=132.90545d0, atnumber=55)
Atom(56) = Element(name="barium", sp="Ba", atmass=137.327d0, atnumber=56)
! Lanthanides
Atom(57) = Element(name="lanthanum", sp="La", atmass=138.9055d0, atnumber=57)
Atom(58) = Element(name="cerium", sp="Ce", atmass=140.116d0, atnumber=58)
Atom(59) = Element(name="praseodymium", sp="Pr", atmass=140.90765d0, atnumber=59)
Atom(60) = Element(name="neodymium", sp="Nd", atmass=144.24d0, atnumber=60)
Atom(61) = Element(name="promethium", sp="Pm", atmass=145.d0, atnumber=61)
Atom(62) = Element(name="samarium", sp="Sm", atmass=150.36d0, atnumber=62)
Atom(63) = Element(name="europium", sp="Eu", atmass=151.964d0, atnumber=63)
Atom(64) = Element(name="gadolinium", sp="Gd", atmass=157.25d0, atnumber=64)
Atom(65) = Element(name="terbium", sp="Tb", atmass=158.92534d0, atnumber=65)
Atom(66) = Element(name="dysprosium", sp="Dy", atmass=162.5d0, atnumber=66)
Atom(67) = Element(name="holmium", sp="Ho", atmass=164.93032d0, atnumber=67)
Atom(68) = Element(name="erbium", sp="Er", atmass=167.259d0, atnumber=68)
Atom(69) = Element(name="thulium", sp="Tm", atmass=168.93421d0, atnumber=69)
Atom(70) = Element(name="ytterbium", sp="Yb", atmass=173.04d0, atnumber=70)
Atom(71) = Element(name="lutetium", sp="Lu", atmass=174.967d0, atnumber=71)
! End of lanthanides
Atom(72) = Element(name="hafnium", sp="Hf", atmass=178.49d0, atnumber=72)
Atom(73) = Element(name="tantalum", sp="Ta", atmass=180.9479d0, atnumber=73)
Atom(74) = Element(name="tungsten", sp="W", atmass=183.84d0, atnumber=74)
Atom(75) = Element(name="rhenium", sp="Re", atmass=186.207d0, atnumber=75)
Atom(76) = Element(name="osmium", sp="Os", atmass=190.23d0, atnumber=76)
Atom(77) = Element(name="iridium", sp="Ir", atmass=192.217d0, atnumber=77)
Atom(78) = Element(name="platinum", sp="Pt", atmass=195.078d0, atnumber=78)
Atom(79) = Element(name="gold", sp="Au", atmass=196.96655d0, atnumber=79)
Atom(80) = Element(name="mercury", sp="Hg", atmass=200.59d0, atnumber=80)
Atom(81) = Element(name="thallium", sp="Tl", atmass=204.3833d0, atnumber=81)
Atom(82) = Element(name="lead", sp="Pb", atmass=207.2d0, atnumber=82)
Atom(83) = Element(name="bismuth", sp="Bi", atmass=208.98038d0, atnumber=83)
Atom(84) = Element(name="polonium", sp="Po", atmass=209.d0, atnumber=84)
Atom(85) = Element(name="astatine", sp="At", atmass=210.d0, atnumber=85)
Atom(86) = Element(name="radon", sp="Rn", atmass=222.d0, atnumber=86)
!
!n=7
Atom(87) = Element(name="francium", sp="Fr", atmass=223.d0, atnumber=87)
Atom(88) = Element(name="radium", sp="Ra", atmass=226.d0, atnumber=88)
! Actinides
Atom(89) = Element(name="actinium", sp="Ac", atmass=227.d0, atnumber=89)
Atom(90) = Element(name="thorium", sp="Th", atmass=232.0381d0, atnumber=90)
Atom(91) = Element(name="protactinium", sp="Pa", atmass=231.03588d0, atnumber=91)
Atom(92) = Element(name="uranium", sp="U", atmass=238.02891d0, atnumber=92)
Atom(93) = Element(name="neptunium", sp="Np", atmass=237.d0, atnumber=93)
Atom(94) = Element(name="plutonium", sp="Pu", atmass=244.d0, atnumber=94)
Atom(95) = Element(name="americium", sp="Am", atmass=243.d0, atnumber=95)
Atom(96) = Element(name="curium", sp="Cm", atmass=247.d0, atnumber=96)
Atom(97) = Element(name="berkelium", sp="Bk", atmass=247.d0, atnumber=97)
Atom(98) = Element(name="californium", sp="Cf", atmass=251.d0, atnumber=98)
Atom(99) = Element(name="einsteinium", sp="Es", atmass=252.d0, atnumber=99)
Atom(100) = Element(name="fermium", sp="Fm", atmass=257.d0, atnumber=100)
Atom(101) = Element(name="mendelevium", sp="Md", atmass=258.d0, atnumber=101)
Atom(102) = Element(name="nobelium", sp="No", atmass=259.d0, atnumber=102)
Atom(103) = Element(name="lawrencium", sp="Lr", atmass=262.d0, atnumber=103)
! End of actinides
Atom(104) = Element(name="rutherfordium", sp="Rf", atmass=261.d0, atnumber=104)
Atom(105) = Element(name="dubnium", sp="Db", atmass=262.d0, atnumber=105)
Atom(106) = Element(name="seaborgium", sp="Sg", atmass=266.d0, atnumber=106)
Atom(107) = Element(name="bohrium", sp="Bh", atmass=264.d0, atnumber=107)
Atom(108) = Element(name="hassium", sp="Hs", atmass=277.d0, atnumber=108)
Atom(109) = Element(name="meitnerium", sp="Mt", atmass=268.d0, atnumber=109)
Atom(110) = Element(name="darmstadtium", sp="Ds", atmass=281.d0, atnumber=110)
Atom(111) = Element(name="roentgenium", sp="Rg", atmass=272.d0, atnumber=111)
Atom(112) = Element(name="copernicium", sp="Cn", atmass=285.d0, atnumber=112)
Atom(113) = Element(name="nihonium", sp="Nh", atmass=286.d0, atnumber=113)
Atom(114) = Element(name="flerovium", sp="Fl", atmass=289.d0, atnumber=114)
Atom(115) = Element(name="moscovium", sp="Mc", atmass=289.d0, atnumber=115)
Atom(116) = Element(name="livermorium", sp="Lv", atmass=293.d0, atnumber=116)
Atom(117) = Element(name="tennessine", sp="Ts", atmass=294.d0, atnumber=117)
Atom(118) = Element(name="oganesson", sp="Og", atmass=294.d0, atnumber=118)
!
END SUBROUTINE INIT_ELEMENTS
!
!
!
!********************************************************
!  ATOMNUMBER
!  This subroutine sets the atomic number of an atom
!  depending on its species.
!********************************************************
!
SUBROUTINE ATOMNUMBER(species,snumber)
!
IMPLICIT NONE
CHARACTER(LEN=2):: species, species2
REAL(dp),INTENT(OUT):: snumber
!
!
!Make sure that the first letter is uppercase, the second one lowercase
species2(1:1) = StrUpCase(species(1:1))
species2(2:2) = StrDnCase(species(2:2))
!
SELECT CASE(species2)
! n=1
CASE('H','D')
      snumber=1.d0
CASE('He')
      snumber=2.d0
!
! n=2
CASE('Li')
      snumber=3.d0
CASE('Be')
      snumber=4.d0
CASE('B')
      snumber=5.d0
CASE('C')
      snumber=6.d0
CASE('N')
      snumber=7.d0
CASE('O')
      snumber=8.d0
CASE('F')
      snumber=9.d0
CASE('Ne')
      snumber=10.d0
!
! n=3
CASE('Na')
      snumber=11.d0
CASE('Mg')
      snumber=12.d0
CASE('Al')
      snumber=13.d0
CASE('Si')
      snumber=14.d0
CASE('P')
      snumber=15.d0
CASE('S')
      snumber=16.d0
CASE('Cl')
      snumber=17.d0
CASE('Ar')
      snumber=18.d0
!
! n=4
CASE('K')
      snumber=19.d0
CASE('Ca')
      snumber=20.d0
CASE('Sc')
      snumber=21.d0
CASE('Ti')
      snumber=22.d0
CASE('V')
      snumber=23.d0
CASE('Cr')
      snumber=24.d0
CASE('Mn')
      snumber=25.d0
CASE('Fe')
      snumber=26.d0
CASE('Co')
      snumber=27.d0
CASE('Ni')
      snumber=28.d0
CASE('Cu')
      snumber=29.d0
CASE('Zn')
      snumber=30.d0
CASE('Ga')
      snumber=31.d0
CASE('Ge')
      snumber=32.d0
CASE('As')
      snumber=33.d0
CASE('Se')
      snumber=34.d0
CASE('Br')
      snumber=35.d0
CASE('Kr')
      snumber=36.d0
!
! n=5
CASE('Rb')
      snumber=37.d0
CASE('Sr')
      snumber=38.d0
CASE('Y')
      snumber=39.d0
CASE('Zr') 
      snumber=40.d0
CASE('Nb') 
      snumber=41.d0
CASE('Mo') 
      snumber=42.d0
CASE('Tc') 
      snumber=43.d0
CASE('Ru') 
      snumber=44.d0
CASE('Rh') 
      snumber=45.d0
CASE('Pd') 
      snumber=46.d0
CASE('Ag') 
      snumber=47.d0
CASE('Cd') 
      snumber=48.d0
CASE('In') 
      snumber=49.d0
CASE('Sn') 
      snumber=50.d0
CASE('Sb') 
      snumber=51.d0
CASE('Te') 
      snumber=52.d0
CASE('I') 
      snumber=53.d0
CASE('Xe') 
      snumber=54.d0
!
! n=6
CASE('Cs')
      snumber=55.d0
CASE('Ba') 
      snumber=56.d0
! Lanthanides
CASE('La') 
      snumber=57.d0
CASE('Ce') 
      snumber=58.d0
CASE('Pr') 
      snumber=59.d0
CASE('Nd') 
      snumber=60.d0
CASE('Pm') 
      snumber=61.d0
CASE('Sm') 
      snumber=62.d0
CASE('Eu') 
      snumber=63.d0
CASE('Gd') 
      snumber=64.d0
CASE('Tb') 
      snumber=65.d0
CASE('Dy') 
      snumber=66.d0
CASE('Ho') 
      snumber=67.d0
CASE('Er') 
      snumber=68.d0
CASE('Tm') 
      snumber=69.d0
CASE('Yb') 
      snumber=70.d0
CASE('Lu') 
      snumber=71.d0
! End of Lanthanides
CASE('Hf') 
      snumber=72.d0
CASE('Ta') 
      snumber=73.d0
CASE('W') 
      snumber=74.d0
CASE('Re') 
      snumber=75.d0
CASE('Os') 
      snumber=76.d0
CASE('Ir') 
      snumber=77.d0
CASE('Pt') 
      snumber=78.d0
CASE('Au') 
      snumber=79.d0
CASE('Hg') 
      snumber=80.d0
CASE('Tl') 
      snumber=81.d0
CASE('Pb') 
      snumber=82.d0
CASE('Bi') 
      snumber=83.d0
CASE('Po') 
      snumber=84.d0
CASE('At') 
      snumber=85.d0
CASE('Rn') 
      snumber=86.d0
!
! n=7
CASE('Fr') 
      snumber=87.d0
CASE('Ra') 
      snumber=88.d0
! Actinides
CASE('Ac') 
      snumber=89.d0
CASE('Th') 
      snumber=90.d0
CASE('Pa') 
      snumber=91.d0
CASE('U') 
      snumber=92.d0
CASE('Np') 
      snumber=93.d0
CASE('Pu') 
      snumber=94.d0
CASE('Am') 
      snumber=95.d0
CASE('Cm') 
      snumber=96.d0
CASE('Bk') 
      snumber=97.d0
CASE('Cf') 
      snumber=98.d0
CASE('Es') 
      snumber=99.d0
CASE('Fm') 
      snumber=100.d0
CASE('Md') 
      snumber=101.d0
CASE('No') 
      snumber=102.d0
CASE('Lr') 
      snumber=103.d0
! End of actinides
CASE('Rf') 
      snumber=104.d0
CASE('Db') 
      snumber=105.d0
CASE('Sg') 
      snumber=106.d0
CASE('Bh') 
      snumber=107.d0
CASE('Hs') 
      snumber=108.d0
CASE('Mt') 
      snumber=109.d0
CASE('Ds') 
      snumber=110.d0
CASE('Rg') 
      snumber=111.d0
CASE('Cn')
      snumber=112.d0
CASE('Nh')
      snumber=113.d0
CASE('Fl')
      snumber=114.d0
CASE('Mc')
      snumber=115.d0
CASE('Lv')
      snumber=116.d0
CASE('Ts')
      snumber=117.d0
CASE('Og')
      snumber=118.d0
!
CASE DEFAULT
  !If the species is not recognized
  snumber=0.d0
END SELECT
!
END SUBROUTINE ATOMNUMBER
!
!
!
!********************************************************
!  ATOMMASS
!  This subroutine sets the mass of an atom
!  depending on its species. The masses are from the
!  National Institute of Standards and Technology (NIST):
!  http://www.nist.gov/pml/data/periodic.cfm
!********************************************************
!
SUBROUTINE ATOMMASS(species,smass)
!
IMPLICIT NONE
CHARACTER(LEN=2),INTENT(IN):: species
CHARACTER(LEN=2):: species2
REAL(dp):: smass
!
!
!Make sure that the first letter is uppercase, the second one lowercase
species2(1:1) = StrUpCase(species(1:1))
species2(2:2) = StrDnCase(species(2:2))
!
!
SELECT CASE(species2)
! n=1
CASE('H')
      smass=1.008d0
CASE('D')
      smass=2.014101777d0
CASE('He')
      smass=4.002602d0
!
! n=2
CASE('Li')
      smass=6.94d0
CASE('Be')
      smass=9.012182d0
CASE('B')
      smass=10.81d0
CASE('C')
      smass=12.011d0
CASE('N')
      smass=14.007d0
CASE('O')
      smass=15.999d0
CASE('F')
      smass=18.9984032d0
CASE('Ne')
      smass=20.1797d0
!
! n=3
CASE('Na') 
      smass=22.98976928d0
CASE('Mg') 
      smass=24.305d0
CASE('Al') 
      smass=26.9815386d0
CASE('Si') 
      smass=28.085d0
CASE('P') 
      smass=30.973762d0
CASE('S') 
      smass=32.06d0
CASE('Cl') 
      smass=35.45d0
CASE('Ar') 
      smass=39.948d0
!
! n=4
CASE('K') 
      smass=39.0983d0
CASE('Ca') 
      smass=40.078d0
CASE('Sc') 
      smass=44.955912d0
CASE('Ti') 
      smass=47.867d0
CASE('V') 
      smass=50.9415d0
CASE('Cr') 
      smass=51.9961d0
CASE('Mn') 
      smass=54.938045d0
CASE('Fe') 
      smass=55.845d0
CASE('Co') 
      smass=58.933195d0
CASE('Ni') 
      smass=58.6934d0
CASE('Cu') 
      smass=63.546d0
CASE('Zn') 
      smass=65.38d0
CASE('Ga') 
      smass=69.723d0
CASE('Ge') 
      smass=72.63d0
CASE('As') 
      smass=74.9216d0
CASE('Se') 
      smass=78.96d0
CASE('Br') 
      smass=79.904d0
CASE('Kr') 
      smass=83.798d0
!
! n=5
CASE('Rb') 
      smass=85.4678d0
CASE('Sr') 
      smass=87.62d0
CASE('Y') 
      smass=88.90585d0
CASE('Zr') 
      smass=91.224d0
CASE('Nb') 
      smass=92.90638d0
CASE('Mo') 
      smass=95.96d0
CASE('Tc') 
      smass=98.906d0
CASE('Ru') 
      smass=101.07d0
CASE('Rh') 
      smass=102.90550d0
CASE('Pd') 
      smass=106.42d0
CASE('Ag') 
      smass=107.8682d0
CASE('Cd') 
      smass=112.411d0
CASE('In') 
      smass=114.818d0
CASE('Sn') 
      smass=118.71d0
CASE('Sb') 
      smass=121.76d0
CASE('Te') 
      smass=127.60d0
CASE('I') 
      smass=126.90447d0
CASE('Xe') 
      smass=131.293d0
!
! n=6
CASE('Cs') 
      smass=132.9054519d0
CASE('Ba') 
      smass=137.327d0
! Lanthanides
CASE('La') 
      smass=138.90547d0
CASE('Ce') 
      smass=140.116d0
CASE('Pr') 
      smass=140.90765d0
CASE('Nd') 
      smass=144.242d0
CASE('Pm') 
      smass=144.91d0
CASE('Sm') 
      smass=150.36d0
CASE('Eu') 
      smass=151.964d0
CASE('Gd') 
      smass=157.25d0
CASE('Tb') 
      smass=158.92535d0
CASE('Dy') 
      smass=162.50d0
CASE('Ho') 
      smass=164.93032d0
CASE('Er') 
      smass=167.259d0
CASE('Tm') 
      smass=168.93421d0
CASE('Yb') 
      smass=173.054d0
CASE('Lu') 
      smass=174.9668d0
! End of Lanthanides
CASE('Hf') 
      smass=178.49d0
CASE('Ta') 
      smass=180.94788d0
CASE('W') 
      smass=183.84d0
CASE('Re') 
      smass=186.207d0
CASE('Os') 
      smass=190.23d0
CASE('Ir') 
      smass=192.217d0
CASE('Pt') 
      smass=195.084d0
CASE('Au') 
      smass=196.966569d0
CASE('Hg') 
      smass=200.59d0
CASE('Tl') 
      smass=204.38d0
CASE('Pb') 
      smass=207.2d0
CASE('Bi') 
      smass=208.9804d0
CASE('Po') 
      smass=209.98d0
CASE('At') 
      smass=209.99d0
CASE('Rn') 
      smass=222.02d0
!
! n=7
CASE('Fr') 
      smass=233.d0
CASE('Ra') 
      smass=226.d0
! Actinides
CASE('Ac') 
      smass=227.d0
CASE('Th') 
      smass=232.03806d0
CASE('Pa') 
      smass=231.03588d0
CASE('U') 
      smass=238.02891d0
CASE('Np') 
      smass=237.d0
CASE('Pu') 
      smass=244.d0
CASE('Am') 
      smass=243.d0
CASE('Cm') 
      smass=247.d0
CASE('Bk') 
      smass=247.d0
CASE('Cf') 
      smass=251.d0
CASE('Es') 
      smass=252.d0
CASE('Fm') 
      smass=257.d0
CASE('Md') 
      smass=258.d0
CASE('No') 
      smass=259.d0
CASE('Lr') 
      smass=262.d0
! End of actinides
CASE('Rf') 
      smass=265.d0
CASE('Db') 
      smass=268.d0
CASE('Sg') 
      smass=271.d0
CASE('Bh') 
      smass=270.d0
CASE('Hs') 
      smass=277.d0
CASE('Mt') 
      smass=276.d0
CASE('Ds') 
      smass=281.d0
CASE('Rg') 
      smass=280.d0
CASE('Cn')
      smass=285.17d0
CASE('Uu')
      smass=284.d0
CASE('Fl')
      smass=289.d0
CASE('Mc')
      smass=288.d0
CASE('Lv')
      smass=293.d0
CASE('Ts')
      smass=294.d0
CASE('Og')
      smass=294.d0
!
CASE DEFAULT
  !If the species is not recognized
  smass=0.d0
END SELECT
!
END SUBROUTINE ATOMMASS
!
!
!
!********************************************************
!  ATOMSPECIES
!  This subroutine sets the species of an atom
!  depending on its atomic number
!********************************************************
!
SUBROUTINE ATOMSPECIES(snumber,species)
!
IMPLICIT NONE
REAL(dp),INTENT(IN):: snumber
CHARACTER(LEN=2):: species
INTEGER:: smint
!
smint = INT(snumber)
!
SELECT CASE(smint)
! n=1
CASE(1)
      species='H'
CASE(2)
      species='He'
!
! n=2
CASE(3)
      species='Li'
CASE(4)
      species='Be'
CASE(5)
      species='B'
CASE(6)
      species='C'
CASE(7)
      species='N'
CASE(8)
      species='O'
CASE(9)
      species='F'
CASE(10)
      species='Ne'
!
! n=3
CASE(11)
      species='Na'
CASE(12)
      species='Mg'
CASE(13)
      species='Al'
CASE(14)
      species='Si'
CASE(15)
      species='P'
CASE(16)
      species='S'
CASE(17)
      species='Cl'
CASE(18)
      species='Ar'
!
! n=4
CASE(19)
      species='K'
CASE(20)
      species='Ca'
CASE(21)
      species='Sc'
CASE(22)
      species='Ti'
CASE(23)
      species='V'
CASE(24)
      species='Cr'
CASE(25)
      species='Mn'
CASE(26)
      species='Fe'
CASE(27)
      species='Co'
CASE(28)
      species='Ni'
CASE(29)
      species='Cu'
CASE(30)
      species='Zn'
CASE(31)
      species='Ga'
CASE(32)
      species='Ge'
CASE(33)
      species='As'
CASE(34)
      species='Se'
CASE(35)
      species='Br'
CASE(36)
      species='Kr'
!
! n=5
CASE(37)
      species='Rb'
CASE(38)
      species='Sr'
CASE(39)
      species='Y'
CASE(40)
      species='Zr'
CASE(41)
      species='Nb'
CASE(42)
      species='Mo'
CASE(43)
      species='Tc'
CASE(44)
      species='Ru'
CASE(45)
      species='Rh'
CASE(46)
      species='Pd'
CASE(47)
      species='Ag'
CASE(48)
      species='Cd'
CASE(49)
      species='In'
CASE(50)
      species='Sn'
CASE(51)
      species='Sb'
CASE(52)
      species='Te'
CASE(53)
      species='I'
CASE(54)
      species='Xe'
!
! n=6
CASE(55)
      species='Cs'
CASE(56)
      species='Ba'
! Lanthanides
CASE(57)
      species='La'
CASE(58)
      species='Ce'
CASE(59)
      species='Pr'
CASE(60)
      species='Nd'
CASE(61)
      species='Pm'
CASE(62)
      species='Sm'
CASE(63)
      species='Eu'
CASE(64)
      species='Gd'
CASE(65)
      species='Tb'
CASE(66)
      species='Dy'
CASE(67)
      species='Ho'
CASE(68)
      species='Er'
CASE(69)
      species='Tm'
CASE(70)
      species='Yb'
CASE(71)
      species='Lu'
! End of lanthanides
CASE(72)
      species='Hf'
CASE(73)
      species='Ta'
CASE(74)
      species='W'
CASE(75)
      species='Re'
CASE(76)
      species='Os'
CASE(77)
      species='Ir'
CASE(78)
      species='Pt'
CASE(79)
      species='Au'
CASE(80)
      species='Hg'
CASE(81)
      species='Tl'
CASE(82)
      species='Pb'
CASE(83)
      species='Bi'
CASE(84)
      species='Po'
CASE(85)
      species='At'
CASE(86)
      species='Rn'
!
! n=7
CASE(87)
      species='Fr'
CASE(88)
      species='Ra'
! Actinides
CASE(89)
      species='Ac'
CASE(90)
      species='Th'
CASE(91)
      species='Pa'
CASE(92)
      species='U'
CASE(93)
      species='Np'
CASE(94)
      species='Pu'
CASE(95)
      species='Am'
CASE(96)
      species='Cm'
CASE(97)
      species='Bk'
CASE(98)
      species='Cf'
CASE(99)
      species='Es'
CASE(100)
      species='Fm'
CASE(101)
      species='Md'
CASE(102)
      species='No'
CASE(103)
      species='Lr'
! End of actinides
CASE(104)
      species='Rf'
CASE(105)
      species='Db'
CASE(106)
      species='Sg'
CASE(107)
      species='Bh'
CASE(108)
      species='Hs'
CASE(109)
      species='Mt'
CASE(110)
      species='Ds'
CASE(111)
      species='Rg'
CASE(112)
      species='Cn'
CASE(113)
      species='Nh'
CASE(114)
      species='Fl'
CASE(115)
      species='Mc'
CASE(116)
      species='Lv'
CASE(117)
      species='Ts'
CASE(118)
      species='Og'
!
CASE DEFAULT
  species='XX'
END SELECT
!
END SUBROUTINE ATOMSPECIES
!
!
!
!********************************************************
!  ATOMMASSSPECIES
!  This subroutine sets the species of an atom,
!  depending on its mass. The provided mass is rounded
!  to the nearest integer to avoid numerical imprecision.
!********************************************************
!
SUBROUTINE ATOMMASSSPECIES(smass,species)
!
IMPLICIT NONE
CHARACTER(LEN=2),INTENT(OUT):: species
INTEGER:: mass
REAL(dp),INTENT(IN):: smass
!
!
!Round mass to nearest integer
mass = NINT(smass)
!
!
SELECT CASE(mass)
! n=1
CASE(1)
      species = 'H'
CASE(2)
      species = 'D'
CASE(4)
      species = 'He'
!
! n=2
CASE(6,7)
      species = 'Li'
CASE(9)
      species = 'Be'
CASE(10,11)
      species = 'B'
CASE(12)
      species = 'C'
CASE(14)
      species = 'N'
CASE(16)
      species = 'O'
CASE(19)
      species = 'F'
CASE(20)
      species = 'Ne'
!
! n=3
CASE(22,23) 
      species = 'Na'
CASE(24,25) 
      species = 'Mg'
CASE(26,27) 
      species = 'Al'
CASE(28) 
      species = 'Si'
CASE(30,31)
      species = 'P'
CASE(32)
      species = 'S'
CASE(35,36) 
      species = 'Cl'
!
! n=4
CASE(39)
      species = 'K'
CASE(40)
    !Can be ambiguous: mass(Ar)=39.948; mass(Ca)=40.078
    IF( smass<40.d0 ) THEN
      species = 'Ar'
    ELSE
      species = 'Ca'
    ENDIF
CASE(44,45)
      species = 'Sc'
CASE(47,48)
      species = 'Ti'
CASE(50,51)
      species = 'V'
CASE(52)
      species = 'Cr'
CASE(54,55)
      species = 'Mn'
CASE(56)
      species = 'Fe'
CASE(58)
      species = 'Ni'
CASE(59)
    !Can be ambiguous: mass(Ni)=58.6934; mass(Co)=58.93319500
    IF( smass>58.8d0 ) THEN
      species = 'Co'
    ELSE
      species = 'Ni'
    ENDIF
CASE(63,64)
      species = 'Cu'
CASE(65,66)
      species = 'Zn'
CASE(69,70)
      species = 'Ga'
CASE(72,73)
      species = 'Ge'
CASE(74,75)
      species = 'As'
CASE(78,79)
      species = 'Se'
CASE(80)
      species = 'Br'
CASE(83,84)
      species = 'Kr'
!
! n=5
CASE(85,86)
      species = 'Rb'
CASE(87,88)
      species = 'Sr'
CASE(89)
      species = 'Y'
CASE(91)
      species = 'Zr'
CASE(93)
      species = 'Nb'
CASE(95,96)
      species = 'Mo'
CASE(98,99)
      species = 'Tc'
CASE(101)
      species = 'Ru'
CASE(102,103)
      species = 'Rh'
CASE(106,107)
      species = 'Pd'
CASE(108)
      species = 'Ag'
CASE(112,113)
      species = 'Cd'
CASE(114,115)
      species = 'In'
CASE(118,119)
      species = 'Sn'
CASE(121,122)
      species = 'Sb'
CASE(128)
      species = 'Te'
CASE(126,127)
      species = 'I'
CASE(131)
      species = 'Xe'
!
! n=6
CASE(133)
      species = 'Cs'
CASE(137)
      species = 'Ba'
! Lanthanides
CASE(138,139)
      species = 'La'
CASE(140)
      species = 'Ce'
CASE(141)
      species = 'Pr'
CASE(144)
      species = 'Nd'
CASE(145)
      species = 'Pm'
CASE(150)
      species = 'Sm'
CASE(151,152)
      species = 'Eu'
CASE(157)
      species = 'Gd'
CASE(158,159)
      species = 'Tb'
CASE(162,163)
      species = 'Dy'
CASE(164,165)
      species = 'Ho'
CASE(167)
      species = 'Er'
CASE(168,169)
      species = 'Tm'
CASE(173)
      species = 'Yb'
CASE(175)
      species = 'Lu'
! End of Lanthanides
CASE(178,179)
      species = 'Hf'
CASE(180,181)
      species = 'Ta'
CASE(183,184)
      species = 'W'
CASE(186)
      species = 'Re'
CASE(190)
      species = 'Os'
CASE(192)
      species = 'Ir'
CASE(195)
      species = 'Pt'
CASE(196,197)
      species = 'Au'
CASE(200,201)
      species = 'Hg'
CASE(204)
      species = 'Tl'
CASE(207)
      species = 'Pb'
CASE(208,209)
      species = 'Bi'
CASE(210)
      species = 'Po'
CASE(222)
      species = 'Rn'
!
! n=7
CASE(233)
      species = 'Fr'
CASE(226)
      species = 'Ra'
! Actinides
CASE(227)
      species = 'Ac'
CASE(232)
      species = 'Th'
CASE(231)
      species = 'Pa'
CASE(238)
      species = 'U'
CASE(237)
      species = 'Np'
CASE(244)
      species = 'Pu'
CASE(243)
      species = 'Am'
CASE(247)
      species = 'Cm'
CASE(251)
      species = 'Cf'
CASE(252)
      species = 'Es'
CASE(257)
      species = 'Fm'
CASE(258)
      species = 'Md'
CASE(259)
      species = 'No'
CASE(262)
      species = 'Lr'
! End of actinides
CASE(265)
      species = 'Rf'
CASE(268)
      species = 'Db'
CASE(271)
      species = 'Sg'
CASE(270)
      species = 'Bh'
CASE(277)
      species = 'Hs'
CASE(276)
      species = 'Mt'
CASE(281)
      species = 'Ds'
CASE(280)
      species = 'Rg'
CASE(285)
      species = 'Cn'
CASE(284)
      species = 'Uu'
CASE(289)
      species = 'Fl'
CASE(288)
      species = 'Mc'
CASE(293)
      species = 'Lv'
CASE(294)
      species = 'Ts'
!
CASE DEFAULT
  !If the species is not recognized
  species = "Xx"
END SELECT
!
END SUBROUTINE ATOMMASSSPECIES
!
!
!
END MODULE atoms
