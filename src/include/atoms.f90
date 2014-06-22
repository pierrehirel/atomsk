MODULE atoms
!
!**********************************************************************************
!*  ATOMS                                                                         *
!**********************************************************************************
!* This module contains subroutines concerning atoms.                             *
!**********************************************************************************
!* (C) Feb. 2014 - Pierre Hirel                                                   *
!*     Unité Matériaux Et Transformations (UMET),                                 *
!*     Université de Lille 1, Bâtiment C6, F-59655 Villeneuve D'Ascq (FRANCE)     *
!*     pierre.hirel@univ-lille1.fr                                                *
!* Last modification: P. Hirel - 19 Feb. 2014                                     *
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
!**********************************************************************************
!
!
USE comv
!
!
CONTAINS
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
CHARACTER(LEN=2):: species
REAL(dp),INTENT(OUT):: snumber
!
!
SELECT CASE(species)
! n=1
CASE('H')
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
CASE('Uu')
      snumber=113.d0
CASE('Fl')
      snumber=114.d0
CASE('Lv')
      snumber=116.d0
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
CHARACTER(LEN=2):: species
REAL(dp):: smass
!
!
SELECT CASE(species)
! n=1
CASE('H')
      smass=1.008d0
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
CASE('Lv')
      smass=293.d0
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
CHARACTER(LEN=2):: species
INTEGER:: smint
REAL(dp):: snumber
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
      species='Uu'
CASE(114)
      species='Fl'
CASE(116)
      species='Lv'
!
CASE DEFAULT
  species='XX'
END SELECT
!
END SUBROUTINE ATOMSPECIES
!
!
!
END MODULE atoms