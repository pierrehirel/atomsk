MODULE constants
!
!**********************************************************************************
!*  CONSTANTS                                                                     *
!**********************************************************************************
!* This module contains mathematical and physical constants                       *
!**********************************************************************************
!* (C) Feb. 2010 - Pierre Hirel                                                   *
!*     Université de Lille, Sciences et Technologies                              *
!*     UMR CNRS 8207, UMET - C6, F-59655 Villeneuve D'Ascq, France                *
!*     pierre.hirel@univ-lille.fr                                                 *
!* Last modification: P. Hirel - 03 Oct. 2025                                     *
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
USE comv  !dp is defined there
!
! MATHS
!
REAL(dp),PARAMETER:: zero  = 0.d0
REAL(dp),PARAMETER:: sqrt2 = DSQRT(2.d0)
REAL(dp),PARAMETER:: pi    = 4.d0*DATAN(1.d0)
REAL(dp),DIMENSION(3,3),PARAMETER:: Id_Matrix = &
        & RESHAPE( (/ 1.d0,0.d0,0.d0, 0.d0,1.d0,0.d0, 0.d0,0.d0,1.d0/) , (/3,3/) )
        ! 3x3 identity matrix
!
!
!
! PHYSICS-MATERIALS
!
!General
REAL(dp),PARAMETER:: kB       = 1.3806504d-23            !Boltzmann's constant (J/K)
REAL(dp),PARAMETER:: h_planck = 6.62607015d-34           !Plank's constant (kg.m²/s)
REAL(dp),PARAMETER:: h_bar    = h_planck/(2.d0*pi)       !Reduced Plank's constant (kg.m²/s)
REAL(dp),PARAMETER:: Navo     = 6.022140857d23           !Avogadro constant (/mol)
REAL(dp),PARAMETER:: Temp_0K  = -273.15d0                !Absolute zero of temperature (°C)
!
!Electromagnetism
REAL(dp),PARAMETER:: e_charge = 1.60217646d-19           !Elementary charge (C)
REAL(dp),PARAMETER:: c_light  = 299792458.d0             !Velocity of light (m/s)
REAL(dp),PARAMETER:: mu_0     = 4.d-7*pi                 !vacuum permeability (kg.m/A²/s²)
REAL(dp),PARAMETER:: eps_0    = 1.d0/(mu_0*c_light**2)   !vacuum permittivity (A².s4/kg/m3)
!
!Materials
REAL(dp),PARAMETER:: mass_amu = 1.660538782d-27          !Atomic mass unit (kg)
REAL(dp),PARAMETER:: a_bohr   = 5.291772108d-11          !Bohr radius (m)
REAL(dp),PARAMETER:: E_hartree= 4.35974417d-18           !1 Hartree (J)
REAL(dp),PARAMETER:: E_rydberg= 2.1798741d-18            !1 Rydberg (J)
REAL(dp),PARAMETER:: alpha_fine=mu_0*c_light*e_charge**2                                  &
                 &            /(2.d0*h_planck)           !Fine-structure constant
!
!
!
END MODULE constants
