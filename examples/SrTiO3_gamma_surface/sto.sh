#!/bin/bash
set -e
export LC_NUMERIC=C

# This script computes generalized stacking fault (GSF) energy densities
# in the (1-10) plane of strontium titanate (SrTiO3 or STO), using
# LAMMPS and the interatomic potential by Thomas et al.,
# Nucl. Instrum. Methods Phys. Res. B 228 (2005) 288.
# The system has the orientation X=[001], Y=[110], Z=[1-10].
# A supercell of size 1x1xN is built from it and a SF is introduced
# by shifting the upper half of the supercell by a vector tau.
# Note that in this setup there are 2 stacking faults / supercell.
# The supercell parameters are kept fixed, and ions are relaxed only
# along the normal of the stacking fault (i.e. along the Z axis).

# REQUIREMENTS:
# - programs atomsk, lmp_atom2charge.sh and LAMMPS available in $PATH
# - LAMMPS input script "sto.lin"
# - for plotting the 3-D energy surface, gnuplot is recommended

# Define command aliases
lammps="lmp_local"
if [ "$(which $lammps)" = "" ]; then
  # Try with a different name
  lammps="lammps"
  if [ "$(which $lammps)" = "" ]; then
    printf "X!X ERROR: LAMMPS executable not found or undefined.\n"
    exit
  fi
fi

# Check for additional files in current directory and for executables
if [ ! -e "sto.lin" ]; then
  printf "X!X ERROR: LAMMPS script 'sto.lin' is missing.\n"
  exit
fi
if [ "$(which atomsk)" = "" ]; then
  printf "X!X ERROR: atomsk executable not found or undefined.\n"
  exit
fi
if [ "$(which lmp_atom2charge.sh)" = "" ]; then
  printf "X!X ERROR: script lmp_atom2charge.sh not found or undefined.\n"
  exit
fi
if [ ! -e $(which gnuplot) ]; then
  printf "/!\ WARNING: gnuplot executable not found, no graph will be produced.\n"
fi

# Output files
Gamma_tau="Gamma_tau.dat"    #Whole (1-10) gamma-surface
Gamma_001="Gamma_001.dat"    #gamma profile along [001]
Gamma_110="Gamma_110.dat"    #gamma profile along [110]

rm -f STO_*.lm* STO_*.out log.lammps $Gamma_tau $Gamma_001 $Gamma_110

printf "#tauX \t tauY \t Etot_NoRelax(J/m²) \t Etot_ionRelaxed(J/m²) \n" > $Gamma_tau
printf "#tauX \t Etot_NoRelax(J/m²) \t Etot_ionRelaxed(J/m²) \n" > $Gamma_001
printf "#tauY \t Etot_NoRelax(J/m²) \t Etot_ionRelaxed(J/m²) \n" > $Gamma_110

# Structural parameters
alat=3.9051     #lattice constant a0 of STO (A) (Thomas potential)
S=$(echo "2.0*$alat*$alat*sqrt(2.0)" | bc -l)  #surface of stacking faults
supersize=8    #Number of times the cell is repeated along Z
# Number of steps along X=[001] and Y=[110] directions
# These numbers can be increased to improve accuracy
stepX=10
stepY=16

# Use a loop to produce the different stacking faults and compute their energy

# Loop for shift along X=[001]
for ((i=0;i<=$stepX;i++)) do
  
  # Loop for shift along Y=[110]
  for ((j=0;j<=$stepY;j++)) do

    # Define current shift vector tau
    tauX=$(echo "$i*$alat/$stepX" | bc -l)
    tauY=$(echo "$j*$alat*sqrt(2.0)/$stepY" | bc -l)
    printf ">>> tau=(%.4f,%.4f). " "$tauX" "$tauY"

    # Use atomsk to create the system with stacking fault:
    # Create unit cell oriented X=[001], Y=[110], Z=[1-10]
    # Duplicate it along Z to form a supercell
    # Use the option '-shift' to build the stacking faults
    # Also use the option '-wrap' to ensure all atoms are in the box
    # Output each structure in LAMMPS data format (*.lmp)
    # Run in silent mode (-v 0)
    printf "Building system... "
    atomsk --create perovskite $alat Sr Ti O orient [001] [110] [1-10]  \
            -dup 1 1 $supersize                                         \
            -shift above 0.501*BOX z $tauX $tauY 0.0                    \
            -wrap                                                       \
            STO_${i}_${j}.lmp                                           \
            -v 0 >/dev/null 2>&1

    # Convert data file for "atom_style charge"
    lmp_atom2charge.sh STO_${i}_${j}.lmp >/dev/null 2>&1

    # Set up the LAMMPS input script
    cp sto.lin sto_curr.lin
    sed -i "/read_data / c\read_data   STO_${i}_${j}.lmp" sto_curr.lin
    sed -i "/dump / c\dump   1 all custom 100 STO_${i}_${j}.lmc id type x y z" sto_curr.lin

    # Run LAMMPS
    printf "Running LAMMPS... "
    $lammps  <sto_curr.lin >/dev/null 2>&1

    # Collect energies, write them to file
    energyNR=$(grep -A 1 "Energy initial" log.lammps | tail -n 1 | awk '{print $1}')
    energyR=$(grep -A 1 "Energy initial" log.lammps | tail -n 1 | awk '{print $3}')
    if [[ ( "$i" == "0" && "$j" == "0" ) ]] ; then
      #This defines the zero of energy
      refenergyNR=$(echo "$energyNR" | bc -l)
      refenergyR=$(echo "$energyR" | bc -l)
    fi
    # Compute SF energy density. The factor 16.022 is to convert from eV/A² to J/m²
    ENR=$(echo "16.022*(($energyNR)-($refenergyNR))/$S" | bc -l)
    ER=$(echo "16.022*(($energyR)-($refenergyR))/$S" | bc -l)

    # Output results
    printf "%.4f \t %.4f \t %.4f \t %.4f \n" "$tauX" "$tauY" "$ENR" "$ER" >> $Gamma_tau
    if [[ "$j" == "0" ]] ; then
      printf "%.4f \t %.4f \t %.4f \n" "$tauX" "$ENR" "$ER" >> $Gamma_001
    fi
    if [[ "$i" == "0" ]] ; then
      printf "%.4f \t %.4f \t %.4f \n" "$tauY" "$ENR" "$ER" >> $Gamma_110
    fi

    printf "Done.\n"

    # Remove temporary files
    rm -f sto_curr.lin log.lammps
  done

  printf "\n" >> $Gamma_tau

done

echo ">>> Plotting graph in EPS file 'gamma.eps'..."
gnuplot sto.gp

echo "\o/ Finished."

