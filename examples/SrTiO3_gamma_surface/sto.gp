#!/usr/bin/gnuplot

### This script aims at displaying a 3-D graph of the gamma-surface.

set encoding iso_8859_1
### Uncomment the 2 following lines to output to an EPS file
set terminal postscript eps enhanced
set output 'gamma.eps'

### Graph parameters
set dgrid3d 60,60 splines
#set pm3d
set samples 40
set isosamples 40
set hidden3d
set view 60,300
set title "GSF energy density of (1@^{/=18-}10) SrTiO_3\n \n LAMMPS, Thomas potential"
set xlabel "[001]" font "Helvetica,20"
set ylabel "[110]" font "Helvetica,20"
set label 1 "{/Symbol g} (J/m^{2})" font "Helvetica,20" at graph 1.6, graph 2.1
set ytics offset 0.0,-0.8
set xyplane 0
set xrange [0:3.905]
set yrange [0:5.5225]
#set zrange [0:2]
set nokey

### The plot itself
file = 'Gamma_tau.dat'
splot file using 1:2:4 with lines

#pause -1

