#!/bin/bash 

pHrange=(5.5 7.5)


for pH in ${pHrange[*]}; do

  rm -rf ${pH}
  mkdir ${pH}
  cd ${pH}
  
  sed "s/<value>/${pH}/g" ../system.top > system.top 
  
  mkdir min
  mkdir eq
  mkdir NpT

  cd min 
  
  gmx grompp -f ../../../mdp_files/min.mdp -c ../../start.gro -p ../system.top -maxwarn 2 -v 
  gmx mdrun -v 
  
  cd ../eq
  
  gmx grompp -f ../../../mdp_files/eq.mdp -c ../min/confout.gro -p ../system.top -v -maxwarn 2
  gmx mdrun  -v -nsteps 200000  -nt 10
  gmx grompp -f ../../../mdp_files/eq-p.mdp -c ../eq/confout.gro -p ../system.top -v -maxwarn 2
  gmx mdrun  -v -nsteps 200000  -nt 10
  
  
  cd ../NpT
  
  gmx grompp -f ../../../mdp_files/NpT.mdp -c ../eq/confout.gro -p ../system.top -v  -maxwarn 1
  gmx mdrun -nsteps 1000000 -v -nt 10
  cd ../..

done
