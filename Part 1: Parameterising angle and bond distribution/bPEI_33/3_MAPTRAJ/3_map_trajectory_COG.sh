#!/bin/bash


read -p "Whats the name of the index file?: " molecule 
echo 0 | gmx trjconv -f ../1_AA/3-run.xtc -o AA-traj.whole.xtc -s ../1_AA/3-run.tpr -pbc whole

gmx grompp -p system-COG.top -f ../1_AA/run.mdp -c ../1_AA/2-eq.gro -o AA-COG.tpr
rm mdout.mdp # clean-up

no_of_beads=$(grep "\[" $molecule.ndx | wc -l)
no_of_beads_minus_1=$( python3 -c "print( $no_of_beads - 1)" )

seq 0 $no_of_beads_minus_1 | gmx traj -f AA-traj.whole.xtc -s AA-COG.tpr -oxt mapped.xtc  -n $molecule.ndx  -ng ${no_of_beads} -com


rm \#*
