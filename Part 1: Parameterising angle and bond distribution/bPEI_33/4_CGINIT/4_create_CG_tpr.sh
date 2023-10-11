declare -a StringArray=("bPEI_33")
for val in ${StringArray[@]}; do
   cp $val/2_ATB/"$val"CG.gro $val/4_CGINIT/
   cd  $val/4_CGINIT/
   no_of_beads=$(grep "\[" ../3_MAPTRAJ/$val.ndx | wc -l)
   no_of_beads_minus_1=$( python3 -c "print( $no_of_beads - 1)" )

   seq 0 $no_of_beads_minus_1 | gmx traj -f ../3_MAPTRAJ/AA-traj.whole.xtc -s ../3_MAPTRAJ/AA-COG.tpr \
                                      -n ../3_MAPTRAJ/$val.ndx -oxt molecule.gro -ng ${no_of_beads} -com -e 0

# Generate the CG tpr
   gmx grompp -f martini_v2.x_new_run.mdp -c "$val"-CG.gro -p system.top -o CG.tpr -maxwarn 1

   done