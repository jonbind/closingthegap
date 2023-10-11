pHrange=(3.0 3.5 4.0 4.5 5.0 5.5 6.0 6.5 7.0 7.5)



for pH in ${pHrange[*]}; do
    cp gyranalysis.py $pH/NpT
    cd $pH/NpT
    echo 2|gmx gyrate -f traj_comp.xtc -s topol.tpr -o gyr.xvg 
    python3 gyranalysis.py
    cd ../..
    done 
