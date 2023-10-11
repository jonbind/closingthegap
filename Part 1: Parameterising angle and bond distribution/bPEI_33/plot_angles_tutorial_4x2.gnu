#!/usr/bin/gnuplot -persist
set term pdf dashed enhanced font "Helvetica,12" size 8.25,3.3
set output 'AAvsCG-ang-tutorial-4x2.pdf'
set bar 1.000000
set border 31 lt -1 lw 1.500
set boxwidth
set style fill empty border
set angles radians
unset grid
unset key
set fit noerrorvariables
#-------------------------------#
set multiplot layout 4,4
#-------------------------------#
set title "{/Helvetica-Italic a}_1, 1 2 3"
set xrange[80:180]
set yrange[0:0.1]
set xtics 20
set ytics 0.025
plot '5_TRGTDISTR/angles_mapped/distr_ang_0.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/angles_mapped/distr_ang_0.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"

set title "{/Helvetica-Italic a}_2, 3 4 5"
set xrange[80:180]
set yrange[0:0.1]
set xtics 20
set ytics 0.025
plot '5_TRGTDISTR/angles_mapped/distr_ang_1.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/angles_mapped/distr_ang_1.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
#
set title "{/Helvetica-Italic a}_3, 3 6 7 "
set xrange[80:180]
set yrange[0:0.1]
set xtics 20
set ytics 0.025
plot '5_TRGTDISTR/angles_mapped/distr_ang_2.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/angles_mapped/distr_ang_2.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
#
set title "{/Helvetica-Italic a}_4, 7 8 9"
set xrange[80:180]
set yrange[0:0.1]
set xtics 20
set ytics 0.025
plot '5_TRGTDISTR/angles_mapped/distr_ang_3.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/angles_mapped/distr_ang_3.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
#
set title "{/Helvetica-Italic a}_5, 8 10 11"
set xrange[80:180]
set yrange[0:0.1]
set xtics 20
set ytics 0.025
plot '5_TRGTDISTR/angles_mapped/distr_ang_4.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/angles_mapped/distr_ang_4.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
#
set title "{/Helvetica-Italic a}_6, 7 12 13"
set xrange[80:180]
set yrange[0:0.1]
set xtics 20
set ytics 0.025
plot '5_TRGTDISTR/angles_mapped/distr_ang_5.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/angles_mapped/distr_ang_5.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
#
set title "{/Helvetica-Italic a}_7, 14 13 15"
set xrange[80:180]
set yrange[0:0.1]
set xtics 20
set ytics 0.025
plot '5_TRGTDISTR/angles_mapped/distr_ang_6.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/angles_mapped/distr_ang_6.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"

set title "{/Helvetica-Italic a}_8, 7 12 16"
set xrange[50:240]
set yrange[0:0.1]
set xtics 20
set ytics 0.025
plot '5_TRGTDISTR/angles_mapped/distr_ang_7.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/angles_mapped/distr_ang_7.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
set title "{/Helvetica-Italic a}_9, 17 16 18"
set xrange[80:180]
set yrange[0:0.1]
set xtics 20
set ytics 0.025
plot '5_TRGTDISTR/angles_mapped/distr_ang_8.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/angles_mapped/distr_ang_8.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"


##-------------------------------#
#unset multiplot        # exit multiplot mode (prompt changes back to 'gnuplot')
##-------------------------------#
#