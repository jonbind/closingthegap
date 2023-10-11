#!/usr/bin/gnuplot -persist
set term pdf dashed enhanced font "Helvetica,12" size 8.25,3.3
set output 'Bond1.pdf'
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
set title "{/Helvetica-Italic b}_1, P-S"
set yrange[0:50]
set xrange[0.2:0.9]
set xtics 0.1
set ytics 25
plot '5_TRGTDISTR/bonds_mapped/distr_bond_0.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/bonds_mapped/distr_bond_0.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"

set title "{/Helvetica-Italic b}_2,  S-T"
set yrange[0:50]
set xrange[0.2:0.9]
set xtics 0.1
set ytics 25
plot '5_TRGTDISTR/bonds_mapped/distr_bond_1.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/bonds_mapped/distr_bond_1.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"

set title "{/Helvetica-Italic b}_3,  P-T"
set yrange[0:50]
set xrange[0.2:0.9]
set xtics 0.1
set ytics 25
plot '5_TRGTDISTR/bonds_mapped/distr_bond_2.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/bonds_mapped/distr_bond_2.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
#
set title "{/Helvetica-Italic b}_4,  S-T"
set yrange[0:50]
set xrange[0.2:0.9]
set xtics 0.1
set ytics 25
plot '5_TRGTDISTR/bonds_mapped/distr_bond_3.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/bonds_mapped/distr_bond_3.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
#
set title "{/Helvetica-Italic b}_5,  S-S"
set yrange[0:50]
set xrange[0.2:0.9]
set xtics 0.1
plot '5_TRGTDISTR/bonds_mapped/distr_bond_4.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/bonds_mapped/distr_bond_4.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
#
set title "{/Helvetica-Italic b}_6, S-S"
set yrange[0:50]
set xrange[0.2:0.9]
set xtics 0.1
plot '5_TRGTDISTR/bonds_mapped/distr_bond_5.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/bonds_mapped/distr_bond_5.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
#
set title "{/Helvetica-Italic b}_7, S-T"
set yrange[0:50]
set xrange[0.2:0.9]
set xtics 0.1
set ytics 25
plot '5_TRGTDISTR/bonds_mapped/distr_bond_6.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/bonds_mapped/distr_bond_6.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"

set title "{/Helvetica-Italic b}_8, S-T"
set yrange[0:50]
set xrange[0.2:0.9]
set xtics 0.1
set ytics 25
plot '5_TRGTDISTR/bonds_mapped/distr_bond_7.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/bonds_mapped/distr_bond_7.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
set title "{/Helvetica-Italic b}_9, P-S"
set yrange[0:50]
set xrange[0.2:0.9]
set xtics 0.1
set ytics 25
plot '5_TRGTDISTR/bonds_mapped/distr_bond_8.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/bonds_mapped/distr_bond_8.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
set title "{/Helvetica-Italic b}_1_0, S-T"
set yrange[0:50]
set xrange[0.2:0.9]
set xtics 0.1
set ytics 25
plot '5_TRGTDISTR/bonds_mapped/distr_bond_9.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/bonds_mapped/distr_bond_9.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
set title "{/Helvetica-Italic b}_1_1, S-S"
set yrange[0:50]
set xrange[0.2:0.9]
set xtics 0.1
plot '5_TRGTDISTR/bonds_mapped/distr_bond_10.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/bonds_mapped/distr_bond_10.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
set title "{/Helvetica-Italic b}_1_2, S-S"
set yrange[0:50]
set xrange[0.2:0.9]
set xtics 0.1
plot '5_TRGTDISTR/bonds_mapped/distr_bond_11.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/bonds_mapped/distr_bond_11.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
set title "{/Helvetica-Italic b}_1_3, P-S"
set yrange[0:50]
set xrange[0.2:0.9]
set xtics 0.1
plot '5_TRGTDISTR/bonds_mapped/distr_bond_12.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/bonds_mapped/distr_bond_12.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
set title "{/Helvetica-Italic b}_1_4,  P-S"
set yrange[0:50]
set xrange[0.2:0.9]
set xtics 0.1
plot '5_TRGTDISTR/bonds_mapped/distr_bond_13.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/bonds_mapped/distr_bond_13.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"

set title "{/Helvetica-Italic b}_1_5, S-S"
set yrange[0:50]
set xrange[0.2:0.9]
set xtics 0.1
plot '5_TRGTDISTR/bonds_mapped/distr_bond_14.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/bonds_mapped/distr_bond_14.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"

set title "{/Helvetica-Italic b}_1_6, S-S"
set yrange[0:50]
set xrange[0.2:0.9]
set xtics 0.1
plot '5_TRGTDISTR/bonds_mapped/distr_bond_15.xvg'  w l lc rgb "#1E90FF" lw 2 title "COG-mapped", \
     '6_CGSIM/bonds_mapped/distr_bond_15.xvg'  w l lc rgb "#DC143C" lw 2 title "Martini"
     

##-------------------------------#
#unset multiplot        # exit multiplot mode (prompt changes back to 'gnuplot')
##-------------------------------#