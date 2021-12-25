set terminal postscript eps enhanced color \
    font 'Helvetica,42' linewidth 1.5

unset key
unset border
set tic scale 0
set boxwidth 1.0 relative
set xrange[-1:12]
set xtics 0,2,10
set grid ytics ls-1 lc rgb "#808080"
set xlabel "Coordination Number"
#set ylabel "Count"
set ytics 20
set ytics format " "
set style fill transparent pattern 4 border

set output "output/statistics/hist_pores_connectivity.eps"
plot 'output/statistics/hist_pores_connectivity' w boxes lc rgb "#404040" notitle, '../P28/hist_connectivity_new' w boxes lc rgb "#404040" notitle,\
'output/statistics/two_hists_connectivity' u 1:($2-$3) w boxes fs solid noborder lc rgb "#404040" notitle

set autoscale x
set xtics autofreq
set boxwidth 1.0 relative
set xlabel "Pore Radius"

set output "output/statistics/hist_pores_r.eps"
plot "output/statistics/hist_pores_r" w boxes lc rgb "#404040" notitle, '../P28/hist_pore_r_new' w boxes lc rgb "#404040" notitle,\
'output/statistics/two_hists_pore_r' u 1:($2-$3) w boxes fs solid noborder lc rgb "#404040" notitle

set boxwidth 1.0 relative
set xlabel "Throat Radius"

set output "output/statistics/hist_throats_r.eps"
plot "output/statistics/hist_throats_r" w boxes lc rgb "#404040" notitle, '../P28/hist_throat_r_new' w boxes lc rgb "#404040" notitle,\
'output/statistics/two_hists_throat_r' u 1:($2-$3) w boxes fs solid noborder lc rgb "#404040" notitle

