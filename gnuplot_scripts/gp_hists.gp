#set terminal postscript eps enhanced color \
#    font 'Helvetica,42' linewidth 1.5

set terminal pngcairo size 800,600 crop enhanced font 'Verdana,28' linewidth 2

unset key
unset border
set ytic scale 0
set boxwidth 1.0 relative
set xtics nomirror
set grid ytics ls-1 lc rgb "#808080"
set ytics 0.2
set ytics format " "
set format x "%2.1f"
set yrange [-0.02:1]

#set xlabel "Coordination Number"
#set output "output/statistics/hist_pores_connectivity.eps"
#plot 'output/statistics/hist_pores_connectivity' w boxes fs solid lc rgb "#404040" notitle

#set autoscale x
#set xtics autofreq nomirror
set xlabel "Pore Radius"
set xrange [25.35:26.25]
set xtics 25.4,0.2,26.2
set output "output/statistics/hist_pores_r.png"
plot "output/statistics/hist_pores_r" w boxes fs solid lc rgb "#404040" notitle


set xlabel "Throat Radius"
#set xrange[0:10]
#set xtics 0,1,13
set xrange[7.9:8.9]
set xtics 8.0,0.2,8.8
set output "output/statistics/hist_throats_r.png"
plot "output/statistics/hist_throats_r" w boxes fs solid lc rgb "#404040" notitle
