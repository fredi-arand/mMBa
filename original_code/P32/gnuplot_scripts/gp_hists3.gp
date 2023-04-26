set terminal postscript eps enhanced color \
    font 'Helvetica,42' linewidth 1.5

unset key
unset border
set tic scale 0
set boxwidth 1.0 relative
set grid ytics ls-1 lc rgb "#808080"
set ytics format " "
set boxwidth 1.0 relative

set autoscale x
set xtics 0,1,3
set xlabel "Difference: Coordination Numbers"
set ytics 100

set output "output/statistics/hist_d_pore_connectivity.eps"
plot "output/statistics/hist_d_pore_connectivity" w boxes fs solid lc rgb "#404040" notitle

set xlabel "Distance: Pore Centers"
set xrange[-0.05:1.45]
set xtics 0,0.3,1.2
set ytics 20
set format x "%2.1f"

set output "output/statistics/hist_d_pore_centers.eps"
plot 'output/statistics/hist_d_pore_centers' w boxes fs solid lc rgb "#404040" notitle

set xtics autofreq
set boxwidth 1.0 relative
set xlabel "Difference: Pore Radii"
set ytics 20
set xrange[-0.95:0.05]
set xtics -0.9,0.3,0

set output "output/statistics/hist_d_pore_radii.eps"
plot "output/statistics/hist_d_pore_radii" w boxes fs solid lc rgb "#404040" notitle

set xlabel "Distance: Throat Centers"
set autoscale x
set xtics autofreq

set output "output/statistics/hist_d_throat_centers.eps"
plot 'output/statistics/hist_d_throat_centers' w boxes fs solid lc rgb "#404040" notitle

set xlabel "Difference: Throat Radii"
set xrange[-0.25:4.25]
set xtics 0,1,4

set output "output/statistics/hist_d_throat_radii.eps"
plot "output/statistics/hist_d_throat_radii" w boxes fs solid lc rgb "#404040" notitle

