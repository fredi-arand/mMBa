set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set output "fig5.png"

unset key
set yrange [sqrt(0.1):]
set xrange [1:300]
set logscale x 10
set logscale y
# set tic scale 1e-6
set grid xtics ytics ls-1 lc rgb "#808080"
# set mxtics 10
set xlabel "Effective Radius [μm]"
set ylabel "Counts"

set arrow from 6.0,graph(0,0) to 6.0,graph(1,1) nohead front lc rgb "red"

plot '~/build/P41/poreVolumeDistribution.txt' u (($1+$2)/2):3:(($2-$1)/1.9) w boxes fs solid lc rgb "#404040"
