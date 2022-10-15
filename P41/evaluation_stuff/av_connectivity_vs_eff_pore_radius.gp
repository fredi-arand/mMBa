set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set output "fig6.png"

unset key
set yrange [sqrt(0.1):]
set xrange [2:300]
set logscale x 10
set logscale y
set xtic scale .5
set ytic scale .5
set xlabel "Effective Radius [μm]"
set ylabel "Average Pore Connectivity"
set grid ytics ls-1 lc rgb "#C0C0C0"
set mxtics 10
set mytics 10
#set ylabel "Counts"
set format y "10^{%T}"

set arrow from 6.0,graph(0,0) to 6.0,graph(1,1) nohead lc rgb "#404040"
set arrow from 2,4.0 to 300,4.0 nohead lc rgb "#404040"

plot '~/build/P41/connectivity_vs_radius.txt' u 1:3:(($2-$1)/2.0):(0) w vectors nohead lw 1.5 lc rgb "black"#,\
# '~/build/P41/connectivity_vs_radius.txt.special' u (($1+$2)/2):3:(($2-$1)/1.9) w boxes fs solid lc rgb "red"
