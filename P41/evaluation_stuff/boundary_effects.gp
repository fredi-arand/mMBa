set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set output "fig39.png"

set logscale x
set xrange [0.3:300]
set yrange [0:1]
set xlabel "r [μm]"
set ylabel "p_{L} (r)"
set grid ls 1 lc rgb "#808080"

set key at graph .45,.9 Left box opaque
#set border back
#set key height 1.0
#set key width 1.5

plot ((750-x)/750)**3 lw 2 lc rgb "red" title " L = 1500 μm", ((325-x)/325)**3 lw 2 lc rgb "black" title " L = 750 μm"
