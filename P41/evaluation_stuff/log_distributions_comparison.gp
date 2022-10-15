set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set output "fig12.png"

unset key
set yrange [sqrt(0.1):]
set xrange [0.3:300]
set logscale x 10
set logscale y
set tic scale .5
set grid ytics ls-1 lc rgb "#808080"
set mxtics 10
set mytics 10
set xlabel "Effective Radius [μm]"
set ylabel "Counts"
set format y "10^{%T}"

small_large_seperation = 7.0

set arrow from small_large_seperation,graph(0,0) to small_large_seperation,graph(1,1) nohead front lc rgb "black"

plot '~/build/P41/poreRadiusDistributionLog.txt' u (($1+$2)/2):3:(($2-$1)) w boxes fs transparent solid 0.5 noborder lc rgb "black",\
'~/build/P41/poreRadiusDistributionDownsampledLog.txt' u (($1+$2)):3:(($2-$1)*2) w boxes fs transparent solid 0.5 noborder lc rgb "#FF0000",\
#'~/build/P41/poreRadiusDistributionUpsampledLog.txt' u (($1+$2)/4):3:(($2-$1)/2) w boxes fs transparent solid 0.5 noborder lc rgb "#000080",\
#'~/build/P41/poreRadiusDistributionUpsampled2Log.txt' u (($1+$2)/8):3:(($2-$1)/4) w boxes fs transparent solid 0.5 noborder lc rgb "#800080",\

set output "fig40.png"
plot '~/build/P41/poreRadiusDistributionLog.txt' u (($1+$2)/2):($3/(((325-($1+$2)/2)/325)**3)):($2-$1) w boxes fs transparent solid 0.5 noborder lc rgb "black",\
'~/build/P41/poreRadiusDistributionDownsampledLog.txt' u (($1+$2)):($3/(((750-($1+$2))/750)**3)):(($2-$1)*2) w boxes fs transparent solid 0.5 noborder lc rgb "#FF0000",\
