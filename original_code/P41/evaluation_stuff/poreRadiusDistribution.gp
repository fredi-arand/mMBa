set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set output "fig13.png"

unset key
set xlabel "Effective Radius [μm]"
set ylabel "Counts"
set xrange [sqrt(0.1):]
set yrange [sqrt(0.1):]
set logscale x
set logscale y

pores = 14884
binWidth = 3.1835-1.83779

plot '~/build/P41/poreRadiusDistribution.txt' u (($1+$2)/2):3 w boxes fs transparent solid 0.5 noborder lc rgb "#404040",\
'~/build/P41/poreRadiusDistributionDownsampled.txt' u (($1+$2)):3 w boxes fs transparent solid 0.5 noborder lc rgb "#404040",\
'~/build/P41/poreRadiusDistributionUpsampled.txt' u (($1+$2)/4):3 w boxes fs transparent solid 0.5 noborder lc rgb "#404040",\

