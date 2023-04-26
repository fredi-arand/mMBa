set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set output "~/Dropbox/Promotion/my_papers/digital_carbon_foam/fig46.png"

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

verticalLine = 7.0
set arrow from verticalLine,graph(0,0) to verticalLine,graph(1,1) nohead front lc rgb "black"

plot '~/build/P41/histogram_r_eff_eroded_dilated.txt' u (($1+$2)/2):3 w boxes fs transparent solid 0.5 noborder lc rgb "blue",\
'~/build/P41/artificial_pores_logarithmic.txt' u (($1+$2)/2):3 w boxes fs transparent solid 0.5 noborder lc rgb "black"
