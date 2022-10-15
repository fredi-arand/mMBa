set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set output "fig11.png"

unset key
set xlabel "Effective Radius [μm]"
set ylabel "Counts"
#set logscale x
set xrange [1.0:7.0]
set mxtics

A = 332.620511032
mu = 0.952824556098
sigma = -0.308947933861
normal_distribution(x) = A/sqrt(2*pi*sigma**2)*exp(-(x-mu)**2/(2*sigma**2))

plot '~/build/P41/pores_smaller_7um.txt' u (($1+$2)/2):($3/((325-($1+$2)/2)/325)**3) w boxes fs transparent solid 0.5 lc rgb "red",\
 '~/build/P41/pores_smaller_7um.txt' u (($1+$2)/2):($3) w boxes fs transparent solid 0.5 lc rgb "#404040",\
normal_distribution(log(x)) lc rgb "black"#,\
