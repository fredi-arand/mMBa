set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set output "fig10.png"

unset key
set xrange [0:65]
set xlabel "Effective Radius [μm]"
set ylabel "Counts"

#alpha = 1013.18925668
alpha = 790.4

pdf(x) = (x>6 && x<60) ? alpha/(x) : 1/0

plot '~/build/P41/mediumPoreRadii.txt' u (($1+$2)/2):3:(($2-$1)) w boxes fs solid lc rgb "#404040",\
pdf(x) w lines lc rgb "red"
