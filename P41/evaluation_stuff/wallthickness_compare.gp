set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set output "fig38.png"

unset key
set xrange [0.0:]
set xlabel "Wall Thickness [μm]"
set ylabel "Relative Frequency"
unset ytics

factor = 0.000250501/0.000137943

plot 'wallthickness_sample.txt' u ($1*1000):($2*factor) w boxes fs transparent solid 0.5 noborder lc rgb "#FF0000",\
'wallthickness_model.txt' u ($1*1000):2 w boxes fs transparent solid 0.5 noborder lc rgb "black",
