set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set output "fig15.png"

unset key
set xrange [.3:300]
set yrange [0:2.5]
set logscale x 10
set xtic scale .5
set ytic scale .00001
set xlabel "r_{eff} [μm]"
set ylabel "<{/Symbol s}_3/{/Symbol s}_2> and <{/Symbol s}_1/{/Symbol s}_2>"
set grid ytics ls-1 lc rgb "#808080"

set arrow from 6.0,graph(0,0) to 6.0,graph(1,1) nohead lc rgb "#C0C0C0"
set arrow from 3.0,graph(0,0) to 3.0,graph(1,1) nohead front lc rgb "black"
set arrow from 30.0,graph(0,0) to 30.0,graph(1,1) nohead front lc rgb "black"

average_ratio = 1.57644
last_bin_end = 2*82.8283

correctionA = 10.0
correctionB = 2.0
correctionC = 10.0
#set arrow 1 from 6,average_ratio to last_bin_end,average_ratio nohead front lw 1.5 lc rgb "black"
plot '~/build/P41/anisotropy_comparison_small_pores.txt' u ($1/2):3:(($2-$1)/2/correctionB):(0) w vectors nohead lw 1.5 lc rgb "#FF0000",\
'~/build/P41/anisotropy_comparison_medium_pores.txt' u 1:3:(($2-$1)/correctionC):(0) w vectors nohead lw 1.5 lc rgb "#C00000",\
'~/build/P41/anisotropy_comparison_large_pores.txt' u ($1*2):3:(($2-$1)*2/correctionA):(0) w vectors nohead lw 1.5 lc rgb "#800000",\
#unset arrow 1

#set key left box
#set key height 0.8
#set key width 0.5

set output "fig16.png"
set yrange [0.0:1.0]
set ylabel "Components of <v_3>"
#set arrow 1 from 6,1 to last_bin_end,1 nohead front lw 1.5 lc rgb "black"
set format y "%g"
correction = 10.0
correction2 = 2.0
correction3 = 10.0

plot '~/build/P41/anisotropy_comparison_small_pores.txt' u ($1/2):(abs($4)):((2*$2-2*$1)/4/correction2):(0) title "(v_{3,x})^2" w vectors nohead lw 1.5 lc rgb "#FF0000",\
'~/build/P41/anisotropy_comparison_small_pores.txt' u ($1/2):(abs($5)):((2*$2-2*$1)/4/correction2):(0) title "(v_{3,y})^2" w vectors nohead lw 1.5 lc rgb "#00FF00",\
'~/build/P41/anisotropy_comparison_small_pores.txt' u ($1/2):(abs($6)):((2*$2-2*$1)/4/correction2):(0) title "(v_{3,z})^2" w vectors nohead lw 1.5 lc rgb "#0000FF",\
'~/build/P41/anisotropy_comparison_medium_pores.txt' u 1:(abs($4)):(($2-$1)/correction3):(0) title "(v_{3,x})^2" w vectors nohead lw 1.5 lc rgb "#C00000",\
'~/build/P41/anisotropy_comparison_medium_pores.txt' u 1:(abs($5)):(($2-$1)/correction3):(0) title "(v_{3,y})^2" w vectors nohead lw 1.5 lc rgb "#00C000",\
'~/build/P41/anisotropy_comparison_medium_pores.txt' u 1:(abs($6)):(($2-$1)/correction3):(0) title "(v_{3,z})^2" w vectors nohead lw 1.5 lc rgb "#0000C0",\
'~/build/P41/anisotropy_comparison_large_pores.txt' u (2*$1):(abs($4)):((2*$2-2*$1)/correction):(0) title "(v_{3,x})^2" w vectors nohead lw 1.5 lc rgb "#800000",\
'~/build/P41/anisotropy_comparison_large_pores.txt' u (2*$1):(abs($5)):((2*$2-2*$1)/correction):(0) title "(v_{3,y})^2" w vectors nohead lw 1.5 lc rgb "#008000",\
'~/build/P41/anisotropy_comparison_large_pores.txt' u (2*$1):(abs($6)):((2*$2-2*$1)/correction):(0) title "(v_{3,z})^2" w vectors nohead lw 1.5 lc rgb "#000080"

