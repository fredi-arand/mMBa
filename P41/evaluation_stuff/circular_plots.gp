set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set size ratio -1
unset key
unset border
unset tics

set output "fig17.png"
set xrange[-9:9]
set yrange[-9:9]
set cbrange [0:1]
set palette defined (0 1 1 1, 1 0.5 0 0)
unset colorbox

set multiplot

plot 'circles.txt' u 1:2:3 w circles lc rgb "black",\
'circle.txt'u 1:2:3 w circles lc rgb "#808080"

plot '~/build/P41/anisotropy_comparison_small_pores.txt' u (($1+$2)/4*$4/($4*$4+$5*$5)):(($1+$2)/4*$5/($4*$4+$5*$5)):(.1):($4*$4+$5*$5) w circles palette fs solid noborder,\
'~/build/P41/anisotropy_comparison_medium_pores.txt' u (($1+$2)/2*$4/($4*$4+$5*$5)):(($1+$2)/2*$5/($4*$4+$5*$5)):(.1):($4*$4+$5*$5) w circles palette fs solid noborder,\
#'~/build/P41/anisotropy_comparison_large_pores.txt' u (($1+$2)/1*$4/($4*$4+$5*$5)):(($1+$2)/1*$5/($4*$4+$5*$5)):(.1):($4*$4+$5*$5) w circles palette fs solid noborder

set palette defined (0 1 1 1, 1 0 0.5 0)

plot '~/build/P41/anisotropy_comparison_small_pores.txt' u (($1+$2)/4*$4/($4*$4+$6*$6)):(($1+$2/4)*$6/($4*$4+$6*$6)):(.1):($4*$4+$6*$6) w circles palette fs solid noborder,\
'~/build/P41/anisotropy_comparison_medium_pores.txt' u (($1+$2)/2*$4/($4*$4+$6*$6)):(($1+$2)/2*$6/($4*$4+$6*$6)):(.1):($4*$4+$6*$6) w circles palette fs solid noborder,\
#'~/build/P41/anisotropy_comparison_large_pores.txt' u (($1+$2)/1*$4/($4*$4+$6*$6)):(($1+$2)/1*$6/($4*$4+$6*$6)):(.1):($4*$4+$6*$6) w circles palette fs solid noborder

set palette defined (0 1 1 1, 1 0 0 0.5)

plot '~/build/P41/anisotropy_comparison_small_pores.txt' u (($1+$2)/4*$5/($5*$5+$6*$6)):(($1+$2)/4*$6/($5*$5+$6*$6)):(.1):($5*$5+$6*$6) w circles palette fs solid noborder,\
'~/build/P41/anisotropy_comparison_medium_pores.txt' u (($1+$2)/2*$5/($5*$5+$6*$6)):(($1+$2)/2*$6/($5*$5+$6*$6)):(.1):($5*$5+$6*$6) w circles palette fs solid noborder,\
#'~/build/P41/anisotropy_comparison_large_pores.txt' u (($1+$2)/1*$5/($5*$5+$6*$6)):(($1+$2)/1*$6/($5*$5+$6*$6)):(.1):($5*$5+$6*$6) w circles palette fs solid noborder

unset multiplot


