set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set size ratio -1
unset key
unset border

set xzeroaxis
#set xtics axis
set yzeroaxis
#set ytics axis
unset tics

set xrange[-10:10]
set yrange[-10:10]

set output "fig21.png"

plot '~/build/P41/anisotropy_comparison_9.0-11.0_pores.txt' u (-($1+$2)/2*$4):(-($1+$2)/2*$5):(2*($1+$2)/2*$4):(2*($1+$2)/2*$5) w vectors heads lc rgb "#800000",\
'~/build/P41/anisotropy_comparison_2.7-3.3_pores.txt' u (-($1+$2)/2*$4):(-($1+$2)/2*$5):(2*($1+$2)/2*$4):(2*($1+$2)/2*$5) w vectors heads lc rgb "red",

set output "fig22.png"

plot '~/build/P41/anisotropy_comparison_9.0-11.0_pores.txt' u (-($1+$2)/2*$4):(-($1+$2)/2*$6):(2*($1+$2)/2*$4):(2*($1+$2)/2*$6) w vectors heads lc rgb "#008000",\
'~/build/P41/anisotropy_comparison_2.7-3.3_pores.txt' u (-($1+$2)/2*$4):(-($1+$2)/2*$6):(2*($1+$2)/2*$4):(2*($1+$2)/2*$6) w vectors heads lc rgb "#00FF00"

set output "fig23.png"

plot '~/build/P41/anisotropy_comparison_9.0-11.0_pores.txt' u (-($1+$2)/2*$5):(-($1+$2)/2*$6):(2*($1+$2)/2*$5):(2*($1+$2)/2*$6) w vectors heads lc rgb "#000080",\
'~/build/P41/anisotropy_comparison_2.7-3.3_pores.txt' u (-($1+$2)/2*$5):(-($1+$2)/2*$6):(2*($1+$2)/2*$5):(2*($1+$2)/2*$6) w vectors heads lc rgb "#0000FF"


set xrange[-9:9]
set yrange[-9:9]
unset colorbox

set output "fig24.png"
set palette defined (0 1 0 0, 1 0 0 0)
plot '~/build/P41/anisotropy_comparison_3.0-9.0_pores.txt' u (-($1+$2)/2*$4):(-($1+$2)/2*$5):(2*($1+$2)/2*$4):(2*($1+$2)/2*$5):(($1+$2)/2) w vectors heads palette

set output "fig25.png"
set palette defined (0 0 1 0, 1 0 0 0)
plot '~/build/P41/anisotropy_comparison_3.0-9.0_pores.txt' u (-($1+$2)/2*$4):(-($1+$2)/2*$6):(2*($1+$2)/2*$4):(2*($1+$2)/2*$6):(($1+$2)/2) w vectors heads palette

set output "fig26.png"
set palette defined (0 0 0 1, 1 0 0 0)
plot '~/build/P41/anisotropy_comparison_3.0-9.0_pores.txt' u (-($1+$2)/2*$5):(-($1+$2)/2*$6):(2*($1+$2)/2*$5):(2*($1+$2)/2*$6):(($1+$2)/2) w vectors heads palette
