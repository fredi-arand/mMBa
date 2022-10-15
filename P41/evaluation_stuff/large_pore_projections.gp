set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set output "fig7.png"
set size ratio -1

unset key
unset border
unset tics
set palette gray
unset colorbox

arrowstart_x = 8
arrowstart_y = 10
arrowlength = 20
set arrow 1 from arrowstart_x,arrowstart_y to arrowstart_x,(arrowstart_y+arrowlength) front lc rgb "white" lw 1.5 nohead
set arrow 2 from arrowstart_x,arrowstart_y to arrowstart_x+arrowlength,arrowstart_y front lc rgb "white" lw 1.5 nohead

letter_xx = -8.2
letter_xy = 3.5
letter_yx = 2.5
letter_yy = -4.2

set view map
# yz
set label 1 "y" at arrowstart_x+arrowlength+4,arrowstart_y front tc rgb "white"
set label 2 "z" at arrowstart_x-4,arrowstart_y+arrowlength+10 front tc rgb "white"
set multiplot
set origin 0,0
set size 1,1
splot '~/build/P41/projection0.mat' matrix w image
set origin 0.19,0.58
set size 0.3,0.3
unset label 1
unset label 2
unset arrow 1
unset arrow 2
splot '~/build/P41/small_pore_projections/projection0.mat' matrix w image
unset multiplot

# xz
set label 1 "x" at arrowstart_x+arrowlength+4,arrowstart_y front tc rgb "white"
set label 2 "z" at arrowstart_x-4,arrowstart_y+arrowlength+10 front tc rgb "white"
set arrow 1 from arrowstart_x,arrowstart_y to arrowstart_x,(arrowstart_y+arrowlength) front lc rgb "white" lw 1.5 nohead
set arrow 2 from arrowstart_x,arrowstart_y to arrowstart_x+arrowlength,arrowstart_y front lc rgb "white" lw 1.5 nohead
set origin 0,0
set size 1,1
set output"fig8.png"
set multiplot
splot '~/build/P41/projection1.mat' matrix w image
set origin 0.19,0.58
set size 0.3,0.3
unset label 1
unset label 2
unset arrow 1
unset arrow 2
splot '~/build/P41/small_pore_projections/projection1.mat' matrix w image
unset multiplot

unset label 1
unset label 2

# xy
set label 1 "x" at arrowstart_x+arrowlength+4,arrowstart_y front tc rgb "white"
set label 2 "y" at arrowstart_x-4,arrowstart_y+arrowlength+10 front tc rgb "white"
set arrow 1 from arrowstart_x,arrowstart_y to arrowstart_x,(arrowstart_y+arrowlength) front lc rgb "white" lw 1.5 nohead
set arrow 2 from arrowstart_x,arrowstart_y to arrowstart_x+arrowlength,arrowstart_y front lc rgb "white" lw 1.5 nohead
set origin 0,0
set size 1,1
set output"fig9.png"
set multiplot
splot '~/build/P41/projection2.mat' matrix w image
set origin 0.19,0.58
set size 0.3,0.3
unset label 1
unset label 2
unset arrow 1
unset arrow 2
splot '~/build/P41/small_pore_projections/projection2.mat' matrix w image
unset multiplot

unset label 1
unset label 2
