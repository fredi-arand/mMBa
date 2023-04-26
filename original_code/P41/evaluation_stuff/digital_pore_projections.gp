set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set size ratio -1

unset key
unset border
unset tics
set palette gray
unset colorbox


letter_xx = -8.2
letter_xy = 3.5
letter_yx = 2.5
letter_yy = -4.2

set view map
set output "fig34.png"
splot '~/build/P41/projections_digital_foam/projection0.mat' matrix w image

# xz
set output"fig35.png"
splot '~/build/P41/projections_digital_foam/projection1.mat' matrix w image

# xy
set output"fig36.png"
splot '~/build/P41/projections_digital_foam/projection2.mat' matrix w image

unset label 1
unset label 2
