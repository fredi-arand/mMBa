#set terminal postscript portrait enhanced color \
#    font 'Helvetica,20' linewidth 1.5
#set output "asdfasdf.eps"

set terminal pngcairo size 800,600 crop enhanced font 'Verdana,24' linewidth 2
set output "asdfasdf.png
set size ratio -1
unset key
unset tics
set tics scale 0
unset border
unset colorbox
set format x ""
set format y ""

set yrange [-0.5:23.5]
set xrange [-0.5:15.5]

#set pm3d

p 'output/morphology/detailed/asdf22.png' binary filetype=png w rgbimage, 'output/morphology/detailed/center_radius_22' w circles ls 1 lc rgb "white" lw 2
