#set terminal postscript portrait enhanced color \
#    font 'Helvetica,20' linewidth 1.5
#set output "asdfasdf.eps"

set terminal pngcairo size 800,600 crop enhanced font 'Verdana,24' linewidth 2
set output "asdfasdf.png

unset key
set size ratio -1
set palette gray
set cbrange[0:3.5]

set ytics scale 0
unset xtics

set yrange [-0.5:10.5]
set xrange [-0.5:8.5]

#unset border
set border lw 3 lc rgb "black"
unset colorbox
set format x ""
set format y ""

p 'output/distanceField' w image, 'output/distanceField' w circles ls 1 lc rgb "red" lw 2
