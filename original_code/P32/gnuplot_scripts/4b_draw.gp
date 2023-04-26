#set terminal postscript portrait enhanced color \
#    font 'Helvetica,20' linewidth 1.5
#set output "asdfasdf.eps"

set terminal pngcairo size 800,600 crop enhanced font 'Verdana,24' linewidth 2
set output "asdfasdf.png"

set size ratio -1

set yrange [-0.5:10.5]
set xrange [-0.5:8.5]

set ytics scale 0 0,1,10

set grid ytics ls -1 lw 2 lc rgb '#808080'

#unset border
set border lw 3 lc rgb "#FEFEFE"
unset xtics
unset key

set format x ""
set format y ""

#set object rectangle from screen 0,0 to screen 1,1 behind fillcolor rgb 'black' fillstyle solid noborder

p 'output/hierarchy_splines' smooth bezier lt 1 lw 4 lc rgb "black", 'output/hierarchy_points' w p pt 7 ps 4 lc rgb "red" notitle
