#set terminal postscript eps size 3.5,3.5 enhanced color \
#    font 'Helvetica,20' linewidth 1.5
#set output "asdfasdf.eps"

set terminal pngcairo size 800,600 crop enhanced font 'Verdana,24' linewidth 2
set output "asdfasdf.png

set contour base
set cntrparam level incremental 0.0, 0.5, 3.5
unset surface
set table 'output/cont.dat'
splot 'output/distanceField'
unset table
reset
unset key
set palette gray
set cbrange[0:3.5]
set size ratio -1

set yrange [-0.5:10.5]
set xrange [-0.5:8.5]

set ytics 0,1,10
set xtics 0,1,8

#set grid xtics ytics ls -1 lw 1.5

#unset border

#p 'output/distanceField' w p pt 7 palette ps 3, 'output/myCircles' w circles ls -1 lw 1.5, 'output/myCircles0' w circles lt 1 lw 2.5, 'output/myCircles0' w p pt 7 ps 1 lc rgb "red" notitle
p 'output/distanceField' w image, 'output/myCircles' w circles ls 1 lc rgb "white" lw 1.5, 'output/myCircles0' w circles lt 1 lc rgb "red" lw 1.5, 'output/myCircles0' w p pt 7 ps 1 lc rgb "red" notitle
