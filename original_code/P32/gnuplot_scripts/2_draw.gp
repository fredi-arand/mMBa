#set terminal postscript eps size 3.5,3.5 enhanced color \
#    font 'Helvetica,20' linewidth 1.5
#
#set output "asdfasdf.eps"

set terminal pngcairo size 800,600 crop enhanced font 'Verdana,24' linewidth 2
set output "asdfasdf.png"

unset key
set palette gray
set cbrange[-0.5:255.5]
set size ratio -1

set cbtics offset -0.8

set yrange [-0.5:10.5]
set xrange [-0.5:8.5]

set ytics 0,1,10
set xtics 0,1,8
set tic scale 0

p 'output/greyValues' w image, 'output/myCircles0' w circles ls 1 lw 2.5, 'output/myCircles0' w p pt 7 ps 1 lc rgb "red" notitle
