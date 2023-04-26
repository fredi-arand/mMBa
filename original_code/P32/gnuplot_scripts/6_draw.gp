#set terminal postscript portrait enhanced color \
#    font 'Helvetica,20' linewidth 1.5
#set output "asdfasdf.eps"

set terminal pngcairo size 800,600 crop enhanced font 'Verdana,24' linewidth 2
set output "asdfasdf.png

unset key
set palette gray
set cbrange[-0.5:255.5]
set size ratio -1

set yrange [-0.5:17.5]
set xrange [-0.5:10.5]

#set ytics 0,1,1
#set xtics 0,1,8

#set y2tics -0.5,1,10.5
#set x2tics -0.5,1,8.5

#set grid x2tics y2tics front ls -1 lw 1.5

set format x2 ""
set format y2 ""
set format x ""
set format y ""
unset colorbox

p 'output/greyValues' w image, 'output/myCircles0' w circles ls 1 lw 2.5, 'output/myCircles0' w p pt 7 ps 1 lc rgb "red" notitle
