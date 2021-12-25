#set terminal postscript eps size 3.5,2.62 enhanced color \
#    font 'Helvetica,20' linewidth 2
#set output "asdfasdf.eps"

set terminal pngcairo size 800,600 crop enhanced font 'Verdana,30' linewidth 3
set output "asdfasdf.png"
unset key
set tics scale 0

#set lmargin 0
#set rmargin 0

set palette gray
set palette maxcolors 4
set cbrange[-0.5:3.5]
set cbtics offset -0.8
set cbtics  0,1,3
set size ratio -1

set xrange [-0.5:3.5]
set yrange [-0.5:3.5]

set xtics 0,1,3
set ytics 0,1,3

set x2tics -0.5,1,3.5
set y2tics -0.5,1,3.5

set format x2 ""
set format y2 ""

set grid xtics ytics ls -1 lw 1.5

unset border

p 'output/cont.dat' w points pt 7 ps 3.5 notitle, 'output/blub' w p pt 7 palette ps 3.5
