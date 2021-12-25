#set terminal postscript eps size 3.5,3.5 enhanced color \
#    font 'Helvetica,20' linewidth 1.5

#set output "asdfasdf.eps"

set terminal pngcairo size 800,600 crop enhanced font 'Verdana,24' linewidth 2
set output "asdfasdf.png"

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
set cbtics offset -0.6

set yrange [-0.5:10.5]
set xrange [-0.5:8.5]
set tic scale 0

set ytics 0,1,10
set xtics 0,1,8
set format cb "%2.1f"

#p 'output/distanceField' w p pt 7 palette ps 3, 'output/cont.dat' w lines lt 1 lw 2.5
p 'output/distanceField' w image, 'output/cont.dat' w lines lt 1 lw 2.5
