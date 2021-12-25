#set terminal postscript eps enhanced color \
#    font 'Helvetica,32' linewidth 2
#set output "asdfasdf.eps"

set terminal pngcairo size 800,600 crop enhanced font 'Verdana,24' linewidth 2
set output "asdfasdf.png"

set output "asdfasdf.png"
set size ratio -1
unset key
unset tics
set tics scale 0
unset border
#unset colorbox
set format x ""
set format y ""
set palette gray
set palette maxcolors 17
set cbrange[-0.5:16.5]
set cbtics offset -1.0

set xrange [-0.5:2.5]
set yrange [-0.5:2.5]

p "image_greyval" w image, "circles" w lines ls -1 lw 2 lc rgb "red", "gridpoints" w points pt 7 ps 1.5 lc rgb "#808080" notitle
