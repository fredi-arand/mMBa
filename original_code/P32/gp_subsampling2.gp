set terminal postscript eps enhanced color \
    font 'Helvetica,32' linewidth 2

set output "asdfasdf.eps"

unset border
unset key
set palette gray
set palette maxcolors 16
set cbrange[-0.5:16.5]

set size ratio -1

unset xtics
unset ytics

set xrange [-0.5:2.5]
set yrange [-0.5:2.5]

plot "image_greyval" w image
