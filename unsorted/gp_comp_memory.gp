#set terminal postscript eps enhanced color \
#    font 'Helvetica,32' linewidth 2
#set output "asdfasdf.eps"

set terminal pngcairo size 2000,1200 crop enhanced font 'Verdana,50' linewidth 4
set output "asdfasdf.png"

set border lw 2
set key at 39.5,4.75
set key box lw 2
set xtics nomirror
#unset ytics
cMBi(x) = 8/13.0 + 84*x/100/13.0
cMBii(x) = 8/13.0 + 30*x/100/13.0
mMB(x) = 13/13.0 + 8*x/100/13.0
set xrange [0:95]
set yrange [0:5]
set xtics 0,10
set ytics 0,1
set ytic scale 0
set ytics format " "
set grid ytics ls-1 lw 2 lc rgb "#808080"
set xlabel "Porosity [%]"
set ylabel "Total Memory" offset 1.5
plot cMBi(x) title "cMBa (i)" w l lt 1 lw 4 lc rgb "black",\
cMBii(x) title "cMBa (ii)" w l lt 1 lw 4 lc rgb "#808080",\
mMB(x) title "mMBa" w l lt 1 lw 4 lc rgb "red"
