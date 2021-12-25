#set terminal postscript eps enhanced color \
#    font 'Helvetica,42' linewidth 1.5

set terminal pngcairo size 1000,600 crop enhanced font 'Verdana,36' linewidth 2

unset key
unset border
set boxwidth 1.0 relative
set style fill solid 0.75
set xtics nomirror
set grid ytics ls-1 lc rgb "#808080"
set format y "10^{%L}"
#set ytics axis
set ytic font ",24"
set ytic offset 0.8,0
set logscale y 10
set ytic scale 0.0001
set yrange [0.31:3100]
set xtic scale 0.0001

#set logscale x

#http://stackoverflow.com/questions/31810947/position-of-tic-marks-in-gnuplot
set xlabel "Pore Radius"
set xrange [-5:175]
set xtics 0,30,150
NXtics=5
Xmax=150
Xmin=0
dX=(Xmax-Xmin)/NXtics
do for [i=0:NXtics] {
  posX = Xmin+i*dX
  set arrow from posX,.215 to posX,.464 nohead front    # bottom
}
set output "hist_pores_r.png"
plot 'statistics/hist_pores_r' w boxes fs solid lc rgb "#404040" notitle
