set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set output "fig27.png"

set xrange [1:300]
set yrange [1:2]
set logscale x 10
set xtic scale .5
set ytic scale .00001
set xlabel "r_{eff} [μm]"
#set ylabel "<{/Symbol s}_3/{/Symbol s}_2> and <{/Symbol s}_1/{/Symbol s}_2>"
set grid ytics ls-1 lc rgb "#808080"

verticalLinePosition = 7.0
set arrow from verticalLinePosition,graph(0,0) to verticalLinePosition,graph(1,1) nohead lc rgb "#black"

average_ratio = 1.57644
last_bin_end = 2*82.8283

set key left box opaque
set key height 1.0
set key width 1.5
correctionA = 3.0
set ylabel "Ratios"
A = 0.179329119456   
B = 0.892206454705

plot '~/build/P41/anisotropy_comparison_2.0-large_pores.txt' u 1:3:(($2-$1)/correctionA):(0) title "<{/Symbol s}_3/{/Symbol s}_2>" w vectors nohead lw 1.5 lc rgb "black",\
'~/build/P41/anisotropy_comparison_2.0-large_pores.txt' u 1:(1/$7):(($2-$1)/correctionA):(0) title "<{/Symbol s}_2/{/Symbol s}_1>" w vectors nohead lw 1.5 lc rgb "red",\
x>=verticalLinePosition && x<200 ? A*log(x)+B : 1/0 w lines lt 3 lc rgb "#808080" notitle
#unset arrow 1
unset ylabel

set output "for_kai.png"
set key height 1.0
set key width 1.5
correctionA = 3.0
set autoscale y
#set arrow 1 from 6,average_ratio to last_bin_end,average_ratio nohead front lw 1.5 lc rgb "black"
plot '~/build/P41/anisotropy_comparison_2.0-large_pores.txt' u 1:(($1*$3/$7*$3)**(1.0/3.0)):(($2-$1)/correctionA):(0) title "<{/Symbol s}_3>" w vectors nohead lw 1.5 lc rgb "black",\
'~/build/P41/anisotropy_comparison_2.0-large_pores.txt' u 1:(($1/$3*$7)**(1.0/3.0)):(($2-$1)/correctionA):(0) title "<{/Symbol s}_2>" w vectors nohead lw 1.5 lc rgb "red",\
'~/build/P41/anisotropy_comparison_2.0-large_pores.txt' u 1:($1/($1/$3*$7)**(1.0/3.0)/($1*$3/$7*$3)**(1.0/3.0)):(($2-$1)/correctionA):(0) title "<{/Symbol s}_1>" w vectors nohead lw 1.5 lc rgb "blue",\
#unset arrow 1



unset key
unset border
set size ratio -1
unset logscale x

set xzeroaxis
set xtics axis
set yzeroaxis
set ytics axis
unset tics
set xrange[-log10(201):log10(201)]
set yrange[-log10(201):log10(201)]
unset colorbox
set cbrange[0:200]
#set cbtics 0,0.2
set cbtics
set format cb "10^{%1.1f}"
unset xlabel

set output "fig28.png"
r_small_large = 7.0
delta_pm = 2.0
set palette defined (0 0 0 .5, (r_small_large-delta_pm)/200.0 0 0 .5, (r_small_large-delta_pm)/200.0 .5 .5 .5, (r_small_large+delta_pm)/200.0 .5 .5 .5, (r_small_large+delta_pm)/200.0 .5 0 0, 1 .5 0 0)

set label 1 at 2,0.1 "x"
set label 2 at 0.1,2 "y"
plot '~/build/P41/anisotropy_comparison_2.0-large_pores.txt' u (-log10(1+($1+$2)/2)*$4):(-log10(1+($1+$2)/2)*$5):(2*log10(1+($1+$2)/2)*$4):(2*log10(1+($1+$2)/2)*$5):(($1+$2)/2) w vectors heads palette #, -log10(100)<x && x<log10(100) ? 0.576062166804*x : 1/0 w lines lt 3 lw 2 lc rgb "#808080" notitle

unset label 1
unset label 2
set label 1 at 2,0.1 "x"
set label 2 at 0.1,2 "z"
set output "fig29.png"
plot '~/build/P41/anisotropy_comparison_2.0-large_pores.txt' u (-log10(1+($1+$2)/2)*$4):(-log10(1+($1+$2)/2)*$6):(2*log10(1+($1+$2)/2)*$4):(2*log10(1+($1+$2)/2)*$6):(($1+$2)/2) w vectors heads palette

unset label 1
unset label 2
set label 1 at 2,0.1 "y"
set label 2 at 0.1,2 "z"
set output "fig30.png"
plot '~/build/P41/anisotropy_comparison_2.0-large_pores.txt' u (-log10(1+($1+$2)/2)*$5):(-log10(1+($1+$2)/2)*$6):(2*log10(1+($1+$2)/2)*$5):(2*log10(1+($1+$2)/2)*$6):(($1+$2)/2) w vectors heads palette

set output "fig31.png"
set palette defined (0 .5 0 .5, r_small_large/200.0 .5 0 .5, r_small_large/200.0 1 .5 0, 1 1 .5 0)

set label 1 at 2,0.1 "x"
set label 2 at 0.1,2 "y"
plot '~/build/P41/anisotropy_comparison_2.0-large_pores.txt' u (-log10(1+($1+$2)/2)*$8):(-log10(1+($1+$2)/2)*$9):(2*log10(1+($1+$2)/2)*$8):(2*log10(1+($1+$2)/2)*$9):(($1+$2)/2) w vectors heads palette

unset label 1
unset label 2
set label 1 at 2,0.1 "x"
set label 2 at 0.1,2 "z"
set output "fig32.png"
#set palette defined (0 1 0 0, 0.311 1 0 0, 0.311 0 0 1, 1 0 0 1)
plot '~/build/P41/anisotropy_comparison_2.0-large_pores.txt' u (-log10(1+($1+$2)/2)*$8):(-log10(1+($1+$2)/2)*$10):(2*log10(1+($1+$2)/2)*$8):(2*log10(1+($1+$2)/2)*$10):(($1+$2)/2) w vectors heads palette

unset label 1
unset label 2
set label 1 at 2,0.1 "y"
set label 2 at 0.1,2 "z"
set output "fig33.png"
#set palette defined (0 0 1 0, 0.311 0 1 0, 0.311 0 0 1, 1 0 0 1)
plot '~/build/P41/anisotropy_comparison_2.0-large_pores.txt' u (-log10(1+($1+$2)/2)*$9):(-log10(1+($1+$2)/2)*$10):(2*log10(1+($1+$2)/2)*$9):(2*log10(1+($1+$2)/2)*$10):(($1+$2)/2) w vectors heads palette



set output "fig33a.png"
#set palette defined (0 1 0 0, 0.311 1 0 0, 0.311 0 1 0, 1 0 1 0)
set palette defined (0 .5 0 .5, r_small_large/200.0 .5 0 .5, r_small_large/200.0 1 .5 0, 1 1 .5 0)

set label 1 at 2,0.1 "x"
set label 2 at 0.1,2 "y"
plot '~/build/P41/anisotropy_comparison_2.0-large_pores.txt' u (-log10(1+($1+$2)/2)*$11):(-log10(1+($1+$2)/2)*$12):(2*log10(1+($1+$2)/2)*$11):(2*log10(1+($1+$2)/2)*$12):(($1+$2)/2) w vectors heads palette

unset label 1
unset label 2
set label 1 at 2,0.1 "x"
set label 2 at 0.1,2 "z"
set output "fig33b.png"
#set palette defined (0 1 0 0, 0.311 1 0 0, 0.311 0 0 1, 1 0 0 1)
plot '~/build/P41/anisotropy_comparison_2.0-large_pores.txt' u (-log10(1+($1+$2)/2)*$11):(-log10(1+($1+$2)/2)*$13):(2*log10(1+($1+$2)/2)*$11):(2*log10(1+($1+$2)/2)*$13):(($1+$2)/2) w vectors heads palette

unset label 1
unset label 2
set label 1 at 2,0.1 "y"
set label 2 at 0.1,2 "z"
set output "fig33c.png"
#set palette defined (0 0 1 0, 0.311 0 1 0, 0.311 0 0 1, 1 0 0 1)
plot '~/build/P41/anisotropy_comparison_2.0-large_pores.txt' u (-log10(1+($1+$2)/2)*$12):(-log10(1+($1+$2)/2)*$13):(2*log10(1+($1+$2)/2)*$12):(2*log10(1+($1+$2)/2)*$13):(($1+$2)/2) w vectors heads palette

