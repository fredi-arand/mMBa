set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set output "fig13.png"

unset key
set xlabel "Effective Radius [μm]"
set ylabel "Counts"
#set xrange [0.0:6.0]
set mxtics

numberOfPores = 11151.0
binWidth = abs(3.00796-3.40777)*2
mu = 2.75896
sigma = 0.97306
mu1 = 8
sigma1 = 3
mu2 = 30
sigma2 = 15

plot '~/build/P41/pores_larger_6um_downsampled.txt' u ($1+$2):($3/(binWidth*numberOfPores)) w boxes fs transparent solid 0.5 lc rgb "#404040",\
x>6 && x<90 ? 1.0/x/(log(90.0)-log(6.0)) : 1/0 lc rgb "black"
#x>0 && x< 60 ? 0.5/sqrt(2*pi*sigma1**2)*exp(-(x-mu1)**2/(2*sigma1**2)) + 0.5/sqrt(2*pi*sigma2**2)*exp(-(x-mu2)**2/(2*sigma2**2)) : 1/0 lc rgb "black"
