set terminal pngcairo size 1200,800 crop enhanced font 'Verdana,24' linewidth 3
set encoding utf8
set output "fig14.png"

unset key
set xrange [3:300]
set logscale x 10
set tic scale .5
set mxtics 10
set xlabel "Effective Radius [μm]"
set ylabel "Counts"

A = 185.0308382   
mu1 = 3.61959556133   
sigma1 = -0.396208471467
B = 260.463150459   
mu2 = 2.24940045021   
sigma2 = -0.586577678754

C = 195.074814477   
mu3 = 3.51508101605   
mu4 = 2.20957691813   
sigma34 = 0.446259975929

normal_distribution(x,mu,sigma) = 1.0/sqrt(2.0*pi*sigma**2)*exp(-(x-mu)**2/(2*sigma**2))
lognormal_distribution(x,mu,sigma) = 1.0/(x*sigma*sqrt(2*pi))*exp(-(log(x)-mu)**2/(2*sigma**2))
boundary_correction(x) = ((1500-2*x)/1500)**3

plot '~/build/P41/pores_larger_7um_downsampled_log.txt' u ($1+$2):($3/boundary_correction($1+$2)) w boxes fs transparent solid 0.5 lc rgb "red",\
'~/build/P41/pores_larger_7um_downsampled_log.txt' u ($1+$2):3 w boxes fs transparent solid 0.5 lc rgb "#404040",\
x>7.0 ? A*normal_distribution(log(x),mu1,sigma1)+B*normal_distribution(log(x),mu2,sigma2) : 1/0 lc rgb "black"
#x>6.0 ? C*(normal_distribution(log(x),mu3,sigma34)+normal_distribution(log(x),mu4,sigma34)) : 1/0 lc rgb "black",\

