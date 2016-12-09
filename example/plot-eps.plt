
set term tikz standalone color size 3.5in, 3.5in

set style line 1 lc rgb '#4DAF4A'      lw 2 pt 6  # green
set style line 2 lc rgb '#F781BF' dt 2 lw 2 pt 1  # pink
set style line 3 lc rgb '#A65628' dt 4 lw 2 pt 10 # brown
set style line 4 lc rgb '#E41A1C' dt 3 lw 2 pt 4 # red
set style line 5 lc rgb '#984EA3' dt 5 lw 2 pt 3 # purple
set style line 6 lc rgb '#FF7F00' dt 6 lw 2 pt 2 # orange
set style line 7 lc rgb '#377EB8' dt 7 lw 2 pt 12 # blue
set style line 8 lc rgb '#FFFF33' dt 8 # yellow

set style line 9 lc rgb "#000000" lw 1.5 dt 2
set style line 10 lc rgb "#000000" lw 2 pt 4

set out "eps.tex"

m = 1.1
b = 0.1

f(x) = m*x + b

set mxtics
set mytics
set xtics nomirror out
set ytics out
fit f(x) 'eps1.dat' using 1:2 via m, b
unset key
stats 'eps1.dat' u 1:2 prefix "A"
set label 1 sprintf("y = %12.3fx + %12.3f",A_slope,A_intercept) at -5,18
set label 2 sprintf("R$^{2}$ = %12.3f",A_correlation**2) at -5,16
plot 'eps.dat' u 1:2 ls 1 notitle, \
     'eps1.dat' u 1:2 ls 2 notitle, \
    f(x) ls 9
