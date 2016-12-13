# MATLAB jet color pallete

# line styles
set style line 1 lt 1 lc rgb '#000080' #
set style line 2 lt 1 lc rgb '#0000ff' #
set style line 3 lt 1 lc rgb '#0080ff' #
set style line 4 lt 1 lc rgb '#00ffff' #
set style line 5 lt 1 lc rgb '#80ff80' #
set style line 6 lt 1 lc rgb '#ffff00' #
set style line 7 lt 1 lc rgb '#ff8000' #
set style line 8 lt 1 lc rgb '#ff0000' #
set style line 9 lt 1 lc rgb '#800000' #
# line style used together with jet (<2014b)
set style line 11 lt 1 lc rgb '#0000ff' # blue
set style line 12 lt 1 lc rgb '#007f00' # green
set style line 13 lt 1 lc rgb '#ff0000' # red
set style line 14 lt 1 lc rgb '#00bfbf' # cyan
set style line 15 lt 1 lc rgb '#bf00bf' # pink
set style line 16 lt 1 lc rgb '#bfbf00' # yellow
set style line 17 lt 1 lc rgb '#3f3f3f' # black

# palette
set palette defined (0  0.0 0.0 0.5, \
                     1  0.0 0.0 1.0, \
                     2  0.0 0.5 1.0, \
                     3  0.0 1.0 1.0, \
                     4  0.5 1.0 0.5, \
                     5  1.0 1.0 0.0, \
                     6  1.0 0.5 0.0, \
                     7  1.0 0.0 0.0, \
                     8  0.5 0.0 0.0 )

set term tikz standalone color size 3.5in, 3.5in

unset key
unset colorbox
set xlabel 'x'
set ylabel 'y'
set zlabel 'z'
set view equal xyz
set tics out
set mxtics
set mytics
set mztics
set out 'cluster.tex'
splot 'cluster.dat' u 1:2:3:4 w points palette

set xlabel '$PC_{1}$'
set ylabel '$PC_{2}$'
set out 'pca.tex'
plot 'pca.dat' u 2:3:1 w points palette

set xlabel '$\lambda_{2}$'
set ylabel '$\lambda_{3}$'
set title "time = 0, bandwidth = 1"
set out 'dmap-1.tex'
plot 'dmap-1.dat' u 2:3:1 w points palette

set title "time = 0, bandwidth = 10"
set out 'dmap-10.tex'
plot 'dmap-10.dat' u 2:3:1 w points palette

set title "time = 0, bandwidth = 100"
set out 'dmap-100.tex'
plot 'dmap-100.dat' u 2:3:1 w points palette
