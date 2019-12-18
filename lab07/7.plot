set term png 

set out "u-1000.png"
set xrange [0:2]
set yrange [0:0.9]
plot "n_s_-1000.000000.dat" u 1:2:5 with image notitle

set out "v-1000.png"
set xrange [0:2]
set yrange [0:0.9]
plot "n_s_-1000.000000.dat" u 1:2:6 with image notitle


set out "u-4000.png"
set xrange [0:2]
set yrange [0:0.9]
plot "n_s_-4000.000000.dat" u 1:2:5 with image notitle

set out "v-4000.png"
set xrange [0:2]
set yrange [0:0.9]
plot "n_s_-4000.000000.dat" u 1:2:6 with image notitle

# set out "2.png"
# set title "nx=ny=50, e1 = e2 = 1"

# plot "V_2.dat" u 1:2:3 with image