set term png 

set out "1.png"
plot "V_1.dat" u 1:2:3 with image

set out "2.png"
set title "nx=ny=50, e1 = e2 = 1"
set xrange [0:5]
set yrange [0:5]
plot "V_2.dat" u 1:2:3 with image

set out "3.png"
set title "nx=ny=100, e1 = e2 = 1"
set xrange [0:10]
set yrange [0:10]
plot "V_3.dat" u 1:2:3 with image

set out "4.png"
set title "nx=ny=200, e1 = e2 = 1"
set xrange [0:20]
set yrange [0:20]
plot "V_4.dat" u 1:2:3 with image

set out "5.png"
set title "nx=ny=100, e1 = e2 = 1"
set xrange [0:10]
set yrange [0:10]
plot "V_5.dat" u 1:2:3 with image

set out "6.png"
set title "nx=ny=100, e1 = e2 = 2"
set xrange [0:10]
set yrange [0:10]
plot "V_6.dat" u 1:2:3 with image

set out "7.png"
set title "nx=ny=100, e1 = e2 = 10"
set xrange [0:10]
set yrange [0:10]
plot "V_7.dat" u 1:2:3 with image


# set logscale x

# set out "z1.png"

# set xl "nr iteracji"
# set yl "S"
# plot "1_0.600000.dat" w l, "1_1.000000.dat" w l

# set out "z2.png"

# set xl "nr iteracji"
# set yl "S"
# plot "2_1.000000.dat" w l , "2_1.400000.dat" w l , "2_1.800000.dat" w l , "2_1.900000.dat" w l


