set term png 

set logscale x

set out "z1.png"

set xl "nr iteracji"
set yl "S"
plot "1_0.600000.dat" w l, "1_1.000000.dat" w l

set out "z2.png"

set xl "nr iteracji"
set yl "S"
plot "2_1.000000.dat" w l , "2_1.400000.dat" w l , "2_1.800000.dat" w l , "2_1.900000.dat" w l

