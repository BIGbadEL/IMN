set term png 

set out "error_z1_0_6.png"
plot "error_1_0.600000.dat" u 1:2:3 with image


set logscale x

set out "z1.png"

set xl "nr iteracji"
set yl "S"
plot "1_0.600000.dat" w l, "1_1.000000.dat" w l

set out "z2.png"

set xl "nr iteracji"
set yl "S"
plot "2_1.000000.dat" w l , "2_1.400000.dat" w l , "2_1.800000.dat" w l , "2_1.900000.dat" w l


