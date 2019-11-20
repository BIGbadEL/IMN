set term png 

set out "z1_a.png"

set xl "t"
set yl "x(t)"

p  "zad1.dat" u 1:3 w l t "TOL=10^-2","zad1b.dat" u 1:3 w l t "TOL=10^-5"

set out "z1_v.png"

set xl "t"
set yl "v(t)"

p  "zad1.dat" u 1:4 w l t "TOL=10^-2","zad1b.dat" u 1:4 w l t "TOL=10^-5"

set out "z1_t.png"

set xl "t"
set yl "dt(t)"
p  "zad1.dat" u 1:2 w l t "TOL=10^-2","zad1b.dat" u 1:2 w p t "TOL=10^-5"

set out "z1_xv.png"

set xl "x"
set yl "v"
p  "zad1.dat" u 3:4 w l t "TOL=10^-2","zad1b.dat" u 3:4 w l t "TOL=10^-5"
################################################################### RK2 ##########################################

set out "z2_a.png"

set xl "t"
set yl "x(t)"

p  "zad2.dat" u 1:3 w l t "TOL=10^-2","zad2b.dat" u 1:3 w l t "TOL=10^-5"


set out "z2_v.png"

set xl "t"
set yl "v(t)"

p  "zad2.dat" u 1:4 w l t "TOL=10^-2","zad2b.dat" u 1:4 w l t "TOL=10^-5"

set out "z2_t.png"

set xl "t"
set yl "dt(t)"
p  "zad2.dat" u 1:2 w l t "TOL=10^-2","zad2b.dat" u 1:2 w p t "TOL=10^-5"

set out "z2_xv.png"

set xl "x"
set yl "v"
p  "zad2.dat" u 3:4 w l t "TOL=10^-2","zad2b.dat" u 3:4 w l t "TOL=10^-5"
