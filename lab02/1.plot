set term png 

set out "z1_a.png"

set xl "t"
set yl "y(t)"

p  "u(t).dat","z(t).dat"

set out "z1_b.png"

set xl "t"
set yl "y(t)"

p  "u(t)_new.dat", "z(t)_new.dat"

