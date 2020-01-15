set term png 

set out "1.png"
plot "out_0.000000 0.000057.dat" u 1:2, "out_0.100000 0.000057.dat" u 1:2

set out "2.png"
plot "out_0.000000 0.000057.dat" u 1:3, "out_0.100000 0.000057.dat" u 1:3

set out "3.png"
#set title "nx=ny=50, e1 = e2 = 1"
#set xrange [0:4]
#set yrange [0:1]
plot "vx.dat" u 1:2:3 with image

set out "4.png"
#set title "nx=ny=50, e1 = e2 = 1"
#set xrange [0:4]
#set yrange [0:1]
plot "vy.dat" u 1:2:3 with image

reset
set term gif size 800,300 animate delay 10
set output "anim.gif"
n=49    #liczba klatek
set view map # widok z gory
set size ratio -1
set cbr [0:]

do for [i=0:n] {
  file = sprintf("out/zad0.000000_it=%i.txt",i)
  splot file u 1:2:3 w pm3d  title sprintf("t=%i",i)
} 

reset
set term gif size 800,300 animate delay 10
set output "anim1.gif"
n=49    #liczba klatek
set view map # widok z gory
set size ratio -1
set cbr [0:]

do for [i=0:n] {
  file = sprintf("out/zad0.100000_it=%i.txt",i)
  splot file u 1:2:3 w pm3d  title sprintf("t=%i",i)
} 

#set out "2.png"
#set title "nx=ny=50, e1 = e2 = 1"
#set xrange [0:5]
#set yrange [0:5]
#plot "V_2.dat" u 1:2:3 with image

#set out "3.png"
#set title "nx=ny=100, e1 = e2 = 1"
#set xrange [0:10]
#set yrange [0:10]
#plot "V_3.dat" u 1:2:3 with image
#
#set out "4.png"
#set title "nx=ny=200, e1 = e2 = 1"
#set xrange [0:20]
#set yrange [0:20]
#plot "V_4.dat" u 1:2:3 with image
#
#set out "5.png"
#set title "nx=ny=100, e1 = e2 = 1"
#set xrange [0:10]
#set yrange [0:10]
#plot "V_5.dat" u 1:2:3 with image
#
#set out "6.png"
#set title "nx=ny=100, e1 = 1, e2 = 2"
#set xrange [0:10]
#set yrange [0:10]
#plot "V_6.dat" u 1:2:3 with image
#
#set out "7.png"
#set title "nx=ny=100, e1 = 1, e2 = 10"
#set xrange [0:10]
#set yrange [0:10]
#plot "V_7.dat" u 1:2:3 with image