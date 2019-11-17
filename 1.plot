set term png 

set out "z1.png"

set xl "t"
set yl "y(t)"


set title "Metoda Eulera"
p  "1c.dat" u 1:2 w p pt 7 t "dt=1", "1b.dat" u 1:2 w p pt 7 t "dt=0.01", "1a.dat" u 1:2 w p pt 7  lc rgb "red" t "dt=0.1"
set out "z1_errors.png"
p "1a_error.dat" w l t "dt=0.1", "1b_error.dat" w l t "dt=0.01", "1c_error.dat" w l t "dt=1"
set out "z2.png"
p  "2c.dat" u 1:2 w p pt 7 t "dt=1", "2b.dat" u 1:2 w p pt 7 t "dt=0.01", "2a.dat" u 1:2 w p pt 7  lc rgb "red" t "dt=0.1"
set out "z2_errors.png"
p "2a_error.dat" w l t "dt=0.1", "2b_error.dat" w l t "dt=0.01", "2c_error.dat" w l t "dt=1"
set out "z3.png"
p  "3b.dat" u 1:2 w p pt 7 t "dt=0.01", "3a.dat" u 1:2 w p pt 7  lc rgb "red" t "dt=0.1", "3c.dat" u 1:2 w p pt 7 t "dt=1"
set out "z3_errors.png"
p "3a_error.dat" w l t "dt=0.1", "3b_error.dat" w l t "dt=0.01", "3c_error.dat" w l t "dt=1"
set xl "t"
set yl "Q(t)"
set out "z4Q.png"
p  "4aq.dat" u 1:2 w l lc rgb "blue" t "wv=0.5*w0", "4bq.dat" u 1:2 w l lc rgb "red" t "wv=0.8*w0", "4cq.dat" u 1:2 w l lc rgb "green" t "wv=1.0*w0", "4dq.dat" u 1:2 w l lc rgb "violet" t "wv=1.2*w0"
set xl "t"
set yl "I(t)"
set out "z4I.png"
p  "4ai.dat" u 1:2 w l lc rgb "blue" t "wv=0.5*w0", "4bi.dat" u 1:2 w l lc rgb "red" t "wv=0.8*w0", "4ci.dat" u 1:2 w l lc rgb "green" t "wv=1.0*w0", "4di.dat" u 1:2 w l lc rgb "violet" t "wv=1.2*w0"
