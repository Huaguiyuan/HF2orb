reset

set terminal pos eps enhanced color butt "Helvetica" 25##
set pm3d map
set size square
set xtics (-1.0,-0.5,0,0.5,1.0) font ",25"
set ytics (-1.0,-0.5,0,0.5,1.0) font ",25"
set xlabel "k_x/{/Symbol p}" font ",25"
set ylabel "k_y/{/Symbol p}" font ",25"
set xrange [-1:1]
set yrange [-1:1]
set palette defined (0 "white", 30 "orange", 70 "red", 150 "blue", 200 "cyan", 250 "black") 

#set output 'bandsTBC1.eps'
#splot 'bandsTBC1.dat' u 1:2:4 noti
#set output 'bandsTBC2.eps'
#splot 'bandsTBC2.dat' u 1:2:4 noti
#set output 'bandsTBC3.eps'
#splot 'bandsTBC3.dat' u 1:2:4 noti
set output 'bandsTBC4.eps'
splot 'bandsTBC4.dat' u 1:2:4 noti
#set output 'bandsTBC5.eps'
#splot 'bandsTBC5.dat' u 1:2:4 noti

#set yrange [-2:2]
#set output 'bandsTBC2.eps'
#splot 'bandsTBC2.dat' u 1:3:4 noti

set xlabel "k_x/{/Symbol p}, ky = 0" font ",25"
set ylabel "E - {/Symbol m}" font ",25"
set xrange [-1:1]
set yrange [-2:2]
set ytics (-2,-1,0,1,2) font ",20"
set output 'bandsTBC5.eps'
splot 'bandsTBC5.dat' u 1:3:4 noti

reset
