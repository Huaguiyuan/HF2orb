reset

set terminal pos eps enhanced color butt "Helvetica" 25##
#set pm3d map
#set size square
#set xtics (-1.0,-0.5,0,0.5,1.0) font ",25"
#set ytics (-1.0,-0.5,0,0.5,1.0) font ",25"
set xlabel "{/Symbol w} - {/Symbol m}" font ",25"
set ylabel "DOS" font ",25"
set xrange [-10:4]
#set yrange [1:2]

set output 'dos2.eps'
plot 'outputTest2' u 1:2 w l noti

reset
