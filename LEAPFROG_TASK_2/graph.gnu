#! /usr/bin/gnuplot -persist

set terminal png
set colorbox

stats 'orbital_elements.dat' nooutput

set xlabel "t, лет"
set ylabel "a, а.е."
set output "Simimajor_axis_Earth.png"
plot 'orbital_elements.dat'  index 0:(STATS_blocks-1):2 using 1:2 w l notitle

set xlabel "t, лет"
set ylabel "e"
set output 'Eccentricity_Earth.png'
plot 'orbital_elements.dat'  index 0:(STATS_blocks-1):2 using 1:3 w l notitle

set xlabel "t, лет"
set ylabel "i, градусы"
set output 'Inclination_Earth.png'
plot 'orbital_elements.dat'  index 0:(STATS_blocks-1):2 using 1:4 w l notitle

set xlabel "t, лет"
set ylabel "Periapsis distance, а.е."
set output 'Periapsis_distance_Earth.png'
plot 'orbital_elements.dat'  index 0:(STATS_blocks-1):2 using 1:5 w l notitle



set xlabel "t, лет"
set ylabel "a, а.е."
set output "Simimajor_axis_Jupiter.png"
plot 'orbital_elements.dat'  index 1:(STATS_blocks-1):2 using 1:2 w l notitle

set xlabel "t, лет"
set ylabel "e"
set output 'Eccentricity_Jupiter.png'
plot 'orbital_elements.dat'  index 1:(STATS_blocks-1):2 using 1:3 w l notitle

set xlabel "t, лет"
set ylabel "i, градусы"
set output 'Inclination_Jupiter.png'
plot 'orbital_elements.dat'  index 1:(STATS_blocks-1):2 using 1:4 w l notitle

set xlabel "t, лет"
set ylabel "Periapsis distance, а.е."
set output 'Periapsis_distance_Jupiter.png'
plot 'orbital_elements.dat'  index 1:(STATS_blocks-1):2 using 1:5 w l notitle

set xlabel "t, лет"
set ylabel "Distance between pairs, а.е."
set output 'Distance_at_Sun_to_Earth.png'
plot 'distance_between_pairs.dat'  index 0:(STATS_blocks-1) using 1:2 w l lt 1 notitle

set xlabel "t, лет"
set ylabel "Distance between pairs, а.е."
set output 'Distance_at_Sun_to_Jupiter.png'
plot 'distance_between_pairs.dat'  index 0:(STATS_blocks-1) using 1:3 w l lt 2 notitle

