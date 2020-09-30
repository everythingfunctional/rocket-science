set terminal png

set xlabel "time"

set output 'p.png'
set ylabel 'p'
plot 'legacy_rocket.out' using 1:2

set output 't.png'
set ylabel 't'
plot 'legacy_rocket.out' using 1:3

set output 'mdotos.png'
set ylabel 'mdotos'
plot 'legacy_rocket.out' using 1:4

set output 'thrust.png'
set ylabel 'thrust'
plot 'legacy_rocket.out' using 1:5

set output 'drag.png'
set ylabel 'drag'
plot 'legacy_rocket.out' using 1:6

set output 'netthrust.png'
set ylabel 'netthrust'
plot 'legacy_rocket.out' using 1:7

set output 'vol.png'
set ylabel 'vol'
plot 'legacy_rocket.out' using 1:8

set output 'accel.png'
set ylabel 'accel'
plot 'legacy_rocket.out' using 1:9

set output 'vel.png'
set ylabel 'vel'
plot 'legacy_rocket.out' using 1:10

set output 'altitude.png'
set ylabel 'altitude'
plot 'legacy_rocket.out' using 1:11
