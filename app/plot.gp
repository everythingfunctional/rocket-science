set terminal png

set style line 1 linewidth 2
set style line 2 linewidth 2

set format y "%.2g"

set xlabel "Time (s)"

set output 'pressure.png'
set ylabel 'Pressure (Pa)'
plot 'legacy_rocket.out' using 1:2 title 'legacy' with lines linestyle 1, \
  'refurbished_rocket.out' using 1:2 title 'refurbished' with lines linestyle 2

set output 'temperature.png'
set ylabel 'Temperature (K)'
plot 'legacy_rocket.out' using 1:3 title 'legacy' with lines linestyle 1, \
  'refurbished_rocket.out' using 1:3 title 'refurbished' with lines linestyle 2

set output 'mass_flow.png'
set ylabel 'Mass Flow Out of Nozzle (kg/s)'
plot 'legacy_rocket.out' using 1:4 title 'legacy' with lines linestyle 1, \
  'refurbished_rocket.out' using 1:4 title 'refurbished' with lines linestyle 2

set output 'thrust.png'
set ylabel 'Thrust (N)'
plot 'legacy_rocket.out' using 1:5 title 'legacy' with lines linestyle 1, \
  'refurbished_rocket.out' using 1:5 title 'refurbished' with lines linestyle 2

set output 'drag.png'
set ylabel 'Drag (N)'
plot 'legacy_rocket.out' using 1:6 title 'legacy' with lines linestyle 1, \
  'refurbished_rocket.out' using 1:6 title 'refurbished' with lines linestyle 2

set output 'netthrust.png'
set ylabel 'Net Thrust (N)'
plot 'legacy_rocket.out' using 1:7 title 'legacy' with lines linestyle 1, \
  'refurbished_rocket.out' using 1:7 title 'refurbished' with lines linestyle 2

set output 'volume.png'
set ylabel 'Volume of Gas In Chamber (m^3)'
plot 'legacy_rocket.out' using 1:8 title 'legacy' with lines linestyle 1, \
  'refurbished_rocket.out' using 1:8 title 'refurbished' with lines linestyle 2

set output 'acceleration.png'
set ylabel 'Acceleration (m/s^2)'
plot 'legacy_rocket.out' using 1:9 title 'legacy' with lines linestyle 1, \
  'refurbished_rocket.out' using 1:9 title 'refurbished' with lines linestyle 2

set output 'velocity.png'
set ylabel 'Velocity (m/s)'
plot 'legacy_rocket.out' using 1:10 title 'legacy' with lines linestyle 1, \
  'refurbished_rocket.out' using 1:10 title 'refurbished' with lines linestyle 2

set output 'altitude.png'
set ylabel 'Altitude (m)'
plot 'legacy_rocket.out' using 1:11 title 'legacy' with lines linestyle 1, \
  'refurbished_rocket.out' using 1:11 title 'refurbished' with lines linestyle 2
