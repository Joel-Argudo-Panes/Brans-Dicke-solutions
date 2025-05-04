# A and 1/B plots.

set term pngcairo
set output "AB_plots.png"

set xrange [0.5:4]
set yrange [-0.5:2]

set xlabel "x"
set ylabel "Metric components"

set key right bottom

f(x) = 1.-1./x # Schwarzschild solution

plot "BD_solutions.dat" index 0 using 1:3 with lines lt rgb "dark-cyan" title "A(x), ω = 1",\
     "BD_solutions.dat" index 0 using 1:4 with lines lt rgb "blue" title "1/B(x), ω = 1",\
     "BD_solutions.dat" index 5 using 1:3 with lines lt rgb "salmon" title "A(x), ω = 6",\
     "BD_solutions.dat" index 5 using 1:4 with lines lt rgb "red" title "1/B(x), ω = 6",\
     "BD_solutions.dat" index 38 using 1:3 with lines lt rgb "light-green"  title "A(x), ω = 300",\
     "BD_solutions.dat" index 38 using 1:4 with lines lt rgb "sea-green" title "1/B(x), ω = 300",\
     f(x) dashtype 2 lt rgb "black" title ""