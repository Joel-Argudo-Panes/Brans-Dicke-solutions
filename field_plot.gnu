# Scalar field for different values of ω.
# Comparison with the weak field solution.

set term pngcairo
set output "field_plot.png"

# Weak field solutions
f(x) = 1./((2.*1.+3.)*x) + (2.*1.+4.)/(2.*1.+3.)
g(x) = 1./((2.*10.+3.)*x) + (2.*10.+4.)/(2.*10.+3.)
h(x) = 1./((2.*100.+3.)*x) + (2.*100.+4.)/(2.*100.+3.)

set xrange [0.5:4]

set xlabel "x"
set ylabel "Gɸ (x)"


plot "BD_solutions.dat" index 0 using 1:(exp($2))  with lines lt rgb "red" title "ω = 1",\
     f(x) dashtype 2 lt rgb  "red" title "weak field, ω = 1",\
     "BD_solutions.dat" index 9 using 1:(exp($2)) with lines lt rgb "blue" title "ω = 10",\
     g(x) dashtype 2 lt rgb  "blue" title "weak field, ω = 10",\
     "BD_solutions.dat" index 18 using 1:(exp($2)) with lines lt rgb "dark-green" title "ω = 100",\
     h(x) dashtype 2 lt rgb  "dark-green" title "weak field, ω = 100"