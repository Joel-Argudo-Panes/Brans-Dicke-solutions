# Brans-Dicke-solutions
This Fortran 90 code solves the Brans-Dicke field equations in vacuum for a static and spherically symmetric spacetime.

- Dimensionless variables are used : x is the radial distance in units of the Schwarzschild radius and y = ln(GΦ), where Φ is the scalar field.

- The program uses the RK4 algorithm and employs the weak field solution to set the initial conditions at x0.

- In order to choose a value for x0, the code solves the equations for increasing values of x0 until the variation of the solution near x = 1 (at x_var) with respect to the previous iteration falls below a desired precision (error). In addition to this, two integration steps are used: step1 for x>xlim and step2 for x<xlim, where the functions change more rapidly.

- The equations are solved for diferent values of the BD parameter w in the range [1,300]. The solutions (x, y, the metric components A(x) and 1/B(x), as well as the values of w and x0) are saved in a file and can be plotted with the gnuplot scripts provided.
