# OSQPSolver
An implementation of the OSQP algorithm (Stellato, 2020).

`OSQPSolver.jl` focuses on reliability and the reuse of the problem container, to reduce memory allocation for repeated usage of this solver as part of an iterative subroutine.

This algorithm yields low-medium accuracy solutions. Solution polishing to get higher accuracy solutions is not implemented presently.

Work in progress.

# Reference
Stellato et. al, OSQP: an operator splitting solver for quadratic programs, 2020. [DOI: 10.1007/s12532-020-00179-2](https://doi.org/10.1007/s12532-020-00179-2).

Also see the official OSQP implementation: [OSQP.jl](https://github.com/osqp/OSQP.jl). Note its Julia interface's update mechanism might have a problem, which is the reason I'm working on this repository.