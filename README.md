# GetRootsF
A Fortran module containing routines for finding the roots of linear/non-linear equations. The methods contained in module "roots" are:
    
    
   root_searchD - Uses exhausive searching for bracketing all the roots within a user-defined range,
                  and then calls the hybridD routine to find the root in each bracket.  The step-size
                  for the search is an optional user input.  If none is supplied, 0.1 is used. The
                  maximum number of iterations is also optional. If none is supplied, 30 iterations
                  are allowed.
    
    
   root_search -  Uses exhausive searching for bracketing all the roots within a user-defined range,
                  and then calls the hybrid routine to find the root in each bracket.  The step-size
                  for the search is an optional user input.  If none is supplied, 0.1 is used. The
                  maximum number of iterations is also optional. If none is supplied, 30 iterations
                  are allowed.


   hybridD - A hybrid Bisection/Newton-Raphson method for finding the root of a
             function F, with derivative Fprime. The root is known to be between "lbound"
             and "rbound", and the result is returned in "root". If the next NR guess
             is within the known bounds, the step is accepted; otherwise, a bisection
             step is taken.


   hybrid -  A hybrid Bisection/Secant method for finding the root of a
             function F when the derivative of F is not available. The Secant method
             uses acceleration to speed convergence. The root is known
             to be between "lbound" and "rbound", and the result is returned in "root".
             If the next SA guess is within the known bounds, the step is accepted;
             otherwise, a bisection step is taken.


This repository contains:
1) roots.f90 - a module containing the root-finding methods.
2) test_roots.f90 - a simple main program to test the module.
3) test_functions.f90 - a library of functions to use with the root finder
4) makefile - a makefie to compile everything.

To run the code:
1)   Clone repository into a directory.
2)  Then run "make".
3) Run "./troots"
